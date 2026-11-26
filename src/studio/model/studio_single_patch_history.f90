submodule (studio_single_patch) studio_single_patch_history
    use, intrinsic :: iso_fortran_env, only: iostat_end
    use studio_history_storage, only: create_history_blob,open_history_blob,release_history_blob
    use studio_patch_stream, only: write_patch,read_patch
    use studio_project_stream, only: write_text,read_text
    implicit none(type, external)

contains

    module procedure next_patch_revision
        revision_clock=revision_clock+1_int64
        revision=revision_clock
    end procedure next_patch_revision

    module procedure advance_patch_revision

        patch%revision=next_patch_revision()
        patch%mesh_revision=patch%revision
        if(present(mesh_revision))patch%mesh_revision=mesh_revision
    end procedure advance_patch_revision

    module procedure single_patch_can_undo

        single_patch_can_undo=single_patch_undo_stamp(catalog)>0
    end procedure single_patch_can_undo

    module procedure single_patch_undo_stamp

        if(present(catalog))then
            stamp=catalog_undo_stamp(catalog)
        else
            stamp=catalog_undo_stamp(default_catalog)
        end if
    end procedure single_patch_undo_stamp

    integer(int64) function catalog_undo_stamp(catalog) result(stamp)
        type(studio_patch_catalog), intent(in) :: catalog

        stamp=0_int64
        if(catalog%undo_count>0)stamp=usable_undo_stamp(catalog%undo_history(catalog%undo_count)%stamp,catalog)
    end function catalog_undo_stamp

    module procedure usable_undo_stamp

        usable=0_int64
        if(present(catalog))then
            if(stamp>catalog%undo_floor)usable=stamp
        else
            if(stamp>default_catalog%undo_floor)usable=stamp
        end if
    end procedure usable_undo_stamp

    module procedure discard_undo_through
        type(studio_patch_catalog), pointer :: owner

        owner=>edit_catalog(catalog)
        owner%undo_floor=max(owner%undo_floor,stamp)
    end procedure discard_undo_through

    module procedure clear_patch_history
        integer :: i

        do i=1,catalog%undo_count
            call release_patch_undo(catalog%undo_history(i))
        end do
        if(allocated(catalog%undo_history))deallocate(catalog%undo_history)
        catalog%undo_count=0
    end procedure clear_patch_history

    module procedure release_patch_undo

        call release_history_blob(edit%before)
        edit=patch_edit()
    end procedure release_patch_undo

    module procedure release_patch_replay

        call release_patch_undo(replay%edit)
        call release_history_blob(replay%after)
        replay=studio_patch_replay()
    end procedure release_patch_replay

    module procedure get_patch_undo_targets
        integer :: before_counts(3),current_counts(3),rank,id,n,at

        if(.not.single_patch_can_undo(catalog))then
            allocate(targets(2,0))
            return
        end if
        associate(edit=>catalog%undo_history(catalog%undo_count))
            current_counts=[(catalog%patch_count(rank),rank=1,3)]
            before_counts=undo_counts(edit,current_counts)
            n=size(edit%targets,2)
            allocate(targets(2,n+sum(current_counts-before_counts)))
            targets(:,:n)=edit%targets
            at=n
            do rank=1,3
                do id=before_counts(rank)+1,current_counts(rank)
                    at=at+1
                    targets(:,at)=[rank,id]
                end do
            end do
        end associate
    end procedure get_patch_undo_targets

    function undo_counts(edit,current) result(counts)
        type(patch_edit), intent(in) :: edit
        integer, intent(in) :: current(3)
        integer :: counts(3)

        counts=current
        if(all(edit%batch_counts>=0))then
            counts=edit%batch_counts
        else if(edit%created_patch_id>0)then
            counts(edit%kind)=edit%created_patch_id-1
        end if
    end function undo_counts

    module procedure undo_single_patch
        type(studio_patch_catalog), pointer :: owner
        type(patch_entry), allocatable :: entries(:)
        integer, allocatable :: targets(:,:)
        integer :: latest,counts(3),selection(3),rank,status
        logical :: names_only

        kind=0
        error=''
        if(present(operation))operation=''
        if(present(axis))axis=0
        if(present(multipatch_id))multipatch_id=0
        owner=>edit_catalog(catalog)
        if(.not.single_patch_can_undo(owner))return
        latest=owner%undo_count
        associate(edit=>owner%undo_history(latest))
            names_only=edit%operation=='rename'
            error='Could not read geometry Undo data.'
            call read_history_state(edit%before,edit%targets,names_only,entries,status)
            if(status/=0)return
            counts=[(owner%patch_count(rank),rank=1,3)]
            if(present(replay))then
                call get_patch_undo_targets(owner,targets)
                error='Could not write geometry Redo data.'
                call write_current_state(owner,targets,names_only,replay%after,status)
                if(status/=0)return
                replay%counts=counts
                replay%selection=owner%selection
                call move_alloc(targets,replay%targets)
                replay%edit=edit
            end if

            selection=owner%selection
            if(all(edit%batch_counts>=0))then
                selection=edit%batch_selection
            else if(edit%previous_selection>0)then
                selection(edit%kind)=edit%previous_selection
            else if(edit%select_restored)then
                selection(edit%kind)=edit%targets(2,1)
            end if

            call apply_history_state(owner,undo_counts(edit,counts),selection,entries,names_only)
            kind=edit%kind
            if(present(operation))operation=edit%operation
            if(present(axis))axis=edit%axis
            if(present(multipatch_id))multipatch_id=edit%multipatch_id
            if(.not.present(replay))call release_history_blob(edit%before)
            edit=patch_edit()
        end associate
        owner%undo_count=latest-1
        error=''
    end procedure undo_single_patch

    module procedure redo_single_patch
        type(studio_patch_catalog), pointer :: owner
        type(patch_entry), allocatable :: entries(:)
        logical :: names_only
        integer :: status

        error='There is no geometry Redo entry.'
        if(replay%edit%stamp==0_int64)return
        owner=>edit_catalog(catalog)
        names_only=replay%edit%operation=='rename'
        error='Could not read geometry Redo data.'
        call read_history_state(replay%after,replay%targets,names_only,entries,status)
        if(status/=0)return
        call apply_history_state(owner,replay%counts,replay%selection,entries,names_only)
        call append_patch_undo(owner,replay%edit)
        call release_patch_replay(replay)
        error=''
    end procedure redo_single_patch

    module procedure record_undo
        type(studio_patch_catalog), pointer :: owner
        type(patch_edit) :: stored

        owner=>edit_catalog(catalog)
        call prepare_patch_undo(edit,owner,stored,error)
        if(len(error)>0)return
        call append_patch_undo(owner,stored)
    end procedure record_undo

    module procedure prepare_patch_undo
        integer :: count,i,unit,status,closed
        logical :: names_only

        stored=patch_edit_metadata(edit)
        names_only=edit%operation=='rename'
        count=0
        if(edit%previous%patch_id>0)then
            count=1
            stored%select_restored=.true.
        end if
        if(allocated(edit%previous_set))count=size(edit%previous_set)
        allocate(stored%targets(2,count))
        error=''
        if(count==0)return

        error='Could not save geometry Undo data.'
        call create_history_blob(stored%before,unit,status,patch_history_project(catalog))
        if(status/=0)return
        if(allocated(edit%previous_set))then
            do i=1,count
                associate(entry=>edit%previous_set(i))
                    stored%targets(:,i)=[entry%dimension,entry%patch_id]
                    call write_history_patch(unit,entry,names_only,status)
                end associate
                if(status/=0)exit
            end do
        else
            stored%targets(:,1)=[edit%previous%dimension,edit%previous%patch_id]
            call write_history_patch(unit,edit%previous,names_only,status)
        end if
        close(unit,iostat=closed)
        if(status==0.and.closed==0)then
            error=''
        else
            call release_patch_undo(stored)
        end if
    end procedure prepare_patch_undo

    function patch_edit_metadata(edit) result(stored)
        type(patch_edit), intent(in) :: edit
        type(patch_edit) :: stored

        stored=patch_edit(kind=edit%kind,axis=edit%axis,created_patch_id=edit%created_patch_id, &
            previous_selection=edit%previous_selection,multipatch_id=edit%multipatch_id, &
            operation=edit%operation,batch_counts=edit%batch_counts,batch_selection=edit%batch_selection)
    end function patch_edit_metadata

    module procedure prepare_patch_set_undo
        integer :: status

        stored=patch_edit_metadata(edit)
        stored%targets=targets
        call write_current_state(catalog,targets,.false.,stored%before,status)
        error=''
        if(status/=0)error='Could not save geometry Undo data.'
    end procedure prepare_patch_set_undo

    module procedure prepare_patch_snapshot_undo
        integer :: i,status,unit,closed

        stored=patch_edit_metadata(edit)
        stored%targets=targets
        error='Could not save geometry Undo data.'
        call create_history_blob(stored%before,unit,status,patch_history_project(catalog))
        if(status/=0)return
        do i=1,size(entries)
            call write_history_patch(unit,entries(i),.false.,status)
            if(status/=0)exit
        end do
        close(unit,iostat=closed)
        if(status==0.and.closed==0)then
            error=''
        else
            call release_patch_undo(stored)
        end if
    end procedure prepare_patch_snapshot_undo

    module procedure append_patch_undo
        type(patch_edit), allocatable :: grown(:)
        integer :: capacity

        capacity=0
        if(allocated(catalog%undo_history))capacity=size(catalog%undo_history)
        if(catalog%undo_count==capacity)then
            allocate(grown(max(8,2*capacity)))
            if(capacity>0)grown(:capacity)=catalog%undo_history
            call move_alloc(grown,catalog%undo_history)
        end if
        if(edit%stamp==0_int64)edit%stamp=next_patch_revision()
        catalog%undo_count=catalog%undo_count+1
        catalog%undo_history(catalog%undo_count)=edit
        edit=patch_edit()
    end procedure append_patch_undo

    subroutine write_current_state(catalog,targets,names_only,blob,status)
        type(studio_patch_catalog), intent(in) :: catalog
        integer, intent(in) :: targets(:,:)
        logical, intent(in) :: names_only
        type(history_blob), intent(out) :: blob
        integer, intent(out) :: status
        integer :: i,rank,id,unit,closed
        type(patch_entry) :: template

        status=0
        if(size(targets,2)==0)return
        call create_history_blob(blob,unit,status,patch_history_project(catalog))
        if(status/=0)return
        do i=1,size(targets,2)
            rank=targets(1,i)
            id=targets(2,i)
            if(allocated(catalog%catalogs(rank)%entries))then
                associate(entry=>catalog%catalogs(rank)%entries(id))
                    if(entry%removed)then
                        template=patch_entry(dimension=rank,patch_id=id,removed=.true.)
                        call write_history_patch(unit,template,names_only,status)
                        if(status/=0)exit
                        cycle
                    end if
                    if(allocated(entry%geometry))then
                        call write_history_patch(unit,entry,names_only,status)
                        if(status/=0)exit
                        cycle
                    end if
                end associate
            end if
            status=1
            if(.not.read_catalog_patch(rank,id,template,catalog))exit
            call write_history_patch(unit,template,names_only,status)
            if(status/=0)exit
        end do
        close(unit,iostat=closed)
        if(status/=0.or.closed/=0)then
            call release_history_blob(blob)
            status=1
        end if
    end subroutine write_current_state

    subroutine read_history_state(blob,targets,names_only,entries,status)
        type(history_blob), intent(in) :: blob
        integer, intent(in) :: targets(:,:)
        logical, intent(in) :: names_only
        type(patch_entry), allocatable, intent(out) :: entries(:)
        integer, intent(out) :: status
        integer :: i,unit,closed
        integer :: trailing
        character :: extra

        allocate(entries(size(targets,2)),stat=status)
        if(status/=0.or.size(targets,2)==0)return
        call open_history_blob(blob,unit,status)
        if(status/=0)return
        do i=1,size(entries)
            call read_history_patch(unit,entries(i),targets(:,i),names_only,status)
            if(status/=0)exit
        end do
        if(status==0)then
            read(unit,iostat=trailing)extra
            if(trailing/=iostat_end)status=1
        end if
        close(unit,iostat=closed)
        if(closed/=0)status=1
    end subroutine read_history_state

    subroutine write_history_patch(unit,entry,names_only,status)
        integer, intent(in) :: unit
        type(patch_entry), intent(in) :: entry
        logical, intent(in) :: names_only
        integer, intent(out) :: status

        write(unit,iostat=status)entry%dimension,entry%patch_id,entry%revision,entry%mesh_revision,entry%resolution,entry%removed
        if(status/=0.or.entry%removed)return
        status=1
        if(.not.allocated(entry%name))return
        call write_text(unit,entry%name,status)
        if(status/=0.or.names_only)return
        status=1
        if(.not.allocated(entry%geometry))return
        call write_patch(unit,entry%dimension,entry%geometry,status)
    end subroutine write_history_patch

    subroutine read_history_patch(unit,entry,target,names_only,status)
        integer, intent(in) :: unit,target(2)
        logical, intent(in) :: names_only
        type(patch_entry), intent(out) :: entry
        integer, intent(out) :: status

        read(unit,iostat=status)entry%dimension,entry%patch_id,entry%revision,entry%mesh_revision,entry%resolution,entry%removed
        if(status/=0)return
        status=1
        if(any([entry%dimension,entry%patch_id]/=target))return
        status=0
        if(entry%removed)return
        call read_text(unit,entry%name,status)
        if(status/=0.or.names_only)return
        call read_patch(unit,entry%dimension,entry%geometry,status)
        if(status/=0)return
        if(.not.refresh_patch_metadata(entry))status=1
    end subroutine read_history_patch

    subroutine apply_history_state(catalog,counts,selection,entries,names_only)
        type(studio_patch_catalog), intent(inout) :: catalog
        integer, intent(in) :: counts(3),selection(3)
        type(patch_entry), intent(inout) :: entries(:)
        logical, intent(in) :: names_only
        type(patch_entry), allocatable :: resized(:)
        integer :: i,rank,id

        do rank=1,3
            if(counts(rank)==catalog%patch_count(rank))cycle
            allocate(resized(counts(rank)))
            if(allocated(catalog%catalogs(rank)%entries))then
                do id=1,min(counts(rank),size(catalog%catalogs(rank)%entries))
                    call move_patch_entry(catalog%catalogs(rank)%entries(id),resized(id))
                end do
            end if
            call move_alloc(resized,catalog%catalogs(rank)%entries)
        end do
        do i=1,size(entries)
            rank=entries(i)%dimension
            id=entries(i)%patch_id
            if(names_only)then
                call move_alloc(entries(i)%name,catalog%catalogs(rank)%entries(id)%name)
            else
                call move_patch_entry(entries(i),catalog%catalogs(rank)%entries(id))
            end if
        end do
        catalog%selection=selection
    end subroutine apply_history_state

    module procedure move_patch_entry
        integer :: axis

        destination%dimension=source%dimension
        destination%patch_id=source%patch_id
        destination%degree=source%degree
        destination%controls=source%controls
        destination%resolution=source%resolution
        destination%revision=source%revision
        destination%mesh_revision=source%mesh_revision
        destination%removed=source%removed
        destination%bounds=source%bounds
        call move_alloc(source%geometry,destination%geometry)
        call move_alloc(source%name,destination%name)
        do axis=1,size(source%axes)
            call move_alloc(source%axes(axis)%knots,destination%axes(axis)%knots)
        end do
    end procedure move_patch_entry

end submodule studio_single_patch_history
