submodule (studio_single_patch) studio_single_patch_catalog
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use studio_forcad, only: describe_patch, copy_native_patch, translate_patch
    use studio_sample_geometry, only: sample_name, make_sample_patch
    implicit none(type, external)

contains

    module procedure initialize_studio_patch_catalog

        if(catalog%template_revision==0_int64)catalog%template_revision=next_patch_revision()
    end procedure initialize_studio_patch_catalog

    module procedure allocate_studio_patch_catalog
        integer :: kind,variant

        staged%selection=selection
        do kind=1,3
            allocate(staged%catalogs(kind)%entries(counts(kind)),stat=status)
            if(status/=0)return
            do variant=1,counts(kind)
                staged%catalogs(kind)%entries(variant)=patch_entry(patch_id=variant,removed=.true.)
            end do
        end do
    end procedure allocate_studio_patch_catalog

    module procedure move_studio_patch_catalog
        integer :: kind

        do kind=1,size(source%catalogs)
            call move_alloc(source%catalogs(kind)%entries,destination%catalogs(kind)%entries)
        end do
        call move_alloc(source%undo_history,destination%undo_history)
        call move_alloc(source%history_project,destination%history_project)
        destination%selection=source%selection
        destination%undo_count=source%undo_count
        destination%undo_floor=source%undo_floor
        destination%template_revision=source%template_revision
        source%selection=builtin_samples
        source%undo_count=0
        source%undo_floor=0_int64
        source%template_revision=0_int64
    end procedure move_studio_patch_catalog

    module procedure stage_studio_patch
        type(patch_entry) :: candidate

        error='The project contains invalid geometry for '//name//'.'
        candidate%patch_id=variant
        candidate%dimension=kind
        candidate%name=name
        candidate%resolution=resolution
        candidate%resolution(:kind)=min(resolution(:kind),merge(curve_samples,tensor_samples,kind==1))
        call move_alloc(geometry,candidate%geometry)
        if(.not.refresh_patch_metadata(candidate))return
        if(.not.can_sample_patch(candidate,sampled))return
        call move_patch_entry(candidate,staged%catalogs(kind)%entries(variant))
        error=''
    end procedure stage_studio_patch

    module procedure studio_patch_is_staged

        available=.false.
        if(kind<1.or.kind>3.or.variant<1)return
        if(variant>size(staged%catalogs(kind)%entries))return
        available=.not.staged%catalogs(kind)%entries(variant)%removed
    end procedure studio_patch_is_staged

    module procedure install_studio_patch_catalog
        type(studio_patch_catalog), pointer :: owner
        integer :: kind,variant

        owner=>edit_catalog(catalog)
        ! A pending edit must not match a template in a different document.
        owner%template_revision=next_patch_revision()
        do kind=1,3
            call move_alloc(staged%catalogs(kind)%entries,owner%catalogs(kind)%entries)
            if(.not.allocated(owner%catalogs(kind)%entries))cycle
            do variant=1,size(owner%catalogs(kind)%entries)
                if(owner%catalogs(kind)%entries(variant)%removed)cycle
                call advance_patch_revision(owner%catalogs(kind)%entries(variant))
            end do
        end do
        owner%selection=staged%selection
        call clear_patch_history(owner)
        owner%undo_floor=revision_clock
    end procedure install_studio_patch_catalog

    module procedure studio_patch_resolution

        if(present(catalog))then
            resolution=catalog_patch_resolution(catalog,kind,variant)
        else
            resolution=catalog_patch_resolution(default_catalog,kind,variant)
        end if
    end procedure studio_patch_resolution

    function catalog_patch_resolution(catalog,kind,variant) result(resolution)
        type(studio_patch_catalog), intent(in) :: catalog
        integer, intent(in) :: kind,variant
        integer :: resolution(3)

        resolution=0
        if(.not.catalog%patch_available(kind,variant))return
        if(allocated(catalog%catalogs(kind)%entries))resolution=catalog%catalogs(kind)%entries(variant)%resolution
    end function catalog_patch_resolution

    module procedure set_active_patch
        type(studio_patch_catalog), pointer :: owner

        owner=>edit_catalog(catalog)
        if (.not.owner%patch_available(kind,patch_id)) return
        owner%selection(kind)=patch_id
    end procedure set_active_patch

    module procedure active_patch_id

        patch_id=0
        if (kind<1.or.kind>size(default_catalog%selection)) return
        if(present(catalog))then
            patch_id=catalog%selection(kind)
        else
            patch_id=default_catalog%selection(kind)
        end if
    end procedure active_patch_id

    module procedure single_patch_count

        if(present(catalog))then
            number=catalog%patch_count(kind)
        else
            number=default_catalog%patch_count(kind)
        end if
    end procedure single_patch_count

    module procedure catalog_patch_count

        number=0
        if (kind<1.or.kind>size(catalog%catalogs)) return
        number=builtin_samples
        if (allocated(catalog%catalogs(kind)%entries)) number=size(catalog%catalogs(kind)%entries)
    end procedure catalog_patch_count

    module procedure single_patch_available

        if(present(catalog))then
            available=catalog%patch_available(kind,variant)
        else
            available=default_catalog%patch_available(kind,variant)
        end if
    end procedure single_patch_available

    module procedure catalog_patch_available

        available=.false.
        if(variant<1.or.variant>catalog%patch_count(kind))return
        available=.true.
        if(allocated(catalog%catalogs(kind)%entries))available=.not.catalog%catalogs(kind)%entries(variant)%removed
    end procedure catalog_patch_available

    module procedure single_patch_revision


        if(present(catalog))then
            revision=catalog%patch_revision(kind,variant)
        else
            revision=default_catalog%patch_revision(kind,variant)
        end if
    end procedure single_patch_revision

    module procedure single_patch_mesh_revision

        if(present(catalog))then
            revision=catalog%patch_revision(kind,variant,mesh=.true.)
        else
            revision=default_catalog%patch_revision(kind,variant,mesh=.true.)
        end if
    end procedure single_patch_mesh_revision

    module procedure catalog_patch_revision

        revision=-1_int64
        if(.not.catalog%patch_available(kind,variant))return
        revision=0_int64
        if(allocated(catalog%catalogs(kind)%entries))then
            associate(patch=>catalog%catalogs(kind)%entries(variant))
                revision=patch%revision
                if(present(mesh))then
                    if(mesh)revision=patch%mesh_revision
                end if
            end associate
        end if
        if(revision==0_int64)revision=catalog%template_revision
    end procedure catalog_patch_revision

    module procedure single_patch_name

        if(present(catalog))then
            name=catalog%patch_name(kind,variant)
        else
            name=default_catalog%patch_name(kind,variant)
        end if
    end procedure single_patch_name

    module procedure catalog_patch_name

        name=''
        if (.not.catalog%patch_available(kind,variant)) return
        if (allocated(catalog%catalogs(kind)%entries))then
            if(allocated(catalog%catalogs(kind)%entries(variant)%name))then
                name=catalog%catalogs(kind)%entries(variant)%name
                return
            end if
        end if
        name=sample_name(kind,variant)
    end procedure catalog_patch_name

    module procedure rename_single_patch
        type(studio_patch_catalog), pointer :: owner
        character(:), allocatable :: proposed
        type(patch_entry) :: previous
        integer :: rank,index

        proposed=trim(adjustl(name))
        error='Enter a patch name.'
        if(len(proposed)==0)return
        error='The patch is no longer available.'
        owner=>edit_catalog(catalog)
        if(.not.read_catalog_patch(kind,variant,previous,owner))return
        if(proposed==previous%name)then
            error=''
            return
        end if
        do rank=1,size(owner%catalogs)
            do index=1,owner%patch_count(rank)
                if(rank==kind.and.index==variant)cycle
                if(owner%patch_name(rank,index)/=proposed)cycle
                error='A patch already uses this name.'
                return
            end do
        end do
        call store_patch_template(previous,owner)
        call record_undo(patch_edit(kind=kind,operation='rename',previous=previous),owner,error)
        if(len(error)>0)return
        owner%catalogs(kind)%entries(variant)%name=proposed
        error=''
    end procedure rename_single_patch

    module procedure remove_single_patch
        type(studio_patch_catalog), pointer :: owner
        type(patch_entry) :: previous

        error='Built-in templates remain in the project.'
        if(variant<=builtin_samples)return
        error='The patch is no longer available.'
        owner=>edit_catalog(catalog)
        if(.not.read_catalog_patch(kind,variant,previous,owner))return
        call record_undo(patch_edit(kind=kind,operation='patch removal',previous=previous),owner,error)
        if(len(error)>0)return
        ! Keep patch IDs stable for scene entries, multipatches and earlier edits.
        owner%catalogs(kind)%entries(variant)=patch_entry(patch_id=variant,removed=.true.)
        if(owner%selection(kind)==variant)owner%selection(kind)=builtin_samples
        error=''
    end procedure remove_single_patch

    module procedure remove_single_patches
        type(patch_edit) :: edit,stored
        integer :: i,rank,id

        error=''
        if(size(targets,2)==0)return
        error='A patch is no longer available for deletion.'
        if(size(targets,1)/=2)return
        do i=1,size(targets,2)
            rank=targets(1,i)
            id=targets(2,i)
            if(id<=builtin_samples.or..not.single_patch_available(rank,id,catalog))return
        end do
        edit%kind=targets(1,1)
        edit%operation='patch removal'
        edit%batch_counts=[(catalog%patch_count(rank),rank=1,3)]
        edit%batch_selection=catalog%selection
        call prepare_patch_set_undo(edit,catalog,targets,stored,error)
        if(len(error)>0)return
        call append_patch_undo(catalog,stored)
        do i=1,size(targets,2)
            rank=targets(1,i)
            id=targets(2,i)
            catalog%catalogs(rank)%entries(id)=patch_entry(patch_id=id,removed=.true.)
            if(catalog%selection(rank)==id)catalog%selection(rank)=builtin_samples
        end do
    end procedure remove_single_patches

    module procedure add_single_patch
        type(studio_patch_catalog), pointer :: owner
        type(patch_entry) :: candidate
        character(:), allocatable :: patch_name
        integer :: rank,index,new_variant
        real(rk) :: translation(3)

        variant=0
        patch_name=trim(adjustl(name))
        error='Enter a patch name.'
        if (len(patch_name)==0) return
        owner=>edit_catalog(catalog)
        error='Choose a valid source patch.'
        if (.not.owner%patch_available(kind,source_variant)) return
        do rank=1,size(owner%catalogs)
            do index=1,owner%patch_count(rank)
                if (owner%patch_name(rank,index)/=patch_name) cycle
                error='A patch already uses this name.'
                return
            end do
        end do
        error='ForCAD could not load the source patch.'
        if (.not.read_catalog_patch(kind,source_variant,candidate,owner)) return
        if(present(offset))then
            error='Enter finite lengths for the X, Y and Z placement offsets.'
            if(.not.all(ieee_is_finite(offset)))return
            translation=real(offset,rk)
            if(.not.all(ieee_is_finite(translation)))return
            if(any(translation/=0._rk))then
                call translate_patch(candidate%geometry,translation)
            end if
        end if
        new_variant=owner%patch_count(kind)+1
        candidate%patch_id=new_variant
        candidate%name=patch_name
        error='ForCAD could not sample the new patch.'
        if (.not.can_sample_patch(candidate)) return

        call record_undo(patch_edit(kind=kind,created_patch_id=new_variant, &
            previous_selection=owner%selection(kind),operation='patch creation'),owner,error)
        if(len(error)>0)return
        call advance_patch_revision(candidate)
        call grow_patch_catalog(owner%catalogs(kind),new_variant)
        call move_patch_entry(candidate,owner%catalogs(kind)%entries(new_variant))
        owner%selection(kind)=new_variant
        variant=new_variant
        error=''
    end procedure add_single_patch

    module procedure publish_native_patches
        type(patch_entry), allocatable :: batch(:)
        type(studio_patch_catalog), pointer :: owner
        type(patch_edit) :: edit
        integer :: i,j,rank,counts(3),selection(3)
        logical :: store_undo

        kind=0; variant=0
        error='No native patches were supplied.'
        if(size(geometries)==0.or.size(geometries)/=size(names))return
        if(present(resolutions))then
            error='The display sample counts do not match the imported geometry.'
            if(size(resolutions,1)/=3.or.size(resolutions,2)/=size(geometries))return
        end if
        owner=>edit_catalog(catalog)
        counts=[(owner%patch_count(i),i=1,3)]
        selection=owner%selection
        edit%batch_counts=counts
        edit%batch_selection=selection
        allocate(batch(size(geometries)))
        do i=1,size(geometries)
            error='Missing native geometry for '//trim(names(i))//'.'
            if(.not.allocated(geometries(i)%geometry))return
            rank=native_patch_kind(geometries(i)%geometry)
            if(rank==0)return
            counts(rank)=counts(rank)+1
            selection(rank)=counts(rank)
            batch(i)%dimension=rank
            batch(i)%patch_id=counts(rank)
            batch(i)%name=available_patch_name(trim(names(i)),batch(:i-1),owner)
            if(present(resolutions))batch(i)%resolution(:rank)=resolutions(:rank,i)
            call copy_native_patch(geometries(i)%geometry,batch(i)%geometry)
            error='ForCAD could not validate '//batch(i)%name//'.'
            if(.not.refresh_patch_metadata(batch(i)))return
            if(.not.can_sample_patch(batch(i)))return
        end do

        store_undo=.true.
        if(present(record))store_undo=record
        kind=batch(size(batch))%dimension
        variant=batch(size(batch))%patch_id
        if(store_undo)then
            edit%kind=kind
            edit%operation='geometry import'
            if(present(operation))edit%operation=operation
            call record_undo(edit,owner,error)
            if(len(error)>0)then
                kind=0
                variant=0
                return
            end if
        end if
        do rank=1,3
            if(counts(rank)==edit%batch_counts(rank))cycle
            call grow_patch_catalog(owner%catalogs(rank),counts(rank))
            do j=1,size(batch)
                if(batch(j)%dimension/=rank)cycle
                call advance_patch_revision(batch(j))
                call move_patch_entry(batch(j),owner%catalogs(rank)%entries(batch(j)%patch_id))
            end do
        end do
        owner%selection=selection
        error=''
    end procedure publish_native_patches

    subroutine grow_patch_catalog(group,count)
        type(patch_catalog), intent(inout) :: group
        integer, intent(in) :: count
        type(patch_entry), allocatable :: grown(:)
        integer :: i

        allocate(grown(count))
        if(allocated(group%entries))then
            do i=1,size(group%entries)
                call move_patch_entry(group%entries(i),grown(i))
            end do
        end if
        call move_alloc(grown,group%entries)
    end subroutine grow_patch_catalog

    module procedure available_patch_name
        character(:), allocatable :: base
        character(16) :: suffix
        integer :: serial,rank,index

        base=requested
        if(len(base)==0)base='Imported patch'
        serial=1
        candidate: do
            name=base
            if(serial>1)then
                write(suffix,'(i0)')serial
                name=base//' ('//trim(suffix)//')'
            end if
            serial=serial+1
            do rank=1,3
                do index=1,single_patch_count(rank,catalog)
                    if(name==single_patch_name(rank,index,catalog))cycle candidate
                end do
            end do
            do index=1,size(batch)
                if(name==batch(index)%name)cycle candidate
            end do
            return
        end do candidate
    end procedure available_patch_name

    module procedure resolve_patch_snapshot
        integer::index

        index=active_patch_id(kind,catalog)
        if(present(patch_id))index=patch_id
        valid=copy_patch_snapshot(kind,index,patch,catalog)
    end procedure resolve_patch_snapshot

    module procedure copy_patch_snapshot

        if (present(catalog)) then
            valid=read_catalog_patch(kind,variant,patch,catalog)
        else
            valid=read_catalog_patch(kind,variant,patch,default_catalog)
        end if
    end procedure copy_patch_snapshot

    module procedure read_catalog_patch

        valid=.false.
        if (.not.catalog%patch_available(kind,variant)) return
        if (allocated(catalog%catalogs(kind)%entries)) then
            if (allocated(catalog%catalogs(kind)%entries(variant)%geometry)) then
                patch=catalog%catalogs(kind)%entries(variant)
                valid=.true.
                return
            end if
        end if
        valid=make_template_patch(kind,variant,patch)
    end procedure read_catalog_patch

    module procedure store_patch_template

        associate(rank_catalog=>catalog%catalogs(patch%dimension))
            if (.not.allocated(rank_catalog%entries)) allocate(rank_catalog%entries(builtin_samples))
            if (.not.allocated(rank_catalog%entries(patch%patch_id)%geometry)) rank_catalog%entries(patch%patch_id)=patch
        end associate
    end procedure store_patch_template

    module procedure edit_catalog

        if(present(catalog))then
            owner=>catalog
        else
            owner=>default_catalog
        end if
    end procedure edit_catalog

    module procedure copy_single_patch_geometry
        type(patch_entry) :: patch

        if(.not.copy_patch_snapshot(kind,variant,patch,catalog))return
        call move_alloc(patch%geometry,geometry)
    end procedure copy_single_patch_geometry

    logical function make_template_patch(kind,sample_index,patch) result(valid)
        integer, intent(in) :: kind,sample_index
        type(patch_entry), intent(out) :: patch

        patch%dimension=kind
        patch%patch_id=sample_index
        patch%name=sample_name(kind,sample_index)
        call make_sample_patch(kind,sample_index,patch%geometry)
        valid=refresh_patch_metadata(patch)
    end function make_template_patch

    module procedure refresh_patch_metadata

        valid=.false.
        if (.not.allocated(patch%geometry)) return
        valid=describe_patch(patch%geometry,patch%dimension,patch%degree,patch%controls,patch%axes,patch%bounds)
    end procedure refresh_patch_metadata

end submodule studio_single_patch_catalog
