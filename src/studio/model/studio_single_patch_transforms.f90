submodule (studio_single_patch) studio_single_patch_transforms
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use studio_forcad, only: rotation_matrix, transform_patch, transform_assembly
    use studio_units, only: valid_unit_id, model_length
    implicit none(type, external)

contains

    module procedure commit_patch_transform
        type(patch_entry), allocatable :: candidates(:)
        type(studio_native_patch), allocatable :: members(:)
        type(affine_transform_t), parameter :: identity=affine_transform_t()
        character(:), allocatable :: operation
        integer :: i,variant

        changed=.false.
        error='Select geometry to transform.'
        if(kind<1.or.kind>3.or.size(variants)==0)return
        error='Enter finite transform values.'
        if(.not.all(ieee_is_finite(transform%linear)).or..not.all(ieee_is_finite(transform%offset)).or. &
            .not.all(ieee_is_finite(transform%pivot)))return
        if(all(transform%linear==identity%linear).and.all(transform%offset==0._c_double))then
            error=''
            return
        end if
        allocate(candidates(size(variants)))
        do i=1,size(variants)
            variant=variants(i)
            error='An object in this transform is unavailable.'
            if(any(variants(:i-1)==variant))return
            if(.not.copy_patch_snapshot(kind,variant,candidates(i),catalog))return
        end do
        error='ForCAD could not apply the transformation.'
        if(present(native_model))then
            if(any(candidates%dimension/=kind))return
            allocate(members(size(candidates)))
            if(.not.transform_assembly(native_model,kind,transform,members))return
            do i=1,size(candidates)
                call move_alloc(members(i)%geometry,candidates(i)%geometry)
            end do
            operation='multipatch transform'
        else
            do i=1,size(candidates)
                if(.not.transform_patch(candidates(i)%geometry,transform))return
            end do
            operation='transform'
        end if
        call publish_patch_set_edit(kind,variants,candidates,operation,changed,error,multipatch_id,catalog)
    end procedure commit_patch_transform

    module procedure commit_native_patch_edits
        type(patch_entry), allocatable :: candidates(:)
        integer :: i

        changed=.false.
        error='The staged edit does not match the member patches.'
        if(kind<1.or.kind>3.or.size(variants)==0.or.size(patches)/=size(variants))return
        allocate(candidates(size(variants)))
        do i=1,size(variants)
            if(any(variants(:i-1)==variants(i)))return
            if(.not.allocated(patches(i)%geometry))return
            if(.not.copy_patch_snapshot(kind,variants(i),candidates(i),catalog))return
            call move_alloc(patches(i)%geometry,candidates(i)%geometry)
            error='ForCAD could not read the edited patch parameters.'
            if(.not.refresh_patch_metadata(candidates(i)))return
        end do
        call publish_patch_set_edit(kind,variants,candidates,operation,changed,error,multipatch_id,catalog,preserve_mesh)
    end procedure commit_native_patch_edits

    subroutine publish_patch_set_edit(kind,variants,candidates,operation,changed,error,multipatch_id,catalog,preserve_mesh)
        integer, intent(in) :: kind,variants(:)
        type(patch_entry), intent(inout) :: candidates(:)
        character(*), intent(in) :: operation
        logical, intent(out) :: changed
        character(:), allocatable, intent(out) :: error
        integer, intent(in), optional :: multipatch_id
        type(studio_patch_catalog), intent(inout), target, optional :: catalog
        type(studio_patch_catalog), pointer :: owner
        logical, intent(in), optional :: preserve_mesh
        type(patch_edit) :: edit,stored
        integer :: targets(2,size(variants)),i,variant
        integer(int64) :: mesh_revision
        logical :: keep_lattice

        changed=.false.
        error='ForCAD could not sample the edited geometry.'
        do i=1,size(candidates)
            if(.not.can_sample_patch(candidates(i)))return
            targets(:,i)=[kind,variants(i)]
        end do
        owner=>edit_catalog(catalog)
        edit%kind=kind
        edit%operation=operation
        if(present(multipatch_id))edit%multipatch_id=multipatch_id
        edit%previous_selection=owner%selection(kind)
        call prepare_patch_set_undo(edit,owner,targets,stored,error)
        if(len(error)>0)return
        call append_patch_undo(owner,stored)
        if(.not.allocated(owner%catalogs(kind)%entries))allocate(owner%catalogs(kind)%entries(builtin_samples))
        keep_lattice=.false.
        if(present(preserve_mesh))keep_lattice=preserve_mesh
        do i=1,size(candidates)
            variant=variants(i)
            if(keep_lattice)then
                mesh_revision=single_patch_mesh_revision(kind,variant,owner)
                call advance_patch_revision(candidates(i),mesh_revision)
            else
                call advance_patch_revision(candidates(i))
            end if
            call move_patch_entry(candidates(i),owner%catalogs(kind)%entries(variant))
        end do
        changed=.true.
        error=''
    end subroutine publish_patch_set_edit

    module procedure commit_scene_transforms
        ! Rows identify each distinct catalog entry as (kind, variant).
        type(patch_entry), allocatable :: candidates(:)
        type(studio_patch_catalog), pointer :: owner
        type(patch_edit) :: edit
        type(patch_edit), allocatable :: stored(:)
        logical, allocatable :: seen(:,:)
        type(affine_transform_t), allocatable :: transforms(:)
        real(c_double) :: angles(3)
        logical, allocatable :: active(:)
        integer :: i,j,k,axis,kind,variant,total
        character(12) :: step_number

        changed=.false.
        error='There is no geometry to transform.'
        total=size(targets,2)
        if(size(targets,1)/=2.or.total==0)return
        error='A scene transform needs at least one step and valid project units.'
        if(size(steps)==0.or..not.valid_unit_id(model_unit))return
        allocate(transforms(size(steps)),active(size(steps)),stored(size(steps)))
        active=.false.
        do j=1,size(steps)
            associate(step=>steps(j))
                axis=step%axis
                error='Every transform needs a finite value and an X, Y or Z axis.'
                if(axis<1.or.axis>3.or..not.ieee_is_finite(step%value))return
                select case(step%op)
                case(scene_rotate)
                    error='Specify rotation angles in degrees.'
                    if(step%unit/=unitless)return
                    angles=0._c_double
                    angles(axis)=modulo(step%value,360._c_double)
                    active(j)=angles(axis)/=0._c_double
                    transforms(j)%linear=rotation_matrix(angles)
                case(scene_translate)
                    error='The move has an invalid length or unit.'
                    if(.not.valid_unit_id(step%unit))return
                    transforms(j)%offset(axis)=model_length(step%value,model_unit,step%unit)
                    if(.not.ieee_is_finite(transforms(j)%offset(axis)))return
                    active(j)=transforms(j)%offset(axis)/=0._c_double
                case default
                    error='The scene transform contains an unknown operation.'
                    return
                end select
            end associate
        end do
        if(.not.any(active))then
            error=''
            return
        end if

        allocate(candidates(total))
        allocate(seen(maxval([(single_patch_count(kind,catalog),kind=1,3)]),3),source=.false.)
        do i=1,total
            kind=targets(1,i)
            variant=targets(2,i)
            error='An object in this transform no longer exists.'
            if(.not.single_patch_available(kind,variant,catalog))return
            error='An object was listed more than once in this transform.'
            if(seen(variant,kind))return
            seen(variant,kind)=.true.
            error='ForCAD could not load an object for the transform.'
            if(.not.copy_patch_snapshot(kind,variant,candidates(i),catalog))return
        end do

        do j=1,size(steps)
            if(.not.active(j))cycle
            edit%kind=targets(1,1)
            edit%operation='scene transform'
            if(steps(j)%op==scene_rotate)edit%operation='scene rotation'
            call prepare_patch_snapshot_undo(edit,candidates,targets,stored(j),error,catalog)
            if(len(error)>0)exit
            do i=1,total
                if(.not.transform_patch(candidates(i)%geometry,transforms(j)))then
                    write(step_number,'(i0)')j
                    error='ForCAD could not apply step '//trim(step_number)//' to '//candidates(i)%name//'.'
                    exit
                end if
                if(.not.can_sample_patch(candidates(i)))then
                    write(step_number,'(i0)')j
                    error='ForCAD could not sample '//candidates(i)%name//' after step '//trim(step_number)//'.'
                    exit
                end if
                call advance_patch_revision(candidates(i))
            end do
            if(len(error)>0)exit
        end do
        if(len(error)>0)then
            do k=1,j
                call release_patch_undo(stored(k))
            end do
            return
        end if

        owner=>edit_catalog(catalog)
        do j=1,size(stored)
            if(active(j))call append_patch_undo(owner,stored(j))
        end do
        do i=1,total
            kind=targets(1,i)
            variant=targets(2,i)
            if(.not.allocated(owner%catalogs(kind)%entries))allocate(owner%catalogs(kind)%entries(builtin_samples))
            call move_patch_entry(candidates(i),owner%catalogs(kind)%entries(variant))
        end do
        changed=.true.
        error=''
    end procedure commit_scene_transforms

end submodule studio_single_patch_transforms
