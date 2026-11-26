module studio_document_deletion
    use studio_document, only: studio_document_t, document_has_patch, document_patch_in_assembly, document_mesh_count, &
        document_mesh_source_group, record_document_mesh_edit, remove_document_mesh, commit_document_deletion
    use studio_multipatch_collection, only: count_multipatches, remove_multipatch
    use studio_single_patch, only: single_patch_name, remove_single_patch, single_patch_count
    use studio_selection, only: selection_target_t,selection_patch,selection_assembly,selection_mesh
    implicit none(type, external)
    private

    public :: can_delete_document_patch, delete_document_patch, delete_document_assembly, delete_document_mesh
    public :: document_deletion_error,delete_document_selection

    type, public :: deleted_mesh_t
        integer :: kind=0, variant=0, source_group=0
        character(:), allocatable :: name
    end type deleted_mesh_t

    type :: deletion_targets
        logical, allocatable :: patches(:,:),groups(:),meshes(:)
    end type deletion_targets

contains

    function document_deletion_error(document,targets) result(error)
        type(studio_document_t), intent(in) :: document
        type(selection_target_t), intent(in) :: targets(:)
        character(:), allocatable :: error
        type(deletion_targets) :: deleted

        call collect_deletion_targets(document,targets,deleted,error)
    end function document_deletion_error

    subroutine collect_deletion_targets(document,targets,deleted,error)
        type(studio_document_t), intent(in) :: document
        type(selection_target_t), intent(in) :: targets(:)
        type(deletion_targets), intent(out) :: deleted
        character(:), allocatable, intent(out) :: error
        integer :: i,kind,id,capacity,group

        error='Select patches, multipatches or meshes to delete.'
        if(size(targets)==0)return
        capacity=maxval([(single_patch_count(i,catalog=document%patches),i=1,3)])
        allocate(deleted%patches(capacity,3),source=.false.)
        allocate(deleted%groups(count_multipatches(document%multipatches)),source=.false.)
        allocate(deleted%meshes(document_mesh_count(document)),source=.false.)
        do i=1,size(targets)
            associate(target=>targets(i))
                error='Select whole patches, multipatches or meshes to delete them.'
                if(target%part/=0.or.target%item/=0.or.target%feature/=0)return
                id=target%id
                kind=target%rank
                select case(target%category)
                case(selection_patch)
                    error='A selected patch is unavailable in this project.'
                    if(.not.document_has_patch(document,kind,id))return
                    error='Built-in templates remain available.'
                    if(id<=2)return
                    deleted%patches(id,kind)=.true.
                case(selection_assembly)
                    error='A selected multipatch is unavailable in this project.'
                    if(id<1.or.id>size(deleted%groups))return
                    deleted%groups(id)=.true.
                case(selection_mesh)
                    error='A selected mesh is unavailable in this project.'
                    if(id<1.or.id>size(deleted%meshes))return
                    deleted%meshes(id)=.true.
                case default
                    return
                end select
            end associate
        end do

        error='Select the containing multipatch too, or detach its member before deleting it.'
        do group=1,size(deleted%groups)
            if(deleted%groups(group))cycle
            associate(definition=>document%multipatches%groups(group)%definition)
                kind=definition%kind
                if(kind<1.or.kind>3.or..not.allocated(definition%members))cycle
                do i=1,size(definition%members)
                    id=definition%members(i)
                    if(id<1.or.id>capacity)cycle
                    if(deleted%patches(id,kind))return
                end do
            end associate
        end do
        error=''
    end subroutine collect_deletion_targets

    subroutine delete_document_selection(document,targets,mesh_only,removed,error)
        type(studio_document_t), intent(inout) :: document
        type(selection_target_t), intent(in) :: targets(:)
        logical, intent(in) :: mesh_only
        integer, intent(out) :: removed
        character(:), allocatable, intent(out) :: error
        type(deletion_targets) :: deleted
        integer, allocatable :: patches(:,:),groups(:),meshes(:)
        integer :: i,n,kind,id,group

        removed=0
        call collect_deletion_targets(document,targets,deleted,error)
        if(len(error)>0)return
        do i=1,size(deleted%meshes)
            if(deleted%meshes(i))cycle
            group=document_mesh_source_group(document,i)
            if(group>0)then
                deleted%meshes(i)=deleted%groups(group)
            else if(.not.allocated(document%meshes(i)%source_group))then
                kind=document%meshes(i)%kind
                id=document%meshes(i)%variant
                if(kind<1.or.kind>3.or.id<1.or.id>size(deleted%patches,1))cycle
                deleted%meshes(i)=deleted%patches(id,kind)
            end if
        end do
        allocate(patches(2,count(deleted%patches)))
        n=0
        do kind=1,3
            do id=1,size(deleted%patches,1)
                if(.not.deleted%patches(id,kind))cycle
                n=n+1
                patches(:,n)=[kind,id]
            end do
        end do
        groups=pack([(i,i=1,size(deleted%groups))],deleted%groups)
        meshes=pack([(i,i=1,size(deleted%meshes))],deleted%meshes)
        call commit_document_deletion(document,patches,groups,meshes,mesh_only,error)
        if(len(error)==0)removed=n+size(groups)+size(meshes)
    end subroutine delete_document_selection

    logical function can_delete_document_patch(document,kind,variant) result(allowed)
        type(studio_document_t), intent(in) :: document
        integer, intent(in) :: kind,variant

        allowed=.false.
        if(variant<=2.or..not.document_has_patch(document,kind,variant))return
        allowed=.not.document_patch_in_assembly(document,kind,variant)
    end function can_delete_document_patch

    subroutine delete_document_patch(document,kind,variant,name,error)
        type(studio_document_t), intent(inout) :: document
        integer, intent(in) :: kind,variant
        character(:), allocatable, intent(out) :: name,error

        error='That patch is unavailable in this project.'
        if(.not.document_has_patch(document,kind,variant))return
        if(.not.can_delete_document_patch(document,kind,variant))then
            if(variant<=2)then
                error='Built-in templates remain available.'
            else
                error='Detach this patch from its assembly before deleting it.'
            end if
            return
        end if
        name=single_patch_name(kind,variant,catalog=document%patches)
        call remove_single_patch(kind,variant,error,catalog=document%patches)
    end subroutine delete_document_patch

    subroutine delete_document_assembly(document,id,error)
        type(studio_document_t), intent(inout) :: document
        integer, intent(in) :: id
        character(:), allocatable, intent(out) :: error

        error='That assembly is unavailable in this project.'
        if(id<1.or.id>count_multipatches(document%multipatches))return
        call remove_multipatch(document%multipatches,id,catalog=document%patches,error=error)
    end subroutine delete_document_assembly

    subroutine delete_document_mesh(document,id,selected_mesh,mesh_only,removed,error)
        type(studio_document_t), intent(inout) :: document
        integer, intent(in) :: id,selected_mesh
        logical, intent(in) :: mesh_only
        type(deleted_mesh_t), intent(out) :: removed
        character(:), allocatable, intent(out) :: error
        logical :: saved

        error='That mesh is unavailable in this project.'
        if(id<1.or.id>document_mesh_count(document))return
        removed%name=document%meshes(id)%name
        removed%kind=document%meshes(id)%kind
        removed%variant=document%meshes(id)%variant
        removed%source_group=document_mesh_source_group(document,id)
        call record_document_mesh_edit(document,[id],0,selected_mesh,mesh_only,saved)
        error='Could not write the mesh history.'
        if(.not.saved)return
        call remove_document_mesh(document,id)
        document%mesh_generation=document%mesh_generation+1
        error=''
    end subroutine delete_document_mesh

end module studio_document_deletion
