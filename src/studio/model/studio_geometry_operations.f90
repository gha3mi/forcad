module studio_geometry_operations
    use, intrinsic :: iso_c_binding, only: c_double
    use studio_forcad, only: rk, affine_transform_t, extract_native_slice, studio_native_patch
    use studio_document, only: studio_document_t, document_has_patch, document_patch_in_assembly, document_mesh_count, &
        record_document_properties
    use studio_single_patch, only: commit_patch_transform, single_patch_name, &
        rename_single_patch, commit_single_patch_split, commit_single_patch_extract, &
        copy_single_patch_geometry, publish_native_patches
    use studio_multipatch_collection, only: count_multipatches, refresh_multipatch
    use studio_scene, only: refresh_studio_scene, studio_scene_is_visible, studio_scene_set_visible, &
        get_studio_material, set_studio_material
    use studio_presentation_text, only: number_text
    implicit none(type, external)
    private

    public :: read_geometry_target, commit_document_transform, commit_patch_range, commit_patch_section, rename_document_patch
    public :: commit_document_visibility

contains

    subroutine read_geometry_target(document,kind,patch_id,assembly_id,patch_kind,patch_ids,name)
        type(studio_document_t), intent(in) :: document
        integer, intent(in) :: kind,patch_id,assembly_id
        integer, intent(out) :: patch_kind
        integer, allocatable, intent(out) :: patch_ids(:)
        character(:), allocatable, intent(out) :: name

        name=''
        patch_kind=kind
        allocate(patch_ids(0))
        if(assembly_id>0)then
            if(assembly_id>count_multipatches(document%multipatches))return
            associate(definition=>document%multipatches%groups(assembly_id)%definition)
                patch_kind=definition%kind
                patch_ids=definition%members
                name=definition%name
            end associate
            return
        end if
        if(assembly_id<0.or.kind<1.or.kind>3.or.patch_id<1)return
        if(.not.document_has_patch(document,kind,patch_id))return
        patch_ids=[patch_id]
        name=document%patches%patch_name(kind,patch_id)
    end subroutine read_geometry_target

    subroutine commit_document_transform(document,patch_kind,patch_ids,assembly_id,transform,changed,error)
        type(studio_document_t), intent(inout) :: document
        integer, intent(in) :: patch_kind,patch_ids(:),assembly_id
        type(affine_transform_t), intent(in) :: transform
        logical, intent(out) :: changed
        character(:), allocatable, intent(out) :: error

        changed=.false.
        error='Select a curve, surface, volume or assembly to transform.'
        if(size(patch_ids)==0.or.assembly_id<0)return
        if(assembly_id>0)then
            if(assembly_id>count_multipatches(document%multipatches))return
            call refresh_multipatch(document%multipatches,assembly_id,.false.,catalog=document%patches)
            error='Validate the multipatch and resolve its interfaces before transforming it.'
            associate(assembly=>document%multipatches%groups(assembly_id))
                if(.not.allocated(assembly%native))return
                error='The transform targets do not match the assembly members.'
                if(patch_kind/=assembly%definition%kind.or.size(patch_ids)/=size(assembly%definition%members))return
                if(any(patch_ids/=assembly%definition%members))return
                call commit_patch_transform(patch_kind,patch_ids,transform,changed,error, &
                    assembly%native,multipatch_id=assembly_id,catalog=document%patches)
            end associate
        else
            call commit_patch_transform(patch_kind,patch_ids,transform,changed,error,catalog=document%patches)
        end if
    end subroutine commit_document_transform

    subroutine commit_patch_range(document,kind,source,axis,lower,upper,splitting,created,error)
        type(studio_document_t), intent(inout) :: document
        integer, intent(in) :: kind,source,axis
        real(c_double), intent(in) :: lower,upper
        logical, intent(in) :: splitting
        integer, intent(out) :: created
        character(:), allocatable, intent(out) :: error
        integer :: rgb(3),opacity
        logical :: source_visible

        created=0
        error='Select a patch and one of its parameter directions.'
        if(axis<1.or.axis>kind.or..not.document_has_patch(document,kind,source))return
        if(splitting)then
            error='Detach this patch from its multipatch before splitting it.'
            if(document_patch_in_assembly(document,kind,source))return
        end if
        source_visible=studio_scene_is_visible(document%scene,kind,source,catalog=document%patches)
        call get_studio_material(document%scene,kind,source,rgb,opacity)
        if(splitting)then
            call commit_single_patch_split(kind,axis,lower,created,error,catalog=document%patches,patch_id=source)
        else
            call commit_single_patch_extract(kind,axis,lower,upper,created,error,catalog=document%patches,patch_id=source)
        end if
        if(len(error)>0.or.created==0)return
        call refresh_studio_scene(document%scene,catalog=document%patches)
        call set_studio_material(document%scene,kind,created,rgb,opacity)
        if(splitting)call studio_scene_set_visible(document%scene,kind,created,source_visible,catalog=document%patches)
    end subroutine commit_patch_range

    subroutine commit_patch_section(document,kind,source,axis,parameter,created,error)
        type(studio_document_t), intent(inout) :: document
        integer, intent(in) :: kind,source,axis
        real(c_double), intent(in) :: parameter
        integer, intent(out) :: created
        character(:), allocatable, intent(out) :: error
        type(studio_native_patch) :: sections(1)
        class(*), allocatable :: geometry
        character(*), parameter :: axes='UVW'
        character(:), allocatable :: name
        integer :: result_kind

        created=0
        error='Select a surface or volume and one of its parameter directions.'
        if(kind<2.or.kind>3.or.axis<1.or.axis>kind)return
        if(.not.document_has_patch(document,kind,source))return
        call copy_single_patch_geometry(kind,source,geometry,catalog=document%patches)
        if(.not.allocated(geometry))return
        call extract_native_slice(geometry,axis,real(parameter,rk),sections(1)%geometry,error)
        if(len(error)>0)return
        name=single_patch_name(kind,source,catalog=document%patches)//' · '//axes(axis:axis)//' = '//number_text(parameter)
        call publish_native_patches(sections,[name],result_kind,created,error,operation='section extraction',catalog=document%patches)
        if(len(error)>0.or.created==0)return
        call refresh_studio_scene(document%scene,catalog=document%patches)
    end subroutine commit_patch_section

    subroutine rename_document_patch(document,kind,variant,name,changed,error)
        type(studio_document_t), intent(inout) :: document
        integer, intent(in) :: kind,variant
        character(*), intent(in) :: name
        logical, intent(out) :: changed
        character(:), allocatable, intent(out) :: error

        changed=.false.
        error='That patch is unavailable in this project.'
        if(.not.document_has_patch(document,kind,variant))return
        error=''
        if(trim(adjustl(name))==single_patch_name(kind,variant,catalog=document%patches))return
        call rename_single_patch(kind,variant,name,error,catalog=document%patches)
        changed=len(error)==0
    end subroutine rename_document_patch

    subroutine commit_document_visibility(document,visible,changed,error,patches,meshes)
        type(studio_document_t), intent(inout) :: document
        logical, intent(in) :: visible
        logical, intent(out) :: changed
        character(:), allocatable, intent(out) :: error
        integer, intent(in), optional :: patches(:,:),meshes(:)
        integer, allocatable :: patch_changes(:,:),mesh_changes(:)
        integer :: i,kind,id,npatches,nmeshes

        changed=.false.
        error=''
        npatches=0
        nmeshes=0
        if(present(patches))then
            allocate(patch_changes(2,size(patches,2)))
            do i=1,size(patches,2)
                kind=patches(1,i)
                id=patches(2,i)
                if(.not.document_has_patch(document,kind,id))cycle
                if(studio_scene_is_visible(document%scene,kind,id,catalog=document%patches).eqv.visible)cycle
                npatches=npatches+1
                patch_changes(:,npatches)=[kind,id]
            end do
        else
            allocate(patch_changes(2,0))
        end if
        if(present(meshes))then
            allocate(mesh_changes(size(meshes)))
            do i=1,size(meshes)
                id=meshes(i)
                if(id<1.or.id>document_mesh_count(document))cycle
                if(document%meshes(id)%visible.eqv.visible)cycle
                nmeshes=nmeshes+1
                mesh_changes(nmeshes)=id
            end do
        else
            allocate(mesh_changes(0))
        end if
        if(npatches+nmeshes==0)return

        if(.not.record_document_properties(document,'Visibility',patches=patch_changes(:,:npatches),meshes=mesh_changes(:nmeshes)))then
            error='Could not record visibility in the edit history.'
            return
        end if
        if(npatches>0)then
            call refresh_studio_scene(document%scene,catalog=document%patches)
            do i=1,npatches
                kind=patch_changes(1,i)
                id=patch_changes(2,i)
                document%scene%catalog(kind)%entries(id)%visible=visible
            end do
            call refresh_studio_scene(document%scene,catalog=document%patches)
        end if
        do i=1,nmeshes
            document%meshes(mesh_changes(i))%visible=visible
        end do
        changed=.true.
    end subroutine commit_document_visibility
end module studio_geometry_operations
