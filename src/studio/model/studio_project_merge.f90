module studio_project_merge
    use, intrinsic :: iso_c_binding, only: c_double
    use studio_forcad, only: rk, affine_transform_t, transform_patch
    use studio_document, only: studio_document_t, document_has_patch, document_mesh_count, document_mesh_valid, &
        document_mesh_source_group, append_document_mesh, record_document_merge
    use studio_mesh, only: studio_mesh_t
    use studio_multipatch, only: multipatch_draft
    use studio_multipatch_collection, only: multipatch_collection, count_multipatches, restore_multipatch_definitions
    use studio_scene, only: refresh_studio_scene, frame_studio_scene
    use studio_selection, only: selection_target_t, selection_patch, selection_control, selection_assembly, selection_mesh
    use studio_single_patch, only: studio_patch_catalog, studio_native_patch, single_patch_count, single_patch_available, &
        single_patch_name, single_patch_mesh_revision, studio_patch_resolution, copy_single_patch_geometry, publish_native_patches
    use studio_units, only: display_length
    use studio_presentation_text, only: integer_text
    implicit none(type, external)
    private

    public :: merge_studio_document

    type :: merge_targets
        logical,allocatable :: patches(:,:),selected_patches(:,:),assemblies(:),selected_assemblies(:),meshes(:)
    end type merge_targets

contains

    logical function merge_studio_document(destination,source,message,selection,added) result(merged)
        type(studio_document_t),intent(inout) :: destination
        type(studio_document_t),intent(in) :: source
        character(:),allocatable,intent(out) :: message
        type(selection_target_t),intent(in),optional :: selection(:)
        type(selection_target_t),allocatable,intent(out),optional :: added(:)
        type(merge_targets) :: targets
        type(studio_patch_catalog) :: patches
        type(multipatch_collection) :: assemblies
        type(studio_mesh_t),allocatable :: meshes(:)
        integer,allocatable :: patch_ids(:,:),assembly_ids(:)
        integer :: prior_assemblies,prior_meshes,i
        real(c_double) :: scale
        logical :: was_framed
        character(:),allocatable :: operation

        merged=.false.
        if(present(added))allocate(added(0))
        call collect_merge_targets(source,selection,targets,message)
        if(len(message)>0)return
        scale=display_length(1._c_double,source%model_unit,destination%model_unit)
        operation='Project import'
        if(present(selection))operation='Paste'
        prior_assemblies=count_multipatches(destination%multipatches)
        prior_meshes=document_mesh_count(destination)
        patches=destination%patches
        call stage_patches(patches,source,targets%patches,scale,patch_ids,message,operation)
        if(len(message)>0)return
        call stage_assemblies(assemblies,source,destination%multipatches,targets%assemblies,patch_ids,patches,assembly_ids,message)
        if(len(message)>0)return
        call stage_meshes(meshes,source,destination,targets%meshes,patch_ids,assembly_ids,assemblies,patches,prior_assemblies,scale,message)
        if(len(message)>0)return

        destination%patches=patches
        if(size(assemblies%groups)>0)then
            if(prior_assemblies==0)then
                call move_alloc(assemblies%groups,destination%multipatches%groups)
            else
                destination%multipatches%groups=[destination%multipatches%groups,assemblies%groups]
            end if
        end if
        do i=1,size(meshes)
            call append_document_mesh(destination,meshes(i))
        end do
        if(size(meshes)>0)destination%mesh_generation=destination%mesh_generation+1
        call record_document_merge(destination,prior_assemblies,prior_meshes,operation)
        was_framed=destination%scene%framed
        call refresh_studio_scene(destination%scene,catalog=destination%patches)
        call restore_transferred_appearance(destination,source,patch_ids)
        if(.not.was_framed)call frame_studio_scene(destination%scene)
        if(present(added))call list_added_objects(targets,patch_ids,assembly_ids,prior_assemblies,prior_meshes,added)
        merged=.true.
    end function merge_studio_document

    subroutine collect_merge_targets(source,selection,targets,message)
        type(studio_document_t),intent(in) :: source
        type(selection_target_t),intent(in),optional :: selection(:)
        type(merge_targets),intent(out) :: targets
        character(:),allocatable,intent(out) :: message
        integer :: counts(3),assemblies,meshes,i,rank,id,group

        counts=[(single_patch_count(rank,catalog=source%patches),rank=1,3)]
        assemblies=count_multipatches(source%multipatches)
        meshes=document_mesh_count(source)
        allocate(targets%patches(maxval(counts),3),source=.false.)
        allocate(targets%selected_patches(maxval(counts),3),source=.false.)
        allocate(targets%assemblies(assemblies),targets%selected_assemblies(assemblies),targets%meshes(meshes),source=.false.)
        if(present(selection))then
            message='A selected object is no longer available to copy.'
            do i=1,size(selection)
                id=selection(i)%id
                rank=selection(i)%rank
                select case(selection(i)%category)
                case(selection_patch,selection_control)
                    if(.not.document_has_patch(source,rank,id))return
                    targets%selected_patches(id,rank)=.true.
                case(selection_assembly)
                    if(id<1.or.id>assemblies)return
                    targets%selected_assemblies(id)=.true.
                case(selection_mesh)
                    if(id<1.or.id>meshes)return
                    targets%meshes(id)=.true.
                case default
                    return
                end select
            end do
        else
            do rank=1,3
                do id=1,counts(rank)
                    targets%selected_patches(id,rank)=document_has_patch(source,rank,id)
                end do
            end do
            targets%selected_assemblies=.true.
            do id=1,meshes
                targets%meshes(id)=document_mesh_valid(source,id)
            end do
        end if
        targets%patches=targets%selected_patches
        targets%assemblies=targets%selected_assemblies
        do id=1,meshes
            if(.not.targets%meshes(id))cycle
            if(.not.document_mesh_valid(source,id))then
                message='Regenerate the selected mesh before copying it.'
                return
            end if
            group=document_mesh_source_group(source,id)
            if(group>0)then
                targets%assemblies(group)=.true.
            else
                rank=source%meshes(id)%kind
                targets%patches(source%meshes(id)%variant,rank)=.true.
            end if
        end do
        message='A referenced assembly member is missing from the source project.'
        do id=1,assemblies
            if(.not.targets%assemblies(id))cycle
            associate(definition=>source%multipatches%groups(id)%definition)
                rank=definition%kind
                if(.not.allocated(definition%members).or.rank<1.or.rank>3)return
                do i=1,size(definition%members)
                    group=definition%members(i)
                    if(.not.single_patch_available(rank,group,catalog=source%patches))return
                    targets%patches(group,rank)=.true.
                end do
            end associate
        end do
        message='There are no objects to transfer.'
        if(any(targets%patches))message=''
    end subroutine collect_merge_targets

    subroutine stage_patches(patches,source,selected,scale,patch_ids,message,operation)
        type(studio_patch_catalog),intent(inout) :: patches
        type(studio_document_t),intent(in) :: source
        logical,intent(in) :: selected(:,:)
        real(c_double),intent(in) :: scale
        integer,allocatable,intent(out) :: patch_ids(:,:)
        character(:),allocatable,intent(out) :: message
        character(*),intent(in) :: operation
        type(studio_native_patch),allocatable :: geometry(:)
        character(:),allocatable :: names(:)
        type(affine_transform_t) :: transform
        integer,allocatable :: resolutions(:,:)
        integer :: rank,id,index,name_length,n,counts(3),kind,variant

        n=count(selected)
        allocate(geometry(n),resolutions(3,n))
        allocate(patch_ids(size(selected,1),size(selected,2)),source=0)
        name_length=1
        do rank=1,3
            do id=1,size(selected,1)
                if(.not.selected(id,rank))cycle
                name_length=max(name_length,len(single_patch_name(rank,id,catalog=source%patches)))
            end do
        end do
        allocate(character(name_length)::names(n))
        counts=[(single_patch_count(rank,catalog=patches),rank=1,3)]
        transform%linear=scale*transform%linear
        index=0
        do rank=1,3
            do id=1,size(selected,1)
                if(.not.selected(id,rank))cycle
                index=index+1
                names(index)=single_patch_name(rank,id,catalog=source%patches)
                call copy_single_patch_geometry(rank,id,geometry(index)%geometry,catalog=source%patches)
                message='The native geometry for '//trim(names(index))//' could not be copied.'
                if(.not.allocated(geometry(index)%geometry))return
                if(scale/=1._c_double)then
                    if(.not.transform_patch(geometry(index)%geometry,transform))return
                end if
                counts(rank)=counts(rank)+1
                patch_ids(id,rank)=counts(rank)
                resolutions(:,index)=studio_patch_resolution(rank,id,catalog=source%patches)
            end do
        end do
        call publish_native_patches(geometry,names,kind,variant,message,resolutions=resolutions,operation=operation,catalog=patches)
    end subroutine stage_patches

    subroutine stage_assemblies(assemblies,source,existing,selected,patch_ids,patches,assembly_ids,message)
        type(multipatch_collection),intent(out) :: assemblies
        type(studio_document_t),intent(in) :: source
        type(multipatch_collection),intent(in) :: existing
        logical,intent(in) :: selected(:)
        integer,intent(in) :: patch_ids(:,:)
        type(studio_patch_catalog),intent(in) :: patches
        integer,allocatable,intent(out) :: assembly_ids(:)
        character(:),allocatable,intent(out) :: message
        type(multipatch_draft),allocatable :: definitions(:)
        integer :: id,n,rank,i,first,second

        allocate(assembly_ids(size(selected)),source=0)
        allocate(definitions(count(selected)))
        message='A source assembly join references missing geometry.'
        n=0
        do id=1,size(selected)
            if(.not.selected(id))cycle
            n=n+1
            assembly_ids(id)=n
            definitions(n)=source%multipatches%groups(id)%definition
            associate(definition=>definitions(n))
                rank=definition%kind
                definition%name=unique_assembly_name(definition%name,existing,definitions(:n-1))
                if(allocated(definition%joins))then
                    do i=1,size(definition%joins)
                        first=definition%joins(i)%patch_a
                        second=definition%joins(i)%patch_b
                        if(.not.any(definition%members==first).or..not.any(definition%members==second))return
                        definition%joins(i)%patch_a=patch_ids(first,rank)
                        definition%joins(i)%patch_b=patch_ids(second,rank)
                    end do
                end if
                definition%members=patch_ids(definition%members,rank)
            end associate
        end do
        call restore_multipatch_definitions(assemblies,definitions,catalog=patches)
        message=''
    end subroutine stage_assemblies

    function unique_assembly_name(base,existing,staged) result(name)
        character(*),intent(in) :: base
        type(multipatch_collection),intent(in) :: existing
        type(multipatch_draft),intent(in) :: staged(:)
        character(:),allocatable :: name
        integer :: serial,i

        name=base
        serial=2
        do
            if(.not.any([(name==existing%groups(i)%definition%name,i=1,count_multipatches(existing))]).and. &
                .not.any([(name==staged(i)%name,i=1,size(staged))]))return
            name=base//' ('//integer_text(serial)//')'
            serial=serial+1
        end do
    end function unique_assembly_name

    subroutine stage_meshes(meshes,source,destination,selected,patch_ids,assembly_ids,assemblies,patches,prior_assemblies,scale,message)
        type(studio_mesh_t),allocatable,intent(out) :: meshes(:)
        type(studio_document_t),intent(in) :: source,destination
        logical,intent(in) :: selected(:)
        integer,intent(in) :: patch_ids(:,:),assembly_ids(:),prior_assemblies
        type(multipatch_collection),intent(in) :: assemblies
        type(studio_patch_catalog),intent(in) :: patches
        real(c_double),intent(in) :: scale
        character(:),allocatable,intent(out) :: message
        integer,allocatable :: members(:)
        integer :: id,n,kind,group,i
        real(rk) :: mesh_scale

        allocate(meshes(count(selected)))
        mesh_scale=real(scale,rk)
        message=''
        n=0
        do id=1,size(selected)
            if(.not.selected(id))cycle
            n=n+1
            associate(original=>source%meshes(id),mesh=>meshes(n))
                kind=original%kind
                mesh=studio_mesh_t(kind=kind,samples=original%samples,nodes=original%nodes,elements=original%elements, &
                    target_size=original%target_size*mesh_scale,visible=original%visible, &
                    per_patch_sizing=original%per_patch_sizing,synchronize=original%synchronize,parts=original%parts)
                group=document_mesh_source_group(source,id)
                if(group>0)then
                    group=assembly_ids(group)
                    associate(assembly=>assemblies%groups(group))
                        if(.not.assembly%validation%valid)then
                            message='The copied assembly could not be validated for its mesh.'
                            return
                        end if
                        mesh%group=prior_assemblies+group
                        mesh%source_group=assembly%definition%name
                        members=assembly%definition%members
                    end associate
                else
                    mesh%variant=patch_ids(original%variant,kind)
                    members=[mesh%variant]
                end if
                mesh%name=unique_mesh_name(original%name,destination,meshes(:n-1))
                do i=1,size(mesh%parts)
                    mesh%parts(i)%variant=members(i)
                    mesh%parts(i)%revision=single_patch_mesh_revision(kind,members(i),catalog=patches)
                    if(mesh_scale/=1._rk)mesh%parts(i)%points=mesh_scale*mesh%parts(i)%points
                end do
                if(allocated(original%connected_points))then
                    mesh%connected_points=original%connected_points
                    mesh%connected_cells=original%connected_cells
                    if(mesh_scale/=1._rk)mesh%connected_points=mesh_scale*mesh%connected_points
                end if
            end associate
        end do
    end subroutine stage_meshes

    function unique_mesh_name(base,destination,staged) result(name)
        character(*),intent(in) :: base
        type(studio_document_t),intent(in) :: destination
        type(studio_mesh_t),intent(in) :: staged(:)
        character(:),allocatable :: name
        integer :: serial,i

        name=base
        serial=2
        do
            if(.not.any([(name==destination%meshes(i)%name,i=1,document_mesh_count(destination))]).and. &
                .not.any([(name==staged(i)%name,i=1,size(staged))]))return
            name=base//' ('//integer_text(serial)//')'
            serial=serial+1
        end do
    end function unique_mesh_name

    subroutine restore_transferred_appearance(destination,source,patch_ids)
        type(studio_document_t),intent(inout) :: destination
        type(studio_document_t),intent(in) :: source
        integer,intent(in) :: patch_ids(:,:)
        integer :: rank,id,variant

        do rank=1,3
            if(.not.allocated(source%scene%catalog(rank)%entries))cycle
            do id=1,min(size(patch_ids,1),size(source%scene%catalog(rank)%entries))
                variant=patch_ids(id,rank)
                if(variant==0)cycle
                associate(old=>source%scene%catalog(rank)%entries(id),new=>destination%scene%catalog(rank)%entries(variant))
                    new%visible=old%visible
                    new%rgb=old%rgb
                    new%opacity=old%opacity
                end associate
            end do
        end do
    end subroutine restore_transferred_appearance

    subroutine list_added_objects(targets,patch_ids,assembly_ids,prior_assemblies,prior_meshes,added)
        type(merge_targets),intent(in) :: targets
        integer,intent(in) :: patch_ids(:,:),assembly_ids(:),prior_assemblies,prior_meshes
        type(selection_target_t),allocatable,intent(out) :: added(:)
        integer :: id,rank,n,mesh

        allocate(added(count(targets%selected_patches)+count(targets%selected_assemblies)+count(targets%meshes)))
        n=0
        do id=1,size(assembly_ids)
            if(.not.targets%selected_assemblies(id))cycle
            n=n+1
            added(n)=selection_target_t(selection_assembly,prior_assemblies+assembly_ids(id))
        end do
        mesh=prior_meshes
        do id=1,size(targets%meshes)
            if(.not.targets%meshes(id))cycle
            n=n+1
            mesh=mesh+1
            added(n)=selection_target_t(selection_mesh,mesh)
        end do
        do rank=1,3
            do id=1,size(patch_ids,1)
                if(.not.targets%selected_patches(id,rank))cycle
                n=n+1
                added(n)=selection_target_t(selection_patch,patch_ids(id,rank),rank)
            end do
        end do
    end subroutine list_added_objects
end module studio_project_merge
