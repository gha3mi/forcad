module studio_workspace_selection
    use studio_geometry_operations, only: read_geometry_target
    use studio_multipatch_collection, only: count_multipatches, multipatch_undo_stamp, multipatch_contains
    use qt, only: Qt_ControlModifier, Qt_ShiftModifier
    use studio_selection, only: selection_target_t, selection_patch, selection_assembly, selection_control, selection_mesh, &
        selection_replace, selection_add, selection_toggle, pick_vertex
    use studio_single_patch, only: active_patch_id, single_patch_count, single_patch_revision
    use studio_state, only: studio_t, studio_mesh_count, studio_mesh_source_group, studio_mesh_valid
    implicit none(type, external)
    private

    public :: selected_group, selected_mesh, workspace_selection_key, primary_workspace_target, pointer_selection_mode
    public :: patch_in_selection, assembly_member_in_selection, mesh_in_selection, mesh_cell_in_selection, control_in_selection
    public :: selection_measurement_targets, selection_mesh_records

contains

    function selection_mesh_records(s,targets) result(ids)
        type(studio_t), intent(in) :: s
        type(selection_target_t), intent(in), optional :: targets(:)
        integer, allocatable :: ids(:)
        logical, allocatable :: explicit_meshes(:), assemblies(:), patches(:)
        integer :: offsets(0:3),mesh_count,rank,i,matching
        logical :: has_assemblies

        mesh_count=studio_mesh_count(s)
        if(mesh_count==0)then
            allocate(ids(0))
            return
        end if
        offsets(0)=0
        do rank=1,3
            offsets(rank)=offsets(rank-1)+single_patch_count(rank,catalog=s%document%patches)
        end do
        allocate(explicit_meshes(mesh_count),assemblies(count_multipatches(s%document%multipatches)),patches(offsets(3)),source=.false.)
        if(present(targets))then
            do i=1,size(targets)
                call mark_target(targets(i))
            end do
        else if(allocated(s%selection%items))then
            do i=1,s%selection%count()
                call mark_target(s%selection%items(i))
            end do
        else
            call mark_target(primary_workspace_target(s))
        end if
        has_assemblies=any(assemblies)
        allocate(ids(mesh_count))
        matching=0
        do i=1,mesh_count
            if(.not.explicit_meshes(i).and..not.source_selected(i))cycle
            if(.not.studio_mesh_valid(s,i))cycle
            matching=matching+1
            ids(matching)=i
        end do
        ids=ids(:matching)

    contains
        subroutine mark_target(target)
            type(selection_target_t), intent(in) :: target
            integer :: id,dimension

            id=target%id
            if(id<1)return
            select case(target%category)
            case(selection_mesh)
                if(id<=size(explicit_meshes))explicit_meshes(id)=.true.
            case(selection_assembly)
                if(id<=size(assemblies))assemblies(id)=.true.
            case(selection_patch,selection_control)
                dimension=target%rank
                if(dimension<1.or.dimension>3)return
                if(id<=offsets(dimension)-offsets(dimension-1))patches(offsets(dimension-1)+id)=.true.
            case default
                continue
            end select
        end subroutine mark_target

        logical function source_selected(id) result(selected)
            integer, intent(in) :: id
            integer :: group,dimension,part,variant

            selected=.false.
            associate(mesh=>s%document%meshes(id))
                if(mesh%group>0.and.has_assemblies)then
                    group=studio_mesh_source_group(s,id)
                    if(group>0)selected=assemblies(group)
                    if(selected)return
                end if
                dimension=mesh%kind
                if(dimension<1.or.dimension>3)return
                if(.not.allocated(mesh%parts))return
                if(size(mesh%parts)==0)return
                do part=1,size(mesh%parts)
                    variant=mesh%parts(part)%variant
                    if(variant<1.or.variant>offsets(dimension)-offsets(dimension-1))return
                    if(.not.patches(offsets(dimension-1)+variant))return
                end do
                selected=.true.
            end associate
        end function source_selected
    end function selection_mesh_records

    function selection_measurement_targets(s) result(targets)
        type(studio_t), intent(in) :: s
        type(selection_target_t), allocatable :: targets(:)
        integer :: i

        allocate(targets(0))
        if(.not.allocated(s%selection%items))then
            call append(primary_workspace_target(s))
            return
        end if
        do i=1,s%selection%count()
            call append(s%selection%items(i))
        end do
    contains
        subroutine append(target)
            type(selection_target_t), intent(in) :: target
            integer, allocatable :: ids(:)
            character(:), allocatable :: name
            integer :: rank,j

            select case(target%category)
            case(selection_patch)
                if(target%feature/=pick_vertex)targets=[targets,target]
            case(selection_assembly)
                call read_geometry_target(s%document,0,0,target%id,rank,ids,name)
                targets=[targets,(selection_target_t(selection_patch,ids(j),rank),j=1,size(ids))]
            case default
                continue
            end select
        end subroutine append
    end function selection_measurement_targets

    integer function pointer_selection_mode(modifiers) result(mode)
        integer, intent(in) :: modifiers

        mode=selection_replace
        if(iand(modifiers,Qt_ShiftModifier)/=0)mode=selection_add
        if(iand(modifiers,Qt_ControlModifier)/=0)mode=selection_toggle
    end function pointer_selection_mode

    function primary_workspace_target(s) result(target)
        type(studio_t), intent(in) :: s
        type(selection_target_t) :: target

        if(allocated(s%selection%items))then
            if(s%selection%count()>0)target=s%selection%items(s%selection%count())
            return
        end if
        if(selected_mesh(s)>0)then
            target=selection_target_t(selection_mesh,s%selected_mesh,part=s%mesh_cell(1),item=s%mesh_cell(2),feature=s%mesh_cell(3))
        else if(selected_group(s)>0)then
            target=selection_target_t(selection_assembly,s%selected_multipatch)
        else if(s%live_scene%valid)then
            if(s%selected_control>0)then
                target=selection_target_t(selection_control,s%live_scene%variant,s%entity_kind,item=s%selected_control)
            else
                target=selection_target_t(selection_patch,s%live_scene%variant,s%entity_kind)
            end if
        end if
    end function primary_workspace_target

    logical function patch_in_selection(s,rank,id,include_controls) result(selected)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: rank,id
        logical, intent(in), optional :: include_controls
        integer :: i

        if(.not.allocated(s%selection%items))then
            selected=s%selected_multipatch==0.and.s%selected_mesh==0.and.s%live_scene%valid.and. &
                s%live_scene%dimension==rank.and.s%live_scene%variant==id
            return
        end if
        selected=.false.
        do i=1,s%selection%count()
            associate(target=>s%selection%items(i))
                if(target%rank/=rank.or.target%id/=id)cycle
                selected=target%category==selection_patch
                if(present(include_controls))selected=selected.or.(include_controls.and.target%category==selection_control)
                if(selected)return
            end associate
        end do
    end function patch_in_selection

    logical function assembly_member_in_selection(s,rank,id) result(selected)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: rank,id
        integer :: i

        if(.not.allocated(s%selection%items))then
            selected=multipatch_contains(s%document%multipatches,s%selected_multipatch,rank,id)
            return
        end if
        selected=.false.
        do i=1,s%selection%count()
            associate(target=>s%selection%items(i))
                if(target%category/=selection_assembly)cycle
                selected=multipatch_contains(s%document%multipatches,target%id,rank,id)
                if(selected)return
            end associate
        end do
    end function assembly_member_in_selection

    logical function mesh_in_selection(s,id) result(selected)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: id
        integer :: i

        if(.not.allocated(s%selection%items))then
            selected=id==s%selected_mesh
            return
        end if
        selected=.false.
        do i=1,s%selection%count()
            associate(target=>s%selection%items(i))
                selected=target%category==selection_mesh.and.target%id==id
                if(selected)return
            end associate
        end do
    end function mesh_in_selection

    logical function mesh_cell_in_selection(s,id,part,cell) result(selected)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: id,part,cell

        if(.not.allocated(s%selection%items))then
            selected=id==s%selected_mesh.and.part==s%mesh_cell(1).and.cell==s%mesh_cell(2)
        else
            selected=s%selection%contains(selection_target_t(selection_mesh,id,part=part,item=cell))
        end if
    end function mesh_cell_in_selection

    logical function control_in_selection(s,point) result(selected)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: point

        if(.not.allocated(s%selection%items))then
            selected=point==s%selected_control
        else
            selected=s%selection%contains(selection_target_t(selection_control,s%live_scene%variant,s%entity_kind,item=point))
        end if
    end function control_in_selection

    integer function selected_group(s) result(index)
        type(studio_t), intent(in) :: s

        index=0
        if(.not.allocated(s%document%multipatches%groups))return
        if(s%selected_multipatch<1.or.s%selected_multipatch>size(s%document%multipatches%groups))return
        index=s%selected_multipatch
    end function selected_group

    integer function selected_mesh(s) result(index)
        type(studio_t), intent(in) :: s

        index=0
        if(s%selected_mesh>=1.and.s%selected_mesh<=studio_mesh_count(s))index=s%selected_mesh
    end function selected_mesh

    function workspace_selection_key(s) result(key)
        type(studio_t), intent(in) :: s
        character(:), allocatable :: key
        character(120) :: numbers
        integer :: group,mesh,variant,i

        mesh=selected_mesh(s)
        group=selected_group(s)
        if(mesh>0)then
            associate(item=>s%document%meshes(mesh))
                write(numbers,'(i0,6(1x,i0),1x,l1)') mesh,item%kind,item%nodes,item%elements,item%samples, &
                    studio_mesh_valid(s,mesh)
                key='mesh '//trim(numbers)//' '//item%name
                if(allocated(item%parts))then
                    do i=1,size(item%parts)
                        write(numbers,'(i0,1x,i0)') item%parts(i)%variant,item%parts(i)%revision
                        key=key//' '//trim(numbers)
                    end do
                end if
            end associate
        else if(group>0)then
            associate(definition=>s%document%multipatches%groups(group)%definition)
                write(numbers,'(i0,2(1x,i0))') group,definition%kind,multipatch_undo_stamp(s%document%multipatches,catalog=s%document%patches)
                key='group '//trim(numbers)//' '//definition%name
                if(allocated(definition%members))then
                    do i=1,size(definition%members)
                        variant=definition%members(i)
                        write(numbers,'(i0,1x,i0)') variant,single_patch_revision(definition%kind,variant,catalog=s%document%patches)
                        key=key//' '//trim(numbers)
                    end do
                end if
            end associate
        else if(s%live_scene%valid)then
            variant=active_patch_id(s%entity_kind,catalog=s%document%patches)
            write(numbers,'(i0,2(1x,i0))') s%entity_kind,variant,single_patch_revision(s%entity_kind,variant,catalog=s%document%patches)
            key='patch '//trim(numbers)
        else
            key='workspace'
        end if
    end function workspace_selection_key

end module studio_workspace_selection
