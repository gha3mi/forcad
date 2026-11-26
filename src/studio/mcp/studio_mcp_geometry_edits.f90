module studio_mcp_geometry_edits
    use, intrinsic :: iso_fortran_env, only: int64
    use qt, only: QJsonArray, QJsonObject, QJsonValue_Array, QJsonValue_Bool, QJsonValue_Object, QJsonValue_String
    use studio_state, only: studio_t, dp
    use studio_scene, only: project_patch_exists
    use studio_single_patch, only: single_patch_scene, active_patch_id, set_active_patch, add_single_patch, &
        read_single_patch_scene, commit_single_patch_knot_values, &
        commit_single_patch_knot, commit_single_patch_degree, commit_single_patch_refinement, commit_single_patch_interval, &
        commit_single_patch_reverse
    use studio_multipatch, only: multipatch_draft
    use studio_multipatch_refinement, only: commit_multipatch_refinement
    use studio_mcp_refinement, only: read_mcp_refinement
    use studio_multipatch_collection, only: count_multipatches, save_multipatch
    use studio_mcp_controls, only: edit_mcp_controls
    use studio_mesh_operations, only: mesh_source_t, read_document_mesh_source, estimate_document_mesh_samples
    use studio_edit_history, only: document_undo_stamp
    use studio_document, only: document_patch_in_assembly
    use studio_geometry_operations, only: commit_patch_range, commit_patch_section
    use studio_ui, only: rename_project_patch, on_assistant_model_created, on_assistant_geometry_changed, &
        generate_workspace_mesh, apply_patch_range, rebuild_project_tree, delete_project_patch, delete_project_assembly, &
        delete_project_mesh, delete_project_selection
    use studio_selection, only: selection_target_t,selection_patch,selection_assembly,selection_mesh
    use studio_mcp_arguments, only: only_fields, bounded_integer, number_arg, vector_arg, real_array_arg, integer_array_arg, &
        optional_boolean, axis_number, read_geometry_identity
    use studio_mcp_inspection, only: describe_selection, patch_identity, assembly_identity, mesh_identity, active_rank, patch_kind
    use studio_presentation_text, only: integer_text
    implicit none(type, external)
    private

    public :: edit_patch, copy_patch, rename_patch, delete_existing, assemble_patches, mesh_selection, preview_mesh

contains

    subroutine edit_patch(s,args,result,error)
        type(studio_t),intent(inout),target :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        type(single_patch_scene) :: patch
        real(dp),allocatable :: knots(:)
        type(QJsonObject) :: item
        character(:),allocatable :: op,direction,kind_name
        integer :: kind,id,axis,repeats,applied,divisions,steps,created,face,previous
        real(dp) :: tolerance,error_bound,first,last
        logical :: changed,reducing

        op=args%stringValue('op')
        if(op=='control')then
            call edit_mcp_controls(s,args,result,error,kind,id,changed)
            if(len(error)>0)return
            if(changed)call on_assistant_geometry_changed(s,kind,id)
            call describe_selection(s,result)
            return
        end if
        if(args%contains('target'))then
            error='Supply an available curve, surface, volume, or assembly target with its kind and ID.'
            if(.not.read_geometry_identity(args,kind_name,id))return
            if(kind_name=='assembly'.or.kind_name=='multipatch')then
                call edit_assembly(s,args,result,error,id,op)
                return
            end if
            kind=patch_kind(kind_name)
            if(kind==0)return
            if(.not.project_patch_exists(s%document%show_templates,kind,id,patches=s%document%patches))return
        else
            if(s%selected_multipatch>0)then
                error='Select one multipatch or supply its assembly ID.'
                if(s%selection%count()>1)return
                call edit_assembly(s,args,result,error,s%selected_multipatch,op)
                return
            end if
            kind=active_rank(s)
            error='Select a curve, surface or volume patch before editing.'
            if(kind==0)return
            id=active_patch_id(kind,catalog=s%document%patches)
        end if
        if(op=='split'.or.op=='reverse')then
            error='Detach this patch from its multipatch before splitting or reversing it.'
            if(document_patch_in_assembly(s%document,kind,id))return
        end if
        direction=args%stringValue('axis')
        axis=axis_number(direction,'UVW')
        error='Supply the fields required by this patch edit.'
        changed=.false.
        select case(op)
        case('move_knot')
            if(.not.only_fields(args,'|op|target|axis|value|to|'))return
            if(axis<1.or.axis>kind)return
            if(.not.number_arg(args,'value',first).or..not.number_arg(args,'to',last))return
            call read_single_patch_scene(kind,patch,variant=id,catalog=s%document%patches)
            error='ForCAD could not read the current knot vector.'
            if(.not.patch%valid)return
            knots=patch%axes(axis)%knots
            error='That value is absent from the current knot vector.'
            if(.not.any(knots==first))return
            where(knots==first)knots=last
            call commit_single_patch_knot_values(kind,axis,knots,changed,error,catalog=s%document%patches,patch_id=id)
        case('set_knot_values')
            if(.not.only_fields(args,'|op|target|axis|values|'))return
            if(axis<1.or.axis>kind)return
            if(.not.real_array_arg(args,'values',knots))return
            call commit_single_patch_knot_values(kind,axis,knots,changed,error,catalog=s%document%patches,patch_id=id)
        case('insert_knot','remove_knot')
            if(.not.only_fields(args,'|op|target|axis|value|count|'))return
            if(axis<1.or.axis>kind)return
            repeats=bounded_integer(args,'count',1,huge(1),fallback=1)
            if(repeats<1.or..not.number_arg(args,'value',first))return
            call commit_single_patch_knot(kind,axis,first,repeats,op=='remove_knot',applied,error,catalog=s%document%patches,patch_id=id)
            if(len(error)==0)call result%insert('applied',applied)
            changed=applied>0
        case('elevate_degree','reduce_degree')
            if(.not.only_fields(args,'|op|target|axis|count|tolerance|'))return
            if(axis<1.or.axis>kind)return
            repeats=bounded_integer(args,'count',1,huge(1),fallback=1)
            if(repeats<1)return
            tolerance=0._dp
            reducing=op=='reduce_degree'
            if(reducing)then
                if(args%contains('tolerance'))then
                    if(.not.number_arg(args,'tolerance',tolerance))return
                    if(tolerance<0._dp)return
                end if
            else if(args%contains('tolerance'))then
                return
            end if
            call commit_single_patch_degree(kind,axis,repeats,reducing,tolerance,applied,error_bound,error,catalog=s%document%patches,patch_id=id)
            if(len(error)==0)then
                call result%insert('applied',applied)
                call result%insert('error_bound',error_bound)
            end if
            changed=applied>0
        case('refine')
            if(.not.only_fields(args,'|op|target|axis|divisions|degree_steps|'))return
            if(direction/='all'.and.args%contains('axis'))then
                if(axis<1.or.axis>kind)return
            end if
            divisions=bounded_integer(args,'divisions',1,huge(1),fallback=1)
            steps=bounded_integer(args,'degree_steps',0,huge(1),fallback=0)
            if(divisions<1.or.steps<0)return
            call commit_single_patch_refinement(kind,axis,divisions,steps,changed,error,catalog=s%document%patches,patch_id=id)
        case('set_interval')
            if(.not.only_fields(args,'|op|target|axis|lower|upper|'))return
            if(axis<1.or.axis>kind)return
            if(.not.number_arg(args,'lower',first).or..not.number_arg(args,'upper',last))return
            call commit_single_patch_interval(kind,axis,first,last,changed,error,catalog=s%document%patches,patch_id=id)
        case('split','extract')
            if(axis<1.or.axis>kind)return
            if(op=='split')then
                if(.not.only_fields(args,'|op|target|axis|value|'))return
                if(.not.number_arg(args,'value',first))return
                last=first
            else
                if(.not.only_fields(args,'|op|target|axis|lower|upper|'))return
                if(.not.number_arg(args,'lower',first).or..not.number_arg(args,'upper',last))return
            end if
            if(active_rank(s)==kind.and.s%live_scene%variant==id)then
                call apply_patch_range(s,axis,first,last,op=='split',created,error)
            else
                previous=active_patch_id(kind,catalog=s%document%patches)
                call commit_patch_range(s%document,kind,id,axis,first,last,op=='split',created,error)
                if(len(error)>0)return
                if(previous>0)call set_active_patch(kind,previous,catalog=s%document%patches)
                call rebuild_project_tree(s)
                call on_assistant_geometry_changed(s,kind,id)
            end if
            if(len(error)>0)return
            item=patch_identity(s,kind,created)
            call result%insert('created',item)
            call item%delete()
            call result%insert('changed',.true.)
            call describe_selection(s,result)
            return
        case('extract_section','extract_face')
            error='Choose a surface or volume to extract a section.'
            if(kind<2)return
            if(op=='extract_face')then
                error='Choose a volume face numbered 1–6 in U min, U max, V min, V max, W min, W max order.'
                if(kind/=3.or..not.only_fields(args,'|op|target|face|'))return
                face=bounded_integer(args,'face',1,6)
                if(face<1)return
                call read_single_patch_scene(kind,patch,variant=id,catalog=s%document%patches)
                if(.not.patch%valid)return
                axis=(face+1)/2
                first=patch%bounds(2-mod(face,2),axis)
            else
                error='Supply an active U, V, or W direction and a literal parameter value.'
                if(.not.only_fields(args,'|op|target|axis|value|').or.axis<1.or.axis>kind)return
                if(.not.number_arg(args,'value',first))return
            end if
            call commit_patch_section(s%document,kind,id,axis,first,created,error)
            if(len(error)>0)return
            call on_assistant_model_created(s,0,kind-1,created,0)
            item=patch_identity(s,kind-1,created)
            call result%insert('created',item)
            call item%delete()
            call result%insert('changed',.true.)
            call describe_selection(s,result)
            return
        case('reverse')
            if(.not.only_fields(args,'|op|target|axis|'))return
            if(axis<1.or.axis>kind)return
            call commit_single_patch_reverse(kind,axis,changed,error,catalog=s%document%patches,patch_id=id)
        case default
            error='Unknown patch edit.'
        end select
        if(len(error)>0)return
        if(changed)call on_assistant_geometry_changed(s,kind,id)
        item=patch_identity(s,kind,id)
        call result%insert('target',item)
        call item%delete()
        call result%insert('changed',changed)
        call describe_selection(s,result)
    end subroutine edit_patch

    subroutine edit_assembly(s,args,result,error,id,op)
        type(studio_t),intent(inout),target :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        integer,intent(in) :: id
        character(*),intent(in) :: op
        type(QJsonObject) :: item
        integer :: kind,member,axis,steps
        real(dp) :: knot
        logical :: changed

        error='Select an available multipatch.'
        if(id<1.or.id>count_multipatches(s%document%multipatches))return
        kind=s%document%multipatches%groups(id)%definition%kind
        call read_mcp_refinement(args,op,member,axis,steps,knot,error)
        if(len(error)>0)return
        if(op=='insert_knot')then
            call commit_multipatch_refinement(s%document,id,member,axis,steps,error,changed,knot=knot)
        else
            call commit_multipatch_refinement(s%document,id,member,axis,steps,error,changed)
        end if
        if(len(error)>0)return

        if(changed)call on_assistant_geometry_changed(s,kind,member)
        item=assembly_identity(s,id)
        call result%insert('target',item)
        call item%delete()
        call result%insert('changed',changed)
        call describe_selection(s,result)
    end subroutine edit_assembly

    subroutine copy_patch(s,args,result,error)
        type(studio_t),intent(inout),target :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        character(:),allocatable :: kind_name,name
        real(dp) :: offset(3)
        integer :: kind,id,created
        type(QJsonObject) :: item

        error='Copy a patch by kind and ID, with a name and an optional XYZ offset in model units.'
        if(.not.only_fields(args,'|kind|id|name|offset|').or.args%valueType('name')/=QJsonValue_String)return
        kind_name=args%stringValue('kind')
        kind=patch_kind(kind_name)
        id=bounded_integer(args,'id',1,huge(1))
        name=trim(adjustl(args%stringValue('name')))
        if(kind<1.or.id<1.or.len(name)==0)return
        if(.not.project_patch_exists(s%document%show_templates,kind,id,patches=s%document%patches))return
        offset=0._dp
        if(args%contains('offset'))then
            if(.not.vector_arg(args,'offset',offset))return
        end if
        call add_single_patch(kind,id,name,created,error,offset,catalog=s%document%patches)
        if(len(error)>0)return
        if(created<1)then
            error='ForCAD Studio could not copy the patch.'
            return
        end if
        call on_assistant_model_created(s,0,kind,created,0)
        item=patch_identity(s,kind,created)
        call result%insert('created',item)
        call item%delete()
        call result%insert('created_kind',kind_name)
        call result%insert('created_id',created)
        call describe_selection(s,result)
    end subroutine copy_patch

    subroutine rename_patch(s,args,result,error)
        type(studio_t),intent(inout) :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        type(QJsonObject) :: item
        character(:),allocatable :: name
        integer :: kind,id
        logical :: changed

        error='Rename a patch with its kind and ID and a name of 1 to 96 characters.'
        if(.not.only_fields(args,'|kind|id|name|'))return
        kind=patch_kind(args%stringValue('kind'))
        id=bounded_integer(args,'id',1,huge(1))
        name=trim(adjustl(args%stringValue('name')))
        if(kind<1.or.id<1.or.len(name)==0.or.len(name)>96)return
        call rename_project_patch(s,kind,id,name,error,changed)
        if(len(error)>0)return
        item=patch_identity(s,kind,id)
        call result%insert('target',item)
        call item%delete()
        call result%insert('changed',changed)
        call result%insert('undo_available',document_undo_stamp(s%document)>0_int64)
    end subroutine rename_patch

    subroutine delete_existing(s,args,result,error)
        type(studio_t),intent(inout),target :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        type(selection_target_t), allocatable :: targets(:)
        type(QJsonArray) :: entries
        type(QJsonObject) :: entry
        integer :: i,removed
        logical :: valid

        if(args%contains('targets'))then
            error='Supply an array of exact object kinds and IDs from inspect.'
            if(.not.only_fields(args,'|targets|').or.args%valueType('targets')/=QJsonValue_Array)return
            entries=args%arrayValue('targets')
            allocate(targets(entries%size()))
            do i=1,size(targets)
                valid=entries%typeAt(i-1)==QJsonValue_Object
                if(valid)then
                    entry=entries%objectAt(i-1)
                    valid=read_delete_target(entry,targets(i))
                    call entry%delete()
                end if
                if(valid)cycle
                call entries%delete()
                return
            end do
            call entries%delete()
            call delete_project_selection(s,targets,removed,error)
        else if(args%contains('selection'))then
            error='Set selection to true to delete the selected whole objects.'
            if(.not.only_fields(args,'|selection|').or..not.args%boolValue('selection'))return
            if(s%selection%count()==0)return
            targets=s%selection%items
            call delete_project_selection(s,targets,removed,error)
        else
            call delete_one_object(s,args,error)
            removed=1
        end if
        if(len(error)>0)return
        call result%insert('deleted',.true.)
        call result%insert('deleted_count',removed)
        call describe_selection(s,result)
    end subroutine delete_existing

    logical function read_delete_target(args,target) result(valid)
        type(QJsonObject), intent(in) :: args
        type(selection_target_t), intent(out) :: target
        character(:), allocatable :: name
        integer :: id,rank

        valid=.false.
        if(.not.only_fields(args,'|kind|id|'))return
        id=bounded_integer(args,'id',1,huge(1))
        if(id<1)return
        name=args%stringValue('kind')
        select case(name)
        case('assembly')
            target=selection_target_t(selection_assembly,id)
        case('mesh')
            target=selection_target_t(selection_mesh,id)
        case default
            rank=patch_kind(name)
            if(rank==0)return
            target=selection_target_t(selection_patch,id,rank)
        end select
        valid=.true.
    end function read_delete_target

    subroutine delete_one_object(s,args,error)
        type(studio_t),intent(inout),target :: s
        type(QJsonObject),intent(in) :: args
        character(:), allocatable, intent(out) :: error
        type(selection_target_t) :: target

        error='Delete requires an object kind and a positive ID from inspect.'
        if(.not.read_delete_target(args,target))return
        select case(target%category)
        case(selection_assembly)
            call delete_project_assembly(s,target%id,error)
        case(selection_mesh)
            call delete_project_mesh(s,target%id,error)
        case(selection_patch)
            call delete_project_patch(s,target%rank,target%id,error)
        case default
            return
        end select
    end subroutine delete_one_object

    subroutine assemble_patches(s,args,result,error)
        use studio_single_patch, only: single_patch_count
        type(studio_t),intent(inout),target :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        type(multipatch_draft) :: draft
        type(QJsonObject) :: join,summary
        type(QJsonArray) :: list,flags
        integer :: i,id,available
        real(dp) :: order
        logical :: valid,valid_order

        error='Supply a name, one patch kind, member IDs, and explicit joins.'
        if(.not.only_fields(args,'|id|name|kind|members|joins|'))return
        draft%name=trim(args%stringValue('name'))
        if(len(draft%name)==0.or.index(draft%name,achar(0))/=0)return
        draft%kind=patch_kind(args%stringValue('kind'))
        if(draft%kind==0)return
        available=single_patch_count(draft%kind,catalog=s%document%patches)
        if(.not.integer_array_arg(args,'members',draft%members,2,available,available))return
        do i=1,size(draft%members)
            if(.not.project_patch_exists(s%document%show_templates,draft%kind,draft%members(i),patches=s%document%patches))return
        end do
        if(args%valueType('joins')/=QJsonValue_Array)return
        list=args%arrayValue('joins')
        allocate(draft%joins(list%size()))
        do i=1,size(draft%joins)
            valid=list%typeAt(i-1)==QJsonValue_Object
            if(valid)then
                join=list%objectAt(i-1)
                valid=only_fields(join,'|patch_a|side_a|patch_b|side_b|continuity|reverse|swap|flip|geometric|reparameterization|')
                associate(spec=>draft%joins(i))
                    spec%patch_a=bounded_integer(join,'patch_a',1,huge(1))
                    spec%patch_b=bounded_integer(join,'patch_b',1,huge(1))
                    spec%side_a=bounded_integer(join,'side_a',1,2*draft%kind)
                    spec%side_b=bounded_integer(join,'side_b',1,2*draft%kind)
                    valid=valid.and.min(spec%patch_a,spec%patch_b,spec%side_a,spec%side_b)>0
                    if(join%contains('continuity'))then
                        valid_order=number_arg(join,'continuity',order)
                        if(valid_order)valid_order=order>=-1._dp.and.order<=real(huge(1),dp).and.order==aint(order)
                        if(valid_order)spec%continuity=int(order)
                        valid=valid.and.valid_order
                    end if
                    valid=valid.and.optional_boolean(join,'reverse').and.optional_boolean(join,'swap').and. &
                        optional_boolean(join,'geometric')
                    spec%reverse=join%boolValue('reverse')
                    spec%swap=join%boolValue('swap')
                    spec%geometric=join%boolValue('geometric')
                    if(spec%geometric.and.spec%continuity>0)then
                        valid_order=real_array_arg(join,'reparameterization',spec%reparameterization,spec%continuity)
                        valid=valid.and.valid_order
                    else
                        valid=valid.and..not.join%contains('reparameterization')
                    end if
                    if(join%contains('flip'))then
                        if(join%valueType('flip')==QJsonValue_Array)then
                            flags=join%arrayValue('flip')
                            valid=valid.and.flags%size()==2.and.flags%typeAt(0)==QJsonValue_Bool.and. &
                                flags%typeAt(1)==QJsonValue_Bool
                            spec%flip=[flags%boolAt(0),flags%boolAt(1)]
                            call flags%delete()
                        else
                            valid=.false.
                        end if
                    end if
                end associate
                call join%delete()
            end if
            if(valid)cycle
            call list%delete()
            error='Invalid assembly join '//integer_text(i)//'.'
            return
        end do
        call list%delete()
        id=bounded_integer(args,'id',1,count_multipatches(s%document%multipatches),fallback=0)
        if(id<0)then
            error='The assembly ID is unavailable.'
            return
        end if
        call save_multipatch(s%document%multipatches,id,draft,error,catalog=s%document%patches)
        if(len(error)>0)return
        call on_assistant_model_created(s,id,draft%kind,0,0)
        summary=assembly_identity(s,id)
        call result%insert('assembly',summary)
        call summary%delete()
    end subroutine assemble_patches

    subroutine mesh_selection(s,args,result,error)
        type(studio_t),intent(inout),target :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        type(QJsonObject) :: item
        type(mesh_source_t) :: source
        integer,allocatable :: cells(:),patch_cells(:,:)
        integer :: id
        real(dp) :: target_size
        logical :: conforming,synchronize

        error='Supply cells, patch_cells, or a positive target_size in model units.'
        if(.not.only_fields(args,'|target|cells|patch_cells|target_size|synchronize|conforming|'))return
        if(.not.optional_boolean(args,'conforming').or..not.optional_boolean(args,'synchronize'))return
        if(count([args%contains('cells'),args%contains('patch_cells'),args%contains('target_size')])/=1)return
        if(args%contains('synchronize').and..not.args%contains('target_size'))return
        if(.not.resolve_mesh_source(s,args,source,error))return
        conforming=.false.
        if(source%mesh_id>0)conforming=allocated(s%document%meshes(source%mesh_id)%connected_cells)
        if(args%contains('conforming'))conforming=args%boolValue('conforming')
        if(args%contains('target_size'))then
            error='Enter a positive finite target size in model units.'
            if(.not.number_arg(args,'target_size',target_size))return
            if(target_size<=0._dp)return
            synchronize=.true.
            if(args%contains('synchronize'))synchronize=args%boolValue('synchronize')
            call generate_workspace_mesh(s,conforming=conforming,mesh_id=id,failure=error,source_override=source, &
                target_size=target_size,synchronize=synchronize)
        else if(args%contains('cells'))then
            error='Supply one positive cell count for each active parameter direction.'
            if(.not.integer_array_arg(args,'cells',cells,1,3,huge(1)-1))return
            call generate_workspace_mesh(s,requested_cells=cells,conforming=conforming, &
                mesh_id=id,failure=error,source_override=source)
        else
            if(.not.read_member_cell_counts(source,args,patch_cells,error))return
            call generate_workspace_mesh(s,patch_cells=patch_cells,conforming=conforming, &
                mesh_id=id,failure=error,source_override=source)
        end if
        if(len(error)>0)return
        item=mesh_identity(s,id)
        call result%insert('mesh',item)
        call item%delete()
    end subroutine mesh_selection

    subroutine preview_mesh(s,args,result,error)
        type(studio_t),intent(in) :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        type(mesh_source_t) :: source
        type(QJsonArray) :: parts,counts
        type(QJsonObject) :: part
        integer,allocatable :: cells(:,:)
        integer :: member
        real(dp) :: target_size
        logical :: synchronize

        error='Supply a positive target_size in model units and optionally synchronize.'
        if(.not.only_fields(args,'|target|target_size|synchronize|'))return
        if(.not.optional_boolean(args,'synchronize'))return
        if(.not.number_arg(args,'target_size',target_size))return
        if(target_size<=0._dp)return
        if(.not.resolve_mesh_source(s,args,source,error))return
        synchronize=.true.
        if(args%contains('synchronize'))synchronize=args%boolValue('synchronize')
        call estimate_document_mesh_samples(s%document,source,target_size,synchronize,cells,error)
        if(len(error)>0)return

        parts=QJsonArray()
        do member=1,size(source%members)
            part=QJsonObject()
            counts=QJsonArray()
            call counts%appendInts(cells(:,member))
            call part%insert('patch_id',source%members(member))
            call part%insert('cells',counts)
            call parts%append(part)
            if(source%assembly_id==0)call result%insert('cells',counts)
            call counts%delete()
            call part%delete()
        end do
        call result%insert('patch_cells',parts)
        call parts%delete()
        call result%insert('target_size',target_size)
        call result%insert('synchronize',synchronize)
    end subroutine preview_mesh

    logical function resolve_mesh_source(s,args,source,error) result(valid)
        type(studio_t),intent(in) :: s
        type(QJsonObject),intent(in) :: args
        type(mesh_source_t),intent(out) :: source
        character(:),allocatable,intent(out) :: error
        character(:),allocatable :: kind_name
        integer :: id,kind,patch

        valid=.false.
        if(args%contains('target'))then
            error='Supply an available patch, assembly, or mesh target with its kind and ID.'
            if(.not.read_geometry_identity(args,kind_name,id))return
            select case(kind_name)
            case('assembly')
                source=read_document_mesh_source(s%document,0,0,id,0)
            case('mesh')
                source=read_document_mesh_source(s%document,0,0,0,id)
            case default
                kind=patch_kind(kind_name)
                if(kind==0)return
                source=read_document_mesh_source(s%document,kind,id,0,0)
            end select
        else
            error='Select a patch or assembly with an available mesh source.'
            patch=0
            if(s%live_scene%valid)patch=s%live_scene%variant
            source=read_document_mesh_source(s%document,s%entity_kind,patch,s%selected_multipatch,s%selected_mesh)
        end if
        if(.not.allocated(source%members))return
        valid=size(source%members)>0
        if(valid)error=''
    end function resolve_mesh_source

    logical function read_member_cell_counts(source,args,counts,error) result(valid)
        type(mesh_source_t),intent(in) :: source
        type(QJsonObject),intent(in) :: args
        integer,allocatable,intent(out) :: counts(:,:)
        character(:),allocatable,intent(out) :: error
        type(QJsonArray) :: rows
        type(QJsonObject) :: row
        integer,allocatable :: values(:)
        integer :: i,id,column

        valid=.false.
        error='Select a multipatch or its mesh to set per-patch cell counts.'
        if(source%assembly_id<1.or..not.allocated(source%members))return
        if(size(source%members)==0.or.source%kind<1.or.source%kind>3)return
        error='Supply one patch_cells row per member, with its patch_id and one positive count per active direction.'
        if(args%valueType('patch_cells')/=QJsonValue_Array)return
        rows=args%arrayValue('patch_cells')
        if(rows%size()==size(source%members))then
            allocate(counts(source%kind,size(source%members)),source=0)
            do i=0,rows%size()-1
                if(rows%typeAt(i)/=QJsonValue_Object)exit
                row=rows%objectAt(i)
                id=bounded_integer(row,'patch_id',1,huge(1))
                column=findloc(source%members,id,dim=1)
                valid=only_fields(row,'|patch_id|cells|').and.column>0
                if(valid)valid=all(counts(:,column)==0)
                if(valid)valid=integer_array_arg(row,'cells',values,source%kind,source%kind,huge(1)-1)
                call row%delete()
                if(.not.valid)exit
                counts(:,column)=values
            end do
            valid=all(counts>0)
        end if
        call rows%delete()
        if(valid)error=''
    end function read_member_cell_counts

end module studio_mcp_geometry_edits
