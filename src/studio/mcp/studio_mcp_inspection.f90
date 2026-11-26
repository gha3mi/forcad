module studio_mcp_inspection
    use, intrinsic :: iso_fortran_env, only: int64
    use qt, only: QJsonObject, QJsonArray
    use studio_diagnostics, only: studio_diagnostic_summary
    use studio_edit_history, only: document_undo_stamp, document_can_redo
    use studio_mcp_arguments, only: only_fields, bounded_integer
    use studio_multipatch_collection, only: count_multipatches
    use studio_multipatch_controls, only: multipatch_controls_t, read_document_multipatch_controls
    use studio_geometry_operations, only: read_geometry_target
    use studio_measurements, only: measure_document_selection, measurement_t, measurement_label, measurement_unit, measurement_warning
    use studio_scene, only: project_patch_exists, project_geometry_targets, studio_scene_is_visible, get_studio_material
    use studio_selection, only: selection_target_t, selection_patch, selection_assembly, selection_control, selection_mesh
    use studio_single_patch, only: single_patch_scene, single_patch_controls, single_patch_name, single_patch_revision, &
        read_single_patch_scene, read_single_patch_controls
    use studio_state, only: studio_t, dp, studio_mesh_count, studio_mesh_valid, studio_mesh_source_group
    use studio_units, only: unit_symbol
    use studio_viewport_scene, only: pick_face, pick_edge, pick_vertex
    use studio_workspace_selection, only: primary_workspace_target, selection_measurement_targets
    implicit none(type, external)
    private

    public :: inspect_scene, inspect_patch, inspect_assembly, measure_geometry, describe_selection, describe_probe
    public :: patch_identity, assembly_identity, mesh_identity, active_rank, patch_kind

    integer, parameter, public :: default_scene_page=64, max_scene_page=128
    integer, parameter :: inline_member_count=16
    character(8), parameter :: rank_names(3)=[character(8)::'curve','surface','volume']

contains

    subroutine inspect_scene(s,result,offset,limit)
        type(studio_t),intent(in) :: s
        type(QJsonObject),intent(inout) :: result
        integer,optional,intent(in) :: offset,limit
        type(QJsonObject) :: item
        type(QJsonArray) :: objects,assemblies,meshes
        integer,allocatable :: targets(:,:)
        integer :: i,kind,id,first,page,object_count,assembly_count,mesh_count,total,last
        character(:),allocatable :: diagnostic

        first=0
        page=default_scene_page
        if(present(offset))first=offset
        if(present(limit))page=limit
        object_count=0
        assembly_count=0
        mesh_count=0
        call result%insert('project_open',.not.s%home_active)
        call result%insert('unit',unit_name(s%document%model_unit))
        call result%insert('document_revision',revision_text(s%document%document_revision))
        objects=QJsonArray()
        assemblies=QJsonArray()
        meshes=QJsonArray()
        if(.not.s%home_active)then
            targets=project_geometry_targets(s%document%show_templates,patches=s%document%patches)
            object_count=size(targets,2)
            assembly_count=count_multipatches(s%document%multipatches)
            mesh_count=studio_mesh_count(s)
            do i=min(first,object_count)+1,page_end(first,page,object_count)
                kind=targets(1,i)
                id=targets(2,i)
                item=patch_identity(s,kind,id)
                call objects%append(item)
                call item%delete()
            end do
            do id=min(first,assembly_count)+1,page_end(first,page,assembly_count)
                item=assembly_identity(s,id,summary=.true.)
                call assemblies%append(item)
                call item%delete()
            end do
            do id=min(first,mesh_count)+1,page_end(first,page,mesh_count)
                item=mesh_identity(s,id,summary=.true.)
                call meshes%append(item)
                call item%delete()
            end do
        end if
        call result%insert('objects',objects)
        call result%insert('assemblies',assemblies)
        call result%insert('meshes',meshes)
        call objects%delete()
        call assemblies%delete()
        call meshes%delete()
        call describe_selection(s,result,offset=first,limit=page)
        total=max(object_count,assembly_count,mesh_count,result%intValue('selected_count'))
        last=page_end(first,page,total)
        call result%insert('object_count',object_count)
        call result%insert('assembly_count',assembly_count)
        call result%insert('mesh_count',mesh_count)
        call result%insert('offset',first)
        call result%insert('next_offset',last)
        call result%insert('has_more',last<total)
        if(.not.s%home_active)then
            kind=active_rank(s)
            if(kind>0)call describe_probe(s,kind,result)
            call result%insert('undo_available',document_undo_stamp(s%document)>0_int64)
            call result%insert('redo_available',document_can_redo(s%document))
            diagnostic=studio_diagnostic_summary(s)
            if(len(diagnostic)>0)call result%insert('diagnostic',diagnostic)
        end if
    end subroutine inspect_scene

    integer function page_end(offset,limit,total) result(last)
        integer,intent(in) :: offset,limit,total
        integer :: first

        first=min(offset,total)
        last=first+min(limit,total-first)
    end function page_end

    subroutine describe_selection(s,result,offset,limit)
        type(studio_t),intent(in) :: s
        type(QJsonObject),intent(inout) :: result
        integer,optional,intent(in) :: offset,limit
        type(QJsonObject) :: selection
        type(QJsonArray) :: selected
        type(selection_target_t) :: target
        integer :: i,first,page,total
        logical :: shared_controls,paged

        first=0
        page=huge(0)
        total=0
        paged=present(limit)
        if(present(offset))first=offset
        if(present(limit))page=limit
        selected=QJsonArray()
        if(.not.s%home_active)then
            shared_controls=s%selected_multipatch>0.and.s%shared_controls%valid.and.s%selected_control>0
            if(shared_controls.and.allocated(s%selected_shared_controls))then
                total=size(s%selected_shared_controls)
                do i=min(first,total)+1,page_end(first,page,total)
                    target=selection_target_t(selection_assembly,s%selected_multipatch,item=s%selected_shared_controls(i))
                    call append_selection_target(selected,target)
                end do
            else if(shared_controls)then
                total=1
                target=selection_target_t(selection_assembly,s%selected_multipatch,item=s%selected_control)
                if(first==0)call append_selection_target(selected,target)
            else if(allocated(s%selection%items))then
                total=s%selection%count()
                do i=min(first,total)+1,page_end(first,page,total)
                    call append_selection_target(selected,s%selection%items(i))
                end do
            else
                target=primary_workspace_target(s)
                if(target%id>0)then
                    total=1
                    if(first==0)call append_selection_target(selected,target)
                end if
            end if
        end if
        if(paged)call result%insert('selected_count',total)
        call result%insert('selected_items',selected)
        call selected%delete()
        if(s%home_active)then
            call result%insertNull('selection')
            return
        end if
        if(s%selected_mesh>0)then
            selection=mesh_identity(s,s%selected_mesh,summary=paged)
        else if(s%selected_multipatch>0)then
            selection=assembly_identity(s,s%selected_multipatch,summary=paged)
        else if(s%live_scene%valid)then
            selection=patch_identity(s,s%entity_kind,s%live_scene%variant)
        else
            call result%insertNull('selection')
            return
        end if
        call result%insert('selection',selection)
        call selection%delete()
    end subroutine describe_selection

    subroutine append_selection_target(items,target)
        type(QJsonArray), intent(inout) :: items
        type(selection_target_t), intent(in) :: target
        type(QJsonObject) :: item
        character(:), allocatable :: kind,feature

        feature=''
        select case(target%category)
        case(selection_assembly)
            kind='assembly'
            if(target%item>0)feature='control'
        case(selection_mesh)
            kind='mesh'
            if(target%item>0)feature='cell'
        case(selection_patch,selection_control)
            if(target%rank<1.or.target%rank>size(rank_names))return
            kind=trim(rank_names(target%rank))
            if(target%category==selection_control)then
                feature='control'
            else
                select case(target%feature)
                case(pick_face); feature='face'
                case(pick_edge); feature='edge'
                case(pick_vertex); feature='vertex'
                case default; continue
                end select
            end if
        case default; return
        end select
        item=QJsonObject()
        call item%insert('kind',kind)
        call item%insert('id',target%id)
        if(len(feature)>0)then
            call item%insert('feature',feature)
            call item%insert('index',target%item)
            if(target%part>0)call item%insert('part',target%part)
        end if
        call items%append(item)
        call item%delete()
    end subroutine append_selection_target

    function patch_identity(s,kind,id) result(item)
        type(studio_t),intent(in) :: s
        integer,intent(in) :: kind,id
        type(QJsonObject) :: item
        type(QJsonArray) :: color
        integer :: rgb(3),opacity

        item=QJsonObject()
        call item%insert('kind',trim(rank_names(kind)))
        call item%insert('id',id)
        call item%insert('name',single_patch_name(kind,id,catalog=s%document%patches))
        call item%insert('visible',studio_scene_is_visible(s%document%scene,kind,id,catalog=s%document%patches))
        call get_studio_material(s%document%scene,kind,id,rgb,opacity)
        color=QJsonArray()
        call color%appendInts(rgb)
        call item%insert('rgb',color)
        call item%insert('opacity',opacity)
        call color%delete()
        call item%insert('revision',revision_text(single_patch_revision(kind,id,catalog=s%document%patches)))
    end function patch_identity

    function assembly_identity(s,id,summary) result(item)
        type(studio_t),intent(in) :: s
        integer,intent(in) :: id
        logical,optional,intent(in) :: summary
        type(QJsonObject) :: item
        type(QJsonArray) :: members
        integer :: i,j,visible,count_members
        logical :: include_members

        item=QJsonObject()
        associate(group=>s%document%multipatches%groups(id))
            call item%insert('kind','assembly')
            call item%insert('id',id)
            call item%insert('name',group%definition%name)
            call item%insert('patch_kind',trim(rank_names(group%definition%kind)))
            call item%insert('valid',group%validation%valid)
            visible=count([(studio_scene_is_visible(s%document%scene,group%definition%kind,group%definition%members(i),catalog=s%document%patches), &
                i=1,size(group%definition%members))])
            call item%insert('visible',visible>0)
            call item%insert('all_members_visible',visible==size(group%definition%members))
            count_members=size(group%definition%members)
            call item%insert('member_count',count_members)
            include_members=.true.
            if(present(summary))include_members=.not.summary.or.count_members<=inline_member_count
            if(include_members)then
                members=QJsonArray()
                call members%appendInts(group%definition%members)
                call item%insert('members',members)
                call members%delete()
            else
                call item%insert('members_in_inspect_assembly',.true.)
            end if
            if(group%validation%valid.and.allocated(group%revisions))then
                do i=1,size(group%definition%members)
                    j=group%definition%members(i)
                    if(single_patch_revision(group%definition%kind,j,catalog=s%document%patches)==group%revisions(i))cycle
                    call item%insert('valid',.false.)
                    call item%insert('diagnostic','A member changed since assembly validation.')
                    return
                end do
            end if
            if(allocated(group%validation%message))call item%insert('diagnostic',group%validation%message)
        end associate
    end function assembly_identity

    subroutine inspect_assembly(s,args,result,error)
        type(studio_t),intent(in) :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        type(QJsonObject) :: item,join
        type(QJsonArray) :: joins,flip,derivatives
        integer :: id,offset,limit,total,last,i,control_offset,control_limit

        error='Supply an assembly ID and optional join or control page offsets and limits.'
        if(.not.only_fields(args,'|id|join_offset|join_limit|control_offset|control_limit|'))return
        id=bounded_integer(args,'id',1,count_multipatches(s%document%multipatches))
        offset=bounded_integer(args,'join_offset',0,huge(1),fallback=0)
        limit=bounded_integer(args,'join_limit',1,128,fallback=64)
        control_offset=bounded_integer(args,'control_offset',0,huge(1),fallback=0)
        control_limit=bounded_integer(args,'control_limit',1,128,fallback=64)
        if(id<1.or.min(offset,control_offset)<0.or.min(limit,control_limit)<1)return
        joins=QJsonArray()
        associate(definition=>s%document%multipatches%groups(id)%definition)
            total=0
            if(allocated(definition%joins))total=size(definition%joins)
            last=min(offset,total)+min(limit,total-min(offset,total))
            do i=min(offset,total)+1,last
                associate(spec=>definition%joins(i))
                    join=QJsonObject()
                    flip=QJsonArray()
                    call flip%append(spec%flip(1))
                    call flip%append(spec%flip(2))
                    call join%insert('patch_a',spec%patch_a)
                    call join%insert('side_a',spec%side_a)
                    call join%insert('patch_b',spec%patch_b)
                    call join%insert('side_b',spec%side_b)
                    call join%insert('continuity',spec%continuity)
                    call join%insert('geometric',spec%geometric)
                    call join%insert('reverse',spec%reverse)
                    call join%insert('swap',spec%swap)
                    call join%insert('flip',flip)
                    if(allocated(spec%reparameterization))then
                        derivatives=QJsonArray()
                        call derivatives%appendDoubles(spec%reparameterization)
                        call join%insert('reparameterization',derivatives)
                        call derivatives%delete()
                    end if
                    call joins%append(join)
                    call flip%delete()
                    call join%delete()
                end associate
            end do
        end associate
        item=assembly_identity(s,id)
        call result%insert('assembly',item)
        call result%insert('joins',joins)
        call result%insert('join_count',total)
        call result%insert('next_join_offset',last)
        call result%insert('joins_have_more',last<total)
        call item%delete()
        call joins%delete()
        error=''
        if(args%contains('control_offset').or.args%contains('control_limit'))then
            call describe_assembly_controls(s,id,control_offset,control_limit,result,error)
        end if
    end subroutine inspect_assembly

    subroutine describe_assembly_controls(s,id,offset,limit,result,error)
        type(studio_t),intent(in) :: s
        integer,intent(in) :: id,offset,limit
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        type(multipatch_controls_t) :: controls
        type(QJsonArray) :: items,position
        type(QJsonObject) :: item
        integer :: first,last,count,i

        call read_document_multipatch_controls(s%document,id,controls)
        error=controls%message
        if(.not.controls%valid)return
        count=size(controls%weights)
        first=min(offset,count)
        last=first+min(limit,count-first)
        items=QJsonArray()
        do i=first+1,last
            item=QJsonObject()
            position=QJsonArray()
            call position%appendDoubles(controls%points(:,i))
            call item%insert('point',i)
            call item%insert('position',position)
            call item%insert('weight',controls%weights(i))
            call item%insert('copies',controls%copies(i))
            call items%append(item)
            call position%delete()
            call item%delete()
        end do
        call result%insert('controls',items)
        call result%insert('control_count',count)
        call result%insert('control_offset',first)
        call result%insert('next_control_offset',last)
        call result%insert('controls_have_more',last<count)
        call result%insert('unit',unit_name(s%document%model_unit))
        call items%delete()
    end subroutine describe_assembly_controls

    subroutine measure_geometry(s,args,result,error)
        type(studio_t),intent(in) :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        type(selection_target_t),allocatable :: targets(:)
        character(:),allocatable :: kind_name,name
        integer,allocatable :: members(:)
        integer :: kind,rank,patch,assembly,id,i

        error='Open a project before measuring geometry.'
        if(s%home_active)return
        error='Measure the selection, or supply an available patch or assembly kind and ID.'
        if(.not.only_fields(args,'|kind|id|'))return
        name='Selection'
        if(args%isEmpty())then
            targets=selection_measurement_targets(s)
            if(size(targets)==1)name=single_patch_name(targets(1)%rank,targets(1)%id,catalog=s%document%patches)
        else
            kind_name=args%stringValue('kind')
            id=bounded_integer(args,'id',1,huge(1))
            if(id<1)return
            kind=0
            patch=0
            assembly=0
            if(kind_name=='assembly')then
                assembly=id
            else
                kind=patch_kind(kind_name)
                if(kind==0)return
                if(.not.project_patch_exists(s%document%show_templates,kind,id,patches=s%document%patches))return
                patch=id
            end if
            call read_geometry_target(s%document,kind,patch,assembly,rank,members,name)
            targets=[(selection_target_t(selection_patch,members(i),rank),i=1,size(members))]
        end if
        call measure_items(s,targets,name,result,error)
    end subroutine measure_geometry

    subroutine measure_items(s,targets,name,result,error)
        type(studio_t), intent(in) :: s
        type(selection_target_t), intent(in) :: targets(:)
        character(*), intent(in) :: name
        type(QJsonObject), intent(inout) :: result
        character(:), allocatable, intent(out) :: error
        type(QJsonArray) :: measurements
        type(QJsonObject) :: item
        type(measurement_t) :: values(3)
        integer :: rank

        error='Select a curve, edge, surface, face, volume, or assembly to measure.'
        if(size(targets)==0)return
        call measure_document_selection(s%document,targets,values)
        if(all(values%count==0))return
        do rank=1,3
            if(values(rank)%count==0)cycle
            error=values(rank)%error
            if(len(error)>0)return
        end do
        measurements=QJsonArray()
        do rank=1,3
            if(values(rank)%count==0)cycle
            item=QJsonObject()
            call write_measurement(s,rank,values(rank),item)
            call measurements%append(item)
            call item%delete()
        end do
        call result%insert('name',name)
        call result%insert('measurements',measurements)
        call measurements%delete()
        if(count(values%count>0)==1)then
            rank=findloc(values%count>0,.true.,dim=1)
            call write_measurement(s,rank,values(rank),result)
        end if
    end subroutine measure_items

    subroutine write_measurement(s,rank,value,result)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: rank
        type(measurement_t), intent(in) :: value
        type(QJsonObject), intent(inout) :: result
        character(:),allocatable :: warning

        call result%insert('quantity',measurement_label(rank))
        call result%insert('value',value%value)
        call result%insert('unit',measurement_unit(s%document,rank))
        warning=measurement_warning(value)
        if(len(warning)>0)call result%insert('warning',warning)
    end subroutine write_measurement

    function mesh_identity(s,id,summary) result(item)
        type(studio_t),intent(in) :: s
        integer,intent(in) :: id
        logical,optional,intent(in) :: summary
        type(QJsonObject) :: item
        type(QJsonArray) :: counts,parts
        type(QJsonObject) :: patch
        integer :: part
        logical :: include_parts

        item=QJsonObject()
        associate(mesh=>s%document%meshes(id))
            call item%insert('kind','mesh')
            call item%insert('id',id)
            call item%insert('name',mesh%name)
            call item%insert('visible',mesh%visible)
            call item%insert('valid',studio_mesh_valid(s,id))
            call item%insert('conforming',allocated(mesh%connected_cells))
            call item%insert('nodes',mesh%nodes)
            call item%insert('elements',mesh%elements)
            call item%insert('source_patch',mesh%variant)
            call item%insert('source_assembly',studio_mesh_source_group(s,id))
            call item%insert('source_kind',trim(rank_names(mesh%kind)))
            counts=QJsonArray()
            call counts%appendInts(mesh%samples(:mesh%kind)-1)
            call item%insert('cells',counts)
            call counts%delete()
            call item%insert('per_patch_sizing',mesh%per_patch_sizing)
            if(mesh%target_size>0._dp)then
                call item%insert('target_size',real(mesh%target_size,dp))
                call item%insert('synchronize',mesh%synchronize)
            end if
            if(mesh%per_patch_sizing.or.any([(any(mesh%parts(part)%samples(:mesh%kind)/=mesh%samples(:mesh%kind)),part=1,size(mesh%parts))]))then
                include_parts=.true.
                if(present(summary))include_parts=.not.summary.or.size(mesh%parts)<=8
                if(.not.include_parts)then
                    call item%insert('patch_cell_count',size(mesh%parts))
                    call item%insert('patch_cells_in_select',.true.)
                    return
                end if
                parts=QJsonArray()
                do part=1,size(mesh%parts)
                    patch=QJsonObject()
                    counts=QJsonArray()
                    call counts%appendInts(mesh%parts(part)%samples(:mesh%kind)-1)
                    call patch%insert('patch_id',mesh%parts(part)%variant)
                    call patch%insert('cells',counts)
                    call parts%append(patch)
                    call counts%delete()
                    call patch%delete()
                end do
                call item%insert('patch_cells',parts)
                call parts%delete()
            end if
        end associate
    end function mesh_identity

    subroutine inspect_patch(s,args,result,error)
        type(studio_t),intent(in) :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        type(QJsonObject) :: item,control
        type(QJsonArray) :: axes,knots,points,coordinates,degree
        type(single_patch_scene) :: scene
        type(single_patch_controls) :: controls
        integer :: kind,id,offset,limit,i,upper,knot_offset,knot_limit,knot_count,lo,hi

        error='Use a curve, surface or volume ID with optional control and knot page offsets and limits.'
        if(.not.only_fields(args, '|kind|id|offset|limit|knot_offset|knot_limit|'))return
        if(s%home_active)return
        kind=patch_kind(args%stringValue('kind'))
        id=bounded_integer(args,'id',1,huge(1))
        offset=bounded_integer(args,'offset',0,huge(1),fallback=0)
        limit=bounded_integer(args,'limit',1,128,fallback=64)
        knot_offset=bounded_integer(args,'knot_offset',0,huge(1),fallback=0)
        knot_limit=bounded_integer(args,'knot_limit',1,128,fallback=128)
        if(min(kind,id,limit,knot_limit)<1.or.min(offset,knot_offset)<0)return
        error='That patch is unavailable in this project.'
        if(.not.project_patch_exists(s%document%show_templates,kind,id,patches=s%document%patches))return
        call read_single_patch_scene(kind,scene,variant=id,catalog=s%document%patches)
        if(.not.scene%valid)then
            if(allocated(scene%diagnostic))error=scene%diagnostic
            return
        end if
        item=patch_identity(s,kind,id)
        call result%insert('patch',item)
        call item%delete()
        call result%insert('unit',unit_name(s%document%model_unit))
        degree=QJsonArray()
        call degree%appendInts(scene%degree(:kind))
        call result%insert('degree',degree)
        call degree%delete()
        axes=QJsonArray()
        do i=1,kind
            item=QJsonObject()
            call item%insert('min',scene%bounds(1,i))
            call item%insert('max',scene%bounds(2,i))
            call item%insert('control_count',scene%control_counts(i))
            knot_count=size(scene%axes(i)%knots)
            lo=min(knot_offset,knot_count)
            hi=lo+min(knot_limit,knot_count-lo)
            call item%insert('knot_count',knot_count)
            call item%insert('knot_offset',lo)
            call item%insert('next_knot_offset',hi)
            call item%insert('knots_have_more',hi<knot_count)
            knots=QJsonArray()
            call knots%appendDoubles(scene%axes(i)%knots(lo+1:hi))
            call item%insert('knots',knots)
            call knots%delete()
            call axes%append(item)
            call item%delete()
        end do
        call result%insert('axes',axes)
        call axes%delete()
        ! Scene controls are in native tensor order. Read the weights for the same patch by ID.
        call read_single_patch_controls(kind,controls,variant=id,catalog=s%document%patches)
        if(.not.controls%valid)then
            error='ForCAD could not read the patch controls.'
            return
        end if
        points=QJsonArray()
        upper=min(offset,size(controls%weights))+min(limit,max(0,size(controls%weights)-offset))
        do i=min(offset,size(controls%weights))+1,upper
            control=QJsonObject()
            coordinates=QJsonArray()
            call coordinates%appendDoubles(controls%points(:,i))
            call control%insert('point',i)
            call control%insert('position',coordinates)
            call control%insert('weight',controls%weights(i))
            call points%append(control)
            call coordinates%delete()
            call control%delete()
        end do
        call result%insert('controls',points)
        call result%insert('control_count',size(controls%weights))
        call result%insert('offset',offset)
        call result%insert('next_offset',upper)
        call result%insert('has_more',upper<size(controls%weights))
        call points%delete()
        error=''
    end subroutine inspect_patch

    subroutine describe_probe(s,kind,result)
        type(studio_t),intent(in) :: s
        integer,intent(in) :: kind
        type(QJsonObject),intent(inout) :: result
        type(QJsonObject) :: item
        type(QJsonArray) :: values
        real(dp) :: normalized(3)

        item=QJsonObject()
        call item%insert('valid',s%live_probe%valid)
        normalized=[s%u,s%v,s%w]
        values=QJsonArray()
        call values%appendDoubles(normalized(:kind))
        call item%insert('normalized',values)
        call values%delete()
        if(s%live_probe%valid)then
            values=QJsonArray()
            call values%appendDoubles(s%live_probe%parameters(:kind))
            call item%insert('parameters',values)
            call values%delete()
            values=QJsonArray()
            call values%appendDoubles(s%live_probe%point)
            call item%insert('position',values)
            call values%delete()
        else if(allocated(s%live_probe%message))then
            call item%insert('diagnostic',s%live_probe%message)
        end if
        call result%insert('probe',item)
        call item%delete()
    end subroutine describe_probe

    integer function active_rank(s) result(kind)
        type(studio_t),intent(in) :: s

        kind=0
        if(s%selected_multipatch==0.and.s%selected_mesh==0.and.s%live_scene%valid)kind=s%entity_kind
    end function active_rank

    integer function patch_kind(name) result(kind)
        character(*),intent(in) :: name
        kind=findloc(rank_names,name,dim=1)
    end function patch_kind

    function unit_name(unit) result(name)
        integer,intent(in) :: unit
        character(:),allocatable :: name

        name=trim(unit_symbol(unit))
        if(name=='')name='unitless'
    end function unit_name

    function revision_text(revision) result(text)
        integer(int64),intent(in) :: revision
        character(:),allocatable :: text
        character(24) :: buffer

        write(buffer,'(i0)')revision
        text=trim(buffer)
    end function revision_text

end module studio_mcp_inspection
