module studio_ai_cursor
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_associated
    use, intrinsic :: iso_fortran_env, only: int64
    use qt, only: QBrush, QColor, QFont, QFontMetrics, QJsonArray, QJsonObject, QJsonValue_Array, &
        QJsonValue_Object, QPainter, QPainter_Antialiasing, QPen, Qt_DashLine, Qt_ElideRight, Qt_NoBrush
    use studio_state, only: studio_t, dp, scene_patch_drawn, scene_mesh_drawn, studio_mesh_valid
    use studio_multipatch_collection, only: count_multipatches
    use studio_multipatch_controls, only: multipatch_controls_t,read_document_multipatch_controls
    use studio_single_patch, only: single_patch_revision
    use studio_scene, only: studio_scene_control_position
    use studio_live_viewport, only: live_viewport_point, workspace_view_bounds, workspace_view_hit
    use studio_viewport_camera, only: in_drawable_view
    use studio_viewport_scene, only: corner_sample, edge_samples
    use studio_knot_interaction, only: domain_knot_anchor
    use studio_mcp_arguments, only: axis_number, number_arg, vector_arg
    use studio_selection, only: selection_target_t, selection_patch, selection_control, selection_assembly, selection_mesh
    use studio_workspace_selection, only: primary_workspace_target
    implicit none(type, external)
    private

    public :: hide_studio_ai_cursor, paint_studio_ai_cursor, studio_ai_cursor_position, studio_ai_cursor_note
    public :: set_studio_ai_step_cursor, set_studio_ai_tool_cursor, set_studio_ai_build_cursor
    integer, public, parameter :: ai_cursor_proposed=1, ai_cursor_executing=2, ai_cursor_applied=3
    integer, parameter :: no_cursor=0
    character(8), parameter :: target_kinds(5)=[character(8)::'curve','surface','volume','assembly','mesh']
    character(7), parameter :: feature_kinds(5)=[character(7)::'face','edge','vertex','control','cell']

    type :: cursor_target
        integer :: rank=0, variant=0, mesh=0, assembly=0, point=0, axis=0
        integer :: feature=0, feature_index=0, part=0, cell=0, step=0, total=0
        integer(int64) :: document_revision=-1_int64, target_revision=-1_int64
        real(dp) :: anchor(3)=0._dp, knot=0._dp, destination(3)=0._dp
        logical :: has_anchor=.false., has_destination=.false.
        character(:), allocatable :: activity
    end type cursor_target

    type :: cursor_marker
        type(cursor_target) :: target
        integer :: phase=no_cursor
        character(:), allocatable :: label
        type(c_ptr) :: canvas=c_null_ptr
    end type cursor_marker

    type(cursor_marker) :: marker,last_completed

contains

    subroutine hide_studio_ai_cursor(s)
        type(studio_t), intent(in) :: s

        if(c_associated(last_completed%canvas,s%surface%ptr))last_completed=cursor_marker()
        if(.not.c_associated(marker%canvas,s%surface%ptr))return
        marker=cursor_marker()
        call request_repaint(s)
    end subroutine hide_studio_ai_cursor

    subroutine restore_completed_cursor(s)
        type(studio_t), intent(in) :: s

        if(marker%phase/=ai_cursor_executing.or..not.c_associated(marker%canvas,s%surface%ptr))return
        marker=cursor_marker()
        if(c_associated(last_completed%canvas,s%surface%ptr))marker=last_completed
        call request_repaint(s)
    end subroutine restore_completed_cursor

    subroutine set_studio_ai_step_cursor(s,new_phase,step,total,label,point,axis,knot,to_position,to_knot,rank,variant)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: new_phase,step,total,rank,variant
        character(*), intent(in) :: label
        integer, optional, intent(in) :: point,axis
        real(dp), optional, intent(in) :: knot,to_position(3),to_knot
        type(cursor_target) :: target

        target%rank=rank
        target%variant=variant
        target%step=step
        target%total=total
        if(present(point))then
            target%point=point
            if(present(to_position))then
                target%destination=to_position
                target%has_destination=all(ieee_is_finite(to_position))
            end if
        else if(present(axis).and.present(knot))then
            target%axis=axis
            target%knot=knot
            if(present(to_knot))then
                target%destination(1)=to_knot
                target%has_destination=ieee_is_finite(to_knot)
            end if
        end if
        call show_cursor(s,new_phase,target,label)
    end subroutine set_studio_ai_step_cursor

    subroutine set_studio_ai_build_cursor(s,label,group,kind,variant,mesh_id)
        type(studio_t), intent(in) :: s
        character(*), intent(in) :: label
        integer, intent(in) :: group,kind,variant,mesh_id
        type(cursor_target) :: target

        if(scene_mesh_drawn(s,mesh_id))then
            target%mesh=mesh_id
        else
            target%assembly=group
            target%rank=kind
            target%variant=variant
        end if
        call show_cursor(s,ai_cursor_applied,target,label)
    end subroutine set_studio_ai_build_cursor

    subroutine set_studio_ai_tool_cursor(s,name,args,result)
        type(studio_t), intent(in) :: s
        character(*), intent(in) :: name
        type(QJsonObject), intent(in) :: args
        type(QJsonObject), optional, intent(in) :: result
        type(cursor_target) :: target
        character(:), allocatable :: label
        integer :: i,new_phase

        new_phase=ai_cursor_executing
        if(present(result))then
            if(.not.result%boolValue('ok',.false.))then
                call restore_completed_cursor(s)
                return
            end if
            new_phase=ai_cursor_applied
        end if
        call identify_tool_target(s,name,args,result,target)
        if(present(result))then
            select case(name)
            case('set_selection','select_region')
                if(target%rank==0.and.target%mesh==0.and.target%assembly==0)then
                    call hide_studio_ai_cursor(s)
                    return
                end if
            case default
            end select
        end if
        label=name
        select case(name)
        case('edit_patch','construct','transform_selection')
            label=args%stringValue('op',name)
            if(name=='edit_patch'.and.label=='control')then
                label='Edit control point'
                if(.not.args%contains('position').and.args%contains('weight'))label='Edit weight'
            end if
        case('execute')
            label=args%stringValue('command',name)
            if(present(result))label=result%stringValue('description',label)
        case('set_visibility')
            label='Hide object'
            if(args%boolValue('visible'))label='Show object'
        case default
        end select
        do i=1,len(label)
            if(index('_.-',label(i:i))>0)label(i:i)=' '
        end do
        if(len(label)>0)then
            i=index('abcdefghijklmnopqrstuvwxyz',label(1:1))
            if(i>0)label(1:1)=achar(iachar('A')+i-1)
        end if
        call show_cursor(s,new_phase,target,label,immediate=.not.present(result))
    end subroutine set_studio_ai_tool_cursor

    subroutine identify_tool_target(s,name,args,result,target)
        type(studio_t), intent(in) :: s
        character(*), intent(in) :: name
        type(QJsonObject), intent(in) :: args
        type(QJsonObject), optional, intent(in) :: result
        type(cursor_target), intent(out) :: target
        type(QJsonObject) :: item
        character(:), allocatable :: operation
        integer :: axis

        select case(name)
        case('inspect_patch','select','measure','set_visibility','rename_patch')
            if(name=='measure'.and.args%isEmpty())then
                call active_target(s,target)
            else
                call read_target(args,target)
            end if
        case('inspect_assembly')
            target%assembly=args%intValue('id')
        case('delete')
            if(present(result))then
                target%activity='Object removed'
            else
                call read_target(args,target)
            end if
        case('transform_selection','set_material')
            if(present(result))call read_selected_target(result,'selected_items',target)
            if(.not.has_identity(target))call active_target(s,target)
        case('transform_scene')
            call active_target(s,target)
            if(present(result).and.marker_in_document(s))then
                if(marker%phase==ai_cursor_executing.and.allocated(marker%target%activity))then
                    if(marker%target%activity=='Entire scene')target=marker%target
                end if
                target%has_anchor=.false.
            end if
            if(target%mesh>0)target=cursor_target()
            target%activity='Entire scene'
        case('edit_patch')
            if(present(result))then
                if(result%valueType('created')==QJsonValue_Object)then
                    call read_member(result,'created',target)
                    return
                end if
                call read_member(result,'target',target)
            end if
            operation=args%stringValue('op')
            if(.not.has_identity(target))then
                if(args%contains('target'))then
                    call read_member(args,'target',target)
                else if(operation=='control'.and.s%selected_multipatch>0.and.s%selection%count()<=1)then
                    target%assembly=s%selected_multipatch
                else
                    call active_patch(s,target)
                end if
            end if
            select case(operation)
            case('control')
                target%point=args%intValue('point')
                if(target%point<1)then
                    target=cursor_target()
                    return
                end if
                if(.not.present(result))target%has_destination=vector_arg(args,'position',target%destination)
            case('insert_knot','move_knot','remove_knot')
                axis=axis_number(args%stringValue('axis'),'UVW')
                if(axis<1.or.axis>target%rank.or..not.number_arg(args,'value',target%knot))then
                    target%variant=0
                    return
                end if
                target%axis=axis
                if(operation=='move_knot')then
                    if(present(result))then
                        target%knot=args%doubleValue('to')
                    else
                        target%has_destination=number_arg(args,'to',target%destination(1))
                    end if
                end if
            case default
            end select
        case('mesh','preview_mesh')
            if(present(result))call read_member(result,'mesh',target)
            if(.not.has_identity(target))then
                if(args%contains('target'))then
                    call read_member(args,'target',target)
                else
                    call active_target(s,target)
                end if
            end if
            target%activity='Mesh source'
        case('assemble')
            if(present(result))then
                call read_member(result,'assembly',target)
            else
                target%assembly=args%intValue('id')
                if(target%assembly==0)call source_patch(args,'members',target)
            end if
            target%activity='Multipatch'
        case('build_curved_volume')
            if(present(result))then
                call read_member(result,'mesh',target)
                if(target%mesh==0)call read_member(result,'assembly',target)
            end if
            target%activity='New multipatch'
        case('construct','build_model')
            if(present(result))then
                call read_member(result,'created',target)
                if(.not.has_identity(target))call read_member(result,'selection',target)
            else if(name=='construct')then
                call source_patch(args,'profiles',target)
            end if
            target%activity='New geometry'
        case('copy_patch')
            if(present(result))then
                call read_member(result,'created',target)
                if(.not.has_identity(target))then
                    target%rank=findloc(target_kinds(:3),result%stringValue('created_kind'),dim=1)
                    if(target%rank==0)target%rank=result%intValue('created_kind')
                    target%variant=result%intValue('created_id')
                end if
            else
                call read_target(args,target)
            end if
            target%activity='Copy geometry'
        case('import_file')
            if(present(result))then
                target%rank=findloc(target_kinds(:3),result%stringValue('selected_kind'),dim=1)
                target%variant=result%intValue('selected_id')
            end if
            target%activity='Project file'
        case('pick_mesh_cell')
            if(present(result))then
                target%mesh=result%intValue('mesh_id')
                target%part=result%intValue('part')
                target%cell=result%intValue('part_cell')
            else
                target%mesh=args%intValue('id',s%selected_mesh)
            end if
        case('set_probe')
            target%activity='Parameter Space'
            if(present(result))then
                if(result%valueType('probe')/=QJsonValue_Object)return
                item=result%objectValue('probe')
                if(item%boolValue('valid'))target%has_anchor=vector_arg(item,'position',target%anchor)
                call item%delete()
                if(target%has_anchor)call active_patch(s,target)
            end if
        case('set_selection')
            if(present(result))then
                call read_selected_target(result,'selected_items',target)
            else
                call read_selected_target(args,'items',target)
            end if
        case('select_region')
            if(present(result))call read_selected_target(result,'selected_items',target)
            target%activity='Viewport selection'
        case('set_camera')
            target%activity='3D view'
        case('set_mesh_view')
            target%activity='Mesh display'
        case('set_parameter_view')
            target%activity='Parameter Space'
        case('set_workspace')
            target%activity='Workspace layout'
        case('undo','redo')
            target%activity='Project history'
        case('save_project','open_project','new_project','export_file')
            target%activity='Project file'
        case('set_project_units')
            target%activity='Project units'
        case('execute')
            call command_target(s,args%stringValue('command'),result,target)
        case default
        end select
    end subroutine identify_tool_target

    subroutine command_target(s,command,result,target)
        type(studio_t), intent(in) :: s
        character(*), intent(in) :: command
        type(QJsonObject), optional, intent(in) :: result
        type(cursor_target), intent(inout) :: target

        if(len(command)==0)return
        select case(command)
        case('workspace.curve','workspace.surface','workspace.volume','mesh.select-source','model.duplicate')
            if(present(result))call read_member(result,'selection',target)
            if(.not.has_identity(target))call active_target(s,target)
        case('camera.fit-selection','scene.reveal','display.isolate','display.translucency','display.toggle-visibility', &
             'mesh.toggle-visibility','multipatch.validate','view.move-gizmo','view.rotate-gizmo','view.scale-gizmo')
            if(present(result))call read_selected_target(result,'selected_items',target)
            if(.not.has_identity(target))call active_target(s,target)
        case('patch.inspect','layout.domain')
            target%activity='Parameter Space'
        case('layout.basis')
            target%activity='Basis Plot'
        case('layout.inspector','inspector.show')
            target%activity='Inspector'
        case('layout.scene','scene.expand','scene.collapse')
            target%activity='Scene tree'
        case default
            target%activity='3D view'
            if(index(command,'layout.')==1)target%activity='Workspace layout'
            if(index(command,'mesh.')==1)target%activity='Mesh display'
            if(index(command,'basis.')==1)target%activity='Basis Plot'
            if(index(command,'domain.')==1.or.index(command,'knots.')==1.or.index(command,'probe.')==1)then
                target%activity='Parameter Space'
            end if
        end select
    end subroutine command_target

    subroutine source_patch(args,key,target)
        type(QJsonObject), intent(in) :: args
        character(*), intent(in) :: key
        type(cursor_target), intent(inout) :: target
        type(QJsonArray) :: items

        if(args%valueType(key)/=QJsonValue_Array)return
        items=args%arrayValue(key)
        if(items%size()>0)then
            target%rank=findloc(target_kinds(:3),args%stringValue('kind'),dim=1)
            target%variant=items%intAt(0)
        end if
        call items%delete()
    end subroutine source_patch

    subroutine read_selected_target(result,key,target)
        type(QJsonObject), intent(in) :: result
        character(*), intent(in) :: key
        type(cursor_target), intent(inout) :: target
        type(QJsonArray) :: items
        type(QJsonObject) :: item
        integer :: last

        if(result%valueType(key)/=QJsonValue_Array)return
        items=result%arrayValue(key)
        last=items%size()-1
        if(last>=0)then
            if(items%typeAt(last)==QJsonValue_Object)then
                item=items%objectAt(last)
                call read_target(item,target)
                call item%delete()
            end if
        end if
        call items%delete()
    end subroutine read_selected_target

    subroutine active_target(s,target)
        type(studio_t), intent(in) :: s
        type(cursor_target), intent(inout) :: target
        type(selection_target_t) :: selected
        integer :: count

        if(allocated(s%selection%items))then
            count=s%selection%count()
            if(count==0)return
            selected=s%selection%items(count)
        else
            selected=primary_workspace_target(s)
        end if
        select case(selected%category)
        case(selection_patch,selection_control)
            target%rank=selected%rank
            target%variant=selected%id
            if(selected%category==selection_control)then
                target%point=selected%item
            else
                target%feature=selected%feature
                target%feature_index=selected%item
            end if
        case(selection_mesh)
            target%mesh=selected%id
            target%part=selected%part
            target%cell=selected%item
        case(selection_assembly)
            target%assembly=selected%id
        case default
        end select
    end subroutine active_target

    logical function has_identity(target) result(identified)
        type(cursor_target), intent(in) :: target

        identified=(target%rank>0.and.target%variant>0).or.target%mesh>0.or.target%assembly>0
    end function has_identity

    subroutine active_patch(s,target)
        type(studio_t), intent(in) :: s
        type(cursor_target), intent(inout) :: target

        if(.not.s%live_scene%valid.or.s%selected_mesh/=0.or.s%selected_multipatch/=0)return
        target%rank=s%entity_kind
        target%variant=s%live_scene%variant
    end subroutine active_patch

    subroutine read_member(result,key,target)
        type(QJsonObject), intent(in) :: result
        character(*), intent(in) :: key
        type(cursor_target), intent(inout) :: target
        type(QJsonObject) :: item

        if(result%valueType(key)/=QJsonValue_Object)return
        item=result%objectValue(key)
        call read_target(item,target)
        call item%delete()
    end subroutine read_member

    subroutine read_target(item,target)
        type(QJsonObject), intent(in) :: item
        type(cursor_target), intent(inout) :: target
        integer :: category,id

        category=findloc(target_kinds,item%stringValue('kind'),dim=1)
        if(item%stringValue('kind')=='multipatch')category=4
        id=item%intValue('id')
        if(id<1)return
        target%feature=findloc(feature_kinds,item%stringValue('feature'),dim=1)
        target%feature_index=item%intValue('index')
        select case(category)
        case(1:3)
            target%rank=category
            target%variant=id
            if(target%feature==4)target%point=target%feature_index
        case(4)
            target%assembly=id
        case(5)
            target%mesh=id
            target%part=item%intValue('part')
            if(target%feature==5)target%cell=target%feature_index
        case default
        end select
    end subroutine read_target

    subroutine show_cursor(s,new_phase,target,label,immediate)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: new_phase
        type(cursor_target), intent(inout) :: target
        character(*), intent(in) :: label
        logical, optional, intent(in) :: immediate
        integer :: panel
        real(dp) :: x,y

        if(.not.resolve_target(s,target))then
            call restore_completed_cursor(s)
            return
        end if
        marker=cursor_marker(target,new_phase,label,s%surface%ptr)
        if(new_phase==ai_cursor_applied)last_completed=marker
        call request_repaint(s)
        if(present(immediate))then
            if(immediate)then
                if(studio_ai_cursor_position(s,panel,x,y))then
                    if(panel==5)then
                        call s%domain%repaint()
                        return
                    end if
                end if
                call s%surface%repaint()
            end if
        end if
    end subroutine show_cursor

    logical function resolve_target(s,target) result(valid)
        type(studio_t), intent(in) :: s
        type(cursor_target), intent(inout) :: target

        target%document_revision=s%document%document_revision
        valid=.false.
        if(.not.has_identity(target))then
            valid=allocated(target%activity)
            return
        end if
        if(target%mesh>0)then
            if(.not.studio_mesh_valid(s,target%mesh))return
            call mesh_anchor(s,target)
            if(.not.target%has_anchor)return
            target%target_revision=s%document%mesh_generation
        else
            if(target%assembly>0)then
                if(target%point>0)then
                    call assembly_control(s,target)
                else
                    call assembly_member(s,target)
                end if
            end if
            if(target%rank<1.or.target%variant<1)return
            target%target_revision=single_patch_revision(target%rank,target%variant,catalog=s%document%patches)
            if(target%target_revision<0)return
            call patch_anchor(s,target)
            if(.not.target%has_anchor.and.target%axis==0)then
                if(scene_patch_drawn(s,target%rank,target%variant))return
            end if
        end if
        valid=.true.
    end function resolve_target

    subroutine assembly_member(s,target)
        type(studio_t), intent(in) :: s
        type(cursor_target), intent(inout) :: target
        integer :: i

        target%variant=0
        if(target%assembly>count_multipatches(s%document%multipatches))return
        associate(group=>s%document%multipatches%groups(target%assembly)%definition)
            target%rank=group%kind
            if(.not.allocated(group%members))return
            if(size(group%members)==0)return
            target%variant=group%members(1)
            do i=1,size(group%members)
                if(.not.scene_patch_drawn(s,target%rank,group%members(i)))cycle
                target%variant=group%members(i)
                return
            end do
        end associate
    end subroutine assembly_member

    subroutine assembly_control(s,target)
        type(studio_t), intent(in) :: s
        type(cursor_target), intent(inout) :: target
        type(multipatch_controls_t) :: controls
        integer :: i,point,first,last,chosen

        target%variant=0
        call read_document_multipatch_controls(s%document,target%assembly,controls)
        if(.not.controls%valid)return
        target%rank=controls%kind
        chosen=0
        do i=1,size(controls%members)
            first=controls%offsets(i)+1
            last=controls%offsets(i+1)
            point=findloc(controls%map(first:last),target%point,dim=1)
            if(point<1)cycle
            if(target%variant==0.or.scene_patch_drawn(s,target%rank,controls%members(i)))then
                target%variant=controls%members(i)
                chosen=point
                if(scene_patch_drawn(s,target%rank,target%variant))exit
            end if
        end do
        if(target%variant>0)target%point=chosen
    end subroutine assembly_control

    subroutine patch_anchor(s,target)
        type(studio_t), intent(in) :: s
        type(cursor_target), intent(inout) :: target
        integer, allocatable :: nodes(:)
        integer :: sample,side,axis,counts(3),index(3),sides(2),n

        if(target%axis>0)return
        if(.not.scene_patch_drawn(s,target%rank,target%variant))then
            target%has_anchor=.false.
            return
        end if
        if(target%point>0)then
            target%has_anchor=studio_scene_control_position(s%document%scene,target%rank,target%variant,target%point,target%anchor, &
                catalog=s%document%patches)
            return
        end if
        if(target%has_anchor)return
        associate(patch=>s%document%scene%catalog(target%rank)%entries(target%variant)%patch)
            if(.not.patch%valid.or..not.allocated(patch%points))return
            counts=patch%sample_counts
            index=(counts+1)/2
            side=target%feature_index
            select case(target%feature)
            case(0)
                sample=sample_index(counts,index)
            case(1)
                if(side<1.or.target%rank<2)return
                if(target%rank==3)then
                    if(side>6)return
                    axis=(side+1)/2
                    index(axis)=merge(counts(axis),1,mod(side,2)==0)
                else if(side/=1)then
                    return
                end if
                sample=sample_index(counts,index)
            case(2)
                allocate(nodes(maxval(counts)))
                call edge_samples(counts,target%rank,side,nodes,n,sides)
                if(n<1)return
                sample=nodes((n+1)/2)
            case(3)
                if(side<1.or.side>2**target%rank)return
                sample=corner_sample(counts,target%rank,side)
            case default
                return
            end select
            if(sample<1.or.sample>size(patch%points,2))return
            target%anchor=patch%points(:,sample)
            target%has_anchor=all(ieee_is_finite(target%anchor))
        end associate
    end subroutine patch_anchor

    pure integer function sample_index(counts,index) result(node)
        integer, intent(in) :: counts(3),index(3)

        node=index(1)+(index(2)-1)*counts(1)+(index(3)-1)*product(counts(:2))
    end function sample_index

    subroutine mesh_anchor(s,target)
        type(studio_t), intent(in) :: s
        type(cursor_target), intent(inout) :: target
        integer :: part,node,dimensions

        if(.not.studio_mesh_valid(s,target%mesh))return
        associate(mesh=>s%document%meshes(target%mesh))
            if(.not.allocated(mesh%parts))return
            if(size(mesh%parts)==0)return
            part=target%part
            if(part==0.and.target%cell==0)part=(size(mesh%parts)+1)/2
            if(part<1.or.part>size(mesh%parts))return
            associate(data=>mesh%parts(part))
                if(.not.allocated(data%points))return
                dimensions=size(data%points,2)
                if(size(data%points,1)==0.or.dimensions<1.or.dimensions>size(target%anchor))return
                if(target%cell>0)then
                    if(.not.allocated(data%cells))return
                    if(target%cell>size(data%cells,1).or.size(data%cells,2)==0)return
                    associate(nodes=>data%cells(target%cell,:))
                        if(any(nodes<1).or.any(nodes>size(data%points,1)))return
                        target%anchor(:dimensions)=sum(data%points(nodes,:),dim=1)/real(size(nodes),dp)
                    end associate
                else
                    node=(size(data%points,1)+1)/2
                    target%anchor(:dimensions)=data%points(node,:)
                end if
                target%has_anchor=all(ieee_is_finite(target%anchor))
            end associate
        end associate
    end subroutine mesh_anchor

    subroutine request_repaint(s)
        type(studio_t), intent(in) :: s

        if(.not.s%surface%isNull())call s%surface%update()
        if(.not.s%domain%isNull())call s%domain%update()
    end subroutine request_repaint

    logical function studio_ai_cursor_position(s,panel,x,y) result(visible)
        type(studio_t), intent(in) :: s
        integer, intent(out) :: panel
        real(dp), intent(out) :: x,y
        type(cursor_target) :: target
        real(dp) :: bounds(4)

        visible=.false.
        panel=0
        x=0._dp
        y=0._dp
        if(.not.marker_in_document(s))return
        target=marker%target
        if(.not.has_identity(target))return
        if(.not.target_current(s,target))return
        if(target%axis>0)then
            if(s%domain%isNull().or..not.s%live_scene%valid)return
            if(.not.s%domain%isVisible().or.s%knot_axis/=target%axis)return
            if(s%entity_kind/=target%rank.or.s%live_scene%variant/=target%variant)return
            if(.not.domain_knot_anchor(s,target%axis,target%knot,x,y))return
            panel=5
            bounds=[0._dp,0._dp,real(s%domain%width(),dp),real(s%domain%height(),dp)]
        else
            if(.not.target_drawn(s,target))return
            if(.not.target%has_anchor)call patch_anchor(s,target)
            if(.not.target%has_anchor)return
            panel=4
            bounds=workspace_view_bounds(s)
            call live_viewport_point(s,target%anchor,x,y)
            if(.not.workspace_view_hit(s,x,y))return
        end if
        visible=within_view(x,y,bounds)
    end function studio_ai_cursor_position

    logical function marker_in_document(s) result(active)
        type(studio_t), intent(in) :: s

        active=marker%phase/=no_cursor.and.c_associated(marker%canvas,s%surface%ptr).and. &
            marker%target%document_revision==s%document%document_revision
    end function marker_in_document

    logical function target_current(s,target) result(current)
        type(studio_t), intent(in) :: s
        type(cursor_target), intent(in) :: target

        if(target%mesh>0)then
            current=target%target_revision==s%document%mesh_generation.and.studio_mesh_valid(s,target%mesh)
        else
            current=target%target_revision==single_patch_revision(target%rank,target%variant,catalog=s%document%patches)
        end if
    end function target_current

    logical function target_drawn(s,target) result(drawn)
        type(studio_t), intent(in) :: s
        type(cursor_target), intent(in) :: target

        if(target%mesh>0)then
            drawn=scene_mesh_drawn(s,target%mesh)
        else
            drawn=scene_patch_drawn(s,target%rank,target%variant)
        end if
    end function target_drawn

    function studio_ai_cursor_note(s) result(note)
        type(studio_t), intent(in) :: s
        character(:), allocatable :: note
        integer :: panel
        real(dp) :: x,y

        note=''
        if(.not.marker_in_document(s))return
        if(.not.has_identity(marker%target))then
            if(allocated(marker%target%activity))note=marker%target%activity
            return
        end if
        if(studio_ai_cursor_position(s,panel,x,y))return
        if(.not.target_current(s,marker%target))then
            note='Target changed'
        else if(marker%target%axis>0)then
            note='Target in Parameter Space'
        else if(.not.target_drawn(s,marker%target))then
            note='Target hidden'
        else
            note='Target outside view'
            if(in_drawable_view(real(s%surface%width(),dp),real(s%surface%height(),dp),x,y))then
                if(.not.workspace_view_hit(s,x,y))note='Target behind a panel'
            end if
        end if
    end function studio_ai_cursor_note

    logical function within_view(x,y,bounds) result(visible)
        real(dp), intent(in) :: x,y,bounds(4)

        visible=.false.
        if(.not.all(ieee_is_finite([x,y])))return
        visible=x>=bounds(1)+4._dp.and.x<=bounds(3)-4._dp.and.y>=bounds(2)+4._dp.and.y<=bounds(4)-4._dp
    end function within_view

    logical function workspace_badge_visible(s,x,y,w,h) result(visible)
        type(studio_t), intent(in) :: s
        real(dp), intent(in) :: x,y,w,h

        visible=workspace_view_hit(s,x,y).and.workspace_view_hit(s,x+w,y).and. &
            workspace_view_hit(s,x,y+h).and.workspace_view_hit(s,x+w,y+h)
    end function workspace_badge_visible

    subroutine place_workspace_badge(s,bounds,w,h,x,y)
        type(studio_t), intent(in) :: s
        real(dp), intent(in) :: bounds(4),w,h
        real(dp), intent(inout) :: x,y
        real(dp) :: candidates(2,4)
        integer :: i

        if(workspace_badge_visible(s,x,y,w,h))return
        candidates(:,1)=[bounds(1)+12._dp,bounds(2)+12._dp]
        candidates(:,2)=[bounds(1)+12._dp,bounds(4)-h-12._dp]
        candidates(:,3)=[bounds(3)-w-12._dp,bounds(2)+12._dp]
        candidates(:,4)=[bounds(3)-w-12._dp,bounds(4)-h-12._dp]
        do i=1,size(candidates,2)
            if(.not.workspace_badge_visible(s,candidates(1,i),candidates(2,i),w,h))cycle
            x=candidates(1,i)
            y=candidates(2,i)
            return
        end do
    end subroutine place_workspace_badge

    subroutine paint_studio_ai_cursor(painter,s,panel)
        type(QPainter), intent(inout) :: painter
        type(studio_t), intent(in) :: s
        integer, intent(in) :: panel
        type(QFontMetrics) :: metrics
        character(:), allocatable :: title,detail,note
        character(32) :: sequence
        integer :: target_panel,available,space_left,space_right
        real(dp) :: x,y,badge_x,badge_y,badge_width,badge_height,font_height,ascent,end_x,end_y,bounds(4)
        logical :: has_path,badge_to_right,anchored
        real(dp), parameter :: tip_x(7)=[0._dp,0._dp,4._dp,7._dp,11._dp,8._dp,14._dp]
        real(dp), parameter :: tip_y(7)=[0._dp,17._dp,13._dp,21._dp,19._dp,12._dp,11._dp]

        anchored=studio_ai_cursor_position(s,target_panel,x,y)
        if(anchored)then
            if(panel/=target_panel)return
        else
            if(panel/=4)return
            note=studio_ai_cursor_note(s)
            if(len(note)==0)return
        end if
        if(panel==4)then
            bounds=workspace_view_bounds(s)
        else
            bounds=[0._dp,0._dp,real(s%domain%width(),dp),real(s%domain%height(),dp)]
        end if
        select case(marker%phase)
        case(ai_cursor_proposed)
            title='AI · Review'
        case(ai_cursor_applied)
            title='AI · Completed'
        case default
            title='AI · Working'
        end select
        if(marker%target%step>0)then
            write(sequence,'(i0,a,i0)') marker%target%step,'/',marker%target%total
            title=title//' '//trim(sequence)
        end if
        has_path=.false.
        if(anchored.and.marker%target%has_destination)then
            if(panel==4)then
                call live_viewport_point(s,marker%target%destination,end_x,end_y)
                has_path=within_view(end_x,end_y,bounds).and.workspace_view_hit(s,end_x,end_y)
            else if(domain_knot_anchor(s,marker%target%axis,marker%target%destination(1),end_x,end_y))then
                has_path=within_view(end_x,end_y,bounds)
            end if
        end if
        detail=marker%label
        if(anchored)then
            space_left=int(x-bounds(1))-27
            space_right=int(bounds(3)-x)-34
            badge_to_right=space_right>=space_left
            available=max(32,max(space_left,space_right))
        else
            available=min(260,int(bounds(3)-bounds(1))-32)
            if(available<48)return
            detail=detail//' · '//note
        end if
        call painter%save()
        call painter%setRenderHint(QPainter_Antialiasing)
        call painter%setFont(QFont(family='Noto Sans',pointSize=9._dp))
        metrics=painter%fontMetrics()
        title=metrics%elidedText(title,Qt_ElideRight,available)
        detail=metrics%elidedText(detail,Qt_ElideRight,available)
        font_height=real(metrics%height(),dp)
        ascent=real(metrics%ascent(),dp)
        badge_width=real(max(metrics%horizontalAdvance(title),metrics%horizontalAdvance(detail)),dp)+12._dp
        badge_height=font_height+6._dp
        if(len(detail)>0)badge_height=badge_height+font_height
        call metrics%delete()
        if(anchored)then
            badge_x=merge(x+16._dp,x-9._dp-badge_width,badge_to_right)
            badge_x=max(bounds(1)+4._dp,min(badge_x,bounds(3)-badge_width-4._dp))
            badge_y=max(bounds(2)+4._dp,min(y+3._dp,bounds(4)-badge_height-4._dp))
        else
            badge_x=bounds(1)+12._dp
            badge_y=bounds(2)+12._dp
        end if
        if(panel==4)call place_workspace_badge(s,bounds,badge_width,badge_height,badge_x,badge_y)
        if(has_path)then
            call painter%setPen(QPen(QColor(195,150,250,200),1.2_dp,Qt_DashLine))
            call painter%drawLine(x,y,end_x,end_y)
            call painter%setPen(QPen(QColor(218,189,255),1.3_dp))
            call painter%setBrush(QBrush(style=Qt_NoBrush))
            call painter%drawEllipse(end_x-5._dp,end_y-5._dp,10._dp,10._dp)
        end if
        if(anchored)then
            call painter%setPen(QPen(QColor(213,172,254),1.3_dp))
            call painter%setBrush(QColor(160,113,240,44))
            call painter%drawEllipse(x-6._dp,y-6._dp,12._dp,12._dp)
            call painter%setPen(QPen(QColor(30,21,46),1._dp))
            call painter%setBrush(QColor(204,172,255))
            call painter%drawPolygon(x+tip_x,y+tip_y)
        end if
        call painter%setPen(QPen(QColor(175,129,237),1._dp))
        call painter%setBrush(QColor(37,29,56,228))
        call painter%drawRect(badge_x,badge_y,badge_width,badge_height)
        call painter%setPen(QPen(QColor(242,230,255),1._dp))
        call painter%drawText(badge_x+6._dp,badge_y+3._dp+ascent,title)
        if(len(detail)>0)then
            call painter%setPen(QPen(QColor(210,190,235),1._dp))
            call painter%drawText(badge_x+6._dp,badge_y+3._dp+font_height+ascent,detail)
        end if
        call painter%restore()
    end subroutine paint_studio_ai_cursor

end module studio_ai_cursor
