module studio_mcp_selection
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use qt, only: QJsonArray, QJsonObject, QJsonValue_Array, QJsonValue_Double, QJsonValue_Object, QJsonValue_String
    use studio_mcp_arguments, only: bounded_integer, only_fields
    use studio_mcp_inspection, only: describe_selection, patch_kind
    use studio_multipatch_collection, only: count_multipatches
    use studio_scene, only: project_patch_exists
    use studio_selection, only: selection_target_t, selection_patch, selection_assembly, selection_control, selection_mesh, &
        selection_replace, selection_add, selection_toggle
    use studio_single_patch, only: single_patch_controls, read_single_patch_controls
    use studio_multipatch_controls, only: multipatch_controls_t, read_document_multipatch_controls
    use studio_canvas_input, only: change_shared_selection, picking_shared_controls
    use studio_commands, only: pick_controls_command
    use studio_state, only: studio_t, studio_mesh_valid, scene_mesh_drawn, dp, area_box, area_lasso
    use studio_ui, only: apply_workspace_selection_batch, invoke_command
    use studio_viewport_scene, only: edge_total, pick_face, pick_edge, pick_vertex
    use studio_viewport_selection, only: collect_viewport_region
    implicit none(type, external)
    private
    public :: set_selection, select_region

contains

    subroutine set_selection(s,args,result,error)
        type(studio_t), intent(inout), target :: s
        type(QJsonObject), intent(in) :: args
        type(QJsonObject), intent(inout) :: result
        character(:), allocatable, intent(out) :: error
        type(selection_target_t), allocatable :: targets(:)
        type(QJsonArray) :: entries
        type(QJsonObject) :: entry
        integer :: i,mode
        logical :: valid

        error='Supply items and an optional selection mode: replace, add, or toggle.'
        if(.not.only_fields(args,'|items|mode|'))return
        mode=read_mode(args)
        if(mode<0.or.args%valueType('items')/=QJsonValue_Array)return
        entries=args%arrayValue('items')
        allocate(targets(entries%size()))
        do i=1,size(targets)
            valid=entries%typeAt(i-1)==QJsonValue_Object
            if(valid)then
                entry=entries%objectAt(i-1)
                valid=read_target(s,entry,targets(i))
                call entry%delete()
            end if
            if(valid)cycle
            call entries%delete()
            error='A selection item identifies an unavailable object or subfeature. Read the scene before selecting.'
            return
        end do
        call entries%delete()
        call apply_mcp_selection(s,targets,mode,error)
        if(len(error)==0)call describe_selection(s,result)
    end subroutine set_selection

    subroutine apply_mcp_selection(s,targets,mode,error)
        type(studio_t), intent(inout), target :: s
        type(selection_target_t), intent(in) :: targets(:)
        integer, intent(in) :: mode
        character(:), allocatable, intent(out) :: error
        type(multipatch_controls_t) :: controls
        integer :: group

        error=''
        if(.not.any(targets%category==selection_assembly.and.targets%item>0))then
            call apply_workspace_selection_batch(s,targets,mode)
            return
        end if
        group=targets(1)%id
        error='Select shared controls from one assembly per request.'
        if(any(targets%category/=selection_assembly.or.targets%id/=group.or.targets%item==0))return
        error='Use replace to select shared controls in another assembly.'
        if(mode/=selection_replace.and.s%selected_multipatch/=group)return
        error='Show CAD geometry before selecting shared controls.'
        if(s%mesh_only)return
        call read_document_multipatch_controls(s%document,group,controls)
        error='A shared control ID is unavailable. Inspect the assembly controls for their exact identifiers.'
        if(.not.controls%valid)return
        if(any(targets%item>size(controls%weights)))return
        if(s%selected_multipatch/=group.or.s%selection%count()>1)then
            call apply_workspace_selection_batch(s,[selection_target_t(selection_assembly,group)],selection_replace)
        end if
        if(.not.picking_shared_controls(s).or..not.s%layers(4))call invoke_command(s,pick_controls_command)
        call change_shared_selection(s,targets%item,mode)
        error=''
    end subroutine apply_mcp_selection

    logical function read_target(s,entry,target) result(valid)
        type(studio_t), intent(in) :: s
        type(QJsonObject), intent(in) :: entry
        type(selection_target_t), intent(out) :: target
        type(single_patch_controls) :: controls
        character(:), allocatable :: kind,feature
        integer :: id,rank,item,part,tag,limit

        valid=.false.
        if(.not.only_fields(entry,'|kind|id|feature|index|part|'))return
        id=bounded_integer(entry,'id',1,huge(1))
        if(id<1)return
        kind=entry%stringValue('kind')
        feature=entry%stringValue('feature')
        if(entry%contains('feature').and.len(feature)==0)return
        item=bounded_integer(entry,'index',1,huge(1),fallback=0)
        part=bounded_integer(entry,'part',1,huge(1),fallback=0)
        if(min(item,part)<0)return
        if(len(feature)==0.and.(item>0.or.part>0))return
        select case(kind)
        case('assembly')
            if(id>count_multipatches(s%document%multipatches).or.part>0)return
            if(len(feature)>0.and.(feature/='control'.or.item<1))return
            target=selection_target_t(selection_assembly,id,item=item)
        case('mesh')
            if(.not.studio_mesh_valid(s,id))return
            target=selection_target_t(selection_mesh,id)
            if(len(feature)>0)then
                if(feature/='cell'.or.min(item,part)<1.or..not.scene_mesh_drawn(s,id))return
                associate(mesh=>s%document%meshes(id))
                    if(part>size(mesh%parts))return
                    if(.not.allocated(mesh%parts(part)%cells))return
                    if(item>size(mesh%parts(part)%cells,1))return
                end associate
                target%part=part
                target%item=item
            end if
        case default
            rank=patch_kind(kind)
            if(rank==0.or.part>0)return
            if(.not.project_patch_exists(s%document%show_templates,rank,id,patches=s%document%patches))return
            target=selection_target_t(selection_patch,id,rank)
            if(len(feature)>0)then
                tag=0
                select case(feature)
                case('control')
                    call read_single_patch_controls(rank,controls,variant=id,catalog=s%document%patches)
                    if(.not.allocated(controls%points))return
                    limit=size(controls%points,2)
                    target%category=selection_control
                case('face')
                    if(rank==1)return
                    tag=pick_face
                    limit=merge(6,1,rank==3)
                case('edge')
                    tag=pick_edge
                    limit=edge_total(rank)
                case('vertex')
                    tag=pick_vertex
                    limit=2**rank
                case default
                    return
                end select
                if(item<1.or.item>limit)return
                target%item=item
                target%feature=tag
            end if
        end select
        valid=.true.
    end function read_target

    subroutine select_region(s,args,result,error)
        type(studio_t), intent(inout), target :: s
        type(QJsonObject), intent(in) :: args
        type(QJsonObject), intent(inout) :: result
        character(:), allocatable, intent(out) :: error
        type(selection_target_t), allocatable :: targets(:)
        type(QJsonArray) :: entries,row
        real(dp), allocatable :: points(:,:)
        integer :: i,mode,tool,count
        logical :: valid

        error='Supply box or lasso with canvas-relative pixel points and an optional selection mode.'
        if(.not.only_fields(args,'|shape|points|mode|'))return
        mode=read_mode(args)
        if(mode<0.or.args%valueType('points')/=QJsonValue_Array)return
        select case(args%stringValue('shape'))
        case('box'); tool=area_box
        case('lasso'); tool=area_lasso
        case default; return
        end select
        entries=args%arrayValue('points')
        count=entries%size()
        if((tool==area_box.and.count/=2).or.(tool==area_lasso.and.count<3))then
            call entries%delete()
            error='Box needs two diagonal points. Lasso needs at least three points.'
            return
        end if
        allocate(points(2,count))
        do i=1,count
            valid=entries%typeAt(i-1)==QJsonValue_Array
            if(valid)then
                row=entries%arrayAt(i-1)
                valid=row%size()==2.and.row%typeAt(0)==QJsonValue_Double.and.row%typeAt(1)==QJsonValue_Double
                if(valid)points(:,i)=row%toDoubles()
                call row%delete()
            end if
            if(.not.valid)exit
            valid=all(ieee_is_finite(points(:,i))).and.all(points(:,i)>=0._dp).and. &
                all(points(:,i)<=real([s%surface%width(),s%surface%height()],dp))
            if(.not.valid)exit
        end do
        call entries%delete()
        if(.not.valid)then
            error='Each region point must be an XY position inside the canvas.'
            return
        end if
        if(tool==area_box.and.any(abs(points(:,2)-points(:,1))<1._dp))then
            error='A selection box must span at least one pixel on both axes.'
            return
        end if
        call collect_viewport_region(s,points,tool,targets)
        if(picking_shared_controls(s))then
            call change_shared_selection(s,targets%item,mode)
        else
            call apply_workspace_selection_batch(s,targets,mode)
        end if
        call describe_selection(s,result)
        error=''
    end subroutine select_region

    integer function read_mode(args) result(mode)
        type(QJsonObject), intent(in) :: args

        mode=-1
        if(args%contains('mode').and.args%valueType('mode')/=QJsonValue_String)return
        select case(args%stringValue('mode','replace'))
        case('replace'); mode=selection_replace
        case('add'); mode=selection_add
        case('toggle'); mode=selection_toggle
        case default; mode=-1
        end select
    end function read_mode

end module studio_mcp_selection
