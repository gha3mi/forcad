module studio_mcp_construction
    use qt, only: QJsonArray, QJsonObject, QJsonValue_Array, QJsonValue_Object
    use studio_state, only: studio_t, dp, studio_mesh_count, append_studio_mesh
    use studio_scene, only: project_patch_exists
    use studio_single_patch, only: studio_native_patch, single_patch_count, active_patch_id, publish_native_patches
    use studio_model_plan, only: model_step, build_model_plan, model_plan_max_steps
    use studio_construction, only: construction_actions, publish_construction
    use studio_curved_volume_plan, only: curved_volume_level, build_curved_volume_plan, &
        curved_volume_max_levels, curved_volume_sectors
    use studio_forcad, only: primitive_actions
    use studio_multipatch, only: multipatch_join
    use studio_multipatch_collection, only: publish_multipatch_volume
    use studio_mesh, only: studio_mesh_t
    use studio_ui, only: on_assistant_model_created
    use studio_mcp_arguments, only: only_fields, bounded_integer, number_arg, vector_arg, integer_array_arg
    use studio_mcp_inspection, only: describe_selection, patch_identity, assembly_identity, mesh_identity, patch_kind
    use studio_presentation_text, only: integer_text
    implicit none(type, external)
    private

    public :: build_model, construct_existing, build_curved_volume

contains

    subroutine build_model(s,args,result,error)
        type(studio_t),intent(inout),target :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        type(QJsonArray) :: list,created
        type(QJsonObject) :: item
        type(model_step),allocatable :: steps(:)
        type(studio_native_patch),allocatable :: patches(:)
        character(96),allocatable :: names(:)
        integer :: i,kind,id,before(3),rank
        logical :: valid

        error='Supply 1 to 12 ordered modeling steps.'
        if(.not.only_fields(args,'|steps|').or.args%valueType('steps')/=QJsonValue_Array)return
        list=args%arrayValue('steps')
        if(list%size()<1.or.list%size()>model_plan_max_steps)then
            call list%delete()
            return
        end if
        allocate(steps(list%size()))
        do i=1,size(steps)
            valid=list%typeAt(i-1)==QJsonValue_Object
            if(valid)then
                item=list%objectAt(i-1)
                valid=parse_model_step(item,steps(i))
                call item%delete()
            end if
            if(valid)cycle
            call list%delete()
            error='Invalid modeling arguments for step '//integer_text(i)//'.'
            return
        end do
        call list%delete()
        call build_model_plan(steps,patches,names,error)
        if(len(error)>0)return
        before=[(single_patch_count(i,catalog=s%document%patches),i=1,3)]
        call publish_native_patches(patches,names,kind,id,error,operation='construction',catalog=s%document%patches)
        if(len(error)>0)return
        call on_assistant_model_created(s,0,kind,id,0)
        created=QJsonArray()
        do rank=1,3
            do i=before(rank)+1,single_patch_count(rank,catalog=s%document%patches)
                item=patch_identity(s,rank,i)
                call created%append(item)
                call item%delete()
            end do
        end do
        call result%insert('created',created)
        call created%delete()
        call describe_selection(s,result)
    end subroutine build_model

    subroutine construct_existing(s,args,result,error)
        type(studio_t),intent(inout),target :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        type(QJsonObject) :: item
        integer,allocatable :: profile_ids(:)
        integer :: operation,kind,i,created_kind,created_id,spine,degree,available
        real(dp) :: vector(3),center(3),angle
        character(:),allocatable :: name

        operation=findloc(construction_actions,args%stringValue('op'),dim=1)
        error='Choose extrude, revolve, sweep or loft.'
        if(operation==0)return
        name=trim(adjustl(args%stringValue('name')))
        error='Enter a name for the result.'
        if(len(name)==0.or.index(name,achar(0))/=0)return
        kind=patch_kind(args%stringValue('kind'))
        error='Choose existing curve or surface profiles by ID.'
        if(kind<1.or.kind>2)return
        available=single_patch_count(kind,catalog=s%document%patches)
        if(operation==4)then
            if(.not.integer_array_arg(args,'profiles',profile_ids,2,available,available))return
        else
            if(.not.integer_array_arg(args,'profiles',profile_ids,1,1,available))return
        end if
        if(.not.all([(project_patch_exists(s%document%show_templates,kind,profile_ids(i),patches=s%document%patches),i=1,size(profile_ids))]))return

        error='Supply the fields required for this construction.'
        select case(operation)
        case(1)
            if(.not.only_fields(args,'|op|kind|profiles|name|vector|'))return
            if(.not.vector_arg(args,'vector',vector))return
            call publish_construction(operation,kind,profile_ids,name,created_kind,error,vector=vector,catalog=s%document%patches)
        case(2)
            if(.not.only_fields(args,'|op|kind|profiles|name|vector|axis_point|angle_degrees|'))return
            if(.not.vector_arg(args,'vector',vector))return
            if(.not.vector_arg(args,'axis_point',center))return
            if(.not.number_arg(args,'angle_degrees',angle))return
            call publish_construction(operation,kind,profile_ids,name,created_kind,error, &
                vector=vector,axis_point=center,angle_degrees=angle,catalog=s%document%patches)
        case(3)
            if(.not.only_fields(args,'|op|kind|profiles|name|spine|'))return
            spine=bounded_integer(args,'spine',1,single_patch_count(1,catalog=s%document%patches))
            error='Choose an existing curve ID for the sweep path.'
            if(.not.project_patch_exists(s%document%show_templates,1,spine,patches=s%document%patches))return
            call publish_construction(operation,kind,profile_ids,name,created_kind,error,spine=spine,catalog=s%document%patches)
        case(4)
            if(.not.only_fields(args,'|op|kind|profiles|name|degree|'))return
            degree=bounded_integer(args,'degree',0,size(profile_ids)-1,fallback=0)
            error='Choose zero for automatic degree or a positive loft degree lower than the section count.'
            if(degree<0)return
            if(degree>0)then
                call publish_construction(operation,kind,profile_ids,name,created_kind,error,degree=degree,catalog=s%document%patches)
            else
                call publish_construction(operation,kind,profile_ids,name,created_kind,error,catalog=s%document%patches)
            end if
        case default
            return
        end select
        if(len(error)>0)return
        created_id=active_patch_id(created_kind,catalog=s%document%patches)
        call on_assistant_model_created(s,0,created_kind,created_id,0)
        item=patch_identity(s,created_kind,created_id)
        call result%insert('created',item)
        call item%delete()
        call describe_selection(s,result)
    end subroutine construct_existing

    subroutine build_curved_volume(s,args,result,error)
        type(studio_t),intent(inout),target :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        character(:),allocatable :: name
        type(curved_volume_level),allocatable :: levels(:)
        type(studio_native_patch),allocatable :: patches(:)
        character(96),allocatable :: names(:)
        type(multipatch_join),allocatable :: joins(:)
        type(studio_mesh_t) :: mesh
        type(QJsonArray) :: list
        type(QJsonObject) :: item
        integer,allocatable :: cells(:)
        integer :: i,group,mesh_id
        logical :: valid

        error='Supply a name, 2 to 12 axial levels, and optional mesh_cells.'
        if(.not.only_fields(args,'|name|levels|mesh_cells|'))return
        name=trim(args%stringValue('name'))
        if(len(name)==0.or.len(name)>96.or.index(name,achar(0))/=0)return
        if(args%valueType('levels')/=QJsonValue_Array)return
        if(args%contains('mesh_cells'))then
            error='Enter three positive integer mesh cell counts.'
            if(.not.integer_array_arg(args,'mesh_cells',cells,3,3,huge(1)-1))return
        end if
        list=args%arrayValue('levels')
        if(list%size()<2.or.list%size()>curved_volume_max_levels)then
            call list%delete()
            error='A curved volume requires 2 to 12 axial levels.'
            return
        end if
        allocate(levels(list%size()))
        do i=1,size(levels)
            valid=list%typeAt(i-1)==QJsonValue_Object
            if(valid)then
                item=list%objectAt(i-1)
                valid=parse_curved_level(item,levels(i),i==1)
                call item%delete()
            end if
            if(valid)cycle
            call list%delete()
            error='Invalid curved volume level '//integer_text(i)// &
                '. Supply z and radius, with control_radius or sectors only after the first level.'
            return
        end do
        call list%delete()
        call build_curved_volume_plan(levels,patches,names,joins,error)
        if(len(error)>0)return

        mesh_id=0
        if(allocated(cells))then
            call publish_multipatch_volume(s%document%multipatches,name,patches,names,joins,group,error,cells+1,mesh,catalog=s%document%patches)
            if(len(error)>0)return
            call append_studio_mesh(s,mesh)
            mesh_id=studio_mesh_count(s)
            item=mesh_identity(s,mesh_id)
            call result%insert('mesh',item)
            call item%delete()
        else
            call publish_multipatch_volume(s%document%multipatches,name,patches,names,joins,group,error,catalog=s%document%patches)
            if(len(error)>0)return
        end if
        call on_assistant_model_created(s,group,3,0,mesh_id)
        item=assembly_identity(s,group)
        call result%insert('assembly',item)
        call item%delete()
        call result%insert('patch_count',size(patches))
        call result%insert('interface_count',size(joins))
        call describe_selection(s,result)
    end subroutine build_curved_volume

    logical function parse_model_step(args,step) result(valid)
        type(QJsonObject),intent(in) :: args
        type(model_step),intent(out) :: step
        character(:),allocatable :: name,op
        real(dp) :: vector(3),angle

        valid=.false.
        op=args%stringValue('op')
        name=trim(args%stringValue('name'))
        if(len(name)==0.or.len(name)>len(step%name).or.index(name,achar(0))/=0)return
        step%name=name
        if(any(op==primitive_actions))then
            if(.not.only_fields(args,'|op|name|origin|values|'))return
            if(.not.vector_arg(args,'origin',vector))return
            step%origin=vector
            if(op/='point'.or.args%contains('values'))then
                if(.not.vector_arg(args,'values',vector))return
                step%values=vector
            end if
        else
            select case(op)
            case('extrude')
                if(.not.only_fields(args,'|op|name|source|values|'))return
                step%source=bounded_integer(args,'source',1,model_plan_max_steps)
                if(step%source<1.or..not.vector_arg(args,'values',vector))return
                step%values=vector
            case('revolve')
                if(.not.only_fields(args,'|op|name|source|origin|values|angle|'))return
                step%source=bounded_integer(args,'source',1,model_plan_max_steps)
                if(step%source<1.or..not.vector_arg(args,'origin',vector))return
                step%origin=vector
                if(.not.vector_arg(args,'values',vector))return
                step%values=vector
                if(.not.number_arg(args,'angle',angle))return
                step%angle=angle
            case('sweep')
                if(.not.only_fields(args,'|op|name|source|spine|'))return
                step%source=bounded_integer(args,'source',1,model_plan_max_steps)
                step%spine=bounded_integer(args,'spine',1,model_plan_max_steps)
                if(min(step%source,step%spine)<1)return
            case('loft')
                if(.not.only_fields(args,'|op|name|sections|degree|'))return
                step%degree=bounded_integer(args,'degree',0,model_plan_max_steps,fallback=0)
                if(step%degree<0)return
                if(.not.integer_array_arg(args,'sections',step%sections,2,model_plan_max_steps,model_plan_max_steps))return
            case default
                return
            end select
        end if
        step%op=op
        valid=.true.
    end function parse_model_step

    logical function parse_curved_level(args,level,first) result(valid)
        type(QJsonObject),intent(in) :: args
        type(curved_volume_level),intent(out) :: level
        logical,intent(in) :: first

        valid=.false.
        if(first)then
            if(.not.only_fields(args,'|z|radius|'))return
        else
            if(.not.only_fields(args,'|z|radius|control_radius|sectors|'))return
        end if
        if(.not.number_arg(args,'z',level%z))return
        if(.not.number_arg(args,'radius',level%radius))return
        if(args%contains('control_radius'))then
            if(.not.number_arg(args,'control_radius',level%control_radius))return
        end if
        if(args%contains('sectors'))then
            if(.not.integer_array_arg(args,'sectors',level%sectors,1,curved_volume_sectors,curved_volume_sectors))return
        end if
        valid=.true.
    end function parse_curved_level

end module studio_mcp_construction
