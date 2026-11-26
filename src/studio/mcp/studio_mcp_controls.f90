module studio_mcp_controls
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use qt, only: QJsonArray,QJsonObject,QJsonValue_Array
    use studio_state, only: studio_t,dp
    use studio_forcad, only: rk,replace_native_controls
    use studio_scene, only: project_patch_exists
    use studio_single_patch, only: single_patch_controls,active_patch_id,read_single_patch_controls, &
        commit_single_patch_controls,copy_single_patch_geometry
    use studio_multipatch_controls, only: multipatch_controls_t,multipatch_control_preview_t,joined_patch_assembly, &
        assembly_control_id,read_document_multipatch_controls,preview_document_multipatch_controls, &
        commit_document_multipatch_controls,scale_control_weight
    use studio_mcp_inspection, only: active_rank,patch_kind,patch_identity,assembly_identity
    use studio_mcp_arguments, only: bounded_integer,number_arg,only_fields,read_geometry_identity,vector_arg
    use studio_presentation_text, only: integer_text,number_text
    implicit none(type,external)
    private
    public :: edit_mcp_controls,preview_mcp_controls,control_approval_summary

    integer,parameter :: summary_limit=6, page_limit=128

    type :: control_request
        integer :: kind=0,patch=0,assembly=0
        type(single_patch_controls) :: local
        type(multipatch_controls_t) :: shared
        integer,allocatable :: ids(:)
        real(dp),allocatable :: old_points(:,:),old_weights(:)
    end type control_request
contains
    subroutine edit_mcp_controls(s,args,result,error,kind,active,changed)
        type(studio_t),intent(inout),target :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        integer,intent(out) :: kind,active
        logical,intent(out) :: changed
        type(control_request) :: request
        type(QJsonArray) :: ids
        integer :: i

        call read_request(s,args,request,error)
        if(len(error)>0)return
        kind=request%kind
        active=request%patch
        if(request%assembly>0)then
            call commit_document_multipatch_controls(s%document,request%assembly,request%shared,changed,error)
            active=active_patch_id(request%kind,catalog=s%document%patches)
            if(.not.any(request%shared%members==active))active=request%shared%members(1)
        else
            call commit_single_patch_controls(request%kind,request%local,changed,error, &
                patch_id=request%patch,catalog=s%document%patches)
        end if
        if(len(error)>0)return
        call describe_target(s,request,result)
        ids=QJsonArray()
        do i=1,size(request%ids)
            call ids%append(request%ids(i))
        end do
        call result%insert('control_ids',ids)
        if(size(request%ids)==1)call result%insert('control_id',request%ids(1))
        call ids%delete()
        call result%insert('changed',changed)
    end subroutine edit_mcp_controls

    subroutine preview_mcp_controls(s,args,result,error)
        type(studio_t),intent(in) :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        type(control_request) :: request
        type(multipatch_control_preview_t) :: preview
        type(QJsonArray) :: rows
        type(QJsonObject) :: item
        integer :: offset,limit,last,n,i,id

        error='Choose a zero-based offset and a page size from 1 to '//integer_text(page_limit)//'.'
        offset=bounded_integer(args,'offset',0,huge(1),fallback=0)
        limit=bounded_integer(args,'limit',1,page_limit,fallback=64)
        if(offset<0.or.limit<1)return
        call read_request(s,args,request,error,paginated=.true.)
        if(len(error)>0)return
        call plan_request(s,request,preview,error)
        if(len(error)>0)return
        n=size(preview%ids)
        if(offset>n)then
            error='The offset exceeds the available dependent controls.'
            return
        end if
        last=offset+min(limit,n-offset)
        call describe_target(s,request,result)
        rows=QJsonArray()
        do i=offset+1,last
            id=preview%ids(i)
            item=QJsonObject()
            call item%insert('control_id',id)
            call item%insert('position_changed',preview%moved(i))
            call item%insert('weight_changed',preview%reweighted(i))
            call insert_position(item,'before_position',request%old_points(:,id))
            call insert_position(item,'after_position',preview%points(:,i))
            call item%insert('before_weight',request%old_weights(id))
            call item%insert('after_weight',preview%weights(i))
            call rows%append(item)
            call item%delete()
        end do
        call result%insert('dependents',rows)
        call rows%delete()
        call result%insert('requested_controls',size(request%ids))
        call result%insert('dependent_controls',n)
        call result%insert('offset',offset)
        call result%insert('next_offset',last)
        call result%insert('has_more',last<n)
        call result%insert('summary',summarize(request,preview,offset+1,last))
    end subroutine preview_mcp_controls

    subroutine control_approval_summary(s,args,summary,error)
        type(studio_t),intent(in) :: s
        type(QJsonObject),intent(in) :: args
        character(:),allocatable,intent(out) :: summary,error
        type(control_request) :: request
        type(multipatch_control_preview_t) :: preview

        summary=''
        call read_request(s,args,request,error)
        if(len(error)>0)return
        call plan_request(s,request,preview,error)
        if(len(error)==0)summary=summarize(request,preview,1,size(preview%ids))
    end subroutine control_approval_summary

    subroutine read_request(s,args,request,error,paginated)
        type(studio_t),intent(in) :: s
        type(QJsonObject),intent(in) :: args
        type(control_request),intent(out) :: request
        character(:),allocatable,intent(out) :: error
        logical,optional,intent(in) :: paginated
        type(QJsonArray) :: batch
        type(QJsonObject) :: row
        character(:),allocatable :: fields,name
        integer :: id,n,i
        logical,allocatable :: handled(:)

        fields='|op|target|point|position|weight|controls|'
        if(present(paginated))then
            if(paginated)fields=fields//'offset|limit|'
        end if
        error='Supply a control point and position or weight, or a controls array.'
        if(.not.only_fields(args,fields))return
        if(args%contains('op'))then
            if(args%stringValue('op')/='control')return
        end if
        error='Select one patch or multipatch, or supply its kind and ID.'
        if(args%contains('target'))then
            if(.not.read_geometry_identity(args,name,id))return
            if(name=='assembly'.or.name=='multipatch')then
                request%assembly=id
            else
                request%kind=patch_kind(name)
                request%patch=id
            end if
        else if(s%selected_multipatch>0)then
            if(s%selection%count()>1)return
            request%assembly=s%selected_multipatch
        else
            request%kind=active_rank(s)
            if(request%kind<1)return
            request%patch=active_patch_id(request%kind,catalog=s%document%patches)
        end if
        if(request%patch>0)then
            if(request%kind<1)return
            if(.not.project_patch_exists(s%document%show_templates,request%kind,request%patch,patches=s%document%patches))return
            call read_single_patch_controls(request%kind,request%local,variant=request%patch,catalog=s%document%patches)
            error='ForCAD could not read the control points.'
            if(.not.request%local%valid)return
            request%assembly=joined_patch_assembly(s%document,request%kind,request%patch)
        end if
        if(request%assembly>0)then
            call read_document_multipatch_controls(s%document,request%assembly,request%shared)
            error=request%shared%message
            if(.not.request%shared%valid)return
            request%kind=request%shared%kind
            request%old_points=request%shared%points
            request%old_weights=request%shared%weights
        else
            if(.not.request%local%valid)return
            request%old_points=request%local%points
            request%old_weights=request%local%weights
        end if
        n=size(request%old_weights)
        if(request%patch>0)n=size(request%local%weights)
        allocate(handled(size(request%old_weights)),source=.false.)
        if(args%contains('controls'))then
            error='Provide one or more unique controls in the controls array.'
            if(args%contains('point').or.args%contains('position').or.args%contains('weight'))return
            if(args%valueType('controls')/=QJsonValue_Array)return
            batch=args%arrayValue('controls')
            if(batch%size()>0.and.batch%size()<=n)then
                allocate(request%ids(batch%size()))
                do i=1,batch%size()
                    row=batch%objectAt(i-1)
                    error='Each controls entry accepts point, position, and weight.'
                    if(only_fields(row,'|point|position|weight|'))call read_update(row,request,i,n,handled,error)
                    call row%delete()
                    if(len(error)>0)exit
                end do
            end if
            call batch%delete()
        else
            allocate(request%ids(1))
            call read_update(args,request,1,n,handled,error)
        end if
    end subroutine read_request

    subroutine read_update(args,request,row,n,handled,error)
        type(QJsonObject),intent(in) :: args
        type(control_request),intent(inout) :: request
        integer,intent(in) :: row,n
        logical,intent(inout) :: handled(:)
        character(:),allocatable,intent(out) :: error
        integer :: point,id
        real(dp) :: original,weight,position(3)

        error='Each edit needs a unique point ID and a finite position or positive weight.'
        point=bounded_integer(args,'point',1,n)
        if(point<1)return
        if(.not.args%contains('position').and..not.args%contains('weight'))return
        id=point
        if(request%assembly>0.and.request%patch>0)id=assembly_control_id(request%shared,request%patch,point)
        if(id<1)return
        if(handled(id))return
        position=request%old_points(:,id)
        original=request%old_weights(id)
        if(request%patch>0)original=request%local%weights(point)
        weight=request%old_weights(id)
        if(args%contains('position'))then
            if(.not.vector_arg(args,'position',position))return
        end if
        if(args%contains('weight'))then
            if(.not.number_arg(args,'weight',weight))return
            if(weight<=0._dp)return
            if(weight==original)then
                weight=request%old_weights(id)
            else
                weight=real(scale_control_weight(real(request%old_weights(id),rk),real(original,rk),real(weight,rk)),dp)
            end if
            error='That linked weight is outside the representable range.'
            if(.not.ieee_is_finite(weight).or.weight<=0._dp)return
        end if
        if(request%assembly>0)then
            request%shared%points(:,id)=position
            request%shared%weights(id)=weight
        else
            request%local%points(:,id)=position
            request%local%weights(id)=weight
        end if
        handled(id)=.true.
        request%ids(row)=id
        error=''
    end subroutine read_update

    subroutine plan_request(s,request,preview,error)
        type(studio_t),intent(in) :: s
        type(control_request),intent(in) :: request
        type(multipatch_control_preview_t),intent(out) :: preview
        character(:),allocatable,intent(out) :: error
        class(*),allocatable :: geometry
        logical :: moved,reweighted

        if(request%assembly>0)then
            call preview_document_multipatch_controls(s%document,request%assembly,request%shared,preview,error)
            return
        end if
        moved=any(request%local%points/=request%old_points)
        reweighted=any(request%local%weights/=request%old_weights)
        if(moved.or.reweighted)then
            error='ForCAD rejected the requested control points or weights.'
            call copy_single_patch_geometry(request%kind,request%patch,geometry,catalog=s%document%patches)
            if(.not.allocated(geometry))return
            if(.not.replace_native_controls(geometry,transpose(real(request%local%points,rk)), &
                real(request%local%weights,rk),moved,reweighted))return
        end if
        allocate(preview%ids(0))
        error=''
    end subroutine plan_request

    subroutine describe_target(s,request,result)
        type(studio_t),intent(in) :: s
        type(control_request),intent(in) :: request
        type(QJsonObject),intent(inout) :: result
        type(QJsonObject) :: item

        if(request%assembly>0)then
            item=assembly_identity(s,request%assembly)
            if(request%patch>0)then
                call result%insert('assembly',item)
                call item%delete()
            end if
        end if
        if(request%patch>0)item=patch_identity(s,request%kind,request%patch)
        call result%insert('target',item)
        call item%delete()
    end subroutine describe_target

    function summarize(request,preview,first,last) result(summary)
        type(control_request),intent(in) :: request
        type(multipatch_control_preview_t),intent(in) :: preview
        integer,intent(in) :: first,last
        character(:),allocatable :: summary
        character,parameter :: line=new_line('a')
        real(dp) :: position(3),weight
        integer :: i,id

        summary='Requested controls · '//integer_text(size(request%ids))//line
        do i=1,min(summary_limit,size(request%ids))
            id=request%ids(i)
            if(request%assembly>0)then
                position=request%shared%points(:,id)
                weight=request%shared%weights(id)
            else
                position=request%local%points(:,id)
                weight=request%local%weights(id)
            end if
            summary=summary//'#'//integer_text(id)//' · XYZ '//xyz_text(position)//' · Weight '//number_text(weight)//line
        end do
        if(size(request%ids)>summary_limit)summary=summary//'…'//line
        summary=summary//'Dependent controls · '//integer_text(size(preview%ids))
        do i=first,min(last,first+summary_limit-1)
            id=preview%ids(i)
            summary=summary//line//'#'//integer_text(id)
            if(preview%moved(i))summary=summary//' · XYZ '//xyz_text(request%old_points(:,id))//' → '//xyz_text(preview%points(:,i))
            if(preview%reweighted(i))summary=summary//' · Weight '//number_text(request%old_weights(id))//' → '//number_text(preview%weights(i))
        end do
        if(size(preview%ids)>summary_limit)summary=summary//line//'Use preview_controls to read all dependent values.'
    end function summarize

    function xyz_text(position) result(text)
        real(dp),intent(in) :: position(3)
        character(:),allocatable :: text

        text='('//number_text(position(1))//', '//number_text(position(2))//', '//number_text(position(3))//')'
    end function xyz_text

    subroutine insert_position(item,key,position)
        type(QJsonObject),intent(inout) :: item
        character(*),intent(in) :: key
        real(dp),intent(in) :: position(3)
        type(QJsonArray) :: array

        array=QJsonArray()
        call array%appendDoubles(position)
        call item%insert(key,array)
        call array%delete()
    end subroutine insert_position
end module studio_mcp_controls
