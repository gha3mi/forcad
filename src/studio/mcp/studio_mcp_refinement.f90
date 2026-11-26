module studio_mcp_refinement
    use qt, only: QJsonArray, QJsonObject
    use studio_state, only: studio_t, dp
    use studio_mcp_arguments, only: only_fields, bounded_integer, number_arg, axis_number, read_geometry_identity
    use studio_multipatch_collection, only: count_multipatches
    use studio_multipatch_refinement, only: multipatch_refinement_plan_t, plan_multipatch_refinement
    use studio_single_patch, only: single_patch_name
    use studio_presentation_text, only: integer_text, number_text
    implicit none(type, external)
    private

    public :: read_mcp_refinement, preview_refinement, refinement_approval_summary

    character(*), parameter :: axes='UVW'
    integer, parameter :: summary_rows=6, default_page_size=64, max_page_size=128

contains

    subroutine read_mcp_refinement(args,op,member,axis,steps,knot,error,paginated)
        type(QJsonObject), intent(in) :: args
        character(*), intent(in) :: op
        integer, intent(out) :: member,axis,steps
        real(dp), intent(out) :: knot
        character(:), allocatable, intent(out) :: error
        logical, intent(in), optional :: paginated
        character(:), allocatable :: allowed

        member=0
        axis=0
        steps=0
        knot=0._dp
        error='Use insert_knot or elevate_degree to refine a multipatch.'
        if(op/='insert_knot'.and.op/='elevate_degree')return
        allowed='|op|target|member|axis|count|'
        if(op=='insert_knot')allowed=allowed//'value|'
        if(present(paginated))then
            if(paginated)allowed=allowed//'offset|limit|'
        end if
        error='Supply a starting member patch ID, parameter direction and positive count.'
        if(.not.only_fields(args,allowed))return
        member=bounded_integer(args,'member',1,huge(1))
        axis=axis_number(args%stringValue('axis'),axes)
        steps=bounded_integer(args,'count',1,huge(1),fallback=1)
        if(member<1.or.axis<1.or.steps<1)return
        if(op=='insert_knot')then
            error='Supply a finite interior knot value in the starting member’s interval.'
            if(.not.number_arg(args,'value',knot))return
        end if
        error=''
    end subroutine read_mcp_refinement

    subroutine preview_refinement(s,args,result,error)
        type(studio_t), intent(in) :: s
        type(QJsonObject), intent(in) :: args
        type(QJsonObject), intent(inout) :: result
        character(:), allocatable, intent(out) :: error
        type(multipatch_refinement_plan_t) :: plan
        type(QJsonArray) :: directions
        type(QJsonObject) :: target,row
        character(:), allocatable :: op
        integer :: assembly,offset,limit,first,last,n,i,axis
        logical :: is_assembly

        call refinement_target(s,args,assembly,is_assembly)
        error='Select one multipatch or supply its assembly ID.'
        if(.not.is_assembly.or.assembly<1)return
        offset=bounded_integer(args,'offset',0,huge(1),fallback=0)
        limit=bounded_integer(args,'limit',1,max_page_size,fallback=default_page_size)
        error='Supply a zero-based offset and a page size from 1 to '//integer_text(max_page_size)//'.'
        if(offset<0.or.limit<1)return

        op=args%stringValue('op')
        call prepare_refinement(s,args,assembly,op,plan,error,paginated=.true.)
        if(len(error)>0)return
        n=size(plan%rows)
        if(offset>=n)then
            error='The offset exceeds the available refinement directions.'
            return
        end if
        first=offset+1
        last=offset+min(limit,n-offset)
        target=QJsonObject()
        call target%insert('kind','assembly')
        call target%insert('id',assembly)
        call target%insert('name',s%document%multipatches%groups(assembly)%definition%name)
        call result%insert('target',target)
        call target%delete()

        directions=QJsonArray()
        do i=first,last
            row=QJsonObject()
            axis=plan%rows(i)%axis
            call row%insert('member',plan%rows(i)%member)
            call row%insert('name',member_name(s,plan%kind,plan%rows(i)%member))
            call row%insert('axis',axes(axis:axis))
            if(plan%inserting)call row%insert('knot',plan%rows(i)%knot)
            call directions%append(row)
            call row%delete()
        end do
        call result%insert('op',op)
        call result%insert('count',plan%steps)
        call result%insert('affected_patches',affected_patches(plan))
        call result%insert('total_directions',n)
        call result%insert('offset',offset)
        call result%insert('next_offset',last)
        call result%insert('has_more',last<n)
        call result%insert('directions',directions)
        call result%insert('summary',refinement_summary(s,plan,first,last))
        call directions%delete()
    end subroutine preview_refinement

    subroutine refinement_approval_summary(s,args,summary,error)
        type(studio_t), intent(in) :: s
        type(QJsonObject), intent(in) :: args
        character(:), allocatable, intent(out) :: summary,error
        type(multipatch_refinement_plan_t) :: plan
        character(:), allocatable :: op
        integer :: assembly
        logical :: is_assembly

        summary=''
        error=''
        op=args%stringValue('op')
        if(op/='insert_knot'.and.op/='elevate_degree')return
        call refinement_target(s,args,assembly,is_assembly)
        if(.not.is_assembly)return

        call prepare_refinement(s,args,assembly,op,plan,error,paginated=.false.)
        if(len(error)==0)summary=refinement_summary(s,plan,1,size(plan%rows))
    end subroutine refinement_approval_summary

    subroutine refinement_target(s,args,assembly,is_assembly)
        type(studio_t), intent(in) :: s
        type(QJsonObject), intent(in) :: args
        integer, intent(out) :: assembly
        logical, intent(out) :: is_assembly
        character(:), allocatable :: kind
        logical :: valid

        assembly=0
        is_assembly=.false.
        if(args%contains('target'))then
            valid=read_geometry_identity(args,kind,assembly)
            if(.not.allocated(kind))return
            is_assembly=kind=='assembly'.or.kind=='multipatch'
            if(.not.valid)assembly=0
        else
            is_assembly=s%selected_multipatch>0
            if(s%selection%count()<=1)assembly=s%selected_multipatch
        end if
    end subroutine refinement_target

    subroutine prepare_refinement(s,args,assembly,op,plan,error,paginated)
        type(studio_t), intent(in) :: s
        type(QJsonObject), intent(in) :: args
        integer, intent(in) :: assembly
        character(*), intent(in) :: op
        type(multipatch_refinement_plan_t), intent(out) :: plan
        character(:), allocatable, intent(out) :: error
        logical, intent(in) :: paginated
        integer :: member,axis,steps
        real(dp) :: knot

        error='Select an available multipatch.'
        if(assembly<1.or.assembly>count_multipatches(s%document%multipatches))return
        call read_mcp_refinement(args,op,member,axis,steps,knot,error,paginated)
        if(len(error)>0)return
        if(op=='insert_knot')then
            call plan_multipatch_refinement(s%document,assembly,member,axis,steps,plan,error,knot=knot)
        else
            call plan_multipatch_refinement(s%document,assembly,member,axis,steps,plan,error)
        end if
    end subroutine prepare_refinement

    integer function affected_patches(plan) result(patches)
        type(multipatch_refinement_plan_t), intent(in) :: plan
        integer :: i,previous

        patches=0
        previous=0
        do i=1,size(plan%rows)
            if(plan%rows(i)%member/=previous)patches=patches+1
            previous=plan%rows(i)%member
        end do
    end function affected_patches

    function refinement_summary(s,plan,first,last) result(summary)
        type(studio_t), intent(in) :: s
        type(multipatch_refinement_plan_t), intent(in) :: plan
        integer, intent(in) :: first,last
        character(:), allocatable :: summary,operation,name,line
        character, parameter :: newline=new_line('a')
        integer :: i,shown,n

        operation='Elevate degree by '//integer_text(plan%steps)
        if(plan%inserting)operation='Insert knot · Count '//integer_text(plan%steps)
        name=summary_name(s%document%multipatches%groups(plan%assembly)%definition%name)
        n=size(plan%rows)
        summary=operation//' in '//name//newline//integer_text(affected_patches(plan))//' affected patches · '// &
            integer_text(n)//' parameter directions'
        shown=min(last,first+summary_rows-1)
        do i=first,shown
            associate(row=>plan%rows(i))
                line=summary_name(member_name(s,plan%kind,row%member))//' (#'//integer_text(row%member)//') · '//axes(row%axis:row%axis)
                if(plan%inserting)line=line//' · Knot '//number_text(row%knot)
                summary=summary//newline//line
            end associate
        end do
        if(first>1.or.shown<n)then
            if(shown>=first)summary=summary//newline//'Showing directions '//integer_text(first)//'–'//integer_text(shown)//' of '//integer_text(n)
            summary=summary//newline//'The preview_refinement tool can return all direction rows in pages.'
        end if
    end function refinement_summary

    function member_name(s,kind,member) result(name)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: kind,member
        character(:), allocatable :: name

        name=single_patch_name(kind,member,catalog=s%document%patches)
    end function member_name

    function summary_name(name) result(display)
        character(*), intent(in) :: name
        character(:), allocatable :: display
        integer, parameter :: name_bytes=56
        integer :: last

        display=name
        if(len(display)<=name_bytes)return
        last=name_bytes
        do while(iand(iachar(display(last+1:last+1)),192)==128)
            last=last-1
        end do
        display=display(:last)//'…'
    end function summary_name

end module studio_mcp_refinement
