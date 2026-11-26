module studio_selection
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none(type, external)
    private

    integer, parameter, public :: selection_patch=1, selection_assembly=2, selection_control=3, selection_mesh=4
    integer, parameter, public :: selection_replace=0, selection_add=1, selection_toggle=2
    integer, parameter, public :: pick_automatic=0, pick_face=1, pick_edge=2, pick_vertex=3, pick_control_points=4
    integer, parameter, public :: pick_multipatch=5
    public :: patch_edge_boundaries, change_selection_batch

    type, public :: selection_target_t
        integer :: category=0, id=0, rank=0
        integer :: part=0, item=0, feature=0
    end type selection_target_t

    type, public :: selection_set_t
        type(selection_target_t), allocatable :: items(:)
    contains
        procedure :: contains => contains_target
        procedure :: change => change_selection
        procedure :: clear => clear_selection
        procedure :: count => selection_count
    end type selection_set_t
contains

    pure function patch_edge_boundaries(rank,edge) result(sides)
        integer, intent(in) :: rank,edge
        integer :: sides(2),axis,flags

        sides=0
        select case(rank)
        case(2)
            if(edge>=1.and.edge<=4)sides(1)=edge
        case(3)
            if(edge<1.or.edge>12)return
            axis=(edge+3)/4
            flags=mod(edge-1,4)
            sides=[2*mod(axis,3)+1+ibits(flags,0,1),2*mod(axis+1,3)+1+ibits(flags,1,1)]
        case default
            continue
        end select
    end function patch_edge_boundaries


    integer function selection_count(selection) result(count)
        class(selection_set_t), intent(in) :: selection

        count=0
        if(allocated(selection%items))count=size(selection%items)
    end function selection_count

    logical function contains_target(selection,target) result(found)
        class(selection_set_t), intent(in) :: selection
        type(selection_target_t), intent(in) :: target

        found=target_index(selection,target)>0
    end function contains_target

    integer function target_index(selection,target) result(index)
        class(selection_set_t), intent(in) :: selection
        type(selection_target_t), intent(in) :: target
        integer :: i

        index=0
        do i=1,selection%count()
            if(.not.same_target(selection%items(i),target))cycle
            index=i
            return
        end do
    end function target_index

    pure logical function same_target(first,second) result(same)
        type(selection_target_t), intent(in) :: first,second

        same=first%category==second%category.and.first%id==second%id.and.first%rank==second%rank.and. &
            first%part==second%part.and.first%item==second%item.and. &
            (first%category==selection_mesh.or.first%feature==second%feature)
    end function same_target

    subroutine change_selection(selection,target,mode)
        class(selection_set_t), intent(inout) :: selection
        type(selection_target_t), intent(in) :: target
        integer, intent(in) :: mode
        integer :: index,count

        if(mode==selection_replace)call selection%clear()
        if(target%category==0.or.target%id<1)return
        if(.not.allocated(selection%items))call selection%clear()
        count=selection%count()
        index=target_index(selection,target)
        if(index==0)then
            selection%items=[selection%items,target]
        else if(mode==selection_toggle)then
            selection%items=[selection%items(:index-1),selection%items(index+1:)]
        else if(index<count)then
            selection%items=[selection%items(:index-1),selection%items(index+1:),target]
        end if
    end subroutine change_selection

    subroutine change_selection_batch(selection,targets,mode)
        class(selection_set_t), intent(inout) :: selection
        type(selection_target_t), intent(in) :: targets(:)
        integer, intent(in) :: mode
        type(selection_target_t), allocatable :: entries(:),ordered(:)
        integer, allocatable :: slots(:),previous(:),next(:)
        integer(int64), parameter :: weights(6)=[13_int64,109_int64,997_int64,2179_int64,4129_int64,8189_int64]
        integer :: old_count,capacity,i,n,count,head,tail,position,slot,feature,item_mode,before,after
        type(selection_target_t) :: target

        if(size(targets)==0)then
            if(mode==selection_replace)call selection%clear()
            return
        end if
        if(size(targets)==1)then
            call selection%change(targets(1),mode)
            return
        end if
        old_count=0
        if(mode/=selection_replace)old_count=selection%count()
        capacity=old_count+size(targets)
        allocate(entries(capacity),previous(capacity),next(capacity))
        allocate(slots(shiftl(1_int64,bit_size(capacity)-leadz(capacity)+1)),source=0)
        n=0; count=0; head=0; tail=0
        do i=1,capacity
            if(i<=old_count)then
                target=selection%items(i)
                item_mode=selection_add
            else
                target=targets(i-old_count)
                item_mode=mode
            end if
            if(target%category==0.or.target%id<1)cycle
            feature=merge(0,target%feature,target%category==selection_mesh)
            slot=1+int(modulo(sum(int([target%category,target%id,target%rank,target%part,target%item,feature],int64)*weights), &
                size(slots,kind=int64)))
            do
                position=slots(slot)
                if(position==0)exit
                if(same_target(entries(position),target))exit
                slot=1+modulo(slot,size(slots))
            end do
            if(position==0)then
                n=n+1
                position=n
                slots(slot)=position
            else if(head==position.or.previous(position)>0)then
                if(item_mode/=selection_toggle.and.tail==position)cycle
                before=previous(position)
                after=next(position)
                if(before==0)then
                    head=after
                else
                    next(before)=after
                end if
                if(after==0)then
                    tail=before
                else
                    previous(after)=before
                end if
                count=count-1
                previous(position)=0
                if(item_mode==selection_toggle)cycle
            end if
            entries(position)=target
            previous(position)=tail
            next(position)=0
            if(tail==0)then
                head=position
            else
                next(tail)=position
            end if
            tail=position
            count=count+1
        end do
        if(mode/=selection_replace.and..not.allocated(selection%items).and.n==0)return
        allocate(ordered(count))
        position=head
        do i=1,count
            ordered(i)=entries(position)
            position=next(position)
        end do
        call move_alloc(ordered,selection%items)
    end subroutine change_selection_batch

    subroutine clear_selection(selection)
        class(selection_set_t), intent(inout) :: selection

        if(allocated(selection%items))deallocate(selection%items)
        allocate(selection%items(0))
    end subroutine clear_selection

end module studio_selection
