! Prepare face-connected NURBS volume patches from a bounded construction plan.
module studio_volume_plan
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use studio_forcad, only: rk, multipatch_join, studio_native_patch, stage_box_volume_patches
    implicit none(type, external)
    private
    public :: volume_block, build_volume_plan, volume_plan_max_blocks, volume_plan_max_patches

    integer, parameter :: volume_plan_max_blocks=20, volume_plan_max_patches=160

    type :: volume_block
        character(64) :: name=''
        real(rk) :: lower(3)=0._rk, upper(3)=0._rk
    end type volume_block

    type :: indexed_box
        integer :: lower(3), upper(3)
    end type indexed_box

contains

    subroutine build_volume_plan(blocks,patches,names,joins,error)
        type(volume_block), intent(in) :: blocks(:)
        type(studio_native_patch), allocatable, intent(out) :: patches(:)
        character(96), allocatable, intent(out) :: names(:)
        ! Joins refer to one-based local patch IDs until the caller assigns catalog IDs.
        type(multipatch_join), allocatable, intent(out) :: joins(:)
        character(:), allocatable, intent(out) :: error
        type(indexed_box), allocatable :: bounds(:)
        type(indexed_box), allocatable :: cells(:)
        type(multipatch_join), allocatable :: proposed(:)
        real(rk), allocatable :: cuts(:,:), lower(:,:), upper(:,:)
        integer, allocatable :: contacts(:,:), source_count(:)
        integer :: ncuts(3)
        integer :: nblocks, part, i, j, axis, shared_axis, shared_span(3), njoins
        logical, allocatable :: block_cuts(:,:,:)
        logical :: tangential(3)
        character(96), allocatable :: part_names(:)
        character(16) :: number

        error='A multipatch plan needs 2 to 20 volume blocks.'
        nblocks=size(blocks)
        if(nblocks<2.or.nblocks>volume_plan_max_blocks)return

        do i=1,nblocks
            error='Every volume block needs a distinct name.'
            if(len_trim(blocks(i)%name)==0.or.any(blocks(:i-1)%name==blocks(i)%name))return
            error='Use finite coordinates for '//trim(blocks(i)%name)//'.'
            if(.not.all(ieee_is_finite(blocks(i)%lower)).or..not.all(ieee_is_finite(blocks(i)%upper)))return
            error='Upper bounds must exceed lower bounds on all axes for '//trim(blocks(i)%name)//'.'
            if(any(blocks(i)%upper<=blocks(i)%lower))return
        end do

        allocate(bounds(nblocks),cuts(2*nblocks,3),contacts(nblocks,nblocks),source_count(nblocks))
        allocate(block_cuts(2*nblocks,3,nblocks))
        block_cuts=.false.
        do axis=1,3
            call unique_axis_cuts([blocks%lower(axis),blocks%upper(axis)],cuts(:,axis),ncuts(axis))
            do i=1,nblocks
                bounds(i)%lower(axis)=findloc(cuts(:ncuts(axis),axis),blocks(i)%lower(axis),dim=1)
                bounds(i)%upper(axis)=findloc(cuts(:ncuts(axis),axis),blocks(i)%upper(axis),dim=1)
                block_cuts(bounds(i)%lower(axis),axis,i)=.true.
                block_cuts(bounds(i)%upper(axis),axis,i)=.true.
            end do
        end do

        contacts=0
        do j=2,nblocks
            do i=1,j-1
                shared_span=min(bounds(i)%upper,bounds(j)%upper)-max(bounds(i)%lower,bounds(j)%lower)
                if(all(shared_span>0))then
                    error='Volume blocks overlap: '//trim(blocks(i)%name)//' and '//trim(blocks(j)%name)//'.'
                    return
                end if
                contacts(i,j)=contact_axis(shared_span)
                contacts(j,i)=contacts(i,j)
            end do
        end do
        error='Every volume block must connect to the others through faces with shared area.'
        if(.not.connected_blocks(contacts))return

        ! A touching face may cover only part of a block. Propagate its tangential
        ! cuts to both blocks until each interface is made of matching whole faces.
        call align_face_cuts(bounds,contacts,block_cuts)
        do i=1,nblocks
            source_count(i)=product(count(block_cuts(:,:,i),dim=1)-1)
        end do
        error='The plan needs more than 160 patches after matching interfaces. Simplify the block layout.'
        if(sum(source_count)>volume_plan_max_patches)return

        allocate(cells(sum(source_count)))
        allocate(part_names(size(cells)),proposed(3*size(cells)))
        part=0
        do i=1,nblocks
            call split_block(block_cuts(:,:,i),cells(part+1:part+source_count(i)))
            do j=1,source_count(i)
                part_names(part+j)=blocks(i)%name
                if(source_count(i)==1)cycle
                write(number,'(i0)')j
                part_names(part+j)=trim(blocks(i)%name)//' · '//trim(number)
            end do
            part=part+source_count(i)
        end do

        njoins=0
        do j=2,size(cells)
            do i=1,j-1
                shared_span=min(cells(i)%upper,cells(j)%upper)-max(cells(i)%lower,cells(j)%lower)
                shared_axis=contact_axis(shared_span)
                if(shared_axis==0)cycle
                tangential=[1,2,3]/=shared_axis
                if(any(tangential.and.cells(i)%lower/=cells(j)%lower).or. &
                    any(tangential.and.cells(i)%upper/=cells(j)%upper))then
                    error='The volume plan could not produce matching faces.'
                    return
                end if
                njoins=njoins+1
                associate(join=>proposed(njoins))
                    join%patch_a=i
                    join%patch_b=j
                    join%side_a=2*shared_axis-1
                    join%side_b=2*shared_axis
                    if(cells(i)%upper(shared_axis)==cells(j)%lower(shared_axis))then
                        join%side_a=2*shared_axis
                        join%side_b=2*shared_axis-1
                    end if
                end associate
            end do
        end do
        allocate(lower(3,size(cells)),upper(3,size(cells)))
        do axis=1,3
            lower(axis,:)=cuts(cells%lower(axis),axis)
            upper(axis,:)=cuts(cells%upper(axis),axis)
        end do
        call stage_box_volume_patches(lower,upper,proposed(:njoins),patches,error)
        if(len(error)>0)return
        joins=proposed(:njoins)
        call move_alloc(part_names,names)
    end subroutine build_volume_plan

    subroutine unique_axis_cuts(input,cuts,count)
        real(rk), intent(in) :: input(:)
        real(rk), intent(out) :: cuts(:)
        integer, intent(out) :: count
        real(rk) :: sorted(size(input)), value
        integer :: i,j

        sorted=input
        do i=2,size(sorted)
            value=sorted(i); j=i-1
            do while(j>=1)
                if(sorted(j)<=value)exit
                sorted(j+1)=sorted(j)
                j=j-1
            end do
            sorted(j+1)=value
        end do
        count=1; cuts(1)=sorted(1)
        do i=2,size(sorted)
            if(sorted(i)==cuts(count))cycle
            count=count+1
            cuts(count)=sorted(i)
        end do
    end subroutine unique_axis_cuts

    integer function contact_axis(span) result(axis)
        integer, intent(in) :: span(3)

        axis=0
        if(count(span>0)==2.and.count(span==0)==1)axis=findloc(span,0,dim=1)
    end function contact_axis

    logical function connected_blocks(contacts) result(connected)
        integer, intent(in) :: contacts(:,:)
        logical :: visited(size(contacts,1))
        integer :: queue(size(contacts,1)), head, tail, j

        visited=.false.; visited(1)=.true.
        queue(1)=1; tail=1; head=1
        do while(head<=tail)
            do j=1,size(visited)
                if(visited(j).or.contacts(queue(head),j)==0)cycle
                visited(j)=.true.
                tail=tail+1
                queue(tail)=j
            end do
            head=head+1
        end do
        connected=all(visited)
    end function connected_blocks

    subroutine align_face_cuts(bounds,contacts,block_cuts)
        type(indexed_box), intent(in) :: bounds(:)
        integer, intent(in) :: contacts(:,:)
        logical, intent(inout) :: block_cuts(:,:,:)
        logical :: unchanged, keep
        integer :: a,b,axis,lo,hi,k,previous

        do
            previous=count(block_cuts)
            do b=2,size(bounds)
                do a=1,b-1
                    if(contacts(a,b)==0)cycle
                    do axis=1,3
                        if(axis==contacts(a,b))cycle
                        lo=max(bounds(a)%lower(axis),bounds(b)%lower(axis))
                        hi=min(bounds(a)%upper(axis),bounds(b)%upper(axis))
                        do k=lo,hi
                            keep=block_cuts(k,axis,a).or.block_cuts(k,axis,b)
                            block_cuts(k,axis,a)=keep
                            block_cuts(k,axis,b)=keep
                        end do
                    end do
                end do
            end do
            unchanged=count(block_cuts)==previous
            if(unchanged)exit
        end do
    end subroutine align_face_cuts

    subroutine split_block(block_cuts,cells)
        logical, intent(in) :: block_cuts(:,:)
        type(indexed_box), intent(out) :: cells(:)
        integer :: indices(size(block_cuts,1),3), sizes(3), axis, i,j,k, part

        do axis=1,3
            sizes(axis)=count(block_cuts(:,axis))
            indices(:sizes(axis),axis)=pack([(i,i=1,size(block_cuts,1))],block_cuts(:,axis))
        end do
        part=0
        do k=1,sizes(3)-1
            do j=1,sizes(2)-1
                do i=1,sizes(1)-1
                    part=part+1
                    cells(part)%lower=[indices(i,1),indices(j,2),indices(k,3)]
                    cells(part)%upper=[indices(i+1,1),indices(j+1,2),indices(k+1,3)]
                end do
            end do
        end do
    end subroutine split_block
end module studio_volume_plan
