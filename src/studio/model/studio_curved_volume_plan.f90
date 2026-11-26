! Build ForCAD volume patches from axial stations and independent per-span sector selections.
module studio_curved_volume_plan
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use studio_forcad, only: rk, multipatch_join, studio_native_patch, stage_curved_volume_patches, &
        core_count=>curved_core_count, curved_volume_sectors
    implicit none(type, external)
    private

    public :: curved_volume_level, build_curved_volume_plan
    public :: curved_volume_sectors
    integer, parameter, public :: curved_volume_max_levels=12
    integer, parameter :: section_count=core_count+curved_volume_sectors
    integer, parameter :: core_neighbour(8)=[4,4,3,3,1,1,2,2]
    integer, parameter :: core_face(8)=[2,4,4,1,1,3,3,2]
    logical, parameter :: core_reversed(8)=[.false.,.true.,.true.,.true.,.true.,.false.,.false.,.false.]

    type :: curved_volume_level
        real(rk) :: z=0._rk
        real(rk) :: radius=0._rk
        ! Radius of the middle quadratic control section for the span ending here.
        ! Zero uses the linear midpoint. The first level has no preceding span.
        real(rk) :: control_radius=0._rk
        ! An omitted mask builds the full disk. Otherwise build only these annular sectors.
        integer,allocatable :: sectors(:)
    end type curved_volume_level

contains

    subroutine build_curved_volume_plan(levels,patches,names,joins,error)
        type(curved_volume_level),intent(in) :: levels(:)
        type(studio_native_patch),allocatable,intent(out) :: patches(:)
        character(96),allocatable,intent(out) :: names(:)
        type(multipatch_join),allocatable,intent(out) :: joins(:)
        character(:),allocatable,intent(out) :: error
        type(multipatch_join),allocatable :: proposed(:)
        integer,allocatable :: slot(:,:)
        real(rk),allocatable :: heights(:,:),radii(:,:)
        integer :: layers,level,sector,next,id,i,part,njoins

        error='A curved volume requires 2 to 12 axial levels.'
        if(size(levels)<2.or.size(levels)>curved_volume_max_levels)return
        error='Use finite coordinates and radii.'
        if(.not.all(ieee_is_finite([levels%z,levels%radius,levels%control_radius])))return
        error='Each axial radius must be positive and each control radius nonnegative.'
        if(minval(levels%radius)<=0._rk.or.minval(levels%control_radius)<0._rk)return
        error='Axial levels must be in strictly increasing Z order.'
        if(any(levels(2:)%z<=levels(:size(levels)-1)%z))return
        error='The first level has no preceding span, control radius, or sector selection.'
        if(levels(1)%control_radius/=0._rk.or.allocated(levels(1)%sectors))return

        layers=size(levels)-1
        allocate(slot(section_count,layers))
        slot=0
        part=0
        do level=1,layers
            if(allocated(levels(level+1)%sectors))then
                associate(selection=>levels(level+1)%sectors)
                    error='A sector selection must contain distinct indices from 1 to 8.'
                    if(size(selection)==0.or.size(selection)>curved_volume_sectors)return
                    if(any(selection<1).or.any(selection>curved_volume_sectors))return
                    if(any([(any(selection(:i-1)==selection(i)),i=1,size(selection))]))return
                    do sector=1,curved_volume_sectors
                        if(.not.any(selection==sector))cycle
                        part=part+1
                        slot(core_count+sector,level)=part
                    end do
                end associate
            else
                slot(:,level)=[(part+i,i=1,section_count)]
                part=part+section_count

            end if
        end do
        error='A curved multipatch needs at least two volume patches.'
        if(part<2)return

        allocate(proposed(3*part))
        njoins=0
        do level=1,layers
            call add_join(slot(1,level),2,slot(2,level),1)
            call add_join(slot(3,level),2,slot(4,level),1)
            call add_join(slot(1,level),4,slot(3,level),3)
            call add_join(slot(2,level),4,slot(4,level),3)
            do sector=1,curved_volume_sectors
                id=slot(core_count+sector,level)
                next=slot(core_count+mod(sector,curved_volume_sectors)+1,level)
                call add_join(id,4,next,3)
                call add_join(slot(core_neighbour(sector),level),core_face(sector),id,1,core_reversed(sector))
            end do
            if(level==1)cycle
            do i=1,section_count
                call add_join(slot(i,level-1),6,slot(i,level),5)
            end do
        end do
        error='All spans and selected sectors must connect through matching whole faces.'
        if(.not.face_connected(part,proposed(:njoins)))return

        allocate(heights(3,layers),radii(3,layers))
        do level=1,layers
            heights(:,level)=[levels(level)%z,.5_rk*(levels(level)%z+levels(level+1)%z),levels(level+1)%z]
            radii(:,level)=[levels(level)%radius,levels(level+1)%control_radius,levels(level+1)%radius]
            if(radii(2,level)==0._rk)radii(2,level)=.5_rk*(radii(1,level)+radii(3,level))
        end do
        call stage_curved_volume_patches(heights,radii,slot,proposed(:njoins),patches,error)
        if(len(error)>0)return

        allocate(names(part))
        do level=1,layers
            do i=1,core_count
                id=slot(i,level)
                if(id>0)write(names(id),'(a,i0,a,i0)')'Span ',level,' · Core ',i
            end do
            do sector=1,curved_volume_sectors
                id=slot(core_count+sector,level)
                if(id>0)write(names(id),'(a,i0,a,i0)')'Span ',level,' · Sector ',sector
            end do
        end do
        joins=proposed(:njoins)
        error=''

    contains

        subroutine add_join(a,side_a,b,side_b,reverse)
            integer,intent(in) :: a,side_a,b,side_b
            logical,optional,intent(in) :: reverse

            if(a==0.or.b==0)return
            njoins=njoins+1
            proposed(njoins)=multipatch_join(patch_a=a,side_a=side_a,patch_b=b,side_b=side_b)
            if(present(reverse))proposed(njoins)%flip(1)=reverse
        end subroutine add_join
    end subroutine build_curved_volume_plan

    logical function face_connected(part,joins) result(valid)
        integer,intent(in) :: part
        type(multipatch_join),intent(in) :: joins(:)
        logical :: visited(part)
        integer :: previous,i,a,b

        visited=.false.
        visited(1)=.true.
        do
            previous=count(visited)
            do i=1,size(joins)
                a=joins(i)%patch_a
                b=joins(i)%patch_b
                if(.not.(visited(a).or.visited(b)))cycle
                visited(a)=.true.
                visited(b)=.true.
            end do
            if(count(visited)==previous)exit
        end do
        valid=all(visited)
    end function face_connected
end module studio_curved_volume_plan
