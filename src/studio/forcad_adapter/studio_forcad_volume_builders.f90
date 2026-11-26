module studio_forcad_volume_builders
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use forcad, only: rk, nurbs_volume, nurbs_multipatch_volume
    use studio_forcad_assembly, only: multipatch_join, multipatch_side_label, append_assembly_patch
    use studio_forcad_transforms, only: studio_native_patch
    implicit none(type, external)
    private

    public :: stage_box_volume_patches, stage_curved_volume_patches
    integer, parameter, public :: curved_core_count=4, curved_volume_sectors=8

    real(rk), parameter :: quadratic(6)=[0._rk,0._rk,0._rk,1._rk,1._rk,1._rk]
    real(rk), parameter :: linear(4)=[0._rk,0._rk,1._rk,1._rk]
    real(rk), parameter :: half_diagonal=sqrt(.5_rk)
    real(rk), parameter :: square_scale=half_diagonal*.5_rk
    real(rk), parameter :: arc_weight=sqrt(.5_rk*(1._rk+half_diagonal))
    ! The 45 degree sectors start at +X and progress counterclockwise.
    real(rk), parameter :: direction(2,8)=reshape([1._rk,0._rk, half_diagonal,half_diagonal, &
        0._rk,1._rk, -half_diagonal,half_diagonal, -1._rk,0._rk, &
        -half_diagonal,-half_diagonal, 0._rk,-1._rk, half_diagonal,-half_diagonal],[2,8])
    real(rk), parameter :: square(2,8)=square_scale*reshape([1._rk,0._rk, 1._rk,1._rk, &
        0._rk,1._rk, -1._rk,1._rk, -1._rk,0._rk, -1._rk,-1._rk, 0._rk,-1._rk, 1._rk,-1._rk],[2,8])

contains

    subroutine stage_box_volume_patches(lower,upper,joins,patches,error)
        real(rk), intent(in) :: lower(:,:),upper(:,:)
        type(multipatch_join), intent(in) :: joins(:)
        type(studio_native_patch), allocatable, intent(out) :: patches(:)
        character(:), allocatable, intent(out) :: error
        type(nurbs_volume) :: volumes(size(lower,2))
        type(nurbs_multipatch_volume) :: native
        integer :: i

        error='ForCAD rejected a volume cell. Check the dimensions and use larger intervals.'
        if(size(lower,1)/=3.or.any(shape(lower)/=shape(upper)))return
        do i=1,size(volumes)
            call volumes(i)%set_hexahedron(upper(:,i)-lower(:,i),[2,2,2])
            call volumes(i)%translate_Xc(lower(:,i))
            if(.not.positive_volume(volumes(i)))return
            if(.not.add_volume(native,volumes(i),i))return
        end do

        error='ForCAD could not validate the connected volume interfaces.'
        if(.not.connected_volume(native,joins))return

        allocate(patches(size(volumes)))
        call copy_volume_geometry(volumes,patches)
        error=''
    end subroutine stage_box_volume_patches

    subroutine stage_curved_volume_patches(heights,radii,slots,joins,patches,error)
        integer, intent(in) :: slots(:,:)
        real(rk), intent(in) :: heights(3,size(slots,2)),radii(3,size(slots,2))
        type(multipatch_join), intent(in) :: joins(:)
        type(studio_native_patch), allocatable, intent(out) :: patches(:)
        character(:), allocatable, intent(out) :: error
        type(nurbs_volume), allocatable :: volumes(:)
        type(nurbs_multipatch_volume) :: native
        integer :: level,sector,id,i

        error='ForCAD rejected a curved volume patch. Increase the spacing between axial levels.'
        if(size(slots,1)/=curved_core_count+curved_volume_sectors.or.size(slots,2)==0)return
        allocate(volumes(maxval(slots)))
        do level=1,size(slots,2)
            do i=1,curved_core_count
                id=slots(i,level)
                if(id>0)call make_core(i,heights(:,level),radii(:,level),volumes(id))
            end do
            do sector=1,curved_volume_sectors
                id=slots(curved_core_count+sector,level)
                if(id>0)call make_sector(sector,heights(:,level),radii(:,level),volumes(id))
            end do
        end do

        do i=1,size(volumes)
            if(.not.positive_volume(volumes(i)))return
            if(.not.add_volume(native,volumes(i),i))return
        end do

        error='ForCAD could not verify the curved C0 interfaces.'
        if(.not.connected_volume(native,joins))return

        allocate(patches(size(volumes)))
        call copy_volume_geometry(volumes,patches)
        error=''
    end subroutine stage_curved_volume_patches

    logical function positive_volume(volume) result(valid)
        type(nurbs_volume), intent(inout) :: volume
        real(rk) :: measure

        valid=.false.
        if(.not.volume%is_valid())return
        call volume%cmp_volume(measure,strict=.true.)
        valid=volume%err%ok.and.ieee_is_finite(measure).and.measure>0._rk
    end function positive_volume

    logical function add_volume(native,volume,expected_id) result(added)
        type(nurbs_multipatch_volume), intent(inout) :: native
        type(nurbs_volume), intent(in) :: volume
        integer, intent(in) :: expected_id
        integer :: id

        added=.false.
        if(.not.volume%err%ok)return
        call append_assembly_patch(native,volume,id)
        added=native%err%ok.and.id==expected_id
    end function add_volume

    logical function connected_volume(native,joins) result(valid)
        type(nurbs_multipatch_volume), intent(inout) :: native
        type(multipatch_join), intent(in) :: joins(:)
        integer :: i

        valid=.false.
        do i=1,size(joins)
            associate(join=>joins(i))
                call native%connect(join%patch_a,multipatch_side_label(3,join%side_a), &
                    join%patch_b,multipatch_side_label(3,join%side_b),join%continuity,flip=join%flip)
            end associate
            if(.not.native%err%ok.or.native%get_nconnection()/=i)return
        end do
        valid=native%is_valid()
    end function connected_volume

    subroutine copy_volume_geometry(volumes,patches)
        type(nurbs_volume), intent(in) :: volumes(:)
        type(studio_native_patch), intent(out) :: patches(size(volumes))
        integer :: i

        do i=1,size(volumes)
            allocate(patches(i)%geometry,source=volumes(i))
        end do
    end subroutine copy_volume_geometry

    subroutine make_core(core,z,radius,volume)
        integer,intent(in) :: core
        real(rk),intent(in) :: z(3),radius(3)
        type(nurbs_volume),intent(out) :: volume
        real(rk) :: x(3),y(3),points(27,3)
        integer :: u,v,w,row

        x=half_square(mod(core-1,2))
        y=half_square((core-1)/2)
        row=0
        do w=1,3
            do v=1,3
                do u=1,3
                    row=row+1
                    points(row,:)=[x(u)*radius(w),y(v)*radius(w),z(w)]
                end do
            end do
        end do
        call volume%set(knot1=quadratic,knot2=quadratic,knot3=quadratic,Xc=points, &
            Wc=spread(1._rk,1,size(points,1)),degree=[2,2,2])
    end subroutine make_core

    pure function half_square(side) result(knots)
        integer,intent(in) :: side
        real(rk) :: knots(3)

        knots=square_scale*(real(side,rk)+[-1._rk,-.5_rk,0._rk])
    end function half_square

    subroutine make_sector(sector,z,radius,volume)
        integer,intent(in) :: sector
        real(rk),intent(in) :: z(3),radius(3)
        type(nurbs_volume),intent(out) :: volume
        real(rk) :: inside(2,3),outside(2,3),points(18,3),weights(18)
        integer :: next,v,w,row

        ! The outside uses the exact quadratic rational arc on a circle.
        ! Unit weights on the inner edge match the planar square core.
        next=mod(sector,curved_volume_sectors)+1
        inside(:,1)=square(:,sector)
        inside(:,3)=square(:,next)
        inside(:,2)=.5_rk*(inside(:,1)+inside(:,3))
        outside(:,1)=direction(:,sector)
        outside(:,3)=direction(:,next)
        outside(:,2)=(outside(:,1)+outside(:,3))/(2._rk*arc_weight**2)
        weights=1._rk
        do w=1,3
            do v=1,3
                row=6*(w-1)+2*(v-1)+1
                points(row,:)=[radius(w)*inside(:,v),z(w)]
                points(row+1,:)=[radius(w)*outside(:,v),z(w)]
                if(v==2)weights(row+1)=arc_weight
            end do
        end do
        call volume%set(knot1=linear,knot2=quadratic,knot3=quadratic,Xc=points,Wc=weights,degree=[1,2,2])
    end subroutine make_sector

end module studio_forcad_volume_builders
