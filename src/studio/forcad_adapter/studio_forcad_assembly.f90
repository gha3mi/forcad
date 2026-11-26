module studio_forcad_assembly
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fordebug, only: debug
    use studio_forcad_constraint_edit, only: complete_constraint_edit
    use forcad, only: rk, nurbs_curve, nurbs_surface, nurbs_volume, &
        nurbs_multipatch_curve, nurbs_multipatch_surface, nurbs_multipatch_volume
    use forcad_utils, only: boundary_axes
    implicit none(type, external)
    private

    public :: multipatch_join, multipatch_side_label, multipatch_join_directions
    public :: initialize_assembly, append_assembly_patch, append_assembly_join, inspect_native_assembly
    public :: native_assembly_control_map, native_assembly_control_correction

    type :: multipatch_join
        ! Patch identifiers are catalog variants, resolved to native ids when assembling.
        integer :: patch_a=0, side_a=0, patch_b=0, side_b=0, continuity=0
        logical :: reverse=.false., swap=.false., flip(2)=.false., geometric=.false.
        real(c_double), allocatable :: reparameterization(:)
    end type multipatch_join

contains

    pure subroutine multipatch_join_directions(kind,join,a,b)
        integer, intent(in) :: kind
        type(multipatch_join), intent(in) :: join
        integer, intent(out) :: a(3),b(3)

        a=boundary_axes(join%side_a)
        b=boundary_axes(join%side_b,join%swap,join%flip,join%reverse)
        if(kind<3)then
            a(kind:)=0
            b(kind:)=0
        end if
    end subroutine multipatch_join_directions

    subroutine append_assembly_join(group,kind,spec,a,b,recorded,valid)
        class(*), intent(inout) :: group
        integer, intent(in) :: kind,a,b
        type(multipatch_join), intent(in) :: spec
        logical, intent(out) :: recorded,valid

        call connect_assembly_join(group,kind,spec,a,b,recorded)
        valid=.false.
        if(recorded)valid=validate_join(group,kind,spec,a,b)
    end subroutine append_assembly_join

    subroutine inspect_native_assembly(group,valid,rational,controls,shared_controls)
        class(*), intent(in) :: group
        logical, intent(out) :: valid,rational
        integer, intent(out) :: controls,shared_controls
        integer, allocatable :: map(:)

        valid=native_is_valid(group)
        rational=.false.
        controls=0
        shared_controls=0
        if(.not.valid)return
        select type(group)
        type is(nurbs_multipatch_curve)
            rational=group%is_rational()
        type is(nurbs_multipatch_surface)
            rational=group%is_rational()
        type is(nurbs_multipatch_volume)
            rational=group%is_rational()
        end select
        call native_assembly_control_map(group,map)
        if(allocated(map))then
            controls=size(map)
            if(size(map)>0)shared_controls=maxval(map)
        end if
    end subroutine inspect_native_assembly

    subroutine native_assembly_control_map(group,map,offsets)
        class(*), intent(in) :: group
        integer, allocatable, intent(out) :: map(:)
        integer, allocatable, optional, intent(out) :: offsets(:)

        select type(group)
        type is(nurbs_multipatch_curve)
            map=group%cmp_dof_map()
            if(present(offsets))offsets=group%cmp_dof_offsets()
        type is(nurbs_multipatch_surface)
            map=group%cmp_dof_map()
            if(present(offsets))offsets=group%cmp_dof_offsets()
        type is(nurbs_multipatch_volume)
            map=group%cmp_dof_map()
            if(present(offsets))offsets=group%cmp_dof_offsets()
        end select
    end subroutine native_assembly_control_map

    subroutine native_assembly_control_correction(group,map,points,weights,new_points,new_weights,pinned,delta,origin,error)
        class(*), intent(in) :: group
        integer, intent(in) :: map(:)
        real(rk), intent(in) :: points(:,:),weights(:),new_points(:,:),new_weights(:)
        logical, intent(in) :: pinned(:)
        real(rk), allocatable, intent(out) :: delta(:,:)
        real(rk), intent(out) :: origin(3)
        character(:), allocatable, intent(out) :: error
        type(debug) :: status
        integer, allocatable :: rowptr(:),columns(:)
        real(rk), allocatable :: coefficients(:)
        real(rk) :: relative
        integer :: i

        error='ForCAD could not assemble the continuity equations for this control edit.'
        select type(group)
        type is(nurbs_multipatch_curve)
            call group%cmp_dof_constraint(rowptr,columns,coefficients,field=.true.,err=status)
        type is(nurbs_multipatch_surface)
            call group%cmp_dof_constraint(rowptr,columns,coefficients,field=.true.,err=status)
        type is(nurbs_multipatch_volume)
            call group%cmp_dof_constraint(rowptr,columns,coefficients,field=.true.,err=status)
        class default
            return
        end select
        if(.not.status%ok.or.size(rowptr)<2.or.size(columns)==0)return
        if(any(columns<1).or.any(columns>size(map)))return
        columns=map(columns)
        error='The requested control edit is outside the representable coordinate or weight range.'
        origin=points(:,findloc(pinned,.true.,dim=1))
        allocate(delta(size(pinned),4),source=0._rk)
        do i=1,size(pinned)
            if(.not.pinned(i))cycle
            relative=new_weights(i)/weights(i)
            if(.not.ieee_is_finite(relative).or.relative<=0._rk)return
            delta(i,4)=relative-1._rk
            delta(i,:3)=new_points(:,i)-points(:,i)+delta(i,4)*(new_points(:,i)-origin)
        end do
        if(.not.all(ieee_is_finite(delta)))return
        call complete_constraint_edit(rowptr,columns,coefficients,pinned,delta,error)
    end subroutine native_assembly_control_correction

    pure function multipatch_side_label(kind,side) result(label)
        integer, intent(in) :: kind,side
        character(:), allocatable :: label
        character(5), parameter :: sides(6)=[character(5)::'u_min','u_max','v_min','v_max','w_min','w_max']

        label=''
        if(kind<1.or.kind>3.or.side<1.or.side>2*kind)return
        if(kind==1)then
            label='left'
            if(side==2)label='right'
        else
            label=sides(side)
        end if
    end function multipatch_side_label

    subroutine initialize_assembly(kind,group)
        integer, intent(in) :: kind
        class(*), allocatable, intent(out) :: group

        select case(kind)
        case(1)
            allocate(nurbs_multipatch_curve::group)
        case(2)
            allocate(nurbs_multipatch_surface::group)
        case(3)
            allocate(nurbs_multipatch_volume::group)
        case default
            return
        end select
    end subroutine initialize_assembly

    subroutine append_assembly_patch(group,patch,id)
        class(*), intent(inout) :: group
        class(*), intent(in) :: patch
        integer, intent(out) :: id

        id=0
        select type(group)
        type is(nurbs_multipatch_curve)
            select type(patch)
            type is(nurbs_curve)
                call group%add_patch(patch,id)
            end select
        type is(nurbs_multipatch_surface)
            select type(patch)
            type is(nurbs_surface)
                call group%add_patch(patch,id)
            end select
        type is(nurbs_multipatch_volume)
            select type(patch)
            type is(nurbs_volume)
                call group%add_patch(patch,id)
            end select
        end select
    end subroutine append_assembly_patch

    subroutine connect_assembly_join(group,kind,spec,a,b,recorded)
        class(*), intent(inout) :: group
        integer, intent(in) :: kind,a,b
        type(multipatch_join), intent(in) :: spec
        logical, intent(out) :: recorded
        character(:), allocatable :: sa,sb
        real(rk), allocatable :: jet(:)
        integer :: existing

        recorded=.false.
        sa=multipatch_side_label(kind,spec%side_a)
        sb=multipatch_side_label(kind,spec%side_b)
        if(allocated(spec%reparameterization))jet=real(spec%reparameterization,rk)
        select type(group)
        type is(nurbs_multipatch_curve)
            existing=group%get_nconnection()
            call group%connect(a,sa,b,sb,spec%continuity,reverse=spec%reverse, &
                geometric=spec%geometric,reparameterization=jet)
            recorded=group%err%ok.and.group%get_nconnection()==existing+1
            if(.not.recorded)call group%err%reset()
        type is(nurbs_multipatch_surface)
            existing=group%get_nconnection()
            call group%connect(a,sa,b,sb,spec%continuity,reverse=spec%reverse, &
                geometric=spec%geometric,reparameterization=jet)
            recorded=group%err%ok.and.group%get_nconnection()==existing+1
            if(.not.recorded)call group%err%reset()
        type is(nurbs_multipatch_volume)
            existing=group%get_nconnection()
            call group%connect(a,sa,b,sb,spec%continuity,reverse=spec%reverse,swap=spec%swap,flip=spec%flip, &
                geometric=spec%geometric,reparameterization=jet)
            recorded=group%err%ok.and.group%get_nconnection()==existing+1
            if(.not.recorded)call group%err%reset()
        end select
    end subroutine connect_assembly_join

    logical function validate_join(group,kind,spec,a,b) result(valid)
        class(*), intent(in) :: group
        integer, intent(in) :: kind,a,b
        type(multipatch_join), intent(in) :: spec
        class(*), allocatable :: pair
        integer :: local_a,local_b
        logical :: recorded

        ! Validate this join in isolation so an unrelated failed seam does not alter its status.
        call initialize_assembly(kind,pair)
        call copy_member(group,a,pair,local_a)
        local_b=local_a
        if(b/=a)call copy_member(group,b,pair,local_b)
        valid=.false.
        if(min(local_a,local_b)<1)return
        call connect_assembly_join(pair,kind,spec,local_a,local_b,recorded)
        if(recorded)valid=native_is_valid(pair)
    end function validate_join

    subroutine copy_member(group,index,pair,id)
        class(*), intent(in) :: group
        integer, intent(in) :: index
        class(*), intent(inout) :: pair
        integer, intent(out) :: id

        id=0
        select type(group)
        type is(nurbs_multipatch_curve)
            call append_assembly_patch(pair,group%get_patch(index),id)
        type is(nurbs_multipatch_surface)
            call append_assembly_patch(pair,group%get_patch(index),id)
        type is(nurbs_multipatch_volume)
            call append_assembly_patch(pair,group%get_patch(index),id)
        end select
    end subroutine copy_member

    logical function native_is_valid(group) result(valid)
        class(*), intent(in) :: group

        valid=.false.
        select type(group)
        type is(nurbs_multipatch_curve)
            valid=group%is_valid()
        type is(nurbs_multipatch_surface)
            valid=group%is_valid()
        type is(nurbs_multipatch_volume)
            valid=group%is_valid()
        end select
    end function native_is_valid

end module studio_forcad_assembly
