module studio_multipatch_refinement
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: iso_fortran_env, only: int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use studio_document, only: studio_document_t
    use studio_forcad, only: rk,studio_native_patch,multipatch_join_directions
    use studio_multipatch, only: multipatch_draft,multipatch_result,inspect_multipatch
    use studio_multipatch_collection, only: count_multipatches,commit_multipatch_patch_edits,validate_multipatch_patch_edits
    use studio_single_patch, only: read_single_patch_bounds,stage_native_patch_refinement,single_patch_name,single_patch_revision
    implicit none(type, external)
    private

    public :: commit_multipatch_refinement,multipatch_refinement_bounds,plan_multipatch_refinement

    type, public :: multipatch_refinement_row_t
        integer :: member=0,axis=0
        real(c_double) :: knot=0._c_double
    end type multipatch_refinement_row_t

    type, public :: multipatch_refinement_plan_t
        integer :: kind=0,assembly=0,steps=0
        logical :: inserting=.false.
        type(multipatch_refinement_row_t), allocatable :: rows(:)
    end type multipatch_refinement_plan_t

contains

    logical function multipatch_refinement_bounds(document,assembly,member,axis,bounds) result(valid)
        type(studio_document_t), intent(in) :: document
        integer, intent(in) :: assembly,member,axis
        real(c_double), intent(out) :: bounds(2)
        real(rk) :: ranges(2,3)

        valid=.false.
        bounds=0._c_double
        if(assembly<1.or.assembly>count_multipatches(document%multipatches))return
        associate(definition=>document%multipatches%groups(assembly)%definition)
            if(axis<1.or.axis>definition%kind)return
            if(.not.any(definition%members==member))return
            if(.not.read_single_patch_bounds(definition%kind,member,document%patches,ranges))return
        end associate
        bounds=real(ranges(:,axis),c_double)
        valid=all(ieee_is_finite(bounds)).and.bounds(1)<bounds(2)
    end function multipatch_refinement_bounds

    subroutine plan_multipatch_refinement(document,assembly,member,axis,steps,plan,error,knot)
        type(studio_document_t), intent(in) :: document
        integer, intent(in) :: assembly,member,axis,steps
        type(multipatch_refinement_plan_t), intent(out) :: plan
        character(:), allocatable, intent(out) :: error
        real(c_double), intent(in), optional :: knot
        type(multipatch_refinement_plan_t) :: staged
        type(studio_native_patch), allocatable :: candidates(:)
        integer, allocatable :: targets(:)

        call stage_multipatch_refinement(document,assembly,member,axis,steps,staged,targets,candidates,error,knot)
        if(len(error)>0)return
        call validate_multipatch_patch_edits(document%multipatches,assembly,staged%kind,targets,candidates,document%patches,error)
        if(len(error)==0)plan=staged
    end subroutine plan_multipatch_refinement

    subroutine commit_multipatch_refinement(document,assembly,member,axis,steps,error,changed,knot)
        type(studio_document_t), intent(inout) :: document
        integer, intent(in) :: assembly,member,axis,steps
        character(:), allocatable, intent(out) :: error
        logical, intent(out) :: changed
        real(c_double), intent(in), optional :: knot
        type(multipatch_refinement_plan_t) :: plan
        type(studio_native_patch), allocatable :: candidates(:)
        integer, allocatable :: targets(:)

        changed=.false.
        call stage_multipatch_refinement(document,assembly,member,axis,steps,plan,targets,candidates,error,knot)
        if(len(error)>0)return
        call commit_multipatch_patch_edits(document%multipatches,assembly,plan%kind,targets,candidates, &
            trim(merge('multipatch knot insertion  ','multipatch degree elevation',plan%inserting)), &
            document%patches,changed,error,preserve_mesh=.true.)
    end subroutine commit_multipatch_refinement

    subroutine stage_multipatch_refinement(document,assembly,member,axis,steps,plan,targets,candidates,error,knot)
        type(studio_document_t), intent(in) :: document
        integer, intent(in) :: assembly,member,axis,steps
        type(multipatch_refinement_plan_t), intent(out) :: plan
        integer, allocatable, intent(out) :: targets(:)
        type(studio_native_patch), allocatable, intent(out) :: candidates(:)
        character(:), allocatable, intent(out) :: error
        real(c_double), intent(in), optional :: knot
        integer, allocatable :: orientations(:,:)
        integer :: kind,i,j,first,index,row
        real(rk) :: fraction,value,width,bounds(2,3),mapped_fraction,knot_values(3)
        logical :: directions(3),mapped(3),check_orientation

        error='Select a multipatch member and one of its parameter directions.'
        if(assembly<1.or.assembly>count_multipatches(document%multipatches))return
        associate(definition=>document%multipatches%groups(assembly)%definition)
            kind=definition%kind
            if(axis<1.or.axis>kind.or..not.allocated(definition%members))return
            first=findloc(definition%members,member,dim=1)
            if(first<1)return
        end associate
        error='Enter a positive refinement count.'
        if(steps<1)return
        error='Validate the multipatch and resolve its interfaces before refining it.'
        if(.not.current_multipatch_valid(document,assembly))return

        value=0._rk
        fraction=0.5_rk
        if(present(knot))then
            error='Enter a finite knot strictly inside the member parameter interval.'
            if(.not.ieee_is_finite(knot))return
            if(.not.read_single_patch_bounds(kind,member,document%patches,bounds))return
            value=real(knot,rk)
            if(.not.ieee_is_finite(value))return
            if(value<=bounds(1,axis).or.value>=bounds(2,axis))return
            width=bounds(2,axis)-bounds(1,axis)
            if(.not.ieee_is_finite(width))return
            fraction=(value-bounds(1,axis))/width
            if(fraction<=0._rk.or.fraction>=1._rk)return
        end if
        check_orientation=present(knot).and.fraction/=0.5_rk
        plan%assembly=assembly
        plan%kind=kind
        plan%steps=steps
        plan%inserting=present(knot)

        associate(definition=>document%multipatches%groups(assembly)%definition)
            if(.not.refinement_directions(definition,first,axis,check_orientation,orientations,error))return
            targets=pack(definition%members,any(orientations/=0,dim=1))
            allocate(candidates(size(targets)),plan%rows(count(orientations/=0)))
            index=0
            row=0
            do i=1,size(definition%members)
                directions=orientations(:,i)/=0
                if(.not.any(directions))cycle
                index=index+1
                knot_values=0._rk
                if(plan%inserting)then
                    error='Could not read the parameter interval of '//single_patch_name(kind,targets(index),catalog=document%patches)//'.'
                    if(.not.read_single_patch_bounds(kind,targets(index),document%patches,bounds))return
                    mapped=directions
                    do j=1,kind
                        if(.not.directions(j))cycle
                        width=bounds(2,j)-bounds(1,j)
                        if(.not.ieee_is_finite(width).or.width<=0._rk)return
                        mapped_fraction=merge(fraction,1._rk-fraction,orientations(j,i)>0)
                        knot_values(j)=bounds(1,j)+width*mapped_fraction
                    end do
                    if(i==first)then
                        knot_values(axis)=value
                        mapped(axis)=.false.
                    end if
                    call stage_native_patch_refinement(kind,targets(index),directions,steps,document%patches, &
                        candidates(index),error,knot_values,mapped)
                else
                    call stage_native_patch_refinement(kind,targets(index),directions,steps,document%patches, &
                        candidates(index),error,mapped=directions)
                end if
                if(len(error)>0)then
                    error=single_patch_name(kind,targets(index),catalog=document%patches)//' · '//error
                    return
                end if
                do j=1,kind
                    if(.not.directions(j))cycle
                    row=row+1
                    plan%rows(row)=multipatch_refinement_row_t(targets(index),j,real(knot_values(j),c_double))
                end do
            end do
        end associate
    end subroutine stage_multipatch_refinement

    logical function current_multipatch_valid(document,assembly) result(valid)
        type(studio_document_t), intent(in) :: document
        integer, intent(in) :: assembly
        type(multipatch_result) :: checked
        integer(int64), allocatable :: revisions(:)
        integer :: i

        associate(group=>document%multipatches%groups(assembly),definition=>document%multipatches%groups(assembly)%definition)
            if(allocated(group%revisions))then
                revisions=[(single_patch_revision(definition%kind,definition%members(i),catalog=document%patches), &
                    i=1,size(definition%members))]
                if(size(group%revisions)==size(revisions))then
                    if(all(group%revisions==revisions))then
                        valid=group%validation%valid
                        return
                    end if
                end if
            end if
            call inspect_multipatch(definition,checked,catalog=document%patches)
            valid=checked%valid
        end associate
    end function current_multipatch_valid

    logical function refinement_directions(definition,member,axis,check_orientation,directions,error) result(valid)
        type(multipatch_draft), intent(in) :: definition
        integer, intent(in) :: member,axis
        logical, intent(in) :: check_orientation
        integer, allocatable, intent(out) :: directions(:,:)
        character(:), allocatable, intent(out) :: error
        integer, allocatable :: endpoints(:,:),orientation(:)
        integer :: join,a,b,i,n,left_axes(3),right_axes(3),left,right,left_sign,right_sign,reflection
        logical :: expanded

        valid=.false.
        error='An interface does not identify valid member directions.'
        allocate(directions(3,size(definition%members)),source=0)
        directions(axis,member)=1
        if(allocated(definition%joins).and.definition%kind>1)then
            allocate(endpoints(2,(definition%kind-1)*size(definition%joins)))
            allocate(orientation(size(endpoints,dim=2)))
            n=0
            do join=1,size(definition%joins)
                associate(seam=>definition%joins(join))
                    if(seam%continuity<0)cycle
                    a=findloc(definition%members,seam%patch_a,dim=1)
                    b=findloc(definition%members,seam%patch_b,dim=1)
                    if(min(a,b)<1)return
                    call multipatch_join_directions(definition%kind,seam,left_axes,right_axes)
                    do i=1,definition%kind-1
                        left=left_axes(i)
                        right=abs(right_axes(i))
                        if(min(left,right)<1.or.max(left,right)>definition%kind)return
                        n=n+1
                        endpoints(:,n)=[3*(a-1)+left,3*(b-1)+right]
                        orientation(n)=sign(1,right_axes(i))
                    end do
                end associate
            end do
            error='This knot maps to both orientations around a closed chain of joins. Use the midpoint or separate the joins.'
            do
                expanded=.false.
                do join=1,n
                    left=endpoints(1,join)
                    right=endpoints(2,join)
                    reflection=orientation(join)
                    left_sign=direction_at(left)
                    right_sign=direction_at(right)
                    if(left_sign==0.and.right_sign==0)cycle
                    if(left_sign==0)then
                        call set_direction(left,reflection*right_sign)
                        expanded=.true.
                    else if(right_sign==0)then
                        call set_direction(right,reflection*left_sign)
                        expanded=.true.
                    else if(check_orientation.and.left_sign/=reflection*right_sign)then
                        return
                    end if
                end do
                if(.not.expanded)exit
            end do
        end if
        valid=.true.
        error=''
    contains
        integer function direction_at(node) result(value)
            integer, intent(in) :: node

            value=directions(mod(node-1,3)+1,(node-1)/3+1)
        end function direction_at

        subroutine set_direction(node,value)
            integer, intent(in) :: node,value

            directions(mod(node-1,3)+1,(node-1)/3+1)=value
        end subroutine set_direction
    end function refinement_directions

end module studio_multipatch_refinement
