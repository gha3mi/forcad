module studio_multipatch_controls
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: iso_fortran_env, only: int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use studio_document, only: studio_document_t
    use studio_forcad, only: rk, studio_native_patch, native_assembly_control_map, native_assembly_control_correction, &
        get_native_controls, replace_native_controls
    use studio_multipatch, only: multipatch_draft, multipatch_join, multipatch_result, inspect_multipatch
    use studio_multipatch_collection, only: count_multipatches, commit_multipatch_patch_edits, validate_multipatch_patch_edits
    use studio_single_patch, only: single_patch_controls, read_single_patch_controls, single_patch_revision, &
        copy_single_patch_geometry
    implicit none(type, external)
    private

    public :: multipatch_controls_t, read_document_multipatch_controls, commit_document_multipatch_controls
    public :: multipatch_control_preview_t, preview_document_multipatch_controls
    public :: joined_patch_assembly, assembly_control_id, set_multipatch_control_displacement, scale_control_weight

    type :: multipatch_controls_t
        logical :: valid=.false., smooth_joins=.false.
        integer :: kind=0
        real(c_double), allocatable :: points(:,:),weights(:)
        integer, allocatable :: copies(:),members(:),offsets(:),map(:)
        integer(int64), allocatable :: revisions(:)
        character(:), allocatable :: message
        type(multipatch_draft), private :: definition
    end type multipatch_controls_t

    type :: multipatch_control_preview_t
        integer, allocatable :: ids(:)
        real(c_double), allocatable :: points(:,:),weights(:)
        logical, allocatable :: moved(:),reweighted(:)
    end type multipatch_control_preview_t

    type :: control_edit
        real(rk), allocatable :: delta(:,:)
        real(rk) :: origin(3)=0._rk
        logical, allocatable :: moved(:),reweighted(:),adjusted(:)
    end type control_edit

contains

    integer function joined_patch_assembly(document,kind,patch_id) result(id)
        type(studio_document_t), intent(in) :: document
        integer, intent(in) :: kind,patch_id
        integer :: group

        id=0
        if(kind<1.or.kind>3.or.patch_id<1)return
        do group=1,count_multipatches(document%multipatches)
            associate(definition=>document%multipatches%groups(group)%definition)
                if(definition%kind/=kind)cycle
                if(.not.allocated(definition%joins))cycle
                if(.not.any(definition%joins%continuity>=0.and. &
                    (definition%joins%patch_a==patch_id.or.definition%joins%patch_b==patch_id)))cycle
                id=group
                return
            end associate
        end do
    end function joined_patch_assembly

    subroutine read_document_multipatch_controls(document,id,controls)
        type(studio_document_t), intent(in) :: document
        integer, intent(in) :: id
        type(multipatch_controls_t), intent(out) :: controls
        type(single_patch_controls) :: local
        type(multipatch_result) :: checked
        class(*), allocatable :: native
        integer :: i,j,first,last,global,count
        logical :: current

        controls%message='Select one multipatch to edit its control points.'
        if(id<1.or.id>count_multipatches(document%multipatches))return
        associate(group=>document%multipatches%groups(id),definition=>document%multipatches%groups(id)%definition)
            controls%kind=definition%kind
            controls%members=definition%members
            controls%definition=definition
            if(allocated(definition%joins))controls%smooth_joins=any(definition%joins%continuity>0)
            controls%revisions=[(single_patch_revision(controls%kind,controls%members(i),catalog=document%patches), &
                i=1,size(controls%members))]
            current=control_cache_matches(group%revisions,controls%revisions)
            controls%message='Resolve invalid multipatch interfaces before editing its controls.'
            if(current)then
                if(.not.allocated(group%native))return
                call native_assembly_control_map(group%native,controls%map,controls%offsets)
            else
                call inspect_multipatch(definition,checked,native,catalog=document%patches)
                if(.not.allocated(native))return
                call native_assembly_control_map(native,controls%map,controls%offsets)
            end if
        end associate
        controls%message='ForCAD could not read the multipatch controls.'
        if(.not.allocated(controls%map).or..not.allocated(controls%offsets))return
        if(size(controls%map)==0.or.size(controls%offsets)/=size(controls%members)+1)return
        if(controls%offsets(1)/=0.or.controls%offsets(size(controls%offsets))/=size(controls%map))return
        if(any(controls%map<1).or.any(controls%map>size(controls%map)))return
        count=maxval(controls%map)
        allocate(controls%points(3,count),controls%weights(count))
        allocate(controls%copies(count),source=0)
        do i=1,size(controls%members)
            first=controls%offsets(i)+1
            last=controls%offsets(i+1)
            if(first<1.or.last>size(controls%map).or.first>last)return
            call read_single_patch_controls(controls%kind,local,variant=controls%members(i),catalog=document%patches)
            if(.not.local%valid)return
            if(size(local%weights)/=last-first+1)return
            do j=first,last
                global=controls%map(j)
                if(controls%copies(global)==0)then
                    controls%points(:,global)=local%points(:,j-first+1)
                    controls%weights(global)=local%weights(j-first+1)
                end if
                controls%copies(global)=controls%copies(global)+1
            end do
        end do
        if(any(controls%copies==0))return
        controls%valid=.true.
        controls%message=''
    end subroutine read_document_multipatch_controls

    subroutine preview_document_multipatch_controls(document,id,controls,preview,error,interactive,baseline)
        type(studio_document_t), intent(in) :: document
        integer, intent(in) :: id
        type(multipatch_controls_t), intent(in) :: controls
        type(multipatch_control_preview_t), intent(out) :: preview
        character(:), allocatable, intent(out) :: error
        logical, optional, intent(in) :: interactive
        type(multipatch_controls_t), optional, target, intent(in) :: baseline
        type(multipatch_controls_t), target :: fresh
        type(multipatch_controls_t), pointer :: current
        type(control_edit) :: edit
        type(studio_native_patch), allocatable :: candidates(:)
        integer, allocatable :: targets(:)
        logical :: live

        live=.false.
        if(present(interactive))live=interactive
        if(live.and.present(baseline))then
            error='The multipatch changed. Start a new control edit.'
            if(.not.snapshot_is_current(document,id,baseline))return
            current=>baseline
        else
            call read_document_multipatch_controls(document,id,fresh)
            current=>fresh
        end if
        call prepare_control_edit(document,id,controls,current,edit,error)
        if(len(error)>0)return
        if(.not.live)then
            call stage_control_candidates(document,controls,current,edit,targets,candidates,error)
            if(len(error)>0)return
            if(allocated(targets))then
                call validate_multipatch_patch_edits(document%multipatches,id,current%kind,targets,candidates,document%patches,error)
                if(len(error)>0)return
            end if
        end if
        call describe_control_dependencies(current,edit,preview,error)
    end subroutine preview_document_multipatch_controls

    subroutine commit_document_multipatch_controls(document,id,controls,changed,error)
        type(studio_document_t), intent(inout) :: document
        integer, intent(in) :: id
        type(multipatch_controls_t), intent(in) :: controls
        logical, intent(out) :: changed
        character(:), allocatable, intent(out) :: error
        type(multipatch_controls_t) :: current
        type(control_edit) :: edit
        type(studio_native_patch), allocatable :: candidates(:)
        integer, allocatable :: targets(:)

        changed=.false.
        call read_document_multipatch_controls(document,id,current)
        call prepare_control_edit(document,id,controls,current,edit,error)
        if(len(error)>0)return
        call stage_control_candidates(document,controls,current,edit,targets,candidates,error)
        if(len(error)>0.or..not.allocated(targets))return
        call commit_multipatch_patch_edits(document%multipatches,id,current%kind,targets,candidates, &
            'multipatch control edit',document%patches,changed,error)
    end subroutine commit_document_multipatch_controls

    subroutine prepare_control_edit(document,id,controls,current,edit,error)
        type(studio_document_t), intent(in) :: document
        integer, intent(in) :: id
        type(multipatch_controls_t), intent(in) :: controls,current
        type(control_edit), intent(out) :: edit
        character(:), allocatable, intent(out) :: error
        logical, allocatable :: pinned(:)
        type(multipatch_result) :: checked
        class(*), allocatable :: native

        error='The multipatch controls are unavailable. Reopen the editor.'
        if(.not.controls%valid.or..not.allocated(controls%points).or..not.allocated(controls%weights))return
        if(.not.current%valid)then
            error=current%message
            return
        end if
        error='The multipatch changed. Reopen the control editor before applying edits.'
        if(.not.same_control_source(controls,current))return
        if(any(shape(controls%points)/=shape(current%points)).or.size(controls%weights)/=size(current%weights))return

        error='Enter finite coordinates and weights greater than zero.'
        if(.not.all(ieee_is_finite(real(controls%points,rk))).or..not.all(ieee_is_finite(real(controls%weights,rk))))return
        if(any(real(controls%weights,rk)<=0._rk))return
        edit%moved=any(real(controls%points,rk)/=real(current%points,rk),dim=1)
        edit%reweighted=real(controls%weights,rk)/=real(current%weights,rk)
        pinned=edit%moved.or.edit%reweighted
        allocate(edit%adjusted(size(pinned)),source=.false.)
        error=''
        if(.not.current%smooth_joins.or..not.any(pinned))return

        associate(group=>document%multipatches%groups(id))
            error='Resolve invalid multipatch interfaces before editing its controls.'
            if(control_cache_matches(group%revisions,current%revisions))then
                if(.not.allocated(group%native))return
                call correct(group%native)
            else
                call inspect_multipatch(group%definition,checked,native,catalog=document%patches)
                if(.not.allocated(native))return
                call correct(native)
            end if
        end associate
        if(len(error)>0)return
        edit%adjusted=any(edit%delta/=0._rk,dim=2).and..not.pinned
        error='Keeping the joins continuous would produce a nonpositive or unrepresentable dependent weight.'
        if(any(edit%adjusted.and.1._rk+edit%delta(:,4)<=0._rk))return
        error=''
    contains
        subroutine correct(group)
            class(*), intent(in) :: group

            call native_assembly_control_correction(group,current%map,real(current%points,rk),real(current%weights,rk), &
                real(controls%points,rk),real(controls%weights,rk),pinned,edit%delta,edit%origin,error)
        end subroutine correct
    end subroutine prepare_control_edit

    subroutine describe_control_dependencies(current,edit,preview,error)
        type(multipatch_controls_t), intent(in) :: current
        type(control_edit), intent(in) :: edit
        type(multipatch_control_preview_t), intent(out) :: preview
        character(:), allocatable, intent(out) :: error
        type(multipatch_control_preview_t) :: rows
        real(rk) :: point(3),weight
        logical, allocatable :: changed(:)
        integer, allocatable :: indices(:)
        integer :: i,n

        n=count(edit%adjusted)
        rows%ids=pack([(i,i=1,size(edit%adjusted))],edit%adjusted)
        allocate(rows%points(3,n),rows%weights(n),rows%moved(n),rows%reweighted(n))
        do i=1,n
            point=real(current%points(:,rows%ids(i)),rk)
            weight=real(current%weights(rows%ids(i)),rk)
            call adjust_dependent_control(point,weight,edit%delta(rows%ids(i),:),edit%origin)
            rows%points(:,i)=real(point,c_double)
            rows%weights(i)=real(weight,c_double)
            rows%moved(i)=any(rows%points(:,i)/=current%points(:,rows%ids(i)))
            rows%reweighted(i)=rows%weights(i)/=current%weights(rows%ids(i))
        end do
        error='The continuity edit exceeds the representable coordinate or weight range.'
        if(.not.all(ieee_is_finite(rows%points)).or..not.all(ieee_is_finite(rows%weights)).or.any(rows%weights<=0._c_double))return
        changed=rows%moved.or.rows%reweighted
        if(all(changed))then
            preview=rows
        else
            indices=pack([(i,i=1,n)],changed)
            preview%ids=rows%ids(indices)
            preview%points=rows%points(:,indices)
            preview%weights=rows%weights(indices)
            preview%moved=rows%moved(indices)
            preview%reweighted=rows%reweighted(indices)
        end if
        error=''
    end subroutine describe_control_dependencies

    subroutine stage_control_candidates(document,controls,current,edit,targets,candidates,error)
        type(studio_document_t), intent(in) :: document
        type(multipatch_controls_t), intent(in) :: controls,current
        type(control_edit), intent(in) :: edit
        integer, allocatable, intent(out) :: targets(:)
        type(studio_native_patch), allocatable, intent(out) :: candidates(:)
        character(:), allocatable, intent(out) :: error
        integer, allocatable :: patch_map(:)
        real(rk), allocatable :: points(:,:),weights(:)
        logical, allocatable :: affected(:),active(:)
        logical :: points_changed,weights_changed
        integer :: i,j,global,n

        error=''
        active=edit%moved.or.edit%reweighted.or.edit%adjusted
        if(.not.any(active))return
        allocate(affected(size(current%members)))
        do i=1,size(affected)
            patch_map=current%map(current%offsets(i)+1:current%offsets(i+1))
            affected(i)=any(active(patch_map))
        end do
        targets=pack(current%members,affected)
        allocate(candidates(size(targets)))
        n=0
        do i=1,size(affected)
            if(.not.affected(i))cycle
            n=n+1
            patch_map=current%map(current%offsets(i)+1:current%offsets(i+1))
            error='A member patch no longer has the expected controls.'
            call copy_single_patch_geometry(current%kind,targets(n),candidates(n)%geometry,catalog=document%patches)
            if(.not.allocated(candidates(n)%geometry))return
            if(.not.get_native_controls(candidates(n)%geometry,size(patch_map),points,weights))return
            points_changed=any(edit%moved(patch_map).or.edit%adjusted(patch_map))
            weights_changed=any(edit%reweighted(patch_map))
            do j=1,size(patch_map)
                global=patch_map(j)
                if(edit%adjusted(global))then
                    call adjust_dependent_control(points(j,:),weights(j),edit%delta(global,:),edit%origin)
                    weights_changed=weights_changed.or.1._rk+edit%delta(global,4)/=1._rk
                else
                    if(edit%moved(global))points(j,:)=real(controls%points(:,global),rk)
                    if(edit%reweighted(global))weights(j)=scale_control_weight(weights(j),real(current%weights(global),rk), &
                        real(controls%weights(global),rk))
                end if
            end do
            error='The continuity edit exceeds the representable coordinate or weight range.'
            if(.not.all(ieee_is_finite(points)).or..not.all(ieee_is_finite(weights)).or.any(weights<=0._rk))return
            error='ForCAD rejected the control points or weights.'
            if(.not.replace_native_controls(candidates(n)%geometry,points,weights,points_changed,weights_changed))return
        end do
        error=''
    end subroutine stage_control_candidates

    pure subroutine adjust_dependent_control(point,weight,delta,origin)
        real(rk), intent(inout) :: point(3),weight
        real(rk), intent(in) :: delta(4),origin(3)
        real(rk) :: factor

        factor=1._rk+delta(4)
        point=point+(delta(:3)-delta(4)*(point-origin))/factor
        weight=weight*factor
    end subroutine adjust_dependent_control

    logical function snapshot_is_current(document,id,baseline) result(current)
        type(studio_document_t), intent(in) :: document
        integer, intent(in) :: id
        type(multipatch_controls_t), intent(in) :: baseline
        integer :: i

        current=.false.
        if(id<1.or.id>count_multipatches(document%multipatches).or..not.baseline%valid)return
        if(.not.same_control_definition(baseline%definition,document%multipatches%groups(id)%definition))return
        do i=1,size(baseline%members)
            if(single_patch_revision(baseline%kind,baseline%members(i),catalog=document%patches)/=baseline%revisions(i))return
        end do
        current=.true.
    end function snapshot_is_current

    logical function control_cache_matches(cached,revisions) result(current)
        integer(int64), allocatable, intent(in) :: cached(:)
        integer(int64), intent(in) :: revisions(:)

        current=.false.
        if(.not.allocated(cached))return
        if(size(cached)==size(revisions))current=all(cached==revisions)
    end function control_cache_matches

    integer function assembly_control_id(controls,patch_id,point) result(id)
        type(multipatch_controls_t), intent(in) :: controls
        integer, intent(in) :: patch_id,point
        integer :: member

        id=0
        if(.not.controls%valid.or.point<1)return
        member=findloc(controls%members,patch_id,dim=1)
        if(member<1)return
        if(point>controls%offsets(member+1)-controls%offsets(member))return
        id=controls%map(controls%offsets(member)+point)
    end function assembly_control_id

    subroutine set_multipatch_control_displacement(controls,point,position,selected)
        type(multipatch_controls_t), intent(inout) :: controls
        integer, intent(in) :: point
        real(c_double), intent(in) :: position(3)
        integer, optional, intent(in) :: selected(:)
        real(c_double) :: displacement(3)
        integer :: i,id

        if(.not.controls%valid.or.point<1)return
        if(point>size(controls%points,2))return
        displacement=position-controls%points(:,point)
        if(present(selected))then
            if(any(selected==point))then
                do i=1,size(selected)
                    id=selected(i)
                    if(id<1.or.id>size(controls%points,2).or.id==point)cycle
                    controls%points(:,id)=controls%points(:,id)+displacement
                end do
            end if
        end if
        controls%points(:,point)=position
    end subroutine set_multipatch_control_displacement

    logical function same_control_source(edited,current) result(same)
        type(multipatch_controls_t), intent(in) :: edited,current
        same=.false.
        if(.not.allocated(edited%revisions).or..not.allocated(edited%map))return
        if(size(edited%revisions)/=size(current%revisions).or.size(edited%map)/=size(current%map))return
        if(any(edited%revisions/=current%revisions).or.any(edited%map/=current%map))return
        same=same_control_definition(edited%definition,current%definition)
    end function same_control_source

    logical function same_control_definition(a,b) result(same)
        type(multipatch_draft), intent(in) :: a,b
        integer :: i

        same=.false.
        if(a%name/=b%name.or.a%kind/=b%kind)return
        if(size(a%members)/=size(b%members))return
        if(any(a%members/=b%members))return
        if(allocated(a%joins).neqv.allocated(b%joins))return
        if(allocated(a%joins))then
            if(size(a%joins)/=size(b%joins))return
            do i=1,size(a%joins)
                if(.not.same_join(a%joins(i),b%joins(i)))return
            end do
        end if
        same=.true.
    end function same_control_definition

    logical function same_join(a,b) result(same)
        type(multipatch_join), intent(in) :: a,b

        same=.false.
        if(any([a%patch_a,a%side_a,a%patch_b,a%side_b,a%continuity]/= &
            [b%patch_a,b%side_a,b%patch_b,b%side_b,b%continuity]))return
        if(any([a%reverse,a%swap,a%flip,a%geometric].neqv.[b%reverse,b%swap,b%flip,b%geometric]))return
        if(allocated(a%reparameterization).neqv.allocated(b%reparameterization))return
        if(allocated(a%reparameterization))then
            if(size(a%reparameterization)/=size(b%reparameterization))return
            if(any(a%reparameterization/=b%reparameterization))return
        end if
        same=.true.
    end function same_join

    pure real(rk) function scale_control_weight(local,representative,updated) result(weight)
        real(rk), intent(in) :: local,representative,updated
        real(rk) :: ratio
        integer :: power

        if(local==representative)then
            weight=updated
            return
        end if
        ratio=fraction(local)*fraction(updated)/fraction(representative)
        power=exponent(local)+exponent(updated)-exponent(representative)+exponent(ratio)
        weight=0._rk
        if(power>maxexponent(weight).or.power<minexponent(weight)-digits(weight))return
        weight=scale(fraction(ratio),power)
    end function scale_control_weight
end module studio_multipatch_controls
