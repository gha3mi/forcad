submodule (studio_single_patch) studio_single_patch_edits
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use studio_forcad, only: parameter_parts, get_native_controls, replace_native_controls, replace_native_knot_values, &
        edit_native_knots, insert_native_knots, edit_native_degree, reparameterize_patch, reverse_patch
    implicit none(type, external)

contains

    module procedure commit_single_patch_controls
        type(patch_entry) :: candidate
        real(rk), allocatable :: old_points(:,:),old_weights(:),points(:,:),weights(:)
        logical :: points_changed,weights_changed
        changed=.false.
        error='The controls do not match the selected patch.'
        if (kind<1.or.kind>3) return
        if (.not.controls%valid.or.controls%dimension/=kind) return
        if (.not.allocated(controls%points).or..not.allocated(controls%weights)) return
        if (.not.resolve_patch_snapshot(kind,candidate,patch_id,catalog)) return
        if (.not.get_native_controls(candidate%geometry,product(candidate%controls),old_points,old_weights)) return
        if (any(controls%counts/=candidate%controls)) return
        if (any(shape(controls%points)/=[3,size(old_points,1)]).or.size(controls%weights)/=size(old_weights)) return
        error='Enter finite coordinates and weights greater than zero.'
        if (.not.all(ieee_is_finite(controls%points)).or..not.all(ieee_is_finite(controls%weights))) return
        if (any(controls%weights<=0.0_c_double)) return
        points=transpose(real(controls%points,rk))
        weights=real(controls%weights,rk)
        if (.not.all(ieee_is_finite(points)).or..not.all(ieee_is_finite(weights))) return
        if (any(weights<=0._rk)) return
        points_changed=any(points/=old_points)
        weights_changed=any(weights/=old_weights)
        error=''
        if (.not.points_changed.and..not.weights_changed) return

        if (.not.replace_native_controls(candidate%geometry,points,weights,points_changed,weights_changed)) then
            error='ForCAD rejected the control edit.'
            return
        end if
        error='ForCAD could not sample the edited patch.'
        changed=publish_edit(candidate,'control edit',0,catalog,error)
    end procedure commit_single_patch_controls

    module procedure commit_single_patch_knot
        type(patch_entry) :: candidate
        real(rk) :: value,lo,hi
        integer :: multiplicity,controls,applied,difference
        character(:), allocatable :: operation

        accepted=0
        error='Choose a parameter direction of the current patch.'
        if (kind<1.or.kind>3.or.axis<1.or.axis>kind) return
        error='Enter a finite knot value and a positive count.'
        if (.not.ieee_is_finite(knot).or.repeats<1) return
        error='ForCAD could not load the patch.'
        if (.not.resolve_patch_snapshot(kind,candidate,patch_id,catalog)) return
        value=real(knot,rk)
        lo=candidate%bounds(1,axis)
        hi=candidate%bounds(2,axis)
        error='Enter a knot in the active parameter interval.'
        if (.not.ieee_is_finite(value).or.value<lo.or.value>hi) return
        multiplicity=count(candidate%axes(axis)%knots==value)
        controls=candidate%controls(axis)

        if (remove) then
            operation='knot removal'
            error='Select a stored knot strictly inside the active interval.'
            if (value<=lo.or.value>=hi.or.multiplicity==0) return
            error='The removal count must fit the knot multiplicity.'
            if (repeats>multiplicity) return
        else
            operation='knot insertion'
            error='The requested insertion exceeds the knot multiplicity limit.'
            if (repeats>candidate%degree(axis)+1-multiplicity) return
        end if

        error='ForCAD rejected the '//operation//'.'
        if (.not.edit_native_knots(candidate%geometry,axis,value,repeats,remove,applied)) return
        if (remove.and.applied==0) then
            error='ForCAD could not remove this knot while preserving the patch geometry.'
            return
        end if
        if (applied<1.or.applied>repeats) return
        if (.not.refresh_patch_metadata(candidate)) return

        difference=merge(-applied,applied,remove)
        if (candidate%controls(axis)/=controls+difference) return
        if (count(candidate%axes(axis)%knots==value)/=multiplicity+difference) return
        error='ForCAD could not sample the edited patch.'
        if (.not.publish_edit(candidate,operation,axis,catalog,error)) return
        accepted=applied
        error=''
    end procedure commit_single_patch_knot

    module procedure commit_single_patch_degree
        type(patch_entry) :: candidate
        integer :: applied,expected_degree(3)
        real(rk) :: tol,bound,bounds(2,3)
        character(:), allocatable :: operation

        accepted=0
        error_bound=0.0_c_double
        error='Choose a parameter direction of the current patch.'
        if (kind<1.or.kind>3.or.axis<1.or.axis>kind) return
        error='Enter a positive number of degree steps.'
        if (steps<1) return
        error='ForCAD could not load the patch.'
        if (.not.resolve_patch_snapshot(kind,candidate,patch_id,catalog)) return
        expected_degree=candidate%degree
        bounds=candidate%bounds

        tol=0._rk
        if (reduce) then
            operation='degree reduction'
            error='The reduction cannot go below degree zero.'
            if (steps>expected_degree(axis)) return
            error='Enter a finite, nonnegative tolerance in model units.'
            if (.not.ieee_is_finite(tolerance).or.tolerance<0.0_c_double) return
            tol=real(tolerance,rk)
            if (.not.ieee_is_finite(tol)) return
        else
            operation='degree elevation'
            error='The requested degree is too large.'
            if (steps>=huge(steps)-expected_degree(axis)) return
            error='The elevated patch would exceed the supported index range.'
            if (.not.refinement_fits_indices(candidate,axis,1,steps)) return
        end if

        error='ForCAD rejected the '//operation//'.'
        if (.not.edit_native_degree(candidate%geometry,axis,steps,reduce,tol,applied,bound)) return
        if (reduce.and.applied==0) then
            error='ForCAD could not reduce the degree within this tolerance.'
            return
        end if
        if (applied<1.or.applied>steps) return
        if (.not.ieee_is_finite(bound).or.bound<0._rk.or.bound>tol) return
        if (.not.ieee_is_finite(real(bound,c_double))) return
        if (.not.refresh_patch_metadata(candidate)) return
        expected_degree(axis)=expected_degree(axis)+merge(-applied,applied,reduce)
        if (any(candidate%degree/=expected_degree).or.any(candidate%bounds/=bounds)) return
        error='ForCAD could not sample the edited patch.'
        if (.not.publish_edit(candidate,operation,axis,catalog,error)) return
        accepted=applied
        error_bound=real(bound,c_double)
        error=''
    end procedure commit_single_patch_degree

    module procedure read_single_patch_bounds
        type(patch_entry) :: snapshot

        bounds=0._rk
        valid=copy_patch_snapshot(kind,variant,snapshot,catalog)
        if(valid)bounds=snapshot%bounds
    end procedure read_single_patch_bounds

    module procedure stage_native_patch_refinement
        type(patch_entry) :: candidate
        real(rk) :: bounds(2,3)
        integer :: axis
        character(3), parameter :: axes='UVW'

        error='ForCAD could not load the member patch.'
        if(.not.copy_patch_snapshot(kind,variant,candidate,catalog))return
        bounds=candidate%bounds
        do axis=1,kind
            if(.not.directions(axis))cycle
            if(present(knot_values))then
                call stage_member_knot(candidate,axis,knot_values(axis),mapped(axis),steps,error)
            else
                call stage_member_elevation(candidate,axis,steps,error)
            end if
            if(len(error)>0)then
                error=axes(axis:axis)//' · '//error
                return
            end if
        end do
        error='The refinement changed the active parameter interval.'
        if(any(candidate%bounds/=bounds))return
        call move_alloc(candidate%geometry,staged%geometry)
        error=''
    end procedure stage_native_patch_refinement

    subroutine stage_member_knot(candidate,axis,knot,mapped,steps,error)
        type(patch_entry), intent(inout) :: candidate
        integer, intent(in) :: axis,steps
        real(rk), intent(inout) :: knot
        logical, intent(in) :: mapped
        character(:), allocatable, intent(out) :: error
        real(rk) :: value,nearest
        integer(int64) :: projected
        integer :: multiplicity,controls,remaining_count,direction

        value=knot
        error='The mapped knot must be inside the active interval.'
        if(.not.ieee_is_finite(value))return
        if(value<=candidate%bounds(1,axis).or.value>=candidate%bounds(2,axis))return
        associate(knots=>candidate%axes(axis)%knots)
            if(mapped)then
                nearest=knots(minloc(abs(knots-value),dim=1))
                if(abs(nearest-value)<=2._rk*max(spacing(nearest),spacing(value)))then
                    if(nearest>candidate%bounds(1,axis).and.nearest<candidate%bounds(2,axis))value=nearest
                end if
            end if
            multiplicity=count(knots==value)
        end associate
        error='The requested insertion exceeds the knot multiplicity limit.'
        if(steps>candidate%degree(axis)+1-multiplicity)return
        controls=candidate%controls(axis)
        error='The refined patch would exceed the supported index range.'
        remaining_count=huge(remaining_count)
        do direction=1,candidate%dimension
            projected=int(candidate%controls(direction),int64)+merge(steps,0,direction==axis)
            if(projected>remaining_count.or.projected+candidate%degree(direction)+1_int64>huge(controls))return
            remaining_count=remaining_count/int(projected)
        end do
        error='ForCAD rejected the knot insertion.'
        if(.not.insert_native_knots(candidate%geometry,axis,[value],[steps]))return
        if(.not.refresh_patch_metadata(candidate))return
        if(candidate%controls(axis)/=controls+steps)return
        if(count(candidate%axes(axis)%knots==value)/=multiplicity+steps)return
        knot=value
        error=''
    end subroutine stage_member_knot

    subroutine stage_member_elevation(candidate,axis,steps,error)
        type(patch_entry), intent(inout) :: candidate
        integer, intent(in) :: axis,steps
        character(:), allocatable, intent(out) :: error
        integer :: applied,expected_degree(3)
        real(rk) :: bound

        error='The elevated patch would exceed the supported index range.'
        if(.not.refinement_fits_indices(candidate,axis,1,steps))return
        expected_degree=candidate%degree
        expected_degree(axis)=expected_degree(axis)+steps
        error='ForCAD rejected the degree elevation.'
        if(.not.edit_native_degree(candidate%geometry,axis,steps,.false.,0._rk,applied,bound))return
        if(.not.refresh_patch_metadata(candidate))return
        if(applied/=steps.or.any(candidate%degree/=expected_degree))return
        error=''
    end subroutine stage_member_elevation

    module procedure commit_single_patch_refinement
        type(patch_entry) :: candidate
        real(rk), allocatable :: insertions(:)
        real(rk) :: bounds(2,3),unused_bound
        integer :: direction,last,old_degree,old_count,applied,expected_degree
        character(3), parameter :: axes='UVW'
        character :: name

        changed=.false.
        error='Choose one direction or all directions of the current patch.'
        if (kind<1.or.kind>3.or.axis<0.or.axis>kind) return
        error='Enter a positive span division count and a nonnegative degree increase.'
        if (subdivisions<1.or.degree_steps<0) return
        error=''
        if (subdivisions==1.and.degree_steps==0) return
        error='ForCAD could not load the patch.'
        if (.not.resolve_patch_snapshot(kind,candidate,patch_id,catalog)) return
        error='The refined patch would exceed the supported index range.'
        if (.not.refinement_fits_indices(candidate,axis,subdivisions,degree_steps)) return
        bounds=candidate%bounds
        last=merge(kind,axis,axis==0)

        ! Refine a private copy so all selected directions share one transaction.
        do direction=max(1,axis),last
            name=axes(direction:direction)
            old_degree=candidate%degree(direction)
            error='The requested '//name//' degree is too large.'
            if (degree_steps>=huge(old_degree)-old_degree) return
            expected_degree=old_degree+degree_steps
            if (degree_steps>0) then
                error='ForCAD could not elevate degree along '//name//'.'
                if (.not.edit_native_degree(candidate%geometry,direction,degree_steps,.false.,0._rk,applied,unused_bound)) return
                if (.not.refresh_patch_metadata(candidate)) return
                if (applied/=degree_steps.or.candidate%degree(direction)/=expected_degree) return
            end if
            if (subdivisions>1) then
                error='The '//name//' spans cannot represent the requested subdivision.'
                if (.not.span_subdivisions(candidate,direction,subdivisions,insertions)) return
                error='ForCAD could not subdivide the '//name//' spans.'
                old_count=candidate%controls(direction)
                if (.not.insert_native_knots(candidate%geometry,direction,insertions,spread(1,1,size(insertions)))) return
                if (.not.refresh_patch_metadata(candidate)) return
                if (candidate%controls(direction)/=old_count+size(insertions)) return
                if (candidate%degree(direction)/=expected_degree) return
            end if
        end do
        error='The refinement changed the active parameter interval.'
        if (any(candidate%bounds/=bounds)) return
        error='ForCAD could not sample the refined patch.'
        changed=publish_edit(candidate,'refinement',axis,catalog,error)
        if (changed) error=''
    end procedure commit_single_patch_refinement

    module procedure single_patch_refinement_limits
        type(patch_entry) :: patch

        if(present(max_divisions))max_divisions=1
        if(present(max_steps))max_steps=0
        if(kind<1.or.kind>3.or.axis<0.or.axis>kind)return
        if(.not.resolve_patch_snapshot(kind,patch,patch_id,catalog))return
        if(present(max_divisions))max_divisions=1+refinement_limit(patch,axis,degree_increase=.false.)
        if(present(max_steps))max_steps=refinement_limit(patch,axis,degree_increase=.true.)
    end procedure single_patch_refinement_limits

    integer function refinement_limit(patch,axis,degree_increase) result(limit)
        type(patch_entry), intent(in) :: patch
        integer, intent(in) :: axis
        logical, intent(in) :: degree_increase
        integer :: upper,midpoint,divisions,steps

        limit=0
        upper=huge(upper)
        do while(upper-limit>1)
            midpoint=limit+(upper-limit)/2
            divisions=1+merge(0,midpoint,degree_increase)
            steps=merge(midpoint,0,degree_increase)
            if(refinement_fits_indices(patch,axis,divisions,steps))then
                limit=midpoint
            else
                upper=midpoint
            end if
        end do
    end function refinement_limit

    logical function refinement_fits_indices(patch,axis,divisions,steps) result(fits)

        type(patch_entry), intent(in) :: patch
        integer, intent(in) :: axis,divisions,steps
        integer :: direction,remaining,count_,degree,spans,elevation_spans
        integer(int64) :: projected

        fits=.false.
        remaining=huge(remaining)
        do direction=1,patch%dimension
            count_=patch%controls(direction)
            if(count_<1)return
            degree=patch%degree(direction)
            projected=count_
            if (axis==0.or.axis==direction) then
                associate(knots=>patch%axes(direction)%knots)
                    spans=count(knots(degree+1:count_)<knots(degree+2:count_+1))
                    elevation_spans=spans
                    if (any(knots(1:degree+1)/=knots(degree+1)).or. &
                        any(knots(count_+1:)/=knots(count_+1))) elevation_spans=spans+1
                end associate
                ! Periodic elevation also extends the repeated control rows by the degree increment.
                projected=projected+int(spans,int64)*(divisions-1)+int(elevation_spans,int64)*steps
            end if
            if (projected>remaining) return
            if (projected+degree+1_int64+merge(steps,0,axis==0.or.axis==direction)>huge(count_))return
            remaining=remaining/int(projected)
        end do
        fits=.true.
    end function refinement_fits_indices

    logical function span_subdivisions(patch,axis,divisions,insertions) result(valid)
        type(patch_entry), intent(in) :: patch
        integer, intent(in) :: axis,divisions
        real(rk), allocatable, intent(out) :: insertions(:)
        integer :: first,last,span,part,nspan,n,maximum,status
        real(rk) :: left,right,fraction,previous

        valid=.false.
        ! Propose evenly spaced values only inside the stored active spans.
        first=patch%degree(axis)+1
        last=patch%controls(axis)+1
        associate(knots=>patch%axes(axis)%knots)
            nspan=count(knots(first:last-1)<knots(first+1:last))
            if (nspan==0) return
            maximum=huge(n)-size(knots)
            if (divisions-1>maximum/nspan) return
            allocate(insertions(nspan*(divisions-1)),stat=status)
            if(status/=0)return
            n=0
            do span=first,last-1
                left=knots(span)
                right=knots(span+1)
                if (left==right) cycle
                previous=left
                do part=1,divisions-1
                    fraction=real(part,rk)/real(divisions,rk)
                    n=n+1
                    insertions(n)=(1._rk-fraction)*left+fraction*right
                    if (.not.ieee_is_finite(insertions(n))) return
                    if (insertions(n)<=previous.or.insertions(n)>=right) return
                    previous=insertions(n)
                end do
            end do
        end associate
        valid=.true.
    end function span_subdivisions

    module procedure commit_single_patch_interval
        type(patch_entry) :: candidate
        real(rk) :: lo,hi

        changed=.false.
        error='Choose a parameter direction of the current patch.'
        if (kind<1.or.kind>3.or.axis<1.or.axis>kind) return
        error='Enter finite bounds with the lower value below the upper value.'
        if (.not.all(ieee_is_finite([lower,upper])).or.lower>=upper) return
        lo=real(lower,rk)
        hi=real(upper,rk)
        if (.not.all(ieee_is_finite([lo,hi])).or.lo>=hi) return
        if (.not.ieee_is_finite(hi-lo)) return
        error='ForCAD could not load the patch.'
        if (.not.resolve_patch_snapshot(kind,candidate,patch_id,catalog)) return
        error=''
        if (all(candidate%bounds(:,axis)==[lo,hi])) return

        call reparameterize_patch(candidate%geometry,axis,lo,hi)
        error='ForCAD could not represent this parameter interval.'
        if (.not.refresh_patch_metadata(candidate)) return
        if (any(candidate%bounds(:,axis)/=[lo,hi])) return
        error='ForCAD could not sample the reparameterized patch.'
        if (.not.publish_edit(candidate,'parameter interval edit',axis,catalog,error)) return
        changed=.true.
        error=''
    end procedure commit_single_patch_interval

    module procedure commit_single_patch_knot_values
        type(patch_entry) :: candidate
        real(rk), allocatable :: values(:)
        integer :: n

        changed=.false.
        error='Choose a parameter direction of the current patch.'
        if (kind<1.or.kind>3.or.axis<1.or.axis>kind) return
        error='ForCAD could not load the patch.'
        if (.not.resolve_patch_snapshot(kind,candidate,patch_id,catalog)) return
        error='Supply the complete knot vector with finite values.'
        n=size(candidate%axes(axis)%knots)
        if (size(knots)/=n.or..not.all(ieee_is_finite(knots))) return
        values=real(knots,rk)
        if (.not.all(ieee_is_finite(values))) return

        associate(old=>candidate%axes(axis)%knots,bounds=>candidate%bounds(:,axis))
        error='Use the parameter interval editor to change the bounds.'
        if (any((old<=bounds(1).or.old>=bounds(2)).and.values/=old)) return
        error='Keep the knot groups in order and preserve their multiplicities.'
        if (any(values(2:)<values(:n-1))) return
        if (any((old(2:)==old(:n-1)).neqv.(values(2:)==values(:n-1)))) return
        error=''
        if (all(values==old)) return
        end associate

        candidate%axes(axis)%knots=values
        error='ForCAD rejected the edited knot vector.'
        if (.not.replace_native_knot_values(candidate%geometry,candidate%axes,candidate%degree)) return
        if (.not.refresh_patch_metadata(candidate)) return
        if (any(candidate%axes(axis)%knots/=values)) return
        error='ForCAD could not sample the edited patch.'
        if (.not.publish_edit(candidate,'knot value edit',axis,catalog,error)) return
        changed=.true.
        error=''
    end procedure commit_single_patch_knot_values

    module procedure commit_single_patch_split


        call publish_parameter_parts(kind,axis,value,value,created,error,split=.true.,catalog=catalog,patch_id=patch_id)
    end procedure commit_single_patch_split

    module procedure commit_single_patch_extract


        call publish_parameter_parts(kind,axis,lower,upper,created,error,split=.false.,catalog=catalog,patch_id=patch_id)
    end procedure commit_single_patch_extract

    subroutine publish_parameter_parts(kind,axis,lower,upper,created,error,split,catalog,patch_id)
        integer,intent(in) :: kind,axis
        real(c_double),intent(in) :: lower,upper
        integer,intent(out) :: created
        character(:),allocatable,intent(out) :: error
        logical,intent(in) :: split
        type(studio_patch_catalog), intent(inout), target, optional :: catalog
        integer, intent(in), optional :: patch_id
        type(studio_patch_catalog), pointer :: owner
        type(patch_entry) :: source,left,right
        type(patch_entry),allocatable :: resized(:)
        type(patch_edit) :: edit
        real(rk) :: lo,hi,start,finish
        integer :: count,variant
        character(*),parameter :: directions='UVW'
        character(:),allocatable :: suffix

        created=0
        error='Choose a parameter direction of the current patch.'
        if(kind<1.or.kind>3.or.axis<1.or.axis>kind)return
        error='ForCAD could not load the patch.'
        owner=>edit_catalog(catalog)
        if(.not.resolve_patch_snapshot(kind,source,patch_id,owner))return
        lo=real(lower,rk); hi=real(upper,rk)
        start=source%bounds(1,axis); finish=source%bounds(2,axis)
        if(split)then
            error='Enter a finite split parameter strictly inside the active interval.'
            if(.not.ieee_is_finite(lo).or.lo<=start.or.lo>=finish)return
            suffix=' · Upper '//directions(axis:axis)
            edit%operation='split'
            error='ForCAD could not split the patch at this parameter.'
        else
            error='Enter finite increasing bounds inside the active interval.'
            if(.not.all(ieee_is_finite([lo,hi])).or.lo<start.or.hi>finish.or.lo>=hi)return
            suffix=' · '//directions(axis:axis)//' extract'
            edit%operation='parameter extraction'
            error='ForCAD could not extract geometry in this range.'
        end if
        call parameter_parts(source%geometry,axis,lo,hi,split,left%geometry,right%geometry)
        right%dimension=kind
        right%resolution=source%resolution
        if(.not.refresh_patch_metadata(right))return
        if(.not.can_sample_patch(right))return
        if(split)then
            left%dimension=kind
            left%patch_id=source%patch_id
            left%name=source%name
            left%resolution=source%resolution
            if(.not.refresh_patch_metadata(left))return
            if(.not.can_sample_patch(left))return
        end if

        ! Stage sampled results before changing the catalog. Both split outputs
        ! share one history record, and the lower piece keeps the source identity.
        count=owner%patch_count(kind)
        variant=count+1
        allocate(resized(variant))
        if(allocated(owner%catalogs(kind)%entries))resized(:count)=owner%catalogs(kind)%entries
        right%patch_id=variant
        right%name=available_patch_name(source%name//suffix,resized(:0),owner)
        call advance_patch_revision(right)
        resized(variant)=right
        if(split)then
            allocate(edit%previous_set(1))
            edit%previous_set(1)=source
            call advance_patch_revision(left)
            resized(source%patch_id)=left
        end if
        edit%kind=kind
        edit%axis=axis
        edit%created_patch_id=variant
        edit%previous_selection=owner%selection(kind)
        call record_undo(edit,owner,error)
        if(len(error)>0)return
        call move_alloc(resized,owner%catalogs(kind)%entries)
        owner%selection(kind)=variant
        created=variant
        error=''
    end subroutine publish_parameter_parts

    module procedure commit_single_patch_reverse
        type(patch_entry) :: candidate

        changed=.false.
        error='Choose a parameter direction of the current patch.'
        if(kind<1.or.kind>3.or.axis<1.or.axis>kind)return
        error='ForCAD could not load the patch.'
        if(.not.resolve_patch_snapshot(kind,candidate,patch_id,catalog))return
        call reverse_patch(candidate%geometry,axis)
        error='ForCAD could not reverse the patch in this direction.'
        if(.not.refresh_patch_metadata(candidate))return
        changed=publish_edit(candidate,'reverse direction',axis,catalog,error)
        if(changed)error=''
    end procedure commit_single_patch_reverse

    logical function publish_edit(candidate,operation,axis,catalog,error) result(published)
        type(patch_entry), intent(in) :: candidate
        character(*), intent(in) :: operation
        integer, intent(in) :: axis
        type(studio_patch_catalog), intent(inout), target, optional :: catalog
        character(:), allocatable, intent(inout) :: error
        type(studio_patch_catalog), pointer :: owner
        type(patch_entry) :: previous
        integer :: kind,variant

        published=can_sample_patch(candidate)
        if (.not.published) return

        kind=candidate%dimension
        variant=candidate%patch_id
        owner=>edit_catalog(catalog)
        published=read_catalog_patch(kind,variant,previous,owner)
        if (.not.published) return
        call store_patch_template(previous,owner)
        call record_undo(patch_edit(kind=kind,axis=axis,operation=operation,previous=previous),owner,error)
        published=len(error)==0
        if(.not.published)return
        owner%catalogs(kind)%entries(variant)=candidate
        call advance_patch_revision(owner%catalogs(kind)%entries(variant))
    end function publish_edit

end submodule studio_single_patch_edits
