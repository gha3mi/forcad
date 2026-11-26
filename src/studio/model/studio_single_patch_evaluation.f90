submodule (studio_single_patch) studio_single_patch_evaluation
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use studio_forcad, only: evaluate_patch, evaluate_patch_basis, sample_patch_basis, check_patch_jacobian, &
        get_native_controls, sample_patch_geometry, sample_volume_slice
    implicit none(type, external)

contains

    module procedure read_single_patch
        type(patch_entry) :: patch
        real(rk) :: parameters(3)

        allocate(info%knots(0))
        info%name=''
        info%message='Choose a valid patch and knot direction.'
        if (kind<1.or.kind>3.or.axis<1.or.axis>kind) return
        info%message='Probe coordinates must be finite and between 0 and 1.'
        if (.not.valid_normalized(kind,normalized)) return
        info%message='ForCAD could not initialize the reference patch.'
        if (.not.resolve_patch_snapshot(kind,patch,patch_id,catalog)) then
            if (allocated(patch%geometry)) then
                if (studio_native_problem(patch%geometry)) then
                    info%diagnostic='ForCAD reported an initialization problem. Review the knots and control data.'
                end if
            end if
            return
        end if
        parameters=mapped_parameters(patch,normalized)
        call inspect_patch(patch,parameters,axis,info)
    end procedure read_single_patch

    module procedure read_single_patch_scene
        type(patch_entry) :: patch
        logical :: valid
        allocate(scene%points(3,0),scene%controls(3,0))
        if (kind<1.or.kind>3) return
        if (.not.resolve_patch_snapshot(kind,patch,variant,catalog)) then
            if (allocated(patch%geometry)) then
                if (studio_native_problem(patch%geometry)) then
                    scene%diagnostic='ForCAD reported an initialization problem. Review the knots and control data.'
                end if
            end if
            return
        end if
        valid=sample_scene(patch,scene)
        scene%valid=valid
    end procedure read_single_patch_scene

    module procedure read_single_patch_plot
        type(patch_entry) :: patch
        real(rk) :: parameters(3),lo,hi,samples(plot_samples)
        real(rk), allocatable :: samples_u(:),samples_v(:),samples_w(:),values(:,:)
        integer :: i
        allocate(plot%parameters(0),plot%values(0,0))
        if (kind<1.or.kind>3.or.axis<1.or.axis>kind) return
        if (.not.valid_normalized(kind,normalized)) return
        if (.not.resolve_patch_snapshot(kind,patch,patch_id,catalog)) return
        lo=patch%bounds(1,axis)
        hi=patch%bounds(2,axis)
        samples=[(lo+(hi-lo)*(real(i-1,rk)/real(plot_samples-1,rk)),i=1,plot_samples)]
        samples(plot_samples)=hi
        parameters=mapped_parameters(patch,normalized)
        samples_u=[parameters(1)]
        samples_v=[parameters(2)]
        samples_w=[parameters(3)]
        select case(axis)
        case(1)
            samples_u=samples
        case(2)
            samples_v=samples
        case(3)
            samples_w=samples
        case default
            return
        end select
        if (.not.sample_patch_basis(patch%geometry,samples_u,samples_v,samples_w,axis,derivative,values)) return
        if (.not.valid_matrix(values,plot_samples,product(patch%controls))) return
        plot%parameters=real(samples,c_double)
        plot%values=real(values,c_double)
        plot%valid=.true.
    end procedure read_single_patch_plot

    module procedure read_single_patch_slice
        type(patch_entry) :: patch
        real(rk) :: parameters(3)
        real(rk), allocatable :: native_points(:,:)
        allocate(points(3,0))
        if (.not.scene%valid.or.scene%dimension/=3) return
        if (scene%variant/=active_patch_id(3,catalog)) return
        if (.not.valid_normalized(3,normalized)) return
        if (any(scene%sample_counts(:2)<2)) return
        if (.not.resolve_patch_snapshot(3,patch,scene%variant,catalog)) return
        parameters=mapped_parameters(patch,normalized)
        if (.not.sample_volume_slice(patch%geometry,parameters(3),scene%sample_counts(:2),native_points)) return
        if (.not.valid_matrix(native_points,product(scene%sample_counts(:2)),3)) return
        points=transpose(real(native_points,c_double))
    end procedure read_single_patch_slice

    module procedure read_single_patch_controls
        type(patch_entry) :: patch
        real(rk), allocatable :: points(:,:),weights(:)
        controls%name=''
        allocate(controls%points(3,0),controls%weights(0))
        if (kind<1.or.kind>3) return
        if (.not.resolve_patch_snapshot(kind,patch,variant,catalog)) return
        if (.not.get_native_controls(patch%geometry,product(patch%controls),points,weights)) return
        controls%points=transpose(real(points,c_double))
        controls%weights=real(weights,c_double)
        if (.not.all(ieee_is_finite(controls%points)).or..not.all(ieee_is_finite(controls%weights))) return
        controls%dimension=kind
        controls%counts=patch%controls
        controls%name=patch%name
        controls%valid=.true.
    end procedure read_single_patch_controls

    logical function valid_normalized(kind,normalized) result(valid)
        integer, intent(in) :: kind
        real(c_double), intent(in) :: normalized(3)
        valid=.false.
        if (.not.all(ieee_is_finite(normalized(:kind)))) return
        valid=all(normalized(:kind)>=0.0_c_double.and.normalized(:kind)<=1.0_c_double)
    end function valid_normalized

    module procedure can_sample_patch
        type(patch_entry) :: validation

        ! Keep saved geometry and undo snapshots independent of display caches.
        validation=candidate
        valid=sample_scene(validation,scene)
        if(present(scene))scene%valid=valid
    end procedure can_sample_patch

    logical function sample_scene(patch,scene) result(valid)
        type(patch_entry), intent(inout) :: patch
        type(single_patch_scene), intent(inout), optional :: scene
        real(rk), allocatable :: native_points(:,:),native_controls(:,:)
        integer :: kind,direction,counts(3)

        valid=.false.
        kind=patch%dimension
        counts=1
        counts(:kind)=tensor_samples
        if(kind==1)counts(1)=curve_samples
        where(patch%resolution(:kind)>=2)counts(:kind)=patch%resolution(:kind)
        if(present(scene))then
            scene%sample_counts=counts
            scene%diagnostic='ForCAD reported a geometry sampling problem. Review the source patch and its resolution.'
        end if
        if(.not.sample_patch_geometry(patch%geometry,counts,native_points,native_controls))return
        if(present(scene))scene%diagnostic='ForCAD sampling returned unusable viewport data. Check the controls, weights and resolution.'
        if(.not.valid_matrix(native_points,product(counts),3))return
        if(.not.valid_matrix(native_controls,product(patch%controls),3))return
        if(rk/=c_double)then
            if(.not.all(ieee_is_finite(real(native_points,c_double))))return
            if(.not.all(ieee_is_finite(real(native_controls,c_double))))return
        end if
        valid=.true.
        if(.not.present(scene))return
        scene%points=transpose(real(native_points,c_double))
        scene%controls=transpose(real(native_controls,c_double))
        scene%dimension=kind
        scene%variant=patch%patch_id
        scene%degree=patch%degree
        scene%control_counts=patch%controls
        scene%bounds=real(patch%bounds,c_double)
        do direction=1,kind
            scene%axes(direction)%knots=real(patch%axes(direction)%knots,c_double)
        end do
        deallocate(scene%diagnostic)
    end function sample_scene

    function mapped_parameters(patch,normalized) result(parameters)
        type(patch_entry), intent(in) :: patch
        real(c_double), intent(in) :: normalized(3)
        real(rk) :: parameters(3)
        integer :: n
        parameters=0._rk
        n=patch%dimension
        parameters(:n)=patch%bounds(1,:n)+(patch%bounds(2,:n)-patch%bounds(1,:n))*real(normalized(:n),rk)
        where (normalized(:n)>=1.0_c_double) parameters(:n)=patch%bounds(2,:n)
    end function mapped_parameters

    subroutine inspect_patch(patch,parameters,axis,info)
        type(patch_entry), intent(inout) :: patch
        real(rk), intent(in) :: parameters(3)
        integer, intent(in) :: axis
        type(single_patch_info), intent(inout) :: info
        real(rk), allocatable :: value(:), weights(:),knots(:)
        logical :: inverted,regular
        integer :: direction,order(3),strongest,control,span,i,p,nc
        info%name=patch%name
        info%dimension=patch%dimension
        info%degree=patch%degree
        info%controls(:patch%dimension)=patch%controls(:patch%dimension)
        info%bounds=real(patch%bounds,c_double)
        info%parameters=real(parameters,c_double)
        info%knots=real(patch%axes(axis)%knots,c_double)
        info%message='ForCAD could not evaluate the reference patch at this probe.'
        info%diagnostic='ForCAD returned invalid point or derivative data at this probe. Try a nearby parameter.'
        value=evaluate_patch(patch%geometry,parameters,[0,0,0])
        if (size(value)/=3.or..not.all(ieee_is_finite(value))) return
        info%point=real(value,c_double)
        do direction=1,patch%dimension
            order=0
            order(direction)=1
            value=evaluate_patch(patch%geometry,parameters,order)
            if (size(value)/=3.or..not.all(ieee_is_finite(value))) return
            info%jacobian(:,direction)=real(value,c_double)
        end do

        call evaluate_patch_basis(patch%geometry,parameters,weights)
        if (studio_native_problem(patch%geometry)) then
            info%diagnostic='ForCAD reported a basis evaluation problem at this probe. Try a nearby parameter.'
            return
        end if
        info%diagnostic='ForCAD returned unusable basis or span data. Review the knot vector at this probe.'
        if (.not.allocated(weights)) return
        if (size(weights)/=product(patch%controls).or..not.all(ieee_is_finite(weights))) return
        strongest=maxloc(weights,dim=1)-1
        do direction=1,patch%dimension
            nc=patch%controls(direction)
            p=patch%degree(direction)
            knots=patch%axes(direction)%knots
            span=0
            do i=p+1,nc
                if (knots(i+1)<=knots(i).or.parameters(direction)<knots(i)) cycle
                if (parameters(direction)>=knots(i+1).and.knots(i+1)<patch%bounds(2,direction)) cycle
                span=i
                exit
            end do
            if (span==0) return
            info%cell(:,direction)=real(knots(span:span+1),c_double)
            ! Convert the strongest native tensor-basis index to its axis index.
            control=mod(strongest,nc)+1
            strongest=strongest/nc
            info%support(:,direction)=real([max(knots(control),patch%bounds(1,direction)), &
                min(knots(control+p+1),patch%bounds(2,direction))],c_double)
        end do
        info%message=''
        info%valid=.true.
        deallocate(info%diagnostic)
        ! Use ForCAD's geometry frame check, scoped to this probe rather than the whole patch.
        call check_patch_jacobian(real(info%jacobian(:,:patch%dimension),rk),inverted,regular)
        if (.not.regular) then
            info%diagnostic='ForCAD probe warning: the local Jacobian is singular or cannot be represented. '// &
                'Inspect the controls or try a nearby parameter.'
        else if (inverted) then
            info%diagnostic='ForCAD probe warning: the volume orientation is reversed at this point. '// &
                'Check the intended orientation before analysis.'
        end if
    end subroutine inspect_patch

    logical function valid_matrix(values,rows,columns) result(valid)
        real(rk), intent(in) :: values(:,:)
        integer, intent(in) :: rows,columns
        valid=all(shape(values)==[rows,columns]).and.all(ieee_is_finite(values))
    end function valid_matrix

end submodule studio_single_patch_evaluation
