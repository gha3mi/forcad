module studio_forcad_patch_mutations
    use forcad, only: rk, nurbs_curve, nurbs_surface, nurbs_volume
    use studio_forcad_patch_data, only: native_knot_axis
    implicit none(type, external)
    private

    public :: replace_native_controls, replace_native_knot_values, edit_native_knots, insert_native_knots, edit_native_degree
    public :: reparameterize_patch, reverse_patch, translate_patch

contains

    logical function replace_native_controls(geometry,points,weights,points_changed,weights_changed) result(valid)
        class(*), intent(inout) :: geometry
        real(rk), intent(in), contiguous :: points(:,:),weights(:)
        logical, intent(in) :: points_changed,weights_changed
        valid=.false.
        select type (geometry)
        type is(nurbs_curve)
            if (points_changed) call geometry%modify_Xc(Xc=points)
            if (weights_changed) call geometry%modify_Wc(Wc=weights)
            valid=geometry%err%ok.and.geometry%is_valid()
        type is(nurbs_surface)
            if (points_changed) call geometry%modify_Xc(Xc=points)
            if (weights_changed) call geometry%modify_Wc(Wc=weights)
            valid=geometry%err%ok.and.geometry%is_valid()
        type is(nurbs_volume)
            if (points_changed) call geometry%modify_Xc(Xc=points)
            if (weights_changed) call geometry%modify_Wc(Wc=weights)
            valid=geometry%err%ok.and.geometry%is_valid()
        end select
    end function replace_native_controls

    logical function replace_native_knot_values(geometry,axes,degree) result(valid)
        class(*), intent(inout) :: geometry
        type(native_knot_axis), intent(in) :: axes(3)
        integer, intent(in) :: degree(3)
        real(rk), allocatable :: points(:,:),weights(:)
        logical :: wrapping(3)

        valid=.false.
        select type (geometry)
        type is(nurbs_curve)
            points=geometry%get_Xc()
            weights=geometry%get_Wc()
            if (size(weights)==0) deallocate(weights)
            wrapping(1)=geometry%get_parameter_wrapping()
            call geometry%set(axes(1)%knots,points,weights,degree(1),wrapping(1))
            valid=geometry%err%ok
        type is(nurbs_surface)
            points=geometry%get_Xc()
            weights=geometry%get_Wc()
            if (size(weights)==0) deallocate(weights)
            wrapping(:2)=geometry%get_parameter_wrapping([1,2])
            call geometry%set(axes(1)%knots,axes(2)%knots,points,weights,degree(:2),wrapping(:2))
            valid=geometry%err%ok
        type is(nurbs_volume)
            points=geometry%get_Xc()
            weights=geometry%get_Wc()
            if (size(weights)==0) deallocate(weights)
            wrapping=geometry%get_parameter_wrapping([1,2,3])
            call geometry%set(axes(1)%knots,axes(2)%knots,axes(3)%knots,points,weights, &
                degree,wrapping)
            valid=geometry%err%ok
        end select
    end function replace_native_knot_values

    logical function edit_native_knots(geometry,axis,value,repeats,remove,applied) result(valid)
        class(*), intent(inout) :: geometry
        integer, intent(in) :: axis,repeats
        real(rk), intent(in) :: value
        logical, intent(in) :: remove
        integer, intent(out) :: applied
        integer, allocatable :: removed(:)

        valid=.false.
        applied=0
        if (.not.remove) then
            valid=insert_native_knots(geometry,axis,[value],[repeats])
            if (valid) applied=repeats
            return
        end if
        select type (geometry)
        type is(nurbs_curve)
            call geometry%remove_knots(Xth=[value],r=[repeats],removed=removed)
            valid=geometry%err%ok
        type is(nurbs_surface)
            call geometry%remove_knots(dir=axis,Xth=[value],r=[repeats],removed=removed)
            valid=geometry%err%ok
        type is(nurbs_volume)
            call geometry%remove_knots(dir=axis,Xth=[value],r=[repeats],removed=removed)
            valid=geometry%err%ok
        end select
        if (.not.valid) return
        valid=allocated(removed)
        if (.not.valid) return
        valid=size(removed)==1
        if (valid) applied=removed(1)
    end function edit_native_knots

    logical function insert_native_knots(geometry,axis,values,repeats) result(valid)
        class(*), intent(inout) :: geometry
        integer, intent(in) :: axis
        real(rk), contiguous, intent(in) :: values(:)
        integer, contiguous, intent(in) :: repeats(:)

        valid=.false.
        select type (geometry)
        type is(nurbs_curve)
            call geometry%insert_knots(Xth=values,r=repeats)
            valid=geometry%err%ok
        type is(nurbs_surface)
            call geometry%insert_knots(dir=axis,Xth=values,r=repeats)
            valid=geometry%err%ok
        type is(nurbs_volume)
            call geometry%insert_knots(dir=axis,Xth=values,r=repeats)
            valid=geometry%err%ok
        end select
    end function insert_native_knots

    logical function edit_native_degree(geometry,axis,steps,reduce,tolerance,applied,bound) result(valid)
        class(*), intent(inout) :: geometry
        integer, intent(in) :: axis,steps
        logical, intent(in) :: reduce
        real(rk), intent(in) :: tolerance
        integer, intent(out) :: applied
        real(rk), intent(out) :: bound

        applied=steps
        bound=0._rk
        valid=.false.
        select type (geometry)
        type is(nurbs_curve)
            if (reduce) then
                call geometry%reduce_degree(t=steps,tol=tolerance,reduced=applied,error_bound=bound)
            else
                call geometry%elevate_degree(t=steps)
            end if
            valid=geometry%err%ok
        type is(nurbs_surface)
            if (reduce) then
                call geometry%reduce_degree(dir=axis,t=steps,tol=tolerance,reduced=applied,error_bound=bound)
            else
                call geometry%elevate_degree(dir=axis,t=steps)
            end if
            valid=geometry%err%ok
        type is(nurbs_volume)
            if (reduce) then
                call geometry%reduce_degree(dir=axis,t=steps,tol=tolerance,reduced=applied,error_bound=bound)
            else
                call geometry%elevate_degree(dir=axis,t=steps)
            end if
            valid=geometry%err%ok
        end select
    end function edit_native_degree

    subroutine reparameterize_patch(geometry,axis,lower,upper)
        class(*), intent(inout) :: geometry
        integer, intent(in) :: axis
        real(rk), intent(in) :: lower,upper

        select type (geometry)
        type is(nurbs_curve)
            geometry=geometry%reparameterize(a=lower,b=upper)
        type is(nurbs_surface)
            geometry=geometry%reparameterize(dir=axis,a=lower,b=upper)
        type is(nurbs_volume)
            geometry=geometry%reparameterize(dir=axis,a=lower,b=upper)
        end select
    end subroutine reparameterize_patch

    subroutine reverse_patch(geometry,axis)
        class(*), intent(inout) :: geometry
        integer, intent(in) :: axis

        select type(geometry)
        type is(nurbs_curve)
            geometry=geometry%reverse()
        type is(nurbs_surface)
            geometry=geometry%reverse(axis)
        type is(nurbs_volume)
            geometry=geometry%reverse(axis)
        end select
    end subroutine reverse_patch

    subroutine translate_patch(geometry,translation)
        class(*), intent(inout) :: geometry
        real(rk), intent(in) :: translation(3)

        select type(geometry)
        type is(nurbs_curve)
            call geometry%translate_Xc(translation)
        type is(nurbs_surface)
            call geometry%translate_Xc(translation)
        type is(nurbs_volume)
            call geometry%translate_Xc(translation)
        end select
    end subroutine translate_patch
end module studio_forcad_patch_mutations
