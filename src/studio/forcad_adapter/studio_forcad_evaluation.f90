module studio_forcad_evaluation
    use forcad, only: rk, nurbs_curve, nurbs_surface, nurbs_volume, &
        nurbs_multipatch_curve, nurbs_multipatch_surface, nurbs_multipatch_volume
    use forcad_utils, only: jacobian_frame
    implicit none(type, external)
    private

    public :: evaluate_patch, evaluate_patch_basis, sample_patch_basis, check_patch_jacobian, studio_native_problem

contains

    function evaluate_patch(geometry,parameters,order) result(value)
        class(*), intent(in) :: geometry
        real(rk), intent(in) :: parameters(3)
        integer, intent(in) :: order(3)
        real(rk), allocatable :: value(:)

        select type(geometry)
        type is(nurbs_curve)
            value=geometry%cmp_Xg(parameters(1),order=order(1))
        type is(nurbs_surface)
            value=geometry%cmp_Xg(parameters(:2),order=order(:2))
        type is(nurbs_volume)
            value=geometry%cmp_Xg(parameters,order=order)
        class default
            allocate(value(0))
        end select
    end function evaluate_patch

    subroutine evaluate_patch_basis(geometry,parameters,weights)
        class(*), intent(inout) :: geometry
        real(rk), intent(in) :: parameters(3)
        real(rk), allocatable, intent(out) :: weights(:)

        select type(geometry)
        type is(nurbs_curve)
            call geometry%basis(Xt=parameters(1),Tgc=weights)
        type is(nurbs_surface)
            call geometry%basis(Xt=parameters(:2),Tgc=weights)
        type is(nurbs_volume)
            call geometry%basis(Xt=parameters,Tgc=weights)
        class default
            return
        end select
    end subroutine evaluate_patch_basis

    logical function sample_patch_basis(geometry,samples_u,samples_v,samples_w,axis,derivative,values) result(valid)
        class(*), intent(inout) :: geometry
        real(rk), contiguous, intent(in) :: samples_u(:),samples_v(:),samples_w(:)
        integer, intent(in) :: axis
        logical, intent(in) :: derivative
        real(rk), allocatable, intent(out) :: values(:,:)
        real(rk), allocatable :: partials(:,:,:)
        integer :: dimension

        valid=.false.
        select type(geometry)
        type is(nurbs_curve)
            dimension=1
            if (derivative) then
                call geometry%derivative(Xt=samples_u,dTgc=values)
            else
                call geometry%basis(Xt=samples_u,Tgc=values)
            end if
        type is(nurbs_surface)
            dimension=2
            if (derivative) then
                call geometry%derivative(Xt1=samples_u,Xt2=samples_v,dTgc=partials)
            else
                call geometry%basis(Xt1=samples_u,Xt2=samples_v,Tgc=values)
            end if
        type is(nurbs_volume)
            dimension=3
            if (derivative) then
                call geometry%derivative(Xt1=samples_u,Xt2=samples_v,Xt3=samples_w,dTgc=partials)
            else
                call geometry%basis(Xt1=samples_u,Xt2=samples_v,Xt3=samples_w,Tgc=values)
            end if
        class default
            return
        end select
        if (studio_native_problem(geometry)) return
        if (derivative.and.dimension>1) then
            if (axis<1.or.axis>dimension.or..not.allocated(partials)) return
            if (size(partials,3)/=dimension) return
            values=partials(:,:,axis)
        end if
        valid=allocated(values)
    end function sample_patch_basis

    subroutine check_patch_jacobian(jacobian,inverted,regular)
        real(rk), intent(in) :: jacobian(:,:)
        logical, intent(out) :: inverted,regular
        real(rk) :: frame(3,size(jacobian,2)),measure

        call jacobian_frame(jacobian,frame,measure,inverted,regular)
    end subroutine check_patch_jacobian

    logical function studio_native_problem(geometry) result(reported)
        class(*), intent(in) :: geometry

        select type(geometry)
        type is(nurbs_curve)
            reported=.not.geometry%err%ok
        type is(nurbs_surface)
            reported=.not.geometry%err%ok
        type is(nurbs_volume)
            reported=.not.geometry%err%ok
        type is(nurbs_multipatch_curve)
            reported=.not.geometry%err%ok
        type is(nurbs_multipatch_surface)
            reported=.not.geometry%err%ok
        type is(nurbs_multipatch_volume)
            reported=.not.geometry%err%ok
        class default
            reported=.false.
        end select
    end function studio_native_problem

end module studio_forcad_evaluation
