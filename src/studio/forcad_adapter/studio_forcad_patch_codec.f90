module studio_forcad_patch_codec
    use forcad, only: rk, nurbs_curve, nurbs_surface, nurbs_volume
    use studio_forcad_patch_data, only: native_knot_axis
    implicit none(type, external)
    private


    type, public :: native_patch_data
        integer :: degree(3)=0, counts(3)=1
        logical :: wrap(3)=.false.
        type(native_knot_axis) :: axes(3)
        real(rk), allocatable :: points(:,:), weights(:)
    end type native_patch_data

    public :: collect_native_patch_data, restore_native_patch_data

contains

    logical function collect_native_patch_data(geometry,data) result(available)
        class(*), intent(in) :: geometry
        type(native_patch_data), intent(out) :: data
        integer :: axis

        available=.false.
        select type(geometry)
        type is(nurbs_curve)
            data%degree(1)=geometry%get_degree()
            data%counts(1)=geometry%get_nc()
            data%wrap(1)=geometry%get_parameter_wrapping()
            data%axes(1)%knots=geometry%get_knot()
            data%points=geometry%get_Xc()
            data%weights=geometry%get_Wc()
        type is(nurbs_surface)
            data%degree(:2)=geometry%get_degree()
            data%counts(:2)=geometry%get_nc()
            data%wrap(:2)=geometry%get_parameter_wrapping([1,2])
            do axis=1,2
                data%axes(axis)%knots=geometry%get_knot(axis)
            end do
            data%points=geometry%get_Xc()
            data%weights=geometry%get_Wc()
        type is(nurbs_volume)
            data%degree=geometry%get_degree()
            data%counts=geometry%get_nc()
            data%wrap=geometry%get_parameter_wrapping([1,2,3])
            do axis=1,3
                data%axes(axis)%knots=geometry%get_knot(axis)
            end do
            data%points=geometry%get_Xc()
            data%weights=geometry%get_Wc()
        class default
            return
        end select
        available=.true.
    end function collect_native_patch_data

    subroutine restore_native_patch_data(kind,data,geometry)
        integer, intent(in) :: kind
        type(native_patch_data), intent(in) :: data
        class(*), allocatable, intent(out) :: geometry

        select case(kind)
        case(1)
            allocate(nurbs_curve::geometry)
            select type(geometry)
            type is(nurbs_curve)
                call geometry%set(data%axes(1)%knots,data%points,data%weights,data%degree(1),data%wrap(1))
            end select
        case(2)
            allocate(nurbs_surface::geometry)
            select type(geometry)
            type is(nurbs_surface)
                call geometry%set(data%axes(1)%knots,data%axes(2)%knots,data%points,data%weights,data%degree(:2),data%wrap(:2))
            end select
        case(3)
            allocate(nurbs_volume::geometry)
            select type(geometry)
            type is(nurbs_volume)
                call geometry%set(data%axes(1)%knots,data%axes(2)%knots,data%axes(3)%knots,data%points,data%weights,data%degree,data%wrap)
            end select
        case default
            return
        end select
    end subroutine restore_native_patch_data

end module studio_forcad_patch_codec
