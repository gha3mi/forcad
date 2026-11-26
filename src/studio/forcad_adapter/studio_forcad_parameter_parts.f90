module studio_forcad_parameter_parts
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use forcad, only: rk, nurbs_curve, nurbs_surface, nurbs_volume
    implicit none(type, external)
    private

    public :: parameter_parts, extract_native_slice

contains

    subroutine parameter_parts(geometry,axis,lo,hi,split,left,right)
        class(*), intent(in) :: geometry
        integer,intent(in) :: axis
        real(rk),intent(in) :: lo,hi
        logical,intent(in) :: split
        class(*), allocatable, intent(out) :: left,right
        type(nurbs_curve) :: lower_curve,upper_curve
        type(nurbs_surface) :: lower_surface,upper_surface
        type(nurbs_volume) :: lower_volume,upper_volume

        select type(geometry)
        type is(nurbs_curve)
            if(split)then
                call geometry%split(lo,lower_curve,upper_curve)
                allocate(left,source=lower_curve)
            else
                upper_curve=geometry%restrict(lo,hi)
            end if
            allocate(right,source=upper_curve)
        type is(nurbs_surface)
            if(split)then
                call geometry%split(axis,lo,lower_surface,upper_surface)
                allocate(left,source=lower_surface)
            else
                upper_surface=geometry%restrict(axis,lo,hi)
            end if
            allocate(right,source=upper_surface)
        type is(nurbs_volume)
            if(split)then
                call geometry%split(axis,lo,lower_volume,upper_volume)
                allocate(left,source=lower_volume)
            else
                upper_volume=geometry%restrict(axis,lo,hi)
            end if
            allocate(right,source=upper_volume)
        end select
    end subroutine parameter_parts

    subroutine extract_native_slice(geometry,axis,parameter,piece,error)
        class(*), intent(in) :: geometry
        integer, intent(in) :: axis
        real(rk), intent(in) :: parameter
        class(*), allocatable, intent(out) :: piece
        character(:), allocatable, intent(out) :: error
        type(nurbs_curve) :: curve
        type(nurbs_surface) :: surface,bounded_surface
        type(nurbs_volume) :: bounded_volume
        real(rk) :: lower,upper

        error='Enter a finite section parameter.'
        if(.not.ieee_is_finite(parameter))return
        error='Select a valid surface or volume and one of its parameter directions.'
        select type(geometry)
        type is(nurbs_surface)
            if(axis<1.or.axis>2.or..not.geometry%is_valid())return
            lower=geometry%get_knot(axis,geometry%get_degree(axis)+1)
            upper=geometry%get_knot(axis,geometry%get_nc(axis)+1)
            error='The section parameter must lie within the active patch domain.'
            if(parameter<lower.or.parameter>upper)return
            error='ForCAD could not extract the surface section.'
            if(parameter==upper.and.geometry%get_parameter_wrapping(axis))then
                bounded_surface=geometry%restrict(axis,lower,upper)
                if(.not.bounded_surface%err%ok)return
                curve=bounded_surface%slice(axis,parameter)
            else
                curve=geometry%slice(axis,parameter)
            end if
            if(.not.curve%err%ok.or..not.curve%is_valid())return
            allocate(piece,source=curve)
        type is(nurbs_volume)
            if(axis<1.or.axis>3.or..not.geometry%is_valid())return
            lower=geometry%get_knot(axis,geometry%get_degree(axis)+1)
            upper=geometry%get_knot(axis,geometry%get_nc(axis)+1)
            error='The section parameter must lie within the active patch domain.'
            if(parameter<lower.or.parameter>upper)return
            error='ForCAD could not extract the volume section.'
            if(parameter==upper.and.geometry%get_parameter_wrapping(axis))then
                bounded_volume=geometry%restrict(axis,lower,upper)
                if(.not.bounded_volume%err%ok)return
                surface=bounded_volume%slice(axis,parameter)
            else
                surface=geometry%slice(axis,parameter)
            end if
            if(.not.surface%err%ok.or..not.surface%is_valid())return
            allocate(piece,source=surface)
        class default
            return
        end select
        error=''
    end subroutine extract_native_slice
end module studio_forcad_parameter_parts
