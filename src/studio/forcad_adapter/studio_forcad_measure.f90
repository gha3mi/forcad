module studio_forcad_measure
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use forcad, only: rk, nurbs_curve, nurbs_surface, nurbs_volume
    implicit none(type, external)
    private
    public :: measure_native_patch
contains
    subroutine measure_native_patch(geometry,measure,valid,orientation_ok)
        class(*), intent(inout) :: geometry
        real(rk), intent(out) :: measure
        logical, intent(out) :: valid
        logical, intent(out), optional :: orientation_ok

        if(present(orientation_ok))orientation_ok=.true.
        measure=0._rk
        valid=.false.
        select type(geometry)
        type is(nurbs_curve)
            call geometry%err%reset()
            call geometry%cmp_length(measure)
            valid=geometry%err%ok
        type is(nurbs_surface)
            call geometry%err%reset()
            call geometry%cmp_area(measure)
            valid=geometry%err%ok
        type is(nurbs_volume)
            call geometry%err%reset()
            if(present(orientation_ok))then
                call geometry%cmp_volume(measure,strict=.true.)
                orientation_ok=geometry%err%ok
                if(orientation_ok)then
                    valid=ieee_is_finite(measure).and.measure>=0._rk
                    return
                end if
                call geometry%err%reset()
            end if
            call geometry%cmp_volume(measure)
            valid=geometry%err%ok
        end select
        valid=valid.and.ieee_is_finite(measure).and.measure>=0._rk
    end subroutine measure_native_patch
end module studio_forcad_measure
