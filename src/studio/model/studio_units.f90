! Length units describe ForCAD Studio coordinates without changing the stored geometry.
module studio_units
    use, intrinsic :: iso_c_binding, only: c_double
    implicit none(type, external)
    private
    public :: unitless, unit_mm, unit_cm, unit_m, unit_in, unit_ft, unit_labels
    public :: valid_unit_id, unit_symbol, length_caption, display_length, model_length

    integer, parameter :: unitless=0, unit_mm=1, unit_cm=2, unit_m=3, unit_in=4, unit_ft=5
    character(20), parameter :: unit_labels(0:5)=[character(20) :: &
        'Unitless', 'Millimeters (mm)', 'Centimeters (cm)', 'Meters (m)', 'Inches (in)', 'Feet (ft)']
    character(2), parameter :: symbols(0:5)=[character(2) :: '', 'mm', 'cm', 'm', 'in', 'ft']
    real(c_double), parameter :: metres_per_unit(1:5)=[.001_c_double, .01_c_double, 1._c_double, .0254_c_double, .3048_c_double]

contains

    pure elemental logical function valid_unit_id(unit_id)
        integer, intent(in) :: unit_id
        valid_unit_id=unit_id>=unitless.and.unit_id<=unit_ft
    end function valid_unit_id

    pure function unit_symbol(unit_id) result(symbol)
        integer, intent(in), optional :: unit_id
        character(2) :: symbol

        symbol=''
        if(.not.present(unit_id))return
        if(valid_unit_id(unit_id))symbol=symbols(unit_id)
    end function unit_symbol

    pure function length_caption(caption,unit_id) result(text)
        character(*), intent(in) :: caption
        integer, intent(in), optional :: unit_id
        character(:), allocatable :: text
        character(2) :: symbol

        symbol=unit_symbol(unit_id)
        text=trim(caption)
        if(symbol/='')text=text//' ('//trim(symbol)//')'
    end function length_caption

    pure elemental real(c_double) function display_length(value,model_unit,display_unit) result(display_value)
        real(c_double), intent(in) :: value
        integer, intent(in), optional :: model_unit,display_unit
        display_value=value*unit_scale(model_unit,display_unit)
    end function display_length

    pure elemental real(c_double) function model_length(value,model_unit,display_unit) result(model_value)
        real(c_double), intent(in) :: value
        integer, intent(in), optional :: model_unit,display_unit
        model_value=value*unit_scale(display_unit,model_unit)
    end function model_length

    pure real(c_double) function unit_scale(from_unit,to_unit) result(scale)
        integer, intent(in), optional :: from_unit,to_unit

        scale=1._c_double
        if(.not.present(from_unit).or..not.present(to_unit))return
        if(.not.valid_unit_id(from_unit).or..not.valid_unit_id(to_unit))return
        if(from_unit==unitless.or.to_unit==unitless.or.from_unit==to_unit)return
        scale=metres_per_unit(from_unit)/metres_per_unit(to_unit)
    end function unit_scale

end module studio_units
