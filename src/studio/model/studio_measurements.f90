module studio_measurements
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use studio_document, only: studio_document_t, document_has_patch
    use studio_forcad, only: rk, measure_native_patch, describe_patch, extract_native_slice
    use studio_single_patch, only: copy_single_patch_geometry, single_patch_name
    use studio_selection, only: selection_target_t, selection_set_t, selection_patch, selection_add, &
        pick_face, pick_edge, pick_vertex, patch_edge_boundaries
    use studio_units, only: display_length, unit_symbol, unitless
    use studio_presentation_text, only: integer_text
    implicit none(type, external)
    private
    public :: measure_document_geometry, measure_document_selection, measurement_label, measurement_unit, measurement_warning
    character(6), parameter :: labels(3)=[character(6) :: 'Length','Area','Volume']
    character(*), parameter :: missing_geometry='Select a curve, edge, surface, face, volume, or assembly to measure.'

    type, public :: measurement_t
        real(c_double) :: value=0._c_double
        integer :: count=0
        logical :: reversed=.false.
        character(:), allocatable :: error
    end type measurement_t

contains
    subroutine measure_document_geometry(document,kind,patch_ids,measure,error,reversed)
        type(studio_document_t), intent(in) :: document
        integer, intent(in) :: kind,patch_ids(:)
        real(c_double), intent(out) :: measure
        character(:), allocatable, intent(out) :: error
        logical, intent(out), optional :: reversed
        type(measurement_t) :: results(3)
        integer :: i

        measure=0._c_double
        if(present(reversed))reversed=.false.
        error=missing_geometry
        if(kind<1.or.kind>3.or.size(patch_ids)==0)return
        if(any(patch_ids<1))then
            error='A selected patch is no longer available.'
            return
        end if
        call measure_document_selection(document,[(selection_target_t(selection_patch,patch_ids(i),kind),i=1,size(patch_ids))],results)
        if(results(kind)%count==0)return
        measure=results(kind)%value
        error=results(kind)%error
        if(present(reversed))reversed=results(kind)%reversed
    end subroutine measure_document_geometry

    subroutine measure_document_selection(document,targets,results)
        type(studio_document_t), intent(in) :: document
        type(selection_target_t), intent(in) :: targets(:)
        type(measurement_t), intent(out) :: results(3)
        type(selection_set_t) :: measured
        type(selection_target_t) :: target
        real(rk) :: sums(3),part
        character(:), allocatable :: error
        integer :: i,kind,sides(2)
        logical :: orientation_ok

        sums=0._rk
        do i=1,size(targets)
            target=targets(i)
            if(target%category/=selection_patch.or.target%rank<1.or.target%rank>3.or.target%feature==pick_vertex)cycle
            kind=target%rank
            sides=0
            select case(target%feature)
            case(pick_face)
                if(kind==3.and.target%item>=1.and.target%item<=6) sides(1)=target%item
            case(pick_edge)
                sides=patch_edge_boundaries(kind,target%item)
            case default
                continue
            end select
            if(all(sides==0))then
                target%feature=0
                target%item=0
            end if
            if(measured%contains(target))cycle
            call measured%change(target,selection_add)
            kind=kind-count(sides>0)
            associate(result=>results(kind))
                if(result%count==0)result%error=''
                result%count=result%count+1
                if(len(result%error)>0)cycle
                call measure_patch(document,target,sides,part,orientation_ok,error)
                if(len(error)>0)then
                    result%error=error
                    cycle
                end if
                if(part>huge(part)-sums(kind))then
                    result%error='The total measure exceeds the numeric range.'
                    cycle
                end if
                sums(kind)=sums(kind)+part
                result%reversed=result%reversed.or..not.orientation_ok
            end associate
        end do

        do kind=1,3
            associate(result=>results(kind))
                if(result%count==0)cycle
                if(len(result%error)>0)cycle
                result%value=display_length(1._c_double,document%model_unit,document%display_unit)**kind*real(sums(kind),c_double)
                if(.not.ieee_is_finite(result%value))then
                    result%value=0._c_double
                    result%error='The measure in display units exceeds the numeric range.'
                end if
            end associate
        end do
    end subroutine measure_document_selection

    subroutine measure_patch(document,target,sides,value,orientation_ok,error)
        use studio_forcad, only: native_knot_axis
        type(studio_document_t), intent(in) :: document
        type(selection_target_t), intent(in) :: target
        integer, intent(in) :: sides(2)
        real(rk), intent(out) :: value
        logical, intent(out) :: orientation_ok
        character(:), allocatable, intent(out) :: error
        class(*), allocatable :: geometry,piece
        type(native_knot_axis) :: axes(3)
        integer :: degree(3),controls(3),axis,end_point,i,remaining(3)
        real(rk) :: bounds(2,3)
        logical :: valid

        orientation_ok=.true.
        error='A selected patch is no longer available.'
        if(.not.document_has_patch(document,target%rank,target%id))return
        call copy_single_patch_geometry(target%rank,target%id,geometry,catalog=document%patches)
        if(.not.allocated(geometry))return
        remaining=[1,2,3]
        do i=1,size(sides)
            if(sides(i)==0)exit
            axis=findloc(remaining,(sides(i)+1)/2,dim=1)
            end_point=2-mod(sides(i),2)
            error='ForCAD could not read the selected boundary.'
            if(.not.describe_patch(geometry,target%rank-i+1,degree,controls,axes,bounds))return
            call extract_native_slice(geometry,axis,bounds(end_point,axis),piece,error)
            if(len(error)>0)return
            call move_alloc(piece,geometry)
            if(axis<3)remaining(axis:2)=remaining(axis+1:3)
            remaining(3)=0
        end do
        call measure_native_patch(geometry,value,valid,orientation_ok)
        if(.not.valid)then
            error='ForCAD could not measure '//single_patch_name(target%rank,target%id,catalog=document%patches)//'.'
            return
        end if
        error=''
    end subroutine measure_patch

    function measurement_warning(result) result(text)
        type(measurement_t), intent(in) :: result
        character(:), allocatable :: text

        text=''
        if(result%count>1)text='Sum of '//integer_text(result%count)//' selected measures. Overlaps contribute to the sum.'
        if(result%reversed)then
            if(len(text)>0)text=text//new_line('a')
            text=text//'ForCAD detected reversed orientation at volume quadrature points. The displayed volume is unsigned.'
        end if
    end function measurement_warning

    pure function measurement_label(kind) result(label)
        integer, intent(in) :: kind
        character(:), allocatable :: label

        label='Measure'
        if(kind>=1.and.kind<=size(labels))label=trim(labels(kind))
    end function measurement_label

    pure function measurement_unit(document,kind) result(units)
        type(studio_document_t), intent(in) :: document
        integer, intent(in) :: kind
        character(:), allocatable :: units

        units=''
        if(document%model_unit==unitless.or.document%display_unit==unitless)return
        units=trim(unit_symbol(document%display_unit))
        if(kind==2)units=units//'²'
        if(kind==3)units=units//'³'
    end function measurement_unit
end module studio_measurements
