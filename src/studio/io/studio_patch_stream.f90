module studio_patch_stream
    use, intrinsic :: iso_fortran_env, only: int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use studio_forcad, only: native_patch_data,collect_native_patch_data,restore_native_patch_data
    use studio_project_stream, only: remaining_file_units,real_units
    implicit none(type, external)
    private
    public :: write_patch,read_patch,indexable_shape

contains

    subroutine write_patch(unit,kind,geometry,status)
        integer, intent(in) :: unit,kind
        class(*), intent(in) :: geometry
        integer, intent(out) :: status
        type(native_patch_data) :: data
        integer :: axis,io

        status=1
        if(.not.collect_native_patch_data(geometry,data))return
        if(size(data%points,2)/=3)return
        if(.not.valid_patch_shape(kind,data%counts,data%degree,merge(1,0,data%wrap)))return
        write(unit,iostat=status)data%counts,data%degree,merge(1,0,data%wrap),size(data%weights)
        if(status/=0)return
        do axis=1,kind
            write(unit,iostat=status)data%axes(axis)%knots
            if(status/=0)return
        end do
        write(unit,iostat=io)data%points,data%weights
        status=merge(0,1,io==0)
    end subroutine write_patch

    subroutine read_patch(unit,kind,geometry,status)
        integer, intent(in) :: unit,kind
        class(*), allocatable, intent(out) :: geometry
        integer, intent(out) :: status
        type(native_patch_data) :: data
        integer :: wrap(3),weight_count,control_count,axis
        integer(int64) :: values

        read(unit,iostat=status)data%counts,data%degree,wrap,weight_count
        if(status/=0)return
        status=1
        if(.not.valid_patch_shape(kind,data%counts,data%degree,wrap))return
        control_count=product(data%counts(:kind))
        if(weight_count/=0.and.weight_count/=control_count)return
        values=3_int64*control_count+weight_count+sum(int(data%counts(:kind),int64)+data%degree(:kind)+1)
        if(values>remaining_file_units(unit)/real_units)return
        do axis=1,kind
            allocate(data%axes(axis)%knots(data%counts(axis)+data%degree(axis)+1),stat=status)
            if(status/=0)return
            read(unit,iostat=status)data%axes(axis)%knots
            if(status/=0)return
        end do
        allocate(data%points(control_count,3),stat=status)
        if(status/=0)return
        if(weight_count>0)then
            allocate(data%weights(weight_count),stat=status)
            if(status/=0)return
        end if
        if(allocated(data%weights))then
            read(unit,iostat=status)data%points,data%weights
        else
            read(unit,iostat=status)data%points
        end if
        if(status/=0)return
        status=1
        if(.not.all(ieee_is_finite(data%points)))return
        do axis=1,kind
            if(.not.all(ieee_is_finite(data%axes(axis)%knots)))return
        end do
        if(allocated(data%weights))then
            if(.not.all(ieee_is_finite(data%weights)))return
        end if
        status=0
        data%wrap=wrap==1
        call restore_native_patch_data(kind,data,geometry)
    end subroutine read_patch

    logical function valid_patch_shape(kind,counts,degree,wrap) result(valid)
        integer, intent(in) :: kind,counts(3),degree(3),wrap(3)

        valid=.false.
        if(kind<1.or.kind>size(counts))return
        if(any(counts<1).or..not.indexable_shape(counts(:kind)))return
        if(any(degree<0).or.any(degree>=counts))return
        if(any(int(counts,int64)+degree+1>huge(1)))return
        if(any(wrap<0).or.any(wrap>1))return
        valid=.true.
    end function valid_patch_shape

    pure logical function indexable_shape(counts) result(fits)
        integer, intent(in) :: counts(:)
        integer :: axis, capacity

        fits=.false.
        capacity=huge(capacity)/3
        do axis=1,size(counts)
            if(counts(axis)<1.or.counts(axis)>capacity)return
            capacity=capacity/counts(axis)
        end do
        fits=.true.
    end function indexable_shape

end module studio_patch_stream
