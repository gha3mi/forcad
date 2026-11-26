module studio_project_stream
    use, intrinsic :: iso_c_binding, only: c_null_char
    use, intrinsic :: iso_fortran_env, only: int64,file_storage_size
    use studio_forcad, only: rk
    implicit none(type, external)
    private
    public :: write_text,read_text,remaining_file_units,integer_units,real_units

    integer, parameter :: integer_units=storage_size(0)/file_storage_size, real_units=storage_size(0._rk)/file_storage_size
    integer, parameter :: character_units=storage_size('a')/file_storage_size

contains

    subroutine write_text(unit,value,status)
        integer, intent(in) :: unit
        character(*), intent(in) :: value
        integer, intent(out) :: status
        integer :: io

        status=1
        if(len_trim(value)==0.or.index(value,c_null_char)/=0)return
        write(unit,iostat=io)len(value),value
        status=merge(0,1,io==0)
    end subroutine write_text

    subroutine read_text(unit,value,status)
        integer, intent(in) :: unit
        character(:), allocatable, intent(out) :: value
        integer, intent(out) :: status
        integer :: length

        read(unit,iostat=status)length
        if(status/=0)return
        status=1
        if(length<1.or.length>remaining_file_units(unit)/character_units)return
        allocate(character(length)::value,stat=status)
        if(status/=0)return
        read(unit,iostat=status)value
        if(status/=0)return
        if(len_trim(value)==0.or.index(value,c_null_char)/=0)status=1
    end subroutine read_text
    integer(int64) function remaining_file_units(unit) result(remaining)
        integer, intent(in) :: unit
        integer(int64) :: pos,length
        integer :: status

        remaining=0
        inquire(unit=unit,pos=pos,size=length,iostat=status)
        if(status/=0)return
        if(pos<1.or.pos>length)return
        remaining=length-pos+1
    end function remaining_file_units

end module studio_project_stream
