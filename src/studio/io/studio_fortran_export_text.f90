module studio_fortran_export_text
    use studio_forcad, only: rk
    use studio_presentation_text, only: integer_text
    implicit none(type,external)
    private
    public :: source_writer, integer_list, real_literal, logical_literal

    type :: source_writer
        integer :: unit=0, status=0
    contains
        procedure :: line => write_line
        procedure :: integers => write_integers
        procedure :: reals => write_reals
    end type source_writer

contains

    subroutine write_line(source,text)
        class(source_writer),intent(inout) :: source
        character(*),intent(in) :: text
        integer :: io

        if(source%status/=0)return
        write(source%unit,'(a)',iostat=io)text
        if(io/=0)source%status=io
    end subroutine write_line

    subroutine write_integers(source,target,values)
        class(source_writer),intent(inout) :: source
        character(*),intent(in) :: target
        integer,intent(in) :: values(:)
        integer,parameter :: values_per_line=8, values_per_assignment=384
        integer :: first,last,row,last_column

        do first=1,size(values),values_per_assignment
            if(source%status/=0)return
            last=min(size(values),first+values_per_assignment-1)
            call start_assignment(source,target,first,last)
            do row=first,last,values_per_line
                last_column=min(last,row+values_per_line-1)
                call finish_source_row(source,integer_list(values(row:last_column)),last_column==last)
            end do
        end do
    end subroutine write_integers

    subroutine write_reals(source,target,values)
        class(source_writer),intent(inout) :: source
        character(*),intent(in) :: target
        real(rk),intent(in) :: values(:)
        integer,parameter :: values_per_line=2, values_per_assignment=384
        character(:),allocatable :: row_text
        integer :: first,last,row,last_column,i

        do first=1,size(values),values_per_assignment
            if(source%status/=0)return
            last=min(size(values),first+values_per_assignment-1)
            call start_assignment(source,target,first,last)
            do row=first,last,values_per_line
                last_column=min(last,row+values_per_line-1)
                row_text=real_literal(values(row))
                do i=row+1,last_column
                    row_text=row_text//', '//real_literal(values(i))
                end do
                call finish_source_row(source,row_text,last_column==last)
            end do
        end do
    end subroutine write_reals

    subroutine start_assignment(source,target,first,last)
        type(source_writer),intent(inout) :: source
        character(*),intent(in) :: target
        integer,intent(in) :: first,last

        call source%line('        '//target//'('//integer_text(first)//':'//integer_text(last)//') = [ &')
    end subroutine start_assignment

    subroutine finish_source_row(source,text,last_row)
        type(source_writer),intent(inout) :: source
        character(*),intent(in) :: text
        logical,intent(in) :: last_row

        if(last_row)then
            call source%line('            '//text//' ]')
        else
            call source%line('            '//text//', &')
        end if
    end subroutine finish_source_row

    function integer_list(values) result(text)
        integer,intent(in) :: values(:)
        character(:),allocatable :: text
        integer :: i

        text=''
        do i=1,size(values)
            if(i>1)text=text//', '
            text=text//integer_text(values(i))
        end do
    end function integer_list

    function real_literal(value) result(text)
        real(rk),intent(in) :: value
        character(:),allocatable :: text
        integer,parameter :: decimal_places=precision(0._rk)+2
        character(24) :: format
        character(decimal_places+11) :: buffer

        if(value==0._rk)then
            text='0.0_rk'
            return
        end if
        if(value==1._rk)then
            text='1.0_rk'
            return
        end if
        write(format,'("(es",i0,".",i0,"e4)")')len(buffer),decimal_places
        write(buffer,format)value
        text=trim(adjustl(buffer))//'_rk'
    end function real_literal

    function logical_literal(value) result(text)
        logical,intent(in) :: value
        character(:),allocatable :: text

        text=trim(merge('.true. ','.false.',value))
    end function logical_literal

end module studio_fortran_export_text
