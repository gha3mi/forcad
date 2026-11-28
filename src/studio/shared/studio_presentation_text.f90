module studio_presentation_text
    use, intrinsic :: iso_c_binding, only: c_double
    implicit none(type, external)
    private

    public :: entity_name, integer_text, number_text, without_nul
    character(*), parameter, public :: studio_product_name='ForCAD Studio'
    character(*), parameter, public :: studio_version = '&
        &FORCAD_RELEASE_VERSION'

contains

    function integer_text(value) result(text)
        integer,intent(in)::value
        character(len=:),allocatable::text
        character(len=16)::buffer
        write(buffer,'(i0)')value
        text=trim(buffer)
    end function integer_text

    function entity_name(kind) result(name)
        integer,intent(in)::kind
        character(len=:),allocatable::name
        select case(kind)
        case(1); name='Curve'
        case(3); name='Volume'
        case default; name='Surface'
        end select
    end function entity_name

    function number_text(value) result(text)
        real(c_double), intent(in) :: value
        character(:), allocatable :: text
        character(len=40) :: buffer
        real(c_double) :: decoded
        integer :: status,exponent,last

        if (value==0._c_double) then
            text='0'
            return
        end if
        if (abs(value)>=1.e98_c_double.or.abs(value)<1.e-96_c_double) then
            write(buffer,'(es30.16e3)') value
        else
            write(buffer,'(g0.12)') value
            read(buffer,*,iostat=status,decimal='point') decoded
            if (status/=0) decoded=0._c_double
            if (decoded/=value) write(buffer,'(g0.17)') value
        end if
        text=trim(adjustl(buffer))
        exponent=index(text,'E')
        last=len(text)
        if (exponent>0) last=exponent-1
        if (index(text,'.')==0) return
        do while (last>1)
            if (text(last:last)/='0') exit
            last=last-1
        end do
        if (text(last:last)=='.') last=last-1
        if (exponent>0) then
            text=text(:last)//text(exponent:)
        else
            text=text(:last)
        end if
    end function number_text

    function without_nul(input) result(safe)
        character(*), intent(in) :: input
        character(len=len(input)) :: safe
        integer :: i

        safe=input
        do i=1,len(safe)
            if(safe(i:i)==achar(0))safe(i:i)=' '
        end do
    end function without_nul

end module studio_presentation_text
