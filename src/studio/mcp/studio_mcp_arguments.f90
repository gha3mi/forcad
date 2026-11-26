module studio_mcp_arguments
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_c_binding, only: dp=>c_double
    use qt, only: QJsonObject, QJsonArray, QJsonValue_Array, QJsonValue_Bool, QJsonValue_Double, QJsonValue_Object, qt_string
    implicit none(type, external)
    private

    public :: read_geometry_identity
    public :: only_fields, bounded_integer, number_arg, vector_arg, real_array_arg, integer_array_arg, optional_boolean, axis_number

contains

    logical function only_fields(args,allowed) result(valid)
        type(QJsonObject),intent(in) :: args
        character(*),intent(in) :: allowed
        type(qt_string),allocatable :: keys(:)
        integer :: i

        keys=args%keys()
        valid=.false.
        do i=1,size(keys)
            if(len(keys(i)%s)==0.or.index(keys(i)%s,'|')>0)return
            if(index(allowed,'|'//keys(i)%s//'|')==0)return
        end do
        valid=.true.
    end function only_fields

    integer function bounded_integer(args,key,low,high,fallback) result(value)
        type(QJsonObject),intent(in) :: args
        character(*),intent(in) :: key
        integer,intent(in) :: low,high
        integer,intent(in),optional :: fallback

        value=-1
        if(present(fallback).and..not.args%contains(key))then
            value=fallback
            return
        end if
        if(args%valueType(key)/=QJsonValue_Double)return
        value=args%intValue(key,fallback=-1)
        if(value<low.or.value>high)value=-1
    end function bounded_integer

    logical function number_arg(args,key,value) result(valid)
        type(QJsonObject),intent(in) :: args
        character(*),intent(in) :: key
        real(dp),intent(out) :: value

        value=0._dp
        valid=args%valueType(key)==QJsonValue_Double
        if(.not.valid)return
        value=args%doubleValue(key)
        valid=ieee_is_finite(value)
    end function number_arg

    logical function vector_arg(args,key,vector) result(valid)
        type(QJsonObject),intent(in) :: args
        character(*),intent(in) :: key
        real(dp),intent(out) :: vector(3)
        real(dp), allocatable :: values(:)

        vector=0._dp
        valid=real_array_arg(args,key,values,3)
        if(valid)vector=values
    end function vector_arg

    logical function real_array_arg(args,key,values,required_size) result(valid)
        type(QJsonObject), intent(in) :: args
        character(*), intent(in) :: key
        real(dp), allocatable, intent(out) :: values(:)
        integer, intent(in), optional :: required_size
        type(QJsonArray) :: array
        integer :: i,n

        valid=.false.
        if(args%valueType(key)/=QJsonValue_Array)return
        array=args%arrayValue(key)
        n=array%size()
        if(present(required_size))then
            if(n/=required_size)then
                call array%delete()
                return
            end if
        end if
        valid=n>0.and.all([(array%typeAt(i-1)==QJsonValue_Double,i=1,n)])
        if(valid)then
            values=array%toDoubles()
            valid=all(ieee_is_finite(values))
        end if
        call array%delete()
    end function real_array_arg

    logical function integer_array_arg(args,key,values,min_size,max_size,high,minimum) result(valid)
        type(QJsonObject),intent(in) :: args
        character(*),intent(in) :: key
        integer,allocatable,intent(out) :: values(:)
        integer,intent(in) :: min_size,max_size,high
        integer,intent(in),optional :: minimum
        type(QJsonArray) :: array
        real(dp),allocatable :: numbers(:)
        integer :: n,i,low

        valid=.false.
        low=1
        if(present(minimum))low=minimum
        if(args%valueType(key)/=QJsonValue_Array)return
        array=args%arrayValue(key)
        n=array%size()
        if(n>=min_size.and.n<=max_size)then
            if(all([(array%typeAt(i-1)==QJsonValue_Double,i=1,n)]))then
                numbers=array%toDoubles()
                if(all(ieee_is_finite(numbers)))then
                    if(all(numbers>=real(low,dp).and.numbers<=real(high,dp)))then
                        values=int(numbers)
                        valid=all(real(values,dp)==numbers)
                    end if
                end if
            end if
        end if
        call array%delete()
    end function integer_array_arg

    logical function optional_boolean(args,key) result(valid)
        type(QJsonObject),intent(in) :: args
        character(*),intent(in) :: key
        valid=.not.args%contains(key).or.args%valueType(key)==QJsonValue_Bool
    end function optional_boolean

    integer function axis_number(name,axes) result(axis)
        character(*),intent(in) :: name,axes
        axis=0
        if(len(name)==1)axis=index(axes,name)
    end function axis_number

    logical function read_geometry_identity(args,kind,id) result(valid)
        type(QJsonObject), intent(in) :: args
        character(:), allocatable, intent(out) :: kind
        integer, intent(out) :: id
        type(QJsonObject) :: target

        valid=.false.
        id=0
        if(args%valueType('target')/=QJsonValue_Object)return
        target=args%objectValue('target')
        kind=target%stringValue('kind')
        id=bounded_integer(target,'id',1,huge(1))
        valid=only_fields(target,'|kind|id|').and.id>0.and.len(kind)>0
        call target%delete()
    end function read_geometry_identity

end module studio_mcp_arguments
