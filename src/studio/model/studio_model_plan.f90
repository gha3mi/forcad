! Validate and stage ordered native modeling operations without changing the open project.
module studio_model_plan
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use studio_forcad, only: rk, studio_native_patch, primitive_actions, primitive_names, stage_studio_primitive, &
        construction_actions, construction_names, stage_studio_construction
    use studio_presentation_text, only: integer_text, number_text
    implicit none(type, external)
    private
    public :: model_step, build_model_plan, model_plan_text, model_plan_max_steps

    integer, parameter :: model_plan_max_steps=12

    type :: model_step
        character(12) :: op=''
        character(64) :: name=''
        real(rk) :: origin(3)=0._rk, values(3)=0._rk, angle=0._rk
        integer :: source=0, spine=0, degree=0
        integer, allocatable :: sections(:)
    end type model_step

contains

    subroutine build_model_plan(steps,patches,names,error)
        type(model_step), intent(in) :: steps(:)
        type(studio_native_patch), allocatable, intent(out) :: patches(:)
        character(96), allocatable, intent(out) :: names(:)
        character(:), allocatable, intent(out) :: error
        type(studio_native_patch), allocatable :: objects(:)
        integer :: i, primitive, construction

        error='A modeling plan needs 1 to 12 operations.'
        if(size(steps)==0.or.size(steps)>model_plan_max_steps)return
        allocate(objects(size(steps)))
        do i=1,size(steps)
            associate(step=>steps(i))
                error='Every construction step needs a distinct name.'
                if(len_trim(step%name)==0.or.any(steps(:i-1)%name==step%name))return
                error='Use finite coordinates and dimensions.'
                if(.not.all(ieee_is_finite([step%origin,step%values,step%angle])))return
                primitive=findloc(primitive_actions,step%op,dim=1)
                construction=findloc(construction_actions,step%op,dim=1)
                if(primitive>0)then
                    call stage_studio_primitive(primitive,step%origin,step%values,objects(i),error)
                else if(construction>0)then
                    call stage_derived_step(step,construction,objects(:i-1),objects(i),error)
                else
                    error='This modeling operation is unavailable.'
                end if
                if(len(error)>0)then
                    error='Step '//integer_text(i)//' · '//trim(step%name)//new_line('a')//error
                    return
                end if
            end associate
        end do
        call move_alloc(objects,patches)
        allocate(names(size(steps)))
        names=steps%name
    end subroutine build_model_plan

    subroutine stage_derived_step(step,construction,earlier,result,error)
        type(model_step), intent(in) :: step
        integer, intent(in) :: construction
        type(studio_native_patch), intent(in) :: earlier(:)
        type(studio_native_patch), intent(out) :: result
        character(:), allocatable, intent(out) :: error
        integer :: i

        if(step%op=='loft')then
            error='Loft needs at least two distinct earlier step numbers as sections.'
            if(.not.allocated(step%sections))return
            if(size(step%sections)<2.or.size(step%sections)>size(earlier))return
            if(any(step%sections<1).or.any(step%sections>size(earlier)))return
            do i=2,size(step%sections)
                if(any(step%sections(:i-1)==step%sections(i)))return
            end do
            if(step%degree==0)then
                call stage_studio_construction(construction,earlier(step%sections),result,error)
            else
                call stage_studio_construction(construction,earlier(step%sections),result,error,degree=step%degree)
            end if
            return
        end if

        error='The profile must refer to an earlier step.'
        if(step%source<1.or.step%source>size(earlier))return
        if(step%op=='sweep')then
            error='The sweep path must refer to an earlier curve step distinct from the profile.'
            if(step%spine<1.or.step%spine>size(earlier).or.step%spine==step%source)return
            call stage_studio_construction(construction,earlier(step%source:step%source),result,error, &
                spine=earlier(step%spine))
        else
            call stage_studio_construction(construction,earlier(step%source:step%source),result,error, &
                vector=step%values,axis_point=step%origin,angle_degrees=step%angle)
        end if
    end subroutine stage_derived_step

    function model_plan_text(steps) result(description)
        type(model_step), intent(in) :: steps(:)
        character(:), allocatable :: description, details
        character, parameter :: lf=new_line('a')
        integer :: i, j, primitive, construction

        description='Create '//integer_text(size(steps))//' objects in model units.'//lf// &
            'Intermediate profiles stay editable. Undo removes the entire construction.'
        do i=1,size(steps)
            associate(step=>steps(i))
                primitive=findloc(primitive_actions,step%op,dim=1)
                construction=findloc(construction_actions,step%op,dim=1)
                description=description//lf//lf//integer_text(i)//'. '//trim(step%name)
                if(primitive>0)description=description//' · '//trim(primitive_names(primitive))
                if(construction>0)description=description//' · '//trim(construction_names(construction))
                select case(step%op)
                case('point')
                    details='Position '//tuple(step%origin)
                case('line')
                    details='Start '//tuple(step%origin)//' → End '//tuple(step%values)
                case('circle','cylinder')
                    details='Center '//tuple(step%origin)//' · Radius '//number_text(step%values(1))
                    if(step%op=='cylinder')details=details//' · Height '//number_text(step%values(2))
                case('arc')
                    details='Center '//tuple(step%origin)//' · Radius '//number_text(step%values(1))// &
                        ' · Start '//number_text(step%values(2))//'° · Sweep '//number_text(step%values(3))//'°'
                case('rectangle','box')
                    details='Origin '//tuple(step%origin)//' · Size '//tuple(step%values(:merge(2,3,step%op=='rectangle')))
                case('tube')
                    details='Base center '//tuple(step%origin)//' · Inner radius '//number_text(step%values(1))// &
                        ' · Outer radius '//number_text(step%values(2))//' · Height '//number_text(step%values(3))
                case('extrude')
                    details='Profile '//integer_text(step%source)//' · Displacement '//tuple(step%values)
                case('revolve')
                    details='Profile '//integer_text(step%source)//' · Axis point '//tuple(step%origin)// &
                        ' · Direction '//tuple(step%values)//' · Sweep '//number_text(step%angle)//'°'
                case('sweep')
                    details='Profile '//integer_text(step%source)//' · Path '//integer_text(step%spine)
                case('loft')
                    details='Sections'
                    if(allocated(step%sections))then
                        do j=1,size(step%sections)
                            if(j>1)details=details//','
                            details=details//' '//integer_text(step%sections(j))
                        end do
                    end if
                    if(step%degree==0)then
                        details=details//' · Automatic degree'
                    else
                        details=details//' · Degree '//integer_text(step%degree)
                    end if
                case default
                    details='Unavailable operation'
                end select
                description=description//lf//details
            end associate
        end do
    end function model_plan_text

    function tuple(values) result(text)
        real(rk), intent(in) :: values(:)
        character(:), allocatable :: text
        integer :: i

        text='('
        do i=1,size(values)
            if(i>1)text=text//', '
            text=text//number_text(values(i))
        end do
        text=text//')'
    end function tuple
end module studio_model_plan
