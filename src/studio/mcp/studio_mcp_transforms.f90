module studio_mcp_transforms
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: int64
    use qt, only: QJsonArray, QJsonObject, QJsonValue_Array, QJsonValue_Double, QJsonValue_Object
    use studio_forcad, only: affine_transform_t, rotation_matrix
    use studio_state, only: studio_t, dp
    use studio_scene, only: project_geometry_targets
    use studio_document, only: document_patch_in_assembly
    use studio_single_patch, only: scene_transform_step, scene_rotate, scene_translate, &
        active_patch_id, commit_scene_transforms
    use studio_edit_history, only: document_undo_stamp
    use studio_ui, only: apply_geometry_transform, refresh_rotated_scene
    use studio_mcp_arguments, only: only_fields, number_arg, vector_arg, axis_number
    use studio_mcp_inspection, only: describe_selection, active_rank
    use studio_presentation_text, only: integer_text
    use studio_units, only: unitless, unit_ft, unit_symbol, model_length
    implicit none(type, external)
    private

    public :: transform_scene, transform_selection

contains

    subroutine transform_scene(s,args,result,error)
        type(studio_t),intent(inout),target :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        type(QJsonArray) :: list
        type(QJsonObject) :: item
        type(scene_transform_step),allocatable :: steps(:)
        integer :: i
        logical :: valid,changed

        error='Supply at least one rotate or move step.'
        if(.not.only_fields(args,'|steps|').or.args%valueType('steps')/=QJsonValue_Array)return
        list=args%arrayValue('steps')
        if(list%size()<1)then
            call list%delete()
            return
        end if
        allocate(steps(list%size()))
        do i=1,size(steps)
            valid=list%typeAt(i-1)==QJsonValue_Object
            if(valid)then
                item=list%objectAt(i-1)
                valid=parse_transform(item,steps(i))
                call item%delete()
            end if
            if(valid)cycle
            call list%delete()
            error='Transform step '//integer_text(i)//' needs an X, Y or Z axis and a valid rotate or move value.'
            return
        end do
        call list%delete()
        call commit_scene_transforms(project_geometry_targets(s%document%show_templates,patches=s%document%patches), &
            steps,s%document%model_unit,changed,error,catalog=s%document%patches)
        if(len(error)>0)return
        if(changed)call refresh_rotated_scene(s)
        call result%insert('changed',changed)
        call result%insert('steps',size(steps))
        call result%insert('undo_available',document_undo_stamp(s%document)>0_int64)
    end subroutine transform_scene

    subroutine transform_selection(s,args,result,error)
        type(studio_t),intent(inout),target :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        type(affine_transform_t) :: transform
        character(:),allocatable :: operation, units
        real(dp) :: values(3)
        integer :: unit,axis,kind,id
        logical :: changed

        error='Choose move, rotate, scale or affine and supply the requested vectors or matrix.'
        operation=args%stringValue('op')
        select case(operation)
        case('move')
            if(.not.only_fields(args,'|op|offset|unit|'))return
            if(.not.vector_arg(args,'offset',transform%offset))return
        case('rotate')
            if(.not.only_fields(args,'|op|angles|pivot|unit|'))return
            if(.not.vector_arg(args,'angles',values))return
            transform%linear=rotation_matrix(values)
        case('scale')
            if(.not.only_fields(args,'|op|factors|pivot|unit|'))return
            if(.not.vector_arg(args,'factors',values))return
            if(any(values==0._dp))then
                error='Scale factors must be nonzero.'
                return
            end if
            do axis=1,3
                transform%linear(axis,axis)=values(axis)
            end do
        case('affine')
            if(.not.only_fields(args,'|op|matrix|offset|pivot|unit|'))return
            if(.not.transform_matrix_arg(args,transform%linear))return
            if(args%contains('offset'))then
                if(.not.vector_arg(args,'offset',transform%offset))return
            end if
        case default
            return
        end select

        if(args%contains('pivot'))then
            if(.not.vector_arg(args,'pivot',transform%pivot))return
        end if
        if(args%contains('unit'))then
            units=args%stringValue('unit')
            do unit=unitless+1,unit_ft
                if(units==trim(unit_symbol(unit)))exit
            end do
            if(unit>unit_ft)then
                error='Use mm, cm, m, in or ft for the offset and pivot.'
                return
            end if
            transform%offset=model_length(transform%offset,s%document%model_unit,unit)
            transform%pivot=model_length(transform%pivot,s%document%model_unit,unit)
        end if

        ! Transform joined geometry as its assembly so native interface validation stays in effect.
        kind=active_rank(s)
        if(kind>0)then
            id=active_patch_id(kind,catalog=s%document%patches)
            if(document_patch_in_assembly(s%document,kind,id))then
                error='This patch belongs to a multipatch. Select that assembly before transforming it.'
                return
            end if
        end if

        call apply_geometry_transform(s,transform,fit_after=.false.,changed_out=changed,error_out=error,reveal=.false.)
        if(len(error)>0)return
        call result%insert('changed',changed)
        call result%insert('operation',operation)
        call result%insert('undo_available',document_undo_stamp(s%document)>0_int64)
        call describe_selection(s,result)
    end subroutine transform_selection

    logical function parse_transform(args,step) result(valid)
        type(QJsonObject),intent(in) :: args
        type(scene_transform_step),intent(out) :: step
        character(:),allocatable :: name
        integer :: unit

        valid=.false.
        if(.not.only_fields(args,'|op|axis|value|unit|'))return
        step%axis=axis_number(args%stringValue('axis'),'XYZ')
        if(step%axis==0)return
        if(.not.number_arg(args,'value',step%value))return
        name=args%stringValue('op')
        select case(name)
        case('rotate')
            if(args%contains('unit'))return
            step%op=scene_rotate
        case('move')
            if(args%contains('unit'))then
                name=args%stringValue('unit')
                do unit=unitless+1,unit_ft
                    if(name==trim(unit_symbol(unit)))exit
                end do
                if(unit>unit_ft)return
                step%unit=unit
            end if
            step%op=scene_translate
        case default
            return
        end select
        valid=.true.
    end function parse_transform

    logical function transform_matrix_arg(args,matrix) result(valid)
        type(QJsonObject),intent(in) :: args
        real(dp),intent(out) :: matrix(3,3)
        type(QJsonArray) :: rows,row
        real(dp) :: values(3)
        integer :: i,j

        valid=.false.
        matrix=0._dp
        if(args%valueType('matrix')/=QJsonValue_Array)return
        rows=args%arrayValue('matrix')
        if(rows%size()==3)then
            valid=.true.
            do i=1,3
                if(rows%typeAt(i-1)/=QJsonValue_Array)then
                    valid=.false.
                    exit
                end if
                row=rows%arrayAt(i-1)
                valid=row%size()==3
                if(valid)valid=all([(row%typeAt(j-1)==QJsonValue_Double,j=1,3)])
                if(valid)then
                    values=row%toDoubles()
                    valid=all(ieee_is_finite(values))
                    if(valid)matrix(i,:)=values
                end if
                call row%delete()
                if(.not.valid)exit
            end do
        end if
        call rows%delete()
    end function transform_matrix_arg

end module studio_mcp_transforms
