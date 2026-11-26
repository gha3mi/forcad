! Create simple scene objects with ForCAD's native constructors.
module studio_primitive_dialog
    use studio_status_hints, only: set_status_hint, publish_window_status
    use studio_forcad, only: rk, primitive_names, &
        point_shape, line_shape, circle_shape, arc_shape, rectangle_shape, box_shape, cylinder_shape, tube_shape
    use qt, only: QComboBox, QDialog, QDialog_Accepted, QDialogButtonBox, QDialogButtonBox_AcceptRole, QDialogButtonBox_Cancel, &
        QGridLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QRegularExpressionValidator, QVBoxLayout, QWidget
    use studio_icons, only: line_icon
    use studio_construction, only: construct_studio_primitive
    use studio_single_patch, only: studio_patch_catalog, single_patch_count, single_patch_name
    use studio_parameter_editor, only: parse_parameter_number,parameter_number_pattern
    use studio_units, only: unitless,model_length,length_caption
    use studio_theme, only: studio_style_primary
    implicit none(type, external)
    private
    public :: create_studio_primitive

    type :: primitive_editor
        type(QDialog) :: dialog
        type(QComboBox) :: choice
        type(QLineEdit) :: name,coordinates(3),inputs(3)
        type(QLabel) :: heading, position_caption, input_labels(3)
        type(QRegularExpressionValidator) :: decimal
        type(studio_patch_catalog), pointer :: catalog=>null()
        character(:), allocatable :: suggested_name
        integer :: created_kind=0,model_unit=unitless,display_unit=unitless
    end type primitive_editor

contains

    integer function create_studio_primitive(parent,primitive,model_unit,display_unit,catalog) result(kind)
        class(QWidget), intent(in) :: parent
        integer, intent(in) :: primitive
        integer, intent(in), optional :: model_unit,display_unit
        type(studio_patch_catalog), intent(inout), target, optional :: catalog
        type(primitive_editor), target :: editor
        type(QVBoxLayout) :: body
        type(QHBoxLayout) :: positions
        type(QGridLayout) :: form
        type(QLabel) :: label, shape_caption
        type(QPushButton) :: add
        type(QDialogButtonBox) :: buttons
        character(len=3), parameter :: axes='XYZ'
        integer :: i

        kind=0
        if(primitive<0.or.primitive>size(primitive_names))return
        if(present(catalog))editor%catalog=>catalog
        if(present(model_unit))editor%model_unit=model_unit
        if(present(display_unit))editor%display_unit=display_unit
        editor%dialog=QDialog(parent)
        call editor%dialog%setObjectName('studioPrimitiveDialog')
        call editor%dialog%setWindowTitle('Create geometry')
        call editor%dialog%setMinimumWidth(440)
        editor%decimal=QRegularExpressionValidator(parameter_number_pattern,editor%dialog)
        body=QVBoxLayout(editor%dialog)
        call body%setContentsMargins(18,16,18,16)
        call body%setSpacing(10)
        editor%heading=QLabel('',editor%dialog)
        call editor%heading%setObjectName('studioDialogTitle')
        call body%addWidget(editor%heading)

        form=QGridLayout()
        call form%setColumnStretch(1,1)
        call form%setSpacing(8)
        shape_caption=QLabel('Shape',editor%dialog)
        call shape_caption%setVisible(primitive==0)
        call form%addWidget(shape_caption,0,0)
        editor%choice=QComboBox(editor%dialog)
        call editor%choice%setObjectName('studioPrimitiveKind')
        call editor%choice%setProperty('accessibleName','Primitive type')
        do i=1,size(primitive_names)
            call editor%choice%addItem(trim(primitive_names(i)))
        end do
        call editor%choice%setCurrentIndex(max(0,primitive-1))
        call editor%choice%setVisible(primitive==0)
        call form%addWidget(editor%choice,0,1)

        label=QLabel('Name',editor%dialog)
        call form%addWidget(label,1,0)
        editor%name=QLineEdit(parent=editor%dialog)
        call editor%name%setObjectName('studioPrimitiveName')
        call editor%name%setProperty('accessibleName','Object name')
        call form%addWidget(editor%name,1,1)

        editor%position_caption=QLabel('',editor%dialog)
        call form%addWidget(editor%position_caption,2,0)
        positions=QHBoxLayout()
        call positions%setSpacing(5)
        do i=1,3
            label=QLabel(axes(i:i),editor%dialog)
            call label%setProperty('dialogRole','detail')
            call positions%addWidget(label)
            call create_input(editor%coordinates(i),editor%dialog,'studioPrimitiveOrigin'//axes(i:i),axes(i:i)//' coordinate')
            call editor%coordinates(i)%setText('0')
            call positions%addWidget(editor%coordinates(i),1)
        end do
        call form%addLayout(positions,2,1)

        do i=1,3
            editor%input_labels(i)=QLabel('',editor%dialog)
            call form%addWidget(editor%input_labels(i),2+i,0)
            call create_input(editor%inputs(i),editor%dialog,'studioPrimitiveValue'//achar(iachar('0')+i),'Shape value')
            call form%addWidget(editor%inputs(i),2+i,1)
        end do
        call body%addLayout(form)

        buttons=QDialogButtonBox(QDialogButtonBox_Cancel,editor%dialog)
        add=buttons%addButton('Create',QDialogButtonBox_AcceptRole)
        call add%setIcon(line_icon(30))
        call add%setIconSize(17,17)
        call add%setObjectName('studioPrimitiveApply')
        call add%setDefault(.true.)
        call studio_style_primary(add)
        call body%addWidget(buttons)
        call update_editor(editor,max(1,primitive))
        call editor%choice%currentIndexChanged%connect(on_shape_changed,editor)
        call editor%name%textChanged%connect(input_changed,editor)
        do i=1,3
            call editor%coordinates(i)%textChanged%connect(input_changed,editor)
            call editor%inputs(i)%textChanged%connect(input_changed,editor)
        end do
        call buttons%accepted%connect(on_create,editor)
        call buttons%rejected%connect(on_cancel,editor)
        if(primitive/=0)call editor%name%setFocus()
        if(editor%dialog%exec()==QDialog_Accepted)kind=editor%created_kind
        call editor%dialog%delete()
    end function create_studio_primitive

    subroutine create_input(field,parent,name,description)
        type(QLineEdit), intent(out) :: field
        class(QWidget), intent(in) :: parent
        character(*), intent(in) :: name,description

        field=QLineEdit(parent=parent)
        call field%setObjectName(name)
        call field%setMaxLength(128)
        call field%setMinimumHeight(27)
        call field%setProperty('accessibleName',description)
    end subroutine create_input

    subroutine update_editor(editor,primitive)
        type(primitive_editor), intent(inout) :: editor
        integer, intent(in) :: primitive
        character(len=20) :: captions(3)
        character(len=6) :: defaults(3)
        character(:), allocatable :: position,note,name
        character(:), allocatable :: caption
        integer :: i
        logical :: new_name,visible

        new_name=.not.allocated(editor%suggested_name)
        if(allocated(editor%suggested_name))new_name=editor%name%text()==editor%suggested_name
        name=next_primitive_name(trim(primitive_names(primitive)),editor%catalog)
        editor%suggested_name=name
        if(new_name)call editor%name%setText(name)

        captions=''
        defaults='1'
        position='Origin'
        select case(primitive)
        case(point_shape)
            position='Position'
            note='Place a point at an exact XYZ coordinate.'
        case(line_shape)
            position='Start'
            note='Set the endpoints of a straight segment.'
            captions=[character(len=20)::'End X','End Y','End Z']
            defaults=[character(len=6)::'1','0','0']
        case(circle_shape)
            position='Center'
            note='Create a circle in the XY plane at the center’s Z coordinate.'
            captions(1)='Radius'
        case(arc_shape)
            position='Center'
            note='Create an arc in a plane parallel to XY. Angles are measured from +X.'
            captions=[character(len=20)::'Radius','Start angle (°)','Sweep (°)']
            defaults=[character(len=6)::'1','0','90']
        case(rectangle_shape)
            note='Place a rectangular surface parallel to XY.'
            captions(1:2)=[character(len=20)::'Width (X)','Height (Y)']
        case(box_shape)
            note='The origin is the minimum XYZ corner.'
            captions=[character(len=20)::'Width (X)','Depth (Y)','Height (Z)']
        case(cylinder_shape)
            position='Base center'
            note='Extrude a circle along +Z to make a cylindrical surface.'
            captions(1:2)=[character(len=20)::'Radius','Height (Z)']
        case(tube_shape)
            position='Base center'
            note='Create an annular volume along +Z.'
            captions=[character(len=20)::'Inner radius','Outer radius','Height (Z)']
            defaults=[character(len=6)::'0.5','1','1']
        case default
            return
        end select
        call editor%heading%setText('Create '//trim(primitive_names(primitive)))
        call set_status_hint(editor%dialog,note)
        call editor%position_caption%setText(length_caption(position,editor%display_unit))
        do i=1,3
            visible=len_trim(captions(i))>0
            call editor%input_labels(i)%setVisible(visible)
            call editor%inputs(i)%setVisible(visible)
            if(.not.visible)cycle
            caption=trim(captions(i))
            if(spatial_input(primitive,i))caption=length_caption(caption,editor%display_unit)
            call editor%input_labels(i)%setText(caption)
            call editor%inputs(i)%setText(trim(defaults(i)))
            call editor%inputs(i)%setProperty('accessibleName',caption)
        end do
        call publish_window_status(editor%dialog,'',blocking=.true.)
    end subroutine update_editor

    pure logical function spatial_input(primitive,column)
        integer, intent(in) :: primitive,column
        spatial_input=primitive/=arc_shape.or.column==1
    end function spatial_input

    function next_primitive_name(prefix,catalog) result(name)
        character(*), intent(in) :: prefix
        type(studio_patch_catalog), intent(in), optional :: catalog
        character(:), allocatable :: name
        character(len=20) :: number
        integer :: rank,item,serial

        serial=0
        names: do
            serial=serial+1
            write(number,'(i0)')serial
            name=prefix//' '//trim(number)
            do rank=1,3
                do item=1,single_patch_count(rank,catalog)
                    if(single_patch_name(rank,item,catalog)==name)cycle names
                end do
            end do
            return
        end do names
    end function next_primitive_name

    subroutine on_shape_changed(index,data)
        integer, intent(in) :: index
        class(*), intent(inout), optional :: data

        select type(editor=>data)
        type is(primitive_editor)
            if(index>=0)call update_editor(editor,index+1)
        end select
    end subroutine on_shape_changed

    subroutine input_changed(text,data)
        character(*), intent(in) :: text
        class(*), optional, intent(inout) :: data

        select type(editor=>data)
        type is(primitive_editor)
            call publish_window_status(editor%dialog,'',blocking=.true.)
        end select
    end subroutine input_changed

    subroutine on_create(data)
        class(*), intent(inout), optional :: data
        character(:), allocatable :: error
        real(rk) :: origin(3),values(3)
        integer :: i,variant,primitive

        select type(editor=>data)
        type is(primitive_editor)
            call publish_window_status(editor%dialog,'',blocking=.true.)
            primitive=editor%choice%currentIndex()+1
            values=0._rk
            do i=1,3
                if(.not.parse_parameter_number(editor%decimal,editor%coordinates(i)%text(),origin(i)))then
                    call show_error(editor,'Enter a finite coordinate.',editor%coordinates(i))
                    return
                end if
                origin(i)=model_length(origin(i),editor%model_unit,editor%display_unit)
                if(.not.editor%inputs(i)%isVisible())cycle
                if(.not.parse_parameter_number(editor%decimal,editor%inputs(i)%text(),values(i)))then
                    call show_error(editor,'Enter a finite '//editor%input_labels(i)%text()//' value.',editor%inputs(i))
                    return
                end if
                if(spatial_input(primitive,i))values(i)=model_length(values(i),editor%model_unit,editor%display_unit)
            end do
            call construct_studio_primitive(primitive,editor%name%text(),origin,values, &
                editor%created_kind,variant,error,editor%catalog)
            if(editor%created_kind>0)then
                call editor%dialog%accept()
            else
                call show_error(editor,error)
            end if
        end select
    end subroutine on_create

    subroutine show_error(editor,error,field)
        type(primitive_editor), intent(inout) :: editor
        character(*), intent(in) :: error
        type(QLineEdit), intent(inout), optional :: field

        call publish_window_status(editor%dialog,error,warning=.true.,blocking=.true.)
        if(present(field))call field%setFocus()
    end subroutine show_error

    subroutine on_cancel(data)
        class(*), intent(inout), optional :: data

        select type(editor=>data)
        type is(primitive_editor)
            call editor%dialog%reject()
        end select
    end subroutine on_cancel


end module studio_primitive_dialog
