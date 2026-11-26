! Exact numeric controls for ForCAD's native affine transform.
module studio_transform_dialog
    use studio_status_hints, only: set_status_hint, publish_window_status
    use qt, only: QComboBox, QDialog, QDialog_Accepted, QDialogButtonBox, QDialogButtonBox_AcceptRole, &
        QDialogButtonBox_Cancel, QDoubleSpinBox, QGridLayout, QHBoxLayout, QLabel, QPushButton, &
        QSizePolicy_Fixed, QSizePolicy_Ignored, QStackedWidget, Qt_PlainText, QVBoxLayout, QWidget
    use studio_icons, only: line_icon
    use, intrinsic :: iso_c_binding, only: dp => c_double
    use studio_forcad, only: affine_transform_t, rotation_matrix
    use studio_units, only: display_length, model_length, length_caption
    use studio_theme, only: studio_style_primary
    implicit none(type, external)
    private
    public :: edit_studio_transform

    character(3), parameter :: axes='XYZ'
    character(16), parameter :: operation_names(0:3)=[character(16) :: 'Move', 'Rotate', 'Scale', 'Transform matrix']
    character(86), parameter :: operation_hints(0:3)=[character(86) :: &
        'Enter a displacement along each world axis.', &
        'Angles are applied about world X, Y, then Z.', &
        'Enter a scale factor per world axis.', &
        'Enter the linear map by row. Translation is applied after the map about its pivot.']
    type :: transform_editor
        type(QDialog) :: dialog
        type(QComboBox) :: operation
        type(QWidget) :: pivot_row
        logical :: centered_pivot=.false.
        type(QStackedWidget) :: forms
        type(QDoubleSpinBox) :: offset(3),pivot(3),angles(3),scale(3),matrix(3,3)
    end type transform_editor
contains
    logical function edit_studio_transform(parent,name,transform,operation,initial_pivot,model_unit,display_unit) result(accepted)
        class(QWidget), intent(in) :: parent
        character(*), intent(in) :: name
        type(affine_transform_t), intent(out) :: transform
        integer, intent(in), optional :: operation,model_unit,display_unit
        real(dp), intent(in), optional :: initial_pivot(3)
        type(transform_editor), target :: editor
        type(QWidget) :: page
        type(QVBoxLayout) :: layout
        type(QGridLayout) :: grid
        type(QLabel) :: label
        type(QDialogButtonBox) :: buttons
        type(QPushButton) :: apply
        integer :: i,j,index
        logical :: show_translation
        real(dp) :: pivot_start(3)
        character(:), allocatable :: apply_text

        index=0
        apply_text='Apply transform'
        if(present(operation))then
            accepted=.false.
            if(operation<0.or.operation>3)return
            index=operation
            apply_text=trim(operation_names(index))
            if(index==3)apply_text='Apply matrix'
        end if
        show_translation=.not.present(operation).or.index==0.or.index==3
        editor%centered_pivot=present(initial_pivot)

        editor%dialog=QDialog(parent)
        call editor%dialog%setObjectName('studioTransformDialog')
        call editor%dialog%setWindowTitle('Transform geometry')
        call editor%dialog%setMinimumWidth(460)
        if(present(operation))call editor%dialog%setWindowTitle(trim(operation_names(index)))
        layout=QVBoxLayout(editor%dialog)
        call layout%setContentsMargins(18,16,18,16)
        call layout%setSpacing(10)
        label=QLabel(name,editor%dialog)
        call label%setTextFormat(Qt_PlainText)
        call label%setWordWrap(.true.)
        call label%setObjectName('studioDialogTitle')
        call layout%addWidget(label)
        editor%operation=QComboBox(editor%dialog)
        call editor%operation%setObjectName('studioTransformMode')
        call editor%operation%setAccessibleName('Transform operation')
        do i=0,3
            call editor%operation%addItem(trim(operation_names(i)))
        end do
        call editor%operation%setCurrentIndex(index)
        call editor%operation%setVisible(.not.present(operation))
        call layout%addWidget(editor%operation)
        editor%forms=QStackedWidget(editor%dialog)
        page=QWidget(); call editor%forms%addWidget(page)
        call add_vector(editor%forms,'Rotation (°)','Angle',0._dp,editor%angles)
        call add_vector(editor%forms,'Scale factors','Scale',1._dp,editor%scale)
        page=QWidget()
        grid=QGridLayout(page)
        call grid%setContentsMargins(0,0,0,0)
        call grid%setSpacing(6)
        do i=1,3
            label=QLabel(axes(i:i))
            call grid%addWidget(label,0,i)
            label=QLabel(axes(i:i))
            call grid%addWidget(label,i,0)
            do j=1,3
                call set_input(editor%matrix(i,j),'Matrix'//axes(i:i)//axes(j:j),transform%linear(i,j))
                call grid%addWidget(editor%matrix(i,j),i,j)
            end do
        end do
        call editor%forms%addWidget(page)
        call layout%addWidget(editor%forms)
        editor%pivot_row=vector_row(length_caption('Pivot',display_unit),'Pivot',0._dp,editor%pivot)
        call layout%addWidget(editor%pivot_row)
        page=vector_row(length_caption('Translation',display_unit),'Offset',0._dp,editor%offset)
        call page%setVisible(show_translation)
        call layout%addWidget(page)
        do i=1,3
            call editor%angles(i)%setSingleStep(5._dp)
            call editor%scale(i)%setSingleStep(.05_dp)
            call editor%pivot(i)%setSingleStep(.25_dp)
            call editor%offset(i)%setSingleStep(.25_dp)
            if(present(initial_pivot))call editor%pivot(i)%setValue(display_length(initial_pivot(i),model_unit,display_unit))
            pivot_start(i)=editor%pivot(i)%value()
        end do
        buttons=QDialogButtonBox(QDialogButtonBox_Cancel,editor%dialog)
        apply=buttons%addButton(apply_text,QDialogButtonBox_AcceptRole)
        call apply%setObjectName('studioTransformApply')
        call apply%setIcon(line_icon(35))
        call apply%setIconSize(16,16)
        call studio_style_primary(apply)
        call apply%setDefault(.true.)
        call layout%addWidget(buttons)
        call buttons%accepted%connect(accept_editor,editor)
        call buttons%rejected%connect(cancel_editor,editor)
        call editor%operation%currentIndexChanged%connect(mode_changed,editor)
        call mode_changed(index,editor)
        select case(index)
        case(0)
            call editor%offset(1)%setFocus()
        case(1)
            call editor%angles(1)%setFocus()
        case(2)
            call editor%scale(1)%setFocus()
        case(3)
            call editor%matrix(1,1)%setFocus()
        end select
        accepted=editor%dialog%exec()==QDialog_Accepted
        if(accepted)then
            if(show_translation)transform%offset=model_length([(editor%offset(i)%value(),i=1,3)],model_unit,display_unit)
            if(editor%operation%currentIndex()>0)then
                transform%pivot=model_length([(editor%pivot(i)%value(),i=1,3)],model_unit,display_unit)
                if(present(initial_pivot))then
                    ! An untouched pivot retains its exact model coordinates.
                    where([(editor%pivot(i)%value(),i=1,3)]==pivot_start)transform%pivot=initial_pivot
                end if
            end if
            select case(editor%operation%currentIndex())
            case(1)
                transform%linear=rotation_matrix([(editor%angles(i)%value(),i=1,3)])
            case(2)
                do i=1,3
                    transform%linear(i,i)=editor%scale(i)%value()
                end do
            case(3)
                do i=1,3
                    do j=1,3
                        transform%linear(i,j)=editor%matrix(i,j)%value()
                    end do
                end do
            end select
        end if
        call editor%dialog%delete()
    end function edit_studio_transform

    subroutine add_vector(stack,title,key,initial,fields)
        type(QStackedWidget), intent(inout) :: stack
        character(*),intent(in) :: title,key
        real(dp),intent(in) :: initial
        type(QDoubleSpinBox),intent(out) :: fields(3)
        type(QWidget) :: row
        row=vector_row(title,key,initial,fields)
        call stack%addWidget(row)
    end subroutine add_vector

    function vector_row(title,key,initial,fields) result(row)
        character(*),intent(in) :: title,key
        real(dp),intent(in) :: initial
        type(QDoubleSpinBox),intent(out) :: fields(3)
        type(QWidget) :: row
        type(QHBoxLayout) :: layout
        type(QLabel) :: label
        integer :: i
        row=QWidget()
        layout=QHBoxLayout(row)
        call layout%setContentsMargins(0,0,0,0)
        call layout%setSpacing(6)
        label=QLabel(title)
        call label%setMinimumWidth(104)
        call layout%addWidget(label)
        do i=1,3
            call set_input(fields(i),key//axes(i:i),initial)
            call fields(i)%setAccessibleName(title//' '//axes(i:i))
            call fields(i)%setPrefix(axes(i:i)//'  ')
            call layout%addWidget(fields(i),1)
        end do
    end function vector_row

    subroutine set_input(field,key,initial)
        type(QDoubleSpinBox),intent(out) :: field
        character(*),intent(in) :: key
        real(dp),intent(in) :: initial
        field=QDoubleSpinBox()
        call field%setObjectName('studioTransform'//key)
        call field%setAccessibleName(key)
        call field%setDecimals(6)
        call field%setRange(-huge(1._dp),huge(1._dp))
        call field%setValue(real(initial,kind(1.d0)))
        call field%setMinimumHeight(31)
        call field%setSizePolicy(QSizePolicy_Ignored,QSizePolicy_Fixed)
    end subroutine set_input

    subroutine mode_changed(index,data)
        integer,intent(in) :: index
        class(*),optional,intent(inout) :: data
        character(:),allocatable :: hint

        if(index<0.or.index>ubound(operation_hints,1))return
        select type(editor=>data)
        type is(transform_editor)
            hint=trim(operation_hints(index))
            if(index>0)then
                if(editor%centered_pivot)then
                    hint=hint//' The initial pivot is the center of the selection.'
                else
                    hint=hint//' The initial pivot is the origin.'
                end if
            end if
            call set_status_hint(editor%dialog,hint)
            call publish_window_status(editor%dialog,hint)
            call editor%forms%setCurrentIndex(index)
            call editor%forms%setVisible(index>0)
            call editor%pivot_row%setVisible(index>0)
        end select
    end subroutine mode_changed

    subroutine accept_editor(data)
        class(*),optional,intent(inout) :: data
        select type(editor=>data)
        type is(transform_editor)
            call editor%dialog%accept()
        end select
    end subroutine accept_editor

    subroutine cancel_editor(data)
        class(*),optional,intent(inout) :: data
        select type(editor=>data)
        type is(transform_editor)
            call editor%dialog%reject()
        end select
    end subroutine cancel_editor
end module studio_transform_dialog
