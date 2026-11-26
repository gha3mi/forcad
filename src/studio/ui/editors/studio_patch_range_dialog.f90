! Select a literal parameter cut or bounded range for a native modeling operation.
module studio_patch_range_dialog
    use studio_status_hints, only: set_status_hint, publish_window_status
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use qt, only: QAbstractButton, QButtonGroup, QDialog, QDialog_Accepted, QDialogButtonBox, QDialogButtonBox_AcceptRole, &
        QDialogButtonBox_Cancel, QFormLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QRegularExpressionValidator, &
        QSizePolicy_Fixed, QSizePolicy_Ignored, Qt_PlainText, QVBoxLayout, QWidget
    use studio_icons, only: line_icon
    use studio_parameter_editor, only: parse_parameter_number, number_text, parameter_number_pattern, &
        parameter_direction_buttons
    use studio_theme, only: studio_style_primary
    implicit none(type, external)
    private
    public :: edit_studio_patch_range, edit_studio_patch_section, edit_studio_parameter_interval, choose_studio_reverse_direction

    integer, parameter :: split_mode=1, range_mode=2, section_mode=3, interval_mode=4, reverse_mode=5

    type :: range_draft
        character(:), allocatable :: lower, upper
    end type range_draft

    type :: patch_range_editor
        type(QDialog) :: dialog
        type(range_draft) :: drafts(3)
        type(QLineEdit) :: lower_field, upper_field
        type(QLabel) :: domain
        type(QRegularExpressionValidator) :: decimal
        real(c_double) :: bounds(2,3), probe(3)
        real(c_double) :: lower, upper
        integer :: mode,kind,axis=0
    end type patch_range_editor

contains

    logical function edit_studio_patch_range(parent,title,kind,axis,bounds,probe,split,lower,upper) result(accepted)
        class(QWidget), intent(in) :: parent
        character(*), intent(in) :: title
        integer, intent(in) :: kind
        integer, intent(inout) :: axis
        real(c_double), intent(in) :: bounds(2,3), probe(3)
        logical, intent(in) :: split
        real(c_double), intent(out) :: lower, upper

        accepted=edit_parameter_extent(parent,title,kind,axis,bounds,probe,merge(split_mode,range_mode,split),lower,upper)
    end function edit_studio_patch_range

    logical function edit_studio_patch_section(parent,title,kind,axis,bounds,probe,parameter) result(accepted)
        class(QWidget), intent(in) :: parent
        character(*), intent(in) :: title
        integer, intent(in) :: kind
        integer, intent(inout) :: axis
        real(c_double), intent(in) :: bounds(2,3), probe(3)
        real(c_double), intent(out) :: parameter
        real(c_double) :: unused_upper

        accepted=edit_parameter_extent(parent,title,kind,axis,bounds,probe,section_mode,parameter,unused_upper)
    end function edit_studio_patch_section

    logical function edit_studio_parameter_interval(parent,title,kind,axis,bounds,lower,upper) result(accepted)
        class(QWidget), intent(in) :: parent
        character(*), intent(in) :: title
        integer, intent(in) :: kind
        integer, intent(inout) :: axis
        real(c_double), intent(in) :: bounds(2,3)
        real(c_double), intent(out) :: lower, upper

        accepted=edit_parameter_extent(parent,title,kind,axis,bounds,bounds(1,:),interval_mode,lower,upper)
    end function edit_studio_parameter_interval

    logical function choose_studio_reverse_direction(parent,title,kind,axis,bounds) result(accepted)
        class(QWidget), intent(in) :: parent
        character(*), intent(in) :: title
        integer, intent(in) :: kind
        integer, intent(inout) :: axis
        real(c_double), intent(in) :: bounds(2,3)
        real(c_double) :: lower,upper

        accepted=edit_parameter_extent(parent,title,kind,axis,bounds,bounds(1,:),reverse_mode,lower,upper)
    end function choose_studio_reverse_direction

    logical function edit_parameter_extent(parent,title,kind,axis,bounds,probe,mode,lower,upper) result(accepted)
        class(QWidget), intent(in) :: parent
        character(*), intent(in) :: title
        integer, intent(in) :: kind,mode
        integer, intent(inout) :: axis
        real(c_double), intent(in) :: bounds(2,3), probe(3)
        real(c_double), intent(out) :: lower, upper

        type(patch_range_editor), target :: editor
        type(QVBoxLayout) :: layout
        type(QHBoxLayout) :: fields
        type(QFormLayout) :: form
        type(QDialogButtonBox) :: buttons
        type(QPushButton) :: apply
        type(QLabel) :: label
        type(QWidget) :: direction_choices,bounds_row
        type(QButtonGroup) :: directions
        type(QAbstractButton) :: selected_direction
        character(:), allocatable :: heading, action, explanation, lower_caption, upper_caption, object_name
        integer :: icon_kind
        logical :: has_upper

        accepted=.false.
        lower=0._c_double
        upper=0._c_double
        if(kind<1.or.kind>3.or.axis<1.or.axis>kind)return
        if(.not.all(ieee_is_finite(bounds(:,:kind))))return
        if(any(bounds(1,:kind)>=bounds(2,:kind)))return

        editor%bounds=bounds
        editor%probe=probe
        editor%mode=mode
        editor%kind=kind
        has_upper=mode==range_mode.or.mode==interval_mode
        upper_caption='Upper'
        object_name='studioPatchRangeDialog'
        icon_kind=91
        select case(mode)
        case(split_mode)
            heading='Split patch'
            lower_caption='Split at'
            action='Split'
            icon_kind=90
            explanation='Split at an interior value. The lower piece keeps the name.'
        case(section_mode)
            if(kind<2)return
            heading='Extract section'
            lower_caption='Section at'
            if(kind==2)then
                action='Extract curve'
                explanation='Create an exact curve at a constant surface parameter.'
            else
                action='Extract surface'
                explanation='Create an exact surface at a constant volume parameter.'
            end if
            object_name='studioPatchSectionDialog'
        case(interval_mode)
            heading='Set domain'
            lower_caption='New lower'
            upper_caption='New upper'
            action='Set domain'
            explanation='Change the bounds while preserving the patch shape.'
            object_name='studioParameterIntervalDialog'
        case(reverse_mode)
            heading='Reverse direction'
            action='Reverse'
            explanation='Reverses the chosen parameter direction while preserving the patch shape.'
            object_name='studioReverseDirectionDialog'
            icon_kind=92
        case default
            heading='Extract range'
            lower_caption='Lower'
            action='Extract range'
            explanation='Extract a new patch and keep the source.'
        end select

        editor%dialog=QDialog(parent)
        call editor%dialog%setObjectName(object_name)
        call editor%dialog%setWindowTitle(heading)
        call editor%dialog%setMinimumWidth(450)
        layout=QVBoxLayout(editor%dialog)
        call layout%setContentsMargins(18,16,18,16)
        call layout%setSpacing(10)
        label=QLabel(heading,editor%dialog)
        call label%setObjectName('studioDialogTitle')
        call layout%addWidget(label)
        label=QLabel(title,editor%dialog)
        call label%setObjectName('muted')
        call label%setTextFormat(Qt_PlainText)
        call label%setWordWrap(.true.)
        call layout%addWidget(label)

        form=QFormLayout()
        call form%setSpacing(8)
        direction_choices=parameter_direction_buttons(editor%dialog,kind,axis,'studioPatchRangeAxis', &
            directions,direction_changed,editor)
        call form%addRow('Direction',direction_choices)
        if(mode/=reverse_mode)then
            editor%domain=QLabel('',editor%dialog)
            call editor%domain%setObjectName('studioPatchRangeDomain')
            call editor%domain%setTextFormat(Qt_PlainText)
            call editor%domain%setWordWrap(.true.)
            call editor%domain%setProperty('dialogRole','detail')
            call form%addRow('Current range',editor%domain)
            if(has_upper)then
                bounds_row=QWidget(editor%dialog)
                fields=QHBoxLayout(bounds_row)
                call fields%setContentsMargins(0,0,0,0)
                call fields%setSpacing(8)
                editor%lower_field=QLineEdit(parent=bounds_row)
                editor%upper_field=QLineEdit(parent=bounds_row)
                call configure_field(editor%lower_field,'studioPatchRangeLower',lower_caption//' bound')
                call configure_field(editor%upper_field,'studioPatchRangeUpper',upper_caption//' bound')
                call editor%lower_field%setSizePolicy(QSizePolicy_Ignored,QSizePolicy_Fixed)
                call editor%upper_field%setSizePolicy(QSizePolicy_Ignored,QSizePolicy_Fixed)
                call editor%lower_field%setMinimumWidth(100)
                call editor%upper_field%setMinimumWidth(100)
                call fields%addWidget(QLabel(lower_caption,bounds_row))
                call fields%addWidget(editor%lower_field,1)
                call fields%addSpacing(6)
                call fields%addWidget(QLabel(upper_caption,bounds_row))
                call fields%addWidget(editor%upper_field,1)
                call form%addRow(bounds_row)
            else
                editor%lower_field=QLineEdit(parent=editor%dialog)
                call configure_field(editor%lower_field,'studioPatchRangeLower',lower_caption)
                call form%addRow(lower_caption,editor%lower_field)
            end if
        end if
        call layout%addLayout(form)

        call set_status_hint(editor%dialog,explanation)
        if(mode/=reverse_mode)then
            editor%decimal=QRegularExpressionValidator(parameter_number_pattern,editor%dialog)
        end if
        buttons=QDialogButtonBox(QDialogButtonBox_Cancel,editor%dialog)
        apply=buttons%addButton(action,QDialogButtonBox_AcceptRole)
        call apply%setObjectName('studioPatchRangeApply')
        call apply%setIcon(line_icon(icon_kind))
        call apply%setIconSize(17,17)
        call apply%setMinimumWidth(88)
        call apply%setDefault(.true.)
        call studio_style_primary(apply)
        call apply%setAccessibleDescription(explanation)
        call layout%addWidget(buttons)
        call buttons%accepted%connect(apply_range,editor)
        call buttons%rejected%connect(cancel_range,editor)
        call direction_changed(axis-1,editor)
        if(mode==reverse_mode)then
            if(kind>1)then
                selected_direction=directions%button(axis-1)
                call selected_direction%setFocus()
            end if
        else
            call editor%lower_field%textChanged%connect(clear_feedback,editor)
            if(has_upper)call editor%upper_field%textChanged%connect(clear_feedback,editor)
            call editor%lower_field%setFocus()
        end if

        if(editor%dialog%exec()==QDialog_Accepted)then
            axis=editor%axis
            if(mode/=reverse_mode)then
                lower=editor%lower
                upper=editor%upper
            end if
            accepted=.true.
        end if
        call editor%dialog%delete()
    end function edit_parameter_extent

    subroutine configure_field(field,object_name,label)
        type(QLineEdit), intent(inout) :: field
        character(*), intent(in) :: object_name, label

        call field%setObjectName(object_name)
        call field%setAccessibleName(label)
        call field%setMinimumHeight(30)
        call field%setMaxLength(128)
    end subroutine configure_field

    subroutine direction_changed(index,data)
        integer, intent(in) :: index
        class(*), optional, intent(inout) :: data
        integer :: axis
        real(c_double) :: first, last, value
        logical :: has_upper

        select type(editor=>data)
        type is(patch_range_editor)
            axis=index+1
            if(axis<1.or.axis>editor%kind.or.axis==editor%axis)return
            if(editor%mode==reverse_mode)then
                editor%axis=axis
                return
            end if

            has_upper=.not.editor%upper_field%isNull()
            if(editor%axis>0)then
                editor%drafts(editor%axis)%lower=editor%lower_field%text()
                if(has_upper)editor%drafts(editor%axis)%upper=editor%upper_field%text()
            end if
            editor%axis=axis
            first=editor%bounds(1,axis)
            last=editor%bounds(2,axis)
            call editor%domain%setText(number_text(first)//' to '//number_text(last))
            if(allocated(editor%drafts(axis)%lower))then
                call editor%lower_field%setText(editor%drafts(axis)%lower)
                if(has_upper)call editor%upper_field%setText(editor%drafts(axis)%upper)
            else if(has_upper)then
                call editor%lower_field%setText(number_text(first))
                call editor%upper_field%setText(number_text(last))
            else
                value=editor%probe(axis)
                if(editor%mode==split_mode)then
                    if(.not.(value>first.and.value<last))value=0.5_c_double*first+0.5_c_double*last
                else
                    if(.not.(value>=first.and.value<=last))value=0.5_c_double*first+0.5_c_double*last
                end if
                call editor%lower_field%setText(number_text(value))
            end if
            call publish_window_status(editor%dialog,'',blocking=.true.)
        end select
    end subroutine direction_changed

    subroutine clear_feedback(input,data)
        character(*), intent(in) :: input
        class(*), optional, intent(inout) :: data

        select type(editor=>data)
        type is(patch_range_editor)
            call publish_window_status(editor%dialog,'',blocking=.true.)
        end select
    end subroutine clear_feedback

    subroutine apply_range(data)
        class(*), optional, intent(inout) :: data
        integer :: axis
        real(c_double) :: first, last, lower, upper

        select type(editor=>data)
        type is(patch_range_editor)
            if(editor%mode==reverse_mode)then
                call editor%dialog%accept()
                return
            end if
            axis=editor%axis
            first=editor%bounds(1,axis)
            last=editor%bounds(2,axis)
            if(.not.parse_parameter_number(editor%decimal,editor%lower_field%text(),lower))then
                if(editor%upper_field%isNull())then
                    call show_error(editor,'Enter a finite parameter value.')
                else
                    call show_error(editor,'Enter a finite lower parameter bound.')
                end if
                call editor%lower_field%setFocus()
                return
            end if
            upper=lower
            select case(editor%mode)
            case(split_mode)
                if(lower<=first.or.lower>=last)then
                    call show_error(editor,'The split parameter must be strictly inside the patch domain.')
                    call editor%lower_field%setFocus()
                    return
                end if
            case(section_mode)
                if(lower<first.or.lower>last)then
                    call show_error(editor,'The section parameter must lie within the patch domain.')
                    call editor%lower_field%setFocus()
                    return
                end if
            case default
                if(.not.parse_parameter_number(editor%decimal,editor%upper_field%text(),upper))then
                    call show_error(editor,'Enter a finite upper parameter bound.')
                    call editor%upper_field%setFocus()
                    return
                end if
                if(lower>=upper)then
                    call show_error(editor,'The lower bound must be less than the upper bound.')
                    call editor%upper_field%setFocus()
                    return
                end if
                if(editor%mode==range_mode.and.(lower<first.or.upper>last))then
                    call show_error(editor,'Both bounds must lie within the patch domain.')
                    if(lower<first)then
                        call editor%lower_field%setFocus()
                    else
                        call editor%upper_field%setFocus()
                    end if
                    return
                end if
            end select
            editor%lower=lower
            editor%upper=upper
            call editor%dialog%accept()
        end select
    end subroutine apply_range

    subroutine show_error(editor,message)
        type(patch_range_editor), intent(inout) :: editor
        character(*), intent(in) :: message

        call publish_window_status(editor%dialog,message,warning=.true.,blocking=.true.)
    end subroutine show_error

    subroutine cancel_range(data)
        class(*), optional, intent(inout) :: data

        select type(editor=>data)
        type is(patch_range_editor)
            call editor%dialog%reject()
        end select
    end subroutine cancel_range

end module studio_patch_range_dialog
