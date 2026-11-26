! Native degree elevation and reduction for the selected ForCAD direction.
module studio_degree_editor
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use qt, only: QAbstractButton, QButtonGroup, QDialog, QDialog_Accepted, QDialogButtonBox, QDialogButtonBox_AcceptRole, &
        QDialogButtonBox_Cancel, QFrame, QGridLayout, QLabel, QLineEdit, QPushButton, &
        QRegularExpressionValidator, QSpinBox, Qt_AlignHCenter, Qt_AlignLeft, QValidator_Acceptable, &
        QVBoxLayout, QWidget
    use studio_single_patch, only: studio_patch_catalog, single_patch_info, read_single_patch, commit_single_patch_degree, &
        single_patch_refinement_limits

    use studio_theme, only: studio_style_primary
    use studio_parameter_editor, only: parameter_direction_buttons
    use studio_status_hints, only: set_status_hint, publish_window_status
    implicit none(type, external)
    private
    public :: edit_single_patch_degree

    type :: degree_editor
        type(QDialog) :: dialog
        type(QButtonGroup) :: direction
        type(QSpinBox) :: amount
        type(QLineEdit) :: tolerance
        type(QPushButton) :: apply
        type(QRegularExpressionValidator) :: decimal
        type(QLabel) :: current_degree, target_degree
        integer :: kind, axis, patch_id, degrees(3), limits(3)=0, accepted=0, requested=0
        real(c_double) :: error_bound=0._c_double
        type(studio_patch_catalog), pointer :: catalog=>null()
        logical :: reducing
    end type degree_editor

contains

    logical function edit_single_patch_degree(parent,kind,axis,probe,reducing,accepted,requested,error_bound, &
        catalog,patch_id) result(changed)
        class(QWidget), intent(in) :: parent
        integer, intent(in) :: kind,patch_id
        integer, intent(inout) :: axis
        real(c_double), intent(in) :: probe(3)
        logical, intent(in) :: reducing
        integer, intent(out) :: accepted,requested
        real(c_double), intent(out) :: error_bound
        type(studio_patch_catalog), intent(inout), target :: catalog
        type(degree_editor), target :: editor
        type(single_patch_info) :: patch
        type(QVBoxLayout) :: layout
        type(QGridLayout) :: form
        type(QLabel) :: label
        type(QDialogButtonBox) :: buttons
        type(QPushButton) :: cancel
        type(QWidget) :: direction_choices
        type(QAbstractButton) :: axis_button
        character(len=:), allocatable :: action,detail,amount_label
        integer :: direction

        changed=.false.
        accepted=0
        requested=0
        error_bound=0._c_double
        if (kind<1.or.kind>3.or.axis<1.or.axis>kind) return
        call read_single_patch(kind,probe,axis,patch,patch_id=patch_id,catalog=catalog)
        ! Degree metadata stays usable when a probe derivative is unavailable.
        if(patch%dimension/=kind.or..not.allocated(patch%knots))return
        if(size(patch%knots)==0)return
        editor%kind=kind
        editor%axis=0
        editor%patch_id=patch_id
        editor%catalog=>catalog
        editor%degrees=patch%degree
        editor%reducing=reducing
        if (reducing) then
            action='Reduce'
            amount_label='Degrees to remove'
            editor%limits=editor%degrees
            detail='Tolerance limits the total shape change. Fewer steps may meet it.'
        else
            action='Elevate'
            amount_label='Degrees to add'
            do direction=1,kind
                call single_patch_refinement_limits(kind,direction,max_steps=editor%limits(direction),patch_id=patch_id,catalog=catalog)
            end do
            detail='The patch keeps its shape. Higher degrees add control points.'
        end if

        editor%dialog=QDialog(parent)
        call editor%dialog%setObjectName('studioDegreeEditor')
        call editor%dialog%setWindowTitle(action//' degree')
        call editor%dialog%setMinimumWidth(472)
        layout=QVBoxLayout(editor%dialog)
        call layout%setContentsMargins(18,16,18,16)
        call layout%setSpacing(10)

        label=QLabel(action//' degree',editor%dialog)
        call label%setObjectName('studioDialogTitle')
        call layout%addWidget(label)
        label=QLabel(patch%name,editor%dialog)
        call label%setObjectName('studioDegreeDirection')
        call label%setProperty('textFormat',0)
        call label%setWordWrap(.true.)
        call label%setProperty('dialogRole','detail')
        call layout%addWidget(label)

        form=QGridLayout()
        call form%setSpacing(8)
        call form%setColumnStretch(2,1)
        label=QLabel('Direction',editor%dialog)
        call form%addWidget(label,0,0)
        direction_choices=parameter_direction_buttons(editor%dialog,kind,axis,'studioDegreeAxis', &
            editor%direction,direction_changed,editor)
        call form%addWidget(direction_choices,0,1,1,2)
        label=QLabel(amount_label,editor%dialog)
        call form%addWidget(label,1,0)
        editor%amount=QSpinBox(editor%dialog)
        call editor%amount%setObjectName('studioDegreeAmount')
        call editor%amount%setAccessibleName(amount_label)
        call editor%amount%setMinimumHeight(30)
        call editor%amount%setFixedWidth(96)
        call form%addWidget(editor%amount,1,1,alignment=Qt_AlignLeft)
        if (reducing) then
            label=QLabel('Error tolerance',editor%dialog)
            call form%addWidget(label,2,0)
            editor%tolerance=QLineEdit('1e-10',editor%dialog)
            call editor%tolerance%setObjectName('studioDegreeTolerance')
            call editor%tolerance%setAccessibleName('Physical error tolerance in model units')
            call set_status_hint(editor%tolerance,'Maximum permitted shape change in model units.')
            call editor%tolerance%setMinimumHeight(30)
            call editor%tolerance%setMaxLength(48)
            call form%addWidget(editor%tolerance,2,1)
            label=QLabel('model units',editor%dialog)
            call label%setProperty('dialogRole','detail')
            call form%addWidget(label,2,2)
            editor%decimal=QRegularExpressionValidator( &
                '[+-]?(?:[0-9]+(?:[.,][0-9]*)?|[.,][0-9]+)(?:[eE][+-]?[0-9]+)?',editor%dialog)
            call editor%tolerance%textChanged%connect(tolerance_changed,editor)
        end if
        call layout%addLayout(form)
        call add_degree_summary(layout,editor)
        call set_status_hint(editor%dialog,detail)

        buttons=QDialogButtonBox(QDialogButtonBox_Cancel,editor%dialog)
        editor%apply=buttons%addButton(action,QDialogButtonBox_AcceptRole)
        call editor%apply%setObjectName('studioDegreeApply')
        call studio_style_primary(editor%apply)
        call editor%apply%setMinimumWidth(80)
        cancel=buttons%button(QDialogButtonBox_Cancel)
        call cancel%setObjectName('studioDegreeCancel')
        call cancel%setMinimumWidth(80)
        call layout%addWidget(buttons)
        call buttons%accepted%connect(apply_degree,editor)
        call buttons%rejected%connect(cancel_degree,editor)
        call editor%amount%valueChanged%connect(amount_changed,editor)
        call direction_changed(axis-1,editor)
        if (editor%limits(axis)>0) then
            call editor%amount%setFocus()
        else if (kind>1) then
            axis_button=editor%direction%button(axis-1)
            call axis_button%setFocus()
        end if

        if (editor%dialog%exec()==QDialog_Accepted) then
            accepted=editor%accepted
            requested=editor%requested
            error_bound=editor%error_bound
            changed=accepted>0
            if(changed)axis=editor%axis
        end if
        call editor%dialog%delete()
    end function edit_single_patch_degree

    subroutine add_degree_summary(layout,editor)
        type(QVBoxLayout), intent(inout) :: layout
        type(degree_editor), intent(inout) :: editor
        type(QFrame) :: summary
        type(QGridLayout) :: grid
        type(QLabel) :: label

        summary=QFrame(editor%dialog)
        call summary%setObjectName('studioDegreeSummary')
        call summary%setStyleSheet('QFrame#studioDegreeSummary { background:#16212a; '// &
            'border:1px solid #344755; }')
        grid=QGridLayout(summary)
        call grid%setContentsMargins(10,6,10,6)
        call grid%setSpacing(3)
        call grid%setColumnStretch(0,1)
        call grid%setColumnStretch(2,1)
        label=QLabel('Current degree',summary)
        call label%setProperty('dialogRole','detail')
        call grid%addWidget(label,0,0,alignment=Qt_AlignHCenter)
        label=QLabel('Requested degree',summary)
        call label%setProperty('dialogRole','detail')
        call grid%addWidget(label,0,2,alignment=Qt_AlignHCenter)
        editor%current_degree=QLabel('',summary)
        call editor%current_degree%setObjectName('studioDegreeCurrent')
        call editor%current_degree%setStyleSheet('QLabel { font-size:20px; font-weight:650; }')
        call grid%addWidget(editor%current_degree,1,0,alignment=Qt_AlignHCenter)
        label=QLabel('→',summary)
        call label%setStyleSheet('QLabel { color:#a9bdca; font-size:18px; }')
        call grid%addWidget(label,1,1,alignment=Qt_AlignHCenter)
        editor%target_degree=QLabel('',summary)
        call editor%target_degree%setObjectName('studioDegreeTarget')
        call editor%target_degree%setStyleSheet('QLabel { color:#7ee6dc; font-size:20px; font-weight:650; }')
        call grid%addWidget(editor%target_degree,1,2,alignment=Qt_AlignHCenter)
        call layout%addWidget(summary)
    end subroutine add_degree_summary

    function integer_text(value) result(text)
        integer, intent(in) :: value
        character(len=:), allocatable :: text
        character(len=20) :: buffer
        write(buffer,'(i0)') value
        text=trim(buffer)
    end function integer_text

    logical function parse_tolerance(editor,tolerance) result(valid)
        type(degree_editor), intent(in) :: editor
        real(c_double), intent(out) :: tolerance
        character(len=:), allocatable :: text
        integer :: comma,status

        valid=.false.
        tolerance=0._c_double
        text=trim(adjustl(editor%tolerance%text()))
        if (editor%decimal%validate(text)/=QValidator_Acceptable) return
        comma=index(text,',')
        if (comma>0) text(comma:comma)='.'
        read(text,*,iostat=status,decimal='point') tolerance
        if (status/=0.or..not.ieee_is_finite(tolerance)) return
        valid=tolerance>=0._c_double
    end function parse_tolerance

    subroutine direction_changed(index,data)
        integer, intent(in) :: index
        class(*), intent(inout), optional :: data
        integer :: axis,limit

        select type(editor=>data)
        type is(degree_editor)
            axis=index+1
            if(axis<1.or.axis>editor%kind.or.axis==editor%axis)return
            editor%axis=axis
            limit=editor%limits(axis)
            call editor%current_degree%setText(integer_text(editor%degrees(axis)))
            call editor%amount%setRange(min(1,limit),limit)
            call editor%amount%setValue(min(1,limit))
            call editor%amount%setEnabled(limit>0)
            if(editor%reducing)call editor%tolerance%setEnabled(limit>0)
            call editor%apply%setEnabled(limit>0)
            call editor%apply%setDefault(limit>0)
            call amount_changed(editor%amount%value(),editor)
            if(limit==0)then
                if(editor%reducing)then
                    call show_error(editor,'This direction is already at degree zero.')
                else
                    call show_error(editor,'Elevation would exceed the supported index range in this direction.')
                end if
            end if
        end select
    end subroutine direction_changed

    subroutine amount_changed(amount,data)
        integer, intent(in) :: amount
        class(*), intent(inout), optional :: data
        integer :: degree

        select type(editor=>data)
        type is(degree_editor)
            degree=editor%degrees(editor%axis)+merge(-amount,amount,editor%reducing)
            call editor%target_degree%setText(integer_text(degree))
            call publish_window_status(editor%dialog,'',blocking=.true.)
        end select
    end subroutine amount_changed

    subroutine tolerance_changed(text,data)
        character(len=*), intent(in) :: text
        class(*), intent(inout), optional :: data

        select type(editor=>data)
        type is(degree_editor)
            call publish_window_status(editor%dialog,'',blocking=.true.)
        end select
    end subroutine tolerance_changed

    subroutine show_error(editor,message)
        type(degree_editor), intent(inout) :: editor
        character(len=*), intent(in) :: message

        call publish_window_status(editor%dialog,message,warning=.true.,blocking=.true.)
    end subroutine show_error

    subroutine apply_degree(data)
        class(*), intent(inout), optional :: data
        character(len=:), allocatable :: error
        real(c_double) :: tolerance
        integer :: steps

        select type(editor=>data)
        type is(degree_editor)
            steps=editor%amount%value()
            if (steps<=0.or.steps>editor%limits(editor%axis)) return
            tolerance=0._c_double
            if (editor%reducing) then
                if (.not.parse_tolerance(editor,tolerance)) then
                    call show_error(editor,'Enter a finite, nonnegative error tolerance.')
                    call editor%tolerance%setFocus()
                    return
                end if
            end if
            call commit_single_patch_degree(editor%kind,editor%axis,steps,editor%reducing,tolerance, &
                editor%accepted,editor%error_bound,error,patch_id=editor%patch_id,catalog=editor%catalog)
            if (len(error)>0) then
                call show_error(editor,error)
                return
            end if
            if (editor%accepted==0) then
                call show_error(editor,'ForCAD could not apply this degree change.')
                return
            end if
            editor%requested=steps
            call editor%dialog%accept()
        end select
    end subroutine apply_degree

    subroutine cancel_degree(data)
        class(*), intent(inout), optional :: data

        select type(editor=>data)
        type is(degree_editor)
            call editor%dialog%reject()
        end select
    end subroutine cancel_degree

end module studio_degree_editor
