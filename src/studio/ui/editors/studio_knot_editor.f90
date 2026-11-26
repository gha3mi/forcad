! Native knot editing for the current ForCAD single-patch direction.
module studio_knot_editor
    use studio_status_hints, only: set_status_hint, publish_window_status
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use qt, only: QButtonGroup, QComboBox, QDialog, QDialog_Accepted, QDialogButtonBox, QDialogButtonBox_AcceptRole, &
        QDialogButtonBox_Cancel, QGridLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QRegularExpressionValidator, &
        QSpinBox, Qt_PlainText, QValidator_Acceptable, QVBoxLayout, QWidget
    use studio_single_patch, only: studio_patch_catalog, single_patch_info, read_single_patch, commit_single_patch_knot
    use studio_theme, only: studio_style_primary
    use studio_parameter_editor, only: parameter_direction_buttons, check_parameter_direction
    implicit none(type, external)
    private
    public :: edit_single_patch_knots

    character(len=3), parameter :: directions='UVW'

    type :: knot_editor
        type(QDialog) :: dialog
        type(QLineEdit) :: value_input
        type(QButtonGroup) :: direction_input
        type(QComboBox) :: group_input
        type(QSpinBox) :: count_input
        type(QPushButton) :: apply
        type(QLabel) :: interval, context
        type(QRegularExpressionValidator) :: decimal
        type(single_patch_info) :: patch
        type(studio_patch_catalog), pointer :: catalog=>null()
        real(c_double), allocatable :: groups(:)
        real(c_double) :: probe(3), edited_value=0._c_double
        integer, allocatable :: multiplicities(:)
        integer :: kind, axis, patch_id, accepted=0, requested=0, group_count=0
        logical :: removing, loading=.false.
    end type knot_editor

contains

    logical function edit_single_patch_knots(parent,kind,axis,probe,removing,accepted,requested,selected_value,edited_value, &
        catalog,patch_id) result(changed)
        class(QWidget), intent(in) :: parent
        integer, intent(in) :: kind,patch_id
        integer, intent(inout) :: axis
        real(c_double), intent(in) :: probe(3)
        logical, intent(in) :: removing
        integer, intent(out) :: accepted,requested
        real(c_double), optional, intent(in) :: selected_value
        real(c_double), optional, intent(out) :: edited_value
        type(studio_patch_catalog), intent(inout), target :: catalog
        type(knot_editor), target :: editor
        type(QVBoxLayout) :: layout
        type(QHBoxLayout) :: direction_row
        type(QGridLayout) :: form
        type(QLabel) :: heading,label
        type(QDialogButtonBox) :: buttons
        type(QPushButton) :: cancel
        character(len=:), allocatable :: action,copies_caption
        real(c_double) :: initial
        type(QWidget) :: direction_choices

        changed=.false.
        accepted=0
        requested=0
        if (present(edited_value)) edited_value=0._c_double
        if (kind<1.or.kind>3.or.axis<1.or.axis>kind) return
        ! Knot metadata is populated before ForCAD evaluates the probe and its derivatives.
        call read_single_patch(kind,probe,axis,editor%patch,patch_id=patch_id,catalog=catalog)
        if (.not.allocated(editor%patch%knots)) return
        if (editor%patch%dimension/=kind.or.size(editor%patch%knots)==0) return
        editor%kind=kind
        editor%axis=axis
        editor%patch_id=patch_id
        editor%catalog=>catalog
        editor%removing=removing
        editor%probe=probe
        initial=editor%patch%parameters(axis)
        if (present(selected_value)) then
            if (interior_knot(editor,selected_value)) initial=selected_value
        end if
        if (removing) then
            action='Remove'
            copies_caption='Copies to remove'
        else
            action='Insert'
            copies_caption='Copies to insert'
        end if

        editor%dialog=QDialog(parent)
        call editor%dialog%setObjectName('studioKnotEditor')
        call editor%dialog%setWindowTitle(action//' knot')
        call editor%dialog%setMinimumWidth(472)
        layout=QVBoxLayout(editor%dialog)
        call layout%setContentsMargins(18,16,18,16)
        call layout%setSpacing(10)

        heading=QLabel(action//' knot',editor%dialog)
        call heading%setObjectName('studioDialogTitle')
        call layout%addWidget(heading)
        editor%interval=QLabel('',editor%dialog)
        call editor%interval%setObjectName('studioKnotRange')
        call editor%interval%setTextFormat(Qt_PlainText)
        call editor%interval%setWordWrap(.true.)
        call editor%interval%setProperty('dialogRole','detail')
        call layout%addWidget(editor%interval)

        direction_row=QHBoxLayout()
        call direction_row%setSpacing(10)
        label=QLabel('Direction',editor%dialog)
        call direction_row%addWidget(label)
        direction_choices=parameter_direction_buttons(editor%dialog,kind,axis,'studioKnotDirection', &
            editor%direction_input,direction_changed,editor)
        call direction_row%addWidget(direction_choices,1)
        call layout%addLayout(direction_row)

        form=QGridLayout()
        call form%setSpacing(7)
        call form%setColumnStretch(0,1)
        label=QLabel('Knot value',editor%dialog)
        call form%addWidget(label,0,0)
        label=QLabel('Copies',editor%dialog)
        call form%addWidget(label,0,1)
        if (removing) then
            editor%group_input=QComboBox(editor%dialog)
            call editor%group_input%setObjectName('studioKnotGroup')
            call editor%group_input%setAccessibleName('Interior knot value')
            call editor%group_input%setFixedHeight(32)
            call form%addWidget(editor%group_input,1,0)
        else
            editor%value_input=QLineEdit('',editor%dialog)
            call editor%value_input%setObjectName('studioKnotValue')
            call editor%value_input%setAccessibleName('Knot value in native parameter units')
            call editor%value_input%setAccessibleDescription('Enter a value strictly inside the active range.')
            call editor%value_input%setFixedHeight(32)
            call editor%value_input%setMaxLength(128)
            editor%decimal=QRegularExpressionValidator( &
                '[+-]?(?:[0-9]+(?:[.,][0-9]*)?|[.,][0-9]+)(?:[eE][+-]?[0-9]+)?',editor%dialog)
            call form%addWidget(editor%value_input,1,0)
        end if
        editor%count_input=QSpinBox(editor%dialog)
        call editor%count_input%setObjectName('studioKnotCount')
        call editor%count_input%setAccessibleName(copies_caption)
        call editor%count_input%setFixedHeight(32)
        call editor%count_input%setMinimumWidth(102)
        call editor%count_input%setRange(1,1)
        call form%addWidget(editor%count_input,1,1)
        call layout%addLayout(form)

        editor%context=QLabel('',editor%dialog)
        call editor%context%setObjectName('studioKnotContext')
        call editor%context%setWordWrap(.true.)
        call editor%context%setProperty('dialogRole','detail')
        call layout%addWidget(editor%context)

        buttons=QDialogButtonBox(QDialogButtonBox_Cancel,editor%dialog)
        editor%apply=buttons%addButton(action,QDialogButtonBox_AcceptRole)
        call editor%apply%setObjectName('studioKnotApply')
        if(removing)then
            call editor%apply%setAccessibleDescription('Remove copies that preserve the patch shape.')
        else
            call editor%apply%setAccessibleDescription('Insert knots while preserving the patch shape.')
        end if
        call studio_style_primary(editor%apply)
        call editor%apply%setMinimumWidth(80)
        call editor%apply%setDefault(.true.)
        cancel=buttons%button(QDialogButtonBox_Cancel)
        call cancel%setObjectName('studioKnotCancel')
        call cancel%setMinimumWidth(80)
        call layout%addWidget(buttons)
        call buttons%accepted%connect(apply_edit,editor)
        call buttons%rejected%connect(cancel_edit,editor)
        if (removing) then
            call editor%group_input%currentIndexChanged%connect(group_changed,editor)
        else
            call editor%value_input%textChanged%connect(value_changed,editor)
        end if
        call editor%count_input%valueChanged%connect(copies_changed,editor)
        call refresh_direction(editor,initial)
        if (removing) then
            if (editor%group_count>0) call editor%group_input%setFocus()
        else
            call editor%value_input%setFocus()
        end if

        if (editor%dialog%exec()==QDialog_Accepted) then
            accepted=editor%accepted
            requested=editor%requested
            if (present(edited_value)) edited_value=editor%edited_value
            changed=accepted>0
            if (changed) axis=editor%axis
        end if
        call editor%dialog%delete()
    end function edit_single_patch_knots

    subroutine refresh_direction(editor,initial)
        type(knot_editor), intent(inout) :: editor
        real(c_double), intent(in) :: initial
        real(c_double) :: lower,upper,value
        integer :: axis

        editor%loading=.true.
        axis=editor%axis
        lower=editor%patch%bounds(1,axis)
        upper=editor%patch%bounds(2,axis)
        call editor%interval%setText(editor%patch%name//' · '//directions(axis:axis)//' range ['// &
            knot_text(lower)//', '//knot_text(upper)//']')
        call editor%count_input%setValue(1)
        if (editor%removing) then
            call load_interior_groups(editor,initial)
        else
            value=initial
            if (.not.interior_knot(editor,value)) value=0.5_c_double*lower+0.5_c_double*upper
            call editor%value_input%setText(knot_text(value))
        end if
        editor%loading=.false.
        call update_context(editor)
    end subroutine refresh_direction

    recursive subroutine direction_changed(index_,data)
        integer, intent(in) :: index_
        class(*), intent(inout), optional :: data
        type(single_patch_info) :: patch
        integer :: axis

        select type(editor=>data)
        type is(knot_editor)
            if (editor%loading) return
            axis=index_+1
            if (axis<1.or.axis>editor%kind.or.axis==editor%axis) return
            call read_single_patch(editor%kind,editor%probe,axis,patch,patch_id=editor%patch_id,catalog=editor%catalog)
            if (allocated(patch%knots)) then
                if (patch%dimension==editor%kind.and.size(patch%knots)>0) then
                    editor%axis=axis
                    editor%patch=patch
                    call refresh_direction(editor,patch%parameters(axis))
                    return
                end if
            end if
            call check_parameter_direction(editor%direction_input,editor%axis)
            call show_error(editor,'Could not read the knots in the selected direction.')
        end select
    end subroutine direction_changed

    logical function interior_knot(editor,knot) result(interior)
        type(knot_editor), intent(in) :: editor
        real(c_double), intent(in) :: knot
        interior=knot>editor%patch%bounds(1,editor%axis).and.knot<editor%patch%bounds(2,editor%axis)
    end function interior_knot

    subroutine load_interior_groups(editor,initial)
        type(knot_editor), intent(inout) :: editor
        real(c_double), intent(in) :: initial
        integer :: first,last,total,multiplicity
        real(c_double) :: knot

        total=size(editor%patch%knots)
        editor%group_count=0
        if (allocated(editor%groups)) deallocate(editor%groups,editor%multiplicities)
        allocate(editor%groups(total),editor%multiplicities(total))
        call editor%group_input%clear()
        call editor%group_input%setAccessibleDescription('')
        first=1
        do while (first<=total)
            knot=editor%patch%knots(first)
            last=first
            do while (last<total)
                if (editor%patch%knots(last+1)/=knot) exit
                last=last+1
            end do
            if (interior_knot(editor,knot)) then
                editor%group_count=editor%group_count+1
                multiplicity=last-first+1
                editor%groups(editor%group_count)=knot
                editor%multiplicities(editor%group_count)=multiplicity
                call editor%group_input%addItem(knot_text(knot))
            end if
            first=last+1
        end do
        call editor%group_input%setEnabled(editor%group_count>0)
        if (editor%group_count>0) then
            call editor%group_input%setCurrentIndex(minloc(abs(editor%groups(:editor%group_count)-initial),dim=1)-1)
        else
            call editor%group_input%addItem('No interior knots')
        end if
    end subroutine load_interior_groups

    function knot_text(value) result(text)
        real(c_double), intent(in) :: value
        character(len=:), allocatable :: text
        character(len=48) :: buffer
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
    end function knot_text

    function integer_text(value) result(text)
        integer, intent(in) :: value
        character(len=:), allocatable :: text
        character(len=20) :: buffer
        write(buffer,'(i0)') value
        text=trim(buffer)
    end function integer_text

    logical function read_value(editor,knot) result(valid)
        type(knot_editor), intent(in) :: editor
        real(c_double), intent(out) :: knot
        character(len=:), allocatable :: text
        integer :: status,comma

        knot=0._c_double
        valid=.false.
        text=trim(adjustl(editor%value_input%text()))
        if (editor%decimal%validate(text)/=QValidator_Acceptable) return
        comma=index(text,',')
        if (comma>0) text(comma:comma)='.'
        read(text,*,iostat=status,decimal='point') knot
        if (status==0) valid=ieee_is_finite(knot)
    end function read_value

    subroutine update_context(editor)
        type(knot_editor), intent(inout) :: editor
        integer :: index_,multiplicity,available
        real(c_double) :: knot
        character(len=:), allocatable :: message

        if (editor%loading) return
        call publish_window_status(editor%dialog,'',blocking=.true.)
        if (editor%removing) then
            index_=editor%group_input%currentIndex()+1
            if (index_<1.or.index_>editor%group_count) then
                call editor%count_input%setEnabled(.false.)
                call editor%apply%setEnabled(.false.)
                call editor%context%setText('No interior knot groups')
                call set_status_hint(editor%dialog,'Insert a knot in this direction to make an interior group available.')
                return
            end if
            multiplicity=editor%multiplicities(index_)
            available=multiplicity
            call editor%group_input%setAccessibleDescription('Current multiplicity '//integer_text(multiplicity))
            message='Current multiplicity '//integer_text(multiplicity)
            call set_status_hint(editor%dialog,'Remove copies that preserve the patch shape.')
        else
            available=editor%count_input%maximum()
            message='Multiplicity unavailable'
            call set_status_hint(editor%dialog,'Choose a finite value strictly inside the active interval.')
            if (read_value(editor,knot)) then
                if(interior_knot(editor,knot))then
                    multiplicity=count(editor%patch%knots==knot)
                    available=editor%patch%degree(editor%axis)+1-multiplicity
                    message='Current multiplicity '//integer_text(multiplicity)//' · Limit '//integer_text(editor%patch%degree(editor%axis)+1)
                    if(available<=0)call set_status_hint(editor%dialog,'Insertion limit reached for this knot value.')
                end if
            end if
        end if
        call editor%count_input%setEnabled(available>0)
        call editor%apply%setEnabled(available>0)
        call editor%count_input%setRange(1,max(1,available))
        call editor%context%setText(message)
    end subroutine update_context

    subroutine group_changed(index_,data)
        integer, intent(in) :: index_
        class(*), intent(inout), optional :: data
        select type(editor=>data)
        type is(knot_editor)
            call update_context(editor)
        end select
    end subroutine group_changed

    subroutine copies_changed(value,data)
        integer, intent(in) :: value
        class(*), intent(inout), optional :: data
        select type(editor=>data)
        type is(knot_editor)
            if (.not.editor%loading) call publish_window_status(editor%dialog,'',blocking=.true.)
        end select
    end subroutine copies_changed

    subroutine value_changed(text,data)
        character(len=*), intent(in) :: text
        class(*), intent(inout), optional :: data
        select type(editor=>data)
        type is(knot_editor)
            call update_context(editor)
        end select
    end subroutine value_changed

    subroutine show_error(editor,message)
        type(knot_editor), intent(inout) :: editor
        character(len=*), intent(in) :: message
        call publish_window_status(editor%dialog,message,warning=.true.,blocking=.true.)
    end subroutine show_error

    subroutine apply_edit(data)
        class(*), intent(inout), optional :: data
        character(len=:), allocatable :: error
        real(c_double) :: knot
        integer :: index_,accepted,requested

        select type(editor=>data)
        type is(knot_editor)
            if (.not.editor%apply%isEnabled()) return
            if (editor%removing) then
                index_=editor%group_input%currentIndex()+1
                if (index_<1.or.index_>editor%group_count) return
                knot=editor%groups(index_)
            else
                if (.not.read_value(editor,knot)) then
                    call show_error(editor,'Enter a finite numeric knot value.')
                    call editor%value_input%setFocus()
                    return
                end if
                if (.not.interior_knot(editor,knot)) then
                    call show_error(editor,'Enter a value strictly inside the active parameter interval.')
                    call editor%value_input%setFocus()
                    return
                end if
            end if
            requested=editor%count_input%value()
            call commit_single_patch_knot(editor%kind,editor%axis,knot,requested,editor%removing,accepted,error, &
                patch_id=editor%patch_id,catalog=editor%catalog)
            if (len(error)>0) then
                call show_error(editor,error)
                return
            end if
            if (accepted==0) then
                call show_error(editor,'ForCAD could not apply any of the requested copies.')
                return
            end if
            editor%accepted=accepted
            editor%requested=requested
            editor%edited_value=knot
            call editor%dialog%accept()
        end select
    end subroutine apply_edit

    subroutine cancel_edit(data)
        class(*), intent(inout), optional :: data
        select type(editor=>data)
        type is(knot_editor)
            call editor%dialog%reject()
        end select
    end subroutine cancel_edit

end module studio_knot_editor
