module studio_parameter_editor
    use studio_status_hints, only: set_status_hint, publish_window_status
    use studio_presentation_text, only: number_text
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use qt, only: QAbstractButton, QAbstractItemView_AnyKeyPressed, QAbstractItemView_DoubleClicked, &
        QAbstractItemView_EditKeyPressed, QAbstractItemView_SelectItems, QAbstractItemView_SingleSelection, &
        QButtonGroup, QColor, QDialog, QDialog_Accepted, QDialogButtonBox, QDialogButtonBox_AcceptRole, &
        QDialogButtonBox_Cancel, QFormLayout, QHBoxLayout, QHeaderView, QHeaderView_Stretch, QHeaderView_ResizeToContents, &
        QLabel, QPushButton, QRegularExpressionValidator, Qt_AlignCenter, Qt_ItemIsEditable, Qt_ItemIsEnabled, &
        Qt_ItemIsSelectable, Qt_ScrollBarAlwaysOff, QTableWidget, QTableWidgetItem, QValidator_Acceptable, &
        QVBoxLayout, QWidget, qt_slot_int
    use studio_single_patch, only: studio_patch_catalog, single_patch_info, read_single_patch, &
        commit_single_patch_knot_values
    use studio_theme, only: studio_active_edge, studio_selected_fill, studio_style_primary
    implicit none(type, external)
    private
    public :: edit_single_patch_knot_values
    public :: parse_parameter_number, number_text, parameter_number_pattern
    public :: parameter_direction_buttons, check_parameter_direction

    character(len=3), parameter :: directions='UVW'
    integer, parameter :: value_column=1
    character(len=*), parameter :: parameter_number_pattern= &
        '[+-]?(?:[0-9]+(?:[.,][0-9]*)?|[.,][0-9]+)(?:[eE][+-]?[0-9]+)?'

    type :: parameter_editor
        type(QDialog) :: dialog
        type(QPushButton) :: apply
        type(QLabel) :: context
        type(QButtonGroup) :: direction_input
        type(QTableWidget) :: table
        type(QRegularExpressionValidator) :: decimal
        type(single_patch_info) :: patch
        type(studio_patch_catalog), pointer :: catalog=>null()
        integer, allocatable :: group_start(:)
        integer :: kind,axis,patch_id
        logical :: changed=.false., loading=.false.
    end type parameter_editor

contains

    function parameter_direction_buttons(parent,kind,axis,object_name,group,changed,data,include_all) result(selector)
        class(QWidget), intent(in) :: parent
        integer, intent(in) :: kind,axis
        character(*), intent(in) :: object_name
        type(QButtonGroup), intent(out) :: group
        procedure(qt_slot_int) :: changed
        class(*), intent(inout), target :: data
        logical, optional, intent(in) :: include_all
        type(QWidget) :: selector
        type(QHBoxLayout) :: row
        type(QPushButton) :: button
        character(:), allocatable :: label,accessible_name
        integer :: i,choice_count

        selector=QWidget(parent)
        call selector%setObjectName(object_name)
        call selector%setAccessibleName('Parameter direction')
        row=QHBoxLayout(selector)
        call row%setContentsMargins(0,0,0,0)
        call row%setSpacing(4)
        group=QButtonGroup(selector)
        call group%setExclusive(.true.)
        choice_count=kind
        if(present(include_all))then
            if(include_all.and.kind>1)choice_count=kind+1
        end if
        do i=1,choice_count
            if(i<=kind)then
                label=directions(i:i)
                accessible_name=label//' parameter direction'
            else
                label='All'
                accessible_name='All parameter directions'
            end if
            button=QPushButton(label,selector)
            call button%setObjectName(object_name//label)
            call button%setAccessibleName(accessible_name)
            call button%setProperty('autoDefault',.false.)
            call button%setCheckable(.true.)
            call button%setChecked(i==axis)
            call button%setEnabled(kind>1)
            call button%setMinimumWidth(36)
            call row%addWidget(button)
            call group%addButton(button,i-1)
        end do
        call row%addStretch(1)
        call group%idClicked%connect(changed,data)
    end function parameter_direction_buttons

    subroutine check_parameter_direction(group,axis)
        type(QButtonGroup), intent(in) :: group
        integer, intent(in) :: axis
        type(QAbstractButton) :: button

        button=group%button(axis-1)
        call button%setChecked(.true.)
    end subroutine check_parameter_direction

    logical function edit_single_patch_knot_values(parent,kind,axis,catalog,patch_id) result(changed)
        class(QWidget), intent(in) :: parent
        integer, intent(in) :: kind,patch_id
        integer, intent(inout) :: axis
        type(studio_patch_catalog), intent(inout), target :: catalog
        type(parameter_editor), target :: editor
        type(QVBoxLayout) :: layout
        type(QFormLayout) :: form
        type(QLabel) :: heading
        type(QDialogButtonBox) :: buttons
        type(QPushButton) :: cancel
        type(QWidget) :: direction_choices
        character(*), parameter :: prefix='studioKnotValues'

        changed=.false.
        if (kind<1.or.kind>3.or.axis<1.or.axis>kind) return
        ! The adapter sets metadata before evaluating the probe and its derivatives.
        call read_single_patch(kind,[.5_c_double,.5_c_double,.5_c_double],axis,editor%patch, &
            patch_id=patch_id,catalog=catalog)
        if (.not.allocated(editor%patch%knots)) return
        if (editor%patch%dimension/=kind.or.size(editor%patch%knots)<2) return
        editor%kind=kind
        editor%axis=axis
        editor%patch_id=patch_id
        editor%catalog=>catalog
        editor%dialog=QDialog(parent)
        call editor%dialog%setObjectName(prefix//'Editor')
        call editor%dialog%setWindowTitle('Knot values')
        call editor%dialog%setMinimumSize(510,330)
        call editor%dialog%resize(600,420)
        call editor%dialog%setStyleSheet( &
            'QTableWidget#studioKnotValuesTable { background:#16212a; alternate-background-color:#1b2933;'// &
            'selection-background-color:'//studio_selected_fill//'; selection-color:#f5fafc; gridline-color:#334653; }'// &
            'QTableWidget#studioKnotValuesTable::item { padding:4px 8px; }'// &
            'QTableWidget#studioKnotValuesTable QLineEdit { background:#102a32; border:1px solid '//studio_active_edge//';'// &
            'padding:2px 6px; selection-background-color:#386d7a; }')
        layout=QVBoxLayout(editor%dialog)
        call layout%setContentsMargins(18,16,18,16)
        call layout%setSpacing(10)

        heading=QLabel('Knot values',editor%dialog)
        call heading%setObjectName('studioDialogTitle')
        call layout%addWidget(heading)
        editor%context=QLabel('',editor%dialog)
        call editor%context%setProperty('textFormat',0)
        call editor%context%setWordWrap(.true.)
        call editor%context%setProperty('dialogRole','detail')
        call layout%addWidget(editor%context)

        form=QFormLayout()
        direction_choices=parameter_direction_buttons(editor%dialog,kind,axis,prefix//'Direction', &
            editor%direction_input,direction_changed,editor)
        call form%addRow('Direction',direction_choices)
        call layout%addLayout(form)
        call set_status_hint(editor%dialog,'Edit an interior value to move all its copies. This can change the patch shape.')
        call add_knot_table(editor,layout)

        buttons=QDialogButtonBox(QDialogButtonBox_Cancel,editor%dialog)
        editor%apply=buttons%addButton('Apply',QDialogButtonBox_AcceptRole)
        call editor%apply%setObjectName(prefix//'Apply')
        call studio_style_primary(editor%apply)
        call editor%apply%setMinimumWidth(80)
        call editor%apply%setDefault(.true.)
        cancel=buttons%button(QDialogButtonBox_Cancel)
        call cancel%setObjectName(prefix//'Cancel')
        call cancel%setMinimumWidth(80)
        call layout%addWidget(buttons)

        editor%decimal=QRegularExpressionValidator(parameter_number_pattern,editor%dialog)
        call buttons%accepted%connect(apply_edit,editor)
        call buttons%rejected%connect(cancel_edit,editor)
        call editor%table%cellChanged%connect(knot_draft_changed,editor)
        call refresh_direction(editor)
        if (editor%apply%isEnabled()) call editor%table%setFocus()
        if (editor%dialog%exec()==QDialog_Accepted) then
            changed=editor%changed
            if (changed) axis=editor%axis
        end if
        call editor%dialog%delete()
    end function edit_single_patch_knot_values

    subroutine add_knot_table(editor,layout)
        type(parameter_editor), intent(inout) :: editor
        type(QVBoxLayout), intent(inout) :: layout
        type(QHeaderView) :: header

        editor%table=QTableWidget(0,3,editor%dialog)
        call editor%table%setObjectName('studioKnotValuesTable')
        call editor%table%setProperty('accessibleName','Knot group values')
        call set_status_hint(editor%table,'Double-click or type to edit an interior value. Endpoint and extension rows are locked.')
        call editor%table%setVerticalHeaderVisible(.false.)
        call editor%table%setWordWrap(.false.)
        call editor%table%setAlternatingRowColors(.true.)
        call editor%table%setSelectionMode(QAbstractItemView_SingleSelection)
        call editor%table%setSelectionBehavior(QAbstractItemView_SelectItems)
        call editor%table%setEditTriggers(ior(QAbstractItemView_DoubleClicked, &
            ior(QAbstractItemView_EditKeyPressed,QAbstractItemView_AnyKeyPressed)))
        call editor%table%setProperty('horizontalScrollBarPolicy',Qt_ScrollBarAlwaysOff)
        call editor%table%setHorizontalHeaderLabels([character(len=16)::'Indices','Parameter value','Multiplicity'])
        header=editor%table%horizontalHeader()
        call header%setSectionResizeMode(0,QHeaderView_ResizeToContents)
        call header%setSectionResizeMode(1,QHeaderView_Stretch)
        call header%setSectionResizeMode(2,QHeaderView_ResizeToContents)
        call layout%addWidget(editor%table,1)
    end subroutine add_knot_table

    subroutine refresh_direction(editor)
        type(parameter_editor), intent(inout) :: editor
        integer :: axis,n,i,row,first,last,first_editable
        real(c_double) :: value
        character(:), allocatable :: indices,text

        axis=editor%axis
        call editor%context%setText(editor%patch%name//' · '//directions(axis:axis)//' direction · Active interval ['// &
            number_text(editor%patch%bounds(1,axis))//', '//number_text(editor%patch%bounds(2,axis))//']')
        n=size(editor%patch%knots)
        editor%group_start=[1,pack([(i,i=2,n)],editor%patch%knots(2:n)/=editor%patch%knots(:n-1)),n+1]
        call editor%table%setRowCount(0)
        call editor%table%setRowCount(size(editor%group_start)-1)
        first_editable=-1
        do row=1,size(editor%group_start)-1
            first=editor%group_start(row)
            last=editor%group_start(row+1)-1
            value=editor%patch%knots(first)
            indices=integer_text(first-1)
            if (last>first) indices=indices//'–'//integer_text(last-1)
            text=number_text(value)
            call set_cell(editor%table,row-1,0,indices,.false.,'Zero-based knot indices')
            call set_cell(editor%table,row-1,2,integer_text(last-first+1),.false.,'Copies at this value')
            if (interior_knot(editor,value)) then
                call set_cell(editor%table,row-1,value_column,text,.true.,text)
                if (first_editable<0) first_editable=row-1
            else
                call set_cell(editor%table,row-1,value_column,text,.false.,'Endpoint or extension value · '//text)
            end if
        end do
        call editor%table%resizeRowsToContents()
        call publish_window_status(editor%dialog,'',blocking=.true.)
        call editor%apply%setEnabled(first_editable>=0)
        if (first_editable>=0) then
            call editor%table%setCurrentCell(first_editable,value_column)
        else
            call show_error(editor,'Insert an interior knot in this direction to edit its value.')
        end if
    end subroutine refresh_direction

    recursive subroutine direction_changed(index_,data)
        integer, intent(in) :: index_
        class(*), intent(inout), optional :: data
        type(single_patch_info) :: patch
        integer :: axis

        select type(editor=>data)
        type is(parameter_editor)
            if (editor%loading) return
            axis=index_+1
            if (axis<1.or.axis>editor%kind.or.axis==editor%axis) return
            call read_single_patch(editor%kind,[.5_c_double,.5_c_double,.5_c_double],axis,patch, &
                patch_id=editor%patch_id,catalog=editor%catalog)
            if (allocated(patch%knots)) then
                if (patch%dimension==editor%kind.and.size(patch%knots)>=2) then
                    editor%axis=axis
                    editor%patch=patch
                    call refresh_direction(editor)
                    return
                end if
            end if
            call check_parameter_direction(editor%direction_input,editor%axis)
            call show_error(editor,'Could not read the knots in the selected direction.')
        end select
    end subroutine direction_changed

    subroutine set_cell(table,row,column,text,editable,tooltip)
        type(QTableWidget), intent(inout) :: table
        integer, intent(in) :: row,column
        character(len=*), intent(in) :: text,tooltip
        logical, intent(in) :: editable
        type(QTableWidgetItem) :: item
        integer :: flags

        item=QTableWidgetItem(text)
        flags=ior(Qt_ItemIsEnabled,Qt_ItemIsSelectable)
        if (editable) flags=ior(flags,Qt_ItemIsEditable)
        call item%setFlags(flags)
        if (.not.editable) call item%setForeground(QColor(169,189,202))
        call item%setTextAlignment(Qt_AlignCenter)
        call set_status_hint(item,tooltip)
        call table%setItem(row,column,item)
    end subroutine set_cell

    logical function interior_knot(editor,value) result(interior)
        type(parameter_editor), intent(in) :: editor
        real(c_double), intent(in) :: value

        interior=value>editor%patch%bounds(1,editor%axis).and.value<editor%patch%bounds(2,editor%axis)
    end function interior_knot

    logical function parse_parameter_number(decimal,input,value) result(valid)
        type(QRegularExpressionValidator), intent(in) :: decimal
        character(len=*), intent(in) :: input
        real(c_double), intent(out) :: value
        character(:), allocatable :: text
        integer :: status,comma

        valid=.false.
        value=0._c_double
        text=trim(adjustl(input))
        if (len(text)==0.or.len(text)>128) return
        if (decimal%validate(text)/=QValidator_Acceptable) return
        comma=index(text,',')
        if (comma>0) text(comma:comma)='.'
        read(text,*,iostat=status,decimal='point') value
        if (status==0) valid=ieee_is_finite(value)
    end function parse_parameter_number

    logical function apply_knots(editor) result(accepted)
        type(parameter_editor), intent(inout) :: editor
        real(c_double), allocatable :: candidate(:)
        real(c_double) :: value
        type(QTableWidgetItem) :: item
        character(:), allocatable :: text,error
        integer :: row,first,last,previous,focus

        accepted=.false.
        candidate=editor%patch%knots
        do row=1,size(editor%group_start)-1
            first=editor%group_start(row)
            last=editor%group_start(row+1)-1
            if (.not.interior_knot(editor,candidate(first))) cycle
            item=editor%table%item(row-1,value_column)
            text=item%text()
            ! Use the exact stored value for every untouched cell.
            if (text==number_text(candidate(first))) cycle
            if (.not.parse_parameter_number(editor%decimal,text,value)) then
                call show_error(editor,'Enter a finite number for the highlighted knot value.')
                call focus_knot(editor,row-1)
                return
            end if
            if (.not.interior_knot(editor,value)) then
                call show_error(editor,'Keep interior knots strictly inside the active interval.')
                call focus_knot(editor,row-1)
                return
            end if
            candidate(first:last)=value
        end do
        do row=2,size(editor%group_start)-1
            previous=editor%group_start(row-1)
            first=editor%group_start(row)
            if (candidate(first)>candidate(previous)) cycle
            call show_error(editor,'Keep consecutive knot groups in strictly ascending order.')
            focus=row-1
            if (candidate(first)==editor%patch%knots(first)) focus=row-2
            call focus_knot(editor,focus)
            return
        end do
        call commit_single_patch_knot_values(editor%kind,editor%axis,candidate,editor%changed,error, &
            patch_id=editor%patch_id,catalog=editor%catalog)
        accepted=len(error)==0
        if (.not.accepted) call show_error(editor,error)
    end function apply_knots

    subroutine apply_edit(data)
        class(*), intent(inout), optional :: data

        select type(editor=>data)
        type is(parameter_editor)
            if (.not.editor%apply%isEnabled()) return
            ! Focus transfer commits an open native table cell before reading it.
            call editor%apply%setFocus()
            call publish_window_status(editor%dialog,'',blocking=.true.)
            if (.not.apply_knots(editor)) return
            call editor%dialog%accept()
        end select
    end subroutine apply_edit

    subroutine cancel_edit(data)
        class(*), intent(inout), optional :: data

        select type(editor=>data)
        type is(parameter_editor)
            call editor%dialog%reject()
        end select
    end subroutine cancel_edit

    subroutine show_error(editor,message)
        type(parameter_editor), intent(inout) :: editor
        character(len=*), intent(in) :: message

        call publish_window_status(editor%dialog,message,warning=.true.,blocking=.true.)
    end subroutine show_error

    subroutine focus_knot(editor,row)
        type(parameter_editor), intent(inout) :: editor
        integer, intent(in) :: row

        call editor%table%setCurrentCell(row,value_column)
        call editor%table%setFocus()
    end subroutine focus_knot

    subroutine knot_draft_changed(row,column,data)
        integer, intent(in) :: row,column
        class(*), intent(inout), optional :: data

        select type(editor=>data)
        type is(parameter_editor)
            if (column==value_column) call publish_window_status(editor%dialog,'',blocking=.true.)
        end select
    end subroutine knot_draft_changed

    function integer_text(value) result(text)
        integer, intent(in) :: value
        character(:), allocatable :: text
        character(len=20) :: buffer

        write(buffer,'(i0)') value
        text=trim(buffer)
    end function integer_text

end module studio_parameter_editor
