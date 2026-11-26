module studio_control_editor
    use studio_status_hints, only: set_status_hint, publish_window_status
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: iso_fortran_env, only: int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use qt, only: QAbstractItemView_AnyKeyPressed, QAbstractItemView_DoubleClicked, QAbstractItemView_NoEditTriggers, &
        QAbstractItemView_EditKeyPressed, QAbstractItemView_SelectItems, QAbstractItemView_SingleSelection, &
        QColor, QDialog, QDialog_Accepted, QDialogButtonBox, QDialogButtonBox_AcceptRole, &
        QDialogButtonBox_ActionRole, QDialogButtonBox_Cancel, QHBoxLayout, QHeaderView, QHeaderView_ResizeToContents, &
        QHeaderView_Stretch, QLabel, QPushButton, QRegularExpressionValidator, &
        Qt_AlignCenter, Qt_AlignRight, Qt_AlignVCenter, Qt_ItemIsEditable, Qt_ItemIsEnabled, &
        Qt_ItemIsSelectable, Qt_PlainText, QTableWidget, QTableWidgetItem, &
        QValidator_Acceptable, QVBoxLayout, QWidget
    use studio_single_patch, only: studio_patch_catalog, single_patch_controls, read_single_patch_controls, &
        commit_single_patch_controls, active_patch_id
    use studio_document, only: studio_document_t
    use studio_multipatch_controls, only: multipatch_controls_t, multipatch_control_preview_t, read_document_multipatch_controls, &
        preview_document_multipatch_controls, commit_document_multipatch_controls
    use studio_units, only: unitless,display_length,model_length,length_caption
    use studio_theme, only: studio_active_edge, studio_selected_fill, studio_style_primary
    implicit none(type, external)
    private
    public :: edit_single_patch_controls, edit_multipatch_controls, index_text

    integer, parameter :: weight_column=4
    character(len=6), parameter :: component_names(4)=[character(len=6)::'X','Y','Z','Weight']

    type :: control_editor
        type(QDialog) :: dialog
        type(QTableWidget) :: table, preview_table
        type(QPushButton) :: apply, preview
        type(QRegularExpressionValidator) :: decimal
        type(single_patch_controls) :: original
        type(studio_patch_catalog), pointer :: catalog=>null()
        type(studio_document_t), pointer :: document=>null()
        type(multipatch_controls_t) :: shared
        integer :: assembly_id=0,kind,variant,invalid_row=-1,invalid_column=-1
        integer :: model_unit=unitless,display_unit=unitless
        integer(int64) :: revision
        logical :: changed=.false., updating_feedback=.false.
    end type control_editor

contains

    logical function edit_single_patch_controls(parent,kind,selected_point,model_unit,display_unit,catalog,patch_id) result(changed)
        class(QWidget), intent(in) :: parent
        integer, intent(in) :: kind,patch_id
        integer, optional, intent(in) :: selected_point
        integer, optional, intent(in) :: model_unit,display_unit
        type(studio_patch_catalog), intent(inout), target :: catalog
        type(control_editor), target :: editor

        editor%kind=kind
        if(present(model_unit))editor%model_unit=model_unit
        if(present(display_unit))editor%display_unit=display_unit
        editor%catalog=>catalog
        editor%variant=patch_id
        editor%revision=catalog%patch_revision(kind,patch_id)
        call read_single_patch_controls(kind,editor%original,patch_id,catalog)
        changed=run_control_editor(parent,editor,selected_point)
    end function edit_single_patch_controls

    logical function edit_multipatch_controls(parent,document,id,selected_point,display_unit) result(changed)
        class(QWidget), intent(in) :: parent
        type(studio_document_t), intent(inout), target :: document
        integer, intent(in) :: id,selected_point,display_unit
        type(control_editor), target :: editor

        changed=.false.
        call read_document_multipatch_controls(document,id,editor%shared)
        if(.not.editor%shared%valid)return
        editor%document=>document
        editor%assembly_id=id
        editor%model_unit=document%model_unit
        editor%display_unit=display_unit
        editor%kind=1
        editor%original%valid=.true.
        editor%original%dimension=1
        editor%original%counts=[size(editor%shared%weights),1,1]
        editor%original%name=document%multipatches%groups(id)%definition%name
        editor%original%points=editor%shared%points
        editor%original%weights=editor%shared%weights
        changed=run_control_editor(parent,editor,selected_point)
    end function edit_multipatch_controls

    logical function run_control_editor(parent,editor,selected_point) result(changed)
        class(QWidget), intent(in) :: parent
        type(control_editor), intent(inout), target :: editor
        integer, optional, intent(in) :: selected_point
        type(QVBoxLayout) :: layout
        type(QHBoxLayout) :: header
        type(QLabel) :: title, count_label
        type(QDialogButtonBox) :: buttons
        type(QPushButton) :: cancel
        type(QHeaderView) :: header_view
        character(len=:), allocatable :: name
        integer :: count_rows,initial_row,height

        count_rows=0
        if(valid_controls(editor%original,editor%kind))count_rows=size(editor%original%weights)
        name='Control points'
        if (allocated(editor%original%name)) then
            if (len(editor%original%name)>0) name=editor%original%name
        end if

        editor%dialog=QDialog(parent)
        call editor%dialog%setObjectName('studioControlEditor')
        call editor%dialog%setWindowTitle('Points + weights')
        if(editor%assembly_id>0)call editor%dialog%setWindowTitle('Multipatch points + weights')
        call editor%dialog%setMinimumSize(510,294)
        call editor%dialog%setStyleSheet( &
            'QTableWidget#studioControlTable { background:#16212a; alternate-background-color:#1b2933;'// &
            'selection-background-color:'//studio_selected_fill//'; selection-color:#f5fafc; gridline-color:#334653; }'// &
            'QTableWidget#studioControlTable::item { padding:4px 8px; }'// &
            'QTableWidget#studioControlTable QLineEdit { background:#102a32; border:1px solid '//studio_active_edge//';'// &
            'padding:2px 6px; selection-background-color:#386d7a; }')
        layout=QVBoxLayout(editor%dialog)
        call layout%setContentsMargins(18,16,18,16)
        call layout%setSpacing(10)

        header=QHBoxLayout()
        title=QLabel(name,editor%dialog)
        call title%setObjectName('studioControlTitle')
        call title%setProperty('dialogRole','title')
        call title%setTextFormat(Qt_PlainText)
        call title%setWordWrap(.true.)
        call header%addWidget(title,1)
        if (count_rows>0) then
            count_label=QLabel(integer_text(count_rows)//' control points',editor%dialog)
            call count_label%setObjectName('studioControlCount')
            call count_label%setProperty('dialogRole','detail')
            call header%addWidget(count_label)
        end if
        call layout%addLayout(header)

        call set_status_hint(editor%dialog,'Edit X, Y, Z, or Weight. Weights must be positive.')
        if(editor%assembly_id>0)call set_status_hint(editor%dialog,'Shared controls update across their connected patches.')
        if(editor%shared%smooth_joins)call set_status_hint(editor%dialog,'Edit coordinates or weights, then preview the dependent controls.')

        editor%table=QTableWidget(count_rows,5,editor%dialog)
        call editor%table%setObjectName('studioControlTable')
        call editor%table%setAccessibleName('Control point coordinates and weights')
        call set_status_hint(editor%table,'Indices are one-based. Edit the X, Y, Z or Weight cells.')
        if(editor%assembly_id>0)call set_status_hint(editor%table,'Each ID is one control. The weight shown is from its first patch.')
        call editor%table%setVerticalHeaderVisible(.false.)
        call editor%table%setWordWrap(.false.)
        call editor%table%setAlternatingRowColors(.true.)
        call editor%table%setShowGrid(.false.)
        call editor%table%setSelectionMode(QAbstractItemView_SingleSelection)
        call editor%table%setSelectionBehavior(QAbstractItemView_SelectItems)
        call editor%table%setEditTriggers(ior(QAbstractItemView_DoubleClicked, &
            ior(QAbstractItemView_EditKeyPressed,QAbstractItemView_AnyKeyPressed)))
        call add_controls(editor,count_rows)
        header_view=editor%table%horizontalHeader()
        call header_view%setMinimumSectionSize(108)
        call editor%table%resizeColumnsToContents()
        height=min(540,max(304,184+29*min(count_rows,13)))
        if(editor%shared%smooth_joins)height=540
        call editor%dialog%resize(max(640,min(900,header_view%length()+58)),height)
        call editor%table%setStretchLastSection(.true.)
        call layout%addWidget(editor%table,3)
        if(editor%shared%smooth_joins)call add_preview_table(editor,layout)


        buttons=QDialogButtonBox(QDialogButtonBox_Cancel,editor%dialog)
        editor%apply=buttons%addButton('Apply',QDialogButtonBox_AcceptRole)
        call editor%apply%setObjectName('studioControlApply')
        call studio_style_primary(editor%apply)
        call editor%apply%setDefault(.true.)
        call editor%apply%setMinimumWidth(80)
        call editor%apply%setEnabled(count_rows>0.and..not.editor%shared%smooth_joins)
        if(editor%shared%smooth_joins)then
            editor%preview=buttons%addButton('Preview',QDialogButtonBox_ActionRole)
            call editor%preview%setObjectName('studioControlPreview')
            call editor%preview%setMinimumWidth(80)
            call editor%preview%setEnabled(count_rows>0)
            call editor%preview%setDefault(.true.)
            call set_status_hint(editor%preview,'Check the joins and list the control point and weight adjustments.')
            call editor%preview%clicked%connect(preview_controls,editor)
        end if
        cancel=buttons%button(QDialogButtonBox_Cancel)
        call cancel%setObjectName('studioControlCancel')
        call cancel%setMinimumWidth(80)
        call layout%addWidget(buttons)

        editor%decimal=QRegularExpressionValidator( &
            '[+-]?(?:[0-9]+(?:[.,][0-9]*)?|[.,][0-9]+)(?:[eE][+-]?[0-9]+)?',editor%dialog)
        call buttons%accepted%connect(apply_controls,editor)
        call buttons%rejected%connect(cancel_controls,editor)
        call editor%table%cellChanged%connect(control_cell_changed,editor)
        if (count_rows>0) then
            initial_row=0
            if(present(selected_point))initial_row=max(0,min(count_rows-1,selected_point-1))
            call editor%table%setCurrentCell(initial_row,1)
            call editor%table%setFocus()
        else
            call show_error(editor,'Control points are unavailable for this patch.')
        end if
        changed=.false.
        if (editor%dialog%exec()==QDialog_Accepted) changed=editor%changed
        call editor%dialog%delete()
    end function run_control_editor

    subroutine add_preview_table(editor,layout)
        type(control_editor), intent(inout) :: editor
        type(QVBoxLayout), intent(inout) :: layout
        type(QHeaderView) :: header


        editor%preview_table=QTableWidget(0,4,editor%dialog)
        call editor%preview_table%setObjectName('studioControlPreviewTable')
        call editor%preview_table%setAccessibleName('Dependent control point and weight adjustments')
        call editor%preview_table%setHorizontalHeaderLabels([character(12)::'Control ID','Field','Before','After'])
        call editor%preview_table%setMinimumHeight(112)
        call editor%preview_table%setVerticalHeaderVisible(.false.)
        call editor%preview_table%setAlternatingRowColors(.true.)
        call editor%preview_table%setWordWrap(.false.)
        call editor%preview_table%setEditTriggers(QAbstractItemView_NoEditTriggers)
        header=editor%preview_table%horizontalHeader()
        call header%setSectionResizeMode(0,QHeaderView_ResizeToContents)
        call header%setSectionResizeMode(1,QHeaderView_ResizeToContents)
        call header%setSectionResizeMode(2,QHeaderView_Stretch)
        call header%setSectionResizeMode(3,QHeaderView_Stretch)
        call layout%addWidget(editor%preview_table,2)
    end subroutine add_preview_table

    logical function valid_controls(controls,kind) result(valid)
        type(single_patch_controls), intent(in) :: controls
        integer, intent(in) :: kind
        valid=.false.
        if (.not.controls%valid.or.controls%dimension/=kind.or.kind<1.or.kind>3) return
        if (.not.allocated(controls%points).or..not.allocated(controls%weights)) return
        if (any(controls%counts<1).or.any(controls%counts(kind+1:)/=1)) return
        valid=size(controls%points,1)==3.and.size(controls%weights)==product(controls%counts).and. &
            size(controls%points,2)==size(controls%weights)
    end function valid_controls

    subroutine add_controls(editor,count_rows)
        type(control_editor), intent(inout) :: editor
        integer, intent(in) :: count_rows
        type(QTableWidgetItem) :: item
        character(len=16) :: index_header
        character(len=:), allocatable :: text
        integer :: row,column,read_flags,edit_flags

        select case(editor%kind)
        case(1); index_header='Index (U)'
        case(2); index_header='Index (U, V)'
        case default; index_header='Index (U, V, W)'
        end select
        if(editor%assembly_id>0)index_header='Control ID'
        call editor%table%setHorizontalHeaderLabels([character(len=16)::index_header, &
            (length_caption(component_names(column),editor%display_unit),column=1,3),'Weight'])
        read_flags=ior(Qt_ItemIsEnabled,Qt_ItemIsSelectable)
        edit_flags=ior(read_flags,Qt_ItemIsEditable)
        do row=1,count_rows
            item=QTableWidgetItem(index_text(row,editor%original%counts,editor%kind))
            call item%setFlags(read_flags)
            if(editor%assembly_id>0)call set_status_hint(item,integer_text(editor%shared%copies(row))//' linked control entries')
            call item%setForeground(QColor(167,192,205))
            call item%setTextAlignment(Qt_AlignCenter)
            call editor%table%setItem(row-1,0,item)
            do column=1,weight_column
                text=display_value(component(editor,row,column))
                item=QTableWidgetItem(text)
                call item%setFlags(edit_flags)
                call item%setTextAlignment(ior(Qt_AlignRight,Qt_AlignVCenter))
                call set_status_hint(item,text)
                call editor%table%setItem(row-1,column,item)
            end do
        end do
        call editor%table%resizeRowsToContents()
    end subroutine add_controls

    function index_text(row,counts,kind) result(text)
        integer, intent(in) :: row,counts(3),kind
        character(len=:), allocatable :: text
        integer :: index_u,index_v,index_w
        index_u=mod(row-1,counts(1))+1
        index_v=mod((row-1)/counts(1),counts(2))+1
        index_w=(row-1)/(counts(1)*counts(2))+1
        text=integer_text(index_u)
        if (kind>1) text=text//', '//integer_text(index_v)
        if (kind>2) text=text//', '//integer_text(index_w)
    end function index_text

    function integer_text(value) result(text)
        integer, intent(in) :: value
        character(len=:), allocatable :: text
        character(len=20) :: buffer
        write(buffer,'(i0)') value
        text=trim(buffer)
    end function integer_text

    real(c_double) function component(editor,row,column) result(value)
        type(control_editor), intent(in) :: editor
        integer, intent(in) :: row,column
        if (column==weight_column) then
            value=editor%original%weights(row)
        else
            value=display_length(editor%original%points(column,row),editor%model_unit,editor%display_unit)
        end if
    end function component

    function display_value(value) result(text)
        real(c_double), intent(in) :: value
        character(len=:), allocatable :: text
        character(len=48) :: buffer
        real(c_double) :: roundtrip
        integer :: status,last,exponent

        if (value==0._c_double) then
            text='0'
            return
        end if
        if (abs(value)>=1.e98_c_double.or.abs(value)<1.e-96_c_double) then
            write(buffer,'(es30.16e3)') value
        else
            write(buffer,'(g0.10)') value
            read(buffer,*,iostat=status,decimal='point') roundtrip
            if (status/=0) roundtrip=0._c_double
            if (roundtrip/=value) write(buffer,'(g0.17)') value
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
    end function display_value

    logical function decimal_value(editor,input,value) result(valid)
        type(control_editor), intent(in) :: editor
        character(len=*), intent(in) :: input
        real(c_double), intent(out) :: value
        character(len=:), allocatable :: text
        integer :: comma,status

        value=0._c_double
        valid=.false.
        text=trim(adjustl(input))
        if (len(text)==0.or.len(text)>128) return
        if (editor%decimal%validate(text)/=QValidator_Acceptable) return
        comma=index(text,',')
        if (comma>0) text(comma:comma)='.'
        read(text,*,iostat=status,decimal='point') value
        if (status/=0) return
        valid=ieee_is_finite(value)
    end function decimal_value

    logical function read_control_input(editor,candidate) result(valid)
        type(control_editor), intent(inout) :: editor
        type(single_patch_controls), intent(out) :: candidate
        type(QTableWidgetItem) :: item
        character(len=:), allocatable :: text
        integer :: row,column
        real(c_double) :: value

        valid=.false.
        candidate=editor%original
        do row=1,size(candidate%weights)
            do column=1,weight_column
                item=editor%table%item(row-1,column)
                text=item%text()
                if (text==display_value(component(editor,row,column))) cycle
                if (.not.decimal_value(editor,text,value)) then
                    call mark_invalid(editor,row-1,column,'Enter a finite number.')
                    return
                end if
                if (column==weight_column) then
                    if (value<=0._c_double) then
                        call mark_invalid(editor,row-1,column,'Enter a weight greater than zero.')
                        return
                    end if
                    candidate%weights(row)=value
                else
                    if(value==component(editor,row,column))cycle
                    value=model_length(value,editor%model_unit,editor%display_unit)
                    if(.not.ieee_is_finite(value))then
                        call mark_invalid(editor,row-1,column,'This coordinate exceeds the model range.')
                        return
                    end if
                    candidate%points(column,row)=value
                end if
            end do
        end do
        valid=.true.
    end function read_control_input

    subroutine invalidate_control_preview(editor)
        type(control_editor), intent(inout) :: editor

        call editor%apply%setEnabled(.false.)
        call editor%preview%setDefault(.true.)
        call editor%preview_table%setRowCount(0)
    end subroutine invalidate_control_preview

    subroutine preview_controls(data)
        class(*), intent(inout), optional :: data
        type(single_patch_controls) :: candidate
        type(multipatch_controls_t) :: requested
        type(multipatch_control_preview_t) :: preview
        character(len=:), allocatable :: error

        select type(editor=>data)
        type is(control_editor)
            if(.not.editor%preview%isEnabled())return
            call editor%preview%setFocus()
            call invalidate_control_preview(editor)
            call clear_error(editor)
            if(.not.read_control_input(editor,candidate))return
            if(all(candidate%points==editor%original%points).and.all(candidate%weights==editor%original%weights))then
                call publish_window_status(editor%dialog,'Edit a coordinate or weight to preview an adjustment.')
                return
            end if
            requested=editor%shared
            requested%points=candidate%points
            requested%weights=candidate%weights
            call preview_document_multipatch_controls(editor%document,editor%assembly_id,requested,preview,error)
            if(len(error)>0)then
                call show_error(editor,error)
                return
            end if
            call show_control_preview(editor,preview)
            editor%shared%points=candidate%points
            editor%shared%weights=candidate%weights
            call editor%apply%setEnabled(.true.)
            call editor%apply%setDefault(.true.)
        end select
    end subroutine preview_controls

    subroutine show_control_preview(editor,preview)
        type(control_editor), intent(inout) :: editor
        type(multipatch_control_preview_t), intent(in) :: preview
        character(len=:), allocatable :: summary,field
        real(c_double) :: before,after
        integer :: i,id,column,row

        if(size(preview%ids)==0)then
            call publish_window_status(editor%dialog,'No dependent controls need to adjust.')
            return
        end if
        summary=integer_text(size(preview%ids))//' dependent control points will adjust.'
        if(size(preview%ids)==1)summary='1 dependent control point will adjust.'
        call publish_window_status(editor%dialog,summary)
        call editor%preview_table%setRowCount(count(preview%reweighted)+ &
            count(preview%points/=editor%original%points(:,preview%ids)))
        row=0
        do i=1,size(preview%ids)
            id=preview%ids(i)
            do column=1,3
                if(preview%points(column,i)==editor%original%points(column,id))cycle
                field=length_caption(component_names(column),editor%display_unit)
                before=component(editor,id,column)
                after=display_length(preview%points(column,i),editor%model_unit,editor%display_unit)
                call set_preview_change(editor%preview_table,row,id,field,before,after)
                row=row+1
            end do
            if(.not.preview%reweighted(i))cycle
            call set_preview_change(editor%preview_table,row,id,'Weight',editor%original%weights(id),preview%weights(i))
            row=row+1
        end do
    end subroutine show_control_preview

    subroutine set_preview_change(table,row,id,field,before,after)
        type(QTableWidget), intent(inout) :: table
        integer, intent(in) :: row,id
        character(len=*), intent(in) :: field
        real(c_double), intent(in) :: before,after
        type(QTableWidgetItem) :: item
        integer :: column
        character(len=:), allocatable :: text

        call table%setItem(row,0,QTableWidgetItem(integer_text(id)))
        call table%setItem(row,1,QTableWidgetItem(field))
        do column=2,3
            text=display_value(merge(before,after,column==2))
            item=QTableWidgetItem(text)
            call item%setTextAlignment(ior(Qt_AlignRight,Qt_AlignVCenter))
            call set_status_hint(item,text)
            call table%setItem(row,column,item)
        end do
    end subroutine set_preview_change

    subroutine apply_controls(data)
        class(*), intent(inout), optional :: data
        type(single_patch_controls) :: candidate
        character(len=:), allocatable :: error

        select type(editor=>data)
        type is(control_editor)
            if (.not.editor%apply%isEnabled()) return
            call editor%apply%setFocus()
            if (.not.editor%apply%isEnabled()) return
            call clear_error(editor)
            if(editor%assembly_id==0)then
                if(active_patch_id(editor%kind,catalog=editor%catalog)/=editor%variant.or. &
                   editor%catalog%patch_revision(editor%kind,editor%variant)/=editor%revision)then
                    call show_error(editor,'This patch changed while the editor was open. Cancel and reopen Controls to edit it.')
                    return
                end if
            end if
            if(.not.read_control_input(editor,candidate))return
            if(editor%assembly_id>0)then
                if(editor%shared%smooth_joins)then
                    if(any(editor%shared%points/=candidate%points).or.any(editor%shared%weights/=candidate%weights))then
                        call invalidate_control_preview(editor)
                        return
                    end if
                end if
                editor%shared%points=candidate%points
                editor%shared%weights=candidate%weights
                call commit_document_multipatch_controls(editor%document,editor%assembly_id,editor%shared,editor%changed,error)
            else
                call commit_single_patch_controls(editor%kind,candidate,editor%changed,error, &
                    patch_id=editor%variant,catalog=editor%catalog)
            end if
            if (allocated(error)) then
                if (len(error)>0) then
                    if(editor%shared%smooth_joins)call invalidate_control_preview(editor)
                    call show_error(editor,error)
                    return
                end if
            end if
            call editor%dialog%accept()
        end select
    end subroutine apply_controls

    subroutine cancel_controls(data)
        class(*), intent(inout), optional :: data
        select type(editor=>data)
        type is(control_editor)
            call editor%dialog%reject()
        end select
    end subroutine cancel_controls

    subroutine show_error(editor,message)
        type(control_editor), intent(inout) :: editor
        character(len=*), intent(in) :: message
        call publish_window_status(editor%dialog,message,warning=.true.,blocking=.true.)
    end subroutine show_error

    subroutine mark_invalid(editor,row,column,message)
        type(control_editor), intent(inout) :: editor
        integer, intent(in) :: row,column
        character(len=*), intent(in) :: message
        type(QTableWidgetItem) :: item
        character(len=:), allocatable :: point

        point=index_text(row+1,editor%original%counts,editor%kind)
        call show_error(editor,'Point ('//point//') · '//trim(component_names(column))//' — '//message)
        editor%invalid_row=row
        editor%invalid_column=column
        editor%updating_feedback=.true.
        item=editor%table%item(row,column)
        call item%setBackground(QColor(108,54,60))
        call set_status_hint(item,message)
        editor%updating_feedback=.false.
        call editor%table%setCurrentCell(row,column)
        call editor%table%setFocus()
    end subroutine mark_invalid

    subroutine clear_error(editor)
        type(control_editor), intent(inout) :: editor
        type(QTableWidgetItem) :: item
        integer :: row
        row=editor%invalid_row
        if (row>=0) then
            editor%updating_feedback=.true.
            item=editor%table%item(row,editor%invalid_column)
            if (mod(row,2)==0) then
                call item%setBackground(QColor(22,33,42))
            else
                call item%setBackground(QColor(27,41,51))
            end if
            call set_status_hint(item,item%text())
            editor%updating_feedback=.false.
        end if
        editor%invalid_row=-1
        editor%invalid_column=-1
        call publish_window_status(editor%dialog,'',blocking=.true.)
    end subroutine clear_error

    recursive subroutine control_cell_changed(row,column,data)
        integer, intent(in) :: row,column
        class(*), intent(inout), optional :: data
        type(QTableWidgetItem) :: item
        real(c_double) :: value
        select type(editor=>data)
        type is(control_editor)
            if (editor%updating_feedback.or.column<1.or.column>weight_column) return
            item=editor%table%item(row,column)
            if (item%isNull()) return
            if(editor%shared%smooth_joins.and.editor%apply%isEnabled())call invalidate_control_preview(editor)
            if (row==editor%invalid_row.and.column==editor%invalid_column) then
                if (.not.decimal_value(editor,item%text(),value)) return
                if (column==weight_column.and.value<=0._c_double) return
                call clear_error(editor)
            else
                editor%updating_feedback=.true.
                call set_status_hint(item,item%text())
                editor%updating_feedback=.false.
                if (editor%invalid_row<0) call clear_error(editor)
            end if
        end select
    end subroutine control_cell_changed

end module studio_control_editor
