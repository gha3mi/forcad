module studio_assembly_refinement_dialog
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: iso_fortran_env, only: int64
    use qt, only: QAbstractItemView_NoEditTriggers,QButtonGroup,QComboBox,QDialog,QDialog_Accepted,QDialogButtonBox, &
        QDialogButtonBox_AcceptRole,QDialogButtonBox_ActionRole,QDialogButtonBox_Cancel,QGridLayout,QHeaderView, &
        QHeaderView_Stretch,QHeaderView_ResizeToContents,QLabel,QLineEdit,QPushButton,QRegularExpressionValidator, &
        QSizePolicy_Ignored,QSizePolicy_Fixed,QSpinBox,QTableWidget,QTableWidgetItem,Qt_PlainText,QVBoxLayout,QWidget
    use studio_state, only: studio_t,scene_patch_drawn
    use studio_presentation_text, only: integer_text
    use studio_edit_history, only: document_undo_stamp
    use studio_multipatch_collection, only: count_multipatches
    use studio_multipatch_refinement, only: commit_multipatch_refinement,multipatch_refinement_bounds, &
        multipatch_refinement_plan_t,plan_multipatch_refinement
    use studio_parameter_editor, only: number_text,parameter_number_pattern,parse_parameter_number,parameter_direction_buttons
    use studio_single_patch, only: single_patch_refinement_limits,single_patch_name
    use studio_status_hints, only: set_status_hint, publish_window_status
    use studio_theme, only: studio_style_primary
    implicit none(type,external)
    private
    public :: edit_assembly_refinement

    character(*),parameter :: directions='UVW'

    type :: refinement_editor
        type(studio_t),pointer :: workspace=>null()
        type(QDialog) :: dialog
        type(QComboBox) :: member
        type(QButtonGroup) :: direction
        type(QLineEdit) :: knot
        type(QSpinBox) :: steps
        type(QLabel) :: bounds
        type(QPushButton) :: apply,preview
        type(QTableWidget) :: table
        type(QRegularExpressionValidator) :: decimal
        integer,allocatable :: members(:)
        integer :: assembly,kind,last_member=0,last_axis=0
        integer(int64) :: revision,stamp
        real(c_double) :: interval(2)
        logical :: inserting,changed=.false.
    end type refinement_editor

contains

    logical function edit_assembly_refinement(workspace,assembly,inserting) result(changed)
        type(studio_t),intent(inout),target :: workspace
        integer,intent(in) :: assembly
        logical,intent(in) :: inserting
        type(refinement_editor),target :: editor
        type(QVBoxLayout) :: layout
        type(QGridLayout) :: form
        type(QHeaderView) :: header
        type(QWidget) :: direction_choices
        type(QLabel) :: label
        type(QDialogButtonBox) :: buttons
        character(:),allocatable :: title,action
        integer :: i,row

        changed=.false.
        if(assembly<1.or.assembly>count_multipatches(workspace%document%multipatches))return
        associate(group=>workspace%document%multipatches%groups(assembly))
            if(.not.allocated(group%native).or..not.allocated(group%definition%members))return
            if(size(group%definition%members)==0)return
            editor%members=group%definition%members
            editor%kind=group%definition%kind
        end associate
        editor%workspace=>workspace
        editor%assembly=assembly
        editor%inserting=inserting
        editor%revision=workspace%document%document_revision
        editor%stamp=document_undo_stamp(workspace%document)
        action='Elevate'
        title='Elevate multipatch degree'
        if(inserting)then
            action='Insert'
            title='Insert multipatch knot'
        end if

        editor%dialog=QDialog(workspace%window)
        call editor%dialog%resize(570,450)
        call editor%dialog%setObjectName('studioAssemblyRefinement')
        call editor%dialog%setWindowTitle(title)
        call editor%dialog%setMinimumWidth(480)
        layout=QVBoxLayout(editor%dialog)
        call layout%setContentsMargins(18,16,18,16)
        call layout%setSpacing(10)

        label=QLabel(title,editor%dialog)
        call label%setObjectName('studioDialogTitle')
        call layout%addWidget(label)
        label=QLabel(workspace%document%multipatches%groups(assembly)%definition%name,editor%dialog)
        call label%setTextFormat(Qt_PlainText)
        call label%setProperty('dialogRole','detail')
        call label%setWordWrap(.true.)
        call layout%addWidget(label)

        form=QGridLayout()
        call form%setSpacing(8)
        call form%setColumnStretch(1,1)
        label=QLabel('Starting patch',editor%dialog)
        call form%addWidget(label,0,0)
        editor%member=QComboBox(editor%dialog)
        call editor%member%setObjectName('studioAssemblyRefineMember')
        call editor%member%setAccessibleName('Starting patch')
        call editor%member%setSizePolicy(QSizePolicy_Ignored,QSizePolicy_Fixed)
        call editor%member%setMinimumWidth(220)
        do i=1,size(editor%members)
            call editor%member%addItem(single_patch_name(editor%kind,editor%members(i),catalog=workspace%document%patches))
        end do
        call form%addWidget(editor%member,0,1)

        label=QLabel('Direction',editor%dialog)
        call form%addWidget(label,1,0)
        direction_choices=parameter_direction_buttons(editor%dialog,editor%kind,1,'studioAssemblyRefineAxis', &
            editor%direction,refresh_direction,editor)
        call form%addWidget(direction_choices,1,1)

        row=2
        if(inserting)then
            label=QLabel('Knot value',editor%dialog)
            call form%addWidget(label,row,0)
            editor%knot=QLineEdit(parent=editor%dialog)
            call editor%knot%setObjectName('studioAssemblyRefineKnot')
            call editor%knot%setAccessibleName('Knot value')
            editor%decimal=QRegularExpressionValidator(parameter_number_pattern,editor%dialog)
            call editor%knot%setValidator(editor%decimal)
            call form%addWidget(editor%knot,row,1)
            row=row+1
        end if

        label=QLabel(merge('Insertions','Increment ',inserting),editor%dialog)
        call form%addWidget(label,row,0)
        editor%steps=QSpinBox(editor%dialog)
        call editor%steps%setObjectName('studioAssemblyRefineSteps')
        call editor%steps%setAccessibleName(trim(merge('Knot insertions ','Degree increment',inserting)))
        call editor%steps%setRange(1,huge(1))
        call form%addWidget(editor%steps,row,1)
        call layout%addLayout(form)

        editor%bounds=QLabel('',editor%dialog)
        call editor%bounds%setObjectName('studioAssemblyRefineRange')
        call editor%bounds%setProperty('dialogRole','detail')
        call editor%bounds%setWordWrap(.true.)
        if(inserting)call layout%addWidget(editor%bounds)

        editor%table=QTableWidget(0,merge(4,3,inserting),editor%dialog)
        call editor%table%setObjectName('studioAssemblyRefinePreview')
        call editor%table%setAccessibleName('Affected patches and parameter directions')
        call editor%table%setMinimumHeight(134)
        call editor%table%setVerticalHeaderVisible(.false.)
        call editor%table%setAlternatingRowColors(.true.)
        call editor%table%setWordWrap(.false.)
        call editor%table%setEditTriggers(QAbstractItemView_NoEditTriggers)
        if(inserting)then
            call editor%table%setHorizontalHeaderLabels([character(12)::'Patch','Axis','Knot','Insertions'])
        else
            call editor%table%setHorizontalHeaderLabels([character(12)::'Patch','Axis','Increment'])
        end if
        header=editor%table%horizontalHeader()
        call header%setSectionResizeMode(0,QHeaderView_Stretch)
        do i=1,editor%table%columnCount()-1
            call header%setSectionResizeMode(i,QHeaderView_ResizeToContents)
        end do
        call layout%addWidget(editor%table,1)


        buttons=QDialogButtonBox(QDialogButtonBox_Cancel,editor%dialog)
        editor%preview=buttons%addButton('Preview',QDialogButtonBox_ActionRole)
        call editor%preview%setObjectName('studioAssemblyRefineCheck')
        call editor%preview%setMinimumWidth(80)
        call editor%preview%setDefault(.true.)
        call set_status_hint(editor%preview,'Check the joins and highlight every visible patch the edit will affect.')
        call editor%preview%clicked%connect(preview_refinement,editor)
        editor%apply=buttons%addButton(action,QDialogButtonBox_AcceptRole)
        call editor%apply%setObjectName('studioAssemblyRefineApply')
        call editor%apply%setMinimumWidth(80)
        call studio_style_primary(editor%apply)
        call set_status_hint(editor%apply,'Refine all affected members in one Undo step.')
        call buttons%accepted%connect(apply_refinement,editor)
        call buttons%rejected%connect(cancel_refinement,editor)
        call layout%addWidget(buttons)
        call editor%member%currentIndexChanged%connect(refresh_direction,editor)
        if(inserting)call editor%knot%textChanged%connect(knot_changed,editor)
        call editor%steps%valueChanged%connect(steps_changed,editor)
        call refresh_direction(0,editor)
        if(editor%dialog%exec()==QDialog_Accepted)changed=editor%changed
        call clear_viewport_preview(editor)
        call editor%dialog%delete()
    end function edit_assembly_refinement

    subroutine refresh_direction(index,data)
        integer,intent(in) :: index
        class(*),intent(inout),optional :: data
        integer :: member,axis,limit
        logical :: available

        if(index<0)return
        select type(editor=>data)
        type is(refinement_editor)
            member=editor%member%currentIndex()+1
            axis=editor%direction%checkedId()+1
            if(member<1.or.member>size(editor%members).or.axis<1)return
            if(member==editor%last_member.and.axis==editor%last_axis)return
            editor%last_member=member
            editor%last_axis=axis
            member=editor%members(member)
            call invalidate_preview(editor)
            available=multipatch_refinement_bounds(editor%workspace%document,editor%assembly,member,axis,editor%interval)
            if(editor%inserting.and.available)then
                call editor%bounds%setText('Interior range · '//number_text(editor%interval(1))//' to '//number_text(editor%interval(2)))
                call editor%knot%setText(number_text(0.5_c_double*editor%interval(1)+0.5_c_double*editor%interval(2)))
            else if(available)then
                call single_patch_refinement_limits(editor%kind,axis,max_steps=limit, &
                    catalog=editor%workspace%document%patches,patch_id=member)
                call editor%steps%setRange(1,max(1,limit))
                available=limit>0
            end if
            call editor%preview%setEnabled(available)
            if(.not.available)call show_error(editor,'Refinement is unavailable in this direction.')
        end select
    end subroutine refresh_direction

    subroutine knot_changed(text,data)
        character(*),intent(in) :: text
        class(*),intent(inout),optional :: data

        select type(editor=>data)
        type is(refinement_editor)
            call invalidate_preview(editor)
            call editor%preview%setEnabled(len_trim(text)>0.and.editor%interval(1)<editor%interval(2))
        end select
    end subroutine knot_changed

    subroutine steps_changed(count,data)
        integer,intent(in) :: count
        class(*),intent(inout),optional :: data

        if(count<1)return
        select type(editor=>data)
        type is(refinement_editor)
            call invalidate_preview(editor)
        end select
    end subroutine steps_changed

    subroutine invalidate_preview(editor)
        type(refinement_editor),intent(inout) :: editor

        call editor%apply%setEnabled(.false.)
        call editor%preview%setDefault(.true.)
        call editor%table%setRowCount(0)
        call publish_window_status(editor%dialog,'',blocking=.true.)
        call clear_viewport_preview(editor)
    end subroutine invalidate_preview

    subroutine clear_viewport_preview(editor)
        type(refinement_editor),intent(in) :: editor

        associate(s=>editor%workspace)
            if(.not.allocated(s%refinement_preview_patches))return
            deallocate(s%refinement_preview_patches)
            s%refinement_preview_kind=0
            call s%surface%update()
        end associate
    end subroutine clear_viewport_preview

    logical function current_editor(editor) result(current)
        type(refinement_editor),intent(in) :: editor

        associate(document=>editor%workspace%document)
            current=document%document_revision==editor%revision.and.document_undo_stamp(document)==editor%stamp
        end associate
    end function current_editor

    subroutine preview_refinement(data)
        class(*),intent(inout),optional :: data
        type(multipatch_refinement_plan_t) :: plan
        real(c_double) :: value
        character(:),allocatable :: error
        integer :: member,axis,steps

        select type(editor=>data)
        type is(refinement_editor)
            if(.not.editor%preview%isEnabled())return
            call invalidate_preview(editor)
            if(.not.current_editor(editor))then
                call show_error(editor,'The project changed. Reopen refinement to use the current interfaces.')
                return
            end if
            member=editor%members(editor%member%currentIndex()+1)
            axis=editor%direction%checkedId()+1
            steps=editor%steps%value()
            associate(document=>editor%workspace%document)
                if(editor%inserting)then
                    if(.not.parse_parameter_number(editor%decimal,editor%knot%text(),value))then
                        call show_error(editor,'Enter a finite knot value inside the parameter range.')
                        return
                    end if
                    call plan_multipatch_refinement(document,editor%assembly,member,axis,steps,plan,error,knot=value)
                else
                    call plan_multipatch_refinement(document,editor%assembly,member,axis,steps,plan,error)
                end if
            end associate
            if(len(error)>0)then
                call show_error(editor,error)
                return
            end if
            call show_preview(editor,plan)
            call editor%apply%setEnabled(.true.)
            call editor%apply%setDefault(.true.)
        end select
    end subroutine preview_refinement

    subroutine show_preview(editor,plan)
        type(refinement_editor),intent(inout) :: editor
        type(multipatch_refinement_plan_t),intent(in) :: plan
        type(QTableWidgetItem) :: item
        integer :: i,id,axis,visible,total
        character(:),allocatable :: name,text

        call editor%table%setRowCount(size(plan%rows))
        visible=0
        associate(s=>editor%workspace)
            s%refinement_preview_kind=plan%kind
            allocate(s%refinement_preview_patches(maxval(plan%rows%member)),source=.false.)
            do i=1,size(plan%rows)
                id=plan%rows(i)%member
                axis=plan%rows(i)%axis
                if(.not.s%refinement_preview_patches(id))then
                    s%refinement_preview_patches(id)=.true.
                    if(scene_patch_drawn(s,plan%kind,id))visible=visible+1
                end if
                name=single_patch_name(plan%kind,id,catalog=s%document%patches)
                item=QTableWidgetItem(name)
                call set_status_hint(item,name//' · Patch '//integer_text(id))
                call editor%table%setItem(i-1,0,item)
                call editor%table%setItem(i-1,1,QTableWidgetItem(directions(axis:axis)))
                if(plan%inserting)then
                    text=number_text(plan%rows(i)%knot)
                    item=QTableWidgetItem(text)
                    call set_status_hint(item,text)
                    call editor%table%setItem(i-1,2,item)
                    call editor%table%setItem(i-1,3,QTableWidgetItem(integer_text(plan%steps)))
                else
                    call editor%table%setItem(i-1,2,QTableWidgetItem('+'//integer_text(plan%steps)))
                end if
            end do
            total=count(s%refinement_preview_patches)
            call s%surface%update()
        end associate
        text=integer_text(total)//' affected patches · '//integer_text(size(plan%rows))//' directions'
        if(visible>0)text=text//' · Visible patches outlined in amber'
        if(visible<total)text=text//' · '//integer_text(total-visible)//' hidden'
        call publish_window_status(editor%dialog,text)
    end subroutine show_preview

    subroutine apply_refinement(data)
        class(*),intent(inout),optional :: data
        real(c_double) :: value
        character(:),allocatable :: error
        integer :: member,axis,steps

        select type(editor=>data)
        type is(refinement_editor)
            if(.not.editor%apply%isEnabled())return
            if(.not.current_editor(editor))then
                call invalidate_preview(editor)
                call show_error(editor,'The project changed. Reopen refinement to use the current interfaces.')
                return
            end if
            member=editor%members(editor%member%currentIndex()+1)
            axis=editor%direction%checkedId()+1
            steps=editor%steps%value()
            if(editor%inserting)then
                if(.not.parse_parameter_number(editor%decimal,editor%knot%text(),value))then
                    call show_error(editor,'Enter a finite knot value inside the parameter range.')
                    call editor%knot%setFocus()
                    return
                end if
                call commit_multipatch_refinement(editor%workspace%document,editor%assembly,member,axis,steps,error,editor%changed,knot=value)
            else
                call commit_multipatch_refinement(editor%workspace%document,editor%assembly,member,axis,steps,error,editor%changed)
            end if
            if(len(error)>0)then
                call show_error(editor,error)
                return
            end if
            call editor%dialog%accept()
        end select
    end subroutine apply_refinement

    subroutine cancel_refinement(data)
        class(*),intent(inout),optional :: data

        select type(editor=>data)
        type is(refinement_editor)
            call editor%dialog%reject()
        end select
    end subroutine cancel_refinement

    subroutine show_error(editor,text)
        type(refinement_editor),intent(in) :: editor
        character(*),intent(in) :: text

        call publish_window_status(editor%dialog,text,warning=.true.,blocking=.true.)
    end subroutine show_error
end module studio_assembly_refinement_dialog
