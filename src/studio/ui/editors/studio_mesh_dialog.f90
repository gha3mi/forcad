module studio_mesh_dialog
    use studio_status_hints, only: set_status_hint, publish_window_status
    use, intrinsic :: iso_c_binding, only: c_double
    use qt, only: QHeaderView, QHeaderView_Stretch, QHeaderView_ResizeToContents, QHBoxLayout, QLineEdit, QRegularExpressionValidator, QTableWidget, QTableWidgetItem, &
        QAbstractButton, QAbstractItemView_NoEditTriggers, QButtonGroup, QCheckBox, QDialog, QDialog_Accepted, QDialogButtonBox, QDialogButtonBox_AcceptRole, &
        QDialogButtonBox_Cancel, QGridLayout, QLabel, QPushButton, QScrollArea, QSize, QSizePolicy_Expanding, &
        QSizePolicy_Ignored, QSpinBox, QStackedWidget, Qt_AlignCenter, &
        Qt_PlainText, QVBoxLayout, QWidget
    use studio_presentation_text, only: integer_text
    use studio_parameter_editor, only: number_text, parse_parameter_number, parameter_number_pattern
    use studio_units, only: model_length, display_length, length_caption
    use studio_document, only: studio_document_t
    use studio_mesh_operations, only: mesh_source_t, estimate_document_mesh_samples
    use studio_theme, only: studio_style_primary
    implicit none(type, external)
    private

    public :: edit_mesh_patch_counts
    character(*), parameter :: directions='UVW'

    type :: mesh_size_editor
        type(QDialog), pointer :: dialog=>null()
        type(studio_document_t), pointer :: document=>null()
        type(mesh_source_t) :: source
        type(QStackedWidget) :: pages
        type(QCheckBox) :: synchronize
        type(QLineEdit) :: target
        type(QRegularExpressionValidator) :: validator
        type(QTableWidget) :: table
        type(QPushButton) :: apply
        integer :: target_page=0,display_unit=0,method=-1
        logical :: generate_on_accept=.false.
        real(c_double) :: estimated_size=0._c_double
        integer, allocatable :: counts(:,:)
        character(:), allocatable :: labels(:)
    end type mesh_size_editor

contains

    logical function edit_mesh_patch_counts(parent,title,labels,cells,shared_limit,patch_limit,allow_per_patch, &
        per_patch,can_join,conforming,document,source,target_size,synchronize,display_unit,start_with_target,operation) result(accepted)
        class(QWidget), intent(in) :: parent
        character(*), intent(in) :: title,labels(:)
        integer, intent(inout) :: cells(:,:)
        integer, intent(in) :: shared_limit,patch_limit,display_unit
        logical, intent(in) :: allow_per_patch,can_join,start_with_target
        logical, intent(inout) :: per_patch,conforming,synchronize
        character(*), optional, intent(in) :: operation
        type(studio_document_t), intent(in), target :: document
        type(mesh_source_t), intent(in) :: source
        real(c_double), intent(inout) :: target_size
        type(QDialog), target :: dialog
        type(mesh_size_editor), target :: editor
        type(QVBoxLayout) :: layout
        type(QHBoxLayout) :: method_row
        type(QLabel) :: heading, context
        type(QButtonGroup) :: methods
        type(QPushButton) :: method
        type(QAbstractButton) :: selected_method
        type(QCheckBox) :: join_faces
        type(QDialogButtonBox) :: buttons
        type(QSpinBox), allocatable :: shared_inputs(:),patch_inputs(:,:)
        integer :: dimension,row,axis,index,initial_method
        character(:), allocatable :: submit_label,submit_hint,dialog_title,method_label,method_name,method_hint

        accepted=.false.
        dimension=size(cells,1)
        if(dimension<1.or.dimension>3.or.size(cells,2)/=size(labels).or.size(labels)==0)return
        editor%document=>document
        editor%source=source
        editor%labels=labels
        editor%display_unit=display_unit
        editor%target_page=merge(2,1,allow_per_patch)
        editor%generate_on_accept=present(operation)

        dialog_title='Mesh sizing'
        submit_label='Apply counts'
        submit_hint='Use these counts the next time you generate the mesh.'
        if(present(operation))then
            dialog_title=operation//' mesh'
            submit_label=operation
            submit_hint='Build the mesh with this sizing.'
            if(source%mesh_id>0)then
                dialog_title='Remesh'
                submit_hint='Replace the existing mesh with this sizing.'
            end if
        end if

        dialog=QDialog(parent)
        editor%dialog=>dialog
        call dialog%setObjectName('studioPatchMeshDialog')
        call dialog%setWindowTitle(dialog_title)
        call dialog%setMinimumWidth(merge(568,448,allow_per_patch))
        layout=QVBoxLayout(dialog)
        call layout%setContentsMargins(18,16,18,16)
        call layout%setSpacing(10)
        heading=QLabel(dialog_title,dialog)
        call heading%setObjectName('studioDialogTitle')
        call layout%addWidget(heading)
        context=QLabel(title,dialog)
        call context%setTextFormat(Qt_PlainText)
        call context%setWordWrap(.true.)
        call context%setProperty('dialogRole','detail')
        call layout%addWidget(context)

        initial_method=merge(1,0,allow_per_patch.and.per_patch)
        if(start_with_target)initial_method=editor%target_page
        method_row=QHBoxLayout()
        call method_row%setSpacing(4)
        methods=QButtonGroup(dialog)
        call methods%setExclusive(.true.)
        do index=0,editor%target_page
            if(index==0)then
                method_label=trim(merge('Shared counts','Cell counts  ',allow_per_patch))
                method_name='studioMeshSizingCounts'
                method_hint='Set the cell counts for all selected patches.'
            else if(index==editor%target_page)then
                method_label='Target size'
                method_name='studioMeshSizingTarget'
                method_hint='Estimate cell counts from a physical target size.'
            else
                method_label='Per patch'
                method_name='studioMeshSizingPerPatch'
                method_hint='Set cell counts separately for each patch.'
            end if
            method=QPushButton(method_label,dialog)
            call method%setObjectName(method_name)
            call method%setAccessibleName(method_label//' mesh sizing')
            call method%setProperty('autoDefault',.false.)
            call method%setCheckable(.true.)
            call method%setChecked(index==initial_method)
            call set_status_hint(method,method_hint)
            call methods%addButton(method,index)
            call method_row%addWidget(method)
        end do
        call method_row%addStretch(1)
        call layout%addLayout(method_row)

        editor%pages=QStackedWidget(dialog)

        allocate(shared_inputs(dimension))
        call add_shared_counts(editor%pages,cells(:,1),shared_limit,shared_inputs)
        if(allow_per_patch)then
            allocate(patch_inputs(dimension,size(labels)))
            call add_patch_counts(editor%pages,labels,cells,patch_limit,patch_inputs)
        end if
        call add_target_sizing(editor,dialog,target_size,synchronize,allow_per_patch.and.dimension>1)
        call layout%addWidget(editor%pages,1)

        if(can_join)then
            join_faces=QCheckBox('Conforming hex mesh',dialog)
            call join_faces%setObjectName('studioPatchMeshConforming')
            call join_faces%setChecked(conforming)
            call set_status_hint(join_faces,'Share nodes on joined faces. Cell counts must match along both mapped face directions.')
            call layout%addWidget(join_faces)
        end if

        buttons=QDialogButtonBox(QDialogButtonBox_Cancel,parent=dialog)
        editor%apply=buttons%addButton(submit_label,QDialogButtonBox_AcceptRole)
        call editor%apply%setObjectName('studioPatchMeshApply')
        call studio_style_primary(editor%apply)
        call editor%apply%setDefault(.true.)
        call set_status_hint(editor%apply,submit_hint)
        call buttons%accepted%connect(accept_mesh_dialog,editor)
        call buttons%rejected%connect(cancel_mesh_dialog,dialog)
        call methods%idClicked%connect(change_count_scope,editor)
        call change_count_scope(initial_method,editor)
        call layout%addWidget(buttons)
        if(start_with_target)then
            call editor%target%setFocus()
        else
            selected_method=methods%button(initial_method)
            call selected_method%setFocus()
        end if

        call dialog%adjustSize()
        accepted=dialog%exec()==QDialog_Accepted
        if(accepted)then
            index=editor%pages%currentIndex()
            if(index==editor%target_page)then
                accepted=allocated(editor%counts)
                if(accepted)then
                    cells=editor%counts
                    per_patch=allow_per_patch
                    target_size=editor%estimated_size
                    synchronize=editor%synchronize%isChecked()
                end if
            else
                target_size=0._c_double
                per_patch=allow_per_patch.and.index==1
                do row=1,size(labels)
                    do axis=1,dimension
                        if(per_patch)then
                            cells(axis,row)=patch_inputs(axis,row)%value()
                        else
                            cells(axis,row)=shared_inputs(axis)%value()
                        end if
                    end do
                end do
            end if
            conforming=.false.
            if(can_join)conforming=join_faces%isChecked()
        end if
        call dialog%delete()
    end function edit_mesh_patch_counts

    subroutine add_target_sizing(editor,dialog,target_size,synchronize,can_synchronize)
        type(mesh_size_editor), intent(inout), target :: editor
        type(QDialog), intent(in) :: dialog
        real(c_double), intent(in) :: target_size
        logical, intent(in) :: synchronize,can_synchronize
        type(QWidget) :: page
        type(QVBoxLayout) :: layout
        type(QHBoxLayout) :: field
        type(QLabel) :: label
        type(QPushButton) :: estimate
        type(QHeaderView) :: header
        character(8) :: headings(4)
        integer :: axis

        page=QWidget(dialog)
        layout=QVBoxLayout(page)
        call layout%setContentsMargins(0,0,0,0)
        call layout%setSpacing(7)
        field=QHBoxLayout()
        label=QLabel(length_caption('Target cell size',editor%display_unit),page)
        editor%target=QLineEdit(parent=page)
        call editor%target%setObjectName('studioMeshTargetSize')
        call editor%target%setAccessibleName(label%text())
        call editor%target%setPlaceholderText('Positive physical size')
        editor%validator=QRegularExpressionValidator(parameter_number_pattern,page)
        call editor%target%setValidator(editor%validator)
        if(target_size>0._c_double)then
            call editor%target%setText(number_text(display_length(target_size,editor%document%model_unit,editor%display_unit)))
        end if
        estimate=QPushButton('Preview counts',page)
        call estimate%setObjectName('studioMeshEstimateCounts')
        call set_status_hint(estimate,'Estimate counts from the patch control geometry.')
        call set_status_hint(editor%target,'Enter a positive physical size, then preview the estimated counts.')
        call field%addWidget(label)
        call field%addWidget(editor%target,1)
        call field%addWidget(estimate)
        call layout%addLayout(field)
        editor%synchronize=QCheckBox('Synchronize joined directions',page)
        call editor%synchronize%setObjectName('studioMeshTargetSynchronize')
        call editor%synchronize%setChecked(synchronize)
        call editor%synchronize%setVisible(can_synchronize)
        call set_status_hint(editor%synchronize,'Use matching counts along connected edges and faces, including swapped axes.')
        call layout%addWidget(editor%synchronize)

        headings(1)='Patch'
        do axis=1,editor%source%kind
            headings(axis+1)=directions(axis:axis)//' cells'
        end do
        editor%table=QTableWidget(0,editor%source%kind+1,page)
        call editor%table%setObjectName('studioMeshTargetPreview')
        call editor%table%setAccessibleName('Estimated cells per patch')
        call editor%table%setEditTriggers(QAbstractItemView_NoEditTriggers)
        call editor%table%setVerticalHeaderVisible(.false.)
        call editor%table%setAlternatingRowColors(.true.)
        call editor%table%setStyleSheet('QTableWidget#studioMeshTargetPreview { alternate-background-color:#202f3b; }')
        call editor%table%setMinimumHeight(126)
        call editor%table%setHorizontalHeaderLabels(headings(:editor%source%kind+1))
        header=editor%table%horizontalHeader()
        call header%setSectionResizeMode(0,QHeaderView_Stretch)
        do axis=1,editor%source%kind
            call header%setSectionResizeMode(axis,QHeaderView_ResizeToContents)
        end do
        call layout%addWidget(editor%table,1)
        call editor%pages%addWidget(page)

        call editor%target%textChanged%connect(target_text_changed,editor)
        call editor%synchronize%toggled%connect(target_sync_changed,editor)
        call estimate%clicked%connect(preview_target_counts,editor)
    end subroutine add_target_sizing

    subroutine preview_target_counts(data)
        class(*), optional, intent(inout) :: data
        real(c_double) :: value
        integer :: row,axis,total
        type(QTableWidgetItem) :: item
        character(:), allocatable :: error
        logical :: valid

        if(.not.present(data))return
        select type(editor=>data)
        type is(mesh_size_editor)
            call clear_target_preview(editor)
            valid=parse_parameter_number(editor%validator,editor%target%text(),value)
            if(.not.valid.or.value<=0._c_double)then
                call publish_window_status(editor%dialog,'Enter a positive target size.',warning=.true.,blocking=.true.)
                call editor%target%setFocus()
                return
            end if
            editor%estimated_size=model_length(value,editor%document%model_unit,editor%display_unit)
            call estimate_document_mesh_samples(editor%document,editor%source,editor%estimated_size, &
                editor%synchronize%isChecked(),editor%counts,error)
            if(len(error)>0)then
                call publish_window_status(editor%dialog,error,warning=.true.)
                return
            end if
            call editor%table%setRowCount(size(editor%counts,2))
            total=0
            do row=1,size(editor%counts,2)
                item=QTableWidgetItem(trim(editor%labels(row)))
                call set_status_hint(item,trim(editor%labels(row)))
                call editor%table%setItem(row-1,0,item)
                do axis=1,size(editor%counts,1)
                    item=QTableWidgetItem(integer_text(editor%counts(axis,row)))
                    call item%setTextAlignment(Qt_AlignCenter)
                    call editor%table%setItem(row-1,axis,item)
                end do
                total=total+product(editor%counts(:,row))
            end do
            call publish_window_status(editor%dialog,integer_text(total)//' estimated cells')
            call editor%apply%setEnabled(.true.)
        end select
    end subroutine preview_target_counts

    subroutine target_text_changed(text,data)
        character(*), intent(in) :: text
        class(*), optional, intent(inout) :: data

        if(.not.present(data))return
        select type(editor=>data)
        type is(mesh_size_editor)
            call clear_target_preview(editor)
        end select
    end subroutine target_text_changed

    subroutine target_sync_changed(checked,data)
        logical, intent(in) :: checked
        class(*), optional, intent(inout) :: data

        if(.not.present(data))return
        select type(editor=>data)
        type is(mesh_size_editor)
            call clear_target_preview(editor)
        end select
    end subroutine target_sync_changed

    subroutine clear_target_preview(editor)
        type(mesh_size_editor), intent(inout) :: editor

        if(allocated(editor%counts))deallocate(editor%counts)
        call editor%table%setRowCount(0)
        call publish_window_status(editor%dialog,'',blocking=.true.)
        call update_submit_state(editor)
    end subroutine clear_target_preview

    subroutine update_submit_state(editor)
        type(mesh_size_editor), intent(inout) :: editor
        logical :: available
        real(c_double) :: value

        available=editor%pages%currentIndex()/=editor%target_page.or.allocated(editor%counts)
        if(.not.available.and.editor%generate_on_accept)then
            available=parse_parameter_number(editor%validator,editor%target%text(),value)
            available=available.and.value>0._c_double
        end if
        call editor%apply%setEnabled(available)
    end subroutine update_submit_state

    subroutine add_shared_counts(pages,counts,cell_limit,inputs)
        type(QStackedWidget), intent(inout) :: pages
        integer, intent(in) :: counts(:),cell_limit
        type(QSpinBox), intent(out) :: inputs(:)
        type(QWidget) :: page
        type(QGridLayout) :: row
        type(QLabel) :: label
        integer :: axis

        page=QWidget()
        row=QGridLayout(page)
        call row%setContentsMargins(0,0,0,0)
        call row%setSpacing(7)
        do axis=1,size(counts)
            label=QLabel(directions(axis:axis)//' cells',page)
            call row%addWidget(label,0,axis-1)
            call mesh_count_input(inputs(axis),page,'studioMeshCells'//directions(axis:axis), &
                directions(axis:axis)//' cells',counts(axis),cell_limit)
            call row%addWidget(inputs(axis),1,axis-1)
        end do
        call row%setRowStretch(2,1)
        call pages%addWidget(page)
    end subroutine add_shared_counts

    subroutine add_patch_counts(pages,labels,cells,cell_limit,inputs)
        type(QStackedWidget), intent(inout) :: pages
        character(*), intent(in) :: labels(:)
        integer, intent(in) :: cells(:,:),cell_limit
        type(QSpinBox), intent(out) :: inputs(:,:)
        type(QScrollArea) :: scroll
        type(QWidget) :: table
        type(QGridLayout) :: rows
        type(QLabel) :: label
        integer :: row,axis

        scroll=QScrollArea()
        call scroll%setObjectName('studioPatchMeshCounts')
        call scroll%setAccessibleName('Cell counts per patch')
        call scroll%setWidgetResizable(.true.)
        table=QWidget()
        rows=QGridLayout(table)
        call rows%setContentsMargins(10,8,10,8)
        call rows%setSpacing(7)
        call rows%setColumnStretch(0,1)
        label=QLabel('Patch',table)
        call rows%addWidget(label,0,0)
        do axis=1,size(cells,1)
            label=QLabel(directions(axis:axis)//' cells',table)
            call label%setAlignment(Qt_AlignCenter)
            call rows%addWidget(label,0,axis)
        end do
        do row=1,size(labels)
            label=QLabel(trim(labels(row)),table)
            call label%setTextFormat(Qt_PlainText)
            call label%setWordWrap(.true.)
            call rows%addWidget(label,row,0)
            do axis=1,size(cells,1)
                call mesh_count_input(inputs(axis,row),table,'studioPatchMeshCells'//directions(axis:axis)//integer_text(row), &
                    trim(labels(row))//' · '//directions(axis:axis)//' cells',cells(axis,row),cell_limit)
                call rows%addWidget(inputs(axis,row),row,axis)
            end do
        end do
        call rows%setRowStretch(size(labels)+1,1)
        call scroll%setWidget(table)
        call pages%addWidget(scroll)
    end subroutine add_patch_counts

    subroutine mesh_count_input(input,parent,object_name,label,value,limit)
        type(QSpinBox), intent(out) :: input
        type(QWidget), intent(in) :: parent
        character(*), intent(in) :: object_name,label
        integer, intent(in) :: value,limit

        input=QSpinBox(parent)
        call input%setObjectName(object_name)
        call input%setAccessibleName(label)
        call input%setRange(1,max(1,limit))
        call input%setValue(value)
        call input%setMinimumWidth(92)
    end subroutine mesh_count_input

    subroutine change_count_scope(index,data)
        integer, intent(in) :: index
        class(*), optional, intent(inout) :: data
        type(QWidget) :: page
        type(QSize) :: preferred
        integer :: i

        if(.not.present(data))return
        select type(editor=>data)
        type is(mesh_size_editor)
            if(index==editor%method)return
            call publish_window_status(editor%dialog,'',blocking=.true.)
            editor%method=index
            do i=0,editor%pages%count()-1
                page=editor%pages%widget(i)
                call page%setSizePolicy(QSizePolicy_Expanding,merge(QSizePolicy_Expanding,QSizePolicy_Ignored,i==index))
            end do
            call editor%pages%setCurrentIndex(index)
            call update_submit_state(editor)
            if(editor%dialog%isVisible())then
                preferred=editor%dialog%sizeHint()
                call editor%dialog%resize(editor%dialog%width(),preferred%height)
            end if
        end select
    end subroutine change_count_scope

    subroutine accept_mesh_dialog(data)
        class(*), optional, intent(inout) :: data

        if(.not.present(data))return
        select type(editor=>data)
        type is(mesh_size_editor)
            if(editor%pages%currentIndex()==editor%target_page)then
                if(editor%generate_on_accept.and..not.allocated(editor%counts))call preview_target_counts(editor)
                if(.not.allocated(editor%counts))return
            end if
            call editor%dialog%accept()
        end select
    end subroutine accept_mesh_dialog

    subroutine cancel_mesh_dialog(data)
        class(*), optional, intent(inout) :: data

        if(.not.present(data))return
        select type(dialog=>data)
        type is(QDialog)
            call dialog%reject()
        end select
    end subroutine cancel_mesh_dialog

end module studio_mesh_dialog
