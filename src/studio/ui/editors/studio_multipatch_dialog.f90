! Edit assembly membership and native interface declarations in Qt.
module studio_multipatch_dialog
    use studio_status_hints, only: set_status_hint, publish_window_status
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use qt, only: QAbstractItemView_NoEditTriggers, QAbstractItemView_SelectRows, &
        QAbstractItemView_SingleSelection, QCheckBox, QColor, QComboBox, QDialog, QDialog_Accepted, &
        QDialogButtonBox, QDialogButtonBox_AcceptRole, QDialogButtonBox_Cancel, QFrame, QGridLayout, &
        QHeaderView, QHeaderView_ResizeToContents, QHeaderView_Stretch, &
        QHBoxLayout, QLabel, QLineEdit, QPushButton, QScrollArea, QSizePolicy_Fixed, QSizePolicy_Ignored, &
        QSpinBox, QStackedWidget, Qt_ItemIsEnabled, Qt_ItemIsSelectable, Qt_PlainText, Qt_ScrollBarAlwaysOff, &
        QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget
    use studio_icons, only: line_icon
    use studio_multipatch, only: multipatch_draft, multipatch_join, multipatch_result, inspect_multipatch
    use studio_single_patch, only: single_patch_count, single_patch_name, studio_patch_catalog
    use studio_theme, only: studio_divider, studio_secondary_text, studio_style_primary
    implicit none(type, external)
    private
    public :: edit_studio_multipatch

    character(7), parameter :: ranks(3)=[character(7)::'Curve','Surface','Volume']
    character(8), parameter :: interface_types(3)=[character(8)::'endpoints','edges','faces']
    character(*), parameter :: dialog_style=&
        'QLabel#assemblyHeading { font-size:15px; font-weight:650; color:#edf5fb; }'// &
        'QLabel#assemblyHint { color:'//studio_secondary_text//'; }'// &
        'QFrame#assemblyCard { background:#151f29; border:1px solid '//studio_divider//'; }'// &
        'QFrame#assemblyCard QScrollArea, QFrame#assemblyCard QWidget#assemblyMemberList { background:transparent; }'// &
        'QFrame#assemblyCard QTableView { background:#131d26; alternate-background-color:#182630; }'// &
        'QFrame#assemblyCard QHeaderView::section { padding:7px 8px; }'// &
        'QFrame#assemblyCard QTableView::item { padding-left:6px; }'// &
        'QFrame#assemblyCard QCheckBox { padding:3px 6px; }'// &
        'QFrame#assemblyCard QCheckBox:hover { background:#213746; }'

    type :: membership_page
        type(QCheckBox), allocatable :: checks(:)
        type(QLabel) :: empty_label
    end type membership_page

    type :: multipatch_editor
        type(studio_patch_catalog), pointer :: catalog=>null()
        type(QDialog) :: dialog
        type(QLineEdit) :: name, filter
        type(QComboBox) :: kind_box
        type(QStackedWidget) :: member_stack
        type(QTableWidget) :: table
        type(QLabel) :: count_label, interface_detail
        type(QPushButton) :: add_button, edit_button, remove_button, check_button, save_button
        type(membership_page) :: pages(3)
        type(multipatch_draft) :: drafts(3)
        type(multipatch_result) :: report
        character(96), allocatable :: reserved_names(:)
        integer :: kind=2, first_variant=1
        logical :: checked=.false.
    end type multipatch_editor

    type :: interface_editor
        type(studio_patch_catalog), pointer :: catalog=>null()
        type(QDialog) :: dialog
        type(QComboBox) :: patches(2), sides(2), continuity
        type(QSpinBox) :: order
        type(QLabel) :: order_label, jet_label
        type(QPushButton) :: preview
        type(QLineEdit) :: jet
        type(QCheckBox) :: reverse, swap, flip(2)
        type(multipatch_join) :: result
        integer, allocatable :: members(:)
        integer :: kind=2
    end type interface_editor

contains

    logical function edit_studio_multipatch(parent,initial_rank,initial_member,draft,reserved_names,first_variant,catalog,initial_members) &
        result(accepted)
        class(QWidget), intent(in) :: parent
        integer, intent(in) :: initial_rank,initial_member
        type(multipatch_draft), intent(inout) :: draft
        character(*), intent(in), optional :: reserved_names(:)
        integer,intent(in),optional::first_variant,initial_members(:)
        type(studio_patch_catalog), intent(in), target :: catalog
        type(multipatch_editor), target :: editor
        type(QVBoxLayout) :: layout, member_layout, interface_layout
        type(QHBoxLayout) :: row, details, footer
        type(QLabel) :: label
        type(QFrame) :: member_card, interface_card
        type(QHeaderView) :: header
        type(QDialogButtonBox) :: buttons
        logical :: editing
        integer :: rank
        character(:), allocatable :: title

        accepted=.false.
        editor%catalog=>catalog
        editing=draft%kind>=1.and.draft%kind<=3
        if(present(reserved_names))editor%reserved_names=reserved_names
        if(present(first_variant))editor%first_variant=first_variant
        editor%kind=max(1,min(3,initial_rank))
        if(editing)editor%kind=draft%kind
        do rank=1,3
            editor%drafts(rank)%kind=rank
            editor%drafts(rank)%name=trim(ranks(rank))//' assembly'
            allocate(editor%drafts(rank)%members(0), editor%drafts(rank)%joins(0))
        end do
        if(editing)then
            editor%drafts(editor%kind)=draft
            if(.not.allocated(editor%drafts(editor%kind)%name))editor%drafts(editor%kind)%name='Assembly'
            if(.not.allocated(editor%drafts(editor%kind)%members))allocate(editor%drafts(editor%kind)%members(0))
            if(.not.allocated(editor%drafts(editor%kind)%joins))allocate(editor%drafts(editor%kind)%joins(0))
        else if(present(initial_members))then
            editor%drafts(editor%kind)%members=pack(initial_members,initial_members>=editor%first_variant.and. &
                initial_members<=single_patch_count(editor%kind,catalog=editor%catalog))
        else if(initial_rank==editor%kind.and.initial_member>=editor%first_variant.and. &
                initial_member<=single_patch_count(editor%kind,catalog=editor%catalog))then
            editor%drafts(editor%kind)%members=[initial_member]
        end if

        title='New multipatch'
        if(editing)title='Edit multipatch'
        editor%dialog=QDialog(parent)
        call editor%dialog%setObjectName('studioMultipatchDialog')
        call editor%dialog%setWindowTitle(title)
        call editor%dialog%setMinimumSize(790,464)
        call editor%dialog%resize(852,520)
        call editor%dialog%setStyleSheet(dialog_style)
        layout=QVBoxLayout(editor%dialog)
        call layout%setContentsMargins(18,16,18,16)
        call layout%setSpacing(10)
        label=QLabel(title,editor%dialog)
        call label%setObjectName('studioDialogTitle')
        call layout%addWidget(label)
        call set_status_hint(editor%dialog,'Choose members, then connect their shared boundaries.')

        row=QHBoxLayout()
        call row%setSpacing(8)
        label=QLabel('Name',editor%dialog)
        call row%addWidget(label)
        editor%name=QLineEdit(editor%drafts(editor%kind)%name,editor%dialog)
        call name_control(editor%name,'studioMultipatchName','Multipatch name')
        call editor%name%setMaxLength(96)
        call editor%name%setPlaceholderText('Enter a name for the assembly')
        call row%addWidget(editor%name,1)
        label=QLabel('Patch type',editor%dialog)
        call row%addWidget(label)
        editor%kind_box=QComboBox(editor%dialog)
        call name_control(editor%kind_box,'studioMultipatchKind','Patch type')
        call editor%kind_box%setMinimumWidth(128)
        do rank=1,size(ranks)
            call editor%kind_box%addItem(trim(ranks(rank)))
        end do
        call editor%kind_box%setCurrentIndex(editor%kind-1)
        call editor%kind_box%setEnabled(.not.editing)
        call row%addWidget(editor%kind_box)
        call layout%addLayout(row)

        details=QHBoxLayout()
        call details%setSpacing(12)
        call build_card(editor%dialog,member_card,member_layout)
        call member_card%setFixedWidth(222)
        row=QHBoxLayout()
        label=heading('Members',member_card)
        call row%addWidget(label,1)
        editor%count_label=hint('',member_card)
        call editor%count_label%setObjectName('studioMultipatchMemberCount')
        call row%addWidget(editor%count_label)
        call member_layout%addLayout(row)
        editor%filter=QLineEdit(parent=member_card)
        call name_control(editor%filter,'studioMultipatchMemberFilter','Filter members by name')
        call editor%filter%setPlaceholderText('Filter by name')
        call editor%filter%setProperty('clearButtonEnabled',.true.)
        call member_layout%addWidget(editor%filter)
        editor%member_stack=QStackedWidget(member_card)
        call editor%member_stack%setObjectName('studioMultipatchMembers')
        call editor%member_stack%setAccessibleName('Patches in this assembly')
        do rank=1,3
            call build_member_page(editor,rank)
        end do
        call editor%member_stack%setCurrentIndex(editor%kind-1)
        call member_layout%addWidget(editor%member_stack,1)
        call details%addWidget(member_card)

        call build_card(editor%dialog,interface_card,interface_layout)
        row=QHBoxLayout()
        call row%setSpacing(7)
        label=heading('Interfaces',interface_card)
        call row%addWidget(label,1)
        editor%add_button=QPushButton('Add',interface_card)
        call editor%add_button%setIcon(line_icon(30))
        call name_control(editor%add_button,'studioMultipatchAddJoin','Add interface')
        call set_status_hint(editor%add_button,'Define a connection between two members')
        call row%addWidget(editor%add_button)
        editor%edit_button=QPushButton('Edit',interface_card)
        call editor%edit_button%setIcon(line_icon(88))
        call name_control(editor%edit_button,'studioMultipatchEditJoin','Edit selected interface')
        call row%addWidget(editor%edit_button)
        editor%remove_button=QPushButton('Remove',interface_card)
        call editor%remove_button%setIcon(line_icon(36))
        call name_control(editor%remove_button,'studioMultipatchRemoveJoin','Remove selected interface')
        call row%addWidget(editor%remove_button)
        call interface_layout%addLayout(row)
        editor%table=QTableWidget(0,4,interface_card)
        call editor%table%setObjectName('studioMultipatchInterfaces')
        call editor%table%setAccessibleName('Interfaces and their validation status')
        call editor%table%setHorizontalHeaderLabels([character(10)::'Boundary A','Boundary B','Continuity','Status'])
        header=editor%table%horizontalHeader()
        call header%setSectionResizeMode(0,QHeaderView_Stretch)
        call header%setSectionResizeMode(1,QHeaderView_Stretch)
        call header%setSectionResizeMode(2,QHeaderView_ResizeToContents)
        call header%setSectionResizeMode(3,QHeaderView_ResizeToContents)
        call editor%table%setVerticalHeaderVisible(.false.)
        call editor%table%setSelectionMode(QAbstractItemView_SingleSelection)
        call editor%table%setSelectionBehavior(QAbstractItemView_SelectRows)
        call editor%table%setEditTriggers(QAbstractItemView_NoEditTriggers)
        call editor%table%setAlternatingRowColors(.true.)
        call editor%table%setWordWrap(.false.)
        call editor%table%setShowGrid(.false.)
        call interface_layout%addWidget(editor%table,1)
        editor%interface_detail=hint('',interface_card)
        call editor%interface_detail%setObjectName('studioMultipatchInterfaceDetail')
        call editor%interface_detail%setMinimumHeight(40)
        call interface_layout%addWidget(editor%interface_detail)
        call details%addWidget(interface_card,1)
        call layout%addLayout(details,1)

        footer=QHBoxLayout()
        editor%check_button=QPushButton('Check interfaces',editor%dialog)
        call editor%check_button%setIcon(line_icon(89))
        call name_control(editor%check_button,'studioMultipatchCheck','Check assembly interfaces')
        call footer%addWidget(editor%check_button)
        call footer%addStretch(1)
        buttons=QDialogButtonBox(QDialogButtonBox_Cancel,editor%dialog)
        editor%save_button=buttons%addButton('Create multipatch',QDialogButtonBox_AcceptRole)
        call editor%save_button%setIcon(line_icon(87))
        if(editing)call editor%save_button%setText('Save multipatch')
        call name_control(editor%save_button,'studioMultipatchSave','Save multipatch')
        call studio_style_primary(editor%save_button)
        call editor%save_button%setMinimumWidth(142)
        call editor%save_button%setDefault(.true.)
        call footer%addWidget(buttons)
        call layout%addLayout(footer)

        call editor%name%textChanged%connect(name_changed,editor)
        call editor%filter%textChanged%connect(filter_members,editor)
        call editor%kind_box%currentIndexChanged%connect(change_kind,editor)
        call editor%table%cellDoubleClicked%connect(double_click_join,editor)
        call editor%table%currentCellChanged%connect(change_join_row,editor)
        call editor%add_button%clicked%connect(add_join,editor)
        call editor%edit_button%clicked%connect(edit_join,editor)
        call editor%remove_button%clicked%connect(remove_join,editor)
        call editor%check_button%clicked%connect(check_draft,editor)
        call buttons%accepted%connect(save_draft,editor)
        call buttons%rejected%connect(cancel_draft,editor)
        call refresh_editor(editor)
        call filter_members('',editor)
        if(editing)then
            call check_draft(editor)
        else
            call editor%name%setFocus()
        end if

        if(editor%dialog%exec()==QDialog_Accepted)then
            accepted=.true.
            draft=editor%drafts(editor%kind)
        end if
        call editor%dialog%delete()
    end function edit_studio_multipatch

    subroutine build_card(parent,card,layout)
        class(QWidget), intent(in) :: parent
        type(QFrame), intent(out) :: card
        type(QVBoxLayout), intent(out) :: layout
        card=QFrame(parent)
        call card%setObjectName('assemblyCard')
        layout=QVBoxLayout(card)
        call layout%setContentsMargins(12,11,12,10)
        call layout%setSpacing(9)
    end subroutine build_card

    function heading(text,parent) result(label)
        character(*), intent(in) :: text
        class(QWidget), intent(in) :: parent
        type(QLabel) :: label
        label=QLabel(text,parent)
        call label%setObjectName('assemblyHeading')
    end function heading

    function hint(text,parent) result(label)
        character(*), intent(in) :: text
        class(QWidget), intent(in) :: parent
        type(QLabel) :: label
        label=QLabel(text,parent)
        call label%setObjectName('assemblyHint')
        call label%setTextFormat(Qt_PlainText)
        call label%setWordWrap(.true.)
    end function hint

    subroutine name_control(widget,id,label)
        class(QWidget), intent(in) :: widget
        character(*), intent(in) :: id,label
        call widget%setObjectName(id)
        call widget%setAccessibleName(label)
        call widget%setMinimumHeight(30)
    end subroutine name_control

    subroutine build_member_page(editor,kind)
        type(multipatch_editor), intent(inout), target :: editor
        integer, intent(in) :: kind
        type(QScrollArea) :: scroller
        type(QWidget) :: body
        type(QVBoxLayout) :: layout
        character(:), allocatable :: patch_name
        integer :: i

        scroller=QScrollArea(editor%member_stack)
        call scroller%setProperty('frameShape',0)
        call scroller%setWidgetResizable(.true.)
        call scroller%setHorizontalScrollBarPolicy(Qt_ScrollBarAlwaysOff)
        body=QWidget(scroller)
        call body%setObjectName('assemblyMemberList')
        layout=QVBoxLayout(body)
        call layout%setContentsMargins(1,0,2,0)
        call layout%setSpacing(4)
        allocate(editor%pages(kind)%checks(single_patch_count(kind,catalog=editor%catalog)))
        do i=1,size(editor%pages(kind)%checks)
            patch_name=single_patch_name(kind,i,catalog=editor%catalog)
            if(i<editor%first_variant)patch_name=''
            associate(check=>editor%pages(kind)%checks(i))
                check=QCheckBox(patch_name,body)
                call name_control(check,'studioMultipatchMember'//integer_text(kind)//'_'//integer_text(i),patch_name)
                call check%setSizePolicy(QSizePolicy_Ignored,QSizePolicy_Fixed)
                call set_status_hint(check,patch_name)
                call check%setChecked(any(editor%drafts(kind)%members==i))
                call check%toggled%connect(member_changed,editor)
                if(len(patch_name)>0)then
                    call layout%addWidget(check)
                else
                    ! Preserve the native variant number in the checkbox array.
                    call check%hide()
                end if
            end associate
        end do
        editor%pages(kind)%empty_label=hint('',body)
        call layout%addWidget(editor%pages(kind)%empty_label)
        call layout%addStretch(1)
        call scroller%setWidget(body)
        call editor%member_stack%addWidget(scroller)
    end subroutine build_member_page

    subroutine filter_members(query,data)
        character(*), intent(in) :: query
        class(*), optional, intent(inout) :: data
        character(:), allocatable :: name,pattern
        integer :: i,shown
        logical :: matches
        select type(editor=>data)
        type is(multipatch_editor)
            pattern=fold_name(trim(adjustl(query)))
            shown=0
            associate(page=>editor%pages(editor%kind))
                do i=1,size(page%checks)
                    name=single_patch_name(editor%kind,i,catalog=editor%catalog)
                    matches=i>=editor%first_variant.and.len(name)>0.and.index(fold_name(name),pattern)>0
                    call page%checks(i)%setVisible(matches)
                    if(matches)shown=shown+1
                end do
                call page%empty_label%setText('No patches of this type in the project.')
                if(len(pattern)>0)call page%empty_label%setText('No patches match this search.')
                call page%empty_label%setVisible(shown==0)
            end associate
        end select
    end subroutine filter_members

    pure function fold_name(name) result(folded)
        character(*), intent(in) :: name
        character(len(name)) :: folded
        integer :: i,code
        folded=name
        do i=1,len(name)
            code=iachar(name(i:i))
            if(code>=iachar('A').and.code<=iachar('Z'))folded(i:i)=achar(code+iachar('a')-iachar('A'))
        end do
    end function fold_name

    subroutine member_changed(checked,data)
        logical, intent(in) :: checked
        class(*), optional, intent(inout) :: data
        integer :: i
        logical, allocatable :: members(:)
        select type(editor=>data)
        type is(multipatch_editor)
            members=[(editor%pages(editor%kind)%checks(i)%isChecked(),i=1,size(editor%pages(editor%kind)%checks))]
            editor%drafts(editor%kind)%members=pack([(i,i=1,size(members))],members)
            call dirty(editor)
        end select
    end subroutine member_changed

    subroutine name_changed(name,data)
        character(*), intent(in) :: name
        class(*), optional, intent(inout) :: data
        select type(editor=>data)
        type is(multipatch_editor)
            editor%drafts(editor%kind)%name=trim(adjustl(name))
            call dirty(editor)
        end select
    end subroutine name_changed

    subroutine change_kind(index,data)
        integer, intent(in) :: index
        class(*), optional, intent(inout) :: data
        select type(editor=>data)
        type is(multipatch_editor)
            if(index<0.or.index>=size(ranks))return
            editor%kind=index+1
            call editor%member_stack%setCurrentIndex(index)
            call editor%name%setText(editor%drafts(editor%kind)%name)
            call filter_members(editor%filter%text(),editor)
            call dirty(editor)
        end select
    end subroutine change_kind

    subroutine change_join_row(row,column,previous_row,previous_column,data)
        integer, intent(in) :: row,column,previous_row,previous_column
        class(*), optional, intent(inout) :: data
        select type(editor=>data)
        type is(multipatch_editor)
            call update_buttons(editor)
            call show_interface_detail(editor)
        end select
    end subroutine change_join_row

    subroutine double_click_join(row,column,data)
        integer, intent(in) :: row,column
        class(*), optional, intent(inout) :: data
        call edit_join(data)
    end subroutine double_click_join

    subroutine add_join(data)
        class(*), optional, intent(inout) :: data
        type(multipatch_join) :: draft
        select type(editor=>data)
        type is(multipatch_editor)
            associate(assembly=>editor%drafts(editor%kind))
                if(size(assembly%members)<2)return
                draft%patch_a=assembly%members(1)
                draft%patch_b=assembly%members(2)
                draft%side_a=2
                draft%side_b=1
                if(.not.edit_interface(editor%dialog,editor%kind,assembly%members,draft,.false.,editor%catalog))return
                assembly%joins=[assembly%joins,draft]
                call dirty(editor)
                call editor%table%setCurrentCell(size(assembly%joins)-1,0)
            end associate
        end select
    end subroutine add_join

    subroutine edit_join(data)
        class(*), optional, intent(inout) :: data
        type(multipatch_join) :: draft
        integer :: row
        select type(editor=>data)
        type is(multipatch_editor)
            row=editor%table%currentRow()+1
            associate(assembly=>editor%drafts(editor%kind))
                if(row<1.or.row>size(assembly%joins).or.size(assembly%members)<2)return
                draft=assembly%joins(row)
                if(.not.edit_interface(editor%dialog,editor%kind,assembly%members,draft,.true.,editor%catalog))return
                assembly%joins(row)=draft
                call dirty(editor)
                call editor%table%setCurrentCell(row-1,0)
            end associate
        end select
    end subroutine edit_join

    subroutine remove_join(data)
        class(*), optional, intent(inout) :: data
        integer :: row
        select type(editor=>data)
        type is(multipatch_editor)
            row=editor%table%currentRow()+1
            associate(assembly=>editor%drafts(editor%kind))
                if(row<1.or.row>size(assembly%joins))return
                assembly%joins=[assembly%joins(:row-1),assembly%joins(row+1:)]
                call dirty(editor)
                if(size(assembly%joins)>0)call editor%table%setCurrentCell(min(row,size(assembly%joins))-1,0)
            end associate
        end select
    end subroutine remove_join

    subroutine dirty(editor)
        type(multipatch_editor), intent(inout) :: editor
        editor%checked=.false.
        call refresh_editor(editor)
    end subroutine dirty

    subroutine refresh_editor(editor)
        type(multipatch_editor), intent(inout) :: editor
        integer :: row,current
        character(:), allocatable :: feedback,label,detail

        current=editor%table%currentRow()
        associate(assembly=>editor%drafts(editor%kind))
            call editor%count_label%setText(integer_text(size(assembly%members))//' selected')
            call editor%table%setRowCount(size(assembly%joins))
            do row=1,size(assembly%joins)
                associate(join=>assembly%joins(row))
                    call table_cell(editor%table,row,1,boundary_text(editor%kind,join%patch_a,join%side_a,editor%catalog))
                    call table_cell(editor%table,row,2,boundary_text(editor%kind,join%patch_b,join%side_b,editor%catalog))
                    call table_cell(editor%table,row,3,continuity_text(join),orientation_text(join))
                    call interface_status(editor,row,label,detail)
                    call table_cell(editor%table,row,4,label,detail)
                end associate
            end do
            if(current>=0.and.current<size(assembly%joins))call editor%table%setCurrentCell(current,0)
            feedback='Run Check interfaces to validate the current draft.'
            if(size(assembly%members)<2)then
                feedback='Select at least two patches to form an assembly.'
            else if(size(assembly%joins)==0)then
                feedback='You can save the assembly now or add an interface between its patches.'
            end if
        end associate
        call update_buttons(editor)
        call show_interface_detail(editor)
        if(editor%checked)then
            if(allocated(editor%report%message))feedback=editor%report%message
            call publish_window_status(editor%dialog,feedback,warning=.not.editor%report%valid)
        else
            call publish_window_status(editor%dialog,'',blocking=.true.)
            call set_status_hint(editor%dialog,feedback)
        end if
    end subroutine refresh_editor

    subroutine interface_status(editor,row,label,detail)
        type(multipatch_editor), intent(in) :: editor
        integer, intent(in) :: row
        character(:), allocatable, intent(out) :: label,detail
        associate(assembly=>editor%drafts(editor%kind))
            associate(join=>assembly%joins(row))
                if(.not.any(assembly%members==join%patch_a).or..not.any(assembly%members==join%patch_b))then
                    label='Missing'
                    detail='Select both members or remove this interface.'
                    return
                end if
            end associate
        end associate
        label='Pending'
        detail='Use Check interfaces to validate this connection.'
        if(.not.editor%checked.or..not.allocated(editor%report%joins))return
        if(row>size(editor%report%joins))return
        associate(report=>editor%report%joins(row))
            if(.not.report%recorded)then
                label='Rejected'
            else if(report%valid)then
                label='Passed'
            else
                label='Failed'
            end if
            if(allocated(report%message))detail=report%message
        end associate
    end subroutine interface_status

    subroutine show_interface_detail(editor)
        type(multipatch_editor), intent(inout) :: editor
        integer :: row,n
        character(:), allocatable :: label,detail
        row=editor%table%currentRow()+1
        n=size(editor%drafts(editor%kind)%joins)
        detail='Select an interface to see its orientation and check result.'
        if(n==0)detail='An interface pairs two '//trim(interface_types(editor%kind))//'.'
        if(row>0.and.row<=n)then
            call interface_status(editor,row,label,detail)
            associate(join=>editor%drafts(editor%kind)%joins(row))
                detail=boundary_text(editor%kind,join%patch_a,join%side_a,editor%catalog)//' ↔ '// &
                    boundary_text(editor%kind,join%patch_b,join%side_b,editor%catalog)//new_line('a')// &
                    detail//new_line('a')//orientation_text(join)
            end associate
        end if
        call editor%interface_detail%setText(detail)
    end subroutine show_interface_detail

    subroutine update_buttons(editor)
        type(multipatch_editor), intent(inout) :: editor
        logical :: can_join,selected
        associate(assembly=>editor%drafts(editor%kind))
            selected=editor%table%currentRow()>=0.and.editor%table%currentRow()<size(assembly%joins)
            can_join=size(assembly%members)>=2
            call editor%add_button%setEnabled(can_join)
            call editor%edit_button%setEnabled(can_join.and.selected)
            call editor%remove_button%setEnabled(selected)
            call editor%check_button%setEnabled(can_join.and.len_trim(assembly%name)>0)
            call editor%save_button%setEnabled(can_join.and.len_trim(assembly%name)>0)
        end associate
    end subroutine update_buttons

    subroutine check_draft(data)
        class(*), optional, intent(inout) :: data
        select type(editor=>data)
        type is(multipatch_editor)
            call inspect_multipatch(editor%drafts(editor%kind),editor%report,catalog=editor%catalog)
            editor%checked=.true.
            call refresh_editor(editor)
        end select
    end subroutine check_draft

    subroutine save_draft(data)
        class(*), optional, intent(inout) :: data
        select type(editor=>data)
        type is(multipatch_editor)
            if(allocated(editor%reserved_names))then
                if(any(editor%reserved_names==editor%drafts(editor%kind)%name))then
                    call publish_window_status(editor%dialog,'A multipatch already uses this name. Enter a different name.', &
                        warning=.true.,blocking=.true.)
                    call editor%name%setFocus()
                    return
                end if
            end if
            call check_draft(editor)
            if(editor%report%valid)then
                call editor%dialog%accept()
            else if(allocated(editor%report%message))then
                call publish_window_status(editor%dialog,editor%report%message,warning=.true.,blocking=.true.)
            end if
        end select
    end subroutine save_draft

    subroutine cancel_draft(data)
        class(*), optional, intent(inout) :: data
        select type(editor=>data)
        type is(multipatch_editor)
            call editor%dialog%reject()
        end select
    end subroutine cancel_draft

    subroutine table_cell(table,row,column,text,tooltip)
        type(QTableWidget), intent(in) :: table
        integer, intent(in) :: row,column
        character(*), intent(in) :: text
        character(*), optional, intent(in) :: tooltip
        type(QTableWidgetItem) :: cell
        cell=QTableWidgetItem(text)
        call cell%setFlags(ior(Qt_ItemIsEnabled,Qt_ItemIsSelectable))
        call set_status_hint(cell,text)
        if(present(tooltip))call set_status_hint(cell,tooltip)
        if(column==4)then
            select case(text)
            case('Passed')
                call cell%setForeground(QColor(146,223,194))
            case('Failed','Rejected','Missing')
                call cell%setForeground(QColor(255,173,162))
            end select
        end if
        call table%setItem(row-1,column-1,cell)
    end subroutine table_cell

    function side_text(side) result(label)
        integer, intent(in) :: side
        character(:), allocatable :: label
        character(3), parameter :: axes='UVW'
        integer :: axis
        label='—'
        if(side<1.or.side>6)return
        axis=(side+1)/2
        if(mod(side,2)==1)then
            label=axes(axis:axis)//' min'
        else
            label=axes(axis:axis)//' max'
        end if
    end function side_text

    function member_text(kind,patch,catalog) result(label)
        integer, intent(in) :: kind,patch
        type(studio_patch_catalog), intent(in) :: catalog
        character(:), allocatable :: label
        label=single_patch_name(kind,patch,catalog=catalog)
        if(len(label)==0)label='Unavailable patch'
    end function member_text

    function boundary_text(kind,patch,side,catalog) result(label)
        integer, intent(in) :: kind,patch,side
        type(studio_patch_catalog), intent(in) :: catalog
        character(:), allocatable :: label
        label=side_text(side)//' · '//member_text(kind,patch,catalog)
    end function boundary_text

    function orientation_text(join) result(label)
        type(multipatch_join), intent(in) :: join
        character(:), allocatable :: label
        label='Boundary orientation'
        if(join%reverse)label=label//' · Reversed'
        if(join%swap)label=label//' · Axes swapped'
        if(join%flip(1))label=label//' · Axis 1 flipped'
        if(join%flip(2))label=label//' · Axis 2 flipped'
        if(.not.(join%reverse.or.join%swap.or.any(join%flip)))label=label//' · Aligned'
    end function orientation_text

    function continuity_text(join) result(label)
        type(multipatch_join), intent(in) :: join
        character(:), allocatable :: label
        label='Free'
        if(join%continuity>=0)label=merge('G','C',join%geometric)//integer_text(join%continuity)
    end function continuity_text

    function integer_text(value) result(text)
        integer, intent(in) :: value
        character(:), allocatable :: text
        character(22) :: buffer
        write(buffer,'(i0)')value
        text=trim(buffer)
    end function integer_text

    logical function edit_interface(parent,kind,members,draft,editing,catalog) result(accepted)
        class(QWidget), intent(in) :: parent
        integer, intent(in) :: kind,members(:)
        type(multipatch_join), intent(inout) :: draft
        logical, intent(in) :: editing
        type(studio_patch_catalog), intent(in), target :: catalog
        type(interface_editor), target :: editor
        type(QVBoxLayout) :: layout
        type(QGridLayout) :: form
        type(QHBoxLayout) :: orientation,footer
        type(QLabel) :: label
        type(QDialogButtonBox) :: buttons
        type(QPushButton) :: save
        character(:), allocatable :: title
        character(1) :: endpoint
        integer :: a,i,initial(2),sides(2),index

        accepted=.false.
        editor%catalog=>catalog
        editor%members=members
        editor%kind=kind
        editor%result=draft
        initial=[draft%patch_a,draft%patch_b]
        sides=[draft%side_a,draft%side_b]
        title='Add interface'
        if(editing)title='Edit interface'
        editor%dialog=QDialog(parent)
        call editor%dialog%setObjectName('studioInterfaceDialog')
        call editor%dialog%setWindowTitle(title)
        call editor%dialog%setMinimumWidth(564)
        call editor%dialog%resize(620,320)
        call editor%dialog%setStyleSheet(dialog_style)
        layout=QVBoxLayout(editor%dialog)
        call layout%setContentsMargins(18,16,18,16)
        call layout%setSpacing(10)
        label=QLabel(title,editor%dialog)
        call label%setObjectName('studioDialogTitle')
        call layout%addWidget(label)
        call set_status_hint(editor%dialog,'Pair the '//trim(interface_types(kind))//' and choose the continuity to check.')

        form=QGridLayout()
        call form%setSpacing(8)
        call form%setColumnStretch(1,1)
        call form%setColumnStretch(2,1)
        do a=1,2
            endpoint=achar(iachar('A')+a-1)
            label=heading('Patch '//endpoint,editor%dialog)
            call form%addWidget(label,0,a)
            editor%patches(a)=QComboBox(editor%dialog)
            call name_control(editor%patches(a),'studioInterfacePatch'//endpoint,'Patch '//endpoint)
            call editor%patches(a)%setSizePolicy(QSizePolicy_Ignored,QSizePolicy_Fixed)
            index=-1
            do i=1,size(members)
                call editor%patches(a)%addItem(single_patch_name(kind,members(i),catalog=editor%catalog))
                if(initial(a)==members(i))index=i-1
            end do
            call editor%patches(a)%setCurrentIndex(index)
            call form%addWidget(editor%patches(a),1,a)
            editor%sides(a)=QComboBox(editor%dialog)
            call name_control(editor%sides(a),'studioInterfaceSide'//endpoint,'Boundary of patch '//endpoint)
            do i=1,2*kind
                call editor%sides(a)%addItem(side_text(i))
            end do
            call editor%sides(a)%setCurrentIndex(sides(a)-1)
            call form%addWidget(editor%sides(a),2,a)
        end do
        label=QLabel('Member',editor%dialog)
        call form%addWidget(label,1,0)
        label=QLabel('Boundary',editor%dialog)
        call form%addWidget(label,2,0)
        label=QLabel('Continuity',editor%dialog)
        call form%addWidget(label,3,0)
        editor%continuity=QComboBox(editor%dialog)
        call name_control(editor%continuity,'studioInterfaceContinuity','Continuity type')
        call editor%continuity%addItem('Parametric · Cⁿ')
        call editor%continuity%addItem('Geometric · Gⁿ')
        call editor%continuity%addItem('Unconstrained')
        index=merge(1,0,draft%geometric)
        if(draft%continuity<0)index=2
        call editor%continuity%setCurrentIndex(index)
        call form%addWidget(editor%continuity,3,1)
        editor%order=QSpinBox(editor%dialog)
        call name_control(editor%order,'studioInterfaceOrder','Continuity order')
        call editor%order%setPrefix('Order  ')
        call editor%order%setRange(0,huge(1))
        call editor%order%setValue(max(0,draft%continuity))
        call form%addWidget(editor%order,3,2)

        editor%jet_label=QLabel('Transition jet',editor%dialog)
        call form%addWidget(editor%jet_label,4,0)
        editor%jet=QLineEdit(parent=editor%dialog)
        call name_control(editor%jet,'studioInterfaceTransition','Normalized transition derivatives')
        call editor%jet%setPlaceholderText('Enter derivatives separated by spaces')
        call editor%jet%setMaxLength(huge(1))
        if(allocated(draft%reparameterization))call editor%jet%setText(jet_text(draft%reparameterization))
        call form%addWidget(editor%jet,4,1,1,2)
        call layout%addLayout(form)

        orientation=QHBoxLayout()
        call orientation%setSpacing(12)
        editor%reverse=QCheckBox('Reverse tangent',editor%dialog)
        if(kind==3)call editor%reverse%setText('Reverse both axes')
        call name_control(editor%reverse,'studioInterfaceReverse','Reverse patch B tangential direction')
        call set_status_hint(editor%reverse,'Reverse the orientation along the selected boundary of patch B')
        call editor%reverse%setChecked(draft%reverse)
        call editor%reverse%setVisible(kind>=2)
        call orientation%addWidget(editor%reverse)
        editor%swap=QCheckBox('Swap axes',editor%dialog)
        call name_control(editor%swap,'studioInterfaceSwap','Swap patch B tangential axes')
        call set_status_hint(editor%swap,'Exchange the two tangential axes on the face of patch B')
        call editor%swap%setChecked(draft%swap)
        call editor%swap%setVisible(kind==3)
        call orientation%addWidget(editor%swap)
        do i=1,2
            editor%flip(i)=QCheckBox('Flip axis '//integer_text(i),editor%dialog)
            call name_control(editor%flip(i),'studioInterfaceFlip'//integer_text(i), &
                'Flip patch B tangential axis '//integer_text(i))
            call set_status_hint(editor%flip(i),'Reverse patch B tangential axis '//integer_text(i)//' before swapping axes')
            call editor%flip(i)%setChecked(draft%flip(i))
            call editor%flip(i)%setVisible(kind==3)
            call orientation%addWidget(editor%flip(i))
        end do
        call orientation%addStretch(1)
        if(kind>=2)call layout%addLayout(orientation)

        buttons=QDialogButtonBox(QDialogButtonBox_Cancel,editor%dialog)
        save=buttons%addButton('Add to draft',QDialogButtonBox_AcceptRole)
        call save%setIcon(line_icon(87))
        if(editing)call save%setText('Save to draft')
        call name_control(save,'studioInterfaceSave','Confirm interface in draft')
        call studio_style_primary(save)
        call save%setMinimumWidth(108)
        call save%setDefault(.true.)
        footer=QHBoxLayout()
        editor%preview=QPushButton('Check this interface',editor%dialog)
        call editor%preview%setIcon(line_icon(89))
        call name_control(editor%preview,'studioInterfacePreview','Check this interface')
        call set_status_hint(editor%preview,'Check the selected boundaries with their current orientation and continuity')
        call footer%addWidget(editor%preview)
        call footer%addStretch(1)
        call footer%addWidget(buttons)
        call layout%addLayout(footer)
        call editor%preview%clicked%connect(preview_interface,editor)
        call buttons%accepted%connect(accept_interface,editor)
        call buttons%rejected%connect(cancel_interface,editor)
        call editor%continuity%currentIndexChanged%connect(interface_option_changed,editor)
        call editor%order%valueChanged%connect(interface_option_changed,editor)
        call editor%jet%textChanged%connect(interface_text_changed,editor)
        call editor%reverse%toggled%connect(interface_orientation_changed,editor)
        call editor%swap%toggled%connect(interface_orientation_changed,editor)
        do a=1,2
            call editor%patches(a)%currentIndexChanged%connect(interface_option_changed,editor)
            call editor%sides(a)%currentIndexChanged%connect(interface_option_changed,editor)
            call editor%flip(a)%toggled%connect(interface_orientation_changed,editor)
        end do
        call update_interface(editor)
        call editor%patches(1)%setFocus()

        if(editor%dialog%exec()==QDialog_Accepted)then
            accepted=.true.
            draft=editor%result
        end if
        call editor%dialog%delete()
    end function edit_interface

    subroutine interface_option_changed(index,data)
        integer, intent(in) :: index
        class(*), optional, intent(inout) :: data
        select type(editor=>data)
        type is(interface_editor)
            call publish_window_status(editor%dialog,'',blocking=.true.)
            call update_interface(editor)
        end select
    end subroutine interface_option_changed

    subroutine interface_text_changed(text,data)
        character(*), intent(in) :: text
        class(*), optional, intent(inout) :: data
        select type(editor=>data)
        type is(interface_editor)
            call publish_window_status(editor%dialog,'',blocking=.true.)
        end select
    end subroutine interface_text_changed

    subroutine interface_orientation_changed(checked,data)
        logical, intent(in) :: checked
        class(*), optional, intent(inout) :: data
        select type(editor=>data)
        type is(interface_editor)
            call publish_window_status(editor%dialog,'',blocking=.true.)
        end select
    end subroutine interface_orientation_changed

    subroutine update_interface(editor)
        type(interface_editor), intent(inout) :: editor
        integer :: mode,order
        logical :: show_jet
        mode=editor%continuity%currentIndex()
        order=editor%order%value()
        show_jet=mode==1.and.order>0
        call editor%order%setEnabled(mode<2)
        call editor%jet_label%setVisible(show_jet)
        call editor%jet%setVisible(show_jet)
        if(show_jet)call set_status_hint(editor%jet,'Enter '//integer_text(order)//' normalized derivatives in order.')
    end subroutine update_interface

    subroutine accept_interface(data)
        class(*), optional, intent(inout) :: data
        select type(editor=>data)
        type is(interface_editor)
            if(read_interface(editor))call editor%dialog%accept()
        end select
    end subroutine accept_interface

    subroutine preview_interface(data)
        class(*), optional, intent(inout) :: data
        type(multipatch_draft) :: draft
        type(multipatch_result) :: report
        character(:), allocatable :: message
        select type(editor=>data)
        type is(interface_editor)
            if(.not.read_interface(editor))return
            draft%kind=editor%kind
            draft%members=[editor%result%patch_a,editor%result%patch_b]
            draft%joins=[editor%result]
            call inspect_multipatch(draft,report,catalog=editor%catalog)
            message=report%message
            if(report%patches==2.and.allocated(report%joins))message=report%joins(1)%message
            call publish_window_status(editor%dialog,message,warning=.not.report%valid)
        end select
    end subroutine preview_interface

    logical function read_interface(editor) result(valid)
        type(interface_editor), intent(inout) :: editor
        integer :: a,b,mode,order,side_a,side_b
        real(c_double), allocatable :: derivatives(:)
        valid=.false.
        a=editor%patches(1)%currentIndex()+1
        b=editor%patches(2)%currentIndex()+1
        side_a=editor%sides(1)%currentIndex()+1
        side_b=editor%sides(2)%currentIndex()+1
        if(min(a,b,side_a,side_b)<1)then
            call interface_error(editor,'Choose a member and boundary for each side.')
            if(a<1)then
                call editor%patches(1)%setFocus()
            else if(b<1)then
                call editor%patches(2)%setFocus()
            else if(side_a<1)then
                call editor%sides(1)%setFocus()
            else
                call editor%sides(2)%setFocus()
            end if
            return
        end if
        if(a==b)then
            call interface_error(editor,'Choose two different patches.')
            call editor%patches(2)%setFocus()
            return
        end if
        mode=editor%continuity%currentIndex()
        order=editor%order%value()
        if(mode==1.and.order>0)then
            allocate(derivatives(order))
            if(.not.read_jet(editor%jet%text(),derivatives))then
                call interface_error(editor,'Enter exactly '//integer_text(order)// &
                    ' finite derivatives separated by spaces, such as 1.0.')
                call editor%jet%setFocus()
                return
            end if
        end if
        associate(join=>editor%result)
            join%patch_a=editor%members(a)
            join%patch_b=editor%members(b)
            join%side_a=side_a
            join%side_b=side_b
            join%continuity=merge(-1,order,mode==2)
            join%geometric=mode==1
            join%reverse=editor%reverse%isChecked()
            join%swap=editor%swap%isChecked()
            join%flip=[editor%flip(1)%isChecked(),editor%flip(2)%isChecked()]
            call move_alloc(derivatives,join%reparameterization)
        end associate
        valid=.true.
    end function read_interface

    subroutine cancel_interface(data)
        class(*), optional, intent(inout) :: data
        select type(editor=>data)
        type is(interface_editor)
            call editor%dialog%reject()
        end select
    end subroutine cancel_interface

    subroutine interface_error(editor,text)
        type(interface_editor), intent(in) :: editor
        character(*), intent(in) :: text
        call publish_window_status(editor%dialog,text,warning=.true.,blocking=.true.)
    end subroutine interface_error

    function jet_text(values) result(text)
        real(c_double), intent(in) :: values(:)
        character(:), allocatable :: text
        character(30) :: buffer
        integer :: i
        text=''
        do i=1,size(values)
            write(buffer,'(es25.16e3)')values(i)
            if(i>1)text=text//' '
            text=text//trim(adjustl(buffer))
        end do
    end function jet_text

    logical function read_jet(text,values) result(valid)
        character(*), intent(in) :: text
        real(c_double), intent(out) :: values(:)
        character(:), allocatable :: token
        integer :: first,last,i,status
        character(*), parameter :: whitespace=' '//achar(9)

        valid=.false.
        first=1
        do i=1,size(values)
            last=verify(text(first:),whitespace)
            if(last==0)return
            first=first+last-1
            last=scan(text(first:),whitespace)
            if(last==0)last=len(text)-first+2
            last=first+last-2
            token=text(first:last)
            if(verify(token,'0123456789.eEdD+-')/=0)return
            read(token,*,iostat=status,decimal='point')values(i)
            if(status/=0)return
            if(.not.ieee_is_finite(values(i)))return
            first=last+1
        end do
        valid=verify(text(first:),whitespace)==0
    end function read_jet

end module studio_multipatch_dialog
