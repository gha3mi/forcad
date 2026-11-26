module studio_status_hints
    use, intrinsic :: iso_c_binding, only: c_associated
    use qt, only: QAction, QApplication_activeModalWidget, QObject, QPoint, QRect, QCursor_pos, QHeaderView, &
        QMainWindow, QMenu, QStatusBar, QTabBar, QTabWidget, QTableWidget, QTableWidgetItem, QTimer, &
        QToolButton, QToolTip_hideText, QToolTip_isVisible, QTreeWidget, QTreeWidgetItem, QWidget, QtfEventFilter, &
        Qt_BacktabFocusReason, Qt_ShortcutFocusReason, Qt_TabFocusReason, qt_slot_string
    use studio_status_bar, only: studio_status_t
    implicit none(type,external)
    private

    public :: register_status_hints, register_status_window, set_status_hint, watch_status_hint, status_hint_text, &
        set_tab_status_hint, register_status_menu, prepare_tree_item_hints, publish_window_status

    character(*), parameter :: hint_property='studioStatusHint'
    character(*), parameter :: validation_property='studioStatusValidation'
    character(*), parameter :: preview_property='studioStatusEditorPreview'
    character(*), parameter :: warning_property='studioStatusEditorWarning'
    character(*), parameter :: tree_hint_offset='studioTreeHintColumnOffset'
    character(*), parameter :: watched_property='studioStatusHintWatched'
    integer, parameter :: hint_interval_ms=50

    type :: hint_site
        type(QWidget) :: window,focused
        type(QStatusBar) :: menu_status
        type(studio_status_t), pointer :: status=>null()
        type(QWidget), allocatable :: floating(:)
        type(QPoint) :: cursor
        procedure(qt_slot_string), pointer, nopass :: publish=>null()
        class(*), pointer :: target=>null()
        character(:), allocatable :: raw, plain, shown
        type(hint_site), pointer :: next=>null()
    end type hint_site

    type :: hint_focus
        type(QWidget) :: widget
        integer :: serial=0
        type(hint_focus), pointer :: next=>null()
    end type hint_focus

    type(hint_site), pointer :: sites=>null()
    type(hint_focus), pointer :: watches=>null()

    interface set_status_hint
        module procedure set_widget_hint, set_action_hint, set_cell_hint, set_tree_item_hint
    end interface set_status_hint

contains

    subroutine register_status_hints(window,publish,data,status)
        class(QWidget), intent(in) :: window
        procedure(qt_slot_string) :: publish
        class(*), optional, intent(inout), target :: data
        type(studio_status_t), optional, intent(inout), target :: status
        type(hint_site), pointer :: site
        type(QTimer) :: timer
        type(QMainWindow) :: main

        site=>owner_site(window)
        if(associated(site))return
        allocate(site)
        site%window=window
        site%cursor=QCursor_pos()
        site%raw=''
        site%plain=''
        site%shown=''
        if(window%inherits('QMainWindow'))then
            call main%attach(window%ptr)
            site%menu_status=main%statusBar()
        end if
        if(present(status))site%status=>status
        site%publish=>publish
        if(present(data))site%target=>data
        site%next=>sites
        sites=>site
        allocate(site%floating(0))
        timer=QTimer(window)
        call timer%timeout%connect(update_registered_hint,site)
        call window%destroyed%connect(unregister_hints,site)
        call timer%start(hint_interval_ms)
    end subroutine register_status_hints

    subroutine publish_window_status(window,text,warning,blocking)
        class(QWidget), intent(in) :: window
        character(*), intent(in) :: text
        logical, optional, intent(in) :: warning,blocking
        type(hint_site), pointer :: site
        type(QWidget) :: modal
        logical :: warn,validate

        warn=.false.
        if(present(warning))warn=warning
        validate=.false.
        if(present(blocking))validate=blocking
        if(validate)then
            call window%setProperty(validation_property,text)
            if(len_trim(text)==0)then
                call window%setProperty(preview_property,'')
                call window%setProperty(warning_property,.false.)
            end if
        else
            call window%setProperty(preview_property,text)
            call window%setProperty(warning_property,warn)
        end if
        site=>owner_site(window)
        if(.not.associated(site))return
        call watch_status_hint(window)
        modal=QApplication_activeModalWidget()
        if(same_widget(window,modal))then
            call update_editor_status(site,modal)
            if(.not.associated(site%status))call publish_site_hint(site,text)
        end if
    end subroutine publish_window_status

    subroutine update_editor_status(site,modal)
        type(hint_site), intent(inout) :: site
        type(QWidget), intent(in) :: modal

        if(.not.associated(site%status))return
        if(modal%isNull())then
            call site%status%setValidation('')
            call site%status%setEditorPreview('',.false.)
        else
            call site%status%setValidation(modal%stringProperty(validation_property))
            call site%status%setEditorPreview(modal%stringProperty(preview_property),modal%boolProperty(warning_property))
        end if
    end subroutine update_editor_status

    subroutine register_status_window(owner,window)
        class(QWidget), intent(in) :: owner,window
        type(hint_site), pointer :: site
        type(QWidget) :: wrapped
        integer :: i

        site=>owner_site(owner)
        if(.not.associated(site))return
        site%floating=pack(site%floating,[(.not.site%floating(i)%isNull(),i=1,size(site%floating))])
        do i=1,size(site%floating)
            if(same_widget(site%floating(i),window))return
        end do
        wrapped=window
        site%floating=[site%floating,wrapped]
    end subroutine register_status_window

    subroutine register_status_menu(menu)
        type(QMenu), intent(in) :: menu

        call menu%setProperty('toolTipsVisible',.false.)
        call watch_status_hint(menu)
    end subroutine register_status_menu

    subroutine unregister_hints(data)
        class(*), optional, intent(inout) :: data
        type(hint_site), pointer :: site,previous

        select type(removed=>data)
        type is(hint_site)
            nullify(previous)
            site=>sites
            do while(associated(site))
                if(site%window%serial==removed%window%serial)then
                    if(associated(previous))then
                        previous%next=>site%next
                    else
                        sites=>site%next
                    end if
                    deallocate(site)
                    return
                end if
                previous=>site
                site=>site%next
            end do
        end select
    end subroutine unregister_hints

    subroutine set_widget_hint(widget,text)
        class(QWidget), intent(in) :: widget
        character(*), intent(in) :: text
        character(:), allocatable :: previous, description

        previous=widget%stringProperty(hint_property)
        if(previous==text.and.len(widget%toolTip())==0)return
        description=widget%accessibleDescription()
        call widget%setProperty(hint_property,text)
        call widget%setToolTip('')
        if(len(description)==0.or.description==plain_hint(previous))call widget%setAccessibleDescription(plain_hint(text))
        if(len_trim(text)==0)return
        call watch_status_hint(widget)
    end subroutine set_widget_hint

    subroutine set_action_hint(action,text)
        type(QAction), intent(in) :: action
        character(*), intent(in) :: text

        call action%setProperty(hint_property,text)
        call action%setStatusTip(plain_hint(text))
    end subroutine set_action_hint

    subroutine set_cell_hint(item,text)
        type(QTableWidgetItem), intent(in) :: item
        character(*), intent(in) :: text

        call item%setData(text)
    end subroutine set_cell_hint

    subroutine prepare_tree_item_hints(tree)
        type(QTreeWidget), intent(in) :: tree
        type(QHeaderView) :: header
        integer :: columns, index

        columns=tree%columnCount()
        call tree%setProperty(tree_hint_offset,columns)
        call tree%setColumnCount(2*columns)
        header=tree%header()
        do index=columns,2*columns-1
            call header%hideSection(index)
        end do
    end subroutine prepare_tree_item_hints

    subroutine set_tree_item_hint(tree,item,column,text)
        type(QTreeWidget), intent(in) :: tree
        type(QTreeWidgetItem), intent(in) :: item
        integer, intent(in) :: column
        character(*), intent(in) :: text
        integer :: offset

        offset=tree%intProperty(tree_hint_offset)
        if(offset>0)call item%setData(offset+column,text)
    end subroutine set_tree_item_hint

    subroutine set_tab_status_hint(tabs,index,text)
        type(QTabWidget), intent(in) :: tabs
        integer, intent(in) :: index
        character(*), intent(in) :: text
        type(QTabBar) :: bar

        call tabs%findChild('qt_tabwidget_tabbar',bar)
        call tabs%setTabToolTip(index,'')
        if(bar%isNull())return
        call bar%setProperty(tab_hint_property(index),text)
        call watch_status_hint(bar)
    end subroutine set_tab_status_hint

    subroutine watch_status_hint(widget)
        class(QWidget), intent(in) :: widget
        type(hint_focus), pointer :: watch
        type(QtfEventFilter) :: filter

        call clear_widget_tooltip(widget)
        if(widget%boolProperty(watched_property))return
        call widget%setProperty(watched_property,.true.)
        allocate(watch)
        watch%widget=widget
        watch%serial=widget%serial
        watch%next=>watches
        watches=>watch
        filter=QtfEventFilter(widget)
        call filter%focusInEvent%connect(focus_hint,watch)
        call filter%focusOutEvent%connect(clear_focus_hint,watch)
        if(widget%inherits('QMenu'))call filter%showEvent%connect(show_menu,watch)
        if(widget%inherits('QDialog'))then
            call filter%showEvent%connect(show_editor_status,watch)
            call filter%hideEvent%connect(hide_editor_status,watch)
        end if
        if(widget%inherits('QToolButton'))call filter%enterEvent%connect(clear_watched_tooltip,watch)
        call widget%destroyed%connect(unwatch_hint,watch)
    end subroutine watch_status_hint

    subroutine clear_widget_tooltip(widget)
        class(QWidget), intent(in) :: widget

        if(len(widget%toolTip())>0)call widget%setToolTip('')
    end subroutine clear_widget_tooltip

    subroutine clear_watched_tooltip(data)
        class(*), optional, intent(inout) :: data

        select type(watch=>data)
        type is(hint_focus)
            call clear_widget_tooltip(watch%widget)
        end select
    end subroutine clear_watched_tooltip

    subroutine unwatch_hint(data)
        class(*), optional, intent(inout) :: data
        type(hint_focus), pointer :: watch,previous

        select type(removed=>data)
        type is(hint_focus)
            nullify(previous)
            watch=>watches
            do while(associated(watch))
                if(watch%serial==removed%serial)then
                    if(associated(previous))then
                        previous%next=>watch%next
                    else
                        watches=>watch%next
                    end if
                    deallocate(watch)
                    return
                end if
                previous=>watch
                watch=>watch%next
            end do
        end select
    end subroutine unwatch_hint

    subroutine show_menu(data)
        class(*), optional, intent(inout) :: data
        type(hint_site), pointer :: site

        select type(watch=>data)
        type is(hint_focus)
            site=>owner_site(watch%widget)
            if(associated(site))call register_status_window(site%window,watch%widget)
        end select
    end subroutine show_menu

    subroutine show_editor_status(data)
        class(*), optional, intent(inout) :: data
        type(hint_site), pointer :: site

        select type(watch=>data)
        type is(hint_focus)
            site=>owner_site(watch%widget)
            if(associated(site))call update_editor_status(site,watch%widget)
        end select
    end subroutine show_editor_status

    subroutine hide_editor_status(data)
        class(*), optional, intent(inout) :: data
        type(hint_site), pointer :: site
        type(QWidget) :: empty

        select type(watch=>data)
        type is(hint_focus)
            site=>owner_site(watch%widget)
            if(associated(site))call update_editor_status(site,empty)
        end select
    end subroutine hide_editor_status

    subroutine focus_hint(reason,data)
        integer, intent(in) :: reason
        class(*), optional, intent(inout) :: data
        type(hint_site), pointer :: site

        if(reason/=Qt_TabFocusReason.and.reason/=Qt_BacktabFocusReason.and.reason/=Qt_ShortcutFocusReason)return
        select type(watch=>data)
        type is(hint_focus)
            site=>owner_site(watch%widget)
            if(.not.associated(site))return
            site%focused=watch%widget
            site%cursor=QCursor_pos()
            call update_site_hint(site)
        end select
    end subroutine focus_hint

    subroutine clear_focus_hint(reason,data)
        integer, intent(in) :: reason
        class(*), optional, intent(inout) :: data
        type(hint_site), pointer :: site
        type(QWidget) :: empty

        select type(watch=>data)
        type is(hint_focus)
            site=>owner_site(watch%widget)
            if(.not.associated(site))return
            if(same_widget(site%focused,watch%widget))site%focused=empty
        end select
    end subroutine clear_focus_hint

    subroutine update_registered_hint(data)
        class(*), optional, intent(inout) :: data

        select type(site=>data)
        type is(hint_site)
            call update_site_hint(site)
        end select
    end subroutine update_registered_hint

    subroutine update_site_hint(site)
        type(hint_site), intent(inout) :: site
        type(hint_site), pointer :: modal_owner
        type(QWidget) :: modal,host,child,empty
        type(QPoint) :: position,global
        character(:), allocatable :: detail
        integer :: i

        if(.not.site%window%isVisible())return
        global=QCursor_pos()
        if(global%x/=site%cursor%x.or.global%y/=site%cursor%y)then
            site%focused=empty
            site%cursor=global
        end if
        modal=QApplication_activeModalWidget()
        if(.not.modal%isNull())then
            modal_owner=>owner_site(modal)
            if(.not.associated(modal_owner))then
                call update_editor_status(site,empty)
                call publish_site_hint(site,'')
                return
            end if
            if(.not.same_widget(modal_owner%window,site%window))then
                call update_editor_status(site,empty)
                call publish_site_hint(site,'')
                return
            end if
            host=modal
            if(.not.within_widget(site%focused,modal))site%focused=empty
        else
            host=site%window
            do i=size(site%floating),1,-1
                if(site%floating(i)%isNull())cycle
                if(.not.site%floating(i)%isVisible())cycle
                if(under_cursor(site%floating(i),global))then
                    host=site%floating(i)
                    exit
                end if
            end do
        end if

        call update_editor_status(site,modal)
        detail=''
        if(under_cursor(host,global))then
            position=host%mapFromGlobal(global)
            child=host%childAt(position)
            if(child%isNull())child=host
            call clear_widget_tooltip(child)
        end if
        if(host%inherits('QMenu').and..not.site%menu_status%isNull())then
            detail=site%menu_status%currentMessage()
        end if
        if(len(detail)==0.and..not.site%focused%isNull())then
            if(site%focused%hasFocus())detail=widget_hint(site%focused,global,keyboard=.true.,boundary=modal)
        end if
        if(.not.child%isNull())then
            if(len(detail)==0)detail=widget_hint(child,global,keyboard=.false.,boundary=modal)
        end if
        if(len(detail)==0.and..not.modal%isNull())detail=modal%stringProperty(hint_property)
        if(QToolTip_isVisible())call QToolTip_hideText()
        if(detail/=site%raw)then
            site%raw=detail
            site%plain=plain_hint(detail)
        end if
        call publish_site_hint(site,site%plain)
    end subroutine update_site_hint

    function status_hint_text(widget) result(text)
        class(QWidget), intent(in) :: widget
        character(:), allocatable :: text
        type(QWidget) :: control

        control=widget
        text=plain_hint(widget_hint(control,QCursor_pos(),keyboard=.true.))
    end function status_hint_text

    function tab_hint_property(tab) result(name)
        integer, intent(in) :: tab
        character(:), allocatable :: name
        character(16) :: number

        write(number,'(i0)')tab
        name=hint_property//trim(number)
    end function tab_hint_property

    function widget_hint(widget,global,keyboard,boundary) result(text)
        type(QWidget), intent(in) :: widget
        type(QPoint), intent(in) :: global
        logical, intent(in) :: keyboard
        type(QWidget), optional, intent(in) :: boundary
        character(:), allocatable :: text
        type(QObject) :: parent
        type(QWidget) :: current
        type(QToolButton) :: tool
        type(QAction) :: action
        type(QTabBar) :: tabs
        type(QPoint) :: pos
        integer :: tab

        current=widget
        text=''
        do while(.not.current%isNull())
            if(current%inherits('QTabBar'))then
                call tabs%attach(current%ptr)
                if(keyboard)then
                    tab=tabs%currentIndex()
                else
                    pos=current%mapFromGlobal(global)
                    tab=tabs%tabAt(pos%x,pos%y)
                end if
                if(tab>=0)then
                    text=tabs%stringProperty(tab_hint_property(tab))
                end if
            end if
            if(len(text)==0)text=item_view_hint(current,global,keyboard)
            if(len(text)==0)text=current%stringProperty(hint_property)
            if(len(text)==0.and.current%inherits('QToolButton'))then
                call tool%attach(current%ptr)
                action=tool%defaultAction()
                if(.not.action%isNull())text=action%stringProperty(hint_property)
            end if
            call clear_widget_tooltip(current)
            if(len(text)>0)return
            if(present(boundary))then
                if(same_widget(current,boundary))return
            end if
            parent=current%parent()
            if(parent%isNull())return
            if(.not.parent%inherits('QWidget'))return
            call current%attach(parent%ptr)
        end do
    end function widget_hint

    function item_view_hint(view,global,keyboard) result(text)
        type(QWidget), intent(in) :: view
        type(QPoint), intent(in) :: global
        logical, intent(in) :: keyboard
        character(:), allocatable :: text
        type(QTableWidget) :: table
        type(QTableWidgetItem) :: cell
        type(QTreeWidget) :: tree
        type(QTreeWidgetItem) :: item
        type(QHeaderView) :: horizontal,vertical
        type(QWidget) :: viewport
        type(QPoint) :: position
        integer :: row,column,offset

        text=''
        if(view%inherits('QTableWidget'))then
            call table%attach(view%ptr)
            if(keyboard)then
                row=table%currentRow()
                column=table%currentColumn()
            else
                viewport=table%viewport()
                if(.not.under_cursor(viewport,global))return
                position=viewport%mapFromGlobal(global)
                vertical=table%verticalHeader()
                horizontal=table%horizontalHeader()
                row=vertical%logicalIndexAt(position%y)
                column=horizontal%logicalIndexAt(position%x)
            end if
            if(row<0.or.column<0)return
            cell=table%item(row,column)
            if(.not.cell%isNull())text=cell%stringData()
        else if(view%inherits('QTreeWidget'))then
            offset=view%intProperty(tree_hint_offset)
            if(offset<1)return
            call tree%attach(view%ptr)
            column=0
            if(keyboard)then
                item=tree%currentItem()
            else
                viewport=tree%viewport()
                if(.not.under_cursor(viewport,global))return
                position=viewport%mapFromGlobal(global)
                horizontal=tree%header()
                column=horizontal%logicalIndexAt(position%x)
                if(column<0)return
                item=tree%itemAt(position)
            end if
            if(item%isNull())return
            text=item%stringData(offset+column)
            if(len(text)==0.and.column/=0)text=item%stringData(offset)
        end if
    end function item_view_hint

    subroutine publish_site_hint(site,text)
        type(hint_site), intent(inout) :: site
        character(*), intent(in) :: text

        if(text==site%shown)return
        site%shown=text
        if(associated(site%target))then
            call site%publish(text,site%target)
        else
            call site%publish(text)
        end if
    end subroutine publish_site_hint

    function owner_site(object) result(site)
        class(QObject), intent(in) :: object
        type(hint_site), pointer :: site
        type(QObject) :: parent,node
        integer :: i

        node=object
        do while(.not.node%isNull())
            site=>sites
            do while(associated(site))
                if(same_widget(node,site%window))return
                do i=1,size(site%floating)
                    if(same_widget(node,site%floating(i)))return
                end do
                site=>site%next
            end do
            parent=node%parent()
            node=parent
        end do
        nullify(site)
    end function owner_site

    logical function within_widget(widget,ancestor) result(inside)
        type(QWidget), intent(in) :: widget,ancestor
        type(QObject) :: current

        inside=.false.
        current=widget%QObject
        do while(.not.current%isNull())
            if(same_widget(current,ancestor))then
                inside=.true.
                return
            end if
            current=current%parent()
        end do
    end function within_widget

    logical function under_cursor(widget,global) result(hit)
        type(QWidget), intent(in) :: widget
        type(QPoint), intent(in) :: global
        type(QRect) :: bounds
        type(QPoint) :: local

        bounds=widget%rect()
        local=widget%mapFromGlobal(global)
        hit=bounds%contains(local)
    end function under_cursor

    logical function same_widget(first,second) result(same)
        class(QObject), intent(in) :: first,second

        same=first%serial==second%serial.and.c_associated(first%ptr,second%ptr)
    end function same_widget

    function plain_hint(text) result(plain)
        character(*), intent(in) :: text
        character(:), allocatable :: plain,line
        integer :: i,separator

        plain=''
        i=1
        do while(i<=len(text))
            separator=scan(text(i:),achar(10)//achar(13)//achar(9))
            if(separator==0)separator=len(text)-i+2
            line=trim(adjustl(text(i:i+separator-2)))
            if(len(line)>0)then
                if(len(plain)>0)plain=plain//' · '
                plain=plain//line
            end if
            i=i+separator
        end do
    end function plain_hint

end module studio_status_hints
