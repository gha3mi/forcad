module studio_status_bar
    use qt, only: QFont,QFontMetrics,QHBoxLayout,QVBoxLayout,QLabel,QMainWindow,QPlainTextEdit,QKeyEvent,QResizeEvent, &
        QRect,QSize,QStatusBar,QToolButton,QtfEventFilter,QTimer,QWidget,Qt_Key_Escape,Qt_ElideRight,Qt_ElideMiddle, &
        Qt_PlainText,QSizePolicy_Ignored,QSizePolicy_Preferred,Qt_ToolButtonTextOnly,Qt_ToolButtonTextBesideIcon,Qt_StrongFocus,QScrollBar
    use studio_icons, only: line_icon,icon_close,icon_history
    use studio_presentation_text, only: integer_text
    implicit none(type,external)
    private
    public :: studio_status_t,studio_status_log

    type :: status_notice
        character(:),allocatable :: summary,detail
        character(8) :: time
        logical :: warning
    end type status_notice

    type :: studio_status_log
        private
        type(status_notice),allocatable :: entries(:)
        integer :: unread=0
    end type studio_status_log

    type :: studio_status_t
        type(QStatusBar) :: widget
        type(QToolButton) :: context,assistant,unit_button
        type(QWidget),private :: content,workspace,drawer
        type(QLabel),private :: guidance,elapsed,warning_symbol
        type(QToolButton),private :: messages,close_drawer
        type(QPlainTextEdit),private :: message_list
        type(QTimer),private :: expiry
        type(studio_status_log),private :: log
        character(:),allocatable,private :: message,hint,activity,validation,selection,brief_selection,units,default_guidance
        character(:),allocatable,private :: assistant_activity,assistant_elapsed,assistant_chat,scene_summary,editor_preview
        logical,private :: warning=.false.,painted_warning=.false.,approval=.false.,editor_warning=.false.
    contains
        procedure :: build => build_status_bar
        procedure :: showMessage => show_message
        procedure :: clearMessage => clear_message
        procedure :: currentMessage => current_message
        procedure :: setHint => set_hint
        procedure :: setActivity => set_activity
        procedure :: setValidation => set_validation
        procedure :: setEditorPreview => set_editor_preview
        procedure :: setAssistantActivity => set_assistant_activity
        procedure :: setSceneSummary => set_scene_summary
        procedure :: setContext => set_context
        procedure :: setSelection => set_selection
        procedure :: takeLog => take_log
        procedure :: restoreLog => restore_log
        procedure,private :: refresh,refresh_messages,position_messages
    end type studio_status_t

contains

    subroutine build_status_bar(self,window,workspace)
        class(studio_status_t),intent(inout),target :: self
        type(QMainWindow),intent(in) :: window
        type(QWidget),intent(in) :: workspace
        type(QHBoxLayout) :: layout,heading
        type(QVBoxLayout) :: drawer_layout
        type(QLabel) :: title
        type(QtfEventFilter) :: events,drawer_events,workspace_events
        type(QSize) :: content_size

        self%workspace=workspace
        self%widget=window%statusBar()
        call self%widget%setObjectName('studioStatusBar')
        self%content=QWidget(self%widget)
        call self%content%setObjectName('studioStatusContent')
        layout=QHBoxLayout(self%content)
        call layout%setContentsMargins(10,2,4,2)
        call layout%setSpacing(4)
        self%warning_symbol=QLabel('⚠')
        call self%warning_symbol%setObjectName('studioStatusWarning')
        call self%warning_symbol%setAccessibleName('Warning')
        call layout%addWidget(self%warning_symbol)
        call self%warning_symbol%hide()
        self%guidance=QLabel('')
        call self%guidance%setObjectName('studioStatusGuidance')
        call self%guidance%setTextFormat(Qt_PlainText)
        call self%guidance%setAccessibleName('Workspace feedback')
        call self%guidance%setSizePolicy(QSizePolicy_Ignored,QSizePolicy_Preferred)
        call layout%addWidget(self%guidance,1)

        self%context=status_button(self%content,'studioStatusContext','Selection details')
        self%assistant=status_button(self%content,'studioStatusAssistant','Assistant activity')
        self%unit_button=status_button(self%content,'studioStatusUnits','Project units')
        self%messages=status_button(self%content,'studioStatusMessages','Messages')
        call self%messages%setToolButtonStyle(Qt_ToolButtonTextBesideIcon)
        call self%messages%setIcon(line_icon(icon_history))
        call self%messages%setIconSize(15,15)
        call self%messages%setMinimumWidth(38)
        call self%messages%clicked%connect(toggle_messages,self)
        self%elapsed=QLabel('')
        call self%elapsed%setObjectName('studioAssistantElapsed')
        call self%elapsed%setTextFormat(Qt_PlainText)
        call self%elapsed%setAccessibleName('Elapsed Assistant time')
        call layout%addWidget(self%context)
        call layout%addWidget(self%assistant)
        call layout%addWidget(self%elapsed)
        call layout%addWidget(self%unit_button)
        call layout%addWidget(self%messages)

        self%drawer=QWidget(workspace)
        call self%drawer%setObjectName('studioMessagesDrawer')
        call self%drawer%setAccessibleName('Project messages')
        drawer_layout=QVBoxLayout(self%drawer)
        call drawer_layout%setContentsMargins(10,6,10,10)
        call drawer_layout%setSpacing(6)
        heading=QHBoxLayout()
        title=QLabel('Messages',self%drawer)
        call title%setObjectName('studioMessagesTitle')
        call heading%addWidget(title,1)
        self%close_drawer=QToolButton(self%drawer)
        call self%close_drawer%setIcon(line_icon(icon_close))
        call self%close_drawer%setFocusPolicy(Qt_StrongFocus)
        call self%close_drawer%setAccessibleName('Close Messages')
        call self%close_drawer%clicked%connect(close_messages_slot,self)
        call heading%addWidget(self%close_drawer)
        call drawer_layout%addLayout(heading)
        self%message_list=QPlainTextEdit('',self%drawer)
        call self%message_list%setObjectName('studioMessagesList')
        call self%message_list%setReadOnly(.true.)
        call self%message_list%setAccessibleName('Recent messages for the active project')
        call drawer_layout%addWidget(self%message_list,1)
        call self%drawer%hide()
        drawer_events=QtfEventFilter(self%drawer)
        call drawer_events%keyPressEvent%connect(on_drawer_key,self)
        drawer_events=QtfEventFilter(self%message_list)
        call drawer_events%keyPressEvent%connect(on_drawer_key,self)
        drawer_events=QtfEventFilter(self%close_drawer)
        call drawer_events%keyPressEvent%connect(on_drawer_key,self)

        self%message=''
        self%hint=''
        self%activity=''
        self%validation=''
        self%editor_preview=''
        self%assistant_activity=''
        self%assistant_elapsed=''
        self%assistant_chat=''
        self%selection=''
        self%brief_selection=''
        self%scene_summary=''
        self%units=''
        self%default_guidance='Ready'
        self%expiry=QTimer(self%widget)
        call self%expiry%setSingleShot(.true.)
        call self%expiry%timeout%connect(on_expiry,self)
        events=QtfEventFilter(self%widget)
        call events%resizeEvent%connect(on_resize,self)
        call events%showEvent%connect(on_status_show,self)
        workspace_events=QtfEventFilter(workspace)
        call workspace_events%resizeEvent%connect(on_workspace_resize,self)
        content_size=self%content%sizeHint()
        call self%widget%setFixedHeight(max(32,content_size%height))
        call self%content%show()
        call on_status_show(self)
    end subroutine build_status_bar

    function status_button(parent,name,accessible) result(button)
        type(QWidget),intent(in) :: parent
        character(*),intent(in) :: name,accessible
        type(QToolButton) :: button

        button=QToolButton(parent)
        call button%setObjectName(name)
        call button%setToolButtonStyle(Qt_ToolButtonTextOnly)
        call button%setAccessibleName(accessible)
        call button%setFocusPolicy(Qt_StrongFocus)
    end function status_button

    subroutine show_message(self,text,timeout,warning,detail,record)
        class(studio_status_t),intent(inout) :: self
        character(*),intent(in) :: text
        integer,intent(in),optional :: timeout
        logical,intent(in),optional :: warning,record
        character(*),intent(in),optional :: detail
        logical :: keep
        character(:),allocatable :: full

        if(self%widget%isNull())return
        self%message=status_line(text)
        self%warning=.false.
        if(present(warning))self%warning=warning
        call self%expiry%stop()
        if(present(timeout))then
            if(timeout>0)call self%expiry%start(timeout)
        end if
        keep=.true.
        if(present(record))keep=record
        if(keep.and.len(self%message)>0)then
            full=trim(text)
            if(present(detail))full=detail
            call record_notice(self,full)
        end if
        call self%refresh()
    end subroutine show_message

    subroutine record_notice(self,detail)
        class(studio_status_t),intent(inout) :: self
        character(*),intent(in) :: detail
        type(status_notice) :: entry
        character(10) :: clock
        integer :: n

        if(.not.allocated(self%log%entries))allocate(self%log%entries(0))
        n=size(self%log%entries)
        if(n>0)then
            associate(previous=>self%log%entries(n))
                if((self%warning.eqv.previous%warning).and.self%message==previous%summary.and.detail==previous%detail)return
            end associate
        end if
        call date_and_time(time=clock)
        entry%time=clock(:2)//':'//clock(3:4)//':'//clock(5:6)
        entry%summary=self%message
        entry%detail=detail
        entry%warning=self%warning
        self%log%entries=[self%log%entries,entry]
        if(self%warning.and..not.self%drawer%isVisible())self%log%unread=self%log%unread+1
        if(self%drawer%isVisible())then
            call refresh_messages(self)
            call self%position_messages()
        end if
    end subroutine record_notice

    subroutine clear_message(self)
        class(studio_status_t),intent(inout) :: self
        call self%showMessage('')
    end subroutine clear_message

    function current_message(self) result(message)
        class(studio_status_t),intent(in) :: self
        character(:),allocatable :: message

        message=''
        if(allocated(self%message))message=self%message
        if(allocated(self%editor_preview))then
            if(len(self%editor_preview)>0)message=self%editor_preview
        end if
        if(allocated(self%validation))then
            if(len(self%validation)>0)message=self%validation
        end if
    end function current_message

    subroutine set_hint(self,text)
        class(studio_status_t),intent(inout) :: self
        character(*),intent(in) :: text

        if(self%widget%isNull())return
        if(text==self%hint)return
        self%hint=text
        call self%refresh()
    end subroutine set_hint

    subroutine set_activity(self,text)
        class(studio_status_t),intent(inout) :: self
        character(*),intent(in) :: text

        if(self%widget%isNull())return
        if(text==self%activity)return
        self%activity=text
        call self%refresh()
    end subroutine set_activity

    subroutine set_validation(self,text)
        class(studio_status_t),intent(inout) :: self
        character(*),intent(in) :: text

        if(self%widget%isNull())return
        if(text==self%validation)return
        self%validation=status_line(text)
        call self%refresh()
    end subroutine set_validation

    subroutine set_editor_preview(self,text,warning)
        class(studio_status_t),intent(inout) :: self
        character(*),intent(in) :: text
        logical,intent(in) :: warning

        if(self%widget%isNull())return
        if(text==self%editor_preview.and.(warning.eqv.self%editor_warning))return
        self%editor_preview=status_line(text)
        self%editor_warning=warning
        call self%refresh()
    end subroutine set_editor_preview

    subroutine set_assistant_activity(self,text,elapsed,approval,chat)
        class(studio_status_t),intent(inout) :: self
        character(*),intent(in),optional :: text,elapsed,chat
        logical,intent(in),optional :: approval

        if(self%widget%isNull())return
        if(present(text))self%assistant_activity=status_line(text)
        if(present(elapsed))self%assistant_elapsed=elapsed
        if(present(approval))then
            if(approval.neqv.self%approval)then
                self%approval=approval
                if(approval)then
                    call self%assistant%setStyleSheet('color:#ffcd89; background:#3b3026')
                else
                    call self%assistant%setStyleSheet('')
                end if
            end if
        end if
        if(present(chat))self%assistant_chat=chat
        call self%refresh()
    end subroutine set_assistant_activity

    subroutine set_scene_summary(self,text)
        class(studio_status_t),intent(inout) :: self
        character(*),intent(in) :: text

        if(self%widget%isNull())return
        if(text==self%scene_summary)return
        self%scene_summary=text
        call self%refresh()
    end subroutine set_scene_summary

    function status_line(text) result(line)
        character(*),intent(in) :: text
        character(:),allocatable :: line,part
        integer :: first,last,break

        line=''
        first=1
        do while(first<=len(text))
            break=scan(text(first:),achar(10)//achar(13))
            last=len(text)
            if(break>0)last=first+break-2
            part=trim(adjustl(text(first:last)))
            if(len(part)>0)then
                if(len(line)>0)line=line//' · '
                line=line//part
            end if
            first=last+2
        end do
    end function status_line

    subroutine set_context(self,selection,units,guidance,compact)
        class(studio_status_t),intent(inout) :: self
        character(*),intent(in) :: selection,units,guidance
        character(*),intent(in),optional :: compact

        if(self%widget%isNull())return
        self%selection=selection
        self%brief_selection=selection
        if(present(compact))self%brief_selection=compact
        self%units=units
        self%default_guidance=guidance
        call self%refresh()
    end subroutine set_context

    subroutine set_selection(self,text,compact)
        class(studio_status_t),intent(inout) :: self
        character(*),intent(in) :: text
        character(*),intent(in),optional :: compact

        if(self%widget%isNull())return
        self%selection=text
        self%brief_selection=text
        if(present(compact))self%brief_selection=compact
        call self%refresh()
    end subroutine set_selection

    integer function preferred_width(widget) result(width)
        class(QWidget),intent(in) :: widget
        type(QSize) :: size

        size=widget%sizeHint()
        width=size%width
    end function preferred_width

    subroutine refresh(self)
        class(studio_status_t),intent(inout) :: self
        type(QFont) :: font
        type(QFontMetrics) :: metrics
        character(:),allocatable :: left,selection,full_selection,activity,notice_count
        integer :: width,remaining,selection_width
        logical :: warn,busy,compact

        width=max(1,self%widget%width()-14)
        compact=width<1280
        warn=.false.
        left=self%validation
        if(len(left)>0)then
            warn=.true.
        else if(len(self%editor_preview)>0)then
            left=self%editor_preview
            warn=self%editor_warning
        else if(len(self%activity)>0)then
            left=self%activity
        else if(len(self%message)>0)then
            left=self%message
            warn=self%warning
        else if(len(self%hint)>0)then
            left=self%hint
        else
            left=self%default_guidance
        end if
        if(warn.neqv.self%painted_warning)then
            call self%guidance%setStyleSheet(merge('color:#ffcd89','color:#d4e4ef',warn))
            call self%warning_symbol%setVisible(warn)
            self%painted_warning=warn
        end if

        busy=len(self%assistant_activity)>0.or.self%approval
        activity=self%assistant_chat
        if(busy)then
            activity=self%assistant_activity
            if(len(activity)==0.or.activity=='Assistant · Working…')activity='Assistant'
            font=self%assistant%font()
            metrics=QFontMetrics(font)
            activity=metrics%elidedText(activity,Qt_ElideRight,max(110,width/6))
            call metrics%delete()
        end if
        if(self%approval)activity='Approval needed'
        call self%assistant%setText(activity)
        call self%assistant%setAccessibleDescription(self%assistant_activity)
        call self%assistant%setVisible(len(activity)>0)
        call self%elapsed%setText(self%assistant_elapsed)
        call self%elapsed%setVisible(len(self%assistant_elapsed)>0)
        call self%unit_button%setText(self%units)
        call self%unit_button%setVisible(len(self%units)>0)
        notice_count=''
        if(.not.compact)notice_count='Messages'
        if(self%log%unread>0)then
            if(len(notice_count)>0)notice_count=notice_count//' · '
            notice_count=notice_count//integer_text(self%log%unread)
        end if
        call self%messages%setText(notice_count)
        call self%messages%setAccessibleDescription(integer_text(self%log%unread)//' unread warnings')
        remaining=width-5*4
        if(warn)remaining=remaining-preferred_width(self%warning_symbol)-4
        remaining=remaining-max(38,preferred_width(self%messages))
        if(len(self%units)>0)remaining=remaining-preferred_width(self%unit_button)
        if(len(activity)>0)remaining=remaining-preferred_width(self%assistant)
        if(len(self%assistant_elapsed)>0)remaining=remaining-preferred_width(self%elapsed)

        full_selection=self%selection
        if(len(full_selection)==0.and.len(self%units)>0)full_selection=self%scene_summary
        selection=full_selection
        call self%context%setVisible(len(selection)>0)
        if(len(selection)>0)then
            font=self%context%font()
            metrics=QFontMetrics(font)
            selection_width=max(16,remaining/2)
            if(compact.and.warn)selection_width=max(16,remaining/3)
            if(compact.and.metrics%horizontalAdvance(selection)>selection_width-18.and.len(self%brief_selection)>0)then
                selection=self%brief_selection
            end if
            call self%context%setText(metrics%elidedText(selection,Qt_ElideMiddle,max(0,selection_width-18)))
            call self%context%setAccessibleDescription(full_selection)
            call self%context%setMaximumWidth(selection_width)
            remaining=remaining-min(selection_width,preferred_width(self%context))
            call metrics%delete()
        end if
        font=self%guidance%font()
        metrics=QFontMetrics(font)
        call self%guidance%setText(metrics%elidedText(left,Qt_ElideRight,max(0,remaining)))
        call self%guidance%setAccessibleDescription(left)
        call metrics%delete()
    end subroutine refresh

    subroutine take_log(self,log)
        class(studio_status_t),intent(inout) :: self
        type(studio_status_log),intent(out) :: log

        call move_alloc(self%log%entries,log%entries)
        log%unread=self%log%unread
        self%log%unread=0
        call self%setSceneSummary('')
        call self%setValidation('')
        call self%setEditorPreview('',.false.)
        call close_messages(self,restore_focus=.false.)
        call self%clearMessage()
    end subroutine take_log

    subroutine restore_log(self,log)
        class(studio_status_t),intent(inout) :: self
        type(studio_status_log),intent(inout) :: log

        call move_alloc(log%entries,self%log%entries)
        self%log%unread=log%unread
        log%unread=0
        call self%refresh()
    end subroutine restore_log

    subroutine toggle_messages(data)
        class(*),intent(inout),optional :: data

        select type(self=>data)
        type is(studio_status_t)
            if(self%drawer%isVisible())then
                call close_messages(self)
                return
            end if
            self%log%unread=0
            call self%refresh()
            call self%refresh_messages()
            call self%position_messages()
            call self%drawer%show()
            call self%drawer%raise()
            call self%close_drawer%setFocus()
        end select
    end subroutine toggle_messages

    subroutine close_messages_slot(data)
        class(*),intent(inout),optional :: data
        select type(self=>data)
        type is(studio_status_t)
            call close_messages(self)
        end select
    end subroutine close_messages_slot

    subroutine close_messages(self,restore_focus)
        class(studio_status_t),intent(inout) :: self
        logical,intent(in),optional :: restore_focus
        logical :: focus

        if(self%drawer%isNull())return
        focus=self%drawer%isVisible()
        if(present(restore_focus))focus=focus.and.restore_focus
        call self%drawer%hide()
        if(focus)call self%messages%setFocus()
    end subroutine close_messages

    subroutine refresh_messages(self)
        class(studio_status_t),intent(inout) :: self
        integer :: i
        character(:),allocatable :: text
        type(QScrollBar) :: scrollbar

        call self%message_list%clear()
        if(allocated(self%log%entries))then
            do i=size(self%log%entries),1,-1
                associate(entry=>self%log%entries(i))
                    text=entry%time//' · '
                    if(entry%warning)text=text//'Warning · '
                    text=text//entry%summary
                    if(status_line(entry%detail)/=entry%summary.or.index(entry%detail,new_line('a'))/=0)then
                        text=text//new_line('a')//entry%detail
                    end if
                    if(i>1)text=text//new_line('a')
                    call self%message_list%appendPlainText(text)
                end associate
            end do
            if(size(self%log%entries)>0)then
                scrollbar=self%message_list%verticalScrollBar()
                call scrollbar%setValue(0)
                return
            end if
        end if
        call self%message_list%setPlainText('No recent project messages.')
    end subroutine refresh_messages

    subroutine position_messages(self)
        class(studio_status_t),intent(inout) :: self
        integer :: width,height,lines
        type(QFontMetrics) :: metrics

        width=min(520,max(1,self%workspace%width()-16))
        lines=self%message_list%blockCount()
        metrics=QFontMetrics(self%message_list%font())
        height=min(max(136,(lines+2)*metrics%height()+40),264,max(1,self%workspace%height()-16))
        call metrics%delete()
        call self%drawer%setGeometry(QRect(self%workspace%width()-width-8,max(0,self%workspace%height()-height-4),width,height))
    end subroutine position_messages

    subroutine on_drawer_key(event,data)
        type(QKeyEvent),intent(in) :: event
        class(*),intent(inout),optional :: data

        if(event%key()/=Qt_Key_Escape)return
        call close_messages_slot(data)
        call event%accept()
    end subroutine on_drawer_key

    subroutine on_expiry(data)
        class(*),intent(inout),optional :: data

        select type(self=>data)
        type is(studio_status_t)
            call self%clearMessage()
        end select
    end subroutine on_expiry

    subroutine on_status_show(data)
        class(*),intent(inout),optional :: data

        select type(self=>data)
        type is(studio_status_t)
            call self%content%setGeometry(QRect(0,0,self%widget%width(),self%widget%height()))
            call self%content%raise()
            call self%refresh()
        end select
    end subroutine on_status_show

    subroutine on_resize(event,data)
        type(QResizeEvent),intent(in) :: event
        class(*),intent(inout),optional :: data

        if(event%width()<=0)return
        call on_status_show(data)
    end subroutine on_resize

    subroutine on_workspace_resize(event,data)
        type(QResizeEvent),intent(in) :: event
        class(*),intent(inout),optional :: data

        if(event%width()<=0)return
        select type(self=>data)
        type is(studio_status_t)
            if(self%drawer%isVisible())call self%position_messages()
        end select
    end subroutine on_workspace_resize
end module studio_status_bar
