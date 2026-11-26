module studio_assistant_transcript
    use qt, only: QClipboard, QCoreApplication_sendEvent, QFontMetrics, QFrame_NoFrame, QGuiApplication_clipboard, QHBoxLayout, QLabel, &
        QIcon, QPushButton, QResizeEvent, QSize, &
        QScrollArea, QScrollBar, QSizePolicy_Expanding, QSizePolicy_Ignored, QSizePolicy_Preferred, QTextEdit, QTimer, QVBoxLayout, QWidget, &
        QWheelEvent, Qt_PlainText, Qt_ScrollBarAlwaysOff, QtfEventFilter
    use studio_icons, only: line_icon, icon_copy
    use studio_presentation_text, only: without_nul, studio_product_name
    use studio_theme, only: studio_divider, studio_secondary_text, studio_selected_fill, &
        studio_active_edge
    use studio_status_hints, only: set_status_hint
    implicit none(type, external)
    private

    public :: assistant_transcript_t, chat_body

    character, parameter :: lf=new_line('a')
    character(*), parameter :: gap=lf//lf, storage_header='forcad-timeline-v1'//lf, digits='0123456789'
    character(*), parameter :: agent_storage_header='forcad-timeline-v2'//lf
    integer, parameter :: agent_length_size=2, max_agent_label=48
    integer, parameter :: timestamp_size=19, length_size=9, record_size=1+timestamp_size+length_size
    integer, parameter :: measurement_height=20, timeline_margin=3, bubble_padding=10
    integer, parameter :: max_tool_summary=320, max_event_detail=8000, visible_detail_size=700

    type :: timeline_entry
        character(:), allocatable :: role, body, agent
        character(timestamp_size) :: timestamp=''
        type(timeline_entry), pointer :: next=>null()
        type(QTextEdit) :: editor
        type(QPushButton) :: disclosure
        type(QWidget) :: card, scroll_viewport
        type(QScrollBar) :: measurement
        type(QTimer) :: measurement_timer
        integer :: measured_width=0, preferred_width=0, card_width=0
        logical :: needs_measure=.true., needs_render=.false., measuring=.false., collapsed=.false.
    end type timeline_entry

    type :: assistant_transcript_t
        private
        type(QScrollArea) :: view
        type(QWidget) :: host
        type(QVBoxLayout) :: cards
        type(QLabel) :: empty_state
        type(QScrollBar) :: scrollbar
        type(QTimer) :: tail_timer
        type(QPushButton) :: jump_button
        type(QIcon) :: copy_icon
        type(timeline_entry), pointer :: first=>null(), last=>null(), active_reply=>null()
        character(:), allocatable :: empty_message
        integer :: scroll_value=0
        logical :: moving_scrollbar=.false., layout_pending=.false., follow_tail=.true.
    contains
        procedure :: attach
        procedure :: release
        procedure :: set_empty_message
        procedure :: append_message
        procedure :: begin_reply
        procedure :: clear_reply
        procedure :: update_reply
        procedure :: append_reply_chunk
        procedure :: append_activity
        procedure :: append_reasoning
        procedure :: append_plan
        procedure :: append_tool
        procedure :: serialized
        procedure :: plain_text
        procedure :: restore
        procedure :: defer_update
        procedure :: at_end
        procedure :: follow
        procedure :: finish_follow
    end type assistant_transcript_t

contains

    subroutine attach(transcript,view,scrollbar,tail_timer,jump_button)
        class(assistant_transcript_t), intent(inout), target :: transcript
        type(QScrollArea), intent(in) :: view
        type(QScrollBar), intent(in) :: scrollbar
        type(QTimer), intent(in) :: tail_timer
        type(QPushButton), intent(in), optional :: jump_button
        type(QPushButton) :: empty_button
        type(QtfEventFilter) :: resize_events

        transcript%view=view
        transcript%scrollbar=scrollbar
        transcript%tail_timer=tail_timer
        transcript%jump_button=empty_button
        if(present(jump_button))transcript%jump_button=jump_button
        call transcript%view%setWidgetResizable(.true.)
        call transcript%view%setHorizontalScrollBarPolicy(Qt_ScrollBarAlwaysOff)
        call transcript%scrollbar%valueChanged%connect(on_scroll_changed,transcript)
        call transcript%scrollbar%rangeChanged%connect(on_range_changed,transcript)
        resize_events=QtfEventFilter(view%viewport())
        call resize_events%resizeEvent%connect(on_resize,transcript)
        call transcript%restore('')
    end subroutine attach

    subroutine release(transcript)
        class(assistant_transcript_t), intent(inout) :: transcript
        type(timeline_entry), pointer :: next,entry

        transcript%layout_pending=.false.
        transcript%follow_tail=.true.
        transcript%scroll_value=0
        transcript%moving_scrollbar=.true.
        entry=>transcript%first
        do while(associated(entry))
            entry%needs_render=.false.
            entry%measuring=.true.
            entry=>entry%next
        end do
        if(.not.transcript%host%isNull())call transcript%host%delete()
        if(.not.transcript%copy_icon%isNull())call transcript%copy_icon%delete()
        nullify(transcript%active_reply,transcript%last)
        do while(associated(transcript%first))
            next=>transcript%first%next
            deallocate(transcript%first)
            transcript%first=>next
        end do
        transcript%moving_scrollbar=.false.
    end subroutine release

    subroutine build_host(transcript)
        class(assistant_transcript_t), intent(inout) :: transcript
        type(QWidget) :: timeline
        type(QVBoxLayout) :: stack

        transcript%copy_icon=line_icon(icon_copy)
        transcript%host=QWidget()
        call transcript%host%setObjectName('studioAssistantTimeline')
        call transcript%host%setStyleSheet('QWidget#studioAssistantTimeline {background:#16212a;}')
        stack=QVBoxLayout(transcript%host)
        call stack%setContentsMargins(timeline_margin,6,timeline_margin,8)
        call stack%setSpacing(0)
        timeline=QWidget()
        transcript%cards=QVBoxLayout(timeline)
        call transcript%cards%setContentsMargins(0,0,0,0)
        call transcript%cards%setSpacing(8)
        if(.not.allocated(transcript%empty_message))transcript%empty_message='Choose an agent to start a chat.'
        transcript%empty_state=QLabel(transcript%empty_message)
        call transcript%empty_state%setObjectName('studioAssistantEmptyConversation')
        call transcript%empty_state%setTextFormat(Qt_PlainText)
        call transcript%empty_state%setWordWrap(.true.)
        call transcript%empty_state%setStyleSheet('padding:8px; color:'//studio_secondary_text//'; font-size:14px;')
        call transcript%cards%addWidget(transcript%empty_state)
        call stack%addWidget(timeline)
        call stack%addStretch(1)
        call transcript%view%setWidget(transcript%host)
    end subroutine build_host

    subroutine set_empty_message(transcript,message)
        class(assistant_transcript_t), intent(inout) :: transcript
        character(*), intent(in) :: message

        transcript%empty_message=message
        if(.not.transcript%empty_state%isNull())call transcript%empty_state%setText(message)
    end subroutine set_empty_message

    subroutine append_message(transcript,role,message,agent)
        class(assistant_transcript_t), intent(inout) :: transcript
        character(*), intent(in) :: role,message
        character(*), intent(in), optional :: agent
        character(:), allocatable :: body

        call transcript%clear_reply()
        body=chat_body(message)
        if(len(body)>0)call append_entry(transcript,role,body,local_timestamp(),agent_label(agent))
    end subroutine append_message

    subroutine begin_reply(transcript)
        class(assistant_transcript_t), intent(inout) :: transcript

        call transcript%clear_reply()
        call transcript%follow()
    end subroutine begin_reply

    subroutine clear_reply(transcript)
        class(assistant_transcript_t), intent(inout) :: transcript

        nullify(transcript%active_reply)
    end subroutine clear_reply

    subroutine update_reply(transcript,text)
        class(assistant_transcript_t), intent(inout) :: transcript
        character(*), intent(in) :: text
        character(:), allocatable :: body

        body=chat_body(text)
        if(len(body)==0)return
        if(.not.associated(transcript%active_reply))then
            call append_entry(transcript,'Assistant',body,local_timestamp())
            transcript%active_reply=>transcript%last
        else if(body/=transcript%active_reply%body)then
            call queue_layout(transcript)
            transcript%active_reply%body=body
            call render_entry(transcript%active_reply)
        end if
    end subroutine update_reply

    subroutine append_reply_chunk(transcript,chunk,agent)
        class(assistant_transcript_t), intent(inout) :: transcript
        character(*), intent(in) :: chunk
        character(*), intent(in), optional :: agent
        character(:), allocatable :: body,label

        label=agent_label(agent)
        if(len(chunk)==0)return
        body=without_nul(chunk)
        if(associated(transcript%active_reply,transcript%last))then
            if(transcript%last%agent==label)then
                call queue_layout(transcript)
                transcript%active_reply%body=transcript%active_reply%body//body
                call render_entry(transcript%active_reply)
                return
            end if
        end if
        if(len(chat_body(body))>0)then
            call append_entry(transcript,'Assistant',body,local_timestamp(),label)
            transcript%active_reply=>transcript%last
        end if
    end subroutine append_reply_chunk

    subroutine append_activity(transcript,summary,agent)
        class(assistant_transcript_t), intent(inout) :: transcript
        character(*), intent(in) :: summary
        character(*), intent(in), optional :: agent
        character(timestamp_size) :: stamp
        character(:), allocatable :: item,label

        label=agent_label(agent)
        item=chat_body(summary)
        if(len(item)==0)return
        stamp=local_timestamp()
        item=stamp(12:19)//'  ·  '//item
        if(associated(transcript%last))then
            if(transcript%last%role=='Activity'.and.transcript%last%agent==label.and. &
                transcript%last%timestamp(:10)==stamp(:10))then
                call queue_layout(transcript)
                transcript%last%body=transcript%last%body//lf//item
                call render_entry(transcript%last)
                return
            end if
        end if
        call append_entry(transcript,'Activity',item,stamp,label)
    end subroutine append_activity

    subroutine append_reasoning(transcript,chunk,agent)
        class(assistant_transcript_t), intent(inout) :: transcript
        character(*), intent(in) :: chunk
        character(*), intent(in), optional :: agent
        character(:), allocatable :: label

        label=agent_label(agent)
        if(len(chunk)==0)return
        if(associated(transcript%last))then
            if(transcript%last%role=='Reasoning'.and.transcript%last%agent==label)then
                call queue_layout(transcript)
                transcript%last%body=transcript%last%body//chunk
                call render_entry(transcript%last)
                return
            end if
        end if
        if(len(chat_body(chunk))>0)call append_entry(transcript,'Reasoning',chunk,local_timestamp(),label)
    end subroutine append_reasoning

    subroutine append_plan(transcript,body,agent)
        class(assistant_transcript_t), intent(inout) :: transcript
        character(*), intent(in) :: body
        character(*), intent(in), optional :: agent
        character(:), allocatable :: detail

        detail=chat_body(body)
        if(len(detail)>0)call append_entry(transcript,'Plan',bounded_text(detail,max_event_detail),local_timestamp(),agent_label(agent))
    end subroutine append_plan

    subroutine append_tool(transcript,summary,detail,agent)
        class(assistant_transcript_t), intent(inout) :: transcript
        character(*), intent(in) :: summary
        character(*), intent(in), optional :: detail,agent
        character(:), allocatable :: body,output
        integer :: i

        body=without_nul(summary)
        do i=1,len(body)
            if(iachar(body(i:i))<32)body(i:i)=' '
        end do
        body=trim(adjustl(body))
        if(len(body)==0)return
        body=bounded_text(body,max_tool_summary)
        if(present(detail))then
            output=chat_body(detail)
            if(len(output)>0)body=body//lf//bounded_text(output,max_event_detail)
        end if
        call append_entry(transcript,'Tool',body,local_timestamp(),agent_label(agent))
    end subroutine append_tool

    function bounded_text(value,maximum) result(text)
        character(*), intent(in) :: value
        integer, intent(in) :: maximum
        character(:), allocatable :: text
        integer :: last

        if(len(value)<=maximum)then
            text=value
            return
        end if
        last=maximum+1
        do while(last>1)
            if(iand(iachar(value(last:last)),192)/=128)exit
            last=last-1
        end do
        text=value(:last-1)//' …'
    end function bounded_text

    function agent_label(agent) result(label)
        character(*), intent(in), optional :: agent
        character(:), allocatable :: label
        integer :: i

        label=''
        if(.not.present(agent))return
        label=agent
        do i=1,len(label)
            if(iachar(label(i:i))<32.or.iachar(label(i:i))==127)label(i:i)=' '
        end do
        label=trim(adjustl(label))
        if(len(label)>max_agent_label)label=bounded_text(label,max_agent_label-4)
    end function agent_label

    function timeline_title(entry) result(title)
        type(timeline_entry), intent(in) :: entry
        character(:), allocatable :: title

        title=entry%role
        if(title=='Studio')title=studio_product_name
        if(len(entry%agent)>0)then
            title=entry%agent//' · '//title
        else if(title=='Activity')then
            title=studio_product_name//' activity'
        end if
    end function timeline_title

    subroutine append_entry(transcript,role,body,stamp,agent)
        class(assistant_transcript_t), intent(inout) :: transcript
        character(*), intent(in) :: role,body
        character(timestamp_size), intent(in) :: stamp
        character(*), intent(in), optional :: agent
        type(timeline_entry), pointer :: entry

        call queue_layout(transcript)
        allocate(entry)
        entry%role=role
        entry%body=body
        entry%agent=''
        if(present(agent))entry%agent=agent
        entry%timestamp=stamp
        entry%scroll_viewport=transcript%view%viewport()
        entry%measurement_timer=transcript%tail_timer
        call make_card(transcript%cards,entry,transcript%copy_icon)
        if(associated(transcript%last))then
            transcript%last%next=>entry
        else
            transcript%first=>entry
        end if
        transcript%last=>entry
        call transcript%empty_state%hide()
    end subroutine append_entry

    subroutine make_card(cards,entry,copy_icon)
        type(QVBoxLayout), intent(inout) :: cards
        type(timeline_entry), target, intent(inout) :: entry
        type(QIcon), intent(in) :: copy_icon
        type(QWidget) :: row, metadata
        type(QHBoxLayout) :: row_layout,header
        type(QVBoxLayout) :: inside
        type(QLabel) :: role_label,clock,summary_label
        type(QPushButton) :: copy
        type(QtfEventFilter) :: selection_events,wheel_events
        type(QSize) :: metadata_size
        character(:), allocatable :: bubble_style,accent,text_style,title,summary
        integer :: detail_start
        logical :: expandable
        character(*), parameter :: action_style='QPushButton {background:transparent;color:#c7d7e2;border:1px solid transparent;'// &
            'padding:0 3px;min-width:0;font-size:11px;} QPushButton:hover {background:#314550;color:#edf3f7;}'// &
            'QPushButton:focus {background:#314550;border-color:'//studio_active_edge//';color:#edf3f7;}'

        select case(entry%role)
        case('You')
            bubble_style='background:'//studio_selected_fill//';border:1px solid #3c5c6a;'
            accent='#bcdfeb'
            text_style='color:#edf3f7;font-size:10.5pt;'
        case('Assistant')
            bubble_style='background:#16212a;border:0;'
            accent='#a1e2df'
            text_style='color:#edf3f7;font-size:14px;'
        case('Activity')
            bubble_style='background:#1b2833;border:1px solid '//studio_divider//';border-left:2px solid #527d86;'
            accent='#b1ced9'
            text_style='color:#b6ccd9;font-size:12px;'
        case('Reasoning')
            bubble_style='background:#1b2833;border:1px solid '//studio_divider//';border-left:2px solid #676782;'
            accent='#c2bfd8'
            text_style='color:#bdc9d7;font-size:13px;'
        case('Plan')
            bubble_style='background:#1b2833;border:1px solid '//studio_divider//';border-left:2px solid #5d7fa9;'
            accent='#b9d7fa'
            text_style='color:#d4e2ef;font-size:13px;'
        case('Tool')
            bubble_style='background:#1b2833;border:1px solid '//studio_divider//';border-left:2px solid #8b784e;'
            accent='#ebd6ac'
            text_style='color:#d2dce7;font-size:12px;'
        case default
            bubble_style='background:#1b2833;border:1px solid '//studio_divider//';'
            accent='#b1ced9'
            text_style='color:#c8d9e5;font-size:13px;'
        end select
        title=timeline_title(entry)
        detail_start=0
        expandable=entry%role=='Reasoning'.or.entry%role=='Plan'
        if(entry%role=='Tool')then
            detail_start=index(entry%body,lf)
            expandable=detail_start>0
            entry%collapsed=.not.expandable.or.len(entry%body)-detail_start>visible_detail_size
        end if
        row=QWidget()
        row_layout=QHBoxLayout(row)
        call row_layout%setContentsMargins(0,0,0,0)
        call row_layout%setSpacing(0)
        entry%card=QWidget()
        call entry%card%setObjectName('assistantTimelineCard')
        call entry%card%setStyleSheet('QWidget#assistantTimelineCard {'//bubble_style//'}')
        if(entry%role=='You')call row_layout%addStretch(1)
        call row_layout%addWidget(entry%card)
        if(entry%role/='You')call row_layout%addStretch(1)
        inside=QVBoxLayout(entry%card)
        call inside%setContentsMargins(bubble_padding,5,bubble_padding,7)
        call inside%setSpacing(2)
        metadata=QWidget(entry%card)
        header=QHBoxLayout(metadata)
        call header%setContentsMargins(4,0,0,0)
        call header%setSpacing(6)
        role_label=QLabel(title)
        call role_label%setObjectName('assistantTimelineRole')
        call role_label%setTextFormat(Qt_PlainText)
        call role_label%setStyleSheet('background:transparent;color:'//accent//';font-size:12px;font-weight:600;')
        if(len(entry%agent)>0)then
            call role_label%setWordWrap(.true.)
            call role_label%setSizePolicy(QSizePolicy_Ignored,QSizePolicy_Preferred)
            call header%addWidget(role_label,1)
        else
            call header%addWidget(role_label)
            call header%addStretch(1)
        end if
        if(len_trim(entry%timestamp)>0)then
            clock=QLabel(entry%timestamp(12:16))
            call clock%setObjectName('assistantTimelineTime')
            call clock%setTextFormat(Qt_PlainText)
            call clock%setAccessibleName('Local time '//entry%timestamp)
            call set_status_hint(clock,entry%timestamp)
            call clock%setStyleSheet('background:transparent;color:'//studio_secondary_text//';font-size:11px;')
            call header%addWidget(clock)
        end if
        copy=QPushButton('')
        call copy%setObjectName('assistantTimelineCopy'//trim(entry%role))
        call copy%setAccessibleName('Copy text')
        call copy%setAccessibleDescription(title)
        call set_status_hint(copy,'Copy text')
        call copy%setIcon(copy_icon)
        call copy%setIconSize(14,14)
        call copy%setFixedSize(24,21)
        call copy%setStyleSheet(action_style)
        call copy%clicked%connect(copy_card_text,entry)
        call header%addWidget(copy)
        if(expandable)then
            entry%disclosure=QPushButton('')
            call entry%disclosure%setObjectName('assistant'//entry%role//'Disclosure')
            call update_disclosure(entry)
            call entry%disclosure%setFixedHeight(21)
            call entry%disclosure%setStyleSheet(action_style)
            call entry%disclosure%clicked%connect(toggle_detail,entry)
            call header%addWidget(entry%disclosure)
        end if
        call inside%addWidget(metadata)

        if(entry%role=='Tool')then
            summary=entry%body
            if(detail_start>0)summary=summary(:detail_start-1)
            summary_label=QLabel(summary)
            call summary_label%setObjectName('assistantToolSummary')
            call summary_label%setTextFormat(Qt_PlainText)
            call summary_label%setWordWrap(.true.)
            call summary_label%setSizePolicy(QSizePolicy_Ignored,QSizePolicy_Preferred)
            call summary_label%setStyleSheet('background:transparent;padding:1px 4px;color:#e2ebf4;font-size:13px;')
            call inside%addWidget(summary_label)
        end if

        entry%editor=QTextEdit()
        call entry%editor%setObjectName('assistantTimeline'//trim(entry%role))
        call entry%editor%setAccessibleName(title//' text')
        call entry%editor%setReadOnly(.true.)
        call entry%editor%setFrameShape(QFrame_NoFrame)
        call entry%editor%setSizePolicy(QSizePolicy_Ignored,QSizePolicy_Expanding)
        call entry%editor%setHorizontalScrollBarPolicy(Qt_ScrollBarAlwaysOff)
        call entry%editor%setVerticalScrollBarPolicy(Qt_ScrollBarAlwaysOff)
        call entry%editor%setStyleSheet('QTextEdit {background:transparent;border:none;padding:0;min-height:0;'// &
            'selection-background-color:#416783;'//text_style//'}')
        entry%measurement=entry%editor%verticalScrollBar()
        call entry%measurement%rangeChanged%connect(on_card_range_changed,entry)
        call inside%addWidget(entry%editor)
        if(entry%collapsed)call entry%editor%hide()
        call entry%editor%setFixedHeight(measurement_height)
        selection_events=QtfEventFilter(entry%editor)
        call selection_events%focusOutEvent%connect(on_editor_focus_lost,entry)
        wheel_events=QtfEventFilter(entry%editor%viewport())
        call wheel_events%wheelEvent%connect(on_editor_wheel,entry)
        call cards%addWidget(row)
        if(entry%role=='You')then
            metadata_size=metadata%sizeHint()
            entry%preferred_width=max(metadata_size%width,question_width(entry))+2*bubble_padding
        end if
        call resize_card(entry)
        call render_entry(entry)
    end subroutine make_card

    integer function question_width(entry) result(width)
        type(timeline_entry), intent(in) :: entry
        type(QFontMetrics) :: metrics
        integer :: first,last

        metrics=QFontMetrics(entry%editor%font())
        width=0
        first=1
        do while(first<=len(entry%body))
            last=index(entry%body(first:),lf)
            if(last==0)last=len(entry%body)-first+2
            width=max(width,metrics%horizontalAdvance(entry%body(first:first+last-2)))
            first=first+last
        end do
        width=width+12
        call metrics%delete()
    end function question_width

    subroutine resize_card(entry)
        type(timeline_entry), intent(inout) :: entry
        integer :: available,width

        available=max(1,entry%scroll_viewport%width()-2*timeline_margin)
        select case(entry%role)
        case('You')
            width=min(entry%preferred_width,520,max(1,available*5/6))
        case('Assistant')
            width=min(620,available)
        case default
            width=min(540,max(1,available-14))
        end select
        if(entry%card_width==width)return
        entry%card_width=width
        call entry%card%setFixedWidth(width)
        entry%needs_measure=.true.
    end subroutine resize_card

    subroutine update_disclosure(entry)
        type(timeline_entry), intent(inout) :: entry
        character(:), allocatable :: label

        select case(entry%role)
        case('Plan'); label=' the plan'
        case('Tool'); label=' tool details'
        case default; label=' agent reasoning'
        end select
        call entry%disclosure%setText(trim(merge('Show','Hide',entry%collapsed)))
        call entry%disclosure%setAccessibleName(entry%disclosure%text()//label)
    end subroutine update_disclosure

    subroutine toggle_detail(data)
        class(*), intent(inout), optional :: data

        select type(entry=>data)
        type is(timeline_entry)
            entry%collapsed=.not.entry%collapsed
            call entry%editor%setVisible(.not.entry%collapsed)
            call update_disclosure(entry)
            entry%needs_measure=.true.
            call entry%measurement_timer%start(0)
        end select
    end subroutine toggle_detail

    subroutine render_entry(entry,refresh)
        type(timeline_entry), intent(inout) :: entry
        logical, intent(in), optional :: refresh
        logical :: keep_selection
        integer :: detail_start

        keep_selection=entry%editor%hasFocus()
        if(present(refresh))keep_selection=keep_selection.and..not.refresh
        if(keep_selection)then
            entry%needs_render=.true.
            return
        end if
        entry%needs_render=.false.
        if(entry%role=='Assistant')then
            call entry%editor%setHtml(assistant_html(entry%body))
        else if(entry%role=='Tool')then
            detail_start=index(entry%body,lf)
            if(detail_start>0)call entry%editor%setPlainText(entry%body(detail_start+1:))
        else
            call entry%editor%setPlainText(entry%body)
        end if
        entry%needs_measure=.true.
    end subroutine render_entry

    subroutine copy_card_text(data)
        class(*), intent(inout), optional :: data
        type(QClipboard) :: clipboard

        select type(entry=>data)
        type is(timeline_entry)
            clipboard=QGuiApplication_clipboard()
            call clipboard%setText(entry%body)
        end select
    end subroutine copy_card_text

    subroutine on_editor_focus_lost(reason,data)
        integer, intent(in) :: reason
        class(*), intent(inout), optional :: data

        select type(entry=>data)
        type is(timeline_entry)
            if(.not.entry%needs_render)return
            call render_entry(entry,refresh=.true.)
            call measure_card(entry)
        end select
    end subroutine on_editor_focus_lost

    subroutine on_editor_wheel(event,data)
        type(QWheelEvent), intent(in) :: event
        class(*), intent(inout), optional :: data

        select type(entry=>data)
        type is(timeline_entry)
            if(QCoreApplication_sendEvent(entry%scroll_viewport,event))call event%accept()
        end select
    end subroutine on_editor_wheel

    subroutine on_card_range_changed(minimum,maximum,data)
        integer, intent(in) :: minimum,maximum
        class(*), intent(inout), optional :: data

        select type(entry=>data)
        type is(timeline_entry)
            if(entry%measuring.or.maximum<=minimum)return
            entry%needs_measure=.true.
            call entry%measurement_timer%start(0)
        end select
    end subroutine on_card_range_changed

    function local_timestamp() result(stamp)
        character(timestamp_size) :: stamp
        integer :: clock(8)

        call date_and_time(values=clock)
        write(stamp,'(i4.4,2(a,i2.2),a,i2.2,2(a,i2.2))') clock(1),'-',clock(2),'-',clock(3),' ',clock(5),':',clock(6),':',clock(7)
    end function local_timestamp

    function serialized(transcript) result(text)
        class(assistant_transcript_t), intent(in) :: transcript
        character(:), allocatable :: text
        character(length_size) :: size_field
        character(agent_length_size) :: agent_field
        type(timeline_entry), pointer :: entry
        integer :: bytes,agent_bytes,first,last
        logical :: attributed

        text=''
        if(.not.associated(transcript%first))return
        bytes=len(storage_header)
        agent_bytes=0
        attributed=.false.
        entry=>transcript%first
        do while(associated(entry))
            bytes=bytes+record_size+len(entry%body)
            agent_bytes=agent_bytes+agent_length_size+len(entry%agent)
            if(len(entry%agent)>0)attributed=.true.
            entry=>entry%next
        end do
        if(attributed)bytes=bytes+agent_bytes
        text=repeat(' ',bytes)
        text(:len(storage_header))=merge(agent_storage_header,storage_header,attributed)
        first=len(storage_header)+1
        entry=>transcript%first
        do while(associated(entry))
            write(size_field,'(i9.9)')len(entry%body)
            last=first+record_size-1
            text(first:last)=role_code(entry%role)//entry%timestamp//size_field
            first=last+1
            if(attributed)then
                write(agent_field,'(i2.2)')len(entry%agent)
                last=first+agent_length_size+len(entry%agent)-1
                text(first:last)=agent_field//entry%agent
                first=last+1
            end if
            last=first+len(entry%body)-1
            text(first:last)=entry%body
            first=last+1
            entry=>entry%next
        end do
    end function serialized

    function plain_text(transcript) result(text)
        class(assistant_transcript_t), intent(in) :: transcript
        character(:), allocatable :: text,title
        type(timeline_entry), pointer :: entry

        text=''
        entry=>transcript%first
        do while(associated(entry))
            if(len(text)>0)text=text//gap
            title=timeline_title(entry)
            if(len_trim(entry%timestamp)>0)title=title//' · '//entry%timestamp
            if(entry%role=='Conversation'.and.len_trim(entry%timestamp)==0)then
                text=text//entry%body
            else
                text=text//title//lf//entry%body
            end if
            entry=>entry%next
        end do
    end function plain_text

    subroutine restore(transcript,text)
        class(assistant_transcript_t), intent(inout) :: transcript
        character(*), intent(in) :: text
        character(:), allocatable :: role,agent
        character(timestamp_size) :: stamp
        integer :: first,last,length,agent_length,status
        logical :: attributed

        call transcript%release()
        call build_host(transcript)
        attributed=index(text,agent_storage_header)==1
        if(attributed.or.index(text,storage_header)==1)then
            first=len(storage_header)+1
            do while(first<=len(text))
                last=first+record_size-1
                if(last>len(text))exit
                role=record_role(text(first:first))
                if(len(role)==0)exit
                stamp=text(first+1:first+timestamp_size)
                if(verify(text(last-length_size+1:last),digits)/=0)exit
                read(text(last-length_size+1:last),*,iostat=status)length
                if(status/=0)exit
                agent=''
                if(attributed)then
                    last=last+agent_length_size
                    if(last>len(text))exit
                    if(verify(text(last-agent_length_size+1:last),digits)/=0)exit
                    read(text(last-agent_length_size+1:last),*,iostat=status)agent_length
                    if(status/=0)exit
                    if(agent_length>max_agent_label.or.agent_length>len(text)-last)exit
                    agent=agent_label(text(last+1:last+agent_length))
                    last=last+agent_length
                end if
                if(length<1.or.length>len(text)-last)exit
                call append_entry(transcript,role,text(last+1:last+length),stamp,agent)
                first=last+length+1
            end do
            if(first<=len(text).and..not.associated(transcript%first))call restore_legacy(transcript,text)
        else
            call restore_legacy(transcript,text)
        end if
        call transcript%follow()
    end subroutine restore

    subroutine restore_legacy(transcript,text)
        class(assistant_transcript_t), intent(inout) :: transcript
        character(*), intent(in) :: text
        character(:), allocatable :: role,next_role,line
        integer :: i,head,body,next_break

        if(len(text)==0)return
        role='Conversation'
        head=1
        body=1
        do i=1,len(text)
            if(i/=1)then
                if(i<3)cycle
                if(text(i-2:i-1)/=gap)cycle
            end if
            next_break=index(text(i:),lf)
            if(next_break==0)cycle
            line=text(i:i+next_break-2)
            select case(line)
            case('You','Assistant')
                next_role=line
            case('Studio','ForCAD Studio')
                next_role='Studio'
            case('Studio activity','ForCAD Studio activity')
                next_role='Activity'
            case default
                cycle
            end select
            if(i>head)call append_entry(transcript,role,text(body:i-3),repeat(' ',timestamp_size))
            role=next_role
            head=i
            body=i+next_break
        end do
        if(body<=len(text))call append_entry(transcript,role,text(body:),repeat(' ',timestamp_size))
    end subroutine restore_legacy

    function role_code(role) result(code)
        character(*), intent(in) :: role
        character :: code

        select case(role)
        case('You'); code='Y'
        case('Assistant'); code='A'
        case('Studio'); code='S'
        case('Activity'); code='E'
        case('Reasoning'); code='R'
        case('Plan'); code='P'
        case('Tool'); code='T'
        case default; code='L'
        end select
    end function role_code

    function record_role(code) result(role)
        character, intent(in) :: code
        character(:), allocatable :: role

        select case(code)
        case('Y'); role='You'
        case('A'); role='Assistant'
        case('S'); role='Studio'
        case('E'); role='Activity'
        case('R'); role='Reasoning'
        case('P'); role='Plan'
        case('T'); role='Tool'
        case('L'); role='Conversation'
        case default; role=''
        end select
    end function record_role

    subroutine queue_layout(transcript)
        class(assistant_transcript_t), intent(inout) :: transcript

        transcript%layout_pending=.true.
        call transcript%tail_timer%start(0)
    end subroutine queue_layout

    logical function defer_update(transcript) result(deferred)
        class(assistant_transcript_t), intent(in) :: transcript

        deferred=.false.
        if(transcript%scrollbar%isNull())return
        deferred=transcript%scrollbar%isSliderDown()
    end function defer_update

    logical function at_end(transcript)
        class(assistant_transcript_t), intent(in) :: transcript

        at_end=.true.
        if(transcript%scrollbar%isNull())return
        at_end=transcript%scrollbar%value()>=transcript%scrollbar%maximum()
    end function at_end

    subroutine follow(transcript)
        class(assistant_transcript_t), intent(inout) :: transcript

        if(transcript%scrollbar%isNull())return
        transcript%follow_tail=.true.
        call queue_layout(transcript)
        call scroll_to_end(transcript)
    end subroutine follow

    subroutine finish_follow(transcript)
        class(assistant_transcript_t), intent(inout) :: transcript
        type(timeline_entry), pointer :: entry
        logical :: remeasure

        if(transcript%scrollbar%isNull())return
        transcript%layout_pending=.true.
        transcript%moving_scrollbar=.true.
        remeasure=.false.
        entry=>transcript%first
        do while(associated(entry))
            call resize_card(entry)
            if(entry%needs_measure.or.entry%measured_width/=entry%editor%width())then
                call measure_card(entry)
                remeasure=.true.
            end if
            entry=>entry%next
        end do
        transcript%moving_scrollbar=.false.
        if(.not.transcript%scrollbar%isSliderDown())then
            if(transcript%follow_tail)then
                call scroll_to_end(transcript)
            else
                transcript%moving_scrollbar=.true.
                call transcript%scrollbar%setValue(transcript%scroll_value)
                transcript%moving_scrollbar=.false.
            end if
        end if
        if(remeasure)then
            call transcript%tail_timer%start(0)
        else
            transcript%layout_pending=.false.
        end if
        call sync_jump(transcript)
    end subroutine finish_follow

    subroutine measure_card(entry)
        type(timeline_entry), intent(inout) :: entry
        integer :: desired_height

        entry%measuring=.true.
        if(.not.entry%collapsed)then
            call entry%editor%setFixedHeight(measurement_height)
            desired_height=entry%editor%height()+entry%measurement%maximum()
            call entry%editor%setFixedHeight(desired_height)
        end if
        entry%measured_width=entry%editor%width()
        entry%needs_measure=.false.
        entry%measuring=.false.
    end subroutine measure_card

    subroutine scroll_to_end(transcript)
        class(assistant_transcript_t), intent(inout) :: transcript

        transcript%moving_scrollbar=.true.
        if(.not.transcript%jump_button%isNull())call transcript%jump_button%hide()
        call transcript%scrollbar%setValue(transcript%scrollbar%maximum())
        transcript%scroll_value=transcript%scrollbar%value()
        transcript%moving_scrollbar=.false.
    end subroutine scroll_to_end

    subroutine sync_jump(transcript)
        class(assistant_transcript_t), intent(inout) :: transcript
        logical :: show_jump

        if(transcript%moving_scrollbar.or.transcript%jump_button%isNull())return
        show_jump=.not.transcript%at_end().and..not.transcript%follow_tail
        if(transcript%jump_button%isVisible().eqv.show_jump)return
        transcript%moving_scrollbar=.true.
        call transcript%jump_button%setVisible(show_jump)
        transcript%moving_scrollbar=.false.
    end subroutine sync_jump

    subroutine on_scroll_changed(value,data)
        integer, intent(in) :: value
        class(*), intent(inout), optional :: data

        select type(transcript=>data)
        type is(assistant_transcript_t)
            if(transcript%moving_scrollbar)return
            if(transcript%layout_pending.and..not.transcript%scrollbar%isSliderDown())return
            transcript%follow_tail=value>=transcript%scrollbar%maximum()
            transcript%scroll_value=value
            call sync_jump(transcript)
        end select
    end subroutine on_scroll_changed

    subroutine on_range_changed(minimum,maximum,data)
        integer, intent(in) :: minimum,maximum
        class(*), intent(inout), optional :: data

        select type(transcript=>data)
        type is(assistant_transcript_t)
            if(transcript%moving_scrollbar)return
            call queue_layout(transcript)
        end select
    end subroutine on_range_changed

    subroutine on_resize(event,data)
        type(QResizeEvent), intent(in) :: event
        class(*), intent(inout), optional :: data

        select type(transcript=>data)
        type is(assistant_transcript_t)
            if(transcript%moving_scrollbar)return
            call queue_layout(transcript)
        end select
    end subroutine on_resize

    function assistant_html(body) result(html)
        character(*), intent(in) :: body
        character(:), allocatable :: html,line,display,buffer
        integer :: first,i,cursor,indent,heading_end
        logical :: fenced

        allocate(character(8*len(body)+320) :: buffer)
        cursor=0
        call append_html(buffer,cursor,'<html><head><style>code {color:#a8e3e8;background-color:#142330;}'// &
            '</style></head><body><p style="margin:0">')
        first=1
        fenced=.false.
        do i=1,len(body)+1
            if(i<=len(body))then
                if(body(i:i)/=lf)cycle
            end if
            line=body(first:i-1)
            display=adjustl(line)
            if(first>1)call append_html(buffer,cursor,'<br/>')
            if(index(display,'```')==1.or.index(display,'~~~')==1)then
                if(.not.fenced)call append_html(buffer,cursor,'<code>')
                call append_html(buffer,cursor,escape_html(trim(display)))
                if(fenced)call append_html(buffer,cursor,'</code>')
                fenced=.not.fenced
            else if(fenced)then
                indent=verify(line,' ')-1
                call append_html(buffer,cursor,repeat('&nbsp;',max(0,indent))//escape_html(trim(display)))
            else
                heading_end=verify(display,'#')
                if(heading_end>1.and.heading_end<len(display))then
                    if(display(heading_end:heading_end)==' ')then
                        call append_html(buffer,cursor,'<b>'//escape_html(trim(display(heading_end+1:)))//'</b>')
                        first=i+1
                        cycle
                    end if
                end if
                if(index(display,'- ')==1.or.index(display,'* ')==1)display='• '//display(3:)
                call append_inline(buffer,cursor,trim(display))
            end if
            first=i+1
        end do
        if(fenced)call append_html(buffer,cursor,'</code>')
        call append_html(buffer,cursor,'</p></body></html>')
        html=buffer(:cursor)
    end function assistant_html

    subroutine append_inline(buffer,cursor,line)
        character(:), allocatable, intent(inout) :: buffer
        integer, intent(inout) :: cursor
        character(*), intent(in) :: line
        character(:), allocatable :: marker,opening,closing
        integer :: start,position,ending

        start=1
        do while(start<=len(line))
            position=scan(line(start:),'*`')+start-1
            if(position<start)exit
            call append_html(buffer,cursor,escape_html(line(start:position-1)))
            marker='`'
            opening='<code>'
            closing='</code>'
            if(line(position:position)=='*')then
                marker='**'
                opening='<b>'
                closing='</b>'
            end if
            start=position+len(marker)
            if(index(line(position:),marker)==1)then
                ending=index(line(start:),marker)+start-1
                if(ending>=start)then
                    call append_html(buffer,cursor,opening//escape_html(line(start:ending-1))//closing)
                    start=ending+len(marker)
                    cycle
                end if
            end if
            call append_html(buffer,cursor,escape_html(line(position:position)))
            start=position+1
        end do
        if(start<=len(line))call append_html(buffer,cursor,escape_html(line(start:)))
    end subroutine append_inline

    function escape_html(input) result(html)
        character(*), intent(in) :: input
        character(:), allocatable :: html,buffer
        integer :: i,cursor

        allocate(character(5*len(input)) :: buffer)
        cursor=0
        do i=1,len(input)
            select case(input(i:i))
            case('&'); call append_html(buffer,cursor,'&amp;')
            case('<'); call append_html(buffer,cursor,'&lt;')
            case('>'); call append_html(buffer,cursor,'&gt;')
            case default; call append_html(buffer,cursor,input(i:i))
            end select
        end do
        html=buffer(:cursor)
    end function escape_html

    subroutine append_html(buffer,cursor,fragment)
        character(:), allocatable, intent(inout) :: buffer
        integer, intent(inout) :: cursor
        character(*), intent(in) :: fragment
        integer :: last

        last=cursor+len(fragment)
        buffer(cursor+1:last)=fragment
        cursor=last
    end subroutine append_html

    function chat_body(value) result(text)
        character(*), intent(in) :: value
        character(:), allocatable :: text,input,line
        character(len(value)) :: buffer
        integer :: i,first,last,used,width,blank_lines
        logical :: in_code

        input=without_nul(value)
        first=1
        used=0
        blank_lines=0
        in_code=.false.
        do i=1,len(input)+1
            if(i<=len(input))then
                if(input(i:i)/=lf)cycle
            end if
            last=i-1
            if(last>=first)then
                if(input(last:last)==achar(13))last=last-1
            end if
            line=input(first:last)
            if(verify(line,' '//achar(9))==0)then
                blank_lines=blank_lines+1
            else
                if(used>0)then
                    width=1+min(blank_lines,1)
                    if(in_code)width=1+blank_lines
                    buffer(used+1:used+width)=repeat(lf,width)
                    used=used+width
                end if
                width=len(line)
                buffer(used+1:used+width)=line
                used=used+width
                blank_lines=0
                line=adjustl(line)
                if(index(line,'```')==1.or.index(line,'~~~')==1)in_code=.not.in_code
            end if
            first=i+1
        end do
        text=buffer(:used)
    end function chat_body

end module studio_assistant_transcript
