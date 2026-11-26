module studio_assistant_panel
    use qt, only: QWidget, QMainWindow, QDockWidget, QVBoxLayout, QHBoxLayout, QLabel, QComboBox, &
        QtfPlainTextEdit, QtfEventFilter, QResizeEvent, QPushButton, QToolButton, QStackedWidget, &
        QAction, QMenu, QScrollArea, QScrollBar, QDialog, QListWidget, QDialogButtonBox, QDateTime, QDateTime_fromString, &
        QIcon, QDockWidget_DockWidgetClosable, QDockWidget_DockWidgetMovable, QDockWidget_DockWidgetFloatable, &
        QRectF, QSize, QFontMetrics, Qt_LeftDockWidgetArea, Qt_RightDockWidgetArea, Qt_PlainText, &
        Qt_AlignHCenter, Qt_StrongFocus, Qt_ScrollBarAlwaysOff, QSizePolicy_Ignored, QSizePolicy_Fixed, QSizePolicy_Preferred, &
        QPlainTextEdit_WidgetWidth, Qt_ToolButtonIconOnly, QToolButton_InstantPopup, QDialogButtonBox_Close, &
        QDialogButtonBox_AcceptRole, QDialogButtonBox_DestructiveRole, &
        QDialog_Accepted, QMessageBox, QAbstractButton, QMessageBox_Icon_Question, QMessageBox_Discard, QMessageBox_Cancel, &
        QMessageBox_Save, QMessageBox_Yes, QMessageBox_warning
    use studio_icons, only: line_icon, stylesheet_icon, icon_close, icon_float, icon_more, icon_refresh, &
        icon_new_chat, icon_history, icon_send, icon_stop, icon_fit_selection, icon_chevron, icon_export, icon_box_selection, &
        icon_delete, icon_copy, icon_screenshot, icon_settings

    use studio_chat_history, only: chat_summary_t, delete_chat
    use studio_theme, only: studio_frame_edge, studio_divider, studio_heading, studio_style_primary, &
        studio_secondary_text, studio_selected_fill, studio_active_edge, studio_tooltip
    use studio_status_hints, only: set_status_hint
    implicit none(type, external)
    private

    public :: assistant_panel_t, build_assistant_panel, choose_saved_chat, confirm_chat_replacement

    integer, parameter :: composer_min_height=48, composer_max_height=116, send_button_height=32
    character(*), parameter :: saved_chats_label='Saved chats', save_chat_label='Save chat'

    type :: assistant_panel_t
        type(QDockWidget) :: dock
        type(QLabel) :: provider_label
        type(QListWidget) :: attachment_list
        type(QComboBox) :: model, agent_model, agent_effort, approval_mode
        type(QWidget) :: effort_row
        type(QtfPlainTextEdit) :: prompt
        type(QScrollArea) :: transcript
        type(QPushButton) :: jump_to_latest
        type(QScrollBar) :: transcript_scrollbar
        type(QStackedWidget) :: send_or_stop
        type(QPushButton) :: ask, stop, new_chat, refresh_models, detach
        type(QPushButton) :: studio_tools, clear_context, continue_scene, team_mode, team_settings
        type(QAction) :: save_button, history_button, copy_chat, export_chat, attach, attach_selection, attach_view
        type(QWidget) :: permission_box, attachment_bar, context_bar
        type(QLabel) :: permission_title, permission_summary, attached_context
        type(QPushButton) :: allow_once, deny_request
    contains
        procedure :: sync => sync_assistant_panel
    end type assistant_panel_t

    type :: history_picker
        type(QDialog) :: dialog
        type(QListWidget) :: list
        type(QPushButton) :: open_button, delete_button
        type(QLabel) :: empty_state
        type(chat_summary_t), allocatable :: items(:)
        character(:), allocatable :: active_id,file_name
        logical :: deleted_active=.false., deleted=.false.
    end type history_picker

contains

    subroutine build_assistant_panel(parent,overlay,panel,dock_host,content)
        class(QWidget), intent(in) :: parent
        type(QWidget), intent(out) :: overlay,content
        type(assistant_panel_t), target, intent(out) :: panel
        type(QMainWindow), intent(out) :: dock_host
        type(QVBoxLayout) :: layout, conversation_layout
        type(QWidget) :: conversation

        call create_assistant_dock(parent,overlay,dock_host,panel%dock)
        content=QWidget()
        call content%setObjectName('studioAssistantBody')
        call content%setMinimumWidth(310)
        layout=QVBoxLayout(content)
        call layout%setContentsMargins(10,8,10,8)
        call layout%setSpacing(8)
        call build_assistant_header(layout,content,panel)
        conversation=QWidget(content)
        call conversation%setObjectName('studioAssistantConversation')
        conversation_layout=QVBoxLayout(conversation)
        call conversation_layout%setContentsMargins(0,0,0,0)
        call conversation_layout%setSpacing(8)
        call build_transcript_area(conversation_layout,panel)
        call build_message_composer(conversation_layout,conversation,panel)
        call layout%addWidget(conversation,1)
        call build_continue_scene(layout,panel)
    end subroutine build_assistant_panel

    subroutine create_assistant_dock(parent,overlay,dock_host,dock)
        class(QWidget), intent(in) :: parent
        type(QWidget), intent(out) :: overlay
        type(QMainWindow), intent(out) :: dock_host
        type(QDockWidget), intent(out) :: dock
        type(QVBoxLayout) :: overlay_layout
        integer :: features,areas

        ! The nested host keeps native detach/redock while owning no canvas layout space.
        overlay=QWidget(parent)
        call overlay%setObjectName('studioAssistantOverlay')
        overlay_layout=QVBoxLayout(overlay)
        call overlay_layout%setContentsMargins(1,1,1,1)
        dock_host=QMainWindow()
        call dock_host%setObjectName('studioAssistantHost')
        call overlay_layout%addWidget(dock_host)
        call overlay%setStyleSheet( &
            'QWidget#studioAssistantOverlay { background:#16212a; border:1px solid '//studio_frame_edge//'; }'// &
            'QMainWindow#studioAssistantHost { background:transparent; }')
        dock=QDockWidget('Assistant',dock_host)
        call dock%setObjectName('studioAssistantDock')
        call dock%setAccessibleName('Assistant')
        features=ior(QDockWidget_DockWidgetClosable, &
            ior(QDockWidget_DockWidgetMovable,QDockWidget_DockWidgetFloatable))
        call dock%setFeatures(features)
        areas=ior(Qt_LeftDockWidgetArea,Qt_RightDockWidgetArea)
        call dock%setAllowedAreas(areas)
        call dock%setStyleSheet( &
            'QDockWidget#studioAssistantDock { background:transparent; color:#edf3f7; font-size:14px; font-weight:700;'// &
            stylesheet_icon(icon_close,16,16,'titlebar-close-icon')// &
            stylesheet_icon(icon_float,16,16,'titlebar-normal-icon')//' }'// &
            'QDockWidget#studioAssistantDock::title { text-align:left; padding:8px 12px;'// &
            'background:'//studio_heading//'; border-bottom:1px solid '//studio_divider//'; }'// &
            'QDockWidget#studioAssistantDock::float-button,QDockWidget#studioAssistantDock::close-button {'// &
            'subcontrol-origin:margin; subcontrol-position:center right; position:absolute;'// &
            'width:24px; height:24px; icon-size:16px; border:1px solid transparent; }'// &
            'QDockWidget#studioAssistantDock::float-button { right:32px; }'// &
            'QDockWidget#studioAssistantDock::close-button { right:4px; }'// &
            'QDockWidget#studioAssistantDock::float-button:hover { background:#30414e; border-color:#536b79; }'// &
            'QDockWidget#studioAssistantDock::close-button:hover { background:#593844; border-color:#805160; }'// &
            'QWidget#studioAssistantBody { background:#16212a; }'// &
            'QLabel[assistantOptionLabel="true"] { color:'//studio_secondary_text//'; font-size:12px; font-weight:600; }'// &
            'QWidget#studioAssistantConversation { background:#16212a; border:0; }'// &
            'QScrollArea#studioAssistantTranscript { background:#16212a; color:#edf3f7; border:0; }'// &
            'QPushButton#studioAssistantJumpToLatest { background:#22303b; border:1px solid '//studio_frame_edge//';'// &
            'padding:0; }'// &
            'QPushButton#studioAssistantJumpToLatest:hover { background:#30414e; }'// &
            'QWidget#studioAssistantComposer { background:#1b2833; border:1px solid '//studio_frame_edge//'; }'// &
            'QPlainTextEdit#studioAssistantPrompt { background:transparent; color:#edf3f7; padding:3px 4px; font-size:14px;'// &
            'border:1px solid transparent; placeholder-text-color:'//studio_secondary_text//'; }'// &
            'QPlainTextEdit#studioAssistantPrompt:focus { background:#16212a; border-color:'//studio_active_edge//'; }'// &
            'QWidget[assistantAttachment="true"] { background:#223440; border:1px solid '//studio_frame_edge//'; }'// &
            'QWidget[assistantAttachment="true"] QLabel { color:#d8e7ef; font-size:12px; }'// &
            'QListWidget#studioAssistantAttachments { color:#d8e7ef; font-size:12px; background:transparent; border:none; }'// &
            'QListWidget#studioAssistantAttachments::item { padding:3px 4px; }'// &
            'QListWidget#studioAssistantAttachments::item:selected { background:'//studio_selected_fill//'; color:#edf3f7; }'// &
            'QWidget[assistantAttachment="true"] QPushButton { padding:0; border:1px solid transparent; }'// &
            'QWidget[assistantAttachment="true"] QPushButton:hover { background:#354b5a; }'// &
            'QWidget[assistantAttachment="true"] QPushButton:focus { border-color:#8de5e1; }'// &
            'QToolButton[assistantAttachmentAction="true"] { background:#22303b; padding:0; border:1px solid '//studio_frame_edge//'; }'// &
            'QToolButton[assistantAttachmentAction="true"]:hover { background:#30414e; border-color:#536b79; }'// &
            'QToolButton[assistantAttachmentAction="true"]:focus { border-color:'//studio_active_edge//'; }'// &
            'QToolButton[assistantAttachmentAction="true"]:disabled { background:#1c2933; }'// &
            'QPushButton[assistantAction="true"] { font-size:12px; padding:2px 6px; }'// &
            'QPushButton#studioAssistantRefreshModels,QPushButton#studioAssistantNewChat {'// &
            'background:#22303b; border-color:'//studio_frame_edge//'; }'// &
            'QPushButton#studioAssistantRefreshModels:hover,QPushButton#studioAssistantNewChat:hover {'// &
            'background:#30414e; border-color:#536b79; }'// &
            'QPushButton#studioAssistantRefreshModels:focus,QPushButton#studioAssistantNewChat:focus {'// &
            'border-color:'//studio_active_edge//'; }'// &
            'QToolButton[assistantMenu="true"] { padding:2px 4px; background:#22303b; border:1px solid '//studio_frame_edge//'; }'// &
            'QToolButton[assistantMenu="true"]::menu-indicator { width:0; image:none; }'// &
            'QToolButton[assistantMenu="true"]:hover { background:#30414e; border-color:#536b79; }'// &
            'QToolButton[assistantMenu="true"]:focus,QPushButton#studioAssistantJumpToLatest:focus {'// &
            'border-color:'//studio_active_edge//'; }'// &
            'QPushButton#studioAssistantStudioTools:checked { background:'//studio_selected_fill//'; color:#d2f4f1; border-color:#557e85; }'// &
            'QPushButton#studioAssistantStudioTools:checked:hover { background:#30535e; }'// &
            'QPushButton#studioAssistantStudioTools:disabled { color:#9eb1be; border-color:#344653; }'// &
            'QPushButton#studioAssistantAsk,QPushButton#studioAssistantStop {'// &
            'padding:3px 7px; font-size:13px; font-weight:600; }'// &
            'QPushButton#studioAssistantAsk:enabled { background:'//studio_selected_fill//'; border:1px solid '//studio_active_edge//'; color:#edf3f7; }'// &
            'QPushButton#studioAssistantAsk:enabled:hover { background:#305966; }'// &
            'QPushButton#studioAssistantAsk:disabled { background:'//studio_heading//'; border:1px solid #344653; color:#90a6b6; }'// &
            'QPushButton#studioAssistantStop { background:#4c3039; color:#ffede9; border:1px solid #a66675; }'// &
            'QPushButton#studioAssistantStop:enabled:hover { background:#65404a; }'// &
            'QPushButton#studioAssistantAsk:focus,QPushButton#studioAssistantStop:focus { border-color:#a5eeeb; }')

    end subroutine create_assistant_dock

    subroutine build_assistant_header(layout,content,panel)
        type(QVBoxLayout), intent(inout) :: layout
        type(QWidget), intent(in) :: content
        type(assistant_panel_t), intent(inout) :: panel
        type(QHBoxLayout) :: row
        type(QToolButton) :: more
        type(QMenu) :: chats

        row=QHBoxLayout()
        call row%setSpacing(6)
        panel%provider_label=QLabel('Agent')
        call panel%provider_label%setObjectName('assistantProvider')
        call panel%provider_label%setProperty('assistantOptionLabel',.true.)
        call panel%provider_label%setMinimumWidth(38)
        panel%model=QComboBox()
        call panel%model%setObjectName('studioAssistantModel')
        call panel%model%setAccessibleName('Assistant agent')
        call panel%model%setPlaceholderText('Choose an agent')
        call panel%model%setSizePolicy(QSizePolicy_Ignored,QSizePolicy_Fixed)
        call panel%model%setMinimumWidth(80)
        call panel%model%setFixedHeight(30)
        call row%addWidget(panel%provider_label)
        call row%addWidget(panel%model,1)

        panel%refresh_models=QPushButton('')
        call panel%refresh_models%setObjectName('studioAssistantRefreshModels')
        call panel%refresh_models%setAccessibleName('Refresh available agents and models')
        call set_status_hint(panel%refresh_models,'Refresh available agents and models')
        call panel%refresh_models%setFixedSize(30,30)
        call panel%refresh_models%setProperty('assistantAction',.true.)
        call set_assistant_button_icon(panel%refresh_models,icon_refresh,16)
        call row%addWidget(panel%refresh_models)

        panel%new_chat=QPushButton('')
        call panel%new_chat%setObjectName('studioAssistantNewChat')
        call panel%new_chat%setAccessibleName('New chat')
        call set_status_hint(panel%new_chat,'New chat')
        call panel%new_chat%setProperty('assistantAction',.true.)
        call panel%new_chat%setFixedSize(30,30)
        call set_assistant_button_icon(panel%new_chat,icon_new_chat,16)
        call row%addWidget(panel%new_chat)

        more=assistant_menu_button(content,'studioAssistantMore','Chat options',icon_more)
        chats=QMenu('Chat options',more)
        panel%save_button=assistant_action(chats,'studioAssistantSave',save_chat_label,3)
        panel%history_button=assistant_action(chats,'studioAssistantHistory',saved_chats_label//'…',icon_history)
        call chats%addSeparator()
        panel%copy_chat=assistant_action(chats,'studioAssistantCopyChat','Copy chat',icon_copy)
        panel%export_chat=assistant_action(chats,'studioAssistantExportChat','Export chat…',icon_export)
        call more%setMenu(chats)
        call row%addWidget(more)
        call layout%addLayout(row)
        call build_model_options_row(layout,panel)

        call build_studio_access_row(layout,panel)
    end subroutine build_assistant_header

    subroutine build_model_options_row(layout,panel)
        type(QVBoxLayout), intent(inout) :: layout
        type(assistant_panel_t), intent(inout) :: panel
        type(QHBoxLayout) :: row,effort_options
        type(QLabel) :: model_label,effort_label

        row=QHBoxLayout()
        call row%setSpacing(6)
        model_label=QLabel('Model')
        call model_label%setProperty('assistantOptionLabel',.true.)
        call model_label%setMinimumWidth(38)
        call row%addWidget(model_label)

        panel%agent_model=QComboBox()
        call panel%agent_model%setObjectName('studioAssistantAgentModel')
        call panel%agent_model%setAccessibleName('Assistant model')
        call panel%agent_model%setSizePolicy(QSizePolicy_Ignored,QSizePolicy_Fixed)
        call panel%agent_model%setFixedHeight(30)
        call panel%agent_model%setMinimumWidth(80)
        call panel%agent_model%setPlaceholderText('Choose a model')
        call row%addWidget(panel%agent_model,1)
        call layout%addLayout(row)

        panel%effort_row=QWidget()
        row=QHBoxLayout(panel%effort_row)
        call row%setContentsMargins(0,0,0,0)
        call row%setSpacing(6)
        effort_label=QLabel('Effort')
        call effort_label%setProperty('assistantOptionLabel',.true.)
        call effort_label%setMinimumWidth(38)
        call row%addWidget(effort_label)
        panel%agent_effort=QComboBox()
        call panel%agent_effort%setObjectName('studioAssistantAgentEffort')
        call panel%agent_effort%setAccessibleName('Model reasoning effort')
        call panel%agent_effort%setSizePolicy(QSizePolicy_Ignored,QSizePolicy_Fixed)
        call panel%agent_effort%setFixedHeight(30)
        call panel%agent_effort%setMinimumWidth(78)
        call panel%agent_effort%addItem('Default')
        call row%addWidget(panel%agent_effort,1)
        call panel%effort_row%hide()

        effort_options=QHBoxLayout()
        call effort_options%setSpacing(6)
        call effort_options%addWidget(panel%effort_row,1)
        call effort_options%addStretch()
        panel%team_mode=QPushButton('Team')
        call panel%team_mode%setObjectName('studioAssistantTeamMode')
        call panel%team_mode%setAccessibleName('Parallel agent team')
        call set_status_hint(panel%team_mode,'Run agents in parallel while the lead applies approved ForCAD Studio edits in order.')
        call panel%team_mode%setProperty('assistantAction',.true.)
        call panel%team_mode%setCheckable(.true.)
        call panel%team_mode%setFixedSize(58,30)
        call effort_options%addWidget(panel%team_mode)
        panel%team_settings=QPushButton('')
        call panel%team_settings%setObjectName('studioAssistantTeamSettings')
        call panel%team_settings%setAccessibleName('Team settings')
        call set_status_hint(panel%team_settings,studio_tooltip('Team settings','Choose an agent, model and effort for each team role.'))
        call panel%team_settings%setProperty('assistantAction',.true.)
        call panel%team_settings%setFixedSize(30,30)
        call set_assistant_button_icon(panel%team_settings,icon_settings,16)
        call panel%team_settings%hide()
        call effort_options%addWidget(panel%team_settings)
        call layout%addLayout(effort_options)
    end subroutine build_model_options_row

    subroutine build_studio_access_row(layout,panel)
        type(QVBoxLayout), intent(inout) :: layout
        type(assistant_panel_t), intent(inout) :: panel
        type(QHBoxLayout) :: row

        row=QHBoxLayout()
        call row%setSpacing(6)
        panel%studio_tools=QPushButton('Tools: Off')
        call panel%studio_tools%setObjectName('studioAssistantStudioTools')
        call panel%studio_tools%setAccessibleName('ForCAD Studio tool access')
        call panel%studio_tools%setProperty('assistantAction',.true.)
        call panel%studio_tools%setCheckable(.true.)
        call panel%studio_tools%setFixedHeight(30)
        call set_assistant_button_icon(panel%studio_tools,8,15)
        call row%addWidget(panel%studio_tools)

        panel%approval_mode=QComboBox()
        call panel%approval_mode%setObjectName('studioAssistantApprovalMode')
        call panel%approval_mode%setAccessibleName('ForCAD Studio tool approval mode')
        call set_status_hint(panel%approval_mode,'Choose whether the Assistant asks before changing the project or writing an export.')
        call panel%approval_mode%setFixedHeight(30)
        call panel%approval_mode%addItem('Ask for approval')
        call panel%approval_mode%addItem('Approve for me')
        call row%addWidget(panel%approval_mode,1)
        call layout%addLayout(row)
    end subroutine build_studio_access_row

    subroutine sync_assistant_panel(panel,active,show_effort)
        class(assistant_panel_t), intent(inout) :: panel
        logical, intent(in) :: active,show_effort

        call panel%send_or_stop%setCurrentIndex(merge(1,0,active))
        call panel%effort_row%setVisible(show_effort)
    end subroutine sync_assistant_panel

    subroutine build_transcript_area(layout,panel)
        type(QVBoxLayout), intent(inout) :: layout
        type(assistant_panel_t), intent(inout) :: panel

        panel%transcript=QScrollArea()
        call panel%transcript%setObjectName('studioAssistantTranscript')
        call panel%transcript%setAccessibleName('Assistant chat')
        call panel%transcript%setWidgetResizable(.true.)
        panel%transcript_scrollbar=panel%transcript%verticalScrollBar()
        call layout%addWidget(panel%transcript,1)
    end subroutine build_transcript_area

    subroutine build_message_composer(layout,content,panel)
        type(QVBoxLayout), intent(inout) :: layout
        type(QWidget), intent(in) :: content
        type(assistant_panel_t), target, intent(inout) :: panel
        type(QWidget) :: composer
        type(QVBoxLayout) :: input
        type(QHBoxLayout) :: footer
        type(QtfEventFilter) :: resize_events

        call build_permission_box(layout,content,panel)
        composer=QWidget(content)
        call composer%setObjectName('studioAssistantComposer')
        input=QVBoxLayout(composer)
        call input%setContentsMargins(8,8,8,8)
        call input%setSpacing(5)
        call build_attachment_list(input,panel)
        call build_attachment_chip(input,panel%context_bar,panel%attached_context,panel%clear_context, &
            'studioAssistantContextBar','studioAssistantAttachedContext','studioAssistantClearContext','Remove selection attachment')

        panel%prompt=QtfPlainTextEdit()
        call panel%prompt%setObjectName('studioAssistantPrompt')
        call panel%prompt%setAccessibleName('Message the Assistant')
        call panel%prompt%setPlaceholderText('Message the Assistant…')
        call set_status_hint(panel%prompt,'Enter to send · Shift+Enter for a new line')
        call panel%prompt%setLineWrapMode(QPlainTextEdit_WidgetWidth)
        call panel%prompt%setFixedHeight(composer_min_height)
        call panel%prompt%textChanged%connect(resize_composer,panel)
        resize_events=QtfEventFilter(panel%prompt)
        call resize_events%resizeEvent%connect(on_composer_resize,panel)
        call input%addWidget(panel%prompt)

        footer=QHBoxLayout()
        call footer%setSpacing(5)
        call add_attachment_action(footer,composer,panel%attach,'studioAssistantAttach','Attach files…', &
            'Choose one or more files to attach.',1)
        call add_attachment_action(footer,composer,panel%attach_selection,'studioAssistantAttachSelection','Attach selection', &
            'Include the current geometry or feature selection.',icon_box_selection)
        call add_attachment_action(footer,composer,panel%attach_view,'studioAssistantAttachView','Attach screenshot', &
            'Capture the current 3D view for this message.',icon_screenshot)

        panel%jump_to_latest=QPushButton('',composer)
        call panel%jump_to_latest%setObjectName('studioAssistantJumpToLatest')
        call panel%jump_to_latest%setAccessibleName('Jump to the newest message')
        call set_status_hint(panel%jump_to_latest,'Jump to the newest message')
        call panel%jump_to_latest%setFixedSize(send_button_height,send_button_height)
        call set_assistant_button_icon(panel%jump_to_latest,icon_chevron,16)
        call footer%addWidget(panel%jump_to_latest)
        call panel%jump_to_latest%hide()
        call build_send_controls(footer,panel)
        call input%addLayout(footer)
        call layout%addWidget(composer)
    end subroutine build_message_composer

    subroutine add_attachment_action(footer,composer,action,name,title,help,icon_kind)
        type(QHBoxLayout), intent(inout) :: footer
        type(QWidget), intent(in) :: composer
        type(QAction), intent(out) :: action
        character(*), intent(in) :: name,title,help
        integer, intent(in) :: icon_kind
        type(QToolButton) :: button
        type(QIcon) :: icon

        button=QToolButton(composer)
        call button%setObjectName(name//'Button')
        action=QAction(title,button)
        call action%setObjectName(name)
        call set_status_hint(action,studio_tooltip(title,help))
        icon=line_icon(icon_kind)
        call action%setIcon(icon)
        call icon%delete()
        call button%setDefaultAction(action)
        call button%setAccessibleName(title)
        call button%setAccessibleDescription(help)
        call set_status_hint(button,studio_tooltip(title,help))
        call button%setToolButtonStyle(Qt_ToolButtonIconOnly)
        call button%setFocusPolicy(Qt_StrongFocus)
        call button%setProperty('assistantAttachmentAction',.true.)
        call button%setIconSize(18,18)
        call button%setSizePolicy(QSizePolicy_Ignored,QSizePolicy_Fixed)
        call button%setMinimumWidth(34)
        call button%setFixedHeight(send_button_height)
        call footer%addWidget(button,1)
    end subroutine add_attachment_action

    subroutine build_send_controls(footer,panel)
        type(QHBoxLayout), intent(inout) :: footer
        type(assistant_panel_t), intent(inout) :: panel

        panel%ask=QPushButton('Send')
        call panel%ask%setObjectName('studioAssistantAsk')
        call panel%ask%setAccessibleName('Send message')
        call set_status_hint(panel%ask,'Send message · Enter')
        call panel%ask%setFixedSize(76,send_button_height)
        call panel%ask%setProperty('assistantAction',.true.)
        call set_assistant_button_icon(panel%ask,icon_send,15)
        panel%stop=QPushButton('Stop')
        call panel%stop%setObjectName('studioAssistantStop')
        call panel%stop%setAccessibleName('Stop the current assistant request')
        call set_status_hint(panel%stop,'Stop the current request')
        call panel%stop%setFixedSize(76,send_button_height)
        call panel%stop%setProperty('assistantAction',.true.)
        call set_assistant_button_icon(panel%stop,icon_stop,15)
        panel%send_or_stop=QStackedWidget()
        call panel%send_or_stop%setFixedSize(76,send_button_height)
        call panel%send_or_stop%addWidget(panel%ask)
        call panel%send_or_stop%addWidget(panel%stop)
        call panel%send_or_stop%setCurrentIndex(0)
        call footer%addWidget(panel%send_or_stop)
    end subroutine build_send_controls

    subroutine build_continue_scene(layout,panel)
        type(QVBoxLayout), intent(inout) :: layout
        type(assistant_panel_t), intent(inout) :: panel

        panel%continue_scene=QPushButton('Continue from current scene')
        call panel%continue_scene%setObjectName('studioAssistantContinueScene')
        call set_status_hint(panel%continue_scene,'Draft a follow-up using the geometry already committed')
        call panel%continue_scene%setProperty('assistantAction',.true.)
        call set_assistant_button_icon(panel%continue_scene,icon_fit_selection,15)
        call layout%addWidget(panel%continue_scene,alignment=Qt_AlignHCenter)
        call panel%continue_scene%hide()
    end subroutine build_continue_scene

    subroutine build_permission_box(layout,content,panel)
        type(QVBoxLayout), intent(inout) :: layout
        type(QWidget), intent(in) :: content
        type(assistant_panel_t), intent(inout) :: panel
        type(QVBoxLayout) :: permission_layout
        type(QHBoxLayout) :: permission_buttons

        panel%permission_box=QWidget(content)
        call panel%permission_box%setObjectName('studioAssistantPermission')
        call panel%permission_box%setAccessibleName('Pending tool approval')
        call panel%permission_box%setStyleSheet( &
            'QWidget#studioAssistantPermission { background:#2d2925; border:1px solid #8d6f4e; }'// &
            'QLabel#studioAssistantPermissionTitle { color:#f1c997; font-weight:600; }')
        permission_layout=QVBoxLayout(panel%permission_box)
        call permission_layout%setContentsMargins(8,6,8,6)
        call permission_layout%setSpacing(5)
        panel%permission_title=QLabel('Tool approval')
        call panel%permission_title%setObjectName('studioAssistantPermissionTitle')
        call permission_layout%addWidget(panel%permission_title)
        panel%permission_summary=QLabel('')
        call panel%permission_summary%setObjectName('studioAssistantPermissionSummary')
        call panel%permission_summary%setTextFormat(Qt_PlainText)
        call panel%permission_summary%setWordWrap(.true.)
        call permission_layout%addWidget(panel%permission_summary)
        permission_buttons=QHBoxLayout()
        call permission_buttons%addStretch()
        panel%deny_request=QPushButton('Deny')
        call panel%deny_request%setObjectName('studioAssistantPermissionDeny')
        call permission_buttons%addWidget(panel%deny_request)
        panel%allow_once=QPushButton('Allow once')
        call panel%allow_once%setObjectName('studioAssistantPermissionAllow')
        call permission_buttons%addWidget(panel%allow_once)
        call permission_layout%addLayout(permission_buttons)
        call layout%addWidget(panel%permission_box)
        call panel%permission_box%hide()
    end subroutine build_permission_box

    subroutine build_attachment_list(layout,panel)
        type(QVBoxLayout),intent(inout) :: layout
        type(assistant_panel_t),intent(inout) :: panel
        type(QHBoxLayout) :: row

        panel%attachment_bar=QWidget()
        call panel%attachment_bar%setObjectName('studioAssistantAttachmentBar')
        call panel%attachment_bar%setProperty('assistantAttachment',.true.)
        row=QHBoxLayout(panel%attachment_bar)
        call row%setContentsMargins(3,2,3,2)
        call row%setSpacing(5)

        panel%attachment_list=QListWidget(panel%attachment_bar)
        call panel%attachment_list%setObjectName('studioAssistantAttachments')
        call panel%attachment_list%setAccessibleName('Files attached to this message')
        call set_status_hint(panel%attachment_list,'Select an attachment to remove it.')
        call panel%attachment_list%setSizePolicy(QSizePolicy_Ignored,QSizePolicy_Fixed)
        call panel%attachment_list%setHorizontalScrollBarPolicy(Qt_ScrollBarAlwaysOff)
        call row%addWidget(panel%attachment_list,1)

        panel%detach=QPushButton('')
        call panel%detach%setObjectName('studioAssistantDetach')
        call panel%detach%setAccessibleName('Remove selected attachment')
        call set_status_hint(panel%detach,'Remove the selected attachment.')
        call panel%detach%setFixedSize(24,22)
        call set_assistant_button_icon(panel%detach,icon_close,14)
        call row%addWidget(panel%detach)
        call layout%addWidget(panel%attachment_bar)
        call panel%attachment_bar%hide()
    end subroutine build_attachment_list

    subroutine build_attachment_chip(layout,bar,label,remove,bar_name,label_name,button_name,button_tip)
        type(QVBoxLayout), intent(inout) :: layout
        type(QWidget), intent(out) :: bar
        type(QLabel), intent(out) :: label
        type(QPushButton), intent(out) :: remove
        character(*), intent(in) :: bar_name,label_name,button_name,button_tip
        type(QHBoxLayout) :: row

        bar=QWidget()
        call bar%setObjectName(bar_name)
        call bar%setProperty('assistantAttachment',.true.)
        row=QHBoxLayout(bar)
        call row%setContentsMargins(7,2,3,2)
        call row%setSpacing(5)
        label=QLabel('')
        call label%setObjectName(label_name)
        call label%setTextFormat(Qt_PlainText)
        call label%setSizePolicy(QSizePolicy_Ignored,QSizePolicy_Preferred)
        call row%addWidget(label,1)
        remove=QPushButton('')
        call remove%setObjectName(button_name)
        call remove%setAccessibleName(button_tip)
        call set_status_hint(remove,button_tip)
        call remove%setFixedSize(24,22)
        call set_assistant_button_icon(remove,icon_close,14)
        call row%addWidget(remove)
        call layout%addWidget(bar)
        call bar%hide()
    end subroutine build_attachment_chip

    subroutine resize_composer(data)
        class(*), intent(inout), optional :: data
        type(QFontMetrics) :: metrics
        type(QRectF) :: bounds
        integer :: block,line_height,height

        select type(panel=>data)
        type is(assistant_panel_t)
            metrics=QFontMetrics(panel%prompt%font())
            line_height=metrics%lineSpacing()
            call metrics%delete()
            height=18
            do block=0,min(4,panel%prompt%blockCount()-1)
                bounds=panel%prompt%blockBoundingGeometry(block)
                height=height+max(line_height,ceiling(bounds%height))
                if(height>=composer_max_height)exit
            end do
            height=max(composer_min_height,min(height,composer_max_height))
            if(panel%prompt%height()/=height)call panel%prompt%setFixedHeight(height)
        end select
    end subroutine resize_composer

    subroutine on_composer_resize(event,data)
        type(QResizeEvent), intent(in) :: event
        class(*), intent(inout), optional :: data
        type(QSize) :: previous,current

        previous=event%oldSize()
        current=event%size()
        if(previous%width/=current%width)call resize_composer(data)
    end subroutine on_composer_resize

    function assistant_menu_button(parent,name,tooltip,icon_kind) result(button)
        type(QWidget), intent(in) :: parent
        character(*), intent(in) :: name,tooltip
        integer, intent(in) :: icon_kind
        type(QToolButton) :: button

        button=QToolButton(parent)
        call button%setObjectName(name)
        call button%setAccessibleName(tooltip)
        call set_status_hint(button,tooltip)
        call button%setFixedSize(30,30)
        call button%setPopupMode(QToolButton_InstantPopup)
        call button%setProperty('assistantMenu',.true.)
        call set_assistant_button_icon(button,icon_kind,16)
    end function assistant_menu_button

    function assistant_action(menu,name,title,icon_kind) result(action)
        type(QMenu), intent(inout) :: menu
        character(*), intent(in) :: name,title
        integer, intent(in) :: icon_kind
        type(QAction) :: action
        type(QIcon) :: icon

        action=QAction(title,menu)
        call action%setObjectName(name)
        call set_status_hint(action,title)
        icon=line_icon(icon_kind)
        call action%setIcon(icon)
        call icon%delete()
        call menu%addAction(action)
    end function assistant_action

    integer function confirm_chat_replacement(parent,dirty,unsent) result(choice)
        class(QWidget), intent(in) :: parent
        logical, intent(in) :: dirty,unsent
        type(QMessageBox) :: dialog
        type(QAbstractButton) :: button
        integer :: buttons

        dialog=QMessageBox(parent)
        call dialog%setObjectName('studioAssistantUnsavedDialog')
        call dialog%setWindowTitle('Unsaved changes')
        call dialog%setIcon(QMessageBox_Icon_Question)
        buttons=ior(QMessageBox_Discard,QMessageBox_Cancel)
        if (dirty) then
            buttons=ior(buttons,QMessageBox_Save)
            call dialog%setText('Save this chat?')
            call dialog%setInformativeText('Reopen it later from '//saved_chats_label//'.')
            if (unsent) call dialog%setInformativeText('Continuing clears the unsent message.')
        else
            call dialog%setText('Discard the unsent message?')
        end if
        call dialog%setStandardButtons(buttons)
        call dialog%setDefaultButton(QMessageBox_Cancel)
        button=dialog%button(QMessageBox_Discard)
        call button%setObjectName('studioAssistantDiscard')
        button=dialog%button(QMessageBox_Cancel)
        call button%setObjectName('studioAssistantCancelReplace')
        if (dirty) then
            button=dialog%button(QMessageBox_Save)
            call button%setText(save_chat_label)
            call set_assistant_button_icon(button,3,16)
            call button%setObjectName('studioAssistantSaveBeforeReplace')
        end if
        choice=dialog%exec()
        call dialog%delete()
    end function confirm_chat_replacement

    function choose_saved_chat(parent,items,file_name,active_id,deleted_active,deleted) result(id)
        class(QWidget), intent(in) :: parent
        character(len=:), allocatable :: id,caption,row_text
        type(chat_summary_t), intent(in) :: items(:)
        character(*), intent(in) :: file_name,active_id
        logical, intent(out) :: deleted_active,deleted
        type(history_picker), target :: picker
        type(QVBoxLayout) :: layout
        type(QDialogButtonBox) :: buttons
        type(QHBoxLayout) :: heading
        type(QLabel) :: title,history_icon
        type(QDateTime) :: updated
        integer :: i,row

        picker%items=items
        picker%active_id=active_id
        picker%file_name=file_name
        picker%dialog=QDialog(parent)
        call picker%dialog%setObjectName('studioAssistantHistoryDialog')
        call picker%dialog%setWindowTitle(saved_chats_label)
        call picker%dialog%setMinimumSize(390,270)
        call picker%dialog%resize(490,360)
        call picker%dialog%setStyleSheet( &
            'QListWidget#studioAssistantHistoryList {'// &
            'background:#151f29; color:#edf3f7; border:1px solid '//studio_frame_edge//'; padding:5px; }'// &
            'QListWidget#studioAssistantHistoryList:focus { border-color:'//studio_active_edge//'; }'// &
            'QListWidget#studioAssistantHistoryList::item { padding:7px 9px; }'// &
            'QListWidget#studioAssistantHistoryList::item:hover:!selected { background:#283b49; }'// &
            'QListWidget#studioAssistantHistoryList::item:selected { background:'//studio_selected_fill//'; color:#f1f6ff; }')
        layout=QVBoxLayout(picker%dialog)
        call layout%setContentsMargins(18,16,18,16)
        call layout%setSpacing(10)

        heading=QHBoxLayout()
        call heading%setSpacing(8)
        history_icon=QLabel('',picker%dialog)
        call set_assistant_label_icon(history_icon,icon_history,21)
        title=QLabel(saved_chats_label,picker%dialog)
        call title%setObjectName('studioAssistantHistoryTitle')
        call title%setProperty('dialogRole','title')
        call heading%addWidget(history_icon)
        call heading%addWidget(title,1)
        call layout%addLayout(heading)

        picker%list=QListWidget(picker%dialog)
        call picker%list%setObjectName('studioAssistantHistoryList')
        call picker%list%setAccessibleName(saved_chats_label)
        call picker%list%setWordWrap(.true.)
        call picker%list%setSpacing(3)
        do i=1,size(items)
            caption=items(i)%updated
            updated=QDateTime_fromString(caption,'yyyy-MM-dd HH:mm:ss')
            if (updated%isValid()) caption=updated%toString('d MMM yyyy · HH:mm')
            row_text=items(i)%title//new_line('a')//caption
            call picker%list%addItem(row_text)
        end do
        call layout%addWidget(picker%list,1)
        picker%empty_state=QLabel('No saved chats.',picker%dialog)
        call picker%empty_state%setObjectName('studioAssistantHistoryEmpty')
        call picker%empty_state%setAlignment(Qt_AlignHCenter)
        call layout%addWidget(picker%empty_state,1)
        call picker%empty_state%hide()

        buttons=QDialogButtonBox(QDialogButtonBox_Close,picker%dialog)
        picker%delete_button=buttons%addButton('Delete chat…',QDialogButtonBox_DestructiveRole)
        call picker%delete_button%setObjectName('studioAssistantHistoryDelete')
        call set_assistant_button_icon(picker%delete_button,icon_delete,16)
        call picker%delete_button%clicked%connect(delete_current_chat,picker)
        picker%open_button=buttons%addButton('Open chat',QDialogButtonBox_AcceptRole)
        call picker%open_button%setObjectName('studioAssistantHistoryOpen')
        call studio_style_primary(picker%open_button)
        call picker%open_button%setDefault(.true.)
        call set_assistant_button_icon(picker%open_button,2,16)
        call layout%addWidget(buttons)
        call buttons%accepted%connect(open_current_chat,picker)
        call buttons%rejected%connect(cancel_chat_picker,picker)
        call picker%list%itemDoubleClicked%connect(open_chat_row,picker)
        call picker%list%currentRowChanged%connect(update_history_actions,picker)
        call picker%list%setCurrentRow(0)
        call picker%list%setFocus()

        id=''
        if (picker%dialog%exec()==QDialog_Accepted) then
            row=picker%list%currentRow()
            if (row>=0.and.row<size(picker%items)) id=picker%items(row+1)%id
        end if
        deleted_active=picker%deleted_active
        deleted=picker%deleted
        call picker%dialog%delete()
    end function choose_saved_chat

    subroutine update_history_actions(row,data)
        integer, intent(in) :: row
        class(*), intent(inout), optional :: data
        logical :: selected

        select type(picker=>data)
        type is(history_picker)
            selected=row>=0.and.row<size(picker%items)
            call picker%open_button%setEnabled(selected)
            call picker%delete_button%setEnabled(selected)
        end select
    end subroutine update_history_actions

    subroutine delete_current_chat(data)
        class(*), intent(inout), optional :: data
        character(:), allocatable :: id,error
        type(QMessageBox) :: confirmation
        type(QAbstractButton) :: button
        integer :: row,result
        logical :: active

        select type(picker=>data)
        type is(history_picker)
            row=picker%list%currentRow()
            if (row<0.or.row>=size(picker%items)) return
            id=picker%items(row+1)%id
            active=id==picker%active_id
            confirmation=QMessageBox(picker%dialog)
            call confirmation%setObjectName('studioAssistantDeleteChatDialog')
            call confirmation%setWindowTitle('Delete saved chat')
            call confirmation%setText('Delete the selected chat from this project’s saved chats?')
            if (active) call confirmation%setInformativeText('The open chat stays in the Assistant as an unsaved chat.')
            call confirmation%setIcon(QMessageBox_Icon_Question)
            call confirmation%setStandardButtons(ior(QMessageBox_Yes,QMessageBox_Cancel))
            call confirmation%setDefaultButton(QMessageBox_Cancel)
            button=confirmation%button(QMessageBox_Yes)
            call button%setText('Delete chat')
            call set_assistant_button_icon(button,icon_delete,16)
            call button%setObjectName('studioAssistantConfirmDeleteChat')
            result=confirmation%exec()
            call confirmation%delete()
            if (result/=QMessageBox_Yes) return

            call delete_chat(id,error,picker%file_name)
            if (len(error)>0) then
                result=QMessageBox_warning('Could not delete chat',error,parent=picker%dialog)
                return
            end if
            picker%deleted=.true.
            picker%deleted_active=picker%deleted_active.or.active
            picker%items=[picker%items(:row),picker%items(row+2:)]
            call picker%list%takeItem(row)
            if (size(picker%items)==0) then
                call picker%list%hide()
                call picker%empty_state%show()
                call update_history_actions(-1,picker)
            else
                call picker%list%setCurrentRow(min(row,size(picker%items)-1))
                call update_history_actions(picker%list%currentRow(),picker)
                call picker%list%setFocus()
            end if
        end select
    end subroutine delete_current_chat

    subroutine open_chat_row(row,data)
        integer, intent(in) :: row
        class(*), intent(inout), optional :: data
        select type(picker=>data)
        type is(history_picker)
            if (row>=0.and.row<picker%list%count()) then
                call picker%list%setCurrentRow(row)
                call picker%dialog%accept()
            end if
        end select
    end subroutine open_chat_row

    subroutine open_current_chat(data)
        class(*), intent(inout), optional :: data
        select type(picker=>data)
        type is(history_picker)
            call open_chat_row(picker%list%currentRow(),picker)
        end select
    end subroutine open_current_chat

    subroutine cancel_chat_picker(data)
        class(*), intent(inout), optional :: data
        select type(picker=>data)
        type is(history_picker)
            call picker%dialog%reject()
        end select
    end subroutine cancel_chat_picker

    subroutine set_assistant_button_icon(button,kind,extent)
        class(QAbstractButton), intent(inout) :: button
        integer, intent(in) :: kind,extent
        type(QIcon) :: icon

        icon=line_icon(kind)
        call button%setIcon(icon)
        call button%setIconSize(extent,extent)
        call icon%delete()
    end subroutine set_assistant_button_icon

    subroutine set_assistant_label_icon(label,kind,extent)
        type(QLabel), intent(inout) :: label
        integer, intent(in) :: kind,extent
        type(QIcon) :: icon

        icon=line_icon(kind,label)
        call label%setScaledContents(.true.)
        call label%setFixedSize(extent,extent)
        call icon%delete()
    end subroutine set_assistant_label_icon
end module studio_assistant_panel
