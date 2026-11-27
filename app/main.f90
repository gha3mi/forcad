program forcad_studio
    use studio_presentation_text, only: studio_product_name
    use qt, only: QApplication,QApplication_setOrganizationName,QApplication_setApplicationName,QGuiApplication_primaryScreen,QScreen
    use, intrinsic :: iso_fortran_env, only: error_unit
    use studio_display_refinement, only: finish_display_refinement
    use studio_state, only: studio_t,initialize_state
    use studio_ui, only: build_ui,on_studio_close,mark_project_saved,show_home,open_project_path, &
        studio_min_width,studio_min_height
    use studio_history_storage, only: finish_history_storage
    use studio_cli, only: launch_options,read_launch_options
    use studio_process_mode, only: dispatch_studio_process_mode
    use studio_capture, only: schedule_preview
    use studio_preview, only: apply_preview
    use studio_mcp_actions, only: execute_studio_mcp
    use studio_mcp_bridge, only: bind_studio_mcp,finish_studio_mcp,start_studio_mcp
    use studio_assistant, only: finish_studio_assistant,begin_studio_assistant_mcp_event,append_studio_assistant_mcp_event
    use studio_assistant_transport, only: register_assistant_worker
    implicit none(type, external)
    type(QApplication) :: app
    type(QScreen) :: screen
    type(studio_t),target :: workspace
    type(launch_options) :: options
    character(:),allocatable :: session_path
    integer :: status

    call dispatch_studio_process_mode()
    call read_launch_options(options)
    call initialize_state(workspace)
    app=QApplication()
    call register_assistant_worker()
    call QApplication_setOrganizationName('ForCAD')
    call QApplication_setApplicationName(studio_product_name)
    call build_ui(workspace,empty_document=.true.)
    call bind_studio_mcp(workspace,execute_studio_mcp,append_studio_assistant_mcp_event,begin_studio_assistant_mcp_event)

    if(.not.options%fixed_size.and..not.allocated(options%screenshot))then
        screen=QGuiApplication_primaryScreen()
        if(.not.screen%isNull())then
            options%width=max(studio_min_width,min(options%width,screen%availableWidth()-32))
            options%height=max(studio_min_height,min(options%height,screen%availableHeight()-52))
        end if
    end if
    call workspace%window%resize(options%width,options%height)
    workspace%remember_workspace=options%remember_workspace
    call workspace%window%closeEvent%connect(on_studio_close,workspace)
    if(options%home)call show_home(workspace)
    call workspace%window%show()
    if(allocated(options%project_path))then
        call open_project_path(workspace,options%project_path)
    else if(.not.workspace%home_active)then
        call apply_preview(workspace,options%preview)
    end if
    if(.not.workspace%home_active.and..not.allocated(options%screenshot))then
        workspace%document%document_started=.true.
        call mark_project_saved(workspace)
    end if
    if(options%enable_mcp)then
        session_path=start_studio_mcp()
        if(len(session_path)>0)then
            write(error_unit,'(a)') 'FORCAD_MCP_SESSION='//session_path
        else
            write(error_unit,'(a)') 'Could not enable local MCP control.'
        end if
        flush(error_unit)
    end if

    call schedule_preview(workspace,options%preview,options%screenshot)
    status=app%exec()
    call workspace%viewport_frame%delete()
    call finish_display_refinement(workspace)
    call finish_studio_assistant()
    call finish_studio_mcp()
    call finish_history_storage()
    if(status/=0)stop 1
end program forcad_studio
