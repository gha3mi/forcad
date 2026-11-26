module studio_preview
    use qt, only: QApplication_processEvents,QTabWidget
    use studio_state, only: studio_t,dp
    use studio_ui, only: invoke_command,new_studio_project,select_entity,show_home,show_ready
    use studio_live_viewport, only: fit_around_panels,zoom_live_viewport
    use studio_commands, only: sample_project_command,canvas_command,toggle_parameter_panels_command, &
        top_view_command,wireframe_display_command,show_inspector_command,show_assistant_command,clean_display_command, &
        project_ribbon_tab,create_ribbon_tab,model_ribbon_tab,patch_ribbon_tab,multipatch_ribbon_tab,view_ribbon_tab
    implicit none(type,external)
    private
    public :: validate_preview,is_dialog_preview,apply_preview

    character(17),parameter :: ribbon_previews(*)=[character(17):: &
        'ribbon-project','ribbon-create','ribbon-model','ribbon-parameters','ribbon-assembly','ribbon-view']
    integer,parameter :: ribbon_tabs(*)=[project_ribbon_tab,create_ribbon_tab,model_ribbon_tab, &
        patch_ribbon_tab,multipatch_ribbon_tab,view_ribbon_tab]

contains

    subroutine validate_preview(name)
        character(*),intent(in) :: name

        if(any(name==ribbon_previews))return
        select case(name)
        case('home','default','chess','chess-close','chess-parameters','canvas','empty','units','commands', &
            'top','wireframe','inspector','assistant','clean')
            return
        case default
            error stop 'Unknown preview state.'
        end select
    end subroutine validate_preview

    logical function is_dialog_preview(name)
        character(*),intent(in) :: name

        is_dialog_preview=name=='commands'.or.name=='units'
    end function is_dialog_preview

    subroutine apply_preview(s,name)
        type(studio_t),intent(inout),target :: s
        character(*),intent(in) :: name
        type(QTabWidget) :: ribbon
        integer :: index

        select case(name)
        case('home')
            call show_home(s)
            return
        case('empty','units')
            call new_studio_project(s,'')
            return
        end select

        call invoke_command(s,sample_project_command)
        select case(name)
        case('canvas','chess-close')
            if(.not.s%canvas_only)call invoke_command(s,canvas_command)
            if(name=='chess-close')then
                call QApplication_processEvents()
                call fit_around_panels(s)
                call zoom_live_viewport(s,4._dp)
            end if
        case('chess-parameters','ribbon-parameters')
            call select_entity(s,3)
            if(s%focused_surface)call invoke_command(s,toggle_parameter_panels_command)
            call QApplication_processEvents()
            call fit_around_panels(s)
        case('top')
            call invoke_command(s,top_view_command)
        case('wireframe')
            call invoke_command(s,wireframe_display_command)
        case('inspector')
            call select_entity(s,3)
            call invoke_command(s,show_inspector_command)
        case('assistant')
            call invoke_command(s,show_assistant_command)
        case('clean')
            call invoke_command(s,clean_display_command)
        end select

        index=findloc(ribbon_previews,name,dim=1)
        if(index>0)then
            call s%workspace%findChild('studioRibbonTabs',ribbon)
            call ribbon%setCurrentIndex(ribbon_tabs(index))
        end if
        call show_ready(s)
    end subroutine apply_preview
end module studio_preview
