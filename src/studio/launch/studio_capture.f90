module studio_capture
    use qt, only: QDialog,QLineEdit,QPixmap,QRect,QTimer_singleShot,QApplication_processEvents,QApplication_quit
    use studio_display_refinement, only: display_refinement_pending,request_display_refinement
    use studio_state, only: studio_t
    use studio_ui, only: invoke_command
    use studio_preview, only: is_dialog_preview
    use studio_commands, only: units_command,find_command
    implicit none(type, external)
    private
    public :: schedule_preview

    character(:),allocatable :: preview_name,screenshot_path
    integer :: refinement_waits=0
contains
    subroutine schedule_preview(workspace,preview,screenshot)
        type(studio_t),target,intent(inout) :: workspace
        character(*),intent(in) :: preview
        character(:),allocatable,intent(in) :: screenshot

        preview_name=preview
        refinement_waits=0
        if(allocated(screenshot))then
            screenshot_path=screenshot
            call QTimer_singleShot(400,capture_preview,workspace)
        else if(is_dialog_preview(preview_name))then
            call QTimer_singleShot(0,open_tool_dialog,workspace)
        end if
    end subroutine schedule_preview

    subroutine capture_tool_dialog(data)
        class(*),intent(inout),optional :: data
        type(QDialog) :: dialog
        type(QLineEdit) :: query
        type(QPixmap) :: pixmap
        logical :: saved
        select type(s=>data)
        type is(studio_t)
            if(preview_name=='units')then
                call s%window%findChild('studioUnitsDialog',dialog)
                if(dialog%isNull())error stop 'Could not find the Project units dialog'
            else
                call s%window%findChild('studioCommandFinder',dialog)
                if(dialog%isNull())error stop 'Could not find the command dialog'
                call dialog%findChild('studioCommandQuery',query)
                if(query%isNull())error stop 'Could not find the command search field'
                call query%setText('camera')
            end if
            pixmap=dialog%grab()
            saved=pixmap%save(screenshot_path)
            call pixmap%delete()
            if(.not.saved)error stop 'Could not save tool dialog screenshot'
            print '(a)', 'Saved '//screenshot_path
            call dialog%reject()
        end select
    end subroutine capture_tool_dialog

    subroutine open_tool_dialog(data)
        class(*),intent(inout),optional :: data
        select type(s=>data)
        type is(studio_t)
            select case(trim(preview_name))
            case('units')
                call invoke_command(s,units_command)
            case default
                call invoke_command(s,find_command)
            end select
        end select
    end subroutine open_tool_dialog

    subroutine capture_preview(data)
        class(*),intent(inout),optional :: data
        type(QPixmap) :: pix
        type(QRect) :: frame,client
        logical :: saved

        select type(s=>data)
        type is(studio_t)
            if(is_dialog_preview(preview_name))then
                call QTimer_singleShot(100,capture_tool_dialog,s)
                call open_tool_dialog(s)
            else
                call QApplication_processEvents()
                if(preview_name=='chess-close'.and.refinement_waits<100)then
                    call request_display_refinement(s)
                    if(display_refinement_pending(s))then
                        refinement_waits=refinement_waits+1
                        call QTimer_singleShot(40,capture_preview,s)
                        return
                    end if
                end if
                frame=s%window%frameGeometry()
                client=s%window%geometry()
                print '(a,4(i0,1x))', 'Outer frame extents (left top right bottom): ', &
                    client%x-frame%x,client%y-frame%y, &
                    frame%width-client%width-(client%x-frame%x), &
                    frame%height-client%height-(client%y-frame%y)
                pix=s%window%grab()
                saved=pix%save(screenshot_path)
                call pix%delete()
                if(.not.saved)error stop 'Could not save screenshot'
                print '(a)', 'Saved '//screenshot_path
            end if
            call QApplication_quit()
        end select
    end subroutine capture_preview
end module studio_capture
