submodule(studio_assistant) studio_assistant_draft_ui
    use studio_chat_attachment, only: read_studio_attachment
    use studio_platform_files, only: file_metadata,inspect_path,restrict_file
    use qt, only: QFontMetrics
    implicit none(type, external)
contains

    module procedure on_attach_selection
        character(:), allocatable :: label,details

        select type(state=>data)
        type is(assistant_t)
            if(state%busy.or..not.associated(state%workspace))return
            call capture_assistant_selection(state%workspace,label,details)
            if(len(details)==0)then
                call show_assistant_status(state,'Make a selection in Scene, the viewport, or Parameter Space first.')
                return
            end if
            state%attached_selection=details
            state%attached_selection_label=label
            state%attached_selection_document=state%workspace%document%document_revision
            call state%panel%attached_context%setText('Selection · '//label)
            call set_status_hint(state%panel%attached_context,details)
            call state%panel%context_bar%show()
            call show_assistant_status(state,'Selection attached.')
        end select
    end procedure on_attach_selection

    module procedure on_clear_draft_selection
        select type(state=>data)
        type is(assistant_t)
            if(.not.state%busy)call clear_draft_selection(state)
        end select
    end procedure on_clear_draft_selection

    module procedure clear_draft_selection
        if(allocated(state%attached_selection))deallocate(state%attached_selection)
        if(allocated(state%attached_selection_label))deallocate(state%attached_selection_label)
        state%attached_selection_document=-1_int64
        if(.not.state%panel%context_bar%isNull())call state%panel%context_bar%hide()
    end procedure clear_draft_selection

    module procedure on_attach_view
        type(studio_attachment) :: candidate
        type(file_metadata) :: copy_directory
        character(:),allocatable :: path,error,project,folder
        logical :: removed
        integer :: added

        select type(state=>data)
        type is(assistant_t)
            if(state%busy)return
            if(.not.ensure_project_folder(state))return
            project=state%workspace%document%project_directory
            folder=temporary_private_directory(studio_project_directory(project,'screenshots'),'viewport-')
            if(len(folder)>0)then
                path=folder//'/Viewport.png'
                if(inspect_path(folder,copy_directory).and.save_viewport_image(state%workspace,path,max_dimension=1440))then
                    if(restrict_file(path//c_null_char)==0)then
                        call read_studio_attachment(path,candidate,error)
                        if(allocated(candidate%name))then
                            candidate%project_copy=folder
                            candidate%copy_directory=copy_directory
                            added=append_draft_attachments(state,[candidate],project)
                            if(added>0)then
                                call show_assistant_status(state,'Screenshot attached.')
                                call update_controls(state)
                                call state%panel%attachment_list%setCurrentRow(size(state%attachments)-1)
                                return
                            end if
                        end if
                    end if
                end if
                removed=QDir_removeRecursively(folder)
            end if
            call show_assistant_status(state,'Could not attach a screenshot. Check that the project folder is writable.',warning=.true.)
        end select
    end procedure on_attach_view

    module procedure sync_draft_attachment
        type(QFontMetrics) :: metrics
        character(:),allocatable :: mode
        integer :: i,vision,row,count

        if(state%panel%attachment_list%isNull())return
        if(.not.allocated(state%attachments))then
            call state%panel%attachment_list%clear()
            call state%panel%attachment_bar%hide()
            return
        end if
        row=state%panel%attachment_list%currentRow()
        count=size(state%attachments)
        if(state%panel%attachment_list%count()/=count)call state%panel%attachment_list%clear()
        vision=attachment_vision(state)
        do i=1,count
            associate(file=>state%attachments(i))
                mode='File reference'
                if(allocated(file%text))mode='Text input'
                if(allocated(file%image).and.vision/=0)then
                    mode='Image input'
                    if(vision<0)mode='Image'
                end if
                if(i>state%panel%attachment_list%count())then
                    call state%panel%attachment_list%addItem(file%name//' · '//mode)
                else if(state%panel%attachment_list%itemText(i-1)/=file%name//' · '//mode)then
                    call state%panel%attachment_list%setItemText(i-1,file%name//' · '//mode)
                end if
            end associate
        end do
        call state%panel%attachment_list%setCurrentRow(max(0,min(row,count-1)))
        metrics=QFontMetrics(state%panel%attachment_list%font())
        call state%panel%attachment_list%setFixedHeight(min(3,count)*(metrics%height()+6)+4)
        call metrics%delete()
        call state%panel%attachment_bar%show()
    end procedure sync_draft_attachment

end submodule studio_assistant_draft_ui
