module studio_diagnostics
    use studio_presentation_text, only: integer_text, without_nul
    use studio_scene, only: project_patch_exists
    use studio_single_patch, only: active_patch_id
    use studio_state, only: studio_t, workspace_diagnostic_t, studio_mesh_valid
    use studio_workspace_selection, only: selected_group, selected_mesh, workspace_selection_key
    implicit none(type, external)
    private

    public :: studio_diagnostic_summary, record_studio_diagnostic

contains

    subroutine record_studio_diagnostic(s,message)
        type(studio_t), intent(inout) :: s
        character(*), intent(in) :: message

        if(len_trim(message)==0)then
            s%diagnostic=workspace_diagnostic_t()
            return
        end if
        s%diagnostic%message=trim(message)
        s%diagnostic%selection_key=workspace_selection_key(s)
        s%diagnostic%document_revision=s%document%document_revision
    end subroutine record_studio_diagnostic

    function studio_diagnostic_summary(s) result(text)
        type(studio_t),intent(in) :: s
        character(:),allocatable :: text
        integer :: mesh,group,i

        text=''
        mesh=selected_mesh(s)
        group=selected_group(s)
        if(mesh>0)then
            if(.not.studio_mesh_valid(s,mesh)) then
                text='The mesh no longer matches its source. Check the source and remesh it.'
            end if
        else if(group>0)then
            associate(validation=>s%document%multipatches%groups(group)%validation)
                if(.not.validation%valid.and.allocated(validation%message))text=validation%message
                if(allocated(validation%joins))then
                    do i=1,size(validation%joins)
                        if(validation%joins(i)%valid.or..not.allocated(validation%joins(i)%message))cycle
                        if(len(text)>650)exit
                        call add_issue('Interface '//integer_text(i)//' · '//validation%joins(i)%message)
                    end do
                end if
            end associate
        else if(s%live_scene%valid)then
            if(allocated(s%live_probe%diagnostic))then
                call add_issue(s%live_probe%diagnostic)
            else if(.not.s%live_probe%valid.and.allocated(s%live_probe%message))then
                call add_issue(s%live_probe%message)
            end if
        else if(allocated(s%live_scene%diagnostic))then
            if(project_patch_exists(s%document%show_templates,s%entity_kind,active_patch_id(s%entity_kind,catalog=s%document%patches),patches=s%document%patches)) then
                call add_issue(s%live_scene%diagnostic)
            end if
        end if
        if(allocated(s%diagnostic%message))then
            if(s%diagnostic%document_revision==s%document%document_revision)then
                if(s%diagnostic%selection_key==workspace_selection_key(s))call add_issue(s%diagnostic%message)
            end if
        end if
        text=without_nul(text)
    contains
        subroutine add_issue(issue)
            character(*),intent(in) :: issue
            if(len_trim(issue)==0)return
            if(index(text,trim(issue))>0)return
            if(len(text)>0)text=text//new_line('a')
            text=text//trim(issue)
        end subroutine add_issue
    end function studio_diagnostic_summary

end module studio_diagnostics
