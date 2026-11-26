module studio_workspace_preferences
    use studio_commands, only: workspace_panels
    use studio_preferences, only: workspace_preferences_t, read_preferences, write_preferences
    use studio_single_patch, only: active_patch_id, set_active_patch
    use studio_state, only: studio_t
    use studio_viewport_camera, only: set_interactive_quality
    implicit none(type, external)
    private
    public :: restore_preferences, save_preferences
contains
    subroutine restore_preferences(workspace,file_name,restored)
        type(studio_t), intent(inout) :: workspace
        character(*), optional, intent(in) :: file_name
        logical, optional, intent(out) :: restored
        type(workspace_preferences_t) :: values
        integer :: i,panel,initial_width,initial_height
        logical :: has_preferences,visible,geometry_restored,left_visible,basis_visible

        values=workspace_settings(workspace)
        call read_preferences(values,has_preferences,file_name)
        if(present(restored))restored=has_preferences
        if(.not.has_preferences)return

        workspace%entity_kind=values%entity_kind
        do i=1,size(values%sample_ids)
            call set_active_patch(i,values%sample_ids(i),catalog=workspace%document%patches)
        end do
        workspace%knot_axis=values%knot_axis
        workspace%basis_derivative=values%basis_derivative
        workspace%surface_tool=values%surface_tool
        workspace%layers=values%layers
        workspace%domain_layers=values%domain_layers
        workspace%basis_legend=values%basis_legend
        workspace%basis_probe_guide=values%basis_probe_guide
        workspace%world_grid=values%world_grid
        workspace%viewport_style=values%viewport_style
        call set_interactive_quality(workspace,values%interactive_quality)
        left_visible=.false.
        basis_visible=.false.
        do i=1,size(workspace_panels)
            panel=workspace_panels(i)
            visible=values%visible_panels(i).or.panel==4
            call workspace%panels(panel)%setVisible(visible)
            if(panel==1)left_visible=visible
            if(panel==6)basis_visible=visible
            if(panel==4)call workspace%panel_actions(panel)%setChecked(visible)
            if(panel==8)workspace%widgets%inspector_requested=visible
        end do
        call workspace%left_column%setVisible(left_visible.and..not.workspace%canvas_only)
        call workspace%lower_row%setVisible(basis_visible.and..not.workspace%focused_surface)
        if(allocated(values%window_geometry))then
            initial_width=workspace%window%width()
            initial_height=workspace%window%height()
            geometry_restored=workspace%window%restoreGeometry(values%window_geometry)
            if(.not.geometry_restored)call workspace%window%resize(initial_width,initial_height)
        end if
        call workspace%columns%setSizes(values%columns)
        call workspace%center_rows%setSizes(values%center_rows)
    end subroutine restore_preferences

    subroutine save_preferences(workspace,file_name)
        type(studio_t), intent(in) :: workspace
        character(*), optional, intent(in) :: file_name
        type(workspace_preferences_t) :: values
        integer :: i,panel

        values=workspace_settings(workspace)
        values%window_geometry=workspace%window%saveGeometry()
        values%save_sample=values%sample_ids<=2
        do i=1,size(workspace_panels)
            panel=workspace_panels(i)
            values%visible_panels(i)=workspace%panel_actions(panel)%isChecked()
            if(panel==8)values%visible_panels(i)=workspace%widgets%inspector_requested
            if(workspace%focused_surface)values%visible_panels(i)=workspace%saved_panels(panel)
        end do
        if(workspace%focused_surface)then
            values%columns=workspace%normal_columns
            values%center_rows=workspace%normal_center
        else
            values%columns=workspace%columns%sizes()
            values%center_rows=workspace%center_rows%sizes()
        end if
        call write_preferences(values,file_name)
    end subroutine save_preferences

    function workspace_settings(workspace) result(values)
        type(studio_t), intent(in) :: workspace
        type(workspace_preferences_t) :: values
        integer :: i

        values%entity_kind=workspace%entity_kind
        values%sample_ids=[(active_patch_id(i,catalog=workspace%document%patches),i=1,3)]
        values%knot_axis=workspace%knot_axis
        values%basis_derivative=workspace%basis_derivative
        values%surface_tool=workspace%surface_tool
        values%layers=workspace%layers
        values%domain_layers=workspace%domain_layers
        values%basis_legend=workspace%basis_legend
        values%basis_probe_guide=workspace%basis_probe_guide
        values%world_grid=workspace%world_grid
        values%viewport_style=workspace%viewport_style
        values%interactive_quality=workspace%interactive_quality
    end function workspace_settings
end module studio_workspace_preferences
