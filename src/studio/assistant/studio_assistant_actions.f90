! Bridge from named assistant suggestions to the same commands as the ForCAD Studio ribbon.
module studio_assistant_actions
    use studio_commands, only: &
        fit_scene_command, fit_selection_command, reset_controls_command, reset_presentation_command, refine_command, translate_geometry_command, &
        affine_geometry_command, primitive_command, project_isolate_command, project_show_all_command, &
        split_patch_command, extract_patch_command, reverse_patch_command, snap_off_command, &
        snap_vertex_command, snap_closest_command, generate_mesh_command, mesh_edges_command, &
        mesh_pick_command, mesh_only_command, physical_knots_command, mesh_shaded_command, &
        edit_controls_command, undo_patch_command, redo_patch_command, insert_knot_command, &
        remove_knot_command, elevate_degree_command, reduce_degree_command, edit_knot_values_command, &
        edit_interval_command, new_multipatch_command, edit_multipatch_command, validate_multipatch_command, &
        color_command, translucent_command, transform_command, snap_command, &
        canvas_command, pick_mode_base, back_view_command, bottom_view_command, &
        left_view_command, right_view_command, toggle_patch_visibility_command, isometric_view_command, &
        zoom_in_command, zoom_out_command, restore_layout_command, front_view_command, &
        top_view_command, select_tool_command, orbit_tool_command, reset_probe_command, &
        basis_u_command, basis_v_command, basis_w_command, basis_du_command, &
        basis_dv_command, basis_dw_command, knots_u_command, knots_v_command, &
        knots_w_command, show_knots_command, show_basis_command, show_inspector_command, &
        toggle_patch_inspection_command, select_curve_command, select_surface_command, select_volume_command, &
        clean_display_command, control_structure_command, wireframe_display_command, fast_shaded_display_command, &
        shaded_edges_display_command, shaded_display_command, rotate_geometry_command, &
        scale_geometry_command, extrude_command, revolve_command, sweep_command, &
        loft_command, fit_around_panels_command, box_selection_command, lasso_selection_command, pick_controls_command, pick_multipatch_command, &
        clear_selection_command, scene_select_mesh_source_command, domain_reset_command, domain_probe_command, &
        domain_zoom_in_command, domain_zoom_out_command, domain_details_command, basis_functions_command, &
        basis_derivative_command, basis_legend_command, basis_probe_command, toggle_panel_command_base, &
        shaded_layer_command, wireframe_layer_command, control_net_layer_command, project_duplicate_command, &
        world_grid_command, viewport_layer_base, parameter_lines_layer, project_expand_command, project_collapse_command, interactive_quality_command, &
        scene_reveal_command, scene_exit_isolation_command, inspect_patch_command, parameter_grid_command, basis_support_command, &
        domain_select_command, domain_insert_command, domain_remove_command, transform_done_command, extract_section_command, &
        extract_face_command, mesh_patch_counts_command, mesh_target_size_command, remove_multipatch_command, &
        insert_assembly_knot_command,elevate_assembly_degree_command,copy_geometry_command,paste_geometry_command
    implicit none(type, external)
    private
    public :: assistant_action_context, assistant_catalog, resolve_assistant_action
    public :: assistant_action, available_assistant_actions

    ! The context describes the current workspace and active parameter direction.
    type :: assistant_action_context
        logical :: has_curve=.false., has_surface=.false., has_volume=.false.
        logical :: has_visible_geometry=.false., has_selection=.false., can_fit_selection=.false.
        logical :: can_copy=.false., can_paste=.false.
        logical :: has_multipatch=.false., can_create_multipatch=.false., can_edit_assembly_controls=.false.
        logical :: has_mesh=.false., has_selected_mesh=.false., has_visible_mesh=.false.
        logical :: has_shaded_mesh=.false., mesh_only=.false.
        logical :: can_edit=.false., can_remove_knot=.false., can_reduce_degree=.false.
        logical :: can_restructure_patch=.false., can_extract_patch=.false., can_undo=.false., can_redo=.false.
        logical :: can_sweep=.false., can_loft=.false.
        logical :: has_isolation=.false., domain_knots_active=.false., can_remove_active_knot=.false.
        logical :: can_transform_gizmo=.false., transform_active=.false.
    end type assistant_action_context

    type :: assistant_action
        character(len=26) :: key
        character(len=54) :: label
        integer, private :: command
        integer, private :: min_dimension = 0
        logical :: dialog_free = .false.
    end type assistant_action

    type(assistant_action), parameter :: actions(*) = [ &
        assistant_action('workspace.reset','Reset workspace controls',reset_controls_command,dialog_free=.true.), &
        assistant_action('layout.reset-defaults','Reset the saved workspace layout and controls',reset_presentation_command,dialog_free=.true.), &
        assistant_action('workspace.curve','Select a curve in the scene',select_curve_command,dialog_free=.true.), &
        assistant_action('workspace.surface','Select a surface in the scene',select_surface_command,dialog_free=.true.), &
        assistant_action('workspace.volume','Select a volume in the scene',select_volume_command,dialog_free=.true.), &
        assistant_action('camera.isometric','Reset to isometric view',isometric_view_command,dialog_free=.true.), &
        assistant_action('camera.fit','Frame visible objects',fit_scene_command,dialog_free=.true.), &
        assistant_action('camera.fit-selection','Frame the visible selection',fit_selection_command,dialog_free=.true.), &
        assistant_action('camera.fit-panels','Frame the scene around open panels',fit_around_panels_command,dialog_free=.true.), &
        assistant_action('camera.front','Use front view',front_view_command,dialog_free=.true.), &
        assistant_action('camera.back','Use back view',back_view_command,dialog_free=.true.), &
        assistant_action('camera.left','Use left view',left_view_command,dialog_free=.true.), &
        assistant_action('camera.right','Use right view',right_view_command,dialog_free=.true.), &
        assistant_action('camera.top','Use top view',top_view_command,dialog_free=.true.), &
        assistant_action('camera.bottom','Use bottom view',bottom_view_command,dialog_free=.true.), &
        assistant_action('camera.zoom-in','Zoom in',zoom_in_command,dialog_free=.true.), &
        assistant_action('camera.zoom-out','Zoom out',zoom_out_command,dialog_free=.true.), &
        assistant_action('layout.toggle-focus','Toggle parameter panels',toggle_patch_inspection_command,1,dialog_free=.true.), &
        assistant_action('layout.canvas','Toggle Canvas only',canvas_command,dialog_free=.true.), &
        assistant_action('layout.restore','Show all workspace panels',restore_layout_command,dialog_free=.true.), &
        assistant_action('layout.scene','Toggle the Scene panel',toggle_panel_command_base+1,dialog_free=.true.), &
        assistant_action('layout.domain','Toggle Parameter Space',toggle_panel_command_base+5,1,dialog_free=.true.), &
        assistant_action('layout.basis','Toggle the basis plot panel',toggle_panel_command_base+6,1,dialog_free=.true.), &
        assistant_action('layout.inspector','Toggle the Inspector',toggle_panel_command_base+8,dialog_free=.true.), &
        assistant_action('patch.inspect','Show Parameter Space and Basis Plot',inspect_patch_command,1,dialog_free=.true.), &
        assistant_action('scene.reveal','Reveal the current selection in the Scene tree',scene_reveal_command,dialog_free=.true.), &
        assistant_action('scene.expand','Expand all Scene tree branches',project_expand_command,dialog_free=.true.), &
        assistant_action('scene.collapse','Collapse all Scene tree branches',project_collapse_command,dialog_free=.true.), &
        assistant_action('scene.exit-isolation','Exit Scene isolation',scene_exit_isolation_command,dialog_free=.true.), &
        assistant_action('display.clean','Use clean preview',clean_display_command,dialog_free=.true.), &
        assistant_action('display.structure','Show the control structure',control_structure_command,1,dialog_free=.true.), &
        assistant_action('display.shaded-faces','Show filled faces',shaded_display_command,dialog_free=.true.), &
        assistant_action('display.shaded-edges','Show filled faces and geometry edges',shaded_edges_display_command,dialog_free=.true.), &
        assistant_action('display.fast-shaded','Use Fast Shaded mode',fast_shaded_display_command,dialog_free=.true.), &
        assistant_action('display.wireframe','Show geometry and mesh edges without fills',wireframe_display_command,dialog_free=.true.), &
        assistant_action('display.shaded','Toggle shaded geometry',shaded_layer_command,dialog_free=.true.), &
        assistant_action('display.wire-layer','Toggle the wire overlay',wireframe_layer_command,dialog_free=.true.), &
        assistant_action('display.control-net','Toggle control net visibility',control_net_layer_command,1,dialog_free=.true.), &
        assistant_action('display.control-points','Toggle control point visibility',viewport_layer_base+3,dialog_free=.true.), &
        assistant_action('display.parameter-lines','Toggle parameter guides on the geometry',viewport_layer_base+parameter_lines_layer-1,dialog_free=.true.), &
        assistant_action('display.ground-grid','Toggle the 3D ground grid',world_grid_command,dialog_free=.true.), &
        assistant_action('display.interactive','Toggle Interactive Quality during zoom',interactive_quality_command,dialog_free=.true.), &
        assistant_action('display.color','Open the object color picker',color_command), &
        assistant_action('display.translucency','Toggle object translucency',translucent_command,dialog_free=.true.), &
        assistant_action('display.toggle-visibility','Show or hide the selection',toggle_patch_visibility_command,dialog_free=.true.), &
        assistant_action('display.isolate','Isolate the selection',project_isolate_command,dialog_free=.true.), &
        assistant_action('display.show-all','Show all geometry and meshes',project_show_all_command,dialog_free=.true.), &
        assistant_action('display.physical-knots','Toggle knot boundaries',physical_knots_command,dialog_free=.true.), &
        assistant_action('view.select','Use the selection tool',select_tool_command,dialog_free=.true.), &
        assistant_action('view.box','Use box selection in the viewport',box_selection_command,dialog_free=.true.), &
        assistant_action('view.lasso','Use lasso selection in the viewport',lasso_selection_command,dialog_free=.true.), &
        assistant_action('view.orbit','Use the orbit tool',orbit_tool_command,dialog_free=.true.), &
        assistant_action('view.move-gizmo','Toggle the Move handles in the 3D view',translate_geometry_command,dialog_free=.true.), &
        assistant_action('view.rotate-gizmo','Toggle the Rotate handles in the 3D view',rotate_geometry_command,dialog_free=.true.), &
        assistant_action('view.scale-gizmo','Toggle the Scale handles in the 3D view',scale_geometry_command,dialog_free=.true.), &
        assistant_action('view.end-transform','Close the 3D transform handles',transform_done_command,dialog_free=.true.), &
        assistant_action('view.snap','Toggle snapping',snap_command,dialog_free=.true.), &
        assistant_action('view.snap-off','Turn snapping off',snap_off_command,dialog_free=.true.), &
        assistant_action('view.snap-vertex','Snap to vertices',snap_vertex_command,dialog_free=.true.), &
        assistant_action('view.snap-closest','Snap to the closest geometry point',snap_closest_command,dialog_free=.true.), &
        assistant_action('pick.automatic','Use automatic viewport picking',pick_mode_base,dialog_free=.true.), &
        assistant_action('pick.faces','Pick faces in the viewport',pick_mode_base+1,dialog_free=.true.), &
        assistant_action('pick.edges','Pick edges in the viewport',pick_mode_base+2,dialog_free=.true.), &
        assistant_action('pick.vertices','Pick vertices in the viewport',pick_mode_base+3,dialog_free=.true.), &
        assistant_action('pick.controls','Pick control points in the viewport',pick_controls_command,dialog_free=.true.), &
        assistant_action('pick.multipatches','Pick whole multipatches in the viewport',pick_multipatch_command,dialog_free=.true.), &
        assistant_action('clipboard.copy','Copy selected objects from this project',copy_geometry_command,dialog_free=.true.), &
        assistant_action('clipboard.paste','Paste objects copied in the assigned folder',paste_geometry_command,dialog_free=.true.), &
        assistant_action('selection.clear','Clear the viewport selection',clear_selection_command,dialog_free=.true.), &
        assistant_action('model.transform','Open the numeric Transform editor',transform_command), &
        assistant_action('model.move','Open the Move tool',translate_geometry_command), &
        assistant_action('model.rotate','Open the Rotate tool',rotate_geometry_command), &
        assistant_action('model.scale','Open the Scale tool',scale_geometry_command), &
        assistant_action('model.affine','Open the numeric Affine Transform editor',affine_geometry_command), &
        assistant_action('model.create','Open the geometry creation dialog',primitive_command), &
        assistant_action('model.point','Open the Point creation dialog',primitive_command+1), &
        assistant_action('model.line','Open the Line creation dialog',primitive_command+2), &
        assistant_action('model.circle','Open the Circle creation dialog',primitive_command+3), &
        assistant_action('model.arc','Open the Arc creation dialog',primitive_command+4), &
        assistant_action('model.rectangle','Open the Rectangle creation dialog',primitive_command+5), &
        assistant_action('model.box','Open the Box creation dialog',primitive_command+6), &
        assistant_action('model.cylinder','Open the Cylinder surface dialog',primitive_command+7), &
        assistant_action('model.tube','Open the Tube volume dialog',primitive_command+8), &
        assistant_action('model.extrude','Open the Extrude construction dialog',extrude_command), &
        assistant_action('model.revolve','Open the Revolve construction dialog',revolve_command), &
        assistant_action('model.sweep','Open the Sweep construction dialog',sweep_command), &
        assistant_action('model.loft','Open the Loft construction dialog',loft_command), &
        assistant_action('model.split','Open the Split Patch editor',split_patch_command,1), &
        assistant_action('model.extract','Open the Extract Parameter Range editor',extract_patch_command,1), &
        assistant_action('model.extract-section','Open the Extract Section editor',extract_section_command,2), &
        assistant_action('model.extract-face','Open the Extract Face editor',extract_face_command,3), &
        assistant_action('model.reverse','Open the Reverse Parameter Direction editor',reverse_patch_command,1), &
        assistant_action('model.refine','Open the Refine Patch editor',refine_command,1), &
        assistant_action('model.duplicate','Duplicate the selected patch',project_duplicate_command,1,dialog_free=.true.), &
        assistant_action('model.undo','Undo the latest modeling edit',undo_patch_command), &
        assistant_action('model.redo','Redo the latest undone edit',redo_patch_command), &
        assistant_action('mesh.generate','Generate the structured mesh',generate_mesh_command), &
        assistant_action('mesh.cell-counts','Open the per-patch Mesh Cell Counts editor',mesh_patch_counts_command), &
        assistant_action('mesh.target-size','Open the Mesh Target Size editor',mesh_target_size_command), &
        assistant_action('mesh.toggle-visibility','Show or hide meshes that belong to the selection',mesh_edges_command,dialog_free=.true.), &
        assistant_action('mesh.shaded','Toggle shaded mesh faces',mesh_shaded_command,dialog_free=.true.), &
        assistant_action('mesh.only','Toggle Meshes only',mesh_only_command,dialog_free=.true.), &
        assistant_action('mesh.pick','Toggle cell picking on visible meshes',mesh_pick_command,dialog_free=.true.), &
        assistant_action('mesh.select-source','Select the source of the selected mesh',scene_select_mesh_source_command,dialog_free=.true.), &
        assistant_action('multipatch.create','Open the Multipatch assembly editor',new_multipatch_command), &
        assistant_action('multipatch.edit','Edit the selected multipatch assembly',edit_multipatch_command), &
        assistant_action('multipatch.dissolve','Open the Dissolve Multipatch confirmation',remove_multipatch_command), &
        assistant_action('multipatch.validate','Validate multipatch interfaces',validate_multipatch_command,dialog_free=.true.), &
        assistant_action('multipatch.insert-knot','Open multipatch knot insertion',insert_assembly_knot_command), &
        assistant_action('multipatch.elevate-degree','Open multipatch degree elevation',elevate_assembly_degree_command), &
        assistant_action('inspector.show','Show the Inspector',show_inspector_command,dialog_free=.true.), &
        assistant_action('controls.edit','Open the control point and weight editor',edit_controls_command), &
        assistant_action('probe.reset','Reset probe coordinates',reset_probe_command,1,dialog_free=.true.), &
        assistant_action('domain.probe','Switch Parameter Space to the probe tool',domain_probe_command,1,dialog_free=.true.), &
        assistant_action('domain.reset-view','Reset the parameter view',domain_reset_command,1,dialog_free=.true.), &
        assistant_action('domain.zoom-in','Zoom in on Parameter Space',domain_zoom_in_command,1,dialog_free=.true.), &
        assistant_action('domain.zoom-out','Zoom out of Parameter Space',domain_zoom_out_command,1,dialog_free=.true.), &
        assistant_action('domain.grid','Toggle the Parameter Space grid',parameter_grid_command,1,dialog_free=.true.), &
        assistant_action('domain.support','Toggle basis support in Parameter Space',basis_support_command,1,dialog_free=.true.), &
        assistant_action('knots.show','Edit knots in parameter space',show_knots_command,1,dialog_free=.true.), &
        assistant_action('knots.groups','Toggle the knot groups table',domain_details_command,1,dialog_free=.true.), &
        assistant_action('knots.select-tool','Activate the Select Knot tool',domain_select_command,1,dialog_free=.true.), &
        assistant_action('knots.insert-tool','Activate the Insert Knot tool',domain_insert_command,1,dialog_free=.true.), &
        assistant_action('knots.remove-tool','Activate the Remove Knot tool',domain_remove_command,1,dialog_free=.true.), &
        assistant_action('knots.insert','Open the Insert Knot editor',insert_knot_command,1), &
        assistant_action('knots.remove','Open the Remove Knot editor',remove_knot_command,1), &
        assistant_action('knots.values','Open the Knot Values editor',edit_knot_values_command,1), &
        assistant_action('knots.domain','Open the Parameter Interval editor',edit_interval_command,1), &
        assistant_action('degree.elevate','Open the Elevate Degree editor',elevate_degree_command,1), &
        assistant_action('degree.reduce','Open the Reduce Degree editor',reduce_degree_command,1), &
        assistant_action('basis.show','Show the basis plot panel',show_basis_command,1,dialog_free=.true.), &
        assistant_action('basis.functions','Show basis functions on the active axis',basis_functions_command,1,dialog_free=.true.), &
        assistant_action('basis.derivative','Show the first derivative on the active axis',basis_derivative_command,1,dialog_free=.true.), &
        assistant_action('basis.legend','Toggle the basis plot legend',basis_legend_command,1,dialog_free=.true.), &
        assistant_action('basis.probe','Toggle the basis plot probe guide',basis_probe_command,1,dialog_free=.true.), &
        assistant_action('basis.u','Plot the U basis functions',basis_u_command,1,dialog_free=.true.), &
        assistant_action('basis.v','Plot the V basis functions',basis_v_command,2,dialog_free=.true.), &
        assistant_action('basis.w','Plot the W basis functions',basis_w_command,3,dialog_free=.true.), &
        assistant_action('basis.du','Plot the U basis derivatives',basis_du_command,1,dialog_free=.true.), &
        assistant_action('basis.dv','Plot the V basis derivatives',basis_dv_command,2,dialog_free=.true.), &
        assistant_action('basis.dw','Plot the W basis derivatives',basis_dw_command,3,dialog_free=.true.), &
        assistant_action('knots.u','Select the U knot vector',knots_u_command,1,dialog_free=.true.), &
        assistant_action('knots.v','Select the V knot vector',knots_v_command,2,dialog_free=.true.), &
        assistant_action('knots.w','Select the W knot vector',knots_w_command,3,dialog_free=.true.) ]
contains
    function available_assistant_actions(entity_kind,mesh_available,context) result(entries)
        integer, intent(in) :: entity_kind
        logical, intent(in), optional :: mesh_available
        type(assistant_action_context), intent(in), optional :: context
        type(assistant_action), allocatable :: entries(:)
        integer :: i

        entries=pack(actions,[(available_for_context(actions(i),entity_kind,mesh_available,context),i=1,size(actions))])
    end function available_assistant_actions

    function assistant_catalog(entity_kind,mesh_available,context) result(text)
        integer, intent(in) :: entity_kind
        logical, intent(in), optional :: mesh_available
        type(assistant_action_context), intent(in), optional :: context
        character(:), allocatable :: text
        type(assistant_action), allocatable :: entries(:)
        integer :: i

        entries=available_assistant_actions(entity_kind,mesh_available,context)
        text=''
        do i=1,size(entries)
            text=text//trim(entries(i)%key)//' — '//trim(entries(i)%label)//new_line('a')
        end do
    end function assistant_catalog

    subroutine resolve_assistant_action(key,entity_kind,command,label,mesh_available,context,dialog_free_only)
        character(len=*), intent(in) :: key
        integer, intent(in) :: entity_kind
        integer, intent(out) :: command
        character(:), allocatable, intent(out) :: label
        logical, intent(in), optional :: mesh_available
        type(assistant_action_context), intent(in), optional :: context
        logical, intent(in), optional :: dialog_free_only
        integer :: i

        command=0
        label=''
        do i=1,size(actions)
            if(key/=trim(actions(i)%key))cycle
            if(present(dialog_free_only))then
                if(dialog_free_only.and..not.actions(i)%dialog_free)return
            end if
            if(.not.available_for_context(actions(i),entity_kind,mesh_available,context))return
            command=actions(i)%command
            label=trim(actions(i)%label)
            return
        end do
    end subroutine resolve_assistant_action

    pure logical function available_for_context(action,entity_kind,mesh_available,context) result(available)
        type(assistant_action), intent(in) :: action
        integer, intent(in) :: entity_kind
        logical, intent(in), optional :: mesh_available
        type(assistant_action_context), intent(in), optional :: context
        type(assistant_action_context) :: state
        logical :: patch, profile, geometry

        available=entity_kind>=action%min_dimension.and.entity_kind<=3
        if(.not.available)return
        if(present(context))state=context
        patch=state%can_edit.and..not.state%has_multipatch
        profile=state%has_curve.or.state%has_surface
        geometry=profile.or.state%has_volume
        select case(action%command)
        case(select_curve_command)
            if(present(context))available=state%has_curve
        case(select_surface_command)
            if(present(context))available=state%has_surface
        case(select_volume_command)
            if(present(context))available=state%has_volume
        case(fit_scene_command,fit_around_panels_command)
            if(present(context))available=(state%has_visible_geometry.and..not.state%mesh_only).or.state%has_visible_mesh
        case(fit_selection_command)
            available=state%can_fit_selection
        case(scene_reveal_command)
            available=state%has_selection
        case(scene_exit_isolation_command)
            available=state%has_isolation
        case(copy_geometry_command)
            available=state%can_copy
        case(paste_geometry_command)
            available=state%can_paste
        case(project_duplicate_command)
            available=patch
        case(scene_select_mesh_source_command)
            if(present(context))available=state%has_selected_mesh
        case(domain_probe_command,domain_reset_command,domain_zoom_in_command,domain_zoom_out_command,domain_details_command, &
             basis_functions_command,basis_derivative_command,basis_legend_command,basis_probe_command, &
             parameter_grid_command,basis_support_command,toggle_panel_command_base+5,toggle_panel_command_base+6, &
             inspect_patch_command,toggle_patch_inspection_command)
            if(present(context))available=patch
        case(domain_select_command,domain_insert_command)
            available=patch.and.state%domain_knots_active
        case(domain_remove_command)
            available=patch.and.state%domain_knots_active.and.state%can_remove_active_knot
        case(color_command,transform_command)
            if(present(context))available=patch.or.state%has_multipatch
        case(translucent_command,affine_geometry_command)
            available=patch.or.state%has_multipatch
        case(translate_geometry_command,rotate_geometry_command,scale_geometry_command)
            available=patch.or.state%has_multipatch
            if(action%dialog_free)available=state%can_transform_gizmo
        case(transform_done_command)
            available=state%transform_active
        case(toggle_patch_visibility_command,project_isolate_command)
            available=state%has_selection
        case(project_show_all_command)
            available=geometry.or.state%has_mesh.or.state%has_visible_mesh
        case(physical_knots_command)
            available=state%has_visible_geometry.and..not.state%mesh_only
        case(reset_probe_command,basis_u_command:basis_w_command,knots_u_command:knots_w_command, &
             show_knots_command,show_basis_command,control_structure_command)
            if(present(context))available=patch
        case(edit_controls_command)
            available=patch.or.state%can_edit_assembly_controls
        case(insert_assembly_knot_command,elevate_assembly_degree_command)
            available=state%can_edit_assembly_controls
        case(pick_multipatch_command)
            if(present(context))available=.not.state%mesh_only
        case(basis_du_command:basis_dw_command,insert_knot_command,elevate_degree_command,edit_interval_command,refine_command)
            available=patch
        case(remove_knot_command,edit_knot_values_command)
            available=patch.and.state%can_remove_knot
        case(reduce_degree_command)
            available=patch.and.state%can_reduce_degree
        case(split_patch_command,reverse_patch_command)
            available=patch.and.state%can_restructure_patch
        case(extract_patch_command,extract_section_command,extract_face_command)
            available=patch.and.state%can_extract_patch
        case(undo_patch_command)
            available=state%can_undo
        case(redo_patch_command)
            available=state%can_redo
        case(extrude_command,revolve_command)
            available=profile
        case(sweep_command)
            available=state%can_sweep
        case(loft_command)
            available=state%can_loft
        case(new_multipatch_command)
            available=state%can_create_multipatch
        case(edit_multipatch_command,remove_multipatch_command,mesh_patch_counts_command)
            available=state%has_multipatch
        case(validate_multipatch_command)
            available=entity_kind==0
            if(present(context))available=state%has_multipatch
        case(generate_mesh_command,mesh_target_size_command)
            available=entity_kind>=1
            if(present(context))available=patch.or.state%has_multipatch
            if(present(mesh_available))available=mesh_available
        case(mesh_edges_command)
            available=state%has_mesh
        case(mesh_pick_command)
            available=state%has_visible_mesh
        case(mesh_only_command)
            available=state%has_visible_mesh.or.state%mesh_only
        case(mesh_shaded_command)
            available=state%has_shaded_mesh
        case default
            continue
        end select
    end function available_for_context
end module studio_assistant_actions
