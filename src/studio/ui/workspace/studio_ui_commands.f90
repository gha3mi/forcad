submodule(studio_ui) studio_ui_commands
    use, intrinsic :: iso_c_binding, only: c_int,c_null_char
    use qt, only: QCoreApplication_applicationFilePath,QTemporaryDir
    use studio_platform_files, only: resolved_file_path,last_path_separator,open_directory,current_directory, &
        publish_file,close_file
    use studio_storage_paths, only: studio_project_directory
    use studio_selection, only: selection_set_t
    use studio_document, only: set_document_project_directory
    use studio_commands, only: about_command, head_sample_command
    use studio_help_dialog, only: show_interaction_guide
    implicit none(type, external)
contains
    function sample_path(filename) result(path)
        character(*),intent(in) :: filename
        character(:),allocatable :: path,executable,sample_file
        integer :: separator

        sample_file='samples/'//filename
        path=''
        executable=QCoreApplication_applicationFilePath()
        separator=last_path_separator(executable)
        if(separator>0)then
            path=resolved_file_path(executable(:separator)//'../share/forcad-studio/'//filename)
            if(len(path)==0)path=resolved_file_path(executable(:separator)//'../../../'//sample_file)
        end if
        if(len(path)==0)path=resolved_file_path(sample_file)
    end function sample_path

    subroutine save_workspace_screenshot(s)
        type(studio_t),intent(inout) :: s
        character(*),parameter :: staging_prefix='.capture-'
        type(QTemporaryDir) :: staging
        character(:),allocatable :: folder,path,leaf,staging_path
        character(8) :: date
        character(10) :: time
        integer(c_int) :: parent,closed
        logical :: saved

        if(.not.allocated(s%document%project_directory))then
            folder=choose_project_folder(s%window)
            if(len(folder)==0)return
            call set_document_project_directory(s%document,folder)
            call sync_studio_assistant(s)
        end if
        folder=studio_project_directory(s%document%project_directory,'screenshots')
        if(len(folder)==0)then
            call s%status%showMessage('Could not open the project screenshot folder',6500,warning=.true.)
            return
        end if
        parent=open_directory(current_directory,folder//c_null_char)
        if(parent<0)then
            call s%status%showMessage('Could not open the project screenshot folder',6500,warning=.true.,detail=folder)
            return
        end if
        staging=QTemporaryDir(folder//'/'//staging_prefix//'XXXXXX')
        saved=.false.
        if(staging%isValid())then
            path=staging%filePath('Screenshot.png')
            staging_path=staging%path()
            call date_and_time(date,time)
            leaf='Screenshot-'//date//'-'//time(:6)//'-'//staging_path(last_path_separator(staging_path)+len(staging_prefix)+1:)//'.png'
            if(capture_workspace(s,path))saved=publish_file(path//c_null_char,parent,leaf//c_null_char,0_c_int)==0
        end if
        call staging%delete()
        closed=close_file(parent)
        if(saved)then
            call s%status%showMessage('Screenshot saved · '//leaf,5000,detail=folder//'/'//leaf)
        else
            call s%status%showMessage('Could not save screenshot',6500,warning=.true.,detail=folder)
        end if
    end subroutine save_workspace_screenshot

    module procedure invoke_command
        character(len=:),allocatable::path,text,sample_title
        integer::i,result,axis
        integer(int64) :: selected_multipatch_stamp
        type(selection_target_t) :: focused
        type(selection_set_t) :: saved_selection
        logical::visible,changed
        if(allocated(s%selection%items))then
            if(s%selection%count()/=1)then
                select case(command)
                case(edit_controls_command,insert_knot_command,remove_knot_command,edit_knot_values_command, &
                    edit_interval_command,edit_interval_axis_base+1:edit_interval_axis_base+3,apply_interval_command, &
                    elevate_degree_command,reduce_degree_command, &
                    translate_geometry_command:affine_geometry_command,transform_command,refine_command, &
                    split_patch_command,extract_patch_command,reverse_patch_command,extract_section_command,extract_face_command, &
                    project_rename_command, &
                    project_duplicate_command,edit_multipatch_command,validate_multipatch_command, &
                    remove_multipatch_command,insert_assembly_knot_command,elevate_assembly_degree_command, &
                    color_command,translucent_command,export_iges_command,export_vtk_command)
                    call s%status%showMessage('Select one active item for this edit',3500)
                    return
                end select
            end if
        end if
        if(s%canvas_only)then
            select case(command)
            case(restore_layout_command,show_knots_command,show_basis_command,show_inspector_command, &
                 toggle_parameter_panels_command,toggle_patch_inspection_command,show_assistant_command,inspect_patch_command, &
                 domain_details_command,domain_probe_command)
                call toggle_canvas(s)
                if(command==toggle_parameter_panels_command.and..not.s%focused_surface)return
            case default
                continue
            end select
        end if
        if(s%selected_multipatch>0)then
            select case(command)
            case(22,25,27,30:42,50:51,show_knots_command,show_basis_command, &
                 toggle_panel_command_base+5:toggle_panel_command_base+6, &
                 toggle_parameter_panels_command,toggle_patch_inspection_command, &
                 domain_details_command,domain_probe_command,domain_reset_command,basis_functions_command,basis_derivative_command, &
                 elevate_degree_command,reduce_degree_command,edit_knot_values_command,edit_interval_command, &
                 refine_command,apply_interval_command,domain_zoom_out_command,domain_zoom_in_command, &
                 split_patch_command,extract_patch_command,reverse_patch_command,extract_section_command,extract_face_command, &
                 domain_select_command:domain_remove_command, &
                 200:202,inspect_patch_command,basis_legend_command,basis_probe_command)
                call s%status%showMessage('Select a member to inspect or edit its parameters',3500)
                return
            case default
                continue
            end select
        end if
        select case(command)
        case(reset_controls_command,reset_presentation_command)
            axis=s%entity_kind
            if(allocated(s%selection%items))then
                saved_selection=s%selection
            else
                call saved_selection%change(primary_workspace_target(s),selection_replace)
            end if
            call cancel_area_drag(s)
            call reset_domain_view(s)
            call initialize_state(s)
            s%entity_kind=axis
            s%selected_multipatch=0
            do i=1,size(s%inspector%headers)
                call s%inspector%headers(i)%setChecked(i==1)
                call update_inspector_section(s%inspector,i)
            end do
            if(allocated(saved_selection%items))call apply_workspace_selection_batch(s,saved_selection%items,selection_replace)
            if(command==reset_presentation_command)then
                if(s%remember_workspace)call clear_preferences()
                call restore_workspace_layout(s)
                call s%status%showMessage('Workspace reset to presentation defaults',4000)
            else
                call s%status%showMessage('Workspace controls reset',3000)
            end if
        case(screenshot_command)
            call save_workspace_screenshot(s)
        case(home_command)
            call show_home(s)
            return
        case(resume_project_command)
            call enter_workspace(s,.false.)
            if(s%resume_assistant)call invoke_command(s,show_assistant_command)
            call s%surface%setFocus()
            call show_ready(s)
            return
        case(new_project_command)
            text=choose_project_folder(s%window)
            if(len(text)==0)return
            call start_project_tab(s)
            call new_studio_project(s,text)
            return
        case(open_project_command)
            call choose_project(s)
            return
        case(close_project_tab_command)
            call close_active_project_tab(s)
            return
        case(import_project_command)
            call run_project_import(s)
            return
        case(copy_geometry_command)
            call copy_project_selection(s)
            return
        case(paste_geometry_command)
            call paste_project_selection(s)
            return
        case(save_project_command,save_as_command)
            visible=write_project(s,command==save_as_command)
            return
        case(project_folder_command)
            if(allocated(s%document%project_directory))then
                text=choose_project_folder(s%window,s%document%project_directory)
            else
                text=choose_project_folder(s%window)
            end if
            if(len(text)>0)then
                call set_document_project_directory(s%document,text)
                call sync_studio_assistant(s)
                call s%status%showMessage('Project folder · '//text,6000)
            end if
            return
        case(measure_geometry_command)
            call measure_selection(s)
            return
        case(units_command)
            call choose_units(s)

            return
        case(sample_project_command,head_sample_command)
            if(command==sample_project_command)then
                path=sample_path('chess.fcad')
                sample_title='Chess'
            else
                path=sample_path('head.fcad')
                sample_title='Head'
            end if
            if(len(path)==0)then
                call s%status%showMessage(sample_title//' sample is unavailable',6500,warning=.true.)
                return
            end if
            call open_project_tab(s,path,text,sample_title=sample_title)
            if(len(text)>0)call s%status%showMessage('Could not open the '//sample_title//' sample',6500,warning=.true.,detail=text)
            return
        case(import_iges_command)
            call run_geometry_import(s)
            return
        case(export_iges_command,export_vtk_command)
            call run_geometry_export(s,command==export_vtk_command)
            return
        case(export_fortran_command)
            call run_fortran_export(s)
            return
        case(mesh_patch_counts_command,mesh_target_size_command)
            call run_patch_mesh_dialog(s,command==mesh_target_size_command)
            return
        case(open_mesh_dialog_command)
            call run_patch_mesh_dialog(s,.false.,generate=.true.)
            return
        case(generate_mesh_command)
            call generate_workspace_mesh(s)
            return
        case(mesh_edges_command)
            call toggle_selected_mesh_visibility(s)
            return
        case(mesh_shaded_command)
            if(.not.has_visible_meshes(s,2))return
            s%mesh_shaded=.not.s%mesh_shaded
            call s%widgets%ribbon%mesh_shaded_action%setChecked(s%mesh_shaded)
            call s%surface%update()
            return
        case(physical_knots_command)
            if(s%mesh_only)return
            s%physical_knots=.not.s%physical_knots
            call s%widgets%ribbon%physical_knot_action%setChecked(s%physical_knots)
            call s%surface%update()
            return
        case(mesh_only_command)
            if(.not.s%mesh_only.and..not.has_visible_meshes(s))then
                call s%status%showMessage('Show a mesh in the current view to use Meshes only',3500)
                return
            end if
            s%mesh_only=.not.s%mesh_only
            if(s%mesh_only)then
                if(s%transform_mode>0)call on_transform_exit(s)
                call cancel_transform_drag(s)
                s%control_dragging=.false.
                s%control_snapped=.false.
                s%viewport_feature=0
                call s%status%showMessage('Meshes only',3500)
            else
                call s%status%showMessage('CAD geometry and meshes',3500)
            end if
            call sync_patch_visibility_button(s)
            call refresh_edit_actions(s)
            call s%surface%update()
            return
        case(mesh_pick_command)
            if(.not.has_visible_meshes(s))return
            s%mesh_pick=.not.s%mesh_pick
            if(s%mesh_pick)then
                call cancel_area_drag(s)
                s%area_selection_tool=area_click
                s%surface_tool=1
                if(s%transform_mode>0)call on_transform_exit(s)
                call cancel_transform_drag(s)
                s%viewport_feature=0
                s%control_dragging=.false.
                s%control_snapped=.false.
                call s%status%showMessage('Click a visible mesh cell',3500)
            else
                call s%status%showMessage('Cell picking off',2500)
            end if
            call update_mesh_ribbon(s)
            call repaint_all(s)
            return
        case(export_mesh_command)
            call save_selected_mesh(s)
            return
        case(transform_command)
            call run_geometry_transform(s)
            return
        case(translate_geometry_command:affine_geometry_command-1)
            call activate_transform(s,command-translate_geometry_command+1)
            return
        case(affine_geometry_command)
            call run_geometry_transform(s,3)
            return
        case(transform_exact_base)
            if(s%transform_mode>0)call run_geometry_transform(s,s%transform_mode-1)
            return
        case(transform_exact_base+1:transform_exact_base+3)
            call run_geometry_transform(s,command-transform_exact_base-1)
            return
        case(transform_apply_command)
            call apply_transform_values(s)
            return
        case(transform_done_command)
            call on_transform_exit(s)
            return
        case(color_command,translucent_command)
            call change_material(s,command==color_command)
        case(snap_command,snap_off_command:snap_closest_command)
            call change_snap_mode(s,command)
            return
        case(canvas_command)
            if(s%home_active)return
            call toggle_canvas(s)
        case(clear_selection_command)
            call clear_workspace_selection(s)
            return
        case(pick_mode_base:pick_mode_base+3,pick_controls_command,pick_multipatch_command)
            if(s%transform_mode>0)call on_transform_exit(s)
            call cancel_transform_drag(s)
            s%surface_tool=1
            s%feature_pick_mode=pick_action_index(command)-1
            if(command==pick_controls_command)s%layers(4)=.true.
            if(command/=pick_controls_command.and.s%selected_multipatch>0)then
                s%selected_control=0
                s%canvas_gesture%viewport_shared_hover=0
                if(allocated(s%selected_shared_controls))deallocate(s%selected_shared_controls)
            end if
            call refresh_shared_controls(s)
            s%viewport_feature=0
            s%mesh_pick=.false.
            s%mesh_cell=0
            call on_control_selection(s)
            text=s%widgets%pick_actions(s%feature_pick_mode+1)%text()//' · Ctrl toggles · Shift adds · Esc clears'
            if(command==pick_multipatch_command)text='Click a member to select its multipatch · Ctrl toggles · Shift adds'
            if(command==pick_controls_command.and.s%selected_multipatch>0)then
                text='Drag points · Shift adds · Ctrl toggles · Ctrl+E opens exact values'
                if(.not.s%shared_controls%valid)text='Resolve the multipatch interfaces before editing controls'
            end if
            call s%status%showMessage(text,4000)
        case(isometric_view_command,front_view_command,top_view_command,back_view_command,bottom_view_command,left_view_command,right_view_command)
            i=findloc(camera_commands,command,dim=1)
            call set_viewport_orientation(s,i)
            s%zoom=1._dp
            s%pan_x=0._dp; s%pan_y=0._dp
            call frame_workspace_scene(s)
            call s%status%showMessage(trim(camera_names(i)),2500)
        case(fit_scene_command)
            s%zoom=1._dp; s%pan_x=0._dp; s%pan_y=0._dp
            call frame_workspace_scene(s)
            call s%status%showMessage('Scene fitted',2500)
        case(fit_around_panels_command)
            call fit_around_panels(s)
            call s%status%showMessage('Scene fitted around panels',2500)
        case(fit_selection_command)
            if(s%isolation%count()>0.or.s%mesh_only)then
                focused=primary_workspace_target(s)
                if(allocated(s%selection%items))then
                    changed=fit_workspace_selection(s)
                else
                    changed=fit_workspace_selection(s,[focused])
                end if
                if(changed)then
                    call s%status%showMessage('Selection fitted',2500)
                    call repaint_all(s)
                else
                    call s%status%showMessage('Select a visible scene item to fit',3500)
                end if
                return
            end if
            if(allocated(s%selection%items))then
                changed=s%selection%count()/=1
                if(s%selection%count()==1)then
                    associate(target=>s%selection%items(1))
                        changed=target%category==selection_control.or.(target%category==selection_mesh.and.target%item>0)
                    end associate
                end if
                if(changed)then
                    if(fit_workspace_selection(s))then
                        call s%status%showMessage('Selection fitted',2500)
                        call repaint_all(s)
                    else
                        call s%status%showMessage('Select visible geometry, controls, or mesh cells to fit',3500)
                    end if
                    return
                end if
            end if
            if(s%selected_mesh>0)then
                if(fit_mesh_selection(s,s%selected_mesh))then
                    call s%status%showMessage('Mesh fitted',2500)
                    call repaint_all(s)
                else
                    call s%status%showMessage('Show a current mesh before framing it',3500)
                end if
                return
            end if
            if(s%selected_multipatch>0)then
                associate(definition=>s%document%multipatches%groups(s%selected_multipatch)%definition)
                    if(.not.any([(studio_scene_is_visible(s%document%scene,definition%kind,definition%members(i),catalog=s%document%patches), &
                        i=1,size(definition%members))]))then
                        call s%status%showMessage('Show the multipatch members before framing them',3500)
                        return
                    end if
                    call frame_studio_scene(s%document%scene,s%yaw,s%pitch,kind=definition%kind,members=definition%members, &
                        include_controls=any(s%layers(3:4)))
                end associate
                call fit_framed_scene_in_clear_area(s)
                call s%status%showMessage('Multipatch fitted',2500)
                call repaint_all(s)
                return
            end if
            i=active_patch_id(s%entity_kind,catalog=s%document%patches)
            if(.not.s%live_scene%valid.or..not.studio_scene_is_visible(s%document%scene,s%entity_kind,i,catalog=s%document%patches))then
                call s%status%showMessage('Show the selected object before framing it',3500)
                return
            end if
            call frame_studio_scene(s%document%scene,s%yaw,s%pitch,s%entity_kind,i,include_controls=any(s%layers(3:4)))
            call fit_framed_scene_in_clear_area(s)
            call s%status%showMessage('Selection fitted',2500)
        case(zoom_in_command)
            call zoom_live_viewport(s,1._dp)
        case(zoom_out_command)
            call zoom_live_viewport(s,-1._dp)
        case(restore_layout_command)
            call restore_workspace_layout(s)
        case(select_tool_command,orbit_tool_command)
            call cancel_area_drag(s)
            if(s%transform_mode>0)call on_transform_exit(s)
            s%surface_tool=command-15
            s%area_selection_tool=area_click
            s%mesh_pick=.false.
            s%mesh_cell=0
            if(command==select_tool_command)then
                call s%status%showMessage('Select · Pick geometry or drag a control point',3500)
            else
                call s%status%showMessage('Orbit · Drag to rotate the scene',3500)
            end if
        case(box_selection_command,lasso_selection_command)
            call cancel_area_drag(s)
            if(s%transform_mode>0)call on_transform_exit(s)
            s%surface_tool=1
            s%area_selection_tool=merge(area_box,area_lasso,command==box_selection_command)
            if(command==box_selection_command)then
                call s%status%showMessage('Box · Left to right encloses · Right to left crosses',4500)
            else
                call s%status%showMessage('Lasso · Draw around or across projected objects',4500)
            end if
        case(shaded_layer_command,wireframe_layer_command)
            s%layers(command-shaded_layer_command+1)=.not.s%layers(command-shaded_layer_command+1)
        case(reset_probe_command)
            s%u=.35_dp; s%v=.70_dp; s%w=.45_dp
            call s%status%showMessage('Selection coordinates reset',2500)
        case(control_net_layer_command)
            s%layers(3)=.not.s%layers(3)
        case(parameter_grid_command)
            s%domain_layers(1)=.not.s%domain_layers(1)
        case(basis_support_command)
            s%domain_layers(4)=.not.s%domain_layers(4)
        case(basis_functions_command,basis_derivative_command)
            s%basis_derivative=command==basis_derivative_command
        case(basis_u_command:basis_dw_command,knots_u_command:knots_w_command)
            if(command<=basis_dw_command)then
                i=mod(command-basis_u_command,3)+1
                if(i>s%entity_kind)return
                s%basis_derivative=command>=basis_du_command
            else
                i=command-knots_u_command+1
                if(i>s%entity_kind)return
            end if
            if(s%knot_axis/=i)then
                s%knot_axis=i
                call refresh_knot_table(s)
            end if
        case(insert_knot_command,remove_knot_command)
            call run_knot_editor(s,command==remove_knot_command)
            return
        case(split_patch_command,extract_patch_command,reverse_patch_command)
            call run_patch_tool(s,command)
            return
        case(extract_section_command,extract_face_command)
            call run_patch_section(s,command==extract_face_command)
            return
        case(elevate_degree_command,reduce_degree_command)
            call run_degree_editor(s,command==reduce_degree_command)
            return
        case(add_patch_command)
            i=create_studio_patch(s%window,s%entity_kind,s%document%model_unit,active_length_unit(s),catalog=s%document%patches)
            if(i>0)then
                call rebuild_project_tree(s)
                call select_entity(s,i)
                s%zoom=1._dp; s%pan_x=0._dp; s%pan_y=0._dp
                call frame_workspace_scene(s)
                call s%surface%update()
                call s%status%showMessage('Added '//single_patch_name(i,active_patch_id(i,catalog=s%document%patches),catalog=s%document%patches),5000)
            end if
            return
        case(primitive_command:primitive_command+size(primitive_names))
            i=create_studio_primitive(s%window,command-primitive_command,s%document%model_unit,active_length_unit(s),catalog=s%document%patches)
            if(i>0)call show_created_patch(s,i)
            return
        case(construction_command:construction_command+size(construction_names)-1)
            axis=0
            if(s%selected_multipatch==0.and.s%live_scene%valid)axis=s%entity_kind
            i=create_studio_construction(s%window,command-construction_command+1,axis,active_patch_id(axis,catalog=s%document%patches), &
                s%document%model_unit,active_length_unit(s),merge(1,3,s%document%show_templates),catalog=s%document%patches)
            if(i>0)call show_created_patch(s,i)
            return
        case(toggle_patch_visibility_command)
            call toggle_scene_selection_visibility(s)
            return
        case(project_duplicate_command,project_rename_command)
            call run_project_edit(s,command)
            return
        case(project_remove_command)
            call run_project_delete(s)
            return
        case(project_hide_selection_command)
            call hide_scene_selection(s)
            return
        case(project_isolate_command,project_show_all_command)
            call update_project_visibility(s,show_all=command==project_show_all_command)
            return
        case(scene_exit_isolation_command)
            call exit_scene_isolation(s)
            return
        case(scene_reveal_command)
            call reveal_scene_selection(s)
            return
        case(scene_select_mesh_source_command)
            call select_scene_mesh_source(s)
            return
        case(project_expand_command)
            call s%widgets%project_browser%tree%expandAll()
        case(project_collapse_command)
            call s%widgets%project_browser%tree%collapseAll()
        case(new_multipatch_command)
            call run_multipatch_editor(s,0)
            return
        case(edit_multipatch_command)
            if(s%selected_multipatch>0)call run_multipatch_editor(s,s%selected_multipatch)
            return
        case(insert_assembly_knot_command,elevate_assembly_degree_command)
            i=s%selected_multipatch
            if(i<1.or.i>count_multipatches(s%document%multipatches))return
            if(edit_assembly_refinement(s,i,command==insert_assembly_knot_command))then
                s%live_scene%valid=.false.
                call refresh_inspector(s)
                call refresh_knot_table(s)
                call repaint_all(s)
                call s%status%showMessage('Multipatch refinement applied · Ctrl+Z to undo',5000)
            end if
            return
        case(validate_multipatch_command)
            i=s%selected_multipatch
            if(i==0)return
            call refresh_multipatch(s%document%multipatches,i,.true.,catalog=s%document%patches)
            call refresh_multipatch_card(s)
            call refresh_workspace_status(s)
            call sync_studio_assistant(s)
            call s%status%showMessage(s%document%multipatches%groups(i)%validation%message,4500)
            call s%surface%update()
            return
        case(remove_multipatch_command)
            i=s%selected_multipatch
            if(i<1.or.i>count_multipatches(s%document%multipatches))return
            selected_multipatch_stamp=multipatch_undo_stamp(s%document%multipatches,catalog=s%document%patches)
            result=QMessageBox_question('Dissolve multipatch?', &
                'Member patches stay in the scene. You can undo this change.', &
                buttons=ior(QMessageBox_Yes,QMessageBox_No),defaultButton=QMessageBox_No,parent=s%window)
            if(result/=QMessageBox_Yes.or.i/=s%selected_multipatch)return
            if(selected_multipatch_stamp/=multipatch_undo_stamp(s%document%multipatches,catalog=s%document%patches))return
            call delete_project_assembly(s,i,text)
            if(len(text)>0)call s%status%showMessage(text,5000)
            return
        case(edit_knot_values_command)
            call run_knot_values_editor(s)
            return
        case(edit_interval_command,edit_interval_axis_base+1:edit_interval_axis_base+3)
            i=s%knot_axis
            if(command/=edit_interval_command)i=command-edit_interval_axis_base
            call run_patch_interval(s,i)
            return
        case(apply_interval_command)
            call apply_inline_interval(s)
            return
        case(refine_command)
            call run_refinement(s)
            return
        case(about_command)
            result=QMessageBox_information('About '//studio_product_name, &
                studio_product_name//new_line('a')//'Parametric geometry and structured meshing'//new_line('a')//new_line('a')// &
                'Powered by ForCAD and Qt for Fortran.',parent=s%window)
        case(interaction_guide_command)
            call show_interaction_guide(s%window)
        case(show_knots_command,domain_probe_command,domain_details_command)
            visible=s%panels(5)%isVisible().and.s%domain_edit_knots.and.s%domain_knot_details
            call set_panel_visibility(s,5,.true.)
            s%domain_edit_knots=command/=domain_probe_command
            if(.not.s%domain_edit_knots)s%domain_knot_tool=knot_select_tool
            if(command==domain_details_command)s%domain_knot_details=.not.visible
            call s%domain%setFocus()
        case(domain_select_command:domain_remove_command)
            s%domain_knot_tool=command-domain_select_command
            call s%domain%setFocus()
        case(domain_reset_command)
            call reset_domain_view(s)
            return
        case(domain_zoom_out_command,domain_zoom_in_command)
            call zoom_domain_view(s,merge(1.2_dp,1._dp/1.2_dp,command==domain_zoom_in_command))
            return
        case(show_basis_command)
            call set_panel_visibility(s,6,.true.)
            call s%basis_buttons(merge(2,1,s%basis_derivative))%setFocus()
        case(show_inspector_command)
            if(s%selected_multipatch==0)s%focused_surface=.false.
            if(assistant_shares_inspector_space(s))call hide_studio_assistant(s)
            call set_panel_visibility(s,8,.true.)
            call s%inspector%headers(1)%setFocus()
        case(inspect_patch_command,toggle_patch_inspection_command)
            if(.not.s%live_scene%valid.or.s%selected_mesh>0.or.s%selection%count()>1)return
            visible=command==inspect_patch_command.or..not.(s%panels(5)%isVisible().and.s%panels(6)%isVisible())
            call set_panel_visibility(s,5,visible)
            call set_panel_visibility(s,6,visible)
            if(visible)call s%domain%setFocus()
        case(close_panel_command_base+1:close_panel_command_base+8)
            i=command-close_panel_command_base
            if(i==4.or..not.any(workspace_panels==i))return
            call set_panel_visibility(s,i,.false.)
            call s%status%showMessage('Panel closed · Reopen it from the View ribbon',3500)
        case(toggle_inspector_command)
            if(s%canvas_only.or..not.s%panels(8)%isVisible())then
                call invoke_command(s,show_inspector_command)
            else
                call invoke_command(s,close_panel_command_base+8)
            end if
            return
        case(toggle_panel_command_base+1:toggle_panel_command_base+7)
            i=command-toggle_panel_command_base
            if(i==4.or..not.any(workspace_panels==i))return
            visible=.not.s%panels(i)%isVisible()
            if(visible.and.s%canvas_only)call toggle_canvas(s)
            call set_panel_visibility(s,i,visible)
        case(toggle_parameter_panels_command)
            s%focused_surface=.not.s%focused_surface
            if(s%focused_surface)then
                s%normal_columns=s%columns%sizes(); s%normal_center=s%center_rows%sizes()
                do axis=1,size(workspace_panels)
                    i=workspace_panels(axis)
                    s%saved_panels(i)=s%panels(i)%isVisible()
                end do
                s%saved_panels(8)=s%widgets%inspector_requested
                call set_panel_visibility(s,5,.false.)
                call set_panel_visibility(s,6,.false.)
                call set_panel_visibility(s,8,.false.)
            else
                call leave_viewport_focus(s)
            end if
        case(select_curve_command:select_volume_command)
            i=command-select_curve_command+1
            call select_patch(s,i,active_patch_id(i,catalog=s%document%patches),reveal=.false.)
            return
        case(inspector_section_command_base+1:inspector_section_command_base+2)
            call update_inspector_section(s%inspector,command-inspector_section_command_base)
        case(clean_display_command,control_structure_command)
            i=command-clean_display_command+1
            s%layers=display_layer_presets(:,i)
            s%viewport_style=viewport_edges
            call sync_controls(s)
            call s%status%showMessage('Display preset · '//s%viewport_presets(i)%text(),2500)
            call s%surface%update()
            return
        case(shaded_display_command,shaded_edges_display_command,wireframe_display_command,fast_shaded_display_command)
            s%viewport_style=findloc(viewport_style_commands,command,dim=1)
            if(s%viewport_style/=viewport_wireframe)s%layers(1)=.true.
            call sync_controls(s)
            call s%status%showMessage('3D display · '//trim(viewport_style_names(s%viewport_style)),2500)
            call s%surface%update()
            return
        case(viewport_layer_base:viewport_layer_base+viewport_layer_count-1)
            i=command-viewport_layer_base+1
            s%layers(i)=.not.s%layers(i)
            if(i==4)call refresh_shared_controls(s)
            call sync_controls(s)
            call s%surface%update()
            return
        case(interactive_quality_command)
            call set_interactive_quality(s,.not.s%interactive_quality)
            call s%widgets%ribbon%interactive_quality_action%setChecked(s%interactive_quality)
            if(s%interactive_quality)then
                call s%status%showMessage('Interactive Quality on',2500)
            else
                call s%status%showMessage('Interactive Quality off',2500)
            end if
            return
        case(world_grid_command)
            s%world_grid=.not.s%world_grid
            call sync_controls(s)
            call s%surface%update()
            return
        case(basis_legend_command,basis_probe_command)
            if(command==basis_legend_command)s%basis_legend=.not.s%basis_legend
            if(command==basis_probe_command)s%basis_probe_guide=.not.s%basis_probe_guide
            call sync_controls(s)
            call s%basis%update()
            return
        case(find_command)
            if(s%home_active)return
            focused=primary_workspace_target(s)
            result=find_studio_command(s%window,s%entity_kind,s%widgets%edit_controls_action%isEnabled(), &
                document_undo_stamp(s%document)>0,can_remove_knot(s,any_direction=.true.), &
                can_reduce_degree(s),has_multipatch=s%selected_multipatch>0,can_restructure_patch=can_restructure_patch(s), &
                can_redo=document_can_redo(s%document), &
                can_reveal=focused%category>0.and.(.not.allocated(s%selection%items).or.s%selection%count()>0), &
                scene_isolated=s%isolation%count()>0,has_mesh=s%selected_mesh>0, &
                has_shape_target=s%widgets%transform_action%isEnabled(),has_project=s%document%document_started, &
                mesh_available=[s%widgets%ribbon%mesh_generate_action%isEnabled(),s%widgets%ribbon%mesh_edges_action%isEnabled(), &
                s%widgets%ribbon%mesh_export_action%isEnabled(),s%widgets%ribbon%mesh_pick_action%isEnabled(),s%widgets%ribbon%mesh_only_action%isEnabled(), &
                s%widgets%ribbon%physical_knot_action%isEnabled(),s%widgets%ribbon%mesh_shaded_action%isEnabled()])
            if(result>0)call invoke_command(s,result)
            return
        case(toggle_assistant_command)
            if(s%canvas_only.or..not.studio_assistant_is_visible(s))then
                call invoke_command(s,show_assistant_command)
            else
                call hide_studio_assistant(s)
                call s%surface%setFocus()
            end if
            return
        case(show_assistant_command,diagnose_with_ai_command)
            call show_studio_assistant(s,sync_assistant_layout)
            if(command==diagnose_with_ai_command)call ask_studio_diagnostic(s)
            call sync_assistant_layout(s)
            return
        case(edit_controls_command)
            call edit_workspace_controls(s)
            return
        case(undo_patch_command)
            call undo_project_edit(s)
            return
        case(redo_patch_command)
            call redo_project_edit(s)
            return
        case(close_window_command)
            call s%window%close(); return
        case(minimize_window_command)
            call s%window%showMinimized(); return
        case(fullscreen_command)
            call toggle_fullscreen(s)
            return
        case(maximize_window_command)
            if(s%window%boolProperty('fullScreen'))then
                call toggle_fullscreen(s)
            else
                if(s%window%isMaximized())then
                    call s%window%showNormal()
                else
                    call s%window%showMaximized()
                end if
                call sync_window_controls(s)
            end if
            return
        case default
            return
        end select
        call sync_assistant_layout(s)
        call sync_controls(s)
        call refresh_inspector(s); call repaint_all(s)
    end procedure invoke_command

    subroutine toggle_fullscreen(s)
        type(studio_t),intent(inout) :: s

        if(s%window%boolProperty('fullScreen'))then
            if(s%widgets%restore_maximized)then
                call s%window%showMaximized()
            else
                call s%window%showNormal()
            end if
        else
            s%widgets%restore_maximized=s%window%isMaximized()
            call s%window%showFullScreen()
        end if
        call sync_window_controls(s)
    end subroutine toggle_fullscreen

    subroutine measure_selection(s)
        use studio_workspace_selection, only: selection_measurement_targets
        use studio_measurements, only: measurement_t, measure_document_selection, measurement_warning
        use studio_viewport_scene, only: viewport_feature_text
        type(studio_t), intent(inout) :: s
        type(selection_target_t), allocatable :: targets(:)
        type(measurement_t) :: measures(3)
        character(:), allocatable :: title,feature

        targets=selection_measurement_targets(s)
        title='Selection'
        if(size(targets)==1)then
            title=single_patch_name(targets(1)%rank,targets(1)%id,catalog=s%document%patches)
            if(targets(1)%feature>0)then
                feature=viewport_feature_text(s)
                if(len(feature)>0)title=feature//' · '//title
            end if
        end if
        call measure_document_selection(s%document,targets,measures)
        if(measures(3)%reversed)then
            call record_studio_diagnostic(s,measurement_warning(measures(3)))
            call refresh_studio_diagnostics(s)
        end if
        call show_geometry_measure(s%window,s%document,measures,title)
    end subroutine measure_selection

end submodule studio_ui_commands
