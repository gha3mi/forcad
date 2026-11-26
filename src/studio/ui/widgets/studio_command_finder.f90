! Native command finder for presentation actions already provided by ForCAD Studio.
module studio_command_finder
    use studio_command_presentation, only: command_visual
    use studio_theme, only: studio_frame_edge, studio_active_edge, studio_selected_fill, studio_style_primary
    use qt, only: QDialog, QDialog_Accepted, QDialogButtonBox, QDialogButtonBox_AcceptRole, &
        QDialogButtonBox_Cancel, QHBoxLayout, QKeyEvent, QLabel, QLineEdit, QListWidget, QPushButton, &
        QStackedWidget, Qt_AlignCenter, Qt_Key_Down, Qt_Key_Enter, Qt_Key_Return, Qt_Key_Up, Qt_NoModifier, &
        QtfEventFilter, QVBoxLayout, QWidget
    use studio_commands, only: &
        toggle_assistant_command, fit_scene_command, fit_selection_command, refine_command, &
        translate_geometry_command, affine_geometry_command, primitive_command, extrude_command, &
        project_duplicate_command, project_rename_command, project_isolate_command, project_show_all_command, &
        split_patch_command, extract_patch_command, reverse_patch_command, snap_off_command, &
        extract_section_command,extract_face_command, &
        snap_vertex_command, snap_closest_command, generate_mesh_command, mesh_edges_command, &
        export_mesh_command, mesh_pick_command, mesh_only_command, physical_knots_command, &
        mesh_shaded_command, home_command, new_project_command, open_project_command, &
        save_project_command, save_as_command, import_project_command, units_command, edit_controls_command, &
        undo_patch_command, redo_patch_command, insert_knot_command, remove_knot_command, elevate_degree_command, &
        reduce_degree_command, add_patch_command, edit_knot_values_command, edit_interval_command, &
        new_multipatch_command, edit_multipatch_command, validate_multipatch_command, remove_multipatch_command, &
        insert_assembly_knot_command,elevate_assembly_degree_command, &
        export_iges_command, export_fortran_command, color_command, translucent_command, transform_command, &
        canvas_command, pick_mode_base, pick_controls_command, pick_multipatch_command, clear_selection_command, export_vtk_command, back_view_command, &
        bottom_view_command, left_view_command, right_view_command, toggle_patch_visibility_command, &
        reset_controls_command, screenshot_command, import_iges_command, isometric_view_command, &
        zoom_in_command, zoom_out_command, restore_layout_command, front_view_command, interactive_quality_command, &
        top_view_command, select_tool_command, orbit_tool_command, reset_probe_command, &
        parameter_grid_command, basis_support_command, basis_u_command, basis_v_command, &
        basis_w_command, basis_du_command, basis_dv_command, basis_dw_command, &
        knots_u_command, knots_v_command, knots_w_command, interaction_guide_command, about_command, fullscreen_command, &
        show_knots_command, show_basis_command, show_inspector_command, toggle_patch_inspection_command, &
        select_curve_command, select_surface_command, select_volume_command, clean_display_command, &
        control_structure_command, wireframe_display_command, fast_shaded_display_command, &
        shaded_edges_display_command, shaded_display_command, rotate_geometry_command, scale_geometry_command, &
        revolve_command, sweep_command, loft_command, measure_geometry_command, project_folder_command, &
        open_mesh_dialog_command, fit_around_panels_command, domain_details_command, &
        scene_reveal_command, scene_exit_isolation_command, scene_select_mesh_source_command, &
        viewport_layer_base, parameter_lines_layer, world_grid_command
    implicit none(type, external)
    private
    public :: find_studio_command
    type :: command_entry
        integer :: id
        character(len=40) :: label
        character(len=64) :: keywords
        character(len=14) :: shortcut = ''
        integer :: min_dimension = 1
    end type command_entry

    type(command_entry), parameter :: commands(*) = [ &
        command_entry(select_curve_command,'Recall last curve','select curve geometry scene patch','Ctrl+1'), &
        command_entry(select_surface_command,'Recall last surface','select surface geometry scene patch','Ctrl+2'), &
        command_entry(select_volume_command,'Recall last volume','select volume geometry scene patch','Ctrl+3'), &
        command_entry(isometric_view_command,'Isometric view','camera reset home','R'), &
        command_entry(fit_scene_command,'Fit scene','camera frame visible geometry','F'), &
        command_entry(fit_selection_command,'Fit selection','camera frame all selected geometry','Shift+F'), &
        command_entry(fit_around_panels_command,'Fit around panels','camera fit clear work area'), &
        command_entry(front_view_command,'Front view','camera orientation'), &
        command_entry(back_view_command,'Back view','camera orientation rear'), &
        command_entry(left_view_command,'Left view','camera orientation side'), &
        command_entry(right_view_command,'Right view','camera orientation side'), &
        command_entry(top_view_command,'Top view','camera orientation'), &
        command_entry(bottom_view_command,'Bottom view','camera orientation underside'), &
        command_entry(zoom_in_command,'Zoom in','camera magnify','Ctrl++'), &
        command_entry(zoom_out_command,'Zoom out','camera reduce','Ctrl+-'), &
        command_entry(toggle_patch_inspection_command,'Parameter panels','show hide domain basis plot inspection','Shift+F11'), &
        command_entry(restore_layout_command,'Restore all panels','window layout workspace'), &
        command_entry(toggle_assistant_command,'Toggle Assistant','open close show hide agent chat OpenCode Codex Claude Ollama','Ctrl+Shift+A'), &
        command_entry(toggle_patch_visibility_command,'Show or hide selection','scene object visibility eye','Ctrl+Shift+H'), &
        command_entry(import_project_command,'Import FCAD project…','add merge another project geometry meshes'), &
        command_entry(import_iges_command,'Import IGES…','open read igs curves surfaces','Ctrl+I'), &
        command_entry(export_iges_command,'Export geometry (IGES)…','save geometry igs exchange'), &
        command_entry(export_vtk_command,'Export geometry (VTK)…','save sampled geometry visualization vtk exchange'), &
        command_entry(export_fortran_command,'Export Fortran…','save project rebuild source program f90 exchange'), &
        command_entry(open_mesh_dialog_command,'Generate mesh…','remesh regenerate structured sizing cell counts target size'), &
        command_entry(mesh_edges_command,'Show or hide selected meshes','visibility generated mesh cells viewport'), &
        command_entry(export_mesh_command,'Export mesh (VTK)…','save generated structured cells connectivity'), &
        command_entry(mesh_pick_command,'Pick mesh cells','select line quadrilateral hexahedron cell viewport'), &
        command_entry(mesh_only_command,'Meshes only','hide CAD geometry display mesh viewport filter'), &
        command_entry(mesh_shaded_command,'Filled mesh faces','shade solid faces wireframe mesh display appearance'), &
        command_entry(physical_knots_command,'Knot lines','geometry viewport physical knot boundaries spans grid'), &
        command_entry(viewport_layer_base+1,'Face grid','display wire overlay surface volume parameter lattice'), &
        command_entry(viewport_layer_base+parameter_lines_layer-1,'Parameter lines','display guides isocurves sample selection'), &
        command_entry(world_grid_command,'Ground grid','display viewport work plane'), &
        command_entry(color_command,'Object color…','appearance shading material rgb'), &
        command_entry(translucent_command,'Translucent','material color opaque transparent'), &
        command_entry(transform_command,'Exact transform…','numeric move translate rotation scale matrix affine','Ctrl+T'), &
        command_entry(translate_geometry_command,'Move','move shift offset'), &
        command_entry(rotate_geometry_command,'Rotate','orientation angles pivot'), &
        command_entry(scale_geometry_command,'Scale','resize factors pivot'), &
        command_entry(affine_geometry_command,'Affine matrix…','linear matrix offset pivot'), &
        command_entry(primitive_command,'Create geometry','primitive point line circle arc rectangle box cylinder tube'), &
        command_entry(extrude_command,'Extrude','construction profile direction displacement surface volume'), &
        command_entry(revolve_command,'Revolve','construction profile axis angle surface volume'), &
        command_entry(sweep_command,'Sweep','construction profile path spine translation'), &
        command_entry(loft_command,'Loft','construction interpolate ordered sections surface volume'), &
        command_entry(project_duplicate_command,'Duplicate selected patch','copy object project'), &
        command_entry(project_rename_command,'Rename selected patch','name label object project'), &
        command_entry(project_isolate_command,'Isolate selection','hide other objects visibility scene'), &
        command_entry(scene_exit_isolation_command,'Exit isolation','restore previous scene visibility unisolate'), &
        command_entry(project_show_all_command,'Show all objects','reveal unhide geometry and meshes visibility scene'), &
        command_entry(scene_reveal_command,'Reveal in Scene','locate browser tree selection'), &
        command_entry(scene_select_mesh_source_command,'Select mesh source','locate original geometry in scene browser'), &
        command_entry(snap_off_command,'Turn snapping off','disable control placement snapping'), &
        command_entry(snap_vertex_command,'Snap to vertices','snap corner point control placement'), &
        command_entry(snap_closest_command,'Snap to closest point','snap geometry curve surface volume control placement'), &
        command_entry(split_patch_command,'Split at parameter…','divide partition parameter modeling'), &
        command_entry(extract_patch_command,'Extract range…','restrict trim bounded patch interval modeling'), &
        command_entry(extract_section_command,'Extract section…','exact isoparametric surface curve slice','',2), &
        command_entry(extract_face_command,'Extract face…','exact volume boundary surface slice','',3), &
        command_entry(reverse_patch_command,'Reverse parameter direction','flip orientation active axis U V W modeling'), &
        command_entry(canvas_command,'Toggle canvas only','full viewport ribbon hide panels','Ctrl+Shift+F'), &
        command_entry(fullscreen_command,'Toggle fullscreen','full screen window display','F11'), &
        command_entry(pick_mode_base,'Automatic picking','select body face edge vertex'), &
        command_entry(pick_mode_base+1,'Pick faces','select face boundary'), &
        command_entry(pick_mode_base+2,'Pick edges','select edges boundary'), &
        command_entry(pick_mode_base+3,'Pick vertices','select corners points'), &
        command_entry(pick_controls_command,'Pick control points','select spline controls vertices'), &
        command_entry(pick_multipatch_command,'Pick multipatches','select assembly groups'), &
        command_entry(clear_selection_command,'Clear selection','deselect none'), &
        command_entry(select_tool_command,'Select tool','pick probe cursor'), &
        command_entry(orbit_tool_command,'Orbit tool','rotate navigate camera'), &
        command_entry(measure_geometry_command,'Measure geometry…','inspect length area volume dimensions'), &
        command_entry(show_inspector_command,'Geometry Inspector','selection probe panel'), &
        command_entry(home_command,'Home','welcome start samples'), &
        command_entry(new_project_command,'New project','blank empty document','Ctrl+N'), &
        command_entry(open_project_command,'Open project…','load file studio fcad','Ctrl+O'), &
        command_entry(save_project_command,'Save project','write file studio fcad','Ctrl+S'), &
        command_entry(save_as_command,'Save project as…','write copy file studio fcad','Ctrl+Shift+S'), &
        command_entry(project_folder_command,'Working folder…','project directory files attachments assistant'), &
        command_entry(units_command,'Project units…','millimeters meters inches display scale'), &
        command_entry(add_patch_command,'Copy a starting patch…','new create duplicate existing patch'), &
        command_entry(new_multipatch_command,'New multipatch…','create group assembly connect multi-patch interfaces', &
                      'Ctrl+Shift+M'), &
        command_entry(edit_multipatch_command,'Edit assembly…','manage assembly group membership connection interface'), &
        command_entry(insert_assembly_knot_command,'Insert multipatch knot…','refine assembly connected knot split span'), &
        command_entry(elevate_assembly_degree_command,'Elevate multipatch degree…','refine assembly connected increase degree'), &
        command_entry(validate_multipatch_command,'Check interfaces','verify connection interface continuity diagnostics'), &
        command_entry(remove_multipatch_command,'Dissolve multipatch','delete assembly group'), &
        command_entry(edit_controls_command,'Control points & weights…','edit coordinates xyz weight','Ctrl+E'), &
        command_entry(undo_patch_command,'Undo modeling edit','restore geometry point knot degree','Ctrl+Z'), &
        command_entry(redo_patch_command,'Redo modeling edit','restore last undone geometry change','Ctrl+Shift+Z'), &
        command_entry(show_knots_command,'Edit knots in domain','knot vectors topology handle editor'), &
        command_entry(domain_details_command,'Show knot groups','parameter knot editing domain table'), &
        command_entry(insert_knot_command,'Insert knot','add knot refine multiplicity'), &
        command_entry(remove_knot_command,'Remove knot','delete knot multiplicity'), &
        command_entry(edit_knot_values_command,'Edit knot values…','move modify knot vector'), &
        command_entry(edit_interval_command,'Set parameter range…','domain bounds range reparameterize'), &
        command_entry(elevate_degree_command,'Elevate degree','increase degree refine order'), &
        command_entry(refine_command,'Refine patch…','split spans subdivide h p k degree elevation'), &
        command_entry(reduce_degree_command,'Reduce degree','lower degree tolerance order'), &
        command_entry(show_basis_command,'Show basis plot','basis inspector derivative panel'), &
        command_entry(reset_probe_command,'Reset probe','parameter coordinates selection'), &
        command_entry(parameter_grid_command,'Toggle parameter grid','domain overlay'), &
        command_entry(basis_support_command,'Toggle basis support','domain overlay'), &
        command_entry(interactive_quality_command,'Interactive Quality','viewport display antialiasing zoom performance'), &
        command_entry(clean_display_command,'Clean preview','display preset shaded'), &
        command_entry(control_structure_command,'Control structure','display preset net lattice polygon'), &
        command_entry(shaded_display_command,'Shaded','display geometry faces'), &
        command_entry(shaded_edges_display_command,'Shaded with Edges','display geometry faces full quality'), &
        command_entry(fast_shaded_display_command,'Fast Shaded','viewport display performance zoom'), &
        command_entry(wireframe_display_command,'Wireframe','display geometry mesh edges'), &
        command_entry(basis_u_command,'Basis in U','function chart plot'), &
        command_entry(basis_v_command,'Basis in V','function chart plot',min_dimension=2), &
        command_entry(basis_w_command,'Basis in W','function chart plot',min_dimension=3), &
        command_entry(basis_du_command,'Derivative in U','basis chart plot'), &
        command_entry(basis_dv_command,'Derivative in V','basis chart plot',min_dimension=2), &
        command_entry(basis_dw_command,'Derivative in W','basis chart plot',min_dimension=3), &
        command_entry(knots_u_command,'Knot axis U','parameter space knot vector'), &
        command_entry(knots_v_command,'Knot axis V','parameter space knot vector',min_dimension=2), &
        command_entry(knots_w_command,'Knot axis W','parameter space knot vector',min_dimension=3), &
        command_entry(screenshot_command,'Capture screenshot','capture export image png save','Ctrl+Alt+S'), &
        command_entry(reset_controls_command,'Reset workspace controls','restore camera defaults','Ctrl+R'), &
        command_entry(about_command,'About ForCAD Studio','application information'), &
        command_entry(interaction_guide_command,'Interaction guide','help shortcuts keyboard','F1') ]

    type :: finder_state
        type(QDialog) :: dialog
        type(QListWidget) :: results
        type(QStackedWidget) :: views
        type(QPushButton) :: run
        integer :: entity_kind
        logical :: can_edit=.false., can_undo=.false., can_redo=.false., can_remove_knot=.false., can_reduce_degree=.false.
        logical :: has_multipatch=.false., can_restructure_patch=.false.
        logical :: can_reveal=.false., scene_isolated=.false., has_mesh=.false., has_shape_target=.false., has_project=.false.
        logical :: mesh_available(generate_mesh_command:mesh_shaded_command)=.true.
        integer :: result_ids(size(commands)) = 0
        integer :: selected_id = 0
    end type finder_state
contains
    integer function find_studio_command(parent,entity_kind,can_edit,can_undo,can_remove_knot,can_reduce_degree, &
                                          has_multipatch,can_restructure_patch,mesh_available,can_redo, &
                                          can_reveal,scene_isolated,has_mesh,has_shape_target,has_project) result(command)
        class(QWidget), intent(in) :: parent
        integer, intent(in) :: entity_kind
        logical,optional,intent(in) :: can_edit,can_undo,can_remove_knot,can_reduce_degree,can_redo
        logical,optional,intent(in) :: has_multipatch,can_restructure_patch
        logical,optional,intent(in) :: can_reveal,scene_isolated,has_mesh,has_shape_target,has_project
        logical,optional,intent(in) :: mesh_available(generate_mesh_command:mesh_shaded_command)
        type(finder_state), target :: state
        type(QVBoxLayout) :: layout
        type(QHBoxLayout) :: footer
        type(QLineEdit) :: query
        type(QLabel) :: empty, keys, title
        type(QDialogButtonBox) :: buttons
        type(QtfEventFilter) :: query_keys, result_keys

        state%entity_kind = entity_kind
        if(present(can_edit))state%can_edit=can_edit
        if(present(can_undo))state%can_undo=can_undo
        if(present(can_redo))state%can_redo=can_redo
        if(present(can_remove_knot))state%can_remove_knot=can_remove_knot
        if(present(can_reduce_degree))state%can_reduce_degree=can_reduce_degree
        if(present(has_multipatch))state%has_multipatch=has_multipatch
        state%can_restructure_patch=state%can_edit.and..not.state%has_multipatch
        if(present(can_restructure_patch))state%can_restructure_patch=can_restructure_patch
        if(present(mesh_available))state%mesh_available=mesh_available
        if(present(can_reveal))state%can_reveal=can_reveal
        if(present(scene_isolated))state%scene_isolated=scene_isolated
        if(present(has_mesh))state%has_mesh=has_mesh
        if(present(has_shape_target))state%has_shape_target=has_shape_target
        if(present(has_project))state%has_project=has_project
        state%dialog = QDialog(parent)
        call state%dialog%setObjectName('studioCommandFinder')
        call state%dialog%setWindowTitle('Find command')
        call state%dialog%setMinimumSize(480,330)
        call state%dialog%resize(530,410)
        call state%dialog%setStyleSheet( &
            'QLineEdit#studioCommandQuery { min-height:32px; padding:3px 10px; }'// &
            'QLineEdit#studioCommandQuery:focus { border-color:'//studio_active_edge//'; }'// &
            'QListWidget#studioCommandResults, QLabel#studioCommandEmpty {'// &
            'background:#151f29; border:1px solid '//studio_frame_edge//'; padding:5px; }'// &
            'QListWidget#studioCommandResults::item { padding:7px 9px; }'// &
            'QListWidget#studioCommandResults::item:hover:!selected { background:#243541; }'// &
            'QListWidget#studioCommandResults::item:selected { background:'//studio_selected_fill//'; color:#f0f8fa; }')
        layout = QVBoxLayout(state%dialog)
        call layout%setContentsMargins(18,16,18,16)
        call layout%setSpacing(10)

        title=QLabel('Find command',state%dialog)
        call title%setObjectName('studioDialogTitle')
        call layout%addWidget(title)
        query = QLineEdit(parent=state%dialog)
        call query%setObjectName('studioCommandQuery')
        call query%setAccessibleName('Search commands')
        call query%setPlaceholderText('Search by command name or shortcut')
        call query%setClearButtonEnabled(.true.)
        call layout%addWidget(query)

        state%views = QStackedWidget(state%dialog)
        state%results = QListWidget(state%views)
        call state%results%setObjectName('studioCommandResults')
        call state%results%setAccessibleName('Matching commands')
        call state%results%setMouseTracking(.true.)
        empty = QLabel('No available commands match this search.',state%views)
        call empty%setObjectName('studioCommandEmpty')
        call empty%setAlignment(Qt_AlignCenter)
        call state%views%addWidget(state%results)
        call state%views%addWidget(empty)
        call layout%addWidget(state%views,1)

        footer = QHBoxLayout()
        keys = QLabel('↑ ↓ Select  ·  Enter Run  ·  Esc Close')
        call keys%setProperty('dialogRole','detail')
        call footer%addWidget(keys)
        call footer%addStretch()
        buttons = QDialogButtonBox(QDialogButtonBox_Cancel,state%dialog)
        state%run = buttons%addButton('Run',QDialogButtonBox_AcceptRole)
        call state%run%setObjectName('studioCommandRun')
        call studio_style_primary(state%run)
        call state%run%setDefault(.true.)
        call footer%addWidget(buttons)
        call layout%addLayout(footer)

        call query%textChanged%connect(filter_commands,state)
        call state%results%currentRowChanged%connect(update_run_button,state)
        call state%results%itemDoubleClicked%connect(pick_row,state)
        call buttons%accepted%connect(pick_current,state)
        call buttons%rejected%connect(cancel,state)
        query_keys = QtfEventFilter(query)
        result_keys = QtfEventFilter(state%results)
        call query_keys%keyPressEvent%connect(on_query_key,state)
        call result_keys%keyPressEvent%connect(on_result_key,state)
        call filter_commands('',state)
        call query%setFocus()

        command = 0
        if (state%dialog%exec()==QDialog_Accepted) command=state%selected_id
        call state%dialog%delete()
    end function find_studio_command

    subroutine filter_commands(input,data)
        character(len=*), intent(in) :: input
        class(*), intent(inout), optional :: data
        character(:), allocatable :: query, label, shortcut, haystack, row_text
        integer :: i, count, rank, position
        integer :: ranks(size(commands))
        select type(state=>data)
        type is(finder_state)
            query=lowercase(trim(adjustl(input)))
            call state%results%clear()
            state%result_ids=0
            count=0
            do i=1,size(commands)
                if (state%entity_kind<commands(i)%min_dimension) cycle
                select case(commands(i)%id)
                case(insert_assembly_knot_command,elevate_assembly_degree_command)
                    if(.not.(state%has_multipatch.and.state%can_edit))cycle
                case(edit_multipatch_command,validate_multipatch_command,remove_multipatch_command)
                    if(.not.(state%has_multipatch.and.state%has_shape_target))cycle
                case(22,25,27,30:35,40:42,show_knots_command,show_basis_command,toggle_patch_inspection_command)
                    if(.not.state%can_edit.or.state%has_multipatch)cycle
                case(open_mesh_dialog_command)
                    if(.not.state%mesh_available(generate_mesh_command))cycle
                case(scene_reveal_command,project_isolate_command)
                    if(.not.state%can_reveal)cycle
                case(scene_exit_isolation_command)
                    if(.not.state%scene_isolated)cycle
                case(scene_select_mesh_source_command)
                    if(.not.state%has_mesh)cycle
                case(domain_details_command)
                    if(.not.state%can_edit.or.state%has_multipatch)cycle
                case(generate_mesh_command:mesh_shaded_command)
                    if(.not.state%mesh_available(commands(i)%id))cycle
                case(transform_command,translate_geometry_command,rotate_geometry_command,scale_geometry_command, &
                    affine_geometry_command,export_vtk_command)
                    if(.not.state%has_shape_target)cycle
                case(export_iges_command)
                    if(.not.state%can_edit.or.state%entity_kind==3.or.state%has_multipatch)cycle
                case(export_fortran_command)
                    if(.not.state%has_project)cycle
                case(edit_controls_command)
                    if(.not.state%can_edit)cycle
                case(insert_knot_command,elevate_degree_command,edit_interval_command,254,290,291)
                    if(.not.state%can_edit.or.state%has_multipatch)cycle
                case(split_patch_command,reverse_patch_command)
                    if(.not.state%can_restructure_patch)cycle
                case(extract_patch_command,extract_section_command,extract_face_command)
                    if(.not.state%can_edit.or.state%has_multipatch)cycle
                case(remove_knot_command,edit_knot_values_command)
                    if(.not.state%can_remove_knot.or.state%has_multipatch)cycle
                case(reduce_degree_command)
                    if(.not.state%can_reduce_degree.or.state%has_multipatch)cycle
                case(undo_patch_command)
                    if(.not.state%can_undo)cycle
                case(redo_patch_command)
                    if(.not.state%can_redo)cycle
                case default
                    continue
                end select
                associate(visual=>command_visual(commands(i)%id,trim(commands(i)%label),0))
                    row_text=visual%caption
                end associate
                label=lowercase(row_text)
                shortcut=lowercase(trim(commands(i)%shortcut))
                if (len(shortcut)>0) row_text=row_text//'  ·  '//trim(commands(i)%shortcut)
                haystack=label//' '//shortcut//' '//lowercase(trim(commands(i)%keywords)//' '//trim(commands(i)%label))
                if (.not.matches(query,haystack)) cycle
                rank=0
                if(len(query)>0)then
                    if(matches(query,label))rank=1
                    if(index(label,query)==1)rank=2
                    if(query==label.or.query==shortcut)rank=3
                end if
                position=count+1
                do while(position>1)
                    if(ranks(position-1)>=rank)exit
                    ranks(position)=ranks(position-1)
                    state%result_ids(position)=state%result_ids(position-1)
                    position=position-1
                end do
                count=count+1
                ranks(position)=rank
                state%result_ids(position)=commands(i)%id
                call state%results%insertItem(position-1,row_text)
            end do
            call state%views%setCurrentIndex(merge(1,0,count==0))
            if (count>0.and.len_trim(adjustl(query))>0) call state%results%setCurrentRow(0)
            call update_run_button(state%results%currentRow(),state)
        end select
    end subroutine filter_commands

    pure logical function matches(query,haystack)
        character(len=*), intent(in) :: query, haystack
        character(:), allocatable :: term
        integer :: first, last, position
        first=1
        do while (first<=len(query))
            if (index(' '//achar(9),query(first:first))>0) then
                first=first+1
                cycle
            end if
            last=first
            do while (last<=len(query))
                if (index(' '//achar(9),query(last:last))>0) exit
                last=last+1
            end do
            term=query(first:last-1)
            if (len(term)==1.and.index('uvw',term)>0.and.len_trim(adjustl(query))>1) then
                ! Axis initials are words in a multiword query. W must not match "show" or "view".
                position=index(' '//haystack//' ',' '//term//' ')
            else
                position=index(haystack,term)
            end if
            if (position==0) then
                matches=.false.
                return
            end if
            first=last
        end do
        matches=.true.
    end function matches

    pure function lowercase(input) result(output)
        character(len=*), intent(in) :: input
        character(len=len(input)) :: output
        integer :: i, letter
        output=input
        do i=1,len(input)
            letter=iachar(input(i:i))
            if (letter>=iachar('A').and.letter<=iachar('Z')) output(i:i)=achar(letter+32)
        end do
    end function lowercase

    subroutine update_run_button(row,data)
        integer, intent(in) :: row
        class(*), intent(inout), optional :: data
        select type(state=>data)
        type is(finder_state)
            call state%run%setEnabled(row>=0.and.row<state%results%count())
        end select
    end subroutine update_run_button

    subroutine pick_row(row,data)
        integer, intent(in) :: row
        class(*), intent(inout), optional :: data
        select type(state=>data)
        type is(finder_state)
            if (row<0.or.row>=state%results%count()) return
            state%selected_id=state%result_ids(row+1)
            if (state%selected_id>0) call state%dialog%accept()
        end select
    end subroutine pick_row

    subroutine pick_current(data)
        class(*), intent(inout), optional :: data
        select type(state=>data)
        type is(finder_state)
            call pick_row(state%results%currentRow(),state)
        end select
    end subroutine pick_current

    subroutine cancel(data)
        class(*), intent(inout), optional :: data
        select type(state=>data)
        type is(finder_state)
            call state%dialog%reject()
        end select
    end subroutine cancel

    subroutine on_query_key(event,data)
        type(QKeyEvent), intent(in) :: event
        class(*), intent(inout), optional :: data
        integer :: row, count
        select type(state=>data)
        type is(finder_state)
            if (event%key()==Qt_Key_Return.or.event%key()==Qt_Key_Enter) then
                call on_result_key(event,state)
                return
            end if
            if (event%modifiers()/=Qt_NoModifier) return
            count=state%results%count()
            if (count==0) return
            row=state%results%currentRow()
            select case(event%key())
            case(Qt_Key_Up)
                call state%results%setCurrentRow(max(0,row-1))
            case(Qt_Key_Down)
                call state%results%setCurrentRow(min(count-1,row+1))
            case default
                return
            end select
            call event%accept()
        end select
    end subroutine on_query_key

    subroutine on_result_key(event,data)
        type(QKeyEvent), intent(in) :: event
        class(*), intent(inout), optional :: data
        if (event%key()/=Qt_Key_Return.and.event%key()/=Qt_Key_Enter) return
        call event%accept()
        call pick_current(data)
    end subroutine on_result_key
end module studio_command_finder
