submodule(studio_ui) studio_ui_inspection
    use qt, only: QHeaderView, Qt_NoButton
    use studio_viewport_scene, only: selected_controls_available
    use studio_selection_status, only: probe_status_text,knot_status_text
    implicit none(type, external)
contains
    function fixed(value,places) result(text)
        real(dp),intent(in)::value
        integer,intent(in),optional::places
        character(len=:),allocatable::text
        character(len=32)::buffer,fmt
        integer::n
        n=4; if(present(places))n=places
        write(fmt,'("(F20.",I0,")")')n
        write(buffer,fmt)value
        text=trim(adjustl(buffer))
    end function fixed

    function html_row(key,value) result(text)
        character(len=*),intent(in)::key,value
        character(len=:),allocatable::text
        text='<tr><td style="border-right:1px solid #364451;border-bottom:1px solid #33414d;'// &
            'color:#cce3f2;">'//key//'</td><td style="border-bottom:1px solid #33414d;">'//value//'</td></tr>'
    end function html_row

    module procedure compact_scalar
        character(len=34)::buffer
        integer::last
        if(value==0._dp)then
            text='0'
            return
        end if
        write(buffer,'(g0.7)')value
        text=trim(buffer)
        if(index(text,'E')/=0.or.index(text,'.')==0)return
        last=len(text)
        do while(last>1)
            if(text(last:last)/='0')exit
            last=last-1
        end do
        if(text(last:last)=='.')last=last-1
        text=text(:last)
    end procedure compact_scalar

    function coordinates_text(values,precise) result(text)
        real(dp),intent(in)::values(3)
        logical,intent(in)::precise
        character(len=:),allocatable::text
        character(len=26)::buffer
        real(dp)::value
        integer::i
        text='('
        do i=1,3
            if(i>1)text=text//', '
            value=values(i)
            if(precise)then
                write(buffer,'(es26.17e3)')value
                text=text//trim(adjustl(buffer))
            else
                if(abs(value)<4._dp*epsilon(1._dp))value=0._dp
                text=text//compact_scalar(value)
            end if
        end do
        text=text//')'
    end function coordinates_text

    module procedure select_entity

        call select_patch(s,kind,variant,reveal=.true.)
    end procedure select_entity

    module procedure select_patch
        logical::was_updating,new_target

        if(kind<1.or.kind>3)return
        if(present(variant))then
            if(.not.document_patch(s,kind,variant))return
            call set_active_patch(kind,variant,catalog=s%document%patches)
            if(reveal)then
                call studio_scene_set_visible(s%document%scene,kind,variant,.true.,catalog=s%document%patches)
                call sync_scene_visibility_icons(s%widgets%project_browser,s,kind,variant)
            end if
        end if
        call forget_workspace_selection(s)
        new_target=s%selected_mesh>0.or.s%selected_multipatch>0.or.kind/=s%entity_kind.or. &
            active_patch_id(kind,catalog=s%document%patches)/=s%live_scene%variant
        call clear_mesh_selection(s)
        s%selected_multipatch=0
        s%viewport_feature=0
        s%entity_kind=kind
        if(s%knot_axis>kind)s%knot_axis=1
        was_updating=s%updating
        s%updating=.true.
        call select_project_patch(s,kind)
        s%updating=was_updating
        call refresh_inspector(s)
        call refresh_knot_table(s)
        if(new_target)call reset_transform_values(s)
        call sync_controls(s)
        call show_ready(s)
        call repaint_all(s)
    end procedure select_patch

    module procedure refresh_inspector
        character(len=:),allocatable::parameter_summary
        integer::i,variant
        logical::was_updating
        integer,parameter::icons(3)=[23,24,8]
        select case(s%entity_kind)
        case(1)
            parameter_summary='1 parameter · U'
            call s%panel_titles(5)%setText('Parameter Space (U)')
            call s%domain%setAccessibleName('Curve parameter domain U')
        case(3)
            parameter_summary='3 parameters · U / V / W'
            call s%panel_titles(5)%setText('Parameter Space (U–V–W)')
            call s%domain%setAccessibleName('Three-dimensional volume parameter domain U, V and W')
        case default
            parameter_summary='2 parameters · U / V'
            call s%panel_titles(5)%setText('Parameter Space (U–V)')
            call s%domain%setAccessibleName('Surface parameter domain U and V')
        end select
        call s%inspector%context_dimensions%setText(parameter_summary)
        if(s%selected_mesh==0.and.s%selected_multipatch==0)call set_inspector_context_icon(s,s%entity_kind,icons(s%entity_kind))
        call s%probe_v%setEnabled(s%entity_kind>=2)
        call s%probe_w%setEnabled(s%entity_kind==3)
        call s%inspector%probe_v_label%setEnabled(s%entity_kind>=2)
        call s%probe_v%setVisible(s%entity_kind>=2); call s%inspector%probe_v_label%setVisible(s%entity_kind>=2)
        call s%probe_w%setVisible(s%entity_kind==3); call s%inspector%probe_w_label%setVisible(s%entity_kind==3)
        do i=1,3
            call s%knot_buttons(i)%setEnabled(i<=s%entity_kind)
            call s%knot_buttons(i)%setVisible(i<=s%entity_kind)
            call s%widgets%entity_menu_actions(i)%setEnabled(any([(document_patch(s,i,variant),variant=1,single_patch_count(i,catalog=s%document%patches))]))
        end do
        call refresh_live_data(s)
        call refresh_inline_interval(s)
        call refresh_live_inspector(s,parameter_summary)
        call refresh_multipatch_card(s)
        call refresh_shared_controls(s)
        call on_control_selection(s)
        call sync_project_meshes(s)
        call sync_patch_visibility_button(s)
        call refresh_edit_actions(s)
        if(s%entity_kind==1)then
            call s%domain_overlays(3)%setText('Span preview')
            call s%domain_overlays(5)%setText('Selected span')
        else
            call s%domain_overlays(3)%setText('Element preview')
            call s%domain_overlays(5)%setText('Selected element')
        end if
        call refresh_domain_readout(s)
        call refresh_basis_readout(s)
        was_updating=s%updating; s%updating=.true.
        call s%probe_u%setValue(s%u); call s%probe_v%setValue(s%v); call s%probe_w%setValue(s%w)
        s%updating=was_updating
        call refresh_workspace_status(s)
        call sync_studio_assistant(s)
        call sync_project_title(s)
    end procedure refresh_inspector

    subroutine refresh_basis_readout(s)
        type(studio_t),intent(in)::s
        character(len=:),allocatable::title,hint,fixed
        character::axis
        integer::i

        axis=parameter_directions(s%knot_axis:s%knot_axis)
        title='Basis Functions ('//axis//')'
        if(s%basis_derivative)title='Basis Derivative ('//axis//')'
        call s%basis_title%setText(title)
        call s%basis%setAccessibleName(title)
        hint=axis//' follows Parameter Space'
        if(s%live_probe%valid)then
            fixed=''
            do i=1,s%entity_kind
                if(i==s%knot_axis)cycle
                if(len(fixed)>0)fixed=fixed//' · '
                fixed=fixed//parameter_directions(i:i)//' '//compact_scalar(s%live_probe%parameters(i))
            end do
            if(len(fixed)>0)hint='Probe slice · '//fixed
        end if
        call set_status_hint(s%basis,hint)
    end subroutine refresh_basis_readout

    subroutine refresh_live_data(s)
        type(studio_t),intent(inout)::s
        real(dp)::normalized(3)
        integer::axis,patch_id,previous_counts(3)
        logical::moved,same_patch,slice_moved

        patch_id=active_patch_id(s%entity_kind,catalog=s%document%patches)
        if((allocated(s%selection%items).and.s%selection%count()==0).or. &
            s%selected_mesh>0.or..not.document_patch(s,s%entity_kind,patch_id))then
            s%live_scene%valid=.false.
            s%live_probe%valid=.false.
            s%live_plot%valid=.false.
            return
        end if
        same_patch=s%live_scene%dimension==s%entity_kind.and.s%live_scene%variant==patch_id
        if(.not.same_patch.or..not.s%live_scene%valid)then
            call refresh_studio_scene(s%document%scene,catalog=s%document%patches)
            previous_counts=s%live_scene%control_counts
            associate(entry=>s%document%scene%catalog(s%entity_kind)%entries(patch_id))
                if(entry%revision==single_patch_revision(s%entity_kind,patch_id,catalog=s%document%patches))then
                    s%live_scene=entry%patch
                else
                    call read_single_patch_scene(s%entity_kind,s%live_scene,patch_id,catalog=s%document%patches)
                end if
            end associate
            if(.not.same_patch.or.any(previous_counts/=s%live_scene%control_counts))then
                s%selected_control=0
                s%control_dragging=.false.
            end if
            s%live_normalized=-1._dp
            s%live_plot_axis=0
            s%live_slice_w=-1._dp
            if(allocated(s%live_slice))deallocate(s%live_slice)
        end if
        if(.not.s%live_scene%valid)then
            s%live_probe%valid=.false.
            s%live_probe%message='Patch data unavailable.'
            s%live_plot%valid=.false.
            return
        end if

        normalized=[s%u,s%v,s%w]
        moved=any(normalized(:s%entity_kind)/=s%live_normalized(:s%entity_kind))
        if(moved.or.s%live_probe_axis/=s%knot_axis)then
            call read_single_patch(s%entity_kind,normalized,s%knot_axis,s%live_probe,patch_id,catalog=s%document%patches)
            s%live_probe_axis=s%knot_axis
        end if
        axis=s%knot_axis
        ! A change along the sampled direction moves the guide only.
        slice_moved=any(normalized(:axis-1)/=s%live_normalized(:axis-1)).or. &
            any(normalized(axis+1:s%entity_kind)/=s%live_normalized(axis+1:s%entity_kind))
        if(s%live_plot_axis/=axis.or.(s%live_plot_derivative.neqv.s%basis_derivative).or.slice_moved)then
            call read_single_patch_plot(s%entity_kind,normalized,axis,s%basis_derivative,s%live_plot,patch_id,catalog=s%document%patches)
            s%live_plot_axis=axis
            s%live_plot_derivative=s%basis_derivative
        end if
        if(s%entity_kind==3.and.s%w/=s%live_slice_w)then
            call read_single_patch_slice(normalized,s%live_scene,s%live_slice,catalog=s%document%patches)
            s%live_slice_w=s%w
        end if
        s%live_normalized=normalized
    end subroutine refresh_live_data

    module procedure refresh_shared_controls
        integer :: id,previous
        integer,allocatable :: members(:),map(:)
        integer(int64),allocatable :: revisions(:)
        logical :: same_controls,unchanged

        id=s%selected_multipatch
        previous=s%shared_controls_assembly
        if(s%shared_controls%valid.and.id==previous)then
            members=s%shared_controls%members
            map=s%shared_controls%map
            revisions=s%shared_controls%revisions
        end if
        if(id>0.and..not.s%mesh_only.and.(s%layers(4).or.s%feature_pick_mode==pick_control_points))then
            call read_document_multipatch_controls(s%document,id,s%shared_controls)
        else
            s%shared_controls=multipatch_controls_t()
        end if
        s%shared_controls_assembly=id

        same_controls=.false.
        unchanged=.false.
        if(allocated(map).and.s%shared_controls%valid)then
            if(size(members)==size(s%shared_controls%members).and.size(map)==size(s%shared_controls%map))then
                same_controls=all(members==s%shared_controls%members).and.all(map==s%shared_controls%map)
                if(same_controls)unchanged=all(revisions==s%shared_controls%revisions)
            end if
        end if
        if(.not.unchanged)then
            s%canvas_gesture%viewport_shared_hover=0
            if(previous>0.and.s%control_dragging)then
                s%control_dragging=.false.
                s%control_snapped=.false.
                s%dragging=.false.
                s%canvas_gesture%pointer_button=Qt_NoButton
                call s%status%setActivity('')
            end if
        end if
        if(same_controls)return
        if(allocated(s%selected_shared_controls))deallocate(s%selected_shared_controls)
        if(id>0)s%selected_control=0
    end procedure refresh_shared_controls

    module procedure on_control_selection
        type(single_patch_controls)::controls
        integer::point
        character(:),allocatable::name

        point=s%selected_control
        if(.not.s%updating.and..not.s%applying_selection)call refresh_workspace_status(s)
        if(s%selected_multipatch>0)then
            if(s%shared_controls%valid.and.point>0)then
                if(point<=size(s%shared_controls%weights))then
                    name='Control '//integer_text(point)
                    if(s%shared_controls%copies(point)>1)name='Shared control '//integer_text(point)
                    if(allocated(s%selected_shared_controls))then
                        if(size(s%selected_shared_controls)>1)name=integer_text(size(s%selected_shared_controls))// &
                            ' selected · Active '//integer_text(point)
                    end if
                    if(active_length_unit(s)/=unitless)name=length_caption(name,active_length_unit(s))
                    call s%inspector%selected_control_label%setText(name//new_line('a')// &
                        coordinates_text(display_length(s%shared_controls%points(:,point),s%document%model_unit, &
                            active_length_unit(s)),.false.)//' · Weight '//compact_scalar(s%shared_controls%weights(point)))
                    call s%inspector%selected_control_row%show()
                    return
                end if
            end if
            call s%inspector%selected_control_row%hide()
            return
        end if
        if(point>0.and.s%live_scene%valid)then
            call read_single_patch_controls(s%entity_kind,controls,variant=s%live_scene%variant,catalog=s%document%patches)
            if(controls%valid)then
                if(point<=size(controls%weights))then
                    name='Control point ('//control_index_text(point,controls%counts,s%entity_kind)//')'
                    if(active_length_unit(s)/=unitless)name=length_caption(name,active_length_unit(s))
                    call s%inspector%selected_control_label%setText(name//new_line('a')// &
                        coordinates_text(display_length(controls%points(:,point),s%document%model_unit,active_length_unit(s)),.false.)// &
                        ' · Weight '//compact_scalar(controls%weights(point)))
                    call s%inspector%selected_control_row%show()
                    return
                end if
            end if
        end if
        call s%inspector%selected_control_row%hide()
    end procedure on_control_selection

    module procedure refresh_edit_actions
        integer::i,kind,variant,rgb(3),opacity
        logical::can_undo,can_redo,can_remove,can_edit,has_group,has_shape,one_selected,can_measure,has_guides,can_edit_group
        type(assistant_action_context) :: context

        has_group=s%selected_multipatch>0
        one_selected=.true.
        if(allocated(s%selection%items))one_selected=s%selection%count()==1
        can_edit=one_selected.and.s%live_scene%valid.and..not.has_group.and.s%selected_mesh==0
        call prune_document_redo(s%document)
        call refresh_project_tab_actions(s)
        can_undo=document_undo_stamp(s%document)>0
        can_remove=can_edit.and.can_remove_knot(s,any_direction=.true.)
        can_edit_group=.false.
        if(has_group.and.one_selected)can_edit_group=allocated(s%document%multipatches%groups(s%selected_multipatch)%native)
        call s%widgets%edit_controls_action%setEnabled(can_edit.or.can_edit_group)
        call s%widgets%undo_patch_action%setEnabled(can_undo)
        can_redo=document_can_redo(s%document)
        call s%toolbar_actions(5)%setEnabled(can_redo)
        call s%widgets%ribbon%alternate_redo_action%setEnabled(can_redo)
        call s%widgets%insert_knot_action%setEnabled(can_edit)
        call s%widgets%insert_knot_button%setEnabled(can_edit)
        call s%widgets%remove_knot_action%setEnabled(can_remove)
        call s%widgets%remove_knot_button%setEnabled(can_remove)
        call s%widgets%elevate_degree_action%setEnabled(can_edit)
        call s%widgets%reduce_degree_action%setEnabled(can_edit.and.can_reduce_degree(s))
        call s%widgets%domain_knots_button%setEnabled(can_edit)
        call s%widgets%edit_knot_values_action%setEnabled(can_remove)
        call s%widgets%domain_values_button%setEnabled(can_remove)
        call s%widgets%domain_refine_button%setEnabled(can_edit)
        call s%widgets%ribbon%refinement_action%setEnabled(can_edit)
        call s%widgets%ribbon%patch_tool_actions(1)%setEnabled(can_edit.and.can_restructure_patch(s))
        call s%widgets%ribbon%patch_tool_actions(2)%setEnabled(can_edit)
        call s%widgets%ribbon%patch_tool_actions(3)%setEnabled(can_edit.and.can_restructure_patch(s))
        call s%widgets%ribbon%patch_tool_actions(4)%setEnabled(can_edit.and.s%entity_kind>=2)
        call s%widgets%ribbon%patch_tool_actions(5)%setEnabled(can_edit.and.s%entity_kind==3)
        call s%widgets%domain_knot_buttons(3)%setEnabled(can_edit.and.can_remove_knot(s))
        call s%widgets%edit_interval_action%setEnabled(can_edit)
        call s%widgets%ribbon%parameters_ribbon%setEnabled(can_edit)
        do i=1,size(s%basis_buttons)
            call s%basis_buttons(i)%setEnabled(can_edit)
        end do
        call update_assembly_ribbon(s,has_group.and.one_selected)
        call update_mesh_ribbon(s)
        has_guides=selected_controls_available(s)
        call s%viewport_presets(2)%setEnabled(has_guides)
        do i=3,size(s%viewport_layers)
            call s%viewport_layers(i)%setEnabled(has_guides)
        end do
        has_shape=one_selected.and.s%selected_mesh==0.and.(has_group.or.s%live_scene%valid)
        can_measure=has_measurable_selection(s)
        call s%widgets%ribbon%measure_action%setEnabled(can_measure)
        call s%panel_actions(5)%setEnabled(.not.has_group)
        call s%panel_actions(6)%setEnabled(.not.has_group)
        context=available_assistant_tools(s)
        call s%widgets%ribbon%profile_actions(1)%setEnabled(context%has_curve.or.context%has_surface)
        call s%widgets%ribbon%profile_actions(2)%setEnabled(context%has_curve.or.context%has_surface)
        call s%widgets%ribbon%profile_actions(3)%setEnabled(context%can_sweep)
        call s%widgets%ribbon%profile_actions(4)%setEnabled(context%can_loft)
        kind=s%entity_kind
        variant=active_patch_id(kind,catalog=s%document%patches)
        if(has_group)then
            kind=s%document%multipatches%groups(s%selected_multipatch)%definition%kind
            variant=s%document%multipatches%groups(s%selected_multipatch)%definition%members(1)
        end if
        call get_studio_material(s%document%scene,kind,variant,rgb,opacity)
        call s%widgets%translucent_action%setEnabled(has_shape)
        call s%widgets%translucent_action%setChecked(has_shape.and.opacity<255)
        call s%widgets%material_color_action%setEnabled(has_shape)
        call s%widgets%transform_action%setEnabled(has_shape)
        do i=1,size(s%widgets%ribbon%transform_tools)
            call s%widgets%ribbon%transform_tools(i)%setEnabled(has_shape)
        end do
        call s%widgets%export_iges_action%setEnabled(can_edit.and.kind<3)
        call s%widgets%export_vtk_action%setEnabled(has_shape)
        if(can_remove)then
            call set_status_hint(s%widgets%remove_knot_button,'Choose a direction and remove an interior knot')
        else
            call set_status_hint(s%widgets%remove_knot_button,'This patch has no interior knots')
        end if
    end procedure refresh_edit_actions

    logical function has_measurable_selection(s) result(available)
        use studio_selection, only: pick_vertex
        type(studio_t), intent(in) :: s
        type(selection_target_t) :: target
        integer :: i

        if(.not.allocated(s%selection%items))then
            target=primary_workspace_target(s)
            available=measurable(target)
            return
        end if
        available=.false.
        do i=1,s%selection%count()
            available=measurable(s%selection%items(i))
            if(available)return
        end do
    contains
        logical function measurable(target) result(matches)
            type(selection_target_t),intent(in) :: target

            matches=target%category==selection_patch.and.target%feature/=pick_vertex
            if(target%category/=selection_assembly)return
            if(target%id<1.or.target%id>count_multipatches(s%document%multipatches))return
            matches=size(s%document%multipatches%groups(target%id)%definition%members)>0
        end function measurable
    end function has_measurable_selection

    module procedure can_reduce_degree

        available=s%live_scene%valid.and.any(s%live_scene%degree(:s%entity_kind)>0)
    end procedure can_reduce_degree

    module procedure can_remove_knot
        integer::axis,first,last

        available=.false.
        if(.not.s%live_scene%valid)return
        first=s%knot_axis
        last=first
        if(present(any_direction))then
            if(any_direction)then
                first=1
                last=s%entity_kind
            end if
        end if
        do axis=first,last
            if(.not.allocated(s%live_scene%axes(axis)%knots))cycle
            associate(knots=>s%live_scene%axes(axis)%knots,bounds=>s%live_scene%bounds(:,axis))
                available=any(knots>bounds(1).and.knots<bounds(2))
            end associate
            if(available)return
        end do
    end procedure can_remove_knot

    module procedure refresh_studio_diagnostics
        character(:),allocatable :: issue

        issue=studio_diagnostic_summary(s)
        call s%inspector%repair_diagnostic%setVisible(len(issue)>0)
        if(len(issue)==0)issue='No issues reported for this selection.'
        call s%diagnostic_label%setText(issue)
    end procedure refresh_studio_diagnostics

    subroutine refresh_live_inspector(s,parameter_summary)

        type(studio_t),intent(inout)::s
        character(len=*),intent(in)::parameter_summary
        character(len=:),allocatable::text,parameters,degrees,controls,note
        character(len=3),parameter::axes='UVW'
        character(len=3),parameter::derivative_axes='uvw'
        integer::i

        call refresh_studio_diagnostics(s)
        call s%inspector%sections(1)%setVisible(s%selected_mesh==0.and.s%selected_multipatch==0)
        if(s%selected_multipatch>0)return
        associate(info=>s%live_probe,scene=>s%live_scene)
        if(s%selected_mesh>0)then
            associate(mesh=>s%document%meshes(s%selected_mesh))
                call s%inspector%context_identity%setText(mesh%name)
                text=entity_name(mesh%kind)//' mesh · '//integer_text(mesh%nodes)//' nodes · '//integer_text(mesh%elements)//' cells'
                if(.not.studio_mesh_valid(s,s%selected_mesh))text=text//new_line('a')//'Needs update · Regenerate from Mesh'
                call s%inspector%context_dimensions%setText(text)
            end associate
            call set_inspector_context_icon(s,5,icon_mesh)
            call set_probe_readout(s%probe_label,'Select a patch to inspect parameters.')
            call set_status_hint(s%probe_label,'')
            return
        end if
        if(.not.scene%valid)then
            call s%inspector%context_identity%setText('No patch selected')
            call s%inspector%context_dimensions%setText('Create or select geometry to inspect it.')
            call set_probe_readout(s%probe_label,'Select a patch to inspect parameters.')
            call set_status_hint(s%probe_label,'')
            return
        end if
        call s%inspector%context_identity%setText(single_patch_name(s%entity_kind,active_patch_id(s%entity_kind,catalog=s%document%patches),catalog=s%document%patches))

        degrees=''
        controls=''
        do i=1,scene%dimension
            if(i>1)then
                degrees=degrees//' × '
                controls=controls//' × '
            end if
            degrees=degrees//integer_text(scene%degree(i))
            controls=controls//integer_text(scene%control_counts(i))
        end do
        call s%inspector%context_dimensions%setText(parameter_summary//new_line('a')//'Degree '//degrees// &
            new_line('a')//'Control points '//controls)
        if(.not.info%valid)then
            note='Probe data unavailable.'
            if(allocated(info%message))then
                if(len(info%message)>0)note=info%message
            end if
            call set_probe_readout(s%probe_label,'Probe unavailable.')
            call set_status_hint(s%probe_label,note)
            return
        end if
        parameters=''
        do i=1,info%dimension
            if(i>1)parameters=parameters//' · '
            parameters=parameters//axes(i:i)//' '//fixed(info%parameters(i))
        end do
        text='<p style="margin:3px 0;color:#bdd2de;font-size:12px"><b>Patch parameters</b> · '//parameters//'</p>'// &
            '<table cellspacing="0" cellpadding="3" width="100%">'// &
            html_row(length_caption('Point',active_length_unit(s)), &
                coordinates_text(display_length(info%point,s%document%model_unit,active_length_unit(s)),.false.))
        note=length_caption('Point',active_length_unit(s))//' '// &
            coordinates_text(display_length(info%point,s%document%model_unit,active_length_unit(s)),.true.)
        do i=1,info%dimension
            text=text//html_row('∂P/∂'//derivative_axes(i:i), &
                coordinates_text(display_length(info%jacobian(:,i),s%document%model_unit,active_length_unit(s)),.false.))
            note=note//new_line('a')//'∂P/∂'//derivative_axes(i:i)//': '// &
                coordinates_text(display_length(info%jacobian(:,i),s%document%model_unit,active_length_unit(s)),.true.)
        end do
        call set_probe_readout(s%probe_label,text//'</table>')
        call set_status_hint(s%probe_label,note)
        end associate
    end subroutine refresh_live_inspector

    subroutine set_probe_readout(label,text)
        type(QLabel), intent(inout) :: label
        character(*), intent(in) :: text

        call label%setText(text)
        associate(hint=>label%minimumSizeHint())
            call label%setMinimumHeight(hint%height)
        end associate
    end subroutine set_probe_readout

    module procedure show_multipatch_inspector_summary
        character(:),allocatable :: details
        integer :: id

        id=s%selected_multipatch
        if(id<1.or.id>count_multipatches(s%document%multipatches))return
        associate(group=>s%document%multipatches%groups(id))
            call s%inspector%context_identity%setText(group%definition%name)
            details=entity_name(group%definition%kind)//' multipatch · '//integer_text(size(group%definition%members))//' members · '// &
                integer_text(group%validation%connections)//' interface'
            if(group%validation%connections/=1)details=details//'s'
            if(group%validation%valid.and.group%validation%controls>0)details=details//new_line('a')// &
                integer_text(group%validation%shared_controls)//' unique controls of '//integer_text(group%validation%controls)
            if(allocated(group%validation%message))then
                if(len(group%validation%message)>0)details=details//new_line('a')//group%validation%message
            end if
            call s%inspector%context_dimensions%setText(details)
        end associate
        call set_inspector_context_icon(s,4,87)
    end procedure show_multipatch_inspector_summary

    subroutine set_inspector_context_icon(s,kind,glyph)
        type(studio_t), intent(inout) :: s
        integer, intent(in) :: kind,glyph
        type(QIcon) :: icon

        if(s%inspector%context_kind==kind)return
        icon=line_icon(glyph,s%inspector%context_icon)
        call icon%delete()
        s%inspector%context_kind=kind
    end subroutine set_inspector_context_icon

    module procedure refresh_inline_interval
        integer::axis,context(3)
        integer(int64)::revision
        logical::available

        axis=s%knot_axis
        available=s%live_scene%valid.and.s%selected_multipatch==0
        call s%widgets%interval_lower%setEnabled(available)
        call s%widgets%interval_upper%setEnabled(available)
        call s%widgets%interval_apply%setEnabled(available)
        call s%widgets%interval_label%setText(parameter_directions(axis:axis)//' range')
        if(.not.available)then
            s%widgets%interval_context=0
            call s%widgets%interval_lower%setText('')
            call s%widgets%interval_upper%setText('')
            return
        end if
        context=[s%entity_kind,s%live_scene%variant,axis]
        revision=s%document%patches%patch_revision(context(1),context(2))
        if(all(s%widgets%interval_context==context).and.s%widgets%interval_revision==revision)return
        s%widgets%interval_context=context
        s%widgets%interval_revision=revision
        call s%widgets%interval_lower%setText(number_text(s%live_scene%bounds(1,axis)))
        call s%widgets%interval_upper%setText(number_text(s%live_scene%bounds(2,axis)))
    end procedure refresh_inline_interval

    module procedure refresh_domain_readout
        character(:),allocatable::hint,selected_knot
        character::axis

        if(.not.s%live_probe%valid.or.s%selected_mesh>0.or.s%selected_multipatch>0)then
            call set_status_hint(s%domain,'Select a patch to explore its parameter space')
            return
        end if
        hint='Drag to probe · Scroll to zoom · Middle-drag to pan'
        if(s%live_probe%dimension==3)hint=hint//' · Shift+middle or right drag to orbit'
        hint=hint//' · Shift to edit knots'
        if(s%domain_edit_knots)then
            axis=parameter_directions(s%knot_axis:s%knot_axis)
            select case(s%domain_knot_tool)
            case(knot_insert_tool)
                hint='Click a gap on the '//axis//' ruler to insert a knot'
            case(knot_remove_tool)
                hint='Click an interior '//axis//' knot to remove it'
            case default
                hint='Drag a knot · Double-click a gap to insert · Delete to remove'
                selected_knot=knot_status_text(s)
                if(len(selected_knot)>0)hint=selected_knot//' · '//hint
            end select
        end if
        call set_status_hint(s%domain,hint)
    end procedure refresh_domain_readout

    module procedure refresh_knot_table
        character::axis_name
        type(QHeaderView)::header,rows
        integer::row_height,visible_groups
        character(*),parameter::indices_title='Indices',multiplicity_title='Mult.'

        axis_name=achar(iachar('U')+s%knot_axis-1)
        call s%knot_table%setHorizontalHeaderLabels([character(len=12)::indices_title,'Knot '//axis_name,multiplicity_title])
        call refresh_live_knots(s,axis_name)
        header=s%knot_table%horizontalHeader()
        call header%setFixedHeight(max(22,header%lineHeight()+6))
        call header%resizeSection(0,max(56,header%textWidth(indices_title)+12))
        call header%resizeSection(2,max(56,header%textWidth(multiplicity_title)+12))
        rows=s%knot_table%verticalHeader()
        row_height=max(23,s%knot_table%lineHeight()+6)
        call rows%setMinimumSectionSize(row_height)
        call rows%setDefaultSectionSize(row_height)
        visible_groups=min(3,s%knot_table%rowCount())
        call s%knot_table%setFixedHeight(header%height()+2+visible_groups*row_height)
        call s%knot_table%scrollToTop()
        call s%knot_table%setCurrentCell(-1,-1)
    end procedure refresh_knot_table

    subroutine refresh_live_knots(s,axis_name)
        type(studio_t),intent(inout)::s
        character(len=1),intent(in)::axis_name
        character(len=:),allocatable::indices,detail
        character(len=26)::precise_knot
        integer::n,last,first,row

        n=0
        if(s%live_scene%valid.and.s%live_scene%dimension==s%entity_kind)then
            if(allocated(s%live_scene%axes(s%knot_axis)%knots))n=size(s%live_scene%axes(s%knot_axis)%knots)
        end if
        if(n==0)then
            s%live_knots=[real(dp)::]
            s%live_knot_groups=[real(dp)::]
            call s%domain%update()
            call s%knot_table%setRowCount(0)
            detail='Knot data unavailable.'
            call set_status_hint(s%knot_table,detail)
            return
        end if
        s%live_knots=s%live_scene%axes(s%knot_axis)%knots
        associate(knots=>s%live_knots)
        s%live_knot_groups=pack(knots,[.true.,knots(2:n)/=knots(:n-1)])
        call s%domain%update()
        call set_status_hint(s%knot_table,integer_text(n)//' stored '//axis_name// &
            ' knots. Identical values are grouped. Double-click an interior group to edit knot values.')
        call s%knot_table%setRowCount(size(s%live_knot_groups))
        first=1
        row=0
        do last=1,n
            if(last<n)then
                if(knots(last)==knots(last+1))cycle
            end if
            indices=integer_text(first-1)
            if(last>first)indices=indices//'–'//integer_text(last-1)
            write(precise_knot,'(es26.17e3)')knots(first)
            call set_knot_item(s%knot_table,row,0,indices)
            call set_knot_item(s%knot_table,row,1,compact_scalar(knots(first)),trim(adjustl(precise_knot)))
            call set_knot_item(s%knot_table,row,2,integer_text(last-first+1))
            first=last+1
            row=row+1
        end do
        end associate
    end subroutine refresh_live_knots

    subroutine set_knot_item(table,row,column,text,detail)
        type(QTableWidget),intent(in)::table
        integer,intent(in)::row,column
        character(len=*),intent(in)::text
        character(len=*),intent(in),optional::detail
        type(QTableWidgetItem)::item
        item=QTableWidgetItem(text)
        call item%setFlags(ior(Qt_ItemIsEnabled,Qt_ItemIsSelectable))
        call item%setTextAlignment(Qt_AlignCenter)
        if(present(detail))then
            call set_status_hint(item,detail)
        else
            call set_status_hint(item,text)
        end if
        call table%setItem(row,column,item)
    end subroutine set_knot_item

    module procedure on_knot_double_click
        real(dp)::value

        if(column/=1)return
        select type(s=>data)
        type is(studio_t)
            if(.not.allocated(s%live_knot_groups))return
            if(row<0.or.row>=size(s%live_knot_groups))return
            value=s%live_knot_groups(row+1)
            associate(bounds=>s%live_scene%bounds(:,s%knot_axis))
                if(value>bounds(1).and.value<bounds(2))call run_knot_values_editor(s)
            end associate
        end select
    end procedure on_knot_double_click

    module procedure on_knot_table_selection
        character(:),allocatable :: text

        if(row==previous_row)return
        select type(s=>data)
        type is(studio_t)
            call s%domain%update()
            call refresh_domain_readout(s)
            if(s%updating.or..not.(s%domain%hasFocus().or.s%knot_table%hasFocus()))return
            text=knot_status_text(s)
            if(len(text)>0)call s%status%showMessage(text,3800,record=.false.)
        end select
    end procedure on_knot_table_selection

    module procedure on_domain
        integer::i
        select type(s=>data)
        type is(studio_t)
            if(s%updating)return
            do i=1,5
                s%domain_layers(i)=s%domain_overlays(i)%isChecked()
            end do
            call sync_controls(s); call s%domain%update()
        end select
    end procedure on_domain

    module procedure on_probe
        select type(s=>data)
        type is(studio_t)
            if(s%updating)return
            call set_parameter_probe(s,[s%probe_u%value(),s%probe_v%value(),s%probe_w%value()])
        end select
    end procedure on_probe

    module procedure set_parameter_probe
        character(:),allocatable :: text

        s%u=normalized(1); s%v=normalized(2); s%w=normalized(3)
        call refresh_inspector(s)
        call repaint_all(s)
        text=probe_status_text(s)
        if(len(text)>0)call s%status%showMessage(text,3500,record=.false.)
    end procedure set_parameter_probe
end submodule studio_ui_inspection
