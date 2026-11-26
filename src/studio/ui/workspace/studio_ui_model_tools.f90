submodule(studio_ui) studio_ui_model_tools
    use studio_commands, only: mesh_ribbon_tab
    use studio_state, only: scene_mesh_drawn,scene_patch_drawn,scene_drawn_patch_count
    use studio_mesh_operations, only: estimate_document_mesh_samples
    use studio_multipatch_collection, only: publish_multipatch_volume
    use studio_workspace_selection, only: selection_mesh_records
    implicit none(type, external)
contains
    module procedure update_assembly_ribbon
        integer::i
        integer,allocatable::mesh_ids(:)
        logical::can_connect,can_fit,members_shown,can_refine
        character(:),allocatable::fit_caption,fit_hint,visibility_caption

        can_connect=any(eligible_assembly_ranks(s))
        call s%widgets%ribbon%assembly_connect_action%setEnabled(can_connect)
        can_fit=.false.
        members_shown=.false.
        can_refine=.false.
        fit_caption='Fit members'
        fit_hint='Frame the visible members of the selected multipatch.'
        if(s%mesh_only)then
            fit_caption='Fit mesh'
            fit_hint='Frame the visible meshes of the selected multipatch.'
        end if
        if(has_group)then
            can_refine=allocated(s%document%multipatches%groups(s%selected_multipatch)%native)
            if(s%mesh_only)then
                mesh_ids=selection_mesh_records(s)
                can_fit=any([(scene_mesh_drawn(s,mesh_ids(i)),i=1,size(mesh_ids))])
                if(.not.can_fit)fit_hint='Show a mesh for this multipatch first.'
            else
                associate(group=>s%document%multipatches%groups(s%selected_multipatch)%definition)
                    can_fit=any([(scene_patch_drawn(s,group%kind,group%members(i)),i=1,size(group%members))])
                end associate
                if(.not.can_fit)fit_hint='Show a member of this multipatch first.'
            end if
            members_shown=multipatch_members_shown(s)
        else
            fit_hint='Select one multipatch to fit.'
        end if
        associate(action=>s%widgets%ribbon%assembly_fit_action)
            call action%setEnabled(can_fit)
            call action%setText(fit_caption)
            call action%setIconText(fit_caption)
            call set_status_hint(action,studio_tooltip(fit_caption,fit_hint))
        end associate
        visibility_caption='Show members'
        if(members_shown)visibility_caption='Hide members'
        associate(action=>s%widgets%ribbon%assembly_visibility_action)
            call action%setEnabled(has_group)
            call action%setText(visibility_caption)
            call action%setIconText(visibility_caption)
            call set_status_hint(action,studio_tooltip(visibility_caption,'Change visibility for every member of the selected multipatch.'))
            if(members_shown)then
                call action%setIcon(s%widgets%project_browser%visible_icon)
            else
                call action%setIcon(s%widgets%project_browser%hidden_icon)
            end if
        end associate
        do i=1,size(s%widgets%ribbon%assembly_refine_actions)
            call s%widgets%ribbon%assembly_refine_actions(i)%setEnabled(can_refine)
        end do
        do i=1,size(s%widgets%multipatch_actions)
            call s%widgets%multipatch_actions(i)%setEnabled(has_group)
        end do
        if(can_connect)then
            call set_status_hint(s%widgets%ribbon%assembly_connect_action,studio_tooltip('New multipatch', &
                'Group patches of one type and define interfaces where their boundaries meet.','Ctrl+Shift+M'))
        else
            call set_status_hint(s%widgets%ribbon%assembly_connect_action,studio_tooltip('New multipatch', &
                'The project needs at least two patches of the same type.'))
        end if
    end procedure update_assembly_ribbon

    function eligible_assembly_ranks(s) result(eligible)
        type(studio_t), intent(in) :: s
        logical :: eligible(3)
        integer :: rank,i

        do rank=1,3
            eligible(rank)=count([(document_patch(s,rank,i),i=1,single_patch_count(rank,catalog=s%document%patches))])>=2
        end do
    end function eligible_assembly_ranks

    function initial_assembly_members(s,rank) result(members)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: rank
        integer, allocatable :: members(:)
        logical :: selected(single_patch_count(rank,catalog=s%document%patches))
        type(selection_target_t) :: target
        integer :: i

        selected=.false.
        if(allocated(s%selection%items))then
            do i=1,s%selection%count()
                target=s%selection%items(i)
                if(target%category/=selection_patch.or.target%rank/=rank)cycle
                if(document_patch(s,rank,target%id))selected(target%id)=.true.
            end do
        else
            target=primary_workspace_target(s)
            if(target%category==selection_patch.and.target%rank==rank)then
                if(document_patch(s,rank,target%id))selected(target%id)=.true.
            end if
        end if
        members=pack([(i,i=1,size(selected))],selected)
    end function initial_assembly_members

    module procedure update_mesh_ribbon
        integer :: id
        integer, allocatable :: mesh_ids(:)
        type(mesh_source_t) :: source
        character(:),allocatable :: detail,build_caption,pick_hint,visibility_hint,only_hint
        logical :: available,current,shown,can_pick,one_selected

        source=mesh_source_selection(s)
        id=source%mesh_id
        one_selected=.true.
        if(allocated(s%selection%items))one_selected=s%selection%count()==1
        available=one_selected.and.mesh_source_available(s,source)
        current=one_selected.and.studio_mesh_valid(s,id)
        mesh_ids=selection_mesh_records(s)
        shown=.false.
        visibility_hint='Select a mesh, cell, or its source geometry.'
        if(size(mesh_ids)>0)then
            shown=all(s%document%meshes(mesh_ids)%visible)
            visibility_hint='Show the meshes of the selection.'
            if(shown)visibility_hint='Hide the meshes of the selection.'
        end if
        can_pick=has_visible_meshes(s)
        call s%widgets%ribbon%mesh_generate_action%setEnabled(available)
        if(id>0)then
            build_caption='Remesh'
            detail='Set the sizing and replace the existing mesh.'
        else
            build_caption='Generate mesh'
            detail='Set the sizing and generate a structured mesh.'
        end if
        if(.not.available)detail='Select one patch, valid multipatch, or mesh with an available source.'
        call s%widgets%ribbon%mesh_generate_action%setIconText(build_caption//'…')
        call s%widgets%ribbon%mesh_generate_action%setText(build_caption//'…')
        call set_status_hint(s%widgets%ribbon%mesh_generate_action,studio_tooltip(build_caption,detail))
        call set_status_hint(s%widgets%ribbon%mesh_edges_action,studio_tooltip('Show meshes',visibility_hint))
        call s%widgets%ribbon%mesh_edges_action%setEnabled(size(mesh_ids)>0)
        call s%widgets%ribbon%mesh_edges_action%setChecked(shown)
        call s%widgets%ribbon%mesh_export_action%setEnabled(current)
        call s%widgets%ribbon%mesh_pick_action%setEnabled(can_pick)
        pick_hint='Show a mesh to enable cell picking.'
        if(can_pick)pick_hint='Click a visible cell to inspect it in the viewport.'
        if(s%mesh_pick)pick_hint='Click a cell to inspect it. Turn off Pick cells to select whole objects.'
        call set_status_hint(s%widgets%ribbon%mesh_pick_action,studio_tooltip('Pick cells',pick_hint))
        call s%widgets%ribbon%mesh_pick_action%setChecked(s%mesh_pick)
        call sync_selection_filter(s)
        call s%widgets%ribbon%mesh_only_action%setEnabled(can_pick.or.s%mesh_only)
        call s%widgets%ribbon%mesh_only_action%setChecked(s%mesh_only)
        only_hint='Show a mesh in the current view first.'
        if(can_pick)only_hint='Show all visible meshes without CAD geometry.'
        if(s%mesh_only)only_hint='Return to the geometry and mesh view.'
        call set_status_hint(s%widgets%ribbon%mesh_only_action,studio_tooltip('Meshes only',only_hint))
        call s%widgets%ribbon%mesh_shaded_action%setEnabled(has_visible_meshes(s,2))
        call s%widgets%ribbon%mesh_shaded_action%setChecked(s%mesh_shaded)
        call s%widgets%ribbon%physical_knot_action%setEnabled(.not.s%mesh_only)
        call s%widgets%ribbon%physical_knot_action%setChecked(s%physical_knots)
        call sync_volume_mesh_controls(s)
    end procedure update_mesh_ribbon

    function mesh_source_revisions(s,source) result(revisions)
        type(studio_t), intent(in) :: s
        type(mesh_source_t), intent(in) :: source
        integer(int64), allocatable :: revisions(:)
        integer :: i,offset

        revisions=[(single_patch_revision(source%kind,source%members(i),catalog=s%document%patches),i=1,size(source%members))]
        if(source%assembly_id==0)return
        associate(assembly=>s%document%multipatches%groups(source%assembly_id)%definition)
            if(.not.allocated(assembly%joins))return
            offset=size(revisions)
            revisions=[revisions,spread(0_int64,1,8*size(assembly%joins))]
            do i=1,size(assembly%joins)
                associate(join=>assembly%joins(i))
                    revisions(offset+1:offset+8)=[int([join%patch_a,join%side_a,join%patch_b,join%side_b],int64), &
                        merge(1_int64,0_int64,[join%reverse,join%swap,join%flip])]
                    offset=offset+8
                end associate
            end do
        end associate
    end function mesh_source_revisions

    logical function mesh_preview_current(s,source,recorded) result(current)
        type(studio_t), intent(in) :: s
        type(mesh_source_t), intent(in) :: source
        integer(int64), intent(in) :: recorded(:)
        integer(int64), allocatable :: revisions(:)

        current=.false.
        if(.not.mesh_source_available(s,source))return
        revisions=mesh_source_revisions(s,source)
        if(size(revisions)/=size(recorded))return
        current=all(revisions==recorded)
    end function mesh_preview_current

    logical function mesh_source_available(s,source) result(available)
        type(studio_t), intent(in) :: s
        type(mesh_source_t), intent(in) :: source

        available=.false.
        if(source%kind<1.or.source%kind>3.or..not.allocated(source%members))return
        if(size(source%members)==0)return
        if(source%assembly_id>0)then
            if(source%assembly_id>count_multipatches(s%document%multipatches))return
            available=allocated(s%document%multipatches%groups(source%assembly_id)%native)
        else
            available=.true.
        end if
    end function mesh_source_available

    function mesh_sizing(s,source) result(settings)
        type(studio_t), intent(in) :: s
        type(mesh_source_t), intent(in) :: source
        type(mesh_sizing_t) :: settings
        integer :: slot,rank,part,i,limit

        slot=mesh_draft_index(s,source)
        if(slot>0)then
            settings=s%widgets%mesh_drafts(slot)
            return
        end if
        rank=source%kind
        if(rank<1.or.rank>3.or.size(source%members)==0)then
            allocate(settings%cells(0,0))
            return
        end if
        settings%rank=rank
        settings%assembly_id=source%assembly_id
        settings%document_revision=s%document%document_revision
        settings%members=source%members
        limit=max(1,mesh_axis_cell_limit(rank,size(source%members)))
        settings%cells=spread(min(s%widgets%mesh_divisions(:rank,rank),limit),2,size(source%members))
        if(source%mesh_id==0)return
        associate(mesh=>s%document%meshes(source%mesh_id))
            settings%conforming=allocated(mesh%connected_cells)
            settings%per_patch=mesh%per_patch_sizing
            settings%target_size=real(mesh%target_size,dp)
            settings%synchronize=mesh%synchronize
            if(.not.allocated(mesh%parts))return
            do i=1,size(source%members)
                part=0
                if(i<=size(mesh%parts))then
                    if(mesh%parts(i)%variant==source%members(i))part=i
                end if
                if(part==0)part=findloc(mesh%parts%variant,source%members(i),dim=1)
                if(part>0)settings%cells(:,i)=mesh%parts(part)%samples(:rank)-1
            end do
        end associate
        settings%per_patch=settings%per_patch.or.any(settings%cells/=spread(settings%cells(:,1),2,size(source%members)))
    end function mesh_sizing

    integer function mesh_draft_index(s,source) result(index)
        type(studio_t), intent(in) :: s
        type(mesh_source_t), intent(in) :: source
        integer :: i

        index=0
        if(.not.allocated(s%widgets%mesh_drafts))return
        do i=1,size(s%widgets%mesh_drafts)
            associate(draft=>s%widgets%mesh_drafts(i))
                if(draft%document_revision/=s%document%document_revision.or.draft%rank/=source%kind.or. &
                    draft%assembly_id/=source%assembly_id.or..not.allocated(draft%members))cycle
                if(size(draft%members)/=size(source%members))cycle
                if(any(draft%members/=source%members))cycle
                index=i
                return
            end associate
        end do
    end function mesh_draft_index

    subroutine set_mesh_draft(s,source,settings)
        type(studio_t), intent(inout) :: s
        type(mesh_source_t), intent(in) :: source
        type(mesh_sizing_t), intent(in) :: settings
        integer :: slot

        if(.not.allocated(s%widgets%mesh_drafts))allocate(s%widgets%mesh_drafts(0))
        slot=mesh_draft_index(s,source)
        if(slot==0)slot=findloc(s%widgets%mesh_drafts%document_revision/=s%document%document_revision,.true.,dim=1)
        if(slot==0)then
            s%widgets%mesh_drafts=[s%widgets%mesh_drafts,settings]
        else
            s%widgets%mesh_drafts(slot)=settings
        end if
    end subroutine set_mesh_draft

    subroutine sync_volume_mesh_controls(s)
        type(studio_t),intent(inout)::s
        integer::extent

        extent=volume_mesh_extent(s,s%mesh_layer_axis)
        s%widgets%syncing_volume_mesh_controls=.true.
        call s%widgets%ribbon%volume_mesh_controls%setEnabled(extent>0)
        call s%widgets%ribbon%volume_mesh_modes(s%mesh_view+1)%setChecked(.true.)
        call s%widgets%ribbon%volume_mesh_layer_controls%setVisible(extent>0.and.s%mesh_view==mesh_layer)
        call s%widgets%ribbon%volume_mesh_axis%setCurrentIndex(s%mesh_layer_axis-1)
        call s%widgets%ribbon%volume_mesh_layer%setRange(1,max(1,extent))
        call s%widgets%ribbon%volume_mesh_layer%setValue(s%mesh_layer_index)
        call sync_canvas_mesh_controls(s,extent)
        s%widgets%syncing_volume_mesh_controls=.false.
    end subroutine sync_volume_mesh_controls

    subroutine sync_canvas_mesh_controls(s,extent)
        type(studio_t), intent(inout) :: s
        integer, intent(in) :: extent
        integer :: layer,offset

        associate(hud=>s%widgets%canvas_mesh)
            if(hud%panel%isNull())return
            call hud%panel%setVisible(extent>0)
            if(extent<=0)return
            layer=min(s%mesh_layer_index,extent)
            offset=nint(100._dp*s%mesh_section_fraction)
            call hud%presentation%setCurrentIndex(s%mesh_view)
            call hud%layer_controls%setVisible(s%mesh_view==mesh_layer)
            call hud%layer_axis%setCurrentIndex(s%mesh_layer_axis-1)
            call hud%layer_index%setRange(1,extent)
            call hud%layer_index%setValue(layer)
            call hud%layer_value%setText(integer_text(layer)//' / '//integer_text(extent))
            call hud%section_axis%setCurrentIndex(s%mesh_section_axis)
            call hud%section_position%setVisible(s%mesh_section_axis>0)
            call hud%section_offset%setValue(offset)
            call hud%section_value%setText(integer_text(offset)//'%')
        end associate
    end subroutine sync_canvas_mesh_controls

    integer function volume_mesh_extent(s,axis) result(extent)
        type(studio_t),intent(in)::s
        integer,intent(in)::axis
        integer::id,part

        extent=0
        do id=1,studio_mesh_count(s)
            if(s%document%meshes(id)%kind/=3.or..not.scene_mesh_drawn(s,id))cycle
            associate(mesh=>s%document%meshes(id))
                do part=1,size(mesh%parts)
                    extent=max(extent,mesh%parts(part)%samples(axis)-1)
                end do
            end associate
        end do
    end function volume_mesh_extent

    module procedure volume_mesh_style_changed

        select type(s=>data)
        type is(studio_t)
            if(s%widgets%syncing_volume_mesh_controls.or.index<mesh_exterior.or.index>mesh_layer)return
            if(index==s%mesh_view)return
            s%mesh_view=index
            call refresh_volume_mesh_view(s)
        end select
    end procedure volume_mesh_style_changed

    module procedure volume_mesh_axis_changed

        select type(s=>data)
        type is(studio_t)
            if(s%widgets%syncing_volume_mesh_controls.or.index<0.or.index>2)return
            s%mesh_layer_axis=index+1
            s%mesh_layer_index=min(s%mesh_layer_index,max(1,volume_mesh_extent(s,index+1)))
            call refresh_volume_mesh_view(s)
        end select
    end procedure volume_mesh_axis_changed

    module procedure volume_mesh_layer_changed

        select type(s=>data)
        type is(studio_t)
            if(s%widgets%syncing_volume_mesh_controls.or.index<1)return
            s%mesh_layer_index=index
            call refresh_volume_mesh_view(s)
        end select
    end procedure volume_mesh_layer_changed

    module procedure mesh_section_axis_changed
        character(:), allocatable :: error

        select type(s=>data)
        type is(studio_t)
            if(s%widgets%syncing_volume_mesh_controls)return
            call set_workspace_mesh_section(s,index,nint(100._dp*s%mesh_section_fraction),error)
        end select
    end procedure mesh_section_axis_changed

    module procedure mesh_section_offset_changed
        character(:), allocatable :: error

        select type(s=>data)
        type is(studio_t)
            if(s%widgets%syncing_volume_mesh_controls)return
            call set_workspace_mesh_section(s,s%mesh_section_axis,index,error)
        end select
    end procedure mesh_section_offset_changed

    module procedure set_workspace_mesh_section

        error='Show a current volume mesh before changing its section.'
        if(.not.has_visible_meshes(s,3))return
        error='Choose off, X, Y or Z and a section percentage from 0 to 100.'
        if(axis<0.or.axis>3.or.percent<0.or.percent>100)return
        s%mesh_section_axis=axis
        s%mesh_section_fraction=real(percent,dp)/100._dp
        call refresh_volume_mesh_view(s,preserve_cell=.true.)
        error=''
    end procedure set_workspace_mesh_section

    module procedure set_workspace_mesh_view
        integer :: direction,extent,index

        error='Choose the exterior, lattice, or layer volume mesh view.'
        if(mode<mesh_exterior.or.mode>mesh_layer)return
        direction=s%mesh_layer_axis
        if(present(axis))direction=axis
        error='Choose a valid volume mesh direction U, V, or W.'
        if(direction<1.or.direction>3)return
        extent=volume_mesh_extent(s,direction)
        error='Show a current volume mesh before changing its presentation.'
        if(extent<1)return
        index=min(s%mesh_layer_index,extent)
        if(present(layer))then
            error='Choose a volume mesh layer between 1 and '//integer_text(extent)//'.'
            if(layer<1.or.layer>extent)return
            index=layer
        end if
        s%mesh_view=mode
        s%mesh_layer_axis=direction
        s%mesh_layer_index=index
        call refresh_volume_mesh_view(s)
        error=''
    end procedure set_workspace_mesh_view

    subroutine refresh_volume_mesh_view(s,preserve_cell)
        type(studio_t),intent(inout)::s
        logical,intent(in),optional :: preserve_cell
        character(:),allocatable::description
        integer::axis
        logical :: clear_cell
        character(*),parameter :: world_axes='XYZ'

        clear_cell=.true.
        if(present(preserve_cell))clear_cell=.not.preserve_cell
        if(clear_cell.and.s%selected_mesh>0)then
            if(s%document%meshes(s%selected_mesh)%kind==3)s%mesh_cell=0
        end if
        if(clear_cell)then
            call update_mesh_ribbon(s)
        else
            call sync_volume_mesh_controls(s)
        end if
        description=trim(mesh_presentation_names(s%mesh_view))
        if(s%mesh_view==mesh_layer)then
            axis=s%mesh_layer_axis
            description=parameter_directions(axis:axis)//' layer '//integer_text(s%widgets%ribbon%volume_mesh_layer%value())
        end if
        axis=s%mesh_section_axis
        if(axis>0)then
            description=description//' · Section '//world_axes(axis:axis)//' · '// &
                integer_text(nint(100._dp*s%mesh_section_fraction))//'%'
        end if
        call s%status%showMessage('Volume mesh · '//description,3500)
        call s%surface%update()
    end subroutine refresh_volume_mesh_view

    module procedure has_visible_meshes
        integer::id,rank
        shown=.false.
        rank=1
        if(present(minimum_rank))rank=minimum_rank
        do id=1,studio_mesh_count(s)
            if(s%document%meshes(id)%kind<rank)cycle
            if(.not.scene_mesh_drawn(s,id))cycle
            shown=.true.
            return
        end do
    end procedure has_visible_meshes

    module procedure sync_snap_controls
        integer::active,i

        active=1
        if(s%snap_enabled)active=merge(3,2,s%snap_closest)
        do i=1,size(s%widgets%snap_mode_actions)
            call s%widgets%snap_mode_actions(i)%setChecked(i==active)
        end do
        call s%widgets%viewport_snap_action%setChecked(s%snap_enabled)
        if(s%snap_enabled)then
            call set_status_hint(s%widgets%viewport_snap,'Snap on · '//s%widgets%snap_mode_actions(active)%text()//' · Click to turn off')
        else
            call set_status_hint(s%widgets%viewport_snap,'Snap off · Click to turn on · Hold Shift to snap temporarily')
        end if
    end procedure sync_snap_controls

    module procedure sync_transform_controls
        integer::i
        logical::visible

        visible=s%transform_mode>0.and.transform_gizmo_available(s)
        call s%widgets%transform_strip%setVisible(visible)
        call s%widgets%transform_snap_button%setChecked(s%snap_enabled)
        do i=1,3
            call s%widgets%ribbon%transform_tools(i)%setChecked(visible.and.s%transform_mode==i)
        end do
    end procedure sync_transform_controls

    module procedure reset_transform_values
        integer::i
        real(dp),parameter::steps(3)=[.1_dp,15._dp,.1_dp]
        character(len=11),parameter::units(3)=[character(len=11)::'World units','Degrees','Scale factor']

        if(s%transform_mode<1.or.s%transform_mode>3)return
        if(s%transform_mode==1)then
            call s%widgets%transform_caption%setText('Move · '//project_unit_status(s))
        else
            call s%widgets%transform_caption%setText(trim(transform_captions(s%transform_mode))//' · '//trim(units(s%transform_mode)))
        end if
        do i=1,3
            call s%widgets%transform_inputs(i)%setRange(-huge(1._dp),huge(1._dp))
            call s%widgets%transform_inputs(i)%setSingleStep(steps(s%transform_mode))
            call s%widgets%transform_inputs(i)%setValue(merge(1._dp,0._dp,s%transform_mode==3))
            call set_status_hint(s%widgets%transform_inputs(i),trim(units(s%transform_mode)))
        end do
    end procedure reset_transform_values

    module procedure activate_transform
        real(dp)::bounds(4),x,y,center(2)

        if(s%transform_mode==mode)then
            call on_transform_exit(s)
            return
        end if
        if(.not.transform_gizmo_available(s))then
            call on_transform_exit(s)
            call run_geometry_transform(s,mode-1)
            return
        end if
        call cancel_area_drag(s)
        s%area_selection_tool=area_click
        s%surface_tool=1
        s%mesh_pick=.false.
        s%transform_mode=mode
        call reset_transform_values(s)
        call sync_transform_controls(s)
        call sync_selection_filter(s)
        call s%surface%setFocus()
        call live_viewport_point(s,transform_gizmo_pivot(s),x,y)
        bounds=workspace_view_bounds(s)
        if(.not.workspace_view_hit(s,x,y).or.x<bounds(1)+55._dp.or.x>bounds(3)-55._dp.or. &
            y<bounds(2)+230._dp.or.y>bounds(4)-55._dp)then
            if(min(x,y)<0._dp)then
                call fit_around_panels(s)
            else
                center=.5_dp*(bounds(:2)+bounds(3:))
                center(2)=max(center(2),bounds(2)+min(234._dp,.77_dp*(bounds(4)-bounds(2))))
                s%pan_x=s%pan_x+center(1)-x
                s%pan_y=s%pan_y+center(2)-y
            end if
        end if
        call s%surface%update()
        call s%status%showMessage(trim(transform_captions(mode))//' · Drag a handle or enter exact values',4000)
    end procedure activate_transform

    module procedure on_transform_exit

        s%transform_mode=0
        call sync_transform_controls(s)
        call sync_selection_filter(s)
        call s%surface%setFocus()
        call s%surface%update()
    end procedure on_transform_exit

    module procedure apply_transform_values
        type(affine_transform_t)::transform
        integer::i
        real(dp)::values(3)

        if(s%transform_mode==0.or..not.transform_gizmo_available(s))return
        values=[(s%widgets%transform_inputs(i)%value(),i=1,3)]
        transform%pivot=transform_gizmo_pivot(s)
        select case(s%transform_mode)
        case(1)
            transform%offset=model_length(values,s%document%model_unit,active_length_unit(s))
        case(2)
            transform%linear=rotation_matrix(values)
        case(3)
            do i=1,3
                transform%linear(i,i)=values(i)
            end do
        end select
        call apply_geometry_transform(s,transform,fit_after=.false.)
    end procedure apply_transform_values

    module procedure on_transform_drag
        type(affine_transform_t)::transform

        transform=affine_transform_t(linear,offset,pivot)
        call apply_geometry_transform(s,transform,fit_after=.false.)
    end procedure on_transform_drag

    module procedure on_assistant_geometry_changed

        if(kind==s%entity_kind.and.variant==active_patch_id(kind,catalog=s%document%patches))then
            s%live_scene%valid=.false.
            call refresh_inspector(s)
            call refresh_knot_table(s)
            call repaint_all(s)
        else
            call refresh_studio_scene(s%document%scene,catalog=s%document%patches)
            call sync_project_meshes(s)
            call sync_patch_visibility_button(s)
            call refresh_multipatch_card(s)
            call refresh_shared_controls(s)
            call refresh_edit_actions(s)
            if(s%selected_multipatch>0)call on_control_selection(s)
            call sync_project_title(s)
            call s%surface%update()
        end if
    end procedure on_assistant_geometry_changed

    module procedure on_assistant_model_created
        type(QTabWidget) :: tabs

        call cancel_transform_drag(s)
        s%transform_mode=0
        s%mesh_only=.false.
        s%mesh_pick=.false.
        if(mesh_id>0)then
            call rebuild_project_tree(s)
            s%mesh_only=.true.
            call select_project_mesh(s,mesh_id)
            call s%workspace%findChild('studioRibbonTabs',tabs)
            call tabs%setCurrentIndex(mesh_ribbon_tab)
        else if(group>0)then
            call show_scene_multipatch(s,group)
        else
            call rebuild_project_tree(s)
            call select_entity(s,kind,variant)
        end if
        call fit_around_panels(s)
        call repaint_all(s)
    end procedure on_assistant_model_created

    module procedure available_assistant_tools
        integer :: objects(3),kind,i
        integer, allocatable :: mesh_ids(:)
        logical :: one_selected
        type(selection_target_t) :: focused

        do kind=1,3
            objects(kind)=count([(document_patch(s,kind,i),i=1,single_patch_count(kind,catalog=s%document%patches))])
        end do
        context%has_curve=objects(1)>0
        context%has_surface=objects(2)>0
        context%has_volume=objects(3)>0
        context%has_visible_geometry=scene_drawn_patch_count(s)>0
        context%has_multipatch=s%selected_multipatch>0
        context%can_create_multipatch=any(objects>=2)
        focused=primary_workspace_target(s)
        one_selected=.not.allocated(s%selection%items).or.s%selection%count()==1
        context%has_selection=focused%category>0.and.(.not.allocated(s%selection%items).or.s%selection%count()>0)
        context%can_copy=context%has_selection.and.allocated(s%document%project_directory)
        context%can_paste=project_clipboard_in_scope(s)
        context%has_isolation=s%isolation%count()>0
        mesh_ids=selection_mesh_records(s)
        context%has_mesh=size(mesh_ids)>0
        context%has_selected_mesh=studio_mesh_valid(s,s%selected_mesh)
        context%has_visible_mesh=has_visible_meshes(s)
        context%has_shaded_mesh=has_visible_meshes(s,2)
        context%mesh_only=s%mesh_only
        context%can_edit=one_selected.and.s%live_scene%valid.and..not.context%has_multipatch.and.s%selected_mesh==0
        if(one_selected.and.context%has_multipatch)then
            context%can_edit_assembly_controls=allocated(s%document%multipatches%groups(s%selected_multipatch)%native)
        end if
        context%can_transform_gizmo=one_selected.and.transform_gizmo_available(s)
        context%transform_active=s%transform_mode>0
        context%domain_knots_active=s%panels(5)%isVisible().and.s%domain_edit_knots
        context%can_remove_active_knot=context%can_edit.and.context%domain_knots_active.and.can_remove_knot(s)
        context%can_remove_knot=context%can_edit.and.can_remove_knot(s,any_direction=.true.)
        context%can_reduce_degree=context%can_edit.and.can_reduce_degree(s)
        context%can_restructure_patch=context%can_edit.and.can_restructure_patch(s)
        context%can_extract_patch=context%can_edit
        context%can_undo=document_undo_stamp(s%document)>0
        context%can_redo=document_can_redo(s%document)
        context%can_sweep=objects(1)>=2.or.(context%has_curve.and.context%has_surface)
        context%can_loft=any(objects(:2)>=2)
        if(.not.s%widgets%ribbon%view_fit_action%isNull())then
            context%can_fit_selection=s%widgets%ribbon%view_fit_action%isEnabled()
        end if
    end procedure available_assistant_tools

    module procedure refresh_rotated_scene
        integer :: group

        call cancel_transform_drag(s)
        s%transform_mode=0
        s%control_dragging=.false.
        s%control_snapped=.false.
        s%viewport_feature=0
        s%live_scene%valid=.false.
        s%mesh_cell=0
        do group=1,count_multipatches(s%document%multipatches)
            call refresh_multipatch(s%document%multipatches,group,.false.,catalog=s%document%patches)
        end do
        call refresh_studio_scene(s%document%scene,catalog=s%document%patches)
        if(.not.has_visible_meshes(s))then
            s%mesh_pick=.false.
            s%mesh_only=.false.
        end if
        call refresh_inspector(s)
        call refresh_knot_table(s)
        call sync_controls(s)
        call fit_around_panels(s)
        call repaint_all(s)
    end procedure refresh_rotated_scene

    module procedure run_multipatch_editor
        type(multipatch_draft)::draft
        character(:),allocatable::error
        character(96),allocatable::reserved_names(:)
        integer,allocatable :: members(:)
        integer :: current,kind,i,member,first_variant,rank,selected_counts(3)
        logical :: eligible(3)

        if(id<0.or.id>count_multipatches(s%document%multipatches))return
        current=id
        kind=s%entity_kind
        members=[integer::]
        allocate(reserved_names(count_multipatches(s%document%multipatches)))
        reserved_names=''
        do i=1,size(reserved_names)
            if(i/=id)reserved_names(i)=s%document%multipatches%groups(i)%definition%name
        end do
        if(id>0)then
            draft=s%document%multipatches%groups(id)%definition
            kind=draft%kind
        else
            eligible=eligible_assembly_ranks(s)
            if(.not.any(eligible))then
                call s%status%showMessage('Add two patches of the same type to create a multipatch.',4500)
                return
            end if
            do rank=1,3
                members=initial_assembly_members(s,rank)
                selected_counts(rank)=size(members)
            end do
            rank=maxloc(selected_counts,dim=1,mask=eligible)
            if(.not.eligible(kind).or.(selected_counts(rank)>=2.and.selected_counts(rank)>selected_counts(kind)))kind=rank
            members=initial_assembly_members(s,kind)
        end if
        member=0
        if(size(members)>0)member=members(1)
        first_variant=merge(1,3,s%document%show_templates)
        if(.not.edit_studio_multipatch(s%window,kind,member,draft,reserved_names,first_variant, &
            catalog=s%document%patches,initial_members=members))return
        call save_multipatch(s%document%multipatches,current,draft,error,catalog=s%document%patches)
        if(len(error)>0)then
            call s%status%showMessage(error,5000,warning=.true.)
            return
        end if
        call rebuild_project_tree(s)
        call select_multipatch_group(s,current,reveal=.false.)
        call invoke_command(s,fit_selection_command)
        call s%status%showMessage('Multipatch saved · '//draft%name,3500)
    end procedure run_multipatch_editor

    module procedure run_geometry_transform
        type(affine_transform_t)::transform
        integer::kind,id,i,current_kind
        integer,allocatable::variants(:),current_variants(:)
        integer(int64),allocatable::revisions(:)
        character(:),allocatable::name,current_name
        logical::accepted

        call selected_geometry(s,kind,variants,name)
        if(size(variants)==0)return
        id=s%selected_multipatch
        revisions=[(single_patch_revision(kind,variants(i),catalog=s%document%patches),i=1,size(variants))]
        if(s%transform_mode>0.and.transform_gizmo_available(s))then
            accepted=edit_studio_transform(s%window,name,transform,operation,transform_gizmo_pivot(s), &
                s%document%model_unit,active_length_unit(s))
        else
            accepted=edit_studio_transform(s%window,name,transform,operation, &
                model_unit=s%document%model_unit,display_unit=active_length_unit(s))
        end if
        if(.not.accepted)return
        call selected_geometry(s,current_kind,current_variants,current_name)
        if(id==s%selected_multipatch.and.kind==current_kind.and.size(variants)==size(current_variants))then
            if(all(variants==current_variants).and. &
                all(revisions==[(single_patch_revision(kind,variants(i),catalog=s%document%patches),i=1,size(variants))]))then
                call apply_geometry_transform(s,transform,fit_after=s%transform_mode==0)
                return
            end if
        end if
        call s%status%showMessage('The transform target changed. Reopen the editor for the current selection.',5500)
    end procedure run_geometry_transform

    module procedure apply_geometry_transform
        integer::kind,id
        integer,allocatable::variants(:)
        character(:),allocatable::name,error
        logical::changed,show_target

        call selected_geometry(s,kind,variants,name)
        id=s%selected_multipatch
        call commit_document_transform(s%document,kind,variants,id,transform,changed,error)
        if(present(changed_out))changed_out=changed
        if(present(error_out))error_out=error
        if(len(error)>0)then
            if(size(variants)>0)call s%status%showMessage(error,5000,warning=.true.)
            return
        end if
        if(.not.changed)return
        s%live_scene%valid=.false.
        s%viewport_feature=0
        show_target=.false.
        if(present(reveal))show_target=reveal
        if(id>0)then
            call select_multipatch_group(s,id,reveal=show_target)
        else
            call select_patch(s,kind,variants(1),reveal=show_target)
        end if
        if(fit_after)call invoke_command(s,fit_selection_command)
        call reset_transform_values(s)
        call s%status%showMessage('Transformed '//name,3500)
    end procedure apply_geometry_transform

    module procedure set_workspace_material
        integer :: kind
        integer,allocatable :: variants(:)
        character(:),allocatable :: name
        logical :: saved

        changed=.false.
        call selected_geometry(s,kind,variants,name)
        error='Select a patch or assembly to set its appearance.'
        if(size(variants)==0)return
        changed=commit_studio_material(s%document%scene,kind,variants,rgb,opacity,catalog=s%document%patches,saved=saved)
        if(.not.saved)then
            error='Could not record the appearance edit in the history.'
            call s%status%showMessage(error,5000,warning=.true.)
            return
        end if
        if(changed)then
            call sync_controls(s)
            call refresh_edit_actions(s)
            call sync_project_title(s)
            call repaint_all(s)
            call s%status%showMessage('Updated appearance · '//name,3500)
        end if
        error=''
    end procedure set_workspace_material

    module procedure select_workspace_mesh_cell
        integer :: i,count,counts(3),stride(3),corner(3),axis,j,node

        part=0
        cell=number
        error='Select a current, visible mesh to pick its cells.'
        if(.not.studio_mesh_valid(s,id))return
        if(.not.s%document%meshes(id)%visible)return
        error='That cell number is outside the mesh.'
        if(number<1.or.number>s%document%meshes(id)%elements)return
        do i=1,size(s%document%meshes(id)%parts)
            count=size(s%document%meshes(id)%parts(i)%cells,1)
            if(cell<=count)then
                part=i
                if(s%document%meshes(id)%kind==3)then
                    counts=s%document%meshes(id)%parts(part)%samples
                    stride=[1,counts(1),product(counts(:2))]
                    node=s%document%meshes(id)%parts(part)%cells(cell,1)-1
                    corner=[(mod(node/stride(j),counts(j))+1,j=1,3)]
                    axis=s%mesh_layer_axis
                    if(s%mesh_view==mesh_layer.or.all(corner>1.and.corner<counts-1))then
                        s%mesh_view=mesh_layer
                        s%mesh_layer_index=corner(axis)
                    end if
                end if
                s%mesh_pick=.true.
                call apply_workspace_selection(s,selection_target_t(selection_mesh,id,part=part,item=cell),selection_replace)
                error=''
                return
            end if
            cell=cell-count
        end do
    end procedure select_workspace_mesh_cell

    module procedure change_material
        integer::kind,rgb(3),opacity
        integer,allocatable::variants(:)
        character(:),allocatable::name
        type(QColor)::color
        logical :: changed,saved

        call selected_geometry(s,kind,variants,name)
        if(size(variants)==0)return
        call get_studio_material(s%document%scene,kind,variants(1),rgb,opacity)
        if(change_color)then
            color=QColorDialog_getColor(QColor(rgb(1),rgb(2),rgb(3)),title='Object color',parent=s%window)
            if(.not.color%isValid())return
            changed=commit_studio_material(s%document%scene,kind,variants,rgb=[color%r,color%g,color%b],catalog=s%document%patches,saved=saved)
        else
            opacity=merge(255,112,opacity<255)
            changed=commit_studio_material(s%document%scene,kind,variants,opacity=opacity,catalog=s%document%patches,saved=saved)
        end if
        if(.not.saved)then
            call s%status%showMessage('Could not record the appearance edit in the history.',5000,warning=.true.)
            return
        end if
        if(.not.changed)return
        call refresh_edit_actions(s)
        call sync_project_title(s)
        call s%status%showMessage('Updated appearance · '//name,3500)
    end procedure change_material

    module procedure generate_workspace_mesh
        type(mesh_source_t) :: source
        type(mesh_publication_t) :: published
        character(:),allocatable :: error,caption
        type(mesh_sizing_t) :: sizing
        integer,allocatable :: cells(:),member_counts(:,:)
        integer :: i,slot
        logical :: connected,default_sizing

        if(present(mesh_id))mesh_id=0
        if(present(failure))failure=''
        if(present(source_override))then
            source=source_override
        else
            source=mesh_source_selection(s)
        end if
        if(source%assembly_id>0)then
            call refresh_multipatch(s%document%multipatches,source%assembly_id,.false.,catalog=s%document%patches)
        end if
        connected=.false.
        if(source%mesh_id>0.and.source%mesh_id<=studio_mesh_count(s))then
            connected=allocated(s%document%meshes(source%mesh_id)%connected_cells)
        end if
        default_sizing=.not.(present(requested_cells).or.present(patch_cells).or.present(target_size))
        if(present(prepared_sizing))then
            sizing=prepared_sizing
            connected=sizing%conforming
        else if(default_sizing)then
            if(allocated(s%selection%items))then
                if(s%selection%count()/=1)then
                    error='Select one patch, multipatch, or mesh to generate.'
                    if(present(failure))failure=error
                    call show_mesh_diagnostic(s,error)
                    return
                end if
            end if
            sizing=mesh_sizing(s,source)
            connected=sizing%conforming
        end if
        if(present(conforming))connected=conforming
        if(present(target_size))sizing%target_size=target_size
        if(present(synchronize))sizing%synchronize=synchronize
        if(present(target_size).or.sizing%target_size>0._dp)then
            error=''
            if(.not.allocated(sizing%preview_revisions))then
                call estimate_document_mesh_samples(s%document,source,sizing%target_size,sizing%synchronize,sizing%cells,error)
            else if(.not.mesh_preview_current(s,source,sizing%preview_revisions))then
                error='The source changed after sizing. Open Generate mesh and try again.'
            end if
            if(len(error)>0)then
                if(present(failure))failure=error
                call show_mesh_diagnostic(s,error)
                return
            end if
            sizing%per_patch=source%assembly_id>0
        end if
        if(present(requested_cells))then
            cells=requested_cells
        else if(present(patch_cells))then
            allocate(cells(0))
        else if(allocated(sizing%cells))then
            if(size(sizing%cells,2)>0)cells=sizing%cells(:,1)
            if(sizing%per_patch)member_counts=sizing%cells
        end if
        if(.not.allocated(cells))allocate(cells(0))
        if(present(patch_cells))member_counts=patch_cells
        call generate_document_mesh(s%document,source,cells,connected,s%selected_mesh,s%mesh_only,published,error,member_counts, &
            sizing%target_size,sizing%synchronize)
        if(len(error)>0)then
            if(present(failure))failure=error
            select case(published%failure)
            case(mesh_source_missing)
                if(s%selected_mesh>0)call s%status%showMessage('Mesh source is unavailable',5000)
            case(mesh_build_failed)
                call show_mesh_diagnostic(s,error)
            case default
                if(.not.present(failure))call show_mesh_diagnostic(s,error)
            end select
            return
        end if

        caption='Generated '
        if(published%remeshed)then
            caption='Remeshed · '
            s%mesh_cell=0
            call reconcile_remeshed_selection(s,source%mesh_id,published)
            do i=size(published%removed_ids),1,-1
                call adjust_mesh_selection_after_removal(s,published%removed_ids(i))
            end do
            if(size(published%removed_ids)>0)call rebuild_project_tree(s)
        end if
        call record_studio_diagnostic(s,'')
        call refresh_studio_diagnostics(s)
        call clamp_mesh_layer(s)
        call sync_project_meshes(s)
        call s%widgets%project_browser%mesh_branch%setExpanded(.true.)
        call sync_patch_visibility_button(s)
        slot=mesh_draft_index(s,source)
        if(slot>0)s%widgets%mesh_drafts(slot)=mesh_sizing_t()
        if(.not.allocated(member_counts).and.sizing%target_size<=0._dp)then
            s%widgets%mesh_divisions(:source%kind,source%kind)=cells
        end if
        call update_mesh_ribbon(s)
        call refresh_edit_actions(s)
        call sync_project_title(s)
        if(present(mesh_id))mesh_id=published%id
        call refresh_workspace_status(s)
        call s%status%showMessage(caption//integer_text(published%cells)//' cells · Scene → Meshes',6000)
        call repaint_all(s)
    end procedure generate_workspace_mesh

    subroutine reconcile_remeshed_selection(s,previous_id,published)
        type(studio_t), intent(inout) :: s
        integer, intent(in) :: previous_id
        type(mesh_publication_t), intent(in) :: published
        logical, allocatable :: keep(:)
        integer :: i,latest

        if(s%selection%count()==0)return
        latest=0
        associate(items=>s%selection%items)
            allocate(keep(size(items)),source=.true.)
            do i=1,size(items)
                if(items(i)%category/=selection_mesh)cycle
                if(items(i)%id==previous_id.or.any(published%removed_ids==items(i)%id))then
                    items(i)=selection_target_t(selection_mesh,published%id)
                    if(latest>0)keep(latest)=.false.
                    latest=i
                else
                    items(i)%id=items(i)%id-count(published%removed_ids<items(i)%id)
                end if
            end do
        end associate
        if(.not.all(keep))s%selection%items=pack(s%selection%items,keep)
    end subroutine reconcile_remeshed_selection

    module procedure run_patch_mesh_dialog
        type(mesh_source_t) :: source,current
        type(mesh_sizing_t) :: settings
        character(:), allocatable :: labels(:),operation
        integer :: group,i,rank,label_width
        integer(int64) :: mesh_generation
        integer(int64), allocatable :: preview_revisions(:)
        logical :: can_join,start_with_target,build_now

        source=mesh_source_selection(s)
        if(source%assembly_id>0)then
            call refresh_multipatch(s%document%multipatches,source%assembly_id,.false.,catalog=s%document%patches)
        end if
        if(.not.mesh_source_available(s,source))return
        if(allocated(s%selection%items))then
            if(s%selection%count()/=1)then
                call show_mesh_diagnostic(s,'Select one patch, multipatch, or mesh to edit the sizing.')
                return
            end if
        end if
        settings=mesh_sizing(s,source)
        build_now=.false.
        if(present(generate))build_now=generate
        start_with_target=by_target_size
        if(build_now)then
            operation='Generate'
            if(source%mesh_id>0)operation='Remesh'
            start_with_target=start_with_target.or.settings%target_size>0._dp
        end if
        group=source%assembly_id
        rank=source%kind
        mesh_generation=s%document%mesh_generation
        preview_revisions=mesh_source_revisions(s,source)
        label_width=maxval([(len(s%document%patches%patch_name(rank,source%members(i))),i=1,size(source%members))])
        allocate(character(max(1,label_width)) :: labels(size(source%members)))
        do i=1,size(source%members)
            labels(i)=s%document%patches%patch_name(rank,source%members(i))
        end do

        can_join=.false.
        if(rank==3.and.group>0)then
            associate(definition=>s%document%multipatches%groups(group)%definition)
                if(allocated(definition%joins))can_join=size(definition%joins)>0
            end associate
        end if
        if(.not.edit_mesh_patch_counts(s%window,source%name,labels,settings%cells, &
            mesh_axis_cell_limit(rank,size(labels)),mesh_axis_cell_limit(rank,1),group>0,settings%per_patch, &
            can_join,settings%conforming,s%document,source,settings%target_size,settings%synchronize,active_length_unit(s), &
            start_with_target,operation))return

        current=mesh_source_selection(s)
        if(s%document%document_revision/=settings%document_revision.or.s%document%mesh_generation/=mesh_generation.or. &
            current%kind/=rank.or.current%assembly_id/=group.or.current%mesh_id/=source%mesh_id.or. &
            size(current%members)/=size(source%members))then
            call show_mesh_diagnostic(s,'The mesh source changed. Select it and reopen the mesh dialog.')
            return
        end if
        if(any(current%members/=source%members))then
            call show_mesh_diagnostic(s,'The multipatch changed. Select it and reopen the mesh dialog.')
            return
        end if
        if(settings%target_size>0._dp)then
            if(.not.mesh_preview_current(s,source,preview_revisions))then
                call show_mesh_diagnostic(s,'The source changed after sizing. Reopen the mesh dialog.')
                return
            end if
            settings%preview_revisions=preview_revisions
        else if(allocated(settings%preview_revisions))then
            deallocate(settings%preview_revisions)
        end if
        if(build_now)then
            call generate_workspace_mesh(s,source_override=source,prepared_sizing=settings)
            return
        end if
        call set_mesh_draft(s,source,settings)
        call update_mesh_ribbon(s)
        if(settings%target_size>0._dp)then
            call s%status%showMessage('Target-size counts ready to generate',3500)
        else
            call s%status%showMessage('Cell counts ready to generate',3500)
        end if
    end procedure run_patch_mesh_dialog

    subroutine show_mesh_diagnostic(s,message)
        type(studio_t),intent(inout) :: s
        character(*),intent(in) :: message

        call record_studio_diagnostic(s,message)
        call refresh_studio_diagnostics(s)
        call s%status%showMessage(message,7000)
    end subroutine show_mesh_diagnostic

    module procedure clamp_mesh_layer


        s%mesh_layer_index=min(s%mesh_layer_index,max(1,volume_mesh_extent(s,s%mesh_layer_axis)))
    end procedure clamp_mesh_layer

    module procedure save_selected_mesh
        character(:),allocatable::path,error
        integer::id
        integer(int64)::tree_stamp

        id=mesh_context_record(s)
        if(.not.studio_mesh_valid(s,id))return
        tree_stamp=s%widgets%project_browser%stamp
        path=choose_export_path(s,'Save generated mesh','mesh','VTK mesh (*.vtk)','.vtk')
        if(len(path)==0)return
        if(tree_stamp/=s%widgets%project_browser%stamp.or..not.studio_mesh_valid(s,id))then
            call s%status%showMessage('The model changed. Generate a new mesh before exporting.',6000)
            return
        end if
        call export_studio_mesh(s%document%meshes(id),path,error)
        if(len(error)>0)then
            call s%status%showMessage(error,7000,warning=.true.)
        else
            call s%status%showMessage('Exported mesh · '//path,6000)
        end if
    end procedure save_selected_mesh

    module procedure toggle_selected_mesh_visibility
        integer, allocatable :: ids(:)
        logical :: visible,changed

        ids=selection_mesh_records(s)
        if(size(ids)==0)return
        visible=.not.all(s%document%meshes(ids)%visible)
        call set_meshes_visibility(s,ids,visible,changed)
    end procedure toggle_selected_mesh_visibility

    module procedure set_mesh_visibility

        call set_meshes_visibility(s,[id],visible,changed,error)
    end procedure set_mesh_visibility

    subroutine set_meshes_visibility(s,ids,visible,changed,error)
        type(studio_t), intent(inout) :: s
        integer, intent(in) :: ids(:)
        logical, intent(in) :: visible
        logical, intent(out) :: changed
        character(:), allocatable, intent(out), optional :: error
        character(:), allocatable :: message

        call commit_document_visibility(s%document,visible,changed,message,meshes=ids)
        if(present(error))error=message
        if(len(message)>0)call s%status%showMessage(message,5000,warning=.true.)
        if(.not.changed)return
        if(.not.has_visible_meshes(s))s%mesh_only=.false.
        call sync_project_meshes(s)
        if(.not.has_visible_meshes(s))s%mesh_pick=.false.
        call sync_patch_visibility_button(s)
        call update_mesh_ribbon(s)
        call refresh_edit_actions(s)
        call sync_project_title(s)
        call repaint_all(s)
    end subroutine set_meshes_visibility

    module procedure on_mesh_cell_selected

        call forget_workspace_selection(s)
        if(id==0)then
            s%mesh_cell=0
            if(s%mesh_pick)call s%status%showMessage('Click a visible mesh cell',3500)
        else
            if(.not.studio_mesh_valid(s,id))return
            if(.not.s%document%meshes(id)%visible)return
            if(id/=s%selected_mesh)call select_project_mesh(s,id)
            s%mesh_cell=[part,cell,face]
        end if
        call update_mesh_ribbon(s)
        call refresh_workspace_status(s)
        call s%surface%update()
    end procedure on_mesh_cell_selected

    module procedure show_created_patch

        call on_transform_exit(s)
        call rebuild_project_tree(s)
        call select_entity(s,kind)
        call fit_around_panels(s)
        call s%surface%update()
        call s%status%showMessage('Created '//single_patch_name(kind,active_patch_id(kind,catalog=s%document%patches),catalog=s%document%patches),5000)
    end procedure show_created_patch

    module procedure change_snap_mode
        character(:),allocatable::message

        select case(command)
        case(snap_command)
            s%snap_enabled=.not.s%snap_enabled
        case(snap_off_command)
            s%snap_enabled=.false.
        case(snap_vertex_command,snap_closest_command)
            s%snap_enabled=.true.
            s%snap_closest=command==snap_closest_command
        end select
        if(.not.s%snap_enabled)then
            message='Snap off · Hold Shift to snap temporarily'
        else if(s%transform_mode>0.and.command==snap_command)then
            message='Snap on · Move 0.1 model units · Rotate 15° · Scale 0.1×'
        else if(s%snap_closest)then
            message='Closest-point snap · Drag a control point near another visible patch'
        else
            message='Vertex snap · Drag a control point near a vertex'
        end if
        call sync_controls(s)
        call s%status%showMessage(message,5000)
        call s%surface%update()
    end procedure change_snap_mode
end submodule studio_ui_model_tools
