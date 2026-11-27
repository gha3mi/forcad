module studio_canvas_input
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use studio_state, only: studio_t, dp, repaint_all, area_click, area_box, area_lasso, clear_control_dependency_preview
    use studio_presentation_text, only: integer_text
    use studio_selection, only: selection_target_t, selection_replace, selection_toggle, &
        selection_control, selection_mesh, pick_multipatch
    use studio_workspace_selection, only: pointer_selection_mode
    use studio_domain_projection, only: volume_domain_projection, parameter_track, domain_axis_track, live_volume_projection, &
        ruler_at, volume_probe_at, planar_probe_at, track_coordinates, cube_pitch, min_ruler_length, domain_ruler_gap
    use studio_viewport_camera, only: zoom_live_viewport,reset_live_viewport_zoom,orbit_live_viewport,viewport_zoom_active, &
        frame_workspace_scene, workspace_view_hit, workspace_view_bounds, set_viewport_orientation
    use studio_viewport_scene, only: pick_live_viewport, begin_live_control_drag, update_live_control_drag, pick_scene_patch, &
        pick_automatic, pick_control_points, pick_selected_control, viewport_patch_target, multipatch_control_scope, &
        pick_shared_control, viewport_control_position
    use studio_viewport_transform, only: begin_transform_drag, update_transform_drag, finish_transform_drag, &
        cancel_transform_drag, transform_dragging, transform_gizmo_hint
    use studio_live_viewport, only: pick_mesh_cell
    use studio_viewport_selection, only: collect_viewport_region
    use studio_knot_interaction, only: knot_edit_in_progress, handle_domain_knot_press, handle_domain_knot_move, handle_domain_knot_release
    use studio_selection_status, only: selection_status_text
    use qt, only: QMouseEvent, QKeyEvent, QWheelEvent, QTimer, QPoint, Qt_LeftButton, Qt_MiddleButton, Qt_RightButton, &
        Qt_NoButton, Qt_ShiftModifier, Qt_Key_Escape, Qt_StrongFocus, Qt_ClosedHandCursor, Qt_CrossCursor
    implicit none(type, external)
    private

    public :: connect_canvas_input, zoom_domain_view, reset_domain_view, cancel_area_drag
    public :: min_domain_zoom, max_domain_zoom, domain_pitch_limit
    public :: control_move_handler, control_selection_handler, selection_handler, transform_handler, viewport_context_handler
    public :: region_selection_handler, change_shared_selection, picking_shared_controls

    integer, parameter :: volume_orbit=-1, domain_pan=-2, domain_interior=4
    real(dp), parameter :: min_domain_zoom=.65_dp, max_domain_zoom=24._dp, domain_zoom_step=1.16_dp
    real(dp), parameter :: orbit_drag_threshold=4._dp
    real(dp), parameter :: domain_pitch_limit=1.30_dp

    abstract interface
        subroutine control_move_handler(s,point,position)
            import studio_t,dp
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            integer, intent(in) :: point
            real(dp), intent(in) :: position(3)
        end subroutine control_move_handler

        subroutine control_selection_handler(s)
            import studio_t
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
        end subroutine control_selection_handler

        subroutine selection_handler(s,item,mode)
            import studio_t,selection_target_t
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            type(selection_target_t), intent(in) :: item
            integer, intent(in) :: mode
        end subroutine selection_handler

        subroutine region_selection_handler(s,items,mode)
            import studio_t,selection_target_t
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            type(selection_target_t), intent(in) :: items(:)
            integer, intent(in) :: mode
        end subroutine region_selection_handler

        subroutine transform_handler(s,linear,offset,pivot)
            import studio_t,dp
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            real(dp), intent(in) :: linear(3,3),offset(3),pivot(3)
        end subroutine transform_handler

        subroutine viewport_context_handler(s,x,y)
            import studio_t,dp
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            real(dp), intent(in) :: x,y
        end subroutine viewport_context_handler

    end interface

    procedure(control_move_handler), pointer :: commit_control_move=>null()
    procedure(control_selection_handler), pointer :: notify_control_selection=>null(), open_selected_control=>null()
    procedure(selection_handler), pointer :: notify_selection=>null()
    procedure(region_selection_handler), pointer :: notify_region_selection=>null()
    procedure(transform_handler), pointer :: commit_transform=>null()
    procedure(viewport_context_handler), pointer :: open_viewport_context=>null()
    procedure(control_selection_handler), pointer :: notify_transform_exit=>null()

contains

    subroutine connect_canvas_input(s,control_move,control_selected,control_open,selection_changed, &
                                transform_committed,viewport_context,transform_exited,region_selected)
        type(studio_t), intent(inout), target :: s
        procedure(control_move_handler), optional :: control_move
        procedure(control_selection_handler), optional :: control_selected,control_open
        procedure(selection_handler), optional :: selection_changed
        procedure(region_selection_handler), optional :: region_selected
        procedure(transform_handler), optional :: transform_committed
        procedure(viewport_context_handler), optional :: viewport_context
        procedure(control_selection_handler), optional :: transform_exited

        nullify(commit_control_move,notify_control_selection,open_selected_control,notify_selection,notify_region_selection)
        nullify(commit_transform,open_viewport_context,notify_transform_exit)
        if(present(control_move))commit_control_move=>control_move
        if(present(control_selected))notify_control_selection=>control_selected
        if(present(control_open))open_selected_control=>control_open
        if(present(selection_changed))notify_selection=>selection_changed
        if(present(region_selected))notify_region_selection=>region_selected
        if(present(transform_committed))commit_transform=>transform_committed
        if(present(viewport_context))open_viewport_context=>viewport_context
        if(present(transform_exited))notify_transform_exit=>transform_exited
        call cancel_transform_drag(s)
        s%canvas_gesture%domain_pointer=0
        s%canvas_gesture%domain_button=Qt_NoButton
        call s%surface%mousePressEvent%connect(surface_press,s)
        call s%surface%mouseMoveEvent%connect(surface_move,s)
        call s%surface%leaveEvent%connect(surface_leave,s)
        call s%surface%mouseReleaseEvent%connect(surface_release,s)
        call s%surface%mouseDoubleClickEvent%connect(surface_reset,s)
        call s%surface%wheelEvent%connect(surface_wheel,s)
        call s%surface%keyPressEvent%connect(surface_key,s)
        if(s%viewport_hover_timer%isNull())then
            s%viewport_hover_timer=QTimer(s%surface)
            call s%viewport_hover_timer%setSingleShot(.true.)
            call s%viewport_hover_timer%timeout%connect(settle_viewport_hover,s)
        end if
        call s%surface%setFocusPolicy(Qt_StrongFocus)
        call s%domain%mousePressEvent%connect(domain_press,s)
        call s%domain%mouseMoveEvent%connect(domain_move,s)
        call s%domain%mouseReleaseEvent%connect(domain_release,s)
        call s%domain%mouseDoubleClickEvent%connect(domain_reset,s)
        call s%domain%wheelEvent%connect(domain_wheel,s)
        call s%surface%setAccessibleDescription('Click to select, Ctrl to toggle, Shift to add, Esc to clear. '// &
            'Right-drag or Shift+middle-drag to orbit, middle-drag to pan, and scroll to zoom.')
        call s%domain%setAccessibleDescription('Drag to probe, scroll to zoom, and middle-drag to pan. '// &
            'Right-drag or Shift+middle-drag to orbit a volume.')
    end subroutine connect_canvas_input

    logical function shifted_middle_button(event) result(shifted)
        type(QMouseEvent),intent(in) :: event

        shifted=event%button()==Qt_MiddleButton.and.iand(event%modifiers(),Qt_ShiftModifier)/=0
    end function shifted_middle_button

    subroutine surface_press(event,data)
        type(QMouseEvent), intent(in) :: event
        class(*), intent(inout), optional :: data
        logical :: hit
        integer :: kind,variant
        if(.not.any(event%button()==[Qt_LeftButton,Qt_MiddleButton,Qt_RightButton]))return
        select type(s=>data)
        type is(studio_t)
            if(s%dragging)return
            if(.not.workspace_view_hit(s,event%x(),event%y()))return
            call clear_viewport_hover(s)
            s%canvas_gesture%pointer_button=event%button()
            s%canvas_gesture%middle_orbit=shifted_middle_button(event)
            s%canvas_gesture%pointer_displaced=.false.
            s%canvas_gesture%probe_dragging=.false.
            s%dragging=.true.
            s%last_x=event%x()
            s%last_y=event%y()
            call s%surface%setFocus()
            if(event%button()==Qt_LeftButton.and.s%area_selection_tool/=area_click.and.s%transform_mode==0.and. &
                associated(notify_region_selection))then
                call begin_area_drag(s,event%x(),event%y())
                return
            end if
            if((s%mesh_pick.or.s%mesh_only).and.s%canvas_gesture%pointer_button==Qt_LeftButton)return
            if(s%transform_mode>0.and.s%canvas_gesture%pointer_button==Qt_LeftButton)then
                if(associated(commit_transform))then
                    hit=begin_transform_drag(s,event%x(),event%y())
                    if(hit)call s%status%setActivity(transform_gizmo_hint(s))
                end if
                return
            end if
            if(picking_controls(s).and.s%canvas_gesture%pointer_button==Qt_LeftButton.and. &
                pointer_selection_mode(event%modifiers())==selection_replace)then
                hit=select_control_at_pointer(s,event%x(),event%y(),selection_replace,drag=.true.)
                if(.not.hit.and.s%feature_pick_mode==pick_automatic.and..not.multipatch_control_scope(s))then
                    call pick_scene_patch(s,event%x(),event%y(),kind,variant)
                    s%canvas_gesture%probe_dragging=kind==s%entity_kind.and.variant==s%live_scene%variant
                    if(s%canvas_gesture%probe_dragging)call pick_live_viewport(s,event%x(),event%y())
                end if
                call s%surface%update()
            end if
        end select
    end subroutine surface_press

    subroutine surface_move(event,data)
        type(QMouseEvent), intent(in) :: event
        class(*), intent(inout), optional :: data
        real(dp) :: dx,dy
        select type(s=>data)
        type is(studio_t)
            if(.not.s%dragging)then
                if(viewport_zoom_active(s))return
                if(s%canvas_gesture%frame_time_ms>24._dp)then
                    s%canvas_gesture%hover_position=[event%x(),event%y()]
                    call s%viewport_hover_timer%start(42)
                else
                    call s%viewport_hover_timer%stop()
                    call update_viewport_hover(s,event%x(),event%y())
                end if
                return
            end if
            if(s%canvas_gesture%area_count>0)then
                call track_area_drag(s,event%x(),event%y(),.false.)
                return
            end if
            if(transform_dragging(s))then
                call update_transform_drag(s,event%x(),event%y(),iand(event%modifiers(),Qt_ShiftModifier)/=0)
                call s%status%setActivity(transform_gizmo_hint(s))
                return
            end if
            if(s%control_dragging)then
                call update_live_control_drag(s,event%x(),event%y(),iand(event%modifiers(),Qt_ShiftModifier)/=0)
                call s%status%setActivity(control_drag_activity(s))
                return
            end if
            dx=event%x()-s%last_x; dy=event%y()-s%last_y
            if(.not.s%canvas_gesture%pointer_displaced.and.hypot(dx,dy)<orbit_drag_threshold)return
            s%canvas_gesture%pointer_displaced=.true.
            if(s%canvas_gesture%probe_dragging)then
                call pick_live_viewport(s,event%x(),event%y())
            else if(s%canvas_gesture%pointer_button==Qt_MiddleButton.and..not.s%canvas_gesture%middle_orbit)then
                s%pan_x=s%pan_x+dx; s%pan_y=s%pan_y+dy
            else if(s%canvas_gesture%middle_orbit.or.s%canvas_gesture%pointer_button==Qt_RightButton.or. &
                (s%surface_tool==2.and.s%transform_mode==0.and..not.s%mesh_pick))then
                call orbit_live_viewport(s,dx,dy)
            end if
            s%last_x=event%x()
            s%last_y=event%y()
            call s%surface%update()
        end select
    end subroutine surface_move

    function control_drag_activity(s) result(message)
        type(studio_t), intent(in) :: s
        character(:), allocatable :: message
        integer :: n

        if(s%control_snapped)then
            message=trim(merge('Snapped to closest point','Snapped to vertex       ',s%snap_closest))
        else if(multipatch_control_scope(s))then
            message='Moving selected multipatch controls'
        else
            message='Moving control point'
        end if
        associate(source=>s%control_dependency_source,dependent=>s%control_dependency_preview)
            if(.not.source%attempted)return
            if(allocated(source%error))then
                if(len(source%error)>0)then
                    message=source%error
                    return
                end if
            end if
            if(.not.allocated(dependent%ids))return
            n=size(dependent%ids)
            if(n==0)then
                message=message//' · No dependent controls'
            else if(n==1)then
                message=message//' · 1 dependent control in amber'
            else
                message=message//' · '//integer_text(n)//' dependent controls in amber'
            end if
        end associate
    end function control_drag_activity

    subroutine surface_leave(data)
        class(*),intent(inout),optional :: data

        select type(s=>data)
        type is(studio_t)
            call clear_viewport_hover(s)
        end select
    end subroutine surface_leave

    subroutine clear_viewport_hover(s)
        type(studio_t),intent(inout),target :: s
        logical :: had_shared_hover

        if(.not.s%viewport_hover_timer%isNull())call s%viewport_hover_timer%stop()
        had_shared_hover=s%canvas_gesture%viewport_shared_hover>0
        if(.not.had_shared_hover.and.all(s%canvas_gesture%viewport_hover==0).and.all(s%canvas_gesture%viewport_mesh_hover==0))return
        s%canvas_gesture%viewport_hover=0
        s%canvas_gesture%viewport_mesh_hover=0
        s%canvas_gesture%viewport_shared_hover=0
        if(had_shared_hover.and.associated(notify_control_selection))call notify_control_selection(s)
        call s%surface%update()
    end subroutine clear_viewport_hover

    subroutine settle_viewport_hover(data)
        class(*),intent(inout),optional :: data

        select type(s=>data)
        type is(studio_t)
            if(s%dragging.or.viewport_zoom_active(s))return
            associate(position=>s%canvas_gesture%hover_position)
                call update_viewport_hover(s,position(1),position(2))
            end associate
        end select
    end subroutine settle_viewport_hover

    subroutine update_viewport_hover(s,x,y)
        type(studio_t),intent(inout),target :: s
        real(dp),intent(in) :: x,y
        integer :: target(4),mesh_target(3),face,shared_point
        logical :: shared_changed

        target=0
        mesh_target=0
        shared_point=0
        if(workspace_view_hit(s,x,y))then
            if(s%mesh_pick.or.s%mesh_only)then
                call pick_mesh_cell(s,x,y,mesh_target(1),mesh_target(2),mesh_target(3),face)
            else if(s%transform_mode==0.and.picking_controls(s).and.multipatch_control_scope(s))then
                shared_point=pick_shared_control(s,x,y)
            else if(s%surface_tool==1.and.s%transform_mode==0.and..not.picking_controls(s))then
                call pick_scene_patch(s,x,y,target(1),target(2),target(3),target(4))
            end if
        end if
        shared_changed=shared_point/=s%canvas_gesture%viewport_shared_hover
        if(.not.shared_changed.and.all(target==s%canvas_gesture%viewport_hover).and. &
            all(mesh_target==s%canvas_gesture%viewport_mesh_hover))return
        s%canvas_gesture%viewport_hover=target
        s%canvas_gesture%viewport_mesh_hover=mesh_target
        s%canvas_gesture%viewport_shared_hover=shared_point
        if(shared_changed.and.associated(notify_control_selection))call notify_control_selection(s)
        call s%surface%update()
    end subroutine update_viewport_hover

    subroutine surface_release(event,data)
        type(QMouseEvent), intent(in) :: event
        class(*), intent(inout), optional :: data
        integer :: point,kind,variant,element,number,mesh,part,cell,face,mode
        real(dp) :: position(3),linear(3,3),pivot(3),origin(3)
        logical :: changed
        type(selection_target_t) :: item
        select type(s=>data)
        type is(studio_t)
            if(.not.s%dragging.or.event%button()/=s%canvas_gesture%pointer_button)return
            s%canvas_gesture%pointer_button=Qt_NoButton
            s%canvas_gesture%middle_orbit=.false.
            s%dragging=.false.
            call s%status%setActivity('')
            if(s%canvas_gesture%area_count>0)then
                if(finish_area_drag(s,event))return
            end if
            call s%surface%update()
            if(transform_dragging(s))then
                call finish_transform_drag(s,event%x(),event%y(),iand(event%modifiers(),Qt_ShiftModifier)/=0, &
                    linear,position,pivot,changed)
                if(changed.and.associated(commit_transform))call commit_transform(s,linear,position,pivot)
                return
            end if
            if(event%button()==Qt_RightButton.and..not.s%canvas_gesture%pointer_displaced.and.associated(open_viewport_context))then
                if(hypot(event%x()-s%last_x,event%y()-s%last_y)<orbit_drag_threshold)then
                    if(workspace_view_hit(s,event%x(),event%y()))call open_viewport_context(s,event%x(),event%y())
                end if
            end if
            if(.not.s%control_dragging)then
                if(event%button()/=Qt_LeftButton.or.s%canvas_gesture%pointer_displaced)return
                if(hypot(event%x()-s%last_x,event%y()-s%last_y)>=orbit_drag_threshold)return
                if(.not.associated(notify_selection).or..not.workspace_view_hit(s,event%x(),event%y()))return
                mode=pointer_selection_mode(event%modifiers())
                if(s%mesh_pick.or.s%mesh_only)then
                    call pick_mesh_cell(s,event%x(),event%y(),mesh,part,cell,face)
                    item=selection_target_t()
                    if(mesh>0)then
                        item=selection_target_t(selection_mesh,mesh)
                        if(s%mesh_pick)item=selection_target_t(selection_mesh,mesh,part=part,item=cell,feature=face)
                    end if
                    if(mesh>0.or.mode==selection_replace)call notify_selection(s,item,mode)
                    return
                end if
                if(picking_controls(s))then
                    if(select_control_at_pointer(s,event%x(),event%y(),mode))return
                end if
                if(s%feature_pick_mode==pick_control_points)then
                    if(mode==selection_replace)then
                        if(multipatch_control_scope(s))then
                            call change_shared_selection(s,[integer ::],selection_replace)
                        else
                            call notify_selection(s,selection_target_t(),mode)
                        end if
                    end if
                    return
                end if
                call pick_scene_patch(s,event%x(),event%y(),kind,variant,element,number)
                if(kind==0)then
                    if(mode==selection_replace)call notify_selection(s,selection_target_t(),mode)
                    return
                end if
                item=viewport_patch_target(s,kind,variant)
                if(s%feature_pick_mode==pick_multipatch)then
                    call notify_selection(s,item,mode)
                    return
                end if
                if(s%feature_pick_mode/=pick_automatic)then
                    item%feature=element
                    item%item=number
                end if
                call notify_selection(s,item,mode)
                if(mode==selection_replace)then
                    if(s%transform_mode==0.and.s%surface_tool==1.and.s%feature_pick_mode==pick_automatic) then
                        call pick_live_viewport(s,event%x(),event%y())
                    end if
                    s%viewport_feature=[kind,variant,element,number]
                    call s%status%setSelection(selection_status_text(s),compact=selection_status_text(s,compact=.true.))
                end if
                call s%surface%update()
                return
            end if
            call update_live_control_drag(s,event%x(),event%y(),iand(event%modifiers(),Qt_ShiftModifier)/=0,with_dependencies=.false.)
            s%control_dragging=.false.
            s%control_snapped=.false.
            call clear_control_dependency_preview(s)
            call s%surface%update()
            point=s%selected_control
            if(.not.associated(commit_control_move))return
            if(.not.viewport_control_position(s,point,origin))return
            position=s%control_preview
            if(any(position/=origin))then
                call commit_control_move(s,point,position)
            else if(multipatch_control_scope(s))then
                if(hypot(event%x()-s%last_x,event%y()-s%last_y)<orbit_drag_threshold)then
                    call change_shared_selection(s,[point],selection_replace)
                end if
            end if
        end select
    end subroutine surface_release

    subroutine surface_reset(event,data)
        type(QMouseEvent), intent(in) :: event
        class(*), intent(inout), optional :: data
        select type(s=>data)
        type is(studio_t)
            if(event%button()==Qt_RightButton)return
            if(shifted_middle_button(event))then
                call surface_press(event,s)
                return
            end if
            if(s%transform_mode>0.and.event%button()==Qt_LeftButton)return
            if(s%area_selection_tool/=area_click.and.event%button()==Qt_LeftButton)then
                call cancel_area_drag(s)
                return
            end if
            if(.not.workspace_view_hit(s,event%x(),event%y()))return
            if(picking_controls(s).and.event%button()==Qt_LeftButton)then
                if(select_control_at_pointer(s,event%x(),event%y(),selection_replace))then
                    s%dragging=.false.
                    s%control_dragging=.false.
                    call clear_control_dependency_preview(s)
                    call s%surface%update()
                    if(associated(open_selected_control))call open_selected_control(s)
                    return
                end if
            end if
            call set_viewport_orientation(s,1)
            call reset_live_viewport_zoom(s)
            call frame_workspace_scene(s)
            s%dragging=.false.
            s%control_dragging=.false.
            s%control_snapped=.false.
            s%canvas_gesture%pointer_button=Qt_NoButton
            call clear_control_dependency_preview(s)
            call s%surface%update()
        end select
    end subroutine surface_reset

    subroutine surface_key(event,data)
        type(QKeyEvent), intent(in) :: event
        class(*), intent(inout), optional :: data
        if(event%key()/=Qt_Key_Escape)return
        select type(s=>data)
        type is(studio_t)
            if(.not.(s%surface%hasFocus().or.s%dragging.or.s%control_dragging.or.transform_dragging(s)))return
            if(s%canvas_gesture%area_count>0)then
                call cancel_area_drag(s)
                call event%accept()
                return
            end if
            if(s%transform_mode>0.or.transform_dragging(s))then
                if(transform_dragging(s))then
                    call cancel_transform_drag(s)
                else
                    s%transform_mode=0
                    if(associated(notify_transform_exit))call notify_transform_exit(s)
                end if
                s%dragging=.false.
                s%canvas_gesture%pointer_button=Qt_NoButton
                call s%status%setActivity('')
                call s%surface%update()
                call event%accept()
                return
            end if
            if(multipatch_control_scope(s).and.(s%control_dragging.or.s%selected_control>0))then
                if(.not.s%control_dragging)call change_shared_selection(s,[integer ::],selection_replace)
                s%dragging=.false.
                s%control_dragging=.false.
                s%control_snapped=.false.
                s%canvas_gesture%pointer_button=Qt_NoButton
                call clear_control_dependency_preview(s)
                call s%status%setActivity('')
                call s%surface%update()
                call event%accept()
                return
            end if
            s%dragging=.false.
            s%control_dragging=.false.
            s%control_snapped=.false.
            call clear_control_dependency_preview(s)
            s%viewport_feature=0
            call s%status%setActivity('')
            s%canvas_gesture%pointer_button=Qt_NoButton
            if(associated(notify_selection))call notify_selection(s,selection_target_t(),selection_replace)
            call s%surface%update()
            call event%accept()
        end select
    end subroutine surface_key

    subroutine begin_area_drag(s,x,y)
        type(studio_t), intent(inout) :: s
        real(dp), intent(in) :: x,y

        associate(gesture=>s%canvas_gesture)
            if(.not.allocated(gesture%area_points))allocate(gesture%area_points(2,16))
            gesture%area_count=1
            gesture%area_points(:,1)=[x,y]
        end associate
        call s%surface%setCursor(Qt_CrossCursor)
    end subroutine begin_area_drag

    subroutine track_area_drag(s,x,y,released)
        type(studio_t), intent(inout) :: s
        real(dp), intent(in) :: x,y
        logical, intent(in) :: released
        real(dp), allocatable :: grown(:,:)
        real(dp) :: bounds(4),point(2)
        integer :: n

        bounds=workspace_view_bounds(s)
        point=min(max([x,y],bounds(:2)),bounds(3:4))
        associate(gesture=>s%canvas_gesture)
            if(.not.gesture%pointer_displaced)then
                if(norm2(point-gesture%area_points(:,1))<orbit_drag_threshold)return
                gesture%pointer_displaced=.true.
            end if
            n=gesture%area_count
            if(s%area_selection_tool==area_box)then
                n=2
            else
                if(.not.released.and.norm2(point-gesture%area_points(:,n))<2._dp)return
                if(n==size(gesture%area_points,2))then
                    allocate(grown(2,2*n))
                    grown(:,:n)=gesture%area_points
                    call move_alloc(grown,gesture%area_points)
                end if
                n=n+1
            end if
            gesture%area_points(:,n)=point
            gesture%area_count=n
        end associate
        if(.not.released)call s%workspace%update()
    end subroutine track_area_drag

    logical function finish_area_drag(s,event) result(consumed)
        type(studio_t), intent(inout), target :: s
        type(QMouseEvent), intent(in) :: event
        type(selection_target_t), allocatable :: targets(:)
        integer :: n,tool,mode

        call track_area_drag(s,event%x(),event%y(),.true.)
        consumed=s%canvas_gesture%pointer_displaced
        n=s%canvas_gesture%area_count
        tool=s%area_selection_tool
        if(consumed)then
            associate(points=>s%canvas_gesture%area_points(:,:n))
                if((tool==area_box.and.all(abs(points(:,2)-points(:,1))>=1._dp)).or.(tool==area_lasso.and.n>=3))then
                    call collect_viewport_region(s,points,tool,targets)
                end if
            end associate
        end if
        call cancel_area_drag(s)
        if(.not.allocated(targets))return
        mode=pointer_selection_mode(event%modifiers())
        if(picking_shared_controls(s))then
            call change_shared_selection(s,targets%item,mode)
        else
            call notify_region_selection(s,targets,mode)
        end if
    end function finish_area_drag

    subroutine cancel_area_drag(s)
        type(studio_t), intent(inout) :: s

        if(s%canvas_gesture%area_count==0)return
        s%canvas_gesture%area_count=0
        s%canvas_gesture%pointer_button=Qt_NoButton
        s%dragging=.false.
        call s%surface%unsetCursor()
        call s%workspace%update()
    end subroutine cancel_area_drag

    logical function picking_controls(s) result(enabled)
        type(studio_t), intent(in) :: s

        enabled=s%layers(4).and..not.(s%mesh_only.or.s%mesh_pick).and. &
            ((s%surface_tool==1.and.s%feature_pick_mode==pick_automatic).or.s%feature_pick_mode==pick_control_points)
    end function picking_controls

    logical function picking_shared_controls(s) result(enabled)
        type(studio_t),intent(in) :: s

        enabled=s%feature_pick_mode==pick_control_points.and.multipatch_control_scope(s).and..not.(s%mesh_pick.or.s%mesh_only)
    end function picking_shared_controls

    logical function select_control_at_pointer(s,x,y,mode,drag) result(hit)
        type(studio_t), intent(inout), target :: s
        real(dp), intent(in) :: x,y
        integer, intent(in) :: mode
        logical, intent(in), optional :: drag
        integer :: kind,id,point
        logical :: start_drag,already_selected

        if(multipatch_control_scope(s))then
            point=pick_shared_control(s,x,y)
            hit=point>0
            if(.not.hit)return
            start_drag=.false.
            if(present(drag))start_drag=drag
            already_selected=.false.
            if(allocated(s%selected_shared_controls))already_selected=any(s%selected_shared_controls==point)
            if(start_drag.and.already_selected)then
                s%selected_control=point
                if(associated(notify_control_selection))call notify_control_selection(s)
            else
                call change_shared_selection(s,[point],mode)
            end if
            if(start_drag)hit=begin_live_control_drag(s,x,y,point)
            return
        end if
        call pick_selected_control(s,x,y,kind,id,point)
        hit=point>0
        if(.not.hit)return
        if(associated(notify_selection))then
            call notify_selection(s,selection_target_t(selection_control,id,kind,item=point),mode)
        else if(kind==s%entity_kind.and.id==s%live_scene%variant)then
            s%selected_control=point
            if(associated(notify_control_selection))call notify_control_selection(s)
        else
            hit=.false.
            return
        end if
        if(.not.present(drag))return
        if(.not.drag)return
        hit=.false.
        if(kind/=s%entity_kind.or.id/=s%live_scene%variant.or.point/=s%selected_control)return
        hit=begin_live_control_drag(s,x,y,point)
    end function select_control_at_pointer

    subroutine change_shared_selection(s,points,mode)
        type(studio_t),intent(inout),target :: s
        integer,intent(in) :: points(:),mode
        logical,allocatable :: selected(:)
        integer :: point,i

        if(.not.s%shared_controls%valid)return
        allocate(selected(size(s%shared_controls%points,2)),source=.false.)
        if(mode/=selection_replace)then
            if(allocated(s%selected_shared_controls))then
                do i=1,size(s%selected_shared_controls)
                    point=s%selected_shared_controls(i)
                    if(point>=1.and.point<=size(selected))selected(point)=.true.
                end do
            else
                point=s%selected_control
                if(point>=1.and.point<=size(selected))selected(point)=.true.
            end if
        end if
        do i=1,size(points)
            point=points(i)
            if(point<1.or.point>size(selected))cycle
            if(mode==selection_toggle)then
                selected(point)=.not.selected(point)
            else
                selected(point)=.true.
            end if
        end do
        s%selected_shared_controls=pack([(i,i=1,size(selected))],selected)
        s%selected_control=0
        if(size(s%selected_shared_controls)>0)then
            s%selected_control=s%selected_shared_controls(size(s%selected_shared_controls))
            do i=size(points),1,-1
                point=points(i)
                if(point<1.or.point>size(selected))cycle
                if(.not.selected(point))cycle
                s%selected_control=point
                exit
            end do
        end if
        if(associated(notify_control_selection))call notify_control_selection(s)
        call s%surface%update()
    end subroutine change_shared_selection

    subroutine surface_wheel(event,data)
        type(QWheelEvent), intent(in) :: event
        class(*), intent(inout), optional :: data
        type(QPoint) :: delta
        real(dp) :: steps
        select type(s=>data)
        type is(studio_t)
            if(s%control_dragging.or.transform_dragging(s))return
            if(.not.workspace_view_hit(s,event%x(),event%y()))return
            steps=real(event%angleDeltaY(),dp)/120._dp
            if(steps==0._dp)then
                delta=event%pixelDelta()
                steps=real(delta%y,dp)/40._dp
            end if
            if(steps==0._dp)return
            call clear_viewport_hover(s)
            call zoom_live_viewport(s,steps,event%x(),event%y())
            call event%accept()
        end select
    end subroutine surface_wheel

    subroutine zoom_domain_view(s,factor,x,y)
        type(studio_t), intent(inout), target :: s
        real(dp), intent(in) :: factor
        real(dp), intent(in), optional :: x,y
        real(dp) :: next_zoom,ratio,center(2),anchor(2)
        integer :: axis
        type(parameter_track) :: track

        if(.not.s%live_scene%valid.or.knot_edit_in_progress(s))return
        if(.not.ieee_is_finite(factor).or.factor<=0._dp)return
        next_zoom=max(min_domain_zoom,min(max_domain_zoom,s%domain_zoom*min(factor,max_domain_zoom/min_domain_zoom)))
        if(next_zoom==s%domain_zoom)return
        center=.5_dp*real([s%domain%width(),s%domain%height()],dp)
        anchor=center
        if(present(x).and.present(y))then
            anchor=[x,y]
            if(s%entity_kind>1)then
                axis=ruler_at(s,x,y)
                ! Rulers keep a constant pixel gap outside the geometry.
                if(axis>0)then
                    track=domain_axis_track(s,axis)
                    anchor=anchor-domain_ruler_gap*track%outside
                end if
            end if
        end if
        ratio=next_zoom/s%domain_zoom
        s%domain_pan_x=ratio*s%domain_pan_x+(1._dp-ratio)*(anchor(1)-center(1))
        if(s%entity_kind>1) s%domain_pan_y=ratio*s%domain_pan_y+(1._dp-ratio)*(anchor(2)-center(2))
        s%domain_zoom=next_zoom
        call s%domain%update()
    end subroutine zoom_domain_view

    subroutine domain_wheel(event,data)
        type(QWheelEvent), intent(in) :: event
        class(*), intent(inout), optional :: data
        type(QPoint) :: delta
        real(dp) :: steps

        select type(s=>data)
        type is(studio_t)
            if(.not.s%live_scene%valid.or.event%buttons()/=Qt_NoButton)return
            steps=real(event%angleDeltaY(),dp)/120._dp
            if(steps==0._dp)then
                delta=event%pixelDelta()
                steps=real(delta%y,dp)/40._dp
            end if
            if(steps==0._dp)return
            call zoom_domain_view(s,domain_zoom_step**max(-64._dp,min(64._dp,steps)),event%x(),event%y())
            call event%accept()
        end select
    end subroutine domain_wheel

    subroutine domain_press(event,data)
        type(QMouseEvent), intent(in) :: event
        class(*), intent(inout), optional :: data
        real(dp) :: point(3)
        type(volume_domain_projection) :: plot
        logical :: orbit

        select type(s=>data)
        type is(studio_t)
            if(.not.s%live_scene%valid)return
            if(event%button()==Qt_MiddleButton.or.(event%button()==Qt_RightButton.and.s%entity_kind==3))then
                if(s%canvas_gesture%domain_pointer<0.or.iand(event%buttons(),Qt_LeftButton)/=0.or.knot_edit_in_progress(s))return
                orbit=s%entity_kind==3.and.(event%button()==Qt_RightButton.or.shifted_middle_button(event))
                s%canvas_gesture%domain_pointer=merge(volume_orbit,domain_pan,orbit)
                s%canvas_gesture%domain_button=event%button()
                s%canvas_gesture%domain_last=[event%x(),event%y()]
                if(s%canvas_gesture%domain_pointer==domain_pan)call s%domain%setCursor(Qt_ClosedHandCursor)
                call event%accept()
                return
            end if
            if(event%button()/=Qt_LeftButton)return
            if(iand(event%buttons(),ior(Qt_MiddleButton,Qt_RightButton))/=0)return
            if(s%canvas_gesture%domain_pointer==domain_pan)call s%domain%unsetCursor()
            s%canvas_gesture%domain_pointer=0
            if(handle_domain_knot_press(event,s))return
            s%canvas_gesture%domain_pointer=ruler_at(s,event%x(),event%y())
            if(s%canvas_gesture%domain_pointer==0.and.s%entity_kind==3)then
                plot=live_volume_projection(s)
                s%canvas_gesture%volume_probe_depth=dot_product([s%u,s%v,s%w]-.5_dp,plot%depth)
                if(volume_probe_at(plot,event%x(),event%y(),s%canvas_gesture%volume_probe_depth,point))s%canvas_gesture%domain_pointer=domain_interior
            else if(s%canvas_gesture%domain_pointer==0)then
                if(planar_probe_at(s,event%x(),event%y(),point))s%canvas_gesture%domain_pointer=domain_interior
            end if
            call set_domain_probe(s,event%x(),event%y())
        end select
    end subroutine domain_press

    subroutine domain_move(event,data)
        type(QMouseEvent), intent(in) :: event
        class(*), intent(inout), optional :: data
        real(dp) :: pointer(2),delta(2)
        select type(s=>data)
        type is(studio_t)
            if(.not.s%live_scene%valid)return
            if(s%canvas_gesture%domain_pointer<0.and.iand(event%buttons(),s%canvas_gesture%domain_button)/=0)then
                if(s%canvas_gesture%domain_pointer==volume_orbit.and.s%entity_kind/=3)return
                pointer=[event%x(),event%y()]
                delta=pointer-s%canvas_gesture%domain_last
                if(s%canvas_gesture%domain_pointer==domain_pan)then
                    s%domain_pan_x=s%domain_pan_x+delta(1)
                    s%domain_pan_y=s%domain_pan_y+delta(2)
                else
                    s%domain_yaw=s%domain_yaw+.012_dp*delta(1)
                    s%domain_pitch=max(-domain_pitch_limit-cube_pitch,min(domain_pitch_limit-cube_pitch,s%domain_pitch+.012_dp*delta(2)))
                end if
                s%canvas_gesture%domain_last=pointer
                call s%domain%update()
                call event%accept()
                return
            end if
            if(iand(event%buttons(),Qt_LeftButton)==0)return
            if(handle_domain_knot_move(event,s))return
            call set_domain_probe(s,event%x(),event%y())
        end select
    end subroutine domain_move

    subroutine domain_release(event,data)
        type(QMouseEvent), intent(in) :: event
        class(*), intent(inout), optional :: data
        select type(s=>data)
        type is(studio_t)
            if(s%canvas_gesture%domain_pointer<0.and.event%button()==s%canvas_gesture%domain_button)then
                if(s%canvas_gesture%domain_pointer==domain_pan)call s%domain%unsetCursor()
                s%canvas_gesture%domain_pointer=0
                s%canvas_gesture%domain_button=Qt_NoButton
                call event%accept()
                return
            end if
            if(event%button()/=Qt_LeftButton)return
            if(s%canvas_gesture%domain_pointer<0)return
            s%canvas_gesture%domain_pointer=0
            if(handle_domain_knot_release(event,s))return
        end select
    end subroutine domain_release

    subroutine domain_reset(event,data)
        type(QMouseEvent), intent(in) :: event
        class(*), intent(inout), optional :: data
        select type(s=>data)
        type is(studio_t)
            if(event%button()/=Qt_MiddleButton.and.(event%button()/=Qt_RightButton.or.s%entity_kind/=3))return
            if(s%entity_kind==3.and.shifted_middle_button(event))then
                call domain_press(event,s)
                return
            end if
            call reset_domain_view(s)
            call event%accept()
        end select
    end subroutine domain_reset

    subroutine reset_domain_view(s)
        type(studio_t), intent(inout), target :: s

        if(knot_edit_in_progress(s))return
        s%domain_zoom=1._dp
        s%domain_pan_x=0._dp; s%domain_pan_y=0._dp
        s%domain_yaw=0._dp; s%domain_pitch=0._dp
        if(s%canvas_gesture%domain_pointer==domain_pan)call s%domain%unsetCursor()
        s%canvas_gesture%domain_pointer=0
        s%canvas_gesture%domain_button=Qt_NoButton
        call s%domain%update()
    end subroutine reset_domain_view

    subroutine set_domain_probe(s,x,y)
        type(studio_t), intent(inout) :: s
        real(dp), intent(in) :: x,y
        real(dp) :: point(3),fraction,distance,edge_length

        if(s%canvas_gesture%domain_pointer<=0)return
        if(s%canvas_gesture%domain_pointer==domain_interior)then
            if(s%entity_kind==3)then
                if(.not.volume_probe_at(live_volume_projection(s),x,y,s%canvas_gesture%volume_probe_depth,point))return
            else
                if(.not.planar_probe_at(s,x,y,point))return
            end if
            call s%probe_u%setValue(point(1))
            if(s%entity_kind>=2)call s%probe_v%setValue(point(2))
            if(s%entity_kind==3)call s%probe_w%setValue(point(3))
        else if(s%canvas_gesture%domain_pointer<=s%entity_kind)then
            call track_coordinates(domain_axis_track(s,s%canvas_gesture%domain_pointer),x,y,fraction,distance,edge_length)
            if(edge_length<min_ruler_length)return
            fraction=max(0._dp,min(1._dp,fraction))
            select case(s%canvas_gesture%domain_pointer)
            case(1)
                call s%probe_u%setValue(fraction)
            case(2)
                call s%probe_v%setValue(fraction)
            case(3)
                call s%probe_w%setValue(fraction)
            case default
                return
            end select
        end if
        call repaint_all(s)
    end subroutine set_domain_probe
end module studio_canvas_input
