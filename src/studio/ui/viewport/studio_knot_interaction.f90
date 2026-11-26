! Native pointer and keyboard knot editing for the parameter canvas.
module studio_knot_interaction
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: int64
    use qt, only: QAction, QContextMenuEvent, QContextMenuEvent_Keyboard, QKeyEvent, QMenu, QMouseEvent, &
        QPoint, Qt_Key_Delete, Qt_Key_Escape, Qt_LeftButton, Qt_ShiftModifier, Qt_StrongFocus
    use studio_state, only: dp, studio_t, knot_select_tool, knot_insert_tool, knot_remove_tool
    use studio_single_patch, only: commit_single_patch_knot, commit_single_patch_knot_values
    use studio_status_hints, only: register_status_menu,set_status_hint
    implicit none(type, external)
    private
    public :: connect_knot_interactions, domain_knot_anchor, knot_selection, knot_drag_position
    public :: handle_domain_knot_press, handle_domain_knot_move, handle_domain_knot_release
    public :: knot_edit_in_progress

    real(dp), parameter :: pick_radius=8._dp, drag_distance=3._dp

    abstract interface
        subroutine edit_report(s,value,operation,error,changed)
            import :: studio_t,dp
            implicit none(type, external)
            type(studio_t), intent(inout), target :: s
            real(dp), intent(in) :: value
            character(len=*), intent(in) :: operation,error
            logical, intent(in) :: changed
        end subroutine edit_report

        logical function domain_map(s,axis,x,y,fraction,pixels_per_unit)
            import :: studio_t,dp
            implicit none(type, external)
            type(studio_t), intent(in) :: s
            integer, intent(in) :: axis
            real(dp), intent(in) :: x,y
            real(dp), intent(out) :: fraction,pixels_per_unit
        end function domain_map

        logical function domain_anchor(s,axis,value,x,y)
            import :: studio_t,dp
            implicit none(type, external)
            type(studio_t), intent(in) :: s
            integer, intent(in) :: axis
            real(dp), intent(in) :: value
            real(dp), intent(out) :: x,y
        end function domain_anchor
    end interface

    type :: axis_location
        real(dp) :: lower,upper,span,fraction,pixels
    end type axis_location

    type :: knot_gesture
        real(dp), allocatable :: original(:)
        real(dp) :: pointer_fraction=0._dp,knot_at_start=0._dp,preview=0._dp
        real(dp) :: yaw=0._dp,pitch=0._dp
        real(dp) :: zoom=1._dp,pan_x=0._dp,pan_y=0._dp
        integer :: kind=0,axis=0,variant=0,width=0,height=0
        integer(int64) :: revision=-1_int64
        logical :: preview_visible=.false.
    end type knot_gesture

    type(studio_t), pointer :: workspace=>null()
    procedure(edit_report), pointer :: report=>null()
    procedure(domain_map), pointer :: map_domain=>null()
    procedure(domain_anchor), pointer :: anchor_domain=>null()
    type(knot_gesture) :: drag
    logical :: domain_selection=.false., editing_pointer=.false.

contains

    subroutine connect_knot_interactions(s,callback,mapper,anchor)
        type(studio_t), intent(inout), target :: s
        procedure(edit_report) :: callback
        procedure(domain_map), optional :: mapper
        procedure(domain_anchor), optional :: anchor

        workspace=>s
        report=>callback
        nullify(map_domain,anchor_domain)
        if (present(mapper)) map_domain=>mapper
        if (present(anchor)) anchor_domain=>anchor
        drag=knot_gesture()
        domain_selection=.false.
        editing_pointer=.false.
        if (.not.associated(map_domain)) return
        call s%domain%setFocusPolicy(Qt_StrongFocus)
        call s%domain%mouseDoubleClickEvent%connect(domain_double_click,s)
        call s%domain%keyPressEvent%connect(domain_key_press,s)
        call s%domain%contextMenuEvent%connect(domain_context_menu,s)
    end subroutine connect_knot_interactions

    logical function domain_knot_anchor(s,axis,value,x,y) result(valid)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: axis
        real(dp), intent(in) :: value
        real(dp), intent(out) :: x,y

        valid=.false.
        x=0._dp
        y=0._dp
        if (associated(anchor_domain)) valid=anchor_domain(s,axis,value,x,y)
    end function domain_knot_anchor

    logical function knots_loaded(s) result(loaded)
        type(studio_t), intent(in) :: s

        loaded=.false.
        if (.not.s%live_scene%valid.or.s%knot_axis<1.or.s%knot_axis>s%entity_kind) return
        if (.not.allocated(s%live_knots).or..not.allocated(s%live_knot_groups)) return
        loaded=size(s%live_knots)>=2.and.size(s%live_knot_groups)>=2
    end function knots_loaded

    logical function interior(s,value) result(editable)
        type(studio_t), intent(in) :: s
        real(dp), intent(in) :: value

        editable=value>s%live_scene%bounds(1,s%knot_axis).and.value<s%live_scene%bounds(2,s%knot_axis)
    end function interior

    logical function knot_selection(s,value) result(selected)
        type(studio_t), intent(in) :: s
        real(dp), intent(out) :: value
        integer :: row

        selected=.false.
        value=0._dp
        if (.not.knots_loaded(s)) return
        if (s%knot_table%isNull()) return
        row=s%knot_table%currentRow()+1
        if (row<1.or.row>size(s%live_knot_groups)) return
        value=s%live_knot_groups(row)
        selected=.true.
    end function knot_selection

    logical function knot_drag_position(s,value) result(visible)
        type(studio_t), intent(in), target :: s
        real(dp), intent(out) :: value

        value=0._dp
        visible=drag_matches(s).and.drag%preview_visible
        if (visible) value=drag%preview
    end function knot_drag_position

    logical function knot_edit_in_progress(s) result(active)
        type(studio_t), intent(in), target :: s
        active=drag_matches(s)
    end function knot_edit_in_progress

    logical function pointer_location(s,x,y,location) result(valid)
        type(studio_t), intent(in) :: s
        real(dp), intent(in) :: x,y
        type(axis_location), intent(out) :: location

        valid=.false.
        if (.not.knots_loaded(s).or..not.associated(map_domain)) return
        if (.not.map_domain(s,s%knot_axis,x,y,location%fraction,location%pixels)) return
        location%lower=s%live_scene%bounds(1,s%knot_axis)
        location%upper=s%live_scene%bounds(2,s%knot_axis)
        location%span=location%upper-location%lower
        valid=ieee_is_finite(location%span).and.location%span>0._dp.and.location%pixels>0._dp
    end function pointer_location

    logical function nearest_knot(s,location,value) result(found)
        type(studio_t), intent(in) :: s
        type(axis_location), intent(in) :: location
        real(dp), intent(out) :: value
        integer :: i,nearest
        real(dp) :: distance,closest

        nearest=0
        closest=pick_radius
        value=0._dp
        do i=1,size(s%live_knot_groups)
            if (s%live_knot_groups(i)<location%lower) cycle
            if (s%live_knot_groups(i)>location%upper) cycle
            distance=location%pixels*abs((s%live_knot_groups(i)-location%lower)/location%span-location%fraction)
            if (distance>closest) cycle
            closest=distance
            nearest=i
        end do
        found=nearest>0
        if (found) value=s%live_knot_groups(nearest)
    end function nearest_knot

    subroutine begin_drag(s,x,y)
        type(studio_t), intent(inout), target :: s
        real(dp), intent(in) :: x,y
        type(axis_location) :: location
        real(dp) :: value

        call clear_drag(s)
        call s%domain%setFocus()
        if (.not.pointer_location(s,x,y,location)) return
        if (.not.nearest_knot(s,location,value)) return
        domain_selection=.true.
        call report(s,value,'Select knot','',.false.)
        if (.not.interior(s,value)) return
        drag%original=s%live_knots
        drag%pointer_fraction=location%fraction
        drag%knot_at_start=value
        drag%kind=s%entity_kind
        drag%axis=s%knot_axis
        drag%variant=s%live_scene%variant
        drag%revision=s%document%patches%patch_revision(drag%kind,drag%variant)
        drag%width=s%domain%width()
        drag%height=s%domain%height()
        drag%yaw=s%domain_yaw
        drag%pitch=s%domain_pitch
        drag%zoom=s%domain_zoom
        drag%pan_x=s%domain_pan_x
        drag%pan_y=s%domain_pan_y
    end subroutine begin_drag

    logical function drag_matches(s) result(matches)
        type(studio_t), intent(in), target :: s

        matches=.false.
        if (.not.allocated(drag%original).or..not.associated(workspace,s)) return
        if (.not.knots_loaded(s).or..not.s%domain%isVisible()) return
        if (drag%kind/=s%entity_kind.or.drag%axis/=s%knot_axis.or.drag%variant/=s%live_scene%variant) return
        if (drag%revision/=s%document%patches%patch_revision(drag%kind,drag%variant)) return
        if (drag%kind==3.and.(s%domain_yaw/=drag%yaw.or.s%domain_pitch/=drag%pitch)) return
        if (s%domain_zoom/=drag%zoom.or.s%domain_pan_x/=drag%pan_x.or.s%domain_pan_y/=drag%pan_y) return
        matches=s%domain%width()==drag%width.and.s%domain%height()==drag%height
    end function drag_matches

    subroutine update_drag(s,x,y)
        type(studio_t), intent(inout), target :: s
        real(dp), intent(in) :: x,y
        type(axis_location) :: location

        if (.not.drag_matches(s)) then
            call clear_drag(s)
            return
        end if
        drag%preview_visible=.false.
        if (pointer_location(s,x,y,location)) then
            if (abs(location%fraction-drag%pointer_fraction)*location%pixels>=drag_distance) then
                drag%preview=drag%knot_at_start+(location%fraction-drag%pointer_fraction)*location%span
                drag%preview_visible=ieee_is_finite(drag%preview)
            end if
        end if
        call s%domain%update()
    end subroutine update_drag

    subroutine finish_drag(s,x,y)
        type(studio_t), intent(inout), target :: s
        real(dp), intent(in) :: x,y
        real(dp), allocatable :: candidate(:)
        real(dp) :: value,original
        logical :: changed
        integer :: patch_id
        character(:), allocatable :: error

        call update_drag(s,x,y)
        if (.not.drag%preview_visible) then
            call clear_drag(s)
            return
        end if
        candidate=drag%original
        value=drag%preview
        original=drag%knot_at_start
        patch_id=drag%variant
        call clear_drag(s)
        where(candidate==original) candidate=value
        call commit_single_patch_knot_values(s%entity_kind,s%knot_axis,candidate,changed,error, &
            patch_id=patch_id,catalog=s%document%patches)
        call report(s,merge(value,original,changed),'Move knot',error,changed)
    end subroutine finish_drag

    subroutine clear_drag(s)
        type(studio_t), intent(inout) :: s

        if (drag%preview_visible) call s%domain%update()
        drag=knot_gesture()
    end subroutine clear_drag

    subroutine insert_at_pointer(s,x,y)
        type(studio_t), intent(inout), target :: s
        real(dp), intent(in) :: x,y
        type(axis_location) :: location
        real(dp) :: value

        call clear_drag(s)
        if (.not.pointer_location(s,x,y,location)) return
        if (nearest_knot(s,location,value)) return
        value=location%lower+location%span*location%fraction
        if (.not.interior(s,value)) return
        call change_knot_count(s,value,.false.)
    end subroutine insert_at_pointer

    subroutine change_knot_count(s,value,remove)
        type(studio_t), intent(inout), target :: s
        real(dp), intent(in) :: value
        logical, intent(in) :: remove
        integer :: accepted
        character(:), allocatable :: error,operation

        if (remove) then
            operation='Remove knot'
        else
            operation='Insert knot'
        end if
        call commit_single_patch_knot(s%entity_kind,s%knot_axis,value,1,remove,accepted,error, &
            patch_id=s%live_scene%variant,catalog=s%document%patches)
        call report(s,value,operation,error,accepted>0)
    end subroutine change_knot_count

    logical function editing_knots(s,modifiers)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: modifiers

        editing_knots=associated(map_domain).and.(s%domain_edit_knots.or.iand(modifiers,Qt_ShiftModifier)/=0)
    end function editing_knots

    logical function handle_domain_knot_press(event,s) result(handled)
        type(QMouseEvent), intent(in) :: event
        type(studio_t), intent(inout), target :: s
        type(axis_location) :: location
        real(dp) :: value

        handled=.false.
        if (event%button()/=Qt_LeftButton) return
        domain_selection=.false.
        editing_pointer=editing_knots(s,event%modifiers())
        handled=editing_pointer
        if (handled) then
            call s%domain%setFocus()
            select case(s%domain_knot_tool)
            case(knot_insert_tool)
                call insert_at_pointer(s,event%x(),event%y())
            case(knot_remove_tool)
                call clear_drag(s)
                if (.not.pointer_location(s,event%x(),event%y(),location)) return
                if (.not.nearest_knot(s,location,value)) return
                if (interior(s,value)) call change_knot_count(s,value,.true.)
            case default
                call begin_drag(s,event%x(),event%y())
            end select
        else
            call clear_drag(s)
        end if
    end function handle_domain_knot_press

    logical function handle_domain_knot_move(event,s) result(handled)
        type(QMouseEvent), intent(in) :: event
        type(studio_t), intent(inout), target :: s

        handled=editing_pointer.and.iand(event%buttons(),Qt_LeftButton)/=0
        if (handled.and.allocated(drag%original)) call update_drag(s,event%x(),event%y())
    end function handle_domain_knot_move

    logical function handle_domain_knot_release(event,s) result(handled)
        type(QMouseEvent), intent(in) :: event
        type(studio_t), intent(inout), target :: s

        handled=.false.
        if (event%button()/=Qt_LeftButton) return
        handled=editing_pointer
        editing_pointer=.false.
        if (handled.and.allocated(drag%original)) call finish_drag(s,event%x(),event%y())
    end function handle_domain_knot_release

    subroutine domain_double_click(event,data)
        type(QMouseEvent), intent(in) :: event
        class(*), intent(inout), optional :: data

        if (event%button()/=Qt_LeftButton) return
        select type(s=>data)
        type is(studio_t)
            if (.not.editing_knots(s,event%modifiers())) return
            if(s%domain_knot_tool==knot_select_tool)call insert_at_pointer(s,event%x(),event%y())
            call event%accept()
        end select
    end subroutine domain_double_click

    subroutine domain_key_press(event,data)
        type(QKeyEvent), intent(in) :: event
        class(*), intent(inout), optional :: data
        real(dp) :: value

        select type(s=>data)
        type is(studio_t)
            select case(event%key())
            case(Qt_Key_Escape)
                if (allocated(drag%original)) then
                    call clear_drag(s)
                else if (s%domain_knot_tool/=knot_select_tool) then
                    s%domain_knot_tool=knot_select_tool
                    call report(s,0._dp,'Select knot tool','',.false.)
                else
                    return
                end if
            case(Qt_Key_Delete)
                if (.not.(s%domain_edit_knots.or.domain_selection)) return
                if (event%isAutoRepeat()) return
                if (.not.knot_selection(s,value)) return
                if (.not.interior(s,value)) return
                call clear_drag(s)
                call change_knot_count(s,value,.true.)
            case default
                return
            end select
            call event%accept()
        end select
    end subroutine domain_key_press

    subroutine domain_context_menu(event,data)
        type(QContextMenuEvent), intent(in) :: event
        class(*), intent(inout), optional :: data
        type(QMenu) :: menu
        type(QAction) :: remove,chosen
        type(QPoint) :: at
        type(axis_location) :: location
        real(dp) :: value,x,y
        integer :: kind,axis,variant
        integer(int64) :: revision

        select type(s=>data)
        type is(studio_t)
            if (event%reason()==QContextMenuEvent_Keyboard) then
                if (.not.(s%domain_edit_knots.or.domain_selection)) return
                if (.not.knot_selection(s,value)) return
                if (.not.domain_knot_anchor(s,s%knot_axis,value,x,y)) return
                at=s%domain%mapToGlobal(QPoint(nint(x),nint(y)))
            else
                ! Right-drag belongs to the volume camera.
                if (s%entity_kind==3.or..not.editing_knots(s,event%modifiers())) return
                if (.not.pointer_location(s,real(event%x(),dp),real(event%y(),dp),location)) return
                if (.not.nearest_knot(s,location,value)) return
                at=event%globalPos()
            end if
            if (.not.interior(s,value)) return
            call report(s,value,'Select knot','',.false.)
            kind=s%entity_kind
            axis=s%knot_axis
            variant=s%live_scene%variant
            revision=s%document%patches%patch_revision(kind,variant)
            menu=QMenu('Knot',s%domain)
            call register_status_menu(menu)
            remove=QAction('Remove one copy',menu)
            call set_status_hint(remove,'Remove one occurrence of the selected knot. Undo with Ctrl+Z.')
            call menu%addAction(remove)
            chosen=menu%exec(at%x,at%y)
            if (.not.chosen%isNull().and.s%entity_kind==kind.and.s%knot_axis==axis.and.s%live_scene%variant==variant) then
                if (s%document%patches%patch_revision(kind,variant)==revision) call change_knot_count(s,value,.true.)
            end if
            call menu%delete()
            call event%accept()
        end select
    end subroutine domain_context_menu

end module studio_knot_interaction
