module studio_viewport_camera
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use qt, only: QIcon,QPoint,QTimer,QWidget,Qt_LeftButton,Qt_MiddleButton,Qt_RightButton
    use studio_icons, only: line_icon,icon_camera_front,icon_camera_back,icon_camera_left,icon_camera_right, &
        icon_camera_top,icon_camera_bottom
    use studio_status_hints, only: set_status_hint
    use studio_state, only: studio_t, dp, studio_mesh_count, mesh_layer, mesh_lattice, scene_patch_drawn, scene_mesh_drawn, &
        viewport_fast_shaded
    use studio_selection, only: selection_target_t, selection_patch, selection_assembly, selection_control, selection_mesh
    use studio_workspace_selection, only: selection_mesh_records
    use studio_scene, only: scene_yaw, scene_pitch, frame_studio_scene, studio_scene_control_position
    implicit none(type, external)
    private

    public :: viewport_camera, live_yaw_offset, live_pitch_offset
    public :: workspace_view_bounds, workspace_view_hit, prepare_camera, rotate_for_view, project_vertex, in_drawable_view
    public :: mesh_drawable, mesh_volume_limits, fit_workspace_selection, fit_mesh_selection, fit_around_panels, fit_framed_scene_in_clear_area
    public :: zoom_live_viewport, reset_live_viewport_zoom, frame_workspace_scene, live_viewport_point
    public :: orbit_live_viewport, wrap_view_angle, viewport_zoom_range
    public :: camera_names, camera_icons, set_viewport_orientation, sync_viewport_orientation
    public :: connect_viewport_zoom_quality,set_interactive_quality,viewport_zoom_active,viewport_multisampling,viewport_raster_step
    public :: viewport_box_visible,project_viewport_box,viewport_clip_intersects,viewport_preview_raster,viewport_navigation_active

    real(dp), parameter :: live_yaw_offset=scene_yaw, live_pitch_offset=scene_pitch
    real(dp), parameter :: scene_top=80._dp, scene_bottom=38._dp, scene_margin=12._dp
    real(dp), parameter :: pi=acos(-1._dp), orbit_rate=.006_dp, zoom_step=1.16_dp
    real(dp), parameter :: camera_yaw(7)=[live_yaw_offset,0._dp,pi,.5_dp*pi,-.5_dp*pi,0._dp,0._dp]
    real(dp), parameter :: camera_pitch(7)=[live_pitch_offset,0._dp,0._dp,0._dp,0._dp,.5_dp*pi,-.5_dp*pi]
    real(dp), parameter :: orientation_tolerance=64._dp*epsilon(1._dp)
    integer, parameter :: free_orbit_icon=28
    integer, parameter :: camera_icons(7)=[25,icon_camera_front,icon_camera_back,icon_camera_left, &
        icon_camera_right,icon_camera_top,icon_camera_bottom]
    character(14), parameter :: camera_names(7)=[character(14):: &
        'Isometric view','Front view','Back view','Left view','Right view','Top view','Bottom view']
    integer, parameter :: zoom_quality_idle_ms=200, interactive_raster_area=300000

    type :: viewport_camera
        real(dp) :: center(3), screen(2), radius, scale, floor
        real(dp) :: yaw_sin, yaw_cos, pitch_sin, pitch_cos
    end type viewport_camera

contains

    function workspace_view_bounds(s) result(bounds)
        type(studio_t), intent(in) :: s
        real(dp) :: bounds(4)
        type(QPoint) :: origin

        bounds=0._dp
        if(s%surface%isNull())return
        bounds(3:4)=real([s%surface%width(),s%surface%height()],dp)
        if(s%viewport_aperture%isNull())return
        if(.not.s%viewport_aperture%isVisible())return
        origin=s%viewport_aperture%mapTo(s%surface,QPoint(0,0))
        bounds(:2)=real([origin%x,origin%y],dp)
        bounds(3:4)=bounds(:2)+real([s%viewport_aperture%width(),s%viewport_aperture%height()],dp)
    end function workspace_view_bounds

    logical function workspace_view_hit(s,x,y) result(available)
        type(studio_t), intent(in) :: s
        real(dp), intent(in) :: x,y
        integer :: panel

        available=.false.
        if(.not.in_drawable_view(real(s%surface%width(),dp),real(s%surface%height(),dp),x,y))return
        if(covers(s%command_ribbon).or.covers(s%viewport_controls).or.covers(s%lower_row))return
        if(covers(s%assistant_overlay))return
        do panel=1,size(s%panels)
            if(panel==4.or.panel==6.or.panel==7)cycle
            if(covers(s%panels(panel)))return
        end do
        available=.true.
    contains
        logical function covers(widget) result(hit)
            class(QWidget),intent(in)::widget
            type(QPoint)::origin

            hit=.false.
            if(widget%isNull())return
            if(.not.widget%isVisible())return
            origin=widget%mapTo(s%surface,QPoint(0,0))
            hit=x>=origin%x.and.x<origin%x+widget%width().and.y>=origin%y.and.y<origin%y+widget%height()
        end function covers
    end function workspace_view_hit

    logical function prepare_camera(s,camera) result(ready)
        type(studio_t), intent(in) :: s
        type(viewport_camera), intent(out) :: camera
        real(dp) :: yaw,pitch,draw_height,draw_width,canvas_center,span(2),lower(3),upper(3)
        integer :: id

        ready=.false.
        if(s%mesh_only.and.s%document%scene%framed)then
            if(.not.any([(mesh_drawable(s,id),id=1,studio_mesh_count(s))]))return
        end if
        ! Panels only affect hit testing and HUD placement. Keep world-to-screen
        ! coordinates independent of their sizes and visibility.
        canvas_center=.5_dp*real(s%surface%width(),dp)
        draw_width=2._dp*(canvas_center-scene_margin)
        draw_height=real(s%surface%height(),dp)-scene_top-scene_bottom
        if (draw_width<=0._dp.or.draw_height<=0._dp) return
        camera%screen=[canvas_center+s%pan_x,scene_top+.5_dp*draw_height+s%pan_y]
        yaw=live_yaw_offset+s%yaw
        pitch=live_pitch_offset+s%pitch
        camera%yaw_cos=cos(yaw); camera%yaw_sin=sin(yaw)
        camera%pitch_cos=cos(pitch); camera%pitch_sin=sin(pitch)
        if(s%document%scene%framed)then
            camera%center=s%document%scene%center
            camera%radius=s%document%scene%radius
            camera%floor=s%document%scene%floor
            span=s%document%scene%span
        else
            if(.not.visible_mesh_bounds(s,lower,upper))return
            camera%center=lower+.5_dp*(upper-lower)
            camera%radius=max(sqrt(tiny(1._dp)),norm2(.5_dp*(upper-lower)))
            camera%floor=lower(3)-.06_dp*camera%radius
            span=mesh_view_span(camera,.5_dp*(upper-lower))
        end if
        camera%scale=.88_dp*s%zoom*min(draw_width/span(1),draw_height/span(2))
        ready=camera%scale>0._dp.and.ieee_is_finite(camera%scale).and.all(ieee_is_finite(camera%screen))
    end function prepare_camera

    logical function visible_mesh_bounds(s,lower,upper,id) result(found)
        type(studio_t), intent(in) :: s
        real(dp), intent(out) :: lower(3),upper(3)
        integer, intent(in), optional :: id
        integer :: first,last,mesh,part,coords
        real(dp) :: part_lower(3),part_upper(3)

        found=.false.
        lower=huge(1._dp); upper=-huge(1._dp)
        first=1; last=studio_mesh_count(s)
        if(present(id))then
            if(id<first.or.id>last)return
            first=id; last=id
        end if
        do mesh=first,last
            if(.not.mesh_drawable(s,mesh))cycle
            do part=1,size(s%document%meshes(mesh)%parts)
                associate(points=>s%document%meshes(mesh)%parts(part)%points)
                    if(size(points,1)==0)cycle
                    coords=size(points,2)
                    part_lower=0._dp; part_upper=0._dp
                    part_lower(:coords)=minval(points,dim=1)
                    part_upper(:coords)=maxval(points,dim=1)
                    lower=min(lower,part_lower)
                    upper=max(upper,part_upper)
                    found=.true.
                end associate
            end do
        end do
    end function visible_mesh_bounds

    pure function mesh_view_span(camera,half) result(span)
        type(viewport_camera), intent(in) :: camera
        real(dp), intent(in) :: half(3)
        real(dp) :: span(2),corner(3),rotated(3)
        integer :: i,axis

        span=0._dp
        do i=0,7
            corner=half*[(merge(1._dp,-1._dp,btest(i,axis-1)),axis=1,3)]
            rotated=rotate_for_view(camera,corner)
            span=max(span,2._dp*abs(rotated(:2)))
        end do
        span=max(span,max(sqrt(tiny(1._dp)),.12_dp*norm2(half)))
    end function mesh_view_span

    logical function fit_workspace_selection(s,targets,isolated_scene) result(fitted)
        type(studio_t), intent(inout) :: s
        type(selection_target_t), intent(in), optional :: targets(:)
        logical, intent(in), optional :: isolated_scene
        type(viewport_camera) :: camera
        type(selection_target_t), allocatable :: fitting(:)
        integer, allocatable :: mesh_ids(:)
        real(dp) :: lower(2), upper(2), view(4), available(2), span(2), factor, zoom, pan(2), position(3)
        integer :: index, member, id
        logical :: found

        fitted=.false.
        if(.not.prepare_camera(s,camera))return
        if(present(targets))then
            fitting=targets
        else if(present(isolated_scene))then
            if(isolated_scene.and.s%isolation%count()>0)then
                fitting=pack(s%isolation%items,s%isolation%items%category==selection_patch)
                do id=1,studio_mesh_count(s)
                    if(mesh_drawable(s,id))fitting=[fitting,selection_target_t(selection_mesh,id)]
                end do
            end if
        end if
        if(.not.allocated(fitting).and.s%selection%count()>0)fitting=s%selection%items
        if(.not.allocated(fitting))return
        view=clear_fit_bounds(s)
        available=view(3:)-view(:2)
        if(any(available<=0._dp))return

        lower=huge(1._dp)
        upper=-huge(1._dp)
        found=.false.
        do index=1,size(fitting)
            associate(target=>fitting(index))
                select case(target%category)
                case(selection_patch)
                    call include_patch(target%rank,target%id,.true.)
                case(selection_assembly)
                    if(.not.allocated(s%document%multipatches%groups))cycle
                    if(target%id<1.or.target%id>size(s%document%multipatches%groups))cycle
                    associate(definition=>s%document%multipatches%groups(target%id)%definition)
                        do member=1,size(definition%members)
                            call include_patch(definition%kind,definition%members(member),.true.)
                        end do
                    end associate
                case(selection_control)
                    if(.not.s%layers(4).or..not.scene_patch_drawn(s,target%rank,target%id))cycle
                    if(studio_scene_control_position(s%document%scene,target%rank,target%id,target%item,position,s%document%patches))then
                        call include_position(position)
                    end if
                case(selection_mesh)
                    call include_mesh(target)
                case default
                    cycle
                end select
            end associate
        end do
        if(s%mesh_only)then
            mesh_ids=selection_mesh_records(s,fitting)
            do index=1,size(mesh_ids)
                id=mesh_ids(index)
                if(any(fitting%category==selection_mesh.and.fitting%id==id))cycle
                call include_mesh(selection_target_t(selection_mesh,id))
            end do
        end if
        if(.not.found)return

        span=upper-lower
        factor=1._dp
        if(maxval(span)>sqrt(epsilon(1._dp))*camera%radius)then
            span=max(span,.01_dp*maxval(span))
            factor=.88_dp*minval(available/span)/camera%scale
        end if
        zoom=s%zoom*factor
        pan=[s%pan_x,s%pan_y] + .5_dp*(view(:2)+view(3:)) - camera%screen - factor*camera%scale*(lower+.5_dp*(upper-lower))
        if(.not.ieee_is_finite(zoom).or.zoom<=0._dp.or..not.all(ieee_is_finite(pan)))return
        s%zoom=zoom
        s%pan_x=pan(1)
        s%pan_y=pan(2)
        fitted=.true.
    contains
        subroutine include_position(point)
            real(dp), intent(in) :: point(3)
            real(dp) :: projected(3)

            if(.not.all(ieee_is_finite(point)))return
            projected=rotate_for_view(camera,point-camera%center)
            lower=min(lower,projected(:2))
            upper=max(upper,projected(:2))
            found=.true.
        end subroutine include_position

        subroutine include_patch(rank,id,controls)
            integer, intent(in) :: rank,id
            logical, intent(in) :: controls
            integer :: point

            if(.not.scene_patch_drawn(s,rank,id))return
            associate(patch=>s%document%scene%catalog(rank)%entries(id)%patch)
                if(.not.patch%valid)return
                do point=1,size(patch%points,2)
                    call include_position(patch%points(:,point))
                end do
                if(.not.controls.or..not.any(s%layers(3:4)))return
                if(.not.allocated(patch%controls))return
                do point=1,size(patch%controls,2)
                    call include_position(patch%controls(:,point))
                end do
            end associate
        end subroutine include_patch

        subroutine include_mesh(target)
            type(selection_target_t), intent(in) :: target
            integer :: first,last,part,point,coords,axis,layer_min(3),layer_max(3),stride,at,node
            real(dp) :: xyz(3)

            if(.not.mesh_drawable(s,target%id))return
            associate(mesh=>s%document%meshes(target%id))
                first=1
                last=size(mesh%parts)
                if(target%item>0)then
                    if(target%part<first.or.target%part>last)return
                    first=target%part
                    last=target%part
                end if
                do part=first,last
                    associate(data=>mesh%parts(part))
                        coords=size(data%points,2)
                        if(coords<1.or.coords>3)cycle
                        if(target%item>0)then
                            if(target%item>size(data%cells,1))cycle
                            if(mesh%kind==3.and..not.mesh_cell_drawable(s,data%samples,data%cells(target%item,:)))cycle
                            do point=1,size(data%cells,2)
                                xyz=0._dp
                                xyz(:coords)=real(data%points(data%cells(target%item,point),:),dp)
                                call include_position(xyz)
                            end do
                            cycle
                        end if
                        axis=max(1,min(3,s%mesh_layer_axis))
                        call mesh_volume_limits(s,data%samples,layer_min,layer_max)
                        stride=product(data%samples(:axis-1))
                        do node=1,size(data%points,1)
                            if(mesh%kind==3.and.s%mesh_view==mesh_layer)then
                                at=modulo((node-1)/stride,data%samples(axis))
                                if(at<layer_min(axis).or.at>layer_max(axis))cycle
                            end if
                            xyz=0._dp
                            xyz(:coords)=real(data%points(node,:),dp)
                            call include_position(xyz)
                        end do
                    end associate
                end do
            end associate
        end subroutine include_mesh
    end function fit_workspace_selection

    logical function mesh_cell_drawable(s,counts,nodes) result(visible)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: counts(3),nodes(:)
        integer :: axis,at(size(nodes)),lower(3),upper(3),stride
        logical :: face

        visible=.false.
        if(any(counts<2))return
        call mesh_volume_limits(s,counts,lower,upper)
        face=.false.
        do axis=1,3
            stride=product(counts(:axis-1))
            at=modulo((nodes-1)/stride,counts(axis))
            if(any(at<lower(axis)).or.any(at>upper(axis)))return
            face=face.or.any(at==lower(axis)).or.any(at==upper(axis))
        end do
        visible=face.or.s%mesh_view==mesh_lattice
    end function mesh_cell_drawable

    subroutine mesh_volume_limits(s,counts,lower,upper)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: counts(3)
        integer, intent(out) :: lower(3),upper(3)
        integer :: axis,layer

        lower=0
        upper=counts-1
        if(s%mesh_view/=mesh_layer.or.any(counts<2))return
        axis=max(1,min(3,s%mesh_layer_axis))
        layer=max(1,min(counts(axis)-1,s%mesh_layer_index))
        lower(axis)=layer-1
        upper(axis)=layer
    end subroutine mesh_volume_limits

    logical function fit_mesh_selection(s,id) result(fitted)
        type(studio_t), intent(inout) :: s
        integer, intent(in) :: id
        real(dp) :: lower(3),upper(3)

        fitted=.false.
        if(.not.visible_mesh_bounds(s,lower,upper,id))return
        call fit_mesh_bounds(s,lower,upper,.true.,fitted)
    end function fit_mesh_selection

    subroutine fit_mesh_bounds(s,lower,upper,clear_area,fitted)
        type(studio_t), intent(inout) :: s
        real(dp), intent(in) :: lower(3),upper(3)
        logical, intent(in) :: clear_area
        logical, intent(out), optional :: fitted
        type(viewport_camera) :: camera
        real(dp) :: center(3),span(2),bounds(4),available(2),factor,pan(2)

        if(present(fitted))fitted=.false.
        if(.not.prepare_camera(s,camera))return
        if(clear_area)then
            bounds=clear_fit_bounds(s)
        else
            bounds=[scene_margin,scene_top,real(s%surface%width(),dp)-scene_margin,real(s%surface%height(),dp)-scene_bottom]
        end if
        available=bounds(3:)-bounds(:2)
        if(any(available<=0._dp))return
        span=mesh_view_span(camera,.5_dp*(upper-lower))
        center=project_vertex(camera,lower+.5_dp*(upper-lower))
        factor=.88_dp*minval(available/span)/camera%scale
        pan=[s%pan_x,s%pan_y]+.5_dp*(bounds(:2)+bounds(3:))-camera%screen-factor*(center(:2)-camera%screen)
        s%zoom=s%zoom*factor
        s%pan_x=pan(1); s%pan_y=pan(2)
        if(present(fitted))fitted=.true.
    end subroutine fit_mesh_bounds

    subroutine fit_around_panels(s)
        type(studio_t), intent(inout) :: s
        real(dp) :: lower(3),upper(3)

        if(s%isolation%count()>0)then
            if(fit_workspace_selection(s,isolated_scene=.true.))call s%surface%update()
            return
        end if
        if(s%mesh_only)then
            if(visible_mesh_bounds(s,lower,upper))call fit_mesh_bounds(s,lower,upper,.true.)
            return
        end if
        call frame_workspace_scene(s)
        call fit_framed_scene_in_clear_area(s)
    end subroutine fit_around_panels

    function clear_fit_bounds(s) result(bounds)
        type(studio_t), intent(in) :: s
        real(dp) :: bounds(4)
        type(QPoint) :: origin

        bounds=workspace_view_bounds(s)
        bounds(:2)=bounds(:2)+[scene_margin,42._dp]
        if(.not.s%viewport_controls%isNull())then
            if(s%viewport_controls%isVisible())then
                origin=s%viewport_controls%mapTo(s%surface,QPoint(0,0))
                bounds(2)=max(bounds(2),real(origin%y+s%viewport_controls%height()+scene_margin,dp))
            end if
        end if
        bounds(3:)=bounds(3:)-[scene_margin,scene_bottom]
    end function clear_fit_bounds

    subroutine fit_framed_scene_in_clear_area(s)
        type(studio_t), intent(inout) :: s
        type(viewport_camera) :: camera
        real(dp) :: bounds(4),available(2),center(2),span(2),lower(3),upper(3)

        call reset_live_viewport_zoom(s)
        if(.not.prepare_camera(s,camera))return
        bounds=clear_fit_bounds(s)
        available=bounds(3:)-bounds(:2)
        if(any(available<=0._dp))return
        center=.5_dp*(bounds(:2)+bounds(3:))
        span=s%document%scene%span
        if(.not.s%document%scene%framed)then
            if(.not.visible_mesh_bounds(s,lower,upper))return
            span=mesh_view_span(camera,.5_dp*(upper-lower))
        end if
        s%zoom=.88_dp*minval(available/span)/camera%scale
        s%pan_x=center(1)-camera%screen(1)
        s%pan_y=center(2)-camera%screen(2)
    end subroutine fit_framed_scene_in_clear_area

    pure real(dp) function wrap_view_angle(angle) result(wrapped)
        real(dp),intent(in) :: angle

        wrapped=modulo(angle,2._dp*pi)
        if(wrapped>pi)wrapped=wrapped-2._dp*pi
    end function wrap_view_angle

    subroutine set_viewport_orientation(s,orientation)
        type(studio_t),intent(inout) :: s
        integer,intent(in) :: orientation

        if(orientation<1.or.orientation>size(camera_yaw))return
        s%yaw=camera_yaw(orientation)-live_yaw_offset
        s%pitch=camera_pitch(orientation)-live_pitch_offset
        call sync_viewport_orientation(s)
    end subroutine set_viewport_orientation

    subroutine sync_viewport_orientation(s)
        type(studio_t),intent(inout) :: s
        type(QIcon) :: icon
        character(:),allocatable :: hint
        integer :: orientation,i,glyph

        if(s%widgets%viewport_orientation%isNull().or.s%widgets%free_orbit_action%isNull())return
        orientation=viewport_orientation(s%yaw,s%pitch)
        if(orientation==s%widgets%displayed_viewport_orientation)return
        s%widgets%displayed_viewport_orientation=orientation
        do i=1,size(s%widgets%camera_view_actions)
            call s%widgets%camera_view_actions(i)%setChecked(i==orientation)
        end do
        call s%widgets%free_orbit_action%setChecked(orientation==0)
        hint='View orientation · Free orbit'
        glyph=free_orbit_icon
        if(orientation>0)then
            hint='View orientation · '//trim(camera_names(orientation))
            glyph=camera_icons(orientation)
        end if
        icon=line_icon(glyph)
        call s%widgets%viewport_orientation%setIcon(icon)
        call icon%delete()
        call s%widgets%viewport_orientation%setAccessibleName(hint)
        call set_status_hint(s%widgets%viewport_orientation,hint)
    end subroutine sync_viewport_orientation

    pure integer function viewport_orientation(yaw,pitch) result(orientation)
        real(dp),intent(in) :: yaw,pitch
        integer :: i

        orientation=0
        do i=1,size(camera_yaw)
            if(abs(wrap_view_angle(live_yaw_offset+yaw-camera_yaw(i)))>orientation_tolerance)cycle
            if(abs(wrap_view_angle(live_pitch_offset+pitch-camera_pitch(i)))>orientation_tolerance)cycle
            orientation=i
            return
        end do
    end function viewport_orientation

    subroutine orbit_live_viewport(s,dx,dy)
        type(studio_t),intent(inout) :: s
        real(dp),intent(in) :: dx,dy

        s%yaw=wrap_view_angle(live_yaw_offset+s%yaw+orbit_rate*dx)-live_yaw_offset
        s%pitch=wrap_view_angle(live_pitch_offset+s%pitch+orbit_rate*dy)-live_pitch_offset
        call sync_viewport_orientation(s)
    end subroutine orbit_live_viewport

    function viewport_zoom_range(s) result(limits)
        type(studio_t),intent(in) :: s
        real(dp) :: limits(2),precision_at_canvas

        precision_at_canvas=epsilon(1._dp)
        if(.not.s%surface%isNull())then
            precision_at_canvas=spacing(real(max(1,s%surface%width(),s%surface%height()),dp))
        end if
        limits=[precision_at_canvas,1._dp/precision_at_canvas]
    end function viewport_zoom_range

    subroutine connect_viewport_zoom_quality(s)
        type(studio_t),intent(inout),target :: s

        if(.not.s%zoom_quality_timer%isNull())return
        s%zoom_quality_timer=QTimer(s%surface)
        call s%zoom_quality_timer%setSingleShot(.true.)
        call s%zoom_quality_timer%timeout%connect(on_zoom_settled,s)
    end subroutine connect_viewport_zoom_quality

    subroutine on_zoom_settled(data)
        class(*),intent(inout),optional :: data

        select type(s=>data)
        type is(studio_t)
            if(s%interactive_quality)call s%surface%update()
        end select
    end subroutine on_zoom_settled

    subroutine set_interactive_quality(s,enabled)
        type(studio_t),intent(inout) :: s
        logical,intent(in) :: enabled

        s%interactive_quality=enabled
        call s%surface%update()
    end subroutine set_interactive_quality

    logical function viewport_zoom_active(s) result(active)
        type(studio_t),intent(in) :: s

        active=.false.
        if(.not.s%zoom_quality_timer%isNull())active=s%zoom_quality_timer%isActive()
    end function viewport_zoom_active

    logical function viewport_multisampling(s) result(enabled)
        type(studio_t),intent(in) :: s

        enabled=.not.(s%viewport_style==viewport_fast_shaded.or.s%control_dragging.or. &
            (s%dragging.and.s%canvas_gesture%area_count==0))
        if(enabled.and.s%interactive_quality)enabled=.not.viewport_zoom_active(s)
    end function viewport_multisampling

    integer function viewport_raster_step(s,extent) result(step)
        type(studio_t),intent(in) :: s
        integer,intent(in) :: extent(2)

        step=1
        if(product(real(extent,dp))<=interactive_raster_area)return
        if(viewport_preview_raster(s))step=2
    end function viewport_raster_step

    logical function viewport_preview_raster(s) result(preview)
        type(studio_t),intent(in) :: s

        preview=s%viewport_style==viewport_fast_shaded.or.viewport_navigation_active(s)
    end function viewport_preview_raster

    logical function viewport_navigation_active(s) result(active)
        type(studio_t),intent(in) :: s
        integer :: pointer_button

        active=s%interactive_quality.and.viewport_zoom_active(s)
        if(s%dragging.and.s%canvas_gesture%pointer_displaced)then
            pointer_button=s%canvas_gesture%pointer_button
            active=active.or.any(pointer_button==[Qt_RightButton,Qt_MiddleButton])
            if(pointer_button==Qt_LeftButton)then
                active=active.or.(s%surface_tool==2.and.s%transform_mode==0.and..not.s%mesh_pick.and. &
                    s%canvas_gesture%area_count==0.and..not.s%canvas_gesture%probe_dragging.and..not.s%control_dragging)
            end if
        end if
    end function viewport_navigation_active

    subroutine zoom_live_viewport(s,steps,x,y)
        type(studio_t),intent(inout) :: s
        real(dp),intent(in) :: steps
        real(dp),intent(in),optional :: x,y
        type(viewport_camera) :: camera
        real(dp) :: limits(2),zoom,log_zoom,anchor(2),center(2),pan(2)

        if(s%control_dragging.or..not.ieee_is_finite(steps).or.steps==0._dp)return
        if(.not.ieee_is_finite(s%zoom).or.s%zoom<=0._dp)return
        limits=viewport_zoom_range(s)
        log_zoom=log(s%zoom)+log(zoom_step)*steps
        if(log_zoom<=log(limits(1)))then
            zoom=limits(1)
        else if(log_zoom>=log(limits(2)))then
            zoom=limits(2)
        else
            zoom=exp(log_zoom)
        end if
        if(zoom==s%zoom)return
        if(prepare_camera(s,camera))then
            center=camera%screen-[s%pan_x,s%pan_y]
            anchor=center
            if(present(x))anchor(1)=x
            if(present(y))anchor(2)=y
            if(.not.all(ieee_is_finite(anchor)))return
            pan=anchor-center+(zoom/s%zoom)*(camera%screen-anchor)
            if(.not.all(ieee_is_finite(pan)))return
            s%pan_x=pan(1)
            s%pan_y=pan(2)
        end if
        s%zoom=zoom
        if(.not.s%zoom_quality_timer%isNull())call s%zoom_quality_timer%start(zoom_quality_idle_ms)
        call s%surface%update()
    end subroutine zoom_live_viewport

    subroutine reset_live_viewport_zoom(s)
        type(studio_t),intent(inout) :: s

        if(.not.s%zoom_quality_timer%isNull())call s%zoom_quality_timer%stop()
        s%zoom=1._dp
        s%pan_x=0._dp; s%pan_y=0._dp
        call s%surface%update()
    end subroutine reset_live_viewport_zoom

    subroutine frame_workspace_scene(s)
        type(studio_t),intent(inout) :: s
        integer :: target(2)
        real(dp) :: lower(3),upper(3)

        call sync_viewport_orientation(s)
        if(s%isolation%count()>0)then
            if(fit_workspace_selection(s,isolated_scene=.true.))call s%surface%update()
            return
        end if
        if(s%mesh_only)then
            if(visible_mesh_bounds(s,lower,upper))call fit_mesh_bounds(s,lower,upper,.false.)
            return
        end if
        target=[s%live_scene%dimension,s%live_scene%variant]
        call frame_studio_scene(s%document%scene,s%yaw,s%pitch,include_controls=any(s%layers(3:4)).and.s%selected_multipatch==0, &
                                control_target=target)
    end subroutine frame_workspace_scene

    pure function rotate_for_view(camera,xyz) result(rotated)
        type(viewport_camera), intent(in) :: camera
        real(dp), intent(in) :: xyz(3)
        real(dp) :: rotated(3),ahead
        ahead=camera%yaw_sin*xyz(1)+camera%yaw_cos*xyz(2)
        rotated(1)=camera%yaw_cos*xyz(1)-camera%yaw_sin*xyz(2)
        rotated(2)=-camera%pitch_sin*ahead-camera%pitch_cos*xyz(3)
        rotated(3)=-camera%pitch_cos*ahead+camera%pitch_sin*xyz(3)
    end function rotate_for_view

    pure function project_vertex(camera,xyz) result(projected)
        type(viewport_camera), intent(in) :: camera
        real(dp), intent(in) :: xyz(3)
        real(dp) :: projected(3)
        projected=rotate_for_view(camera,xyz-camera%center)
        projected(:2)=camera%screen+camera%scale*projected(:2)
    end function project_vertex

    pure logical function viewport_box_visible(camera,controls,clip) result(visible)
        type(viewport_camera),intent(in) :: camera
        real(dp),intent(in) :: controls(:,:),clip(4)

        visible=.true.
        if(size(controls,1)/=3.or.size(controls,2)==0)return
        visible=viewport_clip_intersects(project_viewport_box(camera,minval(controls,dim=2),maxval(controls,dim=2)),clip)
    end function viewport_box_visible

    pure function project_viewport_box(camera,lower,upper) result(bounds)
        type(viewport_camera),intent(in) :: camera
        real(dp),intent(in) :: lower(3),upper(3)
        real(dp) :: bounds(4),half(3),center(3),reach(2)

        half=.5_dp*upper-.5_dp*lower
        center=project_vertex(camera,lower+half)
        reach=camera%scale*[abs(camera%yaw_cos)*half(1)+abs(camera%yaw_sin)*half(2), &
            abs(camera%pitch_sin)*(abs(camera%yaw_sin)*half(1)+abs(camera%yaw_cos)*half(2))+abs(camera%pitch_cos)*half(3)]
        bounds=[center(:2)-reach,center(:2)+reach]
    end function project_viewport_box

    pure logical function viewport_clip_intersects(bounds,clip) result(visible)
        real(dp),intent(in) :: bounds(4),clip(4)
        real(dp) :: margin(2)

        visible=.true.
        if(.not.all(ieee_is_finite(bounds)))return
        margin=8._dp+16._dp*epsilon(1._dp)*max(abs(bounds(:2)),abs(bounds(3:)))
        visible=all(bounds(3:)+margin>=clip(:2)).and.all(bounds(:2)-margin<=clip(3:))
    end function viewport_clip_intersects

    logical function mesh_drawable(s,id) result(visible)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: id

        visible=scene_mesh_drawn(s,id)
    end function mesh_drawable

    subroutine live_viewport_point(s,position,x,y)
        type(studio_t), intent(in) :: s
        real(dp), intent(in) :: position(3)
        real(dp), intent(out) :: x,y
        type(viewport_camera) :: camera
        real(dp) :: point(3)

        x=-1._dp; y=-1._dp
        if (.not.prepare_camera(s,camera)) return
        point=project_vertex(camera,position)
        x=point(1); y=point(2)
    end subroutine live_viewport_point

    pure logical function in_drawable_view(width,height,x,y) result(inside)
        real(dp), intent(in) :: width,height,x,y

        inside=x>=0._dp.and.x<width.and.y>=0._dp.and.y<height
    end function in_drawable_view
end module studio_viewport_camera
