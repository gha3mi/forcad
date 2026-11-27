module studio_render
    use, intrinsic :: iso_fortran_env, only: int64
    use studio_display_refinement, only: connect_display_refinement,request_display_refinement
    use studio_live_viewport, only: paint_live_viewport,paint_viewport_hud
    use studio_viewport_camera, only: connect_viewport_zoom_quality,viewport_camera,prepare_camera
    use studio_viewport_selection, only: paint_viewport_selection
    use studio_live_plots, only: paint_live_basis
    use studio_ai_cursor, only: paint_studio_ai_cursor
    use studio_knot_interaction, only: knot_selection, knot_drag_position
    use studio_domain_projection, only: domain_rect, curve_domain_line, live_planar_domain, volume_domain_projection, &
        parameter_track, live_volume_projection, project_volume, volume_axis_track, domain_axis_track, volume_domain_point, &
        volume_domain_axis_point, domain_knot_point, domain_knot_fraction, min_ruler_length
    use studio_canvas_input, only: connect_canvas_input, zoom_domain_view, reset_domain_view, control_move_handler, &
        control_selection_handler, selection_handler, transform_handler, viewport_context_handler, region_selection_handler
    use studio_theme, only: studio_active_edge, studio_frame_edge
    use qt, only: QBrush, QColor, QFont, QFontMetrics, QImage, QLinearGradient, QPainter, QPainter_Antialiasing, &
        QPainter_SmoothPixmapTransform, QPainter_TextAntialiasing, QPen, QResizeEvent, QTimer, &
        Qt_DashLine, Qt_DotLine, Qt_NoBrush, Qt_NoPen, Qt_SolidLine
    use studio_state, only: studio_t, dp, axis_rgb
    implicit none(type, external)
    private
    public :: connect_canvases, paint_surface, paint_domain, paint_basis
    public :: domain_rect, curve_domain_line
    public :: volume_domain_point, volume_domain_axis_point, domain_knot_fraction, domain_knot_point
    public :: zoom_domain_view, reset_domain_view

    character(len=3), parameter :: parameter_axes='UVW'
    character(*), parameter :: domain_tick_color='#d2e2ee'
    real(dp), parameter :: domain_tick_size=9.5_dp, domain_axis_size=11._dp
    integer, parameter :: viewport_resize_idle_ms=120
    real(dp), parameter :: viewport_resize_budget_ms=24._dp

contains

    subroutine connect_canvases(s,control_move,control_selected,control_open,selection_changed, &
                                transform_committed,viewport_context,transform_exited,region_selected)
        type(studio_t), intent(inout), target :: s
        procedure(control_move_handler), optional :: control_move
        procedure(control_selection_handler), optional :: control_selected,control_open
        procedure(selection_handler), optional :: selection_changed
        procedure(region_selection_handler), optional :: region_selected
        procedure(transform_handler), optional :: transform_committed
        procedure(viewport_context_handler), optional :: viewport_context
        procedure(control_selection_handler), optional :: transform_exited

        call connect_display_refinement(s)
        call connect_viewport_zoom_quality(s)
        s%viewport_resize_timer=QTimer(s%surface)
        call s%viewport_resize_timer%setSingleShot(.true.)
        call s%viewport_resize_timer%timeout%connect(on_viewport_resize_settled,s)
        call s%workspace%resizeEvent%connect(on_viewport_resize,s)
        call s%surface%paintEvent%connect(paint_surface,s)
        call s%domain%paintEvent%connect(paint_domain,s)
        call s%basis%paintEvent%connect(paint_basis,s)
        call connect_canvas_input(s,control_move,control_selected,control_open,selection_changed, &
            transform_committed,viewport_context,transform_exited,region_selected)
    end subroutine connect_canvases
    subroutine backdrop(p,w,h,viewport)
        type(QPainter), intent(inout) :: p
        real(dp), intent(in) :: w,h
        logical, intent(in), optional :: viewport
        type(QLinearGradient) :: g
        logical :: is_viewport

        is_viewport=.false.
        if(present(viewport))is_viewport=viewport
        g=QLinearGradient(0._dp,0._dp,0._dp,h)
        if(is_viewport)then
            call g%setColorAt(0._dp,QColor(23,33,46))
            call g%setColorAt(.52_dp,QColor(19,29,40))
            call g%setColorAt(1._dp,QColor(12,20,30))
        else
            call g%setColorAt(0._dp,QColor(16,24,31))
            call g%setColorAt(1._dp,QColor(21,30,38))
        end if
        call p%fillRect(0._dp,0._dp,w,h,g)
        call g%delete()
        call p%setRenderHint(QPainter_Antialiasing)
        call p%setRenderHint(QPainter_TextAntialiasing)
        call p%setFont(QFont(family='Noto Sans',pointSize=10._dp))
    end subroutine backdrop

    subroutine label_centered(p,x,y,txt,sz,vertical,color,bold)
        type(QPainter), intent(inout) :: p
        real(dp), intent(in) :: x,y,sz
        character(*), intent(in) :: txt
        logical, intent(in), optional :: vertical,bold
        type(QColor), intent(in), optional :: color
        type(QFontMetrics) :: metrics
        real(dp) :: left,baseline
        logical :: emphasized
        emphasized=.false.
        if(present(bold))emphasized=bold
        call p%setFont(QFont(family='Noto Sans',pointSize=sz,bold=emphasized))
        metrics=p%fontMetrics()
        left=x-.5_dp*real(metrics%horizontalAdvance(txt),dp)
        baseline=y
        if(present(vertical))then
            if(vertical)baseline=y+.5_dp*real(metrics%ascent()-metrics%descent(),dp)
        end if
        call metrics%delete()
        call p%setPen(QColor(237,243,247))
        if(present(color))call p%setPen(color)
        call p%drawText(left,baseline,txt)
    end subroutine label_centered

    subroutine label_right_centered(p,metrics,x,y,txt)
        type(QPainter), intent(inout) :: p
        type(QFontMetrics), intent(in) :: metrics
        real(dp), intent(in) :: x,y
        character(*), intent(in) :: txt
        real(dp) :: left,baseline

        left=x-real(metrics%horizontalAdvance(txt),dp)
        baseline=y+.5_dp*real(metrics%ascent()-metrics%descent(),dp)
        call p%setPen(QColor(domain_tick_color))
        call p%drawText(left,baseline,txt)
    end subroutine label_right_centered

    subroutine paint_surface(painter,data)
        type(QPainter), intent(inout) :: painter
        class(*), intent(inout), optional :: data
        real(dp) :: w,h

        select type(s=>data)
        type is(studio_t)
            w=real(s%surface%width(),dp)
            h=real(s%surface%height(),dp)
            call paint_viewport_frame(painter,s,w,h)
            call painter%setRenderHint(QPainter_Antialiasing)
            call painter%setRenderHint(QPainter_TextAntialiasing)
            call paint_viewport_hud(painter,s)
            call paint_viewport_selection(painter,s)
            call paint_studio_ai_cursor(painter,s,4)
        end select
    end subroutine paint_surface

    subroutine paint_viewport_frame(painter,s,w,h)
        type(QPainter),intent(inout) :: painter
        type(studio_t),intent(inout) :: s
        real(dp),intent(in) :: w,h
        type(QPainter) :: frame_painter
        real(dp) :: scale
        integer :: pixels(2)
        integer(int64) :: started,finished,rate
        logical :: rebuild

        scale=s%surface%devicePixelRatio()
        pixels=nint([w,h]*scale)
        rebuild=s%viewport_frame%isNull().or.s%viewport_frame_scale/=scale
        if(.not.rebuild)rebuild=any(pixels/=[s%viewport_frame%width(),s%viewport_frame%height()])
        if(rebuild)then
            if(paint_resize_preview(painter,s,w,h))return
            call s%viewport_frame%delete()
            s%viewport_frame=QImage(pixels(1),pixels(2))
            call s%viewport_frame%fill(QColor(12,20,30))
            s%viewport_frame_scale=scale
        end if

        call request_display_refinement(s)
        if(rebuild.or.s%surface%boolProperty('studioSceneNeedsPaint'))then
            call system_clock(started,rate)
            if(.not.s%viewport_frame%isNull())frame_painter=QPainter(s%viewport_frame)
            if(frame_painter%isActive())then
                call frame_painter%scale(scale,scale)
                call backdrop(frame_painter,w,h,viewport=.true.)
                call paint_live_viewport(frame_painter,s,w,h,frame_camera=s%viewport_frame_camera)
                call frame_painter%end()
                call s%surface%setProperty('studioSceneNeedsPaint',.false.)
            else
                call frame_painter%end()
                call backdrop(painter,w,h,viewport=.true.)
                call paint_live_viewport(painter,s,w,h)
                call s%viewport_frame%delete()
            end if
            call system_clock(finished)
            if(rate>0)s%canvas_gesture%frame_time_ms=1000._dp*real(finished-started,dp)/real(rate,dp)
        end if
        if(.not.s%viewport_frame%isNull())call painter%drawImage(0._dp,0._dp,w,h,s%viewport_frame)
    end subroutine paint_viewport_frame

    subroutine on_viewport_resize(event,data)
        type(QResizeEvent),intent(in) :: event
        class(*),intent(inout),optional :: data

        select type(s=>data)
        type is(studio_t)
            if(s%interactive_quality.and.s%canvas_gesture%frame_time_ms>=viewport_resize_budget_ms)then
                call s%viewport_resize_timer%start(viewport_resize_idle_ms)
            end if
        end select
    end subroutine on_viewport_resize

    subroutine on_viewport_resize_settled(data)
        class(*),intent(inout),optional :: data

        select type(s=>data)
        type is(studio_t)
            call s%workspace%update()
        end select
    end subroutine on_viewport_resize_settled

    logical function paint_resize_preview(painter,s,w,h) result(painted)
        type(QPainter),intent(inout) :: painter
        type(studio_t),intent(in) :: s
        real(dp),intent(in) :: w,h
        type(viewport_camera) :: camera
        real(dp) :: factor,origin(2),extent(2)

        painted=.false.
        if(.not.s%interactive_quality.or..not.s%viewport_resize_timer%isActive())return
        if(s%viewport_frame%isNull().or.s%viewport_frame_camera(3)<=0._dp)return
        if(s%surface%boolProperty('studioSceneNeedsPaint'))return
        if(.not.prepare_camera(s,camera))return
        factor=camera%scale/s%viewport_frame_camera(3)
        origin=camera%screen-factor*s%viewport_frame_camera(:2)
        extent=factor*real([s%viewport_frame%width(),s%viewport_frame%height()],dp)/s%viewport_frame_scale
        call painter%save()
        call backdrop(painter,w,h,viewport=.true.)
        call painter%setRenderHint(QPainter_SmoothPixmapTransform)
        call painter%drawImage(origin(1),origin(2),extent(1),extent(2),s%viewport_frame)
        call painter%restore()
        painted=.true.
    end function paint_resize_preview

    function parameter_breaks(s,axis) result(breaks)
        type(studio_t),intent(in)::s
        integer,intent(in)::axis
        real(dp),allocatable::breaks(:),values(:)
        real(dp)::lower,upper,knot
        integer::i,n

        ! A live scene carries the native parameter bounds and a complete knot vector for each active axis.
        associate(knots=>s%live_scene%axes(axis)%knots)
            lower=s%live_scene%bounds(1,axis)
            upper=s%live_scene%bounds(2,axis)
            allocate(values(size(knots)+2))
            values(1)=0._dp
            n=1
            do i=1,size(knots)
                knot=knots(i)
                if(knot<=lower.or.knot>=upper)cycle
                knot=(knot-lower)/(upper-lower)
                if(knot<=values(n))cycle
                n=n+1
                values(n)=knot
            end do
            values(n+1)=1._dp
            breaks=values(:n+1)
        end associate
    end function parameter_breaks

    subroutine parameter_overlay(s,layer,lo,hi)
        type(studio_t),intent(in)::s
        integer,intent(in)::layer
        real(dp),intent(out)::lo(3),hi(3)
        real(dp)::lower(3),span(3)
        integer::d
        d=s%entity_kind
        lo=0._dp
        hi=0._dp
        lower=s%live_scene%bounds(1,:)
        span=s%live_scene%bounds(2,:)-lower
        if(layer==4)then
            lo(:d)=(s%live_probe%support(1,:d)-lower(:d))/span(:d)
            hi(:d)=(s%live_probe%support(2,:d)-lower(:d))/span(:d)
        else
            lo(:d)=(s%live_probe%cell(1,:d)-lower(:d))/span(:d)
            hi(:d)=(s%live_probe%cell(2,:d)-lower(:d))/span(:d)
        end if
        lo=max(0._dp,min(1._dp,lo))
        hi=max(lo,min(1._dp,hi))
    end subroutine parameter_overlay

    function parameter_tick(s,axis,fraction)result(text)
        type(studio_t),intent(in)::s
        integer,intent(in)::axis
        real(dp),intent(in)::fraction
        character(len=:),allocatable::text
        real(dp)::lower,upper
        character(24)::buffer
        integer::decimal,last
        lower=s%live_scene%bounds(1,axis)
        upper=s%live_scene%bounds(2,axis)
        write(buffer,'(g0.5)')lower+fraction*(upper-lower)
        text=trim(buffer)
        decimal=index(text,'.')
        if(decimal==0.or.index(text,'E')/=0)return
        last=len(text)
        do while(last>decimal)
            if(text(last:last)/='0')exit
            last=last-1
        end do
        if(last==decimal)last=last-1
        text=text(:last)
    end function parameter_tick

    integer function ruler_divisions(pixel_length,minimum,zoom) result(divisions)
        real(dp), intent(in) :: pixel_length,zoom
        integer, intent(in) :: minimum

        divisions=minimum
        if(zoom>1._dp)divisions=max(divisions,2**min(8,max(0,floor(log(max(1._dp,pixel_length/72._dp))/log(2._dp)))))
    end function ruler_divisions

    subroutine paint_curve_domain(p,s)
        type(QPainter), intent(inout) :: p
        type(studio_t), intent(in) :: s
        real(dp) :: l,y,length,x,a,lo(3),hi(3),ruler_y
        real(dp),allocatable::cuts(:)
        integer :: i,divisions
        character(:), allocatable :: tick_text
        type(parameter_track) :: ruler
        call live_planar_domain(s,l,y,length)
        ruler=domain_axis_track(s,1)
        ruler_y=ruler%start(2)
        cuts=parameter_breaks(s,1)
        call p%setPen(QPen(style=Qt_NoPen))
        call p%fillRect(l,y-26._dp,length,52._dp,QColor(18,35,49))
        if(s%domain_layers(3))then
            do i=1,size(cuts)-1
                a=real(mod(i-1,2),dp)
                call p%fillRect(l+length*cuts(i),y-26._dp,length*(cuts(i+1)-cuts(i)),52._dp, &
                    QColor(22+int(31*a),75-int(14*a),97+int(6*a),125))
            end do
        end if
        if(s%live_probe%valid)then
            if(s%domain_layers(4))then
                call parameter_overlay(s,4,lo,hi)
                call p%fillRect(l+length*lo(1),y-26._dp,length*(hi(1)-lo(1)),52._dp,QColor(77,171,142,40))
            end if
            if(s%domain_layers(5))then
                call parameter_overlay(s,5,lo,hi)
                call p%fillRect(l+length*lo(1),y-26._dp,length*(hi(1)-lo(1)),52._dp,QColor(207,172,79,40))
            end if
        end if
        if(s%domain_layers(1))then
            call p%setPen(QPen(QColor(76,116,135,120),.8_dp))
            do i=0,10
                x=l+length*real(i,dp)/10._dp
                call p%drawLine(x,y-26._dp,x,y+26._dp)
            end do
        end if
        call p%setPen(QPen(QColor(studio_frame_edge),1._dp))
        call p%setBrush(QBrush(style=Qt_NoBrush))
        call p%drawRect(l,y-26._dp,length,52._dp)
        call p%setPen(QPen(QColor(138,170,186),1.1_dp))
        call p%drawLine(l,y,l+length,y)
        if(s%domain_layers(2))then
            call p%setPen(QPen(QColor(139,235,232),2._dp))
            do i=2,size(cuts)-1
                x=l+length*cuts(i)
                call p%drawLine(x,y-32._dp,x,y+32._dp)
            end do
        end if
        call paint_domain_ruler(p,s,ruler,1,s%u)
        divisions=2
        if(length>=200._dp)divisions=4
        divisions=ruler_divisions(length,divisions,s%domain_zoom)
        do i=0,divisions
            a=real(i,dp)/divisions
            x=l+length*a
            call p%setPen(QPen(QColor(194,215,223),1._dp))
            call p%drawLine(x,ruler_y,x,ruler_y+4._dp)
            tick_text=parameter_tick(s,1,a)
            call label_centered(p,x,ruler_y+21._dp,tick_text,domain_tick_size,color=QColor(domain_tick_color))
        end do
        call label_centered(p,l+length*.5_dp,ruler_y+41._dp,'U',domain_axis_size,color=parameter_axis_color(1),bold=.true.)
        if(.not.s%live_probe%valid)return
        x=l+length*s%u
        call p%setPen(QPen(QColor(177,221,222,130),1._dp,Qt_DotLine))
        call p%drawLine(x,y-34._dp,x,y+34._dp)
        call paint_domain_probe(p,x,y)
    end subroutine paint_curve_domain

    subroutine volume_segment(p,plot,start,finish)
        type(QPainter), intent(inout) :: p
        type(volume_domain_projection), intent(in) :: plot
        real(dp), intent(in) :: start(3),finish(3)
        real(dp) :: first(2),last(2)

        first=project_volume(plot,start); last=project_volume(plot,finish)
        call p%drawLine(first(1),first(2),last(1),last(2))
    end subroutine volume_segment

    subroutine volume_face(p,plot,axis,at)
        type(QPainter), intent(inout) :: p
        type(volume_domain_projection), intent(in) :: plot
        integer, intent(in) :: axis
        real(dp), intent(in) :: at
        real(dp) :: point(3),screen(2),x(4),y(4)
        integer :: i,j,k

        j=mod(axis,3)+1; k=mod(axis+1,3)+1
        point(axis)=at
        do i=1,4
            point(j)=real(merge(1,0,i==2.or.i==3),dp)
            point(k)=real(merge(1,0,i>=3),dp)
            screen=project_volume(plot,point)
            x(i)=screen(1); y(i)=screen(2)
        end do
        call p%drawPolygon(x,y)
    end subroutine volume_face

    subroutine volume_outline(p,plot,lo,hi,front,style,color)
        type(QPainter), intent(inout) :: p
        type(volume_domain_projection), intent(in) :: plot
        real(dp), intent(in) :: lo(3),hi(3)
        logical, intent(in) :: front
        integer, intent(in) :: style
        type(QColor), intent(in) :: color
        real(dp) :: start(3),finish(3)
        logical :: near(3),j_high,k_high
        integer :: axis,j,k,corner

        near=plot%depth>=0._dp
        call p%setPen(QPen(color,merge(1.35_dp,.9_dp,front),style))
        do axis=1,3
            j=mod(axis,3)+1; k=mod(axis+1,3)+1
            do corner=0,3
                j_high=mod(corner,2)==1; k_high=corner>=2
                if(front.neqv.((near(j).eqv.j_high).or.(near(k).eqv.k_high)))cycle
                start=lo
                start(j)=merge(hi(j),lo(j),j_high)
                start(k)=merge(hi(k),lo(k),k_high)
                finish=start; finish(axis)=hi(axis)
                call volume_segment(p,plot,start,finish)
            end do
        end do
    end subroutine volume_outline

    subroutine volume_knot_plane(p,plot,axis,at,pen_color,line_width,line_style)
        type(QPainter), intent(inout) :: p
        type(volume_domain_projection), intent(in) :: plot
        integer, intent(in) :: axis,line_style
        real(dp), intent(in) :: at,line_width
        type(QColor), intent(in) :: pen_color

        call p%setBrush(QBrush(style=Qt_NoBrush))
        call p%setPen(QPen(pen_color,line_width,line_style))
        call volume_face(p,plot,axis,at)
    end subroutine volume_knot_plane

    subroutine paint_volume_lattice(p,plot,axis,first,second)
        type(QPainter), intent(inout) :: p
        type(volume_domain_projection), intent(in) :: plot
        integer, intent(in) :: axis
        real(dp), intent(in) :: first(:),second(:)
        real(dp) :: start(3),finish(3)
        integer :: j,k,i,m

        j=mod(axis,3)+1; k=mod(axis+1,3)+1
        call p%setPen(QPen(QColor(axis_rgb(1,axis),axis_rgb(2,axis),axis_rgb(3,axis),70),.8_dp))
        start(axis)=0._dp
        ! Dense grids use representative interior lines, while the knot planes
        ! and active-axis handles still show every native break.
        do m=2,size(second)-1,max(1,size(second)/14)
            start(k)=second(m)
            do i=2,size(first)-1,max(1,size(first)/14)
                start(j)=first(i)
                finish=start; finish(axis)=1._dp
                call volume_segment(p,plot,start,finish)
            end do
        end do
    end subroutine paint_volume_lattice

    subroutine paint_volume_track(p,s,plot,axis,probe_value)
        type(QPainter), intent(inout) :: p
        type(studio_t), intent(in) :: s
        type(volume_domain_projection), intent(in) :: plot
        integer, intent(in) :: axis
        real(dp), intent(in) :: probe_value
        type(parameter_track) :: track
        real(dp) :: position(2),text_point(2),edge_length,fraction
        integer :: i,divisions
        type(QColor) :: color

        track=volume_axis_track(plot,axis)
        edge_length=norm2(track%finish-track%start)
        if(edge_length<min_ruler_length)return
        divisions=ruler_divisions(edge_length,1,s%domain_zoom)
        color=parameter_axis_color(axis)
        call paint_domain_ruler(p,s,track,axis,probe_value)
        do i=0,divisions
            fraction=real(i,dp)/divisions
            position=track%start+fraction*(track%finish-track%start)
            call p%setPen(QPen(color,1._dp))
            call p%drawLine(position(1),position(2),position(1)+4._dp*track%outside(1), &
                         position(2)+4._dp*track%outside(2))
            text_point=position+18._dp*track%outside
            call label_centered(p,text_point(1),text_point(2),parameter_tick(s,axis,fraction), &
                                domain_tick_size,vertical=.true.,color=QColor(domain_tick_color))
        end do
        text_point=.5_dp*(track%start+track%finish)+merge(18._dp,35._dp,divisions==1)*track%outside
        call label_centered(p,text_point(1),text_point(2),parameter_axes(axis:axis),domain_axis_size, &
                            vertical=.true.,color=color,bold=.true.)
    end subroutine paint_volume_track

    function parameter_axis_color(axis) result(color)
        integer, intent(in) :: axis
        type(QColor) :: color

        color=QColor(axis_rgb(1,axis),axis_rgb(2,axis),axis_rgb(3,axis))
    end function parameter_axis_color

    subroutine paint_domain_ruler(p,s,track,axis,probe_value)
        type(QPainter), intent(inout) :: p
        type(studio_t), intent(in) :: s
        type(parameter_track), intent(in) :: track
        integer, intent(in) :: axis
        real(dp), intent(in) :: probe_value
        real(dp) :: position(2)
        type(QColor) :: color

        color=parameter_axis_color(axis)
        call p%setPen(QPen(QColor(90,115,137,125),5._dp))
        call p%drawLine(track%start(1),track%start(2),track%finish(1),track%finish(2))
        call p%setPen(QPen(color,merge(1.8_dp,1.1_dp,s%knot_axis==axis)))
        call p%drawLine(track%start(1),track%start(2),track%finish(1),track%finish(2))
        if(.not.s%live_probe%valid.or.s%domain_edit_knots)return
        position=track%start+probe_value*(track%finish-track%start)
        call p%setPen(QPen(QColor(17,31,44),1.2_dp)); call p%setBrush(color)
        call p%drawRect(position(1)-3.5_dp,position(2)-3.5_dp,7._dp,7._dp)
    end subroutine paint_domain_ruler

    subroutine paint_domain_probe(p,x,y)
        type(QPainter), intent(inout) :: p
        real(dp), intent(in) :: x,y

        call p%setPen(QPen(QColor(studio_active_edge),1.2_dp))
        call p%setBrush(QColor(24,38,48))
        call p%drawRect(x-7.5_dp,y-7.5_dp,15._dp,15._dp)
        call p%setPen(QPen(style=Qt_NoPen)); call p%setBrush(QColor(247,251,254))
        call p%drawRect(x-3.5_dp,y-3.5_dp,7._dp,7._dp)
    end subroutine paint_domain_probe

    subroutine paint_volume_domain(p,s)
        type(QPainter), intent(inout) :: p
        type(studio_t), intent(in) :: s
        type(volume_domain_projection) :: plot
        real(dp), parameter :: unit_lo(3)=0._dp,unit_hi(3)=1._dp
        real(dp) :: point(3),foot(3),position(2),lo(3),hi(3),at
        real(dp), allocatable :: uc(:),vc(:),wc(:),cuts(:)
        integer :: axis,i,opacity

        plot=live_volume_projection(s)
        point=[s%u,s%v,s%w]
        uc=parameter_breaks(s,1); vc=parameter_breaks(s,2); wc=parameter_breaks(s,3)

        call p%setPen(QPen(style=Qt_NoPen))
        do axis=1,3
            call p%setBrush(QColor(axis_rgb(1,axis),axis_rgb(2,axis),axis_rgb(3,axis),22))
            call volume_face(p,plot,axis,merge(1._dp,0._dp,plot%depth(axis)>=0._dp))
        end do
        call volume_outline(p,plot,unit_lo,unit_hi,.false.,Qt_DashLine,QColor(137,166,188,112))
        do axis=1,3
            if(s%domain_layers(1).and.plot%scale>=80._dp)then
                do i=1,3
                    at=.25_dp*real(i,dp)
                    call volume_knot_plane(p,plot,axis,at,QColor(107,148,169,44),.65_dp,Qt_DotLine)
                end do
            end if
            if(.not.s%domain_layers(2))cycle
            select case(axis)
            case(1)
                cuts=uc
            case(2)
                cuts=vc
            case(3)
                cuts=wc
            end select
            opacity=merge(155,102,s%knot_axis==axis)
            do i=2,size(cuts)-1
                call volume_knot_plane(p,plot,axis,cuts(i), &
                    QColor(axis_rgb(1,axis),axis_rgb(2,axis),axis_rgb(3,axis),opacity),1.05_dp,Qt_SolidLine)
            end do
        end do
        if(s%domain_layers(3))then
            call paint_volume_lattice(p,plot,1,vc,wc)
            call paint_volume_lattice(p,plot,2,wc,uc)
            call paint_volume_lattice(p,plot,3,uc,vc)
        end if
        call volume_outline(p,plot,unit_lo,unit_hi,.true.,Qt_SolidLine,QColor(151,190,210,200))

        if(s%live_probe%valid)then
            if(s%domain_layers(4))then
                call parameter_overlay(s,4,lo,hi)
                call volume_outline(p,plot,lo,hi,.true.,Qt_DotLine,QColor(130,229,161,212))
                call volume_outline(p,plot,lo,hi,.false.,Qt_DotLine,QColor(130,229,161,92))
            end if
            if(s%domain_layers(5))then
                call parameter_overlay(s,5,lo,hi)
                call volume_outline(p,plot,lo,hi,.true.,Qt_SolidLine,QColor(247,189,112,223))
                call volume_outline(p,plot,lo,hi,.false.,Qt_DashLine,QColor(247,189,112,109))
            end if
            do axis=1,3
                foot=point
                foot(axis)=merge(1._dp,0._dp,plot%depth(axis)<0._dp)
                call p%setPen(QPen(QColor(axis_rgb(1,axis),axis_rgb(2,axis),axis_rgb(3,axis),180),1._dp,Qt_DashLine))
                call volume_segment(p,plot,foot,point)
            end do
            position=project_volume(plot,point)
            call paint_domain_probe(p,position(1),position(2))
        end if
        do axis=1,3
            call paint_volume_track(p,s,plot,axis,point(axis))
        end do
    end subroutine paint_volume_domain

    subroutine paint_domain(painter,data)
        type(QPainter), intent(inout) :: painter
        class(*), intent(inout), optional :: data
        real(dp) :: l,t,side,w,h,x,y,a,tick_width,lo(3),hi(3),u_y,v_x
        real(dp),allocatable::uc(:),vc(:)
        integer :: i,j,shade,divisions
        character(:), allocatable :: tick_text
        type(QFontMetrics) :: metrics
        type(parameter_track) :: ruler
        select type(s=>data)
        type is(studio_t)
            w=real(s%domain%width(),dp)
            h=real(s%domain%height(),dp)
            call backdrop(painter,w,h)
            if(.not.s%live_scene%valid.or.s%live_scene%dimension/=s%entity_kind)then
                call label_centered(painter,w*.5_dp,h*.5_dp,'Parameter domain unavailable',10._dp)
                call paint_studio_ai_cursor(painter,s,5)
                return
            end if
            if(s%entity_kind==1)then
                call paint_curve_domain(painter,s)
                call paint_domain_interaction(painter,s)
                return
            end if
            if(s%entity_kind==3)then
                call paint_volume_domain(painter,s)
                call paint_domain_interaction(painter,s)
                return
            end if
            call live_planar_domain(s,l,t,side)
            uc=parameter_breaks(s,1); vc=parameter_breaks(s,2)
            call painter%setPen(QPen(style=Qt_NoPen))
            call painter%fillRect(l,t,side,side,QColor(18,35,49))
            if(s%domain_layers(3))then
                do j=1,size(vc)-1
                    do i=1,size(uc)-1
                        shade=mod(i+j,2)
                        call painter%fillRect(l+side*uc(i),t+side*(1._dp-vc(j+1)), &
                            side*(uc(i+1)-uc(i)),side*(vc(j+1)-vc(j)), &
                            QColor(16+shade*22,55-shade*8,74+shade*2))
                    end do
                end do
            end if
            if(s%domain_layers(4).and.s%live_probe%valid) then
                call parameter_overlay(s,4,lo,hi)
                call painter%fillRect(l+side*lo(1),t+side*(1._dp-hi(2)), &
                    side*(hi(1)-lo(1)),side*(hi(2)-lo(2)),QColor(84,226,115,46))
            end if
            if(s%domain_layers(5).and.s%live_probe%valid) then
                call parameter_overlay(s,5,lo,hi)
                call painter%fillRect(l+side*lo(1),t+side*(1._dp-hi(2)), &
                    side*(hi(1)-lo(1)),side*(hi(2)-lo(2)),QColor(245,185,45,56))
            end if
            if(s%domain_layers(1)) then
                call painter%setPen(QPen(QColor(66,108,127,160),.8_dp))
                do i=1,7
                    a=side*i/8._dp
                    call painter%drawLine(l+a,t,l+a,t+side)
                    call painter%drawLine(l,t+a,l+side,t+a)
                end do
            end if
            if(s%domain_layers(2)) then
                call painter%setPen(QPen(parameter_axis_color(1),1.35_dp))
                do i=2,size(uc)-1
                    x=l+side*uc(i)
                    call painter%drawLine(x,t,x,t+side)
                end do
                call painter%setPen(QPen(parameter_axis_color(2),1.35_dp))
                do i=2,size(vc)-1
                    y=t+side*(1._dp-vc(i))
                    call painter%drawLine(l,y,l+side,y)
                end do
            end if
            call painter%setPen(QPen(QColor(147,185,195),1.3_dp))
            call painter%setBrush(QBrush(style=Qt_NoBrush))
            call painter%drawRect(l,t,side,side)
            ruler=domain_axis_track(s,1)
            u_y=ruler%start(2)
            call paint_domain_ruler(painter,s,ruler,1,s%u)
            ruler=domain_axis_track(s,2)
            v_x=ruler%start(1)
            call paint_domain_ruler(painter,s,ruler,2,s%v)
            divisions=4
            if(side<200._dp)divisions=2
            divisions=ruler_divisions(side,divisions,s%domain_zoom)
            call painter%setFont(QFont(family='Noto Sans',pointSize=domain_tick_size))
            metrics=painter%fontMetrics()
            tick_width=0._dp
            do i=0,divisions
                a=real(i,dp)/divisions
                x=l+a*side
                y=t+(1._dp-a)*side
                call painter%setPen(QPen(QColor(192,209,224),1._dp))
                call painter%drawLine(x,u_y,x,u_y+4._dp)
                call painter%drawLine(v_x-4._dp,y,v_x,y)
                tick_text=parameter_tick(s,1,a)
                call label_centered(painter,x,u_y+18._dp,tick_text,domain_tick_size,color=QColor(domain_tick_color))
                tick_text=parameter_tick(s,2,a)
                tick_width=max(tick_width,real(metrics%horizontalAdvance(tick_text),dp))
                call label_right_centered(painter,metrics,v_x-6._dp,y,tick_text)
            end do
            call metrics%delete()
            call label_centered(painter,l+side*.5_dp,u_y+37._dp,'U',domain_axis_size, &
                                color=parameter_axis_color(1),bold=.true.)
            x=v_x-15._dp-tick_width
            y=t+side*.5_dp
            if(x<10._dp)then
                x=v_x
                y=t-13._dp
            end if
            call label_centered(painter,x,y,'V',domain_axis_size,vertical=.true.,color=parameter_axis_color(2),bold=.true.)
            if(.not.s%live_probe%valid)then
                call paint_domain_interaction(painter,s)
                return
            end if
            x=l+side*s%u
            y=t+side*(1._dp-s%v)
            call painter%setPen(QPen(QColor(177,221,222,130),1._dp,Qt_DotLine))
            call painter%drawLine(l,y,x,y)
            call painter%drawLine(x,y,x,t+side)
            call paint_domain_probe(painter,x,y)
            call paint_domain_interaction(painter,s)
        end select
    end subroutine paint_domain

    subroutine paint_domain_interaction(p,s)
        type(QPainter), intent(inout) :: p
        type(studio_t), intent(in) :: s
        real(dp) :: value,selected,preview,x,y
        integer :: i,axis
        logical :: has_selected,has_preview
        type(QColor) :: color

        if(s%domain_layers(2).or.s%domain_edit_knots)then
            axis=s%knot_axis
            color=parameter_axis_color(axis)
            has_selected=knot_selection(s,selected)
            has_preview=knot_drag_position(s,preview)
            if(has_selected)call paint_domain_knot(p,s,selected,Qt_SolidLine)
            if(allocated(s%live_knot_groups))then
                do i=1,size(s%live_knot_groups)
                    value=s%live_knot_groups(i)
                    if(.not.domain_knot_point(s,axis,value,x,y))cycle
                    call draw_knot_handle(p,x,y,color,has_selected.and.value==selected,s%domain_edit_knots)
                end do
            end if
            if(has_preview)then
                call paint_domain_knot(p,s,preview,Qt_DashLine)
                if(domain_knot_point(s,axis,preview,x,y))call draw_knot_handle(p,x,y,color,.true.,.true.)
            end if
        end if
        call paint_studio_ai_cursor(p,s,5)
    end subroutine paint_domain_interaction

    subroutine draw_knot_handle(p,x,y,color,selected,editing)
        type(QPainter), intent(inout) :: p
        real(dp), intent(in) :: x,y
        type(QColor), intent(in) :: color
        logical, intent(in) :: selected,editing
        real(dp) :: radius

        radius=merge(4.6_dp,2.7_dp,editing)
        if(selected)then
            call p%setPen(QPen(QColor(245,252,255,225),1.15_dp))
            call p%setBrush(QBrush(style=Qt_NoBrush))
            call p%drawEllipse(x-8._dp,y-8._dp,16._dp,16._dp)
        end if
        call p%setPen(QPen(QColor(15,29,40),1._dp))
        call p%setBrush(color)
        call p%drawPolygon([x,x+radius,x,x-radius],[y-radius,y,y+radius,y])
    end subroutine draw_knot_handle

    subroutine paint_domain_knot(p,s,value,line_style)
        type(QPainter), intent(inout) :: p
        type(studio_t), intent(in) :: s
        real(dp), intent(in) :: value
        integer, intent(in) :: line_style
        real(dp) :: lo,hi,fraction,x,y,l,t,side
        type(QColor) :: color
        integer :: axis

        axis=s%knot_axis
        lo=s%live_scene%bounds(1,axis)
        hi=s%live_scene%bounds(2,axis)
        if(value<=lo.or.value>=hi)return
        fraction=(value-lo)/(hi-lo)
        color=QColor(axis_rgb(1,axis),axis_rgb(2,axis),axis_rgb(3,axis),235)
        call p%setPen(QPen(color,2._dp,line_style))
        select case(s%entity_kind)
        case(1)
            call live_planar_domain(s,l,y,side)
            x=l+fraction*side
            call p%drawLine(x,y-35._dp,x,y+35._dp)
        case(2)
            call live_planar_domain(s,l,t,side)
            if(axis==1)then
                x=l+fraction*side
                call p%drawLine(x,t,x,t+side)
            else
                y=t+(1._dp-fraction)*side
                call p%drawLine(l,y,l+side,y)
            end if
        case(3)
            call volume_knot_plane(p,live_volume_projection(s),axis,fraction,color,1.85_dp,line_style)
        end select
    end subroutine paint_domain_knot

    subroutine paint_basis(painter,data)
        type(QPainter), intent(inout) :: painter
        class(*), intent(inout), optional :: data
        select type(s=>data)
        type is(studio_t)
            call paint_live_basis(painter,s)
            call paint_studio_ai_cursor(painter,s,6)
        end select
    end subroutine paint_basis

end module studio_render
