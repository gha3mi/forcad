module studio_domain_projection
    use studio_state, only: studio_t, dp
    implicit none(type, external)
    private

    public :: domain_rect, curve_domain_line, live_planar_domain
    public :: volume_domain_projection, parameter_track, live_volume_projection, project_volume, volume_axis_track, domain_axis_track
    public :: volume_domain_point, volume_domain_axis_point, domain_knot_point, domain_knot_fraction
    public :: track_coordinates, ruler_at, volume_probe_at, planar_probe_at
    public :: cube_yaw, cube_pitch, min_ruler_length, domain_ruler_gap

    real(dp), parameter :: cube_yaw=.53_dp, cube_pitch=.47_dp
    real(dp), parameter :: ruler_pick_radius=10._dp, min_ruler_length=14._dp, domain_ruler_gap=12._dp

    type :: volume_domain_projection
        real(dp) :: center(2),scale,screen(2,3),depth(3)
    end type volume_domain_projection

    type :: parameter_track
        real(dp) :: start(2)=0._dp,finish(2)=0._dp,outside(2)=0._dp
    end type parameter_track

contains

    subroutine domain_rect(w,h,l,t,side)
        real(dp), intent(in) :: w,h
        real(dp), intent(out) :: l,t,side
        l=54._dp
        side=max(40._dp,min(w-70._dp,h-82._dp))
        t=max(28._dp,.5_dp*(h-54._dp-side))
        if(h<310._dp) then
            t=24._dp
            side=max(40._dp,min(w-70._dp,h-78._dp))
            l=max(46._dp,.5_dp*(w-side))
        end if
    end subroutine domain_rect

    subroutine curve_domain_line(w,h,l,y,length)
        real(dp), intent(in) :: w,h
        real(dp), intent(out) :: l,y,length
        l=36._dp
        length=max(40._dp,w-72._dp)
        y=max(54._dp,min(.5_dp*(h-45._dp),h-86._dp))
    end subroutine curve_domain_line

    subroutine live_planar_domain(s,left,top,span)
        type(studio_t), intent(in) :: s
        real(dp), intent(out) :: left,top,span
        real(dp) :: width,height

        width=real(s%domain%width(),dp)
        height=real(s%domain%height(),dp)
        if(s%entity_kind==1)then
            call curve_domain_line(width,height,left,top,span)
        else
            call domain_rect(width,height,left,top,span)
            top=.5_dp*height+s%domain_zoom*(top-.5_dp*height)
        end if
        left=.5_dp*width+s%domain_zoom*(left-.5_dp*width)+s%domain_pan_x
        top=top+s%domain_pan_y
        span=span*s%domain_zoom
    end subroutine live_planar_domain

    function volume_projection(width,height,yaw,pitch) result(plot)
        real(dp), intent(in) :: width,height,yaw,pitch
        type(volume_domain_projection) :: plot
        real(dp) :: sy,cy,sp,cp

        sy=sin(cube_yaw+yaw); cy=cos(cube_yaw+yaw)
        sp=sin(cube_pitch+pitch); cp=cos(cube_pitch+pitch)
        plot%screen(1,:)=[cy,sy,0._dp]
        plot%screen(2,:)=[sy*sp,-cy*sp,-cp]
        plot%depth=[sy*cp,-cy*cp,sp]
        plot%scale=min(max(32._dp,width-98._dp)/sum(abs(plot%screen(1,:))), &
                       max(32._dp,height-88._dp)/sum(abs(plot%screen(2,:))))
        plot%center=[width*.5_dp,height*.5_dp]
    end function volume_projection

    function live_volume_projection(s) result(plot)
        type(studio_t), intent(in) :: s
        type(volume_domain_projection) :: plot

        plot=volume_projection(real(s%domain%width(),dp),real(s%domain%height(),dp),s%domain_yaw,s%domain_pitch)
        plot%scale=plot%scale*s%domain_zoom
        plot%center=plot%center+[s%domain_pan_x,s%domain_pan_y]
    end function live_volume_projection

    function project_volume(plot,point) result(screen)
        type(volume_domain_projection), intent(in) :: plot
        real(dp), intent(in) :: point(3)
        real(dp) :: screen(2)

        screen=plot%center+plot%scale*matmul(plot%screen,point-.5_dp)
    end function project_volume

    subroutine volume_domain_point(width,height,u,v,w,x,y)
        real(dp), intent(in) :: width,height,u,v,w
        real(dp), intent(out) :: x,y
        real(dp) :: screen(2)

        screen=project_volume(volume_projection(width,height,0._dp,0._dp),[u,v,w])
        x=screen(1); y=screen(2)
    end subroutine volume_domain_point

    function volume_axis_track(plot,axis) result(track)
        type(volume_domain_projection), intent(in) :: plot
        integer, intent(in) :: axis
        type(parameter_track) :: track
        real(dp), parameter :: preference(2,3)=reshape([0._dp,1._dp, 1._dp,0._dp, -1._dp,0._dp],[2,3])
        real(dp) :: normal(2),direction(2),start(3),finish(3),screen(2),outward,score,best,edge_length
        integer :: face,j,k,corner

        j=mod(axis,3)+1; k=mod(axis+1,3)+1
        direction=plot%screen(:,axis)
        edge_length=norm2(direction)
        normal=0._dp
        if(edge_length>1.e-8_dp)normal=[-direction(2),direction(1)]/edge_length
        if(dot_product(normal,preference(:,axis))<0._dp)normal=-normal
        best=-huge(1._dp)
        track%outside=normal
        do corner=0,3
            start=0._dp
            start(j)=real(mod(corner,2),dp); start(k)=real(corner/2,dp)
            finish=start; finish(axis)=1._dp
            screen=project_volume(plot,.5_dp*(start+finish))
            outward=dot_product(screen-plot%center,normal)
            face=merge(2,1,axis==3)
            score=outward+.01_dp*screen(face)
            if(score<=best)cycle
            best=score
            track%start=project_volume(plot,start)+domain_ruler_gap*normal
            track%finish=project_volume(plot,finish)+domain_ruler_gap*normal
        end do
    end function volume_axis_track

    function domain_axis_track(s,axis) result(track)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: axis
        type(parameter_track) :: track
        real(dp) :: left,top,span

        select case(s%entity_kind)
        case(1)
            call live_planar_domain(s,left,top,span)
            track%start=[left,top+36._dp]
            track%finish=track%start+[span,0._dp]
            track%outside=[0._dp,1._dp]
        case(2)
            call live_planar_domain(s,left,top,span)
            if(axis==1)then
                track%start=[left,top+span+domain_ruler_gap]
                track%finish=track%start+[span,0._dp]
                track%outside=[0._dp,1._dp]
            else
                track%start=[left-domain_ruler_gap,top+span]
                track%finish=track%start+[0._dp,-span]
                track%outside=[-1._dp,0._dp]
            end if
        case(3)
            track=volume_axis_track(live_volume_projection(s),axis)
        case default
            return
        end select
    end function domain_axis_track

    subroutine volume_domain_axis_point(s,axis,fraction,x,y)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: axis
        real(dp), intent(in) :: fraction
        real(dp), intent(out) :: x,y
        type(parameter_track) :: track
        real(dp) :: screen(2)

        track=domain_axis_track(s,axis)
        screen=track%start+fraction*(track%finish-track%start)
        x=screen(1); y=screen(2)
    end subroutine volume_domain_axis_point

    logical function domain_knot_point(s,axis,value,x,y) result(valid)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: axis
        real(dp), intent(in) :: value
        real(dp), intent(out) :: x,y
        real(dp) :: lo,hi,fraction,point(2)
        type(parameter_track) :: track

        valid=.false.; x=0._dp; y=0._dp
        if(.not.s%live_scene%valid.or.axis<1.or.axis>s%entity_kind)return
        lo=s%live_scene%bounds(1,axis); hi=s%live_scene%bounds(2,axis)
        if(hi<=lo.or.value<lo.or.value>hi)return
        fraction=(value-lo)/(hi-lo)
        track=domain_axis_track(s,axis)
        if(norm2(track%finish-track%start)<min_ruler_length)return
        point=track%start+fraction*(track%finish-track%start)
        x=point(1); y=point(2)
        valid=.true.
    end function domain_knot_point

    subroutine track_coordinates(track,x,y,fraction,distance,edge_length)
        type(parameter_track), intent(in) :: track
        real(dp), intent(in) :: x,y
        real(dp), intent(out) :: fraction,distance,edge_length
        real(dp) :: delta(2),unit(2),pointer(2)

        delta=track%finish-track%start
        edge_length=norm2(delta)
        distance=huge(1._dp); fraction=0._dp
        if(edge_length<1.e-8_dp)return
        unit=delta/edge_length
        pointer=[x,y]-track%start
        fraction=dot_product(pointer,unit)/edge_length
        distance=norm2(pointer-max(0._dp,min(1._dp,fraction))*delta)
    end subroutine track_coordinates

    integer function ruler_at(s,x,y) result(axis)
        type(studio_t), intent(in) :: s
        real(dp), intent(in) :: x,y
        real(dp) :: nearest,distance,edge_length,fraction
        integer :: candidate

        nearest=ruler_pick_radius; axis=0
        do candidate=1,s%entity_kind
            call track_coordinates(domain_axis_track(s,candidate),x,y,fraction,distance,edge_length)
            if(distance>nearest.or.edge_length<min_ruler_length)cycle
            nearest=distance
            axis=candidate
        end do
    end function ruler_at

    logical function volume_probe_at(plot,x,y,depth,point) result(hit)
        type(volume_domain_projection), intent(in) :: plot
        real(dp), intent(in) :: x,y,depth
        real(dp), intent(out) :: point(3)
        real(dp) :: origin(3),first,last,near,far
        integer :: axis

        ! A camera ray yields a bounded depth interval through the full cube.
        ! Moving the probe preserves its projection depth within that interval.
        origin=.5_dp+matmul(transpose(plot%screen),([x,y]-plot%center)/plot%scale)
        near=-huge(1._dp); far=huge(1._dp)
        hit=.false.
        do axis=1,3
            if(abs(plot%depth(axis))<1.e-8_dp)then
                if(origin(axis)<0._dp.or.origin(axis)>1._dp)return
                cycle
            end if
            first=-origin(axis)/plot%depth(axis)
            last=(1._dp-origin(axis))/plot%depth(axis)
            near=max(near,min(first,last))
            far=min(far,max(first,last))
        end do
        if(near>far)return
        point=max(0._dp,min(1._dp,origin+max(near,min(far,depth))*plot%depth))
        hit=.true.
    end function volume_probe_at

    logical function domain_knot_fraction(s,axis,x,y,fraction,pixels_per_unit) result(hit)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: axis
        real(dp), intent(in) :: x,y
        real(dp), intent(out) :: fraction,pixels_per_unit
        real(dp) :: distance

        hit=.false.; fraction=0._dp; pixels_per_unit=0._dp
        if(.not.s%live_scene%valid.or.axis<1.or.axis>s%entity_kind)return
        call track_coordinates(domain_axis_track(s,axis),x,y,fraction,distance,pixels_per_unit)
        hit=distance<=ruler_pick_radius.and.pixels_per_unit>=min_ruler_length
        if(hit)fraction=max(0._dp,min(1._dp,fraction))
    end function domain_knot_fraction

    logical function planar_probe_at(s,x,y,point) result(hit)
        type(studio_t), intent(in) :: s
        real(dp), intent(in) :: x,y
        real(dp), intent(out) :: point(3)
        real(dp) :: l,t,side

        call live_planar_domain(s,l,t,side)
        point=[s%u,s%v,s%w]
        if(s%entity_kind==1)then
            hit=abs(y-t)<=26._dp
        else
            hit=y>=t.and.y<=t+side
            point(2)=max(0._dp,min(1._dp,1._dp-(y-t)/side))
        end if
        hit=hit.and.x>=l.and.x<=l+side
        point(1)=max(0._dp,min(1._dp,(x-l)/side))
    end function planar_probe_at
end module studio_domain_projection
