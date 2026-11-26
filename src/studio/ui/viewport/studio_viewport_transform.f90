module studio_viewport_transform
    use qt, only: QColor, QFont, QPainter, QPen
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use studio_forcad, only: rotation_matrix
    use studio_state, only: studio_t, dp, viewport_transform_t, world_axes=>viewport_world_axes, scene_patch_drawn
    use studio_multipatch_collection, only: count_multipatches
    use studio_single_patch, only: single_patch_revision
    use studio_viewport_camera, only: viewport_camera, prepare_camera, project_vertex, rotate_for_view, workspace_view_hit
    use studio_viewport_primitives, only: draw_dot
    use studio_viewport_scene, only: active_patch_visible
    implicit none(type, external)
    private

    public :: transform_gizmo_available, transform_gizmo_pivot, transform_gizmo_point, transform_gizmo_hint
    public :: begin_transform_drag, update_transform_drag, finish_transform_drag, cancel_transform_drag, transform_dragging
    public :: transform_drag_matches, paint_transform_gizmo

    real(dp), parameter :: gizmo_length=67._dp,gizmo_radius=54._dp,gizmo_pick=8._dp,gizmo_slop=3._dp
    real(dp), parameter :: pi=acos(-1._dp)
    integer, parameter :: gizmo_colors(3,3)=reshape([245,128,125,116,222,171,127,177,247],[3,3])
    character(len=3), parameter :: gizmo_axes='XYZ'

contains

    subroutine selected_transform_target(s,kind,variants,ready)
        type(studio_t),intent(in) :: s
        integer,intent(out) :: kind
        integer,allocatable,intent(out) :: variants(:)
        logical,intent(out) :: ready
        integer :: id

        ready=.false.
        kind=s%live_scene%dimension
        if(s%mesh_only)return
        id=s%selected_multipatch
        if(id==0)then
            if(.not.active_patch_visible(s))return
            variants=[s%live_scene%variant]
            ready=.true.
            return
        end if
        if(id<1.or.id>count_multipatches(s%document%multipatches))return
        associate(group=>s%document%multipatches%groups(id))
            if(.not.allocated(group%native))return
            kind=group%definition%kind
            variants=group%definition%members
        end associate
        ready=any([(scene_patch_drawn(s,kind,variants(id)),id=1,size(variants))])
    end subroutine selected_transform_target

    logical function transform_gizmo_available(s) result(ready)
        type(studio_t),intent(in) :: s
        integer :: kind
        integer,allocatable :: variants(:)

        call selected_transform_target(s,kind,variants,ready)
    end function transform_gizmo_available

    function transform_gizmo_pivot(s) result(pivot)
        type(studio_t),intent(in) :: s
        real(dp) :: pivot(3)
        real(dp) :: lower(3),upper(3)
        integer :: kind,i
        integer,allocatable :: variants(:)
        logical :: ready

        pivot=0._dp
        call selected_transform_target(s,kind,variants,ready)
        if(.not.ready)return
        lower=huge(1._dp); upper=-huge(1._dp)
        do i=1,size(variants)
            associate(patch=>s%document%scene%catalog(kind)%entries(variants(i))%patch)
                if(.not.patch%valid.or..not.allocated(patch%points))cycle
                lower=min(lower,minval(patch%points,dim=2))
                upper=max(upper,maxval(patch%points,dim=2))
            end associate
        end do
        if(all(upper>=lower))pivot=.5_dp*lower+.5_dp*upper
    end function transform_gizmo_pivot

    pure function gizmo_screen_axes(camera) result(axes)
        type(viewport_camera),intent(in) :: camera
        real(dp) :: axes(2,3),projected(3)
        integer :: axis

        do axis=1,3
            projected=rotate_for_view(camera,world_axes(:,axis))
            axes(:,axis)=projected(:2)
        end do
    end function gizmo_screen_axes

    pure function gizmo_ring_plane(axes,axis) result(plane)
        real(dp),intent(in) :: axes(2,3)
        integer,intent(in) :: axis
        real(dp) :: plane(2,2)

        plane=axes(:,[mod(axis,3)+1,mod(axis+1,3)+1])
    end function gizmo_ring_plane

    pure real(dp) function determinant(plane)
        real(dp),intent(in) :: plane(2,2)

        determinant=plane(1,1)*plane(2,2)-plane(1,2)*plane(2,1)
    end function determinant

    logical function gizmo_point(axes,center,mode,axis,fraction,point) result(valid)
        real(dp),intent(in) :: axes(2,3),center(2),fraction
        integer,intent(in) :: mode,axis
        real(dp),intent(out) :: point(2)
        real(dp) :: angle,plane(2,2),length

        point=center
        valid=.false.
        if(axis==4)then
            valid=mode==1.or.mode==3
        else if(mode==2)then
            plane=gizmo_ring_plane(axes,axis)
            if(abs(determinant(plane))<.11_dp)return
            angle=2._dp*pi*fraction
            point=center+gizmo_radius*matmul(plane,[cos(angle),sin(angle)])
            valid=.true.
        else
            length=norm2(axes(:,axis))
            if(length<.18_dp)return
            point=center+gizmo_length*fraction*axes(:,axis)/length
            valid=.true.
        end if
    end function gizmo_point

    logical function transform_gizmo_point(s,axis,fraction,x,y) result(valid)
        type(studio_t),intent(in) :: s
        integer,intent(in) :: axis
        real(dp),intent(in) :: fraction
        real(dp),intent(out) :: x,y
        type(viewport_camera) :: camera
        real(dp) :: center(3),point(2)

        valid=.false.; x=0._dp; y=0._dp
        if(s%transform_mode<1.or.s%transform_mode>3.or.axis<1.or.axis>4)return
        if(.not.transform_gizmo_available(s).or..not.prepare_camera(s,camera))return
        center=project_vertex(camera,transform_gizmo_pivot(s))
        if(.not.gizmo_point(gizmo_screen_axes(camera),center(:2),s%transform_mode,axis,fraction,point))return
        x=point(1); y=point(2)
        valid=workspace_view_hit(s,x,y)
    end function transform_gizmo_point

    pure real(dp) function gizmo_segment_distance(click,a,b) result(distance)
        real(dp),intent(in) :: click(2),a(2),b(2)
        real(dp) :: segment(2),length,t

        segment=b-a
        length=sum(segment**2)
        t=0._dp
        if(length>0._dp)t=max(0._dp,min(1._dp,dot_product(click-a,segment)/length))
        distance=norm2(click-a-t*segment)
    end function gizmo_segment_distance

    integer function hit_gizmo(axes,center,mode,click) result(hit)
        real(dp),intent(in) :: axes(2,3),center(2),click(2)
        integer,intent(in) :: mode
        integer,parameter :: divisions=48
        real(dp) :: a(2),b(2),distance,nearest
        integer :: axis,i

        hit=0
        if(mode/=2.and.norm2(click-center)<=gizmo_pick)then
            hit=4
            return
        end if
        nearest=gizmo_pick
        do axis=1,3
            if(mode==2)then
                if(.not.gizmo_point(axes,center,mode,axis,0._dp,a))cycle
                distance=huge(1._dp)
                do i=1,divisions
                    if(.not.gizmo_point(axes,center,mode,axis,real(i,dp)/divisions,b))exit
                    distance=min(distance,gizmo_segment_distance(click,a,b))
                    a=b
                end do
            else
                if(.not.gizmo_point(axes,center,mode,axis,.22_dp,a))cycle
                if(.not.gizmo_point(axes,center,mode,axis,1._dp,b))cycle
                distance=gizmo_segment_distance(click,a,b)
            end if
            if(distance>=nearest)cycle
            nearest=distance
            hit=axis
        end do
    end function hit_gizmo

    pure real(dp) function gizmo_angle(plane,offset) result(angle)
        real(dp),intent(in) :: plane(2,2),offset(2)
        real(dp) :: coordinates(2)

        coordinates=[plane(2,2)*offset(1)-plane(1,2)*offset(2), &
            -plane(2,1)*offset(1)+plane(1,1)*offset(2)]/determinant(plane)
        angle=atan2(coordinates(2),coordinates(1))
    end function gizmo_angle

    pure logical function transform_dragging(s)
        type(studio_t),intent(in) :: s
        transform_dragging=s%transform_drag%axis>0
    end function transform_dragging

    logical function transform_drag_matches(s) result(matches)
        type(studio_t),intent(in) :: s
        integer :: kind,i
        integer,allocatable :: variants(:)
        logical :: ready

        matches=.false.
        if(.not.transform_dragging(s).or.s%transform_mode/=s%transform_drag%mode)return
        if(s%selected_multipatch/=s%transform_drag%group)return
        if(any([s%surface%width(),s%surface%height()]/=s%transform_drag%viewport_size))return
        if(any([s%yaw,s%pitch,s%zoom,s%pan_x,s%pan_y]/=s%transform_drag%view))return
        call selected_transform_target(s,kind,variants,ready)
        if(.not.ready.or.kind/=s%transform_drag%kind)return
        if(size(variants)/=size(s%transform_drag%variants))return
        if(any(variants/=s%transform_drag%variants))return
        do i=1,size(variants)
            if(single_patch_revision(kind,variants(i),catalog=s%document%patches)/=s%transform_drag%revisions(i))return
        end do
        matches=.true.
    end function transform_drag_matches

    logical function begin_transform_drag(s,x,y) result(hit)
        type(studio_t),intent(inout) :: s
        real(dp),intent(in) :: x,y
        type(viewport_camera) :: camera
        integer :: axis,kind,i
        integer,allocatable :: variants(:)
        real(dp) :: axes(2,3),pivot(3),center(3)
        logical :: ready

        hit=.false.
        if(s%transform_mode<1.or.s%transform_mode>3)return
        call selected_transform_target(s,kind,variants,ready)
        if(.not.ready.or..not.prepare_camera(s,camera))return
        pivot=transform_gizmo_pivot(s)
        center=project_vertex(camera,pivot)
        axes=gizmo_screen_axes(camera)
        axis=hit_gizmo(axes,center(:2),s%transform_mode,[x,y])
        if(axis==0)return
        s%transform_drag=viewport_transform_t()
        s%transform_drag%axis=axis
        s%transform_drag%mode=s%transform_mode
        s%transform_drag%kind=kind
        s%transform_drag%group=s%selected_multipatch
        s%transform_drag%viewport_size=[s%surface%width(),s%surface%height()]
        s%transform_drag%view=[s%yaw,s%pitch,s%zoom,s%pan_x,s%pan_y]
        s%transform_drag%variants=variants
        s%transform_drag%revisions=[(single_patch_revision(kind,variants(i),catalog=s%document%patches),i=1,size(variants))]
        s%transform_drag%pivot=pivot
        s%transform_drag%screen=center(:2)
        s%transform_drag%start=[x,y]
        if(s%transform_mode==2)then
            s%transform_drag%plane=gizmo_ring_plane(axes,axis)
            s%transform_drag%spin_last=gizmo_angle(s%transform_drag%plane,[x,y]-center(:2))
        else
            s%transform_drag%basis=transpose(axes)/camera%scale
            s%transform_drag%plane(:,1)=axes(:,min(axis,3))
            s%transform_drag%plane(:,2)=[camera%scale,0._dp]
            if(s%transform_mode==3)s%transform_drag%amount=1._dp
        end if
        call s%surface%update()
        hit=.true.
    end function begin_transform_drag

    subroutine update_transform_drag(s,x,y,temporary_snap)
        type(studio_t),intent(inout) :: s
        real(dp),intent(in) :: x,y
        logical,intent(in) :: temporary_snap
        real(dp) :: displacement(2),projected(2),length,amount,angle,angles(3)
        integer :: axis,i
        logical :: snap

        if(.not.transform_drag_matches(s))return
        if(.not.all(ieee_is_finite([x,y])))return
        axis=s%transform_drag%axis
        displacement=[x,y]-s%transform_drag%start
        if(norm2(displacement)<gizmo_slop)displacement=0._dp
        snap=s%snap_enabled.or.temporary_snap
        s%transform_drag%snapped=snap
        select case(s%transform_drag%mode)
        case(2)
            if(all(displacement==0._dp).and.s%transform_drag%spin_total==0._dp)return
            projected=[x,y]-s%transform_drag%screen
            if(norm2(projected)>9._dp)then
                angle=gizmo_angle(s%transform_drag%plane,projected)
                s%transform_drag%spin_total=s%transform_drag%spin_total+modulo(angle-s%transform_drag%spin_last+pi,2._dp*pi)-pi
                s%transform_drag%spin_last=angle
            end if
            amount=s%transform_drag%spin_total*(180._dp/pi)
            if(snap)amount=anint(amount/15._dp)*15._dp
            angles=0._dp; angles(axis)=amount
            s%transform_drag%linear=rotation_matrix(angles)
        case(1)
            s%transform_drag%offset=0._dp
            if(axis==4)then
                s%transform_drag%offset=matmul(s%transform_drag%basis,displacement)
                if(snap)s%transform_drag%offset=anint(s%transform_drag%offset/.1_dp)*.1_dp
                amount=norm2(s%transform_drag%offset)
            else
                projected=s%transform_drag%plane(:,1)
                length=sum(projected**2)*s%transform_drag%plane(1,2)
                amount=dot_product(displacement,projected)/length
                if(snap)amount=anint(amount/.1_dp)*.1_dp
                s%transform_drag%offset(axis)=amount
            end if
        case(3)
            projected=[.70710678_dp,-.70710678_dp]
            if(axis<4)projected=s%transform_drag%plane(:,1)/norm2(s%transform_drag%plane(:,1))
            amount=1._dp+dot_product(displacement,projected)/gizmo_length
            if(snap)amount=anint(amount/.1_dp)*.1_dp
            amount=max(.05_dp,min(100._dp,amount))
            s%transform_drag%linear=world_axes
            if(axis<4)then
                s%transform_drag%linear(axis,axis)=amount
            else
                do i=1,3
                    s%transform_drag%linear(i,i)=amount
                end do
            end if
        end select
        s%transform_drag%amount=amount
        call s%surface%update()
    end subroutine update_transform_drag

    subroutine finish_transform_drag(s,x,y,temporary_snap,linear,offset,pivot,changed)
        type(studio_t),intent(inout) :: s
        real(dp),intent(in) :: x,y
        logical,intent(in) :: temporary_snap
        real(dp),intent(out) :: linear(3,3),offset(3),pivot(3)
        logical,intent(out) :: changed

        changed=.false.
        if(transform_drag_matches(s))then
            call update_transform_drag(s,x,y,temporary_snap)
            linear=s%transform_drag%linear
            offset=s%transform_drag%offset
            pivot=s%transform_drag%pivot
            changed=any(abs(linear-world_axes)>1.e-12_dp).or.any(abs(offset)>1.e-12_dp)
        end if
        call cancel_transform_drag(s)
    end subroutine finish_transform_drag

    subroutine cancel_transform_drag(s)
        type(studio_t),intent(inout) :: s

        s%transform_drag=viewport_transform_t()
        call s%surface%update()
    end subroutine cancel_transform_drag

    function transform_gizmo_hint(s) result(hint)
        type(studio_t),intent(in) :: s
        character(:),allocatable :: hint,unit
        character(24) :: value
        character(6),parameter :: modes(3)=[character(6)::'Move','Rotate','Scale']
        integer :: axis,mode

        hint=''
        if(s%mesh_only)return
        mode=s%transform_mode
        if(mode<1.or.mode>3)return
        if(.not.transform_drag_matches(s))then
            hint=trim(modes(mode))//' · Drag a handle'
            return
        end if
        axis=s%transform_drag%axis
        hint=trim(modes(mode))
        if(axis<=3)hint=hint//' '//gizmo_axes(axis:axis)
        unit=''
        if(mode==2)unit='°'
        if(mode==3)unit='×'
        write(value,'(f20.3)')s%transform_drag%amount
        hint=hint//' · '//trim(adjustl(value))//unit
        if(s%transform_drag%snapped)hint=hint//' · Snap'
    end function transform_gizmo_hint

    subroutine paint_transform_gizmo(painter,s,camera,preview)
        type(QPainter),intent(inout) :: painter
        type(studio_t),intent(in) :: s
        type(viewport_camera),intent(in) :: camera
        logical,intent(in) :: preview
        integer,parameter :: steps=72
        real(dp) :: pivot(3),projected(3),center(2),axes(2,3),point(2),side(2),direction(2),base(2),path(2,steps+1)
        type(QColor) :: color
        integer :: axis,i,alpha,mode

        if(.not.transform_gizmo_available(s))return
        mode=s%transform_mode
        if(mode<1.or.mode>3)return
        pivot=transform_gizmo_pivot(s)
        if(preview)pivot=s%transform_drag%pivot+s%transform_drag%offset
        projected=project_vertex(camera,pivot)
        center=projected(:2)
        axes=gizmo_screen_axes(camera)
        call painter%setFont(QFont(family='Noto Sans',pointSize=8.5_dp,bold=.true.))
        do axis=1,3
            if(.not.gizmo_point(axes,center,mode,axis,1._dp,point))cycle
            alpha=255
            if(preview.and.axis/=s%transform_drag%axis)alpha=132
            color=QColor(gizmo_colors(1,axis),gizmo_colors(2,axis),gizmo_colors(3,axis),alpha)
            if(mode==2)then
                do i=1,steps+1
                    if(.not.gizmo_point(axes,center,mode,axis,real(i-1,dp)/steps,path(:,i)))exit
                end do
                call painter%setPen(QPen(QColor(11,20,29,165),4.4_dp))
                call painter%drawPolyline(path(1,:),path(2,:))
                call painter%setPen(QPen(color,merge(2.8_dp,1.8_dp,preview.and.axis==s%transform_drag%axis)))
                call painter%drawPolyline(path(1,:),path(2,:))
            else
                direction=(point-center)/gizmo_length
                side=[-direction(2),direction(1)]
                call painter%setPen(QPen(QColor(12,20,29,195),4.4_dp))
                call painter%drawLine(center(1),center(2),point(1),point(2))
                call painter%setPen(QPen(color,2._dp))
                call painter%drawLine(center(1),center(2),point(1),point(2))
                call painter%setPen(QPen(QColor('#16232d'),1.2_dp))
                call painter%setBrush(color)
                if(mode==1)then
                    base=point-11._dp*direction
                    call painter%drawPolygon([point(1),base(1)+5._dp*side(1),base(1)-5._dp*side(1)], &
                        [point(2),base(2)+5._dp*side(2),base(2)-5._dp*side(2)])
                else
                    call painter%drawRect(point(1)-4.5_dp,point(2)-4.5_dp,9._dp,9._dp)
                end if
                call painter%setPen(color)
                point=point+11._dp*direction
                call painter%drawText(point(1)-3._dp,point(2)+3._dp,gizmo_axes(axis:axis))
            end if
        end do
        call painter%setPen(QPen(QColor('#cad9e4'),1.5_dp))
        call painter%setBrush(QColor('#233645'))
        if(mode==3)then
            call painter%drawRect(center(1)-5.5_dp,center(2)-5.5_dp,11._dp,11._dp)
        else
            call draw_dot(painter,center,merge(5._dp,2.7_dp,mode==1))
        end if
    end subroutine paint_transform_gizmo

end module studio_viewport_transform
