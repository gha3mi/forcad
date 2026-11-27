! Project the session's native patch samples into one Qt viewport.
module studio_live_viewport
    use qt, only: QColor, QFont, QImage, QPainter, QPainterPath, QPen, QPainter_Antialiasing, Qt_AlignCenter, Qt_DashLine, Qt_DotLine, &
        Qt_NoPen, Qt_SolidLine
    use studio_state, only: studio_t, dp,mesh_layer,mesh_lattice,viewport_shaded,viewport_wireframe
    use studio_commands, only: parameter_lines_layer
    use studio_selection, only: patch_edge_boundaries,selection_control
    use studio_workspace_selection, only: control_in_selection
    use studio_physical_knots, only: physical_knot_grid, read_physical_knot_grid
    use studio_viewport_camera, only: viewport_camera, live_yaw_offset, live_pitch_offset, workspace_view_bounds, &
        workspace_view_hit, prepare_camera, rotate_for_view, project_vertex, fit_mesh_selection, fit_around_panels, &
        fit_framed_scene_in_clear_area, zoom_live_viewport, reset_live_viewport_zoom, frame_workspace_scene, live_viewport_point, &
        viewport_multisampling,viewport_navigation_active
    use studio_viewport_primitives, only: project_points, preview_transform, draw_dot, paint_tensor_lines,append_tensor_path,volume_face_front
    use studio_viewport_scene, only: projected_patch, project_scene, has_slice, show_volume_slice, &
        projected_control_patch,project_selected_controls,project_shared_controls,viewport_control_position,project_control_dependencies, &
        pick_scene_patch, pick_live_viewport,opaque_patch_depth, begin_live_control_drag, update_live_control_drag, viewport_feature_text, &
        pick_automatic, pick_face, pick_edge, pick_vertex, edge_total, projected_edge, projected_face_axes
    use studio_viewport_mesh, only: projected_mesh_part, project_visible_meshes, paint_mesh_layer_guides, &
        paint_mesh_details, pick_sample_mesh_cell=>pick_mesh_cell,mesh_covered_patches
    use studio_viewport_geometry, only: paint_viewport_geometry
    use studio_viewport_transform, only: transform_gizmo_available, transform_gizmo_pivot, transform_gizmo_point, &
        transform_gizmo_hint, begin_transform_drag, update_transform_drag, finish_transform_drag, cancel_transform_drag, &
        transform_dragging, transform_drag_matches, paint_transform_gizmo
    implicit none(type, external)
    private

    public :: save_viewport_image
    public :: paint_live_viewport, paint_viewport_hud, pick_live_viewport, live_viewport_point, live_yaw_offset, live_pitch_offset
    public :: begin_live_control_drag, update_live_control_drag
    public :: pick_scene_patch, zoom_live_viewport, reset_live_viewport_zoom, viewport_feature_text, frame_workspace_scene
    public :: workspace_view_bounds, workspace_view_hit, fit_around_panels, fit_framed_scene_in_clear_area
    public :: pick_mesh_cell, fit_mesh_selection
    public :: transform_gizmo_available, transform_gizmo_pivot, transform_gizmo_point, transform_gizmo_hint
    public :: begin_transform_drag, update_transform_drag, finish_transform_drag, cancel_transform_drag, transform_dragging
    public :: pick_automatic, pick_face, pick_edge, pick_vertex

contains

    subroutine pick_mesh_cell(s,x,y,mesh,part,cell,face)
        type(studio_t),intent(in) :: s
        real(dp),intent(in) :: x,y
        integer,intent(out) :: mesh,part,cell,face
        real(dp) :: depth

        call pick_sample_mesh_cell(s,x,y,mesh,part,cell,face,view_depth=depth)
        if(mesh<=0.or.s%mesh_only)return
        if(opaque_patch_depth(s,x,y)<=depth)return
        mesh=0; part=0; cell=0; face=0
    end subroutine pick_mesh_cell

    logical function save_viewport_image(s,path,max_dimension) result(saved)
        type(studio_t), intent(in) :: s
        character(*), intent(in) :: path
        integer, optional, intent(in) :: max_dimension
        type(QImage) :: image
        type(QPainter) :: painter
        real(dp) :: bounds(4),factor
        integer :: width,height

        saved=.false.
        bounds=workspace_view_bounds(s)
        width=nint(bounds(3)-bounds(1))
        height=nint(bounds(4)-bounds(2))
        if(width<=0.or.height<=0)return
        factor=1._dp
        if(present(max_dimension))factor=min(factor,real(max_dimension,dp)/max(width,height))
        image=QImage(max(1,nint(width*factor)),max(1,nint(height*factor)))
        if(image%isNull())return
        call image%fill(QColor(17,28,38))
        painter=QPainter(image)
        if(painter%isActive())then
            call painter%scale(factor,factor)
            call painter%translate(-bounds(1),-bounds(2))
            call paint_live_viewport(painter,s,real(s%surface%width(),dp),real(s%surface%height(),dp))
            call paint_viewport_hud(painter,s)
            call painter%end()
            saved=image%save(path)
        end if
        call image%delete()
    end function save_viewport_image

    subroutine paint_live_viewport(painter,s,width,height)
        type(QPainter), intent(inout) :: painter
        type(studio_t), intent(in) :: s
        real(dp), intent(in) :: width,height
        type(viewport_camera) :: camera
        type(projected_patch), allocatable :: objects(:)
        type(projected_mesh_part), allocatable :: mesh_parts(:)
        real(dp), allocatable :: controls(:,:),slice(:,:)
        real(dp) :: marker(3)
        integer :: counts(3), dimension, point, object
        logical,allocatable :: show_guides(:),mesh_coverage(:,:),mesh_sources_selected(:,:)
        logical :: active_visible,preview,navigating

        if(.not.prepare_camera(s,camera))return
        preview=transform_drag_matches(s)
        call project_visible_meshes(s,camera,preview,mesh_parts,clip=[0._dp,0._dp,width,height])
        navigating=dense_mesh_navigation(s,mesh_parts)
        call mesh_covered_patches(s,mesh_coverage)
        if((navigating.or.s%viewport_style/=viewport_wireframe).and.s%selected_multipatch==0.and. &
            .not.(s%physical_knots.or.s%layers(2)).and.any(mesh_coverage))then
            call project_scene(s,camera,objects,preview,clip=[0._dp,0._dp,width,height],display=.true., &
                covered=mesh_coverage,source_selected=mesh_sources_selected)
        else
            call project_scene(s,camera,objects,preview,clip=[0._dp,0._dp,width,height],display=.true.)
        end if
        dimension=s%live_scene%dimension
        counts=s%live_scene%sample_counts
        allocate(show_guides(size(objects)),source=.true.)
        if(s%mesh_section_axis>0)then
            do object=1,size(objects)
                if(objects(object)%kind/=3)cycle
                show_guides(object)=.not.mesh_coverage(objects(object)%variant,3)
            end do
        end if
        active_visible=any(objects%selected.and.objects%selection_marked.and.show_guides)

        call painter%save()
        call painter%setClipRect(0._dp,0._dp,width,height)
        if(s%world_grid)call paint_ground_grid(painter,s,camera,width)
        call painter%save()
        if(s%viewport_style==viewport_wireframe.or..not.viewport_multisampling(s))then
            call painter%setRenderHint(QPainter_Antialiasing,.false.)
        end if
        if(s%mesh_view==mesh_layer)call paint_mesh_layer_guides(painter,s,mesh_parts,navigating=navigating)
        call paint_viewport_geometry(painter,s,camera,objects,mesh_parts,nint(width),nint(height),mesh_coverage,navigating)
        do object=1,size(objects)
            if(show_guides(object))call paint_scene_lines(painter,s,objects(object))
        end do
        if(s%physical_knots.and..not.s%mesh_only)call paint_physical_knot_grid(painter,s,camera,objects,preview,show_guides)
        call paint_mesh_details(painter,s,mesh_parts,navigating=navigating,source_selected=mesh_sources_selected)
        if(.not.s%mesh_only.and.s%selected_multipatch>0)call paint_multipatch_interfaces(painter,s,objects,show_guides)
        if(allocated(s%refinement_preview_patches))call paint_refinement_preview(painter,s,objects,show_guides)
        call painter%restore()

        if (active_visible.and.dimension==3.and.has_slice(s,counts).and.show_volume_slice(s)) then
            slice=project_points(s,camera,s%live_slice,preview)
            call paint_slice(painter,s,slice,[counts(1),counts(2),1])
        end if
        if(any(s%layers(3:4)).and..not.s%mesh_only)then
            call paint_control_guides(painter,s,camera,[0._dp,0._dp,width,height],preview,mesh_coverage,controls)
        end if
        if (active_visible.and.s%live_probe%valid.and.any(s%layers(3:))) then
            marker=s%live_probe%point
            if(preview)marker=preview_transform(s%transform_drag,marker)
            marker=project_vertex(camera,marker)
            call painter%setPen(QPen(QColor(13,43,53),1.8_dp))
            call painter%setBrush(QColor(247,252,255))
            call draw_dot(painter,marker(:2),6._dp)
        end if
        if (allocated(controls).and.s%layers(4).and.s%surface_tool==1) then
            call painter%setPen(QPen(QColor('#b9fbf3'),1.6_dp))
            call painter%setBrush(QColor('#264d53'))
            if(s%selected_multipatch==0)then
                do point=1,size(controls,2)
                    if(point==s%selected_control)cycle
                    if(control_in_selection(s,point))call draw_dot(painter,controls(:2,point),5.4_dp)
                end do
            end if
            point=s%selected_control
            if(point>=1.and.point<=size(controls,2))then
                if(s%control_dragging)then
                    if(viewport_control_position(s,point,marker))then
                        marker=project_vertex(camera,marker)
                        call painter%setPen(QPen(QColor(119,218,211,160),1.3_dp,Qt_DotLine))
                        call painter%drawLine(marker(1),marker(2),controls(1,point),controls(2,point))
                    end if
                end if
                if(s%control_dragging.and.s%control_snapped)then
                    marker(1:2)=controls(:2,point)
                    call painter%setPen(QPen(QColor('#b7f7d1'),1.9_dp))
                    call painter%setBrush(QColor('#2ba980'))
                    call painter%drawPolygon(marker(1)+[0._dp,8._dp,0._dp,-8._dp], &
                        marker(2)+[-8._dp,0._dp,8._dp,0._dp])
                else
                    call painter%setPen(QPen(QColor(185,251,243),1.6_dp))
                    call painter%setBrush(QColor(39,163,160))
                    call draw_dot(painter,controls(:2,point),6._dp)
                end if
            end if
        end if
        call painter%restore()
    end subroutine paint_live_viewport

    subroutine paint_viewport_hud(painter,s)
        type(QPainter),intent(inout) :: painter
        type(studio_t),intent(in) :: s
        type(viewport_camera) :: camera
        real(dp) :: bounds(4)
        character(:),allocatable :: empty_title

        bounds=workspace_view_bounds(s)
        call painter%save()
        if(prepare_camera(s,camera))then
            call paint_triad(painter,camera,bounds)
            if(s%transform_mode>0)call paint_transform_gizmo(painter,s,camera,transform_drag_matches(s))
        else
            empty_title='No visible geometry'
            if(s%mesh_only)empty_title='No visible meshes'
            call painter%setPen(QPen(QColor('#eff9ff')))
            call painter%setFont(QFont(family='Adwaita Sans',pointSize=17._dp,bold=.true.))
            call painter%drawText(bounds(1),.5_dp*(bounds(2)+bounds(4))-15._dp,bounds(3)-bounds(1),30._dp, &
                Qt_AlignCenter,empty_title)
        end if
        call painter%restore()
    end subroutine paint_viewport_hud

    subroutine paint_control_guides(painter,s,camera,clip,preview,mesh_coverage,focused)
        type(QPainter),intent(inout) :: painter
        type(studio_t),intent(in) :: s
        type(viewport_camera),intent(in) :: camera
        real(dp),intent(in) :: clip(4)
        logical,intent(in) :: preview,mesh_coverage(:,:)
        real(dp),allocatable,intent(out) :: focused(:,:)
        type(projected_control_patch),allocatable :: patches(:)
        real(dp),allocatable :: shared(:,:)
        logical,allocatable :: shared_visible(:)
        integer,allocatable :: patch_index(:,:)
        integer :: i,active,index,point
        real(dp) :: device_scale
        logical :: motion

        call project_selected_controls(s,camera,clip,patches,preview,mesh_coverage)
        if(size(patches)==0)return
        motion=viewport_navigation_active(s)
        device_scale=s%surface%devicePixelRatio()
        if(s%layers(3))call paint_control_nets(painter,patches,motion,device_scale==1._dp)
        if(.not.s%layers(4))return
        call project_shared_controls(s,camera,patches,shared,shared_visible,preview)
        call paint_control_points(painter,patches,clip,motion,device_scale,shared,shared_visible,s%shared_controls%copies)
        if(s%control_dragging)call paint_control_dependencies(painter,s,camera,clip,shared_visible,mesh_coverage)
        if(size(shared_visible)>0)then
            if(s%surface_tool==1.and.allocated(s%selected_shared_controls))then
                call painter%setPen(QPen(QColor('#b9fbf3'),1.6_dp))
                call painter%setBrush(QColor('#264d53'))
                do i=1,size(s%selected_shared_controls)
                    point=s%selected_shared_controls(i)
                    if(point<1.or.point>size(shared_visible).or.point==s%selected_control)cycle
                    if(shared_visible(point))call draw_dot(painter,shared(:2,point),5.4_dp)
                end do
            end if
            point=s%canvas_gesture%viewport_shared_hover
            if(point>=1.and.point<=size(shared_visible).and..not.s%control_dragging.and..not.motion)then
                if(shared_visible(point))then
                    call painter%setPen(QPen(QColor('#ffcf87'),1.7_dp))
                    call painter%setBrush(QColor(0,0,0,0))
                    call draw_dot(painter,shared(:2,point),8._dp)
                end if
            end if
            i=s%selected_control
            if(i>=1.and.i<=size(shared_visible))then
                if(shared_visible(i))call move_alloc(shared,focused)
            end if
            return
        end if
        active=0
        do i=1,size(patches)
            if(.not.patches(i)%focused)cycle
            active=i
            exit
        end do
        if(s%surface_tool==1.and.allocated(s%selection%items))then
            if(any(s%selection%items%category==selection_control))then
                allocate(patch_index(maxval(patches%id),3),source=0)
                do i=1,size(patches)
                    patch_index(patches(i)%id,patches(i)%kind)=i
                end do
                call painter%setPen(QPen(QColor('#b9fbf3'),1.6_dp))
                call painter%setBrush(QColor('#264d53'))
                do i=1,size(s%selection%items)
                    associate(target=>s%selection%items(i))
                        if(target%category/=selection_control)cycle
                        if(target%rank<1.or.target%rank>3.or.target%id<1.or.target%id>size(patch_index,1))cycle
                        index=patch_index(target%id,target%rank)
                        if(index==0.or.index==active)cycle
                        if(target%item<1.or.target%item>size(patches(index)%points,2))cycle
                        call draw_dot(painter,patches(index)%points(:2,target%item),5.4_dp)
                    end associate
                end do
            end if
        end if
        if(active>0)call move_alloc(patches(active)%points,focused)
    end subroutine paint_control_guides

    subroutine paint_control_dependencies(painter,s,camera,clip,shared_visible,mesh_coverage)
        type(QPainter), intent(inout) :: painter
        type(studio_t), intent(in) :: s
        type(viewport_camera), intent(in) :: camera
        real(dp), intent(in) :: clip(4)
        logical, intent(in) :: shared_visible(:),mesh_coverage(:,:)
        real(dp), allocatable :: points(:,:)
        logical, allocatable :: visible(:)
        real(dp), parameter :: radius=6.5_dp
        integer :: i

        call project_control_dependencies(s,camera,shared_visible,mesh_coverage,points,visible)
        if(.not.any(visible))return
        call painter%save()
        call painter%setPen(QPen(QColor('#613715'),1.7_dp))
        call painter%setBrush(QColor('#ffd18d'))
        do i=1,size(visible)
            if(.not.visible(i))cycle
            if(any(points(:2,i)<clip(:2)-radius).or.any(points(:2,i)>clip(3:)+radius))cycle
            call draw_dot(painter,points(:2,i),radius)
        end do
        call painter%restore()
    end subroutine paint_control_dependencies

    subroutine paint_control_nets(painter,patches,motion,native_scale)
        type(QPainter),intent(inout) :: painter
        type(projected_control_patch),intent(in) :: patches(:)
        logical,intent(in) :: motion,native_scale
        type(QPainterPath) :: path
        type(QPen) :: pen
        real(dp) :: width
        integer :: i
        logical :: batched,axes(3)

        batched=motion.or.native_scale
        width=merge(1._dp,1.1_dp,native_scale)
        if(motion)width=0._dp
        pen=QPen(QColor(228,188,89,220),width)
        call painter%save()
        if(motion)call painter%setRenderHint(QPainter_Antialiasing,.false.)
        if(batched)then
            path=QPainterPath()
        else
            call painter%setPen(pen)
        end if
        do i=1,size(patches)
            associate(patch=>patches(i))
                axes=[.true.,patch%kind>=2,patch%kind==3]
                if(batched)then
                    call append_tensor_path(path,patch%points,patch%counts,axes)
                    if(.not.motion)call flush_guide_path(painter,path,pen)
                else
                    call paint_tensor_lines(painter,patch%points,patch%counts,axes,maxval(patch%counts),boundary_only=.false.)
                end if
            end associate
        end do
        if(batched)then
            call painter%strokePath(path,pen)
            call path%delete()
        end if
        call painter%restore()
    end subroutine paint_control_nets

    subroutine paint_control_points(painter,patches,clip,motion,device_scale,shared,shared_visible,shared_copies)
        type(QPainter),intent(inout) :: painter
        type(projected_control_patch),intent(in) :: patches(:)
        real(dp),intent(in) :: clip(4)
        logical,intent(in) :: motion
        real(dp),intent(in) :: device_scale,shared(:,:)
        logical,intent(in) :: shared_visible(:)
        integer,optional,intent(in) :: shared_copies(:)
        integer,parameter :: motion_size=9,phase_count=8
        real(dp),parameter :: radius=3.7_dp,margin=5._dp,spacing=4._dp,motion_center=.5_dp*motion_size
        type(QImage) :: markers(2),stamps(0:phase_count-1,0:phase_count-1,2)
        type(QColor) :: border(2),fill(2)
        logical,allocatable :: occupied(:,:)
        real(dp) :: lower(2),upper(2)
        integer :: i,j,k,point,cells(2),stamp_origin,stamp_size,palette

        border=[QColor(61,26,35),QColor(17,57,63)]
        fill=[QColor(248,113,111),QColor(106,226,214)]
        lower=clip(:2)-margin
        upper=clip(3:)+margin
        call painter%save()
        if(motion)then
            cells=floor((upper-lower)/spacing)
            allocate(occupied(0:cells(1),0:cells(2)),source=.false.)
        else
            stamp_origin=ceiling(5._dp*device_scale)
            stamp_size=2*stamp_origin+1
            if(device_scale/=1._dp)call painter%scale(1._dp/device_scale,1._dp/device_scale)
        end if
        do i=1,size(patches)
            associate(patch=>patches(i))
                if(allocated(patch%shared_ids))cycle
                do point=1,size(patch%points,2)
                    call paint_marker(patch%points(:2,point),patch%focused,1)
                end do
            end associate
        end do
        do point=1,size(shared_visible)
            if(.not.shared_visible(point))cycle
            palette=1
            if(present(shared_copies))then
                if(shared_copies(point)>1)palette=2
            end if
            call paint_marker(shared(:2,point),.false.,palette)
        end do
        if(motion)then
            do k=1,size(markers)
                call markers(k)%delete()
            end do
        else
            do k=1,size(stamps,3)
                do j=0,phase_count-1
                    do i=0,phase_count-1
                        call stamps(i,j,k)%delete()
                    end do
                end do
            end do
        end if
        call painter%restore()
    contains
        subroutine paint_marker(pos,focused,palette)
            real(dp),intent(in) :: pos(2)
            logical,intent(in) :: focused
            integer,intent(in) :: palette
            integer :: cell(2),phase(2),pixel(2)

            if(.not.(all(pos>=lower).and.all(pos<=upper)))return
            if(motion)then
                cell=floor((pos-lower)/spacing)
                if(occupied(cell(1),cell(2)).and..not.focused)return
                occupied(cell(1),cell(2))=.true.
                if(markers(palette)%isNull())call render_marker(markers(palette),motion_size,[motion_center,motion_center],1._dp,palette)
                call painter%drawImage(nint(pos(1)-motion_center),nint(pos(2)-motion_center),markers(palette))
            else
                pixel=nint(pos*device_scale*phase_count)
                phase=modulo(pixel,phase_count)
                pixel=(pixel-phase)/phase_count-stamp_origin
                associate(stamp=>stamps(phase(1),phase(2),palette))
                    if(stamp%isNull())call render_marker(stamp,stamp_size,stamp_origin+real(phase,dp)/phase_count,device_scale,palette)
                    call painter%drawImage(pixel(1),pixel(2),stamp)
                end associate
            end if
        end subroutine paint_marker

        subroutine render_marker(image,extent,center,scale,palette)
            type(QImage),intent(out) :: image
            integer,intent(in) :: extent,palette
            real(dp),intent(in) :: center(2),scale
            type(QPainter) :: marker_painter

            image=QImage(extent,extent)
            call image%fill(QColor(0,0,0,0))
            marker_painter=QPainter(image)
            call marker_painter%setRenderHint(QPainter_Antialiasing)
            call marker_painter%setPen(QPen(border(palette),scale))
            call marker_painter%setBrush(fill(palette))
            call draw_dot(marker_painter,center,radius*scale)
            call marker_painter%end()
        end subroutine render_marker
    end subroutine paint_control_points

    logical function dense_mesh_navigation(s,parts) result(active)
        type(studio_t),intent(in) :: s
        type(projected_mesh_part),intent(in) :: parts(:)
        real(dp),parameter :: preview_face_budget=8192._dp
        real(dp) :: span(3),faces
        integer :: i

        active=.false.
        if(size(parts)==0)return
        if(.not.viewport_navigation_active(s))return
        faces=0._dp
        do i=1,size(parts)
            associate(part=>parts(i))
                span=real(part%upper-part%lower,dp)
                select case(s%document%meshes(part%mesh)%kind)
                case(1)
                    faces=faces+span(1)
                case(2)
                    faces=faces+span(1)*span(2)
                case(3)
                    faces=faces+2._dp*(span(1)*span(2)+(span(1)+span(2))*span(3))
                    if(s%mesh_view==mesh_lattice)faces=faces+product(span+1._dp)
                case default
                    continue
                end select
            end associate
            if(faces<preview_face_budget)cycle
            active=.true.
            return
        end do
    end function dense_mesh_navigation

    subroutine paint_ground_grid(painter,s,camera,width)
        type(QPainter), intent(inout) :: painter
        type(studio_t), intent(in) :: s
        type(viewport_camera), intent(in) :: camera
        real(dp), intent(in) :: width
        real(dp) :: extent,step,origin(2),level,decade,scaled,across,ahead
        integer :: i,n

        if(camera%scale<=0._dp)return
        extent=.86_dp*width/camera%scale
        decade=10._dp**floor(log10(extent/12._dp))
        scaled=extent/(12._dp*decade)
        step=decade
        if(scaled>2._dp)step=2._dp*decade
        if(scaled>5._dp)step=5._dp*decade
        level=camera%floor
        ! Center the painted portion of the world grid under the viewport as the camera pans.
        across=-s%pan_x/camera%scale
        ahead=0._dp
        if(abs(camera%pitch_sin)>.15_dp) then
          ahead=(s%pan_y/camera%scale-camera%pitch_cos*(level-camera%center(3)))/camera%pitch_sin
        end if
        origin=camera%center(:2)+[camera%yaw_cos*across+camera%yaw_sin*ahead, &
            -camera%yaw_sin*across+camera%yaw_cos*ahead]
        origin=anint(origin/step)*step
        n=min(32,ceiling(extent/step))
        extent=real(n,dp)*step
        do i=-n,n
            call painter%setPen(QPen(QColor(100,129,150,28),.75_dp))
            if(mod(i,5)==0)call painter%setPen(QPen(QColor(112,144,168,51),.9_dp))
            call draw_grid_line(1,origin(2)+real(i,dp)*step)
            call draw_grid_line(2,origin(1)+real(i,dp)*step)
        end do
        call painter%setPen(QPen(QColor(192,109,104,84),1._dp))
        if(abs(origin(2))<=extent)call draw_grid_line(1,0._dp)
        call painter%setPen(QPen(QColor(100,169,137,86),1._dp))
        if(abs(origin(1))<=extent)call draw_grid_line(2,0._dp)

    contains

        subroutine draw_grid_line(axis,value)
            integer, intent(in) :: axis
            real(dp), intent(in) :: value
            real(dp) :: first(3),last(3),a(3),b(3)
            integer :: fixed

            first=[origin,level]; last=first
            fixed=3-axis
            first(fixed)=value; last(fixed)=value
            first(axis)=origin(axis)-extent
            last(axis)=origin(axis)+extent
            a=project_vertex(camera,first); b=project_vertex(camera,last)
            call painter%drawLine(a(1),a(2),b(1),b(2))
        end subroutine draw_grid_line

    end subroutine paint_ground_grid

    subroutine paint_physical_knot_grid(painter,s,camera,objects,moving,show_guides)
        type(QPainter), intent(inout) :: painter
        type(studio_t), intent(in) :: s
        type(viewport_camera), intent(in) :: camera
        type(projected_patch), intent(in) :: objects(:)
        logical, intent(in) :: moving,show_guides(:)
        type(physical_knot_grid), pointer :: grid
        real(dp), allocatable :: points(:,:)
        real(dp) :: width
        logical :: skin,front(0:6),navigating,batched
        integer :: i,line,point,side,pen_kind,node
        type(QPen) :: line_pen,marker_pen
        type(QPainterPath) :: path

        navigating=viewport_navigation_active(s)
        batched=.not.navigating.and.s%surface%devicePixelRatio()==1._dp
        width=merge(1._dp,1.55_dp,batched)
        if(navigating)width=0._dp
        line_pen=QPen(QColor(151,250,205,235),width)
        if(batched)path=QPainterPath()
        marker_pen=QPen(QColor('#143940'),1.35_dp)
        pen_kind=0
        call painter%save()
        call painter%setBrush(QColor('#9af4c9'))
        do i=1,size(objects)
            if(.not.show_guides(i))cycle
            associate(object=>objects(i))
                if(.not.(object%selection_marked.or.object%group_member).or.object%is_point)cycle
                call read_physical_knot_grid(object%kind,object%variant,grid,catalog=s%document%patches,cached_only=navigating)
                if(.not.associated(grid))cycle
                points=project_points(s,camera,real(grid%points,dp), &
                    moving.and.(object%selected.or.(s%selected_multipatch>0.and.object%group_member)))
                if(object%kind==1)then
                    if(pen_kind/=1)call painter%setPen(marker_pen)
                    pen_kind=1
                    do point=1,size(points,2)
                        call draw_dot(painter,points(:2,point),4.7_dp)
                    end do
                    cycle
                end if
                skin=object%kind==3.and.s%layers(1).and.object%opacity==255
                if(skin)then
                    front(0)=.false.
                    front(1:)=[(volume_face_front(object%points,object%counts,side),side=1,6)]
                end if
                if(.not.batched.and.pen_kind/=2)call painter%setPen(line_pen)
                pen_kind=2
                do line=1,size(grid%lines,1)
                    if(skin)then
                        if(.not.any(front(grid%faces(:,line))))cycle
                    end if
                    if(.not.batched)then
                        call painter%drawPolyline(points(1,grid%lines(line,:)),points(2,grid%lines(line,:)))
                        cycle
                    end if
                    node=grid%lines(line,1)
                    call path%moveTo(points(1,node),points(2,node))
                    do point=2,size(grid%lines,2)
                        node=grid%lines(line,point)
                        call path%lineTo(points(1,node),points(2,node))
                    end do
                end do
            end associate
            if(batched)call flush_guide_path(painter,path,line_pen)
        end do
        if(batched)then
            call painter%strokePath(path,line_pen)
            call path%delete()
        end if
        call painter%restore()
    end subroutine paint_physical_knot_grid

    subroutine flush_guide_path(painter,path,pen)
        type(QPainter), intent(inout) :: painter
        type(QPainterPath), intent(inout) :: path
        type(QPen), intent(in) :: pen
        integer, parameter :: batch_elements=16384

        if(path%elementCount()<batch_elements)return
        call painter%strokePath(path,pen)
        call path%delete()
        path=QPainterPath()
    end subroutine flush_guide_path

    subroutine paint_scene_lines(painter,s,object)
        type(QPainter), intent(inout) :: painter
        type(studio_t), intent(in) :: s
        type(projected_patch), intent(in) :: object
        real(dp), allocatable :: points(:,:)
        logical :: front(0:6),wireframe,guides
        integer :: side,limit,axes(2),counts(2),edge,sides(2)

        guides=(object%selection_marked.or.object%group_member).and.s%layers(parameter_lines_layer)
        if(object%kind==1)then
            if(guides.and..not.object%is_point)call paint_curve_parameter_marks(painter,object%points)
            return
        end if
        wireframe=s%layers(2)
        if(.not.(wireframe.or.guides))return
        limit=merge(8,5,wireframe)
        call painter%setPen(QPen(QColor(124,190,212,merge(109,171,wireframe)),.85_dp))
        front(0)=.false.
        front(1:)=.true.
        if(object%kind==3.and..not.wireframe.and.s%layers(1).and.object%opacity==255)then
            front(1:)=[(volume_face_front(object%points,object%counts,side),side=1,6)]
        end if
        if(object%kind==2)then
            call paint_face_isocurves(painter,object%points,object%counts(:2),limit)
        else
            do side=1,6
                if(.not.front(side))cycle
                call projected_face_axes(object,side,axes,counts)
                call paint_face_isocurves(painter,object%boundary(side)%points,counts,limit)
            end do
        end if
        if(s%viewport_style/=viewport_shaded)return
        do edge=1,edge_total(object%kind)
            sides=patch_edge_boundaries(object%kind,edge)
            if(.not.any(front(sides)))cycle
            points=projected_edge(object,edge,sides)
            call painter%drawPolyline(points(1,:),points(2,:))
        end do
    end subroutine paint_scene_lines

    subroutine paint_face_isocurves(painter,points,counts,limit)
        type(QPainter), intent(inout) :: painter
        real(dp), intent(in) :: points(:,:)
        integer, intent(in) :: counts(2),limit
        integer :: axis,across,steps(2),lines,line,first,last

        steps=[1,counts(1)]
        do axis=1,2
            across=3-axis
            lines=min(limit,counts(across)-1)
            if(counts(axis)<2)cycle
            do line=1,lines-1
                first=1+(line*(counts(across)-1)/lines)*steps(across)
                last=first+(counts(axis)-1)*steps(axis)
                call painter%drawPolyline(points(1,first:last:steps(axis)),points(2,first:last:steps(axis)))
            end do
        end do
    end subroutine paint_face_isocurves

    subroutine paint_refinement_preview(painter,s,objects,show_guides)
        type(QPainter),intent(inout) :: painter
        type(studio_t),intent(in) :: s
        type(projected_patch),intent(in) :: objects(:)
        logical,intent(in) :: show_guides(:)
        real(dp),allocatable :: points(:,:)
        logical :: front(0:6)
        integer :: i,edge,side,sides(2)

        call painter%save()
        call painter%setPen(QPen(QColor('#ffbf64'),2.8_dp))
        do i=1,size(objects)
            if(.not.show_guides(i))cycle
            associate(object=>objects(i))
                if(object%kind/=s%refinement_preview_kind)cycle
                if(object%variant>size(s%refinement_preview_patches))cycle
                if(.not.s%refinement_preview_patches(object%variant))cycle
                if(object%kind==1)then
                    call painter%drawPolyline(object%points(1,:),object%points(2,:))
                    cycle
                end if
                front(0)=.false.
                front(1:)=.true.
                if(object%kind==3.and.s%layers(1).and.object%opacity==255)then
                    front(1:)=[(volume_face_front(object%points,object%counts,side),side=1,6)]
                end if
                do edge=1,edge_total(object%kind)
                    sides=patch_edge_boundaries(object%kind,edge)
                    if(.not.any(front(sides)))cycle
                    points=projected_edge(object,edge,sides)
                    call painter%drawPolyline(points(1,:),points(2,:))
                end do
            end associate
        end do
        call painter%restore()
    end subroutine paint_refinement_preview

    subroutine paint_multipatch_interfaces(painter,s,objects,show_guides)
        type(QPainter), intent(inout) :: painter
        type(studio_t), intent(in) :: s
        type(projected_patch), intent(in) :: objects(:)
        logical,intent(in) :: show_guides(:)
        type(QColor) :: ink
        integer,allocatable :: visible(:)
        integer :: join,object,first,second,kind
        logical :: valid,unconstrained

        associate(group=>s%document%multipatches%groups(s%selected_multipatch))
            if(.not.allocated(group%definition%joins).or..not.allocated(group%validation%joins))return
            if(size(group%definition%joins)==0.or.size(objects)==0)return
            kind=group%definition%kind
            if(.not.allocated(s%document%scene%catalog(kind)%entries))return
            allocate(visible(size(s%document%scene%catalog(kind)%entries)),source=0)
            do object=1,size(objects)
                if(.not.(show_guides(object).and.objects(object)%group_member))cycle
                if(objects(object)%kind/=kind)cycle
                visible(objects(object)%variant)=object
            end do
            do join=1,min(size(group%definition%joins),size(group%validation%joins))
                associate(connection=>group%definition%joins(join))
                    first=0; second=0
                    if(connection%patch_a>0.and.connection%patch_a<=size(visible))first=visible(connection%patch_a)
                    if(connection%patch_b>0.and.connection%patch_b<=size(visible))second=visible(connection%patch_b)
                    if(first==0.and.second==0)cycle
                    valid=group%validation%joins(join)%valid
                    unconstrained=connection%continuity<0
                    ink=QColor('#78efbf')
                    if(unconstrained)ink=QColor('#9dcaf1')
                    if(.not.valid)ink=QColor('#ffc092')
                    call painter%setPen(QPen(ink,2.5_dp,merge(Qt_DashLine,Qt_SolidLine,.not.valid.or.unconstrained)))
                    call painter%setBrush(ink)
                    if(first>0)call paint_interface_boundary(painter,objects(first),connection%side_a)
                    if(second>0.and.(.not.valid.or.unconstrained))call paint_interface_boundary(painter,objects(second),connection%side_b)
                end associate
            end do
        end associate
    end subroutine paint_multipatch_interfaces

    subroutine paint_interface_boundary(painter,object,side)
        type(QPainter), intent(inout) :: painter
        type(projected_patch), intent(in) :: object
        integer, intent(in) :: side
        integer :: endpoint,edge,sides(2)
        real(dp), allocatable :: points(:,:)

        if(side<1.or.side>2*object%kind)return
        if(object%kind==1)then
            endpoint=merge(1,size(object%points,2),side==1)
            call draw_dot(painter,object%points(:2,endpoint),4.5_dp)
            return
        end if
        do edge=1,edge_total(object%kind)
            if(.not.any(patch_edge_boundaries(object%kind,edge)==side))cycle
            points=projected_edge(object,edge,sides)
            call painter%drawPolyline(points(1,:),points(2,:))
        end do
    end subroutine paint_interface_boundary

    subroutine paint_curve_parameter_marks(painter,path)
        type(QPainter), intent(inout) :: painter
        real(dp), intent(in) :: path(:,:)
        real(dp) :: arrow(2),tail(2),wing(2),direction(2),length
        integer :: point,step

        step=max(1,(size(path,2)-1)/10)
        call painter%setPen(QPen(QColor(18,58,73),1._dp))
        call painter%setBrush(QColor(171,226,240))
        do point=1,size(path,2),step
            call draw_dot(painter,path(:2,point),2.3_dp)
        end do

        if(size(path,2)<2)return
        point=max(1,2*(size(path,2)-1)/3)
        arrow=path(:2,point+1)
        direction=arrow-path(:2,point)
        length=norm2(direction)
        if(length<.1_dp)return
        direction=direction/length
        tail=arrow-8._dp*direction
        wing=4._dp*[-direction(2),direction(1)]
        call painter%setPen(QPen(QColor(230,250,253),1.55_dp))
        call painter%drawLine(arrow(1),arrow(2),tail(1)+wing(1),tail(2)+wing(2))
        call painter%drawLine(arrow(1),arrow(2),tail(1)-wing(1),tail(2)-wing(2))
    end subroutine paint_curve_parameter_marks

    subroutine paint_slice(painter,s,slice,counts)
        type(QPainter), intent(inout) :: painter
        type(studio_t), intent(in) :: s
        real(dp), intent(in) :: slice(:,:)
        integer, intent(in) :: counts(3)
        integer :: i,j,base,nu,nv,vertices(4)
        if (.not.any(s%layers([1,2,parameter_lines_layer]))) return
        nu=counts(1); nv=counts(2)
        if(s%viewport_style/=viewport_wireframe)then
            call painter%save()
            call painter%setRenderHint(QPainter_Antialiasing,.false.)
            call painter%setPen(QPen(style=Qt_NoPen))
            call painter%setBrush(QColor(118,236,229,70))
            do j=1,nv-1
                do i=1,nu-1
                    base=(j-1)*nu+i
                    vertices=base+[0,1,nu+1,nu]
                    call painter%drawPolygon(slice(1,vertices),slice(2,vertices))
                end do
            end do
            call painter%restore()
        end if
        if (s%layers(parameter_lines_layer)) then
            call painter%setPen(QPen(QColor(135,244,234,184),.9_dp))
            call paint_tensor_lines(painter,slice,counts,[.true.,.true.,.false.],4,boundary_only=.false.)
        end if
        call painter%setPen(QPen(QColor(156,252,237),1.7_dp))
        call paint_tensor_lines(painter,slice,counts,[.true.,.true.,.false.],1,boundary_only=.false.)
    end subroutine paint_slice

    subroutine paint_triad(painter,camera,bounds)
        type(QPainter), intent(inout) :: painter
        type(viewport_camera), intent(in) :: camera
        real(dp), intent(in) :: bounds(4)
        real(dp), parameter :: unit(3,3)=reshape([1._dp,0._dp,0._dp, 0._dp,1._dp,0._dp, 0._dp,0._dp,1._dp],[3,3])
        integer, parameter :: colors(3,3)=reshape([241,122,121, 123,175,244, 110,210,156],[3,3])
        character(len=3), parameter :: names='XYZ'
        real(dp) :: origin(2),direction(3),tip(2),label(2)
        integer :: axis

        origin=[bounds(1)+33._dp,bounds(4)-25._dp]
        call painter%setFont(QFont(family='Noto Sans',pointSize=8.5_dp))
        do axis=1,3
            direction=rotate_for_view(camera,unit(:,axis))
            if (norm2(direction(:2))<.16_dp) cycle
            tip=origin+22._dp*direction(:2)
            label=origin+30._dp*direction(:2)
            call painter%setPen(QPen(QColor(colors(1,axis),colors(2,axis),colors(3,axis)),1.65_dp))
            call painter%drawLine(origin(1),origin(2),tip(1),tip(2))
            call painter%drawText(label(1)-3._dp,label(2)+3._dp,names(axis:axis))
        end do
        call painter%setPen(QPen(style=Qt_NoPen))
        call painter%setBrush(QColor(224,235,242))
        call draw_dot(painter,origin,2._dp)
    end subroutine paint_triad

end module studio_live_viewport
