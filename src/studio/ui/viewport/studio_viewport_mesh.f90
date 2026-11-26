module studio_viewport_mesh
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    !$ use omp_lib, only: omp_get_max_threads
    use qt, only: QBrush, QColor, QPainter, QPen, Qt_DashLine, Qt_NoBrush
    use studio_selection, only: selection_mesh
    use studio_workspace_selection, only: mesh_in_selection, mesh_cell_in_selection
    use studio_state, only: studio_t, dp, studio_mesh_count, mesh_exterior, mesh_lattice, mesh_layer, viewport_wireframe
    use studio_mesh, only: studio_mesh_part
    use studio_multipatch_collection, only: multipatch_contains
    use studio_viewport_camera, only: viewport_camera, prepare_camera, workspace_view_hit, mesh_drawable, mesh_volume_limits, &
        rotate_for_view,project_vertex,project_viewport_box,viewport_clip_intersects
    use studio_viewport_primitives, only: screen_quad, project_points, draw_dot, paint_tensor_lines, volume_face_front, &
        projected_face_light, nearest_curve_parameter, screen_triangle, preview_transform
    use studio_mesh_section, only: mesh_section_face, build_mesh_section, mesh_boundary_cells, mesh_hex_faces
    implicit none(type, external)
    private

    public :: projected_mesh_part, project_visible_meshes, mesh_replaces_patch_fill, paint_mesh_layer_guides
    public :: mesh_covered_patches
    public :: paint_mesh_details, pick_mesh_cell, add_mesh_faces, mesh_face_color, mesh_section_face, mesh_section_color

    real(dp), parameter :: mesh_line_pick_radius=7._dp, straight_grid_tolerance=1.e-8_dp, straight_preview_tolerance=.6_dp
    ! Project cyclic face outlines from ForCAD's direction-1-fastest cell rows.
    integer, parameter :: mesh_quad_nodes(4)=[1,2,4,3]

    type :: projected_mesh_part
        real(dp), allocatable :: points(:,:)
        type(mesh_section_face), allocatable :: section_faces(:)
        integer :: mesh, part, counts(3), lower(3), upper(3)
    end type projected_mesh_part

contains

    subroutine project_visible_meshes(s,camera,moving,parts,clip)
        type(studio_t), intent(in) :: s
        type(viewport_camera), intent(in) :: camera
        logical, intent(in) :: moving
        type(projected_mesh_part), allocatable, intent(out) :: parts(:)
        real(dp), intent(in), optional :: clip(4)
        integer :: id,used,part_count,section_axis,i,part,counts(3),lower(3),upper(3),point,dimensions,y_column
        logical,allocatable :: visible(:)
        integer,allocatable :: targets(:,:)
        real(dp),allocatable :: world(:,:)
        real(dp) :: section_coordinate,section_normal(3),position(3)
        logical :: section_front,transform

        call mesh_section_position(s,section_axis,section_coordinate)
        section_front=.false.
        if(section_axis>0)then
            section_normal=0._dp
            section_normal(section_axis)=1._dp
            section_normal=rotate_for_view(camera,section_normal)
            section_front=section_normal(3)>0._dp
        end if
        part_count=0
        allocate(visible(studio_mesh_count(s)))
        do id=1,studio_mesh_count(s)
            visible(id)=mesh_drawable(s,id)
            if(visible(id))part_count=part_count+size(s%document%meshes(id)%parts)
        end do
        allocate(targets(2,part_count))
        used=0
        do id=1,studio_mesh_count(s)
            if(id/=s%selected_mesh)call collect_parts(id)
        end do
        if(s%selected_mesh>0)call collect_parts(s%selected_mesh)
        allocate(parts(used))
        !$omp parallel do default(none) schedule(static) if(used>=256.and.section_axis==0) num_threads(min(4,omp_get_max_threads())) &
        !$omp shared(s,camera,moving,used,targets,parts,section_axis,section_coordinate,section_front) &
        !$omp private(id,part,counts,lower,upper,world,transform,dimensions,y_column,point,position)
        do i=1,used
            id=targets(1,i)
            part=targets(2,i)
            associate(mesh=>s%document%meshes(id),piece=>s%document%meshes(id)%parts(part),projected=>parts(i))
                counts=piece%samples
                call mesh_volume_limits(s,counts,lower,upper)
                projected%mesh=id
                projected%part=part
                projected%counts=counts
                projected%lower=lower
                projected%upper=upper

                transform=.false.
                if(moving)then
                    if(s%selected_multipatch>0)then
                        transform=multipatch_contains(s%document%multipatches,s%selected_multipatch,mesh%kind,piece%variant)
                    else
                        transform=s%live_scene%dimension==mesh%kind.and.s%live_scene%variant==piece%variant
                    end if
                end if
                if(section_axis>0.and.mesh%kind==3)then
                    world=mesh_world_points(s,piece,transform)
                    projected%points=project_points(s,camera,world)
                    call build_mesh_section(world,projected%points,piece%cells,counts,lower,upper, &
                        section_axis,section_coordinate,section_front,projected%section_faces)
                else
                    dimensions=size(piece%points,2)
                    y_column=min(dimensions,2)
                    allocate(projected%points(3,size(piece%points,1)))
                    do point=1,size(piece%points,1)
                        position=[real(piece%points(point,1),dp), &
                            merge(real(piece%points(point,y_column),dp),0._dp,dimensions>=2), &
                            merge(real(piece%points(point,dimensions),dp),0._dp,dimensions==3)]
                        if(transform)position=preview_transform(s%transform_drag,position)
                        projected%points(:,point)=project_vertex(camera,position)
                    end do
                end if
            end associate
        end do
        !$omp end parallel do
    contains
        subroutine collect_parts(mesh_id)
            integer, intent(in) :: mesh_id
            integer :: index

            if(mesh_id>size(visible))return
            if(.not.visible(mesh_id))return
            associate(mesh=>s%document%meshes(mesh_id))
                do index=1,size(mesh%parts)
                    if(present(clip).and..not.moving)then
                        if(.not.mesh_part_intersects(camera,mesh%parts(index),clip))cycle
                    end if
                    used=used+1
                    targets(:,used)=[mesh_id,index]
                end do
            end associate
        end subroutine collect_parts
    end subroutine project_visible_meshes

    logical function mesh_part_intersects(camera,part,clip) result(intersects)
        type(viewport_camera), intent(in) :: camera
        type(studio_mesh_part), intent(in) :: part
        real(dp), intent(in) :: clip(4)
        real(dp) :: lower(3),upper(3)
        integer :: axis,n

        lower=0._dp; upper=0._dp
        n=size(part%points,2)
        do axis=1,n
            lower(axis)=minval(part%points(:,axis))
            upper(axis)=maxval(part%points(:,axis))
        end do
        intersects=viewport_clip_intersects(project_viewport_box(camera,lower,upper),clip)
    end function mesh_part_intersects

    subroutine mesh_section_position(s,axis,coordinate)
        type(studio_t), intent(in) :: s
        integer, intent(out) :: axis
        real(dp), intent(out) :: coordinate
        real(dp) :: lower,upper,fraction
        integer :: mesh,part,direction

        axis=0
        coordinate=0._dp
        direction=s%mesh_section_axis
        if(direction<1.or.direction>3)return
        if(.not.ieee_is_finite(s%mesh_section_fraction))return
        lower=huge(1._dp)
        upper=-huge(1._dp)
        do mesh=1,studio_mesh_count(s)
            if(.not.mesh_drawable(s,mesh))cycle
            if(s%document%meshes(mesh)%kind/=3)cycle
            do part=1,size(s%document%meshes(mesh)%parts)
                associate(points=>s%document%meshes(mesh)%parts(part)%points)
                    if(size(points,2)<direction)cycle
                    lower=min(lower,real(minval(points(:,direction)),dp))
                    upper=max(upper,real(maxval(points(:,direction)),dp))
                end associate
            end do
        end do
        if(upper<=lower.or..not.ieee_is_finite(lower).or..not.ieee_is_finite(upper))return
        fraction=max(0._dp,min(1._dp,s%mesh_section_fraction))
        coordinate=(1._dp-fraction)*lower+fraction*upper
        axis=direction
    end subroutine mesh_section_position

    logical function mesh_replaces_patch_fill(s,kind,variant) result(replaces)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: kind,variant
        integer :: id

        replaces=.false.
        if(kind==1)return
        if(kind==2.and..not.s%mesh_shaded)return
        if(kind==3.and..not.s%mesh_shaded.and.s%mesh_view==mesh_exterior.and.s%mesh_section_axis==0)return
        do id=1,studio_mesh_count(s)
            if(s%document%meshes(id)%kind/=kind)cycle
            if(.not.mesh_drawable(s,id))cycle
            if(.not.any(s%document%meshes(id)%parts%variant==variant))cycle
            replaces=.true.
            return
        end do
    end function mesh_replaces_patch_fill

    subroutine mesh_covered_patches(s,covered)
        type(studio_t), intent(in) :: s
        logical, allocatable, intent(out) :: covered(:,:)
        integer :: id,kind,maximum,part,variant

        maximum=0
        do kind=1,size(s%document%scene%catalog)
            if(allocated(s%document%scene%catalog(kind)%entries))then
                maximum=max(maximum,size(s%document%scene%catalog(kind)%entries))
            end if
        end do
        allocate(covered(maximum,3),source=.false.)
        do id=1,studio_mesh_count(s)
            associate(mesh=>s%document%meshes(id))
                kind=mesh%kind
                if(kind<2.or.kind>3)cycle
                if(kind==2.and..not.s%mesh_shaded)cycle
                if(kind==3.and..not.s%mesh_shaded.and.s%mesh_view==mesh_exterior.and.s%mesh_section_axis==0)cycle
                if(.not.mesh_drawable(s,id))cycle
                do part=1,size(mesh%parts)
                    variant=mesh%parts(part)%variant
                    if(variant>=1.and.variant<=maximum)covered(variant,kind)=.true.
                end do
            end associate
        end do
    end subroutine mesh_covered_patches

    function mesh_world_points(s,part,transform) result(world)
        type(studio_t), intent(in) :: s
        type(studio_mesh_part), intent(in) :: part
        logical, intent(in) :: transform
        real(dp), allocatable :: world(:,:)
        integer :: point

        allocate(world(3,size(part%points,1)),source=0._dp)
        world(:size(part%points,2),:)=transpose(real(part%points,dp))
        if(transform)then
            do point=1,size(world,2)
                world(:,point)=preview_transform(s%transform_drag,world(:,point))
            end do
        end if
    end function mesh_world_points

    subroutine paint_mesh_layer_guides(painter,s,parts,navigating)
        type(QPainter), intent(inout) :: painter
        type(studio_t), intent(in) :: s
        type(projected_mesh_part), intent(in) :: parts(:)
        logical, optional, intent(in) :: navigating
        integer :: i

        if(size(parts)==0.or.s%mesh_section_axis/=0)return
        if(present(navigating))then
            if(navigating)return
        end if
        call painter%save()
        call painter%setPen(QPen(QColor(124,173,189,115),.9_dp,Qt_DashLine))
        do i=1,size(parts)
            associate(part=>parts(i))
                if(s%document%meshes(part%mesh)%kind/=3)cycle
                if(all(part%lower==0).and.all(part%upper==part%counts-1))cycle
                call paint_tensor_lines(painter,part%points,part%counts,[.true.,.true.,.true.],1,boundary_only=.true.)
            end associate
        end do
        call painter%restore()
    end subroutine paint_mesh_layer_guides

    subroutine paint_mesh_details(painter,s,parts,navigating,source_selected)
        type(QPainter), intent(inout) :: painter
        type(studio_t), intent(in) :: s
        type(projected_mesh_part), intent(in) :: parts(:)
        logical, optional, intent(in) :: navigating
        logical, optional, intent(in) :: source_selected(:,:)
        type(QPen) :: wire,marker,shared_edge
        integer :: i,node,kind,alpha,cell,j,mesh_id,variant
        logical :: selected,wireframe,preview,new_mesh,source_marked,cyan

        wireframe=s%viewport_style==viewport_wireframe
        preview=.false.
        if(present(navigating))preview=navigating
        selected=.false.
        alpha=190
        mesh_id=0
        cyan=.false.
        call painter%save()
        do i=1,size(parts)
            associate(part=>parts(i),mesh=>s%document%meshes(parts(i)%mesh))
                if(allocated(part%section_faces))cycle
                kind=mesh%kind
                new_mesh=part%mesh/=mesh_id
                if(new_mesh)then
                    mesh_id=part%mesh
                    selected=mesh_in_selection(s,mesh_id)
                    alpha=merge(245,190,selected)
                    wire=QPen(QColor(255,213,139,alpha),merge(1.5_dp,1._dp,selected))
                    shared_edge=QPen(QColor(255,213,139,alpha+(alpha*(255-alpha)+127)/255),wire%width)
                end if
                if(preview)then
                    source_marked=.false.
                    if(present(source_selected))then
                        variant=mesh%parts(part%part)%variant
                        if(variant>=1.and.variant<=size(source_selected,1).and.kind>=1.and.kind<=size(source_selected,2))then
                            source_marked=source_selected(variant,kind)
                        end if
                    end if
                    if(new_mesh.or.(source_marked.neqv.cyan))then
                        cyan=source_marked
                        if(cyan)then
                            call painter%setPen(QPen(QColor(130,235,255,215),0._dp))
                        else
                            call painter%setPen(QPen(QColor(255,213,139,alpha),0._dp))
                        end if
                    end if
                end if
                associate(points=>part%points,counts=>part%counts,lower=>part%lower,upper=>part%upper)
                    if(kind==3.and.s%mesh_view==mesh_lattice.and..not.preview)then
                        call painter%setPen(QPen(QColor(228,168,92,merge(103,70,selected)),.8_dp))
                        call paint_tensor_lines(painter,points,counts,[.true.,.true.,.true.], &
                            maxval(counts),boundary_only=.false.)
                    end if
                    if(s%mesh_shaded.and.kind>1.and..not.wireframe)cycle
                    if(preview)then
                        if(kind==3)then
                            call paint_mesh_perimeter(painter,points,counts,lower,upper,kind,straight_preview_tolerance)
                        else
                            call paint_mesh_perimeter(painter,points,counts,[0,0,0],counts-1,kind,straight_preview_tolerance)
                        end if
                    else
                        if(wireframe.and.kind==3)then
                            call painter%setPen(shared_edge)
                            call paint_mesh_perimeter(painter,points,counts,lower,upper,kind,straight_grid_tolerance)
                            call painter%setPen(wire)
                        else if(new_mesh.or.kind==1.or.s%mesh_view==mesh_lattice)then
                            call painter%setPen(wire)
                        end if
                        select case(kind)
                        case(1)
                            call painter%drawPolyline(points(1,:),points(2,:))
                            marker=QPen(QColor(255,231,181,alpha),1.3_dp)
                            call painter%setPen(marker)
                            call painter%setBrush(QColor('#172732'))
                            do node=1,size(points,2)
                                call draw_dot(painter,points(:2,node),merge(3.2_dp,2.7_dp,selected))
                            end do
                        case(2)
                            call paint_mesh_grid(painter,points,counts,[0,0,0],counts-1,kind,wireframe)
                        case(3)
                            call paint_mesh_grid(painter,points,counts,lower,upper,kind,wireframe)
                        end select
                    end if
                    if(all(s%canvas_gesture%viewport_mesh_hover(:2)==[part%mesh,part%part]))then
                        cell=s%canvas_gesture%viewport_mesh_hover(3)
                        if(cell>=1.and.cell<=size(mesh%parts(part%part)%cells,1))then
                            call paint_mesh_cell(painter,points,counts,mesh%parts(part%part)%cells(cell,:), &
                                kind,lower,upper,wireframe,hovered=.true.)
                        end if
                    end if
                    if(.not.selected)cycle
                    if(allocated(s%selection%items))then
                        do j=1,s%selection%count()
                            associate(target=>s%selection%items(j))
                                if(target%category/=selection_mesh.or.target%id/=part%mesh.or.target%part/=part%part)cycle
                                cell=target%item
                                if(cell<1.or.cell>size(mesh%parts(part%part)%cells,1))cycle
                                call paint_mesh_cell(painter,points,counts,mesh%parts(part%part)%cells(cell,:), &
                                    kind,lower,upper,wireframe)
                            end associate
                        end do
                    else if(part%part==s%mesh_cell(1))then
                        cell=s%mesh_cell(2)
                        if(cell<1.or.cell>size(mesh%parts(part%part)%cells,1))cycle
                        call paint_mesh_cell(painter,points,counts,mesh%parts(part%part)%cells(cell,:),kind,lower,upper,wireframe)
                    end if
                end associate
            end associate
        end do
        call painter%restore()
    end subroutine paint_mesh_details

    subroutine paint_mesh_cell(painter,points,counts,nodes,kind,lower,upper,wireframe,hovered)
        type(QPainter), intent(inout) :: painter
        real(dp), intent(in) :: points(:,:)
        integer, intent(in) :: counts(3),nodes(:),kind,lower(3),upper(3)
        logical, intent(in) :: wireframe
        logical, optional, intent(in) :: hovered
        integer :: side,face(4)
        logical :: hover
        type(QColor) :: outline

        hover=.false.
        if(present(hovered))hover=hovered
        outline=QColor('#fde2fc')
        if(hover)outline=QColor('#c4e6f4')
        call painter%save()
        call painter%setPen(QPen(outline,merge(1.5_dp,2.5_dp,hover)))
        if(wireframe.or.hover)then
            call painter%setBrush(QBrush(style=Qt_NoBrush))
        else
            call painter%setBrush(QColor(204,98,222,164))
        end if
        select case(kind)
        case(1)
            if(.not.hover)then
                call painter%setPen(QPen(QColor(146,77,187,238),7._dp))
                call painter%drawLine(points(1,nodes(1)),points(2,nodes(1)),points(1,nodes(2)),points(2,nodes(2)))
                call painter%setPen(QPen(outline,3._dp))
            end if
            call painter%drawLine(points(1,nodes(1)),points(2,nodes(1)),points(1,nodes(2)),points(2,nodes(2)))
        case(2)
            face=nodes(mesh_quad_nodes)
            call painter%drawPolygon(points(1,face),points(2,face))
        case(3)
            do side=1,6
                if(wireframe)then
                    face=nodes(mesh_hex_faces(:,side))
                else
                    if(.not.mesh_boundary_face(points,counts,nodes,side,lower,upper,face))cycle
                end if
                call painter%drawPolygon(points(1,face),points(2,face))
            end do
        end select
        call painter%restore()
    end subroutine paint_mesh_cell

    subroutine paint_mesh_perimeter(painter,points,counts,lower,upper,kind,tolerance)
        type(QPainter), intent(inout) :: painter
        real(dp), intent(in) :: points(:,:)
        integer, intent(in) :: counts(3),lower(3),upper(3),kind
        real(dp), intent(in) :: tolerance
        integer :: origin,stride(3),span(3),axis,across,a,b,corner,first,last

        stride=[1,counts(1),product(counts(:2))]
        span=upper-lower
        origin=1+sum(lower*stride)
        do axis=1,kind
            if(span(axis)<1)cycle
            do corner=0,2**(kind-1)-1
                first=origin
                b=0
                do across=1,kind
                    if(across==axis)cycle
                    a=merge(span(across),0,btest(corner,b))
                    first=first+a*stride(across)
                    b=b+1
                end do
                last=first+span(axis)*stride(axis)
                call paint_mesh_stroke(painter,points,first,last,stride(axis),tolerance)
            end do
        end do
    end subroutine paint_mesh_perimeter

    subroutine paint_mesh_stroke(painter,points,first,last,stride,tolerance)
        type(QPainter), intent(inout) :: painter
        real(dp), intent(in) :: points(:,:),tolerance
        integer, intent(in) :: first,last,stride

        if(straight_mesh_edge(points,first,last,stride,tolerance))then
            call painter%drawLine(points(1,first),points(2,first),points(1,last),points(2,last))
        else
            call painter%drawPolyline(points(1,first:last:stride),points(2,first:last:stride))
        end if
    end subroutine paint_mesh_stroke

    pure logical function straight_mesh_edge(points,first,last,stride,tolerance) result(straight)
        real(dp), intent(in) :: points(:,:)
        integer, intent(in) :: first,last,stride
        real(dp), intent(in) :: tolerance
        real(dp) :: chord(2),delta(2),length_squared,along,cross,previous,cross_limit
        integer :: node

        straight=.true.
        if(last-first==stride)return
        straight=.false.
        chord=points(:2,last)-points(:2,first)
        length_squared=sum(chord*chord)
        if(.not.ieee_is_finite(length_squared).or.length_squared<.25_dp)return
        previous=0._dp
        cross_limit=tolerance**2*length_squared
        do node=first+stride,last-stride,stride
            delta=points(:2,node)-points(:2,first)
            along=dot_product(chord,delta)
            if(.not.ieee_is_finite(along).or.along<previous.or.along>length_squared)return
            cross=chord(1)*delta(2)-chord(2)*delta(1)
            if(.not.ieee_is_finite(cross).or.cross*cross>cross_limit)return
            previous=along
        end do
        straight=.true.
    end function straight_mesh_edge

    logical function mesh_boundary_face(points,counts,nodes,side,lower,upper,face) result(visible)
        real(dp), intent(in) :: points(:,:)
        integer, intent(in) :: counts(3),nodes(:),side,lower(3),upper(3)
        integer, intent(out) :: face(4)
        integer :: axis,stride(3),corner(3),target,inward,i

        visible=.false.
        face=nodes(mesh_hex_faces(:,side))
        axis=(side+1)/2
        stride=[1,counts(1),product(counts(:2))]
        corner=[(mod((nodes(1)-1)/stride(i),counts(i)),i=1,3)]
        if(any(corner<lower).or.any(corner>=upper))return
        target=lower(axis); inward=nodes(8)
        if(mod(side,2)==0)then
            target=upper(axis)
            inward=nodes(1)
        end if
        if(mod((face(1)-1)/stride(axis),counts(axis))/=target)return
        visible=mesh_face_front(points,face,inward)
    end function mesh_boundary_face

    pure logical function mesh_face_front(points,face,inward) result(visible)
        real(dp),intent(in) :: points(:,:)
        integer,intent(in) :: face(4),inward
        real(dp) :: a(3),b(3),c(3),normal(3)

        a=points(:,face(2))-points(:,face(1))
        b=points(:,face(4))-points(:,face(1))
        c=points(:,inward)-points(:,face(1))
        normal=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
        visible=normal(3)*dot_product(normal,c)<0._dp
    end function mesh_face_front

    subroutine paint_mesh_grid(painter,points,counts,lower,upper,kind,wireframe)
        type(QPainter), intent(inout) :: painter
        real(dp), intent(in) :: points(:,:)
        integer, intent(in) :: counts(3),lower(3),upper(3),kind
        logical, intent(in) :: wireframe
        integer :: stride(3),side,normal,along,across,base,i,first,last,span(3),line,lines,far,offset,border
        real(dp) :: spacing,pixel_span
        logical :: front

        ! Each front-facing side keeps the native direction-1-fastest node order.
        stride=[1,counts(1),product(counts(:2))]
        span=upper-lower
        border=merge(1,0,wireframe.and.kind==3)
        do side=merge(5,1,kind==2),merge(5,6,kind==2)
            front=.true.
            if(kind==3)front=volume_face_front(points,counts,side,lower,upper)
            if(.not.front.and..not.wireframe)cycle
            spacing=merge(3._dp,5._dp,front)
            normal=(side+1)/2
            base=1+sum(lower*stride)
            if(mod(side,2)==0)base=base+span(normal)*stride(normal)
            do along=1,3
                if(along==normal)cycle
                across=6-normal-along
                lines=span(across)
                if(wireframe.and.lines>1)then
                    far=base+span(along)*stride(along)
                    offset=span(across)*stride(across)
                    pixel_span=max(projected_line_span(base,offset,stride(across)), &
                        projected_line_span(far,offset,stride(across)))
                    if(pixel_span<spacing*real(lines,dp))lines=max(1,ceiling(pixel_span/spacing))
                end if
                do line=border,lines-border
                    i=line*span(across)/max(1,lines)
                    first=base+i*stride(across)
                    last=first+span(along)*stride(along)
                    if(wireframe)then
                        call paint_mesh_stroke(painter,points,first,last,stride(along),straight_grid_tolerance)
                    else
                        call painter%drawPolyline(points(1,first:last:stride(along)),points(2,first:last:stride(along)))
                    end if
                end do
            end do
        end do
    contains
        real(dp) function projected_line_span(first,offset,step) result(pixels)
            integer,intent(in) :: first,offset,step
            real(dp) :: bounds(2)

            bounds=maxval(points(:2,first:first+offset:step),dim=2)-minval(points(:2,first:first+offset:step),dim=2)
            pixels=norm2(bounds)
        end function projected_line_span
    end subroutine paint_mesh_grid

    subroutine pick_mesh_cell(s,x,y,mesh_id,part,cell,face,view_depth,opaque_only)
        type(studio_t), intent(in) :: s
        real(dp), intent(in) :: x,y
        integer, intent(out) :: mesh_id,part,cell,face
        real(dp), intent(out), optional :: view_depth
        logical, intent(in), optional :: opaque_only
        type(viewport_camera) :: camera
        type(projected_mesh_part), allocatable :: parts(:)
        real(dp) :: click(2),distance,depth,best_distance,best_depth,radius,parameter
        integer :: id,p,c,kind,side,vertices(4),i,j,k,base,step(2),extent(2)
        logical :: hit,require_opaque

        mesh_id=0; part=0; cell=0; face=0
        if(present(view_depth))view_depth=-huge(1._dp)
        require_opaque=.false.
        if(present(opaque_only))require_opaque=opaque_only
        if(require_opaque.and.(.not.s%mesh_shaded.or.s%viewport_style==viewport_wireframe))return
        if(.not.workspace_view_hit(s,x,y))return
        if(.not.prepare_camera(s,camera))return
        click=[x,y]
        best_distance=huge(1._dp)
        best_depth=-huge(1._dp)
        call project_visible_meshes(s,camera,.false.,parts,clip=[x-mesh_line_pick_radius,y-mesh_line_pick_radius, &
            x+mesh_line_pick_radius,y+mesh_line_pick_radius])
        do k=1,size(parts)
            id=parts(k)%mesh
            p=parts(k)%part
            kind=s%document%meshes(id)%kind
            if(require_opaque.and.kind==1)cycle
            radius=merge(mesh_line_pick_radius,0._dp,kind==1)
            associate(points=>parts(k)%points,counts=>parts(k)%counts,lower=>parts(k)%lower,upper=>parts(k)%upper)
                if(any(click<minval(points(:2,:),dim=2)-radius))cycle
                if(any(click>maxval(points(:2,:),dim=2)+radius))cycle
                if(allocated(parts(k)%section_faces))then
                    do i=1,size(parts(k)%section_faces)
                        associate(clipped=>parts(k)%section_faces(i))
                            if(.not.clipped%front)cycle
                            if(.not.mesh_face_hit(clipped%points(:,:clipped%count),click,depth))cycle
                            distance=0._dp
                            c=clipped%cell
                            call offer_hit(clipped%side)
                        end associate
                    end do
                    cycle
                end if
                associate(cells=>s%document%meshes(id)%parts(p)%cells)
                    select case(kind)
                    case(1)
                        do c=1,size(cells,1)
                            call nearest_curve_parameter(points(:,cells(c,:)),click,parameter,hit,distance,depth)
                            if(hit.and.distance<=radius**2)call offer_hit(0)
                        end do
                    case(2)
                        do c=1,size(cells,1)
                            vertices=cells(c,mesh_quad_nodes)
                            if(mesh_face_hit(points(:,vertices),click,depth))then
                                distance=0._dp
                                call offer_hit(0)
                            end if
                        end do
                    case(3)
                        do side=1,6
                            if(.not.s%mesh_shaded)then
                                if(.not.volume_face_front(points,counts,side,lower,upper))cycle
                            end if
                            call mesh_boundary_cells(counts,lower,upper,side,base,step,extent)
                            do j=0,extent(2)-1
                                do i=0,extent(1)-1
                                    c=base+i*step(1)+j*step(2)
                                    if(require_opaque.and.s%mesh_view==mesh_lattice)then
                                        if(.not.mesh_cell_in_selection(s,id,p,c))cycle
                                    end if
                                    if(.not.mesh_boundary_face(points,counts,cells(c,:),side,lower,upper,vertices))cycle
                                    if(.not.mesh_face_hit(points(:,vertices),click,depth))cycle
                                    distance=0._dp
                                    call offer_hit(side)
                                end do
                            end do
                        end do
                    end select
                end associate
            end associate
        end do
        if(present(view_depth))view_depth=best_depth
    contains
        subroutine offer_hit(side)
            integer, intent(in) :: side

            if(distance>best_distance)return
            if(distance==best_distance)then
                if(depth<best_depth)return
                if(depth==best_depth.and.mesh_id==s%selected_mesh)return
                if(depth==best_depth.and.id/=s%selected_mesh)return
            end if
            best_distance=distance
            best_depth=depth
            mesh_id=id; part=p; cell=c; face=side
        end subroutine offer_hit
    end subroutine pick_mesh_cell

    logical function mesh_face_hit(vertices,click,depth) result(hit)
        real(dp), intent(in) :: vertices(:,:),click(2)
        real(dp), intent(out) :: depth
        integer :: trial,rows(3)
        real(dp) :: weight(3)

        hit=.false.
        depth=-huge(1._dp)
        if(any(click<minval(vertices(:2,:),dim=2)).or.any(click>maxval(vertices(:2,:),dim=2)))return
        do trial=2,size(vertices,2)-1
            rows=[1,trial,trial+1]
            if(.not.screen_triangle(vertices(:2,rows),click,weight))cycle
            depth=max(depth,dot_product(weight,vertices(3,rows)))
            hit=.true.
        end do
    end function mesh_face_hit

    subroutine add_mesh_faces(quads,used,s,part,part_index,selected,preview,spacing,canvas)
        type(screen_quad), intent(inout) :: quads(:)
        integer, intent(inout) :: used
        type(studio_t), intent(in) :: s
        type(projected_mesh_part), intent(in) :: part
        integer, intent(in) :: part_index
        logical, intent(in) :: selected,preview
        real(dp), intent(in) :: spacing
        integer,intent(in) :: canvas(2)
        integer :: side,cell,vertices(4),i,j,base,step(2),extent(2)
        integer :: jump(3)
        real(dp) :: clip(4)
        logical :: clip_faces

        if(allocated(part%section_faces))return
        clip=[-4._dp,-4._dp,real(canvas,dp)+4._dp]
        associate(mesh=>s%document%meshes(part%mesh),points=>part%points)
            clip_faces=any(minval(points(:2,:),dim=2)<clip(:2)).or.any(maxval(points(:2,:),dim=2)>clip(3:))
            if(preview.and.mesh%kind>=2)then
                jump=mesh_overview_stride(part,spacing)
                if(any(jump>1))then
                    call append_mesh_overview(quads,used,part,part_index,mesh%kind,selected,jump,clip,clip_faces)
                    return
                end if
            end if
            associate(cells=>mesh%parts(part%part)%cells)
                select case(mesh%kind)
                case(2)
                    do cell=1,size(cells,1)
                        vertices=cells(cell,mesh_quad_nodes)
                        if(clip_faces)then
                            if(.not.mesh_quad_on_screen(points,vertices,clip))cycle
                        end if
                        used=used+1
                        quads(used)=screen_quad(vertices,0,0,.25_dp*sum(points(3,vertices)),selected,part_index,cell)
                    end do
                case(3)
                    do side=1,6
                        call mesh_boundary_cells(part%counts,part%lower,part%upper,side,base,step,extent)
                        do j=0,extent(2)-1
                            do i=0,extent(1)-1
                                cell=base+i*step(1)+j*step(2)
                                if(.not.mesh_boundary_face(points,part%counts,cells(cell,:),side, &
                                    part%lower,part%upper,vertices))cycle
                                if(clip_faces)then
                                    if(.not.mesh_quad_on_screen(points,vertices,clip))cycle
                                end if
                                used=used+1
                                quads(used)=screen_quad(vertices,side,0,.25_dp*sum(points(3,vertices)), &
                                    selected,part_index,cell)
                            end do
                        end do
                    end do
                end select
            end associate
        end associate
    end subroutine add_mesh_faces

    pure logical function mesh_quad_on_screen(points,face,clip) result(visible)
        real(dp),intent(in) :: points(:,:),clip(4)
        integer,intent(in) :: face(4)

        visible=all(maxval(points(:2,face),dim=2)>=clip(:2)).and.all(minval(points(:2,face),dim=2)<=clip(3:))
    end function mesh_quad_on_screen

    function mesh_overview_stride(part,spacing) result(jump)
        type(projected_mesh_part),intent(in) :: part
        real(dp),intent(in) :: spacing
        integer :: jump(3),span(3),stride(3),from,to,base,budget,axis
        real(dp) :: screen_span,tolerance

        jump=1
        span=part%upper-part%lower
        screen_span=maxval(maxval(part%points(:2,:),dim=2)-minval(part%points(:2,:),dim=2))
        if(.not.ieee_is_finite(screen_span).or.screen_span>=spacing*maxval(span))return
        budget=max(1,ceiling(screen_span/spacing))
        jump=max(1,1+(span-1)/budget)
        stride=[1,part%counts(1),product(part%counts(:2))]
        base=1+sum((part%lower+span/2)*stride)
        do axis=1,3
            if(jump(axis)==1)cycle
            from=base-span(axis)/2*stride(axis)
            to=from+span(axis)*stride(axis)
            tolerance=32._dp*epsilon(1._dp)*max(1._dp,maxval(abs(part%points(:,from))),maxval(abs(part%points(:,to))))
            if(norm2(part%points(:,from)-part%points(:,to))<=tolerance)jump(axis)=min(jump(axis),max(1,span(axis)/4))
        end do
    end function mesh_overview_stride

    subroutine append_mesh_overview(quads,used,part,part_index,kind,selected,jump,clip,clip_faces)
        type(screen_quad),intent(inout) :: quads(:)
        integer,intent(inout) :: used
        type(projected_mesh_part),intent(in) :: part
        integer,intent(in) :: part_index,kind,jump(3)
        logical,intent(in) :: selected,clip_faces
        real(dp),intent(in) :: clip(4)
        integer :: side,first,axis,along,across,span(3),stride(3),base,inward,i,j,start,a,b,face(4)

        stride=[1,part%counts(1),product(part%counts(:2))]
        span=part%upper-part%lower
        first=merge(1,5,kind==3)
        do side=first,merge(6,5,kind==3)
            axis=(side+1)/2
            along=mod(axis,3)+1
            across=mod(axis+1,3)+1
            base=1+sum(part%lower*stride)
            inward=stride(axis)
            if(mod(side,2)==0)then
                base=base+span(axis)*stride(axis)
                inward=-inward
            end if
            do j=0,span(across)-1,jump(across)
                b=min(jump(across),span(across)-j)*stride(across)
                do i=0,span(along)-1,jump(along)
                    a=min(jump(along),span(along)-i)*stride(along)
                    start=base+i*stride(along)+j*stride(across)
                    face=start+[0,a,a+b,b]
                    if(clip_faces)then
                        if(.not.mesh_quad_on_screen(part%points,face,clip))cycle
                    end if
                    if(kind==3)then
                        if(.not.mesh_face_front(part%points,face,start+inward))cycle
                    end if
                    used=used+1
                    quads(used)=screen_quad(face,side,0,.25_dp*sum(part%points(3,face)),selected,part_index,0)
                end do
            end do
        end do
    end subroutine append_mesh_overview

    function mesh_face_color(points,vertices,scale,selected,picked,opacity) result(color)
        real(dp), intent(in) :: points(:,:),scale
        integer, intent(in) :: vertices(4),opacity
        logical, intent(in) :: selected,picked
        type(QColor) :: color

        color=lit_mesh_color(projected_face_light(points,vertices,scale),selected,picked,opacity,.false.)
    end function mesh_face_color

    function mesh_section_color(face,scale,selected,picked,opacity) result(color)
        type(mesh_section_face), intent(in) :: face
        real(dp), intent(in) :: scale
        integer, intent(in) :: opacity
        logical, intent(in) :: selected,picked
        type(QColor) :: color

        color=lit_mesh_color(projected_face_light(face%points,[1,2,3,face%count],scale),selected,picked,opacity,face%side==0)
    end function mesh_section_color

    function lit_mesh_color(light,selected,picked,opacity,cap) result(color)
        real(dp), intent(in) :: light
        integer, intent(in) :: opacity
        logical, intent(in) :: selected,picked,cap
        type(QColor) :: color
        integer :: base(3),rgb(3)

        base=[56,113,128]
        if(selected)base=[62,141,157]
        if(cap)base=[115,173,175]
        if(picked)base=[186,105,208]
        rgb=min(255,nint(real(base,dp)*(.57_dp+.56_dp*light))+nint(15._dp*light))
        color=QColor(rgb(1),rgb(2),rgb(3),opacity)
    end function lit_mesh_color
end module studio_viewport_mesh
