module studio_viewport_primitives
    use qt, only: QPainter, QPainterPath
    use studio_state, only: studio_t, dp, viewport_transform_t
    use studio_viewport_camera, only: viewport_camera, project_vertex
    implicit none(type, external)
    private

    public :: screen_quad, project_points, preview_transform, draw_dot, paint_tensor_lines, append_tensor_path
    public :: volume_face_front, projected_face_light, nearest_curve_parameter, screen_triangle, screen_triangle_nearest
    public :: pick_sample_grid

    real(dp), parameter :: pick_radius=17._dp

    type :: screen_quad
        integer :: vertices(4), face, object
        real(dp) :: depth
        logical :: selected
        integer :: mesh_part=0, mesh_cell=0
    end type screen_quad

contains

    function project_points(s,camera,xyz,preview) result(projected)
        type(studio_t),intent(in) :: s
        type(viewport_camera), intent(in) :: camera
        real(dp), intent(in) :: xyz(:,:)
        logical,intent(in),optional :: preview
        real(dp), allocatable :: projected(:,:)
        integer :: point
        logical :: moving
        real(dp) :: position(3)

        moving=.false.
        if(present(preview))moving=preview
        allocate(projected(3,size(xyz,2)))
        do point=1,size(xyz,2)
            position=xyz(:,point)
            if(moving)position=preview_transform(s%transform_drag,position)
            projected(:,point)=project_vertex(camera,position)
        end do
    end function project_points

    pure function preview_transform(drag,position) result(moved)
        type(viewport_transform_t),intent(in) :: drag
        real(dp),intent(in) :: position(3)
        real(dp) :: moved(3)

        ! Transform only display samples. The edit is published through ForCAD on release.
        moved=drag%pivot+matmul(drag%linear,position-drag%pivot)+drag%offset
    end function preview_transform

    subroutine draw_dot(painter,center,radius)
        type(QPainter), intent(inout) :: painter
        real(dp), intent(in) :: center(2),radius
        call painter%drawEllipse(center(1)-radius,center(2)-radius,2._dp*radius,2._dp*radius)
    end subroutine draw_dot

    subroutine paint_tensor_lines(painter,projected,counts,axes,line_limit,boundary_only,visible_faces)
        type(QPainter), intent(inout) :: painter
        real(dp), intent(in) :: projected(:,:)
        integer, intent(in) :: counts(3),line_limit
        logical, intent(in) :: axes(3),boundary_only
        logical, intent(in), optional :: visible_faces(6)

        call trace_tensor_lines(projected,counts,axes,line_limit,boundary_only,visible_faces,painter=painter)
    end subroutine paint_tensor_lines

    subroutine append_tensor_path(path,projected,counts,axes)
        type(QPainterPath), intent(inout) :: path
        real(dp), intent(in) :: projected(:,:)
        integer, intent(in) :: counts(3)
        logical, intent(in) :: axes(3)

        call trace_tensor_lines(projected,counts,axes,maxval(counts),.false.,path=path)
    end subroutine append_tensor_path

    subroutine trace_tensor_lines(projected,counts,axes,line_limit,boundary_only,visible_faces,painter,path)
        type(QPainter), intent(inout), optional :: painter
        real(dp), intent(in) :: projected(:,:)
        integer, intent(in) :: counts(3),line_limit
        logical, intent(in) :: axes(3),boundary_only
        logical, intent(in), optional :: visible_faces(6)
        type(QPainterPath), intent(inout), optional :: path
        integer :: stride(3),axis,first,second,i,j,lines_along_first,lines_along_second
        integer :: index_first,index_second,base,last,step,point

        stride=[1,counts(1),counts(1)*counts(2)]
        do axis=1,3
            if (.not.axes(axis).or.counts(axis)<2) cycle
            first=mod(axis,3)+1
            second=mod(axis+1,3)+1
            lines_along_first=min(line_limit,counts(first)-1)
            lines_along_second=min(line_limit,counts(second)-1)
            step=stride(axis)
            do j=0,lines_along_second
                index_second=j*(counts(second)-1)/max(1,lines_along_second)
                do i=0,lines_along_first
                    if (boundary_only.and.i>0.and.i<lines_along_first.and.j>0.and.j<lines_along_second) cycle
                    if(present(visible_faces))then
                        if(.not.((i==0.and.visible_faces(2*first-1)).or. &
                            (i==lines_along_first.and.visible_faces(2*first)).or. &
                            (j==0.and.visible_faces(2*second-1)).or. &
                            (j==lines_along_second.and.visible_faces(2*second))))cycle
                    end if
                    index_first=i*(counts(first)-1)/max(1,lines_along_first)
                    base=1+index_first*stride(first)+index_second*stride(second)
                    last=base+(counts(axis)-1)*step
                    if(present(path))then
                        call path%moveTo(projected(1,base),projected(2,base))
                        do point=base+step,last,step
                            call path%lineTo(projected(1,point),projected(2,point))
                        end do
                    else
                        call painter%drawPolyline(projected(1,base:last:step),projected(2,base:last:step))
                    end if
                end do
            end do
        end do
    end subroutine trace_tensor_lines

    pure logical function volume_face_front(points,counts,side,lower,upper) result(front)
        real(dp), intent(in) :: points(:,:)
        integer, intent(in) :: counts(3),side
        integer, intent(in), optional :: lower(3),upper(3)
        integer :: axis,first,second,base,inward,stride(3),origin(3),span(3)
        real(dp) :: a(3),b(3),c(3),normal(3)

        origin=0
        if(present(lower))origin=lower
        span=counts-1-origin
        if(present(upper))span=upper-origin
        stride=[1,counts(1),product(counts(:2))]
        axis=(side+1)/2; first=mod(axis,3)+1; second=mod(axis+1,3)+1
        base=1+sum(origin*stride)+(span(first)-1)/2*stride(first)+(span(second)-1)/2*stride(second)
        inward=stride(axis)
        if(mod(side,2)==0)then
            base=base+span(axis)*stride(axis)
            inward=-inward
        end if
        a=points(:,base+stride(first))-points(:,base)
        b=points(:,base+stride(second))-points(:,base)
        c=points(:,base+inward)-points(:,base)
        normal=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
        front=normal(3)*dot_product(normal,c)<0._dp
    end function volume_face_front

    pure real(dp) function projected_face_light(points,vertices,scale) result(light)
        real(dp), intent(in) :: points(:,:),scale
        integer, intent(in) :: vertices(4)
        real(dp) :: a(3),b(3),normal(3),length

        a=points(:,vertices(2))-points(:,vertices(1))
        b=points(:,vertices(4))-points(:,vertices(1))
        a(:2)=a(:2)/scale; b(:2)=b(:2)/scale
        normal=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
        length=norm2(normal)
        light=0._dp
        if(length>tiny(1._dp))light=abs(dot_product(normal/length,[-.35_dp,-.45_dp,.82_dp]))
    end function projected_face_light

    subroutine nearest_curve_parameter(projected,click,u,hit,distance_squared,view_depth)
        real(dp), intent(in) :: projected(:,:),click(2)
        real(dp), intent(out) :: u
        logical, intent(out) :: hit
        real(dp), intent(out), optional :: distance_squared,view_depth
        real(dp) :: delta(2),offset(2),fraction,distance,nearest,length_squared,depth,nearest_depth
        integer :: segment

        nearest=pick_radius**2
        nearest_depth=-huge(1._dp)
        hit=.false.
        u=0._dp
        if(size(projected,2)==1)then
            nearest=sum((click-projected(:2,1))**2)
            nearest_depth=projected(3,1)
            hit=nearest<=pick_radius**2
        end if
        do segment=1,size(projected,2)-1
            delta=projected(:2,segment+1)-projected(:2,segment)
            length_squared=sum(delta**2)
            offset=click-projected(:2,segment)
            fraction=0._dp
            if(length_squared>=1.e-8_dp)fraction=max(0._dp,min(1._dp,dot_product(offset,delta)/length_squared))
            distance=sum((offset-fraction*delta)**2)
            if (distance>nearest) cycle
            depth=(1._dp-fraction)*projected(3,segment)+fraction*projected(3,segment+1)
            if (distance==nearest.and.depth<=nearest_depth) cycle
            nearest=distance; nearest_depth=depth
            u=(real(segment-1,dp)+fraction)/real(size(projected,2)-1,dp)
            hit=.true.
        end do
        if(present(distance_squared))distance_squared=nearest
        if(present(view_depth))view_depth=nearest_depth
    end subroutine nearest_curve_parameter

    logical function screen_triangle(vertices,click,weight) result(hit)
        real(dp), intent(in) :: vertices(2,3),click(2)
        real(dp), intent(out) :: weight(3)
        real(dp) :: edge1(2),edge2(2),relative(2),determinant
        hit=.false.
        edge1=vertices(:,2)-vertices(:,1)
        edge2=vertices(:,3)-vertices(:,1)
        relative=click-vertices(:,1)
        determinant=edge1(1)*edge2(2)-edge1(2)*edge2(1)
        if (abs(determinant)<1.e-6_dp) return
        weight(2)=(relative(1)*edge2(2)-relative(2)*edge2(1))/determinant
        weight(3)=(edge1(1)*relative(2)-edge1(2)*relative(1))/determinant
        weight(1)=1._dp-weight(2)-weight(3)
        hit=all(weight>=-1.e-8_dp)
    end function screen_triangle

    subroutine screen_triangle_nearest(vertices,click,limit,weight,distance)
        real(dp),intent(in) :: vertices(2,3),click(2),limit
        real(dp),intent(out) :: weight(3),distance
        real(dp) :: edge(2),offset(2),length,fraction,d,radius
        integer :: i,j

        distance=huge(1._dp)
        radius=sqrt(limit)
        if(any(click<minval(vertices,dim=2)-radius).or.any(click>maxval(vertices,dim=2)+radius))return
        if(screen_triangle(vertices,click,weight))then
            distance=0._dp
            return
        end if
        do i=1,3
            j=mod(i,3)+1
            offset=click-vertices(:,i)
            edge=vertices(:,j)-vertices(:,i)
            length=sum(edge**2)
            fraction=0._dp
            if(length>=1.e-8_dp)fraction=max(0._dp,min(1._dp,dot_product(offset,edge)/length))
            d=sum((offset-fraction*edge)**2)
            if(d>limit.or.d>=distance)cycle
            distance=d
            weight=0._dp
            weight(i)=1._dp-fraction
            weight(j)=fraction
        end do
    end subroutine screen_triangle_nearest

    subroutine pick_sample_grid(projected,nu,nv,click,uv,hit)
        real(dp), intent(in) :: projected(:,:),click(2)
        integer, intent(in) :: nu,nv
        real(dp), intent(out) :: uv(2)
        logical, intent(out) :: hit
        integer, parameter :: triangles(3,2)=reshape([1,2,3,1,3,4],[3,2])
        integer :: i,j,start,corners(4),nodes(3),trial,nearest_node
        real(dp) :: depth,best_depth,weight(3),local_u(4),local_v(4),distance,nearest

        hit=.false.
        uv=0._dp
        best_depth=-huge(1._dp)
        do j=1,nv-1
            local_v=real([j-1,j-1,j,j],dp)/real(nv-1,dp)
            do i=1,nu-1
                local_u=real([i-1,i,i,i-1],dp)/real(nu-1,dp)
                start=(j-1)*nu+i
                corners=start+[0,1,nu+1,nu]
                do trial=1,2
                    nodes=triangles(:,trial)
                    if (.not.screen_triangle(projected(:2,corners(nodes)),click,weight)) cycle
                    depth=dot_product(weight,projected(3,corners(nodes)))
                    if (depth<best_depth) cycle
                    best_depth=depth
                    uv=[dot_product(weight,local_u(nodes)),dot_product(weight,local_v(nodes))]
                    hit=.true.
                end do
            end do
        end do
        if (hit) return

        ! Near a boundary or edge-on view, allow snapping to the closest sample.
        nearest=pick_radius**2
        nearest_node=0
        do i=1,size(projected,2)
            distance=sum((projected(:2,i)-click)**2)
            if (distance>=nearest) cycle
            nearest=distance
            nearest_node=i
        end do
        if (nearest_node==0) return
        uv=[real(mod(nearest_node-1,nu),dp)/real(nu-1,dp),real((nearest_node-1)/nu,dp)/real(nv-1,dp)]
        hit=.true.
    end subroutine pick_sample_grid

end module studio_viewport_primitives
