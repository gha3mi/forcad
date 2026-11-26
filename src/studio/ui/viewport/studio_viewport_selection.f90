module studio_viewport_selection
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use qt, only: QColor, QPainter, QPainterPath, QPen, QPoint, Qt_DashLine, Qt_SolidLine
    use studio_selection, only: selection_target_t, selection_patch, selection_control, selection_mesh, selection_assembly, &
        pick_multipatch, pick_face, pick_edge, pick_vertex
    use studio_multipatch_collection, only: count_multipatches
    use studio_state, only: studio_t, dp, studio_mesh_count, area_box, area_lasso
    use studio_viewport_camera, only: viewport_camera, prepare_camera, workspace_view_bounds
    use studio_viewport_mesh, only: projected_mesh_part, project_visible_meshes
    use studio_viewport_scene, only: projected_patch, project_scene, projected_face_axes, projected_edge, edge_total, corner_sample, &
        projected_control_patch, project_selected_controls, project_shared_controls, &
        pick_control_points, viewport_patch_target, multipatch_control_scope
    implicit none(type, external)
    private
    public :: collect_viewport_region, paint_viewport_selection

    type :: viewport_region
        type(QPainterPath) :: path
        real(dp) :: bounds(4)=0._dp
        real(dp) :: excluded(4)=0._dp
        logical :: crossing=.false.,box=.false.,has_excluded=.false.
    end type viewport_region

contains

    subroutine collect_viewport_region(s,points,tool,targets)
        type(studio_t), intent(in) :: s
        real(dp), intent(in) :: points(:,:)
        integer, intent(in) :: tool
        type(selection_target_t), allocatable, intent(out) :: targets(:)
        type(viewport_camera) :: camera
        type(viewport_region) :: region
        type(projected_patch), allocatable :: objects(:)
        type(projected_control_patch),allocatable :: controls(:)
        type(selection_target_t), allocatable :: selected(:)
        type(selection_target_t) :: target
        logical,allocatable :: seen(:),included(:)
        logical :: hit
        integer :: i,j,n,capacity

        allocate(targets(0))
        if(.not.valid_outline(points,tool))return
        if(.not.prepare_camera(s,camera))return
        call visible_region(s,points,tool,region)
        if(.not.region%path%isEmpty())then
            n=0
            if(s%mesh_pick.or.s%mesh_only)then
                call collect_meshes(s,camera,region,targets)
            else if(s%feature_pick_mode==pick_control_points)then
                if(multipatch_control_scope(s))then
                    if(s%layers(4))call collect_shared_controls(s,camera,region,targets)
                else if(s%layers(4))then
                    call project_selected_controls(s,camera,region%bounds,controls)
                    capacity=0
                    do i=1,size(controls)
                        capacity=capacity+size(controls(i)%points,2)
                    end do
                    allocate(selected(capacity))
                    do i=1,size(controls)
                        associate(patch=>controls(i))
                            do j=1,size(patch%points,2)
                                if(.not.region%path%containsPoint(patch%points(1,j),patch%points(2,j)))cycle
                                n=n+1
                                selected(n)=selection_target_t(selection_control,patch%id,patch%kind,item=j)
                            end do
                        end associate
                    end do
                    targets=selected(:n)
                end if
            else if(s%feature_pick_mode>=pick_face.and.s%feature_pick_mode<=pick_vertex)then
                call collect_features(s,camera,region,targets)
            else
                if(s%feature_pick_mode==pick_multipatch.and..not.region%crossing)then
                    call project_scene(s,camera,objects,highlights=.false.)
                else
                    call project_scene(s,camera,objects,clip=region%bounds,highlights=.false.)
                end if
                if(s%feature_pick_mode==pick_multipatch)then
                    allocate(seen(count_multipatches(s%document%multipatches)),source=.false.)
                    allocate(included(size(seen)),source=.not.region%crossing)
                end if
                allocate(selected(size(objects)))
                do i=1,size(objects)
                    associate(obj=>objects(i))
                        target=viewport_patch_target(s,obj%kind,obj%variant)
                        hit=patch_hit(region,obj)
                        if(target%category==selection_assembly)then
                            j=target%id
                            seen(j)=.true.
                            if(hit.eqv.region%crossing)included(j)=hit
                        else if(hit)then
                            n=n+1
                            selected(n)=target
                        end if
                    end associate
                end do
                if(allocated(seen))then
                    do i=1,size(seen)
                        if(.not.(seen(i).and.included(i)))cycle
                        n=n+1
                        selected(n)=selection_target_t(selection_assembly,i)
                    end do
                end if
                targets=selected(:n)
            end if
        end if
        call region%path%delete()
    end subroutine collect_viewport_region

    subroutine collect_features(s,camera,region,targets)
        type(studio_t),intent(in) :: s
        type(viewport_camera),intent(in) :: camera
        type(viewport_region),intent(in) :: region
        type(selection_target_t),allocatable,intent(out) :: targets(:)
        type(projected_patch),allocatable :: objects(:)
        type(selection_target_t),allocatable :: selected(:)
        integer :: feature,i,item,n

        feature=s%feature_pick_mode
        call project_scene(s,camera,objects,clip=region%bounds,highlights=.false.)
        allocate(selected(sum([(feature_count(objects(i),feature),i=1,size(objects))])))
        n=0
        do i=1,size(objects)
            associate(patch=>objects(i))
                do item=1,feature_count(patch,feature)
                    if(.not.feature_hit(region,patch,feature,item))cycle
                    n=n+1
                    selected(n)=selection_target_t(selection_patch,patch%variant,patch%kind,item=item,feature=feature)
                end do
            end associate
        end do
        targets=selected(:n)
    end subroutine collect_features

    pure integer function feature_count(patch,feature) result(count)
        type(projected_patch),intent(in) :: patch
        integer,intent(in) :: feature
        integer,parameter :: faces(3)=[0,1,6]

        select case(feature)
        case(pick_face)
            count=faces(patch%kind)
        case(pick_edge)
            count=0
            if(.not.patch%is_point)count=edge_total(patch%kind)
        case(pick_vertex)
            count=merge(1,2**patch%kind,patch%is_point)
        case default
            count=0
        end select
    end function feature_count

    logical function feature_hit(region,patch,feature,item) result(hit)
        type(viewport_region),intent(in) :: region
        type(projected_patch),intent(in) :: patch
        integer,intent(in) :: feature,item
        real(dp),allocatable :: points(:,:)
        integer :: sides(2),n,node

        select case(feature)
        case(pick_face)
            hit=patch_face_hit(region,patch,item)
        case(pick_edge)
            points=projected_edge(patch,item,sides)
            n=size(points,2)
            hit=grid_hit(region,points,1,[n,1,1],[0,0,0],[n-1,0,0])
        case(pick_vertex)
            node=corner_sample(patch%counts,patch%kind,item)
            hit=region%path%containsPoint(patch%points(1,node),patch%points(2,node))
        case default
            hit=.false.
        end select
    end function feature_hit

    subroutine collect_shared_controls(s,camera,region,targets)
        type(studio_t),intent(in) :: s
        type(viewport_camera),intent(in) :: camera
        type(viewport_region),intent(in) :: region
        type(selection_target_t),allocatable,intent(inout) :: targets(:)
        type(projected_control_patch),allocatable :: patches(:)
        real(dp),allocatable :: projected(:,:)
        logical,allocatable :: included(:)
        integer :: i

        call project_selected_controls(s,camera,region%bounds,patches)
        call project_shared_controls(s,camera,patches,projected,included)
        do i=1,size(included)
            if(included(i))included(i)=region%path%containsPoint(projected(1,i),projected(2,i))
        end do
        targets=pack([(selection_target_t(selection_assembly,s%selected_multipatch,item=i),i=1,size(included))],included)
    end subroutine collect_shared_controls

    logical function patch_hit(region,patch) result(hit)
        type(viewport_region),intent(in) :: region
        type(projected_patch),intent(in) :: patch
        integer :: side
        logical :: side_hit

        if(patch%kind<3)then
            hit=grid_hit(region,patch%points,patch%kind,patch%counts,[0,0,0],patch%counts-1)
            return
        end if
        do side=1,6
            side_hit=patch_face_hit(region,patch,side)
            if(side_hit.eqv.region%crossing)exit
        end do
        hit=side_hit
    end function patch_hit

    logical function patch_face_hit(region,patch,side) result(hit)
        type(viewport_region),intent(in) :: region
        type(projected_patch),intent(in) :: patch
        integer,intent(in) :: side
        integer :: axes(2),counts(2)

        call projected_face_axes(patch,side,axes,counts)
        if(patch%kind==3)then
            hit=grid_hit(region,patch%boundary(side)%points,2,[counts,1],[0,0,0],[counts-1,0])
        else
            hit=grid_hit(region,patch%points,2,[counts,1],[0,0,0],[counts-1,0])
        end if
    end function patch_face_hit

    subroutine collect_meshes(s,camera,region,targets)
        type(studio_t), intent(in) :: s
        type(viewport_camera), intent(in) :: camera
        type(viewport_region), intent(in) :: region
        type(selection_target_t), allocatable, intent(out) :: targets(:)
        type(projected_mesh_part), allocatable :: parts(:)
        logical, allocatable :: included(:),visible(:)
        integer :: i,id,rank
        logical :: hit

        if(region%crossing)then
            call project_visible_meshes(s,camera,.false.,parts,clip=region%bounds)
        else
            call project_visible_meshes(s,camera,.false.,parts)
        end if
        allocate(visible(studio_mesh_count(s)),source=.false.)
        allocate(included(size(visible)),source=.not.region%crossing)
        do i=1,size(parts)
            associate(part=>parts(i))
                id=part%mesh
                visible(id)=.true.
                if(included(id).eqv.region%crossing)cycle
                rank=s%document%meshes(id)%kind
                hit=grid_hit(region,part%points,rank,part%counts,part%lower,part%upper)
                included(id)=hit
            end associate
        end do
        targets=pack([(selection_target_t(selection_mesh,id),id=1,size(visible))],included.and.visible)
    end subroutine collect_meshes

    logical function grid_hit(region,vertices,rank,counts,lower,upper) result(hit)
        type(viewport_region), intent(in) :: region
        real(dp), intent(in) :: vertices(:,:)
        integer, intent(in) :: rank,counts(3),lower(3),upper(3)
        integer :: i,j,k,node,stride(3),axis,step,a,b,plane,position,offset
        integer :: face(4)
        real(dp) :: bounds(4)
        logical :: inside

        hit=.false.
        if(any(counts<1).or.size(vertices,2)/=product(counts))return
        if(any(lower<0).or.any(upper>=counts).or.any(upper<lower))return
        bounds=[minval(vertices(:2,:),dim=2),maxval(vertices(:2,:),dim=2)]
        if(.not.boxes_overlap(region%bounds,bounds))return
        if(region%box.and.all(lower==0).and.all(upper==counts-1))then
            if(all(bounds(:2)>region%bounds(:2)).and.all(bounds(3:)<region%bounds(3:)))then
                if(.not.region%has_excluded.or..not.boxes_overlap(region%excluded,bounds))then
                    hit=.true.
                    return
                end if
            end if
        end if
        stride=[1,counts(1),product(counts(:2))]
        do k=lower(3),upper(3)
            do j=lower(2),upper(2)
                do i=lower(1),upper(1)
                    node=1+i+j*stride(2)+k*stride(3)
                    inside=region%path%containsPoint(vertices(1,node),vertices(2,node))
                    if(region%crossing.eqv.inside)then
                        hit=inside
                        return
                    end if
                end do
            end do
        end do
        if(.not.region%crossing)then
            hit=.true.
            return
        end if
        if(rank==1)then
            do i=lower(1)+1,upper(1)
                if(segment_hit(region,vertices(:2,i),vertices(:2,i+1)))then
                    hit=.true.
                    return
                end if
            end do
            return
        end if
        do axis=3,1,-1
            if(rank==2.and.axis/=3)exit
            a=mod(axis,3)+1
            b=mod(axis+1,3)+1
            step=max(1,upper(axis)-lower(axis))
            do plane=lower(axis),upper(axis),step
                offset=1+plane*stride(axis)
                do j=lower(b),upper(b)-1
                    do i=lower(a),upper(a)-1
                        position=offset+i*stride(a)+j*stride(b)
                        face=position+[0,stride(a),stride(a)+stride(b),stride(b)]
                        if(footprint_hit(region,vertices(:2,face)))then
                            hit=.true.
                            return
                        end if
                    end do
                end do
            end do
        end do
    end function grid_hit

    logical function segment_hit(region,a,b) result(hit)
        type(viewport_region), intent(in) :: region
        real(dp), intent(in) :: a(2),b(2)
        real(dp) :: normal(2),length,footprint(2,4)

        hit=.false.
        length=norm2(b-a)
        if(length<=tiny(1._dp))return
        normal=.5_dp*[a(2)-b(2),b(1)-a(1)]/length
        footprint(:,1)=a+normal
        footprint(:,2)=b+normal
        footprint(:,3)=b-normal
        footprint(:,4)=a-normal
        hit=footprint_hit(region,footprint)
    end function segment_hit

    logical function footprint_hit(region,polygon) result(hit)
        type(viewport_region), intent(in) :: region
        real(dp), intent(in) :: polygon(:,:)
        type(QPainterPath) :: shape,intersection

        hit=.false.
        if(.not.overlaps(region%bounds,polygon))return
        shape=QPainterPath()
        call shape%addPolygon(polygon(1,:),polygon(2,:))
        call shape%closeSubpath()
        intersection=shape%intersected(region%path)
        hit=.not.intersection%isEmpty()
        call intersection%delete()
        call shape%delete()
    end function footprint_hit

    logical function overlaps(bounds,points)
        real(dp), intent(in) :: bounds(4),points(:,:)

        overlaps=all(maxval(points,dim=2)>=bounds(:2)).and.all(minval(points,dim=2)<=bounds(3:4))
    end function overlaps

    pure logical function boxes_overlap(first,second) result(overlap)
        real(dp),intent(in) :: first(4),second(4)

        overlap=all(first(3:)>=second(:2)).and.all(second(3:)>=first(:2))
    end function boxes_overlap

    subroutine visible_region(s,points,tool,region)
        type(studio_t), intent(in) :: s
        real(dp), intent(in) :: points(:,:)
        integer, intent(in) :: tool
        type(viewport_region), intent(out) :: region
        type(QPainterPath) :: shape,aperture,controls,uncovered
        type(QPoint) :: origin
        real(dp) :: bounds(4)

        bounds=workspace_view_bounds(s)
        shape=region_outline(points,tool)
        aperture=QPainterPath()
        call aperture%addRect(bounds(1),bounds(2),bounds(3)-bounds(1),bounds(4)-bounds(2))
        region%path=shape%intersected(aperture)
        call aperture%delete()
        call shape%delete()
        if(.not.s%viewport_controls%isNull())then
            if(s%viewport_controls%isVisible())then
                controls=QPainterPath()
                origin=s%viewport_controls%mapTo(s%surface,QPoint(0,0))
                region%has_excluded=.true.
                region%excluded(:2)=[origin%x,origin%y]
                region%excluded(3:)=region%excluded(:2)+[s%viewport_controls%width(),s%viewport_controls%height()]
                call controls%addRect(origin%x,origin%y,s%viewport_controls%width(),s%viewport_controls%height())
                uncovered=region%path%subtracted(controls)
                call region%path%delete()
                region%path=uncovered
                call controls%delete()
            end if
        end if
        region%bounds(:2)=[region%path%boundingX(),region%path%boundingY()]
        region%bounds(3:4)=region%bounds(:2)+[region%path%boundingWidth(),region%path%boundingHeight()]
        region%crossing=tool==area_lasso.or.points(1,1)>points(1,size(points,2))
        region%box=tool==area_box
    end subroutine visible_region

    logical function valid_outline(points,tool) result(valid)
        real(dp), intent(in) :: points(:,:)
        integer, intent(in) :: tool

        valid=.false.
        if(size(points,1)/=2.or.size(points,2)<2)return
        if(.not.all(ieee_is_finite(points)))return
        select case(tool)
        case(area_box); valid=size(points,2)==2.and.all(abs(points(:,2)-points(:,1))>=1._dp)
        case(area_lasso); valid=size(points,2)>=3
        case default; valid=.false.
        end select
    end function valid_outline

    function region_outline(points,tool) result(path)
        real(dp), intent(in) :: points(:,:)
        integer, intent(in) :: tool
        type(QPainterPath) :: path
        real(dp) :: origin(2),extent(2)

        path=QPainterPath()
        if(tool==area_box)then
            origin=min(points(:,1),points(:,2))
            extent=abs(points(:,2)-points(:,1))
            call path%addRect(origin(1),origin(2),extent(1),extent(2))
        else
            call path%addPolygon(points(1,:),points(2,:))
            call path%closeSubpath()
        end if
    end function region_outline

    subroutine paint_viewport_selection(painter,s)
        type(QPainter), intent(inout) :: painter
        type(studio_t), intent(in) :: s
        type(QPainterPath) :: outline
        type(QColor) :: edge,fill
        real(dp) :: bounds(4)
        integer :: n,style
        logical :: crossing

        n=s%canvas_gesture%area_count
        if(n<2)return
        associate(points=>s%canvas_gesture%area_points(:,:n))
            outline=region_outline(points,s%area_selection_tool)
            crossing=s%area_selection_tool==area_lasso.or.points(1,n)<points(1,1)
        end associate
        edge=QColor(113,177,249)
        fill=QColor(75,148,223,45)
        style=Qt_SolidLine
        if(crossing)then
            edge=QColor(102,228,188)
            fill=QColor(69,185,158,44)
            style=Qt_DashLine
        end if
        bounds=workspace_view_bounds(s)
        call painter%save()
        call painter%setClipRect(bounds(1),bounds(2),bounds(3)-bounds(1),bounds(4)-bounds(2))
        call painter%setPen(QPen(edge,1.4_dp,style))
        call painter%setBrush(fill)
        call painter%drawPath(outline)
        call painter%restore()
        call outline%delete()
    end subroutine paint_viewport_selection

end module studio_viewport_selection
