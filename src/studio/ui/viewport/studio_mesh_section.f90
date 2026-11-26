module studio_mesh_section
    use studio_state, only: dp
    implicit none(type, external)
    private
    public :: mesh_section_face, build_mesh_section, mesh_boundary_cells, mesh_hex_faces

    integer, parameter :: max_polygon_nodes=12
    integer, parameter :: mesh_hex_faces(4,6)=reshape([1,3,7,5, 2,4,8,6, 1,2,6,5, 3,4,8,7, 1,2,4,3, 5,6,8,7],[4,6])
    integer, parameter :: hex_edges(2,12)=reshape([1,2, 3,4, 5,6, 7,8, 1,3, 2,4, 5,7, 6,8, 1,5, 2,6, 3,7, 4,8],[2,12])

    type :: mesh_section_face
        real(dp) :: points(3,max_polygon_nodes)=0._dp
        integer :: count=0, cell=0, side=0
        logical :: front=.false.
    end type mesh_section_face

contains

    subroutine build_mesh_section(world,screen,cells,counts,lower,upper,axis,coordinate,cap_front,faces)
        real(dp), intent(in) :: world(:,:),screen(:,:),coordinate
        integer, intent(in) :: cells(:,:),counts(3),lower(3),upper(3),axis
        logical, intent(in) :: cap_front
        type(mesh_section_face), allocatable, intent(out) :: faces(:)
        type(mesh_section_face), allocatable :: buffer(:)
        type(mesh_section_face) :: face
        real(dp), allocatable :: heights(:)
        real(dp) :: tolerance,spacing
        integer :: strides(3),side,cell,nodes(8),corners(4),inside,i,j,k,base,step(2),extent(2),used

        heights=world(axis,:)-coordinate
        spacing=max(maxval(abs(heights)),tiny(1._dp))
        tolerance=32._dp*epsilon(1._dp)*spacing
        if(all(heights<=tolerance))return
        allocate(buffer(0))
        used=0
        if(all(heights>tolerance))then
            call move_alloc(buffer,faces)
            return
        end if
        where(abs(heights)<=tolerance)heights=0._dp

        do side=1,6
            call mesh_boundary_cells(counts,lower,upper,side,base,step,extent)
            do j=0,extent(2)-1
                do i=0,extent(1)-1
                    cell=base+i*step(1)+j*step(2)
                    nodes=cells(cell,:)
                    corners=nodes(mesh_hex_faces(:,side))
                    inside=nodes(8)
                    if(mod(side,2)==0)inside=nodes(1)
                    call clip_mesh_face(screen(:,corners),heights(corners),face)
                    if(face%count<3)cycle
                    face%cell=cell
                    face%side=side
                    face%front=front_facing(screen(:,corners),screen(:,inside))
                    call append_face(buffer,used,face)
                end do
            end do
        end do

        strides=[1,counts(1)-1,product(counts(:2)-1)]
        do k=lower(3),upper(3)-1
            do j=lower(2),upper(2)-1
                do i=lower(1),upper(1)-1
                    cell=1+dot_product([i,j,k],strides)
                    nodes=cells(cell,:)
                    if(minval(heights(nodes))>=0._dp.or.maxval(heights(nodes))<0._dp)cycle
                    if(all(heights(nodes)<=0._dp))then
                        if(on_visible_boundary(nodes,counts,lower,upper,heights))cycle
                    end if
                    call cut_mesh_cell(world(:,nodes),screen(:,nodes),heights(nodes),axis,face)
                    if(face%count<3)cycle
                    face%cell=cell
                    face%front=cap_front
                    call append_face(buffer,used,face)
                end do
            end do
        end do
        if(used==size(buffer))then
            call move_alloc(buffer,faces)
        else
            faces=buffer(:used)
        end if
    end subroutine build_mesh_section

    pure subroutine clip_mesh_face(screen,heights,face)
        real(dp), intent(in) :: screen(3,4),heights(4)
        type(mesh_section_face), intent(out) :: face
        integer :: i,j

        face=mesh_section_face()
        do i=1,4
            j=mod(i,4)+1
            if(heights(i)<=0._dp)then
                face%count=face%count+1
                face%points(:,face%count)=screen(:,i)
            end if
            if(.not.crosses_plane(heights(i),heights(j)))cycle
            face%count=face%count+1
            face%points(:,face%count)=intersection(screen(:,i),screen(:,j),heights(i),heights(j))
        end do
    end subroutine clip_mesh_face

    pure subroutine cut_mesh_cell(world,screen,heights,axis,face)
        real(dp), intent(in) :: world(3,8),screen(3,8),heights(8)
        integer, intent(in) :: axis
        type(mesh_section_face), intent(out) :: face
        real(dp) :: profile(3,max_polygon_nodes),center(3),point(3),angle(max_polygon_nodes),bearing
        integer :: edge,node,first,last,i,j,axes(2),count

        face=mesh_section_face()
        count=0
        do node=1,8
            if(heights(node)/=0._dp)cycle
            count=count+1
            profile(:,count)=world(:,node)
            face%points(:,count)=screen(:,node)
        end do
        do edge=1,size(hex_edges,2)
            first=hex_edges(1,edge)
            last=hex_edges(2,edge)
            if(.not.crosses_plane(heights(first),heights(last)))cycle
            count=count+1
            profile(:,count)=intersection(world(:,first),world(:,last),heights(first),heights(last))
            face%points(:,count)=intersection(screen(:,first),screen(:,last),heights(first),heights(last))
        end do
        if(count<3)return

        face%count=count
        axes=[mod(axis,3)+1,mod(axis+1,3)+1]
        center=profile(:,1)+sum(profile(:,:count)-spread(profile(:,1),2,count),dim=2)/real(count,dp)
        angle(:count)=atan2(profile(axes(2),:count)-center(axes(2)),profile(axes(1),:count)-center(axes(1)))
        do i=2,count
            bearing=angle(i)
            point=face%points(:,i)
            j=i-1
            do while(j>=1)
                if(angle(j)<=bearing)exit
                angle(j+1)=angle(j)
                face%points(:,j+1)=face%points(:,j)
                j=j-1
            end do
            angle(j+1)=bearing
            face%points(:,j+1)=point
        end do
    end subroutine cut_mesh_cell

    pure logical function crosses_plane(first,last) result(crosses)
        real(dp), intent(in) :: first,last

        crosses=(first<0._dp.and.last>0._dp).or.(first>0._dp.and.last<0._dp)
    end function crosses_plane

    pure function intersection(first,last,first_height,last_height) result(point)
        real(dp), intent(in) :: first(3),last(3),first_height,last_height
        real(dp) :: point(3),a,b,scale

        scale=max(abs(first_height),abs(last_height))
        a=abs(first_height)/scale
        b=abs(last_height)/scale
        point=(b/(a+b))*first+(a/(a+b))*last
    end function intersection

    pure logical function on_visible_boundary(nodes,counts,lower,upper,heights) result(found)
        integer, intent(in) :: nodes(8),counts(3),lower(3),upper(3)
        real(dp), intent(in) :: heights(:)
        integer :: side,axis,plane,corner(3),stride(3),i

        found=.false.
        stride=[1,counts(1),product(counts(:2))]
        corner=[(mod((nodes(1)-1)/stride(i),counts(i)),i=1,3)]
        do side=1,6
            if(any(heights(nodes(mesh_hex_faces(:,side)))/=0._dp))cycle
            axis=(side+1)/2
            plane=merge(upper(axis)-1,lower(axis),mod(side,2)==0)
            if(corner(axis)/=plane)cycle
            found=.true.
            return
        end do
    end function on_visible_boundary

    pure logical function front_facing(face,inside) result(front)
        real(dp), intent(in) :: face(3,4),inside(3)
        real(dp) :: a(3),b(3),inward(3),normal(3)

        a=face(:,2)-face(:,1)
        b=face(:,4)-face(:,1)
        inward=inside-face(:,1)
        normal=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
        front=normal(3)*dot_product(normal,inward)<0._dp
    end function front_facing

    subroutine append_face(buffer,used,face)
        type(mesh_section_face), allocatable, intent(inout) :: buffer(:)
        integer, intent(inout) :: used
        type(mesh_section_face), intent(in) :: face
        type(mesh_section_face), allocatable :: larger(:)

        if(used==size(buffer))then
            allocate(larger(max(64,2*used)))
            larger(:used)=buffer
            call move_alloc(larger,buffer)
        end if
        used=used+1
        buffer(used)=face
    end subroutine append_face

    pure subroutine mesh_boundary_cells(counts,lower,upper,side,base,step,extent)
        integer, intent(in) :: counts(3),lower(3),upper(3),side
        integer, intent(out) :: base,step(2),extent(2)
        integer :: strides(3),span(3),normal,along,across

        strides=[1,counts(1)-1,product(counts(:2)-1)]
        span=upper-lower
        normal=(side+1)/2
        along=mod(normal,3)+1
        across=mod(normal+1,3)+1
        base=1+sum(lower*strides)
        if(mod(side,2)==0)base=base+(span(normal)-1)*strides(normal)
        step=[strides(along),strides(across)]
        extent=[span(along),span(across)]
    end subroutine mesh_boundary_cells
end module studio_mesh_section
