module studio_viewport_geometry
    use, intrinsic :: iso_fortran_env, only: int8
    !$ use omp_lib, only: omp_get_max_threads
    use qt, only: QColor,QPainter,QPen,Qt_NoPen
    use studio_state, only: studio_t,dp,mesh_lattice,viewport_shaded,viewport_wireframe,scene_patch_drawn,studio_mesh_count
    use studio_workspace_selection, only: mesh_in_selection
    use studio_selection, only: selection_mesh,patch_edge_boundaries
    use studio_viewport_camera, only: viewport_camera,viewport_multisampling,viewport_raster_step,viewport_preview_raster, &
        viewport_navigation_active
    use studio_viewport_primitives, only: screen_quad,draw_dot
    use studio_viewport_scene, only: projected_patch,corner_sample,edge_total,projected_edge,projected_face_axes, &
        pick_face,pick_edge,pick_vertex
    use studio_viewport_mesh, only: projected_mesh_part,mesh_covered_patches,add_mesh_faces,mesh_face_color,mesh_section_color
    use studio_viewport_depth, only: viewport_depth_t
    implicit none(type,external)
    private
    public :: paint_viewport_geometry

    integer,parameter :: first_triangle(3)=[1,2,3],second_triangle(3)=[1,3,4]
    real(dp),parameter :: edge_ink(3)=[47._dp,66._dp,80._dp],selection_ink(3)=[130._dp,235._dp,255._dp]
    real(dp),parameter :: hover_ink(3)=[196._dp,230._dp,244._dp]
    real(dp),parameter :: mesh_ink(3)=[205._dp,164._dp,103._dp]
    integer,parameter :: shared_mesh_opacity(*)=[180,233,249,253]
    integer,parameter :: dense_preview_patches=1024
    integer,parameter :: parallel_face_count=256,parallel_pixel_count=150000

    type :: face_colors
        real(dp),allocatable :: rgb(:,:)
    end type face_colors

    type :: mesh_cell_selection
        logical,allocatable :: selected(:)
    end type mesh_cell_selection

contains

    subroutine paint_viewport_geometry(painter,s,camera,objects,parts,width,height,mesh_coverage,navigating)
        type(QPainter),intent(inout) :: painter
        type(studio_t),intent(in) :: s
        type(viewport_camera),intent(in) :: camera
        type(projected_patch),intent(in) :: objects(:)
        type(projected_mesh_part),intent(in) :: parts(:)
        integer,intent(in) :: width,height
        logical,intent(in),optional,target :: mesh_coverage(:,:)
        logical,intent(in),optional :: navigating
        type(viewport_depth_t) :: buffer
        type(screen_quad),allocatable :: quads(:)
        type(face_colors),allocatable :: shading(:,:)
        type(mesh_cell_selection) :: selection(size(parts))
        integer,allocatable :: smooth_edges(:),transparent(:)
        integer,allocatable :: mesh_quads(:,:)
        integer(int8),allocatable :: edge_uses(:,:)
        integer :: capacity,used,transparent_count,object,part,side,primitive,counts(3),lower(2),upper(2),pixel_step,suppressed_edges
        integer :: transparent_index, threads, bands, band, band_height
        real(dp) :: bias,mesh_spacing
        logical :: cad_fill(size(objects)),wireframe,overview,light,motion,translucent_faces
        logical,allocatable :: selected_mesh(:)
        logical,allocatable,target :: coverage(:,:)
        logical,pointer :: replaced(:,:)

        if(size(objects)+size(parts)==0)return
        motion=.false.
        if(present(navigating))motion=navigating
        if(present(mesh_coverage))then
            replaced=>mesh_coverage
        else
            call mesh_covered_patches(s,coverage)
            replaced=>coverage
        end if
        do object=1,size(objects)
            associate(patch=>objects(object))
                cad_fill(object)=.true.
                if(patch%variant>=1.and.patch%variant<=size(replaced,1))then
                    cad_fill(object)=.not.replaced(patch%variant,patch%kind)
                end if
            end associate
        end do
        wireframe=s%viewport_style==viewport_wireframe
        if(wireframe)then
            call painter%save()
            do object=1,size(objects)
                associate(patch=>objects(object))
                    if(s%mesh_section_axis>0.and.patch%kind==3.and..not.cad_fill(object))cycle
                    call draw_patch_lines(painter,s,patch,0._dp,camera%scale,0)
                end associate
            end do
            do part=1,size(parts)
                if(.not.allocated(parts(part)%section_faces))cycle
                if(size(parts(part)%section_faces)==0)cycle
                call index_mesh_cell_selection(s,parts,selection)
                call draw_sections(.true.)
                exit
            end do
            call painter%restore()
            return
        end if
        call index_mesh_cell_selection(s,parts,selection)
        call viewport_geometry_bounds(objects,parts,[width,height],lower,upper)
        if(any(lower>=upper))return
        light=viewport_preview_raster(s)
        mesh_spacing=merge(9._dp,4._dp,motion)
        if(any(objects%kind==2))call index_smooth_joins(s,smooth_edges)
        capacity=0
        allocate(shading(6,size(objects)))
        do object=1,size(objects)
            associate(patch=>objects(object))
                if(.not.cad_fill(object).or.patch%kind==1)cycle
                if(patch%selected.and..not.s%layers(1))cycle
                if(patch%kind==2)then
                    capacity=capacity+product(patch%counts(:2)-1)
                else
                    counts=patch%boundary_counts-1
                    capacity=capacity+2*(counts(1)*counts(2)+(counts(1)+counts(2))*counts(3))
                end if
            end associate
        end do
        if(s%mesh_shaded)then
            do part=1,size(parts)
                if(allocated(parts(part)%section_faces))cycle
                counts=parts(part)%upper-parts(part)%lower
                select case(s%document%meshes(parts(part)%mesh)%kind)
                case(2); capacity=capacity+counts(1)*counts(2)
                case(3); capacity=capacity+2*(counts(1)*counts(2)+(counts(1)+counts(2))*counts(3))
                case default
                    cycle
                end select
            end do
        end if
        allocate(quads(capacity))
        used=0
        do object=1,size(objects)
            associate(patch=>objects(object))
                if(.not.cad_fill(object).or.patch%kind==1)cycle
                if(patch%selected.and..not.s%layers(1))cycle
                if(patch%kind==2)then
                    call stage_face(patch%points,patch,1,object)
                else
                    do side=1,6
                        call stage_face(patch%boundary(side)%points,patch,side,object)
                    end do
                end if
            end associate
        end do
        if(s%mesh_shaded)then
            allocate(selected_mesh(studio_mesh_count(s)))
            allocate(mesh_quads(2,size(parts)))
            do part=1,size(selected_mesh)
                selected_mesh(part)=mesh_in_selection(s,part)
            end do
            do part=1,size(parts)
                associate(piece=>parts(part))
                    overview=light.and.(motion.or.s%mesh_view/=mesh_lattice).and. &
                        .not.allocated(selection(part)%selected).and. &
                        .not.all(s%canvas_gesture%viewport_mesh_hover(:2)==[piece%mesh,piece%part])
                    mesh_quads(1,part)=used+1
                    call add_mesh_faces(quads,used,s,piece,part,selected_mesh(piece%mesh),overview,mesh_spacing,[width,height])
                    mesh_quads(2,part)=used
                end associate
            end do
        end if

        pixel_step=viewport_raster_step(s,upper-lower)
        if(pixel_step==1.and.s%interactive_quality.and..not.motion.and.viewport_navigation_active(s))then
            if(count(cad_fill.and.objects%kind>1)>=dense_preview_patches)pixel_step=2
        end if
        if(motion)pixel_step=pixel_step+1
        call buffer%initialize(upper(1)-lower(1),upper(2)-lower(2),lower(1),lower(2), &
            multisample=viewport_multisampling(s),pixel_step=pixel_step)
        if(.not.buffer%ready())return
        bias=.85_dp/camera%scale
        translucent_faces=any(objects%opacity<255).or.(size(parts)>0.and.s%mesh_view==mesh_lattice)
        threads=1
        !$ threads=min(4,omp_get_max_threads())
        if(threads>1.and.pixel_step==1.and.used>=parallel_face_count.and. &
            product(real(upper-lower,dp))>=parallel_pixel_count)then
            call buffer%cover_full_frame()
            band_height=(upper(2)-lower(2)+2*threads-1)/(2*threads)
            bands=(upper(2)-lower(2)+band_height-1)/band_height
            !$omp parallel do default(none) schedule(dynamic,1) num_threads(threads) shared(bands,band_height)
            do band=1,bands
                call paint_opaque_faces([(band-1)*band_height+1,band*band_height])
            end do
            !$omp end parallel do
        else
            call paint_opaque_faces()
        end if
        if(s%mesh_shaded)call draw_sections(.false.)
        do object=1,size(objects)
            if(cad_fill(object).and.objects(object)%kind==1)then
                call draw_patch_lines(painter,s,objects(object),0._dp,camera%scale,0,buffer)
            end if
        end do
        if(translucent_faces)then
            transparent_count=count([(quad_opacity(quads(primitive))<255,primitive=1,used)])
            allocate(transparent(transparent_count))
            transparent_index=0
            do primitive=1,used
                if(quad_opacity(quads(primitive))>=255)cycle
                transparent_index=transparent_index+1
                transparent(transparent_index)=primitive
            end do
            if(transparent_count>0)call sort_quad_indices(quads,transparent,1,transparent_count)
            do transparent_index=1,transparent_count
                primitive=transparent(transparent_index)
                call paint_quad_fill(quads(primitive),quad_opacity(quads(primitive)))
            end do
        end if
        if(allocated(mesh_quads))then
            do part=1,size(parts)
                if(mesh_quads(1,part)>mesh_quads(2,part))cycle
                call paint_mesh_edges(part,mesh_quads(1,part),mesh_quads(2,part))
            end do
        end if
        call draw_sections(.true.)
        do object=1,size(objects)
            if(.not.cad_fill(object).or.objects(object)%kind==1)cycle
            suppressed_edges=0
            if(objects(object)%kind==2.and.allocated(smooth_edges))suppressed_edges=smooth_edges(objects(object)%variant)
            call draw_patch_lines(painter,s,objects(object),bias,camera%scale,suppressed_edges,buffer)
        end do
        call buffer%paint(painter)
    contains

        subroutine paint_opaque_faces(rows)
            integer,optional,intent(in) :: rows(2)
            integer :: index

            do index=1,used
                if(translucent_faces)then
                    if(quad_opacity(quads(index))<255)cycle
                end if
                call paint_quad_fill(quads(index),255,rows)
            end do
        end subroutine paint_opaque_faces

        subroutine paint_mesh_edges(part,first,last)
            integer,intent(in) :: part,first,last
            integer :: index
            logical :: highlighted

            associate(piece=>parts(part))
                highlighted=allocated(selection(part)%selected).or. &
                    all(s%canvas_gesture%viewport_mesh_hover(:2)==[piece%mesh,piece%part])
                if(motion.and..not.highlighted)return
                if(.not.(highlighted.or.light))then
                    if(paint_shared_mesh_edges(piece,quads(first:last)))return
                end if
                do index=first,last
                    call paint_mesh_quad_edges(quads(index))
                end do
            end associate
        end subroutine paint_mesh_edges

        logical function paint_shared_mesh_edges(piece,faces) result(painted)
            type(projected_mesh_part),intent(in) :: piece
            type(screen_quad),intent(in) :: faces(:)
            integer :: strides(3),dim,nodes,axis,node,index,edge,first,last,status

            painted=.false.
            dim=s%document%meshes(piece%mesh)%kind
            strides=[1,piece%counts(1),product(piece%counts(:2))]
            nodes=size(piece%points,2)
            if(allocated(edge_uses))then
                if(size(edge_uses,2)<nodes)deallocate(edge_uses)
            end if
            if(.not.allocated(edge_uses))then
                allocate(edge_uses(3,nodes),stat=status)
                if(status/=0)return
            end if
            edge_uses(:,:nodes)=0_int8
            do index=1,size(faces)
                do edge=1,4
                    first=faces(index)%vertices(edge)
                    last=faces(index)%vertices(mod(edge,4)+1)
                    axis=findloc(strides(:dim),abs(last-first),dim=1)
                    if(axis==0)return
                    node=min(first,last)
                    edge_uses(axis,node)=min(4_int8,edge_uses(axis,node)+1_int8)
                end do
            end do
            call buffer%grid_edges(piece%points,edge_uses(:,:nodes),strides(:dim),mesh_ink,shared_mesh_opacity,.85_dp,bias)
            painted=.true.
        end function paint_shared_mesh_edges

        subroutine stage_face(points,patch,face,object_index)
            real(dp),intent(in) :: points(:,:)
            type(projected_patch),intent(in) :: patch
            integer,intent(in) :: face,object_index
            integer :: axes(2),face_counts(2),before

            call projected_face_axes(patch,face,axes,face_counts)
            before=used
            call append_face(quads,used,points,face_counts,[width,height],face,object_index,patch%selected)
            if(used>before)shading(face,object_index)%rgb=shade_face(s,patch,points,face_counts,face,camera%scale)
        end subroutine stage_face

        integer function quad_opacity(quad) result(opacity)
            type(screen_quad),intent(in) :: quad
            logical :: picked

            if(quad%mesh_part==0)then
                opacity=objects(quad%object)%opacity
            else
                associate(piece=>parts(quad%mesh_part))
                    picked=cell_selected(quad%mesh_part,quad%mesh_cell)
                    opacity=255
                    if(s%document%meshes(piece%mesh)%kind==3.and.s%mesh_view==mesh_lattice.and..not.picked)opacity=130
                end associate
            end if
        end function quad_opacity

        subroutine paint_quad_fill(quad,opacity,rows)
            type(screen_quad),intent(in) :: quad
            integer,intent(in) :: opacity
            integer,optional,intent(in) :: rows(2)
            real(dp) :: points(3,4),rgb(3,4)
            type(QColor) :: color
            logical :: picked

            if(quad%mesh_part==0)then
                associate(patch=>objects(quad%object))
                    if(patch%kind==3)then
                        points=patch%boundary(quad%face)%points(:,quad%vertices)
                    else
                        points=patch%points(:,quad%vertices)
                    end if
                end associate
            else
                points=parts(quad%mesh_part)%points(:,quad%vertices)
            end if
            if(present(rows))then
                if(maxval(points(2,:))<lower(2)+rows(1)-1.or.minval(points(2,:))>lower(2)+rows(2))return
            end if

            if(quad%mesh_part==0)then
                rgb=shading(quad%face,quad%object)%rgb(:,quad%vertices)
                call buffer%triangle(points(:,first_triangle),rgb(:,first_triangle),opacity,rows)
                call buffer%triangle(points(:,second_triangle),rgb(:,second_triangle),opacity,rows)
            else
                picked=cell_selected(quad%mesh_part,quad%mesh_cell)
                color=mesh_face_color(points,[1,2,3,4],camera%scale,quad%selected,picked,opacity)
                call buffer%flat_quad(points,real([color%r,color%g,color%b],dp),opacity,rows)
            end if
        end subroutine paint_quad_fill

        subroutine paint_mesh_quad_edges(quad)
            type(screen_quad),intent(in) :: quad
            real(dp) :: points(3,4),line_rgb(3),thickness
            integer :: step,next
            logical :: picked,hovered

            associate(piece=>parts(quad%mesh_part))
                picked=cell_selected(quad%mesh_part,quad%mesh_cell)
                hovered=all(s%canvas_gesture%viewport_mesh_hover==[piece%mesh,piece%part,quad%mesh_cell])
                if(motion.and..not.(picked.or.hovered))return
                points=piece%points(:,quad%vertices)
            end associate
            if(light.and..not.(hovered.or.picked).and.quad%mesh_cell>0)then
                if(maxval(maxval(points(:2,:),dim=2)-minval(points(:2,:),dim=2))<=3._dp*pixel_step)return
            end if
            line_rgb=mesh_ink
            thickness=.85_dp
            if(hovered)then
                line_rgb=hover_ink
                thickness=1.5_dp
            end if
            if(picked)then
                line_rgb=[252._dp,211._dp,255._dp]
                thickness=1.85_dp
            end if
            do step=1,4
                next=mod(step,4)+1
                call buffer%line(points(:,step),points(:,next),line_rgb,merge(245,180,picked),thickness,bias)
            end do
        end subroutine paint_mesh_quad_edges

        subroutine draw_sections(edges)
            logical,intent(in) :: edges
            real(dp) :: rgb(3),points(3,3),line_rgb(3),thickness
            type(QColor) :: color
            integer :: part,face,vertex,next
            logical :: selected,picked,cell_hovered

            do part=1,size(parts)
                if(.not.allocated(parts(part)%section_faces))cycle
                associate(piece=>parts(part))
                    selected=piece%mesh==s%selected_mesh
                    do face=1,size(piece%section_faces)
                        associate(poly=>piece%section_faces(face))
                            if(.not.poly%front.and..not.wireframe)cycle
                            picked=cell_selected(part,poly%cell)
                            cell_hovered=all(s%canvas_gesture%viewport_mesh_hover==[piece%mesh,piece%part,poly%cell])
                            if(edges)then
                                line_rgb=[208._dp,161._dp,100._dp]
                                thickness=.85_dp
                                if(poly%side==0)line_rgb=[103._dp,210._dp,205._dp]
                                if(cell_hovered)then
                                    line_rgb=hover_ink
                                    thickness=1.5_dp
                                end if
                                if(picked)then
                                    line_rgb=[252._dp,211._dp,255._dp]
                                    thickness=1.85_dp
                                end if
                                if(wireframe)call painter%setPen(QPen(QColor(nint(line_rgb(1)),nint(line_rgb(2)), &
                                    nint(line_rgb(3)),230),thickness))
                                do vertex=1,poly%count
                                    next=mod(vertex,poly%count)+1
                                    if(wireframe)then
                                        call painter%drawLine(poly%points(1,vertex),poly%points(2,vertex), &
                                            poly%points(1,next),poly%points(2,next))
                                    else
                                        call buffer%line(poly%points(:,vertex),poly%points(:,next),line_rgb,230,thickness,bias)
                                    end if
                                end do
                            else
                                color=mesh_section_color(poly,camera%scale,selected,picked,255)
                                rgb=real([color%r,color%g,color%b],dp)
                                do vertex=2,poly%count-1
                                    points=poly%points(:,[1,vertex,vertex+1])
                                    call buffer%triangle(points,spread(rgb,2,3),255)
                                end do
                            end if
                        end associate
                    end do
                end associate
            end do
        end subroutine draw_sections

        logical function cell_selected(part,cell) result(selected)
            integer,intent(in) :: part,cell

            selected=.false.
            if(.not.allocated(selection(part)%selected))return
            if(cell<1.or.cell>size(selection(part)%selected))return
            selected=selection(part)%selected(cell)
        end function cell_selected

    end subroutine paint_viewport_geometry

    subroutine index_mesh_cell_selection(s,parts,selection)
        type(studio_t),intent(in) :: s
        type(projected_mesh_part),intent(in) :: parts(:)
        type(mesh_cell_selection),intent(inout) :: selection(:)
        integer,allocatable :: mesh_ranges(:,:)
        integer :: i,mesh

        if(size(parts)==0)return
        allocate(mesh_ranges(2,maxval(parts%mesh)),source=0)
        do i=1,size(parts)
            mesh=parts(i)%mesh
            if(mesh_ranges(1,mesh)==0)mesh_ranges(1,mesh)=i
            mesh_ranges(2,mesh)=i
        end do
        if(allocated(s%selection%items))then
            do i=1,s%selection%count()
                associate(target=>s%selection%items(i))
                    if(target%category/=selection_mesh.or.target%item<1.or.target%part<1)cycle
                    call mark_cell(target%id,target%part,target%item)
                end associate
            end do
        else if(s%selected_mesh>0.and.all(s%mesh_cell(:2)>0))then
            call mark_cell(s%selected_mesh,s%mesh_cell(1),s%mesh_cell(2))
        end if
    contains
        subroutine mark_cell(mesh,part,cell)
            integer,intent(in) :: mesh,part,cell
            integer :: lower,upper,index,total

            if(mesh<1.or.mesh>size(mesh_ranges,2))return
            lower=mesh_ranges(1,mesh)
            upper=mesh_ranges(2,mesh)
            if(lower==0)return
            do while(lower<=upper)
                index=lower+(upper-lower)/2
                associate(current=>parts(index))
                    if(current%part==part)then
                        total=size(s%document%meshes(mesh)%parts(part)%cells,1)
                        if(cell>total)return
                        if(.not.allocated(selection(index)%selected))allocate(selection(index)%selected(total),source=.false.)
                        selection(index)%selected(cell)=.true.
                        return
                    end if
                    if(current%part<part)then
                        lower=index+1
                    else
                        upper=index-1
                    end if
                end associate
            end do
        end subroutine mark_cell
    end subroutine index_mesh_cell_selection

    subroutine viewport_geometry_bounds(objects,parts,canvas,lower,upper)
        type(projected_patch),intent(in) :: objects(:)
        type(projected_mesh_part),intent(in) :: parts(:)
        integer,intent(in) :: canvas(2)
        integer,intent(out) :: lower(2),upper(2)
        real(dp) :: lo(2),hi(2)
        integer :: i,side

        lo=real(canvas,dp)
        hi=0._dp
        do i=1,size(objects)
            associate(patch=>objects(i))
                if(patch%kind==3)then
                    do side=1,6
                        lo=min(lo,minval(patch%boundary(side)%points(:2,:),dim=2))
                        hi=max(hi,maxval(patch%boundary(side)%points(:2,:),dim=2))
                    end do
                else
                    lo=min(lo,minval(patch%points(:2,:),dim=2))
                    hi=max(hi,maxval(patch%points(:2,:),dim=2))
                end if
            end associate
        end do
        do i=1,size(parts)
            lo=min(lo,minval(parts(i)%points(:2,:),dim=2))
            hi=max(hi,maxval(parts(i)%points(:2,:),dim=2))
        end do
        lower=floor(min(real(canvas,dp),max(0._dp,lo-6._dp)))
        upper=ceiling(min(real(canvas,dp),max(0._dp,hi+6._dp)))
    end subroutine viewport_geometry_bounds

    subroutine append_face(quads,used,points,counts,canvas,side,object,selected)
        type(screen_quad),intent(inout) :: quads(:)
        integer,intent(inout) :: used
        real(dp),intent(in) :: points(:,:)
        integer,intent(in) :: counts(2),canvas(2),side,object
        logical,intent(in) :: selected
        integer :: i,j,base,vertices(4),width
        real(dp) :: x(4),y(4),right,bottom

        width=counts(1)
        right=real(canvas(1),dp)
        bottom=real(canvas(2),dp)
        do j=1,counts(2)-1
            do i=1,width-1
                base=i+(j-1)*width
                vertices=base+[0,1,width+1,width]
                x=points(1,vertices)
                if(maxval(x)<0._dp.or.minval(x)>right)cycle
                y=points(2,vertices)
                if(maxval(y)<0._dp.or.minval(y)>bottom)cycle
                used=used+1
                quads(used)=screen_quad(vertices,side,object,.25_dp*sum(points(3,vertices)),selected)
            end do
        end do
    end subroutine append_face

    function shade_face(s,patch,points,counts,side,scale) result(rgb)
        type(studio_t),intent(in) :: s
        type(projected_patch),intent(in) :: patch
        integer,intent(in) :: counts(2),side
        real(dp),intent(in) :: points(:,:),scale
        real(dp),allocatable :: rgb(:,:)
        real(dp),parameter :: main_light(3)=[-.38_dp,-.45_dp,.808_dp],fill_light(3)=[.78_dp,.32_dp,.54_dp], &
            halfway(3)=[-.204_dp,-.241_dp,.948_dp]
        real(dp) :: base(3),a(3),b(3),normal(3),length,light,gloss,color(3)
        integer :: i,j,vertex,fa,fb,sa,sb,width
        logical :: selected,hovered

        width=counts(1)
        allocate(rgb(3,product(counts)))
        base=[76._dp,128._dp,154._dp]
        if(patch%kind==3)base=[108._dp,130._dp,155._dp]
        if(all(patch%rgb>=0))base=real(patch%rgb,dp)
        selected=feature_selected(patch,pick_face,side)
        hovered=all(s%canvas_gesture%viewport_hover==[patch%kind,patch%variant,pick_face,side]).and..not.selected
        do j=1,counts(2)
            sa=merge(width,0,j>1)
            sb=merge(width,0,j<counts(2))
            do i=1,width
                fa=merge(1,0,i>1)
                fb=merge(1,0,i<width)
                vertex=i+(j-1)*width
                a=points(:,vertex+fb)-points(:,vertex-fa)
                b=points(:,vertex+sb)-points(:,vertex-sa)
                a(:2)=a(:2)/scale
                b(:2)=b(:2)/scale
                normal=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
                length=norm2(normal)
                if(length>tiny(1._dp))then
                    normal=sign(1._dp,normal(3))*normal/length
                else
                    normal=[0._dp,0._dp,1._dp]
                end if
                light=.56_dp+.40_dp*max(0._dp,dot_product(normal,main_light))+.13_dp*max(0._dp,dot_product(normal,fill_light))
                gloss=24._dp*max(0._dp,dot_product(normal,halfway))**16
                color=min(255._dp,light*base+gloss)
                if(hovered)color=.88_dp*color+.12_dp*hover_ink
                if(selected)color=.70_dp*color+.30_dp*selection_ink
                rgb(:,vertex)=color
            end do
        end do
    end function shade_face

    subroutine draw_patch_lines(painter,s,patch,bias,scale,suppressed_edges,buffer)
        type(QPainter),intent(inout) :: painter
        type(studio_t),intent(in) :: s
        type(projected_patch),intent(in) :: patch
        real(dp),intent(in) :: bias,scale
        integer,intent(in) :: suppressed_edges
        type(viewport_depth_t),intent(inout),optional :: buffer
        real(dp) :: ink(3),width,pen_ink(3),pen_width
        real(dp),allocatable :: edge_points(:,:)
        integer :: edge,step,sides(2),vertex,alpha,hover(2),pen_alpha
        logical :: selected,emphasized,prehighlight,edges,wireframe,closed(3)

        hover=0
        pen_ink=0._dp; pen_width=0._dp; pen_alpha=-1
        if(all(s%canvas_gesture%viewport_hover(:2)==[patch%kind,patch%variant]))hover=s%canvas_gesture%viewport_hover(3:)
        selected=(patch%selection_marked.and.patch%feature==0).or.patch%group_member
        wireframe=s%viewport_style==viewport_wireframe.or.(patch%selected.and..not.s%layers(1))
        edges=s%viewport_style/=viewport_shaded.or.wireframe
        if(patch%kind==1)then
            ink=[242._dp,195._dp,129._dp]
            if(all(patch%rgb>=0))ink=real(patch%rgb,dp)
            if(hover(1)>0)ink=hover_ink
            if(selected)ink=selection_ink
            width=merge(3._dp,2.1_dp,selected)
            if(patch%is_point)then
                call draw_path(patch%points(:,:1),ink,patch%opacity,merge(8._dp,5._dp,selected))
            else
                call draw_path(patch%points,ink,patch%opacity,width)
            end if
            return
        end if
        closed=.false.
        if(edges.and..not.wireframe)then
            closed(:patch%kind)=[(smooth_closed_axis(patch,step,scale),step=1,patch%kind)]
        end if
        do edge=1,edge_total(patch%kind)
            sides=patch_edge_boundaries(patch%kind,edge)
            emphasized=feature_selected(patch,pick_edge,edge)
            if(sides(1)>0)emphasized=emphasized.or.feature_selected(patch,pick_face,sides(1))
            if(sides(2)>0)emphasized=emphasized.or.feature_selected(patch,pick_face,sides(2))
            if(patch%kind==2)emphasized=emphasized.or.feature_selected(patch,pick_face,1)
            prehighlight=all(hover==[pick_edge,edge]).or.(hover(1)==pick_face.and.any(sides==hover(2)))
            if(patch%kind==2)prehighlight=prehighlight.or.hover(1)==pick_face
            if(.not.(edges.or.selected.or.prehighlight.or.emphasized))cycle
            if(.not.(wireframe.or.prehighlight.or.emphasized))then
                if(patch%kind==2)then
                    if(closed((edge+1)/2))cycle
                else
                    if(any(closed((sides+1)/2)))cycle
                end if
                if(.not.selected.and.btest(suppressed_edges,edge-1))cycle
            end if
            ink=edge_ink
            alpha=180
            width=.95_dp
            if(wireframe)ink=[165._dp,203._dp,222._dp]
            if(prehighlight)then
                ink=hover_ink; alpha=235; width=1.55_dp
            end if
            if(selected)then
                ink=selection_ink; alpha=215; width=1.45_dp
            end if
            if(emphasized)then
                ink=selection_ink; alpha=255; width=2.5_dp
            end if
            edge_points=projected_edge(patch,edge,sides)
            call draw_path(edge_points,ink,alpha,width)
        end do
        do vertex=1,2**patch%kind
            emphasized=feature_selected(patch,pick_vertex,vertex)
            if(.not.emphasized.and..not.all(hover==[pick_vertex,vertex]))cycle
            step=corner_sample(patch%counts,patch%kind,vertex)
            ink=hover_ink
            if(emphasized)ink=selection_ink
            call draw_path(patch%points(:,step:step),ink,255,merge(8._dp,6._dp,emphasized))
        end do
    contains
        subroutine draw_path(points,rgb,opacity,thickness)
            real(dp),intent(in) :: points(:,:),rgb(3),thickness
            integer,intent(in) :: opacity
            type(QColor) :: color
            integer :: i,n

            n=size(points,2)
            if(n==0)return
            if(present(buffer))then
                if(n==1)then
                    call buffer%line(points(:,1),points(:,1),rgb,opacity,thickness,bias)
                else
                    do i=1,n-1
                        call buffer%line(points(:,i),points(:,i+1),rgb,opacity,thickness,bias)
                    end do
                end if
            else
                color=QColor(nint(rgb(1)),nint(rgb(2)),nint(rgb(3)),opacity)
                if(n==1)then
                    call painter%setPen(QPen(style=Qt_NoPen))
                    call painter%setBrush(color)
                    call draw_dot(painter,points(:2,1),.5_dp*thickness)
                else
                    if(opacity/=pen_alpha.or.thickness/=pen_width.or.any(rgb/=pen_ink))then
                        call painter%setPen(QPen(color,thickness))
                        pen_ink=rgb; pen_width=thickness; pen_alpha=opacity
                    end if
                    if(n==2)then
                        call painter%drawLine(points(1,1),points(2,1),points(1,2),points(2,2))
                    else
                        call painter%drawPolyline(points(1,:),points(2,:))
                    end if
                end if
            end if
        end subroutine draw_path
    end subroutine draw_patch_lines

    pure logical function feature_selected(patch,feature,item) result(selected)
        type(projected_patch),intent(in) :: patch
        integer,intent(in) :: feature,item

        selected=(patch%feature==feature.and.patch%feature_number==item).or.btest(patch%selected_features(feature),item-1)
    end function feature_selected

    logical function smooth_closed_axis(patch,axis,scale) result(closed)
        type(projected_patch),intent(in) :: patch
        integer,intent(in) :: axis
        real(dp),intent(in) :: scale
        integer :: first,second,i,j,start,finish,stride(3)
        real(dp) :: a(3),b(3),separation(3),size(3),tolerance,length

        closed=.false.
        if(patch%counts(axis)<4)return
        first=mod(axis,3)+1
        second=mod(axis+1,3)+1
        stride=[1,patch%counts(1),product(patch%counts(:2))]
        size=maxval(patch%points,dim=2)-minval(patch%points,dim=2)
        size(:2)=size(:2)/scale
        tolerance=1.e-7_dp*maxval(size)
        do j=0,patch%counts(second)-1
            do i=0,patch%counts(first)-1
                start=1+i*stride(first)+j*stride(second)
                finish=start+(patch%counts(axis)-1)*stride(axis)
                separation=patch%points(:,start)-patch%points(:,finish)
                separation(:2)=separation(:2)/scale
                if(norm2(separation)>tolerance)then
                    closed=.false.
                    return
                end if
                a=patch%points(:,start+stride(axis))-patch%points(:,start)
                b=patch%points(:,finish)-patch%points(:,finish-stride(axis))
                a(:2)=a(:2)/scale; b(:2)=b(:2)/scale
                length=norm2(a)*norm2(b)
                if(length<=tiny(1._dp))cycle
                if(dot_product(a,b)<.92_dp*length)then
                    closed=.false.
                    return
                end if
                closed=.true.
            end do
        end do
    end function smooth_closed_axis

    subroutine index_smooth_joins(s,edges)
        type(studio_t),intent(in) :: s
        integer,allocatable,intent(out) :: edges(:)
        integer :: group,i

        if(.not.allocated(s%document%scene%catalog(2)%entries))return
        if(.not.allocated(s%document%multipatches%groups))return
        allocate(edges(size(s%document%scene%catalog(2)%entries)),source=0)
        do group=1,size(s%document%multipatches%groups)
            associate(assembly=>s%document%multipatches%groups(group))
                if(assembly%definition%kind/=2.or..not.allocated(assembly%definition%joins))cycle
                if(.not.allocated(assembly%validation%joins))cycle
                do i=1,min(size(assembly%definition%joins),size(assembly%validation%joins))
                    associate(join=>assembly%definition%joins(i))
                        if(join%continuity<1.or..not.assembly%validation%joins(i)%valid)cycle
                        call suppress_edge(join%patch_a,join%side_a,join%patch_b)
                        call suppress_edge(join%patch_b,join%side_b,join%patch_a)
                    end associate
                end do
            end associate
        end do
    contains
        subroutine suppress_edge(id,side,neighbor)
            integer,intent(in) :: id,side,neighbor

            if(id<1.or.id>size(edges).or.side<1.or.side>4)return
            if(scene_patch_drawn(s,2,neighbor))edges(id)=ibset(edges(id),side-1)
        end subroutine suppress_edge
    end subroutine index_smooth_joins

    recursive subroutine sort_quad_indices(quads,indices,first,last)
        type(screen_quad),intent(in) :: quads(:)
        integer,intent(inout) :: indices(:)
        integer,intent(in) :: first,last
        real(dp) :: pivot
        integer :: i,j,swapped

        if(first>=last)return
        pivot=quads(indices(first+(last-first)/2))%depth
        i=first; j=last
        do
            do while(quads(indices(i))%depth<pivot); i=i+1; end do
            do while(quads(indices(j))%depth>pivot); j=j-1; end do
            if(i>j)exit
            swapped=indices(i); indices(i)=indices(j); indices(j)=swapped
            i=i+1; j=j-1
            if(i>j)exit
        end do
        if(first<j)call sort_quad_indices(quads,indices,first,j)
        if(i<last)call sort_quad_indices(quads,indices,i,last)
    end subroutine sort_quad_indices

end module studio_viewport_geometry
