module studio_viewport_scene
    use studio_display_cache, only: display_grid
    use studio_display_refinement, only: display_grid_current
    use studio_commands, only: parameter_lines_layer
    use studio_forcad, only: closest_point_on_patch
    use studio_selection, only: pick_automatic, pick_face, pick_edge, pick_vertex, pick_control_points, pick_multipatch, &
        patch_edge_boundaries, selection_target_t, selection_patch,selection_control,selection_assembly
    use studio_state, only: viewport_wireframe, studio_t, dp, scene_patch_drawn, scene_drawn_patch_count, &
        clear_control_dependency_preview
    use studio_multipatch_controls, only: multipatch_controls_t, multipatch_control_preview_t, &
        read_document_multipatch_controls, joined_patch_assembly, assembly_control_id, &
        preview_document_multipatch_controls, set_multipatch_control_displacement
    use studio_multipatch_collection, only: count_multipatches,multipatch_contains
    use studio_single_patch, only: copy_single_patch_geometry,single_patch_scene
    use studio_scene, only: get_studio_material
    use studio_viewport_camera, only: viewport_camera, prepare_camera, project_vertex, workspace_view_hit, &
        in_drawable_view,project_viewport_box,viewport_clip_intersects,viewport_preview_raster
    use studio_viewport_mesh, only: mesh_replaces_patch_fill,mesh_covered_patches,pick_mesh_cell
    use studio_viewport_primitives, only: project_points, nearest_curve_parameter, screen_triangle, &
        screen_triangle_nearest, pick_sample_grid
    implicit none(type, external)
    private

    public :: projected_patch, project_scene, active_patch_visible, has_controls, has_slice, show_volume_slice
    public :: pick_scene_patch, pick_live_viewport, begin_live_control_drag, update_live_control_drag, viewport_feature_text
    public :: corner_sample, edge_total, edge_samples, pick_control_point, opaque_patch_depth
    public :: projected_control_patch, project_selected_controls, pick_selected_control, selected_controls_available
    public :: viewport_patch_target, multipatch_control_scope, project_shared_controls, pick_shared_control
    public :: shared_control_position,viewport_control_position, project_control_dependencies
    public :: projected_edge, projected_face_axes
    public :: pick_automatic, pick_face, pick_edge, pick_vertex, pick_control_points

    real(dp), parameter :: control_hit_radius=6._dp, control_drag_slop=3._dp
    real(dp), parameter :: edge_pick_radius=7._dp, vertex_pick_radius=9._dp, snap_pick_radius=14._dp
    type :: projected_patch
        real(dp), allocatable :: points(:,:)
        type(display_grid) :: boundary(6)
        integer :: boundary_counts(3)=1
        integer :: kind, variant, counts(3)
        integer :: rgb(3),opacity,feature,feature_number
        integer :: selected_features(3)=0
        logical :: selected, selection_marked, group_member, is_point
    end type projected_patch

    type :: projected_control_patch
        real(dp),allocatable :: points(:,:)
        integer,allocatable :: shared_ids(:)
        integer :: kind,id,counts(3)
        logical :: focused
    end type projected_control_patch

    type :: patch_display_selection
        logical :: patch=.false.,member=.false.
        integer :: features(3)=0
    end type patch_display_selection

contains

    subroutine project_scene(s,camera,objects,preview,clip,display,highlights,covered,source_selected)
        type(studio_t), intent(in) :: s
        type(viewport_camera), intent(in) :: camera
        type(projected_patch), allocatable, intent(out) :: objects(:)
        logical,intent(in),optional :: preview
        real(dp),intent(in),optional :: clip(4)
        logical,intent(in),optional :: display,highlights
        logical,intent(in),optional :: covered(:,:)
        logical,allocatable,intent(out),optional :: source_selected(:,:)
        integer,allocatable :: targets(:,:)
        real(dp),allocatable :: screen_spans(:)
        type(patch_display_selection),allocatable :: selection(:,:)
        integer :: kind,variant,n,side,i
        real(dp) :: bounds(4),span,spacing
        logical :: moving,refined,preview_active,rendering,mark_selection,light

        preview_active=.false.
        if(present(preview))preview_active=preview
        rendering=.false.
        if(present(display))rendering=display
        mark_selection=.true.
        if(present(highlights))mark_selection=highlights
        light=rendering.and.(s%viewport_style==viewport_wireframe.or.viewport_preview_raster(s))
        spacing=merge(13._dp,6._dp,light)
        if(present(covered).or.present(source_selected))call build_patch_selection(s,selection)
        if(present(source_selected))source_selected=selection%patch.or.selection%member
        allocate(targets(2,scene_drawn_patch_count(s)))
        allocate(screen_spans(size(targets,2)))
        n=0
        if(.not.s%mesh_only)then
            do kind=1,size(s%document%scene%catalog)
                if(.not.allocated(s%document%scene%catalog(kind)%entries))cycle
                do variant=1,size(s%document%scene%catalog(kind)%entries)
                    if(.not.scene_patch_drawn(s,kind,variant))cycle
                    if(present(covered))then
                        if(covered(variant,kind))then
                            if(s%live_scene%dimension/=kind.or.s%live_scene%variant/=variant)then
                                associate(mark=>selection(variant,kind))
                                    if(.not.any(mark%features/=0))then
                                        if(.not.(rendering.and.s%layers(parameter_lines_layer).and.(mark%patch.or.mark%member)))cycle
                                    end if
                                end associate
                            end if
                        end if
                    end if
                    span=huge(1._dp)
                    if((present(clip).or.rendering).and..not.(preview_active.or.s%control_dragging))then
                        associate(patch=>s%document%scene%catalog(kind)%entries(variant)%patch)
                            if(allocated(patch%controls))then
                                bounds=project_viewport_box(camera,minval(patch%controls,dim=2),maxval(patch%controls,dim=2))
                                if(present(clip))then
                                    if(.not.viewport_clip_intersects(bounds,clip))cycle
                                end if
                                if(rendering)span=maxval(bounds(3:)-bounds(:2))
                            end if
                        end associate
                    end if
                    n=n+1
                    targets(:,n)=[kind,variant]
                    screen_spans(n)=span
                end do
            end do
        end if
        allocate(objects(n))
        if(n==0)return
        if((mark_selection.or.preview_active).and..not.allocated(selection))call build_patch_selection(s,selection)
        do i=1,n
            kind=targets(1,i)
            variant=targets(2,i)
            associate(object=>objects(i),entry=>s%document%scene%catalog(kind)%entries(variant),display=>s%display_cache%ranks(kind))
                object%kind=kind
                object%variant=variant
                object%selected=s%selected_multipatch==0.and.s%live_scene%valid.and.s%live_scene%dimension==kind.and. &
                    s%live_scene%variant==variant
                object%selection_marked=.false.
                object%group_member=.false.
                if(allocated(selection))then
                    object%selection_marked=selection(variant,kind)%patch
                    object%group_member=selection(variant,kind)%member
                    object%selected_features=selection(variant,kind)%features
                end if
                object%counts=entry%patch%sample_counts
                object%is_point=kind==1.and.entry%patch%control_counts(1)==1
                object%rgb=entry%rgb
                object%opacity=entry%opacity
                object%feature=0
                object%feature_number=0
                if(object%selected.and.all(s%viewport_feature(:2)==[kind,variant]))then
                    object%feature=s%viewport_feature(3)
                    object%feature_number=s%viewport_feature(4)
                end if
                moving=preview_active.and.(object%selected.or.(s%selected_multipatch>0.and.object%group_member))
                refined=.not.light.and.display_grid_current(s,kind,variant)
                if(refined.and.kind<3)then
                    object%counts=display%entries(variant)%counts
                    object%points=project_points(s,camera,display%entries(variant)%grids(1)%points,moving)
                else if(rendering.and..not.moving.and..not.refined.and.kind>1)then
                    object%counts=display_grid_counts(entry%patch%points,object%counts,screen_spans(i),spacing)
                    object%points=project_display_grid(s,camera,entry%patch%points,entry%patch%sample_counts,object%counts)
                else
                    object%points=project_points(s,camera,entry%patch%points,moving)
                end if
                if(kind/=3)cycle
                object%boundary_counts=object%counts
                if(refined)object%boundary_counts=display%entries(variant)%counts
                do side=1,6
                    if(refined)then
                        object%boundary(side)%points=project_points(s,camera,display%entries(variant)%grids(side)%points,moving)
                    else
                        call copy_coarse_volume_face(object%points,object%counts,side,object%boundary(side)%points)
                    end if
                end do
            end associate
        end do
    end subroutine project_scene

    subroutine build_patch_selection(s,selection)
        type(studio_t),intent(in) :: s
        type(patch_display_selection),allocatable,intent(out) :: selection(:,:)
        integer :: rank,maximum,i

        maximum=0
        do rank=1,size(s%document%scene%catalog)
            if(allocated(s%document%scene%catalog(rank)%entries))then
                maximum=max(maximum,size(s%document%scene%catalog(rank)%entries))
            end if
        end do
        allocate(selection(maximum,3))
        if(.not.allocated(s%selection%items))then
            call add_group(s%selected_multipatch)
            if(s%selected_multipatch/=0.or.s%selected_mesh/=0.or..not.s%live_scene%valid)return
            rank=s%live_scene%dimension
            i=s%live_scene%variant
            if(i>=1.and.i<=maximum.and.rank>=1.and.rank<=3)selection(i,rank)%patch=.true.
            return
        end if
        do i=1,size(s%selection%items)
            associate(target=>s%selection%items(i))
                select case(target%category)
                case(selection_assembly)
                    call add_group(target%id)
                case(selection_patch,selection_control)
                    if(target%id<1.or.target%id>maximum.or.target%rank<1.or.target%rank>3)cycle
                    associate(entry=>selection(target%id,target%rank))
                        entry%patch=.true.
                        if(target%category/=selection_patch)cycle
                        if(target%feature<pick_face.or.target%feature>pick_vertex)cycle
                        if(target%item<1.or.target%item>bit_size(0))cycle
                        entry%features(target%feature)=ibset(entry%features(target%feature),target%item-1)
                    end associate
                case default
                    continue
                end select
            end associate
        end do
    contains
        subroutine add_group(id)
            integer,intent(in) :: id
            integer :: j,member

            if(id<1.or.id>count_multipatches(s%document%multipatches))return
            associate(group=>s%document%multipatches%groups(id)%definition)
                do j=1,size(group%members)
                    member=group%members(j)
                    if(member<1.or.member>maximum)cycle
                    selection(member,group%kind)%member=.true.
                end do
            end associate
        end subroutine add_group
    end subroutine build_patch_selection

    function viewport_patch_target(s,kind,variant) result(target)
        type(studio_t),intent(in) :: s
        integer,intent(in) :: kind,variant
        type(selection_target_t) :: target
        integer :: group,owner

        if(kind<1.or.variant<1)return
        target=selection_target_t(selection_patch,variant,kind)
        if(s%feature_pick_mode/=pick_multipatch)return
        owner=0
        do group=1,count_multipatches(s%document%multipatches)
            if(.not.multipatch_contains(s%document%multipatches,group,kind,variant))cycle
            if(owner>0)return
            owner=group
        end do
        if(owner>0)target=selection_target_t(selection_assembly,owner)
    end function viewport_patch_target

    logical function multipatch_control_scope(s) result(active)
        type(studio_t),intent(in) :: s

        active=.false.
        if(s%selected_multipatch<1.or.s%selected_multipatch>count_multipatches(s%document%multipatches))return
        active=.true.
        if(.not.allocated(s%selection%items))return
        active=s%selection%count()==1.and.s%selection%contains(selection_target_t(selection_assembly,s%selected_multipatch))
    end function multipatch_control_scope

    logical function shared_control_position(s,point,position) result(available)
        type(studio_t),intent(in) :: s
        integer,intent(in) :: point
        real(dp),intent(out) :: position(3)

        available=.false.
        position=0._dp
        if(.not.multipatch_control_scope(s).or..not.s%shared_controls%valid)return
        if(point<1.or.point>size(s%shared_controls%points,2))return
        position=s%shared_controls%points(:,point)
        available=.true.
    end function shared_control_position

    logical function viewport_control_position(s,point,position) result(available)
        type(studio_t),intent(in) :: s
        integer,intent(in) :: point
        real(dp),intent(out) :: position(3)

        if(multipatch_control_scope(s))then
            available=shared_control_position(s,point,position)
            return
        end if
        available=.false.
        position=0._dp
        if(.not.has_controls(s))return
        if(point<1.or.point>size(s%live_scene%controls,2))return
        position=s%live_scene%controls(:,point)
        available=.true.
    end function viewport_control_position

    logical function selected_controls_available(s) result(available)
        type(studio_t),intent(in) :: s
        type(patch_display_selection),allocatable :: selection(:,:)
        integer :: kind,id

        available=.false.
        if(s%mesh_only)return
        call build_patch_selection(s,selection)
        do kind=1,size(s%document%scene%catalog)
            if(.not.allocated(s%document%scene%catalog(kind)%entries))cycle
            do id=1,size(s%document%scene%catalog(kind)%entries)
                if(.not.(selection(id,kind)%patch.or.selection(id,kind)%member))cycle
                if(.not.scene_patch_drawn(s,kind,id))cycle
                if(.not.valid_controls(s%document%scene%catalog(kind)%entries(id)%patch))cycle
                available=.true.
                return
            end do
        end do
    end function selected_controls_available

    subroutine project_selected_controls(s,camera,clip,patches,preview,mesh_coverage)
        type(studio_t),intent(in) :: s
        type(viewport_camera),intent(in) :: camera
        real(dp),intent(in) :: clip(4)
        type(projected_control_patch),allocatable,intent(out) :: patches(:)
        logical,intent(in),optional :: preview
        logical,target,intent(in),optional :: mesh_coverage(:,:)
        type(patch_display_selection),allocatable :: selection(:,:)
        integer,allocatable :: targets(:,:)
        logical,allocatable,target :: coverage(:,:)
        logical,pointer :: covered(:,:)
        logical,allocatable :: moved(:)
        integer,allocatable :: dependency_rows(:)
        integer :: kind,id,n,i,member,first,last,j
        logical :: moving,focused,with_preview,shared
        real(dp) :: displacement(3)

        if(s%mesh_only.or..not.any(s%layers(3:4)))then
            allocate(patches(0))
            return
        end if
        with_preview=.false.
        if(present(preview))with_preview=preview
        shared=multipatch_control_scope(s).and.s%shared_controls%valid
        if(shared.and.present(preview))call shared_control_preview(s,moved,displacement)
        if(present(preview))call control_dependency_rows(s,dependency_rows)
        nullify(covered)
        if(s%mesh_section_axis>0)then
            if(present(mesh_coverage))then
                covered=>mesh_coverage
            else
                call mesh_covered_patches(s,coverage)
                covered=>coverage
            end if
        end if
        call build_patch_selection(s,selection)
        allocate(targets(2,count(selection%patch.or.selection%member)))
        n=0
        do kind=1,size(s%document%scene%catalog)
            if(.not.allocated(s%document%scene%catalog(kind)%entries))cycle
            do id=1,size(s%document%scene%catalog(kind)%entries)
                associate(mark=>selection(id,kind),patch=>s%document%scene%catalog(kind)%entries(id)%patch)
                    if(.not.(mark%patch.or.mark%member))cycle
                    if(.not.scene_patch_drawn(s,kind,id).or..not.valid_controls(patch))cycle
                    if(kind==3.and.associated(covered))then
                        if(covered(id,kind))cycle
                    end if
                    if(.not.(with_preview.or.s%control_dragging))then
                        if(.not.viewport_clip_intersects(project_viewport_box(camera, &
                            minval(patch%controls,dim=2),maxval(patch%controls,dim=2)),clip))cycle
                    end if
                    n=n+1
                    targets(:,n)=[kind,id]
                end associate
            end do
        end do
        allocate(patches(n))
        do i=1,n
            kind=targets(1,i)
            id=targets(2,i)
            associate(patch=>s%document%scene%catalog(kind)%entries(id)%patch,controls=>patches(i))
                focused=s%selected_multipatch==0.and.s%live_scene%valid.and. &
                    s%live_scene%dimension==kind.and.s%live_scene%variant==id
                moving=with_preview.and.(focused.or.(s%selected_multipatch>0.and.selection(id,kind)%member))
                controls%kind=kind
                controls%id=id
                controls%counts=patch%control_counts
                controls%focused=focused
                if(focused.and.has_controls(s))then
                    controls%points=project_points(s,camera,s%live_scene%controls,moving)
                    controls%counts=s%live_scene%control_counts
                else
                    controls%points=project_points(s,camera,patch%controls,moving)
                end if
                if(shared.and.selection(id,kind)%member)then
                    member=findloc(s%shared_controls%members,id,dim=1)
                    if(member>0)then
                        first=s%shared_controls%offsets(member)+1
                        last=s%shared_controls%offsets(member+1)
                        if(last-first+1==size(controls%points,2))controls%shared_ids=s%shared_controls%map(first:last)
                    end if
                end if
                if(.not.present(preview).or..not.s%control_dragging)cycle
                if(allocated(controls%shared_ids))then
                    if(allocated(dependency_rows))call project_patch_dependency_positions( &
                        s,camera,s%shared_controls,controls%shared_ids,patch%controls,dependency_rows,controls%points)
                    if(.not.allocated(moved))cycle
                    do j=1,size(controls%shared_ids)
                        if(.not.moved(controls%shared_ids(j)))cycle
                        if(controls%shared_ids(j)==s%selected_control)then
                            controls%points(:,j)=project_vertex(camera,s%control_preview)
                        else
                            controls%points(:,j)=project_vertex(camera,s%shared_controls%points(:,controls%shared_ids(j))+displacement)
                        end if
                    end do
                else if(focused)then
                    if(allocated(dependency_rows))then
                        associate(base=>s%control_dependency_source%member_controls)
                            if(base%valid)then
                                member=findloc(base%members,id,dim=1)
                                if(member>0)then
                                    first=base%offsets(member)+1
                                    last=base%offsets(member+1)
                                    if(last-first+1==size(controls%points,2))call project_patch_dependency_positions( &
                                        s,camera,base,base%map(first:last),s%live_scene%controls,dependency_rows,controls%points)
                                end if
                            end if
                        end associate
                    end if
                    if(s%selected_control>=1.and.s%selected_control<=size(controls%points,2))then
                        controls%points(:,s%selected_control)=project_vertex(camera,s%control_preview)
                    end if
                end if
            end associate
        end do
    end subroutine project_selected_controls

    subroutine project_shared_controls(s,camera,patches,points,visible,preview)
        type(studio_t),intent(in) :: s
        type(viewport_camera),intent(in) :: camera
        type(projected_control_patch),intent(in) :: patches(:)
        real(dp),allocatable,intent(out) :: points(:,:)
        logical,allocatable,intent(out) :: visible(:)
        logical,intent(in),optional :: preview
        logical,allocatable :: moved(:)
        real(dp) :: displacement(3)
        integer :: i,id

        if(.not.multipatch_control_scope(s).or..not.s%shared_controls%valid)then
            allocate(points(3,0),visible(0))
            return
        end if
        call visible_shared_controls(s,patches,visible)
        points=project_points(s,camera,s%shared_controls%points,preview)
        if(.not.present(preview))return
        call shared_control_preview(s,moved,displacement)
        if(allocated(moved))then
            do i=1,size(visible)
                if(.not.(moved(i).and.visible(i)))cycle
                if(i==s%selected_control)then
                    points(:,i)=project_vertex(camera,s%control_preview)
                else
                    points(:,i)=project_vertex(camera,s%shared_controls%points(:,i)+displacement)
                end if
            end do
        end if
        if(.not.control_dependencies_available(s))return
        associate(dependent=>s%control_dependency_preview)
            do i=1,size(dependent%ids)
                id=dependent%ids(i)
                if(visible(id))points(:,id)=project_vertex(camera,dependent%points(:,i))
            end do
        end associate
    end subroutine project_shared_controls

    logical function control_dependencies_available(s) result(available)
        type(studio_t), intent(in) :: s

        available=.false.
        if(.not.s%control_dragging.or.s%control_dependency_source%assembly<=0)return
        if(.not.allocated(s%control_dependency_preview%ids))return
        available=size(s%control_dependency_preview%ids)>0
    end function control_dependencies_available

    subroutine control_dependency_rows(s,rows)
        type(studio_t), intent(in) :: s
        integer, allocatable, intent(out) :: rows(:)
        integer :: i,control_count

        if(.not.control_dependencies_available(s))return
        if(s%control_dependency_source%member_controls%valid)then
            control_count=size(s%control_dependency_source%member_controls%weights)
        else
            if(.not.s%shared_controls%valid)return
            control_count=size(s%shared_controls%weights)
        end if
        allocate(rows(control_count),source=0)
        do i=1,size(s%control_dependency_preview%ids)
            rows(s%control_dependency_preview%ids(i))=i
        end do
    end subroutine control_dependency_rows

    subroutine project_patch_dependency_positions(s,camera,base,ids,original,rows,projected)
        type(studio_t), intent(in) :: s
        type(viewport_camera), intent(in) :: camera
        type(multipatch_controls_t), intent(in) :: base
        integer, intent(in) :: ids(:),rows(:)
        real(dp), intent(in) :: original(:,:)
        real(dp), intent(inout) :: projected(:,:)
        integer :: i,id,row

        do i=1,size(ids)
            id=ids(i)
            row=rows(id)
            if(row==0)cycle
            if(.not.s%control_dependency_preview%moved(row))cycle
            projected(:,i)=project_vertex(camera,original(:,i)+(s%control_dependency_preview%points(:,row)-base%points(:,id)))
        end do
    end subroutine project_patch_dependency_positions

    subroutine project_control_dependencies(s,camera,shared_visible,mesh_coverage,points,visible)
        type(studio_t), intent(in) :: s
        type(viewport_camera), intent(in) :: camera
        logical, intent(in) :: shared_visible(:),mesh_coverage(:,:)
        real(dp), allocatable, intent(out) :: points(:,:)
        logical, allocatable, intent(out) :: visible(:)
        integer, allocatable :: rows(:)
        integer :: member,patch,rank,i,row

        if(.not.control_dependencies_available(s))then
            allocate(points(3,0),visible(0))
            return
        end if
        associate(dependent=>s%control_dependency_preview,base=>s%control_dependency_source%member_controls)
            allocate(visible(size(dependent%ids)),source=.false.)
            points=project_points(s,camera,dependent%points)
            if(.not.base%valid)then
                if(size(shared_visible)>0)visible=shared_visible(dependent%ids)
                return
            end if
            call control_dependency_rows(s,rows)
            rank=base%kind
            do member=1,size(base%members)
                patch=base%members(member)
                if(.not.scene_patch_drawn(s,rank,patch))cycle
                if(rank==3.and.s%mesh_section_axis>0)then
                    if(mesh_coverage(patch,rank))cycle
                end if
                do i=base%offsets(member)+1,base%offsets(member+1)
                    row=rows(base%map(i))
                    if(row>0)visible(row)=.true.
                end do
            end do
        end associate
    end subroutine project_control_dependencies

    subroutine visible_shared_controls(s,patches,visible)
        type(studio_t),intent(in) :: s
        type(projected_control_patch),intent(in) :: patches(:)
        logical,allocatable,intent(out) :: visible(:)
        integer :: i,j

        allocate(visible(size(s%shared_controls%points,2)),source=.false.)
        do i=1,size(patches)
            if(.not.allocated(patches(i)%shared_ids))cycle
            do j=1,size(patches(i)%shared_ids)
                visible(patches(i)%shared_ids(j))=.true.
            end do
        end do
    end subroutine visible_shared_controls

    subroutine shared_control_preview(s,moved,displacement)
        type(studio_t),intent(in) :: s
        logical,allocatable,intent(out) :: moved(:)
        real(dp),intent(out) :: displacement(3)
        integer :: i,point

        displacement=0._dp
        if(.not.s%control_dragging)return
        if(.not.shared_control_position(s,s%selected_control,displacement))return
        displacement=s%control_preview-displacement
        allocate(moved(size(s%shared_controls%points,2)),source=.false.)
        moved(s%selected_control)=.true.
        if(.not.allocated(s%selected_shared_controls))return
        if(.not.any(s%selected_shared_controls==s%selected_control))return
        do i=1,size(s%selected_shared_controls)
            point=s%selected_shared_controls(i)
            if(point>=1.and.point<=size(moved))moved(point)=.true.
        end do
    end subroutine shared_control_preview

    pure function display_grid_counts(points,counts,span,spacing) result(samples)
        real(dp),intent(in) :: points(:,:)
        integer,intent(in) :: counts(3)
        real(dp),intent(in) :: span,spacing
        integer :: samples(3),stride(3),middle(3),node,first,last,axis
        real(dp) :: tolerance

        samples=counts
        if(span>=spacing*real(maxval(counts)-1,dp))return
        if(span>=0._dp)samples=min(counts,max(2,ceiling(span/spacing)+1))
        if(all(samples==counts))return
        stride=[1,counts(1),product(counts(:2))]
        middle=(counts-1)/2
        node=1+dot_product(middle,stride)
        do axis=1,3
            if(samples(axis)==counts(axis).or.samples(axis)>=5)cycle
            first=node-middle(axis)*stride(axis)
            last=first+(counts(axis)-1)*stride(axis)
            tolerance=32._dp*epsilon(1._dp)*max(1._dp,maxval(abs(points(:,first))),maxval(abs(points(:,last))))
            if(norm2(points(:,first)-points(:,last))<=tolerance)samples(axis)=min(5,counts(axis))
        end do
    end function display_grid_counts

    function project_display_grid(s,camera,points,counts,samples) result(projected)
        type(studio_t),intent(in) :: s
        type(viewport_camera),intent(in) :: camera
        real(dp),intent(in) :: points(:,:)
        integer,intent(in) :: counts(3),samples(3)
        real(dp),allocatable :: projected(:,:)
        integer :: i,j,k,index,node
        real(dp) :: step(3)

        if(all(counts==samples))then
            projected=project_points(s,camera,points)
            return
        end if
        allocate(projected(3,product(samples)))
        step=real(counts-1,dp)/real(max(1,samples-1),dp)
        index=0
        do k=0,samples(3)-1
            do j=0,samples(2)-1
                do i=0,samples(1)-1
                    node=1+nint(i*step(1))+counts(1)*(nint(j*step(2))+counts(2)*nint(k*step(3)))
                    index=index+1
                    projected(:,index)=project_vertex(camera,points(:,node))
                end do
            end do
        end do
    end function project_display_grid

    pure subroutine projected_face_axes(object,side,axes,counts)
        type(projected_patch),intent(in) :: object
        integer,intent(in) :: side
        integer,intent(out) :: axes(2),counts(2)
        integer :: axis

        if(object%kind==3)then
            axis=(side+1)/2
            axes=[mod(axis,3)+1,mod(axis+1,3)+1]
            counts=object%boundary_counts(axes)
        else
            axes=[1,2]
            counts=object%counts(:2)
        end if
    end subroutine projected_face_axes

    pure subroutine copy_coarse_volume_face(volume,volume_counts,side,points)
        real(dp), intent(in) :: volume(:,:)
        integer,intent(in) :: volume_counts(3),side
        real(dp), allocatable, intent(out) :: points(:,:)
        integer :: axes(2),counts(2),axis,stride(3),base,i,j,row,node

        axis=(side+1)/2
        axes=[mod(axis,3)+1,mod(axis+1,3)+1]
        counts=volume_counts(axes)
        stride=[1,volume_counts(1),product(volume_counts(:2))]
        base=1+mod(side+1,2)*(volume_counts(axis)-1)*stride(axis)
        allocate(points(3,product(counts)))
        do j=0,counts(2)-1
            do i=0,counts(1)-1
                row=1+i+j*counts(1)
                node=base+i*stride(axes(1))+j*stride(axes(2))
                points(:,row)=volume(:,node)
            end do
        end do
    end subroutine copy_coarse_volume_face

    function projected_edge(object,edge,sides) result(points)
        type(projected_patch),intent(in) :: object
        integer,intent(in) :: edge
        integer,intent(out) :: sides(2)
        real(dp),allocatable :: points(:,:)
        integer :: counts(3),n,nodes(max(maxval(object%counts),maxval(object%boundary_counts)))
        integer :: axes(2),face_counts(2),stride(3),rows(size(nodes)),i

        counts=object%counts
        if(object%kind==3)counts=object%boundary_counts
        call edge_samples(counts,object%kind,edge,nodes,n,sides)
        if(object%kind/=3.or.n==0)then
            points=object%points(:,nodes(:n))
            return
        end if
        call projected_face_axes(object,sides(1),axes,face_counts)
        stride=[1,counts(1),product(counts(:2))]
        do i=1,n
            rows(i)=1+mod((nodes(i)-1)/stride(axes(1)),face_counts(1))+ &
                face_counts(1)*mod((nodes(i)-1)/stride(axes(2)),face_counts(2))
        end do
        points=object%boundary(sides(1))%points(:,rows(:n))
    end function projected_edge

    logical function show_volume_slice(s) result(shown)
        type(studio_t),intent(in) :: s
        integer :: rgb(3),opacity

        shown=.false.
        if(s%selected_multipatch>0.or.s%live_scene%dimension/=3)return
        shown=s%layers(2).or..not.s%layers(1)
        if(shown)return
        call get_studio_material(s%document%scene,3,s%live_scene%variant,rgb,opacity)
        shown=opacity<255
    end function show_volume_slice

    logical function active_patch_visible(s) result(visible)
        type(studio_t), intent(in) :: s

        visible=s%live_scene%valid.and.scene_patch_drawn(s,s%live_scene%dimension,s%live_scene%variant)
    end function active_patch_visible

    logical function has_controls(s) result(available)
        type(studio_t), intent(in) :: s

        available=valid_controls(s%live_scene)
    end function has_controls

    logical function valid_controls(patch) result(available)
        type(single_patch_scene),intent(in) :: patch

        available=.false.
        if(.not.allocated(patch%controls))return
        if(any(patch%control_counts<1))return
        available=size(patch%controls,1)==3.and.size(patch%controls,2)==product(patch%control_counts)
    end function valid_controls

    logical function has_slice(s,counts) result(available)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: counts(3)
        available=.false.
        if (.not.allocated(s%live_slice)) return
        available=size(s%live_slice,1)==3.and.size(s%live_slice,2)==counts(1)*counts(2)
    end function has_slice

    pure integer function edge_total(kind) result(count)
        integer,intent(in) :: kind
        count=kind*2**(kind-1)
    end function edge_total

    pure integer function corner_sample(counts,kind,corner) result(node)
        integer,intent(in) :: counts(3),kind,corner
        integer :: axis,stride(3)

        stride=[1,counts(1),product(counts(:2))]
        node=1
        do axis=1,kind
            if(btest(corner-1,axis-1))node=node+(counts(axis)-1)*stride(axis)
        end do
    end function corner_sample

    pure subroutine edge_samples(counts,kind,edge,nodes,n,sides)
        integer,intent(in) :: counts(3),kind,edge
        integer,intent(out) :: nodes(:),n,sides(2)
        integer :: base,i,axis,fixed,other,stride(3)

        n=0; sides=0
        if(edge<1.or.edge>edge_total(kind))return
        stride=[1,counts(1),product(counts(:2))]
        base=1; axis=1
        if(kind==2)then
            fixed=(edge+1)/2
            axis=3-fixed
            sides=patch_edge_boundaries(kind,edge)
            if(mod(edge,2)==0)base=base+(counts(fixed)-1)*stride(fixed)
        else if(kind==3)then
            axis=(edge+3)/4
            fixed=mod(axis,3)+1; other=mod(axis+1,3)+1
            sides=patch_edge_boundaries(kind,edge)
            if(mod(sides(1),2)==0)base=base+(counts(fixed)-1)*stride(fixed)
            if(mod(sides(2),2)==0)base=base+(counts(other)-1)*stride(other)
        end if
        n=counts(axis)
        nodes(:n)=[(base+i*stride(axis),i=0,n-1)]
    end subroutine edge_samples

    function viewport_feature_text(s) result(label)
        type(studio_t), intent(in) :: s
        character(:), allocatable :: label
        character(5), parameter :: bounds(6)=[character(5)::'U min','U max','V min','V max','W min','W max']
        integer :: kind,number,sides(2),dummy(2),n,axis

        label=''
        if(s%selected_multipatch>0.or..not.s%live_scene%valid)return
        kind=s%viewport_feature(1); number=s%viewport_feature(4)
        if(kind/=s%live_scene%dimension.or.s%viewport_feature(2)/=s%live_scene%variant.or.number<1)return
        if(.not.active_patch_visible(s))return
        if(kind==1.and.s%live_scene%control_counts(1)==1)then
            label='Point'
            return
        end if
        select case(s%viewport_feature(3))
        case(pick_face)
            if(kind==2)label='Surface face'
            if(kind==3.and.number<=size(bounds))label='Face · '//bounds(number)
        case(pick_edge)
            if(kind==1)label='Curve edge'
            if(kind==2.and.number<=4)label='Edge · '//bounds(number)
            if(kind==3.and.number<=12)then
                call edge_samples([2,2,2],kind,number,dummy,n,sides)
                label='Edge · '//bounds(sides(1))//' / '//bounds(sides(2))
            end if
        case(pick_vertex)
            label='Vertex'
            if(kind==1)label='Endpoint'
            do axis=1,kind
                if(axis==1)then
                    label=label//' · '
                else
                    label=label//' / '
                end if
                label=label//bounds(2*axis-1+ibits(number-1,axis-1,1))
            end do
        case default
            return
        end select
    end function viewport_feature_text

    integer function pick_control_point(s,x,y,view) result(point)
        type(studio_t), intent(in) :: s
        real(dp), intent(in) :: x,y
        type(viewport_camera), intent(out), optional :: view
        type(viewport_camera) :: camera
        real(dp) :: width,height,projected(3),distance,nearest,nearest_depth
        integer :: i

        point=0
        if(s%mesh_only.or..not.s%layers(4).or..not.active_patch_visible(s).or..not.has_controls(s))return
        width=real(s%surface%width(),dp)
        height=real(s%surface%height(),dp)
        if(.not.in_drawable_view(width,height,x,y))return
        if(.not.prepare_camera(s,camera))return
        if(camera%scale<=0._dp)return
        if(present(view))view=camera

        nearest=control_hit_radius**2
        nearest_depth=-huge(1._dp)
        do i=1,size(s%live_scene%controls,2)
            projected=project_vertex(camera,s%live_scene%controls(:,i))
            if(.not.in_drawable_view(width,height,projected(1),projected(2)))cycle
            distance=(projected(1)-x)**2+(projected(2)-y)**2
            if(distance>nearest)cycle
            if(distance==nearest.and.projected(3)<=nearest_depth)cycle
            nearest=distance
            nearest_depth=projected(3)
            point=i
        end do
    end function pick_control_point

    subroutine pick_selected_control(s,x,y,kind,id,point)
        type(studio_t),intent(in) :: s
        real(dp),intent(in) :: x,y
        integer,intent(out) :: kind,id,point
        type(viewport_camera) :: camera
        type(projected_control_patch),allocatable :: patches(:)
        real(dp) :: width,height,distance,nearest,nearest_depth,projected(3)
        integer :: i,j
        logical :: focused

        kind=0; id=0; point=0
        if(s%mesh_only.or..not.s%layers(4).or.multipatch_control_scope(s))return
        width=real(s%surface%width(),dp)
        height=real(s%surface%height(),dp)
        if(.not.in_drawable_view(width,height,x,y))return
        if(.not.prepare_camera(s,camera))return
        call project_selected_controls(s,camera,[x-control_hit_radius,y-control_hit_radius, &
            x+control_hit_radius,y+control_hit_radius],patches)
        nearest=control_hit_radius**2
        nearest_depth=-huge(1._dp)
        focused=.false.
        do i=1,size(patches)
            associate(patch=>patches(i))
                do j=1,size(patch%points,2)
                    projected=patch%points(:,j)
                    if(.not.in_drawable_view(width,height,projected(1),projected(2)))cycle
                    distance=(projected(1)-x)**2+(projected(2)-y)**2
                    if(distance>nearest)cycle
                    if(distance==nearest)then
                        if(projected(3)<nearest_depth)cycle
                        if(projected(3)==nearest_depth.and.(focused.or..not.patch%focused))cycle
                    end if
                    nearest=distance
                    nearest_depth=projected(3)
                    focused=patch%focused
                    kind=patch%kind; id=patch%id; point=j
                end do
            end associate
        end do
    end subroutine pick_selected_control

    integer function pick_shared_control(s,x,y) result(point)
        type(studio_t),intent(in) :: s
        real(dp),intent(in) :: x,y
        type(viewport_camera) :: camera
        type(projected_control_patch),allocatable :: patches(:)
        real(dp) :: projected(3)
        logical,allocatable :: visible(:)
        real(dp) :: width,height,nearest,depth,distance
        integer :: i

        point=0
        if(s%mesh_only.or..not.s%layers(4).or..not.s%shared_controls%valid)return
        if(.not.multipatch_control_scope(s))return
        width=real(s%surface%width(),dp)
        height=real(s%surface%height(),dp)
        if(.not.in_drawable_view(width,height,x,y))return
        if(.not.prepare_camera(s,camera))return
        call project_selected_controls(s,camera,[x-control_hit_radius,y-control_hit_radius, &
            x+control_hit_radius,y+control_hit_radius],patches)
        call visible_shared_controls(s,patches,visible)
        nearest=control_hit_radius**2
        depth=-huge(1._dp)
        do i=1,size(visible)
            if(.not.visible(i))cycle
            projected=project_vertex(camera,s%shared_controls%points(:,i))
            if(.not.in_drawable_view(width,height,projected(1),projected(2)))cycle
            distance=sum((projected(:2)-[x,y])**2)
            if(distance>nearest.or.(distance==nearest.and.projected(3)<=depth))cycle
            nearest=distance
            depth=projected(3)
            point=i
        end do
    end function pick_shared_control

    logical function begin_live_control_drag(s,x,y,control) result(hit)
        type(studio_t), intent(inout) :: s
        real(dp), intent(in) :: x,y
        integer,intent(in),optional :: control
        type(viewport_camera) :: camera
        real(dp) :: origin(3)

        hit=.false.
        if(s%mesh_only)return
        s%control_dragging=.false.
        s%control_snapped=.false.
        call clear_control_dependency_preview(s)
        if(multipatch_control_scope(s))then
            if(.not.s%layers(4))return
            if(present(control))then
                s%selected_control=control
            else
                s%selected_control=pick_shared_control(s,x,y)
            end if
            if(.not.prepare_camera(s,camera))return
        else if(present(control))then
            if(.not.active_patch_visible(s).or..not.has_controls(s).or..not.s%layers(4))return
            if(control<1.or.control>size(s%live_scene%controls,2))return
            if(.not.prepare_camera(s,camera))return
            s%selected_control=control
        else
            s%selected_control=pick_control_point(s,x,y,camera)
        end if
        if(.not.viewport_control_position(s,s%selected_control,origin))return
        s%control_preview=origin

        ! The screen-plane axes are orthonormal rows of the view rotation.
        ! Their transpose maps a pixel displacement to a world displacement at fixed depth.
        s%control_drag_basis(:,1)=[camera%yaw_cos,-camera%yaw_sin,0._dp]/camera%scale
        s%control_drag_basis(:,2)=[-camera%pitch_sin*camera%yaw_sin, &
            -camera%pitch_sin*camera%yaw_cos,-camera%pitch_cos]/camera%scale
        s%last_x=x
        s%last_y=y
        s%control_dragging=.true.
        call begin_control_dependency_preview(s)
        hit=.true.
    end function begin_live_control_drag

    subroutine begin_control_dependency_preview(s)
        type(studio_t), intent(inout) :: s
        type(multipatch_controls_t) :: controls
        integer :: assembly,point

        associate(source=>s%control_dependency_source)
            if(multipatch_control_scope(s))then
                if(.not.s%shared_controls%smooth_joins)return
                source%assembly=s%selected_multipatch
                source%point=s%selected_control
                return
            end if
            assembly=joined_patch_assembly(s%document,s%entity_kind,s%live_scene%variant)
            if(assembly==0)return
            associate(group=>s%document%multipatches%groups(assembly)%definition)
                if(.not.allocated(group%joins))return
                if(.not.any(group%joins%continuity>0))return
            end associate
            call read_document_multipatch_controls(s%document,assembly,controls)
            if(.not.controls%valid)return
            point=assembly_control_id(controls,s%live_scene%variant,s%selected_control)
            if(point<1)return
            source%assembly=assembly
            source%point=point
            source%member_controls=controls
        end associate
    end subroutine begin_control_dependency_preview

    subroutine update_live_control_drag(s,x,y,temporary_snap,with_dependencies)
        type(studio_t), intent(inout) :: s
        real(dp), intent(in) :: x,y
        logical,intent(in),optional :: temporary_snap,with_dependencies
        real(dp) :: displacement(2),origin(3)
        logical :: snap,show_dependencies

        if(.not.s%control_dragging)return
        if(.not.viewport_control_position(s,s%selected_control,origin))return
        displacement=[x-s%last_x,y-s%last_y]
        s%control_preview=origin
        s%control_snapped=.false.
        if(norm2(displacement)>=control_drag_slop)then
            s%control_preview=s%control_preview+matmul(s%control_drag_basis,displacement)
            snap=s%snap_enabled
            if(present(temporary_snap))snap=snap.or.temporary_snap
            if(snap)then
                if(s%snap_closest)then
                    call snap_native_closest(s)
                else
                    call snap_native_corner(s,origin)
                end if
            end if
        end if
        show_dependencies=.true.
        if(present(with_dependencies))show_dependencies=with_dependencies
        if(show_dependencies)call update_control_dependency_preview(s,origin)
        call s%surface%update()
    end subroutine update_live_control_drag

    subroutine update_control_dependency_preview(s,origin)
        type(studio_t), intent(inout) :: s
        real(dp), intent(in) :: origin(3)
        type(multipatch_controls_t) :: edited

        associate(source=>s%control_dependency_source)
            if(source%assembly==0)return
            if(all(s%control_preview==origin))then
                source%attempted=.false.
                source%error=''
                s%control_dependency_preview=multipatch_control_preview_t()
                return
            end if
            if(source%attempted.and.all(source%position==s%control_preview))return
            source%attempted=.true.
            source%position=s%control_preview
            if(source%member_controls%valid)then
                edited=source%member_controls
                call set_multipatch_control_displacement(edited,source%point,s%control_preview)
                call preview_document_multipatch_controls(s%document,source%assembly,edited,s%control_dependency_preview, &
                    source%error,interactive=.true.,baseline=source%member_controls)
            else
                edited=s%shared_controls
                call set_multipatch_control_displacement(edited,source%point,s%control_preview,s%selected_shared_controls)
                call preview_document_multipatch_controls(s%document,source%assembly,edited,s%control_dependency_preview, &
                    source%error,interactive=.true.,baseline=s%shared_controls)
            end if
            if(len(source%error)>0)s%control_dependency_preview=multipatch_control_preview_t()
        end associate
    end subroutine update_control_dependency_preview

    subroutine snap_native_corner(s,start)
        type(studio_t), intent(inout) :: s
        real(dp),intent(in) :: start(3)
        type(viewport_camera) :: camera
        real(dp) :: nearest,depth,distance,width,height,candidate(3),screen(3),target(3),snapped(3)
        integer :: kind,variant,corner,node

        width=real(s%surface%width(),dp); height=real(s%surface%height(),dp)
        if(.not.prepare_camera(s,camera))return
        target=project_vertex(camera,s%control_preview)
        nearest=snap_pick_radius**2; depth=-huge(1._dp)
        do kind=1,size(s%document%scene%catalog)
            if(.not.allocated(s%document%scene%catalog(kind)%entries))cycle
            do variant=1,size(s%document%scene%catalog(kind)%entries)
                associate(entry=>s%document%scene%catalog(kind)%entries(variant))
                    if(.not.scene_patch_drawn(s,kind,variant))cycle
                    do corner=1,2**kind
                        node=corner_sample(entry%patch%sample_counts,kind,corner)
                        candidate=entry%patch%points(:,node)
                        if(all(candidate==start))then
                            if(kind==s%entity_kind.and.variant==s%live_scene%variant)cycle
                            if(multipatch_contains(s%document%multipatches,s%selected_multipatch,kind,variant))cycle
                        end if
                        screen=project_vertex(camera,candidate)
                        if(.not.in_drawable_view(width,height,screen(1),screen(2)))cycle
                        distance=sum((screen(:2)-target(:2))**2)
                        if(distance>nearest.or.(distance==nearest.and.screen(3)<depth))cycle
                        nearest=distance; depth=screen(3); snapped=candidate
                        s%control_snapped=.true.
                    end do
                end associate
            end do
        end do
        if(s%control_snapped)s%control_preview=snapped
    end subroutine snap_native_corner

    subroutine snap_native_closest(s)
        type(studio_t), intent(inout) :: s
        type(viewport_camera) :: camera
        type(projected_patch), allocatable :: objects(:)
        class(*), allocatable :: geometry
        real(dp) :: position(3),target(3),view_normal(3),query(3),candidate(3),screen(3),normalized(3)
        real(dp) :: depth,distance,nearest,nearest_depth
        integer :: i,kind,variant,face
        logical :: hit

        if (.not.prepare_camera(s,camera))return
        position=s%control_preview
        target=project_vertex(camera,position)
        if (.not.workspace_view_hit(s,target(1),target(2)))return
        view_normal=[-camera%pitch_cos*camera%yaw_sin,-camera%pitch_cos*camera%yaw_cos,camera%pitch_sin]
        call project_scene(s,camera,objects,clip=cursor_clip(target(1),target(2)),highlights=.false.)
        nearest=snap_pick_radius**2
        nearest_depth=-huge(1._dp)
        do i=1,size(objects)
            associate(object=>objects(i))
                if(object%selected.or.object%opacity==0)cycle
                kind=object%kind
                variant=object%variant
                if(multipatch_contains(s%document%multipatches,s%selected_multipatch,kind,variant))cycle
                normalized=0._dp
                face=0
                if(kind==1)then
                    call nearest_curve_parameter(object%points,target(:2),normalized(1),hit,distance,depth)
                else
                    call nearest_patch_face(object,target(:2),hit,distance,depth,face,normalized)
                end if
                if(.not.hit.or.distance>snap_pick_radius**2)cycle
                ! Locate the cursor at the rendered target's depth before exact native projection.
                query=position+(depth-target(3))*view_normal
                call copy_single_patch_geometry(kind,variant,geometry,catalog=s%document%patches)
                if(.not.allocated(geometry))cycle
                if(.not.closest_point_on_patch(geometry,face,s%document%scene%catalog(kind)%entries(variant)%patch%bounds, &
                    normalized,query,candidate))cycle
                screen=project_vertex(camera,candidate)
                if (.not.workspace_view_hit(s,screen(1),screen(2)))cycle
                distance=sum((screen(:2)-target(:2))**2)
                if(distance>snap_pick_radius**2)cycle
                if(s%control_snapped)then
                    ! Within a quarter square pixel, favor the visible front target.
                    if(distance>nearest+.25_dp)cycle
                    if(abs(distance-nearest)<=.25_dp.and.screen(3)<=nearest_depth)cycle
                end if
                nearest=distance
                nearest_depth=screen(3)
                s%control_preview=candidate
                s%control_snapped=.true.
            end associate
        end do
    end subroutine snap_native_closest

    subroutine pick_scene_patch(s,x,y,kind,variant,element,number)
        type(studio_t), intent(in) :: s
        real(dp), intent(in) :: x,y
        integer, intent(out) :: kind,variant
        integer,intent(out),optional :: element,number
        type(viewport_camera) :: camera
        type(projected_patch), allocatable :: objects(:)
        real(dp) :: width,height,distance,depth,best_depth,best_distance,parameter,depth_limit,mesh_depth,depth_bias
        integer :: mesh,part,cell,face,mesh_owner,mesh_kind
        integer :: i,side,best,best_side,which,classification,mode,first,last
        logical :: hit,wants_feature

        kind=0; variant=0
        if(present(element))element=0
        if(present(number))number=0
        if(s%mesh_only)return
        wants_feature=present(element).or.present(number)
        mode=pick_automatic
        if(wants_feature.and.s%feature_pick_mode>=pick_face.and.s%feature_pick_mode<=pick_vertex)mode=s%feature_pick_mode
        width=real(s%surface%width(),dp)
        height=real(s%surface%height(),dp)
        if(.not.in_drawable_view(width,height,x,y))return
        if(.not.prepare_camera(s,camera))return
        call project_scene(s,camera,objects,clip=cursor_clip(x,y),highlights=.false.)
        if(s%viewport_style==viewport_wireframe)then
            depth_limit=-huge(1._dp)
            best=0
            classification=pick_vertex
            if(mode==pick_automatic.or.mode==pick_vertex)then
                call choose_scene_feature(s,objects,1,size(objects),pick_vertex,[x,y],depth_limit,best,which)
            end if
            if(best==0.and.(mode==pick_automatic.or.mode==pick_edge))then
                classification=pick_edge
                call choose_scene_feature(s,objects,1,size(objects),pick_edge,[x,y],depth_limit,best,which)
            end if
            if(best==0)return
            kind=objects(best)%kind; variant=objects(best)%variant
            if(present(element))element=classification
            if(present(number))number=which
            return
        end if
        call pick_mesh_cell(s,x,y,mesh,part,cell,face,view_depth=mesh_depth,opaque_only=.true.)
        mesh_owner=0
        mesh_kind=0
        if(mesh>0)then
            mesh_owner=s%document%meshes(mesh)%parts(part)%variant
            mesh_kind=s%document%meshes(mesh)%kind
        end if
        depth_bias=.85_dp/camera%scale
        best_depth=-huge(1._dp)
        best_distance=huge(1._dp)
        best=0; best_side=0
        do i=1,size(objects)
            associate(object=>objects(i))
                if(.not.pickable_object(s,object))cycle
                if(object%kind==1)then
                    if(mode==pick_face)cycle
                    call nearest_curve_parameter(object%points,[x,y],parameter,hit,distance,depth)
                    if(distance<=25._dp)distance=0._dp
                    if(distance>100._dp)hit=.false.
                    side=1
                else
                    call nearest_patch_face(object,[x,y],hit,distance,depth,side)
                end if
                if(.not.hit.or.distance>best_distance)cycle
                if(depth<mesh_depth-depth_bias)then
                    if(object%kind/=mesh_kind.or.object%variant/=mesh_owner)cycle
                end if
                if(distance==best_distance.and.depth<best_depth)cycle
                if(distance==best_distance.and.depth==best_depth.and..not.object%selected.and.best/=0)cycle
                best=i; best_side=side
                best_depth=depth; best_distance=distance
            end associate
        end do
        if(.not.wants_feature)then
            if(best==0)return
            kind=objects(best)%kind; variant=objects(best)%variant
            return
        end if

        depth_limit=-huge(1._dp)
        if(best_distance==0._dp)depth_limit=best_depth-6._dp/camera%scale
        if(mesh>0)then
            if(best==0)then
                depth_limit=max(depth_limit,mesh_depth-depth_bias)
            else if(objects(best)%kind/=mesh_kind.or.objects(best)%variant/=mesh_owner)then
                depth_limit=max(depth_limit,mesh_depth-depth_bias)
            end if
        end if
        which=0
        first=1; last=size(objects)
        if(mode==pick_automatic.and.best>0)then
            first=best; last=best
        end if
        select case(mode)
        case(pick_automatic)
            call choose_scene_feature(s,objects,first,last,pick_vertex,[x,y],depth_limit,i,which)
            classification=pick_vertex
            if(i==0)then
                call choose_scene_feature(s,objects,first,last,pick_edge,[x,y],depth_limit,i,which)
                classification=pick_edge
            end if
            if(i==0.and.best>0)then
                i=best; which=best_side; classification=merge(pick_edge,pick_face,objects(best)%kind==1)
                if(which==0)classification=0
            end if
        case(pick_face)
            i=best; which=best_side; classification=pick_face
            if(which==0)i=0
        case(pick_edge,pick_vertex)
            call choose_scene_feature(s,objects,first,last,mode,[x,y],depth_limit,i,which)
            classification=mode
        case default
            return
        end select
        if(i==0)return
        kind=objects(i)%kind; variant=objects(i)%variant
        if(present(element))element=classification
        if(present(number))number=which
    end subroutine pick_scene_patch

    real(dp) function opaque_patch_depth(s,x,y) result(depth)
        type(studio_t),intent(in) :: s
        real(dp),intent(in) :: x,y
        type(viewport_camera) :: camera
        type(projected_patch),allocatable :: objects(:)
        real(dp) :: candidate,distance
        logical :: hit
        integer :: i,side

        depth=-huge(1._dp)
        if(s%mesh_only.or.s%viewport_style==viewport_wireframe)return
        if(.not.prepare_camera(s,camera))return
        call project_scene(s,camera,objects,clip=cursor_clip(x,y),highlights=.false.)
        do i=1,size(objects)
            associate(object=>objects(i))
                if(object%kind==1.or.object%opacity<255)cycle
                if(object%selected.and..not.s%layers(1))cycle
                if(mesh_replaces_patch_fill(s,object%kind,object%variant))cycle
                call nearest_patch_face(object,[x,y],hit,distance,candidate,side)
                if(hit.and.distance==0._dp)depth=max(depth,candidate)
            end associate
        end do
        depth=depth-.85_dp/camera%scale
    end function opaque_patch_depth

    pure function cursor_clip(x,y) result(clip)
        real(dp),intent(in) :: x,y
        real(dp) :: clip(4)

        clip=[x-snap_pick_radius,y-snap_pick_radius,x+snap_pick_radius,y+snap_pick_radius]
    end function cursor_clip

    logical function pickable_object(s,object) result(pickable)
        type(studio_t), intent(in) :: s
        type(projected_patch), intent(in) :: object
        pickable=.not.object%selected.or.any(s%layers([1,2,parameter_lines_layer]))
        if(s%mesh_section_axis>0.and.object%kind==3)then
            pickable=pickable.and..not.mesh_replaces_patch_fill(s,object%kind,object%variant)
        end if
    end function pickable_object

    subroutine choose_scene_feature(s,objects,first,last,feature,click,depth_limit,best,number)
        type(studio_t),intent(in) :: s
        type(projected_patch),intent(in) :: objects(:)
        integer,intent(in) :: first,last,feature
        real(dp),intent(in) :: click(2),depth_limit
        integer,intent(out) :: best,number
        integer :: i,candidate
        real(dp) :: distance,depth,best_distance,best_depth

        best=0; number=0
        best_distance=huge(1._dp); best_depth=-huge(1._dp)
        do i=first,last
            if(.not.pickable_object(s,objects(i)))cycle
            if(feature==pick_vertex)then
                call nearest_patch_vertex(objects(i),click,depth_limit,candidate,distance,depth)
            else
                call nearest_patch_edge(objects(i),click,depth_limit,candidate,distance,depth)
            end if
            if(candidate==0.or.distance>best_distance)cycle
            if(distance==best_distance.and.depth<best_depth)cycle
            best=i; number=candidate; best_distance=distance; best_depth=depth
        end do
    end subroutine choose_scene_feature

    pure subroutine nearest_patch_vertex(object,click,depth_limit,vertex,distance,depth)
        type(projected_patch),intent(in) :: object
        real(dp),intent(in) :: click(2),depth_limit
        integer,intent(out) :: vertex
        real(dp),intent(out) :: distance,depth
        real(dp) :: candidate(3),d
        integer :: i,node

        vertex=0; distance=vertex_pick_radius**2; depth=-huge(1._dp)
        do i=1,2**object%kind
            node=corner_sample(object%counts,object%kind,i)
            candidate=object%points(:,node)
            if(candidate(3)<depth_limit)cycle
            d=sum((candidate(:2)-click)**2)
            if(d>distance.or.(d==distance.and.candidate(3)<depth))cycle
            vertex=i; distance=d; depth=candidate(3)
        end do
    end subroutine nearest_patch_vertex

    subroutine nearest_patch_edge(object,click,depth_limit,edge,distance,depth)
        type(projected_patch),intent(in) :: object
        real(dp),intent(in) :: click(2),depth_limit
        integer,intent(out) :: edge
        real(dp),intent(out) :: distance,depth
        real(dp) :: a(3),b(3),delta(2),length,fraction,candidate(3),d
        real(dp),allocatable :: points(:,:)
        integer :: i,j,sides(2)

        edge=0; distance=edge_pick_radius**2; depth=-huge(1._dp)
        do i=1,edge_total(object%kind)
            points=projected_edge(object,i,sides)
            do j=1,size(points,2)-1
                a=points(:,j); b=points(:,j+1)
                delta=b(:2)-a(:2)
                length=sum(delta**2)
                if(length<1.e-8_dp)cycle
                fraction=max(0._dp,min(1._dp,dot_product(click-a(:2),delta)/length))
                candidate=a+fraction*(b-a)
                if(candidate(3)<depth_limit)cycle
                d=sum((candidate(:2)-click)**2)
                if(d>distance.or.(d==distance.and.candidate(3)<depth))cycle
                edge=i; distance=d; depth=candidate(3)
            end do
        end do
    end subroutine nearest_patch_edge

    subroutine nearest_patch_face(object,click,hit,distance,depth,face,normalized)
        type(projected_patch), intent(in) :: object
        real(dp), intent(in) :: click(2)
        logical, intent(out) :: hit
        real(dp), intent(out) :: distance,depth
        integer,intent(out) :: face
        real(dp),intent(out),optional :: normalized(3)
        integer :: side,axes(2),counts(2)
        real(dp) :: padding

        distance=huge(1._dp); depth=-huge(1._dp)
        face=0; hit=.false.
        padding=0._dp
        if(present(normalized))then
            normalized=0._dp
            padding=snap_pick_radius
        end if
        if(object%kind==3)then
            do side=1,6
                call projected_face_axes(object,side,axes,counts)
                call test_face(object%boundary(side)%points,side,axes,counts)
            end do
        else
            call projected_face_axes(object,1,axes,counts)
            call test_face(object%points,1,axes,counts)
        end if
        if(hit.or.present(normalized))return
        if(object%kind==3)then
            do side=1,6
                call test_vertices(object%boundary(side)%points)
            end do
        else
            call test_vertices(object%points)
        end if
    contains
        subroutine test_face(points,side,axes,counts)
            real(dp),intent(in) :: points(:,:)
            integer,intent(in) :: side,axes(2),counts(2)
            integer,parameter :: triangles(3,2)=reshape([1,2,3,1,3,4],[3,2])
            integer :: i,j,base,vertices(4),rows(3),trial,stride,axis,strip_start,strip_end
            real(dp) :: weight(3),d,screen_distance,corners(2,4)

            if(any(counts<2))return
            if(any(click<minval(points(:2,:),dim=2)-padding).or.any(click>maxval(points(:2,:),dim=2)+padding))return
            stride=counts(1)
            do j=1,counts(2)-1
                strip_start=1+(j-1)*stride
                strip_end=(j+1)*stride
                if(any(click<minval(points(:2,strip_start:strip_end),dim=2)-padding))cycle
                if(any(click>maxval(points(:2,strip_start:strip_end),dim=2)+padding))cycle
                do i=1,stride-1
                    base=strip_start+i-1
                    vertices=base+[0,1,stride+1,stride]
                    corners=points(:2,vertices)
                    if(any(click<minval(corners,dim=2)-padding).or.any(click>maxval(corners,dim=2)+padding))cycle
                    do trial=1,2
                        rows=vertices(triangles(:,trial))
                        screen_distance=0._dp
                        if(present(normalized))then
                            call screen_triangle_nearest(corners(:,triangles(:,trial)),click, &
                                min(distance,padding**2),weight,screen_distance)
                            if(screen_distance>padding**2.or.screen_distance>distance)cycle
                        else
                            if(.not.screen_triangle(corners(:,triangles(:,trial)),click,weight))cycle
                        end if
                        d=dot_product(weight,points(3,rows))
                        if(screen_distance==distance.and.d<depth)cycle
                        distance=screen_distance
                        depth=d; face=side
                        hit=.true.
                        if(present(normalized))then
                            normalized(axes(1))=dot_product(weight,real(mod(rows-1,stride),dp))/real(stride-1,dp)
                            normalized(axes(2))=dot_product(weight,real((rows-1)/stride,dp))/real(counts(2)-1,dp)
                            if(object%kind==3)then
                                axis=(side+1)/2
                                normalized(axis)=real(mod(side+1,2),dp)
                            end if
                        end if
                    end do
                end do
            end do
        end subroutine test_face

        subroutine test_vertices(points)
            real(dp),intent(in) :: points(:,:)
            integer :: i
            real(dp) :: d

            if(any(click<minval(points(:2,:),dim=2)-8._dp).or.any(click>maxval(points(:2,:),dim=2)+8._dp))return
            do i=1,size(points,2)
                d=sum((points(:2,i)-click)**2)
                if(d>64._dp.or.d>=distance)cycle
                distance=d; depth=points(3,i)
                hit=.true.
            end do
        end subroutine test_vertices
    end subroutine nearest_patch_face

    subroutine pick_live_viewport(s,x,y)
        type(studio_t), intent(inout) :: s
        real(dp), intent(in) :: x,y
        type(viewport_camera) :: camera
        real(dp), allocatable :: projected(:,:)
        real(dp) :: width,height,uv(2)
        integer :: counts(3),dimension
        logical :: hit

        width=real(s%surface%width(),dp)
        height=real(s%surface%height(),dp)
        if (.not.in_drawable_view(width,height,x,y)) return
        if (.not.active_patch_visible(s)) return
        if (.not.any(s%layers([1,2,parameter_lines_layer]))) return
        if (.not.prepare_camera(s,camera)) return
        counts=s%live_scene%sample_counts
        dimension=s%live_scene%dimension
        if (dimension==3) then
            if (.not.has_slice(s,counts)) return
            if (.not.show_volume_slice(s)) return
            projected=project_points(s,camera,s%live_slice)
        else if(display_grid_current(s,dimension,s%live_scene%variant))then
            associate(entry=>s%display_cache%ranks(dimension)%entries(s%live_scene%variant))
                counts=entry%counts
                projected=project_points(s,camera,entry%grids(1)%points)
            end associate
        else
            projected=project_points(s,camera,s%live_scene%points)
        end if

        if (dimension==1) then
            call nearest_curve_parameter(projected,[x,y],uv(1),hit)
        else
            call pick_sample_grid(projected,counts(1),counts(2),[x,y],uv,hit)
        end if
        if (.not.hit) return
        call s%probe_u%setValue(min(1._dp,max(0._dp,uv(1))))
        if (dimension>1) call s%probe_v%setValue(min(1._dp,max(0._dp,uv(2))))
        ! Volume picking is confined to the evaluated W slice and leaves W unchanged.
    end subroutine pick_live_viewport

end module studio_viewport_scene
