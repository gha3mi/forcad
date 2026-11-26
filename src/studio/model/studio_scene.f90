! Cache native patch samples for the 3D scene without changing the active patch.
module studio_scene
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: iso_fortran_env, only: int64
    use studio_single_patch, only: single_patch_scene, single_patch_count, single_patch_revision, read_single_patch_scene, &
        single_patch_available, next_patch_revision, usable_undo_stamp, studio_patch_catalog,patch_history_project
    use studio_history_storage, only: history_blob,release_history_blob
    use studio_property_stream, only: studio_property_values,save_property_history,load_property_history
    implicit none(type, external)
    private
    integer, parameter :: dp=c_double
    real(dp), parameter, public :: scene_yaw=.65_dp, scene_pitch=.56_dp

    type :: scene_entry
        logical :: visible=.false.
        integer(int64) :: revision=-1_int64
        type(single_patch_scene) :: patch
        integer :: rgb(3)=-1,opacity=255
    end type scene_entry

    type :: scene_catalog
        type(scene_entry), allocatable :: entries(:)
    end type scene_catalog

    type, public :: studio_material_replay
        integer(int64) :: stamp=0
        type(history_blob) :: payload
    end type studio_material_replay

    type, public :: studio_scene_t
        type(scene_catalog) :: catalog(3)
        type(studio_material_replay),allocatable,private :: materials(:)
        logical :: framed=.false.
        real(dp) :: center(3)=0._dp, span(2)=1._dp, radius=1._dp, floor=0._dp
    end type studio_scene_t

    public :: refresh_studio_scene, studio_scene_is_visible, studio_scene_set_visible, frame_studio_scene
    public :: studio_scene_count, studio_scene_control_position
    public :: studio_scene_set_visibility, project_patch_exists, project_geometry_targets
    public :: set_studio_material,get_studio_material
    public :: commit_studio_material,studio_material_undo_stamp,undo_studio_material,redo_studio_material,clear_material_history
    public :: snapshot_scene_properties,restore_scene_properties,move_studio_scene

contains

    subroutine move_studio_scene(source,destination)
        type(studio_scene_t), intent(inout) :: source
        type(studio_scene_t), intent(out) :: destination
        integer :: kind

        do kind=1,size(source%catalog)
            call move_alloc(source%catalog(kind)%entries,destination%catalog(kind)%entries)
        end do
        call move_alloc(source%materials,destination%materials)
        destination%framed=source%framed
        destination%center=source%center
        destination%span=source%span
        destination%radius=source%radius
        destination%floor=source%floor
        source=studio_scene_t()
    end subroutine move_studio_scene

    logical function project_patch_exists(include_templates,kind,variant,patches) result(present_in_project)
        logical, intent(in) :: include_templates
        type(studio_patch_catalog), intent(in), optional :: patches
        integer, intent(in) :: kind,variant

        present_in_project=(include_templates.or.variant>2).and.single_patch_available(kind,variant,patches)
    end function project_patch_exists

    function project_geometry_targets(include_templates,patches) result(targets)
        logical, intent(in) :: include_templates
        type(studio_patch_catalog), intent(in), optional :: patches
        integer, allocatable :: targets(:,:), catalog(:,:)
        integer :: kind,variant,n

        allocate(catalog(2,sum([(single_patch_count(kind,patches),kind=1,3)])))
        n=0
        do kind=1,3
            do variant=1,single_patch_count(kind,patches)
                if(.not.project_patch_exists(include_templates,kind,variant,patches))cycle
                n=n+1
                catalog(:,n)=[kind,variant]
            end do
        end do
        targets=catalog(:,:n)
    end function project_geometry_targets


    subroutine set_studio_material(scene,kind,variant,rgb,opacity)
        type(studio_scene_t), intent(inout) :: scene
        integer, intent(in) :: kind,variant
        integer, intent(in), optional :: rgb(3),opacity
        if(kind<1.or.kind>size(scene%catalog))return
        if(.not.allocated(scene%catalog(kind)%entries))return
        if(variant<1.or.variant>size(scene%catalog(kind)%entries))return
        associate(entry=>scene%catalog(kind)%entries(variant))
            if(present(rgb))entry%rgb=max(0,min(255,rgb))
            if(present(opacity))entry%opacity=max(0,min(255,opacity))
        end associate
    end subroutine set_studio_material

    subroutine get_studio_material(scene,kind,variant,rgb,opacity)
        type(studio_scene_t), intent(in) :: scene
        integer, intent(in) :: kind,variant
        integer, intent(out) :: rgb(3),opacity
        integer, parameter :: defaults(3,3)=reshape([242,192,114,77,165,196,119,108,166],[3,3])
        rgb=defaults(:,max(1,min(3,kind)))
        opacity=255
        if(kind<1.or.kind>size(scene%catalog))return
        if(.not.allocated(scene%catalog(kind)%entries))return
        if(variant<1.or.variant>size(scene%catalog(kind)%entries))return
        associate(entry=>scene%catalog(kind)%entries(variant))
            if(all(entry%rgb>=0))rgb=entry%rgb
            opacity=entry%opacity
        end associate
    end subroutine get_studio_material

    logical function commit_studio_material(scene,kind,variants,rgb,opacity,catalog,saved) result(changed)
        type(studio_scene_t),intent(inout) :: scene
        integer,intent(in) :: kind,variants(:)
        integer,intent(in),optional :: rgb(3),opacity
        type(studio_patch_catalog), intent(inout), optional :: catalog
        logical, intent(out), optional :: saved
        type(studio_material_replay) :: edit
        type(studio_property_values) :: values
        integer :: i,current_rgb(3),current_opacity,targets(2,size(variants))

        changed=.false.
        if(present(saved))saved=.true.
        if(.not.present(rgb).and..not.present(opacity))return
        if(kind<1.or.kind>size(scene%catalog).or.size(variants)==0)return
        if(.not.allocated(scene%catalog(kind)%entries))return
        if(any(variants<1).or.any(variants>size(scene%catalog(kind)%entries)))return
        do i=1,size(variants)
            if(.not.single_patch_available(kind,variants(i),catalog))return
            call get_studio_material(scene,kind,variants(i),current_rgb,current_opacity)
            if(present(rgb))changed=changed.or.any(current_rgb/=rgb)
            if(present(opacity))changed=changed.or.current_opacity/=opacity
        end do
        if(.not.changed)return

        targets(1,:)=kind
        targets(2,:)=variants
        values%patches=snapshot_scene_properties(scene,targets)
        changed=save_property_history(values,edit%payload,patch_history_project(catalog))
        if(present(saved))saved=changed
        if(.not.changed)return
        edit%stamp=next_patch_revision()
        call push_material_edit(scene,edit)
        do i=1,size(variants)
            call set_studio_material(scene,kind,variants(i),rgb,opacity)
        end do
    end function commit_studio_material

    integer(int64) function studio_material_undo_stamp(scene,catalog) result(stamp)
        type(studio_scene_t),intent(in) :: scene
        type(studio_patch_catalog), intent(in), optional :: catalog
        integer :: n

        stamp=0
        if(.not.allocated(scene%materials))return
        n=size(scene%materials)
        if(n>0)stamp=usable_undo_stamp(scene%materials(n)%stamp,catalog=catalog)
    end function studio_material_undo_stamp

    subroutine push_material_edit(scene,edit)
        type(studio_scene_t), intent(inout) :: scene
        type(studio_material_replay), intent(in) :: edit

        if(.not.allocated(scene%materials))allocate(scene%materials(0))
        scene%materials=[scene%materials,edit]
    end subroutine push_material_edit

    subroutine undo_studio_material(scene,catalog,replay)
        type(studio_scene_t), intent(inout) :: scene
        type(studio_patch_catalog), intent(in), optional :: catalog
        type(studio_material_replay), intent(out), optional :: replay
        type(studio_material_replay) :: edit
        type(studio_property_values) :: old,current
        integer :: n

        if(studio_material_undo_stamp(scene,catalog)==0)return
        n=size(scene%materials)
        edit=scene%materials(n)
        if(.not.load_property_history(edit%payload,old))return
        if(present(replay))then
            current%patches=snapshot_scene_properties(scene,old%patches(:2,:))
            if(.not.save_property_history(current,replay%payload,patch_history_project(catalog)))return
            replay%stamp=edit%stamp
        end if
        call restore_scene_properties(scene,old%patches,catalog)
        call release_history_blob(edit%payload)
        scene%materials=scene%materials(:n-1)
    end subroutine undo_studio_material

    logical function redo_studio_material(scene,replay,catalog) result(changed)
        type(studio_scene_t), intent(inout) :: scene
        type(studio_material_replay), intent(inout) :: replay
        type(studio_patch_catalog), intent(in), optional :: catalog
        type(studio_material_replay) :: edit
        type(studio_property_values) :: future,current

        changed=.false.
        if(.not.load_property_history(replay%payload,future))return
        current%patches=snapshot_scene_properties(scene,future%patches(:2,:))
        if(.not.save_property_history(current,edit%payload,patch_history_project(catalog)))return
        edit%stamp=replay%stamp
        call restore_scene_properties(scene,future%patches,catalog)
        call push_material_edit(scene,edit)
        call release_history_blob(replay%payload)
        changed=.true.
    end function redo_studio_material

    subroutine clear_material_history(scene)
        type(studio_scene_t), intent(inout) :: scene
        integer :: i

        if(.not.allocated(scene%materials))return
        do i=1,size(scene%materials)
            call release_history_blob(scene%materials(i)%payload)
        end do
        deallocate(scene%materials)
    end subroutine clear_material_history

    function snapshot_scene_properties(scene,targets) result(values)
        type(studio_scene_t), intent(in) :: scene
        integer, intent(in) :: targets(:,:)
        integer :: values(7,size(targets,2)),i,kind,id

        do i=1,size(targets,2)
            kind=targets(1,i)
            id=targets(2,i)
            values(:,i)=[kind,id,merge(1,0,id>=2),-1,-1,-1,255]
            if(.not.allocated(scene%catalog(kind)%entries))cycle
            if(id>size(scene%catalog(kind)%entries))cycle
            associate(entry=>scene%catalog(kind)%entries(id))
                values(3:,i)=[merge(1,0,entry%visible),entry%rgb,entry%opacity]
            end associate
        end do
    end function snapshot_scene_properties

    subroutine restore_scene_properties(scene,values,catalog)
        type(studio_scene_t), intent(inout) :: scene
        integer, intent(in) :: values(:,:)
        type(studio_patch_catalog), intent(in), optional :: catalog
        integer :: i,kind,id
        logical :: visibility_changed,visible

        do i=1,size(values,2)
            kind=values(1,i)
            id=values(2,i)
            if(.not.single_patch_available(kind,id,catalog))cycle
            if(allocated(scene%catalog(kind)%entries))then
                if(id<=size(scene%catalog(kind)%entries))cycle
            end if
            call refresh_studio_scene(scene,catalog)
            exit
        end do

        visibility_changed=.false.
        do i=1,size(values,2)
            kind=values(1,i)
            id=values(2,i)
            if(.not.single_patch_available(kind,id,catalog))cycle
            visible=values(3,i)==1
            associate(entry=>scene%catalog(kind)%entries(id))
                visibility_changed=visibility_changed.or.(entry%visible.neqv.visible)
                entry%visible=visible
                entry%rgb=values(4:6,i)
                entry%opacity=values(7,i)
            end associate
        end do
        if(visibility_changed)call refresh_studio_scene(scene,catalog)
    end subroutine restore_scene_properties

    subroutine refresh_studio_scene(scene,catalog)
        type(studio_scene_t), intent(inout) :: scene
        type(studio_patch_catalog), intent(in), optional :: catalog
        type(scene_entry), allocatable :: resized(:)
        integer :: kind,variant,n,previous_count,retained
        integer(int64) :: revision

        do kind=1,size(scene%catalog)
            n=single_patch_count(kind,catalog)
            associate(group=>scene%catalog(kind))
                previous_count=0
                if(allocated(group%entries))previous_count=size(group%entries)
                if(n/=previous_count)then
                    allocate(resized(n))
                    retained=min(n,previous_count)
                    do variant=1,retained
                        call move_scene_entry(group%entries(variant),resized(variant))
                    end do
                    do variant=retained+1,n
                        resized(variant)%visible=variant>=2
                    end do
                    call move_alloc(resized,group%entries)
                end if
                do variant=1,n
                    if(.not.single_patch_available(kind,variant,catalog))then
                        group%entries(variant)%patch=single_patch_scene()
                        group%entries(variant)%revision=-1_int64
                        cycle
                    end if
                    if(.not.group%entries(variant)%visible)cycle
                    revision=single_patch_revision(kind,variant,catalog)
                    if(group%entries(variant)%revision==revision)cycle
                    call read_single_patch_scene(kind,group%entries(variant)%patch,variant,catalog)
                    group%entries(variant)%revision=revision
                end do
            end associate
        end do
        if(studio_scene_count(scene)==0)scene%framed=.false.
        if(.not.scene%framed)call frame_studio_scene(scene)
    end subroutine refresh_studio_scene

    subroutine move_scene_entry(source,destination)
        type(scene_entry), intent(inout) :: source
        type(scene_entry), intent(out) :: destination
        integer :: axis

        destination%visible=source%visible
        destination%revision=source%revision
        destination%rgb=source%rgb
        destination%opacity=source%opacity
        associate(from=>source%patch,to=>destination%patch)
            to%valid=from%valid
            to%dimension=from%dimension
            to%variant=from%variant
            to%degree=from%degree
            to%sample_counts=from%sample_counts
            to%control_counts=from%control_counts
            to%bounds=from%bounds
            call move_alloc(from%points,to%points)
            call move_alloc(from%controls,to%controls)
            call move_alloc(from%diagnostic,to%diagnostic)
            do axis=1,size(from%axes)
                call move_alloc(from%axes(axis)%knots,to%axes(axis)%knots)
            end do
        end associate
    end subroutine move_scene_entry

    pure logical function studio_scene_is_visible(scene,kind,variant,catalog) result(visible)
        type(studio_scene_t), intent(in) :: scene
        integer, intent(in) :: kind,variant
        type(studio_patch_catalog), intent(in), optional :: catalog

        visible=.false.
        if(kind<1.or.kind>size(scene%catalog))return
        if(.not.allocated(scene%catalog(kind)%entries))return
        if(variant<1.or.variant>size(scene%catalog(kind)%entries))return
        if(.not.single_patch_available(kind,variant,catalog))return
        visible=scene%catalog(kind)%entries(variant)%visible
    end function studio_scene_is_visible

    subroutine studio_scene_set_visible(scene,kind,variant,visible,catalog)
        type(studio_scene_t), intent(inout) :: scene
        integer, intent(in) :: kind,variant
        logical, intent(in) :: visible
        type(studio_patch_catalog), intent(in), optional :: catalog

        if(.not.single_patch_available(kind,variant,catalog))return
        if(.not.allocated(scene%catalog(kind)%entries))call refresh_studio_scene(scene,catalog)
        if(variant>size(scene%catalog(kind)%entries))call refresh_studio_scene(scene,catalog)
        if(scene%catalog(kind)%entries(variant)%visible.eqv.visible)return
        scene%catalog(kind)%entries(variant)%visible=visible
        call refresh_studio_scene(scene,catalog)
    end subroutine studio_scene_set_visible

    subroutine studio_scene_set_visibility(scene,visible,kind,variants,catalog)
        type(studio_scene_t),intent(inout)::scene
        logical,intent(in)::visible
        integer,intent(in),optional::kind,variants(:)
        type(studio_patch_catalog), intent(in), optional :: catalog
        integer::rank,variant

        call refresh_studio_scene(scene,catalog)
        do rank=1,size(scene%catalog)
            do variant=1,size(scene%catalog(rank)%entries)
                if(.not.single_patch_available(rank,variant,catalog))cycle
                scene%catalog(rank)%entries(variant)%visible=visible
                if(.not.present(kind).or..not.present(variants))cycle
                if(rank==kind.and.any(variants==variant))scene%catalog(rank)%entries(variant)%visible=.true.
            end do
        end do
        call refresh_studio_scene(scene,catalog)
    end subroutine studio_scene_set_visibility

    pure integer function studio_scene_count(scene) result(n)
        type(studio_scene_t), intent(in) :: scene
        integer :: kind

        n=0
        do kind=1,size(scene%catalog)
            if(.not.allocated(scene%catalog(kind)%entries))cycle
            n=n+count(scene%catalog(kind)%entries%visible.and.scene%catalog(kind)%entries%patch%valid)
        end do
    end function studio_scene_count

    logical function studio_scene_control_position(scene,kind,variant,point,position,catalog) result(available)
        type(studio_scene_t), intent(in) :: scene
        integer, intent(in) :: kind,variant,point
        real(dp), intent(out) :: position(3)
        type(studio_patch_catalog), intent(in), optional :: catalog

        available=.false.
        position=0._dp
        if(.not.studio_scene_is_visible(scene,kind,variant,catalog))return
        associate(patch=>scene%catalog(kind)%entries(variant)%patch)
            if(.not.patch%valid.or..not.allocated(patch%controls))return
            if(point<1.or.point>size(patch%controls,2))return
            position=patch%controls(:,point)
        end associate
        available=.true.
    end function studio_scene_control_position

    subroutine frame_studio_scene(scene,yaw,pitch,kind,variant,members,include_controls,control_target)
        type(studio_scene_t), intent(inout) :: scene
        real(dp), intent(in), optional :: yaw,pitch
        integer, intent(in), optional :: kind,variant,members(:),control_target(2)
        logical, intent(in), optional :: include_controls
        real(dp) :: lower(3),upper(3),lower_view(3),upper_view(3),basis(3,3),u,v,half(3)
        integer :: rank,id
        logical :: selected_only,found,controls,frame_controls

        if(present(variant).and.present(members))return
        if(present(kind).neqv.(present(variant).or.present(members)))return
        selected_only=present(kind)
        controls=.false.
        if(present(include_controls))controls=include_controls
        if(.not.selected_only.and.studio_scene_count(scene)==0)then
            scene%framed=.false.
            return
        end if
        u=scene_yaw; v=scene_pitch
        if(present(yaw))u=scene_yaw+yaw
        if(present(pitch))v=scene_pitch+pitch
        basis(1,:)=[cos(u),-sin(u),0._dp]
        basis(2,:)=[-sin(v)*sin(u),-sin(v)*cos(u),-cos(v)]
        basis(3,:)=[-cos(v)*sin(u),-cos(v)*cos(u),sin(v)]
        lower=huge(1._dp); lower_view=huge(1._dp)
        upper=-huge(1._dp); upper_view=-huge(1._dp)
        found=.false.
        do rank=1,size(scene%catalog)
            if(selected_only)then
                if(rank/=kind)cycle
            end if
            if(.not.allocated(scene%catalog(rank)%entries))cycle
            do id=1,size(scene%catalog(rank)%entries)
                if(present(variant))then
                    if(id/=variant)cycle
                else if(present(members))then
                    if(.not.any(members==id))cycle
                end if
                associate(entry=>scene%catalog(rank)%entries(id))
                    if(.not.entry%visible.or..not.entry%patch%valid)cycle
                    found=.true.
                    call include_points(entry%patch%points,basis,lower,upper,lower_view,upper_view)
                    frame_controls=controls
                    if(present(control_target))frame_controls=controls.and.all(control_target==[rank,id])
                    if(frame_controls)call include_points(entry%patch%controls,basis,lower,upper,lower_view,upper_view)
                end associate
            end do
        end do
        if(.not.found)return
        half=.5_dp*(upper-lower)
        scene%center=matmul(transpose(basis),lower_view+.5_dp*(upper_view-lower_view))
        scene%radius=max(norm2(half),sqrt(tiny(1._dp)))
        scene%span=max(upper_view(:2)-lower_view(:2),.12_dp*scene%radius)
        if(.not.selected_only)scene%floor=lower(3)-.06_dp*scene%radius
        scene%framed=.true.
    end subroutine frame_studio_scene

    subroutine include_points(points,basis,lower,upper,lower_view,upper_view)
        real(dp), intent(in) :: points(:,:),basis(3,3)
        real(dp), intent(inout) :: lower(3),upper(3),lower_view(3),upper_view(3)
        real(dp) :: projected(3,size(points,2))

        if(size(points,2)==0)return
        lower=min(lower,minval(points,dim=2))
        upper=max(upper,maxval(points,dim=2))
        projected=matmul(basis,points)
        lower_view=min(lower_view,minval(projected,dim=2))
        upper_view=max(upper_view,maxval(projected,dim=2))
    end subroutine include_points

end module studio_scene
