module studio_state
    use studio_document_tabs, only: studio_document_tabs_t
    use studio_display_cache, only: viewport_display_cache
    use studio_commands, only: viewport_layer_count, default_viewport_layers, &
        viewport_shaded,viewport_edges,viewport_wireframe,viewport_fast_shaded
    use qt, only: QAction, QDoubleSpinBox, QLabel, QMainWindow, QPushButton, QSplitter, &
        QStackedWidget, QTableWidget, QTimer, QToolButton, QWidget, Qt_NoButton, QtfWidget, qt_string
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: iso_fortran_env, only: int64
    use studio_single_patch, only: single_patch_info, single_patch_scene, single_patch_plot, initialize_studio_patch_catalog
    use studio_scene, only: frame_studio_scene, studio_scene_count, studio_scene_is_visible
    use studio_mesh, only: studio_mesh_t
    use studio_multipatch_controls, only: multipatch_controls_t, multipatch_control_preview_t
    use studio_selection, only: selection_set_t, selection_target_t, selection_patch
    use studio_status_bar, only: studio_status_t
    use studio_inspector, only: studio_inspector_t
    use studio_workspace_widgets, only: workspace_widgets_t
    use studio_document, only: studio_document_t, document_mesh_count, document_mesh_source_group, append_document_mesh, document_mesh_valid
    implicit none(type, external)
    private

    public :: studio_t, studio_action_binding_t, workspace_diagnostic_t, canvas_gesture_t, viewport_transform_t
    public :: dp, knot_select_tool, knot_insert_tool, knot_remove_tool, mesh_exterior, mesh_lattice, mesh_layer
    public :: axis_rgb, viewport_world_axes
    public :: studio_mesh_count, append_studio_mesh, adjust_mesh_selection_after_removal
    public :: studio_mesh_source_group, studio_mesh_valid, studio_mesh_cell_label, clear_mesh_selection
    public :: scene_patch_drawn, scene_mesh_drawn, scene_drawn_patch_count
    public :: initialize_state, repaint_all, clear_control_dependency_preview
    public :: area_click, area_box, area_lasso
    integer, parameter, public :: default_assistant_width=328
    integer, parameter :: dp = c_double
    integer, parameter :: knot_select_tool=0,knot_insert_tool=1,knot_remove_tool=2
    integer, parameter :: mesh_exterior=0,mesh_lattice=1,mesh_layer=2
    public :: viewport_shaded, viewport_edges, viewport_wireframe, viewport_fast_shaded
    integer, parameter :: area_click=0,area_box=1,area_lasso=2
    integer, parameter :: axis_rgb(3,3)=reshape([94,221,219, 135,181,247, 234,155,214],[3,3])
    real(dp), parameter :: viewport_world_axes(3,3)=reshape([1._dp,0._dp,0._dp,0._dp,1._dp,0._dp,0._dp,0._dp,1._dp],[3,3])

    type :: canvas_gesture_t
        integer :: pointer_button=Qt_NoButton
        logical :: pointer_displaced=.false., probe_dragging=.false., middle_orbit=.false.
        integer :: domain_pointer=0, domain_button=Qt_NoButton
        real(dp) :: domain_last(2)=0._dp, volume_probe_depth=0._dp
        integer :: area_count=0
        real(dp), allocatable :: area_points(:,:)
        integer :: viewport_hover(4)=0,viewport_mesh_hover(3)=0,viewport_shared_hover=0
        real(dp) :: hover_position(2)=0._dp,frame_time_ms=0._dp
    end type canvas_gesture_t

    type :: control_dependency_source_t
        type(multipatch_controls_t) :: member_controls
        integer :: assembly=0, point=0
        logical :: attempted=.false.
        real(dp) :: position(3)=0._dp
        character(:), allocatable :: error
    end type control_dependency_source_t

    type :: viewport_transform_t
        integer :: axis=0,mode=0,kind=0,group=0,viewport_size(2)=0
        integer,allocatable :: variants(:)
        integer(int64),allocatable :: revisions(:)
        real(dp) :: pivot(3)=0._dp,screen(2)=0._dp,start(2)=0._dp
        real(dp) :: linear(3,3)=viewport_world_axes,offset(3)=0._dp,view(5)=0._dp,basis(3,2)=0._dp
        real(dp) :: plane(2,2)=0._dp,spin_last=0._dp,spin_total=0._dp,amount=0._dp
        logical :: snapped=.false.
    end type viewport_transform_t

    type :: studio_action_binding_t
        type(studio_t), pointer :: state => null()
        integer :: command=0
    end type studio_action_binding_t

    type :: workspace_diagnostic_t
        character(:), allocatable :: message, selection_key
        integer(int64) :: document_revision=0_int64
    end type workspace_diagnostic_t

    type :: studio_t
        type(selection_set_t) :: selection
        type(selection_set_t) :: isolation
        type(qt_string), allocatable :: isolation_mesh_names(:)
        logical :: applying_selection=.false.
        type(workspace_widgets_t) :: widgets
        type(studio_action_binding_t) :: action_bindings(384)
        integer :: action_count=0
        type(studio_document_t) :: document
        type(studio_document_tabs_t) :: project_tabs
        type(workspace_diagnostic_t) :: diagnostic
        type(studio_inspector_t) :: inspector
        type(QMainWindow) :: window
        type(QStackedWidget) :: pages
        type(QWidget) :: home_page
        type(QtfWidget) :: surface, domain, basis, workspace
        type(studio_status_t) :: status
        type(QDoubleSpinBox) :: probe_u, probe_v, probe_w
        type(QLabel) :: probe_label, diagnostic_label, basis_title
        type(QLabel) :: panel_titles(8)
        type(QTableWidget) :: knot_table
        type(QPushButton) :: basis_buttons(2), knot_buttons(3)
        type(QPushButton) :: close_buttons(8)
        type(QToolButton) :: viewport_display,domain_display,basis_display
        type(QAction) :: viewport_presets(2),viewport_layers(viewport_layer_count),domain_overlays(5),basis_options(2)
        type(QAction) :: viewport_grid
        type(QAction) :: toolbar_actions(17), panel_actions(8)
        type(QWidget) :: panels(8), left_column, lower_row, command_ribbon, viewport_controls
        type(QWidget) :: assistant_overlay, assistant_space
        type(QtfWidget) :: viewport_aperture
        type(QSplitter) :: columns, center_rows, assistant_columns
        integer :: assistant_width=default_assistant_width
        real(dp) :: u = 0.35_dp, v = 0.70_dp, w = 0.45_dp
        real(dp) :: domain_yaw=0._dp, domain_pitch=0._dp
        real(dp) :: domain_zoom=1._dp,domain_pan_x=0._dp,domain_pan_y=0._dp
        logical :: domain_edit_knots=.false., domain_knot_details=.false.
        integer :: domain_knot_tool=knot_select_tool
        real(dp) :: yaw = 0.0_dp, pitch = 0.0_dp, zoom = 1.0_dp
        real(dp) :: pan_x = 0.0_dp, pan_y = 0.0_dp
        integer :: feature_pick_mode=0, viewport_feature(4)=0
        integer :: transform_mode=0
        type(viewport_transform_t) :: transform_drag
        type(canvas_gesture_t) :: canvas_gesture
        type(viewport_display_cache) :: display_cache
        integer :: area_selection_tool=area_click
        logical :: snap_enabled=.false., control_snapped=.false.
        logical :: snap_closest=.false.
        logical :: canvas_only=.false., before_canvas_focus=.true., before_canvas_project=.true., before_canvas_ai=.true.
        real(dp) :: last_x = 0.0_dp, last_y = 0.0_dp
        logical :: dragging = .false., updating = .false., focused_surface = .false.
        integer :: selected_control = 0
        logical :: control_dragging = .false.
        real(dp) :: control_preview(3) = 0._dp, control_drag_basis(3,2) = 0._dp
        type(control_dependency_source_t) :: control_dependency_source
        type(multipatch_control_preview_t) :: control_dependency_preview
        logical :: saved_panels(8) = .true.
        integer :: normal_columns(3) = [264,1041,342], normal_center(2) = [488,293]
        logical :: remember_workspace = .false.
        logical :: home_active=.false., resume_assistant=.false.
        integer :: viewport_style=viewport_edges
        logical :: layers(viewport_layer_count) = default_viewport_layers
        logical :: domain_layers(5) = [.true., .true., .true., .false., .false.]
        logical :: basis_legend=.true.,basis_probe_guide=.true.,basis_derivative=.false.
        logical :: world_grid=.true., interactive_quality=.true.
        type(QTimer) :: zoom_quality_timer,viewport_hover_timer
        ! Parameter Space, its knot editor and the basis plot share this direction.
        integer :: entity_kind = 2, knot_axis = 1, surface_tool = 2
        real(dp), allocatable :: live_knots(:), live_knot_groups(:)
        type(single_patch_scene) :: live_scene
        integer :: selected_multipatch=0
        integer :: refinement_preview_kind=0
        logical, allocatable :: refinement_preview_patches(:)
        type(multipatch_controls_t) :: shared_controls
        integer :: shared_controls_assembly=0
        integer, allocatable :: selected_shared_controls(:)
        integer :: selected_mesh=0, mesh_cell(3)=0
        logical :: mesh_pick=.false., mesh_only=.false.
        logical :: mesh_shaded=.true.
        integer :: mesh_view=mesh_exterior, mesh_layer_axis=3, mesh_layer_index=1
        integer :: mesh_section_axis=0
        real(dp) :: mesh_section_fraction=.5_dp
        logical :: physical_knots=.false.
        type(single_patch_info) :: live_probe
        type(single_patch_plot) :: live_plot
        real(dp), allocatable :: live_slice(:,:)
        real(dp) :: live_normalized(3) = -1._dp, live_slice_w = -1._dp
        integer :: live_probe_axis = 0, live_plot_axis = 0
        logical :: live_plot_derivative=.false.
    end type studio_t
contains
    logical function scene_patch_drawn(s,kind,variant) result(drawn)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: kind,variant

        drawn=.false.
        if(s%mesh_only)return
        if(.not.studio_scene_is_visible(s%document%scene,kind,variant,s%document%patches))return
        if(.not.s%document%scene%catalog(kind)%entries(variant)%patch%valid)return
        drawn=s%isolation%count()==0.or.s%isolation%contains(selection_target_t(selection_patch,variant,kind))
    end function scene_patch_drawn

    logical function scene_mesh_drawn(s,id) result(drawn)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: id
        integer :: i

        drawn=.false.
        if(id<1.or.id>studio_mesh_count(s))return
        if(.not.s%document%meshes(id)%visible)return
        if(.not.studio_mesh_valid(s,id))return
        if(s%isolation%count()==0)then
            drawn=.true.
            return
        end if
        associate(mesh=>s%document%meshes(id))
            if(allocated(s%isolation_mesh_names).and.allocated(mesh%name))then
                do i=1,size(s%isolation_mesh_names)
                    if(s%isolation_mesh_names(i)%s/=mesh%name)cycle
                    drawn=.true.
                    return
                end do
            end if
            if(.not.allocated(mesh%parts))return
            if(size(mesh%parts)==0)return
            drawn=all([(s%isolation%contains(selection_target_t(selection_patch,mesh%parts(i)%variant,mesh%kind)),i=1,size(mesh%parts))])
        end associate
    end function scene_mesh_drawn

    integer function scene_drawn_patch_count(s) result(n)
        type(studio_t), intent(in) :: s
        integer :: i

        n=0
        if(s%mesh_only)return
        if(s%isolation%count()==0)then
            n=studio_scene_count(s%document%scene)
            return
        end if
        do i=1,s%isolation%count()
            associate(target=>s%isolation%items(i))
                if(target%category/=selection_patch)cycle
                if(scene_patch_drawn(s,target%rank,target%id))n=n+1
            end associate
        end do
    end function scene_drawn_patch_count

    integer function studio_mesh_count(s) result(n)
        type(studio_t), intent(in) :: s

        n=document_mesh_count(s%document)
    end function studio_mesh_count

    subroutine append_studio_mesh(s,mesh)
        type(studio_t), intent(inout) :: s
        type(studio_mesh_t), intent(inout) :: mesh

        call append_document_mesh(s%document,mesh)
    end subroutine append_studio_mesh

    subroutine adjust_mesh_selection_after_removal(s,id)
        type(studio_t), intent(inout) :: s
        integer, intent(in) :: id

        if(s%selected_mesh==id)call clear_mesh_selection(s)
        if(s%selected_mesh>id)s%selected_mesh=s%selected_mesh-1
    end subroutine adjust_mesh_selection_after_removal

    integer function studio_mesh_source_group(s,id) result(group)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: id

        group=document_mesh_source_group(s%document,id)
    end function studio_mesh_source_group

    logical function studio_mesh_valid(s,id) result(valid)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: id

        valid=document_mesh_valid(s%document,id)
    end function studio_mesh_valid

    function studio_mesh_cell_label(s,target) result(label)
        type(studio_t), intent(in) :: s
        type(selection_target_t), intent(in), optional :: target
        character(:), allocatable :: label
        character(48) :: number
        integer :: id,part,cell,i

        label=''
        id=s%selected_mesh
        part=s%mesh_cell(1)
        cell=s%mesh_cell(2)
        if(present(target))then
            id=target%id
            part=target%part
            cell=target%item
        end if
        if(.not.studio_mesh_valid(s,id))return
        if(.not.s%document%meshes(id)%visible)return
        associate(mesh=>s%document%meshes(id))
            if(part<1.or.part>size(mesh%parts).or.cell<1)return
            if(cell>size(mesh%parts(part)%cells,1))return
            do i=1,part-1
                cell=cell+size(mesh%parts(i)%cells,1)
            end do
            write(number,'(a,i0)') 'Cell ',cell
            label=trim(number)
            if(mesh%group>0)then
                write(number,'(a,i0)') ' · Patch ',part
                label=label//trim(number)
            end if
        end associate
    end function studio_mesh_cell_label

    subroutine clear_mesh_selection(s)
        type(studio_t), intent(inout) :: s
        s%selected_mesh=0
        s%mesh_cell=0
        s%mesh_pick=.false.
    end subroutine clear_mesh_selection

    subroutine initialize_state(s)
        type(studio_t), intent(inout) :: s
        call initialize_studio_patch_catalog(s%document%patches)
        ! Reset presentation controls while retaining native patches and Qt handles.
        s%u=.35_dp; s%v=.70_dp; s%w=.45_dp
        s%domain_yaw=0._dp; s%domain_pitch=0._dp
        s%domain_zoom=1._dp; s%domain_pan_x=0._dp; s%domain_pan_y=0._dp
        s%domain_edit_knots=.false.; s%domain_knot_details=.false.
        s%domain_knot_tool=knot_select_tool
        s%yaw=0._dp; s%pitch=0._dp; s%zoom=1._dp
        s%pan_x=0._dp; s%pan_y=0._dp
        call frame_studio_scene(s%document%scene)
        s%last_x=0._dp; s%last_y=0._dp; s%dragging=.false.
        if(allocated(s%selection%items))deallocate(s%selection%items)
        s%selected_control=0; s%control_dragging=.false.
        call clear_control_dependency_preview(s)
        s%shared_controls=multipatch_controls_t()
        s%shared_controls_assembly=0
        if(allocated(s%selected_shared_controls))deallocate(s%selected_shared_controls)
        s%viewport_feature=0; s%control_snapped=.false.
        s%transform_mode=0
        s%transform_drag=viewport_transform_t()
        s%canvas_gesture=canvas_gesture_t()
        s%area_selection_tool=area_click
        call clear_mesh_selection(s)
        s%mesh_only=.false.
        s%mesh_shaded=.true.
        s%mesh_view=mesh_exterior; s%mesh_layer_axis=3; s%mesh_layer_index=1
        s%mesh_section_axis=0; s%mesh_section_fraction=.5_dp
        s%viewport_style=viewport_edges
        s%physical_knots=.false.
        s%layers=default_viewport_layers
        s%domain_layers=[.true.,.true.,.true.,.false.,.false.]
        s%basis_legend=.true.; s%basis_probe_guide=.true.; s%basis_derivative=.false.
        s%world_grid=.true.
        s%interactive_quality=.true.
        if(.not.s%zoom_quality_timer%isNull())call s%zoom_quality_timer%stop()
        s%entity_kind=2; s%knot_axis=1; s%surface_tool=2
        s%live_scene%valid=.false.; s%live_probe%valid=.false.; s%live_plot%valid=.false.
    end subroutine initialize_state

    subroutine clear_control_dependency_preview(s)
        type(studio_t), intent(inout) :: s

        s%control_dependency_source=control_dependency_source_t()
        s%control_dependency_preview=multipatch_control_preview_t()
    end subroutine clear_control_dependency_preview
    subroutine repaint_all(s)
        type(studio_t), intent(inout) :: s
        call s%surface%update()
        call s%domain%update()
        call s%basis%update()
    end subroutine repaint_all
end module studio_state
