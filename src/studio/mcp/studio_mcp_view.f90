module studio_mcp_view
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use qt, only: QJsonArray, QJsonObject, QJsonValue_Array, QJsonValue_Bool, QJsonValue_Double, &
        QJsonValue_Object, QJsonValue_String, QTabWidget, qt_string
    use studio_state, only: studio_t, dp, studio_mesh_count, studio_mesh_valid, mesh_exterior, mesh_layer, &
        area_click, area_box, area_lasso
    use studio_assistant, only: studio_assistant_is_visible
    use studio_inspector, only: update_inspector_section
    use studio_canvas_input, only: min_domain_zoom,max_domain_zoom,domain_pitch_limit
    use studio_domain_projection, only: cube_yaw,cube_pitch
    use studio_scene, only: project_patch_exists
    use studio_multipatch_collection, only: count_multipatches
    use studio_ui, only: select_patch, select_project_mesh, select_multipatch_group, set_workspace_material, select_workspace_mesh_cell, &
        set_patch_visibility, set_mesh_visibility, set_parameter_probe, set_workspace_mesh_view, set_workspace_mesh_section, &
        has_visible_meshes, invoke_command
    use studio_commands, only: toggle_assistant_command, canvas_command, toggle_panel_command_base, toggle_inspector_command, &
        knots_u_command, viewport_layer_base, viewport_layer_count, viewport_style_names, viewport_style_commands, &
        physical_knots_command, world_grid_command, basis_functions_command, basis_derivative_command, &
        basis_legend_command, basis_probe_command, project_ribbon_tab, create_ribbon_tab, model_ribbon_tab, &
        patch_ribbon_tab, multipatch_ribbon_tab, mesh_ribbon_tab, view_ribbon_tab, interactive_quality_command, &
        mesh_pick_command, mesh_shaded_command, mesh_only_command, select_tool_command, orbit_tool_command, &
        box_selection_command, lasso_selection_command, pick_mode_base, pick_controls_command, pick_multipatch_command, &
        snap_off_command, snap_vertex_command, snap_closest_command, show_knots_command, domain_probe_command, domain_select_command
    use studio_mcp_arguments, only: only_fields, bounded_integer, integer_array_arg, optional_boolean, number_arg, axis_number
    use studio_mcp_inspection, only: describe_selection, describe_probe, patch_identity, assembly_identity, mesh_identity, &
        active_rank, patch_kind
    use studio_viewport_camera, only: viewport_camera, workspace_view_bounds, prepare_camera, frame_workspace_scene, &
        live_yaw_offset, live_pitch_offset, wrap_view_angle, viewport_zoom_range
    use studio_project_tree, only: scene_filter_all,scene_filter_geometry,scene_filter_multipatches,scene_filter_meshes, &
        scene_filter_hidden,scene_filter_stale,scene_filter_curves,scene_filter_surfaces,scene_filter_volumes
    implicit none(type, external)
    private

    public :: select_object, set_material, set_visibility, set_probe, select_mesh_cell
    public :: inspect_view, set_camera, set_parameter_view, set_mesh_view, set_workspace
    real(dp), parameter :: radians_per_degree=acos(-1._dp)/180._dp
    character(8), parameter :: mesh_modes(mesh_exterior:mesh_layer)=[character(8)::'exterior','lattice','layer']
    character(*), parameter :: axes='UVW'
    character(10), parameter :: tab_names(*)=[character(10)::'project','create','model','patch','multipatch','mesh','view']
    integer, parameter :: tab_ids(*)=[project_ribbon_tab,create_ribbon_tab,model_ribbon_tab,patch_ribbon_tab, &
        multipatch_ribbon_tab,mesh_ribbon_tab,view_ribbon_tab]
    character(15), parameter :: panel_names(*)=[character(15)::'scene','inspector','parameter_space','basis_plot','assistant']
    integer, parameter :: panel_ids(*)=[1,8,5,6]
    character(15), parameter :: inspector_names(*)=[character(15)::'selection_probe','diagnostics']
    character(15), parameter :: display_names(*)=[character(15)::'shaded','wireframe','control_net','control_points', &
        'parameter_lines','knot_boundaries','ground_grid','interactive','mesh_shading','mesh_only','mesh_picking']
    character(12), parameter :: display_styles(*)=[character(12)::'shaded','shaded_edges','wireframe','fast_shaded']
    character(4), parameter :: section_axes(0:3)=[character(4)::'off','X','Y','Z']
    character(6), parameter :: viewport_tools(*)=[character(6)::'select','orbit']
    character(5), parameter :: selection_tools(area_click:area_lasso)=[character(5)::'click','box','lasso']
    character(10), parameter :: pick_modes(0:5)=[character(10)::'automatic','face','edge','vertex','control','multipatch']
    character(7), parameter :: snap_modes(*)=[character(7)::'off','vertex','closest']
    character(11), parameter :: domain_tools(*)=[character(11)::'probe','select_knot','insert_knot','remove_knot']
    integer, parameter :: viewport_commands(*)=[select_tool_command,orbit_tool_command]
    integer, parameter :: selection_commands(*)=[select_tool_command,box_selection_command,lasso_selection_command]
    integer, parameter :: pick_commands(*)=[pick_mode_base,pick_mode_base+1,pick_mode_base+2,pick_mode_base+3, &
        pick_controls_command,pick_multipatch_command]
    integer, parameter :: snap_commands(*)=[snap_off_command,snap_vertex_command,snap_closest_command]
    integer, parameter :: display_mesh_shading=9,display_mesh_only=10,display_mesh_picking=11
    character(16), parameter :: domain_names(*)=[character(16)::'grid','knots','elements','basis_support','selected_element']
    character(10), parameter :: basis_names(*)=[character(10)::'derivative','legend','probe']
    character(12), parameter :: scene_filters(*)=[character(12)::'all','geometry','multipatches','meshes','hidden', &
        'needs_update','curves','surfaces','volumes']
    integer, parameter :: scene_filter_ids(*)=[scene_filter_all,scene_filter_geometry,scene_filter_multipatches, &
        scene_filter_meshes,scene_filter_hidden,scene_filter_stale,scene_filter_curves,scene_filter_surfaces,scene_filter_volumes]

    type :: workspace_toggle
        logical :: specified=.false., value=.false.
    end type workspace_toggle

contains

    subroutine select_object(s,args,result,error)
        type(studio_t),intent(inout),target :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        character(:),allocatable :: kind_name
        integer :: kind,id

        error='Select requires an object kind and a positive ID from inspect.'
        if(.not.only_fields(args,'|kind|id|'))return
        kind_name=args%stringValue('kind')
        id=bounded_integer(args,'id',1,huge(1))
        if(id<1)return
        select case(kind_name)
        case('mesh')
            error='That mesh is unavailable.'
            if(id>studio_mesh_count(s))return
            call select_project_mesh(s,id)
        case('assembly')
            error='That assembly is unavailable.'
            if(id>count_multipatches(s%document%multipatches))return
            call select_multipatch_group(s,id,reveal=.false.)
        case default
            kind=patch_kind(kind_name)
            if(kind==0)return
            error='That patch is unavailable.'
            if(.not.project_patch_exists(s%document%show_templates,kind,id,patches=s%document%patches))return
            call select_patch(s,kind,id,reveal=.false.)
        end select
        call describe_selection(s,result)
        error=''
    end subroutine select_object

    subroutine set_material(s,args,result,error)
        type(studio_t),intent(inout),target :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        integer,allocatable :: rgb(:),opacity
        logical :: changed

        error='Supply RGB channels or opacity with integer values from 0 to 255.'
        if(.not.only_fields(args,'|rgb|opacity|').or.args%isEmpty())return
        if(args%contains('rgb'))then
            if(.not.integer_array_arg(args,'rgb',rgb,3,3,255,minimum=0))return
        end if
        if(args%contains('opacity'))then
            opacity=bounded_integer(args,'opacity',0,255)
            if(opacity<0)return
        end if
        call set_workspace_material(s,changed,error,rgb,opacity)
        if(len(error)>0)return
        call result%insert('changed',changed)
        call describe_selection(s,result)
    end subroutine set_material

    subroutine set_visibility(s,args,result,error)
        type(studio_t),intent(inout) :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        type(QJsonObject) :: item
        character(:),allocatable :: kind_name
        logical :: visible,changed
        integer :: kind,id

        error='Set visibility with an object kind, its ID and a Boolean visible value.'
        if(.not.only_fields(args,'|kind|id|visible|'))return
        if(args%valueType('visible')/=QJsonValue_Bool)return
        id=bounded_integer(args,'id',1,huge(1))
        if(id<1)return
        visible=args%boolValue('visible')
        kind_name=args%stringValue('kind')
        select case(kind_name)
        case('assembly')
            error='That assembly is unavailable.'
            if(id>count_multipatches(s%document%multipatches))return
            associate(definition=>s%document%multipatches%groups(id)%definition)
                call set_patch_visibility(s,definition%kind,definition%members,visible,changed,error)
            end associate
            if(len(error)>0)return
            item=assembly_identity(s,id)
        case('mesh')
            error='That mesh is unavailable or stale. Regenerate a stale mesh before showing it.'
            if(.not.studio_mesh_valid(s,id))return
            call set_mesh_visibility(s,id,visible,changed,error)
            if(len(error)>0)return
            item=mesh_identity(s,id)
        case default
            error='That patch is unavailable in this project.'
            kind=patch_kind(kind_name)
            if(.not.project_patch_exists(s%document%show_templates,kind,id,patches=s%document%patches))return
            call set_patch_visibility(s,kind,[id],visible,changed,error)
            if(len(error)>0)return
            item=patch_identity(s,kind,id)
        end select
        call result%insert('changed',changed)
        call result%insert('target',item)
        call item%delete()
        error=''
    end subroutine set_visibility

    subroutine set_probe(s,args,result,error)
        type(studio_t),intent(inout) :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        type(QJsonArray) :: values
        real(dp) :: normalized(3),parameters(3),low(3),high(3)
        integer :: i,kind
        logical :: valid

        kind=active_rank(s)
        error='Select a patch before setting the parameter probe.'
        if(kind<1)return
        error='Supply one parameter per active direction in U, V, W order.'
        if(.not.only_fields(args,'|parameters|normalized|').or..not.optional_boolean(args,'normalized'))return
        if(args%valueType('parameters')/=QJsonValue_Array)return
        values=args%arrayValue('parameters')
        valid=values%size()==kind
        if(valid)then
            do i=1,kind
                if(values%typeAt(i-1)/=QJsonValue_Double)valid=.false.
                parameters(i)=values%doubleAt(i-1)
                if(.not.ieee_is_finite(parameters(i)))valid=.false.
            end do
        end if
        call values%delete()
        if(.not.valid)return
        normalized=[s%u,s%v,s%w]
        if(args%boolValue('normalized'))then
            low=0._dp; high=1._dp
        else
            low=s%live_scene%bounds(1,:); high=s%live_scene%bounds(2,:)
        end if
        error='Every probe parameter must lie within its patch interval.'
        if(any(parameters(:kind)<low(:kind)).or.any(parameters(:kind)>high(:kind)))return
        normalized(:kind)=(parameters(:kind)-low(:kind))/(high(:kind)-low(:kind))
        call set_parameter_probe(s,normalized)
        call describe_probe(s,kind,result)
        error=''
    end subroutine set_probe

    subroutine select_mesh_cell(s,args,result,error)
        type(studio_t),intent(inout),target :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        type(QJsonArray) :: vertices,point
        integer :: id,number,part,cell,i,node

        error='Supply a mesh ID and one-based cell number.'
        if(.not.only_fields(args,'|id|cell|'))return
        id=bounded_integer(args,'id',1,huge(1))
        number=bounded_integer(args,'cell',1,huge(1))
        if(min(id,number)<1)return
        call select_workspace_mesh_cell(s,id,number,part,cell,error)
        if(len(error)>0)return
        call result%insert('mesh_id',id)
        call result%insert('cell',number)
        call result%insert('part',part)
        call result%insert('part_cell',cell)
        vertices=QJsonArray()
        associate(data=>s%document%meshes(id)%parts(part))
            do i=1,size(data%cells,2)
                node=data%cells(cell,i)
                point=QJsonArray()
                call point%appendDoubles(real(data%points(node,:),dp))
                call vertices%append(point)
                call point%delete()
            end do
        end associate
        call result%insert('vertices',vertices)
        call vertices%delete()
    end subroutine select_mesh_cell

    subroutine inspect_view(s,args,result,error)
        type(studio_t),intent(in) :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error

        error='Inspect view takes an empty object.'
        if(.not.args%isEmpty())return
        error='Open a ForCAD Studio project before inspecting its view.'
        if(s%home_active.or.s%surface%isNull())return
        call describe_parameter_view(s,result)
        call describe_view(s,result)
        call describe_workspace(s,result)
        error=''
    end subroutine inspect_view

    subroutine set_camera(s,args,result,error)
        type(studio_t),intent(inout) :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        real(dp) :: yaw,pitch,zoom,pan(2),degrees,limits(2)
        logical :: rotating

        error='Supply at least one camera field: yaw_degrees, pitch_degrees, zoom, pan_x, or pan_y.'
        if(.not.only_fields(args,'|yaw_degrees|pitch_degrees|zoom|pan_x|pan_y|').or.args%isEmpty())return
        yaw=s%yaw
        pitch=s%pitch
        zoom=s%zoom
        pan=[s%pan_x,s%pan_y]
        if(args%contains('yaw_degrees'))then
            if(.not.number_arg(args,'yaw_degrees',degrees))return
            yaw=modulo(degrees,360._dp)*radians_per_degree-live_yaw_offset
        end if
        if(args%contains('pitch_degrees'))then
            error='Camera pitch must be a finite angle in degrees.'
            if(.not.number_arg(args,'pitch_degrees',degrees))return
            pitch=wrap_view_angle(modulo(degrees,360._dp)*radians_per_degree)-live_pitch_offset
        end if
        if(args%contains('zoom'))then
            error='Camera zoom must be positive and within the viewport precision range reported by inspect_view.'
            if(.not.number_arg(args,'zoom',zoom))return
            limits=viewport_zoom_range(s)
            if(zoom<limits(1).or.zoom>limits(2))return
        end if
        error='Camera pan must use finite viewport pixel offsets.'
        if(args%contains('pan_x'))then
            if(.not.number_arg(args,'pan_x',pan(1)))return
        end if
        if(args%contains('pan_y'))then
            if(.not.number_arg(args,'pan_y',pan(2)))return
        end if
        rotating=yaw/=s%yaw.or.pitch/=s%pitch
        s%yaw=yaw
        s%pitch=pitch
        if(rotating)call frame_workspace_scene(s)
        s%zoom=zoom
        s%pan_x=pan(1)
        s%pan_y=pan(2)
        call s%surface%update()
        call describe_view(s,result)
        error=''
    end subroutine set_camera

    subroutine set_parameter_view(s,args,result,error)
        type(studio_t),intent(inout) :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        real(dp) :: yaw,pitch,zoom,pan(2),degrees
        integer :: kind

        kind=active_rank(s)
        error='Select a single patch to set the Parameter Space view.'
        if(kind<1.or.s%selected_multipatch>0.or.s%selected_mesh>0.or.s%selection%count()>1)return
        error='Supply zoom, pan_x, pan_y, yaw_degrees or pitch_degrees for Parameter Space.'
        if(.not.only_fields(args,'|zoom|pan_x|pan_y|yaw_degrees|pitch_degrees|').or.args%isEmpty())return
        if(kind/=3.and.(args%contains('yaw_degrees').or.args%contains('pitch_degrees')))then
            error='Parameter Space orbit requires a volume.'
            return
        end if
        yaw=s%domain_yaw; pitch=s%domain_pitch; zoom=s%domain_zoom
        pan=[s%domain_pan_x,s%domain_pan_y]
        error='Parameter Space camera fields must be finite numbers within the ranges from inspect_view.'
        if(args%contains('yaw_degrees'))then
            if(.not.number_arg(args,'yaw_degrees',degrees))return
            yaw=modulo(degrees,360._dp)*radians_per_degree-cube_yaw
        end if
        if(args%contains('pitch_degrees'))then
            if(.not.number_arg(args,'pitch_degrees',degrees))return
            if(abs(degrees)>domain_pitch_limit/radians_per_degree)return
            pitch=degrees*radians_per_degree-cube_pitch
        end if
        if(args%contains('zoom'))then
            if(.not.number_arg(args,'zoom',zoom))return
            if(zoom<min_domain_zoom.or.zoom>max_domain_zoom)return
        end if
        if(args%contains('pan_x'))then
            if(.not.number_arg(args,'pan_x',pan(1)))return
        end if
        if(args%contains('pan_y'))then
            if(.not.number_arg(args,'pan_y',pan(2)))return
        end if
        s%domain_yaw=yaw; s%domain_pitch=pitch; s%domain_zoom=zoom
        s%domain_pan_x=pan(1); s%domain_pan_y=pan(2)
        call s%domain%update()
        call describe_parameter_view(s,result)
        error=''
    end subroutine set_parameter_view

    subroutine describe_parameter_view(s,result)
        type(studio_t),intent(in) :: s
        type(QJsonObject),intent(inout) :: result
        type(QJsonObject) :: camera
        type(QJsonArray) :: limits
        integer :: kind

        kind=active_rank(s)
        if(kind<1.or.s%selected_mesh>0.or.s%selected_multipatch>0.or.s%selection%count()>1)return
        camera=QJsonObject()
        call camera%insert('dimensions',kind)
        call camera%insert('zoom',s%domain_zoom)
        call camera%insert('pan_x',s%domain_pan_x)
        call camera%insert('pan_y',s%domain_pan_y)
        limits=QJsonArray()
        call limits%appendDoubles([min_domain_zoom,max_domain_zoom])
        call camera%insert('zoom_range',limits)
        call limits%delete()
        if(kind==3)then
            call camera%insert('yaw_degrees',modulo((cube_yaw+s%domain_yaw)/radians_per_degree,360._dp))
            call camera%insert('pitch_degrees',(cube_pitch+s%domain_pitch)/radians_per_degree)
            limits=QJsonArray()
            call limits%appendDoubles([-domain_pitch_limit,domain_pitch_limit]/radians_per_degree)
            call camera%insert('pitch_range',limits)
            call limits%delete()
        end if
        call result%insert('parameter_camera',camera)
        call camera%delete()
    end subroutine describe_parameter_view

    subroutine set_mesh_view(s,args,result,error)
        type(studio_t),intent(inout) :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        integer,allocatable :: axis,layer
        integer :: mode,section,percent
        logical :: set_presentation,set_section

        error='Supply a mesh presentation or section settings.'
        if(.not.only_fields(args,'|mode|axis|layer|section_axis|section_percent|').or.args%isEmpty())return
        section=workspace_choice(args,'section_axis',section_axes)
        percent=bounded_integer(args,'section_percent',0,100,fallback=nint(100._dp*s%mesh_section_fraction))
        error='Choose off, X, Y or Z and a section percentage from 0 to 100.'
        if(section<0.or.percent<0)return
        set_section=section>0.or.args%contains('section_percent')
        if(section==0)section=s%mesh_section_axis+1
        if(set_section)then
            error='Show a current volume mesh before changing its section.'
            if(.not.has_visible_meshes(s,3))return
        end if

        set_presentation=args%contains('mode').or.args%contains('axis').or.args%contains('layer')
        mode=s%mesh_view
        error='Choose exterior, lattice, or layer. For layer mode, optionally supply axis U, V, or W and a positive layer.'
        if(args%contains('mode'))mode=findloc(mesh_modes,args%stringValue('mode'),dim=1)-1
        if(mode<mesh_exterior)return
        if(mode/=mesh_layer.and.(args%contains('axis').or.args%contains('layer')))return
        if(args%contains('axis'))then
            axis=axis_number(args%stringValue('axis'),axes)
            if(axis<1)return
        end if
        if(args%contains('layer'))then
            layer=bounded_integer(args,'layer',1,huge(1))
            if(layer<1)return
        end if
        if(set_presentation)then
            call set_workspace_mesh_view(s,mode,axis,layer,error)
            if(len(error)>0)return
        end if
        if(set_section)then
            call set_workspace_mesh_section(s,section-1,percent,error)
            if(len(error)>0)return
        end if
        call describe_view(s,result)
        error=''
    end subroutine set_mesh_view

    subroutine set_workspace(s,args,result,error)
        type(studio_t),intent(inout),target :: s
        type(QJsonObject),intent(in) :: args
        type(QJsonObject),intent(inout) :: result
        character(:),allocatable,intent(out) :: error
        type(workspace_toggle) :: panels(size(panel_names)), display(size(display_names)), domain(size(domain_names)), &
            basis(size(basis_names)), inspector(size(inspector_names))
        type(QTabWidget) :: tabs
        logical :: canvas, show_panels, mesh_only
        integer :: i, tab, filter, axis, style, viewport, selection, pick, snap, parameter_tool
        character(:),allocatable :: option

        error='Open a ForCAD Studio project before changing its workspace.'
        if(s%home_active.or.s%surface%isNull())return
        error='Supply workspace options for the ribbon, panels, canvas, Scene, display, selection, snap, or parameter tools.'
        if(.not.only_fields(args, '|tab|panels|canvas_only|scene_search|scene_filter|display|domain|basis|parameter_axis|inspector|'// &
            'display_style|viewport_tool|selection_tool|pick_mode|snap|domain_tool|').or.args%isEmpty())return

        style=workspace_choice(args,'display_style',display_styles)
        viewport=workspace_choice(args,'viewport_tool',viewport_tools)
        selection=workspace_choice(args,'selection_tool',selection_tools)
        pick=workspace_choice(args,'pick_mode',pick_modes)
        snap=workspace_choice(args,'snap',snap_modes)
        parameter_tool=workspace_choice(args,'domain_tool',domain_tools)
        error='Choose display, viewport, selection, pick, snap and domain options from the values in the tool description.'
        if(min(style,viewport,selection,pick,snap,parameter_tool)<0)return
        error='Use Select mode for box or lasso selection.'
        if(viewport==2.and.selection>1)return

        tab=0
        if(args%contains('tab'))then
            error='Unknown ribbon tab. Choose project, create, model, patch, multipatch, mesh, or view.'
            option=args%stringValue('tab')
            if(len(option)/=len_trim(option))return
            tab=findloc(tab_names,option,dim=1)
            if(tab==0)return
            error='The ribbon tabs are unavailable.'
            call s%workspace%findChild('studioRibbonTabs',tabs)
            if(tabs%isNull())return
        end if
        error='Canvas-only mode requires a Boolean value.'
        if(.not.optional_boolean(args,'canvas_only'))return
        if(args%contains('scene_search'))then
            error='Scene search requires a string.'
            if(args%valueType('scene_search')/=QJsonValue_String)return
        end if
        filter=0
        if(args%contains('scene_filter'))then
            error='Choose a Scene filter from all, geometry, multipatches, meshes, hidden, needs_update, curves, surfaces or volumes.'
            option=args%stringValue('scene_filter')
            if(len(option)/=len_trim(option))return
            filter=findloc(scene_filters,option,dim=1)
            if(filter==0)return
        end if
        axis=0
        if(args%contains('parameter_axis'))then
            error='Select a patch and choose an available parameter axis (U, V or W).'
            axis=axis_number(args%stringValue('parameter_axis'),axes)
            if(axis<1.or.axis>active_rank(s).or.s%selected_multipatch>0.or.s%selected_mesh>0)return
        end if
        call read_workspace_toggles(args,'panels',panel_names,panels,error)
        if(len(error)>0)return
        call read_workspace_toggles(args,'display',display_names,display,error)
        if(len(error)>0)return
        call read_workspace_toggles(args,'domain',domain_names,domain,error)
        if(len(error)>0)return
        call read_workspace_toggles(args,'basis',basis_names,basis,error)
        if(len(error)>0)return
        call read_workspace_toggles(args,'inspector',inspector_names,inspector,error)
        if(len(error)>0)return

        show_panels=any(panels%specified.and.panels%value)
        canvas=s%canvas_only.and..not.show_panels
        if(args%contains('canvas_only'))canvas=args%boolValue('canvas_only')
        error='Turn off Canvas only to open a workspace panel.'
        if(canvas.and.show_panels)return
        if(s%selected_multipatch>0)then
            error='Select a member patch to open or change its parameter space and basis plot.'
            if(any(panels(3:4)%specified.and.panels(3:4)%value).or.any(domain%specified).or.any(basis%specified))return
        end if
        if(parameter_tool>0)then
            error='Open Parameter Space for a single patch before changing its tools.'
            if(canvas.or.s%selected_multipatch>0.or.s%selected_mesh>0.or.s%selection%count()>1.or.active_rank(s)<1)return
            if(parameter_tool==4)then
                error='Choose an axis with an interior knot before activating the removal tool.'
                if(.not.has_interior_knots(s,merge(axis,s%knot_axis,axis>0)))return
            end if
        end if
        mesh_only=s%mesh_only
        if(display(display_mesh_only)%specified)mesh_only=display(display_mesh_only)%value
        error='Show geometry alongside meshes before changing knot boundaries.'
        if(mesh_only.and.needs_toggle(display(6),s%physical_knots))return
        error='Show a mesh before changing mesh display or picking.'
        if(mesh_only.and..not.s%mesh_only.or.needs_toggle(display(display_mesh_picking),s%mesh_pick))then
            if(.not.has_visible_meshes(s))return
        end if
        error='Show a current surface or volume mesh before changing its shading.'
        if(needs_toggle(display(display_mesh_shading),s%mesh_shaded).and..not.has_visible_meshes(s,2))return
        error='Mesh cell picking uses Select and Click mode.'
        if(display(display_mesh_picking)%specified.and.display(display_mesh_picking)%value)then
            if(viewport==2.or.selection>1)return
        end if

        if(canvas.neqv.s%canvas_only)call invoke_command(s,canvas_command)
        if(tab/=0)call tabs%setCurrentIndex(tab_ids(tab))
        do i=1,size(inspector)
            if(.not.needs_toggle(inspector(i),s%inspector%headers(i)%isChecked()))cycle
            call s%inspector%headers(i)%setChecked(inspector(i)%value)
            call update_inspector_section(s%inspector,i)
        end do
        if(args%contains('scene_search'))call s%widgets%project_browser%search%setText(args%stringValue('scene_search'))
        if(filter>0)then
            if(s%widgets%project_browser%filter_mode/=scene_filter_ids(filter))then
                call s%widgets%project_browser%filters(scene_filter_ids(filter)+1)%trigger()
            end if
        end if
        if(axis>0.and.axis/=s%knot_axis)call invoke_command(s,knots_u_command+axis-1)
        if(style>0.and.style/=s%viewport_style)call invoke_command(s,viewport_style_commands(style))
        if(needs_toggle(display(display_mesh_only),s%mesh_only))call invoke_command(s,mesh_only_command)
        if(pick>0)then
            if(pick-1/=s%feature_pick_mode.or.s%mesh_pick.or.s%surface_tool/=1)call invoke_command(s,pick_commands(pick))
        end if
        if(viewport>0.and.viewport/=s%surface_tool)call invoke_command(s,viewport_commands(viewport))
        if(selection>0.and.viewport/=2)then
            if(selection-1/=s%area_selection_tool.or.s%surface_tool/=1)call invoke_command(s,selection_commands(selection))
        end if
        if(snap>0.and.snap/=current_snap_mode(s))call invoke_command(s,snap_commands(snap))
        if(needs_toggle(display(8),s%interactive_quality))call invoke_command(s,interactive_quality_command)
        if(needs_toggle(display(display_mesh_shading),s%mesh_shaded))call invoke_command(s,mesh_shaded_command)
        if(needs_toggle(display(display_mesh_picking),s%mesh_pick))call invoke_command(s,mesh_pick_command)
        if(parameter_tool>0.and.parameter_tool/=current_domain_tool(s))then
            if(parameter_tool==1)then
                call invoke_command(s,domain_probe_command)
            else
                if(.not.s%domain_edit_knots)call invoke_command(s,show_knots_command)
                call invoke_command(s,domain_select_command+parameter_tool-2)
            end if
        end if
        do i=1,viewport_layer_count
            if(needs_toggle(display(i),s%layers(i)))call invoke_command(s,viewport_layer_base+i-1)
        end do
        if(needs_toggle(display(6),s%physical_knots))call invoke_command(s,physical_knots_command)
        if(needs_toggle(display(7),s%world_grid))call invoke_command(s,world_grid_command)
        do i=1,size(domain)
            if(needs_toggle(domain(i),s%domain_layers(i)))call s%domain_overlays(i)%setChecked(domain(i)%value)
        end do
        if(needs_toggle(basis(1),s%basis_derivative)) then
          call invoke_command(s,merge(basis_derivative_command,basis_functions_command,basis(1)%value))
        end if
        if(needs_toggle(basis(2),s%basis_legend))call invoke_command(s,basis_legend_command)
        if(needs_toggle(basis(3),s%basis_probe_guide))call invoke_command(s,basis_probe_command)
        do i=1,size(panel_ids)
            if(.not.needs_toggle(panels(i),s%panels(panel_ids(i))%isVisible()))cycle
            call invoke_command(s,merge(toggle_inspector_command,toggle_panel_command_base+panel_ids(i),panel_ids(i)==8))
        end do
        if(needs_toggle(panels(size(panels)),studio_assistant_is_visible(s)))call invoke_command(s,toggle_assistant_command)
        call describe_workspace(s,result)
        error=''
    end subroutine set_workspace

    logical function has_interior_knots(s,axis) result(available)
        type(studio_t), intent(in) :: s
        integer, intent(in) :: axis

        available=.false.
        if(.not.allocated(s%live_scene%axes(axis)%knots))return
        associate(knots=>s%live_scene%axes(axis)%knots,bounds=>s%live_scene%bounds(:,axis))
            available=any(knots>bounds(1).and.knots<bounds(2))
        end associate
    end function has_interior_knots

    integer function workspace_choice(args,key,names) result(choice)
        type(QJsonObject), intent(in) :: args
        character(*), intent(in) :: key,names(:)
        character(:), allocatable :: value

        choice=0
        if(.not.args%contains(key))return
        choice=-1
        if(args%valueType(key)/=QJsonValue_String)return
        value=args%stringValue(key)
        if(len(value)/=len_trim(value))return
        choice=findloc(names,value,dim=1)
        if(choice==0)choice=-1
    end function workspace_choice

    integer function current_snap_mode(s) result(mode)
        type(studio_t), intent(in) :: s

        mode=1
        if(s%snap_enabled)mode=merge(3,2,s%snap_closest)
    end function current_snap_mode

    integer function current_domain_tool(s) result(tool)
        type(studio_t), intent(in) :: s

        tool=1
        if(s%domain_edit_knots)tool=2+s%domain_knot_tool
    end function current_domain_tool

    subroutine read_workspace_toggles(args,group,names,toggles,error)
        type(QJsonObject),intent(in) :: args
        character(*),intent(in) :: group,names(:)
        type(workspace_toggle),intent(out) :: toggles(:)
        character(:),allocatable,intent(out) :: error
        type(QJsonObject) :: values
        type(qt_string),allocatable :: keys(:)
        integer :: i, index

        error=''
        if(.not.args%contains(group))return
        if(args%valueType(group)/=QJsonValue_Object)then
            error=group//' options require an object of Boolean values.'
            return
        end if
        values=args%objectValue(group)
        keys=values%keys()
        do i=1,size(keys)
            index=findloc(names,keys(i)%s,dim=1)
            if(index==0.or.len(keys(i)%s)/=len_trim(keys(i)%s))then
                error='Unknown '//group//' option '//keys(i)%s//'.'
                exit
            end if
            if(values%valueType(keys(i)%s)/=QJsonValue_Bool)then
                error=group//' option '//keys(i)%s//' requires a Boolean value.'
                exit
            end if
            toggles(index)=workspace_toggle(.true.,values%boolValue(keys(i)%s))
        end do
        call values%delete()
    end subroutine read_workspace_toggles

    logical function needs_toggle(toggle,current)
        type(workspace_toggle),intent(in) :: toggle
        logical,intent(in) :: current

        needs_toggle=toggle%specified.and.(toggle%value.neqv.current)
    end function needs_toggle

    subroutine describe_workspace(s,result)
        type(studio_t),intent(in) :: s
        type(QJsonObject),intent(inout) :: result
        type(QJsonObject) :: workspace
        type(QTabWidget) :: tabs
        integer :: i, tab, filter, axis

        workspace=QJsonObject()
        call s%workspace%findChild('studioRibbonTabs',tabs)
        if(.not.tabs%isNull())then
            tab=findloc(tab_ids,tabs%currentIndex(),dim=1)
            if(tab>0)call workspace%insert('tab',trim(tab_names(tab)))
        end if
        call workspace%insert('canvas_only',s%canvas_only)
        call workspace%insert('display_style',trim(display_styles(s%viewport_style)))
        call workspace%insert('viewport_tool',trim(viewport_tools(s%surface_tool)))
        call workspace%insert('selection_tool',trim(selection_tools(s%area_selection_tool)))
        call workspace%insert('pick_mode',trim(pick_modes(s%feature_pick_mode)))
        call workspace%insert('snap',trim(snap_modes(current_snap_mode(s))))
        call workspace%insert('domain_tool',trim(domain_tools(current_domain_tool(s))))
        call workspace%insert('knot_groups',s%domain_edit_knots.and.s%domain_knot_details.and.s%panels(5)%isVisible())
        call workspace%insert('scene_search',s%widgets%project_browser%search%text())
        filter=findloc(scene_filter_ids,s%widgets%project_browser%filter_mode,dim=1)
        if(filter>0)call workspace%insert('scene_filter',trim(scene_filters(filter)))
        axis=s%knot_axis
        if(axis>=1.and.axis<=active_rank(s).and.s%selected_multipatch==0.and.s%selected_mesh==0)then
            call workspace%insert('parameter_axis',axes(axis:axis))
        end if
        call insert_workspace_toggles(workspace,'panels',panel_names, &
            [(s%panels(panel_ids(i))%isVisible(),i=1,size(panel_ids)),studio_assistant_is_visible(s)])
        call insert_workspace_toggles(workspace,'inspector',inspector_names, &
            [(s%inspector%headers(i)%isChecked(),i=1,size(inspector_names))])
        call insert_workspace_toggles(workspace,'display',display_names, &
            [s%layers,s%physical_knots,s%world_grid,s%interactive_quality,s%mesh_shaded,s%mesh_only,s%mesh_pick])
        call insert_workspace_toggles(workspace,'domain',domain_names,s%domain_layers)
        call insert_workspace_toggles(workspace,'basis',basis_names,[s%basis_derivative,s%basis_legend,s%basis_probe_guide])
        call result%insert('workspace',workspace)
        call workspace%delete()
    end subroutine describe_workspace

    subroutine insert_workspace_toggles(workspace,group,names,enabled)
        type(QJsonObject),intent(inout) :: workspace
        character(*),intent(in) :: group,names(:)
        logical,intent(in) :: enabled(:)
        type(QJsonObject) :: values
        integer :: i

        values=QJsonObject()
        do i=1,size(names)
            call values%insert(trim(names(i)),enabled(i))
        end do
        call workspace%insert(group,values)
        call values%delete()
    end subroutine insert_workspace_toggles

    subroutine describe_view(s,result)
        type(studio_t),intent(in) :: s
        type(QJsonObject),intent(inout) :: result
        type(QJsonObject) :: view
        type(QJsonArray) :: values
        type(viewport_camera) :: camera

        view=QJsonObject()
        call view%insert('yaw_degrees',modulo((live_yaw_offset+s%yaw)/radians_per_degree,360._dp))
        call view%insert('pitch_degrees',(live_pitch_offset+s%pitch)/radians_per_degree)
        call view%insert('zoom',s%zoom)
        call view%insert('interactive_quality',s%interactive_quality)
        call view%insert('display_style',trim(viewport_style_names(s%viewport_style)))
        call view%insert('pan_x',s%pan_x)
        call view%insert('pan_y',s%pan_y)
        call view%insert('width',s%surface%width())
        call view%insert('height',s%surface%height())
        call view%insert('canvas_only',s%canvas_only)
        call view%insert('mesh_only',s%mesh_only)
        call view%insert('mesh_view',trim(mesh_modes(s%mesh_view)))
        call view%insert('mesh_layer_axis',axes(s%mesh_layer_axis:s%mesh_layer_axis))
        call view%insert('mesh_layer_index',s%mesh_layer_index)
        call view%insert('mesh_section_axis',trim(section_axes(s%mesh_section_axis)))
        call view%insert('mesh_section_percent',nint(100._dp*s%mesh_section_fraction))
        values=QJsonArray()
        call values%appendDoubles(viewport_zoom_range(s))
        call view%insert('zoom_range',values)
        call values%delete()
        values=QJsonArray()
        call values%appendDoubles(workspace_view_bounds(s))
        call view%insert('bounds',values)
        call values%delete()
        if(prepare_camera(s,camera))then
            values=QJsonArray()
            call values%appendDoubles(camera%center)
            call view%insert('model_center',values)
            call values%delete()
            call view%insert('pixels_per_model_unit',camera%scale)
        end if
        call result%insert('camera',view)
        call view%delete()
    end subroutine describe_view

end module studio_mcp_view
