module studio_mcp_schema
    use, intrinsic :: iso_c_binding, only: c_double
    use qt, only: QJsonArray,QJsonDocument,QJsonDocument_fromJson,QJsonObject,QJsonValue,qt_string
    implicit none(type, external)
    private
    public :: studio_mcp_catalog, studio_mcp_tool_schema, studio_mcp_valid_arguments, studio_rpc_token

    character(*), parameter :: lowercase='abcdefghijklmnopqrstuvwxyz', uppercase='ABCDEFGHIJKLMNOPQRSTUVWXYZ', digits='0123456789'

    character(*), parameter :: tool_inspect = &
        '{"name":"inspect","description":"Read the live ForCAD Studio design, object IDs, selection, mesh summaries and '// &
        'diagnostics. Each collection is paged using the same zero-based offset. Use next_offset while has_more is true. '// &
        'Use inspect_assembly for large member lists and select for detailed mesh sizing.",'// &
        '"inputSchema":{"type":"object","properties":{"offset":{"type":"integer","minimum":0,"maximum":2147483647},'// &
        '"limit":{"type":"integer","minimum":1,"maximum":128}},"required":[],"additionalProperties":false},'// &
        '"annotations":{"openWorldHint":false,"readOnlyHint":true}}'

    character(*), parameter :: tool_list_commands = &
        '{"name":"list_commands","description":"List the dialog-free ForCAD Studio commands available in the current '// &
        'design.","inputSchema":{"type":"object","properties":{},"required":[],"additionalProperties":false},'// &
        '"annotations":{"openWorldHint":false,"readOnlyHint":true}}'

    character(*), parameter :: tool_studio_skill = &
        '{"name":"studio_skill","description":"Call with no arguments to discover packaged ForCAD Studio workflows. '// &
        'Pass a listed id for its Markdown guidance on demand. The native tool schemas remain authoritative.",'// &
        '"inputSchema":{"type":"object","properties":{"id":{"type":"string","minLength":1}},'// &
        '"required":[],"additionalProperties":false},"annotations":{"openWorldHint":false,"readOnlyHint":true}}'

    character(*), parameter :: tool_execute = &
        '{"name":"execute","description":"Run a named dialog-free command from list_commands.",'// &
        '"inputSchema":{"type":"object","properties":{"command":{"type":"string","pattern":"^[a-z][a-z0-9_.-]{0,'// &
        '63}$"}},"required":["command"],"additionalProperties":false},"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_select = &
        '{"name":"select","description":"Select a patch, assembly or mesh by kind and ID returned by inspect.",'// &
        '"inputSchema":{"type":"object","properties":{"kind":{"type":"string","enum":["curve","surface","volume",'// &
        '"assembly","mesh"]},"id":{"type":"integer","minimum":1,"maximum":2147483647}},"required":["kind","id"],'// &
        '"additionalProperties":false},"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_set_selection = &
        '{"name":"set_selection","description":"Select, add, or toggle exact objects and subfeatures. IDs come '// &
        'from inspect; feature indexes are one-based. For shared controls use one assembly ID and control indexes from '// &
        'inspect_assembly. A mesh cell also needs its part number. Use an empty items array with replace to clear selection.",'// &
        '"inputSchema":{"type":"object","properties":{"mode":{"type":"s'// &
        'tring","enum":["replace","add","toggle"]},"items":{"type":"array","items":{"type":"object",'// &
        '"properties":{"kind":{"type":"string","enum":["curve","surface","volume","assembly","mesh"]},'// &
        '"id":{"type":"integer","minimum":1,"maximum":2147483647},"feature":{"type":"string","enum":["face","edge",'// &
        '"vertex","control","cell"]},"index":{"type":"integer","minimum":1,"maximum":2147483647},'// &
        '"part":{"type":"integer","minimum":1,"maximum":2147483647}},"required":["kind","id"],'// &
        '"additionalProperties":false}}},"required":["items"],"additionalProperties":false},'// &
        '"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_select_region = &
        '{"name":"select_region","description":"Select from the current viewport using canvas-relative logical XY '// &
        'pixels. Use inspect_view for the drawable bounds. Box takes two diagonal points. Left to right selects '// &
        'enclosed geometry and right to left selects crossing geometry. Lasso takes at least three points and '// &
        'selects touched geometry.","inputSchema":{"type":"object","properties":{"mode":{"type":"string",'// &
        '"enum":["replace","add","toggle"]},"shape":{"type":"string","enum":["box","lasso"]},'// &
        '"points":{"type":"array","minItems":2,"items":{"type":"array","items":{"type":"number"},"minItems":2,'// &
        '"maxItems":2}}},"required":["shape","points"],"additionalProperties":false},'// &
        '"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_inspect_patch = &
        '{"name":"inspect_patch","description":"Read native patch degrees, knot vectors and pages of exact control '// &
        'coordinates and weights. Control and knot offsets are zero-based.","inputSchema":{"type":"object",'// &
        '"properties":{"kind":{"type":"string","enum":["curve","surface","volume"]},"id":{"type":"integer",'// &
        '"minimum":1,"maximum":2147483647},"offset":{"type":"integer","minimum":0,"maximum":2147483647},'// &
        '"limit":{"type":"integer","minimum":1,"maximum":128},"knot_offset":{"type":"integer","minimum":0,'// &
        '"maximum":2147483647},"knot_limit":{"type":"integer","minimum":1,"maximum":128}},"required":["kind","id"],'// &
        '"additionalProperties":false},"annotations":{"openWorldHint":false,"readOnlyHint":true}}'

    character(*), parameter :: tool_inspect_assembly = &
        '{"name":"inspect_assembly","description":"Read an assembly definition, member patch IDs and a page of '// &
        'native interface joins. Pass control_offset or control_limit to read shared control IDs, exact XYZ,'// &
        ' representative weights and copy counts. Offsets are zero-based.","inputSchema":{"type":"object",'// &
        '"properties":{"id":{"type":"integer","minimum":1,"maximum":2147483647},"join_offset":{"type":"integer",'// &
        '"minimum":0,"maximum":2147483647},"join_limit":{"type":"integer","minimum":1,"maximum":128},'// &
        '"control_offset":{"type":"integer","minimum":0,"maximum":2147483647},'// &
        '"control_limit":{"type":"integer","minimum":1,"maximum":128}},'// &
        '"required":["id"],"additionalProperties":false},"annotations":{"openWorldHint":false,"readOnlyHint":true}}'

    character(*), parameter :: tool_inspect_view = &
        '{"name":"inspect_view","description":"Read the exact camera, drawable viewport bounds, ribbon tab, panel '// &
        'visibility, Scene search and display options. Pixel positions are relative to the canvas.",'// &
        '"inputSchema":{"type":"object","properties":{},"required":[],"additionalProperties":false},'// &
        '"annotations":{"openWorldHint":false,"readOnlyHint":true}}'

    character(*), parameter :: tool_set_workspace = &
        '{"name":"set_workspace","description":"Set the ribbon tab, panel visibility, Canvas Only, Scene filter or '// &
        'search, parameter axis or tool, display, pick mode, selection, and snap. Read the current settings with '// &
        'inspect_view.","inputSchema":{"type":"object","properties":{"tab":{"type":"string","enum":["project",'// &
        '"create","model","patch","multipatch","mesh","view"]},"canvas_only":{"type":"boolean"},'// &
        '"scene_search":{"type":"string","minLength":0,"maxLength":1024},"scene_filter":{"type":"string",'// &
        '"enum":["all","geometry","multipatches","meshes","hidden","needs_update","curves","surfaces","volumes"]},'// &
        '"parameter_axis":{"type":"string","enum":["U","V","W"]},'// &
        '"viewport_tool":{"type":"string","enum":["select","orbit"]},"selection_tool":{"type":"string","enum":["click","box","lasso"]},'// &
        '"pick_mode":{"type":"string","enum":["automatic","face","edge","vertex","control","multipatch"]},'// &
        '"snap":{"type":"string","enum":["off","vertex","closest"]},'// &
        '"domain_tool":{"type":"string","enum":["probe","select_knot","insert_knot","remove_knot"]},'// &
        '"display_style":{"type":"string","enum":["shaded","shaded_edges","wireframe","fast_shaded"]},"panels":{"type":"object",'// &
        '"properties":{"scene":{"type":"boolean"},"inspector":{"type":"boolean"},"parameter_space":{"type":"boolean'// &
        '"},"basis_plot":{"type":"boolean"},"assistant":{"type":"boolean"}},"required":[],'// &
        '"additionalProperties":false},"inspector":{"type":"object","properties":{"selection_probe":{"type":"boolea'// &
        'n"},"diagnostics":{"type":"boolean"}},"required":[],"additionalProperties":false},'// &
        '"display":{"type":"object","properties":{"shaded":{"type":"boolean"},"wireframe":{"type":"boolean"},'// &
        '"control_net":{"type":"boolean"},"control_points":{"type":"boolean"},"parameter_lines":{"type":"boolean"},'// &
        '"knot_boundaries":{"type":"boolean"},"ground_grid":{"type":"boolean"},"interactive":{"type":"boolean"},'// &
        '"mesh_shading":{"type":"boolean"},"mesh_only":{"type":"boolean"},"mesh_picking":{"type":"boolean"}},"required":[],'// &
        '"additionalProperties":false},"domain":{"type":"object","properties":{"grid":{"type":"boolean"},'// &
        '"knots":{"type":"boolean"},"elements":{"type":"boolean"},"basis_support":{"type":"boolean"},'// &
        '"selected_element":{"type":"boolean"}},"required":[],"additionalProperties":false},'// &
        '"basis":{"type":"object","properties":{"derivative":{"type":"boolean"},"legend":{"type":"boolean"},'// &
        '"probe":{"type":"boolean"}},"required":[],"additionalProperties":false}},"required":[],'// &
        '"additionalProperties":false,"minProperties":1},"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_set_project_units = &
        '{"name":"set_project_units","description":"Set the model length unit and display length unit, without '// &
        'rescaling the stored geometry. Use same to display in the model unit.","inputSchema":{"type":"object",'// &
        '"properties":{"model":{"type":"string","enum":["unitless","mm","cm","m","in","ft"]},'// &
        '"display":{"type":"string","enum":["same","mm","cm","m","in","ft"]}},"required":[],'// &
        '"additionalProperties":false,"minProperties":1},"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_set_camera = &
        '{"name":"set_camera","description":"Set one or more exact camera fields. Yaw and pitch are absolute '// &
        'degrees, zoom is a native factor, and pan is in canvas pixels. Omitted fields are retained. This does not '// &
        'affect geometry or Undo.","inputSchema":{"type":"object","properties":{"yaw_degrees":{"type":"number"},'// &
        '"pitch_degrees":{"type":"number"},"zoom":{"type":"number","exclusiveMinimum":0},'// &
        '"pan_x":{"type":"number"},"pan_y":{"type":"number"}},"required":[],'// &
        '"additionalProperties":false,"minProperties":1},"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_set_parameter_view = &
        '{"name":"set_parameter_view","description":"Set exact Parameter Space zoom and pan for a single patch, '// &
        'and absolute yaw and pitch in degrees for a volume. Get the native ranges and current camera from '// &
        'inspect_view.","inputSchema":{"type":"object","properties":{"zoom":{"type":"number","minimum":0.65,'// &
        '"maximum":24},"pan_x":{"type":"number"},"pan_y":{"type":"number"},"yaw_degrees":{"type":"number"},'// &
        '"pitch_degrees":{"type":"number"}},"required":[],"additionalProperties":false,"minProperties":1},'// &
        '"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_set_mesh_view = &
        '{"name":"set_mesh_view","description":"Change the visible volume mesh presentation or its section. Use the U/V/W '// &
        'axis and one-based index for a layer, or a world X/Y/Z cut with an integer position percentage from 0 to 100. '// &
        'Set section_axis to off to remove the cut. This changes the view only.",'// &
        '"inputSchema":{"type":"object","properties":{"mode":{"type":"string","enum":["exterior","lattice","layer"]},'// &
        '"axis":{"type":"string","enum":["U","V","W"]},"layer":{"type":"integer","minimum":1,"maximum":2147483647},'// &
        '"section_axis":{"type":"string","enum":["off","X","Y","Z"]},"section_percent":{"type":"integer","minimum":0,"maximum":100}},'// &
        '"required":[],"minProperties":1,"additionalProperties":false},"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_measure = &
        '{"name":"measure","description":"Measure geometry using ForCAD in project display units. Omit arguments '// &
        'to measure the current selection. Multiple selections return measurements grouped by dimension. Supply '// &
        'both kind and ID to measure an exact patch or assembly.","inputSchema":{"type":"object",'// &
        '"properties":{"kind":{"type":"string","enum":["curve","surface","volume","assembly"]},'// &
        '"id":{"type":"integer","minimum":1,"maximum":2147483647}},"required":[],"additionalProperties":false},'// &
        '"annotations":{"openWorldHint":false,"readOnlyHint":true}}'

    character(*), parameter :: tool_transform_scene = &
        '{"name":"transform_scene","description":"Apply ordered rotations in degrees about the world origin and '// &
        'world-axis moves to all scene patches. Undo reverses one transform per call. A move without a unit uses '// &
        'project coordinates.","inputSchema":{"type":"object","properties":{"steps":{"type":"array","minItems":1,'// &
        '"items":{"type":"object","properties":{"op":{"type":"string","enum":["rotate","move"]},'// &
        '"axis":{"type":"string","enum":["X","Y","Z"]},"value":{"type":"number"},"unit":{"type":"string",'// &
        '"enum":["mm","cm","m","in","ft"]}},"required":["op","axis","value"],"additionalProperties":false}}},'// &
        '"required":["steps"],"additionalProperties":false},"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_build_model = &
        '{"name":"build_model","description":"Build an ordered sequence of ForCAD primitives and construction '// &
        'operations. References are one-based indexes of earlier steps. Inspect the design for project units.",'// &
        '"inputSchema":{"type":"object","properties":{"steps":{"type":"array","items":{"type":"object",'// &
        '"properties":{"op":{"type":"string","enum":["point","line","circle","arc","rectangle","box","cylinder",'// &
        '"tube","extrude","revolve","sweep","loft"]},"name":{"type":"string","minLength":1,"maxLength":64},'// &
        '"origin":{"type":"array","items":{"type":"number"},"minItems":3,"maxItems":3},"values":{"type":"array",'// &
        '"items":{"type":"number"},"minItems":3,"maxItems":3},"source":{"type":"integer","minimum":1,"maximum":12},'// &
        '"spine":{"type":"integer","minimum":1,"maximum":12},"angle":{"type":"number","minimum":-360,'// &
        '"maximum":360},"degree":{"type":"integer","minimum":0,"maximum":11},"sections":{"type":"array",'// &
        '"items":{"type":"integer","minimum":1,"maximum":12},"minItems":2,"maxItems":12}},"required":["op","name"],'// &
        '"additionalProperties":false},"minItems":1,"maxItems":12}},"required":["steps"],'// &
        '"additionalProperties":false},"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_build_curved_volume = &
        '{"name":"build_curved_volume","description":"Build an exact round NURBS multipatch solid from 2 to 12 '// &
        'user-defined axial levels around the world Z axis. Use strictly ascending z and positive outer radii in '// &
        'project units. Omit control_radius and sectors on the first level. Each later level describes the span '// &
        'ending there. The optional control_radius sets its middle quadratic control radius, defaulting to the '// &
        'linear midpoint. Omit sectors for a full disk span, or select distinct 45-degree annular sectors 1 to 8 '// &
        'counterclockwise from +X for partial spans. Partial spans must connect through matching full faces. '// &
        'mesh_cells optionally builds one conforming hex mesh with the specified cells per patch direction. The '// &
        'volume and mesh are one undoable edit.","inputSchema":{"type":"object","properties":{"name":{"type":"strin'// &
        'g","minLength":1,"maxLength":96},"levels":{"type":"array","items":{"type":"object",'// &
        '"properties":{"z":{"type":"number"},"radius":{"type":"number","exclusiveMinimum":0},'// &
        '"control_radius":{"type":"number","exclusiveMinimum":0},"sectors":{"type":"array",'// &
        '"items":{"type":"integer","minimum":1,"maximum":8},"minItems":1,"maxItems":8}},"required":["z","radius"],'// &
        '"additionalProperties":false},"minItems":2,"maxItems":12},"mesh_cells":{"type":"array",'// &
        '"items":{"type":"integer","minimum":1,"maximum":2147483646},"minItems":3,"maxItems":3}},'// &
        '"required":["name","levels"],"additionalProperties":false},"annotations":{"openWorldHint":false}}'

    character(*), parameter :: control_target_schema = &
        '"target":{"type":"object","properties":{"kind":{"type":"string","enum":["curve","surface","volume","assembly","multipatch"]},'// &
        '"id":{"type":"integer","minimum":1,"maximum":2147483647}},"required":["kind","id"],"additionalProperties":false}'

    character(*), parameter :: control_fields = &
        '"point":{"type":"integer","minimum":1,"maximum":2147483647},'// &
        '"position":{"type":"array","items":{"type":"number"},"minItems":3,"maxItems":3},'// &
        '"weight":{"type":"number","exclusiveMinimum":0}'

    character(*), parameter :: controls_array_schema = &
        '"controls":{"type":"array","minItems":1,"items":{"type":"object","properties":{'//control_fields// &
        '},"required":["point"],"additionalProperties":false}}'

    character(*), parameter :: tool_preview_controls = &
        '{"name":"preview_controls","description":"Preview a control edit with the same inputs as edit_patch control. '// &
        'Shows exact dependent XYZ and weights before and after preserving smooth joins, without editing. Use point and position '// &
        'or weight, or pass several edits in controls. Positions are in model units. Assembly point IDs come from inspect_assembly. '// &
        'Member point IDs come from inspect_patch. Use next_offset to page through dependent controls.",'// &
        '"inputSchema":{"type":"object","properties":{"op":{"type":"string","enum":["control"]},'//control_target_schema//','// &
        control_fields//','//controls_array_schema//','// &
        '"offset":{"type":"integer","minimum":0,"maximum":2147483647},"limit":{"type":"integer","minimum":1,"maximum":128}},'// &
        '"additionalProperties":false},"annotations":{"openWorldHint":false,"readOnlyHint":true}}'

    character(*), parameter :: tool_edit_patch = &
        '{"name":"edit_patch","description":"Edit a patch. Use target with an exact kind and ID from inspect so '// &
        'the edit follows the intended patch if the user selects another object. Omit target to use the current '// &
        'selection. Control edits use a one-based point ID and exact position and/or weight. Use a global point '// &
        'from inspect_assembly with kind assembly or multipatch. Linked member IDs use patch-local numbering. '// &
        'Pass controls for several point and weight edits in one Undo. Call preview_controls before applying them. '// &
        'On an assembly, insert_knot or elevate_degree take the starting member patch ID in member, its axis and optional count. '// &
        'Call preview_refinement with those same arguments before an assembly refinement. '// &
        'Parameter edits use U, V or W. Move a knot group with value and to, or use set_knot_values with the full values array. Extract '// &
        'an exact section from a surface or volume using extract_section, axis and value. Extract a volume '// &
        'boundary using extract_face and face 1–6 in U min, U max, V min, V max, W min, W max order.",'// &
        '"inputSchema":{"type":"object","properties":{"op":{"type":"string","enum":["control","insert_knot",'// &
        '"remove_knot","move_knot","set_knot_values","elevate_degree","reduce_degree","refine","set_interval",'// &
        '"reverse","split","extract","extract_section","extract_face"]},'//control_target_schema//','// &
        '"member":{"type":"integer","minimum":1,"maximum":2147483647},"axis":{"type":"string","enum":["U","V","W","all"]},'// &
        control_fields//','//controls_array_schema//',"value":{"type":"number"},"to":{"type":"number"},'// &
        '"values":{"type":"array","minItems":1,"items":{"type":"number"}},"count":{"type":"integer","minimum":1,'// &
        '"maximum":2147483647},"tolerance":{"type":"number","minimum":0},"divisions":{"type":"integer","minimum":1,'// &
        '"maximum":2147483647},"degree_steps":{"type":"integer","minimum":0,"maximum":2147483647},'// &
        '"lower":{"type":"number"},"upper":{"type":"number"},"face":{"type":"integer","minimum":1,"maximum":6}},'// &
        '"required":["op"],"additionalProperties":false},"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_assemble = &
        '{"name":"assemble","description":"Join native patches using member IDs and explicit interfaces. Supply an '// &
        'existing assembly ID to replace its definition with one Undo. Pass an empty joins array to save an '// &
        'unjoined assembly, or use continuity -1 to record an unconstrained interface. For geometric continuity '// &
        'of order n greater than zero, set geometric to true and provide n finite normal transition derivatives in '// &
        'reparameterization. Use inspect_assembly and inspect_patch to read existing definitions.",'// &
        '"inputSchema":{"type":"object","properties":{"id":{"type":"integer","minimum":1,"maximum":2147483647},'// &
        '"name":{"type":"string","minLength":1},"kind":{"type":"string","enum":["curve","surface",'// &
        '"volume"]},"members":{"type":"array","items":{"type":"integer","minimum":1,"maximum":2147483647},'// &
        '"minItems":2},"joins":{"type":"array","items":{"type":"object",'// &
        '"properties":{"patch_a":{"type":"integer","minimum":1,"maximum":2147483647},"side_a":{"type":"integer",'// &
        '"minimum":1,"maximum":2147483647},"patch_b":{"type":"integer","minimum":1,"maximum":2147483647},'// &
        '"side_b":{"type":"integer","minimum":1,"maximum":2147483647},"continuity":{"type":"integer","minimum":-1,'// &
        '"maximum":2147483647},"reverse":{"type":"boolean"},"swap":{"type":"boolean"},"flip":{"type":"array",'// &
        '"items":{"type":"boolean"},"minItems":2,"maxItems":2},"geometric":{"type":"boolean"},'// &
        '"reparameterization":{"type":"array","items":{"type":"number"},"minItems":1}},'// &
        '"required":["patch_a","side_a","patch_b","side_b"],"additionalProperties":false}}},'// &
        '"required":["name","kind","members","joins"],"additionalProperties":false},'// &
        '"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_preview_refinement = &
        '{"name":"preview_refinement","description":"Validate an assembly knot insertion or degree elevation without editing. '// &
        'Returns affected member IDs, mapped parameter directions and exact knot values. Use the same op, target, member, '// &
        'axis, count and value as edit_patch. Supply value for insert_knot. Omit target to use the selected multipatch. '// &
        'Use next_offset for additional rows.",'// &
        '"inputSchema":{"type":"object","properties":{"op":{"type":"string","enum":["insert_knot","elevate_degree"]},'// &
        '"target":{"type":"object","properties":{"kind":{"type":"string","enum":["assembly","multipatch"]},'// &
        '"id":{"type":"integer","minimum":1,"maximum":2147483647}},"required":["kind","id"],"additionalProperties":false},'// &
        '"member":{"type":"integer","minimum":1,"maximum":2147483647},"axis":{"type":"string","enum":["U","V","W"]},'// &
        '"count":{"type":"integer","minimum":1,"maximum":2147483647},"value":{"type":"number"},'// &
        '"offset":{"type":"integer","minimum":0,"maximum":2147483647},"limit":{"type":"integer","minimum":1,"maximum":128}},'// &
        '"required":["op","member","axis"],"additionalProperties":false},'// &
        '"annotations":{"openWorldHint":false,"readOnlyHint":true}}'

    character(*), parameter :: tool_preview_mesh = &
        '{"name":"preview_mesh","description":"Estimate per-patch cell counts from a positive target size in model '// &
        'units using the same sizing routine as mesh. Optionally synchronize interface counts. This does not '// &
        'generate or change any mesh.","inputSchema":{"type":"object","properties":{"target":{"type":"object",'// &
        '"properties":{"kind":{"type":"string","enum":["curve","surface","volume","assembly","mesh"]},'// &
        '"id":{"type":"integer","minimum":1,"maximum":2147483647}},"required":["kind","id"],'// &
        '"additionalProperties":false},"target_size":{"type":"number","exclusiveMinimum":0},'// &
        '"synchronize":{"type":"boolean"}},"required":["target_size"],"additionalProperties":false},'// &
        '"annotations":{"openWorldHint":false,"readOnlyHint":true}}'

    character(*), parameter :: tool_mesh = &
        '{"name":"mesh","description":"Generate or replace a mesh. Use target with an exact patch, assembly, or '// &
        'mesh kind and ID from inspect, or omit target to use the current selection. Use cells for uniform sizing '// &
        'or supply patch_cells with one row per assembly member. Counts follow U, V, W order. Alternatively pass a '// &
        'positive target_size in model units and optionally synchronize interface counts. Set conforming to stitch '// &
        'a volume assembly into one hex mesh.","inputSchema":{"type":"object","properties":{"target":{"type":"objec'// &
        't","properties":{"kind":{"type":"string","enum":["curve","surface","volume","assembly","mesh"]},'// &
        '"id":{"type":"integer","minimum":1,"maximum":2147483647}},"required":["kind","id"],'// &
        '"additionalProperties":false},"cells":{"type":"array","items":{"type":"integer","minimum":1,'// &
        '"maximum":2147483646},"minItems":1,"maxItems":3},"patch_cells":{"type":"array","minItems":1,'// &
        '"items":{"type":"object","properties":{"patch_id":{"type":"integer","minimum":1,"maximum":2147483647},'// &
        '"cells":{"type":"array","items":{"type":"integer","minimum":1,"maximum":2147483646},"minItems":1,'// &
        '"maxItems":3}},"required":["patch_id","cells"],"additionalProperties":false}},'// &
        '"target_size":{"type":"number","exclusiveMinimum":0},"synchronize":{"type":"boolean"},'// &
        '"conforming":{"type":"boolean"}},"required":[],"additionalProperties":false,'// &
        '"oneOf":[{"required":["cells"]},{"required":["patch_cells"]},{"required":["target_size"]}]},'// &
        '"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_transform_selection = &
        '{"name":"transform_selection","description":"Apply an exact Move, Rotate, Scale or Affine transform to '// &
        'the selected patch or validated multipatch assembly. One Undo restores the selection. For move supply '// &
        'offset XYZ. For rotate supply angles in degrees applied about world X, Y and Z in that order. For scale '// &
        'supply nonzero XYZ factors. For affine supply a 3×3 matrix by row and optionally an offset. Rotate, scale '// &
        'and affine accept a pivot; it defaults to the world origin. Offset and pivot use project units unless an '// &
        'explicit unit is supplied.","inputSchema":{"type":"object","properties":{"op":{"type":"string",'// &
        '"enum":["move","rotate","scale","affine"]},"offset":{"type":"array","items":{"type":"number"},'// &
        '"minItems":3,"maxItems":3},"angles":{"type":"array","items":{"type":"number"},"minItems":3,"maxItems":3},'// &
        '"factors":{"type":"array","items":{"type":"number"},"minItems":3,"maxItems":3},"matrix":{"type":"array",'// &
        '"items":{"type":"array","items":{"type":"number"},"minItems":3,"maxItems":3},"minItems":3,"maxItems":3},'// &
        '"pivot":{"type":"array","items":{"type":"number"},"minItems":3,"maxItems":3},"unit":{"type":"string",'// &
        '"enum":["mm","cm","m","in","ft"]}},"required":["op"],"additionalProperties":false},'// &
        '"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_construct = &
        '{"name":"construct","description":"Extrude, revolve, sweep or loft existing ForCAD Studio geometry. Use native '// &
        'patch IDs from inspect. Extrude and revolve take one profile and a vector. Revolve also takes axis_point '// &
        'and angle_degrees. Sweep takes one curve or surface profile and a spine curve ID. Loft uses two or more '// &
        'profiles of the given kind. Omit degree or set it to zero for automatic loft degree. One Undo removes the result.",'// &
        '"inputSchema":{"type":"object","properties":{"op":{"type":"string","enum":["extrude","revolve","sweep",'// &
        '"loft"]},"name":{"type":"string","minLength":1},"kind":{"type":"string","enum":["curve",'// &
        '"surface"]},"profiles":{"type":"array","items":{"type":"integer","minimum":1,"maximum":2147483647},'// &
        '"minItems":1},"vector":{"type":"array","items":{"type":"number"},"minItems":3,"maxItems":3},'// &
        '"axis_point":{"type":"array","items":{"type":"number"},"minItems":3,"maxItems":3},'// &
        '"angle_degrees":{"type":"number","minimum":-360,"maximum":360},"spine":{"type":"integer","minimum":1,'// &
        '"maximum":2147483647},"degree":{"type":"integer","minimum":0,"maximum":2147483647}},"required":["op","name",'// &
        '"kind","profiles"],"additionalProperties":false},"annotations":{"openWorldHint":false}}'

    character(*), parameter :: deletion_target_properties = &
        '"kind":{"type":"string","enum":["curve","surface","volume","assembly","mesh"]},'// &
        '"id":{"type":"integer","minimum":1,"maximum":2147483647}'

    character(*), parameter :: tool_delete = &
        '{"name":"delete","description":"Delete whole patches, assemblies or meshes. Use a kind and ID from inspect, '// &
        'an array of exact targets, or selection true to delete the current selection. Batch deletion removes affected meshes '// &
        'and is one Undo step. To delete an assembly with its geometry, get its member IDs from inspect_assembly and include '// &
        'both the assembly and those patches in targets. Include every assembly that contains a patch being deleted.",'// &
        '"inputSchema":{"type":"object","properties":{'//deletion_target_properties//',"selection":{"type":"boolean"},'// &
        '"targets":{"type":"array","minItems":1,"items":{"type":"object","properties":{'//deletion_target_properties// &
        '},"required":["kind","id"],"additionalProperties":false}}},"oneOf":[{"required":["kind","id"]},'// &
        '{"required":["targets"]},{"required":["selection"]}],"additionalProperties":false},'// &
        '"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_set_material = &
        '{"name":"set_material","description":"Set the selected patch or all members of the selected assembly to '// &
        'exact RGB color or opacity. Use integers 0–255. Zero opacity is transparent and 255 is opaque. The edit '// &
        'can be undone.","inputSchema":{"type":"object","properties":{"rgb":{"type":"array",'// &
        '"items":{"type":"integer","minimum":0,"maximum":255},"minItems":3,"maxItems":3},'// &
        '"opacity":{"type":"integer","minimum":0,"maximum":255}},"required":[],"additionalProperties":false},'// &
        '"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_set_visibility = &
        '{"name":"set_visibility","description":"Show or hide a patch, all members of an assembly, or a current '// &
        'mesh by kind and ID. Keeps the current selection.","inputSchema":{"type":"object",'// &
        '"properties":{"kind":{"type":"string","enum":["curve","surface","volume","assembly","mesh"]},'// &
        '"id":{"type":"integer","minimum":1,"maximum":2147483647},"visible":{"type":"boolean"}},"required":["kind",'// &
        '"id","visible"],"additionalProperties":false},"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_copy_patch = &
        '{"name":"copy_patch","description":"Copy a native curve, surface or volume by ID, with a new name and an '// &
        'optional XYZ offset in model units. Selects the new patch and records one Undo.",'// &
        '"inputSchema":{"type":"object","properties":{"kind":{"type":"string","enum":["curve","surface","volume"]},'// &
        '"id":{"type":"integer","minimum":1,"maximum":2147483647},"name":{"type":"string","minLength":1},'// &
        '"offset":{"type":"array","items":{"type":"number"},"minItems":3,"maxItems":3}},"required":["kind","id",'// &
        '"name"],"additionalProperties":false},"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_rename_patch = &
        '{"name":"rename_patch","description":"Rename a curve, surface or volume patch by ID. The name edit can be '// &
        'undone.","inputSchema":{"type":"object","properties":{"kind":{"type":"string","enum":["curve","surface",'// &
        '"volume"]},"id":{"type":"integer","minimum":1,"maximum":2147483647},"name":{"type":"string","minLength":1,'// &
        '"maxLength":96}},"required":["kind","id","name"],"additionalProperties":false},'// &
        '"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_set_probe = &
        '{"name":"set_probe","description":"Move the selected patch parameter probe and return the evaluated '// &
        'physical position. Supply one value per active direction in U, V, W order. Values are in patch parameter '// &
        'coordinates unless normalized is true.","inputSchema":{"type":"object","properties":{"parameters":{"type":'// &
        '"array","items":{"type":"number"},"minItems":1,"maxItems":3},"normalized":{"type":"boolean"}},'// &
        '"required":["parameters"],"additionalProperties":false},"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_pick_mesh_cell = &
        '{"name":"pick_mesh_cell","description":"Highlight a cell of a current visible mesh by its one-based cell '// &
        'number. Returns the cell vertices in project units and native mesh order.","inputSchema":{"type":"object",'// &
        '"properties":{"id":{"type":"integer","minimum":1,"maximum":2147483647},"cell":{"type":"integer",'// &
        '"minimum":1,"maximum":2147483647}},"required":["id","cell"],"additionalProperties":false},'// &
        '"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_import_file = &
        '{"name":"import_file","description":"Import IGES (.igs or .iges) or merge an .fcad project into the open '// &
        'document using a project-relative path. For IGES only, optionally supply source_unit or use keep to preserve '// &
        'numeric coordinates. One Undo removes the import.","inputSchema":{"type":"object",'// &
        '"properties":{"path":{"type":"string","minLength":1,"maxLength":1024},"source_unit":{"type":"string",'// &
        '"enum":["mm","cm","m","in","ft","keep"]}},"required":["path"],"additionalProperties":false},'// &
        '"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_export_file = &
        '{"name":"export_file","description":"Export to a new project-relative file. Use .igs or .iges for a '// &
        'selected standalone curve or surface, .vtk for a selected geometry, assembly or mesh, .fcad for a '// &
        'project snapshot, .f90 for a standalone ForCAD rebuild program, or .forcad/screenshots/name.png '// &
        'to capture the 3D viewport without panels or chat. '// &
        'Existing files are preserved.",'// &
        '"inputSchema":{"type":"object","properties":{"path":{"type":"string","minLength":1,"maxLength":1024}},'// &
        '"required":["path"],"additionalProperties":false},"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_new_project = &
        '{"name":"new_project","description":"Start a blank ForCAD Studio project in the assigned folder. Set new_tab true '// &
        'to keep the current project in its own tab. Save changes before replacing the current project.",'// &
        '"inputSchema":{"type":"object","properties":{"new_tab":{"type":"boolean"}},"required":[],'// &
        '"additionalProperties":false},"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_save_project = &
        '{"name":"save_project","description":"Save the active ForCAD Studio project in its assigned folder. Omit path to '// &
        'save the current project file. Use a new relative .fcad path for Save As. Existing unrelated files are '// &
        'preserved.","inputSchema":{"type":"object","properties":{"path":{"type":"string","minLength":1,'// &
        '"maxLength":1024}},"required":[],"additionalProperties":false},"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_open_project = &
        '{"name":"open_project","description":"Open an existing .fcad from the assigned folder using a relative path. '// &
        'Set new_tab true to keep the current project in its own tab. Save changes before replacing it. Inspect '// &
        'again for current object IDs.","inputSchema":{"type":"object","properties":{'// &
        '"path":{"type":"string","minLength":1,"maxLength":1024},"new_tab":{"type":"boolean"}},'// &
        '"required":["path"],"additionalProperties":false},"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_list_project_tabs = &
        '{"name":"list_project_tabs","description":"List open ForCAD Studio tabs in the assigned project folder, '// &
        'including which is active and which has unsaved changes. Use the returned tab and revision to switch or close.",'// &
        '"inputSchema":{"type":"object","properties":{},"additionalProperties":false},'// &
        '"annotations":{"openWorldHint":false,"readOnlyHint":true}}'

    character(*), parameter :: tab_target_schema = &
        '"inputSchema":{"type":"object","properties":{"tab":{"type":"integer","minimum":1,"maximum":2147483647},'// &
        '"revision":{"type":"string","minLength":1}},"required":["tab","revision"],"additionalProperties":false}'

    character(*), parameter :: tool_switch_project_tab = &
        '{"name":"switch_project_tab","description":"Activate an open project tab within the assigned folder. '// &
        'Use the tab and revision from list_project_tabs. Existing edits stay in their tabs.",'// &
        tab_target_schema//',"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_close_project_tab = &
        '{"name":"close_project_tab","description":"Close an inactive, saved project tab within the assigned folder. '// &
        'Use the tab and revision from list_project_tabs.",'// &
        tab_target_schema//',"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_undo = &
        '{"name":"undo","description":"Undo the most recent geometry, assembly, mesh or appearance edit in '// &
        'ForCAD Studio.","inputSchema":{"type":"object","properties":{},"required":[],"additionalProperties":false},'// &
        '"annotations":{"openWorldHint":false}}'

    character(*), parameter :: tool_redo = &
        '{"name":"redo","description":"Restore the latest undone edit in ForCAD Studio.","inputSchema":{"type":"object",'// &
        '"properties":{},"required":[],"additionalProperties":false},"annotations":{"openWorldHint":false}}'

    character(*), parameter :: catalog_json = '['// &
        tool_inspect//','// &
        tool_list_commands//','// &
        tool_studio_skill//','// &
        tool_execute//','// &
        tool_select//','// &
        tool_set_selection//','// &
        tool_select_region//','// &
        tool_inspect_patch//','// &
        tool_inspect_assembly//','// &
        tool_inspect_view//','// &
        tool_set_workspace//','// &
        tool_set_project_units//','// &
        tool_set_camera//','// &
        tool_set_parameter_view//','// &
        tool_set_mesh_view//','// &
        tool_measure//','// &
        tool_transform_scene//','// &
        tool_build_model//','// &
        tool_build_curved_volume//','// &
        tool_edit_patch//','// &
        tool_preview_controls//','// &
        tool_preview_refinement//','// &
        tool_assemble//','// &
        tool_preview_mesh//','// &
        tool_mesh//','// &
        tool_transform_selection//','// &
        tool_construct//','// &
        tool_delete//','// &
        tool_set_material//','// &
        tool_set_visibility//','// &
        tool_copy_patch//','// &
        tool_rename_patch//','// &
        tool_set_probe//','// &
        tool_pick_mesh_cell//','// &
        tool_import_file//','// &
        tool_export_file//','// &
        tool_new_project//','// &
        tool_save_project//','// &
        tool_open_project//','// &
        tool_list_project_tabs//','// &
        tool_switch_project_tab//','// &
        tool_close_project_tab//','// &
        tool_undo//','// &
        tool_redo//']'

contains

    function studio_mcp_catalog() result(catalog)
        type(QJsonArray) :: catalog
        type(QJsonDocument) :: document

        document=QJsonDocument_fromJson(catalog_json)
        catalog=document%array()
        call document%delete()
    end function studio_mcp_catalog

    function studio_mcp_tool_schema(catalog,name) result(schema)
        type(QJsonArray), intent(in) :: catalog
        character(*), intent(in) :: name
        type(QJsonObject) :: schema,entry
        integer :: i

        do i=0,catalog%size()-1
            entry=catalog%objectAt(i)
            if(entry%stringValue('name')==name)then
                schema=entry%objectValue('inputSchema')
                call entry%delete()
                return
            end if
            call entry%delete()
        end do
        schema=QJsonObject()
    end function studio_mcp_tool_schema

    logical function studio_rpc_token(text,tool) result(valid)
        character(*), intent(in) :: text
        logical, intent(in) :: tool

        valid=.false.
        if(len(text)==0)return
        if(tool)then
            valid=len(text)<=32.and.verify(text,lowercase//'_')==0
        else
            valid=len(text)<=128.and.verify(text(:1),lowercase//uppercase)==0.and. &
                verify(text(2:),lowercase//uppercase//digits//'_/.-')==0
        end if
    end function studio_rpc_token

    logical function studio_mcp_valid_arguments(arguments,schema) result(valid)
        type(QJsonValue), intent(in) :: arguments
        type(QJsonObject), intent(in) :: schema

        valid=matches_schema(arguments,schema)
    end function studio_mcp_valid_arguments

    recursive logical function matches_schema(value,schema) result(valid)
        type(QJsonValue), intent(in) :: value
        type(QJsonObject), intent(in) :: schema
        type(QJsonObject) :: object,child_schema
        type(QJsonArray) :: array
        type(QJsonValue) :: child
        character(:), allocatable :: kind
        real(c_double) :: number
        integer :: i

        valid=.false.
        kind=schema%stringValue('type')
        select case(kind)
        case('object')
            if(.not.value%isObject())return
            object=value%toObject()
            valid=valid_object(object,schema)
            call object%delete()
        case('array')
            if(.not.value%isArray())return
            array=value%toArray()
            valid=array%size()>=schema%intValue('minItems',0).and.array%size()<=schema%intValue('maxItems',huge(0))
            if(valid)then
                child_schema=schema%objectValue('items')
                do i=0,array%size()-1
                    child=array%at(i)
                    valid=matches_schema(child,child_schema)
                    call child%delete()
                    if(.not.valid)exit
                end do
                call child_schema%delete()
            end if
            call array%delete()
        case('number','integer')
            if(.not.value%isDouble())return
            number=value%toDouble()
            if(kind=='integer'.and.number/=aint(number))return
            valid=number>=schema%doubleValue('minimum',-huge(number)).and.number<=schema%doubleValue('maximum',huge(number))
            if(schema%contains('exclusiveMinimum'))valid=valid.and.number>schema%doubleValue('exclusiveMinimum')
        case('string')
            if(value%isString())valid=valid_string(value%toString(),schema)
        case('boolean')
            valid=value%isBool()
        end select
    end function matches_schema

    recursive logical function valid_object(object,schema) result(valid)
        type(QJsonObject), intent(in) :: object,schema
        type(QJsonArray) :: alternatives
        type(QJsonObject) :: candidate,properties,child_schema
        type(QJsonValue) :: child
        type(qt_string), allocatable :: keys(:)
        character(:), allocatable :: key
        integer :: i,matching

        valid=.false.
        if(object%size()<schema%intValue('minProperties',0))return
        if(.not.has_required_keys(object,schema))return
        if(schema%contains('oneOf'))then
            alternatives=schema%arrayValue('oneOf')
            matching=0
            do i=0,alternatives%size()-1
                candidate=alternatives%objectAt(i)
                if(has_required_keys(object,candidate))matching=matching+1
                call candidate%delete()
            end do
            call alternatives%delete()
            if(matching/=1)return
        end if
        keys=object%keys()
        properties=schema%objectValue('properties')
        valid=.true.
        do i=1,size(keys)
            key=keys(i)%s
            if(.not.properties%contains(key))then
                valid=.false.
                exit
            end if
            child=object%value(key)
            child_schema=properties%objectValue(key)
            valid=matches_schema(child,child_schema)
            call child%delete()
            call child_schema%delete()
            if(.not.valid)exit
        end do
        call properties%delete()
    end function valid_object

    logical function has_required_keys(object,schema) result(valid)
        type(QJsonObject), intent(in) :: object,schema
        type(QJsonArray) :: fields
        integer :: i

        fields=schema%arrayValue('required')
        valid=.true.
        do i=0,fields%size()-1
            if(object%contains(fields%stringAt(i)))cycle
            valid=.false.
            exit
        end do
        call fields%delete()
    end function has_required_keys

    logical function valid_string(text,schema) result(valid)
        character(*), intent(in) :: text
        type(QJsonObject), intent(in) :: schema
        type(QJsonArray) :: choices
        integer :: i,b,n

        valid=.false.
        n=0
        do i=1,len(text)
            b=iachar(text(i:i))
            if(b<32.or.b==127)return
            if(iand(b,192)/=128)n=n+1
        end do
        if(schema%contains('pattern'))then
            if(schema%stringValue('pattern')/='^[a-z][a-z0-9_.-]{0,63}$')return
            if(n<1.or.n>64)return
            valid=verify(text(:1),lowercase)==0.and.verify(text(2:),lowercase//digits//'_.-')==0
        else if(schema%contains('enum'))then
            choices=schema%arrayValue('enum')
            do i=0,choices%size()-1
                if(text/=choices%stringAt(i))cycle
                valid=.true.
                exit
            end do
            call choices%delete()
        else
            valid=n>=schema%intValue('minLength',0).and.n<=schema%intValue('maxLength',huge(0))
        end if
    end function valid_string

end module studio_mcp_schema
