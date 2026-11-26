module studio_forcad
    use forcad, only: rk
    use studio_forcad_samples, only: make_sample_patch
    use studio_forcad_transforms, only: affine_transform_t, studio_native_patch, native_patch_kind, rotation_matrix, &
        transform_patch, transform_assembly
    use studio_forcad_queries, only: closest_point_on_patch, physical_knot_grid, sample_physical_knots
    use studio_forcad_measure, only: measure_native_patch
    use studio_forcad_evaluation, only: evaluate_patch, evaluate_patch_basis, sample_patch_basis, check_patch_jacobian, studio_native_problem
    use studio_forcad_construction, only: construction_names, construction_actions, stage_studio_construction, &
        extrude_operation, revolve_operation, sweep_operation, loft_operation
    use studio_forcad_primitives, only: stage_studio_primitive, primitive_names, primitive_actions, &
        point_shape, line_shape, circle_shape, arc_shape, rectangle_shape, box_shape, cylinder_shape, tube_shape
    use studio_forcad_parameter_parts, only: parameter_parts, extract_native_slice
    use studio_forcad_patch_data, only: native_knot_axis, describe_patch, copy_native_patch, &
        get_native_controls, sample_patch_geometry, sample_volume_slice, sample_display_grid
    use studio_forcad_patch_mutations, only: replace_native_controls, replace_native_knot_values, edit_native_knots, &
        insert_native_knots, edit_native_degree, reparameterize_patch, reverse_patch, translate_patch
    use studio_forcad_assembly, only: multipatch_join, multipatch_side_label, initialize_assembly, &
        append_assembly_patch, append_assembly_join, inspect_native_assembly, native_assembly_control_map, &
        native_assembly_control_correction, multipatch_join_directions
    use studio_forcad_mesh, only: forcad_mesh_samples, sample_forcad_patch, sample_forcad_assembly, &
        volume_assembly_size, check_forcad_mesh_joins, join_forcad_volume_mesh, estimate_forcad_mesh_samples
    use studio_forcad_exchange, only: write_iges_worker_output, read_iges_worker_output, write_native_iges, &
        write_iges_units, write_native_vtk, write_native_multipatch_vtk, write_mesh_vtk, native_multipatch_vtk_path_limit
    use studio_forcad_patch_codec, only: native_patch_data, collect_native_patch_data, restore_native_patch_data
    use studio_forcad_volume_builders, only: stage_box_volume_patches, stage_curved_volume_patches, &
        curved_core_count, curved_volume_sectors
    implicit none(type, external)
    private

    public :: rk, make_sample_patch
    public :: affine_transform_t, studio_native_patch, native_patch_kind, rotation_matrix, &
        transform_patch, transform_assembly
    public :: closest_point_on_patch, physical_knot_grid, sample_physical_knots
    public :: measure_native_patch
    public :: evaluate_patch, evaluate_patch_basis, sample_patch_basis, check_patch_jacobian, studio_native_problem
    public :: construction_names, construction_actions, stage_studio_construction
    public :: extrude_operation, revolve_operation, sweep_operation, loft_operation
    public :: stage_studio_primitive, primitive_names, primitive_actions
    public :: point_shape, line_shape, circle_shape, arc_shape, rectangle_shape, box_shape, cylinder_shape, tube_shape
    public :: parameter_parts, extract_native_slice
    public :: native_knot_axis, describe_patch, copy_native_patch, get_native_controls, sample_patch_geometry, sample_volume_slice, sample_display_grid
    public :: replace_native_controls, replace_native_knot_values, edit_native_knots, insert_native_knots, edit_native_degree
    public :: reparameterize_patch, reverse_patch, translate_patch
    public :: multipatch_join, multipatch_side_label, initialize_assembly, append_assembly_patch, append_assembly_join, inspect_native_assembly
    public :: native_assembly_control_map, native_assembly_control_correction, multipatch_join_directions
    public :: forcad_mesh_samples, sample_forcad_patch, sample_forcad_assembly, volume_assembly_size, check_forcad_mesh_joins, join_forcad_volume_mesh, estimate_forcad_mesh_samples
    public :: write_iges_worker_output, read_iges_worker_output, write_native_iges, write_iges_units
    public :: write_native_vtk, write_native_multipatch_vtk, write_mesh_vtk, native_multipatch_vtk_path_limit
    public :: native_patch_data, collect_native_patch_data, restore_native_patch_data
    public :: stage_box_volume_patches, stage_curved_volume_patches
    public :: curved_core_count, curved_volume_sectors
end module studio_forcad
