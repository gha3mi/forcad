module studio_single_patch
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: iso_fortran_env, only: int64
    use studio_forcad, only: rk, affine_transform_t, studio_native_patch, native_patch_kind, native_knot_axis, studio_native_problem
    use studio_units, only: unitless
    use studio_history_storage, only: history_blob
    implicit none(type, external)
    private
    public :: single_patch_info, read_single_patch
    public :: single_patch_axis, single_patch_scene, single_patch_plot
    public :: read_single_patch_scene, read_single_patch_plot, read_single_patch_slice
    public :: studio_native_problem
    public :: single_patch_controls, read_single_patch_controls, commit_single_patch_controls
    public :: commit_single_patch_knot, commit_single_patch_degree
    public :: commit_single_patch_refinement, single_patch_refinement_limits
    public :: commit_single_patch_interval, commit_single_patch_knot_values
    public :: commit_single_patch_split, commit_single_patch_extract, commit_single_patch_reverse
    public :: set_active_patch, active_patch_id
    public :: single_patch_count, single_patch_name, add_single_patch, single_patch_available
    public :: rename_single_patch, remove_single_patch, remove_single_patches
    public :: single_patch_can_undo, undo_single_patch, single_patch_undo_stamp
    public :: studio_patch_replay, redo_single_patch, release_patch_replay, clear_patch_history, get_patch_undo_targets
    public :: usable_undo_stamp, discard_undo_through
    public :: single_patch_revision, single_patch_mesh_revision, next_patch_revision, copy_single_patch_geometry
    public :: studio_native_patch, publish_native_patches, native_patch_kind, commit_patch_transform
    public :: commit_scene_transforms, scene_transform_step, commit_native_patch_edits
    public :: read_single_patch_bounds, stage_native_patch_refinement

    integer, parameter, public :: scene_rotate=1, scene_translate=2

    type :: scene_transform_step
        integer :: op=0, axis=0
        real(c_double) :: value=0._c_double
        integer :: unit=unitless
    end type scene_transform_step

    public :: studio_patch_catalog,allocate_studio_patch_catalog,stage_studio_patch,studio_patch_is_staged,move_studio_patch_catalog
    public :: initialize_studio_patch_catalog,install_studio_patch_catalog,studio_patch_resolution,patch_history_project

    type :: single_patch_info
        logical :: valid=.false.
        integer :: dimension=0, degree(3)=0, controls(3)=0
        real(c_double) :: bounds(2,3)=0.0_c_double, parameters(3)=0.0_c_double
        real(c_double) :: point(3)=0.0_c_double, jacobian(3,3)=0.0_c_double
        real(c_double) :: cell(2,3)=0.0_c_double, support(2,3)=0.0_c_double
        real(c_double), allocatable :: knots(:)
        character(:), allocatable :: message,name,diagnostic
    end type single_patch_info

    type :: single_patch_axis
        real(c_double), allocatable :: knots(:)
    end type single_patch_axis

    type :: single_patch_scene
        logical :: valid=.false.
        integer :: dimension=0, variant=0, degree(3)=0, sample_counts(3)=1, control_counts(3)=1
        real(c_double) :: bounds(2,3)=0.0_c_double
        ! Coordinate rows, with the first parameter direction varying fastest.
        real(c_double), allocatable :: points(:,:), controls(:,:)
        type(single_patch_axis) :: axes(3)
        character(:), allocatable :: diagnostic
    end type single_patch_scene

    type :: single_patch_plot
        logical :: valid=.false.
        real(c_double), allocatable :: parameters(:), values(:,:)
    end type single_patch_plot

    type :: single_patch_controls
        logical :: valid=.false.
        integer :: dimension=0, counts(3)=1
        real(c_double), allocatable :: points(:,:), weights(:)
        character(:), allocatable :: name
    end type single_patch_controls

    integer, parameter :: curve_samples=81, tensor_samples=17, plot_samples=81

    type :: patch_entry
        class(*), allocatable :: geometry
        character(:), allocatable :: name
        integer :: dimension=0, patch_id=0, degree(3)=0, controls(3)=1
        integer :: resolution(3)=0
        integer(int64) :: revision=0_int64, mesh_revision=0_int64
        logical :: removed=.false.
        real(rk) :: bounds(2,3)=0._rk
        type(native_knot_axis) :: axes(3)
    end type patch_entry
    ! ForCAD Studio takes query copies on the GUI thread so native diagnostics and
    ! sample caches stay local to each operation.
    integer, parameter :: builtin_samples=2
    type :: patch_catalog
        type(patch_entry), allocatable :: entries(:)
    end type patch_catalog
    type :: patch_edit
        integer :: kind=0, axis=0, created_patch_id=0, previous_selection=0
        integer :: multipatch_id=0
        integer(int64) :: stamp=0_int64
        character(:), allocatable :: operation
        type(patch_entry) :: previous
        type(patch_entry), allocatable :: previous_set(:)
        integer :: batch_counts(3)=-1, batch_selection(3)=0
        integer, allocatable :: targets(:,:)
        logical :: select_restored=.false.
        type(history_blob) :: before
    end type patch_edit

    type :: studio_patch_replay
        private
        type(patch_edit) :: edit
        type(history_blob) :: after
        integer, allocatable :: targets(:,:)
        integer :: counts(3)=0,selection(3)=0
    end type studio_patch_replay

    type :: studio_patch_catalog
        private
        type(patch_catalog) :: catalogs(3)
        integer :: selection(3)=builtin_samples
        type(patch_edit), allocatable :: undo_history(:)
        character(:), allocatable :: history_project
        integer :: undo_count=0
        integer(int64) :: undo_floor=0_int64, template_revision=0_int64
    contains
        procedure :: patch_count => catalog_patch_count
        procedure :: patch_available => catalog_patch_available
        procedure :: patch_revision => catalog_patch_revision
        procedure :: patch_name => catalog_patch_name
        procedure :: set_history_project
    end type studio_patch_catalog

    type(studio_patch_catalog), target :: default_catalog
    integer(int64) :: revision_clock=0_int64

    interface
        module subroutine initialize_studio_patch_catalog(catalog)
            implicit none(type, external)
            type(studio_patch_catalog), intent(inout) :: catalog
        end subroutine initialize_studio_patch_catalog

        module subroutine allocate_studio_patch_catalog(staged,counts,selection,status)
            implicit none(type, external)
            type(studio_patch_catalog), intent(out) :: staged
            integer, intent(in) :: counts(3),selection(3)
            integer, intent(out) :: status
        end subroutine allocate_studio_patch_catalog

        module subroutine move_studio_patch_catalog(source,destination)
            implicit none(type, external)
            type(studio_patch_catalog), intent(inout) :: source
            type(studio_patch_catalog), intent(out) :: destination
        end subroutine move_studio_patch_catalog

        module subroutine stage_studio_patch(staged,kind,variant,name,geometry,resolution,error,sampled)
            implicit none(type, external)
            type(studio_patch_catalog), intent(inout) :: staged
            integer, intent(in) :: kind,variant,resolution(3)
            character(*), intent(in) :: name
            class(*), allocatable, intent(inout) :: geometry
            character(:), allocatable, intent(out) :: error
            type(single_patch_scene), optional, intent(out) :: sampled
        end subroutine stage_studio_patch

        module logical function studio_patch_is_staged(staged,kind,variant) result(available)
            implicit none(type, external)
            type(studio_patch_catalog), intent(in) :: staged
            integer, intent(in) :: kind,variant
        end function studio_patch_is_staged

        module subroutine install_studio_patch_catalog(staged,catalog)
            implicit none(type, external)
            type(studio_patch_catalog), intent(inout) :: staged
            type(studio_patch_catalog), intent(inout), target, optional :: catalog
        end subroutine install_studio_patch_catalog

        module function studio_patch_resolution(kind,variant,catalog) result(resolution)
            implicit none(type, external)
            integer, intent(in) :: kind,variant
            type(studio_patch_catalog), intent(in), optional :: catalog
            integer :: resolution(3)
        end function studio_patch_resolution

        module subroutine set_active_patch(kind,patch_id,catalog)
            implicit none(type, external)
            integer, intent(in) :: kind,patch_id
            type(studio_patch_catalog), intent(inout), target, optional :: catalog
        end subroutine set_active_patch

        module integer function active_patch_id(kind,catalog) result(patch_id)
            implicit none(type, external)
            integer, intent(in) :: kind
            type(studio_patch_catalog), intent(in), optional :: catalog
        end function active_patch_id

        module pure integer function single_patch_count(kind,catalog) result(number)
            implicit none(type, external)
            integer, intent(in) :: kind
            type(studio_patch_catalog), intent(in), optional :: catalog
        end function single_patch_count

        module pure integer function catalog_patch_count(catalog,kind) result(number)
            implicit none(type, external)
            class(studio_patch_catalog), intent(in) :: catalog
            integer, intent(in) :: kind
        end function catalog_patch_count

        module pure logical function single_patch_available(kind,variant,catalog) result(available)
            implicit none(type, external)
            integer, intent(in) :: kind,variant
            type(studio_patch_catalog), intent(in), optional :: catalog
        end function single_patch_available

        module pure logical function catalog_patch_available(catalog,kind,variant) result(available)
            implicit none(type, external)
            class(studio_patch_catalog), intent(in) :: catalog
            integer, intent(in) :: kind,variant
        end function catalog_patch_available

        module integer(int64) function single_patch_revision(kind,variant,catalog) result(revision)
            implicit none(type, external)
            integer, intent(in) :: kind,variant
            type(studio_patch_catalog), intent(in), optional :: catalog
        end function single_patch_revision

        module integer(int64) function single_patch_mesh_revision(kind,variant,catalog) result(revision)
            implicit none(type, external)
            integer, intent(in) :: kind,variant
            type(studio_patch_catalog), intent(in), optional :: catalog
        end function single_patch_mesh_revision

        module integer(int64) function catalog_patch_revision(catalog,kind,variant,mesh) result(revision)
            implicit none(type, external)
            class(studio_patch_catalog), intent(in) :: catalog
            integer, intent(in) :: kind,variant
            logical, intent(in), optional :: mesh
        end function catalog_patch_revision

        module function single_patch_name(kind,variant,catalog) result(name)
            implicit none(type, external)
            integer, intent(in) :: kind,variant
            type(studio_patch_catalog), intent(in), optional :: catalog
            character(:), allocatable :: name
        end function single_patch_name

        module function catalog_patch_name(catalog,kind,variant) result(name)
            implicit none(type, external)
            class(studio_patch_catalog), intent(in) :: catalog
            integer, intent(in) :: kind,variant
            character(:), allocatable :: name
        end function catalog_patch_name

        module subroutine rename_single_patch(kind,variant,name,error,catalog)
            implicit none(type, external)
            integer, intent(in) :: kind,variant
            character(*), intent(in) :: name
            character(:), allocatable, intent(out) :: error
            type(studio_patch_catalog), intent(inout), target, optional :: catalog
        end subroutine rename_single_patch

        module subroutine remove_single_patch(kind,variant,error,catalog)
            implicit none(type, external)
            integer, intent(in) :: kind,variant
            character(:), allocatable, intent(out) :: error
            type(studio_patch_catalog), intent(inout), target, optional :: catalog
        end subroutine remove_single_patch

        module subroutine remove_single_patches(targets,error,catalog)
            implicit none(type, external)
            integer, intent(in) :: targets(:,:)
            character(:), allocatable, intent(out) :: error
            type(studio_patch_catalog), intent(inout), target :: catalog
        end subroutine remove_single_patches

        module subroutine add_single_patch(kind,source_variant,name,variant,error,offset,catalog)
            implicit none(type, external)
            integer, intent(in) :: kind,source_variant
            character(*), intent(in) :: name
            integer, intent(out) :: variant
            character(:), allocatable, intent(out) :: error
            real(c_double), intent(in), optional :: offset(3)
            type(studio_patch_catalog), intent(inout), target, optional :: catalog
        end subroutine add_single_patch

        module subroutine publish_native_patches(geometries,names,kind,variant,error,record,resolutions,operation,catalog)
            implicit none(type, external)
            type(studio_native_patch), intent(in) :: geometries(:)
            character(*), intent(in) :: names(:)
            integer, intent(out) :: kind,variant
            character(:), allocatable, intent(out) :: error
            logical, intent(in), optional :: record
            integer,intent(in),optional :: resolutions(:,:)
            character(*),intent(in),optional :: operation
            type(studio_patch_catalog), intent(inout), target, optional :: catalog
        end subroutine publish_native_patches

        module function available_patch_name(requested,batch,catalog) result(name)
            implicit none(type, external)
            character(*), intent(in) :: requested
            type(patch_entry), intent(in) :: batch(:)
            type(studio_patch_catalog), intent(in), optional :: catalog
            character(:), allocatable :: name
        end function available_patch_name

        module logical function resolve_patch_snapshot(kind,patch,patch_id,catalog) result(valid)
            implicit none(type, external)
            integer,intent(in)::kind
            type(patch_entry),intent(out)::patch
            integer,intent(in),optional::patch_id
            type(studio_patch_catalog), intent(in), optional :: catalog
        end function resolve_patch_snapshot

        module logical function copy_patch_snapshot(kind,variant,patch,catalog) result(valid)
            implicit none(type, external)
            integer, intent(in) :: kind,variant
            type(patch_entry), intent(out) :: patch
            type(studio_patch_catalog), intent(in), optional :: catalog
        end function copy_patch_snapshot

        module logical function read_catalog_patch(kind,variant,patch,catalog) result(valid)
            implicit none(type, external)
            integer, intent(in) :: kind,variant
            type(patch_entry), intent(out) :: patch
            type(studio_patch_catalog), intent(in) :: catalog
        end function read_catalog_patch

        module subroutine store_patch_template(patch,catalog)
            implicit none(type, external)
            type(patch_entry), intent(in) :: patch
            type(studio_patch_catalog), intent(inout) :: catalog
        end subroutine store_patch_template

        module function edit_catalog(catalog) result(owner)
            implicit none(type, external)
            type(studio_patch_catalog), intent(inout), target, optional :: catalog
            type(studio_patch_catalog), pointer :: owner
        end function edit_catalog

        module subroutine copy_single_patch_geometry(kind,variant,geometry,catalog)
            implicit none(type, external)
            integer, intent(in) :: kind,variant
            class(*), allocatable, intent(out) :: geometry
            type(studio_patch_catalog), intent(in), optional :: catalog
        end subroutine copy_single_patch_geometry

        module logical function refresh_patch_metadata(patch) result(valid)
            implicit none(type, external)
            type(patch_entry), intent(inout) :: patch
        end function refresh_patch_metadata

    end interface

    interface
        module subroutine read_single_patch(kind,normalized,axis,info,patch_id,catalog)
            implicit none(type, external)
            integer, intent(in) :: kind,axis
            real(c_double), intent(in) :: normalized(3)
            type(single_patch_info), intent(out) :: info
            integer, intent(in), optional :: patch_id
            type(studio_patch_catalog), intent(in), optional :: catalog
        end subroutine read_single_patch

        module subroutine read_single_patch_scene(kind,scene,variant,catalog)
            implicit none(type, external)
            integer, intent(in) :: kind
            type(single_patch_scene), intent(out) :: scene
            integer, intent(in), optional :: variant
            type(studio_patch_catalog), intent(in), optional :: catalog
        end subroutine read_single_patch_scene

        module subroutine read_single_patch_plot(kind,normalized,axis,derivative,plot,patch_id,catalog)
            implicit none(type, external)
            integer, intent(in) :: kind,axis
            real(c_double), intent(in) :: normalized(3)
            logical, intent(in) :: derivative
            type(single_patch_plot), intent(out) :: plot
            integer, intent(in), optional :: patch_id
            type(studio_patch_catalog), intent(in), optional :: catalog
        end subroutine read_single_patch_plot

        module subroutine read_single_patch_slice(normalized,scene,points,catalog)
            implicit none(type, external)
            real(c_double), intent(in) :: normalized(3)
            type(single_patch_scene), intent(in) :: scene
            real(c_double), allocatable, intent(out) :: points(:,:)
            type(studio_patch_catalog), intent(in), optional :: catalog
        end subroutine read_single_patch_slice

        module subroutine read_single_patch_controls(kind,controls,variant,catalog)
            implicit none(type, external)
            integer, intent(in) :: kind
            type(single_patch_controls), intent(out) :: controls
            integer, intent(in), optional :: variant
            type(studio_patch_catalog), intent(in), optional :: catalog
        end subroutine read_single_patch_controls

        module logical function can_sample_patch(candidate,scene) result(valid)
            implicit none(type, external)
            type(patch_entry), intent(in) :: candidate
            type(single_patch_scene), optional, intent(out) :: scene
        end function can_sample_patch

    end interface

    interface
        module subroutine commit_single_patch_controls(kind,controls,changed,error,catalog,patch_id)
            implicit none(type, external)
            integer, intent(in) :: kind
            type(single_patch_controls), intent(in) :: controls
            logical, intent(out) :: changed
            character(:), allocatable, intent(out) :: error
            type(studio_patch_catalog), intent(inout), target, optional :: catalog
            integer, intent(in), optional :: patch_id
        end subroutine commit_single_patch_controls

        module subroutine commit_single_patch_knot(kind,axis,knot,repeats,remove,accepted,error,catalog,patch_id)
            implicit none(type, external)
            integer, intent(in) :: kind,axis,repeats
            real(c_double), intent(in) :: knot
            logical, intent(in) :: remove
            integer, intent(out) :: accepted
            character(:), allocatable, intent(out) :: error
            type(studio_patch_catalog), intent(inout), target, optional :: catalog
            integer, intent(in), optional :: patch_id
        end subroutine commit_single_patch_knot

        module subroutine commit_single_patch_degree(kind,axis,steps,reduce,tolerance,accepted,error_bound,error,catalog,patch_id)
            implicit none(type, external)
            integer, intent(in) :: kind,axis,steps
            logical, intent(in) :: reduce
            real(c_double), intent(in) :: tolerance
            integer, intent(out) :: accepted
            real(c_double), intent(out) :: error_bound
            character(:), allocatable, intent(out) :: error
            type(studio_patch_catalog), intent(inout), target, optional :: catalog
            integer, intent(in), optional :: patch_id
        end subroutine commit_single_patch_degree

        module subroutine commit_single_patch_refinement(kind,axis,subdivisions,degree_steps,changed,error,catalog,patch_id)
            implicit none(type, external)
            integer, intent(in) :: kind,axis,subdivisions,degree_steps
            logical, intent(out) :: changed
            character(:), allocatable, intent(out) :: error
            type(studio_patch_catalog), intent(inout), target, optional :: catalog
            integer, intent(in), optional :: patch_id
        end subroutine commit_single_patch_refinement

        module subroutine single_patch_refinement_limits(kind,axis,max_divisions,max_steps,catalog,patch_id)
            implicit none(type, external)
            integer, intent(in) :: kind,axis
            integer, intent(out), optional :: max_divisions,max_steps
            type(studio_patch_catalog), intent(in), optional :: catalog
            integer, intent(in), optional :: patch_id
        end subroutine single_patch_refinement_limits

        module subroutine commit_single_patch_interval(kind,axis,lower,upper,changed,error,catalog,patch_id)
            implicit none(type, external)
            integer, intent(in) :: kind,axis
            real(c_double), intent(in) :: lower,upper
            logical, intent(out) :: changed
            character(:), allocatable, intent(out) :: error
            type(studio_patch_catalog), intent(inout), target, optional :: catalog
            integer, intent(in), optional :: patch_id
        end subroutine commit_single_patch_interval

        module subroutine commit_single_patch_knot_values(kind,axis,knots,changed,error,catalog,patch_id)
            implicit none(type, external)
            integer, intent(in) :: kind,axis
            real(c_double), intent(in) :: knots(:)
            logical, intent(out) :: changed
            character(:), allocatable, intent(out) :: error
            type(studio_patch_catalog), intent(inout), target, optional :: catalog
            integer, intent(in), optional :: patch_id
        end subroutine commit_single_patch_knot_values

        module subroutine commit_single_patch_split(kind,axis,value,created,error,catalog,patch_id)
            implicit none(type, external)
            integer,intent(in) :: kind,axis
            real(c_double),intent(in) :: value
            integer,intent(out) :: created
            character(:),allocatable,intent(out) :: error
            type(studio_patch_catalog), intent(inout), target, optional :: catalog
            integer, intent(in), optional :: patch_id
        end subroutine commit_single_patch_split

        module subroutine commit_single_patch_extract(kind,axis,lower,upper,created,error,catalog,patch_id)
            implicit none(type, external)
            integer,intent(in) :: kind,axis
            real(c_double),intent(in) :: lower,upper
            integer,intent(out) :: created
            character(:),allocatable,intent(out) :: error
            type(studio_patch_catalog), intent(inout), target, optional :: catalog
            integer, intent(in), optional :: patch_id
        end subroutine commit_single_patch_extract

        module subroutine commit_single_patch_reverse(kind,axis,changed,error,catalog,patch_id)
            implicit none(type, external)
            integer,intent(in) :: kind,axis
            logical,intent(out) :: changed
            character(:),allocatable,intent(out) :: error
            type(studio_patch_catalog), intent(inout), target, optional :: catalog
            integer, intent(in), optional :: patch_id
        end subroutine commit_single_patch_reverse

    end interface

    interface
        module subroutine commit_patch_transform(kind,variants,transform,changed,error,native_model,multipatch_id,catalog)
            implicit none(type, external)
            integer, intent(in) :: kind,variants(:)
            type(affine_transform_t), intent(in) :: transform
            logical, intent(out) :: changed
            character(:), allocatable, intent(out) :: error
            class(*), intent(in), optional :: native_model
            integer, intent(in), optional :: multipatch_id
            type(studio_patch_catalog), intent(inout), target, optional :: catalog
        end subroutine commit_patch_transform

        module logical function read_single_patch_bounds(kind,variant,catalog,bounds) result(valid)
            implicit none(type, external)
            integer, intent(in) :: kind,variant
            type(studio_patch_catalog), intent(in) :: catalog
            real(rk), intent(out) :: bounds(2,3)
        end function read_single_patch_bounds

        module subroutine stage_native_patch_refinement(kind,variant,directions,steps,catalog,staged,error,knot_values,mapped)
            implicit none(type, external)
            integer, intent(in) :: kind,variant,steps
            logical, intent(in) :: directions(3),mapped(3)
            type(studio_patch_catalog), intent(in) :: catalog
            type(studio_native_patch), intent(out) :: staged
            character(:), allocatable, intent(out) :: error
            real(rk), intent(inout), optional :: knot_values(3)
        end subroutine stage_native_patch_refinement

        module subroutine commit_native_patch_edits(kind,variants,patches,multipatch_id,changed,error,catalog,operation,preserve_mesh)
            implicit none(type, external)
            integer, intent(in) :: kind,variants(:),multipatch_id
            type(studio_native_patch), intent(inout) :: patches(:)
            logical, intent(out) :: changed
            character(:), allocatable, intent(out) :: error
            type(studio_patch_catalog), intent(inout), target :: catalog
            character(*), intent(in) :: operation
            logical, intent(in), optional :: preserve_mesh
        end subroutine commit_native_patch_edits

        module subroutine commit_scene_transforms(targets,steps,model_unit,changed,error,catalog)
            implicit none(type, external)
            integer, intent(in) :: targets(:,:),model_unit
            type(scene_transform_step), intent(in) :: steps(:)
            logical, intent(out) :: changed
            character(:), allocatable, intent(out) :: error
            type(studio_patch_catalog), intent(inout), target, optional :: catalog
        end subroutine commit_scene_transforms

    end interface

    interface
        module integer(int64) function next_patch_revision() result(revision)
            implicit none(type, external)
        end function next_patch_revision

        module subroutine advance_patch_revision(patch,mesh_revision)
            implicit none(type, external)
            type(patch_entry), intent(inout) :: patch
            integer(int64), intent(in), optional :: mesh_revision
        end subroutine advance_patch_revision

        module logical function single_patch_can_undo(catalog)
            implicit none(type, external)
            type(studio_patch_catalog), intent(in), optional :: catalog
        end function single_patch_can_undo

        module integer(int64) function single_patch_undo_stamp(catalog) result(stamp)
            implicit none(type, external)
            type(studio_patch_catalog), intent(in), optional :: catalog
        end function single_patch_undo_stamp

        module integer(int64) function usable_undo_stamp(stamp,catalog) result(usable)
            implicit none(type, external)
            integer(int64), intent(in) :: stamp
            type(studio_patch_catalog), intent(in), optional :: catalog
        end function usable_undo_stamp

        module subroutine discard_undo_through(stamp,catalog)
            implicit none(type, external)
            integer(int64), intent(in) :: stamp
            type(studio_patch_catalog), intent(inout), target, optional :: catalog
        end subroutine discard_undo_through

        module subroutine undo_single_patch(kind,error,operation,axis,multipatch_id,catalog,replay)
            implicit none(type, external)
            integer, intent(out) :: kind
            character(:), allocatable, intent(out) :: error
            character(:), allocatable, optional, intent(out) :: operation
            integer, optional, intent(out) :: axis,multipatch_id
            type(studio_patch_catalog), intent(inout), target, optional :: catalog
            type(studio_patch_replay), intent(out), optional :: replay
        end subroutine undo_single_patch

        module subroutine redo_single_patch(replay,error,catalog)
            implicit none(type, external)
            type(studio_patch_replay), intent(inout) :: replay
            character(:), allocatable, intent(out) :: error
            type(studio_patch_catalog), intent(inout), target, optional :: catalog
        end subroutine redo_single_patch

        module subroutine release_patch_replay(replay)
            implicit none(type, external)
            type(studio_patch_replay), intent(inout) :: replay
        end subroutine release_patch_replay

        module subroutine clear_patch_history(catalog)
            implicit none(type, external)
            type(studio_patch_catalog), intent(inout) :: catalog
        end subroutine clear_patch_history

        module subroutine get_patch_undo_targets(catalog,targets)
            implicit none(type, external)
            type(studio_patch_catalog), intent(in) :: catalog
            integer, allocatable, intent(out) :: targets(:,:)
        end subroutine get_patch_undo_targets

        module subroutine record_undo(edit,catalog,error)
            implicit none(type, external)
            type(patch_edit), intent(in) :: edit
            type(studio_patch_catalog), intent(inout), target, optional :: catalog
            character(:), allocatable, intent(out) :: error
        end subroutine record_undo

        module subroutine prepare_patch_undo(edit,catalog,stored,error)
            implicit none(type, external)
            type(patch_edit), intent(in) :: edit
            type(studio_patch_catalog), intent(in) :: catalog
            type(patch_edit), intent(out) :: stored
            character(:), allocatable, intent(out) :: error
        end subroutine prepare_patch_undo

        module subroutine prepare_patch_set_undo(edit,catalog,targets,stored,error)
            implicit none(type, external)
            type(patch_edit), intent(in) :: edit
            type(studio_patch_catalog), intent(in) :: catalog
            integer, intent(in) :: targets(:,:)
            type(patch_edit), intent(out) :: stored
            character(:), allocatable, intent(out) :: error
        end subroutine prepare_patch_set_undo

        module subroutine prepare_patch_snapshot_undo(edit,entries,targets,stored,error,catalog)
            implicit none(type, external)
            type(patch_edit), intent(in) :: edit
            type(studio_patch_catalog), intent(in), optional :: catalog
            type(patch_entry), intent(in) :: entries(:)
            integer, intent(in) :: targets(:,:)
            type(patch_edit), intent(out) :: stored
            character(:), allocatable, intent(out) :: error
        end subroutine prepare_patch_snapshot_undo

        module subroutine move_patch_entry(source,destination)
            implicit none(type, external)
            type(patch_entry), intent(inout) :: source
            type(patch_entry), intent(out) :: destination
        end subroutine move_patch_entry

        module subroutine append_patch_undo(catalog,edit)
            implicit none(type, external)
            type(studio_patch_catalog), intent(inout) :: catalog
            type(patch_edit), intent(inout) :: edit
        end subroutine append_patch_undo

        module subroutine release_patch_undo(edit)
            implicit none(type, external)
            type(patch_edit), intent(inout) :: edit
        end subroutine release_patch_undo

    end interface

contains

    subroutine set_history_project(catalog,project)
        class(studio_patch_catalog), intent(inout) :: catalog
        character(*), intent(in), optional :: project

        if(allocated(catalog%history_project))deallocate(catalog%history_project)
        if(present(project))then
            if(len_trim(project)>0)catalog%history_project=project
        end if
    end subroutine set_history_project

    function patch_history_project(catalog) result(project)
        type(studio_patch_catalog), intent(in), optional :: catalog
        character(:), allocatable :: project

        project=''
        if(.not.present(catalog))return
        if(allocated(catalog%history_project))project=catalog%history_project
    end function patch_history_project

end module studio_single_patch
