! Native geometry knot polylines cached by source patch revision.
module studio_physical_knots
    use, intrinsic :: iso_fortran_env, only: int64
    use studio_forcad, only: physical_knot_grid, sample_physical_knots
    use studio_single_patch, only: studio_patch_catalog, copy_single_patch_geometry,single_patch_revision
    implicit none(type, external)
    private
    public :: physical_knot_grid, read_physical_knot_grid

    integer, parameter :: bucket_count=4096, slots_per_bucket=4
    integer(int64), parameter :: cache_budget=64_int64*1024*1024

    type :: grid_cache_entry
        integer :: kind=0,variant=0
        integer(int64) :: revision=-1_int64,touched=0_int64,bytes=0_int64
        type(physical_knot_grid) :: grid
    end type grid_cache_entry
    type(grid_cache_entry), allocatable, target :: cache(:,:)
    type(physical_knot_grid), target :: temporary_grid
    integer(int64) :: use_clock=0_int64, cached_bytes=0_int64
    integer :: next_eviction=0
contains
    subroutine read_physical_knot_grid(kind,variant,grid,catalog,cached_only)
        integer,intent(in) :: kind,variant
        type(physical_knot_grid),pointer,intent(out) :: grid
        type(studio_patch_catalog), intent(in), optional :: catalog
        logical, intent(in), optional :: cached_only
        type(physical_knot_grid) :: pending
        integer(int64) :: revision,bytes
        integer :: bucket,slot,i
        logical :: only_cached

        nullify(grid)
        if(allocated(temporary_grid%points))temporary_grid=physical_knot_grid()
        if(kind<1.or.kind>3.or.variant<1)return
        revision=single_patch_revision(kind,variant,catalog)
        if(revision<0_int64)return
        only_cached=.false.
        if(present(cached_only))only_cached=cached_only
        if(.not.allocated(cache))then
            if(only_cached)return
            allocate(cache(slots_per_bucket,bucket_count))
        end if
        bucket=1+int(modulo(revision+52_int64*variant+199_int64*kind,int(bucket_count,int64)))
        use_clock=use_clock+1_int64
        do i=1,slots_per_bucket
            associate(entry=>cache(i,bucket))
                if(entry%kind/=kind.or.entry%variant/=variant.or.entry%revision/=revision)cycle
                entry%touched=use_clock
                if(allocated(entry%grid%points))grid=>entry%grid
                return
            end associate
        end do
        if(only_cached)return

        call build_grid(kind,variant,pending,catalog)
        bytes=grid_bytes(pending)
        if(bytes>cache_budget)then
            call move_grid(pending,temporary_grid)
            if(allocated(temporary_grid%points))grid=>temporary_grid
            return
        end if
        slot=minloc(cache(:,bucket)%touched,dim=1)
        call evict_entry(slot,bucket)
        call make_room(bytes)
        associate(entry=>cache(slot,bucket))
            call move_grid(pending,entry%grid)
            entry%kind=kind
            entry%variant=variant
            entry%revision=revision
            entry%touched=use_clock
            entry%bytes=bytes
            cached_bytes=cached_bytes+bytes
            if(allocated(entry%grid%points))grid=>entry%grid
        end associate
    end subroutine read_physical_knot_grid

    subroutine make_room(bytes)
        integer(int64), intent(in) :: bytes
        integer :: slot,bucket

        do while(cached_bytes+bytes>cache_budget)
            slot=1+mod(next_eviction,slots_per_bucket)
            bucket=1+next_eviction/slots_per_bucket
            next_eviction=mod(next_eviction+1,size(cache))
            call evict_entry(slot,bucket)
        end do
    end subroutine make_room

    subroutine evict_entry(slot,bucket)
        integer, intent(in) :: slot,bucket

        associate(entry=>cache(slot,bucket))
            if(entry%kind==0)return
            cached_bytes=cached_bytes-entry%bytes
            entry=grid_cache_entry()
        end associate
    end subroutine evict_entry

    pure integer(int64) function grid_bytes(grid) result(bytes)
        type(physical_knot_grid), intent(in) :: grid

        bytes=0_int64
        if(allocated(grid%points))bytes=bytes+size(grid%points,kind=int64)*int(storage_size(grid%points)/8,int64)
        if(allocated(grid%lines))bytes=bytes+size(grid%lines,kind=int64)*int(storage_size(grid%lines)/8,int64)
        if(allocated(grid%faces))bytes=bytes+size(grid%faces,kind=int64)*int(storage_size(grid%faces)/8,int64)
    end function grid_bytes

    subroutine move_grid(source,destination)
        type(physical_knot_grid), intent(inout) :: source
        type(physical_knot_grid), intent(out) :: destination

        call move_alloc(source%points,destination%points)
        call move_alloc(source%lines,destination%lines)
        call move_alloc(source%faces,destination%faces)
    end subroutine move_grid

    subroutine build_grid(kind,variant,grid,catalog)
        integer,intent(in) :: kind,variant
        type(physical_knot_grid),intent(out) :: grid
        type(studio_patch_catalog), intent(in), optional :: catalog
        class(*),allocatable :: native

        call copy_single_patch_geometry(kind,variant,native,catalog)
        if(allocated(native))call sample_physical_knots(native,grid)
    end subroutine build_grid
end module studio_physical_knots
