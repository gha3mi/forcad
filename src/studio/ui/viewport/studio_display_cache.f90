module studio_display_cache
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: iso_fortran_env, only: int64
    use qt, only: QThread,QTimer
    use studio_selection, only: selection_set_t
    implicit none(type,external)
    private

    public :: display_grid,display_cache_entry,display_sample_job,display_view,viewport_display_cache

    type :: display_grid
        real(c_double),allocatable :: points(:,:)
    end type display_grid

    type :: display_cache_entry
        integer(int64) :: revision=-1_int64
        integer :: base(3)=1,counts(3)=1,attempted(3)=0
        real(c_double) :: error(3)=0._c_double
        logical :: needs_pilot=.false.
        type(display_grid) :: grids(6)
    end type display_cache_entry

    type :: display_cache_rank
        type(display_cache_entry),allocatable :: entries(:)
    end type display_cache_rank

    type :: display_sample_job
        integer :: kind=0,variant=0,counts(3)=1,base(3)=1
        integer(int64) :: revision=-1_int64,document=-1_int64
        logical :: valid=.false.,pilot=.false.
        real(c_double) :: bounds(2,3)=0._c_double,error(3)=0._c_double
        class(*),allocatable :: geometry
        type(display_grid) :: grids(6)
    end type display_sample_job

    type :: display_view
        integer(int64) :: revisions(3)=-1_int64
        integer :: mesh_options(3)=0,visible_count=-1
        real(c_double) :: projection(10)=0._c_double, clip(4)=0._c_double
        type(selection_set_t) :: isolation
    end type display_view

    type :: viewport_display_cache
        integer(int64) :: document=-1_int64
        integer :: next_rank=1,next_variant=1
        logical :: dirty=.false.
        type(display_view) :: view
        type(display_cache_rank) :: ranks(3)
        type(QTimer) :: timer
        type(QThread) :: worker
        type(display_sample_job),pointer :: job=>null()
    end type viewport_display_cache
end module studio_display_cache
