module studio_system_memory
    use, intrinsic :: iso_c_binding, only: c_int64_t
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none(type, external)
    private
    public :: mesh_memory_budget

    interface
        integer(c_int64_t) function available_memory() bind(C,name='studio_available_memory')
            import c_int64_t
            implicit none(type,external)
        end function available_memory
    end interface

contains

    integer(int64) function mesh_memory_budget() result(budget)
        integer(int64), parameter :: fallback_bytes=256_int64*1024*1024
        integer(int64) :: available

        available=available_memory()
        budget=fallback_bytes
        if(available>=0)budget=available/4
    end function mesh_memory_budget
end module studio_system_memory
