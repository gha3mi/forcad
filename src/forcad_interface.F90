!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Provides a unified interfaces
!> @note
!> note: Due to a known issue with preprocessing in fpm
!> (see: https://github.com/fortran-lang/fpm/issues/773),
!> conditional compilation is currently unreliable in this context.
!> As a workaround, the implementation must be selected manually.
!> @endnote
module forcad_interface

! #ifdef USE_STDLIB_SOLVE
!    use stdlib_linalg, only: solve
! #else
!    use forcad_utils, only: solve
! #endif
   use forcad_utils, only: solve
   ! use stdlib_linalg, only: solve

   implicit none

   private
   public solve
end module
