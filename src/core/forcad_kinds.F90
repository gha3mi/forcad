!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
!> Floating-point kind policy for the complete ForCAD public API.
!!
!! `rk` is selected at preprocessing time. Define at most one of `REAL32`,
!! `REAL64`, `REALXDP`, or `REAL128`; without a selection, ForCAD requests
!! approximately 15 decimal digits (`REAL64`). If several macros are defined,
!! source-order `#elif` precedence selects the first one, so such a build is
!! unsupported. The actual kind is compiler dependent because
!! `selected_real_kind` expresses a precision requirement, not a storage-size
!! requirement. If the requested precision is unavailable, it returns a
!! negative kind value and the compiler cannot build that precision mode.
!!
!! If \(p\) is the requested decimal precision, the selection is equivalent to
!! `selected_real_kind(p)`. Algorithms derive scale-sensitive tolerances from
!! `epsilon(1.0_rk)` and representability limits from `tiny(1.0_rk)` or
!! `huge(1.0_rk)`; no API assumes an IEEE storage width from the macro name.
!!
!! Client code should declare real data and literals with this public kind:
!!
!! ```fortran
!! real(rk) :: tolerance
!! tolerance = 100.0_rk*epsilon(1.0_rk)
!! ```
module forcad_kinds
   implicit none
   private
   public rk
#ifdef REAL32
   integer, parameter :: rk = selected_real_kind(6) !! Working real kind; approximately 6 decimal digits.
#elif REAL64
   integer, parameter :: rk = selected_real_kind(15) !! Working real kind; approximately 15 decimal digits.
#elif REALXDP
   integer, parameter :: rk = selected_real_kind(18) !! Working real kind; approximately 18 decimal digits.
#elif REAL128
   integer, parameter :: rk = selected_real_kind(33) !! Working real kind; approximately 33 decimal digits.
#else
   integer, parameter :: rk = selected_real_kind(15) !! Default working real kind; approximately 15 decimal digits.
#endif
end module
