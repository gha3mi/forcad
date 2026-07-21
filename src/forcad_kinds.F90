!> author: Seyed Ali Ghasemi
!> license: BSD 3-Clause
module forcad_kinds
   implicit none
   private
   public rk
#ifdef REAL32
   integer, parameter :: rk = selected_real_kind(6)
#elif REAL64
   integer, parameter :: rk = selected_real_kind(15)
#elif REALXDP
   integer, parameter :: rk = selected_real_kind(18)
#elif REAL128
   integer, parameter :: rk = selected_real_kind(33)
#else
   integer, parameter :: rk = selected_real_kind(15)
#endif
end module
