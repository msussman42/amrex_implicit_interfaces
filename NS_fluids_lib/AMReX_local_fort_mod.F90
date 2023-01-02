#include <AMReX_Config.H>

module local_amrex_fort_module

  use iso_c_binding, only : c_char, c_short, c_int, c_long, c_long_long, c_float, c_double, c_size_t, c_ptr

  implicit none

  integer, parameter ::    bl_spacedim = AMREX_SPACEDIM
  integer, parameter :: amrex_spacedim = AMREX_SPACEDIM

#ifdef AMREX_USE_FLOAT
  integer, parameter :: amrex_real = c_float
  ! We could/should use Fortran 2008 c_sizeof here.
  integer (kind=c_size_t), parameter :: amrex_real_size = 4_c_size_t
#else
  integer, parameter :: amrex_real = c_double
  ! We could/should use Fortran 2008 c_sizeof here.
  integer (kind=c_size_t), parameter :: amrex_real_size = 8_c_size_t
#endif

#ifdef AMREX_SINGLE_PRECISION_PARTICLES
  integer, parameter :: amrex_particle_real = c_float
#else
  integer, parameter :: amrex_particle_real = c_double
#endif

#ifdef _WIN32
  integer, parameter :: amrex_long = c_long_long
#else
  integer, parameter :: amrex_long = c_long
#endif

contains

end module local_amrex_fort_module
