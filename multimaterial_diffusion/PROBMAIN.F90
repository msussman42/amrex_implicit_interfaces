#undef BL_LANG_CC
#define BL_LANG_FORT

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_BC_TYPES.H"

#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else  
print *,"dimension bust"
stop
#endif

module probmain_module
implicit none

      INTEGER_T global_nten
      INTEGER_T physbc(SDIM,2)
      REAL_T physbc_value(SDIM,2)

      INTEGER_T normal_probe_size
      INTEGER_T ngrow_expansion
      INTEGER_T override_density(100)
      INTEGER_T freezing_model(100)
      INTEGER_T use_exact_temperature(100)
      REAL_T saturation_temp(100)
      REAL_T latent_heat(100)
      REAL_T reaction_rate(100)
      REAL_T visc_coef
contains

end module probmain_module

