#undef BL_LANG_CC
#define BL_LANG_FORT

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
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

      REAL_T :: transition_region
      INTEGER_T height_function_flag_global
      INTEGER_T, PARAMETER ::  DEBUG_LS_MOVE_INTERFACE = 0
      INTEGER_T global_nten
      INTEGER_T physbc(SDIM,2)
      REAL_T physbc_value(SDIM,2)

      INTEGER_T normal_probe_size
      INTEGER_T ngrow_expansion
      INTEGER_T override_density(100)
      INTEGER_T freezing_model(100)
      INTEGER_T use_exact_temperature(100)
      REAL_T saturation_temp(100)
      REAL_T saturation_temp_curv(100)
      REAL_T saturation_temp_vel(100)
      REAL_T latent_heat(100)
      REAL_T reaction_rate(100)
      REAL_T visc_coef

      INTEGER_T CTML_FSI_numsolids
      INTEGER_T CTML_FSI_init
      INTEGER_T CTML_force_model

      INTEGER_T FSI_touch_flag
      INTEGER_T elements_generated
      INTEGER_T nFSI_all,nFSI_sub
      INTEGER_T ngrowFSI
      INTEGER_T global_nparts

      type MG_type
       REAL_T, pointer :: FSI_MF(D_DECL(:,:,:),:)
       REAL_T, pointer :: MASK_NBR_MF(D_DECL(:,:,:),:)
      end type MG_type

      type(MG_type), dimension(:), allocatable :: MG

      INTEGER_T, dimension(:), allocatable :: im_solid_map ! type: 0..nmat-1

      real(kind=8),dimension(:,:), allocatable :: dxlevel
      INTEGER,dimension(:,:), allocatable :: domlo_level
      INTEGER,dimension(:,:), allocatable :: domhi_level

contains

end module probmain_module

