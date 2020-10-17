#ifndef _PROB_CPP_PARMS_F_H_
#define _PROB_CPP_PARMS_F_H_

#ifdef BL_LANG_FORT
#    define FORT_OVERRIDE    override
#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_OVERRIDE    OVERRIDE
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_OVERRIDE    override
#  else
#    define FORT_OVERRIDE    override_
#  endif
#  endif

#include "AMReX_ArrayLim.H"

namespace amrex{

  extern "C" {

     void FORT_OVERRIDE(
       const int* max_level,
       const int* bfact_space_order,
       const int* bfact_time_order,
       const int* prescribe_temperature_outflow,
       const int* solidheat_flag,
       const int* rz_flag,
       const int* FSI_flag,
       const int* ZEYU_DCA_SELECT,
       const int* invert_solid_levelset,
       const Real* denfact,
       const Real* velfact,
       int* n_sites,
       const Real* nucleation_init_time,
       Real* pos_sites,
       const Real* xblob,const Real* yblob,
       const Real* zblob,const Real* radblob,
       const Real* xblob2,const Real* yblob2,
       const Real* zblob2,const Real* radblob2,
       const Real* xblob3,const Real* yblob3,
       const Real* zblob3,const Real* radblob3,
       const Real* xblob4,const Real* yblob4,
       const Real* zblob4,const Real* radblob4,
       const Real* xblob5,const Real* yblob5,
       const Real* zblob5,const Real* radblob5,
       const Real* xblob6,const Real* yblob6,
       const Real* zblob6,const Real* radblob6,
       const Real* xblob7,const Real* yblob7,
       const Real* zblob7,const Real* radblob7,
       const Real* xblob8,const Real* yblob8,
       const Real* zblob8,const Real* radblob8,
       const Real* xblob9,const Real* yblob9,
       const Real* zblob9,const Real* radblob9,
       const Real* xblob10,const Real* yblob10,
       const Real* zblob10,const Real* radblob10,
       const Real* xactive,
       const Real* yactive,
       const Real* zactive,
       const Real* ractivex,
       const Real* ractivey,
       const Real* ractivez,
       const int* probtype,
       const int* adv_dir,
       const Real* adv_vel,
       const int* axis_dir,
       const Real* rgasinlet,
       const Real* vinletgas,
       const Real* twall,
       const Real* advbot,
       const Real* inflow_pressure,
       const Real* outflow_pressure,
       const Real* period_time,
       Real* problox,Real* probloy,Real* probloz,
       Real* probhix,Real* probhiy,Real* probhiz,
       int* num_species_var,
       int* num_materials_viscoelastic,
       int* num_state_material, 
       int* num_state_base, 
       int* ngeom_raw,
       int* ngeom_recon,
       int* num_materials,
       int* num_materials_vel,
       int* num_materials_scalar_solve,
       int* material_type,
       int* nten,
       Real* drhodt,
       Real* drhodz,
       Real* tempconst,
       Real* initial_temperature,
       Real* tempcutoff,
       Real* tempcutoffmax,
       Real* stiffPINF,
       Real* stiffCP,
       Real* stiffCV,
       Real* stiffGAMMA,
       Real* denconst,
       Real* den_floor,
       Real* den_ceiling,
       Real* cavdenconst,
       Real* viscconst,
       Real* viscconst_eddy,
       int*  viscosity_state_model,
       Real* elastic_viscosity,
       int* store_elastic_data,
       Real* heatviscconst,
       Real* prerecalesce_heatviscconst,
       Real* prerecalesce_viscconst,
       Real* prerecalesce_stiffCP,
       Real* prerecalesce_stiffCV,
       Real* speciesconst,
       Real* speciesviscconst,
       Real* latent_heat,
       Real* saturation_temp,
       Real* molar_mass,
       Real* species_molar_mass,
       Real* tension,
       Real* tension_slope,
       Real* tension_T0,
       Real* tension_min,
       Real* prefreeze_tension,
       Real* FORT_MUSHY_THICK,  // default==2.0
       Real* gravity,
       int* gravity_dir,
       int* invert_gravity,
       Real* fort_stop_time,
       int* ioproc);

  }

}/* namespace amrex */

#endif
#endif 
