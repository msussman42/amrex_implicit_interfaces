#ifndef _DIFFUSION_F_H
#define _DIFFUSION_F_H

/*
** $Id: DIFFUSION_F.H,v 1.12 2001/08/01 21:50:53 lijewski Exp $
*/

#ifdef BL_LANG_FORT
#    define FORT_COMPUTE_NEG_MOM_FORCE compute_neg_mom_force
#    define FORT_HOOPIMPLICIT hoopimplicit
#    define FORT_THERMAL_OFFSET_FORCE thermal_offset_force
#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_COMPUTE_NEG_MOM_FORCE COMPUTE_NEG_MOM_FORCE
#    define FORT_HOOPIMPLICIT HOOPIMPLICIT
#    define FORT_THERMAL_OFFSET_FORCE THERMAL_OFFSET_FORCE
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_COMPUTE_NEG_MOM_FORCE compute_neg_mom_force
#    define FORT_HOOPIMPLICIT hoopimplicit
#    define FORT_THERMAL_OFFSET_FORCE thermal_offset_force
#  else
#    define FORT_COMPUTE_NEG_MOM_FORCE compute_neg_mom_force_
#    define FORT_HOOPIMPLICIT hoopimplicit_
#    define FORT_THERMAL_OFFSET_FORCE thermal_offset_force_
#  endif
#  endif

#include <ArrayLim.H>

extern "C"
{
    void FORT_HOOPIMPLICIT (
      int* override_density,
      Real* gravity_normalized,
      int* gravity_dir,
      Real *force, ARLIM_P(forcelo), ARLIM_P(forcehi),
      Real *tensor, ARLIM_P(tensorlo), ARLIM_P(tensorhi),
      Real *thermal, ARLIM_P(thermallo), ARLIM_P(thermalhi),
      Real *recon, ARLIM_P(reconlo), ARLIM_P(reconhi),
      Real *sol, ARLIM_P(sollo), ARLIM_P(solhi),
      const Real *xlo,const Real* dx,
      Real *uold, ARLIM_P(uoldlo), ARLIM_P(uoldhi),
      Real *unew, ARLIM_P(unewlo), ARLIM_P(unewhi),
      Real *lsnew, ARLIM_P(lsnewlo), ARLIM_P(lsnewhi),
      Real *den, ARLIM_P(denlo), ARLIM_P(denhi),  
      Real *mu, ARLIM_P(mulo), ARLIM_P(muhi),
      const int* tilelo, const int* tilehi,
      const int* fablo, const int* fabhi,
      int* bfact,
      int* level,
      int* finest_level,
      Real* visc_coef,
      Real* angular_velocity,
      int* constant_viscosity,
      int* update_state,
      Real* dt,
      int* rzflag,
      int* nmat,
      int* nparts,
      int* nparts_def,
      int* im_solid_map,
      int* ntensorMM,
      int* nsolveMM);


    void FORT_COMPUTE_NEG_MOM_FORCE (
      Real *force, ARLIM_P(forcelo), ARLIM_P(forcehi),
      const Real *xlo,const Real* dx,
      Real *unew, ARLIM_P(unewlo), ARLIM_P(unewhi),
      Real *den, ARLIM_P(denlo), ARLIM_P(denhi),  
      const int* tilelo, const int* tilehi,
      const int* fablo, const int* fabhi,
      int* bfact,
      int* level,
      int* finest_level,
      int* update_state,
      Real* dt,
      Real* prev_time,
      Real* cur_time,
      int* nmat,
      int* nsolveMM);




    void FORT_THERMAL_OFFSET_FORCE (
      int* override_density,
      Real *force, ARLIM_P(forcelo), ARLIM_P(forcehi),
      Real *thermal, ARLIM_P(thermallo), ARLIM_P(thermalhi),
      Real *recon, ARLIM_P(reconlo), ARLIM_P(reconhi),
      const Real *xlo,const Real* dx,
      Real *uold, ARLIM_P(uoldlo), ARLIM_P(uoldhi),
      Real *snew, ARLIM_P(snewlo), ARLIM_P(snewhi),
      Real *den, ARLIM_P(denlo), ARLIM_P(denhi),   // 1/den
      Real *DEDT, ARLIM_P(DEDTlo), ARLIM_P(DEDThi), // 1/(rho cv)
      const int* tilelo, const int* tilehi,
      const int* fablo, const int* fabhi,
      int* bfact,
      int* level,
      int* finest_level,
      int* update_state,
      Real* dt,
      int* rzflag,
      int* nmat,
      int* nstate,
      int* nsolveMM);

}

#endif
#endif /*_DIFFUSION_F_H*/
