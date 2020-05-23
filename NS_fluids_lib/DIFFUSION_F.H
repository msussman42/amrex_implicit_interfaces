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

#include <AMReX_ArrayLim.H>

namespace amrex{

extern "C"
{
    void FORT_HOOPIMPLICIT (
      int* override_density,
      Real* gravity_normalized,
      int* gravity_dir,
      Real *force,AMREX_ARLIM_P(forcelo),AMREX_ARLIM_P(forcehi),
      Real *tensor,AMREX_ARLIM_P(tensorlo),AMREX_ARLIM_P(tensorhi),
      Real *thermal,AMREX_ARLIM_P(thermallo),AMREX_ARLIM_P(thermalhi),
      Real *recon,AMREX_ARLIM_P(reconlo),AMREX_ARLIM_P(reconhi),
      Real *solxfab,
      AMREX_ARLIM_P(solxfablo),AMREX_ARLIM_P(solxfabhi),
      Real *solyfab,
      AMREX_ARLIM_P(solyfablo),AMREX_ARLIM_P(solyfabhi),
      Real *solzfab,
      AMREX_ARLIM_P(solzfablo),AMREX_ARLIM_P(solzfabhi),
      const Real *xlo,const Real* dx,
      Real *uold,AMREX_ARLIM_P(uoldlo),AMREX_ARLIM_P(uoldhi),
      Real *unew,AMREX_ARLIM_P(unewlo),AMREX_ARLIM_P(unewhi),
      Real *lsnew,AMREX_ARLIM_P(lsnewlo),AMREX_ARLIM_P(lsnewhi),
      Real *den,AMREX_ARLIM_P(denlo),AMREX_ARLIM_P(denhi),  
      Real *mu,AMREX_ARLIM_P(mulo),AMREX_ARLIM_P(muhi),
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
      Real *force,AMREX_ARLIM_P(forcelo),AMREX_ARLIM_P(forcehi),
      const Real *xlo,const Real* dx,
      Real *unew,AMREX_ARLIM_P(unewlo),AMREX_ARLIM_P(unewhi),
      Real *den,AMREX_ARLIM_P(denlo),AMREX_ARLIM_P(denhi),  
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
      Real *force,AMREX_ARLIM_P(forcelo),AMREX_ARLIM_P(forcehi),
      Real *thermal,AMREX_ARLIM_P(thermallo),AMREX_ARLIM_P(thermalhi),
      Real *recon,AMREX_ARLIM_P(reconlo),AMREX_ARLIM_P(reconhi),
      const Real *xlo,const Real* dx,
      Real *uold,AMREX_ARLIM_P(uoldlo),AMREX_ARLIM_P(uoldhi),
      Real *snew,AMREX_ARLIM_P(snewlo),AMREX_ARLIM_P(snewhi),
      Real *den,AMREX_ARLIM_P(denlo),AMREX_ARLIM_P(denhi),   // 1/den
      Real *DEDT,AMREX_ARLIM_P(DEDTlo),AMREX_ARLIM_P(DEDThi), // 1/(rho cv)
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

}/* namespace amrex */

#endif
#endif /*_DIFFUSION_F_H*/
