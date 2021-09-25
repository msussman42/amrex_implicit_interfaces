#ifndef _DIFFUSION_F_H
#define _DIFFUSION_F_H

/*
** $Id: DIFFUSION_F.H,v 1.12 2001/08/01 21:50:53 lijewski Exp $
*/

#ifdef BL_LANG_FORT
#else

#include <AMReX_ArrayLim.H>

namespace amrex{

extern "C"
{
    void fort_hoopimplicit(
      int* override_density,
      int* constant_density_all_time,
      int* nstate,
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
      int* uncoupled_viscosity,
      int* update_state,
      Real* dt,
      int* rzflag,
      int* nmat,
      int* nparts,
      int* nparts_def,
      int* im_solid_map,
      int* ntensor,
      int* nsolve);

}

}/* namespace amrex */

#endif
#endif /*_DIFFUSION_F_H*/
