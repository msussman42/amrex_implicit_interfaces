
#ifndef _DERIVE_F_H
#define _DERIVE_F_H

/*
** $Id: DERIVE_F.H,v 1.13 1999/08/25 18:04:03 propp Exp $
*/

#ifdef BL_LANG_FORT
#    define FORT_MAXPRESVEL maxpresvel 
#    define FORT_INTEGRATE_RECALESCE integrate_recalesce
#    define FORT_RESET_TEMPERATURE reset_temperature
#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_MAXPRESVEL MAXPRESVEL
#    define FORT_INTEGRATE_RECALESCE INTEGRATE_RECALESCE
#    define FORT_RESET_TEMPERATURE RESET_TEMPERATURE
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_MAXPRESVEL maxpresvel 
#    define FORT_INTEGRATE_RECALESCE integrate_recalesce
#    define FORT_RESET_TEMPERATURE reset_temperature
#  else
#    define FORT_MAXPRESVEL maxpresvel_ 
#    define FORT_INTEGRATE_RECALESCE integrate_recalesce_
#    define FORT_RESET_TEMPERATURE reset_temperature_
#  endif
#  endif

#include <AMReX_ArrayLim.H>

namespace amrex{

extern "C" {

     void fort_getshear(
       int* ntensor,
       Real* cellten,AMREX_ARLIM_P(celltenlo),AMREX_ARLIM_P(celltenhi),
       Real* vel,AMREX_ARLIM_P(vello),AMREX_ARLIM_P(velhi),
       const Real* delta, const Real* xlo,
       Real* tensordata,AMREX_ARLIM_P(tendtlo),AMREX_ARLIM_P(tendthi),
       int* iproject,int* onlyscalar,Real* time,
       const int* tilelo, const int* tilehi,
       const int* fablo, const int* fabhi,
       int* bfact,
       int* level,
       int* bcrec,int* ngrow,
       int* nmat);

     void fort_dermagtrace(
       int* level,
       int* finest_level,
       int* im,  // im=0..nmat-1
       int* ntensor,
       Real* cellten,AMREX_ARLIM_P(celltenlo),AMREX_ARLIM_P(celltenhi),
       Real* dest,AMREX_ARLIM_P(destlo),AMREX_ARLIM_P(desthi),
       Real* den,AMREX_ARLIM_P(denlo),AMREX_ARLIM_P(denhi),
       Real* tensor,AMREX_ARLIM_P(tensorlo),AMREX_ARLIM_P(tensorhi),
       Real* vel,AMREX_ARLIM_P(vello),AMREX_ARLIM_P(velhi),
       Real* visc,AMREX_ARLIM_P(visclo),AMREX_ARLIM_P(vischi),
       const int* tilelo, const int* tilehi,
       const int* fablo, const int* fabhi,
       int* bfact,
       int* ngrow,
       const Real* dx, const Real* xlo,Real* time,
       int* bcrec,
       int* ncomp_den,
       int* ncomp_tensor,
       int* ncomp_visc,
       int* n_trace,
       int* nmat,
       Real* polymer_factor,
       Real* etaS,
       Real* etaP,
       Real* Carreau_beta,
       Real* elastic_time,
       int* viscoelastic_model,
       Real* elastic_viscosity);

     void fort_getdrag(
      int* isweep,
      Real* globalsum,Real* localsum,
      Real* gravity_normalized,
      int* gravdir,
      int* ntenvisco,
      Real* tdata,AMREX_ARLIM_P(tdatalo),AMREX_ARLIM_P(tdatahi),
      Real* viscoten,AMREX_ARLIM_P(viscotenlo),AMREX_ARLIM_P(viscotenhi),
      Real* den,AMREX_ARLIM_P(denlo),AMREX_ARLIM_P(denhi),
      Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
      Real* vof,AMREX_ARLIM_P(voflo),AMREX_ARLIM_P(vofhi),
      Real* levelpc,AMREX_ARLIM_P(levelpclo),AMREX_ARLIM_P(levelpchi),
      Real* vol,AMREX_ARLIM_P(vollo),AMREX_ARLIM_P(volhi),
      Real* areax,AMREX_ARLIM_P(areaxlo),AMREX_ARLIM_P(areaxhi),
      Real* areay,AMREX_ARLIM_P(areaylo),AMREX_ARLIM_P(areayhi),
      Real* areaz,AMREX_ARLIM_P(areazlo),AMREX_ARLIM_P(areazhi),
      Real* xface,AMREX_ARLIM_P(xfacelo),AMREX_ARLIM_P(xfacehi),
      Real* yface,AMREX_ARLIM_P(yfacelo),AMREX_ARLIM_P(yfacehi),
      Real* zface,AMREX_ARLIM_P(zfacelo),AMREX_ARLIM_P(zfacehi),
      Real* cvisc,AMREX_ARLIM_P(cvisclo),AMREX_ARLIM_P(cvischi),
      Real* c_mat_visc,
      AMREX_ARLIM_P(c_mat_visclo),AMREX_ARLIM_P(c_mat_vischi),
      const Real* xlo,const Real* dx,
      Real* solxfab,AMREX_ARLIM_P(solxfablo),AMREX_ARLIM_P(solxfabhi),
      Real* solyfab,AMREX_ARLIM_P(solyfablo),AMREX_ARLIM_P(solyfabhi),
      Real* solzfab,AMREX_ARLIM_P(solzfablo),AMREX_ARLIM_P(solzfabhi),
      Real* pres,AMREX_ARLIM_P(preslo),AMREX_ARLIM_P(preshi),
      Real* vel,AMREX_ARLIM_P(vello),AMREX_ARLIM_P(velhi),
      Real* drag,AMREX_ARLIM_P(draglo),AMREX_ARLIM_P(draghi),
      const int* tilelo, const int* tilehi, 
      const int* fablo, const int* fabhi,
      int* bfact, 
      int* rzflag,int* bcrec,
      Real* time,
      Real* visc_coef,
      int* ntensor,
      int* nmat,
      int* nparts,
      int* nparts_def,
      int* im_solid_map);


     void FORT_INTEGRATE_RECALESCE(
      int* isweep,
      Real* globalsum,Real* localsum,
      int* recalesce_material,
      Real* snew,AMREX_ARLIM_P(snewlo),AMREX_ARLIM_P(snewhi),
      Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
      Real* vol,AMREX_ARLIM_P(vollo),AMREX_ARLIM_P(volhi),
      const Real* xlo,const Real* dx,
      const int* tilelo, const int* tilehi, 
      const int* fablo, const int* fabhi,
      int* bfact, 
      int* velbc,
      Real* time,
      int* num_integrate,int* nmat,int* ncomp_state,
      int* level,int* finest_level);


     void FORT_RESET_TEMPERATURE(
      int* im_source,
      Real* TSAT,
      Real* snew,AMREX_ARLIM_P(snewlo),AMREX_ARLIM_P(snewhi),
      const Real* xlo,const Real* dx,
      const int* tilelo, const int* tilehi, 
      const int* fablo, const int* fabhi,
      int* bfact, 
      int* velbc,
      Real* time,
      int* nmat,int* ncomp_state,
      int* level,int* finest_level);


     void FORT_MAXPRESVEL(
       Real* minpres,
       Real* maxpres,
       Real* maxvel,
       Real* maxvel_collide,
       const Real* xlo,const Real* dx,
       Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
       Real* vel,AMREX_ARLIM_P(vello),AMREX_ARLIM_P(velhi),
       Real* velx,AMREX_ARLIM_P(velxlo),AMREX_ARLIM_P(velxhi),
       Real* vely,AMREX_ARLIM_P(velylo),AMREX_ARLIM_P(velyhi),
       Real* velz,AMREX_ARLIM_P(velzlo),AMREX_ARLIM_P(velzhi),
       const int* tilelo, const int* tilehi,
       const int* fablo, const int* fabhi,int* bfact);

     void fort_derviscosity(
      int* level,
      int* finest_level,
      Real* visc_coef,
      int* fortran_im,
      int* nmat,
      Real* dt,
      Real* viscosity_coefficient,
      int* shear_thinning_flag,
      Real* Carreau_alpha,
      Real* Carreau_beta,
      Real* Carreau_n,
      Real* Carreau_mu_inf,
      Real* concentration,
      Real* elastic_time,
      int* viscosity_state_model,
      int* viscoelastic_model,
      Real* elastic_viscosity,
      Real* elastic_regularization,
      Real* etaL,Real* etaP,Real* etaS,
      Real* polymer_factor,
      Real* dat,AMREX_ARLIM_P(datlo),AMREX_ARLIM_P(dathi),
      Real* vel,AMREX_ARLIM_P(vello),AMREX_ARLIM_P(velhi),
      Real* eosdata,AMREX_ARLIM_P(eosdatalo),AMREX_ARLIM_P(eosdatahi),
      Real* tensor,AMREX_ARLIM_P(tensorlo),AMREX_ARLIM_P(tensorhi),
      Real* gammadot,AMREX_ARLIM_P(gammadotlo),AMREX_ARLIM_P(gammadothi),
      const int* tilelo, const int* tilehi,
      const int* fablo, const int* fabhi,
      int* bfact,
      Real* time,
      const Real* dx, const Real* xlo,
      const int* bcrec,int* ngrow,
      int* ncompvisc);

     void fort_derconductivity(
      int* NS_sumdata_size,
      Real* NS_sumdata,
      int* ncomp_sum_int_user1,
      int* ncomp_sum_int_user2,
      int* level,
      int* finest_level,
      int* fortran_im,
      int* nmat,
      Real* dt,
      Real* dat,AMREX_ARLIM_P(datlo),AMREX_ARLIM_P(dathi),
      Real* eosdata,AMREX_ARLIM_P(eosdatalo),AMREX_ARLIM_P(eosdatahi),
      Real* vof,AMREX_ARLIM_P(voflo),AMREX_ARLIM_P(vofhi),
      const int* tilelo, const int* tilehi,
      const int* fablo, const int* fabhi,
      int* bfact,
      Real* time,
      const Real* dx, const Real* xlo,
      int* ngrow);

     void fort_derturbvisc(
      int* les_model,
      int* level,
      int* fortran_im,
      int* nmat,
      Real* dt,
      int* ntensor,
      Real* denstate,AMREX_ARLIM_P(denstatelo),AMREX_ARLIM_P(denstatehi),
      Real* vof,AMREX_ARLIM_P(voflo),AMREX_ARLIM_P(vofhi),
      Real* vel,AMREX_ARLIM_P(vello),AMREX_ARLIM_P(velhi),
      Real* visc,AMREX_ARLIM_P(visclo),AMREX_ARLIM_P(vischi),
      Real* cellten,AMREX_ARLIM_P(celltenlo),AMREX_ARLIM_P(celltenhi),
      const int* tilelo, const int* tilehi,
      const int* fablo, const int* fabhi,
      int* bfact,
      Real* cur_time,
      const Real* dx, 
      const Real* xlo,
      const Real* dx_coarsest, 
      int* ngrow,
      int* ncompvisc);

}
}/* namespace amrex */

#endif
#endif
