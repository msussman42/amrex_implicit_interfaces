
#ifndef _DERIVE_F_H
#define _DERIVE_F_H

/*
** $Id: DERIVE_F.H,v 1.13 1999/08/25 18:04:03 propp Exp $
*/

#ifdef BL_LANG_FORT
#    define FORT_GETSHEAR getshear
#    define FORT_MAXPRESVEL maxpresvel 
#    define FORT_GETDRAG  getdrag
#    define FORT_DERERROR    dererror
#    define FORT_DERVISCOSITY derviscosity
#    define FORT_DERTURBVISC derturbvisc
#    define FORT_DERMAGTRACE dermagtrace
#    define FORT_INTEGRATE_RECALESCE integrate_recalesce
#    define FORT_RESET_TEMPERATURE reset_temperature
#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_GETSHEAR GETSHEAR
#    define FORT_MAXPRESVEL MAXPRESVEL
#    define FORT_GETDRAG  GETDRAG
#    define FORT_DERERROR    DERERROR
#    define FORT_DERVISCOSITY DERVISCOSITY
#    define FORT_DERTURBVISC DERTURBVISC
#    define FORT_DERMAGTRACE DERMAGTRACE
#    define FORT_INTEGRATE_RECALESCE INTEGRATE_RECALESCE
#    define FORT_RESET_TEMPERATURE RESET_TEMPERATURE
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_GETSHEAR getshear
#    define FORT_MAXPRESVEL maxpresvel 
#    define FORT_GETDRAG  getdrag
#    define FORT_DERERROR    dererror
#    define FORT_DERVISCOSITY derviscosity
#    define FORT_DERTURBVISC derturbvisc
#    define FORT_DERMAGTRACE dermagtrace
#    define FORT_INTEGRATE_RECALESCE integrate_recalesce
#    define FORT_RESET_TEMPERATURE reset_temperature
#  else
#    define FORT_GETSHEAR getshear_
#    define FORT_MAXPRESVEL maxpresvel_ 
#    define FORT_GETDRAG  getdrag_
#    define FORT_DERERROR    dererror_
#    define FORT_DERVISCOSITY derviscosity_
#    define FORT_DERTURBVISC derturbvisc_
#    define FORT_DERMAGTRACE dermagtrace_
#    define FORT_INTEGRATE_RECALESCE integrate_recalesce_
#    define FORT_RESET_TEMPERATURE reset_temperature_
#  endif
#  endif

#include <ArrayLim.H>

extern "C" {

     void FORT_GETSHEAR(
       int* im,
       int* ntensor,
       Real* cellten, ARLIM_P(celltenlo),ARLIM_P(celltenhi),
       Real* vof, ARLIM_P(voflo),ARLIM_P(vofhi),
       Real* vel, ARLIM_P(vello),ARLIM_P(velhi),
       const Real* delta, const Real* xlo,
       Real* tensordata, ARLIM_P(tendtlo),ARLIM_P(tendthi),
       int* iproject,int* onlyscalar,Real* time,
       const int* tilelo, const int* tilehi,
       const int* fablo, const int* fabhi,
       int* bfact,
       int* level,
       int* bcrec,int* ngrow,
       int* nmat);

     void FORT_DERMAGTRACE(
       int* level,
       int* finest_level,
       int* im,  // im=0..nmat-1
       int* ntensor,
       Real* cellten, ARLIM_P(celltenlo),ARLIM_P(celltenhi),
       Real* vof, ARLIM_P(voflo),ARLIM_P(vofhi),
       Real* dest, ARLIM_P(destlo),ARLIM_P(desthi),
       Real* den, ARLIM_P(denlo),ARLIM_P(denhi),
       Real* tensor, ARLIM_P(tensorlo),ARLIM_P(tensorhi),
       Real* vel, ARLIM_P(vello),ARLIM_P(velhi),
       Real* visc, ARLIM_P(visclo),ARLIM_P(vischi),
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

     void FORT_DERERROR(
      Real* data, ARLIM_P(dat_lo), ARLIM_P(dat_hi),
      const int* nvar, const Real* compdat,
      ARLIM_P(clo), ARLIM_P(chi), const int* ncomp,
      const int* lo, const int* hi,const int* ngrow,
      const int* domain_lo, const int* domain_hi,
      const Real* delta, const Real* xlo,
      const Real* time,const Real* dt,const int* bcrec,
      const int* level, const int* grid_no);

     void FORT_GETDRAG(
      int* isweep,
      Real* globalsum,Real* localsum,
      Real* gravity_normalized,
      int* gravdir,
      int* ntenvisco,
      Real* tdata,ARLIM_P(tdatalo),ARLIM_P(tdatahi),
      Real* viscoten,ARLIM_P(viscotenlo),ARLIM_P(viscotenhi),
      Real* den,ARLIM_P(denlo),ARLIM_P(denhi),
      Real* mask,ARLIM_P(masklo),ARLIM_P(maskhi),
      Real* vof,ARLIM_P(voflo),ARLIM_P(vofhi),
      Real* levelpc,ARLIM_P(levelpclo),ARLIM_P(levelpchi),
      Real* vol, ARLIM_P(vollo),ARLIM_P(volhi),
      Real* areax, ARLIM_P(areaxlo),ARLIM_P(areaxhi),
      Real* areay, ARLIM_P(areaylo),ARLIM_P(areayhi),
      Real* areaz, ARLIM_P(areazlo),ARLIM_P(areazhi),
      Real* xface,ARLIM_P(xfacelo),ARLIM_P(xfacehi),
      Real* yface,ARLIM_P(yfacelo),ARLIM_P(yfacehi),
      Real* zface,ARLIM_P(zfacelo),ARLIM_P(zfacehi),
      Real* cvisc,ARLIM_P(cvisclo),ARLIM_P(cvischi),
      int* facevisc_index,
      int* faceheat_index,
      int* ncphys,
      const Real* xlo,const Real* dx,
      Real* solid, ARLIM_P(solidlo),ARLIM_P(solidhi),
      Real* pres, ARLIM_P(preslo),ARLIM_P(preshi),
      Real* vel, ARLIM_P(vello),ARLIM_P(velhi),
      Real* drag, ARLIM_P(draglo),ARLIM_P(draghi),
      const int* tilelo, const int* tilehi, 
      const int* fablo, const int* fabhi,
      int* bfact, 
      int* rzflag,int* bcrec,
      Real* time,
      Real* visc_coef,
      int* ntensor,
      int* ntensorMM,
      int* nmat,
      int* nparts,
      int* nparts_def,
      int* im_solid_map);


     void FORT_INTEGRATE_RECALESCE(
      int* isweep,
      Real* globalsum,Real* localsum,
      int* recalesce_material,
      Real* snew,ARLIM_P(snewlo),ARLIM_P(snewhi),
      Real* mask,ARLIM_P(masklo),ARLIM_P(maskhi),
      Real* vol, ARLIM_P(vollo),ARLIM_P(volhi),
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
      Real* snew,ARLIM_P(snewlo),ARLIM_P(snewhi),
      const Real* xlo,const Real* dx,
      const int* tilelo, const int* tilehi, 
      const int* fablo, const int* fabhi,
      int* bfact, 
      int* velbc,
      Real* time,
      int* nmat,int* ncomp_state,
      int* level,int* finest_level);


     void FORT_MAXPRESVEL(
       Real* minpres,Real* maxpres,Real* maxvel,
       const Real* xlo,const Real* dx,
       Real* mask, ARLIM_P(masklo),ARLIM_P(maskhi),
       Real* vel, ARLIM_P(vello),ARLIM_P(velhi),
       const int* tilelo, const int* tilehi,
       const int* fablo, const int* fabhi,int* bfact);

     void FORT_DERVISCOSITY(
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
      Real* etaL,Real* etaP,Real* etaS,
      Real* polymer_factor,
      Real* vof, ARLIM_P(voflo), ARLIM_P(vofhi),
      Real* dat, ARLIM_P(datlo), ARLIM_P(dathi),
      Real* vel, ARLIM_P(vello), ARLIM_P(velhi),
      Real* eosdata, ARLIM_P(eosdatalo), ARLIM_P(eosdatahi),
      Real* tensor, ARLIM_P(tensorlo), ARLIM_P(tensorhi),
      Real* gammadot, ARLIM_P(gammadotlo), ARLIM_P(gammadothi),
      const int* tilelo, const int* tilehi,
      const int* fablo, const int* fabhi,
      int* bfact,
      Real* time,
      const Real* dx, const Real* xlo,
      const int* bcrec,int* ngrow,
      int* ncompvisc);


     void FORT_DERTURBVISC(
      int* level,
      int* fortran_im,
      int* nmat,
      Real* dt,
      int* ntensor,
      Real* denstate, ARLIM_P(denstatelo), ARLIM_P(denstatehi),
      Real* vof, ARLIM_P(voflo), ARLIM_P(vofhi),
      Real* vel, ARLIM_P(vello), ARLIM_P(velhi),
      Real* visc, ARLIM_P(visclo), ARLIM_P(vischi),
      Real* cellten, ARLIM_P(celltenlo), ARLIM_P(celltenhi),
      const int* tilelo, const int* tilehi,
      const int* fablo, const int* fabhi,
      int* bfact,
      Real* cur_time,
      const Real* dx, const Real* xlo,
      int* ngrow,
      int* ncompvisc);

}
#endif
#endif
