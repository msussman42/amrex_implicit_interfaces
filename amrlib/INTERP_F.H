#ifndef _INTERP_NAMES_H
#define _INTERP_NAMES_H

#ifdef BL_LANG_FORT
#else

#include <AMReX_ArrayLim.H>

namespace amrex{

extern "C"
{

    void fort_override_finest_level(int* cc_finest_level);

    void fort_gl_slab(
      Real* time_array,
      int* slab_dt_type,
      int* cc_time_order,
      Real* slablow,Real* slabhigh);

    void fort_refine_density_interp (
     const Real* crse, 
     AMREX_ARLIM_P(clo), AMREX_ARLIM_P(chi),
     const int* crse_bx_lo, 
     const int* crse_bx_hi,
     Real* fine, 
     AMREX_ARLIM_P(flo), AMREX_ARLIM_P(fhi),
     const int* fblo, const int* fbhi,
     const Real* prob_lo,
     const Real* dxf, const Real* dxc,
     const int* nvar, 
     int* levelc,int* levelf,
     int* bfactc,int* bfactf);


    void fort_pcinterp (
     int* grid_type,
     int* zapflag,
     const Real* crse, 
     AMREX_ARLIM_P(clo), AMREX_ARLIM_P(chi),
     const int* crse_bx_lo, 
     const int* crse_bx_hi,
     Real* fine, 
     AMREX_ARLIM_P(flo), AMREX_ARLIM_P(fhi),
     const int* fblo, const int* fbhi,
     const Real* prob_lo,
     const Real* dxf, const Real* dxc,
     const int* nvar, 
     int* levelc,int* levelf,
     int* bfactc,int* bfactf);

    // enable_spectral:
    // 0 - low order
    // 1 - space/time spectral
    // 2 - space spectral only
    // 3 - time spectral only
    void fort_seminterp (
     int* enable_spectral,
     const Real* dxc,const Real* dxf,
     const Real* crse, AMREX_ARLIM_P(clo), AMREX_ARLIM_P(chi),
     Real*  fine, AMREX_ARLIM_P(flo), AMREX_ARLIM_P(fhi),
     const int* fblo, const int* fbhi,
     const int* nvar,
     int* levelc,int* levelf,
     int* bfactc,int* bfactf);

    void fort_maskinterppc (
     const Real* crse, AMREX_ARLIM_P(clo), AMREX_ARLIM_P(chi),
     Real*  fine, AMREX_ARLIM_P(flo), AMREX_ARLIM_P(fhi),
     const int* fblo, const int* fbhi,
     const int* nvar, 
     int* levelc,int* levelf,
     int* bfactc,int* bfactf);

    void fort_lsinterp(
       const Real* crse, AMREX_ARLIM_P(clo), AMREX_ARLIM_P(chi),
       const int* clo, const int* chi,
       Real*  fine, AMREX_ARLIM_P(flo), AMREX_ARLIM_P(fhi),
       const int* fblo, const int* fbhi,
       const Real* prob_lo,
       const Real* dxf, const Real* dxc,
       int* ncomp,
       int* levelc,int* levelf,
       int* bfactc,int* bfactf);

    void fort_multimofinterp(
       Real* time,
       const Real* crse, AMREX_ARLIM_P(clo), AMREX_ARLIM_P(chi),
       const int* clo, const int* chi,
       Real*  fine, AMREX_ARLIM_P(flo), AMREX_ARLIM_P(fhi),
       const int* fblo, const int* fbhi,
       Real*  recon, AMREX_ARLIM_P(reconlo), AMREX_ARLIM_P(reconhi),
       const Real* prob_lo,
       const Real* dxf, const Real* dxc,
       int* ngeom_recon,int* ngeom_raw,
       int* levelc,int* levelf,
       int* bfactc,int* bfactf);


    void fort_multiextmofinterp(
       Real* time,
       const Real* crse, AMREX_ARLIM_P(clo), AMREX_ARLIM_P(chi),
       Real*  fine, AMREX_ARLIM_P(flo), AMREX_ARLIM_P(fhi),
       const int* fblo, const int* fbhi,
       const Real* prob_lo,
       const Real* dxf, const Real* dxc,
       int* ngeom_recon,int* ngeom_raw,
       int* levelc,int* levelf,
       int* bfactc,int* bfactf);

    void fort_ext_burnvel_interp(
       int* velflag,
       Real* time,
       const Real* crse, AMREX_ARLIM_P(clo), AMREX_ARLIM_P(chi),
       Real*  fine, AMREX_ARLIM_P(flo), AMREX_ARLIM_P(fhi),
       const int* fblo, const int* fbhi,
       const Real* prob_lo,
       const Real* dxf, const Real* dxc,
       int* nburning,
       int* levelc,int* levelf,
       int* bfactc,int* bfactf);


    // enable_spectral:
    // 0 - low order
    // 1 - space/time spectral
    // 2 - space spectral only
    // 3 - time spectral only
    void fort_edgeinterp(
      int* enable_spectral,
      int* dir,
      const Real* crse, AMREX_ARLIM_P(clo), AMREX_ARLIM_P(chi),
      const int* cloMAC, const int* chiMAC,
      Real*  fine, AMREX_ARLIM_P(flo), AMREX_ARLIM_P(fhi),
      const int* fblo, const int* fbhi,
      const Real* prob_lo,
      const Real* dxf, const Real* dxc,
      const int* ncomp,
      int* levelc,int* levelf,
      int* bfactc,int* bfactf);
}

}/* namespace amrex */

#endif

#endif /*_INTERP_NAMES_H*/
