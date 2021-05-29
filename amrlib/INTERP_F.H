#ifndef _INTERP_NAMES_H
#define _INTERP_NAMES_H

#ifdef BL_LANG_FORT
#    define FORT_OVERRIDE_FINEST_LEVEL override_finest_level
#    define FORT_GL_SLAB     gl_slab
#    define FORT_PCINTERP    pcinterp
#    define FORT_SEMINTERP   seminterp
#    define FORT_MASKINTERPPC maskinterppc
#    define FORT_LSHOINTERP lshointerp
#    define FORT_MULTIMOFINTERP multimofinterp
#    define FORT_MULTIEXTMOFINTERP multiextmofinterp
#    define FORT_EXT_BURNVEL_INTERP ext_burnvel_interp
#    define FORT_EDGEINTERP  edgeinterp
#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_OVERRIDE_FINEST_LEVEL OVERRIDE_FINEST_LEVEL
#    define FORT_GL_SLAB     GL_SLAB
#    define FORT_PCINTERP    PCINTERP
#    define FORT_SEMINTERP   SEMINTERP
#    define FORT_MASKINTERPPC MASKINTERPPC
#    define FORT_LSHOINTERP LSHOINTERP
#    define FORT_MULTIMOFINTERP MULTIMOFINTERP
#    define FORT_MULTIEXTMOFINTERP MULTIEXTMOFINTERP
#    define FORT_EXT_BURNVEL_INTERP EXT_BURNVEL_INTERP
#    define FORT_EDGEINTERP  EDGEINTERP
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_OVERRIDE_FINEST_LEVEL override_finest_level
#    define FORT_GL_SLAB     gl_slab
#    define FORT_PCINTERP    pcinterp
#    define FORT_SEMINTERP   seminterp
#    define FORT_MASKINTERPPC maskinterppc
#    define FORT_LSHOINTERP lshointerp
#    define FORT_MULTIMOFINTERP multimofinterp
#    define FORT_MULTIEXTMOFINTERP multiextmofinterp
#    define FORT_EXT_BURNVEL_INTERP ext_burnvel_interp
#    define FORT_EDGEINTERP  edgeinterp
#  else
#    define FORT_OVERRIDE_FINEST_LEVEL override_finest_level_
#    define FORT_GL_SLAB     gl_slab_
#    define FORT_PCINTERP    pcinterp_
#    define FORT_SEMINTERP   seminterp_
#    define FORT_MASKINTERPPC maskinterppc_
#    define FORT_LSHOINTERP lshointerp_
#    define FORT_MULTIMOFINTERP multimofinterp_
#    define FORT_MULTIEXTMOFINTERP multiextmofinterp_
#    define FORT_EXT_BURNVEL_INTERP ext_burnvel_interp_
#    define FORT_EDGEINTERP  edgeinterp_
#  endif
#  endif

#include <AMReX_ArrayLim.H>

namespace amrex{

extern "C"
{

    void FORT_OVERRIDE_FINEST_LEVEL(int* cc_finest_level);

    void FORT_GL_SLAB(
      Real* time_array,
      int* slab_dt_type,
      int* cc_time_order,
      Real* slablow,Real* slabhigh);

    void FORT_PCINTERP (
     int* grid_type,
     int* zapflag,
     const Real* crse, AMREX_ARLIM_P(clo), AMREX_ARLIM_P(chi),
     const int* crse_bx_lo, 
     const int* crse_bx_hi,
     Real*  fine, AMREX_ARLIM_P(flo), AMREX_ARLIM_P(fhi),
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
    void FORT_SEMINTERP (
     int* enable_spectral,
     const Real* dxc,const Real* dxf,
     const Real* crse, AMREX_ARLIM_P(clo), AMREX_ARLIM_P(chi),
     Real*  fine, AMREX_ARLIM_P(flo), AMREX_ARLIM_P(fhi),
     const int* fblo, const int* fbhi,
     const int* nvar,
     int* levelc,int* levelf,
     int* bfactc,int* bfactf);

    void FORT_MASKINTERPPC (
     const Real* crse, AMREX_ARLIM_P(clo), AMREX_ARLIM_P(chi),
     Real*  fine, AMREX_ARLIM_P(flo), AMREX_ARLIM_P(fhi),
     const int* fblo, const int* fbhi,
     const int* nvar, 
     int* levelc,int* levelf,
     int* bfactc,int* bfactf);

    void FORT_LSHOINTERP(
       int* LSHOInterp_LO,
       const Real* crse, AMREX_ARLIM_P(clo), AMREX_ARLIM_P(chi),
       const int* clo, const int* chi,
       Real*  fine, AMREX_ARLIM_P(flo), AMREX_ARLIM_P(fhi),
       const int* fblo, const int* fbhi,
       const Real* prob_lo,
       const Real* dxf, const Real* dxc,
       int* nmat,
       int* ncomp,
       int* levelc,int* levelf,
       int* bfactc,int* bfactf);

    void FORT_MULTIMOFINTERP(
       Real* time,
       const Real* crse, AMREX_ARLIM_P(clo), AMREX_ARLIM_P(chi),
       const int* clo, const int* chi,
       Real*  fine, AMREX_ARLIM_P(flo), AMREX_ARLIM_P(fhi),
       const int* fblo, const int* fbhi,
       Real*  recon, AMREX_ARLIM_P(reconlo), AMREX_ARLIM_P(reconhi),
       const Real* prob_lo,
       const Real* dxf, const Real* dxc,
       int* nmat,
       int* ngeom_recon,int* ngeom_raw,
       int* levelc,int* levelf,
       int* bfactc,int* bfactf);


    void FORT_MULTIEXTMOFINTERP(
       Real* time,
       const Real* crse, AMREX_ARLIM_P(clo), AMREX_ARLIM_P(chi),
       Real*  fine, AMREX_ARLIM_P(flo), AMREX_ARLIM_P(fhi),
       const int* fblo, const int* fbhi,
       const Real* prob_lo,
       const Real* dxf, const Real* dxc,
       int* nmat,
       int* ngeom_recon,int* ngeom_raw,
       int* levelc,int* levelf,
       int* bfactc,int* bfactf);

    void FORT_EXT_BURNVEL_INTERP(
       int* velflag,
       Real* time,
       const Real* crse, AMREX_ARLIM_P(clo), AMREX_ARLIM_P(chi),
       Real*  fine, AMREX_ARLIM_P(flo), AMREX_ARLIM_P(fhi),
       const int* fblo, const int* fbhi,
       const Real* prob_lo,
       const Real* dxf, const Real* dxc,
       int* nmat,
       int* nten,
       int* nburning,
       int* levelc,int* levelf,
       int* bfactc,int* bfactf);


    // enable_spectral:
    // 0 - low order
    // 1 - space/time spectral
    // 2 - space spectral only
    // 3 - time spectral only
    void FORT_EDGEINTERP(
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
