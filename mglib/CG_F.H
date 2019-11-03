
#ifndef _CG_F_H_
#define _CG_F_H_

#include <AMReX_REAL.H>

#ifdef BL_LANG_FORT
#define FORT_CGUPDATE  cgupdate3dgen
#define FORT_CGADVCP   cgadvcp3dgen
#define FORT_CGXDOTY   cgxdoty3dgen
#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_CGUPDATE  CGUPDATE3DGEN
#    define FORT_CGADVCP   CGADVCP3DGEN
#    define FORT_CGXDOTY   CGXDOTY3DGEN
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_CGUPDATE  cgupdate3dgen
#    define FORT_CGADVCP   cgadvcp3dgen
#    define FORT_CGXDOTY   cgxdoty3dgen
#  else
#    define FORT_CGUPDATE  cgupdate3dgen_
#    define FORT_CGADVCP   cgadvcp3dgen_
#    define FORT_CGXDOTY   cgxdoty3dgen_
#endif
#endif

#include <AMReX_ArrayLim.H>

extern "C"
{
    void FORT_CGUPDATE (
        Real* phi,  AMREX_ARLIM_P(phi_lo),  AMREX_ARLIM_P(phi_hi),
        const Real *aa,
        const Real* tmpy, AMREX_ARLIM_P(tmpy_lo), AMREX_ARLIM_P(tmpy_hi),
        const Real* tmpp, AMREX_ARLIM_P(tmpp_lo), AMREX_ARLIM_P(tmpp_hi),
        const int* tilelo, const int* tilehi,
        const int* fablo,const int* fabhi,
        int* bfact,int* bfact_top);

    void FORT_CGADVCP (
        Real* tmpp,AMREX_ARLIM_P(tmpp_lo), AMREX_ARLIM_P(tmpp_hi),
        const Real* tmpr, AMREX_ARLIM_P(tmpr_lo), AMREX_ARLIM_P(tmpr_hi),
        const Real* tmpy, AMREX_ARLIM_P(tmpy_lo), AMREX_ARLIM_P(tmpy_hi),
        const Real *b,
        const int* tilelo, const int* tilehi,
        const int* fablo,const int* fabhi,
        int* bfact,int* bfact_top);

    void FORT_CGXDOTY (
      int* ncomp,
      Real *pw,
      const Real* phi,  AMREX_ARLIM_P(phi_lo),  AMREX_ARLIM_P(phi_hi),
      const Real* tmpw, AMREX_ARLIM_P(tmpw_lo), AMREX_ARLIM_P(tmpw_hi),
      const int* tilelo, const int* tilehi,
      const int* fablo,const int* fabhi,
      int* bfact,int* bfact_top);


}
#endif

#endif /*_CG_F_H_*/
