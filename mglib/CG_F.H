
#ifndef _CG_F_H_
#define _CG_F_H_

#include <AMReX_REAL.H>

#ifdef BL_LANG_FORT
#else

#include <AMReX_ArrayLim.H>

namespace amrex{

extern "C"
{
    void fort_cgupdate(
        Real* phi,  AMREX_ARLIM_P(phi_lo),  AMREX_ARLIM_P(phi_hi),
        const Real *aa,
        const Real* tmpy, AMREX_ARLIM_P(tmpy_lo), AMREX_ARLIM_P(tmpy_hi),
        const Real* tmpp, AMREX_ARLIM_P(tmpp_lo), AMREX_ARLIM_P(tmpp_hi),
        const int* tilelo, const int* tilehi,
        const int* fablo,const int* fabhi,
        int* bfact,int* bfact_top);

    void fort_cgadvcp(
        Real* tmpp,AMREX_ARLIM_P(tmpp_lo), AMREX_ARLIM_P(tmpp_hi),
        const Real* tmpr, AMREX_ARLIM_P(tmpr_lo), AMREX_ARLIM_P(tmpr_hi),
        const Real* tmpy, AMREX_ARLIM_P(tmpy_lo), AMREX_ARLIM_P(tmpy_hi),
        const Real *b,
        const int* tilelo, const int* tilehi,
        const int* fablo,const int* fabhi,
        int* bfact,int* bfact_top);

    void fort_cgxdoty(
      int* ncomp,
      Real *pw,
      const Real* phi,  AMREX_ARLIM_P(phi_lo),  AMREX_ARLIM_P(phi_hi),
      const Real* tmpw, AMREX_ARLIM_P(tmpw_lo), AMREX_ARLIM_P(tmpw_hi),
      const int* tilelo, const int* tilehi,
      const int* fablo,const int* fabhi,
      int* bfact,int* bfact_top);


}

}/* namespace amrex */

#endif

#endif /*_CG_F_H_*/
