
#ifndef _MG_F_H_
#define _MG_F_H_

#ifdef BL_LANG_FORT
#define FORT_AVERAGE   average_for_mg
#define FORT_INTERP    interp_for_mg
#define FORT_CHECKERBOARD_RB checkerboard_rb
#else
#  ifdef BL_FORT_USE_UPPERCASE
#define FORT_AVERAGE   AVERAGE_FOR_MG
#define FORT_INTERP    INTERP_FOR_MG
#define FORT_CHECKERBOARD_RB CHECKERBOARD_RB
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#define FORT_AVERAGE   average_for_mg
#define FORT_INTERP    interp_for_mg
#define FORT_CHECKERBOARD_RB checkerboard_rb
#  else
#define FORT_AVERAGE   average_for_mg_
#define FORT_INTERP    interp_for_mg_
#define FORT_CHECKERBOARD_RB checkerboard_rb_
#endif
#endif

#include <AMReX_ArrayLim.H>

namespace amrex{

extern "C"
{
    void FORT_AVERAGE (
        Real* crse,       AMREX_ARLIM_P(crse_lo), AMREX_ARLIM_P(crse_hi),
        const Real* fine, AMREX_ARLIM_P(fine_lo), AMREX_ARLIM_P(fine_hi),
        const int *tlo, const int *thi,
        int* iaverage, 
        int* bfact_coarse,
        int* bfact_fine,
        int* bfact_top);

    void FORT_INTERP (
        int* bfact,int* bfact_f,int* bfact_top,
        Real* fine,       AMREX_ARLIM_P(fine_lo), AMREX_ARLIM_P(fine_hi),
        const Real* crse, AMREX_ARLIM_P(crse_lo), AMREX_ARLIM_P(crse_hi),
        const int *tlo, const int *thi);

    void FORT_CHECKERBOARD_RB (
      int* ncomp,
      Real* v,  AMREX_ARLIM_P(v_lo),  AMREX_ARLIM_P(v_hi),
      const int* tilelo, const int* tilehi,
      const int* fablo,const int* fabhi,
      int* bfact,int* bfact_top);

}

}/* namespace amrex */

#endif

#endif /*_MG_F_H_*/
