
#ifndef _MG_F_H_
#define _MG_F_H_

#ifdef BL_LANG_FORT
#else

#include <AMReX_ArrayLim.H>

namespace amrex{

extern "C"
{
    void fort_average(
        Real* crse,       AMREX_ARLIM_P(crse_lo), AMREX_ARLIM_P(crse_hi),
        const Real* fine, AMREX_ARLIM_P(fine_lo), AMREX_ARLIM_P(fine_hi),
        const int *tlo, const int *thi,
        int* iaverage, 
        int* bfact_coarse,
        int* bfact_fine,
        int* bfact_top);

    void fort_interp(
        int* bfact,int* bfact_f,int* bfact_top,
        Real* fine,       AMREX_ARLIM_P(fine_lo), AMREX_ARLIM_P(fine_hi),
        const Real* crse, AMREX_ARLIM_P(crse_lo), AMREX_ARLIM_P(crse_hi),
        const int *tlo, const int *thi);

}

}/* namespace amrex */

#endif

#endif /*_MG_F_H_*/
