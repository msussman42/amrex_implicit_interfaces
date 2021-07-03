
#ifndef _LO_F_H_
#define _LO_F_H_

#include <AMReX_REAL.H>

#ifdef BL_LANG_FORT
#else

#include <AMReX_ArrayLim.H>

namespace amrex{

extern "C"
{

    void fort_residl(
     int* level,
     int* mg_coarsest_level,
     int* nsolve,
     const Real* masksing , 
     AMREX_ARLIM_P(masksing_lo), AMREX_ARLIM_P(masksing_hi),
     Real* res      , AMREX_ARLIM_P(res_lo), AMREX_ARLIM_P(res_hi),
     const Real* rhs, AMREX_ARLIM_P(rhs_lo), AMREX_ARLIM_P(rhs_hi),
     const Real* phi, AMREX_ARLIM_P(phi_lo), AMREX_ARLIM_P(phi_hi),
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,
     int* bfact,int* bfact_top);

    void fort_buildmat(
        int* level,
        int* veldir,
        int* nsolve,
        int* isweep,
	Real* solvemask, 
        AMREX_ARLIM_P(solvemasklo), AMREX_ARLIM_P(solvemaskhi),
	Real* a, 
        AMREX_ARLIM_P(alo), AMREX_ARLIM_P(ahi),
	Real* bx, AMREX_ARLIM_P(bxlo), AMREX_ARLIM_P(bxhi),
	Real* by, AMREX_ARLIM_P(bylo), AMREX_ARLIM_P(byhi),
	Real* bz, AMREX_ARLIM_P(bzlo), AMREX_ARLIM_P(bzhi),
        Real* diagfab,
        AMREX_ARLIM_P(worklo), AMREX_ARLIM_P(workhi),
        Real* bxleft,Real* bxright,
        Real* byleft,Real* byright,
        Real* bzleft,Real* bzright,
        Real* icbx,Real* icby,Real* icbz,
	Real* icdiag,
        Real* icdiagrb,
        Real* mask,
        AMREX_ARLIM_P(masklo), AMREX_ARLIM_P(maskhi),
	const int* tilelo, const int* tilehi,
        const int* fablo, const int* fabhi,
        int* bfact,int* bfact_top);

    void fort_applybc(
        int* nsolve,
	Real *phi, AMREX_ARLIM_P(phi_lo), AMREX_ARLIM_P(phi_hi),
	Real *bfab, AMREX_ARLIM_P(bfablo), AMREX_ARLIM_P(bfabhi),
	Real *mfab, AMREX_ARLIM_P(mfablo), AMREX_ARLIM_P(mfabhi),
	int* bcpres,
        const int *tilelo, const int *tilehi,
        const int *fablo, const int *fabhi,
        int* bfact,int* bfact_top);

    void fort_averagecc(
        int* nsolve,
        int* ncomp_expect,
	Real* crseX,       AMREX_ARLIM_P(crseX_lo), AMREX_ARLIM_P(crseX_hi),
	const Real* fineX, AMREX_ARLIM_P(fineX_lo), AMREX_ARLIM_P(fineX_hi),
	const int *tlo, const int *thi, int* avg,
        int* ngrow,
        int* bfact_coarse,int* bfact_fine,int* bfact_top);

    void fort_averageec(
        int* nsolve,
	Real* crseX,       AMREX_ARLIM_P(crseX_lo), AMREX_ARLIM_P(crseX_hi),
	const Real* fineX, AMREX_ARLIM_P(fineX_lo), AMREX_ARLIM_P(fineX_hi),
	const int *tlo, const int *thi, 
	const int *axis,int* avg,
        int* bfact_coarse,int* bfact_fine,int* bfact_top);
}

}/* namespace amrex */

#endif

#endif /*_LO_F_H_*/
