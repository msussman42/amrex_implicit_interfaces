
#ifndef _LO_F_H_
#define _LO_F_H_

#include <AMReX_REAL.H>

#ifdef BL_LANG_FORT
#define FORT_AVERAGECC          averagecc
#define FORT_AVERAGEEC          averageec
#define FORT_APPLYBC            applybc
#define FORT_RESIDL             resid
#define FORT_BUILDMAT           buildmat
#else
#  ifdef BL_FORT_USE_UPPERCASE
#define FORT_AVERAGECC          AVERAGECC
#define FORT_AVERAGEEC          AVERAGEEC
#define FORT_APPLYBC            APPLYBC
#define FORT_RESIDL             RESID
#define FORT_BUILDMAT           BUILDMAT
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#define FORT_AVERAGECC          averagecc
#define FORT_AVERAGEEC          averageec
#define FORT_APPLYBC            applybc
#define FORT_RESIDL             resid
#define FORT_BUILDMAT           buildmat
#  else
#define FORT_AVERAGECC          averagecc_
#define FORT_AVERAGEEC          averageec_
#define FORT_APPLYBC            applybc_
#define FORT_RESIDL             resid_
#define FORT_BUILDMAT           buildmat_
#endif
#endif

#include <AMReX_ArrayLim.H>

extern "C"
{

    void FORT_RESIDL (
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

    void FORT_BUILDMAT (
        int* level,
        int* veldir,
        int* nsolve,
        int* isweep,
        Real* offdiag_coeff,
        int* check_for_singular,
        Real* diag_regularization,
	Real* solvemask, 
        AMREX_ARLIM_P(solvemasklo), AMREX_ARLIM_P(solvemaskhi),
	Real* a, 
        AMREX_ARLIM_P(alo), AMREX_ARLIM_P(ahi),
	Real* bx, AMREX_ARLIM_P(bxlo), AMREX_ARLIM_P(bxhi),
	Real* by, AMREX_ARLIM_P(bylo), AMREX_ARLIM_P(byhi),
	Real* bz, AMREX_ARLIM_P(bzlo), AMREX_ARLIM_P(bzhi),
        Real* diagnonsing,
        AMREX_ARLIM_P(worklo), AMREX_ARLIM_P(workhi),
        Real* diagsing,
        Real* bxleft,Real* bxright,
        Real* byleft,Real* byright,
        Real* bzleft,Real* bzright,
        Real* icbx,Real* icby,Real* icbz,Real* icdiag,
        Real* icdiagrb,
        Real* mask,
        AMREX_ARLIM_P(masklo), AMREX_ARLIM_P(maskhi),
	const int* tilelo, const int* tilehi,
        const int* fablo, const int* fabhi,
        int* bfact,int* bfact_top);

    void FORT_APPLYBC(
        int* nsolve,
	Real *phi, AMREX_ARLIM_P(phi_lo), AMREX_ARLIM_P(phi_hi),
	Real *bfab, AMREX_ARLIM_P(bfablo), AMREX_ARLIM_P(bfabhi),
	Real *mfab, AMREX_ARLIM_P(mfablo), AMREX_ARLIM_P(mfabhi),
	int* bcpres,
        const int *tilelo, const int *tilehi,
        const int *fablo, const int *fabhi,
        int* bfact,int* bfact_top);

    void FORT_AVERAGECC (
        int* nsolve,
        int* ncomp_expect,
	Real* crseX,       AMREX_ARLIM_P(crseX_lo), AMREX_ARLIM_P(crseX_hi),
	const Real* fineX, AMREX_ARLIM_P(fineX_lo), AMREX_ARLIM_P(fineX_hi),
	const int *tlo, const int *thi, int* avg,
        int* ngrow,
        int* bfact_coarse,int* bfact_fine,int* bfact_top);

    void FORT_AVERAGEEC (
        int* nsolve,
	Real* crseX,       AMREX_ARLIM_P(crseX_lo), AMREX_ARLIM_P(crseX_hi),
	const Real* fineX, AMREX_ARLIM_P(fineX_lo), AMREX_ARLIM_P(fineX_hi),
	const int *tlo, const int *thi, 
	const int *axis,int* avg,
        int* bfact_coarse,int* bfact_fine,int* bfact_top);
}
#endif

#endif /*_LO_F_H_*/
