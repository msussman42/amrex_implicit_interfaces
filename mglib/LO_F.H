/*
** (c) 1996-2000 The Regents of the University of California (through
** E.O. Lawrence Berkeley National Laboratory), subject to approval by
** the U.S. Department of Energy.  Your use of this software is under
** license -- the license agreement is attached and included in the
** directory as license.txt or you may contact Berkeley Lab's Technology
** Transfer Department at TTD@lbl.gov.  NOTICE OF U.S. GOVERNMENT RIGHTS.
** The Software was developed under funding from the U.S. Government
** which consequently retains certain rights as follows: the
** U.S. Government has been granted for itself and others acting on its
** behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
** Software to reproduce, prepare derivative works, and perform publicly
** and display publicly.  Beginning five (5) years after the date
** permission to assert copyright is obtained from the U.S. Department of
** Energy, and subject to any subsequent five (5) year renewals, the
** U.S. Government is granted for itself and others acting on its behalf
** a paid-up, nonexclusive, irrevocable, worldwide license in the
** Software to reproduce, prepare derivative works, distribute copies to
** the public, perform publicly and display publicly, and to permit
** others to do so.
*/

#ifndef _LO_F_H_
#define _LO_F_H_

/*
** $Id: LO_F.H,v 1.2 1999/08/25 18:04:10 propp Exp $
*/

#include <REAL.H>

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

#include <ArrayLim.H>

extern "C"
{

    void FORT_RESIDL (
        int* nsolve,
	Real* res      , ARLIM_P(res_lo), ARLIM_P(res_hi),
	const Real* rhs, ARLIM_P(rhs_lo), ARLIM_P(rhs_hi),
	const Real* phi, ARLIM_P(phi_lo), ARLIM_P(phi_hi),
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
        ARLIM_P(solvemasklo), ARLIM_P(solvemaskhi),
	Real* a, 
        ARLIM_P(alo), ARLIM_P(ahi),
	Real* bx, ARLIM_P(bxlo), ARLIM_P(bxhi),
	Real* by, ARLIM_P(bylo), ARLIM_P(byhi),
	Real* bz, ARLIM_P(bzlo), ARLIM_P(bzhi),
        Real* diagnonsing,
        ARLIM_P(worklo), ARLIM_P(workhi),
        Real* diagsing,
        Real* bxleft,Real* bxright,
        Real* byleft,Real* byright,
        Real* bzleft,Real* bzright,
        Real* icbx,Real* icby,Real* icbz,Real* icdiag,
        Real* icdiagrb,
        Real* mask,
        ARLIM_P(masklo), ARLIM_P(maskhi),
	const int* tilelo, const int* tilehi,
        const int* fablo, const int* fabhi,
        int* bfact,int* bfact_top);

    void FORT_APPLYBC(
        int* nsolve,
	Real *phi, ARLIM_P(phi_lo), ARLIM_P(phi_hi),
	Real *bfab, ARLIM_P(bfablo), ARLIM_P(bfabhi),
	Real *mfab, ARLIM_P(mfablo), ARLIM_P(mfabhi),
	int* bcpres,
        const int *tilelo, const int *tilehi,
        const int *fablo, const int *fabhi,
        int* bfact,int* bfact_top);

    void FORT_AVERAGECC (
        int* nsolve,
        int* ncomp_expect,
	Real* crseX,       ARLIM_P(crseX_lo), ARLIM_P(crseX_hi),
	const Real* fineX, ARLIM_P(fineX_lo), ARLIM_P(fineX_hi),
	const int *tlo, const int *thi, int* avg,
        int* ngrow,
        int* bfact_coarse,int* bfact_fine,int* bfact_top);

    void FORT_AVERAGEEC (
        int* nsolve,
	Real* crseX,       ARLIM_P(crseX_lo), ARLIM_P(crseX_hi),
	const Real* fineX, ARLIM_P(fineX_lo), ARLIM_P(fineX_hi),
	const int *tlo, const int *thi, 
	const int *axis,int* avg,
        int* bfact_coarse,int* bfact_fine,int* bfact_top);
}
#endif

#endif /*_LO_F_H_*/
