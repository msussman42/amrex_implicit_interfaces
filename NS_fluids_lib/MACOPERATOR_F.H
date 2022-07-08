#ifndef _MACOPERATOR_F_H_
#define _MACOPERATOR_F_H_

#ifdef BL_LANG_FORT

#else

#include <AMReX_ArrayLim.H>

namespace amrex{

extern "C"
{


    void fort_regularize_bx(
     int* nsolve,
     Real* bx,AMREX_ARLIM_P(bxlo),AMREX_ARLIM_P(bxhi),
     Real* min_interior_coeff,
     const int* domlo, const int* domhi,
     const int* tilelo, const int* tilehi,
     const int* fablo, const int* fabhi,
     int* bfact,
     int* level,
     const Real* xlo,const Real* dx,
     int* dir);

   void fort_init_mask_sing (
     int* level,
     int* finest_level,
     int* nsolve,
     int* project_option,
     Real* masksolv,AMREX_ARLIM_P(masksolvlo),AMREX_ARLIM_P(masksolvhi),
     Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi),
     Real* alpha,
     AMREX_ARLIM_P(alphalo),AMREX_ARLIM_P(alphahi),
     Real* offdiagcheck, 
     AMREX_ARLIM_P(offdiagchecklo),AMREX_ARLIM_P(offdiagcheckhi),
     Real* maskdivres,AMREX_ARLIM_P(maskdivreslo),AMREX_ARLIM_P(maskdivreshi),
     Real* maskres,AMREX_ARLIM_P(maskreslo),AMREX_ARLIM_P(maskreshi),
     Real* mdot,AMREX_ARLIM_P(mdotlo),AMREX_ARLIM_P(mdothi),
     Real* bx,AMREX_ARLIM_P(bxlo),AMREX_ARLIM_P(bxhi),
     Real* by,AMREX_ARLIM_P(bylo),AMREX_ARLIM_P(byhi),
     Real* bz,AMREX_ARLIM_P(bzlo),AMREX_ARLIM_P(bzhi),
     Real* fwtx,AMREX_ARLIM_P(fwtxlo),AMREX_ARLIM_P(fwtxhi),
     Real* fwty,AMREX_ARLIM_P(fwtylo),AMREX_ARLIM_P(fwtyhi),
     Real* fwtz,AMREX_ARLIM_P(fwtzlo),AMREX_ARLIM_P(fwtzhi),
     const int* tilelo, const int* tilehi,
     const int* fablo, const int* fabhi,
     int* bfact,
     int* bc);

   void fort_nsgenerate(
     int* level,
     int* finest_level,
     int* nsolve,
     int* project_option,
     Real* alpha,
     AMREX_ARLIM_P(alphalo),AMREX_ARLIM_P(alphahi),
     Real* diag_reg, 
     AMREX_ARLIM_P(diag_reglo),AMREX_ARLIM_P(diag_reghi),
     Real* bx,AMREX_ARLIM_P(bxlo),AMREX_ARLIM_P(bxhi),
     Real* by,AMREX_ARLIM_P(bylo),AMREX_ARLIM_P(byhi),
     Real* bz,AMREX_ARLIM_P(bzlo),AMREX_ARLIM_P(bzhi),
     const int* tilelo, const int* tilehi,
     const int* fablo, const int* fabhi,
     int* bfact);

    void fort_dividedx(
     int* nsolve,
     Real* bx,AMREX_ARLIM_P(bxlo),AMREX_ARLIM_P(bxhi),
     const int* tilelo, const int* tilehi,
     const int* fablo, const int* fabhi,
     int* bfact,
     int* level,
     const Real* xlo,const Real* dx,
     int* dir);


    void fort_mult_facewt(
     int* nsolve,
     Real* bx,AMREX_ARLIM_P(bxlo),AMREX_ARLIM_P(bxhi),
     Real* facewt,AMREX_ARLIM_P(facewtlo),AMREX_ARLIM_P(facewthi),
     const int* tilelo, const int* tilehi,
     const int* fablo, const int* fabhi,
     int* bfact,
     int* level,
     const Real* xlo,const Real* dx,
     int* dir);

    void fort_interpmac(
        int* bfact,int* bfact_f,
        Real* fine,AMREX_ARLIM_P(fine_lo),AMREX_ARLIM_P(fine_hi),
        const Real* crse,AMREX_ARLIM_P(crse_lo),AMREX_ARLIM_P(crse_hi),
        const int *tlo, const int *thi,
        Real* cdiag,AMREX_ARLIM_P(cdiag_lo),AMREX_ARLIM_P(cdiag_hi));

    void fort_scalarcoeff(
       int* nsolve,
       const Real* xlo,
       const Real* dx,
       Real* offdiagcheck,
       AMREX_ARLIM_P(offdiagchecklo),AMREX_ARLIM_P(offdiagcheckhi),
       Real* cterm,AMREX_ARLIM_P(ctermlo),AMREX_ARLIM_P(ctermhi),
       Real* c2,AMREX_ARLIM_P(c2lo),AMREX_ARLIM_P(c2hi),
       Real* DeDT,AMREX_ARLIM_P(DeDTlo),AMREX_ARLIM_P(DeDThi),
       Real* lsnew,AMREX_ARLIM_P(lsnewlo),AMREX_ARLIM_P(lsnewhi),
       Real* den,AMREX_ARLIM_P(denlo),AMREX_ARLIM_P(denhi),
       Real* mu,AMREX_ARLIM_P(mulo),AMREX_ARLIM_P(muhi),
       const int *tilelo, const int *tilehi,
       const int *fablo, const int *fabhi,
       int* bfact,
       int* level,
       int* finest_level,
       Real* visc_coef,
       Real* angular_velocity,
       Real* dt,
       Real* cur_time,
       int* project_option,
       int* rzflag,
       int* solidheat_flag);

    void fort_restore_pres(
        Real* offdiagcheck,
	AMREX_ARLIM_P(offdiagchecklo),AMREX_ARLIM_P(offdiagcheckhi),
        Real* savepres,AMREX_ARLIM_P(savepreslo),AMREX_ARLIM_P(savepreshi),
        Real* newpres,AMREX_ARLIM_P(newpreslo),AMREX_ARLIM_P(newpreshi),
        const int *tilelo, const int *tilehi,
        const int *fablo, const int *fabhi,
        int* bfact,
        int* level,
        int* finest_level);


}

}/* namespace amrex */

#endif

#endif 
