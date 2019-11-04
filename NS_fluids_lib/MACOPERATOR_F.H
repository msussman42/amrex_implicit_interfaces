#ifndef _MACOPERATOR_F_H_
#define _MACOPERATOR_F_H_

#ifdef BL_LANG_FORT
#    define FORT_NSGENERATE    nsgenerate
#    define FORT_DIVIDEDX      dividedx
#    define FORT_MULT_FACEWT   mult_facewt
#    define FORT_INTERPMAC     interpmac
#    define FORT_SCALARCOEFF   scalarcoeff
#    define FORT_RESTORE_PRES  restore_pres
#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_NSGENERATE    NSGENERATE
#    define FORT_DIVIDEDX      DIVIDEDX
#    define FORT_MULT_FACEWT   MULT_FACEWT
#    define FORT_INTERPMAC     INTERPMAC
#    define FORT_SCALARCOEFF   SCALARCOEFF
#    define FORT_RESTORE_PRES  RESTORE_PRES
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_NSGENERATE    nsgenerate
#    define FORT_DIVIDEDX      dividedx
#    define FORT_MULT_FACEWT   mult_facewt
#    define FORT_INTERPMAC     interpmac
#    define FORT_SCALARCOEFF   scalarcoeff
#    define FORT_RESTORE_PRES  restore_pres
#  else
#    define FORT_NSGENERATE    nsgenerate_
#    define FORT_DIVIDEDX      dividedx_
#    define FORT_MULT_FACEWT   mult_facewt_
#    define FORT_INTERPMAC     interpmac_
#    define FORT_SCALARCOEFF   scalarcoeff_
#    define FORT_RESTORE_PRES  restore_pres_
#  endif
#endif

#include <AMReX_ArrayLim.H>

extern "C"
{
   void FORT_NSGENERATE (
     int* num_materials_face,
     int* level,
     int* finest_level,
     int* nsolve,
     int* nsolveMM,
     int* nsolveMM_FACE,
     int* nfacefrac,
     int* ncellfrac,
     int* nmat,
     int* project_option,
     int* ncphys,
     Real* cellmm,AMREX_ARLIM_P(cellmmlo),AMREX_ARLIM_P(cellmmhi), 
     Real* xfacemm,AMREX_ARLIM_P(xfacemmlo),AMREX_ARLIM_P(xfacemmhi), 
     Real* yfacemm,AMREX_ARLIM_P(yfacemmlo),AMREX_ARLIM_P(yfacemmhi), 
     Real* zfacemm,AMREX_ARLIM_P(zfacemmlo),AMREX_ARLIM_P(zfacemmhi), 
     Real* xface,AMREX_ARLIM_P(xfacelo),AMREX_ARLIM_P(xfacehi), 
     Real* yface,AMREX_ARLIM_P(yfacelo),AMREX_ARLIM_P(yfacehi), 
     Real* zface,AMREX_ARLIM_P(zfacelo),AMREX_ARLIM_P(zfacehi), 
     Real* masksolv,AMREX_ARLIM_P(masksolvlo),AMREX_ARLIM_P(masksolvhi),
     Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi),
     Real* alpha,AMREX_ARLIM_P(alphalo),AMREX_ARLIM_P(alphahi),
     Real* offdiagcheck, 
    AMREX_ARLIM_P(offdiagchecklo),AMREX_ARLIM_P(offdiagcheckhi),
     Real* diag_non_sing, 
    AMREX_ARLIM_P(diag_non_singlo),AMREX_ARLIM_P(diag_non_singhi),
     Real* diag_sing, 
    AMREX_ARLIM_P(diag_singlo),AMREX_ARLIM_P(diag_singhi),
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

    void FORT_DIVIDEDX(
     int* num_materials_face,
     int* nsolve,
     int* nsolveMM,
     Real* bx,AMREX_ARLIM_P(bxlo),AMREX_ARLIM_P(bxhi),
     const int* tilelo, const int* tilehi,
     const int* fablo, const int* fabhi,
     int* bfact,
     int* level,
     const Real* xlo,const Real* dx,
     int* dir);


    void FORT_MULT_FACEWT(
     int* num_materials_face,
     int* nsolve,
     int* nsolveMM,
     int* nsolveMM_FACE,
     Real* bx,AMREX_ARLIM_P(bxlo),AMREX_ARLIM_P(bxhi),
     Real* facewt,AMREX_ARLIM_P(facewtlo),AMREX_ARLIM_P(facewthi),
     const int* tilelo, const int* tilehi,
     const int* fablo, const int* fabhi,
     int* bfact,
     int* level,
     const Real* xlo,const Real* dx,
     int* dir);

    void FORT_INTERPMAC(
        int* bfact,int* bfact_f,
        Real* fine,AMREX_ARLIM_P(fine_lo),AMREX_ARLIM_P(fine_hi),
        const Real* crse,AMREX_ARLIM_P(crse_lo),AMREX_ARLIM_P(crse_hi),
        const int *tlo, const int *thi,
        Real* cdiag,AMREX_ARLIM_P(cdiag_lo),AMREX_ARLIM_P(cdiag_hi),
        Real* fdiag,AMREX_ARLIM_P(fdiag_lo),AMREX_ARLIM_P(fdiag_hi));

    void FORT_SCALARCOEFF(
        int* num_materials_face,
        int* nsolve,
        int* nsolveMM,
        int* nmat,
        const Real* xlo,
        const Real* dx,
        Real* offdiagcheck,
       AMREX_ARLIM_P(offdiagchecklo),AMREX_ARLIM_P(offdiagcheckhi),
        Real* cterm,AMREX_ARLIM_P(ctermlo),AMREX_ARLIM_P(ctermhi),
        Real* c2,AMREX_ARLIM_P(c2lo),AMREX_ARLIM_P(c2hi),
        Real* DeDT,AMREX_ARLIM_P(DeDTlo),AMREX_ARLIM_P(DeDThi),
        Real* recon,AMREX_ARLIM_P(reconlo),AMREX_ARLIM_P(reconhi),
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
        Real* dual_time_stepping_coefficient,
        int* project_option,
        int* rzflag,
        int* solidheat_flag);

    void FORT_RESTORE_PRES(
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
#endif

#endif 
