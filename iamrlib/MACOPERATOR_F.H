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

#include <ArrayLim.H>

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
     Real* cellmm,ARLIM_P(cellmmlo),ARLIM_P(cellmmhi), 
     Real* xfacemm,ARLIM_P(xfacemmlo),ARLIM_P(xfacemmhi), 
     Real* yfacemm,ARLIM_P(yfacemmlo),ARLIM_P(yfacemmhi), 
     Real* zfacemm,ARLIM_P(zfacemmlo),ARLIM_P(zfacemmhi), 
     Real* xface,ARLIM_P(xfacelo),ARLIM_P(xfacehi), 
     Real* yface,ARLIM_P(yfacelo),ARLIM_P(yfacehi), 
     Real* zface,ARLIM_P(zfacelo),ARLIM_P(zfacehi), 
     Real* masksolv, ARLIM_P(masksolvlo),ARLIM_P(masksolvhi),
     Real* maskcov, ARLIM_P(maskcovlo),ARLIM_P(maskcovhi),
     Real* alpha, ARLIM_P(alphalo),ARLIM_P(alphahi),
     Real* offdiagcheck, 
     ARLIM_P(offdiagchecklo),ARLIM_P(offdiagcheckhi),
     Real* diag_non_sing, 
     ARLIM_P(diag_non_singlo),ARLIM_P(diag_non_singhi),
     Real* diag_sing, 
     ARLIM_P(diag_singlo),ARLIM_P(diag_singhi),
     Real* maskdivres, ARLIM_P(maskdivreslo),ARLIM_P(maskdivreshi),
     Real* maskres, ARLIM_P(maskreslo),ARLIM_P(maskreshi),
     Real* bx, ARLIM_P(bxlo),ARLIM_P(bxhi),
     Real* by, ARLIM_P(bylo),ARLIM_P(byhi),
     Real* bz, ARLIM_P(bzlo),ARLIM_P(bzhi),
     Real* fwtx, ARLIM_P(fwtxlo),ARLIM_P(fwtxhi),
     Real* fwty, ARLIM_P(fwtylo),ARLIM_P(fwtyhi),
     Real* fwtz, ARLIM_P(fwtzlo),ARLIM_P(fwtzhi),
     const int* tilelo, const int* tilehi,
     const int* fablo, const int* fabhi,
     int* bfact,
     int* bc);

    void FORT_DIVIDEDX(
     int* num_materials_face,
     int* nsolve,
     int* nsolveMM,
     Real* bx, ARLIM_P(bxlo), ARLIM_P(bxhi),
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
     Real* bx, ARLIM_P(bxlo), ARLIM_P(bxhi),
     Real* facewt, ARLIM_P(facewtlo), ARLIM_P(facewthi),
     const int* tilelo, const int* tilehi,
     const int* fablo, const int* fabhi,
     int* bfact,
     int* level,
     const Real* xlo,const Real* dx,
     int* dir);

    void FORT_INTERPMAC(
        int* bfact,int* bfact_f,
        Real* fine,ARLIM_P(fine_lo), ARLIM_P(fine_hi),
        const Real* crse, ARLIM_P(crse_lo), ARLIM_P(crse_hi),
        const int *tlo, const int *thi,
        Real* cdiag,ARLIM_P(cdiag_lo), ARLIM_P(cdiag_hi),
        Real* fdiag,ARLIM_P(fdiag_lo), ARLIM_P(fdiag_hi));

    void FORT_SCALARCOEFF(
        int* num_materials_face,
        int* nsolve,
        int* nsolveMM,
        int* nmat,
        const Real* xlo,
        const Real* dx,
        Real* offdiagcheck,
        ARLIM_P(offdiagchecklo), ARLIM_P(offdiagcheckhi),
        Real* cterm,ARLIM_P(ctermlo), ARLIM_P(ctermhi),
        Real* c2,ARLIM_P(c2lo), ARLIM_P(c2hi),
        Real* DeDT,ARLIM_P(DeDTlo), ARLIM_P(DeDThi),
        Real* recon,ARLIM_P(reconlo), ARLIM_P(reconhi),
        Real* lsnew,ARLIM_P(lsnewlo), ARLIM_P(lsnewhi),
        Real* den,ARLIM_P(denlo), ARLIM_P(denhi),
        Real* mu,ARLIM_P(mulo), ARLIM_P(muhi),
        const int *tilelo, const int *tilehi,
        const int *fablo, const int *fabhi,
        int* bfact,
        int* level,
        int* finest_level,
        Real* visc_coef,Real* angular_velocity,
        Real* dt,
        int* project_option,
        int* rzflag,
        int* solidheat_flag);

    void FORT_RESTORE_PRES(
        Real* offdiagcheck,
	ARLIM_P(offdiagchecklo), ARLIM_P(offdiagcheckhi),
        Real* savepres,ARLIM_P(savepreslo), ARLIM_P(savepreshi),
        Real* newpres,ARLIM_P(newpreslo), ARLIM_P(newpreshi),
        const int *tilelo, const int *tilehi,
        const int *fablo, const int *fabhi,
        int* bfact,
        int* level,
        int* finest_level);


}
#endif

#endif 
