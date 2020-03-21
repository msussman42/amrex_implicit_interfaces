#ifndef _AABBEC_F_H_
#define _AABBEC_F_H_

#include <AMReX_REAL.H>

#ifdef BL_LANG_FORT
#define FORT_GSRB          gsrb
#define FORT_ADOTX         adotx
#define FORT_DIAGSUM       diagsum
#else
#  ifdef BL_FORT_USE_UPPERCASE
#define FORT_GSRB     GSRB
#define FORT_ADOTX    ADOTX
#define FORT_DIAGSUM  DIAGSUM
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#define FORT_GSRB     gsrb
#define FORT_ADOTX    adotx
#define FORT_DIAGSUM  diagsum
#  else
#define FORT_GSRB     gsrb_
#define FORT_ADOTX    adotx_
#define FORT_DIAGSUM  diagsum_
#endif
#endif

#include <AMReX_ArrayLim.H>

namespace amrex{

extern "C"
{

    void FORT_GSRB (
     int* level,
     int* mg_coarsest_level,
     int* isweep,
     int* num_sweeps,
     const Real* masksing , 
     AMREX_ARLIM_P(masksing_lo), AMREX_ARLIM_P(masksing_hi),
     Real* phi,AMREX_ARLIM_P(philo), AMREX_ARLIM_P(phihi),
     const Real* rhs, AMREX_ARLIM_P(rhslo), AMREX_ARLIM_P(rhshi),
     const Real* diagfab, 
     AMREX_ARLIM_P(diagfablo),  AMREX_ARLIM_P(diagfabhi),
     const Real* bXleft,
     const Real* bXright,
     const Real* bYleft,
     const Real* bYright,
     const Real* bZleft,
     const Real* bZright,
     const Real* icbx,
     const Real* icby,
     const Real* icbz,
     const Real* icdiag,
     const Real* icdiagrb,
     const Real* mask,
     const Real* ax,
     const Real* solnsave,
     const Real* rhssave,
     const Real* redsoln,
     const Real* blacksoln,
     const int *tilelo, const int *tilehi,
     const int *fablo, const int *fabhi,
     int* bfact,int* bfact_top,
     const int* smooth_type);

    void FORT_ADOTX(
     int* level,
     int* mg_coarsest_level,
     const Real* masksing , 
     AMREX_ARLIM_P(masksing_lo), AMREX_ARLIM_P(masksing_hi),
     Real *y, AMREX_ARLIM_P(ylo), AMREX_ARLIM_P(yhi),
     const Real *x, AMREX_ARLIM_P(xlo), AMREX_ARLIM_P(xhi),
     const Real *diagfab,
     AMREX_ARLIM_P(diagfablo),AMREX_ARLIM_P(diagfabhi),
     const Real* bXleft,
     const Real* bXright,
     const Real* bYleft,
     const Real* bYright,
     const Real* bZleft,
     const Real* bZright,
     const int *tilelo, const int *tilehi,
     const int *fablo, const int *fabhi,
     int* bfact,int* bfact_top);
   
    void FORT_DIAGSUM(
        Real *y      , AMREX_ARLIM_P(y_lo), AMREX_ARLIM_P(y_hi),
        const Real* bX, AMREX_ARLIM_P(bX_lo), AMREX_ARLIM_P(bX_hi),
        const Real* bY, AMREX_ARLIM_P(bY_lo), AMREX_ARLIM_P(bY_hi),
        const Real* bZ, AMREX_ARLIM_P(bZ_lo), AMREX_ARLIM_P(bZ_hi),
        const int *tilelo, const int *tilehi,
        const int *fablo, const int *fabhi,
        int* bfact,int* bfact_top);

}

}/* namespace amrex */

#endif

#endif /*_AABBEC_F_H_*/
