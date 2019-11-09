#ifndef _TECPLOTUTIL_F_H_
#define _TECPLOTUTIL_F_H_

#ifdef BL_LANG_FORT
#    define FORT_TECPLOTFAB  tecplotfab
#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_TECPLOTFAB  TECPLOTFAB
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_TECPLOTFAB  tecplotfab
#  else
#    define FORT_TECPLOTFAB  tecplotfab_
#  endif
#  endif

#include <AMReX_ArrayLim.H>

namespace amrex{

extern "C" {

    void FORT_TECPLOTFAB(
     Real* time,
     Real* fabdata,AMREX_ARLIM_P(fabdatalo),AMREX_ARLIM_P(fabdatahi),
     const int* growlo,const int* growhi,
     const int* fablo,const int* fabhi,
     int* bfact,
     const Real* xlo,const Real* dx,
     int* dir,int* ncomp,
     int* interior_only,int* nsteps);

}
}/* namespace amrex */
    
#endif
#endif /*_TECPLOTUTIL_F_H_*/
