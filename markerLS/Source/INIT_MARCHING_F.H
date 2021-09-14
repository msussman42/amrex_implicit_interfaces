#ifndef _INIT_MARCHING_F_H_
#define _INIT_MARCHING_F_H_

#ifdef BL_LANG_FORT
#else

#include <AMReX_ArrayLim.H>

namespace amrex{

extern "C" {


    void fort_init_marching(const Real* xlo,const Real* dx,const int* n_cell);

}
}/* namespace amrex */
    
#endif
#endif /*_INIT_MARCHING_F_H_*/
