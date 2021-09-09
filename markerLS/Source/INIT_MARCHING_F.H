#ifndef _MARCHING_TETRA_F_H_
#define _MARCHING_TETRA_F_H_

#ifdef BL_LANG_FORT
#else

#include <AMReX_ArrayLim.H>

namespace amrex{

extern "C" {


    void fort_init_marching(const Real* xlo,const Real* dx,const int* n_cell);

}
}/* namespace amrex */
    
#endif
#endif /*_MARCHING_TETRA_F_H_*/
