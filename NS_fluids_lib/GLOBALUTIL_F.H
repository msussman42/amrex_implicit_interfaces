#ifndef _GLOBALUTIL_F_H_
#define _GLOBALUTIL_F_H_

#ifdef BL_LANG_FORT

#else

#include <AMReX_ArrayLim.H>

namespace amrex{

extern "C"
{

   int fort_is_eulerian_elastic_model(
    Real* elastic_visc_in,
    int* viscoelastic_model_in);

}

}/* namespace amrex */

#endif

#endif 
