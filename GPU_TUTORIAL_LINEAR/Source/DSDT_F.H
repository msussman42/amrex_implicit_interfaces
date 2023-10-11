#include <AMReX_ArrayLim.H>

namespace amrex{

  extern "C" {

    void fort_dsdtfab(
     const Real* dsdtfab,
     AMREX_ARLIM_P(dsdt_fablo),AMREX_ARLIM_P(dsdt_fabhi),
     const Real* s_data_fab,
     AMREX_ARLIM_P(s_data_fablo),AMREX_ARLIM_P(s_data_fabhi),
     const Real* dx,
     const Real* dt,
     const Real* sound_speed,
     const int* lo, const int* hi);

  }

}/* namespace amrex */


