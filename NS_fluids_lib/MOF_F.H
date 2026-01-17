#ifndef _MOF_F_H_
#define _MOF_F_H_ 

#ifdef BL_LANG_FORT

#else

namespace amrex{

  extern "C" {

    void fort_initmof(
     int * order_algorithm_in,
     int * renormalize_order_algorithm_in,
     Real* denconst_local,
     int* is_rigid_local,
     int* MOFITERMAX,
     int* MOFITERMAX_AFTER_PREDICT,
     int* MOF_DEBUG_RECON_in,
     int* MOF_TURN_OFF_LS_in,
     int* nthreads,
     int* nmax_in);
  }
}/* namespace amrex */

#endif

#endif

