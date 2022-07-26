#ifndef _MOF_F_H_
#define _MOF_F_H_ 

#ifdef BL_LANG_FORT

#else

  extern "C" {

    void fort_initmof(
     int * order_algorithm_in,
     int* MOFITERMAX,
     int* MOF_DEBUG_RECON_in,
     int* MOF_TURN_OFF_LS_in,
     int* nthreads,
     int* nmax_in);
  }
#endif

#endif

