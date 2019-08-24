#ifndef _MOF_F_H_
#define _MOF_F_H_ 

#ifdef BL_LANG_FORT
#    define FORT_INITMOF initmof
#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_INITMOF INITMOF
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_INITMOF initmof
#  else
#    define FORT_INITMOF initmof_
#  endif
#  endif

  extern "C" {

    void FORT_INITMOF(
     int * order_algorithm_in,
     int* nmat,
     int* MOFITERMAX,
     int* MOF_DEBUG_RECON_in,
     int* MOF_TURN_OFF_LS_in,
     int* nthreads,
     int* nmax_in);
  };
#endif

#endif

