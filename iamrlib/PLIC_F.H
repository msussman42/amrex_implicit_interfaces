#ifndef _PLIC_F_H_
#define _PLIC_F_H_

#ifdef BL_LANG_FORT
#    define FORT_SLOPERECON   sloperecon
#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_SLOPERECON   SLOPERECON
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_SLOPERECON   sloperecon
#  else
#    define FORT_SLOPERECON   sloperecon_
#  endif
#  endif

#include <ArrayLim.H>

extern "C"
{
      void FORT_SLOPERECON(
        int* tid,
        const int* gridno,
        int* level,
        int* finest_level,
        int* max_level,
        int* ngrow,
        const int* tilelo,const int* tilehi,
        const int* fablo,const int* fabhi,
        int* bfact,
        const Real* xlo,const Real* dx,
        Real* masknbr,ARLIM_P(masknbrlo),ARLIM_P(masknbrhi),
        Real* snew,ARLIM_P(snewlo),ARLIM_P(snewhi),
        Real* vof,ARLIM_P(voflo),ARLIM_P(vofhi),
        Real* LS,ARLIM_P(LSlo),ARLIM_P(LShi),
        Real* recon,ARLIM_P(reconlo),ARLIM_P(reconhi),
        int* nsteps,
        Real* time,
        int* nmat,int* nten,
        Real* latent_heat,int* update_flag,
        int* total_calls,
        int* total_iterations,
        int* continuous_mof,
        int* radius_cutoff);

}
#endif

#endif 
