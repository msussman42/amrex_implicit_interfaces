#ifndef _SHALLOWWATER_F_H_
#define _SHALLOWWATER_F_H_ 

#ifdef BL_LANG_FORT
#    define FORT_FEED_FORWARD_ADVECT feed_forward_advect
#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_FEED_FORWARD_ADVECT FEED_FORWARD_ADVECT
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_FEED_FORWARD_ADVECT feed_forward_advect
#  else
#    define FORT_FEED_FORWARD_ADVECT feed_forward_advect_
#  endif
#  endif

#include <ArrayLim.H>

  extern "C" {

      void FORT_FEED_FORWARD_ADVECT(
        int* init_fluxes,
        int* nsolveMM_FACE,
        int* dir,
        int* nflux_feed_forward,
        int* nstate,
        int* velbc,
        int* slab_step,
        Real* dt_slab,
        Real* prev_time_slab,
        const Real* xlo,
        const Real* dx,
        Real* maskfine,ARLIM_P(maskfinelo),ARLIM_P(maskfinehi),//1=fine/fine
        Real* maskcov,ARLIM_P(maskcovlo),ARLIM_P(maskcovhi),//1=not covered
        Real* ax,ARLIM_P(axlo),ARLIM_P(axhi), 
        Real* vol,ARLIM_P(vollo),ARLIM_P(volhi), 
        Real* xflux,ARLIM_P(xfluxlo),ARLIM_P(xfluxhi), 
        Real* velfab,ARLIM_P(velfablo),ARLIM_P(velfabhi), 
        Real* lsfab,ARLIM_P(lsfablo),ARLIM_P(lsfabhi), 
        Real* snewfab,ARLIM_P(snewfablo),ARLIM_P(snewfabhi), 
        Real* soldfab,ARLIM_P(soldfablo),ARLIM_P(soldfabhi), 
        Real* lsnewfab,ARLIM_P(lsnewfablo),ARLIM_P(lsnewfabhi), 
        Real* lsoldfab,ARLIM_P(lsoldfablo),ARLIM_P(lsoldfabhi), 
        Real* xmacnew,ARLIM_P(xmacnewlo),ARLIM_P(xmacnewhi), 
        const int* tilelo,const int* tilehi,
        const int* fablo,const int* fabhi,
        int* bfact,
        int* level,
        int* finest_level,
        const int* domlo,
        const int* domhi,
        int* nmat);

  };
#endif

#endif

