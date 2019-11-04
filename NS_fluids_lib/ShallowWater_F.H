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

#include <AMReX_ArrayLim.H>

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
        Real* maskfine,AMREX_ARLIM_P(maskfinelo),AMREX_ARLIM_P(maskfinehi),//1=fine/fine
        Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi),//1=not covered
        Real* ax,AMREX_ARLIM_P(axlo),AMREX_ARLIM_P(axhi), 
        Real* vol,AMREX_ARLIM_P(vollo),AMREX_ARLIM_P(volhi), 
        Real* xflux,AMREX_ARLIM_P(xfluxlo),AMREX_ARLIM_P(xfluxhi), 
        Real* velfab,AMREX_ARLIM_P(velfablo),AMREX_ARLIM_P(velfabhi), 
        Real* lsfab,AMREX_ARLIM_P(lsfablo),AMREX_ARLIM_P(lsfabhi), 
        Real* snewfab,AMREX_ARLIM_P(snewfablo),AMREX_ARLIM_P(snewfabhi), 
        Real* soldfab,AMREX_ARLIM_P(soldfablo),AMREX_ARLIM_P(soldfabhi), 
        Real* lsnewfab,AMREX_ARLIM_P(lsnewfablo),AMREX_ARLIM_P(lsnewfabhi), 
        Real* lsoldfab,AMREX_ARLIM_P(lsoldfablo),AMREX_ARLIM_P(lsoldfabhi), 
        Real* xmacnew,AMREX_ARLIM_P(xmacnewlo),AMREX_ARLIM_P(xmacnewhi), 
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

