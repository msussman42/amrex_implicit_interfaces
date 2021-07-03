#ifndef _PLIC_F_H_
#define _PLIC_F_H_

#ifdef BL_LANG_FORT

#else

#include <AMReX_ArrayLim.H>

namespace amrex{

extern "C"
{
      void fort_sloperecon(
        int* tid,
        const int* gridno,
        int* level,
        int* finest_level,
        int* max_level,
        int* ngrow,
        int* vofbc,
        const int* tilelo,const int* tilehi,
        const int* fablo,const int* fabhi,
        int* bfact,
        const Real* xlo,const Real* dx,
        Real* masknbr,AMREX_ARLIM_P(masknbrlo),AMREX_ARLIM_P(masknbrhi),
        Real* snew,AMREX_ARLIM_P(snewlo),AMREX_ARLIM_P(snewhi),
        Real* vof,AMREX_ARLIM_P(voflo),AMREX_ARLIM_P(vofhi),
        Real* LS,AMREX_ARLIM_P(LSlo),AMREX_ARLIM_P(LShi),
        Real* recon,AMREX_ARLIM_P(reconlo),AMREX_ARLIM_P(reconhi),
        int* nsteps,
        Real* time,
        int* nmat,int* nten,
        Real* latent_heat,int* update_flag,
        int* total_calls,
        int* total_iterations,
        int* continuous_mof,
        int* force_cmof_at_triple_junctions, 
        int* partial_cmof_stencil_at_walls, 
        int* radius_cutoff);

}
}/* namespace amrex */
#endif

#endif 
