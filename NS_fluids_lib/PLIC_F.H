#ifndef _PLIC_F_H_
#define _PLIC_F_H_

#ifdef BL_LANG_FORT

#else

#include <AMReX_ArrayLim.H>

namespace amrex{

extern "C"
{

      void fort_MOF_DT_training(
       int* num_samples,
       int* finest_level,
       int* bfact,
       const int* domlo,const int* domhi,
       const Real* dx);

      void fort_MOF_training(
       int* num_samples,
       int* op_training,
       int* cpp_training_lo,
       int* cpp_training_hi,
       int* i,
       int* j,
       int* k,
       int* finest_level,
       int* bfact,
       const int* domlo,const int* domhi,
       const Real* dx, 
       int* continuous_mof);

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
        const Real* xlo,
	const Real* dx,
        Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi),
        Real* masknbr,AMREX_ARLIM_P(masknbrlo),AMREX_ARLIM_P(masknbrhi),
        Real* snew,AMREX_ARLIM_P(snewlo),AMREX_ARLIM_P(snewhi),
        Real* vof,AMREX_ARLIM_P(voflo),AMREX_ARLIM_P(vofhi),
        Real* LS,AMREX_ARLIM_P(LSlo),AMREX_ARLIM_P(LShi),
        Real* recon,AMREX_ARLIM_P(reconlo),AMREX_ARLIM_P(reconhi),
        int* nsteps,
        Real* time,
        int* update_flag,
	int* number_centroid_per_core,
	Real* delta_centroid_per_core,
        int* total_calls,
        int* total_iterations,
        Real* total_errors,
        int* continuous_mof,
        int* partial_cmof_stencil_at_walls,
        Real* growth_angle, 
        int* growth_angle_primary_mat, 
        int* growth_angle_tertiary_mat, 
        int* growth_angle_secondary_mat);

}
}/* namespace amrex */
#endif

#endif 
