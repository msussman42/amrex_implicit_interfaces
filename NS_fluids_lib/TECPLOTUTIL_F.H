#ifndef _TECPLOTUTIL_F_H_
#define _TECPLOTUTIL_F_H_

#ifdef BL_LANG_FORT

#else

#include <AMReX_ArrayLim.H>

namespace amrex{

extern "C" {

    void fort_tecplotfab(
     Real* time,
     Real* fabdata,AMREX_ARLIM_P(fabdatalo),AMREX_ARLIM_P(fabdatahi),
     const int* growlo,const int* growhi,
     const int* fablo,const int* fabhi,
     int* bfact,
     const Real* xlo,const Real* dx,
     int* dir,int* ncomp,
     int* interior_only,int* nsteps);


    void fort_tecplotfab_sanity(
     char* root_char_array,
     int* n_root,
     int* data_dir,
     int* bfact,
     const int* fablo,const int* fabhi,
     Real* fabdata,AMREX_ARLIM_P(fabdatalo),AMREX_ARLIM_P(fabdatahi),
     const Real* problo,
     const Real* probhi,
     const Real* dx,
     int* SDC_outer_sweeps,
     int* slab_step,
     int* data_id,
     int* nsteps,
     Real* time,
     int* visual_revolve,
     int* level, 
     int* finest_level, 
     int* ncomp);


    void fort_cellgrid_sanity(
     int* tid,
     int* data_id,
     int* data_dir,
     int* bfact,
     int* ncomp,
     Real* datafab,AMREX_ARLIM_P(datafablo),AMREX_ARLIM_P(datafabhi),
     const Real* problo,
     const Real* probhi,
     const Real* dx,
     const int* lo,const int* hi,
     const int* level,
     const int* finest_level,
     const int* gridno,
     int* rz_flag);

    void fort_combinezones_sanity(
     char* root_char_array,
     int* n_root,
     int* data_dir,
     int* total_number_grids,
     int* grids_per_level_array,
     int* levels_array,
     int* bfact_array,
     int* gridno_array,
     int* gridlo_array,
     int* gridhi_array,
     int* finest_level,
     int* SDC_outer_sweeps,
     int* slab_step,
     int* data_id,
     int* nsteps,
     int* num_levels,
     Real* time,
     int* visual_revolve,
     int* ncomp);
}
}/* namespace amrex */
    
#endif
#endif /*_TECPLOTUTIL_F_H_*/
