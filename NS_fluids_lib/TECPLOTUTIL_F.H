#ifndef _TECPLOTUTIL_F_H_
#define _TECPLOTUTIL_F_H_

#ifdef BL_LANG_FORT
#    define FORT_TECPLOTFAB  tecplotfab
#    define FORT_CELLGRID_SANITY cellgrid_sanity
#    define FORT_COMBINEZONES_SANITY combinezones_sanity
#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_TECPLOTFAB  TECPLOTFAB
#    define FORT_CELLGRID_SANITY CELLGRID_SANITY
#    define FORT_COMBINEZONES_SANITY COMBINEZONES_SANITY
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_TECPLOTFAB  tecplotfab
#    define FORT_CELLGRID_SANITY cellgrid_sanity
#    define FORT_COMBINEZONES_SANITY combinezones_sanity
#  else
#    define FORT_TECPLOTFAB  tecplotfab_
#    define FORT_CELLGRID_SANITY cellgrid_sanity_
#    define FORT_COMBINEZONES_SANITY combinezones_sanity_
#  endif
#  endif

#include <AMReX_ArrayLim.H>

namespace amrex{

extern "C" {

    void FORT_TECPLOTFAB(
     Real* time,
     Real* fabdata,AMREX_ARLIM_P(fabdatalo),AMREX_ARLIM_P(fabdatahi),
     const int* growlo,const int* growhi,
     const int* fablo,const int* fabhi,
     int* bfact,
     const Real* xlo,const Real* dx,
     int* dir,int* ncomp,
     int* interior_only,int* nsteps);


    void FORT_CELLGRID_SANITY(
     int* tid,
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

    void FORT_COMBINEZONES_SANITY(
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
     int* visual_option,
     int* visual_revolve,
     int* ncomp);
}
}/* namespace amrex */
    
#endif
#endif /*_TECPLOTUTIL_F_H_*/
