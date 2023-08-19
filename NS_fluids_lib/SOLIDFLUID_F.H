#ifndef _SOLIDFLUID_F_H_
#define _SOLIDFLUID_F_H_ 

#ifdef BL_LANG_FORT

#else

#include <AMReX_ArrayLim.H>

namespace amrex{

  extern "C" {

    void fort_init_aux_data(int* ioproc);

    void fort_headermsg(
     int* processor_id,
     int* tid,
     int* tilenum,
     const int* gridno,
     int* nthread_parm,
     int* level,
     int* finest_level,
     const int* max_level,
     Real* FSI_input_flattened,
     Real* FSI_output_flattened,
     int* flatten_size,
     int* local_caller_id,
     int* FSI_operation,
     int* FSI_sub_operation,
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,
     int* bfact,
     const Real* xlo,
     const Real* dx,
     Real* dx_max_level,
     const int* velbc,
     const int* vofbc,
     Real* FSIdata,AMREX_ARLIM_P(FSIdatalo),AMREX_ARLIM_P(FSIdatahi),
     Real* drag,AMREX_ARLIM_P(draglo),AMREX_ARLIM_P(draghi),
     Real* masknbr,AMREX_ARLIM_P(masknbrlo),AMREX_ARLIM_P(masknbrhi),
     Real* maskfiner,AMREX_ARLIM_P(maskfinerlo),AMREX_ARLIM_P(maskfinerhi),
     int* nFSI,
     int* ngrow_make_distance_in,
     int* nparts,
     int* im_solid_map,
     Real* h_small,
     Real* cur_time,
     Real* dt,
     int* FSI_refine_factor,
     int* FSI_bounding_box_ngrow,
     int* touch_flag,
     int* CTML_FSI_INIT,
     int* iter,
     int* current_step,
     int* plot_interval,
     int* ioproc);

    void fort_fillcontainer(
     int* level,
     int* finest_level,
     const int* sci_max_level,
     Real* cur_time,
     Real* dt,
     int* tilelo_array,
     int* tilehi_array,
     Real* xlo_array,
     const Real* dx,
     Real* dx_max_level,
     int* num_grids_on_level,
     int* num_grids_on_level_proc,
     int* gridno_array,
     int* num_tiles_on_thread_proc,
     int* nthread_parm,
     int* max_num_tiles_on_thread_proc,
     int* tile_dim,
     int* nparts,
     int* im_solid_map);

  }

}/* namespace amrex */

#endif

#endif

