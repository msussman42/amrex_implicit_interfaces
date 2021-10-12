#ifndef _SOLIDFLUID_F_H_
#define _SOLIDFLUID_F_H_ 

#ifdef BL_LANG_FORT

#else

#include <AMReX_ArrayLim.H>

namespace amrex{

  extern "C" {

    void fort_headermsg(
     int* tid,
     int* tilenum,
     const int* gridno,
     int* nthread_parm,
     int* level,
     int* finest_level,
     const int* max_level,
     int* im_critical,
     int* num_nodes_list,
     int* num_elements_list,
     int* FSI_input_num_nodes,
     int* FSI_input_num_elements,
     Real* FSI_input_node_list,
     int* FSI_input_element_list,
     Real* FSI_input_displacement_list,
     Real* FSI_input_velocity_halftime_list,
     Real* FSI_input_velocity_list,
     Real* FSI_input_force_list,
     Real* FSI_input_mass_list,
     Real* FSI_input_temperature_list,
     int* FSI_output_num_nodes,
     int* FSI_output_num_elements,
     Real* FSI_output_node_list,
     int* FSI_output_element_list,
     Real* FSI_output_displacement_list,
     Real* FSI_output_velocity_halftime_list,
     Real* FSI_output_velocity_list,
     Real* FSI_output_force_list,
     Real* FSI_output_mass_list,
     Real* FSI_output_temperature_list,
     int* FSI_operation,
     int* FSI_sub_operation,
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,
     int* bfact,
     const Real* xlo,
     const Real* dx,
     Real* dx_max_level,
     Real* problo,
     Real* probhi,
     const int* velbc,
     const int* vofbc,
     Real* FSIdata,AMREX_ARLIM_P(FSIdatalo),AMREX_ARLIM_P(FSIdatahi),
     Real* veldata,AMREX_ARLIM_P(veldatalo),AMREX_ARLIM_P(veldatahi),
     Real* masknbr,AMREX_ARLIM_P(masknbrlo),AMREX_ARLIM_P(masknbrhi),
     Real* maskfiner,AMREX_ARLIM_P(maskfinerlo),AMREX_ARLIM_P(maskfinerhi),
     int* nFSI,
     int* nFSI_sub,
     int* ngrowFSI,
     int* nparts,
     int* im_solid_map,
     Real* h_small,
     Real* cur_time,
     Real* dt,
     int* FSI_refine_factor,
     int* FSI_bounding_box_ngrow,
     int* touch_flag,
     int* CTML_FSI_INIT,
     int* CTML_force_model,
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
     int* nmat,
     int* nparts,
     int* im_solid_map,
     Real* problo,
     Real* probhi);

  };

}/* namespace amrex */

#endif

#endif

