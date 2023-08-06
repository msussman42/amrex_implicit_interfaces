#ifndef _CTMLFSI_F_H_
#define _CTMLFSI_F_H_ 

#ifdef BL_LANG_FORT

#else

#include <AMReX_ArrayLim.H>

namespace amrex{

extern "C" {
			
 void fort_ctml_max_nodes( 
  const int* num_materials,
  const int* FSI_flag,
  int* CTML_FSI_numsolids,
  int* CTML_FSI_num_scalars,
  int* max_num_nodes_list,
  int* max_num_elements_list);

}

}/* namespace amrex */

#endif

#endif
