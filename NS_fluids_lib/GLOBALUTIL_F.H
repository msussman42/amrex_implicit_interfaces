#ifndef _GLOBALUTIL_F_H_
#define _GLOBALUTIL_F_H_

#ifdef BL_LANG_FORT

#else

#include <AMReX_ArrayLim.H>

namespace amrex{

extern "C"
{

   Real get_user_latent_heat(
     int iten, //pass by value
     Real temperature, //pass by value
     int default_flag); //pass by value

   int fort_built_in_elastic_model(
    Real* elastic_visc_in,
    int* viscoelastic_model_in);

   int fort_is_passive_advect_test();

   int fort_is_rigid_base(
        int* FSI_flag_local,
	int* nmat,
	int* im);  // 1<=im<=nmat

   int fort_drag_type(
    int* drag_comp,//0..NDRAG-1
    int* drag_im); //0..nmat-1

   int fort_drag_IQ_type(
    int* drag_comp,//0..NDRAG_IQ-1
    int* drag_im); //0..nmat-1

   int fort_read_from_CAD(
    int* fsi_flag_local);

   void fort_check_operation_flag_MAC(int* operation_flag);

}

}/* namespace amrex */

#endif

#endif 
