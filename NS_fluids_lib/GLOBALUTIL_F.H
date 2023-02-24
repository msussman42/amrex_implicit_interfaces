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

   int project_option_needs_scalingF(int* project_option);
   int project_option_projectionF(int* project_option);
   int project_option_FSI_rigid(int* project_option);
   int project_option_is_static(int* project_option);
   int project_option_singular_possibleF(int* project_option);

   int fort_CTML_FSI_mat_base(
        int* FSI_flag_local,
	int* im);  // 1<=im<=num_materials

   int fort_is_FSI_rigid_base(
        int* FSI_flag_local,
	int* im);  // 1<=im<=num_materials

   int fort_is_rigid_base(
        int* FSI_flag_local,
	int* im);  // 1<=im<=num_materials

   int fort_is_lag_part_base(
        int* FSI_flag_local,
	int* im);  // 1<=im<=num_materials

   int fort_is_ice_base(
        int* FSI_flag_local,
	int* im);  // 1<=im<=num_materials

   int fort_FSI_flag_valid_base(
        int* FSI_flag_local,
	int* im);  // 1<=im<=num_materials

   int fort_drag_type(
    int* drag_comp,//0..NDRAG-1
    int* drag_im); //0..num_materials-1

   int fort_drag_IQ_type(
    int* drag_comp,//0..NDRAG_IQ-1
    int* drag_im); //0..num_materials-1

   int is_ice_or_FSI_rigid_material(int* im); //im=1..num_materials

   int fort_read_from_CAD(
    int* fsi_flag_local);

   void get_use_DCA(int* use_DCA);

   void fort_check_operation_flag_MAC(int* operation_flag);

   void fort_derive_gravity_dir(Real* gravity_vector_in,int* gravity_dir);

}

}/* namespace amrex */

#endif

#endif 
