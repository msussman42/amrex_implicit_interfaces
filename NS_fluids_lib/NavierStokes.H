#ifndef _NavierStokes_H_
#define _NavierStokes_H_

#define NS_profile_solver 1

#include <local_thread_class.H>
#include <AMReX_Geometry.H>
#include <AmrLevel.H>
#include <LevelBld.H>
#include <ABecLaplacian.H>

#include <CPP_EPS.H>
#include <EXTRAP_COMP.H>
#include <TOWER_MF_CODES.H>

namespace {
    constexpr double NS_PI = 3.14159265358979323846264338327950288;
}

namespace amrex{

// State_Type:
//  vel, pres, num_state_material x num_materials, ngeom_raw x num_materials, error ind
// Umac_Type, Vmac_Type, Wmac_Type:
//  mac vel
// LS_Type:
//  LS 1..num_materials, LS_slope sdim x num_materials
// DIV_Type (1 variable):
//  (a) temporary "pressure" for project_option==11, or
//  (b) -(pnew-pold)/(rho c^2 dt) + dt mdot/vol
// Solid_State_Type:
// velocity
// Tensor_Type:
// Tij x num_materials_viscoelastic + sdim

// Q is stored in the code.  A=Q+I   A is sometimes called the
// "configuration tensor" 
// See: Stewart, Lay, Sussman, Ohta (and papers cited, and papers who
// cite this)
// See also earlier work involving Pushpendra Singh (NJIT).

class NSBld
    :
    public LevelBld
{
    virtual void variableSetUp () override;
    virtual void variableCleanUp () override;
    virtual AmrLevel *operator() () override;
    virtual AmrLevel *operator() (AmrCore&        papa,
                                  int             lev,
                                  const Geometry& level_geom,
                                  const BoxArray& ba,
                                  const DistributionMapping& dm,
                                  Real            time) override;
};

class blobclass
{
public:
   // X=(x,y,z)  U=(u,v,w)
   // velocity=sum alpha_i phi_i(X)  i=1..6 in 3D (i=1..4 in 2D)
   // E=sum || u_k - sum alpha_i phi_i(X_k) ||^2 w_k
   // least squares matrix 
   // Aij=sum phi_i(X_k) dot phi_j(X_k) w_k
   // Bi=sum phi_i(X_k) dot u_k w_k
   // ux,uy,uz fluid_interior,fluid,solid
  Real blob_matrix[3*(2*AMREX_SPACEDIM)*(2*AMREX_SPACEDIM)]; 
  Real blob_RHS[3*(2*AMREX_SPACEDIM)]; 
   // alphax,alphay,alphaz,alphaxy alphaxz alphayz  
   //  fluid interior,fluid,solid
   // u=alphaxy(y-y0)+alphaxz(z-z0)+alphax  
   // v=-alphaxy(x-x0)+alphayz(z-z0)+alphay
   // w=-alphayz(y-y0)-alphaxz(x-x0)+alphaz
   // phix=(1 0 0) phiy=(0 1 0) phiz=(0 0 1)
   // phixy=( (y-y0)  -(x-x0)   0   )
   // phixz=( (z-z0)     0  -(x-x0) )
   // phiyz=(   0    (z-z0) -(y-y0) )
  Real blob_velocity[3*(2*AMREX_SPACEDIM)]; 
   // translational momentum integral rho u dv
   // angular momentum integral (x-x0) cross rho u dv
   // | i    j    k  |
   // |x-x0 y-y0 z-z0|
   // | u    v    w  |  = ((y-y0)w-(z-z0)v)ihat+((z-z0)u-(x-x0)w))jhat+
   //                     ((x-x0)v-(y-y0)u)khat
   // first group: sum phi_i(x_k) dot u_k wmom_k
   // second group: sum phi_i(x_k) dot phi_i(x_k) wmom_k
   // note:
   // momentum_x=alphax sum_k w_k=alphax sum phi_1 dot phi_1 wmom_k
   // momentum_y and momentum_z are analogous.
   // momentum_khat=-alphaxy sum ((xk-x0)^2 wk + (yk-y0)^2) wmomk=
   //               -alphaxy sum phixy dot phixy wmomk
   // momentum_jhat=alphaxz sum (zk-z0)^2 + (xk-x0)^2 wmomk=
   //               alphaxz sum phixz dot phixz wmomk
   // momentum_ihat=-alphayz sum (zk-z0)^2 + (yk-y0)^2 wmomk=
   //               -alphayz sum phiyz dot phiyz wmomk
  Real blob_integral_momentum[2*(2*AMREX_SPACEDIM)];
  Real blob_energy;
   // fluid_interior,fluid,solid
  Real blob_mass_for_velocity[3];
  Real blob_volume;
  Real blob_cell_count;
  Real blob_cellvol_count;
  Real blob_mass;
  Real blob_pressure;
  Real blob_center_integral[AMREX_SPACEDIM];
  Real blob_center_actual[AMREX_SPACEDIM];
  Real blob_perim;  // surface area of the blob
   // surface area of region separating the blob from a given material im.
  Vector<Real> blob_perim_mat;  
   // perimeter of contact line between im, im1, im2
  Vector< Vector<Real> > blob_triple_perim;
  int im;  // 1<=im<=num_materials
};

class NavierStokes
    :
    public AmrLevel 
{
public:

    NavierStokes ();

    NavierStokes (AmrCore&        papa,
                  int             lev,
                  const Geometry& level_geom,
                  const BoxArray& bl,
                  const DistributionMapping& dmap_in,
                  Real            time);

    virtual ~NavierStokes ();

    virtual void restart (AmrCore&      papa,
                          std::istream& is,
			  int old_finest_level,
			  int new_finest_level);
    //
    // Define data descriptors.
    //

    static void variableSetUp ();
    //
    // Cleanup data descriptors at end of run.
    //
    static void variableCleanUp ();
    //
    // Init grid data at problem start-up.
    //
    virtual void initData ();

    static int project_option_olddata_needed(int project_option);
    static int project_option_singular_possible(int project_option);
    static int project_option_pressure(int project_option);
    static int project_option_needs_scaling(int project_option);
    static int project_option_projection(int project_option);
    static int project_option_momeqn(int project_option);
    static int project_option_is_valid(int project_option);

    static int is_multi_component_evap(int freezing_model,
		    int evap_flag,Real latent_heat);
    static int is_GFM_freezing_model(int freezing_model);
    static int is_hydrate_freezing_model(int freezing_model);
    static int is_valid_freezing_model(int freezing_model);

    //
    // Set time levels of state data.
    //
    virtual void setTimeLevel (Real time,Real& dt);

    virtual void init (AmrLevel& old,
      const BoxArray& ba_in,
      const DistributionMapping& dmap_in);

    virtual void init (
      const BoxArray& ba_in,
      const DistributionMapping& dmap_in);

    //
    // Error estimation for regridding.
    //
    virtual void errorEst (TagBoxArray& tb,
      int clearval, int tagval, 
      Real time,int n_error_buf=0,int ngrow=0);

    int pattern_test(const std::string& source,
      const std::string& pattern);

    void deallocate_potential_forceALL();
    void process_potential_forceALL(int potgrad_surface_tension_mask,
      int project_option);
    void process_potential_force_face(int potgrad_surface_tension_mask,
      int project_option);
    void init_gravity_potential();

    void GetDragALL();
    void GetDrag(int isweep);

    void volWgtSum(int isweep,int fast_mode);

    void Number_CellsALL(Real& rcells);

    void volWgtSumALL(const std::string& caller_string,int fast_mode);

    void zap_resid_where_singular(int index_MF);

    void dotSumONES_size(int project_option,
      Vector<Real>& result_sum,
      Vector<int>& result_flag);

    void dotSumONES(int project_option,
      int index_MF,
      Vector<Real>& result_sum);

    void dotSum(int project_option,
      MultiFab* mf1, MultiFab* mf2, Real& result,int nsolve);
    void dot_productALL(int project_option,int index1_MF,int index2_MF,
     Real& result,int nsolve);

    void dot_productALL_ones(int project_option,int index_MF,
		    Vector<Real>& result);
    void dot_productALL_ones_size(int project_option);

    void project_right_hand_side(int index_MF,int project_option,
		    int& change_flag);

    void init_boundary();
    void init_boundary_list(Vector<int> scomp,Vector<int> ncomp);

    void MaxPressureVelocity(Real& minpres,Real& maxpres,
     Real& maxvel,Real& maxvel_collide);
    void MaxPressureVelocityALL(Real& minpres,Real& maxpres,
     Real& maxvel,Real& maxvel_collide);

     // data_dir=-1 cell centered data
     // data_dir=0..sdim-1 face centered data.
     // data_dir=3  X,Y node
     // data_dir=4  X,Z node
     // data_dir=5  Y,Z node
     // a tecplot file is output in which the raw data is
     // output as follows (2D cell center):
     //   x x               x x
     //   x x               x x
     //
     //
     //
     //
     //   x x               x x
     //   x x               x x
     //   Interpolation of data, coarse/fine BCs, fine/fine BCs,
     //   and physical BCs
     //   are inconsequential since the data is assumed to be 
     //   piecewise constant in cells.
     //   For postprocessing, adjust visual_ncell, and call writeTECPLOT_file
     //   (data will be ascii with filenames uniform??????.tec)
    void writeSanityCheckData(
	const std::string& root_string,
	const std::string& information_string,
	const std::string& caller_string,
        int tower_mf_id,
	int ncomp,
	int data_mf,
	int state_type_mf,
	int data_dir,
	int nsteps_actual); 

    virtual void writePlotFile (
     int do_plot,int do_slice,
     int SDC_outer_sweeps_in,
     int slab_step_in,
     int divu_outer_sweeps_in);

    void debug_ParallelCopy();

    void writeInterfaceReconstruction();
    void writeTECPLOT_File(int do_plot,int do_slice);
    void tecplot_debug(FArrayBox& fabdata,
      const Real* xlo,const int* fablo,const int* fabhi,
      const Real* dx,int dir,int id,int comp,int ncomp,
      int interior_only);

    //
    // Timestep estimation functions follow ...
    //
   Real estTimeStep (Real local_fixed_dt,
     const std::string& caller_string);

   virtual void computeInitialDt (int    finest_level,
                                  Real&  dt,
                                  Real   stop_time);

   virtual void computeNewDt (int   finest_level,
                              Real& dt,
                              Real  stop_time);

    virtual Real advance (Real time,Real dt);

    virtual void post_timestep (Real stop_time);
    virtual void post_regrid (int lbase,int start_level,
      int new_finest,int initialInit_flag,Real time);


    //
    // Estimate the end of the simulation for amrLevel.
    //
    virtual int okToContinue ();
    //
    // Build any additional data structures after restart.
    //
    virtual void post_restart ();

    //
    // Insure state, and pressure are consistent.
    //
    virtual void post_init (Real stop_time);
    void prepare_post_process(const std::string& caller_string);

    //
    // Advance grids at this level in time.
    //

    void debug_memory();

    void allocate_mdot();

    void restore_active_pressure(int save_mf);

    void allocate_maccoefALL(int project_option,int nsolve,
	int create_hierarchy,Real dt_maccoef); //-1,0,1
    void allocate_maccoef(int project_option,int nsolve,
        int create_hierarchy,Real dt_maccoef); //-1,0,1
     
    void deallocate_maccoefALL(int project_option);
    void deallocate_maccoef(int project_option);

    void prepare_mask_nbr(int ngrow);
    void resize_mask_nbr(int ngrow);
    void prepare_displacement();

    // VAHAB HEAT SOURCE
    void heat_source_term_flux_source();

    void stefan_solver_init(MultiFab* coeffMF,int adjust_temperature,
		    int project_option);

    void prepare_viscous_solver();

    void alloc_DTDtALL(int alloc_flag);
    void alloc_DTDt(int alloc_flag);

    void prepare_advect_vars(Real time);
    void delete_advect_vars();

    void exit_viscous_solver();

    void diffusion_heatingALL(int source_idx,int idx_heat);
    void diffusion_heating(int source_idx,int idx_heat);

    void explicit_fluxes(int state_idx,int scomp,int iupdate,MultiFab* Un);

     // iten>=1  im1,im2>=1
    static void get_iten_cpp(int im1,int im2,int& iten);

    static int CTML_FSI_numsolids;

    static Real inflow_pressure;
    static Real outflow_pressure;
    static Real period_time;
    static int  MOF_TURN_OFF_LS;
    static int  MOF_DEBUG_RECON;

    static Real visc_coef; 

    static Real angular_velocity;  //Omega in Lewis and Nagata (Lappa's book)

    static Real viscconst_max;
    static Real viscconst_min;

     // 1..7 "1" corresponds to the GNBC option.
    static int ZEYU_DCA_SELECT;  // DCA=dynamic contact angle
    static int num_local_aux_grids; //default is 0
    static Real R_Palmore_Desjardins;
    static int prescribe_temperature_outflow;
     // default is solidheat_flag=0
    static int solidheat_flag;  // 0=diffuse in solid 1=dirichlet 2=neumann
    static Real gravity;
    static int invert_gravity;
    static int gravity_dir;
    static Vector<Real> gravity_vector;
    static int incremental_gravity_flag;
    static int segregated_gravity_flag;

    static Real nucleation_init_time; //def 0.0
     // 0 - centroid furthest from uncaptured centroid
     // 1 - use MOF error
    static int mof_error_ordering;

    static Vector<int> store_elastic_data; // viscoelastic parameter

     // this is set according to store_elastic_data
    static int num_materials_viscoelastic;

     // FSI_FLUID=0
     // FSI_PRESCRIBED_PROBF90=1
     // FSI_PRESCRIBED_NODES=2
     // FSI_ICE_PROBF90=3
     // FSI_SHOELE_CTML=4
     // FSI_RIGID_NOTPRESCRIBED=5
     // FSI_ICE_NODES=6
     // FSI_FLUID_NODES_INIT=7
     // FSI_ICE_STATIC=9
     
    static Vector<int> FSI_flag; 

    static Vector<int> CTML_max_num_nodes_list; 
    static int CTML_max_num_elements_list; 
    static int CTML_FSI_num_scalars;

    static Vector<int> material_type; // 0 incomp   
    static Vector<int> material_type_interface; 

    static int num_species_var; 
    static int num_materials;
    static int num_interfaces;

    static int num_state_material;
    static int num_state_base;
    static int ngeom_raw;
    static int ngeom_recon;

    static Vector<Real> grid_stretching_parameter;

protected:

    void SDC_setup_step();
    void SDC_setup();

    void debug_ngrow(int idxMF,int ngrow,
      const std::string& caller_string);
    void debug_ixType(int idxMF,int grid_type,
      const std::string& caller_string);
    void debug_ixType_raw(MultiFab* mf,int grid_type,
      const std::string& caller_string);
    void debug_boxArray(MultiFab* mf,int grid_type,
      const std::string& caller_string);

    void avgDown_tag_localMF(int idxMF);
    void level_avgDown_tag(MultiFab& S_crse,MultiFab& S_fine);

    void avgDownDRAG_MF();
    void level_avgDownDRAG(MultiFab& S_crse,MultiFab& S_fine);

    void avgDownBURNING_localMF(
		int burnvel_MF,
		int TSAT_MF);
    void level_avgDownBURNING(MultiFab& S_crse,MultiFab& S_fine,int velflag);

    void avgDownCURV_localMF(int idxMF);
    void level_avgDownCURV(MultiFab& S_crse,MultiFab& S_fine);

    void avgDown_and_Copy_localMF(
      int idx_den_MF,
      int idx_vel_MF,
      int idx_flux_MF,
      int operation_flag);

    void interp_and_Copy_localMF(
      int idx_den_MF,
      int idx_vel_MF,
      int idx_flux_MF,
      int operation_flag);

    void interp_flux_localMF(
      int coarse_flux_MF,
      int fine_flux_MF);

    void sync_flux_var(
      int dir,
      int flux_MF,
      int ncomp_flux);

     // spectral_override==1 => order derived from "enable_spectral"
     // spectral_override==0 => always low order.
     // synchronizes level+1 and level.
    void avgDown_localMF(int idxMF,int scomp,int ncomp,
     int spectral_override);
    void Plus_localMF(int idx_dest,Real val,int dcomp,int ncomp,int ngrow);
    void Mult_localMF(int idx_dest,int idx_source,
     int scomp,int dcomp,int ncomp,int ngrow);
    void Copy_localMF(int idx_dest,int idx_source,
     int scomp,int dcomp,int ncomp,int ngrow);
    void minus_localMF(int idx_dest,int idx_source,
     int scomp,int ncomp,int ngrow);

      // grid_type=-1 ... 5
    void new_localMF(int idx_MF,int ncomp,int ngrow,int grid_type);
      // grid_type=-1 ... 5
    void new_localMF_if_not_exist(int idx_MF,int ncomp,int ngrow,
		    int grid_type);

    //dir=0..sdim-1
    void getStateMAC_localMF(int idx_MF,int ngrow,int dir,Real time);

    void getStateDen_localMF(int idx_MF,int ngrow,Real time);

    void getStateDist_localMF(int idx_MF,int ngrow,Real time,
       const std::string& caller_string);

    void resize_levelset(int ngrow,int idx);
    void allocate_levelset(int ngrow,int idx);
    void allocate_levelset_ALL(int ngrow,int idx);

    void getState_localMF(int idx_MF,int ngrow,
     int scomp,int ncomp,Real time);

    void getStateTensor_localMF(int idx_MF,int ngrow,
     int scomp,int ncomp,Real time);

    void getState_localMF_list(
     int idx_MF,int ngrow,
     int state_index,
     Vector<int> scomp,
     Vector<int> ncomp);

    void getState_localMF_listALL(
     int idx_MF,int ngrow,
     int state_index,
     Vector<int> scomp,
     Vector<int> ncomp);

    void putState_localMF_list(
     int idx_MF,
     int state_index,
     Vector<int> scomp,
     Vector<int> ncomp);

    void putState_localMF_listALL(
     int idx_MF,
     int state_index,
     Vector<int> scomp,
     Vector<int> ncomp);

    void getBCArray_list(Vector<int>& listbc,int state_index,
     int gridno,Vector<int> scomp,Vector<int> ncomp);

    void resize_maskfiner(int ngrow,int mask_id);
    void maskfiner_localMF(int idx_MF,int ngrow,
     Real tag,int clearbdry);

    void avgDownError_ALL();
    void avgDownError();

    void getStateVISC(const std::string& caller_string);
    void getStateVISC_ALL(const std::string& caller_string);

    void getStateCONDUCTIVITY();
    void getStateCONDUCTIVITY_ALL();

    void delete_localMF(int idx_MF,int ncomp); 
    void delete_localMF_if_exist(int idx_MF,int ncomp); 

    void MaxAdvectSpeedALL(
      Real& dt_min,
      Vector<Real>& vel_max_estdt,
      const std::string& caller_string);

    void MaxAdvectSpeed(
      Real& dt_min,
      Vector<Real>& vel_max_estdt,
      const std::string& caller_string);

    void DumpProcNum();

      // interpolates coarsedata from level-1 to level.
    void AllinterpScalarMAC(
      MultiFab* coarsedata,MultiFab* finedata,
      BoxArray& cgridscen,BoxArray& fgridscen,
      DistributionMapping& fdmap,
      MultiFab* cdiagsing,
      int nsolve,int project_option);
    void interpScalarMAC(MultiFab* coarsedata,MultiFab* finedata,
      int nsolve,int project_option);
    void Allaverage(
      MultiFab* coarsedata,MultiFab* finedata,
      BoxArray& cgridscen,BoxArray& fgridscen,
      DistributionMapping& fdmap,
      int iaverage,
      int scomp,int dcomp);
    void averageRhs(int idx_MF,int nsolve,int project_option);
    int NSnumLevels();

    int ns_thread();
    void ns_reconcile_d_num(int caller_loop_id,
       const std::string& caller_string);

    static int some_materials_compressible();

    void get_mm_scomp_solver(
     int num_materials_combine,
     int project_option,
     int& state_index,
     Vector<int>& scomp,
     Vector<int>& ncomp,
     int& ncomp_check);

// pressure/vel/temp/spec=0
    void zero_independent_variable(int project_option,int nsolve); 
    void zero_independent_vel(int project_option,int idx,int nsolve); 

    static void setup_integrated_quantities ();

    static void read_params ();   // Read input file
    static void read_geometry (); // Set metrics

    static void override_enable_spectral(int enable_spectral_in);
    static void override_enable_spectralGHOST(
        int scomp,int ncomp,int enable_spectral_in);

    static void set_tensor_extrap_components(
      int coord,std::string postfix,int indx,int ibase_tensor);
    static void set_tensor_extrap_components_main(
      int coord,std::string postfix,int indx);

    void sum_integrated_quantities (
     const std::string& caller_string,
     Real stop_time);

    void CopyNewToOldALL();
    void CopyOldToNewALL();

    void Geometry_cleanup ();
    void Geometry_setup ();

    void copybc(Vector<int> dest,Vector<int> source,
                int scomp,int dcomp,int ncomp);

    void combine_state_variable(
     int project_option,
     int idx,
     int combine_flag,
     int hflag,
     int update_flux,
     int interface_cond_avail);

    void viscous_boundary_fluxes(
     int project_option,
     MultiFab* xflux,MultiFab* yflux,MultiFab* zflux,
     int nsolve);

    void tensor_extrapolation();
    void tensor_advection_update();
    void tensor_advection_updateALL();

    void level_init_icemask_and_icefacecut();
    void phase_change_redistributeALL();
    void level_phase_change_redistribute(
      Real expect_mdot_sign,
      int im_source,int im_dest,int indexEXP,
      int isweep);

    void clear_blobdata(int i,Vector<blobclass>& blobdata);
    void copy_from_blobdata(int i,int& counter,
     Vector<Real>& blob_array,Vector<blobclass>& blobdata);
    void sum_blobdata(int i,
     Vector<blobclass>& blobdata,
     Vector<blobclass>& level_blobdata,int sweep_num);
    void copy_to_blobdata(int i,int& counter,
     Vector<Real>& blob_array,Vector<blobclass>& blobdata);
    void copy_blobdata(Vector<blobclass>& dest_blobdata,
       Vector<blobclass>& source_blobdata);

    void build_NRM_FD_MF(int fd_mf,int ls_mf);

    void level_phase_change_convert(
      int im_outer,int im_opp_outer,
      int i_phase_change,int n_phase_change);

    void level_phase_change_convertALL();

    void level_species_reaction(const std::string& caller_string);

    void level_phase_change_rate(Vector<blobclass> blobdata,
       int color_count,int nucleation_flag);
    void level_phase_change_rate_extend();

    void level_DRAG_extend();
   
    void allocate_SDC();
    void deallocate_SDC();

    void init_delta_SDC();
    void SEM_advectALL(int source_term);
    void prelim_alloc();
    void update_SEM_forcesALL(int project_option,
     int idx_source,int update_spectral,int update_stable);
    void update_SEM_forces(int project_option,
     int idx_source,int update_spectral,int update_stable);
    void init_splitting_force_SDC();

    void init_gradu_tensor(
     int homflag,
     int idx_vel,
     int idx_cell,
     int idx_face,
     int simple_AMR_BC_flag_viscosity);

    void doit_gradu_tensor(
     int homflag,
     int idx_vel,
     int idx_cell,
     int idx_face,
     int spectral_loop,
     int itensor_iter,
     MultiFab* mask3,
     int simple_AMR_BC_flag_viscosity);

    void init_gradu_tensor_and_material_visc_ALL(
      const std::string& caller_string);

    void init_gradu_tensorALL(
     int idx,
     int do_alloc,
     int idx_cell,
     int idx_face,
     int simple_AMR_BC_flag_viscosity);

    void make_viscoelastic_tensor(int im);
    void make_viscoelastic_tensorALL(int im);

    void make_viscoelastic_tensorMAC(int im,
	int flux_mf,int flux_grid_type,int fill_state_idx);
    void make_viscoelastic_tensorMACALL(int im,
        int flux_mf,int flux_grid_type,
        int fill_state_idx);

    void make_viscoelastic_heating(int im,int idx);

    // VAHAB HEAT SOURCE
    void make_heat_source();

    void add_perturbation();
    void make_SEM_delta_force(int project_option);

    void getStateMOM_DEN(int idx,int ngrow,Real time);

    void check_grid_places();

    void show_norm2(MultiFab* mf,int scomp,int ncomp,int id);
    void show_norm2_id(int mf_id,int id);
    void check_value_max(int id,int mf_id,int scomp, 
		int ncomp,int ngrow,Real max_value);
    void check_value_max_level(int id,MultiFab* mf,
		int scomp,int ncomp,int ngrow,Real max_value);

    void aggressive_debug(
      int datatype,
      int force_check,
      MultiFab* mf,
      int scomp,int ncomp,
      int ngrow,
      int dir,
      Real warning_cutoff);

    void FillBoundaryTENSOR(MultiFab* mf,int sc,int dir);

    Box growntileboxTENSOR(
     MultiFab* mf,int datatype,int ng,int dir,
     const Box& tilebox_box,const Box& vbx_box);
    bool contains_nanTENSOR(MultiFab* mf,
      int datatype,int scomp,int dir);

    void check_for_NAN_TENSOR_base(int datatype,MultiFab* mf,
      int sc,int dir);
    void check_for_NAN_TENSOR(int datatype,MultiFab* mf);
    void check_for_NAN(MultiFab* mf);

    void split_scalar_advection(); 
    void split_scalar_advectionALL();
 
    void SEM_scalar_advection (int init_fluxes,int source_term,
     int spectral_loop,int tileloop);
    void synchronize_flux_register(int operation_flag,int spectral_loop);

    int end_spectral_loop();
    void allocate_flux_register(int operation_flag);
 
    void PCINTERP_fill_borders(int idx_MF,int ngrow,int scomp,int ncomp,
      int index,Vector<int> scompBC_map);
    void PCINTERP_fill_bordersALL(int idx_MF,int ngrow,int scomp,int ncomp,
      int index,Vector<int> scompBC_map);

    void PCINTERP_fill_coarse_patch(int idx_MF,
      int scomp,int ncomp,
      int index,Vector<int> scompBC_map);

    void GetStateFromLocal(int idx_MF,int ngrow,int scomp,int ncomp,
      int index,Vector<int> scompBC_map);
    void GetStateFromLocalALL(int idx_MF,int ngrow,int scomp,int ncomp,
      int index,Vector<int> scompBC_map);

  void mf_combine(
    int project_option,
    int index_x_MF,
    int index_y_MF,
    Real Beta,
    int index_z_MF,int nsolve);

  void mf_combine_ones(int project_option,int index_MF,Vector<Real> Beta);
  void mf_combine_ones_level(int project_option,int index_MF,
		  Vector<Real> Beta);

  void levelCombine(
   int project_option,
   MultiFab* mfx, MultiFab* mfy, MultiFab* mfz,
   Real beta,int ncomp);

    void build_masksemALL();
    void build_masksem(int mask_sweep);

    void MOF_training();

     // update_flag=RECON_UPDATE_NULL
     // update_flag=RECON_UPDATE_STATE_ERR
     // update_flag=RECON_UPDATE_STATE_CENTROID
     // update_flag=RECON_UPDATE_STATE_ERR_AND_CENTROID
    void VOF_Recon(int ngrow,
      Real time,int update_flag,
      int init_vof_prev_time,
      Real& delta_centroid_level,
      int& number_centroid_level); //output: SLOPE_RECON_MF
    void VOF_Recon_ALL(
      const std::string& caller_string,
      int ngrow,Real time,int update_flag,
      int init_vof_prev_time); //output: SLOPE_RECON_MF
    void VOF_Recon_resize(int ngrow); //output: SLOPE_RECON_MF

    void make_physics_varsALL(int project_option,
	   const std::string& caller_string);

    void allocate_physics_vars();
    void make_physics_vars(int project_option,
      const std::string& caller_string);

    void solid_temperature();

    void nonlinear_advection(const std::string& caller_string);

    void advance_MAC_velocity(int project_option);

    void pressure_gradient_code_segment(
      const std::string& caller_string);

    void no_mass_transfer_code_segment(
      const std::string& caller_string);

    void phase_change_code_segment(
      const std::string& caller_string,
      int& color_count,
      Vector<blobclass>& blobdata);

    void nucleation_code_segment(
      const std::string& caller_string,
      int& color_count,
      Vector<blobclass>& blobdata,
      Vector< Vector<Real> >& mdot_data,
      Vector< Vector<Real> >& mdot_comp_data,
      Vector< Vector<Real> >& mdot_data_redistribute,
      Vector< Vector<Real> >& mdot_comp_data_redistribute,
      Vector<int>& type_flag);

    void do_the_advance(Real timeSEM,Real dtSEM,int& advance_status);

    void correct_colors(
     int idx_color,int base_level,
     Vector<int> domaincolormap,
     int max_colors_level);

    void assign_colors(
     int& fully_covered,
     int idx_color,int idx_type,
     Vector<int>& colormax,Vector<int> type_flag,
     int zero_diag_flag);

    void avgDownColor(int idx_color,int idx_type);
    MultiFab* CopyFineToCoarseColor(int idx_color,int idx_type,
		    int zero_diag_flag);

    void sync_colors(
     int idx_color,int idx_type,
     Vector<int> color_per_grid,
     Vector<int>& colormax,
     int max_colors_grid,
     MultiFab* maskmf,
     int check_corners,
     int zero_diag_flag);

    void color_variable(
     int& coarsest_level,
     int idx_color,int idx_type,
     int* color_count,
     Vector<int> type_flag,
     int zero_diag_flag);

    void SumRegions(int isweep); //isweep = 0 or 1

    void ColorSum(
     int operation_flag, // OP_GATHER_MDOT or OP_SCATTER_MDOT
     int tessellate, // =1 or 3
     int sweep_num,
     int ncomp_mdot_alloc,
     int ncomp_mdot,
     MultiFab* typemf,
     MultiFab* color,
     MultiFab* mdot,
     MultiFab* mdot_complement,
     Vector<blobclass>& level_blobdata,
     Vector<blobclass> cum_blobdata,
     Vector< Vector<Real> >& level_mdot_data,
     Vector< Vector<Real> >& level_mdot_comp_data,
     Vector< Vector<Real> > cum_mdot_data,
     Vector< Vector<Real> > cum_mdot_comp_data,
     Vector< Vector<Real> >& level_mdot_data_redistribute,
     Vector< Vector<Real> >& level_mdot_comp_data_redistribute 
     );

    void ColorSumALL(
     int operation_flag, // OP_GATHER_MDOT or OP_SCATTER_MDOT
     int tessellate,  // =1 or 3
     int coarsest_level,
     int& color_count,
     int idx_type,
     int idx_color,
     int idx_mdot,  // ==-1 if no mdot
     int idx_mdot_complement,  // ==-1 if no mdot_complement
     Vector<int>& type_flag, 
     Vector<blobclass>& blobdata,
     Vector< Vector<Real> >& mdot_data,
     Vector< Vector<Real> >& mdot_comp_data,
     Vector< Vector<Real> >& mdot_data_redistribute,
     Vector< Vector<Real> >& mdot_comp_data_redistribute
     );



    void LowMachDIVU(
     int sweep_num,
     MultiFab* typemf,
     MultiFab* color,
     MultiFab* mdot_local,
     MultiFab* mdot_global,
     Vector<blobclass> cum_blobdata,
     Vector< Vector<Real> >& level_mdot_data,
     Vector< Vector<Real> >& cum_mdot_data
     );

    void LowMachDIVUALL(
     int coarsest_level,
     int& color_count,
     int idx_type,
     int idx_color,
     Vector<int>& type_flag, 
     Vector<blobclass> blobdata);


    void Type_level(
     MultiFab* typemf,Vector<int>& type_flag,
     int zero_diag_flag);

    void TypeALL(
     int idx_type,
     Vector<int>& type_flag,
     int zero_diag_flag);

    virtual void post_init_state ();


    void avgDown(MultiFab& S_crse,MultiFab& S_fine,int scomp,int ncomp,
       int spectral_override);
    void MOFavgDown();
    void avgDownMac(); 
    void avgDownMacState(int spectral_override); 

    void avgDownEdge(int grid_type,MultiFab& S_crse,MultiFab& S_fine,
     int scomp,int ncomp,int spectral_override,
     const std::string& caller_string);

     // spectral_override==1 => order derived from "enable_spectral"
     // spectral_override==0 => always low order.
    void avgDownEdge_localMF(int idxMF,int scomp,int ncomp,
     int start_grid_type,int n_grid_type,
     int spectral_override,
     const std::string& caller_string);
     
    void smearvariable(int scomp,
     MultiFab* mask_extrapolate,int ngrow,Real def_value);

     // spectral_override==1 => order derived from "enable_spectral"
     // spectral_override==0 => always low order.
    void avgDown_localMF_ALL(int idxMF,int scomp,int ncomp,
      int spectral_override); 

    void avgDownALL(int stateidx,int startcomp,int numcomp,
      int spectral_override); 

    void avgDownALL_TENSOR(); 

     // spectral_override==1 => order derived from "enable_spectral"
     // spectral_override==0 => always low order.
    void avgDown(int stateidx,int startcomp,int numcomp,int spectral_override); 
    void avgDown_list(int stateidx,Vector<int> scomp,
      Vector<int> ncomp,int spectral_override); 

    void prescribe_solid_geometryALL(Real time,
      int renormalize_only,
      int local_truncate,
      const std::string& caller_string);

    void prescribe_solid_geometry(Real time,int renormalize_only);
#ifdef AMREX_PARTICLES
    using My_ParticleContainer =
      AmrParticleContainer<N_EXTRA_REAL,N_EXTRA_INT,0,0>;
    void init_particle_containerALL(int append_flag,
      const std::string& caller_string);
    void init_particle_container(int append_flag,
      const std::string& caller_string);
    void move_particles(My_ParticleContainer& localPC,
       const std::string& caller_string);
#endif

    void truncate_VOF(Vector<Real>& delta_mass_all);

    void getState_tracemag_ALL(int idx);
    void getState_tracemag(int idx);

    static void print_project_option(int project_option);

    void output_triangles();

      //declared in: NavierStokes2.cpp
    void output_zones(
      int plot_grid_type,
      FArrayBox& visual_fab_output,
      Box& visual_domain,
      int visual_ncomp,
      MultiFab* velmf,
      MultiFab* presmf,
      MultiFab* divmf,
      MultiFab* div_data,
      MultiFab* denmf,
      MultiFab* mom_denmf,
      MultiFab* viscoelasticmf,
      MultiFab* lsdistmf,
      MultiFab* viscmf,
      MultiFab* conductmf,
      MultiFab* magtracemf,
      MultiFab* elasticforcemf,
      MultiFab* gradvelocitymf,
      int& grids_per_level,
      BoxArray& cgrids_minusBA,
      Real* slice_data,
      int do_plot,int do_slice);

     // data_dir=-1 cell centered data
     // data_dir=0..sdim-1 face centered data.
     // data_dir=3  X,Y node
     // data_dir=4  X,Z node
     // data_dir=5  Y,Z node
    void Sanity_output_zones(
      const std::string& information_string,
      int tower_mf_id,
      int data_dir,
      MultiFab* datamf,
      int ncomp,
      int& grids_per_level,
      BoxArray& cgrids_minusBA);

    void adjust_tolerance(Real error0,Real& error0_max,
      int project_option);

void JacobiALL(
 int idx_resid,int idx_rhs,int idx_xnew,
 int project_option,int nsolve);

void DiagInverse(MultiFab* resid,MultiFab* xnew,int nsolve,int project_option);

void residALL(
  int project_option,
  int idx_rhs,int idx_resid,int idx_phi,int nsolve);

void applyALL(
  int project_option,
  int idx_phi,int idx_Aphi,int nsolve);

void applyGradALL(
  int project_option,int idx_phi,int nsolve);

void applyBC_MGLEVEL(int idx_phi,
 MultiFab* pbdry,int homflag,int nsolve,int project_option);

void applyBC_LEVEL(int project_option,int idx_phi,int nsolve);

void apply_div(
  int project_option,int homflag,
  int idx_phi,
  MultiFab* rhsmf,
  MultiFab* mdotcell,
  int idx_gphi,
  int nsolve,
  Real dt_apply_div);

void mac_project_rhs(int project_option,
 int idx_mac_phi_crse,int idx_mac_rhs_crse,int nsolve);

void mac_update(MultiFab* mac_phi_crse,int project_option,int nsolve);

void diffuse_coefficients();

void user_defined_momentum_forceALL(int idx_vel,int idx_thermal);
void user_defined_momentum_force(int idx_vel,int idx_thermal);

void diffuse_hoopALL(int idx_vel,int idx_thermal,
 int idx_force,int update_state);
void diffuse_hoop(int idx_vel,int idx_thermal,
 int idx_force,int update_state);

void apply_pressure_grad(
 int simple_AMR_BC_flag,
 int simple_AMR_BC_flag_viscosity,
 int homflag,
 int energyflag,
 int gp_mf,
 int pboth_mf,
 int project_option,int nsolve,
 Real dt_pressure_grad);

    void relaxLEVEL(
     MultiFab* rhsmf,
     int idx_rhs,int idx_phi,
     int presmooth,
     int project_option,int nsolve);

    void mg_cycleALL(int presmooth,
      int project_option,
      int idx_rhs,int idx_phi,
      int nsolve);

    void check_outer_solver_convergence(
	Real error_n,Real error0,
        Real save_mac_abs_tol_in,
        int& meets_tol);

    void jacobi_cycles(
     int call_adjust_tolerance,
     int ncycles,
     int update_vel,int project_option,
     int idx_mac_rhs_crse,
     int idx_mac_phi_crse,
     Real& error_at_the_beginning,
     Real& error_after_all_jacobi_sweeps,
     Real& error0,Real& error0_max,
     int krylov_subspace_num_outer_iterSOLVER,int nsolve);

    void updatevelALL(
     int project_option,
     int idx_mac_phi_crse,int nsolve);

    void extend_pressure_to_solidALL();

    void Krylov_checkpoint(
      int vcycle,
      Real krylov_error,
      Real& best_error,
      int& best_iter,
      int idx_phi,
      int idx_umac,
      int& restart_flag);

    void multiphase_project(int project_option);

    static void set_local_tolerances(int project_option);

    void multiphase_preconditioner(
     int project_option,int project_timings,
     int presmooth,int postsmooth,
     int idx_Z,int idx_R,int nsolve);

    void multiphase_SHELL_preconditioner(
     int project_option,int project_timings,
     int presmooth,int postsmooth,
     int idx_Z,int idx_R,int nsolve);

    void Prepare_UMAC_for_solver(int project_option,int nsolve);
    void remove_UMAC_for_solver(int project_option);
    void allocate_MAC_velocityALL(int nsolve,int idx);
    void remove_MAC_velocityALL(int idx);
    void remove_project_variables();
    void remove_pressure_work_vars();
    void allocate_pressure_work_vars(int nsolve,int project_option);

    void allocate_project_variables(int nsolve,int project_option);

    void remove_FACE_WEIGHT_vars();

    void sanity_check_face_wt(int project_option);

    void allocate_FACE_WEIGHT(
      int nsolve,
      int project_option,
      int face_weight_op); //SUB_OP_FOR_MAIN or SUB_OP_FOR_SDC

    void residual_correction_form(
      int homflag,int energyflag,
      int project_option,int nsolve);

    void correct_velocity(
     int project_option,
     int macdest,
     int macsrc,
     int gp,int nsolve);

    void overwrite_outflow();

    void scale_variablesALL();
    void unscale_variablesALL();
    void scale_variables(int scale_flag); // scale_flag=0=>scale, 1=>unscale

    void vel_elastic_ALL(int viscoelastic_force_only);

    void veldiffuseALL();

    void Mass_Energy_Sources_SinksALL();

    static int is_zalesak();

    void zalesakVEL();

    static void grid_type_to_box_type_cpp(int grid_type,int* box_type);

     // iten>=1  im1,im2>=1
    static void get_inverse_iten_cpp(int& im1,int& im2,int iten);

     // unew^f=unew^f + beta * localMF[idx_velcell]
    void increment_face_velocity(
      int operation_flag,
      int project_option,
      int idx_velcell,
      Real beta,
      Vector<blobclass> blobdata);

    void increment_face_velocityALL(
      int operation_flag,
      int project_option,
      int idx_velcell,
      Real beta,
      Vector<blobclass> blobdata);

    void project_to_rigid_velocity(Vector<blobclass> blobdata);

    void project_to_rigid_velocityALL();

    void make_MAC_velocity_consistentALL();
    void make_MAC_velocity_consistent();

    void VELMAC_TO_CELLALL(int dest_idx);

    void VELMAC_TO_CELL(int dest_idx);

    void INCREMENT_REGISTERS_ALL(int source_mf);
    void INCREMENT_REGISTERS(int source_mf);
    void push_back_state_register(int idx_MF,Real time);
    void SET_STOKES_MARK(int idx_MF);

    void APPLY_VISCOUS_HEATING(int source_mf);
    void APPLY_VISCOUS_HEATINGALL(int source_mf);

    MultiFab* derive_EOS_pressure(Vector<int> local_material_type);
    void init_pressure_error_indicator();
    void init_advective_pressure(int project_option);

    void level_getshear(
      MultiFab* shear_output_mf,
      MultiFab* vel_mf,
      int only_scalar,
      int destcomp,
      int ngrow);

    void save_to_macvel_state(int idx_umac);

    void init_divup_cell_vel_cell(
     int project_option,
     int energyflag,
     int idx_pres,  
     int idx_umac);

    void update_SEM_delta_force(
     int project_option,
      //GP_DEST_FACE_MF = grad p (MAC grid) SOLVETYPE_PRES
      //MACDIV_MF=-div(2 mu D) SOLVETYPE_VISC
      //MACDIV_MF=-div(k grad T) SOLVETYPE_HEAT
      //HOOP_FORCE_MARK_MF included in update_SEM_delta_force.
     int update_spectral,
     int update_stable,
     int nsolve);

    void increment_potential_forceALL();
    void increment_potential_force();

    MultiFab* getState_list(int ngrow,Vector<int> strt_comp,
       Vector<int> num_comp,Real time);

    void putState_list(Vector<int> strt_comp,
       Vector<int> num_comp,int idx_MF);

    MultiFab* getState(int ngrow,int strt_comp,
       int num_comp,Real time);

    MultiFab* getStateSolid(int ngrow,int strt_comp,
       int num_comp,Real time);

    MultiFab* getStateTensor(int ngrow,int strt_comp,
       int num_comp,Real time);

    void CELL_GRID_ELASTIC_FORCE(int im_elastic);

    MultiFab* getStateMAC(int ngrow,int dir,Real time);

    void getStateDIV_ALL(int idx_source,int scomp_src,
      int idx_dest,int idx_mask);
    void getStateDIV(int idx_source,int scomp_src,
      int idx_dest,int idx_mask);

    void ADVECT_DIV_ALL();
    void ADVECT_DIV();

    void makeCellFrac(int tessellate,int ngrow,int idx);
    void makeFaceFrac(int tessellate,int ngrow,int idx);
    void makeFaceTest(int tessellate,int ngrow,int idx);
    void ProcessFaceFrac(int tessellate,int idxsrc,int idxdst,
		    int ngrow_dest);
    void makeStateDistALL(int keep_all_interfaces,
       int ngrow_make_distance_accept);
    void makeStateDist(int keep_all_interfaces,
       int ngrow_make_distance_accept);
    void correct_dist_uninit();

    void makeStateCurv(int project_option,
	  const std::string& caller_string);

    void make_marangoni_force();

    void delete_array(int idx_localMF);

    void allocate_array(int ngrow,int ncomp,int grid_type,
      int idx_localMF);

    void allocate_independent_var(int nsolve,int idx);
    void allocate_rhs_var(int nsolve,int idx);

    int get_new_data_Type(int mfab_id);

    void Copy_array(int mfab_dest,int mfab_source,
		    int scomp,int dcomp,int ncomp,int ngrow);

    void setVal_array(int ngrow,int scomp,int ncomp,Real dataval,
      int idx_localMF);

    void setVal_localMF(int idx,Real dataval,
      int scomp,int ncomp,int ngrow);

    void mult_array(int ngrow,int ncomp,Real dataval,
      int idx_localMF);
    void minusALL(int ngrow,int scomp,int ncomp,int idx_dest,int idx_source);
    void plusALL(int ngrow,int scomp,int ncomp,int idx_dest,int idx_source);

    void negateALL(int idx_localMF);
    void zeroALL(int ngrow,int ncomp,int idx_localMF);

    MultiFab* getStateDen(int ngrow,Real time);

    MultiFab* getStateDist(int ngrow,Real time,
       const std::string& caller_string);

    MultiFab* getStateDIV_DATA(int ngrow,int scomp,int ncomp,Real time);

    void putStateDIV_DATA(int scomp,int ncomp,int idx_MF);

    void getStateALL(int ngrow,Real time,int scomp,
      int ncomp,int idx_localMF);

    MultiFab* getStatePres(int ngrow,Real time);

    void resize_metrics(int ngrow);
    void metrics_data(int ngrow);

    void metrics_data_min_max(const std::string& caller_string);

    void metrics_dataALL(int ngrow);

    void metrics_data_min_max_ALL(const std::string& caller_string);

    MultiFab* maskfiner(int ngrow,Real tag,int clear_phys_boundary);
    MultiFab* maskfiner_noghosts(Real tag);

    NavierStokes& getLevel (int lev)
    {
        return *(NavierStokes*) &parent->getLevel(lev);
    }

    void init_regrid_history();

    void init_aux_data();

     //cur_time=prev_time+dt
    void ns_header_msg_level(
      int FSI_operation,
      int FSI_sub_operation,
      Real cur_time,
      Real dt,
      int iter,
      const std::string& caller_string);

    static int read_from_CAD();
    static int is_ice_matC(int im);
    static int is_singular_coeff(int im);
    static int is_FSI_rigid_matC(int im);
    static int CTML_FSI_flagC();
    static int CTML_FSI_matC(int im);
    static int FSI_material_exists_CTML();
    static int FSI_flag_valid(int im); // 0<=im<num_materials

    void init_FSI_GHOST_MAC_MF_ALL_predict();
    void init_FSI_GHOST_MAC_MF_predict();

    void init_FSI_GHOST_MAC_MF_ALL(
      int renormalize_only,
      const std::string& caller_string);

    void init_FSI_GHOST_MAC_MF(int dealloc_history);

    void assimilate_state_data();

    void resize_FSI_MF();

      //cur_time=prev_time+dt
    void FSI_make_distance(Real cur_time,Real dt);
    void regenerate_from_eulerian(Real cur_time);

    void create_fortran_grid_struct(Real cur_time,Real dt);

    void cpp_overridepbc(int homflag_in,int project_option_in);

    static int ns_is_rigid(int im);
    static int ns_is_lag_part(int im);

    void copy_velocity_on_sign(int partid);
    void Transfer_FSI_To_STATE(Real cur_time);

    void copy_old_FSI_to_new_level();
    void copy_old_FSI_to_new();
    int ok_copy_FSI_old_to_new();

    void build_moment_from_FSILS(Real cur_time);

protected:

    static int POLYGON_LIST_MAX;

    static bool ns_tiling; // default is false

    static int profile_debug;

    static int ns_time_order;
     // before: get_old_data  now: get_new_data(slab_step)
     // before: get_new_data  now: get_new_data(slab_step+1)
    static int slab_step; 
    static int SDC_outer_sweeps;
    static int divu_outer_sweeps;
    static int very_last_sweep;
    static Real upper_slab_time;
    static Real lower_slab_time;
    static Real delta_slab_time;
    static Real prev_time_slab;
    static Real cur_time_slab;
    static Real vel_time_slab;
    static Real advect_time_slab;
    static Real prescribed_vel_time_slab;
    static Real dt_slab;
    static int advect_iter;
    static int hydrate_flag;

    static int idx_umac_material_mf;
    static int idx_umac_mask_material_mf;
    static int idx_scalar_mask_material_mf;

    static Vector<int> type_ONES_flag;
    static int color_ONES_count;
    static int coarsest_ONES_level;

    static Vector<Real> outer_error_all_solver_calls;
    static Vector<int> number_vcycles_all_solver_calls;
    static Vector<int> max_lev0_cycles_all_solver_calls;
    static Vector<int> median_lev0_cycles_all_solver_calls;
    static Vector<int> lev0_cycles_list;
    static Vector<int> number_solver_calls;

    static Vector<int> map_forward_direct_split;
    static Vector<int> normdir_direct_split;
    static int dir_absolute_direct_split;
    static int order_direct_split; // base_step mod 2

    static Vector< Vector<Real> > min_face_wt;
    static Vector< Vector<Real> > max_face_wt;

    static Vector< Vector<Real> > delta_mass;

    static Real max_problen;
    static Vector< Vector<Real> > minLS;
    static Vector< Vector<Real> > maxLS;

    static Vector<Real> mdotplus;
    static Vector<Real> mdotminus;
    static Vector<Real> mdotcount;
    static Vector<Real> mdot_sum;
    static Vector<Real> mdot_sum2;
    static Vector<Real> mdot_lost;

    static Vector<Real> mdotplus_complement;
    static Vector<Real> mdotminus_complement;
    static Vector<Real> mdotcount_complement;
    static Vector<Real> mdot_sum_complement;
    static Vector<Real> mdot_sum2_complement;
    static Vector<Real> mdot_lost_complement;

    static Vector<Real> curv_min;
    static Vector<Real> curv_max;

    static int override_bc_to_homogeneous;

    static Real total_advance_time;
    static BCRec phys_bc;
    static BCRec viscosity_phys_bc;
    static BCRec temperature_phys_bc;
    static BCRec species_phys_bc;

    static int NS_geometry_coord;

    ABecLaplacian* mac_op;

      // localMF[mf_id] is a MultiFab defined on the BoxArray of the
      // present level.  The idea behind this data structure is as follows:
      // suppose one wants to interpolate data from level-1 to the present
      // level: then one can pass "mf_id" on the present level, and grab
      // the coarse level as follows:
      // ns_coarse=getLevel(level-1)
      // coarse data is then ns_coarse.localMF[mf_id]
      // The main purpose of localMF structure is to have the same 
      // convenience for non-state-data as state data.
    Array<int,MAX_NUM_LOCAL_MF> localMF_grow;
    Array<MultiFab*,MAX_NUM_LOCAL_MF> localMF;

#ifdef AMREX_PARTICLES

    using My_NBR_ParticleContainer =
      NeighborParticleContainer<N_EXTRA_REAL,N_EXTRA_INT,0,0>;

    My_NBR_ParticleContainer* NBR_Particle_Container;

#endif

     // value_color = sum_{i in omega_color} 1
    static Vector<Real> ones_sum_global;
     // for each given color, singular_patch_flag=
     //   0 if color is masked off 
     //   1 if color is not masked off, no compressible/internal dirichlet 
     //     regions, and not touching a Dirichlet condition wall.
     //   2 if color is not masked off, a compressible/internal dirichlet
     //     region exists or color is touching a Dirichlet condition wall.
    static Vector<int> singular_patch_flag;

    int save_mg_precond;

    Vector<int> bcpres_array;
    Vector<Real> xarray;
    int icnt_xarray;


    int is_first_step_after_regrid;
    BoxArray old_intersect_new;

    //
    // Internal (protected) parameters for timestepping.
    //
    static Real init_shrink; // reduction factor of first esimated timestep
    static Real cfl;         // desired maximum cfl

    static int  MOFITERMAX;
    static int  MOFITERMAX_AFTER_PREDICT;

// default: tessellating fluid => default==1
//          non-tesselating or tesselating solid => default==0
    static Vector<int> truncate_volume_fractions;  

    static int particle_nsubdivide;  
    static int particle_max_per_nsubdivide;  

    static Real truncate_thickness;  // default is 2.0

    static int  enable_spectral;

     // continuous_mof = STANDARD_MOF => regular MOF
     // continuous_mof > STANDARD_MOF => CMOF(n) ("continuous MOF ")
    static int  continuous_mof; // default is CMOF_X
    static int  update_centroid_after_recon; //default is 0
    static int  mof_decision_tree_learning; // default is 10000
    static int  mof_machine_learning; // default is 0
    static Vector<Real> centroid_noise_factor; //default is 0.0

    static int  partial_cmof_stencil_at_walls;  // default is 1
    static Real change_max;    // maximum change in dt over a timestep
    static Real change_max_init; // maximum change in dt over a timestep

    static Real fixed_dt;      // set > 0 to specify dt
    static Real fixed_dt_init; // set > 0 to specify dt
    static Real min_velocity_for_dt; // dt<dx/min_velocity_for_dt
     //uu_stefan=max(uu_stefan,min_stefan_velocity_for_dt) 
    static Real min_stefan_velocity_for_dt; 
    static Real fixed_dt_velocity; // set > 0 to specify u for dt=dx/u
    static Real dt_max;        // maximum dt allowed
    static int  sum_interval;  // number of timesteps for conservation stats
     // refine a lagrangian element if hside>refine_factor*hfinest
     // do not modify element structure if refine_factor=0
    static Vector<int> FSI_refine_factor;  
     // bounding box ngrow for finding box of cells to visit for a 
     // given lagrangian element.
    static Vector<int> FSI_bounding_box_ngrow;

    static Vector<int> ns_max_grid_size;

    //
    // Internal parameters for options.
    //
    static int        show_mem;
    static int        show_timings;
    static int        verbose;
    static int        debug_fillpatch; // default=0
    static int        check_nan;
    static int        fab_verbose; // 1=curv 2=error heat 3=both
    static int        output_drop_distribution;
    static int        extend_pressure_into_solid; // default=0
     // default=diagonal length of the domain tangent to the
     // gravity direction.
    static Real       gravity_reference_wavelen;
    static int        NUM_SCALARS;      // number of non-velocity components

    static int include_viscous_heating; // default=0
     //smooth_type: 0=GSRB 1=weighted Jacobi 2=ILU
    static int smooth_type; // Lp.smooth_type def=2 GSRB
    static int bottom_smooth_type; // Lp.bottom_smooth_type def=2 GSRB
    static int global_presmooth; // default=2
    static int global_postsmooth; // default=2
    static int use_mg_precond_in_mglib; //def=1 Lp.use_mg_precon_in_mglib
     //def=0.001 Lp.bottom_bottom_tol_factor
    static Real bottom_bottom_tol_factor;

    static Vector<int>  FSI_touch_flag;
     // 0=> u=u_solid if phi_solid>=0
     // 1=> u=u_solid_ghost if phi_solid>=0
     // 2=> GNBC (u=u_solid_ghost if phi_solid>=0)
    static Vector<int> law_of_the_wall;  // def=0
    static Vector<Real> wall_model_velocity;  // def=0.0
    static Vector<int> interface_mass_transfer_model;  // def=0
    static Real        wall_slip_weight; // def=0.0

    static int FSI_interval; //default=1

    static int CTML_FSI_init;

    static Vector<int>  shock_timestep;

     // vel, pres, num_state_material x num_materials, 
     // ngeom_raw x num_materials, error ind
    static int State_Type; 
     // mac vel
    static int Umac_Type;
    static int Vmac_Type;
    static int Wmac_Type;
     // LS 1..num_materials, LS_slope sdim x num_materials
    static int LS_Type;
     // -(pnew-pold)/(rho c^2 dt) + dt mdot/vol
    static int DIV_Type;
    static int Solid_State_Type;
    static int Tensor_Type;
    static int TensorX_Type;
    static int TensorY_Type;
    static int TensorZ_Type;
    static int NUM_STATE_TYPE;

    static Vector<Real> tension_min;
    static Vector<Real> tension_T0;
    static Vector<Real> tension_slope; //should be non-positive
    static Vector<Real> tension; // 1-2 1-3 1-4 2-3 2-4 3-4
    static Vector<Real> tension_init; // 1-2 1-3 1-4 2-3 2-4 3-4
    static Real unscaled_min_curvature_radius; //default: 2.0
    static Vector<Real> prefreeze_tension; // 1-2 1-3 1-4 2-3 2-4 3-4

     // (dir,side)  (1,1),(2,1),(3,1),(1,2),(2,2),(3,2)
    static Vector<Real> outflow_velocity_buffer_size;

    static Vector<Real> cap_wave_speed; // num_interfaces components

    static Real real_number_of_cells;

     // default=1.0e+6
    static Real mglib_max_ratio; 
     // default=0.0
    static Real min_interior_coeff;

    static int custom_nucleation_model;

    static int FD_curv_interp; // def==1
    static int vof_height_function; // def==1

     // index: 1..num_materials
    static Vector<Real> cavitation_pressure; 
     // index: 1..num_materials
    static Vector<Real> cavitation_vapor_density; 
    static Vector<Real> cavitation_tension; 


      // 1-2 means "material 1 is converted into material 2"
      // 1-2 1-3 1-4 2-3 2-4 3-4
      // 2-1 3-1 4-1 3-2 4-2 4-3
      // if hardwire_Y_gamma and hardwire_T_gamma greater than zero
      // then use these values at the interface changing phase.
    static Vector<Real> hardwire_Y_gamma; 
    static Vector<Real> hardwire_T_gamma; 
    static Vector<Real> saturation_temp; 
     // "p_boil" in Dodd and Ferrante
     // default: 1.0e+6
    static Vector<Real> reference_pressure; 
    static Vector<Real> accommodation_coefficient; 
    static Vector<Real> saturation_temp_curv; 
    static Vector<Real> saturation_temp_vel; 
    static Vector<Real> saturation_temp_min;  //aka T_I_min
    static Vector<Real> saturation_temp_max;  //aka T_I_max
    static Vector<Real> nucleation_pressure; 
    static Vector<Real> nucleation_pmg; 
    static Vector<Real> nucleation_mach; 
    static Vector<Real> nucleation_temp; 
    static Real nucleation_period; 
    static int pos_sites_random_flag; //default=1
    static int n_sites;
    static Vector<Real> pos_sites;

     // if microlayer_size>0 and LS_substrate>0 in 5x5x5 stencil,
     // then the microlayer model is activated.
     // subroutine get_vel_phasechange returns the larger of the microlayer
     // model rate and the conventional rate. 
    static Vector<Real> microlayer_size; // minimum transition length
    static Vector<Real> macrolayer_size; // maximum transition length

     // minimum transition length for heat conduction.
     // e.g. k_heater=k * max(1,dx/thermal_microlayer_size)
    static Vector<Real> thermal_microlayer_size;
     // minimum transition length for viscosity. 
     // e.g. visc+viscconst_eddy_wall<visc * max(1,dx/shear_microlayer_size)
     // e.g. (u_in - u_law_of_wall)/dx<u_in * max(1,dx/shear_microlayer_size)
    static Vector<Real> shear_microlayer_size;
     // e.g. k_model<k * dx/buoyancy_microlayer_size
    static Vector<Real> buoyancy_microlayer_size;
     // grad T \approx (T_probe-T_I)/phasechange_microlayer_size if material
     // "im" macro scale probe cannot be found.
    static Vector<Real> phasechange_microlayer_size;

     // if the contact line size>max_size, then
     // multiply the microlayer model heat flux by 
     // max_size/size
     // default is 0.0
    static Vector<Real> max_contact_line_size;

    static Vector<Real> microlayer_angle; // radians
    static Vector<int>  microlayer_substrate; // index of substrate

     // if freezing_model==0:
     // given im1,im2 pair:
     // if microlayer_substrate(im1)>0 and
     //    microlayer_temperature_substrate(im1)>0.0 and
     //    solidheat_flag==0 (diffuse in solid) then
     //  T=microlayer_temperature_substrate at im2/substrate boundary.
     //
     // if microlayer_substrate(im1 or im2)>0 and
     //    microlayer_temperature_substrate(im1 or im2)>0.0 and
     //    solidheat_flag==2 (Neumann at solid/fluid interface) then
     //  grad T dot n=(microlayer_temperature_substrate-TSAT)/
     //    macrolayer_size at (im1 or im2)/substrate boundary.

    static Vector<Real> microlayer_temperature_substrate; 

      // 1-2 1-3 1-4 2-3 2-4 3-4
      // 2-1 3-1 4-1 3-2 4-2 4-3
      // latent_heat<0 if condensation or solidification
      // latent_heat>0 if boiling or melting
    static Vector<Real> latent_heat; 
    static Vector<Real> latent_heat_min; 
    static Vector<Real> latent_heat_T0; 
    static Vector<Real> latent_heat_slope; //should be non-positive 

    static Vector<Real> reaction_rate; 
     // freezing_model:
     // 0=sharp interface model (TSAT dirichlet BC)
     // 1=source term model (single equation for T with source term)
     // 2=hydrate model 
     //    a) V=V(P,T,RHO,C)
     //    b) single equation for T with source term
     //    c) single equation for C with source term.
     //       Coefficient for diffusivity "D" is zero at junction
     //       between F_water=1 and F_water<1.
     // 5=evaporation
     // 1-2 1-3 1-4 2-3 2-4 3-4
     // 2-1 3-1 4-1 3-2 4-2 4-3
    static Vector<int> freezing_model;
     // 0=Palmore and Desjardins (same as Villegas et al)
     // 1=Tanasawa  2=Schrage  3=Kassemi
    static Vector<int> Tanasawa_or_Schrage_or_Kassemi; 
     // index: 1..num_materials
     // value: 1..num_species_var
    static Vector<int> rigid_fraction_id;
     // index: 1..2*num_interfaces
     // value: 1..num_species_var
    static Vector<int> mass_fraction_id; 
     // index: 1..num_species_var
     // value: 1..num_materials
    static Vector<int> spec_material_id_LIQUID; 
    static Vector<int> spec_material_id_AMBIENT; 
     // 0 - distribute to the destination material 
     //     V=u_src dot n + mdot/rho_src
     // 1 - distribute to the source material
     //     V=u_dst dot n + mdot/rho_dst
    static Vector<int> distribute_from_target;

     // 0 - mdot concentrated at full cells near interface on 1-side
     // 1 - mdot distributed evenly to all full cells; cellvol weight
     // 2 - mdot distributed evenly to all full cells; constant weight
    static Vector<int> distribute_mdot_evenly; // 1..2*num_interfaces
     // 0 - sum mdot <>0   
     // 1 - distribute sum -mdot to the source
     // -1 - distribute sum -mdot to the dest
    static Vector<int> constant_volume_mdot; // 1..2*num_interfaces
    static Vector<int> constant_density_all_time; // 1..num_materials, def=1
    static int is_phasechange;
    static int ngrow_distance;  // should always be 4
    static int ngrow_make_distance;  // should always be 3

    static int num_elements_blobclass;
    static int BLB_MATRIX;
    static int BLB_RHS;
    static int BLB_VEL;
    static int BLB_INT_MOM;
    static int BLB_ENERGY;
    static int BLB_MASS_VEL;
    static int BLB_VOL;
    static int BLB_CEN_INT;
    static int BLB_CEN_ACT;
    static int BLB_PERIM;
    static int BLB_PERIM_MAT;
    static int BLB_TRIPLE_PERIM;
    static int BLB_CELL_CNT;
    static int BLB_CELLVOL_CNT;
    static int BLB_MASS;
    static int BLB_PRES;

    static Vector<int> im_solid_map; //nparts components, in range 0..num_materials-1
    static Vector<int> im_elastic_map;

     // 0=volume fraction  1=mass fraction 2=impedance fraction
    static int pressure_select_criterion; 

    static int  post_init_pressure_solve;
    static Real denfact;
    static Real velfact;

    static Real xblob;
    static Real yblob;
    static Real zblob;
    static Real radblob;

    static Real xblob2;
    static Real yblob2;
    static Real zblob2;
    static Real radblob2;

    static Real xblob3;
    static Real yblob3;
    static Real zblob3;
    static Real radblob3;

    static Real xblob4;
    static Real yblob4;
    static Real zblob4;
    static Real radblob4;

    static Real xblob5;
    static Real yblob5;
    static Real zblob5;
    static Real radblob5;

    static Real xblob6;
    static Real yblob6;
    static Real zblob6;
    static Real radblob6;

    static Real xblob7;
    static Real yblob7;
    static Real zblob7;
    static Real radblob7;

    static Real xblob8;
    static Real yblob8;
    static Real zblob8;
    static Real radblob8;

    static Real xblob9;
    static Real yblob9;
    static Real zblob9;
    static Real radblob9;

    static Real xblob10;
    static Real yblob10;
    static Real zblob10;
    static Real radblob10;

    static Real xactive;
    static Real yactive;
    static Real zactive;
    static Real ractive;
    static Real ractivex;
    static Real ractivey;
    static Real ractivez;

    static int adapt_quad_depth;
    static int probtype;

    static int perturbation_on_restart;
     // sin(2 pi kx/L)
    static int perturbation_mode;  // number of wavelengths each dir.
     // a fraction of delta T or the radial velocity rmax omega.
    static Real perturbation_eps_temp;
    static Real perturbation_eps_vel;

    static int step_through_data; //def=0

    static int visual_tessellate_vfrac; // def=0
    static int visual_revolve; //def=0
    static int visual_output_raw_State_Type; //def=0
    static int visual_output_raw_mac_Type; //def=0
    static int visual_phase_change_plot_int; //def=0
    static int visual_buoyancy_plot_int; //def=0
    static int visual_divergence_plot_int; //def=0
    static int visual_WALLVEL_plot_int; //def=0
    static int visual_drag_plot_int; //def=0
    //0=tecplot nodes
    //1=plt file cells
    //2=tecplot cells (piecewise constant reconstruction).
    static int visual_nddata_format; //def=0

    static Vector<int> visual_ncell;
    static int visual_compare;
 
    static int slice_dir;  // 0,..,sdim-1
    static Vector<Real> xslice;

    static int nblocks;
    static Vector<Real> xblocks;
    static Vector<Real> yblocks;
    static Vector<Real> zblocks;
    static Vector<Real> rxblocks;
    static Vector<Real> ryblocks;
    static Vector<Real> rzblocks;

    static int tecplot_max_level; //def=max_level
    static int max_level_for_use;

    static int ncoarseblocks;
    static Vector<Real> xcoarseblocks;
    static Vector<Real> ycoarseblocks;
    static Vector<Real> zcoarseblocks;
    static Vector<Real> rxcoarseblocks;
    static Vector<Real> rycoarseblocks;
    static Vector<Real> rzcoarseblocks;

    static int ncomp_sum_int_user1; //def=0
    static int ncomp_sum_int_user2; //def=0
    static int ncomp_sum_int_user12; //def=0

    static Vector<Real> compressible_dt_factor;  //def=1.0
    static int disable_advection;
    static int disable_pressure_solve;

    static Vector<Real> etaL; //def=0
    static Vector<Real> etaP; //def=0
    static Vector<Real> etaS; //def=0
    static Vector<Real> concentration; // def=0

    static Vector<Real> Carreau_alpha; // def=1
    static Vector<Real> Carreau_beta; // def=0
    static Vector<Real> Carreau_n; // def=1
    static Vector<Real> Carreau_mu_inf; // def=0
    static Vector<int> shear_thinning_fluid; // def=0

    static Vector<Real> elastic_time; // viscoelastic parameter

     // 0 => viscoelastic FENE-CR material  
     //    force= (div mu_p(f(A)/lambda) Q)   Q=A-I  mu_p="elastic_viscosity"
     //                                       (no shear thinning case)
     //    relaxation term=(f(A)/lambda)Q     f(A)=1/(1-trac(A)/L^2)
     //     updating Q:
     //       (i) CISL advection Q^advect=Q^n(x-V dt)
     //       (ii) A^advect=Q^advect+I
     //       (iii) X=I+dt gradu
     //       (iv)  Q^{n+1}=X A^adv X^T -I=
     //           Q^advect+dt (2D) + dt gradu Q + dt Q gradu^T+O(dt^2)
     //       (v) lambda'=lambda * (1-tr(A)/L^2)
     //       (vi) Q^n+1=lambda' Q^n+1/(lambda'+dt)
     // 1 => Oldroyd-B
     //       Oldroyd-B=lim_{L\rightarrow\infty} FENE_CR
     //       (i) lambda'=lambda
     // 3=> incremental elastic model 
     //   DS/DT=2 (D0-Dp) - (SW-WS)
     //     mu=Lame coefficient (bulk modulus?)
     //     D0=D-tr(D)Id/DIM 
     //       =D if incompressible
     //     W=(1/2)(grad V - grad V^T)    W^T=-W
     //     updating S:
     //       (i) CISL advection S^advect=S^n(x-V dt)
     //       (ii) A^advect=S^advect+I+dt 2D
     //       (iii) X=I+dt W
     //       (iv) S^{n+1}=X A^advect X^T-I=
     //         S^advect+dt (2D)+dt W S+dt S W^T + O(dt^2)
     // 4=> elastic model handled in the coupler (Lagrangian).
     // 5=> FENE-P 
     //       force=div mu_p (1/lambda)(f(A)A - I)
     //       relaxation term=(1/lambda)(f(A)A-I)
     //       (v) lambda'=lambda * (1-tr(A)/L^2)
     //       (vi) ofs=I*tr(A)/L^2
     //       (vi) Q_t = -(1/lambda')(Q+ofs)  
     //       (vii) Q^{n+1}-Q^star=-(dt/lambda')(Q^{n+1} + ofs)
     //       (viii) (1+dt/lambda')Q^{n+1}=Q^star - dt * ofs/lambda'
     //       (viv) Q^{n+1}=(lambda'/(lambda'+dt))*(Q^star-dt*ofs/lambda')
     //       (x) for incompressible flow, source term is equivalent to 
     //           FENE-CR source term which is tau=Q/lambda'.
     // 6=> Linear PTT
     //       (v) lambda'=lambda
     //       (vi) Q_t = -(1/lambda)(Q+Tr(Q)Q/L^2)
     //       (vii) lambda''=lambda*(1/(1+Tr(Q)/L^2))
     //       (vii) Q^{n+1}-Q^star=-(dt/lambda'')Q^{n+1}
     //       (viii) Q^{n+1}=lambda'' Q^n+1/(lambda''+dt)
     // 7=> Neo-Hookean (using Left Cauchy Green tensor B=F F^{T}
     //     Xia, Lu, Tryggvason 2018
    static Vector<int> viscoelastic_model; // def=0
    static Vector<int> viscosity_state_model; // def=0
    static Vector<int> les_model;

    static int transposegradu; //def=0

    static Vector<Real> elastic_viscosity; // viscoelastic parameter
     // for freezing:
     // 1. nonlinear advection
     // 2. calculate mass flux at the liquid/ice interface
     // 3. update the liquid/ice interface according to the mass flux
     //    (including expansion e.g. water freezing)
     // 4. energy equation
     // 5. viscosity force (viscosity coefficient >> 1 in the ice)
     // 6. (a) pressure gradient:
     //      t^{n+1}=(n+1)\Delta t  if uniform \Delta t.
     //      u^{n+1}=
     //         u^{viscosity}-dt F^{surface tension}-grad p^{n+1}/rho^{n+1}
     //      div u^{n+1} = 0
     //
     //      NOTE:
     //      (i) project the ice velocity u^{viscosity} to a rigid body
     //          velocity u^{rigid}.  u=u^{rigid} in the ice.
     //


    static Vector<Real> static_damping_coefficient; // dv/dt = -mu v

    // (1/L)   eps=0.01
    // viscoelastic_model=0  FENE CR   trac(A)<L^2 lambda(A)>eps/L^2
    // viscoelastic_model=1  Oldroyd B trac(A)<inf lambda(A)>eps/L^2
    // viscoelastic_model=5  FENE P    trac(A)<L^2 lambda(A)>eps/L^2
    // viscoelastic_model=6  linearPTT trac(A)<inf lambda(A)>eps/L^2
    static Vector<Real> polymer_factor; // (1/L)

    static int adv_dir;
    static Real adv_vel;
    static int axis_dir;

    static Vector<int> mof_ordering;

    static Real rgasinlet;
    static Real slipcoeff;
    static Real vinletgas;
    static Real twall;
    static Real walltemp;
     // 1=EILE (default), -1=Weymouth Yue
     // 2=always EI   3=always LE
    static int EILE_flag; 
    static int krylov_subspace_max_num_outer_iter;
    static Real projection_pressure_scale;
    static Real projection_velocity_scale;

    static Vector<Real> vorterr; // 0 to disable
    static int pressure_error_flag;  // 0 (check mag, default) 1=check diff
    static Vector<Real> pressure_error_cutoff; // 0 to disable
    static Vector<Real> temperature_error_cutoff; // 0 to disable
    static Real temperature_source;
    static Vector<Real> temperature_source_cen;
    static Vector<Real> temperature_source_rad;
    static Vector<Real> tempcutoff;
    static Vector<Real> tempcutoffmax;
    static Vector<Real> tempconst;
    static Vector<Real> initial_temperature;
    static Real initial_temperature_diffuse_duration; // def=0.0
    static Vector<Real> density_floor;  // def=0.0
    static Vector<Real> density_ceiling;  // def=1.0e+20
    static Vector<Real> molar_mass;  // def=1
    static Vector<Real> denconst;
     // if stokes_flow==1, do not include u dot grad u term.
    static int stokes_flow;  // def=0
     // if cancel_advection==1, do not include u dot grad S terms for any
     // S corresponding to velocity, interface, or state variables.
    static int cancel_advection;  // def=0
    static Vector<Real> denconst_interface_added;

    static Vector<Real> stiffPINF;
    static Vector<Real> prerecalesce_stiffCP; // def=4.1855E+7
    static Vector<Real> prerecalesce_stiffCV; // def=4.1855E+7
    static Vector<Real> stiffCP; // def=4.1855E+7
    static Vector<Real> stiffCV; // def=4.1855E+7
    static Vector<Real> stiffGAMMA;

    static int uncoupled_viscosity;

    static Real centrifugal_force_factor; //Eady disables centrifugal force.
     
     //Du/Dt=-grad (p-rho0 g dot z)/rho0 - g DrhoDT (T-T0) 
     //DrhoDT has units of 1/(Degrees Kelvin)
    static Vector<Real> DrhoDT; // default 0.0

     //1=>rho=rho(T)
     //2=>Du/Dt=-grad (p-rho0 g dot z)/rho0 - g DrhoDT (T-T0) 
    static Vector<int> override_density; // default 0 
    static Vector<Real> viscconst;
    static Vector<Real> viscconst_eddy_wall;
    static Vector<Real> viscconst_eddy_bulk;
    static Vector<Real> heatviscconst_eddy_wall;
    static Vector<Real> heatviscconst_eddy_bulk;
    static Vector<Real> prerecalesce_viscconst;
    static Vector<Real> viscconst_interface;
    static Vector<Real> heatflux_factor;
    static Vector<Real> heatviscconst;
    static Real heatviscconst_max;
    static Real heatviscconst_min;
    static Vector<Real> prerecalesce_heatviscconst;
    static Vector<Real> heatviscconst_interface;
      // num_species * num_materials size.
    static Vector<Real> speciesreactionrate;
      // num_species * num_materials size.
    static Vector<Real> speciesconst;
      // num_species * num_materials size.
    static Vector<Real> speciesviscconst;
      // num_species * num_interfaces
    static Vector<Real> speciesviscconst_interface;
      // num_species
    static Vector<Real> species_molar_mass;  // def=1

     //0 incomp; material_type_evap needed for the Kassemi model.
    static Vector<int> material_type_evap; 
    static Vector<int> material_type_lowmach; // 0 incomp   
    static Vector<int> material_type_visual; // 0 incomp   
    static Real wait_time;
    static Real advbot;

// mg.bot_atol
// mg.visc_bot_atol
// mg.thermal_bot_atol
//
    static Real save_atol_b;
    static Real save_mac_abs_tol;
    static Real save_min_rel_error;

    static Real minimum_relative_error;
    static Real diffusion_minimum_relative_error;
    static Real mac_abs_tol;
    static Real visc_abs_tol;
    static Real thermal_abs_tol;

    static int  multilevel_maxcycle;
    static int  multilevel_restart_period;
    static int  initial_project_cycles;
    static int  initial_viscosity_cycles;
    static int  initial_thermal_cycles;
     // 0 - MGPCG  1-PCG 2-MINV=I
    static int  project_solver_type;
    static int  initial_cg_cycles;
    static int  num_divu_outer_sweeps;

    static Real divu_outer_sweeps_preconditioner;
    static Real activated_divu_preconditioner;

    static int  debug_dot_product;

    static Vector<Real> NS_coflow_Z;
    static Vector<Real> NS_coflow_R_of_Z;

    static Vector<Real> NS_DRAG_integrated_quantities;
    static Vector<int> NS_DRAG_integrated_quantities_sweep;

    static Vector<Real> NS_sumdata;
    static Vector<int> NS_sumdata_type;
    static Vector<int> NS_sumdata_sweep;

};

}/* namespace amrex */

#endif /*_NavierStokes_H_*/
