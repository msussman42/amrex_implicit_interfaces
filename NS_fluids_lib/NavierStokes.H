#ifndef _NavierStokes_H_
#define _NavierStokes_H_

#include <AMReX_Geometry.H>
#include <AmrLevel.H>
#include <LevelBld.H>
#include <ABecLaplacian.H>

namespace {
    constexpr double NS_PI = 3.14159265358979323846264338327950288;
}

namespace amrex{

// assume 3d for face declarations.

#define OFF_DIAG_CHECK_MF 0  // nsolveMM components
#define DIAG_REGULARIZE_MF 1

#define ADVECT_REGISTER_FACE_MF 2

#define ADVECT_REGISTER_MF 5
#define DRAG_MF 6

#define CELL_ELASTIC_FORCE_MF 7

#define REGISTER_MARK_MF 8

#define HOLD_VELOCITY_DATA_MF 9

#define MAGTRACE_MF 10

#define HOLD_ADVECT_STATE_CELL_MF 11

#define HOLD_ADVECT_STATE_FACE_MF 12

#define RECENT_VELOCITY_MF 15
#define MASK_RESIDUAL_MF 16
#define MASK_DIV_RESIDUAL_MF 17

#define FACETENSOR_MF 18
#define CELLTENSOR_MF 19
#define LSTENSOR_MF 20
#define MASKSOLIDTENSOR_MF 21

#define BEFORE_VISCOUS_VELOCITY_MF 22

#define CONSERVE_FLUXES_MF 23

#define ALPHANOVOLUME_MF 26
#define ALPHACOEF_MF 27

#define BXCOEFNOAREA_MF 28
#define BXCOEF_MF 31

#define VOLUME_MF 34

#define AREA_MF 35

#define MASKCOEF_MF 38

#define HISTORY_MF 39

#define SLOPE_RECON_MF 40
#define MASKSEM_MF 41

#define FACE_VAR_MF 42

#define dt_gradp_over_rho_cell_MF 45

// I-scheme, conduction, viscosity, div(up), grad p
#define stableF_MF 46
#define stableF_GP_MF 47

// used in update_SEM_forces
#define GP_DEST_FACE_MF 50

#define SATURATION_TEMP_MF 54

#define SEM_FLUXREG_MF 55

#define LOCAL_FACETENSOR_MF 56
#define LOCAL_CELLTENSOR_MF 57

#define VISCOTEN_MF 58
#define CELLFRAC_MF 59
#define BURNING_VELOCITY_MF 60
#define FACEFRAC_MF 61
#define STENCIL_MF 62
#define ORIGDIST_MF 63

#define MDOT_LOCAL_MF 64

#define HYDROSTATIC_PRESDEN_MF 65
#define MDOT_MF 66
#define VISCHEAT_MF 67
#define FSI_MF 68
#define VISCHEAT_SOURCE_MF 69
#define JUMP_STRENGTH_MF 70
#define CELL_DEDT_MF 71
#define CELL_DEN_MF 72
#define CELL_SOUND_MF 73
#define CELL_VISC_MF 74
#define CELL_VISC_MATERIAL_MF 75
#define LS_NRM_FD_MF 76
#define HOLD_LS_DATA_MF 77
#define DIFFUSIONRHS_MF 78
#define LEVELPC_MF 79
#define VOF_PREV_TIME_MF 80

#define CELL_CONDUCTIVITY_MATERIAL_MF 81

#define MAC_VELOCITY_MF 82

#define MASK_NBR_MF 85
#define DEN_RECON_MF 86
#define LS_RECON_MF 87
#define OUTER_ITER_PRESSURE_MF 88
#define ONES_MF 89
#define POLDHOLD_MF 90

#define UMACSTAR_MF 91
#define UMAC_MF 94
#define MAC_TEMP_MF 97
#define GRADPEDGE_MF 100
#define PEDGE_MF 103

#define REGISTER_MARK_MAC_MF 106

#define FACE_WEIGHT_MF 109

#define dt_gradp_over_rho_face_MF 112

#define POTENTIAL_FORCE_EDGE_MF 115

#define S_NEW_GROW_MF 119
#define INTERIOR_FLAG_MF 120

#define LS_NRM_FD_GNBC_MF 121

#define TYPE_MF 122
#define COLOR_MF 123
#define PRES_MF 124
#define PRESPC_MF 125
#define PRESPC2_MF 126
#define RESID_MF 127
#define OUTER_RESID_MF 128
#define CGRESID_MF 129
#define Z_MF 130
#define P_MF 131

#define LS_NRM_CP_MF 132

#define MAC_PHI_CRSE_MF 133
#define MAC_RHS_CRSE_MF 134
#define OUTER_MAC_RHS_CRSE_MF 135
#define OUTER_MAC_PHI_CRSE_MF 136
#define REGISTER_CURRENT_MF 137

#define FACETEST_MF 138

#define DINVERSE_MF 139
#define TENSOR_RHS_MF 140

#define HOOP_FORCE_MARK_MF 142
#define BOUSSINESQ_TEMP_MF 143

#define VOF_RECON_MF 144

#define CGRESID_ZERO_MF 146
#define V_MF 147
#define Y_MF 148
#define S_MF 149
#define T_MF 150

#define P_SOLN_MF 151

#define FD_NRM_ND_MF 152
#define FD_CURV_CELL_MF 153

#define FACEFRAC_SOLVE_MM_MF 155
#define FACEFRAC_MM_MF 158
#define CELLFRAC_MM_MF 159

#define VELADVECT_MF 161
#define UMACOLD_MF 162

#define SWEPT_CROSSING_MF 166

#define CELL_VOF_MF 167
#define MACDIV_MF 168
#define DIST_CURV_MF 171

#define bicg_FLUXES_MF 172
#define bicg_FLUX_COEFF_MF 175
#define bicg_G_COEF_MF 180
#define bicg_G_MF 181
#define bicg_U0_MF 182
#define bicg_V0_MF 183
#define bicg_R1_MF 186
#define bicg_V1_MF 188
#define bicg_P1_MF 189
#define bicg_R0hat_MF 190
#define bicg_RHS_MF 192
#define bicg_Y_MF 193
#define bicg_Hvec_MF 194
#define bicg_S_MF 195
#define bicg_T_MF 196
#define bicg_AU_MF 198
#define bicg_RESID_MF 199

#define CELL_VELOCITY_MF 200

// nfluxSEM=sdim + 1
// nstate_SDC= nfluxSEM + 1 (thermal conduction) + sdim (viscosity) 
#define spectralF_MF 201
#define spectralF_GP_MF 202

#define delta_MF 205
#define delta_GP_MF 206

#define nodevel_MF 209
#define LSNEW_MF 210
#define donorflag_MF 211

#define L1_LEVELPC_MF 213

#define DIST_TOUCH_MF 214

#define VISCOTEN_ALL_MAT_MF 215

#define TENSOR_RECON_MF 216

#define DIV_SAVE_MF 218

#define AMRSYNC_VEL_MF 220
#define AMRSYNC_PRES_MF 223
#define COARSE_FINE_FLUX_MF 229

#define DELTA_CELL_VEL_MF 232

#define STATE_FOR_RESID_MF 233

#define LS_COLORSUM_MF 234
#define VEL_COLORSUM_MF 235
#define DEN_COLORSUM_MF 236

#define PRESSURE_SAVE_MF 237

#define MOM_DEN_MF 238

#define SAVE_MAC_PRESCRIBE_MF 239
#define SAVE_SNEW_PRESCRIBE_MF 242
#define SAVE_POLDHOLD_PRESCRIBE_MF 243

#define FSI_GHOST_MAC_MF 247

#define HISTORY_MAC_MF 250

#define MAX_GMRES_BUFFER 256

#define GMRES_BUFFER0_V_MF 253
#define GMRES_BUFFER_LAST_V_MF (GMRES_BUFFER0_V_MF+MAX_GMRES_BUFFER-1)

#define GMRES_BUFFER0_U_MF (GMRES_BUFFER_LAST_V_MF+1)
#define GMRES_BUFFER_LAST_U_MF (GMRES_BUFFER0_U_MF+MAX_GMRES_BUFFER-1)

#define GMRES_BUFFER0_Z_MF (GMRES_BUFFER_LAST_U_MF+1)
#define GMRES_BUFFER_LAST_Z_MF (GMRES_BUFFER0_Z_MF+MAX_GMRES_BUFFER-1)

#define GMRES_BUFFER_W_MF (GMRES_BUFFER_LAST_Z_MF+1)

#define XDISPLACE_MF (GMRES_BUFFER_W_MF+1)
#define XDISP_FLUX_MF (XDISPLACE_MF+1)
#define JUMP_STRENGTH_COMPLEMENT_MF (XDISP_FLUX_MF+3)
#define donorflag_complement_MF (JUMP_STRENGTH_COMPLEMENT_MF+1)
#define CURRENT_CELL_VEL_MF (donorflag_complement_MF+1)
#define XDMACOLD_MF (CURRENT_CELL_VEL_MF+1)

#define T_advect_MF (XDMACOLD_MF+3)
#define DTDt_MF (T_advect_MF+1)

#define TEMPERATURE_SAVE_MF (DTDt_MF+1)
#define TEMPERATURE_SMOOTH_MF (TEMPERATURE_SAVE_MF+1)

#define ONES_GROW_MF (TEMPERATURE_SMOOTH_MF+1)

#define TYPE_ONES_MF (ONES_GROW_MF+1)
#define COLOR_ONES_MF (TYPE_ONES_MF+1)

#define MAC_ELASTIC_FLUX_CC_MF (COLOR_ONES_MF+1)
#define MAC_ELASTIC_FLUX_XY_MF (MAC_ELASTIC_FLUX_CC_MF+1)
#define MAC_ELASTIC_FLUX_XZ_MF (MAC_ELASTIC_FLUX_XY_MF+1)
#define MAC_ELASTIC_FLUX_YZ_MF (MAC_ELASTIC_FLUX_XZ_MF+1)

// holds the cell centered interpolant of the MAC XDISP.
#define VISUAL_XDISP_MAC_CELL_MF (MAC_ELASTIC_FLUX_YZ_MF+1)

#define HOLD_VELOCITY_DATA_MAC_MF (VISUAL_XDISP_MAC_CELL_MF+1)

#define NEXT_MF (HOLD_VELOCITY_DATA_MAC_MF+3)

#define MAX_NUM_LOCAL_MF (NEXT_MF+3)

// State_Type:
// velocity, pressure, state x nmat, interface vars x nmat, error ind
// Solid_State_Type:
// velocity
// Tensor_Type:
// Tij x num_materials_viscoelastic + sdim

enum StateNames {

  DenVar=0, TemperatureVar, SpeciesVar };

// Q is stored in the code.  A=Q+I   A is sometimes called the
// "configuration tensor" 
// See: Stewart, Lay, Sussman, Ohta (and papers cited, and papers who
// cite this)
// See also earlier work involving Pushpendra Singh (NJIT).
#if (AMREX_SPACEDIM == 2)

enum TensorNames {T11=0,T12,T22,T33,NUM_TENSOR_TYPE};

enum VelocityNames {
  Xvel=0, Yvel };

enum MOF_Names  { 
   VoFVar=0,CenxVar,CenyVar,OrderVar,NxVar,NyVar,LSVar,NUM_MOF_VAR };

#else

enum TensorNames {T11=0,T12,T22,T33,T13,T23,NUM_TENSOR_TYPE};

enum VelocityNames {
  Xvel=0, Yvel, Zvel };

enum MOF_Names  { 
   VoFVar=0,CenxVar,CenyVar,CenzVar,OrderVar,
   NxVar,NyVar,NzVar,LSVar, NUM_MOF_VAR };

#endif



class NSBld
    :
    public LevelBld
{
    virtual void variableSetUp () override;
    virtual void variableCleanUp () override;
    virtual AmrLevel *operator() () override;
    virtual AmrLevel *operator() (Amr&            papa,
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
  int im;  // 1<=im<=nmat
};

class NavierStokes
    :
    public AmrLevel 
{
public:

    NavierStokes ();

    NavierStokes (Amr&            papa,
                  int             lev,
                  const Geometry& level_geom,
                  const BoxArray& bl,
                  const DistributionMapping& dmap_in,
                  Real            time);

    virtual ~NavierStokes ();

    virtual void restart (Amr&          papa,
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

    void deallocate_potential_forceALL();
    void process_potential_forceALL();
    void process_potential_force_face();
    void init_gravity_potential();

    void GetDrag(Vector<Real>& integrated_quantities,int isweep);
    void GetDragALL(Vector<Real>& integrated_quantities);

    void process_recalesce_dataALL(
     Vector<int> recalesce_material,
     Vector<Real>& recalesce_state_old,
     Vector<Real>& recalesce_state_new);
    void process_recalesce_data(
     Vector<int> recalesce_material,
     Vector<Real>& recalesce_state_old,
     Vector<Real>& recalesce_state_new,
     Vector<Real>& integrated_quantities,
     int isweep);
    void recalesce_temperature(int im_source);

    void volWgtSum(
     Vector<Real>& result,
     Vector<int>& sumdata_type,
     Vector<int>& sumdata_sweep,
     Vector<Real>& ZZ,Vector<Real>& FF,
     int dirx,int diry,int cut_flag,
     MultiFab* dragmf,int isweep);

    void Number_CellsALL(Real& rcells);

    void volWgtSumALL(
     int post_init_flag,
     Vector<Real>& result,
     Vector<int>& sumdata_type,
     Vector<int>& sumdata_sweep,
     Vector<Real>& ZZ,Vector<Real>& FF,
     int dirx,int diry,int cut_flag,
     int isweep);

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
		const std::string& caller_string,
		int data_id,
		int ncomp,
		int data_mf,
		int state_type_mf,
		int data_dir);

    virtual void writePlotFile (
     int do_plot,int do_slice,
     int SDC_outer_sweeps_in,int slab_step_in);

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
   Real estTimeStep (Real local_fixed_dt,int caller_id);

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
    void prepare_post_process(int post_init_flag);

    void init_particle_container(int append_flag);

    void accumulate_info_no_particles(int im);

    void accumulate_PC_info(
      AmrParticleContainer<N_EXTRA_REAL,0,0,0>& localPC_no_nbr,
      NeighborParticleContainer<N_EXTRA_REAL,0>& localPC_nbr,
      int nnbr);

    void assimilate_vel_from_particles(
      AmrParticleContainer<N_EXTRA_REAL,0,0,0>& localPC_no_nbr,
      NeighborParticleContainer<N_EXTRA_REAL,0>& localPC_nbr,
      int nnbr);

    //
    // Advance grids at this level in time.
    //

    void debug_memory();

    void allocate_mdot();

    void restore_active_pressure(int save_mf);

    void allocate_maccoefALL(int project_option,int nsolve,
		    int create_hierarchy);
    void allocate_maccoef(int project_option,int nsolve,
		    int create_hierarchy);
     
    void deallocate_maccoefALL(int project_option);
    void deallocate_maccoef(int project_option);

    void prepare_mask_nbr(int ngrow);
    void resize_mask_nbr(int ngrow);
    void prepare_displacement(int mac_grow);

    // VAHAB HEAT SOURCE
    void heat_source_term_flux_source();

    void stefan_solver_init(MultiFab* coeffMF,int adjust_temperature,
		    int project_option);

    void prepare_viscous_solver();

    void alloc_gradp_over_rho(int alloc_flag);

    void alloc_DTDtALL(int alloc_flag);
    void alloc_DTDt(int alloc_flag);

    void prepare_advect_vars(Real time);
    void delete_advect_vars();
    void exit_viscous_solver();

    void diffusion_heatingALL(int source_idx,int idx_heat);
    void diffusion_heating(int source_idx,int idx_heat);

    void explicit_fluxes(int state_idx,int scomp,int iupdate,MultiFab* Un);

protected:

    void SDC_setup_step();
    void SDC_setup();

    void debug_ngrow(int idxMF,int ngrow,int counter);
    void debug_ixType(int idxMF,int grid_type,int counter);
    void debug_ixType_raw(MultiFab* mf,int grid_type,int counter);
    void debug_boxArray(MultiFab* mf,int grid_type,int counter);

    void avgDown_tag_localMF(int idxMF);
    void level_avgDown_tag(MultiFab& S_crse,MultiFab& S_fine);

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
    void avgDown_localMF(int idxMF,int scomp,int ncomp,
     int spectral_override);
    void Plus_localMF(int idx_dest,Real val,
     int dcomp,int ncomp,int ngrow);
    void Mult_localMF(int idx_dest,int idx_source,
     int scomp,int dcomp,int ncomp,int ngrow);
    void Copy_localMF(int idx_dest,int idx_source,
     int scomp,int dcomp,int ncomp,int ngrow);
    void minus_localMF(int idx_dest,int idx_source,
     int ncomp,int ngrow);

      // grid_type=-1 ... 5
    void new_localMF(int idx_MF,int ncomp,int ngrow,int grid_type);
      // grid_type=-1 ... 5
    void new_localMF_if_not_exist(int idx_MF,int ncomp,int ngrow,
		    int grid_type);

    //MAC_state_idx==Umac_Type or (mac velocity)
    //MAC_state_idx==XDmac_Type   (mac displacement)
    //dir=0..sdim-1
    void getStateMAC_localMF(int MAC_state_idx,
     int idx_MF,int ngrow,int dir,
     int scomp,int ncomp,Real time);

    void getStateDen_localMF(int idx_MF,int ngrow,Real time);
    void getStateDist_localMF(int idx_MF,int ngrow,Real time,int caller_id);
    void resize_levelsetLO(int ngrow,int idx);
    void allocate_levelsetLO(int ngrow,int idx);
    void allocate_levelsetLO_ALL(int ngrow,int idx);

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

    void getStateVISC();
    void getStateVISC_ALL();

    void getStateCONDUCTIVITY();
    void getStateCONDUCTIVITY_ALL();

    void delete_localMF(int idx_MF,int ncomp); 
    void delete_localMF_if_exist(int idx_MF,int ncomp); 

    void MaxAdvectSpeedALL(Vector<Real>& dt_min,
        Real* vel_max_estdt,
	Real& vel_max_cap_wave,int caller_id);
    void MaxAdvectSpeed(Vector<Real>& dt_min,
        Real* vel_max_estdt,
        Real& vel_max_cap_wave,int caller_id);
    void DumpProcNum();

      // interpolates coarsedata from level-1 to level.
    void AllinterpScalarMAC(
      MultiFab* coarsedata,MultiFab* finedata,
      BoxArray& cgridscen,BoxArray& fgridscen,
      DistributionMapping& fdmap,
      MultiFab* cdiagsing,MultiFab* fdiagsing,
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
    void ns_reconcile_d_num(int caller_id);

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

    static void read_params ();   // Read input file
    static void read_geometry (); // Set metrics

    static void override_enable_spectral(int enable_spectral_in);
    static void override_enable_spectralGHOST(
        int scomp,int ncomp,int enable_spectral_in);
    static void override_LS_HO(int Interp_LO); //0=use normals 1=PC

    static void set_tensor_extrap_components(
      int coord,std::string postfix,int indx,int ibase_tensor);
    static void set_tensor_extrap_components_main(
      int coord,std::string postfix,int indx);

    void sum_integrated_quantities (int post_init_flag);

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

    void correct_xdisplace_with_particles();

    void tensor_advection_update();
    void tensor_advection_updateALL();

    void level_init_icemask();
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

    void build_NRM_FD_MF(int fd_mf,int ls_mf, int ngrow);

    void level_phase_change_convert(
      int im_outer,int im_opp_outer,
      int i_phase_change,int n_phase_change);

    void level_phase_change_convertALL();

    void level_phase_change_rate(Vector<blobclass> blobdata,
       int color_count,int nucleation_flag);
    void level_phase_change_rate_extend();
   
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

    void init_gradu_tensorALL(
     int idx,
     int do_alloc,
     int idx_cell,
     int idx_face,
     int simple_AMR_BC_flag_viscosity);

    void make_viscoelastic_tensor(int im);
    void make_viscoelastic_tensorALL(int im);

    void make_viscoelastic_tensorMAC(int im,int interp_Q_to_flux,
	int flux_mf,int flux_grid_type,int fill_state_idx);
    void make_viscoelastic_tensorMACALL(int im,int interp_Q_to_flux,
        int flux_mf,int flux_grid_type,
        int fill_state_idx);

    void make_viscoelastic_heating(int im,int idx);

    // VAHAB HEAT SOURCE
    void make_heat_source();

    void add_perturbation();
    void make_SEM_delta_force(int project_option);

    void getStateMOM_DEN(int idx,int ngrow,Real time);

    void check_grid_places();

    void SOD_SANITY_CHECK(int id);

    void show_norm2(MultiFab* mf,int scomp,int ncomp,int id);
    void show_norm2_id(int mf_id,int id);
     // mf_check_inf_bounds must be 1.
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
      int dir,int id,
      Real warning_cutoff);

    void FillBoundaryTENSOR(MultiFab* mf,int sc,int dir);

    void check_for_NAN_TENSOR_base(int datatype,MultiFab* mf,
      int sc,int dir,int id);
    void check_for_NAN_TENSOR(int datatype,MultiFab* mf,int id);
    void check_for_NAN(MultiFab* mf,int id);

    void split_scalar_advection(); 
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

     // update_flag=0 do not update
     // update_flag=1 update S_old
     // update_flag=2 update S_new
    void VOF_Recon(int ngrow,
      Real time,int update_flag,
      int init_vof_prev_time,int dest_mf);
    void VOF_Recon_ALL(int ngrow,Real time,int update_flag,
      int init_vof_prev_time,int dest_mf);
    void VOF_Recon_resize(int ngrow,int dest_mf);

    void make_physics_varsALL(int project_option,int post_restart_flag,
		    int caller_id);

    void allocate_physics_vars();
    void make_physics_vars(int project_option);

    void solid_temperature();

    void nonlinear_advection();

    void advance_MAC_velocity(int project_option);

    void do_the_advance(Real timeSEM,Real dtSEM,int& advance_status);

    void correct_colors(
     int idx_color,int base_level,
     Vector<int> domaincolormap,int total_colors,
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
     int operation_flag, // =0 or 1
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
     int operation_flag, // =0 or 1
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
    void avgDownMacState(int MAC_state_idx,int spectral_override); 
    void avgDownEdge(int grid_type,MultiFab& S_crse,MultiFab& S_fine,
     int scomp,int ncomp,int spectral_override,int caller_id);
     // spectral_override==1 => order derived from "enable_spectral"
     // spectral_override==0 => always low order.
    void avgDownEdge_localMF(int idxMF,int scomp,int ncomp,
     int start_grid_type,int n_grid_type,
     int spectral_override,int caller_id);

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

     // particle_id==-1 => no info from particles
    void prescribe_solid_geometryALL(Real time,int renormalize_only,
      int local_truncate);
    void prescribe_solid_geometry(Real time,int renormalize_only);

    void move_particles(
      AmrParticleContainer<N_EXTRA_REAL,0,0,0>& localPC_no_nbr,
      NeighborParticleContainer<N_EXTRA_REAL,0>& localPC_nbr);

    void truncate_VOF(Vector<Real>& delta_mass_all);

    void getState_tracemag_ALL(int idx);
    void getState_tracemag(int idx);

    void output_triangles();

      //declared in: NavierStokes2.cpp
    void output_zones(
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
      int data_id,
      int data_dir,
      MultiFab* datamf,
      int ncomp,
      int& grids_per_level,
      BoxArray& cgrids_minusBA);

    void adjust_tolerance(Real& error0,Real& error0_max,
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
  int nsolve);

void mac_project_rhs(int project_option,
 int idx_mac_phi_crse,int idx_mac_rhs_crse,int nsolve);

void mac_update(MultiFab* mac_phi_crse,int project_option,int nsolve);

void diffuse_coefficients();

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
 int project_option,int nsolve);

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
        Real Ar_error_n,Real Ar_error0,
        Real rAr_error_n,Real rAr_error0,
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
     int bicgstab_num_outer_iterSOLVER,int nsolve);

    void updatevelALL(
     int project_option,
     int idx_mac_phi_crse,int nsolve);

    void extend_pressure_to_solidALL();

    void multiphase_project(int project_option);

    static void set_local_tolerances(int project_option);

    void multiphase_preconditioner(
     int project_option,int project_timings,
     int presmooth,int postsmooth,
     int idx_Z,int idx_R,int nsolve);

    void multiphase_GMRES_preconditioner(
     int gmres_precond_iter,
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

    void allocate_FACE_WEIGHT(
      int nsolve,
      int project_option);

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

    static void get_iten_cpp(int im1,int im2,int& iten,int nmat);
     // iten>=1  im1,im2>=1
    static void get_inverse_iten_cpp(int& im1,int& im2,int iten,int nmat);

     // unew^f=unew^f + beta * localMF[idx_velcell]
    void increment_face_velocity(
      int interp_option,
      int project_option,
      int idx_velcell,
      Real beta,
      Vector<blobclass> blobdata);

    void increment_face_velocityALL(
      int interp_option,
      int project_option,
      int idx_velcell,
      Real beta,
      Vector<blobclass> blobdata);

    void make_MAC_velocity_consistentALL();
    void make_MAC_velocity_consistent();

    void VELMAC_TO_CELLALL(
      int vel_or_disp, // =0 velocity; =1 displacement
      int dest_idx);

    void VELMAC_TO_CELL(
      int vel_or_disp, // =0 velocity; =1 displacement
      int dest_idx);

    void increment_KE(Real beta);
    void increment_KE_ALL(Real beta);

    void INCREMENT_REGISTERS_ALL(int source_mf,int caller_id);
    void INCREMENT_REGISTERS(int source_mf,int caller_id);
    void push_back_state_register(int idx_MF,Real time,int caller_id);
    void SET_STOKES_MARK(int idx_MF,int caller_id);
    void APPLY_VISCOUS_HEATING(int source_mf);
    void APPLY_VISCOUS_HEATINGALL(int source_mf);

    MultiFab* derive_EOS_pressure(Vector<int> local_material_type);
    void init_pressure_error_indicator();
    void init_advective_pressure(int project_option);

    void save_to_macvel_state(int idx_umac);

    void apply_cell_pressure_gradient(
     int project_option,
     int energyflag,
     int idx_pres,  
     int idx_umac,
     int idx_gpcell,
     int idx_divup);

    void update_SEM_delta_force(int project_option,
     int idx_gpmac,
     int idx_div,
     int update_spectral,int update_stable,
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

    void MAC_GRID_ELASTIC_FORCE(int im_elastic);

    MultiFab* getStateMAC(int MAC_state_idx,
	int ngrow,int dir,
	int scomp,int ncomp,Real time);

    void getStateDIV_ALL(int idx,int ngrow);
    void getStateDIV(int idx,int ngrow);

    void ADVECT_DIV_ALL();
    void ADVECT_DIV();

    void makeCellFrac(int tessellate,int ngrow,int idx);
    void makeFaceFrac(int tessellate,int ngrow,int idx,int do_face_decomp);
    void makeFaceTest(int tessellate,int ngrow,int idx);
    void ProcessFaceFrac(int tessellate,int idxsrc,int idxdst,
		    int ngrow_dest);
    void makeStateDistALL(int keep_all_interfaces);
    void makeStateDist(int keep_all_interfaces);
    void correct_dist_uninit();
    void makeStateCurv(int project_option,int post_restart_flag);
    void make_marangoni_force();

    void delete_array(int idx_localMF);

    void allocate_array(int ngrow,int ncomp,int grid_type,
      int idx_localMF);

    void allocate_independent_var(int nsolve,int idx);
    void allocate_rhs_var(int nsolve,int idx);

    void setVal_array(int ngrow,int ncomp,Real dataval,
      int idx_localMF);
    void setVal_localMF(int idx,Real dataval,
      int scomp,int ncomp,int ngrow);
    void mult_array(int ngrow,int ncomp,Real dataval,
      int idx_localMF);
    void copyALL(int ngrow,int ncomp,
      int scomp,int dcomp,int idx_dest,int idx_source);
    void minusALL(int ngrow,int ncomp,int idx_dest,int idx_source);

    void negateALL(int idx_localMF);
    void zeroALL(int ngrow,int ncomp,int idx_localMF);

    MultiFab* getStateDen(int ngrow,Real time);
    MultiFab* getStateDist(int ngrow,Real time,int caller_id);
    MultiFab* getStateDIV_DATA(int ngrow,int scomp,int ncomp,Real time);

    void putStateDIV_DATA(int scomp,int ncomp,int idx_MF);

    void getStateALL(int ngrow,Real time,int scomp,
      int ncomp,int idx_localMF);

    MultiFab* getStatePres(int ngrow,Real time);

    void resize_metrics(int ngrow);
    void metrics_data(int ngrow);
    void metrics_data_min_max(int caller_id);
    void metrics_dataALL(int ngrow);
    void metrics_data_min_max_ALL(int caller_id);

    MultiFab* maskfiner(int ngrow,Real tag,int clear_phys_boundary);

    NavierStokes& getLevel (int lev)
    {
        return *(NavierStokes*) &parent->getLevel(lev);
    }

    void init_regrid_history();

    void ns_header_msg_level(int FSI_operation,int FSI_sub_operation,
      Real time,Real dt,int iter);

    static int read_from_CAD();
    static int is_ice_matC(int im);
    static int is_singular_coeff(int im);
    static int is_FSI_rigid_matC(int im);
    static int CTML_FSI_flagC();
    static int CTML_FSI_matC(int im);
    static int FSI_material_exists();
    static int elastic_material_exists();

    void init_FSI_GHOST_MAC_MF_ALL(int caller_id);
    void init_FSI_GHOST_MAC_MF(int dealloc_history);

    void assimilate_state_data();

    void resize_FSI_MF();

    void FSI_make_distance(Real time,Real dt);

    void create_fortran_grid_struct(Real time,Real dt);

    void cpp_overridepbc(int homflag_in,int project_option_in);

    static int ns_is_rigid(int im);
    static int ns_is_lag_part(int im);

    void copy_velocity_on_sign(int partid);
    void Transfer_FSI_To_STATE(Real time);
    void build_moment_from_FSILS();

    void ctml_fsi_transfer_force();

protected:

    static int POLYGON_LIST_MAX;

    static bool ns_tiling; // default is false

    static int profile_debug;

    static int nfluxSEM;
    static int nstate_SDC;
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
    static BCRec temperature_phys_bc;
    static BCRec species_phys_bc;


    ABecLaplacian* mac_op;

      // localMF defined here MARC DAY JULIA REAM
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
    // Internal parameters for timestepping.
    //
    static Real init_shrink;   // reduction factor of first esimated timestep
    static Real cfl;           // desired maximum cfl
    static int  curv_stencil_height; // should be 4

    static int  MOF_TURN_OFF_LS;
    static int  MOF_DEBUG_RECON;
    static int  MOFITERMAX;
// default: tessellating fluid => default==1
//          non-tesselating or tesselating solid => default==0
    static Vector<int> truncate_volume_fractions;  
     // =0 do not hybridize with particles (default)
     // =1 hybridize with particles
    static int particles_flag;  
    static Real particles_weight_XD;  //=0.0 particles have no effect  
    static Real particles_weight_VEL; //=0.0 particles have no effect  
    static int particle_nsubdivide;  
    static int particle_max_per_nsubdivide;  
    static int particle_min_per_nsubdivide;  
    static Real particle_volume;
    static Real particle_relaxation_time_to_fluid;
    static Real fluid_relaxation_time_to_particle;
    static int particle_interaction_ngrow;

    static Real truncate_thickness;  // default is 2.0

    // 0  low order space and time
    // 1  SEM space and time
    // 2  SEM space 
    // 3  SEM time
    static int  enable_spectral;
    static int  viscous_enable_spectral; // def 0
    static int  projection_enable_spectral;
    static int  SEM_upwind;
     //0=div(uS)-S div(u)    1=u dot grad S  def=0
    static int  SEM_advection_algorithm;
     // continuous_mof = 0 => regular MOF
     // continuous_mof = 2 => CMOF ("continuous MOF ")
     // continuous_mof = 3 => CLSVOF except use MOF in cells with >2 materials.
     // continuous_mof = 4 => CLSVOF except use CMOF in cells w/ >2 materials.
     // continuous_mof = 5 => use CLSVOF everywhere (not recommended)
    static int  continuous_mof; // default is 0
     // =1 EXT_DIR,REFLECT_EVEN,embedded; =2 same as 1 plus triple points.
     // default is 1
    static int  force_cmof_at_triple_junctions; 
    static int  partial_cmof_stencil_at_walls;  // default is 0
    static Real change_max;    // maximum change in dt over a timestep
    static Real change_max_init; // maximum change in dt over a timestep

    static int n_scales;
    static Vector<Real> current_dt_group; //0,...,n_scales-1
    static Vector<Real> hold_dt_factors;  //0,...,n_scales-1
    static Vector<Real> fixed_dt_scales;  //0,...,n_scales-1
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
    static Real       MUSHY_THICK; // default==2.0
    static Real       gravity;
    static Real       gravity_reference_wavelen;
    static int        invert_gravity;
    static int        gravity_dir;
    static int        NUM_SCALARS;      // number of non-velocity components
    static Real       visc_coef; 

    static int gmres_precond_iter_base; // default=4
    static int include_viscous_heating; // default=0
     //smooth_type: 0=GSRB 1=ICRB 2=ILU  3=Jacobi
    static int smooth_type; // Lp.smooth_type def=0 GSRB
    static int bottom_smooth_type; // Lp.bottom_smooth_type def=2 ILU
    static int global_presmooth; // default=2
    static int global_postsmooth; // default=2
    static int use_mg_precond_in_mglib; //def=1 Lp.use_mg_precon_in_mglib
     //def=0.001 Lp.bottom_bottom_tol_factor
    static Real bottom_bottom_tol_factor;

    static Vector<int>  FSI_touch_flag;
    static int         elements_generated;
    static int         invert_solid_levelset;
     // 0=> u=u_solid if phi_solid>=0
     // 1=> u=u_solid_ghost if phi_solid>=0
     // 2=> GNBC (u=u_solid_ghost if phi_solid>=0)
    static int         law_of_the_wall;  // def=0
    static Real        wall_slip_weight; // def=0.0
     // 1..7 "1" corresponds to the GNBC option.
    static int         ZEYU_DCA_SELECT;  // DCA=dynamic contact angle
     // 0=fluid, tessellating 
     // 1=prescribed motion, non-tessellating, PROB.F90 solid 
     // 2=prescribed motion, non-tessellating, sci_clsvof.F90 solid 
     // 3=rigid FSI interior, tessellating, ice 
     //    (expansion/triple point dynamics allowed)
     // 4=link w/Kourosh, FSI, sci_clsvof.F90, fluids tessellate the
     //   computational domain, and the CTML elastic solid is embedded
     //   within the domain (non-tessellating) as if it was a rigid solid.
     // 5=FSI, tessellating, rigid (PROB.F90)
     // 6=FSI ice, represented as a tessellating fluid
     //   (initial geometry: sci_clsvof.F90)
     // 7=fluid, tessellating (initial geometry: sci_clsvof.F90)
     // 8=FSI, non-tessellating link w/user defined fortran code.
    static Vector<int>  FSI_flag; 

    static int         CTML_FSI_numsolids;
    static int         CTML_force_model; // 0=Lag force  1=Lag stress
    static int         CTML_FSI_init;

    static Vector<int>  shock_timestep;

    static int State_Type;
    static int Umac_Type;
    static int Vmac_Type;
    static int Wmac_Type;
    static int LS_Type;
    static int DIV_Type;
    static int Solid_State_Type;
    static int Tensor_Type;
    static int TensorXU_Type;
    static int TensorYU_Type;
    static int TensorZU_Type;
    static int TensorZV_Type;
    static int XDmac_Type;
    static int YDmac_Type;
    static int ZDmac_Type;
    static int NUM_STATE_TYPE;

     // def=1
     // 0=> I(u_GL) dot grad u_GG, no sync project
     // 1=> grad dot (u_GL u_GG) - u_GG div u_GL, no sync project
     // 2=> grad dot (u_GL u_GG) - u_GG div u_GL, sync project
    static int conservative_div_uu; 
    static Vector<Real> tension_min;
    static Vector<Real> tension_T0;
    static Vector<Real> tension_slope;
    static Vector<Real> tension; // 1-2 1-3 1-4 2-3 2-4 3-4
    static Vector<Real> prefreeze_tension; // 1-2 1-3 1-4 2-3 2-4 3-4
    static Vector<Real> recalesce_model_parameters; // 3*nmat

     // (dir,side)  (1,1),(2,1),(3,1),(1,2),(2,2),(3,2)
    static Vector<Real> outflow_velocity_buffer_size;

    static Vector<Real> cap_wave_speed; // nten components

    static Real real_number_of_cells;

     // 1.0/(den_max * mglib_min_coeff_factor) default=1000.0
    static Real mglib_min_coeff_factor; 

    static int custom_nucleation_model;

    static int FD_curv_interp; // def==1
    static int vof_height_function; // def==0

     // index: 1..nmat
    static Vector<Real> cavitation_pressure; 
     // index: 1..nmat
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
    static Real nucleation_init_time; //def 0.0
    static int pos_sites_random_flag; //default=1
    static int n_sites;
    static Vector<Real> pos_sites;

     // if microlayer_size>0 and LS_substrate>0 in 5x5x5 stencil,
     // then the microlayer model is activated.
     // subroutine get_vel_phasechange returns the larger of the microlayer
     // model rate and the conventional rate. 
    static Vector<Real> microlayer_size; // minimum transition length
    static Vector<Real> macrolayer_size; // maximum transition length

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

    static Real R_Palmore_Desjardins;

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
     // 4=Tanasawa model (evaporation and condensation for saturated material)
     // 5=evaporation
     // 1-2 1-3 1-4 2-3 2-4 3-4
     // 2-1 3-1 4-1 3-2 4-2 4-3
    static Vector<int> freezing_model;
     // 0=Palmore and Desjardins (same as Villegas et al)
     // 1=Tanasawa  2=Schrage  3=Kassemi
     // 4=Stefan model in which T_interface=f(Pressure_smooth)
    static Vector<int> Tanasawa_or_Schrage_or_Kassemi; 
     // index: 1..2*nten
     // value: 1..num_species_var
    static Vector<int> mass_fraction_id; 
     // index: 1..num_species_var
     // value: 1..nmat
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
    static Vector<int> distribute_mdot_evenly; // 1..2*nten
     // 0 - sum mdot <>0   
     // 1 - distribute sum -mdot to the source
     // -1 - distribute sum -mdot to the dest
    static Vector<int> constant_volume_mdot; // 1..2*nten
    static Vector<int> constant_density_all_time; // 1..nmat, def=1
    static int is_phasechange;
    static int normal_probe_size;  // should always be 1
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

    static int ngrowFSI;//should be 3
     //should be 12=velocity+LS+temperature+flag+stress (3D)
    static int nFSI_sub;
    static Vector<int> im_solid_map; //nparts components, in range 0..nmat-1
    static Vector<int> im_elastic_map;

    static int ngrow_expansion;  // should always be 2
    static int prescribe_temperature_outflow;

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

    static int visual_tessellate_vfrac; // def=0
    static int visual_revolve;
    static int visual_phase_change_plot_int; //def=0
    static int visual_buoyancy_plot_int; //def=0
    static int visual_divergence_plot_int; //def=0

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

    static int num_species_var; 
    static int num_materials;
    static int use_supermesh; //def=0
    static int ncomp_sum_int_user1; //def=0
    static int ncomp_sum_int_user2; //def=0

     // this is set according to store_elastic_data
    static int num_materials_viscoelastic;

     //num_material_elastic*NUM_TENSOR_TYPE
    static int NUM_CELL_ELASTIC; 

     // den,denA,E,T,P,species,configuration tensor
    static int num_state_material;
    static int num_state_base;
    static int ngeom_raw;
    static int ngeom_recon;
    static int ncghost_elastic;
    static int ncghost_state;

    static int ignore_div_up; // def=0
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

    static Vector<Real> elastic_time; // viscoelastic parameter

     // 0,1 => viscoelastic FENE-CR material  
     //     updating Q:
     //       (i) CISL advection Q^advect=Q^n(x-V dt)
     //       (ii) A^advect=Q^advect+I
     //       (iii) X=I+dt gradu
     //       (iv)  Q^{n+1}=X A^adv X^T -I=
     //           Q^advect+dt (2D) + dt gradu Q + dt Q gradu^T+O(dt^2)
     //       (v) Q^n+1=lambda Q^n+1/(lambda+dt)
     //   
     // 2=> elastic material  Q=bulk_modulus*(grad X + grad X^{T}) (linear ex.)
     // 3=> incremental elastic model 
     //   (a) cell centered formulation
     //   (b) face centered formulation
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
     //
    static Vector<int> viscoelastic_model; // def=0
    static Vector<int> viscosity_state_model; // def=0
    static Vector<int> les_model;
     // def=0 if shock_timestep=1, 1 otherwise
     // applies to both density and temperature.
    static Vector<int> temperature_primitive_variable; 

    static Vector<int> store_elastic_data; // viscoelastic parameter
    static Vector<Real> elastic_viscosity; // viscoelastic parameter
    static Vector<Real> elastic_regularization; // elastic parameter
    static Vector<Real> shear_modulus; // elastic parameter
    static Vector<Real> lame_coefficient; // elastic parameter
    static Vector<int> linear_elastic_model; // elastic parameter
    static Vector<Real> damping_coefficient; // dv/dt = -mu v
    static Vector<Real> polymer_factor;

    static int adv_dir;
    static Real adv_vel;
    static int axis_dir;

     // 0 - centroid furthest from uncaptured centroid
     // 1 - use MOF error
    static int mof_error_ordering;
    static Vector<int> mof_ordering;

    static Real rgasinlet;
    static Real slipcoeff;
    static Real vinletgas;
    static Real twall;
    static Real walltemp;
     // 1=EILE (default), -1=Weymouth Yue
     // 2=always EI   3=always LE
    static int EILE_flag; 
    static int bicgstab_max_num_outer_iter;
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
    static Real denconst_max;
    static Real denconst_min;
     // if stokes_flow==1, do not include u dot grad u term.
    static int stokes_flow;  // def=0
     // if cancel_advection==1, do not include u dot grad S terms for any
     // S corresponding to velocity, interface, or state variables.
    static int cancel_advection;  // def=0
    static Vector<Real> denconst_interface;
    static Vector<Real> added_weight;

    static Vector<Real> stiffPINF;
    static Vector<Real> prerecalesce_stiffCP; // def=4.1855E+7
    static Vector<Real> prerecalesce_stiffCV; // def=4.1855E+7
    static Vector<Real> stiffCP; // def=4.1855E+7
    static Vector<Real> stiffCV; // def=4.1855E+7
    static Vector<Real> stiffGAMMA;

    static int constant_viscosity;

    static Real angular_velocity;
    static Vector<Real> DrhoDT; // default 0.0

     //1=>rho=rho(T)
     //2=>Du/Dt=-grad (p-rho0 g dot z)/rho0 - g DrhoDT (T-T0) 
    static Vector<int> override_density; // default 0 
    static Vector<Real> viscconst;
    static Real viscconst_max;
    static Real viscconst_min;
    static Vector<Real> viscconst_eddy;
    static Vector<Real> prerecalesce_viscconst;
    static Vector<Real> viscconst_interface;
    static Real smoothing_length_scale;
    static Vector<Real> heatviscconst;
    static Real heatviscconst_max;
    static Real heatviscconst_min;
    static Vector<Real> prerecalesce_heatviscconst;
    static Vector<Real> heatviscconst_interface;
      // num_species * num_materials size.
    static Vector<Real> speciesconst;
      // num_species * num_materials size.
    static Vector<Real> speciesviscconst;
      // num_species * nten
    static Vector<Real> speciesviscconst_interface;
      // num_species
    static Vector<Real> species_molar_mass;  // def=1
     // default is solidheat_flag=0
    static int solidheat_flag;  // 0=diffuse in solid 1=dirichlet 2=neumann
    static Vector<int> material_type; // 0 incomp   
     //0 incomp; material_type_evap needed for the Kassemi model.
    static Vector<int> material_type_evap; 
    static Vector<int> material_type_lowmach; // 0 incomp   
    static Vector<int> material_type_visual; // 0 incomp   
    static Real wait_time;
    static Real advbot;
    static Real inflow_pressure;
    static Real outflow_pressure;
    static Real period_time;

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

    static int always_use_bicgstab;
    static int  viscous_maxiter;
    static int  multilevel_maxcycle;
    static int  initial_project_cycles;
    static int  initial_viscosity_cycles;
    static int  initial_thermal_cycles;
     // 0 - MGPCG  1-PCG 2-MINV=I
    static int  project_solver_type;
    static int  initial_cg_cycles;
    static int  num_divu_outer_sweeps;
    static int  debug_dot_product;

    static int curv_index;
    static int pforce_index;
    static int faceden_index;
    static int facecut_index;
    static int icefacecut_index;
    static int icemask_index;
    static int facevisc_index;
    static int faceheat_index;
    static int facevel_index;
    static int facespecies_index;
    static int smoothing_index;
    static int massface_index;
    static int vofface_index;
    static int ncphys;

};

}/* namespace amrex */

#endif /*_NavierStokes_H_*/
