#include <algorithm>
#include <vector>

#include <cstdio>
#include <cmath>

#include <string>
#include <stdlib.h>

#include <AMReX_CoordSys.H>
#include <AMReX_Geometry.H>
#include <AMReX_BoxDomain.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_Utility.H>
#include <AMReX_TagBox.H>

#include <AMReX_PlotFileUtil.H>
#ifdef AMREX_USE_HDF5
#include <AMReX_PlotFileUtilHDF5.H>
#endif

#include <NavierStokes.H>
#include <INTEGRATED_QUANTITY.H>
#include <DRAG_COMP.H>
#include <INTERP_F.H>
#include <MACOPERATOR_F.H>
#include <MARCHING_TETRA_F.H>
#include <NAVIERSTOKES_F.H>
#include <GLOBALUTIL_F.H>
#include <TECPLOTUTIL_F.H>
#include <GODUNOV_F.H>
#include <MASS_TRANSFER_F.H>
#include <PROB_CPP_PARMS_F.H>
#include <PROB_F.H>
#include <MOF_F.H>
#include <PLIC_F.H>
#include <LEVEL_F.H>
#include <MOF_REDIST_F.H>
#include <SOLIDFLUID_F.H>
#include <DERIVE_F.H>
#include <MG_F.H>
#include <INDEX_TYPE_MACROS.H>

#ifdef MVAHABFSI
#include <CTMLFSI_F.H>
#endif

namespace amrex{

#define bogus_value 1.e20
#define show_norm2_flag 0
#define mf_check_inf_bounds 1

#define DEFAULT_MOFITERMAX 15
#define DEFAULT_MOFITERMAX_AFTER_PREDICT 15

#define PCOPY_DEBUG 0
//
// Static objects.
//
BCRec NavierStokes::phys_bc;
BCRec NavierStokes::viscosity_phys_bc;
BCRec NavierStokes::temperature_phys_bc;
BCRec NavierStokes::species_phys_bc;

int  NavierStokes::NS_geometry_coord=-1;

int  NavierStokes::profile_debug=0;
bool NavierStokes::ns_tiling=false;

int NavierStokes::POLYGON_LIST_MAX=1000;

int  NavierStokes::ns_time_order=1; // time_blocking_factor
int  NavierStokes::slab_step=0;
int  NavierStokes::SDC_outer_sweeps=0;
int  NavierStokes::divu_outer_sweeps=0;
int  NavierStokes::very_last_sweep=0;
int  NavierStokes::num_divu_outer_sweeps=1;
Real NavierStokes::prev_time_slab=0.0;
Real NavierStokes::cur_time_slab=0.0;
Real NavierStokes::vel_time_slab=0.0;
Real NavierStokes::advect_time_slab=0.0;
Real NavierStokes::upper_slab_time=0.0;
Real NavierStokes::lower_slab_time=0.0;
Real NavierStokes::delta_slab_time=0.0;
Real NavierStokes::prescribed_vel_time_slab=0.0;
Real NavierStokes::dt_slab=1.0;
int NavierStokes::advect_iter=0;

int  NavierStokes::show_mem = 0;
int  NavierStokes::show_timings = 1;
int  NavierStokes::verbose      = 0;
int  NavierStokes::debug_fillpatch = 0;
int  NavierStokes::check_nan    = 0;
// 1=curv 2=error heat 3=both
int  NavierStokes::fab_verbose  = 0;
// ns.output_drop_distribution
int  NavierStokes::output_drop_distribution = 0;
int  NavierStokes::extend_pressure_into_solid = 0;
Real NavierStokes::cfl = 0.5;
int  NavierStokes::MOF_TURN_OFF_LS=0;
int  NavierStokes::MOF_DEBUG_RECON=0;
int  NavierStokes::MOFITERMAX=DEFAULT_MOFITERMAX;
int  NavierStokes::MOFITERMAX_AFTER_PREDICT=DEFAULT_MOFITERMAX_AFTER_PREDICT;

/*
 continuous_mof=STANDARD_MOF

  regular MOF  minimize E=||x_ij^ref-x_ij^derived||
  subject to the constraint that F_ij^ref=F_ij^derived

   x_ij^ref=reference centroid in cell ij
   x_ij^derived=derived centroid in cell ij for a given slope and
     intercept.
   F_ij^ref=reference volume fraction in cell ij
   F_ij^derived=derived volume fraction in cell ij for a given slope and
     intercept.   

continuous_mof>STANDARD_MOF

  CMOF  minimize E=||xS_ij^ref-xS_ij^derived||  "S"=super cell
  subject to the constraint that F_ij^ref=F_ij^derived

   xS_ij^ref=reference centroid in cell stencil i'=i-1,..,i+1,
     j'=j-1,..,j+1

   xS_ij^derived=derived centroid in cell stencil for a given slope and
     intercept. 
   F_ij^ref=reference volume fraction in cell
   F_ij^derived=derived volume fraction in cell for a given
     slope and intercept.



NOTE: rigid materials are not counted as materials in a cell.  Rigid 
materials are immersed into the fluid(s). 

*/

int  NavierStokes::continuous_mof=CMOF_X;
int  NavierStokes::update_centroid_after_recon=0;

//make MOFITERMAX_AFTER_PREDICT=0 if mof_decision_tree_learning>=100^d

#if (AMREX_SPACEDIM==2)
int  NavierStokes::mof_decision_tree_learning=10*10;//100x100 production runs
#elif (AMREX_SPACEDIM==3)
int  NavierStokes::mof_decision_tree_learning=10*10*10;//100^3 production runs
#else
something wrong
#endif
int  NavierStokes::mof_machine_learning=0;

//stability (or well posedness) theory for linear operator P: Pu=f
//P represents the linearized NavierStokes equations for two phase
//flows with surface tension (see e.g. Yaohong Wang's thesis).
//f represents a delta function source term for the centroid variable which
//mimics centroid initialization.
//the operator P is well posed if P^{-1} exists:
//u=P^{-1}f
//||u||<=||P^{-1}|| ||f||
//if P has a zero eigenvalue, i.e. a non-trivial null space, then
//P^{-1} cannot be defined, and there does not exist a "C" such that
//||u||<=C||f||.
//suppose Pu^{null}=0 and ||u^{null}||>0 then
//if Pu=f, then P(u+alpha u^{null})=f for any alpha.
//||u+alpha u^{null}||>=||alpha u^{null}||-||u||   (||A||<=||A-B||+||B||) 
//pick alpha=2C ||f||/||u^{null}|| which leads to contradiction.
//In practice, instability is manifested by lack of convergence, a
//difference between simulation and experiments, or in some cases,
//solution "blow up."
//CMOF advection fixes the problem because suppose x=x^{checkerboard}
//then MOF advection has x^{n+1}=x^{n} (assuming zero velocity), but
//CMOF advection has x^{n+1}=x^{n}/3 (lim_{n\rightarrow\infty}||x^{n}||=0)
//(centroid relative to
//to cell center).   i.e. CMOF advection automatically removes any checkerboard
//mode from appearing in x.  Albeit, CMOF advection leads to a less accurate
//method than MOF (when MOF converges), but CMOF, as we report here 
//anecdotedly, is shown to be more robust than MOF.
//Remark: if high frequency O(h) noise is added to the centroid, then if there
//is any surface tension and viscosity, then the noise should be dissipated
//by the method (i.e. "exact" solution is smooth).  But, MOF solutions exhibit
//persistant noisy data for many surface tension driven flows with viscosity,
//differing from what is expected from the PDE model and from experiments, 
//whereas CMOF solutions are consistent with both.

Vector<Real> NavierStokes::centroid_noise_factor; 

int  NavierStokes::partial_cmof_stencil_at_walls=1;

int  NavierStokes::enable_spectral=0;
// default: tessellating fluid => default==1
//          non-tessellating or tessellating solid => default==0
Vector<int> NavierStokes::truncate_volume_fractions; 

int NavierStokes::particle_nsubdivide=4; 
int NavierStokes::particle_max_per_nsubdivide=10; 

Real NavierStokes::truncate_thickness=2.0;  

Real NavierStokes::init_shrink  = 1.0;
Real NavierStokes::change_max   = 1.01;
Real NavierStokes::change_max_init = 1.01;

Vector<Real> NavierStokes::NS_coflow_Z; 
Vector<Real> NavierStokes::NS_coflow_R_of_Z; 

Vector<Real> NavierStokes::NS_DRAG_integrated_quantities; 
Vector<int> NavierStokes::NS_DRAG_integrated_quantities_sweep;

Vector<Real> NavierStokes::NS_sumdata; 
Vector<int> NavierStokes::NS_sumdata_type; 
Vector<int> NavierStokes::NS_sumdata_sweep; 

Real NavierStokes::fixed_dt     = 0.0;
Real NavierStokes::fixed_dt_init = 0.0;
Real NavierStokes::min_velocity_for_dt = 1.0e-12;
Real NavierStokes::min_stefan_velocity_for_dt = 1.0e-12;
Real NavierStokes::fixed_dt_velocity = 0.0;
Real NavierStokes::dt_max       = 1.0e+10;
Real NavierStokes::gravity      = 0.0;
// default=diagonal length of the domain tangent to the
// gravity direction.
Real NavierStokes::gravity_reference_wavelen = 0.0;

int NavierStokes::gravity_dir = AMREX_SPACEDIM;
int NavierStokes::invert_gravity = 0;
int NavierStokes::incremental_gravity_flag = 0;
int NavierStokes::segregated_gravity_flag = 0;

Vector<Real> NavierStokes::gravity_vector;

int  NavierStokes::sum_interval = -1;
int  NavierStokes::NUM_SCALARS  = 0;

// these vars used to be in probin
Real NavierStokes::denfact=1.0;
Real NavierStokes::velfact=0.0;

Real NavierStokes::xblob=0.0;
Real NavierStokes::yblob=0.0;
Real NavierStokes::zblob=0.0;
Real NavierStokes::radblob=1.0;

Real NavierStokes::xblob2=0.0;
Real NavierStokes::yblob2=0.0;
Real NavierStokes::zblob2=0.0;
Real NavierStokes::radblob2=0.0;

Real NavierStokes::xblob3=0.0;
Real NavierStokes::yblob3=0.0;
Real NavierStokes::zblob3=0.0;
Real NavierStokes::radblob3=0.0;

Real NavierStokes::xblob4=0.0;
Real NavierStokes::yblob4=0.0;
Real NavierStokes::zblob4=0.0;
Real NavierStokes::radblob4=0.0;

Real NavierStokes::xblob5=0.0;
Real NavierStokes::yblob5=0.0;
Real NavierStokes::zblob5=0.0;
Real NavierStokes::radblob5=0.0;

Real NavierStokes::xblob6=0.0;
Real NavierStokes::yblob6=0.0;
Real NavierStokes::zblob6=0.0;
Real NavierStokes::radblob6=0.0;

Real NavierStokes::xblob7=0.0;
Real NavierStokes::yblob7=0.0;
Real NavierStokes::zblob7=0.0;
Real NavierStokes::radblob7=0.0;

Real NavierStokes::xblob8=0.0;
Real NavierStokes::yblob8=0.0;
Real NavierStokes::zblob8=0.0;
Real NavierStokes::radblob8=0.0;

Real NavierStokes::xblob9=0.0;
Real NavierStokes::yblob9=0.0;
Real NavierStokes::zblob9=0.0;
Real NavierStokes::radblob9=0.0;

Real NavierStokes::xblob10=0.0;
Real NavierStokes::yblob10=0.0;
Real NavierStokes::zblob10=0.0;
Real NavierStokes::radblob10=0.0;

// force tag=0 outside this specified box if ractivex>0
Real NavierStokes::xactive=0.0;
Real NavierStokes::yactive=0.0;
Real NavierStokes::zactive=0.0;
Real NavierStokes::ractive=0.0;
Real NavierStokes::ractivex=0.0;
Real NavierStokes::ractivey=0.0;
Real NavierStokes::ractivez=0.0;

int  NavierStokes::probtype=0;
int  NavierStokes::adapt_quad_depth=1;

int  NavierStokes::step_through_data=0;   

// =0 (solids embed, fluids tessellate), 
// =1 (solids and fluids tessellate)
// =3 (solids and fluids tessellate, if F_solid>1/2, replace with F_solid=1,
//    if F_solid<1/2, replace with F_solid=0.
int  NavierStokes::visual_tessellate_vfrac=0;   
int  NavierStokes::visual_revolve=0;   
int  NavierStokes::visual_output_raw_State_Type=0; 
int  NavierStokes::visual_output_raw_mac_Type=0; 
int  NavierStokes::visual_phase_change_plot_int=0; 
int  NavierStokes::visual_buoyancy_plot_int=0; 
int  NavierStokes::visual_divergence_plot_int=0; 
//set visual_WALLVEL_plot_int>0 in order to generate WALLVEL*.plt files.
int  NavierStokes::visual_WALLVEL_plot_int=0; 
int  NavierStokes::visual_drag_plot_int=0; 
//default: tecplot nodes
//0=tecplot nodes
//1=plt file cells
//2=tecplot cells (piecewise constant reconstruction).
int  NavierStokes::visual_nddata_format=0;  

int NavierStokes::visual_compare=0; 
Vector<int> NavierStokes::visual_ncell;

// 0..sdim-1
int NavierStokes::slice_dir=0;
Vector<Real> NavierStokes::xslice;

Vector<Real> NavierStokes::outer_error_all_solver_calls;
Vector<int> NavierStokes::number_vcycles_all_solver_calls;
Vector<int> NavierStokes::max_lev0_cycles_all_solver_calls;
Vector<int> NavierStokes::median_lev0_cycles_all_solver_calls;
Vector<int> NavierStokes::lev0_cycles_list;
Vector<int> NavierStokes::number_solver_calls;

Vector<int> NavierStokes::type_ONES_flag;
int NavierStokes::color_ONES_count;
int NavierStokes::coarsest_ONES_level;
Vector<int> NavierStokes::singular_patch_flag;
Vector<Real> NavierStokes::ones_sum_global;

Vector< Vector<Real> > NavierStokes::min_face_wt;
Vector< Vector<Real> > NavierStokes::max_face_wt;

Vector< Vector<Real> > NavierStokes::delta_mass;

Real NavierStokes::max_problen=0.0;
Vector< Vector<Real> > NavierStokes::minLS;
Vector< Vector<Real> > NavierStokes::maxLS;

Vector<int> NavierStokes::map_forward_direct_split;
Vector<int> NavierStokes::normdir_direct_split;
int NavierStokes::dir_absolute_direct_split=0;
int NavierStokes::order_direct_split=0; // base_step mod 2

Vector<Real> NavierStokes::mdotplus;
Vector<Real> NavierStokes::mdotminus;
Vector<Real> NavierStokes::mdotcount;
Vector<Real> NavierStokes::mdot_sum;
Vector<Real> NavierStokes::mdot_sum2;
Vector<Real> NavierStokes::mdot_lost;

Vector<Real> NavierStokes::mdotplus_complement;
Vector<Real> NavierStokes::mdotminus_complement;
Vector<Real> NavierStokes::mdotcount_complement;
Vector<Real> NavierStokes::mdot_sum_complement;
Vector<Real> NavierStokes::mdot_sum2_complement;
Vector<Real> NavierStokes::mdot_lost_complement;

Vector<Real> NavierStokes::curv_min;
Vector<Real> NavierStokes::curv_max;


// force AMR tagging within these specified boxes.
int NavierStokes::nblocks=0;
Vector<Real> NavierStokes::xblocks;
Vector<Real> NavierStokes::yblocks;
Vector<Real> NavierStokes::zblocks;
Vector<Real> NavierStokes::rxblocks;
Vector<Real> NavierStokes::ryblocks;
Vector<Real> NavierStokes::rzblocks;

int NavierStokes::tecplot_max_level=-1;
int NavierStokes::max_level_for_use=0;

// tagflag forced to 0 OUTSIDE these specified boxes.
int NavierStokes::ncoarseblocks=0;
Vector<Real> NavierStokes::xcoarseblocks;
Vector<Real> NavierStokes::ycoarseblocks;
Vector<Real> NavierStokes::zcoarseblocks;
Vector<Real> NavierStokes::rxcoarseblocks;
Vector<Real> NavierStokes::rycoarseblocks;
Vector<Real> NavierStokes::rzcoarseblocks;

int  NavierStokes::override_bc_to_homogeneous=0;

int  NavierStokes::num_species_var=0;

// search num_materials,AMREX_SPACEDIM+1,SDIM+1,
//  dcomp,pcomp,tcomp,scomp,
//  get_mm_scomp_solver,dencomp,scomp_tensor,
//  velcomp,flagcomp
// nfacefrac,ncellfrac,mm_,cellmm,facemm
int  NavierStokes::num_materials=0;
int  NavierStokes::num_interfaces=0;

int  NavierStokes::ncomp_sum_int_user1=0;
int  NavierStokes::ncomp_sum_int_user2=0;
int  NavierStokes::ncomp_sum_int_user12=0;

// set using elastic_viscosity, and other criteria
int  NavierStokes::num_materials_viscoelastic=0;

int  NavierStokes::num_state_material=ENUM_SPECIESVAR; // den,T
int  NavierStokes::num_state_base=ENUM_SPECIESVAR; // den,T
int  NavierStokes::ngeom_raw=AMREX_SPACEDIM+1;
int  NavierStokes::ngeom_recon=ENUM_NUM_MOF_VAR;

// vel, pres, num_state_material x num_materials, ngeom_raw x num_materials, error ind
int  NavierStokes::State_Type=0;
// mac vel
int  NavierStokes::Umac_Type=State_Type+1;
int  NavierStokes::Vmac_Type=Umac_Type+1;
int  NavierStokes::Wmac_Type=Vmac_Type+AMREX_SPACEDIM-2;
// LS 1..num_materials, LS_slope sdim x num_materials
int  NavierStokes::LS_Type=Wmac_Type+1;
// -(pnew-pold)/(rho c^2 dt) + dt mdot/vol
int  NavierStokes::DIV_Type=LS_Type+1;
int  NavierStokes::Solid_State_Type=DIV_Type+1;
int  NavierStokes::Tensor_Type=Solid_State_Type+1;
int  NavierStokes::TensorX_Type=Tensor_Type+1;
int  NavierStokes::TensorY_Type=TensorX_Type+1;
int  NavierStokes::TensorZ_Type=TensorY_Type+1;

int  NavierStokes::NUM_STATE_TYPE=TensorZ_Type+1;


Vector<Real> NavierStokes::compressible_dt_factor; 

int  NavierStokes::disable_advection=0;
int  NavierStokes::disable_pressure_solve=0;

Vector<Real> NavierStokes::elastic_time; // def=0

// if viscosity_state_model>=1, density and temperature come from State_Type
//  otherwise density=fort_denconst, temperature=fort_tempconst
// mu=get_user_viscconst(im,density,temperature)
// in: PROB.F90,
// REAL_T function get_user_viscconst(im,density,temperature)
// MITSUHIRO: viscosity_state_model=2
Vector<int> NavierStokes::viscosity_state_model; // def=0
// 0 => viscoelastic FENE-CR material  
//     updating Q:
//       (i) CISL advection Q^advect=Q^n(x-V dt)
//       (ii) A^advect=Q^advect+I
//       (iii) X=I+dt gradu
//       (iv)  Q^{n+1}=X A^adv X^T -I=
//           Q^advect+dt (2D) + dt gradu Q + dt Q gradu^T+O(dt^2)
//       (v) lambda'=lambda * (1-tr(A)/L^2)
//       (vi) Q^n+1=lambda' Q^n+1/(lambda'+dt)
// 1 => Oldroyd-B
//       (i) lambda'=lambda
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
// 5=> FENE-P 
//       (v) lambda'=lambda * (1-tr(A)/L^2)
//       (vi) Q_t = -(1/lambda')(Q+I * tr(A)/L^{2})  
//       (vii) Q^{n+1}-Q^star=-(dt/lambda')(Q^{n+1} + I * tr(A)/L^2)
//       (viii) (1+dt/lambda')Q^{n+1}=Q^star - dt * I * tr(A)/L^2
//       (viv) Q^{n+1}=(lambda'/(lambda'+dt))*(Q^star-dt*I*tr(A)/L^2)
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
Vector<int> NavierStokes::viscoelastic_model; // def=0
Vector<int> NavierStokes::les_model; // def=0

int NavierStokes::transposegradu=0;

Vector<int> NavierStokes::store_elastic_data; // def=0, 0...num_materials-1
Vector<Real> NavierStokes::elastic_viscosity; // def=0
Vector<Real> NavierStokes::static_damping_coefficient; // def=0

Vector<Real> NavierStokes::Carreau_alpha; // def=1
Vector<Real> NavierStokes::Carreau_beta; // def=0
Vector<Real> NavierStokes::Carreau_n; // def=1
Vector<Real> NavierStokes::Carreau_mu_inf; // def=0
Vector<int> NavierStokes::shear_thinning_fluid; // def=0

Vector<Real> NavierStokes::concentration; // def=0
Vector<Real> NavierStokes::etaL; // def=0 (etaL0)
Vector<Real> NavierStokes::etaS; // def=0
Vector<Real> NavierStokes::etaP; // def=0 (etaP0)

// (1/L)   eps=0.01
// viscoelastic_model=0  FENE CR   trac(A)<L^2 lambda(A)>eps/L^2
// viscoelastic_model=1  Oldroyd B trac(A)<inf lambda(A)>eps/L^2
// viscoelastic_model=5  FENE P    trac(A)<L^2 lambda(A)>eps/L^2
// viscoelastic_model=6  linearPTT trac(A)<inf lambda(A)>eps/L^2
Vector<Real> NavierStokes::polymer_factor; // def=0

 // 0 - centroid furthest from uncaptured centroid
 // 1 - smallest MOF error
int  NavierStokes::mof_error_ordering=0;

Vector<int> NavierStokes::mof_ordering; // def=0

// adv_dir=1,..,sdim+1
int  NavierStokes::adv_dir=0;
Real NavierStokes::adv_vel=1.0;
int  NavierStokes::axis_dir=0;
Real NavierStokes::rgasinlet=0.0;
Real NavierStokes::slipcoeff=0.0;
Real NavierStokes::vinletgas=0.0;
Real NavierStokes::twall=0.1;

// 1=EILE (default), -1=Weymouth Yue
// 2=always EI   3=always LE
int NavierStokes::EILE_flag=1;

int  NavierStokes::krylov_subspace_max_num_outer_iter=60;
Real NavierStokes::projection_pressure_scale=1.0;
Real NavierStokes::projection_velocity_scale=1.0;

int NavierStokes::ngrow_distance=4;
int NavierStokes::ngrow_make_distance=3;

int NavierStokes::num_elements_blobclass=-32767;
int NavierStokes::BLB_MATRIX=-32767;
int NavierStokes::BLB_RHS=-32767;
int NavierStokes::BLB_VEL=-32767;
int NavierStokes::BLB_INT_MOM=-32767;
int NavierStokes::BLB_ENERGY=-32767;
int NavierStokes::BLB_MASS_VEL=-32767;
int NavierStokes::BLB_VOL=-32767;
int NavierStokes::BLB_CEN_INT=-32767;
int NavierStokes::BLB_CEN_ACT=-32767;
int NavierStokes::BLB_PERIM=-32767;
int NavierStokes::BLB_PERIM_MAT=-32767;
int NavierStokes::BLB_TRIPLE_PERIM=-32767;
int NavierStokes::BLB_CELL_CNT=-32767;
int NavierStokes::BLB_CELLVOL_CNT=-32767;
int NavierStokes::BLB_MASS=-32767;
int NavierStokes::BLB_PRES=-32767;

Vector<int> NavierStokes::im_solid_map; //nparts components, in range 0..num_materials-1
// 0<=im_elastic_map<num_materials
Vector<int> NavierStokes::im_elastic_map; //0...num_materials_viscoelastic-1

Real NavierStokes::real_number_of_cells=0.0; 

Real NavierStokes::mglib_max_ratio=1.0e+5; 
Real NavierStokes::min_interior_coeff=0.0; 

int NavierStokes::idx_umac_material_mf=-1;
int NavierStokes::idx_umac_mask_material_mf=-1;
int NavierStokes::idx_scalar_mask_material_mf=-1;

int NavierStokes::hydrate_flag=0; 
int NavierStokes::post_init_pressure_solve=1; 

Vector<Real> NavierStokes::tension_slope;
Vector<Real> NavierStokes::tension_min;
Vector<Real> NavierStokes::tension_T0;
Vector<Real> NavierStokes::tension;
Vector<Real> NavierStokes::tension_init;

Real NavierStokes::unscaled_min_curvature_radius=2.0;
Vector<Real> NavierStokes::prefreeze_tension;

Vector<Real> NavierStokes::outflow_velocity_buffer_size;

Vector<Real> NavierStokes::cap_wave_speed;

Vector<Real> NavierStokes::grid_stretching_parameter;

Vector<Real> NavierStokes::hardwire_Y_gamma;
Vector<Real> NavierStokes::hardwire_T_gamma;
// "p_boil" in Dodd and Ferrante
// default: 1.0e+6
Vector<Real> NavierStokes::reference_pressure;
Vector<Real> NavierStokes::accommodation_coefficient;
Vector<Real> NavierStokes::saturation_temp;
Vector<Real> NavierStokes::saturation_temp_curv;
Vector<Real> NavierStokes::saturation_temp_vel;
Vector<Real> NavierStokes::saturation_temp_min; //aka T_I_min
Vector<Real> NavierStokes::saturation_temp_max; //aka T_I_max

Vector<int> NavierStokes::microlayer_substrate;
Vector<Real> NavierStokes::microlayer_angle;
Vector<Real> NavierStokes::microlayer_size;
Vector<Real> NavierStokes::macrolayer_size;
Vector<Real> NavierStokes::max_contact_line_size;

// minimum transition length for heat conduction.
// e.g. k_heater=k * max(1,dx/thermal_microlayer_size)
Vector<Real> NavierStokes::thermal_microlayer_size;
// minimum transition length for viscosity. 
// e.g. visc+viscconst_eddy<visc * max(1,dx/shear_microlayer_size)
// e.g. (u_in - u_law_of_wall)/dx<u_in * max(1,dx/shear_microlayer_size)
Vector<Real> NavierStokes::shear_microlayer_size;
// e.g. k_model<k * dx/buoyancy_microlayer_size
Vector<Real> NavierStokes::buoyancy_microlayer_size;
// grad T \approx (T_probe-T_I)/phasechange_microlayer_size if material
// "im" macro scale probe cannot be found.
Vector<Real> NavierStokes::phasechange_microlayer_size;

// if freezing_model==0: (ambient air is 100 percent saturated)
//
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

Vector<Real> NavierStokes::microlayer_temperature_substrate; 

int NavierStokes::custom_nucleation_model=0;

int NavierStokes::FD_curv_interp=1;
int NavierStokes::vof_height_function=1;

Vector<Real> NavierStokes::cavitation_pressure;
Vector<Real> NavierStokes::cavitation_vapor_density;
Vector<Real> NavierStokes::cavitation_tension;

Vector<Real> NavierStokes::nucleation_pressure;
Vector<Real> NavierStokes::nucleation_pmg;
Vector<Real> NavierStokes::nucleation_mach;
Vector<Real> NavierStokes::nucleation_temp;
Real NavierStokes::nucleation_period=0.0;
Real NavierStokes::nucleation_init_time=0.0;
int NavierStokes::n_sites=0;
int NavierStokes::pos_sites_random_flag=1;
Vector<Real> NavierStokes::pos_sites;

int NavierStokes::perturbation_on_restart=0;
 // sin(2 pi k x /L)
int NavierStokes::perturbation_mode=0; // number of wavelen each dir.
// a fraction of delta T or the radial velocity rmax omega
Real NavierStokes::perturbation_eps_temp=0.0;
Real NavierStokes::perturbation_eps_vel=0.0;

// latent_heat<0 if condensation or solidification
// latent_heat>0 if boiling or melting
Vector<Real> NavierStokes::latent_heat;
Vector<Real> NavierStokes::latent_heat_slope;
Vector<Real> NavierStokes::latent_heat_min;
Vector<Real> NavierStokes::latent_heat_T0;

//ergs/(mol kelvin)
Real NavierStokes::R_Palmore_Desjardins=8.31446261815324e+7;

Vector<Real> NavierStokes::reaction_rate;
// 0=T_interface=TSAT-epsC K -epsV V ambient air is 100 percent
//   saturated.
//  
// 1=source term model (single equation for T with source term).
//   interpolation does not assume T=TSAT at the interface.
// 2=hydrate model 
//    a) V=V(P,T,RHO,C)
//    b) single equation for T with source term
//    c) single equation for C with source term.
// 3=wildfire combustion
//   1=liquid  2=ambient gas  3=solid wall  num_species_var=1
//   ->  12 13 23 21 31 32
//   latent_heat
//   ->  +L 0  0  -L 0  0
//   freezing_model
//   ->  6  0  0  6  0  0
//   mass_fraction_id
//   ->  1  0  0  1  0  0
//   Tanasawa_or_Schrage_or_Kassemi
//   ->  2  0  0  2  0  0
//   speciesviscconst
//   ->  0.0 D 0.0
//   distribute_from_target (distribution of the expansion term)
//   ->  0 0 0 1 0 0
//   molar_mass
//   ->  mL  m_ambient 0.0
//   species_molar_mass
//   ->  m_vapor
//   material_type
//   0 ?? 999
//
// 5=evaporation/condensation (Stefan model speed)
// 6=evaporation/condensation (Palmore and Desjardins, JCP 2019)
// 7=cavitation
Vector<int> NavierStokes::freezing_model;

//0=Palmore and Desjardins (Villegas, Tanguy, Desjardins) 
//1=Tanasawa  2=Schrage 3=Kassemi
Vector<int> NavierStokes::Tanasawa_or_Schrage_or_Kassemi; 

//ispec=rigid_fraction_id[0..num_materials-1]=1..num_species_var
Vector<int> NavierStokes::rigid_fraction_id;

//ispec=mass_fraction_id[0..2 num_interfaces-1]=1..num_species_var
Vector<int> NavierStokes::mass_fraction_id; 
//link diffused material to non-diff. (array 1..num_species_var)
//spec_material_id_LIQUID, spec_material_id_AMBIENT are 
//both input AND derived.
Vector<int> NavierStokes::spec_material_id_LIQUID; 
Vector<int> NavierStokes::spec_material_id_AMBIENT; 
// 0 - distribute to the destination material 
//     V=mdot/den_src
//     source term= dF * Vcell/dt^2 * (-1+ den_src/den_dst)
//
// 1 - distribute to the source material
//     V=mdot/den_dst
//     source term= dF * Vcell/dt^2 * (1-den_dst/den_src)
//
// source term is independent of "distribute_from_target" since the
// magnitude of the jump in velocity is always the same: mdot[1/rho].
//
// Palmore and Desjardins: u_Gamma=uL + mdot/rhoL=uG + mdot/rhoG
// if the source term is in the destination, then interface velocity is
//  u_source+mdot/rho_source  (velocity at interface is u_source)
// if the source term is in the source, then interface velocity is
//  u_dest+mdot/rho_dest (velocity at interface is u_dest)
//
// for boiling and evaporation, den_dst<<den_src => 
//   distribute_from_target=1 (pressure more accurate in the
//   liquid regions, also less reliance on dt being uniform, the
//   destination velocity has smaller magnitude.)
// for freezing and melting need to distribute the mdot into
// the liquid material:
//   distribute_from_target=0 melting 
//   distribute_from_target=1 freezing 

Vector<int> NavierStokes::distribute_from_target; // 1..2*num_interfaces

// 0 - mdot concentrated at full cells near interface on 1-side
// 1 - mdot distributed evenly to all full cells; cellvol weight
// 2 - mdot distributed evenly to all full cells; constant weight
Vector<int> NavierStokes::distribute_mdot_evenly; // 1..2*num_interfaces

//  For freezing in which liquid volume is fixed, the liquid density
//  must increase.  The mass fraction of "heavy liquid" increases.
//  For melting in which liquid volume is fixed, the liquid density
//  must decrease.  The mass fraction of "light liquid" increases.
// 0 - sum mdot <>0   
// 1 - distribute sum -mdot to the source
// -1 - distribute sum -mdot to the dest
Vector<int> NavierStokes::constant_volume_mdot; // 1..2*num_interfaces
Vector<int> NavierStokes::constant_density_all_time; // 1..num_materials, def=1

int NavierStokes::is_phasechange=0;
// 0=dirichlet at inflow
// 1=dirichlet at inflow and outflow
// 2=dirichlet at inflow and walls.
// 3=dirichlet at inflow, outflow, and walls.
int NavierStokes::prescribe_temperature_outflow=0; // default is 0

// 0=volume fraction  1=mass fraction 2=impedance fraction
int  NavierStokes::pressure_select_criterion=0;


Vector<Real> NavierStokes::vorterr;
Vector<Real> NavierStokes::pressure_error_cutoff;

// 0 (check mag, default) 1=check diff
int NavierStokes::pressure_error_flag=0;

Vector<Real> NavierStokes::temperature_error_cutoff;

Vector<Real> NavierStokes::tempcutoff;
Vector<Real> NavierStokes::tempcutoffmax;
Vector<Real> NavierStokes::tempconst;
Vector<Real> NavierStokes::initial_temperature;
Real NavierStokes::initial_temperature_diffuse_duration=0.0;

// default is zero which means "get_local_heat_source" returns a zero source.
Real NavierStokes::temperature_source=0.0;
Vector<Real> NavierStokes::temperature_source_cen;
Vector<Real> NavierStokes::temperature_source_rad;

Vector<Real> NavierStokes::density_floor;  // def=0.0
Vector<Real> NavierStokes::density_ceiling;  // def=1.0e+20
Vector<Real> NavierStokes::molar_mass;  // def=1
Vector<Real> NavierStokes::denconst;
Vector<Real> NavierStokes::denconst_interface_added;
int NavierStokes::stokes_flow=0;
int NavierStokes::cancel_advection=0;

Vector<Real> NavierStokes::stiffPINF;
Vector<Real> NavierStokes::prerecalesce_stiffCP;  // def=4.1855E+7
Vector<Real> NavierStokes::prerecalesce_stiffCV;  // def=4.1855E+7
Vector<Real> NavierStokes::stiffCP;  // def=4.1855E+7
Vector<Real> NavierStokes::stiffCV;  // def=4.1855E+7
Vector<Real> NavierStokes::stiffGAMMA; // def=1.4

// uncoupled_viscosity=0 => div (2 mu D)
// uncoupled_viscosity=1 => div (mu grad U)
int NavierStokes::uncoupled_viscosity=0;

Real NavierStokes::angular_velocity=0.0;
Real NavierStokes::centrifugal_force_factor=1.0;

//Du/Dt=-grad (p-rho0 g dot z)/rho0 - g DrhoDT (T-T0) 
//DrhoDT has units of 1/(Degrees Kelvin)
Vector<Real> NavierStokes::DrhoDT;  // def=0.0

// 1=>rho=rho(T,Y,z)
// 2=>Boussinesq approximation Du/Dt=-grad(p-rho0 g dot z)/rho0-g DrhoDT (T-T0)
// DrhoDT has units of 1/Temperature.
Vector<int> NavierStokes::override_density; // def=0
Vector<Real> NavierStokes::prerecalesce_viscconst;
Vector<Real> NavierStokes::viscconst;
Real NavierStokes::viscconst_max=0.0;
Real NavierStokes::viscconst_min=0.0;
Vector<Real> NavierStokes::viscconst_eddy_wall; //default = 0
Vector<Real> NavierStokes::viscconst_eddy_bulk; //default = 0
Vector<Real> NavierStokes::heatviscconst_eddy_wall; //default = 0
Vector<Real> NavierStokes::heatviscconst_eddy_bulk; //default = 0
Vector<Real> NavierStokes::speciesviscconst;// species mass diffusion coeff.
Vector<Real> NavierStokes::prerecalesce_heatviscconst;
Vector<Real> NavierStokes::heatflux_factor;
Vector<Real> NavierStokes::heatviscconst;
Real NavierStokes::heatviscconst_max=0.0;
Real NavierStokes::heatviscconst_min=0.0;
Vector<Real> NavierStokes::viscconst_interface;
Vector<Real> NavierStokes::heatviscconst_interface;
Vector<Real> NavierStokes::speciesconst;  
Vector<Real> NavierStokes::speciesreactionrate;  
Vector<Real> NavierStokes::speciesviscconst_interface;
// 1..num_species_var
Vector<Real> NavierStokes::species_molar_mass; // def=1
// 0=diffuse in solid 1=dirichlet 2=neumann
int NavierStokes::solidheat_flag=0; 

Vector<int> NavierStokes::material_type;
Vector<int> NavierStokes::material_type_interface;

//0 incomp; material_type_evap needed for the Kassemi model.
Vector<int> NavierStokes::material_type_evap;
Vector<int> NavierStokes::material_type_lowmach;
Vector<int> NavierStokes::material_type_visual;

Real NavierStokes::wait_time=0.0;
Real NavierStokes::advbot=1.0;
Real NavierStokes::inflow_pressure=0.0;
Real NavierStokes::outflow_pressure=0.0;
Real NavierStokes::period_time=0.0;

     // 0 - MGPCG  1-PCG  2-MINV=I
int  NavierStokes::project_solver_type=0;
// number of Jacobi method cycles elliptic solver initially does.
int  NavierStokes::initial_project_cycles=3;
// number of Jacobi method cycles elliptic solver initially does for the
// viscosity equation.
int  NavierStokes::initial_viscosity_cycles=3;
// number of Jacobi method cycles elliptic solver initially does for the
// temperature equation.
int  NavierStokes::initial_thermal_cycles=3;
// default is to do 5 MGPCG cycles, then restart the MGPCG iteration
// and do as many cycles as necessary in order to achieve convergence.
int  NavierStokes::initial_cg_cycles=5;
int  NavierStokes::debug_dot_product=0;

int NavierStokes::smooth_type = 2; // 0=GSRB 1=weighted Jacobi 2=ILU
int NavierStokes::bottom_smooth_type = 2; // 0=GSRB 1=weighted Jacobi 2=ILU
int NavierStokes::global_presmooth = 2;
int NavierStokes::global_postsmooth = 2;
int NavierStokes::use_mg_precond_in_mglib=1;
Real NavierStokes::bottom_bottom_tol_factor=0.001;

// for static contact line algorithm: 
//   Arienti and Sussman, IJMF
// for dynamic contact line algorithm using "get_use_DCA" model,
//   see work involving Yongsheng Lian and Sussman.
//   (velocity is input => dynamic angle is output)
// for dynamic contact line algorithm using GNBC, <unpublished>.
//   (dynamic angle is input => velocity is output)
// 0=> u=u_solid if phi_solid>=0
// 1=> u=u_solid_ghost if phi_solid>=0
// 2=> generalized Navier Boundary condition (GNBC),
//   for conventional contact line dynamics, 
//   modify "get_use_DCA" in GLOBALUTIL.F90.
Vector<int> NavierStokes::law_of_the_wall;
Vector<Real> NavierStokes::wall_model_velocity; //1..num_materials
Vector<int> NavierStokes::interface_mass_transfer_model; //1..2*num_interfaces
Real NavierStokes::wall_slip_weight=0.0;
int NavierStokes::ZEYU_DCA_SELECT=-1; // -1 = static angle

// FSI_FLUID=0
// FSI_PRESCRIBED_PROBF90=1
// FSI_PRESCRIBED_NODES=2
// FSI_ICE_PROBF90=3
// FSI_SHOELE_CTML=4
// FSI_RIGID_NOTPRESCRIBED=5
// FSI_ICE_NODES_INIT=6
// FSI_FLUID_NODES_INIT=7
// FSI_ICE_STATIC=9
Vector<int> NavierStokes::FSI_flag; 

Vector<int> NavierStokes::CTML_max_num_nodes_list;
int NavierStokes::CTML_max_num_elements_list=0;
int NavierStokes::CTML_FSI_num_scalars=0;

int NavierStokes::FSI_interval=1;
int NavierStokes::num_local_aux_grids=0;

Vector<int> NavierStokes::FSI_touch_flag; // 0..nthreads-1
// default: 1
Vector<int> NavierStokes::FSI_refine_factor; 
// default: 3
Vector<int> NavierStokes::FSI_bounding_box_ngrow; 

Vector<int> NavierStokes::ns_max_grid_size; 

int NavierStokes::CTML_FSI_numsolids = 0;

int NavierStokes::CTML_FSI_init = 0;

// 0=take into account sound speed only at t=0 if compressible.
// 1=always take into account sound speed
// 2=never take into account sound speed
Vector<int> NavierStokes::shock_timestep; 

Real NavierStokes::visc_coef=0.0;

int NavierStokes::include_viscous_heating=0;

int NavierStokes::multilevel_maxcycle=200;
int NavierStokes::multilevel_restart_period=20000;

// mg.bot_atol
// mg.visc_bot_atol
// mg.thermal_bot_atol
Real NavierStokes::minimum_relative_error = 1.0e-11;
Real NavierStokes::diffusion_minimum_relative_error = 1.0e-11;

Real NavierStokes::save_atol_b=1.0e-14;
Real NavierStokes::save_mac_abs_tol=1.0e-10;
Real NavierStokes::save_min_rel_error=1.0e-11;

Real NavierStokes::mac_abs_tol = 1.0e-10;
Real NavierStokes::visc_abs_tol = 1.0e-10;
Real NavierStokes::thermal_abs_tol = 1.0e-10;
Real NavierStokes::total_advance_time=0.0;

void extra_circle_parameters(
 Real& xblob2,Real& yblob2,Real& zblob2,Real& radblob2,
 Real& xblob3,Real& yblob3,Real& zblob3,Real& radblob3,
 Real& xblob4,Real& yblob4,Real& zblob4,Real& radblob4,
 Real& xblob5,Real& yblob5,Real& zblob5,Real& radblob5,
 Real& xblob6,Real& yblob6,Real& zblob6,Real& radblob6,
 Real& xblob7,Real& yblob7,Real& zblob7,Real& radblob7,
 Real& xblob8,Real& yblob8,Real& zblob8,Real& radblob8,
 Real& xblob9,Real& yblob9,Real& zblob9,Real& radblob9,
 Real& xblob10,Real& yblob10,Real& zblob10,Real& radblob10 ) {

 xblob2=0.0; 
 yblob2=0.0; 
 zblob2=0.0; 
 radblob2=0.0; 

 xblob3=0.0; 
 yblob3=0.0; 
 zblob3=0.0; 
 radblob3=0.0; 

 xblob4=0.0; 
 yblob4=0.0; 
 zblob4=0.0; 
 radblob4=0.0; 

 xblob5=0.0; 
 yblob5=0.0; 
 zblob5=0.0; 
 radblob5=0.0; 

 xblob6=0.0; 
 yblob6=0.0; 
 zblob6=0.0; 
 radblob6=0.0; 

 xblob7=0.0; 
 yblob7=0.0; 
 zblob7=0.0; 
 radblob7=0.0; 

 xblob8=0.0; 
 yblob8=0.0; 
 zblob8=0.0; 
 radblob8=0.0; 

 xblob9=0.0; 
 yblob9=0.0; 
 zblob9=0.0; 
 radblob9=0.0; 

 xblob10=0.0; 
 yblob10=0.0; 
 zblob10=0.0; 
 radblob10=0.0; 

} // subroutine extra_circle_parameters

void read_geometry_raw(int& geometry_coord,
	Vector<Real>& geometry_prob_lo,
        Vector<Real>& geometry_prob_hi,
	Vector<int>& geometry_is_periodic,
	int& geometry_is_any_periodic) {

 ParmParse pp("geometry");
 pp.get("coord_sys",geometry_coord);
 int geometry_coord_override=geometry_coord;
 pp.queryAdd("coord_sys_override",geometry_coord_override);
 geometry_coord=geometry_coord_override;

 pp.getarr("prob_lo",geometry_prob_lo,0,AMREX_SPACEDIM);
 pp.getarr("prob_hi",geometry_prob_hi,0,AMREX_SPACEDIM);
 geometry_is_periodic.resize(AMREX_SPACEDIM);
 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  geometry_is_periodic[dir]=0;
 pp.queryAdd("is_periodic",geometry_is_periodic,AMREX_SPACEDIM);
 geometry_is_any_periodic=0;
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  if (geometry_is_periodic[dir]==0) {
   // do nothing
  } else if (geometry_is_periodic[dir]==1) {
   geometry_is_any_periodic=1;
  } else
    amrex::Error("geometry_is_periodic[dir] invalid");
 } // dir=0...sdim-1

} // end subroutine read_geometry_raw

// this routine is called from main.cpp after
void fortran_deallocate_parameters() {

 fort_deallocate_module();

} // end subroutine fortran_deallocate_parameters

// this routine is called from main.cpp prior to:
//  1. AmrCore* amrptr = new AmrCore();
//  2. amrptr->init(strt_time,stop_time);
void fortran_parameters() {

 int geometry_coord;
 Vector<Real> geometry_prob_lo;
 Vector<Real> geometry_prob_hi;
 Vector<int> geometry_is_periodic;
 int geometry_is_any_periodic;
   //read_geometry_raw is declared in NavierStokes.cpp
 read_geometry_raw(geometry_coord,geometry_prob_lo,geometry_prob_hi,
		 geometry_is_periodic,geometry_is_any_periodic);

 Real problox=geometry_prob_lo[0];
 Real probloy=geometry_prob_lo[1];
 Real probloz=geometry_prob_lo[AMREX_SPACEDIM-1];
 Real probhix=geometry_prob_hi[0];
 Real probhiy=geometry_prob_hi[1];
 Real probhiz=geometry_prob_hi[AMREX_SPACEDIM-1];

 Real problen_min=probhix-problox;
 if (problen_min>0.0) {
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   Real problen=geometry_prob_hi[dir]-geometry_prob_lo[dir];
   if (problen>0.0) {
    if (problen<problen_min) 
     problen_min=problen;
   } else
    amrex::Error("problen invalid");
  }
 } else
  amrex::Error("problen_min invalid");

 Real microlayer_size_default=problen_min/1.0e+9;

 Real denfact;
 Real velfact=0.0;
 Real xblob;
 Real yblob;
 Real zblob;
 Real radblob;

 Real xblob2;
 Real yblob2;
 Real zblob2;
 Real radblob2;

 Real xblob3;
 Real yblob3;
 Real zblob3;
 Real radblob3;

 Real xblob4;
 Real yblob4;
 Real zblob4;
 Real radblob4;

 Real xblob5;
 Real yblob5;
 Real zblob5;
 Real radblob5;

 Real xblob6;
 Real yblob6;
 Real zblob6;
 Real radblob6;

 Real xblob7;
 Real yblob7;
 Real zblob7;
 Real radblob7;

 Real xblob8;
 Real yblob8;
 Real zblob8;
 Real radblob8;

 Real xblob9;
 Real yblob9;
 Real zblob9;
 Real radblob9;

 Real xblob10;
 Real yblob10;
 Real zblob10;
 Real radblob10;

 Real xactive;
 Real yactive;
 Real zactive;
 Real ractive;
 Real ractivex;
 Real ractivey;
 Real ractivez;

 Real fort_stop_time=-1.0;

 int probtype;
 int adv_dir;
 Real adv_vel;
 int axis_dir;
 Real rgasinlet;
 Real vinletgas;
 Real twall;
 Real advbot;
 Real inflow_pressure=NavierStokes::inflow_pressure;
 Real outflow_pressure=NavierStokes::outflow_pressure;
 Real period_time=NavierStokes::period_time;

 ParmParse ppmain;
 fort_stop_time=-1.0;
 ppmain.queryAdd("stop_time",fort_stop_time);

 int ns_max_level;
 ParmParse ppamr("amr");
 ppamr.get("max_level",ns_max_level);
 Vector<int> ns_space_blocking_factor;
 ns_space_blocking_factor.resize(ns_max_level+1);
 for (int lev=0;lev<=ns_max_level;lev++)
  ns_space_blocking_factor[lev]=1;
 ppamr.queryAdd("space_blocking_factor",
   ns_space_blocking_factor,ns_max_level+1);

 int time_blocking_factor=1;
 ppamr.queryAdd("time_blocking_factor",time_blocking_factor); 

 Vector<int> ns_n_cell(AMREX_SPACEDIM);
 ppamr.getarr("n_cell",ns_n_cell,0,AMREX_SPACEDIM);

 ParmParse pp("ns");
 pp.get("probtype",probtype);
 pp.get("axis_dir",axis_dir);
 pp.get("zblob",zblob);

 extra_circle_parameters(
   xblob2,yblob2,zblob2,radblob2,
   xblob3,yblob3,zblob3,radblob3,
   xblob4,yblob4,zblob4,radblob4,
   xblob5,yblob5,zblob5,radblob5,
   xblob6,yblob6,zblob6,radblob6,
   xblob7,yblob7,zblob7,radblob7,
   xblob8,yblob8,zblob8,radblob8,
   xblob9,yblob9,zblob9,radblob9,
   xblob10,yblob10,zblob10,radblob10 );

 pp.get("denfact",denfact);
 pp.get("velfact",velfact);

 pp.get("xblob",xblob);
 pp.get("yblob",yblob);
 pp.get("zblob",zblob);
 pp.get("radblob",radblob);

 pp.queryAdd("xblob2",xblob2);
 pp.queryAdd("yblob2",yblob2);
 pp.queryAdd("zblob2",zblob2);
 pp.queryAdd("radblob2",radblob2);

 pp.queryAdd("xblob3",xblob3);
 pp.queryAdd("yblob3",yblob3);
 pp.queryAdd("zblob3",zblob3);
 pp.queryAdd("radblob3",radblob3);

 pp.queryAdd("xblob4",xblob4);
 pp.queryAdd("yblob4",yblob4);
 pp.queryAdd("zblob4",zblob4);
 pp.queryAdd("radblob4",radblob4);

 pp.queryAdd("xblob5",xblob5);
 pp.queryAdd("yblob5",yblob5);
 pp.queryAdd("zblob5",zblob5);
 pp.queryAdd("radblob5",radblob5);

 pp.queryAdd("xblob6",xblob6);
 pp.queryAdd("yblob6",yblob6);
 pp.queryAdd("zblob6",zblob6);
 pp.queryAdd("radblob6",radblob6);

 pp.queryAdd("xblob7",xblob7);
 pp.queryAdd("yblob7",yblob7);
 pp.queryAdd("zblob7",zblob7);
 pp.queryAdd("radblob7",radblob7);

 pp.queryAdd("xblob8",xblob8);
 pp.queryAdd("yblob8",yblob8);
 pp.queryAdd("zblob8",zblob8);
 pp.queryAdd("radblob8",radblob8);

 pp.queryAdd("xblob9",xblob9);
 pp.queryAdd("yblob9",yblob9);
 pp.queryAdd("zblob9",zblob9);
 pp.queryAdd("radblob9",radblob9);

 pp.queryAdd("xblob10",xblob10);
 pp.queryAdd("yblob10",yblob10);
 pp.queryAdd("zblob10",zblob10);
 pp.queryAdd("radblob10",radblob10);

 xactive=0.0;
 yactive=0.0;
 zactive=0.0;
 ractive=0.0;
 ractivex=0.0;
 ractivey=0.0;
 ractivez=0.0;

 pp.queryAdd("xactive",xactive);
 pp.queryAdd("yactive",yactive);
 pp.queryAdd("zactive",zactive);
 pp.queryAdd("ractive",ractive);
 if (ractive>0.0) {
  ractivex=ractive;
  ractivey=ractive;
  ractivez=ractive;
 }
 pp.queryAdd("ractivex",ractivex);
 pp.queryAdd("ractivey",ractivey);
 pp.queryAdd("ractivez",ractivez);

 pp.get("adv_dir",adv_dir);
 if ((adv_dir<1)||(adv_dir>2*AMREX_SPACEDIM+1))
  amrex::Error("adv_dir invalid");

 pp.get("adv_vel",adv_vel);
 pp.get("rgasinlet",rgasinlet);
 pp.get("vinletgas",vinletgas);
 pp.get("twall",twall);
 pp.get("advbot",advbot);
 pp.queryAdd("inflow_pressure",inflow_pressure);
 pp.queryAdd("outflow_pressure",outflow_pressure);
 pp.queryAdd("period_time",period_time);

 NavierStokes::num_state_material=ENUM_SPECIESVAR;  // den,T
 NavierStokes::num_state_base=ENUM_SPECIESVAR;  // den,T
 NavierStokes::ngeom_raw=AMREX_SPACEDIM+1;
 NavierStokes::ngeom_recon=ENUM_NUM_MOF_VAR;

 pp.get("num_materials",NavierStokes::num_materials);
 if ((NavierStokes::num_materials<2)||(NavierStokes::num_materials>999))
  amrex::Error("num materials invalid");

 NavierStokes::num_interfaces=
  ( (NavierStokes::num_materials-1)*
    (NavierStokes::num_materials-1)+
    NavierStokes::num_materials-1 )/2;

 if ((NavierStokes::num_interfaces<1)||(NavierStokes::num_interfaces>999))
  amrex::Error("num interfaces invalid");

  // this is local variable, not static variable
 int MOFITERMAX=DEFAULT_MOFITERMAX;  
 pp.queryAdd("MOFITERMAX",MOFITERMAX);
 if ((MOFITERMAX<0)||(MOFITERMAX>MOFITERMAX_LIMIT)) {
  std::cout << "MOFITERMAX= " << MOFITERMAX << '\n';
  std::cout << "MOFITERMAX_LIMIT= " << MOFITERMAX_LIMIT << '\n';
  amrex::Error("mof iter max invalid in navierstokes");
 }

  // this is local variable, not static variable
 int MOFITERMAX_AFTER_PREDICT=DEFAULT_MOFITERMAX_AFTER_PREDICT;  
 pp.queryAdd("MOFITERMAX_AFTER_PREDICT",MOFITERMAX_AFTER_PREDICT);
 if ((MOFITERMAX_AFTER_PREDICT<0)|| 
     (MOFITERMAX_AFTER_PREDICT>MOFITERMAX_LIMIT)||
     (MOFITERMAX_AFTER_PREDICT>MOFITERMAX)) {
  std::cout << "MOFITERMAX_AFTER_PREDICT= " << 
    MOFITERMAX_AFTER_PREDICT << '\n';
  std::cout << "MOFITERMAX_LIMIT= " << MOFITERMAX_LIMIT << '\n';
  amrex::Error("mof iter max after predict invalid in navierstokes");
 }
 int MOF_TURN_OFF_LS=NavierStokes::MOF_TURN_OFF_LS;
 pp.queryAdd("MOF_TURN_OFF_LS",MOF_TURN_OFF_LS);
 if ((MOF_TURN_OFF_LS!=0)&&(MOF_TURN_OFF_LS!=1))
  amrex::Error("mof turn off ls invalid in navierstokes");

 int MOF_DEBUG_RECON=NavierStokes::MOF_DEBUG_RECON; 
 pp.queryAdd("MOF_DEBUG_RECON",MOF_DEBUG_RECON);
 if ((MOF_DEBUG_RECON!=0)&&(MOF_DEBUG_RECON!=1)&&
     (MOF_DEBUG_RECON!=2))
  amrex::Error("mof debug recon invalid in navierstokes");

 pp.get("num_species_var",NavierStokes::num_species_var);
 if (NavierStokes::num_species_var<0)
  amrex::Error("num species var invalid");

 NavierStokes::num_state_base=ENUM_SPECIESVAR;   // den,Temperature
 NavierStokes::num_state_material=ENUM_SPECIESVAR;  // den,Temperature
 NavierStokes::num_state_material+=NavierStokes::num_species_var;

 NavierStokes::grid_stretching_parameter.resize(AMREX_SPACEDIM);
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) { 
  NavierStokes::grid_stretching_parameter[dir]=0.0;
 }
 pp.queryAdd("grid_stretching_parameter",
    NavierStokes::grid_stretching_parameter,
    AMREX_SPACEDIM);

 Vector<Real> elastic_viscosity_temp;
 Vector<Real> elastic_time_temp;
 Vector<int> viscoelastic_model_temp;

 elastic_viscosity_temp.resize(NavierStokes::num_materials);
 elastic_time_temp.resize(NavierStokes::num_materials);
 viscoelastic_model_temp.resize(NavierStokes::num_materials);

 NavierStokes::store_elastic_data.resize(NavierStokes::num_materials);
 for (int im=0;im<NavierStokes::num_materials;im++) {

  elastic_viscosity_temp[im]=0.0;
  elastic_time_temp[im]=0.0;
  viscoelastic_model_temp[im]=0;

  NavierStokes::store_elastic_data[im]=0;
 }
 pp.queryAdd("elastic_viscosity",elastic_viscosity_temp,
   NavierStokes::num_materials);
 pp.queryAdd("elastic_time",elastic_time_temp,NavierStokes::num_materials);
 pp.queryAdd("viscoelastic_model",viscoelastic_model_temp,
	NavierStokes::num_materials);

 for (int im=0;im<NavierStokes::num_materials;im++) {
  if (elastic_viscosity_temp[im]>0.0) {
   if (fort_built_in_elastic_model(&elastic_viscosity_temp[im],
              		         &viscoelastic_model_temp[im])==1) {
    NavierStokes::store_elastic_data[im]=1;
   } else if (fort_built_in_elastic_model(&elastic_viscosity_temp[im],
                                        &viscoelastic_model_temp[im])==0) {
    // do nothing
   } else
    amrex::Error("fort_built_in_elastic_model invalid");
  } else if (elastic_viscosity_temp[im]==0.0) {
   // do nothing
  } else
   amrex::Error("elastic_viscosity_temp[im] invalid");
 } // im=0..NavierStokes::num_materials-1 

 NavierStokes::num_materials_viscoelastic=0;
 for (int im=0;im<NavierStokes::num_materials;im++) {
  if (NavierStokes::store_elastic_data[im]==1) {
   NavierStokes::num_materials_viscoelastic++;
  } else if (NavierStokes::store_elastic_data[im]==0) {
   // do nothing
  } else
   amrex::Error("NavierStokes::store_elastic_data invalid");
 } // im=0..NavierStokes::num_materials-1 

 Vector<Real> denconst_temp(NavierStokes::num_materials);
 Vector<Real> den_ceiling_temp(NavierStokes::num_materials);
 Vector<Real> den_floor_temp(NavierStokes::num_materials);
 Vector<Real> cavdenconst_temp(NavierStokes::num_materials);

 Vector<Real> stiffPINFtemp(NavierStokes::num_materials);
 Vector<Real> stiffCPtemp(NavierStokes::num_materials);
 Vector<Real> stiffCVtemp(NavierStokes::num_materials);
 Vector<Real> stiffGAMMAtemp(NavierStokes::num_materials);

 //Du/Dt=-grad (p-rho0 g dot z)/rho0 - g DrhoDT (T-T0) 
 //DrhoDT has units of 1/(Degrees Kelvin)
 Vector<Real> DrhoDTtemp(NavierStokes::num_materials);
 Vector<Real> tempcutofftemp(NavierStokes::num_materials);
 Vector<Real> tempcutoffmaxtemp(NavierStokes::num_materials);
 Vector<Real> tempconst_temp(NavierStokes::num_materials);
 Vector<Real> initial_temperature_temp(NavierStokes::num_materials);
 Vector<Real> viscconst_temp(NavierStokes::num_materials);

 Vector<Real> viscconst_eddy_wall_temp(NavierStokes::num_materials);
 Vector<Real> viscconst_eddy_bulk_temp(NavierStokes::num_materials);
 Vector<Real> heatviscconst_eddy_wall_temp(NavierStokes::num_materials);
 Vector<Real> heatviscconst_eddy_bulk_temp(NavierStokes::num_materials);

 Vector<Real> thermal_microlayer_size_temp(NavierStokes::num_materials);
 Vector<Real> shear_microlayer_size_temp(NavierStokes::num_materials);
 Vector<Real> buoyancy_microlayer_size_temp(NavierStokes::num_materials);
 Vector<Real> phasechange_microlayer_size_temp(NavierStokes::num_materials);

 Vector<int> viscosity_state_model_temp(NavierStokes::num_materials);
 Vector<Real> heatflux_factor_temp(NavierStokes::num_materials);
 Vector<Real> heatviscconst_temp(NavierStokes::num_materials);
 Vector<Real> speciesconst_temp((NavierStokes::num_species_var+1)*NavierStokes::num_materials);
 Vector<Real> speciesviscconst_temp((NavierStokes::num_species_var+1)*NavierStokes::num_materials);

 NavierStokes::material_type.resize(NavierStokes::num_materials);

 NavierStokes::material_type_interface.resize(NavierStokes::num_interfaces);

 NavierStokes::FSI_flag.resize(NavierStokes::num_materials);

 NavierStokes::CTML_FSI_numsolids=0;

 NavierStokes::CTML_max_num_nodes_list.resize(3);
 for (int dir=0;dir<3;dir++) {
  NavierStokes::CTML_max_num_nodes_list[dir]=0;
 }

 NavierStokes::CTML_max_num_elements_list=0;
 NavierStokes::CTML_FSI_num_scalars=0;

 Vector<Real> Carreau_alpha_temp(NavierStokes::num_materials);
 Vector<Real> Carreau_beta_temp(NavierStokes::num_materials);
 Vector<Real> Carreau_n_temp(NavierStokes::num_materials);
 Vector<Real> Carreau_mu_inf_temp(NavierStokes::num_materials);
 Vector<int> shear_thinning_fluid_temp(NavierStokes::num_materials);
 Vector<Real> polymer_factor_temp(NavierStokes::num_materials);

 Vector<Real> concentration_temp(NavierStokes::num_materials);
 Vector<Real> etaL_temp(NavierStokes::num_materials);
 Vector<Real> etaS_temp(NavierStokes::num_materials);
 Vector<Real> etaP_temp(NavierStokes::num_materials);

 Real visc_coef_temp=NavierStokes::visc_coef;
 pp.get("visc_coef",visc_coef_temp);

 Real angular_velocity_temp=NavierStokes::angular_velocity;
 pp.queryAdd("angular_velocity",angular_velocity_temp);

 pp.getarr("material_type",NavierStokes::material_type,0,
    NavierStokes::num_materials);

 for (int im=0;im<NavierStokes::num_materials;im++) {

  stiffPINFtemp[im]=0.0;
  stiffCPtemp[im]=4.1855e+7;
  stiffCVtemp[im]=4.1855e+7;
  stiffGAMMAtemp[im]=1.4;

  DrhoDTtemp[im]=0.0;
  tempcutofftemp[im]=1.0e-8;
  tempcutoffmaxtemp[im]=1.0e+99;
  NavierStokes::FSI_flag[im]=FSI_FLUID;

  Carreau_alpha_temp[im]=1.0;
  Carreau_beta_temp[im]=0.0;
  Carreau_n_temp[im]=1.0;
  Carreau_mu_inf_temp[im]=0.0;
  shear_thinning_fluid_temp[im]=0;

  polymer_factor_temp[im]=0.0;
 } // im=0..NavierStokes::num_materials-1

 pp.queryAdd("polymer_factor",polymer_factor_temp,NavierStokes::num_materials);

 pp.queryAdd("Carreau_alpha",Carreau_alpha_temp,NavierStokes::num_materials);
 pp.queryAdd("Carreau_beta",Carreau_beta_temp,NavierStokes::num_materials);
 pp.queryAdd("Carreau_n",Carreau_n_temp,NavierStokes::num_materials);
 pp.queryAdd("Carreau_mu_inf",Carreau_mu_inf_temp,NavierStokes::num_materials);

 for (int im=0;
      im<(NavierStokes::num_species_var+1)*NavierStokes::num_materials;
      im++) {
  speciesviscconst_temp[im]=0.0;
  speciesconst_temp[im]=0.0;
 }

 pp.queryAdd("FSI_flag",NavierStokes::FSI_flag,NavierStokes::num_materials);

  // fort_ctml_max_nodes is declared in: CTMLFSI.F90
#ifdef MVAHABFSI
 fort_ctml_max_nodes(
   &NavierStokes::num_materials,
   NavierStokes::FSI_flag.dataPtr(),
   &NavierStokes::CTML_FSI_numsolids,
   &NavierStokes::CTML_FSI_num_scalars,
   NavierStokes::CTML_max_num_nodes_list.dataPtr(),
   &NavierStokes::CTML_max_num_elements_list);
#endif

 int num_local_aux_grids_temp=NavierStokes::num_local_aux_grids;
 pp.queryAdd("num_local_aux_grids",num_local_aux_grids_temp);

 pp.queryAdd("tempcutoff",tempcutofftemp,NavierStokes::num_materials);
 pp.queryAdd("tempcutoffmax",tempcutoffmaxtemp,NavierStokes::num_materials);

 pp.getarr("tempconst",tempconst_temp,0,NavierStokes::num_materials);
 for (int im=0;im<NavierStokes::num_materials;im++)
  initial_temperature_temp[im]=tempconst_temp[im];
 pp.queryAdd("initial_temperature",initial_temperature_temp,
   NavierStokes::num_materials);

  //Du/Dt=-grad (p-rho0 g dot z)/rho0 - g DrhoDT (T-T0) 
  //DrhoDT has units of 1/(Degrees Kelvin)
 pp.queryAdd("DrhoDT",DrhoDTtemp,NavierStokes::num_materials);

 pp.queryAdd("stiffPINF",stiffPINFtemp,NavierStokes::num_materials);

 pp.queryAdd("stiffCP",stiffCPtemp,NavierStokes::num_materials);
 for (int im=0;im<NavierStokes::num_materials;im++)
  stiffCVtemp[im]=stiffCPtemp[im];
 pp.queryAdd("stiffCV",stiffCVtemp,NavierStokes::num_materials);

 Vector<Real> prerecalesce_stiffCP_temp(NavierStokes::num_materials);
 for (int im=0;im<NavierStokes::num_materials;im++)
  prerecalesce_stiffCP_temp[im]=stiffCPtemp[im];
 pp.queryAdd("precalesce_stiffCP",prerecalesce_stiffCP_temp,
    NavierStokes::num_materials);
 Vector<Real> prerecalesce_stiffCV_temp(NavierStokes::num_materials);
 for (int im=0;im<NavierStokes::num_materials;im++)
  prerecalesce_stiffCV_temp[im]=stiffCVtemp[im];
 pp.queryAdd("precalesce_stiffCV",prerecalesce_stiffCV_temp,
  NavierStokes::num_materials);

 pp.queryAdd("stiffGAMMA",stiffGAMMAtemp,NavierStokes::num_materials);

 pp.getarr("denconst",denconst_temp,0,NavierStokes::num_materials);

 for (int im=0;im<NavierStokes::num_materials;im++) {
  cavdenconst_temp[im]=0.0;
  den_ceiling_temp[im]=1.0e+20;
  den_floor_temp[im]=0.0;
 }
 pp.queryAdd("cavitation_vapor_density",cavdenconst_temp,
     NavierStokes::num_materials);
 pp.queryAdd("density_floor",den_floor_temp,NavierStokes::num_materials);
 pp.queryAdd("density_ceiling",den_ceiling_temp,NavierStokes::num_materials);

 pp.getarr("viscconst",viscconst_temp,0,NavierStokes::num_materials);
 NavierStokes::viscconst_min=viscconst_temp[0];
 NavierStokes::viscconst_max=viscconst_temp[0];

 for (int im=0;im<NavierStokes::num_materials;im++) {

  if (viscconst_temp[im]<NavierStokes::viscconst_min)
   NavierStokes::viscconst_min=viscconst_temp[im];
  if (viscconst_temp[im]>NavierStokes::viscconst_max)
   NavierStokes::viscconst_max=viscconst_temp[im];

  viscconst_eddy_wall_temp[im]=0.0;
  viscconst_eddy_bulk_temp[im]=0.0;
  heatviscconst_eddy_wall_temp[im]=0.0;
  heatviscconst_eddy_bulk_temp[im]=0.0;

  thermal_microlayer_size_temp[im]=microlayer_size_default;
  shear_microlayer_size_temp[im]=microlayer_size_default;
  buoyancy_microlayer_size_temp[im]=microlayer_size_default;
  phasechange_microlayer_size_temp[im]=microlayer_size_default;
 } //im=0...num_materials-1

 pp.queryAdd("viscconst_eddy_wall",viscconst_eddy_wall_temp,
   NavierStokes::num_materials);
 pp.queryAdd("viscconst_eddy_bulk",viscconst_eddy_bulk_temp,
   NavierStokes::num_materials);
 pp.queryAdd("heatviscconst_eddy_wall",heatviscconst_eddy_wall_temp,
   NavierStokes::num_materials);
 pp.queryAdd("heatviscconst_eddy_bulk",heatviscconst_eddy_bulk_temp,
   NavierStokes::num_materials);

 pp.queryAdd("thermal_microlayer_size",thermal_microlayer_size_temp,NavierStokes::num_materials);
 pp.queryAdd("shear_microlayer_size",shear_microlayer_size_temp,NavierStokes::num_materials);
 pp.queryAdd("buoyancy_microlayer_size",buoyancy_microlayer_size_temp,NavierStokes::num_materials);
 pp.queryAdd("phasechange_microlayer_size", 
   phasechange_microlayer_size_temp,NavierStokes::num_materials);

 Vector<Real> prerecalesce_viscconst_temp(NavierStokes::num_materials);
 for (int im=0;im<NavierStokes::num_materials;im++)
  prerecalesce_viscconst_temp[im]=viscconst_temp[im];
 pp.queryAdd("precalesce_viscconst",prerecalesce_viscconst_temp,
	NavierStokes::num_materials);

 for (int im=0;im<NavierStokes::num_materials;im++) {
  viscosity_state_model_temp[im]=0;
  heatflux_factor_temp[im]=1.0;
 }
 pp.queryAdd("viscosity_state_model",
  viscosity_state_model_temp,NavierStokes::num_materials);

 pp.queryAdd("heatflux_factor",heatflux_factor_temp,
	NavierStokes::num_materials);
 pp.getarr("heatviscconst",heatviscconst_temp,0,NavierStokes::num_materials);

 Vector<Real> prerecalesce_heatviscconst_temp(NavierStokes::num_materials);
 for (int im=0;im<NavierStokes::num_materials;im++)
  prerecalesce_heatviscconst_temp[im]=heatviscconst_temp[im];
 pp.queryAdd("precalesce_heatviscconst",prerecalesce_heatviscconst_temp,
  NavierStokes::num_materials);

 if (NavierStokes::num_species_var>0) {
  pp.queryAdd("speciesconst",speciesconst_temp,
    NavierStokes::num_species_var*NavierStokes::num_materials);
  pp.queryAdd("speciesviscconst",speciesviscconst_temp,
    NavierStokes::num_species_var*NavierStokes::num_materials);
 }

 Vector<Real> tension_slopetemp(NavierStokes::num_interfaces);
 Vector<Real> tension_T0temp(NavierStokes::num_interfaces);
 Vector<Real> tension_mintemp(NavierStokes::num_interfaces);
 for (int im=0;im<NavierStokes::num_interfaces;im++) {
  tension_slopetemp[im]=0.0;
  tension_T0temp[im]=293.0;
  tension_mintemp[im]=0.0;
 }

 Vector<Real> tensiontemp(NavierStokes::num_interfaces);
 Vector<Real> tension_inittemp(NavierStokes::num_interfaces);

 Vector<Real> prefreeze_tensiontemp(NavierStokes::num_interfaces);

 pp.getarr("tension",tensiontemp,0,NavierStokes::num_interfaces);
 for (int im=0;im<NavierStokes::num_interfaces;im++) {
  tension_inittemp[im]=tensiontemp[im];
 }
 pp.queryarr("tension_init",tension_inittemp,0,NavierStokes::num_interfaces);

 pp.queryAdd("tension_slope",tension_slopetemp,NavierStokes::num_interfaces);
 pp.queryAdd("tension_T0",tension_T0temp,NavierStokes::num_interfaces);
 pp.queryAdd("tension_min",tension_mintemp,NavierStokes::num_interfaces);

 Vector<Real> latent_heat_temp(2*NavierStokes::num_interfaces);
 Vector<Real> latent_heat_slopetemp(2*NavierStokes::num_interfaces);
 Vector<Real> latent_heat_T0temp(2*NavierStokes::num_interfaces);
 Vector<Real> latent_heat_mintemp(2*NavierStokes::num_interfaces);

 Vector<Real> saturation_temp_temp(2*NavierStokes::num_interfaces);
 Vector<Real> reference_pressure_temp(2*NavierStokes::num_interfaces);
 for (int i=0;i<NavierStokes::num_interfaces;i++) { 

  if (tension_slopetemp[i]<=0.0) {
   // do nothing
  } else
   amrex::Error("tension_slope must be non-positive(1)");

  saturation_temp_temp[i]=0.0;
  saturation_temp_temp[i+NavierStokes::num_interfaces]=0.0;
  reference_pressure_temp[i]=1.0e+6;
  reference_pressure_temp[i+NavierStokes::num_interfaces]=1.0e+6;

  latent_heat_temp[i]=0.0;
  latent_heat_temp[i+NavierStokes::num_interfaces]=0.0;
  latent_heat_slopetemp[i]=0.0;
  latent_heat_slopetemp[i+NavierStokes::num_interfaces]=0.0;
  latent_heat_T0temp[i]=0.0;
  latent_heat_T0temp[i+NavierStokes::num_interfaces]=0.0;
  latent_heat_mintemp[i]=0.0;
  latent_heat_mintemp[i+NavierStokes::num_interfaces]=0.0;
 }

 pp.queryAdd("latent_heat",latent_heat_temp,2*NavierStokes::num_interfaces);
 pp.queryAdd("latent_heat_slope",latent_heat_slopetemp,
    2*NavierStokes::num_interfaces);
 pp.queryAdd("latent_heat_T0",latent_heat_T0temp,
    2*NavierStokes::num_interfaces);
 pp.queryAdd("latent_heat_min",latent_heat_mintemp,
    2*NavierStokes::num_interfaces);

 for (int i=0;i<2*NavierStokes::num_interfaces;i++) { 

  if (latent_heat_slopetemp[i]<=0.0) {
   // do nothing
  } else
   amrex::Error("latent_heat_slope must be non-positive(1)");

 } //i=0...2*NavierStokes::num_interfaces-1

 pp.queryAdd("saturation_temp",saturation_temp_temp,
	2*NavierStokes::num_interfaces);
 pp.queryAdd("reference_pressure",reference_pressure_temp,
	2*NavierStokes::num_interfaces);

  // ergs/(mol Kelvin)
 pp.queryAdd("R_Palmore_Desjardins",NavierStokes::R_Palmore_Desjardins);

 Vector<Real> molar_mass_temp(NavierStokes::num_materials);
 Vector<Real> species_molar_mass_temp(NavierStokes::num_species_var+1);
 for (int im=0;im<NavierStokes::num_materials;im++) {
  molar_mass_temp[im]=1.0;
 }
 for (int im=0;im<NavierStokes::num_species_var+1;im++) {
  species_molar_mass_temp[im]=1.0;
 }
 pp.queryAdd("molar_mass",molar_mass_temp,NavierStokes::num_materials);

 pp.queryAdd("species_molar_mass",
   species_molar_mass_temp,NavierStokes::num_species_var);

 for (int im=0;im<NavierStokes::num_interfaces;im++) {
  prefreeze_tensiontemp[im]=tensiontemp[im];
 }
 pp.queryAdd("prefreeze_tension",prefreeze_tensiontemp,
	NavierStokes::num_interfaces);

 Vector<int> preset_flag;
 preset_flag.resize(NavierStokes::num_interfaces);

 for (int im=0;im<NavierStokes::num_materials;im++) {
  for (int im_opp=im+1;im_opp<NavierStokes::num_materials;im_opp++) {
   int iten=0;
   NavierStokes::get_iten_cpp(im+1,im_opp+1,iten);
   if ((iten<1)||(iten>NavierStokes::num_interfaces))
    amrex::Error("iten invalid");
   preset_flag[iten-1]=-1;
   if ((NavierStokes::material_type[im]==999)||
       (NavierStokes::material_type[im_opp]==999)) {
    NavierStokes::material_type_interface[iten-1]=999;
    preset_flag[iten-1]=999;
   } else if ((NavierStokes::material_type[im]==0)||
              (NavierStokes::material_type[im_opp]==0)) {
    NavierStokes::material_type_interface[iten-1]=0;
    preset_flag[iten-1]=0;
   } else if ((latent_heat_temp[iten-1]!=0.0)||
  	      (latent_heat_temp[iten-1+NavierStokes::num_interfaces]!=0.0)) {
    NavierStokes::material_type_interface[iten-1]=0;
    preset_flag[iten-1]=0;
   } else {
    NavierStokes::material_type_interface[iten-1]= 
      NavierStokes::material_type[im];
   }
  } // im_opp=im+1;im_opp<NavierStokes::num_materials;im_opp++
 } // im=0;im<NavierStokes::num_materials;im++

 pp.queryAdd("material_type_interface",
   NavierStokes::material_type_interface,NavierStokes::num_interfaces);

 for (int im=0;im<NavierStokes::num_materials;im++) {
  for (int im_opp=im+1;im_opp<NavierStokes::num_materials;im_opp++) {
   int iten=0;
   NavierStokes::get_iten_cpp(im+1,im_opp+1,iten);
   if ((iten<1)||(iten>NavierStokes::num_interfaces))
    amrex::Error("iten invalid");
   if (preset_flag[iten-1]==-1) {
    if ((NavierStokes::material_type_interface[iten-1]>=0)&&
        (NavierStokes::material_type_interface[iten-1]<=999)) {
     //do nothing
    } else
     amrex::Error("NavierStokes::material_type_interface invalid");
   } else if (preset_flag[iten-1]>=0) {
    if (NavierStokes::material_type_interface[iten-1]==preset_flag[iten-1]) {
     //do nothing
    } else
     amrex::Error("NavierStokes::material_type_interface invalid");
   } else
    amrex::Error("preset flag invalid");
  } // im_opp=im+1;im_opp<NavierStokes::num_materials;im_opp++
 } // im=0;im<NavierStokes::num_materials;im++

 int ioproc=0;
 if (ParallelDescriptor::IOProcessor())
  ioproc=1;

 const int cc_int_size=sizeof(int);

  // declared in PROB_CPP_PARMS.F90
 fort_override_MAIN_GLOBALS(
  &cc_int_size,
  &NavierStokes::num_species_var,
  &NavierStokes::num_materials_viscoelastic,
  &NavierStokes::num_state_material,
  &NavierStokes::num_state_base,
  &NavierStokes::ngeom_raw,
  &NavierStokes::ngeom_recon,
  &NavierStokes::num_materials,
  &NavierStokes::num_interfaces,
  &ioproc);

 ParallelDescriptor::Barrier();

 for (int im=0;im<NavierStokes::num_materials;im++) {

  if (NavierStokes::material_type[im]==999) {

    // non-tessellating cases.
   if ((NavierStokes::FSI_flag[im]==FSI_PRESCRIBED_PROBF90)||
       (NavierStokes::FSI_flag[im]==FSI_PRESCRIBED_NODES)||
       (NavierStokes::FSI_flag[im]==FSI_SHOELE_CTML)) {
    //do nothing
   } else
    amrex::Error("NavierStokes::FSI_flag invalid");

   int imp1=im+1;
   if (fort_is_rigid_base(&NavierStokes::FSI_flag[im],&imp1)==1) {
    //do nothing
   } else if (fort_is_rigid_base(&NavierStokes::FSI_flag[im],&imp1)==0) {
    amrex::Error("NavierStokes::FSI_flag and material_type inconsistent");
   } else
    amrex::Error("fort_is_rigid_base corrupt");

  } else if (NavierStokes::material_type[im]==0) {

   if ((NavierStokes::FSI_flag[im]==FSI_FLUID)||
       (NavierStokes::FSI_flag[im]==FSI_FLUID_NODES_INIT)||
       (NavierStokes::FSI_flag[im]==FSI_ICE_PROBF90)||
       (NavierStokes::FSI_flag[im]==FSI_ICE_STATIC)||
       (NavierStokes::FSI_flag[im]==FSI_ICE_NODES_INIT)||
       (NavierStokes::FSI_flag[im]==FSI_RIGID_NOTPRESCRIBED)) {
    //do nothing
   } else
    amrex::Error("NavierStokes::FSI_flag invalid");

   int imp1=im+1;
   if (fort_is_rigid_base(&NavierStokes::FSI_flag[im],&imp1)==0) {
    //do nothing
   } else if (fort_is_rigid_base(&NavierStokes::FSI_flag[im],&imp1)==1) {
    amrex::Error("NavierStokes::FSI_flag and material_type inconsistent");
   } else
    amrex::Error("fort_is_rigid_base corrupt");

  } else if ((NavierStokes::material_type[im]>0)&& 
             (NavierStokes::material_type[im]<999)) {

   if ((NavierStokes::FSI_flag[im]==FSI_FLUID)||
       (NavierStokes::FSI_flag[im]==FSI_FLUID_NODES_INIT)) {
    //do nothing
   } else
    amrex::Error("NavierStokes::FSI_flag invalid");

  } else {
   amrex::Error("material type invalid");
  }

  int imp1=im+1;
  if (fort_is_rigid_base(&NavierStokes::FSI_flag[im],&imp1)==1) {
   shear_thinning_fluid_temp[im]=0;
  } else if (fort_is_rigid_base(&NavierStokes::FSI_flag[im],&imp1)==0) {
   shear_thinning_fluid_temp[im]=0;
   if ((probtype==2)&&(axis_dir>0)&&(im==0))
    shear_thinning_fluid_temp[im]=1;
   if (Carreau_beta_temp[im]!=0.0)
    shear_thinning_fluid_temp[im]=1;

   if (shear_thinning_fluid_temp[im]==1) {
    // do nothing
   } else if (shear_thinning_fluid_temp[im]==0) {
    // do nothing
   } else
    amrex::Error("shear_thinning_fluid invalid");
  } else 
   amrex::Error("fort_is_rigid_base invalid");

  etaL_temp[im]=viscconst_temp[im];  
  etaP_temp[im]=elastic_viscosity_temp[im]; //eta_P0
  etaS_temp[im]=etaL_temp[im]-etaP_temp[im];  

   // c0=etaP0/etaS=etaP0/(etaL0-etaP0)
  concentration_temp[im]=0.0;
  if (etaL_temp[im]-etaP_temp[im]>0.0)
   concentration_temp[im]=etaP_temp[im]/(etaL_temp[im]-etaP_temp[im]);  

 } //im=0..NavierStokes::num_materials-1

 if (NavierStokes::num_state_base!=2)
  amrex::Error("NavierStokes::num_state_base invalid 9");

 int prescribe_temperature_outflow=NavierStokes::prescribe_temperature_outflow;
 pp.queryAdd("prescribe_temperature_outflow",prescribe_temperature_outflow);
 if ((prescribe_temperature_outflow<0)||
     (prescribe_temperature_outflow>3))
  amrex::Error("prescribe_temperature_outflow invalid (fortran_parameters)");

  // 0=diffuse in solid 1=dirichlet 2=neumann
 int solidheat_flag=NavierStokes::solidheat_flag;
 pp.queryAdd("solidheat_flag",solidheat_flag);
 if ((solidheat_flag<0)||
     (solidheat_flag>2))
  amrex::Error("solidheat_flag invalid (fortran_parameters)");

 Real gravity_temp=NavierStokes::gravity; 
 int gravity_dir_temp=NavierStokes::gravity_dir;
 int invert_gravity_temp=NavierStokes::invert_gravity;
 Vector<Real> gravity_vector_temp(AMREX_SPACEDIM);

 bool gravity_in_table=pp.contains("gravity");
 bool gravity_dir_in_table=pp.contains("gravity_dir");
 bool invert_gravity_in_table=pp.contains("invert_gravity");
 bool gravity_vector_in_table=pp.contains("gravity_vector");

 if (gravity_vector_in_table==true) {

  if ((gravity_in_table==false)&&
      (gravity_dir_in_table==false)&&
      (invert_gravity_in_table==false)) {
   // do nothing
  } else
   amrex::Error("gravity parm conflict");

  pp.getarr("gravity_vector",gravity_vector_temp);

 } else if (gravity_vector_in_table==false) {

  pp.query("gravity",gravity_temp);
  pp.query("gravity_dir",gravity_dir_temp);
  pp.query("invert_gravity",invert_gravity_temp);
  if ((gravity_dir_temp<1)||(gravity_dir_temp>AMREX_SPACEDIM))
   amrex::Error("gravity dir invalid");

  for (int dir=0;dir<AMREX_SPACEDIM;dir++)
   gravity_vector_temp[dir]=0.0;

  if (invert_gravity_temp==0) {
   gravity_vector_temp[gravity_dir_temp-1]=-std::abs(gravity_temp);
  } else if (invert_gravity_temp==1) {
   gravity_vector_temp[gravity_dir_temp-1]=std::abs(gravity_temp);
  } else
   amrex::Error("invert_gravity_temp invalid");

 } else
  amrex::Error("gravity_vector_in_table invalid");

 int n_sites=0;
 pp.queryAdd("n_sites",n_sites);
 Real nucleation_init_time=NavierStokes::nucleation_init_time;
 pp.queryAdd("nucleation_init_time",nucleation_init_time);

 pp.queryAdd("ZEYU_DCA_SELECT",NavierStokes::ZEYU_DCA_SELECT);
 if ((NavierStokes::ZEYU_DCA_SELECT==-1)||  //static
     ((NavierStokes::ZEYU_DCA_SELECT>=1)&&
      (NavierStokes::ZEYU_DCA_SELECT<=8))) {
  // do nothing
 } else
  amrex::Error("NavierStokes::ZEYU_DCA_SELECT invalid");

 Vector<Real> temp_pos_sites(4);
 for (int dir=0;dir<4;dir++)
  temp_pos_sites[dir]=0.0;
 
 if (n_sites>0) {
  temp_pos_sites.resize(4*n_sites);
    // these positions are used only for the initial creation of bubbles;
    // at later times, when bubbles are nucleated, the position is random.
    // The code that sets the bubble positions and times is in 
    // NavierStokes::level_phase_change_rate (NavierStokes.cpp)
    // search "MITSUHIRO BUBBLE POSITIONS" for exact places where bubble 
    // positions are set.
  pp.getarr("pos_sites",temp_pos_sites,0,4*n_sites);
 } else if (n_sites==0) {
  // do nothing
 } else {
  std::cout << "n_sites= " << n_sites << '\n';
  amrex::Error("n_sites invalid(1)");
 }
 double start_initialization = ParallelDescriptor::second();

  // declared in PROB_CPP_PARMS.F90
 fort_override(
  &cc_int_size,
  &ns_max_level,
  ns_n_cell.dataPtr(),
  ns_space_blocking_factor.dataPtr(),
  &time_blocking_factor,
  &prescribe_temperature_outflow,
  &solidheat_flag,
  &geometry_coord,
  NavierStokes::FSI_flag.dataPtr(),
  &num_local_aux_grids_temp,
  &NavierStokes::ZEYU_DCA_SELECT,
  &denfact,
  &velfact,
  &n_sites,
  &nucleation_init_time,
  temp_pos_sites.dataPtr(),
  &xblob,&yblob,&zblob,&radblob,
  &xblob2,&yblob2,&zblob2,&radblob2,
  &xblob3,&yblob3,&zblob3,&radblob3,
  &xblob4,&yblob4,&zblob4,&radblob4,
  &xblob5,&yblob5,&zblob5,&radblob5,
  &xblob6,&yblob6,&zblob6,&radblob6,
  &xblob7,&yblob7,&zblob7,&radblob7,
  &xblob8,&yblob8,&zblob8,&radblob8,
  &xblob9,&yblob9,&zblob9,&radblob9,
  &xblob10,&yblob10,&zblob10,&radblob10,
  &xactive,&yactive,&zactive,
  &ractivex,&ractivey,&ractivez,
  &probtype,&adv_dir,&adv_vel,
  &axis_dir,&rgasinlet,
  &vinletgas,&twall,&advbot, 
  &inflow_pressure,
  &outflow_pressure,
  &period_time,
  &problox,&probloy,&probloz,
  &probhix,&probhiy,&probhiz,
  &NavierStokes::num_species_var,
  &NavierStokes::num_materials_viscoelastic,
  &NavierStokes::num_state_material,
  &NavierStokes::num_state_base,
  &NavierStokes::ngeom_raw,
  &NavierStokes::ngeom_recon,
  &NavierStokes::num_materials,
  NavierStokes::material_type.dataPtr(),
  NavierStokes::material_type_interface.dataPtr(),
  &NavierStokes::num_interfaces,
  DrhoDTtemp.dataPtr(),
  tempconst_temp.dataPtr(),
  initial_temperature_temp.dataPtr(),
  tempcutofftemp.dataPtr(),
  tempcutoffmaxtemp.dataPtr(),
  stiffPINFtemp.dataPtr(),
  &NavierStokes::R_Palmore_Desjardins,
  stiffCPtemp.dataPtr(),
  stiffCVtemp.dataPtr(),
  stiffGAMMAtemp.dataPtr(),
  denconst_temp.dataPtr(),
  den_floor_temp.dataPtr(),
  den_ceiling_temp.dataPtr(),
  cavdenconst_temp.dataPtr(),
  viscconst_temp.dataPtr(),
  viscconst_eddy_wall_temp.dataPtr(),
  viscconst_eddy_bulk_temp.dataPtr(),
  heatviscconst_eddy_wall_temp.dataPtr(),
  heatviscconst_eddy_bulk_temp.dataPtr(),
  thermal_microlayer_size_temp.dataPtr(),
  shear_microlayer_size_temp.dataPtr(),
  buoyancy_microlayer_size_temp.dataPtr(),
  phasechange_microlayer_size_temp.dataPtr(),
  viscosity_state_model_temp.dataPtr(),
  elastic_viscosity_temp.dataPtr(),
  elastic_time_temp.dataPtr(),
  viscoelastic_model_temp.dataPtr(),
  NavierStokes::store_elastic_data.dataPtr(),
  heatflux_factor_temp.dataPtr(),
  heatviscconst_temp.dataPtr(),
  prerecalesce_heatviscconst_temp.dataPtr(),
  prerecalesce_viscconst_temp.dataPtr(),
  prerecalesce_stiffCP_temp.dataPtr(),
  prerecalesce_stiffCV_temp.dataPtr(),
  speciesconst_temp.dataPtr(),
  speciesviscconst_temp.dataPtr(),
  latent_heat_temp.dataPtr(),
  latent_heat_slopetemp.dataPtr(),
  latent_heat_T0temp.dataPtr(),
  latent_heat_mintemp.dataPtr(),
  saturation_temp_temp.dataPtr(),
  reference_pressure_temp.dataPtr(),
  molar_mass_temp.dataPtr(),
  species_molar_mass_temp.dataPtr(),
  tensiontemp.dataPtr(),
  tension_inittemp.dataPtr(),
  tension_slopetemp.dataPtr(),
  tension_T0temp.dataPtr(),
  tension_mintemp.dataPtr(),
  prefreeze_tensiontemp.dataPtr(),
  gravity_vector_temp.dataPtr(),
  &fort_stop_time,
  Carreau_alpha_temp.dataPtr(),
  Carreau_beta_temp.dataPtr(),
  Carreau_n_temp.dataPtr(),
  Carreau_mu_inf_temp.dataPtr(),
  shear_thinning_fluid_temp.dataPtr(),
  polymer_factor_temp.dataPtr(),
  concentration_temp.dataPtr(),
  etaL_temp.dataPtr(),
  etaS_temp.dataPtr(),
  etaP_temp.dataPtr(),
  &visc_coef_temp,
  &angular_velocity_temp,
  NavierStokes::grid_stretching_parameter.dataPtr(),
  &ioproc);

 ParallelDescriptor::Barrier();

 double finish_initialization = ParallelDescriptor::second();
 if (ParallelDescriptor::IOProcessor()) {
  std::cout << "elapsed time in fort_override " << finish_initialization-
        start_initialization << '\n';
 }

 int mof_error_ordering_local=NavierStokes::mof_error_ordering;

 pp.queryAdd("mof_error_ordering",mof_error_ordering_local);
 if ((mof_error_ordering_local!=0)&& //centroid furthest from uncapt centroid
     (mof_error_ordering_local!=1))  //smallest MOF error
  amrex::Error("mof_error_ordering_local invalid");

 Vector<int> mof_ordering_local;
 mof_ordering_local.resize(NavierStokes::num_materials);

  //fort_mof_ordering_override is declared in: PROB_CPP_PARMS.F90
 fort_mof_ordering_override(
  mof_ordering_local.dataPtr(),
  &mof_error_ordering_local,
  NavierStokes::FSI_flag.dataPtr());

 pp.queryAdd("mof_ordering",mof_ordering_local,NavierStokes::num_materials);
 for (int i=0;i<NavierStokes::num_materials;i++) {
  if ((mof_ordering_local[i]<0)||
      (mof_ordering_local[i]>NavierStokes::num_materials+1))
   amrex::Error("mof_ordering_local invalid");
 }

 int temp_POLYGON_LIST_MAX=1000;
 
  //fort_initmof is declared in: MOF.F90
 fort_initmof(
   mof_ordering_local.dataPtr(),
   &MOFITERMAX,
   &MOFITERMAX_AFTER_PREDICT,
   &MOF_DEBUG_RECON,
   &MOF_TURN_OFF_LS,
   &thread_class::nthreads,
   &temp_POLYGON_LIST_MAX);

 if (ioproc==1) {
  std::cout << "in c++ code, after fort_override\n";
  for (int im=0;im<NavierStokes::num_materials;im++) {
   std::cout << "im= " << im << " mof_ordering_local= " <<
    mof_ordering_local[im] << '\n';
  }
 }
}  // end subroutine fortran_parameters()




void
NavierStokes::variableCleanUp ()
{
    desc_lst.clear();
    desc_lstGHOST.clear();
}

void
NavierStokes::read_geometry ()
{
    //
    // Must load coord here because 
    // 1. CoordSys hasn't read it in yet.
    // 2. geometry.coord_sys_override needs to be queried for too.
    //
    int geometry_coord;
    Vector<Real> geometry_prob_lo;
    Vector<Real> geometry_prob_hi;
    Vector<int> geometry_is_periodic;
    int geometry_is_any_periodic;

    read_geometry_raw(geometry_coord,geometry_prob_lo,geometry_prob_hi,
		    geometry_is_periodic,geometry_is_any_periodic);

    if (geometry_coord == COORDSYS_RZ) {
     if (AMREX_SPACEDIM==3)
      amrex::Error("No RZ in 3d");
    }

    if ((geometry_coord == COORDSYS_RZ) && 
        (phys_bc.lo(0) != Symmetry)) {

     phys_bc.setLo(0,Symmetry);
     viscosity_phys_bc.setLo(0,Symmetry);
     temperature_phys_bc.setLo(0,Symmetry);
     species_phys_bc.setLo(0,Symmetry);

     if (ParallelDescriptor::IOProcessor()) {
      std::cout << "\n WARNING: Setting phys_bc at xlo to Symmetry\n\n";
     }

    }

    NS_geometry_coord=geometry_coord;

} // end subroutine read_geometry

void
NavierStokes::setup_integrated_quantities() {

 NS_DRAG_integrated_quantities.resize(N_DRAG_IQ); 
 NS_DRAG_integrated_quantities_sweep.resize(N_DRAG_IQ); 
 for (int iq=0;iq<N_DRAG_IQ;iq++) {
  NS_DRAG_integrated_quantities[iq]=0.0;
  NS_DRAG_integrated_quantities_sweep[iq]=1;  // update second sweep
 }
 for (int iq=DRAGCOMP_IQ_COM;
      iq<DRAGCOMP_IQ_COM+3*num_materials;iq++) {
  NS_DRAG_integrated_quantities_sweep[iq]=0;  // update first sweep
 }
 for (int iq=DRAGCOMP_IQ_MASS;iq<DRAGCOMP_IQ_MASS+num_materials;iq++) {
  NS_DRAG_integrated_quantities_sweep[iq]=0;  // update first sweep
 }

 NS_sumdata.resize(IQ_TOTAL_SUM_COMP);
 NS_sumdata_type.resize(IQ_TOTAL_SUM_COMP);
 NS_sumdata_sweep.resize(IQ_TOTAL_SUM_COMP);

 for (int isum=0;isum<IQ_TOTAL_SUM_COMP;isum++) {
  NS_sumdata[isum]=0.0;
  NS_sumdata_type[isum]=1;  // reduce real sum
  NS_sumdata_sweep[isum]=0;  // update first sweep
 }

 for (int im=0;im<ncomp_sum_int_user1;im++) {
  NS_sumdata_sweep[IQ_USER_SUM_COMP+im]=0; //update 1st sweep
 }
 for (int im=0;im<ncomp_sum_int_user2;im++) {
   //update 2nd sweep
  NS_sumdata_sweep[IQ_USER_SUM_COMP+ncomp_sum_int_user1+im]=1; 
 }

 NS_sumdata_type[IQ_VORT_ERROR_SUM_COMP]=3;  // reduce real max (-1.0E+6)
 NS_sumdata_type[IQ_TEMP_ERROR_SUM_COMP]=3;  // reduce real max (-1.0E+6)
 NS_sumdata_type[IQ_VEL_ERROR_SUM_COMP]=3;   // reduce real max (-1.0E+6)

 for (int idir=0;idir<3;idir++) {
  for (int im=0;im<num_materials;im++) {
   NS_sumdata_type[idir+IQ_MININT_SUM_COMP+3*im]=2; // reduce real min (1.0E+6)
   NS_sumdata_type[idir+IQ_MAXINT_SUM_COMP+3*im]=3; // reduce real max (-1.0E+6)
  }
 }
 for (int im=0;im<num_materials;im++) {
  NS_sumdata_type[IQ_MININT_SLICE_SUM_COMP+im]=2; // reduce real min (1.0E+6)
  NS_sumdata_type[IQ_MAXINT_SLICE_SUM_COMP+im]=3; // reduce real max (-1.0E+6)
 }
 for (int idir=0;idir<2*num_materials;idir++) {
  NS_sumdata_type[idir+IQ_MINSTATE_SUM_COMP]=2;  // reduce real min
  NS_sumdata_type[idir+IQ_MAXSTATE_SUM_COMP]=3;  // reduce real max
 }

 NS_sumdata_type[IQ_XNOT_AMP_SUM_COMP]=3;  // x=0 amplitude  material 1

 for (int idir=0;idir<num_materials;idir++) {
  NS_sumdata_type[idir+IQ_MINCEN_SUM_COMP]=2;  // min dist from centroid
  NS_sumdata_type[idir+IQ_MAXCEN_SUM_COMP]=3;  // max dist from centroid
  NS_sumdata_sweep[idir+IQ_MINCEN_SUM_COMP]=1;  
  NS_sumdata_sweep[idir+IQ_MAXCEN_SUM_COMP]=1; 
 }

 for (int isum=0;isum<IQ_TOTAL_SUM_COMP;isum++) {
  NS_sumdata[isum]=0.0;
  if (NS_sumdata_type[isum]==2) // min
   NS_sumdata[isum]=1.0E+6;
  else if (NS_sumdata_type[isum]==3)  // max
   NS_sumdata[isum]=-1.0E+6;
  else if (NS_sumdata_type[isum]==1)
   NS_sumdata[isum]=0.0;
  else
   amrex::Error("sumdata_type invalid");
 } // isum

} // end subroutine setup_integrated_quantities

// read_params is called from:
// NavierStokes::variableSetUp
void
NavierStokes::read_params ()
{
    ParmParse pp("ns");

    pp.queryAdd("check_nan",check_nan);
    pp.queryAdd("v",verbose);
    pp.queryAdd("fab_verbose",fab_verbose);
    pp.queryAdd("output_drop_distribution",output_drop_distribution);
    pp.queryAdd("show_timings",show_timings);
    pp.queryAdd("show_mem",show_mem);

    pp.queryAdd("slice_dir",slice_dir);
    xslice.resize(AMREX_SPACEDIM);
    if ((slice_dir>=0)&&(slice_dir<AMREX_SPACEDIM)) {
     for (int i=0;i<AMREX_SPACEDIM;i++)
      xslice[i]=0.0;
     pp.queryAdd("xslice",xslice,AMREX_SPACEDIM);
    } else
     amrex::Error("slice_dir invalid");


    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "check_nan " << check_nan << '\n';
     std::cout << "NavierStokes.verbose " << verbose << '\n';
     std::cout << "NavierStokes.fab_verbose " << fab_verbose << '\n';
     std::cout << "slice_dir " << slice_dir << '\n';
     for (int i=0;i<AMREX_SPACEDIM;i++) {
      std::cout << "i=" << i << '\n';
      std::cout << "xslice " << xslice[i] << '\n';
     }
    } 

    pp.queryAdd("nblocks",nblocks);

    int nblocks_size=( (nblocks==0) ? 1 : nblocks );
    xblocks.resize(nblocks_size);
    yblocks.resize(nblocks_size);
    zblocks.resize(nblocks_size);
    rxblocks.resize(nblocks_size);
    ryblocks.resize(nblocks_size);
    rzblocks.resize(nblocks_size);

    if ((nblocks>0)&&(nblocks<10)) {
     pp.getarr("xblocks",xblocks,0,nblocks);
     pp.getarr("yblocks",yblocks,0,nblocks);
     pp.getarr("zblocks",zblocks,0,nblocks);
     pp.getarr("rxblocks",rxblocks,0,nblocks);
     pp.getarr("ryblocks",ryblocks,0,nblocks);
     pp.getarr("rzblocks",rzblocks,0,nblocks);
     if (ParallelDescriptor::IOProcessor()) {
      std::cout << "nblocks " << nblocks << '\n';
      for (int i=0;i<nblocks;i++) {
       std::cout << "i=" << i << '\n';
       std::cout << "xblocks " << xblocks[i] << '\n';
       std::cout << "yblocks " << yblocks[i] << '\n';
       std::cout << "zblocks " << zblocks[i] << '\n';
       std::cout << "rxblocks " << rxblocks[i] << '\n';
       std::cout << "ryblocks " << ryblocks[i] << '\n';
       std::cout << "rzblocks " << rzblocks[i] << '\n';
      }
     }
    } else if (nblocks!=0)
     amrex::Error("nblocks out of range");

    num_materials=0;
    num_interfaces=0;
    pp.get("num_materials",num_materials);
    if ((num_materials<2)||(num_materials>999))
     amrex::Error("num materials invalid");

    num_interfaces=( (num_materials-1)*(num_materials-1)+num_materials-1 )/2;
    if ((num_interfaces<1)||(num_interfaces>999))
     amrex::Error("num interfaces invalid");

    ncomp_sum_int_user1=0;
    pp.queryAdd("ncomp_sum_int_user",ncomp_sum_int_user1);
    if (ncomp_sum_int_user1==0) {
     pp.queryAdd("ncomp_sum_int_user1",ncomp_sum_int_user1);
    }
    if (ncomp_sum_int_user1>=0) {
     // do nothing
    } else
     amrex::Error("ncomp_sum_int_user1 invalid");

    pp.queryAdd("ncomp_sum_int_user2",ncomp_sum_int_user2);
    if (ncomp_sum_int_user2>=0) {
     // do nothing
    } else
     amrex::Error("ncomp_sum_int_user2 invalid");

    ncomp_sum_int_user12=ncomp_sum_int_user1+ncomp_sum_int_user2;

    BLB_MATRIX=0;
    BLB_RHS=BLB_MATRIX+3*(2*AMREX_SPACEDIM)*(2*AMREX_SPACEDIM);
    BLB_VEL=BLB_RHS+3*(2*AMREX_SPACEDIM);
    BLB_INT_MOM=BLB_VEL+3*(2*AMREX_SPACEDIM);
    BLB_ENERGY=BLB_INT_MOM+2*(2*AMREX_SPACEDIM);
    BLB_MASS_VEL=BLB_ENERGY+1;
    BLB_VOL=BLB_MASS_VEL+3;
    BLB_CEN_INT=BLB_VOL+1;
    BLB_CEN_ACT=BLB_CEN_INT+AMREX_SPACEDIM;
    BLB_PERIM=BLB_CEN_ACT+AMREX_SPACEDIM;
    BLB_PERIM_MAT=BLB_PERIM+1;
    BLB_TRIPLE_PERIM=BLB_PERIM_MAT+num_materials;
    BLB_CELL_CNT=BLB_TRIPLE_PERIM+num_materials*num_materials; //F_m>=1/2
    BLB_CELLVOL_CNT=BLB_CELL_CNT+1; //F_m>=1/2
    BLB_MASS=BLB_CELLVOL_CNT+1;
    BLB_PRES=BLB_MASS+1; //F_m>=1/2
    num_elements_blobclass=BLB_PRES+1;

    fort_blb_init(
     &BLB_MATRIX,
     &BLB_RHS,
     &BLB_VEL,
     &BLB_INT_MOM,
     &BLB_ENERGY,
     &BLB_MASS_VEL,
     &BLB_VOL,
     &BLB_CEN_INT,
     &BLB_CEN_ACT,
     &BLB_PERIM,
     &BLB_PERIM_MAT,
     &BLB_TRIPLE_PERIM,
     &BLB_CELL_CNT,
     &BLB_CELLVOL_CNT,
     &BLB_MASS,
     &BLB_PRES,
     &num_elements_blobclass);

    int ns_max_level;
    int cnt_max_grid_size;

    ParmParse ppamr("amr");
    ppamr.get("max_level",ns_max_level);
    Vector<int> ns_n_error_buf;
    ns_n_error_buf.resize(ns_max_level);
    for (int ilev=0;ilev<ns_max_level;ilev++) 
     ns_n_error_buf[ilev]=1;
    ppamr.queryAdd("n_error_buf",ns_n_error_buf,ns_max_level);

    ns_max_grid_size.resize(ns_max_level+1);
    for (int ilev=0;ilev<=ns_max_level;ilev++) {
     if (AMREX_SPACEDIM==2) {
      ns_max_grid_size[ilev]=128;
     } else if (AMREX_SPACEDIM==3) {
      ns_max_grid_size[ilev]=32;
     } else
      amrex::Error("AMREX_SPACEDIM invalid");
    } //ilev=0..max_level

    cnt_max_grid_size=ppamr.countval("max_grid_size");
  
    if (cnt_max_grid_size==0) {
     // do nothing
    } else if (cnt_max_grid_size>0) {
     Vector<int> mgs;
     ppamr.getarr("max_grid_size",mgs);
     int last_mgs = mgs.back();
      // the newly added components to "mgs" are initialized with
      // the last component of "mgs" prior to the resizing.
     mgs.resize(ns_max_level+1,last_mgs);
     for (int ilev=0;ilev<=ns_max_level;ilev++) {
      ns_max_grid_size[ilev]=mgs[ilev];

      if (AMREX_SPACEDIM==2) {
       if (ns_max_grid_size[ilev]<64)
        amrex::Error("expecting max_grid_size>=64 in 2D or 3DRZ");
      } else if (AMREX_SPACEDIM==3) {
       if (ns_max_grid_size[ilev]<32)
        amrex::Error("expecting max_grid_size>=32 in 3D");
      } else
       amrex::Error("AMREX_SPACEDIM invalid");

     }
    } else
     amrex::Error("cnt_max_grid_size invalid");

    int def_n_proper=1;
    ppamr.queryAdd("n_proper",def_n_proper);

    int local_plotfile_on_restart=0;
    ppamr.queryAdd("plotfile_on_restart",local_plotfile_on_restart);
    int local_checkpoint_on_restart=0;
    ppamr.queryAdd("checkpoint_on_restart",local_checkpoint_on_restart);

    tecplot_max_level=ns_max_level;
    max_level_for_use=ns_max_level;
    pp.queryAdd("tecplot_max_level",tecplot_max_level);
    pp.queryAdd("max_level_for_use",max_level_for_use);

    int max_level_two_materials=999;
    pp.queryAdd("max_level_two_materials",max_level_two_materials);
    if (max_level_two_materials==999) {
     // do nothing
    } else {
     std::cout << "max_level_two_materials no longer an option\n";
     amrex::Error("aborting ...");
    }

    Vector<int> radius_cutoff;
    radius_cutoff.resize(num_materials);
    for (int i=0;i<num_materials;i++)
     radius_cutoff[i]=999;
    pp.queryAdd("radius_cutoff",radius_cutoff,num_materials);
    for (int i=0;i<num_materials;i++) {
     if (radius_cutoff[i]==999) {
      // do nothing
     } else {
      std::cout << "radius_cutoff no longer an option\n";
      amrex::Error("aborting ...");
     }
    } // i=0..num_materials-1

    if ((tecplot_max_level<0)||
        (tecplot_max_level>ns_max_level))
     amrex::Error("tecplot_max_level invalid"); 

    if ((max_level_for_use<0)||
        (max_level_for_use>ns_max_level))
     amrex::Error("max_level_for_use invalid"); 

    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "ncomp_sum_int_user1 " << 
       ncomp_sum_int_user1 << '\n';
     std::cout << "ncomp_sum_int_user2 " << 
       ncomp_sum_int_user2 << '\n';
     std::cout << "ncomp_sum_int_user12 " << 
       ncomp_sum_int_user12 << '\n';
     std::cout << "tecplot_max_level " << 
       tecplot_max_level << '\n';
     std::cout << "max_level_for_use " << 
       max_level_for_use << '\n';
    }
    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "num_elements_blobclass= " << 
      num_elements_blobclass << '\n';
    }

    pp.queryAdd("ncoarseblocks",ncoarseblocks);
    int ncoarseblocks_size=( (ncoarseblocks==0) ? 1 : ncoarseblocks );
    xcoarseblocks.resize(ncoarseblocks_size);
    ycoarseblocks.resize(ncoarseblocks_size);
    zcoarseblocks.resize(ncoarseblocks_size);
    rxcoarseblocks.resize(ncoarseblocks_size);
    rycoarseblocks.resize(ncoarseblocks_size);
    rzcoarseblocks.resize(ncoarseblocks_size);
    if ((ncoarseblocks>0)&&(ncoarseblocks<10)) {
     pp.getarr("xcoarseblocks",xcoarseblocks,0,ncoarseblocks);
     pp.getarr("ycoarseblocks",ycoarseblocks,0,ncoarseblocks);
     pp.getarr("zcoarseblocks",zcoarseblocks,0,ncoarseblocks);
     pp.getarr("rxcoarseblocks",rxcoarseblocks,0,ncoarseblocks);
     pp.getarr("rycoarseblocks",rycoarseblocks,0,ncoarseblocks);
     pp.getarr("rzcoarseblocks",rzcoarseblocks,0,ncoarseblocks);

     if (ParallelDescriptor::IOProcessor()) {

      std::cout << "ncoarseblocks " << ncoarseblocks << '\n';
      for (int i=0;i<ncoarseblocks;i++) {
       std::cout << "i=" << i << '\n';
       std::cout << "xcoarseblocks " << xcoarseblocks[i] << '\n';
       std::cout << "ycoarseblocks " << ycoarseblocks[i] << '\n';
       std::cout << "zcoarseblocks " << zcoarseblocks[i] << '\n';
       std::cout << "rxcoarseblocks " << rxcoarseblocks[i] << '\n';
       std::cout << "rycoarseblocks " << rycoarseblocks[i] << '\n';
       std::cout << "rzcoarseblocks " << rzcoarseblocks[i] << '\n';
      }
     }
    } else if (ncoarseblocks!=0)
     amrex::Error("ncoarseblocks out of range");
     
    Vector<int> lo_bc(AMREX_SPACEDIM);
    Vector<int> hi_bc(AMREX_SPACEDIM);
    Vector<int> viscosity_lo_bc(AMREX_SPACEDIM);
    Vector<int> viscosity_hi_bc(AMREX_SPACEDIM);
    Vector<int> temperature_lo_bc(AMREX_SPACEDIM);
    Vector<int> temperature_hi_bc(AMREX_SPACEDIM);
    Vector<int> species_lo_bc(AMREX_SPACEDIM);
    Vector<int> species_hi_bc(AMREX_SPACEDIM);
    pp.getarr("lo_bc",lo_bc,0,AMREX_SPACEDIM);
    pp.getarr("hi_bc",hi_bc,0,AMREX_SPACEDIM);
    for (int i = 0; i < AMREX_SPACEDIM; i++) {
     phys_bc.setLo(i,lo_bc[i]);
     phys_bc.setHi(i,hi_bc[i]);
     viscosity_phys_bc.setLo(i,lo_bc[i]);
     viscosity_phys_bc.setHi(i,hi_bc[i]);
     viscosity_lo_bc[i]=lo_bc[i];
     viscosity_hi_bc[i]=hi_bc[i];
     temperature_phys_bc.setLo(i,lo_bc[i]);
     temperature_phys_bc.setHi(i,hi_bc[i]);
     temperature_lo_bc[i]=lo_bc[i];
     temperature_hi_bc[i]=hi_bc[i];
     species_phys_bc.setLo(i,lo_bc[i]);
     species_phys_bc.setHi(i,hi_bc[i]);
     species_lo_bc[i]=lo_bc[i];
     species_hi_bc[i]=hi_bc[i];
    }
    pp.queryAdd("viscosity_lo_bc",viscosity_lo_bc,AMREX_SPACEDIM);
    pp.queryAdd("viscosity_hi_bc",viscosity_hi_bc,AMREX_SPACEDIM);
    pp.queryAdd("temperature_lo_bc",temperature_lo_bc,AMREX_SPACEDIM);
    pp.queryAdd("temperature_hi_bc",temperature_hi_bc,AMREX_SPACEDIM);
    pp.queryAdd("species_lo_bc",species_lo_bc,AMREX_SPACEDIM);
    pp.queryAdd("species_hi_bc",species_hi_bc,AMREX_SPACEDIM);
    for (int i = 0; i < AMREX_SPACEDIM; i++) {
     viscosity_phys_bc.setLo(i,viscosity_lo_bc[i]);
     viscosity_phys_bc.setHi(i,viscosity_hi_bc[i]);
     temperature_phys_bc.setLo(i,temperature_lo_bc[i]);
     temperature_phys_bc.setHi(i,temperature_hi_bc[i]);
     species_phys_bc.setLo(i,species_lo_bc[i]);
     species_phys_bc.setHi(i,species_hi_bc[i]);
    }

     // call after phys_bc initialized since phys_bc might have to be
     // modified in this routine.
    read_geometry();

    int geometry_coord;
    Vector<Real> geometry_prob_lo;
    Vector<Real> geometry_prob_hi;
    Vector<int> geometry_is_periodic;
    int geometry_is_any_periodic;
    read_geometry_raw(geometry_coord,geometry_prob_lo,geometry_prob_hi,
		 geometry_is_periodic,geometry_is_any_periodic);
    NS_geometry_coord=geometry_coord;

    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "NS_geometry_coord= " << NS_geometry_coord << '\n';
    }

    int geometry_is_all_periodic=1;

    if (geometry_is_any_periodic==1) {
     for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
      if (geometry_is_periodic[dir]==1) {
       if ((lo_bc[dir] != Interior)||
           (viscosity_lo_bc[dir] != Interior)||
           (species_lo_bc[dir] != Interior)||
           (temperature_lo_bc[dir]!=Interior)) {
        std::cerr << "NavierStokes::variableSetUp:periodic in direction "
            << dir << " but low BC is not Interior\n";
        amrex::Abort("NavierStokes::read_params()");
       }
       if ((hi_bc[dir] != Interior)||
           (viscosity_hi_bc[dir] != Interior)||
           (species_hi_bc[dir] != Interior)||
           (temperature_hi_bc[dir]!=Interior)) {
        std::cerr << "NavierStokes::variableSetUp:periodic in direction "
            << dir << " but high BC is not Interior\n";
        amrex::Abort("NavierStokes::read_params()");
       }
      } else if (geometry_is_periodic[dir]==0) {
       geometry_is_all_periodic=0;
      } else
       amrex::Error("geometry_is_periodic[dir] invalid"); 
     } // dir=0..sdim-1
    } else if (geometry_is_any_periodic==0) {
     geometry_is_all_periodic=0;
    } else
     amrex::Error("geometry_is_any_periodic invalid");

    fort_set_periodic_var(geometry_is_periodic.dataPtr());

    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
     if (geometry_is_periodic[dir]==0) {
      if ((lo_bc[dir] == Interior)||
          (viscosity_lo_bc[dir] == Interior)||
          (species_lo_bc[dir] == Interior)||
          (temperature_lo_bc[dir]==Interior)) {
       std::cerr << "NavierStokes::variableSetUp:Interior bc in direction "
                 << dir << " but not defined as periodic\n";
       amrex::Abort("NavierStokes::read_params()");
      }
      if ((hi_bc[dir] == Interior)||
          (viscosity_hi_bc[dir] == Interior)||
          (species_hi_bc[dir] == Interior)||
          (temperature_hi_bc[dir]==Interior)) {
       std::cerr << "NavierStokes::variableSetUp:Interior bc in direction "
                 << dir << " but not defined as periodic\n";
       amrex::Abort("NavierStokes::read_params()");
      }
     }
    } // dir=0.. sdim-1

    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "phys_bc= " << phys_bc << '\n';
     std::cout << "viscosity_phys_bc= " << viscosity_phys_bc << '\n';
     std::cout << "temperature_phys_bc= " << temperature_phys_bc << '\n';
     std::cout << "species_phys_bc= " << species_phys_bc << '\n';
    }

    Real problen_min=geometry_prob_hi[0]-geometry_prob_lo[0];
    if (problen_min>0.0) {
     for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
      Real problen=geometry_prob_hi[dir]-geometry_prob_lo[dir];
      if (problen>0.0) {
       if (problen<problen_min) 
        problen_min=problen;
      } else
       amrex::Error("problen invalid");
     }
    } else
     amrex::Error("problen_min invalid");

    Real microlayer_size_default=problen_min/1.0e+9;

    //
    // Get timestepping parameters.
    //
    pp.get("cfl",cfl);

    if ((cfl>0.0)&&(cfl<=0.95)) {
     // do nothing
    } else {
     std::cout << "cfl out of range (0<cfl<=0.95); cfl=" << cfl << '\n';
     amrex::Error("default value for cfl is 1/2");
    }

    pp.queryAdd("enable_spectral",enable_spectral);

    if (enable_spectral==1) {

     if (geometry_is_all_periodic==1) {
      // do nothing
     } else if (geometry_is_all_periodic==0) {
      amrex::Warning("no slip BC not implemented for space order>2");
     } else
      amrex::Error("geometry_is_all_periodic invalid");

    } else if (enable_spectral==0) {
     // do nothing
    } else
     amrex::Error("enable_spectral invalid");

    pp.queryAdd("continuous_mof",continuous_mof);
    if (continuous_mof==2) {
     amrex::Error("continuous_mof==2 is an anachronism, set to 1");
    } else if (continuous_mof>=STANDARD_MOF) {
     //do nothing
    } else
     amrex::Error("continuous_mof invalid");

    pp.queryAdd("update_centroid_after_recon",update_centroid_after_recon);
    if (update_centroid_after_recon==0) {
     if (continuous_mof>CMOF_X)
      amrex::Error("expecting update_centroid_after_recon=1");
    } else if (update_centroid_after_recon==1) {
     if (continuous_mof==STANDARD_MOF)
      amrex::Error("expecting update_centroid_after_recon=0");
    } else
     amrex::Error("expecting update_centroid_after_recon=0 or 1");

    pp.queryAdd("mof_machine_learning",mof_machine_learning);
    pp.queryAdd("mof_decision_tree_learning",mof_decision_tree_learning);

    centroid_noise_factor.resize(num_materials);
    for (int i=0;i<num_materials;i++) {
     centroid_noise_factor[i]=0.0;
    }
    pp.queryAdd("centroid_noise_factor",centroid_noise_factor,num_materials);

    pp.queryAdd("partial_cmof_stencil_at_walls",
		partial_cmof_stencil_at_walls);

    pp.queryAdd("init_shrink",init_shrink);
    if ((init_shrink>0.0)&&(init_shrink<=1.0)) {
     //do nothing
    } else
     amrex::Error("need to have: 0.0<ns.init_shrink<=1.0");

    pp.queryAdd("dt_max",dt_max);

    pp.queryAdd("change_max",change_max);
    if ((change_max>=1.0)&&(change_max<1.011)) {
     //do nothing
    } else {
     std::cout << "ns.change_max now invalid: " << change_max << '\n';
     amrex::Error("need to have 1.0<=ns.change_max<=1.01");
    }

    change_max_init=change_max;
    pp.queryAdd("change_max_init",change_max_init);

    if ((change_max_init>=1.0)&&(change_max_init<1.011)) {
     //do nothing
    } else if (change_max_init<1.0) {
     amrex::Error("need to have 1.0<=ns.change_max_init<=1.01");
    } else if (change_max_init>=1.011) {
     if ((init_shrink>0.0)&&(init_shrink<1.0)) {
      //do nothing
     } else if (init_shrink>=1.0) {
      amrex::Error("cannot have change_max_init>1.01 and init_shrink>=1.0");
     } else
      amrex::Error("init_shrink invalid");
    } else
     amrex::Error("ns.change_max_init invalid");

    pp.queryAdd("fixed_dt",fixed_dt);
    fixed_dt_init=fixed_dt;
    pp.queryAdd("fixed_dt_init",fixed_dt_init);

    pp.queryAdd("min_velocity_for_dt",min_velocity_for_dt);
    if (min_velocity_for_dt<0.0)
     amrex::Error("min_velocity_for_dt invalid");

    pp.queryAdd("min_stefan_velocity_for_dt",min_stefan_velocity_for_dt);
    if (min_stefan_velocity_for_dt<0.0)
     amrex::Error("min_stefan_velocity_for_dt invalid");

    pp.queryAdd("fixed_dt_velocity",fixed_dt_velocity);
    pp.queryAdd("sum_interval",sum_interval);

    pp.queryAdd("profile_debug",profile_debug);
    if ((profile_debug!=0)&&(profile_debug!=1))
     amrex::Error("profile_debug invalid");

    pp.queryAdd("ns_tiling",ns_tiling);
    if (ns_tiling==true) {
     // do nothing
    } else if (ns_tiling==false) {
     // do nothing
    } else
     amrex::Error("ns_tiling invalid");

    if (thread_class::nthreads<1)
     amrex::Error("thread_class::nthreads invalid NavierStokes::read_params");

    FSI_touch_flag.resize(thread_class::nthreads);
    for (int tid=0;tid<thread_class::nthreads;tid++) {
     FSI_touch_flag[tid]=0;
    }

    gravity_vector.resize(AMREX_SPACEDIM);

    bool gravity_in_table=pp.contains("gravity");
    bool gravity_dir_in_table=pp.contains("gravity_dir");
    bool invert_gravity_in_table=pp.contains("invert_gravity");
    bool gravity_vector_in_table=pp.contains("gravity_vector");

    if (gravity_vector_in_table==true) {

     if ((gravity_in_table==false)&&
         (gravity_dir_in_table==false)&&
         (invert_gravity_in_table==false)) {
      // do nothing
     } else
      amrex::Error("gravity parm conflict");

     pp.getarr("gravity_vector",gravity_vector);

    } else if (gravity_vector_in_table==false) {

     pp.query("gravity",gravity);
     pp.query("gravity_dir",gravity_dir);
     pp.query("invert_gravity",invert_gravity);
     if ((gravity_dir<1)||(gravity_dir>AMREX_SPACEDIM))
      amrex::Error("gravity dir invalid");

     for (int dir=0;dir<AMREX_SPACEDIM;dir++)
      gravity_vector[dir]=0.0;

     if (invert_gravity==0) {
      gravity_vector[gravity_dir-1]=-std::abs(gravity);
     } else if (invert_gravity==1) {
      gravity_vector[gravity_dir-1]=std::abs(gravity);
     } else
      amrex::Error("invert_gravity invalid");

    } else
     amrex::Error("gravity_vector_in_table invalid");

    Real gravity_reference_wavelen_default=0.0;

    int gravity_max_index=0;
    fort_derive_gravity_dir(gravity_vector.dataPtr(),&gravity_max_index);

    if ((gravity_max_index>=1)&&
        (gravity_max_index<=AMREX_SPACEDIM)) {

     for (int local_dir=0;local_dir<AMREX_SPACEDIM;local_dir++) {
      if (local_dir+1!=gravity_max_index) {
       gravity_reference_wavelen_default=
        gravity_reference_wavelen_default+
         (geometry_prob_hi[local_dir]-
	  geometry_prob_lo[local_dir])*
         (geometry_prob_hi[local_dir]-
	  geometry_prob_lo[local_dir]);
      }
     } //local_dir=0 ... sdim-1
     gravity_reference_wavelen_default= 
	  std::sqrt(gravity_reference_wavelen_default);

    } else
     amrex::Error("gravity_max_index invalid; NavierStokes::read_params() ");

    if (gravity_reference_wavelen_default>0.0) {
     gravity_reference_wavelen=gravity_reference_wavelen_default;
     pp.queryAdd("gravity_reference_wavelen",gravity_reference_wavelen);
     if ((gravity_reference_wavelen>0.0)&&
         (gravity_reference_wavelen<=
	  gravity_reference_wavelen_default*(1.0001))) {
      // do nothing
     } else
      amrex::Error("gravity_reference_wavelen out of range: read_params()");
    } else
     amrex::Error("gravity_reference_wavelen_default invalid: read_params");

    pp.queryAdd("incremental_gravity_flag",incremental_gravity_flag);
    if (incremental_gravity_flag==1) {
     if (enable_spectral==0) {
      // do nothing
     } else if (enable_spectral==1) {
      //do nothing
     } else
      amrex::Error("enable_spectral invalid");
    } else if (incremental_gravity_flag==0) {
     // do nothing
    } else
     amrex::Error("incremental_gravity_flag invalid");

    pp.queryAdd("segregated_gravity_flag",segregated_gravity_flag);
    if (segregated_gravity_flag==1) {
     if (enable_spectral==0) {
      // do nothing
     } else if (enable_spectral==1) {
      //do nothing
      //having just SOLVETYPE_PRESGRAVITY low order will not
      //adversely effect the numerical dissipation of an 
      //otherwise spectrally accurate method (only the order
      //will be adversely effected)
     } else
      amrex::Error("enable_spectral invalid");
    } else if (segregated_gravity_flag==0) {
     // do nothing
    } else
     amrex::Error("segregated_gravity_flag invalid");

    pp.get("visc_coef",visc_coef);

    pp.queryAdd("include_viscous_heating",include_viscous_heating);
    if ((include_viscous_heating<0)||
        (include_viscous_heating>1))
     amrex::Error("include_viscous_heating invalid");
   
    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "local_plotfile_on_restart (NS) = " << 
	     local_plotfile_on_restart << '\n';
     std::cout << "local_checkpoint_on_restart (NS) = " << 
	     local_checkpoint_on_restart << '\n';
     std::cout << "def_n_proper (NS) = " << def_n_proper << '\n';
     std::cout << "profile_debug= " << profile_debug << '\n';
     std::cout << "ns_tiling= " << ns_tiling << '\n';
     std::cout << "thread_class::nthreads= " << 
       thread_class::nthreads << '\n';
     std::cout << "sum_interval " << sum_interval << '\n';
     std::cout << "show_timings " << show_timings << '\n';
     std::cout << "show_mem " << show_mem << '\n';
     std::cout << "output_drop_distribution " << 
      output_drop_distribution << '\n';
     std::cout << "visc_coef " << visc_coef << '\n';
     std::cout << "include_viscous_heating " << include_viscous_heating << '\n';

     std::cout << "gravity_vector 0..sdim-1: " << 
       gravity_vector[0] << ' ' <<
       gravity_vector[1] << ' ' <<
       gravity_vector[AMREX_SPACEDIM-1] << '\n';

     std::cout << "gravity_reference_wavelen " << 
	  gravity_reference_wavelen << '\n';

     std::cout << "incremental_gravity_flag " << 
	  incremental_gravity_flag << '\n';
     std::cout << "segregated_gravity_flag " << 
	  segregated_gravity_flag << '\n';

     std::cout << "cfl " << cfl << '\n';
     std::cout << "enable_spectral " << enable_spectral << '\n';
     std::cout << "continuous_mof " << continuous_mof << '\n';
     std::cout << "update_centroid_after_recon " << 
	    update_centroid_after_recon << '\n';
     std::cout << "mof_machine_learning " << mof_machine_learning << '\n';
     std::cout << "mof_decision_tree_learning " << 
	     mof_decision_tree_learning << '\n';
     std::cout << "partial_cmof_stencil_at_walls " << 
	    partial_cmof_stencil_at_walls << '\n';
     for (int i=0;i<num_materials;i++) {
      std::cout << "i= " << i << " centroid_noise_factor[i]=" <<
       centroid_noise_factor[i] << '\n';
     }
     
    }

    pp.queryAdd("FD_curv_interp",FD_curv_interp);
    if ((FD_curv_interp!=0)&&
        (FD_curv_interp!=1))
     amrex::Error("FD_curv_interp invalid");

    pp.queryAdd("vof_height_function",vof_height_function);
    if ((vof_height_function!=0)&&
        (vof_height_function!=1))
     amrex::Error("vof_height_function invalid");

    custom_nucleation_model=0;
    pp.queryAdd("custom_nucleation_model",custom_nucleation_model);
    if ((custom_nucleation_model!=0)&&
        (custom_nucleation_model!=1))
     amrex::Error("custom_nucleation_model invalid");

    n_sites=0;
    pp.queryAdd("n_sites",n_sites);
    pos_sites.resize(4);
    for (int dir=0;dir<4;dir++)
     pos_sites[dir]=0.0;
    if (n_sites>0) {
     pos_sites.resize(4*n_sites);
     pp.getarr("pos_sites",pos_sites,0,4*n_sites);
    } else if (n_sites==0) {
     // do nothing
    } else {
     std::cout << "n_sites= " << n_sites << '\n';
     amrex::Error("n_sites invalid(2)");
    }
   
    pp.queryAdd("pos_sites_random_flag",pos_sites_random_flag);

    pp.get("denfact",denfact);
    pp.get("velfact",velfact);
    pp.get("xblob",xblob);
    pp.get("yblob",yblob);
    pp.get("zblob",zblob);
    pp.get("radblob",radblob);
    pp.get("probtype",probtype);
    pp.get("axis_dir",axis_dir);

    visual_ncell.resize(AMREX_SPACEDIM);
    for (int dir=0;dir<AMREX_SPACEDIM;dir++)
     visual_ncell[dir]=8;
    pp.queryAdd("visual_ncell",visual_ncell,AMREX_SPACEDIM);
    pp.queryAdd("visual_compare",visual_compare);

    pp.queryAdd("step_through_data",step_through_data);

    pp.queryAdd("visual_tessellate_vfrac",visual_tessellate_vfrac);
    pp.queryAdd("visual_revolve",visual_revolve);
    pp.queryAdd("visual_output_raw_State_Type",visual_output_raw_State_Type);
    pp.queryAdd("visual_output_raw_mac_Type",visual_output_raw_mac_Type);
    pp.queryAdd("visual_phase_change_plot_int",visual_phase_change_plot_int);
    pp.queryAdd("visual_buoyancy_plot_int",visual_buoyancy_plot_int);
    pp.queryAdd("visual_divergence_plot_int",visual_divergence_plot_int);
    pp.queryAdd("visual_WALLVEL_plot_int",visual_WALLVEL_plot_int);
    pp.queryAdd("visual_drag_plot_int",visual_drag_plot_int);
    pp.queryAdd("visual_nddata_format",visual_nddata_format);

    if ((visual_tessellate_vfrac!=0)&&
        (visual_tessellate_vfrac!=1)&&
	(visual_tessellate_vfrac!=3))
     amrex::Error("visual_tessellate_vfrac invalid");

    if (visual_revolve<0) {
     std::cout << "visual_revolve: " << visual_revolve << '\n';
     amrex::Error("visual_revolve invalid");
    }

    if (AMREX_SPACEDIM==2) {
     adapt_quad_depth=2;
     if (is_zalesak()==1) {
      adapt_quad_depth=5; 
     }
    } else if (AMREX_SPACEDIM==3) {
     adapt_quad_depth=1;
     if (is_zalesak()==1) {
      adapt_quad_depth=3;
     }
    } else
     amrex::Error("dimension bust");

    pp.queryAdd("adapt_quad_depth",adapt_quad_depth);

    law_of_the_wall.resize(num_materials);
    wall_model_velocity.resize(num_materials);
    interface_mass_transfer_model.resize(2*num_interfaces);

    for (int i=0;i<2*num_interfaces;i++) {
     interface_mass_transfer_model[i]=0;
    }

    for (int i=0;i<num_materials;i++) {
     law_of_the_wall[i]=0;
     wall_model_velocity[i]=0.0;
    }
    pp.queryAdd("law_of_the_wall",law_of_the_wall,num_materials);
    pp.queryAdd("wall_model_velocity",wall_model_velocity,num_materials);
    pp.queryAdd("interface_mass_transfer_model",
       interface_mass_transfer_model,2*num_interfaces);

    for (int i=0;i<num_materials;i++) {
     if ((law_of_the_wall[i]==0)||
         (law_of_the_wall[i]==1)||
  	 (law_of_the_wall[i]==2)) {
      // do nothing
     } else {
      amrex::Error("law_of_the_wall invalid");
     }
    }

    pp.queryAdd("wall_slip_weight",wall_slip_weight);
    if ((wall_slip_weight>=0.0)&&
        (wall_slip_weight<=1.0)) {
     // do nothing
    } else {
     amrex::Error("wall_slip_weight invalid");
    }

    pp.queryAdd("ZEYU_DCA_SELECT",ZEYU_DCA_SELECT);
    if ((ZEYU_DCA_SELECT==-1)|| //static
        ((ZEYU_DCA_SELECT>=1)&&
	 (ZEYU_DCA_SELECT<=8))) {
     // do nothing
    } else
     amrex::Error("ZEYU_DCA_SELECT invalid");

    FSI_flag.resize(num_materials);
    FSI_refine_factor.resize(num_materials);
    FSI_bounding_box_ngrow.resize(num_materials);

    material_type.resize(num_materials);
    material_type_interface.resize(num_interfaces);

    pp.getarr("material_type",material_type,0,num_materials);
    material_type_evap.resize(num_materials);
    material_type_lowmach.resize(num_materials);
    material_type_visual.resize(num_materials);
    for (int i=0;i<num_materials;i++) {
     material_type_evap[i]=material_type[i];
     material_type_lowmach[i]=material_type[i];
     material_type_visual[i]=material_type[i];
    }
    pp.queryAdd("material_type_evap",material_type_evap,num_materials);
    pp.queryAdd("material_type_lowmach",material_type_lowmach,num_materials);
    for (int i=0;i<num_materials;i++) {
     material_type_visual[i]=material_type_lowmach[i];
    }
    pp.queryAdd("material_type_visual",material_type_visual,num_materials);
 
    for (int i=0;i<num_materials;i++) {
     FSI_flag[i]=FSI_FLUID;
     FSI_refine_factor[i]=1;
     FSI_bounding_box_ngrow[i]=3;
    }
    pp.queryAdd("FSI_flag",FSI_flag,num_materials);
    pp.queryAdd("num_local_aux_grids",num_local_aux_grids);
    pp.queryAdd("FSI_interval",FSI_interval);
    pp.queryAdd("FSI_refine_factor",FSI_refine_factor,num_materials);
    pp.queryAdd("FSI_bounding_box_ngrow",FSI_bounding_box_ngrow,num_materials);

    for (int i=0;i<num_materials;i++) {
     if (FSI_flag[i]==FSI_SHOELE_CTML) {
      if (FSI_interval==1) {
       // do nothing
      } else
       amrex::Error("Eulerian data must be regenerated every step");
     } else if (FSI_flag_valid(i)==1) {
      // do nothing
     } else {
      amrex::Error("FSI_flag_valid(i)!=1");
     }
    } // i=0..num_materials-1

    int nparts=0;
    int CTML_FSI_numsolids_test=0;
    for (int im=0;im<num_materials;im++) {

     if (FSI_flag[im]==FSI_SHOELE_CTML) { 
      CTML_FSI_numsolids_test++;
     } else if (FSI_flag_valid(im)==1) {
      //do nothing
     } else
      amrex::Error("FSI_flag_valid invalid");
 
     if (ns_is_lag_part(im)==1) {
      nparts++;
     } else if (ns_is_lag_part(im)==0) {
      // do nothing
     } else
      amrex::Error("ns_is_lag_part invalid");
    }  // im=0..num_materials-1
    im_solid_map.resize(nparts);

    if (CTML_FSI_numsolids!=CTML_FSI_numsolids_test)
     amrex::Error("CTML_FSI_numsolids!=CTML_FSI_numsolids_test");

    nparts=0;
    for (int im=0;im<num_materials;im++) {
     if (ns_is_lag_part(im)==1) {
      im_solid_map[nparts]=im;
      nparts++;
     } else if (ns_is_lag_part(im)==0) {
      // do nothing
     } else
      amrex::Error("ns_is_lag_part invalid");
    }  // im=0..num_materials-1
    if (nparts!=im_solid_map.size())
     amrex::Error("nparts!=im_solid_map.size()");

    if (FSI_material_exists_CTML()==1) {
     extend_pressure_into_solid=1;
    } else if (FSI_material_exists_CTML()==0) {
     // do nothing
    } else
     amrex::Error("FSI_material_exists_CTML() invalid");

    pp.queryAdd("extend_pressure_into_solid",extend_pressure_into_solid);

    if (FSI_material_exists_CTML()==1) {
     if (extend_pressure_into_solid!=1)
      amrex::Error("need extend_pressure_into_solid==1 if CTML coupling");
    } else if (FSI_material_exists_CTML()==0) {
     // do nothing
    } else
     amrex::Error("FSI_material_exists_CTML() invalid");

    elastic_viscosity.resize(num_materials);
    static_damping_coefficient.resize(num_materials);
    store_elastic_data.resize(num_materials);

    for (int im=0;im<num_materials;im++) {
     elastic_viscosity[im]=0.0;
     static_damping_coefficient[im]=0.0;
     store_elastic_data[im]=0;
    }
    pp.queryAdd("elastic_viscosity",elastic_viscosity,num_materials);

    viscoelastic_model.resize(num_materials);
    for (int i=0;i<num_materials;i++)
     viscoelastic_model[i]=0;
    pp.queryAdd("viscoelastic_model",viscoelastic_model,num_materials);

    for (int i=0;i<num_materials;i++) {
     if (fort_built_in_elastic_model(&elastic_viscosity[i],
			           &viscoelastic_model[i])==1) {
      // do nothing
     } else if (fort_built_in_elastic_model(&elastic_viscosity[i],
                                          &viscoelastic_model[i])==0) {
      // do nothing
     } else
      amrex::Error("fort_built_in_elastic_model invalid");
    } // i=0..num_materials-1

    pp.queryAdd("static_damping_coefficient",
		static_damping_coefficient,num_materials);

    for (int im=0;im<num_materials;im++) {
     if (fort_built_in_elastic_model(&elastic_viscosity[im],
       &viscoelastic_model[im])==1) {
      store_elastic_data[im]=1;
     } else if (fort_built_in_elastic_model(&elastic_viscosity[im],
                 &viscoelastic_model[im])==0) {
      // do nothing
     } else
      amrex::Error("fort_built_in_elastic_model invalid");
    } // im=0..num_materials-1 

    num_materials_viscoelastic=0;
    for (int im=0;im<num_materials;im++) {
     if (store_elastic_data[im]==1) {
      num_materials_viscoelastic++;
     } else if (store_elastic_data[im]==0) {
      // do nothing
     } else
      amrex::Error("store_elastic_data invalid");
    } // im=0..num_materials-1 

    im_elastic_map.resize(num_materials_viscoelastic);

    int elastic_partid=0;
    for (int i=0;i<num_materials;i++) {
     if (store_elastic_data[i]==1) {
      im_elastic_map[elastic_partid]=i;
      elastic_partid++;
     } else if (store_elastic_data[i]==0) {
      // do nothing
     } else
      amrex::Error("store_elastic_data invalid");
    } // im=0..num_materials-1 

    if (elastic_partid==num_materials_viscoelastic) {
     // do nothing
    } else
     amrex::Error("elastic_partid==num_materials_viscoelastic failed");

    NUM_STATE_TYPE=DIV_Type+1;

    if (nparts==0) {
     Solid_State_Type=-1;
    } else if ((nparts>=1)&&(nparts<=num_materials)) {
     Solid_State_Type=DIV_Type+1;
     NUM_STATE_TYPE++;
    } else
     amrex::Error("nparts invalid");
 
    if ((num_materials_viscoelastic>=1)&&
        (num_materials_viscoelastic<=num_materials)) {
     Tensor_Type=NUM_STATE_TYPE;

     TensorX_Type=Tensor_Type+1;
     TensorY_Type=TensorX_Type+1;
     TensorZ_Type=TensorY_Type+1;

     NUM_STATE_TYPE=TensorZ_Type+1;
    } else if (num_materials_viscoelastic==0) {
     Tensor_Type=-1;
     TensorX_Type=-1;
     TensorY_Type=-1;
     TensorZ_Type=-1;
    } else
     amrex::Error("num_materials_viscoelastic invalid");

    for (int i=0;i<num_materials;i++) {
     if (material_type[i]==0) {
      if (ns_is_rigid(i)!=0)
       amrex::Error("ns_is_rigid invalid");
     } else if (material_type[i]==999) {
      if (ns_is_rigid(i)!=1)
       amrex::Error("ns_is_rigid invalid");
     } else if ((material_type[i]>0)&&
                (material_type[i]<999)) {
      if (ns_is_rigid(i)!=0)
       amrex::Error("ns_is_rigid invalid");
     } else
      amrex::Error("material_type invalid");

    } // i=0..num_materials-1

     //smooth_type: 0=GSRB 1=weighted Jacobi 2=ILU
    ParmParse pplp("Lp");
    pplp.queryAdd("smooth_type",smooth_type);
    pplp.queryAdd("bottom_smooth_type",bottom_smooth_type);
    pplp.queryAdd("use_mg_precond_in_mglib",use_mg_precond_in_mglib);

    pplp.queryAdd("bottom_bottom_tol_factor",bottom_bottom_tol_factor);

    ParmParse ppmg("mg");
    ppmg.queryAdd("presmooth",global_presmooth);
    ppmg.queryAdd("postsmooth",global_postsmooth);
    if (global_presmooth!=global_postsmooth)
     amrex::Error("global_presmooth!=global_postsmooth");

    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "extend_pressure_into_solid " << 
      extend_pressure_into_solid << '\n';
     std::cout << "smooth_type " << smooth_type << '\n';
     std::cout << "bottom_smooth_type " << bottom_smooth_type << '\n';
     std::cout << "global_presmooth " << global_presmooth << '\n';
     std::cout << "global_postsmooth " << global_postsmooth << '\n';
     std::cout << "use_mg_precond_in_mglib " <<use_mg_precond_in_mglib<<'\n';

     std::cout << "bottom_bottom_tol_factor " <<
       bottom_bottom_tol_factor<<'\n';
     std::cout << "FSI_interval " << FSI_interval << '\n';
     std::cout << "num_local_aux_grids " << num_local_aux_grids  << '\n';
     for (int i=0;i<num_materials;i++) {
      std::cout << "i= " << i << " FSI_flag " << FSI_flag[i] << '\n';
      std::cout << "i= " << i << " FSI_refine_factor " << 
	     FSI_refine_factor[i] << '\n';
      std::cout << "i= " << i << " FSI_bounding_box_ngrow " << 
	     FSI_bounding_box_ngrow[i] << '\n';
     }

     for (int i=0;i<num_materials;i++) {
      std::cout << "law_of_the_wall i=" << i << " " << 
	      law_of_the_wall[i] << '\n';
      std::cout << "wall_model_velocity i=" << i << " " << 
	      wall_model_velocity[i] << '\n';
     }
     for (int i=0;i<2*num_interfaces;i++) {
      std::cout << "interface_mass_transfer_model i=" << i << " " << 
	      interface_mass_transfer_model[i] << '\n';
     }
     std::cout << "wall_slip_weight " << wall_slip_weight << '\n';
     std::cout << "ZEYU_DCA_SELECT " << ZEYU_DCA_SELECT << '\n';
     std::cout << "adapt_quad_depth= " << adapt_quad_depth << '\n';
    }

    pp.get("adv_dir",adv_dir);
    pp.get("adv_vel",adv_vel);
    pp.get("rgasinlet",rgasinlet);
    pp.queryAdd("slipcoeff",slipcoeff);
    pp.get("vinletgas",vinletgas);

    pp.get("twall",twall);

    pp.queryAdd("krylov_subspace_max_num_outer_iter",krylov_subspace_max_num_outer_iter);

    pp.queryAdd("EILE_flag",EILE_flag);

    if ((EILE_flag==1)|| // EILE
        (EILE_flag==2)|| // EI
        (EILE_flag==3)|| // LE
        (EILE_flag==-1)) { // Weymouth and Yue
     // do nothing
    } else
     amrex::Error("EILE_flag invalid");

    num_species_var=0;

    pp.get("num_species_var",num_species_var);
    if ((num_species_var<0)||(num_species_var>999))
     amrex::Error("num species var invalid");

    MOFITERMAX=DEFAULT_MOFITERMAX;
    pp.queryAdd("MOFITERMAX",MOFITERMAX);
    if ((MOFITERMAX<0)||(MOFITERMAX>MOFITERMAX_LIMIT)) {
     std::cout << "MOFITERMAX= " << MOFITERMAX << '\n';
     std::cout << "MOFITERMAX_LIMIT= " << MOFITERMAX_LIMIT << '\n';
     amrex::Error("mof iter max invalid in navierstokes");
    }

    MOFITERMAX_AFTER_PREDICT=DEFAULT_MOFITERMAX_AFTER_PREDICT;  
    pp.queryAdd("MOFITERMAX_AFTER_PREDICT",MOFITERMAX_AFTER_PREDICT);
    if ((MOFITERMAX_AFTER_PREDICT<0)||
	(MOFITERMAX_AFTER_PREDICT>MOFITERMAX_LIMIT)||
	(MOFITERMAX_AFTER_PREDICT>MOFITERMAX)) {
     std::cout << "MOFITERMAX= " << MOFITERMAX << '\n';
     std::cout << "MOFITERMAX_LIMIT= " << MOFITERMAX_LIMIT << '\n';
     std::cout << "MOFITERMAX_AFTER_PREDICT= " << 
       MOFITERMAX_AFTER_PREDICT << '\n';
     amrex::Error("mof iter max after predict invalid in navierstokes");
    }

    MOF_TURN_OFF_LS=0;  
    pp.queryAdd("MOF_TURN_OFF_LS",MOF_TURN_OFF_LS);
    if ((MOF_TURN_OFF_LS!=0)&&(MOF_TURN_OFF_LS!=1))
     amrex::Error("mof turn off ls invalid in navierstokes");

    MOF_DEBUG_RECON=0; 
    pp.queryAdd("MOF_DEBUG_RECON",MOF_DEBUG_RECON);
    if ((MOF_DEBUG_RECON!=0)&&(MOF_DEBUG_RECON!=1)&&
        (MOF_DEBUG_RECON!=2))
     amrex::Error("mof debug recon invalid in navierstokes");

    num_state_base=ENUM_SPECIESVAR;  // den,Temperature
    num_state_material=ENUM_SPECIESVAR;  // den,Temperature
    num_state_material+=num_species_var;

    stiffPINF.resize(num_materials);
    stiffCP.resize(num_materials);
    stiffCV.resize(num_materials);
    stiffGAMMA.resize(num_materials);

    DrhoDT.resize(num_materials);
    override_density.resize(num_materials);

    temperature_source_cen.resize(AMREX_SPACEDIM);
    temperature_source_rad.resize(AMREX_SPACEDIM);

    tempconst.resize(num_materials);
    initial_temperature.resize(num_materials);
    tempcutoff.resize(num_materials);
    tempcutoffmax.resize(num_materials);
    viscconst.resize(num_materials);
    viscconst_eddy_wall.resize(num_materials);
    viscconst_eddy_bulk.resize(num_materials);
    heatviscconst_eddy_wall.resize(num_materials);
    heatviscconst_eddy_bulk.resize(num_materials);
    viscosity_state_model.resize(num_materials);
    les_model.resize(num_materials);
    viscconst_interface.resize(num_interfaces);
    speciesreactionrate.resize((num_species_var+1)*num_materials);
    speciesconst.resize((num_species_var+1)*num_materials);
    speciesviscconst.resize((num_species_var+1)*num_materials);
    speciesviscconst_interface.resize((num_species_var+1)*num_interfaces);
    species_molar_mass.resize(num_species_var+1);

    for (int j=0;j<=num_species_var;j++)
     species_molar_mass[j]=1.0;

    for (int i=0;i<num_interfaces;i++) {
     viscconst_interface[i]=0.0;
     for (int j=0;j<=num_species_var;j++)
      speciesviscconst_interface[j*num_interfaces+i]=0.0;
    }
    pp.queryAdd("viscconst_interface",viscconst_interface,num_interfaces);
    if (num_species_var>0) {
     pp.queryAdd("speciesviscconst_interface",
      speciesviscconst_interface,num_species_var*num_interfaces);

     pp.queryAdd("species_molar_mass",
      species_molar_mass,num_species_var);
    }
     // in: read_params

    spec_material_id_LIQUID.resize(num_species_var+1);
    spec_material_id_AMBIENT.resize(num_species_var+1);
    for (int i=0;i<num_species_var+1;i++) {
     spec_material_id_LIQUID[i]=0;
     spec_material_id_AMBIENT[i]=0;
    }

    if (num_species_var>0) {
     pp.queryAdd("spec_material_id_LIQUID",spec_material_id_LIQUID,
       num_species_var);
     pp.queryAdd("spec_material_id_AMBIENT",spec_material_id_AMBIENT,
       num_species_var);
    }
    
    vorterr.resize(num_materials);
    pressure_error_cutoff.resize(num_materials);
    temperature_error_cutoff.resize(num_materials);

    microlayer_substrate.resize(num_materials);
    microlayer_angle.resize(num_materials);
    microlayer_size.resize(num_materials);
    macrolayer_size.resize(num_materials);
    max_contact_line_size.resize(num_materials);

    thermal_microlayer_size.resize(num_materials);
    shear_microlayer_size.resize(num_materials);
    buoyancy_microlayer_size.resize(num_materials);
    phasechange_microlayer_size.resize(num_materials);
 
    microlayer_temperature_substrate.resize(num_materials);

     // in: read_params
     
    cavitation_pressure.resize(num_materials);
    cavitation_vapor_density.resize(num_materials);
    cavitation_tension.resize(num_materials);
    for (int i=0;i<num_materials;i++) {
     cavitation_pressure[i]=0.0; 
     cavitation_vapor_density[i]=0.0; 
     cavitation_tension[i]=0.0; 
    }
    pp.queryAdd("cavitation_pressure",cavitation_pressure,num_materials);
    pp.queryAdd("cavitation_vapor_density",cavitation_vapor_density,
		num_materials);
    pp.queryAdd("cavitation_tension",cavitation_tension,num_materials);
 
     // in: read_params

    grid_stretching_parameter.resize(AMREX_SPACEDIM);

    hardwire_Y_gamma.resize(2*num_interfaces);
    hardwire_T_gamma.resize(2*num_interfaces);
    accommodation_coefficient.resize(2*num_interfaces);
    reference_pressure.resize(2*num_interfaces);
    saturation_temp.resize(2*num_interfaces);
    saturation_temp_curv.resize(2*num_interfaces);
    saturation_temp_vel.resize(2*num_interfaces);
    saturation_temp_min.resize(2*num_interfaces);
    saturation_temp_max.resize(2*num_interfaces);
    nucleation_temp.resize(2*num_interfaces);
    nucleation_pressure.resize(2*num_interfaces);
    nucleation_pmg.resize(2*num_interfaces);
    nucleation_mach.resize(2*num_interfaces);

    latent_heat.resize(2*num_interfaces);
    latent_heat_slope.resize(2*num_interfaces);
    latent_heat_T0.resize(2*num_interfaces);
    latent_heat_min.resize(2*num_interfaces);

    reaction_rate.resize(2*num_interfaces);
    freezing_model.resize(2*num_interfaces);
    Tanasawa_or_Schrage_or_Kassemi.resize(2*num_interfaces);
    rigid_fraction_id.resize(num_materials);
    mass_fraction_id.resize(2*num_interfaces);
    distribute_from_target.resize(2*num_interfaces);
    distribute_mdot_evenly.resize(2*num_interfaces);
    constant_volume_mdot.resize(2*num_interfaces);

    constant_density_all_time.resize(num_materials);

    tension.resize(num_interfaces);
    tension_init.resize(num_interfaces);
    tension_slope.resize(num_interfaces);
    tension_T0.resize(num_interfaces);
    tension_min.resize(num_interfaces);

    prefreeze_tension.resize(num_interfaces);

     // (dir,side)  (1,1),(2,1),(3,1),(1,2),(2,2),(3,2)
    outflow_velocity_buffer_size.resize(2*AMREX_SPACEDIM);

    cap_wave_speed.resize(num_interfaces);

    prerecalesce_stiffCP.resize(num_materials);
    prerecalesce_stiffCV.resize(num_materials);
    prerecalesce_viscconst.resize(num_materials);
    prerecalesce_heatviscconst.resize(num_materials);

    nucleation_period=0.0;
    nucleation_init_time=0.0;

    for (int i=0;i<num_materials;i++) {
     microlayer_substrate[i]=0;
     microlayer_angle[i]=0.0;
     microlayer_size[i]=0.0;
     macrolayer_size[i]=0.0;
     max_contact_line_size[i]=0.0;
     microlayer_temperature_substrate[i]=0.0;

     thermal_microlayer_size[i]=microlayer_size_default;
     shear_microlayer_size[i]=microlayer_size_default;
     buoyancy_microlayer_size[i]=microlayer_size_default;
     phasechange_microlayer_size[i]=microlayer_size_default;
    }

    for (int i=0;i<AMREX_SPACEDIM;i++) { 
     grid_stretching_parameter[i]=0.0;
    }

    for (int i=0;i<num_interfaces;i++) { 
     hardwire_Y_gamma[i]=0.0;
     hardwire_Y_gamma[i+num_interfaces]=0.0;
     hardwire_T_gamma[i]=0.0;
     hardwire_T_gamma[i+num_interfaces]=0.0;
     accommodation_coefficient[i]=0.0;
     accommodation_coefficient[i+num_interfaces]=0.0;
     reference_pressure[i]=1.0e+6;
     reference_pressure[i+num_interfaces]=1.0e+6;
     saturation_temp[i]=0.0;
     saturation_temp[i+num_interfaces]=0.0;
     saturation_temp_curv[i]=0.0;
     saturation_temp_curv[i+num_interfaces]=0.0;
     saturation_temp_vel[i]=0.0;
     saturation_temp_vel[i+num_interfaces]=0.0;
     saturation_temp_min[i]=0.0;
     saturation_temp_min[i+num_interfaces]=0.0;
     saturation_temp_max[i]=1.0e+20;
     saturation_temp_max[i+num_interfaces]=1.0e+20;
     nucleation_temp[i]=0.0;
     nucleation_temp[i+num_interfaces]=0.0;
     nucleation_pressure[i]=0.0;
     nucleation_pmg[i]=0.0;
     nucleation_mach[i]=0.0;
     nucleation_pressure[i+num_interfaces]=0.0;
     nucleation_pmg[i+num_interfaces]=0.0;
     nucleation_mach[i+num_interfaces]=0.0;

     latent_heat[i]=0.0;
     latent_heat[i+num_interfaces]=0.0;

     latent_heat_slope[i]=0.0;
     latent_heat_slope[i+num_interfaces]=0.0;

     latent_heat_T0[i]=0.0;
     latent_heat_T0[i+num_interfaces]=0.0;

     latent_heat_min[i]=0.0;
     latent_heat_min[i+num_interfaces]=0.0;

     reaction_rate[i]=0.0;
     reaction_rate[i+num_interfaces]=0.0;
     freezing_model[i]=0;
     freezing_model[i+num_interfaces]=0;
     Tanasawa_or_Schrage_or_Kassemi[i]=0;
     Tanasawa_or_Schrage_or_Kassemi[i+num_interfaces]=0;
     mass_fraction_id[i]=0;
     mass_fraction_id[i+num_interfaces]=0;
     distribute_from_target[i]=0;
     distribute_from_target[i+num_interfaces]=0;
     distribute_mdot_evenly[i]=0;
     distribute_mdot_evenly[i+num_interfaces]=0;
     constant_volume_mdot[i]=0;
     constant_volume_mdot[i+num_interfaces]=0;
    } // i=0..num_interfaces-1

    for (int i=0;i<num_materials;i++) {
     constant_density_all_time[i]=1;
     rigid_fraction_id[i]=0;
    }
    molar_mass.resize(num_materials);

    density_floor.resize(num_materials);
    for (int i=0;i<num_materials;i++) {
     density_floor[i]=0.0;
    }
    pp.queryAdd("density_floor",density_floor,num_materials);

    density_ceiling.resize(num_materials);
    for (int i=0;i<num_materials;i++) {
     density_ceiling[i]=1.0e+20;
    }

    pp.queryAdd("density_ceiling",density_ceiling,num_materials);

    denconst.resize(num_materials);
    pp.getarr("denconst",denconst,0,num_materials);

    for (int i=0;i<num_materials;i++) {

     if (density_ceiling[i]<=0.0) {
      amrex::Error("density_ceiling[i]<=0.0");
     } else if (density_ceiling[i]<=denconst[i]) {
      amrex::Error("density_ceiling[i]<=denconst[i]");
     }

     if (density_floor[i]<0.0) {
      amrex::Error("density_floor[i]<0.0");
     } else if (density_floor[i]==0.0) {
      // do nothing
     } else if (density_floor[i]>=denconst[i]) {
      amrex::Error("density_floor[i]>=denconst[i]");
     }

    } // i=0..num_materials-1

    for (int i=0;i<num_materials;i++) {
     molar_mass[i]=1.0;
    }

    pp.queryAdd("molar_mass",molar_mass,num_materials);

    denconst_interface_added.resize(num_interfaces);

    for (int iten=0;iten<num_interfaces;iten++) 
     denconst_interface_added[iten]=0.0;

    pp.queryAdd("denconst_interface_added",
      denconst_interface_added,num_interfaces);

    pp.queryAdd("stokes_flow",stokes_flow);
    pp.queryAdd("cancel_advection",cancel_advection);

    for (int i=0;i<(num_species_var+1)*num_materials;i++) {
     speciesconst[i]=0.0;
     speciesviscconst[i]=0.0;
     speciesreactionrate[i]=0.0;
    }

    for (int i=0;i<num_materials;i++) {

     stiffPINF[i]=0.0;
     stiffCP[i]=4.1855e+7;
     stiffCV[i]=4.1855e+7;
     stiffGAMMA[i]=1.4;

     tempcutoff[i]=1.0e-8;
     tempcutoffmax[i]=1.0e+99;
     DrhoDT[i]=0.0;
     override_density[i]=0;
     temperature_error_cutoff[i]=0.0;
    }

    pp.queryAdd("tempcutoff",tempcutoff,num_materials);
    pp.queryAdd("tempcutoffmax",tempcutoffmax,num_materials);

    pp.queryAdd("stiffPINF",stiffPINF,num_materials);
    pp.queryAdd("stiffCP",stiffCP,num_materials);
    for (int i=0;i<num_materials;i++)
     stiffCV[i]=stiffCP[i];

    pp.queryAdd("stiffCV",stiffCV,num_materials);
    pp.queryAdd("stiffGAMMA",stiffGAMMA,num_materials);

    pp.queryAdd("angular_velocity",angular_velocity);
    pp.queryAdd("centrifugal_force_factor",centrifugal_force_factor);

    pp.queryAdd("uncoupled_viscosity",uncoupled_viscosity);

     //Du/Dt=-grad (p-rho0 g dot z)/rho0 - g DrhoDT (T-T0) 
     //DrhoDT has units of 1/(Degrees Kelvin)
    pp.queryAdd("DrhoDT",DrhoDT,num_materials);

    pp.queryAdd("override_density",override_density,num_materials);

    pp.getarr("vorterr",vorterr,0,num_materials);
    pp.queryAdd("pressure_error_flag",pressure_error_flag);
    pp.getarr("pressure_error_cutoff",pressure_error_cutoff,0,num_materials);
    pp.queryAdd("temperature_error_cutoff",
	temperature_error_cutoff,num_materials);

    pp.queryAdd("temperature_source",temperature_source);
    pp.queryAdd("temperature_source_cen",temperature_source_cen,AMREX_SPACEDIM);
    pp.queryAdd("temperature_source_rad",temperature_source_rad,AMREX_SPACEDIM);

    pp.getarr("tempconst",tempconst,0,num_materials);
    for (int i=0;i<num_materials;i++)
     initial_temperature[i]=tempconst[i]; 
    pp.queryAdd("initial_temperature",initial_temperature,num_materials);
    pp.queryAdd("initial_temperature_diffuse_duration",
     initial_temperature_diffuse_duration);
    if (initial_temperature_diffuse_duration<0.0)
     amrex::Error("initial_temperature_diffuse_duration<0.0");

    pp.getarr("viscconst",viscconst,0,num_materials);

    viscconst_min=viscconst[0];
    viscconst_max=viscconst[0];

    for (int i=0;i<num_materials;i++) {

     if (viscconst[i]<viscconst_min)
      viscconst_min=viscconst[i];
     if (viscconst[i]>viscconst_max)
      viscconst_max=viscconst[i];

     viscconst_eddy_wall[i]=0.0;
     viscconst_eddy_bulk[i]=0.0;
     heatviscconst_eddy_wall[i]=0.0;
     heatviscconst_eddy_bulk[i]=0.0;
    }
    pp.queryAdd("viscconst_eddy_wall",viscconst_eddy_wall,num_materials);
    pp.queryAdd("viscconst_eddy_bulk",viscconst_eddy_bulk,num_materials);
    pp.queryAdd("heatviscconst_eddy_wall",
	heatviscconst_eddy_wall,num_materials);
    pp.queryAdd("heatviscconst_eddy_bulk",
	heatviscconst_eddy_bulk,num_materials);

    for (int i=0;i<num_materials;i++)
     viscosity_state_model[i]=0;
    pp.queryAdd("viscosity_state_model",viscosity_state_model,num_materials);
    heatflux_factor.resize(num_materials);
    
    for (int i=0;i<num_materials;i++) {
     heatflux_factor[i]=1.0;
     les_model[i]=0;
    }
    pp.queryAdd("les_model",les_model,num_materials);

    heatviscconst.resize(num_materials);
    heatviscconst_interface.resize(num_interfaces);
    pp.queryAdd("heatflux_factor",heatflux_factor,num_materials);
    pp.getarr("heatviscconst",heatviscconst,0,num_materials);

    heatviscconst_min=heatviscconst[0];
    heatviscconst_max=heatviscconst[0];

    for (int i=0;i<num_materials;i++) {
     if (heatviscconst[i]<heatviscconst_min)
      heatviscconst_min=heatviscconst[i];
     if (heatviscconst[i]>heatviscconst_max)
      heatviscconst_max=heatviscconst[i];
    }

    for (int i=0;i<num_interfaces;i++) {
     heatviscconst_interface[i]=0.0;
    }

    pp.queryAdd("heatviscconst_interface",heatviscconst_interface,
		num_interfaces);
    if (num_species_var>0) {
     pp.queryAdd("speciesreactionrate",speciesreactionrate,
		 num_species_var*num_materials);
     pp.queryAdd("speciesconst",speciesconst,num_species_var*num_materials);
     pp.queryAdd("speciesviscconst",speciesviscconst,
		 num_species_var*num_materials);
    }

    for (int i=0;i<num_materials;i++) {
     prerecalesce_stiffCP[i]=stiffCP[i];
     prerecalesce_stiffCV[i]=stiffCV[i];
     prerecalesce_viscconst[i]=viscconst[i];
     prerecalesce_heatviscconst[i]=heatviscconst[i];
    }
    pp.queryAdd("prerecalesce_viscconst",prerecalesce_viscconst,num_materials);
    pp.queryAdd("prerecalesce_heatviscconst",
       prerecalesce_heatviscconst,num_materials);
    pp.queryAdd("prerecalesce_stiffCP",prerecalesce_stiffCP,num_materials);
    pp.queryAdd("prerecalesce_stiffCV",prerecalesce_stiffCV,num_materials);

    pp.queryAdd("mglib_max_ratio",mglib_max_ratio);
    if (mglib_max_ratio>1.0) {
     //do nothing
    } else
     amrex::Error("mglib_max_ratio invalid");

    pp.getarr("tension",tension,0,num_interfaces);

    for (int iten=0;iten<num_interfaces;iten++)
     tension_init[iten]=tension[iten];

    pp.queryarr("tension_init",tension_init,0,num_interfaces);

    pp.queryAdd("unscaled_min_curvature_radius",unscaled_min_curvature_radius);
    if (unscaled_min_curvature_radius<2.0) {
     amrex::Error("must have unscaled_min_curvature_radius>=2.0");
    } else if (unscaled_min_curvature_radius>=2.0) {
     //do nothing
    } else
     amrex::Error("unscaled_min_curvature_radius = NaN");

    for (int i=0;i<num_interfaces;i++) 
     prefreeze_tension[i]=tension[i];
    pp.queryAdd("prefreeze_tension",prefreeze_tension,num_interfaces);

    for (int i=0;i<num_interfaces;i++) {
     tension_slope[i]=0.0;
     tension_T0[i]=293.0;
     tension_min[i]=0.0;
     cap_wave_speed[i]=0.0;
    }

    for (int i=0;i<2*AMREX_SPACEDIM;i++) {
     outflow_velocity_buffer_size[i]=0.0;
    }
    pp.queryAdd("outflow_velocity_buffer_size",
      outflow_velocity_buffer_size,2*AMREX_SPACEDIM);

    pp.queryAdd("tension_slope",tension_slope,num_interfaces);
    pp.queryAdd("tension_T0",tension_T0,num_interfaces);
    pp.queryAdd("tension_min",tension_min,num_interfaces);

    pp.queryAdd("grid_stretching_parameter",grid_stretching_parameter,
	    AMREX_SPACEDIM);
    if (ParallelDescriptor::IOProcessor()) {
     for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
      std::cout << "dir,grid_stretching_parameter " <<
        dir << ' ' << grid_stretching_parameter[dir] << '\n';
     }
    }

    pp.queryAdd("hardwire_Y_gamma",hardwire_Y_gamma,2*num_interfaces);
    pp.queryAdd("hardwire_T_gamma",hardwire_T_gamma,2*num_interfaces);
    pp.queryAdd("accommodation_coefficient",
		 accommodation_coefficient,2*num_interfaces);
    pp.queryAdd("reference_pressure",reference_pressure,2*num_interfaces);
    pp.queryAdd("saturation_temp",saturation_temp,2*num_interfaces);
    pp.queryAdd("saturation_temp_curv",saturation_temp_curv,2*num_interfaces);
    pp.queryAdd("saturation_temp_vel",saturation_temp_vel,2*num_interfaces);
    pp.queryAdd("saturation_temp_min",saturation_temp_min,2*num_interfaces);
    pp.queryAdd("saturation_temp_max",saturation_temp_max,2*num_interfaces);

    pp.queryAdd("nucleation_period",nucleation_period);
    pp.queryAdd("nucleation_init_time",nucleation_init_time);

    pp.queryAdd("perturbation_on_restart",perturbation_on_restart);
    pp.queryAdd("perturbation_mode",perturbation_mode);
    pp.queryAdd("perturbation_eps_temp",perturbation_eps_temp);
    pp.queryAdd("perturbation_eps_vel",perturbation_eps_vel);
   
    pp.queryAdd("solidheat_flag",solidheat_flag);
    if ((solidheat_flag<0)||(solidheat_flag>2))
     amrex::Error("solidheat_flag invalid"); 
 
    pp.queryAdd("microlayer_substrate",microlayer_substrate,num_materials);
    pp.queryAdd("microlayer_angle",microlayer_angle,num_materials);
    pp.queryAdd("microlayer_size",microlayer_size,num_materials);
    pp.queryAdd("macrolayer_size",macrolayer_size,num_materials);
    pp.queryAdd("max_contact_line_size",
                max_contact_line_size,num_materials);
    pp.queryAdd("microlayer_temperature_substrate",
     microlayer_temperature_substrate,num_materials);

    pp.queryAdd("thermal_microlayer_size",thermal_microlayer_size,
		num_materials);
    pp.queryAdd("shear_microlayer_size",shear_microlayer_size,num_materials);
    pp.queryAdd("buoyancy_microlayer_size",buoyancy_microlayer_size,
		num_materials);
    pp.queryAdd("phasechange_microlayer_size",
	phasechange_microlayer_size,num_materials);

    for (int i=0;i<num_materials;i++) {
     if (microlayer_temperature_substrate[i]<0.0)
      amrex::Error("microlayer_temperature_substrate[i]<0.0");
     if ((microlayer_substrate[i]<0)||
         (microlayer_substrate[i]>num_materials))
      amrex::Error("microlayer_substrate invalid");
     if ((microlayer_angle[i]<0.0)|| 
         (microlayer_angle[i]>=NS_PI))
      amrex::Error("microlayer_angle invalid");
     if ((microlayer_size[i]<0.0)||
         (macrolayer_size[i]<microlayer_size[i]))
      amrex::Error("microlayer_size invalid");
     if (max_contact_line_size[i]<0.0)
      amrex::Error("max_contact_line_size invalid");
     if ((microlayer_size[i]>0.0)&&(solidheat_flag==0))
      amrex::Error("cannot have microlayer_size>0.0&&solidheat_flag==0");
     if (thermal_microlayer_size[i]<microlayer_size_default)
      amrex::Error("thermal_microlayer_size too small");
     if (shear_microlayer_size[i]<microlayer_size_default)
      amrex::Error("shear_microlayer_size too small");
     if (buoyancy_microlayer_size[i]<microlayer_size_default)
      amrex::Error("buoyancy_microlayer_size too small");
     if (phasechange_microlayer_size[i]<microlayer_size_default)
      amrex::Error("phasechange_microlayer_size too small");
    }  // i=0..num_materials-1

    pp.queryAdd("nucleation_temp",nucleation_temp,2*num_interfaces);
    pp.queryAdd("nucleation_pressure",nucleation_pressure,2*num_interfaces);
    pp.queryAdd("nucleation_pmg",nucleation_pmg,2*num_interfaces);
    pp.queryAdd("nucleation_mach",nucleation_mach,2*num_interfaces);

    pp.queryAdd("latent_heat",latent_heat,2*num_interfaces);
    pp.queryAdd("latent_heat_slope",latent_heat_slope,2*num_interfaces);
    pp.queryAdd("latent_heat_T0",latent_heat_T0,2*num_interfaces);
    pp.queryAdd("latent_heat_min",latent_heat_min,2*num_interfaces);


    pp.queryAdd("reaction_rate",reaction_rate,2*num_interfaces);
    pp.queryAdd("freezing_model",freezing_model,2*num_interfaces);
    pp.queryAdd("Tanasawa_or_Schrage_or_Kassemi",
      Tanasawa_or_Schrage_or_Kassemi,2*num_interfaces);

    pp.queryAdd("rigid_fraction_id",rigid_fraction_id,num_materials);

    pp.queryAdd("mass_fraction_id",mass_fraction_id,2*num_interfaces);

     // set defaults for "distribute_from_target"
    for (int iten=0;iten<num_interfaces;iten++) {
     for (int ireverse=0;ireverse<2;ireverse++) {
      int iten_local=ireverse*num_interfaces+iten;

      if (freezing_model[iten_local]>=0) {

       Real LL=get_user_latent_heat(iten_local+1,293.0,1);
       if (LL!=0.0) {
        int im1=0;
        int im2=0;
        int im_source=0;
        int im_dest=0;
	 // get_inverse_iten_cpp declared in NavierStokes2.cpp
	 // 1<=im1<im2<=num_materials
        get_inverse_iten_cpp(im1,im2,iten+1);
        if (ireverse==0) {
         im_source=im1;  
         im_dest=im2;  
        } else if (ireverse==1) {
         im_source=im2;  
         im_dest=im1;  
        } else
         amrex::Error("ireverse invalid");

        if ((material_type[im_source-1]==0)&&
            (material_type[im_dest-1]==0)) {
         // do nothing
        } else
         amrex::Error("comp. divu source term model for phase change invalid");

	 // if freezing, we want distribute_from_target (from ice) to
	 // be 1 since the ice is modeled as a rigid material; need
	 // to account for expansion in the liquid.
        Real den_source=denconst[im_source-1];
        Real den_dest=denconst[im_dest-1];
        if ((den_source>0.0)&&(den_dest>0.0)) {

         Real min_den=den_dest;
         Real max_den=den_dest;
	 if (den_source<min_den)
	  min_den=den_source;
	 if (den_source>max_den)
	  max_den=den_source;

	 if (max_den/min_den>=1.0) {
 	  // do nothing
	 } else
	  amrex::Error("max_den or min_den bust");

	 if (is_ice_matC(im_dest-1)==1) { // freezing
          distribute_from_target[iten_local]=1;
	 } else if (is_ice_matC(im_source-1)==1) { // melting
          distribute_from_target[iten_local]=0;
         } else if ((is_ice_matC(im_dest-1)==0)&&
   	            (is_ice_matC(im_source-1)==0)) {

          // fixed_parm=-1,0, or 1.
 	  if (max_den/min_den<1.0) {
 	   amrex::Error("max_den or min_den bust");
	  } else if (max_den/min_den<2.0) {
           distribute_from_target[iten_local]=0;
          } else if (den_dest<den_source) {
            // s= n dot u_dest + mdot/rho_dest
           distribute_from_target[iten_local]=1;
          } else if (den_source<den_dest) {
            // s= n dot u_source + mdot/rho_source
           distribute_from_target[iten_local]=0;
          } else
           amrex::Error("den_source or den_dest invalid");
	 } else {
	  amrex::Error("is_ice_matC bust");
	 }

        } else
         amrex::Error("den_source or den_dest invalid");

	if (Tanasawa_or_Schrage_or_Kassemi[iten_local]==0) { //PD 
 	 //do nothing
	} else if (Tanasawa_or_Schrage_or_Kassemi[iten_local]==1) { //Tanasawa
	 //do nothing
	} else if (Tanasawa_or_Schrage_or_Kassemi[iten_local]==2) { //Schrage
	 //do nothing
	} else if (Tanasawa_or_Schrage_or_Kassemi[iten_local]==3) { //Kassemi
	 //do nothing
	} else
	 amrex::Error("Tanasawa_or_Schrage_or_Kassemi[iten_local] invalid");

       } else if (LL==0.0) {
        distribute_from_target[iten_local]=0;
       } else
        amrex::Error("latent_heat (LL) invalid");
      } else
       amrex::Error("freezing_model invalid");
     } // ireverse
    } // iten

    pp.queryAdd("distribute_from_target",
      distribute_from_target,2*num_interfaces);
    pp.queryAdd("distribute_mdot_evenly",
      distribute_mdot_evenly,2*num_interfaces);
    pp.queryAdd("constant_volume_mdot",
      constant_volume_mdot,2*num_interfaces);

    pp.queryAdd("constant_density_all_time",
      constant_density_all_time,num_materials);
    for (int i=0;i<num_materials;i++) {
     if (material_type[i]==0) {
      // do nothing
     } else if (material_type[i]==999) {
      // do nothing
     } else if ((material_type[i]>0)&&
                (material_type[i]<999)) {
      constant_density_all_time[i]=0;
     } else
      amrex::Error("material_type invalid");
    }

    for (int iten=0;iten<num_interfaces;iten++) {
     for (int ireverse=0;ireverse<2;ireverse++) {
      int iten_local=ireverse*num_interfaces+iten;

      int im1=0;
      int im2=0;
      int im_source=0;
      int im_dest=0;

      if (freezing_model[iten_local]>=0) {

       Real LL=get_user_latent_heat(iten_local+1,293.0,1);
       if (LL!=0.0) {
        get_inverse_iten_cpp(im1,im2,iten+1);
        if ((im1>=1)&&(im1<=num_materials)&&(im2>=1)&&(im2<=num_materials)&&
            (im1!=im2)) { 
         // do nothing
        } else
         amrex::Error("im1 or im2 invalid");
        
        if (ireverse==0) {
         im_source=im1;  
         im_dest=im2;  
        } else if (ireverse==1) {
         im_source=im2;  
         im_dest=im1;  
        } else
         amrex::Error("ireverse invalid");

        int im_constant=0;
        if (constant_volume_mdot[iten_local]==0) {
         // do nothing
        } else if (constant_volume_mdot[iten_local]==1) {
         im_constant=im_source;
        } else if (constant_volume_mdot[iten_local]==-1) {
         im_constant=im_dest;
        } else
         amrex::Error("constant_volume_mdot[iten_local] invalid");

        if (im_constant==0) {
         // do nothing
        } else if ((im_constant>=1)&&(im_constant<=num_materials)) {
         constant_density_all_time[im_constant-1]=0;
        } else
         amrex::Error("im_constant invalid");
 
       } else if (LL==0.0) {
        // do nothing
       } else
        amrex::Error("latent_heat (LL) invalid");
      } else
       amrex::Error("freezing_model invalid");

      Real LL=get_user_latent_heat(iten_local+1,293.0,1);
      if (LL==0.0) {

       if (distribute_from_target[iten_local]==0) {
        // do nothing
       } else
        amrex::Error("distribute_from_target[iten_local] invalid");

       if (distribute_mdot_evenly[iten_local]==0) {
        // do nothing
       } else
        amrex::Error("distribute_mdot_evenly[iten_local] invalid");

       if (constant_volume_mdot[iten_local]==0) {
        // do nothing
       } else
        amrex::Error("constant_volume_mdot[iten_local] invalid");

      } else if (LL!=0.0) {

       if ((distribute_from_target[iten_local]==0)||
           (distribute_from_target[iten_local]==1)) {

        if ((distribute_mdot_evenly[iten_local]==0)||
            (distribute_mdot_evenly[iten_local]==1)||
            (distribute_mdot_evenly[iten_local]==2)) {
         // do nothing
        } else
         amrex::Error("distribute_mdot_evenly invalid");

       } else
        amrex::Error("distribute_from_target invalid");

      } else
       amrex::Error("latent_heat (LL) invalid");

     } // ireverse
    } // iten

    Vector<int> preset_flag;
    preset_flag.resize(num_interfaces);

    for (int im=0;im<num_materials;im++) {
     for (int im_opp=im+1;im_opp<num_materials;im_opp++) {
      int iten=0;
      get_iten_cpp(im+1,im_opp+1,iten);
      if ((iten<1)||(iten>num_interfaces))
       amrex::Error("iten invalid");
      preset_flag[iten-1]=-1;
      if ((material_type[im]==999)||
          (material_type[im_opp]==999)) {
       material_type_interface[iten-1]=999;
       preset_flag[iten-1]=999;
      } else if ((material_type[im]==0)||
                 (material_type[im_opp]==0)) {
       material_type_interface[iten-1]=0;
       preset_flag[iten-1]=0;
      } else if ((latent_heat[iten-1]!=0.0)||
    	         (latent_heat[iten-1+num_interfaces]!=0.0)) {
       material_type_interface[iten-1]=0;
       preset_flag[iten-1]=0;
      } else {
       material_type_interface[iten-1]=material_type[im];
      }
     } // im_opp=im+1;im_opp<navierStokes::num_materials;im_opp++
    } // im=0;im<navierstokes::num_materials;im++

    pp.queryAdd("material_type_interface",
      material_type_interface,num_interfaces);

    for (int im=0;im<num_materials;im++) {
     for (int im_opp=im+1;im_opp<num_materials;im_opp++) {
      int iten=0;
      get_iten_cpp(im+1,im_opp+1,iten);
      if ((iten<1)||(iten>NavierStokes::num_interfaces))
       amrex::Error("iten invalid");
      if (preset_flag[iten-1]==-1) {
       if ((material_type_interface[iten-1]>=0)&&
           (material_type_interface[iten-1]<=999)) {
        //do nothing
       } else
        amrex::Error("material_type_interface invalid");
      } else if (preset_flag[iten-1]>=0) {
       if (material_type_interface[iten-1]==preset_flag[iten-1]) {
        //do nothing
       } else
        amrex::Error("material_type_interface invalid");
      } else
       amrex::Error("preset flag invalid");
     } // im_opp=im+1;im_opp<num_materials;im_opp++
    } // im=0;im<num_materials;im++


    for (int im=0;im<num_materials;im++) {

     if (is_ice_matC(im)==1) {

      int ispec=rigid_fraction_id[im];

      if ((ispec>=1)&&(ispec<=num_species_var)) {
       for (int im_opp=0;im_opp<num_materials;im_opp++) {
        if ((speciesconst[(ispec-1)*num_materials+im_opp]>0.0)&&
            (speciesconst[(ispec-1)*num_materials+im_opp]<=1.0)) {
         //do nothing
	} else
	 amrex::Error("speciesconst invalid");
        if (speciesreactionrate[(ispec-1)*num_materials+im_opp]>=0.0) {
  	 //do nothing
	} else
	 amrex::Error("speciesreactionrate invalid");
       } //im_opp=0..num_materials-1
      } else
       amrex::Error("rigid_fraction_id (ispec) invalid");

      for (int im_opp=0;im_opp<num_materials;im_opp++) {
       if (im!=im_opp) {
        if (ns_is_rigid(im_opp)==0) {
         int iten;
         get_iten_cpp(im+1,im_opp+1,iten); //declared in NavierStokes2.cpp
         if ((iten<1)||(iten>num_interfaces))
          amrex::Error("iten invalid");
         Real LL1=get_user_latent_heat(iten,293.0,1);
         Real LL2=get_user_latent_heat(iten+num_interfaces,293.0,1);

         if ((LL1!=0.0)||(LL2!=0.0)) {

  	  int local_distribute=-1;

 	  int im_source=-1;
	  int im_dest=-1;
	  int im_ice=im;
	  if ((LL1!=0.0)&&(LL2==0.0)) {
           local_distribute=distribute_from_target[iten-1];
           if (im<im_opp) {
  	    im_source=im;
	    im_dest=im_opp;
	   } else if (im>im_opp) {
 	    im_source=im_opp;
	    im_dest=im;
           } else
  	    amrex::Error("im or im_opp bust");
	  } else if ((LL1==0.0)&&(LL2!=0.0)) {
           local_distribute=distribute_from_target[iten+num_interfaces-1];
           if (im>im_opp) {
  	    im_source=im;
	    im_dest=im_opp;
	   } else if (im<im_opp) {
 	    im_source=im_opp;
	    im_dest=im;
           } else
  	    amrex::Error("im or im_opp bust");
	  } else if ((LL1!=0.0)&&(LL2!=0.0)) {
           amrex::Error("cannot do both melting and freezing yet...");
	  } else
  	   amrex::Error("LL1 or LL2 bust");

	  if (im_source==im_ice) { //melting
           if (local_distribute==0) {
	    // distribute_from_target=false
	    // distribute_to_target=true 
	    // source=ice  dest (target)=melt
	   } else 
            amrex::Error("distribute_from_target should be 0(melting)");
	  } else if (im_dest==im_ice) {
           if (local_distribute==1) {
	    // distribute_from_target=true
	    // distribute_to_target=false
	    // source=melt  dest (target)=ice
	   } else 
            amrex::Error("distribute_from_target should be 1(freezing)");
	  } else
           amrex::Error("im_ice invalid");

         } else if ((LL1==0.0)&&(LL2==0.0)) {
          // do nothing
         } else
          amrex::Error("LL1 or LL2 bad");
        } else if (ns_is_rigid(im_opp)==1) {
         //do nothing
        } else
         amrex::Error("ns_is_rigid invalid");
       } else if (im==im_opp) { 
        //do nothing
       } else
        amrex::Error("im or im_opp bust");
      } // im_opp=0;im_opp<num_materials;im_opp++
     } else if (is_ice_matC(im)==0) {
      //do nothing
     } else
      amrex::Error("is_ice_matC invalid");

     if (is_FSI_rigid_matC(im)==1) {

      int ispec=rigid_fraction_id[im];

      if ((ispec>=1)&&(ispec<=num_species_var)) {
       for (int im_opp=0;im_opp<num_materials;im_opp++) {
        if (im_opp==im) {
         if (speciesconst[(ispec-1)*num_materials+im]==1.0) {
  	  //do nothing
	 } else
	  amrex::Error("speciesconst invalid");
         if (speciesreactionrate[(ispec-1)*num_materials+im]>0.0) {
  	  //do nothing
	 } else
	  amrex::Error("speciesreactionrate invalid");
	} else if (im_opp!=im) {
         if (speciesconst[(ispec-1)*num_materials+im_opp]==1.0) {
  	  //do nothing
	 } else
	  amrex::Error("speciesconst invalid");
         if (speciesreactionrate[(ispec-1)*num_materials+im_opp]>=0.0) {
  	  //do nothing
	 } else
	  amrex::Error("speciesreactionrate invalid");
	} else 
	 amrex::Error("im_opp or im invalid");
       } //im_opp=0..num_materials-1
      } else
       amrex::Error("rigid_fraction_id (ispec) invalid");

     } else if (is_FSI_rigid_matC(im)==0) {
      //do nothing
     } else
      amrex::Error("is_FSI_rigid_matC invalid");

    } // im=0;im<num_materials;im++

    pp.queryAdd("R_Palmore_Desjardins",R_Palmore_Desjardins);

    for (int im=0;im<num_materials;im++) {

     if ((override_density[im]!=0)&&
         (override_density[im]!=1)&& //rho=rho(T,Y,z)
         (override_density[im]!=2)) {//Boussinesq approximation
      std::cout << "num_materials= " << num_materials << '\n';
      std::cout << "im=" << im << " override_density[im]= " <<
	     override_density[im] << '\n';
      amrex::Error("override_density invalid (1) ");
     }
     if (DrhoDT[im]>0.0)
      amrex::Error("DrhoDT cannot be positive");

     if (material_type[im]==999) {
      if (viscconst[im]<=0.0)
       amrex::Error("solid cannot have 0 viscosity");
      if (ns_is_rigid(im)!=1)
       amrex::Error("ns_is_rigid invalid");
      if (override_density[im]!=0)
       amrex::Error("override_density invalid");
     } else if (material_type[im]==0) {
      if (ns_is_rigid(im)!=0)
       amrex::Error("ns_is_rigid invalid");
     } else if ((material_type[im]>0)&& 
                (material_type[im]<999)) {
      if (override_density[im]!=0)
       amrex::Error("override_density invalid");
      if (ns_is_rigid(im)!=0)
       amrex::Error("ns_is_rigid invalid");
     } else {
      amrex::Error("material type invalid");
     }

    }  // im=0, im<num_materials

    if (num_state_base!=2)
     amrex::Error("num_state_base invalid 10");

    for (int i=0;i<num_interfaces;i++) {

     if (is_valid_freezing_model(freezing_model[i])==1) {
      // do nothing
     } else
      amrex::Error("freezing_model invalid in read_params (i)");
     if (is_valid_freezing_model(freezing_model[i+num_interfaces])==1) {
      // do nothing
     } else
      amrex::Error("freezing_model invalid in read_params (i+num_interfaces)");

     if ((distribute_from_target[i]<0)||(distribute_from_target[i]>1))
      amrex::Error("distribute_from_target invalid in read_params (i)");
     if ((distribute_from_target[i+num_interfaces]<0)||
	 (distribute_from_target[i+num_interfaces]>1))
      amrex::Error("distribute_from_target invalid in read_params (i+num_interfaces)");
     if (mass_fraction_id[i]<0)
      amrex::Error("mass_fraction_id invalid in read_params (i)");
     if (mass_fraction_id[i+num_interfaces]<0)
      amrex::Error("mass_fraction_id invalid in read_params (i+num_interfaces)");
     if ((distribute_mdot_evenly[i]<0)||
         (distribute_mdot_evenly[i]>2))
      amrex::Error("distribute_mdot_evenly invalid in read_params (i)");
     if ((distribute_mdot_evenly[i+num_interfaces]<0)||
         (distribute_mdot_evenly[i+num_interfaces]>2))
      amrex::Error("distribute_mdot_evenly invalid in read_params (i+num_interfaces)");
     if ((constant_volume_mdot[i]<-1)||
         (constant_volume_mdot[i]>1))
      amrex::Error("constant_volume_mdot invalid in read_params (i)");
     if ((constant_volume_mdot[i+num_interfaces]<-1)||
         (constant_volume_mdot[i+num_interfaces]>1))
      amrex::Error("constant_volume_mdot invalid in read_params (i+num_interfaces)");
    }  // i=0..num_interfaces-1


    shock_timestep.resize(num_materials);
    for (int i=0;i<num_materials;i++) {
     if (rigid_fraction_id[i]<0)
      amrex::Error("rigid_fraction_id invalid in read_params");
     shock_timestep[i]=0;
    }
    pp.queryAdd("shock_timestep",shock_timestep,num_materials);

    for (int i=0;i<num_materials;i++) {
     if (!((shock_timestep[i]==1)|| //always take into account acoustic waves
	   (shock_timestep[i]==0)|| //taken into account acoustic waves t=0
           (shock_timestep[i]==2))) //never consider acoustic waves.
      amrex::Error("shock_timestep invalid");
    }

    projection_pressure_scale=1.0;

    if (some_materials_compressible()==1) {
     projection_pressure_scale=1.0e+6;
    }

    pp.queryAdd("projection_pressure_scale",projection_pressure_scale);
    if (projection_pressure_scale<=0.0)
     amrex::Error("projection pressure scale invalid");

    projection_velocity_scale=std::sqrt(projection_pressure_scale);

    num_divu_outer_sweeps=1;

    if (some_materials_compressible()==1) {
     num_divu_outer_sweeps=2;
    } else if (some_materials_compressible()==0) {
     // do nothing
    } else {
     amrex::Error("some_materials_compressible invalid");
    }

    if (FSI_material_exists_CTML()==1) {
     num_divu_outer_sweeps=2;
    } else if (FSI_material_exists_CTML()==0) {
     // do nothing
    } else
     amrex::Error("FSI_material_exists_CTML() invalid");

    pp.queryAdd("num_divu_outer_sweeps",num_divu_outer_sweeps);

    if (some_materials_compressible()==1) {
     if (num_divu_outer_sweeps<2) 
      amrex::Error("num_divu_outer_sweeps must be >1 (compres)");
    } else if (some_materials_compressible()==0) {
     // do nothing
    } else {
     amrex::Error("some_materials_compressible invalid");
    }

    if (FSI_material_exists_CTML()==1) {
     if (num_divu_outer_sweeps<2) 
      amrex::Error("num_divu_outer_sweeps must be >1 (CTML coupling)");
    } else if (FSI_material_exists_CTML()==0) {
     // do nothing
    } else
     amrex::Error("FSI_material_exists_CTML() invalid");
     
    pp.queryAdd("post_init_pressure_solve",post_init_pressure_solve);
    if ((post_init_pressure_solve<0)||(post_init_pressure_solve>1))
     amrex::Error("post_init_pressure_solve out of range");

    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {

     if (geometry_is_periodic[dir]==1) {
      // do nothing
     } else if (geometry_is_periodic[dir]==0) {

      if (lo_bc[dir]==Interior) {
       amrex::Error("cannot have Interior lo_bc if not periodic");
      } else if (lo_bc[dir]==Symmetry) {
       // do nothing
      } else if (lo_bc[dir]==Inflow) {
       // do nothing
      } else if (lo_bc[dir]==SlipWall) {
       // do nothing
      } else if (lo_bc[dir]==NoSlipWall) {
       // do nothing
      } else if (lo_bc[dir]==Outflow) {
       // do nothing
      } else
       amrex::Error("lo_bc[dir] not recognized");

      if (hi_bc[dir]==Interior) {
       amrex::Error("cannot have Interior hi_bc if not periodic");
      } else if (hi_bc[dir]==Symmetry) {
       // do nothing
      } else if (hi_bc[dir]==Inflow) {
       // do nothing
      } else if (hi_bc[dir]==SlipWall) {
       // do nothing
      } else if (hi_bc[dir]==NoSlipWall) {
       // do nothing
      } else if (hi_bc[dir]==Outflow) {
       // do nothing
      } else
       amrex::Error("hi_bc[dir] not recognized");

     } else
      amrex::Error("geometry_is_periodic[dir] invalid"); 
    } // dir=0..sdim-1

    pp.queryAdd("pressure_select_criterion",pressure_select_criterion);
    if ((pressure_select_criterion<0)||
        (pressure_select_criterion>2))
     amrex::Error("pressure_select_criterion invalid");

    elastic_time.resize(num_materials);

    Carreau_alpha.resize(num_materials);
    Carreau_beta.resize(num_materials);
    Carreau_n.resize(num_materials);
    Carreau_mu_inf.resize(num_materials);
    shear_thinning_fluid.resize(num_materials);

    polymer_factor.resize(num_materials);

    compressible_dt_factor.resize(num_materials);
    for (int i=0;i<num_materials;i++) {
     compressible_dt_factor[i]=1.0;
    }

    pp.queryAdd("compressible_dt_factor",
       compressible_dt_factor);

    for (int i=0;i<num_materials;i++) {
     if ((compressible_dt_factor[i]>=1.0)&&
         (compressible_dt_factor[i]<=1.0e+12)) {
      // do nothing
     } else
      amrex::Error("compressible_dt_factor invalid");
    }

    pp.queryAdd("disable_advection",disable_advection);

    if ((disable_advection==0)||(disable_advection==1)) {
     // do nothing
    } else
     amrex::Error("disable_advection invalid 1");

    pp.queryAdd("disable_pressure_solve",disable_pressure_solve);

    if ((disable_pressure_solve==0)||(disable_pressure_solve==1)) {
     // do nothing
    } else
     amrex::Error("disable_pressure_solve invalid 1");

    for (int i=0;i<num_materials;i++) {
     elastic_time[i]=0.0;

     Carreau_alpha[i]=1.0;
     Carreau_beta[i]=0.0;
     Carreau_n[i]=1.0;
     Carreau_mu_inf[i]=0.0;
     shear_thinning_fluid[i]=0;

     polymer_factor[i]=0.0;
    }  // i=0..num_materials-1

    pp.queryAdd("elastic_time",elastic_time,num_materials);

    for (int i=0;i<num_materials;i++) {
     if (elastic_viscosity[i]>=0.0) {
      if (fort_built_in_elastic_model(&elastic_viscosity[i],
                                      &viscoelastic_model[i])==1) {
       if (viscoelastic_model[i]==3) { // incremental elastic model
        if (elastic_time[i]>=1.0e+8) {
         // do nothing
        } else
         amrex::Error("elastic time inconsistent with model");
       } else if (viscoelastic_model[i]==7) { // incremental Neo-Hookean model
        if (elastic_time[i]>=1.0e+8) {
         // do nothing
        } else
         amrex::Error("elastic time inconsistent with model");
       } else if (fort_built_in_elastic_model(&elastic_viscosity[i],
	     		                    &viscoelastic_model[i])==1) {
        // do nothing
       } else
        amrex::Error("fort_built_in_elastic_model invalid");
      } else if (fort_built_in_elastic_model(&elastic_viscosity[i],
		 	                   &viscoelastic_model[i])==0) {
       //do nothing
      } else
       amrex::Error("fort_built_in_elastic_model invalid");
     } else
      amrex::Error("elastic_viscosity invalid");
    } // i=0..num_materials-1

    pp.queryAdd("polymer_factor",polymer_factor,num_materials);

    pp.queryAdd("Carreau_alpha",Carreau_alpha,num_materials);
    pp.queryAdd("Carreau_beta",Carreau_beta,num_materials);
    pp.queryAdd("Carreau_n",Carreau_n,num_materials);
    pp.queryAdd("Carreau_mu_inf",Carreau_mu_inf,num_materials);

    etaL.resize(num_materials);
    etaP.resize(num_materials);
    etaS.resize(num_materials);
    concentration.resize(num_materials);

    for (int i=0;i<num_materials;i++) {

     if (Carreau_n[i]>1.0)
      amrex::Error("Carreau_n[i] invalid");
     if (Carreau_mu_inf[i]<0.0)
      amrex::Error("Carreau_mu_inf[i] invalid");

     if (viscosity_state_model[i]>=0) {
      // do nothing
     } else
      amrex::Error("viscosity state model invalid");

     if (fort_built_in_elastic_model(&elastic_viscosity[i],
			           &viscoelastic_model[i])==1) {
      // do nothing
     } else if (fort_built_in_elastic_model(&elastic_viscosity[i],
	 	                          &viscoelastic_model[i])==0) {
      // do nothing
     } else
      amrex::Error("fort_built_in_elastic_model invalid");

     if (les_model[i]<0)
      amrex::Error("les model invalid");

     if ((elastic_time[i]<0.0)||(elastic_viscosity[i]<0.0))
      amrex::Error("elastic_time/elastic_viscosity invalid read_params");

     // (1/L) eps=0.01
     // viscoelastic_model=0  FENE CR   trac(A)<L^2 lambda(A)>eps/L^2
     // viscoelastic_model=1  Oldroyd B trac(A)<inf lambda(A)>eps/L^2
     // viscoelastic_model=5  FENE P    trac(A)<L^2 lambda(A)>eps/L^2
     // viscoelastic_model=6  linearPTT trac(A)<inf lambda(A)>eps/L^2
     if (polymer_factor[i]>=0.0) {
      if ((viscoelastic_model[i]==0)|| //FENE CR
          (viscoelastic_model[i]==1)|| //Oldroyd B
          (viscoelastic_model[i]==5)|| //FENE-P
          (viscoelastic_model[i]==6)) {//linearPTT

       if (elastic_viscosity[i]>0.0) {
        if (polymer_factor[i]>0.0) {
         // do nothing
        } else {
         std::cout << "----------------------------------\n";
         std::cout << "need polymer_factor > 0 viscoelastic_model=0,1,5,6\n";
         std::cout << "lambda(A)>0.01/L^2 \n";
         amrex::Error("polymer_factor invalid");
        }
       } else if (elastic_viscosity[i]==0.0) {
        // do nothing
       } else
        amrex::Error("elastic_viscosity invalid");

      } else if ((viscoelastic_model[i]==7)|| //Neo-Hookean
   	         (viscoelastic_model[i]==3)) { //incremental model
       // do nothing
      } else
       amrex::Error("viscoelastic_model invalid");
	       
     } else
      amrex::Error("polymer_factor invalid");

     if ((Carreau_beta[i]!=0.0)&&(visc_coef==0.0))
      amrex::Error("Cannot have Carreau_beta!=0 and visc_coef==0 ");

     int ip1=i+1;
     if (fort_is_rigid_base(&FSI_flag[i],&ip1)==1) {
      shear_thinning_fluid[i]=0;
     } else if (fort_is_rigid_base(&FSI_flag[i],&ip1)==0) {
      shear_thinning_fluid[i]=0;
      if ((probtype==2)&&(axis_dir>0)&&(i==0))
       shear_thinning_fluid[i]=1;
      if (Carreau_beta[i]!=0.0)
       shear_thinning_fluid[i]=1;

      if (shear_thinning_fluid[i]==1) {
       // do nothing
      } else if (shear_thinning_fluid[i]==0) {
       // do nothing
      } else
       amrex::Error("shear_thinning_fluid invalid");
     } else 
      amrex::Error("fort_is_rigid_base invalid");

      // if first material and Carreau_beta==0 for first material,
      //  probtype==2, axis_dir>0, then 
      //  "call viscosity(axis_dir,visc(D_DECL(i,j,k)),shear)"
      //  VISC_RAW=visc(D_DECL(i,j,k))
      // otherwise:
      //  if Carreau_beta==0, then
      //   VISC_RAW=viscconst
      //  else if (Carreau_beta>0) then
      //   if visco-elastic then
      //    VISC_RAW=( etaL0-etaP0+ 
      //     etaP0*(1+(beta*gamma_dot)**alpha)**( (n-1)/alpha )  )
      //    =
      //    VISC_RAW=( etaS+ 
      //     etaP0*(1+(beta*gamma_dot)**alpha)**( (n-1)/alpha )  )
      //
      //   else
      //    VISC_RAW=(
      //      mu_inf + (etaL0-mu_inf)* 
      //      (1+(beta*gamma_dot)**alpha)**((n-1)/alpha) )
      //
      //   VISCOSITY=visc_coef * VISC_RAW
      //
      // my "etaL0" is Mitsuhiro's eta_S(1+c_0)= viscconst
      // my "etaP0" is Mitsuhiro's c_0 eta_S   = elastic_viscosity
      // my "etaS" is Mitsuhiro's eta_S too.  = viscconst-elastic_viscosity
      // The coefficient for the viscoelastic force term is:
      //  visc_coef * 
      //  (VISC_RAW-viscconst+elastic_viscosity)/(elastic_time*(1-tr(A)/L^2))
      // =
      //  visc_coef * 
      //  (VISC_RAW-etaL0+etaP0)/(elastic_time*(1-tr(A)/L^2))
      // =
      //  visc_coef * 
      //  (VISC_RAW-etaS)/(elastic_time*(1-tr(A)/L^2))
      // =  (assume visc_coef==1)
      //  etaP0*
      //  (1+(beta*gamma_dot)**alpha)**((n-1)/alpha) 
      //
     etaL[i]=viscconst[i];  
     etaP[i]=elastic_viscosity[i]; // eta_{P0} in the JNNFM paper and above.
     etaS[i]=etaL[i]-etaP[i];  
     if (ParallelDescriptor::IOProcessor()) {
      std::cout << "for material " << i << '\n';

      std::cout << "static_damping_coefficient[i]=" << 
	      static_damping_coefficient[i] << '\n';

      std::cout << "etaL0=viscconst[i]=" << etaL[i] << '\n';

      std::cout << "Carreau_alpha=" << Carreau_alpha[i] << '\n';
      std::cout << "Carreau_beta=" << Carreau_beta[i] << '\n';
      std::cout << "Carreau_n=" << Carreau_n[i] << '\n';
      std::cout << "Carreau_mu_inf=" << Carreau_mu_inf[i] << '\n';
      std::cout << "shear_thinning_fluid=" << shear_thinning_fluid[i] << '\n';
     } // io processor

     if ((num_materials_viscoelastic>=1)&&
         (num_materials_viscoelastic<=num_materials)) {

      if (ParallelDescriptor::IOProcessor()) {
       std::cout << "for material " << i << '\n';
       std::cout << "etaP0=elastic_viscosity=" << etaP[i] << '\n';
       std::cout << "etaS=etaL0-etaP0= " << etaS[i] << '\n';
       std::cout << "elastic_viscosity= " << elastic_viscosity[i] << '\n';
       std::cout << "store_elastic_data= " << store_elastic_data[i] << '\n';
       std::cout << "elastic_time= " << elastic_time[i] << '\n';
      }

       // c0=etaP0/etaS=etaP0/(etaL0-etaP0)
      concentration[i]=0.0;
      if (etaL[i]-etaP[i]>0.0)
       concentration[i]=etaP[i]/(etaL[i]-etaP[i]);  

      if (ParallelDescriptor::IOProcessor()) {
       std::cout << "for material " << i << " c0=" << 
        concentration[i] << '\n';
      }

     } else if (num_materials_viscoelastic==0) {
      // do nothing
     } else
      amrex::Error("num_materials_viscoelastic invalid");

    } // i=0..num_materials-1

    pp.queryAdd("wait_time",wait_time);

    pp.get("advbot",advbot);
    pp.queryAdd("inflow_pressure",inflow_pressure);
    pp.queryAdd("outflow_pressure",outflow_pressure);

    pp.queryAdd("multilevel_maxcycle",multilevel_maxcycle);
    pp.queryAdd("multilevel_restart_period",multilevel_restart_period);
    if ((multilevel_maxcycle>1)&&(multilevel_restart_period>1)) {
     // do nothing
    } else
     amrex::Error("multilevel_maxcycle or multilevel_restart_period bad");


    ParmParse ppmac("mac");
    ParmParse ppcg("cg");

    int cg_restart_period=ABecLaplacian::CG_def_restart_period;
    ppcg.queryAdd("restart_period",cg_restart_period);
    int cg_maxiter=ABecLaplacian::CG_def_maxiter;
    ppcg.queryAdd("maxiter",cg_maxiter);

    ppmac.queryAdd( "mac_abs_tol",mac_abs_tol);
      // mac.visc_abs_tol (not ns.visc_abs_tol)
    ppmac.queryAdd( "visc_abs_tol",visc_abs_tol);
    thermal_abs_tol=visc_abs_tol;
    ppmac.queryAdd( "thermal_abs_tol",thermal_abs_tol);
    pp.queryAdd( "minimum_relative_error",minimum_relative_error);
    pp.queryAdd( "diffusion_minimum_relative_error",
      diffusion_minimum_relative_error);
    if (mac_abs_tol<=0.0)
     amrex::Error("mac_abs_tol must be positive");
    if (visc_abs_tol<=0.0)
     amrex::Error("visc_abs_tol must be positive");
    if (thermal_abs_tol<=0.0)
     amrex::Error("thermal_abs_tol must be positive");

    extra_circle_parameters(
       xblob2,yblob2,zblob2,radblob2,
       xblob3,yblob3,zblob3,radblob3,
       xblob4,yblob4,zblob4,radblob4,
       xblob5,yblob5,zblob5,radblob5,
       xblob6,yblob6,zblob6,radblob6,
       xblob7,yblob7,zblob7,radblob7,
       xblob8,yblob8,zblob8,radblob8,
       xblob9,yblob9,zblob9,radblob9,
       xblob10,yblob10,zblob10,radblob10 );

    pp.queryAdd("xblob2",xblob2);
    pp.queryAdd("yblob2",yblob2);
    pp.queryAdd("zblob2",zblob2);
    pp.queryAdd("radblob2",radblob2);

    pp.queryAdd("xblob3",xblob3);
    pp.queryAdd("yblob3",yblob3);
    pp.queryAdd("zblob3",zblob3);
    pp.queryAdd("radblob3",radblob3);

    pp.queryAdd("xblob4",xblob4);
    pp.queryAdd("yblob4",yblob4);
    pp.queryAdd("zblob4",zblob4);
    pp.queryAdd("radblob4",radblob4);

    pp.queryAdd("xblob5",xblob5);
    pp.queryAdd("yblob5",yblob5);
    pp.queryAdd("zblob5",zblob5);
    pp.queryAdd("radblob5",radblob5);

    pp.queryAdd("xblob6",xblob6);
    pp.queryAdd("yblob6",yblob6);
    pp.queryAdd("zblob6",zblob6);
    pp.queryAdd("radblob6",radblob6);

    pp.queryAdd("xblob7",xblob7);
    pp.queryAdd("yblob7",yblob7);
    pp.queryAdd("zblob7",zblob7);
    pp.queryAdd("radblob7",radblob7);

    pp.queryAdd("xblob8",xblob8);
    pp.queryAdd("yblob8",yblob8);
    pp.queryAdd("zblob8",zblob8);
    pp.queryAdd("radblob8",radblob8);

    pp.queryAdd("xblob9",xblob9);
    pp.queryAdd("yblob9",yblob9);
    pp.queryAdd("zblob9",zblob9);
    pp.queryAdd("radblob9",radblob9);

    pp.queryAdd("xblob10",xblob10);
    pp.queryAdd("yblob10",yblob10);
    pp.queryAdd("zblob10",zblob10);
    pp.queryAdd("radblob10",radblob10);

    xactive=0.0;
    yactive=0.0;
    zactive=0.0;
    ractive=0.0;
    ractivex=0.0;
    ractivey=0.0;
    ractivez=0.0;

    pp.queryAdd("xactive",xactive);
    pp.queryAdd("yactive",yactive);
    pp.queryAdd("zactive",zactive);
    pp.queryAdd("ractive",ractive);
    if (ractive>0.0) {
     ractivex=ractive;
     ractivey=ractive;
     ractivez=ractive;
    }
    pp.queryAdd("ractivex",ractivex);
    pp.queryAdd("ractivey",ractivey);
    pp.queryAdd("ractivez",ractivez);

     // 0 - MGPCG  1-PCG 2-MINV=I
    pp.queryAdd("project_solver_type",project_solver_type);
    pp.queryAdd("initial_cg_cycles",initial_cg_cycles);

    pp.queryAdd("initial_project_cycles",initial_project_cycles);
    if (initial_project_cycles<1)
     amrex::Error("must do at least 1 jacobi sweep at the beginning");
    initial_viscosity_cycles=initial_project_cycles;
    pp.queryAdd("initial_viscosity_cycles",initial_viscosity_cycles);
    if (initial_viscosity_cycles<1)
     amrex::Error("must do at least 1 jacobi sweep at the beginning (visc)");
    initial_thermal_cycles=initial_viscosity_cycles;
    pp.queryAdd("initial_thermal_cycles",initial_thermal_cycles);
    if (initial_thermal_cycles<1)
     amrex::Error("must do at least 1 jacobi sweep at the beginning (therm)");

    if ((project_solver_type>=0)&&(project_solver_type<=2)) {
     // do nothing
    } else
     amrex::Error("project_solver_type invalid");

    prescribe_temperature_outflow=0;
    pp.queryAdd("prescribe_temperature_outflow",prescribe_temperature_outflow);
    if ((prescribe_temperature_outflow<0)||
        (prescribe_temperature_outflow>3))
     amrex::Error("prescribe_temperature_outflow invalid");

    is_phasechange=0;
    for (int i=0;i<2*num_interfaces;i++) {
     Real LL=get_user_latent_heat(i+1,293.0,1);
     if (LL!=0.0) {
      is_phasechange=1;
      if (is_valid_freezing_model(freezing_model[i])==1) {
       // do nothing
      } else
       amrex::Error("freezing_model[i] invalid");
     } else if (LL==0.0) {
      //do nothing
     } else {
      amrex::Error("LL is NaN");
     } 
    }  // i=0;i<2*num_interfaces

    hydrate_flag=0;
    for (int i=0;i<2*num_interfaces;i++) {
     Real LL=get_user_latent_heat(i+1,293.0,1);
     if (LL!=0.0) {
      if (is_hydrate_freezing_model(freezing_model[i])==1) {
       hydrate_flag=1;
      } else if (is_hydrate_freezing_model(freezing_model[i])==0) {
       // do nothing
      } else
       amrex::Error("is_hydrate_freezing_model bust");
     } else if (LL==0.0) {
      // do nothing
     } else
      amrex::Error("latent_heat (LL) is NaN");
    } // i

    truncate_volume_fractions.resize(num_materials);

    for (int i=0;i<num_materials;i++) {

     if ((FSI_flag[i]==FSI_FLUID)|| // tessellating
         (FSI_flag[i]==FSI_FLUID_NODES_INIT))  // fluid, tessellating
      truncate_volume_fractions[i]=1;
     else if (is_ice_matC(i)==1) // ice, tessellating
      truncate_volume_fractions[i]=1;
     else if (FSI_flag[i]==FSI_PRESCRIBED_PROBF90) 
      truncate_volume_fractions[i]=0;
     else if (FSI_flag[i]==FSI_PRESCRIBED_NODES) 
      truncate_volume_fractions[i]=0;
     else if (FSI_flag[i]==FSI_SHOELE_CTML) 
      truncate_volume_fractions[i]=0;
     else if (is_FSI_rigid_matC(i)==1) // FSI rigid solid, tessellating
      truncate_volume_fractions[i]=0;
     else
      amrex::Error("FSI_flag invalid");
    }  // i=0..num_materials-1

     //default=4
    pp.queryAdd("particle_nsubdivide",particle_nsubdivide);
     //default=10
    pp.queryAdd("particle_max_per_nsubdivide",particle_max_per_nsubdivide);

    pp.queryAdd("truncate_volume_fractions",truncate_volume_fractions,
		num_materials);

    for (int i=0;i<num_materials;i++) {
     if ((truncate_volume_fractions[i]<0)||
         (truncate_volume_fractions[i]>1))
      amrex::Error("truncate_volume_fractions invalid");
    }

    pp.queryAdd("truncate_thickness",truncate_thickness);
    if (truncate_thickness<1.0)
     amrex::Error("truncate_thickness too small");

    for (int im=1;im<=num_materials;im++) {
     for (int im_opp=im+1;im_opp<=num_materials;im_opp++) {
      for (int ireverse=0;ireverse<=1;ireverse++) {

       if ((im>num_materials)||(im_opp>num_materials))
        amrex::Error("im or im_opp bust 200cpp");
       int iten;
       get_iten_cpp(im,im_opp,iten); //declared in NavierStokes2.cpp
       if ((iten<1)||(iten>num_interfaces))
        amrex::Error("iten invalid");
       int im_source=im;
       int im_dest=im_opp;
       if (ireverse==1) {
        im_source=im_opp;
        im_dest=im;
       }

       int indexEXP=iten+ireverse*num_interfaces-1;

       Real LL=get_user_latent_heat(indexEXP+1,293.0,1);

       if (LL!=0.0) {
        if ((truncate_volume_fractions[im-1]==0)&&
  	    (truncate_volume_fractions[im_opp-1]==0)) {
	 // do nothing
	} else {
	 std::cout << "WARNING: (if microscopic seeds) \n";
	 std::cout << "all materials at mass transfer interface\n";
	 std::cout << "should have truncate_volume_fractions==0\n";
	 std::cout << "im= " << im << '\n';
	 std::cout << "truncate_volume_fractions[im-1]= " << 
  	   truncate_volume_fractions[im-1] << '\n';
	 std::cout << "im_opp= " << im_opp << '\n';
	 std::cout << "truncate_volume_fractions[im_opp-1]= " << 
  	   truncate_volume_fractions[im_opp-1] << '\n';
 	 amrex::Warning("truncate_volume_fractions==1 for im or im_opp");
	}
       } else if (LL==0.0) {
        // do nothing
       } else
        amrex::Error("LL invalid");
      
       if (is_multi_component_evap(freezing_model[indexEXP],
           Tanasawa_or_Schrage_or_Kassemi[indexEXP],LL)==1) {

        int massfrac_id=mass_fraction_id[indexEXP];

        if ((massfrac_id<1)||(massfrac_id>num_species_var))
         amrex::Error("massfrac_id invalid");
        if (LL>0.0) { //evaporation
          spec_material_id_LIQUID[massfrac_id-1]=im_source;
          spec_material_id_AMBIENT[massfrac_id-1]=im_dest;

          if (material_type[im_dest-1]==0) {
	   // do nothing
	  } else
	   amrex::Error("material_type[im_dest-1] invalid for evaporation");

        } else if (LL<0.0) { // condensation
          spec_material_id_LIQUID[massfrac_id-1]=im_dest;
          spec_material_id_AMBIENT[massfrac_id-1]=im_source;

          if (material_type[im_source-1]==0) {
	   // do nothing
	  } else
	   amrex::Error("material_type[im_source-1] invalid for condensation");

        } else
          amrex::Error("LL invalid");

       } else if (is_multi_component_evap(freezing_model[indexEXP],
                  Tanasawa_or_Schrage_or_Kassemi[indexEXP],LL)==0) {
        // do nothing
       } else
        amrex::Error("is_multi_component_evap invalid");

      } // ireverse=0,1
     } //im_opp=im+1..num_materials
    } // im=1..num_materials

    for (int i=0;i<num_species_var;i++) {
     int im=spec_material_id_AMBIENT[i];
     if (im==0) {
      // check nothing
     } else if ((im>=1)&&(im<=num_materials)) {
      if (material_type[im-1]==0) {
       if ((override_density[im-1]==1)|| //rho=rho(T,Y,z)
           (override_density[im-1]==2)) {//Boussinesq approximation
        // do nothing
       } else if (override_density[im-1]==0) {
	// do nothing
       } else {
	std::cout << "override_density==0,1, or 2 allowed if incomp mat.\n";
        std::cout << "num_materials= " << num_materials << '\n';
        std::cout << "im-1=" << im-1 << " override_density[im-1]= " <<
	     override_density[im-1] << '\n';
        amrex::Error("override_density invalid (2)");
       }
      } else if (material_type[im-1]==999) {
       if (override_density[im]!=0)
        amrex::Error("override_density invalid");
      } else if (material_type[im-1]>=1) {
       if (override_density[im]!=0)
        amrex::Error("override_density invalid");
      } else
       amrex::Error("material_type[im-1] invalid");
     } else
      amrex::Error("im invalid 4542");
    } // i=0..num_species_var-1

    mof_error_ordering=0; 
    pp.queryAdd("mof_error_ordering",mof_error_ordering);
    if ((mof_error_ordering!=0)&& //centroid furthest from uncaptured centroid
        (mof_error_ordering!=1))  //smallest MOF error
     amrex::Error("mof_error_ordering invalid");

    mof_ordering.resize(num_materials);

     //fort_mof_ordering_override is declared in: PROB_CPP_PARMS.F90
    fort_mof_ordering_override(
      mof_ordering.dataPtr(),
      &mof_error_ordering,
      FSI_flag.dataPtr());

    pp.queryAdd("mof_ordering",mof_ordering,num_materials);
    for (int i=0;i<num_materials;i++) {
     if ((mof_ordering[i]<0)||
         (mof_ordering[i]>num_materials+1))
      amrex::Error("mof_ordering invalid");
    }

    for (int i=0;i<num_materials;i++) {

     if (visc_coef*viscconst[i]<0.0) {
      amrex::Error("viscosity coefficients invalid");
     } else if (visc_coef*viscconst[i]==0.0) {
      // do nothing
     } else if (visc_coef*viscconst[i]>0.0) {
      // do nothing
     } else
      amrex::Error("viscconst bust");
 
     if (heatviscconst[i]==0.0) {
      // do nothing
     } else if (heatviscconst[i]>0.0) {
      // do nothing
     } else
      amrex::Error("heatviscconst invalid");

     for (int imspec=0;imspec<num_species_var;imspec++) {
      if (speciesviscconst[imspec*num_materials+i]==0.0) {
       // do nothing
      } else if (speciesviscconst[imspec*num_materials+i]>0.0) {
       // do nothing
      } else
       amrex::Error("speciesviscconst invalid");
     } // imspec

    } // i=0..num_materials-1

    if (ParallelDescriptor::IOProcessor()) {

     for (int i=0;i<num_interfaces;i++) {
      std::cout << "i,tension=" << i << ' ' <<
         tension[i] << '\n';
      std::cout << "i,tension_init=" << i << ' ' <<
         tension_init[i] << '\n';
     }

     std::cout << "temperature_source=" << temperature_source << '\n';

     std::cout << "mglib_max_ratio=" << 
        mglib_max_ratio << '\n';
     for (int i=0;i<AMREX_SPACEDIM;i++) {
      std::cout << "i,temperature_source_cen=" << i << ' ' <<
         temperature_source_cen[i] << '\n';
      std::cout << "i,temperature_source_rad=" << i << ' ' <<
         temperature_source_rad[i] << '\n';
     }

     for (int i=0;i<num_interfaces;i++) {
      std::cout << "i= " << i << " denconst_interface_added "  << 
        denconst_interface_added[i] << '\n';
      std::cout << "i= " << i << " viscconst_interface "  << 
        viscconst_interface[i] << '\n';
      std::cout << "i= " << i << " heatviscconst_interface "  << 
        heatviscconst_interface[i] << '\n';
      for (int j=0;j<num_species_var;j++) {
       std::cout << "i= " << i << " j= " << j << 
         " speciesviscconst_interface "  << 
         speciesviscconst_interface[j*num_interfaces+i] << '\n';
      }

     } // i=0 ... num_interfaces-1

     for (int j=0;j<num_species_var;j++) {
      std::cout << " j= " << j << 
         " species_molar_mass "  <<
         species_molar_mass[j] << '\n';
     }  

     std::cout << "CTML_FSI_numsolids " << CTML_FSI_numsolids << '\n';

     std::cout << "mof_error_ordering " << 
      mof_error_ordering << '\n';

     std::cout << "ngrow_make_distance= " << 
      ngrow_make_distance << '\n';
     std::cout << "ngrow_distance= " << 
      ngrow_distance << '\n';
     std::cout << "prescribe_temperature_outflow= " << 
      prescribe_temperature_outflow << '\n';
     std::cout << "solidheat_flag= " << solidheat_flag << '\n';
     std::cout << "truncate_thickness= " << truncate_thickness << '\n';

     for (int i=0;i<num_materials;i++) {
      std::cout << "i= " << i << " compressible_dt_factor= " <<
        compressible_dt_factor[i] << '\n';
     }

     std::cout << "disable_advection= " << disable_advection << '\n';
     std::cout << "disable_pressure_solve= " << disable_pressure_solve << '\n';
     std::cout << "nparts (im_solid_map.size())= " << 
      im_solid_map.size() << '\n';
     std::cout << "Solid_State_Type= " << Solid_State_Type << '\n';
     std::cout << "Tensor_Type= " << Tensor_Type << '\n';
     std::cout << "TensorX_Type= " << TensorX_Type << '\n';
     std::cout << "TensorY_Type= " << TensorY_Type << '\n';
     std::cout << "TensorZ_Type= " << TensorZ_Type << '\n';
     std::cout << "NUM_CELL_ELASTIC= " << NUM_CELL_ELASTIC << '\n';
     std::cout << "NUM_STATE_TYPE= " << NUM_STATE_TYPE << '\n';

     std::cout << "angular_velocity= " << angular_velocity << '\n';
     std::cout << "centrifugal_force_factor= " << 
       centrifugal_force_factor << '\n';

     std::cout << "uncoupled_viscosity= " << uncoupled_viscosity << '\n';

     std::cout << "pressure_error_flag=" << pressure_error_flag << '\n';

     std::cout << "initial_temperature_diffuse_duration=" << 
      initial_temperature_diffuse_duration << '\n';

     std::cout<<"particle_nsubdivide="<<particle_nsubdivide<<'\n';
     std::cout << "particle_max_per_nsubdivide= " <<
        particle_max_per_nsubdivide << '\n';

     for (int i=0;i<num_materials;i++) {
      std::cout << "mof_ordering i= " << i << ' ' <<
        mof_ordering[i] << '\n';

      std::cout << "truncate_volume_fractions i= " << i << ' ' <<
        truncate_volume_fractions[i] << '\n';

      std::cout << "viscosity_state_model i= " << i << ' ' <<
        viscosity_state_model[i] << '\n';
      std::cout << "viscoelastic_model i= " << i << ' ' <<
        viscoelastic_model[i] << '\n';
      std::cout << "les_model i= " << i << ' ' <<
        les_model[i] << '\n';
      std::cout << "shock_timestep i=" << i << " " << 
          shock_timestep[i] << '\n';
      std::cout << "material_type i=" << i << " " << material_type[i] << '\n';
      std::cout << "material_type_evap i=" << i << " " << 
	      material_type_evap[i] << '\n';
      std::cout << "material_type_lowmach i=" << i << " " << 
	      material_type_lowmach[i] << '\n';
      std::cout << "material_type_visual i=" << i << " " << 
	      material_type_visual[i] << '\n';
      std::cout << "pressure_error_cutoff i=" << i << " " << 
        pressure_error_cutoff[i] << '\n';
      std::cout << "vorterr i=" << i << " " << 
        vorterr[i] << '\n';
      std::cout << "temperature_error_cutoff i=" << i << " " << 
        temperature_error_cutoff[i] << '\n';
      std::cout << "stiffPINF i=" << i << " " << stiffPINF[i] << '\n';
      std::cout << "stiffCP i=" << i << " " << stiffCP[i] << '\n';
      std::cout << "stiffCV i=" << i << " " << stiffCV[i] << '\n';
      std::cout << "stiffGAMMA i=" << i << " " << stiffGAMMA[i] << '\n';
      std::cout << "denconst i=" << i << " " << denconst[i] << '\n';
      std::cout << "density_floor i=" << i << " " << density_floor[i] << '\n';
      std::cout << "density_ceiling i="<<i<<" "<< density_ceiling[i] << '\n';
      std::cout << "molar_mass i="<<i<<" "<< 
        molar_mass[i] << '\n';
      std::cout << "tempconst i=" << i << " " << tempconst[i] << '\n';
      std::cout << "initial_temperature i=" << i << " " << 
        initial_temperature[i] << '\n';
      std::cout << "tempcutoff i=" << i << " " << tempcutoff[i] << '\n';
      std::cout << "tempcutoffmax i=" << i << " " << tempcutoffmax[i] << '\n';
      std::cout << "DrhoDT i=" << i << " " << DrhoDT[i] << '\n';
      std::cout << "override_density i=" << i << " " << 
         override_density[i] << '\n';
      std::cout << "viscconst i=" << i << "  " << viscconst[i] << '\n';

      std::cout << "viscconst_eddy_wall i=" <<i<<"  "<<
	      viscconst_eddy_wall[i]<<'\n';
      std::cout << "viscconst_eddy_bulk i=" <<i<<"  "<<
	      viscconst_eddy_bulk[i]<<'\n';
      std::cout << "heatviscconst_eddy_wall i=" <<i<<"  "<<
	      heatviscconst_eddy_wall[i]<<'\n';
      std::cout << "heatviscconst_eddy_bulk i=" <<i<<"  "<<
	      heatviscconst_eddy_bulk[i]<<'\n';

      std::cout << "heatflux_factor i=" << i << "  " << 
          heatflux_factor[i] << '\n';
      std::cout << "heatviscconst i=" << i << "  " << 
          heatviscconst[i] << '\n';
      std::cout << "prerecalesce_viscconst i=" << i << "  " << 
         prerecalesce_viscconst[i] << '\n';
      std::cout << "prerecalesce_heatviscconst i=" << i << "  " << 
         prerecalesce_heatviscconst[i] << '\n';
      std::cout << "prerecalesce_stiffCP i=" << i << "  " << 
         prerecalesce_stiffCP[i] << '\n';
      std::cout << "prerecalesce_stiffCV i=" << i << "  " << 
         prerecalesce_stiffCV[i] << '\n';
     }  // i=0,..,num_materials

     for (int i=0;i<num_species_var*num_materials;i++) {
      std::cout << "speciesviscconst i=" << i << "  " << 
          speciesviscconst[i] << '\n';
      std::cout << "speciesconst i=" << i << "  " << 
          speciesconst[i] << '\n';
      std::cout << "speciesreactionrate i=" << i << "  " << 
          speciesreactionrate[i] << '\n';
     }

     std::cout << "stokes_flow= " << stokes_flow << '\n';
     std::cout << "cancel_advection= " << cancel_advection << '\n';

     std::cout << "is_phasechange= " << is_phasechange << '\n';

     std::cout << "perturbation_on_restart " << perturbation_on_restart << '\n';
     std::cout << "perturbation_mode " << perturbation_mode << '\n';
     std::cout << "perturbation_eps_temp " << perturbation_eps_temp << '\n';
     std::cout << "perturbation_eps_vel " << perturbation_eps_vel << '\n';

     std::cout << "custom_nucleation_model " << 
       custom_nucleation_model << '\n';

     std::cout << "FD_curv_interp " << FD_curv_interp << '\n';
     std::cout << "vof_height_function " << vof_height_function << '\n';

     std::cout << "hydrate flag " << hydrate_flag << '\n';
     std::cout << "nucleation_period= " << nucleation_period << '\n';
     std::cout << "nucleation_init_time= " << nucleation_init_time << '\n';
     std::cout << "n_sites= " << n_sites << '\n';
     if (n_sites>0) {
      for (int i=0;i<pos_sites.size();i++) {
       std::cout << "i, pos_sites= " << i << ' ' << pos_sites[i] << '\n';
      }
     }
     std::cout << "pos_sites_random_flag= " << pos_sites_random_flag << '\n';
    
     for (int i=0;i<num_materials;i++) {
      std::cout << "microlayer_substrate i=" << i << "  " << 
       microlayer_substrate[i] << '\n';
      std::cout << "microlayer_angle i=" << i << "  " << 
       microlayer_angle[i] << '\n';
      std::cout << "microlayer_size i=" << i << "  " << 
       microlayer_size[i] << '\n';
      std::cout << "macrolayer_size i=" << i << "  " << 
       macrolayer_size[i] << '\n';
      std::cout << "max_contact_line_size i=" << i << "  " << 
       max_contact_line_size[i] << '\n';
      std::cout << "microlayer_temperature_substrate i=" << i << "  " << 
       microlayer_temperature_substrate[i] << '\n';
      std::cout << "thermal_microlayer_size i=" << i << "  " << 
       thermal_microlayer_size[i] << '\n';
      std::cout << "shear_microlayer_size i=" << i << "  " << 
       shear_microlayer_size[i] << '\n';
      std::cout << "buoyancy_microlayer_size i=" << i << "  " << 
       buoyancy_microlayer_size[i] << '\n';
      std::cout << "phasechange_microlayer_size i=" << i << "  " << 
       phasechange_microlayer_size[i] << '\n';
     } // i=0..num_materials-1

     for (int i=0;i<num_species_var;i++) {
      std::cout << "spec_material_id_LIQUID i= " << i << " " <<
       spec_material_id_LIQUID[i] << '\n';
      std::cout << "spec_material_id_AMBIENT i= " << i << " " <<
       spec_material_id_AMBIENT[i] << '\n';
     }

     for (int i=0;i<2*AMREX_SPACEDIM;i++) {
      std::cout << "i= " << i << " outflow_velocity_buffer_size= " <<
       outflow_velocity_buffer_size[i] << '\n';
     }

     std::cout << "R_Palmore_Desjardins " << R_Palmore_Desjardins << '\n';

     std::cout << "unscaled_min_curvature_radius=" << 
	      unscaled_min_curvature_radius << '\n';

     for (int i=0;i<num_interfaces;i++) {
      std::cout << "hardwire_T_gamma i=" << i << "  " << 
       hardwire_T_gamma[i] << '\n';
      std::cout << "hardwire_T_gamma i+num_interfaces=" << 
       i+num_interfaces << "  " << 
       hardwire_T_gamma[i+num_interfaces] << '\n';

      std::cout << "hardwire_Y_gamma i=" << i << "  " << 
       hardwire_Y_gamma[i] << '\n';
      std::cout << "hardwire_Y_gamma i+num_interfaces=" << 
       i+num_interfaces << "  " << 
       hardwire_Y_gamma[i+num_interfaces] << '\n';

      std::cout << "accommodation_coefficient i=" << i << "  " << 
       accommodation_coefficient[i] << '\n';
      std::cout << "accommodation_coefficient i+num_interfaces=" << i+num_interfaces << "  " << 
       accommodation_coefficient[i+num_interfaces] << '\n';

      std::cout << "reference_pressure i=" << i << "  " << 
       reference_pressure[i] << '\n';
      std::cout << "reference_pressure i+num_interfaces=" << i+num_interfaces << "  " << 
       reference_pressure[i+num_interfaces] << '\n';

      std::cout << "saturation_temp i=" << i << "  " << 
       saturation_temp[i] << '\n';
      std::cout << "saturation_temp i+num_interfaces=" << i+num_interfaces << "  " << 
       saturation_temp[i+num_interfaces] << '\n';

      std::cout << "saturation_temp_curv i=" << i << "  " << 
       saturation_temp_curv[i] << '\n';
      std::cout << "saturation_temp_curv i+num_interfaces=" << i+num_interfaces << "  " << 
       saturation_temp_curv[i+num_interfaces] << '\n';

      std::cout << "saturation_temp_vel i=" << i << "  " << 
       saturation_temp_vel[i] << '\n';
      std::cout << "saturation_temp_vel i+num_interfaces=" << i+num_interfaces << "  " << 
       saturation_temp_vel[i+num_interfaces] << '\n';

      std::cout << "saturation_temp_min i=" << i << "  " << 
       saturation_temp_min[i] << '\n';
      std::cout << "saturation_temp_max i+num_interfaces=" << i+num_interfaces << "  " << 
       saturation_temp_max[i+num_interfaces] << '\n';

      std::cout << "nucleation_temp i=" << i << "  " << 
       nucleation_temp[i] << '\n';
      std::cout << "nucleation_temp i+num_interfaces=" << i+num_interfaces << "  " << 
       nucleation_temp[i+num_interfaces] << '\n';

      std::cout << "nucleation_pressure i=" << i << "  " << 
       nucleation_pressure[i] << '\n';
      std::cout << "nucleation_pressure i+num_interfaces=" << i+num_interfaces << "  " << 
       nucleation_pressure[i+num_interfaces] << '\n';

      std::cout << "nucleation_pmg i=" << i << "  " << 
       nucleation_pmg[i] << '\n';
      std::cout << "nucleation_pmg i+num_interfaces=" << i+num_interfaces << "  " << 
       nucleation_pmg[i+num_interfaces] << '\n';

      std::cout << "nucleation_mach i=" << i << "  " << 
       nucleation_mach[i] << '\n';
      std::cout << "nucleation_mach i+num_interfaces=" << 
	     i+num_interfaces << "  " << 
       nucleation_mach[i+num_interfaces] << '\n';

      std::cout << "material_type_interface i=" << i << "  " << 
       material_type_interface[i] << '\n';

      std::cout << "latent_heat i=" << i << "  " << 
       latent_heat[i] << '\n';
      std::cout << "latent_heat i+num_interfaces=" << 
       i+num_interfaces << "  " << latent_heat[i+num_interfaces] << '\n';

      std::cout << "latent_heat_slope i=" << i << "  " << 
       latent_heat_slope[i] << '\n';

      std::cout << "latent_heat_slope i+num_interfaces=" << 
       i+num_interfaces << "  " << 
       latent_heat_slope[i+num_interfaces] << '\n';

      if (latent_heat_slope[i]<=0.0) {
       //do nothing
      } else 
       amrex::Error("need latent_heat_slope[i]<=0.0)");

      if (latent_heat_slope[i+num_interfaces]<=0.0) {
       //do nothing
      } else 
       amrex::Error("need latent_heat_slope[i+num_interfaces]<=0.0)");

      std::cout << "latent_heat_T0 i=" << i << "  " << 
       latent_heat_T0[i] << '\n';
      std::cout << "latent_heat_T0 i+num_interfaces=" << 
       i+num_interfaces << "  " << latent_heat_T0[i+num_interfaces] << '\n';

      std::cout << "latent_heat_min i=" << i << "  " << 
       latent_heat_min[i] << '\n';
      std::cout << "latent_heat_min i+num_interfaces=" << 
       i+num_interfaces << "  " << latent_heat_min[i+num_interfaces] << '\n';

      std::cout << "reaction_rate i=" << i << "  " << 
       reaction_rate[i] << '\n';
      std::cout << "reaction_rate i+num_interfaces=" << 
       i+num_interfaces << "  " << reaction_rate[i+num_interfaces] << '\n';

      std::cout << "freezing_model i=" << i << "  " << 
       freezing_model[i] << '\n';
      std::cout << "freezing_model i+num_interfaces=" << 
       i+num_interfaces << "  " << freezing_model[i+num_interfaces] << '\n';
      std::cout << "Tanasawa_or_Schrage_or_Kassemi i=" << i << "  " << 
       Tanasawa_or_Schrage_or_Kassemi[i] << '\n';
      std::cout << "Tanasawa_or_Schrage_or_Kassemi i+num_interfaces=" << 
       i+num_interfaces << "  " << 
       Tanasawa_or_Schrage_or_Kassemi[i+num_interfaces] << '\n';
      std::cout << "mass_fraction_id i=" << i << "  " << 
       mass_fraction_id[i] << '\n';
      std::cout << "mass_fraction_id i+num_interfaces=" << 
       i+num_interfaces << "  " << mass_fraction_id[i+num_interfaces] << '\n';
      std::cout << "distribute_from_target i=" << i << "  " << 
       distribute_from_target[i] << '\n';
      std::cout << "distribute_from_target i+num_interfaces=" << 
       i+num_interfaces << "  " << 
       distribute_from_target[i+num_interfaces] << '\n';

      std::cout << "tension i=" << i << "  " << tension[i] << '\n';
      std::cout << "tension_slope i=" << i << "  " << tension_slope[i] << '\n';

      if (tension_slope[i]<=0.0) {
       // do nothing
      } else
       amrex::Error("tension_slope must be non-positive(2)");

      std::cout << "tension_T0 i=" << i << "  " << tension_T0[i] << '\n';
      std::cout << "tension_min i=" << i << "  " << tension_min[i] << '\n';
      std::cout << "initial cap_wave_speed i=" << i << "  " << 
        cap_wave_speed[i] << '\n';
      std::cout << "prefreeze_tension i=" << i << "  " << 
       prefreeze_tension[i] << '\n';
      std::cout << "distribute_mdot_evenly i=" << i << "  " << 
       distribute_mdot_evenly[i] << '\n';
      std::cout << "distribute_mdot_evenly i+num_interfaces=" << 
       i+num_interfaces << "  " << 
       distribute_mdot_evenly[i+num_interfaces] << '\n';
      std::cout << "constant_volume_mdot  i=" << i << "  " << 
       constant_volume_mdot[i] << '\n';
      std::cout << "constant_volume_mdot  i+num_interfaces=" << 
       i+num_interfaces << "  " << 
       constant_volume_mdot[i+num_interfaces] << '\n';
     }  // i=0..num_interfaces-1

     for (int i=0;i<num_materials;i++) {
      std::cout << "rigid_fraction_id i=" << i << "  " << 
       rigid_fraction_id[i] << '\n';
      std::cout << "constant_density_all_time i=" << i << "  " << 
       constant_density_all_time[i] << '\n';
      std::cout << "cavitation_pressure i=" << i << "  " << 
       cavitation_pressure[i] << '\n';
      std::cout << "cavitation_vapor_density i=" << i << "  " << 
       cavitation_vapor_density[i] << '\n';
      std::cout << "cavitation_tension i=" << i << "  " << 
       cavitation_tension[i] << '\n';
     } // i=0..num_materials-1

     std::cout << "pressure_select_criterion " << 
       pressure_select_criterion << '\n';

     std::cout << "num_materials_viscoelastic " << 
        num_materials_viscoelastic << '\n';
     std::cout << "num_species_var " << num_species_var << '\n';
     std::cout << "num_materials " << num_materials << '\n';
     std::cout << "num_interfaces " << num_interfaces << '\n';
     std::cout << "MOFITERMAX= " << 
	     MOFITERMAX << '\n';
     std::cout << "MOFITERMAX_AFTER_PREDICT= " << 
	     MOFITERMAX_AFTER_PREDICT << '\n';
     std::cout << "MOF_DEBUG_RECON= " << MOF_DEBUG_RECON << '\n';
     std::cout << "MOF_TURN_OFF_LS= " << MOF_TURN_OFF_LS << '\n';

     std::cout << "post_init_pressure_solve " << 
       post_init_pressure_solve << '\n';

     std::cout << "projection_pressure_scale " << 
       projection_pressure_scale << '\n';
     std::cout << "projection_velocity_scale " << 
       projection_velocity_scale << '\n';

      //in AMReX_MemPool.cpp and AMReX_FArrayBox.cpp: 
      // ParmParse pp("fab");
      // pp.queryAdd("init_snan", init_snan);
#if defined(AMREX_DEBUG) || defined(AMREX_TESTING)
     bool init_snan  = true;
#else
     bool init_snan  = false;
#endif
     ParmParse ppfab("fab");
     ppfab.queryAdd("init_snan",init_snan);

     std::cout << "init_snan= " << init_snan << '\n';

     int do_initval=FArrayBox::get_do_initval();
     std::cout << "do_initval= " << do_initval << '\n';
     Real initval=FArrayBox::get_initval();
     std::cout << "initval= " << initval << '\n';

       // SFC == Space Filling Curve
     DistributionMapping::Strategy local_strategy=
         DistributionMapping::strategy();
     if (local_strategy==DistributionMapping::ROUNDROBIN) {
      std::cout << "DistributionMapping::strategy=ROUNDROBIN\n";
     } else if (local_strategy==DistributionMapping::KNAPSACK) {
      std::cout << "DistributionMapping::strategy=KNAPSACK\n";
     } else if (local_strategy==DistributionMapping::SFC) {
      std::cout << "DistributionMapping::strategy=SFC\n";
     } else if (local_strategy==DistributionMapping::RRSFC) {
      std::cout << "DistributionMapping::strategy=RRSFC\n";
     } else
      amrex::Error("local_strategy invalid");

     std::cout << "krylov_subspace_max_num_outer_iter " << 
       krylov_subspace_max_num_outer_iter << '\n';
     std::cout << "slipcoeff " << slipcoeff << '\n';

     std::cout << "EILE_flag " << EILE_flag << '\n';

     std::cout << "ractive " << ractive << '\n';
     std::cout << "ractivex " << ractivex << '\n';
     std::cout << "ractivey " << ractivey << '\n';
     std::cout << "ractivez " << ractivez << '\n';
     std::cout << "wait_time " << wait_time << '\n';

     std::cout << "multilevel_maxcycle " << multilevel_maxcycle << '\n';

     std::cout << "multilevel_restart_period " << 
	     multilevel_restart_period << '\n';
     if (cg_restart_period>=0) {
      std::cout << "cg.restart_period " << 
	     cg_restart_period << '\n';
     }
     if (cg_maxiter>=0) {
      std::cout << "cg.maxiter " << 
	     cg_maxiter << '\n';
     }

     std::cout << "mac.mac_abs_tol " <<mac_abs_tol<< '\n';
     std::cout << "mac.visc_abs_tol " <<visc_abs_tol<< '\n';
     std::cout << "mac.thermal_abs_tol " <<thermal_abs_tol<< '\n';
     std::cout << "project_solver_type " <<project_solver_type<< '\n';
     std::cout << "initial_cg_cycles " <<initial_cg_cycles<< '\n';
     std::cout << "initial_project_cycles " <<initial_project_cycles<< '\n';
     std::cout << "initial_viscosity_cycles " <<initial_viscosity_cycles<< '\n';
     std::cout << "initial_thermal_cycles " <<initial_thermal_cycles<< '\n';
     std::cout << "visual_tessellate_vfrac " << visual_tessellate_vfrac << '\n';

     std::cout << "step_through_data " << step_through_data << '\n';

     std::cout << "visual_revolve " << visual_revolve << '\n';
     std::cout << "visual_output_raw_State_Type " << 
	     visual_output_raw_State_Type << '\n';
     std::cout << "visual_output_raw_mac_Type " << 
	     visual_output_raw_mac_Type << '\n';
     std::cout << "visual_phase_change_plot_int " << 
	     visual_phase_change_plot_int << '\n';
     std::cout << "visual_buoyancy_plot_int " << 
	     visual_buoyancy_plot_int << '\n';
     std::cout << "visual_divergence_plot_int " << 
	     visual_divergence_plot_int << '\n';
     std::cout << "visual_WALLVEL_plot_int " << 
	     visual_WALLVEL_plot_int << '\n';
     std::cout << "visual_drag_plot_int " << 
	     visual_drag_plot_int << '\n';
     std::cout << "visual_nddata_format " << 
	     visual_nddata_format << '\n';

     std::cout << "visual_compare " << visual_compare << '\n';
     for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
      std::cout << "dir,visual_ncell " << dir << ' ' << 
       visual_ncell[dir] << '\n';
     }
     std::cout << "change_max=" << change_max << '\n';
     std::cout << "change_max_init=" << change_max_init << '\n';
     std::cout << "fixed_dt=" << fixed_dt << '\n';
     std::cout << "fixed_dt_init=" << fixed_dt_init << '\n';
     std::cout << "fixed_dt_velocity=" << fixed_dt_velocity << '\n';
     std::cout << "min_velocity_for_dt=" << min_velocity_for_dt << '\n';
     std::cout << "min_stefan_velocity_for_dt=" << 
        min_stefan_velocity_for_dt << '\n';
     std::cout << "dt_max=" << dt_max << '\n';
     std::cout << "minimum_relative_error=" << minimum_relative_error << '\n';
     std::cout << "diffusion_minimum_relative_error=" << 
      diffusion_minimum_relative_error << '\n';
     std::cout << "num_divu_outer_sweeps=" << num_divu_outer_sweeps << '\n';

     std::cout << "ns_max_level= " << ns_max_level << '\n';
     for (int ilev=0;ilev<ns_max_grid_size.size();ilev++) {
      std::cout << "ilev, max_grid_size " << ilev << ' ' <<
       ns_max_grid_size[ilev] << '\n';
     }
     for (int ilev=0;ilev<ns_max_level;ilev++) {
      std::cout << "ilev, ns_n_error_buf " << ilev << ' ' <<
       ns_n_error_buf[ilev] << '\n';
     }
    }  // if IO processor

    set_local_tolerances(SOLVETYPE_PRES);
    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "project_option_tol(regular project)= " << 
	     SOLVETYPE_PRES << '\n';
     std::cout << "save_mac_abs_tol= " << save_mac_abs_tol << '\n';
     std::cout << "save_atol_b= " << save_atol_b << '\n';
     std::cout << "save_min_rel_error= " << save_min_rel_error << '\n';
    }

    set_local_tolerances(SOLVETYPE_HEAT);
    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "project_option_tol(thermal)= " << 
	     SOLVETYPE_HEAT << '\n';
     std::cout << "save_mac_abs_tol= " << save_mac_abs_tol << '\n';
     std::cout << "save_atol_b= " << save_atol_b << '\n';
     std::cout << "save_min_rel_error= " << save_min_rel_error << '\n';
    }

    set_local_tolerances(SOLVETYPE_VISC);
    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "project_option_tol(viscosity)= " << 
	     SOLVETYPE_VISC << '\n';
     std::cout << "save_mac_abs_tol= " << save_mac_abs_tol << '\n';
     std::cout << "save_atol_b= " << save_atol_b << '\n';
     std::cout << "save_min_rel_error= " << save_min_rel_error << '\n';
    }
     for (int im=1;im<=num_materials;im++) {

      if (is_ice_or_FSI_rigid_material(&im)==1) {

       if (enable_spectral==0) {
        // do nothing
       } else if (enable_spectral==1) {
        //having just SOLVETYPE_PRESGRAVITY low order will not
        //adversely effect the numerical dissipation of an 
        //otherwise spectrally accurate method (only the order
        //will be adversely effected)
       } else
        amrex::Error("enable_spectral invalid");

      } // is_ice_or_FSI_rigid_material==1

    } // im=1..num_materials

    if (enable_spectral==1) {

     // do nothing

    } else if (enable_spectral==0) {

     // do nothing

    } else {
     amrex::Error("enable_spectral invalid");
    }
     
    if (some_materials_compressible()==1) {

     if (num_divu_outer_sweeps<2)
      amrex::Error("num_divu_outer_sweeps>=2 for comp materials");

    } else if (some_materials_compressible()==0) {
     // do nothing
    } else
     amrex::Error("compressible flag bust");

    setup_integrated_quantities();

} // end subroutine read_params()

NavierStokes::NavierStokes ()
{
    Geometry_setup();
}

// constructor
NavierStokes::NavierStokes (AmrCore&        papa,
                            int             lev,
                            const Geometry& level_geom,
                            const BoxArray& bl,
                            const DistributionMapping& dmap_in,
                            Real            time)
    :
    AmrLevel(papa,lev,level_geom,bl,dmap_in,time)
{
    Geometry_setup();
}

NavierStokes::~NavierStokes ()
{
    Geometry_cleanup();

    for (int i=0;i<MAX_NUM_LOCAL_MF;i++)
     if (localMF_grow[i]==-1) {
      // do nothing
     } else {
      std::cout << "i= " << i << " localMF_grow= " <<
       localMF_grow[i] << '\n';
      amrex::Error("localMF_grow invalid");
     }

}

int NavierStokes::ns_is_rigid(int im) {

 if ((im<0)|(im>=num_materials))
  amrex::Error("im invalid50 (ns_is_rigid)");
 
 int imp1=im+1;
 int local_flag=fort_is_rigid_base(&FSI_flag[im],&imp1);

 return local_flag;

} // subroutine ns_is_rigid


int NavierStokes::ns_is_lag_part(int im) {

 if ((im<0)|(im>=num_materials))
  amrex::Error("im invalid50 (ns_is_lag_part)");

 int imp1=im+1;
 int local_flag=fort_is_lag_part_base(&FSI_flag[im],&imp1);

 return local_flag;

} // subroutine ns_is_lag_part

int NavierStokes::is_GFM_freezing_model(int loc_freezing_model) {

 if ((loc_freezing_model==0)||   //fully saturated
     (loc_freezing_model==5)||   //stefan model evap or condensation
     (loc_freezing_model==6)) {  //Palmore and Desjardins
  return 1;
 } else if (is_valid_freezing_model(loc_freezing_model)==1) {
  return 0;
 } else {
  amrex::Error("freezing_model bust");
  return 0;
 }

} //end function is_GFM_freezing_model(int freezing_model) 

int NavierStokes::is_hydrate_freezing_model(int loc_freezing_model) {

 if (loc_freezing_model==2) {
  return 1;
 } else if (is_valid_freezing_model(loc_freezing_model)) {
  return 0;
 } else {
  amrex::Error("freezing_model invalid");
  return 0;
 }
} // end function is_hydrate_freezing_model(int freezing_model)

int NavierStokes::is_valid_freezing_model(int loc_freezing_model) {

 if ((loc_freezing_model==5)|| //Stefan model evaporation or condensation
     (loc_freezing_model==6)|| //Palmore and Desjardins
     (loc_freezing_model==7)) {//cavitation
  return 1;
 } else if ((loc_freezing_model==0)|| //Energy jump model
            (loc_freezing_model==1)|| //source term
            (loc_freezing_model==2)|| //hydrate
            (loc_freezing_model==3)) {//wildfire
  return 1;
 } else {
  amrex::Error("loc_freezing_model invalid");
  return 0;
 }

} // end function is_valid_freezing_model


int NavierStokes::is_multi_component_evap(int loc_freezing_model,
           int loc_evap_flag,Real loc_latent_heat) {

 if (loc_latent_heat==0.0) {
  return 0;
 } else if (loc_latent_heat!=0.0) {

  if ((loc_freezing_model==5)|| //Stefan model evaporation or condensation
      (loc_freezing_model==6)|| //Palmore and Desjardins
      (loc_freezing_model==7)) {//cavitation

   if (loc_evap_flag==0) { //Palmore and Desjardins
    return 1;
   } else if ((loc_evap_flag==1)|| //Tanasawa
              (loc_evap_flag==2)|| //Schrage
              (loc_evap_flag==3)) { //Kassemi
    return 0;
   } else {
    amrex::Error("loc_evap_flag invalid");
    return 0;
   }

  } else if (loc_freezing_model==2) { //hydrate
   return 1;
  } else if ((loc_freezing_model==0)|| //Energy jump model
             (loc_freezing_model==1)|| //source term
             (loc_freezing_model==3)) {//wildfire
   return 0;
  } else {
   amrex::Error("loc_freezing_model invalid");
   return 0;
  }

 } else {
  amrex::Error("loc_latent_heat invalid");
  return 0;
 }

} // end function is_multi_component_evap

int NavierStokes::project_option_is_valid(int project_option) {

 if (project_option_momeqn(project_option)==1) {
  return 1;
 } else if (project_option_momeqn(project_option)==0) {
  return 1;
 } else {
  amrex::Error("project_option not valid");
  return 0;
 }

} // end function project_option_is_valid(project_option)

int NavierStokes::project_option_momeqn(int project_option) {

 if ((project_option==SOLVETYPE_PRES)|| 
     (project_option==SOLVETYPE_INITPROJ)||  
     (project_option==SOLVETYPE_PRESGRAVITY)|| 
     (project_option==SOLVETYPE_PRESEXTRAP)|| 
     (project_option==SOLVETYPE_VISC)) {
  return 1;
 } else if ((project_option==SOLVETYPE_HEAT)||  
            ((project_option>=SOLVETYPE_SPEC)&& 
             (project_option<SOLVETYPE_SPEC+num_species_var))) {
  return 0;
 } else {
  amrex::Error("project_option invalid project_option_momeqn");
  return 0;
 }

}  // static function project_option_momeqn(project_option)


int NavierStokes::project_option_singular_possible(int project_option) {

 int local_flag=project_option_singular_possibleF(&project_option);
 return local_flag;

} // end function project_option_singular_possible(project_option)


int NavierStokes::project_option_olddata_needed(int project_option) {

 if ((project_option==SOLVETYPE_PRES)|| 
     (project_option==SOLVETYPE_INITPROJ)|| 
     (project_option==SOLVETYPE_PRESGRAVITY)||  
     (project_option==SOLVETYPE_PRESEXTRAP)) {
  return 0;
 } else if ((project_option==SOLVETYPE_HEAT)|| 
   	    (project_option==SOLVETYPE_VISC)|| 
            ((project_option>=SOLVETYPE_SPEC)&&
	     (project_option<SOLVETYPE_SPEC+num_species_var))) {
  return 1;
 } else {
  amrex::Error("project_option invalid function project_option_olddata_needed");
  return 0;
 }

} // end function project_option_olddata_needed()


int NavierStokes::project_option_pressure(int project_option) {

 if ((project_option==SOLVETYPE_PRES)||
     (project_option==SOLVETYPE_PRESGRAVITY)||
     (project_option==SOLVETYPE_INITPROJ)||
     (project_option==SOLVETYPE_PRESEXTRAP)) {
  return 1;
 } else if ((project_option==SOLVETYPE_HEAT)|| 
	    (project_option==SOLVETYPE_VISC)|| 
            ((project_option>=SOLVETYPE_SPEC)&&
             (project_option<SOLVETYPE_SPEC+num_species_var))) {
  return 0;
 } else {
  amrex::Error("project_option invalid project_option_pressure");
  return 0;
 }

}  // static function project_option_pressure(project_option)

int NavierStokes::project_option_needs_scaling(int project_option) {

 int local_flag=project_option_needs_scalingF(&project_option);
 return local_flag;

}  // static function project_option_needs_scaling(project_option)


int NavierStokes::project_option_projection(int project_option) {

 int local_flag=project_option_projectionF(&project_option);
 return local_flag;

}  // static function project_option_projection(project_option)



// getState_list needs scomp,ncomp
void
NavierStokes::get_mm_scomp_solver(
  int num_materials_combine,
  int project_option,
  int& state_index,
  Vector<int>& scomp,
  Vector<int>& ncomp,
  int& ncomp_check) {

 int nsolve=1;
 int nlist=1;

 if ((project_option==SOLVETYPE_PRES)||   
     (project_option==SOLVETYPE_PRESGRAVITY)||   
     (project_option==SOLVETYPE_INITPROJ)) { 
  nsolve=1;
  nlist=1;
  if (num_materials_combine!=1)
   amrex::Error("num_materials_combine invalid");
 } else if (project_option==SOLVETYPE_PRESEXTRAP) { 
  nsolve=1;
  nlist=1;
  if (num_materials_combine!=1)
   amrex::Error("num_materials_combine invalid");
 } else if (project_option==SOLVETYPE_HEAT) { 
  nsolve=1;
  nlist=num_materials_combine;
  if ((num_materials_combine!=1)&&
      (num_materials_combine!=num_materials))
   amrex::Error("num_materials_combine invalid");
 } else if ((project_option>=SOLVETYPE_SPEC)&&
            (project_option<SOLVETYPE_SPEC+num_species_var)) { 
  nsolve=1;
  nlist=num_materials_combine;
  if ((num_materials_combine!=1)&& 
      (num_materials_combine!=num_materials))
   amrex::Error("num_materials_combine invalid");
 } else if (project_option==SOLVETYPE_VISC) { 
  nsolve=AMREX_SPACEDIM;
  nlist=1;
  if (num_materials_combine!=1)
   amrex::Error("num_materials_combine invalid");
 } else
  amrex::Error("project_option invalid1");

 if ((num_materials_combine!=1)&&
     (num_materials_combine!=num_materials)) 
  amrex::Error("num_materials_combine invalid");

 scomp.resize(nlist);
 ncomp.resize(nlist);

 int nsolveMM=nsolve*num_materials_combine;

 if (project_option_pressure(project_option)==1) {

  if (num_materials_combine==1) { 
   scomp[0]=STATECOMP_PRES;
   ncomp[0]=nsolveMM; 
   state_index=State_Type;
  } else
   amrex::Error("num_materials_combine invalid");

 } else if (project_option==SOLVETYPE_HEAT) { 

  // u,v,w,p,den1,T1,...
  for (int im=0;im<nlist;im++) {
   scomp[im]=STATECOMP_STATES+im*num_state_material+ENUM_TEMPERATUREVAR;
   ncomp[im]=1;
  }
  state_index=State_Type;

 } else if (project_option==SOLVETYPE_VISC) { 

  if (num_materials_combine==1) { 
   scomp[0]=STATECOMP_VEL;
   ncomp[0]=nsolveMM; 
   state_index=State_Type;
  } else
   amrex::Error("num_materials_combine invalid");

 } else if ((project_option>=SOLVETYPE_SPEC)&&
            (project_option<SOLVETYPE_SPEC+num_species_var)) { 

  for (int im=0;im<nlist;im++) {

   if (ENUM_SPECIESVAR==num_state_base) {
    // do nothing
   } else
    amrex::Error("ENUM_SPECIES_VAR invalid");

   scomp[im]=STATECOMP_STATES+
     im*num_state_material+num_state_base+project_option-SOLVETYPE_SPEC;
   ncomp[im]=1;
  }
  state_index=State_Type;

 } else
  amrex::Error("project_option invalid get_mm_scomp_solver");

 ncomp_check=0;
 for (int im=0;im<nlist;im++)
  ncomp_check+=ncomp[im];

 if (ncomp_check!=nsolveMM)
  amrex::Error("ncomp_check invalid");

} // end subroutine get_mm_scomp_solver

void
NavierStokes::zero_independent_vel(int project_option,int idx,int nsolve) {

 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level corrupt");

 if ((nsolve!=1)&&(nsolve!=AMREX_SPACEDIM))
  amrex::Error("nsolve invalid36");

 if (project_option_momeqn(project_option)==1) {

  // do nothing
  
 } else if (project_option_momeqn(project_option)==0) {

  // do nothing

 } else
  amrex::Error("project_option_momeqn(project_option) invalid");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  if (localMF[idx+dir]->nComp()!=nsolve)
   amrex::Error("localMF[idx+dir] has invalid ncomp");
  setVal_localMF(idx+dir,0.0,0,nsolve,0);
 } // dir

} // end subroutine zero_independent_vel

// u,v,w,p,den1,T1,...,den2,T2,...
void
NavierStokes::zero_independent_variable(int project_option,int nsolve) {

 if (num_state_base!=2) 
  amrex::Error("num_state_base invalid");

 if ((nsolve!=1)&&
     (nsolve!=AMREX_SPACEDIM))
  amrex::Error("nsolve invalid36");

 if (project_option_momeqn(project_option)==1) {
  // do nothing
 } else if (project_option_momeqn(project_option)==0) {
  // do nothing
 } else
  amrex::Error("project_option invalid3");

 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level corrupt");

 Vector<int> scomp;
 Vector<int> ncomp;
 int ncomp_check;
 int state_index;
  //num_materials_combine=1
 get_mm_scomp_solver(
   1,
   project_option,
   state_index,
   scomp,
   ncomp,
   ncomp_check);

 if (ncomp_check!=nsolve)
  amrex::Error("nsolve invalid 2732");

 MultiFab& S_new = get_new_data(state_index,slab_step+1);
 for (int icomp=0;icomp<scomp.size();icomp++) {
  S_new.setVal(0.0,scomp[icomp],ncomp[icomp],0);
 }

} // end subroutine zero_independent_variable

void
NavierStokes::init_regrid_history() {

    is_first_step_after_regrid = 0;
    old_intersect_new          = grids;

} // end subroutine init_regrid_history

void
NavierStokes::restart (AmrCore&      papa,
                       std::istream& is,
		       int old_finest_level,
		       int new_finest_level) {

    AmrLevel::restart(papa,is,old_finest_level,new_finest_level);
    init_regrid_history();

} //end subroutine restart

//
// new_time=time  old_time=time-dt
//
void
NavierStokes::setTimeLevel (Real time,Real& dt)
{
 int nstate=state.size();
 if (nstate!=NUM_STATE_TYPE)
  amrex::Error("nstate invalid");
 for (int k=0;k<nstate;k++) 
  state[k].setTimeLevel(time,dt);
} //end subroutine setTimeLevel


void NavierStokes::debug_ngrow(int idxMF,int ngrow,
  const std::string& caller_string) {

 if ((idxMF<0)||(idxMF>=MAX_NUM_LOCAL_MF))
  amrex::Error("idxMF invalid");

 if (1==0) {
  std::cout << "full check of localMF integrity \n";

  for (int i=0;i<MAX_NUM_LOCAL_MF;i++) {
   if (localMF_grow[i]>=0) {
    MultiFab* mf_temp=localMF[i];
    if (! mf_temp->ok()) {
     amrex::Error("! mf_temp->ok()");
    }
   } else if (localMF_grow[i]==-1) {
    if (localMF[i]==0) {
     // do nothing
    } else {
     std::cout << "level = " << level << '\n';
     std::cout << "i = " << i << '\n';
     amrex::Error("localMF[i] invalid");
    }
   } else {
    amrex::Error("localMF_grow[i] invalid");
   }
  } // i=0 ... MAX_NUM_LOCAL_MF-1   
 } // full check of localMF integrity

 MultiFab* mf=localMF[idxMF];
 int mfgrow=localMF_grow[idxMF];

 if (mfgrow<ngrow) {
  std::cout << "caller_string= " << caller_string << '\n';
  std::cout << "idxMF= " << idxMF << '\n';
  std::cout << "mfgrow= " << mfgrow << " expected grow= " <<
   ngrow << '\n';
 }

 if (! mf->ok()) {
  std::cout << "caller_string= " << caller_string << '\n';
  std::cout << "idxMF= " << idxMF << '\n';
  amrex::Error("mf not ok");
 } else if ((mf->nGrow()<ngrow)||(mfgrow<ngrow)) {
  std::cout << "caller_string= " << caller_string << '\n';
  std::cout << "idxMF= " << idxMF << '\n';
  std::cout << "mf->ngrow= " << mf->nGrow() << " expected grow= " <<
   ngrow << '\n';
  std::cout << "mfgrow= " << mfgrow << " expected grow= " <<
   ngrow << '\n';
  amrex::Error("grow invalid in debug_ngrow");
 }

} // end subroutine debug_ngrow



void NavierStokes::debug_ixType(int idxMF,int grid_type,
  const std::string& caller_string) {

 if ((idxMF<0)||(idxMF>=MAX_NUM_LOCAL_MF))
  amrex::Error("idxMF invalid");

 if (1==0) {
  std::cout << "full check of localMF integrity \n";

  for (int i=0;i<MAX_NUM_LOCAL_MF;i++) {
   if (localMF_grow[i]>=0) {
    MultiFab* mf_temp=localMF[i];
    if (! mf_temp->ok()) {
     std::cout << "caller_string " << caller_string <<'\n';
     amrex::Error("! mf_temp->ok()");
    }
   } else if (localMF_grow[i]==-1) {
    if (localMF[i]==0) {
     // do nothing
    } else {
     std::cout << "caller_string " << caller_string <<'\n';
     std::cout << "level = " << level << '\n';
     std::cout << "i = " << i << '\n';
     amrex::Error("localMF[i] invalid");
    }
   } else {
    std::cout << "caller_string " << caller_string <<'\n';
    amrex::Error("localMF_grow[i] invalid");
   }
  } // i=0 ... MAX_NUM_LOCAL_MF-1   
 } // full check of localMF integrity

 MultiFab* mf=localMF[idxMF];

 if (! mf->ok()) {
  std::cout << "caller_string= " << caller_string << '\n';
  std::cout << "idxMF= " << idxMF << '\n';
  amrex::Error("mf not ok");
 } else if (mf->ok()) {
  debug_ixType_raw(mf,grid_type,caller_string);
 } else {
  std::cout << "caller_string= " << caller_string << '\n';
  amrex::Error("mf->ok corrupt");
 }

} // end subroutine debug_ixType


void NavierStokes::debug_ixType_raw(MultiFab* mf,int grid_type,
 const std::string& caller_string) {

 if (! mf->ok()) {
  std::cout << "caller_string= " << caller_string << '\n';
  amrex::Error("mf not ok");
 } else if (mf->ok()) {

  if (1==0) {
   std::fflush(NULL);
   int proc=ParallelDescriptor::MyProc();
   std::cout << "caller_string= " << caller_string << '\n';
   std::cout << "in debug_ixType_raw grid_type= " << grid_type << '\n';
   std::cout << "in debug_ixType_raw proc= " << proc << '\n';
   std::fflush(NULL);
   ParallelDescriptor::Barrier();
  }

  IndexType compare_typ;
  if (grid_type==-1) {
   compare_typ=IndexType::TheCellType();
  } else if (grid_type==0) {
   compare_typ=TheUMACType;
  } else if (grid_type==1) {
   compare_typ=TheVMACType;
  } else if ((grid_type==2)&&(AMREX_SPACEDIM==3)) {
   compare_typ=TheWMACType;
  } else if (grid_type==3) {
   compare_typ=TheYUMACType;
  } else if ((grid_type==4)&&(AMREX_SPACEDIM==3)) {
   compare_typ=TheZUMACType;
  } else if ((grid_type==5)&&(AMREX_SPACEDIM==3)) {
   compare_typ=TheZVMACType;
  } else
   amrex::Error("grid_type invalid");

  if (mf->boxArray().ixType()==compare_typ) {

   if (1==0) {
    std::fflush(NULL);
    int proc=ParallelDescriptor::MyProc();
    std::cout << "in debug_ixType_raw2 grid_type= " << grid_type << '\n';
    std::cout << "in debug_ixType_raw2 proc= " << proc << '\n';
    std::cout << "caller_string= " << caller_string << '\n';
    std::fflush(NULL);
    ParallelDescriptor::Barrier();
   }

  } else {
   std::cout << "caller_string= " << caller_string << '\n';
   amrex::Error("mf->boxArray().ixType()!=compare_typ");
  }

 } else
  amrex::Error("mf->ok corrupt");

} // end subroutine debug_ixType_raw



void NavierStokes::debug_boxArray(MultiFab* mf,int grid_type,
 const std::string& caller_string) {

 if (! mf->ok()) {
  std::cout << "caller_string= " << caller_string << '\n';
  amrex::Error("mf not ok");
 } else if (mf->ok()) {
  IndexType compare_typ;
  if (grid_type==-1) {
   compare_typ=IndexType::TheCellType();
  } else if (grid_type==0) {
   compare_typ=TheUMACType;
  } else if (grid_type==1) {
   compare_typ=TheVMACType;
  } else if ((grid_type==2)&&(AMREX_SPACEDIM==3)) {
   compare_typ=TheWMACType;
  } else if (grid_type==3) {
   compare_typ=TheYUMACType;
  } else if ((grid_type==4)&&(AMREX_SPACEDIM==3)) {
   compare_typ=TheZUMACType;
  } else if ((grid_type==5)&&(AMREX_SPACEDIM==3)) {
   compare_typ=TheZVMACType;
  } else
   amrex::Error("grid_type invalid");

  if (mf->boxArray().ixType()==compare_typ) {
   // do nothing
  } else
   amrex::Error("mf->boxArray().ixType()!=compare_typ");

   // this is the copy constructor:
   // BoxArray (const BoxArray& rhs);
  BoxArray grids_convert(grids);
  grids_convert.convert(compare_typ);
  if (mf->boxArray()==grids_convert) {
   // do nothing
  } else
   amrex::Error("mf->boxArray()!=grids_convert");

 } else
  amrex::Error("mf->ok corrupt");

} // end subroutine debug_boxArray



int NavierStokes::some_materials_compressible() {

 int comp_flag=0;
 for (int im=0;im<num_materials;im++) {
  int imat_type=material_type[im];
  if (imat_type==999) {

   //do nothing

  } else if (imat_type==0) {

   //do nothing

  } else if ((imat_type>=1)&&(imat_type<999)) {
   comp_flag=1;
  } else
   amrex::Error("material type invalid");
 }
 return comp_flag;
} //end subroutine some_materials_compressible()

// number of levels including the current that are "valid"
int NavierStokes::NSnumLevels() {

 int numLevelsMAX = 1024;

 int lv = numLevelsMAX;
    //
    // The routine `falls through' since coarsening and refining
    // a unit box does not yield the initial box.
    //
 const BoxArray& bs = grids;
 int ng=grids.size();

 for (int i = 0; i < ng; ++i) {
  int llv = 0;
  Box tmp = bs[i];
  for (;;) {
   Box ctmp  = tmp;   ctmp.coarsen(2);
   Box rctmp = ctmp; rctmp.refine(2);
   if (tmp != rctmp || ctmp.numPts() == 1)
    break;
   llv++;
   tmp = ctmp;
  }
  if (lv >= llv)
   lv = llv;
 }

 return lv+1; // Including coarsest.

} // end function int NavierStokes::NSnumLevels()

int NavierStokes::read_from_CAD() {

 int local_read_from_CAD=0;
 for (int im=0;im<num_materials;im++) {
  if (fort_read_from_CAD(&FSI_flag[im])==1) {
   local_read_from_CAD=1;
  } else if (fort_read_from_CAD(&FSI_flag[im])==0) {
   // do nothing
  } else
   amrex::Error("fort_read_from_CAD invalid");
 } // im=0..num_materials-1

 return local_read_from_CAD;

}  // end subroutine read_from_CAD


int NavierStokes::is_ice_matC(int im) {

 if ((im<0)|(im>=num_materials))
  amrex::Error("im invalid50b (is_ice_matC)");

 int imp1=im+1;
 int local_is_ice=fort_is_ice_base(&FSI_flag[im],&imp1);

 return local_is_ice;

}  // end subroutine is_ice_matC()

int NavierStokes::FSI_material_exists_CTML() {

 int local_flag=0;

 for (int im=0;im<num_materials;im++) {
  if (FSI_flag[im]==FSI_SHOELE_CTML) { //Wardlaw
   local_flag=1;
  } else if (FSI_flag_valid(im)==1) {
   // do nothing
  } else {
   amrex::Error("FSI_flag_valid(im)!=1");
  }
 }

 return local_flag;

}  // end subroutine FSI_material_exists_CTML()

int NavierStokes::FSI_flag_valid(int im) {

 if ((im<0)|(im>=num_materials))
  amrex::Error("im invalid50c (FSI_flag_valid)");

 int imp1=im+1;
 int local_is_valid=fort_FSI_flag_valid_base(&FSI_flag[im],&imp1);

 return local_is_valid;

}  //end function FSI_flag_valid

int NavierStokes::is_FSI_rigid_matC(int im) {

 if ((im<0)|(im>=num_materials))
  amrex::Error("im invalid50d (is_FSI_rigid_matC)");

 int imp1=im+1;
 int local_is_FSI_rigid=fort_is_FSI_rigid_base(&FSI_flag[im],&imp1);

 return local_is_FSI_rigid;

}  // end function is_FSI_rigid_matC()

int NavierStokes::is_singular_coeff(int im) {

 int local_is_singular_coeff=0;
 if ((im>=0)&&(im<num_materials)) {
  if (FSI_flag[im]==FSI_RIGID_NOTPRESCRIBED) {  
   local_is_singular_coeff=1; //extend pressure into this region
  } else if (FSI_flag[im]==FSI_PRESCRIBED_PROBF90) { 
   local_is_singular_coeff=1; //extend pressure
  } else if (FSI_flag[im]==FSI_PRESCRIBED_NODES) { 
   local_is_singular_coeff=1; //extend pressure
  } else if ((FSI_flag[im]==FSI_FLUID)||
             (FSI_flag[im]==FSI_FLUID_NODES_INIT)) { 
   local_is_singular_coeff=0;
  } else if ((FSI_flag[im]==FSI_ICE_PROBF90)||
  	     (FSI_flag[im]==FSI_ICE_STATIC)||
	     (FSI_flag[im]==FSI_ICE_NODES_INIT)) { 
   local_is_singular_coeff=1; //extend pressure (after SOLVETYPE_PRES)
  } else if (FSI_flag[im]==FSI_SHOELE_CTML) {  //Wardlaw
   local_is_singular_coeff=1;
  } else
   amrex::Error("FSI_flag invalid");
 } else
  amrex::Error("im invalid");

 return local_is_singular_coeff;

}  // end subroutine is_singular_coeff()


int NavierStokes::CTML_FSI_flagC() {

 int local_CTML_FSI_flag=0;

 for (int im=0;im<num_materials;im++) {
  int imp1=im+1;
  if (fort_CTML_FSI_mat_base(&FSI_flag[im],&imp1)==1) {
#ifdef MVAHABFSI
   local_CTML_FSI_flag=1;
#else
   amrex::Error("CTML(C): define MEHDI_VAHAB_FSI in GNUmakefile");
#endif
  } else if (fort_CTML_FSI_mat_base(&FSI_flag[im],&imp1)==0) {
   // do nothing
  } else
   amrex::Error("fort_CTML_FSI_mat invalid");
 } // im=0..num_materials-1

 return local_CTML_FSI_flag;

}  // end function CTML_FSI_flagC()


int NavierStokes::CTML_FSI_matC(int im) {

 if ((im<0)|(im>=num_materials))
  amrex::Error("im invalid50e (CTML_FSI_matC)");

 int imp1=im+1;
 int local_flag=fort_CTML_FSI_mat_base(&FSI_flag[im],&imp1);

 return local_flag;

}  // end subroutine CTML_FSI_matC(int im)



// passes tile information to sci_clsvof.F90 so that Lagrangian
// elements can be distributed amongst the tiles.
// called from FSI_make_distance and ns_header_msg_level 
// cur_time is passed just in case the actual node position depends on time.
// i.e. for finding target of characteristic given the foot.
void NavierStokes::create_fortran_grid_struct(Real cur_time,Real dt) {

 if (read_from_CAD()!=1)
  amrex::Error("read_from_CAD()!=1");

 const int max_level = parent->maxLevel();
 int finest_level=parent->finestLevel();
 if ((level<0)||(level>max_level))
  amrex::Error("level invalid in create_fortran_grid_struct");

 bool use_tiling=ns_tiling;

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

 const Real* dx = geom.CellSize();

 Real dx_max_level[AMREX_SPACEDIM];
 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  dx_max_level[dir]=dx[dir];
 for (int ilev=level+1;ilev<=max_level;ilev++) 
  for (int dir=0;dir<AMREX_SPACEDIM;dir++)
   dx_max_level[dir]/=2.0;

  // this will store information about the grids stored
  // on this processor.
 Vector<int> tilelo_array;
 Vector<int> tilehi_array;
 Vector<int> gridno_array;
 Vector<Real> xlo_array;
 Vector< int > num_tiles_on_thread_proc;
 Vector< int > num_tiles_on_thread_proc_check;
 int max_num_tiles_on_thread_proc=0;
 int tile_dim=0;

 num_tiles_on_thread_proc.resize(thread_class::nthreads);
 num_tiles_on_thread_proc_check.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  num_tiles_on_thread_proc[tid]=0;
  num_tiles_on_thread_proc_check[tid]=0;
 }

 int num_grids_on_level=grids.size();
 int num_grids_on_level_check=0;
 int num_grids_on_level_proc=0;
 
 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(S_new.boxArray().d_numPts());

 // if this loop is not inside of a "pragma omp parallel" command, then
 // omp_get_num_threads() returns 1 so that all tiles live on thread=1.
 // i.e.
 // The omp_get_num_threads routine returns the number of threads in the team
 // that is executing the parallel region to which the routine region binds.
 // If called from the sequential part of a program, omp_get_num_threads()
 // returns 1.
 for (MFIter mfi(S_new,false); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  const Box& tilegrid = mfi.tilebox();

  int tid_current=ns_thread();
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  const int gridno = mfi.index();
  if ((gridno>=0)&&(gridno<num_grids_on_level)) {
   num_grids_on_level_proc++;
   num_grids_on_level_check++;
  } else
   amrex::Error("gridno invalid");
 } // mfi
 ParallelDescriptor::ReduceIntSum(num_grids_on_level_check);
 ns_reconcile_d_num(LOOP_NUM_GRIDS,"create_fortran_grid_struct");

 if ((num_grids_on_level_proc<0)||
     (num_grids_on_level_proc>num_grids_on_level)) {
  std::cout << "num_grids_on_level_proc= " << num_grids_on_level_proc << '\n';
  std::cout << "num_grids_on_level= " << num_grids_on_level << '\n';
  std::cout << "num_grids_on_level_check= " << 
    num_grids_on_level_check << '\n';
  std::cout << "level= " << level << '\n';
  amrex::Error("num_grids_on_level_proc invalid");
 }

 if (num_grids_on_level_check!=num_grids_on_level) {
  std::cout << "num_grids_on_level_proc= " << num_grids_on_level_proc << '\n';
  std::cout << "num_grids_on_level= " << num_grids_on_level << '\n';
  std::cout << "num_grids_on_level_check= " << 
    num_grids_on_level_check << '\n';
  std::cout << "level= " << level << '\n';
  amrex::Error("num_grids_on_level_check invalid");
 }

 for (int grid_sweep=0;grid_sweep<2;grid_sweep++) {

  if (grid_sweep==0) {
   // do nothing
  } else if (grid_sweep==1) {

   for (int tid=0;tid<thread_class::nthreads;tid++) {

    if (num_tiles_on_thread_proc[tid]!=
	num_tiles_on_thread_proc_check[tid])
     amrex::Error("num_tiles_on_thread_proc[tid] failed check");

    if (num_grids_on_level_proc==0) {
     if (num_tiles_on_thread_proc[tid]!=0)
      amrex::Error("num_tiles_on_thread_proc[tid] invalid");
     if (max_num_tiles_on_thread_proc!=0)
      amrex::Error("max_num_tiles_on_thread_proc invalid");
    } else if (num_grids_on_level_proc>0) {

     if (num_tiles_on_thread_proc[tid]>max_num_tiles_on_thread_proc)
      max_num_tiles_on_thread_proc=num_tiles_on_thread_proc[tid];

    } else
     amrex::Error("num_grids_on_level_proc invalid");

   } // tid=0..thread_class::nthreads-1

   tile_dim=thread_class::nthreads*max_num_tiles_on_thread_proc;

   if (num_grids_on_level_proc==0) {
    if (tile_dim!=0)
     amrex::Error("tile_dim invalid");
    int tile_dim_virtual=1;
    tilelo_array.resize(tile_dim_virtual*AMREX_SPACEDIM);   
    tilehi_array.resize(tile_dim_virtual*AMREX_SPACEDIM);   
    gridno_array.resize(tile_dim_virtual);   // gridno for a given tile
    xlo_array.resize(tile_dim_virtual*AMREX_SPACEDIM);   
    for (int tid=0;tid<thread_class::nthreads;tid++) {
     if (num_tiles_on_thread_proc[tid]!=0)
      amrex::Error("num_tiles_on_thread_proc[tid] invalid");
    }
   } else if (num_grids_on_level_proc>0) {

    if ((tile_dim>=1)&&(tile_dim>=num_grids_on_level_proc)) {
     tilelo_array.resize(tile_dim*AMREX_SPACEDIM);   
     tilehi_array.resize(tile_dim*AMREX_SPACEDIM);   
     gridno_array.resize(tile_dim);   // gridno for a given tile
     xlo_array.resize(tile_dim*AMREX_SPACEDIM);   
    } else
     amrex::Error("tile_dim invalid");
   } else
    amrex::Error("num_grids_on_level_proc invalid");

  } else
   amrex::Error("grid_sweep invalid"); 

  if ((grid_sweep==0)||(grid_sweep==1)) {
   for (int tid=0;tid<thread_class::nthreads;tid++) {
    num_tiles_on_thread_proc[tid]=0;
   }
  } else
   amrex::Error("grid_sweep invalid");

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(S_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(S_new,use_tiling); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const Box& tilegrid = mfi.tilebox();

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();
 
   if (grid_sweep==0) {
    // do nothing
   } else if (grid_sweep==1) {
    const int gridno = mfi.index();
    const int* tilelo=tilegrid.loVect();
    const int* tilehi=tilegrid.hiVect();
    const Real* xlo = grid_loc[gridno].lo();

    if (num_grids_on_level_proc==0) {
     // do nothing
    } else if (num_grids_on_level_proc>0) {

     if (max_num_tiles_on_thread_proc>0) {
      int ibase=max_num_tiles_on_thread_proc*AMREX_SPACEDIM*tid_current+
       AMREX_SPACEDIM*num_tiles_on_thread_proc[tid_current];
      for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
       tilelo_array[ibase+dir]=tilelo[dir]; 
       tilehi_array[ibase+dir]=tilehi[dir]; 
       xlo_array[ibase+dir]=xlo[dir]; 
      } // dir=0..sdim-1
      ibase=max_num_tiles_on_thread_proc*tid_current+
       num_tiles_on_thread_proc[tid_current];

      if ((gridno>=0)&&(gridno<num_grids_on_level)) {
       gridno_array[ibase]=gridno;
      } else
       amrex::Error("gridno invalid");
     } else
      amrex::Error("max_num_tiles_on_thread_proc invalid");
    } else
     amrex::Error("num_grids_on_level_proc invalid");
   } else
    amrex::Error("grid_sweep invalid");

   num_tiles_on_thread_proc[tid_current]++;
   if (grid_sweep==0) {
    // check nothing
   } else if (grid_sweep==1) {
    if (num_tiles_on_thread_proc[tid_current]>max_num_tiles_on_thread_proc)
     amrex::Error("num_tiles_on_thread_proc[tid_current] invalid");
   } else
    amrex::Error("grid_sweep invalid");

  } // mfi
}//omp
  ns_reconcile_d_num(LOOP_GRIDNO_ARRAY,"create_fortran_grid_struct");

  for (int tid=0;tid<thread_class::nthreads;tid++) {
   if (grid_sweep==0) {
    num_tiles_on_thread_proc_check[tid]=num_tiles_on_thread_proc[tid];
   } else if (grid_sweep==1) {
    if (num_tiles_on_thread_proc_check[tid]!=
	num_tiles_on_thread_proc[tid])
     amrex::Error("num_tiles_on_thread_proc[tid] failed check");
   } else
    amrex::Error("grid_sweep invalid");
  }

 } // grid_sweep=0..1

 if (num_grids_on_level_proc==0) {
  if (max_num_tiles_on_thread_proc!=0)
   amrex::Error("max_num_tiles_on_thread_proc invalid");
 } else if (num_grids_on_level_proc>0) {
  if (max_num_tiles_on_thread_proc>0) {
   // do nothing
  } else
   amrex::Error("max_num_tiles_on_thread_proc invalid");
 } else
  amrex::Error("num_grids_on_level_proc invalid"); 

 int nparts=im_solid_map.size();
 if ((nparts<1)||(nparts>num_materials))
  amrex::Error("nparts invalid");

  // declared in: SOLIDFLUID.F90
 fort_fillcontainer(
  &level,
  &finest_level,
  &max_level,
  &cur_time,
  &dt,
  tilelo_array.dataPtr(),
  tilehi_array.dataPtr(),
  xlo_array.dataPtr(),
  dx,
  dx_max_level,
  &num_grids_on_level,
  &num_grids_on_level_proc,
  gridno_array.dataPtr(),
  num_tiles_on_thread_proc.dataPtr(),
  &thread_class::nthreads,
  &max_num_tiles_on_thread_proc,
  &tile_dim, 
  &nparts,
  im_solid_map.dataPtr());

} // end subroutine create_fortran_grid_struct

void NavierStokes::print_project_option(int project_option) {

 if (project_option==SOLVETYPE_PRES) {
  std::cout << "project_option= " << project_option << " (SOLVETYPE_PRES) \n";
 } else if (project_option==SOLVETYPE_PRESGRAVITY) {
  std::cout << "project_option= " << project_option << 
    " (SOLVETYPE_PRESGRAVITY) \n";
 } else if (project_option==SOLVETYPE_INITPROJ) {
  std::cout << "project_option= " << project_option << 
    " (SOLVETYPE_INITPROJ) \n";
 } else if (project_option==SOLVETYPE_HEAT) {
  std::cout << "project_option= " << project_option << 
    " (SOLVETYPE_HEAT) \n";
 } else if (project_option==SOLVETYPE_VISC) {
  std::cout << "project_option= " << project_option << 
    " (SOLVETYPE_VISC) \n";
 } else if (project_option==SOLVETYPE_PRESEXTRAP) {
  std::cout << "project_option= " << project_option << 
    " (SOLVETYPE_PRESEXTRAP) \n";
 } else if ((project_option>=SOLVETYPE_SPEC)&&
            (project_option<SOLVETYPE_SPEC+num_species_var)) {
  std::cout << "project_option= " << project_option << 
    " (SOLVETYPE_SPEC) \n";
 } else {
  std::cout << "project_option= " << project_option << 
    " (INVALID) \n";
  amrex::Error("project_option invalid print_project_option");
 }

} // end subroutine print_project_option

void NavierStokes::init_FSI_GHOST_MAC_MF_ALL_predict() {

 int finest_level=parent->finestLevel();
 if (level!=0)
  amrex::Error("level invalid init_FSI_GHOST_MAC_MF_ALL_predict");

 setup_integrated_quantities();

 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.init_FSI_GHOST_MAC_MF_predict();
 } // ilev=level...finest_level

} // end subroutine init_FSI_GHOST_MAC_MF_ALL_predict

void NavierStokes::init_FSI_GHOST_MAC_MF_predict() {

 int finest_level=parent->finestLevel();
 int nparts=im_solid_map.size();
 bool use_tiling=ns_tiling;

 if (ngrow_distance==4) {
  // do nothing
 } else
  amrex::Error("ngrow_distance invalid");

 int nparts_ghost=nparts;
 if (nparts==0) {
  nparts_ghost=1;
 } else if ((nparts>=1)&&(nparts<=num_materials)) {
  // do nothing
 } else {
  amrex::Error("nparts invalid");
 }

 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) { 

  if (localMF_grow[FSI_GHOST_MAC_MF+data_dir]==0) {
   delete_localMF(FSI_GHOST_MAC_MF+data_dir,1);
  } else if (localMF_grow[FSI_GHOST_MAC_MF+data_dir]==-1) {
   // do nothing
  } else
   amrex::Error("localMF_grow[FSI_GHOST_MAC_MF+data_dir] invalid");

  new_localMF(FSI_GHOST_MAC_MF+data_dir,
    nparts_ghost*AMREX_SPACEDIM,0,data_dir);
 }

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 int nstate=STATE_NCOMP;
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");

 MultiFab* solid_vel_mf;
 if (nparts==0) {
  if (nparts_ghost==1) {
   solid_vel_mf=getState(ngrow_distance,STATECOMP_VEL,
    STATE_NCOMP_VEL,
    cur_time_slab);
  } else
   amrex::Error("nparts_ghost invalid");
 } else if (nparts>0) {
  solid_vel_mf=getStateSolid(ngrow_distance,0,
    nparts*AMREX_SPACEDIM,cur_time_slab);
 } else
  amrex::Error("nparts invalid");

  // velocity and pressure
 MultiFab* fluid_vel_mf=getState(ngrow_distance,STATECOMP_VEL,
   STATE_NCOMP_VEL+STATE_NCOMP_PRES,
   cur_time_slab);

 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) { 

  const Real* dx = geom.CellSize();

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(fluid_vel_mf->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(*fluid_vel_mf,use_tiling); mfi.isValid(); ++mfi) {
    BL_ASSERT(grids[mfi.index()] == mfi.validbox());
    const int gridno = mfi.index();
    const Box& tilegrid = mfi.tilebox();
    const Box& fabgrid = grids[gridno];
    const int* tilelo=tilegrid.loVect();
    const int* tilehi=tilegrid.hiVect();
    const int* fablo=fabgrid.loVect();
    const int* fabhi=fabgrid.hiVect();
    int bfact=parent->Space_blockingFactor(level);

    const Real* xlo = grid_loc[gridno].lo();

    FArrayBox& fluidvelfab=(*fluid_vel_mf)[mfi]; //ngrow_distance=4
    FArrayBox& solidvelfab=(*solid_vel_mf)[mfi]; //ngrow_distance=4
    FArrayBox& ghostsolidvelfab=(*localMF[FSI_GHOST_MAC_MF+data_dir])[mfi]; 

    int tid_current=ns_thread();
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    fort_wallfunction_predict( 
     &data_dir,
     im_solid_map.dataPtr(),
     &level,
     &finest_level,
     &ngrow_distance,
     &nparts,
     &nparts_ghost,
     tilelo,tilehi,
     fablo,fabhi,&bfact,
     xlo,dx,
     &dt_slab,
     &cur_time_slab,
     fluidvelfab.dataPtr(),
     ARLIM(fluidvelfab.loVect()),ARLIM(fluidvelfab.hiVect()),
     solidvelfab.dataPtr(),
     ARLIM(solidvelfab.loVect()),ARLIM(solidvelfab.hiVect()),
     ghostsolidvelfab.dataPtr(),
     ARLIM(ghostsolidvelfab.loVect()),ARLIM(ghostsolidvelfab.hiVect()));
  } // mfi
} // omp
  ns_reconcile_d_num(LOOP_WALLFUNC_PRED,"init_FSI_GHOST_MAC_MF_predict"); 
   //thread_class::sync_tile_d_numPts(),
   //ParallelDescriptor::ReduceRealSum
   //thread_class::reconcile_d_numPts(caller_loop_id,caller_string)

 } // data_dir=0..sdim-1


 delete fluid_vel_mf;
 delete solid_vel_mf;

} // end subroutine init_FSI_GHOST_MAC_MF_predict

int NavierStokes::pattern_test(const std::string& source,
 const std::string& pattern) {

 std::size_t local_found=source.find(pattern);
 int return_val=((local_found<std::string::npos) ? 1 : 0);
 return return_val;
}

void NavierStokes::init_FSI_GHOST_MAC_MF_ALL(
  int renormalize_only,
  const std::string& caller_string) {

 int finest_level=parent->finestLevel();
 if (level!=0)
  amrex::Error("level invalid init_FSI_GHOST_MAC_MF_ALL");

 setup_integrated_quantities();
 int fast_mode=1;

 std::string local_caller_string="init_FSI_GHOST_MAC_MF_ALL";
 local_caller_string=caller_string+local_caller_string;

 volWgtSumALL(local_caller_string,fast_mode);

 getStateCONDUCTIVITY_ALL();

 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  int dealloc_history=0;
  ns_level.init_FSI_GHOST_MAC_MF(dealloc_history);
 } // ilev=level...finest_level

  // GNBC DEBUGGING
  //    --------------
  //    | Image cell |
  //    ------|-------   | = interface cell on MAC grid.
  //    | solid cell |   LS_solid(x_image)<0  LS_solid(x_solid)>0
  //    --------------
  //
  //    angle = angle measured at the solid normal probe in the fluid
  //    region   grad LS_solid dot grad LS_fluid = cos(theta) ?

 if (renormalize_only==0) {
  if (pattern_test(local_caller_string,"nonlinear_advection")==1) {
   if (pattern_test(local_caller_string,"prescribe_solid_geometryALL")==1) {

    for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {

     //mat<stuff>.tec has reconstructed interface.
     // (visit can open ascii surface mesh files)
     //
     //MAC grid law of the wall information.
     //WALLFUNCTION<stuff>.plt (visit can open binary tecplot files)
     if (1==0) {
      writeSanityCheckData(
       "WALLFUNCTION",
       "GNBC DEBUGGING velINT,imgVR,solVR,angle",
       local_caller_string,
       HISTORY_MAC_MF+data_dir, // tower_mf_id
       //velINT,image vel,velsol,image vel raster,velsol raster,angle
       localMF[HISTORY_MAC_MF+data_dir]->nComp(), 
       HISTORY_MAC_MF+data_dir,
       -1,  // State_Type==-1 
       data_dir,
       parent->levelSteps(0)); 
     }

     if (visual_WALLVEL_plot_int>0) {
      if (very_last_sweep==1) {
       int nsteps=parent->levelSteps(0)+1; // nsteps==0 very first step.
       int ratio=nsteps/visual_WALLVEL_plot_int;
       ratio=ratio*visual_WALLVEL_plot_int;
       if (ratio==nsteps) {

        int nparts=im_solid_map.size();
        int nparts_def=nparts;
        if (nparts==0) {
         nparts_def=1;
        } else if ((nparts>=1)&&(nparts<=num_materials)) {
         //do nothing
        } else
         amrex::Error("nparts invalid");

        if (localMF[FSI_GHOST_MAC_MF+data_dir]->nComp()!= 
            nparts_def*AMREX_SPACEDIM)
         amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nComp() bad");
        if (localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0)
         amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0");

        //MAC grid rasterized solid boundary condition.
        //WALLVEL<stuff>.plt (visit can open binary tecplot files)
        writeSanityCheckData(
         "WALLVEL",
         "init_FSI_GHOST_MAC_MF_ALL, FSI_GHOST_MAC_MF",//fictitious sol vel
         local_caller_string,
         FSI_GHOST_MAC_MF+data_dir, //tower_mf_id
         localMF[FSI_GHOST_MAC_MF+data_dir]->nComp(),
         FSI_GHOST_MAC_MF+data_dir,
         -1,  // State_Type==-1 
         data_dir,
         nsteps); 
       } else if (ratio<nsteps) {
        //do nothing
       } else
        amrex::Error("ratio or nsteps invalid");

      } else if (very_last_sweep==0) {
       // do nothing
      } else
       amrex::Error("very_last_sweep invalid");
     } else if (visual_WALLVEL_plot_int==0) {
      // do nothing
     } else
      amrex::Error("visual_WALLVEL_plot_int invalid");
    } //data_dir=0..sdim-1
   } // called from prescribe_solid_geometryALL?
  } // called from nonlinear_advection?
 } else if (renormalize_only==1) {
  //do nothing
 } else
  amrex::Error("renormalize_only invalid");
			 
 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
  delete_array(HISTORY_MAC_MF+data_dir);
 }

} // end subroutine init_FSI_GHOST_MAC_MF_ALL

void NavierStokes::init_FSI_GHOST_MAC_MF(int dealloc_history) {

 int finest_level=parent->finestLevel();
 int nparts=im_solid_map.size();
 bool use_tiling=ns_tiling;

 std::string local_caller_string="init_FSI_GHOST_MAC_MF";

 int nparts_ghost=nparts;
 if (nparts==0) {
  nparts_ghost=1;
 } else if ((nparts>=1)&&(nparts<=num_materials)) {
  // do nothing
 } else {
  amrex::Error("nparts invalid");
 }

 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) { 

  if (localMF_grow[FSI_GHOST_MAC_MF+data_dir]==0) {
   delete_localMF(FSI_GHOST_MAC_MF+data_dir,1);
  } else if (localMF_grow[FSI_GHOST_MAC_MF+data_dir]==-1) {
   // do nothing
  } else
   amrex::Error("localMF_grow[FSI_GHOST_MAC_MF+data_dir] invalid");

  new_localMF(FSI_GHOST_MAC_MF+data_dir,
    nparts_ghost*AMREX_SPACEDIM,0,data_dir);
 }

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 int nstate=STATE_NCOMP;
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");

 if (ngrow_distance==4) {
  // do nothing
 } else
  amrex::Error("ngrow_distance invalid");

  // usolid_law_of_the_wall,
  // uimage raster,usolid raster,angle_ACT_cell
 int nhistory_sub=3*AMREX_SPACEDIM+1;
 int nhistory=nparts_ghost*nhistory_sub;

 MultiFab* solid_vel_mf;
 if (nparts==0) {
  if (nparts_ghost==1) {
   solid_vel_mf=getState(ngrow_distance,0,AMREX_SPACEDIM,
    cur_time_slab);
  } else
   amrex::Error("nparts_ghost invalid");
 } else if (nparts>0) {
  solid_vel_mf=getStateSolid(ngrow_distance,0,
    nparts*AMREX_SPACEDIM,cur_time_slab);
 } else
  amrex::Error("nparts invalid");

  // velocity and pressure
 MultiFab* fluid_vel_mf=getState(ngrow_distance,STATECOMP_VEL,
   STATE_NCOMP_VEL+STATE_NCOMP_PRES,cur_time_slab);

  // temperature and density for all of the materials.
 int nden=num_materials*num_state_material;
 MultiFab* state_var_mf=getStateDen(ngrow_distance,cur_time_slab);
 if (state_var_mf->nComp()==nden) {
  // do nothing
 } else
  amrex::Error("state_var_mf->nComp()!=nden");

 if (ngrow_distance==4) {
  // do nothing
 } else
  amrex::Error("ngrow_distance invalid");

 getStateDist_localMF(LS_NRM_CP_MF,ngrow_distance,cur_time_slab,
  		      local_caller_string);
 if (localMF[LS_NRM_CP_MF]->nGrow()!=ngrow_distance)
  amrex::Error("localMF[LS_NRM_CP_MF]->nGrow()!=ngrow_distance");
 if (localMF[LS_NRM_CP_MF]->nComp()!=num_materials*(AMREX_SPACEDIM+1))
  amrex::Error("localMF[LS_NRM_CP_MF]->nComp()!=num_materials*(AMREX_SPACEDIM+1)");

 new_localMF(LS_NRM_FD_GNBC_MF,num_materials*AMREX_SPACEDIM,ngrow_distance,-1);
 build_NRM_FD_MF(LS_NRM_FD_GNBC_MF,LS_NRM_CP_MF);

 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) { 

  if (localMF_grow[HISTORY_MAC_MF+data_dir]==-1) {
   // do nothing
  } else
   amrex::Error("localMF_grow[HISTORY_MAC_MF+data_dir] invalid");

  new_localMF(HISTORY_MAC_MF+data_dir,nhistory,0,data_dir);

  for (int im=0;im<num_materials;im++) {

   if ((law_of_the_wall[im]==0)||   //just use the solid velocity
       (law_of_the_wall[im]==1)||   //turbulent wall flux
       (law_of_the_wall[im]==2)) {  //GNBC
    // do nothing
   } else
    amrex::Error("law_of_the_wall[im] invalid");

  }

  const Real* dx = geom.CellSize();

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(localMF[LS_NRM_CP_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(*localMF[LS_NRM_CP_MF],use_tiling); mfi.isValid(); ++mfi) {
    BL_ASSERT(grids[mfi.index()] == mfi.validbox());
    const int gridno = mfi.index();
    const Box& tilegrid = mfi.tilebox();
    const Box& fabgrid = grids[gridno];
    const int* tilelo=tilegrid.loVect();
    const int* tilehi=tilegrid.hiVect();
    const int* fablo=fabgrid.loVect();
    const int* fabhi=fabgrid.hiVect();
    int bfact=parent->Space_blockingFactor(level);

    const Real* xlo = grid_loc[gridno].lo();

    FArrayBox& lsCPfab=(*localMF[LS_NRM_CP_MF])[mfi];
    FArrayBox& lsFDfab=(*localMF[LS_NRM_FD_GNBC_MF])[mfi];

    FArrayBox& statefab=(*state_var_mf)[mfi]; //ngrow_distance=4
    FArrayBox& fluidvelfab=(*fluid_vel_mf)[mfi]; //ngrow_distance=4
    FArrayBox& solidvelfab=(*solid_vel_mf)[mfi]; //ngrow_distance=4
    FArrayBox& ghostsolidvelfab=(*localMF[FSI_GHOST_MAC_MF+data_dir])[mfi]; 

    FArrayBox& histfab=(*localMF[HISTORY_MAC_MF+data_dir])[mfi]; 

    int tid_current=ns_thread();
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    int nhistory_local=histfab.nComp();

     // CODY ESTEBE: LAW OF THE WALL
     // "HISTORY_MAC_MF" contains image velocity data and angle data.
     // fab = fortran array block
     // DATA: state data (velocity, level set function, temperature, density,
     // pressure) are data that represent the state of some fluid dynamics
     // system.
     // "localMF" data are temporary data used to assist in updating the
     // state data from time t=t^n to time t=t^{n+1}.
     // state data is not discarded at the end of a time step.
     // At the beginning of new time steps, existing state data is copied
     // to "state data old" and then by the end of a time step, 
     // "state data new" contains the updated state.
     // Generalized Navier Boundary Condition GNBC (call get_use_DCA)
     // 1. copy solid velocity into ghost velocity where phi_solid>0
     // 2. copy solid velocity into ghost velocity where phi_solid<-|cutoff|
     // 3. overwrite ghost velocity with law of the wall or GNBC velocity
     //    where phi_solid > -|cutoff|.  (and solid is rigid)
     //    (tangential velocity in phi_solid>0 regions is replaced by
     //     ghost tangential velocity, and unchanged where
     //     0>phi_solid>-|cutoff|)
     //    ghost normal velocity = solid normal velocity everywhere.
     // declared in: GODUNOV_3D.F90

    int local_sumdata_size=NS_sumdata.size();

    fort_wallfunction( 
     &data_dir,
     law_of_the_wall.dataPtr(),
     &local_sumdata_size,
     NS_sumdata.dataPtr(),
     &ncomp_sum_int_user1,
     &ncomp_sum_int_user2,
     wall_model_velocity.dataPtr(),
     im_solid_map.dataPtr(),
     &level,
     &finest_level,
     &ngrow_distance,
     &nparts,
     &nparts_ghost,
     &nden,
     tilelo,tilehi,
     fablo,fabhi,&bfact,
     xlo,dx,
     &dt_slab,//wall_function
     &cur_time_slab,
     lsCPfab.dataPtr(),ARLIM(lsCPfab.loVect()),ARLIM(lsCPfab.hiVect()),
     lsFDfab.dataPtr(),ARLIM(lsFDfab.loVect()),ARLIM(lsFDfab.hiVect()),
     statefab.dataPtr(),ARLIM(statefab.loVect()),ARLIM(statefab.hiVect()),
     fluidvelfab.dataPtr(),
     ARLIM(fluidvelfab.loVect()),ARLIM(fluidvelfab.hiVect()),
     solidvelfab.dataPtr(),
     ARLIM(solidvelfab.loVect()),ARLIM(solidvelfab.hiVect()),
     ghostsolidvelfab.dataPtr(),
     ARLIM(ghostsolidvelfab.loVect()),ARLIM(ghostsolidvelfab.hiVect()),
     histfab.dataPtr(),
     ARLIM(histfab.loVect()),ARLIM(histfab.hiVect()),
     &nhistory_local,
     &visc_coef);
  } // mfi
} // omp
  ns_reconcile_d_num(LOOP_WALLFUNCTION,"init_FSI_GHOST_MAC_MF"); 
   //thread_class::sync_tile_d_numPts(),
   //ParallelDescriptor::ReduceRealSum
   //thread_class::reconcile_d_numPts(caller_loop_id,caller_string)

  if (dealloc_history==0) {
   // do nothing
  } else if (dealloc_history==1) {
   delete_localMF(HISTORY_MAC_MF+data_dir,1);
  } else 
   amrex::Error("dealloc_history invalid");

 } // data_dir=0..sdim-1

 delete_localMF(LS_NRM_CP_MF,1);
 delete_localMF(LS_NRM_FD_GNBC_MF,1);

 delete state_var_mf;
 delete fluid_vel_mf;
 delete solid_vel_mf;

} // end subroutine init_FSI_GHOST_MAC_MF


// called from veldiffuseALL (NavierStokes3.cpp)
void NavierStokes::assimilate_state_data() {

 int finest_level=parent->finestLevel();
 int nparts=im_solid_map.size();
 bool use_tiling=ns_tiling;

 int nparts_ghost=nparts;
 if (nparts==0) {
  nparts_ghost=1;
 } else if ((nparts>=1)&&(nparts<=num_materials)) {
  // do nothing
 } else {
  amrex::Error("nparts invalid");
 }

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 MultiFab& LS_new=get_new_data(LS_Type,slab_step+1);

 if (LS_new.nComp()==num_materials*(1+AMREX_SPACEDIM)) {
  // do nothing
 } else
  amrex::Error("LS_new.nComp() invalid");

 int nstate=STATE_NCOMP;
 if (nstate!=S_new.nComp()) {
  std::cout << "nstate= " << nstate << '\n';
  amrex::Error("nstate invalid in cpp assimilate");
 }

 if (ngrow_distance==4) {
  // do nothing
 } else
  amrex::Error("ngrow_distance invalid");

 const Real* dx = geom.CellSize();

 MultiFab& Smac_new_x = get_new_data(Umac_Type,slab_step+1);
 MultiFab& Smac_new_y = get_new_data(Umac_Type+1,slab_step+1);
 MultiFab& Smac_new_z = get_new_data(Umac_Type+AMREX_SPACEDIM-1,slab_step+1);

 for (int isweep=0;isweep<2;isweep++) {
  init_boundary(); // init ghost cells on the given level.

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(S_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(S_new,use_tiling); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();
   int bfact=parent->Space_blockingFactor(level);

   const Real* xlo = grid_loc[gridno].lo();

   FArrayBox& solidvelx=(*localMF[FSI_GHOST_MAC_MF])[mfi]; 
   FArrayBox& solidvely=(*localMF[FSI_GHOST_MAC_MF+1])[mfi]; 
   FArrayBox& solidvelz=(*localMF[FSI_GHOST_MAC_MF+AMREX_SPACEDIM-1])[mfi]; 
   FArrayBox& snewfab=S_new[mfi]; 
   FArrayBox& lsnewfab=LS_new[mfi]; 
   FArrayBox& macnewx=Smac_new_x[mfi]; 
   FArrayBox& macnewy=Smac_new_y[mfi]; 
   FArrayBox& macnewz=Smac_new_z[mfi]; 

   int tid_current=ns_thread();
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // isweep==0 => update cell data
    // isweep==1 => update MAC data 
    // in: GODUNOV_3D.F90
   fort_assimilate_statedata( 
     &isweep,
     law_of_the_wall.dataPtr(), //currently unused in this routine.
     &wall_slip_weight,
     static_damping_coefficient.dataPtr(),
     im_solid_map.dataPtr(),
     &level,
     &finest_level,
     &nstate,
     &nparts,
     &nparts_ghost,
     tilelo,tilehi,
     fablo,fabhi,&bfact,
     xlo,dx,
     &dt_slab,
     &cur_time_slab,
     lsnewfab.dataPtr(),
     ARLIM(lsnewfab.loVect()),ARLIM(lsnewfab.hiVect()),
     snewfab.dataPtr(),
     ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()),
     macnewx.dataPtr(),
     ARLIM(macnewx.loVect()),ARLIM(macnewx.hiVect()),
     macnewy.dataPtr(),
     ARLIM(macnewy.loVect()),ARLIM(macnewy.hiVect()),
     macnewz.dataPtr(),
     ARLIM(macnewz.loVect()),ARLIM(macnewz.hiVect()),
     solidvelx.dataPtr(),
     ARLIM(solidvelx.loVect()),ARLIM(solidvelx.hiVect()),
     solidvely.dataPtr(),
     ARLIM(solidvely.loVect()),ARLIM(solidvely.hiVect()),
     solidvelz.dataPtr(),
     ARLIM(solidvelz.loVect()),ARLIM(solidvelz.hiVect()) );
  } // mfi
} // omp
  ns_reconcile_d_num(LOOP_ASSIMILATE,"assimilate_state_data"); 
        //thread_class::sync_tile_d_numPts(),
        //ParallelDescriptor::ReduceRealSum
        //thread_class::reconcile_d_numPts(caller_loop_id,caller_string)
 } // isweep=0..1

} // end subroutine assimilate_state_data()

// get rid of the ghost cells
void NavierStokes::resize_FSI_MF() {

 int nparts=im_solid_map.size();
 if (nparts==0) {
  // do nothing
 } else if ((nparts>=1)&&(nparts<=num_materials)) {
  int nFSI=nparts*NCOMP_FSI;
  if (localMF[FSI_MF]->nComp()!=nFSI)
   amrex::Error("localMF[FSI_MF]->nComp()!=nFSI");
  if (localMF[FSI_MF]->nGrow()==0) {
   // do nothing
  } else if (localMF[FSI_MF]->nGrow()>0) {
   MultiFab* save_FSI=
    new MultiFab(grids,dmap,nFSI,0,
     MFInfo().SetTag("save_FSI"),FArrayBoxFactory()); 
   MultiFab::Copy(*save_FSI,*localMF[FSI_MF],0,0,nFSI,0);
   delete_localMF(FSI_MF,1);
   new_localMF(FSI_MF,nFSI,0,-1);
   MultiFab::Copy(*localMF[FSI_MF],*save_FSI,0,0,nFSI,0);
   delete save_FSI;
  } else
   amrex::Error("localMF[FSI_MF]->nGrow() invalid");

 } else {
  amrex::Error("nparts invalid");
 }

} // end subroutine resize_FSI_MF

void NavierStokes::regenerate_from_eulerian(Real cur_time) {

 bool use_tiling=ns_tiling;

 const Real* dx = geom.CellSize();

 int nparts=im_solid_map.size();

 if ((nparts>=1)&&(nparts<=num_materials)) {

  MultiFab& Solid_new = get_new_data(Solid_State_Type,slab_step+1);
  if (Solid_new.nComp()!=nparts*AMREX_SPACEDIM)
   amrex::Error("Solid_new.nComp()!=nparts*AMREX_SPACEDIM");
  MultiFab& S_new=get_new_data(State_Type,slab_step+1);
  if (S_new.nComp()!=STATE_NCOMP)
   amrex::Error("S_new.nComp()!=STATE_NCOMP");
  MultiFab& LS_new=get_new_data(LS_Type,slab_step+1);
  if (LS_new.nComp()!=num_materials*(1+AMREX_SPACEDIM))
   amrex::Error("LS_new.nComp()!=num_materials*(1+AMREX_SPACEDIM)");

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(Solid_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(Solid_new,use_tiling); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();
   int bfact=parent->Space_blockingFactor(level);
   const Real* xlo = grid_loc[gridno].lo();
   const Real* xhi = grid_loc[gridno].hi();

   int tid_current=ns_thread();
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   FArrayBox& snewfab=S_new[mfi];
   FArrayBox& lsnewfab=LS_new[mfi];
   FArrayBox& solidfab=Solid_new[mfi];

    // updates Solid_new for 
    // FSI_flag(im)==FSI_PRESCRIBED_PROBF90 materials.
    // (prescribed solid from PROB.F90, not from CAD)
    // fort_initdatasolid is declared in PROB.F90
   fort_initdatasolid(
     &nparts,
     &ngrow_make_distance,
     im_solid_map.dataPtr(),
     &cur_time,
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,
     solidfab.dataPtr(),
     ARLIM(solidfab.loVect()),ARLIM(solidfab.hiVect()),
     lsnewfab.dataPtr(),
     ARLIM(lsnewfab.loVect()),ARLIM(lsnewfab.hiVect()),
     snewfab.dataPtr(),
     ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()),
     dx,xlo,xhi);  

  } // mfi
} // omp
  ns_reconcile_d_num(LOOP_INITDATA_SOLID,"regenerate_from_eulerian");

 } else if (nparts==0) {
  // do nothing
 } else
  amrex::Error("nparts invalid");

} // end subroutine regenerate_from_eulerian(cur_time)

// create a distance function (velocity and temperature) on this level.
// calls fill coarse patch if level>0
// called from: NavierStokes::initData ()
//              NavierStokes::nonlinear_advection()
// cur_time=prev_time+dt
void NavierStokes::FSI_make_distance(Real cur_time,Real dt) {

 std::string local_caller_string="FSI_make_distance";

 if (read_from_CAD()==1) {
  // do nothing
 } else
  amrex::Error("expecting read_from_CAD()==1");

 int nparts=im_solid_map.size();

 if (nparts==0) {
  // do nothing
 } else if ((nparts>=1)&&(nparts<=num_materials)) {

  int nFSI=nparts*NCOMP_FSI;

  delete_localMF_if_exist(FSI_MF,1);

  new_localMF(FSI_MF,nFSI,ngrow_make_distance,-1);

  for (int partid=0;partid<nparts;partid++) {
   int ibase=partid*NCOMP_FSI;
   setVal_localMF(FSI_MF,0.0,ibase+FSI_VELOCITY,3,ngrow_make_distance); 
   setVal_localMF(FSI_MF,-99999.0,ibase+FSI_LEVELSET,1,ngrow_make_distance);
     // =0 no hits
     // =1 positive
     // =2 negative
     // =3 inconclusive
   setVal_localMF(FSI_MF,0.0,ibase+FSI_SIGN_CONFLICT,1,ngrow_make_distance);
   setVal_localMF(FSI_MF,0.0,ibase+FSI_TEMPERATURE,1,ngrow_make_distance);
   setVal_localMF(FSI_MF,FSI_NOTHING_VALID,
		  ibase+FSI_EXTRAP_FLAG,1,ngrow_make_distance);
  } // partid=0..nparts-1


   // in: NavierStokes::FSI_make_distance
   // 1. create lagrangian container data structure within the 
   //    fortran part that recognizes tiles. 
   //    (fort_fillcontainer in SOLIDFLUID.F90)
   // 2. fill the containers with the Lagrangian information.
   //    (CLSVOF_FILLCONTAINER called from fort_fillcontainer)
   //    i.e. associate to each tile a set of Lagrangian nodes and elements
   //    that are located in or very near the tile.
  create_fortran_grid_struct(cur_time,dt);

  int iter=0;
  resize_mask_nbr(ngrow_make_distance);
   // in: FSI_make_distance
   // 1.FillCoarsePatch
   // 2.traverse lagrangian elements belonging to each tile and update
   //   cells within "bounding box" of the element.
  ns_header_msg_level(
     OP_FSI_MAKE_DISTANCE,
     SUB_OP_FSI_DEFAULT,
     cur_time,
     dt,
     iter,
     local_caller_string);
  
  do {
 
   ns_header_msg_level(
     OP_FSI_MAKE_SIGN,
     SUB_OP_FSI_DEFAULT,
     cur_time,
     dt,
     iter,
     local_caller_string);

   iter++;
  
  } while (FSI_touch_flag[0]==1);

  build_moment_from_FSILS(cur_time);

   // Solid velocity
  MultiFab& Solid_new = get_new_data(Solid_State_Type,slab_step+1);
  if (Solid_new.nComp()!=nparts*AMREX_SPACEDIM)
   amrex::Error("Solid_new.nComp()!=nparts*AMREX_SPACEDIM");

  for (int partid=0;partid<nparts;partid++) {

   int im_part=im_solid_map[partid];
   if ((im_part<0)||(im_part>=num_materials))
    amrex::Error("im_part invalid");

   int ok_to_modify_EUL=1;
   if ((FSI_flag[im_part]==FSI_ICE_NODES_INIT)||  
       (FSI_flag[im_part]==FSI_FLUID_NODES_INIT)) { 
    if (cur_time==0.0) {
     ok_to_modify_EUL=1;
    } else if (cur_time>0.0) {
     ok_to_modify_EUL=0;
    } else
     amrex::Error("cur_time invalid");
   } else if ((FSI_flag[im_part]==FSI_PRESCRIBED_NODES)|| 
              (FSI_flag[im_part]==FSI_SHOELE_CTML)) {//Wardlaw
    ok_to_modify_EUL=1;
    // do nothing
   } else
    amrex::Error("FSI_flag invalid");

   if (ok_to_modify_EUL==1) {

    if (fort_read_from_CAD(&FSI_flag[im_part])==1) { 
     int ibase=partid*NCOMP_FSI;
     MultiFab::Copy(Solid_new,*localMF[FSI_MF],ibase+FSI_VELOCITY,
      partid*AMREX_SPACEDIM,
      AMREX_SPACEDIM,0);
    } else if (FSI_flag[im_part]==FSI_PRESCRIBED_PROBF90) { 
     // do nothing
    } else
     amrex::Error("FSI_flag invalid");

   } else if (ok_to_modify_EUL==0) {
    // do nothing
   } else
    amrex::Error("ok_to_modify_EUL invalid");

  } // partid=0..nparts-1

 } else {
  amrex::Error("nparts invalid");
 }

}  // end subroutine FSI_make_distance

// called from: Transfer_FSI_To_STATE()
// Transfer_FSI_To_STATE() called from: ns_header_msg_level,initData ()
void NavierStokes::copy_velocity_on_sign(int partid) {

 std::string local_caller_string="copy_velocity_on_sign";

 int nparts=im_solid_map.size();
 if ((nparts<1)||(nparts>num_materials))
  amrex::Error("nparts invalid");
 if ((partid<0)||(partid>=nparts))
  amrex::Error("partid invalid");
 debug_ngrow(FSI_MF,ngrow_make_distance,local_caller_string);

 int im_part=im_solid_map[partid];
 if ((im_part<0)||(im_part>=num_materials))
  amrex::Error("im_part invalid");

 if (fort_read_from_CAD(&FSI_flag[im_part])==1) { 

  if ((FSI_flag[im_part]==FSI_PRESCRIBED_NODES)|| 
      (FSI_flag[im_part]==FSI_ICE_NODES_INIT)|| 
      (FSI_flag[im_part]==FSI_SHOELE_CTML)||
      (FSI_flag[im_part]==FSI_FLUID_NODES_INIT)) {

   MultiFab& S_new=get_new_data(State_Type,slab_step+1);
   int nstate=STATE_NCOMP;
   if (nstate!=S_new.nComp())
    amrex::Error("nstate invalid");

   int nFSI=nparts*NCOMP_FSI;
   if (localMF[FSI_MF]->nComp()!=nFSI)
    amrex::Error("localMF[FSI_MF]->nComp()!=nFSI");
  
   bool use_tiling=ns_tiling;

   if (ns_is_rigid(im_part)==1) {

    const Real* dx = geom.CellSize();

    if (thread_class::nthreads<1)
     amrex::Error("thread_class::nthreads invalid");
    thread_class::init_d_numPts(S_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
    for (MFIter mfi(S_new,use_tiling); mfi.isValid(); ++mfi) {
     BL_ASSERT(grids[mfi.index()] == mfi.validbox());
     const int gridno = mfi.index();
     const Box& tilegrid = mfi.tilebox();
     const Box& fabgrid = grids[gridno];
     const int* tilelo=tilegrid.loVect();
     const int* tilehi=tilegrid.hiVect();
     const int* fablo=fabgrid.loVect();
     const int* fabhi=fabgrid.hiVect();
     int bfact=parent->Space_blockingFactor(level);

     const Real* xlo = grid_loc[gridno].lo();

     FArrayBox& snewfab=S_new[mfi];
     FArrayBox& fsifab=(*localMF[FSI_MF])[mfi];

     int tid_current=ns_thread();
     thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     // declared in: GODUNOV_3D.F90
     // copies velocity from fsifab to snewfab if fsifab_LS>0
     fort_copy_vel_on_sign(
      &im_part, 
      &nparts,
      &partid, 
      &ngrow_make_distance, 
      &nFSI, 
      xlo,dx,
      snewfab.dataPtr(),ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()),
      fsifab.dataPtr(),ARLIM(fsifab.loVect()),ARLIM(fsifab.hiVect()),
      tilelo,tilehi,
      fablo,fabhi,&bfact,
      &nstate);
    }  // mfi  
}//omp
    ns_reconcile_d_num(LOOP_COPY_VEL_ON_SIGN,"copy_velocity_on_sign");

   } else if (ns_is_rigid(im_part)==0) {

    // do nothing
   
   } else {
    amrex::Error("ns_is_rigid invalid");
   }

  } else
   amrex::Error("FSI_flag[im_part] invalid");

 } else
  amrex::Error("FSI_flag[im_part] invalid");

} // end subroutine copy_velocity_on_sign

// called from: FSI_make_distance, initData ()
void NavierStokes::build_moment_from_FSILS(Real cur_time) {

 std::string local_caller_string="build_moment_from_FSILS";

 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid 1");


 if (read_from_CAD()!=1)
  amrex::Error("read_from_CAD invalid");

 if (cur_time>=0.0) {
  // do nothing
 } else
  amrex::Error("cur_time invalid");

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 MultiFab& LS_new=get_new_data(LS_Type,slab_step+1);
 int nstate=STATE_NCOMP;
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");
 if (LS_new.nComp()!=num_materials*(1+AMREX_SPACEDIM))
  amrex::Error("LS_new.nComp()!=num_materials*(1+AMREX_SPACEDIM)");

 if (ngrow_make_distance!=3)
  amrex::Error("ngrow_make_distance!=3");
 int nparts=im_solid_map.size();
 if ((nparts<1)||(nparts>num_materials))
  amrex::Error("nparts invalid");
 int nFSI=nparts*NCOMP_FSI;
 if (localMF[FSI_MF]->nComp()!=nFSI)
  amrex::Error("localMF[FSI_MF]->nComp()!=nFSI");
 debug_ngrow(FSI_MF,ngrow_make_distance,local_caller_string);
  
 bool use_tiling=ns_tiling;

 const Real* dx = geom.CellSize();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(S_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(S_new,use_tiling); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();
  int bfact=parent->Space_blockingFactor(level);

  const Real* xlo = grid_loc[gridno].lo();

  FArrayBox& snewfab=S_new[mfi];
  FArrayBox& lsnewfab=LS_new[mfi];
  FArrayBox& fsifab=(*localMF[FSI_MF])[mfi];

  int tid_current=ns_thread();
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   //fort_build_moment is declared in GODUNOV_3D.F90
   //only copies from FSI_MF to state data if
   //fort_read_from_CAD(FSI_flag(im_part)).eq.1) 
  fort_build_moment(
    &cur_time,
    &level,
    &finest_level,
    &nFSI, 
    &nparts,
    &ngrow_make_distance, 
    im_solid_map.dataPtr(),
    xlo,dx,
    snewfab.dataPtr(),ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()),
    lsnewfab.dataPtr(),ARLIM(lsnewfab.loVect()),ARLIM(lsnewfab.hiVect()),
    fsifab.dataPtr(),ARLIM(fsifab.loVect()),ARLIM(fsifab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    &nstate);
 }  // mfi  
}//omp
 ns_reconcile_d_num(LOOP_BUILD_MOMENT,"build_moment_from_FSILS");

} // subroutine build_moment_from_FSILS


// 1. coarseTimeStep
// 2. timeStep
//     a. advance
//         i. CopyNewToOld
//         ii. setTimeLevel(time+dt_AMR,dt_AMR)
//     b. level_steps++
// 3. cumtime += dt_AMR
int NavierStokes::ok_copy_FSI_old_to_new() {

 if (level!=0)
  amrex::Error("level invalid ok_copy_FSI_old_to_new");

 int local_flag=1;
 if (FSI_interval==1) { //regenerate Eulerian data every time step.
  local_flag=0; 
 } else if (parent->levelSteps(0)==0) { //always regenerate at t=0
  local_flag=0;
 } else if (FSI_interval==0) { //never regenerate t>0
  // do nothing
 } else if (FSI_interval>=2) { //regenerate every FSI_interval steps.
  if (parent->levelSteps(0) % FSI_interval == 0)
   local_flag=0;
 } else
  amrex::Error("FSI_interval invalid");

 return local_flag;

} // end subroutine ok_copy_FSI_old_to_new

void NavierStokes::copy_old_FSI_to_new_level() {

 std::string local_caller_string="copy_old_FSI_to_new_level";

 if ((cur_time_slab>0.0)&&(cur_time_slab>prev_time_slab)) {
  // do nothing
 } else
  amrex::Error("expecting cur_time_slab>0, cur_time_slab>prev_time_slab");

 int nparts=im_solid_map.size();

 MultiFab& S_old=get_new_data(State_Type,slab_step);
 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 MultiFab* vofmf=
   getState(1,STATECOMP_MOF,num_materials*ngeom_raw,prev_time_slab);
 MultiFab& LS_old=get_new_data(LS_Type,slab_step);
 MultiFab& LS_new=get_new_data(LS_Type,slab_step+1);
 MultiFab* lsmf=getStateDist(1,prev_time_slab,local_caller_string);  
 MultiFab& Solid_old=get_new_data(Solid_State_Type,slab_step);
 MultiFab& Solid_new=get_new_data(Solid_State_Type,slab_step+1);
 MultiFab* velmf=getStateSolid(1,0,nparts*AMREX_SPACEDIM,prev_time_slab);

 for (int partid=0;partid<nparts;partid++) {
  int im_part=im_solid_map[partid];
  if ((im_part>=0)&&(im_part<num_materials)) {
   if (FSI_flag[im_part]==FSI_SHOELE_CTML) { //Wardlaw
    amrex::Error("must regenerate Eulerian data each step");
   } else if (FSI_flag[im_part]==FSI_PRESCRIBED_NODES){
    MultiFab::Copy(S_old,*vofmf,
      im_part*ngeom_raw,
      STATECOMP_MOF+im_part*ngeom_raw,
      ngeom_raw,1);
    MultiFab::Copy(S_new,*vofmf,
      im_part*ngeom_raw,
      STATECOMP_MOF+im_part*ngeom_raw,
      ngeom_raw,1);
    MultiFab::Copy(LS_old,*lsmf,
      im_part*(1+AMREX_SPACEDIM),
      im_part*(1+AMREX_SPACEDIM),
      (1+AMREX_SPACEDIM),1);
    MultiFab::Copy(LS_new,*lsmf,
      im_part*(1+AMREX_SPACEDIM),
      im_part*(1+AMREX_SPACEDIM),
      (1+AMREX_SPACEDIM),1);
    MultiFab::Copy(Solid_old,*velmf,
      partid*AMREX_SPACEDIM,
      partid*AMREX_SPACEDIM,
      AMREX_SPACEDIM,1);
    MultiFab::Copy(Solid_new,*velmf,
      partid*AMREX_SPACEDIM,
      partid*AMREX_SPACEDIM,
      AMREX_SPACEDIM,1);
   } else if (FSI_flag_valid(im_part)==1) {
    // do nothing
   } else
    amrex::Error("FSI_flag invalid"); 
  } else
   amrex::Error("im_part invalid");

 } // (int partid=0;partid<nparts;partid++)

 delete vofmf;
 delete lsmf;
 delete velmf;

} // end subroutine copy_old_FSI_to_new_level()

void NavierStokes::copy_old_FSI_to_new() {

 if (level!=0)
  amrex::Error("level invalid copy_old_FSI_to_new");

 if ((cur_time_slab>0.0)&&(cur_time_slab>prev_time_slab)) {
  // do nothing
 } else
  amrex::Error("expecting cur_time_slab>0, cur_time_slab>prev_time_slab");

 int finest_level=parent->finestLevel();
 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.copy_old_FSI_to_new_level();
 }

} // end subroutine copy_old_FSI_to_new

// called from: ns_header_msg_level,initData ()
void NavierStokes::Transfer_FSI_To_STATE(Real cur_time) {

 std::string local_caller_string="Transfer_FSI_To_STATE";

 if (ngrow_make_distance!=3)
  amrex::Error("ngrow_make_distance invalid");

 int nparts=im_solid_map.size();
 if ((nparts<0)||(nparts>num_materials))
  amrex::Error("nparts invalid");

  //FSI_PRESCRIBED_NODES
  //FSI_SHOELE_CTML
  //FSI_ICE_NODES_INIT
  //FSI_FLUID_NODES_INIT
 if (read_from_CAD()==1) {
  // do nothing
 } else
  amrex::Error("expecting read_from_CAD()==1");

 if ((nparts<1)||(nparts>num_materials))
  amrex::Error("nparts invalid");
 debug_ngrow(FSI_MF,ngrow_make_distance,local_caller_string);

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 int nstate=STATE_NCOMP;
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");

 MultiFab& Solid_new = get_new_data(Solid_State_Type,slab_step+1);
 if (Solid_new.nComp()!=nparts*AMREX_SPACEDIM)
  amrex::Error("Solid_new.nComp()!=nparts*AMREX_SPACEDIM");

 MultiFab& LS_new = get_new_data(LS_Type,slab_step+1);
 if (LS_new.nComp()!=num_materials*(AMREX_SPACEDIM+1))
  amrex::Error("LS_new invalid ncomp");
 int nFSI=nparts*NCOMP_FSI;
 if (localMF[FSI_MF]->nComp()!=nFSI)
  amrex::Error("localMF[FSI_MF]->nComp()!=nFSI");

 for (int partid=0;partid<nparts;partid++) {

  int im_part=im_solid_map[partid];
  if ((im_part<0)||(im_part>=num_materials))
   amrex::Error("im_part invalid");

  if (fort_read_from_CAD(&FSI_flag[im_part])==1) { 

   int ibase=partid*NCOMP_FSI;

   int ok_to_modify_EUL=1;
   if ((FSI_flag[im_part]==FSI_ICE_NODES_INIT)||  
       (FSI_flag[im_part]==FSI_FLUID_NODES_INIT)) { 
    if (cur_time==0.0) {
     ok_to_modify_EUL=1;
    } else if (cur_time>0.0) {
     ok_to_modify_EUL=0;
    } else
     amrex::Error("cur_time invalid");
   } else if ((FSI_flag[im_part]==FSI_PRESCRIBED_NODES)|| 
              (FSI_flag[im_part]==FSI_SHOELE_CTML)) {
    ok_to_modify_EUL=1;
    // do nothing
   } else
    amrex::Error("FSI_flag invalid");

   if (ok_to_modify_EUL==1) {

    copy_velocity_on_sign(partid);
    // Solid velocity
    //ngrow=0
    MultiFab::Copy(Solid_new,*localMF[FSI_MF],
     ibase+FSI_VELOCITY,partid*AMREX_SPACEDIM,
     AMREX_SPACEDIM,0);
     // LS
     //ngrow=0
    MultiFab::Copy(LS_new,*localMF[FSI_MF],ibase+FSI_LEVELSET,im_part,1,0);
     // temperature
    if (solidheat_flag==0) { // diffuse in solid
     // do nothing
    } else if ((solidheat_flag==1)||  //dirichlet
               (solidheat_flag==2)) { //neumann
      //ngrow=0
     MultiFab::Copy(S_new,*localMF[FSI_MF],ibase+FSI_TEMPERATURE,
      STATECOMP_STATES+im_part*num_state_material+ENUM_TEMPERATUREVAR,1,0);
    } else
     amrex::Error("solidheat_flag invalid"); 

   } else if (ok_to_modify_EUL==0) {
    // do nothing
   } else
    amrex::Error("ok_to_modify_EUL invalid");

  } else if (FSI_flag[im_part]==FSI_PRESCRIBED_PROBF90) { 
   // do nothing
  } else
   amrex::Error("FSI_flag invalid");

 } // partid=0..nparts-1

}  // subroutine Transfer_FSI_To_STATE

// called from: NavierStokes::post_restart(), NavierStokes::initData() 
void NavierStokes::init_aux_data() {

 int ioproc;
 if (ParallelDescriptor::IOProcessor())
  ioproc=1;
 else
  ioproc=0;

 if (level==0) {
  if (num_local_aux_grids>0) {
    //fort_init_aux_data() is declared in SOLIDFLUID.F90
   fort_init_aux_data(&ioproc);
  } else if (num_local_aux_grids==0) {
   // do nothing
  } else
   amrex::Error("num_local_aux_grids invalid");

 } else
  amrex::Error("expecting level==0 in init_aux_data()");

} // end subroutine init_aux_data()

//FSI_operation=OP_FSI_INITIALIZE_NODES  
//   initialize node locations; generate_new_triangles
//FSI_operation=OP_FSI_UPDATE_NODES  update node locations
//FSI_operation=OP_FSI_MAKE_DISTANCE  make distance in narrow band
//  (nparts x (vel, LS, Temp, flag, force)
//FSI_operation=OP_FSI_MAKE_SIGN  update the sign.
//FSI_operation=OP_FSI_LAG_STRESS  copy Eulerian stress to lag.
//
// note for CTML algorithm (FSI_flag==FSI_SHOELE_CTML):
// 0. Eulerian advection
// 1. Lagrangian advection
// 2. copy Lagrangian velocity and position to Eulerian.
// 3. copy Eulerian force to Lagrangian
//
// ns_header_msg_level is called from:
//  NavierStokes::FSI_make_distance
//  NavierStokes::post_restart
//  NavierStokes::initData
//  NavierStokes::nonlinear_advection
void NavierStokes::ns_header_msg_level(
 int FSI_operation,
 int FSI_sub_operation,
 Real cur_time,
 Real dt,
 int iter,
 const std::string& caller_string) {

 std::string local_caller_string="ns_header_msg_level";
 local_caller_string=caller_string+local_caller_string;

 int local_caller_id=caller_FSI_make_distance;

 const int max_level = parent->maxLevel();
 int finest_level=parent->finestLevel();

 if (pattern_test(local_caller_string,"FSI_make_distance")==1) {
  local_caller_id=caller_FSI_make_distance;
 } else if (pattern_test(local_caller_string,"post_restart")==1) {
  local_caller_id=caller_post_restart;
 } else if (pattern_test(local_caller_string,"initData")==1) {
  local_caller_id=caller_initData;
 } else if (pattern_test(local_caller_string,"nonlinear_advection")==1) {
  local_caller_id=caller_nonlinear_advection;
 } else
  amrex::Error("local_caller_string invalid");

 if (cur_time==cur_time_slab) {
  //do nothing
 } else
  amrex::Error("cur_time<>cur_time_slab");

 Vector< int > num_tiles_on_thread_proc;
 num_tiles_on_thread_proc.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  num_tiles_on_thread_proc[tid]=0;
 }

  //initialize node locations; generate_new_triangles
 if (FSI_operation==OP_FSI_INITIALIZE_NODES) { 

  if (level!=0)
   amrex::Error("level invalid");
  if (iter!=0)
   amrex::Error("iter invalid");
  if (FSI_sub_operation!=SUB_OP_FSI_DEFAULT)
   amrex::Error("FSI_sub_operation!=SUB_OP_FSI_DEFAULT");
  if ((local_caller_id==caller_post_restart)||
      (local_caller_id==caller_initData)) {
   //do nothing
  } else
   amrex::Error("local_caller_id invalid");

  //update node locations
 } else if (FSI_operation==OP_FSI_UPDATE_NODES) { 

  if (level!=0)
   amrex::Error("level invalid");
  if (iter!=0)
   amrex::Error("iter invalid");
  if (FSI_sub_operation!=SUB_OP_FSI_DEFAULT)
   amrex::Error("FSI_sub_operation!=SUB_OP_FSI_DEFAULT");
  if (CTML_FSI_flagC()==1) {
   if (num_divu_outer_sweeps<2)
    amrex::Error("num_divu_outer_sweeps<2, need to iterate w/tick");
  } else if (CTML_FSI_flagC()==0) {
   // do nothing
  } else
   amrex::Error("CTML_FSI_flagC() invalid");

  if (local_caller_id==caller_nonlinear_advection) {
   //do nothing
  } else
   amrex::Error("local_caller_id invalid");

  //make distance in narrow band
 } else if (FSI_operation==OP_FSI_MAKE_DISTANCE) { 

  if (iter!=0)
   amrex::Error("iter invalid");
  if (FSI_sub_operation!=SUB_OP_FSI_DEFAULT)
   amrex::Error("FSI_sub_operation!=SUB_OP_FSI_DEFAULT");
  if (local_caller_id==caller_FSI_make_distance) {
   //do nothing
  } else
   amrex::Error("local_caller_id invalid");

  //update the sign.
 } else if (FSI_operation==OP_FSI_MAKE_SIGN) { 

  if (iter<0)
   amrex::Error("iter invalid");
  if (FSI_sub_operation!=SUB_OP_FSI_DEFAULT)
   amrex::Error("FSI_sub_operation!=SUB_OP_FSI_DEFAULT");
  if (local_caller_id==caller_FSI_make_distance) {
   //do nothing
  } else
   amrex::Error("local_caller_id invalid");

  //copy Eulerian pressure to Lagrangian pressure
 } else if (FSI_operation==OP_FSI_LAG_STRESS) { 

  if (level==finest_level) {
   //do nothing
  } else
   amrex::Error("expecting level==finest_level if OP_FSI_LAG_STRESS");

  if (iter!=0)
   amrex::Error("iter invalid");

  if ((FSI_sub_operation==SUB_OP_FSI_CLEAR_LAG_DATA)||
      (FSI_sub_operation==SUB_OP_FSI_SYNC_LAG_DATA)||
      (FSI_sub_operation==SUB_OP_FSI_COPY_TO_LAG_DATA)) {
   //do nothing
  } else
   amrex::Error("FSI_sub_operation invalid for OP_FSI_LAG_STRESS");

  if (local_caller_id==caller_nonlinear_advection) {
   //do nothing
  } else
   amrex::Error("local_caller_id invalid");

 } else
  amrex::Error("FSI_operation invalid NavierStokes::ns_header_msg_level");

 for (int tid=0;tid<thread_class::nthreads;tid++) {
  FSI_touch_flag[tid]=0;
 }

 if ((level>max_level)||(finest_level>max_level))
  amrex::Error("(level>max_level)||(finest_level>max_level)");
 if ((level>max_level_for_use)||(finest_level>max_level_for_use))
  amrex::Error("(level>max_level_for_use)||(finest_level>max_level_for_use)");

 const Real* dx = geom.CellSize();
 Real h_small=dx[0];
 if (h_small>dx[1])
  h_small=dx[1];
 if (h_small>dx[AMREX_SPACEDIM-1])
  h_small=dx[AMREX_SPACEDIM-1];
 for (int i=level+1;i<=max_level;i++)
  h_small/=2.0;

 Real dx_max_level[AMREX_SPACEDIM];
 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  dx_max_level[dir]=dx[dir];
 for (int ilev=level+1;ilev<=max_level;ilev++) {
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   dx_max_level[dir]/=2.0;
  }
 }

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "ns_header_msg_level START\n";
   std::cout << "level= " << level << " finest_level= " << finest_level <<
    " max_level= " << max_level << '\n';
   std::cout << "FSI_operation= " << FSI_operation <<
    " cur_time = " << cur_time << " dt= " << dt << " iter = " << iter << '\n';
  }
 } else if (verbose==0) {
  // do nothing
 } else
  amrex::Error("verbose invalid");

 int ioproc;
 if (ParallelDescriptor::IOProcessor())
  ioproc=1;
 else
  ioproc=0;

 int current_step = nStep();
 int plot_interval=parent->plotInt();

  //ns_header_msg_level only called when read_from_CAD()==1
 if (read_from_CAD()==1) {

  int nparts=im_solid_map.size();
  if ((nparts<1)||(nparts>num_materials))
   amrex::Error("nparts invalid");

  MultiFab& Solid_new=get_new_data(Solid_State_Type,slab_step+1);
  if (Solid_new.nComp()!=nparts*AMREX_SPACEDIM)
   amrex::Error("Solid_new.nComp()!=nparts*AMREX_SPACEDIM");

  MultiFab& S_new=get_new_data(State_Type,slab_step+1);
  MultiFab& LS_new = get_new_data(LS_Type,slab_step+1);
  if (LS_new.nComp()!=num_materials*(AMREX_SPACEDIM+1))
   amrex::Error("LS_new invalid ncomp");

  bool use_tiling=ns_tiling;
  int bfact=parent->Space_blockingFactor(level);

  Real problo[AMREX_SPACEDIM];
  Real probhi[AMREX_SPACEDIM];
  Real problen[AMREX_SPACEDIM];
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   problo[dir]=geom.ProbLo(dir);
   probhi[dir]=geom.ProbHi(dir);
   problen[dir]=probhi[dir]-problo[dir];
   if (problen[dir]>0.0) {
    // do nothing
   } else
    amrex::Error("problen[dir] invalid");
  }

  FSI_container_class FSI_input;
  FSI_container_class FSI_output;
  Vector< Real > FSI_input_flattened;
  Vector< Real > FSI_output_flattened;

  int max_num_nodes_array[3];
  for (int dir=0;dir<3;dir++) {
   max_num_nodes_array[dir]=CTML_max_num_nodes_list[dir];
  }

  FSI_input.initData_FSI(CTML_FSI_numsolids,max_num_nodes_array,
    CTML_max_num_elements_list,CTML_FSI_num_scalars);
  FSI_output.initData_FSI(CTML_FSI_numsolids,max_num_nodes_array,
    CTML_max_num_elements_list,CTML_FSI_num_scalars);

  NavierStokes& ns_level0=getLevel(0);

   //initialize node locations; generate_new_triangles
  if (FSI_operation==OP_FSI_INITIALIZE_NODES) { 

   if (slab_step==ns_time_order-1) {
    //do nothing
   } else
    amrex::Error("expecting slab_step==ns_time_order-1");

   if (ns_level0.new_data_FSI[slab_step+1].CTML_num_solids==0) {
    //do nothing
   } else if (ns_level0.new_data_FSI[slab_step+1].CTML_num_solids==
 	      CTML_FSI_numsolids) {
    if (ns_level0.new_data_FSI[slab_step+1].max_num_nodes[0]==
        CTML_max_num_nodes_list[0]) {
     FSI_input.copyFrom_FSI(ns_level0.new_data_FSI[slab_step+1]);
     FSI_output.copyFrom_FSI(ns_level0.new_data_FSI[slab_step+1]);
    } else
     amrex::Error("new_data_FSI[slab_step+1].max_num_nodes incorrect");
   } else
    amrex::Error("new_data_FSI[slab_step+1].CTML_num_solids incorrect");

   //update node locations
  } else if (FSI_operation==OP_FSI_UPDATE_NODES) { 
    
   if ((slab_step>=0)&&(slab_step<ns_time_order)) {
    //do nothing
   } else
    amrex::Error("slab_step invalid");

   if (ns_level0.new_data_FSI[slab_step+1].CTML_num_solids==
       CTML_FSI_numsolids) {
    if (ns_level0.new_data_FSI[slab_step+1].max_num_nodes[0]==
        CTML_max_num_nodes_list[0]) {
     FSI_input.copyFrom_FSI(ns_level0.new_data_FSI[slab_step]);
     FSI_output.copyFrom_FSI(ns_level0.new_data_FSI[slab_step+1]);
    } else
     amrex::Error("new_data_FSI[slab_step+1].max_num_nodes incorrect");
   } else
    amrex::Error("new_data_FSI[slab_step+1].CTML_num_solids incorrect");

  } else if ((FSI_operation==OP_FSI_MAKE_DISTANCE)||
	     (FSI_operation==OP_FSI_MAKE_SIGN)||
             (FSI_operation==OP_FSI_LAG_STRESS)) {

   if ((slab_step>=0)&&(slab_step<ns_time_order)) {
    //do nothing
   } else
    amrex::Error("slab_step invalid");

   if (ns_level0.new_data_FSI[slab_step+1].CTML_num_solids==
       CTML_FSI_numsolids) {
    if ((ns_level0.new_data_FSI[slab_step+1].max_num_nodes[0]==
         CTML_max_num_nodes_list[0])&&
        (ns_level0.new_data_FSI[slab_step+1].max_num_nodes[1]==
         CTML_max_num_nodes_list[1])&&
        (ns_level0.new_data_FSI[slab_step+1].max_num_nodes[2]==
         CTML_max_num_nodes_list[2])) {
     FSI_input.copyFrom_FSI(ns_level0.new_data_FSI[slab_step+1]);
     FSI_output.copyFrom_FSI(ns_level0.new_data_FSI[slab_step+1]);
    } else
     amrex::Error("new_data_FSI[slab_step+1].max_num_nodes incorrect");
   } else
    amrex::Error("new_data_FSI[slab_step+1].CTML_num_solids incorrect");

  } else
   amrex::Error("FSI_operation invalid");

  FSI_input.FSI_flatten(FSI_input_flattened);
  FSI_output.FSI_flatten(FSI_output_flattened);
  int flatten_size=FSI_input_flattened.size();

  int nFSI=nparts*NCOMP_FSI;

  if ((FSI_operation==OP_FSI_INITIALIZE_NODES)||  // initialize nodes
      (FSI_operation==OP_FSI_UPDATE_NODES)) { // update node locations

   if (FSI_sub_operation!=SUB_OP_FSI_DEFAULT)
    amrex::Error("FSI_sub_operation!=SUB_OP_FSI_DEFAULT");

   int ngrow_make_distance_unitfab=0;
   IntVect unitlo(D_DECL(0,0,0));
   IntVect unithi(D_DECL(0,0,0));
    // construct cell-centered type box
   Box unitbox(unitlo,unithi);

   const int* tilelo=unitbox.loVect();
   const int* tilehi=unitbox.hiVect();
   const int* fablo=unitbox.loVect();
   const int* fabhi=unitbox.hiVect();

   FArrayBox FSIfab(unitbox,nFSI);

   Vector<int> velbc;
   velbc.resize(AMREX_SPACEDIM*2*AMREX_SPACEDIM);
   for (int i=0;i<velbc.size();i++)
    velbc[i]=0;
   Vector<int> vofbc;
   vofbc.resize(2*AMREX_SPACEDIM);
   for (int i=0;i<vofbc.size();i++)
    vofbc[i]=0;

   int tid=0;
   int gridno=0;
   int proc=ParallelDescriptor::MyProc();

    // fort_headermsg is declared in SOLIDFLUID.F90
   fort_headermsg(
      &proc,
      &tid,
      &num_tiles_on_thread_proc[tid],
      &gridno,
      &thread_class::nthreads,
      &level,
      &finest_level,
      &max_level,
      FSI_input_flattened.dataPtr(),
      FSI_output_flattened.dataPtr(),
      &flatten_size,
      &local_caller_id,
      &FSI_operation, // OP_FSI_INITIALIZE(UPDATE)_NODES
      &FSI_sub_operation, // SUB_OP_FSI_DEFAULT
      tilelo,tilehi,
      fablo,fabhi,
      &bfact,
      problo, //xlo
      problen, //dx
      dx_max_level, 
      velbc.dataPtr(),  
      vofbc.dataPtr(), 
      FSIfab.dataPtr(), // placeholder
      ARLIM(FSIfab.loVect()),ARLIM(FSIfab.hiVect()),
      FSIfab.dataPtr(), // drag spot
      ARLIM(FSIfab.loVect()),ARLIM(FSIfab.hiVect()),
      FSIfab.dataPtr(), // mnbrfab spot
      ARLIM(FSIfab.loVect()),ARLIM(FSIfab.hiVect()),
      FSIfab.dataPtr(), // mfiner spot
      ARLIM(FSIfab.loVect()),ARLIM(FSIfab.hiVect()),
      &nFSI,
      &ngrow_make_distance_unitfab,
      &nparts,
      im_solid_map.dataPtr(),
      &h_small,
      &cur_time, 
      &dt, 
      FSI_refine_factor.dataPtr(),
      FSI_bounding_box_ngrow.dataPtr(),
      &FSI_touch_flag[tid],
      &CTML_FSI_init,
      &iter,
      &current_step,
      &plot_interval,
      &ioproc);

   FSI_input.FSI_unflatten(FSI_input_flattened);
   FSI_output.FSI_unflatten(FSI_output_flattened);

   if (FSI_operation==OP_FSI_INITIALIZE_NODES) { 

    for (int i=0;i<=ns_time_order;i++) {
     ns_level0.new_data_FSI[i].copyFrom_FSI(FSI_output); 
    }

   } else if (FSI_operation==OP_FSI_UPDATE_NODES) { 

    ns_level0.new_data_FSI[slab_step+1].copyFrom_FSI(FSI_output); 

   } else
    amrex::Error("FSI_operation invalid");

   CTML_FSI_init=1;

  } else if ((FSI_operation==OP_FSI_MAKE_DISTANCE)|| 
             (FSI_operation==OP_FSI_MAKE_SIGN)) { 

   if (FSI_sub_operation!=SUB_OP_FSI_DEFAULT)
    amrex::Error("FSI_sub_operation!=SUB_OP_FSI_DEFAULT");

    // FSI_MF allocated in FSI_make_distance
   if (ngrow_make_distance!=3)
    amrex::Error("ngrow_make_distance invalid");
   debug_ngrow(FSI_MF,ngrow_make_distance,local_caller_string);
   if (localMF[FSI_MF]->nComp()!=nFSI)
    amrex::Error("localMF[FSI_MF]->nComp() invalid");

   if (FSI_operation==OP_FSI_MAKE_DISTANCE) { 

     // fill coarse patch 
    if (level>0) {

      //ngrow=0
     MultiFab* S_new_coarse=new MultiFab(grids,dmap,AMREX_SPACEDIM,0,
      MFInfo().SetTag("S_new_coarse"),FArrayBoxFactory());
     int dcomp=0;
     int scomp=0;
     FillCoarsePatch(*S_new_coarse,dcomp,cur_time,State_Type,
		     scomp,AMREX_SPACEDIM,debug_fillpatch);

     if (verbose>0) {
      if (ParallelDescriptor::IOProcessor()) {
       std::cout << "check_for_NAN(S_new_coarse)\n";
      }
      std::fflush(NULL);
      check_for_NAN(S_new_coarse);
     }

      //ngrow=0
     MultiFab* Solid_new_coarse=new MultiFab(grids,dmap,
	nparts*AMREX_SPACEDIM,0,
        MFInfo().SetTag("Solid_new_coarse"),FArrayBoxFactory());
     dcomp=0;
     scomp=0;

     if ((verbose>0)&&(1==0)) {
      if (ParallelDescriptor::IOProcessor()) {
       std::cout << "FillCoarsePatch(*Solid_new_coarse)\n";
      }
      std::fflush(NULL);
     }

     FillCoarsePatch(*Solid_new_coarse,dcomp,cur_time,Solid_State_Type,scomp,
        nparts*AMREX_SPACEDIM,debug_fillpatch);

     if (verbose>0) {
      if (ParallelDescriptor::IOProcessor()) {
       std::cout << "check_for_NAN(Solid_new_coarse)\n";
      }
      std::fflush(NULL);
      check_for_NAN(Solid_new_coarse);
     }

      //ngrow=0
     MultiFab* LS_new_coarse=
      new MultiFab(grids,dmap,num_materials*(AMREX_SPACEDIM+1),0,
                   MFInfo().SetTag("LS_new_coarse"),FArrayBoxFactory());
     dcomp=0;
     scomp=0;
     FillCoarsePatch(*LS_new_coarse,dcomp,cur_time,LS_Type,scomp,
        num_materials*(AMREX_SPACEDIM+1),debug_fillpatch);

     if (verbose>0) {
      if (ParallelDescriptor::IOProcessor()) {
       std::cout << "check_for_NAN(LS_new_coarse)\n";
      }
      std::fflush(NULL);
      check_for_NAN(LS_new_coarse);
     }

     for (int partid=0;partid<nparts;partid++) {

      int im_part=im_solid_map[partid];

      if ((im_part<0)||(im_part>=num_materials))
       amrex::Error("im_part invalid");

      if (fort_read_from_CAD(&FSI_flag[im_part])==1) { 

       int ok_to_modify_EUL=1;
       if ((FSI_flag[im_part]==FSI_ICE_NODES_INIT)||  
           (FSI_flag[im_part]==FSI_FLUID_NODES_INIT)) {
	if (cur_time==0.0) {
         ok_to_modify_EUL=1;
	} else if (cur_time>0.0) {
	 ok_to_modify_EUL=0;
	} else
	 amrex::Error("cur_time invalid");
       } else if ((FSI_flag[im_part]==FSI_PRESCRIBED_NODES)|| 
  	          (FSI_flag[im_part]==FSI_SHOELE_CTML)) { 
        ok_to_modify_EUL=1;
       } else
        amrex::Error("FSI_flag invalid");

       if (ok_to_modify_EUL==1) {

        dcomp=im_part;
        scomp=im_part;
         //ngrow==0 (levelset)
        MultiFab::Copy(LS_new,*LS_new_coarse,scomp,dcomp,1,0);
        dcomp=num_materials+im_part*AMREX_SPACEDIM;
        scomp=dcomp;
         //ngrow==0 (levelset normal)
        MultiFab::Copy(LS_new,*LS_new_coarse,scomp,dcomp,AMREX_SPACEDIM,0);

        dcomp=partid*AMREX_SPACEDIM;
        scomp=partid*AMREX_SPACEDIM;
         //ngrow==0
        MultiFab::Copy(Solid_new,*Solid_new_coarse,
		scomp,dcomp,AMREX_SPACEDIM,0);

        //ngrow==0
        MultiFab* new_coarse_thermal=new MultiFab(grids,dmap,1,0,
	    MFInfo().SetTag("new_coarse_thermal"),FArrayBoxFactory());
        dcomp=0;
        int scomp_thermal=STATECOMP_STATES+im_part*num_state_material+ 
		ENUM_TEMPERATUREVAR;
        //ncomp==1
        FillCoarsePatch(*new_coarse_thermal,dcomp,cur_time,State_Type,
         scomp_thermal,1,debug_fillpatch);

         //ngrow==0
        if (solidheat_flag==0) {  // diffuse in solid
         // do nothing
        } else if ((solidheat_flag==1)||   // dirichlet
                   (solidheat_flag==2)) {  // neumann
          //ngrow==0
         MultiFab::Copy(S_new,*new_coarse_thermal,0,scomp_thermal,1,0);
        } else
         amrex::Error("solidheat_flag invalid");

        delete new_coarse_thermal;

       } else if (ok_to_modify_EUL==0) {
        // do nothing
       } else
        amrex::Error("ok_to_modify_EUL invalid");

      } else if (FSI_flag[im_part]==FSI_PRESCRIBED_PROBF90) {  
       // do nothing
      } else
       amrex::Error("FSI_flag invalid");
     } // partid=0..nparts-1

     delete S_new_coarse;
     delete Solid_new_coarse;
     delete LS_new_coarse;

    } else if (level==0) {
     // do nothing
    } else
     amrex::Error("level invalid 3");

   } else if (FSI_operation==OP_FSI_MAKE_SIGN) { 
    // do not fill coarse patch.
   } else
    amrex::Error("FSI_operation invalid");

   MultiFab* solidmf=getStateSolid(ngrow_make_distance,0,
     nparts*AMREX_SPACEDIM,cur_time);
   MultiFab* denmf=getStateDen(ngrow_make_distance,cur_time);  
   MultiFab* LSMF=getStateDist(ngrow_make_distance,cur_time,
		   local_caller_string);
   if (LSMF->nGrow()!=ngrow_make_distance)
    amrex::Error("LSMF->nGrow()!=ngrow_make_distance");

   // FSI_MF allocated in FSI_make_distance
   // all components of FSI_MF are initialized to zero except for LS.
   // LS component of FSI_MF is init to -99999
   for (int partid=0;partid<nparts;partid++) {

    int im_part=im_solid_map[partid];
    if ((im_part<0)||(im_part>=num_materials))
     amrex::Error("im_part invalid");

    int ibase=partid*NCOMP_FSI;
     // velocity
    MultiFab::Copy(*localMF[FSI_MF],*solidmf,partid*AMREX_SPACEDIM,
      ibase+FSI_VELOCITY,AMREX_SPACEDIM,ngrow_make_distance);
     // LS  
    MultiFab::Copy(*localMF[FSI_MF],*LSMF,im_part,
      ibase+FSI_LEVELSET,1,ngrow_make_distance);
     // temperature
    MultiFab::Copy(*localMF[FSI_MF],*denmf,
      im_part*num_state_material+ENUM_TEMPERATUREVAR,
      ibase+FSI_TEMPERATURE,1,ngrow_make_distance);

     // flag (mask) and sign conflict indicator
    if (FSI_operation==OP_FSI_MAKE_DISTANCE) {

     // =0 no hits
     // =1 positive
     // =2 negative
     // =3 inconclusive
     setVal_localMF(FSI_MF,0.0,ibase+FSI_SIGN_CONFLICT,1,ngrow_make_distance); 

     if ((level>0)||
         ((level==0)&&(cur_time>0.0))) {
      setVal_localMF(FSI_MF,FSI_COARSE_LS_SIGN_VEL_VALID,
		     ibase+FSI_EXTRAP_FLAG,1,ngrow_make_distance); 
     } else if ((level==0)&&(cur_time==0.0)) {
      // do nothing
     } else
      amrex::Error("level or cur_time invalid");

    } else if (FSI_operation==OP_FSI_MAKE_SIGN) {
     // do nothing
    } else
     amrex::Error("FSI_operation invalid");

   } // partid=0..nparts-1

   // (1) =1 interior  =1 fine-fine ghost in domain  =0 otherwise
   // (2) =1 interior  =0 otherwise
   // (3) =1 interior+ngrow-1  =0 otherwise
   // (4) =1 interior+ngrow    =0 otherwise
   resize_mask_nbr(ngrow_make_distance);
   debug_ngrow(MASK_NBR_MF,ngrow_make_distance,local_caller_string);
 
   if (thread_class::nthreads<1)
    amrex::Error("thread_class::nthreads invalid");
   thread_class::init_d_numPts(S_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
   for (MFIter mfi(S_new,use_tiling); mfi.isValid(); ++mfi) {
    BL_ASSERT(grids[mfi.index()] == mfi.validbox());
    const int gridno = mfi.index();
    const Box& tilegrid = mfi.tilebox();
    const Box& fabgrid = grids[gridno];
    const int* tilelo=tilegrid.loVect();
    const int* tilehi=tilegrid.hiVect();
    const int* fablo=fabgrid.loVect();
    const int* fabhi=fabgrid.hiVect();
    const Real* xlo = grid_loc[gridno].lo();
    FArrayBox& FSIfab=(*localMF[FSI_MF])[mfi];
    FArrayBox& mnbrfab=(*localMF[MASK_NBR_MF])[mfi];

    Vector<int> velbc=getBCArray(Solid_State_Type,gridno,0,
     nparts*AMREX_SPACEDIM);
    Vector<int> vofbc=getBCArray(State_Type,gridno,STATECOMP_MOF,1);

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();
    int proc=ParallelDescriptor::MyProc();

    fort_headermsg(
     &proc,
     &tid_current,
     &num_tiles_on_thread_proc[tid_current],
     &gridno,
     &thread_class::nthreads,
     &level,
     &finest_level,
     &max_level,
     FSI_input_flattened.dataPtr(),
     FSI_output_flattened.dataPtr(),
     &flatten_size,
     &local_caller_id,
     &FSI_operation, // OP_FSI_MAKE_DISTANCE(SIGN)
     &FSI_sub_operation, // SUB_OP_FSI_DEFAULT
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,
     xlo,
     dx, 
     dx_max_level, 
     velbc.dataPtr(),  
     vofbc.dataPtr(), 
     FSIfab.dataPtr(),
     ARLIM(FSIfab.loVect()),ARLIM(FSIfab.hiVect()),
     FSIfab.dataPtr(), // drag spot
     ARLIM(FSIfab.loVect()),ARLIM(FSIfab.hiVect()),
     mnbrfab.dataPtr(),
     ARLIM(mnbrfab.loVect()),ARLIM(mnbrfab.hiVect()),
     mnbrfab.dataPtr(), // mfiner spot
     ARLIM(mnbrfab.loVect()),ARLIM(mnbrfab.hiVect()),
     &nFSI,
     &ngrow_make_distance,
     &nparts,
     im_solid_map.dataPtr(),
     &h_small,
     &cur_time, 
     &dt, 
     FSI_refine_factor.dataPtr(),
     FSI_bounding_box_ngrow.dataPtr(),
     &FSI_touch_flag[tid_current],
     &CTML_FSI_init,
     &iter,
     &current_step,
     &plot_interval,
     &ioproc);

    num_tiles_on_thread_proc[tid_current]++;
   } //mfi
}//omp
   ns_reconcile_d_num(LOOP_HEADER_MSG_MASK,"ns_header_msg_level");

   for (int tid=1;tid<thread_class::nthreads;tid++) {
    if (FSI_touch_flag[tid]==1) {
     FSI_touch_flag[0]=1;
    } else if (FSI_touch_flag[tid]==0) {
     // do nothing
    } else
     amrex::Error("FSI_touch_flag[tid] invalid");
   } 
   ParallelDescriptor::ReduceIntMax(FSI_touch_flag[0]);

   // idx,ngrow,scomp,ncomp,index,scompBC_map
   // InterpBordersGHOST is ultimately called.
   // dest_lstGHOST for Solid_State_Type defaults to pc_interp.
   // scompBC_map==0 corresponds to extrap_bc, pc_interp and fort_extrapfill
   // scompBC_map==1,2,3 corresponds to x or y or z vel_extrap_bc, pc_interp 
   //   and fort_extrapfill
   // nFSI=nparts * NCOMP_FSI
   for (int partid=0;partid<nparts;partid++) {
    int ibase=partid*NCOMP_FSI;
    Vector<int> scompBC_map;
    scompBC_map.resize(AMREX_SPACEDIM); 
    for (int dir=0;dir<AMREX_SPACEDIM;dir++)
     scompBC_map[dir]=dir+1;

    // This routine interpolates from coarser levels.
    PCINTERP_fill_borders(FSI_MF,ngrow_make_distance,ibase+FSI_VELOCITY,
     AMREX_SPACEDIM,Solid_State_Type,scompBC_map);

    for (int i=FSI_LEVELSET;i<NCOMP_FSI;i++) {
     scompBC_map.resize(1); 
     scompBC_map[0]=0;

     PCINTERP_fill_borders(FSI_MF,ngrow_make_distance,ibase+i,
      1,Solid_State_Type,scompBC_map);
    } // for (int i=FSI_LEVELSET;i<NCOMP_FSI;i++)
   } // partid=0..nparts-1

    // 1. copy_velocity_on_sign
    // 2. update Solid_new
    // 3. update LS_new
    // 4. update S_new(temperature) (if solidheat_flag==1 or 2)

   Transfer_FSI_To_STATE(cur_time);

   delete solidmf;
   delete denmf;
   delete LSMF;

   //copy Eulerian velocity and stress to Lagrangian velocity and stress
  } else if (FSI_operation==OP_FSI_LAG_STRESS) { 

   if (level==finest_level) {
    //do nothing
   } else
    amrex::Error("expecting level==finest_level");

   if (ngrow_make_distance!=3)
    amrex::Error("ngrow_make_distance invalid");
   if ((FSI_sub_operation!=SUB_OP_FSI_CLEAR_LAG_DATA)&&
       (FSI_sub_operation!=SUB_OP_FSI_COPY_TO_LAG_DATA)&&
       (FSI_sub_operation!=SUB_OP_FSI_SYNC_LAG_DATA))
    amrex::Error("FSI_sub_operation invalid");

   // (1) =1 interior  =1 fine-fine ghost in domain  =0 otherwise
   // (2) =1 interior  =0 otherwise
   // (3) =1 interior+ngrow-1  =0 otherwise
   // (4) =1 interior+ngrow    =0 otherwise
   resize_mask_nbr(ngrow_make_distance);
   debug_ngrow(MASK_NBR_MF,ngrow_make_distance,local_caller_string);
   // mask=1 if not covered or if outside the domain.
   // NavierStokes::maskfiner_localMF
   // NavierStokes::maskfiner
   resize_maskfiner(ngrow_make_distance,MASKCOEF_MF);
   debug_ngrow(MASKCOEF_MF,ngrow_make_distance,local_caller_string);

   if (localMF[DRAG_MF]->nGrow()!=ngrow_make_distance)
    amrex::Error("localMF[DRAG_MF] incorrect ngrow");
   if (localMF[DRAG_MF]->nComp()!=N_DRAG)
    amrex::Error("localMF[DRAG_MF] incorrect ncomp");

    //FSI_sub_operation:
    //  SUB_OP_FSI_CLEAR_LAG_DATA
    //  SUB_OP_FSI_COPY_TO_LAG_DATA
    //  SUB_OP_FSI_SYNC_LAG_DATA
    //
   if ((FSI_sub_operation==SUB_OP_FSI_CLEAR_LAG_DATA)|| 
       (FSI_sub_operation==SUB_OP_FSI_SYNC_LAG_DATA)) {

    if (FSI_sub_operation==SUB_OP_FSI_CLEAR_LAG_DATA) {

      // in: NavierStokes::ns_header_msg_level
     create_fortran_grid_struct(cur_time,dt);

    } else if (FSI_sub_operation==SUB_OP_FSI_SYNC_LAG_DATA) {

     // do nothing
	    
    } else
     amrex::Error("FSI_sub_operation bad: NavierStokes::ns_header_msg_level");

    int ngrow_make_distance_unitfab=0;
    IntVect unitlo(D_DECL(0,0,0));
    IntVect unithi(D_DECL(0,0,0));
     // construct cell-centered type box
    Box unitbox(unitlo,unithi);

    const int* tilelo=unitbox.loVect();
    const int* tilehi=unitbox.hiVect();
    const int* fablo=unitbox.loVect();
    const int* fabhi=unitbox.hiVect();

    FArrayBox FSIfab(unitbox,nFSI);

    Vector<int> velbc;
    velbc.resize(AMREX_SPACEDIM*2*AMREX_SPACEDIM);
    for (int i=0;i<velbc.size();i++)
     velbc[i]=0;
    Vector<int> vofbc;
    vofbc.resize(2*AMREX_SPACEDIM);
    for (int i=0;i<vofbc.size();i++)
     vofbc[i]=0;

    int tid=0;
    int gridno=0;
    int proc=ParallelDescriptor::MyProc();

    fort_headermsg(
     &proc,
     &tid,
     &num_tiles_on_thread_proc[tid],
     &gridno,
     &thread_class::nthreads,
     &level,
     &finest_level,
     &max_level,
     FSI_input_flattened.dataPtr(),
     FSI_output_flattened.dataPtr(),
     &flatten_size,
     &local_caller_id,
     &FSI_operation, //OP_FSI_LAG_STRESS
     &FSI_sub_operation, //SUB_OP_FSI_CLEAR_LAG_DATA or SYNC_LAG_DATA
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,
     problo,
     problen, 
     dx_max_level, 
     velbc.dataPtr(),  
     vofbc.dataPtr(), 
     FSIfab.dataPtr(), // placeholder
     ARLIM(FSIfab.loVect()),ARLIM(FSIfab.hiVect()),
     FSIfab.dataPtr(), // drag spot
     ARLIM(FSIfab.loVect()),ARLIM(FSIfab.hiVect()),
     FSIfab.dataPtr(), // mnbrfab spot
     ARLIM(FSIfab.loVect()),ARLIM(FSIfab.hiVect()),
     FSIfab.dataPtr(), // mfiner spot
     ARLIM(FSIfab.loVect()),ARLIM(FSIfab.hiVect()),
     &nFSI,
     &ngrow_make_distance_unitfab,
     &nparts,
     im_solid_map.dataPtr(),
     &h_small,
     &cur_time, 
     &dt, 
     FSI_refine_factor.dataPtr(),
     FSI_bounding_box_ngrow.dataPtr(),
     &FSI_touch_flag[tid],
     &CTML_FSI_init,
     &iter,
     &current_step,
     &plot_interval,
     &ioproc);

   } else if (FSI_sub_operation==SUB_OP_FSI_COPY_TO_LAG_DATA) {

    if (thread_class::nthreads<1)
     amrex::Error("thread_class::nthreads invalid");
    thread_class::init_d_numPts(localMF[MASKCOEF_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
    for (MFIter mfi(*localMF[MASKCOEF_MF],use_tiling); mfi.isValid(); ++mfi) {
     BL_ASSERT(grids[mfi.index()] == mfi.validbox());
     const int gridno = mfi.index();
     const Box& tilegrid = mfi.tilebox();
     const Box& fabgrid = grids[gridno];
     const int* tilelo=tilegrid.loVect();
     const int* tilehi=tilegrid.hiVect();
     const int* fablo=fabgrid.loVect();
     const int* fabhi=fabgrid.hiVect();
     const Real* xlo = grid_loc[gridno].lo();
     FArrayBox& FSIfab=(*localMF[DRAG_MF])[mfi]; //placeholder
     FArrayBox& dragfab=(*localMF[DRAG_MF])[mfi];//ngrow_make_dist ghost cells
     FArrayBox& mnbrfab=(*localMF[MASK_NBR_MF])[mfi];
     FArrayBox& mfinerfab=(*localMF[MASKCOEF_MF])[mfi];

     Vector<int> velbc=getBCArray(State_Type,gridno,
         STATECOMP_VEL,STATE_NCOMP_VEL);
     Vector<int> vofbc=getBCArray(State_Type,gridno,STATECOMP_MOF,1);

     int tid_current=ns_thread();
     if ((tid_current<0)||(tid_current>=thread_class::nthreads))
      amrex::Error("tid_current invalid");
     thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();
     int proc=ParallelDescriptor::MyProc();

     fort_headermsg(
      &proc,
      &tid_current,
      &num_tiles_on_thread_proc[tid_current],
      &gridno,
      &thread_class::nthreads,
      &level,
      &finest_level,
      &max_level,
      FSI_input_flattened.dataPtr(),
      FSI_output_flattened.dataPtr(),
      &flatten_size,
      &local_caller_id,
      &FSI_operation, //OP_FSI_LAG_STRESS
      &FSI_sub_operation, //SUB_OP_FSI_COPY_TO_LAG_DATA
      tilelo,tilehi,
      fablo,fabhi,
      &bfact,
      xlo,
      dx, 
      dx_max_level, 
      velbc.dataPtr(),  
      vofbc.dataPtr(), 
      FSIfab.dataPtr(), // placeholder
      ARLIM(FSIfab.loVect()),ARLIM(FSIfab.hiVect()),
      dragfab.dataPtr(), // ngrow_make_distance ghost cells DRAG_MF
      ARLIM(dragfab.loVect()),ARLIM(dragfab.hiVect()),
      mnbrfab.dataPtr(),
      ARLIM(mnbrfab.loVect()),ARLIM(mnbrfab.hiVect()),
      mfinerfab.dataPtr(),
      ARLIM(mfinerfab.loVect()),ARLIM(mfinerfab.hiVect()),
      &nFSI,
      &ngrow_make_distance,
      &nparts,
      im_solid_map.dataPtr(),
      &h_small,
      &cur_time, 
      &dt, 
      FSI_refine_factor.dataPtr(),
      FSI_bounding_box_ngrow.dataPtr(),
      &FSI_touch_flag[tid_current],
      &CTML_FSI_init,
      &iter,
      &current_step,
      &plot_interval,
      &ioproc);

     num_tiles_on_thread_proc[tid_current]++;
    } //mfi
}//omp
    ns_reconcile_d_num(LOOP_HEADERMSG_COPY_TO_LAG,"ns_header_msg_level");

   } else 
    amrex::Error("FSI_sub_operation invalid ns_header_msg_level");

  } else
   amrex::Error("FSI_operation invalid ns_header_msg_level");

  FSI_input.clear_FSI();
  FSI_output.clear_FSI();

 } else
  amrex::Error("read_from_CAD invalid");

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "ns_header_msg_level FINISH\n";
   std::cout << "level= " << level << " finest_level= " << finest_level <<
    " max_level= " << max_level << '\n';
   std::cout << "FSI_operation= " << FSI_operation <<
    " FSI_sub_operation= " << FSI_sub_operation <<
    " cur_time = " << cur_time << " dt= " << dt << " iter = " << iter << '\n';
  }
 } else if (verbose==0) {
  // do nothing
 } else
  amrex::Error("verbose invalid");

} // end subroutine ns_header_msg_level

// called from AmrCore::restart 
void NavierStokes::post_restart() {

 std::string local_caller_string="post_restart";

 ParmParse ppmain;
 Real local_stop_time=-1.0;
 ppmain.queryAdd("stop_time",local_stop_time);

 SDC_setup();
 ns_time_order=parent->Time_blockingFactor();
 slab_step=ns_time_order-1;

 SDC_outer_sweeps=0;
 SDC_setup_step();

 if (upper_slab_time==cur_time_slab) {
  // do nothing
 } else
  amrex::Error("expecting upper_slab_time==cur_time_slab");

 if (verbose>0)
  if (ParallelDescriptor::IOProcessor())
   std::cout << "in post_restart: level " << level << '\n';

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  std::cout << "dir mfiter_tile_size " << dir << ' ' <<
    FabArrayBase::mfiter_tile_size[dir] << '\n';
 }

 const int max_level = parent->maxLevel();
 const Real* dx = geom.CellSize();
 const Box& domain = geom.Domain();
 const int* domlo = domain.loVect();
 const int* domhi = domain.hiVect();

 Real problo[AMREX_SPACEDIM];
 Real probhi[AMREX_SPACEDIM];
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  problo[dir]=geom.ProbLo(dir);
  probhi[dir]=geom.ProbHi(dir);
 }

 MultiFab& S_new = get_new_data(State_Type,slab_step+1);
 int nc=S_new.nComp();

 fort_initdata_alloc(&nc,
  freezing_model.dataPtr(),
  distribute_from_target.dataPtr(),
  saturation_temp.dataPtr(),
  dx);

 if (level==0) {

  Vector<int> bfact_space_level(max_level+1);
  Vector<int> bfact_grid_level(max_level+1);
  for (int ilev=0;ilev<=max_level;ilev++) {
   bfact_space_level[ilev]=parent->Space_blockingFactor(ilev);
   bfact_grid_level[ilev]=parent->Old_blockingFactor(ilev);
  }
  int ioproc=ParallelDescriptor::IOProcessor();
   //fort_initgridmap is declared in PROB.F90
  fort_initgridmap(
    &verbose,
    &ioproc,
    &max_level,
    bfact_space_level.dataPtr(),
    bfact_grid_level.dataPtr(),
    domlo,domhi,
    dx,
    problo,probhi);

 } else if ((level>0)&&(level<=max_level)) {
  // do nothing
 } else {
  amrex::Error("level invalid post_restart() ");
 }

 metrics_data(2);  

 Real dt_amr=parent->getDt(); // returns dt_AMR

   // inside of post_restart
 if (level==0) {

  if (read_from_CAD()==1) {
   int iter=0;
   // in post_restart: initialize node locations; generate_new_triangles
   ns_header_msg_level(
    OP_FSI_INITIALIZE_NODES,
    SUB_OP_FSI_DEFAULT,
    upper_slab_time,
    dt_amr,
    iter,
    local_caller_string); 
  } else if (read_from_CAD()==0) {
   // do nothing
  } else
   amrex::Error("read_from_CAD invalid");

  // e.g. chkfile=./chk<nsteps>
  std::string chkfile=parent->theRestartFile();
  std::string Level_string = amrex::Concatenate("Level_", level, 1);
  std::string FullPath = chkfile;
  if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/') {
   FullPath += '/';
  }
  FullPath += Level_string;

  std::string FullPathName=FullPath;

  int time_order=parent->Time_blockingFactor();

  for (int i=0;i<=time_order;i++) {

#ifdef AMREX_PARTICLES

   AmrLevel0_new_dataPC[i] = std::make_unique<My_ParticleContainer>(parent);
    
   std::string Part_name="FusionPart";
   std::stringstream i_string_stream(std::stringstream::in |
      std::stringstream::out);
   i_string_stream << i;
   std::string i_string=i_string_stream.str();
   Part_name+=i_string;

   if (ParallelDescriptor::IOProcessor()) {
    std::cout << "Restarting particle container, time_slab i= " <<
      i << " FullPathName= " <<
      FullPathName << " Part_name= " << Part_name << '\n';
   }

   AmrLevel0_new_dataPC[i]->Restart(FullPathName,Part_name);

   Long num_particles=AmrLevel0_new_dataPC[i]->TotalNumberOfParticles();

   if (ParallelDescriptor::IOProcessor()) {
    std::cout << "TotalNumberOfParticles for slab i= " <<
      i << " is equal to " << num_particles << '\n';
   }
#endif

  }//for (int i=0;i<=time_order;i++) 

  init_aux_data();

  prepare_post_process(local_caller_string);

  if (sum_interval>0) {
   sum_integrated_quantities(local_caller_string,local_stop_time);
  }

 } else if ((level>0)&&(level<=max_level)) {
  // do nothing
 } else {
  amrex::Error("level invalid20");
 } 

}  // end subroutine post_restart


// This routine might be called twice at level 0 if AMR program
// run on more than one processor.  At level=0, the BoxArray used in
// "defbaselevel" can be different from the boxarray that optimizes
// load balancing.
void
NavierStokes::initData () {

 std::string local_caller_string="initData";

 int bfact_space=parent->Space_blockingFactor(level);
 int bfact_grid=parent->Old_blockingFactor(level);

 bool use_tiling=ns_tiling;

 if (ParallelDescriptor::IOProcessor()) {
  std::cout << "initData() at level= " << level << '\n';
  std::cout << "amr.space_blocking_factor= " << bfact_space << '\n';
  std::cout << "amr.blocking_factor= " << bfact_grid << '\n';
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   std::cout << "dir mfiter_tile_size " << dir << ' ' <<
     FabArrayBase::mfiter_tile_size[dir] << '\n';
  }
 }

 if (PCOPY_DEBUG==1) {
  debug_ParallelCopy();
 }

 SDC_setup();
  //AMReX_AmrCore.cpp; AmrCore::Time_blockingFactor () const
  //Time_blockingFactor returns "time_blocking_factor"
 ns_time_order=parent->Time_blockingFactor();
 slab_step=ns_time_order-1;

 SDC_outer_sweeps=0;
 SDC_setup_step();

 if (ngeom_raw!=AMREX_SPACEDIM+1)
  amrex::Error("ngeom_raw bust");

 int max_level = parent->maxLevel();

 const Real* dx = geom.CellSize();
 const Box& domain = geom.Domain();
 const int* domlo = domain.loVect();
 const int* domhi = domain.hiVect();
 Real problo[AMREX_SPACEDIM];
 Real probhi[AMREX_SPACEDIM];

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  problo[dir]=geom.ProbLo(dir);
  probhi[dir]=geom.ProbHi(dir);
 }

 if (upper_slab_time!=0.0)
  amrex::Error("upper_slab_time should be zero at the very beginning");
 if (cur_time_slab!=0.0)
  amrex::Error("cur_time_slab should be zero at the very beginning");

 if (level==0) {

  Vector<int> bfact_space_level(max_level+1);
  Vector<int> bfact_grid_level(max_level+1);
  for (int ilev=0;ilev<=max_level;ilev++) {
   bfact_space_level[ilev]=parent->Space_blockingFactor(ilev);
   bfact_grid_level[ilev]=parent->Old_blockingFactor(ilev);
  }
  int ioproc=ParallelDescriptor::IOProcessor();
   //fort_initgridmap is declared in PROB.F90
  fort_initgridmap(
   &verbose,
   &ioproc,
   &max_level,
   bfact_space_level.dataPtr(),
   bfact_grid_level.dataPtr(),
   domlo,domhi,
   dx,
   problo,probhi);

 } else if ((level>0)&&(level<=max_level)) {
  // do nothing
 } else {
  amrex::Error("level invalid 4");
 }

 metrics_data(1);

 if (level==0) {

  init_aux_data();

 } else if ((level>0)&&(level<=max_level)) {
  // do nothing
 } else {
  amrex::Error("level invalid 5");
 } 

 Real dt_amr=parent->getDt(); // returns dt_AMR

  // velocity,pres,state x num_materials,
  // interface variables x num_materials, error ind
 MultiFab& S_new = get_new_data(State_Type,slab_step+1);
 int nc=S_new.nComp();
 int nc_expect=STATE_NCOMP;
 if (nc!=nc_expect)
  amrex::Error("nc invalid in initdata");

 MultiFab& LS_new = get_new_data(LS_Type,slab_step+1);
 if (LS_new.nComp()!=num_materials*(AMREX_SPACEDIM+1))
  amrex::Error("LS_new invalid ncomp");

 MultiFab& DIV_new = get_new_data(DIV_Type,slab_step+1);
 if (DIV_new.nComp()!=1)
  amrex::Error("DIV_new.nComp()!=1");

 int nparts=im_solid_map.size();

 if ((nparts>=1)&&(nparts<=num_materials)) {  
  MultiFab& Solid_new = get_new_data(Solid_State_Type,slab_step+1);
  if (Solid_new.nComp()!=nparts*AMREX_SPACEDIM)
   amrex::Error("Solid_new.nComp()!=nparts*AMREX_SPACEDIM");
  Solid_new.setVal(0.0,0,nparts*AMREX_SPACEDIM,1);
 } else if (nparts==0) {
  // do nothing
 } else 
  amrex::Error("nparts invalid");

 int nparts_tensor=im_elastic_map.size();

 if ((nparts_tensor>=1)&&(nparts_tensor<=num_materials)) {  
  MultiFab& Tensor_new = get_new_data(Tensor_Type,slab_step+1);
  if (Tensor_new.nComp()!=NUM_CELL_ELASTIC)
   amrex::Error("Tensor_new.nComp()!=NUM_CELL_ELASTIC");
  Tensor_new.setVal(0.0,0,NUM_CELL_ELASTIC,1);
 } else if (nparts_tensor==0) {
  // do nothing
 } else 
  amrex::Error("nparts_tensor invalid");

 DIV_new.setVal(0.0);

 S_new.setVal(0.0,0,nc,1);
 LS_new.setVal(-99999.0,0,num_materials,1);
 LS_new.setVal(0.0,num_materials,num_materials*AMREX_SPACEDIM,1); // slopes

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  MultiFab& Smac_new = get_new_data(Umac_Type+dir,slab_step+1);

  if (Smac_new.nComp()!=1) {
   std::cout << "num_materials = " << num_materials << '\n';
   amrex::Error("Smac_new.nComp() invalid in initData");
  }
  Smac_new.setVal(0.0,0,1,0);

 }  // dir=0..sdim-1

 prepare_mask_nbr(1);

  //FSI_PRESCRIBED_NODES
  //FSI_SHOELE_CTML
  //FSI_ICE_NODES_INIT
  //FSI_FLUID_NODES_INIT
 if (read_from_CAD()==1) {
  int iter=0; 
  // in initData: initialize node locations; generate_new_triangles
  if (level==0) {
   ns_header_msg_level(
    OP_FSI_INITIALIZE_NODES,
    SUB_OP_FSI_DEFAULT,
    upper_slab_time,
    dt_amr,
    iter,
    local_caller_string); 
  }

  // create a distance function (velocity and temperature) on this level.
  // calls ns_header_msg_level with 
  //   FSI_operation==OP_FSI_MAKE_DISTANCE(SIGN)
  // ns_header_msg_level calls NavierStokes::Transfer_FSI_To_STATE
  FSI_make_distance(upper_slab_time,dt_amr);

  //FSI_FLUID
  //FSI_PRESCRIBED_PROBF90
  //FSI_ICE_PROBF90
  //FSI_ICE_STATIC
  //FSI_RIGID_NOTPRESCRIBED
 } else if (read_from_CAD()==0) {
  // do nothing
 } else
  amrex::Error("read_from_CAD invalid");

  //regenerate_from_eulerian calls fort_initdatasolid
 regenerate_from_eulerian(upper_slab_time);

 fort_initdata_alloc(&nc,
  freezing_model.dataPtr(),
  distribute_from_target.dataPtr(),
  saturation_temp.dataPtr(),
  dx);

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "initdata loop follows->level=" << level << 
     " max_level= " << max_level << '\n';
  }
 }
 
 MultiFab* lsmf=getStateDist(1,upper_slab_time,local_caller_string);
 MultiFab::Copy(LS_new,*lsmf,0,0,num_materials*(1+AMREX_SPACEDIM),1);
 delete lsmf;

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(S_new.boxArray().d_numPts());
 
#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(S_new,use_tiling); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();
  const Real* xlo = grid_loc[gridno].lo();
  const Real* xhi = grid_loc[gridno].hi();

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  fort_initdata(
   &tid_current,
   &adapt_quad_depth,
   &level,&max_level,
   &upper_slab_time,
   tilelo,tilehi,
   fablo,fabhi,
   &bfact_space,
   &nc,
   saturation_temp.dataPtr(),
   S_new[mfi].dataPtr(),
   ARLIM(S_new[mfi].loVect()),
   ARLIM(S_new[mfi].hiVect()),
   LS_new[mfi].dataPtr(),
   ARLIM(LS_new[mfi].loVect()),
   ARLIM(LS_new[mfi].hiVect()),
   dx,xlo,xhi,
   centroid_noise_factor.dataPtr()); 

  if (1==0) {
   FArrayBox& snewfab=S_new[mfi];
   int interior_only=1;
   tecplot_debug(snewfab,xlo,fablo,fabhi,dx,-1,0,STATECOMP_STATES, 
    num_materials*num_state_material,interior_only); 
  }
 } //mfi
}//omp
 ns_reconcile_d_num(LOOP_INITDATA,"initData");

 if (read_from_CAD()==1) {
  build_moment_from_FSILS(upper_slab_time);
 } else if (read_from_CAD()==0) {
  //do nothing
 } else
  amrex::Error("read_from_CAD() invalid");

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(S_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(S_new,use_tiling); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());

  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();

  const Real* xlo = grid_loc[gridno].lo();
  const Real* xhi = grid_loc[gridno].hi();

  const int* s_lo = S_new[mfi].loVect();
  const int* s_hi = S_new[mfi].hiVect();

  Real Reynolds=visc_coef*viscconst[0];
  if (Reynolds>0.0) Reynolds=1.0/Reynolds;
  Real Weber=tension[0];
  if (Weber>0.0) Weber=1.0/Weber;
  Real RGASRWATER=probhi[0];
  if (xblob>0.0) RGASRWATER/=xblob;

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  fort_initvelocity(
   &level,&upper_slab_time,
   tilelo,tilehi,
   fablo,fabhi,&bfact_space,
   S_new[mfi].dataPtr(),
   ARLIM(s_lo),ARLIM(s_hi),
   dx,xlo,xhi);

 } // mfi
}//omp
 ns_reconcile_d_num(LOOP_INITVELOCITY,"initData");

 if (read_from_CAD()==1) {
  // 1. copy_velocity_on_sign
  // 2. update Solid_new
  // 3. update LS_new
  // 4. update S_new(temperature) (if solidheat_flag==1 or 2)
  Transfer_FSI_To_STATE(upper_slab_time);
 } else if (read_from_CAD()==0) {
  // do nothing
 } else
  amrex::Error("read_from_CAD invalid");


 level_species_reaction(local_caller_string);

  // if nparts>0,
  //  Initialize FSI_GHOST_MAC_MF from Solid_State_Type
  // Otherwise initialize FSI_GHOST_MAC_MF with the fluid velocity.
  // No law of the wall modeling.
  // FSI_GHOST_MAC_MF needed in NavierStokes::post_init_state
 init_FSI_GHOST_MAC_MF_predict();

 init_regrid_history();
 is_first_step_after_regrid=-1;

 int nstate=state.size();
 if (nstate!=NUM_STATE_TYPE)
  amrex::Error("nstate invalid");

 for (int k=0;k<nstate;k++) {
  state[k].CopyNewToOld(level,max_level);  // olddata=newdata 
   // time_array[0]=upper_slab_time-dt_amr  
   // time_array[bfact_time_order]=upper_slab_time (=0.0)
  state[k].setTimeLevel(upper_slab_time,dt_amr); 
 }

#ifdef AMREX_PARTICLES
 NavierStokes& ns_level0=getLevel(0);
 int lev_min=0;
 int lev_max=level;

 My_ParticleContainer& current_PC=ns_level0.newDataPC(ns_time_order);
 int nGrow_Redistribute=0;
 int local_Redistribute=0; //redistribute "from scratch"
 current_PC.Redistribute(lev_min,lev_max,nGrow_Redistribute, 
    local_Redistribute);
 ns_level0.CopyNewToOldPC(lev_max); //copy and redistribute up to lev_max
#endif

}  // end subroutine initData

void NavierStokes::init_boundary_list(Vector<int> scomp,
  Vector<int> ncomp) {

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

 int ncomp_list=0;
 for (int ilist=0;ilist<scomp.size();ilist++) {
  MultiFab* state_mf=getState(1,scomp[ilist],ncomp[ilist],cur_time_slab);
  MultiFab::Copy(S_new,*state_mf,0,scomp[ilist],ncomp[ilist],1);
  delete state_mf;
  ncomp_list+=ncomp[ilist];
 }
 if (ncomp_list<=0)
  amrex::Error("ncomp_list invalid");

} // subroutine init_boundary_list

void NavierStokes::init_boundary() {

 std::string local_caller_string="init_boundary";
 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int nstate=state.size();
 if (nstate!=NUM_STATE_TYPE)
  amrex::Error("nstate invalid");

 int nden=num_materials*num_state_material;

 for (int k=0;k<nstate;k++) {

  if (k==State_Type) {
   MultiFab& S_new=get_new_data(State_Type,slab_step+1);
   MultiFab* vofmf=getState(1,STATECOMP_MOF,
      num_materials*ngeom_raw,cur_time_slab);
   MultiFab::Copy(S_new,*vofmf,0,STATECOMP_MOF,num_materials*ngeom_raw,1);
   delete vofmf;
   MultiFab* velmf=getState(1,STATECOMP_VEL,
	STATE_NCOMP_VEL+STATE_NCOMP_PRES,cur_time_slab);
   MultiFab::Copy(S_new,*velmf,0,0,STATE_NCOMP_VEL+STATE_NCOMP_PRES,1);
   delete velmf;
   MultiFab* denmf=getStateDen(1,cur_time_slab);  
   MultiFab::Copy(S_new,*denmf,0,STATECOMP_STATES,nden,1);
   delete denmf;
  } else if (k==Umac_Type) {
   // do nothing
  } else if (k==Vmac_Type) {
   // do nothing
  } else if (k==Wmac_Type) {
   // do nothing
  } else if (k==LS_Type) {
   MultiFab& LS_new=get_new_data(LS_Type,slab_step+1);
   MultiFab* lsmf=getStateDist(1,cur_time_slab,local_caller_string); 
   MultiFab::Copy(LS_new,*lsmf,0,0,num_materials*(1+AMREX_SPACEDIM),1);
   delete lsmf;
  } else if (k==DIV_Type) {
   // do nothing
  } else if ((k==Solid_State_Type)&&
	     (im_solid_map.size()>=1)&&
	     (im_solid_map.size()<=num_materials)) {
   int nparts=im_solid_map.size();
   if ((nparts<1)||(nparts>num_materials))
    amrex::Error("nparts invalid");
   MultiFab& Solid_new=get_new_data(Solid_State_Type,slab_step+1);
   MultiFab* velmf=getStateSolid(1,0,nparts*AMREX_SPACEDIM,cur_time_slab);
   MultiFab::Copy(Solid_new,*velmf,0,0,nparts*AMREX_SPACEDIM,1);
   delete velmf;
  } else if ((k==Tensor_Type)&&
             (num_materials_viscoelastic>=1)&&
             (num_materials_viscoelastic<=num_materials)) {
   int nparts=im_elastic_map.size();
   if ((nparts<=0)||(nparts>num_materials))
    amrex::Error("nparts invalid");
   if (nparts!=num_materials_viscoelastic)
    amrex::Error("nparts!=num_materials_viscoelastic");
   MultiFab& Tensor_new=get_new_data(Tensor_Type,slab_step+1);
     // ngrow=1 scomp=0
   MultiFab* tensormf=getStateTensor(1,0,
     NUM_CELL_ELASTIC,cur_time_slab);
   MultiFab::Copy(Tensor_new,*tensormf,0,0,
     NUM_CELL_ELASTIC,1);
   delete tensormf;
  } else if ((k==TensorX_Type)&&
             (num_materials_viscoelastic>=1)&&
             (num_materials_viscoelastic<=num_materials)) {
   // do nothing
  } else if ((k==TensorY_Type)&&
             (num_materials_viscoelastic>=1)&&
             (num_materials_viscoelastic<=num_materials)) {
   // do nothing
  } else if ((k==TensorZ_Type)&&
             (num_materials_viscoelastic>=1)&&
             (num_materials_viscoelastic<=num_materials)) {
   // do nothing
  } else 
   amrex::Error("k invalid");

 } // k=0..nstate-1

}  // subroutine init_boundary()


//  AmrLevel* a = (*levelbld)(*this,lev,geom[lev],
//    new_grid_places[lev],new_dmap[lev],cumtime);
void
NavierStokes::init(
  AmrLevel & old,
  const BoxArray& ba_in,  // BoxArray of "this" (new amr_level)
  const DistributionMapping& dmap_in) { // dmap of "this" (new amr_level)
 
 const int max_level = parent->maxLevel();

 NavierStokes* oldns     = (NavierStokes*) &old;

 SDC_setup();
 ns_time_order=parent->Time_blockingFactor();
 slab_step=ns_time_order-1;

 SDC_outer_sweeps=0;

 upper_slab_time=oldns->state[State_Type].slabTime(ns_time_order);
 lower_slab_time=oldns->state[State_Type].slabTime(0);
 cur_time_slab=upper_slab_time;
 prev_time_slab=oldns->state[State_Type].slabTime(ns_time_order-1);

  // in: init(&old,ba_in,dmap_in)
 dt_slab=cur_time_slab-prev_time_slab;
 delta_slab_time=upper_slab_time-lower_slab_time;

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "init(old)  level,cur_time,dt " << level << ' ' <<
     upper_slab_time << ' ' << delta_slab_time << '\n';
  }
 }

  // new time: upper_slab_time   old time: upper_slab_time-delta_slab_time
 setTimeLevel(upper_slab_time,delta_slab_time); 

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 if (ngeom_raw!=AMREX_SPACEDIM+1)
  amrex::Error("ngeom_raw bust");

 int nstate=state.size();  // cell centered, vel MAC, LS, etc
 if (nstate!=NUM_STATE_TYPE)
  amrex::Error("nstate invalid");

 for (int local_index=0;local_index<nstate;local_index++) {

  int state_holds_data=get_desc_lst()[local_index].get_state_holds_data();

  if (state_holds_data==1) {

    // data with respect to new grid structure.
   MultiFab &S_new = get_new_data(local_index,ns_time_order);

   if (S_new.DistributionMap()==dmap_in) {
    // do nothing
   } else {
    amrex::Error("dmap_in invalid in init(old)");
   }
    //  Are the BoxArrays equal after conversion to cell-centered?
   if (S_new.boxArray().CellEqual(ba_in)) {
    // do nothing
   } else {
    amrex::Error("S_new.boxArray().CellEqual(ba_in) failed init(old)");
   }

   int ncomp=S_new.nComp();

   int numparts=1;
   int ncomp_part[4];
   int scomp_part[4];
   
   if (local_index==State_Type) {

    int test_ncomp=0;

    numparts=4;
    scomp_part[0]=0;
    scomp_part[1]=STATE_NCOMP_VEL+STATE_NCOMP_PRES;
    scomp_part[2]=scomp_part[1]+num_materials*num_state_material;
    scomp_part[3]=scomp_part[2]+num_materials*ngeom_raw;

    ncomp_part[0]=scomp_part[1];
    test_ncomp+=ncomp_part[0];
    ncomp_part[1]=scomp_part[2]-scomp_part[1];
    test_ncomp+=ncomp_part[1];
    ncomp_part[2]=scomp_part[3]-scomp_part[2];
    test_ncomp+=ncomp_part[2];
    ncomp_part[3]=1;
    test_ncomp+=ncomp_part[3];
    if (test_ncomp!=ncomp)
     amrex::Error("test ncomp invalid");

   } else if (local_index!=State_Type) {

    numparts=1;
    scomp_part[0]=0;
    ncomp_part[0]=ncomp;

   } else {
    amrex::Error("local_index bust");
   }
    
   for (int part_iter=0;part_iter<numparts;part_iter++) { 
     //FillPatch is declared in amrlib/AmrLevel.cpp
    FillPatch(
       old,
       S_new,
       scomp_part[part_iter],
       upper_slab_time,
       local_index,
       scomp_part[part_iter],
       ncomp_part[part_iter],
       debug_fillpatch);
   }

  } else if (state_holds_data==0) {
   // do nothing
  } else
   amrex::Error("state_holds_data invalid");

 } // local_index=0..nstate-1

 if (level==0) {
  for (int i=0;i<=ns_time_order;i++) {
   new_data_FSI[i].copyFrom_FSI(oldns->new_data_FSI[ns_time_order]);
  }
 }

#ifdef AMREX_PARTICLES

 int lev_min=0;
 int lev_max=level;
 int nGrow_Redistribute=0;
 bool local_copy=true; //do not redistribute inside of copyParticles
 int local_redistribute=0;

  // if level==0, we must copy from the old Amr_level (level==0) to the
  // new Amr_level prior to deleting the old level==0 structure.
 if (level==0) {

  My_ParticleContainer& new_PC=newDataPC(ns_time_order);
  My_ParticleContainer& old_PC=oldns->newDataPC(ns_time_order);

  Long old_num_particles=old_PC.TotalNumberOfParticles();
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "OLD: TotalNumberOfParticles for slab ns_time_order= " <<
     ns_time_order << " is equal to " << old_num_particles << '\n';
  }

  new_PC.clearParticles();
  //make sure hierarchy is initialized.
  new_PC.Redistribute();

  new_PC.copyParticles(old_PC,local_copy);

  Long new_num_particles=new_PC.TotalNumberOfParticles();
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "NEW: TotalNumberOfParticles for slab ns_time_order= " <<
     ns_time_order << " is equal to " << new_num_particles << '\n';
  }

  new_PC.Redistribute(lev_min,lev_max,nGrow_Redistribute, 
   local_redistribute);

 } else if ((level>0)&&(level<=max_level)) {

  NavierStokes& ns_level0=getLevel(0);
  My_ParticleContainer& new_PC=ns_level0.newDataPC(ns_time_order);
  new_PC.Redistribute(lev_min,lev_max,nGrow_Redistribute, 
   local_redistribute);
 
 } else
  amrex::Error("level invalid");
   
#endif

 old_intersect_new = amrex::intersect(grids,oldns->boxArray());
 is_first_step_after_regrid = 1;

 debug_fillpatch=0;

}  // end subroutine init(old)

// init a new level that did not exist on the previous step.
// NavierStokes::init is called from: AmrCore::regrid
//  AmrLevel* a = (*levelbld)(*this,lev,geom[lev],
//    new_grid_places[lev],new_dmap[lev],cumtime);
void
NavierStokes::init(
  const BoxArray& ba_in,  // BoxArray of "this" (new amr_level)
  const DistributionMapping& dmap_in) { // dmap of "this" (new amr_level)

 if (level==0)
  amrex::Error("this init only called for level>0");

 Real dt_amr=parent->getDt();
 Real dt_new=dt_amr;
 parent->setDt(dt_new);

 NavierStokes& old = getLevel(level-1);

 SDC_setup();
 ns_time_order=parent->Time_blockingFactor();
 slab_step=ns_time_order-1;

 SDC_outer_sweeps=0;

 upper_slab_time = old.state[State_Type].slabTime(ns_time_order);
 lower_slab_time = old.state[State_Type].slabTime(0);
 cur_time_slab=upper_slab_time;
 prev_time_slab=old.state[State_Type].slabTime(ns_time_order-1);

  // in: init (ba_in,dmap_in)
 dt_slab=cur_time_slab-prev_time_slab;
 delta_slab_time=upper_slab_time-lower_slab_time;

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "init()  level,cur_time,dt " << level << ' ' <<
     upper_slab_time << ' ' << delta_slab_time << '\n';
  }
 }
  
  // new time: upper_slab_time   old time: upper_slab_time-delta_slab_time
 setTimeLevel(upper_slab_time,delta_slab_time);

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 if (ngeom_raw!=AMREX_SPACEDIM+1)
  amrex::Error("ngeom_raw bust");

 int nstate=state.size(); // cell centered, vel MAC, LS, etc
 if (nstate!=NUM_STATE_TYPE)
  amrex::Error("nstate invalid");

 for (int local_index=0;local_index<nstate;local_index++) {

  int state_holds_data=get_desc_lst()[local_index].get_state_holds_data();

  if (state_holds_data==1) {

   MultiFab &S_new = get_new_data(local_index,ns_time_order);

   if (S_new.DistributionMap()==dmap_in) {
    // do nothing
   } else {
    amrex::Error("dmap_in invalid in init()");
   }

    //  Are the BoxArrays equal after conversion to cell-centered?
   if (S_new.boxArray().CellEqual(ba_in)) {
    // do nothing
   } else {
    amrex::Error("S_new.boxArray().CellEqual(ba_in) failed init()");
   }

   int ncomp=S_new.nComp();

   int numparts=1;
   int ncomp_part[4];
   int scomp_part[4];
   
   if (local_index==State_Type) {

    int test_ncomp=0;

    numparts=4;
    scomp_part[0]=0;
    scomp_part[1]=STATE_NCOMP_VEL+STATE_NCOMP_PRES;
    scomp_part[2]=scomp_part[1]+num_materials*num_state_material;
    scomp_part[3]=scomp_part[2]+num_materials*ngeom_raw;

    ncomp_part[0]=scomp_part[1];
    test_ncomp+=ncomp_part[0];
    ncomp_part[1]=scomp_part[2]-scomp_part[1];
    test_ncomp+=ncomp_part[1];
    ncomp_part[2]=scomp_part[3]-scomp_part[2];
    test_ncomp+=ncomp_part[2];
    ncomp_part[3]=1;
    test_ncomp+=ncomp_part[3];
    if (test_ncomp!=ncomp)
     amrex::Error("test ncomp invalid");

   } else if (local_index!=State_Type) {

    numparts=1;
    scomp_part[0]=0;
    ncomp_part[0]=ncomp;

   } else {
    amrex::Error("local_index bust");
   }

   for (int part_iter=0;part_iter<numparts;part_iter++) { 
    FillCoarsePatch(
      S_new,scomp_part[part_iter],
      upper_slab_time,
      local_index,
      scomp_part[part_iter],
      ncomp_part[part_iter],
      debug_fillpatch);
   }

  } else if (state_holds_data==0) {
   // do nothing
  } else
   amrex::Error("state_holds_data invalid");

 } // local_index=0..nstate-1

#ifdef AMREX_PARTICLES

 int lev_min=0;
 int lev_max=level;
 int nGrow_Redistribute=0;
 int local_Redistribute=0;

 if (level==0) {
  amrex::Error("level==0 cannot be made from nothing");
 } else if (level>0) {
  // do nothing
 } else
  amrex::Error("level invalid");
   
 NavierStokes& ns_level0=getLevel(0);
 My_ParticleContainer& new_PC=ns_level0.newDataPC(ns_time_order);

 new_PC.Redistribute(lev_min,lev_max,nGrow_Redistribute, 
   local_Redistribute);

#endif

 init_regrid_history();
 is_first_step_after_regrid = 2;

}  // end subroutine init()

void NavierStokes::CopyNewToOldALL() {

 int finest_level=parent->finestLevel();
 const int max_level = parent->maxLevel();

 if (level!=0)
  amrex::Error("expecting level=0 in CopyNewToOldALL");

 if (ns_time_order==parent->Time_blockingFactor()) {
  //do nothing
 } else
  amrex::Error("expecting ns_time_order==parent->Time_blockingFactor()");

 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
   // oldtime=newtime newtime+=dt
  int nstate=ns_level.state.size();
  if (nstate!=NUM_STATE_TYPE)
   amrex::Error("nstate invalid");
   //copy bfact_time_order component to the
   //components: 0..bfact_time_order-1
  for (int k = 0; k < nstate; k++) {
   ns_level.state[k].CopyNewToOld(level,max_level);  
  }
 }// for (int ilev=level;ilev<=finest_level;ilev++) 

 for (int i=0;i<ns_time_order;i++) {
  new_data_FSI[i].copyFrom_FSI(new_data_FSI[ns_time_order]); 
 }

#ifdef AMREX_PARTICLES
 NavierStokes& ns_level0=getLevel(0);
 int lev_min=0;
 int lev_max=finest_level;

 My_ParticleContainer& current_PC=ns_level0.newDataPC(ns_time_order);
 int nGrow_Redistribute=0;
 int local_Redistribute=0; //redistribute "from scratch"
 current_PC.Redistribute(lev_min,lev_max,nGrow_Redistribute, 
    local_Redistribute);
 ns_level0.CopyNewToOldPC(lev_max); 

 Long num_particles=current_PC.TotalNumberOfParticles();

 if (ParallelDescriptor::IOProcessor()) {
  std::cout<<"COPYNEWOLD: TotalNumberOfParticles for slab ns_time_order= "<<
   ns_time_order << " is equal to " << num_particles << '\n';
 }
#endif

}  // subroutine CopyNewToOldALL


void NavierStokes::CopyOldToNewALL() {

 int finest_level=parent->finestLevel();
 const int max_level = parent->maxLevel();

 if (level!=0)
  amrex::Error("level invalid CopyOldToNewALL");

 if (ns_time_order==parent->Time_blockingFactor()) {
  //do nothing
 } else
  amrex::Error("expecting ns_time_order==parent->Time_blockingFactor()");

 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  int nstate=ns_level.state.size();
  if (nstate!=NUM_STATE_TYPE)
   amrex::Error("nstate invalid");
   //copy component "0" to components 1..bfact_time_order.
  for (int k = 0; k < nstate; k++) {
   ns_level.state[k].CopyOldToNew(level,max_level);
  }
 }// for (int ilev=level;ilev<=finest_level;ilev++)

 for (int i=1;i<=ns_time_order;i++) {
  new_data_FSI[i].copyFrom_FSI(new_data_FSI[0]); 
 }

#ifdef AMREX_PARTICLES
 NavierStokes& ns_level0=getLevel(0);
 int lev_min=0;
 int lev_max=finest_level;

 My_ParticleContainer& current_PC=ns_level0.newDataPC(0);
 int nGrow_Redistribute=0;
 int local_Redistribute=0; //redistribute "from scratch"
 current_PC.Redistribute(lev_min,lev_max,nGrow_Redistribute, 
    local_Redistribute);
 ns_level0.CopyOldToNewPC(lev_max);

 Long num_particles=current_PC.TotalNumberOfParticles();

 if (ParallelDescriptor::IOProcessor()) {
  std::cout << "OLD: TotalNumberOfParticles for slab zero= " <<
      " is equal to " << num_particles << '\n';
 }

#endif

} // subroutine CopyOldToNewALL


void NavierStokes::Number_CellsALL(Real& rcells) {

 int finest_level=parent->finestLevel();
 if (level!=0)
  amrex::Error("level invalid Number_CellsALL");
 rcells=0.0;
 for (int ilev=0;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  int ngrids=ns_level.grids.size();
  for (int gridno=0;gridno<ngrids;gridno++) {
   const Box& gridbox=ns_level.grids[gridno];
   rcells+=gridbox.numPts();
  }
 }
   
} // subroutine Number_CellsALL  

// called from:
//  NavierStokes::do_the_advance  (near beginning)
//  NavierStokes::prepare_post_process (near beginning)
void NavierStokes::allocate_mdot() {

 int nsolve=1;

 delete_localMF_if_exist(MDOT_MF,1);

  // MDOT has nsolve components.
 new_localMF(MDOT_MF,nsolve,0,-1);
 setVal_localMF(MDOT_MF,0.0,0,nsolve,0); //val,scomp,ncomp,ngrow

} // end subroutine allocate_mdot()

// slab_step and SDC_outer_sweeps are set before calling this routine.
void
NavierStokes::SDC_setup_step() {

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 if ((num_materials<1)||(num_materials>1000))
  amrex::Error("num_materials out of range");

 ns_time_order=parent->Time_blockingFactor();

 if ((SDC_outer_sweeps<0)||
     (SDC_outer_sweeps>=ns_time_order))
  amrex::Error("SDC_outer_sweeps invalid");

 divu_outer_sweeps=0;

 lower_slab_time=state[State_Type].slabTime(0);
 upper_slab_time=state[State_Type].slabTime(ns_time_order);
 delta_slab_time=upper_slab_time-lower_slab_time;
 if (delta_slab_time>=0.0) {
  // do nothing
 } else {
  std::cout << "lower_slab_time= " << lower_slab_time << '\n';
  std::cout << "upper_slab_time= " << upper_slab_time << '\n';
  std::cout << "delta_slab_time= " << delta_slab_time << '\n';
  amrex::Error("delta_slab_time invalid");
 }
 if (delta_slab_time>0.0) {
  // do nothing
 } else if (delta_slab_time==0.0) {
  if (upper_slab_time==0.0) {
   // do nothing
  } else {
   std::cout << "lower_slab_time= " << lower_slab_time << '\n';
   std::cout << "upper_slab_time= " << upper_slab_time << '\n';
   std::cout << "delta_slab_time= " << delta_slab_time << '\n';
   amrex::Error("delta_slab_time or upper_slab_time invalid");
  }
 } else {
  std::cout << "lower_slab_time= " << lower_slab_time << '\n';
  std::cout << "upper_slab_time= " << upper_slab_time << '\n';
  std::cout << "delta_slab_time= " << delta_slab_time << '\n';
  amrex::Error("delta_slab_time invalid");
 }

 if ((slab_step>=0)&&(slab_step<ns_time_order)) {

  prev_time_slab=state[State_Type].slabTime(slab_step);
  cur_time_slab=state[State_Type].slabTime(slab_step+1);
  vel_time_slab=cur_time_slab;
  prescribed_vel_time_slab=cur_time_slab;

   // in: NavierStokes::SDC_setup_step()
  dt_slab=cur_time_slab-prev_time_slab;
  if (dt_slab>=0.0) {
   // do nothing
  } else
   amrex::Error("dt_slab invalid1");

  if ((dt_slab==0.0)&&(cur_time_slab==0.0)) {
   // do nothing
  } else if (dt_slab>0.0) {
   // do nothing
  } else {
   amrex::Error("dt_slab invalid2");
  }

 } else if ((slab_step==-1)&&(ns_time_order>=2)) {
  prev_time_slab=lower_slab_time;
  cur_time_slab=lower_slab_time;
  dt_slab=1.0;
  vel_time_slab=cur_time_slab;
  prescribed_vel_time_slab=cur_time_slab;
 } else if ((slab_step==ns_time_order)&&(ns_time_order>=2)) {
  prev_time_slab=upper_slab_time;
  cur_time_slab=upper_slab_time;
  dt_slab=1.0;
  vel_time_slab=cur_time_slab;
  prescribed_vel_time_slab=cur_time_slab;
 } else
  amrex::Error("slab_step invalid"); 

} // SDC_setup_step()


void
NavierStokes::SDC_setup() {

 SDC_outer_sweeps=0;
 divu_outer_sweeps=0;
 prev_time_slab=0.0;
 cur_time_slab=0.0;
 vel_time_slab=0.0;
 prescribed_vel_time_slab=0.0;
 dt_slab=1.0;

 upper_slab_time=0.0;
 lower_slab_time=0.0;
 delta_slab_time=1.0;

} // subroutine SDC_setup()

void 
NavierStokes::Geometry_setup() {

 std::cout.precision(20);

 for (int i=0;i<MAX_NUM_LOCAL_MF;i++) {
  localMF_grow[i]=-1;
  localMF[i]=0;  // null pointer
 }

 SDC_setup();
 ns_time_order=1;
 slab_step=0;

}

void 
NavierStokes::Geometry_cleanup() {

 for (int i=0;i<MAX_NUM_LOCAL_MF;i++) {
  if (localMF_grow[i]>=0) {
   delete localMF[i];
   ParallelDescriptor::Barrier();
   localMF_grow[i]=-1;
   localMF[i]=0;
   ParallelDescriptor::Barrier();
  } else if (localMF_grow[i]==-1) {
   if (localMF[i]==0) {
    // do nothing
   } else {
    amrex::Error("localMF[i] invalid");
   }
  } else {
   amrex::Error("localMF_grow[i] invalid");
  }
 } // i=0 ... MAX_NUM_LOCAL_MF-1   

} // subroutine Geometry_cleanup()


void NavierStokes::make_viscoelastic_tensorMACALL(int im,
  int flux_mf,int flux_grid_type,
  int fill_state_idx) {

 int finest_level=parent->finestLevel();

 std::string local_caller_string="make_viscoelastic_tensorMACALL";

 if ((fill_state_idx==TensorX_Type)||
     (fill_state_idx==TensorY_Type)||
     ((fill_state_idx==TensorZ_Type)&&(AMREX_SPACEDIM==3))) {
  //do nothing
 } else
  amrex::Error("fill_state_idx invalid");

 if (level==0) {
  // do nothing
 } else
  amrex::Error("level invalid");

 if ((flux_grid_type==0)|| //X MAC
     (flux_grid_type==1)|| //Y MAC
     ((flux_grid_type==2)&&(AMREX_SPACEDIM==3))) { //Z MAC
  // do nothing
 } else
  amrex::Error("flux_grid_type invalid");

 int push_enable_spectral=enable_spectral;
 int elastic_enable_spectral=0;
 override_enable_spectral(elastic_enable_spectral);

 if (localMF_grow[VISCOTEN_MF]==1) {
  // do nothing
 } else 
  amrex::Error("VISCOTEN_MF should be allocated");

 if (localMF_grow[flux_mf]==-1) {
  // do nothing
 } else 
  amrex::Error("flux_mf should not be allocated");

  //ngrow=0  localMF[flux_mf] initialized to 0.0.
 allocate_array(0,ENUM_NUM_TENSOR_TYPE,flux_grid_type,flux_mf);

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.make_viscoelastic_tensorMAC(im,
     flux_mf,flux_grid_type,fill_state_idx);
 }

 if (localMF_grow[VISCOTEN_MF]==1) {
  // do nothing
 } else 
  amrex::Error("VISCOTEN_MF has incorrect Ngrow");

 if (localMF_grow[flux_mf]==0) {
  // do nothing
 } else 
  amrex::Error("flux_mf has incorrect Ngrow");

 if (ENUM_NUM_TENSOR_TYPE==2*AMREX_SPACEDIM) {
  // do nothing
 } else
  amrex::Error("ENUM_NUM_TENSOR_TYPE invalid");

 if (localMF[VISCOTEN_MF]->nComp()==ENUM_NUM_TENSOR_TYPE) {
  // do nothing
 } else 
  amrex::Error("VISCOTEN_MF has incorrect nComp");

 if (localMF[flux_mf]->nComp()==ENUM_NUM_TENSOR_TYPE) {
  // do nothing
 } else 
  amrex::Error("flux_mf has incorrect nComp");

 if (localMF[flux_mf]->nGrow()==0) {
  // do nothing
 } else 
  amrex::Error("flux_mf has incorrect nGrow");

 int Q_grid_type=-1;
 debug_boxArray(localMF[VISCOTEN_MF],Q_grid_type,local_caller_string);
 debug_boxArray(localMF[flux_mf],flux_grid_type,local_caller_string);

  // spectral_override==0 => always low order.
 for (int ilev=finest_level-1;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  int spectral_override=LOW_ORDER_AVGDOWN;//always do low order average down
   //declared in NavierStokes2.cpp
  ns_level.avgDownEdge_localMF(flux_mf,
    0,ENUM_NUM_TENSOR_TYPE,
    flux_grid_type,-1,
    spectral_override,
    local_caller_string);
 } // ilev=finest_level-1 ... level

 for (int scomp_extrap=0;scomp_extrap<ENUM_NUM_TENSOR_TYPE;scomp_extrap++) {
  Vector<int> scompBC_map;
   // desc_lstGHOST.setComponent(Tensor_Type, ...
   // "set_tensor_bc", tensor_pc_interp 
   // fort_extrapfill
   // (i.e. the coarse/fine BC and physical BC will be low order)
  scompBC_map.resize(1);
  scompBC_map[0]=scomp_extrap;
   // idx,ngrow,scomp,ncomp,index,scompBC_map
  PCINTERP_fill_bordersALL(flux_mf,0,scomp_extrap,1,
	fill_state_idx,scompBC_map);
 } // scomp_extrap=0..ENUM_NUM_TENSOR_TYPE-1

 if (1==0) {
  writeSanityCheckData(
   "FLUX_MF",
   "T11,T12,T22,T33,T13,T23",
   local_caller_string,
   flux_mf, //tower_mf_id
   localMF[flux_mf]->nComp(), 
   flux_mf,
   -1, //State_Type==-1
   flux_grid_type,
   parent->levelSteps(0)); 
 }

 override_enable_spectral(push_enable_spectral);

} // end subroutine make_viscoelastic_tensorMACALL

void NavierStokes::make_viscoelastic_tensorMAC(int im,
  int flux_mf,int flux_grid_type,int fill_state_idx) {

 std::string local_caller_string="make_viscoelastic_tensorMAC";

 int finest_level=parent->finestLevel();
 bool use_tiling=ns_tiling;

 if ((flux_grid_type==0)|| //X MAC
     (flux_grid_type==1)|| //Y MAC
     ((flux_grid_type==2)&&(AMREX_SPACEDIM==3))) { //Z MAC
  // do nothing
 } else
  amrex::Error("flux_grid_type invalid");

 IndexType local_typ(get_desc_lstGHOST()[fill_state_idx].getType());
 int flux_box_type[AMREX_SPACEDIM];
  // NavierStokes::grid_type_to_box_type_cpp declared in NavierStokes2.cpp
 grid_type_to_box_type_cpp(flux_grid_type,flux_box_type);
 for (int dir_local=0;dir_local<AMREX_SPACEDIM;dir_local++) {
  if ((flux_box_type[dir_local]==0)&&
      (local_typ.cellCentered(dir_local)==true)) {
   // do nothing
  } else if ((flux_box_type[dir_local]==1)&&
             (local_typ.nodeCentered(dir_local)==true)) {
   // do nothing
  } else
   amrex::Error("flux_box_type and local_typ are inconsistent");
 } // dir_local=0..sdim-1

 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=num_materials)) {
  // do nothing
 } else
  amrex::Error("num_materials_viscoelastic bad:make_viscoelastic_tensorMAC");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 if (localMF_grow[VISCOTEN_MF]==1) {
  // do nothing
 } else 
  amrex::Error("VISCOTEN_MF should be allocated");

 if (localMF[VISCOTEN_MF]->nGrow()==1) {
  // do nothing
 } else 
  amrex::Error("VISCOTEN_MF invalid nGrow");

 if (localMF[VISCOTEN_MF]->nComp()==ENUM_NUM_TENSOR_TYPE) {
  // do nothing
 } else
  amrex::Error("VISCOTEN_MF invalid nComp");

 if (localMF_grow[flux_mf]==0) {
  // do nothing
 } else 
  amrex::Error("flux_mf should be allocated");

 if (localMF[flux_mf]->nGrow()==0) {
  // do nothing
 } else
  amrex::Error("flux_mf invalid nGrow");

 if (localMF[flux_mf]->nComp()==ENUM_NUM_TENSOR_TYPE) {
  // do nothing
 } else
  amrex::Error("flux_mf invalid nComp");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  debug_ngrow(FACE_VAR_MF+dir,0,local_caller_string);
 }
 debug_ngrow(CELL_VISC_MATERIAL_MF,1,local_caller_string);

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

 // 1. viscosity coefficient - 1..num_materials
 // 2. viscoelastic coefficient - 1..num_materials
 // 3. relaxation time - 1..num_materials
 // the viscous and viscoelastic forces should both be multiplied by
 // visc_coef.  
 if (localMF[CELL_VISC_MATERIAL_MF]->nComp()!=3*num_materials) {
  std::cout << "ncomp= " <<
   localMF[CELL_VISC_MATERIAL_MF]->nComp() << 
   " num_materials= " << num_materials << '\n';
  amrex::Error("cell_visc_material ncomp invalid(1)");
 }
 if (localMF[CELL_VISC_MATERIAL_MF]->nGrow()<1)
  amrex::Error("cell_visc_material ngrow invalid");

 int nstate=STATE_NCOMP;
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");

 if ((im<0)||(im>=num_materials))
  amrex::Error("im invalid52");

 if (ns_is_rigid(im)==0) {

  if (store_elastic_data[im]==1) {

   int partid=0;
   while ((im_elastic_map[partid]!=im)&&(partid<im_elastic_map.size())) {
    partid++;
   }

   if (partid<im_elastic_map.size()) {

    if (ENUM_NUM_TENSOR_TYPE!=2*AMREX_SPACEDIM)
     amrex::Error("ENUM_NUM_TENSOR_TYPE invalid");

    const Real* dx = geom.CellSize();

    if (thread_class::nthreads<1)
     amrex::Error("thread_class::nthreads invalid");
    thread_class::init_d_numPts(localMF[VISCOTEN_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
    for (MFIter mfi(*localMF[VISCOTEN_MF],use_tiling); mfi.isValid(); ++mfi) {
     BL_ASSERT(grids[mfi.index()] == mfi.validbox());
     const int gridno = mfi.index();
     const Box& tilegrid = mfi.tilebox();
     const Box& fabgrid = grids[gridno];
     const int* tilelo=tilegrid.loVect();
     const int* tilehi=tilegrid.hiVect();
     const int* fablo=fabgrid.loVect();
     const int* fabhi=fabgrid.hiVect();
     int bfact=parent->Space_blockingFactor(level);

     const Real* xlo = grid_loc[gridno].lo();

     FArrayBox& tenfab=(*localMF[VISCOTEN_MF])[mfi];
     FArrayBox& tenMACfab=(*localMF[flux_mf])[mfi];
     if (tenMACfab.box().ixType()==local_typ) {
      // do nothing
     } else
      amrex::Error("tenMACfab.box().ixType()!=local_typ");

     // 1. maketensor: TQ_{m}=alpha_{m} Q_{m}
     // 2. tensor force: F= div (H_{m} TQ_{m})
     //    H=H(phi_biased)

     FArrayBox& viscfab=(*localMF[CELL_VISC_MATERIAL_MF])[mfi];
     int ncomp_visc=viscfab.nComp();

     int tid_current=ns_thread();
     if ((tid_current<0)||(tid_current>=thread_class::nthreads))
      amrex::Error("tid_current invalid");
     thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

       // fort_maketensor_mac is declared in: GODUNOV_3D.F90
       // viscoelastic_model==0 => FENE-CR
       // viscoelastic_model==1 => Oldroyd B
       // viscoelastic_model==3 => incremental elastic model
       // viscoelastic_model==5 => FENE-P
       // viscoelastic_model==6 => linear PTT
       // viscoelastic_model==7 => incremental Neo-Hookean model
     if (fort_built_in_elastic_model(&elastic_viscosity[im],
			           &viscoelastic_model[im])==1) {
      fort_maketensor_mac(
       &flux_grid_type,
       &partid,
       &level,
       &finest_level,
       &ncomp_visc,
       &im,  // 0..num_materials-1
       xlo,dx,
       viscfab.dataPtr(),ARLIM(viscfab.loVect()),ARLIM(viscfab.hiVect()),
       tenfab.dataPtr(),ARLIM(tenfab.loVect()),ARLIM(tenfab.hiVect()),
       tenMACfab.dataPtr(),
       ARLIM(tenMACfab.loVect()),ARLIM(tenMACfab.hiVect()),
       tilelo,tilehi,
       fablo,fabhi,
       &bfact,
       &elastic_viscosity[im],
       &etaS[im],
       &elastic_time[im],
       &viscoelastic_model[im],
       &polymer_factor[im],
       &NS_geometry_coord);
     } else
      amrex::Error("fort_built_in_elastic_model invalid");
    }  // mfi  
}//omp
    ns_reconcile_d_num(LOOP_MAKETENSOR_MAC,"make_viscoelastic_tensorMAC");

    check_for_NAN(localMF[VISCOTEN_MF]);

   } else
    amrex::Error("partid could not be found: make_viscoelastic_tensor");

  } else 
   amrex::Error("expecting (store_elastic_data[im]==1) ");

 } else if (ns_is_rigid(im)==1) {
  amrex::Error("expecting ns_is_rigid(im)==0");
 } else
  amrex::Error("ns_is_rigid invalid");

}  // subroutine make_viscoelastic_tensorMAC

// called from:
//  NavierStokes::GetDragALL
//  NavierStokes::veldiffuseALL
//  NavierStokes::vel_elastic_ALL
void NavierStokes::make_viscoelastic_tensorALL(int im) {

 int finest_level=parent->finestLevel();

 if (level==0) {
  // do nothing
 } else
  amrex::Error("level invalid");


 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=num_materials)) {
  // do nothing
 } else
  amrex::Error("num_materials_viscoelastic bad:make_viscoelastic_tensorALL");

 int push_enable_spectral=enable_spectral;
 int elastic_enable_spectral=0;
 override_enable_spectral(elastic_enable_spectral);

 if (localMF_grow[VISCOTEN_MF]==-1) {
  // do nothing
 } else 
  amrex::Error("VISCOTEN_MF should not be allocated");

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.make_viscoelastic_tensor(im);
 }
 if (localMF_grow[VISCOTEN_MF]==1) {
  // do nothing
 } else 
  amrex::Error("VISCOTEN_MF has incorrect Ngrow");

 if (ENUM_NUM_TENSOR_TYPE==2*AMREX_SPACEDIM) {
  // do nothing
 } else
  amrex::Error("ENUM_NUM_TENSOR_TYPE invalid");

 if (localMF[VISCOTEN_MF]->nComp()==ENUM_NUM_TENSOR_TYPE) {
  // do nothing
 } else 
  amrex::Error("VISCOTEN_MF has incorrect nComp");

  // spectral_override==0 => always low order.
 avgDown_localMF_ALL(VISCOTEN_MF,0,ENUM_NUM_TENSOR_TYPE,LOW_ORDER_AVGDOWN);
 for (int scomp_extrap=0;scomp_extrap<ENUM_NUM_TENSOR_TYPE;scomp_extrap++) {
  Vector<int> scompBC_map;
   // desc_lstGHOST.setComponent(Tensor_Type, ...
   // "set_tensor_bc", tensor_pc_interp 
   // fort_extrapfill
   // (i.e. the coarse/fine BC and physical BC will be low order)
  scompBC_map.resize(1);
  scompBC_map[0]=scomp_extrap;
   // idx,ngrow,scomp,ncomp,index,scompBC_map
  PCINTERP_fill_bordersALL(VISCOTEN_MF,1,scomp_extrap,1,
	Tensor_Type,scompBC_map);
 } // scomp_extrap=0..ENUM_NUM_TENSOR_TYPE-1

 override_enable_spectral(push_enable_spectral);

} // end subroutine make_viscoelastic_tensorALL

void NavierStokes::make_viscoelastic_tensor(int im) {

 std::string local_caller_string="make_viscoelastic_tensor";

 int finest_level=parent->finestLevel();
 bool use_tiling=ns_tiling;

 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=num_materials)) {
  // do nothing
 } else
  amrex::Error("num_materials_viscoelastic invalid:make_viscoelastic_tensor");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 if (localMF_grow[VISCOTEN_MF]==-1) {
  // do nothing
 } else 
  amrex::Error("VISCOTEN_MF should not be allocated");

 debug_ngrow(CELL_VISC_MATERIAL_MF,1,local_caller_string);

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

 // 1. viscosity coefficient - 1..num_materials
 // 2. viscoelastic coefficient - 1..num_materials
 // 3. relaxation time - 1..num_materials
 // the viscous and viscoelastic forces should both be multiplied by
 // visc_coef.  
 if (localMF[CELL_VISC_MATERIAL_MF]->nComp()!=3*num_materials) {
  std::cout << "ncomp= " <<
   localMF[CELL_VISC_MATERIAL_MF]->nComp() << 
   " num_materials= " << num_materials << '\n';
  amrex::Error("cell_visc_material ncomp invalid(1)");
 }
 if (localMF[CELL_VISC_MATERIAL_MF]->nGrow()<1)
  amrex::Error("cell_visc_material ngrow invalid");

 int nstate=STATE_NCOMP;
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");

 if ((im<0)||(im>=num_materials))
  amrex::Error("im invalid52");

 if (ns_is_rigid(im)==0) {

  if (store_elastic_data[im]==1) {

   int partid=0;
   while ((im_elastic_map[partid]!=im)&&(partid<im_elastic_map.size())) {
    partid++;
   }

   if (partid<im_elastic_map.size()) {

    int scomp_tensor=partid*ENUM_NUM_TENSOR_TYPE;

    if (ENUM_NUM_TENSOR_TYPE!=2*AMREX_SPACEDIM)
     amrex::Error("ENUM_NUM_TENSOR_TYPE invalid");

     // VISCOTEN_MF will be used by NavierStokes::make_viscoelastic_heating
     //  or
     // VISCOTEN_MF will be used by NavierStokes::GetDrag
     //  or
     // VISCOTEN_MF will be used by NavierStokes::make_viscoelastic_force
     //
    getStateTensor_localMF(VISCOTEN_MF,1,scomp_tensor,ENUM_NUM_TENSOR_TYPE,
     cur_time_slab);

    const Real* dx = geom.CellSize();

    if (thread_class::nthreads<1)
     amrex::Error("thread_class::nthreads invalid");
    thread_class::init_d_numPts(localMF[VISCOTEN_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
    for (MFIter mfi(*localMF[VISCOTEN_MF],use_tiling); mfi.isValid(); ++mfi) {
     BL_ASSERT(grids[mfi.index()] == mfi.validbox());
     const int gridno = mfi.index();
     const Box& tilegrid = mfi.tilebox();
     const Box& fabgrid = grids[gridno];
     const int* tilelo=tilegrid.loVect();
     const int* tilehi=tilegrid.hiVect();
     const int* fablo=fabgrid.loVect();
     const int* fabhi=fabgrid.hiVect();
     int bfact=parent->Space_blockingFactor(level);

     const Real* xlo = grid_loc[gridno].lo();

     FArrayBox& tenfab=(*localMF[VISCOTEN_MF])[mfi];
     // 1. maketensor: TQ_{m}=alpha_{m} Q_{m}
     // 2. tensor force: F= div (H_{m} TQ_{m})
     //    H=H(phi_biased)

     FArrayBox& viscfab=(*localMF[CELL_VISC_MATERIAL_MF])[mfi];
     int ncomp_visc=viscfab.nComp();

     int tid_current=ns_thread();
     if ((tid_current<0)||(tid_current>=thread_class::nthreads))
      amrex::Error("tid_current invalid");
     thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     if (fort_built_in_elastic_model(&elastic_viscosity[im],
			           &viscoelastic_model[im])==1) {
       // declared in: GODUNOV_3D.F90
      fort_maketensor(
       &partid,
       &level,
       &finest_level,
       &ncomp_visc,
       &im,  // 0..num_materials-1
       xlo,dx,
       viscfab.dataPtr(),ARLIM(viscfab.loVect()),ARLIM(viscfab.hiVect()),
       tenfab.dataPtr(),ARLIM(tenfab.loVect()),ARLIM(tenfab.hiVect()),
       tilelo,tilehi,
       fablo,fabhi,
       &bfact,
       &elastic_viscosity[im],
       &etaS[im],
       &elastic_time[im],
       &viscoelastic_model[im],
       &polymer_factor[im],
       &NS_geometry_coord);
     } else
      amrex::Error("fort_built_in_elastic_model invalid");
    }  // mfi  
}//omp
    ns_reconcile_d_num(LOOP_MAKETENSOR,"make_viscoelastic_tensor");

    check_for_NAN(localMF[VISCOTEN_MF]);

   } else
    amrex::Error("partid could not be found: make_viscoelastic_tensor");

  } else if (store_elastic_data[im]==0) {

   if (viscoelastic_model[im]!=0)
    amrex::Error("viscoelastic_model[im]!=0");

  } else
   amrex::Error("elastic_time/elastic_viscosity invalid");


 } else if (ns_is_rigid(im)==1) {
  // do nothing
 } else
  amrex::Error("ns_is_rigid invalid");

}  // subroutine make_viscoelastic_tensor


void NavierStokes::make_viscoelastic_heating(int im,int idx) {

 std::string local_caller_string="make_viscoelastic_heating";

 bool use_tiling=ns_tiling;

 int nden=num_materials*num_state_material;

 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=num_materials)) {
  // do nothing
 } else
  amrex::Error("num_materials_viscoelastic bad:make_viscoelastic_heating");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(FACE_VAR_MF+dir,0,local_caller_string);

 debug_ngrow(idx,0,local_caller_string);
 if (localMF[idx]->nComp()!=1)
  amrex::Error("localMF[idx]->nComp() invalid");

 debug_ngrow(CELLTENSOR_MF,1,local_caller_string);
 if (localMF[CELLTENSOR_MF]->nComp()!=AMREX_SPACEDIM_SQR)
  amrex::Error("localMF[CELLTENSOR_MF]->nComp() invalid");

 debug_ngrow(CELL_VISC_MATERIAL_MF,1,local_caller_string);

 debug_ngrow(CELL_DEN_MF,1,local_caller_string); 

 debug_ngrow(CELL_DEDT_MF,1,local_caller_string); 

 if (localMF[CELL_DEN_MF]->nComp()!=1)
  amrex::Error("localMF[CELL_DEN_MF]->nComp() invalid");

 if (localMF[CELL_DEDT_MF]->nComp()!=1)
  amrex::Error("localMF[CELL_DEDT_MF]->nComp() invalid");

 // 1. viscosity coefficient - 1..num_materials
 // 2. viscoelastic coefficient - 1..num_materials
 // 3. relaxation time - 1..num_materials
 // the viscous and viscoelastic forces should both be multiplied by
 // visc_coef.  
 if (localMF[CELL_VISC_MATERIAL_MF]->nComp()!=3*num_materials) {
  std::cout << "ncomp= " <<
   localMF[CELL_VISC_MATERIAL_MF]->nComp() << " num_materials= " << num_materials << '\n';
  amrex::Error("cell_visc_material ncomp invalid(2)");
 }
 if (localMF[CELL_VISC_MATERIAL_MF]->nGrow()<1) {
  std::cout << "ngrow= " <<
   localMF[CELL_VISC_MATERIAL_MF]->nGrow() << " num_materials= " << num_materials << '\n';
  amrex::Error("cell_visc_material ngrow invalid(2)");
 }
 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

 int nstate=STATE_NCOMP;
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");

 if ((im<0)||(im>=num_materials))
  amrex::Error("im invalid53");

 if (ns_is_rigid(im)==0) {

  if (store_elastic_data[im]==1) {

   debug_ngrow(VISCOTEN_MF,1,local_caller_string);
   if (localMF[VISCOTEN_MF]->nComp()!=ENUM_NUM_TENSOR_TYPE)
    amrex::Error("localMF[VISCOTEN_MF] invalid");

   int ncomp_visc=localMF[CELL_VISC_MATERIAL_MF]->nComp();
   if (ncomp_visc!=3*num_materials) {
    std::cout << "ncomp= " <<
     localMF[CELL_VISC_MATERIAL_MF]->nComp() << 
     " num_materials= " << num_materials << '\n';
    amrex::Error("cell_visc_material ncomp invalid (3)");
   }

   resize_levelset(2,LEVELPC_MF);
   debug_ngrow(LEVELPC_MF,2,local_caller_string);
   if (localMF[LEVELPC_MF]->nComp()!=num_materials*(AMREX_SPACEDIM+1))
    amrex::Error("localMF[LEVELPC_MF]->nComp()!=num_materials*(AMREX_SPACEDIM+1)");

   const Real* dx = geom.CellSize();

   if (thread_class::nthreads<1)
    amrex::Error("thread_class::nthreads invalid");
   thread_class::init_d_numPts(localMF[VISCOTEN_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
   for (MFIter mfi(*localMF[VISCOTEN_MF],use_tiling); mfi.isValid(); ++mfi) {
    BL_ASSERT(grids[mfi.index()] == mfi.validbox());
    const int gridno = mfi.index();
    const Box& tilegrid = mfi.tilebox();
    const Box& fabgrid = grids[gridno];
    const int* tilelo=tilegrid.loVect();
    const int* tilehi=tilegrid.hiVect();
    const int* fablo=fabgrid.loVect();
    const int* fabhi=fabgrid.hiVect();
    int bfact=parent->Space_blockingFactor(level);

    const Real* xlo = grid_loc[gridno].lo();

    FArrayBox& tenfab=(*localMF[VISCOTEN_MF])[mfi];
    if (tenfab.nComp()!=ENUM_NUM_TENSOR_TYPE)
     amrex::Error("tenfab.nComp invalid");

    FArrayBox& DeDTinversefab=(*localMF[CELL_DEDT_MF])[mfi]; // 1/(rho cv)
    if (DeDTinversefab.nComp()!=1)
     amrex::Error("DeDTinversefab.nComp() invalid");

    FArrayBox& gradufab=(*localMF[CELLTENSOR_MF])[mfi];
    if (gradufab.nComp()!=AMREX_SPACEDIM_SQR)
     amrex::Error("gradufab.nComp() invalid");

    FArrayBox& heatfab=(*localMF[idx])[mfi];
    if (heatfab.nComp()!=1)
     amrex::Error("heatfab.nComp() invalid");

    FArrayBox& lsfab=(*localMF[LEVELPC_MF])[mfi];
    FArrayBox& viscfab=(*localMF[CELL_VISC_MATERIAL_MF])[mfi];
    if (viscfab.nComp()!=3*num_materials)
     amrex::Error("viscfab.nComp() invalid");

    FArrayBox& xface=(*localMF[FACE_VAR_MF])[mfi];
    FArrayBox& yface=(*localMF[FACE_VAR_MF+1])[mfi];
    FArrayBox& zface=(*localMF[FACE_VAR_MF+AMREX_SPACEDIM-1])[mfi];

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    if (fort_built_in_elastic_model(&elastic_viscosity[im],
		                  &viscoelastic_model[im])==1) {
     // declared in: GODUNOV_3D.F90
     fort_tensorheat(
      &nstate,
      xlo,dx,
      xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()), 
      yface.dataPtr(),ARLIM(yface.loVect()),ARLIM(yface.hiVect()), 
      zface.dataPtr(),ARLIM(zface.loVect()),ARLIM(zface.hiVect()), 
      lsfab.dataPtr(),ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
      DeDTinversefab.dataPtr(),
      ARLIM(DeDTinversefab.loVect()),ARLIM(DeDTinversefab.hiVect()),
      heatfab.dataPtr(),
      ARLIM(heatfab.loVect()),ARLIM(heatfab.hiVect()),
      tenfab.dataPtr(),
      ARLIM(tenfab.loVect()),ARLIM(tenfab.hiVect()),
      gradufab.dataPtr(),
      ARLIM(gradufab.loVect()),ARLIM(gradufab.hiVect()),
      tilelo,tilehi,
      fablo,fabhi,&bfact,&level,
      &dt_slab, //tensorheat
      &NS_geometry_coord,
      &im,&nden);
    } else
     amrex::Error("fort_built_in_elastic_model invalid");
   }  // mfi  
}//omp

   ns_reconcile_d_num(LOOP_TENSORHEAT,"make_viscoelastic_heating");

  } else if (store_elastic_data[im]==0) {

   if (viscoelastic_model[im]!=0)
    amrex::Error("viscoelastic_model[im]!=0");

  } else
   amrex::Error("elastic_time/elastic_viscosity invalid");

 } else if (ns_is_rigid(im)==1) {
  // do nothing
 } else
  amrex::Error("ns_is_rigid invalid");

}   // subroutine make_viscoelastic_heating

void NavierStokes::make_marangoni_force() {

 std::string local_caller_string="make_marangoni_force";

 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();

 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid make_marangoni_force");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 resize_levelset(2,LEVELPC_MF);
 debug_ngrow(LEVELPC_MF,2,local_caller_string);
 if (localMF[LEVELPC_MF]->nComp()!=num_materials*(AMREX_SPACEDIM+1))
  amrex::Error("localMF[LEVELPC_MF]->nComp() invalid");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(FACE_VAR_MF+dir,0,local_caller_string);

 debug_ngrow(DIST_CURV_MF,1,local_caller_string);

 debug_ngrow(CELL_DEN_MF,1,local_caller_string);
 if (localMF[CELL_DEN_MF]->nComp()!=1)
  amrex::Error("localMF[CELL_DEN_MF]->nComp() invalid");

  // mask=1 if not covered or if outside the domain.
  // NavierStokes::maskfiner_localMF
  // NavierStokes::maskfiner
 resize_maskfiner(2,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,2,local_caller_string);
 resize_mask_nbr(2);
 debug_ngrow(MASK_NBR_MF,2,local_caller_string);

 resize_metrics(2);

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

 int nstate=STATE_NCOMP;
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");

 const Real* dx = geom.CellSize();

  // div sigma delta (I-nn^T)/rho=
  //   -sigma kappa grad H/rho + (I-nn^T) (grad sigma) delta/rho 
  // DIST_CURV_MF is calculated in:
  // NavierStokes::makeStateCurv
 int num_curv=num_interfaces*CURVCOMP_NCOMP;
 if (localMF[DIST_CURV_MF]->nComp()!=num_curv)
  amrex::Error("DIST_CURV invalid ncomp");
  
 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(localMF[DIST_CURV_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*localMF[DIST_CURV_MF],use_tiling); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();
  int bfact=parent->Space_blockingFactor(level);
  int bfact_grid=parent->Old_blockingFactor(level);

  const Real* xlo = grid_loc[gridno].lo();

   // DIST_CURV_MF is calculated in:
   // NavierStokes::makeStateCurv
  FArrayBox& curvfab=(*localMF[DIST_CURV_MF])[mfi];
  if (curvfab.nComp()==num_interfaces*CURVCOMP_NCOMP) {
   //do nothing
  } else
   amrex::Error("(curvfab.nComp()!=num_interfaces*CURVCOMP_NCOMP)");

  FArrayBox& rhoinversefab=(*localMF[CELL_DEN_MF])[mfi];
  FArrayBox& lsfab=(*localMF[LEVELPC_MF])[mfi];

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  // force=dt * div ((I-nn^T)(grad sigma) delta) / rho
  // declared in: GODUNOV_3D.F90
  fort_marangoniforce(
   &nstate,
   &num_curv,
   xlo,dx,
   lsfab.dataPtr(),ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
   rhoinversefab.dataPtr(),
   ARLIM(rhoinversefab.loVect()),ARLIM(rhoinversefab.hiVect()),
   curvfab.dataPtr(),ARLIM(curvfab.loVect()),ARLIM(curvfab.hiVect()),
   S_new[mfi].dataPtr(),
   ARLIM(S_new[mfi].loVect()),ARLIM(S_new[mfi].hiVect()),
   tilelo,tilehi,
   fablo,fabhi,
   &bfact,
   &bfact_grid,
   &level,
   &finest_level,
   &dt_slab, //marangoniforce
   &cur_time_slab);
 }  // mfi  
} // omp
 ns_reconcile_d_num(LOOP_MARANGONIFORCE,"make_marangoni_force");

}   // end subroutine make_marangoni_force

void NavierStokes::ns_reconcile_d_num(int caller_loop_id,
		const std::string& caller_string) {

 thread_class::sync_tile_d_numPts();
 ParallelDescriptor::ReduceRealSum(thread_class::tile_d_numPts[0]);
 thread_class::reconcile_d_numPts(caller_loop_id,caller_string);

} // end subroutine ns_reconcile_d_num

int NavierStokes::ns_thread() {
 
 int tid=0;
#ifdef _OPENMP
 tid = omp_get_thread_num();
#endif
 if ((tid<0)||(tid>=thread_class::nthreads))
  amrex::Error("tid invalid");

 return tid;
}  // end subroutine ns_thread()

// add correction term to velocity and/or temperature
// called from:
//  NavierStokes::multiphase_project(before "scale_variablesALL()" is called)
//  NavierStokes::veldiffuseALL
void NavierStokes::make_SEM_delta_force(int project_option) {

 std::string local_caller_string="make_SEM_delta_force";

 bool use_tiling=ns_tiling;

 if ((slab_step<0)||(slab_step>=ns_time_order))
  amrex::Error("slab_step invalid");

 if ((SDC_outer_sweeps<=0)||
     (SDC_outer_sweeps>=ns_time_order))
  amrex::Error("SDC_outer_sweeps invalid");

 if (ns_time_order<=1)
  amrex::Error("ns_time_order invalid");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(FACE_VAR_MF+dir,0,local_caller_string);

 debug_ngrow(delta_MF,0,local_caller_string);
 debug_ngrow(MASKSEM_MF,1,local_caller_string); 

 debug_ngrow(CELL_DEN_MF,1,local_caller_string); 

 debug_ngrow(CELL_DEDT_MF,1,local_caller_string); 

 if (localMF[CELL_DEN_MF]->nComp()!=1)
  amrex::Error("localMF[CELL_DEN_MF]->nComp() invalid");

 if (localMF[CELL_DEDT_MF]->nComp()!=1)
  amrex::Error("localMF[CELL_DEDT_MF]->nComp() invalid");

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

 int nstate=STATE_NCOMP;
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");

 const Real* dx = geom.CellSize();
 int bfact=parent->Space_blockingFactor(level);

 if ((project_option==SOLVETYPE_VISC)|| 
     (project_option==SOLVETYPE_HEAT)) {

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(localMF[delta_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(*localMF[delta_MF],use_tiling); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();

   const Real* xlo = grid_loc[gridno].lo();

   // I-scheme,thermal conduction,viscosity (-div(2 mu D)-force)
   FArrayBox& deltafab=(*localMF[delta_MF])[mfi];
   int deltacomp=0;
   if (project_option==SOLVETYPE_VISC) { 
    deltacomp=slab_step*NSTATE_SDC+SEMDIFFUSE_U;
   } else if (project_option==SOLVETYPE_HEAT) { 
    deltacomp=slab_step*NSTATE_SDC+SEMDIFFUSE_T;
   } else if (project_option==SOLVETYPE_PRES) { 
    amrex::Error("SEM pressure gradient correction on MAC grid");
   } else
    amrex::Error("project_option invalid4");

   FArrayBox& rhoinversefab=(*localMF[CELL_DEN_MF])[mfi];
   FArrayBox& DeDTinversefab=(*localMF[CELL_DEDT_MF])[mfi]; // 1/(rho cv)
   FArrayBox& maskSEMfab=(*localMF[MASKSEM_MF])[mfi];

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // declared in: GODUNOV_3D.F90
   fort_semdeltaforce(
    &nstate,
    &project_option,
    xlo,dx,
    deltafab.dataPtr(deltacomp),
    ARLIM(deltafab.loVect()),ARLIM(deltafab.hiVect()),
    maskSEMfab.dataPtr(),
    ARLIM(maskSEMfab.loVect()),ARLIM(maskSEMfab.hiVect()),
    rhoinversefab.dataPtr(),
    ARLIM(rhoinversefab.loVect()),ARLIM(rhoinversefab.hiVect()),
    DeDTinversefab.dataPtr(),
    ARLIM(DeDTinversefab.loVect()),ARLIM(DeDTinversefab.hiVect()),
    S_new[mfi].dataPtr(),
    ARLIM(S_new[mfi].loVect()),ARLIM(S_new[mfi].hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact,&level,
    &dt_slab); //semdeltaforce
  }  // mfi  
} // omp
  ns_reconcile_d_num(LOOP_SEMDELTAFORCE,"make_SEM_delta_force");

  // pressure gradient at faces.
 } else if (project_option==SOLVETYPE_PRES) {

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

   MultiFab& Umac_new=get_new_data(Umac_Type+dir,slab_step+1);

   if (thread_class::nthreads<1)
    amrex::Error("thread_class::nthreads invalid");
   thread_class::init_d_numPts(localMF[delta_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
   for (MFIter mfi(*localMF[delta_MF],use_tiling); mfi.isValid(); ++mfi) {
    BL_ASSERT(grids[mfi.index()] == mfi.validbox());
    const int gridno = mfi.index();
    const Box& tilegrid = mfi.tilebox();
    const Box& fabgrid = grids[gridno];
    const int* tilelo=tilegrid.loVect();
    const int* tilehi=tilegrid.hiVect();
    const int* fablo=fabgrid.loVect();
    const int* fabhi=fabgrid.hiVect();

    const Real* xlo = grid_loc[gridno].lo();

    FArrayBox& deltafab=(*localMF[delta_GP_MF+dir])[mfi];
    int deltacomp=slab_step;
    FArrayBox& xfacefab=(*localMF[FACE_VAR_MF+dir])[mfi];
    FArrayBox& macfab=Umac_new[mfi];
    FArrayBox& maskSEMfab=(*localMF[MASKSEM_MF])[mfi];

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     // declared in: GODUNOV_3D.F90
    fort_semdeltaforce_face(
     &dir,
     xlo,dx,
     deltafab.dataPtr(deltacomp),
     ARLIM(deltafab.loVect()),ARLIM(deltafab.hiVect()),
     maskSEMfab.dataPtr(),
     ARLIM(maskSEMfab.loVect()),ARLIM(maskSEMfab.hiVect()),
     xfacefab.dataPtr(),
     ARLIM(xfacefab.loVect()),ARLIM(xfacefab.hiVect()),
     macfab.dataPtr(),
     ARLIM(macfab.loVect()),ARLIM(macfab.hiVect()),
     tilelo,tilehi,
     fablo,fabhi,&bfact,&level,
     &dt_slab); //semdeltaforce_face
   }  // mfi  
} // omp
   ns_reconcile_d_num(LOOP_SEMDELTAFORCE_FACE,"make_SEM_delta_force");
  } // dir=0..sdim-1

 } else {
  amrex::Error("project_option invalid in make_SEM_delta_force");
 } 

}   // end subroutine make_SEM_delta_force


// MEHDI VAHAB HEAT SOURCE
// called from veldiffuseALL in NavierStokes3.cpp
void NavierStokes::make_heat_source() {

 std::string local_caller_string="make_heat_source";

 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();

 if ((slab_step<0)||(slab_step>=ns_time_order))
  amrex::Error("slab_step invalid");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(FACE_VAR_MF+dir,0,local_caller_string);

 VOF_Recon_resize(1); //output:SLOPE_RECON_MF
 debug_ngrow(SLOPE_RECON_MF,1,local_caller_string);
 resize_levelset(2,LEVELPC_MF);
 debug_ngrow(LEVELPC_MF,2,local_caller_string);
 if (localMF[LEVELPC_MF]->nComp()!=num_materials*(AMREX_SPACEDIM+1))
  amrex::Error("localMF[LEVELPC_MF]->nComp() invalid");

 resize_metrics(1);  
 debug_ngrow(VOLUME_MF,1,local_caller_string); 

 debug_ngrow(CELL_DEN_MF,1,local_caller_string);
 if (localMF[CELL_DEN_MF]->nComp()!=1)
  amrex::Error("localMF[CELL_DEN_MF]->nComp() invalid");

 debug_ngrow(CELL_DEDT_MF,1,local_caller_string);
 if (localMF[CELL_DEDT_MF]->nComp()!=1)
  amrex::Error("localMF[CELL_DEDT_MF]->nComp() invalid");

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

 int nstate=STATE_NCOMP;
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");

 int nden=num_materials*num_state_material;

 const Real* dx = geom.CellSize();
 int bfact=parent->Space_blockingFactor(level);

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(S_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(S_new,use_tiling); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();

  const Real* xlo = grid_loc[gridno].lo();

  FArrayBox& DeDTinversefab=(*localMF[CELL_DEDT_MF])[mfi]; // 1/(rho cv)
  FArrayBox& lsfab=(*localMF[LEVELPC_MF])[mfi];
  FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi];
  FArrayBox& volfab=(*localMF[VOLUME_MF])[mfi];
  FArrayBox& snewfab=S_new[mfi];

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // MEHDI VAHAB HEAT SOURCE
   // e_{t} + u dot grad e = div k grad T + other terms
   // e = rho cv T (compressible)  rho cp T (incompressible)
   //
   // T^{n+1}=T^{n}+dt * (heat_source)/(rho cv)
   // MKS units:
   // T: Kelvin
   // rho : kg/m^3
   // cv : J/(kg Kelvin)
   // rho cv : J/(m^3 Kelvin)
   // 1/(rho cv) : (m^3 Kelvin)/J
   // heat_source: J/(m^3 s)
   // declared in: GODUNOV_3D.F90
  fort_heatsource(
   &nstate,
   &nden,
   xlo,dx,
   &temperature_source,
   temperature_source_cen.dataPtr(),
   temperature_source_rad.dataPtr(),
   DeDTinversefab.dataPtr(),
   ARLIM(DeDTinversefab.loVect()),ARLIM(DeDTinversefab.hiVect()),
   snewfab.dataPtr(STATECOMP_STATES),
   ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()),
   lsfab.dataPtr(),
   ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
   reconfab.dataPtr(),
   ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
   volfab.dataPtr(),
   ARLIM(volfab.loVect()),ARLIM(volfab.hiVect()),
   tilelo,tilehi,
   fablo,fabhi,&bfact,
   &level,
   &finest_level,
   &dt_slab,//heatsource:time step within a slab SDC, otherwise dt if not SDC.
   &prev_time_slab);
 }  // mfi  
} // omp
 ns_reconcile_d_num(LOOP_HEATSOURCE,"make_heat_source");

}   // end subroutine make_heat_source



void NavierStokes::add_perturbation() {

 std::string local_caller_string="add_perturbation";

 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();

 if (slab_step!=ns_time_order-1)
  amrex::Error("slab_step invalid");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 MultiFab& LS_new = get_new_data(LS_Type,slab_step+1);

 int nstate=STATE_NCOMP;
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");

 const Real* dx = geom.CellSize();
 int bfact=parent->Space_blockingFactor(level);

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(S_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(S_new,use_tiling); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();

  const Real* xlo = grid_loc[gridno].lo();

  FArrayBox& snewfab=S_new[mfi];
  FArrayBox& lsnewfab=LS_new[mfi];

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

   MultiFab& Umac_new=get_new_data(Umac_Type+dir,slab_step+1);
   FArrayBox& macfab=Umac_new[mfi];

   fort_addnoise(
    &dir,
    &angular_velocity,  //parameter for fort_addnoise
    &perturbation_mode, //inputs parameter
    &perturbation_eps_temp, //inputs parameter
    &perturbation_eps_vel,  //inputs parameter
    &nstate,
    xlo,dx,
    snewfab.dataPtr(),
    ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()),
    lsnewfab.dataPtr(),
    ARLIM(lsnewfab.loVect()),ARLIM(lsnewfab.hiVect()),
    macfab.dataPtr(),
    ARLIM(macfab.loVect()),ARLIM(macfab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,
    &bfact,
    &level,
    &finest_level);
  }  // dir
 }  // mfi  
} // omp

 ns_reconcile_d_num(LOOP_ADDNOISE,"add_perturbation");

}   // end subroutine add_perturbation

// called from: update_SEM_forces
// update_SEM_forces is called from: update_SEM_forcesALL
// update_SEM_forcesALL is called from:
//   NavierStokes::do_the_advance
//   NavierStokes::veldiffuseALL
void NavierStokes::update_SEM_delta_force(
 int project_option,
 //GP_DEST_FACE_MF=grad p (MAC grid) SOLVETYPE_PRES
 //MACDIV_MF=-div(k grad T) SOLVETYPE_HEAT
 //MACDIV_MF=-div(2 mu D) SOLVETYPE_VISC
 //HOOP_FORCE_MARK_MF included in this routine (update_SEM_delta_force)
 int update_spectral,
 int update_stable,
 int nsolve) {

 std::string local_caller_string="update_SEM_delta_force";

 bool use_tiling=ns_tiling;

 if ((ns_time_order>=2)&&(ns_time_order<=32)) {
  // do nothing
 } else
  amrex::Error("ns_time_order invalid");

 if (enable_spectral==1) {
  // do nothing
 } else {
  std::cout << "ns_time_order= " << ns_time_order << '\n';
  print_project_option(project_option);
  std::cout << "update_spectral= " << update_spectral << '\n';
  std::cout << "update_stable= " << update_stable << '\n';
  amrex::Error("enable_spectral invalid in update_SEM_delta_force");
 }

 if ((nsolve!=1)&&(nsolve!=AMREX_SPACEDIM))
  amrex::Error("nsolve invalid37");

 if (project_option==SOLVETYPE_PRES) { 
  //do nothing
 } else if (project_option==SOLVETYPE_HEAT) { 
  //do nothing
 } else if (project_option==SOLVETYPE_VISC) { 
  //do nothing
 } else
  amrex::Error("project_option invalid5");

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid update_SEM_forces");

 if ((slab_step<-1)||(slab_step>ns_time_order))
  amrex::Error("slab_step invalid(4)");

 if ((update_spectral==0)||
     (update_spectral==1)) {
  // do nothing
 } else
  amrex::Error("update_spectral invalid update_SEM_delta_force");

 if (update_stable==1) {
  if ((slab_step<0)||(slab_step>=ns_time_order))
   amrex::Error("slab_step invalid");
 } else if (update_stable==0) {
  // do nothing
 } else
  amrex::Error("update_stable invalid update sem delta force");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 debug_ngrow(MACDIV_MF,0,local_caller_string);

 int idx_hoop=MACDIV_MF;

 if (project_option==SOLVETYPE_PRES) { // grad p  (face)
  if (nsolve!=1)
   amrex::Error("nsolve invalid SOLVETYPE_PRES case");
 } else if (project_option==SOLVETYPE_HEAT) { // -div(k grad T)
  if (nsolve!=1)
   amrex::Error("nsolve invalid SOLVETYPE_HEAT case");
  if (localMF[MACDIV_MF]->nComp()!=1) {
   print_project_option(project_option);
   std::cout << "MACDIV_MF = " << MACDIV_MF << '\n';
   std::cout << "nsolve = " << nsolve << '\n';
   std::cout << "localMF[MACDIV_MF] ncomp= " <<
     localMF[MACDIV_MF]->nComp() << '\n';
   amrex::Error("localMF[MACDIV_MF]->nComp() invalid");
  }
 } else if (project_option==SOLVETYPE_VISC) {//-div(2 mu D)-HOOP_FORCE_MARK_MF
  idx_hoop=HOOP_FORCE_MARK_MF;
  if (nsolve!=AMREX_SPACEDIM)
   amrex::Error("nsolve invalid SOLVETYPE_VISC case");
  if (localMF[MACDIV_MF]->nComp()!=nsolve) {
   print_project_option(project_option);
   std::cout << "MACDIV_MF = " << MACDIV_MF << '\n';
   std::cout << "nsolve = " << nsolve << '\n';
   std::cout << "localMF[MACDIV_MF] ncomp= " <<
     localMF[MACDIV_MF]->nComp() << '\n';
   amrex::Error("localMF[MACDIV_MF]->nComp() invalid");
  }
  if (localMF[idx_hoop]->nComp()!=nsolve) {
   print_project_option(project_option);
   std::cout << "idx_hoop = " << idx_hoop << '\n';
   std::cout << "nsolve = " << nsolve << '\n';
   std::cout << "localMF[idx_hoop] ncomp= " <<
     localMF[idx_hoop]->nComp() << '\n';
   amrex::Error("localMF[idx_hoop]->nComp() invalid");
  }
 } else
  amrex::Error("project_option invalid6");

 if (project_option==SOLVETYPE_VISC) { 
  debug_ngrow(idx_hoop,0,local_caller_string);
  if (localMF[idx_hoop]->nComp()!=localMF[MACDIV_MF]->nComp())
   amrex::Error("localMF[idx_hoop]->nComp() invalid");
 } else if (project_option==SOLVETYPE_HEAT) { 
  // check nothing
 } else if (project_option==SOLVETYPE_PRES) {
  // check nothing
 } else
  amrex::Error("project_option invalid 11955");

 if ((project_option==SOLVETYPE_PRES)||  //pressure gradient (face)
     (project_option==SOLVETYPE_HEAT)||  //thermal conductivity
     (project_option==SOLVETYPE_VISC)) { //viscosity

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   debug_ngrow(GP_DEST_FACE_MF+dir,0,local_caller_string);

   if (project_option==SOLVETYPE_PRES) {
    if (nsolve!=1)
     amrex::Error("expecting nsolve==1 if project_option==SOLVETYPE_PRES");
    if (localMF[GP_DEST_FACE_MF+dir]->nComp()!=nsolve)
     amrex::Error("localMF[GP_DEST_FACE_MF+dir]->nComp() invalid");
   } else if ((project_option==SOLVETYPE_HEAT)||
              (project_option==SOLVETYPE_VISC)) {
    // do nothing
   } else
    amrex::Error("project_option invalid7");
 
   debug_ngrow(spectralF_GP_MF+dir,0,local_caller_string);
   debug_ngrow(stableF_GP_MF+dir,0,local_caller_string);
   debug_ngrow(delta_GP_MF+dir,0,local_caller_string);
  } // dir=0..sdim-1

 } else
  amrex::Error("project_option invalid8");

 debug_ngrow(spectralF_MF,0,local_caller_string);
 debug_ngrow(stableF_MF,0,local_caller_string);
 debug_ngrow(delta_MF,0,local_caller_string);

 debug_ngrow(MASKSEM_MF,1,local_caller_string); 

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

 int nstate=STATE_NCOMP;
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");

 const Real* dx = geom.CellSize();
 int bfact=parent->Space_blockingFactor(level);

 if ((project_option==SOLVETYPE_HEAT)||  
     (project_option==SOLVETYPE_VISC)) { 

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(localMF[delta_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(*localMF[delta_MF],use_tiling); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();
 
   const Real* xlo = grid_loc[gridno].lo();
 
   FArrayBox& divfab=(*localMF[MACDIV_MF])[mfi];
   FArrayBox& hoopfab=(*localMF[idx_hoop])[mfi];
   FArrayBox& HOfab=(*localMF[spectralF_MF])[mfi];
   FArrayBox& LOfab=(*localMF[stableF_MF])[mfi];
   FArrayBox& maskSEMfab=(*localMF[MASKSEM_MF])[mfi];

   int LOcomp=0;
   int HOcomp=0;
   if (slab_step==-1)
    HOcomp=0;
   else if ((slab_step>=0)&&(slab_step<ns_time_order))
    HOcomp=NSTATE_SDC*(slab_step+1);
   else
    amrex::Error("slab_step invalid");

   if (update_stable==1) {
    if ((slab_step>=0)&&(slab_step<ns_time_order))
     LOcomp=NSTATE_SDC*slab_step;
    else
     amrex::Error("slab_step invalid");
   } else if (update_stable==0) {
    // do nothing
   } else
    amrex::Error("update_stable invalid");

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();
 
   // declared in: GODUNOV_3D.F90
   fort_updatesemforce(
    &ns_time_order,
    &slab_step,
    &nsolve,
    &update_spectral,
    &update_stable,
    &nstate,
    &project_option,
    xlo,dx,
    divfab.dataPtr(),
    ARLIM(divfab.loVect()),ARLIM(divfab.hiVect()),
    hoopfab.dataPtr(),
    ARLIM(hoopfab.loVect()),ARLIM(hoopfab.hiVect()),
    HOfab.dataPtr(HOcomp),
    ARLIM(HOfab.loVect()),ARLIM(HOfab.hiVect()),
    LOfab.dataPtr(LOcomp),
    ARLIM(LOfab.loVect()),ARLIM(LOfab.hiVect()),
    maskSEMfab.dataPtr(),
    ARLIM(maskSEMfab.loVect()),ARLIM(maskSEMfab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact,&level,
    &dt_slab); //updatesemforce
  }  // mfi  
} // omp
  ns_reconcile_d_num(LOOP_UPDATE_SEMFORCE,"update_SEM_delta_force");

 } else if (project_option==SOLVETYPE_PRES) {

  if (nsolve==1) {
   // do nothing
  } else 
   amrex::Error("nsolve invalid");

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

   if (thread_class::nthreads<1)
    amrex::Error("thread_class::nthreads invalid");
   thread_class::init_d_numPts(localMF[delta_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
   for (MFIter mfi(*localMF[delta_MF],use_tiling); mfi.isValid(); ++mfi) {
    BL_ASSERT(grids[mfi.index()] == mfi.validbox());
    const int gridno = mfi.index();
    const Box& tilegrid = mfi.tilebox();
    const Box& fabgrid = grids[gridno];
    const int* tilelo=tilegrid.loVect();
    const int* tilehi=tilegrid.hiVect();
    const int* fablo=fabgrid.loVect();
    const int* fabhi=fabgrid.hiVect();

    const Real* xlo = grid_loc[gridno].lo();

    FArrayBox& gpfab=(*localMF[GP_DEST_FACE_MF+dir])[mfi];
    FArrayBox& HOfab=(*localMF[spectralF_GP_MF+dir])[mfi];
    FArrayBox& LOfab=(*localMF[stableF_GP_MF+dir])[mfi];
    FArrayBox& maskSEMfab=(*localMF[MASKSEM_MF])[mfi];

    int LOcomp=0;
    int HOcomp=0;
    if (update_spectral==1) {
     if (slab_step==-1)
      HOcomp=0;
     else if ((slab_step>=0)&&(slab_step<ns_time_order))
      HOcomp=(slab_step+1);
     else
      amrex::Error("slab_step invalid");
    } else if (update_spectral==0) {
     // do nothing
    } else
     amrex::Error("update_spectral invalid");

    if (update_stable==1) {
     if ((slab_step>=0)&&(slab_step<ns_time_order))
      LOcomp=slab_step;
     else
      amrex::Error("slab_step invalid");
    } else if (update_stable==0) {
     // do nothing
    } else
     amrex::Error("update_stable invalid");

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     //declared in: GODUNOV_3D.F90
    fort_updatesemforce_face(
     &project_option,
     &ns_time_order,
     &dir,
     &slab_step,
     &update_spectral,
     &update_stable,
     xlo,dx,
     gpfab.dataPtr(),
     ARLIM(gpfab.loVect()),ARLIM(gpfab.hiVect()),
     HOfab.dataPtr(HOcomp),
     ARLIM(HOfab.loVect()),ARLIM(HOfab.hiVect()),
     LOfab.dataPtr(LOcomp),
     ARLIM(LOfab.loVect()),ARLIM(LOfab.hiVect()),
     maskSEMfab.dataPtr(),
     ARLIM(maskSEMfab.loVect()),ARLIM(maskSEMfab.hiVect()),
     tilelo,tilehi,
     fablo,fabhi,&bfact,&level,
     &dt_slab);//updatesemforce_face
   }  // mfi  
} // omp
   ns_reconcile_d_num(LOOP_UPDATESEMFORCE_FACE,"update_SEM_delta_force");

  } // dir=0..sdim-1

 } else {
  amrex::Error("project_option invalid9");
 } 

} // end subroutine update_SEM_delta_force

// called from:
//  NavierStokes::tensor_advection_updateALL()  (NavierStokes3.cpp)
void NavierStokes::tensor_advection_update() {

 std::string local_caller_string="tensor_advection_update";

 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();

 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=num_materials)) {
  // do nothing
 } else
  amrex::Error("num_materials_viscoelastic invalid:tensor_advection_update");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 debug_ngrow(CELLTENSOR_MF,1,local_caller_string);

 debug_ngrow(CELL_DEN_MF,1,local_caller_string);

 debug_ngrow(HOLD_VELOCITY_DATA_MF,1,local_caller_string);
 if (localMF[HOLD_VELOCITY_DATA_MF]->nComp()!=STATE_NCOMP_VEL)
  amrex::Error("localMF[HOLD_VELOCITY_DATA_MF]->nComp()!=STATE_NCOMP_VEL");

 MultiFab& Tensor_new=get_new_data(Tensor_Type,slab_step+1);

 const Real* dx = geom.CellSize();

 int partid_test=0;

 for (int im=0;im<num_materials;im++) {

  if (ns_is_rigid(im)==0) {

   if (store_elastic_data[im]==1) {

    int partid=0;
    while ((im_elastic_map[partid]!=im)&&
	   (partid<im_elastic_map.size())) {
     partid++;
    }

    if (partid==partid_test) {
     //do nothing
    } else
     amrex::Error("partid invalid");

    partid_test++;

    if (partid<im_elastic_map.size()) {

     if (fort_built_in_elastic_model(&elastic_viscosity[im],
      	                             &viscoelastic_model[im])==1) {

      int scomp_tensor=partid*ENUM_NUM_TENSOR_TYPE;

       //CELL_VISC_MATERIAL_MF is build in NavierStokes::getStateVISC_ALL()
       //  1. CELL_VISC_MATERIAL(im)=viscconst(im)  (def) im=1..num_materials
       //  2. CELL_VISC_MATERIAL(num_materials+im)=viscoelastic_coeff *
       //     visc_coef  (def) im=1..num_materials
       //     e.g. viscoelastic_coef=elastic_viscosity/(modtime+dt)
       //  3. CELL_VISC_MATERIAL(2*num_materials+im)=modtime
       //getStateVISC_ALL is called from:
       //  NavierStokes::make_physics_varsALL
       //  NavierStokes::writeTECPLOT_File
       //NavierStokes::init_gradu_tensor_and_material_visc_ALL
      int ncomp_visc=localMF[CELL_VISC_MATERIAL_MF]->nComp();
      if (ncomp_visc!=3*num_materials) {
       std::cout << "ncomp= " <<
        localMF[CELL_VISC_MATERIAL_MF]->nComp() << " num_materials= " << 
	   num_materials << '\n';
       amrex::Error("cell_visc_material ncomp invalid(6)");
      }

      MultiFab* tensor_source_mf=
       getStateTensor(0,scomp_tensor,ENUM_NUM_TENSOR_TYPE,cur_time_slab);

      debug_ngrow(HOLD_GETSHEAR_DATA_MF,0,local_caller_string);

      MultiFab* tendata_mf=localMF[HOLD_GETSHEAR_DATA_MF];
      if (tendata_mf->nGrow()==0) {
       // do nothing
      } else
       amrex::Error("tendata_mf invalid nGrow()");
      if (tendata_mf->nComp()==DERIVE_TENSOR_NCOMP) {
       // do nothing
      } else
       amrex::Error("tendata_mf invalid nComp()");
     
      if (thread_class::nthreads<1)
       amrex::Error("thread_class::nthreads invalid");
      thread_class::init_d_numPts(tensor_source_mf->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
      for (MFIter mfi(*tensor_source_mf,use_tiling); mfi.isValid(); ++mfi) {

       BL_ASSERT(grids[mfi.index()] == mfi.validbox());
       const int gridno = mfi.index();
       const Box& tilegrid = mfi.tilebox();
       const Box& fabgrid = grids[gridno];
       const int* tilelo=tilegrid.loVect();
       const int* tilehi=tilegrid.hiVect();
       const int* fablo=fabgrid.loVect();
       const int* fabhi=fabgrid.hiVect();
       int bfact=parent->Space_blockingFactor(level);

       const Real* xlo = grid_loc[gridno].lo();

       FArrayBox& viscfab=(*localMF[CELL_VISC_MATERIAL_MF])[mfi];
       FArrayBox& one_over_den_fab=(*localMF[CELL_DEN_MF])[mfi];
       FArrayBox& velfab=(*localMF[HOLD_VELOCITY_DATA_MF])[mfi];
       FArrayBox& tensor_new_fab=Tensor_new[mfi];
       FArrayBox& tensor_source_mf_fab=(*tensor_source_mf)[mfi];
       FArrayBox& tendata=(*tendata_mf)[mfi];

       Vector<int> velbc=getBCArray(State_Type,gridno,
        STATECOMP_VEL,STATE_NCOMP_VEL);

       int tid_current=ns_thread();
       if ((tid_current<0)||(tid_current>=thread_class::nthreads))
        amrex::Error("tid_current invalid");
       thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

       if (fort_built_in_elastic_model(&elastic_viscosity[im],
 			            &viscoelastic_model[im])==1) {
        // declared in: GODUNOV_3D.F90
        fort_updatetensor(
         &level,
         &finest_level,
         &im,
         &ncomp_visc,
         viscfab.dataPtr(),ARLIM(viscfab.loVect()),ARLIM(viscfab.hiVect()),
         one_over_den_fab.dataPtr(),
	 ARLIM(one_over_den_fab.loVect()),
	 ARLIM(one_over_den_fab.hiVect()),
         tendata.dataPtr(),ARLIM(tendata.loVect()),ARLIM(tendata.hiVect()),
         dx,xlo,
         velfab.dataPtr(),
         ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
         tensor_new_fab.dataPtr(scomp_tensor),
         ARLIM(tensor_new_fab.loVect()),ARLIM(tensor_new_fab.hiVect()),
         tensor_source_mf_fab.dataPtr(),
         ARLIM(tensor_source_mf_fab.loVect()),
         ARLIM(tensor_source_mf_fab.hiVect()),
         tilelo,tilehi,
         fablo,fabhi,
         &bfact, 
         &dt_slab,//updatetensor
         &elastic_time[im],
         &viscoelastic_model[im],
         &polymer_factor[im],
         &elastic_viscosity[im],
         &NS_geometry_coord,
         velbc.dataPtr(),
         &transposegradu);
       } else {
	std::cout << "illegal: store_elastic_data==1 and visc_model==4\n";
        amrex::Error("fort_built_in_elastic_model invalid");
       }
      }  // mfi
} // omp
      ns_reconcile_d_num(LOOP_UPDATETENSOR,"tensor_advection_update");

      delete tensor_source_mf;
     } else {
      std::cout << "illegal: store_elastic_data==1 and \n";
      std::cout << "fort_built_in_elastic_model!=1";
      amrex::Error("fort_built_in_elastic_model or store_elastic_data invalid");
     }
    } else
     amrex::Error("partid could not be found: tensor_advection_update");

   } else if (store_elastic_data[im]==0) {

    if (viscoelastic_model[im]==0) {
     // do nothing
    } else {
     amrex::Error("viscoelastic_model[im] invalid");
    }

   } else
    amrex::Error("store_elastic_data[im] invalid");

  } else if (ns_is_rigid(im)==1) {

   // do nothing

  } else
   amrex::Error("ns_is_rigid invalid");

 } // im=0..num_materials-1

 if (partid_test==num_materials_viscoelastic) {
  // do nothing
 } else
  amrex::Error("partid_test invalid");

}   // end subroutine tensor_advection_update


// called from:
//  NavierStokes::tensor_advection_updateALL()  (NavierStokes3.cpp)
void NavierStokes::tensor_extrapolation() {

 std::string local_caller_string="tensor_extrapolation";

 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();

 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=num_materials)) {
  // do nothing
 } else
  amrex::Error("num_materials_viscoelastic invalid:tensor_extrapolation");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 MultiFab& Tensor_new=get_new_data(Tensor_Type,slab_step+1);

 const Real* dx = geom.CellSize();

 int partid_test=0;

 for (int im=0;im<num_materials;im++) {

  if (ns_is_rigid(im)==0) {

   if (store_elastic_data[im]==1) {

    int partid=0;
    while ((im_elastic_map[partid]!=im)&&
	   (partid<im_elastic_map.size())) {
     partid++;
    }

    if (partid==partid_test) {
     //do nothing
    } else
     amrex::Error("partid invalid");

    partid_test++;

    if (partid<im_elastic_map.size()) {

     if (fort_built_in_elastic_model(&elastic_viscosity[im],
      	                             &viscoelastic_model[im])==1) {

      int scomp_tensor=partid*ENUM_NUM_TENSOR_TYPE;

      MultiFab* tensor_source_mf=
       getStateTensor(2,scomp_tensor,ENUM_NUM_TENSOR_TYPE,cur_time_slab);

      //LEVELPC_MF is up to date since "allocate_levelset_ALL" was
      //called from "make_physics_varsALL" which was called after
      //the phase change update and before this routine was called.
      resize_levelset(2,LEVELPC_MF);

      if (thread_class::nthreads<1)
       amrex::Error("thread_class::nthreads invalid");
      thread_class::init_d_numPts(tensor_source_mf->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
      for (MFIter mfi(*tensor_source_mf,use_tiling); mfi.isValid(); ++mfi) {

       BL_ASSERT(grids[mfi.index()] == mfi.validbox());
       const int gridno = mfi.index();
       const Box& tilegrid = mfi.tilebox();
       const Box& fabgrid = grids[gridno];
       const int* tilelo=tilegrid.loVect();
       const int* tilehi=tilegrid.hiVect();
       const int* fablo=fabgrid.loVect();
       const int* fabhi=fabgrid.hiVect();
       int bfact=parent->Space_blockingFactor(level);

       const Real* xlo = grid_loc[gridno].lo();

       FArrayBox& LSfab=(*localMF[LEVELPC_MF])[mfi];
       FArrayBox& tensor_new_fab=Tensor_new[mfi];
       FArrayBox& tensor_source_mf_fab=(*tensor_source_mf)[mfi];

       int tid_current=ns_thread();
       if ((tid_current<0)||(tid_current>=thread_class::nthreads))
        amrex::Error("tid_current invalid");
       thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

       if (fort_built_in_elastic_model(&elastic_viscosity[im],
 			            &viscoelastic_model[im])==1) {
        // declared in: GODUNOV_3D.F90
        fort_extrapolate_tensor(
         &level,
         &finest_level,
         &im,
         dx,xlo,
         LSfab.dataPtr(),ARLIM(LSfab.loVect()),ARLIM(LSfab.hiVect()),
         tensor_new_fab.dataPtr(scomp_tensor),
         ARLIM(tensor_new_fab.loVect()),ARLIM(tensor_new_fab.hiVect()),
         tensor_source_mf_fab.dataPtr(),
         ARLIM(tensor_source_mf_fab.loVect()),
         ARLIM(tensor_source_mf_fab.hiVect()),
         tilelo,tilehi,
         fablo,fabhi,
         &bfact);
       } else {
	std::cout << "illegal: store_elastic_data==1 and visc_model==4\n";
        amrex::Error("fort_built_in_elastic_model invalid");
       }
      }  // mfi
} // omp
      ns_reconcile_d_num(LOOP_EXTRAPOLATE_TENSOR,"tensor_extrapolation");

      delete tensor_source_mf;
     } else {
      std::cout << "illegal: store_elastic_data==1 and visc_model==4\n";
      amrex::Error("fort_built_in_elastic_model invalid");
     }
    } else
     amrex::Error("partid could not be found: tensor_extrapolation");

   } else if (store_elastic_data[im]==0) {

    if (viscoelastic_model[im]==0) {
     // do nothing
    } else { 
     amrex::Error("viscoelastic_model[im] invalid");
    }

   } else
    amrex::Error("store_elastic_data[im] invalid");

  } else if (ns_is_rigid(im)==1) {

   // do nothing

  } else
   amrex::Error("ns_is_rigid invalid");

 } // im=0..num_materials-1

 if (partid_test==num_materials_viscoelastic) {
  // do nothing
 } else
  amrex::Error("partid_test invalid");

}   // end subroutine tensor_extrapolation


// non-conservative correction to density.
// if override_density(im)==1,
// rho_im=rho(z)+rho0 * DrhoDT * (T_im - T0_im)
// if override_density(im)=0 or 2, density field is not changed.
// if override_density==2:
// Du/Dt=-grad (p-rho0 g dot z)/rho0 - g DrhoDT (T-T0)
void 
NavierStokes::getStateMOM_DEN(int idx,int ngrow,Real time) {

 std::string local_caller_string="getStateMOM_DEN";

 delete_localMF_if_exist(idx,1);

 bool use_tiling=ns_tiling;

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");
 if (dt_slab>0.0) {
  // do nothing
 } else
  amrex::Error("dt_slab invalid3");

 int finest_level=parent->finestLevel();

 if (ngrow>=1) {
  // do nothing
 } else
  amrex::Error("ngrow>=1 required");

 new_localMF(idx,num_materials,ngrow,-1); // sets values to 0.0
 MultiFab* EOSdata=getStateDen(ngrow,time);

  // NavierStokes::resize_metrics declared in NavierStokes3.cpp
 resize_metrics(ngrow);
  // NavierStokes::resize_maskfiner declared in NavierStokes3.cpp
 resize_maskfiner(ngrow,MASKCOEF_MF);
  // NavierStokes::resize_mask_nbr declared in NavierStokes3.cpp
 resize_mask_nbr(ngrow);

 debug_ngrow(VOLUME_MF,ngrow,local_caller_string); 
 debug_ngrow(MASKCOEF_MF,ngrow,local_caller_string); 
 debug_ngrow(MASK_NBR_MF,ngrow,local_caller_string); 

  // vof,ref centroid,order,slope,intercept  x num_materials
 VOF_Recon_resize(ngrow); //output:SLOPE_RECON_MF
 debug_ngrow(SLOPE_RECON_MF,ngrow,local_caller_string);
 if (localMF[SLOPE_RECON_MF]->nComp()!=num_materials*ngeom_recon)
  amrex::Error("slope_recon_mf has incorrect ncomp");

 const Real* dx = geom.CellSize();

 for (int im=0;im<num_materials;im++) {

  if (DrhoDT[im]==0.0) {
   // check nothing
  } else if (DrhoDT[im]!=0.0) {
   if (override_density[im]==0) {
    amrex::Error("DrhoDT mismatch"); 
   } else if ((override_density[im]==1)||
              (override_density[im]==2)) {
    // do nothing
   } else
    amrex::Error("override_density[im] invalid");
  } else
   amrex::Error("DrhoDT[im] invalid");

  if ((override_density[im]==0)||
      (override_density[im]==1)|| //rho=rho(T,Y,z)
      (override_density[im]==2)) {//Boussinesq approximation
   // do nothing
  } else 
   amrex::Error("override_density[im] invalid");

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(EOSdata->boxArray().d_numPts());
 
#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(*EOSdata,use_tiling); mfi.isValid(); ++mfi) {
    BL_ASSERT(grids[mfi.index()] == mfi.validbox());
    const int gridno = mfi.index();
    const Box& tilegrid = mfi.tilebox();
    const Box& fabgrid = grids[gridno];
    const int* tilelo=tilegrid.loVect();
    const int* tilehi=tilegrid.hiVect();
    const int* fablo=fabgrid.loVect();
    const int* fabhi=fabgrid.hiVect();
    int bfact=parent->Space_blockingFactor(level);

    const Real* xlo = grid_loc[gridno].lo();
    Vector<int> presbc=getBCArray(State_Type,gridno,STATECOMP_PRES,1);
    FArrayBox& maskfab=(*localMF[MASKCOEF_MF])[mfi];
    FArrayBox& masknbrfab=(*localMF[MASK_NBR_MF])[mfi];
    FArrayBox& volfab=(*localMF[VOLUME_MF])[mfi];
    FArrayBox& eosfab=(*EOSdata)[mfi];
    FArrayBox& momfab=(*localMF[idx])[mfi];
    FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi];

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     // declared in: GODUNOV_3D.F90
     // if override_density[im]==1, then rho_im=rho(T) 
    int fort_im=im+1;
    fort_derive_mom_den(
     &fort_im,
     &ngrow,
     constant_density_all_time.dataPtr(),
     spec_material_id_AMBIENT.dataPtr(),
     presbc.dataPtr(),
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,
     &dt_slab,//fort_derive_mom_den
     maskfab.dataPtr(),
     ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
     masknbrfab.dataPtr(),
     ARLIM(masknbrfab.loVect()),ARLIM(masknbrfab.hiVect()),
     volfab.dataPtr(),
     ARLIM(volfab.loVect()),ARLIM(volfab.hiVect()),
     eosfab.dataPtr(),
     ARLIM(eosfab.loVect()),ARLIM(eosfab.hiVect()),
     momfab.dataPtr(),
     ARLIM(momfab.loVect()),ARLIM(momfab.hiVect()),
     reconfab.dataPtr(),
     ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
     xlo,dx,
     DrhoDT.dataPtr(),
     override_density.dataPtr(),
     &level,&finest_level);
  }  // mfi
}  // omp
  ns_reconcile_d_num(LOOP_DERIVE_MOM_DEN,"getStateMOM_DEN");
  
 } // im=0..num_materials-1

 delete EOSdata;

} // end subroutine getStateMOM_DEN

// check that fine grids align with coarse elements and that
// fine grid dimensions are perfectly divisible by the fine order.
// fine grid dimensions are perfectly divisible by the blocking factor.
void NavierStokes::check_grid_places() {

 int finest_level=parent->finestLevel();
 int bfact_SEM_coarse=0;
 int bfact_SEM=parent->Space_blockingFactor(level);
 int bfact_grid=parent->Old_blockingFactor(level);

 int bfact_fine_min=0;
 if ((level>=0)&&(level<finest_level)) {
  bfact_fine_min=((bfact_SEM<4) ? 4 : bfact_SEM);
  if (bfact_grid<4)
   amrex::Error("we must have blocking factor at least 4(1)");
  if (bfact_fine_min<4)
   amrex::Error("bfact_fine_min<4");
 } else if (level==finest_level) {
  bfact_fine_min=((bfact_SEM<4) ? 4 : bfact_SEM);
  if (bfact_grid<4)
   amrex::Error("we must have blocking factor at least 4(1)");
  if (bfact_fine_min<4)
   amrex::Error("bfact_fine_min<4");
 } else
  amrex::Error("level invalid");

 if (bfact_fine_min<bfact_grid)
  bfact_fine_min=bfact_grid;

  // coarse elements must align with fine grid patches.
 if (level>0) {
  bfact_SEM_coarse=parent->Space_blockingFactor(level-1);
  if (bfact_SEM_coarse>bfact_fine_min)
   bfact_fine_min=bfact_SEM_coarse;

  if (2*bfact_SEM_coarse>bfact_fine_min)
   bfact_fine_min=2*bfact_SEM_coarse;
 }

 if ((level>=0)&&(level<finest_level)) {
  if (bfact_fine_min<4)
   amrex::Error("bfact_fine_min<4");
 } else if (level==finest_level) {
  if (bfact_fine_min<4)
   amrex::Error("bfact_fine_min<4");
 } else
  amrex::Error("level invalid");

 for (int gridno=0;gridno<grids.size();gridno++) {
  const Box& fabgrid = grids[gridno];
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   if ((fablo[dir]/bfact_fine_min)*bfact_fine_min!=fablo[dir])
    amrex::Error("fablo index failed bfact_fine_min test");
   if ((fablo[dir]/bfact_grid)*bfact_grid!=fablo[dir])
    amrex::Error("fablo index failed bfact_grid test");
   int testhi=fabhi[dir]+1;
   if ((testhi/bfact_fine_min)*bfact_fine_min!=testhi)
    amrex::Error("fabhi index failed bfact_fine_min test");
   if ((testhi/bfact_grid)*bfact_grid!=testhi)
    amrex::Error("fabhi index failed bfact_grid test");
  } // dir=0..sdim-1

 } // gridno

} // end subroutine check_grid_places


void 
NavierStokes::prepare_mask_nbr(int ngrow) {

 if (ngrow_distance!=4)
  amrex::Error("ngrow_distance invalid");

 if ((ngrow<1)||(ngrow>ngrow_distance))
  amrex::Error("ngrow invalid");

 delete_localMF_if_exist(MASK_NBR_MF,1);

   // mask_nbr:
   // (1) =1 interior  =1 fine-fine ghost in domain  =0 otherwise
   // (2) =1 interior  =0 otherwise
   // (3) =1 interior+ngrow-1  =0 otherwise
   // (4) =1 interior+ngrow    =0 otherwise
 new_localMF(MASK_NBR_MF,4,ngrow,-1);
  // value,scomp,ncomp,ngrow
 localMF[MASK_NBR_MF]->setVal(0.0,0,4,ngrow);
 localMF[MASK_NBR_MF]->setVal(1.0,0,2,0);
   // scomp,ncomp, periodicity
 localMF[MASK_NBR_MF]->FillBoundary(0,1,geom.periodicity());  
 localMF[MASK_NBR_MF]->setVal(1.0,2,1,ngrow-1);
 localMF[MASK_NBR_MF]->setVal(1.0,3,1,ngrow);

} // end subroutine prepare_mask_nbr(int ngrow)

void 
NavierStokes::prepare_displacement() {
 
 bool use_tiling=ns_tiling;

  // 2 ghost cells needed in order to define the displacement MAC
  // velocity for left/right strip of a MAC grid control volume
  // being transported into a target MAC grid control volume.
 int mac_grow=2;

 int finest_level=parent->finestLevel();

 for (int normdir=0;normdir<AMREX_SPACEDIM;normdir++) {

   //Umac_Type
  MultiFab* temp_mac_velocity=getStateMAC(mac_grow,normdir,vel_time_slab); 

   // MAC_VELOCITY_MF deleted towards the end of 
   //   NavierStokes::nonlinear_advection
   // velocity * dt_slab
  new_localMF(MAC_VELOCITY_MF+normdir,1,mac_grow,normdir);

  const Real* dx = geom.CellSize();
  MultiFab& S_new=get_new_data(State_Type,slab_step+1);

  // 1. multiply velocity by dt_slab.  
  // 2. adjust velocity if RZ.  
  // 3. override velocity if it is a passive advection problem.
  // 4. copy into mac_velocity
  // 5. repeat for cell_velocity
  if (thread_class::nthreads<1)
    amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(S_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(S_new,use_tiling); mfi.isValid(); ++mfi) {
    BL_ASSERT(grids[mfi.index()] == mfi.validbox());
    const int gridno = mfi.index();
    const Box& tilegrid = mfi.tilebox();
    const Box& fabgrid = grids[gridno];
    const int* tilelo=tilegrid.loVect();
    const int* tilehi=tilegrid.hiVect();
    const int* fablo=fabgrid.loVect();
    const int* fabhi=fabgrid.hiVect();
    int bfact=parent->Space_blockingFactor(level);

    const Real* xlo = grid_loc[gridno].lo();

    Vector<int> velbc=getBCArray(State_Type,gridno,normdir,1);

    FArrayBox& umactemp=(*temp_mac_velocity)[mfi]; // macgrow
    FArrayBox& umac_displace=
	 (*localMF[MAC_VELOCITY_MF+normdir])[mfi]; // macgrow

    prescribed_vel_time_slab=0.5*(prev_time_slab+cur_time_slab);

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     // declared in: GODUNOV_3D.F90
    fort_velmac_override(
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,
     velbc.dataPtr(),
     &dt_slab,//fort_velmac_override
     &prev_time_slab,
     &prescribed_vel_time_slab,
     &vel_time_slab,
     &dir_absolute_direct_split,
     &normdir,
     umactemp.dataPtr(),
     ARLIM(umactemp.loVect()),ARLIM(umactemp.hiVect()),
     umac_displace.dataPtr(),
     ARLIM(umac_displace.loVect()),ARLIM(umac_displace.hiVect()),
     xlo,dx,
     &mac_grow,
     &map_forward_direct_split[normdir],
     &level,
     &finest_level, 
     &SDC_outer_sweeps,
     &ns_time_order,
     &divu_outer_sweeps,
     &num_divu_outer_sweeps);

  }  // mfi
} // omp

  ns_reconcile_d_num(LOOP_VELMAC_OVERRIDE,"prepare_displacement");

  delete temp_mac_velocity;

  if (cancel_advection==0) {
   // do nothing
  } else if (cancel_advection==1) {
   localMF[MAC_VELOCITY_MF+normdir]->setVal(0.0);
  } else {
   amrex::Error("cancel_advection invalid");
  }

  localMF[MAC_VELOCITY_MF+normdir]->FillBoundary(geom.periodicity());
 } // normdir=0..sdim-1

}  // end subroutine prepare_displacement

//if nucleation_flag==0, then prior to call:
// allocate_levelset_ALL(ngrow_distance,HOLD_LS_DATA_MF);
void
NavierStokes::level_phase_change_rate(Vector<blobclass> blobdata,
	int color_count,int nucleation_flag) {

 std::string local_caller_string="level_phase_change_rate";

 Real problo[AMREX_SPACEDIM];
 Real probhi[AMREX_SPACEDIM];
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  problo[dir]=geom.ProbLo(dir);
  probhi[dir]=geom.ProbHi(dir);
 }

 if (ngrow_distance==4) {
  // do nothing
 } else
  amrex::Error("expecting (ngrow_distance==4)");

 if (ngrow_make_distance==3) {
  // do nothing
 } else
  amrex::Error("expecting (ngrow_make_distance==3)");

 bool use_tiling=ns_tiling;
 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid level_phase_change_rate");

 int nden=num_materials*num_state_material;
 int nburning=EXTRAP_NCOMP_BURNING;
 int ntsat=EXTRAP_NCOMP_TSAT;
 int nstate=STATE_NCOMP;

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int n_normal=(num_materials+num_interfaces)*(AMREX_SPACEDIM+1);

 const Real* dx = geom.CellSize();

 Real dxmax=0.0;
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  if (dx[dir]>0.0) {
   if (dxmax<dx[dir])
    dxmax=dx[dir];
  } else
   amrex::Error("dx[dir] must be positive");
 } // dir=0..sdim-1

  // in: level_phase_change_rate()
 getStateDen_localMF(DEN_RECON_MF,ngrow_distance,cur_time_slab);

 if (localMF[DEN_RECON_MF]->nComp()!=nden)
  amrex::Error("DEN_RECON_MF invalid ncomp");
 debug_ixType(DEN_RECON_MF,-1,local_caller_string);

 Vector<Real> blob_array;
 int blob_arraysize=num_elements_blobclass;
 blob_array.resize(blob_arraysize);

 if (nucleation_flag==0) {

  if (color_count!=blobdata.size())
   amrex::Error("color_count!=blobdata.size()");
  blob_arraysize=color_count*num_elements_blobclass;
  blob_array.resize(blob_arraysize);

  int counter=0;
  for (int i=0;i<color_count;i++) {
   copy_from_blobdata(i,counter,blob_array,blobdata);
  } // i=0..color_count-1

  if (counter!=blob_arraysize)
   amrex::Error("counter invalid");

  if (localMF[COLOR_MF]->nGrow()!=1)
   amrex::Error("localMF[COLOR_MF]->nGrow()!=1");
  if (localMF[TYPE_MF]->nGrow()!=1)
   amrex::Error("localMF[TYPE_MF]->nGrow()!=1");
  debug_ixType(COLOR_MF,-1,local_caller_string);
  debug_ixType(TYPE_MF,-1,local_caller_string);

 } else if (nucleation_flag==1) {
  if (color_count!=1)
   amrex::Error("color_count!=1");
 } else
  amrex::Error("nucleation_flag invalid");

 MultiFab* presmf=getState(ngrow_distance,STATECOMP_PRES,
           1,cur_time_slab);
 debug_ixType_raw(presmf,-1,local_caller_string);


 MultiFab* pres_eos_mf=derive_EOS_pressure(material_type_evap);
 if (pres_eos_mf->nGrow()!=1)
  amrex::Error("pres_eos_mf->nGrow()!=1");
 debug_ixType_raw(pres_eos_mf,-1,local_caller_string);

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,1,local_caller_string); 
 debug_ixType(MASKCOEF_MF,-1,local_caller_string);

 if ((prev_time_slab<0.0)||
     (cur_time_slab<=0.0)||
     (cur_time_slab<=prev_time_slab))
  amrex::Error("prev_time_slab or cur_time_slab invalid");

 int do_the_nucleate;
 Vector<Real> nucleate_pos;
 nucleate_pos.resize(4);
 if (n_sites>0)
  nucleate_pos.resize(4*n_sites);

 do_the_nucleate=0;
 for (int i=0;i<nucleate_pos.size();i++) 
  nucleate_pos[i]=0.0;

 if (n_sites==0) {
  // do nothing
 } else if (n_sites>0) {

  int first_time_nucleate=0;
  if (nucleation_init_time==0.0) {
   if (prev_time_slab==0.0) 
    first_time_nucleate=1;
  } else if (nucleation_init_time>0.0) {
   if ((prev_time_slab<=nucleation_init_time)&&
       (cur_time_slab>nucleation_init_time)) 
    first_time_nucleate=1;
  } else
   amrex::Error("nucleation_init_time invalid");
 
  if (nucleation_period==0.0) { // just nucleate the bubble(s) once
   if (first_time_nucleate==1) {
    do_the_nucleate=1;
    int pos_comp=0;
    for (int i_site=0;i_site<n_sites;i_site++) {
     for (int dir=0;dir<3;dir++) {
	 // MITSUHIRO BUBBLE POSITIONS and RADII
      nucleate_pos[pos_comp]=pos_sites[pos_comp]; 
      pos_comp++;
     }
     double rr=pos_sites[pos_comp];
     if (rr>0.0) {
      if (rr<2.0*dxmax)
       rr=2.0*dxmax;
     } else
      amrex::Error("rr must be positive");
	 // MITSUHIRO BUBBLE POSITIONS and RADII
     nucleate_pos[pos_comp]=rr;
     pos_comp++;
    } // i_site=0..n_sites-1
   } else if (first_time_nucleate==0) {
    // do nothing
   } else
    amrex::Error("first_time_nucleate invalid");

  } else if (nucleation_period>0.0) { // positions random t>t_nucleate_init

   if (first_time_nucleate==1) {

    do_the_nucleate=1;
    int pos_comp=0;
    for (int i_site=0;i_site<n_sites;i_site++) {
     for (int dir=0;dir<3;dir++) {
	 // MITSUHIRO BUBBLE POSITIONS and RADII
      nucleate_pos[pos_comp]=pos_sites[pos_comp]; 
      pos_comp++;
     }
     double rr=pos_sites[pos_comp];
     if (rr>0.0) {
      if (rr<2.0*dxmax)
       rr=2.0*dxmax;
     } else
      amrex::Error("rr must be positive");
	 // MITSUHIRO BUBBLE POSITIONS and RADII
     nucleate_pos[pos_comp]=rr;
     pos_comp++;
    } // i_site=0..n_sites-1

   } else if ((first_time_nucleate==0)&&
              (prev_time_slab>nucleation_init_time)) {
  
    if (level==finest_level) {
 
     int num_periods=0;
     Real mult_period=nucleation_init_time;

     while (mult_period<prev_time_slab) {
      num_periods++;
      mult_period=nucleation_init_time+ 
        num_periods*nucleation_period;
     }

     if (1==0) {
      std::cout << "num_periods= " << num_periods << '\n';
      std::cout << "nucleation_period= " << nucleation_period << '\n';
      std::cout << "prev_time_slab= " << prev_time_slab << '\n';
      std::cout << "cur_time_slab= " << cur_time_slab << '\n';
      std::cout << "mult_period= " << mult_period << '\n';
      std::cout << "nucleation_init_time= " << 
       nucleation_init_time << '\n';
     }

     if (mult_period<cur_time_slab) {

      do_the_nucleate=1;

      for (int nc=0;nc<n_sites;nc++) {

       double rr=pos_sites[nc*4+3];  // radius
       if (rr>0.0) {
        if (rr<2.0*dxmax)
         rr=2.0*dxmax;
       } else
        amrex::Error("rr must be positive");

       Vector<Real> xnucleate(AMREX_SPACEDIM);
       for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
        xnucleate[dir]=-1.0e+99;

       if (ParallelDescriptor::IOProcessor()) {
        for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
         Real save_random=amrex::Random();
         if ((save_random<0.0)||(save_random>1.0)) {
          std::cout << "save_random invalid save_random= " << 
           save_random << '\n';
          amrex::Error("save_random bust");
         }
         
         xnucleate[dir]=problo[dir]+save_random*(probhi[dir]-problo[dir]);
          // (int nc=0;nc<n_sites;nc++)
	 if (pos_sites_random_flag==0) {

          xnucleate[dir]=pos_sites[nc*4+dir];

	 } else if (pos_sites_random_flag==1) {

	  if ((NS_geometry_coord==COORDSYS_RZ)&&(dir==0)) {
           xnucleate[dir]=0.0;
          } else if (dir==AMREX_SPACEDIM-1) {  // vertical coordinate
           xnucleate[dir]=pos_sites[nc*4+dir];
          } else {       
           if (xnucleate[dir]-rr<problo[dir])
            xnucleate[dir]=problo[dir]+rr;
           if (xnucleate[dir]+rr>probhi[dir])
            xnucleate[dir]=probhi[dir]-rr;
          }
	 } else
          amrex::Error("pos_sites_random_flag invalid");
        } // dir
       } // io proc?

       for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
        ParallelDescriptor::ReduceRealMax(xnucleate[dir]);

	 // MITSUHIRO BUBBLE POSITIONS and RADII
	 // (THIS IS THE RANDOM CASE)
       for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
        nucleate_pos[nc*4+dir]=xnucleate[dir];
       nucleate_pos[nc*4+3]=rr;  // radius

      } // nc=0..n_sites-1

     } // mult_period<cur_time_slab

    } else if ((level>=0)&&(level<finest_level)) {
     // do nothing
    } else
     amrex::Error("level invalid nucleate_bubbles 2");

   } else if ((first_time_nucleate==0)&&
              (prev_time_slab<=nucleation_init_time)) {
    // do nothing
   } else
    amrex::Error("first_time_nucleate or prev_time_slab invalid");
  } else 
   amrex::Error("nucleation_period invalid");
 } else {
  std::cout << "n_sites= " << n_sites << '\n';
  amrex::Error("n_sites invalid(3)");
 }

 int nucleate_pos_size=nucleate_pos.size();

 MultiFab& LS_new = get_new_data(LS_Type,slab_step+1);
 if (LS_new.nComp()!=num_materials*(1+AMREX_SPACEDIM)) 
  amrex::Error("LS_new invalid ncomp");

 MultiFab& S_new = get_new_data(State_Type,slab_step+1);
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");

 if (ngrow_distance==ngrow_make_distance+1) {
  // do nothing
 } else
  amrex::Error("ngrow_distance!=ngrow_make_distance+1");

 if (nucleation_flag==0) {

  if (localMF[BURNING_VELOCITY_MF]->nComp()!=nburning)
   amrex::Error("localMF[BURNING_VELOCITY_MF] incorrect ncomp");
  if (localMF[BURNING_VELOCITY_MF]->nGrow()!=ngrow_distance)
   amrex::Error("localMF[BURNING_VELOCITY_MF] incorrect ngrow");
  debug_ixType(BURNING_VELOCITY_MF,-1,local_caller_string);

  if (localMF[SATURATION_TEMP_MF]->nComp()!=ntsat)
   amrex::Error("localMF[SATURATION_TEMP_MF]->nComp()!=ntsat");
  if (localMF[SATURATION_TEMP_MF]->nGrow()!=ngrow_distance)
   amrex::Error("localMF[SATURATION_TEMP_MF] incorrect ngrow");
  debug_ixType(SATURATION_TEMP_MF,-1,local_caller_string);

  if (localMF[FD_NRM_ND_MF]->nComp()!=n_normal)
   amrex::Error("localMF[FD_NRM_ND_MF]->nComp()!=n_normal");
  if (localMF[FD_NRM_ND_MF]->nGrow()!=ngrow_distance)
   amrex::Error("localMF[FD_NRM_ND_MF] incorrect ngrow");
  debug_ixType(FD_NRM_ND_MF,-1,local_caller_string);
  
  if (localMF[FD_CURV_CELL_MF]->nComp()!=2*(num_materials+num_interfaces))
   amrex::Error(
    "localMF[FD_CURV_CELL_MF]->nComp()!=2*(num_materials+num_interfaces)");
  if (localMF[FD_CURV_CELL_MF]->nGrow()!=ngrow_make_distance)
   amrex::Error("localMF[FD_CURV_CELL_MF] incorrect ngrow");
  debug_ixType(FD_CURV_CELL_MF,-1,local_caller_string);

  debug_ngrow(HOLD_LS_DATA_MF,ngrow_distance,local_caller_string);
  if (localMF[HOLD_LS_DATA_MF]->nComp()!=num_materials*(1+AMREX_SPACEDIM)) 
   amrex::Error("localMF[HOLD_LS_DATA_MF]->nComp() invalid");
  debug_ixType(HOLD_LS_DATA_MF,-1,local_caller_string);

  debug_ngrow(MDOT_MF,0,local_caller_string);
  debug_ixType(MDOT_MF,-1,local_caller_string);

  VOF_Recon_resize(ngrow_distance); //output:SLOPE_RECON_MF
  debug_ixType(SLOPE_RECON_MF,-1,local_caller_string);

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(LS_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(LS_new,use_tiling); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();
   int bfact=parent->Space_blockingFactor(level);

   const Real* xlo = grid_loc[gridno].lo();

   FArrayBox& lsfab=(*localMF[HOLD_LS_DATA_MF])[mfi];
   FArrayBox& FD_NRM_ND_fab=(*localMF[FD_NRM_ND_MF])[mfi];

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // fort_fd_node_normal is declared in: MOF_REDIST_3D.F90
    // The output from this routine is used for finding the curvature
    // which is used in ratemasschange.
   fort_fd_node_normal( 
    &level,
    &finest_level,
    lsfab.dataPtr(),
    ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
    FD_NRM_ND_fab.dataPtr(),
    ARLIM(FD_NRM_ND_fab.loVect()),
    ARLIM(FD_NRM_ND_fab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    xlo,dx,
    &n_normal,
    &ngrow_make_distance);
  } // mfi
} // omp
  ns_reconcile_d_num(LOOP_FD_NODE_NORMAL,"level_phase_change_rate");

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(LS_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(LS_new,use_tiling); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();
   int bfact=parent->Space_blockingFactor(level);

   const Real* xlo = grid_loc[gridno].lo();

   FArrayBox& F_fab=(*localMF[SLOPE_RECON_MF])[mfi];
   FArrayBox& lsfab=(*localMF[HOLD_LS_DATA_MF])[mfi];
   FArrayBox& FD_NRM_ND_fab=(*localMF[FD_NRM_ND_MF])[mfi];
   FArrayBox& curvfab=(*localMF[FD_CURV_CELL_MF])[mfi];

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // the normals fab has ngrow_make_distance+1 ghost cells
    //    growntileboxNODE(ngrow_make_distance)
    // the curvature fab should have ngrow_make_distance ghost cells.
    // fort_node_to_cell is declared in: MOF_REDIST_3D.F90
   int height_function_flag=1;
   fort_node_to_cell( 
    &tid_current,
    &level,
    &finest_level,
    &height_function_flag,
    F_fab.dataPtr(),
    ARLIM(F_fab.loVect()),ARLIM(F_fab.hiVect()),
    lsfab.dataPtr(),
    ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
    FD_NRM_ND_fab.dataPtr(),
    ARLIM(FD_NRM_ND_fab.loVect()),
    ARLIM(FD_NRM_ND_fab.hiVect()),
    curvfab.dataPtr(),
    ARLIM(curvfab.loVect()),ARLIM(curvfab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,
    &bfact,
    xlo,dx,
    &n_normal,
    &ngrow_make_distance);
  } // mfi
} // omp
  ns_reconcile_d_num(LOOP_NODE_TO_CELL,"level_phase_change_rate");

  localMF[FD_CURV_CELL_MF]->FillBoundary(geom.periodicity());

 } else if (nucleation_flag==1) {
  // do nothing
 } else
  amrex::Error("nucleation_flag invalid");

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(LS_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(LS_new,use_tiling); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();
   int bfact=parent->Space_blockingFactor(level);

   const Real* xlo = grid_loc[gridno].lo();

    // mask=tag if not covered by level+1 or outside the domain.
   FArrayBox& maskcov=(*localMF[MASKCOEF_MF])[mfi];
   FArrayBox& lsnewfab=LS_new[mfi];
   FArrayBox& snewfab=S_new[mfi];
   FArrayBox& eosfab=(*localMF[DEN_RECON_MF])[mfi];

   FArrayBox& presfab=(*presmf)[mfi]; 
   FArrayBox& pres_eos_fab=(*pres_eos_mf)[mfi]; 

   FArrayBox& thermal_conductivity_fab=
	  (*localMF[CELL_CONDUCTIVITY_MATERIAL_MF])[mfi];
   if (thermal_conductivity_fab.nComp()!=num_materials)
    amrex::Error("thermal_conductivity_fab.nComp()!=num_materials");

   Vector<int> use_exact_temperature(2*num_interfaces);
   for (int im=0;im<2*num_interfaces;im++)
    use_exact_temperature[im]=0;

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   if (nucleation_flag==0) {

    FArrayBox& lsfab=(*localMF[HOLD_LS_DATA_MF])[mfi];

    FArrayBox& burnvelfab=(*localMF[BURNING_VELOCITY_MF])[mfi];
    if (burnvelfab.nComp()!=nburning)
     amrex::Error("burnvelfab.nComp() incorrect");

     // ntsat=num_interfaces*(ncomp_per_tsat+1)
     // e.g. for interface 12,
     //  component 1=0 if T_gamma,Y_gamma not defined
     //             =1 if T_gamma,Y_gamma defined in fort_ratemasschange
     //             =2 if T_gamma,Y_gamma defined after extrapolation
     //             =-1 or -2 for condensation case.
     //  component 2=T_gamma
     //  component 3=Y_gamma
     //  repeats ....
    FArrayBox& Tsatfab=(*localMF[SATURATION_TEMP_MF])[mfi];
    if (Tsatfab.nComp()!=ntsat) {
     std::cout << "Tsatfab.nComp()=" << Tsatfab.nComp() << 
       " ntsat=" << ntsat << '\n';
     amrex::Error("Tsatfab.nComp()!=ntsat 2");
    }

    FArrayBox& colorfab=(*localMF[COLOR_MF])[mfi];
    FArrayBox& typefab=(*localMF[TYPE_MF])[mfi];
     // fluids tessellate; solids overlap
    FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi]; 

    FArrayBox& curvfab=(*localMF[FD_CURV_CELL_MF])[mfi];

     // lsnewfab and burnvelfab are updated.
     // lsfab is not updated.
     // burnvelfab=BURNING_VELOCITY_MF is cell centered velocity.
    fort_ratemasschange( 
     &tid_current,
     &nucleation_flag,
     &level,
     &finest_level,
     &ngrow_distance,
     &nstate,
     &nburning,
     &ntsat,
     &nden,
     &custom_nucleation_model,
     &do_the_nucleate,
     nucleate_pos.dataPtr(),
     &nucleate_pos_size, 
     nucleation_temp.dataPtr(), 
     nucleation_pressure.dataPtr(), 
     nucleation_pmg.dataPtr(), 
     nucleation_mach.dataPtr(), 
     cavitation_pressure.dataPtr(), 
     cavitation_vapor_density.dataPtr(), 
     cavitation_tension.dataPtr(), 
     microlayer_substrate.dataPtr(),
     microlayer_angle.dataPtr(),
     microlayer_size.dataPtr(),
     macrolayer_size.dataPtr(),
     max_contact_line_size.dataPtr(),
     &R_Palmore_Desjardins,
     use_exact_temperature.dataPtr(),
     reaction_rate.dataPtr(),
     hardwire_Y_gamma.dataPtr(),
     hardwire_T_gamma.dataPtr(),
     accommodation_coefficient.dataPtr(),
     reference_pressure.dataPtr(),
     saturation_temp.dataPtr(),
     saturation_temp_curv.dataPtr(),
     saturation_temp_vel.dataPtr(),
     saturation_temp_min.dataPtr(),
     saturation_temp_max.dataPtr(),
     freezing_model.dataPtr(),
     Tanasawa_or_Schrage_or_Kassemi.dataPtr(),
     interface_mass_transfer_model.dataPtr(),
     distribute_from_target.dataPtr(),
     mass_fraction_id.dataPtr(),
     constant_density_all_time.dataPtr(),
     material_type_evap.dataPtr(),
     molar_mass.dataPtr(),
     species_molar_mass.dataPtr(),
     tilelo,tilehi,
     fablo,fabhi,&bfact,
     xlo,dx,
     &prev_time_slab,
     &dt_slab,//ratemasschange
     &blob_arraysize,
     blob_array.dataPtr(),
     &color_count,
     colorfab.dataPtr(),
     ARLIM(colorfab.loVect()),ARLIM(colorfab.hiVect()),
     typefab.dataPtr(),
     ARLIM(typefab.loVect()),ARLIM(typefab.hiVect()),
     maskcov.dataPtr(),
     ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
     thermal_conductivity_fab.dataPtr(), //num_materials components
     ARLIM(thermal_conductivity_fab.loVect()),
     ARLIM(thermal_conductivity_fab.hiVect()),
     burnvelfab.dataPtr(),
     ARLIM(burnvelfab.loVect()),ARLIM(burnvelfab.hiVect()),
     Tsatfab.dataPtr(),
     ARLIM(Tsatfab.loVect()),ARLIM(Tsatfab.hiVect()),
     lsfab.dataPtr(),ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
     lsnewfab.dataPtr(),ARLIM(lsnewfab.loVect()),ARLIM(lsnewfab.hiVect()),
     snewfab.dataPtr(),ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()),
     eosfab.dataPtr(),ARLIM(eosfab.loVect()),ARLIM(eosfab.hiVect()),
     reconfab.dataPtr(),ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
     presfab.dataPtr(),ARLIM(presfab.loVect()),ARLIM(presfab.hiVect()),
     pres_eos_fab.dataPtr(),
     ARLIM(pres_eos_fab.loVect()),ARLIM(pres_eos_fab.hiVect()),
     curvfab.dataPtr(),
     ARLIM(curvfab.loVect()),
     ARLIM(curvfab.hiVect()));

   } else if (nucleation_flag==1) {

    fort_ratemasschange( 
     &tid_current,
     &nucleation_flag,
     &level,
     &finest_level,
     &ngrow_distance,
     &nstate,
     &nburning,
     &ntsat,
     &nden,
     &custom_nucleation_model,
     &do_the_nucleate,
     nucleate_pos.dataPtr(),
     &nucleate_pos_size, 
     nucleation_temp.dataPtr(), 
     nucleation_pressure.dataPtr(), 
     nucleation_pmg.dataPtr(), 
     nucleation_mach.dataPtr(), 
     cavitation_pressure.dataPtr(), 
     cavitation_vapor_density.dataPtr(), 
     cavitation_tension.dataPtr(), 
     microlayer_substrate.dataPtr(),
     microlayer_angle.dataPtr(),
     microlayer_size.dataPtr(),
     macrolayer_size.dataPtr(),
     max_contact_line_size.dataPtr(),
     &R_Palmore_Desjardins,
     use_exact_temperature.dataPtr(),
     reaction_rate.dataPtr(),
     hardwire_Y_gamma.dataPtr(),
     hardwire_T_gamma.dataPtr(),
     accommodation_coefficient.dataPtr(),
     reference_pressure.dataPtr(),
     saturation_temp.dataPtr(),
     saturation_temp_curv.dataPtr(),
     saturation_temp_vel.dataPtr(),
     saturation_temp_min.dataPtr(),
     saturation_temp_max.dataPtr(),
     freezing_model.dataPtr(),
     Tanasawa_or_Schrage_or_Kassemi.dataPtr(),
     interface_mass_transfer_model.dataPtr(),
     distribute_from_target.dataPtr(),
     mass_fraction_id.dataPtr(),
     constant_density_all_time.dataPtr(),
     material_type_evap.dataPtr(),
     molar_mass.dataPtr(),
     species_molar_mass.dataPtr(),
     tilelo,tilehi,
     fablo,fabhi,&bfact,
     xlo,dx,
     &prev_time_slab,
     &dt_slab,//ratemasschange
     &blob_arraysize,
     blob_array.dataPtr(),
     &color_count,
     lsnewfab.dataPtr(), //colorfab
     ARLIM(lsnewfab.loVect()),ARLIM(lsnewfab.hiVect()),
     lsnewfab.dataPtr(), //typefab
     ARLIM(lsnewfab.loVect()),ARLIM(lsnewfab.hiVect()),
     maskcov.dataPtr(),
     ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
     thermal_conductivity_fab.dataPtr(), //num_materials components
     ARLIM(thermal_conductivity_fab.loVect()),
     ARLIM(thermal_conductivity_fab.hiVect()),
     lsnewfab.dataPtr(), //burnvelfab
     ARLIM(lsnewfab.loVect()),ARLIM(lsnewfab.hiVect()),
     lsnewfab.dataPtr(), //Tsatfab
     ARLIM(lsnewfab.loVect()),ARLIM(lsnewfab.hiVect()),
     lsnewfab.dataPtr(), //lsfab
     ARLIM(lsnewfab.loVect()),ARLIM(lsnewfab.hiVect()),
     lsnewfab.dataPtr(),
     ARLIM(lsnewfab.loVect()),ARLIM(lsnewfab.hiVect()),
     snewfab.dataPtr(),
     ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()),
     eosfab.dataPtr(),ARLIM(eosfab.loVect()),ARLIM(eosfab.hiVect()),
     lsnewfab.dataPtr(), //reconfab
     ARLIM(lsnewfab.loVect()),ARLIM(lsnewfab.hiVect()),
     presfab.dataPtr(),
     ARLIM(presfab.loVect()),ARLIM(presfab.hiVect()),
     pres_eos_fab.dataPtr(),
     ARLIM(pres_eos_fab.loVect()),ARLIM(pres_eos_fab.hiVect()),
     lsnewfab.dataPtr(), //curvfab
     ARLIM(lsnewfab.loVect()),ARLIM(lsnewfab.hiVect()));
   } else
    amrex::Error("nucleation_flag invalid");

 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_RATEMASSCHANGE,"level_phase_change_rate");

 delete_localMF(DEN_RECON_MF,1);

 delete presmf;
 delete pres_eos_mf;

} // end subroutine level_phase_change_rate



void
NavierStokes::level_phase_change_rate_extend() {

 std::string local_caller_string="level_phase_change_rate_extend";

 bool use_tiling=ns_tiling;
 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid level_phase_change_rate_extend");


 int ncomp_per_burning=EXTRAP_PER_BURNING;
 int ncomp_per_tsat=EXTRAP_PER_TSAT;

  // flag 1 .. num_interfaces, vel 1 .. num_interfaces
 int nburning=EXTRAP_NCOMP_BURNING;
 int ntsat=EXTRAP_NCOMP_TSAT;

 const Real* dx = geom.CellSize();

 MultiFab& LS_new = get_new_data(LS_Type,slab_step+1);
 if (LS_new.nComp()!=num_materials*(AMREX_SPACEDIM+1)) 
  amrex::Error("LS_new invalid ncomp");

 if (localMF[BURNING_VELOCITY_MF]->nComp()!=nburning)
  amrex::Error("localMF[BURNING_VELOCITY_MF] incorrect ncomp");
 if (localMF[SATURATION_TEMP_MF]->nComp()!=ntsat)
  amrex::Error("localMF[SATURATION_TEMP_MF] incorrect ncomp");

 if (ngrow_make_distance!=3)
  amrex::Error("expecting ngrow_make_distance==3");
 if (ngrow_distance!=4)
  amrex::Error("expecting ngrow_distance==4");

 if (localMF[BURNING_VELOCITY_MF]->nGrow()!=ngrow_distance)
  amrex::Error("localMF[BURNING_VELOCITY_MF] incorrect ngrow");
 if (localMF[SATURATION_TEMP_MF]->nGrow()!=ngrow_distance)
  amrex::Error("localMF[SATURATION_TEMP_MF] incorrect ngrow");

 debug_ngrow(HOLD_LS_DATA_MF,ngrow_distance,local_caller_string);
 if (localMF[HOLD_LS_DATA_MF]->nComp()!=num_materials*(1+AMREX_SPACEDIM)) 
  amrex::Error("localMF[HOLD_LS_DATA_MF]->nComp() invalid");

 for (int velflag=0;velflag<=1;velflag++) {

  int ncomp=0;
  int ncomp_per_interface=0;
  if (velflag==0) {
   ncomp_per_interface=ncomp_per_tsat; // interface temperature, mass fraction
  } else if (velflag==1) {
   ncomp_per_interface=ncomp_per_burning;
  } else
   amrex::Error("velflag invalid");

  ncomp=num_interfaces+num_interfaces*ncomp_per_interface;

  Vector<int> scompBC_map;
  scompBC_map.resize(ncomp);
   // extrap, u_extrap, v_extrap, w_extrap
   // mof recon extrap
   // maskSEMextrap
  int burnvel_start_pos_base=EXTRAPCOMP_BURNVEL;
  int extend_start_pos=burnvel_start_pos_base;
  if (velflag==0) {
   extend_start_pos=EXTRAPCOMP_TSAT;
  } else if (velflag==1) { 
   extend_start_pos=burnvel_start_pos_base;
  } else
   amrex::Error("velflag invalid");

  for (int imdest=0;imdest<ncomp;imdest++)
   scompBC_map[imdest]=extend_start_pos+imdest;

  int local_mf=0;
  if (velflag==0) {
   local_mf=SATURATION_TEMP_MF;
  } else if (velflag==1) {
   local_mf=BURNING_VELOCITY_MF;
  } else
   amrex::Error("velflag invalid");

   //scomp=0
  PCINTERP_fill_borders(local_mf,ngrow_distance,
   0,ncomp,State_Type,scompBC_map);

  if (1==0) {
   int gridno=0;
   const Box& fabgrid = grids[gridno];
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();
   const Real* xlo = grid_loc[gridno].lo();
   int interior_only=0;
   FArrayBox& burnvelfab=(*localMF[local_mf])[0];
   const Real* dxplot = geom.CellSize();
   int scomp=0;
   int dirplot=-1;
   int id=0;
   tecplot_debug(burnvelfab,xlo,fablo,fabhi,dxplot,dirplot,id,
     scomp,ncomp,interior_only);
  }

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(LS_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(LS_new,use_tiling); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();
   int bfact=parent->Space_blockingFactor(level);

   const Real* xlo = grid_loc[gridno].lo();
   FArrayBox& lsfab=(*localMF[HOLD_LS_DATA_MF])[mfi];
   FArrayBox& burnvelfab=(*localMF[local_mf])[mfi];
   if (burnvelfab.nComp()==ncomp) {
    // do nothing
   } else {
    amrex::Error("burnvelfab.nComp() invalid");
   }

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // burnvelfab=BURNING_VELOCITY_MF or 
    // burnvelfab=SATURATION_TEMP_MF is cell centered.
    // sets the burning velocity/saturation temp flag from 0 to 2 if
    // foot of characteristic within range.
    // in: MASS_TRANSFER_3D.F90
   fort_extend_burning_vel( 
    &velflag,
    &level,
    &finest_level,
    xlo,dx,
    &ncomp,
    &ngrow_distance,
    tilelo,tilehi,
    fablo,fabhi,
    &bfact,
    burnvelfab.dataPtr(),
    ARLIM(burnvelfab.loVect()),ARLIM(burnvelfab.hiVect()),
    lsfab.dataPtr(),ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()));
  } // mfi
} // omp
  ns_reconcile_d_num(LOOP_EXTEND_BURNING_VEL,"level_phase_change_rate_extend");

  scompBC_map.resize(ncomp);

  for (int imdest=0;imdest<ncomp;imdest++)
   scompBC_map[imdest]=extend_start_pos+imdest;

    //scomp=0
  PCINTERP_fill_borders(local_mf,ngrow_distance,
   0,ncomp,State_Type,scompBC_map);

  if (1==0) {
   int gridno=0;
   const Box& fabgrid = grids[gridno];
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();
   const Real* xlo = grid_loc[gridno].lo();
   int interior_only=0;
   FArrayBox& burnvelfab=(*localMF[local_mf])[0];
   const Real* dxplot = geom.CellSize();
   int scomp=0;
   int dirplot=-1;
   int id=0;
   tecplot_debug(burnvelfab,xlo,fablo,fabhi,dxplot,dirplot,id,
     scomp,ncomp,interior_only);
  }

 } // velflag=0,1 (velflag==0: interface temp/massfrac;velflag==1: burnvel

} // end subroutine level_phase_change_rate_extend


void
NavierStokes::level_DRAG_extend() {

 std::string local_caller_string="level_DRAG_extend";

 bool use_tiling=ns_tiling;
 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid level_DRAG_extend");

 const Real* dx = geom.CellSize();

 MultiFab& LS_new = get_new_data(LS_Type,slab_step+1);
 if (LS_new.nComp()!=num_materials*(AMREX_SPACEDIM+1)) 
  amrex::Error("LS_new invalid ncomp");

 if (localMF[DRAG_MF]->nComp()!=N_DRAG)
  amrex::Error("localMF[DRAG_MF] incorrect ncomp");

 if (ngrow_make_distance!=3)
  amrex::Error("expecting ngrow_make_distance==3");

 if (localMF[DRAG_MF]->nGrow()!=ngrow_make_distance)
  amrex::Error("localMF[DRAG_MF] incorrect ngrow");

 debug_ngrow(HOLD_LS_DATA_MF,ngrow_distance,local_caller_string);
 if (localMF[HOLD_LS_DATA_MF]->nComp()!=num_materials*(1+AMREX_SPACEDIM)) 
  amrex::Error("localMF[HOLD_LS_DATA_MF]->nComp() invalid");

 int ncomp=N_DRAG;

 Vector<int> scompBC_map;
 scompBC_map.resize(ncomp);
 int extend_start_pos=EXTRAPCOMP_DRAG;

 for (int imdest=0;imdest<ncomp;imdest++)
  scompBC_map[imdest]=extend_start_pos+imdest;

  //scomp=0
 PCINTERP_fill_borders(DRAG_MF,ngrow_make_distance,
   0,ncomp,State_Type,scompBC_map);

 if (1==0) {
  int gridno=0;
  const Box& fabgrid = grids[gridno];
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();
  const Real* xlo = grid_loc[gridno].lo();
  int interior_only=0;
  FArrayBox& dragfab=(*localMF[DRAG_MF])[0];
  const Real* dxplot = geom.CellSize();
  int scomp=0;
  int dirplot=-1;
  int id=0;
  tecplot_debug(dragfab,xlo,fablo,fabhi,dxplot,dirplot,id,
     scomp,ncomp,interior_only);
 }

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(LS_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(LS_new,use_tiling); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();
  int bfact=parent->Space_blockingFactor(level);

  const Real* xlo = grid_loc[gridno].lo();
  FArrayBox& lsfab=(*localMF[HOLD_LS_DATA_MF])[mfi];
  FArrayBox& dragfab=(*localMF[DRAG_MF])[mfi];
  if (dragfab.nComp()==ncomp) {
   // do nothing
  } else {
   amrex::Error("dragfab.nComp() invalid");
  }

  int ngrow=ngrow_distance;
  if (ngrow!=4)
   amrex::Error("expecting ngrow==4");

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // DRAG_MF is cell centered.
   // sets the flag from 0 to 2 if
   // foot of characteristic within range.
   // declared in: MASS_TRANSFER_3D.F90
  fort_extend_drag( 
   &level,
   &finest_level,
   xlo,dx,
   &ncomp,
   &ngrow,
   tilelo,tilehi,
   fablo,fabhi,
   &bfact,
   dragfab.dataPtr(),
   ARLIM(dragfab.loVect()),ARLIM(dragfab.hiVect()),
   lsfab.dataPtr(),ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()));
 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_EXTEND_DRAG,"level_DRAG_extend");

 scompBC_map.resize(ncomp);

 for (int imdest=0;imdest<ncomp;imdest++)
  scompBC_map[imdest]=extend_start_pos+imdest;

   //scomp=0
 PCINTERP_fill_borders(DRAG_MF,ngrow_make_distance,
   0,ncomp,State_Type,scompBC_map);

 if (1==0) {
   int gridno=0;
   const Box& fabgrid = grids[gridno];
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();
   const Real* xlo = grid_loc[gridno].lo();
   int interior_only=0;
   FArrayBox& dragfab=(*localMF[DRAG_MF])[0];
   const Real* dxplot = geom.CellSize();
   int scomp=0;
   int dirplot=-1;
   int id=0;
   tecplot_debug(dragfab,xlo,fablo,fabhi,dxplot,dirplot,id,
     scomp,ncomp,interior_only);
 }


} // end subroutine level_DRAG_extend


void
NavierStokes::level_phase_change_convertALL() {

 std::string local_caller_string="level_phase_change_convertALL";

 int finest_level=parent->finestLevel();
 if (level==0) {
  // do nothing
 } else
  amrex::Error("level must be 0");

 int nden=num_materials*num_state_material;

 int iten;
 int im;
 int im_opp;

 int n_phase_change=0;
 for (im=1;im<=num_materials-1;im++) {
  for (im_opp=im+1;im_opp<=num_materials;im_opp++) {
   get_iten_cpp(im,im_opp,iten); //declared in NavierStokes2.cpp
   if ((iten<1)||(iten>num_interfaces))
    amrex::Error("iten invalid");
   Real LL0=get_user_latent_heat(iten,293.0,1);
   Real LL1=get_user_latent_heat(iten+num_interfaces,293.0,1);
   if ((LL0!=0.0)||(LL1!=0.0)) {
    n_phase_change++;
   } else if ((LL0==0.0)&&(LL1==0.0)) {
    // do nothing
   } else
    amrex::Error("LL0 or LL1 invalid");
  } // im_opp=im+1..num_materials
 } // im=1..num_materials-1

 debug_ngrow(HOLD_LS_DATA_MF,ngrow_distance,local_caller_string);
 if (localMF[HOLD_LS_DATA_MF]->nComp()!=num_materials*(1+AMREX_SPACEDIM)) 
  amrex::Error("localMF[HOLD_LS_DATA_MF]->nComp() invalid");

 debug_ngrow(DEN_RECON_MF,1,local_caller_string);
 if (localMF[DEN_RECON_MF]->nComp()!=nden)
  amrex::Error("DEN_RECON_MF invalid ncomp");

 int i_phase_change=0;
 while (i_phase_change<n_phase_change) {

  for (im=1;im<=num_materials-1;im++) {
   for (im_opp=im+1;im_opp<=num_materials;im_opp++) {
    get_iten_cpp(im,im_opp,iten); //declared in NavierStokes2.cpp
    if ((iten<1)||(iten>num_interfaces))
     amrex::Error("iten invalid");
    Real LL0=get_user_latent_heat(iten,293.0,1);
    Real LL1=get_user_latent_heat(iten+num_interfaces,293.0,1);
    if ((LL0!=0.0)||(LL1!=0.0)) {

     for (int ilev=finest_level;ilev>=level;ilev--) {
      NavierStokes& ns_level=getLevel(ilev);
       // unsplit advection using:
       // BURNING_VELOCITY_MF (interpolated to the nodes)
       // (i) updates centroid and volume fraction for
       //   im and im_opp materials.
       // (ii) updates: JUMP_STRENGTH_MF  (rho_1/rho_2  - 1) expansion factor
       // if i_phase_change+1<n_phase_change
       //  a) copies LS_new[im_opp-1],LS_new[im-1] to HOLD_LS_DATA.
       //  b) copies S_new to DEN_RECON_MF for density, temperature, Y
       //     (im, im_opp)
      ns_level.level_phase_change_convert(im,im_opp,
       i_phase_change,n_phase_change);

       // spectral_override==0 => always low order.
      ns_level.avgDown(LS_Type,0,num_materials,LOW_ORDER_AVGDOWN);
      ns_level.MOFavgDown();
      ns_level.avgDown(State_Type,STATECOMP_STATES,
	num_state_material*num_materials,SPECTRAL_ORDER_AVGDOWN);
     } // ilev=finest_level ... level

     if (i_phase_change+1<n_phase_change) {

      for (int im_count=0;im_count<2;im_count++) {
       int im_current=0;
       if (im_count==0) {
        im_current=im;
       } else if (im_count==1) {
        im_current=im_opp;
       } else
        amrex::Error("im_count invalid");

       Vector<int> scompBC_map;
       scompBC_map.resize(1);
       debug_ngrow(DEN_RECON_MF,1,local_caller_string);

       int dstcomp=(im_current-1)*num_state_material;
       int srccomp=STATECOMP_STATES+(im_current-1)*num_state_material;

       // density
       // spectral_override==0 => always low order
       avgDown_localMF_ALL(DEN_RECON_MF,dstcomp+ENUM_DENVAR,1,
	  SPECTRAL_ORDER_AVGDOWN);
       scompBC_map[0]=srccomp+ENUM_DENVAR;
        //ngrow=1 scomp=dstcomp+ENUM_DENVAR ncomp=1
       GetStateFromLocalALL(DEN_RECON_MF,1,
	  dstcomp+ENUM_DENVAR,1,State_Type,scompBC_map);

       // temperature
       avgDown_localMF_ALL(DEN_RECON_MF,dstcomp+ENUM_TEMPERATUREVAR,1,1);
       scompBC_map[0]=srccomp+ENUM_TEMPERATUREVAR;
        //ngrow=1 scomp=dstcomp+ENUM_TEMPERATUREVAR ncomp=1
       GetStateFromLocalALL(DEN_RECON_MF,1,
	  dstcomp+ENUM_TEMPERATUREVAR,1,State_Type,scompBC_map);

       int ispec=mass_fraction_id[iten-1];
       if (ispec==0) {
        // do nothing
       } else if ((ispec>=1)&&(ispec<=num_species_var)) {
        avgDown_localMF_ALL(DEN_RECON_MF,dstcomp+ENUM_SPECIESVAR+ispec-1,1,1);
        scompBC_map[0]=srccomp+ENUM_SPECIESVAR+ispec-1;
        GetStateFromLocalALL(DEN_RECON_MF,1,
 	   dstcomp+ENUM_SPECIESVAR+ispec-1,1,State_Type,scompBC_map);
       } else
        amrex::Error("ispec invalid");
      } // im_count=0,1

      Vector<int> scompBC_map_LS;
      scompBC_map_LS.resize(num_materials*(AMREX_SPACEDIM+1));
       // spectral_override==0 => always low order
      avgDown_localMF_ALL(HOLD_LS_DATA_MF,0,num_materials*(AMREX_SPACEDIM+1),
		  LOW_ORDER_AVGDOWN);
      for (int im_group=0;im_group<num_materials*(AMREX_SPACEDIM+1);im_group++)
       scompBC_map_LS[im_group]=im_group;
      debug_ngrow(HOLD_LS_DATA_MF,ngrow_distance,local_caller_string);
       //scomp=0
      GetStateFromLocalALL(HOLD_LS_DATA_MF,ngrow_distance,
         0,num_materials*(AMREX_SPACEDIM+1),LS_Type,scompBC_map_LS);

      int update_flag=RECON_UPDATE_STATE_CENTROID;
      int init_vof_prev_time=0;
        // Fluids tessellate; solids overlay; output:SLOPE_RECON_MF
      VOF_Recon_ALL(1,cur_time_slab,update_flag,init_vof_prev_time);
     } else if (i_phase_change+1==n_phase_change) {
      // do nothing
     } else {
      amrex::Error("i_phase_change invalid");
     }
     i_phase_change++;
    } else if ((LL0==0.0)&&(LL1==0.0)) {
     // do nothing
    } else
     amrex::Error("LL0 or LL1 invalid");
   } // im_opp
  } // im
 } // while i_phase_change<n_phase_change

 if (i_phase_change==n_phase_change) {
  // do nothing
 } else
  amrex::Error("i_phase_change invalid");

} // end subroutine level_phase_change_convertALL

// 1. initialize node velocity from BURNING_VELOCITY_MF
// 2. unsplit advection of materials changing phase
// 3. update volume fractions, jump strength, temperature,
//    species
void
NavierStokes::level_phase_change_convert(
  int im_outer,int im_opp_outer,
  int i_phase_change,int n_phase_change) {

 std::string local_caller_string="level_phase_change_convert";

 bool use_tiling=ns_tiling;
 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid level_phase_change_convert");

 if (ngrow_distance==4) {
  // do nothing
 } else
  amrex::Error("expecting ngrow_distance==4");

 int iten;
 get_iten_cpp(im_outer,im_opp_outer,iten); //declared in NavierStokes2.cpp

 if ((im_outer>=1)&&(im_outer<=num_materials)&&
     (im_opp_outer>=1)&&(im_opp_outer<=num_materials)&&
     (im_outer<im_opp_outer)&&
     (i_phase_change>=0)&&
     (i_phase_change<num_interfaces)&&
     (i_phase_change<n_phase_change)&&
     (iten>=1)&&(iten<=num_interfaces)) {

  Real LL0=get_user_latent_heat(iten,293.0,1);
  Real LL1=get_user_latent_heat(iten+num_interfaces,293.0,1);
  if ((LL0!=0.0)||(LL1!=0.0)) {
   // do nothing
  } else
   amrex::Error("LL0 or LL1 invalid");
 } else
  amrex::Error("level_phase_change_convert: invalid parameters");
 
  // first num_interfaces components are the status
 int nburning=EXTRAP_NCOMP_BURNING;
 int ntsat=EXTRAP_NCOMP_TSAT;

 int nden=num_materials*num_state_material;
 int nstate=STATE_NCOMP;

 // mask=1 if not covered or if outside the domain.
 // NavierStokes::maskfiner_localMF
 // NavierStokes::maskfiner
 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,1,local_caller_string); 

 if (localMF[JUMP_STRENGTH_MF]->nGrow()!=ngrow_distance)
  amrex::Error("jump strength invalid ngrow level_phase_change_conv");
 if (localMF[JUMP_STRENGTH_MF]->nComp()!=2*num_interfaces)
  amrex::Error("localMF[JUMP_STRENGTH_MF]->nComp() invalid");

 if (localMF[SATURATION_TEMP_MF]->nComp()!=ntsat)
  amrex::Error("localMF[SATURATION_TEMP_MF]->nComp()!=ntsat");
 if (localMF[SATURATION_TEMP_MF]->nGrow()!=ngrow_distance)
  amrex::Error("localMF[SATURATION_TEMP_MF] incorrect ngrow");

 if (localMF[nodevel_MF]->nGrow()!=1)
  amrex::Error("localMF[nodevel_MF]->nGrow()  invalid");
 if (localMF[nodevel_MF]->nComp()!=2*num_interfaces*AMREX_SPACEDIM)
  amrex::Error("localMF[nodevel_MF]->nComp()  invalid");
 debug_ixType(nodevel_MF,-1,local_caller_string);

 if (localMF[BURNING_VELOCITY_MF]->nComp()!=nburning)
  amrex::Error("burning vel invalid ncomp");
 debug_ixType(BURNING_VELOCITY_MF,-1,local_caller_string);

  // this routine: level_phase_change_convert
  // DEN_RECON_MF is initialized prior to the call
  // to this routine.
 debug_ngrow(DEN_RECON_MF,1,local_caller_string);
 if (localMF[DEN_RECON_MF]->nComp()!=nden)
  amrex::Error("DEN_RECON_MF invalid ncomp");

 debug_ngrow(HOLD_LS_DATA_MF,ngrow_distance,local_caller_string);
 if (localMF[HOLD_LS_DATA_MF]->nComp()!=num_materials*(1+AMREX_SPACEDIM)) 
  amrex::Error("localMF[HOLD_LS_DATA_MF]->nComp() invalid");

 MultiFab& LS_new = get_new_data(LS_Type,slab_step+1);
 if (LS_new.nComp()!=num_materials*(AMREX_SPACEDIM+1)) 
  amrex::Error("LS_new invalid ncomp");

 MultiFab& S_new = get_new_data(State_Type,slab_step+1);
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");

 const Real* dx = geom.CellSize();

 Vector< Vector<Real> > delta_mass_local;
 delta_mass_local.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  // source 1..num_materials  dest 1..num_materials
  delta_mass_local[tid].resize(2*num_materials); 
  for (int im=0;im<2*num_materials;im++)
   delta_mass_local[tid][im]=0.0;
 } // tid


 if (i_phase_change==0) {

  if (level==finest_level) {

   for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
    parent->AMR_max_phase_change_rate[dir]=0.0;
    parent->AMR_min_phase_change_rate[dir]=0.0;
   }

   for (int iten_local=0;iten_local<num_interfaces;iten_local++) {
    int scomp=num_interfaces+iten_local*AMREX_SPACEDIM;
    for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
     Real local_max=
      localMF[BURNING_VELOCITY_MF]->max(scomp+dir); //def nghost=0
     parent->AMR_max_phase_change_rate[dir]=
      max(parent->AMR_max_phase_change_rate[dir],local_max);

     Real local_min=
      localMF[BURNING_VELOCITY_MF]->min(scomp+dir); //def nghost=0
     parent->AMR_min_phase_change_rate[dir]=
      min(parent->AMR_min_phase_change_rate[dir],local_min);
    } //dir=0..sdim-1
   } //iten_local=0..num_interfaces-1

  }  // level==finest_level

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(S_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(S_new,use_tiling); mfi.isValid(); ++mfi) {
    BL_ASSERT(grids[mfi.index()] == mfi.validbox());
    const int gridno = mfi.index();
    const Box& tilegrid = mfi.tilebox();
    const Box& fabgrid = grids[gridno];
    const int* tilelo=tilegrid.loVect();
    const int* tilehi=tilegrid.hiVect();
    const int* fablo=fabgrid.loVect();
    const int* fabhi=fabgrid.hiVect();

    const Real* xlo = grid_loc[gridno].lo();
    Vector<int> velbc=getBCArray(State_Type,gridno,
        STATECOMP_VEL,STATE_NCOMP_VEL);

    FArrayBox& burnvelfab=(*localMF[BURNING_VELOCITY_MF])[mfi];
    FArrayBox& nodevelfab=(*localMF[nodevel_MF])[mfi];
    if (burnvelfab.nComp()==nburning) {
     // do nothing
    } else 
     amrex::Error("burnvelfab.nComp() invalid");

    if (nodevelfab.nComp()==2*num_interfaces*AMREX_SPACEDIM) {
     // do nothing
    } else 
     amrex::Error("nodevelfab.nComp() invalid");

    FArrayBox& olddistfab=(*localMF[HOLD_LS_DATA_MF])[mfi];

    int bfact=parent->Space_blockingFactor(level);

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

      // burnvelfab=BURNING_VELOCITY_MF is cell centered (interior: lo to hi)
      // nodevelfab=nodevel is at the nodes. (interior: lo to hi+1)
      // fort_nodedisplace is declared in: MASS_TRANSFER_3D.F90
    fort_nodedisplace(
     &nburning,
     tilelo,tilehi,
     fablo,fabhi,
     &bfact, 
     velbc.dataPtr(),
     &dt_slab, //nodedisplace
     nodevelfab.dataPtr(),
     ARLIM(nodevelfab.loVect()),ARLIM(nodevelfab.hiVect()),
     burnvelfab.dataPtr(),
     ARLIM(burnvelfab.loVect()),ARLIM(burnvelfab.hiVect()),
     olddistfab.dataPtr(),
     ARLIM(olddistfab.loVect()),ARLIM(olddistfab.hiVect()),
     xlo,dx, 
     &level,&finest_level);
  } // mfi
} // omp
  ns_reconcile_d_num(LOOP_NODEDISPLACE,"level_phase_change_convert");

  if (1==0) {
    int gridno=0;
    const Box& fabgrid = grids[gridno];
    const int* fablo=fabgrid.loVect();
    const int* fabhi=fabgrid.hiVect();
    const Real* xlo = grid_loc[gridno].lo();
    int interior_only=0;
    FArrayBox& nodevelfab=(*localMF[nodevel_MF])[0];
    const Real* dxplot = geom.CellSize();
    int scomp=0;
    int ncomp=2*num_interfaces*AMREX_SPACEDIM;
    int dirplot=-1;
    int id=0;
    std::cout << "dt_slab = " << dt_slab << '\n';
    tecplot_debug(nodevelfab,xlo,fablo,fabhi,dxplot,dirplot,id,
      scomp,ncomp,interior_only);
  }

 } else if ((i_phase_change>=1)&&(i_phase_change<n_phase_change)) {
  // do nothing
 } else
  amrex::Error("i_phase_change invalid");
  

 VOF_Recon_resize(ngrow_distance); //output:SLOPE_RECON_MF

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(S_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(S_new,use_tiling); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();

   const Real* xlo = grid_loc[gridno].lo();

   Vector<int> vofbc=getBCArray(State_Type,gridno,STATECOMP_MOF,1);

    // mask=tag if not covered by level+1 or outside the domain.
   FArrayBox& maskcov=(*localMF[MASKCOEF_MF])[mfi];

   FArrayBox& thermal_conductivity_fab=
	  (*localMF[CELL_CONDUCTIVITY_MATERIAL_MF])[mfi];
   if (thermal_conductivity_fab.nComp()!=num_materials)
    amrex::Error("thermal_conductivity_fab.nComp()!=num_materials");

   FArrayBox& nodevelfab=(*localMF[nodevel_MF])[mfi];
   if (nodevelfab.nComp()==2*num_interfaces*AMREX_SPACEDIM) {
    // do nothing
   } else 
    amrex::Error("nodevelfab.nComp() invalid");

   FArrayBox& olddistfab=(*localMF[HOLD_LS_DATA_MF])[mfi];
   FArrayBox& sweptfab=(*localMF[SWEPT_CROSSING_MF])[mfi];
   FArrayBox& lsnewfab=LS_new[mfi];
   FArrayBox& snewfab=S_new[mfi];

   FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi]; 

   FArrayBox& eosfab=(*localMF[DEN_RECON_MF])[mfi];

   FArrayBox& JUMPfab=(*localMF[JUMP_STRENGTH_MF])[mfi];

   FArrayBox& Tsatfab=(*localMF[SATURATION_TEMP_MF])[mfi];
   if (Tsatfab.nComp()!=ntsat) {
    std::cout << "Tsatfab.nComp()=" << Tsatfab.nComp() << 
      " ntsat=" << ntsat << '\n';
    amrex::Error("Tsatfab.nComp()!=ntsat 3");
   }

   int bfact=parent->Space_blockingFactor(level);

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // nodevelfab (=nodevel) is node based (interior: lo..hi+1)
    // declared in: MASS_TRANSFER_3D.F90
   fort_convertmaterial( 
    &tid_current,
    &im_outer,
    &im_opp_outer,
    &level,&finest_level,
    &nden,
    &nstate,
    &ntsat,
    saturation_temp.dataPtr(),
    freezing_model.dataPtr(),
    Tanasawa_or_Schrage_or_Kassemi.dataPtr(),
    mass_fraction_id.dataPtr(),
    distribute_from_target.dataPtr(),
    constant_density_all_time.dataPtr(),
    tilelo,tilehi,
    fablo,fabhi,
    &bfact, 
    &min_stefan_velocity_for_dt,
    vofbc.dataPtr(),
    xlo,dx,
    &dt_slab,//convertmaterial
    delta_mass_local[tid_current].dataPtr(),
    maskcov.dataPtr(),
    ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
    thermal_conductivity_fab.dataPtr(), //num_materials components
    ARLIM(thermal_conductivity_fab.loVect()),
    ARLIM(thermal_conductivity_fab.hiVect()),
    nodevelfab.dataPtr(),
    ARLIM(nodevelfab.loVect()),ARLIM(nodevelfab.hiVect()),
    JUMPfab.dataPtr(),
    ARLIM(JUMPfab.loVect()),ARLIM(JUMPfab.hiVect()),
    Tsatfab.dataPtr(),
    ARLIM(Tsatfab.loVect()),ARLIM(Tsatfab.hiVect()),
    olddistfab.dataPtr(),
    ARLIM(olddistfab.loVect()),ARLIM(olddistfab.hiVect()),
    lsnewfab.dataPtr(),
    ARLIM(lsnewfab.loVect()),ARLIM(lsnewfab.hiVect()),
    reconfab.dataPtr(),
    ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
    snewfab.dataPtr(),
    ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()),
    eosfab.dataPtr(),
    ARLIM(eosfab.loVect()),ARLIM(eosfab.hiVect()),
    sweptfab.dataPtr(),
    ARLIM(sweptfab.loVect()),ARLIM(sweptfab.hiVect()));
 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_CONVERTMATERIAL,"level_phase_change_convert");

 for (int tid=1;tid<thread_class::nthreads;tid++) {
  for (int im=0;im<2*num_materials;im++) {
   delta_mass_local[0][im]+=delta_mass_local[tid][im];
  }
 } // tid

 ParallelDescriptor::Barrier();

 for (int im=0;im<2*num_materials;im++) {
  ParallelDescriptor::ReduceRealSum(delta_mass_local[0][im]);
  delta_mass[0][im]+=delta_mass_local[0][im];
 }

// if spectral_override==0, then always low order average down.
// JUMP_STRENGTH_MF=0.0 outside the computational domain.
 int spectral_override=LOW_ORDER_AVGDOWN;
 localMF[JUMP_STRENGTH_MF]->FillBoundary(
    iten-1,1,geom.periodicity());
 avgDown_localMF(JUMP_STRENGTH_MF,iten-1,1,spectral_override);

 localMF[JUMP_STRENGTH_MF]->FillBoundary(
    num_interfaces+iten-1,1,geom.periodicity());
 avgDown_localMF(JUMP_STRENGTH_MF,num_interfaces+iten-1,1,spectral_override);

 if (i_phase_change+1<n_phase_change) {

  for (int im_count=0;im_count<2;im_count++) {
   int im_current=0;
   if (im_count==0) {
    im_current=im_outer;
   } else if (im_count==1) {
    im_current=im_opp_outer;
   } else
    amrex::Error("im_count invalid");

   // dest,src,srccomp,dstcomp,ncomp,ngrow
   MultiFab::Copy(*localMF[HOLD_LS_DATA_MF],LS_new,
		  im_current-1,im_current-1,1,0);
   int dstcomp=(im_current-1)*num_state_material;
   int srccomp=STATECOMP_STATES+(im_current-1)*num_state_material;
   // density
   MultiFab::Copy(*localMF[DEN_RECON_MF],S_new,
		  srccomp,dstcomp,1,0);
   // temperature
   MultiFab::Copy(*localMF[DEN_RECON_MF],S_new,
     srccomp+1,dstcomp+1,1,0);

   int ispec=mass_fraction_id[iten-1];
   if (ispec==0) {
    // do nothing
   } else if ((ispec>=1)&&(ispec<=num_species_var)) {
    MultiFab::Copy(*localMF[DEN_RECON_MF],S_new,
      srccomp+1+ispec,dstcomp+1+ispec,1,0);
   } else
    amrex::Error("ispec invalid");
  } // im_count=0,1

 } else if (i_phase_change+1==n_phase_change) {
  // do nothing
 } else {
  amrex::Error("i_phase_change invalid");
 }

} // subroutine level_phase_change_convert


void
NavierStokes::level_species_reaction(const std::string& caller_string) {

 std::string local_caller_string="level_species_reaction";
 local_caller_string=caller_string+local_caller_string;

 int initialize_flag=0;
 if (pattern_test(local_caller_string,"initData")==1) {
  initialize_flag=1;
 } else if (pattern_test(local_caller_string,"veldiffuseALL")==1) {
  initialize_flag=0;
 } else
  amrex::Error("local_caller_string invalid in level_species_reaction");

 bool use_tiling=ns_tiling;
 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid level_species_reaction");

 if (ngrow_distance==4) {
  // do nothing
 } else
  amrex::Error("expecting ngrow_distance==4");

 int nstate=STATE_NCOMP;

 MultiFab& S_new = get_new_data(State_Type,slab_step+1);
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");

 MultiFab& LS_new=get_new_data(LS_Type,slab_step+1);
 if (LS_new.nComp()!=num_materials*(AMREX_SPACEDIM+1))
  amrex::Error("LS_new ncomp invalid");

 MultiFab* local_mask;
 Real local_time;
 Real local_dt;

 if (initialize_flag==1) {
  local_mask=&(S_new);
  local_time=0.0;
  local_dt=0.0;
 } else if (initialize_flag==0) {
  // mask=1 if not covered or if outside the domain.
  // NavierStokes::maskfiner_localMF
  // NavierStokes::maskfiner
  resize_maskfiner(1,MASKCOEF_MF);
  debug_ngrow(MASKCOEF_MF,1,local_caller_string); 
  local_mask=localMF[MASKCOEF_MF];
  local_time=cur_time_slab;
  local_dt=dt_slab; //apply_reaction
 } else
  amrex::Error("initialize_flag invalid");

 const Real* dx = geom.CellSize();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(S_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(S_new,use_tiling); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();

   const Real* xlo = grid_loc[gridno].lo();

   int bfact=parent->Space_blockingFactor(level);

    // mask=tag if not covered by level+1 or outside the domain.
   FArrayBox& maskcov=(*local_mask)[mfi];
   FArrayBox& snewfab=S_new[mfi];
   FArrayBox& LSnewfab=LS_new[mfi];

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    //fort_apply_reaction is declared in: MASS_TRANSFER_3D.F90
   fort_apply_reaction(
    &tid_current,
    &level,&finest_level,
    &nstate,
    speciesreactionrate.dataPtr(),
    rigid_fraction_id.dataPtr(),
    tilelo,tilehi,
    fablo,fabhi,
    &bfact, 
    xlo,dx,
    &initialize_flag,
    &local_dt,
    &local_time,
    maskcov.dataPtr(),
    ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
    LSnewfab.dataPtr(),
    ARLIM(LSnewfab.loVect()),ARLIM(LSnewfab.hiVect()),
    snewfab.dataPtr(),
    ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()));
 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_SPECIES_REACTION,"level_species_reaction");

} // subroutine level_species_reaction


void
NavierStokes::phase_change_redistributeALL() {

 std::string local_caller_string="phase_change_redistributeALL";

 if (level!=0)
  amrex::Error("level invalid phase_change_redistributeALL");

 if (ngrow_distance==4) {
  // do nothing
 } else
  amrex::Error("expecting ngrow_distance==4");

 int finest_level=parent->finestLevel();
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.getStateDist_localMF(LSNEW_MF,ngrow_distance,cur_time_slab,
		  local_caller_string);
 }

 mdotplus.resize(thread_class::nthreads);
 mdotminus.resize(thread_class::nthreads);
 mdotcount.resize(thread_class::nthreads);
 mdot_lost.resize(thread_class::nthreads);
 mdot_sum.resize(thread_class::nthreads);
 mdot_sum2.resize(thread_class::nthreads);

 mdotplus_complement.resize(thread_class::nthreads);
 mdotminus_complement.resize(thread_class::nthreads);
 mdotcount_complement.resize(thread_class::nthreads);
 mdot_lost_complement.resize(thread_class::nthreads);
 mdot_sum_complement.resize(thread_class::nthreads);
 mdot_sum2_complement.resize(thread_class::nthreads);

 allocate_array(ngrow_distance,2*num_interfaces,-1,
		JUMP_STRENGTH_COMPLEMENT_MF); 
 Copy_array(JUMP_STRENGTH_COMPLEMENT_MF,JUMP_STRENGTH_MF,
	    0,0,2*num_interfaces,ngrow_distance);

 for (int im=1;im<=num_materials;im++) {
  for (int im_opp=im+1;im_opp<=num_materials;im_opp++) {
   for (int ireverse=0;ireverse<=1;ireverse++) {
    if ((im>num_materials)||(im_opp>num_materials))
     amrex::Error("im or im_opp bust 200cpp");
    int iten;
    get_iten_cpp(im,im_opp,iten); //declared in NavierStokes2.cpp
    if ((iten<1)||(iten>num_interfaces))
     amrex::Error("iten invalid");

    int indexEXP=iten+ireverse*num_interfaces-1;

    Real LL=get_user_latent_heat(indexEXP+1,293.0,1);
    int distribute_from_targ=distribute_from_target[indexEXP];

    int im_source=-1;
    int im_dest=-1;

    if ((ns_is_rigid(im-1)==1)||
        (ns_is_rigid(im_opp-1)==1)) { 
     // do nothing
    } else if (LL!=0.0) {

     if (ireverse==0) {
      im_source=im;im_dest=im_opp;
     } else if (ireverse==1) {
      im_source=im_opp;im_dest=im;
     } else
      amrex::Error("ireverse invalid");

     Real expect_mdot_sign=0.0;
     if ((distribute_from_targ==0)||
         (distribute_from_targ==1)) {
      expect_mdot_sign=( 
       (denconst[im_source-1]>denconst[im_dest-1]) ? 1.0 : -1.0 );
     } else
      amrex::Error("distribute_from_targ invalid");

     for (int tid=0;tid<thread_class::nthreads;tid++) {
      mdotplus[tid]=0.0;
      mdotminus[tid]=0.0;
      mdotcount[tid]=0.0;
      mdot_sum[tid]=0.0;
      mdot_sum2[tid]=0.0;
      mdot_lost[tid]=0.0;

      mdotplus_complement[tid]=0.0;
      mdotminus_complement[tid]=0.0;
      mdotcount_complement[tid]=0.0;
      mdot_sum_complement[tid]=0.0;
      mdot_sum2_complement[tid]=0.0;
      mdot_lost_complement[tid]=0.0;
     }

     allocate_array(ngrow_distance,1,-1,donorflag_MF);
      //ngrow,scomp,ncomp
     setVal_array(ngrow_distance,0,1,0.0,donorflag_MF);

     allocate_array(ngrow_distance,1,-1,accept_weight_MF);
      //ngrow,scomp,ncomp
     setVal_array(ngrow_distance,0,1,0.0,accept_weight_MF);

     allocate_array(ngrow_distance,1,-1,donorflag_complement_MF);
      //ngrow,scomp,ncomp
     setVal_array(ngrow_distance,0,1,0.0,donorflag_complement_MF);

     allocate_array(ngrow_distance,1,-1,accept_weight_complement_MF);
      //ngrow,scomp,ncomp
     setVal_array(ngrow_distance,0,1,0.0,accept_weight_complement_MF);

      // isweep==0: fort_tagexpansion
      // isweep==1: fort_accept_weight
      // isweep==2: fort_distributeexpansion
      // isweep==3: fort_initjumpterm
     for (int isweep_redistribute=0;isweep_redistribute<=2;
	  isweep_redistribute++) {

      for (int ilev=finest_level;ilev>=level;ilev--) {
       NavierStokes& ns_level=getLevel(ilev);
       ns_level.level_phase_change_redistribute(
        expect_mdot_sign,im_source,im_dest,indexEXP,
        isweep_redistribute);
      } // ilev=finest_level ... level

      // idx,ngrow,scomp,ncomp,index,scompBC_map
      Vector<int> scompBC_map;
      scompBC_map.resize(1);
      scompBC_map[0]=0; //set_extrap_bc, fort_extrapfill

      if (isweep_redistribute==0) {
       PCINTERP_fill_bordersALL(donorflag_MF,1,0,
         1,State_Type,scompBC_map);
       PCINTERP_fill_bordersALL(donorflag_complement_MF,1,0,
         1,State_Type,scompBC_map);
      } else if (isweep_redistribute==1) {
       PCINTERP_fill_bordersALL(accept_weight_MF,1,0,
         1,State_Type,scompBC_map);
       PCINTERP_fill_bordersALL(accept_weight_complement_MF,1,0,
         1,State_Type,scompBC_map);
      } else if (isweep_redistribute==2) { //fort_distributeexpansion
       // do nothing
      } else
       amrex::Error("isweep_redistribute invalid");

      if (ParallelDescriptor::IOProcessor()) {

       if (isweep_redistribute==0) { //fort_tagexpansion
        std::cout << "before:imsrc,imdst,mdot_sum " <<
         im_source << ' ' << im_dest << ' ' << mdot_sum[0] << '\n';
        std::cout << "before:imsrc,imdst,mdot_sum_complement " <<
         im_source << ' ' << im_dest << ' ' << mdot_sum_complement[0] << '\n';
       } else if (isweep_redistribute==1) { //fort_accept_weight
        // do nothing
       } else if (isweep_redistribute==2) { //fort_distributeexpansion
        std::cout << "after:imsrc,imdst,mdot_sum2 " <<   
         im_source << ' ' << im_dest << ' ' << mdot_sum2[0] << '\n';
        std::cout << "after:imsrc,imdst,mdot_lost " <<   
         im_source << ' ' << im_dest << ' ' << mdot_lost[0] << '\n';
        std::cout << "imsrc,imdst,mdot_sum2+mdot_lost " <<   
         im_source << ' ' << im_dest << ' ' <<   
         mdot_sum2[0]+mdot_lost[0] << '\n';

        std::cout << "after:imsrc,imdst,mdot_sum2_complement " <<   
         im_source << ' ' << im_dest << ' ' << mdot_sum2_complement[0] << '\n';
        std::cout << "after:imsrc,imdst,mdot_lost_complement " <<   
         im_source << ' ' << im_dest << ' ' << mdot_lost_complement[0] << '\n';
        std::cout << 
         "imsrc,imdst,mdot_sum2_complement+mdot_lost_complement " <<   
         im_source << ' ' << im_dest << ' ' <<   
         mdot_sum2_complement[0]+mdot_lost_complement[0] << '\n';

       } else
        amrex::Error("isweep_redistribute invalid");

      } // if ParallelDescriptor::IOProcessor()

     } // isweep_redistribute=0,1,2

     delete_array(donorflag_MF);
     delete_array(accept_weight_MF);
     delete_array(donorflag_complement_MF);
     delete_array(accept_weight_complement_MF);

    } else if (LL==0.0) {
     // do nothing
    } else {
     amrex::Error("LL is NaN");
    } 
   } // ireverse
  } // im_opp
 } // im=1..num_materials


 Vector<blobclass> blobdata;
 Vector< Vector<Real> > mdot_data;
 Vector< Vector<Real> > mdot_comp_data;
 Vector< Vector<Real> > mdot_data_redistribute;
 Vector< Vector<Real> > mdot_comp_data_redistribute;
 Vector<int> type_flag;

 int color_count=0;
 int coarsest_level=0;
  //idx_mdot,idx_mdot_complement >=0 signals ColorSumALL to compute
  //auxiliary sums.
 int idx_mdot=JUMP_STRENGTH_MF;
 int idx_mdot_complement=JUMP_STRENGTH_COMPLEMENT_MF;
 int tessellate=3;
 int operation_flag=OP_GATHER_MDOT; 
 ColorSumALL(
  operation_flag, // =OP_GATHER_MDOT
  tessellate,  //=3
  coarsest_level,
  color_count,
  TYPE_MF,
  COLOR_MF,
  idx_mdot,
  idx_mdot_complement,
  type_flag,
  blobdata,
  mdot_data,
  mdot_comp_data,
  mdot_data_redistribute,
  mdot_comp_data_redistribute
  );

 operation_flag=OP_SCATTER_MDOT; // scatter to mdot or density

 ColorSumALL(
  operation_flag, //=OP_SCATTER_MDOT
  tessellate,  //=3
  coarsest_level,
  color_count,
  TYPE_MF,
  COLOR_MF,
  idx_mdot,
  idx_mdot_complement,
  type_flag,
  blobdata,
  mdot_data,
  mdot_comp_data,
  mdot_data_redistribute,
  mdot_comp_data_redistribute
  );

 if (mdot_data.size()==mdot_data_redistribute.size()) {
  if (mdot_data.size()==color_count) {
   // do nothing
  } else
   amrex::Error("mdot_data.size() invalid");
  if (blobdata.size()==color_count) {
   // do nothing
  } else
   amrex::Error("blobdata.size() invalid");
 } else
  amrex::Error("mdot_data or mdot_data_redistribute have wrong size");

 if (mdot_comp_data.size()==mdot_comp_data_redistribute.size()) {
  if (mdot_comp_data.size()==color_count) {
   // do nothing
  } else
   amrex::Error("mdot_comp_data.size() invalid");
  if (blobdata.size()==color_count) {
   // do nothing
  } else
   amrex::Error("blobdata.size() invalid");
 } else
  amrex::Error("mdot_comp_data or mdot_comp_data_redistribute wrong size");


 for (int i=0;i<mdot_data.size();i++) {
  if ((mdot_data[i].size()==2*num_interfaces)&&
      (mdot_data_redistribute[i].size()==2*num_interfaces)&&
      (mdot_data[i].size()==mdot_data_redistribute[i].size())) {
   // do nothing
  } else
   amrex::Error("mdot_data[i] or mdot_data_redistribute[i] have wrong size");
 }

 for (int i=0;i<mdot_comp_data.size();i++) {
  if ((mdot_comp_data[i].size()==2*num_interfaces)&&
      (mdot_comp_data_redistribute[i].size()==2*num_interfaces)&&
      (mdot_comp_data[i].size()==mdot_comp_data_redistribute[i].size())) {
   // do nothing
  } else
   amrex::Error("mdot_comp_data[i] or mdot_comp_data_redistribute[i] bad size");
 }

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout<<"color_count=" << color_count << '\n';
   std::cout<<"i=0..mdot_data.size()-1 (color_count)\n";
   std::cout<<"j=0..mdot_data_redistribute[i].size()-1 (2 * num_interfaces)\n";

   for (int i=0;i<mdot_data.size();i++) {

    int j=0;
    for (j=0;j<mdot_data_redistribute[i].size();j++) {
     std::cout << "i=" << i << " j=" << j << " im=" <<
      blobdata[i].im << 
      " blobdata[i].blob_cell_count=" <<
      blobdata[i].blob_cell_count << 
      " blobdata[i].blob_cellvol_count=" <<
      blobdata[i].blob_cellvol_count << 
      " blobdata[i].blob_mass=" <<
      blobdata[i].blob_mass << 
      " mdot_data[i][j]=" << mdot_data[i][j] << 
      " mdot_data_redistribute[i][j]=" <<
      mdot_data_redistribute[i][j] << '\n';
    } // j=0..2*num_interfaces-1
    if (j==2*num_interfaces) {
     // do nothing
    } else
     amrex::Error("j!=2*num_interfaces");

   } // i=0..color_count-1

   for (int i=0;i<mdot_comp_data.size();i++) {

    int j=0;
    for (j=0;j<mdot_comp_data_redistribute[i].size();j++) {
     std::cout << "i=" << i << " j=" << j << " im=" <<
      blobdata[i].im << 
      " blobdata[i].blob_cell_count=" <<
      blobdata[i].blob_cell_count << 
      " blobdata[i].blob_cellvol_count=" <<
      blobdata[i].blob_cellvol_count << 
      " blobdata[i].blob_mass=" <<
      blobdata[i].blob_mass << 
      " mdot_comp_data[i][j]=" << mdot_comp_data[i][j] << 
      " mdot_comp_data_redistribute[i][j]=" <<
      mdot_comp_data_redistribute[i][j] << '\n';
    } // j=0..2*num_interfaces-1
    if (j==2*num_interfaces) {
     // do nothing
    } else
     amrex::Error("j!=2*num_interfaces");

   } // i=0..color_count-1
  } else if (! ParallelDescriptor::IOProcessor()) {
   // do nothing
  } else
   amrex::Error("ParallelDescriptor::IOProcessor() invalid");
 } else if (verbose==0) {
  // do nothing
 } else
  amrex::Error("verbose invalid");

 delete_array(TYPE_MF);
 delete_array(COLOR_MF);

 // copy contributions from all materials changing phase to a single
 // source term.
 int isweep_combine=3;

 for (int tid=0;tid<thread_class::nthreads;tid++) {
  mdotplus[tid]=0.0;
  mdotminus[tid]=0.0;
  mdotcount[tid]=0.0;
  mdot_sum[tid]=0.0;
  mdot_sum2[tid]=0.0;
  mdot_lost[tid]=0.0;
 }

 Real expect_mdot_sign_filler=0.0;
 int im_source_filler=-1;
 int im_dest_filler=-1;
 int indexEXP_filler=-1;
  
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.level_phase_change_redistribute(
   expect_mdot_sign_filler,
   im_source_filler,im_dest_filler,
   indexEXP_filler,
   isweep_combine); // ==3 (fort_initjumpterm)
 } // ilev=finest_level ... level

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "after fort_initjumpterm\n";
   std::cout << "mdotplus = " << mdotplus[0] << '\n';
   std::cout << "mdotminus = " << mdotminus[0] << '\n';
   std::cout << "mdotcount = " << mdotcount[0] << '\n';
  } // IOProc?
 } // verbose>0

 delete_array(LSNEW_MF);
 delete_array(HOLD_LS_DATA_MF);
 delete_array(JUMP_STRENGTH_COMPLEMENT_MF);

} // end subroutine phase_change_redistributeALL

// isweep==0: fort_tagexpansion
// isweep==1: fort_accept_weight
// isweep==2: fort_distributeexpansion
// isweep==3: fort_initjumpterm
void
NavierStokes::level_phase_change_redistribute(
 Real expect_mdot_sign,
 int im_source,int im_dest,int indexEXP,
 int isweep) {

 std::string local_caller_string="level_phase_change_redistribute";

 bool use_tiling=ns_tiling;
 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid level_phase_change_redistribute");

 if (ngrow_distance==4) {
  // do nothing
 } else
  amrex::Error("expecting ngrow_distance==4");

 debug_ngrow(JUMP_STRENGTH_MF,ngrow_distance,local_caller_string);
 if (localMF[JUMP_STRENGTH_MF]->nComp()!=2*num_interfaces)
  amrex::Error("localMF[JUMP_STRENGTH_MF]->nComp()!=2*num_interfaces");

 debug_ngrow(JUMP_STRENGTH_COMPLEMENT_MF,ngrow_distance,local_caller_string);
 if (localMF[JUMP_STRENGTH_COMPLEMENT_MF]->nComp()!=2*num_interfaces)
  amrex::Error("localMF[JUMP_STRENGTH_COMPLEMENT_MF]->nComp()!=2*num_int");

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,1,local_caller_string);
 
 if (localMF[LSNEW_MF]->nComp()!=num_materials*(1+AMREX_SPACEDIM))
  amrex::Error("localMF[LSNEW_MF]->nComp() invalid");
 debug_ngrow(LSNEW_MF,ngrow_distance,local_caller_string);

 const Real* dx = geom.CellSize();
 const Box& domain = geom.Domain();
 const int* domlo = domain.loVect(); 
 const int* domhi = domain.hiVect();

  // tags for redistribution of source term
  // 1=> donor  2=> receiver  0=> neither

 Real LL=0.0;

 if ((isweep==0)|| //fort_tagexpansion
     (isweep==1)|| //fort_accept_weight
     (isweep==2)) { //fort_distributeexpansion

  if (localMF[donorflag_MF]->nGrow()!=ngrow_distance)
   amrex::Error("localMF[donorflag_MF]->ngrow() invalid");
  if (localMF[donorflag_MF]->nComp()!=1)
   amrex::Error("localMF[donorflag_MF]->nComp() invalid");

  if (localMF[donorflag_complement_MF]->nGrow()!=ngrow_distance)
   amrex::Error("localMF[donorflag_complement_MF]->ngrow() invalid");
  if (localMF[donorflag_complement_MF]->nComp()!=1)
   amrex::Error("localMF[donorflag_complement_MF]->nComp() invalid");

  if (localMF[accept_weight_MF]->nGrow()!=ngrow_distance)
   amrex::Error("localMF[accept_weight_MF]->ngrow() invalid");
  if (localMF[accept_weight_MF]->nComp()!=1)
   amrex::Error("localMF[accept_weight_MF]->nComp() invalid");

  if (localMF[accept_weight_complement_MF]->nGrow()!=ngrow_distance)
   amrex::Error("localMF[accept_weight_complement_MF]->ngrow() invalid");
  if (localMF[accept_weight_complement_MF]->nComp()!=1)
   amrex::Error("localMF[accept_weight_complement_MF]->nComp() invalid");

  if ((indexEXP>=0)&&(indexEXP<2*num_interfaces)) {
   LL=get_user_latent_heat(indexEXP+1,293.0,1);
  } else
   amrex::Error("indexEXP invalid");

 } else if (isweep==3) { //fort_initjumpterm

   // indexEXP is a filler when isweep==3.
  if (indexEXP==-1) {
   LL=0.0;
  } else
   amrex::Error("indexEXP invalid");

 } else
  amrex::Error("isweep invalid");

 VOF_Recon_resize(1); //output:SLOPE_RECON_MF
  
 if (isweep==0) { //fort_tagexpansion

  if (LL!=0.0) {
   // do nothing
  } else if (LL==0.0) {
   amrex::Error("LL invalid");
  } else
   amrex::Error("LL is NaN");

  if (std::abs(expect_mdot_sign)!=1.0)
   amrex::Error("expect_mdot_sign invalid");
  if ((im_source<1)||(im_source>num_materials))
   amrex::Error("im_source invalid");
  if ((im_dest<1)||(im_dest>num_materials))
   amrex::Error("im_dest invalid");
  if ((indexEXP<0)||(indexEXP>=2*num_interfaces))
   amrex::Error("indexEXP invalid");

  int nden=num_materials*num_state_material;
  MultiFab* state_var_mf=getStateDen(1,cur_time_slab);
  if (state_var_mf->nComp()!=nden)
   amrex::Error("state_var_mf->nComp()!=nden");
  
  Vector< Real > mdot_sum_local;
  mdot_sum_local.resize(thread_class::nthreads);
  for (int tid=0;tid<thread_class::nthreads;tid++) {
   mdot_sum_local[tid]=0.0;
  }

  Vector< Real > mdot_sum_complement_local;
  mdot_sum_complement_local.resize(thread_class::nthreads);
  for (int tid=0;tid<thread_class::nthreads;tid++) {
   mdot_sum_complement_local[tid]=0.0;
  }

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(localMF[donorflag_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(*localMF[donorflag_MF],false); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();
   const Real* xlo = grid_loc[gridno].lo();
   Vector<int> vofbc=getBCArray(State_Type,gridno,STATECOMP_MOF,1);

   FArrayBox& maskcov=(*localMF[MASKCOEF_MF])[mfi];

   FArrayBox& donorfab=(*localMF[donorflag_MF])[mfi];
   FArrayBox& JUMPfab=(*localMF[JUMP_STRENGTH_MF])[mfi];

   FArrayBox& donor_comp_fab=(*localMF[donorflag_complement_MF])[mfi];
   FArrayBox& JUMP_comp_fab=(*localMF[JUMP_STRENGTH_COMPLEMENT_MF])[mfi];


   FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi]; 
   FArrayBox& newdistfab=(*localMF[LSNEW_MF])[mfi];
   FArrayBox& denstatefab=(*state_var_mf)[mfi];

   int bfact=parent->Space_blockingFactor(level);
   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // in: GODUNOV_3D.F90
    // isweep==0 
    // A cell that is dominated by an is_rigid(num_materials,im)=1
    // material is neither a donor or a receiver.
    // donorfab is modified.
   fort_tagexpansion( 
    rigid_fraction_id.dataPtr(),
    &nden,
    freezing_model.dataPtr(),
    distribute_from_target.dataPtr(),
    &cur_time_slab,
    vofbc.dataPtr(),
    &expect_mdot_sign,
    &mdot_sum_local[tid_current],
    &mdot_sum_complement_local[tid_current],
    &im_source,
    &im_dest,
    &indexEXP,
    &level,
    &finest_level,
    tilelo,tilehi,
    fablo,fabhi,
    &bfact, 
    xlo,dx,
    &dt_slab,//tagexpansion
    maskcov.dataPtr(),
    ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
    donorfab.dataPtr(),
    ARLIM(donorfab.loVect()),ARLIM(donorfab.hiVect()),
    donor_comp_fab.dataPtr(),
    ARLIM(donor_comp_fab.loVect()),ARLIM(donor_comp_fab.hiVect()),
    JUMPfab.dataPtr(),
    ARLIM(JUMPfab.loVect()),ARLIM(JUMPfab.hiVect()),
    JUMP_comp_fab.dataPtr(),
    ARLIM(JUMP_comp_fab.loVect()),ARLIM(JUMP_comp_fab.hiVect()),
    denstatefab.dataPtr(),
    ARLIM(denstatefab.loVect()),ARLIM(denstatefab.hiVect()),
    newdistfab.dataPtr(),
    ARLIM(newdistfab.loVect()),ARLIM(newdistfab.hiVect()),
    reconfab.dataPtr(),
    ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()));
 
  } // mfi
} // omp
  ns_reconcile_d_num(LOOP_TAGEXPANSION,"level_phase_change_redistribute");

  for (int tid=1;tid<thread_class::nthreads;tid++) {
   mdot_sum_local[0]+=mdot_sum_local[tid];
  }
  ParallelDescriptor::ReduceRealSum(mdot_sum_local[0]);
  mdot_sum[0]+=mdot_sum_local[0];

  avgDown_tag_localMF(donorflag_MF);
  avgDown_tag_localMF(donorflag_complement_MF);

  delete state_var_mf;

 } else if (isweep==1) { //fort_accept_weight

   //accept_weights

  if (LL!=0.0) {
   // do nothing
  } else if (LL==0.0) {
   amrex::Error("LL invalid");
  } else
   amrex::Error("LL is NaN");

  if (std::abs(expect_mdot_sign)!=1.0)
   amrex::Error("expect_mdot_sign invalid");
  if ((im_source<1)||(im_source>num_materials))
   amrex::Error("im_source invalid");
  if ((im_dest<1)||(im_dest>num_materials))
   amrex::Error("im_dest invalid");
  if ((indexEXP<0)||(indexEXP>=2*num_interfaces))
   amrex::Error("indexEXP invalid");

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(localMF[donorflag_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(*localMF[donorflag_MF],false); mfi.isValid(); ++mfi) {
    BL_ASSERT(grids[mfi.index()] == mfi.validbox());
    const int gridno = mfi.index();
    const Box& tilegrid = mfi.tilebox();
    const Box& fabgrid = grids[gridno];
    const int* tilelo=tilegrid.loVect();
    const int* tilehi=tilegrid.hiVect();
    const int* fablo=fabgrid.loVect();
    const int* fabhi=fabgrid.hiVect();
    const Real* xlo = grid_loc[gridno].lo();
    Vector<int> vofbc=getBCArray(State_Type,gridno,STATECOMP_MOF,1);

    FArrayBox& maskcov=(*localMF[MASKCOEF_MF])[mfi];

    FArrayBox& donorfab=(*localMF[donorflag_MF])[mfi];
    FArrayBox& donor_comp_fab=(*localMF[donorflag_complement_MF])[mfi];

    FArrayBox& weightfab=(*localMF[accept_weight_MF])[mfi];
    FArrayBox& weight_comp_fab=(*localMF[accept_weight_complement_MF])[mfi];

    FArrayBox& JUMPfab=(*localMF[JUMP_STRENGTH_MF])[mfi];
    FArrayBox& JUMP_comp_fab=(*localMF[JUMP_STRENGTH_COMPLEMENT_MF])[mfi];
    FArrayBox& newdistfab=(*localMF[LSNEW_MF])[mfi];

    int bfact=parent->Space_blockingFactor(level);

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     // isweep==1
     // declared in: GODUNOV_3D.F90
     // weightfab and weight_comp_fab are modified.
    fort_accept_weight( 
     &im_source,
     &im_dest,
     &indexEXP,
     &level,&finest_level,
     domlo,domhi, 
     tilelo,tilehi,
     fablo,fabhi,
     &bfact, 
     xlo,dx,
     &dt_slab,//fort_accept_weight
     maskcov.dataPtr(),
     ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
     newdistfab.dataPtr(),
     ARLIM(newdistfab.loVect()),ARLIM(newdistfab.hiVect()),
     donorfab.dataPtr(),
     ARLIM(donorfab.loVect()),ARLIM(donorfab.hiVect()),
     donor_comp_fab.dataPtr(),
     ARLIM(donor_comp_fab.loVect()),
     ARLIM(donor_comp_fab.hiVect()),
     weightfab.dataPtr(),
     ARLIM(weightfab.loVect()),ARLIM(weightfab.hiVect()),
     weight_comp_fab.dataPtr(),
     ARLIM(weight_comp_fab.loVect()),
     ARLIM(weight_comp_fab.hiVect()),
     JUMPfab.dataPtr(),
     ARLIM(JUMPfab.loVect()),
     ARLIM(JUMPfab.hiVect()),
     JUMP_comp_fab.dataPtr(),
     ARLIM(JUMP_comp_fab.loVect()),
     ARLIM(JUMP_comp_fab.hiVect()) );
  } // mfi
} //omp
  ns_reconcile_d_num(LOOP_ACCEPT_WEIGHT,"level_phase_change_redistribute");

   // spectral_override==0 => always low order.
  avgDown_localMF(accept_weight_MF,0,1,LOW_ORDER_AVGDOWN);
  avgDown_localMF(accept_weight_complement_MF,0,1,LOW_ORDER_AVGDOWN);

 } else if (isweep==2) { //fort_distributeexpansion

   // redistribution.

  if (LL!=0.0) {
   // do nothing
  } else if (LL==0.0) {
   amrex::Error("LL invalid");
  } else
   amrex::Error("LL is NaN");

  if (std::abs(expect_mdot_sign)!=1.0)
   amrex::Error("expect_mdot_sign invalid");
  if ((im_source<1)||(im_source>num_materials))
   amrex::Error("im_source invalid");
  if ((im_dest<1)||(im_dest>num_materials))
   amrex::Error("im_dest invalid");
  if ((indexEXP<0)||(indexEXP>=2*num_interfaces))
   amrex::Error("indexEXP invalid");

  Vector< Real > mdot_lost_local;
  Vector< Real > mdot_sum2_local;
  mdot_lost_local.resize(thread_class::nthreads);
  mdot_sum2_local.resize(thread_class::nthreads);
  for (int tid=0;tid<thread_class::nthreads;tid++) {
   mdot_sum2_local[tid]=0.0;
   mdot_lost_local[tid]=0.0;
  }

  Vector< Real > mdot_lost_complement_local;
  Vector< Real > mdot_sum2_complement_local;
  mdot_lost_complement_local.resize(thread_class::nthreads);
  mdot_sum2_complement_local.resize(thread_class::nthreads);
  for (int tid=0;tid<thread_class::nthreads;tid++) {
   mdot_sum2_complement_local[tid]=0.0;
   mdot_lost_complement_local[tid]=0.0;
  }

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(localMF[donorflag_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(*localMF[donorflag_MF],false); mfi.isValid(); ++mfi) {
    BL_ASSERT(grids[mfi.index()] == mfi.validbox());
    const int gridno = mfi.index();
    const Box& tilegrid = mfi.tilebox();
    const Box& fabgrid = grids[gridno];
    const int* tilelo=tilegrid.loVect();
    const int* tilehi=tilegrid.hiVect();
    const int* fablo=fabgrid.loVect();
    const int* fabhi=fabgrid.hiVect();
    const Real* xlo = grid_loc[gridno].lo();
    Vector<int> vofbc=getBCArray(State_Type,gridno,STATECOMP_MOF,1);

    FArrayBox& maskcov=(*localMF[MASKCOEF_MF])[mfi];

    FArrayBox& donorfab=(*localMF[donorflag_MF])[mfi];
    FArrayBox& donor_comp_fab=(*localMF[donorflag_complement_MF])[mfi];

    FArrayBox& weightfab=(*localMF[accept_weight_MF])[mfi];
    FArrayBox& weight_comp_fab=(*localMF[accept_weight_complement_MF])[mfi];

    FArrayBox& JUMPfab=(*localMF[JUMP_STRENGTH_MF])[mfi];
    FArrayBox& JUMP_comp_fab=(*localMF[JUMP_STRENGTH_COMPLEMENT_MF])[mfi];
    FArrayBox& newdistfab=(*localMF[LSNEW_MF])[mfi];

    int bfact=parent->Space_blockingFactor(level);

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     // isweep==2
     // declared in: GODUNOV_3D.F90
     // JUMPfab is modified.
    fort_distributeexpansion( 
     &mdot_sum2_local[tid_current],
     &mdot_lost_local[tid_current],
     &mdot_sum2_complement_local[tid_current],
     &mdot_lost_complement_local[tid_current],
     &im_source,
     &im_dest,
     &indexEXP,
     &level,&finest_level,
     domlo,domhi, 
     tilelo,tilehi,
     fablo,fabhi,
     &bfact, 
     xlo,dx,
     &dt_slab, //fort_distributeexpansion
     maskcov.dataPtr(),
     ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
     newdistfab.dataPtr(),
     ARLIM(newdistfab.loVect()),ARLIM(newdistfab.hiVect()),
     donorfab.dataPtr(),
     ARLIM(donorfab.loVect()),ARLIM(donorfab.hiVect()),
     donor_comp_fab.dataPtr(),
     ARLIM(donor_comp_fab.loVect()),
     ARLIM(donor_comp_fab.hiVect()),
     weightfab.dataPtr(),
     ARLIM(weightfab.loVect()),ARLIM(weightfab.hiVect()),
     weight_comp_fab.dataPtr(),
     ARLIM(weight_comp_fab.loVect()),
     ARLIM(weight_comp_fab.hiVect()),
     JUMPfab.dataPtr(),
     ARLIM(JUMPfab.loVect()),
     ARLIM(JUMPfab.hiVect()),
     JUMP_comp_fab.dataPtr(),
     ARLIM(JUMP_comp_fab.loVect()),
     ARLIM(JUMP_comp_fab.hiVect()) );
  } // mfi
} //omp
  ns_reconcile_d_num(LOOP_DISTRIBUTEEXPANSION,
     "level_phase_change_redistribute");

  for (int tid=1;tid<thread_class::nthreads;tid++) {
   mdot_sum2_local[0]+=mdot_sum2_local[tid];
   mdot_lost_local[0]+=mdot_lost_local[tid];
   mdot_sum2_complement_local[0]+=mdot_sum2_complement_local[tid];
   mdot_lost_complement_local[0]+=mdot_lost_complement_local[tid];
  } // tid
  ParallelDescriptor::ReduceRealSum(mdot_sum2_local[0]);
  ParallelDescriptor::ReduceRealSum(mdot_lost_local[0]);
  mdot_sum2[0]+=mdot_sum2_local[0];
  mdot_lost[0]+=mdot_lost_local[0];

  ParallelDescriptor::ReduceRealSum(mdot_sum2_complement_local[0]);
  ParallelDescriptor::ReduceRealSum(mdot_lost_complement_local[0]);
  mdot_sum2_complement[0]+=mdot_sum2_complement_local[0];
  mdot_lost_complement[0]+=mdot_lost_complement_local[0];

  // isweep==0: fort_tagexpansion
  // isweep==1: fort_accept_weight
  // isweep==2: fort_distributeexpansion
  // isweep==3: fort_initjumpterm
 } else if (isweep==3) {

  Vector<Real> mdotplus_local;
  Vector<Real> mdotminus_local;
  Vector<Real> mdotcount_local;
  mdotplus_local.resize(thread_class::nthreads);
  mdotminus_local.resize(thread_class::nthreads);
  mdotcount_local.resize(thread_class::nthreads);

  for (int tid=0;tid<thread_class::nthreads;tid++) {
   mdotplus_local[tid]=0.0;
   mdotminus_local[tid]=0.0;
   mdotcount_local[tid]=0.0;
  }
 
  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(localMF[MDOT_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel 
#endif
{
  for (MFIter mfi(*localMF[MDOT_MF],use_tiling); mfi.isValid(); ++mfi) {
    BL_ASSERT(grids[mfi.index()] == mfi.validbox());
    const int gridno = mfi.index();
    const Box& tilegrid = mfi.tilebox();
    const Box& fabgrid = grids[gridno];
    const int* tilelo=tilegrid.loVect();
    const int* tilehi=tilegrid.hiVect();
    const int* fablo=fabgrid.loVect();
    const int* fabhi=fabgrid.hiVect();
    const Real* xlo = grid_loc[gridno].lo();

    FArrayBox& maskcov=(*localMF[MASKCOEF_MF])[mfi];
    FArrayBox& mdotfab=(*localMF[MDOT_MF])[mfi];
    FArrayBox& JUMPfab=(*localMF[JUMP_STRENGTH_MF])[mfi];
    FArrayBox& newdistfab=(*localMF[LSNEW_MF])[mfi];
    FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi];

    int bfact=parent->Space_blockingFactor(level);

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // NavierStokes::allocate_mdot() called at the beginning of
    //  NavierStokes::do_the_advance
    // mdot initialized in NavierStokes::prelim_alloc()
    // declared in: GODUNOV_3D.F90 (distribute_from_target==0)
    //   a)  jump_strength=
    //       JUMPFAB(D_DECL(i,j,k),iten+ireverse*num_interfaces)
    //      dF * volgrid * (den_source/den_dest-1)/ dt^2 
    //   b)  divu_material=jump_strength  cm^3/s^2
    //   c)  mdot(D_DECL(i,j,k))=mdot(D_DECL(i,j,k))+divu_material
    //   V = U*  - dt grad p/rho
    //   0 = div U* - dt div grad p/rho
    //   -div U*/dt = -div grad p/rho
    //   mdot - vol div U*/dt = -vol div grad p/rho
    //   cm^3  * (1/cm) (cm/s) (1/s) = cm^3/s^2
    //   vol div V/dt = mdot     div V= dt mdot / vol  
    //   dt div V=dt^2 mdot/vol=dF * (den_source/den_dst-1)
    //   if compressible, then instead of increasing the volume, the 
    //   mass is increased instead: rho^expand - rho = -dt rho^expand div V
    //   rho=rho^expand (1+ dt div V)
    //   1+dt div V=1+dt^2 mdot/vol = 1+ dF * (den_source/den_dest-1)
    //   rho=den^dest * (1 + dF *(den_source/den_dest-1))=
    //    (1-dF)*den^dest + dF * den_source
    fort_initjumpterm( 
     &mdotplus_local[tid_current],
     &mdotminus_local[tid_current],
     &mdotcount_local[tid_current],
     &cur_time_slab,
     &level,
     &finest_level,
     saturation_temp.dataPtr(),
     freezing_model.dataPtr(),
     distribute_from_target.dataPtr(),
     constant_volume_mdot.dataPtr(),
     constant_density_all_time.dataPtr(),
     tilelo,tilehi,
     fablo,fabhi,
     &bfact, 
     xlo,dx,
     &dt_slab, //initjumpterm
     maskcov.dataPtr(),
     ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
     JUMPfab.dataPtr(),ARLIM(JUMPfab.loVect()),ARLIM(JUMPfab.hiVect()),
      // mdotfab is incremented.
     mdotfab.dataPtr(),ARLIM(mdotfab.loVect()),ARLIM(mdotfab.hiVect()),
     newdistfab.dataPtr(),
     ARLIM(newdistfab.loVect()),ARLIM(newdistfab.hiVect()),
     reconfab.dataPtr(),
     ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()));

  } // mfi
} // omp
  ns_reconcile_d_num(LOOP_INITJUMPTERM,"level_phase_change_redistribute");

  for (int tid=1;tid<thread_class::nthreads;tid++) {
   mdotplus_local[0]+=mdotplus_local[tid];
   mdotminus_local[0]+=mdotminus_local[tid];
   mdotcount_local[0]+=mdotcount_local[tid];
  }
  ParallelDescriptor::ReduceRealSum(mdotplus_local[0]);
  ParallelDescriptor::ReduceRealSum(mdotminus_local[0]);
  ParallelDescriptor::ReduceRealSum(mdotcount_local[0]);

  mdotplus[0]+=mdotplus_local[0];
  mdotminus[0]+=mdotminus_local[0];
  mdotcount[0]+=mdotcount_local[0];

 } else
  amrex::Error("isweep invalid");

} // subroutine level_phase_change_redistribute

// called from: NavierStokes::make_physics_varsALL
void
NavierStokes::level_init_icemask_and_icefacecut() {

 std::string local_caller_string="level_init_icemask_and_icefacecut";

 bool use_tiling=ns_tiling;
 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid level_init_icemask_and_icefacecut");

 resize_maskfiner(1,MASKCOEF_MF);
 VOF_Recon_resize(1); //output:SLOPE_RECON_MF

 debug_ngrow(SLOPE_RECON_MF,1,local_caller_string);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(FACE_VAR_MF+dir,0,local_caller_string);

 debug_ngrow(MASKCOEF_MF,1,local_caller_string);

 getStateDist_localMF(LSNEW_MF,1,cur_time_slab,local_caller_string);

 debug_ngrow(LSNEW_MF,1,local_caller_string);
 if (localMF[LSNEW_MF]->nComp()!=num_materials*(1+AMREX_SPACEDIM))
  amrex::Error("localMF[LSNEW_MF]->nComp() invalid");

 int nden=num_materials*num_state_material;
 MultiFab* state_var_mf=getStateDen(1,cur_time_slab);
 if (state_var_mf->nComp()!=nden)
  amrex::Error("state_var_mf->nComp()!=nden");

 const Real* dx = geom.CellSize();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(localMF[LSNEW_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel 
#endif
{
  for (MFIter mfi(*localMF[LSNEW_MF],use_tiling); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();
   const Real* xlo = grid_loc[gridno].lo();

   FArrayBox& maskcov=(*localMF[MASKCOEF_MF])[mfi];

   FArrayBox& xface=(*localMF[FACE_VAR_MF])[mfi];
   FArrayBox& yface=(*localMF[FACE_VAR_MF+1])[mfi];
   FArrayBox& zface=(*localMF[FACE_VAR_MF+AMREX_SPACEDIM-1])[mfi];
   FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi]; 
   FArrayBox& newdistfab=(*localMF[LSNEW_MF])[mfi];
   FArrayBox& denstatefab=(*state_var_mf)[mfi];

   int bfact=parent->Space_blockingFactor(level);

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();
   
    // declared in: GODUNOV_3D.F90
   fort_init_icemask_and_icefacecut( 
    rigid_fraction_id.dataPtr(),
    &nden,
    &cur_time_slab,
    &level,&finest_level,
    saturation_temp.dataPtr(),
    freezing_model.dataPtr(),
    distribute_from_target.dataPtr(),
    tilelo,tilehi,
    fablo,fabhi,
    &bfact, 
    xlo,dx,
    &dt_slab, //fort_init_icemask_and_icefacecut
    maskcov.dataPtr(),
    ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
    xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()),
    yface.dataPtr(),ARLIM(yface.loVect()),ARLIM(yface.hiVect()),
    zface.dataPtr(),ARLIM(zface.loVect()),ARLIM(zface.hiVect()),
    denstatefab.dataPtr(),
    ARLIM(denstatefab.loVect()),ARLIM(denstatefab.hiVect()),
    newdistfab.dataPtr(),
    ARLIM(newdistfab.loVect()),ARLIM(newdistfab.hiVect()),
    reconfab.dataPtr(),ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()));

  } // mfi
} // omp
  ns_reconcile_d_num(LOOP_INIT_ICEMASK_AND_ICEFACECUT,
	"level_init_icemask_and_icefacecut");

  delete_localMF(LSNEW_MF,1);
  delete state_var_mf;

} // subroutine level_init_icemask_and_icefacecut


// 1. called if "is_GFM_freezing_model"
//
// 2. multiphase_project->allocate_project_variables->stefan_solver_init
//    (project_option=SOLVETYPE_HEAT or SOLVETYPE_SPEC)
//    (adjust_temperature==1)
//    coeffMF==localMF[OUTER_ITER_PRESSURE_MF]
//
// 3. multiphase_project->allocate_maccoef->stefan_solver_init
//    (project_option=SOLVETYPE_HEAT or SOLVETYPE_SPEC)
//    (adjust_temperature==0)
//    coeffMF==localMF[ALPHANOVOLUME_MF]
//
// 4. multiphase_project->allocate_FACE_WEIGHT->stefan_solver_init
//    (project_option=SOLVETYPE_HEAT or SOLVETYPE_SPEC)
//    (adjust_temperature==-1)
//    coeffMF==localMF[CELL_DEN_MF]
//
// 5. update_SEM_forcesALL->allocate_maccoef->
//    (create_hierarchy=-1)
//
// 6. update_SEM_forcesALL->allocate_FACE_WEIGHT->
//      (face_weight_op==SUB_OP_FOR_SDC)
//
// 7. diffusion_heatingALL->allocate_FACE_WEIGHT->
//      (project_option==SOLVETYPE_VISC)
//
// if adjust_temperature==1,
//  Snew=(c1 Tn + c2 TSAT)/(c1+c2)
//  coeffMF=(c1 Tn + c2 TSAT)/(c1+c2)
//  c1=rho cv/(dt*sweptfactor)  c2=(1/vol) sum_face Aface k_m/(theta dx)
// else if adjust_temperature==0,
//  coeffMF=c1+c2
// else if adjust_temperature==-1,
//  FACECOMP_FACEHEAT component of FACE_VAR_MF (c++)
//
// for mass fraction:
// (rho Y)_t + div (rho u Y) = div rho D grad Y
// since rho_t + div (rho u)=0,
// rho (Y_t + u dot grad Y)=div rho D grad Y
// "rho D" coefficient is in FACE_VAR_MF,
//   FACECOMP_FACESPEC ... FACECOMP_FACESPEC+num_species_var-1 (c++)
// "1/rho" coefficient is in CELL_DEN_MF 

void
NavierStokes::stefan_solver_init(MultiFab* coeffMF,
		int adjust_temperature,
		int project_option) {
 
 std::string local_caller_string="stefan_solver_init";

 int finest_level=parent->finestLevel();

 bool use_tiling=ns_tiling;

 if (adjust_temperature==1) {
  // do nothing
 } else if (adjust_temperature==0) {
  // do nothing
 } else if (adjust_temperature==-1) {
  use_tiling=false;
 } else
  amrex::Error("adjust_temperature invalid");

 int nstate=STATE_NCOMP;

 int nsolve=1;

 int ntsat=EXTRAP_NCOMP_TSAT;

 if (localMF[SATURATION_TEMP_MF]->nComp()!=ntsat)
  amrex::Error("localMF[SATURATION_TEMP_MF]->nComp()!=ntsat");
 if (localMF[SATURATION_TEMP_MF]->nGrow()!=ngrow_distance)
  amrex::Error("localMF[SATURATION_TEMP_MF] incorrect ngrow");

 debug_ngrow(VOLUME_MF,1,local_caller_string);
 debug_ngrow(SWEPT_CROSSING_MF,0,local_caller_string);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(FACE_VAR_MF+dir,0,local_caller_string);

 debug_ngrow(CELL_DEN_MF,1,local_caller_string); 

 debug_ngrow(CELL_DEDT_MF,1,local_caller_string); 

 if (localMF[CELL_DEN_MF]->nComp()!=1)
  amrex::Error("localMF[CELL_DEN_MF]->nComp() invalid");

 if (localMF[CELL_DEDT_MF]->nComp()!=1)
  amrex::Error("localMF[CELL_DEDT_MF]->nComp() invalid");

 if (dt_slab>0.0) {
  // do nothing
 } else
  amrex::Error("dt_slab must be positive");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int GFM_flag=0;
 int face_comp_index=0;

 if (is_phasechange==1) {
  if (project_option==SOLVETYPE_HEAT) { 
   face_comp_index=FACECOMP_FACEHEAT;
   for (int im=0;im<2*num_interfaces;im++) {
    Real LL=get_user_latent_heat(im+1,293.0,1);
    if (LL!=0.0) {
     if (is_GFM_freezing_model(freezing_model[im])==1) {
      GFM_flag=1;
     } else if (is_GFM_freezing_model(freezing_model[im])==0) {
      // do nothing
     } else 
      amrex::Error("is_GFM_freezing_model bust");
    } else if (LL==0.0) {
     // do nothing
    } else
     amrex::Error("latent_heat[im] (LL) invalid");
   } // im=0..2 num_interfaces-1
  } else if ((project_option>=SOLVETYPE_SPEC)&&  
	     (project_option<SOLVETYPE_SPEC+num_species_var)) { //mass fraction
   face_comp_index=FACECOMP_FACESPEC+project_option-SOLVETYPE_SPEC;
   for (int im=0;im<2*num_interfaces;im++) {
    Real LL=get_user_latent_heat(im+1,293.0,1);
    if (LL!=0.0) {
     if (is_GFM_freezing_model(freezing_model[im])==1) {

      if (is_multi_component_evap(freezing_model[im],
	   Tanasawa_or_Schrage_or_Kassemi[im],
           LL)==1) {
       int ispec=mass_fraction_id[im];
       if ((ispec>=1)&&(ispec<=num_species_var)) {
        if (ispec==project_option-SOLVETYPE_SPEC+1) {
         GFM_flag=1;
        }
       } else
        amrex::Error("ispec invalid");
      } else if (is_multi_component_evap(freezing_model[im],
           Tanasawa_or_Schrage_or_Kassemi[im],
           LL)==0) {
       // do nothing
      } else 
       amrex::Error("is_multi_component_evap invalid");
     } else if (is_GFM_freezing_model(freezing_model[im])==0) {
      // do nothing
     } else 
      amrex::Error("is_GFM_freezing_model bust");
    } else if (LL==0.0) {
     // do nothing
    } else
     amrex::Error("latent_heat[im] (LL) invalid");
   } // im=0.. 2 num_interfaces -1
  } else
   amrex::Error("project_option invalid 15679");

 } else if (is_phasechange==0) {
  amrex::Error("is_phasechange invalid in stefan_solver_init (1)");
 } else
  amrex::Error("is_phasechange invalid in stefan_solver_init (2)");

 if (GFM_flag!=1)
  amrex::Error("Gibou et al algorithm only used for GFM");

 const Real* dx = geom.CellSize();

 MultiFab* LSmf=getStateDist(1,cur_time_slab,local_caller_string);
 if (LSmf->nComp()!=num_materials*(1+AMREX_SPACEDIM))
  amrex::Error("LSmf invalid ncomp");
 if (LSmf->nGrow()!=1)
  amrex::Error("LSmf->nGrow()!=1");

  // temperature and density for all of the materials.
 int nden=num_materials*num_state_material;
 MultiFab* state_var_mf=getStateDen(1,cur_time_slab);
 if (state_var_mf->nComp()!=nden)
  amrex::Error("state_var_mf->nComp()!=nden");

 MultiFab& S_new = get_new_data(State_Type,slab_step+1);
 if (S_new.nComp()!=nstate)
  amrex::Error("S_new invalid ncomp");

 int num_materials_combine=num_materials;
 int state_index;
 Vector<int> scomp;
 Vector<int> ncomp;
 int ncomp_check;

 get_mm_scomp_solver(
   num_materials_combine,
   SOLVETYPE_HEAT,
   state_index,
   scomp,ncomp,ncomp_check);

 if ((ncomp_check!=num_materials)||(state_index!=State_Type))
  amrex::Error("(ncomp_check!=num_materials)||(state_index!=State_Type)");

 MultiFab* T_list_mf=getState_list(1,scomp,ncomp,cur_time_slab);

 get_mm_scomp_solver(
   num_materials_combine,
   project_option,
   state_index,
   scomp,ncomp,ncomp_check);

 if ((ncomp_check!=num_materials)||(state_index!=State_Type))
  amrex::Error("(ncomp_check!=num_materials)||(state_index!=State_Type)");

 MultiFab* TorY_list_mf=getState_list(1,scomp,ncomp,cur_time_slab);

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,1,local_caller_string); 
 debug_ixType(MASKCOEF_MF,-1,local_caller_string);

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(LSmf->boxArray().d_numPts());
 
#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*LSmf,use_tiling); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();
   int bfact=parent->Space_blockingFactor(level);

   const Real* xlo = grid_loc[gridno].lo();

   // mask=tag if not covered by level+1 or outside the domain.
   // mask=1-tag if covered by level+1 and inside the domain.
   // NavierStokes::maskfiner  (clear_phys_boundary==0)
   FArrayBox& maskfab=(*localMF[MASKCOEF_MF])[mfi];

   FArrayBox& denstatefab=(*state_var_mf)[mfi];

   FArrayBox& lsfab=(*LSmf)[mfi];
   FArrayBox& T_fab=(*T_list_mf)[mfi];
   FArrayBox& TorY_fab=(*TorY_list_mf)[mfi];

   FArrayBox& snewfab=S_new[mfi];
   FArrayBox& DeDTfab=(*localMF[CELL_DEDT_MF])[mfi];  // 1/(rho cv)
   FArrayBox& denfab=(*localMF[CELL_DEN_MF])[mfi];  // 1/rho
    // localMF[ALPHANOVOLUME_MF] if adjust_temperature==0
    // localMF[CELL_DEN_MF] if adjust_temperature==-1
    // localMF[OUTER_ITER_PRESSURE_MF] if adjust_temperature==1
   FArrayBox& coefffab=(*coeffMF)[mfi];  

   if ((adjust_temperature==0)||
       (adjust_temperature==1)) {
    if (coefffab.nComp()!=nsolve) {
     std::cout << "coefffab.nComp()= " << coefffab.nComp() << '\n';
     std::cout << "nsolve= " << nsolve << '\n';
     std::cout << "adjust_temperature= " << adjust_temperature << '\n';
     amrex::Error("coefffab.nComp() invalid");
    }
   } else if (adjust_temperature==-1) {
    if (coefffab.nComp()!=1) {
     std::cout << "coefffab.nComp()= " << coefffab.nComp() << '\n';
     std::cout << "nsolve= " << nsolve << '\n';
     std::cout << "adjust_temperature= " << adjust_temperature << '\n';
     amrex::Error("coefffab.nComp() invalid");
    }
   } else
    amrex::Error("adjust_temperature invalid");

   FArrayBox& heatx=(*localMF[FACE_VAR_MF])[mfi];
   FArrayBox& heaty=(*localMF[FACE_VAR_MF+1])[mfi];
   FArrayBox& heatz=(*localMF[FACE_VAR_MF+AMREX_SPACEDIM-1])[mfi];
   FArrayBox& areax=(*localMF[AREA_MF])[mfi];
   FArrayBox& areay=(*localMF[AREA_MF+1])[mfi];
   FArrayBox& areaz=(*localMF[AREA_MF+AMREX_SPACEDIM-1])[mfi];
   FArrayBox& volfab=(*localMF[VOLUME_MF])[mfi];
   FArrayBox& sweptfab = (*localMF[SWEPT_CROSSING_MF])[mfi];

   FArrayBox& Tsatfab=(*localMF[SATURATION_TEMP_MF])[mfi];
   if (Tsatfab.nComp()!=ntsat) {
    std::cout << "Tsatfab.nComp()=" << Tsatfab.nComp() << 
      " ntsat=" << ntsat << '\n';
    amrex::Error("Tsatfab.nComp()!=ntsat 4");
   }

   FArrayBox& thermal_conductivity_fab=
	 (*localMF[CELL_CONDUCTIVITY_MATERIAL_MF])[mfi];
   if (thermal_conductivity_fab.nComp()!=num_materials)
    amrex::Error("thermal_conductivity_fab.nComp()!=num_materials");

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // FUTURE: 
    // 1. if an interface changing phase exists, but the raw distance
    // function does not change sign, then assume the interface is an
    // idealized point source at the center of the cell and force
    // T=TSAT there.
    // 2. The nucleation routine only seeds a tiny volume fraction in
    // "fluid" cells.
    // fort_stefansolver declared in: GODUNOV_3D.F90
   fort_stefansolver( 
     //SOLVETYPE_HEAT or SOLVETYPE_SPEC...SOLVETYPE_SPEC+num_species_var-1
    &project_option, 
    &solidheat_flag, //0=diffuse in solid 1=dirichlet 2=Neumann
    microlayer_size.dataPtr(), 
    microlayer_substrate.dataPtr(), 
    microlayer_temperature_substrate.dataPtr(), 
    &adjust_temperature,
    &nstate,
    &ntsat, // num_interfaces*(ncomp_per_tsat+1)
    &nden,  // num_materials*num_state_material
    freezing_model.dataPtr(),
    distribute_from_target.dataPtr(),
    saturation_temp.dataPtr(),
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    &level,
    &finest_level,
    xlo,dx,
    &dt_slab, //stefansolver
    maskfab.dataPtr(),
    ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
    thermal_conductivity_fab.dataPtr(), //num_materials components
    ARLIM(thermal_conductivity_fab.loVect()),
    ARLIM(thermal_conductivity_fab.hiVect()),
    denstatefab.dataPtr(),
    ARLIM(denstatefab.loVect()),ARLIM(denstatefab.hiVect()),
    Tsatfab.dataPtr(),
    ARLIM(Tsatfab.loVect()),ARLIM(Tsatfab.hiVect()),
    sweptfab.dataPtr(),ARLIM(sweptfab.loVect()),ARLIM(sweptfab.hiVect()),
    lsfab.dataPtr(),ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
    T_fab.dataPtr(),
    ARLIM(T_fab.loVect()),ARLIM(T_fab.hiVect()),
    TorY_fab.dataPtr(),
    ARLIM(TorY_fab.loVect()),ARLIM(TorY_fab.hiVect()),
    snewfab.dataPtr(),ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()),
    DeDTfab.dataPtr(),ARLIM(DeDTfab.loVect()),ARLIM(DeDTfab.hiVect()),
    denfab.dataPtr(),
    ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
    coefffab.dataPtr(),ARLIM(coefffab.loVect()),ARLIM(coefffab.hiVect()),
    volfab.dataPtr(),ARLIM(volfab.loVect()),ARLIM(volfab.hiVect()),
    heatx.dataPtr(face_comp_index),ARLIM(heatx.loVect()),ARLIM(heatx.hiVect()),
    heaty.dataPtr(face_comp_index),ARLIM(heaty.loVect()),ARLIM(heaty.hiVect()),
    heatz.dataPtr(face_comp_index),ARLIM(heatz.loVect()),ARLIM(heatz.hiVect()),
    areax.dataPtr(),ARLIM(areax.loVect()),ARLIM(areax.hiVect()),
    areay.dataPtr(),ARLIM(areay.loVect()),ARLIM(areay.hiVect()),
    areaz.dataPtr(),ARLIM(areaz.loVect()),ARLIM(areaz.hiVect()) );
 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_STEFANSOLVER,"stefan_solver_init");

 delete state_var_mf;
 delete LSmf;
 delete T_list_mf;
 delete TorY_list_mf;
 
}  // end subroutine stefan_solver_init

// MEHDI VAHAB HEAT SOURCE
// T^new=T^* += dt A Q/(rho cv V) 
// called from NavierStokes::allocate_project_variables when 
//   project_option==SOLVETYPE_HEAT
void
NavierStokes::heat_source_term_flux_source() {
 
 std::string local_caller_string="heat_source_term_flux_source";

 int finest_level=parent->finestLevel();

 bool use_tiling=ns_tiling;

 resize_metrics(1);

 debug_ngrow(VOLUME_MF,1,local_caller_string);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(FACE_VAR_MF+dir,0,local_caller_string);

 if (dt_slab>0.0) {
  //do nothing
 } else
  amrex::Error("dt_slab must be positive");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int nstate=STATE_NCOMP;

 const Real* dx = geom.CellSize();

 MultiFab* LSmf=getStateDist(1,cur_time_slab,local_caller_string); 
 if (LSmf->nComp()!=num_materials*(1+AMREX_SPACEDIM))
  amrex::Error("LSmf invalid ncomp");
 if (LSmf->nGrow()!=1)
  amrex::Error("LSmf->nGrow()!=1");

 MultiFab& S_new = get_new_data(State_Type,slab_step+1);
 if (S_new.nComp()!=nstate)
  amrex::Error("S_new invalid ncomp");

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(LSmf->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*LSmf,use_tiling); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();
   int bfact=parent->Space_blockingFactor(level);

   const Real* xlo = grid_loc[gridno].lo();
   FArrayBox& lsfab=(*LSmf)[mfi];
   FArrayBox& snewfab=S_new[mfi];
   FArrayBox& DeDTfab=(*localMF[CELL_DEDT_MF])[mfi];  // 1/(rho cv)
   if (DeDTfab.nComp()!=1)
    amrex::Error("DeDTfab.nComp() invalid");

   FArrayBox& denfab=(*localMF[CELL_DEN_MF])[mfi];  // 1/rho
   if (denfab.nComp()!=1)
    amrex::Error("denfab.nComp() invalid");

   FArrayBox& heatx=(*localMF[FACE_VAR_MF])[mfi];
   FArrayBox& heaty=(*localMF[FACE_VAR_MF+1])[mfi];
   FArrayBox& heatz=(*localMF[FACE_VAR_MF+AMREX_SPACEDIM-1])[mfi];
   FArrayBox& areax=(*localMF[AREA_MF])[mfi];
   FArrayBox& areay=(*localMF[AREA_MF+1])[mfi];
   FArrayBox& areaz=(*localMF[AREA_MF+AMREX_SPACEDIM-1])[mfi];
   FArrayBox& volfab=(*localMF[VOLUME_MF])[mfi];

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // declared in: GODUNOV_3D.F90
   fort_heatsource_face( 
    &nstate,
    saturation_temp.dataPtr(),
    tilelo,tilehi,
    fablo,fabhi,
    &bfact,
    xlo,dx,
    &dt_slab, //fort_heatsource_face
    &prev_time_slab,
    &level,
    &finest_level,
    lsfab.dataPtr(),ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
    snewfab.dataPtr(),ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()),
    DeDTfab.dataPtr(),ARLIM(DeDTfab.loVect()),ARLIM(DeDTfab.hiVect()),
    denfab.dataPtr(),
    ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
    volfab.dataPtr(),ARLIM(volfab.loVect()),ARLIM(volfab.hiVect()),
    heatx.dataPtr(FACECOMP_FACEHEAT), 
    ARLIM(heatx.loVect()),ARLIM(heatx.hiVect()),
    heaty.dataPtr(FACECOMP_FACEHEAT),
    ARLIM(heaty.loVect()),ARLIM(heaty.hiVect()),
    heatz.dataPtr(FACECOMP_FACEHEAT),
    ARLIM(heatz.loVect()),ARLIM(heatz.hiVect()),
    areax.dataPtr(),ARLIM(areax.loVect()),ARLIM(areax.hiVect()),
    areay.dataPtr(),ARLIM(areay.loVect()),ARLIM(areay.hiVect()),
    areaz.dataPtr(),ARLIM(areaz.loVect()),ARLIM(areaz.hiVect()) );
 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_HEATSOURCE_FACE,"heat_source_term_flux_source");

 delete LSmf;
 
}  // end subroutine heat_source_term_flux_source

void NavierStokes::show_norm2_id(int mf_id,int id) {

 if (show_norm2_flag==1) {
  int finest_level=parent->finestLevel();
  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   MultiFab* mf=ns_level.localMF[mf_id];
   int scomp=0;
   int ncomp=mf->nComp();
   ns_level.show_norm2(mf,scomp,ncomp,id);
  } // ilev
 } else if (show_norm2_flag==0) {
  // do nothing
 } else
  amrex::Error("show_norm2_flag invalid");
 
} // subroutine show_norm2_id

void NavierStokes::show_norm2(MultiFab* mf,int scomp,int ncomp,int id) {

 if (show_norm2_flag==1) {
  for (int nc=0;nc<ncomp;nc++) {
   Real nm2=MultiFab::Dot(*mf,scomp+nc,*mf,scomp+nc,1,0,0);
   std::cout << "\n norm2 id= " << id << "scomp+nc= " << 
    scomp+nc << " nm2= " << nm2 << '\n';
  } // nc
 } else if (show_norm2_flag==0) {
  // do nothing
 } else
  amrex::Error("show_norm2_flag invalid");
  
} // subroutine show_norm2


void NavierStokes::check_value_max(int id,int mf_id,int scomp, 
		int ncomp,int ngrow,Real max_value) {

 if (mf_check_inf_bounds==1) {
  int finest_level=parent->finestLevel();
  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   MultiFab* mf=ns_level.localMF[mf_id];
   if ((scomp>=0)&&(scomp+ncomp<=mf->nComp())&&
       (ngrow>=0)&&(ngrow<=mf->nGrow())&&
       (max_value>=0.0)) {
    ns_level.check_value_max_level(id,mf,scomp,ncomp,ngrow,max_value);
   } else
    amrex::Error("scomp,ncomp, ngrow or max_value invalid");

  } // ilev=finest_level ... level
 } else if (mf_check_inf_bounds==0) {
  // do nothing
 } else
  amrex::Error("mf_check_inf_bounds invalid");
 
} // subroutine check_value_max

void NavierStokes::check_value_max_level(int id,MultiFab* mf,
		int scomp,int ncomp,int ngrow,Real max_value) {

 int finest_level=parent->finestLevel();

 if (mf_check_inf_bounds==1) {
  for (int nc=0;nc<ncomp;nc++) {
   Real nm_max=mf->norm0(scomp+nc,ngrow);
   if (nm_max>max_value) {
    std::cout << "id= " << id << " scomp+nc= " <<
     scomp+nc << " ncomp= " << ncomp <<
     " ngrow= " << ngrow << " nm_max= " <<
     nm_max << " max_value= " <<
     max_value << '\n';
    std::cout << "level= " << level << " finest_level= " << 
	    finest_level << '\n';
    amrex::Error("nm_max>max_value");
   }
  } // nc
 } else if (mf_check_inf_bounds==0) {
  // do nothing
 } else
  amrex::Error("mf_check_inf_bounds invalid");
  
} // subroutine check_value_max_level


// datatype=0 scalar or vector
// datatype=1 tensor face
// datatype=2 tensor cell
void NavierStokes::aggressive_debug(
  int datatype,
  int force_check,
  MultiFab* mf,
  int scomp,int ncomp,
  int ngrow,
  int dir,
  Real warning_cutoff) {

 if (((verbose==0)||(verbose==1))&&(force_check==0)) {
  // do nothing
 } else if ((verbose>=2)||(force_check==1)) {
  std::fflush(NULL);
  int finest_level=parent->finestLevel();
  const Real* dx = geom.CellSize();
  int ndefined=mf->nComp();
  if (ndefined<scomp+ncomp) 
   amrex::Error("scomp,ncomp invalid in aggressive debug");
  if (mf->nGrow()<ngrow)
   amrex::Error("ngrow invalid in aggressive debug");

  if (verbose>=2) {
   std::cout << "AGGRESSIVE DEBUG scomp= " << scomp << '\n';
   std::cout << "AGGRESSIVE DEBUG ncomp= " << ncomp << '\n';
   std::cout << "AGGRESSIVE DEBUG ngrow= " << ngrow << '\n';
   std::cout << "AGGRESSIVE DEBUG dir= " << dir << '\n';
  }

  const BoxArray mfBA=mf->boxArray();
  int ngrid=mfBA.size();

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(mf->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(*mf,false); mfi.isValid(); ++mfi) {
   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = mfi.validbox();

   if (fabgrid!=mfBA[gridno])
    amrex::Error("fabgrid!=mfBA[gridno]");

   const Box& growntilegrid = growntileboxTENSOR(mf,
     datatype,ngrow,dir,tilegrid,fabgrid);

   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* growlo=growntilegrid.loVect();
   const int* growhi=growntilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();
   int bfact=parent->Space_blockingFactor(level);

   FArrayBox& mffab=(*mf)[mfi];

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // declared in GODUNOV_3D.F90
   fort_aggressive(
    &datatype,
    &warning_cutoff,
    tilelo,tilehi,
    fablo,fabhi,
    growlo,growhi,
    &bfact,
    dx,
    &scomp,
    &ncomp,
    &ndefined,
    &ngrow,
    &dir,
    &verbose,
    &force_check,
    &gridno,&ngrid,&level,&finest_level,
    mffab.dataPtr(),ARLIM(mffab.loVect()),ARLIM(mffab.hiVect()));
  } // mfi
} // omp
  ns_reconcile_d_num(LOOP_AGGRESSIVE,"aggressive_debug");

  std::fflush(NULL);
 } else
  amrex::Error("verbose or force_check invalid");

} // subroutine aggressive_debug

void
NavierStokes::synchronize_flux_register(int operation_flag,
 int spectral_loop) {

 fort_check_operation_flag_MAC(&operation_flag);
 
 if (spectral_loop+1==end_spectral_loop()) {
  delete_localMF(SEM_FLUXREG_MF,1);
 } else if (spectral_loop==0) {
  localMF[SEM_FLUXREG_MF]->FillBoundary(geom.periodicity());
 } else
  amrex::Error("spectral_loop invalid");

} // subroutine synchronize_flux_register

void 
NavierStokes::allocate_flux_register(int operation_flag) {

 int ncfluxreg=0;

  // unew^{f} = 
  //   (i) unew^{f} in non-solid regions
  //   (ii) u^{f,save} + (unew^{c}-u^{c,save})^{c->f} in spectral regions 
  //   (iii) usolid in solid regions
 if (operation_flag==OP_U_COMP_CELL_MAC_TO_MAC) {
  ncfluxreg=AMREX_SPACEDIM;
 } else if (operation_flag==OP_ISCHEME_MAC) {  // advection
  ncfluxreg=AMREX_SPACEDIM*NFLUXSEM;
 } else if (operation_flag==OP_POTGRAD_TO_MAC) { 
  ncfluxreg=AMREX_SPACEDIM; // (grad pressure_potential)_mac
 } else if (operation_flag==OP_UNEW_CELL_TO_MAC) { // ucell -> umac
  ncfluxreg=AMREX_SPACEDIM;
 } else if (operation_flag==OP_UNEW_USOL_MAC_TO_MAC) { // umac -> umac
  ncfluxreg=AMREX_SPACEDIM;
  // umac -> umac+beta F^cell->mac
 } else if (operation_flag==OP_UMAC_PLUS_VISC_CELL_TO_MAC) { 
  ncfluxreg=AMREX_SPACEDIM;
 } else if (operation_flag==OP_UGRAD_MAC) { // grad U
  ncfluxreg=AMREX_SPACEDIM*AMREX_SPACEDIM;
 } else if (operation_flag==OP_PRESGRAD_MAC) { // grad p
  ncfluxreg=AMREX_SPACEDIM;
 } else if (operation_flag==OP_UGRAD_COUPLING_MAC) { // grad U (for crossterm)
  ncfluxreg=AMREX_SPACEDIM*AMREX_SPACEDIM;
 } else
  amrex::Error("operation_flag invalid3");

 new_localMF(SEM_FLUXREG_MF,ncfluxreg,1,-1);
 setVal_localMF(SEM_FLUXREG_MF,0.0,0,ncfluxreg,1); 

} // end subroutine allocate_flux_register

int 
NavierStokes::end_spectral_loop() {

 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid end_spectral_loop");

 int local_end=1;
 int bfact=parent->Space_blockingFactor(level);

 if (enable_spectral==0) {
  local_end=1;
 } else if (enable_spectral==1) {
  if (bfact==1) {
   local_end=1;
  } else if (bfact>=2) {
   local_end=2;
  } else
   amrex::Error("bfact invalid");
 } else
  amrex::Error("enable_spectral invalid end_spectral_loop() ");

 return local_end;
 
} // end subroutine end_spectral_loop

//source_term==SUB_OP_SDC_LOW_TIME:
//  slab_step (aka "k") = 0,1,2,...,ns_time_order
//source_term==SUB_OP_SDC_ISCHEME:
//  slab_step (aka "k") = 0,1,2,...,ns_time_order-1
//SDC_outer_sweeps=0...ns_time_order-1
//
//source_term==1==SUB_OP_SDC_LOW_TIME => compute F(t^{n+k/order})
//source_term==0==SUB_OP_SDC_ISCHEME => 
// compute F(t^{n+k/order,(0)})  SUB_OP_ISCHEME_PREDICT
// compute F(t^{n+k/order,(1)})  SUB_OP_ISCHEME_CORRECT
//
// first: NavierStokes::nonlinear_advection() is called
// second: NavierStokes::SEM_advectALL is called
// This routine is called from  NavierStokes::SEM_advectALL
//
// for advect_iter=0 ... 1 (source_term==0==SUB_OP_SDC_ISCHEME)
//
//  DEN_RECON_MF, VELADVECT_MF init at time advect_time_slab
// 
//  AMRSYNC_PRES_MF+dir, 
//  CONSERVE_FLUXES_MF+dir, 
//  COARSE_FINE_FLUX_MF+dir init to 1.0E+40
//
//  for spectral_loop=0 ... 1
//  for tileloop=0..1
//  for ilev=finest_level ... level
//   call this routine.
//
//  Average down CONSERVE_FLUXES_MF
//
//  init_fluxes=0
//  for ilev=finest_level ... level
//   call this routine.
// 
// end advect_iter loop
//
// spectral_loop==0 => fine data transferred to coarse in a special way
// spectral_loop==1 => coarse fluxes interpolated to fine level.
//
void 
NavierStokes::SEM_scalar_advection(int init_fluxes,int source_term,
  int spectral_loop,int tileloop) {

 std::string local_caller_string="SEM_scalar_advection";

 int finest_level=parent->finestLevel();
 
 bool use_tiling=ns_tiling;

 int num_colors=0;
 Vector<Real> blob_array;
 blob_array.resize(1);
 int blob_array_size=blob_array.size();

 if ((SDC_outer_sweeps>=0)&&(SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid");

 if (source_term==SUB_OP_SDC_LOW_TIME) {

  if ((slab_step>=0)&&(slab_step<=ns_time_order)) {
   // do nothing
  } else
   amrex::Error("slab_step invalid");

  if (advect_iter!=SUB_OP_ISCHEME_PREDICT)
   amrex::Error("advect_iter invalid");

 } else if (source_term==SUB_OP_SDC_ISCHEME) {

  if ((slab_step>=0)&&(slab_step<ns_time_order)) {
   // do nothing
  } else
   amrex::Error("slab_step invalid");

  if ((advect_iter!=SUB_OP_ISCHEME_PREDICT)&&
      (advect_iter!=SUB_OP_ISCHEME_CORRECT))
   amrex::Error("advect_iter invalid");
  
 } else
  amrex::Error("source_term invalid");

 if ((ns_time_order==1)&&(enable_spectral!=0)) 
  amrex::Error("(ns_time_order==1)&&(enable_spectral!=0)");
 if ((ns_time_order>=2)&&(enable_spectral!=1)) 
  amrex::Error("(ns_time_order>=2)&&(enable_spectral!=1)");

 if (enable_spectral==1) {

  int nparts=im_solid_map.size();
  if ((nparts<0)||(nparts>num_materials))
   amrex::Error("nparts invalid");
  Vector<int> im_solid_map_null;
  im_solid_map_null.resize(1);

  int* im_solid_map_ptr;
  int nparts_def=nparts;
  if (nparts==0) {
   im_solid_map_ptr=im_solid_map_null.dataPtr();
   nparts_def=1;
  } else if ((nparts>=1)&&(nparts<=num_materials)) {
   im_solid_map_ptr=im_solid_map.dataPtr();
  } else
   amrex::Error("nparts invalid");

  for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) { 
   if (localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0)
    amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0");
   if (localMF[FSI_GHOST_MAC_MF+data_dir]->nComp()!=nparts_def*AMREX_SPACEDIM)
    amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nComp() invalid");
  }

  if (localMF[LEVELPC_MF]->nComp()!=num_materials*(1+AMREX_SPACEDIM))
   amrex::Error("localMF[LEVELPC_MF]->nComp()!=num_materials*(1+AMREX_SPACEDIM)");

  debug_ngrow(LEVELPC_MF,1,local_caller_string);
  debug_ngrow(delta_MF,0,local_caller_string);
  debug_ngrow(MASKCOEF_MF,1,local_caller_string);
  debug_ngrow(DEN_RECON_MF,1,local_caller_string);
  debug_ngrow(VELADVECT_MF,1,local_caller_string);

  if (localMF[DEN_RECON_MF]->nComp()!=num_materials*num_state_material)
   amrex::Error(
    "localMF[DEN_RECON_MF]->nComp()!=num_materials*num_state_material");
  if (localMF[VELADVECT_MF]->nComp()!=AMREX_SPACEDIM)
   amrex::Error("localMF[VELADVECT_MF]->nComp()!=AMREX_SPACEDIM");

  int nstate=STATE_NCOMP;

  MultiFab& S_new_test=get_new_data(State_Type,slab_step+1);
  if (S_new_test.nComp()!=nstate)
   amrex::Error("nstate invalid");

  int bfact=parent->Space_blockingFactor(level);
  int bfact_c=bfact;
  int bfact_f=bfact;
  if (level>0)
   bfact_c=parent->Space_blockingFactor(level-1);
  if (level<finest_level)
   bfact_f=parent->Space_blockingFactor(level+1);

  int project_option_visc=SOLVETYPE_VISC;

  const Real* dx = geom.CellSize();
  const Box& domain = geom.Domain();
  const int* domlo = domain.loVect(); 
  const int* domhi = domain.hiVect();

  if (init_fluxes==1) {

   int operation_flag=OP_ISCHEME_MAC;
   if (localMF[SEM_FLUXREG_MF]->nComp()!=AMREX_SPACEDIM*NFLUXSEM)
    amrex::Error("localMF[SEM_FLUXREG_MF]->nComp() invalid8");

// flux variables: average down in the tangential direction 
// to the box face, copy in
// the normal direction.  Since the blocking factor is >=2, it is
// impossible to have a grid box with size of 1 cell width.
   if ((spectral_loop==0)&&(tileloop==0)) {

    if (level<finest_level) {
     avgDown_and_Copy_localMF(
      DEN_RECON_MF,
      VELADVECT_MF,
      AMRSYNC_PRES_MF,
      operation_flag);
    } else if (level==finest_level) {
     // do nothing
    } else
     amrex::Error("level invalid21");

    if ((level>=1)&&(level<=finest_level)) {
     interp_and_Copy_localMF(
      DEN_RECON_MF,
      VELADVECT_MF,
      AMRSYNC_PRES_MF,
      operation_flag);
    } else if (level==0) {
     // do nothing
    } else
     amrex::Error("level invalid22");

   } else if ((spectral_loop==1)&&(tileloop==0)) {

     // interpolate CONSERVE_FLUXES from the coarse level to the fine level:
     // COARSE_FINE_FLUX
     // Since spectral_loop==1 and the fluxes of interest are not "shared",
     // it means that CONSERVE_FLUXES has already been init during the
     // spectral_loop==0 sweep.
    if ((level>=1)&&(level<=finest_level)) {
     interp_flux_localMF(
      CONSERVE_FLUXES_MF,
      COARSE_FINE_FLUX_MF);
    } else if (level==0) {
     // do nothing
    } else
     amrex::Error("level invalid22");

   } else if ((spectral_loop==0)&&(tileloop==1)) {
    // do nothing
   } else if ((spectral_loop==1)&&(tileloop==1)) {
    // do nothing
   } else
    amrex::Error("spectral_loop or tileloop invalid");

   for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

    debug_ngrow(UMAC_MF+dir,0,local_caller_string);

    if (thread_class::nthreads<1)
     amrex::Error("thread_class::nthreads invalid");
    thread_class::init_d_numPts(localMF[DEN_RECON_MF]->boxArray().d_numPts());
  
#ifdef _OPENMP
#pragma omp parallel
#endif
{
    for (MFIter mfi(*localMF[DEN_RECON_MF],use_tiling); 
         mfi.isValid(); ++mfi) {
     BL_ASSERT(grids[mfi.index()] == mfi.validbox());
     int gridno=mfi.index();
     const Box& tilegrid = mfi.tilebox();
     const Box& fabgrid = grids[gridno];
     const int* tilelo=tilegrid.loVect();
     const int* tilehi=tilegrid.hiVect();
     const int* fablo=fabgrid.loVect();
     const int* fabhi=fabgrid.hiVect();

     const Real* xlo = grid_loc[gridno].lo();

     FArrayBox& xface=(*localMF[CONSERVE_FLUXES_MF+dir])[mfi];  

     FArrayBox& xp=(*localMF[AMRSYNC_PRES_MF+dir])[mfi];  
     if (xp.nComp()==NFLUXSEM) {
      // do nothing
     } else
      amrex::Error("xp.nComp invalid");

     FArrayBox& xgp=(*localMF[COARSE_FINE_FLUX_MF+dir])[mfi];  

     FArrayBox& xvel=(*localMF[UMAC_MF+dir])[mfi];  

     FArrayBox& velfab=(*localMF[VELADVECT_MF])[mfi];
     FArrayBox& denfab=(*localMF[DEN_RECON_MF])[mfi];

     FArrayBox& solfab=(*localMF[FSI_GHOST_MAC_MF+dir])[mfi];
     FArrayBox& levelpcfab=(*localMF[LEVELPC_MF])[mfi];

     FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi];

     // mask=tag if not covered by level+1 or outside the domain.
     FArrayBox& maskcoeffab=(*localMF[MASKCOEF_MF])[mfi];

     FArrayBox& semfluxfab=(*localMF[SEM_FLUXREG_MF])[mfi];
     int ncfluxreg=semfluxfab.nComp();

     FArrayBox& maskSEMfab=(*localMF[MASKSEM_MF])[mfi];

     Vector<int> velbc=getBCArray(State_Type,gridno,
        STATECOMP_VEL,STATE_NCOMP_VEL);
     Vector<int> denbc=getBCArray(State_Type,gridno,STATECOMP_STATES,
      num_materials*num_state_material);

     int energyflag=SUB_OP_DEFAULT;
     int local_enable_spectral=enable_spectral;
     int simple_AMR_BC_flag=0;
     int ncomp_xp=NFLUXSEM;
     int ncomp_xgp=NFLUXSEM;
     int ncomp_mgoni=denfab.nComp();

     int tid_current=ns_thread();
     if ((tid_current<0)||(tid_current>=thread_class::nthreads))
      amrex::Error("tid_current invalid");
     thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     // in: SEM_scalar_advection
     // advect: rho u, rho, temperature 
     // fort_cell_to_mac in LEVELSET_3D.F90
     int nsolve=1; //unused here
     int ncphys_proxy=NFLUXSEM;

     fort_cell_to_mac(
      &ncomp_mgoni,
      &ncomp_xp,
      &ncomp_xgp,
      &simple_AMR_BC_flag,
      &nsolve,
      &tileloop,
      &dir,
      &operation_flag, // OP_ISCHEME_MAC
      &energyflag,
      &visc_coef, //beta
      &visc_coef,
      &local_enable_spectral,
      &ncphys_proxy, // ncphys (NFLUXSEM)
      constant_density_all_time.dataPtr(),
      denbc.dataPtr(),  // presbc
      velbc.dataPtr(),  
      &slab_step,
      &dt_slab,  // CELL_TO_MAC (OP_ISCEME_MAC)
      &prev_time_slab, 
      xlo,dx,
      &spectral_loop,
      &ncfluxreg,
      semfluxfab.dataPtr(),
      ARLIM(semfluxfab.loVect()),ARLIM(semfluxfab.hiVect()),
      maskfab.dataPtr(), // mask=1.0 fine/fine   0.0=coarse/fine
      ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
      maskcoeffab.dataPtr(), // maskcoef: 1=not covered 0=covered
      ARLIM(maskcoeffab.loVect()),ARLIM(maskcoeffab.hiVect()),
      maskSEMfab.dataPtr(),
      ARLIM(maskSEMfab.loVect()),ARLIM(maskSEMfab.hiVect()),
      levelpcfab.dataPtr(),
      ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),
      solfab.dataPtr(),
      ARLIM(solfab.loVect()),ARLIM(solfab.hiVect()),
      xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()), //xcut
      xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()), // xflux
      xgp.dataPtr(),ARLIM(xgp.loVect()),ARLIM(xgp.hiVect()),//holds COARSE_FINE
      xp.dataPtr(),ARLIM(xp.loVect()),ARLIM(xp.hiVect()),//holds AMRSYNC_PRES
      xvel.dataPtr(),ARLIM(xvel.loVect()),ARLIM(xvel.hiVect()), 
      velfab.dataPtr(),
      ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
      denfab.dataPtr(),  // pres
      ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
      denfab.dataPtr(),  // den
      ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
      denfab.dataPtr(),  // mgoni 
      ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
      denfab.dataPtr(),  // color
      ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
      denfab.dataPtr(),  // type
      ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
      tilelo,tilehi,
      fablo,fabhi,
      &bfact,&bfact_c,&bfact_f, 
      &level,&finest_level,
      &NS_geometry_coord,
      domlo,domhi, 
      &nparts,
      &nparts_def,
      im_solid_map_ptr,
      blob_array.dataPtr(),
      &blob_array_size,
      &num_colors,
      &project_option_visc);
    }   // mfi
} // omp
    ns_reconcile_d_num(LOOP_CELL_TO_MAC_ISCHEME_MAC,"SEM_scalar_advection");
   } // dir=0..sdim-1

  } else if (init_fluxes==0) {

   if (bfact>=1) {

    debug_ngrow(UMAC_MF,0,local_caller_string);
    debug_ngrow(UMAC_MF+1,0,local_caller_string);
    debug_ngrow(UMAC_MF+AMREX_SPACEDIM-1,0,local_caller_string);

    MultiFab& S_new=get_new_data(State_Type,slab_step+1);
    MultiFab& S_old=get_new_data(State_Type,slab_step);

    MultiFab* rhs=new MultiFab(grids,dmap,NFLUXSEM,0,
		  MFInfo().SetTag("rhs"),FArrayBoxFactory());

    if (thread_class::nthreads<1)
     amrex::Error("thread_class::nthreads invalid");
    thread_class::init_d_numPts(S_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
    for (MFIter mfi(S_new,use_tiling); mfi.isValid(); ++mfi) {
     BL_ASSERT(grids[mfi.index()] == mfi.validbox());
     const int gridno = mfi.index();
     const Box& tilegrid = mfi.tilebox();
     const Box& fabgrid = grids[gridno];
     const int* tilelo=tilegrid.loVect();
     const int* tilehi=tilegrid.hiVect();
     const int* fablo=fabgrid.loVect();
     const int* fabhi=fabgrid.hiVect();

     const Real* xlo = grid_loc[gridno].lo();

     FArrayBox& ax = (*localMF[AREA_MF])[mfi];
     FArrayBox& ay = (*localMF[AREA_MF+1])[mfi];
     FArrayBox& az = (*localMF[AREA_MF+AMREX_SPACEDIM-1])[mfi];

      // in: SEM_scalar_advection
     FArrayBox& xface=(*localMF[CONSERVE_FLUXES_MF])[mfi];  
     FArrayBox& yface=(*localMF[CONSERVE_FLUXES_MF+1])[mfi];  
     FArrayBox& zface=(*localMF[CONSERVE_FLUXES_MF+AMREX_SPACEDIM-1])[mfi];  

     FArrayBox& xvel=(*localMF[UMAC_MF])[mfi];  
     FArrayBox& yvel=(*localMF[UMAC_MF+1])[mfi];  
     FArrayBox& zvel=(*localMF[UMAC_MF+AMREX_SPACEDIM-1])[mfi];  

     FArrayBox& vol = (*localMF[VOLUME_MF])[mfi];

     // mask=1 if fine/fine
     FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi];
     // mask=tag if not covered by level+1 or outside the domain.
     FArrayBox& maskcoeffab=(*localMF[MASKCOEF_MF])[mfi];

     FArrayBox& maskSEMfab = (*localMF[MASKSEM_MF])[mfi];
     FArrayBox& rhsfab = (*rhs)[mfi];
     FArrayBox& snewfab = S_new[mfi];
     FArrayBox& S_old_fab = S_old[mfi];

     FArrayBox& velfab=(*localMF[VELADVECT_MF])[mfi];
     FArrayBox& denfab=(*localMF[DEN_RECON_MF])[mfi];

     FArrayBox& solxfab=(*localMF[FSI_GHOST_MAC_MF])[mfi];
     FArrayBox& solyfab=(*localMF[FSI_GHOST_MAC_MF+1])[mfi];
     FArrayBox& solzfab=(*localMF[FSI_GHOST_MAC_MF+AMREX_SPACEDIM-1])[mfi];

     FArrayBox& levelpcfab=(*localMF[LEVELPC_MF])[mfi];

     FArrayBox& deltafab=(*localMF[delta_MF])[mfi];
     int deltacomp=0;
     if (source_term==SUB_OP_SDC_LOW_TIME) {
      deltacomp=0;
     } else if (source_term==SUB_OP_SDC_ISCHEME) {
      if ((slab_step>=0)&&(slab_step<ns_time_order)) {
       deltacomp=slab_step*NSTATE_SDC+SEMADVECT;
      } else
       amrex::Error("slab_step invalid");
     } else
      amrex::Error("source_term invalid");

     Vector<int> velbc=getBCArray(State_Type,gridno,
        STATECOMP_VEL,STATE_NCOMP_VEL);
     Vector<int> denbc=getBCArray(State_Type,gridno,STATECOMP_STATES,
      num_materials*num_state_material);

     int operation_flag=OP_ISCHEME_CELL; // advection
     int energyflag=advect_iter;
     int nsolve=NFLUXSEM;
     int homflag=source_term;
     int local_enable_spectral=enable_spectral;

     int ncomp_denold=num_materials*num_state_material;
     int ncomp_veldest=snewfab.nComp();
     int ncomp_dendest=snewfab.nComp()-STATECOMP_STATES;

     int tid_current=ns_thread();
     if ((tid_current<0)||(tid_current>=thread_class::nthreads))
      amrex::Error("tid_current invalid");
     thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     int ncphys_proxy=NFLUXSEM;

      // called from: NavierStokes::SEM_scalar_advection
     fort_mac_to_cell(
      &ns_time_order,
      &divu_outer_sweeps,
      &num_divu_outer_sweeps,
      &operation_flag, // 107=OP_ISCHEME_CELL=advection
      &energyflag,
      constant_density_all_time.dataPtr(),
      &nparts,
      &nparts_def,
      im_solid_map_ptr,
      &level, 
      &finest_level,
      &project_option_visc,
      &local_enable_spectral,
      &ncphys_proxy, // ncphys (NFLUXSEM)
      velbc.dataPtr(),
      denbc.dataPtr(),  // presbc
      &prev_time_slab, 
      &slab_step,
      &dt_slab,  // calling fort_mac_to_cell (OP_ISCHEME_CELL)
      xlo,dx,
      tilelo,tilehi,
      fablo,fabhi,
      &bfact,
      xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()), //xp
      yface.dataPtr(),ARLIM(yface.loVect()),ARLIM(yface.hiVect()), //yp
      zface.dataPtr(),ARLIM(zface.loVect()),ARLIM(zface.hiVect()), //zp
      xvel.dataPtr(),ARLIM(xvel.loVect()),ARLIM(xvel.hiVect()), //xvel
      yvel.dataPtr(),ARLIM(yvel.loVect()),ARLIM(yvel.hiVect()), //xvel
      zvel.dataPtr(),ARLIM(zvel.loVect()),ARLIM(zvel.hiVect()), //xvel
      xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()), //xflux
      yface.dataPtr(),ARLIM(yface.loVect()),ARLIM(yface.hiVect()), //yflux
      zface.dataPtr(),ARLIM(zface.loVect()),ARLIM(zface.hiVect()), //zflux
      ax.dataPtr(),ARLIM(ax.loVect()),ARLIM(ax.hiVect()),
      ay.dataPtr(),ARLIM(ay.loVect()),ARLIM(ay.hiVect()),
      az.dataPtr(),ARLIM(az.loVect()),ARLIM(az.hiVect()),
      vol.dataPtr(),ARLIM(vol.loVect()),ARLIM(vol.hiVect()),
      rhsfab.dataPtr(),ARLIM(rhsfab.loVect()),ARLIM(rhsfab.hiVect()),
      snewfab.dataPtr(), // veldest
      ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()),
      snewfab.dataPtr(STATECOMP_STATES), // dendest
      ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()),
      maskfab.dataPtr(), // mask=1.0 fine/fine   0.0=coarse/fine
      ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
      maskcoeffab.dataPtr(), // maskcoef: 1=not covered 0=covered
      ARLIM(maskcoeffab.loVect()),ARLIM(maskcoeffab.hiVect()),
      maskSEMfab.dataPtr(), 
      ARLIM(maskSEMfab.loVect()),ARLIM(maskSEMfab.hiVect()),
      levelpcfab.dataPtr(),
      ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),
      solxfab.dataPtr(),
      ARLIM(solxfab.loVect()),ARLIM(solxfab.hiVect()),
      solyfab.dataPtr(),
      ARLIM(solyfab.loVect()),ARLIM(solyfab.hiVect()),
      solzfab.dataPtr(),
      ARLIM(solzfab.loVect()),ARLIM(solzfab.hiVect()),
      deltafab.dataPtr(deltacomp), // cterm
      ARLIM(deltafab.loVect()),ARLIM(deltafab.hiVect()),
      rhsfab.dataPtr(), // pold
      ARLIM(rhsfab.loVect()),ARLIM(rhsfab.hiVect()),
      S_old_fab.dataPtr(STATECOMP_STATES), // denold
      ARLIM(S_old_fab.loVect()),ARLIM(S_old_fab.hiVect()),
      S_old_fab.dataPtr(), // ustar
      ARLIM(S_old_fab.loVect()),ARLIM(S_old_fab.hiVect()),
      velfab.dataPtr(), // OP_ISCHEME_CELL, mdot
      ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
      denfab.dataPtr(), // OP_ISCHEME_CELL, maskdivres
      ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
      S_old_fab.dataPtr(), // maskres
      ARLIM(S_old_fab.loVect()),ARLIM(S_old_fab.hiVect()),
      &SDC_outer_sweeps,
      &homflag,
      &nsolve,
      &ncomp_denold,
      &ncomp_veldest,
      &ncomp_dendest);

     if (1==0) {
      std::cout << "SEM_scalar_advect c++ level,finest_level " << 
       level << ' ' << finest_level << '\n';
      int interior_only=1;
      tecplot_debug(snewfab,xlo,fablo,fabhi,dx,-1,0,0,
       AMREX_SPACEDIM,interior_only);
     }

    } // mfi
}// omp
    ns_reconcile_d_num(LOOP_MAC_TO_CELL_ISCHEME_CELL,"SEM_scalar_advection");

    // rhs=div(uF)  
    if (source_term==SUB_OP_SDC_ISCHEME) {

     if ((slab_step>=0)&&(slab_step<ns_time_order)) {

      if ((advect_iter==SUB_OP_ISCHEME_CORRECT)&&
          (divu_outer_sweeps+1==num_divu_outer_sweeps)) {
       int deltacomp=slab_step*NSTATE_SDC+SEMADVECT;
        // dest,src,srccomp,dstcomp,ncomp,ngrow
       MultiFab::Copy(*localMF[stableF_MF],*rhs,0,deltacomp,NFLUXSEM,0);
      } else if ((advect_iter==SUB_OP_ISCHEME_CORRECT)&&
                 (divu_outer_sweeps+1>=1)&&
                 (divu_outer_sweeps+1<num_divu_outer_sweeps)) {
       // do nothing
      } else if (advect_iter==SUB_OP_ISCHEME_PREDICT) {
       // do nothing
      } else
       amrex::Error("advect_iter invalid");

     } else
      amrex::Error("slab_step invalid");

    } else if (source_term==SUB_OP_SDC_LOW_TIME) {

     if ((slab_step==0)&&(SDC_outer_sweeps==0)) {
      // this is ok: F_advect(t^n)
     } else if ((slab_step==0)&&(SDC_outer_sweeps!=0)) {
      amrex::Error("F_advect(t^n) already init when SDC_outer_sweeps==0");
     } 

     if ((slab_step>=0)&&(slab_step<=ns_time_order)) {
      int deltacomp=slab_step*NSTATE_SDC+SEMADVECT;
      MultiFab::Copy(*localMF[spectralF_MF],*rhs,0,deltacomp,NFLUXSEM,0);
     } else
      amrex::Error("slab_step invalid");

    } else 
     amrex::Error("source_term invalid");
 
    delete rhs;

   } else
    amrex::Error("bfact invalid");

  } else
   amrex::Error("init_fluxes invalid");

 } else
  amrex::Error("enable_specral invalid");

} // subroutine SEM_scalar_advection

// Lagrangian solid info lives at t=t^n 
// order_direct_split=base_step mod 2
// must go from finest level to coarsest.
void 
NavierStokes::split_scalar_advection() { 
 
 std::string local_caller_string="split_scalar_advection";

 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();
 int normdir_here=normdir_direct_split[dir_absolute_direct_split];
 if ((normdir_here<0)||(normdir_here>=AMREX_SPACEDIM))
  amrex::Error("normdir_here invalid");

 int bfact=parent->Space_blockingFactor(level);
 int bfact_f=bfact;
 if (level<finest_level)
  bfact_f=parent->Space_blockingFactor(level+1);

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 if ((order_direct_split!=0)&&
     (order_direct_split!=1))
  amrex::Error("order_direct_split invalid");

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else {
  std::cout << "SDC_outer_sweeps= " << SDC_outer_sweeps << '\n';
  amrex::Error("SDC_outer_sweeps invalid");
 }

 if ((dir_absolute_direct_split<0)||
     (dir_absolute_direct_split>=AMREX_SPACEDIM))
  amrex::Error("dir_absolute_direct_split invalid");

  // 2 ghost cells needed in order to form (rho u)^{mac,old} in the
  // mac grid ghost cell.
 int ngrow_mass=2;
 int ngrow_scalar=1;

 int ngrow_mac_old=2;

 if (NUM_CELL_ELASTIC==num_materials_viscoelastic*ENUM_NUM_TENSOR_TYPE) {
  // do nothing
 } else
  amrex::Error("NUM_CELL_ELASTIC invalid");

  // vof,ref centroid,order,slope,intercept  x num_materials
 VOF_Recon_resize(ngrow_mass); //output:SLOPE_RECON_MF
 debug_ngrow(SLOPE_RECON_MF,ngrow_mass,local_caller_string);
 resize_maskfiner(ngrow_mass,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,ngrow_mass,local_caller_string);
 debug_ngrow(VOF_PREV_TIME_MF,ngrow_scalar,local_caller_string);
 if (localMF[VOF_PREV_TIME_MF]->nComp()!=num_materials)
  amrex::Error("vof prev time invalid ncomp");

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 int ncomp_state=S_new.nComp();
 if (ncomp_state!=STATECOMP_STATES+
     num_materials*(num_state_material+ngeom_raw)+1)
  amrex::Error("ncomp_state invalid");

 MultiFab& LS_new=get_new_data(LS_Type,slab_step+1);
 if (LS_new.nComp()!=num_materials*(AMREX_SPACEDIM+1))
  amrex::Error("LS_new ncomp invalid");

 const Real* dx = geom.CellSize();
 const Box& domain = geom.Domain();
 const int* domlo = domain.loVect();
 const int* domhi = domain.hiVect();

 Vector<int> dombc(2*AMREX_SPACEDIM);
 const BCRec& descbc = get_desc_lst()[State_Type].getBC(STATECOMP_MOF);
 const int* b_rec=descbc.vect();
 for (int m=0;m<2*AMREX_SPACEDIM;m++)
  dombc[m]=b_rec[m];

  // in: split_scalar_advection
 getStateDen_localMF(DEN_RECON_MF,ngrow_mass,advect_time_slab);

 //  (rho Y)_t + div(rho u Y)=div(D rho grad Y)

 // getStateMOM_DEN declared in: NavierStokes.cpp
 // if override_density(im)==1,
 // rho_im=rho(z)+rho0* DrhoDT * (T_im - T0_im)
 // if override_density(im)=0 or 2, density is not modified:
 // Du/Dt=-grad (p-rho0 g dot z)/rho0 - g DrhoDT (T-T0)
 getStateMOM_DEN(MOM_DEN_MF,ngrow_mass,advect_time_slab);

 int TENSOR_RECON_MF_local=-1;
 int Tensor_Type_local=-1;

 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=num_materials)) {
  TENSOR_RECON_MF_local=TENSOR_RECON_MF;
  Tensor_Type_local=Tensor_Type;
   //ngrow_scalar=1
  getStateTensor_localMF(
   TENSOR_RECON_MF_local,
   ngrow_scalar,0,
   NUM_CELL_ELASTIC,
   advect_time_slab);
 } else if (num_materials_viscoelastic==0) {
  Tensor_Type_local=State_Type;
  TENSOR_RECON_MF_local=DEN_RECON_MF;
 } else
  amrex::Error("num_materials_viscoelastic invalid");

 MultiFab& Tensor_new=get_new_data(Tensor_Type_local,slab_step+1);

   //ngrow_scalar=1
 getStateDist_localMF(LS_RECON_MF,ngrow_scalar,advect_time_slab,
     local_caller_string);

   // the pressure from before will be copied to the new pressure.
 getState_localMF(VELADVECT_MF,ngrow_mass,
  STATECOMP_VEL,
  STATE_NCOMP_VEL+STATE_NCOMP_PRES,
  advect_time_slab); 

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  getStateMAC_localMF(UMACOLD_MF+dir,ngrow_mac_old,dir,advect_time_slab);
 } // dir = 0..sdim-1

 if ((dir_absolute_direct_split<0)||
     (dir_absolute_direct_split>=AMREX_SPACEDIM))
  amrex::Error("dir_absolute_direct_split invalid");

 if (dir_absolute_direct_split==0) {

  // do nothing

 } else if ((dir_absolute_direct_split>=1)&&
            (dir_absolute_direct_split<AMREX_SPACEDIM)) {
  // do nothing
 } else
  amrex::Error("dir_absolute_direct_split invalid");

 int vofrecon_ncomp=localMF[SLOPE_RECON_MF]->nComp();
 if (vofrecon_ncomp!=num_materials*ngeom_recon)
   amrex::Error("recon ncomp bust");

 int den_recon_ncomp=localMF[DEN_RECON_MF]->nComp();
 if (den_recon_ncomp!=num_state_material*num_materials)
   amrex::Error("den_recon invalid");
 if (localMF[MOM_DEN_MF]->nComp()!=num_materials)
  amrex::Error("MOM_DEN_MF invalid nComp()");

 int LS_recon_ncomp=localMF[LS_RECON_MF]->nComp();
 if (LS_recon_ncomp!=num_materials*(1+AMREX_SPACEDIM))
   amrex::Error("LS_recon invalid");

 debug_ngrow(LS_RECON_MF,ngrow_scalar,local_caller_string);

 resize_mask_nbr(ngrow_mass);
 debug_ngrow(MASK_NBR_MF,ngrow_mass,local_caller_string); 

 MultiFab* umac_new[AMREX_SPACEDIM];
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  umac_new[dir]=&get_new_data(Umac_Type+dir,slab_step+1);
 }

 int ngrid=grids.size();

 int nc_conserve=CISLCOMP_CONS_NCOMP;
 MultiFab* conserve=new MultiFab(grids,dmap,nc_conserve,
    ngrow_mass,
    MFInfo().SetTag("conserve"),FArrayBoxFactory());

 int nc_bucket=CISLCOMP_NCOMP;

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(localMF[SLOPE_RECON_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*localMF[SLOPE_RECON_MF],use_tiling); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();

  FArrayBox& consfab=(*conserve)[mfi];
  FArrayBox& denfab=(*localMF[DEN_RECON_MF])[mfi];
  FArrayBox& mom_denfab=(*localMF[MOM_DEN_MF])[mfi];
  FArrayBox& velfab=(*localMF[VELADVECT_MF])[mfi];

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // declared in: GODUNOV_3D.F90
  fort_build_conserve( 
   constant_density_all_time.dataPtr(),
   consfab.dataPtr(),ARLIM(consfab.loVect()),ARLIM(consfab.hiVect()),
   denfab.dataPtr(),
   ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
   mom_denfab.dataPtr(),
   ARLIM(mom_denfab.loVect()),ARLIM(mom_denfab.hiVect()),
   velfab.dataPtr(),ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
   tilelo,tilehi,
   fablo,fabhi,&bfact,
   &ngrow_mass,
   &normdir_here,
   &nc_conserve,
   &den_recon_ncomp);
 }  // mfi
} // omp
 ns_reconcile_d_num(LOOP_BUILD_CONSERVE,"split_scalar_advection");

 MultiFab* xvel[AMREX_SPACEDIM]; 
 MultiFab* side_bucket_mom[AMREX_SPACEDIM]; // 2 components
 MultiFab* side_bucket_mass[AMREX_SPACEDIM]; // 2 components

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

   //ncomp=1
  xvel[dir]=new MultiFab(state[Umac_Type+dir].boxArray(),dmap,
    1,ngrow_mac_old,MFInfo().SetTag("xvel"),FArrayBoxFactory());

    //ncomp=2 ngrow=1
  side_bucket_mom[dir]=new MultiFab(grids,dmap,
     2,1,MFInfo().SetTag("side_bucket_mom"),FArrayBoxFactory());
    //scomp=0 ncomp=2 ngrow=1
  side_bucket_mom[dir]->setVal(0.0,0,2,1);

   //ncomp=2 ngrow=1
  side_bucket_mass[dir]=new MultiFab(grids,dmap,
     2,1,MFInfo().SetTag("side_bucket_mass"),FArrayBoxFactory());
   //scomp=0 ncomp=2 ngrow=1
  side_bucket_mass[dir]->setVal(0.0,0,2,1);
 }  // dir = 0..sdim-1

 for (int veldir=1;veldir<=AMREX_SPACEDIM;veldir++) {

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(localMF[SLOPE_RECON_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(*localMF[SLOPE_RECON_MF],use_tiling); 
       mfi.isValid(); ++mfi) {
    BL_ASSERT(grids[mfi.index()] == mfi.validbox());
    const int gridno = mfi.index();
    const Box& tilegrid = mfi.tilebox();
    const Box& fabgrid = grids[gridno];
    const int* tilelo=tilegrid.loVect();
    const int* tilehi=tilegrid.hiVect();
    const int* fablo=fabgrid.loVect();
    const int* fabhi=fabgrid.hiVect();

    const Real* xlo = grid_loc[gridno].lo();

    FArrayBox& xmac_old=(*localMF[UMACOLD_MF+veldir-1])[mfi];
    FArrayBox& xvelfab=(*xvel[veldir-1])[mfi]; 

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     // declared in GODUNOV_3D.F90
    fort_build_macvof( 
     &level,
     &finest_level,
     &normdir_here,
     xmac_old.dataPtr(),ARLIM(xmac_old.loVect()),ARLIM(xmac_old.hiVect()),
     xvelfab.dataPtr(),ARLIM(xvelfab.loVect()),ARLIM(xvelfab.hiVect()),
     xlo,dx,
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,
     &ngrow_mass,  //=2
     &ngrow_mac_old, //=2
     &veldir);
  }  // mfi
}// omp
  ns_reconcile_d_num(LOOP_BUILD_MACVOF,"split_scalar_advection");

 } // veldir=1..sdim

 Vector<int> nprocessed;
 nprocessed.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  nprocessed[tid]=0.0;
 }


 double profile_time_start=0.0;
 if (profile_debug==1) {
  profile_time_start=ParallelDescriptor::second();
 }

 // in: split_scalar_advection
 // initialize selected state variables 

 S_new.setVal(0.0,STATECOMP_VEL,STATE_NCOMP_VEL,1);

 if (divu_outer_sweeps==0) {
  S_new.setVal(0.0,STATECOMP_PRES,STATE_NCOMP_PRES,1);
 } else if ((divu_outer_sweeps>=1)&&
            (divu_outer_sweeps<num_divu_outer_sweeps)) {
  //do nothing
 } else
  amrex::Error("divu_outer_sweeps invalid");

 for (int im=0;im<num_materials;im++) {
  if (ns_is_rigid(im)==0) {
   S_new.setVal(0.0,
      STATECOMP_STATES+im*num_state_material,
      num_state_material,1);
   S_new.setVal(0.0,
      STATECOMP_MOF+im*ngeom_raw,
      ngeom_raw,1);
   LS_new.setVal(0.0,im,1,1);
  } else if (ns_is_rigid(im)==1) {
   if (solidheat_flag==0) {  // thermal diffuse in solid (default)
    S_new.setVal(0.0,
     STATECOMP_STATES+im*num_state_material+ENUM_TEMPERATUREVAR,1,1);
   } else if (solidheat_flag==2) { // Neumann
    // do nothing
   } else if (solidheat_flag==1) { // dirichlet
    // do nothing
   } else
    amrex::Error("solidheat_flag invalid");
  } else
   amrex::Error("ns_is_rigid(im) invalid");
 } // im=0..num_materials-1

 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=num_materials)) {

  if (NUM_CELL_ELASTIC==num_materials_viscoelastic*ENUM_NUM_TENSOR_TYPE) {
   // do nothing
  } else
   amrex::Error("NUM_CELL_ELASTIC invalid");

  if (Tensor_new.nComp()==NUM_CELL_ELASTIC) {
   // do nothing
  } else
   amrex::Error("(Tensor_new.nComp()==NUM_CELL_ELASTIC) failed");

  Tensor_new.setVal(0.0,0,NUM_CELL_ELASTIC,1);
 } else if (num_materials_viscoelastic==0) {
  // do nothing
 } else
  amrex::Error("num_materials_viscoelastic invalid:split_scalar_advection");

 if (dir_absolute_direct_split==0) {

   // initialize the error indicator to be 0.0
  S_new.setVal(0.0,ncomp_state-1,1,1);

 } else if ((dir_absolute_direct_split==1)||
            (dir_absolute_direct_split==AMREX_SPACEDIM-1)) {
  // do nothing
 } else {
  amrex::Error("dir_absolute_direct_split invalid");
 }

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(S_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(S_new,use_tiling); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());

  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();

  const Real* xlo = grid_loc[gridno].lo();

    // mask=tag if not covered by level+1 or outside the domain.
    // mask=1-tag if covered by level+1 and inside the domain.
    // NavierStokes::maskfiner  (clear_phys_boundary==0)
  FArrayBox& maskfab=(*localMF[MASKCOEF_MF])[mfi];
   // mask_nbr:
   // (1) =1 interior  =1 fine-fine ghost in domain  =0 otherwise
   // (2) =1 interior  =0 otherwise
   // (3) =1 interior+ngrow-1  =0 otherwise
   // (4) =1 interior+ngrow    =0 otherwise
  FArrayBox& masknbrfab=(*localMF[MASK_NBR_MF])[mfi];

   // velocity * dt
  FArrayBox& umac_displace=(*localMF[MAC_VELOCITY_MF+normdir_here])[mfi];
  if (umac_displace.nComp()!=1)
   amrex::Error("umac_displace has invalid ncomp");

    // this is the original data
  FArrayBox& LSfab=(*localMF[LS_RECON_MF])[mfi];
  FArrayBox& denfab=(*localMF[DEN_RECON_MF])[mfi];
  FArrayBox& mom_denfab=(*localMF[MOM_DEN_MF])[mfi];
  FArrayBox& tenfab=(*localMF[TENSOR_RECON_MF_local])[mfi];
  FArrayBox& velfab=(*localMF[VELADVECT_MF])[mfi];

    // this is the slope data
  FArrayBox& vofslopefab=(*localMF[SLOPE_RECON_MF])[mfi];

  FArrayBox& vof0fab=(*localMF[VOF_PREV_TIME_MF])[mfi];

  Vector<int> velbc=getBCArray(State_Type,gridno,normdir_here,1);

     // this is the result
  FArrayBox& destfab=S_new[mfi];
  FArrayBox& tennewfab=Tensor_new[mfi];
  FArrayBox& LSdestfab=LS_new[mfi];

  FArrayBox& consfab=(*conserve)[mfi];

  FArrayBox& xvelfab=(*xvel[0])[mfi]; 
  FArrayBox& yvelfab=(*xvel[1])[mfi];
  FArrayBox& zvelfab=(*xvel[AMREX_SPACEDIM-1])[mfi];

  FArrayBox& xmomside=(*side_bucket_mom[0])[mfi];
  FArrayBox& ymomside=(*side_bucket_mom[1])[mfi];
  FArrayBox& zmomside=(*side_bucket_mom[AMREX_SPACEDIM-1])[mfi];

  FArrayBox& xmassside=(*side_bucket_mass[0])[mfi];
  FArrayBox& ymassside=(*side_bucket_mass[1])[mfi];
  FArrayBox& zmassside=(*side_bucket_mass[AMREX_SPACEDIM-1])[mfi];

  FArrayBox& xmac_new=(*umac_new[0])[mfi];
  FArrayBox& ymac_new=(*umac_new[1])[mfi];
  FArrayBox& zmac_new=(*umac_new[AMREX_SPACEDIM-1])[mfi];

  FArrayBox& xmac_old=(*localMF[UMACOLD_MF])[mfi];
  FArrayBox& ymac_old=(*localMF[UMACOLD_MF+1])[mfi];
  FArrayBox& zmac_old=(*localMF[UMACOLD_MF+AMREX_SPACEDIM-1])[mfi];

  prescribed_vel_time_slab=0.5*(prev_time_slab+cur_time_slab);

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // solid distance function and solid moments are not modified.
   // solid temperature is modified only if solidheat_flag==0.
  fort_vfrac_split(
   &nprocessed[tid_current],
   &tid_current,
   density_floor.dataPtr(),
   density_ceiling.dataPtr(),
   &solidheat_flag, //0==diffuse in solid 1==dirichlet 2==neumann
   freezing_model.dataPtr(),
   distribute_from_target.dataPtr(),
   constant_density_all_time.dataPtr(),
   velbc.dataPtr(),
   &divu_outer_sweeps,
   &num_divu_outer_sweeps,
   &EILE_flag,
   &dir_absolute_direct_split,
   &normdir_here,
   tilelo,tilehi,
   fablo,fabhi,
   &bfact,
   &bfact_f,
   &dt_slab, // fort_vfrac_split
   &prev_time_slab,
   &cur_time_slab,
   &prescribed_vel_time_slab,
     // this is the original data
   LSfab.dataPtr(), //LS_RECON_MF, ngrow_scalar
   ARLIM(LSfab.loVect()),ARLIM(LSfab.hiVect()),
   denfab.dataPtr(),
   ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
   mom_denfab.dataPtr(),
   ARLIM(mom_denfab.loVect()),ARLIM(mom_denfab.hiVect()),
   tenfab.dataPtr(),
   ARLIM(tenfab.loVect()),ARLIM(tenfab.hiVect()),
   velfab.dataPtr(), //VELADVECT_MF
   ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
     // slope data
   vofslopefab.dataPtr(),
   ARLIM(vofslopefab.loVect()),ARLIM(vofslopefab.hiVect()),
     // this is the result
   destfab.dataPtr(),
   ARLIM(destfab.loVect()),ARLIM(destfab.hiVect()),
   tennewfab.dataPtr(),
   ARLIM(tennewfab.loVect()),ARLIM(tennewfab.hiVect()),
   LSdestfab.dataPtr(),
   ARLIM(LSdestfab.loVect()),ARLIM(LSdestfab.hiVect()),
    // other vars.
   vof0fab.dataPtr(),ARLIM(vof0fab.loVect()),ARLIM(vof0fab.hiVect()),
   maskfab.dataPtr(),ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
   masknbrfab.dataPtr(),
   ARLIM(masknbrfab.loVect()),ARLIM(masknbrfab.hiVect()),
   umac_displace.dataPtr(),
   ARLIM(umac_displace.loVect()),
   ARLIM(umac_displace.hiVect()),
   xlo,dx,
    // local variables
   consfab.dataPtr(),ARLIM(consfab.loVect()),ARLIM(consfab.hiVect()),
   xvelfab.dataPtr(),ARLIM(xvelfab.loVect()),ARLIM(xvelfab.hiVect()),
   yvelfab.dataPtr(),ARLIM(yvelfab.loVect()),ARLIM(yvelfab.hiVect()),
   zvelfab.dataPtr(),ARLIM(zvelfab.loVect()),ARLIM(zvelfab.hiVect()),
   xmomside.dataPtr(),ARLIM(xmomside.loVect()),ARLIM(xmomside.hiVect()),
   ymomside.dataPtr(),ARLIM(ymomside.loVect()),ARLIM(ymomside.hiVect()),
   zmomside.dataPtr(),ARLIM(zmomside.loVect()),ARLIM(zmomside.hiVect()),
   xmassside.dataPtr(),ARLIM(xmassside.loVect()),ARLIM(xmassside.hiVect()),
   ymassside.dataPtr(),ARLIM(ymassside.loVect()),ARLIM(ymassside.hiVect()),
   zmassside.dataPtr(),ARLIM(zmassside.loVect()),ARLIM(zmassside.hiVect()),
     //umac_new[0..sdim-1]
   xmac_new.dataPtr(),ARLIM(xmac_new.loVect()),ARLIM(xmac_new.hiVect()),
   ymac_new.dataPtr(),ARLIM(ymac_new.loVect()),ARLIM(ymac_new.hiVect()),
   zmac_new.dataPtr(),ARLIM(zmac_new.loVect()),ARLIM(zmac_new.hiVect()),
    //UMACOLD_MF
   xmac_old.dataPtr(),ARLIM(xmac_old.loVect()),ARLIM(xmac_old.hiVect()),
   ymac_old.dataPtr(),ARLIM(ymac_old.loVect()),ARLIM(ymac_old.hiVect()),
   zmac_old.dataPtr(),ARLIM(zmac_old.loVect()),ARLIM(zmac_old.hiVect()),
   &stokes_flow,
   denconst_interface_added.dataPtr(), //unused in fort_vfrac_split
   &ngrow_mass, //=2
   &ngrow_mac_old,
   &nc_conserve,
   &map_forward_direct_split[normdir_here],
   &vofrecon_ncomp,
   &den_recon_ncomp,
   &ncomp_state,
   &nc_bucket,
   &verbose,
   &gridno,&ngrid,
   &level,
   &finest_level,
   dombc.dataPtr(), 
   domlo,domhi);

 }  // mfi
} // omp
 ns_reconcile_d_num(LOOP_VFRAC_SPLIT,"split_scalar_advection");

 for (int tid=1;tid<thread_class::nthreads;tid++) {
  nprocessed[0]+=nprocessed[tid];
 }
 ParallelDescriptor::ReduceIntSum(nprocessed[0]);

 if (profile_debug==1) {
  double profile_time_end=ParallelDescriptor::second();
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "level= " << level << '\n';
   std::cout << "nprocessed= " << nprocessed[0] << '\n';
   std::cout << "profile fort_vfrac_split time = " << 
     profile_time_end-profile_time_start << '\n';
  }
 }

 delete conserve;
 
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   delete xvel[dir];
   delete side_bucket_mom[dir];
   delete side_bucket_mass[dir];
 }

 delete_localMF(VELADVECT_MF,1);
 delete_localMF(DEN_RECON_MF,1);
 delete_localMF(MOM_DEN_MF,1);

 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=num_materials)) {
  delete_localMF(TENSOR_RECON_MF_local,1);
 } else if (num_materials_viscoelastic==0) {
  // do nothing
 } else
  amrex::Error("num_materials_viscoelastic invalid:split_scalar_advection");
 
 delete_localMF(LS_RECON_MF,1);
 
 delete_localMF(UMACOLD_MF,AMREX_SPACEDIM);
 
 if ((level>=0)&&(level<finest_level)) {

  int spectral_override=1; // order derived from "enable_spectral"
   //Umac_Type
  avgDownMacState(spectral_override);
 
  avgDown(LS_Type,0,num_materials,0);
  MOFavgDown();
     // velocity and pressure
  avgDown(State_Type,STATECOMP_VEL,STATE_NCOMP_VEL+STATE_NCOMP_PRES,1);
  avgDown(State_Type,STATECOMP_STATES,num_state_material*num_materials,1);
  if ((num_materials_viscoelastic>=1)&&
      (num_materials_viscoelastic<=num_materials)) {
    // spectral_override==0 => always low order
   avgDown(Tensor_Type,0,NUM_CELL_ELASTIC,0);
  } else if (num_materials_viscoelastic==0) {
   // do nothing
  } else
   amrex::Error("num_materials_viscoelastic invalid:split_scalar_advection");

 } else if (level==finest_level) {
  // do nothing
 } else
  amrex::Error("level invalid23");

}  // subroutine split_scalar_advection

void
NavierStokes::errorEst (TagBoxArray& tags,int clearval,int tagval,
 Real time,int n_error_buf,int ngrow)
{
 
 const int max_level = parent->maxLevel();
 if (level>=max_level)
  amrex::Error("level too big in errorEst");

 int bfact=parent->Space_blockingFactor(level);

 if (level==max_level-1) {
  if (n_error_buf<2)
   amrex::Error("amr.n_error_buf<2 on level==max_level-1");
 } else if ((level>=0)&&(level<max_level-1)) {
  if (n_error_buf<2)
   amrex::Error("amr.n_error_buf<2 on level<max_level-1");
 } else
  amrex::Error("level invalid in errorEst");

 const int*  domain_lo = geom.Domain().loVect();
 const int*  domain_hi = geom.Domain().hiVect();
 const Real* dx        = geom.CellSize();
 const Real* prob_lo   = geom.ProbLo();

 Vector<int> error_set_count;
 error_set_count.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  error_set_count[tid]=0;
 }

 int local_time_order=parent->Time_blockingFactor();
 Real nudge_time=state[State_Type].slabTime(local_time_order);

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 if (STATECOMP_ERR+1==S_new.nComp()) {
  // do nothing
 } else
  amrex::Error("S_new.nComp() invalid");

 MultiFab* err_mf = getState(0,STATECOMP_ERR,1,nudge_time); 
 if (err_mf->nComp()==1) {
  // do nothing
 } else
  amrex::Error("err_mf->nComp()=!1");

 Real err_norm1=err_mf->norm1();

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "errorEst: clearval,tagval,time,n_error_buf,ngrow " <<
     clearval << ' ' << tagval << ' ' << time << ' ' << n_error_buf << ' ' <<
     ngrow << '\n';
   std::cout << "errorEst: max_level,level,bfact,time order,nudge_time " <<
     max_level << ' ' << level << ' ' << bfact << ' ' <<  
     local_time_order << ' ' << nudge_time << '\n';
   std::cout << "errorEst: slab_step " << slab_step << '\n';
   std::cout << "errorEst: err_norm1 " << err_norm1 << '\n';
  }
 }

 bool use_tiling=ns_tiling;

 if (1==1) {
  use_tiling=false;
 }

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(err_mf->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  // itags: Vector<int> ("domain"=box dimensions (not whole domain))
  // itags=TagBox::CLEAR in "tags"
  // then itags = older tag values.
 Vector<int>  itags;

 for (MFIter mfi(*err_mf,use_tiling); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  const int gridno = mfi.index();

  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];

  if (tilegrid==fabgrid) {
   // do nothing
  } else
   amrex::Error("FIX errorEST for tiling");

  const int* tilelo=tilegrid.loVect(); 
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();
  const Real* xlo = grid_loc[gridno].lo();

  TagBox& tagfab = tags[mfi];

  // We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
  // So we are going to get a temporary integer array.
  tagfab.get_itags(itags, tilegrid);

  int*        tptr  = itags.dataPtr();
  const int*  tlo   = tilegrid.loVect();
  const int*  thi   = tilegrid.hiVect();

  FArrayBox& errfab=(*err_mf)[mfi];
  const Box& errfab_box = errfab.box();
  if (errfab_box==tilegrid) {
   // do nothing
  } else
   amrex::Error("FIX errorEST for tiling");

  Real*       errfab_dat = errfab.dataPtr();
  const int*  errlo   = errfab_box.loVect();
  const int*  errhi   = errfab_box.hiVect();

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // declared in PROB.F90
  fort_vfracerror(
    &tid_current,
    &error_set_count[tid_current],
    tptr, 
    ARLIM(tlo), ARLIM(thi), 
    &tagval, &clearval, 
    errfab_dat, 
    ARLIM(errlo), ARLIM(errhi),
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    domain_lo, domain_hi,
    dx, xlo, prob_lo, 
    &upper_slab_time, 
    &level,
    &max_level,
    &max_level_for_use,
    &nblocks,
    xblocks.dataPtr(),yblocks.dataPtr(),zblocks.dataPtr(),
    rxblocks.dataPtr(),ryblocks.dataPtr(),rzblocks.dataPtr(),
    &ncoarseblocks,
    xcoarseblocks.dataPtr(),ycoarseblocks.dataPtr(),
    zcoarseblocks.dataPtr(),
    rxcoarseblocks.dataPtr(),rycoarseblocks.dataPtr(),
    rzcoarseblocks.dataPtr());

  tagfab.tags_and_untags(itags,tilegrid);
 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_VFRACERROR,"errorEst");
 for (int tid=1;tid<thread_class::nthreads;tid++) {
  error_set_count[0]+=error_set_count[tid];
 }
 ParallelDescriptor::Barrier();

 ParallelDescriptor::ReduceIntSum(error_set_count[0]);
 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "errorEst: level, error_set_count " << level << ' ' <<
    error_set_count[0] << '\n';
  }
 }

 delete err_mf;

} // end subroutine errorEst

void NavierStokes::GetDragALL() {

 std::string local_caller_string="GetDragALL";

//
// compute (-pI+2\mu D)dot n_solid for cut cells.
//
// F=mA  rho du/dt = div sigma + rho g
// integral rho du/dt = integral (div sigma + rho g)
// du/dt=(integral_boundary sigma dot n)/mass  + g
// torque_axis=r_axis x F=r_axis x (m alpha_axis x r_axis)   
// alpha=angular acceleration  (1/s^2)
// torque_axis=I_axis alpha_axis  
// I_axis=integral rho r_axis x (e_axis x r_axis)   I=moment of inertia
// I_axis=integral rho r_axis^{2}
// r_axis=(x-x_{COM})_{axis}  y_axis=y-(y dot e_axis)e_axis
// torque_axis=integral r_axis x (div sigma + rho g)=
//   integral_boundary r_axis x sigma dot n + integral r_axis x rho g
// NS_DRAG_integrated_quantities: see <DRAG_COMP.H>

 int finest_level=parent->finestLevel();

 if (level!=0)
  amrex::Error("it is required that level=0 in GetDragALL");

 debug_ngrow(CELL_VISC_MATERIAL_MF,1,local_caller_string);
 if (localMF[CELL_VISC_MATERIAL_MF]->nComp()==3*num_materials) {
  // do nothing
 } else 
  amrex::Error("GetDragALL: CELL_VISC_MATERIAL_MF invalid ncomp");
 
  //ngrow_make_distance=3
  //ngrow_distance=4
 debug_ngrow(DRAG_MF,ngrow_make_distance,local_caller_string);
 debug_ixType(DRAG_MF,-1,local_caller_string);
 if (localMF[DRAG_MF]->nComp()==N_DRAG) {
  // do nothing
 } else 
  amrex::Error("GetDragALL: DRAG_MF invalid ncomp");
 

 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=num_materials)) {
  debug_ngrow(VISCOTEN_ALL_MAT_MF,1,local_caller_string);
  if (localMF[VISCOTEN_ALL_MAT_MF]->nComp()==
      num_materials_viscoelastic*ENUM_NUM_TENSOR_TYPE) {
   // do nothing
  } else 
   amrex::Error("GetDragALL: VISCOTEN_ALL_MAT_MF invalid ncomp");
 } else if (num_materials_viscoelastic==0) {
  // do nothing
 } else
  amrex::Error("num_materials_viscoelastic invalid:GetDragALL");

 for (int isweep_drag=0;isweep_drag<2;isweep_drag++) {
  for (int ilev=level;ilev<=finest_level;ilev++) {
   NavierStokes& ns_level=getLevel(ilev);
   ns_level.GetDrag(isweep_drag);
  }
 }

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
   // declared in: NavierStokes2.cpp
  ns_level.avgDownDRAG_MF();
 }

 allocate_levelset_ALL(ngrow_distance,HOLD_LS_DATA_MF);
 if (localMF[HOLD_LS_DATA_MF]->nComp()!=num_materials*(AMREX_SPACEDIM+1))
  amrex::Error("hold_LS_DATA_MF (nComp()) !=num_materials*(AMREX_SPACEDIM+1)");
 debug_ngrow(HOLD_LS_DATA_MF,ngrow_distance,local_caller_string);

 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.level_DRAG_extend();
 }
 delete_array(HOLD_LS_DATA_MF);

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "DRAGCOMP num_materials=" << num_materials << '\n';
   for (int iq=0;iq<N_DRAG_IQ;iq++) {
    if (iq==DRAGCOMP_IQ_BODYFORCE) {
     std::cout << "DRAGCOMP_IQ_BODYFORCE 3xnum_materials\n";
    } else if (iq==DRAGCOMP_IQ_FORCE) {
     std::cout << "DRAGCOMP_IQ_FORCE 3xnum_materials\n";
    } else if (iq==DRAGCOMP_IQ_PFORCE) {
     std::cout << "DRAGCOMP_IQ_PFORCE 3xnum_materials\n";
    } else if (iq==DRAGCOMP_IQ_VISCOUSFORCE) {
     std::cout << "DRAGCOMP_IQ_VISCOUSFORCE 3xnum_materials\n";
    } else if (iq==DRAGCOMP_IQ_VISCOUS0FORCE) {
     std::cout << "DRAGCOMP_IQ_VISCOUS0FORCE 3xnum_materials\n";
    } else if (iq==DRAGCOMP_IQ_VISCOFORCE) {
     std::cout << "DRAGCOMP_IQ_VISCOFORCE 3xnum_materials\n";
    } else if (iq==DRAGCOMP_IQ_BODYTORQUE) {
     std::cout << "DRAGCOMP_IQ_BODYTORQUE 3xnum_materials\n";
    } else if (iq==DRAGCOMP_IQ_TORQUE) {
     std::cout << "DRAGCOMP_IQ_TORQUE 3xnum_materials\n";
    } else if (iq==DRAGCOMP_IQ_PTORQUE) {
     std::cout << "DRAGCOMP_IQ_PTORQUE 3xnum_materials\n";
    } else if (iq==DRAGCOMP_IQ_VISCOUSTORQUE) {
     std::cout << "DRAGCOMP_IQ_VISCOUSTORQUE 3xnum_materials\n";
    } else if (iq==DRAGCOMP_IQ_VISCOUS0TORQUE) {
     std::cout << "DRAGCOMP_IQ_VISCOUS0TORQUE 3xnum_materials\n";
    } else if (iq==DRAGCOMP_IQ_VISCOTORQUE) {
     std::cout << "DRAGCOMP_IQ_VISCOTORQUE 3xnum_materials\n";
    } else if (iq==DRAGCOMP_IQ_COM) {
     std::cout << "DRAGCOMP_IQ_COM 3xnum_materials\n";
    } else if (iq==DRAGCOMP_IQ_MOMINERTIA) {
     std::cout << "DRAGCOMP_IQ_MOMINERTIA 3xnum_materials\n";
    } else if (iq==DRAGCOMP_IQ_MASS) {
     std::cout << "DRAGCOMP_MASS num_materials\n";
    } else if (iq==DRAGCOMP_IQ_PERIM) {
     std::cout << "DRAGCOMP_PERIM num_materials\n";
    } else {
     //do nothing
    }
    int drag_im=-1;
    int drag_type=fort_drag_IQ_type(&iq,&drag_im);

    std::cout << "GetDragALL  iq= " << iq << " drag_im(0..num_materials-1)= " <<
	    drag_im << " drag_type= " << drag_type <<
     " NS_DRAG_integrated_quantities= " <<
     NS_DRAG_integrated_quantities[iq] << '\n';
   }
   for (int im=0;im<num_materials;im++) {
    Real mass=NS_DRAG_integrated_quantities[DRAGCOMP_IQ_MASS+im];
    if (mass>0.0) {
     for (int idir=0;idir<3;idir++) {
      std::cout << "COM im,idir= " << im << ' ' << idir << " COM= " <<
        NS_DRAG_integrated_quantities[DRAGCOMP_IQ_COM+im*3+idir]/mass << '\n';
     } //idir=0,1,2
    } //mass>0.0

   } //im=0;im<num_materials
  } //IOproc
 } //verbose>0

} // end subroutine GetDragALL

// sweep=0: integral (rho x), integral (rho) 
// sweep=1: find force, torque, moments of inertia, center of mass,mass
void
NavierStokes::GetDrag(int isweep) {

 std::string local_caller_string="GetDrag";

 bool use_tiling=ns_tiling;
 int bfact=parent->Space_blockingFactor(level);
 int finest_level=parent->finestLevel();
 if ((level<=finest_level)&&(level>=0)) {
  // do nothing
 } else
  amrex::Error("level or finest_level invalid");

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 debug_ngrow(CELLTENSOR_MF,1,local_caller_string);
 debug_ngrow(FACETENSOR_MF,1,local_caller_string);
 resize_levelset(2,LEVELPC_MF);
 debug_ngrow(LEVELPC_MF,2,local_caller_string);
 if (localMF[LEVELPC_MF]->nComp()!=num_materials*(1+AMREX_SPACEDIM))
  amrex::Error("levelpc mf has incorrect ncomp");
 VOF_Recon_resize(1); //output:SLOPE_RECON_MF
 debug_ngrow(SLOPE_RECON_MF,1,local_caller_string);

 debug_ngrow(CELL_VISC_MATERIAL_MF,1,local_caller_string);
 debug_ngrow(CELL_VISC_MF,1,local_caller_string);

 if (localMF[CELL_VISC_MATERIAL_MF]->nComp()==3*num_materials) {
  // do nothing
 } else {
  amrex::Error("GetDrag: CELL_VISC_MATERIAL_MF invalid ncomp");
 }

 if (localMF[CELL_VISC_MF]->nComp()==1) {
  // do nothing
 } else {
  amrex::Error("GetDrag: CELL_VISC_MF invalid ncomp");
 }

 resize_metrics(1);
 debug_ngrow(VOLUME_MF,1,local_caller_string);
 debug_ngrow(DRAG_MF,ngrow_make_distance,local_caller_string);
 debug_ixType(DRAG_MF,-1,local_caller_string);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(FACE_VAR_MF+dir,0,local_caller_string);

 int nparts=im_solid_map.size();
 if ((nparts<0)||(nparts>num_materials))
  amrex::Error("nparts invalid");
 Vector<int> im_solid_map_null;
 im_solid_map_null.resize(1);

 int* im_solid_map_ptr;
 int nparts_def=nparts;
 if (nparts==0) {
  im_solid_map_ptr=im_solid_map_null.dataPtr();
  nparts_def=1;
 } else if ((nparts>=1)&&(nparts<=num_materials)) {
  im_solid_map_ptr=im_solid_map.dataPtr();
 } else
  amrex::Error("nparts invalid");

 if (localMF[CELLTENSOR_MF]->nComp()!=AMREX_SPACEDIM_SQR)
  amrex::Error("localMF[CELLTENSOR_MF]->nComp() invalid");
 if (localMF[FACETENSOR_MF]->nComp()!=AMREX_SPACEDIM_SQR)
  amrex::Error("localMF[FACETENSOR_MF]->nComp() invalid");

 int nstate=STATE_NCOMP;
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");

 if (NS_DRAG_integrated_quantities.size()!=N_DRAG_IQ)
  amrex::Error("NS_DRAG_integrated_quantities invalid size");
 if (NS_DRAG_integrated_quantities_sweep.size()!=N_DRAG_IQ)
  amrex::Error("NS_DRAG_integrated_quantities_sweep invalid size");
 
 Vector< Vector<Real> > local_integrated_quantities;
 local_integrated_quantities.resize(thread_class::nthreads);

 for (int tid=0;tid<thread_class::nthreads;tid++) {
  local_integrated_quantities[tid].resize(N_DRAG_IQ);
  for (int iq=0;iq<N_DRAG_IQ;iq++)
   local_integrated_quantities[tid][iq]=0.0;
 } // tid

// gear problem: probtype=563, axis_dir=2, 3D
// scale torque by 2 pi vinletgas/60

 if (localMF[DRAG_MF]->nComp()!=N_DRAG)
  amrex::Error("drag ncomp invalid");

 const Real* dx = geom.CellSize();

 int combine_flag=2;
 int hflag=0;
 int combine_idx=-1;  // update state variables
 int update_flux=0;
 int interface_cond_avail=0;

 combine_state_variable(
  SOLVETYPE_VISC,
  combine_idx,
  combine_flag,
  hflag,
  update_flux,
  interface_cond_avail);
 update_flux=1;
 combine_state_variable(
  SOLVETYPE_PRES,
  combine_idx,
  combine_flag,
  hflag,
  update_flux,
  interface_cond_avail);

  // p(rho,T) in compressible parts
  // projection pressure in incompressible parts.
 MultiFab* pres=getStatePres(1,cur_time_slab); 

 if (pres->nComp()!=1)
  amrex::Error("pres->nComp() invalid");
 
 MultiFab* vel=getState(2,0,AMREX_SPACEDIM,cur_time_slab);
  // mask=tag if not covered by level+1 or outside the domain.
 int ngrowmask=2;
 int clear_phys_boundary=0;
 Real tag=1.0;
 MultiFab* mask=maskfiner(ngrowmask,tag,clear_phys_boundary);  
 MultiFab* den_recon=getStateDen(1,cur_time_slab);  
 clear_phys_boundary=3;
  // mask=tag at exterior fine/fine border.
  // mask=1-tag at other exterior boundaries.
 MultiFab* mask3=maskfiner(ngrowmask,tag,clear_phys_boundary);

 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
  if (localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0)
   amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0");
  if (localMF[FSI_GHOST_MAC_MF+data_dir]->nComp()!=nparts_def*AMREX_SPACEDIM)
   amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nComp() invalid");
 }

 if (ENUM_NUM_TENSOR_TYPE==2*AMREX_SPACEDIM) {
  // do nothing
 } else
  amrex::Error("ENUM_NUM_TENSOR_TYPE invalid");

 int VISCOTEN_ALL_MAT_MF_local=0;
 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=num_materials)) {
  VISCOTEN_ALL_MAT_MF_local=VISCOTEN_ALL_MAT_MF;

  debug_ngrow(VISCOTEN_ALL_MAT_MF_local,1,local_caller_string);
  if (localMF[VISCOTEN_ALL_MAT_MF_local]->nComp()==NUM_CELL_ELASTIC) {
   //do nothing
  } else {
   amrex::Error("VISCOTEN_ALL_MAT_MF_local invalid nComp");
  }
 } else if (num_materials_viscoelastic==0) {
  VISCOTEN_ALL_MAT_MF_local=CELLTENSOR_MF;
 } else
  amrex::Error("num_materials_viscoelastic invalid:GetDrag");

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(localMF[DRAG_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*localMF[DRAG_MF],use_tiling); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();

  const Real* xlo = grid_loc[gridno].lo();
  Vector<int> velbc=getBCArray(State_Type,gridno,
     STATECOMP_VEL,STATE_NCOMP_VEL);

  FArrayBox& maskfab=(*mask)[mfi];
  FArrayBox& volfab=(*localMF[VOLUME_MF])[mfi];
  FArrayBox& areax=(*localMF[AREA_MF])[mfi];
  FArrayBox& areay=(*localMF[AREA_MF+1])[mfi];
  FArrayBox& areaz=(*localMF[AREA_MF+AMREX_SPACEDIM-1])[mfi];

  FArrayBox& mufab=(*localMF[CELL_VISC_MF])[mfi];
  FArrayBox& mu_mat_fab=(*localMF[CELL_VISC_MATERIAL_MF])[mfi];
  if (mu_mat_fab.nComp()==3*num_materials) {
   // do nothing
  } else {
   amrex::Error("mu_mat_fab.nComp()==3*num_materials is false");
  }

  FArrayBox& xface=(*localMF[FACE_VAR_MF])[mfi];
  FArrayBox& yface=(*localMF[FACE_VAR_MF+1])[mfi];
  FArrayBox& zface=(*localMF[FACE_VAR_MF+AMREX_SPACEDIM-1])[mfi];

  FArrayBox& solxfab=(*localMF[FSI_GHOST_MAC_MF])[mfi]; 
  FArrayBox& solyfab=(*localMF[FSI_GHOST_MAC_MF+1])[mfi];
  FArrayBox& solzfab=(*localMF[FSI_GHOST_MAC_MF+AMREX_SPACEDIM-1])[mfi];

  FArrayBox& denfab=(*den_recon)[mfi];
  FArrayBox& presfab=(*pres)[mfi];
  FArrayBox& velfab=(*vel)[mfi];
  FArrayBox& dragfab=(*localMF[DRAG_MF])[mfi];

  FArrayBox& voffab=(*localMF[SLOPE_RECON_MF])[mfi];
  FArrayBox& levelpcfab=(*localMF[LEVELPC_MF])[mfi];

  FArrayBox& tensor_data=(*localMF[CELLTENSOR_MF])[mfi];
  FArrayBox& elastic_tensor_data=
     (*localMF[VISCOTEN_ALL_MAT_MF_local])[mfi];

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // hoop stress, centripetal force, coriolis effect still not
   // considered.
   // fort_getdrag is declared in: DERIVE_3D.F90
  fort_getdrag(
   &tid_current,
   &level,
   &finest_level,
   &isweep,
   NS_DRAG_integrated_quantities.dataPtr(),
   NS_DRAG_integrated_quantities_sweep.dataPtr(),
   local_integrated_quantities[tid_current].dataPtr(),
   tensor_data.dataPtr(),
   ARLIM(tensor_data.loVect()),ARLIM(tensor_data.hiVect()),
   elastic_tensor_data.dataPtr(),
   ARLIM(elastic_tensor_data.loVect()),
   ARLIM(elastic_tensor_data.hiVect()),
   denfab.dataPtr(),ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
   maskfab.dataPtr(),ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
   voffab.dataPtr(),ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),
   levelpcfab.dataPtr(),
   ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),
   volfab.dataPtr(),ARLIM(volfab.loVect()),ARLIM(volfab.hiVect()),
   areax.dataPtr(),ARLIM(areax.loVect()),ARLIM(areax.hiVect()),
   areay.dataPtr(),ARLIM(areay.loVect()),ARLIM(areay.hiVect()),
   areaz.dataPtr(),ARLIM(areaz.loVect()),ARLIM(areaz.hiVect()),
   xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()),
   yface.dataPtr(),ARLIM(yface.loVect()),ARLIM(yface.hiVect()),
   zface.dataPtr(),ARLIM(zface.loVect()),ARLIM(zface.hiVect()),
   mufab.dataPtr(),ARLIM(mufab.loVect()),ARLIM(mufab.hiVect()),
   mu_mat_fab.dataPtr(),
   ARLIM(mu_mat_fab.loVect()),ARLIM(mu_mat_fab.hiVect()),
   xlo,dx,
   solxfab.dataPtr(),
   ARLIM(solxfab.loVect()),ARLIM(solxfab.hiVect()),
   solyfab.dataPtr(),
   ARLIM(solyfab.loVect()),ARLIM(solyfab.hiVect()),
   solzfab.dataPtr(),
   ARLIM(solzfab.loVect()),ARLIM(solzfab.hiVect()),
   presfab.dataPtr(),ARLIM(presfab.loVect()),ARLIM(presfab.hiVect()),
   velfab.dataPtr(),ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
   dragfab.dataPtr(),ARLIM(dragfab.loVect()),ARLIM(dragfab.hiVect()),
   tilelo,tilehi,
   fablo,fabhi,&bfact,
   &NS_geometry_coord,
   velbc.dataPtr(),&cur_time_slab,
   &visc_coef,
   &nparts,
   &nparts_def,
   im_solid_map_ptr);

 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_GETDRAG,"GetDrag");

 for (int tid=1;tid<thread_class::nthreads;tid++) {
  for (int iq=0;iq<N_DRAG_IQ;iq++) {
   if ((isweep==0)&&(NS_DRAG_integrated_quantities_sweep[iq]==0)) {
    local_integrated_quantities[0][iq]+=local_integrated_quantities[tid][iq];
   } else if ((isweep==1)&&(NS_DRAG_integrated_quantities_sweep[iq]==1)) {
    local_integrated_quantities[0][iq]+=local_integrated_quantities[tid][iq];
   } else if ((isweep==0)&&(NS_DRAG_integrated_quantities_sweep[iq]==1)) {
    // do nothing
   } else if ((isweep==1)&&(NS_DRAG_integrated_quantities_sweep[iq]==0)) {
    // do nothing
   } else
    amrex::Error("isweep or drag isweep error");
  } // iq
 } // tid

 ParallelDescriptor::Barrier();

 for (int iq=0;iq<N_DRAG_IQ;iq++) {
  if ((isweep==0)&&(NS_DRAG_integrated_quantities_sweep[iq]==0)) {
   ParallelDescriptor::ReduceRealSum(local_integrated_quantities[0][iq]);
   NS_DRAG_integrated_quantities[iq]+=local_integrated_quantities[0][iq];
  } else if ((isweep==1)&&(NS_DRAG_integrated_quantities_sweep[iq]==1)) {
   ParallelDescriptor::ReduceRealSum(local_integrated_quantities[0][iq]);
   NS_DRAG_integrated_quantities[iq]+=local_integrated_quantities[0][iq];
  } else if ((isweep==0)&&(NS_DRAG_integrated_quantities_sweep[iq]==1)) {
   // do nothing
  } else if ((isweep==1)&&(NS_DRAG_integrated_quantities_sweep[iq]==0)) {
   // do nothing
  } else
   amrex::Error("isweep or drag isweep error");
 } // iq

 combine_flag=2;
 hflag=0;
 combine_idx=-1; 
 update_flux=0;

 combine_state_variable(
  SOLVETYPE_VISC,
  combine_idx,
  combine_flag,
  hflag,
  update_flux,
  interface_cond_avail);

 update_flux=1;
 combine_state_variable(
  SOLVETYPE_PRES,
  combine_idx,
  combine_flag,
  hflag,
  update_flux,
  interface_cond_avail);
 
 delete pres; 
 delete vel; 
 delete den_recon; 
 delete mask; 
 delete mask3; 

} // end subroutine GetDrag

//called from:
//NavierStokes::multiphase_SHELL_preconditioner
//NavierStokes::multiphase_preconditioner
//NavierStokes::multiphase_project
//NavierStokes::JacobiALL
//NavierStokes::residALL
//NavierStokes::applyALL
void NavierStokes::project_right_hand_side(
  int index_MF,int project_option,int& change_flag) {

 int finest_level=parent->finestLevel();
 if ((level==0)&&(finest_level>=level)) {
  // do nothing
 } else
  amrex::Error("level or finest_level invalid");

 change_flag=0;

  //SOLVETYPE_PRES, SOLVETYPE_PRESGRAVITY, SOLVETYPE_INITPROJ,
  //SOLVETYPE_PRESEXTRAP
 if (project_option_singular_possible(project_option)==1) {

  int at_least_one_active=0;
  int at_least_one_zap=0;

  for (int icolor=0;icolor<color_ONES_count;icolor++) {

   if (ones_sum_global[icolor]>=1.0) {

    if (verbose>0) {
     if (ParallelDescriptor::IOProcessor()) {
      std::cout << "prj_rhs,icolor="<<icolor<<" ones_sum_global="<<
       ones_sum_global[icolor]<<" singular_patch_flag="<<
       singular_patch_flag[icolor]<<'\n';
     }
    } //verbose>0

   } else
    amrex::Error("ones_sum_global[icolor] invalid");

    //patch contains prescribed solid.
   if (singular_patch_flag[icolor]==0) {

    at_least_one_zap=1;

     //patch contains no compressible or dirichlet cells.
   } else if (singular_patch_flag[icolor]==1) {

    at_least_one_active=1;

     //patch contains a compressible or dirichlet cell.
   } else if (singular_patch_flag[icolor]==2) {

    at_least_one_active=1;

   } else
    amrex::Error("singular_patch_flag invalid");

  } // icolor=0..color_ONES_count-1

  if (at_least_one_zap==1) {
   change_flag=1;
   zap_resid_where_singular(index_MF); // multiply by ONES_MF
  } else if (at_least_one_zap==0) {
   // do nothing
  } else
   amrex::Error("at_least_one_zap invalid");

  if (at_least_one_active==1) {

   change_flag=1;

    // rhsnew=rhs H-alpha H
    // 0 =sum rhs H-alpha sum H
    // alpha=sum rhs H / sum H
   Vector<Real> coef;
   coef.resize(color_ONES_count);
   dot_productALL_ones(project_option,index_MF,coef);

   for (int icolor=0;icolor<color_ONES_count;icolor++) {

     //patch is in a prescribed solid region.    
    if (singular_patch_flag[icolor]==0) {

     if (ones_sum_global[icolor]>=1.0) {
      // do nothing
     } else
      amrex::Error("ones_sum_global[icolor] invalid");

     //patch contains no compressible or dirichlet cells.
    } else if (singular_patch_flag[icolor]==1) {

     if (ones_sum_global[icolor]>=1.0) {

      if (verbose>0) {
       if (ParallelDescriptor::IOProcessor()) {
        std::cout << "prj_rhs, icolor=" << icolor << '\n';
        std::cout << "prj_rhs, coef=" << coef[icolor] << '\n';
        std::cout << "prj_rhs, ones_sum_global=" << 
         ones_sum_global[icolor] << '\n';
       } 
      } // verbose>0

      coef[icolor]=-coef[icolor]/ones_sum_global[icolor];

     } else
      amrex::Error("ones_sum_global[icolor] invalid");

     //patch contains a compressible or dirichlet cell.
    } else if (singular_patch_flag[icolor]==2) {

     if (ones_sum_global[icolor]>=1.0) {
      // do nothing
     } else
      amrex::Error("ones_sum_global[icolor] invalid");

     if ((coef[icolor]>=0.0)||(coef[icolor]<=0.0)) {
      coef[icolor]=0.0;
     } else {
      for (int icolor2=0;icolor2<color_ONES_count;icolor2++) {
       std::cout << "icolor2=" << icolor2 << " ones_sum_global[icolor2]=" <<
        ones_sum_global[icolor2] << " singular_patch_flag[icolor2]=" <<
        singular_patch_flag[icolor2] << '\n';
      }
      std::cout << "icolor= " << icolor << " coef[icolor]="
       << coef[icolor] << '\n';
      amrex::Error("coef[icolor] invalid");
     }

    } else
     amrex::Error("singular_patch_flag[icolor] invalid");

    if (verbose>0) {
     if (ParallelDescriptor::IOProcessor()) {
      std::cout << "project_right_hand_side, coef=" << coef[icolor] << '\n';
      std::cout << "project_right_hand_side, denom=" << 
         ones_sum_global[icolor] << '\n';
     } 
    } // verbose>0

   } // icolor=0..color_ONES_count-1

    //if singular_patch_flag(icolor)==2, then index_MF is not incremented.
   mf_combine_ones(project_option,index_MF,coef);
   zap_resid_where_singular(index_MF);

  } else if (at_least_one_active==0) {
   std::cout << "index_MF = " << index_MF << '\n';
   print_project_option(project_option);
   std::cout << "color_ONES_count = " << color_ONES_count << '\n';
   for (int icolor=0;icolor<color_ONES_count;icolor++) {
    std::cout << "prj_rhs,icolor="<<icolor<<" ones_sum_global="<<
      ones_sum_global[icolor]<<" singular_patch_flag="<<
      singular_patch_flag[icolor]<<'\n';
   }
   std::cout << "at_least_one_active = " << at_least_one_active << '\n';
   std::cout << "verify rigid body level set functions init. correctly\n";
   amrex::Error("at_least_one_active is 0");
  } else 
   amrex::Error("at least_one_active neither 0 nor 1");

 } else if (project_option_singular_possible(project_option)==0) {

  // do nothing

 } else
  amrex::Error("project_option_singular_possible inv project_right_hand_side");

} // end subroutine project_right_hand_side

// called from: NavierStokes::project_right_hand_side
void NavierStokes::dot_productALL_ones(int project_option,
  int index_MF,Vector<Real>& coef) {

 int finest_level=parent->finestLevel();
 if (level!=0)
  amrex::Error("level=0 in dot_productALL_ones");
 if (coarsest_ONES_level>finest_level)
  amrex::Error("coarsest_ONES_level invalid");

 if (color_ONES_count>0) {
  // do nothing
 } else {
  std::cout << "color_ONES_count=" << color_ONES_count << '\n';
  amrex::Error("color_ONES_count invalid 2");
 }

 for (int icolor=0;icolor<color_ONES_count;icolor++) {
  coef[icolor]=0.0;
 }

 Vector<Real> temp_sum;

 temp_sum.resize(color_ONES_count);

 for (int k = coarsest_ONES_level; k <= finest_level; k++) {
  NavierStokes& ns_level = getLevel(k);
  ns_level.dotSumONES(
   project_option,
   index_MF,
   temp_sum);

  for (int icolor=0;icolor<color_ONES_count;icolor++) {
   coef[icolor]+=temp_sum[icolor];
  }
 } // for (int k = coarsest_ONES_level; k <= finest_level; k++)

} // end subroutine dot_productALL_ones

// index_MF = index_MF + coef
// called from:
// NavierStokes::project_right_hand_side
void NavierStokes::mf_combine_ones(int project_option,
  int index_MF,Vector<Real> coef) {

 int finest_level=parent->finestLevel();
 if (level!=0)
  amrex::Error("level=0 in mf_combine_ones");
 if (coarsest_ONES_level>finest_level)
  amrex::Error("coarsest_ONES_level invalid");

 if (color_ONES_count>0) {
  // do nothing
 } else {
  std::cout << "color_ONES_count=" << color_ONES_count << '\n';
  amrex::Error("color_ONES_count invalid 3");
 }

 for (int k = coarsest_ONES_level; k <= finest_level; k++) {
  NavierStokes& ns_level = getLevel(k);
  ns_level.mf_combine_ones_level(
   project_option,
   index_MF,
   coef);

 } // for (int k = coarsest_ONES_level; k <= finest_level; k++)

} // end subroutine mf_combine_ones

void NavierStokes::zap_resid_where_singular(int index_MF) {

 int finest_level=parent->finestLevel();
 if (level!=0)
  amrex::Error("level=0 in zap_resid_where_singular");

 for (int k = 0; k <= finest_level; k++) {
  NavierStokes& ns_level = getLevel(k);
  MultiFab* dest=ns_level.localMF[index_MF];
  int nsolve=dest->nComp();
  if (nsolve==1) {
   //do nothing
  } else
   amrex::Error("nsolve invalid in zap_resid_where_singular");

  MultiFab::Multiply(*dest,*ns_level.localMF[ONES_MF],0,0,1,0);
 }  // k=0..finest_level

}  // subroutine zap_resid_where_singular

//called from: NavierStokes::multiphase_project
void NavierStokes::dot_productALL_ones_size(int project_option) {

 int finest_level=parent->finestLevel();
 if (level!=0)
  amrex::Error("level=0 in dot_productALL_ones");
 if (coarsest_ONES_level>finest_level)
  amrex::Error("coarsest_ONES_level invalid");

 if (color_ONES_count>0) {
  // do nothing
 } else {
  std::cout << "color_ONES_count=" << color_ONES_count << '\n';
  amrex::Error("color_ONES_count invalid 4");
 }

 Vector<Real> temp_sum;
 Vector<int> temp_flag;
 temp_sum.resize(color_ONES_count);
 temp_flag.resize(color_ONES_count);

 for (int icolor=0;icolor<color_ONES_count;icolor++) {
  singular_patch_flag[icolor]=0;
  ones_sum_global[icolor]=0.0;
 }

 for (int k = coarsest_ONES_level; k <= finest_level; k++) {
  NavierStokes& ns_level = getLevel(k);
  ns_level.dotSumONES_size(
   project_option,
   temp_sum,temp_flag);

  for (int icolor=0;icolor<color_ONES_count;icolor++) {
   singular_patch_flag[icolor]=
     std::max(singular_patch_flag[icolor],
	      temp_flag[icolor]);
   ones_sum_global[icolor]+=temp_sum[icolor];
  }
 } // for (int k = coarsest_ONES_level; k <= finest_level; k++)

} // end subroutine dot_productALL_ones_size


void NavierStokes::dot_productALL(int project_option,
 int index1_MF,int index2_MF,
 Real& result,int nsolve) {

 int finest_level=parent->finestLevel();
 if (level!=0)
  amrex::Error("level=0 in dot_productALL");

 Real tempsum=0.0;
 Real total_sum=0.0;

 for (int k = 0; k <= finest_level; k++) {
  NavierStokes& ns_level = getLevel(k);
  ns_level.dotSum( 
     project_option,
     ns_level.localMF[index1_MF],
     ns_level.localMF[index2_MF], tempsum,nsolve);
  total_sum+=tempsum;
 }

 result=total_sum;
}

void
NavierStokes::dotSum(int project_option,
  MultiFab* mf1, MultiFab* mf2, Real& result,int nsolve) {
 
 std::string local_caller_string="dotSum";

 bool use_tiling=ns_tiling;

 if ((nsolve!=1)&&(nsolve!=AMREX_SPACEDIM))
  amrex::Error("nsolve invalid");

 if (project_option_momeqn(project_option)==1) {
  //do nothing
 } else if (project_option_momeqn(project_option)==0) {
  //do nothing
 } else
  amrex::Error("project_option_momeqn(project_option) invalid2");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(FACE_VAR_MF+dir,0,local_caller_string);

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,0,local_caller_string);

 if (mf1->nComp()!=nsolve) {
  std::cout << "mf1_ncomp nsolve " << mf1->nComp() << ' ' << nsolve << '\n';
  amrex::Error("dotSum: mf1 invalid ncomp");
 }
 if (mf2->nComp()!=nsolve)
  amrex::Error("mf2 invalid ncomp");

 Vector<Real> sum;
 sum.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  sum[tid]=0.0;
 }

 int finest_level=parent->finestLevel();
 if (level>finest_level)
  amrex::Error("level too big");

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(mf1->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel 
#endif
{
 for (MFIter mfi(*mf1,use_tiling); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());

  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();
  int bfact=parent->Space_blockingFactor(level);

  FArrayBox& fab = (*mf1)[mfi];
  FArrayBox& fab2 = (*mf2)[mfi];

  Real tsum=0.0;
  FArrayBox& mfab=(*localMF[MASKCOEF_MF])[mfi];
  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // in: NAVIERSTOKES_3D.F90
  fort_sumdot(&tsum,
    fab.dataPtr(),ARLIM(fab.loVect()), ARLIM(fab.hiVect()),
    fab2.dataPtr(),ARLIM(fab2.loVect()), ARLIM(fab2.hiVect()),
    mfab.dataPtr(),ARLIM(mfab.loVect()),ARLIM(mfab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    &debug_dot_product,
    &level,&gridno,
    &nsolve);
  sum[tid_current] += tsum;
 } // mfi1
} // omp
 ns_reconcile_d_num(LOOP_SUMDOT,"dotSum");
 for (int tid=1;tid<thread_class::nthreads;tid++) {
  sum[0]+=sum[tid];
 }
 ParallelDescriptor::ReduceRealSum(sum[0]);

 result=sum[0];
} // end subroutine dotSum

//called from: NavierStokes::dot_productALL_ones_size
void
NavierStokes::dotSumONES_size(int project_option,
  Vector<Real>& result_sum,
  Vector<int>& result_flag) {

 std::string local_caller_string="dotSumONES_size";

 bool use_tiling=ns_tiling;

 if ((result_sum.size()==color_ONES_count)&&
     (result_flag.size()==color_ONES_count)) {
  // do nothing
 } else
  amrex::Error("result_sum or result_flag have invalid size");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(FACE_VAR_MF+dir,0,local_caller_string); 

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,0,local_caller_string);

 debug_ngrow(ONES_MF,0,local_caller_string);
 debug_ngrow(ONES_GROW_MF,0,local_caller_string);
 debug_ngrow(TYPE_ONES_MF,0,local_caller_string); //=1 if diagonal=0.0; =2 if diagonal>0.0
 debug_ngrow(COLOR_ONES_MF,0,local_caller_string);

 if (type_ONES_flag.size()==2) {
  // do nothing
 } else
  amrex::Error("type_ONES_flag.size() invalid");

 debug_ngrow(ALPHACOEF_MF,0,local_caller_string);
 int nsolve_test=localMF[ALPHACOEF_MF]->nComp();
 int nsolve_expect=1;
 if (project_option==SOLVETYPE_VISC) { 
  nsolve_expect=AMREX_SPACEDIM;
 } else if (project_option_is_valid(project_option)==1) {
  nsolve_expect=1;
 } else
  amrex::Error("project_option invalid 18580");

 if (nsolve_expect==nsolve_test) {
  // do nothing
 } else
  amrex::Error("nsolve_test invalid");

 Vector< Vector<Real> > local_sum;
 Vector< Vector<int> > local_flag;

  // for each given color, singular_patch_flag=
  //   0 if color is masked off 
  //   1 if color is not masked off, no compressible/internal dirichlet 
  //     regions, and not touching a Dirichlet condition wall.
  //   2 if color is not masked off, a compressible/internal dirichlet
  //     region exists or color is touching a Dirichlet condition wall.
 local_sum.resize(thread_class::nthreads);
 local_flag.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  local_sum[tid].resize(color_ONES_count);
  local_flag[tid].resize(color_ONES_count);
  for (int icolor=0;icolor<color_ONES_count;icolor++) {
   local_sum[tid][icolor]=0.0;
   local_flag[tid][icolor]=0;
  }
 } // tid=0..nthreads-1

 int finest_level=parent->finestLevel();
 if (level>finest_level)
  amrex::Error("level too big");

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(localMF[ONES_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel 
#endif
{
 for (MFIter mfi(*localMF[ONES_MF],use_tiling); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());

  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();
  int bfact=parent->Space_blockingFactor(level);

  FArrayBox& ones_fab = (*localMF[ONES_MF])[mfi];
  FArrayBox& type_fab = (*localMF[TYPE_ONES_MF])[mfi];
  FArrayBox& color_fab = (*localMF[COLOR_ONES_MF])[mfi];
  FArrayBox& alpha_fab = (*localMF[ALPHACOEF_MF])[mfi];
  FArrayBox& mfab=(*localMF[MASKCOEF_MF])[mfi];

  Vector<int> presbc=getBCArray(State_Type,gridno,STATECOMP_PRES,1);

  Vector<Real> fab_sum;
  fab_sum.resize(color_ONES_count);
  Vector<int> fab_flag;
  fab_flag.resize(color_ONES_count);
  for (int icolor=0;icolor<color_ONES_count;icolor++) {
   fab_sum[icolor]=0.0;
   fab_flag[icolor]=0;
  }

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // in: NAVIERSTOKES_3D.F90
  fort_sumdot_ones_size(
   fab_sum.dataPtr(),
   fab_flag.dataPtr(),
   ones_fab.dataPtr(),
   ARLIM(ones_fab.loVect()),ARLIM(ones_fab.hiVect()),
   type_fab.dataPtr(),
   ARLIM(type_fab.loVect()),ARLIM(type_fab.hiVect()),
   color_fab.dataPtr(),
   ARLIM(color_fab.loVect()),ARLIM(color_fab.hiVect()),
   alpha_fab.dataPtr(),
   ARLIM(alpha_fab.loVect()),ARLIM(alpha_fab.hiVect()),
   mfab.dataPtr(),
   ARLIM(mfab.loVect()),ARLIM(mfab.hiVect()),
   tilelo,tilehi,
   fablo,fabhi,&bfact,
   &level,
   &gridno,
   &nsolve_expect, 
   presbc.dataPtr(),
   type_ONES_flag.dataPtr(),
   &color_ONES_count,
   &project_option);

  for (int icolor=0;icolor<color_ONES_count;icolor++) {
   local_sum[tid_current][icolor] += fab_sum[icolor];
   local_flag[tid_current][icolor]=std::max(
     local_flag[tid_current][icolor],
     fab_flag[icolor]);
  } // icolor=0;icolor<color_ONES_count
 } // mfi1
} // omp
 ns_reconcile_d_num(LOOP_SUMDOT_ONES_SIZE,"dotSumONES_size");
 for (int tid=1;tid<thread_class::nthreads;tid++) {
  for (int icolor=0;icolor<color_ONES_count;icolor++) {
   local_sum[0][icolor]+=local_sum[tid][icolor];
   local_flag[0][icolor]=std::max(
     local_flag[0][icolor],
     local_flag[tid][icolor]);
  } // icolor=0;icolor<color_ONES_count
 }
 for (int icolor=0;icolor<color_ONES_count;icolor++) {
  ParallelDescriptor::ReduceRealSum(local_sum[0][icolor]);
  ParallelDescriptor::ReduceIntMax(local_flag[0][icolor]);
 }

 for (int icolor=0;icolor<color_ONES_count;icolor++) {
  result_sum[icolor]=local_sum[0][icolor];
  result_flag[icolor]=local_flag[0][icolor];
 }

} // end subroutine dotSumONES_size

//called by: NavierStokes::dot_productALL_ones
void
NavierStokes::dotSumONES(int project_option,
  int index_MF,
  Vector<Real>& result_sum) {
 
 std::string local_caller_string="dotSumONES";

 bool use_tiling=ns_tiling;

 if (result_sum.size()==color_ONES_count) {
  // do nothing
 } else
  amrex::Error("result_sum has invalid size");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(FACE_VAR_MF+dir,0,local_caller_string); 

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(index_MF,0,local_caller_string);
 debug_ngrow(MASKCOEF_MF,0,local_caller_string);

 debug_ngrow(ONES_MF,0,local_caller_string);
 debug_ngrow(ONES_GROW_MF,0,local_caller_string);
 debug_ngrow(TYPE_ONES_MF,0,local_caller_string); //=1 if diagonal=0.0; =2 if diagonal>0.0
 debug_ngrow(COLOR_ONES_MF,0,local_caller_string);

 if (type_ONES_flag.size()==2) {
  // do nothing
 } else
  amrex::Error("type_ONES_flag.size() invalid");

 debug_ngrow(ALPHACOEF_MF,0,local_caller_string);
 int nsolve_test=localMF[ALPHACOEF_MF]->nComp();
 int nsolve_expect=1;
 if (project_option==SOLVETYPE_VISC) { 
  nsolve_expect=AMREX_SPACEDIM;
 } else if (project_option_is_valid(project_option)==1) {
  nsolve_expect=1;
 } else
  amrex::Error("project_option invalid nsolve_expect nsolve_test");

 if (nsolve_expect==nsolve_test) {
  // do nothing
 } else
  amrex::Error("nsolve_test invalid");

 Vector< Vector<Real> > local_sum;

  // for each given color, singular_patch_flag=
  //   0 if color is masked off 
  //   1 if color is not masked off, no compressible/internal dirichlet 
  //     regions, and not touching a Dirichlet condition wall.
  //   2 if color is not masked off, a compressible/internal dirichlet
  //     region exists or color is touching a Dirichlet condition wall.
 local_sum.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  local_sum[tid].resize(color_ONES_count);
  for (int icolor=0;icolor<color_ONES_count;icolor++) {
   local_sum[tid][icolor]=0.0;
  }
 } // tid=0..nthreads-1

 int finest_level=parent->finestLevel();
 if (level>finest_level)
  amrex::Error("level too big");

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(localMF[ONES_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel 
#endif
{
 for (MFIter mfi(*localMF[ONES_MF],use_tiling); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());

  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();
  int bfact=parent->Space_blockingFactor(level);

  FArrayBox& data_fab = (*localMF[index_MF])[mfi];
  FArrayBox& ones_fab = (*localMF[ONES_MF])[mfi];
  FArrayBox& type_fab = (*localMF[TYPE_ONES_MF])[mfi];
  FArrayBox& color_fab = (*localMF[COLOR_ONES_MF])[mfi];
  FArrayBox& alpha_fab = (*localMF[ALPHACOEF_MF])[mfi];
  FArrayBox& mfab=(*localMF[MASKCOEF_MF])[mfi];

  Vector<int> presbc=getBCArray(State_Type,gridno,STATECOMP_PRES,1);

  Vector<Real> fab_sum;
  fab_sum.resize(color_ONES_count);
  Vector<int> fab_flag;
  fab_flag.resize(color_ONES_count);
  for (int icolor=0;icolor<color_ONES_count;icolor++) {
   fab_sum[icolor]=0.0;
   fab_flag[icolor]=singular_patch_flag[icolor];
  }

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // in: NAVIERSTOKES_3D.F90
  fort_sumdot_ones(
   fab_sum.dataPtr(),
   fab_flag.dataPtr(),
   data_fab.dataPtr(),
   ARLIM(data_fab.loVect()),ARLIM(data_fab.hiVect()),
   ones_fab.dataPtr(),
   ARLIM(ones_fab.loVect()),ARLIM(ones_fab.hiVect()),
   type_fab.dataPtr(),
   ARLIM(type_fab.loVect()),ARLIM(type_fab.hiVect()),
   color_fab.dataPtr(),
   ARLIM(color_fab.loVect()),ARLIM(color_fab.hiVect()),
   alpha_fab.dataPtr(),
   ARLIM(alpha_fab.loVect()),ARLIM(alpha_fab.hiVect()),
   mfab.dataPtr(),
   ARLIM(mfab.loVect()),ARLIM(mfab.hiVect()),
   tilelo,tilehi,
   fablo,fabhi,&bfact,
   &level,
   &gridno,
   &nsolve_expect, 
   presbc.dataPtr(),
   type_ONES_flag.dataPtr(),
   &color_ONES_count,
   &project_option);

  for (int icolor=0;icolor<color_ONES_count;icolor++) {
   local_sum[tid_current][icolor] += fab_sum[icolor];
  } // icolor=0;icolor<color_ONES_count
 } // mfi1
} // omp
 ns_reconcile_d_num(LOOP_SUMDOT_ONES,"dotSumONES");
 for (int tid=1;tid<thread_class::nthreads;tid++) {
  for (int icolor=0;icolor<color_ONES_count;icolor++) {
   local_sum[0][icolor]+=local_sum[tid][icolor];
  } // icolor=0;icolor<color_ONES_count
 }
 for (int icolor=0;icolor<color_ONES_count;icolor++) {
  ParallelDescriptor::ReduceRealSum(local_sum[0][icolor]);
 }

 for (int icolor=0;icolor<color_ONES_count;icolor++) {
  result_sum[icolor]=local_sum[0][icolor];
 }

} // end subroutine dotSumONES

//called from: NavierStokes::mf_combine_ones
void
NavierStokes::mf_combine_ones_level(int project_option,
  int index_MF,
  Vector<Real> beta) {
 
 std::string local_caller_string="mf_combine_ones_level";

 bool use_tiling=ns_tiling;

 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(FACE_VAR_MF+dir,0,local_caller_string); 

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(index_MF,0,local_caller_string);
 debug_ngrow(MASKCOEF_MF,0,local_caller_string);

 debug_ngrow(ONES_MF,0,local_caller_string);
 debug_ngrow(ONES_GROW_MF,0,local_caller_string);
 debug_ngrow(TYPE_ONES_MF,0,local_caller_string); //=1 if diagonal=0.0; =2 if diagonal>0.0
 debug_ngrow(COLOR_ONES_MF,0,local_caller_string);

 if (type_ONES_flag.size()==2) {
  // do nothing
 } else
  amrex::Error("type_ONES_flag.size() invalid");

 debug_ngrow(ALPHACOEF_MF,0,local_caller_string);
 int nsolve_test=localMF[ALPHACOEF_MF]->nComp();
 int nsolve_expect=1;
 if (project_option==SOLVETYPE_VISC) { 
  nsolve_expect=AMREX_SPACEDIM;
 } else if (project_option_is_valid(project_option)==1) {
  nsolve_expect=1;
 } else
  amrex::Error("project_option invalid nsolve_expect");

 if (nsolve_expect==nsolve_test) {
  // do nothing
 } else
  amrex::Error("nsolve_test invalid");

  // for each given color, singular_patch_flag=
  //   0 if color is masked off 
  //   1 if color is not masked off, no compressible/internal dirichlet 
  //     regions, and not touching a Dirichlet condition wall.
  //   2 if color is not masked off, a compressible/internal dirichlet
  //     region exists or color is touching a Dirichlet condition wall.

 int finest_level=parent->finestLevel();
 if (level>finest_level)
  amrex::Error("level too big");

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(localMF[ONES_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel 
#endif
{
 for (MFIter mfi(*localMF[ONES_MF],use_tiling); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());

  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();
  int bfact=parent->Space_blockingFactor(level);

  FArrayBox& data_fab = (*localMF[index_MF])[mfi];
  FArrayBox& ones_fab = (*localMF[ONES_MF])[mfi];
  FArrayBox& type_fab = (*localMF[TYPE_ONES_MF])[mfi];
  FArrayBox& color_fab = (*localMF[COLOR_ONES_MF])[mfi];
  FArrayBox& alpha_fab = (*localMF[ALPHACOEF_MF])[mfi];
  FArrayBox& mfab=(*localMF[MASKCOEF_MF])[mfi];

  Vector<int> presbc=getBCArray(State_Type,gridno,STATECOMP_PRES,1);

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // in: NAVIERSTOKES_3D.F90
  fort_fabcom_ones(
   beta.dataPtr(),
   singular_patch_flag.dataPtr(),
   data_fab.dataPtr(),
   ARLIM(data_fab.loVect()),ARLIM(data_fab.hiVect()),
   ones_fab.dataPtr(),
   ARLIM(ones_fab.loVect()),ARLIM(ones_fab.hiVect()),
   type_fab.dataPtr(),
   ARLIM(type_fab.loVect()),ARLIM(type_fab.hiVect()),
   color_fab.dataPtr(),
   ARLIM(color_fab.loVect()),ARLIM(color_fab.hiVect()),
   alpha_fab.dataPtr(),
   ARLIM(alpha_fab.loVect()),ARLIM(alpha_fab.hiVect()),
   mfab.dataPtr(),
   ARLIM(mfab.loVect()),ARLIM(mfab.hiVect()),
   tilelo,tilehi,
   fablo,fabhi,&bfact,
   &level,
   &gridno,
   &nsolve_expect, 
   presbc.dataPtr(),
   type_ONES_flag.dataPtr(),
   &color_ONES_count,
   &project_option);

 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_FABCOM_ONES,"mf_combine_ones_level");

} // end subroutine mf_combine_ones_level


void NavierStokes::mf_combine(
  int project_option,
  int index_x_MF,
  int index_y_MF,
  Real Beta,
  int index_z_MF,int nsolve) {
  // amf_z = amf_x + Beta amf_y

  int finest_level=parent->finestLevel();
  for (int k = 0; k <= finest_level; ++k) {
    NavierStokes& ns_level = getLevel(k);
    ns_level.levelCombine(
     project_option,
     ns_level.localMF[index_x_MF],
     ns_level.localMF[index_y_MF],
     ns_level.localMF[index_z_MF], Beta,nsolve);
  }

}

  // mfz = mfx + Beta mfy
void NavierStokes::levelCombine(
 int project_option,
 MultiFab* mfx, MultiFab* mfy, MultiFab* mfz,
 Real beta,int nsolve) {
 
 std::string local_caller_string="levelCombine";

 bool use_tiling=ns_tiling;

 if ((nsolve!=1)&&(nsolve!=AMREX_SPACEDIM))
  amrex::Error("nsolve invalid");

 if (project_option_momeqn(project_option)==1) {
  //do nothing
 } else if (project_option_momeqn(project_option)==0) {
  //do nothing
 } else
  amrex::Error("project_option_momeqn(project_option) invalid2");

 int finest_level=parent->finestLevel();
 if (level > finest_level)
  amrex::Error("level too big");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(FACE_VAR_MF+dir,0,local_caller_string);

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,0,local_caller_string);

 if (mfx->nComp()!= nsolve)
  amrex::Error("mfx invalid ncomp");
 if (mfy->nComp()!= nsolve)
  amrex::Error("mfy invalid ncomp");
 if (mfz->nComp()!= nsolve)
  amrex::Error("mfz invalid ncomp");

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(mfx->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel 
#endif
{
 for (MFIter mfi(*mfx,use_tiling); mfi.isValid(); ++mfi) {
  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();
  int bfact=parent->Space_blockingFactor(level);

  FArrayBox& fabx = (*mfx)[mfi];
  FArrayBox& faby = (*mfy)[mfi];
  FArrayBox& mfab = (*localMF[MASKCOEF_MF])[mfi];
  FArrayBox& fabz = (*mfz)[mfi];
  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  //fabz = fabx + beta * faby
  //in: NAVIERSTOKES_3D.F90
  fort_fabcom(
   fabx.dataPtr(),ARLIM(fabx.loVect()),ARLIM(fabx.hiVect()),
   faby.dataPtr(),ARLIM(faby.loVect()),ARLIM(faby.hiVect()),
   mfab.dataPtr(),ARLIM(mfab.loVect()),ARLIM(mfab.hiVect()),
   fabz.dataPtr(),ARLIM(fabz.loVect()),ARLIM(fabz.hiVect()),
   &beta,
   tilelo,tilehi,
   fablo,fabhi,&bfact,
   &nsolve);
 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_FABCOM,"levelCombine");

} // subroutine levelCombine

void NavierStokes::volWgtSum(int isweep,int fast_mode) {
 
 std::string local_caller_string="volWgtSum";

 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();
 NavierStokes& ns_fine = getLevel(finest_level);

 if ((level<=finest_level)&&(level>=0)) {
  // do nothing
 } else
  amrex::Error("level or finest_level invalid");


 if (IQ_TOTAL_SUM_COMP!=NS_sumdata.size())
  amrex::Error("(IQ_TOTAL_SUM_COMP!=NS_sumdata.size())");
 if (IQ_TOTAL_SUM_COMP!=NS_sumdata_type.size())
  amrex::Error("(IQ_TOTAL_SUM_COMP!=NS_sumdata_type.size())");
 if (IQ_TOTAL_SUM_COMP!=NS_sumdata_sweep.size())
  amrex::Error("(IQ_TOTAL_SUM_COMP!=NS_sumdata_sweep.size())");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 if (localMF[CELLTENSOR_MF]->nComp()!=AMREX_SPACEDIM_SQR)
  amrex::Error("localMF[CELLTENSOR_MF]->nComp() invalid");

 MultiFab* den_recon=getStateDen(1,upper_slab_time);  
 int den_ncomp=den_recon->nComp();

 MultiFab* error_heat_map_mf=new MultiFab(grids,dmap,num_materials,0,
	MFInfo().SetTag("error_heat_map_mf"),FArrayBoxFactory());
 error_heat_map_mf -> setVal(0.0);

 int combine_flag=2;
 int hflag=0;
 int combine_idx=-1;  // update state variables
 int update_flux=0;
 int interface_cond_avail=0; // T_I,Y_I not available.

 combine_state_variable(
  SOLVETYPE_VISC,
  combine_idx,
  combine_flag,
  hflag,
  update_flux,
  interface_cond_avail);

 if (fast_mode==0) {

  update_flux=1;
  combine_state_variable(
   SOLVETYPE_PRES,
   combine_idx,
   combine_flag,
   hflag,
   update_flux,
   interface_cond_avail);

 } else if (fast_mode==1) {

  // do nothing
	
 } else
  amrex::Error("fast_mode invalid");

 resize_maskfiner(2,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,2,local_caller_string);

 VOF_Recon_resize(2); //output:SLOPE_RECON_MF
 debug_ngrow(SLOPE_RECON_MF,2,local_caller_string);
 debug_ngrow(CELLTENSOR_MF,1,local_caller_string);

 debug_ngrow(DRAG_MF,ngrow_make_distance,local_caller_string);
 debug_ixType(DRAG_MF,-1,local_caller_string);
 if (localMF[DRAG_MF]->nComp()!=N_DRAG)
  amrex::Error("drag ncomp invalid");

  // velocity and pressure
 MultiFab* vel=getState(1,STATECOMP_VEL,
   STATE_NCOMP_VEL+STATE_NCOMP_PRES,upper_slab_time);

 if (NUM_CELL_ELASTIC==num_materials_viscoelastic*ENUM_NUM_TENSOR_TYPE) {
  //do nothing
 } else
  amrex::Error("NUM_CELL_ELASTIC invalid");

 MultiFab* viscoelastic_tensor=nullptr;

 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=num_materials)) {
  viscoelastic_tensor=getStateTensor(1,0,NUM_CELL_ELASTIC,
      upper_slab_time);
 } else if (num_materials_viscoelastic==0) {
  viscoelastic_tensor=vel;
 } else
  amrex::Error("num_materials_viscoelastic invalid"); 

 const Real* dx = geom.CellSize();

 int resultsize=NS_sumdata.size();
 if (resultsize!=IQ_TOTAL_SUM_COMP)
  amrex::Error("resultsize invalid");
 int num_cells=0;
 int Z_dir=1;
 int R_dir=0;
 const Box& fdomain = ns_fine.geom.Domain();
 const int* fdomlo = fdomain.loVect();
 const int* fdomhi = fdomain.hiVect();

 if (fast_mode==1) {
  // do nothing
 } else if (fast_mode==0) {
  fort_coflow(
   &upper_slab_time, 
   fdomlo,
   fdomhi,
   &Z_dir,
   &R_dir,
   &num_cells,
   NS_coflow_Z.dataPtr(),
   NS_coflow_R_of_Z.dataPtr());
 } else
  amrex::Error("fast_mode invalid");

 if (num_cells+1==NS_coflow_Z.size()) {
  // do nothing
 } else
  amrex::Error("NS_coflow_Z.size() invalid");

 if (num_cells+1==NS_coflow_R_of_Z.size()) {
  // do nothing
 } else
  amrex::Error("NS_coflow_R_of_Z.size() invalid");

 Vector< Vector<Real> > local_result;
 Vector< Vector<Real> > local_coflow_Z;
 Vector< Vector<Real> > local_coflow_R_of_Z;

 local_result.resize(thread_class::nthreads);
 local_coflow_Z.resize(thread_class::nthreads);
 local_coflow_R_of_Z.resize(thread_class::nthreads);

 for (int tid=0;tid<thread_class::nthreads;tid++) {
  local_result[tid].resize(resultsize);

  for (int isum=0;isum<resultsize;isum++) {
   local_result[tid][isum]=0.0;
   if (NS_sumdata_type[isum]==2) // min
    local_result[tid][isum]=1.0E+15;
   else if (NS_sumdata_type[isum]==3)  // max
    local_result[tid][isum]=-1.0E+15;
   else if (NS_sumdata_type[isum]==1)
    local_result[tid][isum]=0.0;
   else
    amrex::Error("NS_sumdata_type invalid");
  } // isum

  local_coflow_Z[tid].resize(num_cells+1);
  local_coflow_R_of_Z[tid].resize(num_cells+1);
  for (int iz=0;iz<=num_cells;iz++) {
   local_coflow_Z[tid][iz]=0.0;
   local_coflow_R_of_Z[tid][iz]=0.0;
  }
 }  // tid

 MultiFab* lsmf=getStateDist(2,upper_slab_time,local_caller_string);

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(localMF[MASKCOEF_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*localMF[MASKCOEF_MF],use_tiling); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   int gridno=mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();

   const Real* xlo = grid_loc[gridno].lo();

   FArrayBox& mfab=(*localMF[MASKCOEF_MF])[mfi];
   FArrayBox& maskSEMfab=(*localMF[MASKSEM_MF])[mfi];
   FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi];
   FArrayBox& denfab=(*den_recon)[mfi];
   FArrayBox& velfab=(*vel)[mfi];
   FArrayBox& viscofab=(*viscoelastic_tensor)[mfi];

   FArrayBox& cellten=(*localMF[CELLTENSOR_MF])[mfi];
   if (cellten.nComp()!=AMREX_SPACEDIM_SQR)
    amrex::Error("cellten invalid ncomp");

   FArrayBox& dragfab=(*localMF[DRAG_MF])[mfi];
   Real problo[AMREX_SPACEDIM];
   Real probhi[AMREX_SPACEDIM];
   for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
    problo[dir]=geom.ProbLo(dir);
    probhi[dir]=geom.ProbHi(dir);
   }
   FArrayBox& lsfab=(*lsmf)[mfi];
   int bfact=parent->Space_blockingFactor(level);

   int local_adapt_quad_depth=adapt_quad_depth;
   if (fast_mode==1) {
    local_adapt_quad_depth=1;
   } else if (fast_mode==0) {
    // do nothing
   } else
    amrex::Error("fast_mode invalid");

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // declared in: NAVIERSTOKES_3D.F90
   fort_summass(
    &tid_current,
    &level,
    &finest_level,
    &ncomp_sum_int_user1,
    &ncomp_sum_int_user2,
    &local_adapt_quad_depth,
    &slice_dir,
    xslice.dataPtr(),
    problo,probhi, 
    xlo,dx,
    cellten.dataPtr(),ARLIM(cellten.loVect()),ARLIM(cellten.hiVect()),
    lsfab.dataPtr(), //lsmf=getStateDist
    ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
    maskSEMfab.dataPtr(),
    ARLIM(maskSEMfab.loVect()),ARLIM(maskSEMfab.hiVect()),
    mfab.dataPtr(),ARLIM(mfab.loVect()),ARLIM(mfab.hiVect()),
    dragfab.dataPtr(),ARLIM(dragfab.loVect()),ARLIM(dragfab.hiVect()),
    reconfab.dataPtr(), //localMF[SLOPE_RECON_MF])
    ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
    denfab.dataPtr(), //den_recon=getStateDen
    ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
    velfab.dataPtr(), //vel=getState
    ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
    viscofab.dataPtr(), //viscoelastic_tensor=getStateTensor
    ARLIM(viscofab.loVect()),ARLIM(viscofab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,
    &bfact,
    &upper_slab_time,
    local_result[tid_current].dataPtr(),
    NS_sumdata.dataPtr(),
    NS_sumdata_type.dataPtr(),
    NS_sumdata_sweep.dataPtr(),
    &resultsize,
    &num_cells,
    local_coflow_Z[tid_current].dataPtr(),
    local_coflow_R_of_Z[tid_current].dataPtr(),
    &Z_dir,
    &R_dir,
    &den_ncomp,
    &isweep);

 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_SUMMASS,"volWgtSum");

 for (int tid=1;tid<thread_class::nthreads;tid++) {

  for (int idest=1;idest<IQ_TOTAL_SUM_COMP;idest++) {

   if (((NS_sumdata_sweep[idest]==0)&&(isweep==0))|| //def (update 1st sweep)
       ((NS_sumdata_sweep[idest]==1)&&(isweep==1))) {//(update 2nd sweep) 

    if (NS_sumdata_type[idest]==1) { // reduce real sum (default)
     local_result[0][idest]+=local_result[tid][idest];
    } else if (NS_sumdata_type[idest]==2) { // reduce real min 
     if (local_result[tid][idest]<local_result[0][idest])
      local_result[0][idest]=local_result[tid][idest];
    } else if (NS_sumdata_type[idest]==3) { // reduce real max
     if (local_result[tid][idest]>local_result[0][idest])
      local_result[0][idest]=local_result[tid][idest];
    } else
     amrex::Error("sumdata_type invalid");

   } else if (NS_sumdata_sweep[idest]==0) {
    // do nothing
   } else if (NS_sumdata_sweep[idest]==1) {
    // do nothing
   } else
    amrex::Error("NS_sumdata_sweep invalid");

  } // idest
 
  if (isweep==0) { 
   for (int iz=0;iz<=num_cells;iz++) {

    if (local_coflow_Z[tid][iz]>local_coflow_Z[0][iz])
     local_coflow_Z[0][iz]=local_coflow_Z[tid][iz];

    if (local_coflow_R_of_Z[tid][iz]>local_coflow_R_of_Z[0][iz])
     local_coflow_R_of_Z[0][iz]=local_coflow_R_of_Z[tid][iz];
   } // iz
  }  // isweep==0
 } // tid

 ParallelDescriptor::Barrier();

 NS_sumdata[IQ_FILLER_SUM_COMP]=0.0;

 for (int idest=1;idest<IQ_TOTAL_SUM_COMP;idest++) {

  if (((NS_sumdata_sweep[idest]==0)&&(isweep==0))||
      ((NS_sumdata_sweep[idest]==1)&&(isweep==1))) { 

   if (NS_sumdata_type[idest]==1) { // reduce real sum
    ParallelDescriptor::ReduceRealSum(local_result[0][idest]);
    NS_sumdata[idest]=NS_sumdata[idest]+local_result[0][idest];
   } else if (NS_sumdata_type[idest]==2) { // reduce real min
    ParallelDescriptor::ReduceRealMin(local_result[0][idest]);
    if (local_result[0][idest]<NS_sumdata[idest])
     NS_sumdata[idest]=local_result[0][idest];
   } else if (NS_sumdata_type[idest]==3) { // reduce real max
    ParallelDescriptor::ReduceRealMax(local_result[0][idest]);
    if (local_result[0][idest]>NS_sumdata[idest])
     NS_sumdata[idest]=local_result[0][idest];
   } else
    amrex::Error("NS_sumdata_type invalid");

  } else if (NS_sumdata_sweep[idest]==0) {
   // do nothing
  } else if (NS_sumdata_sweep[idest]==1) {
   // do nothing
  } else
   amrex::Error("NS_sumdata_sweep invalid");

 } // idest
 
 if (isweep==0) { 
  for (int iz=0;iz<=num_cells;iz++) {
   ParallelDescriptor::ReduceRealMax(local_coflow_Z[0][iz]);
   ParallelDescriptor::ReduceRealMax(local_coflow_R_of_Z[0][iz]);

   if (local_coflow_Z[0][iz]>NS_coflow_Z[iz])
    NS_coflow_Z[iz]=local_coflow_Z[0][iz];

   if (local_coflow_R_of_Z[0][iz]>NS_coflow_R_of_Z[iz])
    NS_coflow_R_of_Z[iz]=local_coflow_R_of_Z[0][iz];
  }
 }  // isweep==0

 if ((fab_verbose==2)||(fab_verbose==3)) {
  std::cout << "c++ level,finest_level " << level << ' ' <<
    finest_level << '\n';

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(error_heat_map_mf->boxArray().d_numPts());

  for (MFIter mfi(*error_heat_map_mf,false); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const Box& tilegrid = mfi.tilebox();

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   const int gridno = mfi.index();
   const Box& fabgrid = grids[gridno];
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();
   const Real* xlo = grid_loc[gridno].lo();
   std::cout << "gridno= " << gridno << '\n';
   std::cout << "output of error heatmap " << '\n';
   int interior_only=0;
   FArrayBox& errfab=(*error_heat_map_mf)[mfi];
   tecplot_debug(errfab,xlo,fablo,fabhi,dx,-1,0,0,
     num_materials,interior_only);
  }// mfi
  ns_reconcile_d_num(LOOP_TECPLOT_DEBUG,"volWgtSum");
 } // fab_verbose=2 or 3

 delete error_heat_map_mf;
 delete lsmf;
 delete den_recon;
 delete vel;

 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=num_materials)) {
  delete viscoelastic_tensor;
 } else if (num_materials_viscoelastic==0) {
  // do nothing
 } else
  amrex::Error("num_materials_viscoelastic invalid"); 

}  // end subroutine volWgtSum

//put "ns.show_mem=1" in the inputs file to activate this.
//called from NavierStokes::do_the_advance
void 
NavierStokes::debug_memory() {

 if (level!=0)
  amrex::Error("level invalid debug_memory");

 if (show_mem==1) {
  std::fflush(NULL);
  ParallelDescriptor::Barrier();
  int proc=ParallelDescriptor::MyProc();
  std::cout << "calling memory status on processor=" << proc << '\n';
  fort_memstatus(&proc);
  std::cout << "after calling memory status on processor=" << proc << '\n';
  std::fflush(NULL);
  ParallelDescriptor::Barrier();
 } else if (show_mem!=0)
  amrex::Error("show_mem invalid");

} // end subroutine debug_memory()

void NavierStokes::writeInterfaceReconstruction() {

 std::string local_caller_string="writeInterfaceReconstruction";

 std::string path1="./temptecplot";
 UtilCreateDirectoryDestructive(path1);

 debug_ngrow(SLOPE_RECON_MF,1,local_caller_string);
 if (level!=0)
  amrex::Error("level should be zero");

 int finest_level = parent->finestLevel();
 Vector<int> grids_per_level;
 grids_per_level.resize(finest_level+1);
 for (int ilev=finest_level;ilev>=0;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  grids_per_level[ilev]=ns_level.grids.size();
   // NavierStokes2.cpp: num_materials materials at once
  ns_level.output_triangles();  
 }
 int nsteps=parent->levelSteps(0);
 int arrdim=finest_level+1;
 int plotint=parent->plotInt();

 ParallelDescriptor::Barrier();
 if (ParallelDescriptor::IOProcessor()) {

  for (int im=1;im<=num_materials;im++) {
    // in: MARCHING_TETRA_3D.F90
   fort_combinetriangles(grids_per_level.dataPtr(),
    &finest_level,
    &nsteps,
    &im,
    &arrdim,
    &cur_time_slab,
    &plotint);
  } // im=1..num_materials

#ifdef AMREX_PARTICLES

  // fort_combine_particles is declared in: NAVIERSTOKES_3D.F90
  // output to: ./PARCON_pos<nsteps>.tec
  fort_combine_particles(
   grids_per_level.dataPtr(),
   &finest_level,
   &nsteps,
   &arrdim,
   &cur_time_slab,
   &plotint);

#endif
 }
 ParallelDescriptor::Barrier();

 std::string path2="./temptecplot";
 UtilCreateDirectoryDestructive(path2);

}  // end subroutine writeInterfaceReconstruction

// in the future, FillPatchUtil should check if 
// new fine grid data gets initialized properly by 
// placing 1's everywhere except the data that needs to be 
// filled.
// 
void NavierStokes::debug_ParallelCopy() {

 ParallelDescriptor::Barrier();
 int s_nghost=0;
 int d_nghost=0;
 int scomp=0;
 int dcomp=0;
 int ncomp=1;
 int src_ncomp=1;
 int dest_ncomp=1;
 int src_ngrow=0;
 int dest_ngrow=0;
 BoxArray src_box_array;
 BoxArray dest_box_array;
 DistributionMapping src_dmap;
 DistributionMapping dest_dmap;
 Vector<int> src_pmap;
 Vector<int> dest_pmap;
//int readFrom (std::istream& is);
//std::ostream& writeOn (std::ostream&) const;
//ofstream myfile; (or ifstream)
//myfile.open("example.txt");
//myfile.close();
 std::cout << "READING FROM PCOPY_DATA_SRC\n";
 std::ifstream is_src;
 is_src.open("PCOPY_DATA_SRC");
 src_box_array.readFrom(is_src);

 std::cout << "READING FROM PCOPY_DATA_DEST\n";
 std::ifstream is_dest;
 is_dest.open("PCOPY_DATA_DEST");
 dest_box_array.readFrom(is_dest);

 std::cout << "src_box_array= " << src_box_array << '\n';
 std::cout << "dest_box_array= " << dest_box_array << '\n';

 src_pmap.resize(src_box_array.size());
 dest_pmap.resize(dest_box_array.size());
 for (int i=0;i<src_pmap.size();i++) {
  src_pmap[i]=0;  
 }
 if (amrex::ParallelDescriptor::NProcs()==4) {
  src_pmap[0]=3;
  src_pmap[1]=1;
  src_pmap[2]=2;
  src_pmap[3]=0;
 }
 for (int i=0;i<dest_pmap.size();i++) {
  dest_pmap[i]=0;
 }
 if (amrex::ParallelDescriptor::NProcs()==4) {
  dest_pmap[0]=1;
  dest_pmap[1]=1;
 }
 src_dmap.define(src_pmap);
 dest_dmap.define(dest_pmap);
 
 std::cout << "src_dmap= " << src_dmap << '\n';
 std::cout << "dest_dmap= " << dest_dmap << '\n';

 MultiFab* src_mf=new MultiFab(src_box_array,src_dmap,src_ncomp,src_ngrow, 
    MFInfo().SetTag("src_mf"),FArrayBoxFactory());
 MultiFab* dest_mf=new MultiFab(dest_box_array,dest_dmap,dest_ncomp,dest_ngrow, 
    MFInfo().SetTag("dest_mf"),FArrayBoxFactory());
 IntVect vec_period(AMREX_SPACEDIM);

  // vec_period[dir]=domain.length(dir) * is_periodic[dir]
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
  vec_period[dir]=0;

 std::cout << "vec_period= " << vec_period << '\n';

 Periodicity my_periodicity(vec_period);
 src_mf->setVal(1.0);
 dest_mf->setVal(1.0);
 ParallelDescriptor::Barrier();
     // scomp,dcomp,ncomp,s_nghost,d_nghost
 dest_mf->ParallelCopy(*src_mf,scomp,dcomp,ncomp,s_nghost,d_nghost,
   my_periodicity);
 ParallelDescriptor::Barrier();
 amrex::Error("debug_ParallelCopy() completed");

} // end subroutine debug_ParallelCopy()

// VOF_Recon_ALL called before this routine is called.
// init_FSI_GHOST_MAC_MF() called for all relevant 
// levels prior to this routine.
void NavierStokes::writeTECPLOT_File(int do_plot,int do_slice) {

 std::string local_caller_string="writeTECPLOT_File";

 std::string path1="./temptecplot";
 UtilCreateDirectoryDestructive(path1);

 if (level!=0)
  amrex::Error("level invalid writeTECPLOT_File");

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid");

 int nsteps=parent->levelSteps(0);


 int finest_level = parent->finestLevel();

 ParallelDescriptor::Barrier();

 debug_ngrow(SLOPE_RECON_MF,1,local_caller_string);
 if (localMF[SLOPE_RECON_MF]->nComp()!=num_materials*ngeom_recon)
  amrex::Error("localMF[SLOPE_RECON_MF]->nComp() invalid");

 int plot_sdim_macro=AMREX_SPACEDIM;

  // uses "slope_recon" 
 if (do_plot==1) {
  writeInterfaceReconstruction();

  if (visual_output_raw_State_Type==1) {
   MultiFab& S_new_temp=get_new_data(State_Type,slab_step+1);
   writeSanityCheckData(
    "RawStateType",
    "RawStateType: vel,pres,den_temp_spec,mofvars,error ind",
    local_caller_string,
    State_Type+GET_NEW_DATA_OFFSET, //tower_mf_id
    S_new_temp.nComp(),
    -1,  //data_mf=-1
    State_Type, //state_type_mf==State_Type
    -1, //data_dir=-1
    nsteps); 
  } else if (visual_output_raw_State_Type==0) {
   // do nothing
  } else
   amrex::Error("visual_output_raw_State_Type invalid");

  if (visual_output_raw_mac_Type==1) {

   for (int dir_mac=0;dir_mac<AMREX_SPACEDIM;dir_mac++) {
    MultiFab& mac_new_temp=get_new_data(Umac_Type+dir_mac,slab_step+1);
    writeSanityCheckData(
     "RawMacType",
     "RawMacType: vel",
     local_caller_string,
     Umac_Type+dir_mac+GET_NEW_DATA_OFFSET, //tower_mf_id
     mac_new_temp.nComp(),
     -1,  //data_mf=-1
     Umac_Type+dir_mac, //state_type_mf
     dir_mac, //data_dir==dir_mac
     nsteps); 
   } //dir_mac=0,..,sdim-1

  } else if (visual_output_raw_mac_Type==0) {
   // do nothing
  } else
   amrex::Error("visual_output_raw_mac_Type invalid");

 } else if (do_plot==0) {
  // do nothing
 } else
  amrex::Error("do_plot invalid");

 Vector<int> grids_per_level_array;
 grids_per_level_array.resize(finest_level+1);
 Vector<BoxArray> cgrids_minusBA_array;
 cgrids_minusBA_array.resize(finest_level+1);

 Vector<Real> slice_data;
 NavierStokes& ns_finest=getLevel(finest_level);
 const Box& domain_finest = ns_finest.geom.Domain();
 const int* domlo_finest = domain_finest.loVect();
 const int* domhi_finest = domain_finest.hiVect();

  // in order to get domain dimensions:
  //geom.ProbLo(dir);
  //geom.ProbHi(dir);
 
 int nslice=domhi_finest[slice_dir]-domlo_finest[slice_dir]+3;
    // x,y,z,xvel,yvel,zvel,PMG,PEOS,DIV,den,Temp,KE
    // (value of material with LS>0)
 int nstate_slice=SLICECOMP_NCOMP;
 slice_data.resize(nslice*nstate_slice);
 for (int i=0;i<nslice*nstate_slice;i++)
  slice_data[i]=-1.0e+30;

  // in: NavierStokes::writeTECPLOT_File
  // (plot_sdim_macro declared up above)
 allocate_array(1,PLOTCOMP_NCOMP,-1,MULTIFAB_TOWER_PLT_MF);

 allocate_levelset_ALL(1,LEVELPC_MF);

// HOLD_VELOCITY_DATA_MF not already allocated,
// so "init_gradu_tensorALL" needs to allocate HOLD_VELOCITY_DATA_MF
// internally for its' own uses then delete it after this call.
 if (localMF_grow[HOLD_VELOCITY_DATA_MF]!=-1)
  amrex::Error("localMF_grow[HOLD_VELOCITY_DATA_MF] invalid");

  //Algorithm for getting the elastic force for plotting purposes:
  //a) save the state velocity
  //b) call vel_elastic_ALL:  u=u+ dt F_elastic
  //c) F_elastic=(current_velocity-saved_velocity)/dt
  //d) restore the saved velocity
 getStateALL(1,cur_time_slab,0,
   AMREX_SPACEDIM,HOLD_VELOCITY_DATA_MF);

 int viscoelastic_force_only=1;
 vel_elastic_ALL(viscoelastic_force_only);

 getStateALL(1,cur_time_slab,STATECOMP_VEL,
   STATE_NCOMP_VEL,CELL_ELASTIC_FORCE_MF);

  //ngrow,scomp,ncomp
 minusALL(1,0,AMREX_SPACEDIM,CELL_ELASTIC_FORCE_MF,HOLD_VELOCITY_DATA_MF);
 if (dt_slab>0.0) {
  Real over_dt=1.0/dt_slab;
  mult_array(1,AMREX_SPACEDIM,over_dt,CELL_ELASTIC_FORCE_MF);
 } else
  amrex::Error("cannot have dt_slab<=0.0 in writeTECPLOT_File");

 Copy_array(GET_NEW_DATA_OFFSET+State_Type,HOLD_VELOCITY_DATA_MF,
   0,STATECOMP_VEL,STATE_NCOMP_VEL,1);
 delete_array(HOLD_VELOCITY_DATA_MF);

 int simple_AMR_BC_flag_viscosity=1;
 int do_alloc=1; 
 init_gradu_tensorALL(
   HOLD_VELOCITY_DATA_MF,//alloc and delete since do_alloc==1
   do_alloc,
   CELLTENSOR_MF,
   FACETENSOR_MF,
   simple_AMR_BC_flag_viscosity);

 if (localMF_grow[HOLD_VELOCITY_DATA_MF]!=-1)
  amrex::Error("localMF_grow[HOLD_VELOCITY_DATA_MF] invalid");

 if (localMF[CELLTENSOR_MF]->nComp()!=AMREX_SPACEDIM_SQR)
  amrex::Error("localMF[CELLTENSOR_MF]->nComp() invalid");
 if (localMF[FACETENSOR_MF]->nComp()!=AMREX_SPACEDIM_SQR)
  amrex::Error("localMF[FACETENSOR_MF]->nComp() invalid");

  //localMF[CELL_VISC_MATERIAL_MF] is deleted in ::Geometry_cleanup()
 getStateVISC_ALL(); //we are in writeTECPLOT_file
 debug_ngrow(CELL_VISC_MATERIAL_MF,1,local_caller_string);
 if (localMF[CELL_VISC_MATERIAL_MF]->nComp()!=3*num_materials)
  amrex::Error("viscmf invalid ncomp");

 debug_ngrow(CELL_CONDUCTIVITY_MATERIAL_MF,1,local_caller_string);
 if (localMF[CELL_CONDUCTIVITY_MATERIAL_MF]->nComp()!=num_materials)
  amrex::Error("thermal_conductivity_data invalid ncomp");

  // getStateDIV_ALL is declared in: MacProj.cpp
  // getStateDIV_ALL computes the derived divergence of the MAC
  // velocity field.
 int idx_source=-1;
 int scomp_src=0;
 int idx_mask=-1;
 getStateDIV_ALL(idx_source,scomp_src,MACDIV_MF,idx_mask);
 if (localMF[MACDIV_MF]->nComp()!=1)
  amrex::Error("localMF[MACDIV_MF]->nComp() invalid");
 if (localMF[MACDIV_MF]->nGrow()!=1)
  amrex::Error("localMF[MACDIV_MF]->nGrow() invalid");

 // if FENE-CR+Carreau,
 // liquid viscosity=etaS+etaP ( 1+ (beta gamma_dot)^alpha )^((n-1)/alpha)
 //
 // for each material, there are 5 components:
 // 1. \dot{gamma}=std::sqrt(2 * D:D)  D=(grad U + grad U^T)/2 (plot label: DT)
 // 2. Tr(A) if viscoelastic (plot label: TR)
 //    \dot{gamma} o.t.
 // 3. Tr(A) (liquid viscosity - etaS)/etaP  if FENE-CR+Carreau
 //    Tr(A) if FENE-CR
 //    \dot{gamma} o.t.
 // 4. (3) * f(A)  if viscoelastic
 //    \dot{gamma} o.t.
 // 5. vorticity magnitude.

 // calls fort_getshear and DERMAGTRACE
 int ntrace=5*num_materials;
 getState_tracemag_ALL(MAGTRACE_MF); //ngrow==1
 if (localMF[MAGTRACE_MF]->nComp()!=ntrace)
  amrex::Error("localMF[MAGTRACE_MF]->nComp() invalid");

  // idx,ngrow,scomp,ncomp,index,scompBC_map
 Vector<int> scompBC_map;
 scompBC_map.resize(1);
 scompBC_map[0]=0; //set_extrap_bc, fort_extrapfill
 PCINTERP_fill_bordersALL(MACDIV_MF,1,0,
   1,State_Type,scompBC_map);

 for (int i=0;i<ntrace;i++) {
  scompBC_map.resize(1);
  scompBC_map[0]=0;
  PCINTERP_fill_bordersALL(MAGTRACE_MF,1,i,1,State_Type,scompBC_map);
 }
 
  // save a copy of the State_Type cell velocity since it will be
  // overwritten by the mass weighted MAC velocity interpolant.
 getStateALL(1,cur_time_slab,STATECOMP_VEL,
   STATE_NCOMP_VEL,HOLD_VELOCITY_DATA_MF);

 int dest_idx=-1; // we put the interpolant in State_Type so that the
                  // command MultiFab* velmf=ns_level.getState( ... 
                  // gets the interpolated data.  We have to restore
                  // HOLD_VELOCITY_DATA_MF at the end.  Note: this should
		  // be done after getStateVISC_ALL() since the WALE model
		  // for eddy viscosity depends on the velocity.
 VELMAC_TO_CELLALL(dest_idx);

 int tecplot_finest_level=finest_level;
 if ((tecplot_max_level<tecplot_finest_level)&&
     (tecplot_max_level>=0)) {
  tecplot_finest_level=tecplot_max_level;
 } else if (tecplot_max_level>=tecplot_finest_level) {
  // do nothing
 } else
  amrex::Error("tecplot_max_level invalid");

 IntVect visual_fab_lo(IntVect::TheZeroVector()); 
 IntVect visual_fab_hi(visual_ncell); 
 Box visual_node_box(visual_fab_lo,visual_fab_hi);
 visual_fab_hi-=IntVect::TheUnitVector();
 Box visual_domain(visual_fab_lo,visual_fab_hi);
  // x,u,p,den,T,Y1..Yn,mag vort,LS
 int visual_ncomp=VISUALCOMP_NCOMP;  
 FArrayBox visual_fab_output(visual_node_box,visual_ncomp);
 FArrayBox visual_fab_input(visual_node_box,visual_ncomp); 

 Array4<Real> const& visual_output_array=visual_fab_output.array();
 Array4<Real> const& visual_input_array=visual_fab_input.array();
 const Dim3 lo3=amrex::lbound(visual_node_box);
 const Dim3 hi3=amrex::ubound(visual_node_box);
 for (int n=0;n<visual_ncomp;++n) {
 for (int z=lo3.z;z<=hi3.z;++z) {
 for (int y=lo3.y;y<=hi3.y;++y) {
 for (int x=lo3.x;x<=hi3.x;++x) {
  visual_output_array(x,y,z,n)=-1.0e+20;
  visual_input_array(x,y,z,n)=-1.0e+20;
 }
 }
 }
 }

 Vector<int> gridlo(3);
 Vector<int> gridhi(3);
 for (int dir=0;dir<3;dir++) {
  gridlo[dir]=0;
  gridhi[dir]=0;
 }
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  gridlo[dir]=0;
  gridhi[dir]=visual_ncell[dir];
 }

 if (visual_compare==1) {

  if (ParallelDescriptor::IOProcessor()) {
   int do_input=1; //read from COARSEDATA.tec and put in visual_fab_input
    //fort_io_compare is declared in: NAVIERSTOKES_3D.F90
   fort_io_compare(
    &nsteps,
    &do_input,
    &visual_compare,
    &cur_time_slab,
    visual_fab_input.dataPtr(),
    ARLIM(visual_fab_input.loVect()), 
    ARLIM(visual_fab_input.hiVect()), 
    visual_fab_output.dataPtr(),
    ARLIM(visual_fab_output.loVect()), 
    ARLIM(visual_fab_output.hiVect()), 
    visual_domain.loVect(),
    visual_domain.hiVect(),
    &visual_ncomp);
  } else if (! ParallelDescriptor::IOProcessor()) {
   // do nothing
  } else
   amrex::Error("ParallelDescriptor::IOProcessor() invalid");

  ParallelDescriptor::Barrier();

   // communicate visual_fab_input from the IO proc to all of the
   // other processors.
  for (int n=0;n<visual_ncomp;n++) {
  for (int k=gridlo[2];k<=gridhi[2];k++) {
  for (int j=gridlo[1];j<=gridhi[1];j++) {
  for (int i=gridlo[0];i<=gridhi[0];i++) {
   Real local_data=visual_input_array(i,j,k,n);
   ParallelDescriptor::ReduceRealMax(local_data);
   ParallelDescriptor::Barrier();
   visual_input_array(i,j,k,n)=local_data;
   ParallelDescriptor::Barrier();
  } // i
  } // j
  } // k
  } // n 

 } else if (visual_compare==0) {
  // do nothing
 } else
  amrex::Error("visual_compare invalid");

 std::string path2="./temptecplot";
 UtilCreateDirectoryDestructive(path2);

 for (int ilev=tecplot_finest_level;ilev>=0;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);

  ns_level.debug_ngrow(MACDIV_MF,1,local_caller_string);

  MultiFab* velmf=ns_level.getState(1,STATECOMP_VEL,
	STATE_NCOMP_VEL+STATE_NCOMP_PRES,cur_time_slab);
  MultiFab* presmf=ns_level.derive_EOS_pressure(material_type_visual); 
  if (presmf->nComp()!=1)
   amrex::Error("presmf has invalid ncomp");

  MultiFab* denmf=ns_level.getStateDen(1,cur_time_slab); 
  ns_level.getStateMOM_DEN(MOM_DEN_MF,1,cur_time_slab); 
  MultiFab* lsdist=ns_level.getStateDist(1,cur_time_slab,local_caller_string);
   //getStateDIV_DATA is declared in: NavierStokes.cpp
   //ngrow=1 scomp=0 ncomp=1  state[DIV_Type]
   //this data goes in the "PLOTCOMP_STATE_DIV" component.
  MultiFab* div_data=ns_level.getStateDIV_DATA(1,0,1,cur_time_slab);
  if (1==0) {
   std::cout << "level= " << ilev << " div_data norm0= " << 
    div_data->norm0() << '\n'; 
   std::cout << "level= " << ilev << " div_datanorm0+1grow= " << 
    div_data->norm0(0,1) << '\n'; 
  }
  MultiFab* viscoelasticmf=nullptr;
  if (NUM_CELL_ELASTIC==num_materials_viscoelastic*ENUM_NUM_TENSOR_TYPE) {
   // do nothing
  } else
   amrex::Error("NUM_CELL_ELASTIC invalid");

  if ((num_materials_viscoelastic>=1)&&
      (num_materials_viscoelastic<=num_materials)) {

   viscoelasticmf=ns_level.getStateTensor(1,0,NUM_CELL_ELASTIC,cur_time_slab);

  } else if (num_materials_viscoelastic==0) {
   viscoelasticmf = lsdist; //placeholder
  } else
   amrex::Error("num_materials_viscoelastic invalid:writeTECPLOT_File");

   //plot_grid_type==0 data interpolated to nodes.
   //plot_grid_type==1 data lives at the cells.
  for (int plot_grid_type=0;plot_grid_type<=1;plot_grid_type++) {

   ParallelDescriptor::Barrier();

   //plot_grid_type==0 data interpolated to nodes.
   //plot_grid_type==1 data lives at the cells.
   ns_level.output_zones(
    plot_grid_type,
    visual_fab_output,
    visual_domain,
    visual_ncomp,
    velmf,
    presmf,
    ns_level.localMF[MACDIV_MF],
    div_data,
    denmf,
    ns_level.localMF[MOM_DEN_MF],
    viscoelasticmf,
    lsdist,
    ns_level.localMF[CELL_VISC_MATERIAL_MF],
    ns_level.localMF[CELL_CONDUCTIVITY_MATERIAL_MF],
    ns_level.localMF[MAGTRACE_MF],
    ns_level.localMF[CELL_ELASTIC_FORCE_MF],
    ns_level.localMF[CELLTENSOR_MF],
    grids_per_level_array[ilev],
    cgrids_minusBA_array[ilev],
    slice_data.dataPtr(), 
    do_plot,do_slice);

   ParallelDescriptor::Barrier();

  } //for (int plot_grid_type=0;plot_grid_type<=1;plot_grid_type++) 

  if ((slice_dir>=0)&&(slice_dir<AMREX_SPACEDIM)) {
   for (int i=0;i<nslice*nstate_slice;i++)
    ParallelDescriptor::ReduceRealMax(slice_data[i]);
  } else
   amrex::Error("slice_dir invalid");

  if ((num_materials_viscoelastic>=1)&&
      (num_materials_viscoelastic<=num_materials)) {
   delete viscoelasticmf;
  } else if (num_materials_viscoelastic==0) {
   // do nothing
  } else
   amrex::Error("num_materials_viscoelastic invalid:writeTECPLOT_File(2)");
   
  delete div_data;
  delete velmf;
  delete denmf;
  delete presmf;
  delete lsdist;
  ns_level.delete_localMF(MOM_DEN_MF,1);
 }  // ilev=tecplot_finest_level ... 0

 ParallelDescriptor::Barrier();

 for (int n=0;n<visual_ncomp;n++) {
 for (int k=gridlo[2];k<=gridhi[2];k++) {
 for (int j=gridlo[1];j<=gridhi[1];j++) {
 for (int i=gridlo[0];i<=gridhi[0];i++) {
   Real local_data=visual_output_array(i,j,k,n);
   ParallelDescriptor::ReduceRealMax(local_data);
   ParallelDescriptor::Barrier();
   visual_output_array(i,j,k,n)=local_data;
   ParallelDescriptor::Barrier();
 } // i
 } // j
 } // k
 } // n 

 int total_number_grids=0;
 for (int ilev=0;ilev<=tecplot_finest_level;ilev++) {
  total_number_grids+=grids_per_level_array[ilev];
 }

 Vector<int> levels_array(total_number_grids);
 Vector<int> bfact_array(total_number_grids);
 Vector<int> gridno_array(total_number_grids);
 Vector<int> gridlo_array(AMREX_SPACEDIM*total_number_grids);
 Vector<int> gridhi_array(AMREX_SPACEDIM*total_number_grids);

 int temp_number_grids=0;
 for (int ilev=0;ilev<=tecplot_finest_level;ilev++) {
  BoxArray cgrids_minusBA;
  cgrids_minusBA=cgrids_minusBA_array[ilev];
  int bfact=parent->Space_blockingFactor(ilev);
  for (int igrid=0;igrid<cgrids_minusBA.size();igrid++) {
   levels_array[temp_number_grids]=ilev; 
   bfact_array[temp_number_grids]=bfact; 
   gridno_array[temp_number_grids]=igrid; 
   const Box& fabgrid = cgrids_minusBA[igrid];
   const int* lo=fabgrid.loVect();
   const int* hi=fabgrid.hiVect();
   for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
    gridlo_array[AMREX_SPACEDIM*temp_number_grids+dir]=lo[dir];
    gridhi_array[AMREX_SPACEDIM*temp_number_grids+dir]=hi[dir];
   }
   temp_number_grids++;
  } // igrid
 } // ilev

 if (temp_number_grids!=total_number_grids)
  amrex::Error("temp_number_grids invalid");
 
 int num_levels=tecplot_finest_level+1;
 int plotint=parent->plotInt();
 int sliceint=parent->sliceInt();

 int nparts=im_solid_map.size();
 if ((nparts<0)||(nparts>num_materials))
  amrex::Error("nparts invalid");
 Vector<int> im_solid_map_null;
 im_solid_map_null.resize(1);
 im_solid_map_null[0]=0;

 int* im_solid_map_ptr;
 int nparts_def=nparts;
 if (nparts==0) {
  im_solid_map_ptr=im_solid_map_null.dataPtr();
  nparts_def=1;
 } else if ((nparts>=1)&&(nparts<=num_materials)) {
  im_solid_map_ptr=im_solid_map.dataPtr();
 } else
  amrex::Error("nparts invalid");

 if (do_plot==1) {
 
   //0=tecplot nodes
  if (visual_nddata_format==0) {

   ParallelDescriptor::Barrier();

   if (ParallelDescriptor::IOProcessor()) {

    fort_combinezones(
     &total_number_grids,
     grids_per_level_array.dataPtr(),
     levels_array.dataPtr(),
     bfact_array.dataPtr(),
     gridno_array.dataPtr(),
     gridlo_array.dataPtr(),
     gridhi_array.dataPtr(),
     &tecplot_finest_level,
     &nsteps,
     &num_levels,
     &cur_time_slab,
     &visual_revolve,
     &plotint,
     &nparts,
     &nparts_def,
     im_solid_map_ptr);

   } else if (!ParallelDescriptor::IOProcessor()) {
    // do nothing
   } else
    amrex::Error("ParallelDescriptor::IOProcessor() corrupt");

   ParallelDescriptor::Barrier();

   //1=plt file cells
  } else if (visual_nddata_format==1) {

   ParallelDescriptor::Barrier();
   std::fflush(NULL);
   if (1==0) {
    std::cout << 
     "getting ready to call WriteMultiLevelPlotfile on proc " << 
       amrex::ParallelDescriptor::MyProc() << "\n";
    std::cout <<
     "tecplot_max_level= " << tecplot_max_level << "\n";
   }
   std::fflush(NULL);
   ParallelDescriptor::Barrier();

   int ncomp_plot_MOF=num_materials*ngeom_recon;
   debug_ngrow(SLOPE_RECON_MF,1,local_caller_string);
   if (localMF[SLOPE_RECON_MF]->nComp()!=ncomp_plot_MOF)
    amrex::Error("localMF[SLOPE_RECON_MF]->nComp() invalid");

   Vector<const MultiFab*> mf_tower_MOF;
   mf_tower_MOF.resize(tecplot_finest_level+1);
   Vector<std::string> varnames_MOF;
   varnames_MOF.resize(ncomp_plot_MOF);

   const Vector<Geometry>& ns_geom=parent->Geom();
   Vector<IntVect> ref_ratio;
   ref_ratio.resize(tecplot_finest_level+1);
   Vector<int> level_steps;
   level_steps.resize(tecplot_finest_level+1);

   for (int ilev=0;ilev<=tecplot_finest_level;ilev++) {
    NavierStokes& ns_level=getLevel(ilev);
    mf_tower_MOF[ilev]=ns_level.localMF[SLOPE_RECON_MF];
   }

   for (int ilev=0;ilev<=tecplot_finest_level;ilev++) {
    for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
     ref_ratio[ilev][dir]=2;
    }
    level_steps[ilev]=nsteps;
   } 

   std::stringstream steps_string_stream(std::stringstream::in |
     std::stringstream::out);
   steps_string_stream << std::setw(8) << std::setfill('0') << nsteps;
   std::string steps_string=steps_string_stream.str();

   std::string plotfilename_MOF="MOF_PLT"; 
   plotfilename_MOF+=steps_string;

   int icomp_MOF=0;

    // vfrac,centroid,order,slope,intercept x num_materials
   for (int im=0;im<num_materials;im++) {
    std::stringstream im_string_stream(std::stringstream::in |
     std::stringstream::out);
    im_string_stream << std::setw(2) << std::setfill('0') << im+1;
    std::string im_string=im_string_stream.str();
    varnames_MOF[icomp_MOF]="F"+im_string;
    icomp_MOF++;
    varnames_MOF[icomp_MOF]="XCEN"+im_string;
    icomp_MOF++;
    varnames_MOF[icomp_MOF]="YCEN"+im_string;
    icomp_MOF++;
    if (AMREX_SPACEDIM==3) {
     varnames_MOF[icomp_MOF]="ZCEN"+im_string;
     icomp_MOF++;
    }
    varnames_MOF[icomp_MOF]="Order"+im_string;
    icomp_MOF++;
    varnames_MOF[icomp_MOF]="XSLOPE"+im_string;
    icomp_MOF++;
    varnames_MOF[icomp_MOF]="YSLOPE"+im_string;
    icomp_MOF++;
    if (AMREX_SPACEDIM==3) {
     varnames_MOF[icomp_MOF]="ZSLOPE"+im_string;
     icomp_MOF++;
    }
    varnames_MOF[icomp_MOF]="INTERCEPT"+im_string;
    icomp_MOF++;
   }

   if (icomp_MOF==ncomp_plot_MOF) {
    // do nothing
   } else
    amrex::Error("icomp_MOF!=ncomp_plot_MOF");

#ifdef AMREX_USE_HDF5
   WriteMultiLevelPlotfileHDF5(plotfilename_MOF,
     tecplot_finest_level+1, //nlevels
     mf_tower_MOF,
     varnames_MOF,
     ns_geom,   
     cur_time_slab,
     level_steps,
     ref_ratio);
#else
   WriteMultiLevelPlotfile(plotfilename_MOF,
     tecplot_finest_level+1, //nlevels
     mf_tower_MOF,
     varnames_MOF,
     ns_geom,   
     cur_time_slab,
     level_steps,
     ref_ratio);
#endif

   ParallelDescriptor::Barrier();

   int ncomp_plot=PLOTCOMP_NCOMP;
   if (localMF[MULTIFAB_TOWER_PLT_MF]->nComp()==ncomp_plot) {
    // do nothing
   } else
    amrex::Error("localMF[MULTIFAB_TOWER_PLT_MF]->nComp()!=PLOTCOMP_NCOMP");

   Vector<const MultiFab*> mf_tower;
   mf_tower.resize(tecplot_finest_level+1);
   Vector<std::string> varnames;
   varnames.resize(ncomp_plot);

   for (int ilev=0;ilev<=tecplot_finest_level;ilev++) {
    NavierStokes& ns_level=getLevel(ilev);
    mf_tower[ilev]=ns_level.localMF[MULTIFAB_TOWER_PLT_MF];
   }

   std::string plotfilename="nddataPLT"; 
   plotfilename+=steps_string;

   int icomp=0;
   varnames[icomp]="X";
   icomp++;
   varnames[icomp]="Y";
   if (AMREX_SPACEDIM==3) {
    icomp++;
    varnames[icomp]="Z";
   }
   icomp++;
   varnames[icomp]="x_velocity";
   icomp++;
   varnames[icomp]="y_velocity";
   if (AMREX_SPACEDIM==3) {
    icomp++;
    varnames[icomp]="z_velocity";
   }
   icomp++;
   varnames[icomp]="PRES_MG";
   icomp++;
   varnames[icomp]="PRES_EOS";
   icomp++;
   varnames[icomp]="DIV_DERIVED";
   icomp++;
   varnames[icomp]="DIV_EXPECT";
   icomp++;
   varnames[icomp]="MACH";

   for (int im=0;im<num_materials;im++) {
    std::stringstream im_string_stream(std::stringstream::in |
     std::stringstream::out);
    im_string_stream << std::setw(2) << std::setfill('0') << im+1;
    std::string im_string=im_string_stream.str();
    icomp++;
    varnames[icomp]="F"+im_string;
   }

   for (int im=0;im<num_materials;im++) {
    std::stringstream im_string_stream(std::stringstream::in |
     std::stringstream::out);
    im_string_stream << std::setw(2) << std::setfill('0') << im+1;
    std::string im_string=im_string_stream.str();
    std::string im_opp_string=im_string_stream.str();
    icomp++;
    varnames[icomp]="L"+im_string+im_opp_string;
   }

   for (int im=0;im<num_materials;im++) {
    std::stringstream im_string_stream(std::stringstream::in |
     std::stringstream::out);
    im_string_stream << std::setw(2) << std::setfill('0') << im+1;
    std::string im_string=im_string_stream.str();
    std::string im_opp_string=im_string_stream.str();
    icomp++;
    varnames[icomp]="x_normal"+im_string+im_opp_string;
    icomp++;
    varnames[icomp]="y_normal"+im_string+im_opp_string;
    if (AMREX_SPACEDIM==3) {
     icomp++;
     varnames[icomp]="z_normal"+im_string+im_opp_string;
    }
   }

   if (icomp+1==PLOTCOMP_SCALARS) {
    // do nothing
   } else
    amrex::Error("icomp+1!=PLOTCOMP_SCALARS");

   for (int im=0;im<num_materials;im++) {
    std::stringstream im_string_stream(std::stringstream::in |
     std::stringstream::out);
    im_string_stream << std::setw(2) << std::setfill('0') << im+1;
    std::string im_string=im_string_stream.str();
    icomp++;
    varnames[icomp]="D"+im_string;
    icomp++;
    varnames[icomp]="T"+im_string;

    for (int ispec=0;ispec<num_species_var;ispec++) {
     std::stringstream ispec_string_stream(std::stringstream::in |
      std::stringstream::out);
     ispec_string_stream << std::setw(2) << std::setfill('0') << ispec+1;
     std::string ispec_string=ispec_string_stream.str();
     icomp++;
     varnames[icomp]="S"+ispec_string+"-"+im_string;
    } //ispec=0;ispec<num_species_var;ispec++

   } //im=0..num_materials-1

   if (icomp+1==PLOTCOMP_SCALARS_MERGE) {
    // do nothing
   } else
    amrex::Error("icomp+1!=PLOTCOMP_SCALARS_MERGE");

   icomp++;
   varnames[icomp]="DMERGE";
   icomp++;
   varnames[icomp]="TMERGE";

   for (int ispec=0;ispec<num_species_var;ispec++) {
    std::stringstream ispec_string_stream(std::stringstream::in |
      std::stringstream::out);
    ispec_string_stream << std::setw(2) << std::setfill('0') << ispec+1;
    std::string ispec_string=ispec_string_stream.str();
    icomp++;
    varnames[icomp]="SMERGE"+ispec_string;
   }  // for (int ispec=0;ispec<num_species_var;ispec++)

   if (icomp+1==PLOTCOMP_MOMDEN) {
    // do nothing
   } else
    amrex::Error("icomp+1!=PLOTCOMP_MOMDEN");

   for (int im=0;im<num_materials;im++) {
    std::stringstream im_string_stream(std::stringstream::in |
     std::stringstream::out);
    im_string_stream << std::setw(2) << std::setfill('0') << im+1;
    std::string im_string=im_string_stream.str();
    icomp++;
    varnames[icomp]="MOMDEN"+im_string;
   }  // for (int im=0;im<num_materials;im++)

   for (int partid=0;partid<num_materials_viscoelastic;partid++) {
    int im=im_elastic_map[partid];
    if ((im>=0)&&(im<num_materials)) {
     std::stringstream im_string_stream(std::stringstream::in |
      std::stringstream::out);
     im_string_stream << std::setw(2) << std::setfill('0') << im+1;
     std::string im_string=im_string_stream.str();
     for (int ispec=0;ispec<ENUM_NUM_TENSOR_TYPE;ispec++) {
      std::stringstream ispec_string_stream(std::stringstream::in |
       std::stringstream::out);
      ispec_string_stream << std::setw(2) << std::setfill('0') << ispec+1;
      std::string ispec_string=ispec_string_stream.str();
      icomp++;
      varnames[icomp]="CT"+ispec_string+"-"+im_string;
     } //ispec
    } else
     amrex::Error("im invalid");
   } //partid

   if (icomp+1==PLOTCOMP_VISC) {
    // do nothing
   } else
    amrex::Error("icomp+1!=PLOTCOMP_VISC");

   for (int im=0;im<num_materials;im++) {
    std::stringstream im_string_stream(std::stringstream::in |
     std::stringstream::out);
    im_string_stream << std::setw(2) << std::setfill('0') << im+1;
    std::string im_string=im_string_stream.str();
    icomp++;
    varnames[icomp]="MU"+im_string;
   }
   for (int im=0;im<num_materials;im++) {
    std::stringstream im_string_stream(std::stringstream::in |
     std::stringstream::out);
    im_string_stream << std::setw(2) << std::setfill('0') << im+1;
    std::string im_string=im_string_stream.str();
    icomp++;
    varnames[icomp]="K_THERMAL"+im_string;
   }
   for (int im=0;im<num_materials;im++) {
    std::stringstream im_string_stream(std::stringstream::in |
     std::stringstream::out);
    im_string_stream << std::setw(2) << std::setfill('0') << im+1;
    std::string im_string=im_string_stream.str();
    icomp++;
    varnames[icomp]="DT"+im_string;
    icomp++;
    varnames[icomp]="TR"+im_string;
    icomp++;
    varnames[icomp]="TRT"+im_string;
    icomp++;
    varnames[icomp]="TRTF"+im_string;
    icomp++;
    varnames[icomp]="VORT"+im_string;
   }

   if (icomp+1==PLOTCOMP_F_ELASTIC_X) {
    // do nothing
   } else
    amrex::Error("icomp+1!=PLOTCOMP_F_ELASTIC_X");

   icomp++;
   varnames[icomp]="X_ELSTCFORCE";
   icomp++;
   varnames[icomp]="Y_ELSTCFORCE";
   if (AMREX_SPACEDIM==3) {
    icomp++;
    varnames[icomp]="Z_ELSTCFORCE";
   }

   if (icomp+1==PLOTCOMP_GRAD_VELOCITY) {
    // do nothing
   } else
    amrex::Error("icomp+1!=PLOTCOMP_GRAD_VELOCITY");

   for (int igrad=0;igrad<AMREX_SPACEDIM_SQR;igrad++) {
    std::stringstream igrad_string_stream(std::stringstream::in |
     std::stringstream::out);
    igrad_string_stream << std::setw(2) << std::setfill('0') << igrad+1;
    std::string igrad_string=igrad_string_stream.str();
    icomp++;
    varnames[icomp]="GRADVEL"+igrad_string;
   }

   if (icomp+1==PLOTCOMP_NCOMP) {
    // do nothing
   } else
    amrex::Error("icomp+1!=PLOTCOMP_NCOMP");

#ifdef AMREX_USE_HDF5
   WriteMultiLevelPlotfileHDF5(plotfilename,
     tecplot_finest_level+1, //nlevels
     mf_tower,
     varnames,
     ns_geom,   
     cur_time_slab,
     level_steps,
     ref_ratio);
#else
   WriteMultiLevelPlotfile(plotfilename,
     tecplot_finest_level+1, //nlevels
     mf_tower,
     varnames,
     ns_geom,   
     cur_time_slab,
     level_steps,
     ref_ratio);
#endif

   ParallelDescriptor::Barrier();

   //2=tecplot cells (piecewise constant reconstruction).
  } else if (visual_nddata_format==2) {

   ParallelDescriptor::Barrier();

    //plot_sdim_macro declared up above.
   int ncomp_plot=PLOTCOMP_NCOMP;
   if (localMF[MULTIFAB_TOWER_PLT_MF]->nComp()==ncomp_plot) {
    // do nothing
   } else
    amrex::Error("localMF[MULTIFAB_TOWER_PLT_MF]->nComp()!=PLOTCOMP_NCOMP");

   int data_dir=-1;
   writeSanityCheckData(
    "nddata_piecewise_const",
    "nddata_piecewise_const",
    local_caller_string,
    MULTIFAB_TOWER_PLT_MF, //tower_mf_id
    ncomp_plot,
    MULTIFAB_TOWER_PLT_MF,
    -1, //State_Type==-1
    data_dir,
    nsteps);

   ParallelDescriptor::Barrier();

  } else {
   amrex::Error("visual_nddata_format invalid");
  }
 } else if (do_plot==0) {
  // do nothing
 } else
  amrex::Error("do_plot invalid");

 if (do_slice==1) {
  if ((slice_dir>=0)&&(slice_dir<AMREX_SPACEDIM)) {
   if (ParallelDescriptor::IOProcessor()) {
    fort_outputslice(&cur_time_slab,&nsteps,&sliceint,
     slice_data.dataPtr(),&nslice,&nstate_slice);
   } else if (!ParallelDescriptor::IOProcessor()) {
    // do nothing
   } else
    amrex::Error("ParallelDescriptor::IOProcessor() corrupt");
   
  } else
   amrex::Error("slice_dir invalid");
 } else if (do_slice==0) {
  // do nothing
 } else
  amrex::Error("do_slice invalid");

  // in: NAVIERSTOKES_3D.F90
  // just output data, no comparison with other data.
  // "visual_compare" unused here since do_input==0.
  // the data file name is uniform??????.tec
  // X,Y,Z,U,V,W,pres,density,T,Y1,..,Yn,MGVORT,LS01,...,LS_m
 int do_input=0; 

 if (ParallelDescriptor::IOProcessor()) {

   //fort_io_compare is declared in: NAVIERSTOKES_3D.F90
  fort_io_compare(
   &nsteps,
   &do_input,
   &visual_compare,//if visual_compare==1 then compare fab_input to fab_output
   &cur_time_slab,
   visual_fab_input.dataPtr(), //initialized from COARSEDATA.tec
   ARLIM(visual_fab_input.loVect()), 
   ARLIM(visual_fab_input.hiVect()), 
   visual_fab_output.dataPtr(), //initialized by "ns_level.output_zones"
   ARLIM(visual_fab_output.loVect()), 
   ARLIM(visual_fab_output.hiVect()), 
   visual_domain.loVect(),
   visual_domain.hiVect(),
   &visual_ncomp);

 } else if (!ParallelDescriptor::IOProcessor()) {
  // do nothing
 } else
  amrex::Error("ParallelDescriptor::IOProcessor() corrupt");

 ParallelDescriptor::Barrier();

 Copy_array(GET_NEW_DATA_OFFSET+State_Type,HOLD_VELOCITY_DATA_MF,
   0,STATECOMP_VEL,STATE_NCOMP_VEL,1);
 delete_array(HOLD_VELOCITY_DATA_MF);

 delete_array(MULTIFAB_TOWER_PLT_MF);

 delete_array(MACDIV_MF);
 delete_array(MAGTRACE_MF); 
 delete_array(CELL_ELASTIC_FORCE_MF); 
 delete_array(CELLTENSOR_MF);
 delete_array(FACETENSOR_MF);

 std::string path3="./temptecplot";
 UtilCreateDirectoryDestructive(path3);

} // subroutine writeTECPLOT_File



void NavierStokes::writeSanityCheckData(
	const std::string& root_string,
	const std::string& information_string,
	const std::string& caller_string,
	int tower_mf_id,
        int ncomp,
        int data_mf, 
	int state_type_mf,
        int data_dir,
	int nsteps_actual) {

 if (tower_mf_id>=0) {
  //do nothing
 } else if (tower_mf_id-GET_NEW_DATA_OFFSET>=0) {
   //do nothing
 } else
  amrex::Error("tower_mf_id out of range");

 std::string path1="./temptecplot";
 UtilCreateDirectoryDestructive(path1);

 if (ParallelDescriptor::IOProcessor()) {
  std::cout << "in: writeSanityCheckData, root_string= " <<
    root_string << '\n';
  std::cout << "in: writeSanityCheckData, information_string= " <<
    information_string << '\n';
  std::cout << "in: writeSanityCheckData, caller_string= " <<
    caller_string << '\n';
  std::cout << "in: writeSanityCheckData, data_mf= " <<
    data_mf << " state_type_mf=" << state_type_mf << '\n';
  std::cout << "in: writeSanityCheckData, tower_mf_id= " <<
    tower_mf_id << '\n';
  std::cout << "in: writeSanityCheckData, tower_mf_id-ofs= " <<
    tower_mf_id-GET_NEW_DATA_OFFSET << '\n';
 }

 if (nsteps_actual>=0) {
  // do nothing
 } else
  amrex::Error("nsteps_actual invalid");

 if (level!=0)
  amrex::Error("level invalid writeSanityCheckData");

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid");

 int finest_level = parent->finestLevel();

 ParallelDescriptor::Barrier();

 Vector<int> grids_per_level_array;
 grids_per_level_array.resize(finest_level+1);
 Vector<BoxArray> cgrids_minusBA_array;
 cgrids_minusBA_array.resize(finest_level+1);

 MultiFab* raw_data_lev0_mf=nullptr;

 if (data_mf>=0) {
  if (state_type_mf==-1) {
   // do nothing
  } else
   amrex::Error("state_type_mf invalid");

  if (localMF_grow[data_mf]>=0) {
   // do nothing
  } else
   amrex::Error("localMF_grow[data_mf] invalid");

  raw_data_lev0_mf=localMF[data_mf];

 } else if (data_mf==-1) {

  if (state_type_mf>=0) {
   raw_data_lev0_mf=&get_new_data(state_type_mf,slab_step+1);
  } else {
   raw_data_lev0_mf=nullptr;
   amrex::Error("state_type_mf invalid");
  }

  if (raw_data_lev0_mf->nGrow()>=0) {
   // do nothing
  } else
   amrex::Error("raw_data_lev0_mf->nGrow() invalid");

 } else
  amrex::Error("data_mf invalid");

 int tecplot_finest_level=finest_level;

 if ((tecplot_max_level<tecplot_finest_level)&&
     (tecplot_max_level>=0)) {
  tecplot_finest_level=tecplot_max_level;
 } else if (tecplot_max_level>=tecplot_finest_level) {
  // do nothing
 } else
  amrex::Error("tecplot_max_level invalid");

 for (int ilev=tecplot_finest_level;ilev>=0;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);

  MultiFab* raw_data_mf=nullptr;

  if (data_mf>=0) {
   if (state_type_mf==-1) {
    // do nothing
   } else
    amrex::Error("state_type_mf invalid");

   if (ns_level.localMF_grow[data_mf]>=0) {
    // do nothing
   } else
    amrex::Error("ns_level.localMF_grow[data_mf] invalid");

   raw_data_mf=ns_level.localMF[data_mf];
   ns_level.debug_ngrow(data_mf,0,caller_string);

  } else if (data_mf==-1) {

   if (state_type_mf>=0) {
    raw_data_mf=&ns_level.get_new_data(state_type_mf,slab_step+1);
   } else {
    raw_data_mf=nullptr;
    amrex::Error("state_type_mf invalid");
   }

   if (raw_data_mf->nGrow()>=0) {
    // do nothing
   } else
    amrex::Error("raw_data_mf->nGrow() invalid");

  } else
   amrex::Error("data_mf invalid");

   // data_dir=-1 cell centered data
   // data_dir=0..sdim-1 face centered data.
   // data_dir=3  X,Y node
   // data_dir=4  X,Z node
   // data_dir=5  Y,Z node
  ns_level.Sanity_output_zones(
   information_string, 
   tower_mf_id,
   data_dir,
   raw_data_mf,
   ncomp,
   grids_per_level_array[ilev],
   cgrids_minusBA_array[ilev]);

 }  // ilev=tecplot_finest_level ... 0

 ParallelDescriptor::Barrier();

 int total_number_grids=0;
 for (int ilev=0;ilev<=tecplot_finest_level;ilev++) {
  total_number_grids+=grids_per_level_array[ilev];
 }

 Vector<int> levels_array(total_number_grids);
 Vector<int> bfact_array(total_number_grids);
 Vector<int> gridno_array(total_number_grids);
 Vector<int> gridlo_array(AMREX_SPACEDIM*total_number_grids);
 Vector<int> gridhi_array(AMREX_SPACEDIM*total_number_grids);

 int temp_number_grids=0;
 for (int ilev=0;ilev<=tecplot_finest_level;ilev++) {
  BoxArray cgrids_minusBA;
  cgrids_minusBA=cgrids_minusBA_array[ilev];
  int bfact=parent->Space_blockingFactor(ilev);
  for (int igrid=0;igrid<cgrids_minusBA.size();igrid++) {
   levels_array[temp_number_grids]=ilev; 
   bfact_array[temp_number_grids]=bfact; 
   gridno_array[temp_number_grids]=igrid; 
   const Box& fabgrid = cgrids_minusBA[igrid];
   const int* lo=fabgrid.loVect();
   const int* hi=fabgrid.hiVect();
   for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
    gridlo_array[AMREX_SPACEDIM*temp_number_grids+dir]=lo[dir];
    gridhi_array[AMREX_SPACEDIM*temp_number_grids+dir]=hi[dir];
   }
   temp_number_grids++;
  } // igrid
 } // ilev=0...tecplot_finest_level

 int n_root=root_string.length();
 Vector<char> root_char_array;
 root_char_array.resize(n_root+1);
 for (int i_string=0;i_string<n_root;i_string++)
  root_char_array[i_string]=root_string[i_string];

 if (temp_number_grids!=total_number_grids)
  amrex::Error("temp_number_grids invalid");
 
 int num_levels=tecplot_finest_level+1;

 ParallelDescriptor::Barrier();

 if (ParallelDescriptor::IOProcessor()) {

  // declared in: TECPLOTUTIL.F90
  fort_combinezones_sanity(
   root_char_array.dataPtr(),
   &n_root,
   &data_dir,
   &total_number_grids,
   grids_per_level_array.dataPtr(),
   levels_array.dataPtr(),
   bfact_array.dataPtr(),
   gridno_array.dataPtr(),
   gridlo_array.dataPtr(),
   gridhi_array.dataPtr(),
   &tecplot_finest_level,
   &SDC_outer_sweeps,
   &slab_step,
   &tower_mf_id,
   &nsteps_actual,
   &num_levels,
   &cur_time_slab,
   &visual_revolve,
   &ncomp);

 } else if (!ParallelDescriptor::IOProcessor()) {
  // do nothing
 } else
  amrex::Error("ParallelDescriptor::IOProcessor() corrupt");

 ParallelDescriptor::Barrier();

 std::string path2="./temptecplot";
 UtilCreateDirectoryDestructive(path2);

} // end subroutine writeSanityCheckData

// 1. coarseTimeStep
// 2. timeStep
//     a. advance
//         i. CopyNewToOld
//         ii. setTimeLevel(time+dt_AMR,dt_AMR)
//     b. level_steps++
// 3. cumtime += dt_AMR
// 4. if (level_steps[0] % plot_int == 0)
//     writePlotFile
void
NavierStokes::writePlotFile (
  int do_plot,int do_slice,
  int SDC_outer_sweeps_in,
  int slab_step_in,
  int divu_outer_sweeps_in) {

 std::string path1="./temptecplot";
 UtilCreateDirectoryDestructive(path1);

 std::string local_caller_string="writePlotFile";

 SDC_setup();
 ns_time_order=parent->Time_blockingFactor();

 SDC_outer_sweeps=SDC_outer_sweeps_in;
 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid in writePlotFile");

 slab_step=slab_step_in; 

 SDC_setup_step();

 divu_outer_sweeps=divu_outer_sweeps_in;
 if ((divu_outer_sweeps>=0)&&
     (divu_outer_sweeps<num_divu_outer_sweeps)) {
  // do nothing
 } else
  amrex::Error("divu_outer_sweeps invalid in writePlotFile");

  // metrics_dataALL
  // MASKCOEF_MF
  // MASK_NBR_MF
  // init_FSI_GHOST_MAC_MF_ALL
  // LEVELPC_MF
  // MASKSEM_MF
  // VOF_Recon_ALL
  // make_physics_varsALL
 if (level==0) {
  prepare_post_process(local_caller_string);
 } // level==0

  // output tecplot zonal files  x,y,z,u,v,w,phi,psi

 if (level==0) {
  writeTECPLOT_File(do_plot,do_slice);
 }

 std::string path2="./temptecplot";
 UtilCreateDirectoryDestructive(path2);

} // end subroutine writePlotFile

void NavierStokes::DumpProcNum() {

 
 bool use_tiling=ns_tiling;
 MultiFab& unew = get_new_data(State_Type,slab_step+1);

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(unew.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(unew,use_tiling); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  const Box& tilegrid = mfi.tilebox();

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  const int i = mfi.index();
  if (verbose>0)
   std::cout << "level=" << level << " grid=" << i << " proc=" <<
    ParallelDescriptor::MyProc() << " thread= " << ns_thread() << "\n";
 }// mfi
}// omp
 ns_reconcile_d_num(LOOP_DUMPPROCNUM,"DumpProcNum");
}

// called from: estTimeStep 
//              sum_integrated_quantities 
void NavierStokes::MaxAdvectSpeedALL(
  Real& dt_min,
  Vector<Real>& vel_max_estdt,
  const std::string& caller_string) {

 if (vel_max_estdt.size()!=AMREX_SPACEDIM+1)
  amrex::Error("vel_max_estdt has invalid size");

 std::string local_caller_string="MaxAdvectSpeedALL";
 local_caller_string=caller_string+local_caller_string;

 if (pattern_test(local_caller_string,"estTimeStep")==1) {
  // do nothing
 } else if (pattern_test(local_caller_string,"sum_integrated_quantities")==1) {
  // do nothing
 } else
  amrex::Error("local_caller_string invalid");

 int finest_level=parent->finestLevel();
 if (level!=0)
  amrex::Error("level invalid MaxAdvectSpeedALL");

  // last component is max|c|^2
 Vector<Real> local_vel_max_estdt(AMREX_SPACEDIM+1);  

 Real local_dt_min;

 for (int dir=0;dir<AMREX_SPACEDIM+1;dir++) {
  vel_max_estdt[dir]=0.0;
 }

 dt_min=1.0E+30;
 if (dt_min<dt_max) 
  dt_min=dt_max;

 if (localMF_grow[FSI_GHOST_MAC_MF]==-1) {
   int filler_renormalize_only=1;
   init_FSI_GHOST_MAC_MF_ALL(filler_renormalize_only,local_caller_string);
 } else if (localMF_grow[FSI_GHOST_MAC_MF]>=0) {
  // do nothing
 } else {
  amrex::Error("localMF_grow[FSI_GHOST_MAC_MF] invalid");
 }

 for (int ilev=finest_level;ilev>=0;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);

  for (int dir=0;dir<AMREX_SPACEDIM+1;dir++) {
   local_vel_max_estdt[dir]=0.0;
  }

  local_dt_min=1.0E+30;
  if (local_dt_min<dt_max) 
   local_dt_min=dt_max;
  
  ns_level.MaxAdvectSpeed(
    local_dt_min,
    local_vel_max_estdt,
    local_caller_string); 

  for (int dir=0;dir<AMREX_SPACEDIM+1;dir++) {
   vel_max_estdt[dir] = std::max(vel_max_estdt[dir],local_vel_max_estdt[dir]);
  }
  dt_min=std::min(dt_min,local_dt_min);

 } // ilev = finest_level downto 0

} // end subroutine MaxAdvectSpeedALL

// MaxAdvectSpeedALL called from: 
//     estTimeStep 
//     sum_integrated_quantities
// vel_max_estdt= max vel in direction + extra cell centered velocity
//   considerations if spectral element method.
//   vel_max_estdt[sdim]=max c^2
void NavierStokes::MaxAdvectSpeed(
 Real& dt_min, // minimum for the given level.
 Vector<Real>& vel_max_estdt, // maximum for the given level.
 const std::string& caller_string) {

 if (vel_max_estdt.size()!=AMREX_SPACEDIM+1)
  amrex::Error("vel_max_estdt has invalid size");

 std::string local_caller_string="MaxAdvectSpeed";
 local_caller_string=caller_string+local_caller_string;

 if (pattern_test(local_caller_string,"estTimeStep")==1) {
  // do nothing
 } else if (pattern_test(local_caller_string,"sum_integrated_quantities")==1) {
  // do nothing
 } else
  amrex::Error("local_caller_string invalid");

 int finest_level=parent->finestLevel();

 int nparts=im_solid_map.size();
 if ((nparts<0)||(nparts>num_materials))
  amrex::Error("nparts invalid");
 Vector<int> im_solid_map_null;
 im_solid_map_null.resize(1);

 int* im_solid_map_ptr;
 int nparts_def=nparts;
 if (nparts==0) {
  im_solid_map_ptr=im_solid_map_null.dataPtr();
  nparts_def=1;
 } else if ((nparts>=1)&&(nparts<=num_materials)) {
  im_solid_map_ptr=im_solid_map.dataPtr();
 } else
  amrex::Error("nparts invalid");

 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
   if (localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0)
    amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0");
   if (localMF[FSI_GHOST_MAC_MF+data_dir]->nComp()!=nparts_def*AMREX_SPACEDIM)
    amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nComp() invalid");
 }

 MultiFab* distmf=getStateDist(2,cur_time_slab,local_caller_string);
  // num_materials*num_state_material
 MultiFab* denmf=getStateDen(1,cur_time_slab);  
 MultiFab* vofmf=
    getState(1,STATECOMP_MOF,num_materials*ngeom_raw,cur_time_slab);

 const Real* dx = geom.CellSize();

 dt_min=1.0E+30;
 if (dt_min<dt_max) 
  dt_min=dt_max;

 for (int dir=0;dir<AMREX_SPACEDIM+1;dir++) {
  vel_max_estdt[dir]=0.0;
 }

 Vector< Vector<Real> > local_cap_wave_speed;
 Vector< Vector<Real> > local_vel_max_estdt;
 Vector< Real > local_dt_min;
 local_cap_wave_speed.resize(thread_class::nthreads);
 local_vel_max_estdt.resize(thread_class::nthreads);

 local_dt_min.resize(thread_class::nthreads);

 for (int tid=0;tid<thread_class::nthreads;tid++) {

  local_cap_wave_speed[tid].resize(num_interfaces); 
  for (int iten=0;iten<num_interfaces;iten++) {
   local_cap_wave_speed[tid][iten]=cap_wave_speed[iten];
  }

  local_vel_max_estdt[tid].resize(AMREX_SPACEDIM+1);//last component max|c|^2
  for (int dir=0;dir<AMREX_SPACEDIM+1;dir++) {
   local_vel_max_estdt[tid][dir]=0.0;
  }

  local_dt_min[tid]=1.0E+30;
  if (local_dt_min[tid]<dt_max) 
   local_dt_min[tid]=dt_max;
 }  // tid 

 MultiFab* velcell=getState(1,0,AMREX_SPACEDIM,cur_time_slab);
  
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

  if (pattern_test(local_caller_string,"computeInitialDt")==1) {
   // do nothing
  } else if (pattern_test(local_caller_string,"sum_integrated_quantities")==1) {
   // do nothing
  } else if 
   ((pattern_test(local_caller_string,"computeNewDt")==1)|| 
    (pattern_test(local_caller_string,"do_the_advance")==1)) {

   if ((cur_time_slab>prev_time_slab)&&
       (upper_slab_time>lower_slab_time)&&
       (upper_slab_time>0.0)&&
       (lower_slab_time>=0.0)&&
       (upper_slab_time>=cur_time_slab)&&
       (lower_slab_time<=prev_time_slab)) {

    //do nothing

   } else
    amrex::Error("slab times are incorrect");

  } else
   amrex::Error("local_caller_string invalid");

   //ngrow=0
   //Umac_Type
  MultiFab* velmac=getStateMAC(0,dir,cur_time_slab);

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(denmf->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(*denmf); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();
   int bfact=parent->Space_blockingFactor(level);

   const Real* xlo = grid_loc[gridno].lo();
   FArrayBox& Umac=(*velmac)[mfi];
   FArrayBox& Ucell=(*velcell)[mfi];
   FArrayBox& distfab=(*distmf)[mfi];
   FArrayBox& voffab=(*vofmf)[mfi];
   FArrayBox& denfab=(*denmf)[mfi];

   int local_enable_spectral=enable_spectral;

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   Real local_dt_min_thread=local_dt_min[tid_current];

   FArrayBox& solidfab=(*localMF[FSI_GHOST_MAC_MF+dir])[mfi];

   // declared in: GODUNOV_3D.F90
   fort_estdt(
    interface_mass_transfer_model.dataPtr(),
    &tid_current,
    &local_enable_spectral,
    parent->AMR_min_phase_change_rate.dataPtr(),
    parent->AMR_max_phase_change_rate.dataPtr(),
    elastic_time.dataPtr(),
    microlayer_substrate.dataPtr(),
    microlayer_angle.dataPtr(),
    microlayer_size.dataPtr(),
    macrolayer_size.dataPtr(),
    reaction_rate.dataPtr(),
    freezing_model.dataPtr(),
    Tanasawa_or_Schrage_or_Kassemi.dataPtr(),
    distribute_from_target.dataPtr(),
    saturation_temp.dataPtr(),
    mass_fraction_id.dataPtr(),
    molar_mass.dataPtr(),
    species_molar_mass.dataPtr(),
    denconst_interface_added.dataPtr(),
    Umac.dataPtr(),ARLIM(Umac.loVect()),ARLIM(Umac.hiVect()),
    Ucell.dataPtr(),ARLIM(Ucell.loVect()),ARLIM(Ucell.hiVect()),
    solidfab.dataPtr(),ARLIM(solidfab.loVect()),ARLIM(solidfab.hiVect()),
    denfab.dataPtr(),ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
    voffab.dataPtr(),ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),
    distfab.dataPtr(),ARLIM(distfab.loVect()),ARLIM(distfab.hiVect()),
    xlo,dx,
    tilelo,tilehi,
    fablo,fabhi,
    &bfact,
    &min_stefan_velocity_for_dt,
    local_cap_wave_speed[tid_current].dataPtr(),
    local_vel_max_estdt[tid_current].dataPtr(),
    &local_dt_min_thread,
    &NS_geometry_coord,
    denconst.dataPtr(),
    &visc_coef,
    &gravity_reference_wavelen,
    &dir,
    &nparts,
    &nparts_def,
    im_solid_map_ptr,
    material_type.dataPtr(),
    &cur_time_slab,
    shock_timestep.dataPtr(),
    &cfl,
    &EILE_flag, 
    &level,
    &finest_level);

   local_dt_min[tid_current]=local_dt_min_thread;

  }  //mfi
} // omp
  ns_reconcile_d_num(LOOP_ESTDT,"MaxAdvectSpeed");

  for (int tid=1;tid<thread_class::nthreads;tid++) {

   for (int iten=0;iten<num_interfaces;iten++)
    if (local_cap_wave_speed[tid][iten]>local_cap_wave_speed[0][iten])
     local_cap_wave_speed[0][iten]=local_cap_wave_speed[tid][iten];

   if (local_dt_min[tid]<local_dt_min[0])
    local_dt_min[0]=local_dt_min[tid];

   if (local_vel_max_estdt[tid][dir]>local_vel_max_estdt[0][dir])
    local_vel_max_estdt[0][dir]=local_vel_max_estdt[tid][dir];
   if (local_vel_max_estdt[tid][AMREX_SPACEDIM]>
       local_vel_max_estdt[0][AMREX_SPACEDIM])
    local_vel_max_estdt[0][AMREX_SPACEDIM]=
	    local_vel_max_estdt[tid][AMREX_SPACEDIM];

  } // tid=1..nthreads

  for (int iten=0;iten<num_interfaces;iten++) {
   ParallelDescriptor::ReduceRealMax(local_cap_wave_speed[0][iten]);
  }

  ParallelDescriptor::ReduceRealMax(local_vel_max_estdt[0][dir]);
  ParallelDescriptor::ReduceRealMax(local_vel_max_estdt[0][AMREX_SPACEDIM]);

  ParallelDescriptor::ReduceRealMin(local_dt_min[0]);

  delete velmac;
 }  // dir=0..sdim-1

 delete velcell;

 for (int iten=0;iten<num_interfaces;iten++) {
  cap_wave_speed[iten]=local_cap_wave_speed[0][iten];
 }

 for (int dir=0;dir<=AMREX_SPACEDIM;dir++) {
  vel_max_estdt[dir]=local_vel_max_estdt[0][dir];
 }

 dt_min=local_dt_min[0];

 delete denmf;
 delete vofmf;
 delete distmf;

} // subroutine MaxAdvectSpeed

// called from: do_the_advance 
//              computeNewDt (nsteps>0) 
//              computeInitialDt 
// fixed_dt==0.0 if dt not prescribed.
// fixed_dt>0.0 if dt prescribed.
//
Real NavierStokes::estTimeStep (Real local_fixed_dt,
  const std::string& caller_string) {

 std::string local_caller_string="estTimeStep";
 local_caller_string=caller_string+local_caller_string;

 Real return_dt=0.0;

 if (level!=0)
  amrex::Error("estTimeStep only called at level=0");

 const int finest_level = parent->finestLevel();
 NavierStokes& ns_fine = getLevel(finest_level);
 const Real* dxfine = ns_fine.geom.CellSize();
 Real smallest_dx=dxfine[0];
 for (int dir=1;dir<AMREX_SPACEDIM;dir++) {
  if (smallest_dx>dxfine[dir])
   smallest_dx=dxfine[dir];
 }

 Real dt_min=0.0;
  // last component is max|c|^2
 Vector<Real> u_max_estdt(AMREX_SPACEDIM+1);  

 if (local_fixed_dt>0.0) {

  return_dt=local_fixed_dt;

 } else if (local_fixed_dt==0.0) {

  if (fixed_dt_velocity > 0.0) {

   return_dt=smallest_dx/fixed_dt_velocity;

  } else if (fixed_dt_velocity==0.0) {

   MaxAdvectSpeedALL(dt_min,u_max_estdt,local_caller_string);

   if (verbose>0) {
    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "after MaxAdvectSpeedALL " << '\n';
     std::cout << "local_caller_string= " << local_caller_string << '\n';
     std::cout << "dt_min "<<dt_min<< '\n';
     for (int dir=0;dir<AMREX_SPACEDIM+1;dir++)
      std::cout << "dir u_max_estdt "<<dir<<' '<<u_max_estdt[dir]<< '\n';
    }
   }

   if (min_velocity_for_dt>0.0) {
    Real local_dt_max=smallest_dx/min_velocity_for_dt;
    if (dt_min>local_dt_max)
     dt_min=local_dt_max;
   }

   return_dt=cfl*dt_min;

  } else {

   amrex::Error("fixed_dt_velocity invalid");

  }

  Vector<Real> time_array;
  time_array.resize(ns_time_order+1);
  Real slablow=0.0;
  Real slabhigh=1.0;
  int slab_dt_type=parent->get_slab_dt_type();
  fort_gl_slab(time_array.dataPtr(),
               &slab_dt_type,
               &ns_time_order,
               &slablow,&slabhigh);
  Real max_sub_dt=0.0;
  for (int i=0;i<ns_time_order;i++) {
   Real test_sub_dt=time_array[i+1]-time_array[i];
   if (test_sub_dt<=0.0)
    amrex::Error("test_sub_dt invalid");
   if (test_sub_dt>max_sub_dt)
    max_sub_dt=test_sub_dt;
  } // i
  if ((max_sub_dt<=0.0)||(max_sub_dt>1.0))
   amrex::Error("max_sub_dt invalid");

  if (verbose>0)
   if (ParallelDescriptor::IOProcessor())
    std::cout << "ns_time_order= " << ns_time_order << 
     " time slab factor= " << 1.0/max_sub_dt << '\n';

  return_dt=return_dt/max_sub_dt;

 } else {
  amrex::Error("local_fixed_dt invalid");
 }

 if (return_dt>dt_max)
  return_dt=dt_max;

 return return_dt;

} // subroutine estTimeStep

// post_regrid is called from either:
// 1. AmrCore::initialInit, 
// 2. AmrCore::regrid_level_0_on_restart, or
// 3. AmrCore::regrid 
void NavierStokes::post_regrid (int lbase,
  int start_level,int new_finest,int initialInit_flag,Real time) {

  const int max_level = parent->maxLevel();

  if ((lbase>=0)&&(lbase<=max_level)) {
   //do nothing
  } else
   amrex::Error("lbase invalid");

  if ((start_level>=0)&&(start_level<=max_level)) {
   //do nothing
  } else
   amrex::Error("start_level invalid");

  Real dt_amr=parent->getDt(); // returns dt_AMR
  int nstate=state.size();
  if (nstate!=NUM_STATE_TYPE)
   amrex::Error("nstate invalid");

  if (initialInit_flag==1) {
   // do nothing
  } else if (initialInit_flag==0) {
   // do nothing
  } else
   amrex::Error("initialInit_flag invalid");

  NavierStokes& ns_level0=getLevel(0);

#ifdef AMREX_PARTICLES

  int lev_min=0;
  int lev_max=new_finest;
  int nGrow_Redistribute=0;
  int local_Redistribute=0;

  My_ParticleContainer& current_PC=ns_level0.newDataPC(ns_time_order);
  current_PC.Redistribute(lev_min,lev_max,nGrow_Redistribute, 
     local_Redistribute);

  Long num_particles=current_PC.TotalNumberOfParticles();

  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "TotalNumberOfParticles for slab ns_time_order= " <<
     ns_time_order << " is equal to " << num_particles << '\n';
  }

#endif

   // olddata=newdata  
  for (int k=0;k<nstate;k++) {
   state[k].CopyNewToOld(level,max_level); 
   state[k].setTimeLevel(time,dt_amr);
  }
#ifdef AMREX_PARTICLES
  ns_level0.CopyNewToOldPC(lev_max);
#endif

} // end subroutine post_regrid

void NavierStokes::computeNewDt (int finest_level,
  Real& dt,Real stop_time) {

 std::string local_caller_string="computeNewDt";

 int nsteps=parent->levelSteps(0);

 if (nsteps>0) {

  Real local_fixed_dt;
  Real local_change_max;
  if (nsteps==1) {
   local_fixed_dt=fixed_dt;
   local_change_max=change_max_init;
  } else if (nsteps>1) {
   local_fixed_dt=fixed_dt;
   local_change_max=change_max;
  } else
   amrex::Error("nsteps invalid");
   
  if (verbose>0) {
   if (ParallelDescriptor::IOProcessor()) {
    std::cout << "start: computeNewDt nsteps=" << nsteps << 
     " local_fixed_dt= " << local_fixed_dt << " local_change_max= " <<
     local_change_max << '\n';
   }
  }

  int max_level = parent->maxLevel();

  if ((finest_level>=0)&&(finest_level<=max_level)) {
   // do nothing
  } else
   amrex::Error("finest_level invalid");


  if (level==0) {

   Real newdt=estTimeStep(local_fixed_dt,local_caller_string);

   if ((local_fixed_dt==0.0)&&(fixed_dt_velocity==0.0)) {
    if  (newdt>local_change_max*dt)
     newdt=local_change_max*dt;
   } else if ((local_fixed_dt>0.0)||(fixed_dt_velocity>0.0)) {
    // do nothing
   } else {
    amrex::Error("local_fixed_dt or fixed_dt_velocity invalid");
   }

   Real dt_0=newdt;

   const Real eps      = 0.0001*dt_0;
   const Real eps2     = 0.000001*dt_0;
   upper_slab_time = state[State_Type].slabTime(ns_time_order);
   if (stop_time >= 0.0) {
      if ((upper_slab_time + dt_0) > (stop_time - eps))
          dt_0 = stop_time - upper_slab_time + eps2;
   }

   const Real check_per = parent->checkPer();
   if (check_per > 0.0) {
      int a = int((upper_slab_time + eps ) / check_per);
      int b = int((upper_slab_time + dt_0) / check_per);
      if (a != b)
          dt_0 = b * check_per - upper_slab_time;
   }

   const Real plot_per = parent->plotPer();
   if (plot_per > 0.0) {
      int a = int((upper_slab_time + eps ) / plot_per);
      int b = int((upper_slab_time + dt_0) / plot_per);
      if (a != b)
          dt_0 = b * plot_per - upper_slab_time;
   }

   dt=dt_0;
  } else if ((level>0)&&(level<=max_level)) {
   // do nothing
  } else
   amrex::Error("level invalid computeNewDt");

  if (verbose>0) {
   if (ParallelDescriptor::IOProcessor()) {
    std::cout << "end: computeNewDt \n";
   }
  }

 } else
  amrex::Error("nsteps invalid");

} //end subroutine computeNewDt

void NavierStokes::computeInitialDt (int finest_level,
   Real& dt,Real stop_time) {

 std::string local_caller_string="computeInitialDt";

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "start: computeInitialDt \n";
  }
 }

 int max_level = parent->maxLevel();

 if ((finest_level>=0)&&(finest_level<=max_level)) {
   // do nothing
 } else
   amrex::Error("finest_level invalid");

 if (level!=0)
  amrex::Error("level invalid computeInitialDt");

  // The first time step does not have grad p initialized, so we
  // must take a very small initial time step in order to maintain
  // spectral accuracy.
 Real shrink_factor=1.0;
 if ((ns_time_order>=2)&&(ns_time_order<=32)) {
  for (int i=0;i<ns_time_order;i++)
   shrink_factor*=2.0;
 } else if (ns_time_order==1) {
  // do nothing
 } else {
  amrex::Error("ns_time_order invalid");
 }

 Real newdt=
   init_shrink*estTimeStep(fixed_dt_init,local_caller_string)/shrink_factor;

 Real dt_0 = newdt;

 const Real eps      = 0.0001*dt_0;
 upper_slab_time = state[State_Type].slabTime(ns_time_order);
 if (stop_time >= 0.0) {
     if ((upper_slab_time + dt_0) > (stop_time - eps))
         dt_0 = stop_time - upper_slab_time;
 }

 dt=dt_0;

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "end: computeInitialDt \n";
  }
 }

} // subroutine computeInitialDt

//
// Fills in amrLevel okToContinue.
//

int
NavierStokes::okToContinue ()
{
    return true;
}


void
NavierStokes::post_timestep (Real stop_time) {

 std::string local_caller_string="post_timestep";

 int max_level = parent->maxLevel();

 SDC_outer_sweeps=0;
 slab_step=ns_time_order-1;
 SDC_setup_step();

 if (level==0) {
  if ((sum_interval>0)||
      (visual_drag_plot_int>0)) {

   int sum_interval_trigger=0;
   if (sum_interval>0) {
    if (parent->levelSteps(0)%sum_interval == 0) {
     sum_interval_trigger=1;
    }
   } else if ((sum_interval==0)||(sum_interval==-1)) {
    // do nothing
   } else
    amrex::Error("sum_interval invalid");

   int visual_drag_plot_int_trigger=0;

   if (visual_drag_plot_int>0) {
    if (parent->levelSteps(0)%visual_drag_plot_int == 0) {
     visual_drag_plot_int_trigger=1;
    }
   } else if ((visual_drag_plot_int==0)||(visual_drag_plot_int==-1)) {
    // do nothing
   } else
    amrex::Error("visual_drag_plot_int invalid");

   if ( (sum_interval_trigger==1)||
        (visual_drag_plot_int_trigger==1)||
        (stop_time-upper_slab_time<1.0E-8) ) {
    sum_integrated_quantities(local_caller_string,stop_time);
   }
  } else if (((sum_interval==0)||
	      (sum_interval==-1))&&
	     ((visual_drag_plot_int==0)||
	      (visual_drag_plot_int==-1))) {
   // do nothing
  } else {
   amrex::Error("sum_interval or visual_drag_plot_int invalid");
  }
 } else if ((level<=max_level)&&(level>0)) {
  // do nothing
 } else {
  amrex::Error("level invalid");
 } 

 init_regrid_history();

} // end subroutine post_timestep

//
// Ensure state, and pressure are consistent.
//
void
NavierStokes::post_init (Real stop_time)
{

 std::string local_caller_string="post_init";

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "start: post_init \n";
  }
 }

 SDC_setup();    
 ns_time_order=parent->Time_blockingFactor();
 slab_step=ns_time_order-1;

 const int finest_level = parent->finestLevel();

 if (level==0) {

   // in post_init: delete FSI_MF used by initData
  for (int ilev = 0; ilev <= finest_level; ilev++) {
   NavierStokes& ns_level = getLevel(ilev);
   ns_level.delete_localMF_if_exist(FSI_MF,1);
  }

  Real strt_time = state[State_Type].slabTime(ns_time_order);

   // bogus value for dt to use while doing initial project - 
   // differentiate between prev_time and cur_time.
   // newtime=strt_time, oldtime=strt_time-dt_save

  Real dt_save=1.0;

  for (int ilev = 0; ilev <= finest_level; ilev++) {
   NavierStokes& ns_level = getLevel(ilev);
   ns_level.setTimeLevel(strt_time,dt_save);
  }

  parent->setDt(dt_save);

  // Ensure state is consistent, i.e. velocity field is non-divergent,
  // Coarse levels are fine level averages, 
  //
  post_init_state();

  //
  // Re-Estimate the initial timestepping. 

  computeInitialDt(finest_level,dt_save,stop_time);

   // newtime=strt_time, oldtime=strt_time-dt_save
  for (int ilev = 0; ilev <= finest_level; ilev++) {
   NavierStokes& ns_level = getLevel(ilev);
   ns_level.setTimeLevel(strt_time,dt_save);
  }

  parent->setDt(dt_save);

  //
  // Compute the initial estimate of conservation.
  //
  int sum_interval_local=sum_interval;

  if (sum_interval_local > 0) {
     sum_integrated_quantities(local_caller_string,stop_time);
  }

  ParallelDescriptor::Barrier();
 } else if ((level>=1)&&(level<=finest_level)) {
  // do nothing
 } else {
  amrex::Error("level invalid");
 } 

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "end: post_init \n";
  }
 }

}  // end subroutine post_init


// status==1 success
// status==0 failure
void matrix_solveCPP(Real** AA,Real* xx,Real* bb,
 int& status,int numelem) {
        
 Real alpha,holdvalue;
 int i,j,k,holdj;

 status=1;
 for (i=1;i<=numelem-1;i++) {
  holdj=i;
  holdvalue=std::abs(AA[i-1][i-1]);
  for (j=i+1;j<=numelem;j++) {
   if (std::abs(AA[j-1][i-1])>holdvalue) {
    holdj=j;
    holdvalue=std::abs(AA[j-1][i-1]);
   }
  }
  if (holdj!=i) {
   for (j=i;j<=numelem;j++) {
    holdvalue=AA[i-1][j-1];
    AA[i-1][j-1]=AA[holdj-1][j-1];
    AA[holdj-1][j-1]=holdvalue;
   }
  }
  holdvalue=bb[i-1];
  bb[i-1]=bb[holdj-1];
  bb[holdj-1]=holdvalue;
  if (std::abs(AA[i-1][i-1])<1.0E-32) 
   status=0;
  else {
   for (j=i+1;j<=numelem;j++) {
    alpha=AA[j-1][i-1]/AA[i-1][i-1];
    for (k=i;k<=numelem;k++)
     AA[j-1][k-1]=AA[j-1][k-1]-alpha*AA[i-1][k-1];
    bb[j-1]=bb[j-1]-alpha*bb[i-1];
   }
  }
 }

 for (i=numelem;i>=1;i--) {
  if (status!=0) {
   holdvalue=bb[i-1];
   for (j=i+1;j<=numelem;j++)
    holdvalue=holdvalue-AA[i-1][j-1]*xx[j-1];
   if (std::abs(AA[i-1][i-1])<1.0E-32) 
    status=0;
   else
    xx[i-1]=holdvalue/AA[i-1][i-1];
  }
 }

} // end subroutine matrix_solveCPP

//called from: NavierStokes::sum_integrated_quantities (NS_setup.cpp)
//   sum_integrated_quantities called from post_timestep or
//   if sum_integrated_quantities called from post_init or
//   if sum_integrated_quantities called from post_restart or
//   if called from init_FSI_GHOST_MAC_MF_ALL or
//   if called from nonlinear_advection
void
NavierStokes::volWgtSumALL(
  const std::string& caller_string,
  int fast_mode) {

 std::string local_caller_string="volWgtSumALL";
 local_caller_string=caller_string+local_caller_string;

 int finest_level=parent->finestLevel();
 NavierStokes& ns_fine = getLevel(finest_level);

 if (level!=0)
  amrex::Error("it is required that level=0 in volWgtSumALL");

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else {
  std::cout << "SDC_outer_sweeps= " << SDC_outer_sweeps << '\n';
  amrex::Error("SDC_outer_sweeps invalid");
 }

 if (fast_mode==1) {
  init_FSI_GHOST_MAC_MF_ALL_predict();
 } else if (fast_mode==0) {
  // do nothing
 } else
  amrex::Error("fast_mode invalid");

    // vof,ref cen, order,slope,int
 int update_flag=RECON_UPDATE_NULL;
 int init_vof_prev_time=0;
  //output: SLOPE_RECON_MF
 VOF_Recon_ALL(1,cur_time_slab,update_flag,
   init_vof_prev_time); 

  // need to initialize viscosity and density temporary 
  // variables.
  // in: volWgtSumALL
  //
 if (pattern_test(local_caller_string,"post_timestep")==1) {
  // do nothing
 } else if (pattern_test(local_caller_string,"post_init")==1) {
  // do nothing
 } else if (pattern_test(local_caller_string,"post_restart")==1) {
  // do nothing
 } else if (pattern_test(local_caller_string,"init_FSI_GHOST_MAC_MF_ALL")==1) {

  if (fast_mode==1) {
   // do nothing
  } else
   amrex::Error("expecting fast_mode==1");

 } else if (pattern_test(local_caller_string,"nonlinear_advection")==1) {
  // do nothing
 } else
  amrex::Error("local_caller_string invalid 21614");

 allocate_levelset_ALL(2,LEVELPC_MF);

 if (fast_mode==0) {
  //make_physics_varsALL calls "init_gradu_tensor_and_material_visc_ALL"
  make_physics_varsALL(SOLVETYPE_INITPROJ,local_caller_string);
 } else if (fast_mode==1) {
  //localMF[CELL_VISC_MATERIAL_MF] is deleted in ::Geometry_cleanup()
  //responsibility of caller to issue commands,
  // delete_array(CELLTENSOR_MF);
  // delete_array(FACETENSOR_MF);
  //
  init_gradu_tensor_and_material_visc_ALL(local_caller_string);
 } else
  amrex::Error("fast_mode invalid");

  // see <DRAG_COMP.H>
 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.delete_localMF_if_exist(DRAG_MF,1);
 }

  // initializes DRAG_MF to 0.0.
  // ngrow_make_distance=3
  // ngrow_distance=4
 allocate_array(ngrow_make_distance,N_DRAG,-1,DRAG_MF);

 debug_ngrow(CELL_VISC_MATERIAL_MF,1,local_caller_string);
 if (localMF[CELL_VISC_MATERIAL_MF]->nComp()==3*num_materials) {
  // do nothing
 } else {
  amrex::Error("volWgtSumALL: CELL_VISC_MATERIAL_MF invalid ncomp");
 }

 if (ENUM_NUM_TENSOR_TYPE==2*AMREX_SPACEDIM) {
  // do nothing
 } else
  amrex::Error("expecting ENUM_NUM_TENSOR_TYPE==2*AMREX_SPACEDIM");

 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=num_materials)) {

  //ngrow,ncomp,grid_type,mf id
  //VISCOTEN_ALL_MAT_MF initialized to 0.0
  allocate_array(1,num_materials_viscoelastic*ENUM_NUM_TENSOR_TYPE,-1,
	 VISCOTEN_ALL_MAT_MF);

  for (int im=0;im<num_materials;im++) {

   if (ns_is_rigid(im)==0) {

    if (store_elastic_data[im]==1) {
     if (elastic_viscosity[im]>0.0) {

      int partid=0;
      while ((im_elastic_map[partid]!=im)&&(partid<im_elastic_map.size())) {
       partid++;
      }

      if (partid<im_elastic_map.size()) {
       // we are currently in "volWgtSumALL"
       make_viscoelastic_tensorALL(im); // (mu_p/lambda)(f(A)A-I) if FENE-P
       Copy_array(VISCOTEN_ALL_MAT_MF,VISCOTEN_MF,
         0,partid*ENUM_NUM_TENSOR_TYPE,ENUM_NUM_TENSOR_TYPE,1);
       delete_array(VISCOTEN_MF);
      } else
       amrex::Error("partid could not be found: volWgtSumALL");
     } else
      amrex::Error("elastic_viscosity invalid");
    } else if (store_elastic_data[im]==0) {

     if (viscoelastic_model[im]!=0)
      amrex::Error("viscoelastic_model[im]!=0");

    } else
     amrex::Error("elastic_time/elastic_viscosity invalid");

   } else if (ns_is_rigid(im)==1) {
    // do nothing
   } else
    amrex::Error("ns_is_rigid invalid");

  } // im=0..num_materials-1
 } else if (num_materials_viscoelastic==0) {
  // do nothing
 } else
  amrex::Error("num_materials_viscoelastic invalid:volWgtSumALL");

 if (fast_mode==0) {
  GetDragALL();
 } else if (fast_mode==1) {
  // do nothing (need, e.g., maximum temperature in the wall which can be
  // computed without drag information)
 } else
  amrex::Error("fast_mode invalid");

 int num_cells=0;
 int Z_dir=1;
 int R_dir=0;
 const Box& fdomain = ns_fine.geom.Domain();
 const int* fdomlo = fdomain.loVect();
 const int* fdomhi = fdomain.hiVect();
 NS_coflow_Z.resize(num_cells+1);
 NS_coflow_R_of_Z.resize(num_cells+1);

 if (fast_mode==1) {
  // do nothing
 } else if (fast_mode==0) {
  fort_coflow(
   &upper_slab_time, 
   fdomlo,
   fdomhi,
   &Z_dir,
   &R_dir,
   &num_cells,
   NS_coflow_Z.dataPtr(),
   NS_coflow_R_of_Z.dataPtr());

  NS_coflow_Z.resize(num_cells+1);
  NS_coflow_R_of_Z.resize(num_cells+1);
 } else
  amrex::Error("fast_mode invalid");

 for (int jfine=0;jfine<=num_cells;jfine++) {
  NS_coflow_Z[jfine]=0.0;
  NS_coflow_R_of_Z[jfine]=0.0;
 }

 for (int isweep=0;isweep<2;isweep++) {

  for (int ilev = 0; ilev <= finest_level; ilev++) {

   NavierStokes& ns_level = getLevel(ilev);
   ns_level.volWgtSum(isweep,fast_mode);

  }  // ilev=0..finest_level 

  if (isweep==0) {

   for (int im=0;im<num_materials;im++) {
    Real volmat=NS_sumdata[IQ_FE_SUM_COMP+2*im];
    Real LSvolmat=NS_sumdata[IQ_LS_F_SUM_COMP+im];
    if (volmat>0.0) {
     for (int dir=0;dir<AMREX_SPACEDIM;dir++)
      NS_sumdata[3*im+IQ_CEN_SUM_COMP+dir]=
	      NS_sumdata[3*im+IQ_CEN_SUM_COMP+dir]/volmat;
    }
    if (LSvolmat>0.0) {
     for (int dir=0;dir<AMREX_SPACEDIM;dir++)
      NS_sumdata[3*im+IQ_LS_CEN_SUM_COMP+dir]=
       NS_sumdata[3*im+IQ_LS_CEN_SUM_COMP+dir]/LSvolmat;
    }
   } // im=0..num_materials-1
  }  // isweep=0
 } //for (int isweep=0;isweep<2;isweep++) 

 int local_comp=0;
 for (int im=0;im<num_materials;im++) {
  for (int dir=0;dir<3;dir++) {
   int idest=IQ_BODYDRAG_SUM_COMP+local_comp;
   int isource=DRAGCOMP_IQ_BODYFORCE+local_comp;
   NS_sumdata[idest]=NS_DRAG_integrated_quantities[isource];
   idest=IQ_BODYTORQUE_SUM_COMP+local_comp;
   isource=DRAGCOMP_IQ_BODYTORQUE+local_comp;
   NS_sumdata[idest]=NS_DRAG_integrated_quantities[isource];
   idest=IQ_DRAG_SUM_COMP+local_comp;
   isource=DRAGCOMP_IQ_FORCE+local_comp;
   NS_sumdata[idest]=NS_DRAG_integrated_quantities[isource];
   idest=IQ_PDRAG_SUM_COMP+local_comp;
   isource=DRAGCOMP_IQ_PFORCE+local_comp;
   NS_sumdata[idest]=NS_DRAG_integrated_quantities[isource];
   idest=IQ_VISCOUSDRAG_SUM_COMP+local_comp;
   isource=DRAGCOMP_IQ_VISCOUSFORCE+local_comp;
   NS_sumdata[idest]=NS_DRAG_integrated_quantities[isource];
   idest=IQ_VISCOUS0DRAG_SUM_COMP+local_comp;
   isource=DRAGCOMP_IQ_VISCOUS0FORCE+local_comp;
   NS_sumdata[idest]=NS_DRAG_integrated_quantities[isource];
   idest=IQ_VISCODRAG_SUM_COMP+local_comp;
   isource=DRAGCOMP_IQ_VISCOFORCE+local_comp;
   NS_sumdata[idest]=NS_DRAG_integrated_quantities[isource];

   idest=IQ_TORQUE_SUM_COMP+local_comp;
   isource=DRAGCOMP_IQ_TORQUE+local_comp;
   NS_sumdata[idest]=NS_DRAG_integrated_quantities[isource];
   idest=IQ_PTORQUE_SUM_COMP+local_comp;
   isource=DRAGCOMP_IQ_PTORQUE+local_comp;
   NS_sumdata[idest]=NS_DRAG_integrated_quantities[isource];
   idest=IQ_VISCOUSTORQUE_SUM_COMP+local_comp;
   isource=DRAGCOMP_IQ_VISCOUSTORQUE+local_comp;
   NS_sumdata[idest]=NS_DRAG_integrated_quantities[isource];
   idest=IQ_VISCOUS0TORQUE_SUM_COMP+local_comp;
   isource=DRAGCOMP_IQ_VISCOUS0TORQUE+local_comp;
   NS_sumdata[idest]=NS_DRAG_integrated_quantities[isource];
   idest=IQ_VISCOTORQUE_SUM_COMP+local_comp;
   isource=DRAGCOMP_IQ_VISCOTORQUE+local_comp;
   NS_sumdata[idest]=NS_DRAG_integrated_quantities[isource];

   local_comp++;
  } // dir=0 ... 2
 } // im=0 .. num_materials-1

 for (int im=0;im<num_materials;im++) {

  int idest=IQ_STEP_PERIM_SUM_COMP+im;
  int isource=DRAGCOMP_IQ_PERIM+im;
  NS_sumdata[idest]=NS_DRAG_integrated_quantities[isource];

 } // im=0 .. num_materials-1

 if (num_cells>0) {
  if (ParallelDescriptor::IOProcessor()) {
   fort_coflow(
    &upper_slab_time, 
    fdomlo,
    fdomhi,
    &Z_dir,
    &R_dir,
    &num_cells,
    NS_coflow_Z.dataPtr(),
    NS_coflow_R_of_Z.dataPtr());
  }
 } else if (num_cells==0) {
  // do nothing
 } else
  amrex::Error("num_cells invalid");

 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=num_materials)) {
  delete_array(VISCOTEN_ALL_MAT_MF);
 } else if (num_materials_viscoelastic==0) {
  // do nothing
 } else
  amrex::Error("num_materials_viscoelastic invalid:volWgtSumALL");

 delete_array(CELLTENSOR_MF);
 delete_array(FACETENSOR_MF);

}  // end subroutine volWgtSumALL

void
NavierStokes::MaxPressureVelocityALL(
   Real& minpres,Real& maxpres,
   Real& maxvel,Real& maxvel_collide) {

 int finest_level=parent->finestLevel();
 if (level!=0)
  amrex::Error("level=0 in MaxPressureVelocityALL");

 Real local_minpres;
 Real local_maxpres;
 Real local_maxvel;
 Real local_maxvel_collide;
 maxpres=-1.0e+99;
 minpres=1.0e+99;
 maxvel=0.0;
 maxvel_collide=0.0;
 for (int k = 0; k <= finest_level; k++) {
  NavierStokes& ns_level = getLevel(k);
  ns_level.MaxPressureVelocity(local_minpres,local_maxpres,
   local_maxvel,local_maxvel_collide);
  if (local_minpres<minpres)
   minpres=local_minpres;
  if (local_maxpres>maxpres)
   maxpres=local_maxpres;
  if (local_maxvel>maxvel)
   maxvel=local_maxvel;
  if (local_maxvel_collide>maxvel_collide)
   maxvel_collide=local_maxvel_collide;
 }

} // end subroutine MaxPressureVelocityALL

void 
NavierStokes::MaxPressureVelocity(Real& minpres,Real& maxpres,
 Real& maxvel,Real& maxvel_collide) {
 
 bool use_tiling=ns_tiling;

  // ngrow=0  
 MultiFab* vel=getState(0,STATECOMP_VEL,
   STATE_NCOMP_VEL+STATE_NCOMP_PRES,cur_time_slab);
 MultiFab* velmac[AMREX_SPACEDIM];
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   //Umac_Type
  velmac[dir]=getStateMAC(0,dir,cur_time_slab);
 }
 
  // mask=tag if not covered by level+1 
  // ngrow=0 so no need to worry about ghost cells.
 int ngrowmask=0;
 Real tag=1.0;
 int clear_phys_boundary=2;
 MultiFab* mask=maskfiner(ngrowmask,tag,clear_phys_boundary);

 const Real* dx = geom.CellSize();
 maxpres=-1.0e+99;
 minpres=1.0e+99;
 maxvel=0.0;
 maxvel_collide=0.0;

 Vector<Real> minpresA;
 Vector<Real> maxpresA;
 Vector<Real> maxvelA;
 Vector<Real> maxvel_collideA;
 minpresA.resize(thread_class::nthreads);
 maxpresA.resize(thread_class::nthreads);
 maxvelA.resize(thread_class::nthreads);
 maxvel_collideA.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  minpresA[tid]=1.0e+99;
  maxpresA[tid]=-1.0e+99;
  maxvelA[tid]=0.0;
  maxvel_collideA[tid]=0.0;
 }

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(mask->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*mask,use_tiling); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();
  int bfact=parent->Space_blockingFactor(level);

  const Real* xlo = grid_loc[gridno].lo();
  FArrayBox& maskfab=(*mask)[mfi];
  FArrayBox& velfab=(*vel)[mfi];
  FArrayBox& velx=(*velmac[0])[mfi];
  FArrayBox& vely=(*velmac[1])[mfi];
  FArrayBox& velz=(*velmac[AMREX_SPACEDIM-1])[mfi];
  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // declared in: DERIVE_3D.F90
  fort_maxpresvel(
   &minpresA[tid_current],
   &maxpresA[tid_current],
   &maxvelA[tid_current],
   &maxvel_collideA[tid_current],
   xlo,dx,
   maskfab.dataPtr(),ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
   velfab.dataPtr(),ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
   velx.dataPtr(),ARLIM(velx.loVect()),ARLIM(velx.hiVect()),
   vely.dataPtr(),ARLIM(vely.loVect()),ARLIM(vely.hiVect()),
   velz.dataPtr(),ARLIM(velz.loVect()),ARLIM(velz.hiVect()),
   tilelo,tilehi,
   fablo,fabhi,&bfact);
 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_MAXPRESVEL,"MaxPressureVelocity");

 for (int tid=0;tid<thread_class::nthreads;tid++) {
  if (minpres>minpresA[tid])
   minpres=minpresA[tid];
  if (maxpres<maxpresA[tid])
   maxpres=maxpresA[tid];
  if (maxvel<maxvelA[tid])
   maxvel=maxvelA[tid];
  if (maxvel_collide<maxvel_collideA[tid])
   maxvel_collide=maxvel_collideA[tid];
 }
 ParallelDescriptor::ReduceRealMin(minpres);
 ParallelDescriptor::ReduceRealMax(maxpres);
 ParallelDescriptor::ReduceRealMax(maxvel);
 ParallelDescriptor::ReduceRealMax(maxvel_collide);

 delete mask;
 delete vel;
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  delete velmac[dir];
 }
} // end subroutine MaxPressureVelocity

// called from 
// writePlotFile   
// post_init_state 
// post_restart   
void
NavierStokes::prepare_post_process(const std::string& caller_string) {

 if (level!=0)
  amrex::Error("level invalid prepare_post_process");

 std::string local_caller_string="prepare_post_process";
 local_caller_string=caller_string+local_caller_string;

 int max_level = parent->maxLevel();
 int finest_level = parent->finestLevel();
 if ((max_level>=0)&&(finest_level<=max_level)) {
  //do nothing
 } else
  amrex::Error("max_level invalid");
 if ((finest_level>=0)&&(finest_level<=max_level)) {
  //do nothing
 } else
  amrex::Error("finest_level invalid");

 //init VOLUME_MF and AREA_MF; metrics_dataALL is declared in NavierStokes2.cpp
 metrics_dataALL(1);

//note: fort_initgridmap is called from:
// (i) NavierStokes::initData ()
// (ii) NavierStokes::post_restart()

 if (pattern_test(local_caller_string,"post_init_state")==1) {
   // called from post_init_state
  MOF_training();
 } else if (pattern_test(local_caller_string,"post_restart")==1) {
  // called from post_restart
  MOF_training();
 } else if (pattern_test(local_caller_string,"writePlotFile")==1) {
  // called from writePlotFile
  // do nothing
 } else
  amrex::Error("local_caller_string invalid");

 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);

  ns_level.allocate_mdot(); //MDOT_MF=0.0

    // mask=tag if not covered by level+1 or outside the domain.
  Real tag=1.0;
  int clearbdry=0; 
   // ngrow=1
  ns_level.maskfiner_localMF(MASKCOEF_MF,1,tag,clearbdry);
  ns_level.prepare_mask_nbr(1);

  if (pattern_test(local_caller_string,"post_init_state")==1) {
   // called from post_init_state
   // do nothing
  } else if (pattern_test(local_caller_string,"post_restart")==1) {
   // called from post_restart
   // do nothing
  } else if (pattern_test(local_caller_string,"writePlotFile")==1) {
    // called from writePlotFile
    // in: NavierStokes::prepare_post_process
   ns_level.allocate_levelset(1,LEVELPC_MF);
  } else
   amrex::Error("local_caller_string invalid 22064");
   
 } // ilev=level ... finest_level

 build_masksemALL();

 int filler_renormalize_only=1;
 init_FSI_GHOST_MAC_MF_ALL(filler_renormalize_only,local_caller_string);

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  if (ilev<finest_level) {
   ns_level.MOFavgDown();
  }
 } // ilev=finest_level ... level

 int init_vof_prev_time=0;
 int error_update_flag=0;
 int renormalize_only=0; // init:solid TEMP,VEL,LS,extend LSfluid into solid.
 int local_truncate=0; // do not force removal of flotsam.

 if (pattern_test(local_caller_string,"writePlotFile")==1) {
  // called from writePlotFile, do not update S_new
  error_update_flag=RECON_UPDATE_NULL;  
 } else if (pattern_test(local_caller_string,"post_init_state")==1) {
  // called from post_init_state, update S_new
  error_update_flag=RECON_UPDATE_STATE_ERR_AND_CENTROID;  
 } else if (pattern_test(local_caller_string,"post_restart")==1) {
  // called from post_restart, update S_new
  error_update_flag=RECON_UPDATE_STATE_ERR_AND_CENTROID;  
 } else
  amrex::Error("local_caller_string invalid 22091");
	
  //output:SLOPE_RECON_MF
 VOF_Recon_ALL(1,cur_time_slab,error_update_flag,
  init_vof_prev_time);

 if (pattern_test(local_caller_string,"post_init_state")==1) {

  int keep_all_interfaces=1;
  int ngrow_make_distance_accept=ngrow_make_distance;
  makeStateDistALL(keep_all_interfaces,ngrow_make_distance_accept);
  prescribe_solid_geometryALL(cur_time_slab,
		  renormalize_only,
		  local_truncate,
		  local_caller_string);

 } else if (pattern_test(local_caller_string,"writePlotFile")==1) {
  // called from writePlotFile

 } else if (pattern_test(local_caller_string,"post_restart")==1) {
  // called from post_restart

  if (1==0) {
   int keep_all_interfaces=0;
   int ngrow_make_distance_accept=ngrow_make_distance;
   makeStateDistALL(keep_all_interfaces,ngrow_make_distance_accept);
   prescribe_solid_geometryALL(cur_time_slab,renormalize_only,
		   local_truncate,
		   local_caller_string);
  }

 } else
  amrex::Error("local_caller_string invalid 22125");

 make_physics_varsALL(SOLVETYPE_INITPROJ,local_caller_string);
 delete_array(CELLTENSOR_MF);
 delete_array(FACETENSOR_MF);

}  // end subroutine prepare_post_process

#ifdef AMREX_PARTICLES

// This routine called from:
// 1. post_init_state() and
// 2. nonlinear_advection()
// 3. phase_change_code_segment()
// 4. nucleation_code_segment()
// 5. do_the_advance()
void
NavierStokes::init_particle_containerALL(int append_flag,
	const std::string& caller_string) {

 int max_level = parent->maxLevel();
 int finest_level=parent->finestLevel();

 if (finest_level<=max_level) {
  // do nothing
 } else
  amrex::Error("max_level invalid");

 if ((slab_step>=0)&&(slab_step<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("expecting 0<=slab_step<ns_time_order");

 if (level==0) {
  // do nothing
 } else 
  amrex::Error("expecting level==0 init_particle_containerALL");

 std::string local_caller_string="init_particle_containerALL";
 local_caller_string=caller_string+local_caller_string;

 int num_neighbors=1;
 if (append_flag==OP_PARTICLE_INIT) {
  num_neighbors=0;

  if (slab_step==ns_time_order-1) {
   // do nothing
  } else
   amrex::Error("expecting slab_step==ns_time_order-1");

 } else if (append_flag==OP_PARTICLE_ADD) {
  num_neighbors=0;
 } else if (append_flag==OP_PARTICLE_ASSIMILATE) {
  num_neighbors=1;
 } else
  amrex::Error("append_flag invalid");

 My_ParticleContainer& localPC=newDataPC(slab_step+1);
 Long num_particles=localPC.TotalNumberOfParticles();
 if (ParallelDescriptor::IOProcessor()) {
  std::cout << local_caller_string << '\n';
  std::cout << "append_flag= " << append_flag << '\n';
  std::cout << "TotalNumberOfParticles for slab_step+1 = " <<
     slab_step+1 << " is equal to " << num_particles << '\n';
 }

  //m_num_neighbor_cells=nneighbor=ncells=num_neighbors
 NBR_Particle_Container=
  new My_NBR_ParticleContainer(parent->GetParGDB(),num_neighbors);

 int lev_min=0;
 int lev_max=finest_level;
 int nGrow_Redistribute=0;
 bool local_copy=true; //do not redistribute inside of copyParticles
 int local_redistribute=0;

 NBR_Particle_Container->clearParticles();
 NBR_Particle_Container->Redistribute();
 NBR_Particle_Container->copyParticles(localPC,local_copy);
 NBR_Particle_Container->Redistribute(lev_min,lev_max,nGrow_Redistribute, 
    local_redistribute);
 NBR_Particle_Container->fillNeighbors();

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.init_particle_container(append_flag,local_caller_string);
 }

 NBR_Particle_Container->clearNeighbors();
 delete NBR_Particle_Container;

 if ((append_flag==OP_PARTICLE_INIT)||
     (append_flag==OP_PARTICLE_ADD)) {
  localPC.Redistribute(lev_min,lev_max,nGrow_Redistribute, 
    local_redistribute);
 } else if (append_flag==OP_PARTICLE_ASSIMILATE) {
  // do nothing
 } else {
  amrex::Error("append_flag invalid");
 }

} // end subroutine init_particle_containerALL
 
void
NavierStokes::init_particle_container(int append_flag,
   const std::string& caller_string) {

 std::string local_caller_string="init_particle_container";
 local_caller_string=caller_string+local_caller_string;

 bool use_tiling=ns_tiling;
 int max_level = parent->maxLevel();
 int finest_level=parent->finestLevel();

 if (finest_level<=max_level) {
  // do nothing
 } else
  amrex::Error("max_level invalid");

 if ((slab_step>=0)&&(slab_step<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("expecting 0<=slab_step<ns_time_order");

 if ((level>=0)&&(level<=finest_level)) {
  // do nothing
 } else 
  amrex::Error("0<=level<=finest_level failed");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 const Real* dx = geom.CellSize();

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,1,local_caller_string);

 MultiFab* lsmf=getStateDist(1,cur_time_slab,local_caller_string); 

 if (lsmf->nComp()==num_materials*(AMREX_SPACEDIM+1)) {
  //do nothing
 } else
  amrex::Error("lsmf->nComp() invalid");

 MultiFab& LS_new=get_new_data(LS_Type,slab_step+1);

 int ncomp_state=LS_new.nComp();

 if (ncomp_state==num_materials*(1+AMREX_SPACEDIM)) {
  // do nothing
 } else
  amrex::Error("LS_new.nComp() invalid");

 NavierStokes& ns_level0=getLevel(0);
 My_ParticleContainer& localPC=ns_level0.newDataPC(slab_step+1);

 int number_sweeps=1;

 if (append_flag==OP_PARTICLE_INIT) {
  number_sweeps=2;

  if (slab_step==ns_time_order-1) {
   // do nothing
  } else
   amrex::Error("expecting slab_step==ns_time_order-1");

 } else if (append_flag==OP_PARTICLE_ADD) {
  number_sweeps=2;
 } else if (append_flag==OP_PARTICLE_ASSIMILATE) {
  number_sweeps=1;
 } else
  amrex::Error("append_flag invalid");

 Vector<int> dombc(2*AMREX_SPACEDIM);
 const BCRec& descbc = get_desc_lst()[State_Type].getBC(STATECOMP_MOF);
 const int* b_rec=descbc.vect();
 for (int m=0;m<2*AMREX_SPACEDIM;m++)
  dombc[m]=b_rec[m];

 const Box& domain = geom.Domain();
 const int* domlo = domain.loVect();
 const int* domhi = domain.hiVect();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(lsmf->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*lsmf,use_tiling); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  const int gridno = mfi.index();

  const Box& tilegrid = mfi.tilebox();
  Box tilegrid_grow(grow(tilegrid,1));

  const Box& fabgrid = grids[gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();
  int bfact=parent->Space_blockingFactor(level);

  const Real* xlo = grid_loc[gridno].lo();

  FArrayBox& mfinerfab=(*localMF[MASKCOEF_MF])[mfi];

  FArrayBox& lsfab=(*lsmf)[mfi];

  FArrayBox& lsnewfab=LS_new[mfi];

   // component 1: number of cell (i,j,k) particles.
   // component 2: link to the first particle in the list of (i,j,k) particles.
   //  cell_particle_count(i,j,k,2)=0 => no particles
   //  1<=cell_particle_count(i,j,k,2)<=Number particles
  BaseFab<int> cell_particle_count(tilegrid_grow,2);
  Array4<int> const& cell_particle_count_array=cell_particle_count.array();
  const Dim3 lo3=amrex::lbound(tilegrid_grow);
  const Dim3 hi3=amrex::ubound(tilegrid_grow);
  for (int n=0;n<2;++n) {
  for (int z=lo3.z;z<=hi3.z;++z) {
  for (int y=lo3.y;y<=hi3.y;++y) {
  for (int x=lo3.x;x<=hi3.x;++x) {
   cell_particle_count_array(x,y,z,n)=0;
  }
  }
  }
  }

  if (N_EXTRA_REAL==0) {
   // do nothing
  } else 
   amrex::Error("N_EXTRA_REAL invalid");

  if (N_EXTRA_INT==1) {
   // do nothing
  } else 
   amrex::Error("N_EXTRA_INT invalid");

  // allocate for just one particle for now.
  // after the first sweep, this command is given:
  // new_particle_data.resize(Np_append*single_particle_size);
  int single_particle_size=AMREX_SPACEDIM+N_EXTRA_REAL+N_EXTRA_INT;
  Vector< Real > new_particle_data;
  new_particle_data.resize(single_particle_size);

    // this is an object with a pointer to AoS data
  auto& particles_grid_tile = localPC.GetParticles(level)
    [std::make_pair(mfi.index(),mfi.LocalTileIndex())];
 
  auto& particles_AoS = particles_grid_tile.GetArrayOfStructs();
  unsigned int Np=particles_AoS.size();

  auto& NBR_particles_grid_tile = 
   ns_level0.NBR_Particle_Container->GetParticles(level)
    [std::make_pair(mfi.index(),mfi.LocalTileIndex())];
  auto& NBR_particles_AoS = NBR_particles_grid_tile.GetArrayOfStructs();
  unsigned int NBR_Np=NBR_particles_AoS.size();

  AMREX_ALWAYS_ASSERT(Np<=NBR_Np);

   // The link index will start at 1.
  Vector< int > particle_link_data;
   // i_particle_link_1,i1,j1,k1,   (child link, parent link)
   // i_particle_link_2,i2,j2,k2,  ...
   // if (i,j,k) has particles then:
   // 1. 1<=cell_particle_count(i,j,k,2)<=num particles is link to top of list
   // 2. particle_link_data( (pt_top-1)(1+sdim)+1 )=link to 2nd in list if <>0
  particle_link_data.resize(NBR_Np*(1+AMREX_SPACEDIM));

  for (unsigned int i_link=0;i_link<NBR_Np*(1+AMREX_SPACEDIM);i_link++) {
   particle_link_data[i_link]=0;
  }

  // 1 if particle should be deleted.
  Vector< int > particle_delete_flag; 
  particle_delete_flag.resize(Np);
  for (unsigned int i_delete=0;i_delete<Np;i_delete++) {
   particle_delete_flag[i_delete]=0;
  }

  int Np_append=0;  // number of particles to append

  Vector<int> velbc=getBCArray(State_Type,gridno,
    STATECOMP_VEL,STATE_NCOMP_VEL);

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  for (int isweep=0;isweep<number_sweeps;isweep++) {

   int new_Pdata_size=new_particle_data.size();

   // declared in: LEVELSET_3D.F90
   // 1. subdivide each cell with 
   //    "particle_nsubdivide[bulk|narrow|curvature]" divisions.
   //    e.g. if particle_nsubdivide_bulk=2 => 4 pieces in 2D.
   //                 "         "        =4 => 64 pieces in 2D.
   // 2. for each small sub-box, add a particle at the sub-box center.
   fort_init_particle_container( 
     local_caller_string.c_str(),
     local_caller_string.size(),
     &tid_current,
     velbc.dataPtr(),
     dombc.dataPtr(),
     domlo,domhi,
     &single_particle_size,
     &isweep,
     &number_sweeps,
     &append_flag,
     &particle_nsubdivide,
     &particle_max_per_nsubdivide,
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,
     &level,
     &finest_level,
     &cur_time_slab,
     &dt_slab, //init_particle_container
     xlo,dx,
     &ncomp_state,
     particles_AoS.data(), // existing particles
     NBR_particles_AoS.data(), 
     Np,  // pass by value
     NBR_Np,  // pass by value
     new_particle_data.dataPtr(), // size is "new_Pdata_size"
     &new_Pdata_size,
     &Np_append,  // Np_append number of new particles to add.
     particle_link_data.dataPtr(),
     particle_delete_flag.dataPtr(),
     cell_particle_count.dataPtr(),
     ARLIM(cell_particle_count.loVect()),
     ARLIM(cell_particle_count.hiVect()),
     lsfab.dataPtr(),
     ARLIM(lsfab.loVect()),
     ARLIM(lsfab.hiVect()),
     lsnewfab.dataPtr(),
     ARLIM(lsnewfab.loVect()),ARLIM(lsnewfab.hiVect()),
     mfinerfab.dataPtr(),
     ARLIM(mfinerfab.loVect()),ARLIM(mfinerfab.hiVect()));
    
   if (isweep==0) {
    new_particle_data.resize(Np_append*single_particle_size);
   }
  } // isweep=0,...,number_sweeps-1

  unsigned int Np_delete=0;
  for (unsigned int i_delete=0;i_delete<Np;i_delete++) {
   if (particle_delete_flag[i_delete]==1) {
    Np_delete++;
   } else if (particle_delete_flag[i_delete]==0) {
    // do nothing
   } else
    amrex::Error("particle_delete_flag[i_delete] invalid");
  }

  AMREX_ALWAYS_ASSERT(Np_delete<=Np);

  if ((append_flag==OP_PARTICLE_INIT)||
      (append_flag==OP_PARTICLE_ADD)) {

   Vector< My_ParticleContainer::ParticleType > mirrorPC_AoS;
   unsigned int Np_mirror_AoS=Np-Np_delete+Np_append;
   mirrorPC_AoS.resize(Np_mirror_AoS);

   //save the existing particle data to:
   // mirrorPC_AoS
   unsigned int i_mirror=0;
   for (unsigned int i_delete=0;i_delete<Np;i_delete++) {
    if (particle_delete_flag[i_delete]==1) {
     // do nothing
    } else if (particle_delete_flag[i_delete]==0) {
     mirrorPC_AoS[i_mirror]=particles_AoS[i_delete];
     i_mirror++;
    } else
     amrex::Error("particle_delete_flag[i_delete] invalid");
   } // for (unsigned int i_delete=0;i_delete<Np;i_delete++) 

   AMREX_ALWAYS_ASSERT(i_mirror==Np-Np_delete);

   for (int i_append=0;i_append<Np_append;i_append++) {

    My_ParticleContainer::ParticleType p;
    p.id() = My_ParticleContainer::ParticleType::NextID();
 
    p.cpu() = ParallelDescriptor::MyProc();

    int ibase=i_append*single_particle_size;
    //pos(AMREX_SPACEDIM)
    for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
     p.pos(dir) = new_particle_data[ibase+dir];
    }
    for (int dir=0;dir<N_EXTRA_REAL;dir++) {
     p.rdata(dir) = new_particle_data[ibase+AMREX_SPACEDIM+dir];
    }
    //material id.
    for (int dir=0;dir<N_EXTRA_INT;dir++) {
     p.idata(dir) = 
      (int) new_particle_data[ibase+AMREX_SPACEDIM+N_EXTRA_REAL+dir];
    }
    mirrorPC_AoS[i_mirror]=p;
    i_mirror++;
   } // i_append=0..Np_append-1

   AMREX_ALWAYS_ASSERT(i_mirror==Np_mirror_AoS);

   particles_grid_tile.resize(0);

   for (i_mirror=0;i_mirror<Np_mirror_AoS;i_mirror++) {
    particles_grid_tile.push_back(mirrorPC_AoS[i_mirror]);
   }

  } else if (append_flag==OP_PARTICLE_ASSIMILATE) {
   //do nothing
  } else
   amrex::Error("append_flag invalid");

 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_INIT_PARTICLE_CONTAINER,"init_particle_container");

 delete lsmf;

}  // end subroutine init_particle_container()

#endif


// should be cur_time=0 and prev_time=-1
// called from post_init
void
NavierStokes::post_init_state () {
    
 std::string local_caller_string="post_init_state";

 if (level==0) {
  // do nothing
 } else
  amrex::Error("require level==0 in post_init_state");

 SDC_setup();
 ns_time_order=parent->Time_blockingFactor();
 slab_step=ns_time_order-1; 

 SDC_outer_sweeps=0;
 SDC_setup_step();
 dt_slab=1.0;
 delta_slab_time=1.0;

 Number_CellsALL(real_number_of_cells);

 int combine_flag=2;  
 int hflag=0;
 int combine_idx=-1;
 int update_flux=0;
 int interface_cond_avail=0;

 const int finest_level = parent->finestLevel();

   // inside of post_init_state

   // metrics_data
   // allocate_mdot (MDOT_MF=0.0)
   // MASKCOEF
   // init_FSI_GHOST_MAC_MF
   // VOF_Recon_ALL (RECON_UPDATE_STATE_ERR_AND_CENTROID)
   // makeStateDistALL
   // prescribe_solid_geometryALL
   // make_physics_varsALL
 prepare_post_process(local_caller_string);

#ifdef AMREX_PARTICLES

 init_particle_containerALL(OP_PARTICLE_INIT,local_caller_string);

 NavierStokes& ns_level0=getLevel(0);
 My_ParticleContainer& localPC=ns_level0.newDataPC(slab_step+1);

 int lev_min=0;
 int lev_max=-1;
 int nGrow_Redistribute=0;
  //particles are being redistributed for the first time.
 int local_Redistribute=0;
 localPC.Redistribute(lev_min,lev_max,nGrow_Redistribute,local_Redistribute);

#endif

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);

  combine_flag=2; //combine if vfrac<VOFTOL 
  hflag=0;
  combine_idx=-1;
  update_flux=0;
  interface_cond_avail=0;

  ns_level.combine_state_variable(
    SOLVETYPE_HEAT,
    combine_idx,
    combine_flag,
    hflag,
    update_flux,
    interface_cond_avail); 

  for (int ns=0;ns<num_species_var;ns++) {
   ns_level.combine_state_variable(
    SOLVETYPE_SPEC+ns,
    combine_idx,
    combine_flag,
    hflag,
    update_flux,
    interface_cond_avail); 
  }

 } // ilev=finest_level ... level

 Vector<blobclass> blobdata;
 Vector< Vector<Real> > mdot_data;
 Vector< Vector<Real> > mdot_comp_data;
 Vector< Vector<Real> > mdot_data_redistribute;
 Vector< Vector<Real> > mdot_comp_data_redistribute;
 Vector<int> type_flag;

 int color_count=0;
 int coarsest_level=0;
 int idx_mdot=-1; //idx_mdot==-1 => do not collect auxiliary data.
 int tessellate=1;
 int operation_flag=OP_GATHER_MDOT;

 ColorSumALL(
  operation_flag, //=OP_GATHER_MDOT
  tessellate,  //=1
  coarsest_level,
  color_count,
  TYPE_MF,
  COLOR_MF,
  idx_mdot,
  idx_mdot,
  type_flag,
  blobdata,
  mdot_data,
  mdot_comp_data,
  mdot_data_redistribute,
  mdot_comp_data_redistribute
  );

 if (color_count!=blobdata.size())
  amrex::Error("color_count!=blobdata.size()");

 if (is_zalesak()) {
  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   // init cell center gas/liquid velocity for all materials.
   ns_level.zalesakVEL(); 
  }
 } else if (! is_zalesak()) {
   // do nothing
 } else
  amrex::Error("is_zalesak() invalid");

  // unew^{f} = unew^{c->f}
  // in post_init_state 
  // if project_option==SOLVETYPE_INITPROJ, then the velocity in the ice
  // is overwritten with a projected rigid body velocity.
 operation_flag=OP_UNEW_CELL_TO_MAC;
 int idx_velcell=-1;
 Real beta=0.0;

 increment_face_velocityALL(
   operation_flag,
   SOLVETYPE_INITPROJ,
   idx_velcell,
   beta,blobdata); 

 delete_array(TYPE_MF);
 delete_array(COLOR_MF);

 if (step_through_data==1) {
  int basestep_debug=nStep();
  parent->writeDEBUG_PlotFile(
    basestep_debug,
    SDC_outer_sweeps,
    slab_step,
    divu_outer_sweeps);
  std::cout << "press any number then enter (prior post_init_pressure) \n";
  int n_input;
  std::cin >> n_input;
 }

 if ((post_init_pressure_solve==1)&&
     (!is_zalesak())) { 

  double start_pressure_solve = ParallelDescriptor::second();

   // U^CELL and U^MAC; assimilates the solid velocity (MAC and CELL)
   // and  ice velocity (MAC)
  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   combine_flag=2;
   hflag=0;
   combine_idx=-1;
   update_flux=0;
   interface_cond_avail=0;
   ns_level.combine_state_variable(
    SOLVETYPE_VISC,
    combine_idx,
    combine_flag,
    hflag,
    update_flux, 
    interface_cond_avail);
   update_flux=1;
   ns_level.combine_state_variable(
    SOLVETYPE_PRES,
    combine_idx,
    combine_flag,
    hflag,
    update_flux,
    interface_cond_avail);

   ns_level.make_MAC_velocity_consistent();
  }

  multiphase_project(SOLVETYPE_INITPROJ); 

   // U^CELL and U^MAC
  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);

   combine_flag=2;
   hflag=0;
   combine_idx=-1;
   update_flux=0;
   interface_cond_avail=0;
   ns_level.combine_state_variable(
    SOLVETYPE_VISC,
    combine_idx,
    combine_flag,
    hflag,
    update_flux,
    interface_cond_avail);
   update_flux=1;
   ns_level.combine_state_variable(
    SOLVETYPE_PRES,
    combine_idx,
    combine_flag,
    hflag,
    update_flux,
    interface_cond_avail);

   ns_level.make_MAC_velocity_consistent();
  }

  double end_pressure_solve = ParallelDescriptor::second();
  if (verbose>0)
   if (ParallelDescriptor::IOProcessor())
    std::cout << "initial pressure solve time " << end_pressure_solve-
       start_pressure_solve << '\n';

 } else if ((post_init_pressure_solve==0)||(is_zalesak())) {
  // do nothing
 } else {
  amrex::Error("post_init_pressure_solve or is_zalesak() invalid");
 } 

 for (int ilev=finest_level;ilev>=0;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.avgDown(State_Type,STATECOMP_PRES,1,1);
  ns_level.avgDown(State_Type,STATECOMP_STATES,
     num_state_material*num_materials,1);
 } // ilev

 delete_array(MASKCOEF_MF);
 delete_array(MASK_NBR_MF);

 CopyNewToOldALL();

 if (step_through_data==1) {
  int basestep_debug=nStep();
  parent->writeDEBUG_PlotFile(
    basestep_debug,
    SDC_outer_sweeps,
    slab_step,
    divu_outer_sweeps);
  std::cout << "press any number then enter: post_init_state\n";
  int n_input;
  std::cin >> n_input;
 }

}  // end subroutine post_init_state

void
NavierStokes::level_avgDown_tag(MultiFab& S_crse,MultiFab& S_fine) {

 int scomp=0;
 int ncomp=1;

 int finest_level=parent->finestLevel();
 if (level>=finest_level)
  amrex::Error("level invalid in level_avgDown_tag");

 int f_level=level+1;
 NavierStokes&   fine_lev = getLevel(f_level);
 const BoxArray& fgrids=fine_lev.grids;
 const DistributionMapping& fdmap=fine_lev.dmap;

 if (grids!=S_crse.boxArray())
  amrex::Error("S_crse invalid level_avgDown_tag");
 if (fgrids!=S_fine.boxArray())
  amrex::Error("S_fine invalid");
 if (S_crse.nComp()!=S_fine.nComp())
  amrex::Error("nComp mismatch");
 if (S_crse.nComp()<scomp+ncomp)
  amrex::Error("S_crse.nComp() invalid level_avgDown_tag");

 BoxArray crse_S_fine_BA(fgrids.size());
 for (int i = 0; i < fgrids.size(); ++i) {
  crse_S_fine_BA.set(i,amrex::coarsen(fgrids[i],2));
 }
 DistributionMapping crse_dmap=fdmap;
 MultiFab crse_S_fine(crse_S_fine_BA,crse_dmap,ncomp,0,
	MFInfo().SetTag("crse_S_fine"),FArrayBoxFactory());

 ParallelDescriptor::Barrier();

 const Real* dx = geom.CellSize();
 const Real* dxf = fine_lev.geom.CellSize();
 const Real* prob_lo   = geom.ProbLo();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(S_fine.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(S_fine,false); mfi.isValid(); ++mfi) {
  BL_ASSERT(fgrids[mfi.index()] == mfi.validbox());
  const Box& tilegrid = mfi.tilebox();

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  const int gridno = mfi.index();
  const Real* xlo_fine = fine_lev.grid_loc[gridno].lo();
  const Box& ovgrid = crse_S_fine_BA[gridno];
  const int* ovlo=ovgrid.loVect();
  const int* ovhi=ovgrid.hiVect();

  FArrayBox& fine_fab=S_fine[gridno];
  const Box& fgrid=fine_fab.box();
  const int* flo=fgrid.loVect();
  const int* fhi=fgrid.hiVect();
  const Real* f_dat=fine_fab.dataPtr(scomp);

  FArrayBox& coarse_fab=crse_S_fine[gridno];
  const Box& cgrid = coarse_fab.box();
  const int* clo=cgrid.loVect();
  const int* chi=cgrid.hiVect();
  const Real* c_dat=coarse_fab.dataPtr();

  int bfact_c=parent->Space_blockingFactor(level);
  int bfact_f=parent->Space_blockingFactor(f_level);

  const Box& fine_fabgrid = fine_lev.grids[gridno];
  const int* fine_fablo=fine_fabgrid.loVect();
  const int* fine_fabhi=fine_fabgrid.hiVect();

  fort_avgdown_tag(
   prob_lo,
   dxf,
   &level,&f_level,
   &bfact_c,&bfact_f,
   xlo_fine,dx,
   &ncomp,
   c_dat,ARLIM(clo),ARLIM(chi),
   f_dat,ARLIM(flo),ARLIM(fhi),
   ovlo,ovhi,
   fine_fablo,fine_fabhi);

 }// mfi
} //omp
 ns_reconcile_d_num(LOOP_AVGDOWN_TAG,"level_avgDown_tag");
 S_crse.ParallelCopy(crse_S_fine,0,scomp,ncomp);
 ParallelDescriptor::Barrier();
} // subroutine level_avgDown_tag


void
NavierStokes::level_avgDownBURNING(MultiFab& S_crse,MultiFab& S_fine, 
		int velflag) {

 int nburning=EXTRAP_NCOMP_BURNING;
 int ntsat=EXTRAP_NCOMP_TSAT;
 int scomp=0;
 int ncomp=0;
 if (velflag==1) {
  ncomp=nburning;
 } else if (velflag==0) {
  ncomp=ntsat;
 } else
  amrex::Error("velflag invalid");

 int finest_level=parent->finestLevel();
 if (level>=finest_level)
  amrex::Error("level invalid in level_avgDownBURNING");

 int f_level=level+1;
 NavierStokes&   fine_lev = getLevel(f_level);
 const BoxArray& fgrids=fine_lev.grids;
 const DistributionMapping& fdmap=fine_lev.dmap;

 if (grids!=S_crse.boxArray())
  amrex::Error("S_crse invalid level_avgDownBURNING");
 if (fgrids!=S_fine.boxArray())
  amrex::Error("S_fine invalid");
 if (S_crse.nComp()!=S_fine.nComp())
  amrex::Error("nComp mismatch");
 if (S_crse.nComp()!=scomp+ncomp)
  amrex::Error("S_crse.nComp() invalid level_avgDownBurning");

 BoxArray crse_S_fine_BA(fgrids.size());
 for (int i = 0; i < fgrids.size(); ++i) {
  crse_S_fine_BA.set(i,amrex::coarsen(fgrids[i],2));
 }
 DistributionMapping crse_dmap=fdmap;
 MultiFab crse_S_fine(crse_S_fine_BA,crse_dmap,ncomp,0,
   MFInfo().SetTag("crse_S_fine"),FArrayBoxFactory());

 ParallelDescriptor::Barrier();

 const Real* dx = geom.CellSize();
 const Real* dxf = fine_lev.geom.CellSize();
 const Real* prob_lo   = geom.ProbLo();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(S_fine.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(S_fine,false); mfi.isValid(); ++mfi) {
  BL_ASSERT(fgrids[mfi.index()] == mfi.validbox());
  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  const Real* xlo_fine = fine_lev.grid_loc[gridno].lo();
  const Box& ovgrid = crse_S_fine_BA[gridno];
  const int* ovlo=ovgrid.loVect();
  const int* ovhi=ovgrid.hiVect();

  FArrayBox& fine_fab=S_fine[gridno];
  const Box& fgrid=fine_fab.box();
  const int* flo=fgrid.loVect();
  const int* fhi=fgrid.hiVect();
  const Real* f_dat=fine_fab.dataPtr(scomp);

  FArrayBox& coarse_fab=crse_S_fine[gridno];
  const Box& cgrid = coarse_fab.box();
  const int* clo=cgrid.loVect();
  const int* chi=cgrid.hiVect();
  const Real* c_dat=coarse_fab.dataPtr();

  int bfact_c=parent->Space_blockingFactor(level);
  int bfact_f=parent->Space_blockingFactor(f_level);

  const Box& fine_fabgrid = fine_lev.grids[gridno];
  const int* fine_fablo=fine_fabgrid.loVect();
  const int* fine_fabhi=fine_fabgrid.hiVect();

   // in: NAVIERSTOKES_3D.F90
  fort_avgdown_burning(
   &velflag,
   prob_lo,
   dxf,
   &level,&f_level,
   &bfact_c,&bfact_f,
   xlo_fine,dx,
   &ncomp,
   c_dat,ARLIM(clo),ARLIM(chi),
   f_dat,ARLIM(flo),ARLIM(fhi),
   ovlo,ovhi,
   fine_fablo,fine_fabhi);

 }// mfi
} //omp
 ns_reconcile_d_num(LOOP_AVGDOWN_BURNING,"level_avgDownBURNING");
 S_crse.ParallelCopy(crse_S_fine,0,scomp,ncomp);
 ParallelDescriptor::Barrier();

} // end subroutine level_avgDownBURNING

void
NavierStokes::level_avgDownDRAG(MultiFab& S_crse,MultiFab& S_fine) {

 int scomp=0;
 int ncomp=N_DRAG;

 int finest_level=parent->finestLevel();
 if (level>=finest_level)
  amrex::Error("level invalid in level_avgDownDRAG");

 int f_level=level+1;
 NavierStokes&   fine_lev = getLevel(f_level);
 const BoxArray& fgrids=fine_lev.grids;
 const DistributionMapping& fdmap=fine_lev.dmap;

 if (grids!=S_crse.boxArray())
  amrex::Error("S_crse invalid level_avgDownDRAG");
 if (fgrids!=S_fine.boxArray())
  amrex::Error("S_fine invalid");
 if (S_crse.nComp()!=S_fine.nComp())
  amrex::Error("nComp mismatch");
 if (S_crse.nComp()!=scomp+ncomp)
  amrex::Error("S_crse.nComp() invalid level_avgDownDRAG");

 BoxArray crse_S_fine_BA(fgrids.size());
 for (int i = 0; i < fgrids.size(); ++i) {
  crse_S_fine_BA.set(i,amrex::coarsen(fgrids[i],2));
 }
 DistributionMapping crse_dmap=fdmap;
 MultiFab crse_S_fine(crse_S_fine_BA,crse_dmap,ncomp,0,
   MFInfo().SetTag("crse_S_fine"),FArrayBoxFactory());

 ParallelDescriptor::Barrier();

 const Real* dx = geom.CellSize();
 const Real* dxf = fine_lev.geom.CellSize();
 const Real* prob_lo   = geom.ProbLo();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(S_fine.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(S_fine,false); mfi.isValid(); ++mfi) {
  BL_ASSERT(fgrids[mfi.index()] == mfi.validbox());
  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  const Real* xlo_fine = fine_lev.grid_loc[gridno].lo();
  const Box& ovgrid = crse_S_fine_BA[gridno];
  const int* ovlo=ovgrid.loVect();
  const int* ovhi=ovgrid.hiVect();

  FArrayBox& fine_fab=S_fine[gridno];
  const Box& fgrid=fine_fab.box();
  const int* flo=fgrid.loVect();
  const int* fhi=fgrid.hiVect();
  const Real* f_dat=fine_fab.dataPtr(scomp);

  FArrayBox& coarse_fab=crse_S_fine[gridno];
  const Box& cgrid = coarse_fab.box();
  const int* clo=cgrid.loVect();
  const int* chi=cgrid.hiVect();
  const Real* c_dat=coarse_fab.dataPtr();

  int bfact_c=parent->Space_blockingFactor(level);
  int bfact_f=parent->Space_blockingFactor(f_level);

  const Box& fine_fabgrid = fine_lev.grids[gridno];
  const int* fine_fablo=fine_fabgrid.loVect();
  const int* fine_fabhi=fine_fabgrid.hiVect();

   // in: NAVIERSTOKES_3D.F90
  fort_avgdown_drag(
   prob_lo,
   dxf,
   &level,&f_level,
   &bfact_c,&bfact_f,
   xlo_fine,dx,
   &ncomp,
   c_dat,ARLIM(clo),ARLIM(chi),
   f_dat,ARLIM(flo),ARLIM(fhi),
   ovlo,ovhi,
   fine_fablo,fine_fabhi);

 }// mfi
} //omp
 ns_reconcile_d_num(LOOP_AVGDOWN_DRAG,"level_avgDownDRAG");
 S_crse.ParallelCopy(crse_S_fine,0,scomp,ncomp);
 ParallelDescriptor::Barrier();

} // end subroutine level_avgDownDRAG


void
NavierStokes::level_avgDownCURV(MultiFab& S_crse,MultiFab& S_fine) {

 int scomp=0;
 int ncomp_curv_total=num_interfaces*CURVCOMP_NCOMP;

 int finest_level=parent->finestLevel();
 if (level>=finest_level)
  amrex::Error("level invalid in level_avgDownCURV");

 int f_level=level+1;
 NavierStokes&   fine_lev = getLevel(f_level);
 const BoxArray& fgrids=fine_lev.grids;
 const DistributionMapping& fdmap=fine_lev.dmap;

 if (grids!=S_crse.boxArray())
  amrex::Error("S_crse invalid level_avgDownCURV");
 if (fgrids!=S_fine.boxArray())
  amrex::Error("S_fine invalid");
 if (S_crse.nComp()!=S_fine.nComp())
  amrex::Error("nComp mismatch");
 if (S_crse.nComp()!=scomp+ncomp_curv_total)
  amrex::Error("S_crse.nComp() invalid level_avgDownCURV");

 BoxArray crse_S_fine_BA(fgrids.size());
 for (int i = 0; i < fgrids.size(); ++i) {
  crse_S_fine_BA.set(i,amrex::coarsen(fgrids[i],2));
 }
 DistributionMapping crse_dmap=fdmap;
 MultiFab crse_S_fine(crse_S_fine_BA,crse_dmap,ncomp_curv_total,0,
   MFInfo().SetTag("crse_S_fine"),FArrayBoxFactory());

 ParallelDescriptor::Barrier();

 const Real* dx = geom.CellSize();
 const Real* dxf = fine_lev.geom.CellSize();
 const Real* prob_lo   = geom.ProbLo();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(S_fine.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(S_fine,false); mfi.isValid(); ++mfi) {
  BL_ASSERT(fgrids[mfi.index()] == mfi.validbox());
  const Box& tilegrid = mfi.tilebox();

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  const int gridno = mfi.index();
  const Real* xlo_fine = fine_lev.grid_loc[gridno].lo();
  const Box& ovgrid = crse_S_fine_BA[gridno];
  const int* ovlo=ovgrid.loVect();
  const int* ovhi=ovgrid.hiVect();

  FArrayBox& fine_fab=S_fine[gridno];
  const Box& fgrid=fine_fab.box();
  const int* flo=fgrid.loVect();
  const int* fhi=fgrid.hiVect();
  const Real* f_dat=fine_fab.dataPtr(scomp);

  FArrayBox& coarse_fab=crse_S_fine[gridno];
  const Box& cgrid = coarse_fab.box();
  const int* clo=cgrid.loVect();
  const int* chi=cgrid.hiVect();
  const Real* c_dat=coarse_fab.dataPtr();

  int bfact_c=parent->Space_blockingFactor(level);
  int bfact_f=parent->Space_blockingFactor(f_level);

  const Box& fine_fabgrid = fine_lev.grids[gridno];
  const int* fine_fablo=fine_fabgrid.loVect();
  const int* fine_fabhi=fine_fabgrid.hiVect();

  fort_avgdown_curv(
   prob_lo,
   dxf,
   &level,&f_level,
   &bfact_c,&bfact_f,
   xlo_fine,dx,
   &ncomp_curv_total,
   c_dat,ARLIM(clo),ARLIM(chi),
   f_dat,ARLIM(flo),ARLIM(fhi),
   ovlo,ovhi,
   fine_fablo,fine_fabhi);

 }// mfi
} //omp
 ns_reconcile_d_num(LOOP_AVGDOWN_CURV,"level_avgDownCURV");
 S_crse.ParallelCopy(crse_S_fine,0,scomp,ncomp_curv_total);
 ParallelDescriptor::Barrier();

} // subroutine level_avgDownCURV


// spectral_override==0 => always low order.
void
NavierStokes::avgDown(MultiFab& S_crse,MultiFab& S_fine,
  int scomp,int ncomp,int spectral_override) {

 int finest_level=parent->finestLevel();
 if (level>=finest_level)
  amrex::Error("level invalid in avgDown");

 int f_level=level+1;
 NavierStokes&   fine_lev = getLevel(f_level);
 const BoxArray& fgrids=fine_lev.grids;
 const DistributionMapping& fdmap=fine_lev.dmap;

 if (grids!=S_crse.boxArray()) {
  std::cout << "scomp=" << scomp << "ncomp=" << ncomp <<
   "spectral_override=" << spectral_override << '\n';
  amrex::Error("S_crse invalid avgDown");
 }
 if (fgrids!=S_fine.boxArray())
  amrex::Error("S_fine invalid");
 if (S_crse.nComp()!=S_fine.nComp())
  amrex::Error("nComp mismatch");

 BoxArray crse_S_fine_BA(fgrids.size());
 for (int i = 0; i < fgrids.size(); ++i) {
  crse_S_fine_BA.set(i,amrex::coarsen(fgrids[i],2));
 }
 DistributionMapping crse_dmap=fdmap;
 MultiFab crse_S_fine(crse_S_fine_BA,crse_dmap,ncomp,0,
  MFInfo().SetTag("crse_S_fine"),FArrayBoxFactory());

 ParallelDescriptor::Barrier();

 const Real* dx = geom.CellSize();
 const Real* dxf = fine_lev.geom.CellSize();
 const Real* prob_lo   = geom.ProbLo();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(S_fine.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(S_fine,false); mfi.isValid(); ++mfi) {
  BL_ASSERT(fgrids[mfi.index()] == mfi.validbox());
  const Box& tilegrid = mfi.tilebox();

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  const int gridno = mfi.index();
  const Real* xlo_fine = fine_lev.grid_loc[gridno].lo();
  const Box& ovgrid = crse_S_fine_BA[gridno];
  const int* ovlo=ovgrid.loVect();
  const int* ovhi=ovgrid.hiVect();

  FArrayBox& fine_fab=S_fine[gridno];
  const Box& fgrid=fine_fab.box();
  const int* flo=fgrid.loVect();
  const int* fhi=fgrid.hiVect();
  const Real* f_dat=fine_fab.dataPtr(scomp);

  FArrayBox& coarse_fab=crse_S_fine[gridno];
  const Box& cgrid = coarse_fab.box();
  const int* clo=cgrid.loVect();
  const int* chi=cgrid.hiVect();
  const Real* c_dat=coarse_fab.dataPtr();

  int bfact_c=parent->Space_blockingFactor(level);
  int bfact_f=parent->Space_blockingFactor(f_level);

  const Box& fine_fabgrid = fine_lev.grids[gridno];
  const int* fine_fablo=fine_fabgrid.loVect();
  const int* fine_fabhi=fine_fabgrid.hiVect();

  if (spectral_override==1) {

   FArrayBox& maskfab=(*fine_lev.localMF[MASKSEM_MF])[mfi];
   fort_avgdown( 
    &enable_spectral,
    &finest_level,
    &spectral_override,
    prob_lo,
    dxf,
    &level,&f_level,
    &bfact_c,&bfact_f,     
    xlo_fine,dx,
    &ncomp,
    c_dat,ARLIM(clo),ARLIM(chi),
    f_dat,ARLIM(flo),ARLIM(fhi),
    maskfab.dataPtr(),
    ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
    ovlo,ovhi, 
    fine_fablo,fine_fabhi);

  } else if (spectral_override==0) {

   fort_avgdown_low(
    prob_lo,
    dxf,
    &level,&f_level,
    &bfact_c,&bfact_f,
    xlo_fine,dx,
    &ncomp,
    c_dat,ARLIM(clo),ARLIM(chi),
    f_dat,ARLIM(flo),ARLIM(fhi),
    ovlo,ovhi,
    fine_fablo,fine_fabhi);

  } else
   amrex::Error("spectral_override invalid");

 }// mfi
} //omp
 ns_reconcile_d_num(LOOP_AVGDOWN_OR_AVGDOWN_LOW,"avgDown");
 S_crse.ParallelCopy(crse_S_fine,0,scomp,ncomp);
 ParallelDescriptor::Barrier();
}  // subroutine avgDown

// spectral_override==1 => order derived from "enable_spectral"
// spectral_override==0 => always low order.
void 
NavierStokes::avgDown_list(int stateidx,Vector<int> scomp,
   Vector<int> ncomp,int spectral_override) {

 int ncomp_list=0;
 for (int ilist=0;ilist<scomp.size();ilist++) {
  avgDown(stateidx,scomp[ilist],ncomp[ilist],spectral_override);
  ncomp_list+=ncomp[ilist];
 }
 if (ncomp_list<=0)
  amrex::Error("ncomp_list invalid");

} // subroutine avgDown_list

// spectral_override==1 => order derived from "enable_spectral"
// spectral_override==0 => always low order.
void
NavierStokes::avgDown(int stateidx,int startcomp,int numcomp,
 int spectral_override) {

 if ((stateidx!=LS_Type)&&
     (stateidx!=State_Type)&&
     (stateidx!=DIV_Type)&&
     (stateidx!=Tensor_Type)) {
  std::cout << "stateidx= " << stateidx << '\n';
  std::cout << "startcomp= " << startcomp << '\n';
  std::cout << "numcomp= " << numcomp << '\n';
  amrex::Error("stateidx invalid");
 }

 int finest_level=parent->finestLevel();

 if (level==finest_level) {
  // do nothing
 } else if ((level>=0)&&(level<finest_level)) {

  NavierStokes&   fine_lev = getLevel(level+1);
  MultiFab& S_crse = get_new_data(stateidx,slab_step+1);
  MultiFab& S_fine = fine_lev.get_new_data(stateidx,slab_step+1);
  avgDown(S_crse,S_fine,startcomp,numcomp,spectral_override);

 } else
  amrex::Error("level invalid");

} // subroutine avgDown



void NavierStokes::MOFavgDown() {

 std::string local_caller_string="MOFavgDown";

 int finest_level=parent->finestLevel();

 if (level == finest_level)
  return;

 int f_level=level+1;
 NavierStokes&   fine_lev = getLevel(f_level);
 const BoxArray& fgrids=fine_lev.grids;
 const DistributionMapping& fdmap=fine_lev.dmap;
 resize_metrics(1);
 debug_ngrow(VOLUME_MF,0,local_caller_string);
 fine_lev.resize_metrics(1);
 fine_lev.debug_ngrow(VOLUME_MF,1,local_caller_string);

 MultiFab& S_fine=fine_lev.get_new_data(State_Type,slab_step+1);
 MultiFab& S_crse = get_new_data(State_Type,slab_step+1);

 const Real* dxf = fine_lev.geom.CellSize();
 const Real* dxc = geom.CellSize();
 const Real* prob_lo   = geom.ProbLo();

 if (grids!=S_crse.boxArray())
  amrex::Error("S_crse invalid MOFavgDown");
 if (fgrids!=S_fine.boxArray())
  amrex::Error("S_fine invalid");
 if (S_crse.nComp()!=S_fine.nComp())
  amrex::Error("nComp mismatch");

 BoxArray crse_S_fine_BA(fgrids.size());
 for (int i = 0; i < fgrids.size(); ++i) {
  crse_S_fine_BA.set(i,amrex::coarsen(fgrids[i],2));
 }

 DistributionMapping crse_dmap=fdmap;
 MultiFab crse_S_fine(crse_S_fine_BA,crse_dmap,num_materials*ngeom_raw,0,
   MFInfo().SetTag("crse_S_fine"),FArrayBoxFactory());

 ParallelDescriptor::Barrier();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(S_fine.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(S_fine,false); mfi.isValid(); ++mfi) {
  BL_ASSERT(fgrids[mfi.index()] == mfi.validbox());
  const Box& tilegrid = mfi.tilebox();

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  const int gridno = mfi.index();
  const Box& ovgrid = crse_S_fine_BA[gridno];
  const int* ovlo=ovgrid.loVect();
  const int* ovhi=ovgrid.hiVect();
  
  FArrayBox& finefab=S_fine[gridno];
  const Box& fgrid=finefab.box();
  const int* flo=fgrid.loVect();
  const int* fhi=fgrid.hiVect();
  const Real* f_dat=finefab.dataPtr(STATECOMP_MOF);

  FArrayBox& coarsefab=crse_S_fine[gridno];
  const Box& cgrid = coarsefab.box();
  const int* clo=cgrid.loVect();
  const int* chi=cgrid.hiVect();
  const Real* c_dat=coarsefab.dataPtr();

  int bfact_c=parent->Space_blockingFactor(level);
  int bfact_f=parent->Space_blockingFactor(f_level);

  fort_mofavgdown(
   &cur_time_slab,
   prob_lo,
   dxc,
   dxf,
   &bfact_c,&bfact_f,
   c_dat,ARLIM(clo),ARLIM(chi),
   f_dat,ARLIM(flo),ARLIM(fhi),
   ovlo,ovhi);
 } // mfi
} //omp
 ns_reconcile_d_num(LOOP_MOFAVGDOWN,"MOFavgDown");
 S_crse.ParallelCopy(crse_S_fine,0,STATECOMP_MOF,num_materials*ngeom_raw);
 ParallelDescriptor::Barrier();
}


void NavierStokes::avgDownError() {

 std::string local_caller_string="avgDownError";

 int finest_level=parent->finestLevel();

 if (level == finest_level)
  return;

 int f_level=level+1;
 NavierStokes&   fine_lev = getLevel(f_level);
 const BoxArray& fgrids=fine_lev.grids;
 const DistributionMapping& fdmap=fine_lev.dmap;
 resize_metrics(1);
 debug_ngrow(VOLUME_MF,0,local_caller_string);
 fine_lev.resize_metrics(1);
 fine_lev.debug_ngrow(VOLUME_MF,1,local_caller_string);

 MultiFab& S_fine=fine_lev.get_new_data(State_Type,slab_step+1);
 MultiFab& S_crse = get_new_data(State_Type,slab_step+1);

 const Real* dxf = fine_lev.geom.CellSize();
 const Real* prob_lo   = geom.ProbLo();

 if (grids!=S_crse.boxArray())
  amrex::Error("S_crse invalid avgDownError()");
 if (fgrids!=S_fine.boxArray())
  amrex::Error("S_fine invalid");
 if (S_crse.nComp()!=S_fine.nComp())
  amrex::Error("nComp mismatch");

 BoxArray crse_S_fine_BA(fgrids.size());
 for (int i = 0; i < fgrids.size(); ++i) {
  crse_S_fine_BA.set(i,amrex::coarsen(fgrids[i],2));
 }

 if (S_crse.nComp()!=STATE_NCOMP)
  amrex::Error("S_crse.nComp()!=STATE_NCOMP");

 DistributionMapping crse_dmap=fdmap;
 MultiFab crse_S_fine(crse_S_fine_BA,crse_dmap,1,0,
	MFInfo().SetTag("crse_S_fine"),FArrayBoxFactory());

 ParallelDescriptor::Barrier();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(S_fine.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(S_fine,false); mfi.isValid(); ++mfi) {
  BL_ASSERT(fgrids[mfi.index()] == mfi.validbox());
  const Box& tilegrid = mfi.tilebox();

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  const int gridno = mfi.index();
  const Box& ovgrid = crse_S_fine_BA[gridno];
  const int* ovlo=ovgrid.loVect();
  const int* ovhi=ovgrid.hiVect();
  
  FArrayBox& finefab=S_fine[gridno];
  const Box& fgrid=finefab.box();
  const int* flo=fgrid.loVect();
  const int* fhi=fgrid.hiVect();
  const Real* f_dat=finefab.dataPtr(STATECOMP_ERR);

  FArrayBox& coarsefab=crse_S_fine[gridno];
  const Box& cgrid = coarsefab.box();
  const int* clo=cgrid.loVect();
  const int* chi=cgrid.hiVect();
  const Real* c_dat=coarsefab.dataPtr();

  int bfact_c=parent->Space_blockingFactor(level);
  int bfact_f=parent->Space_blockingFactor(f_level);

  fort_erroravgdown(
   prob_lo,
   dxf,
   &bfact_c,&bfact_f,
   c_dat,ARLIM(clo),ARLIM(chi),
   f_dat,ARLIM(flo),ARLIM(fhi),
   ovlo,ovhi);
 } // mfi
} //omp
 ns_reconcile_d_num(LOOP_ERRORAVGDOWN,"avgDownError");
 S_crse.ParallelCopy(crse_S_fine,0,STATECOMP_ERR,1);
 ParallelDescriptor::Barrier();
} // subroutine avgDownError


void NavierStokes::getBCArray_list(Vector<int>& listbc,int state_index,
     int gridno,Vector<int> scomp,Vector<int> ncomp) {

 int ncomp_list=0;
 for (int ilist=0;ilist<scomp.size();ilist++)
  ncomp_list+=ncomp[ilist];
 if (ncomp_list<=0)
  amrex::Error("ncomp_list invalid");

 listbc.resize(ncomp_list*AMREX_SPACEDIM*2);

 int dcomp=0;
 for (int ilist=0;ilist<scomp.size();ilist++) {
  Vector<int> bc_single=getBCArray(state_index,gridno,
    scomp[ilist],ncomp[ilist]);

  int bc_scomp=0;
  for (int nn=0;nn<ncomp[ilist];nn++) {
   for (int side=0;side<=1;side++) {
    for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
     listbc[dcomp]=bc_single[bc_scomp];
     bc_scomp++;
     dcomp++;
    } // dir
   } // side
  } // nn
 } // ilist
 if (dcomp!=ncomp_list*AMREX_SPACEDIM*2)
  amrex::Error("dcomp invalid");

}  // subroutine getBCArray_list

MultiFab* NavierStokes::getState_list(
 int ngrow,Vector<int> scomp,Vector<int> ncomp,
 Real time) {

 if ((ngrow>=0)&&(ngrow<=8)) {
  // do nothing
 } else
  amrex::Error("ngrow out of range in getState_list");

 int ncomp_list=0;
 for (int ilist=0;ilist<scomp.size();ilist++)
  ncomp_list+=ncomp[ilist];

 if (ncomp_list<=0)
  amrex::Error("ncomp_list invalid");
 
 MultiFab* mf = new MultiFab(state[State_Type].boxArray(),dmap,ncomp_list,
   ngrow,MFInfo().SetTag("mf get state list"),FArrayBoxFactory());

 int dcomp=0;
 for (int ilist=0;ilist<scomp.size();ilist++) {
  MultiFab* mf_single=getState(ngrow,scomp[ilist],ncomp[ilist],time);
  MultiFab::Copy(*mf,*mf_single,0,dcomp,ncomp[ilist],ngrow); 
  dcomp+=ncomp[ilist];
  delete mf_single;
 }
 if (dcomp!=ncomp_list)
  amrex::Error("dcomp invalid");

 return mf;

} // subroutine getState_list

// copy localMF[idx_MF] to s_new
void NavierStokes::putState_list(
 Vector<int> scomp,Vector<int> ncomp,
 int idx_MF) {

 int ncomp_list=0;
 for (int ilist=0;ilist<scomp.size();ilist++)
  ncomp_list+=ncomp[ilist];

 if (ncomp_list>=1) {
  // do nothing
 } else
  amrex::Error("ncomp_list invalid");
 
 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

 int scomp_localMF=0;
 for (int ilist=0;ilist<scomp.size();ilist++) {
   // dst,src,scomp,dcomp,ncomp,ngrow
  MultiFab::Copy(S_new,*localMF[idx_MF],scomp_localMF,scomp[ilist],
		 ncomp[ilist],0); 
  scomp_localMF+=ncomp[ilist];
 }
 if (scomp_localMF!=ncomp_list)
  amrex::Error("scomp_localMF invalid");

} // subroutine putState_list


MultiFab* NavierStokes::getState (
  int ngrow, int  scomp,
  int ncomp, Real time) {

 if ((scomp<STATECOMP_ERR)&&(scomp+ncomp-1>=STATECOMP_MOF)) {

  if ((scomp<=STATECOMP_MOF)&&(scomp+ncomp>=STATECOMP_ERR)) {
   // do nothing
  } else {
   std::cout << "VOF: scomp,ncomp " << scomp << ' ' << ncomp << '\n';
   amrex::Error("must get all vof data at once getState");
  }

 }

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 int ntotal=S_new.nComp();
 if (scomp<0)
  amrex::Error("scomp invalid getState"); 
 if (ncomp<=0)
  amrex::Error("ncomp invalid in get state"); 
 if (scomp+ncomp>ntotal)
  amrex::Error("scomp,ncomp invalid");

 MultiFab* mf = new MultiFab(state[State_Type].boxArray(),dmap,ncomp,
   ngrow,MFInfo().SetTag("mf getState"),FArrayBoxFactory());

  //FillPatch is declared in amrlib/AmrLevel.cpp
 FillPatch(*this,*mf,0,time,State_Type,scomp,ncomp,debug_fillpatch);

 ParallelDescriptor::Barrier();

 return mf;
} // end subroutine getState

MultiFab* NavierStokes::getStateSolid (
  int ngrow, int  scomp,
  int ncomp, Real time) {


  // nparts x (velocity + LS + temperature + flag)
 int nparts=im_solid_map.size();
 if ((nparts<1)||(nparts>num_materials))
  amrex::Error("nparts invalid");
 if (ncomp%AMREX_SPACEDIM!=0)
  amrex::Error("ncomp invalid");
 int partid=scomp/AMREX_SPACEDIM;
 if ((partid<0)||(partid>=nparts))
  amrex::Error("partid invalid");

 MultiFab& Solid_new=get_new_data(Solid_State_Type,slab_step+1);
 int ntotal=Solid_new.nComp();
 if (scomp<0)
  amrex::Error("scomp invalid getStateSolid"); 
 if (ncomp<=0)
  amrex::Error("ncomp invalid in getstateSolid"); 
 if (scomp+ncomp>ntotal)
  amrex::Error("scomp,ncomp invalid");

 MultiFab* mf = new MultiFab(state[Solid_State_Type].boxArray(),dmap,ncomp,
   ngrow,MFInfo().SetTag("mf getStateSolid"),FArrayBoxFactory());

 FillPatch(*this,*mf,0,time,Solid_State_Type,scomp,ncomp,debug_fillpatch);

 ParallelDescriptor::Barrier();

 return mf;
} // end subroutine getStateSolid



MultiFab* NavierStokes::getStateTensor (
  int ngrow, int  scomp,
  int ncomp, Real time) {


 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=num_materials)) {

   // 0<=im_elastic_map[i]<num_materials
  int nparts=im_elastic_map.size();

  if (nparts==num_materials_viscoelastic) {

   if (ENUM_NUM_TENSOR_TYPE==2*AMREX_SPACEDIM) {
    // do nothing
   } else
    amrex::Error("ENUM_NUM_TENSOR_TYPE became corrupted");

    // Tensor_Type:
    //   nparts * ENUM_NUM_TENSOR_TYPE
   int ntotal_test=NUM_CELL_ELASTIC;

   if (NUM_CELL_ELASTIC==num_materials_viscoelastic*ENUM_NUM_TENSOR_TYPE) {
    // do nothing
   } else
    amrex::Error("NUM_CELL_ELASTIC invalid");

   if ((ncomp==ntotal_test)&&(scomp==0)) {
    // do nothing
   } else if (ncomp%ENUM_NUM_TENSOR_TYPE==0) {
    int partid=scomp/ENUM_NUM_TENSOR_TYPE;
    if ((partid<0)||(partid>=nparts))
     amrex::Error("partid invalid");
   } else {
    std::cout << "ncomp= " << ncomp << 
      " scomp=" << scomp << 
      " num_materials_viscoelastic= " << num_materials_viscoelastic << 
      " ENUM_NUM_TENSOR_TYPE= " << ENUM_NUM_TENSOR_TYPE << 
      " NUM_CELL_ELASTIC= " << NUM_CELL_ELASTIC << '\n';
    amrex::Error("ncomp or scomp invalid");
   }

   MultiFab& Tensor_new=get_new_data(Tensor_Type,slab_step+1);
   int ntotal=Tensor_new.nComp();
   if (ntotal==ntotal_test) {
    // do nothing
   } else
    amrex::Error("ntotal invalid");

   if (scomp<0)
    amrex::Error("scomp invalid getStateTensor"); 
   if (ncomp<=0)
    amrex::Error("ncomp invalid in getstateTensor"); 
   if (scomp+ncomp>ntotal)
    amrex::Error("scomp,ncomp invalid");

   MultiFab* mf = new MultiFab(state[Tensor_Type].boxArray(),dmap,ncomp,
    ngrow,MFInfo().SetTag("mf getStateTensor"),FArrayBoxFactory());

   FillPatch(*this,*mf,0,time,Tensor_Type,scomp,ncomp, debug_fillpatch);

   ParallelDescriptor::Barrier();

   return mf;
  } else
   amrex::Error("nparts!=num_materials_viscoelastic");
 } else
  amrex::Error("num_materials_viscoelastic bad:getStateTensor");

 return nullptr;

} // end subroutine getStateTensor

MultiFab* NavierStokes::getStateDist (int ngrow,Real time,
   const std::string& caller_string) {

 if (verbose>0)
  if (ParallelDescriptor::IOProcessor())
   std::cout << "getStateDist: time,caller_string " << time << ' ' << 
     caller_string <<'\n';

 if ((ngrow<0)||(ngrow>ngrow_distance))
  amrex::Error("ngrow invalid");

 
 MultiFab& S_new=get_new_data(LS_Type,slab_step+1);
 int ntotal=S_new.nComp();
 if (ntotal!=num_materials*(AMREX_SPACEDIM+1))
  amrex::Error("ntotal invalid");

 MultiFab* mf = new MultiFab(state[State_Type].boxArray(),dmap,
   num_materials*(AMREX_SPACEDIM+1),
   ngrow,MFInfo().SetTag("mf getStateDist"),FArrayBoxFactory());

  // scomp=0
 FillPatch(*this,*mf,0,time,LS_Type,0,num_materials*(AMREX_SPACEDIM+1),debug_fillpatch);

 ParallelDescriptor::Barrier();

 return mf;
} // subroutine getStateDist

void NavierStokes::cpp_overridepbc(int homflag_in,int project_option_in) {

 if ((homflag_in==0)||(homflag_in==1)) {
  override_bc_to_homogeneous=homflag_in;
   //fort_overridepbc is declared in: PROB.F90
  fort_overridepbc(&homflag_in,&project_option_in);
 } else
  amrex::Error("homflag_in invalid");

}  // subroutine cpp_overridepbc

MultiFab* NavierStokes::getStateDIV_DATA(int ngrow,
		int scomp,int ncomp,Real time) {

 int save_bc_status=override_bc_to_homogeneous;
  //cpp_overridepbc is declared in: NavierStokes.cpp
  //homflag_in=1
 cpp_overridepbc(1,SOLVETYPE_PRES);

 MultiFab& S_new=get_new_data(DIV_Type,slab_step+1);
 int ntotal=S_new.nComp();
 if (ntotal!=1)
  amrex::Error("ntotal invalid");

 MultiFab* mf = new MultiFab(state[DIV_Type].boxArray(),dmap,ncomp,
   ngrow,MFInfo().SetTag("mf getStateDIV_DATA"),FArrayBoxFactory());

 FillPatch(*this,*mf,0,time,DIV_Type,scomp,ncomp,debug_fillpatch);

 ParallelDescriptor::Barrier();

 cpp_overridepbc(save_bc_status,SOLVETYPE_PRES);

 return mf;
} // subroutine getStateDIV_DATA


void NavierStokes::putStateDIV_DATA(
 int scomp,int ncomp,int idx_MF) {

 MultiFab& S_new=get_new_data(DIV_Type,slab_step+1);
 int ntotal=S_new.nComp();
 if (ntotal!=1)
  amrex::Error("ntotal invalid");

 int scomp_localMF=0;
   // dst,src,scomp,dcomp,ncomp,ngrow
 MultiFab::Copy(S_new,*localMF[idx_MF],scomp_localMF,scomp,
  	        ncomp,0); 

} // subroutine putStateDIV_DATA



// called from:
// NavierStokes::prepare_post_process if via "post_init_state"
// NavierStokes::do_the_advance
void
NavierStokes::makeStateDistALL(int keep_all_interfaces,
   int ngrow_make_distance_accept) {

 if (level!=0)
  amrex::Error("level invalid in makeStateDistALL");

 int finest_level=parent->finestLevel();

 Real problo[AMREX_SPACEDIM];
 Real probhi[AMREX_SPACEDIM];
 Real problen[AMREX_SPACEDIM];
 max_problen=0.0;
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  problo[dir]=geom.ProbLo(dir);
  probhi[dir]=geom.ProbHi(dir);
  problen[dir]=probhi[dir]-problo[dir];
  if (problen[dir]>0.0) {
   max_problen+=(problen[dir]*problen[dir]);
  } else
   amrex::Error("problen[dir] invalid");
 }
 max_problen=std::sqrt(max_problen);
 if (max_problen>0.0) {
  // do nothing
 } else
  amrex::Error("max_problen invalid");

 minLS.resize(thread_class::nthreads);
 maxLS.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  minLS[tid].resize(num_materials);
  maxLS[tid].resize(num_materials);
  for (int im=0;im<num_materials;im++) {
   minLS[tid][im]=max_problen;
   maxLS[tid][im]=-max_problen;
  } // tid
 }

  // traverse from coarsest to finest so that
  // coarse data normals will be available for filling in
  // ghost values.
  // Also do a fillCoarsePatch for the level set functions and DIST_TOUCH_MF
  // (using piecewise constant interp) in order to init. some level set
  // function values.
 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.makeStateDist(keep_all_interfaces,ngrow_make_distance_accept);
 }
  // fort_correct_uninit is in MOF_REDIST_3D.F90
 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.correct_dist_uninit();
 }

 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.delete_localMF(FACEFRAC_MF,1);
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   ns_level.delete_localMF(FACEFRAC_SOLVE_MM_MF+dir,1);
  }
  ns_level.delete_localMF(FACETEST_MF,1);
  ns_level.delete_localMF(ORIGDIST_MF,1);
  ns_level.delete_localMF(STENCIL_MF,1);
  ns_level.delete_localMF(DIST_TOUCH_MF,1);
 }


 if (verbose>0)
  if (ParallelDescriptor::IOProcessor())
   for (int im=0;im<num_materials;im++) {
    std::cout << "im= " << im << '\n';
    std::cout << "minLS = " << minLS[0][im] << '\n';
    std::cout << "maxLS = " << maxLS[0][im] << '\n';
   }

} // subroutine makeStateDistALL()

// called from: NavierStokes::do_the_advance 
// (prior to level_phase_change_rate) and
// called from: NavierStokes::init_FSI_GHOST_MAC_MF
// fd_mf used by the GNBC algorithm.
void 
NavierStokes::build_NRM_FD_MF(int fd_mf,int ls_mf) {

 bool use_tiling=ns_tiling;
 int finest_level=parent->finestLevel();
 const Real* dx = geom.CellSize();

 if (ngrow_distance==4) {
  // do nothing
 } else
  amrex::Error("ngrow_distance invalid");

 if ((localMF[fd_mf]->nGrow()>=ngrow_distance)&&
     (localMF[ls_mf]->nGrow()>=ngrow_distance)) {
  if ((localMF[fd_mf]->nComp()==num_materials*AMREX_SPACEDIM)&&
      (localMF[ls_mf]->nComp()==num_materials*(AMREX_SPACEDIM+1))) {
   MultiFab::Copy(*localMF[fd_mf],*localMF[ls_mf],num_materials,0,
	    num_materials*AMREX_SPACEDIM,ngrow_distance);

   if (thread_class::nthreads<1)
    amrex::Error("thread_class::nthreads invalid");
   thread_class::init_d_numPts(localMF[fd_mf]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
   for (MFIter mfi(*localMF[fd_mf],use_tiling); mfi.isValid(); ++mfi) {
     BL_ASSERT(grids[mfi.index()] == mfi.validbox());
     const int gridno = mfi.index();
     const Box& tilegrid = mfi.tilebox();
     const Box& fabgrid = grids[gridno];
     const int* tilelo=tilegrid.loVect();
     const int* tilehi=tilegrid.hiVect();
     const int* fablo=fabgrid.loVect();
     const int* fabhi=fabgrid.hiVect();
     int bfact=parent->Space_blockingFactor(level);

     const Real* xlo = grid_loc[gridno].lo();

     FArrayBox& lsfab=(*localMF[ls_mf])[mfi];
     FArrayBox& nrmfab=(*localMF[fd_mf])[mfi];
  
     int tid_current=ns_thread();
     if ((tid_current<0)||(tid_current>=thread_class::nthreads))
      amrex::Error("tid_current invalid");
     thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

      // fort_dt_normal is declared in: MOF_REDIST_3D.F90
     fort_fd_normal( 
      &level,
      &finest_level,
      lsfab.dataPtr(),
      ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
      nrmfab.dataPtr(),
      ARLIM(nrmfab.loVect()),ARLIM(nrmfab.hiVect()),
      tilelo,tilehi,
      fablo,fabhi,
      &bfact,
      xlo,dx);
   } // mfi
} // omp
   ns_reconcile_d_num(LOOP_FD_NORMAL,"build_NRM_FD_MF");

   localMF[fd_mf]->FillBoundary(geom.periodicity()); 
  } else
    amrex::Error("fd_mf or ls_mf nComp() invalid");
 } else
  amrex::Error("fd_mf or ls_mf nGrow() invalid");
			 
} // end subroutine build_NRM_FD_MF

void
NavierStokes::makeStateDist(int keep_all_interfaces,
		int ngrow_make_distance_accept) {

 std::string local_caller_string="makeStateDist";

 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();

 int profile_dist=0;

 int dist_timings=0;

 if (verbose>0) {
  dist_timings=1;
 } else if (verbose==0) {
  // do nothing
 } else
  amrex::Error("verbose invalid");

 double before_dist;
 double after_dist;
 double before_profile;
 double after_profile;

 if (dist_timings==1)
  before_dist = ParallelDescriptor::second();

 const Real* dx = geom.CellSize();

 MultiFab& LS_new = get_new_data(LS_Type,slab_step+1);

 getStateDist_localMF(ORIGDIST_MF,ngrow_distance,cur_time_slab,
   local_caller_string);

 resize_mask_nbr(ngrow_distance);
 debug_ngrow(MASK_NBR_MF,ngrow_distance,local_caller_string);
 if (localMF[MASK_NBR_MF]->nComp()!=4)
  amrex::Error("invalid ncomp for mask nbr");
 VOF_Recon_resize(ngrow_distance); //output:SLOPE_RECON_MF
 debug_ngrow(SLOPE_RECON_MF,ngrow_distance,local_caller_string);
 debug_ngrow(ORIGDIST_MF,ngrow_distance,local_caller_string);
 if (localMF[ORIGDIST_MF]->nComp()!=num_materials*(1+AMREX_SPACEDIM))
  amrex::Error("invalid ncomp for origdist");

 if (profile_dist==1)
  before_profile = ParallelDescriptor::second();

 new_localMF(DIST_TOUCH_MF,num_materials,0,-1);
 localMF[DIST_TOUCH_MF]->setVal(0.0);

 MultiFab* dist_coarse_mf;
 MultiFab* dist_touch_coarse_mf;

 if ((level>0)&&(level<=finest_level)) {
   //ngrow=0
  dist_coarse_mf=new MultiFab(grids,dmap,num_materials*(AMREX_SPACEDIM+1),0,
	MFInfo().SetTag("dist_coarse_mf"),FArrayBoxFactory());
  int dcomp=0;
  int scomp=0;
  FillCoarsePatch(*dist_coarse_mf,dcomp,cur_time_slab,LS_Type,scomp,
        num_materials*(AMREX_SPACEDIM+1),debug_fillpatch);

  // idx,scomp,ncomp,index,scompBC_map
  // FillCoarsePatchGHOST is ultimately called.
  // dest_lstGHOST for State_Type defaults to pc_interp.
  // scompBC_map==0 corresponds to pc_interp and fort_extrapfill
  dist_touch_coarse_mf=
   new MultiFab(grids,dmap,num_materials,0,
    MFInfo().SetTag("dist_touch_coarse_mf"),FArrayBoxFactory());

  for (int i=0;i<num_materials;i++) {
   Vector<int> scompBC_map;
   scompBC_map.resize(1);
   scompBC_map[0]=0;
    //scomp=i ncomp=1
   PCINTERP_fill_coarse_patch(DIST_TOUCH_MF,i,
     1,State_Type,scompBC_map);
  } // i=0..num_materials-1
  MultiFab::Copy(*dist_touch_coarse_mf,*localMF[DIST_TOUCH_MF],0,0,num_materials,0);
  localMF[DIST_TOUCH_MF]->setVal(0.0);

 } else if (level==0) {

  dist_coarse_mf=localMF[ORIGDIST_MF];
  dist_touch_coarse_mf=localMF[DIST_TOUCH_MF];

 } else
  amrex::Error("level invalid 6"); 

 if (profile_dist==1) {
  after_profile = ParallelDescriptor::second();
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "level= " << level << '\n';
   std::cout << "coarse dist time " << after_profile-before_profile << '\n';
  }
 }

 if (profile_dist==1)
  before_profile = ParallelDescriptor::second();

 int nstar=9;
 if (AMREX_SPACEDIM==3)
  nstar*=3;
 int nface=num_materials*AMREX_SPACEDIM*2; 

  // (num_materials,num_materials,2)  left material, right material, 
  // frac_pair+dist_pair
 int nface_dst=num_materials*num_materials*2;

 new_localMF(STENCIL_MF,nstar,ngrow_distance,-1);
 localMF[STENCIL_MF]->setVal(0.0);

 int tessellate=0;
  // fort_faceinit is in: MOF_REDIST_3D.F90
 makeFaceFrac(tessellate,ngrow_distance,FACEFRAC_MF);
  // fort_faceprocess is in: MOF_REDIST_3D.F90
 ProcessFaceFrac(tessellate,FACEFRAC_MF,FACEFRAC_SOLVE_MM_MF,ngrow_distance);

 if (profile_dist==1) {
  after_profile = ParallelDescriptor::second();
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "level= " << level << '\n';
   std::cout << "makeFaceFrac time " << after_profile-before_profile << '\n';
  }
 }

 if (profile_dist==1)
  before_profile = ParallelDescriptor::second();

  // fort_faceinittest is in MOF_REDIST_3D.F90
  // FACETEST_MF has num_materials * sdim components.
  // tessellate=0
 makeFaceTest(tessellate,ngrow_distance,FACETEST_MF);

 if (profile_dist==1) {
  after_profile = ParallelDescriptor::second();
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "level= " << level << '\n';
   std::cout << "makeFaceTest time " << after_profile-before_profile << '\n';
  }
 }

 if (profile_dist==1)
  before_profile = ParallelDescriptor::second();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(LS_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(LS_new,use_tiling); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();
   int bfact=parent->Space_blockingFactor(level);

   const Real* xlo = grid_loc[gridno].lo();

   FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi];
   FArrayBox& voffab=(*localMF[SLOPE_RECON_MF])[mfi];
   FArrayBox& stencilfab=(*localMF[STENCIL_MF])[mfi];

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // in: MOF_REDIST_3D.F90
    // internally sets tessellate=0
    // 3x3 stencil for each cell in 2D
    // 3x3x3 stencil for each cell in 3D
    // fluid material id for each cell edge point is initialized.
   fort_steninit( 
    &level,
    &finest_level,
    stencilfab.dataPtr(),
    ARLIM(stencilfab.loVect()),ARLIM(stencilfab.hiVect()),
    maskfab.dataPtr(),
    ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
    voffab.dataPtr(),
    ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    &NS_geometry_coord,
    xlo,dx,
    &cur_time_slab,
    &ngrow_distance,
    &nstar);
 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_STENINIT,"makeStateDist");

 localMF[STENCIL_MF]->FillBoundary(geom.periodicity());

 if (profile_dist==1) {
  after_profile = ParallelDescriptor::second();
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "level= " << level << '\n';
   std::cout << "STENINIT time " << after_profile-before_profile << '\n';
  }
 }

 debug_ngrow(FACETEST_MF,ngrow_distance,local_caller_string);
 if (localMF[FACETEST_MF]->nComp()!=num_materials*AMREX_SPACEDIM)
  amrex::Error("localMF[FACETEST_MF]->nComp() invalid");

 debug_ngrow(FACEFRAC_MF,ngrow_distance,local_caller_string);
 if (localMF[FACEFRAC_MF]->nComp()!=nface)
  amrex::Error("localMF[FACEFRAC_MF]->nComp() invalid");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  debug_ngrow(FACEFRAC_SOLVE_MM_MF+dir,ngrow_distance,local_caller_string);
  if (localMF[FACEFRAC_SOLVE_MM_MF+dir]->nComp()!=nface_dst)
   amrex::Error("localMF[FACEFRAC_SOLVE_MM_MF+dir]->nComp()!=nface_dst");
 }

 debug_ngrow(STENCIL_MF,ngrow_distance,local_caller_string);

 Vector<int> nprocessed;
 nprocessed.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  nprocessed[tid]=0.0;
 }

 Real profile_time_start=0.0;
 if ((profile_debug==1)||(profile_dist==1)) {
  profile_time_start=ParallelDescriptor::second();
 }

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(LS_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(LS_new,use_tiling); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();
   int bfact=parent->Space_blockingFactor(level);

   const Real* xlo = grid_loc[gridno].lo();

   FArrayBox& voffab=(*localMF[SLOPE_RECON_MF])[mfi];
   FArrayBox& lsfab=LS_new[mfi];
   FArrayBox& origdist=(*localMF[ORIGDIST_MF])[mfi];

   FArrayBox& stencilfab=(*localMF[STENCIL_MF])[mfi];

   FArrayBox& facepairXfab=(*localMF[FACEFRAC_SOLVE_MM_MF])[mfi];
   FArrayBox& facepairYfab=(*localMF[FACEFRAC_SOLVE_MM_MF+1])[mfi];
   FArrayBox& facepairZfab=
      (*localMF[FACEFRAC_SOLVE_MM_MF+AMREX_SPACEDIM-1])[mfi];

   FArrayBox& facetestfab=(*localMF[FACETEST_MF])[mfi];
   FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi];

   FArrayBox& touchfab=(*localMF[DIST_TOUCH_MF])[mfi];
   FArrayBox& crsetouch=(*dist_touch_coarse_mf)[mfi];
   FArrayBox& crsedist=(*dist_coarse_mf)[mfi];

   Vector<int> vofbc=getBCArray(State_Type,gridno,STATECOMP_MOF,1);

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // in: MOF_REDIST_3D.F90
   fort_levelstrip( 
    &keep_all_interfaces,
    &nprocessed[tid_current],
    minLS[tid_current].dataPtr(),
    maxLS[tid_current].dataPtr(),
    &max_problen,
    &level,
    &finest_level,
    truncate_volume_fractions.dataPtr(),
    maskfab.dataPtr(),
    ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
    facepairXfab.dataPtr(),
    ARLIM(facepairXfab.loVect()),ARLIM(facepairXfab.hiVect()),
    facepairYfab.dataPtr(),
    ARLIM(facepairYfab.loVect()),ARLIM(facepairYfab.hiVect()),
    facepairZfab.dataPtr(),
    ARLIM(facepairZfab.loVect()),ARLIM(facepairZfab.hiVect()),
    facetestfab.dataPtr(),
    ARLIM(facetestfab.loVect()),ARLIM(facetestfab.hiVect()),
    stencilfab.dataPtr(),
    ARLIM(stencilfab.loVect()),ARLIM(stencilfab.hiVect()),
    voffab.dataPtr(),
    ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),
    origdist.dataPtr(),
    ARLIM(origdist.loVect()),ARLIM(origdist.hiVect()),
    lsfab.dataPtr(),
    ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
    touchfab.dataPtr(),
    ARLIM(touchfab.loVect()),ARLIM(touchfab.hiVect()),
    crsetouch.dataPtr(),
    ARLIM(crsetouch.loVect()),ARLIM(crsetouch.hiVect()),
    crsedist.dataPtr(),
    ARLIM(crsedist.loVect()),ARLIM(crsedist.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    vofbc.dataPtr(),
    &NS_geometry_coord,
    xlo,dx,
    &cur_time_slab,
    &ngrow_distance,
    &ngrow_make_distance_accept,
    &nstar,
    &nface_dst);
 } // mfi
} // omp

 ns_reconcile_d_num(LOOP_LEVELSTRIP,"makeStateDist");

 for (int tid=1;tid<thread_class::nthreads;tid++) {
  for (int im=0;im<num_materials;im++) {
   if (minLS[tid][im]<minLS[0][im])
    minLS[0][im]=minLS[tid][im];
   if (maxLS[tid][im]>maxLS[0][im])
    maxLS[0][im]=maxLS[tid][im];
  } // im=0..num_materials-1
  nprocessed[0]+=nprocessed[tid];
 }
 ParallelDescriptor::ReduceIntSum(nprocessed[0]);
 for (int im=0;im<num_materials;im++) {
  ParallelDescriptor::ReduceRealMax(maxLS[0][im]);
  ParallelDescriptor::ReduceRealMin(minLS[0][im]);
 }

 if ((level>0)&&(level<=finest_level)) {
  delete dist_coarse_mf;
  delete dist_touch_coarse_mf;
 } else if (level==0) {
  // do nothing
 } else
  amrex::Error("level invalid 7");
 
 if ((profile_debug==1)||(profile_dist==1)) {
  double profile_time_end=ParallelDescriptor::second();
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "level= " << level << '\n';
   std::cout << "nprocessed= " << nprocessed[0] << '\n';
   std::cout << "profile LEVELSTRIP time = " << 
    profile_time_end-profile_time_start << '\n';
  }
 }

 if (dist_timings==1) {
  after_dist = ParallelDescriptor::second();
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "level= " << level << '\n';
   std::cout << "dist cpu time " << after_dist-before_dist << '\n';
  }
 }

} // subroutine makeStateDist


void
NavierStokes::correct_dist_uninit() {

 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();

 const Real* dx = geom.CellSize();

 MultiFab& LS_new = get_new_data(LS_Type,slab_step+1);
 if (localMF[DIST_TOUCH_MF]->nComp()!=num_materials)
  amrex::Error("localMF[DIST_TOUCH_MF]->nComp()!=num_materials");

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(LS_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(LS_new,use_tiling); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();
   int bfact=parent->Space_blockingFactor(level);

   const Real* xlo = grid_loc[gridno].lo();

   FArrayBox& lsfab=LS_new[mfi];
   FArrayBox& touchfab=(*localMF[DIST_TOUCH_MF])[mfi];

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   fort_correct_uninit( 
    minLS[0].dataPtr(),
    maxLS[0].dataPtr(),
    &max_problen,
    &level,
    &finest_level,
    lsfab.dataPtr(),
    ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
    touchfab.dataPtr(),
    ARLIM(touchfab.loVect()),ARLIM(touchfab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    xlo,dx,
    &cur_time_slab);
 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_CORRECT_UNINIT,"correct_dist_uninit");

} // subroutine correct_dist_uninit


// WARNING:  allocates, but does not delete.
// called from NavierStokes::ColorSum
// called from NavierStokes::makeStateDist
void
NavierStokes::ProcessFaceFrac(int tessellate,int idxsrc,int idxdst,
		int ngrow_dest) {
  
 std::string local_caller_string="ProcessFaceFrac";

 bool use_tiling=ns_tiling;

  // (num_materials,sdim,2,sdim+1) area on each face of a cell.
 int nface_src=num_materials*AMREX_SPACEDIM*2; 
  // (num_materials,num_materials,2)  left material, right material, frac_pair+dist_pair
 int nface_dst=num_materials*num_materials*2;

 int finest_level=parent->finestLevel();

 if ((tessellate!=0)&&(tessellate!=1)&&(tessellate!=3))
  amrex::Error("tessellate invalid60");

 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid ProcessFaceFrac");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  delete_localMF_if_exist(idxdst+dir,1);
  new_localMF(idxdst+dir,nface_dst,ngrow_dest,dir);
  localMF[idxdst+dir]->setVal(0.0);
 }

 int ngrow_source=ngrow_dest;
 if (ngrow_dest==0)
  ngrow_source=1;

 debug_ngrow(idxsrc,ngrow_source,local_caller_string);
 if (localMF[idxsrc]->nComp()!=nface_src)
  amrex::Error("idxsrc has invalid ncomp");

 VOF_Recon_resize(ngrow_source); //output:SLOPE_RECON_MF
 debug_ngrow(SLOPE_RECON_MF,ngrow_source,local_caller_string);

   // (1) =1 interior  =1 fine-fine ghost in domain  =0 otherwise
   // (2) =1 interior  =0 otherwise
   // (3) =1 interior+ngrow-1  =0 otherwise
   // (4) =1 interior+ngrow    =0 otherwise
 resize_mask_nbr(ngrow_source);
 debug_ngrow(MASK_NBR_MF,ngrow_source,local_caller_string);

 const Real* dx = geom.CellSize();

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(localMF[idxsrc]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(*localMF[idxsrc],use_tiling); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();
   int bfact=parent->Space_blockingFactor(level);

   const Real* xlo = grid_loc[gridno].lo();

   FArrayBox& voffab=(*localMF[SLOPE_RECON_MF])[mfi];
   FArrayBox& facefab=(*localMF[idxsrc])[mfi];
   FArrayBox& dstfab=(*localMF[idxdst+dir])[mfi];

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // in: MOF_REDIST_3D.F90
   fort_faceprocess( 
    &ngrow_source,
    &ngrow_dest,
    &tid_current,
    &dir,
    &tessellate, // =0,1, or 3
    &level,
    &finest_level,
    dstfab.dataPtr(),
    ARLIM(dstfab.loVect()),ARLIM(dstfab.hiVect()),
    facefab.dataPtr(),
    ARLIM(facefab.loVect()),ARLIM(facefab.hiVect()),
    voffab.dataPtr(),
    ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    &NS_geometry_coord,
    xlo,dx,
    &cur_time_slab,
    &nface_src,&nface_dst);
  } // mfi
} // omp
  ns_reconcile_d_num(LOOP_FACEPROCESS,"ProcessFaceFrac");

  localMF[idxdst+dir]->FillBoundary(geom.periodicity());
 } //dir=0..sdim-1

} // subroutine ProcessFaceFrac



// WARNING: makeFaceFrac allocates, but does not delete.
// if caller from: makeStateDist (tessellate==0)
// if caller from: ColorSum (tessellate==1)
void
NavierStokes::makeFaceFrac(
 int tessellate,int ngrow,int idx) {

 std::string local_caller_string="makeFaceFrac";
 
 bool use_tiling=ns_tiling;

  // (num_materials,sdim,2) area on each face of a cell.
 int nface=num_materials*AMREX_SPACEDIM*2; 

 int finest_level=parent->finestLevel();

 delete_localMF_if_exist(idx,1);
 new_localMF(idx,nface,ngrow,-1);
 localMF[idx]->setVal(0.0);
 VOF_Recon_resize(ngrow); //output:SLOPE_RECON_MF
 debug_ngrow(SLOPE_RECON_MF,ngrow,local_caller_string);
 resize_mask_nbr(ngrow);
 debug_ngrow(MASK_NBR_MF,ngrow,local_caller_string);

 const Real* dx = geom.CellSize();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(localMF[idx]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*localMF[idx],use_tiling); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();
   int bfact=parent->Space_blockingFactor(level);

   const Real* xlo = grid_loc[gridno].lo();

   FArrayBox& voffab=(*localMF[SLOPE_RECON_MF])[mfi];
   FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi];
   FArrayBox& facefab=(*localMF[idx])[mfi];

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // in: MOF_REDIST_3D.F90
   fort_faceinit( 
    &tid_current,
    &tessellate,  // 0,1, or 3
    &level,
    &finest_level,
    facefab.dataPtr(),
    ARLIM(facefab.loVect()),ARLIM(facefab.hiVect()),
    maskfab.dataPtr(),
    ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
    voffab.dataPtr(),
    ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    &NS_geometry_coord,
    xlo,dx,
    &cur_time_slab,
    &ngrow,
    &nface);
 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_FACEINIT,"makeFaceFrac");

 localMF[idx]->FillBoundary(geom.periodicity());

} // subroutine makeFaceFrac


// WARNING: makeFaceTest allocates, but does not delete.
void
NavierStokes::makeFaceTest(int tessellate,int ngrow,int idx) {

 std::string local_caller_string="makeFaceTest";
 
 bool use_tiling=ns_tiling;

  // (im,dir,side)  area  
 int nface=num_materials*AMREX_SPACEDIM*2; 

 int finest_level=parent->finestLevel();

 if ((tessellate!=0)&&(tessellate!=1)&&(tessellate!=3))
  amrex::Error("tessellate invalid61");

 if (localMF_grow[idx]==-1) {
  // do nothing
 } else
  amrex::Error("makeFaceTest: forgot to delete");

 new_localMF(idx,num_materials*AMREX_SPACEDIM,ngrow,-1);
 localMF[idx]->setVal(0.0);

 VOF_Recon_resize(ngrow); //output:SLOPE_RECON_MF
 debug_ngrow(SLOPE_RECON_MF,ngrow,local_caller_string);
 debug_ngrow(FACEFRAC_MF,ngrow,local_caller_string);
 resize_mask_nbr(ngrow);
 debug_ngrow(MASK_NBR_MF,ngrow,local_caller_string);

 if (localMF[FACEFRAC_MF]->nComp()!=nface)
  amrex::Error("localMF[FACEFRAC_MF]->nComp() invalid");

 const Real* dx = geom.CellSize();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(localMF[idx]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*localMF[idx],use_tiling); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();
   int bfact=parent->Space_blockingFactor(level);

   const Real* xlo = grid_loc[gridno].lo();

   FArrayBox& voffab=(*localMF[SLOPE_RECON_MF])[mfi];
   FArrayBox& facefab=(*localMF[FACEFRAC_MF])[mfi];
   FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi];
   FArrayBox& facetest=(*localMF[idx])[mfi];

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   fort_faceinittest( 
    &tid_current,
    &tessellate,
    &level,
    &finest_level,
    facefab.dataPtr(),
    ARLIM(facefab.loVect()),ARLIM(facefab.hiVect()),
    facetest.dataPtr(),
    ARLIM(facetest.loVect()),ARLIM(facetest.hiVect()),
    maskfab.dataPtr(),
    ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
    voffab.dataPtr(),
    ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    &NS_geometry_coord,
    xlo,dx,
    &cur_time_slab,
    &ngrow,
    &nface);
 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_FACEINITTEST,"makeFaceTest");

} // subroutine makeFaceTest


// WARNING: makeCellFrac allocates, but does not delete.
void
NavierStokes::makeCellFrac(
  int tessellate, // = 0,1, or 3
  int ngrow,int idx) {
 
 std::string local_caller_string="makeCellFrac";

 bool use_tiling=ns_tiling;

  // (num_materials,num_materials,3+sdim)
  // im_inside,im_outside,3+sdim --> area, dist_to_line, dist, line normal.
 int ncellfrac=num_materials*num_materials*(3+AMREX_SPACEDIM); 
 int finest_level=parent->finestLevel();

 delete_localMF_if_exist(idx,1);
 new_localMF(idx,ncellfrac,ngrow,-1);
 localMF[idx]->setVal(0.0);

 int ngrow_resize=ngrow;
 if (ngrow_resize==0)
  ngrow_resize++;

 VOF_Recon_resize(ngrow_resize); //output:SLOPE_RECON_MF
 debug_ngrow(SLOPE_RECON_MF,ngrow_resize,local_caller_string);
 resize_mask_nbr(ngrow_resize);
 debug_ngrow(MASK_NBR_MF,ngrow_resize,local_caller_string);

 const Real* dx = geom.CellSize();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(localMF[idx]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*localMF[idx],use_tiling); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();
   int bfact=parent->Space_blockingFactor(level);

   const Real* xlo = grid_loc[gridno].lo();

   FArrayBox& voffab=(*localMF[SLOPE_RECON_MF])[mfi];
   FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi];
   FArrayBox& facefab=(*localMF[idx])[mfi];

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // in: LEVELSET_3D.F90
   fort_cellfaceinit( 
    &tid_current,
    &tessellate,  // = 0,1, or 3
    &level,
    &finest_level,
    facefab.dataPtr(),
    ARLIM(facefab.loVect()),ARLIM(facefab.hiVect()),
    maskfab.dataPtr(),
    ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
    voffab.dataPtr(),
    ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    &NS_geometry_coord,
    xlo,dx,
    &cur_time_slab,
    &ngrow,
    &ncellfrac);
 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_CELLFACEINIT,"makeCellFrac");

 localMF[idx]->FillBoundary(geom.periodicity());

} // makeCellFrac


// called from make_physics_varsALL
void
NavierStokes::makeStateCurv(int project_option,
   const std::string& caller_string) {
 
 std::string local_caller_string="makeStateCurv";
 local_caller_string=caller_string+local_caller_string;

 if (pattern_test(local_caller_string,"volWgtSumALL")==1) {
  //do nothing
 } else if (pattern_test(local_caller_string,"prepare_post_process")==1) {
  //do nothing
 } else if (pattern_test(local_caller_string,"do_the_advance")==1) {
  //do nothing
 } else {
  std::cout << "local_caller_string=" << local_caller_string << '\n';
  amrex::Error("local_caller_string invalid");
 }

 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid makeStateCurv");

 if (ngrow_distance!=4)
  amrex::Error("ngrow_distance invalid");

 Vector< Real > curv_min_local;
 Vector< Real > curv_max_local;
 curv_min_local.resize(thread_class::nthreads);
 curv_max_local.resize(thread_class::nthreads);

 for (int tid=0;tid<thread_class::nthreads;tid++) {
  curv_min_local[tid]=1.0e+99;
  curv_max_local[tid]=-1.0e+99;
 } // tid

 int num_curv=num_interfaces*CURVCOMP_NCOMP;

 resize_metrics(1);

 resize_levelset(ngrow_distance,LEVELPC_MF);

 if (localMF[LEVELPC_MF]->nComp()!=num_materials*(1+AMREX_SPACEDIM))
  amrex::Error("localMF[LEVELPC_MF]->nComp() invalid");
 if (localMF[LEVELPC_MF]->nGrow()!=ngrow_distance)
  amrex::Error("localMF[LEVELPC_MF]->nGrow() invalid");

 int nhistory=localMF[HISTORY_MF]->nComp();
 if (nhistory==num_interfaces*2) {
  // do nothing
 } else
  amrex::Error("nhistory invalid");

 delete_localMF_if_exist(DIST_CURV_MF,1);
 new_localMF(DIST_CURV_MF,num_curv,1,-1);
 localMF[DIST_CURV_MF]->setVal(0.0);

 VOF_Recon_resize(ngrow_distance); //output:SLOPE_RECON_MF
 debug_ngrow(SLOPE_RECON_MF,ngrow_distance,local_caller_string);
 if (localMF[SLOPE_RECON_MF]->nComp()==num_materials*ngeom_recon) {
  // do nothing
 } else
  amrex::Error("localMF[SLOPE_RECON_MF]->nComp() invalid");

 if ((project_option==SOLVETYPE_PRES)||
     (project_option==SOLVETYPE_INITPROJ)) {

  const Real* dx = geom.CellSize();

  Real cl_time=prev_time_slab;

  if (project_option==SOLVETYPE_PRES)  
   cl_time=prev_time_slab;
  else if (project_option==SOLVETYPE_INITPROJ) 
   cl_time=cur_time_slab;
  else
   amrex::Error("project_option invalid makeStateCurv");

  MultiFab* CL_velocity=getState(2,STATECOMP_VEL,
     STATE_NCOMP_VEL+STATE_NCOMP_PRES,cl_time);
  MultiFab* den=getStateDen(2,cl_time);
  if (den->nComp()!=num_materials*num_state_material)
   amrex::Error("invalid ncomp for den");

   // mask=1 if not covered or if outside the domain.
   // NavierStokes::maskfiner_localMF
   // NavierStokes::maskfiner
  resize_maskfiner(1,MASKCOEF_MF);
  debug_ngrow(MASKCOEF_MF,1,local_caller_string); 
  resize_mask_nbr(1);
  debug_ngrow(MASK_NBR_MF,1,local_caller_string);

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(CL_velocity->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(*CL_velocity,use_tiling); mfi.isValid(); ++mfi) {
    BL_ASSERT(grids[mfi.index()] == mfi.validbox());
    const int gridno = mfi.index();
    const Box& tilegrid = mfi.tilebox();
    const Box& fabgrid = grids[gridno];
    const int* tilelo=tilegrid.loVect();
    const int* tilehi=tilegrid.hiVect();
    const int* fablo=fabgrid.loVect();
    const int* fabhi=fabgrid.hiVect();
    int bfact_space=parent->Space_blockingFactor(level);
    int bfact_grid=parent->Old_blockingFactor(level);

    const Real* xlo = grid_loc[gridno].lo();

    FArrayBox& histfab=(*localMF[HISTORY_MF])[mfi];

    // mask=tag if not covered by level+1 or outside the domain.
    FArrayBox& maskcov=(*localMF[MASKCOEF_MF])[mfi];

    FArrayBox& lsfab=(*localMF[LEVELPC_MF])[mfi];
    FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi];

    FArrayBox& curvfab=(*localMF[DIST_CURV_MF])[mfi];
    FArrayBox& velfab=(*CL_velocity)[mfi];
    FArrayBox& denfab=(*den)[mfi];
    FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi];

    FArrayBox& areax=(*localMF[AREA_MF])[mfi];
    FArrayBox& areay=(*localMF[AREA_MF+1])[mfi];
    FArrayBox& areaz=(*localMF[AREA_MF+AMREX_SPACEDIM-1])[mfi];
    FArrayBox& volfab=(*localMF[VOLUME_MF])[mfi];

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     // declared in: LEVELSET_3D.F90
    fort_curvstrip(
     local_caller_string.c_str(),
     local_caller_string.size(),
     &vof_height_function,
     &level,
     &finest_level,
     &curv_min_local[tid_current],
     &curv_max_local[tid_current],
     &nhistory,
     histfab.dataPtr(),
     ARLIM(histfab.loVect()),ARLIM(histfab.hiVect()),
     maskcov.dataPtr(),
     ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
     volfab.dataPtr(),ARLIM(volfab.loVect()),ARLIM(volfab.hiVect()),
     areax.dataPtr(),ARLIM(areax.loVect()),ARLIM(areax.hiVect()),
     areay.dataPtr(),ARLIM(areay.loVect()),ARLIM(areay.hiVect()),
     areaz.dataPtr(),ARLIM(areaz.loVect()),ARLIM(areaz.hiVect()),
     maskfab.dataPtr(),
     ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
     lsfab.dataPtr(),
     ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
     reconfab.dataPtr(),
     ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
     curvfab.dataPtr(),
     ARLIM(curvfab.loVect()),ARLIM(curvfab.hiVect()),
     velfab.dataPtr(),
     ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
     denfab.dataPtr(),
     ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
     tilelo,tilehi,
     fablo,fabhi, 
     &bfact_space,
     &bfact_grid,
     &NS_geometry_coord,
     xlo,dx,
     &cur_time_slab,
     &visc_coef,
     &unscaled_min_curvature_radius,
     &num_curv,
     &ngrow_distance);
  } // mfi
} //omp
  ns_reconcile_d_num(LOOP_CURVSTRIP,"makeStateCurv");

  for (int tid=1;tid<thread_class::nthreads;tid++) {
    if (curv_min_local[tid]<curv_min_local[0])
     curv_min_local[0]=curv_min_local[tid];
    if (curv_max_local[tid]>curv_max_local[0])
     curv_max_local[0]=curv_max_local[tid];
  } // tid

  ParallelDescriptor::ReduceRealMin(curv_min_local[0]);
  ParallelDescriptor::ReduceRealMax(curv_max_local[0]);
  if (curv_min_local[0]<curv_min[0])
    curv_min[0]=curv_min_local[0];
  if (curv_max_local[0]>curv_max[0])
    curv_max[0]=curv_max_local[0];

  localMF[DIST_CURV_MF]->FillBoundary(geom.periodicity());

  if ((fab_verbose==1)||(fab_verbose==3)) {

    std::cout << "c++ level,finest_level " << level << ' ' <<
     finest_level << '\n';
    std::cout << "c++ ngrow_distance " << ngrow_distance << '\n';
    std::cout << "curv_min_local(HT)= " << curv_min_local[0] << '\n';
    std::cout << "curv_max_local(HT)= " << curv_max_local[0] << '\n';

    if (thread_class::nthreads<1)
     amrex::Error("thread_class::nthreads invalid");
    thread_class::init_d_numPts(CL_velocity->boxArray().d_numPts());

    for (MFIter mfi(*CL_velocity,false); mfi.isValid(); ++mfi) {
     BL_ASSERT(grids[mfi.index()] == mfi.validbox());
     const Box& tilegrid = mfi.tilebox();

     int tid_current=ns_thread();
     if ((tid_current<0)||(tid_current>=thread_class::nthreads))
      amrex::Error("tid_current invalid");
     thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     const int gridno = mfi.index();
     const Box& fabgrid = grids[gridno];
     const int* fablo=fabgrid.loVect();
     const int* fabhi=fabgrid.hiVect();
     const Real* xlo = grid_loc[gridno].lo();
     std::cout << "gridno= " << gridno << '\n';
     int interior_only=0;

     std::cout << "output of curvfab" << '\n';
     FArrayBox& curvfab=(*localMF[DIST_CURV_MF])[mfi];
     tecplot_debug(curvfab,xlo,fablo,fabhi,dx,-1,0,0,num_curv,interior_only);

     std::cout << "output of lsfab (LEVELPC)" << '\n';
     FArrayBox& lsfab=(*localMF[LEVELPC_MF])[mfi];
     tecplot_debug(lsfab,xlo,fablo,fabhi,dx,-1,0,0,num_materials,interior_only);

    } // mfi
    ns_reconcile_d_num(LOOP_TECPLOT_DEBUG_CURV,"makeStateCurv");

  } // ((fab_verbose==1)||(fab_verbose==3))

  delete CL_velocity;
  delete den;

 } else
   amrex::Error("project_option invalid10");

}  // subroutine makeStateCurv

//dir=0..sdim-1
MultiFab* NavierStokes::getStateMAC(int ngrow,int dir,Real time) {

 if ((dir<0)||(dir>=AMREX_SPACEDIM))
  amrex::Error("dir invalid get state mac");

 MultiFab& S_new=get_new_data(Umac_Type+dir,slab_step+1);
 int ntotal=S_new.nComp();

 if (ntotal!=1)
  amrex::Error("ntotal invalid");

  //ncomp=1
 MultiFab* mf = new MultiFab(state[Umac_Type+dir].boxArray(),dmap,1,
  ngrow,MFInfo().SetTag("mf getStateMAC"),FArrayBoxFactory());

  //scomp_dest=0 scomp_source=0 ncomp=1
 FillPatch(*this,*mf,0,time,Umac_Type+dir,0,1,debug_fillpatch);

 ParallelDescriptor::Barrier();

 return mf;

}  // subroutine getStateMAC

}/* namespace amrex */

