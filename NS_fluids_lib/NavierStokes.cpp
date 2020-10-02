// tensor, num_materials_viscoelastic, Tensor_Type, im_elastic_map,
// Tensor_new
//#include <winstd.H>

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

/* 
  Narrow Band WENO LEVEL SET METHOD:
  1. t=0 distance function is given.
  2. traverse grid and add (i,j,k) indices for extended 
     narrow band points and regular
     narrow band points.
     For this step, 
     a) two Particle Container objects are created: extended_narrow_band_pc and
     regular_narrow_band_pc 
       NStructReal=0
       NStructInt=0
       NArrayReal=0
       NArrayInt=sdim   (i,j,k)
     b) traverse grid and add to the respective particle container objects.

      for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {

       // ``particles'' starts off empty
       auto& particles = GetParticles(lev)[std::make_pair(mfi.index(),
                                        mfi.LocalTileIndex())];

       (optional) for all indices on tile:

       ParticleType p;
       p.id()   = ParticleType::NextID();
       p.cpu()  = ParallelDescriptor::MyProc();
       p.pos(0) = ...
       etc...

       // AoS real data
       p.rdata(0) = ...
       p.rdata(1)  = ...

       // AoS int data
       p.idata(0) = ...
       p.idata(1) = ...

       // Particle real attributes (SoA)
       std::array<double, 2> real_attribs;
       real_attribs[0] = ...
       real_attribs[1] = ...

       // Particle int attributes (SoA)
       std::array<int, 2> int_attribs;
       int_attribs[0] = ...
       int_attribs[1]  = ...

       particles.push_back(p);
       particles.push_back_real(real_attribs);
       particles.push_back_int(int_attribs);

       // ... add more particles if desired ...
      }

  3. using MyParIter = ParIter<0, 0, 0, sdim>;
       // loop through every grid (or tile if tiling enabled) that has particles.
     for (MyParIter pti(regular_narrow_band_pc, lev); pti.isValid(); ++pti) {
      auto& particle_attributes = pti.GetStructOfArrays();
      // Vector<Real>& real_comp0 = particle_attributes.GetRealData(0);
      for (int i=0;i<SDIM;i++) 
       Vector<int>& int_comp[i]  = particle_attributes.GetIntData(i);
      for (int i = 0; i < pti.numParticles; ++i) {
        // do stuff with your SoA data... (int_comp[j], j=0..sdim-1)
	advect LS data with WENO.
      }
     }

     for (MyParIter pti(extended_narrow_band_pc, lev); pti.isValid(); ++pti) {
      auto& particle_attributes = pti.GetStructOfArrays();
      // Vector<Real>& real_comp0 = particle_attributes.GetRealData(0);
      for (int i=0;i<SDIM;i++)
       Vector<int>& int_comp[i]  = particle_attributes.GetIntData(i);
      for (int i = 0; i < pti.numParticles(); ++i) {
        // do stuff with your SoA data... (int_comp[j], j=0..sdim-1)
        reinitialize LS data with WENO.
      }
      FORT_DOSTUFF_WITH_SOA(
			real_comp0.size(),  // =pti.numParticles()
			real_comp0.dataPtr(),
			int_comp[0].dataPtr(), 
			int_comp[1].dataPtr(), 
			int_comp[AMREX_SPACEDIM-1].dataPtr(), .... );

         or
      Array<int> int_compALL(AMREX_SPACEDIM * pti.numParticles())
      for (int i=0;i<SDIM;i++) {
       for (int j = 0; j < pti.numParticles; ++j) {
	       int k=i*pti.numParticles()+j;
	       int_compALL[k]=int_comp[i][j];
       }
      }

     }


  https://amrex-codes.github.io/amrex/docs_html/Particle.html

  if particles:
  #include <AMReX_Particles.H>
  NStructReal=number of extra Real variables (not including particle position)  
  NStructInt=number of extra int variables (not including cpu and id)

  Array-of-Structs: particle1, particle2, particle3, ....
  Struct-of-Arrays: foo1,foo2,foo3, ...  NArrayReal=2
                    bar1,bar2,bar3, ... 
		    l1,l2,l3, ....   NArrayInt=2
		    n1,n2,n3, ....
  ParticleContainer<NStructReal,NStructInt,NArrayReal,NArrayInt> mypc

    (see AMReX_Particles.H)
    rr[n]=refinement ratio between levels n and n+1
  Particle(const Vector<Geometry> &geom,
           const Vector<DistributionMapping> &dmap,
	   const Vector<BoxArray> &ba,
	   const Vector<int> &rr);  
 
  Redistribute()

  amrex/Tutorials/Particles/NeighborList


  for Cody,
    1. declare mypc (object of type ParticleContainer)
    2. fill the particle container
    3. Redistribute()
    4. advect LS
    5. advect particles
    6. Redistribute()
    7. correct LS

    (how to use the GPUs?)
*/

#include <NavierStokes.H>
#include <INTERP_F.H>
#include <MACOPERATOR_F.H>
#include <NAVIERSTOKES_F.H>
#include <TECPLOTUTIL_F.H>
#include <GODUNOV_F.H>
#include <MASS_TRANSFER_F.H>
#include <PROB_F.H>
#include <MOF_F.H>
#include <PLIC_F.H>
#include <LEVEL_F.H>
#include <MOF_REDIST_F.H>
#include <ShallowWater_F.H>
#include <SOLIDFLUID_F.H>
#include <DERIVE_F.H>
#include <MG_F.H>

#ifdef MVAHABFSI
#include <CTMLFSI_F.H>
#endif

namespace amrex{

#define bogus_value 1.e20
#define show_norm2_flag 0
#define mf_check_inf_bounds 1

#define DEFAULT_MOFITERMAX 15

#define debug_PC 1
//
// Static objects.
//
BCRec NavierStokes::phys_bc;
BCRec NavierStokes::temperature_phys_bc;
BCRec NavierStokes::species_phys_bc;

int  NavierStokes::profile_debug=0;
bool NavierStokes::ns_tiling=false;

int NavierStokes::POLYGON_LIST_MAX=1000;

int  NavierStokes::nfluxSEM=0;
int  NavierStokes::nstate_SDC=0;
int  NavierStokes::ns_time_order=1; // time_blocking_factor
int  NavierStokes::slab_step=0;
int  NavierStokes::SDC_outer_sweeps=0;
int  NavierStokes::divu_outer_sweeps=0;
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
int  NavierStokes::check_nan    = 0;
// 1=curv 2=error heat 3=both
int  NavierStokes::fab_verbose  = 0;
int  NavierStokes::output_drop_distribution = 0;
int  NavierStokes::extend_pressure_into_solid = 0;
Real NavierStokes::cfl          = 0.5;
int  NavierStokes::MOF_TURN_OFF_LS=0;
int  NavierStokes::MOF_DEBUG_RECON=0;
int  NavierStokes::MOFITERMAX=DEFAULT_MOFITERMAX;

/*
 continuous_mof=0

  regular MOF  minimize E=||x_ij^ref-x_ij^derived||
  subject to the constraint that F_ij^ref=F_ij^derived

   x_ij^ref=reference centroid in cell ij
   x_ij^derived=derived centroid in cell ij for a given slope and
     intercept.
   F_ij^ref=reference volume fraction in cell ij
   F_ij^derived=derived volume fraction in cell ij for a given slope and
     intercept.   


continuous_mof=2 (if same number of materials in center cell as in stencil)

  CMOF  minimize E=||xS_ij^ref-xS_ij^derived||  "S"=super cell
  subject to the constraint that F_ij^ref=F_ij^derived

   xS_ij^ref=reference centroid in cell stencil i'=i-1,..,i+1,
     j'=j-1,..,j+1

   xS_ij^derived=derived centroid in cell stencil for a given slope and
     intercept. 
   F_ij^ref=reference volume fraction in cell
   F_ij^derived=derived volume fraction in cell for a given
     slope and intercept.

continuous_mof=3 

  use CLSVOF in 2 material cells and MOF in >2 mat cells.

continuous_mof=4 

  use CLSVOF in 2 material cells and CMOF in >2 mat cells.

continuous_mof=5 

  use CLSVOF everywhere.


NOTE: rigid materials are not counted as materials in a cell.  Rigid 
materials are immersed into the fluid(s). 

future work:
use machine learning techniques in order to quickly derive an initial
guess for the MOF reconstruction.  Also, getting the intercept given
the volume fraction and slope.

calibrate sub-scale model parameters for any given problem, so that 
certain basic benchmark tests give the expected results. (e.g. drag
on an interface, growth rate of perturbations, pressure drop)
*/

int  NavierStokes::VOF_reflux=0;
int  NavierStokes::continuous_mof=0;

// 0  low order space and time
// 1  SEM space and time
// 2  SEM space 
// 3  SEM time
int  NavierStokes::enable_spectral=0;
int  NavierStokes::viscous_enable_spectral=0;
int  NavierStokes::projection_enable_spectral=0;
int  NavierStokes::SEM_upwind=1;
//0=div(uS)-S div(u)    1=u dot grad S
int  NavierStokes::SEM_advection_algorithm=0;
// default: tessellating fluid => default==1
//          non-tesselating or tesselating solid => default==0
Vector<int> NavierStokes::truncate_volume_fractions; 

// default=1
Vector<int> NavierStokes::particle_nsubdivide; 
Vector<int> NavierStokes::particle_max_per_nsubdivide; 
Vector<int> NavierStokes::particleLS_flag; 
int NavierStokes::NS_ncomp_particles=0;

Real NavierStokes::truncate_thickness=2.0;  
Real NavierStokes::init_shrink  = 1.0;
Real NavierStokes::change_max   = 1.1;
Real NavierStokes::change_max_init = 1.1;
Real NavierStokes::fixed_dt     = 0.0;
Real NavierStokes::fixed_dt_init = 0.0;
Real NavierStokes::min_velocity_for_dt = 1.0e-12;
Real NavierStokes::fixed_dt_velocity = 0.0;
Real NavierStokes::dt_max       = 1.0e+10;
Real NavierStokes::MUSHY_THICK  = 2.0;
Real NavierStokes::gravity      = 0.0;
// terminal_velocity_dt==1 =>
// use the terminal velocity for CFL condition instead 
// of the default condition: u dt < dx  u=g dt  g dt^2 < dx   dt<sqrt(dx/g)
int NavierStokes::terminal_velocity_dt = 0;
int NavierStokes::gravity_dir = AMREX_SPACEDIM;
int NavierStokes::invert_gravity = 0;
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

int  NavierStokes::visual_tessellate_vfrac=0;   
int  NavierStokes::visual_revolve=0;   
int  NavierStokes::visual_option=-2; // -2 zonal tecplot,-1 plot files (visit)

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

Vector< Vector<Real> > NavierStokes::min_face_wt;
Vector< Vector<Real> > NavierStokes::max_face_wt;

Vector< Vector<Real> > NavierStokes::DVOF;
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

int NavierStokes::tecplot_max_level=0;
int NavierStokes::max_level_two_materials=0;

// default=0.
// 0=> never adapt  -1=> tag cell for AMR if owned by material in question.
// otherwise, if radius<radius_cutoff * dx then adapt.
// fluid-fluid interfaces are always adapted.
Vector<int> NavierStokes::radius_cutoff;

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

// search num_materials,AMREX_SPACEDIM+1,SDIM+1,idenbase,ipres_base,
//  iden_base,scomp_mofvars,nstate=,scomp_den,pressure_comp,dcomp,
//  pcomp,tcomp,scomp,scomp_pres,get_mm_scomp_solver,dencomp,scomp_tensor,
//  im_pres,velcomp,prescomp,flagcomp
int  NavierStokes::num_materials=0;
int  NavierStokes::num_materials_vel=1;
int  NavierStokes::num_materials_scalar_solve=1;

int  NavierStokes::use_supermesh=0;
int  NavierStokes::ncomp_sum_int_user=0;

// set using elastic_viscosity
int  NavierStokes::num_materials_viscoelastic=0;

int  NavierStokes::num_state_material=SpeciesVar; // den,T
int  NavierStokes::num_state_base=SpeciesVar; // den,T
int  NavierStokes::ngeom_raw=AMREX_SPACEDIM+1;
int  NavierStokes::ngeom_recon=NUM_MOF_VAR;

int  NavierStokes::State_Type=0;
int  NavierStokes::Umac_Type=1;
int  NavierStokes::Vmac_Type=2;
int  NavierStokes::Wmac_Type=AMREX_SPACEDIM;
int  NavierStokes::LS_Type=AMREX_SPACEDIM+1;
int  NavierStokes::DIV_Type=AMREX_SPACEDIM+2;
int  NavierStokes::Solid_State_Type=AMREX_SPACEDIM+3;
int  NavierStokes::Tensor_Type=AMREX_SPACEDIM+4;
int  NavierStokes::NUM_STATE_TYPE=AMREX_SPACEDIM+5;

// 0=velocity stored at cells  
// 1=velocity stored at faces
int  NavierStokes::face_flag=0;

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
// 0,1 => viscoelastic FENE-CR material  2=> elastic material
Vector<int> NavierStokes::viscoelastic_model; // def=0
Vector<int> NavierStokes::les_model; // def=0
// temperature_primitive_variable defaults to 0 (conservative) for
// inviscid compressible materials, and 1 (non conservative) for other
// materials.
// E=u dot u/2 + e
// e=cv T (compressible)   cp T (incompressible)
// 0=> conservative advection of temperature:
//    (rho E)_t + div(rho u E + up)=div(k grad T) 
//    it is assumed the fluid is inviscid and 
//    compressible for this option.
// 1=> non-conservative advection of temperature:
//    (rho e)_t + div(rho u e)+div(u)p=div(k grad T)
// if u is changed, then T must be changed if energy is to be conserved:
//  u_new^2/2 + cv T_new = u_old^2/2 + cv T_old
//  Tnew=Told+(1/cv)(u_old^2/2 - u_new^2/2)
Vector<int> NavierStokes::temperature_primitive_variable; 

Vector<Real> NavierStokes::elastic_viscosity; // def=0

Vector<Real> NavierStokes::Carreau_alpha; // def=1
Vector<Real> NavierStokes::Carreau_beta; // def=0
Vector<Real> NavierStokes::Carreau_n; // def=1
Vector<Real> NavierStokes::Carreau_mu_inf; // def=0

Vector<Real> NavierStokes::concentration; // def=0
Vector<Real> NavierStokes::etaL; // def=0 (etaL0)
Vector<Real> NavierStokes::etaS; // def=0
Vector<Real> NavierStokes::etaP; // def=0 (etaP0)
Vector<Real> NavierStokes::polymer_factor; // def=0

 // 0 - centroid furthest from uncaptured centroid
 // 1 - use MOF error
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
// if make_interface_incomp==1:
// density equation is Drho/Dt=0 in the incompressible zone.
// if make_interface_incomp==2:
// density equation is rhot_t+div(rho u)=0 in the incompressible zone.
int  NavierStokes::make_interface_incomp=0; 
Vector<int> NavierStokes::advection_order; // def=1
// def=advection_order
Vector<int> NavierStokes::density_advection_order; 

// 0=Sussman and Puckett algorithm 
// 1=EILE (default), -1=Weymouth Yue
// 2=always EI   3=always LE
int NavierStokes::EILE_flag=1;

// this should always be 0.  Directionally split is better.
// In future: implement super mesh algorithm and extend the projected velocity
// to the faces.
int NavierStokes::unsplit_flag=0;

// 0=no limiter 1=minmod 2=minmod with slope=0 at interface
// 3=no limit, but slope=0 at interface
int  NavierStokes::slope_limiter_option=2; 
int  NavierStokes::bicgstab_max_num_outer_iter=60;
Real NavierStokes::projection_pressure_scale=1.0;
Real NavierStokes::projection_velocity_scale=1.0;

int NavierStokes::curv_stencil_height=4; 
int NavierStokes::ngrow_distance=4;
int NavierStokes::ngrow_make_distance=3;

// blob_matrix,blob_RHS,blob_velocity,
// blob_integral_momentum,blob_energy,
// blob_mass_for_velocity (3 components)
// blob_volume, 
// blob_center_integral,blob_center_actual
// blob_perim, blob_perim_mat, blob_triple_perim, 
int NavierStokes::num_elements_blobclass=0;

int NavierStokes::ngrowFSI=3;
int NavierStokes::nFSI_sub=12; //velocity+LS+temperature+flag+stress (3D)
Vector<int> NavierStokes::im_solid_map; //nparts components, in range 0..nmat-1
Vector<int> NavierStokes::im_elastic_map; 

int NavierStokes::ngrow_expansion=2;
 
Real NavierStokes::real_number_of_cells=0.0; 

// 1.0/(den_max * mglib_min_coeff_factor) default=1000.0
Real NavierStokes::mglib_min_coeff_factor=1000.0; 

int NavierStokes::hydrate_flag=0; 
int NavierStokes::singular_possible=0; 
int NavierStokes::solvability_projection=0; 
int NavierStokes::local_solvability_projection=0; 
int NavierStokes::post_init_pressure_solve=1; 

int NavierStokes::conservative_tension_force=0;
// 0=> I(u_GL) dot grad u_GG, no sync project
// 1=> grad dot (u_GL u_GG) - u_GG div u_GL, no sync project
// 2=> grad dot (u_GL u_GG) - u_GG div u_GL, sync project
int NavierStokes::conservative_div_uu=1;

Vector<Real> NavierStokes::tension_slope;
Vector<Real> NavierStokes::tension_min;
Vector<Real> NavierStokes::tension_T0;
Vector<Real> NavierStokes::tension;
Vector<Real> NavierStokes::prefreeze_tension;
Vector<Real> NavierStokes::recalesce_model_parameters;

Vector<Real> NavierStokes::outflow_velocity_buffer_size;

Vector<Real> NavierStokes::cap_wave_speed;

Vector<Real> NavierStokes::hardwire_Y_gamma;
Vector<Real> NavierStokes::hardwire_T_gamma;
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

Vector<Real> NavierStokes::cavitation_pressure;
Vector<Real> NavierStokes::cavitation_vapor_density;
Vector<Real> NavierStokes::cavitation_tension;

// 1.. num_species_var
Vector<Real> NavierStokes::species_evaporation_density;

Vector<Real> NavierStokes::nucleation_pressure;
Vector<Real> NavierStokes::nucleation_pmg;
Vector<Real> NavierStokes::nucleation_mach;
Vector<Real> NavierStokes::nucleation_temp;
Real NavierStokes::nucleation_period=0.0;
Real NavierStokes::nucleation_init_time=0.0;
int NavierStokes::n_sites=0;
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
// 4=source term model (single equation for T with source term).
//   Tanasawa model or Schrage is used for evaporation and condensation.
//   TSAT used to determine if phase change happens.
//   expansion source, and offsetting sink evenly distributed.
//   For Tannasawa model implemented here, it is assumed that Y=1
//   at the interface.  Schrage Model does not assume Y=1 at interface.
// MEHDI EVAPORATION: freezing_model=4
//   1=liquid  2=ambient gas  3=solid wall  num_species_var=1
//   ->  12 13 23 21 31 32
//   latent_heat
//   ->  +L 0  0  -L 0  0
//   freezing_model
//   ->  4  0  0  4  0  0
//   mass_fraction_id
//   ->  1  0  0  1  0  0
//   Tanasawa_or_Schrage
//   ->  2  0  0  2  0  0
//   speciesviscconst
//   ->  0.0 D 0.0
//   distribute_from_target (distribution of the expansion term)
//   ->  0 0 0 1 0 0
//   molar_mass
//   ->  mL  m_ambient 0.0
//   species_molar_mass
//   ->  m_vapor
//   species_evaporation_density
//   ->  density_vapor
//   solvability_projection=0
//   material_type
//   0 ?? 999
//   temperature_primitive_variable
//   1 1 1
//
// 5=evaporation/condensation (Stefan model speed)
// 6=evaporation/condensation (Palmore and Desjardins, JCP 2019)
// 7=cavitation
Vector<int> NavierStokes::freezing_model;
Vector<int> NavierStokes::Tanasawa_or_Schrage; //1=Tanasawa  2=Schrage
//ispec=mass_fraction_id[0..2 nten-1]=1..num_species_var
Vector<int> NavierStokes::mass_fraction_id; 
//link diffused material to non-diff. (array 1..num_species_var)
//spec_material_id_LIQUID, spec_material_id_AMBIENT are 
//both input AND derived.
Vector<int> NavierStokes::spec_material_id_LIQUID; 
Vector<int> NavierStokes::spec_material_id_AMBIENT; 
// 0 - distribute to the destination material (default)
//     V=mdot/rho_src
// 1 - distribute to the source material
//     V=mdot/rho_dst
Vector<int> NavierStokes::distribute_from_target;
int NavierStokes::is_phasechange=0;
int NavierStokes::normal_probe_size=1;
// 0=dirichlet at inflow
// 1=dirichlet at inflow and outflow
// 2=dirichlet at inflow and walls.
// 3=dirichlet at inflow, outflow, and walls.
int NavierStokes::prescribe_temperature_outflow=0; // default is 0

int  NavierStokes::use_lsa=0;
Real NavierStokes::Uref=0.0;
Real NavierStokes::Lref=0.0;

Real NavierStokes::pgrad_dt_factor=1.0;
// 0=volume fraction  1=mass fraction 2=impedance fraction
int  NavierStokes::pressure_select_criterion=0;

int  NavierStokes::last_finest_level=-1;

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
Vector<Real> NavierStokes::density_floor_expansion;  // def=denconst
Vector<Real> NavierStokes::density_ceiling_expansion;  // def=denconst
Vector<Real> NavierStokes::molar_mass;  // def=1
Vector<Real> NavierStokes::denconst;
Real NavierStokes::denconst_max=0.0;
Real NavierStokes::denconst_min=0.0;
Vector<Real> NavierStokes::denconst_interface;
Vector<Real> NavierStokes::denconst_gravity; // def=1.0
int NavierStokes::stokes_flow=0;
int NavierStokes::cancel_advection=0;

// passed to MAC_TO_CELL, CELL_TO_MAC, 
//  VFRAC_SPLIT, VFRAC_UNSPLIT
Vector<Real> NavierStokes::added_weight; // def=1.0

Vector<Real> NavierStokes::stiffPINF;
Vector<Real> NavierStokes::prerecalesce_stiffCP;  // def=4.1855E+7
Vector<Real> NavierStokes::prerecalesce_stiffCV;  // def=4.1855E+7
Vector<Real> NavierStokes::stiffCP;  // def=4.1855E+7
Vector<Real> NavierStokes::stiffCV;  // def=4.1855E+7
Vector<Real> NavierStokes::stiffGAMMA;

int NavierStokes::constant_viscosity=0;

Real NavierStokes::angular_velocity=0.0;
Vector<Real> NavierStokes::DrhoDT;  // def=0.0
Vector<Real> NavierStokes::DrhoDz;  // def=0.0

 // 1=>rho=rho(T,Y,z)
 // 2=>P_hydro=P_hydro(rho(T,Y,z)) (Boussinesq like approximation)
Vector<int> NavierStokes::override_density; // def=0
Vector<Real> NavierStokes::prerecalesce_viscconst;
Vector<Real> NavierStokes::viscconst;
Real NavierStokes::viscconst_max=0.0;
Real NavierStokes::viscconst_min=0.0;
Vector<Real> NavierStokes::viscconst_eddy;
Vector<Real> NavierStokes::speciesviscconst;// species mass diffusion coeff.
Vector<Real> NavierStokes::prerecalesce_heatviscconst;
Vector<Real> NavierStokes::heatviscconst;
Real NavierStokes::heatviscconst_max=0.0;
Real NavierStokes::heatviscconst_min=0.0;
Vector<Real> NavierStokes::viscconst_interface;
Vector<Real> NavierStokes::heatviscconst_interface;
Vector<Real> NavierStokes::speciesconst;  // unused currently
Vector<Real> NavierStokes::speciesviscconst_interface;
// 1..num_species_var
Vector<Real> NavierStokes::species_molar_mass; // def=1
// 0=diffuse in solid 1=dirichlet 2=neumann
int NavierStokes::solidheat_flag=0; 
int NavierStokes::diffusionface_flag=1; // 0=use LS  1=use VOF
int NavierStokes::elasticface_flag=1; // 0=use LS  1=use VOF
int NavierStokes::temperatureface_flag=1; // 0=use LS  1=use VOF

Vector<int> NavierStokes::material_type;

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

int NavierStokes::gmres_precond_iter_base=4; 

int NavierStokes::smooth_type = 2; // 0=GSRB 1=ICRB 2=ILU  3=Jacobi
int NavierStokes::bottom_smooth_type = 2; // 0=GSRB 1=ICRB 2=ILU 3=Jacobi
int NavierStokes::use_mg_precond_in_mglib=1;
Real NavierStokes::bottom_bottom_tol_factor=0.01;

// 0=> u=u_solid if phi_solid>=0
// 1=> u=u_solid_ghost if phi_solid>=0
// 2=> generalized Navier Boundary condition (GNBC),
//   for conventional contact line dynamics, 
//   modify "get_use_DCA" in PROB.F90.
int NavierStokes::law_of_the_wall=0;
int NavierStokes::ZEYU_DCA_SELECT=-1; // -1 = static angle

// 0 fluid, tessellating (default)
// 1 prescribed rigid solid, non-tessellating (PROB.F90)
// 2 prescribed rigid solid, non-tessellating (sci_clsvof.F90)
// 3 FSI ice,tessellating
// 4 FSI, non-tessellating link w/Kourosh Shoele
// 5 FSI rigid solid, tessellating (PROB.F90)
// 6 FSI ice, tessellating (initial geometry: sci_clsvof.F90)
// 7 fluid, tessellating (initial geometry: sci_clsvof.F90)
Vector<int> NavierStokes::FSI_flag; 
Vector<int> NavierStokes::FSI_touch_flag; // 0..nthreads-1
// default: 1
Vector<int> NavierStokes::FSI_refine_factor; 
// default: 3
Vector<int> NavierStokes::FSI_bounding_box_ngrow; 

int NavierStokes::CTML_FSI_numsolids = 0;
int NavierStokes::CTML_force_model = 0; // 0=Lag force 1=Lag stress
int NavierStokes::CTML_FSI_init = 0;

int NavierStokes::invert_solid_levelset = 0; 
int NavierStokes::elements_generated = 0; 

// 0=take into account sound speed only at t=0 if compressible.
// 1=always take into account sound speed
// 2=never take into account sound speed
Vector<int> NavierStokes::shock_timestep; 

Real NavierStokes::visc_coef=0.0;

int NavierStokes::include_viscous_heating=0;

int NavierStokes::multilevel_maxcycle=200;

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
int NavierStokes::viscous_maxiter = 1;
int NavierStokes::always_use_bicgstab = 0;
Real NavierStokes::total_advance_time=0.0;

int NavierStokes::curv_index=0;
int NavierStokes::pforce_index=1;
int NavierStokes::faceden_index=2;
int NavierStokes::facecut_index=3;
int NavierStokes::icefacecut_index=4;
int NavierStokes::icemask_index=5;
int NavierStokes::facevisc_index=6;
int NavierStokes::faceheat_index=7;
int NavierStokes::facevel_index=8;
int NavierStokes::facespecies_index=9;
int NavierStokes::massface_index=10;
int NavierStokes::vofface_index=11;
int NavierStokes::ncphys=12;

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

// ns.mof_ordering overrides this.
void mof_ordering_override(Vector<int>& mof_ordering_local,
 int nmat,int probtype,int axis_dir,Real radblob3,
 Real radblob4,Real radblob7,
 int mof_error_ordering_local,
 Vector<int> FSI_flag_temp) {

 if (nmat!=mof_ordering_local.size())
  amrex::Error("mof_ordering_local invalid size");
 if (nmat!=FSI_flag_temp.size())
  amrex::Error("FSI_flag_temp invalid size");
 if (nmat<1)
  amrex::Error("nmat out of range");

 for (int im=0;im<nmat;im++) {
  mof_ordering_local[im]=0;

  if ((FSI_flag_temp[im]==0)||
      (FSI_flag_temp[im]==7)) { // fluid
   // do nothing, tessellating
  } else if (FSI_flag_temp[im]==1) { // prescribed rigid solid (PROB.F90)
   mof_ordering_local[im]=1;  // non-tessellating
  } else if (FSI_flag_temp[im]==2) { // prescribed rigid solid (sci_clsvof.F90)
   mof_ordering_local[im]=1;  // non-tessellating
  } else if ((FSI_flag_temp[im]==3)||
             (FSI_flag_temp[im]==6)) { // ice (PROB.F90),ice (sci_clsvof.F90)
   // do nothing, tessellating
   
   // FSI elastic link w/Kourosh (sci_clsvof.F90)
  } else if (FSI_flag_temp[im]==4) {  
   mof_ordering_local[im]=1; // non-tessellating
  } else if (FSI_flag_temp[im]==5) { // FSI rigid solid (PROB.F90)
   mof_ordering_local[im]=1; // tessellating
  } else
   amrex::Error("FSI_flag_temp invalid");

 }

  // default: centroid farthest from uncaptured centroid.
 if (mof_error_ordering_local==0) { 

  for (int im=0;im<nmat;im++) {

   if ((FSI_flag_temp[im]==0)||
       (FSI_flag_temp[im]==7)) { // fluid, tessellating
    mof_ordering_local[im]=nmat;
   } else if (FSI_flag_temp[im]==1) { //prescribed rigid solid (PROB.F90)
    mof_ordering_local[im]=1; // non-tessellating
   } else if (FSI_flag_temp[im]==2) { //prescribed rigid solid (sci_clsvof.F90)
    mof_ordering_local[im]=1; // non-tessellating
   } else if ((FSI_flag_temp[im]==3)||
	      (FSI_flag_temp[im]==6)) { // ice (PROB.F90),ice (sci_clsvof.F90)
    mof_ordering_local[im]=nmat; // tessellating
   } else if (FSI_flag_temp[im]==4) { // FSI link w/Kourosh (sci_clsvof.F90)
    mof_ordering_local[im]=1;  // non-tessellating
   } else if (FSI_flag_temp[im]==5) { // FSI rigid solid (PROB.F90)
    mof_ordering_local[im]=1;  // tessellating
   } else
    amrex::Error("FSI_flag_temp invalid");

  } // im

   // impinge jets unlike material
  if ((probtype==530)&&(AMREX_SPACEDIM==3)) {
   if (axis_dir==1)
    mof_ordering_local[1]=nmat+1;  // make gas have low priority
   else if (axis_dir!=0)
    amrex::Error("axis_dir invalid probtype=530");
  }

   // ns.mof_ordering overrides this.

   // 2d colliding droplets, boiling, freezing problems
  if ((probtype==55)&&(AMREX_SPACEDIM==2)) {
   if (radblob7>0.0)
    mof_ordering_local[1]=nmat+1;  // make gas have low priority
   if (axis_dir==0) {
    // do nothing
   } else if (axis_dir==1) {
    // 0=water 1=gas 2=ice 3=cold plate
    mof_ordering_local[2]=1;
   } else if (axis_dir==5) {
    // 0=water 1=gas 2=ice 3=cold plate
    mof_ordering_local[2]=1;
    // 0=water 1=vapor 2=hot plate or
    // 0=water 1=vapor 2=gas 3=hot plate 
   } else if (axis_dir==6) {  // nucleate boiling incompressible
    mof_ordering_local[nmat-1]=1;
    // 0=water 1=vapor 2=hot plate or
    // 0=water 1=vapor 2=gas 3=hot plate 
   } else if (axis_dir==7) {  // nucleate boiling compressible
    mof_ordering_local[nmat-1]=1;
   } else
    amrex::Error("axis_dir invalid probtype==55");
  }

  if (probtype==540) {
   if ((radblob4>0.0)&&(radblob3>0.0)) {
    amrex::Error("conflict of parametrs for 540");
   } else if (radblob3>0.0) {  
    mof_ordering_local[1]=nmat+1;  // make gas have low priority
   } else if (radblob4>0.0) {
    mof_ordering_local[2]=nmat+1;  // make filament gas have low priority
   }
  }

  if (probtype==202) {  // liquidlens
   mof_ordering_local[1]=1; // make (circle) material 2 have high priority
  }

  if ((probtype==17)&&(nmat==3)&&(1==0)) {  // droplet impact 3 materials
   mof_ordering_local[1]=1; // make gas material 2 have high priority
  }

 } else if (mof_error_ordering_local==1) {

  // do nothing, order=0 
  // (except if FSI_flag_temp[im]==1,2,4,5 )
  // FSI_flag=1,2,4 non-tessellating
  // FSI_flag=0,3,5,6,7  tessellating

 } else
  amrex::Error("mof_error_ordering invalid");

} // end subroutine mof_ordering_override

void read_geometry_raw(int& geometry_coord,
	Vector<Real>& geometry_prob_lo,
        Vector<Real>& geometry_prob_hi,
	Vector<int>& geometry_is_periodic,
	int& geometry_is_any_periodic) {

 ParmParse pp("geometry");
 pp.get("coord_sys",geometry_coord);
 pp.getarr("prob_lo",geometry_prob_lo,0,AMREX_SPACEDIM);
 pp.getarr("prob_hi",geometry_prob_hi,0,AMREX_SPACEDIM);
 geometry_is_periodic.resize(AMREX_SPACEDIM);
 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  geometry_is_periodic[dir]=0;
 pp.queryarr("is_periodic",geometry_is_periodic,0,AMREX_SPACEDIM);
 geometry_is_any_periodic=0;
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  if (geometry_is_periodic[dir]==0) {
   // do nothing
  } else if (geometry_is_periodic[dir]==1) {
   geometry_is_any_periodic=1;
  } else
    amrex::Error("geometry_is_periodic[dir] invalid");
 } // dir=0...sdim-1

} // subroutine read_geometry_raw

void fortran_parameters() {

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
 Real inflow_pressure=0.0;
 Real outflow_pressure=0.0;
 Real period_time=0.0;

 ParmParse ppmain;
 fort_stop_time=-1.0;
 ppmain.query("stop_time",fort_stop_time);

 int ns_max_level;
 ParmParse ppamr("amr");
 ppamr.get("max_level",ns_max_level);
 Vector<int> ns_space_blocking_factor;
 ns_space_blocking_factor.resize(ns_max_level+1);
 for (int lev=0;lev<=ns_max_level;lev++)
  ns_space_blocking_factor[lev]=2;
 ppamr.queryarr("space_blocking_factor",
   ns_space_blocking_factor,0,ns_max_level+1);

 int time_blocking_factor=1;
 ppamr.query("time_blocking_factor",time_blocking_factor); 

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

 pp.query("xblob2",xblob2);
 pp.query("yblob2",yblob2);
 pp.query("zblob2",zblob2);
 pp.query("radblob2",radblob2);

 pp.query("xblob3",xblob3);
 pp.query("yblob3",yblob3);
 pp.query("zblob3",zblob3);
 pp.query("radblob3",radblob3);

 pp.query("xblob4",xblob4);
 pp.query("yblob4",yblob4);
 pp.query("zblob4",zblob4);
 pp.query("radblob4",radblob4);

 pp.query("xblob5",xblob5);
 pp.query("yblob5",yblob5);
 pp.query("zblob5",zblob5);
 pp.query("radblob5",radblob5);

 pp.query("xblob6",xblob6);
 pp.query("yblob6",yblob6);
 pp.query("zblob6",zblob6);
 pp.query("radblob6",radblob6);

 pp.query("xblob7",xblob7);
 pp.query("yblob7",yblob7);
 pp.query("zblob7",zblob7);
 pp.query("radblob7",radblob7);

 pp.query("xblob8",xblob8);
 pp.query("yblob8",yblob8);
 pp.query("zblob8",zblob8);
 pp.query("radblob8",radblob8);

 pp.query("xblob9",xblob9);
 pp.query("yblob9",yblob9);
 pp.query("zblob9",zblob9);
 pp.query("radblob9",radblob9);

 pp.query("xblob10",xblob10);
 pp.query("yblob10",yblob10);
 pp.query("zblob10",zblob10);
 pp.query("radblob10",radblob10);

 xactive=0.0;
 yactive=0.0;
 zactive=0.0;
 ractive=0.0;
 ractivex=0.0;
 ractivey=0.0;
 ractivez=0.0;

 pp.query("xactive",xactive);
 pp.query("yactive",yactive);
 pp.query("zactive",zactive);
 pp.query("ractive",ractive);
 if (ractive>0.0) {
  ractivex=ractive;
  ractivey=ractive;
  ractivez=ractive;
 }
 pp.query("ractivex",ractivex);
 pp.query("ractivey",ractivey);
 pp.query("ractivez",ractivez);

 pp.get("adv_dir",adv_dir);
 if ((adv_dir<1)||(adv_dir>2*AMREX_SPACEDIM+1))
  amrex::Error("adv_dir invalid");

 pp.get("adv_vel",adv_vel);
 pp.get("rgasinlet",rgasinlet);
 pp.get("vinletgas",vinletgas);
 pp.get("twall",twall);
 pp.get("advbot",advbot);
 pp.query("inflow_pressure",inflow_pressure);
 pp.query("outflow_pressure",outflow_pressure);
 pp.query("period_time",period_time);

 int invert_solid_levelset=0;
 pp.query("invert_solid_levelset",invert_solid_levelset);
 if (!((invert_solid_levelset==1)||(invert_solid_levelset==0)))
  amrex::Error("invert_solid_levelset invalid");

 int num_species_var=0;
 int num_materials=0;
 int num_materials_vel=1;
 int num_materials_scalar_solve=1;
 int num_materials_viscoelastic=0;

 int num_state_material=SpeciesVar;  // den,T
 int num_state_base=SpeciesVar;  // den,T
 int ngeom_raw=AMREX_SPACEDIM+1;
 int ngeom_recon=NUM_MOF_VAR;

 pp.get("num_materials",num_materials);
 if ((num_materials<2)||(num_materials>999))
  amrex::Error("num materials invalid");

 int nmat=num_materials;

 pp.query("num_materials_vel",num_materials_vel);
 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel==1 required");

 pp.query("num_materials_scalar_solve",num_materials_scalar_solve);
 if ((num_materials_scalar_solve!=1)&&(num_materials_scalar_solve!=nmat))
  amrex::Error("num_materials_scalar_solve invalid");

  // this is local variable, not static variable
 int MOFITERMAX=DEFAULT_MOFITERMAX;  
 pp.query("MOFITERMAX",MOFITERMAX);
 if ((MOFITERMAX<0)||(MOFITERMAX>50))
  amrex::Error("mof iter max invalid in navierstokes");

 int MOF_TURN_OFF_LS=0;  // this is local variable, not static variable
 pp.query("MOF_TURN_OFF_LS",MOF_TURN_OFF_LS);
 if ((MOF_TURN_OFF_LS!=0)&&(MOF_TURN_OFF_LS!=1))
  amrex::Error("mof turn off ls invalid in navierstokes");

 int MOF_DEBUG_RECON=0;  // this is local variable, not static variable
 pp.query("MOF_DEBUG_RECON",MOF_DEBUG_RECON);
 if ((MOF_DEBUG_RECON!=0)&&(MOF_DEBUG_RECON!=1)&&
     (MOF_DEBUG_RECON!=2))
  amrex::Error("mof debug recon invalid in navierstokes");

 pp.get("num_species_var",num_species_var);
 if (num_species_var<0)
  amrex::Error("num species var invalid");

 num_state_base=SpeciesVar;  // den,Temperature
 num_state_material=SpeciesVar;  // den,Temperature
 num_state_material+=num_species_var;

 Vector<Real> elastic_viscosity_temp;
 elastic_viscosity_temp.resize(nmat);
 for (int im=0;im<nmat;im++) 
  elastic_viscosity_temp[im]=0.0;
 pp.queryarr("elastic_viscosity",elastic_viscosity_temp,0,nmat);

 num_materials_viscoelastic=0;
 for (int im=0;im<nmat;im++) {
  if (elastic_viscosity_temp[im]>0.0) {
   num_materials_viscoelastic++;
  } else if (elastic_viscosity_temp[im]==0.0) {
   // do nothing
  } else
   amrex::Error("elastic_viscosity_temp invalid");
 } // im=0..nmat-1 

 Vector<Real> denconst_temp(nmat);
 Vector<Real> den_ceiling_temp(nmat);
 Vector<Real> den_floor_temp(nmat);
 Vector<Real> cavdenconst_temp(nmat);

 Vector<Real> stiffPINFtemp(nmat);
 Vector<Real> stiffCPtemp(nmat);
 Vector<Real> stiffCVtemp(nmat);
 Vector<Real> stiffGAMMAtemp(nmat);

 Vector<Real> DrhoDTtemp(nmat);
 Vector<Real> DrhoDztemp(nmat);
 Vector<Real> tempcutofftemp(nmat);
 Vector<Real> tempcutoffmaxtemp(nmat);
 Vector<Real> tempconst_temp(nmat);
 Vector<Real> initial_temperature_temp(nmat);
 Vector<Real> viscconst_temp(nmat);
 Vector<Real> viscconst_eddy_temp(nmat);
 Vector<int> viscosity_state_model_temp(nmat);
 Vector<Real> heatviscconst_temp(nmat);
 Vector<Real> speciesconst_temp((num_species_var+1)*nmat);
 Vector<Real> speciesviscconst_temp((num_species_var+1)*nmat);
 Vector<int> material_type_temp(nmat);
 Vector<int> FSI_flag_temp(nmat);

 int ZEYU_DCA_SELECT_temp=-1;  // -1=static angle


 pp.getarr("material_type",material_type_temp,0,nmat);

 for (int im=0;im<nmat;im++) {

  stiffPINFtemp[im]=0.0;
  stiffCPtemp[im]=4.1855e+7;
  stiffCVtemp[im]=4.1855e+7;
  stiffGAMMAtemp[im]=0.0;

  DrhoDTtemp[im]=0.0;
  DrhoDztemp[im]=0.0;
  tempcutofftemp[im]=1.0e-8;
  tempcutoffmaxtemp[im]=1.0e+99;
  FSI_flag_temp[im]=0;
 }
 for (int im=0;im<(num_species_var+1)*nmat;im++) {
  speciesviscconst_temp[im]=0.0;
  speciesconst_temp[im]=0.0;
 }

 pp.queryarr("FSI_flag",FSI_flag_temp,0,nmat);

 pp.queryarr("tempcutoff",tempcutofftemp,0,nmat);
 pp.queryarr("tempcutoffmax",tempcutoffmaxtemp,0,nmat);

 pp.getarr("tempconst",tempconst_temp,0,nmat);
 for (int im=0;im<nmat;im++)
  initial_temperature_temp[im]=tempconst_temp[im];
 pp.queryarr("initial_temperature",initial_temperature_temp,0,nmat);

 pp.queryarr("DrhoDT",DrhoDTtemp,0,nmat);
 pp.queryarr("DrhoDz",DrhoDztemp,0,nmat);

 pp.queryarr("stiffPINF",stiffPINFtemp,0,nmat);

 pp.queryarr("stiffCP",stiffCPtemp,0,nmat);
 for (int im=0;im<nmat;im++)
  stiffCVtemp[im]=stiffCPtemp[im];
 pp.queryarr("stiffCV",stiffCVtemp,0,nmat);

 Vector<Real> prerecalesce_stiffCP_temp(nmat);
 for (int im=0;im<nmat;im++)
  prerecalesce_stiffCP_temp[im]=stiffCPtemp[im];
 pp.queryarr("precalesce_stiffCP",prerecalesce_stiffCP_temp,0,nmat);
 Vector<Real> prerecalesce_stiffCV_temp(nmat);
 for (int im=0;im<nmat;im++)
  prerecalesce_stiffCV_temp[im]=stiffCVtemp[im];
 pp.queryarr("precalesce_stiffCV",prerecalesce_stiffCV_temp,0,nmat);

 pp.queryarr("stiffGAMMA",stiffGAMMAtemp,0,nmat);

 pp.getarr("denconst",denconst_temp,0,nmat);

 for (int im=0;im<nmat;im++) {
  cavdenconst_temp[im]=0.0;
  den_ceiling_temp[im]=1.0e+20;
  den_floor_temp[im]=0.0;
 }
 pp.queryarr("cavitation_vapor_density",cavdenconst_temp,0,nmat);
 pp.queryarr("density_floor",den_floor_temp,0,nmat);
 pp.queryarr("density_ceiling",den_ceiling_temp,0,nmat);

 pp.getarr("viscconst",viscconst_temp,0,nmat);
 for (int im=0;im<nmat;im++)
  viscconst_eddy_temp[im]=0.0;
 pp.queryarr("viscconst_eddy",viscconst_eddy_temp,0,nmat);

 Vector<Real> prerecalesce_viscconst_temp(nmat);
 for (int im=0;im<nmat;im++)
  prerecalesce_viscconst_temp[im]=viscconst_temp[im];
 pp.queryarr("precalesce_viscconst",prerecalesce_viscconst_temp,0,nmat);

 for (int im=0;im<nmat;im++)
  viscosity_state_model_temp[im]=0;
 pp.queryarr("viscosity_state_model",
  viscosity_state_model_temp,0,nmat);

 pp.getarr("heatviscconst",heatviscconst_temp,0,nmat);

 Vector<Real> prerecalesce_heatviscconst_temp(nmat);
 for (int im=0;im<nmat;im++)
  prerecalesce_heatviscconst_temp[im]=heatviscconst_temp[im];
 pp.queryarr("precalesce_heatviscconst",prerecalesce_heatviscconst_temp,0,nmat);

 if (num_species_var>0) {
  pp.queryarr("speciesconst",speciesconst_temp,0,num_species_var*nmat);
  pp.queryarr("speciesviscconst",speciesviscconst_temp,0,num_species_var*nmat);
 }

 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;
 Vector<Real> tension_slopetemp(nten);
 Vector<Real> tension_T0temp(nten);
 Vector<Real> tension_mintemp(nten);
 for (int im=0;im<nten;im++) {
  tension_slopetemp[im]=0.0;
  tension_T0temp[im]=293.0;
  tension_mintemp[im]=0.0;
 }
 Vector<Real> tensiontemp(nten);
 Vector<Real> prefreeze_tensiontemp(nten);
 pp.getarr("tension",tensiontemp,0,nten);
 pp.queryarr("tension_slope",tension_slopetemp,0,nten);
 pp.queryarr("tension_T0",tension_T0temp,0,nten);
 pp.queryarr("tension_min",tension_mintemp,0,nten);

 Vector<Real> latent_heat_temp(2*nten);
 Vector<Real> saturation_temp_temp(2*nten);
 for (int i=0;i<nten;i++) { 
  saturation_temp_temp[i]=0.0;
  saturation_temp_temp[i+nten]=0.0;
  latent_heat_temp[i]=0.0;
  latent_heat_temp[i+nten]=0.0;
 }
 pp.queryarr("latent_heat",latent_heat_temp,0,2*nten);
 pp.queryarr("saturation_temp",saturation_temp_temp,0,2*nten);

 Vector<Real> molar_mass_temp(nmat);
 Vector<Real> species_molar_mass_temp(num_species_var+1);
 for (int im=0;im<nmat;im++) {
  molar_mass_temp[im]=1.0;
 }
 for (int im=0;im<num_species_var+1;im++) {
  species_molar_mass_temp[im]=1.0;
 }
 pp.queryarr("molar_mass",molar_mass_temp,0,nmat);

 pp.queryarr("species_molar_mass",
   species_molar_mass_temp,0,num_species_var);

 for (int im=0;im<nten;im++)
  prefreeze_tensiontemp[im]=tensiontemp[im];
 pp.queryarr("prefreeze_tension",prefreeze_tensiontemp,0,nten);


 for (int im=0;im<nmat;im++) {

  if (material_type_temp[im]==999) {

    // non-tessellating cases.
   if ((FSI_flag_temp[im]!=1)&& // prescribed PROB.F90 rigid solid
       (FSI_flag_temp[im]!=2)&& // prescribed sci_clsvof.F90 rigid solid
       (FSI_flag_temp[im]!=4))  // FSI CTML solid
    amrex::Error("FSI_flag_temp invalid");

  } else if (material_type_temp[im]==0) {

   if ((FSI_flag_temp[im]!=0)&& // fluid tessellating
       (FSI_flag_temp[im]!=7)&& // fluid tessellating
       (FSI_flag_temp[im]!=3)&& // ice   tessellating
       (FSI_flag_temp[im]!=6)&& // ice   tessellating
       (FSI_flag_temp[im]!=5))  // FSI PROB.F90 rigid solid. tessellating
    amrex::Error("FSI_flag_temp invalid");

  } else if ((material_type_temp[im]>0)&& 
             (material_type_temp[im]<999)) {

   if ((FSI_flag_temp[im]!=0)&&   // tessallating
       (FSI_flag_temp[im]!=7))    // fluid, tessellating
    amrex::Error("FSI_flag_temp invalid");

  } else {
   amrex::Error("material type invalid");
  }
 }

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid 9");

 int geometry_coord;
 Vector<Real> geometry_prob_lo;
 Vector<Real> geometry_prob_hi;
 Vector<int> geometry_is_periodic;
 int geometry_is_any_periodic;
 read_geometry_raw(geometry_coord,geometry_prob_lo,geometry_prob_hi,
		 geometry_is_periodic,geometry_is_any_periodic);

 int rz_flag=0;
 if ((CoordSys::CoordType) geometry_coord == CoordSys::RZ)  
  rz_flag=1;
 else if ((CoordSys::CoordType) geometry_coord == CoordSys::cartesian)  
  rz_flag=0;
 else if ((CoordSys::CoordType) geometry_coord == CoordSys::CYLINDRICAL)  
  rz_flag=3;
 else
  amrex::Error("CoordSys bust 1");

 Real problox=geometry_prob_lo[0];
 Real probloy=geometry_prob_lo[1];
 Real probloz=geometry_prob_lo[AMREX_SPACEDIM-1];
 Real probhix=geometry_prob_hi[0];
 Real probhiy=geometry_prob_hi[1];
 Real probhiz=geometry_prob_hi[AMREX_SPACEDIM-1];

 int ioproc=0;
 if (ParallelDescriptor::IOProcessor())
  ioproc=1;

 int prescribe_temperature_outflow=0;
 pp.query("prescribe_temperature_outflow",prescribe_temperature_outflow);
 if ((prescribe_temperature_outflow<0)||
     (prescribe_temperature_outflow>3))
  amrex::Error("prescribe_temperature_outflow invalid (fortran_parameters)");

 Real MUSHY_THICK=2.0;
 pp.query("MUSHY_THICK",MUSHY_THICK);

 Real gravity=0.0; 
 int gravity_dir=AMREX_SPACEDIM;
 int invert_gravity=0;
 pp.query("gravity",gravity);
 pp.query("gravity_dir",gravity_dir);
 pp.query("invert_gravity",invert_gravity);
 if ((gravity_dir<1)||(gravity_dir>AMREX_SPACEDIM))
  amrex::Error("gravity dir invalid");

 int n_sites=0;
 pp.query("n_sites",n_sites);
 Real nucleation_init_time=0.0;
 pp.query("nucleation_init_time",nucleation_init_time);

 pp.query("ZEYU_DCA_SELECT",ZEYU_DCA_SELECT_temp);
 if ((ZEYU_DCA_SELECT_temp==-1)||
     ((ZEYU_DCA_SELECT_temp>=1)&&
      (ZEYU_DCA_SELECT_temp<=8))) {
  // do nothing
 } else
  amrex::Error("ZEYU_DCA_SELECT_temp invalid");

 Vector<Real> temp_pos_sites(4);
 for (int dir=0;dir<4;dir++)
  temp_pos_sites[dir]=0.0;
 
 if (n_sites>0) {
  temp_pos_sites.resize(4*n_sites);
  pp.getarr("pos_sites",temp_pos_sites,0,4*n_sites);
 } else if (n_sites==0) {
  // do nothing
 } else
  amrex::Error("n_sites invalid");

 FORT_OVERRIDE(
  &ns_max_level,
  ns_space_blocking_factor.dataPtr(),
  &time_blocking_factor,
  &prescribe_temperature_outflow,
  &rz_flag,
  FSI_flag_temp.dataPtr(),
  &ZEYU_DCA_SELECT_temp,
  &invert_solid_levelset,
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
  &num_species_var,
  &num_materials_viscoelastic,
  &num_state_material,
  &num_state_base,
  &ngeom_raw,
  &ngeom_recon,
  &nmat,
  &num_materials_vel,
  &num_materials_scalar_solve,
  material_type_temp.dataPtr(),
  &nten,
  DrhoDTtemp.dataPtr(),
  DrhoDztemp.dataPtr(),
  tempconst_temp.dataPtr(),
  initial_temperature_temp.dataPtr(),
  tempcutofftemp.dataPtr(),
  tempcutoffmaxtemp.dataPtr(),
  stiffPINFtemp.dataPtr(),
  stiffCPtemp.dataPtr(),
  stiffCVtemp.dataPtr(),
  stiffGAMMAtemp.dataPtr(),
  denconst_temp.dataPtr(),
  den_floor_temp.dataPtr(),
  den_ceiling_temp.dataPtr(),
  cavdenconst_temp.dataPtr(),
  viscconst_temp.dataPtr(),
  viscconst_eddy_temp.dataPtr(),
  viscosity_state_model_temp.dataPtr(),
  elastic_viscosity_temp.dataPtr(),
  heatviscconst_temp.dataPtr(),
  prerecalesce_heatviscconst_temp.dataPtr(),
  prerecalesce_viscconst_temp.dataPtr(),
  prerecalesce_stiffCP_temp.dataPtr(),
  prerecalesce_stiffCV_temp.dataPtr(),
  speciesconst_temp.dataPtr(),
  speciesviscconst_temp.dataPtr(),
  latent_heat_temp.dataPtr(),
  saturation_temp_temp.dataPtr(),
  molar_mass_temp.dataPtr(),
  species_molar_mass_temp.dataPtr(),
  tensiontemp.dataPtr(),
  tension_slopetemp.dataPtr(),
  tension_T0temp.dataPtr(),
  tension_mintemp.dataPtr(),
  prefreeze_tensiontemp.dataPtr(),
  &MUSHY_THICK,
  &gravity,
  &gravity_dir,
  &invert_gravity,
  &fort_stop_time,
  &ioproc);

 int mof_error_ordering_local=0;
 pp.query("mof_error_ordering",mof_error_ordering_local);
 if ((mof_error_ordering_local!=0)&&
     (mof_error_ordering_local!=1))
  amrex::Error("mof_error_ordering_local invalid");
 Vector<int> mof_ordering_local;
 mof_ordering_local.resize(nmat);

 mof_ordering_override(mof_ordering_local,
  nmat,probtype,
  axis_dir,radblob3,
  radblob4,radblob7,
  mof_error_ordering_local,
  FSI_flag_temp);

 pp.queryarr("mof_ordering",mof_ordering_local,0,nmat);
 for (int i=0;i<nmat;i++) {
  if ((mof_ordering_local[i]<0)||
      (mof_ordering_local[i]>nmat+1))
   amrex::Error("mof_ordering_local invalid");
 }

 int temp_POLYGON_LIST_MAX=1000;
 
 FORT_INITMOF(
   mof_ordering_local.dataPtr(),
   &nmat,&MOFITERMAX,
   &MOF_DEBUG_RECON,
   &MOF_TURN_OFF_LS,
   &thread_class::nthreads,
   &temp_POLYGON_LIST_MAX);

 if (ioproc==1) {
  std::cout << "in c++ code, after fort_override\n";
  for (int im=0;im<nmat;im++) {
   std::cout << "im= " << im << " mof_ordering_local= " <<
    mof_ordering_local[im] << '\n';
  }
 }
}  // subroutine fortran_parameters()




void
NavierStokes::variableCleanUp ()
{
    desc_lst.clear();
}

void
NavierStokes::read_geometry ()
{
    //
    // Must load coord here because CoordSys hasn't read it in yet.
    //
    int geometry_coord;
    Vector<Real> geometry_prob_lo;
    Vector<Real> geometry_prob_hi;
    Vector<int> geometry_is_periodic;
    int geometry_is_any_periodic;

    read_geometry_raw(geometry_coord,geometry_prob_lo,geometry_prob_hi,
		    geometry_is_periodic,geometry_is_any_periodic);

    if ((CoordSys::CoordType) geometry_coord == CoordSys::RZ)  
     if (AMREX_SPACEDIM==3)
      amrex::Error("No RZ in 3d");

    if (((CoordSys::CoordType) geometry_coord == CoordSys::RZ) && 
        (phys_bc.lo(0) != Symmetry)) {
        phys_bc.setLo(0,Symmetry);
        temperature_phys_bc.setLo(0,Symmetry);
        species_phys_bc.setLo(0,Symmetry);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "\n WARNING: Setting phys_bc at xlo to Symmetry\n\n";
    }


} // end subroutine read_geometry

void
NavierStokes::read_params ()
{
    ParmParse pp("ns");

    pp.query("check_nan",check_nan);
    pp.query("v",verbose);
    pp.query("fab_verbose",fab_verbose);
    pp.query("output_drop_distribution",output_drop_distribution);
    pp.query("extend_pressure_into_solid",extend_pressure_into_solid);
    pp.query("show_timings",show_timings);
    pp.query("show_mem",show_mem);

    pp.query("slice_dir",slice_dir);
    xslice.resize(AMREX_SPACEDIM);
    if ((slice_dir>=0)&&(slice_dir<AMREX_SPACEDIM)) {
     for (int i=0;i<AMREX_SPACEDIM;i++)
      xslice[i]=0.0;
     pp.queryarr("xslice",xslice,0,AMREX_SPACEDIM);
    } else
     amrex::Error("slice_dir invalid");


    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "check_nan " << check_nan << '\n';
     std::cout << "NavierStokes.verbose " << verbose << '\n';
     std::cout << "slice_dir " << slice_dir << '\n';
     for (int i=0;i<AMREX_SPACEDIM;i++) {
      std::cout << "i=" << i << '\n';
      std::cout << "xslice " << xslice[i] << '\n';
     }
    } 

    pp.query("nblocks",nblocks);

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
    num_materials_vel=1;
    num_materials_scalar_solve=1;
    pp.get("num_materials",num_materials);
    if ((num_materials<2)||(num_materials>999))
     amrex::Error("num materials invalid");

    int nmat=num_materials;
    int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

    pp.query("num_materials_vel",num_materials_vel);
    if (num_materials_vel!=1)
     amrex::Error("num_materials_vel==1 required");

    pp.query("num_materials_scalar_solve",num_materials_scalar_solve);
    if ((num_materials_scalar_solve!=1)&&
        (num_materials_scalar_solve!=nmat))
     amrex::Error("num_materials_scalar_solve invalid");

    pp.query("use_supermesh",use_supermesh);
    if ((use_supermesh!=0)&&
        (use_supermesh!=1))
     amrex::Error("use_supermesh invalid");

    pp.query("ncomp_sum_int_user",ncomp_sum_int_user);
    if (ncomp_sum_int_user>=0) {
     // do nothing
    } else
     amrex::Error("ncomp_sum_int_user invalid");

    // blob_matrix,blob_RHS,blob_velocity,
    // blob_integral_momentum,blob_energy,
    // blob_mass_for_velocity (3 components)
    // blob_volume, 
    // blob_center_integral,blob_center_actual
    // blob_perim, blob_perim_mat, blob_triple_perim, 
    num_elements_blobclass=
     3*(2*AMREX_SPACEDIM)*(2*AMREX_SPACEDIM)+
     3*(2*AMREX_SPACEDIM)+
     3*(2*AMREX_SPACEDIM)+
     2*(2*AMREX_SPACEDIM)+
     1+
     3+1+
     2*AMREX_SPACEDIM+1+nmat+nmat*nmat;

    int ns_max_level;
    Vector<int> ns_max_grid_size;
    int the_max_grid_size=0;
    int cnt_max_grid_size;

    ParmParse ppamr("amr");
    ppamr.get("max_level",ns_max_level);
    Vector<int> ns_n_error_buf;
    ns_n_error_buf.resize(ns_max_level);
    for (int ilev=0;ilev<ns_max_level;ilev++) 
     ns_n_error_buf[ilev]=1;
    ppamr.queryarr("n_error_buf",ns_n_error_buf,0,ns_max_level);

    cnt_max_grid_size=ppamr.countval("max_grid_size");
  
    if (cnt_max_grid_size==0) {
     ns_max_grid_size.resize(1);
     ns_max_grid_size[0]=0;
    } else if (cnt_max_grid_size==1) {
     ppamr.get("max_grid_size",the_max_grid_size);
     ns_max_grid_size.resize(1);
     ns_max_grid_size[0]=the_max_grid_size;
    } else if (cnt_max_grid_size>1) {
     ns_max_grid_size.resize(cnt_max_grid_size);
     ppamr.getarr("max_grid_size",ns_max_grid_size,0,cnt_max_grid_size);
    } else
     amrex::Error("cnt_max_grid_size invalid");

    int def_n_proper=1;
    ppamr.query("n_proper",def_n_proper);

    int local_plotfile_on_restart=0;
    ppamr.query("plotfile_on_restart",local_plotfile_on_restart);
    int local_checkpoint_on_restart=0;
    ppamr.query("checkpoint_on_restart",local_checkpoint_on_restart);
    int local_regrid_on_restart=0;
    ppamr.query("regrid_on_restart",local_regrid_on_restart);

    tecplot_max_level=ns_max_level;
    max_level_two_materials=ns_max_level;
    pp.query("tecplot_max_level",tecplot_max_level);
    pp.query("max_level_two_materials",max_level_two_materials);

    radius_cutoff.resize(nmat);
    for (int i=0;i<nmat;i++)
     radius_cutoff[i]=0;

     // default=0.
     // 0=>never adapt  -1=>tag cell for AMR if owned by material in question.
     // otherwise, if radius<radius_cutoff * dx then adapt.
     // fluid-fluid interfaces are always adapted.
    pp.queryarr("radius_cutoff",radius_cutoff,0,nmat);
    for (int i=0;i<nmat;i++)
     if (radius_cutoff[i]<-1)
      amrex::Error("radius_cutoff invalid");

    if ((tecplot_max_level<0)||
        (tecplot_max_level>ns_max_level))
     amrex::Error("tecplot_max_level invalid"); 

    if ((max_level_two_materials<0)||
        (max_level_two_materials>ns_max_level))
     amrex::Error("max_level_two_materials invalid"); 

    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "use_supermesh " << 
       use_supermesh << '\n';
     std::cout << "ncomp_sum_int_user " << 
       ncomp_sum_int_user << '\n';
     std::cout << "tecplot_max_level " << 
       tecplot_max_level << '\n';
     std::cout << "max_level_two_materials " << 
       max_level_two_materials << '\n';
     for (int i=0;i<nmat;i++)
      std::cout << "im=" << i << " radius_cutoff= " << 
        radius_cutoff[i] << '\n';
    }
    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "num_elements_blobclass= " << 
      num_elements_blobclass << '\n';
    }

    pp.query("ncoarseblocks",ncoarseblocks);
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
    Vector<int> temperature_lo_bc(AMREX_SPACEDIM);
    Vector<int> temperature_hi_bc(AMREX_SPACEDIM);
    Vector<int> species_lo_bc(AMREX_SPACEDIM);
    Vector<int> species_hi_bc(AMREX_SPACEDIM);
    pp.getarr("lo_bc",lo_bc,0,AMREX_SPACEDIM);
    pp.getarr("hi_bc",hi_bc,0,AMREX_SPACEDIM);
    for (int i = 0; i < AMREX_SPACEDIM; i++) {
     phys_bc.setLo(i,lo_bc[i]);
     phys_bc.setHi(i,hi_bc[i]);
     temperature_phys_bc.setLo(i,lo_bc[i]);
     temperature_phys_bc.setHi(i,hi_bc[i]);
     temperature_lo_bc[i]=lo_bc[i];
     temperature_hi_bc[i]=hi_bc[i];
     species_phys_bc.setLo(i,lo_bc[i]);
     species_phys_bc.setHi(i,hi_bc[i]);
     species_lo_bc[i]=lo_bc[i];
     species_hi_bc[i]=hi_bc[i];
    }
    pp.queryarr("temperature_lo_bc",temperature_lo_bc,0,AMREX_SPACEDIM);
    pp.queryarr("temperature_hi_bc",temperature_hi_bc,0,AMREX_SPACEDIM);
    pp.queryarr("species_lo_bc",species_lo_bc,0,AMREX_SPACEDIM);
    pp.queryarr("species_hi_bc",species_hi_bc,0,AMREX_SPACEDIM);
    for (int i = 0; i < AMREX_SPACEDIM; i++) {
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

    int geometry_is_all_periodic=1;

    if (geometry_is_any_periodic==1) {
     for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
      if (geometry_is_periodic[dir]==1) {
       if ((lo_bc[dir] != Interior)||
           (temperature_lo_bc[dir]!=Interior)) {
        std::cerr << "NavierStokes::variableSetUp:periodic in direction "
            << dir << " but low BC is not Interior\n";
        amrex::Abort("NavierStokes::read_params()");
       }
       if ((hi_bc[dir] != Interior)||
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

    FORT_SET_PERIODIC_VAR(geometry_is_periodic.dataPtr());

    for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
     if (geometry_is_periodic[dir]==0) {
      if ((lo_bc[dir] == Interior)||
          (temperature_lo_bc[dir]==Interior)) {
       std::cerr << "NavierStokes::variableSetUp:Interior bc in direction "
                 << dir << " but not defined as periodic\n";
       amrex::Abort("NavierStokes::read_params()");
      }
      if ((hi_bc[dir] == Interior)||
          (temperature_hi_bc[dir]==Interior)) {
       std::cerr << "NavierStokes::variableSetUp:Interior bc in direction "
                 << dir << " but not defined as periodic\n";
       amrex::Abort("NavierStokes::read_params()");
      }
     }
    } // dir=0.. sdim-1

    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "phys_bc= " << phys_bc << '\n';
     std::cout << "temperature_phys_bc= " << temperature_phys_bc << '\n';
     std::cout << "species_phys_bc= " << species_phys_bc << '\n';
    }

    //
    // Get timestepping parameters.
    //
    pp.get("cfl",cfl);
    if (cfl>0.9)
     amrex::Warning("WARNING: cfl should be less than or equal to 0.9");

    pp.query("enable_spectral",enable_spectral);

    if (enable_spectral==0) {
     viscous_enable_spectral=0;
    } else if (enable_spectral==1) {
     viscous_enable_spectral=0;
    } else if (enable_spectral==2) {
     viscous_enable_spectral=0;
    } else if (enable_spectral==3) {
     viscous_enable_spectral=0;
    } else
     amrex::Error("enable_spectral invalid");
	
    projection_enable_spectral=enable_spectral;

    pp.query("viscous_enable_spectral",viscous_enable_spectral);
    pp.query("projection_enable_spectral",projection_enable_spectral);

    if ((viscous_enable_spectral==1)||  //SEM space and time
        (viscous_enable_spectral==2)) { //SEM space

     if (geometry_is_all_periodic==1) {
      // do nothing
     } else if (geometry_is_all_periodic==0) {
      amrex::Error("no slip BC not implemented for space order>2");
     } else
      amrex::Error("geometry_is_all_periodic invalid");

    } else if ((viscous_enable_spectral==0)||  //2nd space, 1st time
	       (viscous_enable_spectral==3)) { // 2nd order space, SEM time
     // do nothing
    } else
     amrex::Error("viscous_enable_spectral invalid");

    pp.query("SEM_upwind",SEM_upwind);
    if (SEM_upwind==1) {
     // do nothing
    } else if (SEM_upwind==0) {
     // do nothing
    } else
     amrex::Error("SEM_upwind invalid");

    pp.query("SEM_advection_algorithm",SEM_advection_algorithm);
    if (SEM_advection_algorithm==1) {
     // do nothing
    } else if (SEM_advection_algorithm==0) {
     // do nothing
    } else
     amrex::Error("SEM_advection_algorithm invalid");

    pp.query("continuous_mof",continuous_mof);
    pp.query("VOF_reflux",VOF_reflux);

    pp.query("init_shrink",init_shrink);
    pp.query("dt_max",dt_max);

    pp.query("change_max",change_max);
    change_max_init=change_max;
    pp.query("change_max_init",change_max_init);

    pp.query("fixed_dt",fixed_dt);
    fixed_dt_init=fixed_dt;
    pp.query("fixed_dt_init",fixed_dt_init);

    pp.query("min_velocity_for_dt",min_velocity_for_dt);
    if (min_velocity_for_dt<0.0)
     amrex::Error("min_velocity_for_dt invalid");

    pp.query("fixed_dt_velocity",fixed_dt_velocity);
    pp.query("sum_interval",sum_interval);

    pp.query("profile_debug",profile_debug);
    if ((profile_debug!=0)&&(profile_debug!=1))
     amrex::Error("profile_debug invalid");

    pp.query("ns_tiling",ns_tiling);
    if (ns_tiling==true) {
     // do nothing
    } else if (ns_tiling==false) {
     // do nothing
    } else
     amrex::Error("ns_tiling invalid");

    if (thread_class::nthreads<1)
     amrex::Error("thread_class::nthreads invalid ns init");


    MUSHY_THICK=2.0;
    pp.query("MUSHY_THICK",MUSHY_THICK);

    pp.query("gravity",gravity);
    pp.query("gravity_dir",gravity_dir);
    pp.query("terminal_velocity_dt",terminal_velocity_dt);
    pp.query("invert_gravity",invert_gravity);
    if ((gravity_dir<1)||(gravity_dir>AMREX_SPACEDIM))
     amrex::Error("gravity dir invalid");
    if ((terminal_velocity_dt<0)||
        (terminal_velocity_dt>1))
     amrex::Error("terminal_velocity_dt invalid");

    pp.get("visc_coef",visc_coef);

    pp.query("include_viscous_heating",include_viscous_heating);
    if ((include_viscous_heating<0)||
        (include_viscous_heating>1))
     amrex::Error("include_viscous_heating invalid");
   
    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "local_plotfile_on_restart (NS) = " << 
	     local_plotfile_on_restart << '\n';
     std::cout << "local_checkpoint_on_restart (NS) = " << 
	     local_checkpoint_on_restart << '\n';
     std::cout << "local_regrid_on_restart (NS) = " << 
	     local_regrid_on_restart << '\n';
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
     std::cout << "extend_pressure_into_solid " << 
      extend_pressure_into_solid << '\n';
     std::cout << "visc_coef " << visc_coef << '\n';
     std::cout << "include_viscous_heating " << include_viscous_heating << '\n';

     std::cout << "MUSHY_THICK " << MUSHY_THICK << '\n';

     std::cout << "gravity " << gravity << '\n';
     std::cout << "invert_gravity " << invert_gravity << '\n';
     std::cout << "gravity_dir " << gravity_dir << '\n';
     std::cout << "terminal_velocity_dt " << terminal_velocity_dt << '\n';
     std::cout << "cfl " << cfl << '\n';
     std::cout << "enable_spectral " << enable_spectral << '\n';
     std::cout << "viscous_enable_spectral " << 
       viscous_enable_spectral << '\n';
     std::cout << "projection_enable_spectral " << 
       projection_enable_spectral << '\n';
     std::cout << "SEM_upwind " << SEM_upwind << '\n';
     std::cout << "SEM_advection_algorithm " << 
       SEM_advection_algorithm << '\n';
     std::cout << "continuous_mof " << continuous_mof << '\n';
     std::cout << "VOF_reflux " << VOF_reflux << '\n';
    }

    pp.query("FD_curv_interp",FD_curv_interp);
    if ((FD_curv_interp!=0)&&
        (FD_curv_interp!=1))
     amrex::Error("FD_curv_interp invalid");

    custom_nucleation_model=0;
    pp.query("custom_nucleation_model",custom_nucleation_model);
    if ((custom_nucleation_model!=0)&&
        (custom_nucleation_model!=1))
     amrex::Error("custom_nucleation_model invalid");

    n_sites=0;
    pp.query("n_sites",n_sites);
    pos_sites.resize(4);
    for (int dir=0;dir<4;dir++)
     pos_sites[dir]=0.0;
    if (n_sites>0) {
     pos_sites.resize(4*n_sites);
     pp.getarr("pos_sites",pos_sites,0,4*n_sites);
    } else if (n_sites==0) {
     // do nothing
    } else
     amrex::Error("n_sites invalid");
   
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
    pp.queryarr("visual_ncell",visual_ncell,0,AMREX_SPACEDIM);
    pp.query("visual_compare",visual_compare);

    pp.query("visual_tessellate_vfrac",visual_tessellate_vfrac);
    pp.query("visual_revolve",visual_revolve);
    pp.query("visual_option",visual_option);

    if ((visual_tessellate_vfrac!=0)&&
        (visual_tessellate_vfrac!=1))
     amrex::Error("visual_tessellate_vfrac invalid");

    if (visual_revolve<0) {
     std::cout << "visual_revolve: " << visual_revolve << '\n';
     amrex::Error("visual_revolve invalid");
    }
    if ((visual_option<-2)||(visual_option>-1))
     amrex::Error("visual_option invalid try -1 or -2");

    if (AMREX_SPACEDIM==2) {
     adapt_quad_depth=2;
     if ((probtype==28)||(probtype==29)||(probtype==31))
      adapt_quad_depth=5; 
    } else if (AMREX_SPACEDIM==3) {
     adapt_quad_depth=1;
     if ((probtype==28)||(probtype==29)||(probtype==31))
      adapt_quad_depth=3;
    } else
     amrex::Error("dimension bust");

    pp.query("adapt_quad_depth",adapt_quad_depth);

    pp.query("invert_solid_levelset",invert_solid_levelset);
    if (!((invert_solid_levelset==1)||(invert_solid_levelset==0)))
     amrex::Error("invert_solid_levelset invalid");

    pp.query("law_of_the_wall",law_of_the_wall);
    if ((law_of_the_wall==0)||
        (law_of_the_wall==1)||
	(law_of_the_wall==2)) {
     // do nothing
    } else {
     amrex::Error("law_of_the_wall invalid");
    }


    pp.query("ZEYU_DCA_SELECT",ZEYU_DCA_SELECT);
    if ((ZEYU_DCA_SELECT==-1)||
        ((ZEYU_DCA_SELECT>=1)&&
	 (ZEYU_DCA_SELECT<=8))) {
     // do nothing
    } else
     amrex::Error("ZEYU_DCA_SELECT invalid");

    FSI_flag.resize(nmat);
    FSI_refine_factor.resize(nmat);
    FSI_bounding_box_ngrow.resize(nmat);

    material_type.resize(nmat);
    pp.getarr("material_type",material_type,0,nmat);
 
    for (int i=0;i<nmat;i++) {
     FSI_flag[i]=0;
     FSI_refine_factor[i]=1;
     FSI_bounding_box_ngrow[i]=3;
    }
    pp.queryarr("FSI_flag",FSI_flag,0,nmat);
    pp.queryarr("FSI_refine_factor",FSI_refine_factor,0,nmat);
    pp.queryarr("FSI_bounding_box_ngrow",FSI_bounding_box_ngrow,0,nmat);

    pp.query("CTML_FSI_numsolids",CTML_FSI_numsolids);
    pp.query("CTML_force_model",CTML_force_model);
    if ((CTML_force_model!=0)&&(CTML_force_model!=1))
     amrex::Error("CTML_force_model invalid");

    int nparts=0;
    int CTML_FSI_numsolids_test=0;
    for (int im=0;im<nmat;im++) {

     if (FSI_flag[im]==4) // non-tessellating
      CTML_FSI_numsolids_test++;
 
     if (ns_is_lag_part(im)==1) {
      nparts++;
     } else if (ns_is_lag_part(im)==0) {
      // do nothing
     } else
      amrex::Error("ns_is_lag_part invalid");
    }  // im=0..nmat-1
    im_solid_map.resize(nparts);

    if (CTML_FSI_numsolids!=CTML_FSI_numsolids_test)
     amrex::Error("CTML_FSI_numsolids!=CTML_FSI_numsolids_test");

    nparts=0;
    for (int im=0;im<nmat;im++) {
     if (ns_is_lag_part(im)==1) {
      im_solid_map[nparts]=im;
      nparts++;
     } else if (ns_is_lag_part(im)==0) {
      // do nothing
     } else
      amrex::Error("ns_is_lag_part invalid");
    }  // im=0..nmat-1
    if (nparts!=im_solid_map.size())
     amrex::Error("nparts!=im_solid_map.size()");

    elastic_viscosity.resize(nmat);
    for (int im=0;im<nmat;im++) {
     elastic_viscosity[im]=0.0;
    }
    pp.queryarr("elastic_viscosity",elastic_viscosity,0,nmat);

    num_materials_viscoelastic=0;
    for (int i=0;i<nmat;i++) {
     if (elastic_viscosity[i]>0.0) {
      im_elastic_map.resize(num_materials_viscoelastic+1);
      im_elastic_map[num_materials_viscoelastic]=i;
      num_materials_viscoelastic++;
     } else if (elastic_viscosity[i]==0.0) {
      // do nothing
     } else
      amrex::Error("elastic_viscosity invalid");
    } // im=0..nmat-1 

    NUM_STATE_TYPE=DIV_Type+1;

    if (nparts==0) {
     Solid_State_Type=-1;
    } else if ((nparts>=1)&&(nparts<=nmat)) {
     Solid_State_Type=NUM_STATE_TYPE;
     NUM_STATE_TYPE++;
    } else
     amrex::Error("nparts invalid");
  
    if ((num_materials_viscoelastic>=0)&&
        (num_materials_viscoelastic<=nmat)) {
     Tensor_Type=NUM_STATE_TYPE;
     NUM_STATE_TYPE++;
    } else
     amrex::Error("num_materials_viscoelastic invalid");

    for (int i=0;i<nmat;i++) {
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

    } // i=0..nmat-1

    pp.query("gmres_precond_iter_base",gmres_precond_iter_base);

    ParmParse pplp("Lp");
    pplp.query("smooth_type",smooth_type);
    pplp.query("bottom_smooth_type",bottom_smooth_type);
    if (smooth_type!=2)
     amrex::Warning("WARNING: ILU smoother is best");
    pplp.query("use_mg_precond_in_mglib",use_mg_precond_in_mglib);

    pplp.query("bottom_bottom_tol_factor",bottom_bottom_tol_factor);

    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "smooth_type " << smooth_type << '\n';
     std::cout << "bottom_smooth_type " << bottom_smooth_type << '\n';
     std::cout << "use_mg_precond_in_mglib " <<use_mg_precond_in_mglib<<'\n';

     std::cout << "bottom_bottom_tol_factor " <<
       bottom_bottom_tol_factor<<'\n';
     for (int i=0;i<nmat;i++) {
      std::cout << "i= " << i << " FSI_flag " << FSI_flag[i] << '\n';
      std::cout << "i= " << i << " FSI_refine_factor " << 
	     FSI_refine_factor[i] << '\n';
      std::cout << "i= " << i << " FSI_bounding_box_ngrow " << 
	     FSI_bounding_box_ngrow[i] << '\n';
     }
     std::cout << "CTML_FSI_numsolids " << CTML_FSI_numsolids << '\n';

     std::cout << "invert_solid_levelset " << invert_solid_levelset << '\n';
     std::cout << "law_of_the_wall " << law_of_the_wall << '\n';
     std::cout << "ZEYU_DCA_SELECT " << ZEYU_DCA_SELECT << '\n';
     std::cout << "adapt_quad_depth= " << adapt_quad_depth << '\n';
    }

    pp.get("adv_dir",adv_dir);
    pp.get("adv_vel",adv_vel);
    pp.get("rgasinlet",rgasinlet);
    pp.query("slipcoeff",slipcoeff);
    pp.get("vinletgas",vinletgas);

    pp.get("twall",twall);

    pp.query("curv_stencil_height",curv_stencil_height);
    if (curv_stencil_height!=4)
     amrex::Error("must have curv_stencil_height==4");


    pp.query("bicgstab_max_num_outer_iter",bicgstab_max_num_outer_iter);
    pp.query("slope_limiter_option",slope_limiter_option);

    pp.query("EILE_flag",EILE_flag);

    if ((EILE_flag==0)|| // Sussman and Puckett
        (EILE_flag==1)|| // EILE
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

    massface_index=facespecies_index+num_species_var;
    vofface_index=massface_index+2*num_materials;
    ncphys=vofface_index+2*num_materials;

    MOFITERMAX=DEFAULT_MOFITERMAX;
    pp.query("MOFITERMAX",MOFITERMAX);
    if ((MOFITERMAX<0)||(MOFITERMAX>50))
     amrex::Error("mof iter max invalid in navierstokes");

    MOF_TURN_OFF_LS=0;  
    pp.query("MOF_TURN_OFF_LS",MOF_TURN_OFF_LS);
    if ((MOF_TURN_OFF_LS!=0)&&(MOF_TURN_OFF_LS!=1))
     amrex::Error("mof turn off ls invalid in navierstokes");

    MOF_DEBUG_RECON=0; 
    pp.query("MOF_DEBUG_RECON",MOF_DEBUG_RECON);
    if ((MOF_DEBUG_RECON!=0)&&(MOF_DEBUG_RECON!=1)&&
        (MOF_DEBUG_RECON!=2))
     amrex::Error("mof debug recon invalid in navierstokes");

    num_state_base=SpeciesVar;  // den,Temperature
    num_state_material=SpeciesVar;  // den,Temperature
    num_state_material+=num_species_var;

    pp.query("make_interface_incomp",make_interface_incomp);
    if ((make_interface_incomp!=0)&&
        (make_interface_incomp!=1)&&
        (make_interface_incomp!=2))
     amrex::Error("make_interface_incomp invalid");

    advection_order.resize(nmat);
    density_advection_order.resize(nmat);
    for (int i=0;i<nmat;i++) {
     advection_order[i]=1;
     density_advection_order[i]=1;
    }
    pp.queryarr("advection_order",advection_order,0,nmat);
    for (int i=0;i<nmat;i++) {
     density_advection_order[i]=advection_order[i];
     if (!((advection_order[i]==1)||
           (advection_order[i]==2)))
      amrex::Error("advection_order invalid");
    }
    pp.queryarr("density_advection_order",density_advection_order,0,nmat);
    for (int i=0;i<nmat;i++) {
     if (!((density_advection_order[i]==1)||
           (density_advection_order[i]==2)))
      amrex::Error("density_advection_order invalid");
    }

    stiffPINF.resize(nmat);
    stiffCP.resize(nmat);
    stiffCV.resize(nmat);
    stiffGAMMA.resize(nmat);

    DrhoDT.resize(nmat);
    DrhoDz.resize(nmat);
    override_density.resize(nmat);

    temperature_source_cen.resize(AMREX_SPACEDIM);
    temperature_source_rad.resize(AMREX_SPACEDIM);

    tempconst.resize(nmat);
    initial_temperature.resize(nmat);
    tempcutoff.resize(nmat);
    tempcutoffmax.resize(nmat);
    viscconst.resize(nmat);
    viscconst_eddy.resize(nmat);
    viscosity_state_model.resize(nmat);
    viscoelastic_model.resize(nmat);
    les_model.resize(nmat);
    viscconst_interface.resize(nten);
    speciesconst.resize((num_species_var+1)*nmat);
    speciesviscconst.resize((num_species_var+1)*nmat);
    speciesviscconst_interface.resize((num_species_var+1)*nten);
    species_molar_mass.resize(num_species_var+1);

    for (int j=0;j<=num_species_var;j++)
     species_molar_mass[j]=1.0;

    for (int i=0;i<nten;i++) {
     viscconst_interface[i]=0.0;
     for (int j=0;j<=num_species_var;j++)
      speciesviscconst_interface[j*nten+i]=0.0;
    }
    pp.queryarr("viscconst_interface",viscconst_interface,0,nten);
    if (num_species_var>0) {
     pp.queryarr("speciesviscconst_interface",
      speciesviscconst_interface,0,num_species_var*nten);

     pp.queryarr("species_molar_mass",
      species_molar_mass,0,num_species_var);
    }
     // in: read_params

    species_evaporation_density.resize(num_species_var+1);
    for (int i=0;i<num_species_var+1;i++)
     species_evaporation_density[i]=0.0;

    if (num_species_var>0)
     pp.queryarr("species_evaporation_density",species_evaporation_density,
      0,num_species_var);

    spec_material_id_LIQUID.resize(num_species_var+1);
    spec_material_id_AMBIENT.resize(num_species_var+1);
    for (int i=0;i<num_species_var+1;i++) {
     spec_material_id_LIQUID[i]=0;
     spec_material_id_AMBIENT[i]=0;
    }

    if (num_species_var>0) {
     pp.queryarr("spec_material_id_LIQUID",spec_material_id_LIQUID,
       0,num_species_var);
     pp.queryarr("spec_material_id_AMBIENT",spec_material_id_AMBIENT,
       0,num_species_var);
    }
    
    vorterr.resize(nmat);
    pressure_error_cutoff.resize(nmat);
    temperature_error_cutoff.resize(nmat);

    recalesce_model_parameters.resize(3*nmat);

    microlayer_substrate.resize(nmat);
    microlayer_angle.resize(nmat);
    microlayer_size.resize(nmat);
    macrolayer_size.resize(nmat);
    max_contact_line_size.resize(nmat);
 
    microlayer_temperature_substrate.resize(nmat);

     // in: read_params
     
    cavitation_pressure.resize(nmat);
    cavitation_vapor_density.resize(nmat);
    cavitation_tension.resize(nmat);
    for (int i=0;i<nmat;i++) {
     cavitation_pressure[i]=0.0; 
     cavitation_vapor_density[i]=0.0; 
     cavitation_tension[i]=0.0; 
    }
    pp.queryarr("cavitation_pressure",cavitation_pressure,0,nmat);
    pp.queryarr("cavitation_vapor_density",cavitation_vapor_density,0,nmat);
    pp.queryarr("cavitation_tension",cavitation_tension,0,nmat);
 
     // in: read_params

    hardwire_Y_gamma.resize(2*nten);
    hardwire_T_gamma.resize(2*nten);
    saturation_temp.resize(2*nten);
    saturation_temp_curv.resize(2*nten);
    saturation_temp_vel.resize(2*nten);
    saturation_temp_min.resize(2*nten);
    saturation_temp_max.resize(2*nten);
    nucleation_temp.resize(2*nten);
    nucleation_pressure.resize(2*nten);
    nucleation_pmg.resize(2*nten);
    nucleation_mach.resize(2*nten);
    latent_heat.resize(2*nten);
    reaction_rate.resize(2*nten);
    freezing_model.resize(2*nten);
    Tanasawa_or_Schrage.resize(2*nten);
    mass_fraction_id.resize(2*nten);
    distribute_from_target.resize(2*nten);
    tension.resize(nten);
    tension_slope.resize(nten);
    tension_T0.resize(nten);
    tension_min.resize(nten);
    prefreeze_tension.resize(nten);

     // (dir,side)  (1,1),(2,1),(3,1),(1,2),(2,2),(3,2)
    outflow_velocity_buffer_size.resize(2*AMREX_SPACEDIM);

    cap_wave_speed.resize(nten);

    prerecalesce_stiffCP.resize(nmat);
    prerecalesce_stiffCV.resize(nmat);
    prerecalesce_viscconst.resize(nmat);
    prerecalesce_heatviscconst.resize(nmat);

    for (int i=0;i<3*nmat;i++) { 
     recalesce_model_parameters[i]=0.0;
    }

    nucleation_period=0.0;
    nucleation_init_time=0.0;

    for (int i=0;i<nmat;i++) {
     microlayer_substrate[i]=0;
     microlayer_angle[i]=0.0;
     microlayer_size[i]=0.0;
     macrolayer_size[i]=0.0;
     max_contact_line_size[i]=0.0;
     microlayer_temperature_substrate[i]=0.0;
    }

    for (int i=0;i<nten;i++) { 
     hardwire_Y_gamma[i]=0.0;
     hardwire_Y_gamma[i+nten]=0.0;
     hardwire_T_gamma[i]=0.0;
     hardwire_T_gamma[i+nten]=0.0;
     saturation_temp[i]=0.0;
     saturation_temp[i+nten]=0.0;
     saturation_temp_curv[i]=0.0;
     saturation_temp_curv[i+nten]=0.0;
     saturation_temp_vel[i]=0.0;
     saturation_temp_vel[i+nten]=0.0;
     saturation_temp_min[i]=0.0;
     saturation_temp_min[i+nten]=0.0;
     saturation_temp_max[i]=1.0e+20;
     saturation_temp_max[i+nten]=1.0e+20;
     nucleation_temp[i]=0.0;
     nucleation_temp[i+nten]=0.0;
     nucleation_pressure[i]=0.0;
     nucleation_pmg[i]=0.0;
     nucleation_mach[i]=0.0;
     nucleation_pressure[i+nten]=0.0;
     nucleation_pmg[i+nten]=0.0;
     nucleation_mach[i+nten]=0.0;
     latent_heat[i]=0.0;
     reaction_rate[i]=0.0;
     latent_heat[i+nten]=0.0;
     reaction_rate[i+nten]=0.0;
     freezing_model[i]=0;
     freezing_model[i+nten]=0;
     Tanasawa_or_Schrage[i]=0;
     Tanasawa_or_Schrage[i+nten]=0;
     mass_fraction_id[i]=0;
     mass_fraction_id[i+nten]=0;
     distribute_from_target[i]=0;
     distribute_from_target[i+nten]=0;
    } // i=0..nten-1

    density_floor.resize(nmat);
    density_floor_expansion.resize(nmat);
    for (int i=0;i<nmat;i++)
     density_floor[i]=0.0;
    pp.queryarr("density_floor",density_floor,0,nmat);
    density_ceiling.resize(nmat);
    density_ceiling_expansion.resize(nmat);
    molar_mass.resize(nmat);
    for (int i=0;i<nmat;i++)
     density_ceiling[i]=1.0e+20;
    pp.queryarr("density_ceiling",density_ceiling,0,nmat);

    denconst.resize(nmat);
    pp.getarr("denconst",denconst,0,nmat);

    denconst_min=denconst[0];
    denconst_max=denconst[0];

    for (int i=0;i<nmat;i++) {
     if (denconst[i]<denconst_min)
      denconst_min=denconst[i];
     if (denconst[i]>denconst_max)
      denconst_max=denconst[i];
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
     density_ceiling_expansion[i]=denconst[i];
     density_floor_expansion[i]=denconst[i];
     molar_mass[i]=1.0;
    } // i=0..nmat-1

    pp.queryarr("density_floor_expansion",density_floor_expansion,0,nmat);
    pp.queryarr("density_ceiling_expansion",density_ceiling_expansion,0,nmat);
    pp.queryarr("molar_mass",molar_mass,0,nmat);

    for (int i=0;i<nmat;i++) {
     if (density_ceiling_expansion[i]<=0.0) {
      amrex::Error("density_ceiling_expansion[i]<=0.0");
     } else if (density_ceiling_expansion[i]<denconst[i]) {
      amrex::Error("density_ceiling_expansion[i]<denconst[i]");
     }
     if (density_floor_expansion[i]<=0.0) {
      amrex::Error("density_floor_expansion[i]<=0.0");
     } else if (density_floor_expansion[i]>denconst[i]) {
      amrex::Error("density_floor_expansion[i]>denconst[i]");
     }
    } // i=0..nmat-1


    denconst_interface.resize(nten);
    for (int i=0;i<nten;i++) 
     denconst_interface[i]=0.0;
    pp.queryarr("denconst_interface",denconst_interface,0,nten);

    denconst_gravity.resize(nmat);
    for (int i=0;i<nmat;i++) 
     denconst_gravity[i]=1.0;
    pp.queryarr("denconst_gravity",denconst_gravity,0,nmat);

    pp.query("stokes_flow",stokes_flow);
    pp.query("cancel_advection",cancel_advection);

    added_weight.resize(nmat);
    for (int i=0;i<nmat;i++) 
     added_weight[i]=1.0;
    pp.queryarr("added_weight",added_weight,0,nmat);

    for (int i=0;i<(num_species_var+1)*nmat;i++) {
     speciesconst[i]=0.0;
     speciesviscconst[i]=0.0;
    }

    for (int i=0;i<nmat;i++) {

     stiffPINF[i]=0.0;
     stiffCP[i]=4.1855e+7;
     stiffCV[i]=4.1855e+7;
     stiffGAMMA[i]=0.0;

     tempcutoff[i]=1.0e-8;
     tempcutoffmax[i]=1.0e+99;
     DrhoDT[i]=0.0;
     DrhoDz[i]=0.0;
     override_density[i]=0;
     temperature_error_cutoff[i]=0.0;
    }

    pp.queryarr("tempcutoff",tempcutoff,0,nmat);
    pp.queryarr("tempcutoffmax",tempcutoffmax,0,nmat);

    pp.queryarr("stiffPINF",stiffPINF,0,nmat);
    pp.queryarr("stiffCP",stiffCP,0,nmat);
    for (int i=0;i<nmat;i++)
     stiffCV[i]=stiffCP[i];

    pp.queryarr("stiffCV",stiffCV,0,nmat);
    pp.queryarr("stiffGAMMA",stiffGAMMA,0,nmat);

    pp.query("angular_velocity",angular_velocity);

    pp.query("constant_viscosity",constant_viscosity);

    pp.queryarr("DrhoDT",DrhoDT,0,nmat);
    pp.queryarr("DrhoDz",DrhoDz,0,nmat);
    pp.queryarr("override_density",override_density,0,nmat);

    pp.getarr("vorterr",vorterr,0,nmat);
    pp.query("pressure_error_flag",pressure_error_flag);
    pp.getarr("pressure_error_cutoff",pressure_error_cutoff,0,nmat);
    pp.queryarr("temperature_error_cutoff",temperature_error_cutoff,0,nmat);

    pp.query("temperature_source",temperature_source);
    pp.queryarr("temperature_source_cen",temperature_source_cen,0,AMREX_SPACEDIM);
    pp.queryarr("temperature_source_rad",temperature_source_rad,0,AMREX_SPACEDIM);

    pp.getarr("tempconst",tempconst,0,nmat);
    for (int i=0;i<nmat;i++)
     initial_temperature[i]=tempconst[i]; 
    pp.queryarr("initial_temperature",initial_temperature,0,nmat);
    pp.query("initial_temperature_diffuse_duration",
     initial_temperature_diffuse_duration);
    if (initial_temperature_diffuse_duration<0.0)
     amrex::Error("initial_temperature_diffuse_duration<0.0");

    pp.getarr("viscconst",viscconst,0,nmat);

    viscconst_min=viscconst[0];
    viscconst_max=viscconst[0];

    for (int i=0;i<nmat;i++) {

     if (viscconst[i]<viscconst_min)
      viscconst_min=viscconst[i];
     if (viscconst[i]>viscconst_max)
      viscconst_max=viscconst[i];

     viscconst_eddy[i]=0.0;
    }
    pp.queryarr("viscconst_eddy",viscconst_eddy,0,nmat);

    for (int i=0;i<nmat;i++)
     viscosity_state_model[i]=0;
    pp.queryarr("viscosity_state_model",viscosity_state_model,0,nmat);
    
    for (int i=0;i<nmat;i++)
     viscoelastic_model[i]=0;

    pp.queryarr("viscoelastic_model",viscoelastic_model,0,nmat);

    for (int i=0;i<nmat;i++) {
     if (viscoelastic_model[i]==0) {
      // do nothing
     } else if (viscoelastic_model[i]==1) {
      // do nothing
     } else if (viscoelastic_model[i]==2) {
      if (radius_cutoff[i]==-1) {
       // do nothing
      } else
       amrex::Error("expecting radius_cutoff==-1 for elastic material");
     } else
      amrex::Error("viscoelastic_model invalid");
    } // i=0..nmat-1

    for (int i=0;i<nmat;i++)
     les_model[i]=0;
    pp.queryarr("les_model",les_model,0,nmat);

    heatviscconst.resize(nmat);
    heatviscconst_interface.resize(nten);
    pp.getarr("heatviscconst",heatviscconst,0,nmat);

    heatviscconst_min=heatviscconst[0];
    heatviscconst_max=heatviscconst[0];

    for (int i=0;i<nmat;i++) {
     if (heatviscconst[i]<heatviscconst_min)
      heatviscconst_min=heatviscconst[i];
     if (heatviscconst[i]>heatviscconst_max)
      heatviscconst_max=heatviscconst[i];
    }

    for (int i=0;i<nten;i++) {
     heatviscconst_interface[i]=0.0;
    }

    pp.queryarr("heatviscconst_interface",heatviscconst_interface,0,nten);
    if (num_species_var>0) {
     pp.queryarr("speciesconst",speciesconst,0,num_species_var*nmat);
     pp.queryarr("speciesviscconst",speciesviscconst,0,num_species_var*nmat);
    }

    for (int i=0;i<nmat;i++) {
     prerecalesce_stiffCP[i]=stiffCP[i];
     prerecalesce_stiffCV[i]=stiffCV[i];
     prerecalesce_viscconst[i]=viscconst[i];
     prerecalesce_heatviscconst[i]=heatviscconst[i];
    }
    pp.queryarr("prerecalesce_viscconst",prerecalesce_viscconst,0,nmat);
    pp.queryarr("prerecalesce_heatviscconst",prerecalesce_heatviscconst,0,nmat);
    pp.queryarr("prerecalesce_stiffCP",prerecalesce_stiffCP,0,nmat);
    pp.queryarr("prerecalesce_stiffCV",prerecalesce_stiffCV,0,nmat);

    pp.query("conservative_tension_force",conservative_tension_force);
    pp.query("conservative_div_uu",conservative_div_uu);

    pp.query("mglib_min_coeff_factor",mglib_min_coeff_factor);

    pp.getarr("tension",tension,0,nten);
    for (int i=0;i<nten;i++) 
     prefreeze_tension[i]=tension[i];
    pp.queryarr("prefreeze_tension",prefreeze_tension,0,nten);
    for (int i=0;i<nten;i++) {
     tension_slope[i]=0.0;
     tension_T0[i]=293.0;
     tension_min[i]=0.0;
     cap_wave_speed[i]=0.0;
    }

    for (int i=0;i<2*AMREX_SPACEDIM;i++) {
     outflow_velocity_buffer_size[i]=0.0;
    }
    pp.queryarr("outflow_velocity_buffer_size",
      outflow_velocity_buffer_size,0,2*AMREX_SPACEDIM);

    pp.queryarr("tension_slope",tension_slope,0,nten);
    pp.queryarr("tension_T0",tension_T0,0,nten);
    pp.queryarr("tension_min",tension_min,0,nten);

    pp.queryarr("recalesce_model_parameters",recalesce_model_parameters,
       0,3*nmat);

    pp.queryarr("hardwire_Y_gamma",hardwire_Y_gamma,0,2*nten);
    pp.queryarr("hardwire_T_gamma",hardwire_T_gamma,0,2*nten);
    pp.queryarr("saturation_temp",saturation_temp,0,2*nten);
    pp.queryarr("saturation_temp_curv",saturation_temp_curv,0,2*nten);
    pp.queryarr("saturation_temp_vel",saturation_temp_vel,0,2*nten);
    pp.queryarr("saturation_temp_min",saturation_temp_min,0,2*nten);
    pp.queryarr("saturation_temp_max",saturation_temp_max,0,2*nten);

    pp.query("nucleation_period",nucleation_period);
    pp.query("nucleation_init_time",nucleation_init_time);

    pp.query("perturbation_on_restart",perturbation_on_restart);
    pp.query("perturbation_mode",perturbation_mode);
    pp.query("perturbation_eps_temp",perturbation_eps_temp);
    pp.query("perturbation_eps_vel",perturbation_eps_vel);
   
    pp.query("solidheat_flag",solidheat_flag);
    if ((solidheat_flag<0)||(solidheat_flag>2))
     amrex::Error("solidheat_flag invalid"); 
 
    pp.queryarr("microlayer_substrate",microlayer_substrate,0,nmat);
    pp.queryarr("microlayer_angle",microlayer_angle,0,nmat);
    pp.queryarr("microlayer_size",microlayer_size,0,nmat);
    pp.queryarr("macrolayer_size",macrolayer_size,0,nmat);
    pp.queryarr("max_contact_line_size",
                max_contact_line_size,0,nmat);
    pp.queryarr("microlayer_temperature_substrate",
     microlayer_temperature_substrate,0,nmat);
    for (int i=0;i<nmat;i++) {
     if (microlayer_temperature_substrate[i]<0.0)
      amrex::Error("microlayer_temperature_substrate[i]<0.0");
     if ((microlayer_substrate[i]<0)||
         (microlayer_substrate[i]>nmat))
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
    }  // i=0..nmat-1

    pp.queryarr("nucleation_temp",nucleation_temp,0,2*nten);
    pp.queryarr("nucleation_pressure",nucleation_pressure,0,2*nten);
    pp.queryarr("nucleation_pmg",nucleation_pmg,0,2*nten);
    pp.queryarr("nucleation_mach",nucleation_mach,0,2*nten);
    pp.queryarr("latent_heat",latent_heat,0,2*nten);
    pp.queryarr("reaction_rate",reaction_rate,0,2*nten);
    pp.queryarr("freezing_model",freezing_model,0,2*nten);
    pp.queryarr("Tanasawa_or_Schrage",Tanasawa_or_Schrage,0,2*nten);
    pp.queryarr("mass_fraction_id",mass_fraction_id,0,2*nten);
    pp.queryarr("distribute_from_target",distribute_from_target,0,2*nten);

    pp.query("R_Palmore_Desjardins",R_Palmore_Desjardins);

    for (int im=0;im<nmat;im++) {

     if ((override_density[im]!=0)&&
         (override_density[im]!=1)&&
         (override_density[im]!=2))
      amrex::Error("override_density invalid");
     if (DrhoDT[im]>0.0)
      amrex::Error("DrhoDT cannot be positive");
     if ((DrhoDz[im]!=-1.0)&&(DrhoDz[im]<0.0))
      amrex::Error("DrhoDz invalid"); 

     if (constant_viscosity==0) {
      // do nothing
     } else if (constant_viscosity==1) {
      if (std::abs(viscconst[im]-viscconst[0])>1.0e-12)
       amrex::Warning("variable visc but constant_viscosity==1");
     } else
      amrex::Error("constant_viscosity invalid");

     if (material_type[im]==999) {
      if (viscconst[im]<=0.0)
       amrex::Error("solid cannot have 0 viscosity");
      if (ns_is_rigid(im)!=1)
       amrex::Error("ns_is_rigid invalid");
     } else if (material_type[im]==0) {
      if (ns_is_rigid(im)!=0)
       amrex::Error("ns_is_rigid invalid");
     } else if ((material_type[im]>0)&& 
                (material_type[im]<999)) {
      if (ns_is_rigid(im)!=0)
       amrex::Error("ns_is_rigid invalid");
     } else {
      amrex::Error("material type invalid");
     }

    }  // im=0, im<nmat

    if (num_state_base!=2)
     amrex::Error("num_state_base invalid 10");

    for (int i=0;i<nten;i++) {
     if ((freezing_model[i]<0)||(freezing_model[i]>7))
      amrex::Error("freezing_model invalid in read_params (i)");
     if ((freezing_model[i+nten]<0)||(freezing_model[i+nten]>7))
      amrex::Error("freezing_model invalid in read_params (i+nten)");
     if ((distribute_from_target[i]<0)||(distribute_from_target[i]>1))
      amrex::Error("distribute_from_target invalid in read_params (i)");
     if ((distribute_from_target[i+nten]<0)||(distribute_from_target[i+nten]>1))
      amrex::Error("distribute_from_target invalid in read_params (i+nten)");
     if (mass_fraction_id[i]<0)
      amrex::Error("mass_fraction_id invalid in read_params (i)");
     if (mass_fraction_id[i+nten]<0)
      amrex::Error("mass_fraction_id invalid in read_params (i+nten)");
    }  // i=0..nten-1

    shock_timestep.resize(nmat);
    for (int i=0;i<nmat;i++) 
     shock_timestep[i]=0;
    pp.queryarr("shock_timestep",shock_timestep,0,nmat);

    int all_advective=1;

    for (int i=0;i<nmat;i++) {
     if (!((shock_timestep[i]==1)||(shock_timestep[i]==0)||
           (shock_timestep[i]==2)))
      amrex::Error("shock_timestep invalid");

     if (shock_timestep[i]==1)
      all_advective=0;
    }
    if ((cfl>1.0)&&(all_advective==1))
     amrex::Error("cfl should be less than or equal to 1");

    projection_pressure_scale=1.0;

    if (some_materials_compressible()==1) {
     projection_pressure_scale=1.0e+6;
    }

    pp.query("projection_pressure_scale",projection_pressure_scale);
    if (projection_pressure_scale<=0.0)
     amrex::Error("projection pressure scale invalid");

    projection_velocity_scale=sqrt(projection_pressure_scale);

    num_divu_outer_sweeps=1;
    if (some_materials_compressible()==1) {
     num_divu_outer_sweeps=2;
    } else if (some_materials_compressible()==0) {
     // do nothing
    } else {
     amrex::Error("some_materials_compressible invalid");
    }
    pp.query("num_divu_outer_sweeps",num_divu_outer_sweeps);
     
    pp.query("post_init_pressure_solve",post_init_pressure_solve);
    if ((post_init_pressure_solve<0)||(post_init_pressure_solve>1))
     amrex::Error("post_init_pressure_solve out of range");

    pp.query("solvability_projection",solvability_projection);

    if (solvability_projection==0) {
     // do nothing
    } else if (solvability_projection==1) {
     // do nothing
    } else
     amrex::Error("solvability_projection invalid");

    pp.query("use_lsa",use_lsa);
    if ((use_lsa<0)||(use_lsa>1))
     amrex::Error("use_lsa out of range");
    pp.query("Uref",Uref);
    pp.query("Lref",Lref);
    if (Uref<0.0)
     amrex::Error("Uref invalid");
    if (Lref<0.0)
     amrex::Error("Lref invalid");

    pp.query("pgrad_dt_factor",pgrad_dt_factor);
    if (pgrad_dt_factor<1.0)
     amrex::Error("pgrad_dt_factor too small");

    pp.query("pressure_select_criterion",pressure_select_criterion);
    if ((pressure_select_criterion<0)||
        (pressure_select_criterion>2))
     amrex::Error("pressure_select_criterion invalid");

    temperature_primitive_variable.resize(nmat);
    elastic_time.resize(nmat);

    Carreau_alpha.resize(nmat);
    Carreau_beta.resize(nmat);
    Carreau_n.resize(nmat);
    Carreau_mu_inf.resize(nmat);

    polymer_factor.resize(nmat);

    pp.query("face_flag",face_flag);

    if ((face_flag==0)||(face_flag==1)) {
     // do nothing
    } else
     amrex::Error("face_flag invalid 1");

    pp.query("disable_advection",disable_advection);

    if ((disable_advection==0)||(disable_advection==1)) {
     // do nothing
    } else
     amrex::Error("disable_advection invalid 1");

    pp.query("disable_pressure_solve",disable_pressure_solve);

    if ((disable_pressure_solve==0)||(disable_pressure_solve==1)) {
     // do nothing
    } else
     amrex::Error("disable_pressure_solve invalid 1");

    for (int i=0;i<nmat;i++) {
     temperature_primitive_variable[i]=0;
     if (is_ice_matC(i)==1) {
      temperature_primitive_variable[i]=1;
     } else if (is_FSI_rigid_matC(i)==1) {
      temperature_primitive_variable[i]=1;
     } else if (CTML_FSI_matC(i)==1) {
      temperature_primitive_variable[i]=1;
     } else if (ns_is_rigid(i)==1) {
      temperature_primitive_variable[i]=1;
     } else if (visc_coef*viscconst[i]>0.0) {
      temperature_primitive_variable[i]=1;
     } else if (material_type[i]==0) {
      temperature_primitive_variable[i]=1;
     } else if (material_type[i]==999) {
      temperature_primitive_variable[i]=1;
     } else if (ns_is_rigid(i)==0) {
      // do nothing
     } else if (visc_coef*viscconst[i]==0.0) {
      // do nothing
     } else if (material_type[i]>0) {
      // do nothing
     } else {
      amrex::Error("parameter bust");
     }

     elastic_time[i]=0.0;

     Carreau_alpha[i]=1.0;
     Carreau_beta[i]=0.0;
     Carreau_n[i]=1.0;
     Carreau_mu_inf[i]=0.0;

     polymer_factor[i]=0.0;
    }  // i=0..nmat-1

    pp.queryarr("temperature_primitive_variable",
     temperature_primitive_variable,0,nmat);

    for (int i=0;i<nmat;i++) {

     if (temperature_primitive_variable[i]==0) {
      if (ns_is_rigid(i)==1) {
       amrex::Error("make temperature_primitive_variable=1 for solids");
      } else if (is_ice_matC(i)==1) {
       amrex::Error("make temperature_primitive_variable=1 for ice");
      } else if (is_FSI_rigid_matC(i)==1) {
       amrex::Error("make temperature_primitive_variable=1 for FSI_rigid");
      } else if (CTML_FSI_matC(i)==1) {
       amrex::Error("make temperature_primitive_variable=1 for CTML");
      } else if (visc_coef*viscconst[i]>0.0) {
       amrex::Error("make temperature_primitive_variable=1 for visc. fluids");
      } else if (material_type[i]==0) {
       amrex::Error("make temperature_primitive_variable=1 for incomp fluids");
      } else if (material_type[i]==999) {
       amrex::Error("make temperature_primitive_variable=1 for solids");
      } else if (ns_is_rigid(i)==0) {
       // do nothing
      } else if (visc_coef*viscconst[i]==0.0) {
       // do nothing
      } else if (material_type[i]>0) {
       // do nothing
      } else {
       amrex::Error("parameter bust");
      }
     } else if (temperature_primitive_variable[i]==1) {
      if ((ns_is_rigid(i)==0)&& 
          (is_ice_matC(i)==0)&&
          (is_FSI_rigid_matC(i)==0)&&
          (CTML_FSI_matC(i)==0)&&
          (visc_coef*viscconst[i]==0)&& 
          (material_type[i]>=1)&&
          (material_type[i]<999)) {
       amrex::Error("make temperature_primitive_variable=1 for inv gases");
      }
     } else {
      amrex::Error("temperature_primitive_variable[i] invalid");
     }

    }  // i=0..nmat-1

    pp.queryarr("elastic_time",elastic_time,0,nmat);

    for (int i=0;i<nmat;i++) {
     if (elastic_viscosity[i]>0.0) {
      if (viscoelastic_model[i]==2) { // elastic model
       if (elastic_time[i]>=1.0e+8) {
        // do nothing
       } else
        amrex::Error("elastic time inconsistent with model");
      } else if ((viscoelastic_model[i]==1)||
   	         (viscoelastic_model[i]==0)) {
       // do nothing
      } else
       amrex::Error("viscoelastic_model invalid");
     } else if (elastic_viscosity[i]==0.0) {
      // do nothing
     } else
      amrex::Error("elastic_viscosity invalid");
    } // i=0..nmat-1

    pp.queryarr("polymer_factor",polymer_factor,0,nmat);

    pp.queryarr("Carreau_alpha",Carreau_alpha,0,nmat);
    pp.queryarr("Carreau_beta",Carreau_beta,0,nmat);
    pp.queryarr("Carreau_n",Carreau_n,0,nmat);
    pp.queryarr("Carreau_mu_inf",Carreau_mu_inf,0,nmat);

    etaL.resize(nmat);
    etaP.resize(nmat);
    etaS.resize(nmat);
    concentration.resize(nmat);

    for (int i=0;i<nmat;i++) {

     if (Carreau_n[i]>1.0)
      amrex::Error("Carreau_n[i] invalid");
     if (Carreau_mu_inf[i]<0.0)
      amrex::Error("Carreau_mu_inf[i] invalid");

     if (viscosity_state_model[i]>=0) {
      // do nothing
     } else
      amrex::Error("viscosity state model invalid");

     if ((viscoelastic_model[i]>=0)&&(viscoelastic_model[i]<=2)) {
      // do nothing
     } else
      amrex::Error("viscoelastic_model invalid");

     if (les_model[i]<0)
      amrex::Error("les model invalid");

     if ((elastic_time[i]<0.0)||(elastic_viscosity[i]<0.0))
      amrex::Error("elastic_time/elastic_viscosity invalid read_params");
     if (polymer_factor[i]<0.0)
      amrex::Error("polymer_factor invalid");

     if ((Carreau_beta[i]!=0.0)&&(visc_coef==0.0))
      amrex::Error("Cannot have Carreau_beta!=0 and visc_coef==0 ");

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

      std::cout << "etaL0=viscconst[i]=" << etaL[i] << '\n';

      std::cout << "Carreau_alpha=" << Carreau_alpha[i] << '\n';
      std::cout << "Carreau_beta=" << Carreau_beta[i] << '\n';
      std::cout << "Carreau_n=" << Carreau_n[i] << '\n';
      std::cout << "Carreau_mu_inf=" << Carreau_mu_inf[i] << '\n';
     } // io processor

     if ((num_materials_viscoelastic>=1)&&(num_materials_viscoelastic<=nmat)) {

      if (visc_coef<=0.0)
       amrex::Error("cannot have no viscosity if viscoelastic");

      if (ParallelDescriptor::IOProcessor()) {
       std::cout << "for material " << i << '\n';
       std::cout << "etaP0=elastic_viscosity=" << etaP[i] << '\n';
       std::cout << "etaS=etaL0-etaP0= " << etaS[i] << '\n';
       std::cout << "elastic_viscosity= " << elastic_viscosity[i] << '\n';
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
      amrex::Error("num_materials_viscoelastic  invalid");

    } // i=0..nmat-1

    if ((face_flag<0)||(face_flag>1))
     amrex::Error("face_flag invalid 2");

    pp.query("wait_time",wait_time);

    pp.get("advbot",advbot);
    pp.query("inflow_pressure",inflow_pressure);
    pp.query("outflow_pressure",outflow_pressure);

    pp.query( "multilevel_maxcycle",multilevel_maxcycle);

    ParmParse ppmac("mac");
    ParmParse ppcg("cg");

    pp.query( "viscous_maxiter",viscous_maxiter);
    if ((viscous_maxiter<1)||(viscous_maxiter>2)) 
     amrex::Error("viscous_maxiter should be 1 or 2");

    pp.query( "always_use_bicgstab",always_use_bicgstab);
    if ((always_use_bicgstab==0)||
        (always_use_bicgstab==1)) {
     // do nothing
    } else
     amrex::Error("always_use_bicgstab invalid");

    int cg_abec_use_bicgstab=0;
    ppcg.query("cg.abec_use_bicgstab", cg_abec_use_bicgstab);

    ppmac.query( "mac_abs_tol",mac_abs_tol);
      // mac.visc_abs_tol (not ns.visc_abs_tol)
    ppmac.query( "visc_abs_tol",visc_abs_tol);
    thermal_abs_tol=visc_abs_tol;
    ppmac.query( "thermal_abs_tol",thermal_abs_tol);
    pp.query( "minimum_relative_error",minimum_relative_error);
    pp.query( "diffusion_minimum_relative_error",
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

    pp.query("xblob2",xblob2);
    pp.query("yblob2",yblob2);
    pp.query("zblob2",zblob2);
    pp.query("radblob2",radblob2);

    pp.query("xblob3",xblob3);
    pp.query("yblob3",yblob3);
    pp.query("zblob3",zblob3);
    pp.query("radblob3",radblob3);

    pp.query("xblob4",xblob4);
    pp.query("yblob4",yblob4);
    pp.query("zblob4",zblob4);
    pp.query("radblob4",radblob4);

    pp.query("xblob5",xblob5);
    pp.query("yblob5",yblob5);
    pp.query("zblob5",zblob5);
    pp.query("radblob5",radblob5);

    pp.query("xblob6",xblob6);
    pp.query("yblob6",yblob6);
    pp.query("zblob6",zblob6);
    pp.query("radblob6",radblob6);

    pp.query("xblob7",xblob7);
    pp.query("yblob7",yblob7);
    pp.query("zblob7",zblob7);
    pp.query("radblob7",radblob7);

    pp.query("xblob8",xblob8);
    pp.query("yblob8",yblob8);
    pp.query("zblob8",zblob8);
    pp.query("radblob8",radblob8);

    pp.query("xblob9",xblob9);
    pp.query("yblob9",yblob9);
    pp.query("zblob9",zblob9);
    pp.query("radblob9",radblob9);

    pp.query("xblob10",xblob10);
    pp.query("yblob10",yblob10);
    pp.query("zblob10",zblob10);
    pp.query("radblob10",radblob10);

    xactive=0.0;
    yactive=0.0;
    zactive=0.0;
    ractive=0.0;
    ractivex=0.0;
    ractivey=0.0;
    ractivez=0.0;

    pp.query("xactive",xactive);
    pp.query("yactive",yactive);
    pp.query("zactive",zactive);
    pp.query("ractive",ractive);
    if (ractive>0.0) {
     ractivex=ractive;
     ractivey=ractive;
     ractivez=ractive;
    }
    pp.query("ractivex",ractivex);
    pp.query("ractivey",ractivey);
    pp.query("ractivez",ractivez);

     // 0 - MGPCG  1-PCG 2-MINV=I
    pp.query("project_solver_type",project_solver_type);
    pp.query("initial_cg_cycles",initial_cg_cycles);

    pp.query("initial_project_cycles",initial_project_cycles);
    if (initial_project_cycles<1)
     amrex::Error("must do at least 1 jacobi sweep at the beginning");
    initial_viscosity_cycles=initial_project_cycles;
    pp.query("initial_viscosity_cycles",initial_viscosity_cycles);
    if (initial_viscosity_cycles<1)
     amrex::Error("must do at least 1 jacobi sweep at the beginning (visc)");
    initial_thermal_cycles=initial_viscosity_cycles;
    pp.query("initial_thermal_cycles",initial_thermal_cycles);
    if (initial_thermal_cycles<1)
     amrex::Error("must do at least 1 jacobi sweep at the beginning (therm)");

    if ((project_solver_type>=0)&&(project_solver_type<=2)) {
     // do nothing
    } else
     amrex::Error("project_solver_type invalid");

    prescribe_temperature_outflow=0;
    pp.query("prescribe_temperature_outflow",prescribe_temperature_outflow);
    if ((prescribe_temperature_outflow<0)||
        (prescribe_temperature_outflow>3))
     amrex::Error("prescribe_temperature_outflow invalid");

    pp.query("diffusionface_flag",diffusionface_flag);
    if ((diffusionface_flag<0)||(diffusionface_flag>1))
     amrex::Error("diffusionface_flag invalid"); 

    pp.query("elasticface_flag",elasticface_flag);
    if ((elasticface_flag<0)||(elasticface_flag>1))
     amrex::Error("elasticface_flag invalid"); 

    pp.query("temperatureface_flag",temperatureface_flag);
    if ((temperatureface_flag!=0)&&
        (temperatureface_flag!=1))
     amrex::Error("temperatureface_flag invalid"); 

    is_phasechange=0;
    for (int i=0;i<2*nten;i++) {
     if (latent_heat[i]!=0.0) {
      is_phasechange=1;
      if ((freezing_model[i]==0)||  // Stefan model for phase change (fully sat)
          (freezing_model[i]==5)||  // Stefan model for saturated evap/cond.
          (freezing_model[i]==6)) { // Palmore and Desjardins
       if (temperatureface_flag!=0)
        amrex::Error("must have temperatureface_flag==0");
      } else if ((freezing_model[i]==1)||
                 (freezing_model[i]==2)||  // hydrate
                 (freezing_model[i]==3)||
                 (freezing_model[i]==4)||  // Tannasawa or Schrage model
		 (freezing_model[i]==7)) { // cavitation
       if (temperatureface_flag!=1)
        amrex::Error("must have temperatureface_flag==1");
      } else
       amrex::Error("freezing_model[i] invalid");

     } // latent_heat<>0
    }  // i=0;i<2*nten

    hydrate_flag=0;
    for (int i=0;i<2*nten;i++) {
     if (latent_heat[i]!=0.0)
      if (freezing_model[i]==2)
       hydrate_flag=1;
    } // i

    truncate_volume_fractions.resize(nmat);
    particle_nsubdivide.resize(nmat);
    particle_max_per_nsubdivide.resize(nmat);
    particleLS_flag.resize(nmat);

    for (int i=0;i<nmat;i++) {

     particle_nsubdivide[i]=1;
     particle_max_per_nsubdivide[i]=3;
     particleLS_flag[i]=0;

     if ((FSI_flag[i]==0)|| // tessellating
         (FSI_flag[i]==7))  // fluid, tessellating
      truncate_volume_fractions[i]=1;
     else if (is_ice_matC(i)==1) // ice, tessellating
      truncate_volume_fractions[i]=1;
     else if (FSI_flag[i]==1) // prescribed PROB.F90 solid, non-tessellating
      truncate_volume_fractions[i]=0;
     else if (FSI_flag[i]==2) // prescribed sci_clsvof.F90 solid, non-tessell.
      truncate_volume_fractions[i]=0;
     else if (FSI_flag[i]==4) // FSI CTML solid, non-tesellating
      truncate_volume_fractions[i]=0;
     else if (is_FSI_rigid_matC(i)==1) // FSI rigid solid, tessellating
      truncate_volume_fractions[i]=0;
     else
      amrex::Error("FSI_flag invalid");
    }  // i=0..nmat-1

    pp.queryarr("particle_nsubdivide",particle_nsubdivide,0,nmat);
    pp.queryarr("particle_max_per_nsubdivide",
	    particle_max_per_nsubdivide,0,nmat);

    pp.queryarr("particleLS_flag",particleLS_flag,0,nmat);

    NS_ncomp_particles=0;
    for (int i=0;i<nmat;i++) {
     if (particleLS_flag[i]==1)
      NS_ncomp_particles++;  // levelset == 0.0 for interface particles
    } // i=0..nmat-1

    pp.queryarr("truncate_volume_fractions",truncate_volume_fractions,0,nmat);
    for (int i=0;i<nmat;i++) {
     if ((truncate_volume_fractions[i]<0)||
         (truncate_volume_fractions[i]>1))
      amrex::Error("truncate_volume_fractions invalid");
     if ((particle_nsubdivide[i]<1)||
         (particle_nsubdivide[i]>6))
      amrex::Error("particle_nsubdivide invalid");
     if ((particle_max_per_nsubdivide[i]<2)||
         (particle_max_per_nsubdivide[i]>100))
      amrex::Error("particle_max_per_nsubdivide invalid");
     if ((particleLS_flag[i]<0)||
         (particleLS_flag[i]>1))
      amrex::Error("particleLS_flag invalid");
    }

    pp.query("truncate_thickness",truncate_thickness);
    if (truncate_thickness<1.0)
     amrex::Error("truncate_thickness too small");

    for (int im=1;im<=nmat;im++) {
     for (int im_opp=im+1;im_opp<=nmat;im_opp++) {
      for (int ireverse=0;ireverse<=1;ireverse++) {

       if ((im>nmat)||(im_opp>nmat))
        amrex::Error("im or im_opp bust 200cpp");
       int iten;
       get_iten_cpp(im,im_opp,iten,nmat);
       if ((iten<1)||(iten>nten))
        amrex::Error("iten invalid");
       int im_source=im;
       int im_dest=im_opp;
       if (ireverse==1) {
        im_source=im_opp;
        im_dest=im;
       }

       int indexEXP=iten+ireverse*nten-1;

       Real LL=latent_heat[indexEXP];

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
       
       if ((freezing_model[indexEXP]==4)||  // Tannasawa or Schrage
           (freezing_model[indexEXP]==5)||  // Stefan model evap/cond.
           (freezing_model[indexEXP]==6)||  // Palmore and Desjardins
	   (freezing_model[indexEXP]==7)) { // cavitation
        if (LL!=0.0) {
         int massfrac_id=mass_fraction_id[indexEXP];
         if ((massfrac_id<1)||(massfrac_id>num_species_var))
          amrex::Error("massfrac_id invalid");
         if (LL>0.0) { //evaporation
          spec_material_id_LIQUID[massfrac_id-1]=im_source;
          spec_material_id_AMBIENT[massfrac_id-1]=im_dest;
         } else if (LL<0.0) { // condensation
          spec_material_id_LIQUID[massfrac_id-1]=im_dest;
          spec_material_id_AMBIENT[massfrac_id-1]=im_source;
         } else
          amrex::Error("LL invalid");
        } else if (LL==0.0) {
         // do nothing
	} else
	 amrex::Error("LL invalid");
       } // if (freezing_model[indexEXP]==4,5 or 6)
      } // ireverse
     } //im_opp
    } // im=1..nmat

    for (int i=0;i<num_species_var;i++) {
     int im=spec_material_id_AMBIENT[i];
     if (im==0) {
      // check nothing
     } else if ((im>=1)&&(im<=nmat)) {
      if (material_type[im-1]==0) {
       if ((override_density[im-1]==1)||
           (override_density[im-1]==2)) {
        // do nothing
       } else
        amrex::Error("override_density invalid");
      } else if (material_type[im-1]==999) {
       // do nothing
      } else if (material_type[im-1]>=1) {
       // do nothing
      } else
       amrex::Error("material_type[im-1] invalid");
     } else
      amrex::Error("im invalid");
    } // i=0..num_species_var-1


    pp.query("normal_probe_size",normal_probe_size);
    if (normal_probe_size!=1)
     amrex::Error("normal_probe_size invalid");
   
    pp.query("ngrow_distance",ngrow_distance);
    if (ngrow_distance!=4)
     amrex::Error("ngrow_distance invalid");

    pp.query("ngrowFSI",ngrowFSI);
    if (ngrowFSI!=3)
     amrex::Error("ngrowFSI invalid");

    pp.query("ngrow_expansion",ngrow_expansion);
    if (ngrow_expansion!=2)
     amrex::Error("ngrow_expansion invalid");

    mof_error_ordering=0; 
    pp.query("mof_error_ordering",mof_error_ordering);
    if ((mof_error_ordering!=0)&&
        (mof_error_ordering!=1))
     amrex::Error("mof_error_ordering invalid");
    mof_ordering.resize(nmat);

    mof_ordering_override(mof_ordering,
      nmat,probtype,
      axis_dir,radblob3,
      radblob4,radblob7,
      mof_error_ordering,
      FSI_flag);

    pp.queryarr("mof_ordering",mof_ordering,0,nmat);
    for (int i=0;i<nmat;i++) {
     if ((mof_ordering[i]<0)||
         (mof_ordering[i]>nmat+1))
      amrex::Error("mof_ordering invalid");
    }


    for (int i=0;i<nmat;i++) {

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
        if (speciesviscconst[imspec*nmat+i]==0.0) {
         // do nothing
        } else if (speciesviscconst[imspec*nmat+i]>0.0) {
         // do nothing
        } else
         amrex::Error("speciesviscconst invalid");
       } // imspec

    } // i=0..nmat-1

    if (ParallelDescriptor::IOProcessor()) {

     std::cout << "temperature_source=" << temperature_source << '\n';

     std::cout << "mglib_min_coeff_factor=" << 
        mglib_min_coeff_factor << '\n';
     for (int i=0;i<AMREX_SPACEDIM;i++) {
      std::cout << "i,temperature_source_cen=" << i << ' ' <<
         temperature_source_cen[i] << '\n';
      std::cout << "i,temperature_source_rad=" << i << ' ' <<
         temperature_source_rad[i] << '\n';
     }
 
     for (int i=0;i<nten;i++) {
      std::cout << "i= " << i << " denconst_interface "  << 
        denconst_interface[i] << '\n';
      std::cout << "i= " << i << " viscconst_interface "  << 
        viscconst_interface[i] << '\n';
      std::cout << "i= " << i << " heatviscconst_interface "  << 
        heatviscconst_interface[i] << '\n';
      for (int j=0;j<num_species_var;j++) {
       std::cout << "i= " << i << " j= " << j << 
         " speciesviscconst_interface "  << 
         speciesviscconst_interface[j*nten+i] << '\n';
      }

     } // i=0 ... nten-1

     for (int j=0;j<num_species_var;j++) {
      std::cout << " j= " << j << 
         " species_evaporation_density "  <<
         species_evaporation_density[j] << '\n';
      std::cout << " j= " << j << 
         " species_molar_mass "  <<
         species_molar_mass[j] << '\n';
     }  

     std::cout << "CTML_FSI_numsolids " << CTML_FSI_numsolids << '\n';
     std::cout << "CTML_force_model " << CTML_force_model << '\n';

     std::cout << "mof_error_ordering " << 
      mof_error_ordering << '\n';

     std::cout << "ngrow_make_distance= " << 
      ngrow_make_distance << '\n';
     std::cout << "ngrow_distance= " << 
      ngrow_distance << '\n';
     std::cout << "ngrowFSI= " << 
      ngrowFSI << '\n';
     std::cout << "ngrow_expansion= " << 
      ngrow_expansion << '\n';
     std::cout << "normal_probe_size= " << 
      normal_probe_size << '\n';
     std::cout << "prescribe_temperature_outflow= " << 
      prescribe_temperature_outflow << '\n';
     std::cout << "solidheat_flag= " << solidheat_flag << '\n';
     std::cout << "diffusionface_flag= " << diffusionface_flag << '\n';
     std::cout << "elasticface_flag= " << elasticface_flag << '\n';
     std::cout << "temperatureface_flag= " << temperatureface_flag << '\n';
     std::cout << "truncate_thickness= " << truncate_thickness << '\n';
     std::cout << "face_flag= " << face_flag << '\n';
     std::cout << "disable_advection= " << disable_advection << '\n';
     std::cout << "disable_pressure_solve= " << disable_pressure_solve << '\n';
     std::cout << "nparts (im_solid_map.size())= " << 
      im_solid_map.size() << '\n';
     std::cout << "Solid_State_Type= " << Solid_State_Type << '\n';
     std::cout << "Tensor_Type= " << Tensor_Type << '\n';
     std::cout << "NUM_STATE_TYPE= " << NUM_STATE_TYPE << '\n';

     std::cout << "angular_velocity= " << angular_velocity << '\n';

     std::cout << "constant_viscosity= " << constant_viscosity << '\n';

     std::cout << "pressure_error_flag=" << pressure_error_flag << '\n';

     std::cout << "initial_temperature_diffuse_duration=" << 
      initial_temperature_diffuse_duration << '\n';

     std::cout << "make_interface_incomp " << make_interface_incomp << '\n';
 
     for (int i=0;i<nmat;i++) {
      std::cout << "mof_ordering i= " << i << ' ' <<
        mof_ordering[i] << '\n';
      std::cout << "truncate_volume_fractions i= " << i << ' ' <<
        truncate_volume_fractions[i] << '\n';

      std::cout << "particle_nsubdivide i= " << i << ' ' <<
        particle_nsubdivide[i] << '\n';
      std::cout << "particle_max_per_nsubdivide i= " << i << ' ' <<
        particle_max_per_nsubdivide[i] << '\n';
      std::cout << "particleLS_flag i= " << i << ' ' <<
        particleLS_flag[i] << '\n';

      std::cout << "NS_ncomp_particles= " << NS_ncomp_particles << '\n';

      std::cout << "viscosity_state_model i= " << i << ' ' <<
        viscosity_state_model[i] << '\n';
      std::cout << "viscoelastic_model i= " << i << ' ' <<
        viscoelastic_model[i] << '\n';
      std::cout << "les_model i= " << i << ' ' <<
        les_model[i] << '\n';
      std::cout << "temperature_primitive_variable i= " << i << ' ' <<
       temperature_primitive_variable[i] << '\n';
      std::cout << "shock_timestep i=" << i << " " << 
          shock_timestep[i] << '\n';
      std::cout << "material_type i=" << i << " " << material_type[i] << '\n';
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
      std::cout << "added_weight i=" << i << " " << added_weight[i] << '\n';
      std::cout << "denconst_gravity i=" << i << " " << 
         denconst_gravity[i] << '\n';
      std::cout << "denconst i=" << i << " " << denconst[i] << '\n';
      std::cout << "density_floor i=" << i << " " << density_floor[i] << '\n';
      std::cout << "density_ceiling i="<<i<<" "<< density_ceiling[i] << '\n';
      std::cout << "density_floor_expansion i=" << i << " " << 
        density_floor_expansion[i] << '\n';
      std::cout << "density_ceiling_expansion i="<<i<<" "<< 
        density_ceiling_expansion[i] << '\n';
      std::cout << "molar_mass i="<<i<<" "<< 
        molar_mass[i] << '\n';
      std::cout << "tempconst i=" << i << " " << tempconst[i] << '\n';
      std::cout << "initial_temperature i=" << i << " " << 
        initial_temperature[i] << '\n';
      std::cout << "tempcutoff i=" << i << " " << tempcutoff[i] << '\n';
      std::cout << "tempcutoffmax i=" << i << " " << tempcutoffmax[i] << '\n';
      std::cout << "DrhoDT i=" << i << " " << DrhoDT[i] << '\n';
      std::cout << "DrhoDz i=" << i << " " << DrhoDz[i] << '\n';
      std::cout << "override_density i=" << i << " " << 
         override_density[i] << '\n';
      std::cout << "viscconst i=" << i << "  " << viscconst[i] << '\n';
      std::cout << "viscconst_eddy i=" <<i<<"  "<<viscconst_eddy[i]<<'\n';
      std::cout << "heatviscconst i=" << i << "  " << 
          heatviscconst[i] << '\n';
      std::cout << "advection_order i=" << i << "  " << 
          advection_order[i] << '\n';
      std::cout << "density_advection_order i=" << i << "  " << 
          density_advection_order[i] << '\n';
      std::cout << "prerecalesce_viscconst i=" << i << "  " << 
         prerecalesce_viscconst[i] << '\n';
      std::cout << "prerecalesce_heatviscconst i=" << i << "  " << 
         prerecalesce_heatviscconst[i] << '\n';
      std::cout << "prerecalesce_stiffCP i=" << i << "  " << 
         prerecalesce_stiffCP[i] << '\n';
      std::cout << "prerecalesce_stiffCV i=" << i << "  " << 
         prerecalesce_stiffCV[i] << '\n';
     }  // i=0,..,nmat

     for (int i=0;i<num_species_var*nmat;i++) {
      std::cout << "speciesviscconst i=" << i << "  " << 
          speciesviscconst[i] << '\n';
      std::cout << "speciesconst i=" << i << "  " << 
          speciesconst[i] << '\n';
     }

     std::cout << "stokes_flow= " << stokes_flow << '\n';
     std::cout << "cancel_advection= " << cancel_advection << '\n';

     std::cout << "is_phasechange= " << is_phasechange << '\n';

     for (int i=0;i<3*nmat;i++) {
      std::cout << "recalesce_model_parameters i=" << i << "  " << 
       recalesce_model_parameters[i] << '\n';
     }

     std::cout << "perturbation_on_restart " << perturbation_on_restart << '\n';
     std::cout << "perturbation_mode " << perturbation_mode << '\n';
     std::cout << "perturbation_eps_temp " << perturbation_eps_temp << '\n';
     std::cout << "perturbation_eps_vel " << perturbation_eps_vel << '\n';

     std::cout << "custom_nucleation_model " << 
       custom_nucleation_model << '\n';

     std::cout << "conservative_div_uu " << 
       conservative_div_uu << '\n';
     std::cout << "conservative_tension_force " << 
       conservative_tension_force << '\n';
     std::cout << "FD_curv_interp " << FD_curv_interp << '\n';

     std::cout << "hydrate flag " << hydrate_flag << '\n';
     std::cout << "nucleation_period= " << nucleation_period << '\n';
     std::cout << "nucleation_init_time= " << nucleation_init_time << '\n';
     std::cout << "n_sites= " << n_sites << '\n';
     if (n_sites>0) {
      for (int i=0;i<pos_sites.size();i++) {
       std::cout << "i, pos_sites= " << i << ' ' << pos_sites[i] << '\n';
      }
     }
    
     for (int i=0;i<nmat;i++) {
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
     } // i=0..nmat-1

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

     for (int i=0;i<nten;i++) {
      std::cout << "hardwire_T_gamma i=" << i << "  " << 
       hardwire_T_gamma[i] << '\n';
      std::cout << "hardwire_T_gamma i+nten=" << i+nten << "  " << 
       hardwire_T_gamma[i+nten] << '\n';

      std::cout << "hardwire_Y_gamma i=" << i << "  " << 
       hardwire_Y_gamma[i] << '\n';
      std::cout << "hardwire_Y_gamma i+nten=" << i+nten << "  " << 
       hardwire_Y_gamma[i+nten] << '\n';

      std::cout << "saturation_temp i=" << i << "  " << 
       saturation_temp[i] << '\n';
      std::cout << "saturation_temp i+nten=" << i+nten << "  " << 
       saturation_temp[i+nten] << '\n';

      std::cout << "saturation_temp_curv i=" << i << "  " << 
       saturation_temp_curv[i] << '\n';
      std::cout << "saturation_temp_curv i+nten=" << i+nten << "  " << 
       saturation_temp_curv[i+nten] << '\n';

      std::cout << "saturation_temp_vel i=" << i << "  " << 
       saturation_temp_vel[i] << '\n';
      std::cout << "saturation_temp_vel i+nten=" << i+nten << "  " << 
       saturation_temp_vel[i+nten] << '\n';

      std::cout << "saturation_temp_min i=" << i << "  " << 
       saturation_temp_min[i] << '\n';
      std::cout << "saturation_temp_max i+nten=" << i+nten << "  " << 
       saturation_temp_max[i+nten] << '\n';

      std::cout << "nucleation_temp i=" << i << "  " << 
       nucleation_temp[i] << '\n';
      std::cout << "nucleation_temp i+nten=" << i+nten << "  " << 
       nucleation_temp[i+nten] << '\n';

      std::cout << "nucleation_pressure i=" << i << "  " << 
       nucleation_pressure[i] << '\n';
      std::cout << "nucleation_pressure i+nten=" << i+nten << "  " << 
       nucleation_pressure[i+nten] << '\n';

      std::cout << "nucleation_pmg i=" << i << "  " << 
       nucleation_pmg[i] << '\n';
      std::cout << "nucleation_pmg i+nten=" << i+nten << "  " << 
       nucleation_pmg[i+nten] << '\n';

      std::cout << "nucleation_mach i=" << i << "  " << 
       nucleation_mach[i] << '\n';
      std::cout << "nucleation_mach i+nten=" << i+nten << "  " << 
       nucleation_mach[i+nten] << '\n';

      std::cout << "latent_heat i=" << i << "  " << 
       latent_heat[i] << '\n';
      std::cout << "latent_heat i+nten=" << i+nten << "  " << 
       latent_heat[i+nten] << '\n';

      std::cout << "reaction_rate i=" << i << "  " << 
       reaction_rate[i] << '\n';
      std::cout << "reaction_rate i+nten=" << i+nten << "  " << 
       reaction_rate[i+nten] << '\n';

      std::cout << "freezing_model i=" << i << "  " << 
       freezing_model[i] << '\n';
      std::cout << "freezing_model i+nten=" << i+nten << "  " << 
       freezing_model[i+nten] << '\n';
      std::cout << "Tanasawa_or_Schrage i=" << i << "  " << 
       Tanasawa_or_Schrage[i] << '\n';
      std::cout << "Tanasawa_or_Schrage i+nten=" << i+nten << "  " << 
       Tanasawa_or_Schrage[i+nten] << '\n';
      std::cout << "mass_fraction_id i=" << i << "  " << 
       mass_fraction_id[i] << '\n';
      std::cout << "mass_fraction_id i+nten=" << i+nten << "  " << 
       mass_fraction_id[i+nten] << '\n';
      std::cout << "distribute_from_target i=" << i << "  " << 
       distribute_from_target[i] << '\n';
      std::cout << "distribute_from_target i+nten=" << i+nten << "  " << 
       distribute_from_target[i+nten] << '\n';


      std::cout << "tension i=" << i << "  " << tension[i] << '\n';
      std::cout << "tension_slope i=" << i << "  " << tension_slope[i] << '\n';
      std::cout << "tension_T0 i=" << i << "  " << tension_T0[i] << '\n';
      std::cout << "tension_min i=" << i << "  " << tension_min[i] << '\n';
      std::cout << "initial cap_wave_speed i=" << i << "  " << 
        cap_wave_speed[i] << '\n';
      std::cout << "prefreeze_tension i=" << i << "  " << 
       prefreeze_tension[i] << '\n';
     }  // i=0..nten-1

     for (int i=0;i<nmat;i++) {
      std::cout << "cavitation_pressure i=" << i << "  " << 
       cavitation_pressure[i] << '\n';
      std::cout << "cavitation_vapor_density i=" << i << "  " << 
       cavitation_vapor_density[i] << '\n';
      std::cout << "cavitation_tension i=" << i << "  " << 
       cavitation_tension[i] << '\n';
     } // i=0..nmat-1

     std::cout << "Uref " << Uref << '\n';
     std::cout << "Lref " << Lref << '\n';
     std::cout << "use_lsa " << use_lsa << '\n';
     std::cout << "pgrad_dt_factor " << pgrad_dt_factor << '\n';
     std::cout << "pressure_select_criterion " << 
       pressure_select_criterion << '\n';

     std::cout << "num_materials_viscoelastic " << 
        num_materials_viscoelastic << '\n';
     std::cout << "num_species_var " << num_species_var << '\n';
     std::cout << "num_materials " << num_materials << '\n';
     std::cout << "num_materials_vel " << num_materials_vel << '\n';
     std::cout << "num_materials_scalar_solve " << 
      num_materials_scalar_solve << '\n';
     std::cout << "MOFITERMAX= " << MOFITERMAX << '\n';
     std::cout << "MOF_DEBUG_RECON= " << MOF_DEBUG_RECON << '\n';
     std::cout << "MOF_TURN_OFF_LS= " << MOF_TURN_OFF_LS << '\n';

     std::cout << "post_init_pressure_solve " << 
       post_init_pressure_solve << '\n';

     std::cout << "solvability_projection " << solvability_projection << '\n';

     std::cout << "curv_stencil_height " << curv_stencil_height << '\n';

     std::cout << "projection_pressure_scale " << 
       projection_pressure_scale << '\n';
     std::cout << "projection_velocity_scale " << 
       projection_velocity_scale << '\n';

     int init_snan=FArrayBox::get_init_snan();
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

     std::cout << "bicgstab_max_num_outer_iter " << 
       bicgstab_max_num_outer_iter << '\n';
     std::cout << "slope_limiter_option " << slope_limiter_option << '\n';
     std::cout << "slipcoeff " << slipcoeff << '\n';

     std::cout << "EILE_flag " << EILE_flag << '\n';
     std::cout << "unsplit_flag " << unsplit_flag << '\n';

     std::cout << "ractive " << ractive << '\n';
     std::cout << "ractivex " << ractivex << '\n';
     std::cout << "ractivey " << ractivey << '\n';
     std::cout << "ractivez " << ractivez << '\n';
     std::cout << "wait_time " << wait_time << '\n';
     std::cout << "multilevel_maxcycle " << multilevel_maxcycle << '\n';

     std::cout << "mac.mac_abs_tol " <<mac_abs_tol<< '\n';
     std::cout << "mac.visc_abs_tol " <<visc_abs_tol<< '\n';
     std::cout << "mac.thermal_abs_tol " <<thermal_abs_tol<< '\n';
     std::cout << "viscous_maxiter " <<viscous_maxiter<< '\n';
     std::cout << "always_use_bicgstab " <<always_use_bicgstab<< '\n';
     std::cout << "cg_abec_use_bicgstab " <<cg_abec_use_bicgstab<< '\n';
     std::cout << "project_solver_type " <<project_solver_type<< '\n';
     std::cout << "initial_cg_cycles " <<initial_cg_cycles<< '\n';
     std::cout << "initial_project_cycles " <<initial_project_cycles<< '\n';
     std::cout << "initial_viscosity_cycles " <<initial_viscosity_cycles<< '\n';
     std::cout << "initial_thermal_cycles " <<initial_thermal_cycles<< '\n';
     std::cout << "visual_tessellate_vfrac " << visual_tessellate_vfrac << '\n';
     std::cout << "visual_revolve " << visual_revolve << '\n';
     std::cout << "visual_option " << visual_option << '\n';

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

    if (some_materials_compressible()==1) {
     if (num_divu_outer_sweeps<2)
      amrex::Warning("WARNING:divu_outer_sweeps>=2 for comp materials");
     if (face_flag==0) {
      if ((make_interface_incomp==0)||
          (make_interface_incomp==1)||
          (make_interface_incomp==2)) {
       // do nothing
      } else 
       amrex::Error("make_interface_incomp invalid 1");
     } else if (face_flag==1) {
      if ((make_interface_incomp==1)||
          (make_interface_incomp==2)) {
       // do nothing
      } else
       amrex::Error("make_interface_incomp invalid 2");
     } else
      amrex::Error("face_flag invalid");
    } else if (some_materials_compressible()==0) {
     // do nothing
    } else
     amrex::Error("compressible flag bust");

} // subroutine read_params()


NavierStokes::NavierStokes ()
{
    Geometry_setup();
}

// constructor
NavierStokes::NavierStokes (Amr&            papa,
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
     if (localMF_grow[i]>=0) {
      std::cout << "i= " << i << " localMF_grow= " <<
       localMF_grow[i] << '\n';
      amrex::Error("forgot to delete localMF variables");
     }

}

int NavierStokes::ns_is_rigid(int im) {

 if ((im<0)|(im>=num_materials))
  amrex::Error("im invalid50");

 int local_flag=-1;

 if ((FSI_flag[im]==0)|| // fluid, tessellating
     (FSI_flag[im]==7)|| // fluid, tessellating
     (is_FSI_rigid_matC(im)==1)|| // FSI PROB.F90 rigid solid, 5 tessellating
     (is_ice_matC(im)==1)) { //3,6 tessellating
  local_flag=0;
 } else if ((FSI_flag[im]==1)|| //prescribed PROB.F90 rigid solid,non-tess.
            (FSI_flag[im]==2)|| //prescribed sci_clsvof.F90 rigid solid,non-tess
            (FSI_flag[im]==4)) {//FSI CTML solid,non-tess.
  local_flag=1;
 } else
  amrex::Error("FSI_flag invalid");

 return local_flag;

} // subroutine ns_is_rigid


int NavierStokes::ns_is_lag_part(int im) {

 if ((im<0)|(im>=num_materials))
  amrex::Error("im invalid50");

 int local_flag=-1;

 if ((FSI_flag[im]==0)||  // fluid
     (FSI_flag[im]==3)||  // ice
     (FSI_flag[im]==5)) { //FSI rigid solid (PROB.F90)
  local_flag=0;
 } else if ((FSI_flag[im]==1)|| // prescribed PROB.F90 rigid solid
            (FSI_flag[im]==2)|| // prescribed sci_clsvof.F90 rigid solid
            (FSI_flag[im]==6)|| // lag ice
            (FSI_flag[im]==7)|| // lag fluid
            (FSI_flag[im]==4)) { // FSI CTML solid
  local_flag=1;
 } else
  amrex::Error("FSI_flag invalid");

 return local_flag;

} // subroutine ns_is_lag_part


// getState_list needs scomp,ncomp
void
NavierStokes::get_mm_scomp_solver(
  int num_materials_combine,
  int project_option,
  int& state_index,
  Vector<int>& scomp,
  Vector<int>& ncomp,
  int& ncomp_check) {

 int nmat=num_materials;

 int nsolve=1;
 int nlist=1;

 if ((project_option==0)||
     (project_option==13)||  // FSI_material_exists (1st project)
     (project_option==1)) { // pressure
  nsolve=1;
  nlist=1;
  if (num_materials_combine!=num_materials_vel)
   amrex::Error("num_materials_combine invalid");
 } else if ((project_option==10)||  // divu pressure
            (project_option==11)||  // FSI_material_exists (2nd project)
	    (project_option==12)) { // pressure extension
  nsolve=1;
  nlist=1;
  if (num_materials_combine!=num_materials_vel)
   amrex::Error("num_materials_combine invalid");
 } else if (project_option==2) { // temperature
  nsolve=1;
  nlist=num_materials_combine;
 } else if ((project_option>=100)&&
            (project_option<100+num_species_var)) { // species
  nsolve=1;
  nlist=num_materials_combine;
 } else if (project_option==3) { // viscosity
  nsolve=AMREX_SPACEDIM;
  nlist=1;
  if (num_materials_combine!=num_materials_vel)
   amrex::Error("num_materials_combine invalid");
 } else
  amrex::Error("project_option invalid1");

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 if ((num_materials_combine!=1)&&
     (num_materials_combine!=nmat)) 
  amrex::Error("num_materials_combine invalid");

 scomp.resize(nlist);
 ncomp.resize(nlist);

 int nsolveMM=nsolve*num_materials_combine;

 if ((project_option==0)||
     (project_option==1)||
     (project_option==13)|| //FSI_material_exists (1st project)
     (project_option==12)) { // pressure extrapolation
 
  scomp[0]=num_materials_vel*AMREX_SPACEDIM;
  ncomp[0]=nsolveMM; 
  state_index=State_Type;

 } else if ((project_option==10)||  // divu pressure
            (project_option==11)) { // FSI_material_exists (2nd project);
	                            // divu pressure is independent var.

  scomp[0]=0;
  ncomp[0]=nsolveMM; 
  state_index=DIV_Type;

 } else if (project_option==2) { // temperature

  // u,v,w,p,den1,T1,...
  for (int im=0;im<nlist;im++) {
   scomp[im]=num_materials_vel*(AMREX_SPACEDIM+1)+
     im*num_state_material+1;
   ncomp[im]=1;
  }
  state_index=State_Type;

 } else if (project_option==3) { // viscosity

  scomp[0]=0;
  ncomp[0]=nsolveMM; 
  state_index=State_Type;

 } else if ((project_option>=100)&&
            (project_option<100+num_species_var)) { // species

  for (int im=0;im<nlist;im++) {
   scomp[im]=num_materials_vel*(AMREX_SPACEDIM+1)+
     im*num_state_material+num_state_base+project_option-100;
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

} // get_mm_scomp_solver

void
NavierStokes::zero_independent_vel(int project_option,int idx,int nsolve) {

 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level corrupt");

 int nmat=num_materials;

 if ((nsolve!=1)&&(nsolve!=AMREX_SPACEDIM))
  amrex::Error("nsolve invalid36");

 int num_materials_face=num_materials_vel;

 if ((project_option==0)||
     (project_option==1)||
     (project_option==10)||
     (project_option==11)|| // FSI_material_exists (2nd project)
     (project_option==13)|| // FSI_material_exists (1st project)
     (project_option==12)||
     (project_option==3)) {  // viscosity
  if (num_materials_face!=1)
   amrex::Error("num_materials_face invalid");
 } else if ((project_option==2)||  // thermal diffusion
            ((project_option>=100)&&
             (project_option<100+num_species_var))) {
  num_materials_face=num_materials_scalar_solve;
 } else
  amrex::Error("project_option invalid2");

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 if ((num_materials_face!=1)&&
     (num_materials_face!=nmat))
  amrex::Error("num_materials_face invalid");

 int nsolveMM_FACE=nsolve*num_materials_face;
 if (num_materials_face==1) {
  // do nothing
 } else if (num_materials_face==nmat) { 
  nsolveMM_FACE*=2;
 } else
  amrex::Error("num_materials_face invalid");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  if (localMF[idx+dir]->nComp()!=nsolveMM_FACE)
   amrex::Error("localMF[idx+dir] has invalid ncomp");
  setVal_localMF(idx+dir,0.0,0,nsolveMM_FACE,0);
 } // dir

} // subroutine zero_independent_vel

// u,v,w,p,den1,T1,...,den2,T2,...
void
NavierStokes::zero_independent_variable(int project_option,int nsolve) {

 if (num_state_base!=2) 
  amrex::Error("num_state_base invalid");

 if ((nsolve!=1)&&
     (nsolve!=AMREX_SPACEDIM))
  amrex::Error("nsolve invalid36");

 int nmat=num_materials;

 int num_materials_face=num_materials_vel;

 if ((project_option==0)||
     (project_option==1)||
     (project_option==10)||
     (project_option==11)|| // FSI_material_exists (2nd project)
     (project_option==13)|| // FSI_material_exists (1st project)
     (project_option==12)||
     (project_option==3)) {  // viscosity
  if (num_materials_face!=1)
   amrex::Error("num_materials_face invalid");
 } else if ((project_option==2)||  // thermal diffusion
            ((project_option>=100)&&
             (project_option<100+num_species_var))) {
  num_materials_face=num_materials_scalar_solve;
 } else
  amrex::Error("project_option invalid3");

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 if ((num_materials_face!=1)&&
     (num_materials_face!=nmat))
  amrex::Error("num_materials_face invalid");

 int nsolveMM=nsolve*num_materials_face;

 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level corrupt");

 Vector<int> scomp;
 Vector<int> ncomp;
 int ncomp_check;
 int state_index;
 get_mm_scomp_solver(
   num_materials_face,
   project_option,
   state_index,
   scomp,
   ncomp,
   ncomp_check);

 if (ncomp_check!=nsolveMM)
  amrex::Error("nsolveMM invalid 2732");

 MultiFab& S_new = get_new_data(state_index,slab_step+1);
 for (int icomp=0;icomp<scomp.size();icomp++) 
  S_new.setVal(0.0,scomp[icomp],ncomp[icomp],0);

} // zero_independent_variable

void
NavierStokes::init_regrid_history() {

    is_first_step_after_regrid = 0;
    old_intersect_new          = grids;

}

void
NavierStokes::restart (Amr&          papa,
                       std::istream& is,
		       int old_finest_level,
		       int new_finest_level) {

    AmrLevel::restart(papa,is,old_finest_level,new_finest_level);
    init_regrid_history();

}

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
}


void NavierStokes::debug_ngrow(int idxMF,int ngrow,int counter) {

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
  std::cout << "counter= " << counter << '\n';
  std::cout << "idxMF= " << idxMF << '\n';
  std::cout << "mfgrow= " << mfgrow << " expected grow= " <<
   ngrow << '\n';
 }

 if (! mf->ok()) {
  std::cout << "counter= " << counter << '\n';
  std::cout << "idxMF= " << idxMF << '\n';
  amrex::Error("mf not ok");
 } else if ((mf->nGrow()<ngrow)||(mfgrow<ngrow)) {
  std::cout << "counter= " << counter << '\n';
  std::cout << "idxMF= " << idxMF << '\n';
  std::cout << "mf->ngrow= " << mf->nGrow() << " expected grow= " <<
   ngrow << '\n';
  std::cout << "mfgrow= " << mfgrow << " expected grow= " <<
   ngrow << '\n';
  amrex::Error("grow invalid in debug_ngrow");
 }

} // subroutine debug_ngrow

int NavierStokes::some_materials_compressible() {

 int comp_flag=0;
 int nmat=num_materials;
 for (int im=0;im<nmat;im++) {
  int imat_type=material_type[im];
  if (imat_type==999) {
    // do nothing
  } else if (imat_type==0) {
    // do nothing
  } else if ((imat_type>=1)&&(imat_type<999)) {
   comp_flag=1;
  } else
   amrex::Error("material type invalid");
 }
 return comp_flag;
}

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

} // function int NavierStokes::NSnumLevels()

int NavierStokes::read_from_CAD() {

 int nmat=num_materials;
 int local_read_from_CAD=0;
 for (int im=0;im<nmat;im++) {
  if ((FSI_flag[im]==2)||   // prescribed sci_clsvof.F90 rigid material 
      (FSI_flag[im]==4)||   // FSI CTML sci_clsvof.F90 material
      (FSI_flag[im]==6)||   // sci_clsvof.F90 ice
      (FSI_flag[im]==7)) {  // sci_clsvof.F90 fluid
   local_read_from_CAD=1;
  } else if ((FSI_flag[im]==0)|| // fluid
             (FSI_flag[im]==1)|| // prescribed PROB.F90 rigid material
	     (FSI_flag[im]==3)|| // ice
	     (FSI_flag[im]==5)) {// FSI rigid
   // do nothing
  } else
   amrex::Error("FSI_flag invalid");
 } // im=0..nmat-1

 return local_read_from_CAD;

}  // read_from_CAD


int NavierStokes::is_ice_matC(int im) {

 int nmat=num_materials;
 int local_is_ice=0;
 if ((im>=0)&&(im<nmat)) {
  if ((FSI_flag[im]==3)||
      (FSI_flag[im]==6)) {  
   local_is_ice=1;
  } else if ((FSI_flag[im]==0)||  // fluid
             (FSI_flag[im]==7)||  // fluid
             (FSI_flag[im]==1)||  // prescribed PROB.F90 rigid material
             (FSI_flag[im]==2)||  // prescribed sci_clsvof.F90 rigid material
             (FSI_flag[im]==4)||  // FSI CTML sci_clsvof.F90 material
	     (FSI_flag[im]==5)) { // FSI PROB.F90 rigid material
   // do nothing
  } else
   amrex::Error("FSI_flag invalid");
 } else
  amrex::Error("im invalid");

 return local_is_ice;

}  // is_ice_matC()


int NavierStokes::FSI_material_exists() {

 int local_flag=0;
 int nmat=num_materials;

 for (int im=0;im<nmat;im++) {
  if ((is_ice_matC(im)==0)&&
      (is_FSI_rigid_matC(im)==0)) {
   // do nothing
  } else if ((is_ice_matC(im)==1)||
             (is_FSI_rigid_matC(im)==1)) {
   local_flag=1;
  } else
   amrex::Error("is_ice_matC or is_FSI_rigid_matC invalid");
 } // im=0..nmat-1
 return local_flag;

}  // FSI_material_exists()


int NavierStokes::is_FSI_rigid_matC(int im) {

 int nmat=num_materials;
 int local_is_FSI_rigid=0;
 if ((im>=0)&&(im<nmat)) {
  if (FSI_flag[im]==5) {  // FSI PROB.F90 rigid material
   local_is_FSI_rigid=1;
  } else if ((FSI_flag[im]==0)||  // fluid
             (FSI_flag[im]==7)||  // fluid
             (FSI_flag[im]==1)||  // prescribed PROB.F90 rigid material
             (FSI_flag[im]==2)||  // prescribed sci_clsvof.F90 rigid material
             (FSI_flag[im]==4)||  // FSI CTML sci_clsvof.F90 material
             (FSI_flag[im]==3)||  // ice
             (FSI_flag[im]==6)) { // ice
   // do nothing
  } else
   amrex::Error("FSI_flag invalid");
 } else
  amrex::Error("im invalid");

 return local_is_FSI_rigid;

}  // is_FSI_rigid_matC()

int NavierStokes::is_singular_coeff(int im) {

 int nmat=num_materials;
 int local_is_singular_coeff=0;
 if ((im>=0)&&(im<nmat)) {
  if (FSI_flag[im]==5) {  // FSI PROB.F90 rigid material
   local_is_singular_coeff=1;
  } else if (FSI_flag[im]==1) { // prescribed PROB.F90 rigid material
   local_is_singular_coeff=1;
  } else if (FSI_flag[im]==2) { // prescribed sci_clsvof.F90 rigid material
   local_is_singular_coeff=1;
  } else if ((FSI_flag[im]==0)||
             (FSI_flag[im]==7)) { // fluid
   local_is_singular_coeff=0;
  } else if ((FSI_flag[im]==3)||
	     (FSI_flag[im]==6)) { // ice
   local_is_singular_coeff=1;
  } else if (FSI_flag[im]==4) { // FSI CTML sci_clsvof.F90 material
   local_is_singular_coeff=0;
  } else
   amrex::Error("FSI_flag invalid");
 } else
  amrex::Error("im invalid");

 return local_is_singular_coeff;

}  // is_singular_coeff()


int NavierStokes::CTML_FSI_flagC() {

 int nmat=num_materials;
 int local_CTML_FSI_flag=0;
 for (int im=0;im<nmat;im++) {
  if (FSI_flag[im]==4) {  // FSI CTML sci_clsvof.F90 
#ifdef MVAHABFSI
   local_CTML_FSI_flag=1;
#else
   amrex::Error("CTML(C): define MEHDI_VAHAB_FSI in GNUmakefile");
#endif
  } else if ((FSI_flag[im]==0)||  // fluid
             (FSI_flag[im]==7)||  // fluid
             (FSI_flag[im]==1)||  // prescribed PROB.F90 rigid material
             (FSI_flag[im]==2)||  // prescribed sci_clsvof.F90 rigid material
             (is_FSI_rigid_matC(im)==1)||  // FSI PROB.F90 rigid material
             (is_ice_matC(im)==1)) { // FSI ice material
   // do nothing
  } else
   amrex::Error("FSI_flag invalid");
 } // im=0..nmat-1

 return local_CTML_FSI_flag;

}  // CTML_FSI_flagC()


int NavierStokes::CTML_FSI_matC(int im) {

 int nmat=num_materials;
 int local_CTML_FSI_flag=0;
 if ((im>=0)&&(im<nmat)) {
  if (FSI_flag[im]==4) {  // FSI CTML sci_clsvof.F90
#ifdef MVAHABFSI
   local_CTML_FSI_flag=1;
#else
   amrex::Error("CTML(C): define MEHDI_VAHAB_FSI in GNUmakefile");
#endif
  } else if ((FSI_flag[im]==0)||  // fluid material
             (FSI_flag[im]==7)||  // fluid
             (FSI_flag[im]==1)||  // prescribed PROB.F90 rigid material
             (FSI_flag[im]==2)||  // prescribed sci_clsvof.F90 rigid material
             (is_FSI_rigid_matC(im)==1)||  // FSI PROB.F90 rigid material
             (is_ice_matC(im)==1)) { // FSI ice material
   // do nothing
  } else
   amrex::Error("FSI_flag invalid");
 } else
  amrex::Error("im invalid51");

 return local_CTML_FSI_flag;

}  // CTML_FSI_matC(int im)



// passes tile information to sci_clsvof.F90 so that Lagrangian
// elements can be distributed amongst the tiles.
// called from FSI_make_distance and ns_header_msg_level 
//  (FSI_operation==4, FSI_sub_operation==0)
// time is used just in case the actual node position depends on time.
// i.e. for finding target of characteristic given the foot.
void NavierStokes::create_fortran_grid_struct(Real time,Real dt) {

 if (read_from_CAD()!=1)
  amrex::Error("read_from_CAD()!=1");

 const int max_level = parent->maxLevel();
 int finest_level=parent->finestLevel();
 if ((level<0)||(level>max_level))
  amrex::Error("level invalid in create_fortran_grid_struct");

 Real problo[AMREX_SPACEDIM];
 Real probhi[AMREX_SPACEDIM];
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  problo[dir]=geom.ProbLo(dir);
  probhi[dir]=geom.ProbHi(dir);
 }

 bool use_tiling=ns_tiling;

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

 int nmat=num_materials;
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
 ns_reconcile_d_num(43);

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
  ns_reconcile_d_num(44);

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
 if ((nparts<1)||(nparts>nmat))
  amrex::Error("nparts invalid");

 FORT_FILLCONTAINER(
  &level,
  &finest_level,
  &max_level,
  &time,
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
  &nmat,
  &nparts,
  im_solid_map.dataPtr(),
  problo,
  probhi);

} // end subroutine create_fortran_grid_struct

// called from:
//  NavierStokes::prescribe_solid_geometryALL (if correcting solid state)
//  NavierStokes::do_the_advance (begin of divu_outer_sweeps loop)
//  NavierStokes::do_the_advance (prior to viscous diffusion)
//  NavierStokes::MaxAdvectSpeedALL
//  NavierStokes::sum_integrated_quantities
//  NavierStokes::prepare_post_process
void NavierStokes::init_FSI_GHOST_MAC_MF_ALL(int caller_id) {

 int finest_level=parent->finestLevel();
 if (level!=0)
  amrex::Error("level invalid init_FSI_GHOST_MAC_MF_ALL");

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
  //    angle = ange measured at the solid normal probe in the fluid
  //    region   grad LS_solid dot grad LS_fluid = cos(theta) ?
 if ((1==0)&&(caller_id==3)) {

  for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {

    //mat<stuff>.tec has reconstructed interface.
    // (visit can open ascii surface mesh files)
    //
    //MAC grid law of the wall information.
    //WALLFUNCTION<stuff>.plt (visit can open binary tecplot files)
   writeSanityCheckData(
    "WALLFUNCTION",
    "GNBC DEBUGGING velINT,imgVR,solVR,angle",
    caller_id,
     //velINT,image vel,velsol,image vel raster,velsol raster,angle
    localMF[HISTORY_MAC_MF+data_dir]->nComp(), 
    HISTORY_MAC_MF+data_dir,
    -1,  // State_Type==-1 
    data_dir); 

    //MAC grid rasterized solid boundary condition.
    //WALLVEL<stuff>.plt (visit can open binary tecplot files)
   writeSanityCheckData(
    "WALLVEL",
    "init_FSI_GHOST_MAC_MF_ALL, FSI_GHOST_MAC_MF",//fictitious solid velocity
    caller_id,
    localMF[FSI_GHOST_MAC_MF+data_dir]->nComp(),
    FSI_GHOST_MAC_MF+data_dir,
    -1,  // State_Type==-1 
    data_dir); 
  }
 }

 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
  delete_array(HISTORY_MAC_MF+data_dir);
 }

} // end subroutine init_FSI_GHOST_MAC_MF_ALL

//    create a ghost solid velocity variable:
//    simple method: ghost solid velocity=solid velocity
//    law of wall  : ghost solid velocity in the solid
//                   is some kind of reflection of
//                   the interior fluid velocity.  ghost 
//                   solid velocity in the fluid=fluid velocity.
// initialize Fluid Structure Interaction Ghost Multifab
// multifab = multiple fortran array blocks.
// called from:
//  NavierStokes::init_FSI_GHOST_MAC_MF_ALL
//  NavierStokes::initData ()
void NavierStokes::init_FSI_GHOST_MAC_MF(int dealloc_history) {

 int finest_level=parent->finestLevel();
 int nmat=num_materials;
 int nparts=im_solid_map.size();
 bool use_tiling=ns_tiling;

 int nparts_ghost=nparts;
 int ghost_state_type=Solid_State_Type;
 if (nparts==0) {
  nparts_ghost=1;
  ghost_state_type=State_Type;
 } else if ((nparts>=1)&&(nparts<=nmat)) {
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
 int nstate=num_materials_vel*(AMREX_SPACEDIM+1)+
  nmat*(num_state_material+ngeom_raw)+1;
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");

  // usolid_law_of_the_wall,
  // uimage raster,usolid raster,angle_ACT_cell
 int nhistory_sub=3*AMREX_SPACEDIM+1;
 int nhistory=nparts_ghost*nhistory_sub;
 int ngrow_law_of_wall=4;

 MultiFab* solid_vel_mf;
 if (nparts==0) {
  if (nparts_ghost==1) {
   solid_vel_mf=getState(ngrow_law_of_wall,0,AMREX_SPACEDIM,
    cur_time_slab);
  } else
   amrex::Error("nparts_ghost invalid");
 } else if (nparts>0) {
  solid_vel_mf=getStateSolid(ngrow_law_of_wall,0,
    nparts*AMREX_SPACEDIM,cur_time_slab);
 } else
  amrex::Error("nparts invalid");

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

  // velocity and pressure
 MultiFab* fluid_vel_mf=getState(ngrow_law_of_wall,0,AMREX_SPACEDIM+1,
    cur_time_slab);

  // temperature and density for all of the materials.
 int nden=nmat*num_state_material;
 MultiFab* state_var_mf=getStateDen(ngrow_law_of_wall,cur_time_slab);
 if (state_var_mf->nComp()==nden) {
  // do nothing
 } else
  amrex::Error("state_var_mf->nComp()!=nden");

 if (ngrow_distance==4) {
  // do nothing
 } else
  amrex::Error("ngrow_distance invalid");

  // caller_id==1
 getStateDist_localMF(LS_NRM_CP_MF,ngrow_distance,cur_time_slab,1);
 if (localMF[LS_NRM_CP_MF]->nGrow()!=ngrow_distance)
  amrex::Error("localMF[LS_NRM_CP_MF]->nGrow()!=ngrow_distance");
 if (localMF[LS_NRM_CP_MF]->nComp()!=nmat*(AMREX_SPACEDIM+1))
  amrex::Error("localMF[LS_NRM_CP_MF]->nComp()!=nmat*(AMREX_SPACEDIM+1)");

 new_localMF(LS_NRM_FD_GNBC_MF,nmat*AMREX_SPACEDIM,ngrow_distance,-1);
 build_NRM_FD_MF(LS_NRM_FD_GNBC_MF,LS_NRM_CP_MF,ngrow_distance);

 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) { 

  if (localMF_grow[HISTORY_MAC_MF+data_dir]==-1) {
   // do nothing
  } else
   amrex::Error("localMF_grow[HISTORY_MAC_MF+data_dir] invalid");

  new_localMF(HISTORY_MAC_MF+data_dir,nhistory,0,data_dir);

  if ((law_of_the_wall==0)||   //just use the solid velocity
      (law_of_the_wall==1)||   //turbulent wall flux
      (law_of_the_wall==2)) {  //GNBC

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

    FArrayBox& statefab=(*state_var_mf)[mfi];
    FArrayBox& fluidvelfab=(*fluid_vel_mf)[mfi]; 
    FArrayBox& solidvelfab=(*solid_vel_mf)[mfi]; 
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
     // in: GODUNOV_3D.F90
    FORT_WALLFUNCTION( 
     &data_dir,
     &law_of_the_wall,
     im_solid_map.dataPtr(),
     &level,
     &finest_level,
     &ngrow_law_of_wall,
     &ngrow_distance,
     &nmat,
     &nparts,
     &nparts_ghost,
     &nden,
     tilelo,tilehi,
     fablo,fabhi,&bfact,
     xlo,dx,
     &dt_slab,
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
   ns_reconcile_d_num(45); //thread_class::sync_tile_d_numPts(),
                           //ParallelDescriptor::ReduceRealSum
			   //thread_class::reconcile_d_numPts(caller_id)

   if (dealloc_history==0) {
    // do nothing
   } else if (dealloc_history==1) {
    delete_localMF(HISTORY_MAC_MF+data_dir,1);
   } else 
    amrex::Error("dealloc_history invalid");

  } else
   amrex::Error("law_of_the_wall invalid");

 } // data_dir=0..sdim-1

 delete_localMF(LS_NRM_CP_MF,1);
 delete_localMF(LS_NRM_FD_GNBC_MF,1);

 delete state_var_mf;
 delete fluid_vel_mf;
 delete solid_vel_mf;

} // end subroutine init_FSI_GHOST_MAC_MF



void NavierStokes::assimilate_state_data() {

 int finest_level=parent->finestLevel();
 int nmat=num_materials;
 int nparts=im_solid_map.size();
 bool use_tiling=ns_tiling;

 int nparts_ghost=nparts;
 int ghost_state_type=Solid_State_Type;
 if (nparts==0) {
  nparts_ghost=1;
  ghost_state_type=State_Type;
 } else if ((nparts>=1)&&(nparts<=nmat)) {
  // do nothing
 } else {
  amrex::Error("nparts invalid");
 }

 init_boundary();

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 int nstate=num_materials_vel*(AMREX_SPACEDIM+1)+
  nmat*(num_state_material+ngeom_raw)+1;
 if (nstate!=S_new.nComp()) {
  std::cout << "nstate= " << nstate << '\n';
  amrex::Error("nstate invalid in cpp assimilate");
 }

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 if (ngrow_distance==4) {
  // do nothing
 } else
  amrex::Error("ngrow_distance invalid");

 const Real* dx = geom.CellSize();

 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {

  MultiFab& Smac_new = get_new_data(Umac_Type+data_dir,slab_step+1);
	
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

   FArrayBox& ghostsolidvelfab=(*localMF[FSI_GHOST_MAC_MF+data_dir])[mfi]; 
   FArrayBox& snewfab=S_new[mfi]; 
   FArrayBox& smacnewfab=Smac_new[mfi]; 

   int tid_current=ns_thread();
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   FORT_ASSIMILATE_STATEDATA( 
     &data_dir,
     &law_of_the_wall,
     im_solid_map.dataPtr(),
     &level,
     &finest_level,
     &nstate,
     &nmat,
     &nparts,
     &nparts_ghost,
     tilelo,tilehi,
     fablo,fabhi,&bfact,
     xlo,dx,
     &dt_slab,
     &cur_time_slab,
     snewfab.dataPtr(),
     ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()),
     smacnewfab.dataPtr(),
     ARLIM(smacnewfab.loVect()),ARLIM(smacnewfab.hiVect()),
     ghostsolidvelfab.dataPtr(),
     ARLIM(ghostsolidvelfab.loVect()),ARLIM(ghostsolidvelfab.hiVect()) );
  } // mfi
} // omp
  ns_reconcile_d_num(45); //thread_class::sync_tile_d_numPts(),
                          //ParallelDescriptor::ReduceRealSum
		          //thread_class::reconcile_d_numPts(caller_id)

 } // data_dir=0..sdim-1

} // end subroutine assimilate_state_data()

// get rid of the ghost cells
void NavierStokes::resize_FSI_MF() {

 int nmat=num_materials;
 int nparts=im_solid_map.size();
 if (nparts==0) {
  // do nothing
 } else if ((nparts>=1)&&(nparts<=nmat)) {
  if (nFSI_sub!=12)
   amrex::Error("nFSI_sub invalid");
  int nFSI=nparts*nFSI_sub;
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



// create a distance function (velocity and temperature) on this level.
// calls fill coarse patch if level>0
// called from: NavierStokes::initData ()
//              NavierStokes::nonlinear_advection()
void NavierStokes::FSI_make_distance(Real time,Real dt) {

 int nmat=num_materials;
 int nparts=im_solid_map.size();

 if (nparts==0) {
  // do nothing
 } else if ((nparts>=1)&&(nparts<=nmat)) {

  // nparts x (velocity + LS + temperature + flag+stress)   3D
  if (nFSI_sub!=12)
   amrex::Error("nFSI_sub invalid");
  int nFSI=nparts*nFSI_sub;
  if (ngrowFSI!=3)
   amrex::Error("ngrowFSI invalid");

  if (localMF_grow[FSI_MF]>=0)
   delete_localMF(FSI_MF,1);

  new_localMF(FSI_MF,nFSI,ngrowFSI,-1);

  for (int partid=0;partid<nparts;partid++) {
   int ibase=partid*nFSI_sub;
   setVal_localMF(FSI_MF,0.0,ibase,3,ngrowFSI); // velocity
   setVal_localMF(FSI_MF,-99999.0,ibase+3,1,ngrowFSI); // LS
   setVal_localMF(FSI_MF,0.0,ibase+4,1,ngrowFSI); // temperature
   setVal_localMF(FSI_MF,0.0,ibase+5,1,ngrowFSI); // mask
   setVal_localMF(FSI_MF,0.0,ibase+6,6,ngrowFSI); //stress 
  } // partid=0..nparts-1

  if (read_from_CAD()==1) {

    // in: NavierStokes::FSI_make_distance
    // 1. create lagrangian container data structure within the 
    //    fortran part that recognizes tiles. (FILLCONTAINER in SOLIDFLUID.F90)
    // 2. fill the containers with the Lagrangian information.
    //    (CLSVOF_FILLCONTAINER called from FILLCONTAINER)
    //    i.e. associate to each tile a set of Lagrangian nodes and elements
    //    that are located in or very near the tile.
   create_fortran_grid_struct(time,dt);

   int iter=0; // touch_flag=0
   int FSI_operation=2;  // make distance in narrow band
   int FSI_sub_operation=0;
   resize_mask_nbr(ngrowFSI);
    // in: FSI_make_distance
    // 1.FillCoarsePatch
    // 2.traverse lagrangian elements belonging to each tile and update
    //   cells within "bounding box" of the element.
   ns_header_msg_level(FSI_operation,FSI_sub_operation,time,dt,iter);
  
   do {
 
    FSI_operation=3; // sign update   
    FSI_sub_operation=0;
    ns_header_msg_level(FSI_operation,FSI_sub_operation,time,dt,iter);
    iter++;
  
   } while (FSI_touch_flag[0]==1);

   build_moment_from_FSILS();

  } else if (read_from_CAD()==0) {
   // do nothing
  } else
   amrex::Error("read_from_CAD invalid");
 
  bool use_tiling=ns_tiling;

  const Real* dx = geom.CellSize();

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(localMF[FSI_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(*localMF[FSI_MF],use_tiling); mfi.isValid(); ++mfi) {
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

   FArrayBox& solidfab=(*localMF[FSI_MF])[mfi];

   int tid_current=ns_thread();
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // updates FSI_MF for FSI_flag(im)==1 type materials.
   FORT_INITDATASOLID(
     &nmat,
     &nparts,
     &nFSI_sub,
     &nFSI,
     &ngrowFSI,
     im_solid_map.dataPtr(),
     &time,
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,
     solidfab.dataPtr(),
     ARLIM(solidfab.loVect()),ARLIM(solidfab.hiVect()),
     dx,xlo,xhi);  
  } // mfi
} // omp
  ns_reconcile_d_num(46);

   // Solid velocity
  MultiFab& Solid_new = get_new_data(Solid_State_Type,slab_step+1);
  if (Solid_new.nComp()!=nparts*AMREX_SPACEDIM)
   amrex::Error("Solid_new.nComp()!=nparts*AMREX_SPACEDIM");
  for (int partid=0;partid<nparts;partid++) {
   int ibase=partid*nFSI_sub;
   MultiFab::Copy(Solid_new,*localMF[FSI_MF],ibase,partid*AMREX_SPACEDIM,
     AMREX_SPACEDIM,0);
  } // partid=0..nparts-1

 } else {
  amrex::Error("nparts invalid");
 }

}  // subroutine FSI_make_distance

// called from: Transfer_FSI_To_STATE()
// Transfer_FSI_To_STATE() called from: ns_header_msg_level,initData ()
void NavierStokes::copy_velocity_on_sign(int partid) {

 int nmat=num_materials;
 int nparts=im_solid_map.size();
 if ((nparts<1)||(nparts>nmat))
  amrex::Error("nparts invalid");
 if ((partid<0)||(partid>=nparts))
  amrex::Error("partid invalid");
 debug_ngrow(FSI_MF,ngrowFSI,1);

 int im_part=im_solid_map[partid];
 if ((im_part<0)||(im_part>=nmat))
  amrex::Error("im_part invalid");

 if ((FSI_flag[im_part]==2)|| // prescribed solid CAD
     (FSI_flag[im_part]==6)|| // ice CAD
     (FSI_flag[im_part]==7)) {// fluid CAD

  MultiFab& S_new=get_new_data(State_Type,slab_step+1);
  int nstate=num_materials_vel*(AMREX_SPACEDIM+1)+
   nmat*(num_state_material+ngeom_raw)+1;
  if (nstate!=S_new.nComp())
   amrex::Error("nstate invalid");

   // nparts x (velocity + LS + temperature + flag + stress) 3D
  if (nFSI_sub!=12)
   amrex::Error("nFSI_sub invalid");

  int nFSI=nparts*nFSI_sub;
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

    FORT_COPY_VEL_ON_SIGN(
     &im_part, 
     &nparts,
     &partid, 
     &ngrowFSI, 
     &nFSI, 
     &nFSI_sub, 
     xlo,dx,
     snewfab.dataPtr(),ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()),
     fsifab.dataPtr(),ARLIM(fsifab.loVect()),ARLIM(fsifab.hiVect()),
     tilelo,tilehi,
     fablo,fabhi,&bfact,
     &nmat,&nstate);
   }  // mfi  
}//omp
   ns_reconcile_d_num(47);

  } else if (ns_is_rigid(im_part)==0) {

   // do nothing
   
  } else {
   amrex::Error("ns_is_rigid invalid");
  }

 } else if (FSI_flag[im_part]==4) { // FSI CTML material
  // do nothing (CTML)
 } else
  amrex::Error("FSI_flag[im_part] invalid");

} // subroutine copy_velocity_on_sign

// called from: FSI_make_distance, initData ()
void NavierStokes::build_moment_from_FSILS() {

 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid 1");

 int nmat=num_materials;

 if (read_from_CAD()!=1)
  amrex::Error("read_from_CAD invalid");

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 MultiFab& LS_new=get_new_data(LS_Type,slab_step+1);
 int nstate=num_materials_vel*(AMREX_SPACEDIM+1)+
   nmat*(num_state_material+ngeom_raw)+1;
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");
 if (LS_new.nComp()!=nmat*(1+AMREX_SPACEDIM))
  amrex::Error("LS_new.nComp()!=nmat*(1+AMREX_SPACEDIM)");

   // nparts x (velocity + LS + temperature + flag+stress)
 if (nFSI_sub!=12)
  amrex::Error("nFSI_sub invalid");
 if (ngrowFSI!=3)
  amrex::Error("ngrowFSI!=3");
 int nparts=im_solid_map.size();
 if ((nparts<1)||(nparts>nmat))
  amrex::Error("nparts invalid");
 int nFSI=nparts*nFSI_sub;
 if (localMF[FSI_MF]->nComp()!=nFSI)
  amrex::Error("localMF[FSI_MF]->nComp()!=nFSI");
 debug_ngrow(FSI_MF,ngrowFSI,1);
  
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

  FORT_BUILD_MOMENT(
    &level,
    &finest_level,
    &nFSI, 
    &nFSI_sub, 
    &nparts,
    &ngrowFSI, 
    im_solid_map.dataPtr(),
    xlo,dx,
    snewfab.dataPtr(),ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()),
    lsnewfab.dataPtr(),ARLIM(lsnewfab.loVect()),ARLIM(lsnewfab.hiVect()),
    fsifab.dataPtr(),ARLIM(fsifab.loVect()),ARLIM(fsifab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    &nmat,&nstate);
 }  // mfi  
}//omp
 ns_reconcile_d_num(48);

} // subroutine build_moment_from_FSILS

// called from: ns_header_msg_level,initData ()
void NavierStokes::Transfer_FSI_To_STATE(Real time) {

 // nparts x (velocity + LS + temperature + flag + stress)
 int nmat=num_materials;
 int nparts=im_solid_map.size();
 if ((nparts<0)||(nparts>nmat))
  amrex::Error("nparts invalid");

 if (read_from_CAD()==1) {

  if ((nparts<1)||(nparts>nmat))
   amrex::Error("nparts invalid");
  debug_ngrow(FSI_MF,ngrowFSI,1);

  MultiFab& S_new=get_new_data(State_Type,slab_step+1);
  int nstate=num_materials_vel*(AMREX_SPACEDIM+1)+
   nmat*(num_state_material+ngeom_raw)+1;
  int dencomp=num_materials_vel*(AMREX_SPACEDIM+1);     
  if (nstate!=S_new.nComp())
   amrex::Error("nstate invalid");

  MultiFab& Solid_new = get_new_data(Solid_State_Type,slab_step+1);
  if (Solid_new.nComp()!=nparts*AMREX_SPACEDIM)
   amrex::Error("Solid_new.nComp()!=nparts*AMREX_SPACEDIM");

  MultiFab& LS_new = get_new_data(LS_Type,slab_step+1);
  if (LS_new.nComp()!=nmat*(AMREX_SPACEDIM+1))
   amrex::Error("LS_new invalid ncomp");
  if (nFSI_sub!=12)
   amrex::Error("nFSI_sub invalid");
  int nFSI=nparts*nFSI_sub;
  if (localMF[FSI_MF]->nComp()!=nFSI)
   amrex::Error("localMF[FSI_MF]->nComp()!=nFSI");

  for (int partid=0;partid<nparts;partid++) {

   int im_part=im_solid_map[partid];
   if ((im_part<0)||(im_part>=nmat))
    amrex::Error("im_part invalid");

   if ((FSI_flag[im_part]==2)|| //prescribed sci_clsvof.F90 rigid solid 
       (FSI_flag[im_part]==4)|| //FSI CTML sci_clsvof.F90 solid
       (FSI_flag[im_part]==6)|| //initial ice from CAD file
       (FSI_flag[im_part]==7)) {//initial fluid from CAD file 

    int ibase=partid*nFSI_sub;

    int ok_to_modify_EUL=1;
    if ((FSI_flag[im_part]==6)||
        (FSI_flag[im_part]==7)) {
     if (time==0.0) {
      // do nothing
     } else if (time>0.0) {
      ok_to_modify_EUL=0;
     } else
      amrex::Error("time invalid");
    } else if ((FSI_flag[im_part]==2)||
               (FSI_flag[im_part]==4)) {
     // do nothing
    } else
     amrex::Error("FSI_flag invalid");

    if (ok_to_modify_EUL==1) {

     copy_velocity_on_sign(partid);
     // Solid velocity
     //ngrow=0
     MultiFab::Copy(Solid_new,*localMF[FSI_MF],ibase,partid*AMREX_SPACEDIM,
      AMREX_SPACEDIM,0);
      // LS
      //ngrow=0
     MultiFab::Copy(LS_new,*localMF[FSI_MF],ibase+3,im_part,1,0);
      // temperature
     if (solidheat_flag==0) { // diffuse in solid
      // do nothing
     } else if ((solidheat_flag==1)||  //dirichlet
                (solidheat_flag==2)) { //neumann
       //ngrow=0
      MultiFab::Copy(S_new,*localMF[FSI_MF],ibase+4,
       dencomp+im_part*num_state_material+1,1,0);
     } else
      amrex::Error("solidheat_flag invalid"); 

    } else if (ok_to_modify_EUL==0) {
     // do nothing
    } else
     amrex::Error("ok_to_modify_EUL invalid");

   } else if (FSI_flag[im_part]==1) { // prescribed PROB.F90 rigid solid
    // do nothing
   } else
    amrex::Error("FSI_flag invalid");

  } // partid=0..nparts-1

 } else if (read_from_CAD()==0) {
  // do nothing
 } else
  amrex::Error("read_from_CAD invalid");

}  // subroutine Transfer_FSI_To_STATE

//FSI_operation=0  initialize node locations; generate_new_triangles
//FSI_operation=1  update node locations
//FSI_operation=2  make distance in narrow band
//  (nparts x (vel, LS, Temp, flag, stress)
//FSI_operation=3  update the sign.
//FSI_operation=4  copy Eulerian vel. to lag.
//
// note for CTML algorithm:
// 1. copy Eulerian velocity to Lagrangian velocity.
// 2. update node locations
// 3. copy Lagrangian Force to Eulerian Force and update Eulerian velocity.
void NavierStokes::ns_header_msg_level(
 int FSI_operation,int FSI_sub_operation,
 Real time,Real dt,int iter) {

 Vector< int > num_tiles_on_thread_proc;
 num_tiles_on_thread_proc.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  num_tiles_on_thread_proc[tid]=0;
 }

 if (FSI_operation==0) { //initialize node locations; generate_new_triangles
  if (iter!=0)
   amrex::Error("iter invalid");
  if (FSI_sub_operation!=0)
   amrex::Error("FSI_sub_operation!=0");
 } else if (FSI_operation==1) { //update node locations
  if (iter!=0)
   amrex::Error("iter invalid");
  if (FSI_sub_operation!=0)
   amrex::Error("FSI_sub_operation!=0");
  if (CTML_FSI_flagC()==1) {
   if (num_divu_outer_sweeps!=1)
    amrex::Error("num_divu_outer_sweeps!=1");
   if (ns_time_order!=1)
    amrex::Error("ns_time_order!=1");
  } else if (CTML_FSI_flagC()==0) {
   // do nothing
  } else
   amrex::Error("CTML_FSI_flagC() invalid");
 } else if (FSI_operation==2) { //make distance in narrow band
  if (iter!=0)
   amrex::Error("iter invalid");
  if (FSI_sub_operation!=0)
   amrex::Error("FSI_sub_operation!=0");
 } else if (FSI_operation==3) { //update the sign.
  if (iter<0)
   amrex::Error("iter invalid");
  if (FSI_sub_operation!=0)
   amrex::Error("FSI_sub_operation!=0");
 } else if (FSI_operation==4) { //copy Eulerian velocity to Lagrangian velocity
  if (iter!=0)
   amrex::Error("iter invalid");
  if ((FSI_sub_operation<0)||
      (FSI_sub_operation>2)) 
   amrex::Error("FSI_sub_operation invalid");
 } else
  amrex::Error("FSI_operation out of range");

 if (iter==0) {
  FSI_touch_flag.resize(thread_class::nthreads);
  for (int tid=0;tid<thread_class::nthreads;tid++) {
   FSI_touch_flag[tid]=0;
  }
 } else if (iter>0) {
  for (int tid=0;tid<thread_class::nthreads;tid++) {
   FSI_touch_flag[tid]=0;
  }
 } else {
  amrex::Error("iter invalid");
 }

 int nmat=num_materials;
 int scomp_mofvars=num_materials_vel*(AMREX_SPACEDIM+1)+
  nmat*num_state_material;
 int dencomp=num_materials_vel*(AMREX_SPACEDIM+1);     
  
 const int max_level = parent->maxLevel();
 int finest_level=parent->finestLevel();

 if ((level>max_level)||(finest_level>max_level))
  amrex::Error("(level>max_level)||(finest_level>max_level)");

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
 for (int ilev=level+1;ilev<=max_level;ilev++) 
  for (int dir=0;dir<AMREX_SPACEDIM;dir++)
   dx_max_level[dir]/=2.0;

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "ns_header_msg_level START\n";
   std::cout << "level= " << level << " finest_level= " << finest_level <<
    " max_level= " << max_level << '\n';
   std::cout << "FSI_operation= " << FSI_operation <<
    " time = " << time << " dt= " << dt << " iter = " << iter << '\n';
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

 if (FSI_operation==0) { // init node locations
  if (level==0) {
   elements_generated=0;
  } else {
   elements_generated=1;
  }
 } else if (FSI_operation==1) { // update node locations
  if (level==0) {
   elements_generated=0;
  } else {
   elements_generated=1;
  }
 } else if ((FSI_operation>=2)&&(FSI_operation<=3)) {
  elements_generated=1;
 } else if (FSI_operation==4) { // copy Eul. fluid vel to Lag. fluid vel.
  elements_generated=1;
 } else
  amrex::Error("FSI_operation invalid");

 int current_step = nStep();
 int plot_interval=parent->plotInt();

 if (read_from_CAD()==1) {

   // nparts x (velocity + LS + temperature + flag)
  int nparts=im_solid_map.size();
  if ((nparts<1)||(nparts>nmat))
   amrex::Error("nparts invalid");

  MultiFab& Solid_new=get_new_data(Solid_State_Type,slab_step+1);
  if (Solid_new.nComp()!=nparts*AMREX_SPACEDIM)
   amrex::Error("Solid_new.nComp()!=nparts*AMREX_SPACEDIM");

  MultiFab& S_new=get_new_data(State_Type,slab_step+1);
  MultiFab& LS_new = get_new_data(LS_Type,slab_step+1);
  if (LS_new.nComp()!=nmat*(AMREX_SPACEDIM+1))
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
   if (problen[dir]<=0.0)
    amrex::Error("problen[dir]<=0.0");
  }

  if (nFSI_sub!=12)
   amrex::Error("nFSI_sub invalid");
  int nFSI=nparts*nFSI_sub;

  if ((FSI_operation==0)||  // initialize nodes
      (FSI_operation==1)) { // update node locations

   if (FSI_sub_operation!=0)
    amrex::Error("FSI_sub_operation!=0");

   if (elements_generated==0) {
    int ngrowFSI_unitfab=0;
    IntVect unitlo(D_DECL(0,0,0));
    IntVect unithi(D_DECL(0,0,0));
     // construct cell-centered type box
    Box unitbox(unitlo,unithi);

    const int* tilelo=unitbox.loVect();
    const int* tilehi=unitbox.hiVect();
    const int* fablo=unitbox.loVect();
    const int* fabhi=unitbox.hiVect();

    FArrayBox FSIfab(unitbox,nFSI);

    if (num_materials_vel!=1)
     amrex::Error("num_materials_vel invalid");

    Vector<int> velbc;
    velbc.resize(num_materials_vel*AMREX_SPACEDIM*2*AMREX_SPACEDIM);
    for (int i=0;i<velbc.size();i++)
     velbc[i]=0;
    Vector<int> vofbc;
    vofbc.resize(2*AMREX_SPACEDIM);
    for (int i=0;i<vofbc.size();i++)
     vofbc[i]=0;

    int tid=0;
    int gridno=0;

    FORT_HEADERMSG(
     &tid,
     &num_tiles_on_thread_proc[tid],
     &gridno,
     &thread_class::nthreads,
     &level,
     &finest_level,
     &max_level,
     &FSI_operation, // 0 or 1 (initialize or update nodes)
     &FSI_sub_operation, // 0
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,
     problo,
     problen, 
     dx_max_level, 
     problo,
     probhi, 
     velbc.dataPtr(),  
     vofbc.dataPtr(), 
     FSIfab.dataPtr(), // placeholder
     ARLIM(FSIfab.loVect()),ARLIM(FSIfab.hiVect()),
     FSIfab.dataPtr(), // velfab spot
     ARLIM(FSIfab.loVect()),ARLIM(FSIfab.hiVect()),
     FSIfab.dataPtr(), // mnbrfab spot
     ARLIM(FSIfab.loVect()),ARLIM(FSIfab.hiVect()),
     FSIfab.dataPtr(), // mfiner spot
     ARLIM(FSIfab.loVect()),ARLIM(FSIfab.hiVect()),
     &nFSI,
     &nFSI_sub,
     &ngrowFSI_unitfab,
     &nparts,
     im_solid_map.dataPtr(),
     &h_small,
     &time, 
     &dt, 
     FSI_refine_factor.dataPtr(),
     FSI_bounding_box_ngrow.dataPtr(),
     &FSI_touch_flag[tid],
     &CTML_FSI_init,
     &CTML_force_model,
     &iter,
     &current_step,
     &plot_interval,
     &ioproc);

    elements_generated=1;
   } else if (elements_generated==1) {
    // do nothing
   } else 
    amrex::Error("elements_generated invalid");

   elements_generated=1;

   CTML_FSI_init=1;

  } else if ((FSI_operation==2)||  // make distance in narrow band
             (FSI_operation==3)) { // update the sign

   if (FSI_sub_operation!=0)
    amrex::Error("FSI_sub_operation!=0");

   elements_generated=1;

    // FSI_MF allocated in FSI_make_distance
   if (ngrowFSI!=3)
    amrex::Error("ngrowFSI invalid");
   debug_ngrow(FSI_MF,ngrowFSI,1);
   if (localMF[FSI_MF]->nComp()!=nFSI)
    amrex::Error("localMF[FSI_MF]->nComp() invalid");

   if (FSI_operation==2) { // make distance in narrow band.

    if (num_materials_vel!=1)
     amrex::Error("num_materials_vel invalid");

     // fill coarse patch 
    if (level>0) {

      //ngrow=0
     MultiFab* S_new_coarse=new MultiFab(grids,dmap,AMREX_SPACEDIM,0,
      MFInfo().SetTag("S_new_coarse"),FArrayBoxFactory());
     int dcomp=0;
     int scomp=0;
     FillCoarsePatch(*S_new_coarse,dcomp,time,State_Type,scomp,AMREX_SPACEDIM);

     if (verbose>0) {
      if (ParallelDescriptor::IOProcessor()) {
       std::cout << "check_for_NAN(S_new_coarse,200)\n";
      }
      std::fflush(NULL);
      check_for_NAN(S_new_coarse,200);
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

     FillCoarsePatch(*Solid_new_coarse,dcomp,time,Solid_State_Type,scomp,
        nparts*AMREX_SPACEDIM);

     if (verbose>0) {
      if (ParallelDescriptor::IOProcessor()) {
       std::cout << "check_for_NAN(Solid_new_coarse,200)\n";
      }
      std::fflush(NULL);
      check_for_NAN(Solid_new_coarse,201);
     }

      //ngrow=0
     MultiFab* LS_new_coarse=new MultiFab(grids,dmap,nmat*(AMREX_SPACEDIM+1),0,
      MFInfo().SetTag("LS_new_coarse"),FArrayBoxFactory());
     dcomp=0;
     scomp=0;
     FillCoarsePatch(*LS_new_coarse,dcomp,time,LS_Type,scomp,
        nmat*(AMREX_SPACEDIM+1));

     if (verbose>0) {
      if (ParallelDescriptor::IOProcessor()) {
       std::cout << "check_for_NAN(LS_new_coarse,200)\n";
      }
      std::fflush(NULL);
      check_for_NAN(LS_new_coarse,202);
     }

     for (int partid=0;partid<nparts;partid++) {

      int im_part=im_solid_map[partid];

      if ((im_part<0)||(im_part>=nmat))
       amrex::Error("im_part invalid");
 
      if ((FSI_flag[im_part]==2)|| //prescribed sci_clsvof.F90 rigid solid 
          (FSI_flag[im_part]==4)|| //FSI CTML sci_clsvof.F90 solid
	  (FSI_flag[im_part]==6)||
	  (FSI_flag[im_part]==7)) { 

       int ok_to_modify_EUL=1;
       if ((FSI_flag[im_part]==6)||
           (FSI_flag[im_part]==7)) {
	if (time==0.0) {
	 // do nothing
	} else if (time>0.0) {
	 ok_to_modify_EUL=0;
	} else
	 amrex::Error("time invalid");
       } else if ((FSI_flag[im_part]==2)||
  	          (FSI_flag[im_part]==4)) {
        // do nothing
       } else
        amrex::Error("FSI_flag invalid");

       if (ok_to_modify_EUL==1) {

        dcomp=im_part;
        scomp=im_part;
         //ngrow==0 (levelset)
        MultiFab::Copy(LS_new,*LS_new_coarse,scomp,dcomp,1,0);
        dcomp=nmat+im_part*AMREX_SPACEDIM;
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
        int scomp_thermal=dencomp+im_part*num_state_material+1;
        //ncomp==1
        FillCoarsePatch(*new_coarse_thermal,dcomp,time,State_Type,
         scomp_thermal,1);

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

      } else if (FSI_flag[im_part]==1) { // prescribed PROB.F90 rigid solid 
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

   } else if (FSI_operation==3) { // update sign
    // do not fill coarse patch.
   } else
    amrex::Error("FSI_operation invalid");

   MultiFab* solidmf=getStateSolid(ngrowFSI,0,
     nparts*AMREX_SPACEDIM,time);
   MultiFab* denmf=getStateDen(ngrowFSI,time);  
   MultiFab* LSMF=getStateDist(ngrowFSI,time,2);
   if (LSMF->nGrow()!=ngrowFSI)
    amrex::Error("LSMF->nGrow()!=ngrow_distance");

   // FSI_MF allocated in FSI_make_distance
   // all components of FSI_MF are initialized to zero except for LS.
   // LS component of FSI_MF is init to -99999
   // nparts x (velocity + LS + temperature + flag + stress)
   for (int partid=0;partid<nparts;partid++) {

    int im_part=im_solid_map[partid];
    if ((im_part<0)||(im_part>=nmat))
     amrex::Error("im_part invalid");

    int ibase=partid*nFSI_sub;
     // velocity
    MultiFab::Copy(*localMF[FSI_MF],*solidmf,partid*AMREX_SPACEDIM,
      ibase,AMREX_SPACEDIM,ngrowFSI);
     // LS  
    MultiFab::Copy(*localMF[FSI_MF],*LSMF,im_part,
      ibase+3,1,ngrowFSI);
     // temperature
    MultiFab::Copy(*localMF[FSI_MF],*denmf,im_part*num_state_material+1,
      ibase+4,1,ngrowFSI);

     // flag (mask)
    if (FSI_operation==2) {

     if ((level>0)||
         ((level==0)&&(time>0.0))) {
      setVal_localMF(FSI_MF,10.0,ibase+5,1,ngrowFSI); 
     } else if ((level==0)&&(time==0.0)) {
      // do nothing
     } else
      amrex::Error("level or time invalid");

    } else if (FSI_operation==3) {
     // do nothing
    } else
     amrex::Error("FSI_operation invalid");

   } // partid=0..nparts-1

   // (1) =1 interior  =1 fine-fine ghost in domain  =0 otherwise
   // (2) =1 interior  =0 otherwise
   // (3) =1 interior+ngrow-1  =0 otherwise
   // (4) =1 interior+ngrow    =0 otherwise
   resize_mask_nbr(ngrowFSI);
   debug_ngrow(MASK_NBR_MF,ngrowFSI,2);
 
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
    Vector<int> vofbc=getBCArray(State_Type,gridno,scomp_mofvars,1);

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    FORT_HEADERMSG(
     &tid_current,
     &num_tiles_on_thread_proc[tid_current],
     &gridno,
     &thread_class::nthreads,
     &level,
     &finest_level,
     &max_level,
     &FSI_operation, // 2 or 3 (make distance or update sign)
     &FSI_sub_operation, // 0
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,
     xlo,
     dx, 
     dx_max_level, 
     problo,
     probhi, 
     velbc.dataPtr(),  
     vofbc.dataPtr(), 
     FSIfab.dataPtr(),
     ARLIM(FSIfab.loVect()),ARLIM(FSIfab.hiVect()),
     FSIfab.dataPtr(), // velfab spot
     ARLIM(FSIfab.loVect()),ARLIM(FSIfab.hiVect()),
     mnbrfab.dataPtr(),
     ARLIM(mnbrfab.loVect()),ARLIM(mnbrfab.hiVect()),
     mnbrfab.dataPtr(), // mfiner spot
     ARLIM(mnbrfab.loVect()),ARLIM(mnbrfab.hiVect()),
     &nFSI,
     &nFSI_sub,
     &ngrowFSI,
     &nparts,
     im_solid_map.dataPtr(),
     &h_small,
     &time, 
     &dt, 
     FSI_refine_factor.dataPtr(),
     FSI_bounding_box_ngrow.dataPtr(),
     &FSI_touch_flag[tid_current],
     &CTML_FSI_init,
     &CTML_force_model,
     &iter,
     &current_step,
     &plot_interval,
     &ioproc);

    num_tiles_on_thread_proc[tid_current]++;
   } //mfi
}//omp
   ns_reconcile_d_num(49);

   for (int tid=1;tid<thread_class::nthreads;tid++) {
    if (FSI_touch_flag[tid]==1) {
     FSI_touch_flag[0]=1;
    } else if (FSI_touch_flag[tid]==0) {
     // do nothing
    } else
     amrex::Error("FSI_touch_flag[tid] invalid");
   } 
   ParallelDescriptor::ReduceIntMax(FSI_touch_flag[0]);

   if (num_materials_vel!=1)
    amrex::Error("num_materials_vel invalid");

   // idx,ngrow,scomp,ncomp,index,scompBC_map
   // InterpBordersGHOST is ultimately called.
   // dest_lstGHOST for Solid_State_Type defaults to pc_interp.
   // scompBC_map==0 corresponds to extrap_bc, pc_interp and FORT_EXTRAPFILL
   // scompBC_map==1,2,3 corresponds to x or y or z vel_extrap_bc, pc_interp 
   //   and FORT_EXTRAPFILL
   // nFSI=nparts * (vel + LS + temp + flag + stress)
   for (int partid=0;partid<nparts;partid++) {
    int ibase=partid*nFSI_sub;
    Vector<int> scompBC_map;
    scompBC_map.resize(AMREX_SPACEDIM); 
    for (int dir=0;dir<AMREX_SPACEDIM;dir++)
     scompBC_map[dir]=dir+1;

    // This routine interpolates from coarser levels.
    PCINTERP_fill_borders(FSI_MF,ngrowFSI,ibase,
     AMREX_SPACEDIM,Solid_State_Type,scompBC_map);

    for (int i=AMREX_SPACEDIM;i<nFSI_sub;i++) {
     scompBC_map.resize(1); 
     scompBC_map[0]=0;
     PCINTERP_fill_borders(FSI_MF,ngrowFSI,ibase+i,
      1,Solid_State_Type,scompBC_map);
    } // i=AMREX_SPACEDIM  ... nFSI_sub-1
   } // partid=0..nparts-1

    // 1. copy_velocity_on_sign
    // 2. update Solid_new
    // 3. update LS_new
    // 4. update S_new(temperature) (if solidheat_flag==1 or 2)

   Transfer_FSI_To_STATE(time);

   delete solidmf;
   delete denmf;
   delete LSMF;

  } else if (FSI_operation==4) { // copy Eul. vel to struct vel.

   elements_generated=1;
   if (ngrowFSI!=3)
    amrex::Error("ngrowFSI invalid");
   if (num_materials_vel!=1)
    amrex::Error("num_materials_vel invalid");
   if ((FSI_sub_operation!=0)&&
       (FSI_sub_operation!=1)&&
       (FSI_sub_operation!=2))
    amrex::Error("FSI_sub_operation invalid");

   // (1) =1 interior  =1 fine-fine ghost in domain  =0 otherwise
   // (2) =1 interior  =0 otherwise
   // (3) =1 interior+ngrow-1  =0 otherwise
   // (4) =1 interior+ngrow    =0 otherwise
   resize_mask_nbr(ngrowFSI);
   debug_ngrow(MASK_NBR_MF,ngrowFSI,2);
   // mask=1 if not covered or if outside the domain.
   // NavierStokes::maskfiner_localMF
   // NavierStokes::maskfiner
   resize_maskfiner(ngrowFSI,MASKCOEF_MF);
   debug_ngrow(MASKCOEF_MF,ngrowFSI,28);

   if ((FSI_sub_operation==0)|| //init VELADVECT_MF, fortran grid structure,...
       (FSI_sub_operation==2)) {//delete VELADVECT_MF

    if (FSI_sub_operation==0) {
     // Two layers of ghost cells are needed if
     // (INTP_CORONA = 1) in UTIL_BOUNDARY_FORCE_FSI.F90
     getState_localMF(VELADVECT_MF,ngrowFSI,0,
      num_materials_vel*AMREX_SPACEDIM,cur_time_slab); 

      // in: NavierStokes::ns_header_msg_level
     create_fortran_grid_struct(time,dt);
    } else if (FSI_sub_operation==2) {
     delete_localMF(VELADVECT_MF,1);
    } else
     amrex::Error("FSI_sub_operation invalid");

    int ngrowFSI_unitfab=0;
    IntVect unitlo(D_DECL(0,0,0));
    IntVect unithi(D_DECL(0,0,0));
     // construct cell-centered type box
    Box unitbox(unitlo,unithi);

    const int* tilelo=unitbox.loVect();
    const int* tilehi=unitbox.hiVect();
    const int* fablo=unitbox.loVect();
    const int* fabhi=unitbox.hiVect();

    FArrayBox FSIfab(unitbox,nFSI);

    if (num_materials_vel!=1)
     amrex::Error("num_materials_vel invalid");

    Vector<int> velbc;
    velbc.resize(num_materials_vel*AMREX_SPACEDIM*2*AMREX_SPACEDIM);
    for (int i=0;i<velbc.size();i++)
     velbc[i]=0;
    Vector<int> vofbc;
    vofbc.resize(2*AMREX_SPACEDIM);
    for (int i=0;i<vofbc.size();i++)
     vofbc[i]=0;

    int tid=0;
    int gridno=0;

    FORT_HEADERMSG(
     &tid,
     &num_tiles_on_thread_proc[tid],
     &gridno,
     &thread_class::nthreads,
     &level,
     &finest_level,
     &max_level,
     &FSI_operation, // 4
     &FSI_sub_operation, // 0 (clear lag data) or 2 (sync lag data)
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,
     problo,
     problen, 
     dx_max_level, 
     problo,
     probhi, 
     velbc.dataPtr(),  
     vofbc.dataPtr(), 
     FSIfab.dataPtr(), // placeholder
     ARLIM(FSIfab.loVect()),ARLIM(FSIfab.hiVect()),
     FSIfab.dataPtr(), // velfab spot
     ARLIM(FSIfab.loVect()),ARLIM(FSIfab.hiVect()),
     FSIfab.dataPtr(), // mnbrfab spot
     ARLIM(FSIfab.loVect()),ARLIM(FSIfab.hiVect()),
     FSIfab.dataPtr(), // mfiner spot
     ARLIM(FSIfab.loVect()),ARLIM(FSIfab.hiVect()),
     &nFSI,
     &nFSI_sub,
     &ngrowFSI_unitfab,
     &nparts,
     im_solid_map.dataPtr(),
     &h_small,
     &time, 
     &dt, 
     FSI_refine_factor.dataPtr(),
     FSI_bounding_box_ngrow.dataPtr(),
     &FSI_touch_flag[tid],
     &CTML_FSI_init,
     &CTML_force_model,
     &iter,
     &current_step,
     &plot_interval,
     &ioproc);

   } else if (FSI_sub_operation==1) {

    if (thread_class::nthreads<1)
     amrex::Error("thread_class::nthreads invalid");
    thread_class::init_d_numPts(localMF[VELADVECT_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
    for (MFIter mfi(*localMF[VELADVECT_MF],use_tiling); mfi.isValid(); ++mfi) {
     BL_ASSERT(grids[mfi.index()] == mfi.validbox());
     const int gridno = mfi.index();
     const Box& tilegrid = mfi.tilebox();
     const Box& fabgrid = grids[gridno];
     const int* tilelo=tilegrid.loVect();
     const int* tilehi=tilegrid.hiVect();
     const int* fablo=fabgrid.loVect();
     const int* fabhi=fabgrid.hiVect();
     const Real* xlo = grid_loc[gridno].lo();
     FArrayBox& FSIfab=(*localMF[VELADVECT_MF])[mfi]; // placeholder
     FArrayBox& velfab=(*localMF[VELADVECT_MF])[mfi]; // ngrowFSI ghost cells
     FArrayBox& mnbrfab=(*localMF[MASK_NBR_MF])[mfi];
     FArrayBox& mfinerfab=(*localMF[MASKCOEF_MF])[mfi];

     Vector<int> velbc=getBCArray(State_Type,gridno,0,
      num_materials_vel*AMREX_SPACEDIM);
     Vector<int> vofbc=getBCArray(State_Type,gridno,scomp_mofvars,1);

     int tid_current=ns_thread();
     if ((tid_current<0)||(tid_current>=thread_class::nthreads))
      amrex::Error("tid_current invalid");
     thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     FORT_HEADERMSG(
      &tid_current,
      &num_tiles_on_thread_proc[tid_current],
      &gridno,
      &thread_class::nthreads,
      &level,
      &finest_level,
      &max_level,
      &FSI_operation, // 4 (copy eul. fluid vel to lag. solid vel)
      &FSI_sub_operation, // 1 
      tilelo,tilehi,
      fablo,fabhi,
      &bfact,
      xlo,
      dx, 
      dx_max_level, 
      problo,
      probhi, 
      velbc.dataPtr(),  
      vofbc.dataPtr(), 
      FSIfab.dataPtr(), // placeholder
      ARLIM(FSIfab.loVect()),ARLIM(FSIfab.hiVect()),
      velfab.dataPtr(), // ngrowFSI ghost cells VELADVECT_MF
      ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
      mnbrfab.dataPtr(),
      ARLIM(mnbrfab.loVect()),ARLIM(mnbrfab.hiVect()),
      mfinerfab.dataPtr(),
      ARLIM(mfinerfab.loVect()),ARLIM(mfinerfab.hiVect()),
      &nFSI,
      &nFSI_sub,
      &ngrowFSI,
      &nparts,
      im_solid_map.dataPtr(),
      &h_small,
      &time, 
      &dt, 
      FSI_refine_factor.dataPtr(),
      FSI_bounding_box_ngrow.dataPtr(),
      &FSI_touch_flag[tid_current],
      &CTML_FSI_init,
      &CTML_force_model,
      &iter,
      &current_step,
      &plot_interval,
      &ioproc);

     num_tiles_on_thread_proc[tid_current]++;
    } //mfi
}//omp
    ns_reconcile_d_num(50);

   } else 
    amrex::Error("FSI_sub_operation invalid");

  } else
   amrex::Error("FSI_operation invalid");

 } else if (read_from_CAD()==0) {
  // do nothing
 } else
  amrex::Error("read_from_CAD invalid");

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "ns_header_msg_level FINISH\n";
   std::cout << "level= " << level << " finest_level= " << finest_level <<
    " max_level= " << max_level << '\n';
   std::cout << "FSI_operation= " << FSI_operation <<
    " FSI_sub_operation= " << FSI_sub_operation <<
    " time = " << time << " dt= " << dt << " iter = " << iter << '\n';
  }
 } else if (verbose==0) {
  // do nothing
 } else
  amrex::Error("verbose invalid");

} // end subroutine ns_header_msg_level

// called from Amr::restart 
void NavierStokes::post_restart() {

 SDC_setup();
 ns_time_order=parent->Time_blockingFactor();
 slab_step=ns_time_order-1;

 SDC_outer_sweeps=0;
 SDC_setup_step();

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
 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;
 int nc=S_new.nComp();

 FORT_INITDATA_ALLOC(&nmat,&nten,&nc,
  latent_heat.dataPtr(),
  freezing_model.dataPtr(),
  distribute_from_target.dataPtr(),
  saturation_temp.dataPtr(),
  dx);

 if (level==0) {

  Vector<int> bfact_space_level(max_level+1);
  Vector<int> bfact_grid_level(max_level+1);
  for (int ilev=0;ilev<=max_level;ilev++) {
   bfact_space_level[ilev]=parent->Space_blockingFactor(ilev);
   bfact_grid_level[ilev]=parent->blockingFactor(ilev);
  }
  FORT_INITGRIDMAP(
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

 int iter=0;
  // in post_restart: initialize node locations; generate_new_triangles
 int FSI_operation=0; 
 int FSI_sub_operation=0; 
 ns_header_msg_level(FSI_operation,FSI_sub_operation,
   upper_slab_time,dt_amr,iter); 

   // inside of post_restart
 if (level==0) {

  int post_init_flag=2; // post_restart
  prepare_post_process(post_init_flag);

  if (sum_interval>0) {
   sum_integrated_quantities(post_init_flag);
  }

 } else if (level>0) {
  // do nothing
 } else {
  amrex::Error("level invalid20");
 } 

}  // subroutine post_restart


// This routine might be called twice at level 0 if AMR program
// run on more than one processor.  At level=0, the BoxArray used in
// "defbaselevel" can be different from the boxarray that optimizes
// load balancing.
void
NavierStokes::initData () {

 Real strt_time=0.0;

 int bfact_space=parent->Space_blockingFactor(level);
 int bfact_grid=parent->blockingFactor(level);

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

 SDC_setup();
 ns_time_order=parent->Time_blockingFactor();
 slab_step=ns_time_order-1;

 SDC_outer_sweeps=0;
 SDC_setup_step();

 int nmat=num_materials;
 if (ngeom_raw!=AMREX_SPACEDIM+1)
  amrex::Error("ngeom_raw bust");

 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

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

 if (level==0) {

  Vector<int> bfact_space_level(max_level+1);
  Vector<int> bfact_grid_level(max_level+1);
  for (int ilev=0;ilev<=max_level;ilev++) {
   bfact_space_level[ilev]=parent->Space_blockingFactor(ilev);
   bfact_grid_level[ilev]=parent->blockingFactor(ilev);
  }
  FORT_INITGRIDMAP(
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

  int at_least_one_ice=0;
  for (int im=1;im<=nmat;im++) {
   if (is_ice_matC(im-1)==1)
    at_least_one_ice=1;
  }

  Vector<int> recalesce_material;
  recalesce_material.resize(nmat);

  int at_least_one=0;
  for (int im=1;im<=nmat;im++) {
   recalesce_material[im-1]=parent->AMR_recalesce_flag(im);
   if (parent->AMR_recalesce_flag(im)>0) {
    if (at_least_one_ice!=1)
     amrex::Error("expecting FSI_flag==3 or 6");
    at_least_one=1;
   }
  }
      
  Vector<Real> recalesce_state_old;
  int recalesce_num_state=6;
  recalesce_state_old.resize(recalesce_num_state*nmat);
  if (at_least_one==1) {
   parent->recalesce_init(nmat);
   parent->recalesce_get_state(recalesce_state_old,nmat);
  } else if (at_least_one==0) {
   for (int im=0;im<recalesce_num_state*nmat;im++) {
    recalesce_state_old[im]=-1.0;
   }
  } else
   amrex::Error("at_least_one invalid");

    // this must be done before the volume fractions, centroids, and
    // level set function are initialized.
  FORT_INITRECALESCE(
   recalesce_material.dataPtr(),
   recalesce_state_old.dataPtr(),
   &recalesce_num_state,&nmat); 

 } else if ((level>0)&&(level<=max_level)) {
  // do nothing
 } else {
  amrex::Error("level invalid 5");
 } 

 Real dt_amr=parent->getDt(); // returns dt_AMR

  // velocity,pres,state x nmat,interface variables x nmat, error ind
 MultiFab& S_new = get_new_data(State_Type,slab_step+1);
 int nc=S_new.nComp();
 int dencomp=num_materials_vel*(AMREX_SPACEDIM+1);
 int nc_expect=dencomp+nmat*num_state_material+nmat*ngeom_raw+1;
 if (nc!=nc_expect)
  amrex::Error("nc invalid in initdata");

 MultiFab& LS_new = get_new_data(LS_Type,slab_step+1);
 if (LS_new.nComp()!=nmat*(AMREX_SPACEDIM+1))
  amrex::Error("LS_new invalid ncomp");

 MultiFab& DIV_new = get_new_data(DIV_Type,slab_step+1);
 if (DIV_new.nComp()!=num_materials_vel)
  amrex::Error("DIV_new.nComp()!=num_materials_vel");

 int nparts=im_solid_map.size();

 if ((nparts>=1)&&(nparts<=nmat)) {  
  MultiFab& Solid_new = get_new_data(Solid_State_Type,slab_step+1);
  if (Solid_new.nComp()!=nparts*AMREX_SPACEDIM)
   amrex::Error("Solid_new.nComp()!=nparts*AMREX_SPACEDIM");
  Solid_new.setVal(0.0,0,nparts*AMREX_SPACEDIM,1);
 } else if (nparts==0) {
  // do nothing
 } else 
  amrex::Error("nparts invalid");

 int nparts_tensor=im_elastic_map.size();

 if ((nparts_tensor>=0)&&(nparts_tensor<=nmat)) {  
  MultiFab& Tensor_new = get_new_data(Tensor_Type,slab_step+1);
  if (Tensor_new.nComp()!=nparts_tensor*NUM_TENSOR_TYPE+AMREX_SPACEDIM)
   amrex::Error("Tensor_new.nComp()!=nparts_tensor*NUM_TENSOR_TYPE");
  Tensor_new.setVal(0.0,0,nparts_tensor*NUM_TENSOR_TYPE+AMREX_SPACEDIM,1);
 } else 
  amrex::Error("nparts_tensor invalid");

 DIV_new.setVal(0.0);

 S_new.setVal(0.0,0,nc,1);
 LS_new.setVal(-99999.0,0,nmat,1);
 LS_new.setVal(0.0,nmat,nmat*AMREX_SPACEDIM,1); // slopes

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  MultiFab& Smac_new = get_new_data(Umac_Type+dir,slab_step+1);

  int nsolve=1;
  int nsolveMM_FACE=nsolve*num_materials_vel;

  if (num_materials_vel!=1)
   amrex::Error("num_materials_vel invalid");

  if (Smac_new.nComp()!=nsolveMM_FACE) {
   std::cout << "nmat = " << nmat << '\n';
   std::cout << "num_materials_vel = " << num_materials_vel << '\n';
   std::cout << "nsolveMM_FACE = " << nsolveMM_FACE << '\n';
   amrex::Error("Smac_new.nComp() invalid in initData");
  }
  Smac_new.setVal(0.0,0,nsolveMM_FACE,0);
 }  // dir=0..sdim-1

 int iter=0; // =>  FSI_touch_flag[tid]=0
  // in initData: initialize node locations; generate_new_triangles
 int FSI_operation=0; 
 int FSI_sub_operation=0; 
 ns_header_msg_level(FSI_operation,FSI_sub_operation,
   upper_slab_time,dt_amr,iter); 

  // create a distance function (velocity and temperature) on this level.
  // calls ns_header_msg_level with FSI_operation==2,3
  // ns_header_msg_level calls NavierStokes::Transfer_FSI_To_STATE
 prepare_mask_nbr(1);
  // in: initData
 FSI_make_distance(upper_slab_time,dt_amr);

 FORT_INITDATA_ALLOC(&nmat,&nten,&nc,
  latent_heat.dataPtr(),
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
 
 // initialize one ghost cell of LS_new so that LS_stencil needed
 // by the AMR error indicator can be initialized for those
 // level set components with FSI_flag==2,4,6,7.
 MultiFab* lsmf=getStateDist(1,cur_time_slab,101);
 MultiFab::Copy(LS_new,*lsmf,0,0,nmat*(1+AMREX_SPACEDIM),1);
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

   // if FSI_flag==2,4,6,7 then LS_new is used.
   // if FSI_flag==2,4,6,7 then S_new (volume fractions and centroids)
   //  is used.
  FORT_INITDATA(
   &tid_current,
   &adapt_quad_depth,
   &level,&max_level,
   &upper_slab_time,
   tilelo,tilehi,
   fablo,fabhi,
   &bfact_space,
   &nc,
   &nmat,&nten,
   latent_heat.dataPtr(),
   saturation_temp.dataPtr(),
   radius_cutoff.dataPtr(),
   S_new[mfi].dataPtr(),
   ARLIM(S_new[mfi].loVect()),
   ARLIM(S_new[mfi].hiVect()),
   LS_new[mfi].dataPtr(),
   ARLIM(LS_new[mfi].loVect()),
   ARLIM(LS_new[mfi].hiVect()),
   dx,xlo,xhi);  

  if (1==0) {
   FArrayBox& snewfab=S_new[mfi];
   int interior_only=1;
   tecplot_debug(snewfab,xlo,fablo,fabhi,dx,-1,0,dencomp, 
    nmat*num_state_material,interior_only); 
  }
 } //mfi
}//omp
 ns_reconcile_d_num(51);

 if (read_from_CAD()==1)
  build_moment_from_FSILS();

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

  if (num_materials_vel!=1)
   amrex::Error("num_materials_vel invalid");

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  FORT_INITVELOCITY(
   &level,&upper_slab_time,
   tilelo,tilehi,
   fablo,fabhi,&bfact_space,
   S_new[mfi].dataPtr(),
   ARLIM(s_lo),ARLIM(s_hi),
   dx,xlo,xhi,
   &Reynolds,&Weber,&RGASRWATER,&use_lsa );

 } // mfi
}//omp
 ns_reconcile_d_num(52);

 // for FSI_flag==2 or 4: (prescribed sci_clsvof.F90 rigid material)
 // 1. copy_velocity_on_sign
 // 2. update Solid_new
 // 3. update LS_new
 // 4. update S_new(temperature) (if solidheat_flag==1 or 2)
 Transfer_FSI_To_STATE(strt_time);

  // if nparts>0,
  //  Initialize FSI_GHOST_MAC_MF from Solid_State_Type
  // Otherwise initialize FSI_GHOST_MAC_MF with the fluid velocity.
 int dealloc_history=1;
 init_FSI_GHOST_MAC_MF(dealloc_history);

 init_regrid_history();
 is_first_step_after_regrid=-1;

 int nstate=state.size();
 if (nstate!=NUM_STATE_TYPE)
  amrex::Error("nstate invalid");

 for (int k=0;k<nstate;k++) {
  if (debug_PC==1) {
   std::cout << "PC: before CopyNewToOld k,nstate,level,max_level " << 
     k << ' ' << nstate << ' ' <<
     level << ' ' << max_level << '\n';
  }
  state[k].CopyNewToOld(level,max_level);  // olddata=newdata 
   // time_array[0]=strt_time-dt_amr  
   // time_array[bfact_time_order]=strt_time
  state[k].setTimeLevel(strt_time,dt_amr); 
 }

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

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int nstate=state.size();
 if (nstate!=NUM_STATE_TYPE)
  amrex::Error("nstate invalid");

 int nmat=num_materials;
 int mofcomp=num_materials_vel*(AMREX_SPACEDIM+1)+
  nmat*num_state_material;
 int dencomp=num_materials_vel*(AMREX_SPACEDIM+1);
 int nden=nmat*num_state_material;

 for (int k=0;k<nstate;k++) {

  if (k==State_Type) {
   MultiFab& S_new=get_new_data(State_Type,slab_step+1);
   MultiFab* vofmf=getState(1,mofcomp,nmat*ngeom_raw,cur_time_slab);
   MultiFab::Copy(S_new,*vofmf,0,mofcomp,nmat*ngeom_raw,1);
   delete vofmf;
   MultiFab* velmf=getState(1,0,
    num_materials_vel*(AMREX_SPACEDIM+1),cur_time_slab);
   MultiFab::Copy(S_new,*velmf,0,0,
    num_materials_vel*(AMREX_SPACEDIM+1),1);
   delete velmf;
   MultiFab* denmf=getStateDen(1,cur_time_slab);  
   MultiFab::Copy(S_new,*denmf,0,dencomp,nden,1);
   delete denmf;
  } else if (k==Umac_Type) {
   // do nothing
  } else if (k==Vmac_Type) {
   // do nothing
  } else if (k==Wmac_Type) {
   // do nothing
  } else if (k==LS_Type) {
   MultiFab& LS_new=get_new_data(LS_Type,slab_step+1);
   MultiFab* lsmf=getStateDist(1,cur_time_slab,3);  
   MultiFab::Copy(LS_new,*lsmf,0,0,nmat*(1+AMREX_SPACEDIM),1);
   delete lsmf;
  } else if (k==DIV_Type) {
   // do nothing
  } else if (k==Solid_State_Type) {
   int nparts=im_solid_map.size();
   if ((nparts<1)||(nparts>nmat))
    amrex::Error("nparts invalid");
   MultiFab& Solid_new=get_new_data(Solid_State_Type,slab_step+1);
   MultiFab* velmf=getStateSolid(1,0,nparts*AMREX_SPACEDIM,cur_time_slab);
   MultiFab::Copy(Solid_new,*velmf,0,0,nparts*AMREX_SPACEDIM,1);
   delete velmf;
  } else if (k==Tensor_Type) {
   int nparts=im_elastic_map.size();
   if ((nparts<0)||(nparts>nmat))
    amrex::Error("nparts invalid");
   if (nparts!=num_materials_viscoelastic)
    amrex::Error("nparts!=num_materials_viscoelastic");
   MultiFab& Tensor_new=get_new_data(Tensor_Type,slab_step+1);
     // ngrow=1 scomp=0
   MultiFab* tensormf=getStateTensor(1,0,nparts*NUM_TENSOR_TYPE+
		   AMREX_SPACEDIM,cur_time_slab);
   MultiFab::Copy(Tensor_new,*tensormf,0,0,nparts*NUM_TENSOR_TYPE+
		   AMREX_SPACEDIM,1);
   delete tensormf;
  } else 
   amrex::Error("k invalid");

 } // k=0..nstate-1

}  // subroutine init_boundary()

void
NavierStokes::init(
  AmrLevel & old,
  const BoxArray& ba_in,  // BoxArray of "this" (new amr_level)
  const DistributionMapping& dmap_in) { // dmap of "this" (new amr_level)
 
 NavierStokes* oldns     = (NavierStokes*) &old;

 int max_level = parent->maxLevel();

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

 int nmat=num_materials;
 if (ngeom_raw!=AMREX_SPACEDIM+1)
  amrex::Error("ngeom_raw bust");

 int nstate=state.size();  // cell centered, vel MAC, LS, etc
 if (nstate!=NUM_STATE_TYPE)
  amrex::Error("nstate invalid");

 for (int k=0;k<nstate;k++) {

  MultiFab &S_new = get_new_data(k,ns_time_order);

  if (S_new.DistributionMap()==dmap_in) {
   // do nothing
  } else {
   amrex::Error("dmap_in invalid");
  }
   //  Are the BoxArrays equal after conversion to cell-centered?
  if (S_new.boxArray().CellEqual(ba_in)) {
   // do nothing
  } else {
   amrex::Error("S_new.boxArray().CellEqual(ba_in) failed");
  }


  int ncomp=S_new.nComp();

  int numparts=1;
  int ncomp_part[4];
  int scomp_part[4];
  
  if (k==State_Type) {
   int test_ncomp=0;

   numparts=4;
   scomp_part[0]=0;
   scomp_part[1]=num_materials_vel*(AMREX_SPACEDIM+1);
   scomp_part[2]=scomp_part[1]+nmat*num_state_material;
   scomp_part[3]=scomp_part[2]+nmat*ngeom_raw;

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

  } else {
   numparts=1;
   scomp_part[0]=0;
   ncomp_part[0]=ncomp;
  }
   
  for (int part_iter=0;part_iter<numparts;part_iter++) { 
   FillPatch(old,S_new,scomp_part[part_iter],
      upper_slab_time,k,
      scomp_part[part_iter],
      ncomp_part[part_iter]);
  }
 }  // k=0..nstate-1

 if (level==0) {
   // old particle data will be deleted, so that the data
   // must be saved here.
  for (int ipart=0;ipart<NS_ncomp_particles;ipart++) {
   int lev_min=0;
   int lev_max=-1;
   int nGrow_Redistribute=0;
   int local_Redistribute=0;

   AmrParticleContainer<N_EXTRA_REAL,0,0,0>& old_PC=
    oldns->get_new_dataPC(State_Type,ns_time_order,ipart);

    // new_PC inherits the old box structure.
   AmrParticleContainer<N_EXTRA_REAL,0,0,0>& new_PC=
    get_new_dataPC(State_Type,ns_time_order,ipart);

    // Redistribute called since level=0 grids might have changed.
   old_PC.Redistribute(lev_min,lev_max,nGrow_Redistribute, 
     local_Redistribute);
    //new_PC cleared prior to copy
   bool local_copy_flag=false; 
   new_PC.copyParticles(old_PC,local_copy_flag);

  } // ipart=0..NS_ncomp_particles-1
 } else if ((level>=1)&&(level<=max_level)) {
  // do nothing
 } else
  amrex::Error("level invalid");

 old_intersect_new = amrex::intersect(grids,oldns->boxArray());
 is_first_step_after_regrid = 1;

}  // subroutine init(old)


// init a new level that did not exist on the previous step.
void
NavierStokes::init (const BoxArray& ba_in,
         const DistributionMapping& dmap_in) {

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int nmat=num_materials;
 if (ngeom_raw!=AMREX_SPACEDIM+1)
  amrex::Error("ngeom_raw bust");

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

 int nstate=state.size();
 if (nstate!=NUM_STATE_TYPE)
  amrex::Error("nstate invalid");

 for (int k=0;k<nstate;k++) {

  MultiFab &S_new = get_new_data(k,ns_time_order);
  int ncomp=S_new.nComp();

  int numparts=1;
  int ncomp_part[4];
  int scomp_part[4];
  
  if (k==State_Type) {
   int test_ncomp=0;

   numparts=4;
   scomp_part[0]=0;
   scomp_part[1]=num_materials_vel*(AMREX_SPACEDIM+1);
   scomp_part[2]=scomp_part[1]+nmat*num_state_material;
   scomp_part[3]=scomp_part[2]+nmat*ngeom_raw;

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

  } else {
   numparts=1;
   scomp_part[0]=0;
   ncomp_part[0]=ncomp;
  }

  for (int part_iter=0;part_iter<numparts;part_iter++) { 
   FillCoarsePatch(S_new,scomp_part[part_iter],upper_slab_time,k,
     scomp_part[part_iter],ncomp_part[part_iter]);
  }
 } // k=0..nstate-1

 init_regrid_history();
 is_first_step_after_regrid = 2;

}  // subroutine init()

void NavierStokes::CopyNewToOldALL() {

 int finest_level=parent->finestLevel();
 const int max_level = parent->maxLevel();
 if (level!=0)
  amrex::Error("level invalid CopyNewToOldALL");
 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
   // oldtime=newtime newtime+=dt
  int nstate=ns_level.state.size();
  if (nstate!=NUM_STATE_TYPE)
   amrex::Error("nstate invalid");
  for (int k = 0; k < nstate; k++) {
   ns_level.state[k].CopyNewToOld(level,max_level);  
  }
 }
 int nmat=num_materials;

 int at_least_one=0;
 for (int im=1;im<=nmat;im++) {
  if (parent->AMR_recalesce_flag(im)>0)
   at_least_one=1;
 }
 if (at_least_one==1) {
  parent->recalesce_copy_new_to_old(nmat);
 } else if (at_least_one==0) {
  // do nothing
 } else
  amrex::Error("at_least_one invalid");

}  // subroutine CopyNewToOldALL


void NavierStokes::CopyOldToNewALL() {

 int finest_level=parent->finestLevel();
 const int max_level = parent->maxLevel();
 if (level!=0)
  amrex::Error("level invalid CopyOldToNewALL");
 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
   // oldtime=newtime=start_time newtime+=dt
  int nstate=ns_level.state.size();
  if (nstate!=NUM_STATE_TYPE)
   amrex::Error("nstate invalid");
  for (int k = 0; k < nstate; k++) {
   ns_level.state[k].CopyOldToNew(level,max_level);
  }
 }

 int nmat=num_materials;

 int at_least_one=0;
 for (int im=1;im<=nmat;im++) {
  if (parent->AMR_recalesce_flag(im)>0)
   at_least_one=1;
 }
 if (at_least_one==1) {
  parent->recalesce_copy_old_to_new(nmat);
 } else if (at_least_one==0) {
  // do nothing
 } else
  amrex::Error("at_least_one invalid");

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

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 int nsolve=1;

 if (localMF_grow[MDOT_MF]>=0) {
  delete_localMF(MDOT_MF,1);
 } 

  // MDOT has nsolve components.
 new_localMF(MDOT_MF,nsolve,0,-1);
 setVal_localMF(MDOT_MF,0.0,0,nsolve,0); //val,scomp,ncomp,ngrow

} // end subroutine allocate_mdot()

// slab_step and SDC_outer_sweeps are set before calling this routine.
void
NavierStokes::SDC_setup_step() {

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int nmat=num_materials;
 if ((nmat<1)||(nmat>1000))
  amrex::Error("nmat out of range");

 nfluxSEM=AMREX_SPACEDIM+num_state_base;
  // I-scheme,thermal conduction,viscosity,div(up),gp,-force
 nstate_SDC=nfluxSEM+1+AMREX_SPACEDIM+1+AMREX_SPACEDIM+AMREX_SPACEDIM;

 ns_time_order=parent->Time_blockingFactor();

 if ((SDC_outer_sweeps<0)||
     (SDC_outer_sweeps>=ns_time_order))
  amrex::Error("SDC_outer_sweeps invalid");

 divu_outer_sweeps=0;

 lower_slab_time=state[State_Type].slabTime(0);
 upper_slab_time=state[State_Type].slabTime(ns_time_order);
 delta_slab_time=upper_slab_time-lower_slab_time;
 if (delta_slab_time<0.0) {
  std::cout << "lower_slab_time= " << lower_slab_time << '\n';
  std::cout << "upper_slab_time= " << upper_slab_time << '\n';
  std::cout << "delta_slab_time= " << delta_slab_time << '\n';
  amrex::Error("delta_slab_time invalid");
 }
 if ((delta_slab_time==0.0)&&(upper_slab_time>0.0)) {
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
  if (dt_slab<0.0)
   amrex::Error("dt_slab invalid1");
  if ((dt_slab==0.0)&&(cur_time_slab>0.0))
   amrex::Error("dt_slab invalid2");

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
   localMF_grow[i]=-1;
   localMF[i]=0;
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

void NavierStokes::SOD_SANITY_CHECK(int id) {

 MultiFab& snew=get_new_data(State_Type,slab_step+1);
 int nc=snew.nComp();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(snew.boxArray().d_numPts());

 for (MFIter mfi(snew,false); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  const Box& tilegrid = mfi.tilebox();
  const int gridno = mfi.index();
  const Box& fabgrid = grids[gridno];
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();
  FArrayBox& snewfab=snew[mfi];

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  FORT_SOD_SANITY(
    &id,&nc,fablo,fabhi,
    snewfab.dataPtr(),ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()));
 }
 ns_reconcile_d_num(53);

} // subroutine SOD_SANITY_CHECK

void NavierStokes::make_viscoelastic_tensor(int im) {

 int nmat=num_materials;
 bool use_tiling=ns_tiling;

 if ((num_materials_viscoelastic>=1)&&(num_materials_viscoelastic<=nmat)) {
  // do nothing
 } else
  amrex::Error("num_materials_viscoelastic invalid");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 if (localMF_grow[VISCOTEN_MF]==-1) {
  // do nothing
 } else 
  amrex::Error("VISCOTEN_MF should not be allocated");

 int ngrow=1;  // number of grow cells for the tensor
 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(FACE_VAR_MF+dir,0,2);
 debug_ngrow(CELL_VISC_MATERIAL_MF,ngrow,3);

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

 // 1. viscosity coefficient - 1..nmat
 // 2. viscoelastic coefficient - 1..nmat
 // 3. relaxation time - 1..nmat
 // the viscous and viscoelastic forces should both be multiplied by
 // visc_coef.  
 if (localMF[CELL_VISC_MATERIAL_MF]->nComp()<nmat)
  amrex::Error("cell_visc_material ncomp invalid");
 if (localMF[CELL_VISC_MATERIAL_MF]->nGrow()<ngrow)
  amrex::Error("cell_visc_material ngrow invalid");

 int nstate=num_materials_vel*(AMREX_SPACEDIM+1)+
   nmat*(num_state_material+ngeom_raw)+1;
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");

 if ((im<0)||(im>=nmat))
  amrex::Error("im invalid52");

 if (ns_is_rigid(im)==0) {

  if ((elastic_time[im]>0.0)&&(elastic_viscosity[im]>0.0)) {

   int partid=0;
   while ((im_elastic_map[partid]!=im)&&(partid<im_elastic_map.size())) {
    partid++;
   }

   if (partid<im_elastic_map.size()) {

    int scomp_tensor=partid*NUM_TENSOR_TYPE;

    if (NUM_TENSOR_TYPE!=2*AMREX_SPACEDIM)
     amrex::Error("NUM_TENSOR_TYPE invalid");

     // VISCOTEN_MF will be used by NavierStokes::make_viscoelastic_heating
     //  or
     // VISCOTEN_MF will be used by NavierStokes::GetDrag
     //  or
     // VISCOTEN_MF will be used by NavierStokes::make_viscoelastic_force
     //
    getStateTensor_localMF(VISCOTEN_MF,1,scomp_tensor,NUM_TENSOR_TYPE,
     cur_time_slab);

    int rzflag=0;
    if (geom.IsRZ())
     rzflag=1;
    else if (geom.IsCartesian())
     rzflag=0;
    else if (geom.IsCYLINDRICAL())
     rzflag=3;
    else
     amrex::Error("CoordSys bust 2");

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
     //    H=H(F-1/2) or H=H(phi)

     FArrayBox& viscfab=(*localMF[CELL_VISC_MATERIAL_MF])[mfi];
     int ncomp_visc=viscfab.nComp();

     int tid_current=ns_thread();
     if ((tid_current<0)||(tid_current>=thread_class::nthreads))
      amrex::Error("tid_current invalid");
     thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

       // in: GODUNOV_3D.F90
       // viscoelastic_model==0 => (eta/lambda_mod)*visc_coef*Q
       // viscoelastic_model==2 => (eta)*visc_coef*Q
     FORT_MAKETENSOR(
      &ncomp_visc,&im, 
      xlo,dx,
      viscfab.dataPtr(),ARLIM(viscfab.loVect()),ARLIM(viscfab.hiVect()),
      tenfab.dataPtr(),ARLIM(tenfab.loVect()),ARLIM(tenfab.hiVect()),
      tilelo,tilehi,
      fablo,fabhi,&bfact,
      &elastic_viscosity[im],
      &etaS[im],
      &elastic_time[im],
      &viscoelastic_model[im],
      &polymer_factor[im],
      &rzflag,&ngrow,&nmat);
    }  // mfi  
}//omp
    ns_reconcile_d_num(54);
   } else
    amrex::Error("partid could not be found: make_viscoelastic_tensor");

  } else if ((elastic_time[im]==0.0)||
             (elastic_viscosity[im]==0.0)) {

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

 bool use_tiling=ns_tiling;

 int nmat=num_materials;
 int nden=nmat*num_state_material;

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 int ntensor=AMREX_SPACEDIM*AMREX_SPACEDIM;
 int ntensorMM=ntensor*num_materials_vel;

 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=nmat)) {
  // do nothing
 } else
  amrex::Error("num_materials_viscoelastic invalid");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int ngrow=1;  // number of grow cells for the tensor

 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(FACE_VAR_MF+dir,0,2);

 debug_ngrow(idx,0,4);
 if (localMF[idx]->nComp()!=num_materials_vel)
  amrex::Error("localMF[idx]->nComp() invalid");

 debug_ngrow(CELLTENSOR_MF,1,6);
 if (localMF[CELLTENSOR_MF]->nComp()!=ntensorMM)
  amrex::Error("localMF[CELLTENSOR_MF]->nComp() invalid");

 debug_ngrow(CELL_VISC_MATERIAL_MF,ngrow,3);

 debug_ngrow(CELL_DEN_MF,1,28); 
 debug_ngrow(CELL_DEDT_MF,1,28); 

 if (localMF[CELL_DEN_MF]->nComp()!=nmat+1)
  amrex::Error("localMF[CELL_DEN_MF]->nComp() invalid");
 if (localMF[CELL_DEDT_MF]->nComp()!=nmat+1)
  amrex::Error("localMF[CELL_DEDT_MF]->nComp() invalid");

 // 1. viscosity coefficient - 1..nmat
 // 2. viscoelastic coefficient - 1..nmat
 // 3. relaxation time - 1..nmat
 // the viscous and viscoelastic forces should both be multiplied by
 // visc_coef.  
 if (localMF[CELL_VISC_MATERIAL_MF]->nComp()<nmat)
  amrex::Error("cell_visc_material ncomp invalid");
 if (localMF[CELL_VISC_MATERIAL_MF]->nGrow()<ngrow)
  amrex::Error("cell_visc_material ngrow invalid");
 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

 int nstate=num_materials_vel*(AMREX_SPACEDIM+1)+
   nmat*(num_state_material+ngeom_raw)+1;
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");

 if ((im<0)||(im>=nmat))
  amrex::Error("im invalid53");

 if (ns_is_rigid(im)==0) {

  if ((elastic_time[im]>0.0)&&(elastic_viscosity[im]>0.0)) {

   debug_ngrow(VISCOTEN_MF,1,5);
   if (localMF[VISCOTEN_MF]->nComp()!=NUM_TENSOR_TYPE)
    amrex::Error("localMF[VISCOTEN_MF] invalid");

   int ncomp_visc=localMF[CELL_VISC_MATERIAL_MF]->nComp();
   if (ncomp_visc!=3*nmat)
    amrex::Error("cell_visc_material ncomp invalid");

   resize_levelsetLO(2,LEVELPC_MF);
   debug_ngrow(LEVELPC_MF,2,5);

   int rzflag=0;
   if (geom.IsRZ())
    rzflag=1;
   else if (geom.IsCartesian())
    rzflag=0;
   else if (geom.IsCYLINDRICAL())
    rzflag=3;
   else
    amrex::Error("CoordSys bust 2");

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
    if (tenfab.nComp()!=NUM_TENSOR_TYPE)
     amrex::Error("tenfab.nComp invalid");

    FArrayBox& DeDTinversefab=(*localMF[CELL_DEDT_MF])[mfi]; // 1/(rho cv)
    if (DeDTinversefab.nComp()!=nmat+1)
     amrex::Error("DeDTinversefab.nComp() invalid");

    FArrayBox& gradufab=(*localMF[CELLTENSOR_MF])[mfi];
    if (gradufab.nComp()!=ntensorMM)
     amrex::Error("gradufab.nComp() invalid");

    FArrayBox& heatfab=(*localMF[idx])[mfi];
    if (heatfab.nComp()!=num_materials_vel)
     amrex::Error("heatfab.nComp() invalid");

    FArrayBox& lsfab=(*localMF[LEVELPC_MF])[mfi];
    FArrayBox& viscfab=(*localMF[CELL_VISC_MATERIAL_MF])[mfi];
    if (viscfab.nComp()<nmat)
     amrex::Error("viscfab.nComp() invalid");

    FArrayBox& xface=(*localMF[FACE_VAR_MF])[mfi];
    FArrayBox& yface=(*localMF[FACE_VAR_MF+1])[mfi];
    FArrayBox& zface=(*localMF[FACE_VAR_MF+AMREX_SPACEDIM-1])[mfi];

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    FORT_TENSORHEAT(
     &elasticface_flag,
     &massface_index,
     &vofface_index,
     &ncphys,
     &ntensor,
     &ntensorMM,
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
     &dt_slab,&rzflag,&im,&nmat,&nden);
   }  // mfi  
}//omp

   ns_reconcile_d_num(55);

  } else if ((elastic_time[im]==0.0)||
             (elastic_viscosity[im]==0.0)) {

   if (viscoelastic_model[im]!=0)
    amrex::Error("viscoelastic_model[im]!=0");

  } else
   amrex::Error("elastic_time/elastic_viscosity invalid");

 } else if (ns_is_rigid(im)==1) {
  // do nothing
 } else
  amrex::Error("ns_is_rigid invalid");

}   // subroutine make_viscoelastic_heating

// TO DO: for elastic materials, the force should be calculated
//  in such a way as to preserve the reversibility property that
//  materials return to their original shape if the external
//  forces are turned off. 
void NavierStokes::make_viscoelastic_force(int im) {

 int nmat=num_materials;
 bool use_tiling=ns_tiling;

 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=nmat)) {
  // do nothing
 } else
  amrex::Error("num_materials_viscoelastic invalid");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int ngrow=1;

 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(FACE_VAR_MF+dir,0,2);

 debug_ngrow(CELL_VISC_MATERIAL_MF,ngrow,3);

 int nden=nmat*num_state_material;

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

 if (localMF[CELL_VISC_MATERIAL_MF]->nComp()<nmat)
  amrex::Error("cell_visc_material ncomp invalid");
 if (localMF[CELL_VISC_MATERIAL_MF]->nGrow()<ngrow)
  amrex::Error("cell_visc_material ngrow invalid");

 int nstate=num_materials_vel*(AMREX_SPACEDIM+1)+
   nmat*(num_state_material+ngeom_raw)+1;
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");

 if ((im<0)||(im>=nmat))
  amrex::Error("im invalid54");

 if (ns_is_rigid(im)==0) {

  if ((elastic_time[im]>0.0)&&(elastic_viscosity[im]>0.0)) {

   debug_ngrow(VISCOTEN_MF,1,5);
    // CELL_VISC_MATERIAL init in getStateVISC_ALL which
    // calls getStateVISC which calls:
    //  FORT_GETSHEAR,FORT_DERVISCOSITY, and
    //  FORT_DERTURBVISC
    //  FORT_DERVISCOSITY is in DERIVE_3D.F90
    //  a. 1..nmat           mu or etaS+etaP*(bterm**pterm)   "VISC_RAW"
    //  b. nmat+1..2*nmat    (i)   visc_coef*(VISC_RAW-etaS)/lambda or
    //                       (ii)  visc_coef*(VISC_RAW-etaS) or
    //                       (iii) visc_coef*elastic_viscosity
    //  c. n*nmat+1..3*nmat  lambda
    //  etaS=viscconst-elastic_viscosity
    //  etaP=etaP0=elastic_viscosity
   int ncomp_visc=localMF[CELL_VISC_MATERIAL_MF]->nComp();
   if (ncomp_visc!=3*nmat)
    amrex::Error("cell_visc_material ncomp invalid");

   resize_levelsetLO(2,LEVELPC_MF);
   debug_ngrow(LEVELPC_MF,2,5);
   if (localMF[LEVELPC_MF]->nComp()!=nmat*(AMREX_SPACEDIM+1))
    amrex::Error("localMF[LEVELPC_MF]->nComp()!=nmat*(AMREX_SPACEDIM+1)");

   int rzflag=0;
   if (geom.IsRZ())
    rzflag=1;
   else if (geom.IsCartesian())
    rzflag=0;
   else if (geom.IsCYLINDRICAL())
    rzflag=3;
   else
    amrex::Error("CoordSys bust 2");

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
    FArrayBox& rhoinversefab=(*localMF[CELL_DEN_MF])[mfi];

    FArrayBox& lsfab=(*localMF[LEVELPC_MF])[mfi];
    FArrayBox& viscfab=(*localMF[CELL_VISC_MATERIAL_MF])[mfi];
    if (viscfab.nComp()<nmat)
     amrex::Error("viscfab.nComp() invalid");

    FArrayBox& xface=(*localMF[FACE_VAR_MF])[mfi];
    FArrayBox& yface=(*localMF[FACE_VAR_MF+1])[mfi];
    FArrayBox& zface=(*localMF[FACE_VAR_MF+AMREX_SPACEDIM-1])[mfi];

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    FORT_TENSORFORCE(
     &elasticface_flag,
     &massface_index,
     &vofface_index,
     &ncphys,
     &nstate,
     xlo,dx,
     xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()), 
     yface.dataPtr(),ARLIM(yface.loVect()),ARLIM(yface.hiVect()), 
     zface.dataPtr(),ARLIM(zface.loVect()),ARLIM(zface.hiVect()), 
     lsfab.dataPtr(),ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
     rhoinversefab.dataPtr(),
     ARLIM(rhoinversefab.loVect()),ARLIM(rhoinversefab.hiVect()),
     S_new[mfi].dataPtr(),
     ARLIM(S_new[mfi].loVect()),ARLIM(S_new[mfi].hiVect()),
     tenfab.dataPtr(),
     ARLIM(tenfab.loVect()),ARLIM(tenfab.hiVect()),
     tilelo,tilehi,
     fablo,fabhi,&bfact,&level,
     &dt_slab,&rzflag,&im,&nmat,&nden);
   }  // mfi  
} // omp

   ns_reconcile_d_num(56);

  } else if ((elastic_time[im]==0.0)||
             (elastic_viscosity[im]==0.0)) {

   if (viscoelastic_model[im]!=0)
    amrex::Error("viscoelastic_model[im]!=0");

  } else
   amrex::Error("elastic_time/elastic_viscosity invalid");


 } else if (ns_is_rigid(im)==1) {
  // do nothing
 } else
  amrex::Error("ns_is_rigid invalid");

}   // subroutine make_viscoelastic_force

void NavierStokes::make_marangoni_force(int isweep) {

 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();

 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid make_marangoni_force");

 int nmat=num_materials;

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");
 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 resize_levelsetLO(2,LEVELPC_MF);
 debug_ngrow(LEVELPC_MF,2,5);
 if (localMF[LEVELPC_MF]->nComp()!=nmat*(AMREX_SPACEDIM+1))
  amrex::Error("localMF[LEVELPC_MF]->nComp() invalid");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(FACE_VAR_MF+dir,0,2);

 VOF_Recon_resize(2,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,2,3);
 debug_ngrow(DIST_CURV_MF,1,3);
 debug_ngrow(CELL_DEN_MF,1,5);
 if (localMF[CELL_DEN_MF]->nComp()!=nmat+1)
  amrex::Error("localMF[CELL_DEN_MF]->nComp() invalid");

  // mask=1 if not covered or if outside the domain.
  // NavierStokes::maskfiner_localMF
  // NavierStokes::maskfiner
 resize_maskfiner(2,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,2,28);
 resize_mask_nbr(2);
 debug_ngrow(MASK_NBR_MF,2,2);

 resize_metrics(2);

 MultiFab* CL_velocity;

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  debug_ngrow(CONSERVE_FLUXES_MF+dir,0,7);
  if (localMF[CONSERVE_FLUXES_MF+dir]->nComp()!=AMREX_SPACEDIM)
   amrex::Error("localMF[CONSERVE_FLUXES_MF+dir]->nComp() invalid");
 }   

 if (isweep==0) {
   // second order coarse-fine interp for LS
  getStateDist_localMF(GHOSTDIST_MF,2,cur_time_slab,16);
  getStateDen_localMF(DEN_RECON_MF,2,prev_time_slab);

  CL_velocity=getState(2,0,num_materials_vel*(AMREX_SPACEDIM+1),prev_time_slab);
  if (CL_velocity->nComp()!=AMREX_SPACEDIM+1)
   amrex::Error("CL_velocity->nComp()!=AMREX_SPACEDIM+1");

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   setVal_localMF(CONSERVE_FLUXES_MF+dir,0.0,0,AMREX_SPACEDIM,0);
   new_localMF(POTENTIAL_EDGE_MF+dir,2*AMREX_SPACEDIM,0,dir); //Ften-,Ften+
   setVal_localMF(POTENTIAL_EDGE_MF+dir,0.0,0,2*AMREX_SPACEDIM,0);
  }
 } else if (isweep==1) {
  CL_velocity=localMF[GHOSTDIST_MF];
  if (CL_velocity->nComp()!=localMF[GHOSTDIST_MF]->nComp())
   amrex::Error("CL_velocity->nComp()!=localMF[GHOSTDIST_MF]->nComp()");
 } else {
  amrex::Error("isweep invalid");
 }
 if (CL_velocity->nGrow()!=2)
  amrex::Error("CL_velocity->nGrow()!=2");

 debug_ngrow(GHOSTDIST_MF,2,5);
 if (localMF[GHOSTDIST_MF]->nComp()!=nmat*(AMREX_SPACEDIM+1))
  amrex::Error("localMF[GHOSTDIST_MF]->nComp() invalid");
 debug_ngrow(DEN_RECON_MF,2,5);
 if (localMF[DEN_RECON_MF]->nComp()!=nmat*num_state_material)
  amrex::Error("den_recon has invalid ncomp");

 int nden=nmat*num_state_material;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

 int nstate=num_materials_vel*(AMREX_SPACEDIM+1)+
   nmat*(num_state_material+ngeom_raw)+1;
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");

 int scomp_mofvars=num_materials_vel*(AMREX_SPACEDIM+1)+
  nmat*num_state_material;

 const Real* dx = geom.CellSize();

  // height function curvature
  // finite difference curvature
  // pforce
  // marangoni force
  // dir/side flag
  // im3
  // x nten
 int num_curv=nten*(5+AMREX_SPACEDIM);
 if (localMF[DIST_CURV_MF]->nComp()!=num_curv)
  amrex::Error("DIST_CURV invalid ncomp");
  
 MultiFab& Umac_new=get_new_data(Umac_Type,slab_step+1);
 MultiFab& Vmac_new=get_new_data(Umac_Type+1,slab_step+1);
 MultiFab& Wmac_new=get_new_data(Umac_Type+AMREX_SPACEDIM-1,slab_step+1);

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
  int bfact_grid=parent->blockingFactor(level);

  const Real* xlo = grid_loc[gridno].lo();

   // mask=tag if not covered by level+1 or outside the domain.
  FArrayBox& maskcov=(*localMF[MASKCOEF_MF])[mfi];

  FArrayBox& curvfab=(*localMF[DIST_CURV_MF])[mfi];
  FArrayBox& voffab=(*localMF[SLOPE_RECON_MF])[mfi];
  FArrayBox& rhoinversefab=(*localMF[CELL_DEN_MF])[mfi];
  FArrayBox& lsfab=(*localMF[LEVELPC_MF])[mfi];
  FArrayBox& lshofab=(*localMF[GHOSTDIST_MF])[mfi];
  FArrayBox& velfab=(*CL_velocity)[mfi];
  FArrayBox& denfab=(*localMF[DEN_RECON_MF])[mfi];
   // mask_nbr:
   // (1) =1 interior  =1 fine-fine ghost in domain  =0 otherwise
   // (2) =1 interior  =0 otherwise
   // (3) =1 interior+ngrow-1  =0 otherwise
   // (4) =1 interior+ngrow    =0 otherwise
  FArrayBox& masknbr=(*localMF[MASK_NBR_MF])[mfi];
  FArrayBox& areax=(*localMF[AREA_MF])[mfi];
  FArrayBox& areay=(*localMF[AREA_MF+1])[mfi];
  FArrayBox& areaz=(*localMF[AREA_MF+AMREX_SPACEDIM-1])[mfi];
  FArrayBox& xflux=(*localMF[CONSERVE_FLUXES_MF])[mfi];
  FArrayBox& yflux=(*localMF[CONSERVE_FLUXES_MF+1])[mfi];
  FArrayBox& zflux=(*localMF[CONSERVE_FLUXES_MF+AMREX_SPACEDIM-1])[mfi];
  FArrayBox& volfab=(*localMF[VOLUME_MF])[mfi];
  FArrayBox& xface=(*localMF[FACE_VAR_MF])[mfi];
  FArrayBox& yface=(*localMF[FACE_VAR_MF+1])[mfi];
  FArrayBox& zface=(*localMF[FACE_VAR_MF+AMREX_SPACEDIM-1])[mfi];
  FArrayBox& xp=(*localMF[POTENTIAL_EDGE_MF])[mfi];
  FArrayBox& yp=(*localMF[POTENTIAL_EDGE_MF+1])[mfi];
  FArrayBox& zp=(*localMF[POTENTIAL_EDGE_MF+AMREX_SPACEDIM-1])[mfi];

  Vector<int> presbc=getBCArray(State_Type,gridno,
   num_materials_vel*AMREX_SPACEDIM,1);
  Vector<int> velbc=getBCArray(State_Type,gridno,0,
   num_materials_vel*AMREX_SPACEDIM);
  Vector<int> vofbc=getBCArray(State_Type,gridno,scomp_mofvars,1);

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  // force=dt * div((I-nn^T)(sigma) delta) / rho
  // in: GODUNOV_3D.F90
  FORT_MARANGONIFORCE(
   &conservative_tension_force,
   &isweep,
   &nstate,
   &nten,
   &num_curv,
   xlo,dx,
   &facecut_index,
   &icefacecut_index,
   &curv_index,
   &pforce_index,
   &faceden_index,
   &icemask_index,
   &massface_index,
   &vofface_index,
   &ncphys,
   xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()), 
   yface.dataPtr(),ARLIM(yface.loVect()),ARLIM(yface.hiVect()), 
   zface.dataPtr(),ARLIM(zface.loVect()),ARLIM(zface.hiVect()), 
   xp.dataPtr(),ARLIM(xp.loVect()),ARLIM(xp.hiVect()),
   yp.dataPtr(),ARLIM(yp.loVect()),ARLIM(yp.hiVect()),
   zp.dataPtr(),ARLIM(zp.loVect()),ARLIM(zp.hiVect()),
   maskcov.dataPtr(),
   ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
   masknbr.dataPtr(),
   ARLIM(masknbr.loVect()),ARLIM(masknbr.hiVect()),
   volfab.dataPtr(),ARLIM(volfab.loVect()),ARLIM(volfab.hiVect()),
   areax.dataPtr(),ARLIM(areax.loVect()),ARLIM(areax.hiVect()),
   areay.dataPtr(),ARLIM(areay.loVect()),ARLIM(areay.hiVect()),
   areaz.dataPtr(),ARLIM(areaz.loVect()),ARLIM(areaz.hiVect()),
   xflux.dataPtr(),ARLIM(xflux.loVect()),ARLIM(xflux.hiVect()),
   yflux.dataPtr(),ARLIM(yflux.loVect()),ARLIM(yflux.hiVect()),
   zflux.dataPtr(),ARLIM(zflux.loVect()),ARLIM(zflux.hiVect()),
   velfab.dataPtr(),
   ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
   denfab.dataPtr(),
   ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
   lsfab.dataPtr(),ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
   lshofab.dataPtr(),
   ARLIM(lshofab.loVect()),ARLIM(lshofab.hiVect()),
   rhoinversefab.dataPtr(),
   ARLIM(rhoinversefab.loVect()),ARLIM(rhoinversefab.hiVect()),
   voffab.dataPtr(),ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),
   curvfab.dataPtr(),ARLIM(curvfab.loVect()),ARLIM(curvfab.hiVect()),
   S_new[mfi].dataPtr(),
   ARLIM(S_new[mfi].loVect()),ARLIM(S_new[mfi].hiVect()),
   Umac_new[mfi].dataPtr(),
   ARLIM(Umac_new[mfi].loVect()),ARLIM(Umac_new[mfi].hiVect()),
   Vmac_new[mfi].dataPtr(),
   ARLIM(Vmac_new[mfi].loVect()),ARLIM(Vmac_new[mfi].hiVect()),
   Wmac_new[mfi].dataPtr(),
   ARLIM(Wmac_new[mfi].loVect()),ARLIM(Wmac_new[mfi].hiVect()),
   tilelo,tilehi,
   fablo,fabhi,
   &bfact,
   &bfact_grid,
   &level,
   &finest_level,
   &dt_slab,
   &cur_time_slab,
   &visc_coef,
   &solvability_projection,  // not used
   presbc.dataPtr(),
   velbc.dataPtr(),
   vofbc.dataPtr(),
   &nmat,&nden);
 }  // mfi  
} // omp
 ns_reconcile_d_num(57);

 if (isweep==0) {
  delete CL_velocity;
  int spectral_override=0;  // always do low order average down
  int ncomp_edge=-1;
  avgDownEdge_localMF(
   CONSERVE_FLUXES_MF,
   0,ncomp_edge,0,AMREX_SPACEDIM,spectral_override,1);
  avgDownEdge_localMF(
   POTENTIAL_EDGE_MF,
   0,ncomp_edge,0,AMREX_SPACEDIM,spectral_override,2);
 } else if (isweep==1) {
  delete_localMF(GHOSTDIST_MF,1);
  delete_localMF(DEN_RECON_MF,1);
  delete_localMF(POTENTIAL_EDGE_MF,AMREX_SPACEDIM);
 } else {
  amrex::Error("isweep invalid");
 }

}   // make_marangoni_force

void NavierStokes::ns_reconcile_d_num(int caller_id) {

 thread_class::sync_tile_d_numPts();
 ParallelDescriptor::ReduceRealSum(thread_class::tile_d_numPts[0]);
 thread_class::reconcile_d_numPts(caller_id);

}

int NavierStokes::ns_thread() {
 
 int tid=0;
#ifdef _OPENMP
 tid = omp_get_thread_num();
#endif
 if ((tid<0)||(tid>=thread_class::nthreads))
  amrex::Error("tid invalid");

 return tid;
}

// add correction term to velocity and/or temperature
// called from:
//  NavierStokes::multiphase_project
//  NavierStokes::veldiffuseALL
void NavierStokes::make_SEM_delta_force(int project_option) {

 bool use_tiling=ns_tiling;

 if ((slab_step<0)||(slab_step>=ns_time_order))
  amrex::Error("slab_step invalid");

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 if ((SDC_outer_sweeps<=0)||
     (SDC_outer_sweeps>=ns_time_order))
  amrex::Error("SDC_outer_sweeps invalid");
 if (ns_time_order<=1)
  amrex::Error("ns_time_order invalid");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int nmat=num_materials;

 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(FACE_VAR_MF+dir,0,2);

 debug_ngrow(delta_MF,0,3);
 debug_ngrow(MASKSEM_MF,1,28); 
 debug_ngrow(CELL_DEN_MF,1,28); 
 debug_ngrow(CELL_DEDT_MF,1,28); 
 if (localMF[CELL_DEN_MF]->nComp()!=nmat+1)
  amrex::Error("localMF[CELL_DEN_MF]->nComp() invalid");
 if (localMF[CELL_DEDT_MF]->nComp()!=nmat+1)
  amrex::Error("localMF[CELL_DEDT_MF]->nComp() invalid");

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

 int nstate=num_materials_vel*(AMREX_SPACEDIM+1)+
   nmat*(num_state_material+ngeom_raw)+1;
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");

 const Real* dx = geom.CellSize();
 int bfact=parent->Space_blockingFactor(level);

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

  // I-scheme,thermal conduction,viscosity,div(up),gp,-force
  FArrayBox& deltafab=(*localMF[delta_MF])[mfi];
  int deltacomp=0;
  if (project_option==3) { // viscosity
   deltacomp=slab_step*nstate_SDC+nfluxSEM+1;
  } else if (project_option==4) { // -momentum force at t^n+1
   deltacomp=slab_step*nstate_SDC+nfluxSEM+1+AMREX_SPACEDIM+1+AMREX_SPACEDIM;
  } else if (project_option==2) { // thermal conduction
   deltacomp=slab_step*nstate_SDC+nfluxSEM;
  } else if (project_option==0) { // div up and gp 
     // advection, thermal conduction,viscosity,div(up),gp
   deltacomp=slab_step*nstate_SDC+nfluxSEM+1+AMREX_SPACEDIM;
  } else
   amrex::Error("project_option invalid4");

  FArrayBox& rhoinversefab=(*localMF[CELL_DEN_MF])[mfi];
  FArrayBox& DeDTinversefab=(*localMF[CELL_DEDT_MF])[mfi]; // 1/(rho cv)
  FArrayBox& maskSEMfab=(*localMF[MASKSEM_MF])[mfi];

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // in: GODUNOV_3D.F90
  FORT_SEMDELTAFORCE(
   &nstate,
   &nfluxSEM,
   &nstate_SDC,
   &nmat,
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
   &dt_slab);
 }  // mfi  
} // omp
 ns_reconcile_d_num(58);

  // pressure gradient at faces.
 if (project_option==0) {

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

      // faceden_index=1/rho
    FORT_SEMDELTAFORCE_FACE(
     &dir,
     &faceden_index,
     &ncphys,
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
     &dt_slab);
   }  // mfi  
} // omp
   ns_reconcile_d_num(59);
  } // dir=0..sdim-1

 } else if ((project_option==2)|| // thermal conduction
	    (project_option==3)|| // viscosity
	    (project_option==4)) { // -momentum force at t^n+1
	 // do nothing
 } else {
  amrex::Error("project_option invalid");
 } 

}   // subroutine make_SEM_delta_force


// MEHDI VAHAB HEAT SOURCE
// called from veldiffuseALL in NavierStokes3.cpp
void NavierStokes::make_heat_source() {

 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();
 int nmat=num_materials;

 if ((slab_step<0)||(slab_step>=ns_time_order))
  amrex::Error("slab_step invalid");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(FACE_VAR_MF+dir,0,2);

 VOF_Recon_resize(1,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,1,3);
 resize_levelsetLO(2,LEVELPC_MF);
 debug_ngrow(LEVELPC_MF,2,5);
 if (localMF[LEVELPC_MF]->nComp()!=nmat*(AMREX_SPACEDIM+1))
  amrex::Error("localMF[LEVELPC_MF]->nComp() invalid");

 resize_metrics(1);  
 debug_ngrow(VOLUME_MF,1,28); 

 debug_ngrow(CELL_DEN_MF,1,5);
 if (localMF[CELL_DEN_MF]->nComp()!=nmat+1)
  amrex::Error("localMF[CELL_DEN_MF]->nComp() invalid");
 debug_ngrow(CELL_DEDT_MF,1,5);
 if (localMF[CELL_DEDT_MF]->nComp()!=nmat+1)
  amrex::Error("localMF[CELL_DEDT_MF]->nComp() invalid");

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

 int nstate=num_materials_vel*(AMREX_SPACEDIM+1)+
   nmat*(num_state_material+ngeom_raw)+1;
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");

 int dencomp=num_materials_vel*(AMREX_SPACEDIM+1);
 int nden=nmat*num_state_material;

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

  FArrayBox& rhoinversefab=(*localMF[CELL_DEN_MF])[mfi];
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
   // in: GODUNOV_3D.F90
  FORT_HEATSOURCE(
   &nstate,
   &nmat,
   &nden,
   xlo,dx,
   &temperature_source,
   temperature_source_cen.dataPtr(),
   temperature_source_rad.dataPtr(),
   rhoinversefab.dataPtr(),
   ARLIM(rhoinversefab.loVect()),ARLIM(rhoinversefab.hiVect()),
   DeDTinversefab.dataPtr(),
   ARLIM(DeDTinversefab.loVect()),ARLIM(DeDTinversefab.hiVect()),
   snewfab.dataPtr(dencomp),
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
   &dt_slab,  // time step within a slab if SDC, otherwise dt if not SDC.
   &prev_time_slab);
 }  // mfi  
} // omp
 ns_reconcile_d_num(60);

}   // subroutine make_heat_source



void NavierStokes::add_perturbation() {

 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();

 if (slab_step!=ns_time_order-1)
  amrex::Error("slab_step invalid");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 MultiFab& LS_new = get_new_data(LS_Type,slab_step+1);

 int nmat=num_materials;

 int nstate=num_materials_vel*(AMREX_SPACEDIM+1)+
   nmat*(num_state_material+ngeom_raw)+1;
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

   FORT_ADDNOISE(
    &dir,
    &angular_velocity,
    &perturbation_mode,
    &perturbation_eps_temp,
    &perturbation_eps_vel,
    &nstate,
    &nmat,
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
 ns_reconcile_d_num(61);

}   // subroutine add_perturbation

// called from: update_SEM_forces
// update_SEM_forces is called from: update_SEM_forcesALL
// update_SEM_forcesALL is called from:
//   NavierStokes::do_the_advance
//   NavierStokes::veldiffuseALL
void NavierStokes::update_SEM_delta_force(
 int project_option,
 int idx_gp,int idx_gpmac,int idx_div,
 int update_spectral,int update_stable,
 int nsolve) {

 bool use_tiling=ns_tiling;

 int nmat=num_materials;

 if ((ns_time_order>=2)&&(ns_time_order<=32)) {
  // do nothing
 } else
  amrex::Error("ns_time_order invalid");

 if ((enable_spectral==1)||(enable_spectral==3)) {
  // do nothing
 } else {
  std::cout << "ns_time_order= " << ns_time_order << '\n';
  std::cout << "project_option= " << project_option << '\n';
  std::cout << "update_spectral= " << update_spectral << '\n';
  std::cout << "update_stable= " << update_stable << '\n';
  amrex::Error("enable_spectral invalid in update_SEM_delta_force");
 }

 if ((nsolve!=1)&&(nsolve!=AMREX_SPACEDIM))
  amrex::Error("nsolve invalid37");

 int num_materials_face=num_materials_vel;

 if (project_option==0) {
  if (num_materials_face!=1)
   amrex::Error("num_materials_face invalid");
 } else if (project_option==2) { // thermal diffusion
  num_materials_face=num_materials_scalar_solve;
 } else if (project_option==3) { // viscosity
  if (num_materials_face!=1)
   amrex::Error("num_materials_face invalid");
 } else if (project_option==4) { // -momentum force
  if (num_materials_face!=1)
   amrex::Error("num_materials_face invalid");
 } else
  amrex::Error("project_option invalid5");

 int nsolveMM=nsolve*num_materials_face;

 int nsolveMM_FACE=nsolve*num_materials_face;
 if (num_materials_face==1) {
  // do nothing
 } else if (num_materials_face==nmat) {
  nsolveMM_FACE*=2;
 } else
  amrex::Error("num_materials_face invalid");

 if (update_stable==1) {
  if ((slab_step<0)||(slab_step>=ns_time_order))
   amrex::Error("slab_step invalid");
 } else if (update_stable==0) {
  // do nothing
 } else
  amrex::Error("update_stable invalid update sem delta force");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 debug_ngrow(idx_div,0,3);
 debug_ngrow(idx_gp,0,3);

 int idx_hoop=idx_div;

 if (project_option==0) { // grad p, div(up)
  idx_hoop=idx_div;
  if (nsolve!=1)
   amrex::Error("nsolve invalid");
  if (localMF[idx_gp]->nComp()!=nsolveMM*AMREX_SPACEDIM)
   amrex::Error("localMF[idx_gp]->nComp() invalid");
  if (localMF[idx_div]->nComp()!=nsolveMM)
   amrex::Error("localMF[idx_div]->nComp() invalid");
 } else if (project_option==2) { // -div(k grad T)-THERMAL_FORCE_MF
  idx_hoop=THERMAL_FORCE_MF;
  if (nsolve!=1)
   amrex::Error("nsolve invalid");
  if (localMF[idx_div]->nComp()!=nsolveMM) {
   std::cout << "project_option = " << project_option << '\n';
   std::cout << "idx_div = " << idx_div << '\n';
   std::cout << "nsolveMM = " << nsolveMM << '\n';
   std::cout << "nsolve = " << nsolve << '\n';
   std::cout << "localMF ncomp= " <<
     localMF[idx_div]->nComp() << '\n';
   amrex::Error("localMF[idx_div]->nComp() invalid");
  }
  if (localMF[idx_hoop]->nComp()!=nsolveMM) {
   std::cout << "project_option = " << project_option << '\n';
   std::cout << "idx_hoop = " << idx_hoop << '\n';
   std::cout << "nsolveMM = " << nsolveMM << '\n';
   std::cout << "nsolve = " << nsolve << '\n';
   std::cout << "localMF ncomp= " << 
     localMF[idx_hoop]->nComp() << '\n';
   amrex::Error("localMF[idx_hoop]->nComp() invalid");
  }
 } else if (project_option==3) { // -div(2 mu D)-HOOP_FORCE_MARK_MF
  idx_hoop=HOOP_FORCE_MARK_MF;
  if (nsolve!=AMREX_SPACEDIM)
   amrex::Error("nsolve invalid");
  if (localMF[idx_div]->nComp()!=nsolveMM) {
   std::cout << "project_option = " << project_option << '\n';
   std::cout << "idx_div = " << idx_div << '\n';
   std::cout << "nsolveMM = " << nsolveMM << '\n';
   std::cout << "nsolve = " << nsolve << '\n';
   std::cout << "localMF ncomp= " <<
     localMF[idx_div]->nComp() << '\n';
   amrex::Error("localMF[idx_div]->nComp() invalid");
  }
  if (localMF[idx_hoop]->nComp()!=nsolveMM) {
   std::cout << "project_option = " << project_option << '\n';
   std::cout << "idx_hoop = " << idx_hoop << '\n';
   std::cout << "nsolveMM = " << nsolveMM << '\n';
   std::cout << "nsolve = " << nsolve << '\n';
   std::cout << "localMF ncomp= " <<
     localMF[idx_hoop]->nComp() << '\n';
   amrex::Error("localMF[idx_hoop]->nComp() invalid");
  }
 } else if (project_option==4) { // -momentum force
  idx_hoop=idx_div;
  if (idx_div!=idx_gp)
   amrex::Error("expecting idx_div==idx_gp");
  if (nsolve!=AMREX_SPACEDIM)
   amrex::Error("nsolve invalid");
  if (localMF[idx_hoop]->nComp()!=nsolveMM)
   amrex::Error("localMF[idx_hoop]->nComp() invalid");
 } else
  amrex::Error("project_option invalid6");

 debug_ngrow(idx_hoop,0,3);
 if (localMF[idx_hoop]->nComp()!=localMF[idx_div]->nComp())
  amrex::Error("localMF[idx_hoop]->nComp() invalid");

 if ((project_option==0)||
     (project_option==2)||
     (project_option==3)) {

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   debug_ngrow(idx_gpmac+dir,0,3);

   if (project_option==0) {
    if (localMF[idx_gpmac]->nComp()!=nsolveMM_FACE)
     amrex::Error("localMF[idx_gpmac]->nComp() invalid");
   } else if ((project_option==2)||
              (project_option==3)) {
    // do nothing
   } else
    amrex::Error("project_option invalid7");
 
   debug_ngrow(spectralF_GP_MF+dir,0,3);
   debug_ngrow(stableF_GP_MF+dir,0,3);
   debug_ngrow(delta_GP_MF+dir,0,3);
  } // dir=0..sdim-1

 } else if (project_option==4) { // -momentum force
  // check nothing
 } else
  amrex::Error("project_option invalid8");

 debug_ngrow(spectralF_MF,0,3);
 debug_ngrow(stableF_MF,0,3);
 debug_ngrow(delta_MF,0,3);

 debug_ngrow(MASKSEM_MF,1,28); 

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

 int nstate=num_materials_vel*(AMREX_SPACEDIM+1)+
   nmat*(num_state_material+ngeom_raw)+1;
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");

 const Real* dx = geom.CellSize();
 int bfact=parent->Space_blockingFactor(level);

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

  FArrayBox& gpfab=(*localMF[idx_gp])[mfi];
  FArrayBox& divfab=(*localMF[idx_div])[mfi];
  FArrayBox& hoopfab=(*localMF[idx_hoop])[mfi];
  FArrayBox& HOfab=(*localMF[spectralF_MF])[mfi];
  FArrayBox& LOfab=(*localMF[stableF_MF])[mfi];
  FArrayBox& maskSEMfab=(*localMF[MASKSEM_MF])[mfi];

  int LOcomp=0;
  int HOcomp=0;
  if (slab_step==-1)
   HOcomp=0;
  else if ((slab_step>=0)&&(slab_step<ns_time_order))
   HOcomp=nstate_SDC*(slab_step+1);
  else
   amrex::Error("slab_step invalid");

  if (update_stable==1) {
   if ((slab_step>=0)&&(slab_step<ns_time_order))
    LOcomp=nstate_SDC*slab_step;
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
 
   // in: GODUNOV_3D.F90
  FORT_UPDATESEMFORCE(
   &ns_time_order,
   &slab_step,
   &nsolve,
   &update_spectral,
   &update_stable,
   &nstate,
   &nfluxSEM,
   &nstate_SDC,
   &nmat,
   &project_option,
   xlo,dx,
   gpfab.dataPtr(),
   ARLIM(gpfab.loVect()),ARLIM(gpfab.hiVect()),
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
   &dt_slab);
 }  // mfi  
} // omp
 ns_reconcile_d_num(62);

 if (project_option==0) {

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

    FArrayBox& gpfab=(*localMF[idx_gpmac+dir])[mfi];
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

    FORT_UPDATESEMFORCE_FACE(
     &project_option,
     &num_materials_face,
     &nsolveMM_FACE,
     &ns_time_order,
     &dir,
     &slab_step,
     &update_spectral,
     &update_stable,
     &nmat,
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
     &dt_slab);
   }  // mfi  
} // omp
   ns_reconcile_d_num(63);

  } // dir=0..sdim-1

 } else if ((project_option==2)||
            (project_option==3)||
            (project_option==4)) {
  // do nothing
 } else {
  amrex::Error("project_option invalid9");
 } 

} // subroutine update_SEM_delta_force

// called from:
//  NavierStokes::tensor_advection_updateALL()  (NavierStokes3.cpp)
void NavierStokes::tensor_advection_update() {

 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();

 int ngrow_zero=0;
 int nmat=num_materials;

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=nmat)) {
  // do nothing
 } else
  amrex::Error("num_materials_viscoelastic invalid");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(FACE_VAR_MF+dir,0,2);

 resize_levelsetLO(2,LEVELPC_MF);
 debug_ngrow(LEVELPC_MF,2,8);
 if (localMF[LEVELPC_MF]->nComp()!=nmat*(AMREX_SPACEDIM+1))
  amrex::Error("(localMF[LEVELPC_MF]->nComp()!=nmat*(AMREX_SPACEDIM+1))");

 VOF_Recon_resize(1,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,1,9);
 debug_ngrow(CELLTENSOR_MF,1,9);

 MultiFab& Tensor_new=get_new_data(Tensor_Type,slab_step+1);

 int rzflag=0;
 if (geom.IsRZ())
  rzflag=1;
 else if (geom.IsCartesian())
  rzflag=0;
 else if (geom.IsCYLINDRICAL())
  rzflag=3;
 else
  amrex::Error("CoordSys bust 3");

 const Real* dx = geom.CellSize();

 for (int im=0;im<nmat;im++) {

  if (ns_is_rigid(im)==0) {

   if ((elastic_time[im]>0.0)&&(elastic_viscosity[im]>0.0)) {

    int partid=0;
    while ((im_elastic_map[partid]!=im)&&(partid<im_elastic_map.size())) {
     partid++;
    }

    if (partid<im_elastic_map.size()) {

     int scomp_tensor=partid*NUM_TENSOR_TYPE;

     int ncomp_visc=localMF[CELL_VISC_MATERIAL_MF]->nComp();
     if (ncomp_visc!=3*nmat)
      amrex::Error("cell_visc_material ncomp invalid");

     MultiFab* tensor_source_mf=
      getStateTensor(0,scomp_tensor,NUM_TENSOR_TYPE,cur_time_slab);

     int scomp_xdisplace=num_materials_viscoelastic*NUM_TENSOR_TYPE;
     MultiFab* xdisplace_mf=getStateTensor(1,scomp_xdisplace,AMREX_SPACEDIM,
       cur_time_slab);

     MultiFab* velmf=getState(1,0,AMREX_SPACEDIM,cur_time_slab);
   
     MultiFab* tendata_mf=new MultiFab(grids,dmap,20,ngrow_zero,
		  MFInfo().SetTag("tendata_mf"),FArrayBoxFactory());

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

      FArrayBox& cellten=(*localMF[CELLTENSOR_MF])[mfi];
      int ntensor=AMREX_SPACEDIM*AMREX_SPACEDIM;
      int ntensorMM=ntensor*num_materials_vel;

      if (cellten.nComp()!=ntensorMM)
       amrex::Error("cellten invalid ncomp");

      FArrayBox& voffab=(*localMF[SLOPE_RECON_MF])[mfi];
      FArrayBox& viscfab=(*localMF[CELL_VISC_MATERIAL_MF])[mfi];
      if (viscfab.nComp()<nmat)
       amrex::Error("viscfab.nComp() invalid");

      FArrayBox& velfab=(*velmf)[mfi];
      FArrayBox& tendata=(*tendata_mf)[mfi];

      Vector<int> velbc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);

      // get |grad U|,D,grad U 
      int iproject=0;
      int onlyscalar=0; 

      int tid_current=ns_thread();
      if ((tid_current<0)||(tid_current>=thread_class::nthreads))
       amrex::Error("tid_current invalid");
      thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

      // 0<=im<=nmat-1
      FORT_GETSHEAR(
       &im,
       &ntensor,
       cellten.dataPtr(),
       ARLIM(cellten.loVect()),ARLIM(cellten.hiVect()),
       voffab.dataPtr(),ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),
       velfab.dataPtr(),
       ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
       dx,xlo,
       tendata.dataPtr(),
       ARLIM(tendata.loVect()),ARLIM(tendata.hiVect()),
       &iproject,&onlyscalar,
       &cur_time_slab,
       tilelo,tilehi,
       fablo,fabhi,
       &bfact, 
       &level, 
       velbc.dataPtr(),
       &ngrow_zero, // 0
       &nmat);
     } // mfi  
} // omp
     ns_reconcile_d_num(64);

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

      FArrayBox& voffab=(*localMF[SLOPE_RECON_MF])[mfi];
      FArrayBox& viscfab=(*localMF[CELL_VISC_MATERIAL_MF])[mfi];
      FArrayBox& velfab=(*velmf)[mfi];
      FArrayBox& tensor_new_fab=Tensor_new[mfi];
      FArrayBox& tensor_source_mf_fab=(*tensor_source_mf)[mfi];
      FArrayBox& xdisplace_mf_fab=(*xdisplace_mf)[mfi];
      FArrayBox& tendata=(*tendata_mf)[mfi];

      Vector<int> velbc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);

      int transposegradu=0;

      int tid_current=ns_thread();
      if ((tid_current<0)||(tid_current>=thread_class::nthreads))
       amrex::Error("tid_current invalid");
      thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

        // in: GODUNOV_3D.F90
	// last step in this routine: (Q^n+1-Q^n)/dt = -Q^n+1/lambda
	// if viscoelastic_model==2, then modtime=elastic_time
      FORT_UPDATETENSOR(
       &level,
       &finest_level,
       &nmat,&im,
       &ncomp_visc,
       voffab.dataPtr(),ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),
       viscfab.dataPtr(),ARLIM(viscfab.loVect()),ARLIM(viscfab.hiVect()),
       tendata.dataPtr(),ARLIM(tendata.loVect()),ARLIM(tendata.hiVect()),
       dx,xlo,
       velfab.dataPtr(),
       ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
       tensor_new_fab.dataPtr(scomp_tensor),
       ARLIM(tensor_new_fab.loVect()),ARLIM(tensor_new_fab.hiVect()),
       tensor_source_mf_fab.dataPtr(),
       ARLIM(tensor_source_mf_fab.loVect()),
       ARLIM(tensor_source_mf_fab.hiVect()),
       xdisplace_mf_fab.dataPtr(),
       ARLIM(xdisplace_mf_fab.loVect()),
       ARLIM(xdisplace_mf_fab.hiVect()),
       tilelo,tilehi,
       fablo,fabhi,
       &bfact, 
       &dt_slab,
       &elastic_time[im],
       &viscoelastic_model[im],
       &polymer_factor[im],
       &rzflag,velbc.dataPtr(),&transposegradu);
     }  // mfi
} // omp
     ns_reconcile_d_num(65);

     if ((AMREX_SPACEDIM==2)&&(rzflag==1)) {

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

       FArrayBox& voffab=(*localMF[SLOPE_RECON_MF])[mfi];
       FArrayBox& tensor_new_fab=Tensor_new[mfi];

       int tid_current=ns_thread();
       if ((tid_current<0)||(tid_current>=thread_class::nthreads))
        amrex::Error("tid_current invalid");
       thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

       FORT_FIX_HOOP_TENSOR(
        &level,
        &finest_level,
        &nmat,&im,
        voffab.dataPtr(),ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),
        dx,xlo,
        tensor_new_fab.dataPtr(scomp_tensor),
        ARLIM(tensor_new_fab.loVect()),ARLIM(tensor_new_fab.hiVect()),
        tilelo,tilehi,
        fablo,fabhi,
        &bfact, 
        &rzflag);
      }  // mfi
} // omp
      ns_reconcile_d_num(66);

     } else if ((AMREX_SPACEDIM==3)||
                (rzflag==0)||
                (rzflag==3)) {
      // do nothing
     } else
      amrex::Error("bl_spacedim or rzflag invalid");
      
     delete tendata_mf;
     delete tensor_source_mf;
     delete xdisplace_mf;
     delete velmf;
    } else
     amrex::Error("partid could not be found: tensor_advection_update");

   } else if ((elastic_time[im]==0.0)||(elastic_viscosity[im]==0.0)) {
    if (viscoelastic_model[im]!=0)
     amrex::Error("viscoelastic_model[im]!=0");
   } else
    amrex::Error("viscoelastic parameter bust");

  } else if (ns_is_rigid(im)==1) {

   // do nothing

  } else
   amrex::Error("ns_is_rigid invalid");

 } // im=0..nmat-1


}   // subroutine tensor_advection_update

// extrapolate where the volume fraction is less than 1/2.  
// note: the gradients used to update the tensor are only valid in cells
// where LS(im_viscoelastic)>=0.
void NavierStokes::tensor_extrapolate() {

 int finest_level=parent->finestLevel();

 int nmat=num_materials;

 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=nmat)) {
  // do nothing
 } else
  amrex::Error("num_materials_viscoelastic invalid");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int ngrow_extrap=2;

 VOF_Recon_resize(ngrow_extrap,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,ngrow_extrap,13);

 MultiFab& Tensor_new=get_new_data(Tensor_Type,slab_step+1);

 const Real* dx = geom.CellSize();

 for (int im=0;im<nmat;im++) {

  if (ns_is_rigid(im)==0) {

   if ((elastic_time[im]>0.0)&&(elastic_viscosity[im]>0.0)) {

    int partid=0;
    while ((im_elastic_map[partid]!=im)&&(partid<im_elastic_map.size())) {
     partid++;
    }

    if (partid<im_elastic_map.size()) {

     int scomp_tensor=partid*NUM_TENSOR_TYPE;

     MultiFab* tensor_source_mf= getStateTensor(ngrow_extrap,scomp_tensor,
      NUM_TENSOR_TYPE,cur_time_slab);
   
     //use_tiling==false 
     //for future: have a separate FAB for each thread?
     if (thread_class::nthreads<1)
      amrex::Error("thread_class::nthreads invalid");
     thread_class::init_d_numPts(tensor_source_mf->boxArray().d_numPts());
 
#ifdef _OPENMP
#pragma omp parallel
#endif
{
     for (MFIter mfi(*tensor_source_mf,false); mfi.isValid(); ++mfi) {

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
      FArrayBox& tensor_new_fab=Tensor_new[mfi];
      FArrayBox& tensor_source_mf_fab=(*tensor_source_mf)[mfi];

      int tid_current=ns_thread();
      if ((tid_current<0)||(tid_current>=thread_class::nthreads))
       amrex::Error("tid_current invalid");
      thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

      FORT_EXTRAPTENSOR(
       &level,
       &finest_level,
       &nmat,&im,
       &ngrow_extrap,
       voffab.dataPtr(),ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),
       dx,xlo,
       tensor_new_fab.dataPtr(scomp_tensor),
       ARLIM(tensor_new_fab.loVect()),ARLIM(tensor_new_fab.hiVect()),
       tensor_source_mf_fab.dataPtr(),
       ARLIM(tensor_source_mf_fab.loVect()),
       ARLIM(tensor_source_mf_fab.hiVect()),
       tilelo,tilehi,
       fablo,fabhi,&bfact);
     }  // mfi
} //omp
     ns_reconcile_d_num(67);

     delete tensor_source_mf;

    } else
     amrex::Error("partid could not be found: tensor_extrapolate");

   } else if ((elastic_time[im]==0.0)||(elastic_viscosity[im]==0.0)) {

    if (viscoelastic_model[im]!=0)
     amrex::Error("viscoelastic_model[im]!=0");

   } else
    amrex::Error("viscoelastic parameter bust");
  } else if (ns_is_rigid(im)==1) {
   // do nothing
  } else
   amrex::Error("ns_is_rigid bust");
 } // im

}   // subroutine tensor_extrapolate


void 
NavierStokes::correct_density() {

 bool use_tiling=ns_tiling;

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");
 if (dt_slab<=0.0)
  amrex::Error("dt_slab invalid3");

 int nmat=num_materials;

 int non_conservative_density=0;

 for (int im=0;im<nmat;im++) {
  if ((DrhoDT[im]!=0.0)&&(override_density[im]==0))
   amrex::Error("DrhoDT mismatch"); 
  if ((DrhoDz[im]!=0.0)&&(override_density[im]==0))
   amrex::Error("DrhoDz mismatch"); 
  if (override_density[im]==1) { // rho=rho(T,Y,z)
   non_conservative_density=1;
  } else if (override_density[im]==0) {
   // do nothing

   // P_hydro=P_hydro(rho(T,Y,z)) (Boussinesq like approximation)
  } else if (override_density[im]==2) { 
   // do nothing
  } else
   amrex::Error("override density invalid");
 } // im

 if (non_conservative_density==0) {
  // do nothing
 } else if (non_conservative_density==1) {

   int finest_level=parent->finestLevel();

   resize_metrics(1);
   resize_maskfiner(1,MASKCOEF_MF);
   resize_mask_nbr(1);

   debug_ngrow(VOLUME_MF,1,28); 
   debug_ngrow(MASKCOEF_MF,1,28); 
   debug_ngrow(MASK_NBR_MF,1,28); 
   debug_ngrow(MASKSEM_MF,1,28); 
   MultiFab& S_new=get_new_data(State_Type,slab_step+1);
   const Real* dx = geom.CellSize();
 
   int scomp_pres=num_materials_vel*AMREX_SPACEDIM;

   Real gravity_normalized=std::abs(gravity);
   if (invert_gravity==1)
    gravity_normalized=-gravity_normalized;
   else if (invert_gravity==0) {
    // do nothing
   } else
    amrex::Error("invert_gravity invalid");

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
    Vector<int> presbc=getBCArray(State_Type,gridno,scomp_pres,1);
    FArrayBox& maskfab=(*localMF[MASKCOEF_MF])[mfi];
    FArrayBox& masknbrfab=(*localMF[MASK_NBR_MF])[mfi];
    FArrayBox& volfab=(*localMF[VOLUME_MF])[mfi];
    int dencomp=num_materials_vel*(AMREX_SPACEDIM+1);

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     // in: GODUNOV_3D.F90
     // if override_density[im]==1, then rho_im=rho(T,Y,z) 
    FORT_DENCOR(
      spec_material_id_AMBIENT.dataPtr(),
      species_evaporation_density.dataPtr(),
      presbc.dataPtr(),
      tilelo,tilehi,
      fablo,fabhi,&bfact,
      &dt_slab,
      maskfab.dataPtr(),
      ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
      masknbrfab.dataPtr(),
      ARLIM(masknbrfab.loVect()),ARLIM(masknbrfab.hiVect()),
      volfab.dataPtr(),ARLIM(volfab.loVect()),ARLIM(volfab.hiVect()),
      S_new[mfi].dataPtr(dencomp),
      ARLIM(S_new[mfi].loVect()),ARLIM(S_new[mfi].hiVect()),
      xlo,dx,
      &gravity_normalized,
      DrhoDT.dataPtr(),
      override_density.dataPtr(),
      &nmat,&level,&finest_level);
   }  // mfi
}  // omp
   ns_reconcile_d_num(68);
  
 } else
   amrex::Error("non_conservative_density invalid");

} // subroutine correct_density

// check that fine grids align with coarse elements and that
// fine grid dimensions are perfectly divisible by the fine order.
// fine grid dimensions are perfectly divisible by the blocking factor.
void NavierStokes::check_grid_places() {

 int bfact_SEM_coarse=0;
 int bfact_SEM=parent->Space_blockingFactor(level);
 int bfact_grid=parent->blockingFactor(level);
 if (bfact_grid<4)
  amrex::Error("we must have blocking factor at least 4");

 int bfact_fine_min=((bfact_SEM<4) ? 4 : bfact_SEM);
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

 if (bfact_fine_min<4)
  amrex::Error("bfact_fine_min<4");

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

} // subroutine check_grid_places


void 
NavierStokes::prepare_mask_nbr(int ngrow) {

 if (ngrow_distance!=4)
  amrex::Error("ngrow_distance invalid");
 if (curv_stencil_height!=ngrow_distance)
  amrex::Error("curv_stencil_height invalid");

 if ((ngrow<1)||(ngrow>ngrow_distance))
  amrex::Error("ngrow invalid");

 if (localMF_grow[MASK_NBR_MF]>=0) {
  delete_localMF(MASK_NBR_MF,1);
 }

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

} // subroutine prepare_mask_nbr()

// unsplit_displacement=0 if called from split_scalar_advection
// unsplit_displacement=1 if called from unsplit_scalar_advection
void 
NavierStokes::prepare_displacement(int mac_grow,int unsplit_displacement) {
 
 bool use_tiling=ns_tiling;

 if (divu_outer_sweeps==0)
  vel_time_slab=prev_time_slab;
 else if (divu_outer_sweeps>0)
  vel_time_slab=cur_time_slab;
 else
  amrex::Error("divu_outer_sweeps invalid");

 int mac_grow_expect=1;
 if (face_flag==0) {
  // do nothing
 } else if (face_flag==1) {
  mac_grow_expect++;
 } else
  amrex::Error("face_flag invalid 3");

 if (mac_grow!=mac_grow_expect)
  amrex::Error("mac_grow invalid in prepare_displacement");

 int finest_level=parent->finestLevel();
 int nmat=num_materials;

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 int nsolveMM_FACE=num_materials_vel;

 getState_localMF(CELL_VELOCITY_MF,mac_grow,0,
  num_materials_vel*AMREX_SPACEDIM,vel_time_slab);

 for (int normdir=0;normdir<AMREX_SPACEDIM;normdir++) {

   // mac_grow+1 for finding slopes
  MultiFab* temp_mac_velocity=getStateMAC(mac_grow+1,normdir,
   0,nsolveMM_FACE,vel_time_slab); 

   // mac_grow+1 for finding slopes
   // MAC_VELOCITY_MF deleted towards the end of 
   //   NavierStokes::nonlinear_advection
   // component 1: the velocity
   // component(s) 2..sdim+1: derivatives
  new_localMF(MAC_VELOCITY_MF+normdir,nsolveMM_FACE*(AMREX_SPACEDIM+1),
		  mac_grow+1,normdir);

  const Real* dx = geom.CellSize();
  MultiFab& S_new=get_new_data(State_Type,slab_step+1);

  // first sweep: 
  // 1. multiply velocity by dt.  
  // 2. adjust velocity if RZ.  
  // 3. override velocity if it is a passive advection problem.
  // 4. copy into mac_velocity
  // 5. repeat for cell_velocity
  // second sweep:
  // 1. calculate mac velocity slopes.

  for (int isweep=0;isweep<=1;isweep++) {

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

    FArrayBox& unodetemp=(*temp_mac_velocity)[mfi]; // macgrow+1
    FArrayBox& unode=(*localMF[MAC_VELOCITY_MF+normdir])[mfi]; // macgrow+1
    FArrayBox& ucell=(*localMF[CELL_VELOCITY_MF])[mfi];

    prescribed_vel_time_slab=0.5*(prev_time_slab+cur_time_slab);

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     // in: GODUNOV_3D.F90
    FORT_VELMAC_OVERRIDE(
     &isweep,
     &unsplit_displacement,
     &nsolveMM_FACE,
     &nmat,
     tilelo,tilehi,
     fablo,fabhi,&bfact,
     velbc.dataPtr(),
     &dt_slab,
     &prev_time_slab,
     &prescribed_vel_time_slab,
     &vel_time_slab,
     &dir_absolute_direct_split,
     &normdir,
     unodetemp.dataPtr(),ARLIM(unodetemp.loVect()),ARLIM(unodetemp.hiVect()),
     unode.dataPtr(),ARLIM(unode.loVect()),ARLIM(unode.hiVect()),
     ucell.dataPtr(),ARLIM(ucell.loVect()),ARLIM(ucell.hiVect()),
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

   ns_reconcile_d_num(69);

  } // isweep=0..1

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

 if (cancel_advection==0) {
  // do nothing
 } else if (cancel_advection==1) {
  localMF[CELL_VELOCITY_MF]->setVal(0.0);
 } else {
  amrex::Error("cancel_advection invalid");
 }

 localMF[CELL_VELOCITY_MF]->FillBoundary(geom.periodicity());

}  // prepare_displacement

void
NavierStokes::level_phase_change_rate(Vector<blobclass> blobdata,
	int color_count,int nucleation_flag) {

 Real problo[AMREX_SPACEDIM];
 Real probhi[AMREX_SPACEDIM];
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  problo[dir]=geom.ProbLo(dir);
  probhi[dir]=geom.ProbHi(dir);
 }

 int rz_flag=0;
 if (geom.IsRZ())
  rz_flag=1;
 else if (geom.IsCartesian())
  rz_flag=0;
 else if (geom.IsCYLINDRICAL())
  rz_flag=3;
 else
  amrex::Error("CoordSys bust 1");

 bool use_tiling=ns_tiling;
 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid level_phase_change_rate");

 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;
 int nden=nmat*num_state_material;
 int ncomp_per_burning=AMREX_SPACEDIM;
 int ncomp_per_tsat=2;
 int nburning=nten*(ncomp_per_burning+1);
 int ntsat=nten*(ncomp_per_tsat+1);
 int nstate=num_materials_vel*(AMREX_SPACEDIM+1)+
  nmat*(num_state_material+ngeom_raw)+1;

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");
 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 int n_normal=(nmat+nten)*(AMREX_SPACEDIM+1);

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
 getStateDen_localMF(DEN_RECON_MF,normal_probe_size+3,cur_time_slab);

 if (localMF[DEN_RECON_MF]->nComp()!=nden)
  amrex::Error("DEN_RECON_MF invalid ncomp");

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

 } else if (nucleation_flag==1) {
  if (color_count!=1)
   amrex::Error("color_count!=1");
 } else
  amrex::Error("nucleation_flag invalid");

 MultiFab* presmf=
  getState(normal_probe_size+3,num_materials_vel*AMREX_SPACEDIM,
           1,cur_time_slab);

 MultiFab* pres_eos_mf=derive_EOS_pressure();
 if (pres_eos_mf->nGrow()!=1)
  amrex::Error("pres_eos_mf->nGrow()!=1");

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,1,28); 

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
      nucleate_pos[pos_comp]=pos_sites[pos_comp]; 
      pos_comp++;
     }
     double rr=pos_sites[pos_comp];
     if (rr>0.0) {
      if (rr<2.0*dxmax)
       rr=2.0*dxmax;
     } else
      amrex::Error("rr must be positive");
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
      nucleate_pos[pos_comp]=pos_sites[pos_comp]; 
      pos_comp++;
     }
     double rr=pos_sites[pos_comp];
     if (rr>0.0) {
      if (rr<2.0*dxmax)
       rr=2.0*dxmax;
     } else
      amrex::Error("rr must be positive");
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
         if ((rz_flag==1)&&(dir==0)) {
          xnucleate[dir]=0.0;
         } else if (dir==AMREX_SPACEDIM-1) {  // vertical coordinate
          xnucleate[dir]=pos_sites[nc*4+dir];
         } else {       
          if (xnucleate[dir]-rr<problo[dir])
           xnucleate[dir]=problo[dir]+rr;
          if (xnucleate[dir]+rr>probhi[dir])
           xnucleate[dir]=probhi[dir]-rr;
         }
        } // dir
       } // io proc?

       for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
        ParallelDescriptor::ReduceRealMax(xnucleate[dir]);

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
 } else
  amrex::Error("n_sites invalid");

 int nucleate_pos_size=nucleate_pos.size();

 MultiFab& LS_new = get_new_data(LS_Type,slab_step+1);
 if (LS_new.nComp()!=nmat*(1+AMREX_SPACEDIM)) 
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
  if (localMF[BURNING_VELOCITY_MF]->nGrow()!=ngrow_make_distance)
   amrex::Error("localMF[BURNING_VELOCITY_MF] incorrect ngrow");

  if (localMF[SATURATION_TEMP_MF]->nComp()!=ntsat)
   amrex::Error("localMF[SATURATION_TEMP_MF]->nComp()!=ntsat");
  if (localMF[SATURATION_TEMP_MF]->nGrow()!=ngrow_make_distance)
   amrex::Error("localMF[SATURATION_TEMP_MF] incorrect ngrow");

  if (localMF[FD_NRM_ND_MF]->nComp()!=n_normal)
   amrex::Error("localMF[FD_NRM_ND_MF]->nComp()!=n_normal");
  if (localMF[FD_NRM_ND_MF]->nGrow()!=ngrow_make_distance+1)
   amrex::Error("localMF[FD_NRM_ND_MF] incorrect ngrow");

  if (localMF[FD_CURV_CELL_MF]->nComp()!=2*(nmat+nten))
   amrex::Error("localMF[FD_CURV_CELL_MF]->nComp()!=2*(nmat+nten)");
  if (localMF[FD_CURV_CELL_MF]->nGrow()!=ngrow_make_distance)
   amrex::Error("localMF[FD_CURV_CELL_MF] incorrect ngrow");

  debug_ngrow(HOLD_LS_DATA_MF,normal_probe_size+3,30);
  debug_ngrow(HOLD_LS_DATA_MF,ngrow_distance,30);
  if (localMF[HOLD_LS_DATA_MF]->nComp()!=nmat*(1+AMREX_SPACEDIM)) 
   amrex::Error("localMF[HOLD_LS_DATA_MF]->nComp() invalid");

  debug_ngrow(LS_NRM_FD_MF,1,30);
  if (localMF[LS_NRM_FD_MF]->nComp()!=nmat*AMREX_SPACEDIM) 
   amrex::Error("localMF[LS_NRM_FD_MF]->nComp() invalid");

  debug_ngrow(MDOT_MF,0,355);

  VOF_Recon_resize(normal_probe_size+3,SLOPE_RECON_MF);

  int ngrow_dest=ngrow_distance-1;

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
   FArrayBox& nrmFDfab=(*localMF[FD_NRM_ND_MF])[mfi];

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   FORT_FD_NODE_NORMAL( 
    &level,
    &finest_level,
    lsfab.dataPtr(),
    ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
    nrmFDfab.dataPtr(),
    ARLIM(nrmFDfab.loVect()),ARLIM(nrmFDfab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    xlo,dx,
    &nmat,
    &nten,
    &n_normal,
    &ngrow_dest);
  } // mfi
} // omp
  ns_reconcile_d_num(70);

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
   FArrayBox& nrmFDfab=(*localMF[FD_NRM_ND_MF])[mfi];
   FArrayBox& curvfab=(*localMF[FD_CURV_CELL_MF])[mfi];

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // the normals fab has ngrow_dest+1 ghost cells
    //    growntileboxNODE(ngrow_dest)
    // the curvature fab should have ngrow_dest ghost cells.
   FORT_NODE_TO_CELL( 
    &level,
    &finest_level,
    lsfab.dataPtr(),
    ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
    nrmFDfab.dataPtr(),
    ARLIM(nrmFDfab.loVect()),ARLIM(nrmFDfab.hiVect()),
    curvfab.dataPtr(),
    ARLIM(curvfab.loVect()),ARLIM(curvfab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    xlo,dx,
    &nmat,
    &nten,
    &n_normal,
    &ngrow_dest);
  } // mfi
} // omp
  ns_reconcile_d_num(70);

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

   int stefan_flag=1;
   Vector<int> use_exact_temperature(2*nten);
   for (int im=0;im<2*nten;im++)
    use_exact_temperature[im]=0;

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   if (nucleation_flag==0) {

    FArrayBox& lsfab=(*localMF[HOLD_LS_DATA_MF])[mfi];
    FArrayBox& nrmFDfab=(*localMF[LS_NRM_FD_MF])[mfi];

    FArrayBox& burnvelfab=(*localMF[BURNING_VELOCITY_MF])[mfi];
    if (burnvelfab.nComp()!=nburning)
     amrex::Error("burnvelfab.nComp() incorrect");

     // ntsat=nten*(ncomp_per_tsat+1)
     // e.g. for interface 12,
     //  component 1=0 if T_gamma,Y_gamma not defined
     //             =1 if T_gamma,Y_gamma defined in RATEMASSCHANGE
     //             =2 if T_gamma,Y_gamma defined after extrapolation
     //             =-1 or -2 for condensation case.
     //  component 2=T_gamma
     //  component 3=Y_gamma
     //  repeats ....
    FArrayBox& Tsatfab=(*localMF[SATURATION_TEMP_MF])[mfi];
    if (Tsatfab.nComp()!=ntsat)
     amrex::Error("Tsatfab.nComp()!=ntsat");

    FArrayBox& colorfab=(*localMF[COLOR_MF])[mfi];
    FArrayBox& typefab=(*localMF[TYPE_MF])[mfi];
     // fluids tessellate; solids overlap
    FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi]; 

    FArrayBox& curvfab=(*localMF[FD_CURV_CELL_MF])[mfi];

     // if stefan_flag==1:
     // lsnewfab and burnvelfab are updated.
     // lsfab is not updated.
     // burnvelfab=BURNING_VELOCITY_MF is cell centered velocity.
    FORT_RATEMASSCHANGE( 
     &tid_current,
     &nucleation_flag,
     &stefan_flag,
     &level,
     &finest_level,
     &normal_probe_size,
     &ngrow_distance,
     &nstate,
     &nmat,
     &nten,
     &nburning,
     &ntsat,
     &nden,
     density_floor_expansion.dataPtr(),
     density_ceiling_expansion.dataPtr(),
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
     latent_heat.dataPtr(),
     use_exact_temperature.dataPtr(),
     reaction_rate.dataPtr(),
     hardwire_Y_gamma.dataPtr(),
     hardwire_T_gamma.dataPtr(),
     saturation_temp.dataPtr(),
     saturation_temp_curv.dataPtr(),
     saturation_temp_vel.dataPtr(),
     saturation_temp_min.dataPtr(),
     saturation_temp_max.dataPtr(),
     freezing_model.dataPtr(),
     Tanasawa_or_Schrage.dataPtr(),
     distribute_from_target.dataPtr(),
     mass_fraction_id.dataPtr(),
     species_evaporation_density.dataPtr(),
     molar_mass.dataPtr(),
     species_molar_mass.dataPtr(),
     &use_supermesh,
     tilelo,tilehi,
     fablo,fabhi,&bfact,
     xlo,dx,
     &prev_time_slab,
     &dt_slab,
     &blob_arraysize,
     blob_array.dataPtr(),
     &num_elements_blobclass,
     &color_count,
     colorfab.dataPtr(),
     ARLIM(colorfab.loVect()),ARLIM(colorfab.hiVect()),
     typefab.dataPtr(),
     ARLIM(typefab.loVect()),ARLIM(typefab.hiVect()),
     maskcov.dataPtr(),
     ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
     burnvelfab.dataPtr(),
     ARLIM(burnvelfab.loVect()),ARLIM(burnvelfab.hiVect()),
     Tsatfab.dataPtr(),
     ARLIM(Tsatfab.loVect()),ARLIM(Tsatfab.hiVect()),
     lsfab.dataPtr(),ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
     lsnewfab.dataPtr(),ARLIM(lsnewfab.loVect()),ARLIM(lsnewfab.hiVect()),
     snewfab.dataPtr(),ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()),
     nrmFDfab.dataPtr(),ARLIM(nrmFDfab.loVect()),ARLIM(nrmFDfab.hiVect()),
     eosfab.dataPtr(),ARLIM(eosfab.loVect()),ARLIM(eosfab.hiVect()),
     reconfab.dataPtr(),ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
     presfab.dataPtr(),ARLIM(presfab.loVect()),ARLIM(presfab.hiVect()),
     pres_eos_fab.dataPtr(),
     ARLIM(pres_eos_fab.loVect()),ARLIM(pres_eos_fab.hiVect()),
     curvfab.dataPtr(),
     ARLIM(curvfab.loVect()),
     ARLIM(curvfab.hiVect()));

   } else if (nucleation_flag==1) {

    FORT_RATEMASSCHANGE( 
     &tid_current,
     &nucleation_flag,
     &stefan_flag,
     &level,
     &finest_level,
     &normal_probe_size,
     &ngrow_distance,
     &nstate,
     &nmat,
     &nten,
     &nburning,
     &ntsat,
     &nden,
     density_floor_expansion.dataPtr(),
     density_ceiling_expansion.dataPtr(),
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
     latent_heat.dataPtr(),
     use_exact_temperature.dataPtr(),
     reaction_rate.dataPtr(),
     hardwire_Y_gamma.dataPtr(),
     hardwire_T_gamma.dataPtr(),
     saturation_temp.dataPtr(),
     saturation_temp_curv.dataPtr(),
     saturation_temp_vel.dataPtr(),
     saturation_temp_min.dataPtr(),
     saturation_temp_max.dataPtr(),
     freezing_model.dataPtr(),
     Tanasawa_or_Schrage.dataPtr(),
     distribute_from_target.dataPtr(),
     mass_fraction_id.dataPtr(),
     species_evaporation_density.dataPtr(),
     molar_mass.dataPtr(),
     species_molar_mass.dataPtr(),
     &use_supermesh,
     tilelo,tilehi,
     fablo,fabhi,&bfact,
     xlo,dx,
     &prev_time_slab,
     &dt_slab,
     &blob_arraysize,
     blob_array.dataPtr(),
     &num_elements_blobclass,
     &color_count,
     lsnewfab.dataPtr(), //colorfab
     ARLIM(lsnewfab.loVect()),ARLIM(lsnewfab.hiVect()),
     lsnewfab.dataPtr(), //typefab
     ARLIM(lsnewfab.loVect()),ARLIM(lsnewfab.hiVect()),
     maskcov.dataPtr(),
     ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
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
     lsnewfab.dataPtr(), //nrmFDfab
     ARLIM(lsnewfab.loVect()),ARLIM(lsnewfab.hiVect()),
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
 ns_reconcile_d_num(70);

 delete_localMF(DEN_RECON_MF,1);

 delete presmf;
 delete pres_eos_mf;

} // subroutine level_phase_change_rate



void
NavierStokes::level_phase_change_rate_extend() {

 bool use_tiling=ns_tiling;
 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid level_phase_change_rate_extend");

 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

 int ncomp_per_burning=AMREX_SPACEDIM;
 int ncomp_per_tsat=2;

 int nburning=nten*(ncomp_per_burning+1);
 int ntsat=nten*(ncomp_per_tsat+1);

 const Real* dx = geom.CellSize();

 MultiFab& LS_new = get_new_data(LS_Type,slab_step+1);
 if (LS_new.nComp()!=nmat*(AMREX_SPACEDIM+1)) 
  amrex::Error("LS_new invalid ncomp");

 if (localMF[BURNING_VELOCITY_MF]->nComp()!=nburning)
  amrex::Error("localMF[BURNING_VELOCITY_MF] incorrect ncomp");
 if (localMF[SATURATION_TEMP_MF]->nComp()!=ntsat)
  amrex::Error("localMF[SATURATION_TEMP_MF] incorrect ncomp");

 if (ngrow_make_distance!=3)
  amrex::Error("expecting ngrow_make_distance==3");

 if (localMF[BURNING_VELOCITY_MF]->nGrow()!=ngrow_make_distance)
  amrex::Error("localMF[BURNING_VELOCITY_MF] incorrect ngrow");
 if (localMF[SATURATION_TEMP_MF]->nGrow()!=ngrow_make_distance)
  amrex::Error("localMF[SATURATION_TEMP_MF] incorrect ngrow");

 debug_ngrow(HOLD_LS_DATA_MF,normal_probe_size+3,30);
 if (localMF[HOLD_LS_DATA_MF]->nComp()!=nmat*(1+AMREX_SPACEDIM)) 
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

  ncomp=nten+nten*ncomp_per_interface;

  Vector<int> scompBC_map;
  scompBC_map.resize(ncomp);
   // extrap, u_extrap, v_extrap, w_extrap
   // mof recon extrap
   // maskSEMextrap
  int burnvel_start_pos_base=1+AMREX_SPACEDIM+nmat*ngeom_recon+1;
  int extend_start_pos=burnvel_start_pos_base;
  if (velflag==0) {
   extend_start_pos=burnvel_start_pos_base+nburning;
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

  PCINTERP_fill_borders(local_mf,ngrow_make_distance,
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

   int ngrow=normal_probe_size+3;
   if (ngrow!=4)
    amrex::Error("expecting ngrow==4");

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // burnvelfab=BURNING_VELOCITY_MF or 
    // burnvelfab=SATURATION_TEMP_MF is cell centered.
    // sets the burning velocity/saturation temp flag from 0 to 2 if
    // foot of characteristic within range.
    // in: MASS_TRANSFER_3D.F90
   FORT_EXTEND_BURNING_VEL( 
    &velflag,
    &level,
    &finest_level,
    xlo,dx,
    &nmat,
    &nten,
    &ncomp,
    &ngrow,
    latent_heat.dataPtr(),
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    burnvelfab.dataPtr(),
    ARLIM(burnvelfab.loVect()),ARLIM(burnvelfab.hiVect()),
    lsfab.dataPtr(),ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()));
  } // mfi
} // omp
  ns_reconcile_d_num(71);

  scompBC_map.resize(ncomp);

  for (int imdest=0;imdest<ncomp;imdest++)
   scompBC_map[imdest]=extend_start_pos+imdest;

  PCINTERP_fill_borders(local_mf,ngrow_make_distance,
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

} // subroutine level_phase_change_rate_extend

// isweep==0:
// 1. initialize node velocity from BURNING_VELOCITY_MF
// 2. unsplit advection of materials changing phase
// 3. determine overall change in volume
// isweep==1:
// 1. scale overall change in volume so that amount evaporated equals the
//    amount condensed if heat pipe problem.
// 2. update volume fractions, jump strength, temperature
void
NavierStokes::level_phase_change_convert(int isweep) {

 bool use_tiling=ns_tiling;
 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid level_phase_change_convert");

 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

 int ncomp_per_burning=AMREX_SPACEDIM;
 int ncomp_per_tsat=2;

 int nburning=nten*(ncomp_per_burning+1);
 int ntsat=nten*(ncomp_per_tsat+1);

 int nden=nmat*num_state_material;
 int nstate=num_materials_vel*(AMREX_SPACEDIM+1)+
  nmat*(num_state_material+ngeom_raw)+1;
 int scomp_mofvars=num_materials_vel*(AMREX_SPACEDIM+1)+
  nmat*num_state_material;

 // mask=1 if not covered or if outside the domain.
 // NavierStokes::maskfiner_localMF
 // NavierStokes::maskfiner
 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,1,28); 

 if (localMF[JUMP_STRENGTH_MF]->nGrow()!=ngrow_expansion)
  amrex::Error("jump strength invalid ngrow level_phase_change_conv");
 if (localMF[JUMP_STRENGTH_MF]->nComp()!=2*nten)
  amrex::Error("localMF[JUMP_STRENGTH_MF]->nComp() invalid");

 if (localMF[SATURATION_TEMP_MF]->nComp()!=ntsat)
  amrex::Error("localMF[SATURATION_TEMP_MF]->nComp()!=ntsat");
 if (localMF[SATURATION_TEMP_MF]->nGrow()!=ngrow_make_distance)
  amrex::Error("localMF[SATURATION_TEMP_MF] incorrect ngrow");

 if (localMF[nodevel_MF]->nGrow()!=1)
  amrex::Error("localMF[nodevel_MF]->nGrow()  invalid");
 if (localMF[nodevel_MF]->nComp()!=2*nten*AMREX_SPACEDIM)
  amrex::Error("localMF[nodevel_MF]->nComp()  invalid");

 if (localMF[deltaVOF_MF]->nComp()!=3*nmat)
  amrex::Error("localMF[deltaVOF_MF]->nComp()  invalid");

 if (localMF[BURNING_VELOCITY_MF]->nComp()!=nburning)
  amrex::Error("burning vel invalid ncomp");

  // in: level_phase_change_convert
  // DEN_RECON_MF is initialized prior to isweep==0
 if (localMF[DEN_RECON_MF]->nComp()!=nden)
  amrex::Error("DEN_RECON_MF invalid ncomp");

 debug_ngrow(HOLD_LS_DATA_MF,ngrow_distance,30);
 if (localMF[HOLD_LS_DATA_MF]->nComp()!=nmat*(1+AMREX_SPACEDIM)) 
  amrex::Error("localMF[HOLD_LS_DATA_MF]->nComp() invalid");

 MultiFab& LS_new = get_new_data(LS_Type,slab_step+1);
 if (LS_new.nComp()!=nmat*(AMREX_SPACEDIM+1)) 
  amrex::Error("LS_new invalid ncomp");

 MultiFab& S_new = get_new_data(State_Type,slab_step+1);
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");

 const Real* dx = geom.CellSize();

 Vector< Vector<Real> > DVOF_local;
 DVOF_local.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  DVOF_local[tid].resize(nmat);
  for (int im=0;im<nmat;im++)
   DVOF_local[tid][im]=0.0;
 } // tid

 Vector< Vector<Real> > delta_mass_local;
 delta_mass_local.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  delta_mass_local[tid].resize(2*nmat); // source 1..nmat  dest 1..nmat
  for (int im=0;im<2*nmat;im++)
   delta_mass_local[tid][im]=0.0;
 } // tid


 if (isweep==0) {

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
   Vector<int> velbc=getBCArray(State_Type,gridno,0,
    num_materials_vel*AMREX_SPACEDIM);

   FArrayBox& burnvelfab=(*localMF[BURNING_VELOCITY_MF])[mfi];
   FArrayBox& nodevelfab=(*localMF[nodevel_MF])[mfi];
   if (burnvelfab.nComp()==nburning) {
    // do nothing
   } else 
    amrex::Error("burnvelfab.nComp() invalid");

   if (nodevelfab.nComp()==2*nten*AMREX_SPACEDIM) {
    // do nothing
   } else 
    amrex::Error("nodevelfab.nComp() invalid");

   int bfact=parent->Space_blockingFactor(level);

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     // burnvelfab=BURNING_VELOCITY_MF is cell centered (interior: lo to hi)
     // nodevelfab=nodevel is at the nodes. (interior: lo to hi+1)
   FORT_NODEDISPLACE(
    &nmat,
    &nten,
    &nburning,
    tilelo,tilehi,
    fablo,fabhi,
    &bfact, 
    velbc.dataPtr(),
    &dt_slab,
    nodevelfab.dataPtr(),
    ARLIM(nodevelfab.loVect()),ARLIM(nodevelfab.hiVect()),
    burnvelfab.dataPtr(),
    ARLIM(burnvelfab.loVect()),ARLIM(burnvelfab.hiVect()),
    xlo,dx, 
    &level,&finest_level);
  } // mfi
} // omp
  ns_reconcile_d_num(72);

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
   int ncomp=2*nten*AMREX_SPACEDIM;
   int dirplot=-1;
   int id=0;
   std::cout << "dt_slab = " << dt_slab << '\n';
   tecplot_debug(nodevelfab,xlo,fablo,fabhi,dxplot,dirplot,id,
     scomp,ncomp,interior_only);
  }

 } else if (isweep==1) {
  // do nothing
 } else
  amrex::Error("isweep invalid");

 if (isweep==0) {
  // do nothing
 } else if (isweep==1) {
  for (int im=0;im<nmat;im++)
   DVOF_local[0][im]=DVOF[0][im];
 } else
  amrex::Error("isweep invalid");

 VOF_Recon_resize(normal_probe_size+3,SLOPE_RECON_MF);

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

   Vector<int> vofbc=getBCArray(State_Type,gridno,scomp_mofvars,1);

    // mask=tag if not covered by level+1 or outside the domain.
   FArrayBox& maskcov=(*localMF[MASKCOEF_MF])[mfi];

   FArrayBox& deltafab=(*localMF[deltaVOF_MF])[mfi];

   FArrayBox& nodevelfab=(*localMF[nodevel_MF])[mfi];
   if (nodevelfab.nComp()==2*nten*AMREX_SPACEDIM) {
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
   if (Tsatfab.nComp()!=ntsat)
    amrex::Error("Tsatfab.nComp()!=ntsat");

   int bfact=parent->Space_blockingFactor(level);

   int tid_update=0;
   if (isweep==0) {
    tid_update=ns_thread();
   } else if (isweep==1) {
    tid_update=0; // only use 0th component for 2nd sweep
   } else
    amrex::Error("isweep invalid");

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // nodevelfab (=nodevel) is node based (interior: lo..hi+1)
    // in: MASS_TRANSFER_3D.F90
   FORT_CONVERTMATERIAL( 
    &tid_current,
    &isweep, 
    &solvability_projection, // if solvability_projection==1 => net growth=0
    &ngrow_expansion,
    &level,&finest_level,
    &normal_probe_size,
    &nmat,
    &nten,
    &nden,
    &nstate,
    &ntsat,
    density_floor_expansion.dataPtr(),
    density_ceiling_expansion.dataPtr(),
    latent_heat.dataPtr(),
    saturation_temp.dataPtr(),
    freezing_model.dataPtr(),
    mass_fraction_id.dataPtr(),
    species_evaporation_density.dataPtr(),
    distribute_from_target.dataPtr(),
    tilelo,tilehi,
    fablo,fabhi,
    &bfact, 
    vofbc.dataPtr(),xlo,dx,
    &dt_slab,
    delta_mass_local[tid_current].dataPtr(),
    DVOF_local[tid_update].dataPtr(),
    maskcov.dataPtr(),
    ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
    deltafab.dataPtr(),
    ARLIM(deltafab.loVect()),ARLIM(deltafab.hiVect()),
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
 ns_reconcile_d_num(73);

 if (isweep==0) {
  for (int tid=1;tid<thread_class::nthreads;tid++) {
   for (int im=0;im<nmat;im++) 
    DVOF_local[0][im]+=DVOF_local[tid][im];
  } // tid
 } else if (isweep==1) {
  for (int tid=1;tid<thread_class::nthreads;tid++) {
   for (int im=0;im<2*nmat;im++) {
    delta_mass_local[0][im]+=delta_mass_local[tid][im];
   }
  } // tid
 } else {
  amrex::Error("isweep invalid");
 }

 ParallelDescriptor::Barrier();

 if (isweep==0) {
  for (int im=0;im<nmat;im++) {
   ParallelDescriptor::ReduceRealSum(DVOF_local[0][im]);
   DVOF[0][im]+=DVOF_local[0][im];
  }  // im
 } else if (isweep==1) {
  for (int im=0;im<2*nmat;im++) {
   ParallelDescriptor::ReduceRealSum(delta_mass_local[0][im]);
   delta_mass[0][im]+=delta_mass_local[0][im];
  }
 } else
  amrex::Error("isweep invalid");

 localMF[JUMP_STRENGTH_MF]->FillBoundary(geom.periodicity());
 avgDown_localMF(JUMP_STRENGTH_MF,0,2*nten,0);

} // subroutine level_phase_change_convert

void
NavierStokes::phase_change_redistributeALL() {

 if (level!=0)
  amrex::Error("level invalid phase_change_redistributeALL");

 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

 int finest_level=parent->finestLevel();
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.getStateDist_localMF(LSNEW_MF,2*ngrow_expansion,cur_time_slab,4);
 }

 mdotplus.resize(thread_class::nthreads);
 mdotminus.resize(thread_class::nthreads);
 mdotcount.resize(thread_class::nthreads);
 mdot_lost.resize(thread_class::nthreads);
 mdot_sum.resize(thread_class::nthreads);
 mdot_sum2.resize(thread_class::nthreads);

 for (int im=1;im<=nmat;im++) {
  for (int im_opp=im+1;im_opp<=nmat;im_opp++) {
   for (int ireverse=0;ireverse<=1;ireverse++) {
    if ((im>nmat)||(im_opp>nmat))
     amrex::Error("im or im_opp bust 200cpp");
    int iten;
    get_iten_cpp(im,im_opp,iten,nmat);
    if ((iten<1)||(iten>nten))
     amrex::Error("iten invalid");

    int indexEXP=iten+ireverse*nten-1;

    Real LL=latent_heat[indexEXP];
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
     }

     allocate_array(2*ngrow_expansion,1,-1,donorflag_MF);
     setVal_array(2*ngrow_expansion,1,0.0,donorflag_MF);

     for (int isweep_redistribute=0;isweep_redistribute<3;
	  isweep_redistribute++) {

      for (int ilev=finest_level;ilev>=level;ilev--) {
       NavierStokes& ns_level=getLevel(ilev);
       ns_level.level_phase_change_redistribute(
        expect_mdot_sign,im_source,im_dest,indexEXP,
        isweep_redistribute);
      } // ilev=finest_level ... level

      if (isweep_redistribute==0) {
       std::cout << "before:imsrc,imdst,mdot_sum " <<
        im_source << ' ' << im_dest << ' ' << mdot_sum[0] << '\n';
      } else if (isweep_redistribute==1) {
       // do nothing
      } else if (isweep_redistribute==2) {
       std::cout << "after:imsrc,imdst,mdot_sum2 " <<   
        im_source << ' ' << im_dest << ' ' << mdot_sum2[0] << '\n';
       std::cout << "after:imsrc,imdst,mdot_lost " <<   
        im_source << ' ' << im_dest << ' ' << mdot_lost[0] << '\n';
       std::cout << "imsrc,imdst,mdot_sum2+mdot_lost " <<   
        im_source << ' ' << im_dest << ' ' <<   
        mdot_sum2[0]+mdot_lost[0] << '\n';
      } else
       amrex::Error("isweep_redistribute invalid");

     } // isweep_redistribute=0,1,2

     delete_array(donorflag_MF);

    } // LL!=0
   } // ireverse
  } // im_opp
 } // im=1..nmat


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
   indexEXP_filler,isweep_combine);
 } // ilev=finest_level ... level

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "mdotplus = " << mdotplus[0] << '\n';
   std::cout << "mdotminus = " << mdotminus[0] << '\n';
   std::cout << "mdotcount = " << mdotcount[0] << '\n';
  } // IOProc?
 } // verbose>0

 if (solvability_projection==1) { // require net growth = 0

  int isweep_solvability=4;

  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   ns_level.level_phase_change_redistribute(
    expect_mdot_sign_filler,
    im_source_filler,im_dest_filler,
    indexEXP_filler,isweep_solvability);
  } // ilev=finest_level ... level

 } else if (solvability_projection==0) {
  // do nothing
 } else
  amrex::Error("solvability_projection invalid");

 delete_array(LSNEW_MF);
 delete_array(HOLD_LS_DATA_MF);

} // subroutine phase_change_redistributeALL

void
NavierStokes::level_phase_change_redistribute(
 Real expect_mdot_sign,
 int im_source,int im_dest,int indexEXP,
 int isweep) {

 bool use_tiling=ns_tiling;
 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid level_phase_change_redistribute");

 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;
 int scomp_mofvars=num_materials_vel*(AMREX_SPACEDIM+1)+
  nmat*num_state_material;

 debug_ngrow(JUMP_STRENGTH_MF,ngrow_expansion,355);
 if (localMF[JUMP_STRENGTH_MF]->nComp()!=2*nten)
  amrex::Error("localMF[JUMP_STRENGTH_MF]->nComp()!=2*nten level_phase ...");

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,1,6001);
 
 if (localMF[LSNEW_MF]->nComp()!=nmat*(1+AMREX_SPACEDIM))
  amrex::Error("localMF[LSNEW_MF]->nComp() invalid");
 debug_ngrow(LSNEW_MF,2*ngrow_expansion,6001);

 const Real* dx = geom.CellSize();

  // tags for redistribution of source term
  // 1=> donor  2=> receiver  0=> neither

 Real LL=0.0;

 if ((isweep==0)||(isweep==1)||(isweep==2)) {
  if (localMF[donorflag_MF]->nGrow()!=2*ngrow_expansion)
   amrex::Error("localMF[donorflag_MF]->ngrow() invalid");
  if (localMF[donorflag_MF]->nComp()!=1)
   amrex::Error("localMF[donorflag_MF]->nComp() invalid");

  if ((indexEXP>=0)&&(indexEXP<2*nten)) {
   LL=latent_heat[indexEXP];
  } else
   amrex::Error("indexEXP invalid");

 } else if ((isweep==3)||(isweep==4)) {

  if (indexEXP==-1) {
   LL=0.0;
  } else
   amrex::Error("indexEXP invalid");

 } else
  amrex::Error("isweep invalid");

 VOF_Recon_resize(1,SLOPE_RECON_MF);
  
 if (isweep==0) {
 
  if (LL==0.0)
   amrex::Error("LL invalid");
  if (std::abs(expect_mdot_sign)!=1.0)
   amrex::Error("expect_mdot_sign invalid");
  if ((im_source<1)||(im_source>nmat))
   amrex::Error("im_source invalid");
  if ((im_dest<1)||(im_dest>nmat))
   amrex::Error("im_dest invalid");
  if ((indexEXP<0)||(indexEXP>=2*nten))
   amrex::Error("indexEXP invalid");
  
  Vector< Real > mdot_sum_local;
  mdot_sum_local.resize(thread_class::nthreads);
  for (int tid=0;tid<thread_class::nthreads;tid++) {
   mdot_sum_local[tid]=0.0;
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
   Vector<int> vofbc=getBCArray(State_Type,gridno,scomp_mofvars,1);

   FArrayBox& maskcov=(*localMF[MASKCOEF_MF])[mfi];

   FArrayBox& donorfab=(*localMF[donorflag_MF])[mfi];
   FArrayBox& JUMPfab=(*localMF[JUMP_STRENGTH_MF])[mfi];
   FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi]; 
   FArrayBox& newdistfab=(*localMF[LSNEW_MF])[mfi];
   int bfact=parent->Space_blockingFactor(level);
   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // isweep==0 
   FORT_TAGEXPANSION( 
    latent_heat.dataPtr(),
    freezing_model.dataPtr(),
    distribute_from_target.dataPtr(),
    &ngrow_expansion,
    &cur_time_slab,
    vofbc.dataPtr(),
    &expect_mdot_sign,
    &mdot_sum_local[tid_current],
    &im_source,
    &im_dest,
    &indexEXP,
    &level,&finest_level,
    &nmat,&nten, 
    tilelo,tilehi,
    fablo,fabhi,
    &bfact, 
    xlo,dx,&dt_slab,
    maskcov.dataPtr(),
    ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
    donorfab.dataPtr(),
    ARLIM(donorfab.loVect()),ARLIM(donorfab.hiVect()),
    JUMPfab.dataPtr(),
    ARLIM(JUMPfab.loVect()),ARLIM(JUMPfab.hiVect()),
    newdistfab.dataPtr(),
    ARLIM(newdistfab.loVect()),ARLIM(newdistfab.hiVect()),
    reconfab.dataPtr(),
    ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()));
 
  } // mfi
} // omp
  ns_reconcile_d_num(74);

  for (int tid=1;tid<thread_class::nthreads;tid++) {
   mdot_sum_local[0]+=mdot_sum_local[tid];
  }
  ParallelDescriptor::ReduceRealSum(mdot_sum_local[0]);
  mdot_sum[0]+=mdot_sum_local[0];

  localMF[donorflag_MF]->FillBoundary(geom.periodicity());
  avgDown_tag_localMF(donorflag_MF);

 } else if (isweep==1) {

   // redistribution.

  if (LL==0.0)
   amrex::Error("LL invalid");
  if (std::abs(expect_mdot_sign)!=1.0)
   amrex::Error("expect_mdot_sign invalid");
  if ((im_source<1)||(im_source>nmat))
   amrex::Error("im_source invalid");
  if ((im_dest<1)||(im_dest>nmat))
   amrex::Error("im_dest invalid");
  if ((indexEXP<0)||(indexEXP>=2*nten))
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
    Vector<int> vofbc=getBCArray(State_Type,gridno,scomp_mofvars,1);

    FArrayBox& maskcov=(*localMF[MASKCOEF_MF])[mfi];
    FArrayBox& donorfab=(*localMF[donorflag_MF])[mfi];
    FArrayBox& JUMPfab=(*localMF[JUMP_STRENGTH_MF])[mfi];
    FArrayBox& newdistfab=(*localMF[LSNEW_MF])[mfi];

    int bfact=parent->Space_blockingFactor(level);

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     // isweep==1
    FORT_DISTRIBUTEEXPANSION( 
     &ngrow_expansion,
     &im_source,
     &im_dest,
     &indexEXP,
     &level,&finest_level,
     &nmat,&nten, 
     tilelo,tilehi,
     fablo,fabhi,
     &bfact, 
     xlo,dx,&dt_slab,
     maskcov.dataPtr(),
     ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
     newdistfab.dataPtr(),
     ARLIM(newdistfab.loVect()),ARLIM(newdistfab.hiVect()),
     donorfab.dataPtr(),
     ARLIM(donorfab.loVect()),ARLIM(donorfab.hiVect()),
     JUMPfab.dataPtr(),
     ARLIM(JUMPfab.loVect()),ARLIM(JUMPfab.hiVect()));
  } // mfi
} //omp
  ns_reconcile_d_num(75);

 } else if (isweep==2) {

   // clear out mdot in the donor cells.

  if (LL==0.0)
   amrex::Error("LL invalid");
  if (std::abs(expect_mdot_sign)!=1.0)
   amrex::Error("expect_mdot_sign invalid");
  if ((im_source<1)||(im_source>nmat))
   amrex::Error("im_source invalid");
  if ((im_dest<1)||(im_dest>nmat))
   amrex::Error("im_dest invalid");
  if ((indexEXP<0)||(indexEXP>=2*nten))
   amrex::Error("indexEXP invalid");

  Vector< Real > mdot_lost_local;
  Vector< Real > mdot_sum2_local;
  mdot_lost_local.resize(thread_class::nthreads);
  mdot_sum2_local.resize(thread_class::nthreads);
  for (int tid=0;tid<thread_class::nthreads;tid++) {
   mdot_sum2_local[tid]=0.0;
   mdot_lost_local[tid]=0.0;
  }
  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(localMF[donorflag_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(*localMF[donorflag_MF],use_tiling);mfi.isValid(); ++mfi) {
    BL_ASSERT(grids[mfi.index()] == mfi.validbox());
    const int gridno = mfi.index();
    const Box& tilegrid = mfi.tilebox();
    const Box& fabgrid = grids[gridno];
    const int* tilelo=tilegrid.loVect();
    const int* tilehi=tilegrid.hiVect();
    const int* fablo=fabgrid.loVect();
    const int* fabhi=fabgrid.hiVect();
    const Real* xlo = grid_loc[gridno].lo();
    Vector<int> vofbc=getBCArray(State_Type,gridno,scomp_mofvars,1);

    FArrayBox& maskcov=(*localMF[MASKCOEF_MF])[mfi];
    FArrayBox& donorfab=(*localMF[donorflag_MF])[mfi];
    FArrayBox& JUMPfab=(*localMF[JUMP_STRENGTH_MF])[mfi];

    int bfact=parent->Space_blockingFactor(level);
    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     // isweep==2
    FORT_CLEAREXPANSION( 
     &ngrow_expansion,
     &mdot_sum2_local[tid_current],
     &mdot_lost_local[tid_current],
     &im_source,
     &im_dest,
     &indexEXP,
     &level,&finest_level,
     &nmat,&nten, 
     tilelo,tilehi,
     fablo,fabhi,
     &bfact, 
     xlo,dx,&dt_slab,
     maskcov.dataPtr(),
     ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
     donorfab.dataPtr(),
     ARLIM(donorfab.loVect()),ARLIM(donorfab.hiVect()),
     JUMPfab.dataPtr(),
     ARLIM(JUMPfab.loVect()),ARLIM(JUMPfab.hiVect()));
  } // mfi
} // omp
  ns_reconcile_d_num(76);

  for (int tid=1;tid<thread_class::nthreads;tid++) {
   mdot_sum2_local[0]+=mdot_sum2_local[tid];
   mdot_lost_local[0]+=mdot_lost_local[tid];
  } // tid
  ParallelDescriptor::ReduceRealSum(mdot_sum2_local[0]);
  ParallelDescriptor::ReduceRealSum(mdot_lost_local[0]);
  mdot_sum2[0]+=mdot_sum2_local[0];
  mdot_lost[0]+=mdot_lost_local[0];

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

    int bfact=parent->Space_blockingFactor(level);

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // NavierStokes::allocate_mdot() called at the beginning of
    //  NavierStokes::do_the_advance
    // mdot initialized in NavierStokes::prelim_alloc()
    FORT_INITJUMPTERM( 
     &mdotplus_local[tid_current],
     &mdotminus_local[tid_current],
     &mdotcount_local[tid_current],
     &ngrow_expansion,
     &cur_time_slab,
     &level,&finest_level,
     &nmat,&nten,
     latent_heat.dataPtr(),
     saturation_temp.dataPtr(),
     freezing_model.dataPtr(),
     distribute_from_target.dataPtr(),
     tilelo,tilehi,
     fablo,fabhi,
     &bfact, 
     xlo,dx,&dt_slab,
     maskcov.dataPtr(),
     ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
     JUMPfab.dataPtr(),ARLIM(JUMPfab.loVect()),ARLIM(JUMPfab.hiVect()),
      // mdotfab is incremented.
     mdotfab.dataPtr(),ARLIM(mdotfab.loVect()),ARLIM(mdotfab.hiVect()),
     newdistfab.dataPtr(),
     ARLIM(newdistfab.loVect()),ARLIM(newdistfab.hiVect()));

  } // mfi
} // omp
  ns_reconcile_d_num(77);

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

 } else if (isweep==4) {

  if (solvability_projection!=1)
   amrex::Error("solvability_projection!=1");

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

   int bfact=parent->Space_blockingFactor(level);

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // NavierStokes::allocate_mdot() called at the beginning of
    //  NavierStokes::do_the_advance
    // mdot initialized in NavierStokes::prelim_alloc()
    // mdot updated in nucleate_bubbles.
   FORT_RENORM_MDOT( 
    &mdotplus[0],&mdotminus[0],&mdotcount[0],
    &level,&finest_level,
    &nmat,&nten,
    tilelo,tilehi,
    fablo,fabhi,
    &bfact, 
    xlo,dx,&dt_slab,
    maskcov.dataPtr(),
    ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
    mdotfab.dataPtr(),
    ARLIM(mdotfab.loVect()),ARLIM(mdotfab.hiVect()));

  } // mfi
} // omp
  ns_reconcile_d_num(78);

 } else
  amrex::Error("isweep invalid");

} // subroutine level_phase_change_redistribute

// called from: NavierStokes::make_physics_varsALL
void
NavierStokes::level_init_icemask() {

 bool use_tiling=ns_tiling;
 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid level_init_icemask");

 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

 resize_maskfiner(1,MASKCOEF_MF);
 VOF_Recon_resize(1,SLOPE_RECON_MF);

 debug_ngrow(SLOPE_RECON_MF,1,3);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(FACE_VAR_MF+dir,0,2);

 debug_ngrow(MASKCOEF_MF,1,6001);

 getStateDist_localMF(LSNEW_MF,1,cur_time_slab,5);

 debug_ngrow(LSNEW_MF,1,6001);
 if (localMF[LSNEW_MF]->nComp()!=nmat*(1+AMREX_SPACEDIM))
  amrex::Error("localMF[LSNEW_MF]->nComp() invalid");

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
   int bfact=parent->Space_blockingFactor(level);

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();
   
    // in: GODUNOV_3D.F90
   FORT_INIT_ICEMASK( 
    &cur_time_slab,
    &facecut_index,
    &icefacecut_index,
    &icemask_index,
    &massface_index,
    &vofface_index,
    &ncphys,
    &level,&finest_level,
    &nmat,&nten,
    latent_heat.dataPtr(),
    saturation_temp.dataPtr(),
    freezing_model.dataPtr(),
    distribute_from_target.dataPtr(),
    tilelo,tilehi,
    fablo,fabhi,
    &bfact, 
    xlo,dx,&dt_slab,
    maskcov.dataPtr(),
    ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
    xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()),
    yface.dataPtr(),ARLIM(yface.loVect()),ARLIM(yface.hiVect()),
    zface.dataPtr(),ARLIM(zface.loVect()),ARLIM(zface.hiVect()),
    newdistfab.dataPtr(),
    ARLIM(newdistfab.loVect()),ARLIM(newdistfab.hiVect()),
    reconfab.dataPtr(),ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()));

  } // mfi
} // omp
  ns_reconcile_d_num(79);

  delete_localMF(LSNEW_MF,1);

} // subroutine level_init_icemask


// 1. called if freezing_model==0,5,6.
// 2. multiphase_project->allocate_project_variables->stefan_solver_init
//    (adjust_temperature==1)
//    coeffMF==localMF[OUTER_ITER_PRESSURE_MF]
// 3. multiphase_project->allocate_maccoef->stefan_solver_init
//    update_SEM_forcesALL->allocate_maccoef->stefan_solver_init
//    (adjust_temperature==0)
//    coeffMF==localMF[ALPHANOVOLUME_MF]
// 4. multiphase_project->allocate_FACE_WEIGHT->stefan_solver_init
//    update_SEM_forcesALL->allocate_FACE_WEIGHT->stefan_solver_init
//    diffusion_heatingALL->allocate_FACE_WEIGHT->stefan_solver_init
//    (adjust_temperature==-1)
//    coeffMF==localMF[CELL_DEN_MF]
// if adjust_temperature==1,
//  Snew=(c1 Tn + c2 TSAT)/(c1+c2)
//  coeffMF=(c1 Tn + c2 TSAT)/(c1+c2)
//  c1=rho cv/(dt*sweptfactor)  c2=(1/vol) sum_face Aface k_m/(theta dx)
// else if adjust_temperature==0,
//  coeffMF=c1+c2
// else if adjust_temperature==-1,
//  faceheat_index component of FACE_VAR_MF
//
// for mass fraction:
// (rho Y)_t + div (rho u Y) = div rho D grad Y
// since rho_t + div (rho u)=0,
// rho (Y_t + u dot grad Y)=div rho D grad Y
// "rho D" coefficient is in FACE_VAR_MF,
//   facespecies_index ... facespecies_index+num_species_var-1
// "1/rho" coefficient is in CELL_DEN_MF 

void
NavierStokes::stefan_solver_init(MultiFab* coeffMF,
		int adjust_temperature,
		int project_option) {
 
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

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;
 int nstate=num_materials_vel*(AMREX_SPACEDIM+1)+
  nmat*(num_state_material+ngeom_raw)+1;

 int nsolve=1;
 int nsolveMM=nsolve*num_materials_scalar_solve;

 int ncomp_per_tsat=2;
 int ntsat=nten*(ncomp_per_tsat+1);

 if (localMF[SATURATION_TEMP_MF]->nComp()!=ntsat)
  amrex::Error("localMF[SATURATION_TEMP_MF]->nComp()!=ntsat");
 if (localMF[SATURATION_TEMP_MF]->nGrow()!=ngrow_make_distance)
  amrex::Error("localMF[SATURATION_TEMP_MF] incorrect ngrow");

 resize_metrics(1);
 VOF_Recon_resize(1,SLOPE_RECON_MF);

 debug_ngrow(VOLUME_MF,1,34);
 debug_ngrow(SWEPT_CROSSING_MF,0,34);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(FACE_VAR_MF+dir,0,2);

 debug_ngrow(SLOPE_RECON_MF,1,31);
 if (localMF[SLOPE_RECON_MF]->nComp()!=nmat*ngeom_recon)
  amrex::Error("localMF[SLOPE_RECON_MF]->nComp() invalid");

 debug_ngrow(CELL_DEN_MF,1,28); 
 debug_ngrow(CELL_DEDT_MF,1,28); 

 if (localMF[CELL_DEN_MF]->nComp()!=nmat+1)
  amrex::Error("localMF[CELL_DEN_MF]->nComp() invalid");
 if (localMF[CELL_DEDT_MF]->nComp()!=nmat+1)
  amrex::Error("localMF[CELL_DEDT_MF]->nComp() invalid");


 if (dt_slab<=0.0)
  amrex::Error("dt_slab must be positive");
 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int GFM_flag=0;
 int face_comp_index=0;

 if (is_phasechange==1) {
  if (project_option==2) { // thermal conduction
   face_comp_index=faceheat_index;
   for (int im=0;im<2*nten;im++) {
    if (latent_heat[im]!=0.0) {
     if ((freezing_model[im]==0)|| //fully saturated
         (freezing_model[im]==5)|| //cavitation
         (freezing_model[im]==6))  //Palmore and Desjardins
      GFM_flag=1;
    } else if (latent_heat[im]==0.0) {
     // do nothing
    } else
     amrex::Error("latent_heat[im] invalid");
   } // im=0..2 nten-1
  } else if ((project_option>=100)&&
	     (project_option<100+num_species_var)) {
   face_comp_index=facespecies_index+project_option-100;
   for (int im=0;im<2*nten;im++) {
    if (latent_heat[im]!=0.0) {
     if (freezing_model[im]==6) {  // Palmore and Desjardins
      int ispec=mass_fraction_id[im];
      if ((ispec>=1)&&(ispec<=num_species_var)) {
       if (ispec==project_option-100+1)
        GFM_flag=1;
      } else
       amrex::Error("ispec invalid");
     }
    } else if (latent_heat[im]==0.0) {
     // do nothing
    } else
     amrex::Error("latent_heat[im] invalid");
   } // im=0.. 2 nten -1
  } else
   amrex::Error("project_option invalid");

 } else if (is_phasechange==0) {
  amrex::Error("is_phasechange invalid in stefan_solver_init (1)");
 } else
  amrex::Error("is_phasechange invalid in stefan_solver_init (2)");

 if (GFM_flag!=1)
  amrex::Error("Gibou et al algorithm only used for GFM");

 const Real* dx = geom.CellSize();

 MultiFab* LSmf=getStateDist(1,cur_time_slab,7);  
 if (LSmf->nComp()!=nmat*(1+AMREX_SPACEDIM))
  amrex::Error("LSmf invalid ncomp");
 if (LSmf->nGrow()!=1)
  amrex::Error("LSmf->nGrow()!=1");

  // temperature and density for all of the materials.
 int nden=nmat*num_state_material;
 MultiFab* state_var_mf=getStateDen(1,cur_time_slab);
 if (state_var_mf->nComp()!=nden)
  amrex::Error("state_var_mf->nComp()!=nden");

 MultiFab& S_new = get_new_data(State_Type,slab_step+1);
 if (S_new.nComp()!=nstate)
  amrex::Error("S_new invalid ncomp");

 int mm_areafrac_index=FACE_VAR_MF;
 int mm_cell_areafrac_index=SLOPE_RECON_MF;
 if (num_materials_scalar_solve==nmat) {
  mm_areafrac_index=FACEFRAC_SOLVE_MM_MF;
  mm_cell_areafrac_index=CELLFRAC_MM_MF;
 } else if (num_materials_scalar_solve==1) {
  // do nothing
 } else
  amrex::Error("num_materials_scalar_solve invalid");

 // (ml,mr,2) frac_pair(ml,mr), dist_pair(ml,mr)  
 int nfacefrac=nmat*nmat*2;
 // im_inside,im_outside,3+sdim -->
 //   area, dist_to_line, dist, line normal.
 int ncellfrac=nmat*nmat*(3+AMREX_SPACEDIM);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(mm_areafrac_index+dir,0,111);
 debug_ngrow(mm_cell_areafrac_index,0,113);

 int num_materials_combine=nmat;
 int state_index;
 Vector<int> scomp;
 Vector<int> ncomp;
 int ncomp_check;

 int project_option_thermal=2;
 get_mm_scomp_solver(
   num_materials_combine,
   project_option_thermal,
   state_index,
   scomp,ncomp,ncomp_check);

 if ((ncomp_check!=nmat)||(state_index!=State_Type))
  amrex::Error("(ncomp_check!=nmat)||(state_index!=State_Type)");

 MultiFab* T_list_mf=getState_list(1,scomp,ncomp,cur_time_slab);

 get_mm_scomp_solver(
   num_materials_combine,
   project_option,
   state_index,
   scomp,ncomp,ncomp_check);

 if ((ncomp_check!=nmat)||(state_index!=State_Type))
  amrex::Error("(ncomp_check!=nmat)||(state_index!=State_Type)");

 MultiFab* TorY_list_mf=getState_list(1,scomp,ncomp,cur_time_slab);

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

   FArrayBox& statefab=(*state_var_mf)[mfi];

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
    if (coefffab.nComp()!=nsolveMM) {
     std::cout << "coefffab.nComp()= " << coefffab.nComp() << '\n';
     std::cout << "nsolveMM= " << nsolveMM << '\n';
     std::cout << "adjust_temperature= " << adjust_temperature << '\n';
     amrex::Error("coefffab.nComp() invalid");
    }
   } else if (adjust_temperature==-1) {
    if (coefffab.nComp()<1) {
     std::cout << "coefffab.nComp()= " << coefffab.nComp() << '\n';
     std::cout << "nsolveMM= " << nsolveMM << '\n';
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

   FArrayBox& xfacemm=(*localMF[mm_areafrac_index])[mfi];
   FArrayBox& yfacemm=(*localMF[mm_areafrac_index+1])[mfi];
   FArrayBox& zfacemm=(*localMF[mm_areafrac_index+AMREX_SPACEDIM-1])[mfi];
   FArrayBox& cellfracmm=(*localMF[mm_cell_areafrac_index])[mfi];

   FArrayBox& Tsatfab=(*localMF[SATURATION_TEMP_MF])[mfi];
   if (Tsatfab.nComp()!=ntsat)
    amrex::Error("Tsatfab.nComp()!=ntsat");

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // in: GODUNOV_3D.F90
   FORT_STEFANSOLVER( 
    &project_option, //2=thermal diffusion  or 100 ... 100+num_species_var-1
    &solidheat_flag, //0=diffuse in solid 1=dirichlet 2=Neumann
    microlayer_size.dataPtr(), 
    microlayer_substrate.dataPtr(), 
    microlayer_temperature_substrate.dataPtr(), 
    &adjust_temperature,
    &nmat,
    &nten,
    &nstate,
    &ntsat, // nten*(ncomp_per_tsat+1)
    &nden,  // nmat*num_state_material
    latent_heat.dataPtr(),
    freezing_model.dataPtr(),
    distribute_from_target.dataPtr(),
    saturation_temp.dataPtr(),
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    &level,
    &finest_level,
    &nfacefrac,
    &ncellfrac,
    xlo,dx,
    &dt_slab,
    statefab.dataPtr(),
    ARLIM(statefab.loVect()),ARLIM(statefab.hiVect()),
    Tsatfab.dataPtr(),
    ARLIM(Tsatfab.loVect()),ARLIM(Tsatfab.hiVect()),
    cellfracmm.dataPtr(),
    ARLIM(cellfracmm.loVect()),ARLIM(cellfracmm.hiVect()),
    xfacemm.dataPtr(),ARLIM(xfacemm.loVect()),ARLIM(xfacemm.hiVect()),
    yfacemm.dataPtr(),ARLIM(yfacemm.loVect()),ARLIM(yfacemm.hiVect()),
    zfacemm.dataPtr(),ARLIM(zfacemm.loVect()),ARLIM(zfacemm.hiVect()),
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
 ns_reconcile_d_num(81);

 delete state_var_mf;
 delete LSmf;
 delete T_list_mf;
 delete TorY_list_mf;
 
}  // stefan_solver_init

// MEHDI VAHAB HEAT SOURCE
// T^new=T^* += dt A Q/(rho cv V) 
// called from NavierStokes::allocate_project_variables when project_option==2
void
NavierStokes::heat_source_term_flux_source() {
 
 int finest_level=parent->finestLevel();

 bool use_tiling=ns_tiling;

 resize_metrics(1);

 debug_ngrow(VOLUME_MF,1,34);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(FACE_VAR_MF+dir,0,2);

 if (dt_slab<=0.0)
  amrex::Error("dt_slab must be positive");
 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;
 int nstate=num_materials_vel*(AMREX_SPACEDIM+1)+
  nmat*(num_state_material+ngeom_raw)+1;

 const Real* dx = geom.CellSize();

 MultiFab* LSmf=getStateDist(1,cur_time_slab,8);  
 if (LSmf->nComp()!=nmat*(1+AMREX_SPACEDIM))
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
   if (DeDTfab.nComp()!=nmat+1)
    amrex::Error("DeDTfab.nComp() invalid");

   FArrayBox& denfab=(*localMF[CELL_DEN_MF])[mfi];  // 1/rho
   if (denfab.nComp()!=nmat+1)
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

    // in: GODUNOV_3D.F90
   FORT_HEATSOURCE_FACE( 
    &nmat,&nten,&nstate,
    latent_heat.dataPtr(),
    saturation_temp.dataPtr(),
    tilelo,tilehi,
    fablo,fabhi,
    &bfact,
    xlo,dx,
    &dt_slab,
    &prev_time_slab,
    &level,
    &finest_level,
    lsfab.dataPtr(),ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
    snewfab.dataPtr(),ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()),
    DeDTfab.dataPtr(),ARLIM(DeDTfab.loVect()),ARLIM(DeDTfab.hiVect()),
    denfab.dataPtr(),
    ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
    volfab.dataPtr(),ARLIM(volfab.loVect()),ARLIM(volfab.hiVect()),
    heatx.dataPtr(faceheat_index),ARLIM(heatx.loVect()),ARLIM(heatx.hiVect()),
    heaty.dataPtr(faceheat_index),ARLIM(heaty.loVect()),ARLIM(heaty.hiVect()),
    heatz.dataPtr(faceheat_index),ARLIM(heatz.loVect()),ARLIM(heatz.hiVect()),
    areax.dataPtr(),ARLIM(areax.loVect()),ARLIM(areax.hiVect()),
    areay.dataPtr(),ARLIM(areay.loVect()),ARLIM(areay.hiVect()),
    areaz.dataPtr(),ARLIM(areaz.loVect()),ARLIM(areaz.hiVect()) );
 } // mfi
} // omp
 ns_reconcile_d_num(82);

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
  int dir,int id,
  Real warning_cutoff) {

 if (((verbose==0)||(verbose==1))&&(force_check==0)) {
  // do nothing
 } else if ((verbose==2)||(force_check==1)) {
  std::fflush(NULL);
  int finest_level=parent->finestLevel();
  const Real* dx = geom.CellSize();
  int ndefined=mf->nComp();
  if (ndefined<scomp+ncomp) 
   amrex::Error("scomp,ncomp invalid in aggressive debug");
  if (mf->nGrow()<ngrow)
   amrex::Error("ngrow invalid in aggressive debug");

  if (verbose==2) {
   std::cout << "AGGRESSIVE DEBUG scomp= " << scomp << '\n';
   std::cout << "AGGRESSIVE DEBUG ncomp= " << ncomp << '\n';
   std::cout << "AGGRESSIVE DEBUG ngrow= " << ngrow << '\n';
   std::cout << "AGGRESSIVE DEBUG dir= " << dir << '\n';
   std::cout << "AGGRESSIVE DEBUG id= " << id << '\n';
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
   const Box& growntilegrid = mfi.growntileboxTENSOR(datatype,ngrow,dir);
   const Box& fabgrid = mfi.validbox();

   if (fabgrid!=mfBA[gridno])
    amrex::Error("fabgrid!=mfBA[gridno]");

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

   FORT_AGGRESSIVE(
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
    &ngrow,&dir,&id,
    &verbose,
    &force_check,
    &gridno,&ngrid,&level,&finest_level,
    mffab.dataPtr(),ARLIM(mffab.loVect()),ARLIM(mffab.hiVect()));
  } // mfi
} // omp
  ns_reconcile_d_num(83);

  std::fflush(NULL);
 } else
  amrex::Error("verbose or force_check invalid");

} // subroutine aggressive_debug

void
NavierStokes::synchronize_flux_register(int operation_flag,
 int spectral_loop) {

 if ((operation_flag<0)||(operation_flag>11))
  amrex::Error("operation_flag invalid2");
 
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
  //   (i) unew^{f} in incompressible non-solid regions
  //   (ii) u^{f,save} + (unew^{c}-u^{c,save})^{c->f} in spectral regions or
  //        compressible regions.
  //   (iii) usolid in solid regions
 if (operation_flag==11) {
  ncfluxreg=AMREX_SPACEDIM;
 } else if (operation_flag==10) { // ucell,umac -> umac
  ncfluxreg=AMREX_SPACEDIM;
 } else if (operation_flag==9) {  // density CELL -> MAC
  ncfluxreg=AMREX_SPACEDIM;
 } else if (operation_flag==7) {  // advection
  ncfluxreg=AMREX_SPACEDIM*nfluxSEM;
 } else if (operation_flag==1) { // interp press from cell to MAC.
  ncfluxreg=AMREX_SPACEDIM;
 } else if (operation_flag==2) { // ppot cell-> mac  grad ppot
  ncfluxreg=2*AMREX_SPACEDIM; // (grad ppot)_mac, ppot_mac
 } else if (operation_flag==3) { // ucell -> umac
  ncfluxreg=AMREX_SPACEDIM;
 } else if (operation_flag==4) { // umac -> umac
  ncfluxreg=AMREX_SPACEDIM;
 } else if (operation_flag==5) { // umac -> umac+beta F^cell->mac
  ncfluxreg=AMREX_SPACEDIM;
 } else if (operation_flag==6) { // grad U
  ncfluxreg=AMREX_SPACEDIM*AMREX_SPACEDIM;
 } else if (operation_flag==0) { // grad p
  ncfluxreg=AMREX_SPACEDIM;
 } else if (operation_flag==8) { // grad U (for crossterm)
  ncfluxreg=AMREX_SPACEDIM*AMREX_SPACEDIM;
 } else
  amrex::Error("operation_flag invalid3");

 new_localMF(SEM_FLUXREG_MF,ncfluxreg,1,-1);
 setVal_localMF(SEM_FLUXREG_MF,0.0,0,ncfluxreg,1); 

} // subroutine allocate_flux_register

int 
NavierStokes::end_spectral_loop() {

 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid end_spectral_loop");

 int local_end=1;
 int bfact=parent->Space_blockingFactor(level);

 if ((enable_spectral==0)||(enable_spectral==3)) {
  local_end=1;
 } else if ((enable_spectral==1)||(enable_spectral==2)) {
  if (bfact==1)
   local_end=1;
  else if (bfact>=2)
   local_end=2;
  else
   amrex::Error("bfact invalid");
 } else
  amrex::Error("enable_spectral invalid end_spectral_loop() ");

 return local_end;
 
} // end_spectral_loop

// first: NavierStokes::nonlinear_advection() is called
// second: NavierStokes::SEM_advectALL is called
// This routine is called from  NavierStokes::SEM_advectALL
//
// for advect_iter=0 ... 1 (source_term==0)
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

 int finest_level=parent->finestLevel();
 
 bool use_tiling=ns_tiling;

 int num_colors=0;
 Vector<Real> blob_array;
 blob_array.resize(1);
 int blob_array_size=blob_array.size();

 if (nfluxSEM!=AMREX_SPACEDIM+num_state_base)
  amrex::Error("nfluxSEM!=AMREX_SPACEDIM+num_state_base");

 if ((SDC_outer_sweeps>=0)&&(SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid");

 if (source_term==1) {

  if (advect_iter!=0)
   amrex::Error("advect_iter invalid");

 } else if (source_term==0) {

  // do nothing
  
 } else
  amrex::Error("source_term invalid");

 int nmat=num_materials;

 int nsolveMM_FACE=num_materials_vel;

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 if ((ns_time_order==1)&&
     (enable_spectral==3))
  amrex::Error("(ns_time_order==1)&&(enable_spectral==3)");
 if ((ns_time_order>=2)&&
     ((enable_spectral==0)||
      (enable_spectral==2)))
  amrex::Error("(ns_time_order>=2)&&(enable_spectral==0 or 2)");

 if ((enable_spectral==1)||
     (enable_spectral==2)||
     (ns_time_order>=2)) {

  int nparts=im_solid_map.size();
  if ((nparts<0)||(nparts>nmat))
   amrex::Error("nparts invalid");
  Vector<int> im_solid_map_null;
  im_solid_map_null.resize(1);

  int* im_solid_map_ptr;
  int nparts_def=nparts;
  if (nparts==0) {
   im_solid_map_ptr=im_solid_map_null.dataPtr();
   nparts_def=1;
  } else if ((nparts>=1)&&(nparts<=nmat)) {
   im_solid_map_ptr=im_solid_map.dataPtr();
  } else
   amrex::Error("nparts invalid");

  for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) { 
   if (localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0)
    amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0");
   if (localMF[FSI_GHOST_MAC_MF+data_dir]->nComp()!=nparts_def*AMREX_SPACEDIM)
    amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nComp() invalid");
  }

  if (localMF[LEVELPC_MF]->nComp()!=nmat*(1+AMREX_SPACEDIM))
   amrex::Error("localMF[LEVELPC_MF]->nComp()!=nmat*(1+AMREX_SPACEDIM)");

  debug_ngrow(LEVELPC_MF,1,37);
  debug_ngrow(SLOPE_RECON_MF,1,37);
  debug_ngrow(delta_MF,0,37);
  debug_ngrow(MASKCOEF_MF,1,36);
  debug_ngrow(DEN_RECON_MF,1,37);
  debug_ngrow(VELADVECT_MF,1,37);

  if (localMF[DEN_RECON_MF]->nComp()!=nmat*num_state_material)
   amrex::Error("localMF[DEN_RECON_MF]->nComp()!=nmat*num_state_material");
  if (localMF[VELADVECT_MF]->nComp()!=num_materials_vel*AMREX_SPACEDIM)
   amrex::Error("localMF[VELADVECT_MF]->nComp()!=materials_vel*AMREX_SPACEDIM");

  int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;
  int nstate=num_materials_vel*(AMREX_SPACEDIM+1)+
   nmat*(num_state_material+ngeom_raw)+1;

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

  int project_option_visc=3;

  int fluxvel_index=0;
  int fluxden_index=AMREX_SPACEDIM;

  const Real* dx = geom.CellSize();
  const Box& domain = geom.Domain();
  const int* domlo = domain.loVect(); 
  const int* domhi = domain.hiVect();

  int rzflag=0;
  if (geom.IsRZ())
   rzflag=1;
  else if (geom.IsCartesian())
   rzflag=0;
  else if (geom.IsCYLINDRICAL())
   rzflag=3;
  else
   amrex::Error("CoordSys bust 20");

  if (init_fluxes==1) {

   int operation_flag=7;
   if (localMF[SEM_FLUXREG_MF]->nComp()!=AMREX_SPACEDIM*nfluxSEM)
    amrex::Error("localMF[SEM_FLUXREG_MF]->nComp() invalid8");

   int mm_areafrac_index=CONSERVE_FLUXES_MF;
   int mm_cell_areafrac_index=SLOPE_RECON_MF;
    //(ml,mr,2) frac_pair(ml,mr),dist_pair(ml,mr)
   int nfacefrac=nmat*nmat*2; 
    // im_inside,im_outside,3+sdim -->
    //   area, dist_to_line, dist, line normal.
   int ncellfrac=nmat*nmat*(3+AMREX_SPACEDIM);

   for (int dirloc=0;dirloc<AMREX_SPACEDIM;dirloc++)
    debug_ngrow(mm_areafrac_index+dirloc,0,111);
   debug_ngrow(mm_cell_areafrac_index,0,114);

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

    debug_ngrow(UMAC_MF+dir,0,115);

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
     FArrayBox& xgp=(*localMF[COARSE_FINE_FLUX_MF+dir])[mfi];  

     FArrayBox& xfacemm=(*localMF[mm_areafrac_index+dir])[mfi];  
     FArrayBox& xcellmm=(*localMF[mm_cell_areafrac_index])[mfi];  

     FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi];  

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

     Vector<int> velbc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);
     int dcomp=num_materials_vel*(AMREX_SPACEDIM+1);
     Vector<int> denbc=getBCArray(State_Type,gridno,dcomp,
      nmat*num_state_material);

     int energyflag=0;
     int local_enable_spectral=enable_spectral;
     int num_materials_face=1;
     int simple_AMR_BC_flag=0;
     int ncomp_xp=nfluxSEM;
     int ncomp_xgp=nfluxSEM;

     int tid_current=ns_thread();
     if ((tid_current<0)||(tid_current>=thread_class::nthreads))
      amrex::Error("tid_current invalid");
     thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     // in SEM_scalar_advection
     // advect: rho u, rho, temperature (non conservatively)
     // FORT_CELL_TO_MAC in LEVELSET_3D.F90
     FORT_CELL_TO_MAC(
      &ncomp_xp,
      &ncomp_xgp,
      &simple_AMR_BC_flag,
      &nsolveMM_FACE,
      &num_materials_face,
      &tileloop,
      &dir,
      &operation_flag, // 7
      &energyflag,
      &visc_coef, //beta
      &visc_coef,
      &face_flag,
      temperature_primitive_variable.dataPtr(),
      &local_enable_spectral,
      &fluxvel_index,
      &fluxden_index,
      &facevel_index,
      &facecut_index,
      &icefacecut_index,
      &curv_index,
      &conservative_tension_force,
      &conservative_div_uu,
      &pforce_index,
      &faceden_index,
      &icemask_index,
      &massface_index,
      &vofface_index,
      &nfluxSEM, // ncphys (nflux for advection)
      &make_interface_incomp,
      override_density.dataPtr(),
      &solvability_projection,
      denbc.dataPtr(),  // presbc
      velbc.dataPtr(),  
      &slab_step,
      &dt_slab,  // CELL_TO_MAC
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
      xfacemm.dataPtr(),ARLIM(xfacemm.loVect()),ARLIM(xfacemm.hiVect()),
      xcellmm.dataPtr(),ARLIM(xcellmm.loVect()),ARLIM(xcellmm.hiVect()),
      reconfab.dataPtr(),ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
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
      &rzflag,domlo,domhi, 
      &nmat,
      &nparts,
      &nparts_def,
      im_solid_map_ptr,
      added_weight.dataPtr(),
      blob_array.dataPtr(),
      &blob_array_size,
      &num_elements_blobclass,
      &num_colors,
      &nten,
      &nfacefrac,
      &ncellfrac,
      &project_option_visc,
      &SEM_upwind,
      &SEM_advection_algorithm);
    }   // mfi
} // omp
    ns_reconcile_d_num(84);
   } // dir=0..sdim-1

  } else if (init_fluxes==0) {

   if (bfact>=1) {

    debug_ngrow(UMAC_MF,0,116);
    debug_ngrow(UMAC_MF+1,0,117);
    debug_ngrow(UMAC_MF+AMREX_SPACEDIM-1,0,118);

    MultiFab& S_new=get_new_data(State_Type,slab_step+1);
    MultiFab& S_old=get_new_data(State_Type,slab_step);

    MultiFab* rhs=new MultiFab(grids,dmap,nfluxSEM,0,
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

     FArrayBox& solxfab=(*localMF[FSI_GHOST_MAC_MF])[mfi];
     FArrayBox& solyfab=(*localMF[FSI_GHOST_MAC_MF+1])[mfi];
     FArrayBox& solzfab=(*localMF[FSI_GHOST_MAC_MF+AMREX_SPACEDIM-1])[mfi];

     FArrayBox& levelpcfab=(*localMF[LEVELPC_MF])[mfi];
     FArrayBox& slopefab=(*localMF[SLOPE_RECON_MF])[mfi];

     FArrayBox& deltafab=(*localMF[delta_MF])[mfi];
     int deltacomp=0;
     if (source_term==1) {
      deltacomp=0;
     } else if (source_term==0) {
      if ((slab_step>=0)&&(slab_step<ns_time_order)) {
       deltacomp=slab_step*nstate_SDC;
      } else
       amrex::Error("slab_step invalid");
     } else
      amrex::Error("source_term invalid");

     Vector<int> velbc=getBCArray(State_Type,gridno,0,
      num_materials_vel*AMREX_SPACEDIM);
     int dcomp=num_materials_vel*(AMREX_SPACEDIM+1);
     Vector<int> denbc=getBCArray(State_Type,gridno,dcomp,
      nmat*num_state_material);

     int operation_flag=6; // advection
     int energyflag=advect_iter;
     int nsolve=nfluxSEM;
     int homflag=source_term;
     int local_enable_spectral=projection_enable_spectral;
     int num_materials_face=1;
     int use_VOF_weight=0;

     int ncomp_denold=nmat*num_state_material;
     int ncomp_veldest=snewfab.nComp();
     int ncomp_dendest=snewfab.nComp()-dcomp;

     int tid_current=ns_thread();
     if ((tid_current<0)||(tid_current>=thread_class::nthreads))
      amrex::Error("tid_current invalid");
     thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

      // called from: NavierStokes::SEM_scalar_advection
     FORT_MAC_TO_CELL(
      &nsolveMM_FACE,
      &num_materials_face,
      &ns_time_order,
      &divu_outer_sweeps,
      &num_divu_outer_sweeps,
      &operation_flag, // 6=advection
      &energyflag,
      temperature_primitive_variable.dataPtr(),
      &nmat,
      &nparts,
      &nparts_def,
      im_solid_map_ptr,
      added_weight.dataPtr(),
      &nten,
      &level, 
      &finest_level,
      &face_flag,
      &make_interface_incomp,
      &solvability_projection,
      &project_option_visc,
      &local_enable_spectral,
      &fluxvel_index,
      &fluxden_index,
      &facevel_index,
      &facecut_index,
      &icefacecut_index,
      &curv_index,
      &conservative_tension_force,
      &conservative_div_uu,
      &pforce_index,
      &faceden_index,
      &icemask_index,
      &massface_index,
      &vofface_index,
      &nfluxSEM, // ncphys
      velbc.dataPtr(),
      denbc.dataPtr(),  // presbc
      &prev_time_slab, 
      &slab_step,
      &dt_slab,  // MAC_TO_CELL
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
      snewfab.dataPtr(dcomp), // dendest
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
      S_old_fab.dataPtr(dcomp), // denold
      ARLIM(S_old_fab.loVect()),ARLIM(S_old_fab.hiVect()),
      S_old_fab.dataPtr(), // ustar
      ARLIM(S_old_fab.loVect()),ARLIM(S_old_fab.hiVect()),
      slopefab.dataPtr(), // recon
      ARLIM(slopefab.loVect()),ARLIM(slopefab.hiVect()),
      velfab.dataPtr(), // mdot
      ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
      S_old_fab.dataPtr(), // maskdivres
      ARLIM(S_old_fab.loVect()),ARLIM(S_old_fab.hiVect()),
      S_old_fab.dataPtr(), // maskres
      ARLIM(S_old_fab.loVect()),ARLIM(S_old_fab.hiVect()),
      &SDC_outer_sweeps,
      &homflag,
      &use_VOF_weight,
      &nsolve,
      &ncomp_denold,
      &ncomp_veldest,
      &ncomp_dendest,
      &SEM_advection_algorithm);

     if (1==0) {
      std::cout << "SEM_scalar_advect c++ level,finest_level " << 
       level << ' ' << finest_level << '\n';
      int interior_only=1;
      tecplot_debug(snewfab,xlo,fablo,fabhi,dx,-1,0,0,
       AMREX_SPACEDIM,interior_only);
     }

    } // mfi
}// omp
    ns_reconcile_d_num(85);

    // rhs=div(uF)  (except for temperature: rhs=u dot grad Theta)
    if (source_term==0) {

     if ((slab_step>=0)&&(slab_step<ns_time_order)) {

      if ((advect_iter==1)&&
          (divu_outer_sweeps+1==num_divu_outer_sweeps)) {
       int deltacomp=slab_step*nstate_SDC;
       MultiFab::Copy(*localMF[stableF_MF],*rhs,0,deltacomp,nfluxSEM,0);
      } else if (advect_iter==0) {
       // do nothing
      } else
       amrex::Error("advect_iter invalid");

     } else
      amrex::Error("source_term invalid");

    } else if (source_term==1) {

     if ((slab_step==0)&&(SDC_outer_sweeps==0)) {
      // this is ok: F_advect(t^n)
     } else if ((slab_step==0)&&(SDC_outer_sweeps!=0)) {
      amrex::Error("F_advect(t^n) already init when SDC_outer_sweeps==0");
     } 

     if ((slab_step>=0)&&(slab_step<=ns_time_order)) {
      int deltacomp=slab_step*nstate_SDC;
      MultiFab::Copy(*localMF[spectralF_MF],*rhs,0,deltacomp,nfluxSEM,0);
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
  amrex::Error("enable_specral or ns_time_order invalid");

} // subroutine SEM_scalar_advection

// Lagrangian solid info lives at t=t^n 
// order_direct_split=base_step mod 2
// must go from finest level to coarsest.
void 
NavierStokes::split_scalar_advection() { 
 
 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();
 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;
 int normdir_here=normdir_direct_split[dir_absolute_direct_split];
 if ((normdir_here<0)||(normdir_here>=AMREX_SPACEDIM))
  amrex::Error("normdir_here invalid");

 int bfact=parent->Space_blockingFactor(level);
 int bfact_f=bfact;
 if (level<finest_level)
  bfact_f=parent->Space_blockingFactor(level+1);

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 int nsolve=1;
 int nsolveMM=nsolve;
 int nsolveMM_FACE=nsolveMM;

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

 int ngrow=2;
 int mac_grow=1; 
 int ngrow_mac_old=0;

 if (face_flag==0) {
  // do nothing
 } else if (face_flag==1) {
  ngrow=3;
  ngrow_mac_old=2;
  mac_grow=2;
 } else
  amrex::Error("face_flag invalid 4");

  // vof,ref centroid,order,slope,intercept  x nmat
 VOF_Recon_resize(ngrow,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,ngrow,36);
 resize_maskfiner(ngrow,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,ngrow,36);
 debug_ngrow(VOF_LS_PREV_TIME_MF,1,38);
 if (localMF[VOF_LS_PREV_TIME_MF]->nComp()!=2*nmat)
  amrex::Error("vof ls prev time invalid ncomp");

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 int ncomp_state=S_new.nComp();
 if (ncomp_state!=num_materials_vel*(AMREX_SPACEDIM+1)+
     nmat*(num_state_material+ngeom_raw)+1)
  amrex::Error("ncomp_state invalid");

 MultiFab& LS_new=get_new_data(LS_Type,slab_step+1);
 if (LS_new.nComp()!=nmat*(AMREX_SPACEDIM+1))
  amrex::Error("LS_new ncomp invalid");

 const Real* dx = geom.CellSize();
 const Box& domain = geom.Domain();
 const int* domlo = domain.loVect();
 const int* domhi = domain.hiVect();

 int scomp_mofvars=num_materials_vel*(AMREX_SPACEDIM+1)+nmat*num_state_material;
 Vector<int> dombc(2*AMREX_SPACEDIM);
 const BCRec& descbc = get_desc_lst()[State_Type].getBC(scomp_mofvars);
 const int* b_rec=descbc.vect();
 for (int m=0;m<2*AMREX_SPACEDIM;m++)
  dombc[m]=b_rec[m];

 vel_time_slab=prev_time_slab;
 if (divu_outer_sweeps==0) 
  vel_time_slab=prev_time_slab;
 else if (divu_outer_sweeps>0)
  vel_time_slab=cur_time_slab;
 else
  amrex::Error("divu_outer_sweeps invalid");

  // in: split_scalar_advection
 getStateDen_localMF(DEN_RECON_MF,ngrow,advect_time_slab);

 getStateTensor_localMF(TENSOR_RECON_MF,1,0,
   num_materials_viscoelastic*NUM_TENSOR_TYPE+AMREX_SPACEDIM,
   advect_time_slab);

 MultiFab& Tensor_new=get_new_data(Tensor_Type,slab_step+1);

 getStateDist_localMF(LS_RECON_MF,1,advect_time_slab,10);

   // the pressure from before will be copied to the new pressure.
 getState_localMF(VELADVECT_MF,ngrow,0,
  num_materials_vel*(AMREX_SPACEDIM+1),
  advect_time_slab); 

   // in: split_scalar_advection
 new_localMF(CONSERVE_FLUXES_MF+normdir_here,nmat,0,normdir_here);
 setVal_localMF(CONSERVE_FLUXES_MF+normdir_here,0.0,0,nmat,0);

  // average down volume fraction fluxes 
 if (level<finest_level) {
  if ((bfact==1)&&(bfact_f==1)) {
   NavierStokes& ns_fine=getLevel(level+1);
   vofflux_sum(normdir_here,
     *localMF[CONSERVE_FLUXES_MF+normdir_here],
     *ns_fine.localMF[CONSERVE_FLUXES_MF+normdir_here]);
  }
 }

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  getStateMAC_localMF(UMACOLD_MF+dir,ngrow_mac_old,dir,
    0,nsolveMM_FACE,advect_time_slab);
 } // dir

 if ((dir_absolute_direct_split<0)||
     (dir_absolute_direct_split>=AMREX_SPACEDIM))
  amrex::Error("dir_absolute_direct_split invalid");

 if (dir_absolute_direct_split==0) {

  int unsplit_displacement=0;
  prepare_displacement(mac_grow,unsplit_displacement);

 } else if ((dir_absolute_direct_split>=1)&&
            (dir_absolute_direct_split<AMREX_SPACEDIM)) {
  // do nothing
 } else
  amrex::Error("dir_absolute_direct_split invalid");

 int vofrecon_ncomp=localMF[SLOPE_RECON_MF]->nComp();
 if (vofrecon_ncomp!=nmat*ngeom_recon)
   amrex::Error("recon ncomp bust");

 int den_recon_ncomp=localMF[DEN_RECON_MF]->nComp();
 if (den_recon_ncomp!=num_state_material*nmat)
   amrex::Error("den_recon invalid");

 int LS_recon_ncomp=localMF[LS_RECON_MF]->nComp();
 if (LS_recon_ncomp!=nmat*(1+AMREX_SPACEDIM))
   amrex::Error("LS_recon invalid");

 debug_ngrow(LS_RECON_MF,1,40);

 resize_mask_nbr(ngrow);
 debug_ngrow(MASK_NBR_MF,ngrow,28); 

 MultiFab* umac_new[AMREX_SPACEDIM];
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  umac_new[dir]=&get_new_data(Umac_Type+dir,slab_step+1);
 }


 int ngrid=grids.size();

 int nc_conserve=AMREX_SPACEDIM+nmat*num_state_material;
 MultiFab* conserve=new MultiFab(grids,dmap,nc_conserve,ngrow,
	MFInfo().SetTag("conserve"),FArrayBoxFactory());

 int iden_base=AMREX_SPACEDIM;
 int itensor_base=iden_base+nmat*num_state_material;
 int imof_base=itensor_base+num_materials_viscoelastic*NUM_TENSOR_TYPE+
	 AMREX_SPACEDIM;
 int iLS_base=imof_base+nmat*ngeom_raw;
 int iFtarget_base=iLS_base+nmat;
 int iden_mom_base=iFtarget_base+nmat;
 int nc_bucket=iden_mom_base+nmat;

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
  FArrayBox& velfab=(*localMF[VELADVECT_MF])[mfi];

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // in: GODUNOV_3D.F90
    // note for evaporation: density_air = total density of ambient gas and
    //   vapor mixture.
    // note: Y_vapor = mass_fraction=mass_vapor/mass_total
    // mass_vapor=V_vapor * den_vapor=F_vapor * V_total * den_vapor
    // mass_total=V_total * den_total
    // so, Y_vapor=F_vapor * V_total * den_vapor/(V_total * den_total)=
    //     F_vapor * den_vapor/den_total 
  FORT_BUILD_CONSERVE( 
   &iden_base,
   override_density.dataPtr(),
   temperature_primitive_variable.dataPtr(),
   consfab.dataPtr(),ARLIM(consfab.loVect()),ARLIM(consfab.hiVect()),
   denfab.dataPtr(),ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
   velfab.dataPtr(),ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
   tilelo,tilehi,
   fablo,fabhi,&bfact,
   &nmat,&ngrow,
   &normdir_here,
   &nc_conserve,
   &den_recon_ncomp);
 }  // mfi
} // omp
 ns_reconcile_d_num(86);

 MultiFab* momslope=new MultiFab(grids,dmap,nc_conserve,1,
  MFInfo().SetTag("momslope"),FArrayBoxFactory());

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

  const Real* xlo = grid_loc[gridno].lo();

  FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi];
  FArrayBox& consfab=(*conserve)[mfi];
  FArrayBox& slopefab=(*momslope)[mfi];
  FArrayBox& masknbrfab=(*localMF[MASK_NBR_MF])[mfi];

  Vector<int> velbc=getBCArray(State_Type,gridno,normdir_here,1);

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  FORT_BUILD_SLOPES( 
   masknbrfab.dataPtr(),
   ARLIM(masknbrfab.loVect()),ARLIM(masknbrfab.hiVect()),
   reconfab.dataPtr(),ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
   consfab.dataPtr(),ARLIM(consfab.loVect()),ARLIM(consfab.hiVect()),
   slopefab.dataPtr(),
   ARLIM(slopefab.loVect()),ARLIM(slopefab.hiVect()),
   &nc_conserve,
   &nmat, 
   tilelo,tilehi,
   fablo,fabhi,&bfact, // slopes conserved vars
   &level,
   &finest_level,
   velbc.dataPtr(),
   xlo,dx,
   &normdir_here,
   &ngrow,
   advection_order.dataPtr(), 
   density_advection_order.dataPtr(), 
   &slope_limiter_option);

 }  // mfi
}// omp
 ns_reconcile_d_num(87);

 MultiFab* xvof[AMREX_SPACEDIM];
 MultiFab* xvel[AMREX_SPACEDIM]; // xvel
 MultiFab* xvelslope[AMREX_SPACEDIM]; // xvelslope,xcen
 MultiFab* side_bucket_mom[AMREX_SPACEDIM];
 MultiFab* side_bucket_mass[AMREX_SPACEDIM];

  // (dir-1)*2*nmat + (side-1)*nmat + im
 int nrefine_vof=2*nmat*AMREX_SPACEDIM;

  // (veldir-1)*2*nmat*sdim + (side-1)*nmat*sdim + (im-1)*sdim+dir
 int nrefine_cen=2*nmat*AMREX_SPACEDIM*AMREX_SPACEDIM;

 if (face_flag==0) {

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   xvof[dir]=localMF[SLOPE_RECON_MF];
   xvel[dir]=localMF[SLOPE_RECON_MF];
   xvelslope[dir]=localMF[SLOPE_RECON_MF];
   side_bucket_mom[dir]=localMF[SLOPE_RECON_MF];
   side_bucket_mass[dir]=localMF[SLOPE_RECON_MF];
  }

 } else if (face_flag==1) {

  MultiFab* vofF=new MultiFab(grids,dmap,nrefine_vof,ngrow,
	MFInfo().SetTag("vofF"),FArrayBoxFactory());

   // linear expansion: 
   // 1. find slopes u_x
   // 2. (rho u)(x)_m = rho_m (u + u_x (x - x_m_centroid))
  MultiFab* cenF=new MultiFab(grids,dmap,nrefine_cen,ngrow,
	MFInfo().SetTag("cenF"),FArrayBoxFactory());
  MultiFab* massF=new MultiFab(grids,dmap,nrefine_vof,ngrow,
	MFInfo().SetTag("massF"),FArrayBoxFactory());

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(vofF->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(*vofF,use_tiling); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();

   const Real* xlo = grid_loc[gridno].lo();

   FArrayBox& slopefab=(*localMF[SLOPE_RECON_MF])[mfi];
   FArrayBox& denstatefab=(*localMF[DEN_RECON_MF])[mfi];

   FArrayBox& vofFfab=(*vofF)[mfi];
   FArrayBox& cenFfab=(*cenF)[mfi];
   FArrayBox& massFfab=(*massF)[mfi];

   int tessellate=0;

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // in: LEVELSET_3D.F90
    // centroid in absolute coordinates.
   FORT_BUILD_SEMIREFINEVOF(
    &tid_current,
    &tessellate,
    &ngrow,
    &nrefine_vof,
    &nrefine_cen,
    &nten,
    spec_material_id_AMBIENT.dataPtr(),
    mass_fraction_id.dataPtr(),
    species_evaporation_density.dataPtr(),
    cavitation_vapor_density.dataPtr(),
    override_density.dataPtr(),
    xlo,dx,
    slopefab.dataPtr(),
    ARLIM(slopefab.loVect()),ARLIM(slopefab.hiVect()),
    denstatefab.dataPtr(),
    ARLIM(denstatefab.loVect()),ARLIM(denstatefab.hiVect()),
    vofFfab.dataPtr(),ARLIM(vofFfab.loVect()),ARLIM(vofFfab.hiVect()),
    cenFfab.dataPtr(),ARLIM(cenFfab.loVect()),ARLIM(cenFfab.hiVect()),
    massFfab.dataPtr(),ARLIM(massFfab.loVect()),ARLIM(massFfab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,
    &bfact,
    &nmat,
    &level,&finest_level);
  }  // mfi
}// omp
  ns_reconcile_d_num(88);

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   xvof[dir]=new MultiFab(state[Umac_Type+dir].boxArray(),dmap,nmat,
     ngrow_mac_old,MFInfo().SetTag("xvof"),FArrayBoxFactory());

   xvel[dir]=new MultiFab(state[Umac_Type+dir].boxArray(),dmap,1,
     ngrow_mac_old,MFInfo().SetTag("xvel"),FArrayBoxFactory());
    // xvelslope,xcen
   xvelslope[dir]=new MultiFab(state[Umac_Type+dir].boxArray(),dmap,1+nmat,
     ngrow_mac_old,MFInfo().SetTag("xvelslope"),FArrayBoxFactory());

    //ncomp=2 ngrow=1
   side_bucket_mom[dir]=new MultiFab(grids,dmap,2,1,
	MFInfo().SetTag("side_bucket_mom"),FArrayBoxFactory());
    //scomp=0 ncomp=2 ngrow=1
   side_bucket_mom[dir]->setVal(0.0,0,2,1);

    //ncomp=2 ngrow=1
   side_bucket_mass[dir]=new MultiFab(grids,dmap,2,1,
	MFInfo().SetTag("side_bucket_mass"),FArrayBoxFactory());
    //scomp=0 ncomp=2 ngrow=1
   side_bucket_mass[dir]->setVal(0.0,0,2,1);
  }  // dir 

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
    FArrayBox& xvoffab=(*xvof[veldir-1])[mfi];
    FArrayBox& xvelfab=(*xvel[veldir-1])[mfi]; // xvelleft,xvelright
    FArrayBox& xvelslopefab=(*xvelslope[veldir-1])[mfi]; // xvelslope,xcen
    FArrayBox& vofFfab=(*vofF)[mfi];
    FArrayBox& cenFfab=(*cenF)[mfi];
    int unsplit_advection=0;

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    FORT_BUILD_MACVOF( 
     &unsplit_advection,
     &nsolveMM_FACE,
     &level,
     &finest_level,
     &normdir_here,
     &nrefine_vof,
     &nrefine_cen,
     vofFfab.dataPtr(),ARLIM(vofFfab.loVect()),ARLIM(vofFfab.hiVect()),
     cenFfab.dataPtr(),ARLIM(cenFfab.loVect()),ARLIM(cenFfab.hiVect()),
     xmac_old.dataPtr(),ARLIM(xmac_old.loVect()),ARLIM(xmac_old.hiVect()),
     xvoffab.dataPtr(),ARLIM(xvoffab.loVect()),ARLIM(xvoffab.hiVect()),
     xvelfab.dataPtr(),ARLIM(xvelfab.loVect()),ARLIM(xvelfab.hiVect()),
     xvelslopefab.dataPtr(),
     ARLIM(xvelslopefab.loVect()),ARLIM(xvelslopefab.hiVect()),
     xlo,dx,
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,
     &nmat,&ngrow, 
     &ngrow_mac_old,&veldir);
   }  // mfi
}// omp
   ns_reconcile_d_num(89);

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

    FArrayBox& xvoffab=(*xvof[veldir-1])[mfi];
    FArrayBox& xvelfab=(*xvel[veldir-1])[mfi]; 
    FArrayBox& xvelslopefab=(*xvelslope[veldir-1])[mfi]; //xvelslope,xcen
    FArrayBox& masknbrfab=(*localMF[MASK_NBR_MF])[mfi];

    Vector<int> velbc=getBCArray(State_Type,gridno,normdir_here,1);
    int slopedir=veldir-1;

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    FORT_BUILD_SLOPES_FACE( 
     masknbrfab.dataPtr(),
     ARLIM(masknbrfab.loVect()),ARLIM(masknbrfab.hiVect()),
      // "vfrac"
     xvoffab.dataPtr(),ARLIM(xvoffab.loVect()),ARLIM(xvoffab.hiVect()),
      // "slsrc"
     xvelfab.dataPtr(),ARLIM(xvelfab.loVect()),ARLIM(xvelfab.hiVect()),
      // "sldst"
     xvelslopefab.dataPtr(),
     ARLIM(xvelslopefab.loVect()),ARLIM(xvelslopefab.hiVect()),
     &nmat, 
     tilelo,tilehi,
     fablo,fabhi,&bfact,
     &level,
     &finest_level,
     velbc.dataPtr(),
     xlo,dx,
     &normdir_here,
     &slopedir,
     &ngrow_mac_old,
     advection_order.dataPtr(), 
     &slope_limiter_option);

   }  // mfi
} // omp
   ns_reconcile_d_num(90);

  } // veldir=1..sdim

  delete vofF;
  delete cenF;
  delete massF;

 } else
  amrex::Error("face_flag invalid 5");

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

  // velocity and pressure
 int scomp_init=0;
 int ncomp_init=num_materials_vel*(AMREX_SPACEDIM+1); 
 S_new.setVal(0.0,scomp_init,ncomp_init,1);

 for (int im=0;im<nmat;im++) {
  if (ns_is_rigid(im)==0) {
   scomp_init=num_materials_vel*(AMREX_SPACEDIM+1)+im*num_state_material;
   ncomp_init=num_state_material;
   S_new.setVal(0.0,scomp_init,ncomp_init,1);
   scomp_init=num_materials_vel*(AMREX_SPACEDIM+1)+nmat*num_state_material+
     im*ngeom_raw;
   ncomp_init=ngeom_raw;
   S_new.setVal(0.0,scomp_init,ncomp_init,1);
   LS_new.setVal(0.0,im,1,1);
  } else if (ns_is_rigid(im)==1) {
   if (solidheat_flag==0) {  // thermal diffuse in solid (default)
    scomp_init=num_materials_vel*(AMREX_SPACEDIM+1)+im*num_state_material+1;
    ncomp_init=1;
    S_new.setVal(0.0,scomp_init,ncomp_init,1);
   } else if (solidheat_flag==2) { // Neumann
    // do nothing
   } else if (solidheat_flag==1) { // dirichlet
    // do nothing
   } else
    amrex::Error("solidheat_flag invalid");
  } else
   amrex::Error("ns_is_rigid(im) invalid");
 } // im=0..nmat-1

 if ((num_materials_viscoelastic>=0)&&
     (num_materials_viscoelastic<=nmat)) {
  Tensor_new.setVal(0.0,0,num_materials_viscoelastic*NUM_TENSOR_TYPE+
		  AMREX_SPACEDIM,1);
 } else
  amrex::Error("num_materials_viscoelastic invalid");

 int ntensor=Tensor_new.nComp();

 MultiFab* cons_cor;
 if (EILE_flag==0) {  // Sussman and Puckett
  cons_cor=localMF[CONS_COR_MF];
 } else if ((EILE_flag==-1)||  // Weymouth and Yue
            (EILE_flag==1)||   // EI-LE
            (EILE_flag==2)||   // always EI
            (EILE_flag==3)) {  // always LE
  cons_cor=localMF[DEN_RECON_MF];
 } else
  amrex::Error("EILE_flag invalid");

 if (dir_absolute_direct_split==0) {
   // initialize the error indicator to be 0.0
  S_new.setVal(0.0,ncomp_state-1,1,1);

  if (EILE_flag==0) {
   if (cons_cor->nComp()!=nmat)
    amrex::Error("cons_cor->nComp()!=nmat");
   if (cons_cor->nGrow()!=1)
    amrex::Error("cons_cor->nGrow()!=1");
   cons_cor->setVal(0.0,0,nmat,1);
  } else if ((EILE_flag==-1)||
             (EILE_flag==1)||
             (EILE_flag==2)||
             (EILE_flag==3)) {
   // do nothing
  } else
   amrex::Error("EILE_flag invalid");
    
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

  FArrayBox& unode=(*localMF[MAC_VELOCITY_MF+normdir_here])[mfi];
  if (unode.nComp()!=nsolveMM_FACE*(AMREX_SPACEDIM+1))
   amrex::Error("unode has invalid ncomp");

    // this is the original data
  FArrayBox& LSfab=(*localMF[LS_RECON_MF])[mfi];
  FArrayBox& denfab=(*localMF[DEN_RECON_MF])[mfi];
  FArrayBox& tenfab=(*localMF[TENSOR_RECON_MF])[mfi];
  FArrayBox& velfab=(*localMF[VELADVECT_MF])[mfi];

    // this is the slope data
  FArrayBox& vofslopefab=(*localMF[SLOPE_RECON_MF])[mfi];

  FArrayBox& vofls0fab=(*localMF[VOF_LS_PREV_TIME_MF])[mfi];
  FArrayBox& corfab=(*cons_cor)[mfi];

  int dencomp=num_materials_vel*(AMREX_SPACEDIM+1);
  int mofcomp=dencomp+nmat*num_state_material;
  int errcomp=mofcomp+nmat*ngeom_raw;

  Vector<int> velbc=getBCArray(State_Type,gridno,normdir_here,1);

     // this is the result
  FArrayBox& destfab=S_new[mfi];
  FArrayBox& tennewfab=Tensor_new[mfi];
  FArrayBox& LSdestfab=LS_new[mfi];

  FArrayBox& consfab=(*conserve)[mfi];

  FArrayBox& xvelfab=(*xvel[0])[mfi]; 
  FArrayBox& yvelfab=(*xvel[1])[mfi];
  FArrayBox& zvelfab=(*xvel[AMREX_SPACEDIM-1])[mfi];

  FArrayBox& xvelslopefab=(*xvelslope[0])[mfi];  // xvelslope,xcen
  FArrayBox& yvelslopefab=(*xvelslope[1])[mfi];
  FArrayBox& zvelslopefab=(*xvelslope[AMREX_SPACEDIM-1])[mfi];

  FArrayBox& slopefab=(*momslope)[mfi];

  FArrayBox& xmomside=(*side_bucket_mom[0])[mfi];
  FArrayBox& ymomside=(*side_bucket_mom[1])[mfi];
  FArrayBox& zmomside=(*side_bucket_mom[AMREX_SPACEDIM-1])[mfi];

  FArrayBox& xmassside=(*side_bucket_mass[0])[mfi];
  FArrayBox& ymassside=(*side_bucket_mass[1])[mfi];
  FArrayBox& zmassside=(*side_bucket_mass[AMREX_SPACEDIM-1])[mfi];

  FArrayBox& ucellfab=(*localMF[CELL_VELOCITY_MF])[mfi];

  FArrayBox& vofflux=(*localMF[CONSERVE_FLUXES_MF+normdir_here])[mfi];

  prescribed_vel_time_slab=0.5*(prev_time_slab+cur_time_slab);

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // solid distance function and solid moments are not modified.
   // solid temperature is modified only if solidheat_flag==0.
  FORT_VFRAC_SPLIT(
   &nsolveMM_FACE,
   &nprocessed[tid_current],
   &tid_current,
   &make_interface_incomp,
   added_weight.dataPtr(),
   density_floor.dataPtr(),
   density_ceiling.dataPtr(),
   &solidheat_flag, //0==diffuse in solid 1==dirichlet 2==neumann
   temperature_primitive_variable.dataPtr(),
   &dencomp,&mofcomp,&errcomp,
   latent_heat.dataPtr(),
   freezing_model.dataPtr(),
   distribute_from_target.dataPtr(),
   &nten,
   &face_flag,
   override_density.dataPtr(),
   velbc.dataPtr(),
   &EILE_flag,
   &VOF_reflux,
   &dir_absolute_direct_split,
   &normdir_here,
   tilelo,tilehi,
   fablo,fabhi,
   &bfact,
   &bfact_f,
   &dt_slab, // VFRAC_SPLIT
   &prev_time_slab,
   &prescribed_vel_time_slab,
   vofflux.dataPtr(),
   ARLIM(vofflux.loVect()),ARLIM(vofflux.hiVect()),
     // this is the original data
   LSfab.dataPtr(),
   ARLIM(LSfab.loVect()),ARLIM(LSfab.hiVect()),
   denfab.dataPtr(),
   ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
   tenfab.dataPtr(),
   ARLIM(tenfab.loVect()),ARLIM(tenfab.hiVect()),
   velfab.dataPtr(),
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
   ucellfab.dataPtr(),ARLIM(ucellfab.loVect()),ARLIM(ucellfab.hiVect()),
   vofls0fab.dataPtr(),ARLIM(vofls0fab.loVect()),ARLIM(vofls0fab.hiVect()),
   corfab.dataPtr(),ARLIM(corfab.loVect()),ARLIM(corfab.hiVect()),
   maskfab.dataPtr(),ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
   masknbrfab.dataPtr(),
   ARLIM(masknbrfab.loVect()),ARLIM(masknbrfab.hiVect()),
   unode.dataPtr(),ARLIM(unode.loVect()),ARLIM(unode.hiVect()),
   xlo,dx,
    // local variables
   consfab.dataPtr(),ARLIM(consfab.loVect()),ARLIM(consfab.hiVect()),
    // xvelleft,xvelright
   xvelfab.dataPtr(),ARLIM(xvelfab.loVect()),ARLIM(xvelfab.hiVect()),
   yvelfab.dataPtr(),ARLIM(yvelfab.loVect()),ARLIM(yvelfab.hiVect()),
   zvelfab.dataPtr(),ARLIM(zvelfab.loVect()),ARLIM(zvelfab.hiVect()),
    // xvelslope,xcen
   xvelslopefab.dataPtr(),
   ARLIM(xvelslopefab.loVect()),ARLIM(xvelslopefab.hiVect()),
   yvelslopefab.dataPtr(),
   ARLIM(yvelslopefab.loVect()),ARLIM(yvelslopefab.hiVect()),
   zvelslopefab.dataPtr(),
   ARLIM(zvelslopefab.loVect()),ARLIM(zvelslopefab.hiVect()),
   slopefab.dataPtr(),ARLIM(slopefab.loVect()),ARLIM(slopefab.hiVect()),
   xmomside.dataPtr(),ARLIM(xmomside.loVect()),ARLIM(xmomside.hiVect()),
   ymomside.dataPtr(),ARLIM(ymomside.loVect()),ARLIM(ymomside.hiVect()),
   zmomside.dataPtr(),ARLIM(zmomside.loVect()),ARLIM(zmomside.hiVect()),
   xmassside.dataPtr(),ARLIM(xmassside.loVect()),ARLIM(xmassside.hiVect()),
   ymassside.dataPtr(),ARLIM(ymassside.loVect()),ARLIM(ymassside.hiVect()),
   zmassside.dataPtr(),ARLIM(zmassside.loVect()),ARLIM(zmassside.hiVect()),
   &ngrow,
   &ngrow_mac_old,
   &nc_conserve,
   &iden_base,
   &nmat,
   &map_forward_direct_split[normdir_here],
   &vofrecon_ncomp,
   &den_recon_ncomp,
   &ncomp_state,
   &ntensor,
   &nc_bucket,
   &nrefine_vof,
   &verbose,
   &gridno,&ngrid,
   &level,
   &finest_level,
   dombc.dataPtr(), 
   domlo,domhi);

 }  // mfi
} // omp
 ns_reconcile_d_num(91);

 for (int tid=1;tid<thread_class::nthreads;tid++) {
  nprocessed[0]+=nprocessed[tid];
 }
 ParallelDescriptor::ReduceIntSum(nprocessed[0]);

 if (profile_debug==1) {
  double profile_time_end=ParallelDescriptor::second();
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "level= " << level << '\n';
   std::cout << "nprocessed= " << nprocessed[0] << '\n';
   std::cout << "profile VFRAC_SPLIT time = " << 
     profile_time_end-profile_time_start << '\n';
  }
 }

 if (face_flag==0) {
  // do nothing
 } else if (face_flag==1) {

  MultiFab* mask_unsplit=new MultiFab(grids,dmap,1,1,
	MFInfo().SetTag("mask_unsplit"),FArrayBoxFactory());
  mask_unsplit->setVal(1.0);

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
   FArrayBox& maskfab=(*localMF[MASKCOEF_MF])[mfi];
   FArrayBox& maskunsplitfab=(*mask_unsplit)[mfi];

   FArrayBox& xmomside=(*side_bucket_mom[0])[mfi];
   FArrayBox& ymomside=(*side_bucket_mom[1])[mfi];
   FArrayBox& zmomside=(*side_bucket_mom[AMREX_SPACEDIM-1])[mfi];

   FArrayBox& xmassside=(*side_bucket_mass[0])[mfi];
   FArrayBox& ymassside=(*side_bucket_mass[1])[mfi];
   FArrayBox& zmassside=(*side_bucket_mass[AMREX_SPACEDIM-1])[mfi];

   FArrayBox& xmac_new=(*umac_new[0])[mfi];
   FArrayBox& ymac_new=(*umac_new[1])[mfi];
   FArrayBox& zmac_new=(*umac_new[AMREX_SPACEDIM-1])[mfi];

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   FORT_BUILD_NEWMAC(
    &nsolveMM_FACE,
    &normdir_here,
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    xmomside.dataPtr(),ARLIM(xmomside.loVect()),ARLIM(xmomside.hiVect()),
    ymomside.dataPtr(),ARLIM(ymomside.loVect()),ARLIM(ymomside.hiVect()),
    zmomside.dataPtr(),ARLIM(zmomside.loVect()),ARLIM(zmomside.hiVect()),
    xmassside.dataPtr(),ARLIM(xmassside.loVect()),ARLIM(xmassside.hiVect()),
    ymassside.dataPtr(),ARLIM(ymassside.loVect()),ARLIM(ymassside.hiVect()),
    zmassside.dataPtr(),ARLIM(zmassside.loVect()),ARLIM(zmassside.hiVect()),
    xmac_new.dataPtr(),ARLIM(xmac_new.loVect()),ARLIM(xmac_new.hiVect()),
    ymac_new.dataPtr(),ARLIM(ymac_new.loVect()),ARLIM(ymac_new.hiVect()),
    zmac_new.dataPtr(),ARLIM(zmac_new.loVect()),ARLIM(zmac_new.hiVect()),
    maskfab.dataPtr(),ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
    maskunsplitfab.dataPtr(),
    ARLIM(maskunsplitfab.loVect()),ARLIM(maskunsplitfab.hiVect()),
    xlo,dx,
    &nmat,
    &level,
    &finest_level);

  }  // mfi
} // omp
  ns_reconcile_d_num(92);

  delete mask_unsplit;
 } else
  amrex::Error("face_flag invalid 6");

 delete conserve;
 delete momslope;
 
 if (face_flag==0) {
  // do nothing
 } else if (face_flag==1) {
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   delete xvof[dir];
   delete xvel[dir];
   delete xvelslope[dir];
   delete side_bucket_mom[dir];
   delete side_bucket_mass[dir];
  }
 } else
  amrex::Error("face_flag invalid 7");

 delete_localMF(VELADVECT_MF,1);
 delete_localMF(DEN_RECON_MF,1);

 if ((num_materials_viscoelastic>=0)&&
     (num_materials_viscoelastic<=nmat)) {
  delete_localMF(TENSOR_RECON_MF,1);
 } else
  amrex::Error("num_materials_viscoelastic invalid");
 
 delete_localMF(LS_RECON_MF,1);
 
 if (stokes_flow==0) {
  // do nothing
 } else if (stokes_flow==1) {
  MultiFab& S_old=get_new_data(State_Type,slab_step);
  MultiFab::Copy(S_new,S_old,0,0,AMREX_SPACEDIM,1);
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   MultiFab::Copy(*umac_new[dir],*localMF[UMACOLD_MF+dir],0,0,1,0);
  }
 } else
  amrex::Error("stokes_flow invalid");

 delete_localMF(UMACOLD_MF,AMREX_SPACEDIM);
 
 if ((level>=0)&&(level<finest_level)) {

  int spectral_override=1;

  if (face_flag==1) {
   avgDownMacState(spectral_override);
  } else if (face_flag==0) {
   // do nothing
  } else
   amrex::Error("face_flag invalid 8");
 
  avgDown(LS_Type,0,nmat,0);
  MOFavgDown();
     // velocity and pressure
  avgDown(State_Type,0,num_materials_vel*(AMREX_SPACEDIM+1),1);
  int scomp_den=num_materials_vel*(AMREX_SPACEDIM+1);
  avgDown(State_Type,scomp_den,num_state_material*nmat,1);
  if ((num_materials_viscoelastic>=0)&&
      (num_materials_viscoelastic<=nmat)) {
    // spectral_override==0 => always low order
   avgDown(Tensor_Type,0,num_materials_viscoelastic*NUM_TENSOR_TYPE+
		   AMREX_SPACEDIM,0);
  } else
   amrex::Error("num_materials_viscoelastic invalid");

 } else if (level==finest_level) {
  // do nothing
 } else
  amrex::Error("level invalid23");

}  // subroutine split_scalar_advection



// Lagrangian solid info lives at t=t^n 
// order_direct_split=base_step mod 2
// must go from finest level to coarsest.
void 
NavierStokes::unsplit_scalar_advection() { 
 
 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();
 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

 int bfact=parent->Space_blockingFactor(level);
 int bfact_f=bfact;
 if (level<finest_level)
  bfact_f=parent->Space_blockingFactor(level+1);

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 int nsolve=1;
 int nsolveMM=nsolve;
 int nsolveMM_FACE=nsolveMM;

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 if ((order_direct_split==0)||
     (order_direct_split==1)) {
  // do nothing
 } else
  amrex::Error("order_direct_split invalid");

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else {
  std::cout << "SDC_outer_sweeps= " << SDC_outer_sweeps << '\n';
  amrex::Error("SDC_outer_sweeps invalid");
 }

 if (dir_absolute_direct_split==0) {
  // do nothing
 } else
  amrex::Error("dir_absolute_direct_split invalid");

 int ngrow=1;
 int mac_grow=1; 
 int ngrow_mac_old=0;

 if (face_flag==0) {
  // do nothing
 } else if (face_flag==1) {
  ngrow++;
  mac_grow++;
  ngrow_mac_old++;
 } else
  amrex::Error("face_flag invalid 4");

  // vof,ref centroid,order,slope,intercept  x nmat
 VOF_Recon_resize(ngrow,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,ngrow,36);
 resize_maskfiner(ngrow,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,ngrow,36);
 debug_ngrow(VOF_LS_PREV_TIME_MF,1,38);
 if (localMF[VOF_LS_PREV_TIME_MF]->nComp()!=2*nmat)
  amrex::Error("vof ls prev time invalid ncomp");

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 int ncomp_state=S_new.nComp();
 if (ncomp_state!=num_materials_vel*(AMREX_SPACEDIM+1)+
     nmat*(num_state_material+ngeom_raw)+1)
  amrex::Error("ncomp_state invalid");

 MultiFab& LS_new=get_new_data(LS_Type,slab_step+1);
 if (LS_new.nComp()!=nmat*(AMREX_SPACEDIM+1))
  amrex::Error("LS_new ncomp invalid");

 const Real* dx = geom.CellSize();
 const Box& domain = geom.Domain();
 const int* domlo = domain.loVect();
 const int* domhi = domain.hiVect();

 int scomp_mofvars=num_materials_vel*(AMREX_SPACEDIM+1)+nmat*num_state_material;
 Vector<int> dombc(2*AMREX_SPACEDIM);
 const BCRec& descbc = get_desc_lst()[State_Type].getBC(scomp_mofvars);
 const int* b_rec=descbc.vect();
 for (int m=0;m<2*AMREX_SPACEDIM;m++)
  dombc[m]=b_rec[m];

 vel_time_slab=prev_time_slab;
 if (divu_outer_sweeps==0) 
  vel_time_slab=prev_time_slab;
 else if (divu_outer_sweeps>0)
  vel_time_slab=cur_time_slab;
 else
  amrex::Error("divu_outer_sweeps invalid");

  // in: unsplit_scalar_advection
 getStateDen_localMF(DEN_RECON_MF,ngrow,advect_time_slab);

 getStateTensor_localMF(TENSOR_RECON_MF,1,0,
   num_materials_viscoelastic*NUM_TENSOR_TYPE+
   AMREX_SPACEDIM,advect_time_slab);

 MultiFab& Tensor_new=get_new_data(Tensor_Type,slab_step+1);

 getStateDist_localMF(LS_RECON_MF,1,advect_time_slab,10);

   // the pressure from before will be copied to the new pressure.
   // VELADVECT_MF and UMACOLD_MF deleted at the end of this routine.
 getState_localMF(VELADVECT_MF,ngrow,0,
  num_materials_vel*(AMREX_SPACEDIM+1),
  advect_time_slab); 

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  getStateMAC_localMF(UMACOLD_MF+dir,ngrow_mac_old,dir,
    0,nsolveMM_FACE,advect_time_slab);
 } // dir

 int unsplit_displacement=1;
  // MAC_VELOCITY_MF deleted towards the end of 
  //   NavierStokes::nonlinear_advection
 prepare_displacement(mac_grow,unsplit_displacement);

 int vofrecon_ncomp=localMF[SLOPE_RECON_MF]->nComp();
 if (vofrecon_ncomp!=nmat*ngeom_recon)
   amrex::Error("recon ncomp bust");

 int den_recon_ncomp=localMF[DEN_RECON_MF]->nComp();
 if (den_recon_ncomp!=num_state_material*nmat)
   amrex::Error("den_recon invalid");

 int LS_recon_ncomp=localMF[LS_RECON_MF]->nComp();
 if (LS_recon_ncomp!=nmat*(1+AMREX_SPACEDIM))
   amrex::Error("LS_recon invalid");

 debug_ngrow(LS_RECON_MF,1,40);

 resize_mask_nbr(ngrow);
 debug_ngrow(MASK_NBR_MF,ngrow,28); 

 MultiFab* umac_new[AMREX_SPACEDIM];
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  umac_new[dir]=&get_new_data(Umac_Type+dir,slab_step+1);
 }

 int ngrid=grids.size();

 int nc_conserve=AMREX_SPACEDIM+nmat*num_state_material;
 MultiFab* conserve=new MultiFab(grids,dmap,nc_conserve,ngrow,
  MFInfo().SetTag("conserve"),FArrayBoxFactory());

 MultiFab* mask_unsplit=new MultiFab(grids,dmap,1,1,
  MFInfo().SetTag("mask_unsplit"),FArrayBoxFactory());

 int iden_base=AMREX_SPACEDIM;
 int itensor_base=iden_base+nmat*num_state_material;
 int imof_base=itensor_base+num_materials_viscoelastic*NUM_TENSOR_TYPE+
	 AMREX_SPACEDIM;
 int iLS_base=imof_base+nmat*ngeom_raw;
 int iFtarget_base=iLS_base+nmat;
 int iden_mom_base=iFtarget_base+nmat;
 int nc_bucket=iden_mom_base+nmat;

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
  FArrayBox& velfab=(*localMF[VELADVECT_MF])[mfi];

  int normdir_unsplit=0; // not used

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  FORT_BUILD_CONSERVE( 
   &iden_base,
   override_density.dataPtr(),
   temperature_primitive_variable.dataPtr(),
   consfab.dataPtr(),ARLIM(consfab.loVect()),ARLIM(consfab.hiVect()),
   denfab.dataPtr(),ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
   velfab.dataPtr(),ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
   tilelo,tilehi,
   fablo,fabhi,&bfact,
   &nmat,&ngrow,
   &normdir_unsplit,
   &nc_conserve,
   &den_recon_ncomp);
 }  // mfi
} // omp
 ns_reconcile_d_num(93);

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

  const Real* xlo = grid_loc[gridno].lo();

  FArrayBox& maskfab=(*mask_unsplit)[mfi];
  FArrayBox& vofls0fab=(*localMF[VOF_LS_PREV_TIME_MF])[mfi];

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  FORT_BUILD_MASK_UNSPLIT( 
   &unsplit_flag,
   &make_interface_incomp,
   xlo,dx,
   maskfab.dataPtr(),ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
   vofls0fab.dataPtr(),ARLIM(vofls0fab.loVect()),ARLIM(vofls0fab.hiVect()),
   tilelo,tilehi,
   fablo,fabhi,
   &bfact,
   &nmat,
   &level,&finest_level);
 }  // mfi
} // omp
 ns_reconcile_d_num(94);

 MultiFab* xvof[AMREX_SPACEDIM];
 MultiFab* xvel[AMREX_SPACEDIM]; // xvel
 MultiFab* side_bucket_mom[AMREX_SPACEDIM];
 MultiFab* side_bucket_mass[AMREX_SPACEDIM];

  // (dir-1)*2*nmat + (side-1)*nmat + im
 int nrefine_vof=2*nmat*AMREX_SPACEDIM;

  // (veldir-1)*2*nmat*sdim + (side-1)*nmat*sdim + (im-1)*sdim+dir
 int nrefine_cen=2*nmat*AMREX_SPACEDIM*AMREX_SPACEDIM;

 if (face_flag==0) {

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   xvof[dir]=localMF[SLOPE_RECON_MF];
   xvel[dir]=localMF[SLOPE_RECON_MF];
   side_bucket_mom[dir]=localMF[SLOPE_RECON_MF];
   side_bucket_mass[dir]=localMF[SLOPE_RECON_MF];
  }

 } else if (face_flag==1) {

  MultiFab* vofF=new MultiFab(grids,dmap,nrefine_vof,ngrow,
	MFInfo().SetTag("vofF"),FArrayBoxFactory());

   // linear expansion: 
   // 1. find slopes u_x
   // 2. (rho u)(x)_m = rho_m (u + u_x (x - x_m_centroid))
  MultiFab* cenF=new MultiFab(grids,dmap,nrefine_cen,ngrow,
	MFInfo().SetTag("cenF"),FArrayBoxFactory());
  MultiFab* massF=new MultiFab(grids,dmap,nrefine_vof,ngrow,
		MFInfo().SetTag("massF"),FArrayBoxFactory());

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(vofF->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(*vofF,use_tiling); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();

   const Real* xlo = grid_loc[gridno].lo();

   FArrayBox& slopefab=(*localMF[SLOPE_RECON_MF])[mfi];
   FArrayBox& denstatefab=(*localMF[DEN_RECON_MF])[mfi];

   FArrayBox& vofFfab=(*vofF)[mfi];
   FArrayBox& cenFfab=(*cenF)[mfi];
   FArrayBox& massFfab=(*massF)[mfi];

   int tessellate=0;

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // centroid in absolute coordinates.
   FORT_BUILD_SEMIREFINEVOF(
    &tid_current,
    &tessellate,
    &ngrow,
    &nrefine_vof,
    &nrefine_cen,
    &nten,
    spec_material_id_AMBIENT.dataPtr(),
    mass_fraction_id.dataPtr(),
    species_evaporation_density.dataPtr(),
    cavitation_vapor_density.dataPtr(),
    override_density.dataPtr(),
    xlo,dx,
    slopefab.dataPtr(),
    ARLIM(slopefab.loVect()),ARLIM(slopefab.hiVect()),
    denstatefab.dataPtr(),
    ARLIM(denstatefab.loVect()),ARLIM(denstatefab.hiVect()),
    vofFfab.dataPtr(),ARLIM(vofFfab.loVect()),ARLIM(vofFfab.hiVect()),
    cenFfab.dataPtr(),ARLIM(cenFfab.loVect()),ARLIM(cenFfab.hiVect()),
    massFfab.dataPtr(),ARLIM(massFfab.loVect()),ARLIM(massFfab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,
    &bfact,
    &nmat,
    &level,&finest_level);
  }  // mfi
}// omp
  ns_reconcile_d_num(95);

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   xvof[dir]=new MultiFab(state[Umac_Type+dir].boxArray(),dmap,nmat,
     ngrow_mac_old,MFInfo().SetTag("xvof"),FArrayBoxFactory());

   xvel[dir]=new MultiFab(state[Umac_Type+dir].boxArray(),dmap,1,
     ngrow_mac_old,MFInfo().SetTag("xvel"),FArrayBoxFactory());

    //ncomp=2 ngrow=1
   side_bucket_mom[dir]=new MultiFab(grids,dmap,2,1,
		MFInfo().SetTag("side_bucket_mom"),FArrayBoxFactory());
    //scomp=0 ncomp=2 ngrow=1
   side_bucket_mom[dir]->setVal(0.0,0,2,1);

    //ncomp=2 ngrow=1
   side_bucket_mass[dir]=new MultiFab(grids,dmap,2,1,
	MFInfo().SetTag("side_bucket_mass"),FArrayBoxFactory());
    //scomp=0 ncomp=2 ngrow=1
   side_bucket_mass[dir]->setVal(0.0,0,2,1);
  }  // dir 

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
    FArrayBox& xvoffab=(*xvof[veldir-1])[mfi];
    FArrayBox& xvelfab=(*xvel[veldir-1])[mfi]; // xvelleft,xvelright
    FArrayBox& vofFfab=(*vofF)[mfi];
    FArrayBox& cenFfab=(*cenF)[mfi];
    int unsplit_advection=1;
    int normdir_unsplit=0; // not used

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    FORT_BUILD_MACVOF( 
     &unsplit_advection,
     &nsolveMM_FACE,
     &level,
     &finest_level,
     &normdir_unsplit,
     &nrefine_vof,
     &nrefine_cen,
     vofFfab.dataPtr(),ARLIM(vofFfab.loVect()),ARLIM(vofFfab.hiVect()),
     cenFfab.dataPtr(),ARLIM(cenFfab.loVect()),ARLIM(cenFfab.hiVect()),
     xmac_old.dataPtr(),ARLIM(xmac_old.loVect()),ARLIM(xmac_old.hiVect()),
     xvoffab.dataPtr(),ARLIM(xvoffab.loVect()),ARLIM(xvoffab.hiVect()),
     xvelfab.dataPtr(),
     ARLIM(xvelfab.loVect()),ARLIM(xvelfab.hiVect()),
     xvelfab.dataPtr(),  // xvelslopefab
     ARLIM(xvelfab.loVect()),ARLIM(xvelfab.hiVect()),
     xlo,dx,
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,
     &nmat,&ngrow, 
     &ngrow_mac_old,&veldir);
   }  // mfi
}// omp
   ns_reconcile_d_num(96);

  } // veldir=1..sdim

  delete vofF;
  delete cenF;
  delete massF;

 } else
  amrex::Error("face_flag invalid 5");

 Vector<int> nprocessed;
 nprocessed.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  nprocessed[tid]=0.0;
 }

 double profile_time_start=0.0;
 if (profile_debug==1) {
  profile_time_start=ParallelDescriptor::second();
 }

 int ntensor=Tensor_new.nComp();

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

  FArrayBox& maskunsplitfab=(*mask_unsplit)[mfi];

  FArrayBox& unode=(*localMF[MAC_VELOCITY_MF])[mfi];
  if (unode.nComp()!=nsolveMM_FACE*(1+AMREX_SPACEDIM))
   amrex::Error("unode has invalid ncomp");
  FArrayBox& vnode=(*localMF[MAC_VELOCITY_MF+1])[mfi];
  if (vnode.nComp()!=nsolveMM_FACE*(1+AMREX_SPACEDIM))
   amrex::Error("vnode has invalid ncomp");
  FArrayBox& wnode=(*localMF[MAC_VELOCITY_MF+AMREX_SPACEDIM-1])[mfi];
  if (wnode.nComp()!=nsolveMM_FACE*(1+AMREX_SPACEDIM))
   amrex::Error("wnode has invalid ncomp");

    // this is the original data
  FArrayBox& LSfab=(*localMF[LS_RECON_MF])[mfi];
  FArrayBox& denfab=(*localMF[DEN_RECON_MF])[mfi];

  FArrayBox& tenfab=(*localMF[TENSOR_RECON_MF])[mfi];

  FArrayBox& velfab=(*localMF[VELADVECT_MF])[mfi];

    // this is the slope data
  FArrayBox& vofslopefab=(*localMF[SLOPE_RECON_MF])[mfi];

  FArrayBox& vofls0fab=(*localMF[VOF_LS_PREV_TIME_MF])[mfi];

  int dencomp=num_materials_vel*(AMREX_SPACEDIM+1);
  int mofcomp=dencomp+nmat*num_state_material;
  int errcomp=mofcomp+nmat*ngeom_raw;

  Vector<int> velbc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);

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

  FArrayBox& ucellfab=(*localMF[CELL_VELOCITY_MF])[mfi];

  prescribed_vel_time_slab=0.5*(prev_time_slab+cur_time_slab);

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // solid distance function and solid moments are not modified.
   // solid temperature is modified only if solidheat_flag==0.
  FORT_VFRAC_UNSPLIT(
   &unsplit_flag,  
   &nsolveMM_FACE,
   &nprocessed[tid_current],
   &tid_current,
   &make_interface_incomp,
   added_weight.dataPtr(),
   density_floor.dataPtr(),
   density_ceiling.dataPtr(),
   &solidheat_flag, //0==diffuse in solid 1==dirichlet 2==neumann
   temperature_primitive_variable.dataPtr(),
   &dencomp,&mofcomp,&errcomp,
   latent_heat.dataPtr(),
   freezing_model.dataPtr(),
   distribute_from_target.dataPtr(),
   &nten,
   &face_flag,
   override_density.dataPtr(),
   velbc.dataPtr(),
   tilelo,tilehi,
   fablo,fabhi,
   &bfact,
   &bfact_f,
   &dt_slab, // VFRAC_UNSPLIT
   &prev_time_slab,
   &prescribed_vel_time_slab,
     // this is the original data
   LSfab.dataPtr(),
   ARLIM(LSfab.loVect()),ARLIM(LSfab.hiVect()),
   denfab.dataPtr(),
   ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
   tenfab.dataPtr(),
   ARLIM(tenfab.loVect()),ARLIM(tenfab.hiVect()),
   velfab.dataPtr(),
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
   ucellfab.dataPtr(),ARLIM(ucellfab.loVect()),ARLIM(ucellfab.hiVect()),
   vofls0fab.dataPtr(),ARLIM(vofls0fab.loVect()),ARLIM(vofls0fab.hiVect()),
   maskfab.dataPtr(),ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
   masknbrfab.dataPtr(),
   ARLIM(masknbrfab.loVect()),ARLIM(masknbrfab.hiVect()),
   maskunsplitfab.dataPtr(),
   ARLIM(maskunsplitfab.loVect()),ARLIM(maskunsplitfab.hiVect()),
   unode.dataPtr(),ARLIM(unode.loVect()),ARLIM(unode.hiVect()),
   vnode.dataPtr(),ARLIM(vnode.loVect()),ARLIM(vnode.hiVect()),
   wnode.dataPtr(),ARLIM(wnode.loVect()),ARLIM(wnode.hiVect()),
   xlo,dx,
    // local variables
   consfab.dataPtr(),ARLIM(consfab.loVect()),ARLIM(consfab.hiVect()),
    // xvelleft,xvelright
   xvelfab.dataPtr(),ARLIM(xvelfab.loVect()),ARLIM(xvelfab.hiVect()),
   yvelfab.dataPtr(),ARLIM(yvelfab.loVect()),ARLIM(yvelfab.hiVect()),
   zvelfab.dataPtr(),ARLIM(zvelfab.loVect()),ARLIM(zvelfab.hiVect()),
   xmomside.dataPtr(),ARLIM(xmomside.loVect()),ARLIM(xmomside.hiVect()),
   ymomside.dataPtr(),ARLIM(ymomside.loVect()),ARLIM(ymomside.hiVect()),
   zmomside.dataPtr(),ARLIM(zmomside.loVect()),ARLIM(zmomside.hiVect()),
   xmassside.dataPtr(),ARLIM(xmassside.loVect()),ARLIM(xmassside.hiVect()),
   ymassside.dataPtr(),ARLIM(ymassside.loVect()),ARLIM(ymassside.hiVect()),
   zmassside.dataPtr(),ARLIM(zmassside.loVect()),ARLIM(zmassside.hiVect()),
   &ngrow,
   &ngrow_mac_old,
   &nc_conserve,
   &iden_base,
   &nmat,
   &vofrecon_ncomp,
   &den_recon_ncomp,
   &ncomp_state,
   &ntensor,
   &nc_bucket,
   &nrefine_vof,
   &verbose,
   &gridno,&ngrid,
   &level,
   &finest_level,
   dombc.dataPtr(), 
   domlo,domhi);

 }  // mfi
} // omp
 ns_reconcile_d_num(97);

 for (int tid=1;tid<thread_class::nthreads;tid++) {
  nprocessed[0]+=nprocessed[tid];
 }
 ParallelDescriptor::ReduceIntSum(nprocessed[0]);

 if (profile_debug==1) {
  double profile_time_end=ParallelDescriptor::second();
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "level= " << level << '\n';
   std::cout << "nprocessed= " << nprocessed[0] << '\n';
   std::cout << "profile VFRAC_SPLIT time = " << 
     profile_time_end-profile_time_start << '\n';
  }
 }

 if (face_flag==0) {
  // do nothing
 } else if (face_flag==1) {

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
   FArrayBox& maskfab=(*localMF[MASKCOEF_MF])[mfi];

   FArrayBox& xmomside=(*side_bucket_mom[0])[mfi];
   FArrayBox& ymomside=(*side_bucket_mom[1])[mfi];
   FArrayBox& zmomside=(*side_bucket_mom[AMREX_SPACEDIM-1])[mfi];

   FArrayBox& xmassside=(*side_bucket_mass[0])[mfi];
   FArrayBox& ymassside=(*side_bucket_mass[1])[mfi];
   FArrayBox& zmassside=(*side_bucket_mass[AMREX_SPACEDIM-1])[mfi];

   FArrayBox& xmac_new=(*umac_new[0])[mfi];
   FArrayBox& ymac_new=(*umac_new[1])[mfi];
   FArrayBox& zmac_new=(*umac_new[AMREX_SPACEDIM-1])[mfi];

   FArrayBox& maskunsplitfab=(*mask_unsplit)[mfi];

   int normdir_unsplit=0; // not used

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   FORT_BUILD_NEWMAC(
    &nsolveMM_FACE,
    &normdir_unsplit,
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    xmomside.dataPtr(),ARLIM(xmomside.loVect()),ARLIM(xmomside.hiVect()),
    ymomside.dataPtr(),ARLIM(ymomside.loVect()),ARLIM(ymomside.hiVect()),
    zmomside.dataPtr(),ARLIM(zmomside.loVect()),ARLIM(zmomside.hiVect()),
    xmassside.dataPtr(),ARLIM(xmassside.loVect()),ARLIM(xmassside.hiVect()),
    ymassside.dataPtr(),ARLIM(ymassside.loVect()),ARLIM(ymassside.hiVect()),
    zmassside.dataPtr(),ARLIM(zmassside.loVect()),ARLIM(zmassside.hiVect()),
    xmac_new.dataPtr(),ARLIM(xmac_new.loVect()),ARLIM(xmac_new.hiVect()),
    ymac_new.dataPtr(),ARLIM(ymac_new.loVect()),ARLIM(ymac_new.hiVect()),
    zmac_new.dataPtr(),ARLIM(zmac_new.loVect()),ARLIM(zmac_new.hiVect()),
    maskfab.dataPtr(),ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
    maskunsplitfab.dataPtr(),
    ARLIM(maskunsplitfab.loVect()),ARLIM(maskunsplitfab.hiVect()),
    xlo,dx,
    &nmat,
    &level,
    &finest_level);

  }  // mfi
} // omp
  ns_reconcile_d_num(98);
 } else
  amrex::Error("face_flag invalid 6");

 delete conserve;
 delete mask_unsplit;
 
 if (face_flag==0) {
  // do nothing
 } else if (face_flag==1) {
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   delete xvof[dir];
   delete xvel[dir];
   delete side_bucket_mom[dir];
   delete side_bucket_mass[dir];
  }
 } else
  amrex::Error("face_flag invalid 7");

 delete_localMF(VELADVECT_MF,1);
 delete_localMF(DEN_RECON_MF,1);

 if ((num_materials_viscoelastic>=0)&&
     (num_materials_viscoelastic<=nmat)) {
  delete_localMF(TENSOR_RECON_MF,1);
 } else
  amrex::Error("num_materials_viscoelastic invalid");
 
 delete_localMF(LS_RECON_MF,1);
 
 if (stokes_flow==0) {
  // do nothing
 } else if (stokes_flow==1) {
  MultiFab& S_old=get_new_data(State_Type,slab_step);
  MultiFab::Copy(S_new,S_old,0,0,AMREX_SPACEDIM,1);
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   MultiFab::Copy(*umac_new[dir],*localMF[UMACOLD_MF+dir],0,0,1,0);
  }
 } else
  amrex::Error("stokes_flow invalid");

 delete_localMF(UMACOLD_MF,AMREX_SPACEDIM);
 
 if ((level>=0)&&(level<finest_level)) {

  int spectral_override=1;

  if (face_flag==1) {
   avgDownMacState(spectral_override);
  } else if (face_flag==0) {
   // do nothing
  } else
   amrex::Error("face_flag invalid 8");
 
  avgDown(LS_Type,0,nmat,0);
  MOFavgDown();
     // velocity and pressure
  avgDown(State_Type,0,num_materials_vel*(AMREX_SPACEDIM+1),1);
  int scomp_den=num_materials_vel*(AMREX_SPACEDIM+1);
  avgDown(State_Type,scomp_den,num_state_material*nmat,1);
  if ((num_materials_viscoelastic>=0)&&
      (num_materials_viscoelastic<=nmat)) {
    // spectral_override==0 => always low order
   avgDown(Tensor_Type,0,num_materials_viscoelastic*NUM_TENSOR_TYPE+
		   AMREX_SPACEDIM,0);
  } else
   amrex::Error("num_materials_viscoelastic invalid");

 } else if (level==finest_level) {
  // do nothing
 } else
  amrex::Error("level invalid23");

}  // subroutine unsplit_scalar_advection


void
NavierStokes::errorEst (TagBoxArray& tags,int clearval,int tagval,
 int n_error_buf,int ngrow)
{
 
 const int max_level = parent->maxLevel();
 if (level>=max_level)
  amrex::Error("level too big in errorEst");

 int bfact=parent->Space_blockingFactor(level);

 if (n_error_buf<2)
  amrex::Error("amr.n_error_buf<2");

 const int*  domain_lo = geom.Domain().loVect();
 const int*  domain_hi = geom.Domain().hiVect();
 const Real* dx        = geom.CellSize();
 const Real* prob_lo   = geom.ProbLo();

 int local_time_order=parent->Time_blockingFactor();
 Real nudge_time=state[State_Type].slabTime(local_time_order);
 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 int nc_error=S_new.nComp()-1;
 MultiFab* mf = getState(0,nc_error,1,nudge_time); 

 bool use_tiling=ns_tiling;

 if (1==1) {
  use_tiling=false;
 }

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(mf->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  // itags: Vector<int> ("domain"=box dimensions (not whole domain))
  // itags=TagBox::CLEAR in "tags"
  // then itags = older tag values.
 Vector<int>  itags;

 for (MFIter mfi(*mf,use_tiling); mfi.isValid(); ++mfi) {
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

  FArrayBox& snewfab=(*mf)[mfi];
  const Box& snew_box = snewfab.box();
  if (snew_box==tilegrid) {
   // do nothing
  } else
   amrex::Error("FIX errorEST for tiling");

  Real*       dat   = snewfab.dataPtr();
  const int*  dlo   = snew_box.loVect();
  const int*  dhi   = snew_box.hiVect();
  const int   ncomp = snewfab.nComp();

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  FORT_VFRACERROR(
    tptr, ARLIM(tlo), ARLIM(thi), 
    &tagval, &clearval, 
    dat, ARLIM(dlo), ARLIM(dhi),
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    &ncomp, 
    domain_lo, domain_hi,
    dx, xlo, prob_lo, 
    &upper_slab_time, 
    &level,
    &max_level,
    &max_level_two_materials,
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
 ns_reconcile_d_num(99);

 delete mf;
} // subroutine errorEst

// called from volWgtSumALL
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
// integrated_quantities:
//  1..3  : F=integral_boundary sigma dot n + integral rho g
//  4..6  : torque=integral_boundary r x sigma dot n + integral r x rho g  
//  7..9  : moments of inertia
//  10..12: integral rho x
//  13    : integral rho
void 
NavierStokes::GetDragALL(Vector<Real>& integrated_quantities) {

 int finest_level=parent->finestLevel();

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid");

 if (integrated_quantities.size()!=13)
  amrex::Error("invalid size for integrated_quantities");

 for (int iq=0;iq<13;iq++)
  integrated_quantities[iq]=0.0;

  // in: NavierStokes::GetDragALL
 allocate_levelsetLO_ALL(2,LEVELPC_MF);

 int simple_AMR_BC_flag_viscosity=1;
 int do_alloc=1;
 init_gradu_tensorALL(HOLD_VELOCITY_DATA_MF,do_alloc,CELLTENSOR_MF,
  FACETENSOR_MF,
  simple_AMR_BC_flag_viscosity);

 for (int isweep=0;isweep<2;isweep++) {
  for (int ilev=level;ilev<=finest_level;ilev++) {
   NavierStokes& ns_level=getLevel(ilev);
   ns_level.GetDrag(integrated_quantities,isweep);
  }
 }

 delete_array(CELLTENSOR_MF);
 delete_array(FACETENSOR_MF);

  // isweep
 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   for (int iq=0;iq<13;iq++)
    std::cout << "GetDrag  iq= " << iq << " integrated_quantities= " <<
     integrated_quantities[iq] << '\n';
   Real mass=integrated_quantities[12];
   if (mass>0.0) {
    for (int iq=9;iq<9+AMREX_SPACEDIM;iq++)
     std::cout << "COM iq= " << iq << " COM= " <<
      integrated_quantities[iq]/mass << '\n';
   }

  }
 }

} // GetDragALL

// sweep=0: integral (rho x), integral (rho) 
// sweep=1: find force, torque, moments of inertia, center of mass,mass
void
NavierStokes::GetDrag(Vector<Real>& integrated_quantities,int isweep) {

 int bfact=parent->Space_blockingFactor(level);

 int nmat=num_materials;
 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

 int rzflag=0;
 if (geom.IsRZ())
  rzflag=1;
 else if (geom.IsCartesian())
  rzflag=0;
 else if (geom.IsCYLINDRICAL())
  rzflag=3;
 else
  amrex::Error("CoordSys bust 4");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 debug_ngrow(CELLTENSOR_MF,1,45);
 debug_ngrow(FACETENSOR_MF,1,45);
 resize_levelsetLO(2,LEVELPC_MF);
 debug_ngrow(LEVELPC_MF,2,45);
 if (localMF[LEVELPC_MF]->nComp()!=nmat*(1+AMREX_SPACEDIM))
  amrex::Error("levelpc mf has incorrect ncomp");
 VOF_Recon_resize(1,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,1,46);
 debug_ngrow(CELL_VISC_MATERIAL_MF,0,47);
 debug_ngrow(CELL_VISC_MF,1,47);
 resize_metrics(1);
 debug_ngrow(VOLUME_MF,1,48);
 debug_ngrow(DRAG_MF,0,50);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(FACE_VAR_MF+dir,0,2);

 int nparts=im_solid_map.size();
 if ((nparts<0)||(nparts>nmat))
  amrex::Error("nparts invalid");
 Vector<int> im_solid_map_null;
 im_solid_map_null.resize(1);

 int* im_solid_map_ptr;
 int nparts_def=nparts;
 if (nparts==0) {
  im_solid_map_ptr=im_solid_map_null.dataPtr();
  nparts_def=1;
 } else if ((nparts>=1)&&(nparts<=nmat)) {
  im_solid_map_ptr=im_solid_map.dataPtr();
 } else
  amrex::Error("nparts invalid");

 int ntensor=AMREX_SPACEDIM*AMREX_SPACEDIM;
 int ntensorMM=ntensor*num_materials_vel;

 if (localMF[CELLTENSOR_MF]->nComp()!=ntensorMM)
  amrex::Error("localMF[CELLTENSOR_MF]->nComp() invalid");
 if (localMF[FACETENSOR_MF]->nComp()!=ntensorMM)
  amrex::Error("localMF[FACETENSOR_MF]->nComp() invalid");

 int nstate=num_materials_vel*(AMREX_SPACEDIM+1)+
  nmat*(num_state_material+ngeom_raw)+1;
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");

 if (integrated_quantities.size()!=13)
  amrex::Error("integrated_quantities invalid size");
 
 Vector< Vector<Real> > local_integrated_quantities;
 local_integrated_quantities.resize(thread_class::nthreads);

 for (int tid=0;tid<thread_class::nthreads;tid++) {
  local_integrated_quantities[tid].resize(13);
  for (int iq=0;iq<13;iq++)
   local_integrated_quantities[tid][iq]=0.0;
 } // tid

// gear problem: probtype=563, axis_dir=2, 3D
// scale torque by 2 pi vinletgas/60


 if (localMF[DRAG_MF]->nComp()!=4*AMREX_SPACEDIM+1)
  amrex::Error("drag ncomp invalid");

 const Real* dx = geom.CellSize();

 int project_option_combine=3;  // velocity in GetDrag
 int combine_flag=2;
 int hflag=0;
 int combine_idx=-1;  // update state variables
 int update_flux=0;
 combine_state_variable(
  project_option_combine,
  combine_idx,combine_flag,hflag,update_flux);
 project_option_combine=0; // mac velocity
 update_flux=1;
 combine_state_variable(
  project_option_combine,
  combine_idx,combine_flag,hflag,update_flux);

  // p(rho,T) in compressible parts
  // projection pressure in incompressible parts.
 MultiFab* pres=getStatePres(1,cur_time_slab); 

 if (pres->nComp()!=num_materials_vel)
  amrex::Error("pres->nComp() invalid");
 
 MultiFab* vel=getState(2,0,AMREX_SPACEDIM*num_materials_vel,cur_time_slab);
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

 MultiFab* elastic_tensor_mf=den_recon;
 int elastic_ntensor=den_recon->nComp();

 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
  if (localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0)
   amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0");
  if (localMF[FSI_GHOST_MAC_MF+data_dir]->nComp()!=nparts_def*AMREX_SPACEDIM)
   amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nComp() invalid");
 }

 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=nmat)) {

  elastic_ntensor=num_materials_viscoelastic*NUM_TENSOR_TYPE;
  elastic_tensor_mf=
   new MultiFab(grids,dmap,elastic_ntensor,1,
	MFInfo().SetTag("elastic_tensor_mf"),FArrayBoxFactory());
  elastic_tensor_mf->setVal(0.0);
    
  for (int im=0;im<nmat;im++) {

   if (ns_is_rigid(im)==0) {

    if ((elastic_time[im]>0.0)&& 
        (elastic_viscosity[im]>0.0)) {

     int partid=0;
     while ((im_elastic_map[partid]!=im)&&(partid<im_elastic_map.size())) {
      partid++;
     }

     if (partid<im_elastic_map.size()) {
       // in: GetDrag
      make_viscoelastic_tensor(im);
      MultiFab::Copy(*elastic_tensor_mf,*localMF[VISCOTEN_MF],0,
       partid*NUM_TENSOR_TYPE,NUM_TENSOR_TYPE,1);
      delete_localMF(VISCOTEN_MF,1);
     } else
      amrex::Error("partid could not be found: GetDrag");
    } else if ((elastic_time[im]==0.0)||
               (elastic_viscosity[im]==0.0)) {

     if (viscoelastic_model[im]!=0)
      amrex::Error("viscoelastic_model[im]!=0");

    } else
     amrex::Error("elastic_time/elastic_viscosity invalid");

   } else if (ns_is_rigid(im)==1) {
    // do nothing
   } else
    amrex::Error("ns_is_rigid invalid");

  } // im=0..nmat-1
 } else if (num_materials_viscoelastic==0) {
  // do nothing
 } else
  amrex::Error("num_materials_viscoelastic invalid");

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(localMF[DRAG_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*localMF[DRAG_MF],false); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();

  const Real* xlo = grid_loc[gridno].lo();
  Vector<int> velbc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);

  FArrayBox& maskfab=(*mask)[mfi];
  FArrayBox& volfab=(*localMF[VOLUME_MF])[mfi];
  FArrayBox& areax=(*localMF[AREA_MF])[mfi];
  FArrayBox& areay=(*localMF[AREA_MF+1])[mfi];
  FArrayBox& areaz=(*localMF[AREA_MF+AMREX_SPACEDIM-1])[mfi];

  FArrayBox& mufab=(*localMF[CELL_VISC_MF])[mfi];
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
  FArrayBox& elastic_tensor_data=(*elastic_tensor_mf)[mfi];

  Real gravity_normalized=std::abs(gravity);
  if (invert_gravity==1)
   gravity_normalized=-gravity_normalized;
  else if (invert_gravity==0) {
   // do nothing
  } else
   amrex::Error("invert_gravity invalid");

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // hoop stress, centripetal force, coriolis effect still not
   // considered.
   // in: DERIVE_3D.F90
  FORT_GETDRAG(
   &isweep,
   integrated_quantities.dataPtr(),
   local_integrated_quantities[tid_current].dataPtr(),
   &gravity_normalized,
   &gravity_dir,
   &elastic_ntensor,
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
   &facevisc_index,
   &faceheat_index,
   &ncphys,
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
   &rzflag,velbc.dataPtr(),&cur_time_slab,
   &visc_coef,
   &ntensor,
   &ntensorMM,
   &nmat,
   &nparts,
   &nparts_def,
   im_solid_map_ptr);

 } // mfi
} // omp
 ns_reconcile_d_num(100);

 int iqstart=0;
 int iqend=13;
 if (isweep==0) 
  iqstart=9;
 else if (isweep==1)
  iqend=9;
 else
  amrex::Error("isweep invalid");

 for (int tid=1;tid<thread_class::nthreads;tid++) {
  for (int iq=iqstart;iq<iqend;iq++) {
   local_integrated_quantities[0][iq]+=local_integrated_quantities[tid][iq];
  }
 } // tid

 ParallelDescriptor::Barrier();

 for (int iq=iqstart;iq<iqend;iq++) {
  ParallelDescriptor::ReduceRealSum(local_integrated_quantities[0][iq]);
  integrated_quantities[iq]+=local_integrated_quantities[0][iq];
 }

 project_option_combine=3; // velocity in GetDrag
 combine_flag=2;
 hflag=0;
 combine_idx=-1; 
 update_flux=0;
 combine_state_variable(
  project_option_combine,
  combine_idx,combine_flag,hflag,update_flux);
 project_option_combine=0; // mac velocity
 update_flux=1;
 combine_state_variable(
  project_option_combine,
  combine_idx,combine_flag,hflag,update_flux);
 
 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=nmat)) {
  delete elastic_tensor_mf;
 } else if (num_materials_viscoelastic==0) {
  // do nothing
 } else
  amrex::Error("num_materials_viscoelastic invalid");

 delete pres; 
 delete vel; 
 delete den_recon; 
 delete mask; 
 delete mask3; 

} // subroutine GetDrag

void NavierStokes::project_right_hand_side(
  int index_MF,int project_option,int& change_flag) {

 change_flag=0;

 if ((project_option==0)||
     (project_option==1)||
     (project_option==10)||  // sync project due to advection
     (project_option==11)||  // FSI_material_exists (2nd project)
     (project_option==13)||  // FSI_material_exists (1st project)
     (project_option==12)) { // pressure extension

  int finest_level=parent->finestLevel();
  if (finest_level>=0) {

   if (level==0) {

    if (ones_sum_global>=1.0) {

     if (local_solvability_projection==0) {
      if (singular_possible==1) {
       zap_resid_where_singular(index_MF); // multiply by ONES_MF
       change_flag=1;
      } else if (singular_possible==0) {
       // do nothing
      } else
       amrex::Error("singular_possible invalid"); 
     } else if (local_solvability_projection==1) {

      if (singular_possible==1) {
       // rhsnew=rhs H-alpha H
       // 0 =sum rhs H-alpha sum H
       // alpha=sum rhs H / sum H
       zap_resid_where_singular(index_MF); // multiply by ONES_MF
       Real coef;
       dot_productALL_ones(project_option,index_MF,coef);
       coef=-coef/ones_sum_global;
       mf_combine_ones(project_option,index_MF,coef);
       zap_resid_where_singular(index_MF);
       change_flag=1;
      } else
       amrex::Error("singular_possible invalid"); 

     } else
      amrex::Error("local_solvability_projection invalid");

    } else {
     std::cout << "index_MF = " << index_MF << '\n';
     std::cout << "project_option = " << project_option << '\n';
     std::cout << "ones_sum_global = " << ones_sum_global << '\n';
     amrex::Error("ones_sum_global invalid");
    }
   } else
    amrex::Error("level invalid");
  } else
   amrex::Error("finest_level invalid");

 } else if (project_option==2) { // thermal diffusion
  if (singular_possible==0) {
   // do nothing
  } else
   amrex::Error("singular_possible invalid");
 } else if (project_option==3) {  // viscosity
  if (singular_possible==0) {
   // do nothing
  } else
   amrex::Error("singular_possible invalid");
 } else if ((project_option>=100)&&(project_option<100+num_species_var)) {
  if (singular_possible==0) {
   // do nothing
  } else
   amrex::Error("singular_possible invalid");
 } else
  amrex::Error("project_option invalid project_right_hand_side");

} // subroutine project_right_hand_side

void NavierStokes::init_checkerboardLEV(
  int index_MF,int project_option,
  int nsolve,int nsolveMM) {

 bool use_tiling=ns_tiling;

 if ((nsolve!=1)&&(nsolve!=AMREX_SPACEDIM))
  amrex::Error("nsolve invalid");

 debug_ngrow(index_MF,0,51);

 if (localMF[index_MF]->nGrow()==0) {
  // do nothing
 } else
  amrex::Error("localMF[index_MF]->nGrow() invalid");
 if (localMF[index_MF]->nComp()==nsolveMM) {
  // do nothing
 } else
  amrex::Error("localMF[index_MF]->nComp() invalid");

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 int finest_level=parent->finestLevel();
 if (level>finest_level)
  amrex::Error("level too big");

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(localMF[index_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel 
#endif
{
 for (MFIter mfi(*localMF[index_MF],use_tiling); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());

  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();
  int bfact=parent->Space_blockingFactor(level);

  FArrayBox& fab = (*localMF[index_MF])[mfi];

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // in: MG_3D.F90
  FORT_CHECKERBOARD_RB(&nsolveMM,
    fab.dataPtr(),ARLIM(fab.loVect()), ARLIM(fab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,
    &bfact,&bfact);
 } // mfi
} // omp
 ns_reconcile_d_num(101);

} // init_checkerboardLEV

void NavierStokes::init_checkerboardALL(
  int index_MF,int project_option,
  int nsolve,int nsolveMM) {

 if ((project_option==0)||
     (project_option==1)||
     (project_option==10)||  // sync project due to advection
     (project_option==11)||  // FSI_material_exists (2nd project)
     (project_option==13)||  // FSI_material_exists (1st project)
     (project_option==12)) { // pressure extension

  int finest_level=parent->finestLevel();
  if (finest_level>=0) {

   if (level==0) {

    if (ones_sum_global>=1.0) {

     if (local_solvability_projection==0) { // system is nonsingular.
      if (singular_possible==1) {  // some regions might be masked off
       // do nothing
      } else if (singular_possible==0) {
       // do nothing
      } else
       amrex::Error("singular_possible invalid"); 
     } else if (local_solvability_projection==1) { // system is singular

      if (singular_possible==1) { // some parts of domain might be masked off.
       for (int ilev=finest_level;ilev>=level;ilev--) {
        NavierStokes& ns_level=getLevel(ilev);
	ns_level.init_checkerboardLEV(index_MF,project_option,
			nsolve,nsolveMM);
       }
      } else
       amrex::Error("singular_possible invalid"); 

     } else
      amrex::Error("local_solvability_projection invalid");

    } else {
     std::cout << "index_MF = " << index_MF << '\n';
     std::cout << "project_option = " << project_option << '\n';
     std::cout << "ones_sum_global = " << ones_sum_global << '\n';
     amrex::Error("ones_sum_global invalid");
    }
   } else
    amrex::Error("level invalid");
  } else
   amrex::Error("finest_level invalid");

 } else if (project_option==2) { // thermal diffusion
  if (singular_possible==0) {
   // do nothing
  } else
   amrex::Error("singular_possible invalid");
 } else if (project_option==3) {  // viscosity
  if (singular_possible==0) {
   // do nothing
  } else
   amrex::Error("singular_possible invalid");
 } else if ((project_option>=100)&&(project_option<100+num_species_var)) {
  if (singular_possible==0) {
   // do nothing
  } else
   amrex::Error("singular_possible invalid");
 } else
  amrex::Error("project_option invalid project_right_hand_side");

} // subroutine init_checkerboardALL



void NavierStokes::dot_productALL_ones(int project_option,
   int index_MF,Real& result) {

 int finest_level=parent->finestLevel();
 if (level!=0)
  amrex::Error("level=0 in dot_productALL_ones");

 Real tempsum=0.0;
 Real total_sum=0.0;
 int nsolve=1;

 for (int k = 0; k <= finest_level; k++) {
  NavierStokes& ns_level = getLevel(k);
  ns_level.dotSum(
    project_option,
    ns_level.localMF[index_MF],
    ns_level.localMF[ONES_MF], 
    tempsum,nsolve);

  total_sum+=tempsum;
 }

 result=total_sum;
}  // subroutine dot_projectALL_ones


void NavierStokes::zap_resid_where_singular(int index_MF) {

 int finest_level=parent->finestLevel();
 if (level!=0)
  amrex::Error("level=0 in dot_productALL_ones");

 for (int k = 0; k <= finest_level; k++) {
  NavierStokes& ns_level = getLevel(k);
  MultiFab::Multiply(*ns_level.localMF[index_MF],
	*ns_level.localMF[ONES_MF],0,0,1,0);
 }  // k=0..finest_level

}  // subroutine zap_resid_where_singular

void NavierStokes::dot_productALL_ones_size(int project_option,Real& result) {

 int finest_level=parent->finestLevel();
 if (level!=0)
  amrex::Error("level=0 in dot_productALL_ones");

 Real tempsum=0.0;
 Real total_sum=0.0;
 int nsolve=1;

 for (int k = 0; k <= finest_level; k++) {
  NavierStokes& ns_level = getLevel(k);
  ns_level.dotSum(
   project_option,
   ns_level.localMF[ONES_MF], 
   ns_level.localMF[ONES_MF], tempsum, nsolve);

  total_sum+=tempsum;
 }

 result=total_sum;

} // subroutine dot_productALL_ones_size




void NavierStokes::mf_combine_ones(
 int project_option,
 int index_MF,Real& Beta) {

  int nsolve=1;

    // amf_x=amf_x+Beta * ones_mf  (where mask=1)
  int finest_level=parent->finestLevel();
  for (int k = 0; k <= finest_level; ++k) {
    NavierStokes& ns_level = getLevel(k);
    ns_level.levelCombine(
     project_option,
     ns_level.localMF[index_MF],
     ns_level.localMF[ONES_MF],
     ns_level.localMF[index_MF], Beta,nsolve);
  }

}


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
 
 bool use_tiling=ns_tiling;

 int num_materials_face=num_materials_vel;

 if ((nsolve!=1)&&(nsolve!=AMREX_SPACEDIM))
  amrex::Error("nsolve invalid");

 if ((project_option==0)||
     (project_option==1)||
     (project_option==10)||
     (project_option==11)|| //FSI_material_exists (2nd project)
     (project_option==13)|| //FSI_material_exists (1st project)
     (project_option==12)|| //pressure extrapolation
     (project_option==3)) { //viscosity
  if (num_materials_face!=1)
   amrex::Error("num_materials_face invalid");
 } else if ((project_option==2)||  // thermal diffusion
            ((project_option>=100)&&
             (project_option<100+num_species_var))) {
  num_materials_face=num_materials_scalar_solve;
 } else
  amrex::Error("project_option invalid2");

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 if ((num_materials_face!=1)&&
     (num_materials_face!=num_materials))
  amrex::Error("num_materials_face invalid");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(FACE_VAR_MF+dir,0,2);

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,0,51);
 debug_ngrow(DOTMASK_MF,0,51);
 if (localMF[DOTMASK_MF]->nComp()!=num_materials_face)
  amrex::Error("localMF[DOTMASK_MF]->nComp()!=num_materials_face");

 int nsolveMM=nsolve*num_materials_face;

 if (mf1->nComp()!=nsolveMM) {
  std::cout << "mf1_ncomp nsolveMM " << mf1->nComp() << ' ' << nsolveMM << '\n';
  amrex::Error("dotSum: mf1 invalid ncomp");
 }
 if (mf2->nComp()!=nsolveMM)
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
  FArrayBox& dotfab=(*localMF[DOTMASK_MF])[mfi];
  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // in: NAVIERSTOKES_3D.F90
  FORT_SUMDOT(&tsum,
    fab.dataPtr(),ARLIM(fab.loVect()), ARLIM(fab.hiVect()),
    fab2.dataPtr(),ARLIM(fab2.loVect()), ARLIM(fab2.hiVect()),
    dotfab.dataPtr(),ARLIM(dotfab.loVect()),ARLIM(dotfab.hiVect()),
    mfab.dataPtr(),ARLIM(mfab.loVect()),ARLIM(mfab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    &debug_dot_product,
    &level,&gridno,
    &nsolve, 
    &nsolveMM, 
    &num_materials_face);
  sum[tid_current] += tsum;
 } // mfi1
} // omp
 ns_reconcile_d_num(101);
 for (int tid=1;tid<thread_class::nthreads;tid++) {
  sum[0]+=sum[tid];
 }
 ParallelDescriptor::ReduceRealSum(sum[0]);

 result=sum[0];
}


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
 
 bool use_tiling=ns_tiling;

 int num_materials_face=num_materials_vel;

 if ((nsolve!=1)&&(nsolve!=AMREX_SPACEDIM))
  amrex::Error("nsolve invalid");

 if ((project_option==0)||
     (project_option==1)||
     (project_option==10)||
     (project_option==11)|| //FSI_material_exists (2nd project)
     (project_option==13)|| //FSI_material_exists (1st project)
     (project_option==12)|| //pressure extrapolation
     (project_option==3)) { //viscosity
  if (num_materials_face!=1)
   amrex::Error("num_materials_face invalid");
 } else if ((project_option==2)||  // thermal diffusion
            ((project_option>=100)&&
             (project_option<100+num_species_var))) {
  num_materials_face=num_materials_scalar_solve;
 } else
  amrex::Error("project_option invalid2");

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 if ((num_materials_face!=1)&&
     (num_materials_face!=num_materials))
  amrex::Error("num_materials_face invalid");

 int finest_level=parent->finestLevel();
 if (level > finest_level)
  amrex::Error("level too big");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(FACE_VAR_MF+dir,0,2);

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,0,51);
 debug_ngrow(DOTMASK_MF,0,51);
 if (localMF[DOTMASK_MF]->nComp()!=num_materials_face)
  amrex::Error("localMF[DOTMASK_MF]->nComp()!=num_materials_face");

 int nsolveMM=nsolve*num_materials_face;

 if (mfx->nComp()!= nsolveMM)
  amrex::Error("mfx invalid ncomp");
 if (mfy->nComp()!= nsolveMM)
  amrex::Error("mfy invalid ncomp");
 if (mfz->nComp()!= nsolveMM)
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
  FArrayBox& dotfab=(*localMF[DOTMASK_MF])[mfi];
  FArrayBox& fabz = (*mfz)[mfi];
  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  //fabz = fabx + beta * faby
  //in: NAVIERSTOKES_3D.F90
  FORT_FABCOM(
   fabx.dataPtr(),ARLIM(fabx.loVect()),ARLIM(fabx.hiVect()),
   faby.dataPtr(),ARLIM(faby.loVect()),ARLIM(faby.hiVect()),
   dotfab.dataPtr(),ARLIM(dotfab.loVect()),ARLIM(dotfab.hiVect()),
   mfab.dataPtr(),ARLIM(mfab.loVect()),ARLIM(mfab.hiVect()),
   fabz.dataPtr(),ARLIM(fabz.loVect()),ARLIM(fabz.hiVect()),
   &beta,
   tilelo,tilehi,
   fablo,fabhi,&bfact,
   &nsolve,
   &nsolveMM,
   &num_materials_face);
 } // mfi
} // omp
 ns_reconcile_d_num(102);

} // subroutine levelCombine



void NavierStokes::volWgtSum(
  Vector<Real>& result,
  Vector<int>& sumdata_type,
  Vector<int>& sumdata_sweep,
  Vector<Real>& ZZ,Vector<Real>& FF,
  int dirx,int diry,int cut_flag,
  MultiFab* dragmf,int isweep) {

 
 bool use_tiling=ns_tiling;

 int nmat=num_materials;
 int ntensor=AMREX_SPACEDIM*AMREX_SPACEDIM;
 int ntensorMM=ntensor*num_materials_vel;

  // 0 empty
  // F,E  2 x nmat
  // drag (3 comp)
  // min interface location 3 x nmat  (x1,y1,z1   x2,y2,z2  ...)
  // max interface location 3 x nmat  (x1,y1,z1   x2,y2,z2  ...)
  // pressure drag (3 comp)
  // min den,denA 2 x nmat
  // max den,denA 2 x nmat
  // x=0 amplitude
  // centroid 3 x nmat (x1,y1,z1  x2,y2,z2  ... )
  // min dist from centroid  nmat
  // max dist from centroid  nmat
  // mass      nmat
  // momentum  3 x nmat
  // energy    nmat
  // left pressure,right pressure, left wt,right wt
  // kinetic energy derived  nmat
  // LS F  nmat
  // LS centroid 3 x nmat (x1,y1,z1  x2,y2,z2 ... )
  // torque (3 comp)
  // pressure torque (3 comp)
  // perimeter (rasterized) (1 comp)
  // min interface extent on slice (nmat comp)
  // max interface extent on slice (nmat comp)
  // integral of vorticity (3 comp)
  // vort_error (1 comp)
  // vel_error (1 comp)
  // energy_moment (1 comp)

 int filler_comp=0;
 int FE_sum_comp=filler_comp+1;
 int drag_sum_comp=FE_sum_comp+2*nmat;
 int minint_sum_comp=drag_sum_comp+3;
 int maxint_sum_comp=minint_sum_comp+3*nmat;
 int pdrag_sum_comp=maxint_sum_comp+3*nmat;
 int minden_sum_comp=pdrag_sum_comp+3;
 int maxden_sum_comp=minden_sum_comp+2*nmat;
 int xnot_amp_sum_comp=maxden_sum_comp+2*nmat;
 int cen_sum_comp=xnot_amp_sum_comp+1;
 int mincen_sum_comp=cen_sum_comp+3*nmat;
 int maxcen_sum_comp=mincen_sum_comp+nmat;
 int mass_sum_comp=maxcen_sum_comp+nmat;
 int mom_sum_comp=mass_sum_comp+nmat;
 int energy_sum_comp=mom_sum_comp+3*nmat;
 int left_pressure_sum=energy_sum_comp+nmat;
 int kinetic_energy_sum_comp=left_pressure_sum+4;
 int LS_F_sum_comp=kinetic_energy_sum_comp+nmat;
 int LS_cen_sum_comp=LS_F_sum_comp+nmat;
 int torque_sum_comp=LS_cen_sum_comp+3*nmat;
 int ptorque_sum_comp=torque_sum_comp+3;
 int step_perim_sum_comp=ptorque_sum_comp+3;
 int minint_slice=step_perim_sum_comp+1;
 int maxint_slice=minint_slice+nmat;
 int vort_sum_comp=maxint_slice+nmat;
 int vort_error=vort_sum_comp+3; 
 int vel_error=vort_error+1; 
 int energy_moment=vel_error+1; 
 int enstrophy=energy_moment+1; // integral of w dot w
 int user_comp=enstrophy+nmat;
 int total_comp=user_comp+ncomp_sum_int_user; 

 if (total_comp!=result.size())
  amrex::Error("result size invalid");
 if (total_comp!=sumdata_type.size())
  amrex::Error("sumdata_type size invalid");
 if (total_comp!=sumdata_sweep.size())
  amrex::Error("sumdata_sweep size invalid");

 int finest_level=parent->finestLevel();

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 if (localMF[CELLTENSOR_MF]->nComp()!=ntensorMM)
  amrex::Error("localMF[CELLTENSOR_MF]->nComp() invalid");

 MultiFab* den_recon=getStateDen(1,upper_slab_time);  
 int den_ncomp=den_recon->nComp();

 MultiFab* error_heat_map_mf=new MultiFab(grids,dmap,nmat,0,
	MFInfo().SetTag("error_heat_map_mf"),FArrayBoxFactory());
 error_heat_map_mf -> setVal(0.0);

 int project_option_combine=3; // velocity in volWgtSum
 int combine_flag=2;
 int hflag=0;
 int combine_idx=-1;  // update state variables
 int update_flux=0;
 combine_state_variable(
  project_option_combine,
  combine_idx,combine_flag,hflag,update_flux);
 project_option_combine=0; // mac velocity
 update_flux=1;
 combine_state_variable(
  project_option_combine,
  combine_idx,combine_flag,hflag,update_flux);

 resize_maskfiner(2,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,2,53);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(FACE_VAR_MF+dir,0,2);

 VOF_Recon_resize(2,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,2,54);
 debug_ngrow(CELLTENSOR_MF,1,54);

  // velocity and pressure
 MultiFab* vel=getState(1,0,num_materials_vel*(AMREX_SPACEDIM+1),
   upper_slab_time);

 const Real* dx = geom.CellSize();

 int resultsize=result.size();
 if (resultsize!=total_comp)
  amrex::Error("resultsize invalid");

 int NN=ZZ.size()-1;
 if (NN!=FF.size()-1)
  amrex::Error("FF and ZZ fail size sanity check");

 Vector< Vector<Real> > local_result;
 Vector< Vector<Real> > local_ZZ;
 Vector< Vector<Real> > local_FF;

 local_result.resize(thread_class::nthreads);
 local_ZZ.resize(thread_class::nthreads);
 local_FF.resize(thread_class::nthreads);

 for (int tid=0;tid<thread_class::nthreads;tid++) {
  local_result[tid].resize(resultsize);

  for (int isum=0;isum<resultsize;isum++) {
   local_result[tid][isum]=0.0;
   if (sumdata_type[isum]==2) // min
    local_result[tid][isum]=1.0E+15;
   else if (sumdata_type[isum]==3)  // max
    local_result[tid][isum]=-1.0E+15;
   else if (sumdata_type[isum]==1)
    local_result[tid][isum]=0.0;
   else
    amrex::Error("sumdata_type invalid");
  } // isum

  local_ZZ[tid].resize(NN+1);
  local_FF[tid].resize(NN+1);
  for (int iz=0;iz<=NN;iz++) {
   local_ZZ[tid][iz]=0.0;
   local_FF[tid][iz]=0.0;
  }
 }  // tid

 MultiFab* lsmf=getStateDist(2,upper_slab_time,11);

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

   FArrayBox& cellten=(*localMF[CELLTENSOR_MF])[mfi];
   if (cellten.nComp()!=ntensorMM)
    amrex::Error("cellten invalid ncomp");

   FArrayBox& dragfab=(*dragmf)[mfi];
   Real problo[AMREX_SPACEDIM];
   Real probhi[AMREX_SPACEDIM];
   for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
    problo[dir]=geom.ProbLo(dir);
    probhi[dir]=geom.ProbHi(dir);
   }
   FArrayBox& lsfab=(*lsmf)[mfi];
   int bfact=parent->Space_blockingFactor(level);

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   FORT_SUMMASS(
    &tid_current,
    &ncomp_sum_int_user,
    &adapt_quad_depth,
    &slice_dir,
    xslice.dataPtr(),
    problo,probhi, 
    xlo,dx,
    cellten.dataPtr(),ARLIM(cellten.loVect()),ARLIM(cellten.hiVect()),
    lsfab.dataPtr(),ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
    maskSEMfab.dataPtr(),
    ARLIM(maskSEMfab.loVect()),ARLIM(maskSEMfab.hiVect()),
    mfab.dataPtr(),ARLIM(mfab.loVect()),ARLIM(mfab.hiVect()),
    dragfab.dataPtr(),ARLIM(dragfab.loVect()),ARLIM(dragfab.hiVect()),
    reconfab.dataPtr(),ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
    denfab.dataPtr(),ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
    velfab.dataPtr(),ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,
    &bfact,
    &upper_slab_time,
    local_result[tid_current].dataPtr(),
    result.dataPtr(),
    sumdata_type.dataPtr(),
    sumdata_sweep.dataPtr(),
    &resultsize,
    &NN,
    local_ZZ[tid_current].dataPtr(),
    local_FF[tid_current].dataPtr(),
    &dirx,&diry,&cut_flag,
    &nmat,
    &ntensorMM,
    &den_ncomp,
    &isweep);

 } // mfi
} // omp
 ns_reconcile_d_num(103);

 for (int tid=1;tid<thread_class::nthreads;tid++) {

  for (int idest=1;idest<total_comp;idest++) {

   if (((sumdata_sweep[idest]==0)&&(isweep==0))|| //default (update 1st sweep)
       ((sumdata_sweep[idest]==1)&&(isweep==1))) {//(update 2nd sweep) 

    if (sumdata_type[idest]==1) { // reduce real sum (default)
     local_result[0][idest]+=local_result[tid][idest];
    } else if (sumdata_type[idest]==2) { // reduce real min 
     if (local_result[tid][idest]<local_result[0][idest])
      local_result[0][idest]=local_result[tid][idest];
    } else if (sumdata_type[idest]==3) { // reduce real max
     if (local_result[tid][idest]>local_result[0][idest])
      local_result[0][idest]=local_result[tid][idest];
    } else
     amrex::Error("sumdata_type invalid");

   } else if (sumdata_sweep[idest]==0) {
    // do nothing
   } else if (sumdata_sweep[idest]==1) {
    // do nothing
   } else
    amrex::Error("sumdata_sweep invalid");

  } // idest
 
  if (isweep==0) { 
   for (int iz=0;iz<=NN;iz++) {
    if (local_ZZ[tid][iz]>local_ZZ[0][iz])
     local_ZZ[0][iz]=local_ZZ[tid][iz];
    if (local_FF[tid][iz]>local_FF[0][iz])
     local_FF[0][iz]=local_FF[tid][iz];
   } // iz
  }  // isweep==0
 } // tid

 ParallelDescriptor::Barrier();

 project_option_combine=3; // velocity in volWgtSum
 combine_flag=2;
 hflag=0;
 combine_idx=-1;  // update state variables
 update_flux=0;
 combine_state_variable(
  project_option_combine,
  combine_idx,combine_flag,hflag,update_flux);
 project_option_combine=0; // mac velocity
 update_flux=1;
 combine_state_variable(
  project_option_combine,
  combine_idx,combine_flag,hflag,update_flux);

 result[filler_comp]=0.0;

 for (int idest=1;idest<total_comp;idest++) {

  if (((sumdata_sweep[idest]==0)&&(isweep==0))||
      ((sumdata_sweep[idest]==1)&&(isweep==1))) { 

   if (sumdata_type[idest]==1) { // reduce real sum
    ParallelDescriptor::ReduceRealSum(local_result[0][idest]);
    result[idest]=result[idest]+local_result[0][idest];
   } else if (sumdata_type[idest]==2) { // reduce real min
    ParallelDescriptor::ReduceRealMin(local_result[0][idest]);
    if (local_result[0][idest]<result[idest])
     result[idest]=local_result[0][idest];
   } else if (sumdata_type[idest]==3) { // reduce real max
    ParallelDescriptor::ReduceRealMax(local_result[0][idest]);
    if (local_result[0][idest]>result[idest])
     result[idest]=local_result[0][idest];
   } else
    amrex::Error("sumdata_type invalid");

  } else if (sumdata_sweep[idest]==0) {
   // do nothing
  } else if (sumdata_sweep[idest]==1) {
   // do nothing
  } else
   amrex::Error("sumdata_sweep invalid");

 } // idest
 
 if (isweep==0) { 
  for (int iz=0;iz<=NN;iz++) {
   ParallelDescriptor::ReduceRealMax(local_ZZ[0][iz]);
   ParallelDescriptor::ReduceRealMax(local_FF[0][iz]);
   if (local_ZZ[0][iz]>ZZ[iz])
    ZZ[iz]=local_ZZ[0][iz];
   if (local_FF[0][iz]>FF[iz])
    FF[iz]=local_FF[0][iz];
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
     nmat,interior_only);
  }// mfi
  ns_reconcile_d_num(104);
 } // verbose

 delete error_heat_map_mf;
 delete lsmf;
 delete den_recon;
 delete vel;
}  // subroutine volWgtSum


void
NavierStokes::setPlotVariables()
{
  AmrLevel::setPlotVariables();
}

std::string
NavierStokes::thePlotFileType () const
{
    //
    // Increment this whenever the writePlotFile() format changes.
    //
    static const std::string the_plot_file_type("NavierStokes-V1.1");

    return the_plot_file_type;
}

void 
NavierStokes::debug_memory() {

 if (level!=0)
  amrex::Error("level invalid debug_memory");

 if (show_mem==1) {
  std::fflush(NULL);
  ParallelDescriptor::Barrier();
  int proc=ParallelDescriptor::MyProc();
  std::cout << "calling memory status on processor=" << proc << '\n';
  FORT_MEMSTATUS(&proc);
  std::cout << "after calling memory status on processor=" << proc << '\n';
  std::fflush(NULL);
  ParallelDescriptor::Barrier();
 } else if (show_mem!=0)
  amrex::Error("show_mem invalid");

}

void NavierStokes::writeInterfaceReconstruction() {

 debug_ngrow(SLOPE_RECON_MF,1,55);
 if (level!=0)
  amrex::Error("level should be zero");

 int finest_level = parent->finestLevel();
 Vector<int> grids_per_level;
 grids_per_level.resize(finest_level+1);
 for (int ilev=finest_level;ilev>=0;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  grids_per_level[ilev]=ns_level.grids.size();
  ns_level.output_triangles();  // NavierStokes2.cpp: nmat materials at once
 }
 ParallelDescriptor::Barrier();
 if (ParallelDescriptor::IOProcessor()) {
  int nsteps=parent->levelSteps(0);
  int arrdim=finest_level+1;

  int plotint=parent->plotInt();
  int nmat=num_materials;
  for (int im=1;im<=nmat;im++) {
    // in: NAVIERSTOKES_3D.F90
   FORT_COMBINETRIANGLES(grids_per_level.dataPtr(),
    &finest_level,
    &nsteps,
    &im,
    &arrdim,
    &cur_time_slab,
    &plotint);
  } // im=1..nmat

  for (int ipart=0;ipart<NS_ncomp_particles;ipart++) {
    // in: NAVIERSTOKES_3D.F90
   fort_combine_particles(grids_per_level.dataPtr(),
    &finest_level,
    &nsteps,
    &ipart,
    &NS_ncomp_particles,
    &arrdim,
    &cur_time_slab,
    &plotint);
  } // ipart=0..NS_ncomp_particles-1

 }
 ParallelDescriptor::Barrier();

}  // writeInterfaceReconstruction


// VOF_Recon_ALL called before this routine is called.
// init_FSI_GHOST_MAC_MF() called for all relevant 
// levels prior to this routine.
void NavierStokes::writeTECPLOT_File(int do_plot,int do_slice) {

 if (level!=0)
  amrex::Error("level invalid writeTECPLOT_File");

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel != 1");

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid");

 int nsteps=parent->levelSteps(0);

 int nmat=num_materials;
 int ntensor=AMREX_SPACEDIM*AMREX_SPACEDIM;
 int ntensorMM=ntensor*num_materials_vel;

 int finest_level = parent->finestLevel();

 ParallelDescriptor::Barrier();

 debug_ngrow(SLOPE_RECON_MF,1,65);
 if (localMF[SLOPE_RECON_MF]->nComp()!=nmat*ngeom_recon)
  amrex::Error("localMF[SLOPE_RECON_MF]->nComp() invalid");

  // uses "slope_recon" 
 if (do_plot==1) {
  writeInterfaceReconstruction();
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
 int nstate_slice=AMREX_SPACEDIM+AMREX_SPACEDIM+6;
 slice_data.resize(nslice*nstate_slice);
 for (int i=0;i<nslice*nstate_slice;i++)
  slice_data[i]=-1.0e+30;

  // in: NavierStokes::writeTECPLOT_File
 allocate_levelsetLO_ALL(1,LEVELPC_MF);

// HOLD_VELOCITY_DATA_MF not already allocated,
// so "init_gradu_tensorALL" needs to allocate HOLD_VELOCITY_DATA_MF
// internally for its' own uses then delete it after this call.
 if (localMF_grow[HOLD_VELOCITY_DATA_MF]!=-1)
  amrex::Error("localMF_grow[HOLD_VELOCITY_DATA_MF] invalid");

 int simple_AMR_BC_flag_viscosity=1;
 int do_alloc=1; 
 init_gradu_tensorALL(HOLD_VELOCITY_DATA_MF,do_alloc,
   CELLTENSOR_MF,FACETENSOR_MF,
   simple_AMR_BC_flag_viscosity);

 if (localMF_grow[HOLD_VELOCITY_DATA_MF]!=-1)
  amrex::Error("localMF_grow[HOLD_VELOCITY_DATA_MF] invalid");

 if (localMF[CELLTENSOR_MF]->nComp()!=ntensorMM)
  amrex::Error("localMF[CELLTENSOR_MF]->nComp() invalid");
 if (localMF[FACETENSOR_MF]->nComp()!=ntensorMM)
  amrex::Error("localMF[FACETENSOR_MF]->nComp() invalid");

 getStateVISC_ALL(CELL_VISC_MATERIAL_MF,1);
 if (localMF[CELL_VISC_MATERIAL_MF]->nComp()<nmat)
  amrex::Error("viscmf invalid ncomp");

 getStateDIV_ALL(MACDIV_MF,1);
 if (localMF[MACDIV_MF]->nComp()!=num_materials_vel)
  amrex::Error("localMF[MACDIV_MF]->nComp() invalid");

   // if FENE-CR+Carreau,
   // liquid viscosity=etaS+etaP ( 1+ (beta gamma_dot)^alpha )^((n-1)/alpha)
   //
   // for each material, there are 5 components:
   // 1. \dot{gamma}
   // 2. Tr(A) if viscoelastic
   //    \dot{gamma} o.t.
   // 3. Tr(A) (liquid viscosity - etaS)/etaP  if FENE-CR+Carreau
   //    Tr(A) if FENE-CR
   //    \dot{gamma} o.t.
   // 4. (3) * f(A)  if viscoelastic
   //    \dot{gamma} o.t.
   // 5. vorticity magnitude.

   // calls GETSHEAR and DERMAGTRACE
 int ntrace=5*nmat;
 getState_tracemag_ALL(MAGTRACE_MF,1);
 if (localMF[MAGTRACE_MF]->nComp()!=ntrace)
  amrex::Error("localMF[MAGTRACE_MF]->nComp() invalid");

  // idx,ngrow,scomp,ncomp,index,scompBC_map
 Vector<int> scompBC_map;
 scompBC_map.resize(num_materials_vel);
 for (int i=0;i<num_materials_vel;i++)
  scompBC_map[i]=0;
 PCINTERP_fill_bordersALL(MACDIV_MF,1,0,
   num_materials_vel,State_Type,scompBC_map);

 for (int i=0;i<ntrace;i++) {
  scompBC_map.resize(1);
  scompBC_map[0]=0;
  PCINTERP_fill_bordersALL(MAGTRACE_MF,1,i,1,State_Type,scompBC_map);
 }
 
 int output_MAC_vel=0;
 if (face_flag==1)
  output_MAC_vel=1;

 if (output_MAC_vel==1) {
  getStateALL(1,cur_time_slab,0,
    num_materials_vel*AMREX_SPACEDIM,HOLD_VELOCITY_DATA_MF);
  for (int ilev=finest_level;ilev>=0;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   int use_VOF_weight=1;
   ns_level.VELMAC_TO_CELL(use_VOF_weight);
  }
 }

 int tecplot_finest_level=finest_level;
 if (tecplot_max_level<tecplot_finest_level)
  tecplot_finest_level=tecplot_max_level;

 IntVect visual_fab_lo(IntVect::TheZeroVector()); 
 IntVect visual_fab_hi(visual_ncell); 
 Box visual_node_box(visual_fab_lo,visual_fab_hi);
 visual_fab_hi-=IntVect::TheUnitVector();
 Box visual_domain(visual_fab_lo,visual_fab_hi);
  // x,u,p,den,T,Y1..Yn,mag vort,LS
 int visual_ncomp=2*AMREX_SPACEDIM+3+num_species_var+1+nmat;  
 FArrayBox visual_fab_output(visual_node_box,visual_ncomp);
 FArrayBox visual_fab_input(visual_node_box,visual_ncomp); 

 visual_fab_output.setVal(-1.0e+20);
 visual_fab_input.setVal(-1.0e+20);

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
   int do_input=1;
   FORT_IO_COMPARE(
    &nmat,
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
  }
  ParallelDescriptor::Barrier();

   // communicate visual_fab_input from the IO proc to all of the
   // other processors.
  for (int i=gridlo[0];i<=gridhi[0];i++) {
   for (int j=gridlo[1];j<=gridhi[1];j++) {
    for (int k=gridlo[2];k<=gridhi[2];k++) {
     for (int n=0;n<visual_ncomp;n++) {
      Vector<int> arr_index(AMREX_SPACEDIM);
      arr_index[0]=i;
      arr_index[1]=j;
      if (AMREX_SPACEDIM==3) {
       arr_index[AMREX_SPACEDIM-1]=k;
      }
      IntVect p(arr_index);
      Real local_data=visual_fab_input(p,n); 
      ParallelDescriptor::ReduceRealMax(local_data);
      ParallelDescriptor::Barrier();
      if (1==0) {
       Box pbox(p,p);
       visual_fab_input.setVal(local_data,pbox,n);
      } else {
       visual_fab_input(p,n)=local_data;
      }
      ParallelDescriptor::Barrier();
     } // n
    } // k
   } // j
  } // i 

 } else if (visual_compare==0) {
  // do nothing
 } else
  amrex::Error("visual_compare invalid");

 for (int ilev=tecplot_finest_level;ilev>=0;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);

  ns_level.debug_ngrow(MACDIV_MF,1,250);

  MultiFab* velmf=ns_level.getState(1,0,
    num_materials_vel*(AMREX_SPACEDIM+1),cur_time_slab);
  MultiFab* presmf=ns_level.derive_EOS_pressure(); 
  if (presmf->nComp()!=num_materials_vel)
   amrex::Error("presmf has invalid ncomp");

  MultiFab* denmf=ns_level.getStateDen(1,cur_time_slab); 
  MultiFab* lsdist=ns_level.getStateDist(1,cur_time_slab,12);
  MultiFab* div_data=ns_level.getStateDIV_DATA(1,0,num_materials_vel,
    cur_time_slab);
  if (1==0) {
   std::cout << "level= " << ilev << " div_data norm0= " << 
    div_data->norm0() << '\n'; 
   std::cout << "level= " << ilev << " div_datanorm0+1grow= " << 
    div_data->norm0(0,1) << '\n'; 
  }
  MultiFab* viscoelasticmf=ns_level.getStateTensor(1,0,
     num_materials_viscoelastic*NUM_TENSOR_TYPE+AMREX_SPACEDIM,cur_time_slab);

  ns_level.output_zones(
   visual_fab_output,
   visual_domain,
   visual_ncomp,
   velmf,
   presmf,
   ns_level.localMF[MACDIV_MF],
   div_data,
   denmf,
   viscoelasticmf,
   lsdist,
   ns_level.localMF[CELL_VISC_MATERIAL_MF],
   ns_level.localMF[MAGTRACE_MF],
   grids_per_level_array[ilev],
   cgrids_minusBA_array[ilev],
   slice_data.dataPtr(), 
   do_plot,do_slice);

  if ((slice_dir>=0)&&(slice_dir<AMREX_SPACEDIM)) {
   for (int i=0;i<nslice*nstate_slice;i++)
    ParallelDescriptor::ReduceRealMax(slice_data[i]);
  } else
   amrex::Error("slice_dir invalid");

  delete viscoelasticmf;

  delete div_data;
  delete velmf;
  delete denmf;
  delete presmf;
  delete lsdist;
 }  // ilev=tecplot_finest_level ... 0

 ParallelDescriptor::Barrier();

 for (int i=gridlo[0];i<=gridhi[0];i++) {
  for (int j=gridlo[1];j<=gridhi[1];j++) {
   for (int k=gridlo[2];k<=gridhi[2];k++) {
    for (int n=0;n<visual_ncomp;n++) {
     Vector<int> arr_index(AMREX_SPACEDIM);
     arr_index[0]=i;
     arr_index[1]=j;
     if (AMREX_SPACEDIM==3) {
      arr_index[AMREX_SPACEDIM-1]=k;
     }
     IntVect p(arr_index);
     Real local_data=visual_fab_output(p,n); 
     ParallelDescriptor::ReduceRealMax(local_data);
     ParallelDescriptor::Barrier();
     if (1==0) {
      Box pbox(p,p);
      visual_fab_output.setVal(local_data,pbox,n);
     } else {
      visual_fab_output(p,n)=local_data;
     }
     ParallelDescriptor::Barrier();
    } // n
   } // k
  } // j
 } // i 

 if (ParallelDescriptor::IOProcessor()) {

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
  if ((nparts<0)||(nparts>nmat))
   amrex::Error("nparts invalid");
  Vector<int> im_solid_map_null;
  im_solid_map_null.resize(1);
  im_solid_map_null[0]=0;

  int* im_solid_map_ptr;
  int nparts_def=nparts;
  if (nparts==0) {
   im_solid_map_ptr=im_solid_map_null.dataPtr();
   nparts_def=1;
  } else if ((nparts>=1)&&(nparts<=nmat)) {
   im_solid_map_ptr=im_solid_map.dataPtr();
  } else
   amrex::Error("nparts invalid");

  if (do_plot==1) {
   FORT_COMBINEZONES(
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
    &visual_option,
    &visual_revolve,
    &plotint,
    &nmat, 
    &nparts,
    &nparts_def,
    im_solid_map_ptr);
  } else if (do_plot==0) {
   // do nothing
  } else
   amrex::Error("do_plot invalid");

  if (do_slice==1) {
   if ((slice_dir>=0)&&(slice_dir<AMREX_SPACEDIM)) {
    FORT_OUTPUTSLICE(&cur_time_slab,&nsteps,&sliceint,
     slice_data.dataPtr(),&nslice,&nstate_slice,
     &visual_option);
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
  FORT_IO_COMPARE(
   &nmat,
   &nsteps,
   &do_input,
   &visual_compare,  // unused since do_input==0
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

 } else if (!ParallelDescriptor::IOProcessor()) {
  // do nothing
 } else
  amrex::Error("ParallelDescriptor::IOProcessor() corrupt");

 ParallelDescriptor::Barrier();

 if (output_MAC_vel==1) {
  for (int ilev=finest_level;ilev>=0;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   MultiFab& S_new=ns_level.get_new_data(State_Type,slab_step+1);
   MultiFab::Copy(S_new,*ns_level.localMF[HOLD_VELOCITY_DATA_MF],
     0,0,num_materials_vel*AMREX_SPACEDIM,1);
   ns_level.delete_localMF(HOLD_VELOCITY_DATA_MF,1);
  }  // ilev
 }

 delete_array(MACDIV_MF);
 delete_array(MAGTRACE_MF); 
 delete_array(CELLTENSOR_MF);
 delete_array(FACETENSOR_MF);

} // subroutine writeTECPLOT_File



void NavierStokes::writeSanityCheckData(
		const std::string& root_string,
		const std::string& caller_string,
		int data_id,
                int ncomp,
                int data_mf, 
		int state_type_mf,
                int data_dir) {

 if (ParallelDescriptor::IOProcessor()) {
  std::cout << "in: writeSanityCheckData, root_string= " <<
    root_string << '\n';
  std::cout << "in: writeSanityCheckData, caller_string= " <<
    caller_string << '\n';
  std::cout << "in: writeSanityCheckData, data_mf= " <<
    data_mf << " state_type_mf=" << state_type_mf << '\n';
 }

 if (level!=0)
  amrex::Error("level invalid writeSanityCheckData");

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid");

 int nsteps=parent->levelSteps(0);

 int finest_level = parent->finestLevel();

 ParallelDescriptor::Barrier();

 Vector<int> grids_per_level_array;
 grids_per_level_array.resize(finest_level+1);
 Vector<BoxArray> cgrids_minusBA_array;
 cgrids_minusBA_array.resize(finest_level+1);

 MultiFab* raw_data_lev0_mf;

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
  } else
   amrex::Error("state_type_mf invalid");

  if (raw_data_lev0_mf->nGrow()>=0) {
   // do nothing
  } else
   amrex::Error("raw_data_lev0_mf->nGrow() invalid");

 } else
  amrex::Error("data_mf invalid");

 int tecplot_finest_level=finest_level;
 if (tecplot_max_level<tecplot_finest_level)
  tecplot_finest_level=tecplot_max_level;

 for (int ilev=tecplot_finest_level;ilev>=0;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);

  MultiFab* raw_data_mf;

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
   ns_level.debug_ngrow(data_mf,0,250);

  } else if (data_mf==-1) {

   if (state_type_mf>=0) {
    raw_data_mf=&ns_level.get_new_data(state_type_mf,slab_step+1);
   } else
    amrex::Error("state_type_mf invalid");

   if (raw_data_mf->nGrow()>=0) {
    // do nothing
   } else
    amrex::Error("raw_data_mf->nGrow() invalid");

  } else
   amrex::Error("data_mf invalid");

   // data_dir=-1 cell centered data
   // data_dir=0..sdim-1 face centered data.
   // data_dir=sdim node data
  ns_level.Sanity_output_zones(
   data_id,
   data_dir,
   raw_data_mf,
   ncomp,
   grids_per_level_array[ilev],
   cgrids_minusBA_array[ilev]);

 }  // ilev=tecplot_finest_level ... 0

 ParallelDescriptor::Barrier();

 if (ParallelDescriptor::IOProcessor()) {

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
  char root_char_array[n_root+1];
  strcpy(root_char_array,root_string.c_str());

  if (temp_number_grids!=total_number_grids)
   amrex::Error("temp_number_grids invalid");
  
  int num_levels=tecplot_finest_level+1;

   // in: NAVIERSTOKES_3D.F90
  FORT_COMBINEZONES_SANITY(
    root_char_array,
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
    &data_id,
    &nsteps,
    &num_levels,
    &cur_time_slab,
    &visual_option,
    &visual_revolve,
    &ncomp);

 } else if (!ParallelDescriptor::IOProcessor()) {
  // do nothing
 } else
  amrex::Error("ParallelDescriptor::IOProcessor() corrupt");

 ParallelDescriptor::Barrier();

} // subroutine writeSanityCheckData


void
NavierStokes::writePlotFile (
  const std::string& dir,
  std::ostream& os,
  int do_plot,int do_slice,
  int SDC_outer_sweeps_in,int slab_step_in) {

 SDC_setup();
 ns_time_order=parent->Time_blockingFactor();

 SDC_outer_sweeps=SDC_outer_sweeps_in;
 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid");

 slab_step=slab_step_in; 

 SDC_setup_step();

 int i, n;

 int f_lev = parent->finestLevel();

  // metrics_dataALL
  // MASKCOEF_MF
  // MASK_NBR_MF
  // init_FSI_GHOST_MAC_MF_ALL
  // LEVELPC_MF
  // MASKSEM_MF
  // VOF_Recon_ALL
  // make_physics_varsALL
 if (level==0) {
  int post_init_flag=0; // in: writePlotFile
  prepare_post_process(post_init_flag);
 } // level==0

  // output tecplot zonal files  x,y,z,u,v,w,phi,psi
 if (visual_option==-2) {

  if (level==0) {
   writeTECPLOT_File(do_plot,do_slice);
  }

  // output isosurface for LevelVapor and LevelSolid
  // output plot files
 } else if (visual_option==-1) {

  if (level==0) { 
   if (do_slice==1) {
    int do_plot_kluge=0;
    writeTECPLOT_File(do_plot_kluge,do_slice);
   } 
  }

  ParallelDescriptor::Barrier();

  if (do_plot==1) {

   if (level == 0) {
    writeInterfaceReconstruction();
   } // level==0

   //
   // The list of indices of State to write to plotfile.
   // first component of pair is state_type,
   // second component of pair is component # within the state_type
   //

   std::vector<std::pair<int,int> > plot_var_map;

   for (int typ=State_Type;typ<NUM_STATE_TYPE;typ++) {
    for (int comp = 0; comp < desc_lst[typ].nComp();comp++) {
     if ((parent->isStatePlotVar(desc_lst[typ].name(comp))) &&
         (desc_lst[typ].getType() == IndexType::TheCellType())) {
      plot_var_map.push_back(std::pair<int,int>(typ,comp));
     }
    }  // comp
   } // typ

   int n_data_items = plot_var_map.size();

   if (level == 0 && ParallelDescriptor::IOProcessor()) {
     //
     // The first thing we write out is the plotfile type.
     //
     os << thePlotFileType() << '\n';

     if (n_data_items == 0)
         amrex::Error("Must specify at least one valid data item to plot");

     os << n_data_items << '\n';

     //
     // Names of variables -- first state, then derived
     //
     for (i =0; i < plot_var_map.size(); i++) {
         int typ  = plot_var_map[i].first;
         int comp = plot_var_map[i].second;
         os << desc_lst[typ].name(comp) << '\n';
     }

     os << AMREX_SPACEDIM << '\n';
     os << parent->cumTime() << '\n';
     os << f_lev << '\n';
     for (i = 0; i < AMREX_SPACEDIM; i++)
         os << geom.ProbLo(i) << ' ';
     os << '\n';
     for (i = 0; i < AMREX_SPACEDIM; i++)
         os << geom.ProbHi(i) << ' ';
     os << '\n';
     for (i = 0; i < f_lev; i++)
         os << 2 << ' ';
     os << '\n';
     for (i = 0; i <= f_lev; i++)
         os << parent->Geom(i).Domain() << ' ';
     os << '\n';
     for (i = 0; i <= f_lev; i++)
         os << parent->levelSteps(i) << ' ';
     os << '\n';
     for (i = 0; i <= f_lev; i++) {
         for (int k = 0; k < AMREX_SPACEDIM; k++)
             os << parent->Geom(i).CellSize()[k] << ' ';
         os << '\n';
     }
     os << (int) geom.Coord() << '\n';
     os << "0\n"; // Write bndry data.
   }
   // Build the directory to hold the MultiFab at this level.
   // The name is relative to the directory containing the Header file.
   //
   static const std::string BaseName = "/Cell";
   char buf[64];
   sprintf(buf, "Level_%d", level);
   std::string Level_str = buf;
   //
   // Now for the full pathname of that directory.
   //
   std::string FullPath = dir;
   if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/')
     FullPath += '/';
   FullPath += Level_str;
   //
   // Only the I/O processor makes the directory if it doesn't already exist.
   //
   if (ParallelDescriptor::IOProcessor())
     if (!amrex::UtilCreateDirectory(FullPath, 0755))
         amrex::CreateDirectoryFailed(FullPath);
   //
   // Force other processors to wait till directory is built.
   //
   ParallelDescriptor::Barrier();

   if (ParallelDescriptor::IOProcessor()) {
     os << level << ' ' << grids.size() << ' ' << cur_time_slab << '\n';
     os << parent->levelSteps(level) << '\n';

     for (i = 0; i < grids.size(); ++i) {
         for (n = 0; n < AMREX_SPACEDIM; n++)
             os << grid_loc[i].lo(n) << ' ' << grid_loc[i].hi(n) << '\n';
     }
     //
     // The full relative pathname of the MultiFabs at this level.
     // The name is relative to the Header file containing this name.
     // It's the name that gets written into the Header.
     //
     if (n_data_items > 0) {
         std::string PathNameInHeader = Level_str;
         PathNameInHeader += BaseName;
         os << PathNameInHeader << '\n';
     }
   } // if IOProc

   //
   // We combine all of the multifabs -- state, etc -- into one
   // multifab -- plotMF.
   // Each state variable has one component.
   // The VOF variables have to be obtained together with Centroid vars.
   int       cnt   = 0;
   int       ncomp = 1;
   const int nGrow = 0;
   MultiFab  plotMF(grids,dmap,n_data_items,nGrow,
		 MFInfo().SetTag("plotMF"),FArrayBoxFactory());
   //
   // Cull data from state variables 
   //
   for (i = 0; i < plot_var_map.size(); i++) {
     int typ  = plot_var_map[i].first;
     int comp = plot_var_map[i].second;
     MultiFab* this_dat;
     ncomp=1;

     if (typ==State_Type) {
      int nmat=num_materials;
      int scomp_mofvars=num_materials_vel*(AMREX_SPACEDIM+1)+
       nmat*num_state_material;
      if (comp==scomp_mofvars) {
       ncomp=nmat*ngeom_raw;
      }
      this_dat=getState(nGrow,comp,ncomp,cur_time_slab);
     } else if ((typ==Solid_State_Type)&&(im_solid_map.size()!=0)) {

      if (comp!=0) {
       std::cout << "comp=" << comp << " ncomp= " <<
        desc_lst[typ].nComp() << '\n';
       amrex::Error("comp invalid for Solid_State_Type");
      }
      int nmat=num_materials;
      int nparts=im_solid_map.size();
      if ((nparts<1)||(nparts>nmat))
       amrex::Error("nparts invalid");
      if (comp==0) {
       ncomp=nparts*AMREX_SPACEDIM;
      } else
       amrex::Error("comp invalid");

      this_dat=getStateSolid(nGrow,comp,ncomp,cur_time_slab);

      if (this_dat->nComp()!=ncomp)
       amrex::Error("this_dat->nComp() invalid");

     } else if ((typ==Tensor_Type)&&
		(im_elastic_map.size()>=0)) {

      if (comp!=0) {
       std::cout << "comp=" << comp << " ncomp= " <<
        desc_lst[typ].nComp() << '\n';
       amrex::Error("comp invalid for Tensor_Type");
      }
      int nmat=num_materials;
      int nparts=im_elastic_map.size();
      if ((nparts<0)||(nparts>nmat)||
          (nparts!=num_materials_viscoelastic))
       amrex::Error("nparts invalid");
      if (comp==0) {
       ncomp=nparts*NUM_TENSOR_TYPE+AMREX_SPACEDIM;
      } else
       amrex::Error("comp invalid");

      this_dat=getStateTensor(nGrow,comp,ncomp,cur_time_slab);

      if (this_dat->nComp()!=ncomp) 
       amrex::Error("this_dat->nComp() invalid");

     } else if (typ==LS_Type) {
      if (comp!=0) {
       std::cout << "comp=" << comp << " ncomp= " <<
        desc_lst[typ].nComp() << '\n';
       amrex::Error("comp invalid for LS_Type");
      }

      int nmat=num_materials;
      if (comp==0) {
       ncomp=nmat*(AMREX_SPACEDIM+1);
      } else
       amrex::Error("comp invalid");

      this_dat=getStateDist(nGrow,cur_time_slab,13);
      if (this_dat->nComp()!=ncomp)
       amrex::Error("this_dat->nComp() invalid");

     } else if (typ==DIV_Type) {
      if (comp!=0)
       amrex::Error("comp invalid for DIV_Type");
      ncomp=num_materials_vel;
      this_dat=getStateDIV_DATA(1,comp,ncomp,cur_time_slab);
     } else
      amrex::Error("typ invalid");

     MultiFab::Copy(plotMF,*this_dat,0,cnt,ncomp,nGrow);
     delete this_dat;
     cnt+= ncomp;
     i+=(ncomp-1);
   }  // i
   //
   // Use the Full pathname when naming the MultiFab.
   //
   std::string TheFullPath = FullPath;
   TheFullPath += BaseName;
   VisMF::Write(plotMF,TheFullPath);
   ParallelDescriptor::Barrier();
  } else if (do_plot==0) {
   // do nothing
  } else
   amrex::Error("do_plot invalid");

 } else
  amrex::Error("visual_option invalid try -1 or -2");

}

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
 ns_reconcile_d_num(105);
}

// called from: estTimeStep
void NavierStokes::MaxAdvectSpeedALL(Real& dt_min,
  Real* vel_max,Real* vel_max_estdt,Real& vel_max_cap_wave) {

 int finest_level=parent->finestLevel();
 if (level!=0)
  amrex::Error("level invalid MaxAdvectSpeedALL");

 last_finest_level=finest_level; 
 
 Real local_vel_max[AMREX_SPACEDIM+1];  // last component is max|c|^2
 Real local_vel_max_estdt[AMREX_SPACEDIM+1];  // last component is max|c|^2
 Real local_vel_max_cap_wave;
 Real local_dt_min;

 for (int dir=0;dir<AMREX_SPACEDIM+1;dir++) {
  vel_max[dir]=0.0;
  vel_max_estdt[dir]=0.0;
 }
 vel_max_cap_wave=0.0;

 dt_min=1.0E+30;
 if (dt_min<dt_max) 
  dt_min=dt_max;

 if (localMF_grow[FSI_GHOST_MAC_MF]<0) {

  init_FSI_GHOST_MAC_MF_ALL(1);

 }

 for (int ilev=finest_level;ilev>=0;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  for (int dir=0;dir<AMREX_SPACEDIM+1;dir++) {
   local_vel_max[dir]=0.0;
   local_vel_max_estdt[dir]=0.0;
  }
  local_vel_max_cap_wave=0.0;
  local_dt_min=1.0E+25;
  if (local_dt_min<dt_max) 
   local_dt_min=dt_max;

  ns_level.MaxAdvectSpeed(local_dt_min,local_vel_max,local_vel_max_estdt,
		  local_vel_max_cap_wave); 
  for (int dir=0;dir<AMREX_SPACEDIM+1;dir++) {
   vel_max[dir] = std::max(vel_max[dir],local_vel_max[dir]);
   vel_max_estdt[dir] = std::max(vel_max_estdt[dir],local_vel_max_estdt[dir]);
  }
  vel_max_cap_wave = std::max(vel_max_cap_wave,local_vel_max_cap_wave);
  dt_min=std::min(dt_min,local_dt_min);
 } // ilev 
} // end subroutine MaxAdvectSpeedALL

// vel_max[0,1,2]= max vel in direction  vel_max[sdim]=max c^2
void NavierStokes::MaxAdvectSpeed(Real& dt_min,Real* vel_max,
 Real* vel_max_estdt,Real& vel_max_cap_wave) {

 int finest_level=parent->finestLevel();
 int nmat=num_materials;
 int scomp_mofvars=num_materials_vel*(AMREX_SPACEDIM+1)+
  nmat*num_state_material;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;
 int nsolve=1;
 int nsolveMM_FACE=nsolve*num_materials_vel;

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 int nparts=im_solid_map.size();
 if ((nparts<0)||(nparts>nmat))
  amrex::Error("nparts invalid");
 Vector<int> im_solid_map_null;
 im_solid_map_null.resize(1);

 int* im_solid_map_ptr;
 int nparts_def=nparts;
 if (nparts==0) {
  im_solid_map_ptr=im_solid_map_null.dataPtr();
  nparts_def=1;
 } else if ((nparts>=1)&&(nparts<=nmat)) {
  im_solid_map_ptr=im_solid_map.dataPtr();
 } else
  amrex::Error("nparts invalid");

 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
  if (localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0)
   amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0");
  if (localMF[FSI_GHOST_MAC_MF+data_dir]->nComp()!=nparts_def*AMREX_SPACEDIM)
   amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nComp() invalid");
 }

 MultiFab* distmf=getStateDist(2,cur_time_slab,14);
 MultiFab* denmf=getStateDen(1,cur_time_slab);  // nmat*num_state_material
 MultiFab* vofmf=getState(1,scomp_mofvars,nmat*ngeom_raw,cur_time_slab);

 int rzflag=0;
 if (geom.IsRZ())
  rzflag=1;
 else if (geom.IsCartesian())
  rzflag=0;
 else if (geom.IsCYLINDRICAL())
  rzflag=3;
 else
  amrex::Error("CoordSys bust 5");

 const Real* dx = geom.CellSize();

 dt_min=1.0E+25;
 if (dt_min<dt_max) 
  dt_min=dt_max;

 for (int dir=0;dir<AMREX_SPACEDIM+1;dir++) {
  vel_max[dir]=0.0;
  vel_max_estdt[dir]=0.0;
 }
 vel_max_cap_wave=0.0;

 Vector< Vector<Real> > local_cap_wave_speed;
 Vector< Vector<Real> > local_vel_max;
 Vector< Vector<Real> > local_vel_max_estdt;
 Vector<Real> local_vel_max_cap_wave;
 Vector< Real > local_dt_min;
 local_cap_wave_speed.resize(thread_class::nthreads);
 local_vel_max.resize(thread_class::nthreads);
 local_vel_max_estdt.resize(thread_class::nthreads);
 local_vel_max_cap_wave.resize(thread_class::nthreads);
 local_dt_min.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {

  local_cap_wave_speed[tid].resize(nten); 
  for (int iten=0;iten<nten;iten++) {
   local_cap_wave_speed[tid][iten]=cap_wave_speed[iten];
  }

  local_vel_max[tid].resize(AMREX_SPACEDIM+1);//last component max|c|^2
  local_vel_max_estdt[tid].resize(AMREX_SPACEDIM+1);//last component max|c|^2
  for (int dir=0;dir<AMREX_SPACEDIM+1;dir++) {
   local_vel_max[tid][dir]=0.0;
   local_vel_max_estdt[tid][dir]=0.0;
  }
  local_vel_max_cap_wave[tid]=0.0;

  local_dt_min[tid]=1.0E+25;
  if (local_dt_min[tid]<dt_max) 
   local_dt_min[tid]=dt_max;
 }  // tid 

 MultiFab* velcell=getState(1,0,num_materials_vel*AMREX_SPACEDIM,cur_time_slab);
  
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
 
  MultiFab* velmac=getStateMAC(0,dir,0,nsolveMM_FACE,cur_time_slab);

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
   FArrayBox& solidfab=(*localMF[FSI_GHOST_MAC_MF+dir])[mfi];

   int local_enable_spectral=enable_spectral;

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // in: GODUNOV_3D.F90
   FORT_ESTDT(
    &nsolveMM_FACE,
    &local_enable_spectral,
    elastic_time.dataPtr(),
    microlayer_substrate.dataPtr(),
    microlayer_angle.dataPtr(),
    microlayer_size.dataPtr(),
    macrolayer_size.dataPtr(),
    latent_heat.dataPtr(),
    reaction_rate.dataPtr(),
    freezing_model.dataPtr(),
    Tanasawa_or_Schrage.dataPtr(),
    distribute_from_target.dataPtr(),
    saturation_temp.dataPtr(),
    mass_fraction_id.dataPtr(),
    molar_mass.dataPtr(),
    species_molar_mass.dataPtr(),
    species_evaporation_density.dataPtr(),
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
    local_cap_wave_speed[tid_current].dataPtr(),
    local_vel_max[tid_current].dataPtr(),
    local_vel_max_estdt[tid_current].dataPtr(),
    &local_vel_max_cap_wave[tid_current],
    &local_dt_min[tid_current],
    &rzflag,
    &Uref,&Lref,
    &nten,
    &use_lsa,
    denconst.dataPtr(),
    denconst_gravity.dataPtr(),
    &visc_coef,
    &gravity,
    &terminal_velocity_dt,
    &dir,
    &nmat,
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
  }  //mfi
} // omp
  ns_reconcile_d_num(106);

  for (int tid=1;tid<thread_class::nthreads;tid++) {

   for (int iten=0;iten<nten;iten++)
    if (local_cap_wave_speed[tid][iten]>local_cap_wave_speed[0][iten])
     local_cap_wave_speed[0][iten]=local_cap_wave_speed[tid][iten];

   if (local_dt_min[tid]<local_dt_min[0])
    local_dt_min[0]=local_dt_min[tid];

   if (local_vel_max[tid][dir]>local_vel_max[0][dir])
    local_vel_max[0][dir]=local_vel_max[tid][dir];
   if (local_vel_max[tid][AMREX_SPACEDIM]>local_vel_max[0][AMREX_SPACEDIM])
    local_vel_max[0][AMREX_SPACEDIM]=local_vel_max[tid][AMREX_SPACEDIM];

   if (local_vel_max_estdt[tid][dir]>local_vel_max_estdt[0][dir])
    local_vel_max_estdt[0][dir]=local_vel_max_estdt[tid][dir];
   if (local_vel_max_estdt[tid][AMREX_SPACEDIM]>
       local_vel_max_estdt[0][AMREX_SPACEDIM])
    local_vel_max_estdt[0][AMREX_SPACEDIM]=
	    local_vel_max_estdt[tid][AMREX_SPACEDIM];

   if (local_vel_max_cap_wave[tid]>local_vel_max_cap_wave[0])
    local_vel_max_cap_wave[0]=local_vel_max_cap_wave[tid];
  } // tid

  for (int iten=0;iten<nten;iten++) {
   ParallelDescriptor::ReduceRealMax(local_cap_wave_speed[0][iten]);
  }

  ParallelDescriptor::ReduceRealMax(local_vel_max[0][dir]);
  ParallelDescriptor::ReduceRealMax(local_vel_max[0][AMREX_SPACEDIM]);
  ParallelDescriptor::ReduceRealMax(local_vel_max_estdt[0][dir]);
  ParallelDescriptor::ReduceRealMax(local_vel_max_estdt[0][AMREX_SPACEDIM]);

  ParallelDescriptor::ReduceRealMax(local_vel_max_cap_wave[0]);

  ParallelDescriptor::ReduceRealMin(local_dt_min[0]);

  delete velmac;
 }  // dir=0..sdim-1

 delete velcell;

 for (int iten=0;iten<nten;iten++)
  cap_wave_speed[iten]=local_cap_wave_speed[0][iten];

 for (int dir=0;dir<=AMREX_SPACEDIM;dir++) {
  vel_max[dir]=local_vel_max[0][dir];
  vel_max_estdt[dir]=local_vel_max_estdt[0][dir];
 }
 vel_max_cap_wave=local_vel_max_cap_wave[0];

 dt_min=local_dt_min[0];

 delete denmf;
 delete vofmf;
 delete distmf;

} // subroutine MaxAdvectSpeed

// called from: computeNewDt, computeInitialDt
Real NavierStokes::estTimeStep (Real local_fixed_dt) {

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

 if (local_fixed_dt>0.0) {

  return_dt=local_fixed_dt;

 } else if (local_fixed_dt==0.0) {

  if (fixed_dt_velocity > 0.0) {

   return_dt=smallest_dx/fixed_dt_velocity;

  } else if (fixed_dt_velocity==0.0) {

   Real dt_min;
   Real u_max[AMREX_SPACEDIM+1];  // last component is max|c|^2
   Real u_max_estdt[AMREX_SPACEDIM+1];  // last component is max|c|^2
   Real u_max_cap_wave;

   MaxAdvectSpeedALL(dt_min,u_max,u_max_estdt,u_max_cap_wave);
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
  FORT_GL_SLAB(time_array.dataPtr(),
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

void NavierStokes::post_regrid (int lbase,
  int start_level,int new_finest,int initialInit_flag,Real time) {

  NavierStokes& ns_level0=getLevel(0);

  const int max_level = parent->maxLevel();
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

  for (int ipart=0;ipart<NS_ncomp_particles;ipart++) {
   AmrParticleContainer<N_EXTRA_REAL,0,0,0>& current_PC=
      ns_level0.get_new_dataPC(State_Type,ns_time_order,ipart);
   int lev_min=0;
   int lev_max=-1;
   int nGrow_Redistribute=0;
   int local_Redistribute=0;
   current_PC.Redistribute(lev_min,lev_max,nGrow_Redistribute, 
     local_Redistribute);
  } // ipart=0..NS_ncomp_particles-1
    
   // olddata=newdata  
  for (int k=0;k<nstate;k++) {
   state[k].CopyNewToOld(level,max_level); 
   state[k].setTimeLevel(time,dt_amr);
  }
}

void NavierStokes::computeNewDt (int finest_level,
  Real& dt,Real stop_time,int post_regrid_flag) {

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

  if (level==0) {

   Real newdt=estTimeStep(local_fixed_dt);

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

}

void NavierStokes::computeInitialDt (int finest_level,
   Real& dt,Real stop_time) {

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "start: computeInitialDt \n";
  }
 }

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

 Real newdt=init_shrink*estTimeStep(fixed_dt_init)/shrink_factor;

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

 SDC_outer_sweeps=0;
 slab_step=ns_time_order-1;
 SDC_setup_step();

 if (level==0) {
  if (sum_interval>0) {
   if ( (parent->levelSteps(0)%sum_interval == 0)||
        (stop_time-upper_slab_time<1.0E-8) ) {
    int post_init_flag=0;
    sum_integrated_quantities(post_init_flag);
   }
  }
 } 

 init_regrid_history();
}

//
// Ensure state, and pressure are consistent.
//
void
NavierStokes::post_init (Real stop_time)
{

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
   if (ns_level.localMF_grow[FSI_MF]>=0)
    ns_level.delete_localMF(FSI_MF,1);
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
  // initialize particles.
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
     int post_init_flag=1;
     sum_integrated_quantities(post_init_flag);
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

}  // subroutine post_init

// HH is m+1 x m
// yy is m x 1
// solve H^T H y = H^T beta e1  (slow, since H is m+1xm  H^T H is mxm
//   Gaussian elimination for mxm system has a cost of O(m^3))
//  Givens' rotation approach has a cost of O(m^2)
//  ( min_y ||Hy - beta e1|| )
// status==1 success
// status==0 poor conditioning
// HT is m x m+1
// HH=  h11 h12 h13 ... h1m
//      h21 h22 h23 ... h2m
//          h32 h33 ... h3m
//              h43 ... h4m
//                  .
//                  .
//               hm,m-1 hmm
//                      hm+1,m
//
// HT= h11 h21
//     h12 h22 h32
//     h13 h23 h33 h43
//         .
//         .
//         .
//     h1m h2m h3m h4m ... hmm hm+1,m
//
// beta_vec_input[m+1]
// HH_input[m+1][m]
//
void GMRES_HELPER(Real** HH_input,Real* beta_vec_input,Real* yy,int m) {

 if (m>=1) {
  // do nothing
 } else
  amrex::Error("m invalid in GMRES_HELPER");

 Real** HH=new Real*[m+1];
 for (int i=0;i<m+1;i++) 
  HH[i]=new Real[m];

 Real* sn=new Real[m];
 Real* cs=new Real[m];
 Real* h=new Real[m+1];
 Real* beta_vec=new Real[m+1];
 Real beta1;
 Real beta2;

 for (int i=0;i<m+1;i++) {
  for (int j=0;j<m;j++) {
   if ((HH_input[i][j]>=0.0)||
       (HH_input[i][j]<=0.0)) {
    HH[i][j]=HH_input[i][j];
   } else
    amrex::Error("HH_input bust");
  }
  if ((beta_vec_input[i]>=0.0)||
      (beta_vec_input[i]<=0.0)) {
   beta_vec[i]=beta_vec_input[i];
  } else
   amrex::Error("beta_vec_input bust");
 } // i=0..m

 for (int i=0;i<m;i++) {
  sn[i]=0.0;
  cs[i]=0.0;
 }
 for (int i=0;i<m+1;i++) {
  h[i]=0.0;
 }

 for (int k=1;k<=m;k++) {

  for (int i=1;i<=m+1;i++)
   h[i-1]=0.0;
  for (int i=1;i<=k+1;i++)
   h[i-1]=HH[i-1][k-1];

  for (int i=1;i<=k-1;i++) {
   Real temp=cs[i-1]*h[i-1]+sn[i-1]*h[i];
   h[i]=-sn[i-1]*h[i-1]+cs[i-1]*h[i];
   h[i-1]=temp;
  }
  Real cs_k=0.0;
  Real sn_k=0.0;
  Real v1=h[k-1];
  Real v2=h[k];
  if (v1==0.0) {
   cs_k=0.0;
   sn_k=1.0;
  } else if (v1!=0.0) {
   Real t=sqrt(v1*v1+v2*v2);
   if (t>0.0) {
    cs_k=std::abs(v1)/t;
    sn_k=cs_k*v2/v1; 
   } else {
    std::cout << "t= " << t << '\n';
    amrex::Error("t invalid");
   }
  } else
   amrex::Error("v1 is corrupt");

    // the rotation matrix:  ( cs  sn 
    //                         -sn cs ) 
  h[k-1]=cs_k*h[k-1]+sn_k*h[k];
  h[k]=0.0;

  for (int i=1;i<=k+1;i++)
   HH[i-1][k-1]=h[i-1];

  cs[k-1]=cs_k;
  sn[k-1]=sn_k;

  beta1=beta_vec[k-1];
  beta2=beta_vec[k];
  beta_vec[k-1]=cs[k-1]*beta1+sn[k-1]*beta2;
  beta_vec[k]=-sn[k-1]*beta1+cs[k-1]*beta2;

 } // k=1..m

 for (int k=m;k>=1;k--) {

  yy[k-1]=beta_vec[k-1];
  for (int j=k+1;j<=m;j++)
   yy[k-1]-=HH[k-1][j-1]*yy[j-1];
  Real hdiag=HH[k-1][k-1];
  if (std::abs(hdiag)>0.0) 
   yy[k-1]/=hdiag;
  else
   amrex::Error("hdiag became 0");

 } // k=m ... 1 

 delete [] sn;
 delete [] cs;
 delete [] h;
 delete [] beta_vec;

 for (int i=0;i<m+1;i++) 
  delete [] HH[i];
 delete [] HH;

} // end subroutine GMRES_HELPER


void GMRES_MIN_CPP(Real** HH,Real beta, Real* yy,
		int m,int m_small,
		int caller_id,int project_option,
                int mg_level,int& status) {

#define profile_gmres 0

#if (profile_gmres==1)
 std::string subname="GMRES_MIN_CPP";
 std::stringstream id_string_stream(std::stringstream::in |
   std::stringstream::out);
 id_string_stream << caller_id;
 std::string profname=subname+id_string_stream.str();

 BLProfiler bprof(profname);
#endif

 status=1; // status==1 means success, status==0 means poor conditioning.

 if ((m_small>=1)&&(m_small<=m)) {
  // do nothing
 } else
  amrex::Error("m_small invalid");

 Real** HH_small=new Real*[m_small+1];
 for (int i=0;i<m_small+1;i++) { 
  HH_small[i]=new Real[m_small];
  for (int j=0;j<m_small;j++) 
   HH_small[i][j]=HH[i][j];
 }

 Real** HCOPY=new Real*[m_small+1];
 for (int i=0;i<m_small+1;i++) 
  HCOPY[i]=new Real[m_small];

 for (int i=0;i<m_small+1;i++) { 
  for (int j=0;j<m_small;j++) {

   if (HH[i][j]<0.0) {
    // do nothing
   } else if (HH[i][j]>0.0) {
    // do nothing
   } else if (HH[i][j]==0.0) {
    // do nothing
   } else {
    amrex::Error("HH[i][j] corrupt");
   }
   if (i>=j+2) {
    if (HH[i][j]==0.0) {
     // do nothing
    } else
     amrex::Error("HH should be 0");
   } else if (i<j+2) {
    // check nothing
   } else
    amrex::Error("i or j became corrupt");

   if (HH_small[i][j]<0.0) {
    // do nothing
   } else if (HH_small[i][j]>0.0) {
    // do nothing
   } else if (HH_small[i][j]==0.0) {
    // do nothing
   } else {
    amrex::Error("HH_small[i][j] corrupt");
   }
   if (i>=j+2) {
    if (HH_small[i][j]==0.0) {
     // do nothing
    } else
     amrex::Error("HH_small should be 0");
   } else if (i<j+2) {
    // check nothing
   } else
    amrex::Error("i or j became corrupt");

   HCOPY[i][j]=HH[i][j];
  } // j=0;j<m_small
 } //  i=0;i<m_small+1

 if (beta>0.0) {
  // do nothing
 } else {
  std::cout << "beta= " << beta << '\n';
  amrex::Error("beta invalid");
 }

 Real* delta_y=new Real[m_small];

 Real* beta_vec=new Real[m_small+1];
 for (int i=0;i<m_small+1;i++) {
  beta_vec[i]=0.0;
 }
 beta_vec[0]=beta;

 GMRES_HELPER(HCOPY,beta_vec,yy,m_small);

 for (int i=0;i<m_small+1;i++) {
  if (i==0) {
   beta_vec[i]=beta;
  } else if ((i>=1)&&(i<m_small+1)) {
   beta_vec[i]=0.0;
  } else
   amrex::Error("i invalid");

  for (int j=0;j<m_small;j++)
   beta_vec[i]-=HCOPY[i][j]*yy[j];
 } // i=0..m_small+1

 GMRES_HELPER(HCOPY,beta_vec,delta_y,m_small);

 Real norm_y=0.0;
 Real norm_delta_y=0.0;
 for (int i=0;i<m_small;i++) { 
  norm_y+=yy[i]*yy[i];
  norm_delta_y+=delta_y[i]*delta_y[i];
 }

 if (norm_delta_y>=0.0) {
  norm_delta_y=sqrt(norm_delta_y);
 } else
  amrex::Error("norm_delta_y invalid");

 if (norm_y>0.0) {
  norm_y=sqrt(norm_y);
  double relative_error=norm_delta_y/norm_y;

  if (relative_error>0.01) {
   status=0;
   if (1==0) {
    std::cout << "caller_id= " << caller_id << '\n';
    std::cout << "project_option= " << project_option << '\n';
    std::cout << "mg_level= " << mg_level << '\n';
    std::cout << "relative_error= " << relative_error << '\n';
    std::cout << "beta= " << beta << '\n';
    std::cout << "norm_y= " << norm_y << '\n';
    std::cout << "norm_delta_y= " << norm_delta_y << '\n';
    std::cout << "m_small= " << m_small << '\n';
    amrex::Error("relative_error large, decrease ns.mglib_min_coeff_factor");
   }
  }
 } else if (norm_y==0.0) {
  status=0;
  if ((norm_delta_y!=0.0)||(1==0)) {
   std::cout << "caller_id= " << caller_id << '\n';
   std::cout << "project_option= " << project_option << '\n';
   std::cout << "mg_level= " << mg_level << '\n';
   std::cout << "beta= " << beta << '\n';
   std::cout << "norm_y= " << norm_y << '\n';
   std::cout << "norm_delta_y= " << norm_delta_y << '\n';
   std::cout << "m_small= " << m_small << '\n';
   amrex::Error("norm_y became zero, and norm_delta_y might be invalid too");
  }

 } else
  amrex::Error("norm_y invalid");

 delete [] beta_vec;

 delete [] delta_y;

 for (int i=0;i<m_small+1;i++) 
  delete [] HCOPY[i];
 delete [] HCOPY;

 for (int i=0;i<m_small+1;i++) 
  delete [] HH_small[i];
 delete [] HH_small;

#if (profile_gmres==1)
 bprof.stop();
#endif 

#undef profile_gmres

} // subroutine GMRES_MIN_CPP

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

} // matrix_solveCPP


void
NavierStokes::volWgtSumALL(
 int post_init_flag,
 Vector<Real>& result,
 Vector<int>& sumdata_type,
 Vector<int>& sumdata_sweep,
 Vector<Real>& ZZ,Vector<Real>& FF,
 int dirx,int diry,int cut_flag,
 int isweep) {

 int finest_level=parent->finestLevel();
 if (level!=0)
  amrex::Error("level=0 in volWgtSumALL");

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else {
  std::cout << "SDC_outer_sweeps= " << SDC_outer_sweeps << '\n';
  amrex::Error("SDC_outer_sweeps invalid");
 }

    // vof,ref cen, order,slope,int
 int update_flag=0;  // do not update the error
 int init_vof_ls_prev_time=0;
 VOF_Recon_ALL(1,cur_time_slab,update_flag,
   init_vof_ls_prev_time,SLOPE_RECON_MF); 

 int project_option=1;  // initial project
  // need to initialize viscosity and density temporary 
  // variables.
  // in: volWgtSumALL
 int post_restart_flag=0;
 if ((post_init_flag==0)||(post_init_flag==1)) {
  // do nothing
 } else if (post_init_flag==2) {
  post_restart_flag=1;
 } else
  amrex::Error("post_init_flag invalid");

 make_physics_varsALL(project_option,post_restart_flag,0);

  // force,torque,moment of inertia,COM,mass
 allocate_array(0,4*AMREX_SPACEDIM+1,-1,DRAG_MF);
 Vector<Real> integrated_quantities;
 integrated_quantities.resize(13); 
 GetDragALL(integrated_quantities);

 int do_alloc=1;
 int simple_AMR_BC_flag_viscosity=1;
 init_gradu_tensorALL(HOLD_VELOCITY_DATA_MF,do_alloc,
   CELLTENSOR_MF,FACETENSOR_MF,
   simple_AMR_BC_flag_viscosity);

 for (int ilev = 0; ilev <= finest_level; ilev++) {

  NavierStokes& ns_level = getLevel(ilev);
  ns_level.volWgtSum(
    result,
    sumdata_type,
    sumdata_sweep,
    ZZ,FF,
    dirx,diry,cut_flag,
    ns_level.localMF[DRAG_MF],isweep);

 }  // ilev 

 delete_array(DRAG_MF);
 delete_array(CELLTENSOR_MF);
 delete_array(FACETENSOR_MF);

}  // subroutine volWgtSumALL

void
NavierStokes::MaxPressureVelocityALL(
   Real& minpres,Real& maxpres,Real& maxvel) {

 int finest_level=parent->finestLevel();
 if (level!=0)
  amrex::Error("level=0 in MaxPressureVelocityALL");

 Real local_minpres;
 Real local_maxpres;
 Real local_maxvel;
 maxpres=-1.0e+99;
 minpres=1.0e+99;
 maxvel=0.0;
 for (int k = 0; k <= finest_level; k++) {
  NavierStokes& ns_level = getLevel(k);
  ns_level.MaxPressureVelocity(local_minpres,local_maxpres,local_maxvel);
  if (local_minpres<minpres)
   minpres=local_minpres;
  if (local_maxpres>maxpres)
   maxpres=local_maxpres;
  if (local_maxvel>maxvel)
   maxvel=local_maxvel;
 }

} // end subroutine MaxPressureVelocityALL

void 
NavierStokes::MaxPressureVelocity(Real& minpres,Real& maxpres,Real& maxvel) {
 
 bool use_tiling=ns_tiling;

  // ngrow=0  
 MultiFab* vel=getState(0,0,
   num_materials_vel*(AMREX_SPACEDIM+1),cur_time_slab);
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

 Vector<Real> minpresA;
 Vector<Real> maxpresA;
 Vector<Real> maxvelA;
 minpresA.resize(thread_class::nthreads);
 maxpresA.resize(thread_class::nthreads);
 maxvelA.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  minpresA[tid]=1.0e+99;
  maxpresA[tid]=-1.0e+99;
  maxvelA[tid]=0.0;
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
  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  FORT_MAXPRESVEL(
   &minpresA[tid_current],
   &maxpresA[tid_current],
   &maxvelA[tid_current],
   xlo,dx,
   maskfab.dataPtr(),ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
   velfab.dataPtr(),ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
   tilelo,tilehi,
   fablo,fabhi,&bfact);
 } // mfi
} // omp
 ns_reconcile_d_num(107);

 for (int tid=0;tid<thread_class::nthreads;tid++) {
  if (minpres>minpresA[tid])
   minpres=minpresA[tid];
  if (maxpres<maxpresA[tid])
   maxpres=maxpresA[tid];
  if (maxvel<maxvelA[tid])
   maxvel=maxvelA[tid];
 }
 ParallelDescriptor::ReduceRealMin(minpres);
 ParallelDescriptor::ReduceRealMax(maxpres);
 ParallelDescriptor::ReduceRealMax(maxvel);

 delete mask;
 delete vel;
} // subroutine MaxPressureVelocity

// called from 
// 1. writePlotFile   (post_init_flag=0)
// 2. post_init_state (post_init_flag=1)
// 3. post_restart    (post_init_flag=2)
void
NavierStokes::prepare_post_process(int post_init_flag) {

 if (level!=0)
  amrex::Error("level invalid prepare_post_process");

 const int finest_level = parent->finestLevel();

  // init VOLUME_MF and AREA_MF
 metrics_dataALL(1);

 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);

  ns_level.allocate_mdot();

    // mask=tag if not covered by level+1 or outside the domain.
  Real tag=1.0;
  int clearbdry=0; 
   // ngrow=1
  ns_level.maskfiner_localMF(MASKCOEF_MF,1,tag,clearbdry);
  ns_level.prepare_mask_nbr(1);

  if (post_init_flag==1) { // called from post_init_state
   // do nothing
  } else if (post_init_flag==2) { // called from post_restart
   // do nothing
  } else if (post_init_flag==0) { // called from writePlotFile
    // in: NavierStokes::prepare_post_process
   ns_level.allocate_levelsetLO(1,LEVELPC_MF);
  } else
   amrex::Error("post_init_flag invalid");
   
 } // ilev=level ... finest_level

 init_FSI_GHOST_MAC_MF_ALL(2);

 build_masksemALL();

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  if (ilev<finest_level) {
   ns_level.MOFavgDown();
  }
 } // ilev=finest_level ... level

   // in prepare_post_process:
   // post_init_flag==0 called from writePlotFile
   // post_init_flag==1 called from post_init_state
   // post_init_flag==2 called from post_restart
   //
 int init_vof_ls_prev_time=0;
 int error_update_flag=0;
 int renormalize_only=0; // init:solid TEMP,VEL,LS,extend LSfluid into solid.
 int local_truncate=0; // do not force removal of flotsam.

 if (post_init_flag==0) { // called from writePlotFile
  // do nothing
 } else if (post_init_flag==1) {
  error_update_flag=1;  // called from post_init_state, update S_old
 } else if (post_init_flag==2) {
  error_update_flag=1;  // called from post_restart, update S_old
 } else
  amrex::Error("post_init_flag invalid");
	
 if (post_init_flag==1) { // called from post_init_state
  VOF_Recon_ALL(1,cur_time_slab,error_update_flag,
   init_vof_ls_prev_time,SLOPE_RECON_MF);
  int keep_all_interfaces=1;
  makeStateDistALL(keep_all_interfaces);
  prescribe_solid_geometryALL(cur_time_slab,renormalize_only,local_truncate);
  int project_option=1;  // initial project
  int post_restart_flag=0;
  make_physics_varsALL(project_option,post_restart_flag,1);
 } else if (post_init_flag==0) { // called from writePlotFile
  if (1==1) {
   VOF_Recon_ALL(1,cur_time_slab,error_update_flag,
    init_vof_ls_prev_time,SLOPE_RECON_MF);
   int project_option=1;  // initial project
   int post_restart_flag=0;
   make_physics_varsALL(project_option,post_restart_flag,2);
  }
 } else if (post_init_flag==2) { // called from post_restart
  VOF_Recon_ALL(1,cur_time_slab,error_update_flag,
    init_vof_ls_prev_time,SLOPE_RECON_MF);
  if (1==0) {
   int keep_all_interfaces=0;
   makeStateDistALL(keep_all_interfaces);
   prescribe_solid_geometryALL(cur_time_slab,renormalize_only,local_truncate);
  }
  int project_option=1;  // initial project
  int post_restart_flag=1;
  make_physics_varsALL(project_option,post_restart_flag,3);
 } else
  amrex::Error("post_init_flag invalid");

}  // subroutine prepare_post_process

// called from: NavierStokes::tensor_advection_updateALL() (NavierStokes3.cpp)
//   (called after calls to ns_level.tensor_advection_update)
void 
NavierStokes::accumulate_PC_info(int im_elastic) {

 NavierStokes& ns_level0=getLevel(0);

 int nmat=num_materials;
 bool use_tiling=ns_tiling;
 int finest_level=parent->finestLevel();
 if (level==finest_level) {
  // do nothing
 } else
  amrex::Error("expecting level==finest_level");

 resize_levelsetLO(2,LEVELPC_MF);
 debug_ngrow(LEVELPC_MF,2,8);
 if (localMF[LEVELPC_MF]->nComp()!=nmat*(AMREX_SPACEDIM+1))
  amrex::Error("(localMF[LEVELPC_MF]->nComp()!=nmat*(AMREX_SPACEDIM+1))");

 if ((im_elastic>=0)&&(im_elastic<nmat)) {
  // do nothing
 } else
  amrex::Error("im_elastic invalid");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");
 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 const Real* dx = geom.CellSize();

 int ipart=0;
 if (particleLS_flag[im_elastic]==1) {
  for (int im_local=0;im_local<im_elastic;im_local++) {
   if (particleLS_flag[im_local]==1) {
    ipart++; 
   } else if (particleLS_flag[im_local]==0) {
    // do nothing
   } else
    amrex::Error("particleLS_flag invalid");
  } // im_local=0..im_elastic-1

  if (NS_ncomp_particles>=ipart+1) {
   // do nothing
  } else
   amrex::Error("NS_ncomp_particles invalid");
 } else if (particleLS_flag[im_elastic]==0) {
  // do nothing
 } else
  amrex::Error("particleLS_flag[im_elastic] invalid");

 MultiFab& Tensor_new=get_new_data(Tensor_Type,slab_step+1);

 int elastic_partid=0;
 int scomp_tensor=0;

 if (ns_is_rigid(im_elastic)==0) {
  if ((elastic_time[im_elastic]>0.0)&&
      (elastic_viscosity[im_elastic]>0.0)) {
   elastic_partid=0;
   while ((im_elastic_map[elastic_partid]!=im_elastic)&&
          (elastic_partid<im_elastic_map.size())) {
    elastic_partid++;
   }
   if (elastic_partid<im_elastic_map.size()) {
    scomp_tensor=elastic_partid*NUM_TENSOR_TYPE;
   } else
    amrex::Error("elastic_partid too large");
  } else
   amrex::Error("elastic_time or elastic_viscosity invalid");
 } else
  amrex::Error("ns_is_rigid invalid");

 int matrix_points=1;  // sum_{xp in Omega_cell} W(xp,x_cell,LS)
 int RHS_points=1;     // sum_{xp in Omega_cell} (X_cell(xp)-X_cell_p)*W
 int ncomp_accumulate=matrix_points+AMREX_SPACEDIM*RHS_points;
 MultiFab* accumulate_mf=new MultiFab(grids,dmap,ncomp_accumulate,0,
	  MFInfo().SetTag("accumulate_mf"),FArrayBoxFactory());
 accumulate_mf->setVal(0.0);

 const Vector<Geometry>& ns_geom=parent->Geom();
 const Vector<DistributionMapping>& ns_dmap=parent->DistributionMap();
 const Vector<BoxArray>& ns_ba=parent->boxArray();

 Vector<int> rr;
 rr.resize(ns_ba.size());
 for (int ilev=0;ilev<rr.size();ilev++)
  rr[ilev]=2;
 int nnbr=1;

 bool local_copy_flag=true; 

 if (localMF_grow[VISCOTEN_MF]==-1) {
  // do nothing
 } else 
  amrex::Error("VISCOTEN_MF should not be allocated");

 int scomp_xdisplace=num_materials_viscoelastic*NUM_TENSOR_TYPE;

 int ncomp_tensor=NUM_TENSOR_TYPE;

 if (particleLS_flag[im_elastic]==0) {

  getStateTensor_localMF(VISCOTEN_MF,2,scomp_xdisplace,AMREX_SPACEDIM,
   cur_time_slab);

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(accumulate_mf->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(*accumulate_mf,use_tiling); mfi.isValid(); ++mfi) {
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

   FArrayBox& TNEWfab=Tensor_new[mfi];
   FArrayBox& XDISP_fab=(*localMF[VISCOTEN_MF])[mfi];
   FArrayBox& levelpcfab=(*localMF[LEVELPC_MF])[mfi];
 
   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   fort_assimilate_tensor_from_xdisplace( 
     &im_elastic, // 0..nmat-1
     &tid_current,
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,
     &level,
     &finest_level,
     xlo,dx,
     &ncomp_tensor,
     &nmat,
     levelpcfab.dataPtr(),
     ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),
     TNEWfab.dataPtr(scomp_tensor),
     ARLIM(TNEWfab.loVect()),ARLIM(TNEWfab.hiVect()),
     XDISP_fab.dataPtr(),
     ARLIM(XDISP_fab.loVect()),ARLIM(XDISP_fab.hiVect()));
  } // mfi
} // omp
  ns_reconcile_d_num(81);

  delete_localMF(VISCOTEN_MF,1);

 } else if (particleLS_flag[im_elastic]==1) {

  // All the particles should live on level==0.
  // particle levelset==0.0 for interface particles.

  AmrParticleContainer<N_EXTRA_REAL,0,0,0>& localPC_no_nbr=
    ns_level0.get_new_dataPC(State_Type,slab_step+1,ipart);

  NeighborParticleContainer<N_EXTRA_REAL,0> 
   localPC(ns_geom,ns_dmap,ns_ba,rr,nnbr);
  // the two PC have same hierarchy, no need to call Redistribute after the
  // copy.
  localPC.copyParticles(localPC_no_nbr,local_copy_flag);

  localPC.fillNeighbors();

  for (int isweep=0;isweep<=1;isweep++) {

   getStateTensor_localMF(VISCOTEN_MF,2,scomp_xdisplace,AMREX_SPACEDIM,
    cur_time_slab);

   if (thread_class::nthreads<1)
    amrex::Error("thread_class::nthreads invalid");
   thread_class::init_d_numPts(accumulate_mf->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
   for (MFIter mfi(*accumulate_mf,use_tiling); mfi.isValid(); ++mfi) {
    BL_ASSERT(grids[mfi.index()] == mfi.validbox());
    const int gridno = mfi.index();
     // std::cout << tilegrid << '\n';
     // std::cout << gridno << '\n';
    const Box& tilegrid = mfi.tilebox();
    const Box& fabgrid = grids[gridno];
    const int* tilelo=tilegrid.loVect();
    const int* tilehi=tilegrid.hiVect();
    const int* fablo=fabgrid.loVect();
    const int* fabhi=fabgrid.hiVect();
    int bfact=parent->Space_blockingFactor(level);

    const Real* xlo = grid_loc[gridno].lo();

    // particles is of type "ParticleTileType::AoS"
    // see AMReX_ArrayOfStructs.H
    // particles is of type:
    //  amrex::ParticleTile<SDIM,0,0,0>
    //

    auto& particles = localPC.GetParticles(level)
     [std::make_pair(mfi.index(),mfi.LocalTileIndex())];
    auto& particles_AoS = particles.GetArrayOfStructs();

    int Np=particles_AoS.size();

     // ParticleVector&
    auto& neighbors_local = 
	localPC.GetNeighbors(level,mfi.index(),mfi.LocalTileIndex());
    int Nn=neighbors_local.size();

    FArrayBox& matrixfab=(*accumulate_mf)[mfi];
    FArrayBox& TNEWfab=Tensor_new[mfi];
    FArrayBox& XDISP_fab=(*localMF[VISCOTEN_MF])[mfi];
    FArrayBox& levelpcfab=(*localMF[LEVELPC_MF])[mfi];
 
    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // in: GODUNOV_3D.F90
    // updates (1) configuration tensor and
    // (2) XDISPLACE data.
    fort_assimilate_tensor_from_particles( 
     &im_elastic, // 0..nmat-1
     &isweep,
     &tid_current,
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,
     &level,
     &finest_level,
     xlo,dx,
     particles_AoS.data(),
     neighbors_local.data(),
     Np,       //pass by value
     Nn,       //pass by value
     &ncomp_tensor,
     &matrix_points,
     &RHS_points,
     &ncomp_accumulate,
     &nmat,
     levelpcfab.dataPtr(),
     ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),
     TNEWfab.dataPtr(scomp_tensor),
     ARLIM(TNEWfab.loVect()),ARLIM(TNEWfab.hiVect()),
     TNEWfab.dataPtr(scomp_xdisplace),
     ARLIM(TNEWfab.loVect()),ARLIM(TNEWfab.hiVect()),
     XDISP_fab.dataPtr(),
     ARLIM(XDISP_fab.loVect()),ARLIM(XDISP_fab.hiVect()),
     matrixfab.dataPtr(),
     ARLIM(matrixfab.loVect()),ARLIM(matrixfab.hiVect()));
   } // mfi
} // omp
   ns_reconcile_d_num(81);

   delete_localMF(VISCOTEN_MF,1);
  } // isweep=0,1

  localPC.clearNeighbors();

 } else {

  amrex::Error("particleLS_flag[im_elastic] invalid");

 }

 delete accumulate_mf;

} // end subroutine accumulate_PC_info


// initialize particles and copy to all "slab locations"
// ONLY LEVEL==max_level STATEDATA PARTICLES GET INITIALIZED: THEY HOLD
// PARTICLES THAT APPEAR ON ALL THE LEVELS.
// ALSO: Only state[State_Type] has the particles.
// DO NOT FORGET TO HAVE CHECKPOINT/RESTART CAPABILITY FOR PARTICLES.
// This routine called for level=finest_level 
void
NavierStokes::init_particle_container(int im_PLS,int ipart,int append_flag) {

 NavierStokes& ns_level0=getLevel(0);

 bool use_tiling=ns_tiling;
 int max_level = parent->maxLevel();
 int finest_level=parent->finestLevel();

 if (slab_step==ns_time_order-1) {
  // do nothing
 } else
  amrex::Error("expecting slab_step==ns_time_order-1");

 if (level<=finest_level) {
  // do nothing
 } else 
  amrex::Error("level<=finest_level failed");

 int nmat=num_materials;
 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");
 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 const Real* dx = geom.CellSize();

 if (NS_ncomp_particles>0) {

  if ((im_PLS>=0)&&(im_PLS<nmat)) {
   // do nothing
  } else
   amrex::Error("im_PLS invalid");

  if (particleLS_flag[im_PLS]==1) {

   MultiFab* LSmf=getStateDist(1,cur_time_slab,7);  
   if (LSmf->nComp()!=nmat*(1+AMREX_SPACEDIM))
    amrex::Error("LSmf invalid ncomp");
   if (LSmf->nGrow()!=1)
    amrex::Error("LSmf->nGrow()!=1");

    // All the particles should live on level==0.
    // particle levelset==0.0 for interface particles.
   AmrParticleContainer<N_EXTRA_REAL,0,0,0>& localPC=
     ns_level0.get_new_dataPC(State_Type,slab_step+1,ipart);

    // ngrow=1
    // scomp=num_materials_viscoelastic*NUM_TENSOR_TYPE
    // ncomp=sdim
   MultiFab* x_displace_mf=getStateTensor(1,
      num_materials_viscoelastic*NUM_TENSOR_TYPE,
      AMREX_SPACEDIM,cur_time_slab);

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
    FArrayBox& xdisplacefab=(*x_displace_mf)[mfi];

     // component 1: number of particles linked to the cell.
     // component 2: the link to the list of particles.
    BaseFab<int> cell_particle_count(tilegrid,2);
    cell_particle_count.setVal(0);

     // allocate for just one particle for now.
    int single_particle_size=AMREX_SPACEDIM+N_EXTRA_REAL;
    Vector< Real > new_particle_data;
    new_particle_data.resize(single_particle_size);

      // this is an object with a pointer to both AoS and
      // SoA data.
    auto& particles = localPC.GetParticles(level)
      [std::make_pair(mfi.index(),mfi.LocalTileIndex())];

    auto& particles_AoS = particles.GetArrayOfStructs();
    int Np=particles_AoS.size();

     // The link index will start at 1.
    Vector< int > particle_link_data;
     // i_particle_link_1,i1,j1,k1,   (child link, parent link)
     // i_particle_link_2,i2,j2,k2,  ...
    particle_link_data.resize(Np*(1+AMREX_SPACEDIM));
    for (int i_link=0;i_link<Np*(1+AMREX_SPACEDIM);i_link++)
     particle_link_data[i_link]=0;

      // 1 if particle should be deleted.
    Vector< int > particle_delete_flag; 
    particle_delete_flag.resize(Np);
    for (int i_delete=0;i_delete<Np;i_delete++)
     particle_delete_flag[i_delete]=0;

    int Np_append=0;  // number of particles to append

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    for (int isweep=0;isweep<2;isweep++) {

     int new_Pdata_size=new_particle_data.size();

     // in: LEVELSET_3D.F90
     // 1. subdivide each cell with "particle_nsubdivide" divisions.
     //    e.g. if particle_nsubdivide=2 => 4 pieces in 2D.
     //                 "         "   =4 => 64 pieces in 2D.
     // 2. for each small sub-box, add a particles at the sub-box center
     //    and add a particle "x-phi grad phi/|grad phi|"
     fort_init_particle_container( 
       &tid_current,
       &single_particle_size,
       &isweep,
       &append_flag,
       particle_nsubdivide.dataPtr(),
       particle_max_per_nsubdivide.dataPtr(),
       particleLS_flag.dataPtr(),
       &im_PLS,
       &nmat,
       tilelo,tilehi,
       fablo,fabhi,&bfact,
       &level,
       &finest_level,
       &cur_time_slab,
       xlo,dx,
       particles_AoS.data(), // existing particles
       Np,  // pass by value
       new_particle_data.dataPtr(), // size is "new_Pdata_size"
       &new_Pdata_size,
       &Np_append,  // Np_append number of new particles to add.
       particle_link_data.dataPtr(),
       particle_delete_flag.dataPtr(),
       cell_particle_count.dataPtr(),
       ARLIM(cell_particle_count.loVect()),
       ARLIM(cell_particle_count.hiVect()),
       xdisplacefab.dataPtr(),
       ARLIM(xdisplacefab.loVect()),ARLIM(xdisplacefab.hiVect()),
       lsfab.dataPtr(),
       ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()) );

     if (isweep==0) {
      new_particle_data.resize(Np_append*single_particle_size);
     }
    } // isweep=0...1

    using My_ParticleContainer =
      AmrParticleContainer<N_EXTRA_REAL,0,0,0>;

    if (debug_PC==1) {
     std::cout << "gridno " << gridno << " Np_append " <<
       Np_append << '\n';
    }

    int Np_delete=0;
    for (int i_delete=0;i_delete<Np;i_delete++) {
     if (particle_delete_flag[i_delete]==1) {
      Np_delete++;
     } else if (particle_delete_flag[i_delete]==0) {
      // do nothing
     } else
      amrex::Error("particle_delete_flag[i_delete] invalid");
    }
    if (Np_delete<=Np) {
     // do nothing
    } else
     amrex::Error("Np_delete invalid");

     // mirrorPC will only contain AoS data.  In the future,
     // SoA data will have to be managed too.
    Vector< My_ParticleContainer::ParticleType > mirrorPC_AoS;
    int Np_mirror_AoS=Np-Np_delete+Np_append;
    mirrorPC_AoS.resize(Np_mirror_AoS);
    int i_mirror=0;
    for (int i_delete=0;i_delete<Np;i_delete++) {
     if (particle_delete_flag[i_delete]==1) {
      // do nothing
     } else if (particle_delete_flag[i_delete]==0) {
      mirrorPC_AoS[i_mirror]=particles_AoS[i_delete];
      i_mirror++;
     } else
      amrex::Error("particle_delete_flag[i_delete] invalid");
    }

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
     mirrorPC_AoS[i_mirror]=p;
     i_mirror++;
    } // i_append=0..Np_append-1
    if (i_mirror==Np_mirror_AoS) {
     particles.resize(0);
     for (i_mirror=0;i_mirror<Np_mirror_AoS;i_mirror++) {
      particles.push_back(mirrorPC_AoS[i_mirror]);
     }
    } else
     amrex::Error("i_mirror <> Np_mirror_AoS");

    if (ipart<NS_ncomp_particles) {
     // do nothing
    } else
     amrex::Error("ipart invalid");

   } // mfi
} // omp
   ns_reconcile_d_num(81);

   delete x_displace_mf;
   delete LSmf;

  } else
   amrex::Error("particleLS_flag[im_PLS] invalid");

 } else
  amrex::Error("NS_ncomp_particles invalid");

}  // end subroutine init_particle_container()

 
// should be cur_time=0 and prev_time=-1
// called from post_init
void
NavierStokes::post_init_state () {
    
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

 int project_option_combine=-1;
 int combine_flag=2;  
 int hflag=0;
 int combine_idx=-1;
 int update_flux=0;

 int nmat=num_materials;

 int project_option=1;  // initial project

 const int finest_level = parent->finestLevel();

 NavierStokes& ns_finest=getLevel(finest_level);

   // inside of post_init_state

   // metrics_data
   // allocate_mdot
   // MASKCOEF
   // init_FSI_GHOST_MAC_MF
   // VOF_Recon_ALL (update_flag==1)
   // makeStateDistALL
   // prescribe_solid_geometryALL
   // make_physics_varsALL
 int post_init_flag=1; // in: post_init_state
 prepare_post_process(post_init_flag);

 if (NS_ncomp_particles>0) {

  int ipart=0;
  for (int im=0;im<nmat;im++) {

   if (particleLS_flag[im]==1) {

    if (debug_PC==1) {
     std::cout << "PC: slab_step, ns_time_order, im, ipart " <<
      slab_step << ' ' << ns_time_order << ' ' << im << ' ' << ipart << '\n';
    }

    NavierStokes& ns_level0=getLevel(0);
    AmrParticleContainer<N_EXTRA_REAL,0,0,0>& localPC=
     ns_level0.get_new_dataPC(State_Type,slab_step+1,ipart);
    
    int append_flag=0;
    ns_finest.init_particle_container(im,ipart,append_flag);

    int lev_min=0;
    int lev_max=-1;
    int nGrow_Redistribute=0;
    int local_Redistribute=0;
    localPC.Redistribute(lev_min,lev_max,nGrow_Redistribute, 
     local_Redistribute);

    ipart++;
   } else if (particleLS_flag[im]==0) {
    // do nothing
   } else
    amrex::Error("particleLS_flag[im] invalid");

  } // im=0..nmat-1

  if (ipart==NS_ncomp_particles) {
   // do nothing
  } else
   amrex::Error("ipart invalid");
 
 } else if (NS_ncomp_particles==0) {
  // do nothing
 } else
  amrex::Error("NS_ncomp_particles invalid");


 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);

  project_option_combine=2;  // temperature in post_init_state
  combine_flag=2;  
  hflag=0;
  combine_idx=-1;
  update_flux=0;

  ns_level.combine_state_variable(
    project_option_combine,
    combine_idx,combine_flag,hflag,update_flux); 

  for (int ns=0;ns<num_species_var;ns++) {
   project_option_combine=100+ns; // species in post_init_state
   ns_level.combine_state_variable(
    project_option_combine,
    combine_idx,combine_flag,hflag,update_flux); 
  }

 } // ilev=finest_level ... level

 Vector<blobclass> blobdata;

 int color_count=0;
 int coarsest_level=0;

  // tessellate==1
 ColorSumALL(coarsest_level,color_count,TYPE_MF,COLOR_MF,blobdata);
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
  // in post_init_state (project_option==1)
  // if project_option==1, then the velocity in the ice
  // is overwritten with a projected rigid body velocity.
 int interp_option=0;
 int idx_velcell=-1;
 Real beta=0.0;

 increment_face_velocityALL(
   interp_option,project_option,
   idx_velcell,beta,blobdata); 

 delete_array(TYPE_MF);
 delete_array(COLOR_MF);

 if ((post_init_pressure_solve==1)&&
     (!is_zalesak())) { 

  double start_pressure_solve = ParallelDescriptor::second();

   // U^CELL and U^MAC; assimilates the solid velocity (MAC and CELL)
   // and  ice velocity (MAC)
  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   project_option_combine=3; // velocity in post_init_state
   combine_flag=2;
   hflag=0;
   combine_idx=-1;
   update_flux=0;
   ns_level.combine_state_variable(
    project_option_combine,
    combine_idx,combine_flag,hflag,update_flux);
   project_option_combine=0; // mac velocity
   update_flux=1;
   ns_level.combine_state_variable(
    project_option_combine,
    combine_idx,combine_flag,hflag,update_flux);

   ns_level.make_MAC_velocity_consistent();
  }

  multiphase_project(project_option); // initial project

   // U^CELL and U^MAC
  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);

   project_option_combine=3; // velocity in post_init_state
   combine_flag=2;
   hflag=0;
   combine_idx=-1;
   update_flux=0;
   ns_level.combine_state_variable(
    project_option_combine,
    combine_idx,combine_flag,hflag,update_flux);
   project_option_combine=0; // mac velocity
   update_flux=1;
   ns_level.combine_state_variable(
    project_option_combine,
    combine_idx,combine_flag,hflag,update_flux);

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

 int scomp_pres=num_materials_vel*AMREX_SPACEDIM;
 int scomp_den=num_materials_vel*(AMREX_SPACEDIM+1);
 for (int ilev=finest_level;ilev>=0;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.avgDown(State_Type,scomp_pres,num_materials_vel,1);
  ns_level.avgDown(State_Type,scomp_den,num_state_material*nmat,1);
 } // ilev

 delete_array(MASKCOEF_MF);
 delete_array(MASK_NBR_MF);

 CopyNewToOldALL();

}  // subroutine post_init_state

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
  amrex::Error("S_crse invalid");
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

  FORT_AVGDOWN_TAG(
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
 ns_reconcile_d_num(108);
 S_crse.copy(crse_S_fine,0,scomp,ncomp);
 ParallelDescriptor::Barrier();
} // subroutine level_avgDown_tag


void
NavierStokes::level_avgDownBURNING(MultiFab& S_crse,MultiFab& S_fine, 
		int velflag) {

 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;
 int ncomp_per_burning=AMREX_SPACEDIM;
 int ncomp_per_tsat=2;
 int nburning=nten*(ncomp_per_burning+1);
 int ntsat=nten*(ncomp_per_tsat+1);
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
  amrex::Error("S_crse invalid");
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
  FORT_AVGDOWN_BURNING(
   &velflag,
   prob_lo,
   dxf,
   &level,&f_level,
   &bfact_c,&bfact_f,
   xlo_fine,dx,
   &ncomp,
   &nmat,
   &nten,
   c_dat,ARLIM(clo),ARLIM(chi),
   f_dat,ARLIM(flo),ARLIM(fhi),
   ovlo,ovhi,
   fine_fablo,fine_fabhi);

 }// mfi
} //omp
 ns_reconcile_d_num(109);
 S_crse.copy(crse_S_fine,0,scomp,ncomp);
 ParallelDescriptor::Barrier();
} // subroutine level_avgDownBURNING


void
NavierStokes::level_avgDownCURV(MultiFab& S_crse,MultiFab& S_fine) {

 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

 int scomp=0;
 int ncomp=nten*(AMREX_SPACEDIM+5);

 int finest_level=parent->finestLevel();
 if (level>=finest_level)
  amrex::Error("level invalid in level_avgDownCURV");

 int f_level=level+1;
 NavierStokes&   fine_lev = getLevel(f_level);
 const BoxArray& fgrids=fine_lev.grids;
 const DistributionMapping& fdmap=fine_lev.dmap;

 if (grids!=S_crse.boxArray())
  amrex::Error("S_crse invalid");
 if (fgrids!=S_fine.boxArray())
  amrex::Error("S_fine invalid");
 if (S_crse.nComp()!=S_fine.nComp())
  amrex::Error("nComp mismatch");
 if (S_crse.nComp()!=scomp+ncomp)
  amrex::Error("S_crse.nComp() invalid level_avgDownCURV");

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

  FORT_AVGDOWN_CURV(
   prob_lo,
   dxf,
   &level,&f_level,
   &bfact_c,&bfact_f,
   xlo_fine,dx,
   &ncomp,&nmat,&nten,
   c_dat,ARLIM(clo),ARLIM(chi),
   f_dat,ARLIM(flo),ARLIM(fhi),
   ovlo,ovhi,
   fine_fablo,fine_fabhi);

 }// mfi
} //omp
 ns_reconcile_d_num(110);
 S_crse.copy(crse_S_fine,0,scomp,ncomp);
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

 if (grids!=S_crse.boxArray())
  amrex::Error("S_crse invalid");
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
   FORT_AVGDOWN( 
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

   FORT_AVGDOWN_LOW(
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
 ns_reconcile_d_num(111);
 S_crse.copy(crse_S_fine,0,scomp,ncomp);
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

 int finest_level=parent->finestLevel();

 if (level == finest_level)
  return;

 int f_level=level+1;
 NavierStokes&   fine_lev = getLevel(f_level);
 const BoxArray& fgrids=fine_lev.grids;
 const DistributionMapping& fdmap=fine_lev.dmap;
 resize_metrics(1);
 debug_ngrow(VOLUME_MF,0,80);
 fine_lev.resize_metrics(1);
 fine_lev.debug_ngrow(VOLUME_MF,1,81);

 MultiFab& S_fine=fine_lev.get_new_data(State_Type,slab_step+1);
 MultiFab& S_crse = get_new_data(State_Type,slab_step+1);

 const Real* dxf = fine_lev.geom.CellSize();
 const Real* dxc = geom.CellSize();
 const Real* prob_lo   = geom.ProbLo();

 if (grids!=S_crse.boxArray())
  amrex::Error("S_crse invalid");
 if (fgrids!=S_fine.boxArray())
  amrex::Error("S_fine invalid");
 if (S_crse.nComp()!=S_fine.nComp())
  amrex::Error("nComp mismatch");

 BoxArray crse_S_fine_BA(fgrids.size());
 for (int i = 0; i < fgrids.size(); ++i) {
  crse_S_fine_BA.set(i,amrex::coarsen(fgrids[i],2));
 }

 int nmat=num_materials;
 int scomp_mofvars=num_materials_vel*(AMREX_SPACEDIM+1)+
   nmat*num_state_material;

 DistributionMapping crse_dmap=fdmap;
 MultiFab crse_S_fine(crse_S_fine_BA,crse_dmap,nmat*ngeom_raw,0,
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
  const Real* f_dat=finefab.dataPtr(scomp_mofvars);

  FArrayBox& coarsefab=crse_S_fine[gridno];
  const Box& cgrid = coarsefab.box();
  const int* clo=cgrid.loVect();
  const int* chi=cgrid.hiVect();
  const Real* c_dat=coarsefab.dataPtr();

  int bfact_c=parent->Space_blockingFactor(level);
  int bfact_f=parent->Space_blockingFactor(f_level);

  FORT_MOFAVGDOWN(
   &cur_time_slab,
   prob_lo,
   dxc,
   dxf,
   &bfact_c,&bfact_f,
   c_dat,ARLIM(clo),ARLIM(chi),
   f_dat,ARLIM(flo),ARLIM(fhi),
   ovlo,ovhi,&nmat);
 } // mfi
} //omp
 ns_reconcile_d_num(112);
 S_crse.copy(crse_S_fine,0,scomp_mofvars,nmat*ngeom_raw);
 ParallelDescriptor::Barrier();
}


void NavierStokes::avgDownError() {

 int finest_level=parent->finestLevel();

 if (level == finest_level)
  return;

 int f_level=level+1;
 NavierStokes&   fine_lev = getLevel(f_level);
 const BoxArray& fgrids=fine_lev.grids;
 const DistributionMapping& fdmap=fine_lev.dmap;
 resize_metrics(1);
 debug_ngrow(VOLUME_MF,0,80);
 fine_lev.resize_metrics(1);
 fine_lev.debug_ngrow(VOLUME_MF,1,81);

 MultiFab& S_fine=fine_lev.get_new_data(State_Type,slab_step+1);
 MultiFab& S_crse = get_new_data(State_Type,slab_step+1);

 const Real* dxf = fine_lev.geom.CellSize();
 const Real* prob_lo   = geom.ProbLo();

 if (grids!=S_crse.boxArray())
  amrex::Error("S_crse invalid");
 if (fgrids!=S_fine.boxArray())
  amrex::Error("S_fine invalid");
 if (S_crse.nComp()!=S_fine.nComp())
  amrex::Error("nComp mismatch");

 BoxArray crse_S_fine_BA(fgrids.size());
 for (int i = 0; i < fgrids.size(); ++i) {
  crse_S_fine_BA.set(i,amrex::coarsen(fgrids[i],2));
 }

 int nmat=num_materials;
 int scomp_error=num_materials_vel*(AMREX_SPACEDIM+1)+
   nmat*num_state_material+nmat*ngeom_raw;
 if (S_crse.nComp()!=scomp_error+1)
  amrex::Error("scomp_error invalid");

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
  const Real* f_dat=finefab.dataPtr(scomp_error);

  FArrayBox& coarsefab=crse_S_fine[gridno];
  const Box& cgrid = coarsefab.box();
  const int* clo=cgrid.loVect();
  const int* chi=cgrid.hiVect();
  const Real* c_dat=coarsefab.dataPtr();

  int bfact_c=parent->Space_blockingFactor(level);
  int bfact_f=parent->Space_blockingFactor(f_level);

  FORT_ERRORAVGDOWN(
   prob_lo,
   dxf,
   &bfact_c,&bfact_f,
   c_dat,ARLIM(clo),ARLIM(chi),
   f_dat,ARLIM(flo),ARLIM(fhi),
   ovlo,ovhi);
 } // mfi
} //omp
 ns_reconcile_d_num(113);
 S_crse.copy(crse_S_fine,0,scomp_error,1);
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

 int nmat=num_materials;
 int scomp_mofvars=num_materials_vel*(AMREX_SPACEDIM+1)+
  nmat*num_state_material;
 int scomp_error=scomp_mofvars+nmat*ngeom_raw;

 if ((scomp<scomp_error)&&(scomp+ncomp-1>=scomp_mofvars)) {

  if ((scomp<=scomp_mofvars)&&(scomp+ncomp>=scomp_error)) {
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

 FillPatch(*this,*mf,0,time,State_Type,scomp,ncomp);

 ParallelDescriptor::Barrier();

 return mf;
} // end subroutine getState

MultiFab* NavierStokes::getStateSolid (
  int ngrow, int  scomp,
  int ncomp, Real time) {

 int nmat=num_materials;

  // nparts x (velocity + LS + temperature + flag)
 int nparts=im_solid_map.size();
 if ((nparts<1)||(nparts>nmat))
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

 FillPatch(*this,*mf,0,time,Solid_State_Type,scomp,ncomp);

 ParallelDescriptor::Barrier();

 return mf;
} // end subroutine getStateSolid



MultiFab* NavierStokes::getStateTensor (
  int ngrow, int  scomp,
  int ncomp, Real time) {

 int nmat=num_materials;

  // nparts x NUM_TENSOR_TYPE
 int nparts=im_elastic_map.size();
 if ((nparts<0)||(nparts>nmat))
  amrex::Error("nparts invalid");
 if ((ncomp==nparts*NUM_TENSOR_TYPE+AMREX_SPACEDIM)&&
     (scomp==0)) {
  // do nothing
 } else if ((ncomp==AMREX_SPACEDIM)&&
            (scomp==nparts*NUM_TENSOR_TYPE)) {
  // do nothing
 } else if (ncomp%NUM_TENSOR_TYPE==0) {
  int partid=scomp/NUM_TENSOR_TYPE;
  if ((partid<0)||(partid>=nparts))
   amrex::Error("partid invalid");
 } else
  amrex::Error("ncomp or scomp invalid");

 MultiFab& Tensor_new=get_new_data(Tensor_Type,slab_step+1);
 int ntotal=Tensor_new.nComp();
 if (scomp<0)
  amrex::Error("scomp invalid getStateTensor"); 
 if (ncomp<=0)
  amrex::Error("ncomp invalid in getstateTensor"); 
 if (scomp+ncomp>ntotal)
  amrex::Error("scomp,ncomp invalid");

 MultiFab* mf = new MultiFab(state[Tensor_Type].boxArray(),dmap,ncomp,
   ngrow,MFInfo().SetTag("mf getStateTensor"),FArrayBoxFactory());

 FillPatch(*this,*mf,0,time,Tensor_Type,scomp,ncomp);

 ParallelDescriptor::Barrier();

 return mf;
} // end subroutine getStateTensor


MultiFab* NavierStokes::getStateDist (int ngrow,Real time,int caller_id) {

 if (verbose>0)
  if (ParallelDescriptor::IOProcessor())
   std::cout << "getStateDist: time,caller_id " << time << ' ' << 
     caller_id <<'\n';

 if ((ngrow<0)||(ngrow>ngrow_distance))
  amrex::Error("ngrow invalid");

 int nmat=num_materials;
 
 MultiFab& S_new=get_new_data(LS_Type,slab_step+1);
 int ntotal=S_new.nComp();
 if (ntotal!=nmat*(AMREX_SPACEDIM+1))
  amrex::Error("ntotal invalid");

 MultiFab* mf = new MultiFab(state[State_Type].boxArray(),dmap,
   nmat*(AMREX_SPACEDIM+1),
   ngrow,MFInfo().SetTag("mf getStateDist"),FArrayBoxFactory());

  // scomp=0
 FillPatch(*this,*mf,0,time,LS_Type,0,nmat*(AMREX_SPACEDIM+1));

 ParallelDescriptor::Barrier();

 return mf;
} // subroutine getStateDist

void NavierStokes::CPP_OVERRIDEPBC(int homflag_in,int project_option_in) {

 if ((homflag_in==0)||(homflag_in==1)) {
  override_bc_to_homogeneous=homflag_in;
  FORT_OVERRIDEPBC(&homflag_in,&project_option_in);
 } else
  amrex::Error("homflag_in invalid");

}  // subroutine CPP_OVERRIDEPBC

MultiFab* NavierStokes::getStateDIV_DATA(int ngrow,
		int scomp,int ncomp,Real time) {

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel!=1");

 int project_option=10;
 int save_bc_status=override_bc_to_homogeneous;
 CPP_OVERRIDEPBC(1,project_option);

 MultiFab& S_new=get_new_data(DIV_Type,slab_step+1);
 int ntotal=S_new.nComp();
 if (ntotal!=num_materials_vel)
  amrex::Error("ntotal invalid");

 MultiFab* mf = new MultiFab(state[DIV_Type].boxArray(),dmap,ncomp,
   ngrow,MFInfo().SetTag("mf getStateDIV_DATA"),FArrayBoxFactory());

 FillPatch(*this,*mf,0,time,DIV_Type,scomp,ncomp);

 ParallelDescriptor::Barrier();

 CPP_OVERRIDEPBC(save_bc_status,project_option);

 return mf;
} // subroutine getStateDIV_DATA


void NavierStokes::putStateDIV_DATA(
 int scomp,int ncomp,int idx_MF) {

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel!=1");

 MultiFab& S_new=get_new_data(DIV_Type,slab_step+1);
 int ntotal=S_new.nComp();
 if (ntotal!=num_materials_vel)
  amrex::Error("ntotal invalid");

 int scomp_localMF=0;
   // dst,src,scomp,dcomp,ncomp,ngrow
 MultiFab::Copy(S_new,*localMF[idx_MF],scomp_localMF,scomp,
  	        ncomp,0); 

} // subroutine putStateDIV_DATA



// called from:
// NavierStokes::prepare_post_process  (post_init_flag==1)
// NavierStokes::do_the_advance
void
NavierStokes::makeStateDistALL(int keep_all_interfaces) {

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
 max_problen=sqrt(max_problen);
 if (max_problen>0.0) {
  // do nothing
 } else
  amrex::Error("max_problen invalid");

 int nmat=num_materials;
 minLS.resize(thread_class::nthreads);
 maxLS.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  minLS[tid].resize(nmat);
  maxLS[tid].resize(nmat);
  for (int im=0;im<nmat;im++) {
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
  ns_level.makeStateDist(keep_all_interfaces);
 }
  // CORRECT_UNINIT is in MOF_REDIST_3D.F90
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
   for (int im=0;im<nmat;im++) {
    std::cout << "im= " << im << '\n';
    std::cout << "minLS = " << minLS[0][im] << '\n';
    std::cout << "maxLS = " << maxLS[0][im] << '\n';
   }

} // subroutine makeStateDistALL()

// called from: NavierStokes::do_the_advance 
// (prior to level_phase_change_rate) and
// called from: NavierStokes::init_FSI_GHOST_MAC_MF
void 
NavierStokes::build_NRM_FD_MF(int fd_mf,int ls_mf,int ngrow) {

 int nmat=num_materials;
 bool use_tiling=ns_tiling;
 int finest_level=parent->finestLevel();
 const Real* dx = geom.CellSize();

 if (ngrow>=0) {
  if ((localMF[fd_mf]->nGrow()>=ngrow)&&
      (localMF[ls_mf]->nGrow()>=ngrow)) {
   if ((localMF[fd_mf]->nComp()==nmat*AMREX_SPACEDIM)&&
       (localMF[ls_mf]->nComp()==nmat*(AMREX_SPACEDIM+1))) {
    MultiFab::Copy(*localMF[fd_mf],*localMF[ls_mf],nmat,0,
		    nmat*AMREX_SPACEDIM,ngrow);

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

      // in: MOF_REDIST_3D.F90
     FORT_FD_NORMAL( 
      &level,
      &finest_level,
      lsfab.dataPtr(),
      ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
      nrmfab.dataPtr(),
      ARLIM(nrmfab.loVect()),ARLIM(nrmfab.hiVect()),
      tilelo,tilehi,
      fablo,fabhi,
      &bfact,
      xlo,dx,
      &nmat);
    } // mfi
} // omp
    ns_reconcile_d_num(114);

    localMF[fd_mf]->FillBoundary(geom.periodicity()); 
   } else
    amrex::Error("fd_mf or ls_mf nComp() invalid");
  } else
   amrex::Error("fd_mf or ls_mf nGrow() invalid");
 } else
  amrex::Error("ngrow invalid");
			 
} // subroutine build_NRM_FD_MF

void
NavierStokes::makeStateDist(int keep_all_interfaces) {

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

 int rzflag=0;
 if (geom.IsRZ())
  rzflag=1;
 else if (geom.IsCartesian())
  rzflag=0;
 else if (geom.IsCYLINDRICAL())
  rzflag=3;
 else
  amrex::Error("CoordSys bust 6");

 const Real* dx = geom.CellSize();
 int nmat=num_materials;
 int scomp_mofvars=num_materials_vel*(AMREX_SPACEDIM+1)+
  nmat*num_state_material;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

 MultiFab& LS_new = get_new_data(LS_Type,slab_step+1);

 getStateDist_localMF(ORIGDIST_MF,ngrow_distance,cur_time_slab,15);

 resize_mask_nbr(ngrow_distance);
 debug_ngrow(MASK_NBR_MF,ngrow_distance,90);
 if (localMF[MASK_NBR_MF]->nComp()!=4)
  amrex::Error("invalid ncomp for mask nbr");
 VOF_Recon_resize(ngrow_distance,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,ngrow_distance,90);
 debug_ngrow(ORIGDIST_MF,ngrow_distance,90);
 if (localMF[ORIGDIST_MF]->nComp()!=nmat*(1+AMREX_SPACEDIM))
  amrex::Error("invalid ncomp for origdist");

 if (profile_dist==1)
  before_profile = ParallelDescriptor::second();

 new_localMF(DIST_TOUCH_MF,nmat,0,-1);
 localMF[DIST_TOUCH_MF]->setVal(0.0);

 MultiFab* dist_coarse_mf;
 MultiFab* dist_touch_coarse_mf;

 if ((level>0)&&(level<=finest_level)) {
  dist_coarse_mf=new MultiFab(grids,dmap,nmat*(AMREX_SPACEDIM+1),0,
	MFInfo().SetTag("dist_coarse_mf"),FArrayBoxFactory());
  int dcomp=0;
  int scomp=0;
  FillCoarsePatch(*dist_coarse_mf,dcomp,cur_time_slab,LS_Type,scomp,
        nmat*(AMREX_SPACEDIM+1));

  // idx,scomp,ncomp,index,scompBC_map
  // FillCoarsePatchGHOST is ultimately called.
  // dest_lstGHOST for State_Type defaults to pc_interp.
  // scompBC_map==0 corresponds to pc_interp and FORT_EXTRAPFILL
  dist_touch_coarse_mf=
   new MultiFab(grids,dmap,nmat,0,
    MFInfo().SetTag("dist_touch_coarse_mf"),FArrayBoxFactory());

  for (int i=0;i<nmat;i++) {
   Vector<int> scompBC_map;
   scompBC_map.resize(1);
   scompBC_map[0]=0;
   PCINTERP_fill_coarse_patch(DIST_TOUCH_MF,i,
     1,State_Type,scompBC_map);
  } // i=0..nmat-1
  MultiFab::Copy(*dist_touch_coarse_mf,*localMF[DIST_TOUCH_MF],0,0,nmat,0);
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
 int nface=nmat*AMREX_SPACEDIM*2*(1+AMREX_SPACEDIM); 

  // (nmat,nmat,2)  left material, right material, frac_pair+dist_pair
 int nface_dst=nmat*nmat*2;

 new_localMF(STENCIL_MF,nstar,ngrow_distance,-1);
 localMF[STENCIL_MF]->setVal(0.0);

 int do_face_decomp=0;
 int tessellate=0;
  // FACEINIT is in: MOF_REDIST_3D.F90
 makeFaceFrac(tessellate,ngrow_distance,FACEFRAC_MF,do_face_decomp);
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

  // FACEINITTEST is in MOF_REDIST_3D.F90
  // FACETEST_MF has nmat * sdim components.
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
   FORT_STENINIT( 
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
    &rzflag,
    xlo,dx,
    &cur_time_slab,
    &ngrow_distance,
    &nmat,&nstar);
 } // mfi
} // omp
 ns_reconcile_d_num(115);

 localMF[STENCIL_MF]->FillBoundary(geom.periodicity());

 if (profile_dist==1) {
  after_profile = ParallelDescriptor::second();
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "level= " << level << '\n';
   std::cout << "STENINIT time " << after_profile-before_profile << '\n';
  }
 }

 debug_ngrow(FACETEST_MF,ngrow_distance,90);
 if (localMF[FACETEST_MF]->nComp()!=nmat*AMREX_SPACEDIM)
  amrex::Error("localMF[FACETEST_MF]->nComp() invalid");

 debug_ngrow(FACEFRAC_MF,ngrow_distance,90);
 if (localMF[FACEFRAC_MF]->nComp()!=nface)
  amrex::Error("localMF[FACEFRAC_MF]->nComp() invalid");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  debug_ngrow(FACEFRAC_SOLVE_MM_MF+dir,ngrow_distance,722);
  if (localMF[FACEFRAC_SOLVE_MM_MF+dir]->nComp()!=nface_dst)
   amrex::Error("localMF[FACEFRAC_SOLVE_MM_MF+dir]->nComp()!=nface_dst");
 }

 debug_ngrow(STENCIL_MF,ngrow_distance,90);

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
   FArrayBox& facefracfab=(*localMF[FACEFRAC_MF])[mfi];

   FArrayBox& facepairXfab=(*localMF[FACEFRAC_SOLVE_MM_MF])[mfi];
   FArrayBox& facepairYfab=(*localMF[FACEFRAC_SOLVE_MM_MF+1])[mfi];
   FArrayBox& facepairZfab=
      (*localMF[FACEFRAC_SOLVE_MM_MF+AMREX_SPACEDIM-1])[mfi];

   FArrayBox& facetestfab=(*localMF[FACETEST_MF])[mfi];
   FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi];

   FArrayBox& touchfab=(*localMF[DIST_TOUCH_MF])[mfi];
   FArrayBox& crsetouch=(*dist_touch_coarse_mf)[mfi];
   FArrayBox& crsedist=(*dist_coarse_mf)[mfi];

   int vofcomp=scomp_mofvars;
   Vector<int> vofbc=getBCArray(State_Type,gridno,vofcomp,1);

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // in: MOF_REDIST_3D.F90
   FORT_LEVELSTRIP( 
    &keep_all_interfaces,
    &nprocessed[tid_current],
    minLS[tid_current].dataPtr(),
    maxLS[tid_current].dataPtr(),
    &max_problen,
    &level,
    &finest_level,
    truncate_volume_fractions.dataPtr(),
    latent_heat.dataPtr(),
    maskfab.dataPtr(),
    ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
    facepairXfab.dataPtr(),
    ARLIM(facepairXfab.loVect()),ARLIM(facepairXfab.hiVect()),
    facepairYfab.dataPtr(),
    ARLIM(facepairYfab.loVect()),ARLIM(facepairYfab.hiVect()),
    facepairZfab.dataPtr(),
    ARLIM(facepairZfab.loVect()),ARLIM(facepairZfab.hiVect()),
    facefracfab.dataPtr(),
    ARLIM(facefracfab.loVect()),ARLIM(facefracfab.hiVect()),
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
    &rzflag,
    xlo,dx,
    &cur_time_slab,
    &ngrow_distance,
    &nmat,&nten,
    &nstar,
    &nface,
    &nface_dst);
 } // mfi
} // omp

 ns_reconcile_d_num(116);

 for (int tid=1;tid<thread_class::nthreads;tid++) {
  for (int im=0;im<nmat;im++) {
   if (minLS[tid][im]<minLS[0][im])
    minLS[0][im]=minLS[tid][im];
   if (maxLS[tid][im]>maxLS[0][im])
    maxLS[0][im]=maxLS[tid][im];
  } // im=0..nmat-1
  nprocessed[0]+=nprocessed[tid];
 }
 ParallelDescriptor::ReduceIntSum(nprocessed[0]);
 for (int im=0;im<nmat;im++) {
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
   std::cout << "dist time " << after_dist-before_dist << '\n';
  }
 }

} // subroutine makeStateDist


void
NavierStokes::correct_dist_uninit() {

 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();

 const Real* dx = geom.CellSize();
 int nmat=num_materials;

 MultiFab& LS_new = get_new_data(LS_Type,slab_step+1);
 if (localMF[DIST_TOUCH_MF]->nComp()!=nmat)
  amrex::Error("localMF[DIST_TOUCH_MF]->nComp()!=nmat");

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

   FORT_CORRECT_UNINIT( 
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
    &cur_time_slab,
    &nmat);
 } // mfi
} // omp
 ns_reconcile_d_num(117);

} // subroutine correct_dist_uninit


// WARNING:  allocates, but does not delete.
// called from NavierStokes::make_physics_varsALL if using supermesh.
// called from NavierStokes::ColorSum
// called from NavierStokes::makeStateDist
void
NavierStokes::ProcessFaceFrac(int tessellate,int idxsrc,int idxdst,
		int ngrow_dest) {
  
 bool use_tiling=ns_tiling;

 int nmat=num_materials;
  // (nmat,sdim,2,sdim+1) area+centroid on each face of a cell.
 int nface_src=nmat*AMREX_SPACEDIM*2*(1+AMREX_SPACEDIM); 
  // (nmat,nmat,2)  left material, right material, frac_pair+dist_pair
 int nface_dst=nmat*nmat*2;

 int finest_level=parent->finestLevel();

 if ((tessellate!=0)&&(tessellate!=1))
  amrex::Error("tessellate invalid");

 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid ProcessFaceFrac");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  if (localMF_grow[idxdst+dir]>=0) {
   delete_localMF(idxdst+dir,1);
  }
  new_localMF(idxdst+dir,nface_dst,ngrow_dest,dir);
  localMF[idxdst+dir]->setVal(0.0);
 }

 int ngrow_source=ngrow_dest;
 if (ngrow_dest==0)
  ngrow_source=1;

 debug_ngrow(idxsrc,ngrow_source,90);
 if (localMF[idxsrc]->nComp()!=nface_src)
  amrex::Error("idxsrc has invalid ncomp");

 VOF_Recon_resize(ngrow_source,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,ngrow_source,90);

   // (1) =1 interior  =1 fine-fine ghost in domain  =0 otherwise
   // (2) =1 interior  =0 otherwise
   // (3) =1 interior+ngrow-1  =0 otherwise
   // (4) =1 interior+ngrow    =0 otherwise
 resize_mask_nbr(ngrow_source);
 debug_ngrow(MASK_NBR_MF,ngrow_source,90);

 int rzflag=0;
 if (geom.IsRZ())
  rzflag=1;
 else if (geom.IsCartesian())
  rzflag=0;
 else if (geom.IsCYLINDRICAL())
  rzflag=3;
 else
  amrex::Error("CoordSys bust 7");

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
   FORT_FACEPROCESS( 
    &ngrow_source,
    &ngrow_dest,
    &tid_current,
    &dir,
    &tessellate,
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
    &rzflag,
    xlo,dx,
    &cur_time_slab,
    &nmat,&nface_src,&nface_dst);
  } // mfi
} // omp
  ns_reconcile_d_num(118);

  localMF[idxdst+dir]->FillBoundary(geom.periodicity());
 } //dir=0..sdim-1

} // subroutine ProcessFaceFrac



// WARNING: makeFaceFrac allocates, but does not delete.
void
NavierStokes::makeFaceFrac(
 int tessellate,int ngrow,int idx,int do_face_decomp) {

 
 bool use_tiling=ns_tiling;

 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

  // (nmat,sdim,2,sdim+1) area+centroid on each face of a cell.
 int nface=nmat*AMREX_SPACEDIM*2*(1+AMREX_SPACEDIM); 

  // inside,outside,area+centroid,dir,side)
  // (nmat,nmat,sdim+1,sdim,2)
 int nface_decomp=0;
 if (do_face_decomp==1) {
  nface_decomp=nmat*nmat*(AMREX_SPACEDIM+1)*AMREX_SPACEDIM*2;
 } else if (do_face_decomp==0) {
  // do nothing
 } else
  amrex::Error("do_face_decomp invalid"); 

 int finest_level=parent->finestLevel();

 if (localMF_grow[idx]>=0) {
  delete_localMF(idx,1);
 }

 new_localMF(idx,nface+nface_decomp,ngrow,-1);
 localMF[idx]->setVal(0.0);
 VOF_Recon_resize(ngrow,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,ngrow,90);
 resize_mask_nbr(ngrow);
 debug_ngrow(MASK_NBR_MF,ngrow,90);

 int rzflag=0;
 if (geom.IsRZ())
  rzflag=1;
 else if (geom.IsCartesian())
  rzflag=0;
 else if (geom.IsCYLINDRICAL())
  rzflag=3;
 else
  amrex::Error("CoordSys bust 7");

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
   FORT_FACEINIT( 
    &tid_current,
    &tessellate,
    &nten,
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
    &rzflag,
    xlo,dx,
    &cur_time_slab,
    &ngrow,
    &nmat,
    &nface,
    &nface_decomp);
 } // mfi
} // omp
 ns_reconcile_d_num(119);

 localMF[idx]->FillBoundary(geom.periodicity());

} // subroutine makeFaceFrac


// WARNING: makeFaceTest allocates, but does not delete.
void
NavierStokes::makeFaceTest(int tessellate,int ngrow,int idx) {

 
 bool use_tiling=ns_tiling;

 int nmat=num_materials;
  // (im,dir,side,dir2)  (dir2==1 => area  dir2==2..sdim+1 => cen)
 int nface=nmat*AMREX_SPACEDIM*2*(1+AMREX_SPACEDIM); 

 int finest_level=parent->finestLevel();

 if ((tessellate!=0)&&(tessellate!=1))
  amrex::Error("tessellate invalid");

 if (localMF_grow[idx]>=0)
  amrex::Error("makeFaceTest: forgot to delete");

 new_localMF(idx,nmat*AMREX_SPACEDIM,ngrow,-1);
 localMF[idx]->setVal(0.0);

 VOF_Recon_resize(ngrow,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,ngrow,90);
 debug_ngrow(FACEFRAC_MF,ngrow,90);
 resize_mask_nbr(ngrow);
 debug_ngrow(MASK_NBR_MF,ngrow,90);

 if (localMF[FACEFRAC_MF]->nComp()!=nface)
  amrex::Error("localMF[FACEFRAC_MF]->nComp() invalid");

 int rzflag=0;
 if (geom.IsRZ())
  rzflag=1;
 else if (geom.IsCartesian())
  rzflag=0;
 else if (geom.IsCYLINDRICAL())
  rzflag=3;
 else
  amrex::Error("CoordSys bust 7");

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

   FORT_FACEINITTEST( 
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
    &rzflag,
    xlo,dx,
    &cur_time_slab,
    &ngrow,
    &nmat,
    &nface);
 } // mfi
} // omp
 ns_reconcile_d_num(120);

} // subroutine makeFaceTest


void
NavierStokes::makeDotMask(int nsolve,int project_option) {
 
 bool use_tiling=ns_tiling;

 int nmat=num_materials;

 int finest_level=parent->finestLevel();

 if (localMF_grow[DOTMASK_MF]>=0)
  amrex::Error("makeDotMask: forgot to delete");

 int num_materials_face=num_materials_vel;
 if ((project_option==0)||
     (project_option==1)||
     (project_option==10)||
     (project_option==11)||  // FSI_material_exists (2nd project)
     (project_option==13)||  // FSI_material_exists (1st project)
     (project_option==12)||  // pressure extension
     (project_option==3)) {  // viscosity
  if (num_materials_face!=1)
   amrex::Error("num_materials_face invalid");
 } else if ((project_option==2)||  // thermal diffusion
            ((project_option>=100)&&
             (project_option<100+num_species_var))) {
  num_materials_face=num_materials_scalar_solve;
 } else
  amrex::Error("project_option invalid9");

 new_localMF(DOTMASK_MF,num_materials_face,0,-1);
 setVal_localMF(DOTMASK_MF,1.0,0,num_materials_face,0);

 VOF_Recon_resize(1,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,1,90);
 resize_mask_nbr(1);
 debug_ngrow(MASK_NBR_MF,1,90);
 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,1,90);

 const Real* dx = geom.CellSize();

 if (num_materials_face==nmat) {

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(localMF[DOTMASK_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(*localMF[DOTMASK_MF],use_tiling); mfi.isValid(); ++mfi) {
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

     // mask=tag if not covered by level+1 or outside the domain.
   FArrayBox& maskfab=(*localMF[MASKCOEF_MF])[mfi];
   FArrayBox& dotmaskfab=(*localMF[DOTMASK_MF])[mfi];

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     // in: LEVELSET_3D.F90
   FORT_DOTMASK_BUILD( 
    &num_materials_face,
    &level,
    &finest_level,
    dotmaskfab.dataPtr(),
    ARLIM(dotmaskfab.loVect()),ARLIM(dotmaskfab.hiVect()),
    maskfab.dataPtr(),
    ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
    voffab.dataPtr(),
    ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    xlo,dx,
    &cur_time_slab,
    &nmat);
  } // mfi
} // omp
  ns_reconcile_d_num(121);

 } else if (num_materials_face==1) {
  // do nothing
 } else
  amrex::Error("num_materials_face invalid");

} // subroutine makeDotMask


// WARNING: makeCellFrac allocates, but does not delete.
void
NavierStokes::makeCellFrac(int tessellate,int ngrow,int idx) {
 
 bool use_tiling=ns_tiling;

 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;
  // (nmat,nmat,3+sdim)
  // im_inside,im_outside,3+sdim --> area, dist_to_line, dist, line normal.
 int ncellfrac=nmat*nmat*(3+AMREX_SPACEDIM); 
 int finest_level=parent->finestLevel();

 if (localMF_grow[idx]>=0) {
  delete_localMF(idx,1);
 }

 new_localMF(idx,ncellfrac,ngrow,-1);
 localMF[idx]->setVal(0.0);

 int ngrow_resize=ngrow;
 if (ngrow_resize==0)
  ngrow_resize++;

 VOF_Recon_resize(ngrow_resize,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,ngrow_resize,90);
 resize_mask_nbr(ngrow_resize);
 debug_ngrow(MASK_NBR_MF,ngrow_resize,90);

 int rzflag=0;
 if (geom.IsRZ())
  rzflag=1;
 else if (geom.IsCartesian())
  rzflag=0;
 else if (geom.IsCYLINDRICAL())
  rzflag=3;
 else
  amrex::Error("CoordSys bust 8");

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
   FORT_CELLFACEINIT( 
    &tid_current,
    &tessellate,
    &nten,
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
    &rzflag,
    xlo,dx,
    &cur_time_slab,
    &ngrow,
    &nmat,
    &ncellfrac);
 } // mfi
} // omp
 ns_reconcile_d_num(122);

 localMF[idx]->FillBoundary(geom.periodicity());

} // makeCellFrac


// called from make_physics_varsALL
void
NavierStokes::makeStateCurv(int project_option,int post_restart_flag) {
 
 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid makeStateCurv");

 int rzflag=0;
 if (geom.IsRZ())
  rzflag=1;
 else if (geom.IsCartesian())
  rzflag=0;
 else if (geom.IsCYLINDRICAL())
  rzflag=3;
 else
  amrex::Error("CoordSys bust 9");

 if (curv_stencil_height!=4)
  amrex::Error("curv_stencil_height invalid");

 if (ngrow_distance!=4)
  amrex::Error("ngrow_distance invalid");

 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

 Vector< Real > curv_min_local;
 Vector< Real > curv_max_local;
 curv_min_local.resize(thread_class::nthreads);
 curv_max_local.resize(thread_class::nthreads);

 for (int tid=0;tid<thread_class::nthreads;tid++) {
  curv_min_local[tid]=1.0e+99;
  curv_max_local[tid]=-1.0e+99;
 } // tid

  // height function curvature
  // finite difference curvature
  // pforce
  // marangoni force
  // dir/side flag
  // im3
  // x nten
 int num_curv=nten*(AMREX_SPACEDIM+5); 

 resize_metrics(1);

 resize_levelsetLO(ngrow_distance,LEVELPC_MF);

 if (localMF[LEVELPC_MF]->nComp()!=nmat*(1+AMREX_SPACEDIM))
  amrex::Error("localMF[LEVELPC_MF]->nComp() invalid");
 if (localMF[LEVELPC_MF]->nGrow()!=ngrow_distance)
  amrex::Error("localMF[LEVELPC_MF]->nGrow() invalid");

 int nhistory=localMF[HISTORY_MF]->nComp();
 if (nhistory==nten*2) {
  // do nothing
 } else
  amrex::Error("nhistory invalid");

 if (localMF_grow[DIST_CURV_MF]>=0) {
  delete_localMF(DIST_CURV_MF,1);
 }
 new_localMF(DIST_CURV_MF,num_curv,1,-1);
 localMF[DIST_CURV_MF]->setVal(0.0);

 getStateDist_localMF(GHOSTDIST_MF,ngrow_distance,cur_time_slab,16);
 debug_ngrow(GHOSTDIST_MF,ngrow_distance,906);
 if (localMF[GHOSTDIST_MF]->nComp()!=(1+AMREX_SPACEDIM)*nmat)
  amrex::Error("localMF[GHOSTDIST_MF]->nComp() invalid");

 if ((project_option==0)||
     (project_option==1)) {

  const Real* dx = geom.CellSize();

  Real cl_time=prev_time_slab;
  if (project_option==0)  // regular project
   cl_time=prev_time_slab;
  else if (project_option==1) // initial project
   cl_time=cur_time_slab;
  else if (project_option==10) // sync project
   cl_time=prev_time_slab;
  else
   amrex::Error("project_option invalid makeStateCurv");

  if (num_materials_vel!=1)
   amrex::Error("num_materials_vel invalid");

  MultiFab* CL_velocity=getState(2,0,
    num_materials_vel*(AMREX_SPACEDIM+1),cl_time);
  MultiFab* den=getStateDen(2,cl_time);
  if (den->nComp()!=nmat*num_state_material)
   amrex::Error("invalid ncomp for den");

   // mask=1 if not covered or if outside the domain.
   // NavierStokes::maskfiner_localMF
   // NavierStokes::maskfiner
  resize_maskfiner(1,MASKCOEF_MF);
  debug_ngrow(MASKCOEF_MF,1,28); 
  resize_mask_nbr(1);
  debug_ngrow(MASK_NBR_MF,1,2);

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
    int bfact_grid=parent->blockingFactor(level);

    const Real* xlo = grid_loc[gridno].lo();

    FArrayBox& histfab=(*localMF[HISTORY_MF])[mfi];

    // mask=tag if not covered by level+1 or outside the domain.
    FArrayBox& maskcov=(*localMF[MASKCOEF_MF])[mfi];

    FArrayBox& lsfab=(*localMF[LEVELPC_MF])[mfi];
    FArrayBox& lshofab=(*localMF[GHOSTDIST_MF])[mfi];

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

     // in: LEVELSET_3D.F90
    FORT_CURVSTRIP(
     &post_restart_flag,
     &conservative_tension_force,
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
     lshofab.dataPtr(),
     ARLIM(lshofab.loVect()),ARLIM(lshofab.hiVect()),
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
     &rzflag,
     xlo,dx,
     &cur_time_slab,
     &visc_coef,
     &nmat,&nten,
     &num_curv,
     &ngrow_distance,
     &curv_stencil_height);
  } // mfi
} //omp
  ns_reconcile_d_num(123);

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
    std::cout << "c++ ngrow,csten " << ngrow_distance << ' ' <<
     curv_stencil_height << ' ' << '\n';

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
     tecplot_debug(lsfab,xlo,fablo,fabhi,dx,-1,0,0,nmat,interior_only);

     std::cout << "output of lshofab (GHOSTDIST)" << '\n';
     FArrayBox& lshofab=(*localMF[GHOSTDIST_MF])[mfi];
     tecplot_debug(lshofab,xlo,fablo,fabhi,dx,-1,0,0,nmat,interior_only);

    } // mfi
    ns_reconcile_d_num(124);

  } // ((fab_verbose==1)||(fab_verbose==3))

  delete CL_velocity;
  delete den;

 } else if (project_option==10) {

   // do nothing

 } else
   amrex::Error("project_option invalid10");

 delete_localMF(GHOSTDIST_MF,1);

}  // subroutine makeStateCurv


MultiFab* NavierStokes::getStateMAC(int ngrow,int dir,
 int scomp,int ncomp,Real time) {

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 int nsolve=1;
 int nsolveMM_FACE=nsolve*num_materials_vel;
 
 if ((dir<0)||(dir>=AMREX_SPACEDIM))
  amrex::Error("dir invalid get state mac");

 MultiFab& S_new=get_new_data(Umac_Type+dir,slab_step+1);
 int ntotal=S_new.nComp();
 if (ntotal!=nsolveMM_FACE)
  amrex::Error("ntotal bust");
 if (scomp+ncomp>ntotal)
  amrex::Error("scomp invalid getStateMAC");

 MultiFab* mf = new MultiFab(state[Umac_Type+dir].boxArray(),dmap,ncomp,
   ngrow,MFInfo().SetTag("mf getStateMAC"),FArrayBoxFactory());

 FillPatch(*this,*mf,0,time,Umac_Type+dir,scomp,ncomp);

 ParallelDescriptor::Barrier();

 return mf;

}  // subroutine getStateMAC



void
NavierStokes::ctml_fsi_transfer_force() {
	
 if (ParallelDescriptor::IOProcessor() && verbose)
  std::cout << "in NavierStokes::ctml_fsi_transfer_force() \n";

 bool use_tiling=ns_tiling;

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

 int nmat=num_materials;

  // nparts x (velocity + LS + temperature + flag + stress)
 int nparts=im_solid_map.size();
 if ((nparts<1)||(nparts>nmat))
  amrex::Error("nparts invalid");

 if (nFSI_sub!=12)
  amrex::Error("nFSI_sub invalid");
 if (ngrowFSI!=3)
  amrex::Error("ngrowFSI invalid");
 int nFSI=nparts*nFSI_sub;

 debug_ngrow(FSI_MF,0,1);
 if (localMF[FSI_MF]->nComp()!=nFSI)
  amrex::Error("localMF[FSI_MF]->nComp() invalid");

 for (int partid=0;partid<nparts;partid++) {

  int im_part=im_solid_map[partid];

  if (ns_is_rigid(im_part)==1) {

   if (CTML_FSI_matC(im_part)==1) {

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

#ifdef MVAHABFSI
     const int gridno = mfi.index();
     const Box& fabgrid = grids[gridno];
     const int* tilelo=tilegrid.loVect();
     const int* tilehi=tilegrid.hiVect();
     const int* fablo=fabgrid.loVect();
     const int* fabhi=fabgrid.hiVect();

     FArrayBox& snewfab=S_new[mfi];
     FArrayBox& forcefab=(*localMF[FSI_MF])[mfi];
      // nparts x (velocity + LS + temperature + flag + stress)
     int ibase=partid*nFSI_sub+6;

     FORT_CTMLTRANSFERFORCE(
      tilelo, tilehi, 
      fablo, fabhi, 
      snewfab.dataPtr(), 
      ARLIM(snewfab.loVect()), ARLIM(snewfab.hiVect()), 
      forcefab.dataPtr(ibase), 
      ARLIM(forcefab.loVect()), ARLIM(forcefab.hiVect()));
#else
     amrex::Error("CTML(C): define MEHDI_VAHAB_FSI in GNUmakefile");
#endif
    }  // mfi  
} // omp
    ns_reconcile_d_num(125);

   } else if (CTML_FSI_matC(im_part)==0) {
    // do nothing
   } else 
    amrex::Error("CTML_FSI_matC(im_part) invalid");

  } else
   amrex::Error("ns_is_rigid(im_part) invalid");

 } // partid=0 ... nparts-1

} // subroutine ctml_fsi_transfer_force()

}/* namespace amrex */

