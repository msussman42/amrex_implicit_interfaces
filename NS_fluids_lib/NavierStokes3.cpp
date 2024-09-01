//#include <winstd.H>

#include <algorithm>
#include <vector>

#include <cstdio>
#include <cmath>

#include <AMReX_CoordSys.H>
#include <AMReX_Geometry.H>
#include <AMReX_BoxDomain.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_Utility.H>
#include <NavierStokes.H>
#include <GLOBALUTIL_F.H>
#include <GODUNOV_F.H>
#include <NAVIERSTOKES_F.H>
#include <MACOPERATOR_F.H>
#include <PROB_F.H>
#include <LEVEL_F.H>
#include <SOLIDFLUID_F.H>
#include <DERIVE_F.H>
#include <DIFFUSION_F.H>

namespace amrex{

#define profile_solver 0

// status==1 success
// status==0 failure
extern void matrix_solveCPP(Real** AA,Real* xx,Real* bb,
 int& status,int numelem);

extern void set_x_vel_bc_NS_setup(BCRec& bc,const BCRec& phys_bc);
extern void set_y_vel_bc_NS_setup(BCRec& bc,const BCRec& phys_bc);
extern void set_z_vel_bc_NS_setup(BCRec& bc,const BCRec& phys_bc);

// if ncomp_input==-1, then ncomp=S_crse.ncomp()
// spectral_override==LOW_ORDER_AVGDOWN => always do low order average down.
// grid_type=-1,..,5
void
NavierStokes::avgDownEdge(
 int grid_type,MultiFab& S_crse,MultiFab& S_fine,
 int scomp,int ncomp_input,int spectral_override,
 const std::string& caller_string) {

 if (1==0) {
  std::cout << "avgDownEdge caller_string= " << caller_string << 
   " grid_type= " << grid_type << '\n';
 }
 std::string local_caller_string="avgDownEdge";
 local_caller_string=caller_string+local_caller_string;

 if ((grid_type>=-1)&&(grid_type<=5)) {
  // do nothing
 } else 
  amrex::Error("grid_type invalid in avgdown edge");

 int finest_level=parent->finestLevel();
 if (level>=finest_level) 
  amrex::Error("level invalid in avgDownEdge");

 int f_level=level+1;
 NavierStokes& fine_lev = getLevel(f_level);
 resize_metrics(1);
 debug_ngrow(VOLUME_MF,0,local_caller_string);
 debug_ixType(VOLUME_MF,-1,local_caller_string);
 fine_lev.resize_metrics(1);
 fine_lev.debug_ngrow(VOLUME_MF,0,local_caller_string);
 fine_lev.debug_ixType(VOLUME_MF,-1,local_caller_string);

 debug_ixType_raw(&S_fine,grid_type,local_caller_string);
 debug_ixType_raw(&S_crse,grid_type,local_caller_string);

 fine_lev.debug_boxArray(&S_fine,grid_type,local_caller_string);
 debug_boxArray(&S_crse,grid_type,local_caller_string);

 int ncomp=S_crse.nComp();
 if ((S_crse.nComp()!=ncomp)||(S_fine.nComp()!=ncomp)) {
  std::cout << "coarse ncomp " << S_crse.nComp() << '\n';
  std::cout << "fine ncomp " << S_fine.nComp() << '\n';
  amrex::Error("ncomp invalid in avgdownedge");
 }
 if ((scomp<0)||(scomp>=ncomp))
  amrex::Error("scomp invalid avgDownEdge");

 if (ncomp_input==-1) {
  // do nothing - ncomp already has the correct value.
 } else if ((ncomp_input>=1)&&(ncomp_input<=ncomp)) {
  ncomp=ncomp_input;
 } else
  amrex::Error("ncomp_input invalid");

 const BoxArray& fgridsMAC=S_fine.boxArray();
 const BoxArray& fgridscen=fine_lev.grids;
 const DistributionMapping& fdmap=fine_lev.dmap;

 if (fgridsMAC.size()==fgridscen.size()) {
  // do nothing
 } else {
  amrex::Error("expecting: fgridsMAC.size()==fgridscen.size()");
 }

 BoxArray crse_S_fine_BA_MAC(fgridsMAC.size());
 BoxArray crse_cen_fine_BA(fgridscen.size());

 for (int i = 0; i < fgridsMAC.size(); ++i) {
  crse_S_fine_BA_MAC.set(i,amrex::coarsen(fgridsMAC[i],2));
  Box cbox=amrex::coarsen(fgridscen[i],2);
  crse_cen_fine_BA.set(i,cbox);
 }

 DistributionMapping crse_dmap=fdmap;
 MultiFab crse_S_fine_MAC(crse_S_fine_BA_MAC,crse_dmap,ncomp,0,
   MFInfo().SetTag("crse_S_fine_MAC"),FArrayBoxFactory());
 debug_ixType_raw(&crse_S_fine_MAC,grid_type,local_caller_string);

 ParallelDescriptor::Barrier();

 int bfact=parent->Space_blockingFactor(level);
 int bfact_f=parent->Space_blockingFactor(f_level);
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
  BL_ASSERT(fgridsMAC[mfi.index()] == mfi.validbox());
  const Box& tilegrid = mfi.tilebox();

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  const int i = mfi.index();
  const Real* xlo_fine = fine_lev.grid_loc[i].lo();

  const Box& cbox = crse_cen_fine_BA[i];
  FArrayBox& crse_fab = crse_S_fine_MAC[mfi];

  const Box& fbox = fgridscen[i];
  const FArrayBox& fine_fab = S_fine[mfi];

  FArrayBox& maskfab=(*fine_lev.localMF[MASKSEM_MF])[mfi];

   // declared in: NAVIERSTOKES_3D.F90
  fort_edgeavgdown(
   &enable_spectral,
   &finest_level,
   &spectral_override,
   prob_lo,
   dxf,
   &level,&f_level,
   &bfact,&bfact_f,
   xlo_fine,dx,
   &grid_type,
   crse_fab.dataPtr(),
   ARLIM(crse_fab.loVect()),ARLIM(crse_fab.hiVect()),
   fine_fab.dataPtr(scomp),
   ARLIM(fine_fab.loVect()),ARLIM(fine_fab.hiVect()),
   maskfab.dataPtr(),
   ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
   cbox.loVect(),cbox.hiVect(),
   fbox.loVect(),fbox.hiVect(),
   &ncomp);
 }  // mfi
} //omp
 ns_reconcile_d_num(LOOP_EDGEAVGDOWN,"avgDownEdge");

 S_crse.ParallelCopy(crse_S_fine_MAC,0,scomp,ncomp);
 ParallelDescriptor::Barrier();

 const Box& domain = geom.Domain();

 int box_type[AMREX_SPACEDIM]; 
 grid_type_to_box_type_cpp(grid_type,box_type);

 for (int local_dir=0;local_dir<AMREX_SPACEDIM;local_dir++) {
  if (box_type[local_dir]==0) {
   // do nothing
  } else if (box_type[local_dir]==1) {
   if (geom.isPeriodic(local_dir)) {
    IntVect pshift=IntVect::TheZeroVector();
    pshift[local_dir]=domain.length(local_dir);
    crse_S_fine_MAC.shift(pshift);

    ParallelDescriptor::Barrier();
    S_crse.ParallelCopy(crse_S_fine_MAC,0,scomp,ncomp);
    ParallelDescriptor::Barrier();

    pshift[local_dir]=-2*domain.length(local_dir);
    crse_S_fine_MAC.shift(pshift);

    S_crse.ParallelCopy(crse_S_fine_MAC,0,scomp,ncomp);
    ParallelDescriptor::Barrier();
   } else if (!geom.isPeriodic(local_dir)) {
    // do nothing
   } else {
    amrex::Error("geom.isPeriodic(local_dir) bust");
   } 
  } else
   amrex::Error("box_type bust");
 } // local_dir=0 ... sdim-1

}  // end subroutine avgDownEdge


// called from updatevelALL,multiphase_project
// interpolate from level+1 to level.
void
NavierStokes::avgDownMac() {

 std::string local_caller_string="avgDownMac";

 int ncomp_edge=-1;
 int scomp=0;
  // spectral_override==0 => always do low order average down.
 int spectral_override=SPECTRAL_ORDER_AVGDOWN;
 avgDownEdge_localMF(UMAC_MF,scomp,ncomp_edge,0,AMREX_SPACEDIM, 
   spectral_override,local_caller_string);

}  // end subroutine avgDownMac

// spectral_override==0 => always do low order average down.
void NavierStokes::avgDownMacState(int spectral_override) {

 std::string local_caller_string="avgDownMacState";

 int finest_level = parent->finestLevel();

 if (level>=finest_level) 
  amrex::Error("level invalid avgDownMacState");

 NavierStokes& fine_lev = getLevel(level+1);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  MultiFab& S_crse = get_new_data(Umac_Type+dir,slab_step+1);
  MultiFab& S_fine = fine_lev.get_new_data(Umac_Type+dir,slab_step+1);
  int scomp=0;
  int ncomp_edge=S_crse.nComp(); 
 
  if (ncomp_edge==1) {
   // do nothing
  } else
   amrex::Error("ncomp_edge invalid in avgDownMacState");
   
  if ((S_crse.nComp()!=ncomp_edge)||
      (S_fine.nComp()!=ncomp_edge))
   amrex::Error("S_crse.nComp() or S_fine.nComp() invalid");

  avgDownEdge(dir,S_crse,S_fine,scomp,ncomp_edge,spectral_override,
    local_caller_string);
 }  // dir=0..sdim-1

}  // end subroutine avgDownMacState

void NavierStokes::init_rest_fraction(const std::string& caller_string) {

 if (rest_fraction==0.0) {
  //do nothing
 } else if ((rest_fraction>0.0)&&(rest_fraction<=4.0)) {
  ParallelDescriptor::Barrier();
  std::fflush(NULL);
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "init_rest_fraction caller_string= " << caller_string << '\n';
  }
  std::fflush(NULL);
  ParallelDescriptor::Barrier();
  start_rest=ParallelDescriptor::second(); 
 } else
  amrex::Error("rest_fraction invalid");

}

void NavierStokes::finalize_rest_fraction(const std::string& caller_string) {

 if (rest_fraction==0.0) {
  //do nothing
 } else if ((rest_fraction>0.0)&&(rest_fraction<=4.0)) {
  ParallelDescriptor::Barrier();
  std::fflush(NULL);
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "finalize_rest_fraction caller_string= " << 
     caller_string << '\n';
  }
  std::fflush(NULL);
  ParallelDescriptor::Barrier();
  end_rest=ParallelDescriptor::second();
  double len_rest=end_rest-start_rest;
  if (len_rest>0.0) {
   double abs_rest=len_rest*rest_fraction;
   if (ParallelDescriptor::IOProcessor()) {
    std::cout << "end_rest-start_rest= " << len_rest << '\n';
   }
   amrex::Sleep(abs_rest);
   ParallelDescriptor::Barrier();
   if (ParallelDescriptor::IOProcessor()) {
    std::cout << "end_rest-start_rest= " << len_rest << '\n';
    std::cout << "abs_rest= " << abs_rest << '\n';
    std::cout << "done sleeping\n";
   }
   std::fflush(NULL);
  } else
   amrex::Error("len_rest invalid");

 } else
  amrex::Error("rest_fraction invalid");

}

void NavierStokes::nonlinear_advection(const std::string& caller_string) {

 std::string local_caller_string="nonlinear_advection";
 local_caller_string=caller_string+local_caller_string;

 init_rest_fraction(local_caller_string);

 if (pattern_test(local_caller_string,"do_the_advance")==1) {
  advect_time_slab=prev_time_slab;

  if (divu_outer_sweeps==0) 
   vel_time_slab=prev_time_slab;
  else if (divu_outer_sweeps>0)
   vel_time_slab=cur_time_slab;
  else
   amrex::Error("divu_outer_sweeps invalid nonlinear_advection");

 } else
  amrex::Error("caller is invalid in nonlinear_advection");

 int renormalize_only=1;

 if (level!=0)
  amrex::Error("level invalid nonlinear_advection");

 if (cur_time_slab>0.0) {
  // do nothing
 } else
  amrex::Error("cur_time_slab invalid");

 if ((slab_step>=0)&&(slab_step<ns_time_order)) {
  //do nothing
 } else
  amrex::Error("slab_step invalid");

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid nonlinear_advection");

 if ((divu_outer_sweeps>=0)&&
     (divu_outer_sweeps<num_divu_outer_sweeps)) {
  // do nothing
 } else
  amrex::Error("divu_outer_sweeps invalid nonlinear_advection");

 int finest_level=parent->finestLevel();

 if (std::abs(cur_time_slab-prev_time_slab-dt_slab)>CPP_EPS_5_3) {
  std::cout << "cur_time_slab " << cur_time_slab << '\n';
  std::cout << "prev_time_slab " << prev_time_slab << '\n';
  std::cout << "dt_slab " << dt_slab << '\n';
  amrex::Error("slab time bust1");
 }

 if (verbose>0)
  if (ParallelDescriptor::IOProcessor())
   std::cout << "nonlinear advect \n";

#if (NS_profile_solver==1)
 BLProfiler bprof(local_caller_string);
#endif

 NavierStokes& ns_fine=getLevel(finest_level);
 int basestep=ns_fine.nStep();

 order_direct_split=basestep-2*(basestep/2);

 if ((order_direct_split!=0)&&
     (order_direct_split!=1))
  amrex::Error("order_direct_split invalid");

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.new_localMF(VOF_PREV_TIME_MF,num_materials,2,-1);
 }  // ilev=finest_level ... level

  // order_direct_split=base_step mod 2
 map_forward_direct_split.resize(AMREX_SPACEDIM);
 normdir_direct_split.resize(AMREX_SPACEDIM);

 int map_forward_start=0;
 if (AMREX_SPACEDIM==2) {
  map_forward_start=0;
 } else if (AMREX_SPACEDIM==3) {
  if (order_direct_split==0)
   map_forward_start=0;
  else if (order_direct_split==1)
   map_forward_start=1;
  else
   amrex::Error("order_direct_split invalid");
 } else
  amrex::Error("dimension bust");
   
 for (dir_absolute_direct_split=0;
      dir_absolute_direct_split<AMREX_SPACEDIM;
      dir_absolute_direct_split++) {

  if (order_direct_split==0)
   normdir_direct_split[dir_absolute_direct_split]=dir_absolute_direct_split;
  else if (order_direct_split==1)
   normdir_direct_split[dir_absolute_direct_split]=
     AMREX_SPACEDIM-(dir_absolute_direct_split+1);
  else
   amrex::Error("order_direct_split invalid");

  int normdir_here=normdir_direct_split[dir_absolute_direct_split];

  if ((EILE_flag==-1)||  // Weymouth and Yue
      (EILE_flag==2)) {  // always EI
   map_forward_direct_split[normdir_here]=0;
  } else if (EILE_flag==3) { // always LE
   map_forward_direct_split[normdir_here]=1;
  } else if (EILE_flag==1) { // EI-LE 2d, EI-LE-EI or LE-EI-LE 3d
   map_forward_direct_split[normdir_here]=map_forward_start;
   map_forward_start=1-map_forward_start;
  } else
   amrex::Error("EILE_flag invalid");
 } // dir_absolute_direct_split=0..sdim-1

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "order_direct_split " << order_direct_split << '\n';
  }
 }

 dir_absolute_direct_split=0;
 int init_vof_prev_time=1;
 int normdir_here=normdir_direct_split[dir_absolute_direct_split];
 if ((normdir_here<0)||(normdir_here>=AMREX_SPACEDIM))
  amrex::Error("normdir_here invalid (prior to loop)");

 // delete_advect_vars() called in NavierStokes::do_the_advance
 // delete_transport_vars() called in NavierStokes::do_the_advance
 // right after increment_face_velocityALL. 
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  // initialize ADVECT_REGISTER_FACE_MF and ADVECT_REGISTER_MF
  ns_level.prepare_advect_vars(prev_time_slab);
  // initialize TRANSPORT_REGISTER_FACE_MF
  ns_level.prepare_transport_vars(vel_time_slab);
 }

#ifdef AMREX_PARTICLES

 if ((slab_step>=0)&&(slab_step<ns_time_order)) {

  int lev_min=0;
  int lev_max=-1;
  int nGrow_Redistribute=0;
  int local_Redistribute=0; 
  bool local_copy=true; //do not redistribute inside of copyParticles
  bool remove_negative=true;

    // level=0
  My_ParticleContainer& prevPC=newDataPC(slab_step);

  prevPC.Redistribute(lev_min,lev_max,
    nGrow_Redistribute,local_Redistribute,remove_negative);

    // level=0
  My_ParticleContainer& localPC=newDataPC(slab_step+1);

  localPC.clearParticles();
  localPC.Redistribute();
  localPC.copyParticles(prevPC,local_copy);

    //prior to advection
  init_particle_containerALL(OP_PARTICLE_ADD,local_caller_string);

  localPC.Redistribute(lev_min,lev_max,
    nGrow_Redistribute,local_Redistribute,remove_negative);

 } else
  amrex::Error("slab_step invalid");

#endif

 interface_touch_flag=1; //nonlinear_advection

 //output:SLOPE_RECON_MF
 VOF_Recon_ALL( 
  local_caller_string, //nonlinear_advection
  advect_time_slab,
  RECON_UPDATE_NULL,
  init_vof_prev_time);

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.prepare_displacement();
 } // ilev

 for (dir_absolute_direct_split=0;
      dir_absolute_direct_split<AMREX_SPACEDIM;
      dir_absolute_direct_split++) {

   normdir_here=normdir_direct_split[dir_absolute_direct_split];
   if ((normdir_here<0)||(normdir_here>=AMREX_SPACEDIM))
    amrex::Error("normdir_here invalid (in loop)");

   init_vof_prev_time=0;

   if (dir_absolute_direct_split==0) {

    //do nothing

   } else if ((dir_absolute_direct_split>0)&&
              (dir_absolute_direct_split<AMREX_SPACEDIM)) {

    advect_time_slab=cur_time_slab;

    interface_touch_flag=1; //nonlinear_advection

     //output::SLOPE_RECON_MF
    VOF_Recon_ALL(
      local_caller_string, //nonlinear_advection
      advect_time_slab,
      RECON_UPDATE_STATE_CENTROID,
     init_vof_prev_time);

   } else
    amrex::Error("dir_absolute_direct_split invalid");

   split_scalar_advectionALL();

   interface_touch_flag=1; //nonlinear_advection

#ifdef AMREX_PARTICLES

   if ((slab_step>=0)&&(slab_step<ns_time_order)) {

    int lev_min=0;
    int lev_max=-1;
    int nGrow_Redistribute=0;
    int local_Redistribute=0; 
    bool remove_negative=true;

      //level==0
    My_ParticleContainer& localPC=newDataPC(slab_step+1);

    for (int ilev=finest_level;ilev>=level;ilev--) {
     NavierStokes& ns_level=getLevel(ilev);
      //move_particles() declared in NavierStokes2.cpp
     ns_level.move_particles(
       normdir_here,
       localPC,
       local_caller_string);
    }

    localPC.Redistribute(lev_min,lev_max,
      nGrow_Redistribute,local_Redistribute,remove_negative);

   } else
    amrex::Error("slab_step invalid");

#endif

   if ((dir_absolute_direct_split>=0)&&
       (dir_absolute_direct_split<AMREX_SPACEDIM-1)) {
     // in: nonlinear_advection
     // calls MOFavgDown, LS_Type avgDown
     // projects volume fractions so that sum F_m_fluid=1.
    renormalize_only=1;
    int local_truncate=0;
    int update_particles=0;
    prescribe_solid_geometryALL(prev_time_slab,renormalize_only,
      local_truncate,local_caller_string,update_particles);

     // velocity and pressure
    avgDownALL(State_Type,STATECOMP_VEL,
       STATE_NCOMP_VEL+STATE_NCOMP_PRES,SPECTRAL_ORDER_AVGDOWN);
     // "state" (all materials)
    avgDownALL(State_Type,STATECOMP_STATES,num_state_material*num_materials,
       SPECTRAL_ORDER_AVGDOWN);

    if ((num_materials_viscoelastic>=1)&&
        (num_materials_viscoelastic<=num_materials)) {
     avgDownALL_TENSOR();
    } else if (num_materials_viscoelastic==0) {
     // do nothing
    } else
     amrex::Error("num_materials_viscoelastic invalid in nonlinear_advection");

    if ((num_materials_compressible>=1)&&
        (num_materials_compressible<=num_materials)) {
     avgDownALL_refine_density();
    } else if (num_materials_compressible==0) {
     // do nothing
    } else
     amrex::Error("num_materials_compressible invalid in nonlinear_advection");

    for (int ilev=finest_level-1;ilev>=level;ilev--) {
     NavierStokes& ns_level=getLevel(ilev);
     ns_level.avgDownMacState(LOW_ORDER_AVGDOWN);
    }

   } else if (dir_absolute_direct_split==AMREX_SPACEDIM-1) {
     // do nothing
   } else
    amrex::Error("parameter bust");

 }  // dir_absolute_direct_split=0..sdim-1

 for (int ilev=level;ilev<=finest_level;ilev++) {
   NavierStokes& ns_level=getLevel(ilev);
   ns_level.delete_localMF(VOF_PREV_TIME_MF,1);

   ns_level.delete_localMF(MAC_VELOCITY_MF,AMREX_SPACEDIM);
   ns_level.delete_localMF(RAW_MAC_VELOCITY_MF,AMREX_SPACEDIM);

    // sanity check
   ns_level.debug_ngrow(MASKCOEF_MF,1,local_caller_string);
 } // ilev=level..finest_level

// 0. if read_from_CAD==1
//   (a) copy eulerian velocity and stress to lagrangian.
//   (b) update Lagrangian node positions
//   (c) convert Lagrangian position, velocity, temperature, and
//       force (if CTML) to Eulerian.
// 1. renormalize variables
// 2. extend from F>0 fluid regions into F=0 regions
// 3. if renormalize_only==0, 
//    a. init F,X,LS for the solid materials.
//    b. init U,T in the solid regions.
//    c. extrapolate F,X,LS from fluid regions into solid regions.

 if (read_from_CAD()==1) {

    //init_FSI_GHOST_MAC_MF_ALL is declared in NavierStokes.cpp
   renormalize_only=1;
   init_FSI_GHOST_MAC_MF_ALL(renormalize_only,local_caller_string);

   int fast_mode=0;
   setup_integrated_quantities();
   volWgtSumALL(local_caller_string,fast_mode);

   if (ok_copy_FSI_old_to_new()==1) { //stationary solid

    copy_old_FSI_to_new();

   } else if (ok_copy_FSI_old_to_new()==0) { //deforming solid

    int iter=0;

     //CTML solid must be wholly contained on the finest level.
    NavierStokes& ns_finest=getLevel(finest_level);

    ns_finest.resize_mask_nbr(ngrow_make_distance);

    ns_finest.ns_header_msg_level(
       OP_FSI_LAG_STRESS,
       SUB_OP_FSI_CLEAR_LAG_DATA,
       cur_time_slab,
       dt_slab,
       iter,
       local_caller_string);

    ns_finest.ns_header_msg_level(
       OP_FSI_LAG_STRESS,
       SUB_OP_FSI_COPY_TO_LAG_DATA,
       cur_time_slab,
       dt_slab,
       iter,
       local_caller_string);

    ns_finest.ns_header_msg_level(
       OP_FSI_LAG_STRESS,
       SUB_OP_FSI_SYNC_LAG_DATA,
       cur_time_slab,
       dt_slab,
       iter,
       local_caller_string);

    if (level==0) {
     // fort_headermsg (SOLIDFLUID.F90)
     // CLSVOF_ReadNodes (sci_clsvof.F90)
     // if FSI_flag==FSI_SHOELE_CTML then
     //  tick_fib is called (in ../StructureCodeShoele/tick.F)
     ns_header_msg_level(
      OP_FSI_UPDATE_NODES,
      SUB_OP_FSI_DEFAULT,
      cur_time_slab,
      dt_slab,
      iter,
      local_caller_string);
    } else
     amrex::Error("expecting level==0");

    // convert Lagrangian position, velocity, temperature, and force to
    // Eulerian.
    // go from coarsest to finest.
    for (int ilev=level;ilev<=finest_level;ilev++) {
     NavierStokes& ns_level=getLevel(ilev);
     ns_level.FSI_make_distance(cur_time_slab,dt_slab);
    } // ilev

    interface_touch_flag=1; //nonlinear_advection

    for (int ilev=level;ilev<=finest_level;ilev++) {
     NavierStokes& ns_level=getLevel(ilev);
     ns_level.resize_FSI_MF();
    }

   } else
    amrex::Error("ok_copy_FSI_old_to_new invalid");

 } else if (read_from_CAD()==0) {
  // do nothing
 } else
  amrex::Error("read_from_CAD() invalid");

 regenerate_from_eulerian(cur_time_slab);

 interface_touch_flag=1; //nonlinear_advection

 if (1==0) {
   // S_new is level 0 data
   MultiFab& S_new=get_new_data(State_Type,slab_step+1);
   // data file name "BEFOREPRESCRIBE<stuff>.plt"
   // xvel,yvel,zvel,pressure,(density, temperature) x num_materials,
   // (VFRAC,centroid) x num_materials, error indicator
   writeSanityCheckData(
    "BEFOREPRESCRIBE",
    "in: NavierStokes::nonlinear_advection, State_Type ", 
    local_caller_string,
    State_Type+GET_NEW_DATA_OFFSET, //tower_mf_id
    S_new.nComp(),
    -1, // data_mf==-1
    State_Type, //state_type_mf==State_Type
    -1, // data_dir==-1
    parent->levelSteps(0)); 
 }

  // in: nonlinear_advection
  // level set function, volume fractions, and centroids are
  // made "consistent" amongst the levels.
  // prescribe_solid_geometryALL is declared in: NavierStokes2.cpp
 renormalize_only=0;
 int local_truncate=1;
 int update_particles=1;
 prescribe_solid_geometryALL(cur_time_slab,renormalize_only,
   local_truncate,local_caller_string,update_particles);

 interface_touch_flag=1; //nonlinear_advection

 avgDownALL(State_Type,STATECOMP_VEL,STATE_NCOMP_VEL+STATE_NCOMP_PRES,1);

 if (1==0) {
    // S_new is level 0 data
  MultiFab& S_new=get_new_data(State_Type,slab_step+1);
   // data file name "AFTERPRESCRIBE<stuff>.plt"
   // xvel,yvel,zvel,pressure,(density, temperature) x num_materials,
   // (VFRAC,centroid) x num_materials, error indicator
  writeSanityCheckData(
   "AFTERPRESCRIBE",
   "in: NavierStokes::nonlinear_advection, State_Type ", 
   local_caller_string,
   State_Type+GET_NEW_DATA_OFFSET, //tower_mf_id
   S_new.nComp(),
   -1, // data_mf==-1
   State_Type, //state_type_mf==State_Type
   -1, // data_dir==-1
   parent->levelSteps(0)); 
 }
#if (NS_profile_solver==1)
 bprof.stop();
#endif

 finalize_rest_fraction(local_caller_string);

}  // end subroutine nonlinear_advection


void NavierStokes::allocate_SDC() {

 if ((ns_time_order==1)&&(enable_spectral!=0))
  amrex::Error("(ns_time_order==1)&&(enable_spectral!=0)");
 if ((ns_time_order>=2)&&(enable_spectral!=1))
  amrex::Error("(ns_time_order>=2)&&(enable_spectral!=1)");

 if ((ns_time_order>=2)&&(enable_spectral==1)) {

  if (localMF_grow[stableF_MF]==-1) {
   new_localMF(stableF_MF,NSTATE_SDC*ns_time_order,0,-1);
  } else
   amrex::Error("localMF_grow[stableF_MF] invalid");

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
   new_localMF(stableF_GP_MF+dir,ns_time_order,0,dir);

  if (localMF_grow[spectralF_MF]==-1) {
   new_localMF(spectralF_MF,NSTATE_SDC*(ns_time_order+1),0,-1);
  } else
   amrex::Error("localMF_grow[spectralF_MF] invalid");

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
   new_localMF(spectralF_GP_MF+dir,(ns_time_order+1),0,dir);

  if (localMF_grow[delta_MF]==-1) {
   new_localMF(delta_MF,NSTATE_SDC*ns_time_order,0,-1);
  } else
   amrex::Error("localMF_grow[delta_MF] invalid");

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
   new_localMF(delta_GP_MF+dir,ns_time_order,0,dir);

 } else if ((ns_time_order==1)&&(enable_spectral==0)) {
  // do nothing
 } else
  amrex::Error("ns_time_order or enable_spectral invalid");

}  // subroutine allocate_SDC


void NavierStokes::deallocate_SDC() {

 if ((ns_time_order==1)&&(enable_spectral!=0))
  amrex::Error("(ns_time_order==1)&&(enable_spectral!=0)");
 if ((ns_time_order>=2)&&(enable_spectral!=1))
  amrex::Error("(ns_time_order>=2)&&(enable_spectral!=1)");

 if ((ns_time_order>=2)&&(enable_spectral==1)) {

  if (localMF_grow[spectralF_MF]==0) {
   // do nothing
  } else
   amrex::Error("localMF_grow[spectralF_MF] invalid");

  delete_localMF(spectralF_MF,1);
  delete_localMF(spectralF_GP_MF,AMREX_SPACEDIM);

  if (localMF_grow[stableF_MF]==0) {
   // do nothing
  } else
   amrex::Error("localMF_grow[stableF_MF] invalid");

  delete_localMF(stableF_MF,1);
  delete_localMF(stableF_GP_MF,AMREX_SPACEDIM);

  if (localMF_grow[delta_MF]==0) {
   // do nothing
  } else
   amrex::Error("localMF_grow[delta_MF] invalid");

  delete_localMF(delta_MF,1);
  delete_localMF(delta_GP_MF,AMREX_SPACEDIM);

 } else if ((ns_time_order==1)&&(enable_spectral==0)) {
  // do nothing
 } else
  amrex::Error("ns_time_order or enable_spectral invalid");

}  // subroutine deallocate_SDC


// called before veldiffuseALL() from NavierStokes::do_the_advance
// Second half of D^{upside down triangle}/Dt
void NavierStokes::tensor_advection_updateALL() {

 int finest_level=parent->finestLevel();

 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=num_materials)) {

  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
    //ngrow=1
   ns_level.getState_localMF(HOLD_VELOCITY_DATA_MF,
     1,STATECOMP_VEL,STATE_NCOMP_VEL,cur_time_slab);
  }

   // init_gradu_tensorALL fills CELLTENSOR_MF using these steps:
   // 1. find all velocity derivatives at faces.
   // 2. interpolate derivatives from faces to cells using 1-sided
   //    interpolation in the case that e.g. lsleft(im_primary)>=0
   //    but lsright(im_primary)<0.
   //    (im_primary is the main material in the cell, lspoint(im_primary)>=0)
   // init_gradu_tensorALL declared in NavierStokes2.cpp
  int do_alloc=0;
  int simple_AMR_BC_flag_viscosity=1;
  init_gradu_tensorALL(
    HOLD_VELOCITY_DATA_MF,
    do_alloc,
    CELLTENSOR_MF,
    FACETENSOR_MF,
    simple_AMR_BC_flag_viscosity);

  allocate_array(0,DERIVE_TENSOR_NCOMP,-1,HOLD_GETSHEAR_DATA_MF);

  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   // get std::sqrt(2 D:D),D,grad U 
   // DERIVE_TENSOR_MAG: std::sqrt(2 * D : D)
   // DERIVE_TENSOR_RATE_DEFORM: D11,D12,D13,D21,D22,D23,D31,D32,D33
   // DERIVE_TENSOR_GRAD_VEL: ux,uy,uz,vx,vy,vz,wx,wy,wz
   int only_scalar=0; 
   int destcomp=0;
   int ngrow_zero=0;
   ns_level.level_getshear(
       ns_level.localMF[HOLD_GETSHEAR_DATA_MF],
       ns_level.localMF[HOLD_VELOCITY_DATA_MF],
       only_scalar,destcomp,ngrow_zero);
  }

   // tensor_advection_update is declared in: NavierStokes.cpp
  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   ns_level.tensor_advection_update();
  }
  avgDownALL_TENSOR();

  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   ns_level.tensor_extrapolation();
  }
  avgDownALL_TENSOR();

  delete_array(HOLD_GETSHEAR_DATA_MF);
  delete_array(HOLD_VELOCITY_DATA_MF);
  delete_array(CELLTENSOR_MF);
  delete_array(FACETENSOR_MF);

 } else if (num_materials_viscoelastic==0) {
  // do nothing
 } else
  amrex::Error("num_materials_viscoelastic bad tensor_advection_updateALL");

} // subroutine tensor_advection_updateALL

void NavierStokes::resize_maskfiner(int ngrow,int mask_id) {

 if (localMF_grow[mask_id]>=0) {
  // do nothing
 } else
  amrex::Error("localMF_grow[mask_id]<0");

 if (localMF[mask_id]->nGrow()!=ngrow) {
  Real tag=1.0;
  int clearbdry=0;
  maskfiner_localMF(mask_id,ngrow,tag,clearbdry);
 }

} // subroutine resize_maskfiner


void NavierStokes::resize_mask_nbr(int ngrow) {

 if (localMF_grow[MASK_NBR_MF]>=0) {
  // do nothing
 } else {
  amrex::Error("localMF_grow[MASK_NBR_MF]<0");
 }

 if (localMF[MASK_NBR_MF]->nGrow()!=ngrow) {
  prepare_mask_nbr(ngrow);
 }

} // subroutine resize_mask_nbr


void NavierStokes::resize_metrics(int ngrow) {

 if (localMF_grow[VOLUME_MF]>=0) {
  // do nothing
 } else
  amrex::Error("localMF_grow[VOLUME_MF]<0");

 if (localMF[VOLUME_MF]->nGrow()!=ngrow) {
   // NavierStokes::metrics_data declared in NavierStokes2.cpp
  metrics_data(ngrow);
 }

} // subroutine resize_metrics

void NavierStokes::init_delta_SDC() {

 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid init_delta_SDC");

 if (ns_time_order>=1) {
  // do nothing
 } else
  amrex::Error("ns_time_order must be >=1");

 if ((ns_time_order==1)&&(enable_spectral!=0))
  amrex::Error("(ns_time_order==1)&&(enable_spectral!=0)");
 if ((ns_time_order>=2)&&(enable_spectral!=1))
  amrex::Error("(ns_time_order>=2)&&(enable_spectral!=1)");

 if ((ns_time_order>=2)&&(enable_spectral==1)) {

  if (localMF[delta_MF]->nComp()==NSTATE_SDC*ns_time_order) {
   // do nothing
  } else
   amrex::Error("localMF[delta_MF]->nComp() invalid");

  if (localMF[stableF_MF]->nComp()==NSTATE_SDC*ns_time_order) {
   // do nothing
  } else
   amrex::Error("localMF[stableF_MF]->nComp() invalid");

  if (localMF[spectralF_MF]->nComp()==NSTATE_SDC*(ns_time_order+1)) {
   // do nothing
  } else
   amrex::Error("localMF[spectralF_MF]->nComp() invalid");

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

   if (localMF[delta_GP_MF+dir]->nComp()==ns_time_order) {
    // do nothing
   } else
    amrex::Error("localMF[delta_GP_MF]->nComp() invalid");

   if (localMF[stableF_GP_MF+dir]->nComp()==ns_time_order) {
    // do nothing
   } else
    amrex::Error("localMF[stableF_GP_MF]->nComp() invalid");

   if (localMF[spectralF_GP_MF+dir]->nComp()==ns_time_order+1) {
    // do nothing
   } else
    amrex::Error("localMF[spectralF_GP_MF]->nComp() invalid");

  } //dir=0..sdim-1

   // NSTATE_SDC is a macro declared in EXTRAP_COMP.H
  if (SDC_outer_sweeps==0) {
   setVal_localMF(delta_MF,0.0,0,NSTATE_SDC*ns_time_order,0);
   for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
    setVal_localMF(delta_GP_MF+dir,0.0,0,ns_time_order,0);
  } else if ((SDC_outer_sweeps>0)&&
             (SDC_outer_sweeps<ns_time_order)) {
   // do nothing
  } else
   amrex::Error("SDC_outer_sweeps invalid");

  setVal_localMF(stableF_MF,0.0,0,NSTATE_SDC*ns_time_order,0);
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
   setVal_localMF(stableF_GP_MF+dir,0.0,0,ns_time_order,0);

  int scomp=0;
  int ncomp=NSTATE_SDC*(ns_time_order+1);
  int scompGP=0;
  int ncompGP=ns_time_order+1;

  if ((SDC_outer_sweeps>0)&&
      (SDC_outer_sweeps<ns_time_order)) {
   // do not zap F(t^n) data.
   scomp=NSTATE_SDC;
   ncomp-=scomp;
   scompGP=1;
   ncompGP-=scompGP;
  } else if (SDC_outer_sweeps==0) {
   // ok to zap F(t^n) data.
  } else
   amrex::Error("SDC_outer_sweeps invalid");

  setVal_localMF(spectralF_MF,0.0,scomp,ncomp,0);
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
   setVal_localMF(spectralF_GP_MF+dir,0.0,scompGP,ncompGP,0);

 } else if ((ns_time_order==1)&&(enable_spectral==0)) {
  // do nothing
 } else
  amrex::Error("ns_time_order or enable_spectral invalid");

}  // subroutine init_delta_SDC

Real NavierStokes::advance(Real time,Real dt) {

 std::string local_caller_string="advance";

 if (ParallelDescriptor::IOProcessor()) 
  std::cout << "advance time= " << time << " dt= " << dt << '\n';

 int finest_level = parent->finestLevel();
 const int max_level = parent->maxLevel();

 if (finest_level<=max_level) {
  // do nothing
 } else
  amrex::Error("it is required that finest_level<=max_level");

 if (finest_level<=max_level_for_use) {
  // do nothing
 } else
  amrex::Error("it is required that finest_level<=max_level_for_use");

 Real dt_new=dt;
 int advance_status=1;

 if (level==0) {

  int nsteps=parent->levelSteps(0); 
  ULong cpu_seed=(ULong) nsteps;
  cpu_seed+=ParallelDescriptor::MyProc()+1; 
  amrex::InitRandom(cpu_seed,ParallelDescriptor::NProcs());

   //do {...} while (advance_status==0);
  do {

   interface_touch_flag=1; //advance

   SDC_outer_sweeps=0;
   slab_step=ns_time_order-1;
   SDC_setup_step(); 

   if ((time>=0.0)&&(time<=1.0)) {
    if (std::abs(upper_slab_time-time)<=CPP_EPS_12_5) {
     // do nothing
    } else {
     amrex::Error("upper_slab_time-time>tol (a)");
    }
   } else if (time>1.0) {
    if (std::abs(upper_slab_time-time)<=CPP_EPS_12_5*time) {
     // do nothing
    } else {
     amrex::Error("upper_slab_time-time>tol(time) (b)");
    }
   } else {
    amrex::Error("time invalid");
   }

   if (verbose>0) {
    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "calling metrics \n";
    }
   }

   if (1==0) {
    std::fflush(NULL);
    int proc=ParallelDescriptor::MyProc();
    std::cout << "prior to metrics_dataALL on processor=" << proc << '\n';
    std::fflush(NULL);
    ParallelDescriptor::Barrier();
   }

    //declared in: NavierStokes2.cpp
   metrics_dataALL(1);  

   if (verbose>0) {
    std::fflush(NULL);
    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "initializing masks \n";
    }
    std::fflush(NULL);
    ParallelDescriptor::Barrier();
   }

   for (int ilev=level;ilev<=finest_level;ilev++) {
    NavierStokes& ns_level=getLevel(ilev);
    //mask=tag if not covered by level+1 or outside the domain.
    Real tag=1.0;
    int clearbdry=0; 
    ns_level.maskfiner_localMF(MASKCOEF_MF,1,tag,clearbdry);
    ns_level.prepare_mask_nbr(1);
   }

   build_masksemALL();

   if (verbose>0) {
    std::fflush(NULL);
    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "prescribe solid geometry (after regridding) \n";
    }
    std::fflush(NULL);
    ParallelDescriptor::Barrier();
   }

   if (perturbation_on_restart==1) {

    perturbation_on_restart=0;

    for (int ilev=finest_level;ilev>=level;ilev--) {
     NavierStokes& ns_level=getLevel(ilev);
     ns_level.add_perturbation();
    } 
   
   } else if (perturbation_on_restart==0) {
    // do nothing
   } else
    amrex::Error("perturbation_on_restart invalid");
  

   // take care of AMR grid change.

   make_MAC_velocity_consistentALL();
  
    // velocity and pressure
   avgDownALL(State_Type,STATECOMP_VEL,STATE_NCOMP_VEL+STATE_NCOMP_PRES,1);
    // "state" (all materials)
   avgDownALL(State_Type,STATECOMP_STATES,num_state_material*num_materials,1);
    // expected "DIV" 
   avgDownALL(DIV_Type,0,1,1);

    // in: advance
    // calls MOFavgDown, LS_Type avgDown
    // projects volume fractions so that sum F_m_fluid=1.
   int renormalize_only=0;
   int local_truncate=0;
   int update_particles=0;
   prescribe_solid_geometryALL(upper_slab_time,renormalize_only,
      local_truncate,local_caller_string,update_particles);

   if (verbose>0) {
    std::fflush(NULL);
    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "copy new to old ... \n";
    }
    std::fflush(NULL);
    ParallelDescriptor::Barrier();
   }

#ifdef AMREX_PARTICLES

   int lev_min=0;
   int lev_max=-1;
   int nGrow_Redistribute=0;
   int local_Redistribute=0;
   bool remove_negative=true;

   NavierStokes& ns_level0=getLevel(0);
   My_ParticleContainer& old_PC=ns_level0.newDataPC(ns_time_order);
   old_PC.Redistribute(lev_min,lev_max,nGrow_Redistribute, 
     local_Redistribute,remove_negative);

#endif

   //copy bfact_time_order component to the
   //components: 0..bfact_time_order-1
   CopyNewToOldALL();

   interface_touch_flag=1; //advance

   // new_time=time+dt_new  old_time=time
   for (int ilev=level;ilev<=finest_level;ilev++) {
    NavierStokes& ns_level=getLevel(ilev);
    ns_level.setTimeLevel(time+dt_new,dt_new);
   }

   SDC_setup_step();

   Real time_scale=1.0;
   if (upper_slab_time>time_scale)
    time_scale=upper_slab_time;
   time_scale*=CPP_EPS_10_5;

   if (std::abs(upper_slab_time-lower_slab_time-dt_new)<=time_scale) {
    // do nothing
   } else
    amrex::Error("SDC_setup_step failed");

   if (std::abs(lower_slab_time-time)<=time_scale) {
    // do nothing
   } else
    amrex::Error("lower_slab_time set improperly");

   do_the_advance(lower_slab_time,dt_new,advance_status);

   if (advance_status==1) {
    // do nothing (success)
   } else if (advance_status==0) { // failure
    dt_new=0.5*dt_new;
     //copy component "0" to components 1..bfact_time_order.
    CopyOldToNewALL();

    interface_touch_flag=1; //advance

    for (int ilev=level;ilev<=finest_level;ilev++) {
     NavierStokes& ns_level=getLevel(ilev);
      // new_time=lower_slab_time  old_time=lower_slab_time-dt_new
     ns_level.setTimeLevel(lower_slab_time,dt_new);
    }
    parent->setDt(dt_new);
   } else
    amrex::Error("advance_status invalid");

   delete_array(MASKCOEF_MF);
   delete_array(MASK_NBR_MF);

  } while (advance_status==0);

 } else
  amrex::Error("advance should only be called at level=0");

 return dt_new;

} // end subroutine advance


// called from: NavierStokes::do_the_advance
// delta=integral_tn^tnp1  f^spectral dt - deltatn F^stable
void NavierStokes::init_splitting_force_SDC() {

 std::string local_caller_string="init_splitting_force_SDC";
 
 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();

 if (divu_outer_sweeps==0) {
  // do nothing
 } else
  amrex::Error("divu_outer_sweeps invalid");

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid");

 if (ns_time_order>=2) {

  if (slab_step!=ns_time_order)
   amrex::Error("slab_step invalid");

  if (num_state_base!=2)
   amrex::Error("num_state_base invalid");

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   debug_ngrow(spectralF_GP_MF+dir,0,local_caller_string);
   debug_ngrow(stableF_GP_MF+dir,0,local_caller_string);
   debug_ngrow(delta_GP_MF+dir,0,local_caller_string);
  }
  debug_ngrow(spectralF_MF,0,local_caller_string);
  debug_ngrow(stableF_MF,0,local_caller_string);
  debug_ngrow(delta_MF,0,local_caller_string);

  MultiFab& S_new=get_new_data(State_Type,slab_step+1);

  int nstate=STATE_NCOMP;
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

   FArrayBox& HOfab=(*localMF[spectralF_MF])[mfi];
   FArrayBox& LOfab=(*localMF[stableF_MF])[mfi];
   FArrayBox& deltafab=(*localMF[delta_MF])[mfi];
   FArrayBox& masksem=(*localMF[MASKSEM_MF])[mfi];

   int HOncomp=HOfab.nComp();
   int LOncomp=LOfab.nComp();
   int delta_ncomp=deltafab.nComp();

   if ((HOncomp!=(ns_time_order+1)*NSTATE_SDC)||
       (LOncomp!=ns_time_order*NSTATE_SDC)||
       (delta_ncomp!=ns_time_order*NSTATE_SDC))
    amrex::Error("HOncomp, LOncomp, or delta_ncomp invalid");

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    //declared in: GODUNOV_3D.F90
   fort_sdc_time_quad(
    &HOncomp,
    &LOncomp,
    &delta_ncomp,
    &nstate,
    xlo,dx,
    deltafab.dataPtr(),
    ARLIM(deltafab.loVect()),ARLIM(deltafab.hiVect()),
    HOfab.dataPtr(),
    ARLIM(HOfab.loVect()),ARLIM(HOfab.hiVect()),
    LOfab.dataPtr(),
    ARLIM(LOfab.loVect()),ARLIM(LOfab.hiVect()),
    masksem.dataPtr(),
    ARLIM(masksem.loVect()),ARLIM(masksem.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    &level,
    &finest_level,
    &delta_slab_time);
  }  // mfi  
} // omp
  ns_reconcile_d_num(LOOP_SDC_TIME_QUAD,"init_splitting_force_SDC");

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

    FArrayBox& HOfab=(*localMF[spectralF_GP_MF+dir])[mfi];
    FArrayBox& LOfab=(*localMF[stableF_GP_MF+dir])[mfi];
    FArrayBox& deltafab=(*localMF[delta_GP_MF+dir])[mfi];
    FArrayBox& masksem=(*localMF[MASKSEM_MF])[mfi];

    int HOncomp=HOfab.nComp();
    int LOncomp=LOfab.nComp();
    int delta_ncomp=deltafab.nComp();

    if ((HOncomp!=(ns_time_order+1))||
        (LOncomp!=ns_time_order)||
        (delta_ncomp!=ns_time_order))
     amrex::Error("HOncomp, LOncomp, or delta_ncomp invalid");

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     //declared in: GODUNOV_3D.F90
    fort_sdc_time_quad_face(
     &dir,
     &HOncomp,
     &LOncomp,
     &delta_ncomp,
     &nstate,
     xlo,dx,
     deltafab.dataPtr(),
     ARLIM(deltafab.loVect()),ARLIM(deltafab.hiVect()),
     HOfab.dataPtr(),
     ARLIM(HOfab.loVect()),ARLIM(HOfab.hiVect()),
     LOfab.dataPtr(),
     ARLIM(LOfab.loVect()),ARLIM(LOfab.hiVect()),
     masksem.dataPtr(),
     ARLIM(masksem.loVect()),ARLIM(masksem.hiVect()),
     tilelo,tilehi,
     fablo,fabhi,&bfact,
     &level,
     &finest_level,
     &delta_slab_time);
   }  // mfi  
} // omp
   ns_reconcile_d_num(LOOP_SDC_TIME_QUAD_FACE,"init_splitting_force_SDC");

  } // dir=0..sdim-1

 } else {
  amrex::Error("ns_time_order invalid init_splitting_force_SDC");
 }

} // subroutine init_splitting_force_SDC

//source_term==SUB_OP_SDC_LOW_TIME:
//  slab_step (aka "k") = 0,1,2,...,ns_time_order
//source_term==SUB_OP_SDC_ISCHEME:
//  slab_step (aka "k") = 0,1,2,...,ns_time_order-1
//SDC_outer_sweeps=0...ns_time_order-1
//
//source_term==1=SUB_OP_SDC_LOW_TIME => compute F(t^{n+k/order})
//source_term==0=SUB_OP_SDC_ISCHEME  => compute 
//  (a) F(t^{n+k/order,(0)})  (SUB_OP_ISCHEME_PREDICT)
//  (b) F(t^{n+k/order,(1)})  (SUB_OP_ISCHEME_CORRECT)
void NavierStokes::SEM_advectALL(int source_term) {

 std::string local_caller_string="SEM_advectALL";

 if ((SDC_outer_sweeps>=0)&&(SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid");

 if (stokes_flow==0) {

  int finest_level=parent->finestLevel();
  if ((ns_time_order==1)&&(enable_spectral!=0))
   amrex::Error("(ns_time_order==1)&&(enable_spectral!=0)");
  if ((ns_time_order>=2)&&(enable_spectral!=1))
   amrex::Error("(ns_time_order>=2)&&(enable_spectral!=1)");

  if ((ns_time_order>=2)&&(enable_spectral==1)) {

   int operation_flag=OP_ISCHEME_MAC;

   int SEM_end_spectral_loop=2;

   prescribed_vel_time_slab=prev_time_slab;
   vel_time_slab=prev_time_slab;

   if (source_term==SUB_OP_SDC_LOW_TIME) {

    if ((slab_step>=0)&&(slab_step<=ns_time_order)) {
     // do nothing
    } else
     amrex::Error("slab_step invalid");

    vel_time_slab=prev_time_slab;
   } else if (source_term==SUB_OP_SDC_ISCHEME) {

    if ((slab_step>=0)&&(slab_step<ns_time_order)) {
     // do nothing
    } else
     amrex::Error("slab_step invalid");

    if (divu_outer_sweeps==0) 
     vel_time_slab=prev_time_slab;
    else if (divu_outer_sweeps>0)
     vel_time_slab=cur_time_slab;
    else
     amrex::Error("divu_outer_sweeps invalid SEM_advectALL");
   } else
    amrex::Error("source_term invalid");

   for (int ilev=finest_level;ilev>=level;ilev--) {
    NavierStokes& ns_level=getLevel(ilev);
    for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
      //ngrow,dir,time
      //Umac_Type
     ns_level.getStateMAC_localMF(UMAC_MF+dir,0,dir,vel_time_slab);
    }  //dir=0,...,sdim-1
   } //ilev=finest_level ... level

   int advect_iter_max=2;
   if (source_term==SUB_OP_SDC_LOW_TIME) {
    advect_iter_max=1;
   } else if (source_term==SUB_OP_SDC_ISCHEME) {
    advect_iter_max=2;
   } else
    amrex::Error("advect_iter_max invalid");
  
   for (advect_iter=0;advect_iter<advect_iter_max;advect_iter++) {

    if (source_term==SUB_OP_SDC_LOW_TIME) { 
     advect_time_slab=prev_time_slab;
    } else if (source_term==SUB_OP_SDC_ISCHEME) {
     if (advect_iter==SUB_OP_ISCHEME_PREDICT) {
      advect_time_slab=prev_time_slab;
     } else if (advect_iter==SUB_OP_ISCHEME_CORRECT) {
      advect_time_slab=cur_time_slab;
     } else
      amrex::Error("advect_iter invalid");
    } else
     amrex::Error("source_term invalid");

    for (int ilev=finest_level;ilev>=level;ilev--) {
     NavierStokes& ns_level=getLevel(ilev);

     ns_level.getStateDen_localMF(DEN_RECON_MF,1,advect_time_slab);
     ns_level.getState_localMF(VELADVECT_MF,1,0,
      AMREX_SPACEDIM,advect_time_slab); 
     for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
      ns_level.new_localMF(AMRSYNC_PRES_MF+dir,NFLUXSEM,0,dir);
      ns_level.setVal_localMF(AMRSYNC_PRES_MF+dir,1.0e+30,0,NFLUXSEM,0);
      ns_level.new_localMF(CONSERVE_FLUXES_MF+dir,NFLUXSEM,0,dir);
      ns_level.setVal_localMF(CONSERVE_FLUXES_MF+dir,1.0e+30,0,NFLUXSEM,0);
      ns_level.new_localMF(COARSE_FINE_FLUX_MF+dir,NFLUXSEM,0,dir);
      ns_level.setVal_localMF(COARSE_FINE_FLUX_MF+dir,1.0e+30,0,NFLUXSEM,0);
     } // dir=0..sdim-1
     ns_level.resize_levelset(2,LEVELPC_MF);
     ns_level.VOF_Recon_resize(1); //output:SLOPE_RECON_MF
     ns_level.resize_maskfiner(1,MASKCOEF_MF);
     ns_level.resize_mask_nbr(1);
     ns_level.resize_metrics(1);
     ns_level.allocate_flux_register(operation_flag);
    } // ilev=finest_level ... level

     // spectral_loop==0 => fine data transferred to coarse in a special way
     // spectral_loop==1 => coarse fluxes interpolated to fine level.
    int init_fluxes=1;
    for (int spectral_loop=0;spectral_loop<SEM_end_spectral_loop;
         spectral_loop++) {
      for (int tileloop=0;tileloop<=1;tileloop++) {
       for (int ilev=finest_level;ilev>=level;ilev--) {
        NavierStokes& ns_level=getLevel(ilev);
        ns_level.SEM_scalar_advection(init_fluxes,source_term, 
	spectral_loop,tileloop);
       } // ilev=finest_level ... level
      } // tileloop=0..1

      for (int ilev=finest_level;ilev>=level;ilev--) {
       NavierStokes& ns_level=getLevel(ilev);
       if (spectral_loop==SEM_end_spectral_loop-1) {
        ns_level.delete_localMF(SEM_FLUXREG_MF,1);
       } else if (spectral_loop==0) {
        ns_level.localMF[SEM_FLUXREG_MF]->FillBoundary(geom.periodicity());
       } else
        amrex::Error("spectral_loop invalid");
      } // ilev=finest_level ... level

    } // spectral_loop=0...SEM_end_spectral_loop-1

    // spectral_override==0 => always do low order average down.
    for (int ilev=finest_level-1;ilev>=level;ilev--) {
     NavierStokes& ns_level=getLevel(ilev);
     int ncomp_edge=-1;
     int scomp=0;
     int spectral_override=SPECTRAL_ORDER_AVGDOWN;
     ns_level.avgDownEdge_localMF(CONSERVE_FLUXES_MF,scomp,ncomp_edge,
       0,AMREX_SPACEDIM,spectral_override,local_caller_string);
    } // ilev=finest_level-1 ... level

    init_fluxes=0;
    int spectral_loop=0;
    int tileloop=0;
    for (int ilev=finest_level;ilev>=level;ilev--) {
     NavierStokes& ns_level=getLevel(ilev);
     ns_level.SEM_scalar_advection(init_fluxes,source_term,
        spectral_loop,tileloop);
    } // ilev=finest_level ... level

    avgDownALL(State_Type,STATECOMP_VEL,STATE_NCOMP_VEL+STATE_NCOMP_PRES,1);
    avgDownALL(State_Type,STATECOMP_STATES,num_state_material*num_materials,1);

    for (int ilev=finest_level;ilev>=level;ilev--) {
     NavierStokes& ns_level=getLevel(ilev);
     ns_level.delete_localMF(AMRSYNC_PRES_MF,AMREX_SPACEDIM);
     ns_level.delete_localMF(COARSE_FINE_FLUX_MF,AMREX_SPACEDIM);
     ns_level.delete_localMF(CONSERVE_FLUXES_MF,AMREX_SPACEDIM);
     ns_level.delete_localMF(DEN_RECON_MF,1);
     ns_level.delete_localMF(VELADVECT_MF,1);
    }
  
   } // advect_iter=0 ... advect_iter_max-1

   for (int ilev=finest_level;ilev>=level;ilev--) {
    NavierStokes& ns_level=getLevel(ilev);
    ns_level.delete_localMF(UMAC_MF,AMREX_SPACEDIM);
   }

  } else
   amrex::Error("enable_spectral invalid");

 } else if (stokes_flow==1) {
  // do nothing
 } else
  amrex::Error("stokes_flow invalid");

} // subroutine SEM_advectALL

// called from NavierStokes::do_the_advance first thing in the loop
// divu_outer_sweeps=0..local_num_divu_outer_sweeps-1
void NavierStokes::prelim_alloc() {

 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid prelim_alloc");


 int nsolve=1;

 allocate_physics_vars();  

 if (localMF[MDOT_MF]->nComp()!=nsolve)
  amrex::Error("localMF[MDOT_MF]->nComp() invalid");

   //val,scomp,ncomp,ngrow
 setVal_localMF(MDOT_MF,0.0,0,nsolve,0); 

} // end subroutine prelim_alloc

// called from: NavierStokes::do_the_advance
void NavierStokes::advance_MAC_velocity(int project_option) {

 if ((project_option==SOLVETYPE_PRES)||  // is_zalesak()==FALSE
     (project_option==SOLVETYPE_INITPROJ)) { //is_zalesak()==TRUE
  // do nothing
 } else
  amrex::Error("project_option invalid advance_MAC_velocity");

 int idx_velcell=-1;
 Real beta=0.0;

 // unew^{f}=
 // (i) unew^{f} in non-solid regions
 // (ii) u^{f,save} + (unew^{c}-u^{c,save})^{c->f} in spectral 
 //      regions 
 //      (u^{c,save} = *localMF[ADVECT_REGISTER_MF])
 //      (u^{f,save} = *localMF[ADVECT_REGISTER_FACE_MF+dir])
 // (iii) usolid in solid regions
 int operation_flag=OP_U_SEM_CELL_MAC_TO_MAC;  

 Vector<blobclass> blobdata;

 increment_face_velocityALL(
   operation_flag,project_option,
   idx_velcell,beta,blobdata);

} // end subroutine advance_MAC_velocity()

void NavierStokes::pressure_gradient_code_segment(
  const std::string& caller_string) {

 std::string local_caller_string="pressure_gradient_code_segment";
 local_caller_string=caller_string+local_caller_string;

#if (NS_profile_solver==1)
 BLProfiler bprof(local_caller_string);
#endif

  // MDOT term included
 multiphase_project(SOLVETYPE_PRES);

 int singular_parts_exist=0;
 for (int im=0;im<num_materials;im++) {
  if (is_singular_coeff(im)==0) {
   // do nothing
  } else if (is_singular_coeff(im)==1) {
   singular_parts_exist=1;
  } else
   amrex::Error("is_singular_coeff invalid");
 } // im=0..num_materials-1

 if (singular_parts_exist==1) {

  if (extend_pressure_into_solid==1) {
   if (FSI_outer_sweeps==num_FSI_outer_sweeps-1) {
    multiphase_project(SOLVETYPE_PRESEXTRAP);
   } else if ((FSI_outer_sweeps>=0)&&
              (FSI_outer_sweeps<num_FSI_outer_sweeps-1)) {
    //do nothing
   } else
    amrex::Error("FSI_outer_sweeps invalid");

  } else if (extend_pressure_into_solid==0) {
   // do nothing
  } else
   amrex::Error("extend_pressure_into_solid invalid");

 } else if (singular_parts_exist==0) {
  // do nothing
 } else
  amrex::Error("singular_parts_exist invalid");

#if (NS_profile_solver==1)
 bprof.stop();
#endif

} // end subroutine pressure_gradient_code_segment

void NavierStokes::phase_change_code_segment(
  const std::string& caller_string,
  int& color_count,
  Vector<blobclass>& blobdata) {

 if ((slab_step>=0)&&(slab_step<ns_time_order)) {
  //do nothing
 } else 
  amrex::Error("slab_step invalid");

 int finest_level=parent->finestLevel();
 if (level==0) {
  //do nothing
 } else
  amrex::Error("level invalid phase_change_code_segment");

 std::string local_caller_string="phase_change_code_segment";
 local_caller_string=caller_string+local_caller_string;

#if (NS_profile_solver==1)
 BLProfiler bprof(local_caller_string);
#endif

 if (ngrow_make_distance!=3)
  amrex::Error("ngrow_make_distance!=3");
 if (ngrow_distance!=4)
  amrex::Error("ngrow_distance!=4");

  // first num_interfaces components correspond to the status.
 int nburning=EXTRAP_NCOMP_BURNING;
 int ntsat=EXTRAP_NCOMP_TSAT;

  // SATURATION_TEMP_MF is passed to the following fortran
  // routines:
  //  RATEMASSCHANGE,
  //  AVGDOWN_BURNING, 
  //  EXT_BURNVEL_INTERP,
  //  EXTEND_BURNING_VEL,
  //  CONVERTMATERIAL
  //  STEFANSOLVER
  // BURNING_VELOCITY_MF is passed to the following fortran
  // routines:
  //  RATEMASSCHANGE,
  //  AVGDOWN_BURNING, 
  //  EXT_BURNVEL_INTERP,
  //  EXTEND_BURNING_VEL,
  //  NODEDISPLACE,
  //  CONVERTMATERIAL
 for (int ilev=level;ilev<=finest_level;ilev++) {

  NavierStokes& ns_level=getLevel(ilev);

  ns_level.new_localMF(BURNING_VELOCITY_MF,nburning,
    ngrow_distance,-1);
  ns_level.setVal_localMF(BURNING_VELOCITY_MF,0.0,0,
    nburning,ngrow_distance);

  int n_normal=(num_materials+num_interfaces)*(AMREX_SPACEDIM+1);

  ns_level.new_localMF(FD_NRM_ND_MF,n_normal,
    ngrow_make_distance+1,-1);
  ns_level.setVal_localMF(FD_NRM_ND_MF,0.0,0,
    n_normal,ngrow_make_distance+1);

   // first num_materials+num_interfaces components are curvature
   // second num_materials+num_interfaces components are 
   // status (0=bad 1=good)
  ns_level.new_localMF(FD_CURV_CELL_MF,2*(num_materials+num_interfaces),
    ngrow_make_distance,-1);
  ns_level.setVal_localMF(FD_CURV_CELL_MF,0.0,0,
    2*(num_materials+num_interfaces),ngrow_make_distance);

  ns_level.new_localMF(SATURATION_TEMP_MF,ntsat,
    ngrow_distance,-1);
  ns_level.setVal_localMF(SATURATION_TEMP_MF,0.0,0,
    ntsat,ngrow_distance);

  ns_level.new_localMF(JUMP_STRENGTH_MF,2*num_interfaces,
  		     ngrow_distance,-1); 
  ns_level.setVal_localMF(JUMP_STRENGTH_MF,0.0,0,
  		        2*num_interfaces,ngrow_distance);

 } // ilev=level ... finest_level

 debug_ngrow(JUMP_STRENGTH_MF,ngrow_distance,local_caller_string);
 debug_ngrow(SWEPT_CROSSING_MF,0,local_caller_string);
 debug_ngrow(BURNING_VELOCITY_MF,ngrow_distance,local_caller_string);
 debug_ixType(BURNING_VELOCITY_MF,-1,local_caller_string);
 debug_ngrow(SATURATION_TEMP_MF,ngrow_distance,local_caller_string);
 debug_ngrow(FD_NRM_ND_MF,ngrow_make_distance+1,local_caller_string);
 debug_ixType(FD_NRM_ND_MF,-1,local_caller_string);
 debug_ngrow(FD_CURV_CELL_MF,ngrow_make_distance,local_caller_string);
 debug_ixType(FD_CURV_CELL_MF,-1,local_caller_string);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  debug_ngrow(AREA_MF+dir,1,local_caller_string);
  debug_ngrow(FACE_VAR_MF+dir,0,local_caller_string);
 }
 debug_ngrow(MDOT_MF,0,local_caller_string);

 zeroALL(ngrow_distance,nburning,BURNING_VELOCITY_MF);
 zeroALL(ngrow_distance,2*num_interfaces,JUMP_STRENGTH_MF);
  //ngrow,scomp,ncomp,val,dest_mf
 setVal_array(0,0,num_materials,1.0,SWEPT_CROSSING_MF);
  // piecewise constant interpolation at coarse/fine borders.
  // fluid LS can be positive in the solid regions.
  // HOLD_LS_DATA_MF is deleted in phase_change_redistributeALL()
 allocate_levelset_ALL(ngrow_distance,HOLD_LS_DATA_MF);
 if (localMF[HOLD_LS_DATA_MF]->nComp()!=num_materials*(AMREX_SPACEDIM+1))
  amrex::Error("hold_LS_DATA_MF (nComp())!=num_materials*(AMREX_SPACEDIM+1)");
 debug_ngrow(HOLD_LS_DATA_MF,ngrow_distance,local_caller_string);

  // BURNING_VELOCITY_MF flag==+ or - 1 if valid rate of phase change.
 for (int ilev=level;ilev<=finest_level;ilev++) {
  int nucleation_flag=0;
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.level_phase_change_rate(blobdata,color_count,
    nucleation_flag);
 }

 delete_array(TYPE_MF);
 delete_array(COLOR_MF);

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
    // in: NavierStokes2.cpp
  ns_level.avgDownBURNING_localMF(
  		BURNING_VELOCITY_MF,
  		SATURATION_TEMP_MF);
  ns_level.avgDown(LS_Type,0,num_materials,0);
 }

   // FIND RATE OF PHASE CHANGE V=[k grad T]/L for fully saturated
   // boiling for example.
   // EVAPORATION (partially saturated), or
   // EVAPORATION (fully saturated), or
   // BOILING (fully saturated), or
   // CONDENSATION (partially saturated), .....
   // traverse from coarsest to finest so that coarse/fine
   // BC are well defined.
   // sets the burning velocity flag from 0 to 2 if
   // foot of characteristic within range.
   // calls fort_extend_burning_vel
 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.level_phase_change_rate_extend();
 }

 if (visual_phase_change_plot_int>0) {
  if (very_last_sweep==1) {
   // nsteps==0 very first step.
   // in: phase_change_code_segment
   int nsteps=parent->levelSteps(0); 
   int ratio=(nsteps+1)/visual_phase_change_plot_int;
   ratio=ratio*visual_phase_change_plot_int;
   if (ratio==nsteps+1) {

     //writeSanityCheckData outputs raw data that exists on the
     //computational domain boundary or within.
     //TY_GAMMA<stuff>.plt (visit can open binary tecplot files)
    writeSanityCheckData(
     "TY_GAMMA",
     "SATURATION_TEMP_MF: flag12,flag13,flag23,T_GAMMA12,Y_GAMMA12, ...",
     local_caller_string,
     SATURATION_TEMP_MF, //tower_mf_id
     localMF[SATURATION_TEMP_MF]->nComp(), 
     SATURATION_TEMP_MF,
     -1,  // State_Type==-1 
     -1, // data_dir==-1 (cell centered)
     parent->levelSteps(0)); 

     //BURNVEL<stuff>.plt (visit can open binary tecplot files)
    writeSanityCheckData(
     "BURNVEL",
     "BURNING_VELOCITY_MF: flag12,flag13,flag23,[xyz]V12,[xyz]V13, ..",
     local_caller_string,
     BURNING_VELOCITY_MF, //tower_mf_id
     localMF[BURNING_VELOCITY_MF]->nComp(), 
     BURNING_VELOCITY_MF,
     -1,  // State_Type==-1 
     -1,  // data_dir==-1 (cell centered)
     parent->levelSteps(0)); 

     //BURNVEL<stuff>.plt (visit can open binary tecplot files)
    writeSanityCheckData(
     "CURV_CELL",
     "FD_CURV_CELL_MF:curv:1..num_materials+num_interfaces stat:num_materials+num_interfaces+1..2(num_materials+nsten)",
     local_caller_string,
     FD_CURV_CELL_MF, //tower_mf_id
     localMF[FD_CURV_CELL_MF]->nComp(), 
     FD_CURV_CELL_MF,
     -1,  // State_Type==-1 
     -1,  // data_dir==-1 (cell centered)
     parent->levelSteps(0)); 

   } // (ratio==nsteps+1) {
  } else if (very_last_sweep==0) {
   // do nothing
  } else
   amrex::Error("very_last_sweep invalid");

 } else if (visual_phase_change_plot_int==0) {
  // do nothing
 } else
  amrex::Error("visual_phase_change_plot_int invalid");

 allocate_array(1,2*num_interfaces*AMREX_SPACEDIM,-1,nodevel_MF);
  //ngrow,scomp,ncomp
 setVal_array(1,0,2*num_interfaces*AMREX_SPACEDIM,0.0,nodevel_MF);

 delta_mass.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  // source 1..num_materials  dest 1..num_materials
  delta_mass[tid].resize(2*num_materials); 
  for (int im=0;im<2*num_materials;im++)
   delta_mass[tid][im]=0.0;
 } // tid

  // in: phase_change_code_segment
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.getStateDen_localMF(DEN_RECON_MF,1,cur_time_slab);
 }

 // 1.initialize node velocity from BURNING_VELOCITY_MF
 // 2.unsplit advection of materials changing phase
 // 3.update volume fractions, jump strength, temperature
 level_phase_change_convertALL();

 delete_array(FD_NRM_ND_MF);
 delete_array(FD_CURV_CELL_MF);

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   for (int im=0;im<num_materials;im++) {
    std::cout << "convert statistics: im,source,dest " << im << ' ' <<
     delta_mass[0][im] << ' ' << delta_mass[0][im+num_materials] << '\n';
   }
  }
 }

  // in: phase_change_code_segment
 delete_array(DEN_RECON_MF);


#ifdef AMREX_PARTICLES

 if ((slab_step>=0)&&(slab_step<ns_time_order)) {

  My_ParticleContainer& localPC=newDataPC(slab_step+1);

  int lev_min=0;
  int lev_max=-1;
  int nGrow_Redistribute=0;
  int local_Redistribute=0; 
  bool remove_negative=true;

  localPC.Redistribute(lev_min,lev_max,
    nGrow_Redistribute,local_Redistribute,remove_negative);

  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   int splitting_dir=-1;
    //move_particles() declared in NavierStokes2.cpp
   ns_level.move_particles(splitting_dir,localPC,local_caller_string);
  }

  localPC.Redistribute(lev_min,lev_max,
    nGrow_Redistribute,local_Redistribute,remove_negative);

 } else
  amrex::Error("slab_step invalid");

#endif

 delete_array(BURNING_VELOCITY_MF);
 delete_array(nodevel_MF);

 interface_touch_flag=1; //phase_change_code_segment

 int init_vof_prev_time=0;
 //output:SLOPE_RECON_MF
 VOF_Recon_ALL(
    local_caller_string, //phase_change_code_segment
    cur_time_slab,
    RECON_UPDATE_STATE_ERR_AND_CENTROID,
    init_vof_prev_time);

  // in: phase_change_code_segment
  // 1. prescribe solid temperature, velocity, and geometry where
  //    appropriate.
  // 2. extend level set functions into the solid.
 int renormalize_only=0;
 int local_truncate=0;
 int update_particles=1;
 prescribe_solid_geometryALL(cur_time_slab,renormalize_only,
  local_truncate,local_caller_string,update_particles);

 makeStateDistALL(update_particles);

#if (NS_profile_solver==1)
 bprof.stop();
#endif

} //end subroutine phase_change_code_segment

void NavierStokes::no_mass_transfer_code_segment(
  const std::string& caller_string) {

 if (level==0) {
  //do nothing
 } else
  amrex::Error("level invalid no_mass_transfer_code_segment");

 std::string local_caller_string="no_mass_transfer_code_segment";
 local_caller_string=caller_string+local_caller_string;

#if (NS_profile_solver==1)
 BLProfiler bprof(local_caller_string);
#endif

 interface_touch_flag=1; //no_mass_transfer_code_segment

 int init_vof_prev_time=0;
 //output:SLOPE_RECON_MF
 VOF_Recon_ALL(
   local_caller_string, //no_mass_transfer_code_segment
   cur_time_slab,
   RECON_UPDATE_STATE_ERR_AND_CENTROID,
   init_vof_prev_time);

 int update_particles=1;
 makeStateDistALL(update_particles);

#if (NS_profile_solver==1)
 bprof.stop();
#endif

} // end subroutine no_mass_transfer_code_segment

void NavierStokes::nucleation_code_segment(
  const std::string& caller_string,
  int& color_count,
  Vector<blobclass>& blobdata,
  Vector< Vector<Real> >& mdot_data,
  Vector< Vector<Real> >& mdot_comp_data,
  Vector< Vector<Real> >& mdot_data_redistribute,
  Vector< Vector<Real> >& mdot_comp_data_redistribute,
  Vector<int>& type_flag) {

 int finest_level=parent->finestLevel();
 if (level==0) {
  //do nothing
 } else
  amrex::Error("level invalid nucleation_code_segment");

 std::string local_caller_string="nucleation_code_segment";
 local_caller_string=caller_string+local_caller_string;

#if (NS_profile_solver==1)
 BLProfiler bprof(local_caller_string);
#endif

 if (1==0) {
  int basestep_debug=nStep();
  parent->writeDEBUG_PlotFile(
    basestep_debug,
    SDC_outer_sweeps,
    slab_step,
    divu_outer_sweeps);
  std::cout << "press any number then enter: before nucleate_bubbles\n";
  int n_input;
  std::cin >> n_input;
 }

  // CREATE SEEDS, NUCLEATION.
 for (int ilev=level;ilev<=finest_level;ilev++) {
  int nucleation_flag=1;
  color_count=1; // filler
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.level_phase_change_rate(blobdata,color_count,
    nucleation_flag);
 }

 delta_mass.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  // source 1..num_materials  dest 1..num_materials
  delta_mass[tid].resize(2*num_materials); 
  for (int im=0;im<2*num_materials;im++)
   delta_mass[tid][im]=0.0;
 } // tid=0 ... nthreads-1

 ParallelDescriptor::Barrier();

 int tessellate=1;
 int idx_mdot=-1; //idx_mdot==-1 => do not collect auxiliary data.
 int operation_flag=OP_GATHER_MDOT;
 int coarsest_level=0;

  //calling from: NavierStokes::nucleation_code_segment()
 ColorSumALL( 
   operation_flag, //=OP_GATHER_MDOT
   tessellate, //=1
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

 ParallelDescriptor::Barrier();

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.avgDown(LS_Type,0,num_materials,0);
  ns_level.MOFavgDown();
  ns_level.avgDown(State_Type,STATECOMP_STATES,
     num_state_material*num_materials,1);
 }  // ilev=finest_level ... level  

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   for (int im=0;im<num_materials;im++) {
    std::cout << "Nucleation stats: im,source,dest " << im << ' ' <<
     delta_mass[0][im] << ' ' << delta_mass[0][im+num_materials] << '\n';
   }
  }
 } 

 interface_touch_flag=1; //nucleation_code_segment

 // generates SLOPE_RECON_MF
 int init_vof_prev_time=0;
  // Fluids tessellate; solids overlay.
  // output:SLOPE_RECON_MF
 VOF_Recon_ALL(
   local_caller_string,  //nucleation_code_segment
   cur_time_slab,
   RECON_UPDATE_STATE_CENTROID,init_vof_prev_time);

 int update_particles=1;
 makeStateDistALL(update_particles);

 make_physics_varsALL(SOLVETYPE_PRES,local_caller_string); 
 delete_array(CELLTENSOR_MF);
 delete_array(FACETENSOR_MF);

 if (1==0) {
  int basestep_debug=nStep();
  parent->writeDEBUG_PlotFile(
    basestep_debug,
    SDC_outer_sweeps,
    slab_step,
    divu_outer_sweeps);
  std::cout << "press any number then enter: after nucleate_bubbles\n";
  int n_input;
  std::cin >> n_input;
 }
#if (NS_profile_solver==1)
 bprof.stop();
#endif

} //end subroutine nucleation_code_segment()

// called from: NavierStokes::advance
void NavierStokes::do_the_advance(Real timeSEM,Real dtSEM,
  int& advance_status) {

 if (ParallelDescriptor::IOProcessor()) 
  std::cout << "do_the_advance timeSEM= " << timeSEM << 
   " dtSEM= " << dtSEM << '\n';

 very_last_sweep=0;

 advance_status=1; // 1=success 0=failure

 std::string local_caller_string="do_the_advance";

 int finest_level = parent->finestLevel();
 const int max_level = parent->maxLevel();

 if (level>0) 
  amrex::Error("level should equal zero in do_the_advance");
 if (finest_level>max_level)
  amrex::Error("max_level or finest_level invalid");
 if (finest_level>max_level_for_use)
  amrex::Error("max_level_for_use or finest_level invalid");
 
 debug_memory();

 double start_advance = ParallelDescriptor::second();

 if ((ns_time_order<1)||(ns_time_order>32))
  amrex::Error("ns_time_order invalid");

 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
   // allocate MDOT_MF (delete it first if it already exists)
   // MDOT_MF has no ghost cells.
   // MDOT_MF is initialized to zero.
  ns_level.allocate_mdot(); 
  if (verbose>0) {
   ns_level.DumpProcNum();
  }

  int bfact_space=parent->Space_blockingFactor(ilev);
  int bfact_grid=parent->Old_blockingFactor(ilev);
  if ((bfact_space<1)||(bfact_space>64))
   amrex::Error("bfact_space out of range");
  if ((ilev>=0)&&(ilev<finest_level)) {
   if (bfact_grid<4)
    amrex::Error("we must have blocking factor at least 4(1)");
  } else if (ilev==finest_level) {
   if (bfact_grid<4)
    amrex::Error("we must have blocking factor at least 4(1)");
  } else
   amrex::Error("ilev invalid");
  ns_level.check_grid_places();
 } // ilev=level ... finest_level

 double after_init = ParallelDescriptor::second();
 if ((verbose>0)||(show_timings==1)) {
  std::fflush(NULL);
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "elapsed time in preliminary allocation " << after_init-
        start_advance << '\n';
  }
  std::fflush(NULL);
  ParallelDescriptor::Barrier();
 }

 // nsteps==0 very first step.
 // in: do_the_advance
 int nsteps=parent->levelSteps(0); 

 SDC_outer_sweeps=0;
 slab_step=0;
 SDC_setup_step();

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.allocate_SDC();
 }

 Number_CellsALL(real_number_of_cells);

 int SDC_outer_sweeps_end=ns_time_order;

 if (lower_slab_time>=0.0) {
  // do nothing
 } else
  amrex::Error("lower_slab_time invalid");

 if (lower_slab_time==0.0)
  SDC_outer_sweeps_end=1;

 for (SDC_outer_sweeps=0;
      ((SDC_outer_sweeps<SDC_outer_sweeps_end)&&
       (advance_status==1));
      SDC_outer_sweeps++) {

  slab_step=0;
  SDC_setup_step();
  interface_touch_flag=1; //do_the_advance

  if ((SDC_outer_sweeps>=0)&&(SDC_outer_sweeps<ns_time_order)) {
   // do nothing
  } else
   amrex::Error("SDC_outer_sweeps invalid");

   // delta_MF=0.0 if SDC_outer_sweeps==0
   // stableF_MF cleared to 0.0
   // spectralF_MF cleared to 0.0 for all entries if
   //  SDC_outer_sweeps==0 and for all except the first entry if
   //  SDC_outer_sweeps>0.
  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   ns_level.init_delta_SDC();
  }

  int slab_step_start=0;
  int slab_step_end=ns_time_order-1;

  if (SDC_outer_sweeps_end>1) {

   if (SDC_outer_sweeps_end==ns_time_order) {
    if (SDC_outer_sweeps==0)
     slab_step_start--; 
    slab_step_end++; 
   } else
    amrex::Error("SDC_outer_sweeps_end invalid");

  } else if (SDC_outer_sweeps_end==1) {

   // do nothing

  } else
   amrex::Error("SDC_outer_sweeps_end invalid");

   //slab_step==-1 is needed to compute grad p^n with spectral accuracy.
   //slab_step==ns_time_order is needed to compute the SDC correction:
   // delta=integral_tn^tnp1  f^spectral dt - deltatn F^stable
  for (slab_step=slab_step_start;
       ((slab_step<=slab_step_end)&&(advance_status==1));
       slab_step++) {

   interface_touch_flag=1; //do_the_advance

   SDC_setup_step();

   int local_num_divu_outer_sweeps=num_divu_outer_sweeps;

   if ((slab_step==-1)||(slab_step==ns_time_order)) {
    local_num_divu_outer_sweeps=1;
   } else if ((slab_step>=0)&&(slab_step<ns_time_order)) {

    if (dt_slab>0.0) {
     // do nothing
    } else
     amrex::Error("dt_slab must be positive");

    local_num_divu_outer_sweeps=num_divu_outer_sweeps;
   } else
    amrex::Error("slab_step invalid");

   if (verbose>0) {
    std::fflush(NULL);
    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "dt_slab= " << dt_slab << '\n';
    }
    std::fflush(NULL);
    ParallelDescriptor::Barrier();
   }

    // ns.num_divu_outer_sweeps
   for (divu_outer_sweeps=0;
        ((divu_outer_sweeps<local_num_divu_outer_sweeps)&&
         (advance_status==1));
        divu_outer_sweeps++) {

    very_last_sweep=0;

    if ((SDC_outer_sweeps+1==SDC_outer_sweeps_end)&&
        (slab_step+1==ns_time_order)&&
        (divu_outer_sweeps+1==num_divu_outer_sweeps))
     very_last_sweep=1;

    // allocate_physics_vars()
    // mdot=0.0
    for (int ilev=finest_level;ilev>=level;ilev--) {
     NavierStokes& ns_level=getLevel(ilev);
     ns_level.prelim_alloc();
    }

    double start_divu_sweep = ParallelDescriptor::second();

    if (1==0) {
     int gridno=0;
     const Box& fabgrid = grids[gridno];
     const int* fablo=fabgrid.loVect();
     const int* fabhi=fabgrid.hiVect();
     const Real* xlo = grid_loc[gridno].lo();
     int interior_only=1;
     MultiFab& S_new=get_new_data(State_Type,slab_step+1);
     FArrayBox& snewfab=S_new[0];
     const Real* dx = geom.CellSize();
     int scomp_plot=STATECOMP_STATES;
     int ncomp_plot=num_state_material*num_materials;
     tecplot_debug(snewfab,xlo,fablo,fabhi,dx,-1,0,
       scomp_plot,ncomp_plot,interior_only);
    }

      // initialize "law of the wall" velocity derived from solid velocity.
      //  or
      // initialize "GNBC" velocity.
      // in: NavierStokes::do_the_advance (prior to nonlinear_advect)
    int renormalize_only=1;
    init_FSI_GHOST_MAC_MF_ALL(renormalize_only,local_caller_string);

    int mass_transfer_active=0;

     // 1. ADVECTION (both Eulerian and Lagrangian materials)
     // 2. IF MDOT <> 0 in previous time step, then there
     //    is expansion/compression
    if ((slab_step>=0)&&(slab_step<ns_time_order)) {

     if (disable_advection==0) {

      nonlinear_advection(local_caller_string);

      if (step_through_data==1) {
       int basestep_debug=nStep();
       parent->writeDEBUG_PlotFile(
	 basestep_debug,
	 SDC_outer_sweeps,
	 slab_step,
	 divu_outer_sweeps);
       std::cout << "press any number then enter: after nonlinear_advection\n";
       std::cout << "timeSEM= " << timeSEM << '\n';
       std::cout << "dtSEM= " << dtSEM << '\n';
       std::cout << "divu_outer_sweeps= " << divu_outer_sweeps << '\n';
       std::cout << "local_num_divu_outer_sweeps= " << 
	       local_num_divu_outer_sweeps << '\n';
       std::cout << "num_divu_outer_sweeps= " << 
	       num_divu_outer_sweeps << '\n';
       std::cout << "slab_step= " << 
	       slab_step << '\n';
       std::cout << "SDC_outer_sweeps= " << 
	       SDC_outer_sweeps << '\n';
       int n_input;
       std::cin >> n_input;
      }

     } else if (disable_advection==1) {
      
      if (enable_spectral==0) {
       //do nothing
      } else if (enable_spectral==1) {

       amrex::Error("expecting enable_spectral==0 if disable_adv==1");

      } else
       amrex::Error("enable_spectral invalid");

     } else
      amrex::Error("disable_advection invalid");

    } else if ((slab_step==-1)||
               (slab_step==ns_time_order)) {
     // do nothing
    } else
     amrex::Error("slab_step invalid");


     // in: NavierStokes::do_the_advance
    allocate_levelset_ALL(1,LEVELPC_MF);

    if ((ns_time_order==1)&&(enable_spectral!=0))
     amrex::Error("(ns_time_order==1)&&(enable_spectral!=0)");
    if ((ns_time_order>=2)&&(enable_spectral!=1))
     amrex::Error("(ns_time_order>=2)&&(enable_spectral!=1)");

    if ((ns_time_order>=2)&&(enable_spectral==1)) {

     if (disable_advection==0) {

      double start_SEMADV_time=ParallelDescriptor::second();

      int source_term=SUB_OP_SDC_LOW_TIME;

      if ((slab_step>=0)&&(slab_step<=ns_time_order)) {

        // SEM_advectALL starts off by using the prev_time_slab data.
       if ((slab_step==0)&&
           (SDC_outer_sweeps>0)&&
           (SDC_outer_sweeps<ns_time_order)) {
        // do nothing: F(t^n) already init when SDC_outer_sweeps==0.
       } else if ((slab_step==0)&&
                  (SDC_outer_sweeps==0)) {
        SEM_advectALL(source_term);
       } else if ((slab_step>0)&&
                  (SDC_outer_sweeps>=0)&&
                  (SDC_outer_sweeps<ns_time_order)) {
        SEM_advectALL(source_term);
       } else
        amrex::Error("slab_step or SDC_outer_sweeps invalid");

       // above: ((slab_step>=0)&&(slab_step<=ns_time_order))
       
      } else if (slab_step==-1) {
       // do nothing
      } else
       amrex::Error("slab_step invalid");

       // SEM_advectALL starts off by using the prev_time_slab data.
      source_term=SUB_OP_SDC_ISCHEME;
      if ((slab_step>=0)&&(slab_step<ns_time_order)) {
       SEM_advectALL(source_term);
      } else if ((slab_step==-1)||(slab_step==ns_time_order)) {
       // do nothing
      } else {
       amrex::Error("slab_step invalid");
      }

      double end_SEMADV_time=ParallelDescriptor::second();
      if ((verbose>0)||(show_timings==1)) {
       if (ParallelDescriptor::IOProcessor()) {
        std::cout << "elapsed time in SEM_advectALL " << 
         end_SEMADV_time-start_SEMADV_time << '\n';
       }
      }

     } else if (disable_advection==1) {
      // do nothing
     } else
      amrex::Error("disable_advection invalid");

    } else if ((ns_time_order==1)&&(enable_spectral==0)) {
     // do nothing
    } else
     amrex::Error("ns_time_order, enable_spectral invalid do the advance");

      // in: NavierStokes::do_the_advance

    debug_memory();

    mass_transfer_active=0;

    if ((slab_step>=0)&&(slab_step<ns_time_order)) {
     if (is_zalesak()) {
      mass_transfer_active=0;
     } else if (!is_zalesak()) {
      if (is_phasechange==1) { //get_user_latent_heat!=0.0?
       mass_transfer_active=1;
      } else if (is_phasechange==0) {
       mass_transfer_active=0;
      } else
       amrex::Error("is_phasechange invalid");
     } else
      amrex::Error("is_zalesak() bust");
    } else if ((slab_step==-1)||
	       (slab_step==ns_time_order)) {
     // do nothing
    } else
     amrex::Error("slab_step invalid");

    Vector<blobclass> blobdata;
    Vector< Vector<Real> > mdot_data;
    Vector< Vector<Real> > mdot_comp_data;
    Vector< Vector<Real> > mdot_data_redistribute;
    Vector< Vector<Real> > mdot_comp_data_redistribute;
    Vector<int> type_flag;

    blobdata.resize(1);
    int color_count=0;

      // 2. If mass transfer
      //    (a) mass transfer rate (nucleation)
      //    (b) slopes/redistance
      //    (c) mass transfer rate (stefan problem) (level_phase_change_rate)
      //    (d) unsplit advection
      //    (e) slopes/redistance
      //    (f) redistribute mass increments
      //    If no mass transfer
      //    (a) slopes/redistance
    if ((slab_step>=0)&&(slab_step<ns_time_order)) {

      if (mass_transfer_active==1) {

       nucleation_code_segment(
	 local_caller_string,
	 color_count,
	 blobdata,
	 mdot_data,
	 mdot_comp_data,
	 mdot_data_redistribute,
	 mdot_comp_data_redistribute,
	 type_flag);

      } else if (mass_transfer_active==0) {

       no_mass_transfer_code_segment(local_caller_string);
	
      } else
       amrex::Error("mass_transfer_active invalid");

      for (int ilev=finest_level;ilev>=level;ilev--) {
       NavierStokes& ns_level=getLevel(ilev);
       int combine_flag=2; // only update if vfrac<VOFTOL
       int hflag=0;
       // combine_idx==-1 => update S_new  
       // combine_idx>=0  => update localMF[combine_idx]
       int combine_idx=-1;  
       int update_flux=0;
       int interface_cond_avail=0;

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

    } else if ((slab_step==-1)||(slab_step==ns_time_order)) {

      interface_touch_flag=1; //do_the_advance

      int init_vof_prev_time=0;
       //output:SLOPE_RECON_MF
      VOF_Recon_ALL(
	 local_caller_string, //do_the_advance
         cur_time_slab,
         RECON_UPDATE_NULL,init_vof_prev_time);

    } else
      amrex::Error("slab_step invalid");


      // velocity and pressure
    avgDownALL(State_Type,STATECOMP_VEL,STATE_NCOMP_VEL+STATE_NCOMP_PRES,1);

     // den,temp,species, ....
    avgDownALL(State_Type,STATECOMP_STATES, 
               num_state_material*num_materials,1);  

    if ((num_materials_compressible>=1)&&
        (num_materials_compressible<=num_materials)) {
     avgDownALL_refine_density();
    } else if (num_materials_compressible==0) {
     // do nothing
    } else
     amrex::Error("num_materials_compressible invalid in do_the_advance");

    debug_memory();

    double start_phys_time=ParallelDescriptor::second();

    if ((slab_step>=0)&&(slab_step<ns_time_order)) {

      if (mass_transfer_active==1) {

       phase_change_code_segment(
	local_caller_string,
	color_count,
	blobdata);

      } else if (mass_transfer_active==0) {
       // do nothing
      } else
       amrex::Error("mass_transfer_active invalid");

    } else if ((slab_step==-1)||(slab_step==ns_time_order)) {
      // do nothing
    } else
      amrex::Error("slab_step invalid");

    if (1==0) {
     int basestep_debug=nStep();
     parent->writeDEBUG_PlotFile(
	basestep_debug,
	SDC_outer_sweeps,
	slab_step,
	divu_outer_sweeps);
     std::cout << "press any number then enter: before make_physics_varsALL\n";
     int n_input;
     std::cin >> n_input;
    }

     // initialize "law of the wall" velocity derived from solid velocity.
     //  or
     // initialize "GNBC" velocity.
     // in: NavierStokes::do_the_advance (prior to viscous force step, and
     //  after reinitialization)
    renormalize_only=1;
    init_FSI_GHOST_MAC_MF_ALL(renormalize_only,local_caller_string);

// At this stage, variables are not scaled, so FACECOMP_FACEVEL component (c++)
// will have to be scaled later.
    debug_memory();
    if (is_zalesak()) {
     make_physics_varsALL(SOLVETYPE_INITPROJ,local_caller_string); 
    } else {
     make_physics_varsALL(SOLVETYPE_PRES,local_caller_string); 
    }
    delete_array(CELLTENSOR_MF);
    delete_array(FACETENSOR_MF);

    if (1==0) {
     int basestep_debug=nStep();
     parent->writeDEBUG_PlotFile(
       basestep_debug,
       SDC_outer_sweeps,
       slab_step,
       divu_outer_sweeps);
     std::cout << "press any number then enter: after make_physics_varsALL\n";
     int n_input;
     std::cin >> n_input;
    }

    double start_velocity_diff = ParallelDescriptor::second();
    if ((verbose>0)||(show_timings==1)) {
      if (ParallelDescriptor::IOProcessor()) {
       std::cout << "elapsed time in phase change, make_phys_vars " << 
        start_velocity_diff-start_phys_time << '\n';
       std::cout << 
        "nonlinear advective terms, phase change (divu_sweep,cpu time):  " << 
        divu_outer_sweeps << ' ' <<
        start_velocity_diff-start_divu_sweep << '\n';
      }
    }

    if (verbose>0) {
      if (ParallelDescriptor::IOProcessor())
       std::cout << "velocity diffusion \n";
    }
    start_velocity_diff = ParallelDescriptor::second();

    debug_memory();


    if ((slab_step>=0)&&(slab_step<ns_time_order)) {

      //  unew^{f}=
      // (i) unew^{f} in non-solid regions
      // (ii) u^{f,save} + (unew^{c}-u^{c,save})^{c->f} in spectral 
      //      regions.
      //      (u^{c,save} = *localMF[ADVECT_REGISTER_MF])
      //      (u^{f,save} = *localMF[ADVECT_REGISTER_FACE_MF+dir])
      // (iii) usolid in solid regions
      if (is_zalesak()) {
       advance_MAC_velocity(SOLVETYPE_INITPROJ);
      } else {
       advance_MAC_velocity(SOLVETYPE_PRES);
      }

      // TENSOR ADVECTION (non-Newtonian materials)
      // second half of D^{upside down triangle}/Dt
      // extrapolates Q at the end.
      tensor_advection_updateALL();

      if (step_through_data==1) {
       int basestep_debug=nStep();
       parent->writeDEBUG_PlotFile(
	 basestep_debug,
	 SDC_outer_sweeps,
	 slab_step,
	 divu_outer_sweeps);
       std::cout << "press any number then enter: after tensor_advection_updateALL\n";
       std::cout << "timeSEM= " << timeSEM << '\n';
       std::cout << "dtSEM= " << dtSEM << '\n';
       std::cout << "divu_outer_sweeps= " << divu_outer_sweeps << '\n';
       std::cout << "local_num_divu_outer_sweeps= " << 
	       local_num_divu_outer_sweeps << '\n';
       std::cout << "num_divu_outer_sweeps= " << 
	       num_divu_outer_sweeps << '\n';
       std::cout << "slab_step= " << 
	       slab_step << '\n';
       std::cout << "SDC_outer_sweeps= " << 
	       SDC_outer_sweeps << '\n';
       int n_input;
       std::cin >> n_input;
      }

      for (int ilev=finest_level;ilev>=level;ilev--) {
       NavierStokes& ns_level=getLevel(ilev);
       // delete ADVECT_REGISTER_FACE_MF and ADVECT_REGISTER_MF
       ns_level.delete_advect_vars();
       ns_level.delete_transport_vars();
       // initialize ADVECT_REGISTER_FACE_MF and ADVECT_REGISTER_MF
       // delete_advect_vars() called twice in NavierStokes::do_the_advance
       ns_level.prepare_advect_vars(cur_time_slab);
      } // ilev=finest_level ... level


      for (int ilev=finest_level;ilev>=level;ilev--) {
       NavierStokes& ns_level=getLevel(ilev);
       int combine_flag=2; // update F_m=0 cells only.
       int hflag=0;
       // combine_idx==-1 => update S_new  
       // combine_idx>=0  => update localMF[combine_idx]
       int combine_idx=-1;  
       int update_flux=0;
       int interface_cond_avail=0;

       ns_level.combine_state_variable(
        SOLVETYPE_HEAT,
        combine_idx,
        combine_flag, //combine_flag==2 (Fm=0 cells only)
        hflag,
        update_flux,
        interface_cond_avail); 
       for (int ns=0;ns<num_species_var;ns++) {
        ns_level.combine_state_variable(
         SOLVETYPE_SPEC+ns,
         combine_idx,
         combine_flag, //combine_flag==2 (Fm=0 cells only)
         hflag,
         update_flux,
         interface_cond_avail); 
       }
       ns_level.combine_state_variable(
        SOLVETYPE_VISC, //cell centered velocity
        combine_idx,
        combine_flag, //combine_flag==2 (Fm=0 cells only)
        hflag,
        update_flux,
        interface_cond_avail);

       update_flux=1;
       ns_level.combine_state_variable(
        SOLVETYPE_PRES, //flux var=mac velocity
        combine_idx,
        combine_flag,
        hflag,
        update_flux,
        interface_cond_avail);
      } // ilev = finest_level ... level

      if (verbose>0) {
       if (ParallelDescriptor::IOProcessor()) {
        std::cout << "after increment_face_velocityALL \n";
        std::cout << "after combine \n";
        std::cout << "after update_total_energy \n";
       }
      }
  
      if (mass_transfer_active==1) {

       // 1. modifies localMF[JUMP_STRENGTH_MF] (size 2 * num_interfaces)
       //    a) fort_tagexpansion
       //    b) fort_distributeexpansion
       //    c) fort_clearexpansion
       //    d) fort_initjumpterm ( modifies localMF[MDOT_MF] )
       //    e) fort_getcolorsum (twice)
       phase_change_redistributeALL();

       delete_array(JUMP_STRENGTH_MF);
   
      } else if (mass_transfer_active==0) {
       // do nothing
      } else
       amrex::Error("mass_transfer_active invalid");

      if (verbose>0) {
       if (ParallelDescriptor::IOProcessor()) {
        std::cout << "after phasechange, level_phase_change_redist \n";
       }
      }
      debug_memory();

      if (is_zalesak()==1) {

       for (int ilev=finest_level;ilev>=level;ilev--) {
        NavierStokes& ns_level=getLevel(ilev);
        // velocity scaled by global_velocity_scale
        // init cell center gas/liquid velocity 
        ns_level.zalesakVEL();  
       } // ilev=finest_level ... level

        // unew^{f} = unew^{c->f}
       int operation_flag=OP_UNEW_CELL_TO_MAC;
       int idx_velcell=-1;

       Real beta_local=0.0;
       Vector<blobclass> local_blobdata;
       increment_face_velocityALL(
         operation_flag,
         SOLVETYPE_INITPROJ,
         idx_velcell,beta_local,local_blobdata);

      } else if (is_zalesak()==0) {

       for (int ilev=finest_level;ilev>=level;ilev--) {
        NavierStokes& ns_level=getLevel(ilev);
        int combine_flag=2; // combine if vfrac<VOFTOL
        int hflag=0; // inhomogeneous option
        int combine_idx=-1;  // update state variables
        int update_flux=0;
        int interface_cond_avail=0;

	  // declared in: Diffusion.cpp
        ns_level.combine_state_variable(
         SOLVETYPE_VISC,  //cell centered velocity
         combine_idx,
         combine_flag,
         hflag,
         update_flux, 
         interface_cond_avail);
        update_flux=1;
        ns_level.combine_state_variable(
         SOLVETYPE_PRES, //mac velocity
         combine_idx,
         combine_flag,
         hflag,
         update_flux,
         interface_cond_avail);
       } // ilev = finest_level ... level

        // T_advect_MF=new temperature
       int alloc_flag=1;
       alloc_DTDtALL(alloc_flag);

       for (FSI_outer_sweeps=0;FSI_outer_sweeps<num_FSI_outer_sweeps; 
            FSI_outer_sweeps++) {

        //project_to_rigid_velocityALL is declared in NavierStokes2.cpp
        //fort_project_to_rigid_velocity is declared in LEVELSET_3D.F90
        //overwrite the MAC and CELL velocity.
        project_to_rigid_velocityALL();

        if (num_FSI_outer_sweeps==1) {
         //do nothing
        } else if ((num_FSI_outer_sweeps>=2)&&
                   (num_FSI_outer_sweeps<=num_materials)) {

         for (int ilev=finest_level;ilev>=level;ilev--) {
          NavierStokes& ns_level=getLevel(ilev);
          ns_level.level_init_elasticmask_and_elasticmaskpart();

          ns_level.avgDownEdge_localMF(FACE_VAR_MF,FACECOMP_ELASTICMASK,1,0,
           AMREX_SPACEDIM,LOW_ORDER_AVGDOWN,local_caller_string);
          ns_level.avgDownEdge_localMF(FACE_VAR_MF,FACECOMP_ELASTICMASKPART,1,0,
           AMREX_SPACEDIM,LOW_ORDER_AVGDOWN,local_caller_string);

          ns_level.avgDownEdge_localMF(FACE_VAR_MF,FACECOMP_FACEDEN,1,0,
 	   AMREX_SPACEDIM,LOW_ORDER_AVGDOWN,local_caller_string);
          ns_level.avgDownEdge_localMF(FACE_VAR_MF,FACECOMP_FACEVISC,1,0,
           AMREX_SPACEDIM,LOW_ORDER_AVGDOWN,local_caller_string);

          ns_level.manage_FSI_data(); 
         }
        } else
         amrex::Error("num_FSI_outer_sweeps invalid");

        // 4. Backwards Euler building block: VISCOSITY, thermal diffusion,
        //    species diffusion
        //   a. hoop stress
        //   b. boussinesq approximation:
        //      rho Du/Dt=-grad p + rho g zhat +
        //                 rho Omega^2 r rhat
        //      rho0 Du/Dt=-grad p + rho g zhat +
        //                 rho Omega^2 r rhat
        //      rho/rho0=(1+beta(T-T0))
        //      Du/Dt=-grad p/rho0+g zhat + g zhat beta (T-T0)+
        //            Omega^{2} r rhat +
        //            Omega^{2} r rhat beta (T-T0)
        //      (beta<0 here)
        //      
        //   c. coriolis effect
        //   d. viscous force
        //   e. viscoelastic force
        //   f. FSI force
        //   g. momentum force
        //   h. Marangoni force 
        //
        veldiffuseALL();  

        debug_memory();
  
        double end_velocity_diff = ParallelDescriptor::second();
        if ((verbose>0)||(show_timings==1)) {
         if (ParallelDescriptor::IOProcessor()) {
          std::cout << "veldiffuse time " << 
           end_velocity_diff-start_velocity_diff << '\n';
         }
        }

        //The viscous solver might have slightly perturbed the 
        //velocity in the rigid solid regions, so the velocity
        //in these regions must be re-prescribed.
        for (int ilev=finest_level;ilev>=level;ilev--) {
         NavierStokes& ns_level=getLevel(ilev);
         int combine_flag=2; // combine if vfrac<VOFTOL
         int hflag=0; // inhomogeneous option
         int combine_idx=-1;  // update state variables
         int update_flux=0;
         int interface_cond_avail=0;

	  // declared in: Diffusion.cpp
         ns_level.combine_state_variable(
          SOLVETYPE_VISC,  //cell centered velocity
          combine_idx,
          combine_flag, //=2 (combine if vfrac<VOFTOL)
          hflag,
          update_flux, 
          interface_cond_avail);

         update_flux=1;
         ns_level.combine_state_variable(
          SOLVETYPE_PRES, //mac velocity
          combine_idx,
          combine_flag,
          hflag,
          update_flux,
          interface_cond_avail);
        } // ilev = finest_level ... level


        if (1==0) {
           // S_new is level 0 data
         MultiFab& S_new=get_new_data(State_Type,slab_step+1);
	 // data file name "VISCSOLVE<stuff>.plt"
	 // after the viscous solve, but before the pressure projection.
	 // cell data in the fluid, next to the solid, should "make sense"
	 // xvel,yvel,zvel,pressure,(density, temperature) x num_materials,
	 // (VFRAC,centroid) x num_materials, error indicator
         writeSanityCheckData(
          "VISCSOLVE",
          "in: NavierStokes::do_the_advance, State_Type after veldiffuseALL", 
          local_caller_string,
          State_Type+GET_NEW_DATA_OFFSET, //tower_mf_id
          S_new.nComp(),
          -1, // data_mf==-1
          State_Type, //state_type_mf==State_Type
          -1,  // data_dir==-1
          parent->levelSteps(0)); 
        }

        if (1==0) {
         int basestep_debug=nStep();
         parent->writeDEBUG_PlotFile(
	  basestep_debug,
	  SDC_outer_sweeps,
	  slab_step,
	  divu_outer_sweeps);
        }

        if (FSI_outer_sweeps==0) {

         //NavierStokes::Mass_Energy_Sources_SinksALL declared in 
	 //this file: NavierStokes3.cpp
         Mass_Energy_Sources_SinksALL();

         // DTDt_MF=T_new - T_advect_MF
         alloc_flag=2;
         alloc_DTDtALL(alloc_flag);

         int is_any_lowmach=0;

         for (int im_low=0;im_low<num_materials;im_low++) {
          if (ns_is_rigid(im_low)==1) {
           // do nothing
          } else if (ns_is_rigid(im_low)==0) {
           if ((material_type[im_low]==0)&&
               (material_type_lowmach[im_low]>0)&&
               (material_type_lowmach[im_low]<999)) {
            is_any_lowmach=1;
           } else if ((material_type[im_low]==0)&&
                      (material_type_lowmach[im_low]==0)) {
            // do nothing
           } else if ((material_type[im_low]>0)&&
                      (material_type[im_low]<999)&&
                      (material_type_lowmach[im_low]==
    		       material_type[im_low])) {
            // do nothing
           } else
            amrex::Error("material_type or material_type_lowmach invalid");
          } else
           amrex::Error("ns_is_rigid invalid");
         }  //im_low=0..num_materials-1

         if (is_any_lowmach==1) {

          Vector<blobclass> local_blobdata;
  
          //local_mdot_data, local_mdot_comp_data,
          //local_mdot_data_redistribute, 
          //local_mdot_comp_data_redistribute,
          //are not used; they are placeholders.
          Vector< Vector<Real> > local_mdot_data;
          Vector< Vector<Real> > local_mdot_comp_data;
          Vector< Vector<Real> > local_mdot_data_redistribute;
          Vector< Vector<Real> > local_mdot_comp_data_redistribute;
          Vector<int> local_type_flag;

          int local_color_count=0;
          int coarsest_level=0;
          int idx_mdot=-1; //idx_mdot==-1 => do not collect auxiliary data.
          int local_tessellate=3;
          int operation_flag=OP_GATHER_MDOT; // allocate TYPE_MF,COLOR_MF

          // for each blob, find sum_{F>=1/2} pressure * vol and
 	  // sum_{F>=1/2} vol.
          // calling from: NavierStokes::do_the_advance()
          ColorSumALL(
           operation_flag, // =OP_GATHER_MDOT
           local_tessellate, //=3
           coarsest_level,
           local_color_count,
           TYPE_MF,COLOR_MF,
           idx_mdot,
           idx_mdot,
           local_type_flag,
           local_blobdata,
           local_mdot_data,
           local_mdot_comp_data,
           local_mdot_data_redistribute,
           local_mdot_comp_data_redistribute 
          );
          ParallelDescriptor::Barrier();

          if (color_count!=blobdata.size())
           amrex::Error("color_count!=blobdata.size()");

          // increment MDOF_MF
          LowMachDIVUALL(
           coarsest_level,
           local_color_count,
           TYPE_MF,
           COLOR_MF,
           local_type_flag, 
           local_blobdata);

          delete_array(TYPE_MF);
          delete_array(COLOR_MF);

         } else if (is_any_lowmach==0) {
          // do nothing
         } else
          amrex::Error("is_any_lowmach invalid");

         // delete T_advect_MF, DTDt_MF
         alloc_flag=0;
         alloc_DTDtALL(alloc_flag);

        } else if ((FSI_outer_sweeps>=1)&&
                   (FSI_outer_sweeps<num_FSI_outer_sweeps)) {
         //do nothing
        } else
         amrex::Error("FSI_outer_sweeps invalid");

        double intermediate_time = ParallelDescriptor::second();

         // 5. PRESSURE PROJECTION 
         //    a. gravity
         //    b. surface tension
         //    c. pressure gradient
        if (disable_pressure_solve==0) {

         pressure_gradient_code_segment(local_caller_string);

        } else if (disable_pressure_solve==1) {
         // do nothing
        } else
         amrex::Error("disable_pressure_solve invalid");


        double end_pressure_solve = ParallelDescriptor::second();

        if ((verbose>0)||(show_timings==1)) {
         if (ParallelDescriptor::IOProcessor()) {
          std::cout << "FSI_outer_sweeps=" << FSI_outer_sweeps << '\n';
          std::cout << "pressure solve time " << end_pressure_solve-
            intermediate_time << '\n';
          std::cout << "number of cells in the pressure solve " <<
             real_number_of_cells << '\n';
         }
        }

       } //FSI_outer_sweeps=0 ... num_FSI_outer_sweeps-1

       FSI_outer_sweeps=0;


        // in: do_the_advance
        // 1. prescribe solid temperature, velocity, and geometry where
        //    appropriate.
        // 2. extend level set functions into the solid.
       renormalize_only=0;
       int local_truncate=0;
       int update_particles=0;
       prescribe_solid_geometryALL(cur_time_slab,renormalize_only,
        local_truncate,local_caller_string,update_particles);

       for (int ilev=finest_level;ilev>=level;ilev--) {
        NavierStokes& ns_level=getLevel(ilev);
        int combine_flag=2; // combine if vfrac<VOFTOL
        int hflag=0;
        int combine_idx=-1;  // update state variables
        int update_flux=0;
        int interface_cond_avail=0;
        ns_level.combine_state_variable(
         SOLVETYPE_VISC, // cell centered velocity
         combine_idx,
         combine_flag, // =2 (combine if vfrac<VOFTOL)
         hflag,
         update_flux,
         interface_cond_avail);

        update_flux=1;
        ns_level.combine_state_variable(
         SOLVETYPE_PRES, //flux var = mac velocity
	 combine_idx,
         combine_flag,
         hflag,
         update_flux,
         interface_cond_avail);
       } // for (int ilev=finest_level;ilev>=level;ilev--) 
       debug_memory();

      } else
       amrex::Error("is_zalesak invalid");

      for (int ilev=finest_level;ilev>=level;ilev--) {
       NavierStokes& ns_level=getLevel(ilev);
       // delete ADVECT_REGISTER_FACE_MF and ADVECT_REGISTER_MF
       ns_level.delete_advect_vars();
      }

      if (mass_transfer_active==1) {
       delete_array(SATURATION_TEMP_MF);
      } else if (mass_transfer_active==0) {
       // do nothing
      } else
       amrex::Error("mass_transfer_active invalid");

      debug_memory();

      if (very_last_sweep==0) {

        // fixed_dt==0.0 if dt not prescribed.
        // fixed_dt>0.0 if dt prescribed.
       Real local_fixed_dt;
       if (nsteps==0) {
        local_fixed_dt=fixed_dt_init;
       } else if (nsteps>0) {
        local_fixed_dt=fixed_dt;
       } else {
        local_fixed_dt=0.0;
        amrex::Error("nsteps invalid");
       }

       Real dt_predict=estTimeStep(local_fixed_dt,local_caller_string);
       Real dt_predict_max=dt_predict;
       Real dt_predict_min=dt_predict;
       ParallelDescriptor::ReduceRealMax(dt_predict_max);
       ParallelDescriptor::ReduceRealMin(dt_predict_min);
       Real dt_error=dt_predict_max-dt_predict_min;
       if ((dt_error>=0.0)&&
           (dt_predict_min>0.0)&&
           (dt_predict_max>0.0)) {
	Real dt_tol=dt_predict_max*CPP_EPS_13_6;
        if (dt_predict_min<1.0)
         dt_tol=CPP_EPS_13_6;
	if (dt_error<=dt_tol) {
 	 // do nothing
	} else 
	 amrex::Error("dt_error>dt_tol");
       } else
        amrex::Error("dt_error, dt_predict_min, or dt_predict_max bad");

       if (dt_predict_min<=0.5*delta_slab_time)
        advance_status=0;

      } else if (very_last_sweep==1) {
       // do nothing
      } else
       amrex::Error("very_last_sweep invalid");

      // get rid of uninit in the boundaries of the state variables.
      for (int ilev=level;ilev<=finest_level;ilev++) {
       NavierStokes& ns_level=getLevel(ilev);
       ns_level.init_boundary();
      }

    } else if ((slab_step==-1)||(slab_step==ns_time_order)) {

      // do nothing
      
    } else
      amrex::Error("slab_step invalid");

    if ((ns_time_order>=2)&&(advance_status==1)) {

      if (enable_spectral==1) {
       // do nothing
      } else
       amrex::Error("enable_spectral invalid 2");

      if ((slab_step>=-1)&&(slab_step<ns_time_order)) {

       if (divu_outer_sweeps+1==local_num_divu_outer_sweeps) {

        int update_spectralF=1;

        // grad p (MAC)
        int update_stableF=1;
        if (slab_step==-1)
         update_stableF=0;

        Vector<int> scomp;  
        Vector<int> ncomp;  
        int ncomp_check;
        int state_index;

	 //num_materials_combine=1
        get_mm_scomp_solver(
         1,
         SOLVETYPE_PRES,
         state_index,
         scomp,
         ncomp,
         ncomp_check);

        int nsolve=1;

        if (state_index!=State_Type)
         amrex::Error("state_index invalid");
        if (ncomp_check!=nsolve)
         amrex::Error("ncomp_check invalid");

         // data at time = cur_time_slab
        getState_localMF_listALL(
          PRESPC2_MF,1,
          state_index,
          scomp,
          ncomp);

	 // currently in: NavierStokes::do_the_advance
         // HOfab=grad p (MAC) if (update_spectralF==1)
         // LOfab=grad p (MAC) if (update_stableF==1)
	 // NavierStokes::update_SEM_forcesALL declared in MacProj.cpp
	 // NavierStokes::update_SEM_forces declared in MacProj.cpp
         // fort_updatesemforce declared in GODUNOV_3D.F90
        update_SEM_forcesALL(
         SOLVETYPE_PRES,
         PRESPC2_MF,
         update_spectralF,update_stableF);

        // end: grad p (MAC)

	 //num_materials_combine=1
        get_mm_scomp_solver(
         1,
         SOLVETYPE_VISC,
         state_index,
         scomp,
         ncomp,
         ncomp_check);

        nsolve=AMREX_SPACEDIM;

        if (state_index!=State_Type)
         amrex::Error("state_index invalid");
        if (ncomp_check!=nsolve)
         amrex::Error("ncomp_check invalid");

          // data at time = cur_time_slab
        getState_localMF_listALL(
          REGISTER_MARK_MF,1,
          state_index,
          scomp,
          ncomp);

         // -div(k grad T)
        update_stableF=0;

	 //num_materials_combine=1
        get_mm_scomp_solver(
         1,
         SOLVETYPE_HEAT,
         state_index,
         scomp,
         ncomp,
         ncomp_check);

        nsolve=1;

        if (state_index!=State_Type)
         amrex::Error("state_index invalid");
        if (ncomp_check!=nsolve)
         amrex::Error("ncomp_check invalid");

          // data at time = cur_time_slab
        getState_localMF_listALL(
          BOUSSINESQ_TEMP_MF,1,
          state_index,
          scomp,
          ncomp);

	// currently in: NavierStokes::do_the_advance
        // HOfab=-div(k grad T)
	// NavierStokes::update_SEM_forcesALL declared in MacProj.cpp
	// NavierStokes::update_SEM_forces declared in MacProj.cpp
        // fort_updatesemforce declared in GODUNOV_3D.F90
        if (enable_spectral==1) {
         update_SEM_forcesALL(SOLVETYPE_HEAT,
          BOUSSINESQ_TEMP_MF,
          update_spectralF,update_stableF);
        } else
         amrex::Error("enable_spectral invalid");

        // HOfab=-div(2 mu D)-HOOP_FORCE_MARK_MF
        // calls: fort_updatesemforce declared in GODUNOV_3D.F90
        update_stableF=0;

	 //num_materials_combine=1
        get_mm_scomp_solver(
         1,
         SOLVETYPE_VISC,
         state_index,
         scomp,
         ncomp,
         ncomp_check);

        nsolve=AMREX_SPACEDIM;

        if (state_index!=State_Type)
         amrex::Error("state_index invalid");
        if (ncomp_check!=nsolve)
         amrex::Error("ncomp_check invalid");

        allocate_array(1,nsolve,-1,HOOP_FORCE_MARK_MF);
         // update_state==OP_HOOP_BOUSSINESQ_IMPLICIT:
         //  unp1(1)=unp1(1)/(one+param2*hoop_force_coef)-dt|g|beta(T-T0)
         // update_state==OP_HOOP_BOUSSINESQ_EXPLICIT:
         //  unp1(1)=unp1(1)-param2*hoop_force_coef*un(1)-dt|g|beta(T-T0)
	 // force=rho(unp1-un)/dt
        int update_state=OP_HOOP_BOUSSINESQ_EXPLICIT;
        diffuse_hoopALL(REGISTER_MARK_MF,BOUSSINESQ_TEMP_MF,
         HOOP_FORCE_MARK_MF,update_state);

	// currently in: NavierStokes::do_the_advance
        // HOfab=-div(2 mu D)-HOOP_FORCE_MARK_MF
	// NavierStokes::update_SEM_forcesALL declared in MacProj.cpp
	// NavierStokes::update_SEM_forces declared in MacProj.cpp
        // fort_updatesemforce declared in GODUNOV_3D.F90
        if (enable_spectral==1) {
         update_SEM_forcesALL(SOLVETYPE_VISC,
          REGISTER_MARK_MF,
          update_spectralF,update_stableF);
        } else
         amrex::Error("enable_spectral invalid");

        delete_array(PRESPC2_MF); // pressure
        delete_array(REGISTER_MARK_MF); // velocity
        delete_array(BOUSSINESQ_TEMP_MF); // temperature
        delete_array(HOOP_FORCE_MARK_MF);

       } else if ((divu_outer_sweeps>=0)&&
                  (divu_outer_sweeps+1<local_num_divu_outer_sweeps)) {
         // do nothing
       } else
        amrex::Error("divu_outer_sweeps invalid do_the_advance");
    
      } else if (slab_step==ns_time_order) {

       if (divu_outer_sweeps==0) {

        if ((SDC_outer_sweeps>=0)&&
            (SDC_outer_sweeps<ns_time_order)) {

         // delta=integral_tn^tnp1  f^spectral dt - deltatn F^stable
         for (int ilev=finest_level;ilev>=level;ilev--) {
          NavierStokes& ns_level=getLevel(ilev);
          ns_level.init_splitting_force_SDC();
         }

	} else
	 amrex::Error("SDC_outer_sweeps invalid");
       } else
	amrex::Error("divu_outer_sweeps invalid");

      } else
       amrex::Error("slab_step invalid");

    } else if ((ns_time_order==1)||(advance_status==0)) {
      // do nothing
    } else
      amrex::Error("ns_time_order or advance_status invalid: do_the..");


    if (step_through_data==1) {
     int basestep_debug=nStep();
     parent->writeDEBUG_PlotFile(
       basestep_debug,
       SDC_outer_sweeps,
       slab_step,
       divu_outer_sweeps);
     std::cout << "press any number then enter: just before close of divu_outer_sweeps loop\n";
     std::cout << "timeSEM= " << timeSEM << '\n';
     std::cout << "dtSEM= " << dtSEM << '\n';
     std::cout << "divu_outer_sweeps= " << divu_outer_sweeps << '\n';
     std::cout << "local_num_divu_outer_sweeps= " << 
             local_num_divu_outer_sweeps << '\n';
     std::cout << "num_divu_outer_sweeps= " << 
             num_divu_outer_sweeps << '\n';
     std::cout << "slab_step= " << 
             slab_step << '\n';
     std::cout << "SDC_outer_sweeps= " << 
             SDC_outer_sweeps << '\n';
     int n_input;
     std::cin >> n_input;
    }

   } // divu_outer_sweeps loop

  } // slab_step loop

 } // SDC_outer_sweeps loop

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.deallocate_SDC();
 }

 if (advance_status==0) {

  // do nothing (failure)

 } else if (advance_status==1) { // success

  slab_step=ns_time_order-1;

  // declared in: MacProj.cpp
  // CELL_SOUND_MF contains: (i) coefficient (1/(rho c^2 dt^2))
  //                         (ii) p_advect
  // MultiFab& DIV_new=get_new_data(DIV_Type,slab_step+1);
  // if compressible: DIV_new=-dt(pnew-padv)/(rho c^2 dt^2)+MDOT_MF dt/vol=
  //                          -(pnew-padv)/(rho c^2 dt)+MDOT_MF dt/vol
  //                          
  // if incompressible: DIV_new=MDOT_MF dt/vol
  ADVECT_DIV_ALL();

  if (visual_divergence_plot_int>0) {

   SDC_outer_sweeps=ns_time_order-1;
   divu_outer_sweeps=num_divu_outer_sweeps-1;

   int ratio=(nsteps+1)/visual_divergence_plot_int;
   ratio=ratio*visual_divergence_plot_int;
   if (ratio==nsteps+1) {

    // declared in: MacProj.cpp
    int idx_source=-1;
    int scomp_src=0;
    int idx_mask=-1;
    getStateDIV_ALL(idx_source,scomp_src,MACDIV_MF,idx_mask);
    if (localMF[MACDIV_MF]->nComp()!=1)
     amrex::Error("localMF[MACDIV_MF]->nComp() invalid");
    if (localMF[MACDIV_MF]->nGrow()!=1)
     amrex::Error("localMF[MACDIV_MF]->nGrow() invalid");

     //MACDIV<stuff>.plt (visit can open binary tecplot files)
    writeSanityCheckData(
      "MACDIV",
      "MACDIV_MF: actual div u",
      local_caller_string,
      MACDIV_MF, //tower_mf_id
      localMF[MACDIV_MF]->nComp(), 
      MACDIV_MF,
      -1,  // State_Type==-1 
      -1,  // data_dir==-1 (cell centered)
      parent->levelSteps(0)); 

    MultiFab& DIV_new=get_new_data(DIV_Type,slab_step+1);
    //DIV_Type<stuff>.plt (visit can open binary tecplot files)
    writeSanityCheckData(
      "DIV_Type",
      "DIV_Type: -(pnew-padv)/(rho c^2 dt)+MDOT_MF dt/vol",
      local_caller_string,
      DIV_Type+GET_NEW_DATA_OFFSET,  // tower_mf_id
      DIV_new.nComp(), 
      -1, //data_mf==-1
      DIV_Type,  // state_type_mf==DIV_TYPE
      -1,  // data_dir==-1 (cell centered)
      parent->levelSteps(0)); 

    delete_array(MACDIV_MF);
   }

  } else if (visual_divergence_plot_int==0) {
   // do nothing
  } else
   amrex::Error("visual_divergence_plot_int invalid");

  debug_memory();

  double end_advance = ParallelDescriptor::second();
  start_advance=end_advance-start_advance;
  total_advance_time+=start_advance;

  if ((verbose>0)||(show_timings==1)) {
   if (ParallelDescriptor::IOProcessor()) {
    std::cout << "advance time " << start_advance << '\n';
    std::cout << "total advance time " << total_advance_time << '\n';
   }
  }
  if (1==0) {
   std::cout << "PROC= " << ParallelDescriptor::MyProc() << 
	   " advance time " << start_advance << '\n';
   std::cout << "PROC= " << ParallelDescriptor::MyProc() <<
	   " total advance time " << total_advance_time << '\n';
  }

  ParallelDescriptor::Barrier();
 } else 
  amrex::Error("advance status invalid");

}  // end subroutine do_the_advance

void NavierStokes::push_stack(
       Vector<int>& stackdata,
       int& stackptr,
       int& data) {

 int i_size=stackdata.size();
 if (i_size>=1) {
  //do nothing
 } else
  amrex::Error("i_size invalid");

 if (stackptr+1>i_size-1)
  amrex::Error("stackptr invalid");

 stackptr++;
 stackdata[stackptr]=data;
}

void NavierStokes::pop_stack(
       Vector<int>& stackdata,
       int& stackptr,
       int& data) {

 int i_size=stackdata.size();
 if (i_size>=1) {
  //do nothing
 } else
  amrex::Error("i_size invalid");

 if (stackptr>i_size-1)
  amrex::Error("stackptr invalid");

 if (stackptr<0)
  amrex::Error("stack is empty");

 data=stackdata[stackptr];

 stackptr--;

}

void NavierStokes::cross_check(
  Vector<int>& levelcolormap,
  Vector<int>& stackdata,
  Vector< Vector<int> >& grid_color,
  int i) {

 int n_unique=grid_color.size();
 int n_assoc=0;
 for (int i_unique=0;i_unique<n_unique;i_unique++) { 
  n_assoc+=grid_color[i_unique].size();
 }
 if (n_assoc==stackdata.size()) {
  //do nothing
 } else
  amrex::Error("stackdata.size() invalid");

 int stackptr=-1;

 if (levelcolormap.size()==grid_color.size()) {
  //do nothing
 } else
  amrex::Error("levelcolormap.size() invalid"); 

 for (int j=0;j<grid_color[i].size();j++) {
  push_stack(stackdata,stackptr,grid_color[i][j]);
 }

 while (stackptr>=0) {
  int j;
  pop_stack(stackdata,stackptr,j); 

  if (levelcolormap[j]==0) {
   levelcolormap[j]=levelcolormap[i];
   if ((levelcolormap[j]<1)||(levelcolormap[j]>levelcolormap.size()))
    amrex::Error("levelcolormap[j] invalid");

   for (int jj=0;jj<grid_color[j].size();jj++) {
    push_stack(stackdata,stackptr,grid_color[j][jj]);
   }

  } else if (levelcolormap[j]!=levelcolormap[i])
   amrex::Error("something wrong in cross_check");
 }  // while stackptr>=0

} // end subroutine cross_check

void NavierStokes::correct_colors(
 int idx_color,int base_level,
 Vector<int> domaincolormap,
 int max_colors_level) {

 int finest_level=parent->finestLevel();

 if ((base_level>=0)&&(base_level<finest_level)) {
  //do nothing
 } else
  amrex::Error("base_level invalid");

 const Real* dx = geom.CellSize();

 MultiFab* colormf=localMF[idx_color];

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(colormf->boxArray().d_numPts());

// COLORING LOOP
#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*colormf,false); mfi.isValid(); ++mfi) {
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

  FArrayBox& colorfab=(*colormf)[mfi];
  int arrsize=domaincolormap.size();

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // in: LEVELSET_3D.F90
  fort_levelrecolor(
   colorfab.dataPtr(),
   ARLIM(colorfab.loVect()),ARLIM(colorfab.hiVect()),
   xlo,dx,
   tilelo,tilehi,
   fablo,fabhi,&bfact,
   domaincolormap.dataPtr(),
   &max_colors_level,
   &level,
   &base_level,&arrsize);
 } // mfi
} //omp
 ns_reconcile_d_num(LOOP_LEVELRECOLOR,"correct_colors");

}  // subroutine correct_colors

// type_flag should have size num_materials.
// type_flag[i]=1 if fluid "i" exists. (note, fictitious solid on the
//  boundaries will show up as existing if ngrow>0)
//
void NavierStokes::assign_colors(
 int& fully_covered,
 int idx_color,int idx_type,
 Vector<int>& colormax,Vector<int> type_flag,
 int zero_diag_flag) {

 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid assign_colors");

 fully_covered=0;

 int ncomp_type=num_materials;
 if (zero_diag_flag==1) {
  ncomp_type=2;
 } else if (zero_diag_flag==0) {
  ncomp_type=num_materials;
 } else
  amrex::Error("zero_diag_flag invalid");

 int typedim=type_flag.size();
 if (typedim!=ncomp_type)
  amrex::Error("typedim invalid");

 colormax[level]=0;
 for (int i=0;i<typedim;i++)
  colormax[level]+=type_flag[i];

 MultiFab* typemf=localMF[idx_type];
 MultiFab* colormf=localMF[idx_color];
 int number_grids=grids.size();

 Vector< Vector<int> > color_per_grid_array;
 Vector<int> max_colors_grid_array;
 color_per_grid_array.resize(thread_class::nthreads);
 max_colors_grid_array.resize(thread_class::nthreads);

 for (int tid=0;tid<thread_class::nthreads;tid++) {
  color_per_grid_array[tid].resize(number_grids);
  for (int i=0;i<number_grids;i++)
   color_per_grid_array[tid][i]=0;
  max_colors_grid_array[tid]=0;
 } // tid

  // mask=tag if not covered by level+1 and at fine-fine ghost cell.
 int ngrowmask=1;
 Real tag=1.0;
 int clear_phys_boundary=2;
 MultiFab* maskmf=maskfiner(ngrowmask,tag,clear_phys_boundary); 

 colormf->setVal(0.0,0,1,1);

  // first pass initialize the colors
  // 1<=color_num<=max_colors_grid
  // 0<=grid_no<num_grids
  // 2nd pass, replace color_num with max_colors_grid*grid_no+color_num
  // no tiling on this loop.
 for (int ipass=0;ipass<2;ipass++) {

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(typemf->boxArray().d_numPts());

// COLORING LOOP
#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(*typemf,false); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const Box& tilegrid = mfi.tilebox();

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   const int gridno = mfi.index();
   const Box& grd=grids[gridno];
   const int* lo = grd.loVect();
   const int* hi = grd.hiVect();

   FArrayBox& maskfab=(*maskmf)[mfi]; // maskmf=1 if not covered
   FArrayBox& typefab=(*typemf)[mfi];
   FArrayBox& colorfab=(*colormf)[mfi];

   int tid_max_colors=tid_current;
   if (ipass==0) {
    // do nothing
   } else if (ipass==1) {
    tid_max_colors=0;
   } else
    amrex::Error("ipass invalid");

   BaseFab<int> ijkstack(grd,AMREX_SPACEDIM);
    // colors are NOT assigned where mask=0
    // colors only assigned where mask=1
    // only valid cells are investigated.
    // in: LEVELSET_3D.F90
   fort_colorfill(
     maskfab.dataPtr(),ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
     typefab.dataPtr(),ARLIM(typefab.loVect()),ARLIM(typefab.hiVect()),
     colorfab.dataPtr(),ARLIM(colorfab.loVect()),ARLIM(colorfab.hiVect()),
     ijkstack.dataPtr(),ARLIM(ijkstack.loVect()),ARLIM(ijkstack.hiVect()),
     lo,hi,
     &ipass,
     &number_grids,
     color_per_grid_array[tid_current].dataPtr(),
     &gridno,
     &max_colors_grid_array[tid_max_colors],
     &typedim);

   if (ipass==0)
    if (color_per_grid_array[tid_current][gridno]>
	max_colors_grid_array[tid_current])
     max_colors_grid_array[tid_current]= 
	     color_per_grid_array[tid_current][gridno];
  } // mfi
} // omp
  ns_reconcile_d_num(LOOP_COLORFILL,"assign_colors");

  if (ipass==0) {

   for (int tid=1;tid<thread_class::nthreads;tid++) {
    if (max_colors_grid_array[tid]>max_colors_grid_array[0])
     max_colors_grid_array[0]=max_colors_grid_array[tid];
    for (int i=0;i<number_grids;i++) {
     if (color_per_grid_array[tid][i]>color_per_grid_array[0][i]) 
      color_per_grid_array[0][i]=color_per_grid_array[tid][i];
    }  // i
   } // tid

  } // ipass==0


  ParallelDescriptor::Barrier();

  if (ipass==0) {
   ParallelDescriptor::ReduceIntMax(max_colors_grid_array[0]);
   for (int i=0;i<number_grids;i++)
    ParallelDescriptor::ReduceIntMax(color_per_grid_array[0][i]);
  } // ipass==0


 } // ipass=0 to 1

   // max_colors_grid_array[0]==0 if this level is completely covered
 if (max_colors_grid_array[0]==0) {
   fully_covered=1;
   if (level==finest_level)
    amrex::Error("max_colors_grid_array[0]==0 on finest_level");
 } else if (max_colors_grid_array[0]>=1) {
   int check_corners=1;
   sync_colors(idx_color,idx_type,
    color_per_grid_array[0],
    colormax,
    max_colors_grid_array[0],
    maskmf,check_corners,
    zero_diag_flag);
 } else
   amrex::Error("max_colors_grid_array[0] invalid");
  
 delete maskmf;
} // subroutine assign_colors


void NavierStokes::avgDownColor(int idx_color,int idx_type) {

 int finest_level=parent->finestLevel();

 if (level>=finest_level)
  amrex::Error("level invalid avgDownColor");

 int f_level=level+1;
 NavierStokes&   fine_lev = getLevel(f_level);
 const BoxArray& fgrids=fine_lev.grids;
 const DistributionMapping& fdmap=fine_lev.dmap;

 MultiFab* S_crse=localMF[idx_color];
 MultiFab* S_fine=fine_lev.localMF[idx_color];

 if (grids!=S_crse->boxArray())
  amrex::Error("S_crse invalid avgDownColor");
 if (fgrids!=S_fine->boxArray())
  amrex::Error("S_fine invalid");

 BoxArray crse_S_fine_BA(fgrids.size());
 for (int i = 0; i < fgrids.size(); ++i) {
  crse_S_fine_BA.set(i,amrex::coarsen(fgrids[i],2));
 }
 DistributionMapping crse_dmap=fdmap;
 MultiFab crse_S_fine(crse_S_fine_BA,crse_dmap,1,0,
   MFInfo().SetTag("crse_S_fine"),FArrayBoxFactory());
 MultiFab type_coarse_fine(crse_S_fine_BA,crse_dmap,1,0,
   MFInfo().SetTag("type_coarse_fine"),FArrayBoxFactory());
 type_coarse_fine.ParallelCopy(*localMF[idx_type],0,0,1);

 ParallelDescriptor::Barrier();

 int bfact=parent->Space_blockingFactor(level);
 int bfact_f=parent->Space_blockingFactor(level+1);
 const Real* dx = geom.CellSize();
 const Real* dxf = fine_lev.geom.CellSize();
 const Real* prob_lo   = geom.ProbLo();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(S_fine->boxArray().d_numPts());

// COLORING LOOP
#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*S_fine,false); mfi.isValid(); ++mfi) {
  BL_ASSERT(fgrids[mfi.index()] == mfi.validbox());
  const Box& tilegrid = mfi.tilebox();

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  const int gridno = mfi.index();
  const Real* xlo_fine = fine_lev.grid_loc[gridno].lo();

  const Box& fabgrid = crse_S_fine_BA[gridno];
  const int* ovlo=fabgrid.loVect();
  const int* ovhi=fabgrid.hiVect();

  FArrayBox& finefab=(*S_fine)[gridno];
  const Box& fgrid=finefab.box();
  const int* flo=fgrid.loVect();
  const int* fhi=fgrid.hiVect();
  const Real* f_dat=finefab.dataPtr();

  FArrayBox& coarsefab=crse_S_fine[gridno];
  const Box& cgrid = coarsefab.box();
  const int* clo=cgrid.loVect();
  const int* chi=cgrid.hiVect();
  const Real* c_dat=coarsefab.dataPtr();

  const FArrayBox& typef=(*fine_lev.localMF[idx_type])[gridno];
  const FArrayBox& typec=type_coarse_fine[gridno];

  fort_avgdowncolor(
    prob_lo,dxf, 
    &bfact_f,&bfact,
    xlo_fine,dx,
    c_dat,ARLIM(clo),ARLIM(chi),
    f_dat,ARLIM(flo),ARLIM(fhi),
    typef.dataPtr(),ARLIM(typef.loVect()),ARLIM(typef.hiVect()),
    typec.dataPtr(),ARLIM(typec.loVect()),ARLIM(typec.hiVect()),
    ovlo,ovhi);
 } // mfi
} //omp
 ns_reconcile_d_num(LOOP_AVGDOWNCOLOR,"avgDownColor");

 S_crse->ParallelCopy(crse_S_fine,0,0,1);
}


// maskfinemf corresponds to level+1
// components are: color1,type1,color2,type2,color3,type3
MultiFab* NavierStokes::CopyFineToCoarseColor(
  int idx_color,int idx_type,int zero_diag_flag) {

 int finest_level=parent->finestLevel();

 if (level>=finest_level)
  amrex::Error("level invalid CopyFineToCoarseColor");

 int f_level=level+1;
 NavierStokes&   fine_lev = getLevel(f_level);

 MultiFab* S_crse=localMF[idx_color];
 MultiFab* S_fine=fine_lev.localMF[idx_color];

 MultiFab* mf=new MultiFab(grids,dmap,6,1,
  MFInfo().SetTag("mf CopyFineToCoarseColor"),FArrayBoxFactory());
 for (int i=0;i<3;i++) {
  MultiFab::Copy(*mf,*S_crse,0,2*i,1,1);
  MultiFab::Copy(*mf,*localMF[idx_type],0,2*i+1,1,1);
 }

 const BoxArray& fgrids=fine_lev.grids;
 const DistributionMapping& fdmap=fine_lev.dmap;

  // mask=tag if not covered by level+1 and at fine/fine ghost cell.
 int ngrowmask=1;
 Real tag=1.0;
 int clear_phys_boundary=2;
 MultiFab* maskfinemf=fine_lev.maskfiner(ngrowmask,tag,clear_phys_boundary);

 if (grids!=S_crse->boxArray())
  amrex::Error("S_crse invalid CopyFineToCoarseColor");
 if (fgrids!=S_fine->boxArray())
  amrex::Error("S_fine invalid");
 if (fgrids!=maskfinemf->boxArray())
  amrex::Error("maskfinemf invalid");

 BoxArray crse_S_fine_BA(fgrids.size());
 for (int i = 0; i < fgrids.size(); ++i) {
  crse_S_fine_BA.set(i,amrex::coarsen(fgrids[i],2));
 }
 DistributionMapping crse_dmap=fdmap;
 MultiFab crse_S_fine(crse_S_fine_BA,crse_dmap,6,0,
   MFInfo().SetTag("crse_S_fine"),FArrayBoxFactory());

 ParallelDescriptor::Barrier();
 int bfact=parent->Space_blockingFactor(level);
 int bfact_f=parent->Space_blockingFactor(level+1);
 const Real* dx = geom.CellSize();
 const Real* dxf = fine_lev.geom.CellSize();
 const Real* prob_lo   = geom.ProbLo();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(S_fine->boxArray().d_numPts());

// COLORING LOOP
#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*S_fine,false); mfi.isValid(); ++mfi) {
  BL_ASSERT(fgrids[mfi.index()] == mfi.validbox());
  const Box& tilegrid = mfi.tilebox();

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  const int gridno = mfi.index();
  const Real* xlo_fine = fine_lev.grid_loc[gridno].lo();

  const Box& fgrd=fgrids[gridno];
  const int* lofine=fgrd.loVect();
  const int* hifine=fgrd.hiVect();

  const Box& fabgrid = crse_S_fine_BA[gridno];
  const int* ovlo=fabgrid.loVect();
  const int* ovhi=fabgrid.hiVect();

  FArrayBox& finefab=(*S_fine)[mfi];
  const Box& fgridbox=finefab.box();
  const int* flo=fgridbox.loVect();
  const int* fhi=fgridbox.hiVect();
  const Real* f_dat=finefab.dataPtr();

  FArrayBox& coarsefab=crse_S_fine[mfi];
  const Box& cgridbox = coarsefab.box();
  const int* clo=cgridbox.loVect();
  const int* chi=cgridbox.hiVect();
  const Real* c_dat=coarsefab.dataPtr();

  const FArrayBox& typef=(*fine_lev.localMF[idx_type])[mfi];

  FArrayBox & maskfinefab=(*maskfinemf)[mfi];
  
  fort_copyfinecoarsecolor(
    prob_lo,dxf,&bfact_f,&bfact,
    xlo_fine,dx,
    c_dat,ARLIM(clo),ARLIM(chi),
    f_dat,ARLIM(flo),ARLIM(fhi),
    typef.dataPtr(),ARLIM(typef.loVect()),ARLIM(typef.hiVect()),
    maskfinefab.dataPtr(),
    ARLIM(maskfinefab.loVect()),ARLIM(maskfinefab.hiVect()),
    ovlo,ovhi,lofine,hifine,
    &zero_diag_flag);
 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_COPYFINECOARSECOLOR,"CopyFineToCoarseColor");

 mf->ParallelCopy(crse_S_fine,0,0,6);

 delete maskfinemf;

 return mf;
} //end subroutine CopyFineToCoarseColor

// c++ example of iterating components of a FAB
// for more examples, check FArrayBox.cpp:
// for (IntVect p = bx.smallEnd(); p <= bx.bigEnd(); bx.next(p)) ...
void NavierStokes::sync_colors(
 int idx_color,int idx_type,
 Vector<int> color_per_grid,
 Vector<int>& colormax,
 int max_colors_grid,
 MultiFab* maskmf,
 int check_corners,
 int zero_diag_flag) {

 int finest_level=parent->finestLevel();

 if ((check_corners!=0)&&(check_corners!=1))
  amrex::Error("check_corners invalid");

 MultiFab* typemf=localMF[idx_type];
 MultiFab* colormf=localMF[idx_color];
 int number_grids=grids.size();

 int total_colors=0;
 colormf->FillBoundary(geom.periodicity());

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "in sync colors level= " << level << '\n';
   std::cout << "ngrids,maxcolorgrid " << number_grids << ' ' <<
    max_colors_grid << '\n';
  }
 }

 int Nside=number_grids*max_colors_grid;

 if (Nside<=0) {
  std::cout << "cannot have the coarse level completely covered by a finer\n";
  std::cout << "level for the sync_colors routine\n";
  std::cout << "set the coarse level resolution to be level 1 and reduce\n";
  std::cout << "max_level by 1\n";
  amrex::Error("cannot have Nside 0");
 }

 Vector< Vector< Vector<int> > > grid_color_array;
 grid_color_array.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  grid_color_array[tid].resize(Nside);
  for (int i=0;i<Nside;i++) {
   grid_color_array[tid][i].resize(0);
  }
  for (int igrid=0;igrid<number_grids;igrid++) {
   for (int icolor=1;icolor<=color_per_grid[igrid];icolor++) {
    int i=max_colors_grid*igrid+icolor-1;
    grid_color_array[tid][i].resize(1);
    grid_color_array[tid][i][0]=i;
   } // icolor
  } // igrid
 } //tid

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(typemf->boxArray().d_numPts());

// grid_color_array[tid][i][j]=1 if color i neighbors color j
// and they are both not covered by a finer level.
//
// COLORING LOOP
#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*typemf,false); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  const Box& tilegrid = mfi.tilebox();

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  const int gridno = mfi.index();
  const Box& fabgrid=grids[gridno];
  const int* lo = fabgrid.loVect();
  const int* hi = fabgrid.hiVect();
  FArrayBox& maskfab=(*maskmf)[mfi];
  FArrayBox& typefab=(*typemf)[mfi];
  FArrayBox& colorfab=(*colormf)[mfi];

#if (AMREX_SPACEDIM==3)
  for (int i=lo[0];i<=hi[0];i++)
  for (int j=lo[1];j<=hi[1];j++)
  for (int k=lo[2];k<=hi[2];k++) {
#endif
#if (AMREX_SPACEDIM==2)
  for (int i=lo[0];i<=hi[0];i++)
  for (int j=lo[1];j<=hi[1];j++) {
#endif
   IntVect p(D_DECL(i,j,k));

   if (maskfab(p)==1.0) {
    int icolor_round=(int) (colorfab(p)+0.5);
    int icolor_trunc=(int) (colorfab(p));
    int primary_type_round=(int) (typefab(p)+0.5);
    int primary_type_trunc=(int) (typefab(p));

    if ((icolor_round==icolor_trunc)&&
        (primary_type_round==primary_type_trunc)) {
     //do nothing
    } else
     amrex::Error("expecting round=trunc");

    if (icolor_round>0) {
#if (AMREX_SPACEDIM==3)
     for (int ii=-1;ii<=1;ii++)
     for (int jj=-1;jj<=1;jj++)
     for (int kk=-1;kk<=1;kk++) {
      int idist=std::abs(ii)+std::abs(jj)+std::abs(kk);
#endif
#if (AMREX_SPACEDIM==2)
     for (int ii=-1;ii<=1;ii++)
     for (int jj=-1;jj<=1;jj++) {
      int idist=std::abs(ii)+std::abs(jj);
#endif
      IntVect pofs(D_DECL(i+ii,j+jj,k+kk));

      if ((check_corners==1)||(idist<=1)) {

       int jcolor_round=(int) (colorfab(pofs)+0.5);
       int jcolor_trunc=(int) (colorfab(pofs));
       int mask2=(int) maskfab(pofs); 
       int secondary_type_round=(int) (typefab(pofs)+0.5);
       int secondary_type_trunc=(int) (typefab(pofs));
       if ((jcolor_round==jcolor_trunc)&&
           (secondary_type_round==secondary_type_trunc)) {
        //do nothing
       } else
        amrex::Error("expecting round=trunc");

       if (mask2==1) {

        if (jcolor_round<=0) {
         std::cout << "level= " << level << '\n';
         std::cout << "ii,jj= " << ii << ' ' << jj << '\n';
         std::cout << "i,j= " << i << ' ' << j << '\n';
         std::cout << "ngrids,maxcolorgrid " << number_grids << ' ' <<
           max_colors_grid << '\n';
         for (int dir2=0;dir2<AMREX_SPACEDIM;dir2++) {
          std::cout << "dir,lo,hi " << dir2 << ' ' <<
           lo[dir2] << ' ' << hi[dir2] << '\n';
         }
         std::cout << "jcolor_round = " << jcolor_round << '\n';
         amrex::Error("jcolor_round invalid"); 
        } // jcolor<=0
	   
        if (secondary_type_round==primary_type_round) {

         if ((icolor_round>Nside)||(jcolor_round>Nside))
          amrex::Error("icolor or jcolor invalid"); 

	 if (icolor_round!=jcolor_round) {

	  int size_i=grid_color_array[tid_current][icolor_round-1].size();
	  int size_j=grid_color_array[tid_current][jcolor_round-1].size();

	  if ((size_i>=1)&&(size_j>=1)) {

  	   int dup_flag=0;
	   for (int idup=0;idup<size_i;idup++) {
 	    if (grid_color_array[tid_current][icolor_round-1][idup]==
			    jcolor_round-1)
	     dup_flag=1;
	   }
	   for (int idup=0;idup<size_j;idup++) {
 	    if (grid_color_array[tid_current][jcolor_round-1][idup]==
			    icolor_round-1)
	     dup_flag=1;
	   }
           if (dup_flag==0) {
            grid_color_array[tid_current][icolor_round-1].resize(size_i+1);
            grid_color_array[tid_current][icolor_round-1][size_i]=
		    jcolor_round-1;
            grid_color_array[tid_current][jcolor_round-1].resize(size_j+1);
            grid_color_array[tid_current][jcolor_round-1][size_j]=
		    icolor_round-1;
	   }
	  } else
	   amrex::Error("size_i or size_j invalid");

	 } else if (icolor_round==jcolor_round) {
	  //do nothing
	 } else
  	  amrex::Error("icolor or jcolor bust");
				 
        } // primary_type==secondary_type
       } else if (mask2==0) {
        //do nothing
       } else
        amrex::Error("mask2 invalid");

      } else if ((check_corners==0)&&(idist>1)) {
       // do nothing
      } else
       amrex::Error("check_corners or idist invalid");

     } // ii,jj,kk
    } else
     amrex::Error("icolor invalid");
   } else if (maskfab(p)==0.0) {
    //do nothing
   } else
    amrex::Error("maskfab(p) invalid");
  } // i,j,k
 } // mfi
} //omp
 ns_reconcile_d_num(LOOP_SYNC_COLORS,"sync_colors");

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "after initializing grid_color_array \n";
  }
 }

  // first reduce grid_color_array from the different threads
 for (int igrid=0;igrid<number_grids;igrid++) {
  for (int icolor=1;icolor<=color_per_grid[igrid];icolor++) {
   int i_index=max_colors_grid*igrid+icolor-1;

   for (int tid=0;tid<thread_class::nthreads;tid++) {
    int i_size=grid_color_array[tid][i_index].size();
    for (int i_trial=0;i_trial<i_size;i_trial++) {
     int j_index=grid_color_array[tid][i_index][i_trial];
     int dup_flag=0;
     int base_size=grid_color_array[0][i_index].size();
     for (int idup=0;idup<base_size;idup++) {
      if (grid_color_array[0][i_index][idup]==j_index)
       dup_flag=1;
     }
     if (dup_flag==0) {
      grid_color_array[0][i_index].resize(base_size+1);
      grid_color_array[0][i_index][base_size]=j_index;
     }
    } //i_trial
   } //tid
  } // icolor
 } // igrid

 ParallelDescriptor::Barrier();

 for (int igrid=0;igrid<number_grids;igrid++) {
  for (int icolor=1;icolor<=color_per_grid[igrid];icolor++) {
   int i_index=max_colors_grid*igrid+icolor-1;

   ParallelDescriptor::Barrier();

   Vector<int> local_array;
   local_array.resize(Nside);
   for (int j_index=0;j_index<Nside;j_index++)
    local_array[j_index]=0;
   int i_size_sync=grid_color_array[0][i_index].size();
   for (int i_trial=0;i_trial<i_size_sync;i_trial++) {
    int j_index=grid_color_array[0][i_index][i_trial];
    local_array[j_index]=1;
   }

   grid_color_array[0][i_index].resize(0);

   ParallelDescriptor::Barrier();

   for (int jgrid=0;jgrid<number_grids;jgrid++) {
    for (int jcolor=1;jcolor<=color_per_grid[jgrid];jcolor++) {

     int j_index=max_colors_grid*jgrid+jcolor-1;

     ParallelDescriptor::ReduceIntMax(local_array[j_index]);
     if (local_array[j_index]==0) {
      //do nothing
     } else if (local_array[j_index]==1) {
      int i_size_grid=grid_color_array[0][i_index].size();
      grid_color_array[0][i_index].resize(i_size_grid+1);
      grid_color_array[0][i_index][i_size_grid]=j_index;
     } else
      amrex::Error("local_array[j_index] invalid");
    }  // jcolor
   } // jgrid
     
  } // icolor
 } // igrid
 ParallelDescriptor::Barrier();

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "after ReduceIntMax \n";
  }
 }

   // levelcolormap[i] associates color "i" in the previous scheme with a
   // new absolute coloring scheme on the level.
 Vector<int> levelcolormap;
 levelcolormap.resize(Nside);
 for (int i=0;i<Nside;i++) 
  levelcolormap[i]=0;

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "before cross_check \n";
  }
 }

 Vector<int> stackdata;
 int n_assoc=0;


  //levelcolormap only needs to be updated on the IOProcessor and
  //initialized to 0 on the other processors.  A "reduceintmax"
  //command synchronizes the IOProcessor to the others.
 if (ParallelDescriptor::IOProcessor()) {

  n_assoc=0;
  for (int igrid=0;igrid<number_grids;igrid++) {
   for (int icolor=1;icolor<=color_per_grid[igrid];icolor++) {
    int i_index=max_colors_grid*igrid+icolor-1;
    n_assoc+=grid_color_array[0][i_index].size();
   }
  }

  stackdata.resize(n_assoc);

  for (int igrid=0;igrid<number_grids;igrid++) {
   for (int icolor=1;icolor<=color_per_grid[igrid];icolor++) {
    int i_index=max_colors_grid*igrid+icolor-1;
    if (levelcolormap[i_index]==0) {
     if (grid_color_array[0][i_index].size()>0) {
      total_colors++;
      levelcolormap[i_index]=total_colors;
      if (levelcolormap[i_index]<1)
       amrex::Error("levelcolormap invalid");
      cross_check(levelcolormap,stackdata,grid_color_array[0],i_index);
     } else
      amrex::Error("expecting grid_color_array[0][i_index].size()>0");
    } 
   } // icolor
  } // igrid

  stackdata.resize(0);
 }  // IOProcessor

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "after first cross_check: level,ngrids,n_assoc " <<
	   level << ' ' << number_grids << ' ' << n_assoc << '\n';
  }
 }

 for (int tid=0;tid<thread_class::nthreads;tid++) {
  grid_color_array[tid].resize(1);
 }

 ParallelDescriptor::Barrier();

 for (int igrid=0;igrid<number_grids;igrid++) {
  for (int icolor=1;icolor<=color_per_grid[igrid];icolor++) {
   int i_index=max_colors_grid*igrid+icolor-1;
   ParallelDescriptor::ReduceIntMax(levelcolormap[i_index]);
  }
 }
 ParallelDescriptor::ReduceIntMax(total_colors);
 ParallelDescriptor::Barrier();

 const Real* dx = geom.CellSize();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(typemf->boxArray().d_numPts());

 // resets one ghost layer.
 // switches from local coloring scheme on the level to a global
 // coloring scheme on the level.

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*typemf,false); mfi.isValid(); ++mfi) {
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

  FArrayBox& maskfab=(*maskmf)[mfi];
  FArrayBox& colorfab=(*colormf)[mfi];
  int arrsize=levelcolormap.size();

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // in: LEVELSET_3D.F90
  fort_gridrecolor(
    maskfab.dataPtr(),ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
    colorfab.dataPtr(),ARLIM(colorfab.loVect()),ARLIM(colorfab.hiVect()),
    xlo,dx,
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    levelcolormap.dataPtr(),
    &max_colors_grid,
    &number_grids,&arrsize);
 }  //mfi
} //omp
 ns_reconcile_d_num(LOOP_GRIDRECOLOR,"sync_colors");

 colormax[level]=total_colors;

 if (level<finest_level) {

  avgDownColor(idx_color,idx_type);
   
   // 1 ghost cell is initialized
   // components are: color1,type1,color2,type2,color3,type3
   // only copies where maskfine=1
  MultiFab* fine_coarse_color=
    CopyFineToCoarseColor(idx_color,idx_type,zero_diag_flag);
  fine_coarse_color->FillBoundary(geom.periodicity());

  int max_colors_level=colormax[level+1];
  if (colormax[level]>max_colors_level)
   max_colors_level=colormax[level];
  int arrsize=2*max_colors_level;

  Vector< Vector< Vector<int> > > level_color_array;
  level_color_array.resize(thread_class::nthreads);

  for (int tid=0;tid<thread_class::nthreads;tid++) {
   level_color_array[tid].resize(arrsize);
   for (int i=0;i<arrsize;i++) {
    level_color_array[tid][i].resize(0);
   }
   for (int ilevel=0;ilevel<=1;ilevel++) {
    for (int icolor=1;icolor<=colormax[level+ilevel];icolor++) {
     int i=max_colors_level*ilevel+icolor-1;
     level_color_array[tid][i].resize(1);
     level_color_array[tid][i][0]=i;
    } // icolor
   } // ilevel
  } //tid

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(typemf->boxArray().d_numPts());

   // level_color_array[i][j]=1 if color i neighbors color j and
   // the mask of one of them is 1.
   // COLORING LOOP
#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(*typemf,false); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* lo = fabgrid.loVect();
   const int* hi = fabgrid.hiVect();

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   FArrayBox& maskfab=(*maskmf)[mfi];

   //components are: color1,type1,color2,type2,color3,type3
   FArrayBox& colorfab=(*fine_coarse_color)[mfi];

#if (AMREX_SPACEDIM==3)
   for (int i=lo[0];i<=hi[0];i++)
   for (int j=lo[1];j<=hi[1];j++)
   for (int k=lo[2];k<=hi[2];k++) {
#endif
#if (AMREX_SPACEDIM==2)
   for (int i=lo[0];i<=hi[0];i++)
   for (int j=lo[1];j<=hi[1];j++) {
#endif
    for (int nbase=1;nbase<=3;nbase++) {

     IntVect p(D_DECL(i,j,k));

     int icolor_round=(int) (colorfab(p,2*nbase-2)+0.5);
     int icolor_trunc=(int) (colorfab(p,2*nbase-2));
     int base_type_round=(int) (colorfab(p,2*nbase-1)+0.5);
     int base_type_trunc=(int) (colorfab(p,2*nbase-1));
     if ((icolor_round==icolor_trunc)&&
         (base_type_round==base_type_trunc)) {
      //do nothing
     } else
      amrex::Error("expecting round=trunc");

     int local_mask=(int) maskfab(p);

     if ((icolor_round>0)||(local_mask==0)) {
      //do nothing
     } else {
      amrex::Error("icolor_round>0 || local_mask==0 failed");
     }

     if (icolor_round>0) {
      if (local_mask==0) {
       icolor_round+=max_colors_level;
      } else if (local_mask==1) {
       //do nothing
      } else
       amrex::Error("local_mask invalid");

#if (AMREX_SPACEDIM==3)
      for (int ii=-1;ii<=1;ii++)
      for (int jj=-1;jj<=1;jj++)
      for (int kk=-1;kk<=1;kk++) {
       int idist=std::abs(ii)+std::abs(jj)+std::abs(kk);
#endif
#if (AMREX_SPACEDIM==2)
      for (int ii=-1;ii<=1;ii++)
      for (int jj=-1;jj<=1;jj++) {
       int idist=std::abs(ii)+std::abs(jj);
#endif
       IntVect pofs(D_DECL(i+ii,j+jj,k+kk));

       if ((check_corners==1)||(idist<=1)) {
        int mask2=(int) maskfab(pofs); 
        for (int nbase2=1;nbase2<=3;nbase2++) {
         int jcolor_round=(int) (colorfab(pofs,2*nbase2-2)+0.5);
         int jcolor_trunc=(int) (colorfab(pofs,2*nbase2-2));
         int near_type_round=(int) (colorfab(pofs,2*nbase2-1)+0.5);
         int near_type_trunc=(int) (colorfab(pofs,2*nbase2-1));
         if ((jcolor_round==jcolor_trunc)&&
             (near_type_round==near_type_trunc)) {
          //do nothing
         } else
          amrex::Error("expecting round=trunc");

         if ((jcolor_round>0)||(mask2==0)) {
          //do nothing
         } else {
          amrex::Error("jcolor_round>0 || mask2==0 failed");
         }
         if (jcolor_round>0) {
          if (near_type_round==base_type_round) {
           if (mask2==0) {
            jcolor_round+=max_colors_level;
           } else if (mask2==1) {
            //do nothing
           } else
            amrex::Error("mask2 invalid");

           if ((icolor_round>arrsize)||(jcolor_round>arrsize))
            amrex::Error("icolor_round>arrsize or jcolor_round>arrsize");

	   if ((mask2==1)||(local_mask==1)) {

            if (icolor_round!=jcolor_round) {
 	     int size_i=level_color_array[tid_current][icolor_round-1].size();
   	     int size_j=level_color_array[tid_current][jcolor_round-1].size();
  	     if ((size_i>=1)&&(size_j>=1)) {

 	      int dup_flag=0;
  	      for (int idup=0;idup<size_i;idup++) {
 	       if (level_color_array[tid_current][icolor_round-1][idup]==
			       jcolor_round-1)
	        dup_flag=1;
	      }
	      for (int idup=0;idup<size_j;idup++) {
 	       if (level_color_array[tid_current][jcolor_round-1][idup]==
			       icolor_round-1)
	        dup_flag=1;
	      }

              if (dup_flag==0) {
               level_color_array[tid_current][icolor_round-1].resize(size_i+1);
               level_color_array[tid_current][icolor_round-1][size_i]=
		       jcolor_round-1;
               level_color_array[tid_current][jcolor_round-1].resize(size_j+1);
               level_color_array[tid_current][jcolor_round-1][size_j]=
		       icolor_round-1;
	      }

	     } else
	      amrex::Error("size_i or size_j invalid");
	    } else if (icolor_round==jcolor_round) {
	     //do nothing
	    } else
	     amrex::Error("icolor_round or jcolor_round bust");

	   } else if ((mask2==0)&&(local_mask==0)) {
	    //do nothing
	   } else
  	    amrex::Error("mask2 or local_mask invalid");
	  } else if (near_type_round!=base_type_round) {
	   //do nothing
	  } else
  	   amrex::Error("near_type or base_type invalid");
	 } else if (jcolor_round==0) {
          //do nothing
	 } else
          amrex::Error("jcolor invalid");
        } //for (int nbase2=1;nbase2<=3;nbase2++) 
       } else if ((check_corners==0)&&(idist>1)) {
        //do nothing
       } else
        amrex::Error("check_corners or idist invalid");
      } //ii,jj,kk
     } else if (icolor_round==0) {
      //do nothing
     } else
      amrex::Error("icolor_round invalid");
    } //for (int nbase=1;nbase<=3;nbase++) 
     
   } //i,j,k

  } //mfi
} //omp
  ns_reconcile_d_num(LOOP_LEVELCOLORINIT,"sync_colors");

  delete fine_coarse_color;

  // first reduce level_color_array from the different threads
  for (int ilevel=0;ilevel<=1;ilevel++) {
   for (int icolor=1;icolor<=colormax[level+ilevel];icolor++) {
    int i_index=max_colors_level*ilevel+icolor-1;

    for (int tid=0;tid<thread_class::nthreads;tid++) {
     int i_size_level2=level_color_array[tid][i_index].size();
     for (int i_trial=0;i_trial<i_size_level2;i_trial++) {
      int j_index=level_color_array[tid][i_index][i_trial];
      int dup_flag=0;
      int base_size=level_color_array[0][i_index].size();
      for (int idup=0;idup<base_size;idup++) {
       if (level_color_array[0][i_index][idup]==j_index)
        dup_flag=1;
      }
      if (dup_flag==0) {
       level_color_array[0][i_index].resize(base_size+1);
       level_color_array[0][i_index][base_size]=j_index;
      }
     } //i_trial
    } //tid
   } // icolor
  } // ilevel

  ParallelDescriptor::Barrier();

  for (int ilevel=0;ilevel<=1;ilevel++) {
   for (int icolor=1;icolor<=colormax[level+ilevel];icolor++) {
    int i_index=max_colors_level*ilevel+icolor-1;

    ParallelDescriptor::Barrier();

    Vector<int> local_array;
    local_array.resize(arrsize);
    for (int j_index=0;j_index<arrsize;j_index++)
     local_array[j_index]=0;

    int i_size_level_sync=level_color_array[0][i_index].size();
    for (int i_trial=0;i_trial<i_size_level_sync;i_trial++) {
     int j_index=level_color_array[0][i_index][i_trial];
     local_array[j_index]=1;
    }

    level_color_array[0][i_index].resize(0);

    ParallelDescriptor::Barrier();

    for (int jlevel=0;jlevel<=1;jlevel++) {
     for (int jcolor=1;jcolor<=colormax[level+jlevel];jcolor++) {
      int j_index=max_colors_level*jlevel+jcolor-1;
      ParallelDescriptor::ReduceIntMax(local_array[j_index]);
      if (local_array[j_index]==0) {
       //do nothing
      } else if (local_array[j_index]==1) {
       int i_size_level3=level_color_array[0][i_index].size();
       level_color_array[0][i_index].resize(i_size_level3+1);
       level_color_array[0][i_index][i_size_level3]=j_index;
      } else
       amrex::Error("local_array[j_index] invalid");
     }  // jcolor
    } // jlevel
     
   } // icolor
  } // ilevel
  ParallelDescriptor::Barrier();

  Vector<int> domaincolormap;
  domaincolormap.resize(arrsize);
  for (int i=0;i<arrsize;i++)
   domaincolormap[i]=0;

  total_colors=0;

  stackdata.resize(0);
  n_assoc=0;
  for (int ilevel=0;ilevel<=1;ilevel++) {
   for (int icolor=1;icolor<=colormax[level+ilevel];icolor++) {
    int i_index=max_colors_level*ilevel+icolor-1;
    n_assoc+=level_color_array[0][i_index].size();
   }
  }
  stackdata.resize(n_assoc);

  for (int ilevel=0;ilevel<=1;ilevel++) {
   for (int icolor=1;icolor<=colormax[level+ilevel];icolor++) {
    int i_index=max_colors_level*ilevel+icolor-1;
    if (domaincolormap[i_index]==0) {
     if (level_color_array[0][i_index].size()>0) {
      total_colors++;
      domaincolormap[i_index]=total_colors;
       cross_check(domaincolormap,stackdata,level_color_array[0],i_index);
     } else
      amrex::Error("expecting level_color_array[0][i_index].size()>0");
    }
   } //icolor
  } //ilevel
  
  stackdata.resize(0);

  for (int tid=0;tid<thread_class::nthreads;tid++) {
   level_color_array[tid].resize(1);
  }

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "after 2nd cross_check: level,max_colors_level,n_assoc " <<
     level << ' ' << max_colors_level << ' ' << n_assoc << '\n';
  }
 }



  for (int ilev=finest_level;ilev>=level;ilev--) {
   colormax[ilev]=total_colors;
   NavierStokes& ns_level=getLevel(ilev);
    
   ns_level.correct_colors(idx_color,level,domaincolormap,max_colors_level);
   // if coarse_type=fine_type, then coarse_color=fine_color otherwise
   // coarse_color=0.
   if (ilev<finest_level)
    ns_level.avgDownColor(idx_color,idx_type);
  }
 } // level<finest_level


} // subroutine sync_colors
 
// type_flag should have size num_materials.
// type_flag[i]=1 if fluid "i" exists. (note, fictitious solid on the
//  boundaries will show up as existing if ngrow>0)
void NavierStokes::color_variable(
 int& coarsest_level,
 int idx_color,int idx_type,
 int* color_count,
 Vector<int> type_flag,
 int zero_diag_flag) {

int ilev;

 int finest_level=parent->finestLevel();

 if (level!=0)
  amrex::Error("level invalid color_variable");

 allocate_array(1,1,-1,idx_color);

 Vector<int> colormax;
 colormax.resize(finest_level+1);

 int fully_covered=0;
 for (ilev=finest_level;((ilev>=0)&&(fully_covered==0));ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  
  ns_level.assign_colors(fully_covered,
   idx_color,idx_type,colormax,type_flag,
   zero_diag_flag); 
  if (fully_covered==0) {
   // do nothing
  } else if (fully_covered==1) {
   coarsest_level=ilev+1;
  } else
   amrex::Error("fully_covered invalid");
 }
 if (fully_covered==0) {
  coarsest_level=0;
 } else if (fully_covered==1) {
  if ((coarsest_level<=0)||(coarsest_level>finest_level))
   amrex::Error("coarsest_level invalid");
 } else
  amrex::Error("fully_covered invalid");

 *color_count=colormax[coarsest_level];

  // fort_extrapfill, pc_interp
 Vector<int> scompBC_map;
 scompBC_map.resize(1);
 scompBC_map[0]=0;

  // ngrow=1 scomp=0 ncomp=1 
  // fort_extrapfill, pc_interp
 PCINTERP_fill_bordersALL(idx_color,1,0,1,State_Type,scompBC_map);

 if (verbose>0)
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "coarsest_level= " << coarsest_level << '\n';
   for (ilev=coarsest_level;ilev<=finest_level;ilev++)
    std::cout << "ilev,colormax_ilev " << ilev << ' ' << 
     colormax[ilev] << '\n';
  }
} // subroutine color_variable

//called from NavierStokes::ColorSumALL
//operation_flag==1 OP_SCATTER_MDOT => scatter data collected when 
//operation_flag==0 (OP_GATHER_MDOT) to mdot or density.
void
NavierStokes::ColorSum(
 int operation_flag, //OP_GATHER_MDOT or OP_SCATTER_MDOT
 int tessellate,  // =1 or 3
 int sweep_num,
 int ncomp_mdot_alloc,
 int ncomp_mdot,
 MultiFab* typemf,
 MultiFab* color,
 MultiFab* mdot, // holds typemf if ncomp_mdot==0
 MultiFab* mdot_complement, // holds typemf if ncomp_mdot==0
 Vector<blobclass>& level_blobdata,
 Vector<blobclass> cum_blobdata,
 Vector< Vector<Real> >& level_mdot_data,
 Vector< Vector<Real> >& level_mdot_comp_data,
 Vector< Vector<Real> > cum_mdot_data,
 Vector< Vector<Real> > cum_mdot_comp_data,
 Vector< Vector<Real> >& level_mdot_data_redistribute,
 Vector< Vector<Real> >& level_mdot_comp_data_redistribute
 ) {

 std::string local_caller_string="ColorSum";

 int finest_level=parent->finestLevel();
 bool use_tiling=ns_tiling;

 if (level>finest_level)
  amrex::Error("level invalid ColorSum");

 if (ngrow_distance!=4)
  amrex::Error("ngrow_distance invalid");

 int num_colors=level_blobdata.size();
 if (num_colors==0) 
  amrex::Error("num_colors should be positive in ColorSum");
 if (num_colors!=cum_blobdata.size())
  amrex::Error("num_colors!=cum_blobdata.size()");

 if (num_colors!=cum_mdot_data.size())
  amrex::Error("num_colors!=cum_mdot_data.size()");
 if (num_colors!=level_mdot_data.size())
  amrex::Error("num_colors!=level_mdot_data.size()");
 if (num_colors!=level_mdot_data_redistribute.size())
  amrex::Error("num_colors!=level_mdot_data_redistribute.size()");

 if (num_colors!=cum_mdot_comp_data.size())
  amrex::Error("num_colors!=cum_mdot_comp_data.size()");
 if (num_colors!=level_mdot_comp_data.size())
  amrex::Error("num_colors!=level_mdot_comp_data.size()");
 if (num_colors!=level_mdot_comp_data_redistribute.size())
  amrex::Error("num_colors!=level_mdot_comp_data_redistribute.size()");

 if (ncomp_mdot==0) {
  if (ncomp_mdot_alloc==1) {
   // do nothing  
  } else
   amrex::Error("ncomp_mdot_alloc invalid");
 } else if (ncomp_mdot>0) {
  int ncomp_mdot_test=mdot->nComp();
  if (ncomp_mdot_test==2*num_interfaces) {
   // do nothing
  } else
   amrex::Error("ncomp_mdot_test invalid");
  if ((ncomp_mdot==ncomp_mdot_test)&&
      (ncomp_mdot==ncomp_mdot_alloc)) {
   // do nothing
  } else
   amrex::Error("ncomp_mdot_test or ncomp_mdot_alloc invalid");
  ncomp_mdot_test=mdot_complement->nComp();
  if (ncomp_mdot_test==2*num_interfaces) {
   // do nothing
  } else
   amrex::Error("ncomp_mdot_test invalid");
 } else
  amrex::Error("ncomp_mdot invalid"); 

 if (operation_flag==OP_GATHER_MDOT) {

  if ((sweep_num==0)||(sweep_num==1)) {

   for (int i=0;i<num_colors;i++) {
    clear_blobdata(i,level_blobdata);
    for (int j=0;j<ncomp_mdot_alloc;j++) {
     level_mdot_data[i][j]=0.0;
     level_mdot_data_redistribute[i][j]=0.0;
     level_mdot_comp_data[i][j]=0.0;
     level_mdot_comp_data_redistribute[i][j]=0.0;
    }
   } // i=0..num_colors-1
     
  } else
   amrex::Error("sweep_num invalid");

 } else if (operation_flag==OP_SCATTER_MDOT) {

  if (sweep_num==0) {

   if (ncomp_mdot>=1) {
    for (int i=0;i<num_colors;i++) {
     int j=0;
     for (j=0;j<ncomp_mdot_alloc;j++) {
      level_mdot_data_redistribute[i][j]=0.0;
      level_mdot_comp_data_redistribute[i][j]=0.0;
     }
     if (j==2*num_interfaces) {
      // do nothing
     } else 
      amrex::Error("expecting j==2*num_interfaces");
    } // i=0..num_colors-1
   } else
    amrex::Error("ncomp_mdot invalid");
  } else
   amrex::Error("sweep_num invalid");

 } else
  amrex::Error("operation_flag invalid ::ColorSum");

 int blob_array_size=num_colors*num_elements_blobclass;

 int mdot_array_size=num_colors*ncomp_mdot_alloc;

 Vector<Real> cum_blob_array;
 cum_blob_array.resize(blob_array_size);

 Vector<Real> cum_mdot_array;
 cum_mdot_array.resize(mdot_array_size);

 Vector<Real> cum_mdot_comp_array;
 cum_mdot_comp_array.resize(mdot_array_size);

 int counter=0;
 int mdot_counter=0;

 for (int i=0;i<num_colors;i++) {
  copy_from_blobdata(i,counter,cum_blob_array,cum_blobdata);
  for (int j=0;j<ncomp_mdot_alloc;j++) {
   cum_mdot_array[mdot_counter]=cum_mdot_data[i][j];
   cum_mdot_comp_array[mdot_counter]=cum_mdot_comp_data[i][j];
   mdot_counter++;
  }
 }  // i=0..num_colors-1

 if (counter!=blob_array_size)
  amrex::Error("counter invalid");
 if (mdot_counter!=mdot_array_size) {
  std::cout << "mdot_counter=" << mdot_counter << '\n';
  std::cout << "mdot_array_size=" << mdot_array_size << '\n';
  amrex::Error("mdot_counter invalid in ColorSum1");
 }

 Vector< Vector<Real> > level_blob_array;
 Vector< Vector<int> > level_blob_type_array;
 level_blob_array.resize(thread_class::nthreads);
 level_blob_type_array.resize(thread_class::nthreads);

 Vector< Vector<Real> > level_mdot_array;
 level_mdot_array.resize(thread_class::nthreads);

 Vector< Vector<Real> > level_mdot_comp_array;
 level_mdot_comp_array.resize(thread_class::nthreads);

 Vector< Vector<Real> > level_mdot_redistribute_array;
 level_mdot_redistribute_array.resize(thread_class::nthreads);

 Vector< Vector<Real> > level_mdot_comp_redistribute_array;
 level_mdot_comp_redistribute_array.resize(thread_class::nthreads);

 for (int tid=0;tid<thread_class::nthreads;tid++) {
  level_blob_type_array[tid].resize(num_colors);
  level_blob_array[tid].resize(blob_array_size);
  for (int i=0;i<blob_array_size;i++) 
   level_blob_array[tid][i]=0.0;
  for (int i=0;i<num_colors;i++) 
   level_blob_type_array[tid][i]=0;

  level_mdot_array[tid].resize(mdot_array_size);
  level_mdot_redistribute_array[tid].resize(mdot_array_size);

  level_mdot_comp_array[tid].resize(mdot_array_size);
  level_mdot_comp_redistribute_array[tid].resize(mdot_array_size);

  for (int i=0;i<mdot_array_size;i++) {
   level_mdot_array[tid][i]=0.0;
   level_mdot_redistribute_array[tid][i]=0.0;
   level_mdot_comp_array[tid][i]=0.0;
   level_mdot_comp_redistribute_array[tid][i]=0.0;
  }

  if (operation_flag==OP_GATHER_MDOT) {

   // do nothing
   
  } else if (operation_flag==OP_SCATTER_MDOT) {
   
   counter=0;
   mdot_counter=0;
 
   for (int i=0;i<num_colors;i++) {
    copy_from_blobdata(i,counter,level_blob_array[tid],level_blobdata);
    if ((level_blobdata[i].im>=1)&&
        (level_blobdata[i].im<=num_materials)) {
     level_blob_type_array[tid][i]=level_blobdata[i].im;
    } else
     amrex::Error("level_blobdata[i].im invalid");

    int j=0;
    for (j=0;j<ncomp_mdot_alloc;j++) {
     level_mdot_array[tid][mdot_counter]=level_mdot_data[i][j];
     level_mdot_comp_array[tid][mdot_counter]=level_mdot_comp_data[i][j];
     mdot_counter++;
    }
    if (j==2*num_interfaces) {
     // do nothing
    } else
     amrex::Error("j invalid");

   } // i=0..num_colors
   if (counter!=blob_array_size)
    amrex::Error("counter invalid");
   if (mdot_counter!=mdot_array_size) {
    std::cout << "mdot_counter=" << mdot_counter << '\n';
    std::cout << "mdot_array_size=" << mdot_array_size << '\n';
    amrex::Error("mdot_counter invalid in ColorSum2");
   }

  } else
   amrex::Error("operation_flag invalid");

 } // tid=0..thread_class::nthreads-1

 resize_metrics(1);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  debug_ngrow(AREA_MF+dir,1,local_caller_string);
  MultiFab& Umac_new=get_new_data(Umac_Type+dir,slab_step+1);
  if (localMF[AREA_MF+dir]->boxArray()!=Umac_new.boxArray())
   amrex::Error("area_mf boxarrays do not match");
 } // dir=0..sdim-1
 
 debug_ngrow(VOLUME_MF,0,local_caller_string);
 VOF_Recon_resize(1); //output:SLOPE_RECON_MF
 debug_ngrow(SLOPE_RECON_MF,1,local_caller_string);
 if (localMF[SLOPE_RECON_MF]->nComp()!=num_materials*ngeom_recon)
  amrex::Error("localMF[SLOPE_RECON_MF]->nComp() invalid");

 if (sweep_num==0) {
  getStateDist_localMF(LS_COLORSUM_MF,1,cur_time_slab,local_caller_string);
  getStateDen_localMF(DEN_COLORSUM_MF,1,cur_time_slab);
   // velocity + pressure
  getState_localMF(VEL_COLORSUM_MF,1,STATECOMP_VEL,
    STATE_NCOMP_VEL+STATE_NCOMP_PRES,cur_time_slab);

  makeFaceFrac(tessellate,ngrow_distance,FACEFRAC_MM_MF);
  ProcessFaceFrac(tessellate,FACEFRAC_MM_MF,FACEFRAC_SOLVE_MM_MF,0);
  makeCellFrac(tessellate,0,CELLFRAC_MM_MF);
 } else if (sweep_num==1) {
  // do nothing
 } else
  amrex::Error("sweep_num invalid");

 if (localMF[LS_COLORSUM_MF]->nGrow()!=1)
  amrex::Error("localMF[LS_COLORSUM_MF]->nGrow()!=1");

 if (localMF[LS_COLORSUM_MF]->nComp()!=(1+AMREX_SPACEDIM)*num_materials)
  amrex::Error("localMF[LS_COLORSUM_MF]->nComp()!=(1+AMREX_SPACEDIM)*num_materials");
 if (localMF[VEL_COLORSUM_MF]->nComp()!=STATE_NCOMP_VEL+STATE_NCOMP_PRES)
  amrex::Error("localMF[VEL_COLORSUM_MF]->nComp()!=NCOMP_VEL+NCOMP_PRES");
 if (localMF[DEN_COLORSUM_MF]->nComp()!=num_state_material*num_materials)
  amrex::Error("localMF[DEN_COLORSUM_MF]->nComp()!=num_state_material*num_materials");

   // (num_materials,sdim,2) area on each face of a cell.
 int nface=num_materials*AMREX_SPACEDIM*2;
  // (num_materials,num_materials,2)  
  // left material, right material, frac_pair+dist_pair
 int nface_dst=num_materials*num_materials*2;
  // (num_materials,num_materials,3+sdim)
  // im_inside,im_outside,3+sdim --> area, dist_to_line, dist, line normal.
 int ncellfrac=num_materials*num_materials*(3+AMREX_SPACEDIM);

 debug_ngrow(FACEFRAC_MM_MF,ngrow_distance,local_caller_string);
 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(FACEFRAC_SOLVE_MM_MF+dir,0,local_caller_string);
 debug_ngrow(CELLFRAC_MM_MF,0,local_caller_string);

 if (localMF[FACEFRAC_MM_MF]->nComp()!=nface)
  amrex::Error("localMF[FACEFRAC_MM_MF]->nComp()!=nface");
 if (localMF[FACEFRAC_SOLVE_MM_MF]->nComp()!=nface_dst)
  amrex::Error("localMF[FACEFRAC_SOLVE_MM_MF]->nComp()!=nface_dst");
 if (localMF[CELLFRAC_MM_MF]->nComp()!=ncellfrac)
  amrex::Error("localMF[CELLFRAC_MM_MF]->nComp()!=ncellfrac");

 if (typemf->nGrow()!=1)
  amrex::Error("typemf->nGrow()!=1");
 if (color->nGrow()!=1)
  amrex::Error("color->nGrow()!=1");
 if (mdot->nGrow()>=0) {
  // do nothing
 } else
  amrex::Error("mdot->nGrow() invalid");
 if (mdot_complement->nGrow()>=0) {
  // do nothing
 } else
  amrex::Error("mdot_complement->nGrow() invalid");

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 int nstate=STATE_NCOMP;
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");


  // mask=tag if not covered by level+1 and at fine/fine ghost cell.
 int ngrowmask=1;
 Real tag=1.0;
 int clear_phys_boundary=2;
 MultiFab* mask=maskfiner(ngrowmask,tag,clear_phys_boundary);  
 const Real* dx = geom.CellSize();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(mask->boxArray().d_numPts());

// COLORING LOOP
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

  FArrayBox& snewfab=S_new[mfi];

  FArrayBox& mdotfab=(*mdot)[mfi];
  FArrayBox& mdot_comp_fab=(*mdot_complement)[mfi];

  FArrayBox& typefab=(*typemf)[mfi];
  FArrayBox& lsfab=(*localMF[LS_COLORSUM_MF])[mfi];
  FArrayBox& velfab=(*localMF[VEL_COLORSUM_MF])[mfi];
  FArrayBox& denfab=(*localMF[DEN_COLORSUM_MF])[mfi];
  FArrayBox& voffab=(*localMF[SLOPE_RECON_MF])[mfi];

  FArrayBox& xfacepair=(*localMF[FACEFRAC_SOLVE_MM_MF])[mfi];
  FArrayBox& yfacepair=(*localMF[FACEFRAC_SOLVE_MM_MF+1])[mfi];
  FArrayBox& zfacepair=(*localMF[FACEFRAC_SOLVE_MM_MF+AMREX_SPACEDIM-1])[mfi];
  FArrayBox& cellfab=(*localMF[CELLFRAC_MM_MF])[mfi];

  FArrayBox& colorfab=(*color)[mfi];
  FArrayBox& maskfab=(*mask)[mfi];

  FArrayBox& areax=(*localMF[AREA_MF])[mfi];
  FArrayBox& areay=(*localMF[AREA_MF+1])[mfi];
  FArrayBox& areaz=(*localMF[AREA_MF+AMREX_SPACEDIM-1])[mfi];

  Vector<int> levelbc=getBCArray(LS_Type,gridno,0,1);
  Vector<int> velbc=getBCArray(State_Type,gridno,
     STATECOMP_VEL,STATE_NCOMP_VEL);

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // declared in: LEVELSET_3D.F90
  fort_getcolorsum(
   &tid_current,
   &operation_flag,//OP_GATHER_MDOT || OP_SCATTER_MDOT
   &sweep_num,
   &tessellate,
   distribute_mdot_evenly.dataPtr(),
   constant_volume_mdot.dataPtr(),
   distribute_from_target.dataPtr(),
   constant_density_all_time.dataPtr(),
   &cur_time_slab,
   &dt_slab, //used if "OP_SCATTER_MDOT"
   dx, 
   xlo, 
   &nstate,
   snewfab.dataPtr(),ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()),
   mdotfab.dataPtr(),
   ARLIM(mdotfab.loVect()),
   ARLIM(mdotfab.hiVect()),
   mdot_comp_fab.dataPtr(),
   ARLIM(mdot_comp_fab.loVect()),
   ARLIM(mdot_comp_fab.hiVect()),
   lsfab.dataPtr(),ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
   velfab.dataPtr(),ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
   denfab.dataPtr(),ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
   voffab.dataPtr(),ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),
   xfacepair.dataPtr(),ARLIM(xfacepair.loVect()),ARLIM(xfacepair.hiVect()),
   yfacepair.dataPtr(),ARLIM(yfacepair.loVect()),ARLIM(yfacepair.hiVect()),
   zfacepair.dataPtr(),ARLIM(zfacepair.loVect()),ARLIM(zfacepair.hiVect()),
   areax.dataPtr(),ARLIM(areax.loVect()),ARLIM(areax.hiVect()),
   areay.dataPtr(),ARLIM(areay.loVect()),ARLIM(areay.hiVect()),
   areaz.dataPtr(),ARLIM(areaz.loVect()),ARLIM(areaz.hiVect()),
   cellfab.dataPtr(),ARLIM(cellfab.loVect()),ARLIM(cellfab.hiVect()),
   typefab.dataPtr(),ARLIM(typefab.loVect()),ARLIM(typefab.hiVect()),
   colorfab.dataPtr(),ARLIM(colorfab.loVect()),ARLIM(colorfab.hiVect()),
   maskfab.dataPtr(),ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
   tilelo,tilehi,
   fablo,fabhi,
   &bfact,
   &level,
   &finest_level,
   &NS_geometry_coord,
   &num_colors,
   cum_blob_array.dataPtr(),
   cum_mdot_array.dataPtr(),
   cum_mdot_comp_array.dataPtr(),
   level_blob_array[tid_current].dataPtr(),
   level_blob_type_array[tid_current].dataPtr(),
   level_mdot_array[tid_current].dataPtr(),
   level_mdot_comp_array[tid_current].dataPtr(),
   level_mdot_redistribute_array[tid_current].dataPtr(),
   level_mdot_comp_redistribute_array[tid_current].dataPtr(),
   &blob_array_size,
   &mdot_array_size,
   &ncomp_mdot_alloc,
   &ncomp_mdot,
   levelbc.dataPtr(),
   velbc.dataPtr(),
   material_type_lowmach.dataPtr(),
   material_type_visual.dataPtr(),
   &nface_dst,
   &ncellfrac);
 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_GETCOLORSUM,"ColorSum");

 if (operation_flag==OP_GATHER_MDOT) {

  for (int tid=1;tid<thread_class::nthreads;tid++) {
   for (int i=0;i<blob_array_size;i++) {
    level_blob_array[0][i]+=level_blob_array[tid][i];
   }
   for (int i=0;i<num_colors;i++) {
    if (level_blob_type_array[tid][i]>level_blob_type_array[0][i])
     level_blob_type_array[0][i]=level_blob_type_array[tid][i];
   }
   for (int i=0;i<mdot_array_size;i++) {
    level_mdot_array[0][i]+=level_mdot_array[tid][i];
    level_mdot_comp_array[0][i]+=level_mdot_comp_array[tid][i];
   }
  } // tid

 } else if (operation_flag==OP_SCATTER_MDOT) {

  for (int tid=1;tid<thread_class::nthreads;tid++) {
   for (int i=0;i<mdot_array_size;i++) {
    level_mdot_redistribute_array[0][i]+=
      level_mdot_redistribute_array[tid][i];
    level_mdot_comp_redistribute_array[0][i]+=
      level_mdot_comp_redistribute_array[tid][i];
   }
  }

 } else
  amrex::Error("operation_flag invalid");
 
 ParallelDescriptor::Barrier();

 if (operation_flag==OP_GATHER_MDOT) {

  for (int i=0;i<blob_array_size;i++) 
   ParallelDescriptor::ReduceRealSum(level_blob_array[0][i]);
  for (int i=0;i<num_colors;i++) 
   ParallelDescriptor::ReduceIntMax(level_blob_type_array[0][i]);

  for (int i=0;i<mdot_array_size;i++) {
   ParallelDescriptor::ReduceRealSum(level_mdot_array[0][i]);
   ParallelDescriptor::ReduceRealSum(level_mdot_comp_array[0][i]);
  }

  counter=0;
  mdot_counter=0;

  for (int i=0;i<num_colors;i++) {
   copy_to_blobdata(i,counter,level_blob_array[0],level_blobdata);
   level_blobdata[i].im=level_blob_type_array[0][i];

   for (int j=0;j<ncomp_mdot_alloc;j++) {
    level_mdot_data[i][j]=level_mdot_array[0][mdot_counter];
    level_mdot_comp_data[i][j]=level_mdot_comp_array[0][mdot_counter];
    mdot_counter++;
   }

  }  // i=0..num_colors-1

  if (counter!=blob_array_size)
   amrex::Error("counter invalid");
  if (mdot_counter!=mdot_array_size) {
   std::cout << "mdot_counter=" << mdot_counter << '\n';
   std::cout << "mdot_array_size=" << mdot_array_size << '\n';
   amrex::Error("mdot_counter invalid in ColorSum3");
  }

 } else if (operation_flag==OP_SCATTER_MDOT) {

  for (int i=0;i<mdot_array_size;i++) {
   ParallelDescriptor::ReduceRealSum(level_mdot_redistribute_array[0][i]);
   ParallelDescriptor::ReduceRealSum(level_mdot_comp_redistribute_array[0][i]);
  }

  mdot_counter=0;
  for (int i=0;i<num_colors;i++) {
   for (int j=0;j<ncomp_mdot_alloc;j++) {
    level_mdot_data_redistribute[i][j]=
      level_mdot_redistribute_array[0][mdot_counter];
    level_mdot_comp_data_redistribute[i][j]=
      level_mdot_comp_redistribute_array[0][mdot_counter];
    mdot_counter++;
   }
  }
  if (mdot_counter!=mdot_array_size) {
   std::cout << "mdot_counter=" << mdot_counter << '\n';
   std::cout << "mdot_array_size=" << mdot_array_size << '\n';
   amrex::Error("mdot_counter invalid in ColorSum4");
  }

 } else
  amrex::Error("operation_flag invalid");

 delete mask;

 if ((sweep_num==0)&&(operation_flag==OP_GATHER_MDOT)) {
  // do nothing
 } else if ((sweep_num==1)||(operation_flag==OP_SCATTER_MDOT)) {
  delete_localMF(LS_COLORSUM_MF,1);
  delete_localMF(DEN_COLORSUM_MF,1);
  delete_localMF(VEL_COLORSUM_MF,1);
 } else
  amrex::Error("sweep_num or operation_flag invalid");

}  // end subroutine ColorSum

void
NavierStokes::SumRegions(
 int isweep) { // isweep=0 or 1

 std::string local_caller_string="SumRegions";

 int finest_level=parent->finestLevel();
 bool use_tiling=ns_tiling;

 if (level>finest_level)
  amrex::Error("level invalid SumRegions");

 if (MDOT_MF>=0) {
  if (localMF[MDOT_MF]->nComp()==1) {
   // do nothing
  } else
   amrex::Error("localMF[MDOT_MF]->nComp() invalid");
  if (localMF[MDOT_MF]->nGrow()>=0) {
   // do nothing
  } else
   amrex::Error("localMF[MDOT_MF]->nGrow() invalid");

 } else
  amrex::Error("MDOT_MF invalid");

 resize_metrics(1);

 debug_ngrow(VOLUME_MF,0,local_caller_string);
 VOF_Recon_resize(1); //output:SLOPE_RECON_MF
 debug_ngrow(SLOPE_RECON_MF,1,local_caller_string);
 if (localMF[SLOPE_RECON_MF]->nComp()!=num_materials*ngeom_recon)
  amrex::Error("localMF[SLOPE_RECON_MF]->nComp() invalid");

 if (isweep==0) {
  getStateDen_localMF(DEN_COLORSUM_MF,1,cur_time_slab);
 } else if (isweep==1) {
  // do nothing
 } else
  amrex::Error("isweep invalid");

 if (localMF[DEN_COLORSUM_MF]->nComp()!=num_state_material*num_materials)
  amrex::Error("localMF[DEN_COLORSUM_MF]->nComp()!=num_state_material*num_materials");

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 MultiFab& Umac_new=get_new_data(Umac_Type,slab_step+1);
 MultiFab& Vmac_new=get_new_data(Umac_Type+1,slab_step+1);
 MultiFab& Wmac_new=get_new_data(Umac_Type+AMREX_SPACEDIM-1,slab_step+1);

 int nstate=STATE_NCOMP;
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");

  // mask=tag if not covered by level+1 and at fine/fine ghost cell.
 int ngrowmask=1;
 Real tag=1.0;
 int clear_phys_boundary=2;
 MultiFab* mask=maskfiner(ngrowmask,tag,clear_phys_boundary);  
 const Real* dx = geom.CellSize();

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

  FArrayBox& snewfab=S_new[mfi];

  FArrayBox& volumefab=(*localMF[VOLUME_MF])[mfi];
  FArrayBox& mdotfab=(*localMF[MDOT_MF])[mfi];
  FArrayBox& denfab=(*localMF[DEN_COLORSUM_MF])[mfi];
  FArrayBox& voffab=(*localMF[SLOPE_RECON_MF])[mfi];

  FArrayBox& maskfab=(*mask)[mfi];

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // in: NAVIERSTOKES_3D.F90
  fort_regionsum(
   &tid_current,
   &isweep,  //isweep=0 or 1
   constant_density_all_time.dataPtr(),
   &cur_time_slab,
   &dt_slab,
   dx,
   xlo,
   &nstate,
   snewfab.dataPtr(),ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()),
   Umac_new[mfi].dataPtr(),
   ARLIM(Umac_new[mfi].loVect()),ARLIM(Umac_new[mfi].hiVect()),
   Vmac_new[mfi].dataPtr(),
   ARLIM(Vmac_new[mfi].loVect()),ARLIM(Vmac_new[mfi].hiVect()),
   Wmac_new[mfi].dataPtr(),
   ARLIM(Wmac_new[mfi].loVect()),ARLIM(Wmac_new[mfi].hiVect()),
   mdotfab.dataPtr(),ARLIM(mdotfab.loVect()),ARLIM(mdotfab.hiVect()),
   denfab.dataPtr(),ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
   voffab.dataPtr(),ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),
   volumefab.dataPtr(),ARLIM(volumefab.loVect()),ARLIM(volumefab.hiVect()),
   maskfab.dataPtr(),ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
   tilelo,tilehi,
   fablo,fabhi,
   &bfact,
   &level,
   &finest_level);
 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_REGIONSUM,"SumRegions");

 ParallelDescriptor::Barrier();

 delete mask;

 if (isweep==0) {
  // do nothing
 } else if (isweep==1) {
  delete_localMF(DEN_COLORSUM_MF,1);
 } else
  amrex::Error("isweep invalid");

}  //end subroutine SumRegions

void
NavierStokes::LowMachDIVU(
 int sweep_num,
 MultiFab* typemf,
 MultiFab* color,
 MultiFab* mdot_local, 
 MultiFab* mdot_global, 
 Vector<blobclass> cum_blobdata,
 Vector< Vector<Real> >& level_mdot_data,
 Vector< Vector<Real> >& cum_mdot_data
 ) {

 std::string local_caller_string="LowMachDIVU";

 int finest_level=parent->finestLevel();
 bool use_tiling=ns_tiling;

 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid LowMachDIVU");

 if (ngrow_distance!=4)
  amrex::Error("ngrow_distance invalid");

 int num_colors=cum_blobdata.size();
 if (num_colors==0) 
  amrex::Error("num_colors should be positive in LowMachDIVU");

 if (num_colors!=cum_mdot_data.size())
  amrex::Error("num_colors!=cum_mdot_data.size()");
 if (num_colors!=level_mdot_data.size())
  amrex::Error("num_colors!=level_mdot_data.size()");

 int ncomp_mdot_test=mdot_local->nComp();
 if (ncomp_mdot_test==1) {
  // do nothing
 } else
  amrex::Error("ncomp_mdot_test invalid");

 if (mdot_global->nComp()==1) {
  // do nothing
 } else
  amrex::Error("mdot_global->nComp() invalid");

 for (int i=0;i<num_colors;i++) {
  for (int j=0;j<2;j++) {
   level_mdot_data[i][j]=0.0;
  }
 } // i=0..num_colors-1

 int blob_array_size=num_colors*num_elements_blobclass;

 int mdot_array_size=num_colors*2;

 Vector<Real> cum_blob_array;
 cum_blob_array.resize(blob_array_size);

 Vector<Real> cum_mdot_array;
 cum_mdot_array.resize(mdot_array_size);

 int counter=0;
 int mdot_counter=0;

 for (int i=0;i<num_colors;i++) {
  copy_from_blobdata(i,counter,cum_blob_array,cum_blobdata);
  for (int j=0;j<2;j++) {
   cum_mdot_array[mdot_counter]=cum_mdot_data[i][j];
   mdot_counter++;
  }
 }  // i=0..num_colors-1

 if (counter!=blob_array_size)
  amrex::Error("counter invalid");

 if (mdot_counter!=mdot_array_size) {
  std::cout << "mdot_counter=" << mdot_counter << '\n';
  std::cout << "mdot_array_size=" << mdot_array_size << '\n';
  amrex::Error("mdot_counter invalid in LowMachDIVU 1");
 }

 Vector< Vector<Real> > level_mdot_array;
 level_mdot_array.resize(thread_class::nthreads);

 for (int tid=0;tid<thread_class::nthreads;tid++) {

   //mdot_array_size=num_colors*2
  level_mdot_array[tid].resize(mdot_array_size);

  for (int i=0;i<mdot_array_size;i++) {
   if (sweep_num==0) {
    level_mdot_array[tid][i]=0.0;
   } else if (sweep_num==1) {
    level_mdot_array[tid][i]=cum_mdot_array[i];
   } else
    amrex::Error("sweep_num invalid");
  } // i=0;i<mdot_array_size;i++

 } // tid=0..thread_class::nthreads-1

 resize_metrics(1);

 debug_ngrow(VOLUME_MF,0,local_caller_string);
 VOF_Recon_resize(1); //output:SLOPE_RECON_MF
 debug_ngrow(SLOPE_RECON_MF,1,local_caller_string);
 if (localMF[SLOPE_RECON_MF]->nComp()!=num_materials*ngeom_recon)
  amrex::Error("localMF[SLOPE_RECON_MF]->nComp() invalid");

 if (sweep_num==0) {
  getStateDist_localMF(LS_COLORSUM_MF,1,cur_time_slab,local_caller_string);
  getStateDen_localMF(DEN_COLORSUM_MF,1,cur_time_slab);
 } else if (sweep_num==1) {
  // do nothing
 } else
  amrex::Error("sweep_num invalid");

 if (localMF[LS_COLORSUM_MF]->nGrow()!=1)
  amrex::Error("localMF[LS_COLORSUM_MF]->nGrow()!=1");

 if (localMF[LS_COLORSUM_MF]->nComp()!=(1+AMREX_SPACEDIM)*num_materials)
  amrex::Error("localMF[LS_COLORSUM_MF]->nComp()!=(1+AMREX_SPACEDIM)*num_materials");
 if (localMF[DEN_COLORSUM_MF]->nComp()!=num_state_material*num_materials)
  amrex::Error("localMF[DEN_COLORSUM_MF]->nComp()!=num_state_material*num_materials");

 if (typemf->nGrow()!=1)
  amrex::Error("typemf->nGrow()!=1");
 if (color->nGrow()!=1)
  amrex::Error("color->nGrow()!=1");
 if (mdot_local->nGrow()>=0) {
  // do nothing
 } else
  amrex::Error("mdot_local->nGrow() invalid");
 if (mdot_global->nGrow()>=0) {
  // do nothing
 } else
  amrex::Error("mdot_global->nGrow() invalid");

  // mask=tag if not covered by level+1 and at fine/fine ghost cell.
 int ngrowmask=1;
 Real tag=1.0;
 int clear_phys_boundary=2;
 MultiFab* mask=maskfiner(ngrowmask,tag,clear_phys_boundary);  
 const Real* dx = geom.CellSize();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(mask->boxArray().d_numPts());

// COLORING LOOP
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

  FArrayBox& mdot_local_fab=(*mdot_local)[mfi];
  FArrayBox& mdot_global_fab=(*mdot_global)[mfi];

  FArrayBox& typefab=(*typemf)[mfi];
  FArrayBox& lsfab=(*localMF[LS_COLORSUM_MF])[mfi];
  FArrayBox& denfab=(*localMF[DEN_COLORSUM_MF])[mfi];
  FArrayBox& DTDtfab=(*localMF[DTDt_MF])[mfi];
  FArrayBox& voffab=(*localMF[SLOPE_RECON_MF])[mfi];

  FArrayBox& colorfab=(*color)[mfi];
  FArrayBox& maskfab=(*mask)[mfi];

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // in: LEVELSET_3D.F90
  fort_get_lowmach_divu(
   &tid_current,
   &sweep_num,
   constant_density_all_time.dataPtr(),
   &dt_slab,
   dx,xlo,
   mdot_local_fab.dataPtr(),
   ARLIM(mdot_local_fab.loVect()),
   ARLIM(mdot_local_fab.hiVect()),
   mdot_global_fab.dataPtr(),
   ARLIM(mdot_global_fab.loVect()),
   ARLIM(mdot_global_fab.hiVect()),
   lsfab.dataPtr(),ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
   denfab.dataPtr(),ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
   DTDtfab.dataPtr(),ARLIM(DTDtfab.loVect()),ARLIM(DTDtfab.hiVect()),
   voffab.dataPtr(),ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),
   typefab.dataPtr(),ARLIM(typefab.loVect()),ARLIM(typefab.hiVect()),
   colorfab.dataPtr(),ARLIM(colorfab.loVect()),ARLIM(colorfab.hiVect()),
   maskfab.dataPtr(),ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
   tilelo,tilehi,
   fablo,fabhi,
   &bfact,
   &level,
   &finest_level,
   &NS_geometry_coord,
   &num_colors,
   cum_blob_array.dataPtr(),
   level_mdot_array[tid_current].dataPtr(),
   &blob_array_size,
   &mdot_array_size,
   material_type_lowmach.dataPtr());
 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_LOWMACH_DIVU,"LowMachDIVU");

if (sweep_num==0) {
 for (int tid=1;tid<thread_class::nthreads;tid++) {
  for (int i=0;i<mdot_array_size;i++) {
   level_mdot_array[0][i]+=level_mdot_array[tid][i];
  }
 } // tid
} else if (sweep_num==1) {
 // do nothing
} else
 amrex::Error("sweep_num invalid");

ParallelDescriptor::Barrier();

if (sweep_num==0) {
 for (int i=0;i<mdot_array_size;i++) {
  ParallelDescriptor::ReduceRealSum(level_mdot_array[0][i]);
 }
} else if (sweep_num==1) {
 // do nothing
} else
 amrex::Error("sweep_num invalid");

counter=0;
mdot_counter=0;

for (int i=0;i<num_colors;i++) {
 for (int j=0;j<2;j++) {
  level_mdot_data[i][j]=level_mdot_array[0][mdot_counter];
  mdot_counter++;
 }
}  // i=0..num_colors-1

if (mdot_counter!=mdot_array_size) {
 std::cout << "mdot_counter=" << mdot_counter << '\n';
 std::cout << "mdot_array_size=" << mdot_array_size << '\n';
 amrex::Error("mdot_counter invalid in LowMachDIVU 2");
}

delete mask;

if (sweep_num==0) {
  // do nothing
} else if (sweep_num==1) {
 delete_localMF(LS_COLORSUM_MF,1);
 delete_localMF(DEN_COLORSUM_MF,1);
} else
 amrex::Error("sweep_num invalid");

}  // end subroutine LowMachDIVU

void
NavierStokes::LowMachDIVUALL(
 int coarsest_level,
 int& color_count,
 int idx_type,
 int idx_color,
 Vector<int>& /* type_flag */, 
 Vector<blobclass> blobdata) {

 std::string local_caller_string="LowMachDIVUALL";

 int finest_level=parent->finestLevel();

 if ((coarsest_level<0)||(coarsest_level>finest_level))
  amrex::Error("coarsest_level invalid");

 if (level!=0)
  amrex::Error("level=0 in LowMachDIVUALL");

 if (color_count==0)
  amrex::Error("num_colors=0 in LowMachDIVUALL");

   // initializes MDOT_LOCAL_MF to 0.0
   // ngrow=0, ncomp=1, grid_type=-1
 allocate_array(0,1,-1,MDOT_LOCAL_MF); 

 if (MDOT_LOCAL_MF>=0) {
  int ncomp_mdot=localMF[MDOT_LOCAL_MF]->nComp();
  if (ncomp_mdot==1) {
   // do nothing
  } else
   amrex::Error("ncomp_mdot invalid");
  int ngrow_mdot=localMF[MDOT_LOCAL_MF]->nGrow();
  if (ngrow_mdot>=0) {
   // do nothing
  } else
   amrex::Error("ngrow_mdot invalid");

 } else
  amrex::Error("MDOT_LOCAL_MF invalid");


 if (MDOT_MF>=0) {
  if (localMF[MDOT_MF]->nComp()==1) {
   // do nothing
  } else
   amrex::Error("localMF[MDOT_MF]->nComp() invalid");
  if (localMF[MDOT_MF]->nGrow()>=0) {
   // do nothing
  } else
   amrex::Error("localMF[MDOT_MF]->nGrow() invalid");

 } else
  amrex::Error("MDOT_MF invalid");

  //mdot_data[icolor][j=0 or 1]
 Vector< Vector<Real> > mdot_data;

 mdot_data.resize(color_count);
 for (int i=0;i<color_count;i++) {
  mdot_data[i].resize(2);
  int j=0;
  for (j=0;j<2;j++) {
   mdot_data[i][j]=0.0;
  }

  if (MDOT_LOCAL_MF>=0) {
   if (j==2) {
     // do nothing
   } else
    amrex::Error("expecting j==2");
  } else
   amrex::Error("MDOT_LOCAL_MF invalid");

 }  // i=0..color_count-1

  //level_mdot_data[icolor][j=0 or 1]
 Vector< Vector<Real> > level_mdot_data;
 level_mdot_data.resize(color_count);

 for (int i=0;i<color_count;i++) {

  level_mdot_data[i].resize(2);

  int j=0;
  for (j=0;j<2;j++) {
   level_mdot_data[i][j]=0.0;
  }

  if (MDOT_LOCAL_MF>=0) {
   if (j==2) {
    // do nothing
   } else
    amrex::Error("expecting j==2");
  } else
   amrex::Error("MDOT_LOCAL_MF invalid");

 } // i=0..color_count-1

 int num_sweeps=2;

 for (int sweep_num=0;sweep_num<num_sweeps;sweep_num++) {

  for (int ilev = coarsest_level; ilev <= finest_level; ilev++) {

   NavierStokes& ns_level = getLevel(ilev);

   MultiFab* mdot_local=ns_level.localMF[MDOT_LOCAL_MF];
   MultiFab* mdot_global=ns_level.localMF[MDOT_MF];

   ns_level.LowMachDIVU(
    sweep_num,
    ns_level.localMF[idx_type],
    ns_level.localMF[idx_color],
    mdot_local,
    mdot_global,
    blobdata,
    level_mdot_data,
    mdot_data
    );

   if (sweep_num==0) {

    for (int i=0;i<color_count;i++) {

     int j=0;
     for (j=0;j<2;j++) {
      mdot_data[i][j]+=level_mdot_data[i][j];
     }
     if (MDOT_LOCAL_MF>=0) {
      if (j==2) {
       // do nothing
      } else
       amrex::Error("expecting j==2");
     } else
      amrex::Error("MDOT_LOCAL_MF invalid");

    }  // i=0..color_count-1

   } else if (sweep_num==1) {

    // do nothing

   } else
     amrex::Error("sweep_num invalid");

  } // ilev=coarsest_level..finest_level

 } // sweep_num=0..1

 if (1==0) {
  writeSanityCheckData(
   "MDOT_LOCAL",
   "in: NavierStokes::LowMachDIVUALL", 
   local_caller_string,
   MDOT_LOCAL_MF,  //tower_mf_id
   localMF[MDOT_LOCAL_MF]->nComp(),
   MDOT_LOCAL_MF, 
   -1, //State_Type==-1
   -1, // data_dir==-1 (cell centered)
   parent->levelSteps(0)); 
 }

 delete_array(MDOT_LOCAL_MF);

} // end subroutine LowMachDIVUALL


void NavierStokes::copy_to_blobdata(int i,int& counter,
  Vector<Real>& blob_array,Vector<blobclass>& blobdata) {

 int num_colors=blob_array.size()/num_elements_blobclass;
 if (blob_array.size()!=num_colors*num_elements_blobclass)
  amrex::Error("num_colors invalid");
 if (num_colors!=blobdata.size())
  amrex::Error("blobdata.size() invalid");
 if ((i<0)||(i>=num_colors))
  amrex::Error("i invalid");
 if (counter!=i*num_elements_blobclass)
  amrex::Error("counter invalid");
 for (int dir=0;dir<3*(2*AMREX_SPACEDIM)*(2*AMREX_SPACEDIM);dir++) {
  blobdata[i].blob_matrix[dir]=blob_array[counter+BLB_MATRIX+dir];
 }
 for (int dir=0;dir<3*(2*AMREX_SPACEDIM);dir++) {
  blobdata[i].blob_RHS[dir]=blob_array[counter+BLB_RHS+dir];
 }
 for (int dir=0;dir<3*(2*AMREX_SPACEDIM);dir++) {
  blobdata[i].blob_velocity[dir]=blob_array[counter+BLB_VEL+dir];
 }
 for (int dir=0;dir<2*(2*AMREX_SPACEDIM);dir++) {
  blobdata[i].blob_integral_momentum[dir]=
	  blob_array[counter+BLB_INT_MOM+dir];
 }
 blobdata[i].blob_energy=blob_array[counter+BLB_ENERGY];

 for (int dir=0;dir<3;dir++) {
  blobdata[i].blob_mass_for_velocity[dir]=
	  blob_array[counter+BLB_MASS_VEL+dir];
 }
 blobdata[i].blob_volume=blob_array[counter+BLB_VOL];

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  blobdata[i].blob_center_integral[dir]=
	  blob_array[counter+BLB_CEN_INT+dir];
 }
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  blobdata[i].blob_center_actual[dir]=
	  blob_array[counter+BLB_CEN_ACT+dir];
 }
 blobdata[i].blob_perim=blob_array[counter+BLB_PERIM];

 for (int imnbr=0;imnbr<num_materials;imnbr++) {
  blobdata[i].blob_perim_mat[imnbr]=blob_array[counter+BLB_PERIM_MAT+imnbr];
 }
 for (int im1=0;im1<num_materials;im1++) {
  for (int im2=0;im2<num_materials;im2++) {
   blobdata[i].blob_triple_perim[im1][im2]=
     blob_array[counter+BLB_TRIPLE_PERIM+num_materials*im1+im2];
  } // im2
 } // im1

 blobdata[i].blob_cell_count=blob_array[counter+BLB_CELL_CNT];
 blobdata[i].blob_cellvol_count=blob_array[counter+BLB_CELLVOL_CNT];
 blobdata[i].blob_mass=blob_array[counter+BLB_MASS];
 blobdata[i].blob_pressure=blob_array[counter+BLB_PRES];

 for (int dir=0;dir<6;dir++) {
  blobdata[i].blob_second_moment[dir]= 
	blob_array[counter+BLB_SECONDMOMENT+dir];
 }

 counter+=num_elements_blobclass;

} // end subroutine copy_to_blobdata


void NavierStokes::copy_blobdata(Vector<blobclass>& dest_blobdata,
  Vector<blobclass>& source_blobdata) {

 int num_colors=source_blobdata.size();
 int num_colors_test=dest_blobdata.size();

 if ((num_colors>0)&&
     (num_colors==num_colors_test)) {

  for (int i=0;i<num_colors;i++) {

   dest_blobdata[i].im=source_blobdata[i].im;

   for (int dir=0;dir<3*(2*AMREX_SPACEDIM)*(2*AMREX_SPACEDIM);dir++) {
    dest_blobdata[i].blob_matrix[dir]=
      source_blobdata[i].blob_matrix[dir];
   }
   for (int dir=0;dir<3*(2*AMREX_SPACEDIM);dir++) {
    dest_blobdata[i].blob_RHS[dir]=source_blobdata[i].blob_RHS[dir];
   }
   for (int dir=0;dir<3*(2*AMREX_SPACEDIM);dir++) {
    dest_blobdata[i].blob_velocity[dir]=
      source_blobdata[i].blob_velocity[dir];
   }
   for (int dir=0;dir<2*(2*AMREX_SPACEDIM);dir++) {
    dest_blobdata[i].blob_integral_momentum[dir]=
      source_blobdata[i].blob_integral_momentum[dir];
   }
   dest_blobdata[i].blob_energy=source_blobdata[i].blob_energy;
   for (int dir=0;dir<3;dir++) {
    dest_blobdata[i].blob_mass_for_velocity[dir]=
      source_blobdata[i].blob_mass_for_velocity[dir];
   }
   dest_blobdata[i].blob_volume=
      source_blobdata[i].blob_volume;
   for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
    dest_blobdata[i].blob_center_integral[dir]=
      source_blobdata[i].blob_center_integral[dir];
   }
   for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
    dest_blobdata[i].blob_center_actual[dir]=
      source_blobdata[i].blob_center_actual[dir];
   }
   dest_blobdata[i].blob_perim=
     source_blobdata[i].blob_perim;
   for (int imnbr=0;imnbr<num_materials;imnbr++) {
    dest_blobdata[i].blob_perim_mat[imnbr]=
      source_blobdata[i].blob_perim_mat[imnbr];
   }
   for (int im1=0;im1<num_materials;im1++) {
    for (int im2=0;im2<num_materials;im2++) {
     dest_blobdata[i].blob_triple_perim[im1][im2]=
       source_blobdata[i].blob_triple_perim[im1][im2];
    } // im2
   } // im1

   dest_blobdata[i].blob_cell_count=
     source_blobdata[i].blob_cell_count;

   dest_blobdata[i].blob_cellvol_count=
     source_blobdata[i].blob_cellvol_count;

   dest_blobdata[i].blob_mass=
     source_blobdata[i].blob_mass;

   dest_blobdata[i].blob_pressure=
     source_blobdata[i].blob_pressure;

   for (int dir=0;dir<6;dir++) {
    dest_blobdata[i].blob_second_moment[dir]=
      source_blobdata[i].blob_second_moment[dir];
   }

  } // i=0..num_colors-1

 } else
  amrex::Error("num_colors invalid");

} // end subroutine copy_blobdata


void NavierStokes::sum_blobdata(int i,
  Vector<blobclass>& blobdata,
  Vector<blobclass>& level_blobdata,int sweep_num) {

 int num_colors=blobdata.size();
 if ((i<0)||(i>=num_colors))
  amrex::Error("i invalid");
 if (blobdata.size()!=level_blobdata.size())
  amrex::Error("level_blobdata.size() invalid");

 if (sweep_num==0) {

  blobdata[i].blob_volume+=level_blobdata[i].blob_volume;
  blobdata[i].blob_cell_count+=level_blobdata[i].blob_cell_count;
  blobdata[i].blob_cellvol_count+=level_blobdata[i].blob_cellvol_count;
  blobdata[i].blob_mass+=level_blobdata[i].blob_mass;
  blobdata[i].blob_pressure+=level_blobdata[i].blob_pressure;

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   blobdata[i].blob_center_integral[dir]+=
    level_blobdata[i].blob_center_integral[dir];
  }
  blobdata[i].blob_perim+=level_blobdata[i].blob_perim;

  for (int imnbr=0;imnbr<num_materials;imnbr++)
   blobdata[i].blob_perim_mat[imnbr]+=
    level_blobdata[i].blob_perim_mat[imnbr];

  for (int im1=0;im1<num_materials;im1++)
   for (int im2=0;im2<num_materials;im2++)
    blobdata[i].blob_triple_perim[im1][im2]+=
     level_blobdata[i].blob_triple_perim[im1][im2];

 } else if (sweep_num==1) {

  blobdata[i].blob_energy+=level_blobdata[i].blob_energy;

  for (int dir=0;dir<3*(2*AMREX_SPACEDIM)*(2*AMREX_SPACEDIM);dir++) 
   blobdata[i].blob_matrix[dir]+=
    level_blobdata[i].blob_matrix[dir];
  for (int dir=0;dir<3*(2*AMREX_SPACEDIM);dir++) 
   blobdata[i].blob_RHS[dir]+=
    level_blobdata[i].blob_RHS[dir];
  for (int dir=0;dir<2*(2*AMREX_SPACEDIM);dir++) 
   blobdata[i].blob_integral_momentum[dir]+=
    level_blobdata[i].blob_integral_momentum[dir];
  for (int dir=0;dir<3;dir++) {
   blobdata[i].blob_mass_for_velocity[dir]+=
    level_blobdata[i].blob_mass_for_velocity[dir];
  }
  for (int dir=0;dir<6;dir++) {
   blobdata[i].blob_second_moment[dir]+=
    level_blobdata[i].blob_second_moment[dir];
  }

 } else
  amrex::Error("sweep_num invalid");

} // end subroutine sum_blobdata

void NavierStokes::copy_from_blobdata(int i,int& counter,
  Vector<Real>& blob_array,Vector<blobclass>& blobdata) {

 int num_colors=blob_array.size()/num_elements_blobclass;
 if (blob_array.size()!=num_colors*num_elements_blobclass)
  amrex::Error("num_colors invalid");
 if (num_colors!=blobdata.size())
  amrex::Error("blobdata.size() invalid");
 if ((i<0)||(i>=num_colors))
  amrex::Error("i invalid");
 if (counter!=i*num_elements_blobclass)
  amrex::Error("counter invalid");

 for (int dir=0;dir<3*(2*AMREX_SPACEDIM)*(2*AMREX_SPACEDIM);dir++) {
  blob_array[counter+BLB_MATRIX+dir]=blobdata[i].blob_matrix[dir];
 }
 for (int dir=0;dir<3*(2*AMREX_SPACEDIM);dir++) {
  blob_array[counter+BLB_RHS+dir]=blobdata[i].blob_RHS[dir];
 }
 for (int dir=0;dir<3*(2*AMREX_SPACEDIM);dir++) {
  blob_array[counter+BLB_VEL+dir]=blobdata[i].blob_velocity[dir];
 }
 for (int dir=0;dir<2*(2*AMREX_SPACEDIM);dir++) {
  blob_array[counter+BLB_INT_MOM+dir]=
	  blobdata[i].blob_integral_momentum[dir];
 }
 blob_array[counter+BLB_ENERGY]=blobdata[i].blob_energy;

 for (int dir=0;dir<3;dir++) {
  blob_array[counter+BLB_MASS_VEL+dir]=
	  blobdata[i].blob_mass_for_velocity[dir];
 }
 blob_array[counter+BLB_VOL]=blobdata[i].blob_volume;

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  blob_array[counter+BLB_CEN_INT+dir]=blobdata[i].blob_center_integral[dir];
 }
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  blob_array[counter+BLB_CEN_ACT+dir]=blobdata[i].blob_center_actual[dir];
 }
 blob_array[counter+BLB_PERIM]=blobdata[i].blob_perim;

 for (int imnbr=0;imnbr<num_materials;imnbr++) {
  blob_array[counter+BLB_PERIM_MAT+imnbr]=blobdata[i].blob_perim_mat[imnbr];
 }
 for (int im1=0;im1<num_materials;im1++) {
  for (int im2=0;im2<num_materials;im2++) {
   blob_array[counter+BLB_TRIPLE_PERIM+num_materials*im1+im2]=
     blobdata[i].blob_triple_perim[im1][im2];
  } // im2
 } // im1
 blob_array[counter+BLB_CELL_CNT]=blobdata[i].blob_cell_count;
 blob_array[counter+BLB_CELLVOL_CNT]=blobdata[i].blob_cellvol_count;
 blob_array[counter+BLB_MASS]=blobdata[i].blob_mass;
 blob_array[counter+BLB_PRES]=blobdata[i].blob_pressure;
 for (int dir=0;dir<6;dir++) {
  blob_array[counter+BLB_SECONDMOMENT+dir]=
     blobdata[i].blob_second_moment[dir];
 }

 counter+=num_elements_blobclass;

}  // end subroutine copy_from_blobdata

void NavierStokes::clear_blobdata(int i,Vector<blobclass>& blobdata) {

 int num_colors=blobdata.size();
 if ((i<0)||(i>=num_colors))
  amrex::Error("i invalid");

 for (int dir=0;dir<3*(2*AMREX_SPACEDIM)*(2*AMREX_SPACEDIM);dir++)
  blobdata[i].blob_matrix[dir]=0.0;
 for (int dir=0;dir<3*(2*AMREX_SPACEDIM);dir++)
  blobdata[i].blob_RHS[dir]=0.0;
 for (int dir=0;dir<3*(2*AMREX_SPACEDIM);dir++)
  blobdata[i].blob_velocity[dir]=0.0;
 for (int dir=0;dir<2*(2*AMREX_SPACEDIM);dir++)
  blobdata[i].blob_integral_momentum[dir]=0.0;
 blobdata[i].blob_energy=0.0;
 for (int dir=0;dir<3;dir++)
  blobdata[i].blob_mass_for_velocity[dir]=0.0;

 blobdata[i].blob_volume=0.0;
 blobdata[i].blob_cell_count=0.0;
 blobdata[i].blob_cellvol_count=0.0;
 blobdata[i].blob_mass=0.0;
 blobdata[i].blob_pressure=0.0;

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  blobdata[i].blob_center_integral[dir]=0.0;
  blobdata[i].blob_center_actual[dir]=0.0;
 }
 blobdata[i].blob_perim=0.0;
 blobdata[i].blob_perim_mat.resize(num_materials);
 blobdata[i].blob_triple_perim.resize(num_materials);
 for (int imnbr=0;imnbr<num_materials;imnbr++) {
  blobdata[i].blob_triple_perim[imnbr].resize(num_materials);
  blobdata[i].blob_perim_mat[imnbr]=0.0;
 }
 for (int im1=0;im1<num_materials;im1++)
  for (int im2=0;im2<num_materials;im2++)
   blobdata[i].blob_triple_perim[im1][im2]=0.0;

 for (int dir=0;dir<6;dir++) {
  blobdata[i].blob_second_moment[dir]=0.0;
 }

 blobdata[i].im=0;

} // end subroutine clear_blobdata

//operation_flag==1 (OP_SCATTER_MDOT) => scatter data collected when 
//operation_flag==0 (OP_GATHER_MDOT) to mdot or density.
void
NavierStokes::ColorSumALL(
 int operation_flag, //OP_GATHER_MDOT or OP_SCATTER_MDOT
 int tessellate,  // 1 or 3
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
 ) {

 int finest_level=parent->finestLevel();

 if ((coarsest_level<0)||(coarsest_level>finest_level))
  amrex::Error("coarsest_level invalid");

 if (level!=0)
  amrex::Error("level=0 in ColorSumALL");

 if (operation_flag==OP_GATHER_MDOT) {

  // type_flag[im]=1 if material im exists in the domain.
  // type_mf(i,j,k)=im if material im dominates cell (i,j,k)
  // updates one ghost cell of TYPE_MF
  // fluid(s) and solid(s) tessellate the domain.
  int zero_diag_flag=0;
  TypeALL(idx_type,type_flag,zero_diag_flag);

  // color_count=number of colors
  // ngrow=1, fort_extrapfill, pc_interp for COLOR_MF
  color_variable(coarsest_level,
   idx_color,idx_type,&color_count,type_flag,zero_diag_flag);

 } else if (operation_flag==OP_SCATTER_MDOT) {

  // do nothing

 } else
  amrex::Error("operation_flag invalid");

 Real problo_array[AMREX_SPACEDIM];
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  problo_array[dir]=geom.ProbLo(dir);
 }

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  if (force_blob_symmetry[dir]==1) {
   // do nothing
  } else if (force_blob_symmetry[dir]==0) {
   // do nothing
  } else
   amrex::Error("force_blob_symmetry[dir] out of range");
 } //dir=0..sdim-1

 if (geom.IsRZ()) {
  if (force_blob_symmetry[0]==1) {
   // do nothing
  } else if (force_blob_symmetry[0]==0) {
   amrex::Error("need force_blob_symmetry[0]==1 if geom.IsRZ()");
  } else
   amrex::Error("force_blob_symmetry[0] out of range");
 }  // IsRZ?

 if (color_count==0)
  amrex::Error("num_colors=0 in ColorSumALL");

 int ncomp_mdot_alloc=1;
 int ncomp_mdot=0;

 if (idx_mdot==-1) {
  // do nothing
 } else if (idx_mdot>=0) {
  ncomp_mdot=localMF[idx_mdot]->nComp();
  ncomp_mdot_alloc=ncomp_mdot;
  if (ncomp_mdot==2*num_interfaces) {
   // do nothing
  } else
   amrex::Error("ncomp_mdot invalid");
  int ngrow_mdot=localMF[idx_mdot]->nGrow();
  if (ngrow_mdot>=0) {
   // do nothing
  } else
   amrex::Error("ngrow_mdot invalid");

 } else
  amrex::Error("idx_mdot invalid");

 if (operation_flag==OP_GATHER_MDOT) {

  blobdata.resize(color_count);
  mdot_data.resize(color_count);
  mdot_data_redistribute.resize(color_count);
  mdot_comp_data.resize(color_count);
  mdot_comp_data_redistribute.resize(color_count);
  for (int i=0;i<color_count;i++) {
   clear_blobdata(i,blobdata);
   mdot_data[i].resize(ncomp_mdot_alloc);
   mdot_data_redistribute[i].resize(ncomp_mdot_alloc);
   mdot_comp_data[i].resize(ncomp_mdot_alloc);
   mdot_comp_data_redistribute[i].resize(ncomp_mdot_alloc);
   int j=0;
   for (j=0;j<ncomp_mdot_alloc;j++) {
    mdot_data[i][j]=0.0;
    mdot_data_redistribute[i][j]=0.0;
    mdot_comp_data[i][j]=0.0;
    mdot_comp_data_redistribute[i][j]=0.0;
   }

   if (idx_mdot>=0) {
    if (j==2*num_interfaces) {
     // do nothing
    } else
     amrex::Error("expecting j==2*num_interfaces");
   } else if (idx_mdot==-1) {
    // check nothing
   } else
    amrex::Error("idx_mdot invalid");

  }  // i=0..color_count-1

 } else if (operation_flag==OP_SCATTER_MDOT) {

  // do nothing

 } else
  amrex::Error("operation_flag invalid");

 Vector<blobclass> level_blobdata;
 level_blobdata.resize(color_count);

 Vector< Vector<Real> > level_mdot_data;
 level_mdot_data.resize(color_count);
 Vector< Vector<Real> > level_mdot_data_redistribute;
 level_mdot_data_redistribute.resize(color_count);

 Vector< Vector<Real> > level_mdot_comp_data;
 level_mdot_comp_data.resize(color_count);
 Vector< Vector<Real> > level_mdot_comp_data_redistribute;
 level_mdot_comp_data_redistribute.resize(color_count);

 for (int i=0;i<color_count;i++) {
  clear_blobdata(i,level_blobdata); 

  level_mdot_data[i].resize(ncomp_mdot_alloc);
  level_mdot_data_redistribute[i].resize(ncomp_mdot_alloc);

  level_mdot_comp_data[i].resize(ncomp_mdot_alloc);
  level_mdot_comp_data_redistribute[i].resize(ncomp_mdot_alloc);

  int j=0;
  for (j=0;j<ncomp_mdot_alloc;j++) {
   level_mdot_data[i][j]=0.0;
   level_mdot_data_redistribute[i][j]=0.0;

   level_mdot_comp_data[i][j]=0.0;
   level_mdot_comp_data_redistribute[i][j]=0.0;
  }

  if (idx_mdot>=0) {
   if (j==2*num_interfaces) {
    // do nothing
   } else
    amrex::Error("expecting j==2*num_interfaces");
  } else if (idx_mdot==-1) {
   // check nothing
  } else
   amrex::Error("idx_mdot invalid");

 } // i=0..color_count-1

 int num_sweeps=2;

 if (operation_flag==OP_GATHER_MDOT) {
  // do nothing
 } else if (operation_flag==OP_SCATTER_MDOT) { 
   // (dest,source)
  copy_blobdata(level_blobdata,blobdata);

  for (int i=0;i<color_count;i++) {
   int j=0;
   for (j=0;j<ncomp_mdot_alloc;j++) {
    level_mdot_data[i][j]=mdot_data[i][j];
    level_mdot_comp_data[i][j]=mdot_comp_data[i][j];
   }
   if (idx_mdot>=0) {
    if (j==2*num_interfaces) {
     // do nothing
    } else
     amrex::Error("expecting j==2*num_interfaces");
   } else if (idx_mdot==-1) {
    amrex::Error("expecting idx_mdot>=0");
   } else
    amrex::Error("idx_mdot invalid");
  }

  num_sweeps=1;

 } else
  amrex::Error("operation_flag invalid");

 for (int sweep_num=0;sweep_num<num_sweeps;sweep_num++) {

  for (int ilev = coarsest_level; ilev <= finest_level; ilev++) {

   NavierStokes& ns_level = getLevel(ilev);

   MultiFab* mdot=nullptr;
   MultiFab* mdot_complement=nullptr;
   if (ncomp_mdot==0) {
    mdot=ns_level.localMF[idx_type];
    mdot_complement=ns_level.localMF[idx_type];
   } else if (ncomp_mdot==2*num_interfaces) {
    mdot=ns_level.localMF[idx_mdot];
    mdot_complement=ns_level.localMF[idx_mdot_complement];
   } else {
    mdot=nullptr;
    mdot_complement=nullptr;
    amrex::Error("ncomp_mdot invalid");
   }

   ns_level.ColorSum(
    operation_flag, // OP_GATHER_MDOT or OP_SCATTER_MDOT
    tessellate,  // =1 or 3
    sweep_num,
    ncomp_mdot_alloc,
    ncomp_mdot,
    ns_level.localMF[idx_type],
    ns_level.localMF[idx_color],
    mdot,
    mdot_complement,
    level_blobdata,
    blobdata,  //cum_blobdata
    level_mdot_data,
    level_mdot_comp_data,
    mdot_data,
    mdot_comp_data,
    level_mdot_data_redistribute,
    level_mdot_comp_data_redistribute
    );

   if (operation_flag==OP_GATHER_MDOT) {

    if (sweep_num==0) {

     for (int i=0;i<color_count;i++) {

      sum_blobdata(i,blobdata,level_blobdata,sweep_num);

      int j=0;
      for (j=0;j<ncomp_mdot_alloc;j++) {
       mdot_data[i][j]+=level_mdot_data[i][j];
       mdot_comp_data[i][j]+=level_mdot_comp_data[i][j];
      }
      if (idx_mdot>=0) {
       if (j==2*num_interfaces) {
        // do nothing
       } else
        amrex::Error("expecting j==2*num_interfaces");
      } else if (idx_mdot==-1) {
       // check nothing
      } else
       amrex::Error("idx_mdot invalid");

      if ((level_blobdata[i].im>=1)&&
          (level_blobdata[i].im<=num_materials)) {
       int im_test=level_blobdata[i].im;
       if ((im_test<1)||(im_test>num_materials))
        amrex::Error("im_test invalid");
       blobdata[i].im=im_test;
      } else if (level_blobdata[i].im==0) {
       // do nothing
      } else
       amrex::Error("level_blobdata[i].im invalid");

     }  // i=0..color_count-1

    } else if (sweep_num==1) {

     for (int i=0;i<color_count;i++) {
      sum_blobdata(i,blobdata,level_blobdata,sweep_num);
     }  // i=0..color_count-1

    } else
     amrex::Error("sweep_num invalid");

   } else if (operation_flag==OP_SCATTER_MDOT) { //blobdata not updated.

    if (sweep_num==0) {

     for (int i=0;i<color_count;i++) {

      int j=0;
      for (j=0;j<ncomp_mdot_alloc;j++) {
       mdot_data_redistribute[i][j]+=
	   level_mdot_data_redistribute[i][j];
       mdot_comp_data_redistribute[i][j]+=
	   level_mdot_comp_data_redistribute[i][j];
      }
      if (j==2*num_interfaces) {
       // do nothing
      } else
       amrex::Error("expecting j==2*nsten");

      if ((level_blobdata[i].im>=1)&&
          (level_blobdata[i].im<=num_materials)) {
       if (level_blobdata[i].im==blobdata[i].im) {
        // do nothing
       } else
        amrex::Error("level_blobdata[i].im!=blobdata[i].im");
      } else
       amrex::Error("level_blobdata[i].im invalid");
     }  // i=0..color_count-1

    } else
     amrex::Error("sweep_num invalid for operation_flag==OP_SCATTER_MDOT");

   } else
    amrex::Error("operation_flag invalid");

  } // ilev=coarsest_level..finest_level

  if (operation_flag==OP_GATHER_MDOT) {

   if (sweep_num==0) {

    for (int i=0;i<color_count;i++) {
     Real blobvol=blobdata[i].blob_volume;
     if (blobvol>0.0) {
      for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
       blobdata[i].blob_center_actual[dir]= 
        blobdata[i].blob_center_integral[dir]/blobvol;
      }
      for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
       if (force_blob_symmetry[dir]==1) {
        blobdata[i].blob_center_actual[dir]=0.0;
       } else if (force_blob_symmetry[dir]==0) {
        // do nothing
       } else {
        amrex::Error("force_blob_symmetry[dir] out of range");
       }
      } //dir=0..sdim-2
     } else if (blobvol==0.0) {
      // do nothing
     } else
      amrex::Error("blobvol invalid");
    } // i=0..color_count-1

   } else if (sweep_num==1) {

    for (int i=0;i<color_count;i++) {

     Real blobvol=blobdata[i].blob_volume;
     if (blobvol>0.0) {
      for (int dir=0;dir<6;dir++) {
       blobdata[i].blob_second_moment[dir]= 
        blobdata[i].blob_second_moment[dir]/blobvol;
      }
     } else if (blobvol==0.0) {
      // do nothing
     } else
      amrex::Error("blobvol invalid");

     for (int veltype=0;veltype<3;veltype++) {

      if (blobdata[i].blob_mass_for_velocity[veltype]>0.0) {
       Real** AA3D = new Real*[2*AMREX_SPACEDIM];
       Real** AA2D = new Real*[AMREX_SPACEDIM+1];
       for (int irow=0;irow<2*AMREX_SPACEDIM;irow++) 
        AA3D[irow]=new Real[2*AMREX_SPACEDIM];
       for (int irow=0;irow<AMREX_SPACEDIM+1;irow++) 
        AA2D[irow]=new Real[AMREX_SPACEDIM+1];

       Real BB3D[2*AMREX_SPACEDIM];
       Real BB2D[AMREX_SPACEDIM+1];
       Real XX3D[2*AMREX_SPACEDIM];
       Real XX2D[AMREX_SPACEDIM+1];
       int matrix_ncomp=(2*AMREX_SPACEDIM)*(2*AMREX_SPACEDIM);
       int imatrix=matrix_ncomp*veltype;
       for (int irow=0;irow<2*AMREX_SPACEDIM;irow++) {
        for (int icol=0;icol<2*AMREX_SPACEDIM;icol++) {
         AA3D[irow][icol]=blobdata[i].blob_matrix[imatrix];
         if (irow!=icol)
	  AA3D[irow][icol]=0.0;  // basis functions form an orthogonal set
         if ((irow<AMREX_SPACEDIM+1)&&(icol<AMREX_SPACEDIM+1)) {
          AA2D[irow][icol]=blobdata[i].blob_matrix[imatrix];
	  if (irow!=icol)
	   AA2D[irow][icol]=0.0;  // basis functions form an orthogonal set
	 }
         imatrix++;
        } //icol=0..2*sdim-1
       } //irow=0..2*sdim-1 
       if (imatrix!=matrix_ncomp*(veltype+1))
        amrex::Error("imatrix invalid");

       for (int irow=0;irow<2*AMREX_SPACEDIM;irow++) {
        BB3D[irow]=blobdata[i].blob_RHS[2*AMREX_SPACEDIM*veltype+irow];
        XX3D[irow]=0.0;
	 //rotational motions 
        if ((irow>=AMREX_SPACEDIM)&&(irow<=2*AMREX_SPACEDIM-1)) {
	 if (force_blob_symmetry[0]==1) {
	  BB3D[irow]=0.0;
	 } else if (force_blob_symmetry[0]==0) {
	  // do nothing
	 } else {
          amrex::Error("force_blob_symmetry[0] out of range");
	 }

         if (veltype==1)
	  BB3D[irow]=0.0;

	  //translational motions
	} else if ((irow>=0)&&(irow<AMREX_SPACEDIM)) {
 	 // do nothing
	} else {
	 amrex::Error("irow invalid");
        }

	 //motions for 2D or 3D
        if ((irow>=0)&&(irow<=AMREX_SPACEDIM)) {

         BB2D[irow]=blobdata[i].blob_RHS[2*AMREX_SPACEDIM*veltype+irow];
         XX2D[irow]=0.0;
	  //rotational motions 
         if ((irow>=AMREX_SPACEDIM)&&(irow<=2*AMREX_SPACEDIM-1)) {

	  if (force_blob_symmetry[0]==1) {
           BB2D[irow]=0.0;
          } else if (force_blob_symmetry[0]==0) {
           // do nothing
          } else
           amrex::Error("force_blob_symmetry[0] invalid");

          if (veltype==1)
  	   BB2D[irow]=0.0;

	  //translational motions
	 } else if ((irow>=0)&&(irow<AMREX_SPACEDIM)) {
 	  // do nothing
	 } else {
	  amrex::Error("irow invalid");
         }

	 //3D motions only.
	} else if ((irow>=AMREX_SPACEDIM+1)&&(irow<=2*AMREX_SPACEDIM-1)) {
 	 // do nothing
	} else {
	 amrex::Error("irow invalid");
        } 
       } //for (int irow=0;irow<2*AMREX_SPACEDIM;irow++)

       int mat_status=0;
       if (AMREX_SPACEDIM==3) {
        matrix_solveCPP(AA3D,XX3D,BB3D,mat_status,2*AMREX_SPACEDIM);
       } else if (AMREX_SPACEDIM==2) {
        matrix_solveCPP(AA2D,XX2D,BB2D,mat_status,AMREX_SPACEDIM+1);
       } else
        amrex::Error("dimension bust");

       if (mat_status==1) {
        if (AMREX_SPACEDIM==3) {
         for (int dir=0;dir<2*AMREX_SPACEDIM;dir++) {
          blobdata[i].blob_velocity[2*AMREX_SPACEDIM*veltype+dir]=XX3D[dir];
         }
	 if (force_blob_symmetry[0]==1) {
          for (int dir=AMREX_SPACEDIM;dir<2*AMREX_SPACEDIM;dir++)
           blobdata[i].blob_velocity[2*AMREX_SPACEDIM*veltype+dir]=0.0;
         } else if (force_blob_symmetry[0]==0) {
          // do nothing
         } else {
          amrex::Error("force_blob_symmetry[0] invalid");
         }
        } else if (AMREX_SPACEDIM==2) {
         for (int dir=0;dir<AMREX_SPACEDIM+1;dir++) {
          blobdata[i].blob_velocity[2*AMREX_SPACEDIM*veltype+dir]=XX2D[dir];
         }
	 if (force_blob_symmetry[0]==1) {
          int dir=AMREX_SPACEDIM;
          blobdata[i].blob_velocity[2*AMREX_SPACEDIM*veltype+dir]=0.0;
         } else if (force_blob_symmetry[0]==0) {
          // do nothing
         } else {
          amrex::Error("force_blob_symmetry[0] invalid");
	 }

         if (geom.IsRZ()) {
          int dir=0;
          blobdata[i].blob_velocity[2*AMREX_SPACEDIM*veltype+dir]=0.0;
	 }

        } else
         amrex::Error("dimension bust");

       } else if (mat_status==0) {
        std::cout << "mat_status==0  for i= " << i << '\n';
        amrex::Error("mat_status==0 error");
       } else
        amrex::Error("mat_status invalid");

       for (int irow=0;irow<2*AMREX_SPACEDIM;irow++)
        delete [] AA3D[irow];
       for (int irow=0;irow<AMREX_SPACEDIM+1;irow++)
        delete [] AA2D[irow];

       delete [] AA3D;
       delete [] AA2D;

       if ((veltype==0)||(veltype==2)) {

        for (int irow=0;irow<2*AMREX_SPACEDIM;irow++) {
	 // integral u rho dV
         Real original_mom=blobdata[i].blob_integral_momentum[irow];
         int ibase=veltype*(2*AMREX_SPACEDIM)+irow;
	 // momentum = velocity * mass
	 Real corrected_velocity=blobdata[i].blob_velocity[ibase];
	 Real proposed_velocity=corrected_velocity;
	 Real mass=blobdata[i].blob_integral_momentum[irow+2*AMREX_SPACEDIM];

         Real proposed_mom=proposed_velocity*mass;

	 if ((original_mom*proposed_mom<=0.0)|| 
	     (std::abs(original_mom)<std::abs(proposed_mom))) {
	  if (original_mom*proposed_mom<=0.0) {
	   corrected_velocity=0.0;
	  } else if (std::abs(original_mom)<std::abs(proposed_mom)) {
	   corrected_velocity*=std::abs(original_mom/proposed_mom);
	  } else
	   amrex::Error("original_mom or proposed_mom became corrupt");

	  blobdata[i].blob_velocity[ibase]=corrected_velocity;
	
          if (verbose>=2) {
           if (ParallelDescriptor::IOProcessor()) {
	    std::cout << " ------------------------------------\n";
            std::cout << " avoid momentum overshoot i= " << i << " im= " <<
              blobdata[i].im << '\n';
            std::cout << " irow= " << irow << " veltype= " << veltype << 
             " orig_mom " << original_mom << " proposed_mom " << 
	     proposed_mom << 
	     " corrected_mom " << corrected_velocity*mass << '\n';
	    std::cout << " irow= " << irow << " veltype= " << veltype <<
	     " proposed_velocity= " << proposed_velocity << '\n';
	    std::cout << " irow= " << irow << " veltype= " << veltype <<
	     " blob_velocity= " << blobdata[i].blob_velocity[ibase] << '\n';
	    std::cout << " ------------------------------------\n";
	   } // ParallelDescriptor::IOProcessor()
 	  }  // verbose>=2
         } else if ((original_mom*proposed_mom>0.0)&& 
	            (std::abs(original_mom)>=std::abs(proposed_mom))) {
 	  // do nothing
	 } else
	  amrex::Error("original_mom or proposed_mom invalid");

        } // irow=0...2*sdim-1
        
        Real original_KE=blobdata[i].blob_energy;
        if (original_KE>=0.0) {
         Real proposed_KE=0.0;
         for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
          int ibase=veltype*(2*AMREX_SPACEDIM)+dir;
	  Real mass=blobdata[i].blob_integral_momentum[dir+2*AMREX_SPACEDIM];
	  // (1/2)*mass*(u^2)
          proposed_KE+=
	     blobdata[i].blob_velocity[ibase]*
  	     blobdata[i].blob_velocity[ibase]*mass;
         } // dir=0..sdim-1
         proposed_KE=0.5*proposed_KE;
         if (proposed_KE>original_KE) {
          for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
           int ibase=veltype*(2*AMREX_SPACEDIM)+dir;
  	   blobdata[i].blob_velocity[ibase]*=std::sqrt(original_KE/proposed_KE);
	  } //dir=0..sdim-1
          if (verbose>0) {
           if (ParallelDescriptor::IOProcessor()) {
	    std::cout << " ------------------------------------\n";
            std::cout << " avoid energy overshoot i= " << i << " im= " <<
              blobdata[i].im << '\n';
            std::cout << " veltype= " << veltype << 
             " orig_KE " << original_KE << " proposed_KE " << proposed_KE << 
	     " orig/proposed " << original_KE/proposed_KE;
            for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
             int ibase=veltype*(2*AMREX_SPACEDIM)+dir;
	     std::cout << " dir= " << dir << " veltype= " << veltype <<
  	      " blob_velocity= " << blobdata[i].blob_velocity[ibase] << '\n';
	    } //dir=0..sdim-1
  	    std::cout << " ------------------------------------\n";
	   } //IOProc
	  } //verbose>0 
         } else if ((proposed_KE>=0.0)&&(proposed_KE<=original_KE)) {
	  // do nothing
	 } else
	  amrex::Error("proposed_KE invalid");
        } else
	 amrex::Error("original_KE invalid");

       } else if (veltype==1) {
        // do nothing
       } else
        amrex::Error("veltype invalid");

      } else if (blobdata[i].blob_mass_for_velocity[veltype]==0.0) {
       // do nothing
      } else
       amrex::Error("blobdata[i].blob_mass_for_velocity[veltype] invalid");

     } // veltype=0,1,2

    } // i=0..color_count-1

   } else
    amrex::Error("sweep_num invalid");

  } else if (operation_flag==OP_SCATTER_MDOT) {
   // do nothing
  } else
   amrex::Error("operation_flag invalid");

 } // sweep_num=0..1

 if (verbose>=2) {
  if (ParallelDescriptor::IOProcessor()) {

   std::cout << "in color sum color_count = " << color_count << '\n';
   std::cout << "force_blob_symmetry[0]= " << force_blob_symmetry[0]  << '\n';

   for (int i=0;i<color_count;i++) {

    Real blobvol=blobdata[i].blob_volume;

    if (blobvol==0.0) {
     std::cout << "NULL COLOR i,im " << i << ' ' << 
       blobdata[i].im << '\n';
    } else if (blobvol>0.0) {
     int imbase=blobdata[i].im;
     if ((imbase<1)||(imbase>num_materials))
      amrex::Error("imbase invalid");

     std::cout << "i, vol,perim,im " << i << ' ' <<
      blobdata[i].blob_volume << ' ' <<
      blobdata[i].blob_perim << ' ' <<
      imbase << '\n';

     for (int imnbr=0;imnbr<num_materials;imnbr++) {
      std::cout << "perim(nbr): i,imnbr,im,perim " << i << ' ' <<
      imnbr+1 << ' ' <<
      imbase << ' ' <<
      blobdata[i].blob_perim_mat[imnbr] << '\n';
     }
     for (int im1=0;im1<num_materials;im1++) {
      for (int im2=0;im2<num_materials;im2++) {
       std::cout << "blob_triple_perim: i,im1,im2,im,perim " << i << ' ' <<
       im1+1 << ' ' << im2+1 << ' ' <<
       imbase << ' ' <<
       blobdata[i].blob_triple_perim[im1][im2] << '\n';
      } // im2=0..num_materials-1
     } // im1=0..num_materials-1

     if (blobvol>0.0) {
      std::cout << "i,x,y,z " << i << ' ' <<
       blobdata[i].blob_center_actual[0] << ' ' <<
       blobdata[i].blob_center_actual[1] << ' ' <<
       blobdata[i].blob_center_actual[AMREX_SPACEDIM-1] << '\n';

      std::cout << "i,xx,xy,xz,yy,yz,zz " << i << ' ' <<
       blobdata[i].blob_second_moment[0] << ' ' <<
       blobdata[i].blob_second_moment[1] << ' ' <<
       blobdata[i].blob_second_moment[2] << ' ' <<
       blobdata[i].blob_second_moment[3] << ' ' <<
       blobdata[i].blob_second_moment[4] << ' ' <<
       blobdata[i].blob_second_moment[5] << '\n';

     } else if (blobvol==0.0) {
      // do nothing
     } else
      amrex::Error("blobvol invalid");

     for (int veltype=0;veltype<3;veltype++) {
      std::cout << " im= " << imbase <<
       " veltype= " << veltype << " blob_mass_for_velocity= " <<
       blobdata[i].blob_mass_for_velocity[veltype] << '\n';
      for (int dir=0;dir<2*AMREX_SPACEDIM;dir++) {
       std::cout << " im= " << imbase <<
        " veltype= " << veltype << " dir= " << dir << " vel= " <<
        blobdata[i].blob_velocity[2*AMREX_SPACEDIM*veltype+dir] << '\n';
      } // dir=0..2 sdim-1
     } // veltype=0..2
     for (int dir=0;dir<2*AMREX_SPACEDIM;dir++) {

      std::string momstring="momentum";
      if (dir==0) 
       momstring="momx";
      else if (dir==1)
       momstring="momy";
      else if ((dir==2)&&(AMREX_SPACEDIM==3))
       momstring="momz";
      else if (dir==AMREX_SPACEDIM)
       momstring="momxy";
      else if (dir==AMREX_SPACEDIM+1)
       momstring="momxz";
      else if (dir==AMREX_SPACEDIM+2)
       momstring="momyz";
      else
       amrex::Error("dir invalid");

      std::string velstring="vel";
      if (dir==0) 
       velstring="velx";
      else if (dir==1)
       velstring="vely";
      else if ((dir==2)&&(AMREX_SPACEDIM==3))
       velstring="velz";
      else if (dir==AMREX_SPACEDIM)
       velstring="velxy";
      else if (dir==AMREX_SPACEDIM+1)
       velstring="velxz";
      else if (dir==AMREX_SPACEDIM+2)
       velstring="velyz";
      else
       amrex::Error("dir invalid");

      std::cout << " im= " << imbase <<
       " dir= " << dir << " " << momstring << "= " <<
       blobdata[i].blob_integral_momentum[dir] << '\n';
      Real numerator=blobdata[i].blob_integral_momentum[dir];
      Real denom=blobdata[i].blob_integral_momentum[2*AMREX_SPACEDIM+dir];
      Real avg_vel=numerator;
      if (denom>0.0) 
       avg_vel/=denom;
      std::cout << " im= " << imbase <<
       " dir= " << dir << " average " << velstring << "= " << avg_vel << '\n';
      std::cout << " im= " << imbase <<
       " dir= " << dir << " " << momstring << " divisor= " << denom << '\n';
     } // dir=0..2 sdim-1
     std::cout << " im= " << imbase <<
      " energy= " <<
      blobdata[i].blob_energy << '\n';

    } else 
     amrex::Error("blobdata[i].blob_volume invalid");
   } // i=0..color_count-1
  } // if IOproc
 } else if ((verbose==0)||(verbose==1)) {
  // do nothing
 } else
  amrex::Error("verbose must be >=0");

} // end subroutine ColorSumALL


void
NavierStokes::Type_level(
  MultiFab* typemf,Vector<int>& type_flag,
  int zero_diag_flag) {

 std::string local_caller_string="Type_level";

 int finest_level=parent->finestLevel();

 int ncomp_type=num_materials;
 if (zero_diag_flag==1) {
  ncomp_type=2;
 } else if (zero_diag_flag==0) {
  ncomp_type=num_materials;
 } else
  amrex::Error("zero_diag_flag invalid");

 int typedim=type_flag.size();
 if (typedim!=ncomp_type)
  amrex::Error("typedim invalid");

 Vector< Vector<int> > local_type_flag;
 local_type_flag.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  local_type_flag[tid].resize(typedim);
  for (int i=0;i<typedim;i++)
   local_type_flag[tid][i]=type_flag[i];
 }  // tid

 if (level>finest_level)
  amrex::Error("level invalid Type_level");

 MultiFab* Type_Source_MF=nullptr;
 int ncomp_source=0;

 if (zero_diag_flag==0) {
  Type_Source_MF=getStateDist(1,cur_time_slab,local_caller_string);
  ncomp_source=Type_Source_MF->nComp();
  if (Type_Source_MF->nComp()!=num_materials*(1+AMREX_SPACEDIM))
   amrex::Error("Type_Source_MF->nComp() invalid");
 } else if (zero_diag_flag==1) {
  Type_Source_MF=localMF[ONES_GROW_MF];
  ncomp_source=Type_Source_MF->nComp();
  if (Type_Source_MF->nComp()!=1)
   amrex::Error("Type_Source_MF->nComp() invalid");
 } else
  amrex::Error("zero_diag_flag invalid");

 if (typemf->nGrow()!=1)
  amrex::Error("typemf->nGrow()!=1");

 const Real* dx = geom.CellSize();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(typemf->boxArray().d_numPts());

// COLORING LOOP
#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*typemf,false); mfi.isValid(); ++mfi) {
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

  FArrayBox& source_fab=(*Type_Source_MF)[mfi];
  FArrayBox& typefab=(*typemf)[mfi];

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // updates one ghost cell.
   // declared in: LEVELSET_3D.F90
   //  for each cell,
   //   if is_rigid(im)==1 and LS>=0 then type=im
  fort_gettypefab(
   source_fab.dataPtr(),
   ARLIM(source_fab.loVect()),ARLIM(source_fab.hiVect()),
   typefab.dataPtr(),ARLIM(typefab.loVect()),ARLIM(typefab.hiVect()),
   xlo,dx,
   tilelo,tilehi,
   fablo,fabhi,&bfact,
   local_type_flag[tid_current].dataPtr(),
   &ncomp_type,
   &ncomp_source,
   &zero_diag_flag);
 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_GETTYPEFAB,"Type_level");

 for (int tid=1;tid<thread_class::nthreads;tid++) {
  for (int im=0;im<ncomp_type;im++) { 
   if (local_type_flag[tid][im]>local_type_flag[0][im])
    local_type_flag[0][im]=local_type_flag[tid][im];
  } // im
 } // tid

 ParallelDescriptor::Barrier();

 for (int im=0;im<ncomp_type;im++) {
  ParallelDescriptor::ReduceIntMax(local_type_flag[0][im]);
  type_flag[im]=local_type_flag[0][im];
 }

 if (zero_diag_flag==0) {
  delete Type_Source_MF;
 } else if (zero_diag_flag==1) {
  // do nothing
 } else
  amrex::Error("zero_diag_flag invalid");
 
}  // subroutine Type_level

// zero_diag_flag=0 => color by material
// zero_diag_flag=1 => color by masked cells
void NavierStokes::TypeALL(int idx_type,Vector<int>& type_flag,
		int zero_diag_flag) {

 int finest_level=parent->finestLevel();

 int ncomp_type=num_materials;
 if (zero_diag_flag==1) {
  ncomp_type=2;
 } else if (zero_diag_flag==0) {
  ncomp_type=num_materials;
 } else
  amrex::Error("zero_diag_flag invalid");

 type_flag.resize(ncomp_type);
 for (int im=0;im<ncomp_type;im++) {
  type_flag[im]=0;
 }
 allocate_array(1,1,-1,idx_type);
 if (level!=0)
  amrex::Error("level=0 in TypeALL");

  // updates one ghost cell.
 for (int k = 0; k <= finest_level; k++) {
  NavierStokes& ns_level = getLevel(k);
  ns_level.Type_level(ns_level.localMF[idx_type],type_flag,zero_diag_flag);
 }
 int color_counter=0;
 for (int im=0;im<ncomp_type;im++) {
  color_counter+=type_flag[im];
 }
 if (verbose>0)
  if (ParallelDescriptor::IOProcessor()) {
   for (int im=0;im<ncomp_type;im++) 
    std::cout << "TypeALL im,type " << im << ' ' <<
     type_flag[im] << ' ' << '\n';
   std::cout << "TypeALL color_counter= " << color_counter << '\n'; 
  }
} // subroutine TypeALL

void NavierStokes::remove_pressure_work_vars() {

 delete_localMF(UMACSTAR_MF,AMREX_SPACEDIM);
 delete_localMF(RESTART_UMACSTAR_MF,AMREX_SPACEDIM);
 delete_localMF(GRADPEDGE_MF,AMREX_SPACEDIM);
 delete_localMF(PEDGE_MF,AMREX_SPACEDIM);
 delete_localMF(AMRSYNC_PRES_MF,AMREX_SPACEDIM);

} // end subroutine remove_pressure_work_vars

// called from: NavierStokes::multiphase_project
void NavierStokes::remove_project_variables() {

 delete_localMF(POLDHOLD_MF,1);
 delete_localMF(ONES_MF,1);
 delete_localMF(ONES_GROW_MF,1);
 delete_localMF(TYPE_ONES_MF,1);
 delete_localMF(COLOR_ONES_MF,1);
 delete_localMF(OUTER_ITER_PRESSURE_MF,1);
}

void NavierStokes::allocate_MAC_velocityALL(int nsolve,int idx) {

 if (level!=0)
  amrex::Error("level invalid allocate_MAC_velocity_ALL");

 if ((nsolve!=1)&&(nsolve!=AMREX_SPACEDIM))
  amrex::Error("nsolve invalid");

 int finest_level=parent->finestLevel();
 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   ns_level.new_localMF(idx+dir,nsolve,0,dir);
   ns_level.setVal_localMF(idx+dir,0.0,0,nsolve,0);
  }
 } // ilev

} // allocate_MAC_velocityALL


void NavierStokes::remove_MAC_velocityALL(int idx) {

 if (level!=0)
  amrex::Error("level invalid remove_MAC_velocityALL");

 int finest_level=parent->finestLevel();
 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.delete_localMF(idx,AMREX_SPACEDIM);
 }
} // end subroutine remove_MAC_velocityALL

// called from:
// NavierStokes::update_SEM_forcesALL
// NavierStokes::multiphase_project
// NavierStokes::diffusion_heatingALL 
void NavierStokes::remove_FACE_WEIGHT_vars() {

 delete_localMF(FACE_WEIGHT_MF,AMREX_SPACEDIM);
 delete_localMF(OFF_DIAG_CHECK_MF,1);

} // end subroutine remove_FACE_WEIGHT_vars()

// called from:
// NavierStokes::update_SEM_forcesALL
// NavierStokes::multiphase_project
// NavierStokes::diffusion_heatingALL 
void NavierStokes::allocate_FACE_WEIGHT(
 int nsolve,
 int project_option,
 int face_weight_op) {
 
 std::string local_caller_string="allocate_FACE_WEIGHT";

 int finest_level=parent->finestLevel();

 bool use_tiling=ns_tiling;

 if (project_option_momeqn(project_option)==1) {
  //do nothing
 } else if (project_option_momeqn(project_option)==0) {
  //do nothing
 } else
  amrex::Error("project_option_momeqn(project_option) invalid30");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");
 Vector<int> scomp;
 Vector<int> ncomp;
 int state_index;
 int ncomp_check;

  //num_materials_combine=1
 get_mm_scomp_solver(
  1,
  project_option,
  state_index,
  scomp,
  ncomp,
  ncomp_check);

 if (ncomp_check!=nsolve)
  amrex::Error("ncomp_check alid");

 resize_mask_nbr(1);
 debug_ngrow(MASK_NBR_MF,1,local_caller_string);

 int local_face_index=FACECOMP_FACEDEN;  // 1/rho
 int local_face_ncomp=FACECOMP_NCOMP;

 if (project_option_is_valid(project_option)==1) {
  // do nothing
 } else
  amrex::Error("project_option invalid 7875");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  debug_ngrow(FACE_VAR_MF+dir,0,local_caller_string);
 }

 debug_ngrow(CELL_VISC_MF,1,local_caller_string);
 debug_ngrow(CELL_DEN_MF,1,local_caller_string);
 if (localMF[CELL_VISC_MF]->nComp()!=1)
  amrex::Error("localMF[CELL_VISC_MF]->nComp() invalid");
 if (localMF[CELL_DEN_MF]->nComp()!=1)
  amrex::Error("localMF[CELL_DEN_MF]->nComp() invalid");

 int bcsize=AMREX_SPACEDIM*2*nsolve*grids.size();
 bcpres_array.resize(bcsize);
 for (int gridno=0;gridno<grids.size();gridno++) {

   // presbc declared as presbc(sdim,2) in fortran
   // components ordered at  1,1  2,1  3,1  1,2  2,2  3,2
  Vector<int> presbc;
  getBCArray_list(presbc,state_index,gridno,scomp,ncomp);
  if (presbc.size()!=AMREX_SPACEDIM*2*nsolve)
   amrex::Error("presbc.size() invalid");

  int ibase=AMREX_SPACEDIM*2*nsolve*gridno;
  int j=0;
  for (int nn=0;nn<nsolve;nn++) {
   for (int side=0;side<=1;side++) {
    for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
     int pbc=presbc[j];
     bcpres_array[ibase+j]=pbc;

     if (project_option_singular_possible(project_option)==1) {
      // do nothing
     } else if (project_option_singular_possible(project_option)==0) {
      // do nothing
     } else
      amrex::Error("project_option_singular_possible bad allocate_FACE_WEIGHT");

     j++;
    }  // dir
   } // side
  } // nn
 } // gridno

 int GFM_flag=0;
 int adjust_temperature=-1;  // adjust FACECOMP_FACEHEAT

  // adjust FACECOMP_FACEHEAT if thermal diffusion
  // and phase change using sharp interface method.
 if (project_option==SOLVETYPE_HEAT) { 

  for (int im=0;im<2*num_interfaces;im++) {
   Real LL=get_user_latent_heat(im+1,293.0,1);
   if (LL!=0.0) {
    if (is_GFM_freezing_model(freezing_model[im])==1) {
     GFM_flag=1; 
    } else if (is_GFM_freezing_model(freezing_model[im])==0) {
     // do nothing
    } else 
     amrex::Error("is_GFM_freezing_model invalid");
   } else if (LL==0.0) {
    // do nothing
   } else
    amrex::Error("latent_heat[im] invalid");
  } // im=0..2 num_interfaces -1

 } else if ((project_option>=SOLVETYPE_SPEC)&&
            (project_option<SOLVETYPE_SPEC+num_species_var)) {

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
     amrex::Error("is_GFM_freezing_model invalid");
   } else if (LL==0.0) {
    // do nothing
   } else
    amrex::Error("latent_heat[im] (LL) invalid");

  } // im=0..2 num_interfaces -1
 } else if (project_option_is_valid(project_option)==1) {
  // do nothing
 } else
  amrex::Error("project_option invalid31");

 if (GFM_flag==1) {
  if (face_weight_op==SUB_OP_FOR_MAIN) {
   stefan_solver_init(
    localMF[CELL_DEN_MF],
    adjust_temperature,
    project_option);
  } else if (face_weight_op==SUB_OP_FOR_SDC) {
   //do nothing
  } else
   amrex::Error("face_weight_op invalid");
 }

 const Real* dx = geom.CellSize();

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  new_localMF(FACE_WEIGHT_MF+dir,nsolve,0,dir);
 } // dir
 new_localMF(OFF_DIAG_CHECK_MF,nsolve,0,-1);

 if (project_option_projection(project_option)==1) {
  local_face_index=FACECOMP_FACEDEN;  // 1/rho_added
 } else if (project_option==SOLVETYPE_PRESEXTRAP) { 
   // 1/rho (only used in eval_face_coeff for sanity check purposes)
  local_face_index=FACECOMP_FACEDEN;  
 } else if (project_option==SOLVETYPE_HEAT) {
  local_face_index=FACECOMP_FACEHEAT; 
 } else if (project_option==SOLVETYPE_VISC) {
  local_face_index=FACECOMP_FACEVISC; 
 } else if ((project_option>=SOLVETYPE_SPEC)&&
            (project_option<SOLVETYPE_SPEC+num_species_var)) { // rho D
  local_face_index=FACECOMP_FACESPEC+project_option-SOLVETYPE_SPEC;
 } else
  amrex::Error("project_option invalid allocate_FACE_WEIGHT");

 setVal_localMF(OFF_DIAG_CHECK_MF,0.0,0,nsolve,0);

 for (int facewt_iter=0;facewt_iter<=1;facewt_iter++) {

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(
    localMF[CELL_DEN_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(*localMF[CELL_DEN_MF],use_tiling); 
       mfi.isValid(); ++mfi) {

   BL_ASSERT(grids[mfi.index()] == mfi.validbox());

   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();
   int bfact=parent->Space_blockingFactor(level);

   FArrayBox& cenden=(*localMF[CELL_DEN_MF])[mfi];  // 1/rho
   FArrayBox& cenvisc=(*localMF[CELL_VISC_MF])[mfi];

   FArrayBox& xfwt=(*localMF[FACE_WEIGHT_MF])[mfi];
   FArrayBox& yfwt=(*localMF[FACE_WEIGHT_MF+1])[mfi];
   FArrayBox& zfwt=(*localMF[FACE_WEIGHT_MF+AMREX_SPACEDIM-1])[mfi];

   FArrayBox& offdiagcheck=(*localMF[OFF_DIAG_CHECK_MF])[mfi];
 
   FArrayBox& xface=(*localMF[FACE_VAR_MF])[mfi];  
   FArrayBox& yface=(*localMF[FACE_VAR_MF+1])[mfi];  
   FArrayBox& zface=(*localMF[FACE_VAR_MF+AMREX_SPACEDIM-1])[mfi];  

   FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi];  // mask=1 at fine-fine bc
   const Real* xlo = grid_loc[gridno].lo();

   Vector<int> presbc;
   getBCArray_list(presbc,state_index,gridno,scomp,ncomp);
   if (presbc.size()==2*AMREX_SPACEDIM*nsolve) {
    // do nothing
   } else
    amrex::Error("presbc.size() invalid");

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // fort_buildfacewt is declared in LEVELSET_3D.F90
   fort_buildfacewt(
    im_elastic_map.dataPtr(),
    &num_FSI_outer_sweeps,
    &FSI_outer_sweeps,
    &facewt_iter,
    &level,
    &finest_level,
    &nsolve,
    &local_face_index,
    &local_face_ncomp,
    xlo,
    dx,
    offdiagcheck.dataPtr(),
    ARLIM(offdiagcheck.loVect()),ARLIM(offdiagcheck.hiVect()),
    cenden.dataPtr(),
    ARLIM(cenden.loVect()),ARLIM(cenden.hiVect()),
    cenvisc.dataPtr(),
    ARLIM(cenvisc.loVect()),ARLIM(cenvisc.hiVect()),
    xfwt.dataPtr(),ARLIM(xfwt.loVect()),ARLIM(xfwt.hiVect()),
    yfwt.dataPtr(),ARLIM(yfwt.loVect()),ARLIM(yfwt.hiVect()),
    zfwt.dataPtr(),ARLIM(zfwt.loVect()),ARLIM(zfwt.hiVect()),
    xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()),
    yface.dataPtr(),ARLIM(yface.loVect()),ARLIM(yface.hiVect()),
    zface.dataPtr(),ARLIM(zface.loVect()),ARLIM(zface.hiVect()),
    maskfab.dataPtr(),ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    min_face_wt[tid_current].dataPtr(),
    max_face_wt[tid_current].dataPtr(),
    presbc.dataPtr(),
    &visc_coef,
    &uncoupled_viscosity,
    &project_option);

  }  // mfi
} // omp
  ns_reconcile_d_num(LOOP_BUILDFACEWT,"allocate_FACE_WEIGHT");

 } // facewt_iter=0..1

 for (int tid=1;tid<thread_class::nthreads;tid++) {

  for (int iwt=0;iwt<NCOMP_FACE_WT;iwt++) {
   if (min_face_wt[tid][iwt]<min_face_wt[0][iwt])
    min_face_wt[0][iwt]=min_face_wt[tid][iwt];
   if (max_face_wt[tid][iwt]>max_face_wt[0][iwt])
    max_face_wt[0][iwt]=max_face_wt[tid][iwt];
  }

 } // tid

 ParallelDescriptor::Barrier();


 for (int iwt=0;iwt<NCOMP_FACE_WT;iwt++) {
  ParallelDescriptor::ReduceRealMin(min_face_wt[0][iwt]);
  ParallelDescriptor::ReduceRealMax(max_face_wt[0][iwt]);
 }
 for (int tid=1;tid<thread_class::nthreads;tid++) {
  for (int iwt=0;iwt<NCOMP_FACE_WT;iwt++) {
   min_face_wt[tid][iwt]=min_face_wt[0][iwt];
   max_face_wt[tid][iwt]=max_face_wt[0][iwt];
  }
 }

} // end subroutine allocate_FACE_WEIGHT 

void NavierStokes::sanity_check_face_wt(int project_option) {

 for (int tid=0;tid<thread_class::nthreads;tid++) {
  for (int iwt=0;iwt<NCOMP_FACE_WT;iwt++) {

   if ((min_face_wt[tid][iwt]==min_face_wt[0][iwt])&&
       (max_face_wt[tid][iwt]==max_face_wt[0][iwt])&&
       (min_face_wt[tid][iwt]<=max_face_wt[tid][iwt])&&
       (min_face_wt[tid][iwt]>=0.0)) {
    // do nothing
   } else
    amrex::Error("min_face_wt or max_face_wt invalid");

  } //iwt=0..NCOMP_FACE_WT-1
    
  if (max_face_wt[tid][CC_COMP_FACE_WT]<=1.0) {
   //do nothing
  } else
   amrex::Error("(max_face_wt[tid][CC_COMP_FACE_WT]<=1.0) violated");

  if (max_face_wt[tid][MERGE_COMP_FACE_WT]<=
      2.0*max_face_wt[tid][DD_COMP_FACE_WT]) {
   //do nothing
  } else
   amrex::Error("expecting max: MERGE_COMP<=2.0*DD_COMP");

 } //tid=0..nthreads-1

 if (project_option_singular_possible(project_option)==1) {

  if ((max_face_wt[0][DD_COMP_FACE_WT]>0.0)&&
      (max_face_wt[0][MERGE_COMP_FACE_WT]>0.0)&&
      (max_face_wt[0][MERGE_COMP_FACE_WT]<=
       2.0*max_face_wt[0][DD_COMP_FACE_WT])) {

   if (mglib_max_ratio>1.0) {

    min_interior_coeff=max_face_wt[0][MERGE_COMP_FACE_WT]/mglib_max_ratio;

   } else
    amrex::Error("mglib_max_ratio invalid");

  } else {
   std::cout << "max_face_wt[0][DD_COMP_FACE_WT] " <<
    max_face_wt[0][DD_COMP_FACE_WT] << '\n';
   std::cout << "max_face_wt[0][MERGE_COMP_FACE_WT] " <<
    max_face_wt[0][MERGE_COMP_FACE_WT] << '\n';
   for (int iwt=0;iwt<NCOMP_FACE_WT;iwt++) {
    std::cout << "iwt,min_face_wt " << iwt << ' ' << 
     min_face_wt[0][iwt] << '\n';
    std::cout << "iwt,max_face_wt " << iwt << ' ' << 
     max_face_wt[0][iwt] << '\n';
   }

   amrex::Error("max_face_wt invalid");
  }

 } else if (project_option_singular_possible(project_option)==0) {

  min_interior_coeff=0.0;

 } else
  amrex::Error("project_option_singular_possible invalid");

} // end subroutine sanity_check_face_wt
  
// called from NavierStokes::multiphase_project
void NavierStokes::allocate_project_variables(int nsolve,int project_option) {
 
 std::string local_caller_string="allocate_project_variables";

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 if (project_option_momeqn(project_option)==1) {
  // do nothing
 } else if (project_option_momeqn(project_option)==0) {
  // do nothing
 } else
  amrex::Error("project_option_momeqn invalid32");

 debug_ngrow(FACE_VAR_MF,0,local_caller_string);

 int state_index;
 Vector<int> scomp;
 Vector<int> ncomp;
 int ncomp_check;
  //num_materials_combine=1
 get_mm_scomp_solver(
  1,
  project_option,
  state_index,
  scomp,
  ncomp,
  ncomp_check);
 if (ncomp_check!=nsolve)
  amrex::Error("nsolve or ncomp_check invalid 3919");
 
 MultiFab& S_new=get_new_data(state_index,slab_step+1);

  // in: allocate_project_variables
  // ONES_MF=1 if diag>0  ONES_MF=0 if diag==0.
 new_localMF(ONES_MF,1,0,-1);
 new_localMF(ONES_GROW_MF,1,1,-1);
 setVal_localMF(ONES_MF,1.0,0,1,0);
 setVal_localMF(ONES_GROW_MF,1.0,0,1,1);

 new_localMF(POLDHOLD_MF,nsolve,0,-1);
 setVal_localMF(POLDHOLD_MF,0.0,0,nsolve,0);

  // S_new and outer_iter_pressure will both hold S^*
  // at the very beginning.
  // when solving, -alpha S + div beta grad S = -alpha S^* ,
  // initially S_init=S^* for the pressure projection problem and
  // temperature solve.
  // For the viscous solve, it might be
  // in the future that S_init<>S^*.  In this case, something else
  // must be copied into ``initial_guess'' below.
  // At the very end of this routine:
  //  OUTER_ITER_PRESSURE = S_init
  //  S_new = S_init
  //  POLDHOLD = S^* - S^init

   // this holds S^*
 new_localMF(OUTER_ITER_PRESSURE_MF,nsolve,0,-1);

 int adjust_temperature=1; 
 int GFM_flag=0;

  // temperature diffusion
 if (project_option==SOLVETYPE_HEAT) {

    // MEHDI VAHAB HEAT SOURCE
    // T^new=T^* += dt A Q/(rho cv V) 
    // in: allocate_project_variables
    // NavierStokes::heat_source_term_flux_source  (in:NavierStokes.cpp)
    // heat_source_term_flux_source calls GODUNOV_3D::fort_heatsource_face
  heat_source_term_flux_source();

  if (is_phasechange==1) { //get_user_latent_heat!=0.0?
    // both S_new and outer_iter_pressure are adjusted.
    // vol*(T-T^n)*(rho cv)/dt-vol*grad dot k grad T = -1/dt vol*div u+
    //   diffusionRHS
    // c1 (T-Tn)+c2(T-TSAT)-LT=F
    // (c1+c2)T - (c1 Tn+c2 TSAT) -LT=F
    // (c1+c2)TN=c1 Tn + c2 TSAT  TN=(c1 Tn + c2 TSAT)/(c1+c2)
    // Snew=outer_iter_pressure=TN
    // c1=rho cv/dt   c2=(1/vol) sum_face Aface k_m/(theta dx)
    //
    //  For the finite volume method:
    //   a) find V=[k grad T]/L
    //   b) unsplit advection to find F^*
    //   c) T_ice=(F T_ice + (F^*-F) Tsat)/F^*
    //      T_water does not change
   for (int im=0;im<2*num_interfaces;im++) {
    Real LL=get_user_latent_heat(im+1,293.0,1);
    if (LL!=0.0) {
     if (is_GFM_freezing_model(freezing_model[im])==1) {
      GFM_flag=1;
     } else if (is_GFM_freezing_model(freezing_model[im])==0) {
      // do nothing
     } else 
      amrex::Error("is_GFM_freezing_model invalid");
    } else if (LL==0.0) {
     // do nothing
    } else
     amrex::Error("latent_heat[im] (LL) invalid");
   } // im=0..2 num_interfaces -1
  } else if (is_phasechange==0) {
   // do nothing
  } else
   amrex::Error("is_phasechange invalid");

 } // project_option==SOLVETYPE_HEAT

 if ((project_option>=SOLVETYPE_SPEC)&&
     (project_option<SOLVETYPE_SPEC+num_species_var)) {

  if (is_phasechange==1) { //get_user_latent_heat!=0.0?

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
      amrex::Error("is_GFM_freezing_model invalid");
    } else if (LL==0.0) {
     // do nothing
    } else
     amrex::Error("latent_heat[im] (LL) invalid");
     
   } // im=0.. 2 num_interfaces -1

  } else if (is_phasechange==0) {
   // do nothing
  } else
   amrex::Error("is_phasechange invalid");
 }

  // both S_new and OUTER_ITER_PRESSURE_MF are modified when 
  // adjust_temperature==1
 if (GFM_flag==1) {
  stefan_solver_init(
   localMF[OUTER_ITER_PRESSURE_MF],
   adjust_temperature,
   project_option);
 }

 MultiFab* current_contents_mf=nullptr;

  // this is S^*
  // alpha(S - S^*) - div beta grad S = 0
 if (state_index==State_Type) {
  current_contents_mf=getState_list(1,scomp,ncomp,cur_time_slab);
 } else {
  current_contents_mf=nullptr;
  amrex::Error("state_index invalid");
 }

  // ``OUTER_ITER_PRESSURE'' = S_new = S^*
 MultiFab::Copy(*localMF[OUTER_ITER_PRESSURE_MF],
		*current_contents_mf,0,0,nsolve,0);

  //project_option_olddata_needed is declared in: NavierStokes.cpp
 if (project_option_olddata_needed(project_option)==1) { 
  
   // if S^initial <> S^*, then put something else into
   // ``initial_guess.''   
   // outer_iter_pressure=S_new=S^*
   // For now, initial_guess=S^*
  MultiFab* initial_guess;
  initial_guess=localMF[OUTER_ITER_PRESSURE_MF]; 
  
  MultiFab* dp=new MultiFab(grids,dmap,nsolve,0,
   MFInfo().SetTag("dp"),FArrayBoxFactory());

  MultiFab::Copy(*dp,*initial_guess,0,0,nsolve,0);
   // dp=initial_guess - S^*   (S_new=S^*)
  MultiFab::Subtract(*dp,*current_contents_mf,0,0,nsolve,0);
   // snew+=(initial_guess - S^*)=initial_guess  (S_new=S^* beforehand)
  MultiFab::Add(*current_contents_mf,*dp,0,0,nsolve,0);

  int scomp_temp=0;
  for (int ilist=0;ilist<scomp.size();ilist++) {
   MultiFab::Copy(S_new,*current_contents_mf,
		  scomp_temp,scomp[ilist],ncomp[ilist],0);
   scomp_temp+=ncomp[ilist];
  }
  if (scomp_temp!=nsolve)
   amrex::Error("scomp_temp invalid");

   // dp=S^init-S^*   dS=S-S^init
   // alpha(S^init + dS - S^*) - div beta grad (S^init + dS) = 0
   // alpha dS - div beta grad dS = -alpha dp + div beta grad S^init 
   //                             = alpha POLDHOLD + div beta grad S^init
   // OUTER_ITER_PRESSURE=S^* + (S^init - S^*)=S^init
  MultiFab::Add(*localMF[OUTER_ITER_PRESSURE_MF],*dp,0,0,nsolve,0);
   // POLDHOLD=0 - (S^init-S^*) = S^* - S^init
  MultiFab::Subtract(*localMF[POLDHOLD_MF],*dp,0,0,nsolve,0);
   // later on, (1) UMAC_MF-=beta grad S^init,  (S_new=S^init)
   //           (2) S_new=0.0
  delete dp;

 } else if (project_option_olddata_needed(project_option)==0) { 
  // do nothing
 } else
  amrex::Error("project_option_olddata_needed invalid33");

 delete current_contents_mf;

} // end subroutine allocate_project_variables


void NavierStokes::allocate_pressure_work_vars(int nsolve,int project_option) {

 if (project_option_momeqn(project_option)==1) {
  // do nothing
 } else if (project_option_momeqn(project_option)==0) {
  // do nothing
 } else
  amrex::Error("project_option invalid34");

 if ((nsolve!=1)&&(nsolve!=AMREX_SPACEDIM))
  amrex::Error("nsolve invalid");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  new_localMF(UMACSTAR_MF+dir,nsolve,0,dir);
  new_localMF(RESTART_UMACSTAR_MF+dir,nsolve,0,dir);
  new_localMF(GRADPEDGE_MF+dir,nsolve,0,dir);

   // PEDGE_MF only used if pressure projection.
   // 0=use_face_pres=VALID_PEDGE  
   // 1=(2nd component) pface=PRESSURE_PEDGE
  new_localMF(PEDGE_MF+dir,NCOMP_PEDGE,0,dir);

  new_localMF(AMRSYNC_PRES_MF+dir,nsolve,0,dir);
  setVal_localMF(AMRSYNC_PRES_MF+dir,1.0e+30,0,nsolve,0);
 } // dir=0..sdim-1

} // subroutine allocate_pressure_work_vars

// called in the setup stage from NavierStokes::multiphase_project when
// project_option==SOLVETYPE_PRES 
void NavierStokes::overwrite_outflow() {
 
 bool use_tiling=ns_tiling;

 const Real* dx = geom.CellSize();
 const Real* prob_lo   = geom.ProbLo();
 const Real* prob_hi   = geom.ProbHi();
 MultiFab& U_new = get_new_data(State_Type,slab_step+1); 

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   MultiFab& Umac_new=get_new_data(Umac_Type+dir,slab_step+1);
   if (Umac_new.nComp()!=1)
    amrex::Error("Umac_new.nComp() invalid");

   if (thread_class::nthreads<1)
    amrex::Error("thread_class::nthreads invalid");
   thread_class::init_d_numPts(U_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
   for (MFIter mfi(U_new,use_tiling); mfi.isValid(); ++mfi) {
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
    FArrayBox& vel=U_new[mfi];
    FArrayBox& velmac=Umac_new[mfi];
 
    Vector<int> presbc=getBCArray(State_Type,gridno,STATECOMP_PRES,1);

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     // fort_forcevelocity is declared in: PROB.F90
    fort_forcevelocity(
      prob_lo,prob_hi,
      vel.dataPtr(),ARLIM(vel.loVect()),ARLIM(vel.hiVect()),
      velmac.dataPtr(),ARLIM(velmac.loVect()),ARLIM(velmac.hiVect()),
      &dir,xlo,dx,
      tilelo,tilehi,
      fablo,fabhi,&bfact,
      &cur_time_slab,
      presbc.dataPtr(),
      outflow_velocity_buffer_size.dataPtr());
   } // mfi
}  // omp
   ns_reconcile_d_num(LOOP_FORCEVELOCITY,"overwrite_outlow");

 } // dir=0..sdim-1

} // end subroutine overwrite_outflow


// called from: mac_update (from updatevelALL from multiphase_project)
//    residual_correction_form (from multiphase_project)
//    relaxLEVEL (from mg_cycleALL from multiphase_preconditioner from
//      multiphase_SHELL_preconditioner from multiphase_project
//      or from multphase_project)
//
// macdest=macsrc+gp
// for pressure projection: 
// macdest=UMAC,UMACSTAR,UMAC
// macsrc =UMAC,UMACSTAR,MAC_TEMP
void NavierStokes::correct_velocity(
  int project_option,
  int macdest,
  int macsrc,
  int gp,int nsolve) {

 std::string local_caller_string="correct_velocity";

 int finest_level=parent->finestLevel();

 if (project_option_momeqn(project_option)==1) {
  //do nothing
 } else if (project_option_momeqn(project_option)==0) {
  //do nothing
 } else
  amrex::Error("project_option invalid35");

 Vector<int> scomp;
 Vector<int> ncomp;
 int state_index;
 int ncomp_check;

  //num_materials_combine=1
 get_mm_scomp_solver(
  1,
  project_option,
  state_index,
  scomp,
  ncomp,
  ncomp_check);

 if (ncomp_check!=nsolve)
  amrex::Error("ncomp_check invalid");

 bool use_tiling=ns_tiling;

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  if ((localMF[macdest+dir]->nComp()!=nsolve)||
      (localMF[macsrc+dir]->nComp()!=nsolve)||
      (localMF[gp+dir]->nComp()!=nsolve))
   amrex::Error("invalid ncomp");
 }

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 if (project_option_is_valid(project_option)==1) {
  // do nothing
 } else
  amrex::Error("project option invalid");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  debug_ngrow(FACE_VAR_MF+dir,0,local_caller_string);
 }

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,1,local_caller_string);

 const Real* dx = geom.CellSize();

 MultiFab& S_new=get_new_data(state_index,slab_step+1);

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(S_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(S_new,use_tiling); mfi.isValid();++mfi) {
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

  FArrayBox& xdest=(*localMF[macdest])[mfi];
  FArrayBox& ydest=(*localMF[macdest+1])[mfi];
  FArrayBox& zdest=(*localMF[macdest+AMREX_SPACEDIM-1])[mfi];

  FArrayBox& xsrc=(*localMF[macsrc])[mfi];
  FArrayBox& ysrc=(*localMF[macsrc+1])[mfi];
  FArrayBox& zsrc=(*localMF[macsrc+AMREX_SPACEDIM-1])[mfi];

  FArrayBox& xgp=(*localMF[gp])[mfi];
  FArrayBox& ygp=(*localMF[gp+1])[mfi];
  FArrayBox& zgp=(*localMF[gp+AMREX_SPACEDIM-1])[mfi];

  FArrayBox& xface=(*localMF[FACE_VAR_MF])[mfi];
  FArrayBox& yface=(*localMF[FACE_VAR_MF+1])[mfi];
  FArrayBox& zface=(*localMF[FACE_VAR_MF+AMREX_SPACEDIM-1])[mfi];

  Vector<int> presbc;
  getBCArray_list(presbc,state_index,gridno,scomp,ncomp);
  if (presbc.size()!=nsolve*AMREX_SPACEDIM*2)
   amrex::Error("presbc.size() invalid");

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  for (int velcomp=0;velcomp<nsolve;velcomp++) {

    // declared in: NAVIERSTOKES_3D.F90
   fort_fluidsolidcor(
    im_elastic_map.dataPtr(),
    &num_FSI_outer_sweeps,
    &FSI_outer_sweeps,
    &level,
    &finest_level,
    &velcomp,
    &nsolve,
    &project_option,
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    presbc.dataPtr(),
    maskcov.dataPtr(),ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
    xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()),
    yface.dataPtr(),ARLIM(yface.loVect()),ARLIM(yface.hiVect()),
    zface.dataPtr(),ARLIM(zface.loVect()),ARLIM(zface.hiVect()),
    xgp.dataPtr(velcomp),ARLIM(xgp.loVect()),ARLIM(xgp.hiVect()),
    ygp.dataPtr(velcomp),ARLIM(ygp.loVect()),ARLIM(ygp.hiVect()),
    zgp.dataPtr(velcomp),ARLIM(zgp.loVect()),ARLIM(zgp.hiVect()),
    xsrc.dataPtr(velcomp),ARLIM(xsrc.loVect()),ARLIM(xsrc.hiVect()),
    ysrc.dataPtr(velcomp),ARLIM(ysrc.loVect()),ARLIM(ysrc.hiVect()),
    zsrc.dataPtr(velcomp),ARLIM(zsrc.loVect()),ARLIM(zsrc.hiVect()),
    xdest.dataPtr(velcomp),ARLIM(xdest.loVect()),ARLIM(xdest.hiVect()),
    ydest.dataPtr(velcomp),ARLIM(ydest.loVect()),ARLIM(ydest.hiVect()),
    zdest.dataPtr(velcomp),ARLIM(zdest.loVect()),ARLIM(zdest.hiVect()),
    xlo,
    dx,
    &dt_slab,
    &cur_time_slab);
  } // velcomp=0..nsolve-1
 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_FLUIDSOLIDCOR,"correct_velocity");

} // end subroutine correct_velocity

//called from: NavierStokes::multiphase_project
void NavierStokes::residual_correction_form(
  int homflag_residual_correction_form,
  int energyflag,
  int project_option,int nsolve) {

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 if (project_option_momeqn(project_option)==1) {
  //do nothing
 } else if (project_option_momeqn(project_option)==0) {
  //do nothing
 } else
  amrex::Error("project_option_momeqn invalid36");

 if (project_option==SOLVETYPE_INITPROJ) { 

  if (homflag_residual_correction_form==1) {
   // do nothing
  } else
   amrex::Error("expecting homflag_residual_correction_form==1");

 } else if (project_option_is_valid(project_option)==1) {

   // -dt grad p face_weight  
   // SEM BC if enable_spectral==1
  int simple_AMR_BC_flag=0;
  int simple_AMR_BC_flag_viscosity=0;
   // apply_pressure_grad declared in NavierStokes2.cpp
  apply_pressure_grad(
   simple_AMR_BC_flag,
   simple_AMR_BC_flag_viscosity,
   homflag_residual_correction_form,//unused except for:SOLVETYPE_VISC,HEAT
   energyflag,
   GRADPEDGE_MF,
   STATE_FOR_RESID_MF,
   project_option,nsolve,
   dt_slab); //calling from residual_correction_form

   // UMAC_MF-=GRADPEDGE_MF
  correct_velocity(project_option,
   UMAC_MF, UMAC_MF, GRADPEDGE_MF,nsolve);

 } else
  amrex::Error("project_option invalid residual_correction_form");

}  // end subroutine residual_correction_form


// local_MF[idx_phi]=0 on all levels on input
void NavierStokes::mg_cycleALL(int presmooth,
 int project_option,
 int idx_rhs,
 int idx_phi,
 int nsolve) {

#if (profile_solver==1)
 std::string subname="NavierStokes::mg_cycleALL";
 std::stringstream popt_string_stream(std::stringstream::in |
      std::stringstream::out);
 popt_string_stream << project_option;
 std::string profname=subname+popt_string_stream.str();

 BLProfiler bprof(profname);
#endif 

 int finest_level=parent->finestLevel();

 if (level==finest_level) {
  // do nothing
 } else
  amrex::Error("level invalid mg_cycleALL");

 if (project_option_momeqn(project_option)==1) {
  //do nothing
 } else if (project_option_momeqn(project_option)==0) {
  //do nothing
 } else
  amrex::Error("project_option_momeqn(project_option) invalid37");

   // residmf represents the current status of:
   // f+ div grad p^init + a(POLDHOLD)  (POLDHOLD=p^*-p^init)
   // originally:
   // a(p-p*)-div grad p = f
   // let dp=p-p^init
   // a(dp+p^init-p^*)- div grad (dp+p^init) = f 
   // a dp - div grad dp = f+a(p^*-p^init)+div grad p^init=f+a POLDHOLD +
   // div grad p^init
 MultiFab* residmf=new MultiFab(grids,dmap,nsolve,0,
	MFInfo().SetTag("residmf"),FArrayBoxFactory());
 MultiFab::Copy(*residmf,*localMF[idx_rhs],0,0,nsolve,0);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  setVal_localMF(UMACSTAR_MF+dir,0.0,0,nsolve,0);
 } 

#if (profile_solver==1)
 bprof.stop();
#endif

 relaxLEVEL(
   residmf,
   idx_rhs,
   idx_phi,
   presmooth,
   project_option,nsolve);

 delete residmf;
} // end subroutine mg_cycleALL

// relaxLEVEL called from mg_cycleALL from multiphase_preconditioner from
//   multiphase_SHELL_preconditioner from multiphase_project
// this recursive routine first called from the finest_level.
// localMF[idx_phi]=0.0 on all levels prior to the first call of this
// routine.
void NavierStokes::relaxLEVEL(
  MultiFab* rhsmf,
  int idx_rhs,
  int idx_phi,
  int presmooth,
  int project_option,int nsolve) {

 std::string local_caller_string="relaxLEVEL";

#if (profile_solver==1)
 std::string subname="NavierStokes::relaxLEVEL";
 std::stringstream popt_string_stream(std::stringstream::in |
      std::stringstream::out);
 popt_string_stream << project_option;
 std::string profname=subname+popt_string_stream.str();
 profname=profname+"_";
 std::stringstream lev_string_stream(std::stringstream::in |
      std::stringstream::out);
 lev_string_stream << level;
 profname=profname+lev_string_stream.str();

 BLProfiler bprof(profname);
#endif

 if (project_option_momeqn(project_option)==1) {
  // do nothing
 } else if (project_option_momeqn(project_option)==0) {
  // do nothing
 } else
  amrex::Error("project_option invalid38");

 MultiFab* pbdry=new MultiFab(grids,dmap,nsolve,1,
	MFInfo().SetTag("pbdry"),FArrayBoxFactory());

  // if mask=1 (uncovered cell), rhsmf=idx_rhs-div u/dt
  // otherwise, rhsmf=rhsmf
  // idx_phi[level]=0.0 on input.
  // rhsmf=( idx_phi )*alpha -
  //       vol div UMACSTAR + idx_rhs=
  //       -vol div UMACSTAR + idx_rhs
  // a p - div grad p = f + a p^adv 
  // p=dp + p^init
  // (a)dp - div grad dp = f - (a)p^init + div grad p^init + a p^adv 
  //                          
  // p^init=-POLDHOLD+p^adv
  // p^last=p^adv-POLDHOLD
  // a dp - div grad dp = f- (a)(-POLDHOLD+p^adv) +
  //       div grad p^init + a p^adv =
  //       f- (a)(-POLDHOLD+p^adv) +
  //       div grad p^init + a p^adv =
  //       f+ div grad p^init + a(POLDHOLD)=
  //       f- div UMACSTAR + a(POLDHOLD)
  // 
  // on the finest level: UMACSTAR = 0 and idx_phi = 0.0
  // on coarser levels: idx_phi=0.0
 int homflag_apply_div=3;
 apply_div(
  project_option, 
  homflag_apply_div,
  idx_phi,
  rhsmf,
  localMF[idx_rhs],
  UMACSTAR_MF,
  nsolve,
  dt_slab); //calling from: relaxLEVEL

  // sets pbdry
  // going down the V-cycle: pdry=localMF[idx_phi]=0.0
 int homflag_down_V1=1; 
 applyBC_MGLEVEL(idx_phi,pbdry,homflag_down_V1,nsolve,project_option);

#if (profile_solver==1)
 bprof.stop();
#endif

  // the smoother uses A_LOW: e.g.
  // z^{k+1}=z^{k}+D_LOW^{-1}(r-A_LOW z^{k})
 if (level>0) {

#if (profile_solver==1)
  bprof.start();
#endif

  NavierStokes& ns_coarse=getLevel(level-1);
  ns_coarse.setVal_localMF(idx_phi,0.0,0,nsolve,1);

  for (int i=presmooth;i>0;i--) {
   int apply_lev_presmooth=0;
   mac_op->smooth(*localMF[idx_phi],*rhsmf,
    apply_lev_presmooth,*pbdry,bcpres_array,smooth_type);
  }
  MultiFab* residmf=new MultiFab(grids,dmap,nsolve,0,
	MFInfo().SetTag("residmf"),FArrayBoxFactory());

  int bfact=parent->Space_blockingFactor(level);
  if ((bfact<1)||(bfact>64))
   amrex::Error("bfact out of range");

   // low order BC
  int apply_lev_BC=0;
  mac_op->applyBC(*localMF[idx_phi],apply_lev_BC,*pbdry,bcpres_array);

   // gradpedge= -dt * grad p * denedgebc * densolidedgebc =
   //            -dt * grad p * face_weight_stable
  int homflag_down_V2=1;
  int energyflag=SUB_OP_FOR_MAIN;

    // we must have 
    // simple_AMR_BC_flag=1 and
    // simple_AMR_BC_flag_viscosity=1 
    // for otherwise inhomogeneous values from level+1 
    // (corresponding to a different rhs) will be
    // used in the stencil for elements neighboring a level+1
    // element.
    // It is ok that this is a "surrogate" gradient since 
    // "UMACSTAR" values are discarded after the preconditioner is
    // finished.
    // NOTE: the "simple" coarse/fine discretization is good enough
    // since localMF[idx_phi]=0.0 uniformly on level-1.
    // GRADPEDGE=-dt gradp/rho
  int simple_AMR_BC_flag=1;
  int simple_AMR_BC_flag_viscosity=1;

  apply_pressure_grad(
   simple_AMR_BC_flag,
   simple_AMR_BC_flag_viscosity,
   homflag_down_V2,
   energyflag,
   GRADPEDGE_MF,
   idx_phi,
   project_option,nsolve,
   dt_slab); //calling from relaxLEVEL

    // residmf=rhsmf-( (alpha+da) * phi - vol div grad phi )
  int apply_lev_resid=0;
  mac_op->residual(*residmf,*rhsmf,*localMF[idx_phi],
    apply_lev_resid,
    *pbdry,bcpres_array);
   // update resid where mask=1.  This step is necessary if
   // the complete operator differs from the simplified "mac_op"
   // operator.
   // residmf=-phi * (alpha+da) - div u + rhsmf=
   //         rhsmf-(phi*(alpha+da)+divu)
  int homflag_down_V3=2;
  apply_div(
   project_option, 
   homflag_down_V3,
   idx_phi,
   residmf,  // called "rhsmf" in apply_div
   rhsmf,    // called "diffusionRHScell" in apply_div
   GRADPEDGE_MF,
   nsolve,
   dt_slab); //calling from relaxLEVEL

   // UMACSTAR=UMACSTAR+GRADPEDGE
  correct_velocity(project_option,
    UMACSTAR_MF, UMACSTAR_MF, GRADPEDGE_MF,nsolve);

  MultiFab* residmf_coarse=
   new MultiFab(ns_coarse.grids,ns_coarse.dmap,nsolve,0,
	MFInfo().SetTag("residmf_coarse"),FArrayBoxFactory());
  residmf_coarse->setVal(0.0,0,nsolve,0);

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
   ns_coarse.setVal_localMF(UMACSTAR_MF+dir,0.0,0,nsolve,0);

  int ncomp_edge=-1;
  int scomp_edge=0;
  int start_dir=0;
   //order determined from enable_spectral
  int spectral_override=SPECTRAL_ORDER_AVGDOWN;
   // average down from level to level-1.
  ns_coarse.avgDownEdge_localMF(
    UMACSTAR_MF,
    scomp_edge,ncomp_edge,
    start_dir,AMREX_SPACEDIM,spectral_override,
    local_caller_string);

  int iavg=0;
  BoxArray& fgridscen=grids;
  DistributionMapping& fdmap=dmap;
  BoxArray& cgridscen=ns_coarse.grids;
  for (int veldir=0;veldir<nsolve;veldir++) {
    // average down from level to level-1.
    // calls fort_average which is low order.
   ns_coarse.Allaverage(  
    residmf_coarse,residmf,
    cgridscen,fgridscen,
    fdmap,
    iavg,
    veldir,veldir);
  }

#if (profile_solver==1)
  bprof.stop();
#endif

   // residmf_coarse becomes "rhsmf" on the coarser level.
  ns_coarse.relaxLEVEL(residmf_coarse,
   idx_rhs,idx_phi,
   presmooth,project_option,nsolve);

  delete residmf;
  delete residmf_coarse;

   // prolongates from the coarse phi_array
   // going up the V-cycle.
  int homflag_up_V1=0; 
  applyBC_MGLEVEL(idx_phi,pbdry,homflag_up_V1,nsolve,project_option);

  for (int i=presmooth;i>0;i--) {
   int apply_lev_post_smooth=0;
   mac_op->smooth(*localMF[idx_phi],*rhsmf,
    apply_lev_post_smooth,*pbdry,bcpres_array,smooth_type);
  }

 } else if (level==0) {

#if (profile_solver==1)
  bprof.start();
#endif

  int is_bottom=0;
  int usecg_at_bottom=1;
  int temp_meets_tol=0;
  Real error0=0.0;
  Real bottom_bottom_tol=save_atol_b*bottom_bottom_tol_factor;

  int lev0_cycles=0;

#if (profile_solver==1)
  bprof.stop();
#endif

  int apply_lev_cg_solve=0;
  mac_op->CG_solve(
   lev0_cycles,
   verbose,is_bottom,
   *localMF[idx_phi],*rhsmf,
   save_atol_b,
   bottom_bottom_tol,
   *pbdry,
   bcpres_array,
   usecg_at_bottom,
   temp_meets_tol,
   smooth_type,
   bottom_smooth_type,
   presmooth,presmooth,  //postsmooth="presmooth"
   error0,
   apply_lev_cg_solve);

#if (profile_solver==1)
  bprof.start();
#endif

  int current_list_size=lev0_cycles_list.size();
  int placeholder=-1;
  for (int i=current_list_size-1;i>=0;i--) {
   if (lev0_cycles>lev0_cycles_list[i])
    placeholder=i;
  }

  if (placeholder==-1)
   placeholder=current_list_size;

  lev0_cycles_list.resize(current_list_size+1);
  for (int i=current_list_size;i>placeholder;i--)
   lev0_cycles_list[i]=lev0_cycles_list[i-1];
  lev0_cycles_list[placeholder]=lev0_cycles;

#if (profile_solver==1)
  bprof.stop();
#endif

 } else
  amrex::Error("level invalid relaxLEVEL");

 delete pbdry;

} // subroutine relaxLEVEL

void NavierStokes::check_outer_solver_convergence(
	Real error_n,Real error0,
        Real save_mac_abs_tol_in,
        int& meets_tol) {

 if (error0>=0.0) {
  // do nothing
 } else
  amrex::Error("error0 must be nonnegative");

 if (error_n<=save_mac_abs_tol_in) {
  meets_tol=1;
 } else if (error_n>=save_mac_abs_tol_in) {
  // do nothing
 } else
  amrex::Error("error_n invalid");

} // subroutine check_outer_solver_convergence

// If called from NavierStokes::multiphase_project,
// update_vel=1 if called at the beginning of each outer_iter sweep.
//
// If called from NavierStokes::multiphase_preconditioner,
// update_vel=0 if this routine used as a preconditioner (instead of MG).
void NavierStokes::jacobi_cycles(
 int call_adjust_tolerance,
 int ncycles,
 int update_vel,int project_option,
 int idx_mac_rhs_crse,
 int idx_mac_phi_crse,
 Real& error_at_the_beginning,
 Real& error_after_all_jacobi_sweeps,
 Real& error0,
 Real& error0_max,
 int krylov_subspace_num_outer_iterSOLVER,
 int nsolve) {

 int finest_level=parent->finestLevel();

 if (level==0) {
  // do nothing
 } else
  amrex::Error("jacobi_cycles should only be called from level==0");

 error_at_the_beginning=0.0;
 error_after_all_jacobi_sweeps=0.0;

 if ((nsolve!=1)&&(nsolve!=AMREX_SPACEDIM))
  amrex::Error("nsolve invalid33");

 if (project_option_momeqn(project_option)==1) {
  //do nothing
 } else if (project_option_momeqn(project_option)==0) {
  //do nothing
 } else
  amrex::Error("project_option invalid39");

 if ((krylov_subspace_num_outer_iterSOLVER>=0)&&
     (krylov_subspace_num_outer_iterSOLVER<=
      krylov_subspace_max_num_outer_iter)) {
  //do nothing
 } else
  amrex::Error("krylov_subspace_num_outer_iterSOLVER invalid");

 allocate_array(0,nsolve,-1,RESID_MF);

 int temp_ncycles=ncycles;
 if (ncycles==0) 
  temp_ncycles=1;

 for (int vcycle_jacobi=0;vcycle_jacobi<temp_ncycles;vcycle_jacobi++) {

   // null space filtered out of residual in this routine.
  residALL(project_option,idx_mac_rhs_crse,
    RESID_MF,idx_mac_phi_crse,nsolve);

  if (update_vel==1) {  // not called as a preconditioner

   Real local_error;
   dot_productALL(project_option,RESID_MF,RESID_MF,local_error,nsolve);
   local_error=std::sqrt(local_error);
   if (verbose>0) {
    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "vcycle_jacobi,error " << vcycle_jacobi << ' ' << 
      local_error << '\n';
    }
    if (1==0) {
     for (int ilev=finest_level;ilev>=level;ilev--) {
      NavierStokes& ns_level=getLevel(ilev);
      const Real* dx = ns_level.geom.CellSize();
      std::cout << "ilev=" << ilev << '\n';
      IntVect test_index=ns_level.localMF[RESID_MF]->maxIndex(0,0);
      Real test_pos[AMREX_SPACEDIM];
      for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
       test_pos[dir]=geom.ProbLo(dir)+(test_index[dir]+0.5)*dx[dir];
       std::cout << "pos(maxIndex)=" << test_pos[dir] << '\n';
      }
      Real test_val=ns_level.localMF[RESID_MF]->max(0,0);
      std::cout << "max_index " << test_index << '\n';
      std::cout << "max_value " << test_val << '\n';
      test_index=ns_level.localMF[RESID_MF]->minIndex(0,0);
      for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
       test_pos[dir]=geom.ProbLo(dir)+(test_index[dir]+0.5)*dx[dir];
       std::cout << "pos(minIndex)=" << test_pos[dir] << '\n';
      }
      test_val=ns_level.localMF[RESID_MF]->min(0,0);
      std::cout << "min_index " << test_index << '\n';
      std::cout << "min_value " << test_val << '\n';
     } // ilev
     std::fflush(NULL);
    } // debugging

   }  // verbose>0

   if (vcycle_jacobi==0)
    error_at_the_beginning=local_error;

    //very first "outer iteration"
   if (krylov_subspace_num_outer_iterSOLVER==0) {

    error0=local_error;

    if (call_adjust_tolerance==1) {
     adjust_tolerance(error0,error0_max,project_option);
    } else if (call_adjust_tolerance==0) {
     // do nothing
    } else
     amrex::Error("call_adjust_tolerance invalid");

   } else if ((krylov_subspace_num_outer_iterSOLVER>0)&&
              (krylov_subspace_num_outer_iterSOLVER<=
               krylov_subspace_max_num_outer_iter)) {
    // do nothing
   } else
    amrex::Error("krylov_subspace_num_outer_iterSOLVER invalid");

  } else if (update_vel==0) { // called as preconditioner
    // do nothing
  } else
   amrex::Error("update_vel invalid");

  if (ncycles>0) {
    // 1. (begin) project_right_hand_side(idx_mac_phi_crse)
    // 2. (end)   project_right_hand_side(idx_mac_phi_crse)
   JacobiALL(RESID_MF,idx_mac_rhs_crse,
     idx_mac_phi_crse,project_option,nsolve); 
  } else if (ncycles==0) {
   // do nothing
  } else
   amrex::Error("ncycles invalid");

 }  // vcycle_jacobi=0..temp_ncycles-1

 if (ncycles>0) {

   // null space filtered out.
  residALL(project_option,idx_mac_rhs_crse,
    RESID_MF,idx_mac_phi_crse,nsolve);

  if (update_vel==1) {  // not called as a preconditioner

   Real local_error;
   dot_productALL(project_option,RESID_MF,RESID_MF,local_error,nsolve);
   local_error=std::sqrt(local_error);
   if (verbose>0) {
    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "jacobi after last cycle,error " << local_error << '\n';
    }
   }
   error_after_all_jacobi_sweeps=local_error;
  } else if (update_vel==0) {
   // do nothing
  } else
   amrex::Error("update_vel invalid");

 } else if (ncycles==0) {
  // do nothing
 } else
  amrex::Error("ncycles invalid");
 
 delete_array(RESID_MF);

}  // end subroutine jacobi_cycles

// called from:
//  NavierStokes::multiphase_project
void NavierStokes::updatevelALL(
 int project_option,
 int idx_mac_phi_crse,int nsolve) {

 int finest_level=parent->finestLevel();

 if (project_option_momeqn(project_option)==1) {
  //do nothing
 } else if (project_option_momeqn(project_option)==0) {
  //do nothing
 } else
  amrex::Error("project_option_momeqn(project_option)  invalid40");

 // calling from updatevelALL from multiphase_project
 // gradpedge=-dt W grad p
 // NavierStokes::applyGradALL is declared in MacProj.cpp
 applyGradALL(project_option,idx_mac_phi_crse,nsolve);

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
     // POLDHOLD-=idx_mac_phi_crse
     // PNEW+=idx_mac_phi_crse
     // UMAC+=GRADPEDGE  (GRADPEDGE=-dt W grad p)
  ns_level.mac_update(ns_level.localMF[idx_mac_phi_crse],
    project_option,nsolve);
  ns_level.setVal_localMF(idx_mac_phi_crse,0.0,0,nsolve,1);

  if (ilev<finest_level) {
   ns_level.avgDownMac();   // works on UMAC_MF
  }
 }

} // end subroutine updatevelALL


void NavierStokes::Prepare_UMAC_for_solver(int project_option,
  int nsolve) {

 if (dt_slab>0.0) {
  // do nothing
 } else
  amrex::Error("dt_slab invalid:Prepare_UMAC_for_solver");

 if ((nsolve!=1)&&(nsolve!=AMREX_SPACEDIM))
  amrex::Error("nsolve invalid");

 int finest_level=parent->finestLevel();

 if ((level>=0)&&(level<=finest_level)) {
  // do nothing
 } else
  amrex::Error("level invalid Prepare_UMAC_for_solver");

 if (project_option_momeqn(project_option)==1) {
  //do nothing
 } else if (project_option_momeqn(project_option)==0) {
  //do nothing
 } else
  amrex::Error("project_option_momeqn invalid41");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  new_localMF(MAC_TEMP_MF+dir,nsolve,0,dir);
  setVal_localMF(MAC_TEMP_MF+dir,0.0,0,nsolve,0);
 } // dir

 new_localMF(DIFFUSIONRHS_MF,nsolve,0,-1);
 if (project_option==SOLVETYPE_PRESEXTRAP) {
  setVal_localMF(DIFFUSIONRHS_MF,0.0,0,nsolve,0);
 } else if (project_option==SOLVETYPE_PRES)  { 
  int scomp=0;
   // MDOT_MF already premultiplied by the cell volume
  Copy_localMF(DIFFUSIONRHS_MF,MDOT_MF,0,scomp,nsolve,0);
 } else if (project_option==SOLVETYPE_INITPROJ) { 
  int scomp=0;
   // MDOT_MF already premultiplied by the cell volume
   // MDOT_MF might not be zero at t=0 if sources/sinks of mass.
   // (but for now, MDOT_MF should be zero)
  Copy_localMF(DIFFUSIONRHS_MF,MDOT_MF,0,scomp,nsolve,0);
  zero_independent_variable(project_option,nsolve);
 } else if (project_option==SOLVETYPE_HEAT) { 

  zero_independent_vel(project_option,UMAC_MF,nsolve);
  setVal_localMF(DIFFUSIONRHS_MF,0.0,0,nsolve,0);

 } else if (project_option==SOLVETYPE_VISC) {

  zero_independent_vel(project_option,UMAC_MF,nsolve);
  setVal_localMF(DIFFUSIONRHS_MF,0.0,0,nsolve,0);

 } else if ((project_option>=SOLVETYPE_SPEC)&&
	    (project_option<SOLVETYPE_SPEC+num_species_var)) {

  zero_independent_vel(project_option,UMAC_MF,nsolve);
  setVal_localMF(DIFFUSIONRHS_MF,0.0,0,nsolve,0);

 } else
  amrex::Error("project_option invalid prepare_UMAC_for_solver");

} // subroutine Prepare_UMAC_for_solver

void NavierStokes::remove_UMAC_for_solver(int project_option) {

 if (project_option_is_valid(project_option)==1) {
  // do nothing
 } else
  amrex::Error("project_option_is_valid(project_option)==1 failed");

 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid remove_UMAC_for_solver");

 delete_localMF(DIFFUSIONRHS_MF,1);
 delete_localMF(MAC_TEMP_MF,AMREX_SPACEDIM);

}

void NavierStokes::multiphase_SHELL_preconditioner(
 int project_option,int project_timings,
 int presmooth,int postsmooth,
 int idx_Z,int idx_R,int nsolve) {

 if (level==0) {
  // do nothing
 } else
  amrex::Error("level invalid");

 if (override_bc_to_homogeneous==1) {
  // do nothing
 } else
  amrex::Error("expecting homogeneous mode in preconditioner");

 if (project_option_momeqn(project_option)==1) {
  //do nothing
 } else if (project_option_momeqn(project_option)==0) {
  //do nothing
 } else
  amrex::Error("project_option invalid44");

 int change_flag=0;
 project_right_hand_side(idx_R,project_option,change_flag);

   // Z=M^{-1}R
   // Z=project(Z)
   // the first command in multiphase_precond: zeroALL(1,nsolve,idx_Z) 
 multiphase_preconditioner(
   project_option,project_timings,
   presmooth,postsmooth,
   idx_Z,idx_R,nsolve);

 project_right_hand_side(idx_Z,project_option,change_flag);

} // end subroutine multiphase_SHELL_preconditioner

void NavierStokes::multiphase_preconditioner(
 int project_option,int project_timings,
 int presmooth,int postsmooth,
 int idx_Z,int idx_R,int nsolve) {

#if (profile_solver==1)
 std::string subname="NavierStokes::multiphase_preconditioner";
 std::stringstream popt_string_stream(std::stringstream::in |
      std::stringstream::out);
 popt_string_stream << project_option;
 std::string profname=subname+popt_string_stream.str();

 BLProfiler bprof(profname);
#endif

 int finest_level=parent->finestLevel();

 if (level==0) {
  // do nothing
 } else
  amrex::Error("multiphase_precond should only be called from level==0");

 double before_cycle=0.0;
 if (project_timings==1)
  before_cycle=ParallelDescriptor::second();

 if ((nsolve!=1)&&(nsolve!=AMREX_SPACEDIM))
  amrex::Error("nsolve invalid");
 
 if (project_option_momeqn(project_option)==1) {
  //do nothing
 } else if (project_option_momeqn(project_option)==0) {
  //do nothing
 } else
  amrex::Error("project_option_momeqn invalid42");

  //declared in:NavierStokes2.cpp
 zeroALL(1,nsolve,idx_Z);

#if (profile_solver==1)
 bprof.stop();
#endif

  // PCG (jacobi preconditioner)
 if (project_solver_type==1) {

#if (profile_solver==1)
  bprof.start();
#endif

  int smooth_cycles=presmooth+postsmooth;
  int update_vel=0;
  Real error_at_the_beginning=0.0;
  Real error_after_all_jacobi_sweeps=0.0;
  Real error_n=0.0;
  Real error0_max=0.0;
  int krylov_subspace_num_outer_iterSOLVER=0;

  int call_adjust_tolerance=0;

    //calling from: multiphase_preconditioner
  jacobi_cycles(
   call_adjust_tolerance, //=0
   smooth_cycles, //smooth_cycles=presmooth+postsmooth
   update_vel,  // =0 (i.e. jacobi_cycles used as a preconditioner)
   project_option,
   idx_R,idx_Z,
   error_at_the_beginning,
   error_after_all_jacobi_sweeps,
   error_n,    // not modified
   error0_max, // not modified
   krylov_subspace_num_outer_iterSOLVER,nsolve);

#if (profile_solver==1)
  bprof.stop();
#endif

   // MGPCG
 } else if (project_solver_type==0) {

  NavierStokes& ns_finest=getLevel(finest_level);
  ns_finest.mg_cycleALL(presmooth,
  project_option,
  idx_R,
  idx_Z,
  nsolve);

   // MINV=I
 } else if (project_solver_type==2) {

#if (profile_solver==1)
  bprof.start();
#endif

  int ncomp=localMF[idx_Z]->nComp();
  int ngrow=localMF[idx_Z]->nGrow();

  if ((ncomp==1)||
      (ncomp==AMREX_SPACEDIM)) {
   // do nothing
  } else
   amrex::Error("expecting ncomp=1 or sdim");

   //ngrow,scomp,ncomp
  setVal_array(ngrow,0,ncomp,0.0,idx_Z);
  Copy_array(idx_Z,idx_R,0,0,ncomp,0);

#if (profile_solver==1)
  bprof.stop();
#endif

 } else
  amrex::Error("project_solver_type invalid");

 int change_flag=0;
 project_right_hand_side(idx_Z,project_option,change_flag);

 double after_cycle=0.0;
 if (project_timings==1) {
  after_cycle=ParallelDescriptor::second();
  if (ParallelDescriptor::IOProcessor())
   std::cout << "cycle time " << after_cycle-before_cycle << '\n';
 }

} // subroutine multiphase_preconditioner

void NavierStokes::set_local_tolerances(int project_option) {

 save_atol_b=CPP_EPS_14_6;
 save_mac_abs_tol=mac_abs_tol;
 save_min_rel_error=minimum_relative_error;
 ParmParse pp("mg");

 if (project_option_singular_possible(project_option)==1) {
  save_mac_abs_tol=mac_abs_tol;
  save_atol_b=0.01*save_mac_abs_tol;
  pp.queryAdd("bot_atol",save_atol_b);
  save_min_rel_error=minimum_relative_error;
 } else if (project_option==SOLVETYPE_HEAT) {
  save_mac_abs_tol=thermal_abs_tol;
  save_atol_b=0.01*save_mac_abs_tol; 
  pp.queryAdd("thermal_bot_atol",save_atol_b);
  save_min_rel_error=diffusion_minimum_relative_error;
 } else if (project_option==SOLVETYPE_VISC) {
  save_mac_abs_tol=visc_abs_tol;
  save_atol_b=0.01*save_mac_abs_tol; 
  pp.queryAdd("visc_bot_atol",save_atol_b);
  save_min_rel_error=diffusion_minimum_relative_error;
 } else if ((project_option>=SOLVETYPE_SPEC)&&
	    (project_option<SOLVETYPE_SPEC+num_species_var)) {
  save_mac_abs_tol=visc_abs_tol;
  save_atol_b=0.01*save_mac_abs_tol; 
  pp.queryAdd("visc_bot_atol",save_atol_b);
  save_min_rel_error=diffusion_minimum_relative_error;
 } else
  amrex::Error("project_option invalid 51");

} // end subroutine set_local_tolerances

void NavierStokes::Krylov_checkpoint(
  int vcycle,
  Real krylov_error,
  Real& best_error,
  int& best_iter,
  int idx_phi,
  int idx_umac,
  int& restart_flag) {

 std::string local_caller_string="Krylov_checkpoint";

 if (level==0) {
  // do nothing
 } else
  amrex::Error("level invalid");

 Real breakdown_factor=1000.0;

 int update_best=0;
 if (vcycle==-1) {
  update_best=1;
  best_error=krylov_error;
  best_iter=vcycle;
  if (krylov_error>=0.0) {
   // do nothing
  } else
   amrex::Error("krylov_error invalid");
 } else if (vcycle>=0) {
  if ((krylov_error>=0.0)&&
      (krylov_error<best_error)) {
   update_best=1;
   best_error=krylov_error;
   best_iter=vcycle;
  } else if ((krylov_error>=best_error)&&
             (krylov_error<=breakdown_factor*best_error)) {
   // do nothing
  } else if (krylov_error>=breakdown_factor*best_error) {
   restart_flag=1;
  } else if (krylov_error<0.0) {
   amrex::Error("krylov_error must be positive");
  } else {
   restart_flag=1; // krylov_error=NaN
  }
 } else
  amrex::Error("vcycle invalid");
 
 int nsolve=localMF[idx_phi]->nComp();
 if ((nsolve==1)||(nsolve==AMREX_SPACEDIM)) {
  if (update_best==1) {
   Copy_array(RESTART_MAC_PHI_CRSE_MF,idx_phi,0,0,nsolve,1);

   for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
    Copy_array(RESTART_UMACSTAR_MF+dir,idx_umac+dir,0,0,nsolve,0);
   }
  } else if (update_best==0) {
   // do nothing
  } else
   amrex::Error("update_best invalid");
 } else
  amrex::Error("nsolve invalid");
    
} // end subroutine NavierStokes::Krylov_checkpoint

// if project_option==SOLVETYPE_PRES,
// RHS= -(1/dt)(a_{i+1/2}u_{i+1/2}-a_{i-1/2}u_{i-1/2}+...)+
//  alpha_{i}poldhold_{i}  + diffusionRHS
//
// alpha_{i}p_{i}-(bx_{i+1/2} (p_{i+1}-p_{i})-bx_{i-1/2} (p_{i}-p_{i-1})) = RHS
//
// if pressure solve,
//   ns_level.init_EOS_pressure() *DOES NOT* overwrite P01 with p(rho,T)
//   Result of the Helmholtz solve is stored instead.
//
void NavierStokes::multiphase_project(int project_option) {

 std::string local_caller_string="multiphase_project";

 int finest_level=parent->finestLevel();
 if (level!=0)
  amrex::Error("level invalid multiphase_project");

 init_rest_fraction(local_caller_string);

 std::fflush(NULL);

#if (profile_solver==1)
 std::string subname="NavierStokes::multiphase_project";
 std::stringstream popt_string_stream(std::stringstream::in |
      std::stringstream::out);
 popt_string_stream << project_option;
 std::string profname=subname+popt_string_stream.str();
 BLProfiler bprof(profname);
#endif

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid multiphase_project");

 if ((slab_step>=0)&&(slab_step<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("slab_step invalid");

 const Real* coarse_dx=geom.CellSize();

  // initial pressure needed for FSI (very first time step)
 if (project_option==SOLVETYPE_INITPROJ) {

   //ngrow,ncomp,grid_type
  allocate_array(1,1,-1,PRESSURE_SAVE_MF);
  Copy_array(PRESSURE_SAVE_MF,GET_NEW_DATA_OFFSET+State_Type,
	  STATECOMP_PRES,0,STATE_NCOMP_PRES,1);

 } else if (project_option==SOLVETYPE_PRESEXTRAP) {

  allocate_array(1,1,-1,PRESSURE_SAVE_MF);
  Copy_array(PRESSURE_SAVE_MF,GET_NEW_DATA_OFFSET+State_Type,
	  STATECOMP_PRES,0,STATE_NCOMP_PRES,1);

 } else if (project_option_is_valid(project_option)==1) {
  // do not save anything
 } else
  amrex::Error("project_option invalid 9692");

 int save_enable_spectral=enable_spectral;

 if (project_option_projection(project_option)==1) {
  //do nothing
 } else if (project_option==SOLVETYPE_PRESEXTRAP) {
  override_enable_spectral(0); // always low order
 } else if (project_option==SOLVETYPE_HEAT) { 
 } else if (project_option==SOLVETYPE_VISC) {
 } else if ((project_option>=SOLVETYPE_SPEC)&&
	    (project_option<SOLVETYPE_SPEC+num_species_var)) { 
  override_enable_spectral(0); // always low order
 } else
  amrex::Error("project_option invalid43");

 int energyflag=SUB_OP_FOR_MAIN;

 if (project_option_momeqn(project_option)==1) {
  //do nothing
 } else if (project_option_momeqn(project_option)==0) {
  //do nothing
 } else
  amrex::Error("project_option invalid44");

 int nsolve=1;

  //SOLVETYPE_INITPROJ, 
  //SOLVETYPE_PRES
 if (project_option_projection(project_option)==1) {

  if (project_option==SOLVETYPE_INITPROJ) { 
   // do nothing
  } else if (project_option==SOLVETYPE_PRES) {
   //do nothing
  } else
   amrex::Error("project_option invalid 45"); 

 } else if (project_option==SOLVETYPE_PRESEXTRAP) { 
  // do nothing
 } else if (project_option==SOLVETYPE_HEAT) { 
  // do nothing
 } else if (project_option==SOLVETYPE_VISC) { 
  nsolve=AMREX_SPACEDIM;
 } else if ((project_option>=SOLVETYPE_SPEC)&&
            (project_option<SOLVETYPE_SPEC+num_species_var)) {
  // do nothing
 } else
  amrex::Error("project_option invalid multiphase_project");

 Vector<int> scomp;
 Vector<int> ncomp;
 int state_index;  
 int ncomp_check;
  //num_materials_combine=1
 get_mm_scomp_solver(
  1,
  project_option,
  state_index,
  scomp,
  ncomp,
  ncomp_check);
 if (ncomp_check!=nsolve)
  amrex::Error("nsolve or ncomp_check invalid 4904");

  // localMF[UMAC_MF] = 0.0
 allocate_MAC_velocityALL(nsolve,UMAC_MF);

 resize_metrics(1);
 debug_ngrow(VOLUME_MF,1,local_caller_string);
 debug_ngrow(FACE_VAR_MF,0,local_caller_string);

 int project_timings=0;

 double begin_project=0.0;
 if (project_timings==1)
  begin_project=ParallelDescriptor::second();


  // if project_option==SOLVETYPE_INITPROJ, then 
  // initial guess is p=0 and homogeneous BC are used
  // for p (the potential function).
  //
  // if project_option==SOLVETYPE_PRES, then:
  //  pressure=padvect  (P01)
  //
  // allocate_project_variables comes after this since 
  // allocate_project_variables does the following:
  // p^adv = outer_iter_pressure = p_new
  // p_init = outer_iter_pressure = p_new
  // POLDHOLD = p^adv - p_init
  // p_new = p_init
  // outer_iter_pressure=p_init

 for (int ilev=level;ilev<=finest_level;ilev++) {
   NavierStokes& ns_level=getLevel(ilev);
    // localMF[MAC_TEMP_MF]=0
    // diffusionRHS=mdot project_option==SOLVETYPE_PRES
    // diffusionRHS=mdot pressure=0 project_option==SOLVETYPE_INITPROJ
    // diffusionRHS=0.0 UMAC=0 project_option==SOLVETYPE_HEAT 
    // diffusionRHS=0.0 UMAC=0 project_option==SOLVETYPE_VISC 
    // diffusionRHS=0.0 UMAC=0 project_option==SOLVETYPE_SPEC ...
    // diffusionRHS=0.0 if project_option==SOLVETYPE_PRESEXTRAP 
   ns_level.Prepare_UMAC_for_solver(project_option,nsolve);
 }  // ilev=level ... finest_level

 Real max_nlevels=0;
 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
   //NavierStokes::NSnumLevels() is declared in: NavierStokes.cpp
   //ns_level.NSnumLevels()=min_{grids \in level} numLevels(grid)
  int nlevels=ns_level.NSnumLevels();
  if (nlevels>max_nlevels)
   max_nlevels=nlevels;
 } // ilev=level ... finest_level

 std::fflush(NULL);

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << " ---------------------------------------------- \n";
   print_project_option(project_option);
   std::cout << " divu_outer_sweeps= " << divu_outer_sweeps << '\n';
   std::cout << " SDC_outer_sweeps= " << SDC_outer_sweeps << '\n';
   std::cout << " ns_time_order= " << ns_time_order << '\n';
   std::cout << " slab_step= " << slab_step << '\n';
   std::cout << " dt_slab= " << dt_slab << '\n';
   std::cout << " lower_slab_time= " << lower_slab_time << '\n';
   std::cout << " upper_slab_time= " << upper_slab_time << '\n';
   std::cout << " prev_time_slab= " << prev_time_slab << '\n';
   std::cout << " cur_time_slab= " << cur_time_slab << '\n';
   std::cout << "max_nlevels= " << max_nlevels << '\n';
   std::cout << " ---------------------------------------------- \n";
  }
 } else if (verbose==0) {
  // do nothing
 } else
  amrex::Error("verbose invalid");

 std::fflush(NULL);


  // automatically initializes mac_phi_crse_array=0.0
 allocate_independent_var(nsolve,MAC_PHI_CRSE_MF);
 allocate_independent_var(nsolve,RESTART_MAC_PHI_CRSE_MF);
  // automatically initializes mac_rhs_crse_array=0.0
 allocate_rhs_var(nsolve,MAC_RHS_CRSE_MF);
 
  // currently in: multiphase_project
 min_face_wt.resize(thread_class::nthreads);
 max_face_wt.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  min_face_wt[tid].resize(NCOMP_FACE_WT);
  max_face_wt[tid].resize(NCOMP_FACE_WT);
  for (int iwt=0;iwt<NCOMP_FACE_WT;iwt++) {
   min_face_wt[tid][iwt]=1.0e+20;
   max_face_wt[tid][iwt]=-1.0e+20;
  }
 } // tid

  // in multiphase_project
 if (project_option==SOLVETYPE_PRES) {

  int potgrad_surface_tension_mask=POTGRAD_NULLOPTION;

  if ((FSI_outer_sweeps>=0)&&
      (FSI_outer_sweeps<num_FSI_outer_sweeps-1)) {

   if (incremental_gravity_flag==1) {
    potgrad_surface_tension_mask=POTGRAD_INCREMENTAL_GRAV;
   } else if (incremental_gravity_flag==0) {
    potgrad_surface_tension_mask=POTGRAD_BASE_GRAV;
   } else
    amrex::Error("incremental_gravity_flag invalid");

  } else if (FSI_outer_sweeps==num_FSI_outer_sweeps-1) {

   if (incremental_gravity_flag==1) {
    potgrad_surface_tension_mask=POTGRAD_SURFTEN_INCREMENTAL_GRAV;
   } else if (incremental_gravity_flag==0) {
    potgrad_surface_tension_mask=POTGRAD_SURFTEN_BASE_GRAV;
   } else
    amrex::Error("incremental_gravity_flag invalid");

  } else
   amrex::Error("FSI_outer_sweeps invalid");


   // 1. init_gravity_potential
   //      output: HYDROSTATIC_PRESDEN_MF
   // 2. process_potential_force_face
   //      output: POTENTIAL_FORCE_EDGE_MF (OP_POTGRAD_TO_MAC)
   //
  if (potgrad_surface_tension_mask==POTGRAD_NULLOPTION) {
   // do nothing
  } else if (potgrad_surface_tension_mask!=POTGRAD_NULLOPTION) {
   process_potential_forceALL(potgrad_surface_tension_mask,project_option);
  } else
   amrex::Error("potgrad_surface_tension_mask invalid");

// 1. overwrites cell/face velocity perhaps
// 2. must be called before adding gravity and surface tension.
// 3. cannot be called after the project because the velocity
//    will then fail to satisfy the discrete divergence condition.
  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   ns_level.overwrite_outflow();  
  }

   // increment_potential_forceALL is declared in NavierStokes2.cpp
   // increment_potential_force is declared in NavierStokes2.cpp
   // fort_addgravity is declared in NAVIERSTOKES_3D.F90
   // FUTURE: E+=dt u dot g + dt^2 g dot g/2
  if (potgrad_surface_tension_mask==POTGRAD_NULLOPTION) {
   // do nothing
  } else if (potgrad_surface_tension_mask!=POTGRAD_NULLOPTION) {
   increment_potential_forceALL(); 
  } else
   amrex::Error("potgrad_surface_tension_mask invalid");

  if (1==0) {
   int basestep_debug=nStep()+1;
   parent->writeDEBUG_PlotFile(
     basestep_debug,
     SDC_outer_sweeps,
     slab_step,
     divu_outer_sweeps);
   std::cout << "press any number then enter AFTER GRAV\n";
   int n_input;
   std::cin >> n_input;
  }  

  if (visual_buoyancy_plot_int>0) {

   // nsteps==0 very first step.
   // in: Multiphase_project
   int nsteps=parent->levelSteps(0); 

   if (very_last_sweep==1) {
    int ratio=(nsteps+1)/visual_buoyancy_plot_int;
    ratio=ratio*visual_buoyancy_plot_int;
    if (ratio==nsteps+1) {

     for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
      // filenames: "FACE_VAR<stuff>.plt" (MAC data)
      // see: " FACECOMP_< ... > " in EXTRAP_COMP.H
      writeSanityCheckData(
       "FACE_VAR",
       "project_option==SOLVETYPE_PRES:FACE_VAR_MF",
       local_caller_string,
       FACE_VAR_MF+dir, //tower_mf_id
       localMF[FACE_VAR_MF+dir]->nComp(),
       FACE_VAR_MF+dir,
       -1, // State_Type==-1
       dir,
       parent->levelSteps(0)); 


       // gravity * dt
      if (potgrad_surface_tension_mask==POTGRAD_NULLOPTION) {
       // do nothing
      } else if (potgrad_surface_tension_mask!=POTGRAD_NULLOPTION) {
       writeSanityCheckData(
        "POTENTIAL_FORCE_EDGE",
        "project_option==SOLVETYPE_PRES:POTENTIAL_FORCE_EDGE",
        local_caller_string,
        POTENTIAL_FORCE_EDGE_MF+dir, //tower_mf_id
        localMF[POTENTIAL_FORCE_EDGE_MF+dir]->nComp(),
        POTENTIAL_FORCE_EDGE_MF+dir,
        -1, // State_Type==-1
        dir,
        parent->levelSteps(0)); 
      } else
       amrex::Error("potgrad_surface_tension_mask invalid");

     } // dir=0..sdim-1

    } // ratio==nsteps+1
   } else if (very_last_sweep==0) {
    // do nothing
   } else
    amrex::Error("very_last_sweep invalid");

  } else if (visual_buoyancy_plot_int==0) {
   // do nothing
  } else
   amrex::Error("visual_buoyancy_plot_int invalid");

  if (potgrad_surface_tension_mask==POTGRAD_NULLOPTION) {
   // do nothing
  } else if (potgrad_surface_tension_mask!=POTGRAD_NULLOPTION) {
   deallocate_potential_forceALL(); 
  } else
   amrex::Error("potgrad_surface_tension_mask invalid");

   // grad p  face
   // u=u-(1/rho)(int gp - dt gp)
  if ((SDC_outer_sweeps>0)&&
      (SDC_outer_sweeps<ns_time_order)&&
      (divu_outer_sweeps+1==num_divu_outer_sweeps)) {

   if (project_option==SOLVETYPE_PRES) {

    for (int ilev=finest_level;ilev>=level;ilev--) {
     NavierStokes& ns_level=getLevel(ilev);
     //NavierStokes::make_SEM_delta_force declared in NavierStokes.cpp
     ns_level.make_SEM_delta_force(SOLVETYPE_PRES); 
    }

   } else
    amrex::Error("SOLVETYPE_PRES allowed only.");

  } else if (SDC_outer_sweeps==0) {
   // do nothing
  } else if (divu_outer_sweeps+1<num_divu_outer_sweeps) {
   // do nothing
  } else
   amrex::Error("SDC_outer_sweeps or divu_outer_sweeps invalid multiphase prj");

 } else if (project_option!=SOLVETYPE_PRES) {
  //do nothing
 } else
  amrex::Error("project_option bust");	 

 if (project_option_needs_scaling(project_option)==1) {

   // fortran pressure and velocity scales
   // dt_slab
   // s_new velocity
   // s_new pressure
   // div_new 
   // umac  velocity
   // face_var: rigid velocity
   // diffusion_rhs
   // solid_var velocity
   // second component of CELL_SOUND (padvect_avg)
  scale_variablesALL();

 } else if (project_option_needs_scaling(project_option)==0) {
  // do nothing
 } else
  amrex::Error("project_option_needs_scaling invalid46");

  //SOLVETYPE_PRES, 
  //SOLVETYPE_INITPROJ, 
 if (project_option_FSI_rigid(&project_option)==1) {

  Vector<blobclass> blobdata;
  Vector< Vector<Real> > mdot_data;
  Vector< Vector<Real> > mdot_comp_data;
  Vector< Vector<Real> > mdot_data_redistribute;
  Vector< Vector<Real> > mdot_comp_data_redistribute;
  Vector<int> type_flag;

  if (verbose>0) {
   if (ParallelDescriptor::IOProcessor()) {
    std::cout << "BEGIN: color_variable, multiphase_project\n";
    print_project_option(project_option);
   }
  }

  int color_count=0;
  int coarsest_level=0;

  int idx_mdot=-1; //idx_mdot==-1 => do not collect auxiliary data.

  int tessellate=1;
  int operation_flag=OP_GATHER_MDOT;
   //calling from: NavierStokes::multiphase_project
  ColorSumALL(
     operation_flag, // =OP_GATHER_MDOT
     tessellate, //=1
     coarsest_level,
     color_count,
     TYPE_MF,COLOR_MF,
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

  if (verbose>0) {
   if (ParallelDescriptor::IOProcessor()) {
    std::cout << "END: color_variable, multiphase_project\n";
    print_project_option(project_option);
   }
  }

  int idx_velcell=-1;

  operation_flag=OP_UNEW_CELL_TO_MAC;
  Real beta=0.0;

  if ((project_option==SOLVETYPE_PRES)||
      (project_option==SOLVETYPE_INITPROJ)) {

   if (nsolve!=1)
    amrex::Error("nsolve invalid");

   if (project_option==SOLVETYPE_PRES) {
    operation_flag=OP_UNEW_USOL_MAC_TO_MAC;
   } else if (project_option==SOLVETYPE_INITPROJ) {
    operation_flag=OP_UNEW_CELL_TO_MAC;
   } else 
    amrex::Error("project_option invalid47");

   increment_face_velocityALL(
    operation_flag,
    project_option,
    idx_velcell,beta,blobdata); 

   for (int ilev=finest_level;ilev>=level;ilev--) {
    NavierStokes& ns_level=getLevel(ilev);

    for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
     MultiFab* macvel=ns_level.getStateMAC(0,dir,cur_time_slab); 
     MultiFab::Copy(
      *ns_level.localMF[MAC_TEMP_MF+dir],
      *macvel,0,0,nsolve,0);
     MultiFab::Copy(
      *ns_level.localMF[UMAC_MF+dir],
      *macvel,0,0,nsolve,0);

     if (1==0) {
      int gridno=0;
      const Box& fabgrid = grids[gridno];
      const int* fablo=fabgrid.loVect();
      const int* fabhi=fabgrid.hiVect();
      const Real* xlo = grid_loc[gridno].lo();
      int interior_only=1;
      FArrayBox& macfab=(*macvel)[0];
      int scomp_debug=0;
      int ncomp_debug=nsolve;
      std::cout << "WARNING: this velocity is scaled\n";
      tecplot_debug(macfab,xlo,fablo,fabhi,coarse_dx,dir,0,
	     scomp_debug,ncomp_debug,interior_only);
     }
     delete macvel;
    }  // dir=0..sdim-1
   } // ilev=finest_level ... level

  } else
   amrex::Error("project_option invalid 10208");

  delete_array(TYPE_MF);
  delete_array(COLOR_MF);

 } else if (project_option_FSI_rigid(&project_option)==0) {
  // do nothing
 } else {
  amrex::Error("project_option_FSI_rigid invalid48");
 } 

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);

  if (project_option==SOLVETYPE_INITPROJ) { 
   // do nothing
  } else if (project_option_projection(project_option)==1) {

   if (project_option==SOLVETYPE_PRES) {
    // do nothing
   } else
    amrex::Error("project_option invalid 10233");

   // updates CELL_SOUND_MF, DIFFUSIONRHS, and S_new. 
   //  State_Type is updated by solver if project_option==SOLVETYPE_PRES. 
   // 
   //  NavierStokes::init_advective_pressure declared in NavierStokes2.cpp
   ns_level.init_advective_pressure(project_option); 

  } else if (project_option_projection(project_option)==0) {
   // do nothing 
  } else
   amrex::Error("project_option_projection invalid 49");

   // currently in: multiphase_project
   // calls fort_buidfacewt
   // fort_buildfacewt updates static variables min_face_wt and max_face_wt
   // max_face_wt[0][DD_COMP_FACE_WT] has max 
   //    of (1/rho) or (visc_coef*mu) or (k) or (D)
  int face_weight_op=SUB_OP_FOR_MAIN;
  ns_level.allocate_FACE_WEIGHT(nsolve,project_option,face_weight_op);

  ns_level.allocate_pressure_work_vars(nsolve,project_option);

   // 1. allocates and initializes ONES_MF, ONES_GROW_MF
   // 2. allocates, and sets POLDHOLD_MF=S^adv - S^init, 
   // 3. allocates, and sets OUTER_ITER_PRESSURE_MF=S^init,
   // 4. snew=S^init
   //
  ns_level.allocate_project_variables(nsolve,project_option);
 }  // ilev=finest_level ... level

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  show_norm2_id(FACE_WEIGHT_MF+dir,10+dir);
 }
 show_norm2_id(OUTER_ITER_PRESSURE_MF,15);
 show_norm2_id(POLDHOLD_MF,16);

 for (int ilist=0;ilist<scomp.size();ilist++) 
  avgDownALL(state_index,scomp[ilist],ncomp[ilist],1);

 Real maxden=denconst[0];
 for (int im=1;im<num_materials;im++) {
  if (denconst[im]>maxden)
   maxden=denconst[im];
 }

 if (maxden>0.0) {

  Real problen_max=0.0;
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   Real problen_local=geom.ProbHi(dir)-geom.ProbLo(dir);
   if (problen_local>0.0) {
    if (problen_local>problen_max)
     problen_max=problen_local;
   } else
    amrex::Error("problen_local invalid");
  } // dir=0..sdim-1

  if (problen_max>0.0) {
   //do nothing
  } else
   amrex::Error("problen_max invalid");
 } else
  amrex::Error("maxden invalid");

 sanity_check_face_wt(project_option);

 int finest_total=0;

#if (profile_solver==1)
 bprof.stop();
#endif

 set_local_tolerances(project_option);

#if (profile_solver==1)
 bprof.start();
 bprof.stop();
#endif

#if (profile_solver==1)
 bprof.start();
#endif

 energyflag=SUB_OP_FOR_MAIN;
 int homflag_residual_correction_form=0; 

 if (project_option==SOLVETYPE_INITPROJ) {
  homflag_residual_correction_form=1; 
 } else if (project_option_is_valid(project_option)==1) {
  homflag_residual_correction_form=0; 
 } else
  amrex::Error("project_option invalid 10327");

 cpp_overridepbc(homflag_residual_correction_form,project_option);

   // STATE_FOR_RESID is an input to 
   //  NavierStokes::residual_correction_form
   
  // ngrow=1
 getState_localMF_listALL(STATE_FOR_RESID_MF,1,
   state_index,scomp,ncomp);

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  if (ns_level.localMF[STATE_FOR_RESID_MF]->nComp()!=nsolve)
   amrex::Error("ns_level.localMF[STATE_FOR_RESID_MF]->nComp()!=nsolve");

  ns_level.resize_levelset(2,LEVELPC_MF);
  ns_level.debug_ngrow(LEVELPC_MF,2,local_caller_string);
 } // ilev=finest_level ... level

 // allocate_maccoefALL is declared in: MacProj.cpp
 //
 // initializes diagsing,mask_div_residual,mask_residual,
 // ONES_MF,ONES_GROW_MF
 //
 //  i.e.
 //  
 // calls:fort_scalarcoeff,fort_mult_facewt, fort_dividedx, fort_nsgenerate
 // initializes arrays holding the diagonal, ONES_MF, ONES_GROW_MF.
 // Regularizes FACE_WEIGHT_MF if necessary.
 // create_hierarchy=-1,0,1
 int create_hierarchy=0;
 allocate_maccoefALL(project_option,nsolve,create_hierarchy,dt_slab);

 if (create_hierarchy==0) {

  int zero_diag_flag=1;
  TypeALL(TYPE_ONES_MF,type_ONES_flag,zero_diag_flag);
  color_variable(coarsest_ONES_level,COLOR_ONES_MF,TYPE_ONES_MF,
   &color_ONES_count, 
   type_ONES_flag,zero_diag_flag);

  ones_sum_global.resize(color_ONES_count);
  // for each given color, singular_patch_flag=
  //   0 if color is masked off 
  //   1 if color is not masked off, no compressible/internal dirichlet 
  //     regions, and not touching a Dirichlet condition wall.
  //   2 if color is not masked off, a compressible/internal dirichlet
  //     region exists or color is touching a Dirichlet condition wall.
  singular_patch_flag.resize(color_ONES_count);

   // rhsnew=rhs-alpha H
   // 0 =sum rhs - alpha sum H
   // alpha=sum rhs/sum H
   // in otherwords:
   // if v in the null space of A,
   // we want rhs dot v =0
   // one_sum_global = v dot v
  dot_productALL_ones_size(project_option);

  if (project_option_singular_possible(project_option)==1) {

   if (nsolve!=1)
    amrex::Error("nsolve invalid34");

  } else if (project_option_singular_possible(project_option)==0) {

   for (int icolor=0;icolor<color_ONES_count;icolor++) {
    if (singular_patch_flag[icolor]==0) {
     // do nothing
    } else if (singular_patch_flag[icolor]==2) {
     // do nothing
    } else 
     amrex::Error("invalid singular_patch_flag[icolor]");
   } // icolor=0..color_ONES_count-1

  } else
   amrex::Error("project_option invalid50");

 } else
  amrex::Error("create_hierarchy invalid");

 int change_flag=0;

 // at very beginning:
 // POLDHOLD_MF=S^adv-S^init
 // OUTER_ITER_PRESSURE_MF=S^init
 // snew=S^init
 // STATE_FOR_RESID=S^init
 //
 project_right_hand_side(POLDHOLD_MF,project_option,change_flag);
 project_right_hand_side(OUTER_ITER_PRESSURE_MF,project_option,change_flag);
 project_right_hand_side(STATE_FOR_RESID_MF,project_option,change_flag);

 if (change_flag==0) {
    // do nothing
 } else if (change_flag==1) {
  putState_localMF_listALL(STATE_FOR_RESID_MF,
    state_index,scomp,ncomp);
  delete_array(STATE_FOR_RESID_MF);
  getState_localMF_listALL(STATE_FOR_RESID_MF,1,
    state_index,scomp,ncomp);
 } else
  amrex::Error("change_flag invalid");

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);

   // residual_correction_form declared in NavierStokes3.cpp
   // UMAC_MF-=beta grad STATE_FOR_RESID
  ns_level.residual_correction_form(
   homflag_residual_correction_form,
   energyflag,
   project_option,nsolve);

  if (ilev<finest_level)
   ns_level.avgDownMac();  // interpolates UMAC_MF from ilev+1

    // mac_temp store the velocity that comes from residual_correction_form
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   ns_level.Copy_localMF(MAC_TEMP_MF+dir,UMAC_MF+dir,0,0,nsolve,0);
  }

 }  // ilev=finest_level ... level

 delete_array(STATE_FOR_RESID_MF);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  show_norm2_id(MAC_TEMP_MF+dir,5+dir);
 }

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "0=DD_COMP_FACE_WT 1=CC_COMP_FACE_WT 2=MERGE_COMP_FACE_WT\n";
   print_project_option(project_option);
   for (int iwt=0;iwt<NCOMP_FACE_WT;iwt++) {
    std::cout << "iwt= " << iwt << " min_face_wt " << 
      min_face_wt[0][iwt] << " max_face_wt " <<
      max_face_wt[0][iwt] << '\n';
   }
   std::cout << "mglib_max_ratio=" << mglib_max_ratio << 
	   " min_interior_coeff=" << min_interior_coeff << '\n';
  } //  (ParallelDescriptor::IOProcessor()) 
 } // verbose>0

 deallocate_maccoefALL(project_option);
   
 int meets_tol=0;

 set_local_tolerances(project_option);

 Real error0_max=0.0;
 Real error0=0.0;
 Real error_n=0.0;
 Real error_at_the_beginning=0.0;
 Real error_after_all_jacobi_sweeps=0.0;
    
 if (verbose>0) {
   if (ParallelDescriptor::IOProcessor()) {
    std::cout << "LEVEL PROJECT from level=" << level << " to level=" <<
      ' ' << finest_level << '\n';
    print_project_option(project_option);
    std::cout << "NSOLVE=" << nsolve << '\n';
   }
 }

 int vcycle;
 int krylov_subspace_num_outer_iterSOLVER=0;
 int outer_iter_done=0;

 int min_krylov_subspace_outer_iter=0;

  // initializes diagsing,mask_div_residual,mask_residual,
  // ONES_MF,ONES_GROW_MF
  // create_hierarchy=-1,0,1
 create_hierarchy=1;
 allocate_maccoefALL(project_option,nsolve,create_hierarchy,dt_slab);

   // this must be done after allocate_maccoef (stefan_solver_init relies on
   // inhomogeneous BCs)
   // set BCs to homogeneous for the outer_iter loop.
 cpp_overridepbc(1,project_option);

 int total_number_vcycles=0;

 krylov_subspace_num_outer_iterSOLVER=0;
 outer_iter_done=0;
 
 lev0_cycles_list.resize(0);

#if (profile_solver==1)
 bprof.stop();
#endif

   // error_n,abs_tol
 Vector< Array<Real,2> > outer_error_history;
 outer_error_history.resize(krylov_subspace_max_num_outer_iter+1);
 for (int ehist=0;ehist<outer_error_history.size();ehist++) {
   outer_error_history[ehist][0]=0.0;
   outer_error_history[ehist][1]=0.0;
 }
 Real outer_error=0.0;

 Real best_error=0.0;
 int best_iter=-1;

 while (outer_iter_done==0) {

#if (profile_solver==1)
    bprof.start();
#endif

     // outer_iter_pressure holds the accumulated pressure now.
     // ext_dir boundary conditions are homogeneous for now on.

    for (int ilev=level;ilev<=finest_level;ilev++) {
      NavierStokes& ns_level=getLevel(ilev);
      ns_level.zero_independent_variable(project_option,nsolve);
    }

    meets_tol=0;

     // UMAC should be populated with velocity to be projected.
     // for viscous first sweep: UMAC=-dt mu grad U0   poldhold=0
    for (int ilev=finest_level;ilev>=level;ilev--) {
     NavierStokes& ns_level=getLevel(ilev);
     if (ilev<finest_level) {
      ns_level.setVal_localMF(MAC_RHS_CRSE_MF,0.0,0,nsolve,0); 
      ns_level.averageRhs(MAC_RHS_CRSE_MF,nsolve,project_option);
      ns_level.avgDownMac();  // works on UMAC_MF
     }
      // mac_phi_crse=0
      // mac_rhs_crse=POLDHOLD * alpha -  (note: POLDHOLD=p^advect-p^init)
      //              vol div UMAC/dt + diffusionRHS
     ns_level.mac_project_rhs(project_option,MAC_PHI_CRSE_MF,
       MAC_RHS_CRSE_MF,nsolve);
    } // ilev=finest_level ... level

      // this assumes that RHS is a discretization to 
      // volume * div (u)/dt
      // if v in the nullspace(A) then
      // 1. U=U0 + c v
      // 2. U dot v = U0 dot v + c(v dot v) = 0
      // 3. c=-(U0 dot v)/(v dot v)
      // 4. U=U0-[(U0 dot v)/(v dot v)] v
    change_flag=0;
     // change_flag set to 1 if MAC_RHS_CRSE_MF is modified.
    project_right_hand_side(MAC_RHS_CRSE_MF,project_option,change_flag);

    if (initial_project_cycles<1)
     amrex::Error("must do at least 1 jacobi cycle");
    if (initial_viscosity_cycles<1)
     amrex::Error("must do at least 1 jacobi cycle");
    if (initial_thermal_cycles<1)
     amrex::Error("must do at least 1 jacobi cycle");

    int jacobi_cycles_count=initial_project_cycles;

    if (project_option_singular_possible(project_option)==1) {
      // do nothing
    } else if (project_option==SOLVETYPE_HEAT) {
     jacobi_cycles_count=initial_thermal_cycles;
    } else if ((project_option>=SOLVETYPE_SPEC)&&
               (project_option<SOLVETYPE_SPEC+num_species_var)) { 
     jacobi_cycles_count=initial_thermal_cycles;
    } else if (project_option==SOLVETYPE_VISC) { 
     jacobi_cycles_count=initial_viscosity_cycles;
    } else
     amrex::Error("project_option invalid52");

    int update_vel=1;//update error0 IF krylov_subspace_num_outer_iterSOLVER==0
    int call_adjust_tolerance=1;

      // NavierStokes::jacobi_cycles declared in NavierStokes3.cpp
      // calling from: multiphase project
    jacobi_cycles(
      call_adjust_tolerance,
      jacobi_cycles_count,
      update_vel, // =1
      project_option,
      MAC_RHS_CRSE_MF,
      MAC_PHI_CRSE_MF, // null space projected out.
      error_at_the_beginning, // error before any Jacobi iterations 
      error_after_all_jacobi_sweeps, //error after jacobi iter. 
                                     //regardless num_outer_iterSOLVER
      error0, // error after Jacobi iterations IF num_outer_iterSOLVER==0
      error0_max, //max error0 during Jacobi Iterations IF 
                  //(num_outer_iterSOLVER==0)&&(call_adjust_tol==1)
      krylov_subspace_num_outer_iterSOLVER,
      nsolve);

    error_n=error0; // error after Jacobi iter. IF num_outer_iterSOLVER==0

    if (verbose>0) {
     if (ParallelDescriptor::IOProcessor()) {
      std::cout << "AFTER jacobi_cycles error0, error_after= " << error0 << 
       ' ' << error_after_all_jacobi_sweeps << '\n';
     }
    }

     // alpha deltap - div beta grad deltap=
     //   -(1/dt)div U + alpha poldhold
     // 
     // alpha dp - div beta grad dp=
     //   -(1/dt)div (U+V) + alpha poldhold 
     //
     // UMAC=UMAC-beta grad mac_phi_crse
     // S_new=S_new+mac_phi_crse
     //
     // POLDHOLD=POLDHOLD-mac_phi_crse
     //
     // mac_phi_crse=0
     //
     // updatevelALL calls mac_update.
     // calling from: multiphase_project.
    updatevelALL(project_option,MAC_PHI_CRSE_MF,nsolve);

    double after_startup=0.0;
    if (project_timings==1) {
     after_startup=ParallelDescriptor::second();
     if (ParallelDescriptor::IOProcessor())
      std::cout << "project start time " 
       << after_startup-begin_project << '\n';
    }

    int cg_loop_max=1;
    if (initial_cg_cycles>0) {
     cg_loop_max=2;
    } else if (initial_cg_cycles==0) {
     // do nothing
    } else
     amrex::Error("initial_cg_cycles invalid");

      // error_n,abs_tol
    Vector< Array<Real,2> > error_history;

#if (profile_solver==1)
    bprof.stop();
#endif

    for (int cg_loop=0;cg_loop<cg_loop_max;cg_loop++) {

#if (profile_solver==1)
     bprof.start();
#endif

     allocate_array(0,nsolve,-1,CGRESID_MF);
     allocate_array(1,nsolve,-1,P_SOLN_MF);
     allocate_array(0,nsolve,-1,bicg_V1_MF);

     for (int ilev=finest_level;ilev>=level;ilev--) {
      NavierStokes& ns_level=getLevel(ilev);
      if (ilev<finest_level) {
            // get rid of uninit.
       ns_level.setVal_localMF(MAC_RHS_CRSE_MF,0.0,0,nsolve,0);
       ns_level.averageRhs(MAC_RHS_CRSE_MF,nsolve,project_option);
       ns_level.avgDownMac(); // works on UMAC_MF
      }
       // mac_phi_crse_mf=0.0
       // mac_rhs_crse=POLDHOLD * alpha - 
       //              vol div UMAC/dt + diffusionRHS
       // NavierStokes::mac_project_rhs is declared in MacProj.cpp
      ns_level.mac_project_rhs(project_option,MAC_PHI_CRSE_MF,
        MAC_RHS_CRSE_MF,nsolve);
     } // ilev=finest_level ... level

       // MAC_PHI_CRSE=0.0 (from above)
       // CGRESID=MAC_RHS_CRSE-( alpha*phi-div grad phi )
       // null space is filtered out.
     residALL(project_option,MAC_RHS_CRSE_MF,
      CGRESID_MF,MAC_PHI_CRSE_MF,nsolve);

     Real local_error_n=0.0;
     dot_productALL(project_option,CGRESID_MF,CGRESID_MF,
      	      local_error_n,nsolve);
     if (local_error_n>=0.0) {
      local_error_n=std::sqrt(local_error_n);
     } else
      amrex::Error("local_error_n invalid");

      // if cg_loop==0 then error_n=error0 
      // error0=error after Jacobi iter. IF num_outer_iterSOLVER==0
     check_outer_solver_convergence(
        error_n,error0,
        save_mac_abs_tol,
        meets_tol);

     if (cg_loop==0) {
      if (error_after_all_jacobi_sweeps<save_mac_abs_tol) {
       meets_tol=1;
      }
     } else if (cg_loop==1) {
      // do nothing
     } else
      amrex::Error("cg_loop invalid");

      // double check that residual still meets the criterion
     if (meets_tol==1) {

      if (verbose>0) {
       if (ParallelDescriptor::IOProcessor()) {
        std::cout << "meets_tol=1 at top of CG or BICGSTAB: error_n=" <<
          error_n << '\n';
        if (cg_loop==0)
         std::cout << "error_after_all_jacobi_sweeps (jacobi method)=" <<
          error_after_all_jacobi_sweeps << '\n';
       } // ioproc
      }  // verbose>0

      error_n=local_error_n;
      Real local_tol=1.1*save_mac_abs_tol;
      meets_tol=0;
      check_outer_solver_convergence(
        error_n,error0,
        local_tol,
        meets_tol);

      if (verbose>0) {
       if (ParallelDescriptor::IOProcessor()) {
        std::cout << "double checking meets_tol=1: error_n=" <<
          error_n << " meets_tol= " << meets_tol << '\n';
       } // ioproc
      }  // verbose>0
     } else if (meets_tol==0) {
      // do nothing
     } else
      amrex::Error("meets_tol invalid");

     int local_presmooth=global_presmooth;
     int local_postsmooth=global_postsmooth;
     Real rho0=1.0;
     Real rho1=1.0;
     Real alpha=1.0;
     Real beta=0.0;
     Real w0=1.0;
     Real w1=0.0;
     Real a1=0.0;
     Real a2=0.0;

     int vcycle_max=multilevel_maxcycle;
     if ((initial_cg_cycles==0)&&(cg_loop==0)) {
      // do nothing
     } else if ((initial_cg_cycles>0)&&(cg_loop==0)) {
      vcycle_max=initial_cg_cycles;
     } else if ((initial_cg_cycles>0)&&(cg_loop==1)) {
      // do nothing
     } else
      amrex::Error("initial_cg_cycles or cg_loop invalid");

     int restart_flag=0;

     Real dnorm=0.0;

       // variables initialized to 0.0
     allocate_array(1,nsolve,-1,Z_MF);
     allocate_array(0,nsolve,-1,P_MF);
     allocate_array(0,nsolve,-1,bicg_R0hat_MF);
     allocate_array(1,nsolve,-1,bicg_U0_MF);
     allocate_array(0,nsolve,-1,bicg_V0_MF);
     allocate_array(0,nsolve,-1,bicg_P1_MF);
     allocate_array(0,nsolve,-1,bicg_R1_MF);
     allocate_array(1,nsolve,-1,bicg_Y_MF);
     allocate_array(1,nsolve,-1,bicg_Hvec_MF);
     allocate_array(0,nsolve,-1,bicg_S_MF);
     allocate_array(0,nsolve,-1,bicg_T_MF);

     Copy_array(bicg_R0hat_MF,CGRESID_MF,0,0,nsolve,0);

      // MAC_PHI_CRSE(U1)=0.0
     zeroALL(1,nsolve,bicg_U0_MF);

     error_history.resize(vcycle_max+1);
     for (int ehist=0;ehist<error_history.size();ehist++) {
      error_history[ehist][0]=0.0;
      error_history[ehist][1]=0.0;
     }

#if (profile_solver==1)
     bprof.stop();
#endif

     int BICGSTAB_ACTIVE=0;

     if (enable_spectral==1) {

      BICGSTAB_ACTIVE=1;

     } else if (enable_spectral==0) {

      if (project_option==SOLVETYPE_VISC) { 

       if (uncoupled_viscosity==1) {
        BICGSTAB_ACTIVE=0;
       } else if (uncoupled_viscosity==0) {
        BICGSTAB_ACTIVE=1;
       } else
        amrex::Error("local_uncoupled_viscosity invalid");

      } else if (project_option_is_valid(project_option)==1) {

       BICGSTAB_ACTIVE=0;

      } else
       amrex::Error("project_option invalid 10832");

     } else
      amrex::Error("enable_spectral invalid 3");

     vcycle=-1;
     best_error=0.0;
     best_iter=-1;

     Krylov_checkpoint(vcycle,error_n,best_error,best_iter,
         MAC_PHI_CRSE_MF,UMACSTAR_MF,restart_flag);

     int multilevel_restart_count=0;

     for (vcycle=0;((vcycle<=vcycle_max)&&(meets_tol==0));vcycle++) {

#if (profile_solver==1)
      bprof.start();
#endif

        // CGRESID(R0) is the residual when using U0
        // MAC_PHI_CRSE(U1) and U0 are the same at this point.
      dot_productALL(project_option,CGRESID_MF,CGRESID_MF,error_n,nsolve);
      if (error_n>=0.0) {
       error_n=std::sqrt(error_n);
      } else
       amrex::Error("error_n invalid");

      adjust_tolerance(error_n,error0_max,project_option);

      error_history[vcycle][0]=error_n;
      error_history[vcycle][1]=save_mac_abs_tol;

      restart_flag=0;

      check_outer_solver_convergence(
	error_n,error0,
        save_mac_abs_tol,
        meets_tol);

      if (verbose>0) {
       if (ParallelDescriptor::IOProcessor()) {
        std::cout << "cg_loop,vcycle,E0,En " << cg_loop << ' ' <<
         vcycle << ' ' << error0 << ' ' << error_n << '\n';
       }
       std::fflush(NULL);
      } else if ((verbose==0)||(verbose==-1)) {
       // do nothing
      } else
       amrex::Error("verbose invalid");

#if (profile_solver==1)
      bprof.stop();
#endif

      if (meets_tol==0) {

#if (profile_solver==1)
       bprof.start();
#endif

       if (BICGSTAB_ACTIVE==0) { //MGPCG

        // Z=M^{-1}CGRESID
        // Z=project(Z)
        multiphase_preconditioner(
         project_option,project_timings,
         local_presmooth,local_postsmooth,
         Z_MF,CGRESID_MF,nsolve);

        // rho1=z dot r=MINV r dot r =r^T MINV^T r
        // if rho1<=0 then safe to assume convergence has been 
        // achieved or something very wrong.
        dot_productALL(project_option,Z_MF,CGRESID_MF,rho1,nsolve);

        if (rho1>=0.0) {
         if (rho0>0.0) {
          beta=rho1/rho0;
          // P_MF=0 initially or on restart.
          // P_MF=Z_MF + beta P_MF
          mf_combine(project_option,
           Z_MF,P_MF,beta,P_MF,nsolve);
          change_flag=0;
          project_right_hand_side(P_MF,project_option,change_flag);

	  Copy_array(P_SOLN_MF,P_MF,0,0,nsolve,0);

          // V1=A P_SOLN
          // 1. (begin)calls project_right_hand_side(P)
          // 2. (end)  calls project_right_hand_side(V1)
          applyALL(project_option,P_SOLN_MF,bicg_V1_MF,nsolve);

          Real pAp=0.0;
          dot_productALL(project_option,P_MF,bicg_V1_MF,pAp,nsolve);

          if (pAp>0.0) {
           alpha=rho1/pAp;

            // mac_phi_crse(U1)=mac_phi_crse+alpha P
           mf_combine(project_option,
             MAC_PHI_CRSE_MF,P_MF,alpha,MAC_PHI_CRSE_MF,nsolve);
           change_flag=0;
           project_right_hand_side(MAC_PHI_CRSE_MF,
             project_option,change_flag);
            // U0=MAC_PHI_CRSE(U1)
	   Copy_array(bicg_U0_MF,MAC_PHI_CRSE_MF,0,0,nsolve,1);

           // CGRESID=RHS-A mac_phi_crse(U1)
           // 1. (start) calls project_right_hand_side(MAC_PHI_CRSE_MF)
           // 2. (end)   calls project_right_hand_side(CGRESID_MF)
           residALL(project_option,MAC_RHS_CRSE_MF,CGRESID_MF,
            MAC_PHI_CRSE_MF,nsolve);

           Real PCG_error_n=0.0;
           dot_productALL(project_option,CGRESID_MF,CGRESID_MF,
             PCG_error_n,nsolve);
           if (PCG_error_n>=0.0) {
            PCG_error_n=std::sqrt(PCG_error_n);
           } else
            amrex::Error("PCG_error_n invalid");

           Krylov_checkpoint(vcycle,PCG_error_n,best_error,best_iter,
             MAC_PHI_CRSE_MF,UMACSTAR_MF,restart_flag);

          } else if (pAp==0.0) {
           meets_tol=1;

           // this case can happen due to round off error when
           // A is a singular (indefinite) matrix
          } else if (pAp<0.0) { 
           meets_tol=1;
          } else {
           std::cout << "pAp= " << pAp << '\n';
           amrex::Error("pAp invalid in main solver");
          }
         } else if (rho0==0.0) {
          meets_tol=1;

          // this case can happen due to round off error when
          // A is a singular (indefinite) matrix
         } else if (rho0<0.0) {
          meets_tol=1;
         } else
          amrex::Error("rho0 invalid");

         // this case can happen due to round off error when
         // A is a singular (indefinite) matrix
        } else if (rho1<0.0) {
         meets_tol=1;
        } else
         amrex::Error("rho1 invalid");

       } else if (BICGSTAB_ACTIVE==1) { //MG PRECOND BICGSTAB

         // rho1=R0hat dot CGRESID(R0)
         //  =(b-A x0) dot (b-A xn) =
         //  b^T b + x0^T A^T A xn - x0^T A^T b - b^T A xn
        dot_productALL(project_option,bicg_R0hat_MF,CGRESID_MF,rho1,nsolve);

        if (vcycle==0) { // R0hat=R when vcycle==0

         if (rho1>0.0) {
          Real sanity_error=std::sqrt(rho1);
          if (sanity_error>0.9*save_mac_abs_tol) {
           // do nothing
          } else
           amrex::Error("sanity_error invalid");
         } else {
          std::cout << "rho1=" << rho1 << " vcycle= " << vcycle << '\n';
          amrex::Error("rho1 invalid");
         }

        } else if (vcycle>0) {
         // check nothing
        } else
         amrex::Error("vcycle invalid");

        if (rho0<=0.0) {
         restart_flag=1;
        } else if (rho0>0.0) {
         // do nothing
        } else
         amrex::Error("rho0 failed Nav3");

        if (w0<=0.0) {
         restart_flag=1;
        } else if (w0>0.0) {
         // do nothing
        } else
         amrex::Error("w0 failed Nav3");

        if (rho1<=0.0) {
         restart_flag=1;
        } else if (rho1>0.0) {
         // do nothing
        } else
         amrex::Error("rho1 invalid mglib");

#if (profile_solver==1)
        bprof.stop();
#endif

        if (restart_flag==0) {

#if (profile_solver==1)
         bprof.start();
#endif

         beta=(rho1/rho0)*(alpha/w0);

	 if (beta>0.0) {
	  // do nothing
	 } else
	  amrex::Error("beta invalid Nav3");

         a1=1.0;
         a2=-w0;

          // P1=P - w0 V0
         mf_combine(project_option,
	  P_MF,bicg_V0_MF,a2,bicg_P1_MF,nsolve); 

	 change_flag=0;
         project_right_hand_side(bicg_P1_MF,project_option,change_flag);

          // P1=CGRESID(R0)+beta P1
         mf_combine(project_option,
          CGRESID_MF,bicg_P1_MF,beta,bicg_P1_MF,nsolve); 

	 change_flag=0;
         project_right_hand_side(bicg_P1_MF,project_option,change_flag);

#if (profile_solver==1)
         bprof.stop();
#endif

          // Y=M^{-1}P1
	  // Y=project(Y)
         multiphase_SHELL_preconditioner(
          project_option,project_timings,
          local_presmooth,local_postsmooth,
          bicg_Y_MF,bicg_P1_MF,nsolve);

#if (profile_solver==1)
         bprof.start();
#endif

          // V1=A Y
          // 1. (begin)calls project_right_hand_side(Y)
          // 2. (end)  calls project_right_hand_side(V1)
         applyALL(project_option,bicg_Y_MF,bicg_V1_MF,nsolve);

         // alpha=R0hat dot V1=R0hat dot A M^-1 P1=
	 //       R0hat dot A M^-1 (R0+beta(P-w0 V0))
         dot_productALL(project_option,bicg_R0hat_MF,bicg_V1_MF,alpha,nsolve);

	 if (alpha<=0.0) {
	  restart_flag=1;
	 } else if (alpha>0.0) {
   	  // do nothing
	 } else
	  amrex::Error("alpha failed Nav3");

#if (profile_solver==1)
         bprof.stop();
#endif

         if (restart_flag==0) {

#if (profile_solver==1)
          bprof.start();
#endif

          alpha=rho1/alpha;

          // Hvec=U0+alpha Y
          a1=1.0;
          a2=alpha;
          mf_combine(project_option,
	   bicg_U0_MF,bicg_Y_MF,a2,bicg_Hvec_MF,nsolve);

	  change_flag=0;
          project_right_hand_side(bicg_Hvec_MF,project_option,change_flag);
          // mac_phi_crse(U1)=Hvec
          Copy_array(MAC_PHI_CRSE_MF,bicg_Hvec_MF,0,0,nsolve,1);

	  change_flag=0;
          project_right_hand_side(MAC_PHI_CRSE_MF,project_option,change_flag);

          // R1=RHS-A mac_phi_crse(U1)
	  // 1. (start) calls project_right_hand_side(MAC_PHI_CRSE_MF)
	  // 2. (end)   calls project_right_hand_side(bicg_R1_MF)
          residALL(project_option,MAC_RHS_CRSE_MF,bicg_R1_MF,
            MAC_PHI_CRSE_MF,nsolve);

	  // dnorm=R1 dot R1
          dot_productALL(project_option,bicg_R1_MF,bicg_R1_MF,dnorm,nsolve);
	  if (dnorm>=0.0) {
           dnorm=std::sqrt(dnorm);
	  } else
           amrex::Error("dnorm invalid Nav3");

          Krylov_checkpoint(vcycle,dnorm,best_error,best_iter,
             MAC_PHI_CRSE_MF,UMACSTAR_MF,restart_flag);

#if (profile_solver==1)
          bprof.stop();
#endif

          if (dnorm>save_mac_abs_tol) {

#if (profile_solver==1)
           bprof.start();
#endif

           // S=CGRESID(R0)-alpha V1
           a1=1.0;
           a2=-alpha;
           mf_combine(project_option,
  	    CGRESID_MF,bicg_V1_MF,a2,bicg_S_MF,nsolve);

	   change_flag=0;
           project_right_hand_side(bicg_S_MF,project_option,change_flag);

#if (profile_solver==1)
           bprof.stop();
#endif

           // Z=M^{-1}S
	   // Z=project(Z)
           multiphase_SHELL_preconditioner(
            project_option,project_timings,
            local_presmooth,local_postsmooth,
            Z_MF,bicg_S_MF,nsolve);

#if (profile_solver==1)
           bprof.start();
#endif

           // T=A Z
           // 1. (begin)calls project_right_hand_side(Z)
           // 2. (end)  calls project_right_hand_side(T)
           applyALL(project_option,Z_MF,bicg_T_MF,nsolve);

	   // a1 = T dot S = AZ dot MZ >=0.0 if A and M are SPD
           dot_productALL(project_option,bicg_T_MF,bicg_S_MF,a1,nsolve);
	   // a2 = T dot T = AZ dot AZ >=0.0 
           dot_productALL(project_option,bicg_T_MF,bicg_T_MF,a2,nsolve);

	   if (a2>0.0) {
   	    // do nothing
	   } else if (a2==0.0) {
            meets_tol=1;
	   } else {
	    amrex::Error("a2 invalid Nav3");
	   }

	   if (a1>0.0) {
	    // do nothing
	   } else if (a1<=0.0) {
            meets_tol=1;
	   } else
	    amrex::Error("a1 invalid");

#if (profile_solver==1)
           bprof.stop();
#endif

	   if (meets_tol==0) {

            if (restart_flag==0) {

#if (profile_solver==1)
             bprof.start();
#endif

             w1=a1/a2;
             // mac_phi_crse(U1)=Hvec+w1 Z
             a1=1.0;
             a2=w1;
             mf_combine(project_option,
   	      bicg_Hvec_MF,Z_MF,a2,MAC_PHI_CRSE_MF,nsolve);
	     change_flag=0;
             project_right_hand_side(MAC_PHI_CRSE_MF,
                 project_option,change_flag);

#if (profile_solver==1)
             bprof.stop();
#endif

            } else if (restart_flag==1) {
   	     // do nothing
	    } else
	     amrex::Error("restart_flag invalid");

	   } else if (meets_tol==1) {
	    // do nothing
	   } else
	    amrex::Error("meets_tol invalid");

          } else if ((dnorm>=0.0)&&(dnorm<=save_mac_abs_tol)) {
	   // do nothing (dnorm=R1 dot R1)
          } else
	   amrex::Error("dnorm invalid");

#if (profile_solver==1)
          bprof.start();
#endif

          // R1=RHS-A mac_phi_crse(U1)
	  // 1. (start) calls project_right_hand_side(MAC_PHI_CRSE_MF)
	  // 2. (end)   calls project_right_hand_side(bicg_R1_MF)
          residALL(project_option,MAC_RHS_CRSE_MF,bicg_R1_MF,
            MAC_PHI_CRSE_MF,nsolve);

	  // dnorm=R1 dot R1
          dot_productALL(project_option,bicg_R1_MF,bicg_R1_MF,dnorm,nsolve);
	  if (dnorm>=0.0) {
           dnorm=std::sqrt(dnorm);
	  } else
	   amrex::Error("dnorm invalid Nav3");

          Krylov_checkpoint(vcycle,dnorm,best_error,best_iter,
             MAC_PHI_CRSE_MF,UMACSTAR_MF,restart_flag);

          w0=w1;
          // CGRESID(R0)=R1
          Copy_array(CGRESID_MF,bicg_R1_MF,0,0,nsolve,0);
          Copy_array(P_MF,bicg_P1_MF,0,0,nsolve,0);
          Copy_array(bicg_V0_MF,bicg_V1_MF,0,0,nsolve,0);
          // U0=MAC_PHI_CRSE(U1)
          Copy_array(bicg_U0_MF,MAC_PHI_CRSE_MF,0,0,nsolve,0);

#if (profile_solver==1)
          bprof.stop();
#endif

         } else if (restart_flag==1) {
          // do nothing
         } else
          amrex::Error("restart_flag invalid");

        } else if (restart_flag==1) {
         // do nothing
        } else
         amrex::Error("restart_flag invalid");

       } else 
        amrex::Error("BICGSTAB_ACTIVE invalid");

       rho0=rho1;

       multilevel_restart_count++;

       if (multilevel_restart_count>=multilevel_restart_period) {
        restart_flag=1;
       } else if (multilevel_restart_count<multilevel_restart_period) {
        // do nothing
       } else
        amrex::Error("multilevel_restart_count invalid");

       if (restart_flag==0) {
        // do nothing
       } else if (restart_flag==1) {

        multilevel_restart_count=0;

#if (profile_solver==1)
        bprof.start();
#endif

	if (verbose>0) {
         if (ParallelDescriptor::IOProcessor()) {
          std::cout << "WARNING:RESTARTING: vcycle: " 
	         << vcycle << '\n';
          std::cout << "RESTARTING: BICGSTAB_ACTIVE=" << 
           BICGSTAB_ACTIVE << '\n';
          std::cout << "RESTARTING: local_presmooth= " << 
           local_presmooth << '\n';
          std::cout << "RESTARTING: local_postsmooth= " << 
           local_postsmooth << '\n';
          std::cout << "RESTARTING: error_history[vcycle][0,1]= " << 
           error_history[vcycle][0] << ' ' <<
	   error_history[vcycle][1] << '\n';
	  std::cout << "RESTARTING: best_error, best_iter = " <<
           best_error << ' ' << best_iter << '\n';
	  print_project_option(project_option);
	  std::cout << "END SOLVER RESTARTING NOTIFICATION\n";
         }
	} else if (verbose==0) {
	 // do nothing
	} else
	 amrex::Error("verbose invalid");

        Copy_array(bicg_U0_MF,RESTART_MAC_PHI_CRSE_MF,0,0,nsolve,1);
        for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
         Copy_array(UMACSTAR_MF+dir,RESTART_UMACSTAR_MF+dir,0,0,nsolve,0);
        }

        // CGRESID(R0)=RHS-A U0
	// 1. (start) calls project_right_hand_side(bicg_U0_MF)
	// 2. (end)   calls project_right_hand_side(CGRESID_MF)
        residALL(project_option,MAC_RHS_CRSE_MF,CGRESID_MF,
          bicg_U0_MF,nsolve);
        // R0hat=CGRESID(R0)
        Copy_array(bicg_R0hat_MF,CGRESID_MF,0,0,nsolve,0);
        // MAC_PHI_CRSE(U1)=U0
        Copy_array(MAC_PHI_CRSE_MF,bicg_U0_MF,0,0,nsolve,1);
        rho0=1.0;
        rho1=1.0;
        alpha=1.0;
        beta=0.0;
        w0=1.0;
        w1=0.0;
        a1=0.0;
        a2=0.0;

	// dnorm=RESID dot RESID
        dot_productALL(project_option,CGRESID_MF,CGRESID_MF,dnorm,nsolve);
	if (dnorm>=0.0) {
         dnorm=std::sqrt(dnorm);
	} else
         amrex::Error("dnorm invalid Nav3");

	if (verbose>0) {
         if (ParallelDescriptor::IOProcessor()) {
          std::cout << "RESTARTING: dnorm after restart: " << dnorm << '\n';
         }
	} else if (verbose==0) {
	 // do nothing
	} else
	 amrex::Error("verbose invalid");

        zeroALL(0,nsolve,bicg_V0_MF);
        zeroALL(0,nsolve,P_MF);

#if (profile_solver==1)
        bprof.stop();
#endif

       } else
        amrex::Error("restart_flag invalid");

      } // meets_tol==0

      // top level: BiCGStab  
      // preconditioner: MG
      // smooth for multigrid: "local_presmooth" iterations of the
      // ILU smoother going down the V-cycle and "local_postsmooth" 
      // iterations going up the V-cycle.
      // See: Pei, Vahab, Sussman, Hussaini JSC Hierarchical Spectral 
      // element method for multiphase flows. 2020
      if (restart_flag==1) {
       local_presmooth++;
       local_postsmooth++;
      } else if (restart_flag==0) {
       // do nothing
      } else
       amrex::Error("restart_flag invalid");
       
      if ((local_presmooth>1000)||(local_postsmooth>1000)) {	
       std::cout << "local_presmooth overflow " << 
          local_presmooth << '\n';
       std::cout << "local_postsmooth overflow " << 
          local_postsmooth << '\n';
       print_project_option(project_option);
       for (int ehist=0;ehist<error_history.size();ehist++) {
        std::cout << "vcycle " << ehist << 
         " error_history[vcycle][0,1] " <<
         error_history[ehist][0] << ' ' <<
         error_history[ehist][1] << '\n';
       }
       for (int ehist=0;ehist<outer_error_history.size();ehist++) {
        std::cout << "outer_iter " << ehist << 
         " outer_error_history[vcycle][0,1] " <<
         outer_error_history[ehist][0] << ' ' <<
         outer_error_history[ehist][1] << '\n';
       }
       amrex::Error("abort:NavierStokes3.cpp,local_pre(post)smooth overflow");
      }

      total_number_vcycles++;

     } // vcycle=0..vcycle_max or meets_tol!=0

#if (profile_solver==1)
     bprof.start();
#endif

     delete_array(Z_MF);
     delete_array(P_MF);
     delete_array(bicg_R0hat_MF);
     delete_array(bicg_U0_MF);
     delete_array(bicg_V0_MF);
     delete_array(bicg_P1_MF);
     delete_array(bicg_R1_MF);
     delete_array(bicg_Y_MF);
     delete_array(bicg_Hvec_MF);
     delete_array(bicg_S_MF);
     delete_array(bicg_T_MF);

     delete_array(bicg_V1_MF);
     delete_array(CGRESID_MF);
     delete_array(P_SOLN_MF);

     // alpha deltap - div beta grad deltap=
     //   -(1/dt)div U + alpha poldhold
     // 
     // alpha dp - div beta grad dp=
     //   -(1/dt)div (U+V) + alpha poldhold 
     // 
     // UMAC=UMAC-grad(mac_phi_crse) 
     // S_new=S_new+mac_phi_crse 
     //
     // POLDHOLD=POLDHOLD-mac_phi_crse
     //
     // mac_phi_crse=0
     //
     // updatevelALL calls mac_update.
     // calling from: multiphase_project
     updatevelALL(project_option,MAC_PHI_CRSE_MF,nsolve);

#if (profile_solver==1)
     bprof.stop();
#endif

    } // cg_loop=0..cg_loop_max-1

#if (profile_solver==1)
    bprof.start();
#endif

    std::fflush(NULL);

    if (vcycle>multilevel_maxcycle) {
     std::cout << "WARNING: vcycle too big, multilevel_maxcycle=" <<
      multilevel_maxcycle << '\n';
     std::cout << "->MyProc()= " << ParallelDescriptor::MyProc() << "\n";
     print_project_option(project_option);
     std::cout << "->error_n= " << error_n << '\n';
     std::cout << "->cg_loop_max= " << cg_loop_max << '\n';
     std::cout << "->ERROR HISTORY " << cg_loop_max << '\n';
     for (int ehist=0;ehist<error_history.size();ehist++) {
      std::cout << "vcycle " << ehist << 
       " error_history[vcycle][0,1] " <<
       error_history[ehist][0] << ' ' <<
       error_history[ehist][1] << '\n';
     }
     for (int ehist=0;ehist<outer_error_history.size();ehist++) {
      std::cout << "outer_iter " << ehist << 
       " outer_error_history[vcycle][0,1] " <<
       outer_error_history[ehist][0] << ' ' <<
       outer_error_history[ehist][1] << '\n';
     }
     std::cout << " best_error, best_iter = " <<
           best_error << ' ' << best_iter << '\n';
    } else if ((vcycle>=0)&&(vcycle<=multilevel_maxcycle)) {
     // do nothing
    } else {
     amrex::Error("vcycle bust");
    }

    std::fflush(NULL);

    for (int ilist=0;ilist<scomp.size();ilist++) 
     avgDownALL(state_index,scomp[ilist],ncomp[ilist],1);

      // override_bc_to_homogeneous=1
      // call fort_overridebc
    cpp_overridepbc(1,project_option);

    for (int ilev=finest_level;ilev>=level;ilev--) {
     NavierStokes& ns_level=getLevel(ilev);

     for (int ilist=0;ilist<scomp.size();ilist++) 
      ns_level.avgDown(state_index,scomp[ilist],ncomp[ilist],1);

     ns_level.getState_localMF_list(
      PRESPC_MF,1,
      state_index,
      scomp,
      ncomp);
    } // ilev=finest_level ... level

    for (int ilev=finest_level;ilev>=level;ilev--) {
     NavierStokes& ns_level=getLevel(ilev);
     int homflag_outer_iter_pressure2=1;
     energyflag=SUB_OP_FOR_MAIN; // energyflag=SUB_OP_FOR_SDC=>GRADPEDGE=gradp
      // GRADPEDGE=-dt gradp/rho
     int simple_AMR_BC_flag=0;
     int simple_AMR_BC_flag_viscosity=0;
     ns_level.apply_pressure_grad(
      simple_AMR_BC_flag,
      simple_AMR_BC_flag_viscosity,
      homflag_outer_iter_pressure2,
      energyflag,
      GRADPEDGE_MF,
      PRESPC_MF,
      project_option,nsolve,
      dt_slab); //calling from multiphase_project

     if (ilev<finest_level) {
      int ncomp_edge=-1;
       //spectral_override=1 => order derived from "enable_spectral"
      ns_level.avgDownEdge_localMF(GRADPEDGE_MF,
       0,ncomp_edge,0,AMREX_SPACEDIM,
       SPECTRAL_ORDER_AVGDOWN,local_caller_string);
     }

      // UMAC=MAC_TEMP+GRADPEDGE=MAC_TEMP-dt gradp/rho
     ns_level.correct_velocity(project_option,
      UMAC_MF, MAC_TEMP_MF,GRADPEDGE_MF,nsolve);
    }  // ilev=finest_level ... level

     // outer_iter_pressure+=S_new
     // S_new=0
     // mac_temp=UMAC
    for (int ilev=finest_level;ilev>=level;ilev--) {
     NavierStokes& ns_level=getLevel(ilev);
     for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
      ns_level.Copy_localMF(MAC_TEMP_MF+dir,UMAC_MF+dir,0,0,nsolve,0);

     MultiFab* snew_mf=nullptr;
     if (state_index==State_Type) {
      snew_mf=ns_level.getState_list(1,scomp,ncomp,cur_time_slab);
     } else {
      snew_mf=nullptr;
      amrex::Error("state_index invalid");
     }

     if (snew_mf->nComp()!=nsolve)
      amrex::Error("snew_mf->nComp() invalid");

     MultiFab::Add(
      *ns_level.localMF[OUTER_ITER_PRESSURE_MF],*snew_mf,0,0,nsolve,0);

     delete snew_mf;

     ns_level.zero_independent_variable(project_option,nsolve);
    }  // ilev=finest_level ... level

    change_flag=0;
    project_right_hand_side(OUTER_ITER_PRESSURE_MF,project_option,
      	     change_flag);

    delete_array(PRESPC_MF);

       // variables initialized to 0.0
    allocate_array(1,nsolve,-1,OUTER_MAC_PHI_CRSE_MF);
    allocate_array(0,nsolve,-1,OUTER_RESID_MF);
    allocate_array(0,nsolve,-1,OUTER_MAC_RHS_CRSE_MF);

    for (int ilev=finest_level;ilev>=level;ilev--) {
     NavierStokes& ns_level=getLevel(ilev);
     if (ilev<finest_level) {
      ns_level.setVal_localMF(OUTER_MAC_RHS_CRSE_MF,0.0,0,nsolve,0); 
      ns_level.averageRhs(OUTER_MAC_RHS_CRSE_MF,nsolve,project_option);
      ns_level.avgDownMac(); // works on UMAC_MF
     }
      // OUTER_MAC_PHI_CRSE=0
      // OUTER_MAC_RHS=POLDHOLD * alpha - 
      //               vol div UMAC/dt + diffusionRHS
     ns_level.mac_project_rhs(project_option,OUTER_MAC_PHI_CRSE_MF,
       OUTER_MAC_RHS_CRSE_MF,nsolve);
    }  // ilev=finest_level ... level

     // OUTER_RESID=OUTER_RHS-
     //          ( alpha * phi + (alpha+da) * phi - div grad phi )
     // null space is filtered
    residALL(project_option,OUTER_MAC_RHS_CRSE_MF,
      OUTER_RESID_MF,OUTER_MAC_PHI_CRSE_MF,nsolve);
    outer_error=0.0;
    dot_productALL(project_option,
      OUTER_RESID_MF,OUTER_RESID_MF,outer_error,nsolve);
    if (outer_error>=0.0) {
     outer_error=std::sqrt(outer_error);
    } else
     amrex::Error("outer_error invalid");

    delete_array(OUTER_MAC_RHS_CRSE_MF);
    delete_array(OUTER_MAC_PHI_CRSE_MF);
    delete_array(OUTER_RESID_MF);

    Real outer_tol=100.0*save_mac_abs_tol;

    outer_error_history[krylov_subspace_num_outer_iterSOLVER][0]=outer_error;
    outer_error_history[krylov_subspace_num_outer_iterSOLVER][1]=outer_tol;

    krylov_subspace_num_outer_iterSOLVER++;

    if (verbose>0) {
     if (ParallelDescriptor::IOProcessor()) {
      print_project_option(project_option);
      std::cout << "krylov_subspace_num_outer_iterSOLVER,E " << 
       krylov_subspace_num_outer_iterSOLVER << ' ' << outer_error << '\n';
     }
    }

    outer_iter_done=0;

    check_outer_solver_convergence(
      outer_error,error0,
      outer_tol,
      outer_iter_done);

    if (krylov_subspace_num_outer_iterSOLVER>krylov_subspace_max_num_outer_iter)
     outer_iter_done=1;

     // We cannot do iterative refinement (see Burden and Faires -
     // iterative techniques in matrix algebra) in this case since:
     // RHSPROJ( div(U - dt grad P) ) <> 
     // RHSPROJ( RHSPROJ(div U)- dt div grad P ) 
    if (krylov_subspace_num_outer_iterSOLVER>1)
     outer_iter_done=1;

    if (krylov_subspace_num_outer_iterSOLVER<min_krylov_subspace_outer_iter)
     outer_iter_done=0;

#if (profile_solver==1)
    bprof.stop();
#endif

 }  //  while outer_iter_done==0

#if (profile_solver==1)
 bprof.start();
#endif

 if (number_vcycles_all_solver_calls.size()==
     number_solver_calls.size()) {

  int prev_size=number_solver_calls.size();
  int new_size=project_option+1;
  if (new_size<=prev_size) {
   // do nothing
  } else if (new_size>prev_size) {
   number_vcycles_all_solver_calls.resize(new_size,0);
   max_lev0_cycles_all_solver_calls.resize(new_size,0);
   median_lev0_cycles_all_solver_calls.resize(new_size,0);
   outer_error_all_solver_calls.resize(new_size,0.0);
   number_solver_calls.resize(new_size,0);
  } else
   amrex::Error("new_size or prev_size invalid");

  int local_max_lev0_cycles=0;
  int local_median_lev0_cycles=0;
  int current_list_size=lev0_cycles_list.size();
  if (current_list_size>0) {
   local_max_lev0_cycles=lev0_cycles_list[0];
   local_median_lev0_cycles=lev0_cycles_list[current_list_size/2];
  }

  number_solver_calls[project_option]++;
  max_lev0_cycles_all_solver_calls[project_option]+=
         local_max_lev0_cycles;
  median_lev0_cycles_all_solver_calls[project_option]+=
         local_median_lev0_cycles;
  number_vcycles_all_solver_calls[project_option]+=total_number_vcycles;
  outer_error_all_solver_calls[project_option]+=outer_error;

  if (ParallelDescriptor::IOProcessor()) {
   Real avg_iter=number_vcycles_all_solver_calls[project_option]/
                number_solver_calls[project_option];
   std::cout << "SOLVER STATISTICS  TIME = " << cur_time_slab << '\n';
   print_project_option(project_option);
   std::cout << "project_option= " << project_option <<
          " SDC_outer_sweeps= " << SDC_outer_sweeps <<
          " slab_step= " << slab_step << 
          " FSI_outer_sweeps= " << FSI_outer_sweeps << '\n';
   std::cout << "project_option= " << project_option <<
	  " error0= " << error0 << '\n';

   Real avg_outer_error=outer_error_all_solver_calls[project_option]/
                number_solver_calls[project_option];

   std::cout << "project_option= " << project_option << 
          " number calls= " << number_solver_calls[project_option] << 
           " solver avg outer_error= " << avg_outer_error << '\n';

   std::cout << "project_option= " << project_option << 
          " number calls= " << number_solver_calls[project_option] << 
           " local outer_error= " << outer_error << '\n';

   std::cout << "project_option= " << project_option << 
          " number calls= " << number_solver_calls[project_option] << 
           " solver avg iter= " << avg_iter << '\n';
   std::cout << "project_option= " << project_option <<
          " TIME= " << cur_time_slab << " Current iterations= " <<
          total_number_vcycles << '\n';

   avg_iter=max_lev0_cycles_all_solver_calls[project_option]/
           number_solver_calls[project_option];
   std::cout << "project_option= " << project_option << 
          " number calls= " << number_solver_calls[project_option] << 
           " precon avg max= " << avg_iter << '\n';

   avg_iter=median_lev0_cycles_all_solver_calls[project_option]/
           number_solver_calls[project_option];
   std::cout << "project_option= " << project_option << 
          " number calls= " << number_solver_calls[project_option] << 
           " precon avg median= " << avg_iter << '\n';

   std::cout << "project_option= " << project_option <<
          " TIME= " << cur_time_slab << " local_max_lev0_cycles= " <<
          local_max_lev0_cycles << '\n';
   std::cout << "project_option= " << project_option <<
          " TIME= " << cur_time_slab << " local_median_lev0_cycles= " <<
          local_median_lev0_cycles << '\n';
   std::cout << "project_option= " << project_option <<
          " TIME= " << cur_time_slab << " best_iter= " << best_iter << '\n';
   std::cout << "project_option= " << project_option <<
          " TIME= " << cur_time_slab << " best_error= " << best_error << '\n';
  }

 } else
  amrex::Error("number_solver_calls.size() invalid");

 deallocate_maccoefALL(project_option);

    // copy OUTER_ITER_PRESSURE to s_new
 putState_localMF_listALL(OUTER_ITER_PRESSURE_MF,state_index,
   scomp,ncomp);

 for (int ilist=0;ilist<scomp.size();ilist++) 
  avgDownALL(state_index,scomp[ilist],ncomp[ilist],1);

 int homflag_dual_time=0;

 if (project_option==SOLVETYPE_INITPROJ) {
  homflag_dual_time=1;
 } else if (project_option==SOLVETYPE_PRESEXTRAP) { 
  homflag_dual_time=0;
 } else if ((project_option==SOLVETYPE_PRES)|| 
            (project_option==SOLVETYPE_HEAT)) { 
  homflag_dual_time=0;
 } else if (project_option==SOLVETYPE_VISC) { 
  homflag_dual_time=0;
 } else if ((project_option>=SOLVETYPE_SPEC)&&
	    (project_option<SOLVETYPE_SPEC+num_species_var)) {
  homflag_dual_time=0;
 } else
  amrex::Error("project_option invalid 53");

 cpp_overridepbc(homflag_dual_time,project_option);

 for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   if (ilev<finest_level)
    ns_level.avgDownMac(); // interpolates UMAC_MF from ilev+1
 } // ilev=finest_level ... level


#if (profile_solver==1)
 bprof.stop();
#endif

#if (profile_solver==1)
 bprof.start();
#endif

  //SOLVETYPE_INITPROJ, 
  //SOLVETYPE_PRES
 if (project_option_projection(project_option)==1) {

  getState_localMF_listALL(
    PRESPC2_MF,1,
    state_index,
    scomp,
    ncomp);

  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);

   // presBILINEAR2 and get_new_data(Umac_type+dir) are inputs.
   // Update cell velocity and temperature (compressible materials).

   // (multiphase_project)
   // If update_energy=SUB_OP_THERMAL_DIVUP_OK, 
   // then temperature might be updated if compressible material.
   // Copies UMAC to Umac_new: "save_to_macvel_state(UMAC_MF)"
   int update_energy=SUB_OP_THERMAL_DIVUP_NULL;
   if (project_option==SOLVETYPE_PRES) {

    if (num_FSI_outer_sweeps==1) {
     // update temperature if compressible material
     if (FSI_outer_sweeps==0) {
      update_energy=SUB_OP_THERMAL_DIVUP_OK; 
     } else
      amrex::Error("FSI_outer_sweeps invalid");
    } else if ((num_FSI_outer_sweeps>1)&&
               (num_FSI_outer_sweeps<=num_materials)) {
     if ((FSI_outer_sweeps>=0)&&
         (FSI_outer_sweeps<num_FSI_outer_sweeps-1)) {
      update_energy=SUB_OP_THERMAL_DIVUP_NULL;
     } else if (FSI_outer_sweeps==num_FSI_outer_sweeps-1) {
      update_energy=SUB_OP_THERMAL_DIVUP_OK; 
     } else
      amrex::Error("FSI_outer_sweeps invalid");
    } else
     amrex::Error("num_FSI_outer_sweeps invalid");

   } else if (project_option==SOLVETYPE_INITPROJ) {

    update_energy=SUB_OP_THERMAL_DIVUP_NULL;

   } else 
    amrex::Error("project_option invalid 11824");

    // NavierStokes::init_divup_cell_vel_cell declared
    // in NavierStokes2.cpp.
    // At the end of init_divup_cell_vel_cell,
    // UMAC_MF is copied to UMAC_new. 
   ns_level.init_divup_cell_vel_cell(project_option,
    update_energy,PRESPC2_MF,UMAC_MF);

   if (project_option==SOLVETYPE_PRES) {

    int combine_flag=2; //combine if vfrac<VOFTOL
    int hflag=0;
     // combine_idx==-1 => update S_new  
     // combine_idx>=0  => update localMF[combine_idx]
    int combine_idx=-1;  
    int update_flux=0;
    int interface_cond_avail=0;

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

      // velocity and pressure
    ns_level.avgDown(State_Type,STATECOMP_VEL,
		STATE_NCOMP_VEL+STATE_NCOMP_PRES,1);
    ns_level.avgDown(State_Type,STATECOMP_STATES,
		    num_state_material*num_materials,1);

   } else if (project_option==SOLVETYPE_INITPROJ) {
    ns_level.avgDown(State_Type,STATECOMP_VEL,STATE_NCOMP_VEL,1);
   } else
    amrex::Error("project_option invalid 54");

  } // ilev=finest_level ... level

  delete_array(PRESPC2_MF);

  if (project_option_needs_scaling(project_option)==1) {

   unscale_variablesALL();

  } else if (project_option==SOLVETYPE_INITPROJ) {
   // do nothing
  } else
   amrex::Error("project_option invalid 11886");

  if (project_option==SOLVETYPE_PRES) {
 
   if ((num_FSI_outer_sweeps==1)||
       ((num_FSI_outer_sweeps>1)&&
        (FSI_outer_sweeps==num_FSI_outer_sweeps-1))) {

    //calculate the vorticity needed by "init_pressure_error_indicator"
    int do_alloc=1;
    int simple_AMR_BC_flag_viscosity=1;
    init_gradu_tensorALL(
     HOLD_VELOCITY_DATA_MF,//deleted at end of init_gradu_tensorALL(do_alloc=1)
     do_alloc,
     CELLTENSOR_MF,
     FACETENSOR_MF,
     simple_AMR_BC_flag_viscosity);

    for (int ilev=finest_level;ilev>=level;ilev--) {
     NavierStokes& ns_level=getLevel(ilev);
     ns_level.init_pressure_error_indicator();  
    }
    avgDownError_ALL();

    delete_array(FACETENSOR_MF);
    delete_array(CELLTENSOR_MF);

   } else if ((num_FSI_outer_sweeps>1)&&
              (FSI_outer_sweeps>=0)&&
              (FSI_outer_sweeps<num_FSI_outer_sweeps-1)) {
    //do nothing
   } else 
    amrex::Error("num_FSI_outer_sweeps bust");

  } else if (project_option==SOLVETYPE_INITPROJ) {
   // do nothing
  } else
   amrex::Error("project_option invalid 55");

 } else if (project_option_needs_scaling(project_option)==1) {

  unscale_variablesALL();

 } else if (project_option_needs_scaling(project_option)==0) {
  // do nothing
 } else 
  amrex::Error("project_option invalid multiphase project 19");

 override_enable_spectral(save_enable_spectral);

  //homflag==0
 cpp_overridepbc(0,project_option);

 if (project_option==SOLVETYPE_PRESEXTRAP) {

  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
    // in: MacProj.cpp (calls fort_restore_pres)
   ns_level.restore_active_pressure(PRESSURE_SAVE_MF);
    // spectral_override==1 => order derived from "enable_spectral"
    // average from ilev+1 to ilev
   ns_level.avgDown(State_Type,STATECOMP_PRES,STATE_NCOMP_PRES,
     SPECTRAL_ORDER_AVGDOWN); 
  }
  delete_array(PRESSURE_SAVE_MF);

 } else if (project_option==SOLVETYPE_INITPROJ) {

  Copy_array(GET_NEW_DATA_OFFSET+State_Type,PRESSURE_SAVE_MF,
	  0,STATECOMP_PRES,STATE_NCOMP_PRES,1);
  delete_array(PRESSURE_SAVE_MF);

 } else if (project_option_is_valid(project_option)==1) {
  // do not restore anything
 } else
  amrex::Error("project_option invalid 11963");

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.remove_project_variables(); //remove poldhold,ones,outer_iter
  ns_level.remove_pressure_work_vars(); //remove umacstar, gradpedge, pedge
  ns_level.remove_FACE_WEIGHT_vars(); //remove FACE_WEIGHT,OFF_DIAG_CHECK
 }

 delete_array(MAC_PHI_CRSE_MF);
 delete_array(RESTART_MAC_PHI_CRSE_MF);
 delete_array(MAC_RHS_CRSE_MF);

 if (verbose>0)
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "FINEST SMOOTH TOTAL " << finest_total << '\n'; 
  }

   // delete diffusionRHS
   // delete mac_temp
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.remove_UMAC_for_solver(project_option);
 }

 remove_MAC_velocityALL(UMAC_MF);

 if (project_option_is_valid(project_option)==1) {
  // do nothing
 } else {
  amrex::Error("project_option invalid11993");
 }

#if (profile_solver==1)
 bprof.stop();
#endif

 std::fflush(NULL);

 finalize_rest_fraction(local_caller_string);

}  // end subroutine multiphase_project


void NavierStokes::diffusion_heatingALL(
  int source_idx,int idx_heat) {

 std::string local_caller_string="diffusion_heatingALL";

 int finest_level=parent->finestLevel();

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.resize_levelset(2,LEVELPC_MF);
  ns_level.debug_ngrow(LEVELPC_MF,2,local_caller_string);
  ns_level.VOF_Recon_resize(1); //output:SLOPE_RECON_MF
  ns_level.debug_ngrow(SLOPE_RECON_MF,1,local_caller_string);
  ns_level.debug_ngrow(source_idx,1,local_caller_string);
  ns_level.debug_ngrow(idx_heat,0,local_caller_string);
  ns_level.debug_ngrow(FACE_VAR_MF,0,local_caller_string);
  ns_level.resize_metrics(1);
  ns_level.debug_ngrow(VOLUME_MF,1,local_caller_string);
 }

 if (NS_geometry_coord==COORDSYS_RZ) {
  if (AMREX_SPACEDIM!=2)
   amrex::Error("dimension bust");
 } else if (NS_geometry_coord==COORDSYS_CARTESIAN) {
  // do nothing
 } else if (NS_geometry_coord==COORDSYS_CYLINDRICAL) {
  // do nothing
 } else
  amrex::Error("NS_geometry_coord bust 61");

  // currently in: diffusion_heatingALL
 min_face_wt.resize(thread_class::nthreads);
 max_face_wt.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  min_face_wt[tid].resize(NCOMP_FACE_WT);
  max_face_wt[tid].resize(NCOMP_FACE_WT);
  for (int iwt=0;iwt<NCOMP_FACE_WT;iwt++) {
   min_face_wt[tid][iwt]=1.0e+20;
   max_face_wt[tid][iwt]=-1.0e+20;
  }
 } // tid

  // currently in: NavierStokes::diffusion_heatingALL
  
 int nsolve=AMREX_SPACEDIM;
 for (int ilev=finest_level;ilev>=level;ilev--) {
  int face_weight_op=SUB_OP_FOR_MAIN;
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.allocate_FACE_WEIGHT(nsolve,SOLVETYPE_VISC,face_weight_op);
 }

 sanity_check_face_wt(SOLVETYPE_VISC);

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
   // ::diffusion_heating declared in Diffusion.cpp
  ns_level.diffusion_heating(source_idx,idx_heat);
 } // ilev

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.remove_FACE_WEIGHT_vars();
 }

}  // subroutine diffusion_heatingALL

void NavierStokes::avgDownALL_refine_density() {

 int finest_level = parent->finestLevel();
 if (level!=0)
  amrex::Error("only call with level=0");

 if ((num_materials_compressible>=1)&&
     (num_materials_compressible<=num_materials)) {

  for (int im_comp=0;im_comp<num_materials_compressible;im_comp++) {

   for (int i=finest_level-1;i>=level;i--) {
    NavierStokes& ns_level=getLevel(i);
    ns_level.avgDown_refine_density(im_comp);
   }

  } //im_comp

 } else if (num_materials_compressible==0) {
  //do nothing
 } else
  amrex::Error("num_materials_compressible invalid");

} //end subroutine avgDownALL_refine_density

void NavierStokes::avgDownALL_TENSOR() {

 int finest_level=parent->finestLevel();

 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=num_materials)) {

  for (int scomp=0;scomp<NUM_CELL_ELASTIC_REFINE;
       scomp+=ENUM_NUM_REFINE_DENSITY_TYPE) {

   for (int i=finest_level-1;i>=level;i--) {
    NavierStokes& ns_level=getLevel(i);
    ns_level.avgDown_refine_tensor(scomp);
   }

  } //scomp

 } else
  amrex::Error("num_materials_viscoelastic invalid avgDownALL_TENSOR");

} // end subroutine avgDownALL_TENSOR

// VISCOELASTIC FORCE
// if viscoelastic_force_only==0, then the following command is given 
// prior to this routine:
//  SET_STOKES_MARK(REGISTER_MARK_MF);
void NavierStokes::vel_elastic_ALL(int viscoelastic_force_only) {

 int elastic_force_mac_grid=1;
 if (viscoelastic_force_only==1) {
  elastic_force_mac_grid=0;
 } else if (viscoelastic_force_only==0) {
  //do nothing
 } else
  amrex::Error("viscoelastic_force_only has invalid value");

 std::string local_caller_string="vel_elastic_ALL";

 int finest_level=parent->finestLevel();

 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=num_materials)) {

  for (int im=0;im<num_materials;im++) {
    if (ns_is_rigid(im)==0) {
     if ((elastic_time[im]>0.0)&&
         (elastic_viscosity[im]>0.0)) {

      if (store_elastic_data[im]==1) {

       int push_enable_spectral=enable_spectral;
       int elastic_enable_spectral=0;
       override_enable_spectral(elastic_enable_spectral);

        // note: tensor_advection_updateALL is called before veldiffuseALL.
        // VISCOTEN_MF initialized in NavierStokes::make_viscoelastic_tensor
	// fort_maketensor called from ::make_viscoelastic_tensor.
	// We are currently in vel_elastic_ALL
       make_viscoelastic_tensorALL(im);

       if (1==0) {
        writeSanityCheckData(
         "VISCOTEN", //root_string
         "VISCOTEN", //information_string
         local_caller_string,
         VISCOTEN_MF, //tower_mf_id
         localMF[VISCOTEN_MF]->nComp(),
         VISCOTEN_MF,
         -1, // State_Type==-1
         -1, // data_dir==-1
         parent->levelSteps(0)); 
       }

       if (localMF[VISCOTEN_MF]->nComp()==
           ENUM_NUM_TENSOR_TYPE_REFINE) {
        // do nothing
       } else
        amrex::Error("localMF[VISCOTEN_MF]->nComp() invalid");

       if (localMF[VISCOTEN_MF]->nGrow()==1) {
        // do nothing
       } else
        amrex::Error("localMF[VISCOTEN_MF]->nGrow() invalid");

       int is_rigid_CL_flag=0;
       int imp1=im+1;

       if (fort_is_ice_base(&FSI_flag[im],&imp1)==1) {
        is_rigid_CL_flag=1;
       } else if (fort_is_FSI_rigid_base(&FSI_flag[im],&imp1)==1) {
        is_rigid_CL_flag=1;
       } else if (fort_is_FSI_elastic_base(&FSI_flag[im],&imp1)==1) {
        is_rigid_CL_flag=1;
       } else if (fort_FSI_flag_valid_base(&FSI_flag[im],&imp1)==1) {
        //do nothing
       } else
        amrex::Error("FSI_flag[im] invalid");

       int im_cutoff=num_materials;
       if (FSI_outer_sweeps==0) {
        im_cutoff=num_materials;
       } else if ((FSI_outer_sweeps>0)&&
                  (FSI_outer_sweeps<num_FSI_outer_sweeps)) {
        im_cutoff=im_elastic_map[FSI_outer_sweeps-1]+1;
       } else
        amrex::Error("FSI_outer_sweeps invalid");

       if ((FSI_outer_sweeps==0)||
           (viscoelastic_force_only==1)||
           ((FSI_outer_sweeps>=1)&&
            (is_rigid_CL_flag==0))||
           ((FSI_outer_sweeps>=1)&&
            (is_rigid_CL_flag==1)&&
            (imp1>im_cutoff))) {

         // find divergence of the X,Y,Z variables.
	 // NavierStokes::CELL_GRID_ELASTIC_FORCE is declared in
	 //    NavierStokes2.cpp
	 // CELL_GRID_ELASTIC_FORCE -> fort_elastic_force ->
	 // tensor_Heaviside
        for (int ilev=finest_level;ilev>=level;ilev--) {
         NavierStokes& ns_level=getLevel(ilev);
         ns_level.CELL_GRID_ELASTIC_FORCE(im,elastic_force_mac_grid);
        }

       } else if ((FSI_outer_sweeps>0)&&
                  (viscoelastic_force_only==0)&&
                  (is_rigid_CL_flag==1)&&
                  (imp1<=im_cutoff)) {
        //do nothing
       } else
        amrex::Error("FSI_outer_sweeps or is_rigid_CL_flag invalid");

       delete_array(VISCOTEN_MF);

       override_enable_spectral(push_enable_spectral);

      } else
       amrex::Error("store_elastic_data invalid");
 
     } else if ((elastic_time[im]==0.0)||
  	        (elastic_viscosity[im]==0.0)) {
      // do nothing
     } else
      amrex::Error("elastic_time[im] or elastic_viscosity[im] invalid");

    } else if (ns_is_rigid(im)==1) {
     // do nothing
    } else
     amrex::Error("ns_is_rigid invalid");
  } // im=0..num_materials-1
   
  // spectral_override==1 => order derived from "enable_spectral"
  avgDownALL(State_Type,STATECOMP_VEL,STATE_NCOMP_VEL+STATE_NCOMP_PRES,
    SPECTRAL_ORDER_AVGDOWN);

  if (elastic_force_mac_grid==1) {

   make_MAC_velocity_consistentALL();

  } else if (elastic_force_mac_grid==0) {
   //do nothing
  } else
   amrex::Error("elastic_force_mac_grid invalid");

   //vel_elastic_ALL called from veldiffuseALL
  if (viscoelastic_force_only==0) { 

   if (localMF[REGISTER_MARK_MF]->nComp()<AMREX_SPACEDIM)
    amrex::Error("REGISTER_MARK_MF invalid ncomp");
   if (localMF[REGISTER_MARK_MF]->nGrow()<1)
    amrex::Error("REGISTER_MARK_MF invalid ngrow");

    // umacnew+=INTERP_TO_MAC(unew-register_mark)  
    // (OP_UMAC_PLUS_VISC_CELL_TO_MAC)
   if (elastic_force_mac_grid==1) {
    //do nothing
   } else if (elastic_force_mac_grid==0) {
    INCREMENT_REGISTERS_ALL(REGISTER_MARK_MF); 
   } else
    amrex::Error("elastic_force_mac_grid invalid");

    // register_mark=unew
   SET_STOKES_MARK(REGISTER_MARK_MF);

   //vel_elastic_ALL called from writeTECPLOT_File
  } else if (viscoelastic_force_only==1) {
   // do nothing
  } else
   amrex::Error("viscoelastic_force_only invalid");

 } else if (num_materials_viscoelastic==0) {
  // do nothing
 } else
  amrex::Error("num_materials_viscoelastic invalid");

} // end subroutine vel_elastic_ALL


void NavierStokes::Mass_Energy_Sources_SinksALL() {

 int finest_level=parent->finestLevel();
 int nthreads_parm=thread_class::nthreads;

 fort_init_regions_list(
   constant_density_all_time.dataPtr(),
   &nthreads_parm);

 for (int isweep=0;isweep<2;isweep++) {
  if (level!=0)
   amrex::Error("it is required that level=0 in sum_regions isweep loop");
  if (level<=finest_level) {
   for (int ilev = 0; ilev <= finest_level; ilev++) {

    NavierStokes& ns_level = getLevel(ilev);
    ns_level.SumRegions(isweep);

   }

  } else
   amrex::Error("level<=finest_level required");

   // multi threading synchronization  
  fort_reduce_sum_regions(&isweep);

  int cpp_number_regions=0;
  fort_get_number_regions(&cpp_number_regions);

  if (cpp_number_regions==0) {
   // do nothing
  } else if (cpp_number_regions>0) {

   Vector<Real> cpp_energy_per_kelvin(cpp_number_regions);
   Vector<Real> cpp_mass(cpp_number_regions);
   Vector<Real> cpp_energy(cpp_number_regions);
   Vector<Real> cpp_volume(cpp_number_regions);
   Vector<Real> cpp_volume_raster(cpp_number_regions);
   Vector<Real> cpp_mass_after(cpp_number_regions);
   Vector<Real> cpp_energy_after(cpp_number_regions);
   Vector<Real> cpp_volume_after(cpp_number_regions);

   fort_get_region_data(&isweep,
    cpp_energy_per_kelvin.dataPtr(),
    cpp_mass.dataPtr(),
    cpp_energy.dataPtr(),
    cpp_volume.dataPtr(),
    cpp_volume_raster.dataPtr(),
    cpp_mass_after.dataPtr(),
    cpp_energy_after.dataPtr(),
    cpp_volume_after.dataPtr());

   for (int iregions=0;iregions<cpp_number_regions;iregions++) {
    if (isweep==0) {
     ParallelDescriptor::ReduceRealSum(cpp_energy_per_kelvin[iregions]);
     ParallelDescriptor::ReduceRealSum(cpp_mass[iregions]);
     ParallelDescriptor::ReduceRealSum(cpp_energy[iregions]);
     ParallelDescriptor::ReduceRealSum(cpp_volume[iregions]);
     ParallelDescriptor::ReduceRealSum(cpp_volume_raster[iregions]);
    } else if (isweep==1) {
     ParallelDescriptor::ReduceRealSum(cpp_mass_after[iregions]);
     ParallelDescriptor::ReduceRealSum(cpp_energy_after[iregions]);
     ParallelDescriptor::ReduceRealSum(cpp_volume_after[iregions]);
    } else
     amrex::Error("isweep invalid");
   }

   fort_put_region_data(&isweep,
    cpp_energy_per_kelvin.dataPtr(),
    cpp_mass.dataPtr(),
    cpp_energy.dataPtr(),
    cpp_volume.dataPtr(),
    cpp_volume_raster.dataPtr(),
    cpp_mass_after.dataPtr(),
    cpp_energy_after.dataPtr(),
    cpp_volume_after.dataPtr());

  } else
   amrex::Error("cpp_number_regions invalid");

 } //isweep=0,1

 int ioproc=0;
 if (ParallelDescriptor::IOProcessor()) {
  ioproc=1;
 }

 fort_delete_regions_list(&ioproc);

} // end subroutine Mass_Energy_Sources_SinksALL()

void NavierStokes::veldiffuseALL() {

 std::string local_caller_string="veldiffuseALL";

#if (NS_profile_solver==1)
 BLProfiler bprof(local_caller_string);
#endif

 // AmrLevel.H, protected:
 // static DescriptorList desc_lst
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  BCRec simulation_bc;
  if (dir==0) {
   set_x_vel_bc_NS_setup(simulation_bc,viscosity_phys_bc);
  } else if (dir==1) {
   set_y_vel_bc_NS_setup(simulation_bc,viscosity_phys_bc);
  } else if ((dir==AMREX_SPACEDIM-1)&&(AMREX_SPACEDIM==3)) {
   set_z_vel_bc_NS_setup(simulation_bc,viscosity_phys_bc);
  } else
   amrex::Error("dir invalid: using viscosity_phys_bc");

  desc_lst.reset_bcrecs(State_Type,STATECOMP_VEL+dir,simulation_bc);
 } //dir=0 .. sdim-1

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid");

 if ((slab_step<0)||(slab_step>=ns_time_order))
  amrex::Error("slab_step invalid");

 int finest_level=parent->finestLevel();

 int nden=num_materials*num_state_material;

 avgDownALL(State_Type,STATECOMP_VEL,STATE_NCOMP_VEL+STATE_NCOMP_PRES,1);

 for (int im=0;im<num_materials;im++) {
   if (heatviscconst[im]>0.0) {
    //do nothing
   } else if (heatviscconst[im]==0.0) {
    //do nothing
   } else
    amrex::Error("heatviscconst invalid");
 } // im 

 for (int im=0;im<num_materials*num_species_var;im++) {
   if (speciesviscconst[im]>0.0) {
    //do nothing
   } else if (speciesviscconst[im]==0.0) { //speciesreactionrate
    //do nothing
   } else
    amrex::Error("speciesviscconst invalid");
 } // im=0..num_materials*num_species_var-1 

 if (FSI_outer_sweeps==0) {

  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);

   int hflag=0;

   ns_level.solid_temperature();  // if solid temperature is prescribed

   ns_level.level_species_reaction(local_caller_string);

   // MEHDI VAHAB HEAT SOURCE
   // NavierStokes.cpp: void NavierStokes::make_heat_source()
   // make_heat_source calls GODUNOV_3D.F90::fort_heatsource which
   // calls PROB.F90::get_local_heat_source
   // The same temperature increment is added to all of the materials.
   ns_level.make_heat_source();  // updates S_new

   ns_level.getStateDen_localMF(save_state_MF,1,cur_time_slab);

   if (ns_level.localMF[save_state_MF]->nComp()==
       num_materials*num_state_material) {
    //do nothing
   } else
    amrex::Error("ns_level.localMF[save_state_MF].nComp() invalid");

   //FVM->GFM if phase change
   //FVM->mass weighted average if no phase change. 
  int combine_flag=0;  
   // combine_idx==-1 => update S_new  
   // combine_idx>=0  => update localMF[combine_idx]
  int combine_idx=-1; 
  int update_flux=0;
  int interface_cond_avail=1;
 
   //FVM->GFM if phase change
   //FVM->mass weighted average if no phase change. 
   combine_flag=0; 

   ns_level.combine_state_variable(
    SOLVETYPE_HEAT,
    combine_idx,
    combine_flag,
    hflag,
    update_flux,
    interface_cond_avail); 

   for (int ns=0;ns<num_species_var;ns++) {

    //FVM->GFM if phase change
    //FVM->mass weighted average if no phase change. 
    combine_flag=0; 

    ns_level.combine_state_variable(
     SOLVETYPE_SPEC+ns,
     combine_idx,
     combine_flag,
     hflag,
     update_flux,
     interface_cond_avail); 
   }  // ns=0;ns<num_species_var;ns++

  }  // ilev=finest_level ... level

  avgDownALL(State_Type,STATECOMP_STATES,nden,1);

 } else if ((FSI_outer_sweeps>0)&&
            (FSI_outer_sweeps<num_FSI_outer_sweeps)) {
  //do nothing
 } else
  amrex::Error("FSI_outer_sweeps invalid");


 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.debug_ngrow(FACE_VAR_MF,0,local_caller_string);
  ns_level.resize_metrics(1);
  ns_level.debug_ngrow(VOLUME_MF,1,local_caller_string);
 }

  // in: veldiffuseALL
  //  allocate scratch variables (including CONSERVE_FLUXES_MF)
  //
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.prepare_viscous_solver();
 }

 Vector<int> scomp;
 Vector<int> ncomp;
 int ncomp_check;
 int state_index;

   //num_materials_combine=1
 get_mm_scomp_solver(
  1,
  SOLVETYPE_HEAT,
  state_index,
  scomp,
  ncomp,
  ncomp_check);

 int nsolve_thermal=1;
 if (ncomp_check!=nsolve_thermal)
  amrex::Error("ncomp_check invalid");

 for (int im=0;im<num_materials;im++) {
  if (override_density[im]==0) { // Drho/DT=-divu rho
   // check nothing
  } else if (override_density[im]==1) { // rho=rho(T,Y)
   // check nothing
   
   //DrhoDT has units of 1/(Degrees Kelvin)
   // Du/Dt=-grad (p-rho0 g dot z)/rho0 - g DrhoDT (T-T0)
  } else if (override_density[im]==2) { 
   // check nothing
  } else
   amrex::Error("override_density invalid");  
 } // im=0..num_materials-1

   //ngrow=1
 getState_localMF_listALL(
   BOUSSINESQ_TEMP_MF,1,
   state_index,
   scomp,
   ncomp);

  // register_mark=unew (1 ghost)
 SET_STOKES_MARK(REGISTER_MARK_MF);

 show_norm2_id(REGISTER_MARK_MF,1);

 // -dt * |g| * beta * (T-T0)  beta<0
 // dt * beta * (T-T0) * omega^2 * r
 // u+=dt * (v^2/r +2 omega v)
 // v+=dt * (-uv/r -2 omega u)
 //
 // uncoupled_viscosity==0:
 // u=u-dt * (1/r) * mu * (3 v_t/r + 2 u/r)/rho
 // v=v+dt * (1/r) * mu * (3 u_t/r - v/r) 
 //
 // uncoupled_viscosity==1:
 // u=u-dt * (1/r) * mu * (2 v_t/r + u/r)/rho
 // v=v+dt * (1/r) * mu * (2 u_t/r - v/r) 
 //
 // U_t = F(U)  U=REGISTER_MARK_MF
 // 1. compute explicit terms: u^* = u^n + dt F1(u^n)
 // 2. compute implicit terms: u^** = u^* + dt F2(u^**)
 //

  // HOOP_FORCE_MARK_MF=(unp1-un)rho/dt
  // update_state==OP_HOOP_BOUSSINESQ_IMPLICIT:
  //  unp1(1)=unp1(1)/(one+param2*hoop_force_coef)-dt |g| beta (t-t0)
  //  unew=unp1
  // update_state==OP_HOOP_BOUSSINESQ_EXPLICIT:
  //  unp1(1)=unp1(1)-param2*hoop_force_coef*un(1)-dt |g| beta (t-t0)
  // The Boussinesq force is always explicit.
  //
  // force(D_DECL(i,j,k),dir)=(unp1(dir)-un(dir))/(inverseden*dt)
 int update_state=OP_HOOP_BOUSSINESQ_IMPLICIT;
 diffuse_hoopALL(REGISTER_MARK_MF,BOUSSINESQ_TEMP_MF,
   HOOP_FORCE_MARK_MF,update_state);

 show_norm2_id(REGISTER_MARK_MF,2);

   // spectral_override==1 => order derived from "enable_spectral"
 avgDownALL(State_Type,STATECOMP_VEL,STATE_NCOMP_VEL+STATE_NCOMP_PRES,
    SPECTRAL_ORDER_AVGDOWN);

  // umacnew+=INTERP_TO_MAC(unew-register_mark)  
  // (OP_UMAC_PLUS_VISC_CELL_TO_MAC)
 INCREMENT_REGISTERS_ALL(REGISTER_MARK_MF); 

 avgDownALL(State_Type,STATECOMP_STATES,nden,1);

  // register_mark=unew (1 ghost cell)
 SET_STOKES_MARK(REGISTER_MARK_MF);

 user_defined_momentum_forceALL(
   REGISTER_MARK_MF, //velocity
   BOUSSINESQ_TEMP_MF); //temperature

   // spectral_override==1 => order derived from "enable_spectral"
 avgDownALL(State_Type,STATECOMP_VEL,STATE_NCOMP_VEL+STATE_NCOMP_PRES,
    SPECTRAL_ORDER_AVGDOWN);

  // umacnew+=INTERP_TO_MAC(unew-register_mark)  
  // (OP_UMAC_PLUS_VISC_CELL_TO_MAC)
 INCREMENT_REGISTERS_ALL(REGISTER_MARK_MF); 

  // register_mark=unew
 SET_STOKES_MARK(REGISTER_MARK_MF);

// -----------veldiffuseALL: viscosity -----------------------------


 if (num_FSI_outer_sweeps-1==FSI_outer_sweeps) {

  if ((SDC_outer_sweeps>0)&&
      (SDC_outer_sweeps<ns_time_order)&&
      (divu_outer_sweeps+1==num_divu_outer_sweeps)) {

   if (ns_time_order>=2) {

    // fort_updatesemforce declared in GODUNOV_3D.F90:
    // HOFAB=-div(2 mu D) - HOOP_FORCE_MARK_MF 
    // (update_state=OP_HOOP_BOUSSINESQ_EXPLICIT at end of 
    //  NavierStokes::do_the_advance)
    // unew=unew-(1/rho)(int (HOFAB) - dt (LOFAB))
    if (enable_spectral==1) {
     for (int ilev=finest_level;ilev>=level;ilev--) {
      NavierStokes& ns_level=getLevel(ilev);
      // calls: fort_semdeltaforce in GODUNOV_3D.F90
      // does not look at: enable_spectral
      ns_level.make_SEM_delta_force(SOLVETYPE_VISC);
     }
    } else if (enable_spectral==0) {
     // do nothing
    } else
     amrex::Error("enable_spectral invalid");

    // spectral_override==1 => order derived from "enable_spectral"
    avgDownALL(State_Type,STATECOMP_VEL,STATE_NCOMP_VEL+STATE_NCOMP_PRES,
       SPECTRAL_ORDER_AVGDOWN);

    // umacnew+=INTERP_TO_MAC(unew-register_mark)
    // (OP_UMAC_PLUS_VISC_CELL_TO_MAC)
    INCREMENT_REGISTERS_ALL(REGISTER_MARK_MF); 

    // register_mark=unew
    SET_STOKES_MARK(REGISTER_MARK_MF);
   } else
    amrex::Error("ns_time_order invalid");

  } else if (SDC_outer_sweeps==0) {
   // do nothing
  } else if ((divu_outer_sweeps>=0)&&
             (divu_outer_sweeps+1<num_divu_outer_sweeps)) {
   // do nothing
  } else 
   amrex::Error("SDC_outer_sweeps or divu_outer_sweeps invalid");

 } else if ((FSI_outer_sweeps>=0)&&
            (FSI_outer_sweeps<num_FSI_outer_sweeps-1)) {
  //do nothing
 } else
  amrex::Error("FSI_outer_sweeps invalid");

 SET_STOKES_MARK(REGISTER_MARK_MF); //register_mark=unew

 show_norm2_id(REGISTER_MARK_MF,4);

  //multigrid precond. BiCGStab viscosity
 multiphase_project(SOLVETYPE_VISC); 
  
  //VISCHEAT_SOURCE_MF is a parameter for:
  //  update_SEM_forcesALL(SOLVETYPE_VISC,VISCHEAT_SOURCE_MF,...)
  //  init_gradu_tensorALL(VISCHEAT_SOURCE_MF,...)
  //  diffusion_heatingALL(VISCHEAT_SOURCE_MF,VISCHEAT_MF);
  //
 SET_STOKES_MARK(VISCHEAT_SOURCE_MF);

  // spectral_override==1 => order derived from "enable_spectral"
 avgDownALL(State_Type,STATECOMP_VEL,STATE_NCOMP_VEL+STATE_NCOMP_PRES,
    SPECTRAL_ORDER_AVGDOWN);

  // umacnew+=INTERP_TO_MAC(unew-register_mark)
  // (OP_UMAC_PLUS_VISC_CELL_TO_MAC)
 INCREMENT_REGISTERS_ALL(REGISTER_MARK_MF); 

  // spectral_override==1 => not always low order
 avgDownALL(State_Type,STATECOMP_STATES,nden,SPECTRAL_ORDER_AVGDOWN);

 SET_STOKES_MARK(REGISTER_MARK_MF); //register_mark=unew

// ---------------- end viscosity ---------------------

 int viscoelastic_force_only=0;
 vel_elastic_ALL(viscoelastic_force_only);

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.make_marangoni_force();
 } // ilev=finest_level ... level

  // spectral_override==1 => order derived from "enable_spectral"
 avgDownALL(State_Type,STATECOMP_VEL,STATE_NCOMP_VEL+STATE_NCOMP_PRES,
		 SPECTRAL_ORDER_AVGDOWN);
  // umacnew+=INTERP_TO_MAC(unew-register_mark)
  // (OP_UMAC_PLUS_VISC_CELL_TO_MAC)
 INCREMENT_REGISTERS_ALL(REGISTER_MARK_MF); 


 if (FSI_outer_sweeps==0) {

  // ---------------- begin thermal diffusion ---------------------

   // HOFAB=-div k grad T 
   // Tnew=Tnew-(1/(rho cv))(int (HOFAB) - dt (LOFAB))
   // call to fort_semdeltaforce with dcomp=slab_step*NSTATE_SDC+SEMDIFFUSE_T
  if ((SDC_outer_sweeps>0)&&
      (SDC_outer_sweeps<ns_time_order)&&
      (divu_outer_sweeps+1==num_divu_outer_sweeps)) {

   if (ns_time_order>=2) {

    if (enable_spectral==1) {

     for (int ilev=finest_level;ilev>=level;ilev--) {
      NavierStokes& ns_level=getLevel(ilev);
      // calls: fort_semdeltaforce in GODUNOV_3D.F90
      ns_level.make_SEM_delta_force(SOLVETYPE_HEAT); 
     }

    } else if (enable_spectral==0) {

     // do nothing
     
    } else
     amrex::Error("enable_spectral invalid");

    avgDownALL(State_Type,STATECOMP_STATES,nden,1);
   } else
    amrex::Error("ns_time_order invalid");

  } else if (SDC_outer_sweeps==0) {
   // do nothing
  } else if ((divu_outer_sweeps>=0)&&
             (divu_outer_sweeps+1<num_divu_outer_sweeps)) {
   // do nothing
  } else 
   amrex::Error("SDC_outer_sweeps or divu_outer_sweeps invalid");

    // why average down the density here?
  avgDownALL(State_Type,STATECOMP_STATES,nden,1);

  multiphase_project(SOLVETYPE_HEAT); 

  // --------------- end thermal diffusion -------------------



  // ---------------begin save stable (i.e. low order)
  //                thermal diffusion and viscous forces

  if ((ns_time_order>=2)&&
      (ns_time_order<=32)&&
      (divu_outer_sweeps+1==num_divu_outer_sweeps)) {

    //num_materials_combine=1
   get_mm_scomp_solver(
    1,
    SOLVETYPE_HEAT,
    state_index,
    scomp,
    ncomp,
    ncomp_check);

   if (ncomp_check!=nsolve_thermal)
    amrex::Error("ncomp_check invalid");

     //localMF[PRESPC2_MF] will hold the latest temperature from the implicit
     //backwards Euler system.
     //ngrow=1
   getState_localMF_listALL(
     PRESPC2_MF,  
     1,
     state_index,
     scomp,
     ncomp);
   
   int update_spectralF=0;
   int update_stableF=1;
    //currently in: NavierStokes::veldiffuseALL
    //LOfab=-div(k grad T)
    //NavierStokes::update_SEM_forcesALL declared in MacProj.cpp
    //NavierStokes::update_SEM_forces declared in MacProj.cpp
    //fort_updatesemforce declared in GODUNOV_3D.F90
   if (enable_spectral==1) {
    update_SEM_forcesALL(SOLVETYPE_HEAT,PRESPC2_MF,
     update_spectralF,update_stableF);
   } else if (enable_spectral==0) {
    // do nothing
   } else
    amrex::Error("enable_spectral invalid");

   delete_array(PRESPC2_MF);


    //currently in: NavierStokes::veldiffuseALL
    //LOfab=-div(2 mu D)-HOOP_FORCE_MARK_MF
    //NavierStokes::update_SEM_forcesALL declared in MacProj.cpp
    //NavierStokes::update_SEM_forces declared in MacProj.cpp
    //fort_updatesemforce declared in GODUNOV_3D.F90
   if (enable_spectral==1) {
    update_SEM_forcesALL(SOLVETYPE_VISC,VISCHEAT_SOURCE_MF,
     update_spectralF,update_stableF);
   } else if (enable_spectral==0) {
    // do nothing
   } else
    amrex::Error("enable_spectral invalid");

  } else if ((ns_time_order==1)||
             ((divu_outer_sweeps>=0)&&
              (divu_outer_sweeps+1<num_divu_outer_sweeps))) {
   // do nothing
  } else
   amrex::Error("ns_time_order or divu_outer_sweeps invalid");

  // ---------------end save stable thermal diffusion and viscous forces

  for (int species_comp=0;species_comp<num_species_var;species_comp++) {

    // (rho Y)_{t} + div(rho u Y) =div (rho D) grad Y
    // [rho_t + div(rho u)]Y+rho DY/DT=div(rho D) grad Y
    // DY/DT=div (rho D) grad Y/rho

    // why average down the density here?
   avgDownALL(State_Type,STATECOMP_STATES,nden,1);

   multiphase_project(SOLVETYPE_SPEC+species_comp); 

  } // species_comp=0..num_species_var-1

  avgDownALL(State_Type,STATECOMP_STATES,nden,1);

  for (int ilev=level;ilev<=finest_level;ilev++) {
   NavierStokes& ns_level=getLevel(ilev);

   int hflag=0;
   ns_level.solid_temperature();

   int combine_flag=1;  // GFM -> FVM
   int combine_idx=-1;  // update state variables
   int update_flux=0;
   int interface_cond_avail=1;

    //GFM->FVM if phase change
    //GFM copied to FVM if no phase change.
   combine_flag=1; 

   ns_level.combine_state_variable(
    SOLVETYPE_HEAT,
    combine_idx,
    combine_flag,
    hflag,
    update_flux,
    interface_cond_avail); 

   for (int ns=0;ns<num_species_var;ns++) {

    //GFM->FVM if phase change
    //GFM copied to FVM if no phase change.
    combine_flag=1;

    ns_level.combine_state_variable(
     SOLVETYPE_SPEC+ns,
     combine_idx,
     combine_flag,
     hflag,
     update_flux,
     interface_cond_avail); 
   }
  } // ilev=level ... finest_level


   // viscous heating is grad u : tau
  if (include_viscous_heating==1) {

   int do_alloc=0;

   int simple_AMR_BC_flag_viscosity=1;
   init_gradu_tensorALL(
     VISCHEAT_SOURCE_MF,
     do_alloc,
     CELLTENSOR_MF,
     FACETENSOR_MF,
     simple_AMR_BC_flag_viscosity);

   if ((num_materials_viscoelastic>=1)&&
       (num_materials_viscoelastic<=num_materials)) {

    for (int im=0;im<num_materials;im++) {
     if (ns_is_rigid(im)==0) {
      if ((elastic_time[im]>0.0)&&
          (elastic_viscosity[im]>0.0)) {
       // initializes VISCOTEN_MF
       // we are currently in "veldiffuseALL"
       make_viscoelastic_tensorALL(im);
       for (int ilev=finest_level;ilev>=level;ilev--) {
        NavierStokes& ns_level=getLevel(ilev);
         // VISCHEAT_MF is initialized to zero in ::prepare_viscous_solver().
         // VISCHEAT_MF is incremented with heating terms due to 
         // viscoelastic heating in this loop.
         // uses VISCOTEN_MF
        ns_level.make_viscoelastic_heating(im,VISCHEAT_MF);
       }  
       delete_array(VISCOTEN_MF);
      } else if ((elastic_time[im]==0.0)||
                 (elastic_viscosity[im]==0.0)) {
       // do nothing
      } else
       amrex::Error("elastic_time or elastic_viscosity invalid");

     } else if (ns_is_rigid(im)==1) {
      // do nothing
     } else
      amrex::Error("ns_is_rigid invalid");
    } // im=0..num_materials-1

   } else if (num_materials_viscoelastic==0) {
    // do nothing
   } else {
    amrex::Error("num_materials_viscoelastic invalid");
   }

    // viscosity heating
    // VISCHEAT_SOURCE is the velocity from the viscosity solver.
    // 1. localMF[CONSERVE_FLUXES_MF]=mu(grad U+grad U^T)
    // 2. localMF[VISCHEAT_MF]+=(dt/(rho cv)) mu(grad U+grad U^T) ddot grad U 
   diffusion_heatingALL(VISCHEAT_SOURCE_MF,VISCHEAT_MF);

    // add viscous (and viscoelastic) heating to T_m  m=1...M
   APPLY_VISCOUS_HEATINGALL(VISCHEAT_MF); // increment

      // overwrite T_m if phi_solid>0   m=1...M
   for (int ilev=finest_level;ilev>=level;ilev--) {
    NavierStokes& ns_level=getLevel(ilev);
    ns_level.solid_temperature();
   }

   delete_array(CELLTENSOR_MF);
   delete_array(FACETENSOR_MF);

  } else if (include_viscous_heating==0) {
    // do nothing
  } else
   amrex::Error("include_viscous_heating invalid");

 } else if ((FSI_outer_sweeps>0)&&
            (FSI_outer_sweeps<num_FSI_outer_sweeps)) {
  //do nothing
 } else
  amrex::Error("FSI_outer_sweeps invalid");

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.assimilate_state_data();
 }

 avgDownALL(State_Type,STATECOMP_VEL,STATE_NCOMP_VEL+STATE_NCOMP_PRES,1);
 avgDownALL(State_Type,STATECOMP_STATES,nden,1);

 // delete scratch variables (including CONSERVE_FLUXES_MF)
 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.exit_viscous_solver();
 }  // ilev
 delete_array(BOUSSINESQ_TEMP_MF);

 // AmrLevel.H, protected:
 // static DescriptorList desc_lst
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  BCRec simulation_bc;
  if (dir==0) {
   set_x_vel_bc_NS_setup(simulation_bc,phys_bc);
  } else if (dir==1) {
   set_y_vel_bc_NS_setup(simulation_bc,phys_bc);
  } else if ((dir==AMREX_SPACEDIM-1)&&(AMREX_SPACEDIM==3)) {
   set_z_vel_bc_NS_setup(simulation_bc,phys_bc);
  } else
   amrex::Error("dir invalid: using phys_bc");

  desc_lst.reset_bcrecs(State_Type,STATECOMP_VEL+dir,simulation_bc);
 } //dir=0 .. sdim-1

 if (FSI_outer_sweeps==0) {

  for (int ilev=level;ilev<=finest_level;ilev++) {
   NavierStokes& ns_level=getLevel(ilev);
   MultiFab& S_new=ns_level.get_new_data(State_Type,slab_step+1);
   MultiFab* save_S_new=ns_level.localMF[save_state_MF];
   for (int im=0;im<num_materials;im++) {

    if (heatviscconst[im]>0.0) {
     //do nothing
    } else if (heatviscconst[im]==0.0) {
     MultiFab::Copy(S_new,*save_S_new,
         num_state_material*im+ENUM_TEMPERATUREVAR,
         STATECOMP_STATES+num_state_material*im+ENUM_TEMPERATUREVAR,1,1);
    } else
     amrex::Error("heatviscconst invalid");

    for (int ns=0;ns<num_species_var;ns++) {

     if (speciesviscconst[ns*num_materials+im]>0.0) {
      //do nothing
     } else if (speciesviscconst[ns*num_materials+im]==0.0) {
      MultiFab::Copy(S_new,*save_S_new,
         num_state_material*im+num_state_base+ns,
         STATECOMP_STATES+num_state_material*im+num_state_base+ns,1,1);
     } else
      amrex::Error("heatviscconst invalid");

    } //ns=0 ... num_species_var-1

   } // im=0..num_materials-1
  } //ilev=level ... finest_level

  delete_array(save_state_MF);
 } else if ((FSI_outer_sweeps>0)&&
            (FSI_outer_sweeps<num_FSI_outer_sweeps)) {
  //do nothing
 } else
  amrex::Error("FSI_outer_sweeps invalid");


#if (NS_profile_solver==1)
 bprof.stop();
#endif


}   // end subroutine veldiffuseALL

void NavierStokes::PCINTERP_fill_bordersALL(int idx_MF,
  int ngrow,int scomp,int ncomp,int index,Vector<int> scompBC_map) {

 int finest_level=parent->finestLevel();
 if (level!=0)
  amrex::Error("level invalid PCINTERP_fill_bordersALL");

 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.PCINTERP_fill_borders(idx_MF,ngrow,scomp,ncomp,
    index,scompBC_map);
 } // ilev

} //end subroutine PCINTERP_fill_bordersALL


void NavierStokes::PCINTERP_fill_borders(int idx_MF,int ngrow,
  int scomp,int ncomp,int index,Vector<int> scompBC_map) {

 int ncompcheck=localMF[idx_MF]->nComp();
 int ngrowcheck=localMF[idx_MF]->nGrow();

 if (scomp+ncomp>ncompcheck)
  amrex::Error("ncomp too big");
 if (ngrowcheck<ngrow)
  amrex::Error("ngrow too big in PCINTERP_fill_borders");
 
 MultiFab* cmf=nullptr;

 if (level==0) {
  cmf=localMF[idx_MF];
 } else if (level>0) {
  NavierStokes& ns_coarse=getLevel(level-1);
  cmf=ns_coarse.localMF[idx_MF];
 } else {
  cmf=nullptr;
  amrex::Error("level invalid PCINTERP_fill_borders");
 }

  // uses desc_lstGHOST[index] instead of dest_lst[index]
  // AmrLevel::InterpBordersGHOST is declared in AmrLevel.cpp.
 InterpBordersGHOST(
    *cmf,
    *localMF[idx_MF],
    cur_time_slab,
    index,
    scomp,
    scompBC_map,
    ncomp,
    debug_fillpatch);   

} //end subroutine PCINTERP_fill_borders


void NavierStokes::PCINTERP_fill_coarse_patch(int idx_MF,
  int scomp,int ncomp,int index,Vector<int> scompBC_map) {

 if (level<=0)
  amrex::Error("level invalid: PCINTERP_fill_coarse_patch");

 int ncompcheck=localMF[idx_MF]->nComp();
 int ngrowcheck=localMF[idx_MF]->nGrow();

 if (scomp+ncomp>ncompcheck)
  amrex::Error("ncomp too big");
 if (ngrowcheck<0)
  amrex::Error("ngrowcheck invalid in PCINTERP_fill_coarse_patch");
 
 NavierStokes& ns_coarse=getLevel(level-1);
 MultiFab* cmf=ns_coarse.localMF[idx_MF];

 FillCoarsePatchGHOST(
    *cmf,
    *localMF[idx_MF],
    cur_time_slab,
    index,
    scomp,
    scompBC_map,
    ncomp,
    debug_fillpatch);   

} //end subroutine PCINTERP_fill_coarse_patch



void NavierStokes::GetStateFromLocalALL(int idx_MF,
  int ngrow,int scomp,int ncomp,int index,Vector<int> scompBC_map) {

 int finest_level=parent->finestLevel();
 if (level!=0)
  amrex::Error("level invalid GetStateFromLocalALL");

 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.GetStateFromLocal(idx_MF,ngrow,scomp,ncomp,
    index,scompBC_map);
 } // ilev

} //end subroutine GetStateFromLocalALL


void NavierStokes::GetStateFromLocal(int idx_MF,int ngrow,
  int scomp,int ncomp,int index,Vector<int> scompBC_map) {

 int ncompcheck=localMF[idx_MF]->nComp();
 int ngrowcheck=localMF[idx_MF]->nGrow();

 if (scomp+ncomp>ncompcheck)
  amrex::Error("ncomp too big");
 if (ngrowcheck<ngrow)
  amrex::Error("ngrow too big in GetStateFromLocal");
 
 MultiFab* cmf=nullptr;

 if (level==0) {
  cmf=localMF[idx_MF];
 } else if (level>0) {
  NavierStokes& ns_coarse=getLevel(level-1);
  cmf=ns_coarse.localMF[idx_MF];
 } else {
  cmf=nullptr;
  amrex::Error("level invalid GetStateFromLocal");
 }

  // uses desc_lst[index] 
 InterpBorders(
    *cmf,
    *localMF[idx_MF],
    cur_time_slab,
    index,
    scomp,
    scompBC_map,
    ncomp,
    debug_fillpatch);   

} //end subroutine GetStateFromLocal



void NavierStokes::delete_advect_vars() {

 delete_localMF(ADVECT_REGISTER_FACE_MF,AMREX_SPACEDIM);
 delete_localMF(ADVECT_REGISTER_MF,1);

} // subroutine delete_advect_vars()

void NavierStokes::delete_transport_vars() {

 delete_localMF(TRANSPORT_REGISTER_FACE_MF,AMREX_SPACEDIM);

} // subroutine delete_transport_vars()

void NavierStokes::exit_viscous_solver() {

 delete_localMF(CONSERVE_FLUXES_MF,AMREX_SPACEDIM);

 delete_localMF(REGISTER_MARK_MF,1);
 delete_localMF(HOOP_FORCE_MARK_MF,1);
 delete_localMF(VISCHEAT_MF,1);
 delete_localMF(VISCHEAT_SOURCE_MF,1);

}

// theta new_m = theta new_m + dtheta  m=1..num_materials  
// dtheta=viscous heating term
void NavierStokes::APPLY_VISCOUS_HEATINGALL(int source_mf) { // increment (du)

 int finest_level=parent->finestLevel();
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.APPLY_VISCOUS_HEATING(source_mf);
 }  // ilev=finest_level ... level

} // end subroutine APPLY_VISCOUS_HEATINGALL

void NavierStokes::APPLY_VISCOUS_HEATING(int source_mf) {

 int finest_level=parent->finestLevel();

 bool use_tiling=ns_tiling;

 int nsolve=1;

 if (localMF[source_mf]->nComp()!=nsolve)
  amrex::Error("diffuse_register invalid ncomp");
 if (localMF[source_mf]->nGrow()<0)
  amrex::Error("diffuse_register invalid ngrow");
 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

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
   amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nComp() bad");
 }

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 MultiFab& LS_new=get_new_data(LS_Type,slab_step+1);

 int nstate=STATE_NCOMP;
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");
 if (LS_new.nComp()!=num_materials*(AMREX_SPACEDIM+1))
  amrex::Error("LS_new.nComp()!=num_materials*(AMREX_SPACEDIM+1)");

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

  FArrayBox& solxfab=(*localMF[FSI_GHOST_MAC_MF])[mfi];
  FArrayBox& solyfab=(*localMF[FSI_GHOST_MAC_MF+1])[mfi];
  FArrayBox& solzfab=(*localMF[FSI_GHOST_MAC_MF+AMREX_SPACEDIM-1])[mfi];

  FArrayBox& snewfab=S_new[mfi];
  FArrayBox& lsfab=LS_new[mfi];
  FArrayBox& dufab=(*localMF[source_mf])[mfi];

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // in: GODUNOV_3D.F90
   // snew=TNEW+du
  fort_heatadvance(
    &level,
    &finest_level,
    &cur_time_slab,
    &nparts,
    &nparts_def,
    im_solid_map_ptr,
    &nsolve, 
    &nstate,
    xlo,dx,
    solxfab.dataPtr(),ARLIM(solxfab.loVect()),ARLIM(solxfab.hiVect()),
    solyfab.dataPtr(),ARLIM(solyfab.loVect()),ARLIM(solyfab.hiVect()),
    solzfab.dataPtr(),ARLIM(solzfab.loVect()),ARLIM(solzfab.hiVect()),
    snewfab.dataPtr(),ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()),
    lsfab.dataPtr(),ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
    dufab.dataPtr(),ARLIM(dufab.loVect()),ARLIM(dufab.hiVect()),
    tilelo,tilehi, 
    fablo,fabhi,&bfact);

 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_HEATADVANCE,"APPLY_VISCOUS_HEATING");

} // subroutine APPLY_VISCOUS_HEATING

//REGISTER_CURRENT_MF=unew-source_mf
//uface+=INTERP_TO_MAC(REGISTER_CURRENT_MF)
// (OP_UMAC_PLUS_VISC_CELL_TO_MAC)
void NavierStokes::INCREMENT_REGISTERS_ALL(int source_mf) {

 if (level!=0)
  amrex::Error("level invalid INCREMENT_REGISTERS_ALL");
 int finest_level=parent->finestLevel();

 if (source_mf==REGISTER_MARK_MF) {
  // do nothing
 } else
  amrex::Error("expecting source_mf==REGISTER_MARK_MF");

  // 1. allocate REGISTER_CURRENT_MF
  // 2. REGISTER_CURRENT_MF=unew-source_mf
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.INCREMENT_REGISTERS(source_mf);
 }

  // unew^f=unew^f+beta * diffuse_register^{c->f}
  // in: INCREMENT_REGISTERS_ALL
 int operation_flag=OP_UMAC_PLUS_VISC_CELL_TO_MAC;
 Real beta=1.0;
 Vector<blobclass> blobdata;

  // operation_flag==OP_UMAC_PLUS_VISC_CELL_TO_MAC 
  // unew^f=unew^f+beta * REGISTER_CURRENT_MF^{c->f}
 increment_face_velocityALL(
   operation_flag,
   SOLVETYPE_VISC,
   REGISTER_CURRENT_MF, //dummy parameter name: "idx_velcell"
   beta,
   blobdata);

 delete_array(REGISTER_CURRENT_MF);

} // end subroutine INCREMENT_REGISTERS_ALL

// 1. allocate REGISTER_CURRENT_MF
// 2. REGISTER_CURRENT_MF=(unew-source)
// (OP_UMAC_PLUS_VISC_CELL_TO_MAC uses this data)
void NavierStokes::INCREMENT_REGISTERS(int source_mf) {

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 if (source_mf==REGISTER_MARK_MF) {
  // do nothing
 } else
  amrex::Error("expecting source_mf==REGISTER_MARK_MF");

 int nsolve=AMREX_SPACEDIM;

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 int nstate=S_new.nComp();
 if (nstate!=STATE_NCOMP)
  amrex::Error("nstate invalid");

 new_localMF(REGISTER_CURRENT_MF,nsolve,1,-1);
 push_back_state_register(REGISTER_CURRENT_MF,cur_time_slab);

 MultiFab::Subtract(
  *localMF[REGISTER_CURRENT_MF],
  *localMF[source_mf],0,0,nsolve,1);

} // INCREMENT_REGISTERS


void NavierStokes::push_back_state_register(int idx_MF,Real time) {

 int nsolve=AMREX_SPACEDIM;

 int state_index;
 Vector<int> scomp;
 Vector<int> ncomp;
 int ncomp_check;
  //num_materials_combine=1
 get_mm_scomp_solver(
  1,
  SOLVETYPE_VISC,
  state_index,
  scomp,
  ncomp,
  ncomp_check);
 if (ncomp_check!=nsolve)
  amrex::Error("nsolve invalid 6613");

 if (state_index!=State_Type)
  amrex::Error("state_index invalid");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");
 if (localMF[idx_MF]->nComp()!=nsolve)
  amrex::Error("cell_register invalid ncomp");
 if (localMF[idx_MF]->nGrow()<1)
  amrex::Error("cell_register invalid ngrow");

 MultiFab* snew_mf;
 snew_mf=getState_list(1,scomp,ncomp,time);

 if (snew_mf->nComp()!=nsolve)
  amrex::Error("snew_mf->nComp() invalid");

 check_for_NAN(snew_mf);

 MultiFab::Copy(*localMF[idx_MF],*snew_mf,0,0,nsolve,1);
 delete snew_mf;

} // subroutine push_back_state_register

// stores the current velocity in localMF[idx_MF]
void NavierStokes::SET_STOKES_MARK(int idx_MF) {

 if (level!=0)
  amrex::Error("level invalid SET_STOKES_MARK");
 int finest_level=parent->finestLevel();

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.push_back_state_register(idx_MF,cur_time_slab);  
 }

}  // SET_STOKES_MARK


//allocates, saves, and extends FSI_CELL|MAC_VELOCITY_MF if sweeps>=0 and
//sweeps<num_FSI_outer_sweeps-1
//copies fluid velocity and deletes if sweeps>=1 and 
//sweeps<=num_FSI_outer_sweeps
void NavierStokes::manage_FSI_data() {

 int finest_level=parent->finestLevel();

 std::string local_caller_string="manage_FSI_data";

 if ((num_FSI_outer_sweeps>=2)&&
     (num_FSI_outer_sweeps<=num_materials)) {

  bool use_tiling=ns_tiling;

  if (num_state_base!=2)
   amrex::Error("num_state_base invalid");

  resize_levelset(3,LEVELPC_MF);
  debug_ngrow(LEVELPC_MF,3,local_caller_string);

  resize_maskfiner(1,MASKCOEF_MF);
  resize_mask_nbr(1);

  const Real* dx = geom.CellSize();

  const Box& domain = geom.Domain();
  const int* domlo = domain.loVect();
  const int* domhi = domain.hiVect();

  MultiFab& S_new=get_new_data(State_Type,slab_step+1);

  if ((FSI_outer_sweeps>=1)&&
      (FSI_outer_sweeps<num_FSI_outer_sweeps)) {

   for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

    MultiFab& Umac_new=get_new_data(Umac_Type+dir,slab_step+1);

    if (thread_class::nthreads<1)
     amrex::Error("thread_class::nthreads invalid");
    thread_class::init_d_numPts(S_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
    for (MFIter mfi(S_new,use_tiling);mfi.isValid(); ++mfi) {
     BL_ASSERT(grids[mfi.index()] == mfi.validbox());
     int gridno=mfi.index();
     const Box& tilegrid = mfi.tilebox();
     const Box& fabgrid = grids[gridno];
     const int* tilelo=tilegrid.loVect();
     const int* tilehi=tilegrid.hiVect();
     const int* fablo=fabgrid.loVect();
     const int* fabhi=fabgrid.hiVect();
     int bfact=parent->Space_blockingFactor(level);
     const Real* xlo = grid_loc[gridno].lo();
    
     FArrayBox& velMAC=Umac_new[mfi];
     FArrayBox& velCELL=S_new[mfi];

     FArrayBox& lsfab=(*localMF[LEVELPC_MF])[mfi];

     FArrayBox& FSIvelMAC=(*localMF[FSI_MAC_VELOCITY_MF+dir])[mfi];
     if (FSIvelMAC.nComp()!=1)
      amrex::Error("FSIvelMAC.nComp() invalid");
     FArrayBox& FSIvelCELL=(*localMF[FSI_CELL_VELOCITY_MF])[mfi];
     if (FSIvelCELL.nComp()!=AMREX_SPACEDIM)
      amrex::Error("FSIvelCELL.nComp() invalid");

     // mask=tag if not covered by level+1 or outside the domain.
     FArrayBox& maskcoeffab=(*localMF[MASKCOEF_MF])[mfi];

     Vector<int> velbc=getBCArray(State_Type,gridno,
       STATECOMP_VEL,STATE_NCOMP_VEL);

     int extend_solid_velocity=0;

     int tid_current=ns_thread();
     if ((tid_current<0)||(tid_current>=thread_class::nthreads))
      amrex::Error("tid_current invalid");
     thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

      // fort_manage_elastic_velocity is declared in: LEVELSET_3D.F90
     fort_manage_elastic_velocity(
      &extend_solid_velocity,
      im_elastic_map.dataPtr(),
      &num_FSI_outer_sweeps,
      &FSI_outer_sweeps,
      &dir, //dir=0,1,2
      velbc.dataPtr(),  
      &slab_step,
      &cur_time_slab, 
      xlo,dx,
      // mask=tag if not covered by level+1 or outside the domain.
      maskcoeffab.dataPtr(),
      ARLIM(maskcoeffab.loVect()),ARLIM(maskcoeffab.hiVect()),
      lsfab.dataPtr(),
      ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
      velMAC.dataPtr(),
      ARLIM(velMAC.loVect()),ARLIM(velMAC.hiVect()), 
      velCELL.dataPtr(STATECOMP_VEL+dir),
      ARLIM(velCELL.loVect()),ARLIM(velCELL.hiVect()),
      FSIvelMAC.dataPtr(),
      ARLIM(FSIvelMAC.loVect()),ARLIM(FSIvelMAC.hiVect()), 
      FSIvelCELL.dataPtr(dir),
      ARLIM(FSIvelCELL.loVect()),ARLIM(FSIvelCELL.hiVect()),
      tilelo,tilehi,
      fablo,fabhi,
      &bfact,
      &level,&finest_level,
      &NS_geometry_coord,
      domlo,domhi);
    } // mfi
} // omp
    ns_reconcile_d_num(LOOP_MANAGE_VEL,"fort_manage_elastic_velcity");

   } // dir=0..sdim-1

   delete_localMF(FSI_MAC_VELOCITY_MF,AMREX_SPACEDIM);
   delete_localMF(FSI_CELL_VELOCITY_MF,1);
  } else if (FSI_outer_sweeps==0) {
   //do nothing
  } else
   amrex::Error("FSI_outer_sweeps invalid");

  if ((FSI_outer_sweeps>=0)&&
      (FSI_outer_sweeps<num_FSI_outer_sweeps-1)) {

   new_localMF(FSI_CELL_VELOCITY_MF,AMREX_SPACEDIM,1,-1);
   for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
    //ngrow=2
    //Umac_Type+dir
    //getStateMAC_localMF is declared in NavierStokes2.cpp
    //localMF[FSI_MAC_VELOCITY_MF+dir allocated in getStateMAC_localMF.
    getStateMAC_localMF(FSI_MAC_VELOCITY_MF+dir,2,dir,cur_time_slab);
   } // dir=0 ... sdim-1

   // advect_register has 1 ghost initialized.
   push_back_state_register(FSI_CELL_VELOCITY_MF,cur_time_slab);

   if (localMF[FSI_CELL_VELOCITY_MF]->nGrow()>=1) {
    //do nothing
   } else
    amrex::Error("localMF[FSI_CELL_VELOCITY_MF]->nGrow()<1");

   if (localMF[FSI_CELL_VELOCITY_MF]->nComp()==AMREX_SPACEDIM) {
    //do nothing
   } else
    amrex::Error("localMF[FSI_CELL_VELOCITY_MF]->nComp()!=sdim");

   for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

    MultiFab& Umac_new=get_new_data(Umac_Type+dir,slab_step+1);

    if (thread_class::nthreads<1)
     amrex::Error("thread_class::nthreads invalid");
    thread_class::init_d_numPts(S_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
    for (MFIter mfi(S_new,use_tiling);mfi.isValid(); ++mfi) {
     BL_ASSERT(grids[mfi.index()] == mfi.validbox());
     int gridno=mfi.index();
     const Box& tilegrid = mfi.tilebox();
     const Box& fabgrid = grids[gridno];
     const int* tilelo=tilegrid.loVect();
     const int* tilehi=tilegrid.hiVect();
     const int* fablo=fabgrid.loVect();
     const int* fabhi=fabgrid.hiVect();
     int bfact=parent->Space_blockingFactor(level);
     const Real* xlo = grid_loc[gridno].lo();
    
     FArrayBox& velMAC=Umac_new[mfi];
     FArrayBox& velCELL=S_new[mfi];

     FArrayBox& lsfab=(*localMF[LEVELPC_MF])[mfi];

     FArrayBox& FSIvelMAC=(*localMF[FSI_MAC_VELOCITY_MF+dir])[mfi];
     if (FSIvelMAC.nComp()!=1)
      amrex::Error("FSIvelMAC.nComp() invalid");
     FArrayBox& FSIvelCELL=(*localMF[FSI_CELL_VELOCITY_MF])[mfi];
     if (FSIvelCELL.nComp()!=AMREX_SPACEDIM)
      amrex::Error("FSIvelCELL.nComp() invalid");

     // mask=tag if not covered by level+1 or outside the domain.
     FArrayBox& maskcoeffab=(*localMF[MASKCOEF_MF])[mfi];

     Vector<int> velbc=getBCArray(State_Type,gridno,
       STATECOMP_VEL,STATE_NCOMP_VEL);

     int extend_solid_velocity=1;

     int tid_current=ns_thread();
     if ((tid_current<0)||(tid_current>=thread_class::nthreads))
      amrex::Error("tid_current invalid");
     thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

      // fort_manage_elastic_velocity is declared in: LEVELSET_3D.F90
     fort_manage_elastic_velocity(
      &extend_solid_velocity,
      im_elastic_map.dataPtr(),
      &num_FSI_outer_sweeps,
      &FSI_outer_sweeps,
      &dir, //dir=0,1,2
      velbc.dataPtr(),  
      &slab_step,
      &cur_time_slab, 
      xlo,dx,
      // mask=tag if not covered by level+1 or outside the domain.
      maskcoeffab.dataPtr(),
      ARLIM(maskcoeffab.loVect()),ARLIM(maskcoeffab.hiVect()),
      lsfab.dataPtr(),
      ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
      velMAC.dataPtr(),
      ARLIM(velMAC.loVect()),ARLIM(velMAC.hiVect()), 
      velCELL.dataPtr(STATECOMP_VEL+dir),
      ARLIM(velCELL.loVect()),ARLIM(velCELL.hiVect()),
      FSIvelMAC.dataPtr(),
      ARLIM(FSIvelMAC.loVect()),ARLIM(FSIvelMAC.hiVect()), 
      FSIvelCELL.dataPtr(dir),
      ARLIM(FSIvelCELL.loVect()),ARLIM(FSIvelCELL.hiVect()),
      tilelo,tilehi,
      fablo,fabhi,
      &bfact,
      &level,&finest_level,
      &NS_geometry_coord,
      domlo,domhi);
    } // mfi
} // omp
    ns_reconcile_d_num(LOOP_MANAGE_VEL,"fort_manage_elastic_velcity");

   } // dir=0..sdim-1

  } else if (FSI_outer_sweeps==num_FSI_outer_sweeps-1) {
   //do nothing
  } else
   amrex::Error("FSI_outer_sweeps invalid");

 } else
  amrex::Error("expecting num_FSI_outer_sweeps>=2 and <= num_materials");

} // end subroutine manage_FSI_data()


void NavierStokes::prepare_advect_vars(Real time) {

 if (time>=0.0) {
  //do nothing
 } else
  amrex::Error("time invalid");

 new_localMF(ADVECT_REGISTER_MF,AMREX_SPACEDIM,1,-1);
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  //ngrow=1
  //Umac_Type+dir
  //getStateMAC_localMF is declared in NavierStokes2.cpp
  //localMF[ADVECT_REGISTER_FACE_MF+dir allocated in getStateMAC_localMF.
  getStateMAC_localMF(ADVECT_REGISTER_FACE_MF+dir,1,dir,time);
 } // dir

 // advect_register has 1 ghost initialized.
 push_back_state_register(ADVECT_REGISTER_MF,time);

 if (localMF[ADVECT_REGISTER_MF]->nGrow()>=1) {
  //do nothing
 } else
  amrex::Error("localMF[ADVECT_REGISTER_MF]->nGrow()<1");

 if (localMF[ADVECT_REGISTER_MF]->nComp()==AMREX_SPACEDIM) {
  //do nothing
 } else
  amrex::Error("localMF[ADVECT_REGISTER_MF]->nComp()!=sdim");

} // end subroutine prepare_advect_vars(Real time)


void NavierStokes::prepare_transport_vars(Real time) {

 if (time>=0.0) {
  //do nothing
 } else
  amrex::Error("time invalid");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  //ngrow=1
  //Umac_Type+dir
  //getStateMAC_localMF is declared in NavierStokes2.cpp
  //localMF[TRANSPORT_REGISTER_FACE_MF+dir allocated in getStateMAC_localMF.
  getStateMAC_localMF(TRANSPORT_REGISTER_FACE_MF+dir,1,dir,time);
 } // dir

} // end subroutine prepare_transport_vars(Real time)



void NavierStokes::alloc_DTDtALL(int alloc_flag) {

 if (level==0) {
  // do nothing
 } else
  amrex::Error("level invalid");

 int finest_level=parent->finestLevel();
 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.alloc_DTDt(alloc_flag);
 }

} // end subroutine alloc_DTDtALL

void NavierStokes::alloc_DTDt(int alloc_flag) {

 Vector<int> scomp;  
 Vector<int> ncomp;  
 int state_index;
 int ncomp_check;
 int ngrow=0;

  // DTdt_MF=T_new - T_advect_MF

 int num_materials_combine=num_materials;
 get_mm_scomp_solver(
   num_materials_combine,
   SOLVETYPE_HEAT,
   state_index,
   scomp,
   ncomp,
   ncomp_check);

  // called after phase change, and before thermal diffusion.
 if (alloc_flag==1) {

  getState_localMF_list(
   T_advect_MF,ngrow,
   state_index,
   scomp,
   ncomp);

  // delete DTDt_MF and T_advect_MF after mdot appended with
  // V_T rho DTDt - average(V_T rho DTDt) 
 } else if (alloc_flag==0) {

  delete_localMF(T_advect_MF,1);
  delete_localMF(DTDt_MF,1);

 } else if (alloc_flag==2) {  // DTDt=T_new - T_advect

  getState_localMF_list(
   DTDt_MF,ngrow,
   state_index,
   scomp,
   ncomp);

  MultiFab::Subtract(*localMF[DTDt_MF],*localMF[T_advect_MF],0,0,
		  num_materials,ngrow);
 } else
  amrex::Error("alloc_flag invalid in alloc_DTDt");

} // end subroutine alloc_DTDt


void NavierStokes::prepare_viscous_solver() {

 std::string local_caller_string="prepare_viscous_solver";

 int nsolve=AMREX_SPACEDIM;

 new_localMF(REGISTER_MARK_MF,nsolve,1,-1);
 new_localMF(HOOP_FORCE_MARK_MF,nsolve,1,-1);
 new_localMF(VISCHEAT_SOURCE_MF,nsolve,1,-1);

 new_localMF(VISCHEAT_MF,1,0,-1);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   // store stress fluxes (-pI+2 mu D) 
  new_localMF(CONSERVE_FLUXES_MF+dir,nsolve,0,dir);
 } // dir

 setVal_localMF(REGISTER_MARK_MF,0.0,0,nsolve,1);
 setVal_localMF(VISCHEAT_SOURCE_MF,0.0,0,nsolve,1);

 localMF[VISCHEAT_MF]->setVal(0.0,0,1,0);

 resize_metrics(1);
 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(FACE_VAR_MF,0,local_caller_string);
 debug_ngrow(VOLUME_MF,1,local_caller_string);
 debug_ngrow(MASKCOEF_MF,1,local_caller_string);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  if (localMF[AREA_MF+dir]->boxArray()!=
      localMF[FACE_VAR_MF+dir]->boxArray())
   amrex::Error("boxarrays do not match");
 }

}  // prepare_viscous_solver

// probtype=31 translating circle or sphere
int NavierStokes::is_zalesak() {

 if (fort_is_passive_advect_test()==1) {
  return 1;
 } else if (fort_is_passive_advect_test()==0) {
  return 0;
 } else {
  amrex::Error("fort_is_passive_advect_test() invalid");
  return 0;
 }

}

// velocity scaled by global_velocity_scale
void NavierStokes::zalesakVEL() {
 
 bool use_tiling=ns_tiling;

 const Real* dx = geom.CellSize();

 const Box& domain = geom.Domain();
 const int* domlo = domain.loVect();
 const int* domhi = domain.hiVect();

 MultiFab& U_new=get_new_data(State_Type,slab_step+1);
 U_new.setVal(0.0,STATECOMP_PRES,1,1);

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(U_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(U_new,use_tiling); mfi.isValid(); ++mfi) {
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
   FArrayBox& velfab = U_new[mfi];

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // takes into consideration global_velocity_scale
   fort_zalesak_cell(
    xlo,dx,
    velfab.dataPtr(),
    domlo,domhi,
    tilelo,tilehi,
    fablo,fabhi,
    &bfact,
    &level,
    ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
    &cur_time_slab);

 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_ZALESAK,"zalesakVEL");
}

#undef profile_solver

}/* namespace amrex */

