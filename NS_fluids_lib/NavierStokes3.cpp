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
#include <GODUNOV_F.H>
#include <NAVIERSTOKES_F.H>
#include <MACOPERATOR_F.H>
#include <PROB_F.H>
#include <LEVEL_F.H>
#include <SOLIDFLUID_F.H>
#include <DERIVE_F.H>
#include <DIFFUSION_F.H>
#include <Zeyu_Matrix_Functions.H>

namespace amrex{

#define profile_solver 0

// status==1 success
// status==0 failure
extern void matrix_solveCPP(Real** AA,Real* xx,Real* bb,
 int& status,int numelem);
extern void GMRES_MIN_CPP(Real** HH,Real beta, Real* yy,
 int m,int m_small,
 int caller_id,int project_option,
 int mg_level,int& status);

// if ncomp_input==-1, then ncomp=S_crse.ncomp()
// spectral_override==0 => always do low order average down.
void
NavierStokes::avgDownEdge(int dir,MultiFab& S_crse,MultiFab& S_fine,
 int scomp,int ncomp_input,int spectral_override,int caller_id) {

 if (1==0) {
  std::cout << "avgDownEdge caller_id= " << caller_id << '\n';
 }

 if ((dir<0)||(dir>=AMREX_SPACEDIM))
  amrex::Error("dir invalid avgdown edge");

 int finest_level=parent->finestLevel();
 if (level>=finest_level) 
  amrex::Error("level invalid in avgDownEdge");

 int f_level=level+1;
 NavierStokes& fine_lev = getLevel(f_level);
 resize_metrics(1);
 debug_ngrow(VOLUME_MF,0,700);
 fine_lev.resize_metrics(1);
 fine_lev.debug_ngrow(VOLUME_MF,0,700);

 if (fine_lev.localMF[AREA_MF+dir]->boxArray()!=S_fine.boxArray()) {
  const BoxArray& farray=(*fine_lev.localMF[AREA_MF+dir]).boxArray();
  const BoxArray& sarray=S_fine.boxArray();
  std::cout << "farray " << farray << '\n';
  std::cout << "sarray " << sarray << '\n';
  std::cout << "dir " << dir << '\n';
  std::cout << "level,finest_level " << level << ' ' << finest_level <<'\n';

  amrex::Error("invalid boxes in avgDownEdge fine level");
 }
 if ((*localMF[AREA_MF+dir]).boxArray()!=S_crse.boxArray())
  amrex::Error("invalid boxes in avgDownEdge crse level");

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

  FORT_EDGEAVGDOWN(
   &enable_spectral,
   &finest_level,
   &spectral_override,
   prob_lo,
   dxf,
   &level,&f_level,
   &bfact,&bfact_f,
   xlo_fine,dx,
   &dir,
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
 ns_reconcile_d_num(169);

 S_crse.copy(crse_S_fine_MAC,0,scomp,ncomp);
 ParallelDescriptor::Barrier();

 const Box& domain = geom.Domain();
 if (geom.isPeriodic(dir)) {
  IntVect pshift=IntVect::TheZeroVector();
  pshift[dir]=domain.length(dir);
  crse_S_fine_MAC.shift(pshift);

  ParallelDescriptor::Barrier();
  S_crse.copy(crse_S_fine_MAC,0,scomp,ncomp);
  ParallelDescriptor::Barrier();

  pshift[dir]=-2*domain.length(dir);
  crse_S_fine_MAC.shift(pshift);

  S_crse.copy(crse_S_fine_MAC,0,scomp,ncomp);
  ParallelDescriptor::Barrier();
 }  // isPeriodic(dir)

}  // avgDownEdge


// called from updatevelALL,multiphase_project
// interpolate from level+1 to level.
void
NavierStokes::avgDownMac() {

 int ncomp_edge=-1;
 int scomp=0;
  // spectral_override==0 => always do low order average down.
 int spectral_override=1;
 avgDownEdge_localMF(UMAC_MF,scomp,ncomp_edge,0,AMREX_SPACEDIM, 
   spectral_override,14);

}  // avgDownMac

// spectral_override==0 => always do low order average down.
void NavierStokes::avgDownMacState(int MAC_state_idx,int spectral_override) {

 int finest_level = parent->finestLevel();

 if (level>=finest_level) 
  amrex::Error("level invalid avgDownMacState");

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 int nsolveMM_FACE=num_materials_vel;

 NavierStokes& fine_lev = getLevel(level+1);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  MultiFab& S_crse = get_new_data(MAC_state_idx+dir,slab_step+1);
  MultiFab& S_fine = fine_lev.get_new_data(MAC_state_idx+dir,slab_step+1);
  int scomp=0;
  int ncomp_edge=S_crse.nComp(); 
 
  if (MAC_state_idx==Umac_Type) {
   if (ncomp_edge==nsolveMM_FACE) {
    // do nothing
   } else
    amrex::Error("ncomp_edge invalid in avgDownMacState");
  } else if (MAC_state_idx==XDmac_Type) {
   if (ncomp_edge==num_materials_viscoelastic) {
    // do nothing
   } else
    amrex::Error("ncomp_edge invalid in avgDownMacState");
   if (spectral_override==0) {
    // do nothing
   } else
    amrex::Error("spectral_override invalid");
  } else
   amrex::Error("MAC_state_idx invalid");
   

  if ((S_crse.nComp()!=ncomp_edge)||
      (S_fine.nComp()!=ncomp_edge))
   amrex::Error("S_crse.nComp() or S_fine.nComp() invalid");

  int caller_id=1;
  avgDownEdge(dir,S_crse,S_fine,scomp,ncomp_edge,spectral_override,caller_id);
 }  // dir=0..sdim-1

}  // subroutine avgDownMacState

void NavierStokes::nonlinear_advection() {

 int renormalize_only=1;

 if (level!=0)
  amrex::Error("level invalid nonlinear_advection");

 int nmat=num_materials;

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

 if (std::abs(cur_time_slab-prev_time_slab-dt_slab)>1.0E-5) {
  std::cout << "cur_time_slab " << cur_time_slab << '\n';
  std::cout << "prev_time_slab " << prev_time_slab << '\n';
  std::cout << "dt_slab " << dt_slab << '\n';
  amrex::Error("slab time bust1");
 }

 if (verbose>0)
  if (ParallelDescriptor::IOProcessor())
   std::cout << "nonlinear advect \n";

 NavierStokes& ns_fine=getLevel(finest_level);
 int basestep=ns_fine.nStep();

 order_direct_split=basestep-2*(basestep/2);

 if ((order_direct_split!=0)&&
     (order_direct_split!=1))
  amrex::Error("order_direct_split invalid");

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.new_localMF(VOF_LS_PREV_TIME_MF,2*nmat,1,-1);
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


 if (verbose>0)
  if (ParallelDescriptor::IOProcessor())
   std::cout << "order_direct_split " << order_direct_split << '\n';

 dir_absolute_direct_split=0;
 int init_vof_ls_prev_time=1;
 int normdir_here=normdir_direct_split[dir_absolute_direct_split];
 if ((normdir_here<0)||(normdir_here>=AMREX_SPACEDIM))
  amrex::Error("normdir_here invalid (prior to loop)");
 advect_time_slab=prev_time_slab;
 int update_flag=0;  // do not update the error. 
 VOF_Recon_ALL(1,advect_time_slab,update_flag,
  init_vof_ls_prev_time,SLOPE_RECON_MF);

 if (face_flag==1) {
   // delete_advect_vars() called in NavierStokes::do_the_advance
   // right after increment_face_velocityALL. 
  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
     // initialize ADVECT_REGISTER_FACE_MF and
     //            ADVECT_REGISTER_MF
   ns_level.prepare_advect_vars(prev_time_slab);
  }

 } else if (face_flag==0) {
  // do nothing
 } else
  amrex::Error("face_flag invalid");

 NavierStokes& ns_finest=getLevel(finest_level);
 int ipart_id=0;
 for (int im=0;im<nmat;im++) {
  if (particleLS_flag[im]==1) {

   // 1. void addParticles (const ParticleContainerType& other, 
   //     bool local=false);  (local==true => do not redistribute at end?)
   // 2. Copy Eulerian level set to Lagrangian particles.
   // 3. advect the particles using RK
   // (note: interface particles distinguished by levelset==0.0)
   // in: NavierStokes2.cpp
   ns_finest.move_particles(im,ipart_id);
   ipart_id++;
  } else if (particleLS_flag[im]==0) {
   // do nothing
  } else
   amrex::Error("particleLS_flag[im] invalid");
 } //im=0..nmat-1

 if (ipart_id==NS_ncomp_particles) {
  // do nothing
 } else
  amrex::Error("ipart_id invalid");

 for (dir_absolute_direct_split=0;
      dir_absolute_direct_split<AMREX_SPACEDIM;
      dir_absolute_direct_split++) {

   normdir_here=normdir_direct_split[dir_absolute_direct_split];
   if ((normdir_here<0)||(normdir_here>=AMREX_SPACEDIM))
    amrex::Error("normdir_here invalid (in loop)");

   init_vof_ls_prev_time=0;

   if (dir_absolute_direct_split==0) {
    advect_time_slab=prev_time_slab;
   } else if ((dir_absolute_direct_split>0)&&
              (dir_absolute_direct_split<AMREX_SPACEDIM)) {
    advect_time_slab=cur_time_slab;
    update_flag=0;  // do not update the error. 
    VOF_Recon_ALL(1,advect_time_slab,update_flag,
     init_vof_ls_prev_time,SLOPE_RECON_MF);

     // uses SLOPE_RECON_MF
     // prior to split_scalar_advection, in 
     //  NavierStokes::nonlinear_advection()
    if ((num_materials_viscoelastic>=1)&&
	(num_materials_viscoelastic<=nmat)) {
     for (int ilev=finest_level;ilev>=level;ilev--) {
      NavierStokes& ns_level=getLevel(ilev);
      ns_level.tensor_extrapolate(); //declared in: NavierStokes.cpp
     }
    } else
     amrex::Error("num_materials_viscoelastic invalid");

    avgDownALL_TENSOR();

   } else
    amrex::Error("dir_absolute_direct_split invalid");

    // order_direct_split=base_step mod 2=0 or 1
    // must go from finest level to coarsest.
   for (int ilev=finest_level;ilev>=level;ilev--) {
    NavierStokes& ns_level=getLevel(ilev);
    ns_level.split_scalar_advection();
   } // ilev

   if ((dir_absolute_direct_split>=0)&&
       (dir_absolute_direct_split<AMREX_SPACEDIM-1)) {
     // in: nonlinear_advection
     // calls MOFavgDown, LS_Type avgDown
     // projects volume fractions so that sum F_m_fluid=1.
    renormalize_only=1;
    int local_truncate=0;
    prescribe_solid_geometryALL(prev_time_slab,renormalize_only,
      local_truncate);

     // velocity and pressure
     // spectral_override==1 => order derived from "enable_spectral"
    avgDownALL(State_Type,0,num_materials_vel*(AMREX_SPACEDIM+1),1);
     // "state" (all materials)
    int scomp_den=num_materials_vel*(AMREX_SPACEDIM+1);
    avgDownALL(State_Type,scomp_den,num_state_material*nmat,1);
   } else if (dir_absolute_direct_split==AMREX_SPACEDIM-1) {
     // do nothing
   } else
    amrex::Error("parameter bust");

 }  // dir_absolute_direct_split=0..sdim-1

 for (int ilev=level;ilev<=finest_level;ilev++) {
   NavierStokes& ns_level=getLevel(ilev);
   ns_level.delete_localMF(VOF_LS_PREV_TIME_MF,1);

   ns_level.delete_localMF(MAC_VELOCITY_MF,AMREX_SPACEDIM);
   ns_level.delete_localMF(CELL_VELOCITY_MF,1);

    // sanity check
   ns_level.debug_ngrow(MASKCOEF_MF,1,6003);
 } // ilev=level..finest_level

// 0. convert Lagrangian description to Eulerian if read_from_CAD==1
// 1. renormalize variables
// 2. extend from F>0 fluid regions into F=0 regions
// 3. if renormalize_only==0, 
//    a. init F,X,LS for the solid materials.
//    b. init U,T in the solid regions.
//    c. extrapolate F,X,LS from fluid regions into solid regions.

 if (read_from_CAD()==1) {
  int iter=0;
  int FSI_operation=4; // eul vel t=cur_time_slab -> structure vel
  int FSI_sub_operation=0;
  for (FSI_sub_operation=0;FSI_sub_operation<3;FSI_sub_operation++) {
   for (int ilev=level;ilev<=finest_level;ilev++) {
    NavierStokes& ns_level=getLevel(ilev);
    ns_level.resize_mask_nbr(ngrowFSI);
    ns_level.ns_header_msg_level(
      FSI_operation,FSI_sub_operation,
      cur_time_slab,dt_slab,iter);
   } // ilev=level..finest_level
  } // FSI_sub_operation=0,1,2

   // FORT_HEADERMSG (SOLIDFLUID.F90)
   // CLSVOF_ReadNodes (sci_clsvof.F90)
   // if FSI_flag==4, then
   //  a) CTML_SOLVE_SOLID is called (in CTMLFSI.F90)
   //  b) tick is called (in ../Vicar3D/distFSI/tick.F)
   //
  FSI_operation=1; // update node locations
  FSI_sub_operation=0;
  ns_header_msg_level(FSI_operation,FSI_sub_operation,
   cur_time_slab,dt_slab,iter);
 } else if (read_from_CAD()==0) {
  // do nothing
 } else
  amrex::Error("read_from_CAD() invalid");

 // convert Lagrangian position, velocity, temperature, and force to
 // Eulerian.
 // go from coarsest to finest.
 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.FSI_make_distance(cur_time_slab,dt_slab);
 } // ilev

 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.resize_FSI_MF();
 }



 if (1==0) {
    // S_new is level 0 data
  MultiFab& S_new=get_new_data(State_Type,slab_step+1);
  int caller_id=1;
   // data file name "BEFOREPRESCRIBE<stuff>.plt"
   // xvel,yvel,zvel,pressure,(density, temperature) x nmat,
   // (VFRAC,centroid) x nmat, error indicator
  writeSanityCheckData(
   "BEFOREPRESCRIBE",
   "in: NavierStokes::nonlinear_advection, State_Type ", 
   caller_id,
   S_new.nComp(),
   -1, // data_mf==-1
   State_Type,
   -1); // data_dir==-1
 }


  // in: nonlinear_advection
  // level set function, volume fractions, and centroids are
  // made "consistent" amongst the levels.
  // in: NavierStokes2.cpp
 renormalize_only=0;
 int local_truncate=1;
 prescribe_solid_geometryALL(cur_time_slab,renormalize_only,
   local_truncate);

 avgDownALL(State_Type,0,
  num_materials_vel*(AMREX_SPACEDIM+1),1);


 if (1==0) {
    // S_new is level 0 data
  MultiFab& S_new=get_new_data(State_Type,slab_step+1);
  int caller_id=1;
   // data file name "AFTERPRESCRIBE<stuff>.plt"
   // xvel,yvel,zvel,pressure,(density, temperature) x nmat,
   // (VFRAC,centroid) x nmat, error indicator
  writeSanityCheckData(
   "AFTERPRESCRIBE",
   "in: NavierStokes::nonlinear_advection, State_Type ", 
   caller_id,
   S_new.nComp(),
   -1, // data_mf==-1
   State_Type,
   -1); // data_dir==-1
 }

}  // subroutine nonlinear_advection


void NavierStokes::allocate_SDC() {

 if ((ns_time_order==1)&&
     (enable_spectral==3))
  amrex::Error("(ns_time_order==1)&&(enable_spectral==3)");
 if ((ns_time_order>=2)&&
     ((enable_spectral==0)||
      (enable_spectral==2)))
  amrex::Error("(ns_time_order>=2)&&(enable_spectral==0 or 2)");

 if ((ns_time_order>=2)||
     (enable_spectral==1)||   // spectral in space and time
     (enable_spectral==2)) {  // spectral in space only

   // I-scheme,thermal conduction, viscosity, div(up),gp, -force
  if (localMF_grow[stableF_MF]==-1) {
   new_localMF(stableF_MF,nstate_SDC*ns_time_order,0,-1);
  } else
   amrex::Error("localMF_grow[stableF_MF] invalid");

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
   new_localMF(stableF_GP_MF+dir,ns_time_order,0,dir);

  if (localMF_grow[spectralF_MF]==-1) {
   new_localMF(spectralF_MF,nstate_SDC*(ns_time_order+1),0,-1);
  } else
   amrex::Error("localMF_grow[spectralF_MF] invalid");

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
   new_localMF(spectralF_GP_MF+dir,(ns_time_order+1),0,dir);

  if (localMF_grow[delta_MF]==-1) {
   new_localMF(delta_MF,nstate_SDC*ns_time_order,0,-1);
  } else
   amrex::Error("localMF_grow[delta_MF] invalid");

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
   new_localMF(delta_GP_MF+dir,ns_time_order,0,dir);

 } else if ((ns_time_order==1)&&
	    ((enable_spectral==0)||
	     (enable_spectral==3))) {
  // do nothing
 } else
  amrex::Error("ns_time_order or enable_spectral invalid");

}  // subroutine allocate_SDC


void NavierStokes::deallocate_SDC() {

 if ((ns_time_order==1)&&
     (enable_spectral==3))
  amrex::Error("(ns_time_order==1)&&(enable_spectral==3)");
 if ((ns_time_order>=2)&&
     ((enable_spectral==0)||
      (enable_spectral==2)))
  amrex::Error("(ns_time_order>=2)&&(enable_spectral==0 or 2)");

 if ((ns_time_order>=2)||
     (enable_spectral==1)||   // spectral in space and time
     (enable_spectral==2)) {  // spectral in space only

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

 } else if ((ns_time_order==1)&&
	    ((enable_spectral==0)||
	     (enable_spectral==3))) {
  // do nothing
 } else
  amrex::Error("ns_time_order or enable_spectral invalid");

}  // subroutine deallocate_SDC

// called before veldiffuseALL() from NavierStokes::do_the_advance
// Second half of D^{upside down triangle}/Dt
void NavierStokes::tensor_advection_updateALL() {

 int finest_level=parent->finestLevel();

 int nmat=num_materials;

 if ((num_materials_viscoelastic>=1)&&(num_materials_viscoelastic<=nmat)) {

   // SLOPE_RECON_MF is a parameter but not used.
   // prior to tensor_advection_update(), in
   //  NavierStokes::tensor_advection_updateALL() 
  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   ns_level.tensor_extrapolate(); //declared in: NavierStokes.cpp
  }
  avgDownALL_TENSOR();

   // init_gradu_tensorALL fills CELLTENSOR_MF using these steps:
   // 1. find all velocity derivatives at faces.
   // 2. interpolate derivatives from faces to cells using 1-sided
   //    interpolation in the case that e.g. lsleft(im_primary)>=0
   //    but lsright(im_primary)<0.
   //    (im_primary is the main material in the cell, lspoint(im_primary)>=0)
  int do_alloc=1;
  int simple_AMR_BC_flag_viscosity=1;
  int im_tensor=-1;
  int idx_elastic_flux=-1;
  init_gradu_tensorALL(
    im_tensor,
    HOLD_VELOCITY_DATA_MF,
    do_alloc,
    CELLTENSOR_MF,
    FACETENSOR_MF,
    idx_elastic_flux,
    simple_AMR_BC_flag_viscosity);

   // in: NavierStokes.cpp
  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   ns_level.tensor_advection_update();
  }

  delete_array(CELLTENSOR_MF);
  delete_array(FACETENSOR_MF);

  for (int im=0;im<nmat;im++) {
   if ((particleLS_flag[im]==1)||
       (particleLS_flag[im]==0)) { 
    if (ns_is_rigid(im)==0) {
     if ((elastic_time[im]>=0.0)&&
         (elastic_viscosity[im]>=0.0)) {
      if (store_elastic_data[im]==1) {
       if (viscoelastic_model[im]==2) {
	  // particles only appear on the finest level.
          // The flexible substrate is wholly contained on
          // the finest level.
        NavierStokes& ns_finest=getLevel(finest_level);
        ns_finest.accumulate_PC_info(im);
       } else if ((viscoelastic_model[im]==1)||
  		  (viscoelastic_model[im]==0)||
		  (viscoelastic_model[im]==3)) {
        // do nothing
       } else
        amrex::Error("viscoelastic_model[im] invalid");
      } else if (store_elastic_data[im]==0) {
       // do nothing
      } else
       amrex::Error("store_elastic_data invalid");
     } else
      amrex::Error("elastic_time or elastic_viscosity invalid");
    } else if (ns_is_rigid(im)==1) {
     // do nothing
    } else
     amrex::Error("ns_is_rigid(im) invalid");
   } else
    amrex::Error("particleLS_flag[im] invalid");

  } // im=0..nmat-1

  avgDownALL_TENSOR();

   // uses SLOPE_RECON_MF
   // after tensor_advection_update(), in
   //  NavierStokes::tensor_advection_updateALL() 
  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   ns_level.tensor_extrapolate(); // in: NavierStokes.cpp
  }
  avgDownALL_TENSOR();

 } else
  amrex::Error("num_materials_viscoelastic invalid");

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

 if ((ns_time_order==1)&&
     (enable_spectral==3))
  amrex::Error("(ns_time_order==1)&&(enable_spectral==3)");
 if ((ns_time_order>=2)&&
     ((enable_spectral==0)||
      (enable_spectral==2)))
  amrex::Error("(ns_time_order>=2)&&(enable_spectral==0 or 2)");

 if ((ns_time_order>=2)||
     (enable_spectral==1)||   // spectral in space and time
     (enable_spectral==2)) {  // spectral in space only

  if (SDC_outer_sweeps==0) {
   setVal_localMF(delta_MF,0.0,0,nstate_SDC*ns_time_order,0);
   for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
    setVal_localMF(delta_GP_MF+dir,0.0,0,ns_time_order,0);
  } else if ((SDC_outer_sweeps>0)&&
             (SDC_outer_sweeps<ns_time_order)) {
   // do nothing
  } else
   amrex::Error("SDC_outer_sweeps invalid");

  setVal_localMF(stableF_MF,0.0,0,nstate_SDC*ns_time_order,0);
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
   setVal_localMF(stableF_GP_MF+dir,0.0,0,ns_time_order,0);

  int scomp=0;
  int ncomp=nstate_SDC*(ns_time_order+1);
  int scompGP=0;
  int ncompGP=ns_time_order+1;

  if ((SDC_outer_sweeps>0)&&
      (SDC_outer_sweeps<ns_time_order)) {
   // do not zap F(t^n) data.
   scomp=nstate_SDC;
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

 } else if ((ns_time_order==1)&&
	    ((enable_spectral==0)||
	     (enable_spectral==3))) {
  // do nothing
 } else
  amrex::Error("ns_time_order or enable_spectral invalid");

}  // subroutine init_delta_SDC

Real NavierStokes::advance(Real time,Real dt) {

 if (ParallelDescriptor::IOProcessor()) 
  std::cout << "advance time= " << time << " dt= " << dt << '\n';

 int finest_level = parent->finestLevel();
 const int max_level = parent->maxLevel();

 if (finest_level==max_level) {
  // do nothing
 } else
  amrex::Error("it is required that finest_level==max_level");

 int nmat=num_materials;
 
 Real dt_new=dt;
 int advance_status=1;

 if (level==0) {

  do {

   SDC_outer_sweeps=0;
   slab_step=ns_time_order-1;
   SDC_setup_step(); 

   if ((time>=0.0)&&(time<=1.0)) {
    if (std::abs(upper_slab_time-time)>1.0e-12)
     amrex::Error("upper_slab_time-time>tol (a)");
   } else if (time>1.0) {
    if (std::abs(upper_slab_time-time)>1.0e-12*time)
     amrex::Error("upper_slab_time-time>tol(time) (b)");
   } else
    amrex::Error("time invalid");

   if (verbose>0) {
    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "calling metrics \n";
    }
   }

   metrics_dataALL(1);  

   if (verbose>0) {
    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "initializing masks \n";
    }
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
    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "prescribe solid geometry (after regridding) \n";
    }
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
   avgDownALL(State_Type,0,
    num_materials_vel*(AMREX_SPACEDIM+1),1);
    // "state" (all materials)
   int scomp_den=num_materials_vel*(AMREX_SPACEDIM+1);
   avgDownALL(State_Type,scomp_den,num_state_material*nmat,1);
    // expected "DIV" 
   avgDownALL(DIV_Type,0,num_materials_vel,1);

    // in: advance
    // calls MOFavgDown, LS_Type avgDown
    // projects volume fractions so that sum F_m_fluid=1.
   int renormalize_only=0;
   int local_truncate=0;
   prescribe_solid_geometryALL(upper_slab_time,renormalize_only,
      local_truncate);

   if (verbose>0) {
    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "copy new to old ... \n";
    }
   }

   for (int ipart=0;ipart<NS_ncomp_particles;ipart++) {
    int lev_min=0;
    int lev_max=-1;
    int nGrow_Redistribute=0;
    int local_Redistribute=0;

    AmrParticleContainer<N_EXTRA_REAL,0,0,0>& old_PC=
       get_new_dataPC(State_Type,ns_time_order,ipart);
    old_PC.Redistribute(lev_min,lev_max,nGrow_Redistribute, 
      local_Redistribute);
   }

   CopyNewToOldALL();
   for (int ilev=level;ilev<=finest_level;ilev++) {
    NavierStokes& ns_level=getLevel(ilev);
    ns_level.setTimeLevel(time+dt_new,dt_new);
   }

   SDC_setup_step();

   Real time_scale=1.0;
   if (upper_slab_time>time_scale)
    time_scale=upper_slab_time;
   time_scale*=1.0E-10;

   if (std::abs(upper_slab_time-lower_slab_time-dt_new)>time_scale)
    amrex::Error("SDC_setup_step failed");
   if (std::abs(lower_slab_time-time)>time_scale)
    amrex::Error("lower_slab_time set improperly");

   do_the_advance(lower_slab_time,dt_new,advance_status);

   if (advance_status==1) {
    // do nothing (success)
   } else if (advance_status==0) { // failure
    dt_new=0.5*dt_new;
    CopyOldToNewALL();
    for (int ilev=level;ilev<=finest_level;ilev++) {
     NavierStokes& ns_level=getLevel(ilev);
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

} // subroutine advance



//   Stage| Description
//   -----+--------------------------------------
//    -1  +-initialization (0)
//     0  | Liquid cooling
//     1  +-Nucleation
//     2  | Recalescence in progress
//     3  +-Recalescence finished, frost starting (1)
//     4  | Frost  
//     5  +-Regular freezing starts (2) 
//     6  | Regular freezing in progress
//     6  | Solid cooling in progress
//

// (0): Change the ice surface tensions to values similar to
//      substarte. It's done in the calling function.
//      Noting to be done here.
// (1): Solid sruface tensions are set to be the same as the substrate
//      Surface tension to emulate ice growing out of the substarte,
//      T_w->T_f, and Z_w->Z_mush, Z={Cp, mu, k}
// (2): Ice surface tensions are set to be the original value

// Right now stage 1, 3, and 5 are internal only and cannot seen by
// the calling function.

#define INVALID_TIME -1000

void advance_recalesce(
	
 int nmat, // number of materials
 int nten, // number of surface tension coef
 Vector<int>& recal_material,  // 0=no 1=static 2=dynamic 
   // treatment, size=nmat
 Vector<Real>& freezing_state_old, 
  // nmat*(0-stage,1-touch_time,2-nucleation_time, 3-frost_time, 
  // 4-reg_freezing_time, 5-mushy_frac). size=6*nmat
 Vector<Real>& freezing_state_new, 
  // nmat*(0-stage,1-touch_time,2-nucleation_time, 3-frost_time,
  // 4-reg_freezing_time, 5-mushy_frac). size=6*nmat
 Vector<Real>& material_state, 
  // nmat*(0-T_average, 1-T_top, 2-dist_top, 3-dist_bottom, 
  // 4-touched, 5-centroid, 6-volume). size=nmat* (6+AMREX_SPACEDIM)
 Vector<Real>& exper_var, 
 // nmat*(0-trigger_temperature, 1-recal_speed, 2-frost_period). size=3*nmat
 Vector<Real>& latent_heat_source, // size=nmat
 Vector<Real>& heat_conduct, // size=nmat
 Vector<Real>& spec_heat, // size=nmat
 Real visc_factor,
 Vector<Real>& visc_coef, // size=nmat
 Vector<Real>& density, // size=nmat
 Vector<Real>& surf_ten, // size=nten
 Vector<Real>& TSAT, // size=2*nten
 Real time)
 {
  int nten_test=( (nmat-1)*(nmat-1)+nmat-1 )/2;
  if (nten_test!=nten)
   amrex::Error("nten invalid");
  int num_integrate=6+AMREX_SPACEDIM;
  if (material_state.size()!=nmat*num_integrate)
   amrex::Error("material_state invalid size");
  int recalesce_num_state=6;
  if (freezing_state_old.size()!=recalesce_num_state*nmat)
   amrex::Error("freezing_state_old invalid size");
  if (freezing_state_new.size()!=recalesce_num_state*nmat)
   amrex::Error("freezing_state_new invalid size");
  if (recal_material.size()!=nmat)
   amrex::Error("recal_material invalid size");
  if (exper_var.size()!=3*nmat)
   amrex::Error("exper_var invalid size");
  if (latent_heat_source.size()!=nmat)
   amrex::Error("latent_heat_source invalid size");
  if (heat_conduct.size()!=nmat)
   amrex::Error("heat_conduct invalid size");
  if (spec_heat.size()!=nmat)
   amrex::Error("spec_heat invalid size");
  if (visc_coef.size()!=nmat)
   amrex::Error("visc_coef invalid size");
  if (density.size()!=nmat)
   amrex::Error("density invalid size");
  if (surf_ten.size()!=nten)
   amrex::Error("surf_ten invalid size");
  if (TSAT.size()!=2*nten)
   amrex::Error("TSAT invalid size");

  Real drop_height, recal_duration;
   // Sanity checks
  for (int im=0; im<nmat; im++) {
   int ibase=im*num_integrate;
   // touched==1 -> touch_time <= time
   if((material_state[ibase+4]==1.0 &&
      (freezing_state_old[im*recalesce_num_state+1]>time))) {
    amrex::Error("multistage_freezing: invalid touch time!");
   }
  } // im

  int verbose=1;
  if (verbose>0){
   if (ParallelDescriptor::IOProcessor()) {

    for (int i=0; i<nmat; i++) {
     if(recal_material[i]!=0){
      std::cout << "TIME= " << time << " MAT= " << i << " stage     = ";
      std::cout << freezing_state_old[i*recalesce_num_state+0] << "\n";
      std::cout << "TIME= " << time << " MAT= " << i << " touch_time= ";
      std::cout << freezing_state_old[i*recalesce_num_state+1] << "\n";
      std::cout << "TIME= " << time << " MAT= " << i << " nucl_time = ";
      std::cout << freezing_state_old[i*recalesce_num_state+2] << "\n";
      std::cout << "TIME= " << time << " MAT= " << i << " fros_time = ";
      std::cout << freezing_state_old[i*recalesce_num_state+3] << "\n";
      std::cout << "TIME= " << time << " MAT= " << i << " regf_time = ";
      std::cout << freezing_state_old[i*recalesce_num_state+4] << "\n";
      std::cout << "TIME= " << time << " MAT= " << i << " mushy_frac= ";
      std::cout << freezing_state_old[i*recalesce_num_state+5] << "\n";
     }
    }

    for (int i=0; i<nmat; i++) {
     if(recal_material[i]!=0){
      std::cout << "TIME= " << time << " MAT= " << i << " T_average = ";
      std::cout << material_state[i*num_integrate+0] << "\n";
      std::cout << "TIME= " << time << " MAT= " << i << " T_top     = ";
      std::cout << material_state[i*num_integrate+1] << "\n";
      std::cout << "TIME= " << time << " MAT= " << i << " dist_top  = ";
      std::cout << material_state[i*num_integrate+2] << "\n";
      std::cout << "TIME= " << time << " MAT= " << i << " dist_botm = ";
      std::cout << material_state[i*num_integrate+3] << "\n";
      std::cout << "TIME= " << time << " MAT= " << i << " touched   = ";
      std::cout << material_state[i*num_integrate+4] << "\n";
      std::cout << "TIME= " << time << " MAT= " << i << " centroid  = ";
      D_TERM(std::cout <<        material_state[i*num_integrate+5];,
         std::cout << ", "<< material_state[i*num_integrate+6];,
      std::cout << ", "<< material_state[i*num_integrate+7];);
      std::cout << "\n";
      std::cout << "TIME= " << time << " MAT= " << i << " volume    = ";
      std::cout << material_state[i*num_integrate+AMREX_SPACEDIM+5] << "\n";
     }
    }
   } // IOProcessor
  } // verbose

  // Copy old stage data to new stage data...
  for (int i=0; i<freezing_state_old.size(); i++) {
   freezing_state_new[i]=freezing_state_old[i];
  }
  // and make nessecary changes for each material
  for (int im_source=0; im_source<nmat; im_source++) {

   if ((recal_material[im_source]==1) || 
       (recal_material[im_source]==2)) {
    Real LL=latent_heat_source[im_source];
    if (LL==0.0)
     amrex::Error("LL invalid");

    int ibase=im_source*recalesce_num_state; 
    int ibaseI=im_source*num_integrate; 
	////// stage -1: initialization
    if (freezing_state_new[ibase+0]==-1.0) {
	freezing_state_new[ibase+0]=0.0;
    } // end of stage -1

	////// stage 0: liquid cooling before impact
    if (freezing_state_new[ibase+0]==0.0) {
	// if no contact
     if (material_state[ibaseI+4]==0){
	// do nothing
      freezing_state_new[ibase+1]=INVALID_TIME;
     } else if (material_state[ibaseI+4]==1) { // if has contact
	// record the contact time if not recorded yet
      if (freezing_state_old[ibase+1]==INVALID_TIME){
	freezing_state_new[ibase+1]=time;
      }	
	// check for the nucleation critrria. Necleation starts at
	// specific top temperature that comes from experiment for
	// static droplet, or an deduced average droplet temperature for
	// the impact case.
      if (((recal_material[im_source]==1) && 
           (freezing_state_new[ibase+1]!=INVALID_TIME) &&
           (material_state[ibaseI+1]<=exper_var[im_source*3+0])) ||
          ((recal_material[im_source]==2) && 
           (freezing_state_new[ibase+1]!=INVALID_TIME) &&
           (material_state[ibaseI+0]<=exper_var[im_source*3+0]))) {
	// tansit to nucleation stage
       freezing_state_new[ibase+0]=1.0;
	// set nuclation_time
       freezing_state_new[ibase+2]=time;
      }
     } else {
      amrex::Error("invalid touch flag");
     }
    } // end of stage 0

	////// stage 1: nucleation
    if (freezing_state_new[ibase+0]==1.0) {
	// tansit to recalescence stage
     freezing_state_new[ibase+0]=2.0;
    } // end of stage 1

      ////// stage 2: recalescence
    if (freezing_state_new[ibase+0]==2.0) {
     drop_height = 
        material_state[ibaseI+2] -
	material_state[ibaseI+3];
     if (drop_height<0.0)
      amrex::Error("drop_height invalid");
	// recal_duration = drop_height / recal_speed
     recal_duration = drop_height / exper_var[im_source*3+1];
	// update end of recal stage
     freezing_state_new[ibase+3]=freezing_state_new[ibase+2] 
		+ recal_duration;
	// check if passed the end of recal stage
     if(time>=freezing_state_new[ibase+3]){
      freezing_state_new[ibase+0]=3.0;
     }
    } // end of stage 2
	
	////// stage 3: end of recal, initiating frost
    if (freezing_state_new[ibase+0]==3.0) {
	// tansit to frost stage
     freezing_state_new[ibase+0]=4.0;
     freezing_state_new[ibase+4]=freezing_state_new[ibase+3] 
	+ exper_var[im_source*3+2] ;
	// droplet temperature is changed to T_m (in the calling
	// function), and water properties are changed to the mushy
	// material properties (in the calling function)
    } // end of stage 3

	////// stage 4: frost
    if (freezing_state_new[ibase+0]==4.0) {
	// wait to pass the frost
     if (time>=freezing_state_new[ibase+4]){
      freezing_state_new[ibase+0]=5.0;
     }
    } // end of stage 4

	////// stage 5: initiation of regular freezing
    if (freezing_state_new[ibase+0]==5.0) {
	// The surface tensions are changed for solid material to
	// the physical values (in the calling function). 
     freezing_state_new[ibase+0]=6.0;
    } // end of stage 5

       ////// stage 6: regular freezung or solid cooling
    if (freezing_state_new[ibase+0]==6.0) {
     // do nothing
    } // end of stage 6

   } else if (recal_material[im_source]==0){
	// do nothing
   } else {
    amrex::Error("multistage_freezing: invalid freezing model!");
   } // recal_material
  } // im_source
 
} // advance_recalesce


void
NavierStokes::recalesce_temperature(int im_source) {

 
 bool use_tiling=ns_tiling;

 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;
 int finest_level=parent->finestLevel();

 if ((im_source<0)||(im_source>=nmat))
  amrex::Error("im_source invalid");

 Real TSAT=0.0;
 for (int im_opp=0;im_opp<nmat;im_opp++) {
  int iten;
  int ireverse;
  if (im_opp!=im_source) {
   get_iten_cpp(im_source+1,im_opp+1,iten,nmat);
   if ((iten<1)||(iten>nten))
    amrex::Error("iten invalid");
   if (im_source<im_opp)
    ireverse=0;
   else 
    ireverse=1;
   TSAT=saturation_temp[iten+ireverse*nten-1];
  }  // im_opp<>im_source
 } // im_opp
   
 if (TSAT<=0.0)
  amrex::Error("saturation temperature invalid");

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 int ncomp_state=S_new.nComp();
 if (ncomp_state!=num_materials_vel*(AMREX_SPACEDIM+1)+
     nmat*(num_state_material+ngeom_raw)+1)
  amrex::Error("ncomp_state invalid");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 const Real* dx = geom.CellSize();

 if (num_materials_vel!=1)
  amrex::Error("this code not ready yet for num_materials_vel!=1");

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
  Vector<int> velbc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);

  FArrayBox& snewfab=S_new[mfi];

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  FORT_RESET_TEMPERATURE(
   &im_source,
   &TSAT,
   snewfab.dataPtr(),
   ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()),
   xlo,dx,
   tilelo,tilehi,
   fablo,fabhi,&bfact,
   velbc.dataPtr(),
   &cur_time_slab,
   &nmat, 
   &ncomp_state,
   &level,&finest_level);

 } // mfi
} // omp
 ns_reconcile_d_num(171);

}  // recalesce_temperature



// sweep=0: integral rho T F, integral x rho F, integral rho F, 
//          dist_farthest, dist_closest, touch_flag 
// sweep=1: temperature_farthest
void
NavierStokes::process_recalesce_data(
 Vector<int> recalesce_material,
 Vector<Real>& recalesce_state_old,
 Vector<Real>& recalesce_state_new,
 Vector<Real>& integrated_quantities,
 int isweep) {
 

 int nmat=num_materials;

 int finest_level=parent->finestLevel();

 int recalesce_num_state=6;
 if (recalesce_material.size()!=nmat) 
  amrex::Error("invalid size for recalesce_material");
 if (recalesce_state_old.size()!=recalesce_num_state*nmat)
  amrex::Error("invalid size for recalesce_state_old");
 if (recalesce_state_new.size()!=recalesce_num_state*nmat)
  amrex::Error("invalid size for recalesce_state_new");

  // average temperature, temperature_farthest, dist_farthest,
  // dist_closest, touch_flag, xcentroid, mass
 int num_integrate=5+AMREX_SPACEDIM+1;
 if (integrated_quantities.size()!=nmat*num_integrate)
  amrex::Error("integrated_quantities invalid size");


 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 int ncomp_state=S_new.nComp();
 if (ncomp_state!=num_materials_vel*(AMREX_SPACEDIM+1)+
     nmat*(num_state_material+ngeom_raw)+1)
  amrex::Error("ncomp_state invalid");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 Vector<Real> local_integrated_quantities;
 local_integrated_quantities.resize(nmat*num_integrate);
 int iq,im;
 for (iq=0;iq<nmat*num_integrate;iq++)
  local_integrated_quantities[iq]=0.0;
 for (im=0;im<nmat;im++) {
  int ibase=im*num_integrate;
  iq=3; // dist_closest
  local_integrated_quantities[ibase+iq]=1.0e+10;
 }

 const Real* dx = geom.CellSize();

 if (num_materials_vel!=1)
  amrex::Error("this code not ready yet for num_materials_vel!=1");

 resize_metrics(1);

  // mask=tag if not covered by level+1 or outside the domain.
 int ngrowmask=0;
 int clear_phys_boundary=0;
 Real tag=1.0;
 MultiFab* mask=maskfiner(ngrowmask,tag,clear_phys_boundary);  

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(S_new.boxArray().d_numPts());

 for (MFIter mfi(S_new,false); mfi.isValid(); ++mfi) {
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
  Vector<int> velbc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);

  FArrayBox& maskfab=(*mask)[mfi];
  FArrayBox& volfab=(*localMF[VOLUME_MF])[mfi];

  FArrayBox& snewfab=S_new[mfi];

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  FORT_INTEGRATE_RECALESCE(
   &isweep,
   integrated_quantities.dataPtr(),
   local_integrated_quantities.dataPtr(),
   recalesce_material.dataPtr(),
   snewfab.dataPtr(),
   ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()),
   maskfab.dataPtr(),ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
   volfab.dataPtr(),ARLIM(volfab.loVect()),ARLIM(volfab.hiVect()),
   xlo,dx,
   tilelo,tilehi,
   fablo,fabhi,&bfact,
   velbc.dataPtr(),
   &cur_time_slab,
   &num_integrate,&nmat, 
   &ncomp_state,
   &level,&finest_level);

 } // mfi
 ns_reconcile_d_num(172);

 for (im=0;im<nmat;im++) {
  int ibase=im*num_integrate;
  Real A,B;

  // average temperature, temperature_farthest, dist_farthest,
  // dist_closest, touch_flag, xcentroid, mass
  if (isweep==0) {
   for (iq=5;iq<=5+AMREX_SPACEDIM;iq++) {
    ParallelDescriptor::ReduceRealSum(local_integrated_quantities[ibase+iq]);
    integrated_quantities[ibase+iq]+=local_integrated_quantities[ibase+iq];
   }
   iq=0; // avg temp
   ParallelDescriptor::ReduceRealSum(local_integrated_quantities[ibase+iq]);
   integrated_quantities[ibase+iq]+=local_integrated_quantities[ibase+iq];

   iq=2;  // dist_farthest
   ParallelDescriptor::ReduceRealMax(local_integrated_quantities[ibase+iq]);
   A=integrated_quantities[ibase+iq];
   B=local_integrated_quantities[ibase+iq];
   if (A<B) A=B;
   integrated_quantities[ibase+iq]=A;

   iq=3;  // dist_closest
   ParallelDescriptor::ReduceRealMin(local_integrated_quantities[ibase+iq]);
   A=integrated_quantities[ibase+iq];
   B=local_integrated_quantities[ibase+iq];
   if (A>B) A=B;
   integrated_quantities[ibase+iq]=A;
   iq=4;  // touch_flag
   ParallelDescriptor::ReduceRealMax(local_integrated_quantities[ibase+iq]);
   A=integrated_quantities[ibase+iq];
   B=local_integrated_quantities[ibase+iq];
   if (A<B) A=B;
   integrated_quantities[ibase+iq]=A;
  } else if (isweep==1) {

   iq=1;  // temp_farthest
   ParallelDescriptor::ReduceRealMax(local_integrated_quantities[ibase+iq]);
   A=integrated_quantities[ibase+iq];
   B=local_integrated_quantities[ibase+iq];
   if (A<B) A=B;
   integrated_quantities[ibase+iq]=A;

  } else
   amrex::Error("isweep invalid");
  
 } // im

 delete mask; 

}  // process_recalesce_data


void NavierStokes::process_recalesce_dataALL(
 Vector<int> recalesce_material,
 Vector<Real>& recalesce_state_old,
 Vector<Real>& recalesce_state_new) {

 int finest_level=parent->finestLevel();

 int recalesce_num_state=6;

 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

 if (recalesce_material.size()!=nmat)
  amrex::Error("invalid size for recalesce_material");
 if (recalesce_state_old.size()!=recalesce_num_state*nmat)
  amrex::Error("invalid size for recalesce_state_old");
 if (recalesce_state_new.size()!=recalesce_num_state*nmat)
  amrex::Error("invalid size for recalesce_state_new");

  // average temperature, temperature_farthest, dist_farthest,
  // dist_closest, touch_flag, xcentroid, mass
 int num_integrate=5+AMREX_SPACEDIM+1;
 Vector<Real> integrated_quantities;
 integrated_quantities.resize(nmat*num_integrate);
 int iq,im;
 for (iq=0;iq<integrated_quantities.size();iq++)
  integrated_quantities[iq]=0.0;
 for (im=0;im<nmat;im++) {
  int ibase=im*num_integrate;
  iq=3; // dist_closest
  integrated_quantities[ibase+iq]=1.0e+10;
 }

 for (int isweep=0;isweep<2;isweep++) {
  for (int ilev=level;ilev<=finest_level;ilev++) {
   NavierStokes& ns_level=getLevel(ilev);
   ns_level.process_recalesce_data(
     recalesce_material,
     recalesce_state_old,
     recalesce_state_new,
     integrated_quantities,isweep);
  } // ilev
 } // isweep
 for (im=0;im<nmat;im++) {
  int ibase=im*num_integrate;
  // average temperature, temperature_farthest, dist_farthest,
  // dist_closest, touch_flag, xcentroid, mass
  int iqvol=5+AMREX_SPACEDIM;
  Real massmat=integrated_quantities[ibase+iqvol];
  if (massmat>0.0) {
   iq=0;  // avg temp
   integrated_quantities[ibase+iq]=integrated_quantities[ibase+iq]/massmat;
   for (iq=5;iq<5+AMREX_SPACEDIM;iq++)
    integrated_quantities[ibase+iq]=integrated_quantities[ibase+iq]/massmat;
  } else if (massmat==0.0) {
   // do nothing
  } else
   amrex::Error("massmat invalid");
 } // im

 Vector<Real> latent_heat_source;
 latent_heat_source.resize(nmat);
 for (int im_source=0;im_source<nmat;im_source++) {
  Real LL=0.0;
  if ((recalesce_material[im_source]==1)||
      (recalesce_material[im_source]==2)) {

   for (int im_opp=0;im_opp<nmat;im_opp++) {
    int iten;
    int ireverse;
    if (im_opp!=im_source) {
     get_iten_cpp(im_source+1,im_opp+1,iten,nmat);
     if ((iten<1)||(iten>nten))
      amrex::Error("iten invalid");
     if (im_source<im_opp)
      ireverse=0;
     else 
      ireverse=1;
     LL=latent_heat[iten+ireverse*nten-1];
    }  // im_opp<>im_source
   } // im_opp
   
   if (LL==0.0)
    amrex::Error("latent_heat invalid");
  } else if (recalesce_material[im_source]==0) {
   // do nothing
  } else 
   amrex::Error("recalesce_material invalid");
  latent_heat_source[im_source]=LL;
 } // im_source

  // recalesce_material=0,1 (static),2 (dynamic)
 advance_recalesce(nmat,nten,
  recalesce_material,
  recalesce_state_old,
  recalesce_state_new,
  integrated_quantities,
  recalesce_model_parameters,
  latent_heat_source,
  heatviscconst,
  stiffCP,
  visc_coef,
  viscconst,
  denconst,
  tension,
  saturation_temp, 
  cur_time_slab); 

 for (im=0;im<nmat;im++) {
  if ((recalesce_material[im]==1)||
      (recalesce_material[im]==2)) {
   int ibase=im*recalesce_num_state;
   double stage=recalesce_state_old[ibase];
   double stage_new=recalesce_state_new[ibase];
   if ((stage<=2.5)&&(stage_new>=2.5)) {

    for (int ilev=level;ilev<=finest_level;ilev++) {
     NavierStokes& ns_level=getLevel(ilev);
     ns_level.recalesce_temperature(im);
    }

   }
  } else if (recalesce_material[im]==0) {
   // do nothing
  } else
   amrex::Error("recalesce material invalid");
 } // im=0..nmat-1

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   for (iq=0;iq<integrated_quantities.size();iq++) {
    std::cout << "recalesce  iq= " << iq << 
     " integrated_quantities= " <<
     integrated_quantities[iq] << '\n';
   }  // iq
  }  // io proc ?
 } // verbose>0?

} // process_recalesce_dataALL

// delta=integral_tn^tnp1  f^spectral dt - deltatn F^stable
void NavierStokes::init_splitting_force_SDC() {

 
 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();

 if (ns_time_order>=2) {

  if (slab_step!=ns_time_order)
   amrex::Error("slab_step invalid");

  if (num_state_base!=2)
   amrex::Error("num_state_base invalid");

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   debug_ngrow(spectralF_GP_MF+dir,0,3);
   debug_ngrow(stableF_GP_MF+dir,0,3);
   debug_ngrow(delta_GP_MF+dir,0,3);
  }
  debug_ngrow(spectralF_MF,0,3);
  debug_ngrow(stableF_MF,0,3);
  debug_ngrow(delta_MF,0,3);

  MultiFab& S_new=get_new_data(State_Type,slab_step+1);

  int nmat=num_materials;
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

   FArrayBox& HOfab=(*localMF[spectralF_MF])[mfi];
   FArrayBox& LOfab=(*localMF[stableF_MF])[mfi];
   FArrayBox& deltafab=(*localMF[delta_MF])[mfi];
   FArrayBox& masksem=(*localMF[MASKSEM_MF])[mfi];

   int HOncomp=HOfab.nComp();
   int LOncomp=LOfab.nComp();
   int delta_ncomp=deltafab.nComp();

   if ((HOncomp!=(ns_time_order+1)*nstate_SDC)||
       (LOncomp!=ns_time_order*nstate_SDC)||
       (delta_ncomp!=ns_time_order*nstate_SDC))
    amrex::Error("HOncomp, LOncomp, or delta_ncomp invalid");

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   FORT_SDC_TIME_QUAD(
    &HOncomp,
    &LOncomp,
    &delta_ncomp,
    &nstate,
    &nfluxSEM,
    &nstate_SDC,
    &nmat,
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
  ns_reconcile_d_num(173);

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

    FORT_SDC_TIME_QUAD_FACE(
     &dir,
     &HOncomp,
     &LOncomp,
     &delta_ncomp,
     &nstate,
     &nfluxSEM,
     &nstate_SDC,
     &nmat,
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
   ns_reconcile_d_num(174);

  } // dir=0..sdim-1

 } else {
  amrex::Error("ns_time_order invalid init_splitting_force_SDC");
 }

} // subroutine init_splitting_force_SDC

void NavierStokes::SEM_advectALL(int source_term) {

 if (stokes_flow==0) {

  int finest_level=parent->finestLevel();
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

   int operation_flag=7;

   int SEM_end_spectral_loop=2;
   if ((enable_spectral==1)||(enable_spectral==2)) {
    SEM_end_spectral_loop=2;
   } else if ((enable_spectral==0)||(enable_spectral==3)) {
    SEM_end_spectral_loop=1;
   } else {
    std::cout << "enable_spectral= " << enable_spectral << '\n';
    amrex::Error("enable_spectral invalid 1");
   }

   int nmat=num_materials;

   prescribed_vel_time_slab=prev_time_slab;
   vel_time_slab=prev_time_slab;

   if (source_term==1) {
    vel_time_slab=prev_time_slab;
   } else if (source_term==0) {
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
    for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
     ns_level.getStateMAC_localMF(
       Umac_Type,UMAC_MF+dir,0,dir,0,1,vel_time_slab);
   } //ilev=finest_level ... level

   int advect_iter_max=2;
   if (source_term==1) {
    advect_iter_max=1;
   } else if (source_term==0) {
    advect_iter_max=2;
   } else
    amrex::Error("advect_iter_max invalid");
  
   for (advect_iter=0;advect_iter<advect_iter_max;advect_iter++) {

    if (source_term==1) { 
     advect_time_slab=prev_time_slab;
    } else if (source_term==0) {
     if (advect_iter==0) {
      advect_time_slab=prev_time_slab;
     } else if (advect_iter==1) {
      advect_time_slab=cur_time_slab;
     } else
      amrex::Error("advect_iter invalid");
    } else
     amrex::Error("source_term invalid");

    for (int ilev=finest_level;ilev>=level;ilev--) {
     NavierStokes& ns_level=getLevel(ilev);

     ns_level.getStateDen_localMF(DEN_RECON_MF,1,advect_time_slab);
     ns_level.getState_localMF(VELADVECT_MF,1,0,
      num_materials_vel*AMREX_SPACEDIM,advect_time_slab); 
     for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
      ns_level.new_localMF(AMRSYNC_PRES_MF+dir,nfluxSEM,0,dir);
      ns_level.setVal_localMF(AMRSYNC_PRES_MF+dir,1.0e+40,0,nfluxSEM,0);
      ns_level.new_localMF(CONSERVE_FLUXES_MF+dir,nfluxSEM,0,dir);
      ns_level.setVal_localMF(CONSERVE_FLUXES_MF+dir,1.0e+40,0,nfluxSEM,0);
      ns_level.new_localMF(COARSE_FINE_FLUX_MF+dir,nfluxSEM,0,dir);
      ns_level.setVal_localMF(COARSE_FINE_FLUX_MF+dir,1.0e+40,0,nfluxSEM,0);
     } // dir=0..sdim-1
     ns_level.resize_levelsetLO(2,LEVELPC_MF);
     ns_level.VOF_Recon_resize(1,SLOPE_RECON_MF);
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
     int spectral_override=1;
     ns_level.avgDownEdge_localMF(CONSERVE_FLUXES_MF,scomp,ncomp_edge,
       0,AMREX_SPACEDIM,spectral_override,3);
    } // ilev=finest_level-1 ... level

    init_fluxes=0;
    int spectral_loop=0;
    int tileloop=0;
    for (int ilev=finest_level;ilev>=level;ilev--) {
     NavierStokes& ns_level=getLevel(ilev);
     ns_level.SEM_scalar_advection(init_fluxes,source_term,
        spectral_loop,tileloop);
    } // ilev=finest_level ... level

    avgDownALL(State_Type,0,num_materials_vel*(AMREX_SPACEDIM+1),1);
    int scomp_den=num_materials_vel*(AMREX_SPACEDIM+1);
    avgDownALL(State_Type,scomp_den,num_state_material*nmat,1);

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
   amrex::Error("enable_spectral or ns_time_order invalid");

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

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 int nsolve=1;
 int nsolveMM=nsolve*num_materials_vel;
 if (nsolveMM!=1)
  amrex::Error("nsolveMM invalid 1486");

 allocate_physics_vars();  

 if (localMF[MDOT_MF]->nComp()!=nsolve)
  amrex::Error("localMF[MDOT_MF]->nComp() invalid");

   //val,scomp,ncomp,ngrow
 setVal_localMF(MDOT_MF,0.0,0,nsolve,0); 

} // subroutine prelim_alloc

void NavierStokes::advance_MAC_velocity(int project_option) {

 int interp_option=0;
 int idx_velcell=-1;
 Real beta=0.0;

 if (num_materials_vel==1) {
  // do nothing
 } else
  amrex::Error("num_materials_vel invalid");

 if (face_flag==0) {
  beta=0.0;
  interp_option=0;  // unew^{f} = unew^{c->f}
 } else if (face_flag==1) {
  beta=0.0;
  // interp_option==4:
  // unew^{f}=
  // (i) unew^{f} in incompressible non-solid regions
  // (ii) u^{f,save} + (unew^{c}-u^{c,save})^{c->f} in spectral 
  //      regions or compressible regions.
  //      (u^{c,save} = *localMF[ADVECT_REGISTER_MF])
  //      (u^{f,save} = *localMF[ADVECT_REGISTER_FACE_MF+dir])
  // (iii) usolid in solid regions
  interp_option=4;  
 } else
  amrex::Error("face_flag invalid 9");

 Vector<blobclass> blobdata;

 increment_face_velocityALL(
   interp_option,project_option,
   idx_velcell,beta,blobdata);

} // subroutine advance_MAC_velocity()

// called from: NavierStokes::advance
void NavierStokes::do_the_advance(Real timeSEM,Real dtSEM,
  int& advance_status) {

 if (ParallelDescriptor::IOProcessor()) 
  std::cout << "do_the_advance timeSEM= " << timeSEM << 
   " dtSEM= " << dtSEM << '\n';

 very_last_sweep=0;

 advance_status=1; // 1=success 0=failure

 int post_restart_flag=0; 
 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

 int finest_level = parent->finestLevel();
 const int max_level = parent->maxLevel();

 int tessellate=1;

 if (level>0) 
  amrex::Error("level should equal zero in do_the_advance");
 if (finest_level>max_level)
  amrex::Error("max_level or finest_level invalid");
 
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
  if (bfact_grid<4)
   amrex::Error("we must have blocking factor at least 4");
  ns_level.check_grid_places();
 } // ilev=level ... finest_level

 double after_init = ParallelDescriptor::second();
 if ((verbose>0)||(show_timings==1)) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "elapsed time in preliminary allocation " << after_init-
        start_advance << '\n';
  }
 }

 int nsteps=parent->levelSteps(0); // nsteps==0 very first step.

 SDC_outer_sweeps=0;
 slab_step=0;
 SDC_setup_step();

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.allocate_SDC();
 }

 Number_CellsALL(real_number_of_cells);

 int SDC_outer_sweeps_end=ns_time_order;
 if (lower_slab_time<0.0)
  amrex::Error("lower_slab_time invalid");
 if (lower_slab_time==0.0)
  SDC_outer_sweeps_end=1;

 for (SDC_outer_sweeps=0;
      ((SDC_outer_sweeps<SDC_outer_sweeps_end)&&
       (advance_status==1));
      SDC_outer_sweeps++) {

  slab_step=0;
  SDC_setup_step();

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

   if (SDC_outer_sweeps==0)
    slab_step_start--; 
   slab_step_end++; 

  } else if (SDC_outer_sweeps_end==1) {

   // do nothing

  } else
   amrex::Error("SDC_outer_sweeps_end invalid");

  for (slab_step=slab_step_start;
       ((slab_step<=slab_step_end)&&(advance_status==1));
       slab_step++) {

   SDC_setup_step();

   int local_num_divu_outer_sweeps=num_divu_outer_sweeps;
   if ((slab_step==-1)||(slab_step==ns_time_order))
    local_num_divu_outer_sweeps=1;
   else if ((slab_step>=0)&&(slab_step<ns_time_order))
    local_num_divu_outer_sweeps=num_divu_outer_sweeps;
   else
    amrex::Error("slab_step invalid");

   for (int ilev=finest_level;ilev>=level;ilev--) {
    NavierStokes& ns_level=getLevel(ilev);
    int alloc_flag=1; // allocate and init to 0.0
    ns_level.alloc_gradp_over_rho(alloc_flag);
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
     int scomp_plot=AMREX_SPACEDIM+1;
     int ncomp_plot=num_state_material*num_materials;
     tecplot_debug(snewfab,xlo,fablo,fabhi,dx,-1,0,
       scomp_plot,ncomp_plot,interior_only);
    }

      // initialize "law of the wall" velocity derived from solid velocity.
      //  or
      // initialize "GNBC" velocity.
      // in: NavierStokes::do_the_advance (prior to nonlinear_advect)
    init_FSI_GHOST_MAC_MF_ALL(4);

    int SEM_VISCOUS_SANITY_CHECK=0;

    if (SEM_VISCOUS_SANITY_CHECK==1) {
     amrex::Warning("SEM_VISCOUS_SANITY_CHECK==1");

     int save_enable_spectral=enable_spectral;
     override_enable_spectral(viscous_enable_spectral);

     int update_placeholder=0;
     int project_option_placeholder=3;
     Vector<int> scomp;  
     Vector<int> ncomp;  
     int ncomp_check;
     int state_index;
     get_mm_scomp_solver(
       num_materials_vel,
       project_option_placeholder,
       state_index,
       scomp,
       ncomp,
       ncomp_check);
     int nsolve=AMREX_SPACEDIM;
     int nsolveMM=nsolve*num_materials_vel;
     if (state_index!=State_Type)
      amrex::Error("state_index invalid");
     if (ncomp_check!=nsolveMM)
      amrex::Error("ncomp_check invalid");

      // data at time = cur_time_slab
     getState_localMF_listALL(
       REGISTER_MARK_MF,1,
       state_index,
       scomp,
       ncomp);

     update_SEM_forcesALL(project_option_placeholder,REGISTER_MARK_MF,
       update_placeholder,update_placeholder);

     override_enable_spectral(save_enable_spectral);
    } else if (SEM_VISCOUS_SANITY_CHECK==0) {
     // do nothing
    } else
     amrex::Error("SEM_VISCOUS_SANITY_CHECK invalid");

    int mass_transfer_active=0;

     // 1. ADVECTION (both Eulerian and Lagrangian materials)
     // 2. IF MDOT <> 0 in previous time step, then there
     //    is expansion/compression
    if ((slab_step>=0)&&(slab_step<ns_time_order)) {

     if (disable_advection==0) {

      nonlinear_advection();

     } else if (disable_advection==1) {
      
      if (face_flag==1) {
       // delete_advect_vars() called in NavierStokes::do_the_advance
       // right after increment_face_velocityALL. 
       for (int ilev=finest_level;ilev>=level;ilev--) {
        NavierStokes& ns_level=getLevel(ilev);
        // initialize ADVECT_REGISTER_FACE_MF and
        //            ADVECT_REGISTER_MF
        ns_level.prepare_advect_vars(prev_time_slab);
       }
      } else if (face_flag==0) {
       // do nothing
      } else
       amrex::Error("face_flag invalid");

     } else
      amrex::Error("disable_advection invalid");

    } else if ((slab_step==-1)||
               (slab_step==ns_time_order)) {
     // do nothing
    } else
     amrex::Error("slab_step invalid");


     // in: NavierStokes::do_the_advance
    allocate_levelsetLO_ALL(1,LEVELPC_MF);

    if ((ns_time_order==1)&&
        (enable_spectral==3))
     amrex::Error("(ns_time_order==1)&&(enable_spectral==3)");
    if ((ns_time_order>=2)&&
        ((enable_spectral==0)||
         (enable_spectral==2)))
     amrex::Error("(ns_time_order>=2)&&(enable_spectral==0 or 2)");

    if ((enable_spectral==1)||  // space-time SEM
	(enable_spectral==2)||  // SEM space
        (ns_time_order>=2)) {   // SEM time

     if (disable_advection==0) {

      double start_SEMADV_time=ParallelDescriptor::second();

      int source_term=1;
      if ((slab_step>=0)&&(slab_step<=ns_time_order)) {

        // SEM_advectALL starts off by using the prev_time_slab data.
       if ((slab_step==0)&&
           (SDC_outer_sweeps>0)&&
           (SDC_outer_sweeps<ns_time_order)) {
        // do nothing: F(t^n) already init.
       } else if ((slab_step==0)&&
                  (SDC_outer_sweeps==0)) {
        SEM_advectALL(source_term);
       } else if ((slab_step>0)&&
                  (SDC_outer_sweeps>=0)&&
                  (SDC_outer_sweeps<ns_time_order)) {
        SEM_advectALL(source_term);
       } else
        amrex::Error("slab_step or SDC_outer_sweeps invalid");

      } // ((slab_step>=0)&&(slab_step<=ns_time_order))

       // SEM_advectALL starts off by using the prev_time_slab data.
      source_term=0;
      if ((slab_step>=0)&&(slab_step<ns_time_order)) {
       SEM_advectALL(source_term);
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

    } else if (((enable_spectral==0)||
  	        (enable_spectral==3))&&
               (ns_time_order==1)) {
     // do nothing
    } else
     amrex::Error("enable_spectral or ns_time_order invalid do the advance");

      // in: NaveriStokes::do_the_advance
    int project_option=0;

    if (is_zalesak()) {
      project_option=1;  // initial project
    } else {
      project_option=0;
    }

      // 0=do not update the error 1=update the error
    int update_flag=0; 

    debug_memory();

    mass_transfer_active=0;

    if ((slab_step>=0)&&(slab_step<ns_time_order)) {
     if (is_phasechange==1) {
      mass_transfer_active=1;
     } else if (is_phasechange==0) {
      mass_transfer_active=0;
     } else
      amrex::Error("is_phasechange invalid");
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
    int coarsest_level=0;

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

	// CREATE SEEDS, NUCLEATION.
       for (int ilev=level;ilev<=finest_level;ilev++) {
        int nucleation_flag=1;
	color_count=1; // filler
        NavierStokes& ns_level=getLevel(ilev);
        ns_level.level_phase_change_rate(blobdata,color_count,
          nucleation_flag);
       }

       int at_least_one_ice=0;
       for (int im=0;im<nmat;im++) {
        if (is_ice_matC(im)==1) {
         at_least_one_ice=1;
        } else if (is_ice_matC(im)==0) {
         // do nothing
        } else
         amrex::Error("is_ice_matC invalid");
       } // im=0..nmat-1 

       Vector<int> recalesce_material;
       recalesce_material.resize(nmat);
       int at_least_one=0;
       for (int im=1;im<=nmat;im++) {
        recalesce_material[im-1]=parent->AMR_recalesce_flag(im);
        if (parent->AMR_recalesce_flag(im)>0) {
         if (at_least_one_ice!=1)
          amrex::Error("expecting at least one material FSI_flag==3 or 6");
         at_least_one=1;
        } 
       } //im=1..nmat
       if (at_least_one==1) {
        Vector<Real> recalesce_state_old;
        Vector<Real> recalesce_state_new;
        int recalesce_num_state=6;
        recalesce_state_old.resize(recalesce_num_state*nmat);
        recalesce_state_new.resize(recalesce_num_state*nmat);
        parent->recalesce_get_state(recalesce_state_old,nmat);

        FORT_INITRECALESCE(
         recalesce_material.dataPtr(),
         recalesce_state_old.dataPtr(),
         &recalesce_num_state,&nmat);
   
        process_recalesce_dataALL(recalesce_material,
         recalesce_state_old,recalesce_state_new);
        parent->recalesce_put_state(recalesce_state_new,nmat);
       } else if (at_least_one==0) {
        // do nothing
       } else
        amrex::Error("at_least_one invalid");

       delta_mass.resize(thread_class::nthreads);
       for (int tid=0;tid<thread_class::nthreads;tid++) {
        delta_mass[tid].resize(2*nmat); // source 1..nmat  dest 1..nmat
       for (int im=0;im<2*nmat;im++)
        delta_mass[tid][im]=0.0;
       } // tid

       if (1==0) {
        int basestep_debug=nStep();
        parent->writeDEBUG_PlotFile(basestep_debug,SDC_outer_sweeps,slab_step);
        std::cout << "press any number then enter: before nucleate_bubbles\n";
        int n_input;
        std::cin >> n_input;
       }

       ParallelDescriptor::Barrier();

       tessellate=1;
       int idx_mdot=-1;
       int operation_flag=0;

       ColorSumALL( 
         operation_flag, //=0
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
        ns_level.avgDown(LS_Type,0,nmat,0);
        ns_level.MOFavgDown();
        int scomp_den=num_materials_vel*(AMREX_SPACEDIM+1);
        ns_level.avgDown(State_Type,scomp_den,num_state_material*nmat,1);
       }  // ilev=finest_level ... level  

       if (1==0) {
        int basestep_debug=nStep();
        parent->writeDEBUG_PlotFile(basestep_debug,SDC_outer_sweeps,slab_step);
        std::cout << "press any number then enter: after nucleate_bubbles\n";
        int n_input;
        std::cin >> n_input;
       }

       if (verbose>0) {
        if (ParallelDescriptor::IOProcessor()) {
         for (int im=0;im<nmat;im++) {
          std::cout << "Nucleation stats: im,source,dest " << im << ' ' <<
           delta_mass[0][im] << ' ' << delta_mass[0][im+nmat] << '\n';
         }
        }
       } 

       // generates SLOPE_RECON_MF
       update_flag=0; // do not update the error indicator
       int init_vof_ls_prev_time=0;
        // Fluids tessellate; solids overlay.
       VOF_Recon_ALL(1,cur_time_slab,update_flag,init_vof_ls_prev_time,
        SLOPE_RECON_MF);
       int keep_all_interfaces=1;
       makeStateDistALL(keep_all_interfaces);

       make_physics_varsALL(project_option,post_restart_flag,5); 

        // if face_flag==0: unew^{f} = unew^{c->f}
	// if face_flag==1: 
	//  unew^{f}=
        // (i) unew^{f} in incompressible non-solid regions
        // (ii) u^{f,save} + (unew^{c}-u^{c,save})^{c->f} in spectral 
        //      regions or compressible regions.
        //      (u^{c,save} = *localMF[ADVECT_REGISTER_MF])
        //      (u^{f,save} = *localMF[ADVECT_REGISTER_FACE_MF+dir])
        // (iii) usolid in solid regions

       advance_MAC_velocity(project_option);
  
      } else if (mass_transfer_active==0) {

       update_flag=1;  // update the error in S_new
       int init_vof_ls_prev_time=0;
       VOF_Recon_ALL(1,cur_time_slab,update_flag,init_vof_ls_prev_time,
        SLOPE_RECON_MF);
       int keep_all_interfaces=0;
       makeStateDistALL(keep_all_interfaces);
      } else
       amrex::Error("mass_transfer_active invalid");

      for (int ilev=finest_level;ilev>=level;ilev--) {
       NavierStokes& ns_level=getLevel(ilev);
       int project_option_combine=2; // temperature in do_the_advance
       int combine_flag=2; // only update if vfrac<VOFTOL
       int hflag=0;
       // combine_idx==-1 => update S_new  
       // combine_idx>=0  => update localMF[combine_idx]
       int combine_idx=-1;  
       int update_flux=0;
       int interface_cond_avail=0;

       ns_level.combine_state_variable(
        project_option_combine,
        combine_idx,
        combine_flag,
        hflag,
        update_flux,
        interface_cond_avail); 
       for (int ns=0;ns<num_species_var;ns++) {
        project_option_combine=100+ns; // species in do_the_advance
        ns_level.combine_state_variable(
         project_option_combine,
         combine_idx,
         combine_flag,
         hflag,
         update_flux,
         interface_cond_avail); 
       }

      } // ilev=finest_level ... level

    } else if ((slab_step==-1)||(slab_step==ns_time_order)) {

      update_flag=0;  // do not update the error in S_new
      int init_vof_ls_prev_time=0;
      VOF_Recon_ALL(1,cur_time_slab,update_flag,init_vof_ls_prev_time,
        SLOPE_RECON_MF);

    } else
      amrex::Error("slab_step invalid");


      // velocity and pressure
    avgDownALL(State_Type,0,num_materials_vel*(AMREX_SPACEDIM+1),1);

    int scomp_den=num_materials_vel*(AMREX_SPACEDIM+1);

     // den,denA,(total E)/rho,temp,pres,...
    avgDownALL(State_Type,scomp_den,num_state_material*nmat,1);  
    debug_memory();

    double start_phys_time=ParallelDescriptor::second();

    if ((slab_step>=0)&&(slab_step<ns_time_order)) {

      if (mass_transfer_active==1) {

       if (ngrow_expansion!=2)
        amrex::Error("ngrow_expansion!=2");
       if (ngrow_make_distance!=3)
        amrex::Error("ngrow_make_distance!=3");

        // first nten components correspond to the status.
       int ncomp_per_burning=AMREX_SPACEDIM;
       int ncomp_per_tsat=2;
       int nburning=nten*(ncomp_per_burning+1);
       int ntsat=nten*(ncomp_per_tsat+1);

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
          ngrow_make_distance,-1);
        ns_level.setVal_localMF(BURNING_VELOCITY_MF,0.0,0,
          nburning,ngrow_make_distance);

	int n_normal=(nmat+nten)*(AMREX_SPACEDIM+1);

        ns_level.new_localMF(FD_NRM_ND_MF,n_normal,
          ngrow_make_distance+1,-1);
        ns_level.setVal_localMF(FD_NRM_ND_MF,0.0,0,
          n_normal,ngrow_make_distance+1);

        ns_level.new_localMF(FD_CURV_CELL_MF,2*(nmat+nten),
          ngrow_make_distance,-1);
        ns_level.setVal_localMF(FD_CURV_CELL_MF,0.0,0,
          2*(nmat+nten),ngrow_make_distance);

        ns_level.new_localMF(SATURATION_TEMP_MF,ntsat,
          ngrow_make_distance,-1);
        ns_level.setVal_localMF(SATURATION_TEMP_MF,0.0,0,
          ntsat,ngrow_make_distance);

        ns_level.new_localMF(JUMP_STRENGTH_MF,2*nten,ngrow_expansion,-1); 
        ns_level.setVal_localMF(JUMP_STRENGTH_MF,0.0,0,2*nten,ngrow_expansion);

       } // ilev=level ... finest_level

       debug_ngrow(JUMP_STRENGTH_MF,ngrow_expansion,30);
       debug_ngrow(SWEPT_CROSSING_MF,0,31);
       debug_ngrow(BURNING_VELOCITY_MF,ngrow_make_distance,31);
       debug_ngrow(SATURATION_TEMP_MF,ngrow_make_distance,31);
       debug_ngrow(FD_NRM_ND_MF,ngrow_make_distance+1,31);
       debug_ngrow(FD_CURV_CELL_MF,ngrow_make_distance,31);

       for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
        debug_ngrow(AREA_MF+dir,1,355);
        debug_ngrow(FACE_VAR_MF+dir,0,355);
       }
       debug_ngrow(MDOT_MF,0,355);

       zeroALL(ngrow_make_distance,nburning,BURNING_VELOCITY_MF);
       zeroALL(ngrow_expansion,2*nten,JUMP_STRENGTH_MF);
        //ngrow,ncomp,val,dest_mf
       setVal_array(0,nmat,1.0,SWEPT_CROSSING_MF);
        // piecewise constant interpolation at coarse/fine borders.
        // fluid LS can be positive in the solid regions.
        // HOLD_LS_DATA_MF is deleted in phase_change_redistributeALL()
       allocate_levelsetLO_ALL(ngrow_distance,HOLD_LS_DATA_MF);
       if (localMF[HOLD_LS_DATA_MF]->nComp()!=nmat*(AMREX_SPACEDIM+1))
        amrex::Error("localMF[idx]->nComp()!=nmat*(AMREX_SPACEDIM+1)");
       debug_ngrow(HOLD_LS_DATA_MF,normal_probe_size+3,30);

       for (int ilev=level;ilev<=finest_level;ilev++) {
        NavierStokes& ns_level=getLevel(ilev);
        ns_level.new_localMF(LS_NRM_FD_MF,nmat*AMREX_SPACEDIM,1,-1);
	ns_level.build_NRM_FD_MF(LS_NRM_FD_MF,HOLD_LS_DATA_MF,1);
       }

       FIX ME SMOOTH TEMP HERE

        // BURNING_VELOCITY_MF flag==+ or - 1 if valid rate of phase change.
       for (int ilev=level;ilev<=finest_level;ilev++) {
        int nucleation_flag=0;
        NavierStokes& ns_level=getLevel(ilev);
        ns_level.level_phase_change_rate(blobdata,color_count,
          nucleation_flag);
       }

       delete_array(TEMPERATURE_SMOOTH_MF);
       delete_array(TYPE_MF);
       delete_array(COLOR_MF);

       for (int ilev=finest_level;ilev>=level;ilev--) {
        NavierStokes& ns_level=getLevel(ilev);
	  // in: NavierStokes2.cpp
        ns_level.avgDownBURNING_localMF(
			BURNING_VELOCITY_MF,
			SATURATION_TEMP_MF);
        ns_level.avgDown(LS_Type,0,nmat,0);
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
	 // calls FORT_EXTEND_BURNING_VEL
       for (int ilev=level;ilev<=finest_level;ilev++) {
        NavierStokes& ns_level=getLevel(ilev);
        ns_level.level_phase_change_rate_extend();
       }

       if (visual_phase_change_plot_int>0) {
        if (very_last_sweep==1) {
         int ratio=(nsteps+1)/visual_phase_change_plot_int;
	 ratio=ratio*visual_phase_change_plot_int;
	 if (ratio==nsteps+1) {

	  int caller_id=1;
           //TY_GAMMA<stuff>.plt (visit can open binary tecplot files)
          writeSanityCheckData(
           "TY_GAMMA",
           "SATURATION_TEMP_MF: flag12,flag13,flag23,T_GAMMA12,Y_GAMMA12, ...",
           caller_id,
           localMF[SATURATION_TEMP_MF]->nComp(), 
           SATURATION_TEMP_MF,
           -1,  // State_Type==-1 
           -1); // data_dir==-1 (cell centered)

	  caller_id=2;
           //BURNVEL<stuff>.plt (visit can open binary tecplot files)
          writeSanityCheckData(
           "BURNVEL",
           "BURNING_VELOCITY_MF: flag12,flag13,flag23,[xyz]V12,[xyz]V13, ..",
           caller_id,
           localMF[BURNING_VELOCITY_MF]->nComp(), 
           BURNING_VELOCITY_MF,
           -1,  // State_Type==-1 
           -1); // data_dir==-1 (cell centered)

	 }
	} else if (very_last_sweep==0) {
	 // do nothing
	} else
	 amrex::Error("very_last_sweep invalid");

       } else if (visual_phase_change_plot_int==0) {
        // do nothing
       } else
        amrex::Error("visual_phase_change_plot_int invalid");

       allocate_array(1,2*nten*AMREX_SPACEDIM,-1,nodevel_MF);
       setVal_array(1,2*nten*AMREX_SPACEDIM,0.0,nodevel_MF);

       delta_mass.resize(thread_class::nthreads);
       for (int tid=0;tid<thread_class::nthreads;tid++) {
        delta_mass[tid].resize(2*nmat); // source 1..nmat  dest 1..nmat
        for (int im=0;im<2*nmat;im++)
         delta_mass[tid][im]=0.0;
       } // tid

        // in: do_the_advance
       for (int ilev=finest_level;ilev>=level;ilev--) {
        NavierStokes& ns_level=getLevel(ilev);
        ns_level.getStateDen_localMF(DEN_RECON_MF,1,cur_time_slab);
       }

       // 1.initialize node velocity from BURNING_VELOCITY_MF
       // 2.unsplit advection of materials changing phase
       // 3.update volume fractions, jump strength, temperature
       level_phase_change_convertALL();

       delete_array(LS_NRM_FD_MF);
       delete_array(BURNING_VELOCITY_MF);
       delete_array(FD_NRM_ND_MF);
       delete_array(FD_CURV_CELL_MF);

       if (verbose>0) {
        if (ParallelDescriptor::IOProcessor()) {
         for (int im=0;im<nmat;im++) {
          std::cout << "convert statistics: im,source,dest " << im << ' ' <<
           delta_mass[0][im] << ' ' << delta_mass[0][im+nmat] << '\n';
         }
        }
       }

        // in: do_the_advance
       delete_array(DEN_RECON_MF);
       delete_array(nodevel_MF);

       update_flag=1;  // update the error in S_new
       int init_vof_ls_prev_time=0;
       VOF_Recon_ALL(1,cur_time_slab,update_flag,init_vof_ls_prev_time,
        SLOPE_RECON_MF);

        // in: do_the_advance
        // 1. prescribe solid temperature, velocity, and geometry where
        //    appropriate.
        // 2. extend level set functions into the solid.
       int renormalize_only=0;
       int local_truncate=0;
       prescribe_solid_geometryALL(cur_time_slab,renormalize_only,
        local_truncate);

       int keep_all_interfaces=0;
       makeStateDistALL(keep_all_interfaces);

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
     parent->writeDEBUG_PlotFile(basestep_debug,SDC_outer_sweeps,slab_step);
     std::cout << "press any number then enter: before make_physics_varsALL\n";
     int n_input;
     std::cin >> n_input;
    }

     // initialize "law of the wall" velocity derived from solid velocity.
     //  or
     // initialize "GNBC" velocity.
     // in: NavierStokes::do_the_advance (prior to viscous force step, and
     //  after reinitialization)
    init_FSI_GHOST_MAC_MF_ALL(4);

// At this stage, variables are not scaled, so facevel_index will have
// to be scaled later.
    debug_memory();
    make_physics_varsALL(project_option,post_restart_flag,6); 

    if (1==0) {
     int basestep_debug=nStep();
     parent->writeDEBUG_PlotFile(basestep_debug,SDC_outer_sweeps,slab_step);
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

       // if face_flag==0: unew^{f} = unew^{c->f}
       // if face_flag==1: 
       //  unew^{f}=
       // (i) unew^{f} in incompressible non-solid regions
       // (ii) u^{f,save} + (unew^{c}-u^{c,save})^{c->f} in spectral 
       //      regions or compressible regions.
       //      (u^{c,save} = *localMF[ADVECT_REGISTER_MF])
       //      (u^{f,save} = *localMF[ADVECT_REGISTER_FACE_MF+dir])
       // (iii) usolid in solid regions
      advance_MAC_velocity(project_option);

      if (face_flag==1) {
       for (int ilev=finest_level;ilev>=level;ilev--) {
        NavierStokes& ns_level=getLevel(ilev);
	 // delete ADVECT_REGISTER_FACE_MF and ADVECT_REGISTER_MF
        ns_level.delete_advect_vars();
       } // ilev=finest_level ... level

       for (int ilev=finest_level;ilev>=level;ilev--) {
        NavierStokes& ns_level=getLevel(ilev);
         //if temperature_primitive_var==0,
         // add beta * (1/cv) * (u dot u/2) to temp
        Real beta=1.0;
        ns_level.increment_KE(beta); 
	int use_VOF_weight=1;
        ns_level.VELMAC_TO_CELL(use_VOF_weight);
        beta=-1.0;
        ns_level.increment_KE(beta);
       }
 
      } else if (face_flag==0) {
       // do nothing
      } else
       amrex::Error("face_flag invalid");

      for (int ilev=finest_level;ilev>=level;ilev--) {
       NavierStokes& ns_level=getLevel(ilev);
       int project_option_combine=2; // temperature in do_the_advance
       int combine_flag=2; // update F_m=0 cells only.
       int hflag=0;
       // combine_idx==-1 => update S_new  
       // combine_idx>=0  => update localMF[combine_idx]
       int combine_idx=-1;  
       int update_flux=0;
       int interface_cond_avail=0;

       ns_level.combine_state_variable(
        project_option_combine,
        combine_idx,
        combine_flag,
        hflag,
        update_flux,
        interface_cond_avail); 
       for (int ns=0;ns<num_species_var;ns++) {
        project_option_combine=100+ns;
        ns_level.combine_state_variable(
         project_option_combine,
         combine_idx,
         combine_flag,
         hflag,
         update_flux,
         interface_cond_avail); 
       }
       project_option_combine=3;  // velocity
       ns_level.combine_state_variable(
        project_option_combine,
        combine_idx,
        combine_flag,
        hflag,
        update_flux,
        interface_cond_avail);
       project_option_combine=0; // mac velocity
       update_flux=1;
       ns_level.combine_state_variable(
        project_option_combine,
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

       // 1. modifies localMF[JUMP_STRENGTH_MF] (size 2 * nten)
       //    a) FORT_TAGEXPANSION
       //    b) FORT_DISTRIBUTEEXPANSION
       //    c) FORT_CLEAREXPANSION
       //    d) FORT_INITJUMPTERM ( modifies localMF[MDOT_MF] )
       //    e) FORT_GETCOLORSUM (twice)
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

      // 3. TENSOR ADVECTION
      //  (non-Newtonian materials)
      // second half of D^{upside down triangle}/Dt
      tensor_advection_updateALL();

      if (is_zalesak()==1) {

       for (int ilev=finest_level;ilev>=level;ilev--) {
        NavierStokes& ns_level=getLevel(ilev);
        // velocity scaled by global_velocity_scale
        // init cell center gas/liquid velocity 
        ns_level.zalesakVEL();  
       } // ilev=finest_level ... level

        // unew^{f} = unew^{c->f}
       int interp_option=0;
       int idx_velcell=-1;

       if (num_materials_vel==1) {
        // do nothing
       } else
        amrex::Error("num_materials_vel invalid");

       Real beta=0.0;
       Vector<blobclass> local_blobdata;
       increment_face_velocityALL(
         interp_option,project_option,
         idx_velcell,beta,local_blobdata);

      } else if (is_zalesak()==0) {

       for (int ilev=finest_level;ilev>=level;ilev--) {
        NavierStokes& ns_level=getLevel(ilev);
        int project_option_combine=3; // velocity in do_the_advance
        int combine_flag=2; // combine if vfrac<VOFTOL
        int hflag=0; // inhomogeneous option
        int combine_idx=-1;  // update state variables
        int update_flux=0;
        int interface_cond_avail=0;

	  // declared in: Diffusion.cpp
        ns_level.combine_state_variable(
         project_option_combine,
         combine_idx,
         combine_flag,
         hflag,
         update_flux, 
         interface_cond_avail);
        project_option_combine=0; // mac velocity
        update_flux=1;
        ns_level.combine_state_variable(
         project_option_combine,
         combine_idx,
         combine_flag,
         hflag,
         update_flux,
         interface_cond_avail);
       } // ilev = finest_level ... level

        // T_advect_MF=new temperature
       int alloc_flag=1;
       alloc_DTDtALL(alloc_flag);

       // 4. Backwards Euler building block: VISCOSITY, thermal diffusion,
       //    species diffusion, conservative surface tension force.
       //   a. hoop stress
       //   b. boussinesq approximation
       //   c. coriolis effect
       //   d. viscous force
       //   e. viscoelastic force
       //   f. FSI force
       //   g. momentum force
       //   h. Marangoni force and conservative surface tension force
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

       if (1==0) {
          // S_new is level 0 data
        MultiFab& S_new=get_new_data(State_Type,slab_step+1);
	int caller_id=1;
	 // data file name "VISCSOLVE<stuff>.plt"
	 // after the viscous solve, but before the pressure projection.
	 // cell data in the fluid, next to the solid, should "make sense"
	 // xvel,yvel,zvel,pressure,(density, temperature) x nmat,
	 // (VFRAC,centroid) x nmat, error indicator
        writeSanityCheckData(
         "VISCSOLVE",
         "in: NavierStokes::do_the_advance, State_Type after veldiffuseALL", 
         caller_id,
         S_new.nComp(),
         -1, // data_mf==-1
         State_Type,
         -1); // data_dir==-1
       }

       if (1==0) {
        int basestep_debug=nStep();
        parent->writeDEBUG_PlotFile(basestep_debug,SDC_outer_sweeps,slab_step);
       }

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
         } else if (((material_type[im_low]>0)&&
                     (material_type[im_low]<999))||
                    (material_type_lowmach[im_low]==0)) {
          // do nothing
         } else
          amrex::Error("material_type or material_type_lowmach invalid");
        } else
         amrex::Error("ns_is_rigid invalid");
       }  //im_low=0..nmat-1

       if (is_any_lowmach==1) {

        Vector<blobclass> local_blobdata;
        Vector< Vector<Real> > local_mdot_data;
        Vector< Vector<Real> > local_mdot_comp_data;
        Vector< Vector<Real> > local_mdot_data_redistribute;
        Vector< Vector<Real> > local_mdot_comp_data_redistribute;
        Vector<int> local_type_flag;

        int local_color_count=0;
        int local_coarsest_level=0;
        int idx_mdot=-1;
        int local_tessellate=3;
        int operation_flag=0; // allocate TYPE_MF,COLOR_MF

         // for each blob, find sum_{F>=1/2} pressure * vol and
	// sum_{F>=1/2} vol.
        ColorSumALL(
         operation_flag, // =0
         local_tessellate, //=3
         local_coarsest_level,
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
         local_coarsest_level,
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

       double intermediate_time = ParallelDescriptor::second();

        // 5. PRESSURE PROJECTION 
        //    a. gravity
        //    b. surface tension
        //    c. pressure gradient
       if (disable_pressure_solve==0) {
 
	 // FSI_flag=3,6 (ice) or FSI_flag=5 (FSI PROB.F90 rigid material)
        if (FSI_material_exists()==1) {

          // MDOT term included
         int rigid_project_option=0;
         multiphase_project(rigid_project_option);

          // MDOT term not included, instead 
          // If compressible: DIV_new=-dt(pnew-padv)/(rho c^2 dt^2)+
	  //                          MDOT_MF dt/vol
          // If incompressible: DIV_new=MDOT_MF dt/vol
          // If one of the adjoining cells of a face are in the 
          // "flexible solid," then the face coefficient = 0. 
          // See: FORT_BUILDFACEWT, FACE_VAR_MF
         rigid_project_option=11; // final project
         multiphase_project(rigid_project_option);

        } else if (FSI_material_exists()==0) {

         multiphase_project(project_option); // pressure

        } else
         amrex::Error("FSI_material_exists invalid");

        int singular_parts_exist=0;
        for (int im=0;im<nmat;im++) {
         if (is_singular_coeff(im)==0) {
          // do nothing
         } else if (is_singular_coeff(im)==1) {
          singular_parts_exist=1;
         } else
          amrex::Error("is_singular_coeff invalid");
        } // im=0..nmat-1
 
        if (singular_parts_exist==1) {
 
         if (extend_pressure_into_solid==1) {
          int project_option_extend=12;
          multiphase_project(project_option_extend);
         } else if (extend_pressure_into_solid==0) {
          // do nothing
         } else
          amrex::Error("extend_pressure_into_solid invalid");
 
        } else if (singular_parts_exist==0) {
         // do nothing
        } else
         amrex::Error("singular_parts_exist invalid");

       } else if (disable_pressure_solve==1) {
        // do nothing
       } else
        amrex::Error("disable_pressure_solve invalid");

        // in: do_the_advance
        // 1. prescribe solid temperature, velocity, and geometry where
        //    appropriate.
        // 2. extend level set functions into the solid.
       int renormalize_only=0;
       int local_truncate=0;
       prescribe_solid_geometryALL(cur_time_slab,renormalize_only,
        local_truncate);

       for (int ilev=finest_level;ilev>=level;ilev--) {
        NavierStokes& ns_level=getLevel(ilev);
        int project_option_combine=3;  // velocity in do_the_advance
        int combine_flag=2; // combine if vfrac<VOFTOL
        int hflag=0;
        int combine_idx=-1;  // update state variables
        int update_flux=0;
        ns_level.combine_state_variable(
         project_option_combine,
         combine_idx,
         combine_flag,
         hflag,
         update_flux,
         interface_cond_avail);
        project_option_combine=0; // mac velocity
        update_flux=1;
        ns_level.combine_state_variable(
         project_option_combine, 
	 combine_idx,
         combine_flag,
         hflag,
         update_flux,
         interface_cond_avail);
       } // ilev 
       debug_memory();

       double end_pressure_solve = ParallelDescriptor::second();

       if ((verbose>0)||(show_timings==1)) {
        if (ParallelDescriptor::IOProcessor()) {
         std::cout << "pressure solve time " << end_pressure_solve-
            intermediate_time << '\n';
         std::cout << "number of cells in the pressure solve " <<
            real_number_of_cells << '\n';
        }
       }

       for (int ilev=finest_level;ilev>=level;ilev--) {
        NavierStokes& ns_level=getLevel(ilev);
        int local_alloc_flag=5; // subtract unew from dt_gradp_over_rho
        ns_level.alloc_gradp_over_rho(local_alloc_flag);
       }

      } else
       amrex::Error("is_zalesak invalid");

      if (mass_transfer_active==1) {
       delete_array(SATURATION_TEMP_MF);
      } else if (mass_transfer_active==0) {
       // do nothing
      } else
       amrex::Error("mass_transfer_active invalid");

      debug_memory();

      if (very_last_sweep==0) {

       Real local_fixed_dt;
       if (nsteps==0) {
        local_fixed_dt=fixed_dt_init;
       } else if (nsteps>0) {
        local_fixed_dt=fixed_dt;
       } else {
        local_fixed_dt=0.0;
        amrex::Error("nsteps invalid");
       }

       Real dt_predict=estTimeStep(local_fixed_dt);
       Real dt_predict_max=dt_predict;
       Real dt_predict_min=dt_predict;
       ParallelDescriptor::ReduceRealMax(dt_predict_max);
       ParallelDescriptor::ReduceRealMin(dt_predict_min);
       Real dt_error=dt_predict_max-dt_predict_min;
       if ((dt_error>=0.0)&&
           (dt_predict_min>0.0)&&
           (dt_predict_max>0.0)) {
	Real dt_tol=dt_predict_max*1.0e-13;
        if (dt_predict_min<1.0)
         dt_tol=1.0e-13;
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

      if ((enable_spectral==1)||
	  (enable_spectral==3)) {
       // do nothing
      } else
       amrex::Error("enable_spectral invalid 2");

      if ((viscous_enable_spectral==0)||  // do low order viscosity in time.
          (viscous_enable_spectral==1)||  // SEM space and time
          (viscous_enable_spectral==3)) { // SEM time
       // do nothing
      } else
       amrex::Error("viscous_enable_spectral invalid");

      if ((projection_enable_spectral==1)||  // SEM space and time
          (projection_enable_spectral==3)) { // SEM time
       // do nothing
      } else
       amrex::Error("projection_enable_spectral invalid");

      if ((slab_step>=-1)&&(slab_step<ns_time_order)) {

       if (divu_outer_sweeps+1==local_num_divu_outer_sweeps) {

        int update_spectralF=1;

        // grad p and div(up)
        int update_stableF=1;
        if (slab_step==-1)
         update_stableF=0;

        Vector<int> scomp;  
        Vector<int> ncomp;  
        int ncomp_check;
        int state_index;
        int project_option_op=0;

        get_mm_scomp_solver(
         num_materials_vel,
         project_option_op,
         state_index,
         scomp,
         ncomp,
         ncomp_check);

        int nsolve=1;
        int nsolveMM=nsolve*num_materials_vel;

        if (state_index!=State_Type)
         amrex::Error("state_index invalid");
        if (ncomp_check!=nsolveMM)
         amrex::Error("ncomp_check invalid");

        int save_enable_spectral=enable_spectral;
        override_enable_spectral(projection_enable_spectral);

         // data at time = cur_time_slab
        getState_localMF_listALL(
          PRESPC2_MF,1,
          state_index,
          scomp,
          ncomp);

         // HOfab=grad p,  div(up)
         // calls: UPDATESEMFORCE in GODUNOV_3D.F90
        update_SEM_forcesALL(project_option_op,PRESPC2_MF,
         update_spectralF,update_stableF);

        override_enable_spectral(save_enable_spectral);
        // end: grad p and div(up)

        save_enable_spectral=enable_spectral;
        override_enable_spectral(viscous_enable_spectral);

         // velocity here to be used later.
        project_option_op=3; // viscosity project_option

        get_mm_scomp_solver(
         num_materials_vel,
         project_option_op,
         state_index,
         scomp,
         ncomp,
         ncomp_check);

        nsolve=AMREX_SPACEDIM;
        nsolveMM=nsolve*num_materials_vel;

        if (state_index!=State_Type)
         amrex::Error("state_index invalid");
        if (ncomp_check!=nsolveMM)
         amrex::Error("ncomp_check invalid");

          // data at time = cur_time_slab
        getState_localMF_listALL(
          REGISTER_MARK_MF,1,
          state_index,
          scomp,
          ncomp);

         // -div(k grad T)-THERMAL_FORCE_MF
        update_stableF=0;
        project_option_op=2; // temperature project_option

        get_mm_scomp_solver(
         num_materials_scalar_solve,
         project_option_op,
         state_index,
         scomp,
         ncomp,
         ncomp_check);

        nsolve=1;
        nsolveMM=nsolve*num_materials_scalar_solve;

        if (state_index!=State_Type)
         amrex::Error("state_index invalid");
        if (ncomp_check!=nsolveMM)
         amrex::Error("ncomp_check invalid");

          // data at time = cur_time_slab
        getState_localMF_listALL(
          BOUSSINESQ_TEMP_MF,1,
          state_index,
          scomp,
          ncomp);

        allocate_array(1,nsolveMM,-1,THERMAL_FORCE_MF);
        int update_state=0;
        thermal_transform_forceALL(REGISTER_MARK_MF,BOUSSINESQ_TEMP_MF,
         THERMAL_FORCE_MF,update_state);

        // HOfab=-div(k grad T)-THERMAL_FORCE_MF
        // calls: UPDATESEMFORCE in GODUNOV_3D.F90
        if ((viscous_enable_spectral==1)||  // SEM space and time
            (viscous_enable_spectral==3)) { // SEM time
         update_SEM_forcesALL(project_option_op,BOUSSINESQ_TEMP_MF,
          update_spectralF,update_stableF);
        } else if (viscous_enable_spectral==0) {
         // do nothing
        } else
         amrex::Error("viscous_enable_spectral invalid");

        override_enable_spectral(save_enable_spectral);

         // -momforce
        update_stableF=0;
        project_option_op=4;  // project_option for prescribed Mom. force

        nsolve=AMREX_SPACEDIM;
        nsolveMM=nsolve*num_materials_vel;

        allocate_array(1,nsolveMM,-1,NEG_MOM_FORCE_MF);
         // force at time = cur_time_slab
        update_state=0;
        mom_forceALL(NEG_MOM_FORCE_MF,update_state);

        // force at time = cur_time_slab
        // HOfab=NEG_MOM_FORCE_MF
        // calls: UPDATESEMFORCE in GODUNOV_3D.F90
        update_SEM_forcesALL(project_option_op,NEG_MOM_FORCE_MF,
         update_spectralF,update_stableF);

        save_enable_spectral=enable_spectral;
        override_enable_spectral(viscous_enable_spectral);

        // HOfab=-div(2 mu D)-HOOP_FORCE_MARK_MF
        // calls: UPDATESEMFORCE in GODUNOV_3D.F90
        update_stableF=0;
        project_option_op=3;  // viscosity project option

        get_mm_scomp_solver(
         num_materials_vel,
         project_option_op,
         state_index,
         scomp,
         ncomp,
         ncomp_check);

        nsolve=AMREX_SPACEDIM;
        nsolveMM=nsolve*num_materials_vel;

        if (state_index!=State_Type)
         amrex::Error("state_index invalid");
        if (ncomp_check!=nsolveMM)
         amrex::Error("ncomp_check invalid");

        allocate_array(1,nsolveMM,-1,HOOP_FORCE_MARK_MF);
         // update_state==1:
         //  unp1(1)=unp1(1)/(one+param2*hoop_force_coef)
         // update_state==0:
         //  unp1(1)=unp1(1)-param2*hoop_force_coef*un(1)
        update_state=0;  
        diffuse_hoopALL(REGISTER_MARK_MF,BOUSSINESQ_TEMP_MF,
         HOOP_FORCE_MARK_MF,update_state);

        // HOfab=-div(2 mu D)-HOOP_FORCE_MARK_MF
        // calls: UPDATESEMFORCE in GODUNOV_3D.F90
        if ((viscous_enable_spectral==1)||   // SEM space and time
            (viscous_enable_spectral==3)) {  // SEM time
         update_SEM_forcesALL(project_option_op,REGISTER_MARK_MF,
          update_spectralF,update_stableF);
        } else if (viscous_enable_spectral==0) {
         // do nothing
        } else
         amrex::Error("viscous_enable_spectral invalid");

        override_enable_spectral(save_enable_spectral);

        delete_array(PRESPC2_MF); // pressure
        delete_array(REGISTER_MARK_MF); // velocity
        delete_array(BOUSSINESQ_TEMP_MF); // temperature
        delete_array(HOOP_FORCE_MARK_MF);
        delete_array(NEG_MOM_FORCE_MF);
        delete_array(THERMAL_FORCE_MF);

       } else if ((divu_outer_sweeps>=0)&&
                  (divu_outer_sweeps+1<local_num_divu_outer_sweeps)) {
         // do nothing
       } else
        amrex::Error("divu_outer_sweeps invalid do_the_advance");
    
      } else if (slab_step==ns_time_order) {

       // delta=integral_tn^tnp1  f^spectral dt - deltatn F^stable
       for (int ilev=finest_level;ilev>=level;ilev--) {
        NavierStokes& ns_level=getLevel(ilev);
        ns_level.init_splitting_force_SDC();
       }

      } else
       amrex::Error("slab_step invalid");

    } else if ((ns_time_order==1)||(advance_status==0)) {
      // do nothing
    } else
      amrex::Error("ns_time_order or advance_status invalid: do_the..");

   } // divu_outer_sweeps loop

   for (int ilev=finest_level;ilev>=level;ilev--) {
    NavierStokes& ns_level=getLevel(ilev);
    int alloc_flag=0; // deallocate
    ns_level.alloc_gradp_over_rho(alloc_flag);
   }

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
  // MultiFab& DIV_new=get_new_data(DIV_Type,slab_step+1);
  // if compressible: DIV_new=-dt(pnew-padv)/(rho c^2 dt^2)+MDOT_MF dt/vol
  // if incompressible: DIV_new=MDOT_MF dt/vol
  ADVECT_DIV_ALL();

  if (visual_divergence_plot_int>0) {

   SDC_outer_sweeps=ns_time_order-1;
   divu_outer_sweeps=num_divu_outer_sweeps-1;

   int ratio=(nsteps+1)/visual_divergence_plot_int;
   ratio=ratio*visual_divergence_plot_int;
   if (ratio==nsteps+1) {

    // declared in: MacProj.cpp
    getStateDIV_ALL(MACDIV_MF,1);
    if (localMF[MACDIV_MF]->nComp()!=num_materials_vel)
     amrex::Error("localMF[MACDIV_MF]->nComp() invalid");

    int caller_id=1;
     //MACDIV<stuff>.plt (visit can open binary tecplot files)
    writeSanityCheckData(
      "MACDIV",
      "MACDIV_MF: actual div u",
      caller_id,
      localMF[MACDIV_MF]->nComp(), 
      MACDIV_MF,
      -1,  // State_Type==-1 
      -1); // data_dir==-1 (cell centered)

    caller_id=2;
    MultiFab& DIV_new=get_new_data(DIV_Type,slab_step+1);
    //DIV_Type<stuff>.plt (visit can open binary tecplot files)
    writeSanityCheckData(
      "DIV_Type",
      "DIV_Type: -dt(pnew-padv)/(rho c^2 dt^2)+MDOT_MF dt/vol",
      caller_id,
      DIV_new.nComp(), 
      -1,
      DIV_Type,  // State_Type==-1 
      -1); // data_dir==-1 (cell centered)

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
  if (1==1) {
   std::cout << "PROC= " << ParallelDescriptor::MyProc() << 
	   " advance time " << start_advance << '\n';
   std::cout << "PROC= " << ParallelDescriptor::MyProc() <<
	   " total advance time " << total_advance_time << '\n';
  }

  ParallelDescriptor::Barrier();
 } else 
  amrex::Error("advance status invalid");

}  // end subroutine do_the_advance


void push_stack(Vector<int>& stackdata,
  Long& stackptr,int& data) {

 if (stackptr==stackdata.size()-1) {
  std::cout << "stackptr, stackdata_size " << stackptr << ' ' <<
   stackdata.size() << '\n';
  amrex::Error("stack overflow");
 }
 stackptr++;
 stackdata[stackptr]=data;

}

void pop_stack(Vector<int>& stackdata,
  Long& stackptr,int& data) {

 if (stackptr<0)
  amrex::Error("stack is empty");
 data=stackdata[stackptr];
 stackptr--;

}

void cross_check(Vector<int>& levelcolormap,
  Vector<int>& stackdata,
  std::vector<bool>& grid_color,int i) {

 if (grid_color.size()>stackdata.size())
  amrex::Error("stackdata too small");
 Long stackptr=-1;
 
 for (int j=0;j<levelcolormap.size();j++) {
  Long k=i*levelcolormap.size()+j;
  if (grid_color[k]==true) 
   push_stack(stackdata,stackptr,j);
 }

 while (stackptr>=0) {
  int j;
  pop_stack(stackdata,stackptr,j); 

  if (levelcolormap[j]==0) {
   levelcolormap[j]=levelcolormap[i];
   if (levelcolormap[j]<1)
    amrex::Error("levelcolormap invalid");

   for (int jj=0;jj<levelcolormap.size();jj++) {
    Long k=j*levelcolormap.size()+jj;
    if (grid_color[k]==true) 
     push_stack(stackdata,stackptr,jj);
   }

  } else if (levelcolormap[j]!=levelcolormap[i])
   amrex::Error("something wrong in cross_check");
 }  // stackptr>=0

} // subroutine cross_check



void NavierStokes::correct_colors(
 int idx_color,int base_level,
 Vector<int> domaincolormap,int total_colors,
 int max_colors_level) {

 int finest_level=parent->finestLevel();
 if (base_level>=finest_level)
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
  FORT_LEVELRECOLOR(
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
 ns_reconcile_d_num(175);

}  // subroutine correct_colors

// type_flag should have size nmat.
// type_flag[i]=1 if fluid "i" exists. (note, fictitious solid on the
//  boundaries will show up as existing if ngrow>0)
//
void NavierStokes::assign_colors(
 int& fully_covered,
 int idx_color,int idx_type,
 Vector<int>& colormax,Vector<int> type_flag) {

 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid assign_colors");

 fully_covered=0;

 int nmat=num_materials;

 int ipass_max=2;  
 int typedim=type_flag.size();
 if (typedim!=nmat)
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
 for (int ipass=0;ipass<ipass_max;ipass++) {

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
   FORT_COLORFILL(
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
  ns_reconcile_d_num(176);

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


 } // ipass

 if (ipass_max==2) {
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
    maskmf,check_corners);
  } else
   amrex::Error("max_colors_grid_array[0] invalid");
 } else
  amrex::Error("ipass_max invalid");
  
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
  amrex::Error("S_crse invalid");
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
 type_coarse_fine.copy(*localMF[idx_type],0,0,1);

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

  FORT_AVGDOWNCOLOR(
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
 ns_reconcile_d_num(177);

 S_crse->copy(crse_S_fine,0,0,1);
}


// maskfinemf corresponds to level+1
// components are: color1,type1,color2,type2,color3,type3
MultiFab* NavierStokes::CopyFineToCoarseColor(int idx_color,int idx_type) {

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
  amrex::Error("S_crse invalid");
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
  
  FORT_COPYFINECOARSECOLOR(
    prob_lo,dxf,&bfact_f,&bfact,
    xlo_fine,dx,
    c_dat,ARLIM(clo),ARLIM(chi),
    f_dat,ARLIM(flo),ARLIM(fhi),
    typef.dataPtr(),ARLIM(typef.loVect()),ARLIM(typef.hiVect()),
    maskfinefab.dataPtr(),
    ARLIM(maskfinefab.loVect()),ARLIM(maskfinefab.hiVect()),
    ovlo,ovhi,lofine,hifine);
 } // mfi
} // omp
 ns_reconcile_d_num(178);

 mf->copy(crse_S_fine,0,0,6);

 delete maskfinemf;

 return mf;
}

// c++ example of iterating components of a FAB
// for more examples, check FArrayBox.cpp:
// for (IntVect p = bx.smallEnd(); p <= bx.bigEnd(); bx.next(p)) ...
void NavierStokes::sync_colors(
 int idx_color,int idx_type,
 Vector<int> color_per_grid,
 Vector<int>& colormax,
 int max_colors_grid,
 MultiFab* maskmf,
 int check_corners) {

 int finest_level=parent->finestLevel();

 if ((check_corners!=0)&&(check_corners!=1))
  amrex::Error("check_corners invalid");

 MultiFab* typemf=localMF[idx_type];
 MultiFab* colormf=localMF[idx_color];
 int number_grids=grids.size();

 int total_colors=0;
 colormf->FillBoundary(geom.periodicity());

 if (verbose>0)
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "in sync colors level= " << level << '\n';
   std::cout << "ngrids,maxcolorgrid " << number_grids << ' ' <<
    max_colors_grid << '\n';
  }

 int Nside=number_grids*max_colors_grid;
 if (Nside<=0) {
  std::cout << "cannot have the coarse level completely covered by a finer\n";
  std::cout << "level for the sync_colors routine\n";
  std::cout << "set the coarse level resolution to be level 1 and reduce\n";
  std::cout << "max_level by 1\n";
  amrex::Error("cannot have Nside 0");
 }

 Long Nside2=Nside*Nside;
 Vector< std::vector<bool> > grid_color_array;
 grid_color_array.resize(thread_class::nthreads);

// COLORING LOOP
#ifdef _OPENMP
#pragma omp parallel
#endif
{
 int tid=ns_thread();

 grid_color_array[tid].resize(Nside2);
 for (Long i=0;i<Nside2;i++)
  grid_color_array[tid][i]=false;
  
   // set the diagonal of grid_color_array
 for (int igrid=0;igrid<number_grids;igrid++) {
  for (int icolor=1;icolor<=color_per_grid[igrid];icolor++) {
   Long i=max_colors_grid*igrid+icolor-1;
   Long k=Nside*i+i;
   grid_color_array[tid][k]=true;
  } // icolor
 } // igrid
} // omp

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(typemf->boxArray().d_numPts());

   // grid_color_array[tid][i][j]=1 if color i neighbors color j
   // and they are both not covered by a finer level.
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
     int icolor=(int) (colorfab(p)+0.5);
     int primary_type=(int) (typefab(p)+0.5);
     if (icolor>0) {
#if (AMREX_SPACEDIM==3)
      for (int ii=-1;ii<=1;ii++)
      for (int jj=-1;jj<=1;jj++)
      for (int kk=-1;kk<=1;kk++) {
       IntVect pofs(i+ii,j+jj,k+kk);
       int idist=std::abs(ii)+std::abs(jj)+std::abs(kk);
#endif
#if (AMREX_SPACEDIM==2)
      for (int ii=-1;ii<=1;ii++)
      for (int jj=-1;jj<=1;jj++) {
       IntVect pofs(i+ii,j+jj);
       int idist=std::abs(ii)+std::abs(jj);
#endif
       if ((check_corners==1)||(idist<=1)) {

        int jcolor=(int) (colorfab(pofs)+0.5);
        Real mask2=maskfab(pofs); 
        int secondary_type=(int) (typefab(pofs)+0.5);
        if (mask2==1.0) {
         if (jcolor<=0) {
          std::cout << "level= " << level << '\n';
          std::cout << "ii,jj= " << ii << ' ' << jj << '\n';
          std::cout << "i,j= " << i << ' ' << j << '\n';
          std::cout << "ngrids,maxcolorgrid " << number_grids << ' ' <<
            max_colors_grid << '\n';
          for (int dir2=0;dir2<AMREX_SPACEDIM;dir2++) {
           std::cout << "dir,lo,hi " << dir2 << ' ' <<
            lo[dir2] << ' ' << hi[dir2] << '\n';
          }
          std::cout << "jcolor = " << jcolor << '\n';
          amrex::Error("jcolor invalid"); 
         } // jcolor<=0
         if (secondary_type==primary_type) {
          if ((icolor>Nside)||(jcolor>Nside))
           amrex::Error("icolor or jcolor invalid"); 
          Long igrid=Nside*(icolor-1)+jcolor-1;
          grid_color_array[tid_current][igrid]=true;
          igrid=Nside*(jcolor-1)+icolor-1;
          grid_color_array[tid_current][igrid]=true;
         } // primary_type==secondary_type
        } else if (mask2!=0.0)
         amrex::Error("mask2 invalid");

       }  // check_corners==1 or idist<=1

      } // ii,jj,kk
     } else
      amrex::Error("icolor invalid");
    } else if (maskfab(p)!=0.0)
     amrex::Error("maskfab invalid");
   } // i,j,k
 } // mfi
} //omp
 ns_reconcile_d_num(179);

 if (verbose>0)
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "after initializing grid_color_array \n";
  }

  // first reduce grid_color_array from the different threads
 for (int igrid=0;igrid<number_grids;igrid++) {
  for (int icolor=1;icolor<=color_per_grid[igrid];icolor++) {
   Long i_index=max_colors_grid*igrid+icolor-1;

   for (int jgrid=0;jgrid<number_grids;jgrid++) {
    for (int jcolor=1;jcolor<=color_per_grid[jgrid];jcolor++) {

     Long j_index=max_colors_grid*jgrid+jcolor-1;

     Long i=Nside*i_index+j_index;  
   
     for (int tid=0;tid<thread_class::nthreads;tid++) {
      int tempbit=grid_color_array[tid][i];
      if ((tempbit!=0)&&(tempbit!=1))
       amrex::Error("bits can only be 0 or 1");
      if (tempbit==1)
       grid_color_array[0][i]=true;
     }  // tid
    }  // jcolor
   } // jgrid
  } // icolor
 } // igrid

 ParallelDescriptor::Barrier();

 for (int igrid=0;igrid<number_grids;igrid++) {
  for (int icolor=1;icolor<=color_per_grid[igrid];icolor++) {
   Long i_index=max_colors_grid*igrid+icolor-1;

   for (int jgrid=0;jgrid<number_grids;jgrid++) {
    for (int jcolor=1;jcolor<=color_per_grid[jgrid];jcolor++) {

     Long j_index=max_colors_grid*jgrid+jcolor-1;

     Long i=Nside*i_index+j_index;  
   
     int tempbit=grid_color_array[0][i];
     if ((tempbit!=0)&&(tempbit!=1))
      amrex::Error("bits can only be 0 or 1");
     ParallelDescriptor::ReduceIntMax(tempbit);
     if (tempbit==0)
      grid_color_array[0][i]=false;
     else if (tempbit==1)
      grid_color_array[0][i]=true;
     else
      amrex::Error("tempbit invalid");
    }  // jcolor
   } // jgrid
  } // icolor
 } // igrid
 ParallelDescriptor::Barrier();

 if (verbose>0)
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "after ReduceIntMax \n";
  }

   // levelcolormap[i] associates color "i" in the previous scheme with a
   // new absolute coloring scheme on the level.
 Vector<int> levelcolormap;
 levelcolormap.resize(Nside);
 for (int i=0;i<Nside;i++) 
   levelcolormap[i]=0;

 if (verbose>0)
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "before cross_check \n";
  }

 if (ParallelDescriptor::IOProcessor()) {
  Vector<int> stackdata;
  stackdata.resize(Nside2);

  for (int igrid=0;igrid<number_grids;igrid++) {
   for (int icolor=1;icolor<=color_per_grid[igrid];icolor++) {
    int i_index=max_colors_grid*igrid+icolor-1;
    if (levelcolormap[i_index]==0) {
     Long k=i_index*Nside+i_index;
     if (grid_color_array[0][k]==true) {
      total_colors++;
      levelcolormap[i_index]=total_colors;
      if (levelcolormap[i_index]<1)
       amrex::Error("levelcolormap invalid");
      cross_check(levelcolormap,stackdata,grid_color_array[0],i_index);
     }
    } 
   } // icolor
  } // igrid

  stackdata.resize(1);
 }  // IOProcessor

 if (verbose>0)
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "after cross_check \n";
  }

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 int tid=ns_thread();
 grid_color_array[tid].resize(1);
} // omp

 ParallelDescriptor::Barrier();

 for (int igrid=0;igrid<number_grids;igrid++) {
  for (int icolor=1;icolor<=color_per_grid[igrid];icolor++) {
   Long i_index=max_colors_grid*igrid+icolor-1;
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

   FArrayBox& maskfab=(*maskmf)[mfi];
   FArrayBox& colorfab=(*colormf)[mfi];
   int arrsize=levelcolormap.size();

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // in: LEVELSET_3D.F90
   FORT_GRIDRECOLOR(
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
 ns_reconcile_d_num(180);

 colormax[level]=total_colors;

 if (level<finest_level) {

   avgDownColor(idx_color,idx_type);
   
    // 1 ghost cell is initialized
    // components are: color1,type1,color2,type2,color3,type3
    // only copies where maskfine=1
   MultiFab* fine_coarse_color=CopyFineToCoarseColor(idx_color,idx_type);
   fine_coarse_color->FillBoundary(geom.periodicity());

   int max_colors_level=colormax[level+1];
   if (colormax[level]>max_colors_level)
    max_colors_level=colormax[level];
   int arrsize=2*max_colors_level;
   Long arrsize2=arrsize*arrsize;

   Vector< Vector<int> > level_color_array;
   level_color_array.resize(thread_class::nthreads);

// COLORING LOOP
#ifdef _OPENMP
#pragma omp parallel
#endif
{
   int tid=ns_thread();
   level_color_array[tid].resize(arrsize2);
   for (Long i=0;i<arrsize2;i++)
    level_color_array[tid][i]=0;  // false

     // set the diagonal of level_color
   for (int ilevel=0;ilevel<=1;ilevel++) {
    for (int icolor=1;icolor<=colormax[level+ilevel];icolor++) {
     Long i=max_colors_level*ilevel+icolor-1;
     Long k=2*max_colors_level*i+i;
     level_color_array[tid][k]=1;  // true
    }
   } // ilevel
} // omp

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
    const int* tilelo=tilegrid.loVect();
    const int* tilehi=tilegrid.hiVect();
    const int* fablo=fabgrid.loVect();
    const int* fabhi=fabgrid.hiVect();
    int bfact=parent->Space_blockingFactor(level);
    const Real* xlo = grid_loc[gridno].lo();

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    FArrayBox& maskfab=(*maskmf)[mfi];
    FArrayBox& colorfab=(*fine_coarse_color)[mfi];
     // in: LEVELSET_3D.F90
    FORT_LEVELCOLORINIT(
     maskfab.dataPtr(),ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
     colorfab.dataPtr(),ARLIM(colorfab.loVect()),ARLIM(colorfab.hiVect()),
     xlo,dx,
     tilelo,tilehi,
     fablo,fabhi,&bfact,
     level_color_array[tid_current].dataPtr(),
     &max_colors_level,
     &arrsize,&check_corners);
   } //mfi
} //omp
   ns_reconcile_d_num(181);

   delete fine_coarse_color;

   for (int tid=1;tid<thread_class::nthreads;tid++) {
    for (Long i=0;i<arrsize2;i++) {
     if (level_color_array[tid][i]>level_color_array[0][i])
      level_color_array[0][i]=level_color_array[tid][i];
    } // i
   } // tid

   ParallelDescriptor::Barrier();

   for (Long i=0;i<arrsize2;i++)
    ParallelDescriptor::ReduceIntMax(level_color_array[0][i]);

   Vector<int> domaincolormap;
   domaincolormap.resize(2*max_colors_level);
   Vector<int> stackdata;
   stackdata.resize(arrsize2);
   for (int i=0;i<domaincolormap.size();i++)
    domaincolormap[i]=0;
   total_colors=0;

   std::vector<bool> level_color_bool;
   level_color_bool.resize(arrsize2);
   for (Long k=0;k<level_color_bool.size();k++) {
    if (level_color_array[0][k]==0)
     level_color_bool[k]=false;
    else if (level_color_array[0][k]==1)
     level_color_bool[k]=true;
    else
     amrex::Error("level_color_array invalid value");
   }

   for (int i=0;i<domaincolormap.size();i++) {
    if (domaincolormap[i]==0) {
     Long k=i*domaincolormap.size()+i;
     if (level_color_bool[k]==true) {
      total_colors++;
      domaincolormap[i]=total_colors;
      cross_check(domaincolormap,stackdata,level_color_bool,i);
     }
    }
   }

   for (int ilev=finest_level;ilev>=level;ilev--) {
    colormax[ilev]=total_colors;
    NavierStokes& ns_level=getLevel(ilev);
    
    ns_level.correct_colors(idx_color,level,domaincolormap,
     total_colors,max_colors_level);
     // if coarse_type=fine_type, then coarse_color=fine_color otherwise
     // coarse_color=0.
    if (ilev<finest_level)
     ns_level.avgDownColor(idx_color,idx_type);
   }
 } // level<finest_level


} // subroutine sync_colors
 
// type_flag should have size nmat.
// type_flag[i]=1 if fluid "i" exists. (note, fictitious solid on the
//  boundaries will show up as existing if ngrow>0)
void NavierStokes::color_variable(
 int& coarsest_level,
 int idx_color,int idx_type,
 int* color_count,Vector<int> type_flag) {

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
   idx_color,idx_type,colormax,type_flag); 
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

  // FORT_EXTRAPFILL, pc_interp
 Vector<int> scompBC_map;
 scompBC_map.resize(1);
 scompBC_map[0]=0;

  // ngrow=1 scomp=0 ncomp=1 
  // FORT_EXTRAPFILL, pc_interp
 PCINTERP_fill_bordersALL(idx_color,1,0,1,State_Type,scompBC_map);

 if (verbose>0)
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "coarsest_level= " << coarsest_level << '\n';
   for (ilev=coarsest_level;ilev<=finest_level;ilev++)
    std::cout << "ilev,colormax_ilev " << ilev << ' ' << 
     colormax[ilev] << '\n';
  }
} // subroutine color_variable


void
NavierStokes::ColorSum(
 int operation_flag, //=0 or 1
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

 int finest_level=parent->finestLevel();
 bool use_tiling=ns_tiling;

 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

 if (level>finest_level)
  amrex::Error("level invalid ColorSum");

 if (ngrow_distance!=4)
  amrex::Error("ngrow_distance invalid");

 int rzflag=0;
 if (geom.IsRZ())
  rzflag=1;
 else if (geom.IsCartesian())
  rzflag=0;
 else if (geom.IsCYLINDRICAL())
  rzflag=3;
 else
  amrex::Error("CoordSys bust 30");

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
  if (ncomp_mdot_test==2*nten) {
   // do nothing
  } else
   amrex::Error("ncomp_mdot_test invalid");
  if ((ncomp_mdot==ncomp_mdot_test)&&
      (ncomp_mdot==ncomp_mdot_alloc)) {
   // do nothing
  } else
   amrex::Error("ncomp_mdot_test or ncomp_mdot_alloc invalid");
  ncomp_mdot_test=mdot_complement->nComp();
  if (ncomp_mdot_test==2*nten) {
   // do nothing
  } else
   amrex::Error("ncomp_mdot_test invalid");
 } else
  amrex::Error("ncomp_mdot invalid"); 

 if (operation_flag==0) {

  for (int i=0;i<num_colors;i++) {
   clear_blobdata(i,level_blobdata);
   for (int j=0;j<ncomp_mdot_alloc;j++) {
    level_mdot_data[i][j]=0.0;
    level_mdot_data_redistribute[i][j]=0.0;
    level_mdot_comp_data[i][j]=0.0;
    level_mdot_comp_data_redistribute[i][j]=0.0;
   }
  } // i=0..num_colors-1

 } else if (operation_flag==1) {

  if (sweep_num==0) {
   if (ncomp_mdot>=1) {
    for (int i=0;i<num_colors;i++) {
     int j=0;
     for (j=0;j<ncomp_mdot_alloc;j++) {
      level_mdot_data_redistribute[i][j]=0.0;
      level_mdot_comp_data_redistribute[i][j]=0.0;
     }
     if (j==2*nten) {
      // do nothing
     } else 
      amrex::Error("expecting j==2*nten");
    } // i=0..num_colors-1
   } else
    amrex::Error("ncomp_mdot invalid");
  } else
   amrex::Error("sweep_num invalid");

 } else
  amrex::Error("operation_flag invalid");

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

  if (operation_flag==0) {
   // do nothing
  } else if (operation_flag==1) {
   
   counter=0;
   mdot_counter=0;
 
   for (int i=0;i<num_colors;i++) {
    copy_from_blobdata(i,counter,level_blob_array[tid],level_blobdata);
    if ((level_blobdata[i].im>=1)&&
        (level_blobdata[i].im<=nmat)) {
     level_blob_type_array[tid][i]=level_blobdata[i].im;
    } else
     amrex::Error("level_blobdata[i].im invalid");

    int j=0;
    for (j=0;j<ncomp_mdot_alloc;j++) {
     level_mdot_array[tid][mdot_counter]=level_mdot_data[i][j];
     level_mdot_comp_array[tid][mdot_counter]=level_mdot_comp_data[i][j];
     mdot_counter++;
    }
    if (j==2*nten) {
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
  debug_ngrow(AREA_MF+dir,1,355);
  MultiFab& Umac_new=get_new_data(Umac_Type+dir,slab_step+1);
  if (localMF[AREA_MF+dir]->boxArray()!=Umac_new.boxArray())
   amrex::Error("area_mf boxarrays do not match");
 } // dir=0..sdim-1
 
 debug_ngrow(VOLUME_MF,0,722);
 VOF_Recon_resize(1,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,1,31);
 if (localMF[SLOPE_RECON_MF]->nComp()!=nmat*ngeom_recon)
  amrex::Error("localMF[SLOPE_RECON_MF]->nComp() invalid");

 int do_face_decomp=0;

 if (sweep_num==0) {
  getStateDist_localMF(LS_COLORSUM_MF,1,cur_time_slab,20);
  getStateDen_localMF(DEN_COLORSUM_MF,1,cur_time_slab);
   // velocity + pressure
  getState_localMF(VEL_COLORSUM_MF,1,0,AMREX_SPACEDIM+1,cur_time_slab);

  makeFaceFrac(tessellate,ngrow_distance,FACEFRAC_MM_MF,do_face_decomp);
  ProcessFaceFrac(tessellate,FACEFRAC_MM_MF,FACEFRAC_SOLVE_MM_MF,0);
  makeCellFrac(tessellate,0,CELLFRAC_MM_MF);
 } else if (sweep_num==1) {
  // do nothing
 } else
  amrex::Error("sweep_num invalid");

 if (localMF[LS_COLORSUM_MF]->nGrow()!=1)
  amrex::Error("localMF[LS_COLORSUM_MF]->nGrow()!=1");

 if (localMF[LS_COLORSUM_MF]->nComp()!=(1+AMREX_SPACEDIM)*nmat)
  amrex::Error("localMF[LS_COLORSUM_MF]->nComp()!=(1+AMREX_SPACEDIM)*nmat");
 if (localMF[VEL_COLORSUM_MF]->nComp()!=AMREX_SPACEDIM+1)
  amrex::Error("localMF[VEL_COLORSUM_MF]->nComp()!=AMREX_SPACEDIM+1");
 if (localMF[DEN_COLORSUM_MF]->nComp()!=num_state_material*nmat)
  amrex::Error("localMF[DEN_COLORSUM_MF]->nComp()!=num_state_material*nmat");

   // (nmat,sdim,2,sdim+1) area+centroid on each face of a cell.
 int nface=nmat*AMREX_SPACEDIM*2*(1+AMREX_SPACEDIM);
  // (nmat,nmat,2)  left material, right material, frac_pair+dist_pair
 int nface_dst=nmat*nmat*2;
  // (nmat,nmat,3+sdim)
  // im_inside,im_outside,3+sdim --> area, dist_to_line, dist, line normal.
 int ncellfrac=nmat*nmat*(3+AMREX_SPACEDIM);

 debug_ngrow(FACEFRAC_MM_MF,ngrow_distance,722);
 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(FACEFRAC_SOLVE_MM_MF+dir,0,722);
 debug_ngrow(CELLFRAC_MM_MF,0,722);

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
 int nstate=num_materials_vel*(AMREX_SPACEDIM+1)+
  nmat*(num_state_material+ngeom_raw)+1;
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

  FArrayBox& facefab=(*localMF[FACEFRAC_MM_MF])[mfi];
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
  Vector<int> velbc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // in: LEVELSET_3D.F90
  FORT_GETCOLORSUM(
   &tid_current,
   &operation_flag,
   &sweep_num,
   &tessellate,
   distribute_mdot_evenly.dataPtr(),
   constant_volume_mdot.dataPtr(),
   latent_heat.dataPtr(),
   distribute_from_target.dataPtr(),
   constant_density_all_time.dataPtr(),
   &dt_slab,
   dx,xlo,
   &nmat,
   &nten,
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
   facefab.dataPtr(),ARLIM(facefab.loVect()),ARLIM(facefab.hiVect()),
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
   &rzflag,
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
   &num_elements_blobclass,
   &ncomp_mdot_alloc,
   &ncomp_mdot,
   levelbc.dataPtr(),
   velbc.dataPtr(),
   material_type_lowmach.dataPtr(),
   &nface,
   &nface_dst,
   &ncellfrac);
 } // mfi
} // omp
 ns_reconcile_d_num(182);

 if (operation_flag==0) {

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

 } else if (operation_flag==1) {

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

 if (operation_flag==0) {

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

 } else if (operation_flag==1) {

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

 if ((sweep_num==0)&&(operation_flag==0)) {
  // do nothing
 } else if ((sweep_num==1)||(operation_flag==1)) {
  delete_localMF(LS_COLORSUM_MF,1);
  delete_localMF(DEN_COLORSUM_MF,1);
  delete_localMF(VEL_COLORSUM_MF,1);
 } else
  amrex::Error("sweep_num or operation_flag invalid");

}  // subroutine ColorSum



void
NavierStokes::LowMachDIVU(
 int sweep_num,
 MultiFab* typemf,
 MultiFab* color,
 MultiFab* mdot_local, 
 MultiFab* mdot_global, 
 Vector<blobclass> cum_blobdata,
 Vector< Vector<Real> >& level_mdot_data,
 Vector< Vector<Real> > cum_mdot_data
 ) {

 int finest_level=parent->finestLevel();
 bool use_tiling=ns_tiling;

 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

 if (level>finest_level)
  amrex::Error("level invalid ColorSum");

 if (ngrow_distance!=4)
  amrex::Error("ngrow_distance invalid");

 int rzflag=0;
 if (geom.IsRZ())
  rzflag=1;
 else if (geom.IsCartesian())
  rzflag=0;
 else if (geom.IsCYLINDRICAL())
  rzflag=3;
 else
  amrex::Error("CoordSys bust 30");

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

 debug_ngrow(VOLUME_MF,0,722);
 VOF_Recon_resize(1,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,1,31);
 if (localMF[SLOPE_RECON_MF]->nComp()!=nmat*ngeom_recon)
  amrex::Error("localMF[SLOPE_RECON_MF]->nComp() invalid");

 if (sweep_num==0) {
  getStateDist_localMF(LS_COLORSUM_MF,1,cur_time_slab,20);
  getStateDen_localMF(DEN_COLORSUM_MF,1,cur_time_slab);
 } else if (sweep_num==1) {
  // do nothing
 } else
  amrex::Error("sweep_num invalid");

 if (localMF[LS_COLORSUM_MF]->nGrow()!=1)
  amrex::Error("localMF[LS_COLORSUM_MF]->nGrow()!=1");

 if (localMF[LS_COLORSUM_MF]->nComp()!=(1+AMREX_SPACEDIM)*nmat)
  amrex::Error("localMF[LS_COLORSUM_MF]->nComp()!=(1+AMREX_SPACEDIM)*nmat");
 if (localMF[DEN_COLORSUM_MF]->nComp()!=num_state_material*nmat)
  amrex::Error("localMF[DEN_COLORSUM_MF]->nComp()!=num_state_material*nmat");

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
  FORT_GET_LOWMACH_DIVU(
   &tid_current,
   &sweep_num,
   constant_density_all_time.dataPtr(),
   &dt_slab,
   dx,xlo,
   &nmat,
   &nten,
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
   &rzflag,
   &num_colors,
   cum_blob_array.dataPtr(),
   level_mdot_array[tid_current].dataPtr(),
   &blob_array_size,
   &mdot_array_size,
   &num_elements_blobclass,
   material_type_lowmach.dataPtr());
 } // mfi
} // omp
 ns_reconcile_d_num(182);

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
 Vector<int>& type_flag, 
 Vector<blobclass> blobdata) {

 int finest_level=parent->finestLevel();

 if ((coarsest_level<0)||(coarsest_level>finest_level))
  amrex::Error("coarsest_level invalid");

 if (level!=0)
  amrex::Error("level=0 in LowMachDIVUALL");

 Real problo_array[AMREX_SPACEDIM];
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
  problo_array[dir]=geom.ProbLo(dir);
 int symmetric_flag=1;
 for (int dir=0;dir<AMREX_SPACEDIM-1;dir++) {
  if (phys_bc.lo(dir)!=Symmetry)
   symmetric_flag=0;
  if (problo_array[dir]!=0.0)
   symmetric_flag=0;
 } // dir=0..sdim-2
 if (geom.IsRZ()) {
  if (symmetric_flag==1) {
   // do nothing
  } else if (symmetric_flag==0) {
   amrex::Error("we force symmetric_flag==1 if geom.IsRZ()");
  } else
   amrex::Error("symmetric_flag error");
 }  // IsRZ?

 if (color_count==0)
  amrex::Error("num_colors=0 in ColorSumALL");

   // initializes MDOT_LOCAL_MF to 0.0
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
  int caller_id=500;
  writeSanityCheckData(
   "MDOT_LOCAL",
   "in: NavierStokes::LowMachDIVUALL", 
   caller_id,
   localMF[MDOT_LOCAL_MF]->nComp(),
   MDOT_LOCAL_MF, 
   -1, //State_Type==-1
   -1); // data_dir==-1 (cell centered)
 }

 delete_array(MDOT_LOCAL_MF);

} // end subroutine LowMachDIVUALL


void NavierStokes::copy_to_blobdata(int i,int& counter,
  Vector<Real>& blob_array,Vector<blobclass>& blobdata) {

 int nmat=num_materials;
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
  blobdata[i].blob_matrix[dir]=blob_array[counter];
  counter++;
 }
 for (int dir=0;dir<3*(2*AMREX_SPACEDIM);dir++) {
  blobdata[i].blob_RHS[dir]=blob_array[counter];
  counter++;
 }
 for (int dir=0;dir<3*(2*AMREX_SPACEDIM);dir++) {
  blobdata[i].blob_velocity[dir]=blob_array[counter];
  counter++;
 }
 for (int dir=0;dir<2*(2*AMREX_SPACEDIM);dir++) {
  blobdata[i].blob_integral_momentum[dir]=blob_array[counter];
  counter++;
 }
 blobdata[i].blob_energy=blob_array[counter];
 counter++;
 for (int dir=0;dir<3;dir++) {
  blobdata[i].blob_mass_for_velocity[dir]=blob_array[counter];
  counter++;
 }
 blobdata[i].blob_volume=blob_array[counter];
 counter++;
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  blobdata[i].blob_center_integral[dir]=blob_array[counter];
  counter++;
 }
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  blobdata[i].blob_center_actual[dir]=blob_array[counter];
  counter++;
 }
 blobdata[i].blob_perim=blob_array[counter];
 counter++;
 for (int imnbr=0;imnbr<nmat;imnbr++) {
  blobdata[i].blob_perim_mat[imnbr]=blob_array[counter];
  counter++;
 }
 for (int im1=0;im1<nmat;im1++) {
  for (int im2=0;im2<nmat;im2++) {
   blobdata[i].blob_triple_perim[im1][im2]=blob_array[counter];
   counter++;
  } // im2
 } // im1

 blobdata[i].blob_cell_count=blob_array[counter];
 counter++;
 blobdata[i].blob_cellvol_count=blob_array[counter];
 counter++;
 blobdata[i].blob_mass=blob_array[counter];
 counter++;
 blobdata[i].blob_pressure=blob_array[counter];
 counter++;

} // end subroutine copy_to_blobdata


void NavierStokes::copy_blobdata(Vector<blobclass>& dest_blobdata,
  Vector<blobclass>& source_blobdata) {

 int nmat=num_materials;
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
   for (int imnbr=0;imnbr<nmat;imnbr++) {
    dest_blobdata[i].blob_perim_mat[imnbr]=
      source_blobdata[i].blob_perim_mat[imnbr];
   }
   for (int im1=0;im1<nmat;im1++) {
    for (int im2=0;im2<nmat;im2++) {
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
  } // i=0..num_colors-1

 } else
  amrex::Error("num_colors invalid");

} // end subroutine copy_blobdata


void NavierStokes::sum_blobdata(int i,
  Vector<blobclass>& blobdata,
  Vector<blobclass>& level_blobdata,int sweep_num) {

 int nmat=num_materials;
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

  for (int imnbr=0;imnbr<nmat;imnbr++)
   blobdata[i].blob_perim_mat[imnbr]+=
    level_blobdata[i].blob_perim_mat[imnbr];

  for (int im1=0;im1<nmat;im1++)
   for (int im2=0;im2<nmat;im2++)
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
  for (int dir=0;dir<3;dir++)
   blobdata[i].blob_mass_for_velocity[dir]+=
    level_blobdata[i].blob_mass_for_velocity[dir];

 } else
  amrex::Error("sweep_num invalid");

} // end subroutine sum_blobdata

void NavierStokes::copy_from_blobdata(int i,int& counter,
  Vector<Real>& blob_array,Vector<blobclass>& blobdata) {

 int nmat=num_materials;
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
  blob_array[counter]=blobdata[i].blob_matrix[dir];
  counter++;
 }
 for (int dir=0;dir<3*(2*AMREX_SPACEDIM);dir++) {
  blob_array[counter]=blobdata[i].blob_RHS[dir];
  counter++;
 }
 for (int dir=0;dir<3*(2*AMREX_SPACEDIM);dir++) {
  blob_array[counter]=blobdata[i].blob_velocity[dir];
  counter++;
 }
 for (int dir=0;dir<2*(2*AMREX_SPACEDIM);dir++) {
  blob_array[counter]=blobdata[i].blob_integral_momentum[dir];
  counter++;
 }
 blob_array[counter]=blobdata[i].blob_energy;
 counter++;
 for (int dir=0;dir<3;dir++) {
  blob_array[counter]=blobdata[i].blob_mass_for_velocity[dir];
  counter++;
 }
 blob_array[counter]=blobdata[i].blob_volume;
 counter++;
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  blob_array[counter]=blobdata[i].blob_center_integral[dir];
  counter++;
 }
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  blob_array[counter]=blobdata[i].blob_center_actual[dir];
  counter++;
 }
 blob_array[counter]=blobdata[i].blob_perim;
 counter++;
 for (int imnbr=0;imnbr<nmat;imnbr++) {
  blob_array[counter]=blobdata[i].blob_perim_mat[imnbr];
  counter++;
 }
 for (int im1=0;im1<nmat;im1++) {
  for (int im2=0;im2<nmat;im2++) {
   blob_array[counter]=blobdata[i].blob_triple_perim[im1][im2];
   counter++;
  } // im2
 } // im1
 blob_array[counter]=blobdata[i].blob_cell_count;
 counter++;
 blob_array[counter]=blobdata[i].blob_cellvol_count;
 counter++;
 blob_array[counter]=blobdata[i].blob_mass;
 counter++;
 blob_array[counter]=blobdata[i].blob_pressure;
 counter++;

}  // end subroutine copy_from_blobdata

void NavierStokes::clear_blobdata(int i,Vector<blobclass>& blobdata) {

 int nmat=num_materials;

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
 blobdata[i].blob_perim_mat.resize(nmat);
 blobdata[i].blob_triple_perim.resize(nmat);
 for (int imnbr=0;imnbr<nmat;imnbr++) {
  blobdata[i].blob_triple_perim[imnbr].resize(nmat);
  blobdata[i].blob_perim_mat[imnbr]=0.0;
 }
 for (int im1=0;im1<nmat;im1++)
  for (int im2=0;im2<nmat;im2++)
   blobdata[i].blob_triple_perim[im1][im2]=0.0;
 blobdata[i].im=0;

} // end subroutine clear_blobdata

void
NavierStokes::ColorSumALL(
 int operation_flag, // =0 or 1
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

 if (operation_flag==0) {

  // type_flag[im]=1 if material im exists in the domain.
  // type_mf(i,j,k)=im if material im dominates cell (i,j,k)
  // updates one ghost cell of TYPE_MF
  // fluid(s) and solid(s) tessellate the domain.
  TypeALL(idx_type,type_flag);

  // color_count=number of colors
  // ngrow=1, FORT_EXTRAPFILL, pc_interp for COLOR_MF
  color_variable(coarsest_level,
   idx_color,idx_type,&color_count,type_flag);

 } else if (operation_flag==1) {

  // do nothing

 } else
  amrex::Error("operation_flag invalid");

 Real problo_array[AMREX_SPACEDIM];
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
  problo_array[dir]=geom.ProbLo(dir);
 int symmetric_flag=1;
 for (int dir=0;dir<AMREX_SPACEDIM-1;dir++) {
  if (phys_bc.lo(dir)!=Symmetry)
   symmetric_flag=0;
  if (problo_array[dir]!=0.0)
   symmetric_flag=0;
 } // dir=0..sdim-2
 if (geom.IsRZ()) {
  if (symmetric_flag==1) {
   // do nothing
  } else if (symmetric_flag==0) {
   amrex::Error("we force symmetric_flag==1 if geom.IsRZ()");
  } else
   amrex::Error("symmetric_flag error");
 }  // IsRZ?

 if (color_count==0)
  amrex::Error("num_colors=0 in ColorSumALL");

 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

 int ncomp_mdot_alloc=1;
 int ncomp_mdot=0;

 if (idx_mdot==-1) {
  // do nothing
 } else if (idx_mdot>=0) {
  ncomp_mdot=localMF[idx_mdot]->nComp();
  ncomp_mdot_alloc=ncomp_mdot;
  if (ncomp_mdot==2*nten) {
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

 if (operation_flag==0) {

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
    if (j==2*nten) {
     // do nothing
    } else
     amrex::Error("expecting j==2*nten");
   } else if (idx_mdot==-1) {
    // check nothing
   } else
    amrex::Error("idx_mdot invalid");

  }  // i=0..color_count-1

 } else if (operation_flag==1) {

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
   if (j==2*nten) {
    // do nothing
   } else
    amrex::Error("expecting j==2*nten");
  } else if (idx_mdot==-1) {
   // check nothing
  } else
   amrex::Error("idx_mdot invalid");

 } // i=0..color_count-1

 int num_sweeps=2;
 if (operation_flag==0) {
  // do nothing
 } else if (operation_flag==1) { 
   // (dest,source)
  copy_blobdata(level_blobdata,blobdata);
  for (int i=0;i<color_count;i++) {
   int j=0;
   for (j=0;j<ncomp_mdot_alloc;j++) {
    level_mdot_data[i][j]=mdot_data[i][j];
    level_mdot_comp_data[i][j]=mdot_comp_data[i][j];
   }
   if (idx_mdot>=0) {
    if (j==2*nten) {
     // do nothing
    } else
     amrex::Error("expecting j==2*nten");
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
   } else if (ncomp_mdot==2*nten) {
    mdot=ns_level.localMF[idx_mdot];
    mdot_complement=ns_level.localMF[idx_mdot_complement];
   } else {
    mdot=nullptr;
    mdot_complement=nullptr;
    amrex::Error("ncomp_mdot invalid");
   }

   ns_level.ColorSum(
    operation_flag, // =0 or 1
    tessellate,  // =1 or 3
    sweep_num,
    ncomp_mdot_alloc,
    ncomp_mdot,
    ns_level.localMF[idx_type],
    ns_level.localMF[idx_color],
    mdot,
    mdot_complement,
    level_blobdata,
    blobdata,
    level_mdot_data,
    level_mdot_comp_data,
    mdot_data,
    mdot_comp_data,
    level_mdot_data_redistribute,
    level_mdot_comp_data_redistribute
    );

   if (operation_flag==0) {

    if (sweep_num==0) {

     for (int i=0;i<color_count;i++) {
      // blob_volume, blob_center_integral, blob_perim, blob_perim_mat,
      // blob_triple_perim, 
      // blob_cell_count,blob_cellvol_count, blob_mass,
      // blob_pressure
      sum_blobdata(i,blobdata,level_blobdata,sweep_num);

      int j=0;
      for (j=0;j<ncomp_mdot_alloc;j++) {
       mdot_data[i][j]+=level_mdot_data[i][j];
       mdot_comp_data[i][j]+=level_mdot_comp_data[i][j];
      }
      if (idx_mdot>=0) {
       if (j==2*nten) {
        // do nothing
       } else
        amrex::Error("expecting j==2*nten");
      } else if (idx_mdot==-1) {
       // check nothing
      } else
       amrex::Error("idx_mdot invalid");

      if ((level_blobdata[i].im>=1)&&
          (level_blobdata[i].im<=nmat)) {
       int im_test=level_blobdata[i].im;
       if ((im_test<1)||(im_test>nmat))
        amrex::Error("im_test invalid");
       blobdata[i].im=im_test;
      } else if (level_blobdata[i].im==0) {
       // do nothing
      } else
       amrex::Error("level_blobdata[i].im invalid");

     }  // i=0..color_count-1

    } else if (sweep_num==1) {

     // blob_energy, blob_matrix, blob_RHS, blob_integral_momentum,
     // blob_mass_for_velocity
     for (int i=0;i<color_count;i++) {
      sum_blobdata(i,blobdata,level_blobdata,sweep_num);
     }  // i=0..color_count-1

    } else
     amrex::Error("sweep_num invalid");

   } else if (operation_flag==1) {

    if (sweep_num==0) {

     for (int i=0;i<color_count;i++) {

      int j=0;
      for (j=0;j<ncomp_mdot_alloc;j++) {
       mdot_data_redistribute[i][j]+=
	   level_mdot_data_redistribute[i][j];
       mdot_comp_data_redistribute[i][j]+=
	   level_mdot_comp_data_redistribute[i][j];
      }
      if (j==2*nten) {
       // do nothing
      } else
       amrex::Error("expecting j==2*nsten");

      if ((level_blobdata[i].im>=1)&&
          (level_blobdata[i].im<=nmat)) {
       if (level_blobdata[i].im==blobdata[i].im) {
        // do nothing
       } else
        amrex::Error("level_blobdata[i].im!=blobdata[i].im");
      } else
       amrex::Error("level_blobdata[i].im invalid");
     }  // i=0..color_count-1

    } else
     amrex::Error("sweep_num invalid for operation_flag==1 case");

   } else
    amrex::Error("operation_flag invalid");

  } // ilev=coarsest_level..finest_level

  if (operation_flag==0) {

   if (sweep_num==0) {

    for (int i=0;i<color_count;i++) {
     Real blobvol=blobdata[i].blob_volume;
     if (blobvol>0.0) {
      for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
       blobdata[i].blob_center_actual[dir]= 
        blobdata[i].blob_center_integral[dir]/blobvol;
      }
      if (symmetric_flag==1) {
       for (int dir=0;dir<AMREX_SPACEDIM-1;dir++)
        blobdata[i].blob_center_actual[dir]=0.0;
      }
     } else if (blobvol==0.0) {
      // do nothing
     } else
      amrex::Error("blobvol invalid");
    } // i=0..color_count-1

   } else if (sweep_num==1) {

    for (int i=0;i<color_count;i++) {

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
        }
       } 
       if (imatrix!=matrix_ncomp*(veltype+1))
        amrex::Error("imatrix invalid");
       for (int irow=0;irow<2*AMREX_SPACEDIM;irow++) {
        BB3D[irow]=blobdata[i].blob_RHS[2*AMREX_SPACEDIM*veltype+irow];
        XX3D[irow]=0.0;
        if (irow>=AMREX_SPACEDIM) {
	 if (symmetric_flag==1) {
	  BB3D[irow]=0.0;
	 } else if (symmetric_flag==0) {
	  // do nothing
	 } else
          amrex::Error("symmetric_flag invalid");

         if (veltype==1)
	  BB3D[irow]=0.0;
        }
        if (irow<AMREX_SPACEDIM+1) {
         BB2D[irow]=blobdata[i].blob_RHS[2*AMREX_SPACEDIM*veltype+irow];
         XX2D[irow]=0.0;
         if (irow>=AMREX_SPACEDIM) {
	  if (symmetric_flag==1) {
           BB2D[irow]=0.0;
          } else if (symmetric_flag==0) {
           // do nothing
          } else
           amrex::Error("symmetric_flag invalid");

          if (veltype==1)
  	   BB2D[irow]=0.0;
         }
        } // irow<sdim+1
       }
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
         if (symmetric_flag==1) {
          for (int dir=0;dir<AMREX_SPACEDIM-1;dir++)
           blobdata[i].blob_velocity[2*AMREX_SPACEDIM*veltype+dir]=0.0;
          for (int dir=AMREX_SPACEDIM;dir<2*AMREX_SPACEDIM;dir++)
           blobdata[i].blob_velocity[2*AMREX_SPACEDIM*veltype+dir]=0.0;
         } else if (symmetric_flag==0) {
          // do nothing
         } else {
          amrex::Error("symmetric_flag invalid");
         }
        } else if (AMREX_SPACEDIM==2) {
         for (int dir=0;dir<AMREX_SPACEDIM+1;dir++) {
          blobdata[i].blob_velocity[2*AMREX_SPACEDIM*veltype+dir]=XX2D[dir];
         }
         if ((geom.IsRZ())||(symmetric_flag==1)) {
          int dir=0;
          blobdata[i].blob_velocity[2*AMREX_SPACEDIM*veltype+dir]=0.0;
          dir=AMREX_SPACEDIM;
          blobdata[i].blob_velocity[2*AMREX_SPACEDIM*veltype+dir]=0.0;
         } else if ((!geom.IsRZ())&&(symmetric_flag==0)) {
          // do nothing
         } else {
          amrex::Error("IsRZ or symmetric_flag invalid");
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
  	  blobdata[i].blob_velocity[ibase]*=sqrt(original_KE/proposed_KE);
	 }
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
	   }
	   std::cout << " ------------------------------------\n";
	  }
	 }
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

  } else if (operation_flag==1) {
   // do nothing
  } else
   amrex::Error("operation_flag invalid");

 } // sweep_num=0..1

 if (verbose>=2) {
  if (ParallelDescriptor::IOProcessor()) {

   std::cout << "in color sum color_count = " << color_count << '\n';
   std::cout << "symmetric_flag= " << symmetric_flag << '\n';

   for (int i=0;i<color_count;i++) {

    Real blobvol=blobdata[i].blob_volume;

    if (blobvol==0.0) {
     std::cout << "NULL COLOR i,im " << i << ' ' << 
       blobdata[i].im << '\n';
    } else if (blobvol>0.0) {
     int imbase=blobdata[i].im;
     if ((imbase<1)||(imbase>nmat))
      amrex::Error("imbase invalid");

     std::cout << "i, vol,perim,im " << i << ' ' <<
      blobdata[i].blob_volume << ' ' <<
      blobdata[i].blob_perim << ' ' <<
      imbase << '\n';

     for (int imnbr=0;imnbr<nmat;imnbr++) {
      std::cout << "perim(nbr): i,imnbr,im,perim " << i << ' ' <<
      imnbr+1 << ' ' <<
      imbase << ' ' <<
      blobdata[i].blob_perim_mat[imnbr] << '\n';
     }
     for (int im1=0;im1<nmat;im1++) {
      for (int im2=0;im2<nmat;im2++) {
       std::cout << "blob_triple_perim: i,im1,im2,im,perim " << i << ' ' <<
       im1+1 << ' ' << im2+1 << ' ' <<
       imbase << ' ' <<
       blobdata[i].blob_triple_perim[im1][im2] << '\n';
      } // im2=0..nmat-1
     } // im1=0..nmat-1

     if (blobvol>0.0) {
      std::cout << "i,x,y,z " << 
       blobdata[i].blob_center_actual[0] << ' ' <<
       blobdata[i].blob_center_actual[1] << ' ' <<
       blobdata[i].blob_center_actual[AMREX_SPACEDIM-1] << '\n';
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

} // subroutine ColorSumALL


void
NavierStokes::Type_level(
  MultiFab* typemf,Vector<int>& type_flag) {

 int finest_level=parent->finestLevel();
 int nmat=num_materials;

 int typedim=type_flag.size();
 if (typedim!=nmat)
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

 MultiFab* LS=getStateDist(1,cur_time_slab,21);
 if (LS->nComp()!=nmat*(1+AMREX_SPACEDIM))
  amrex::Error("LS->nComp() invalid");
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

  FArrayBox& LSfab=(*LS)[mfi];
  FArrayBox& typefab=(*typemf)[mfi];

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // updates one ghost cell.
   // in: LEVELSET_3D.F90
   //  for each cell,
   //   if is_rigid(nmat,im)==1 and LS>=0 then type=im
  FORT_GETTYPEFAB(
   LSfab.dataPtr(),ARLIM(LSfab.loVect()),ARLIM(LSfab.hiVect()),
   typefab.dataPtr(),ARLIM(typefab.loVect()),ARLIM(typefab.hiVect()),
   xlo,dx,
   tilelo,tilehi,
   fablo,fabhi,&bfact,
   local_type_flag[tid_current].dataPtr(),
   &nmat);
 } // mfi
} // omp
 ns_reconcile_d_num(183);

 for (int tid=1;tid<thread_class::nthreads;tid++) {
  for (int im=0;im<nmat;im++) { 
   if (local_type_flag[tid][im]>local_type_flag[0][im])
    local_type_flag[0][im]=local_type_flag[tid][im];
  } // im
 } // tid

 ParallelDescriptor::Barrier();

 for (int im=0;im<nmat;im++) {
  ParallelDescriptor::ReduceIntMax(local_type_flag[0][im]);
  type_flag[im]=local_type_flag[0][im];
 }

 delete LS;
}  // subroutine Type_level


void NavierStokes::TypeALL(int idx_type,Vector<int>& type_flag) {

 int finest_level=parent->finestLevel();
 int nmat=num_materials;

 type_flag.resize(nmat);
 for (int im=0;im<nmat;im++) {
  type_flag[im]=0;
 }
 allocate_array(1,1,-1,idx_type);
 if (level!=0)
  amrex::Error("level=0 in TypeALL");

  // updates one ghost cell.
 for (int k = 0; k <= finest_level; k++) {
  NavierStokes& ns_level = getLevel(k);
  ns_level.Type_level(ns_level.localMF[idx_type],type_flag);
 }
 int color_counter=0;
 for (int im=0;im<nmat;im++) {
  color_counter+=type_flag[im];
 }
 if (verbose>0)
  if (ParallelDescriptor::IOProcessor()) {
   for (int im=0;im<nmat;im++) 
    std::cout << "TypeALL im,type " << im << ' ' <<
     type_flag[im] << ' ' << '\n';
   std::cout << "TypeALL color_counter= " << color_counter << '\n'; 
  }
} // subroutine TypeALL

void NavierStokes::remove_pressure_work_vars() {

 delete_localMF(UMACSTAR_MF,AMREX_SPACEDIM);
 delete_localMF(GRADPEDGE_MF,AMREX_SPACEDIM);
 delete_localMF(PEDGE_MF,AMREX_SPACEDIM);
 delete_localMF(AMRSYNC_PEDGE_MF,AMREX_SPACEDIM);
 delete_localMF(AMRSYNC_PRES_MF,AMREX_SPACEDIM);

} // end subroutine remove_pressure_work_vars

// called from: NavierStokes::multiphase_project
void NavierStokes::remove_project_variables() {

 delete_localMF(POLDHOLD_MF,1);
 delete_localMF(ONES_MF,1);
 delete_localMF(DOTMASK_MF,1);
 delete_localMF(OUTER_ITER_PRESSURE_MF,1);
}

void NavierStokes::allocate_MAC_velocityALL(int nsolve,int idx) {

 if (level!=0)
  amrex::Error("level invalid allocate_MAC_velocity_ALL");

 if ((nsolve!=1)&&(nsolve!=AMREX_SPACEDIM))
  amrex::Error("nsolve invalid");

 int nmat=num_materials;
 int nsolveMM=nsolve*num_materials_vel;
 int nsolveMM_FACE=nsolveMM;
 if (num_materials_vel==1) {
  // do nothing
 } else if (num_materials_vel==nmat) {
  nsolveMM_FACE*=2;
 } else
  amrex::Error("num_materials_vel invalid");

 int finest_level=parent->finestLevel();
 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

   ns_level.new_localMF(idx+dir,nsolveMM_FACE,0,dir);
   ns_level.setVal_localMF(idx+dir,0.0,0,nsolveMM_FACE,0);

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
} // remove_MAC_velocityALL

// called from:
// NavierStokes::update_SEM_forcesALL
// NavierStokes::multiphase_project
// NavierStokes::diffusion_heatingALL 
void NavierStokes::allocate_FACE_WEIGHT(
 int nsolve,
 int project_option) {
 
 int finest_level=parent->finestLevel();

 bool use_tiling=ns_tiling;

 if (dt_slab<=0.0)
  amrex::Error("cannot have dt_slab<=0 in allocate_FACE_WEIGHT");

 int nmat=num_materials;

 int num_materials_face=num_materials_vel;

 if (project_option_momeqn(project_option)==1) {
  if (num_materials_face!=1)
   amrex::Error("num_materials_face invalid");
 } else if (project_option_momeqn(project_option)==0) {
  num_materials_face=num_materials_scalar_solve;
 } else
  amrex::Error("project_option_momeqn(project_option) invalid30");

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 if ((num_materials_face!=1)&&
     (num_materials_face!=nmat))
  amrex::Error("num_materials_face invalid");

 int nsolveMM=nsolve*num_materials_face;
 int nsolveMM_FACE=nsolveMM;
 if (num_materials_face==1) {
  // do nothing
 } else if (num_materials_face==nmat) {
  nsolveMM_FACE*=2;
 } else
  amrex::Error("num_materials_face invalid");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");
 Vector<int> scomp;
 Vector<int> ncomp;
 int state_index;
 int ncomp_check;

 get_mm_scomp_solver(
  num_materials_face,
  project_option,
  state_index,
  scomp,
  ncomp,
  ncomp_check);

 if (ncomp_check!=nsolveMM)
  amrex::Error("ncomp_check alid");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  debug_ngrow(FACE_VAR_MF+dir,0,850);
 }

 resize_mask_nbr(1);
 debug_ngrow(MASK_NBR_MF,1,851);

 VOF_Recon_resize(1,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,1,31);
 if (localMF[SLOPE_RECON_MF]->nComp()!=nmat*ngeom_recon)
  amrex::Error("localMF[SLOPE_RECON_MF]->nComp() invalid");

 debug_ngrow(CELL_VISC_MF,1,47);
 debug_ngrow(CELL_DEN_MF,1,47);
 if (localMF[CELL_VISC_MF]->nComp()!=nmat+1)
  amrex::Error("localMF[CELL_VISC_MF]->nComp() invalid");
 if (localMF[CELL_DEN_MF]->nComp()!=nmat+1)
  amrex::Error("localMF[CELL_DEN_MF]->nComp() invalid");

 solvability_level_flag=1;

 int bcsize=AMREX_SPACEDIM*2*nsolveMM*grids.size();
 bcpres_array.resize(bcsize);
 for (int gridno=0;gridno<grids.size();gridno++) {

   // presbc declared as presbc(sdim,2) in fortran
   // components ordered at  1,1  2,1  3,1  1,2  2,2  3,2
  Vector<int> presbc;
  getBCArray_list(presbc,state_index,gridno,scomp,ncomp);
  if (presbc.size()!=AMREX_SPACEDIM*2*nsolveMM)
   amrex::Error("presbc.size() invalid");

  int ibase=AMREX_SPACEDIM*2*nsolveMM*gridno;
  int j=0;
  for (int nn=0;nn<nsolveMM;nn++) {
   for (int side=0;side<=1;side++) {
    for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
     int pbc=presbc[j];
     bcpres_array[ibase+j]=pbc;

     if (project_option_singular_possible(project_option)==1) {

      if (project_option==12) { // pressure extrapolation
       singular_possible=1; // non-solid regions can have coefficients all zero.
       solvability_level_flag=0;
      } else if (project_option_singular_possible(project_option)==1) {

       singular_possible=1; // solid regions can have coefficients all zero.

       if (local_solvability_projection==0) {
        solvability_level_flag=0;
       } else if (local_solvability_projection==1) {
        if (pbc==EXT_DIR)
         amrex::Error("cannot have outflow if solvability=1");
       } else
        amrex::Error("local_solvability_projection invalid");

      } else
       amrex::Error("project_option became corrupt");

     } else if (project_option_singular_possible(project_option)==0) {

      singular_possible=0; // diagonally dominant everywhere.
      solvability_level_flag=0;

     } else
      amrex::Error("project_option_singular_possible bad allocate_FACE_WEIGHT");

     j++;
    }  // dir
   } // side
  } // nn
 } // gridno

 int GFM_flag=0;
 int adjust_temperature=-1;  // adjust faceheat_index
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

  // adjust faceheat_index if thermal diffusion
  // and phase change using sharp interface method.
 if (project_option==2) { 

  for (int im=0;im<2*nten;im++) {
   if (latent_heat[im]!=0.0) {
    if ((freezing_model[im]==0)|| // fully saturated
        (freezing_model[im]==5)||
        (freezing_model[im]==6))  // Palmore and Desjardins
     GFM_flag=1; 
   } else if (latent_heat[im]==0.0) {
    // do nothing
   } else
    amrex::Error("latent_heat[im] invalid");
  } // im=0..2 nten -1

 } else if ((project_option>=100)&&
            (project_option<100+num_species_var)) {

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
  } // im=0..2 nten -1
 } else if (project_option_is_valid(project_option)==1) {
  // do nothing
 } else
  amrex::Error("project_option invalid31");

 if (GFM_flag==1) {
  stefan_solver_init(
   localMF[CELL_DEN_MF],
   adjust_temperature,
   project_option);
 }

 const Real* dx = geom.CellSize();

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  new_localMF(FACE_WEIGHT_MF+dir,nsolveMM_FACE,0,dir);
 } // dir
 new_localMF(OFF_DIAG_CHECK_MF,nsolveMM,0,-1);

 int local_face_index=faceden_index;  // 1/rho
 if (project_option_projection(project_option)==1) {
  local_face_index=faceden_index;  // 1/rho
 } else if (project_option==12) {  // pressure extension
   // 1/rho (only used in eval_face_coeff for sanity check purposes)
  local_face_index=faceden_index;  
 } else if (project_option==2) {
  local_face_index=faceheat_index; // thermal conductivity "k"
 } else if (project_option==3) {
  local_face_index=facevisc_index; // viscosity "mu"
 } else if ((project_option>=100)&&
            (project_option<100+num_species_var)) { // rho D
  local_face_index=facespecies_index+project_option-100;
 } else if (project_option==200) {
  local_face_index=smoothing_index;
 } else
  amrex::Error("project_option invalid allocate_FACE_WEIGHT");

 setVal_localMF(OFF_DIAG_CHECK_MF,0.0,0,nsolveMM,0);

 int mm_areafrac_index=FACE_VAR_MF;
 int mm_cell_areafrac_index=SLOPE_RECON_MF;
 if (num_materials_face==nmat) {
  mm_areafrac_index=FACEFRAC_SOLVE_MM_MF;
  mm_cell_areafrac_index=CELLFRAC_MM_MF;
 } else if (num_materials_face==1) {
  // do nothing
 } else
  amrex::Error("num_materials_face invalid");

 // (ml,mr,2) frac_pair(ml,mr), dist_pair(ml,mr)  
 int nfacefrac=nmat*nmat*2; 
 // im_inside,im_outside,3+sdim -->
 //   area, dist_to_line, dist, line normal.
 int ncellfrac=nmat*nmat*(3+AMREX_SPACEDIM);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++)
  debug_ngrow(mm_areafrac_index+dir,0,111);
 debug_ngrow(mm_cell_areafrac_index,0,133);

 Vector<int> solvability_level_flag_arr;
 solvability_level_flag_arr.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  solvability_level_flag_arr[tid]=solvability_level_flag;
 }

 for (int facewt_iter=0;facewt_iter<=1;facewt_iter++) {

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
   int bfact=parent->Space_blockingFactor(level);

   FArrayBox& slopefab=(*localMF[SLOPE_RECON_MF])[mfi];
   FArrayBox& cenden=(*localMF[CELL_DEN_MF])[mfi];  // 1/rho
   FArrayBox& cenvisc=(*localMF[CELL_VISC_MF])[mfi];

   FArrayBox& xfwt=(*localMF[FACE_WEIGHT_MF])[mfi];
   FArrayBox& yfwt=(*localMF[FACE_WEIGHT_MF+1])[mfi];
   FArrayBox& zfwt=(*localMF[FACE_WEIGHT_MF+AMREX_SPACEDIM-1])[mfi];

   FArrayBox& offdiagcheck=(*localMF[OFF_DIAG_CHECK_MF])[mfi];
 
   FArrayBox& xface=(*localMF[FACE_VAR_MF])[mfi];  
   FArrayBox& yface=(*localMF[FACE_VAR_MF+1])[mfi];  
   FArrayBox& zface=(*localMF[FACE_VAR_MF+AMREX_SPACEDIM-1])[mfi];  

   FArrayBox& xfacemm=(*localMF[mm_areafrac_index])[mfi];  
   FArrayBox& yfacemm=(*localMF[mm_areafrac_index+1])[mfi];  
   FArrayBox& zfacemm=(*localMF[mm_areafrac_index+AMREX_SPACEDIM-1])[mfi];  
   FArrayBox& cellfracmm=(*localMF[mm_cell_areafrac_index])[mfi];  

   FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi];  // mask=1 at fine-fine bc
   const Real* xlo = grid_loc[gridno].lo();

// want face_frac=0 if presbc<> interior or exterior dirichlet.
// solvability_level_flag=0 if coarse/fine interface found.
   Vector<int> presbc;
   getBCArray_list(presbc,state_index,gridno,scomp,ncomp);
   if (presbc.size()==2*AMREX_SPACEDIM*nsolveMM) {
    // do nothing
   } else
    amrex::Error("presbc.size() invalid");

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // BUILDFACEWT is defined in LEVELSET_3D.F90
   FORT_BUILDFACEWT(
    &facewt_iter,
    &num_materials_face,
    &level,
    &finest_level,
    &nsolve,
    &nsolveMM,
    &nsolveMM_FACE,
    &nfacefrac,
    &ncellfrac,
    &local_face_index,
    &facecut_index,
    &icefacecut_index,
    &ncphys,
    &nmat,
    xlo,dx,
    offdiagcheck.dataPtr(),
    ARLIM(offdiagcheck.loVect()),ARLIM(offdiagcheck.hiVect()),
    slopefab.dataPtr(),ARLIM(slopefab.loVect()),ARLIM(slopefab.hiVect()),
    cenden.dataPtr(),ARLIM(cenden.loVect()),ARLIM(cenden.hiVect()),
    cenvisc.dataPtr(),ARLIM(cenvisc.loVect()),ARLIM(cenvisc.hiVect()),
    cellfracmm.dataPtr(),
    ARLIM(cellfracmm.loVect()),ARLIM(cellfracmm.hiVect()),
    xfacemm.dataPtr(),ARLIM(xfacemm.loVect()),ARLIM(xfacemm.hiVect()),
    yfacemm.dataPtr(),ARLIM(yfacemm.loVect()),ARLIM(yfacemm.hiVect()),
    zfacemm.dataPtr(),ARLIM(zfacemm.loVect()),ARLIM(zfacemm.hiVect()),
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
    &singular_possible,
    &local_solvability_projection, // input to the fortran
    &solvability_level_flag_arr[tid_current], // declared inout in the fortran
    presbc.dataPtr(),
    &visc_coef,
    &constant_viscosity,
    &project_option);

  }  // mfi
} // omp
  ns_reconcile_d_num(184);

 } // facewt_iter=0..1

 for (int tid=1;tid<thread_class::nthreads;tid++) {

  if (solvability_level_flag_arr[tid]<solvability_level_flag_arr[0])
   solvability_level_flag_arr[0]=solvability_level_flag_arr[tid];

  for (int iwt=0;iwt<4;iwt++) {
   if (min_face_wt[tid][iwt]<min_face_wt[0][iwt])
    min_face_wt[0][iwt]=min_face_wt[tid][iwt];
   if (max_face_wt[tid][iwt]>max_face_wt[0][iwt])
    max_face_wt[0][iwt]=max_face_wt[tid][iwt];
  }
 } // tid

 solvability_level_flag=solvability_level_flag_arr[0];
 ParallelDescriptor::Barrier();

 ParallelDescriptor::ReduceIntMin(solvability_level_flag);

 for (int iwt=0;iwt<4;iwt++) {
  ParallelDescriptor::ReduceRealMin(min_face_wt[0][iwt]);
  ParallelDescriptor::ReduceRealMax(max_face_wt[0][iwt]);
 }
 for (int tid=1;tid<thread_class::nthreads;tid++) {
  for (int iwt=0;iwt<4;iwt++) {
   min_face_wt[tid][iwt]=min_face_wt[0][iwt];
   max_face_wt[tid][iwt]=max_face_wt[0][iwt];
  }
 }

} // subroutine allocate_FACE_WEIGHT 

// called from NavierStokes::multiphase_project
void NavierStokes::allocate_project_variables(int nsolve,int project_option) {
 
 int nmat=num_materials;

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int num_materials_face=num_materials_vel;
 if (project_option_momeqn(project_option)==1) {
  if (num_materials_face!=1)
   amrex::Error("num_materials_face invalid");
 } else if (project_option_momeqn(project_option)==0) {
  num_materials_face=num_materials_scalar_solve;
 } else
  amrex::Error("project_option_momeqn invalid32");

 if (dt_slab<=0.0)
  amrex::Error("cannot have dt_slab<=0 in allocate_project_variables");
 debug_ngrow(FACE_VAR_MF,0,850);

 int state_index;
 Vector<int> scomp;
 Vector<int> ncomp;
 int ncomp_check;
 get_mm_scomp_solver(
  num_materials_face,
  project_option,
  state_index,
  scomp,
  ncomp,
  ncomp_check);
 int nsolveMM=nsolve*num_materials_face;
 if (ncomp_check!=nsolveMM)
  amrex::Error("nsolveMM invalid 3919");
 
 MultiFab& S_new=get_new_data(state_index,slab_step+1);

 new_localMF(ONES_MF,num_materials_face,0,-1);
 setVal_localMF(ONES_MF,1.0,0,num_materials_face,0);
 ones_sum_global=0.0;

  // allocates and initializes DOTMASK_MF (will be used for
  // next generation heat equation solver)
 makeDotMask(nsolve,project_option);

 new_localMF(POLDHOLD_MF,nsolveMM,0,-1);
 setVal_localMF(POLDHOLD_MF,0.0,0,nsolveMM,0);

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
 new_localMF(OUTER_ITER_PRESSURE_MF,nsolveMM,0,-1);

 int adjust_temperature=1; 
 int GFM_flag=0;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

  // temperature diffusion
 if (project_option==2) {

    // MEHDI VAHAB HEAT SOURCE
    // T^new=T^* += dt A Q/(rho cv V) 
    // in: allocate_project_variables
    // NavierStokes::heat_source_term_flux_source  (in:NavierStokes.cpp)
    // heat_source_term_flux_source calls GODUNOV_3D::FORT_HEATSOURCE_FACE
  heat_source_term_flux_source();

  if (is_phasechange==1) {
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
   for (int im=0;im<2*nten;im++) {
    if (latent_heat[im]!=0.0) {
     if ((freezing_model[im]==0)|| // fully saturated.
         (freezing_model[im]==5)||
         (freezing_model[im]==6))  // Palmore and Desjardins
      GFM_flag=1;
    } else if (latent_heat[im]==0.0) {
     // do nothing
    } else
     amrex::Error("latent_heat[im] invalid");
   } // im=0..2 nten -1
  } else if (is_phasechange==0) {
   // do nothing
  } else
   amrex::Error("is_phasechange invalid");

 } // project_option==2

 if ((project_option>=100)&&
     (project_option<100+num_species_var)) {

  if (is_phasechange==1) {

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
 } else if (state_index==DIV_Type) {
  if (scomp.size()!=1)
   amrex::Error("scomp.size() invalid");
  if (ncomp[0]!=num_materials_vel)
   amrex::Error("ncomp[0] invalid");
  current_contents_mf=getStateDIV_DATA(1,scomp[0],ncomp[0],cur_time_slab);
 } else {
  current_contents_mf=nullptr;
  amrex::Error("state_index invalid");
 }

  // ``OUTER_ITER_PRESSURE'' = S_new = S^*
 MultiFab::Copy(*localMF[OUTER_ITER_PRESSURE_MF],
		*current_contents_mf,0,0,nsolveMM,0);

 if (project_option_olddata_needed(project_option)==1) { 
  
   // if S^initial <> S^*, then put something else into
   // ``initial_guess.''   
   // outer_iter_pressure=S_new=S^*
   // For now, initial_guess=S^*
  MultiFab* initial_guess;
  initial_guess=localMF[OUTER_ITER_PRESSURE_MF]; 
  
  MultiFab* dp=new MultiFab(grids,dmap,nsolveMM,0,
   MFInfo().SetTag("dp"),FArrayBoxFactory());

  MultiFab::Copy(*dp,*initial_guess,0,0,nsolveMM,0);
   // dp=initial_guess - S^*   (S_new=S^*)
  MultiFab::Subtract(*dp,*current_contents_mf,0,0,nsolveMM,0);
   // snew+=(initial_guess - S^*)=initial_guess  (S_new=S^* beforehand)
  MultiFab::Add(*current_contents_mf,*dp,0,0,nsolveMM,0);

  int scomp_temp=0;
  for (int ilist=0;ilist<scomp.size();ilist++) {
   MultiFab::Copy(S_new,*current_contents_mf,
		  scomp_temp,scomp[ilist],ncomp[ilist],0);
   scomp_temp+=ncomp[ilist];
  }
  if (scomp_temp!=nsolveMM)
   amrex::Error("scomp_temp invalid");

   // dp=S^init-S^*   dS=S-S^init
   // alpha(S^init + dS - S^*) - div beta grad (S^init + dS) = 0
   // alpha dS - div beta grad dS = -alpha dp + div beta grad S^init 
   //                             = alpha POLDHOLD + div beta grad S^init
   // OUTER_ITER_PRESSURE=S^* + (S^init - S^*)=S^init
  MultiFab::Add(*localMF[OUTER_ITER_PRESSURE_MF],*dp,0,0,nsolveMM,0);
   // POLDHOLD=0 - (S^init-S^*) = S^* - S^init
  MultiFab::Subtract(*localMF[POLDHOLD_MF],*dp,0,0,nsolveMM,0);
   // later on, (1) UMAC_MF-=beta grad S^init,  (S_new=S^init)
   //           (2) S_new=0.0
  delete dp;

 } else if (project_option_olddata_needed(project_option)==0) { 
  // do nothing
 } else
  amrex::Error("project_option_olddata_needed invalid33");

 delete current_contents_mf;

} // subroutine allocate_project_variables


void NavierStokes::allocate_pressure_work_vars(int nsolve,int project_option) {

 int nmat=num_materials;

 int num_materials_face=num_materials_vel;
 if (project_option_momeqn(project_option)==1) {
  if (num_materials_face!=1)
   amrex::Error("num_materials_face invalid");
 } else if (project_option_momeqn(project_option)==0) {
  num_materials_face=num_materials_scalar_solve;
 } else
  amrex::Error("project_option invalid34");

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 if ((num_materials_face!=1)&&
     (num_materials_face!=nmat))
  amrex::Error("num_materials_face invalid");

 if ((nsolve!=1)&&(nsolve!=AMREX_SPACEDIM))
  amrex::Error("nsolve invalid");
 int nsolveMM=nsolve*num_materials_face;
 int nsolveMM_FACE=nsolveMM;
 if (num_materials_face==1) {
  // do nothing
 } else if (num_materials_face==nmat) {
  nsolveMM_FACE*=2;
 } else
  amrex::Error("num_materials_face invalid");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  new_localMF(UMACSTAR_MF+dir,nsolveMM_FACE,0,dir);
  new_localMF(GRADPEDGE_MF+dir,nsolveMM_FACE,0,dir);

   // PEDGE_MF only used if pressure projection.
   // 0=use_face_pres  1=grid flag  2+im_vel (im_vel=0..nsolveMM_FACE-1) pface
  new_localMF(PEDGE_MF+dir,2+nsolveMM_FACE,0,dir);

  new_localMF(AMRSYNC_PRES_MF+dir,nsolveMM_FACE,0,dir);
  new_localMF(AMRSYNC_PEDGE_MF+dir,1,0,dir);
  setVal_localMF(AMRSYNC_PRES_MF+dir,1.0e+40,0,nsolveMM_FACE,0);
  setVal_localMF(AMRSYNC_PEDGE_MF+dir,1.0e+40,0,1,0);
 } // dir=0..sdim-1

} // subroutine allocate_pressure_work_vars

// called in the setup stage from NavierStokes::multiphase_project when
// project_option==0 
void NavierStokes::overwrite_outflow() {
 
 bool use_tiling=ns_tiling;

 int nsolveMM_FACE=num_materials_vel;

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 int scomp_pres=num_materials_vel*AMREX_SPACEDIM;

 const Real* dx = geom.CellSize();
 const Real* prob_lo   = geom.ProbLo();
 const Real* prob_hi   = geom.ProbHi();
 MultiFab& U_new = get_new_data(State_Type,slab_step+1); 

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   MultiFab& Umac_new=get_new_data(Umac_Type+dir,slab_step+1);
   if (Umac_new.nComp()!=nsolveMM_FACE)
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
 
    Vector<int> presbc=getBCArray(State_Type,gridno,scomp_pres,
     num_materials_vel);

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     // in: PROB.F90
    FORT_FORCEVELOCITY(
      &nsolveMM_FACE,
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
   ns_reconcile_d_num(185);

 } // dir=0..sdim-1

} // subroutine overwrite_outflow


// macdest=macsrc+gp
// for pressure projection: 
// macdest=UMAC,UMACSTAR,UMAC
// macsrc =UMAC,UMACSTAR,MAC_TEMP
void NavierStokes::correct_velocity(
  int project_option,
  int macdest,
  int macsrc,
  int gp,int nsolve) {

 int finest_level=parent->finestLevel();

 int nmat=num_materials;

 int num_materials_face=num_materials_vel;
 if (project_option_momeqn(project_option)==1) {
  if (num_materials_face!=1)
   amrex::Error("num_materials_face invalid");
 } else if (project_option_momeqn(project_option)==0) {
  num_materials_face=num_materials_scalar_solve;
 } else
  amrex::Error("project_option invalid35");

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 if ((num_materials_face!=1)&&
     (num_materials_face!=nmat))
  amrex::Error("num_materials_face invalid");

 int nsolveMM=nsolve*num_materials_face;
 int nsolveMM_FACE=nsolveMM;
 if (num_materials_face==1) {
  // do nothing
 } else if (num_materials_face==nmat) {
  nsolveMM_FACE*=2;
 } else
  amrex::Error("num_materials_face invalid");

 Vector<int> scomp;
 Vector<int> ncomp;
 int state_index;
 int ncomp_check;

 get_mm_scomp_solver(
  num_materials_face,
  project_option,
  state_index,
  scomp,
  ncomp,
  ncomp_check);

 if (ncomp_check!=nsolveMM)
  amrex::Error("ncomp_check invalid");

 bool use_tiling=ns_tiling;

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  if ((localMF[macdest+dir]->nComp()!=nsolveMM_FACE)||
      (localMF[macsrc+dir]->nComp()!=nsolveMM_FACE)||
      (localMF[gp+dir]->nComp()!=nsolveMM_FACE))
   amrex::Error("invalid ncomp");
 }

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 if (project_option_is_valid(project_option)==1) {
  // do nothing
 } else
  amrex::Error("project option invalid");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  debug_ngrow(FACE_VAR_MF+dir,0,862);
 }

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,1,6001);

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
  if (presbc.size()!=nsolveMM*AMREX_SPACEDIM*2)
   amrex::Error("presbc.size() invalid");

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  for (int velcomp=0;velcomp<nsolveMM_FACE;velcomp++) {

    // in: NAVIERSTOKES_3D.F90
   FORT_FLUIDSOLIDCOR(
    &level,
    &finest_level,
    &velcomp,
    &nsolve,
    &nsolveMM,
    &nsolveMM_FACE,
    &facecut_index,
    &icefacecut_index,
    &ncphys,
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
    xlo,dx,
    &cur_time_slab,
    &nmat);
  } // velcomp=0..nsolveMM_FACE
 } // mfi
} // omp
 ns_reconcile_d_num(186);

} // subroutine correct_velocity


void NavierStokes::residual_correction_form(
  int homflag_residual_correction_form,
  int energyflag,
  int project_option,int nsolve) {

 int nmat=num_materials;

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int num_materials_face=num_materials_vel;
 if (project_option_momeqn(project_option)==1) {
  if (num_materials_face!=1)
   amrex::Error("num_materials_face invalid");
 } else if (project_option_momeqn(project_option)==0) {
  num_materials_face=num_materials_scalar_solve;
 } else
  amrex::Error("project_option_momeqn invalid36");

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 if ((num_materials_face!=1)&&
     (num_materials_face!=nmat))
  amrex::Error("num_materials_face invalid");

 if (project_option==1) { // initial project
  //do nothing
 } else if (project_option_is_valid(project_option)==1) {

   // -dt grad p face_weight  
   // SEM BC if projection_enable_spectral==1,2
  int simple_AMR_BC_flag=0;
  int simple_AMR_BC_flag_viscosity=0;
  apply_pressure_grad(
   simple_AMR_BC_flag,
   simple_AMR_BC_flag_viscosity,
   homflag_residual_correction_form,
   energyflag,
   GRADPEDGE_MF,
   STATE_FOR_RESID_MF,
   project_option,nsolve);

   // UMAC_MF-=GRADPEDGE_MF
  correct_velocity(project_option,
   UMAC_MF, UMAC_MF, GRADPEDGE_MF,nsolve);

 } else
  amrex::Error("project_option invalid residual_correction_form");

}  // residual_correction_form


// local_MF[idx_phi]=0 on all levels on input
void NavierStokes::mg_cycleALL(int presmooth,
 int project_option,
 int idx_rhs,int idx_phi,
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

 int nmat=num_materials;

 int num_materials_face=num_materials_vel;
 if (project_option_momeqn(project_option)==1) {
  if (num_materials_face!=1)
   amrex::Error("num_materials_face invalid");
 } else if (project_option_momeqn(project_option)==0) {
  num_materials_face=num_materials_scalar_solve;
 } else
  amrex::Error("project_option_momeqn(project_option) invalid37");

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 if ((num_materials_face!=1)&&
     (num_materials_face!=nmat))
  amrex::Error("num_materials_face invalid");

 int nsolveMM=nsolve*num_materials_face;
 int nsolveMM_FACE=nsolveMM;
 if (num_materials_face==1) {
  // do nothing
 } else if (num_materials_face==nmat) {
  nsolveMM_FACE*=2;
 } else
  amrex::Error("num_materials_face invalid");

   // residmf represents the current status of:
   // f+ div grad p^init + a(POLDHOLD)  (POLDHOLD=p^*-p^init)
   // originally:
   // a(p-p*)-div grad p = f
   // let dp=p-p^init
   // a(dp+p^init-p^*)- div grad (dp+p^init) = f 
   // a dp - div grad dp = f+a(p^*-p^init)+div grad p^init=f+a POLDHOLD +
   // div grad p^init
 MultiFab* residmf=new MultiFab(grids,dmap,nsolveMM,0,
	MFInfo().SetTag("residmf"),FArrayBoxFactory());
 MultiFab::Copy(*residmf,*localMF[idx_rhs],0,0,nsolveMM,0);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  setVal_localMF(UMACSTAR_MF+dir,0.0,0,nsolveMM_FACE,0);
 } 

#if (profile_solver==1)
 bprof.stop();
#endif

 relaxLEVEL(residmf,idx_rhs,idx_phi,
   presmooth,
   project_option,nsolve);

 delete residmf;
} // subroutine mg_cycleALL

// this recursive routine first called from the finest_level.
// localMF[idx_phi]=0.0 on all levels prior to the first call of this
// routine.
void NavierStokes::relaxLEVEL(
  MultiFab* rhsmf,
  int idx_rhs,int idx_phi,
  int presmooth,
  int project_option,int nsolve) {

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

 int nmat=num_materials;

 int num_materials_face=num_materials_vel;
 if (project_option_momeqn(project_option)==1) {
  if (num_materials_face!=1)
   amrex::Error("num_materials_face invalid");
 } else if (project_option_momeqn(project_option)==0) {
  num_materials_face=num_materials_scalar_solve;
 } else
  amrex::Error("project_option invalid38");

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 if ((num_materials_face!=1)&&
     (num_materials_face!=nmat))
  amrex::Error("num_materials_face invalid");

 int nsolveMM=nsolve*num_materials_face;
 int nsolveMM_FACE=nsolveMM;
 if (num_materials_face==1) {
  // do nothing
 } else if (num_materials_face==nmat) {
  nsolveMM_FACE*=2;
 } else
  amrex::Error("num_materials_face invalid");

 MultiFab* pbdry=new MultiFab(grids,dmap,nsolveMM,1,
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
  nsolve);

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
  ns_coarse.setVal_localMF(idx_phi,0.0,0,nsolveMM,1);

  for (int i=presmooth;i>0;i--) {
   int apply_lev_presmooth=0;
   mac_op->smooth(*localMF[idx_phi],*rhsmf,
    apply_lev_presmooth,*pbdry,bcpres_array,smooth_type);
  }
  MultiFab* residmf=new MultiFab(grids,dmap,nsolveMM,0,
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
  int energyflag=0;

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
   project_option,nsolve);

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
   nsolve);

   // UMACSTAR=UMACSTAR+GRADPEDGE
  correct_velocity(project_option,
    UMACSTAR_MF, UMACSTAR_MF, GRADPEDGE_MF,nsolve);

  MultiFab* residmf_coarse=
   new MultiFab(ns_coarse.grids,ns_coarse.dmap,nsolveMM,0,
	MFInfo().SetTag("residmf_coarse"),FArrayBoxFactory());
  residmf_coarse->setVal(0.0,0,nsolveMM,0);

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
   ns_coarse.setVal_localMF(UMACSTAR_MF+dir,0.0,0,nsolveMM_FACE,0);

  int ncomp_edge=-1;
  int scomp_edge=0;
  int start_dir=0;
  int spectral_override=1; // order determined from enable_spectral
   // average down from level to level-1.
  ns_coarse.avgDownEdge_localMF(
    UMACSTAR_MF,
    scomp_edge,ncomp_edge,
    start_dir,AMREX_SPACEDIM,spectral_override,15);

  int iavg=0;
  BoxArray& fgridscen=grids;
  DistributionMapping& fdmap=dmap;
  BoxArray& cgridscen=ns_coarse.grids;
  for (int veldir=0;veldir<nsolveMM;veldir++) {
    // average down from level to level-1.
    // calls FORT_AVERAGE which is low order.
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
  Real A_error0=0.0;
  Real rAr_error0=0.0;
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
   A_error0,
   rAr_error0,
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
        Real Ar_error_n,Real Ar_error0,
        Real rAr_error_n,Real rAr_error0,
        Real save_mac_abs_tol_in,
        int& meets_tol) {

 if (error_n<save_mac_abs_tol_in) {
  meets_tol=1;
 } else if (Ar_error_n==0.0) {
  meets_tol=1;
 } else if (rAr_error_n<=0.0) {
  meets_tol=1;
 } else if ((error_n>0.0)&&
   	    (Ar_error_n>0.0)&&
	    (rAr_error_n>0.0)&&
	    (error0>0.0)&&
	    (Ar_error0>0.0)&&
	    (rAr_error0>0.0)) {
  Real scale_Ar_error_n=Ar_error_n*error0/Ar_error0;
  Real scale_rAr_error_n=rAr_error_n*error0/rAr_error0;
  if (scale_Ar_error_n<save_mac_abs_tol_in) {
   meets_tol=1;
  } else if (scale_rAr_error_n<save_mac_abs_tol_in) {
   meets_tol=1;
  } else if ((scale_Ar_error_n>=save_mac_abs_tol_in)&&
             (scale_rAr_error_n>=save_mac_abs_tol_in)) {
   // do nothing
  } else
   amrex::Error("scale_Ar_error_n or scale_rAr_error_n invalid");
 } else
  amrex::Error("one of error_r,Ar,rAr invalid");

} // subroutine check_outer_solver_convergence

// update_vel=1 if called at the beginning of each outer_iter sweep.
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
 int bicgstab_num_outer_iterSOLVER,int nsolve) {

 int finest_level=parent->finestLevel();

 if (level==0) {
  // do nothing
 } else
  amrex::Error("jacobi_cycles should only be called from level==0");

 error_at_the_beginning=0.0;
 error_after_all_jacobi_sweeps=0.0;

 if ((nsolve!=1)&&(nsolve!=AMREX_SPACEDIM))
  amrex::Error("nsolve invalid33");

 int nmat=num_materials;

 int num_materials_face=num_materials_vel;
 if (project_option_momeqn(project_option)==1) {
  if (num_materials_face!=1)
   amrex::Error("num_materials_face invalid");
 } else if (project_option_momeqn(project_option)==0) {
  num_materials_face=num_materials_scalar_solve;
 } else
  amrex::Error("project_option invalid39");

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 if ((num_materials_face!=1)&&
     (num_materials_face!=nmat))
  amrex::Error("num_materials_face invalid");

 int nsolveMM=nsolve*num_materials_face;

 if (bicgstab_num_outer_iterSOLVER<0)
  amrex::Error("bicgstab_num_outer_iterSOLVER invalid");

 allocate_array(0,nsolveMM,-1,RESID_MF);

 int temp_ncycles=ncycles;
 if (ncycles==0) 
  temp_ncycles=1;

 for (int vcycle_jacobi=0;vcycle_jacobi<temp_ncycles;vcycle_jacobi++) {

    // if local_solvability_projection, then this
    // routine modifies RESID so that it sums to zero.
    // RESID=project( mac_rhs-(alpha*phi+div (-grad p)/dt) )
    // if singular_possible, then this routine zeros out the
    // residual where the matrix diagonal (prior to dual
    // time stepping modification) is 0.
  residALL(project_option,idx_mac_rhs_crse,
    RESID_MF,idx_mac_phi_crse,nsolve);

  if (update_vel==1) {  // not called as a preconditioner

   Real local_error;
   dot_productALL(project_option,RESID_MF,RESID_MF,local_error,nsolve);
   local_error=sqrt(local_error);
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

   if (bicgstab_num_outer_iterSOLVER==0) {

    error0=local_error;

    if (call_adjust_tolerance==1) {
     adjust_tolerance(error0,error0_max,project_option);
    } else if (call_adjust_tolerance==0) {
     // do nothing
    } else
     amrex::Error("call_adjust_tolerance invalid");

   } else if (bicgstab_num_outer_iterSOLVER>0) {
    // do nothing
   } else
    amrex::Error("bicgstab_num_outer_iterSOLVER invalid");

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

    // if local_solvability_projection, then this
    // routine modifies RESID so that it sums to zero.
    // RESID=project( mac_rhs-(alpha*phi+div (-grad p)/dt) )
    // if singular_possible, then this routine zeros out the
    // residual where the matrix diagonal (prior to dual time stepping
    // modification) is 0.
  residALL(project_option,idx_mac_rhs_crse,
    RESID_MF,idx_mac_phi_crse,nsolve);

  if (update_vel==1) {  // not called as a preconditioner

   Real local_error;
   dot_productALL(project_option,RESID_MF,RESID_MF,local_error,nsolve);
   local_error=sqrt(local_error);
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

}  // jacobi_cycles

// called from:
//  NavierStokes::multiphase_project
void NavierStokes::updatevelALL(
 int project_option,
 int idx_mac_phi_crse,int nsolve) {

 int finest_level=parent->finestLevel();

 int nmat=num_materials;

 int num_materials_face=num_materials_vel;

 if (project_option_momeqn(project_option)==1) {
  if (num_materials_face!=1)
   amrex::Error("num_materials_face invalid");
 } else if (project_option_momeqn(project_option)==0) {
  num_materials_face=num_materials_scalar_solve;
 } else
  amrex::Error("project_option_momeqn(project_option)  invalid40");

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 if ((num_materials_face!=1)&&
     (num_materials_face!=nmat))
  amrex::Error("num_materials_face invalid");

 int nsolveMM=nsolve*num_materials_face;

   // gradpedge=-dt W grad p
 applyGradALL(project_option,idx_mac_phi_crse,nsolve);

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
     // POLDHOLD-=idx_mac_phi_crse
     // PNEW+=idx_mac_phi_crse
     // UMAC+=GRADPEDGE  (GRADPEDGE=-dt W grad p)
  ns_level.mac_update(ns_level.localMF[idx_mac_phi_crse],
    project_option,nsolve);
  ns_level.setVal_localMF(idx_mac_phi_crse,0.0,0,nsolveMM,1);

  if (ilev<finest_level) {
   ns_level.avgDownMac();   // works on UMAC_MF
  }
 }

} // end subroutine updatevelALL



void NavierStokes::Prepare_UMAC_for_solver(int project_option,
  int nsolve) {

 if (dt_slab<=0.0)
  amrex::Error("dt_slab invalid4");
 if ((nsolve!=1)&&(nsolve!=AMREX_SPACEDIM))
  amrex::Error("nsolve invalid");

 int nmat=num_materials;
 int finest_level=parent->finestLevel();

 if ((level>=0)&&(level<=finest_level)) {
  // do nothing
 } else
  amrex::Error("level invalid Prepare_UMAC_for_solver");

 int num_materials_face=num_materials_vel;
 if (project_option_momeqn(project_option)==1) {
  if (num_materials_face!=1)
   amrex::Error("num_materials_face invalid");
 } else if (project_option_momeqn(project_option)==0) {
  num_materials_face=num_materials_scalar_solve;
 } else
  amrex::Error("project_option_momeqn invalid41");

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 if ((num_materials_face!=1)&&
     (num_materials_face!=nmat))
  amrex::Error("num_materials_face invalid");

 int nsolveMM=nsolve*num_materials_face;
 int nsolveMM_FACE=nsolveMM;
 if (num_materials_face==1) {
  // do nothing
 } else if (num_materials_face==nmat) {
  nsolveMM_FACE*=2;
 } else
  amrex::Error("num_materials_face invalid");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  new_localMF(MAC_TEMP_MF+dir,nsolveMM_FACE,0,dir);
  setVal_localMF(MAC_TEMP_MF+dir,0.0,0,nsolveMM_FACE,0);
 } // dir

 new_localMF(DIFFUSIONRHS_MF,nsolveMM,0,-1);
 if ((project_option==11)|| //FSI_material_exists last project
     (project_option==12)) {//pressure extension
  setVal_localMF(DIFFUSIONRHS_MF,0.0,0,nsolveMM,0);
 } else if (project_option==0)  { // regular project
  if (num_materials_face!=1)
   amrex::Error("num_materials_face invalid");
  int scomp=0;
   // MDOT_MF already premultiplied by the cell volume
  Copy_localMF(DIFFUSIONRHS_MF,MDOT_MF,0,scomp,nsolveMM,0);
 } else if (project_option==1) { // initial project
  if (num_materials_face!=1)
   amrex::Error("num_materials_face invalid");
  int scomp=0;
   // MDOT_MF already premultiplied by the cell volume
  Copy_localMF(DIFFUSIONRHS_MF,MDOT_MF,0,scomp,nsolve,0);
  zero_independent_variable(project_option,nsolve);
 } else if (project_option==2) { // thermal conduction

  zero_independent_vel(project_option,UMAC_MF,nsolve);
  setVal_localMF(DIFFUSIONRHS_MF,0.0,0,nsolveMM,0);

 } else if (project_option==3) { // viscosity

  zero_independent_vel(project_option,UMAC_MF,nsolve);
  setVal_localMF(DIFFUSIONRHS_MF,0.0,0,nsolveMM,0);

 } else if ((project_option>=100)&&
	    (project_option<100+num_species_var)) {//species

  zero_independent_vel(project_option,UMAC_MF,nsolve);
  setVal_localMF(DIFFUSIONRHS_MF,0.0,0,nsolveMM,0);

 } else if (project_option==200) {//smoothing

  zero_independent_vel(project_option,UMAC_MF,nsolve);
  setVal_localMF(DIFFUSIONRHS_MF,0.0,0,nsolveMM,0);

 } else
  amrex::Error("project_option invalid prepare_UMAC_for_solver");

} // subroutine Prepare_UMAC_for_solver

void NavierStokes::remove_UMAC_for_solver(int project_option) {

 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid remove_UMAC_for_solver");

 delete_localMF(DIFFUSIONRHS_MF,1);
 delete_localMF(MAC_TEMP_MF,AMREX_SPACEDIM);

}

void NavierStokes::multiphase_GMRES_preconditioner(
 int gmres_precond_iter,
 int project_option,int project_timings,
 int presmooth,int postsmooth,
 int idx_Z,int idx_R,int nsolve) {

 int caller_id=0;

 if (level==0) {
  // do nothing
 } else
  amrex::Error("level invalid");

 if (override_bc_to_homogeneous==1) {
  // do nothing
 } else
  amrex::Error("expecting homogeneous mode in preconditioner");

 int num_materials_face=num_materials_vel;
 if (project_option_momeqn(project_option)==1) {
  if (num_materials_face!=1)
   amrex::Error("num_materials_face invalid");
 } else if (project_option_momeqn(project_option)==0) {
  num_materials_face=num_materials_scalar_solve;
 } else
  amrex::Error("project_option invalid44");

 int nmat=num_materials; 
 if ((num_materials_face!=1)&&
     (num_materials_face!=nmat))
  amrex::Error("num_materials_face invalid");
 int nsolveMM=nsolve*num_materials_face;

 int change_flag=0;
 project_right_hand_side(idx_R,project_option,change_flag);

 if (gmres_precond_iter==0) {
   // Z=M^{-1}R
   // Z=project(Z)
   // first calls: zeroALL(1,nsolveMM,idx_Z) 
  multiphase_preconditioner(
   project_option,project_timings,
   presmooth,postsmooth,
   idx_Z,idx_R,nsolve);
 } else if ((gmres_precond_iter>0)&&
            (gmres_precond_iter*nsolveMM<=MAX_GMRES_BUFFER)) {
  int m=gmres_precond_iter*nsolveMM;

  Real GMRES_tol=KRYLOV_NORM_CUTOFF;
  Real beta=0.0;
  dot_productALL(project_option,idx_R,idx_R,beta,nsolve);

  if (beta>0.0) {

   beta=sqrt(beta);

    // variables initialized to 0.0
   allocate_array(0,nsolveMM,-1,GMRES_BUFFER0_V_MF);
   allocate_array(0,nsolveMM,-1,GMRES_BUFFER_W_MF);

   Real aa=1.0/beta;
    // V0=V0+aa R = R/||R||
   mf_combine(project_option,
    GMRES_BUFFER0_V_MF,idx_R,aa,GMRES_BUFFER0_V_MF,nsolve); 
   project_right_hand_side(GMRES_BUFFER0_V_MF,project_option,change_flag);

    // H_j is a j+2 x j+1 matrix  j=0..m-1
   Real** HH=new Real*[m+1];
   for (int i=0;i<m+1;i++) { 
    HH[i]=new Real[m];
    for (int j=0;j<m;j++) 
     HH[i][j]=0.0;
   }

   Real* yy=new Real[m];
   Real* beta_e1=new Real[2*(m+1)];
   for (int j=0;j<2*(m+1);j++)
    beta_e1[j]=0.0;
   beta_e1[0]=beta;

   int status=1;

   int m_small=m;

    // for Zeyus routine:
    // inputs: 
    //   double** HH
    //   integer HH_rows,HH_cols
    //   integer HH_sub_rows,HH_sub_cols
    // Another routine for Zeyu?
    //   given B=[ H
    //             G ]
    //   find min_y ||By -beta e1||
    // inputs:
    //   double** B
    //   integer B_rows,B_cols
    //   integer B_sub_rows,B_sub_cols
    //  B will be some kind of m+p x m matrix
    //  and if the algorithm is slower than O(m^2), then
    //  no need to code, since a Gaussian elimination solver is already
    //  available.

   if (breakdown_free_flag==0) {

    for (int i=1;i<m;i++) {
     // variables initialized to 0.0
     allocate_array(0,nsolveMM,-1,GMRES_BUFFER0_V_MF+i);
    }
    for (int i=0;i<m;i++) {
     // variables initialized to 0.0
     allocate_array(1,nsolveMM,-1,GMRES_BUFFER0_Z_MF+i);
    }

    for (int j=0;j<m_small;j++) {

     if (verbose>0) {
      if (ParallelDescriptor::IOProcessor()) {
       std::cout << "GMRES project_option= " << project_option <<
	      " nsolveMM= " << nsolveMM <<
	      " loop j= " << j << " m= " << m << '\n';
      }
     }  // verbose>0

      // Z=M^{-1}Vj
      // Z=project(Z)
      // Vj is a normalized vector
     multiphase_preconditioner(
      project_option,project_timings,
      presmooth,postsmooth,
      GMRES_BUFFER0_Z_MF+j,
      GMRES_BUFFER0_V_MF+j,nsolve);

      // W=A Z
      // 1. (begin)calls project_right_hand_side(Z)
      // 2. (end)  calls project_right_hand_side(W)
     applyALL(project_option,
       GMRES_BUFFER0_Z_MF+j,
       GMRES_BUFFER_W_MF,nsolve);

     for (int i=0;i<=j;i++) {
      dot_productALL(project_option,
       GMRES_BUFFER_W_MF,
       GMRES_BUFFER0_V_MF+i,HH[i][j],nsolve);
      aa=-HH[i][j];
       // W=W+aa Vi
      mf_combine(project_option,
       GMRES_BUFFER_W_MF,GMRES_BUFFER0_V_MF+i,aa,GMRES_BUFFER_W_MF,nsolve); 
      project_right_hand_side(GMRES_BUFFER_W_MF,project_option,change_flag);
     } // i=0..j

     dot_productALL(project_option,
       GMRES_BUFFER_W_MF,
       GMRES_BUFFER_W_MF,HH[j+1][j],nsolve);

     if (HH[j+1][j]>=0.0) {
      HH[j+1][j]=sqrt(HH[j+1][j]);
     } else
      amrex::Error("HH[j+1][j] invalid");

     if (HH[j+1][j]>KRYLOV_NORM_CUTOFF) {

      status=1;

      if ((j>=0)&&(j<m_small)) {
       caller_id=11;
       GMRES_MIN_CPP(HH,beta,yy,
	m,j+1,
	caller_id,project_option,
        level,status);
      } else
       amrex::Error("j invalid");

      if (status==1) {

       if ((j>=0)&&(j<m-1)) {
        aa=1.0/HH[j+1][j];
         // V=V+aa W
        mf_combine(project_option,
         GMRES_BUFFER0_V_MF+j+1,
         GMRES_BUFFER_W_MF,aa,
         GMRES_BUFFER0_V_MF+j+1,nsolve); 
        project_right_hand_side(GMRES_BUFFER0_V_MF+j+1,
		       project_option,change_flag);
       } else if (j==m-1) {
        // do nothing
       } else
        amrex::Error("j invalid");

      } else if (status==0) {

       m_small=j;

      } else
       amrex::Error("status invalid");

     } else if ((HH[j+1][j]>=0.0)&&(HH[j+1][j]<=KRYLOV_NORM_CUTOFF)) {
      m_small=j; // valid indexes of HH: i1=0..j  i2=0..j-1
     } else
      amrex::Error("HH[j+1][j] invalid");

    } // j=0..m_small-1
  
    status=1;

    if ((m_small>=1)&&
        (m_small<=m)) {

     caller_id=1;
     GMRES_MIN_CPP(HH,beta,yy,
	m,m_small,
	caller_id,project_option,
        level,status);

     if (status==1) {
      // do nothing
     } else
      amrex::Error("expecting status==1 here");

    } else if (m_small==0) {

      // Z=M^{-1}R
      // first calls: zeroALL(1,nsolveMM,idx_Z) 
      // Z=project(Z)
     multiphase_preconditioner(
      project_option,project_timings,
      presmooth,postsmooth,
      idx_Z,idx_R,nsolve);

    } else {
     std::cout << "m_small= " << m_small << '\n';
     amrex::Error("NavierStokes3.cpp m_small invalid");
    }

    if (status==1) {
     zeroALL(1,nsolveMM,idx_Z);
     for (int j=0;j<m_small;j++) {
      aa=yy[j];
       // Z=Z+aa Zj
      mf_combine(project_option,
       idx_Z,
       GMRES_BUFFER0_Z_MF+j,aa,
       idx_Z,nsolve); 
     }
    } else
     amrex::Error("status invalid");

    for (int i=1;i<m;i++) {
     delete_array(GMRES_BUFFER0_V_MF+i); 
    }
    for (int i=0;i<m;i++) {
     delete_array(GMRES_BUFFER0_Z_MF+i); 
    }

   } else if (breakdown_free_flag==1) {

    Real breakdown_tol=KRYLOV_NORM_CUTOFF;

     // G_j is a p(j)+1 x j+1 matrix   j=0..m-1
    Real** GG=new Real*[m+1];
    for (int i=0;i<m+1;i++) { 
     GG[i]=new Real[m];
     for (int j=0;j<m;j++) 
      GG[i][j]=0.0;
    }

    Real** HHGG=new Real*[2*m+2];
    for (int i=0;i<2*m+2;i++) { 
     HHGG[i]=new Real[2*m+1];
     for (int j=0;j<2*m+1;j++) 
      HHGG[i][j]=0.0;
    }

    int convergence_flag=0;
    int p_local=-1;
    int p_local_allocated=0;
    int p_local_init_allocated=p_local_allocated;

    int j_local=0;

    if (debug_BF_GMRES==1) {
     std::cout << "BF_GMRES: project_option= " << project_option
	    << "gmres_precond_iter= " << gmres_precond_iter << '\n';
     std::cout << "BF_GMRES: beta= " << beta << '\n';
    }

    Vector<Real> error_historyGMRES;
    error_historyGMRES.resize(m);

    int use_previous_iterate=0;

    for (j_local=0;((j_local<m_small)&&(convergence_flag==0));j_local++) {

     use_previous_iterate=0;

     // variables initialized to 0.0
     allocate_array(1,nsolveMM,-1,GMRES_BUFFER0_Z_MF+j_local);
      // Zj=M^{-1}Vj
      // Z=project(Z)
     multiphase_preconditioner(
      project_option,project_timings,
      presmooth,postsmooth,
      GMRES_BUFFER0_Z_MF+j_local,
      GMRES_BUFFER0_V_MF+j_local,nsolve);

      // W=A Z
      // 1. (begin)calls project_right_hand_side(Z)
      // 2. (end)  calls project_right_hand_side(W)
     applyALL(project_option,
       GMRES_BUFFER0_Z_MF+j_local,
       GMRES_BUFFER_W_MF,nsolve);

      // H_j is a j+2 x j+1 matrix  j=0..m-1
     for (int i=0;i<=j_local;i++) {
       // H_ij=W dot Vi
      dot_productALL(project_option,
       GMRES_BUFFER_W_MF,
       GMRES_BUFFER0_V_MF+i,HH[i][j_local],nsolve);
     } // i=0..j_local

      // G_j is a p(j)+1 x j+1 matrix  j=0..m-1
     for (int i=0;i<=p_local;i++) {
       // G_ij=W dot Ui
      dot_productALL(project_option,
       GMRES_BUFFER_W_MF,
       GMRES_BUFFER0_U_MF+i,GG[i][j_local],nsolve);
     } // i=0..p_local

     for (int i=0;i<=j_local;i++) {
      aa=-HH[i][j_local];
       // W=W+aa Vi
      mf_combine(project_option,
       GMRES_BUFFER_W_MF,GMRES_BUFFER0_V_MF+i,aa,GMRES_BUFFER_W_MF,nsolve); 
      project_right_hand_side(GMRES_BUFFER_W_MF,project_option,change_flag);
     }
     for (int i=0;i<=p_local;i++) {
      aa=-GG[i][j_local];
       // W=W+aa Ui
      mf_combine(project_option,
       GMRES_BUFFER_W_MF,GMRES_BUFFER0_U_MF+i,aa,GMRES_BUFFER_W_MF,nsolve); 
      project_right_hand_side(GMRES_BUFFER_W_MF,project_option,change_flag);
     }

      // H_j is a j+2 x j+1 matrix  j=0..m-1
     dot_productALL(project_option,
       GMRES_BUFFER_W_MF,
       GMRES_BUFFER_W_MF,HH[j_local+1][j_local],nsolve);

     if (HH[j_local+1][j_local]>=0.0) {
      HH[j_local+1][j_local]=sqrt(HH[j_local+1][j_local]);
     } else
      amrex::Error("HH[j_local+1][j_local] invalid");

     double local_tol=breakdown_tol;
     if (p_local>=0) {
      local_tol/=pow(10.0,2*p_local);
     } else if (p_local==-1) {
      // do nothing
     } else
      amrex::Error("p_local invalid");

     int condition_number_blowup=0;
     double zeyu_condnum=1.0;

     if (HH[j_local+1][j_local]>KRYLOV_NORM_CUTOFF) {

       // Real** HH   i=0..m  j=0..m-1
       // active region: i=0..j_local+1  j=0..j_local
      if (disable_additional_basis==0) {
       zeyu_condnum=CondNum(HH,m+1,m,j_local+2,j_local+1,local_tol);
      } else if (disable_additional_basis==1) {
       zeyu_condnum=0.0;
      } else
       amrex::Error("disable_additional_basis invalid");

      if (zeyu_condnum>1.0/local_tol) { 
       condition_number_blowup=1;  
      } else if (zeyu_condnum<=1.0/local_tol) {
       condition_number_blowup=0;  
      } else
       amrex::Error("zeyu_condnum NaN");

     } else if ((HH[j_local+1][j_local]>=0.0)&&
		(HH[j_local+1][j_local]<=KRYLOV_NORM_CUTOFF)) {
      condition_number_blowup=2;  
      use_previous_iterate=1;
     } else
      amrex::Error("HH[j_local+1][j_local] invalid");

     if (debug_BF_GMRES==1) {
      std::cout << "BFGMRES: j_local="<< j_local << '\n';
      std::cout << "BFGMRES: zeyu_condnum="<< zeyu_condnum << '\n';
      std::cout << "BFGMRES: condition_number_blowup="<< 
	     condition_number_blowup << '\n';

      std::cout << "BFGMRES: HH[j_local+1][j_local]=" <<
	      HH[j_local+1][j_local] << '\n';
      std::cout << "BFGMRES: HH[j_local+1][j_local]=" <<
	      HH[j_local+1][j_local] << '\n';
     }

     if (condition_number_blowup==1) {

      p_local++;

       // variables initialized to 0.0  dir=-1
      if (p_local>=p_local_allocated) {
       allocate_array(0,nsolveMM,-1,GMRES_BUFFER0_U_MF+p_local);
       p_local_allocated++;
      }

       // ngrow,ncomp,idx_dest,idx_source
      copyALL(0,nsolveMM,GMRES_BUFFER0_U_MF+p_local,
	      GMRES_BUFFER0_V_MF+j_local);

      if (j_local==0) {
	// do nothing
        // G_j-1 is a p(j-1)+1 x j matrix  j=0..m-1
        // H_j-1 is a j+1 x j matrix  j=0..m-1
      } else if ((j_local>=1)&&(j_local<m)) {
       for (int i=0;i<j_local;i++) {
        GG[p_local][i]=HH[j_local][i];
	HH[j_local][i]=0.0;
       }
      } else
       amrex::Error("j_local invalid");

      int max_vhat_sweeps=4;
      int vhat_counter=0;
      for (vhat_counter=0;((vhat_counter<max_vhat_sweeps)&&
			   (condition_number_blowup==1));vhat_counter++) {

       setVal_array(0,nsolveMM,1.0,GMRES_BUFFER0_V_MF+j_local);
       init_checkerboardALL(GMRES_BUFFER0_V_MF+j_local,
		     project_option,nsolve,nsolveMM);
       project_right_hand_side(GMRES_BUFFER0_V_MF+j_local,
		     project_option,change_flag);
        // v_i dot v_j = 0 i=0..j-1
	// u_p dot v_j = 0
	// for i=0..j-1
	//  v_j = v_j - (vj,vi)vi/(vi,vi)

       double vj_dot_vi=0.0;
       double vi_dot_vi=0.0;
       for (int i=0;i<=j_local;i++) {
        int idx_vi;
        if ((i>=0)&&(i<j_local)) {
         idx_vi=GMRES_BUFFER0_V_MF+i;
        } else if (i==j_local) {
         idx_vi=GMRES_BUFFER0_U_MF+p_local;
        } else
         amrex::Error("i invalid");

	 // SANITY CHECK
        dot_productALL(project_option,
         GMRES_BUFFER0_V_MF+j_local,
         GMRES_BUFFER0_V_MF+j_local,vi_dot_vi,nsolve);

        if (vi_dot_vi>0.0) {
	 // do nothing
	} else {	
         std::cout << "vi_dot_vi= " << vi_dot_vi << endl;
	 std::cout << "i= " << i << endl;
         std::cout << "vhat_counter,j_local,p_local,m " <<
          vhat_counter << ' ' << ' ' << j_local << ' ' <<
          p_local << ' ' << m << endl;
         std::cout << "beta= " << beta << endl;
         for (int eh=0;eh<j_local;eh++) {
          std::cout << "eh,error_historyGMRES[eh] " << eh << ' ' <<
           error_historyGMRES[eh] << endl;
         }
         std::cout << "nsolve= " << nsolve << endl;
         std::cout << "project_option= " << project_option << '\n';

         amrex::Error("vi_dot_vi==0.0 (mid loop (main) sanity check)");
	}

        dot_productALL(project_option,
         idx_vi,
         GMRES_BUFFER0_V_MF+j_local,vj_dot_vi,nsolve);
        dot_productALL(project_option,
         idx_vi,
         idx_vi,vi_dot_vi,nsolve);
	if (vi_dot_vi>0.0) {
         aa=-vj_dot_vi/vi_dot_vi;
          // vj=vj+aa vi
         mf_combine(project_option,
           GMRES_BUFFER0_V_MF+j_local,idx_vi,aa,
	   GMRES_BUFFER0_V_MF+j_local,nsolve); 
         project_right_hand_side(GMRES_BUFFER0_V_MF+j_local,
			 project_option,change_flag);
	} else
 	 amrex::Error("vi_dot_vi==0.0");
       } // i=0..j_local
       dot_productALL(project_option,
        GMRES_BUFFER0_V_MF+j_local,
        GMRES_BUFFER0_V_MF+j_local,vi_dot_vi,nsolve);
       if (vi_dot_vi>0.0) {
        aa=1.0/vi_dot_vi;
        mult_array(0,nsolveMM,aa,GMRES_BUFFER0_V_MF+j_local);
        project_right_hand_side(GMRES_BUFFER0_V_MF+j_local, 
			project_option,change_flag);
       } else
        amrex::Error("vi_dot_vi==0.0 (renormalization)");

        // Zj=M^{-1}Vj
        // Z=project(Z)
       multiphase_preconditioner(
         project_option,project_timings,
         presmooth,postsmooth,
         GMRES_BUFFER0_Z_MF+j_local,
         GMRES_BUFFER0_V_MF+j_local,nsolve);

         // W=A Z
         // 1. (begin)calls project_right_hand_side(Z)
         // 2. (end)  calls project_right_hand_side(W)
       applyALL(project_option,
         GMRES_BUFFER0_Z_MF+j_local,
         GMRES_BUFFER_W_MF,nsolve);

        // H_j is a j+2 x j+1 matrix  j=0..m-1
       for (int i=0;i<=j_local;i++) {
        // H_ij=W dot Vi
        dot_productALL(project_option,
          GMRES_BUFFER_W_MF,
          GMRES_BUFFER0_V_MF+i,HH[i][j_local],nsolve);
       } // i=0..j_local

        // G_j is a p(j)+1 x j+1 matrix  j=0..m-1
       for (int i=0;i<=p_local;i++) {
        // G_ij=W dot Ui
        dot_productALL(project_option,
         GMRES_BUFFER_W_MF,
         GMRES_BUFFER0_U_MF+i,GG[i][j_local],nsolve);
       } // i=0..p_local

       for (int i=0;i<=j_local;i++) {
        aa=-HH[i][j_local];
         // W=W+aa Vi
        mf_combine(project_option,
         GMRES_BUFFER_W_MF,GMRES_BUFFER0_V_MF+i,aa,GMRES_BUFFER_W_MF,nsolve); 
        project_right_hand_side(GMRES_BUFFER_W_MF,project_option,change_flag);
       }
       for (int i=0;i<=p_local;i++) {
        aa=-GG[i][j_local];
         // W=W+aa Ui
        mf_combine(project_option,
         GMRES_BUFFER_W_MF,GMRES_BUFFER0_U_MF+i,aa,GMRES_BUFFER_W_MF,nsolve); 
        project_right_hand_side(GMRES_BUFFER_W_MF,project_option,change_flag);
       }

        // H_j is a j+2 x j+1 matrix  j=0..m-1
       dot_productALL(project_option,
         GMRES_BUFFER_W_MF,
         GMRES_BUFFER_W_MF,HH[j_local+1][j_local],nsolve);

       condition_number_blowup=0;
       zeyu_condnum=1.0;

       if (HH[j_local+1][j_local]>=0.0) {
        HH[j_local+1][j_local]=sqrt(HH[j_local+1][j_local]);
       } else
        amrex::Error("HH[j_local+1][j_local] invalid");

       if (HH[j_local+1][j_local]>KRYLOV_NORM_CUTOFF) {
         // Real** HH   i=0..m  j=0..m-1
         // active region: i=0..j_local+1  j=0..j_local
         
        zeyu_condnum=CondNum(HH,m+1,m,j_local+2,j_local+1,local_tol);
        if (zeyu_condnum>1.0/local_tol) { 
         condition_number_blowup=1;  
        } else if (zeyu_condnum<=1.0/local_tol) {
         condition_number_blowup=0;  
        } else
         amrex::Error("zeyu_condnum NaN");

       } else if ((HH[j_local+1][j_local]>=0.0)&&
  		  (HH[j_local+1][j_local]<=KRYLOV_NORM_CUTOFF)) {
        condition_number_blowup=2;  
        use_previous_iterate=1;
	p_local--;
       } else
        amrex::Error("HH[j_local+1][j_local] invalid");

       if (debug_BF_GMRES==1) {
        std::cout << "BFGMRES: j_local="<< j_local << '\n';
        std::cout << "BFGMRES: p_local="<< p_local << '\n';
        std::cout << "BFGMRES: vhat_counter="<< vhat_counter << '\n';
        std::cout << "BFGMRES: zeyu_condnum="<< zeyu_condnum << '\n';
        std::cout << "BFGMRES: condition_number_blowup="<< 
	     condition_number_blowup << '\n';

        std::cout << "BFGMRES: HH[j_local+1][j_local]=" <<
	      HH[j_local+1][j_local] << '\n';
        std::cout << "BFGMRES: HH[j_local+1][j_local]=" <<
	      HH[j_local+1][j_local] << '\n';
       }

      } // vhat_counter=0..max_vhat_sweeps-1 or condition_number_blowup==0

      if (vhat_counter==max_vhat_sweeps) {

       use_previous_iterate=1;

      } else if ((vhat_counter>0)&&(vhat_counter<max_vhat_sweeps)) {
       if (condition_number_blowup==0) {
        // do nothing
       } else if (condition_number_blowup==2) {
        use_previous_iterate=1;
       } else
        amrex::Error("condition_number_blowup invalid");
      } else
       amrex::Error("vhat_counter invalid");

     } else if (condition_number_blowup==0) {
      // do nothing
     } else if (condition_number_blowup==2) {
      use_previous_iterate=1;
     } else
      amrex::Error("condition_number_blowup invalid");

     if (use_previous_iterate==0) {

       // H_j is a j+2 x j+1 matrix j=0..m-1
       // G_j is a p(j)+1 x j+1 matrix j=0..m-1
      for (int i1=0;i1<=j_local+1;i1++) {
       for (int i2=0;i2<=j_local;i2++) {
        HHGG[i1][i2]=HH[i1][i2];
       }
      }
      for (int i1=0;i1<=p_local;i1++) {
       for (int i2=0;i2<=j_local;i2++) {
        HHGG[i1+j_local+2][i2]=GG[i1][i2];
       }
      }
      status=1;
      zeroALL(1,nsolveMM,idx_Z);

      // sub box dimensions: p_local+j_local+3  x j_local+1
      // j_local=0..m-1
      LeastSquaresQR(HHGG,yy,beta_e1,2*m+2,2*m+1,
         	     p_local+j_local+3,j_local+1);
      if (status==1) {
       for (int i2=0;i2<=j_local;i2++) {
        aa=yy[i2];
         // Z=Z+aa Z_{i2}
        mf_combine(project_option,
         idx_Z,
         GMRES_BUFFER0_Z_MF+i2,aa,
         idx_Z,nsolve); 
        project_right_hand_side(idx_Z,project_option,change_flag);
       }
      } else
       amrex::Error("status invalid");

      if (disable_additional_basis==0) {

       // residALL calls applyALL which calls 
       // project_right_hand_side(idx_Z).
       // At end end of residALL, 
       // project_right_hand_side(W) is called.
       residALL(project_option,
        idx_R, // rhs
        GMRES_BUFFER_W_MF, // resid
        idx_Z, // source
        nsolve);

       Real beta_compare=0.0;
       dot_productALL(project_option,
	GMRES_BUFFER_W_MF, 
	GMRES_BUFFER_W_MF,beta_compare,nsolve);
       if (beta_compare>=0.0) {
        beta_compare=sqrt(beta_compare);
        if (beta_compare<=GMRES_tol*beta) {
         convergence_flag=1;
        } else if (beta_compare>GMRES_tol*beta) {
         convergence_flag=0;
        } else
         amrex::Error("beta_compare invalid");
       } else
        amrex::Error("beta_compare invalid");

       if (debug_BF_GMRES==1) {
        std::cout << "BFGMRES: beta_compare=" <<
         beta_compare << " beta=" << beta << '\n';
       }

       error_historyGMRES[j_local]=beta_compare;

       if (j_local==0) {
        if (beta_compare>beta) {
         convergence_flag=1;
        } else if (beta_compare<=beta) {
         // do nothing
        } else
         amrex::Error("beta_compare or beta invalid");
       } else if (j_local>0) {
        if (beta_compare>error_historyGMRES[j_local-1]) {
         convergence_flag=1;
        } else if (beta_compare<=error_historyGMRES[j_local-1]) {
         // do nothing
        } else
         amrex::Error("beta_compare or error_historyGMRES[j_local-1] invalid");
       } else
        amrex::Error("j_local invalid");

      } else if (disable_additional_basis==1) {

       convergence_flag=0;

      } else
       amrex::Error("disable_additional_basis invalid");

      if (convergence_flag==0) {

       if (HH[j_local+1][j_local]>0.0) {
        if ((j_local>=0)&&(j_local<m-1)) {

         // variables initialized to 0.0
         allocate_array(0,nsolveMM,-1,GMRES_BUFFER0_V_MF+j_local+1);

         aa=1.0/HH[j_local+1][j_local];
          // V=V+aa W
         mf_combine(project_option,
          GMRES_BUFFER0_V_MF+j_local+1,
          GMRES_BUFFER_W_MF,aa,
          GMRES_BUFFER0_V_MF+j_local+1,nsolve); 
         project_right_hand_side(GMRES_BUFFER0_V_MF+j_local+1,
		       project_option,change_flag);
        } else if (j_local==m-1) {
         // do nothing
        } else
         amrex::Error("j_local invalid");
       } else if (HH[j_local+1][j_local]==0.0) {
        amrex::Error("HH[j_local+1][j_local] should not be 0");
       } else {
        amrex::Error("HH[j_local+1][j_local] invalid");
       }

      } else if (convergence_flag==1) {
       // do nothing
      } else
       amrex::Error("convergence_flag invalid");

     } else if (use_previous_iterate==1) {
      convergence_flag=1;
      if (j_local==0) {

        // Z=M^{-1}R
        // Z=project(Z)
       multiphase_preconditioner(
        project_option,project_timings,
        presmooth,postsmooth,
        idx_Z,idx_R,nsolve);
      }
     } else
      amrex::Error("use_previous_iterate invalid");

    } // j_local=0..m-1 (and convergence_flag==0)

    if (debug_BF_GMRES==1) {
     for (int i=0;i<j_local;i++) {
      std::cout << "main: i,error_historyGMRES[i] " << i << ' ' <<
	     error_historyGMRES[i] << endl;
     } 
    }

    for (int i=0;i<j_local;i++) {
     delete_array(GMRES_BUFFER0_Z_MF+i); 
    }
      // GMRES_BUFFER0_V_MF+0 deleted after if statement below. 
    for (int i=1;i<j_local;i++) {
     delete_array(GMRES_BUFFER0_V_MF+i); 
    }
    for (int i=p_local_init_allocated;i<p_local_allocated;i++) {
     delete_array(GMRES_BUFFER0_U_MF+i); 
    }

    for (int i=0;i<m+1;i++) 
     delete [] GG[i];
    delete [] GG;

    for (int i=0;i<2*m+2;i++) 
     delete [] HHGG[i];
    delete [] HHGG;

   } else
    amrex::Error("breakdown_free_flag invalid");

   delete [] yy;
   delete [] beta_e1;

   for (int i=0;i<m+1;i++) 
    delete [] HH[i];
   delete [] HH;

   delete_array(GMRES_BUFFER0_V_MF); 
   delete_array(GMRES_BUFFER_W_MF);

  } else if (beta==0.0) {
    // Z=M^{-1}R
    // Z=project(Z)
   multiphase_preconditioner(
    project_option,project_timings,
    presmooth,postsmooth,
    idx_Z,idx_R,nsolve);
  } else
   amrex::Error("beta invalid");

 } else {
  std::cout << "gmres_precond_iter= " << gmres_precond_iter << '\n';
  amrex::Error("NavierStokes3.cpp: gmres_precond_iter invalid");
 }
 project_right_hand_side(idx_Z,project_option,change_flag);

} // end subroutine multiphase_GMRES_preconditioner

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
 
 int nmat=num_materials;

 int num_materials_face=num_materials_vel;
 if (project_option_momeqn(project_option)==1) {
  if (num_materials_face!=1)
   amrex::Error("num_materials_face invalid");
 } else if (project_option_momeqn(project_option)==0) {
  num_materials_face=num_materials_scalar_solve;
 } else
  amrex::Error("project_option_momeqn invalid42");

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 if ((num_materials_face!=1)&&
     (num_materials_face!=nmat))
  amrex::Error("num_materials_face invalid");

 int nsolveMM=nsolve*num_materials_face;

 zeroALL(1,nsolveMM,idx_Z);

#if (profile_solver==1)
 bprof.stop();
#endif

  // PCG
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
  int bicgstab_num_outer_iterSOLVER=0;

  int call_adjust_tolerance=0;

  jacobi_cycles(
   call_adjust_tolerance, //=0
   smooth_cycles,
   update_vel,  // =0 (i.e. jacobi_cycles used as a preconditioner)
   project_option,
   idx_R,idx_Z,
   error_at_the_beginning,
   error_after_all_jacobi_sweeps,
   error_n,    // not modified
   error0_max, // not modified
   bicgstab_num_outer_iterSOLVER,nsolve);

#if (profile_solver==1)
  bprof.stop();
#endif

   // MGPCG
 } else if (project_solver_type==0) {

  NavierStokes& ns_finest=getLevel(finest_level);
  ns_finest.mg_cycleALL(presmooth,
  project_option,
  idx_R,idx_Z,nsolve);

   // MINV=I
 } else if (project_solver_type==2) {

#if (profile_solver==1)
  bprof.start();
#endif

  for (int ilev=level;ilev<=finest_level;ilev++) {
   NavierStokes& ns_level=getLevel(ilev);
   int ncomp=ns_level.localMF[idx_Z]->nComp();
   int ngrow=ns_level.localMF[idx_Z]->nGrow();
   ns_level.setVal_localMF(idx_Z,0.0,0,ncomp,ngrow);
   MultiFab::Copy(*ns_level.localMF[idx_Z],
		  *ns_level.localMF[idx_R],0,0,ncomp,0);
  } // ilev=level ... finest_level

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

 save_atol_b=1.0e-14;
 save_mac_abs_tol=mac_abs_tol;
 save_min_rel_error=minimum_relative_error;
 ParmParse pp("mg");

 if (project_option_singular_possible(project_option)==1) {
  save_mac_abs_tol=mac_abs_tol;
  save_atol_b=0.01*save_mac_abs_tol;
  pp.query("bot_atol",save_atol_b);
  save_min_rel_error=minimum_relative_error;
 } else if (project_option==2) {
  save_mac_abs_tol=thermal_abs_tol;
  save_atol_b=0.01*save_mac_abs_tol; 
  pp.query("thermal_bot_atol",save_atol_b);
  save_min_rel_error=diffusion_minimum_relative_error;
 } else if (project_option==3) {
  save_mac_abs_tol=visc_abs_tol;
  save_atol_b=0.01*save_mac_abs_tol; 
  pp.query("visc_bot_atol",save_atol_b);
  save_min_rel_error=diffusion_minimum_relative_error;
 } else if ((project_option>=100)&&
	    (project_option<100+num_species_var)) {//species
  save_mac_abs_tol=visc_abs_tol;
  save_atol_b=0.01*save_mac_abs_tol; 
  pp.query("visc_bot_atol",save_atol_b);
  save_min_rel_error=diffusion_minimum_relative_error;
 } else if (project_option==200) { //smoothing
  save_mac_abs_tol=visc_abs_tol;
  save_atol_b=0.01*save_mac_abs_tol; 
  pp.query("visc_bot_atol",save_atol_b);
  save_min_rel_error=diffusion_minimum_relative_error;
 } else
  amrex::Error("project_option invalid 51");

} // end subroutine set_local_tolerances

// project_option=0  regular project
// project_option=1  initial project
// project_option=2  temperature solve
// project_option=3  viscous forces.
// project_option=11 FSI_material_exists (last project)
// project_option=12 extend pressure into cells where all coefficients==0.0
// project_option=100,..,100+num_species_var-1  species vars.
// project_option=200 smoothing
// 
// if project_option=0,
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

 int finest_level=parent->finestLevel();
 if (level!=0)
  amrex::Error("level invalid multiphase_project");

 int nmat=num_materials;

 std::fflush(NULL);

#if (profile_solver==1)
 std::string subname="NavierStokes::multiphase_project";
 std::stringstream popt_string_stream(std::stringstream::in |
      std::stringstream::out);
 popt_string_stream << project_option;
 std::string profname=subname+popt_string_stream.str();
 BLProfiler bprof(profname);
#endif


 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid multiphase_project");

 const Real* coarse_dx=geom.CellSize();

 
  // FSI_material_exists last project
  // The independent variable is "DIV_Type"
  // which means the contents must be saved.
 if (project_option==11) { // FSI_material_exists (last project)
  allocate_array(1,1,-1,DIV_SAVE_MF);
  for (int ilev=level;ilev<=finest_level;ilev++) {
   NavierStokes& ns_level=getLevel(ilev);
   MultiFab& DIV_new=ns_level.get_new_data(DIV_Type,slab_step+1);
   MultiFab::Copy(
      *ns_level.localMF[DIV_SAVE_MF],
      DIV_new,0,0,1,1);
  } // ilev=level ... finest_level
  // in: MacProj.cpp
  // if compressible: DIV_new=-dt(pnew-padv)/(rho c^2 dt^2)+MDOT_MF dt/vol
  // if incompressible: DIV_new=MDOT_MF dt/vol
  ADVECT_DIV_ALL();
 } else if (project_option==12) { // pressure extension
  allocate_array(1,1,-1,PRESSURE_SAVE_MF);
  for (int ilev=level;ilev<=finest_level;ilev++) {
   NavierStokes& ns_level=getLevel(ilev);
   MultiFab& P_new=ns_level.get_new_data(State_Type,slab_step+1);
   int pcomp=num_materials_vel*AMREX_SPACEDIM;
   MultiFab::Copy(
      *ns_level.localMF[PRESSURE_SAVE_MF],
      P_new,pcomp,0,1,1);
  } // ilev=level ... finest_level
 } else if (project_option==200) { // smoothing

  if ((localMF_grow[TEMPERATURE_SAVE_MF]==-1)&&
      (localMF_grow[TEMPERATURE_SMOOTH_MF]==-1)) {
   // do nothing
  } else
   amrex::Error("TEMPERATURE_SAVE_MF or TEMPERATURE_SMOOTH_MF not deleted");

  Vector<int> scomp_thermal;
  Vector<int> ncomp_thermal;
  int state_index_thermal;  
  int ncomp_check_thermal;
  get_mm_scomp_solver(
    nmat,
    project_option,
    state_index_thermal,
    scomp_thermal,
    ncomp_thermal,
    ncomp_check_thermal);
  if (ncomp_check_thermal!=nmat)
   amrex::Error("ncomp_check_thermal invalid");

   // data at time = cur_time_slab
  getState_localMF_listALL(
    TEMPERATURE_SAVE_MF,1,
    state_index_thermal,
    scomp_thermal,
    ncomp_thermal);

  int project_option_combine=2; 
  int combine_flag=0;  // FVM -> GFM 
   // combine_idx==-1 => update S_new  
   // combine_idx>=0  => update localMF[combine_idx]
  int combine_idx=-1; 
  int update_flux=0;
  int hflag_combine=0;
  int interface_cond_avail=0;
 
  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   ns_level.combine_state_variable(
    project_option_combine,
    combine_idx,
    combine_flag,
    hflag_combine,
    update_flux,
    interface_cond_avail); 
  } // ilev=finest_level ... level

  for (int ilist=0;ilist<scomp_thermal.size();ilist++) 
   avgDownALL(state_index_thermal,scomp_thermal[ilist],
	ncomp_thermal[ilist],1);

 } else if (project_option_is_valid(project_option)==1) {
  // do not save anything
 } else
  amrex::Error("project_option invalid");

 int save_enable_spectral=enable_spectral;

 if (project_option_projection(project_option)==1) {
  override_enable_spectral(projection_enable_spectral);
 } else if (project_option==12) { // pressure extension
  override_enable_spectral(0); // always low order
 } else if (project_option==2) { // thermal diffusion
  // do nothing
 } else if (project_option==3) { // viscosity
  // do nothing
 } else if ((project_option>=100)&&
	    (project_option<100+num_species_var)) { //species
  // do nothing
 } else if (project_option==200) { //smoothing
  override_enable_spectral(0); // always low order
 } else
  amrex::Error("project_option invalid43");

 int energyflag=0;
 int scomp_den=num_materials_vel*(AMREX_SPACEDIM+1);

 int num_materials_face=num_materials_vel;
 if (project_option_momeqn(project_option)==1) {
  if (num_materials_face!=1)
   amrex::Error("num_materials_face invalid");
 } else if (project_option_momeqn(project_option)==0) {
  num_materials_face=num_materials_scalar_solve;
 } else
  amrex::Error("project_option invalid44");

 if ((num_materials_face!=1)&&
     (num_materials_face!=nmat))
  amrex::Error("num_materials_face invalid");

 int nsolve=1;

 if (project_option_projection(project_option)==1) {

  singular_possible=1; // all zero coefficients possible in solid regions.
  local_solvability_projection=solvability_projection;
  if (project_option==1) { // initial project
   // do nothing
  } else if ((project_option==0)||   //regular project
             (project_option==11)) { //FSI_material_exists last project
   if (some_materials_compressible())
    local_solvability_projection=0;
  } else
   amrex::Error("project_option invalid 45"); 
 } else if (project_option==12) { // pressure extension
  singular_possible=1; //all zero coefficients possible in non-solid regions.
  local_solvability_projection=0;
 } else if (project_option==2) { // thermal conduction
  singular_possible=0; // diagonally dominant everywhere.
  local_solvability_projection=0;
 } else if (project_option==3) { // viscosity
  singular_possible=0; // diagonally dominant everywhere.
  nsolve=AMREX_SPACEDIM;
  local_solvability_projection=0;
 } else if ((project_option>=100)&&
            (project_option<100+num_species_var)) {//species
  singular_possible=0; // diagonally dominant everywhere.
  local_solvability_projection=0;
 } else if (project_option==200) { //smoothing
  singular_possible=0; // diagonally dominant everywhere.
  local_solvability_projection=0;
 } else
  amrex::Error("project_option invalid multiphase_project");

 Vector<int> scomp;
 Vector<int> ncomp;
 int state_index;  
 int ncomp_check;
 get_mm_scomp_solver(
  num_materials_face,
  project_option,
  state_index,
  scomp,
  ncomp,
  ncomp_check);
 int nsolveMM=nsolve*num_materials_face;
 if (ncomp_check!=nsolveMM)
  amrex::Error("nsolveMM invalid 4904");

 int nsolveMM_FACE=nsolveMM;
 if (num_materials_face==1) {
  // do nothing
 } else if (num_materials_face==nmat) {
  nsolveMM_FACE*=2;
 } else
  amrex::Error("num_materials_face invalid");

  // localMF[UMAC_MF] = 0.0
 allocate_MAC_velocityALL(nsolve,UMAC_MF);

 resize_metrics(1);
 debug_ngrow(VOLUME_MF,1,700);
 debug_ngrow(FACE_VAR_MF,0,702);

 int project_timings=0;

 double begin_project=0.0;
 if (project_timings==1)
  begin_project=ParallelDescriptor::second();


  // if project_option==1 (initial project), then 
  // initial guess is p=0 and homogeneous BC are used
  // for p (the potential function).
  //
  // if project_option==0, then:
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
    // diffusionRHS=mdot project_option==0
    // diffusionRHS=mdot pressure=0 project_option==1
    // diffusionRHS=0.0 UMAC=0 project_option==2 (thermal conduction)
    // diffusionRHS=0.0 UMAC=0 project_option==3 (viscosity)
    // diffusionRHS=0.0 UMAC=0 project_option==100.. (species)
    // diffusionRHS=0.0 if project_option=11 (FSI_material_exists last project)
    // diffusionRHS=0.0 if project_option=12 (pressure_extension)
    // (DIV_Type contents cannot be zapped because it is needed for
    //  initializing CELL_SOUND_MF, DIFFUSIONRHS if project_option==11)
   ns_level.Prepare_UMAC_for_solver(project_option,nsolve);
 }  // ilev=level ... finest_level

 if (project_option==11) {
   check_value_max(1,DIFFUSIONRHS_MF,0,1,0,0.0);
 }

 Real max_nlevels=0;
 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  int nlevels=ns_level.NSnumLevels();
  if (nlevels>max_nlevels)
   max_nlevels=nlevels;
 } // ilev=level ... finest_level

 std::fflush(NULL);

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << " ---------------------------------------------- \n";
   std::cout << " ELLIPTIC SOLVE project_option= " <<
    project_option << '\n';
   std::cout << " num_materials_vel= " << num_materials_vel << '\n';
   std::cout << " num_materials_face= " << num_materials_face << '\n';
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
  // automatically initializes mac_rhs_crse_array=0.0
 allocate_rhs_var(nsolve,MAC_RHS_CRSE_MF);
 
 if (project_option==11) {
   check_value_max(2,DIFFUSIONRHS_MF,0,1,0,0.0);
 }

 min_face_wt.resize(thread_class::nthreads);
 max_face_wt.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  min_face_wt[tid].resize(4);
  max_face_wt[tid].resize(4);
  for (int iwt=0;iwt<4;iwt++) {
   min_face_wt[tid][iwt]=1.0e+20;
   max_face_wt[tid][iwt]=-1.0e+20;
  }
 } // tid

  // in multiphase_project
 if (project_option==0) {

   // gravity and surface tension
  process_potential_forceALL();


// 1. overwrites cell/face velocity perhaps
// 2. must be called before adding gravity and surface tension.
// 3. cannot be called after the project because the velocity
//    will then fail to satisfy the discrete divergence condition.
  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   ns_level.overwrite_outflow();  
  }

    // gravity and surface tension cell/face
    // FUTURE: E+=dt u dot g + dt^2 g dot g/2
  increment_potential_forceALL(); 

  if (1==0) {
   int basestep_debug=nStep()+1;
   parent->writeDEBUG_PlotFile(basestep_debug,SDC_outer_sweeps,slab_step);
   std::cout << "press any number then enter AFTER GRAV\n";
   int n_input;
   std::cin >> n_input;
  }  

  if (visual_buoyancy_plot_int>0) {

   int nsteps=parent->levelSteps(0); // nsteps==0 very first step.

   if (very_last_sweep==1) {
    int ratio=(nsteps+1)/visual_buoyancy_plot_int;
    ratio=ratio*visual_buoyancy_plot_int;
    if (ratio==nsteps+1) {

     int caller_id=1;

     for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
      // filenames: "FACE_VAR<stuff>.plt" (MAC data)
      // curv_index,pforce_index (unused for now),
      // faceden_index,facecut_index,
      // icefacecut_index=4,icemask_index,facevisc_index,
      // faceheat_index,facevel_index,facespecies_index,
      // smoothing_index,
      // massface_index,vofface_index
      writeSanityCheckData(
       "FACE_VAR",
       "project_option==0:FACE_VAR_MF",//faceden_index=2 facevisc_index=6
       caller_id,
       localMF[FACE_VAR_MF+dir]->nComp(),
       FACE_VAR_MF+dir,
       -1, // State_Type==-1
       dir);

       // gravity * dt
      writeSanityCheckData(
       "POTENTIAL_FORCE_EDGE",
       "project_option==0:POTENTIAL_FORCE_EDGE",
       caller_id,
       localMF[POTENTIAL_FORCE_EDGE_MF+dir]->nComp(),
       POTENTIAL_FORCE_EDGE_MF+dir,
       -1, // State_Type==-1
       dir);
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


  deallocate_potential_forceALL(); 


   // div up and grad p  cell/face
   // T=T-(1/(rho cv))(int div(up)-dt div(up))
   // u=u-(1/rho)(int gp - dt gp)
  if ((SDC_outer_sweeps>0)&&
      (SDC_outer_sweeps<ns_time_order)&&
      (divu_outer_sweeps+1==num_divu_outer_sweeps)) {
   for (int ilev=finest_level;ilev>=level;ilev--) {
    NavierStokes& ns_level=getLevel(ilev);
    int local_project_option=0;
    ns_level.make_SEM_delta_force(local_project_option); 
   }
  } else if (SDC_outer_sweeps==0) {
   // do nothing
  } else if (divu_outer_sweeps+1<num_divu_outer_sweeps) {
   // do nothing
  } else
   amrex::Error("SDC_outer_sweeps or divu_outer_sweeps invalid multiphase prj");

 }  // project_option==0 


 if (project_option==11) { // FSI_material_exists, last project
   check_value_max(3,DIFFUSIONRHS_MF,0,1,0,0.0);
 }

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

 if (project_option==11) { //FSI_material_exists, last project
   check_value_max(4,DIFFUSIONRHS_MF,0,1,0,0.0);
 }

 if (project_option_projection(project_option)==1) {

  Vector<blobclass> blobdata;
  Vector< Vector<Real> > mdot_data;
  Vector< Vector<Real> > mdot_comp_data;
  Vector< Vector<Real> > mdot_data_redistribute;
  Vector< Vector<Real> > mdot_comp_data_redistribute;
  Vector<int> type_flag;

  int alloc_blobdata=0;
  
  if ((project_option==11)|| //FSI_material_exists last project
      ((project_option==0)&&
       (FSI_material_exists()==1))||
      (project_option==1)) {
   alloc_blobdata=1;
  }

  if (alloc_blobdata==1) {

   if (verbose>0) {
    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "BEGIN: color_variable, multiphase_project " <<
	          "project_option=" << project_option << '\n';
    }
   }

   int color_count=0;
   int coarsest_level=0;

   if (project_option==11) {
    check_value_max(41,DIFFUSIONRHS_MF,0,1,0,0.0);
   }

   int idx_mdot=-1;

   int tessellate=1;
   int operation_flag=0;
   ColorSumALL(
     operation_flag, // =0
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

   if (project_option==11) {
    check_value_max(42,DIFFUSIONRHS_MF,0,1,0,0.0);
   }

   if (verbose>0) {
    if (ParallelDescriptor::IOProcessor()) {
     std::cout << "END: color_variable, multiphase_project " <<
	          "project_option= " << project_option << '\n';
    }
   }

  } else if (alloc_blobdata==0) {
   // do nothing
  } else
   amrex::Error("alloc_blobdata invalid");

  int idx_velcell=-1;

  if (num_materials_vel==1) {
   // do nothing
  } else
   amrex::Error("num_materials_vel invalid");

  int interp_option=0;
  Real beta=0.0;

  if (nsolve!=1)
   amrex::Error("nsolve invalid");

  if ((project_option==0)||
      (project_option==11)) { //FSI_material_exists last project
   // unew^{f} = unew^{f} 
   interp_option=1;
  } else if (project_option==1) {
   //unew^{f} = unew^{c->f}
   interp_option=0;
  } else 
   amrex::Error("project_option invalid47");

   // if project_option==11 
   // then the velocity in the ice
   // is overwritten with a projected rigid body velocity:
   //  a) call get_rigid_velocity
   //  b) uedge(im_vel)=test_current_icefacecut*uedge(im_vel)+ &
   //       (one-test_current_icefacecut)*uedge_rigid
   //
   // In LEVELSET_3D.F90, FORT_CELL_TO_MAC,
   // operation_flag=3,4,5,10,11, num_colors.gt.0, if either
   // left cell or right cell has ice (FSI_flag==3,6)
   // or is FSI_rigid (FSI_flag==5) then velocity is
   // overwritten.
   // In GODUNOV_3D.F90, FORT_INIT_ICEMASK,
   // "icefacecut_index" is initialized by calling get_icemask
   // and if "im_FSI_rigid==im_primary" for one of a faces'
   // adjoining cells for example, then, icefacecut_index=0.0
   //

  increment_face_velocityALL(
    interp_option,project_option,
    idx_velcell,beta,blobdata); 

  if (project_option==11) {
   check_value_max(43,DIFFUSIONRHS_MF,0,1,0,0.0);
  }

  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);

     // 1. MAC_TEMP=Umac_new, UMAC=Umac_new
     // 2. call to residual_correction_form (UMAC=UMAC-beta grad p_init)
     // 3. MAC_TEMP=UMAC
     // 4. .... (UMAC updated)
   for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
    MultiFab* macvel=
      ns_level.getStateMAC(Umac_Type,0,dir,0,nsolveMM_FACE,cur_time_slab); 
    MultiFab::Copy(
      *ns_level.localMF[MAC_TEMP_MF+dir],
      *macvel,0,0,nsolveMM_FACE,0);
    MultiFab::Copy(
      *ns_level.localMF[UMAC_MF+dir],
      *macvel,0,0,nsolveMM_FACE,0);

    if (1==0) {
     int gridno=0;
     const Box& fabgrid = grids[gridno];
     const int* fablo=fabgrid.loVect();
     const int* fabhi=fabgrid.hiVect();
     const Real* xlo = grid_loc[gridno].lo();
     int interior_only=1;
     FArrayBox& macfab=(*macvel)[0];
     int scomp_debug=0;
     int ncomp_debug=nsolveMM_FACE;
     std::cout << "WARNING: this velocity is scaled\n";
     tecplot_debug(macfab,xlo,fablo,fabhi,coarse_dx,dir,0,
	     scomp_debug,ncomp_debug,interior_only);
    }
    delete macvel;
   }  // dir=0..sdim-1
  } // ilev=finest_level ... level

  if (project_option==11) {
   check_value_max(44,DIFFUSIONRHS_MF,0,1,0,0.0);
  }

  if (alloc_blobdata==1) {
   delete_array(TYPE_MF);
   delete_array(COLOR_MF);
  } else if (alloc_blobdata==0) {
   // do nothing
  } else
   amrex::Error("alloc_blobdata invalid");

 } else if (project_option_projection(project_option)==0) {
  // do nothing
 } else {
  amrex::Error("project_option_projection invalid48");
 } 

 if (project_option==11) {
  check_value_max(5,DIFFUSIONRHS_MF,0,1,0,0.0);
 }

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);

  if (project_option==1) { // initial project
   // do nothing
  } else if (project_option_projection(project_option)==1) {

   // updates CELL_SOUND_MF, DIFFUSIONRHS, and 
   // S_new 
   //  State_Type if project_option==0 
   //  DIV_Type if project_option==11
   //   DIV_Type=-dt*(pnew-padv)/(rho c^2 dt^2) + mdot * dt/vol
   //   if (project_option==11) then
   //    if (incomp): csnd(2)=0, DIFFSUIONRHS=(1/dt)( DIV_Type ) * vol
   //    if (comp): csnd(2)=DIV_TYPE
   // 
   //  NavierStokes::init_advective_pressure declared in NavierStokes2.cpp
   ns_level.init_advective_pressure(project_option); 
  } else if (project_option_projection(project_option)==0) {
   // do nothing 
  } else
   amrex::Error("project_option_projection invalid 49");

   // in: multiphase_project
   // calls FORT_BUILDFACEWT
   // BUILDFACEWT updates static variables min_face_wt and max_face_wt
   // max_face_wt[0][1] has max of (1/rho) or (visc_coef*mu) or (k) or (D)
  ns_level.allocate_FACE_WEIGHT(nsolve,project_option);

  ns_level.allocate_pressure_work_vars(nsolve,project_option);

   // updates the following variables:
   // ONES_MF, DOTMASK_MF, 
   // POLDHOLD_MF=S^adv - S^init, 
   // OUTER_ITER_PRESSURE_MF=S^init,
   // snew=S^init
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
 for (int im=1;im<nmat;im++) {
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

  } else
   amrex::Error("problen_max invalid");
 } else
  amrex::Error("maxden invalid");

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

 energyflag=0;
 int homflag_residual_correction_form=0; 

 if ((project_option==11)|| // FSI_material_exists last project
     (project_option==1)) { // initial project
  homflag_residual_correction_form=1; 
 } else if (project_option_is_valid(project_option)==1) {
  homflag_residual_correction_form=0; 
 } else
  amrex::Error("project_option invalid");

 CPP_OVERRIDEPBC(homflag_residual_correction_form,project_option);

   // STATE_FOR_RESID is an input to 
   //  NavierStokes::residual_correction_form
   
  // ngrow=1
 getState_localMF_listALL(STATE_FOR_RESID_MF,1,
   state_index,scomp,ncomp);

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  if (ns_level.localMF[STATE_FOR_RESID_MF]->nComp()!=nsolveMM)
   amrex::Error("ns_level.localMF[STATE_FOR_RESID_MF]->nComp()!=nsolveMM");

  ns_level.resize_levelsetLO(2,LEVELPC_MF);
  ns_level.debug_ngrow(LEVELPC_MF,2,870);
 } // ilev=finest_level ... level

    // initializes diagsing,mask_div_residual,mask_residual,ONES_MF,
    // ones_sum_global
    //
    //  i.e.
    //  
    // calls:FORT_SCALARCOEFF,FORT_MULT_FACEWT,FORT_DIVIDEDX,FORT_NSGENERATE
    // initializes arrays holding the diagonal and ONES_MF.
 int create_hierarchy=0;
 allocate_maccoefALL(project_option,nsolve,create_hierarchy);

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
  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   ns_level.putState_localMF_list(STATE_FOR_RESID_MF,
     state_index,scomp,ncomp);
  } // ilev=finest_level ... level
  delete_array(STATE_FOR_RESID_MF);
  getState_localMF_listALL(STATE_FOR_RESID_MF,1,
    state_index,scomp,ncomp);
 } else
  amrex::Error("change_flag invalid");

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);

     // UMAC_MF-=beta grad STATE_FOR_RESID
  ns_level.residual_correction_form(
   homflag_residual_correction_form,
   energyflag,
   project_option,nsolve);

  if (ilev<finest_level)
   ns_level.avgDownMac();  // interpolates UMAC_MF from ilev+1

    // mac_temp store the velocity that comes from residual_correction_form
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   ns_level.Copy_localMF(MAC_TEMP_MF+dir,UMAC_MF+dir,0,0,nsolveMM_FACE,0);
  }

 }  // ilev=finest_level ... level

 delete_array(STATE_FOR_RESID_MF);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  show_norm2_id(MAC_TEMP_MF+dir,5+dir);
 }

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "iwt:0=denface 1=cutface 2=desing. 3=sing";
   std::cout << "project_option= " << project_option << '\n';
   for (int iwt=0;iwt<4;iwt++) {
    std::cout << "iwt= " << iwt << " min " << 
      min_face_wt[0][iwt] << " max " <<
      max_face_wt[0][iwt] << '\n';
   }
  }
 } // verbose>0

 deallocate_maccoefALL(project_option);
   
 int meets_tol=0;

 set_local_tolerances(project_option);

 Real error0_max=0.0;
 Real error0=0.0;
 Real error_n=0.0;
 Real rAr_error0=0.0;
 Real rAr_error_n=0.0;
 Real Ar_error0=0.0;
 Real Ar_error_n=0.0;
 Real error_at_the_beginning=0.0;
 Real error_after_all_jacobi_sweeps=0.0;
    
 if (verbose>0) {
   if (ParallelDescriptor::IOProcessor()) {
    std::cout << "LEVEL PROJECT from level=" << level << " to level=" <<
      ' ' << finest_level << '\n';
    std::cout << "PROJECT OPTION=" << project_option << '\n';
    std::cout << "NSOLVE=" << nsolve << '\n';
   }
 }

 int vcycle;
 int bicgstab_num_outer_iterSOLVER=0;
 int outer_iter_done=0;

 int min_bicgstab_outer_iter=0;

     // initializes diagsing,mask_div_residual,mask_residual,ONES_MF,
     //  ones_sum_global
 create_hierarchy=1;
 allocate_maccoefALL(project_option,nsolve,create_hierarchy);

   // this must be done after allocate_maccoef (stefan_solver_init relies on
   // inhomogeneous BCs)
   // set BCs to homogeneous for the outer_iter loop.
 CPP_OVERRIDEPBC(1,project_option);

 int total_number_vcycles=0;

 bicgstab_num_outer_iterSOLVER=0;
 outer_iter_done=0;
 
 lev0_cycles_list.resize(0);

#if (profile_solver==1)
 bprof.stop();
#endif

 Vector< Array<Real,4> > outer_error_history;
 outer_error_history.resize(bicgstab_max_num_outer_iter+1);
 for (int ehist=0;ehist<outer_error_history.size();ehist++) {
   outer_error_history[ehist][0]=0.0;
   outer_error_history[ehist][1]=0.0;
   outer_error_history[ehist][2]=0.0;
   outer_error_history[ehist][3]=0.0;
 }
 Real outer_error=0.0;

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
      ns_level.setVal_localMF(MAC_RHS_CRSE_MF,0.0,0,nsolveMM,0); 
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
    } else if (project_option==2) {//thermal diffusion
     jacobi_cycles_count=initial_thermal_cycles;
    } else if ((project_option>=100)&&
               (project_option<100+num_species_var)) { //species
     jacobi_cycles_count=initial_thermal_cycles;
    } else if (project_option==200) { //smoothing
     jacobi_cycles_count=initial_viscosity_cycles;
    } else if (project_option==3) { //viscosity
     jacobi_cycles_count=initial_viscosity_cycles;
    } else
     amrex::Error("project_option invalid52");

    int update_vel=1; // update error0 IF bicgstab_num_outer_iterSOLVER==0
    int call_adjust_tolerance=1;

      // NavierStokes::jacobi_cycles in NavierStokes3.cpp
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
      bicgstab_num_outer_iterSOLVER,
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

    Vector< Array<Real,4> > error_history;

#if (profile_solver==1)
    bprof.stop();
#endif

    for (int cg_loop=0;cg_loop<cg_loop_max;cg_loop++) {

#if (profile_solver==1)
     bprof.start();
#endif

     allocate_array(0,nsolveMM,-1,CGRESID_MF);
     allocate_array(1,nsolveMM,-1,P_SOLN_MF);
     allocate_array(0,nsolveMM,-1,bicg_V1_MF);

     for (int ilev=finest_level;ilev>=level;ilev--) {
      NavierStokes& ns_level=getLevel(ilev);
      if (ilev<finest_level) {
            // get rid of uninit.
       ns_level.setVal_localMF(MAC_RHS_CRSE_MF,0.0,0,nsolveMM,0);
       ns_level.averageRhs(MAC_RHS_CRSE_MF,nsolve,project_option);
       ns_level.avgDownMac(); // works on UMAC_MF
      }
       // mac_phi_crse_mf=0.0
       // mac_rhs_crse=POLDHOLD * alpha - 
       //              vol div UMAC/dt + diffusionRHS
      ns_level.mac_project_rhs(project_option,MAC_PHI_CRSE_MF,
        MAC_RHS_CRSE_MF,nsolve);
     } // ilev=finest_level ... level

       // MAC_PHI_CRSE=0.0 (from above)
       // CGRESID=MAC_RHS_CRSE-( alpha*phi-div grad phi )
       // if local_solvability_projection, then this
       // routine modifies CGRESID so that it sums to zero.
       // if singular_possible, then this routine zeros out the
       // residual where the matrix diagonal (prior to dual time stepping
       // modification) is 0.
     residALL(project_option,MAC_RHS_CRSE_MF,
      CGRESID_MF,MAC_PHI_CRSE_MF,nsolve);

     Real local_error_n=0.0;
     dot_productALL(project_option,CGRESID_MF,CGRESID_MF,
      	      local_error_n,nsolve);
     if (local_error_n>=0.0) {
      local_error_n=sqrt(local_error_n);
     } else
      amrex::Error("local_error_n invalid");

     copyALL(0,nsolveMM,P_SOLN_MF,CGRESID_MF); //P_SOLN=CGRESID_MF
       // V1=A P_SOLN
       // 1. (begin)calls project_right_hand_side(P)
       // 2. (end)  calls project_right_hand_side(V1)
     applyALL(project_option,P_SOLN_MF,bicg_V1_MF,nsolve);
     dot_productALL(project_option,CGRESID_MF,bicg_V1_MF,rAr_error_n,nsolve);
     dot_productALL(project_option,bicg_V1_MF,bicg_V1_MF,Ar_error_n,nsolve);

     if (rAr_error_n<0.0) {
      rAr_error_n=0.0;
     } else if (rAr_error_n>=0.0) {
      rAr_error_n=sqrt(rAr_error_n);
     } else
      amrex::Error("rAr_error_n invalid");

     if (Ar_error_n>=0.0) {
      Ar_error_n=sqrt(Ar_error_n);
     } else
      amrex::Error("Ar_error_n invalid");

     if ((bicgstab_num_outer_iterSOLVER==0)&&(cg_loop==0)) {
      Ar_error0=Ar_error_n;
      rAr_error0=rAr_error_n;
     } else if ((bicgstab_num_outer_iterSOLVER>0)||
      	       (cg_loop==1)) {
      // do nothing
     } else
      amrex::Error("bicgstab_num_outer_iterSOLVER or cg_loop invalid");

      // if cg_loop==0 then error_n=error0 
      // error0=error after Jacobi iter. IF num_outer_iterSOLVER==0
     check_outer_solver_convergence(error_n,error0,
             Ar_error_n,Ar_error0,
             rAr_error_n,rAr_error0,
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
        std::cout << "meets_tol=1 at top of CG or BICGSTAB: Ar_error_n=" <<
          Ar_error_n << '\n';
        std::cout << "meets_tol=1 at top of CG or BICGSTAB: rAr_error_n=" <<
          rAr_error_n << '\n';
        if (cg_loop==0)
         std::cout << "error_after_all_jacobi_sweeps (jacobi method)=" <<
          error_after_all_jacobi_sweeps << '\n';
       } // ioproc
      }  // verbose>0

      error_n=local_error_n;
      Real local_tol=1.1*save_mac_abs_tol;
      meets_tol=0;
      check_outer_solver_convergence(error_n,error0,
             Ar_error_n,Ar_error0,
             rAr_error_n,rAr_error0,
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
     int gmres_precond_iter=gmres_precond_iter_base;

     Real dnorm=0.0;

       // variables initialized to 0.0
     allocate_array(1,nsolveMM,-1,Z_MF);
     allocate_array(0,nsolveMM,-1,P_MF);
     allocate_array(0,nsolveMM,-1,bicg_R0hat_MF);
     allocate_array(1,nsolveMM,-1,bicg_U0_MF);
     allocate_array(0,nsolveMM,-1,bicg_V0_MF);
     allocate_array(0,nsolveMM,-1,bicg_P1_MF);
     allocate_array(0,nsolveMM,-1,bicg_R1_MF);
     allocate_array(1,nsolveMM,-1,bicg_Y_MF);
     allocate_array(1,nsolveMM,-1,bicg_Hvec_MF);
     allocate_array(0,nsolveMM,-1,bicg_S_MF);
     allocate_array(0,nsolveMM,-1,bicg_T_MF);


     copyALL(0,nsolveMM,bicg_R0hat_MF,CGRESID_MF);  // R0hat=CGRESID(R0)
      // MAC_PHI_CRSE(U1)=0.0
     zeroALL(1,nsolveMM,bicg_U0_MF);

     error_history.resize(vcycle_max+1);
     for (int ehist=0;ehist<error_history.size();ehist++) {
      error_history[ehist][0]=0.0;
      error_history[ehist][1]=0.0;
      error_history[ehist][2]=0.0;
      error_history[ehist][3]=0.0;
     }

#if (profile_solver==1)
     bprof.stop();
#endif

     int BICGSTAB_ACTIVE=0;
     if ((enable_spectral==1)||  // SEM space and time
         (enable_spectral==2)) { // SEM space
      BICGSTAB_ACTIVE=1;
     } else if ((enable_spectral==0)||
                (enable_spectral==3)) {  // SEM time

      if (project_option==3) { //viscosity
       BICGSTAB_ACTIVE=1;
      } else if (project_option_is_valid(project_option)==1) {

       if (always_use_bicgstab==1) {
        BICGSTAB_ACTIVE=1;
       } else if (always_use_bicgstab==0) {
        BICGSTAB_ACTIVE=0;
       } else
        amrex::Error("always_use_bicgstab invalid");

      } else
       amrex::Error("project_option invalid");
     } else
      amrex::Error("enable_spectral invalid 3");

     for (vcycle=0;((vcycle<=vcycle_max)&&(meets_tol==0));vcycle++) {

#if (profile_solver==1)
      bprof.start();
#endif

        // CGRESID(R0) is the residual when using U0
        // MAC_PHI_CRSE(U1) and U0 are the same at this point.
      dot_productALL(project_option,CGRESID_MF,CGRESID_MF,error_n,nsolve);
      if (error_n>=0.0) {
       error_n=sqrt(error_n);
      } else
       amrex::Error("error_n invalid");

      adjust_tolerance(error_n,error0_max,project_option);

      copyALL(0,nsolveMM,P_SOLN_MF,CGRESID_MF); //P_SOLN=CGRESID_MF
       // V1=A P_SOLN
       // 1. (begin)calls project_right_hand_side(P)
       // 2. (end)  calls project_right_hand_side(V1)
      applyALL(project_option,P_SOLN_MF,bicg_V1_MF,nsolve);
      dot_productALL(project_option,CGRESID_MF,bicg_V1_MF,rAr_error_n,nsolve);
      dot_productALL(project_option,bicg_V1_MF,bicg_V1_MF,Ar_error_n,nsolve);

      if (rAr_error_n<0.0) {
       rAr_error_n=0.0;
      } else if (rAr_error_n>=0.0) {
       rAr_error_n=sqrt(rAr_error_n);
      } else
       amrex::Error("rAr_error_n invalid");

      if (Ar_error_n>=0.0) {
       Ar_error_n=sqrt(Ar_error_n);
      } else
       amrex::Error("Ar_error_n invalid");

      error_history[vcycle][0]=error_n;
      error_history[vcycle][1]=rAr_error_n;
      error_history[vcycle][1]=Ar_error_n;
      error_history[vcycle][3]=save_mac_abs_tol;

      Real restart_tol=0.0;
      restart_flag=0;

      check_outer_solver_convergence(error_n,error0,
             Ar_error_n,Ar_error0,
             rAr_error_n,rAr_error0,
             save_mac_abs_tol,
             meets_tol);

      if (verbose>0) {
       if (ParallelDescriptor::IOProcessor()) {
        std::cout << "cg_loop,vcycle,E0,En,rArEn,ArEn " << cg_loop << ' ' <<
         vcycle << ' ' << error0 << ' ' << error_n << ' ' <<
         rAr_error_n << ' ' << Ar_error_n << '\n';
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
         if (rho0>restart_tol) {
          beta=rho1/rho0;
          // P_MF=0 initially or on restart.
          // P_MF=Z_MF + beta P_MF
          mf_combine(project_option,
           Z_MF,P_MF,beta,P_MF,nsolve);
          change_flag=0;
          project_right_hand_side(P_MF,project_option,change_flag);

          copyALL(0,nsolveMM,P_SOLN_MF,P_MF); //P_SOLN=P
          // V1=A P_SOLN
          // 1. (begin)calls project_right_hand_side(P)
          // 2. (end)  calls project_right_hand_side(V1)
          applyALL(project_option,P_SOLN_MF,bicg_V1_MF,nsolve);

          Real pAp=0.0;
          dot_productALL(project_option,P_MF,bicg_V1_MF,pAp,nsolve);

          if (pAp>restart_tol) {
           alpha=rho1/pAp;

            // mac_phi_crse(U1)=mac_phi_crse+alpha P
           mf_combine(project_option,
             MAC_PHI_CRSE_MF,P_MF,alpha,MAC_PHI_CRSE_MF,nsolve);
           change_flag=0;
           project_right_hand_side(MAC_PHI_CRSE_MF,
             project_option,change_flag);
            // U0=MAC_PHI_CRSE(U1)
           copyALL(1,nsolveMM,bicg_U0_MF,MAC_PHI_CRSE_MF);

           // CGRESID=RHS-A mac_phi_crse(U1)
           // 1. (start) calls project_right_hand_side(MAC_PHI_CRSE_MF)
           // 2. (end)   calls project_right_hand_side(CGRESID_MF)
           residALL(project_option,MAC_RHS_CRSE_MF,CGRESID_MF,
            MAC_PHI_CRSE_MF,nsolve);

          } else if ((pAp>=0.0)&&(pAp<=restart_tol)) {
           meets_tol=1;

           // this case can happen due to round off error when
           // A is a singular (indefinite) matrix
          } else if (pAp<0.0) { 
           meets_tol=1;
          } else {
           std::cout << "pAp= " << pAp << endl;
           amrex::Error("pAp invalid in main solver");
          }
         } else if ((rho0>=0.0)&&(rho0<=restart_tol)) {
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

       } else if (BICGSTAB_ACTIVE==1) { //MG-GMRES PCG

         // rho1=R0hat dot CGRESID(R0)
         //  =(b-A x0) dot (b-A xn) =
         //  b^T b + x0^T A^T A xn - x0^T A^T b - b^T A xn
        dot_productALL(project_option,bicg_R0hat_MF,CGRESID_MF,rho1,nsolve);

        if (vcycle==0) { // R0hat=R when vcycle==0

         if (rho1>0.0) {
          Real sanity_error=sqrt(rho1);
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

        if (rho0<=restart_tol) {
         restart_flag=1;
        } else if (rho0>restart_tol) {
         // do nothing
        } else
         amrex::Error("rho0 failed Nav3");

        if (w0<=restart_tol) {
         restart_flag=1;
        } else if (w0>restart_tol) {
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
         multiphase_GMRES_preconditioner(
          gmres_precond_iter,
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

	 if (alpha<=restart_tol) {
	  restart_flag=1;
	 } else if (alpha>restart_tol) {
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
          copyALL(1,nsolveMM,MAC_PHI_CRSE_MF,bicg_Hvec_MF);

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
           dnorm=sqrt(dnorm);
	  } else
           amrex::Error("dnorm invalid Nav3");

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
           multiphase_GMRES_preconditioner(
            gmres_precond_iter,
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

	   if (a2>restart_tol) {
   	    // do nothing
	   } else if ((a2>=0.0)&&(a2<=restart_tol)) {
//	    restart_flag=1;
            meets_tol=1;
	   } else {
	    amrex::Error("a2 invalid Nav3");
	   }

	   if (a1>0.0) {
	    // do nothing
	   } else if (a1<=0.0) {
//	    restart_flag=1;
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
           dnorm=sqrt(dnorm);
	  } else
	   amrex::Error("dnorm invalid Nav3");

          w0=w1;
          // CGRESID(R0)=R1
          copyALL(0,nsolveMM,CGRESID_MF,bicg_R1_MF);
          copyALL(0,nsolveMM,P_MF,bicg_P1_MF);
          copyALL(0,nsolveMM,bicg_V0_MF,bicg_V1_MF);
          // U0=MAC_PHI_CRSE(U1)
          copyALL(1,nsolveMM,bicg_U0_MF,MAC_PHI_CRSE_MF);

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

       if (restart_flag==0) {
        // do nothing
       } else if (restart_flag==1) {

#if (profile_solver==1)
        bprof.start();
#endif

	if (verbose>0) {
         if (ParallelDescriptor::IOProcessor()) {
          std::cout << "WARNING:RESTARTING: bicgstab vcycle " 
	         << vcycle << '\n';
          std::cout << "RESTARTING: BICGSTAB_ACTIVE=" << 
           BICGSTAB_ACTIVE << '\n';
          std::cout << "RESTARTING: gmres_precond_iter= " << 
           gmres_precond_iter << '\n';
          std::cout << "RESTARTING: local_presmooth= " << 
           local_presmooth << '\n';
          std::cout << "RESTARTING: local_postsmooth= " << 
           local_postsmooth << '\n';
          std::cout << "RESTARTING: error_history[vcycle][0,1,2,3]= " << 
           error_history[vcycle][0] << ' ' <<
	   error_history[vcycle][1] << ' ' <<
	   error_history[vcycle][2] << ' ' <<
	   error_history[vcycle][3] << '\n';
          std::cout << "RESTARTING: project_option= " << 
	         project_option << '\n';
         }
	} else if (verbose==0) {
	 // do nothing
	} else
	 amrex::Error("verbose invalid");

        // CGRESID(R0)=RHS-A U0
	// 1. (start) calls project_right_hand_side(bicg_U0_MF)
	// 2. (end)   calls project_right_hand_side(CGRESID_MF)
        residALL(project_option,MAC_RHS_CRSE_MF,CGRESID_MF,
          bicg_U0_MF,nsolve);
        // R0hat=CGRESID(R0)
        copyALL(0,nsolveMM,bicg_R0hat_MF,CGRESID_MF);
        // MAC_PHI_CRSE(U1)=U0
        copyALL(1,nsolveMM,MAC_PHI_CRSE_MF,bicg_U0_MF);
        rho0=1.0;
        rho1=1.0;
        alpha=1.0;
        w0=1.0;
	// dnorm=RESID dot RESID
        dot_productALL(project_option,CGRESID_MF,CGRESID_MF,dnorm,nsolve);
	if (dnorm>=0.0) {
         dnorm=sqrt(dnorm);
	} else
         amrex::Error("dnorm invalid Nav3");

        zeroALL(0,nsolveMM,bicg_V0_MF);
        zeroALL(0,nsolveMM,P_MF);

#if (profile_solver==1)
        bprof.stop();
#endif

       } else
        amrex::Error("restart_flag invalid");

      } // meets_tol==0

      // top level: BiCGStab  
      // preconditioner: preconditioned GMRES
      // GMRES, like CG or BiCGStab, is a Krylov subspace method in which
      // the solutions is approximated as
      // sum_i=0^M alpha_i A^{i}r
      // (r is the residual)
      // preconditioner for GMRES: multigrid
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
       std::cout << "->project_option= " << project_option << '\n';
       for (int ehist=0;ehist<error_history.size();ehist++) {
        std::cout << "vcycle " << ehist << 
         " error_history[vcycle][0,1,2,3] " <<
         error_history[ehist][0] << ' ' <<
         error_history[ehist][1] << ' ' <<
         error_history[ehist][2] << ' ' <<
         error_history[ehist][3] << '\n';
       }
       for (int ehist=0;ehist<outer_error_history.size();ehist++) {
        std::cout << "outer_iter " << ehist << 
         " outer_error_history[vcycle][0,1,2,3] " <<
         outer_error_history[ehist][0] << ' ' <<
         outer_error_history[ehist][1] << ' ' <<
         outer_error_history[ehist][2] << ' ' <<
         outer_error_history[ehist][3] << '\n';
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
     std::cout << "->project_option= " << project_option << '\n';
     std::cout << "->error_n= " << error_n << '\n';
     std::cout << "->cg_loop_max= " << cg_loop_max << '\n';
     std::cout << "->ERROR HISTORY " << cg_loop_max << '\n';
     for (int ehist=0;ehist<error_history.size();ehist++) {
      std::cout << "vcycle " << ehist << 
       " error_history[vcycle][0,1,2,3] " <<
       error_history[ehist][0] << ' ' <<
       error_history[ehist][1] << ' ' <<
       error_history[ehist][2] << ' ' <<
       error_history[ehist][3] << '\n';
     }
     for (int ehist=0;ehist<outer_error_history.size();ehist++) {
      std::cout << "outer_iter " << ehist << 
       " outer_error_history[vcycle][0,1,2,3] " <<
       outer_error_history[ehist][0] << ' ' <<
       outer_error_history[ehist][1] << ' ' <<
       outer_error_history[ehist][2] << ' ' <<
       outer_error_history[ehist][3] << '\n';
     }
    } else if (vcycle>=0) {
     // do nothing
    } else {
     amrex::Error("vcycle bust");
    }

    std::fflush(NULL);

    for (int ilist=0;ilist<scomp.size();ilist++) 
     avgDownALL(state_index,scomp[ilist],ncomp[ilist],1);

      // override_bc_to_homogeneous=1
      // call FORT_OVERRIDEBC
    CPP_OVERRIDEPBC(1,project_option);

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
     energyflag=0; // energyflag=2 => GRADPEDGE=gradp
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
      project_option,nsolve);

     if (ilev<finest_level) {
      int ncomp_edge=-1;
      ns_level.avgDownEdge_localMF(GRADPEDGE_MF,
       0,ncomp_edge,0,AMREX_SPACEDIM,1,16);
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
      ns_level.Copy_localMF(MAC_TEMP_MF+dir,UMAC_MF+dir,0,0,nsolveMM_FACE,0);

     MultiFab* snew_mf=nullptr;
     if (state_index==State_Type) {
      snew_mf=ns_level.getState_list(1,scomp,ncomp,cur_time_slab);
     } else if (state_index==DIV_Type) {
      if (scomp.size()!=1)
       amrex::Error("scomp.size() invalid");
      if (ncomp[0]!=num_materials_face)
       amrex::Error("ncomp[0] invalid");
      snew_mf=ns_level.getStateDIV_DATA(1,scomp[0],ncomp[0],cur_time_slab);
     } else {
      snew_mf=nullptr;
      amrex::Error("state_index invalid");
     }

     if (snew_mf->nComp()!=nsolveMM)
      amrex::Error("snew_mf->nComp() invalid");

     MultiFab::Add(
      *ns_level.localMF[OUTER_ITER_PRESSURE_MF],*snew_mf,0,0,nsolveMM,0);

     delete snew_mf;

     ns_level.zero_independent_variable(project_option,nsolve);
    }  // ilev=finest_level ... level

    change_flag=0;
    project_right_hand_side(OUTER_ITER_PRESSURE_MF,project_option,
      	     change_flag);

    delete_array(PRESPC_MF);

       // variables initialized to 0.0
    allocate_array(1,nsolveMM,-1,OUTER_MAC_PHI_CRSE_MF);
    allocate_array(0,nsolveMM,-1,OUTER_RESID_MF);
    allocate_array(0,nsolveMM,-1,OUTER_MAC_RHS_CRSE_MF);

    allocate_array(1,nsolveMM,-1,P_SOLN_MF);
    allocate_array(0,nsolveMM,-1,bicg_V1_MF);

    for (int ilev=finest_level;ilev>=level;ilev--) {
     NavierStokes& ns_level=getLevel(ilev);
     if (ilev<finest_level) {
      ns_level.setVal_localMF(OUTER_MAC_RHS_CRSE_MF,0.0,0,nsolveMM,0); 
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
     // if local_solvability_projection, then this
     // routine modifies OUTER_RESID so that it sums to zero.
     // if singular_possible, then this routine zeros out the
     // residual where the matrix diagonal (prior to dual time stepping
     // modification) is 0.
    residALL(project_option,OUTER_MAC_RHS_CRSE_MF,
      OUTER_RESID_MF,OUTER_MAC_PHI_CRSE_MF,nsolve);
    outer_error=0.0;
    Real rAr_outer_error=0.0;
    Real Ar_outer_error=0.0;
    dot_productALL(project_option,
      OUTER_RESID_MF,OUTER_RESID_MF,outer_error,nsolve);
    if (outer_error>=0.0) {
     outer_error=sqrt(outer_error);
    } else
     amrex::Error("outer_error invalid");

    copyALL(0,nsolveMM,P_SOLN_MF,OUTER_RESID_MF); //P_SOLN=OUTER_RESID_MF
       // V1=A P_SOLN
       // 1. (begin)calls project_right_hand_side(P)
       // 2. (end)  calls project_right_hand_side(V1)
    applyALL(project_option,P_SOLN_MF,bicg_V1_MF,nsolve);
    dot_productALL(project_option,OUTER_RESID_MF,bicg_V1_MF,
      	     rAr_outer_error,nsolve);
    dot_productALL(project_option,bicg_V1_MF,bicg_V1_MF,
      	     Ar_outer_error,nsolve);

    if (rAr_outer_error<0.0) {
     rAr_outer_error=0.0;
    } else if (rAr_outer_error>=0.0) {
     rAr_outer_error=sqrt(rAr_outer_error);
    } else
     amrex::Error("rAr_outer_error invalid");

    if (Ar_outer_error>=0.0) {
     Ar_outer_error=sqrt(Ar_outer_error);
    } else
     amrex::Error("Ar_outer_error invalid");

    delete_array(bicg_V1_MF);
    delete_array(P_SOLN_MF);

    delete_array(OUTER_MAC_RHS_CRSE_MF);
    delete_array(OUTER_MAC_PHI_CRSE_MF);
    delete_array(OUTER_RESID_MF);

    Real outer_tol=100.0*save_mac_abs_tol;

    outer_error_history[bicgstab_num_outer_iterSOLVER][0]=outer_error;
    outer_error_history[bicgstab_num_outer_iterSOLVER][1]=rAr_outer_error;
    outer_error_history[bicgstab_num_outer_iterSOLVER][2]=Ar_outer_error;
    outer_error_history[bicgstab_num_outer_iterSOLVER][3]=outer_tol;

    bicgstab_num_outer_iterSOLVER++;

    if (verbose>0) {
     if (ParallelDescriptor::IOProcessor()) {
      std::cout << "project_option= " << project_option << '\n';
      std::cout << "bicgstab_num_outer_iterSOLVER,E,rAr_E,Ar_E " << 
       bicgstab_num_outer_iterSOLVER << ' ' << outer_error << 
       ' ' << rAr_outer_error << 
       ' ' << Ar_outer_error << '\n';
     }
    }

    outer_iter_done=0;

    check_outer_solver_convergence(outer_error,error0,
             Ar_outer_error,Ar_error0,
             rAr_outer_error,rAr_error0,
             outer_tol,
             outer_iter_done);

    if (bicgstab_num_outer_iterSOLVER>bicgstab_max_num_outer_iter)
     outer_iter_done=1;

     // We cannot do iterative refinement (see Burden and Faires -
     // iterative techniques in matrix algebra) in this case since:
     // RHSPROJ( div(U - dt grad P) ) <> 
     // RHSPROJ( RHSPROJ(div U)- dt div grad P ) 
    if (local_solvability_projection==1) {
     if (bicgstab_num_outer_iterSOLVER>1)
      outer_iter_done=1;
    } else if (local_solvability_projection==0) {
     // do nothing
    } else
     amrex::Error("local solvability_projection invalid");

    if (bicgstab_num_outer_iterSOLVER<min_bicgstab_outer_iter)
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
   std::cout << "SOLVER STATISTICS  TIME, project_option = " <<
          cur_time_slab << " " << project_option << '\n';
   std::cout << "project_option= " << project_option <<
          " SDC_outer_sweeps= " << SDC_outer_sweeps <<
          " slab_step= " << slab_step << '\n';

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
  }

 } else
  amrex::Error("number_solver_calls.size() invalid");

 deallocate_maccoefALL(project_option);

    // copy OUTER_ITER_PRESSURE to s_new
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.putState_localMF_list(OUTER_ITER_PRESSURE_MF,state_index,
   scomp,ncomp);
 }  // ilev=finest_level ... level

 for (int ilist=0;ilist<scomp.size();ilist++) 
  avgDownALL(state_index,scomp[ilist],ncomp[ilist],1);

 int homflag_dual_time=0;

 if ((project_option==1)||   // initial project
     (project_option==11)) { // FSI_material_exists (last project)
  homflag_dual_time=1;
 } else if (project_option==12) { // pressure extrapolation
  homflag_dual_time=0;
 } else if ((project_option==0)||  //regular project
            (project_option==2)) { //thermal conductivity
  homflag_dual_time=0;
 } else if (project_option==3) { // viscosity
  homflag_dual_time=0;
 } else if ((project_option>=100)&&
	    (project_option<100+num_species_var)) {//species
  homflag_dual_time=0;
 } else if (project_option==200) { //smoothing
  homflag_dual_time=0;
 } else
  amrex::Error("project_option invalid 53");

 CPP_OVERRIDEPBC(homflag_dual_time,project_option);

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

 if (project_option_projection(project_option)==1) {

  getState_localMF_listALL(
    PRESPC2_MF,1,
    state_index,
    scomp,
    ncomp);

  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);

     // presBILINEAR2 and get_new_data(Umac_type+dir) are inputs.
     // Update cell velocity, density (if non-conservative),
     // total Energy, and internal energy.
     // initializes ns_level.conservative_energy_mask

     // (multiphase_project)
     // If update_energy=1, then rho (if non-cons),E,T,consmask updated.
     // Copies UMAC to Umac_new: "save_to_macvel_state(UMAC_MF)"
   int idx_gpcell=-1;
   int idx_divup=-1;
   int update_energy=0;
   if (project_option==0) {
    update_energy=1; // update temperature and density (if non-cons)
   } else if ((project_option==1)||   //initial project
	      (project_option==11)) { //FSI_material_exists(last project)
    // do nothing
   } else 
    amrex::Error("project_option invalid");

    // NavierStokes::apply_cell_pressure_gradient declared
    // in NavierStokes2.cpp.
    // At the end of apply_cell_pressure_gradient,
    // UMAC_MF is copied to UMAC_new. 
   ns_level.apply_cell_pressure_gradient(project_option,
    update_energy,PRESPC2_MF,UMAC_MF,idx_gpcell,idx_divup);

   if (project_option==0) {

    int project_option_combine=2;  // temperature in multiphase_project
    int combine_flag=2; //combine if vfrac<VOFTOL
    int hflag=0;
     // combine_idx==-1 => update S_new  
     // combine_idx>=0  => update localMF[combine_idx]
    int combine_idx=-1;  
    int update_flux=0;
    int interface_cond_avail=0;

    ns_level.combine_state_variable(
     project_option_combine,
     combine_idx,
     combine_flag,
     hflag,
     update_flux,
     interface_cond_avail); 
    for (int ns=0;ns<num_species_var;ns++) {
     project_option_combine=100+ns; // species
     ns_level.combine_state_variable(
      project_option_combine,
      combine_idx,
      combine_flag,
      hflag,
      update_flux,
      interface_cond_avail); 
    }

      // velocity and pressure
    ns_level.avgDown(State_Type,0,
     num_materials_vel*(AMREX_SPACEDIM+1),1);
    ns_level.avgDown(State_Type,scomp_den,num_state_material*nmat,1);

   } else if (project_option==1) {
    ns_level.avgDown(State_Type,0,num_materials_vel*AMREX_SPACEDIM,1);
   } else if (project_option==11) { //FSI_material_exists (last project)
    ns_level.avgDown(State_Type,0,num_materials_vel*AMREX_SPACEDIM,1);
   } else
    amrex::Error("project_option invalid 54");

  } // ilev=finest_level ... level

  delete_array(PRESPC2_MF);

  if (project_option_needs_scaling(project_option)==1) {

   unscale_variablesALL();

  } else if (project_option==1) {
   // do nothing
  } else
   amrex::Error("project_option invalid");

  if (project_option==0) {
 
   int do_alloc=1;
   int simple_AMR_BC_flag_viscosity=1;
   int im_tensor=-1;
   int idx_elastic_flux=-1;
   init_gradu_tensorALL(
     im_tensor,
     HOLD_VELOCITY_DATA_MF,
     do_alloc,
     CELLTENSOR_MF,
     FACETENSOR_MF,
     idx_elastic_flux,
     simple_AMR_BC_flag_viscosity);

   for (int ilev=finest_level;ilev>=level;ilev--) {
    NavierStokes& ns_level=getLevel(ilev);
    ns_level.init_pressure_error_indicator();  
   }
   avgDownError_ALL();

   delete_array(FACETENSOR_MF);
   delete_array(CELLTENSOR_MF);

  } else if ((project_option==1)||
             (project_option==11)) { //FSI_material_exists (last project)
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

  // project_option tells whether pressure (0,1), 
  // temperature (2), velocity (3,4,5), or species.
 CPP_OVERRIDEPBC(0,project_option);

 if (project_option==12) {  // pressure extrapolation
  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
    // in: MacProj.cpp (calls FORT_RESTORE_PRES)
   ns_level.restore_active_pressure(PRESSURE_SAVE_MF);
   int pcomp=num_materials_vel*AMREX_SPACEDIM;
   ns_level.avgDown(State_Type,pcomp,1,1); // average from ilev+1 to ilev
  }
  delete_array(PRESSURE_SAVE_MF);
 }

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.remove_project_variables(); // remove poldhold,ones,outer_iter
  ns_level.remove_pressure_work_vars(); // remove umacstar, gradpedge, pedge
  ns_level.delete_localMF(FACE_WEIGHT_MF,AMREX_SPACEDIM);
  ns_level.delete_localMF(OFF_DIAG_CHECK_MF,1);
 }

 delete_array(MAC_PHI_CRSE_MF);
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

 if (project_option==11) { //FSI_material_exists (last project)
  for (int ilev=level;ilev<=finest_level;ilev++) {
   NavierStokes& ns_level=getLevel(ilev);
   MultiFab& DIV_new=ns_level.get_new_data(DIV_Type,slab_step+1);
   MultiFab::Copy(
      DIV_new,
      *ns_level.localMF[DIV_SAVE_MF],
      0,0,1,1);
  } // ilev=level ... finest_level
  delete_array(DIV_SAVE_MF);
 } else if (project_option==200) {  // temperature smoothing

  if ((localMF_grow[TEMPERATURE_SAVE_MF]>=0)&&
      (localMF_grow[TEMPERATURE_SMOOTH_MF]==-1)) {
   // do nothing
  } else
   amrex::Error("TEMPERATURE_SAVE_MF or TEMPERATURE_SMOOTH_MF invalid");

  Vector<int> scomp_thermal;
  Vector<int> ncomp_thermal;
  int state_index_thermal;  
  int ncomp_check_thermal;
  get_mm_scomp_solver(
    nmat,
    project_option,
    state_index_thermal,
    scomp_thermal,
    ncomp_thermal,
    ncomp_check_thermal);
  if (ncomp_check_thermal!=nmat)
   amrex::Error("ncomp_check_thermal invalid");

  int project_option_combine=2; 
  int combine_flag=1;  // GFM -> FVM  
   // combine_idx==-1 => update S_new  
   // combine_idx>=0  => update localMF[combine_idx]
  int combine_idx=-1; 
  int update_flux=0;
  int hflag_combine=0;
  int interface_cond_avail=0;
 
  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   ns_level.combine_state_variable(
    project_option_combine,
    combine_idx,
    combine_flag,
    hflag_combine,
    update_flux,
    interface_cond_avail); 
  } // ilev=finest_level ... level

  for (int ilist=0;ilist<scomp_thermal.size();ilist++) 
   avgDownALL(state_index_thermal,scomp_thermal[ilist],
	ncomp_thermal[ilist],1);

   // data at time = cur_time_slab
  getState_localMF_listALL(
    TEMPERATURE_SMOOTH_MF,
    normal_probe_size+3,
    state_index_thermal,
    scomp_thermal,
    ncomp_thermal);

  FIX ME restore state here and delete the temp var
 } else if (project_option_is_valid(project_option)==1) {
  // do nothing
 } else {
  amrex::Error("project_option invalid");
 }

#if (profile_solver==1)
 bprof.stop();
#endif

 std::fflush(NULL);

}  // end subroutine multiphase_project


void NavierStokes::diffusion_heatingALL(
  int source_idx,int idx_heat) {

 int finest_level=parent->finestLevel();

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.resize_levelsetLO(2,LEVELPC_MF);
  ns_level.debug_ngrow(LEVELPC_MF,2,840);
  ns_level.VOF_Recon_resize(1,SLOPE_RECON_MF);
  ns_level.debug_ngrow(SLOPE_RECON_MF,1,841);
  ns_level.debug_ngrow(source_idx,1,842);
  ns_level.debug_ngrow(idx_heat,0,842);
  ns_level.debug_ngrow(FACE_VAR_MF,0,844);
  ns_level.resize_metrics(1);
  ns_level.debug_ngrow(VOLUME_MF,1,845);
 }

 if (geom.IsRZ()) {
  //rzflag=1
  if (AMREX_SPACEDIM!=2)
   amrex::Error("dimension bust");
 } else if (geom.IsCartesian()) {
  //rzflag=0
 } else if (geom.IsCYLINDRICAL()) {
  //rzflag=3
 } else
  amrex::Error("CoordSys bust 61");

 min_face_wt.resize(thread_class::nthreads);
 max_face_wt.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  min_face_wt[tid].resize(4);
  max_face_wt[tid].resize(4);
  for (int iwt=0;iwt<4;iwt++) {
   min_face_wt[tid][iwt]=1.0e+20;
   max_face_wt[tid][iwt]=-1.0e+20;
  }
 } // tid


 int nsolve=AMREX_SPACEDIM;
 int project_option=3; // viscosity
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.allocate_FACE_WEIGHT(nsolve,project_option);
 }

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.diffusion_heating(source_idx,idx_heat);
 } // ilev

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.delete_localMF(FACE_WEIGHT_MF,AMREX_SPACEDIM);
  ns_level.delete_localMF(OFF_DIAG_CHECK_MF,1);
 }

}  // subroutine diffusion_heatingALL

void NavierStokes::avgDownALL_TENSOR() {

 int nmat=num_materials;
 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=nmat)) {
   // spectral_override==1 => order derived from "enable_spectral"
   // spectral_override==0 => always low order.
  avgDownALL(Tensor_Type,0,
    num_materials_viscoelastic*NUM_CELL_ELASTIC,0);

  if (MAC_grid_displacement==0) {
   // do nothing
  } else if (MAC_grid_displacement==1) {
   if (face_flag==1) {
    avgDownMacState(XDmac_Type,0);
   } else 
    amrex::Error("expecting face_flag==1 if MAC_grid_displacement==1");
  } else
   amrex::Error("MAC_grid_displacement invalid");

 } else
  amrex::Error("num_materials_viscoelastic invalid");

} // subroutine avgDownALL_TENSOR

// VISCOELASTIC, CTML FORCE
void NavierStokes::vel_elastic_ALL() {

 int nmat=num_materials;
 int finest_level=parent->finestLevel();

 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=nmat)) {

  for (int im=0;im<nmat;im++) {
   if ((particleLS_flag[im]==1)||
       (particleLS_flag[im]==0)) { 
    if (ns_is_rigid(im)==0) {
     if ((elastic_time[im]>0.0)&&
         (elastic_viscosity[im]>0.0)) {

      if (store_elastic_data[im]==1) {

       if (MAC_grid_displacement==0) {

        for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
         allocate_array(0,AMREX_SPACEDIM*AMREX_SPACEDIM,
          dir,XDISP_FLUX_MF+dir);
         setVal_array(0,AMREX_SPACEDIM*AMREX_SPACEDIM,
          0.0,XDISP_FLUX_MF+dir);
        }

        if (viscoelastic_model[im]==2) {

 	 int push_enable_spectral=enable_spectral;
	 int elastic_enable_spectral=0;
	 override_enable_spectral(elastic_enable_spectral);

	 int do_alloc=1;
         int simple_AMR_BC_flag_viscosity=1;
	 init_gradu_tensorALL(
          im, // 0<=im<=nmat-1  (signifies displacement input)
	  XDISPLACE_MF, // deleted in init_gradu_tensorALL since do_alloc==1
	  do_alloc,
	  CELLTENSOR_MF,
	  FACETENSOR_MF,
	  XDISP_FLUX_MF, // elastic_idx
          simple_AMR_BC_flag_viscosity);

	 override_enable_spectral(push_enable_spectral);
         delete_array(CELLTENSOR_MF);
         delete_array(FACETENSOR_MF);

        } else if ((viscoelastic_model[im]==1)||
   		   (viscoelastic_model[im]==0)||
		   (viscoelastic_model[im]==3)) { // incremental
         // do nothing
        } else
        amrex::Error("viscoelastic_model[im] invalid");

        for (int ilev=finest_level;ilev>=level;ilev--) {
         NavierStokes& ns_level=getLevel(ilev);
         // note: tensor_advection_updateALL is called before veldiffuseALL.
         // VISCOTEN_MF initialized in NavierStokes::make_viscoelastic_tensor
         ns_level.make_viscoelastic_tensor(im);
         ns_level.make_viscoelastic_force(im);
        }

        delete_array(VISCOTEN_MF);

        for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
         delete_array(XDISP_FLUX_MF+dir);

       } else if (MAC_grid_displacement==1) {

        for (int ilev=finest_level;ilev>=level;ilev--) {
         NavierStokes& ns_level=getLevel(ilev);
         // note: tensor_advection_updateALL is called before veldiffuseALL.
         // VISCOTEN_MF initialized in NavierStokes::make_viscoelastic_tensor
         ns_level.make_viscoelastic_tensor(im);
         ns_level.MAC_GRID_ELASTIC_FORCE(im);
        }

        delete_array(VISCOTEN_MF);
       } else
        amrex::Error("MAC_grid_displacement invalid");

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
   } else
    amrex::Error("particleLS_flag[im] invalid");
  } // im=0..nmat-1
   
  if (MAC_grid_displacement==0) {

   // spectral_override==1 => order derived from "enable_spectral"
   avgDownALL(State_Type,0,
    num_materials_vel*(AMREX_SPACEDIM+1),1);

   // if filter_vel==0 and face_flag==1, then
   //  umacnew+=INTERP_TO_MAC(unew^CELL-register_mark^CELL) 
   // else if filter_vel==1 or face_flag==0, then
   //  umacnew=Interp_from_cell_to_MAC(unew^CELL)
   INCREMENT_REGISTERS_ALL(REGISTER_MARK_MF,1); 

  } else if (MAC_grid_displacement==1) {

   if (face_flag==1) {
     // average down the MAC velocity, set the boundary conditions.
    make_MAC_velocity_consistentALL();
    int use_VOF_weight=1;
    VELMAC_TO_CELLALL(use_VOF_weight);
   } else 
    amrex::Error("expecting face_flag==1 if MAC_grid_displacement==1");

  } else
   amrex::Error("MAC_grid_displacement invalid");

  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   ns_level.tensor_extrapolate(); // in: NavierStokes.cpp
  }
  avgDownALL_TENSOR();

   // register_mark=unew
  SET_STOKES_MARK(REGISTER_MARK_MF,101);

 } else
  amrex::Error("num_materials_viscoelastic invalid");


 if (CTML_FSI_flagC()==1) {

  // Add the solid force term on the right hand side
  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   ns_level.ctml_fsi_transfer_force();
  }

   // spectral_override==1 => order derived from "enable_spectral"
  avgDownALL(State_Type,0,
   num_materials_vel*(AMREX_SPACEDIM+1),1);

   // umacnew+=INTERP_TO_MAC(unew-register_mark)
   // (filter_vel==0 and face_flag==1)
  INCREMENT_REGISTERS_ALL(REGISTER_MARK_MF,2); 

    // register_mark=unew
  SET_STOKES_MARK(REGISTER_MARK_MF,102);

 } else if (CTML_FSI_flagC()==0) {
  // do nothing
 } else
  amrex::Error("CTML_FSI_flagC() invalid");

} // end subroutine vel_elastic_ALL


void NavierStokes::veldiffuseALL() {

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid");

 int nmat=num_materials;
 int finest_level=parent->finestLevel();

 int dencomp=num_materials_vel*(AMREX_SPACEDIM+1);
 int nden=nmat*num_state_material;

 int save_enable_spectral=enable_spectral;
 override_enable_spectral(viscous_enable_spectral);

 avgDownALL(State_Type,0,num_materials_vel*(AMREX_SPACEDIM+1),1);

 int convert_temperature=0;
 int convert_species=0;

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 if (num_materials_scalar_solve==1) {

  for (int im=0;im<nmat;im++) {
   if (heatviscconst[im]>0.0) {
    convert_temperature=1;
   } else if (heatviscconst[im]==0.0) {
    // do nothing
   } else
    amrex::Error("heatviscconst invalid");
  } // im 

  for (int im=0;im<nmat*num_species_var;im++) {
   if (speciesviscconst[im]>0.0) {
    convert_species=1;
   } else if (speciesviscconst[im]==0.0) {
    // do nothing
   } else
    amrex::Error("speciesviscconst invalid");
  } // im 


 } else if (num_materials_scalar_solve==nmat) {

  // do nothing

 } else
  amrex::Error("num_materials_scalar_solve invalid");

 if (verbose>0)
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "convert_temperature= " << convert_temperature << '\n';
   std::cout << "convert_species= " << convert_species << '\n';
  }
 
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  int hflag=0;
  ns_level.solid_temperature();  // if solid temperature is prescribed

  // MEHDI VAHAB HEAT SOURCE
  // NavierStokes.cpp: void NavierStokes::make_heat_source()
  // make_heat_source calls GODUNOV_3D.F90::FORT_HEATSOURCE which
  // calls PROB.F90::get_local_heat_source
  // if not supermesh algorithm, then the same temperature 
  // increment is added to all of the materials.
  ns_level.make_heat_source();  // updates S_new

  int project_option_combine=2;  // temperature in veldiffuseALL
  int combine_flag=0;  // FVM -> GFM 
   // combine_idx==-1 => update S_new  
   // combine_idx>=0  => update localMF[combine_idx]
  int combine_idx=-1; 
  int update_flux=0;
  int interface_cond_avail=1;
  
  if (convert_temperature==1) {
   combine_flag=0; // FVM -> GFM
  } else if (convert_temperature==0) { //thermal diffusivities all zero
   combine_flag=2;
  } else
   amrex::Error("convert_temperature invalid");

   // MEHDI VAHAB: COMBINE TEMPERATURES HERE IF NOT SUPERMESH APPROACH
  ns_level.combine_state_variable(
   project_option_combine,
   combine_idx,
   combine_flag,
   hflag,
   update_flux,
   interface_cond_avail); 

  for (int ns=0;ns<num_species_var;ns++) {
   project_option_combine=100+ns; // species

   if (convert_species==1) {
    combine_flag=0; // FVM -> GFM
   } else if (convert_species==0) {
    combine_flag=2;
   } else
    amrex::Error("convert_species invalid");

   ns_level.combine_state_variable(
    project_option_combine,
    combine_idx,
    combine_flag,
    hflag,
    update_flux,
    interface_cond_avail); 
  }
 }  // ilev=finest_level ... level

 avgDownALL(State_Type,dencomp,nden,1);

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.debug_ngrow(FACE_VAR_MF,0,870);
  ns_level.resize_metrics(1);
  ns_level.debug_ngrow(VOLUME_MF,1,872);
 }

  // in: veldiffuseALL
  // 1. substract dt_gradp_over_rho:
  //   viscosity: (1) (u^* - u^advect) = dt div 2\mu D/rho - dt grad p/rho
  //              u^advect <-- u^advect - dt grad p/rho then 
  //              solve 
  //              (u^* - u^advect) = dt div 2\mu D/rho which is then equivalent
  //              to (1)
  // 2. allocate scratch variables (including CONSERVE_FLUXES_MF)
  //
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  int alloc_flag=3;
  ns_level.alloc_gradp_over_rho(alloc_flag);
  ns_level.prepare_viscous_solver();
 }

 int project_option_temperature=2; // temperature conduction
 int vel_project_option=3; // viscosity
 int mom_force_project_option=4; // neg mom force

 Vector<int> scomp;
 Vector<int> ncomp;
 int ncomp_check;
 int state_index;

 get_mm_scomp_solver(
  num_materials_scalar_solve,
  project_option_temperature,
  state_index,
  scomp,
  ncomp,
  ncomp_check);

 int nsolve_thermal=1;
 int nsolveMM_thermal=nsolve_thermal*num_materials_scalar_solve;
 if (ncomp_check!=nsolveMM_thermal)
  amrex::Error("ncomp_check invalid");

 for (int im=0;im<nmat;im++) {
  if (override_density[im]==0) { // Drho/DT=-divu rho
   // check nothing
  } else if (override_density[im]==1) { // rho=rho(T,Y,z)
   // check nothing
  } else if (override_density[im]==2) { // P_hydro=P_hydro(rho(T,Y,Z))
    // convert_temperature==0 if all thermal diffusivities are zero.
    // P_hydro is expecting a reasonable temperature defined for each
    //  separate material?
   if ((convert_temperature==0)&&
       (num_materials_scalar_solve==1))
    amrex::Error("convert_temperature invalid");
  } else
   amrex::Error("override_density invalid");  
 } // im=0..nmat-1

   //ngrow=1
 getState_localMF_listALL(
   BOUSSINESQ_TEMP_MF,1,
   state_index,
   scomp,
   ncomp);

  // register_mark=unew (1 ghost)
 SET_STOKES_MARK(REGISTER_MARK_MF,103);

 show_norm2_id(REGISTER_MARK_MF,1);

 // -dt * |g| * beta * (T-T0)  beta<0
 // dt * beta * (T-T0) * omega^2 * r
 // u+=dt * (v^2/r +2 omega v)
 // v+=dt * (-uv/r -2 omega u)
 //
 // constant_viscosity==0:
 // u=u-dt * (1/r) * mu * (3 v_t/r + 2 u/r)/rho
 // v=v+dt * (1/r) * mu * (3 u_t/r - v/r) 
 //
 // constant_viscosity==1:
 // u=u-dt * (1/r) * mu * (2 v_t/r + u/r)/rho
 // v=v+dt * (1/r) * mu * (2 u_t/r - v/r) 
 //
 // U_t = F(U)  U=REGISTER_MARK_MF
 // 1. compute explicit terms: u^* = u^n + dt F1(u^n)
 // 2. compute implicit terms: u^** = u^* + dt F2(u^**)
 //

  // HOOP_FORCE_MARK_MF=(unp1-un)rho/dt
  // update_state==1:
  //  unp1(1)=unp1(1)/(one+param2*hoop_force_coef)
  //  unew=unp1
  // update_state==0:
  //  unp1(1)=unp1(1)-param2*hoop_force_coef*un(1)
 int update_state=1;
 diffuse_hoopALL(REGISTER_MARK_MF,BOUSSINESQ_TEMP_MF,
   HOOP_FORCE_MARK_MF,update_state);

  // for annulus problem:
  // rho cv (theta_t + u dot grad theta)= div k grad theta
  // T=theta + T_{1}(r)
  // if update_state==1,
  //  T_new=T_new-dt (-u T_{1}'(r))/(rho cv)
  // thermal_force=-(-u T_{1}'(r))
 thermal_transform_forceALL(REGISTER_MARK_MF,BOUSSINESQ_TEMP_MF,
   THERMAL_FORCE_MF,update_state);

 show_norm2_id(REGISTER_MARK_MF,2);

   // spectral_override==1 => order derived from "enable_spectral"
 avgDownALL(State_Type,0,
   num_materials_vel*(AMREX_SPACEDIM+1),1);

  // unew^MAC+=INTERP_TO_MAC(unew-register_mark)
  //    or
  // unew^MAC=INTERP_TO_MAC(unew) if filter_velocity[im]==1 or face_flag==0
 INCREMENT_REGISTERS_ALL(REGISTER_MARK_MF,3); 

 avgDownALL(State_Type,dencomp,nden,1);

  // register_mark=unew
 SET_STOKES_MARK(REGISTER_MARK_MF,104);

 //vel_elastic_ALL();

// -----------veldiffuseALL: viscosity -----------------------------


 if ((SDC_outer_sweeps>0)&&
     (SDC_outer_sweeps<ns_time_order)&&
     (divu_outer_sweeps+1==num_divu_outer_sweeps)) {

  if (ns_time_order>=2) {

   // UPDATESEMFORCE:
   // HOFAB=-div(2 mu D) - HOOP_FORCE_MARK_MF (update_state=0 at end of
   //                                          NavierStokes::do_the_advance)
   // unew=unew-(1/rho)(int (HOFAB) - dt (LOFAB))
   if ((viscous_enable_spectral==1)||  // SEM space and time
       (viscous_enable_spectral==3)) { // SEM time
    for (int ilev=finest_level;ilev>=level;ilev--) {
     NavierStokes& ns_level=getLevel(ilev);
     // calls: SEMDELTAFORCE in GODUNOV_3D.F90
     // does not look at: enable_spectral
     ns_level.make_SEM_delta_force(vel_project_option);
    }
   } else if (viscous_enable_spectral==0) {
    // do nothing
   } else
    amrex::Error("viscous_enable_spectral invalid");

   // spectral_override==1 => order derived from "enable_spectral"
   avgDownALL(State_Type,0,
    num_materials_vel*(AMREX_SPACEDIM+1),1);

   // umacnew+=INTERP_TO_MAC(unew-register_mark)
   // (filter_vel==0 and face_flag==1)
   INCREMENT_REGISTERS_ALL(REGISTER_MARK_MF,4); 

   // register_mark=unew
   SET_STOKES_MARK(REGISTER_MARK_MF,105);
  } else
   amrex::Error("ns_time_order invalid");

 } else if (SDC_outer_sweeps==0) {
  // do nothing
 } else if ((divu_outer_sweeps>=0)&&
            (divu_outer_sweeps+1<num_divu_outer_sweeps)) {
  // do nothing
 } else 
  amrex::Error("SDC_outer_sweeps or divu_outer_sweeps invalid");

 SET_STOKES_MARK(REGISTER_MARK_MF,106); //register_mark=unew

 show_norm2_id(REGISTER_MARK_MF,4);

  //multigrid-GMRES precond. BiCGStab viscosity
 multiphase_project(vel_project_option); 
 SET_STOKES_MARK(VISCHEAT_SOURCE_MF,107);

  // spectral_override==1 => order derived from "enable_spectral"
 avgDownALL(State_Type,0,
   num_materials_vel*(AMREX_SPACEDIM+1),1);

   // umacnew+=INTERP_TO_MAC(unew-register_mark)
   // (filter_vel==0 and face_flag==1)
 INCREMENT_REGISTERS_ALL(REGISTER_MARK_MF,5); 

 avgDownALL(State_Type,dencomp,nden,1);

 SET_STOKES_MARK(REGISTER_MARK_MF,108); //register_mark=unew

// ---------------- end viscosity ---------------------

 vel_elastic_ALL();

   // back to "enable_spectral" for momentum eqn
 override_enable_spectral(save_enable_spectral);

  // force at time = cur_time_slab
  // NEG_MOM_FORCE_MF=-(unp1-un)rho/dt
  // unew=unp1 if update_state==1
 update_state=1;
 mom_forceALL(NEG_MOM_FORCE_MF,update_state);

 if ((SDC_outer_sweeps>0)&&
     (SDC_outer_sweeps<ns_time_order)&&
     (divu_outer_sweeps+1==num_divu_outer_sweeps)) {

  if (ns_time_order>=2) {
   // UPDATESEMFORCE:
   // HOFAB=NEG_MOM_FORCE_MF (update_state=0 at end of 
   //                         NavierStokes::do_the_advance)
   // unew=unew-(1/rho)(int (HOFAB) - dt (LOFAB))
   for (int ilev=finest_level;ilev>=level;ilev--) {
    NavierStokes& ns_level=getLevel(ilev);
    // calls: SEMDELTAFORCE in GODUNOV_3D.F90
    // does not look at: enable_spectral
    ns_level.make_SEM_delta_force(mom_force_project_option);
   } // ilev=finest_level ... level

  } else
   amrex::Error("ns_time_order invalid");

 } else if (SDC_outer_sweeps==0) {
  // do nothing
 } else if ((divu_outer_sweeps>=0)&&
            (divu_outer_sweeps+1<num_divu_outer_sweeps)) {
  // do nothing
 } else 
  amrex::Error("SDC_outer_sweeps or divu_outer_sweeps invalid");

 save_enable_spectral=enable_spectral;
 override_enable_spectral(viscous_enable_spectral);

  // spectral_override==1 => not always low order
 avgDownALL(State_Type,0,
   num_materials_vel*(AMREX_SPACEDIM+1),1);

   // umacnew+=INTERP_TO_MAC(unew-register_mark)
   // (filter_vel==0 and face_flag==1)
 INCREMENT_REGISTERS_ALL(REGISTER_MARK_MF,6); 

   // register_mark=unew
 SET_STOKES_MARK(REGISTER_MARK_MF,109);

  // CONSERVATIVE SURFACE TENSION (Marangoni) FORCE
  // rhonew unew+=dt F_marangoni
  // 1. initialize CONSERVE_FLUXES_MF to store the force on the MAC grid and
  //    update the MAC velocity.
  //    Average down CONSERVE_FLUXES_MF
  // 2. interpolate F_marangoni from MAC to CELL and update the CELL velocity 
  //   
  // faces on the coarse grid or that neighbor a coarse grid 
  // use the old surface tension algorithm.
 for (int isweep=0;isweep<2;isweep++) {
  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   ns_level.make_marangoni_force(isweep);
   if (isweep==0) {
    // do nothing
   } else if (isweep==1) {
     // (i) avgDownMacState
     // (ii) fill boundary for Umac_new
    ns_level.make_MAC_velocity_consistent();
   } else
    amrex::Error("isweep invalid");
  } // ilev=finest_level ... level
 } // isweep=0,1

 avgDownALL(State_Type,0,
   num_materials_vel*(AMREX_SPACEDIM+1),1);

// ---------------- begin thermal diffusion ---------------------

   // UPDATESEMFORCE:
   // HOFAB=-div k grad T - THERMAL_FORCE_MF
   // Tnew=Tnew-(1/(rho cv))(int (HOFAB) - dt (LOFAB))
   // call to semdeltaforce with dcomp=slab_step*nstate_SDC+nfluxSEM
 if ((SDC_outer_sweeps>0)&&
     (SDC_outer_sweeps<ns_time_order)&&
     (divu_outer_sweeps+1==num_divu_outer_sweeps)) {

  if (ns_time_order>=2) {

   if ((viscous_enable_spectral==1)||   // SEM space and time
       (viscous_enable_spectral==3)) {  // SEM time

    for (int ilev=finest_level;ilev>=level;ilev--) {
     NavierStokes& ns_level=getLevel(ilev);
     // calls: SEMDELTAFORCE in GODUNOV_3D.F90
     ns_level.make_SEM_delta_force(project_option_temperature);
    }
   } else if (viscous_enable_spectral==0) {
    // do nothing
   } else
    amrex::Error("viscous_enable_spectral invalid");

   avgDownALL(State_Type,dencomp,nden,1);
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
 avgDownALL(State_Type,dencomp,nden,1);

 multiphase_project(project_option_temperature); // MGP BiCGStab temperature.

// --------------- end thermal diffusion -------------------



// ---------------begin save stable thermal diffusion and viscous forces

 if ((ns_time_order>=2)&&
     (ns_time_order<=32)&&
     (divu_outer_sweeps+1==num_divu_outer_sweeps)) {

  get_mm_scomp_solver(
   num_materials_scalar_solve,
   project_option_temperature,
   state_index,
   scomp,
   ncomp,
   ncomp_check);

  if (ncomp_check!=nsolveMM_thermal)
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
   // LOfab=-div(k grad T)-THERMAL_FORCE_MF
   // calls: UPDATESEMFORCE in GODUNOV_3D.F90
  if ((viscous_enable_spectral==1)||  // SEM space and time
      (viscous_enable_spectral==3)) { // SEM time
   update_SEM_forcesALL(project_option_temperature,PRESPC2_MF,
    update_spectralF,update_stableF);
  } else if (viscous_enable_spectral==0) {
   // do nothing
  } else
   amrex::Error("viscous_enable_spectral invalid");

  delete_array(PRESPC2_MF);


   // LOfab=-div(2 mu D)-HOOP_FORCE_MARK_MF
   // calls: UPDATESEMFORCE in GODUNOV_3D.F90
  if ((viscous_enable_spectral==1)||   // SEM space and time
      (viscous_enable_spectral==3)) {  // SEM time
   update_SEM_forcesALL(vel_project_option,VISCHEAT_SOURCE_MF,
    update_spectralF,update_stableF);
  } else if (viscous_enable_spectral==0) {
   // do nothing
  } else
   amrex::Error("viscous_enable_spectral invalid");

   // enable spectral for momentum equation.
  override_enable_spectral(save_enable_spectral);

   // LOfab=NEG_MOM_FORCE_MF
   // calls: UPDATESEMFORCE in GODUNOV_3D.F90
  update_SEM_forcesALL(mom_force_project_option,NEG_MOM_FORCE_MF,
   update_spectralF,update_stableF);

  save_enable_spectral=enable_spectral;
  override_enable_spectral(viscous_enable_spectral);

 } else if ((ns_time_order==1)||
            ((divu_outer_sweeps>=0)&&
             (divu_outer_sweeps+1<num_divu_outer_sweeps))) {
  // do nothing
 } else
  amrex::Error("ns_time_order or divu_outer_sweeps invalid");

// ---------------end save stable thermal diffusion and viscous forces

 for (int species_comp=0;species_comp<num_species_var;species_comp++) {

  int project_option_species=100+species_comp;

   // (rho Y)_{t} + div(rho u Y) =div (rho D) grad Y
   // [rho_t + div(rho u)]Y+rho DY/DT=div(rho D) grad Y
   // DY/DT=div (rho D) grad Y/rho
   // UPDATESEMFORCE:
   // HOFAB=-div rho D grad Y 
   // Ynew=Ynew-(1/rho)(int (HOFAB)dt - dt (LOFAB))
   // call to semdeltaforce with dcomp=slab_step*nstate_SDC+nfluxSEM+1+
   // species_comp
  if ((SDC_outer_sweeps>0)&&
      (SDC_outer_sweeps<ns_time_order)&&
      (divu_outer_sweeps+1==num_divu_outer_sweeps)) {

   if (ns_time_order>=2) {

    if ((viscous_enable_spectral==1)||   // SEM space and time
        (viscous_enable_spectral==3)) {  // SEM time

     for (int ilev=finest_level;ilev>=level;ilev--) {
      NavierStokes& ns_level=getLevel(ilev);
      // calls: SEMDELTAFORCE in GODUNOV_3D.F90
      ns_level.make_SEM_delta_force(project_option_species);
     }
    } else if (viscous_enable_spectral==0) {
     // do nothing
    } else
     amrex::Error("viscous_enable_spectral invalid");

    avgDownALL(State_Type,dencomp,nden,1);
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
  avgDownALL(State_Type,dencomp,nden,1);

  multiphase_project(project_option_species); // species

   // ---------------begin save stable species diffusion 

  if ((ns_time_order>=2)&&
      (ns_time_order<=32)&&
      (divu_outer_sweeps+1==num_divu_outer_sweeps)) {

   get_mm_scomp_solver(
    num_materials_scalar_solve,
    project_option_species,
    state_index,
    scomp,
    ncomp,
    ncomp_check);

   if (ncomp_check!=nsolveMM_thermal)
    amrex::Error("ncomp_check invalid");

     //localMF[PRESPC2_MF] will hold the latest species from the implicit
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
    // LOfab=-div(rho D grad Y)
    // calls: UPDATESEMFORCE in GODUNOV_3D.F90
   if ((viscous_enable_spectral==1)||  // SEM space and time
       (viscous_enable_spectral==3)) { // SEM time
    update_SEM_forcesALL(project_option_species,PRESPC2_MF,
      update_spectralF,update_stableF);
   } else if (viscous_enable_spectral==0) {
    // do nothing
   } else
    amrex::Error("viscous_enable_spectral invalid");

   delete_array(PRESPC2_MF);

  } else if ((ns_time_order==1)||
             ((divu_outer_sweeps>=0)&&
              (divu_outer_sweeps+1<num_divu_outer_sweeps))) {
   // do nothing
  } else
   amrex::Error("ns_time_order or divu_outer_sweeps invalid");

  // ---------------end save stable species diffusion forces

 } // species_comp=0..num_species_var-1

 avgDownALL(State_Type,dencomp,nden,1);

 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);

  // if (FSI_flag(im)==1,2,4)
  //  T(im)=TSOLID
  // else if (FSI_flag(im)=0,3,5,6,7)
  //  T(im)=TSOLID if in the solid.
  int hflag=0;
  ns_level.solid_temperature();

  int project_option_combine=2;  // temperature in veldiffuseALL
  int combine_flag=1;  // GFM -> FVM
  int combine_idx=-1;  // update state variables
  int update_flux=0;
  int interface_cond_avail=1;

  if (convert_temperature==1) {
   combine_flag=1; // GFM -> FVM
  } else if (convert_temperature==0) {
   combine_flag=2;
  } else
   amrex::Error("convert_temperature invalid");

  ns_level.combine_state_variable(
   project_option_combine,
   combine_idx,
   combine_flag,
   hflag,
   update_flux,
   interface_cond_avail); 

  for (int ns=0;ns<num_species_var;ns++) {
   project_option_combine=100+ns; // species

   if (convert_species==1) {
    combine_flag=1;  // GFM -> FVM
   } else if (convert_species==0) {
    combine_flag=2;
   } else
    amrex::Error("convert_species invalid");

   ns_level.combine_state_variable(
    project_option_combine,
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
  int im_tensor=-1;
  int idx_elastic_flux=-1;
  init_gradu_tensorALL(
    im_tensor,
    VISCHEAT_SOURCE_MF,
    do_alloc,
    CELLTENSOR_MF,
    FACETENSOR_MF,
    idx_elastic_flux,
    simple_AMR_BC_flag_viscosity);

  if ((num_materials_viscoelastic>=1)&&
      (num_materials_viscoelastic<=nmat)) {

   for (int im=0;im<nmat;im++) {
    if (ns_is_rigid(im)==0) {
     if ((elastic_time[im]>0.0)&&
         (elastic_viscosity[im]>0.0)) {
      for (int ilev=finest_level;ilev>=level;ilev--) {
       NavierStokes& ns_level=getLevel(ilev);
        // initializes VISCOTEN_MF
       ns_level.make_viscoelastic_tensor(im);
        // VISCHEAT_MF is initialized to zero.
        // VISCHEAT_MF is incremented with heating terms due to viscosity
        // and viscoelastic heating.
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
   } // im=0..nmat-1

  } else {
   amrex::Error("num_materials_viscoelastic invalid");
  }

   // viscosity heating
   // VISCHEAT_SOURCE is the velocity from the viscosity solver.
   // 1. localMF[CONSERVE_FLUXES_MF]=mu(grad U+grad U^T)
   // 2. localMF[VISCHEAT_MF]+=(dt/(rho cv)) mu(grad U+grad U^T) ddot grad U 
  diffusion_heatingALL(VISCHEAT_SOURCE_MF,VISCHEAT_MF);

   // add viscous heating term to T_m  m=1...M
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

 for (int im_assimilate=0;im_assimilate<nmat;im_assimilate++) {
  if (particleLS_flag[im_assimilate]==1) {
   if (ns_is_rigid(im_assimilate)==0) {
    if ((elastic_time[im_assimilate]>=0.0)&&
        (elastic_viscosity[im_assimilate]>=0.0)) {
     if (store_elastic_data[im_assimilate]==1) {
      if (viscoelastic_model[im_assimilate]==2) {
       // particles only appear on the finest level.
       // The flexible substrate is wholly contained on
       // the finest level.
       NavierStokes& ns_finest=getLevel(finest_level);
       ns_finest.assimilate_vel_from_particles(im_assimilate);
      } else if ((viscoelastic_model[im_assimilate]==1)||
 	         (viscoelastic_model[im_assimilate]==0)||
		 (viscoelastic_model[im_assimilate]==3)) { //incremental
       // do nothing
      } else
       amrex::Error("viscoelastic_model[im_assimilate] invalid");
     } else if (store_elastic_data[im_assimilate]==0) {
      // do nothing
     } else
      amrex::Error("store_elastic_data invalid");
    } else
     amrex::Error("elastic_time or elastic_viscosity invalid");
   } else if (ns_is_rigid(im_assimilate)==1) {
    // do nothing
   } else
    amrex::Error("ns_is_rigid(im_assimilate) invalid");
  } else if (particleLS_flag[im_assimilate]==0) {
	  // do nothing
  } else
   amrex::Error("particleLS_flag[im_assimilate] invalid");

 } // im_assimilate=0..nmat-1

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.assimilate_state_data();
 }

 avgDownALL(State_Type,0,
   num_materials_vel*(AMREX_SPACEDIM+1),1);
 avgDownALL(State_Type,dencomp,nden,1);

  // 1. add dt_gradp_over_rho:  u^* <-- u^* + dt grad p/rho
  // 2. save resulting velocity, u^*, to dt_gradp_over_rho
  //    dt_gradp_over_rho <-- u^*
  // 3. delete scratch variables (including CONSERVE_FLUXES_MF)
 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  int alloc_flag=4;
  ns_level.alloc_gradp_over_rho(alloc_flag);
    //    dt_gradp_over_rho <-- u^*
    //    after the projection: dt_gradp_over_rho <--- u^*-u^{n+1}
    //    note: in Jemison, Sussman, Arient:
    //    for k=0...num_divu_outer_sweeps-1
    //     S_t + u^{n+1,(k)} =0
    //     u^* = u^advect + dt div (2 mu D)^*/rho
    //     u^{n+1}=u^* - dt grad p^{n+1,(k)}/rho
    //    end
    //    NOW,
    //    for k=0...num_divu_outer_sweeps-1
    //     S_t + u^{n+1,(k)} =0
    //     u^* = u^advect + dt div (2 mu D)^*/rho - dt (grad p/rho)^{n+1,(k)}
    //     u^* = u^* + dt (grad p/rho)^{n+1,(k)}
    //     u^{n+1}=u^* - dt grad p^{n+1,(k+1)}/rho
    //    end
  alloc_flag=2;
  ns_level.alloc_gradp_over_rho(alloc_flag);

  ns_level.exit_viscous_solver();
 }  // ilev

 delete_array(BOUSSINESQ_TEMP_MF);

 override_enable_spectral(save_enable_spectral);

}   // subroutine veldiffuseALL

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

void NavierStokes::exit_viscous_solver() {

 delete_localMF(CONSERVE_FLUXES_MF,AMREX_SPACEDIM);

 delete_localMF(REGISTER_MARK_MF,1);
 delete_localMF(HOOP_FORCE_MARK_MF,1);
 delete_localMF(NEG_MOM_FORCE_MF,1);
 delete_localMF(THERMAL_FORCE_MF,1);
 delete_localMF(VISCHEAT_MF,1);
 delete_localMF(VISCHEAT_SOURCE_MF,1);

}

// theta new_m = theta new_m + dtheta  m=1..nmat  
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

 int nmat=num_materials;

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 int  num_materials_face=num_materials_scalar_solve;

 int nsolve=1;
 int nsolveMM=nsolve*num_materials_face;

 if (localMF[source_mf]->nComp()!=nsolveMM)
  amrex::Error("diffuse_register invalid ncomp");
 if (localMF[source_mf]->nGrow()<0)
  amrex::Error("diffuse_register invalid ngrow");
 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

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
   amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nComp() bad");
 }

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 MultiFab& LS_new=get_new_data(LS_Type,slab_step+1);

 int nstate=num_materials_vel*(AMREX_SPACEDIM+1)+
  nmat*(num_state_material+ngeom_raw)+1;
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");
 if (LS_new.nComp()!=nmat*(AMREX_SPACEDIM+1))
  amrex::Error("LS_new.nComp()!=nmat*(AMREX_SPACEDIM+1)");

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
  FORT_HEATADVANCE(
    &level,
    &finest_level,
    &cur_time_slab,
    &nmat,
    &nparts,
    &nparts_def,
    im_solid_map_ptr,
    &nsolve, 
    &nsolveMM, 
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
 ns_reconcile_d_num(187);

} // subroutine APPLY_VISCOUS_HEATING

//REGISTER_CURRENT_MF=unew-source_mf
//uface+=INTERP_TO_MAC(REGISTER_CURRENT_MF)
// (filter_vel==0 and face_flag==1)
void NavierStokes::INCREMENT_REGISTERS_ALL(int source_mf,int caller_id) {

 if (level!=0)
  amrex::Error("level invalid INCREMENT_REGISTERS_ALL");
 int finest_level=parent->finestLevel();

  // REGISTER_CURRENT_MF=unew-source_mf
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.INCREMENT_REGISTERS(source_mf,caller_id);
 }

 if (num_materials_vel==1) {
  // do nothing
 } else
  amrex::Error("num_materials_vel invalid");

  // unew^f=unew^f+beta * diffuse_register^{c->f}
  // (filter_vel==0 and face_flag==1)
  // in: INCREMENT_REGISTERS_ALL
 int interp_option=2;
 int project_option=3; // viscosity
 Real beta=1.0;
 Vector<blobclass> blobdata;

  // operation_flag==5 (interp_option==2)
  // if filter_velocity[im]==1 or face_flag==0,
  //  unew^f=INTERP_TO_MAC(unew)
 increment_face_velocityALL(
   interp_option,project_option,
   REGISTER_CURRENT_MF,beta,blobdata);

 delete_array(REGISTER_CURRENT_MF);

} // end subroutine INCREMENT_REGISTERS_ALL

// REGISTER_CURRENT_MF=(unew-source)
void NavierStokes::INCREMENT_REGISTERS(int source_mf,int caller_id) {

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int nsolve=AMREX_SPACEDIM;
 int nsolveMM=nsolve*num_materials_vel;

 int nmat=num_materials;
 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 int nstate=S_new.nComp();
 if (nstate!=num_materials_vel*(AMREX_SPACEDIM+1)+
             nmat*(num_state_material+ngeom_raw)+1)
  amrex::Error("nstate invalid");

 new_localMF(REGISTER_CURRENT_MF,nsolveMM,1,-1);
 push_back_state_register(REGISTER_CURRENT_MF,cur_time_slab,caller_id);

 MultiFab::Subtract(
  *localMF[REGISTER_CURRENT_MF],
  *localMF[source_mf],0,0,nsolveMM,1);

} // INCREMENT_REGISTERS


void NavierStokes::push_back_state_register(int idx_MF,Real time,
  int caller_id) {

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 int project_option_visc=3;
 int nsolve=AMREX_SPACEDIM;
 int nsolveMM=nsolve*num_materials_vel;

 int state_index;
 Vector<int> scomp;
 Vector<int> ncomp;
 int ncomp_check;
 get_mm_scomp_solver(
  num_materials_vel,
  project_option_visc,
  state_index,
  scomp,
  ncomp,
  ncomp_check);
 if (ncomp_check!=nsolveMM)
  amrex::Error("nsolveMM invalid 6613");

 if (state_index!=State_Type)
  amrex::Error("state_index invalid");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");
 if (localMF[idx_MF]->nComp()!=nsolveMM)
  amrex::Error("cell_register invalid ncomp");
 if (localMF[idx_MF]->nGrow()<1)
  amrex::Error("cell_register invalid ngrow");

 MultiFab* snew_mf;
 snew_mf=getState_list(1,scomp,ncomp,time);

 if (snew_mf->nComp()!=nsolveMM)
  amrex::Error("snew_mf->nComp() invalid");

 check_for_NAN(snew_mf,caller_id+1000);

 MultiFab::Copy(*localMF[idx_MF],*snew_mf,0,0,nsolveMM,1);
 delete snew_mf;

} // subroutine push_back_state_register

// stores the current velocity in localMF[idx_MF]
void NavierStokes::SET_STOKES_MARK(int idx_MF,int caller_id) {

 if (level!=0)
  amrex::Error("level invalid SET_STOKES_MARK");
 int finest_level=parent->finestLevel();

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.push_back_state_register(idx_MF,cur_time_slab,caller_id);  
 }

}  // SET_STOKES_MARK

void NavierStokes::prepare_advect_vars(Real time) {

 if (time<0.0)
  amrex::Error("time invalid");

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");
 int nsolve=AMREX_SPACEDIM;
 int nsolveMM=nsolve*num_materials_vel;
 int nsolveMM_FACE_MAC=num_materials_vel;

 new_localMF(ADVECT_REGISTER_MF,nsolveMM,1,-1);
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  getStateMAC_localMF(Umac_Type,ADVECT_REGISTER_FACE_MF+dir,0,dir,
    0,nsolveMM_FACE_MAC,time);
 } // dir
  // advect_register has 1 ghost initialized.
 push_back_state_register(ADVECT_REGISTER_MF,time,201);

} // end subroutine prepare_advect_vars(Real time)

// FUTURE: do the same treatment for advection:
//  variable: dt div (-pI + tau)/rho = unp1-u^advect = dt_non_advect_force
void NavierStokes::alloc_gradp_over_rho(int alloc_flag) {

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");
 int nsolve=AMREX_SPACEDIM;
 int nsolveMM=nsolve*num_materials_vel;
 int nsolveMM_FACE_MAC=num_materials_vel;

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

 if (alloc_flag==1) {
  new_localMF(dt_gradp_over_rho_cell_MF,nsolveMM,1,-1);
  setVal_localMF(dt_gradp_over_rho_cell_MF,0.0,0,nsolveMM,1);
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   new_localMF(dt_gradp_over_rho_face_MF+dir,nsolveMM_FACE_MAC,0,dir);
   setVal_localMF(dt_gradp_over_rho_face_MF+dir,0.0,0,
     nsolveMM_FACE_MAC,0);
  }
 } else if (alloc_flag==0) {
  delete_localMF(dt_gradp_over_rho_cell_MF,1);
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   delete_localMF(dt_gradp_over_rho_face_MF+dir,1);
  }
 } else if (alloc_flag==2) {  // save velocity at end of veldiffuseALL
  init_boundary();
  MultiFab::Copy(*localMF[dt_gradp_over_rho_cell_MF],S_new,0,0,nsolveMM,1);
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   MultiFab& Umac_new=get_new_data(Umac_Type+dir,slab_step+1);
   MultiFab::Copy(*localMF[dt_gradp_over_rho_face_MF+dir],Umac_new,0,0,
     nsolveMM_FACE_MAC,0);
  }
 } else if (alloc_flag==3) { //subtract dt grad p/rho at begin of veldiffuseALL
  MultiFab::Subtract(S_new,*localMF[dt_gradp_over_rho_cell_MF],0,0,nsolveMM,0);
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   MultiFab& Umac_new=get_new_data(Umac_Type+dir,slab_step+1);
   MultiFab::Subtract(Umac_new,*localMF[dt_gradp_over_rho_face_MF+dir],0,0,
     nsolveMM_FACE_MAC,0);
  }
  init_boundary();
 } else if (alloc_flag==4) { //add dt grad p/rho at end of veldiffuseALL
                             //(but before saving the velocity)
  MultiFab::Add(S_new,*localMF[dt_gradp_over_rho_cell_MF],0,0,nsolveMM,0);
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   MultiFab& Umac_new=get_new_data(Umac_Type+dir,slab_step+1);
   MultiFab::Add(Umac_new,*localMF[dt_gradp_over_rho_face_MF+dir],0,0,
     nsolveMM_FACE_MAC,0);
  }
  init_boundary();
 } else if (alloc_flag==5) { //subtract unew from dt_gradp_over_rho at 
                             //end of pressure projection.
  init_boundary();
  MultiFab::Subtract(*localMF[dt_gradp_over_rho_cell_MF],S_new,0,0,nsolveMM,1);
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   MultiFab& Umac_new=get_new_data(Umac_Type+dir,slab_step+1);
   MultiFab::Subtract(*localMF[dt_gradp_over_rho_face_MF+dir],Umac_new,0,0,
     nsolveMM_FACE_MAC,0);
  }
 } else
  amrex::Error("alloc_flag invalid in alloc_gradp_over_rho");

} // end subroutine alloc_gradp_over_rho

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

 int nmat=num_materials;
 int project_option_thermal=2; // temperature project_option
 Vector<int> scomp;  
 Vector<int> ncomp;  
 int state_index;
 int ncomp_check;
 int ngrow=0;

  // DTdt_MF=T_new - T_advect_MF

 get_mm_scomp_solver(
   nmat,
   project_option_thermal,
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

  MultiFab::Subtract(*localMF[DTDt_MF],*localMF[T_advect_MF],0,0,nmat,ngrow);
 } else
  amrex::Error("alloc_flag invalid in alloc_DTDt");

} // end subroutine alloc_DTDt


void NavierStokes::prepare_viscous_solver() {

 int nsolve=AMREX_SPACEDIM;
 int nsolveMM=nsolve*num_materials_vel;

 if (num_materials_vel!=1)
  amrex::Error("num_materials_vel invalid");

 int nsolveMM_FACE=nsolveMM;

 new_localMF(REGISTER_MARK_MF,nsolveMM,1,-1);
 new_localMF(HOOP_FORCE_MARK_MF,nsolveMM,1,-1);
 new_localMF(NEG_MOM_FORCE_MF,nsolveMM,1,-1);
 new_localMF(THERMAL_FORCE_MF,num_materials_scalar_solve,1,-1);
 new_localMF(VISCHEAT_SOURCE_MF,nsolveMM,1,-1);

 new_localMF(VISCHEAT_MF,num_materials_scalar_solve,0,-1);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   // store stress fluxes (-pI+2 mu D) 
  new_localMF(CONSERVE_FLUXES_MF+dir,nsolveMM_FACE,0,dir);
 } // dir

 setVal_localMF(REGISTER_MARK_MF,0.0,0,nsolveMM,1);
 setVal_localMF(VISCHEAT_SOURCE_MF,0.0,0,nsolveMM,1);

 localMF[VISCHEAT_MF]->setVal(0.0,0,num_materials_scalar_solve,0);

 resize_metrics(1);
 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(FACE_VAR_MF,0,750);
 debug_ngrow(VOLUME_MF,1,751);
 debug_ngrow(MASKCOEF_MF,1,752);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  if (localMF[AREA_MF+dir]->boxArray()!=
      localMF[FACE_VAR_MF+dir]->boxArray())
   amrex::Error("boxarrays do not match");
 }

}  // prepare_viscous_solver

// probtype=31 translating circle or sphere
int NavierStokes::is_zalesak() {

 return ( ((probtype==28)||(probtype==29)||(probtype==31)) ? 1 : 0 );
}

// velocity scaled by global_velocity_scale
void NavierStokes::zalesakVEL() {
 
 bool use_tiling=ns_tiling;

 int nmat=num_materials;
 if ((num_materials_vel!=1)&&(num_materials_vel!=nmat))
  amrex::Error("num_materials_vel invalid");

 const Real* dx = geom.CellSize();

 const Box& domain = geom.Domain();
 const int* domlo = domain.loVect();
 const int* domhi = domain.hiVect();

 MultiFab& U_new=get_new_data(State_Type,slab_step+1);
 int scomp_pres=num_materials_vel*AMREX_SPACEDIM;
 U_new.setVal(0.0,scomp_pres,num_materials_vel,1);

 for (int im=0;im<num_materials_vel;im++) {

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

   int velcomp=im*AMREX_SPACEDIM;

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // takes into consideration global_velocity_scale
   FORT_ZALESAKNODE(
    xlo,dx,
    velfab.dataPtr(velcomp),
    domlo,domhi,
    tilelo,tilehi,
    fablo,fabhi,
    &bfact,
    &level,
    ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
    &cur_time_slab);

  } // mfi
} // omp
  ns_reconcile_d_num(188);
 } // im
}

#undef profile_solver

}/* namespace amrex */

