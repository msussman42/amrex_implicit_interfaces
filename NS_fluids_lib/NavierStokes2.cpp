// get rid of autoindent   :setl noai nocin nosi inde=
//
// #include <winstd.H>
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
#include <TECPLOTUTIL_F.H>
#include <MARCHING_TETRA_F.H>
#include <NAVIERSTOKES_F.H>
#include <MACOPERATOR_F.H>
#include <PROB_F.H>
#include <GODUNOV_F.H>
#include <PLIC_F.H>
#include <LEVEL_F.H>
#include <SOLIDFLUID_F.H>
#include <DERIVE_F.H>

#define GEOM_GROW   1
#define bogus_value 1.e20

#define BOGUS (1.0E+6)
#define VOFTOL (1.0E-8)

namespace amrex{

void NavierStokes::Mult_localMF(int idx_dest,int idx_source,
  int scomp,int dcomp,int ncomp,int ngrow) {

 debug_ngrow(idx_dest,ngrow,500); 
 debug_ngrow(idx_source,ngrow,501); 
 if (localMF[idx_dest]->boxArray()!=
     localMF[idx_source]->boxArray())
  amrex::Error("Mult_localMF: boxarrays do not match");

 MultiFab::Multiply(*localMF[idx_dest],*localMF[idx_source],
    scomp,dcomp,ncomp,ngrow);

}  // subroutine Mult_localMF

void NavierStokes::Plus_localMF(int idx_dest,Real val,
  int dcomp,int ncomp,int ngrow) {

 debug_ngrow(idx_dest,ngrow,500); 
 if (localMF[idx_dest]->boxArray()!=grids)
  amrex::Error("Plus_localMF: boxarrays do not match");

 localMF[idx_dest]->plus(val,dcomp,ncomp,ngrow);

}  // subroutine Plus_localMF


void NavierStokes::Copy_localMF(int idx_dest,int idx_source,
  int scomp,int dcomp,int ncomp,int ngrow) {

 debug_ngrow(idx_dest,ngrow,500); 
 debug_ngrow(idx_source,ngrow,501); 
 if (localMF[idx_dest]->boxArray()!=
     localMF[idx_source]->boxArray())
  amrex::Error("Copy_localMF: boxarrays do not match");

 MultiFab::Copy(*localMF[idx_dest],*localMF[idx_source],scomp,
   dcomp,ncomp,ngrow);

}  // copy_LocalMF


void NavierStokes::minus_localMF(int idx_dest,int idx_source,
  int ncomp,int ngrow) {

 debug_ngrow(idx_dest,ngrow,500); 
 debug_ngrow(idx_source,ngrow,501); 
 if (localMF[idx_dest]->boxArray()!=
     localMF[idx_source]->boxArray())
  amrex::Error("minus_localMF: boxarrays do not match");

 localMF[idx_dest]->
   minus(*localMF[idx_source],0,ncomp,ngrow);

}  // minus_LocalMF

//grid_type=-1 ... 5
void NavierStokes::new_localMF(int idx_MF,int ncomp,int ngrow,int grid_type) {

 if (1==0) {
  std::cout << "in new_localMF idx_MF= " << idx_MF << '\n';
  std::cout << "in new_localMF ncomp= " << ncomp << '\n';
  std::cout << "in new_localMF ngrow= " << ngrow << '\n';
  std::cout << "in new_localMF grid_type= " << grid_type << '\n';
  std::cout << "in new_localMF level= " << level << '\n';
 }

 if (localMF_grow[idx_MF]==-1) { 
  // do nothing
 } else {
  std::cout << "in new_localMF idx_MF= " << idx_MF << '\n';
  std::cout << "in new_localMF ncomp= " << ncomp << '\n';
  std::cout << "in new_localMF ngrow= " << ngrow << '\n';
  std::cout << "in new_localMF grid_type= " << grid_type << '\n';
  std::cout << "in new_localMF level= " << level << '\n';
  amrex::Error("localMF_grow invalid");
 }
 if (ngrow<0)
  amrex::Error("ngrow invalid");

 BoxArray edge_boxes(grids);
 if (grid_type==-1) {
  // do nothing
 } else if ((grid_type>=0)&&(grid_type<AMREX_SPACEDIM)) {
  edge_boxes.surroundingNodes(grid_type);
 } else if (grid_type==3) {
  edge_boxes.surroundingNodes(0);
  edge_boxes.surroundingNodes(1);
 } else if ((grid_type==4)&&(AMREX_SPACEDIM==3)) {
  edge_boxes.surroundingNodes(0);
  edge_boxes.surroundingNodes(AMREX_SPACEDIM-1);
 } else if ((grid_type==5)&&(AMREX_SPACEDIM==3)) {
  edge_boxes.surroundingNodes(1);
  edge_boxes.surroundingNodes(AMREX_SPACEDIM-1);
 } else
  amrex::Error("grid_type invalid new_localMF");

 localMF[idx_MF]=new MultiFab(edge_boxes,dmap,ncomp,ngrow,
	MFInfo().SetTag("localMF[idx_MF]"),FArrayBoxFactory());
 localMF[idx_MF]->setVal(0.0,0,ncomp,ngrow);
 localMF_grow[idx_MF]=ngrow;
 debug_ixType(idx_MF,grid_type,idx_MF);

} //new_localMF

//grid_type=-1 ... 5
void NavierStokes::new_localMF_if_not_exist(int idx_MF,int ncomp,
 int ngrow,int grid_type) {

 if (localMF_grow[idx_MF]==-1) {
  new_localMF(idx_MF,ncomp,ngrow,grid_type);
 } else if (localMF_grow[idx_MF]>=0) {
  // do nothing
 } else
  amrex::Error("localMF_grow[idx_MF] invalid");

} //new_localMF_if_not_exist


void NavierStokes::getStateMAC_localMF(int MAC_state_idx,
  int idx_MF,int ngrow,int dir,
  int scomp,int ncomp,Real time) {

 if (localMF_grow[idx_MF]==-1) {
  // do nothing
 } else {
  std::cout << "idx_MF= " << idx_MF << " ngrow = " << ngrow << 
   " dir= " << dir << " time= " << time << '\n';
  amrex::Error("localMF_grow invalid getstateMAC_localMF");
 }
 if (ngrow<0)
  amrex::Error("ngrow invalid");

   // declared in NavierStokes.cpp
 localMF[idx_MF]=getStateMAC(MAC_state_idx,ngrow,dir,scomp,ncomp,time);
 localMF_grow[idx_MF]=ngrow;
} // end getStateMAC_localMF
 
void NavierStokes::getStateDen_localMF(int idx_MF,int ngrow,Real time) {

 if (localMF_grow[idx_MF]==-1) {
  // do nothing
 } else 
  amrex::Error("localMF_grow invalid getStateDen_localMF");

 if (ngrow<0)
  amrex::Error("ngrow invalid");

 localMF[idx_MF]=getStateDen(ngrow,time);
 localMF_grow[idx_MF]=ngrow;
}


void NavierStokes::getStateDist_localMF(int idx_MF,int ngrow,
  Real time,int caller_id) {

 if ((ngrow<0)||(ngrow>ngrow_distance))
  amrex::Error("ngrow invalid");

 if (localMF_grow[idx_MF]==-1) {
  // do nothing
 } else {
  std::cout << "idx_MF,time " << idx_MF << ' ' << time << '\n';
  amrex::Error("localMF_grow invalid getStateDist_localMF ");
 }

 localMF[idx_MF]=getStateDist(ngrow,time,caller_id);
 localMF_grow[idx_MF]=localMF[idx_MF]->nGrow();
} // end subroutine getStateDist_localMF

void NavierStokes::getState_localMF_list(
  int idx_MF,int ngrow,
  int state_index,
  Vector<int> scomp,
  Vector<int> ncomp) {

 if ((ngrow>=0)&&(ngrow<=8)) {
  // do nothing
 } else
  amrex::Error("ngrow invalid in getState_localMF_list");

 if (localMF_grow[idx_MF]==-1) {
  // do nothing
 } else
  amrex::Error("localMF_grow invalid getState_localMF_list");

 if (state_index==State_Type) {
  localMF[idx_MF]=getState_list(ngrow,scomp,ncomp,cur_time_slab);
 } else if (state_index==DIV_Type) {
  if (scomp.size()!=1)
   amrex::Error("scomp.size() invalid");
  localMF[idx_MF]=getStateDIV_DATA(ngrow,scomp[0],ncomp[0],cur_time_slab);
 } else
  amrex::Error("state_index invalid");

 localMF_grow[idx_MF]=ngrow;

} // subroutine getState_localMF_list

void NavierStokes::getState_localMF_listALL(
  int idx_MF,int ngrow,
  int state_index,
  Vector<int> scomp,
  Vector<int> ncomp) {

 if ((ngrow>=0)&&(ngrow<=8)) {
  // do nothing
 } else {
  std::cout << "idx_MF=" << idx_MF << " ngrow= " << ngrow <<
   " state_index= " << state_index << " scomp.size()=" <<
   scomp.size() << '\n';
  amrex::Error("ngrow invalid getState_localMF_listALL");
 }
 if (localMF_grow[idx_MF]==-1) {
  // do nothing
 } else
  amrex::Error("localMF_grow invalid getState_localMF_listALL");

 int finest_level=parent->finestLevel();
 if (level==0) {
  // do nothing
 } else
  amrex::Error("expecting level==0");

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.getState_localMF_list(
   idx_MF,ngrow,state_index,scomp,ncomp);
 } // ilev=finest_level ... level

} // subroutine getState_localMF_listALL


// copy localMF[idx_MF] to s_new
void NavierStokes::putState_localMF_list(
  int idx_MF,
  int state_index,
  Vector<int> scomp,
  Vector<int> ncomp) {

 if (localMF_grow[idx_MF]>=0) {

  if (state_index==State_Type) {
   putState_list(scomp,ncomp,idx_MF);
  } else if (state_index==DIV_Type) {
   if (scomp.size()!=1)
    amrex::Error("scomp.size() invalid");
   putStateDIV_DATA(scomp[0],ncomp[0],idx_MF);
  } else
   amrex::Error("state_index invalid");

 } else
  amrex::Error("localMF_grow[idx_MF] invalid");

} // subroutine putState_localMF_list


// copy localMF[idx_MF] to s_new
void NavierStokes::putState_localMF_listALL(
  int idx_MF,
  int state_index,
  Vector<int> scomp,
  Vector<int> ncomp) {

 int finest_level=parent->finestLevel();
 if (level==0) {
  // do nothing
 } else
  amrex::Error("expecting level==0");

 if (localMF_grow[idx_MF]>=0) {

  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   ns_level.putState_localMF_list(
    idx_MF,state_index,scomp,ncomp);
  } // ilev=finest_level ... level

 } else
  amrex::Error("localMF_grow[idx_MF] invalid");

} // subroutine putState_localMF_listALL



void NavierStokes::getState_localMF(int idx_MF,int ngrow,
  int scomp,int ncomp,Real time) {

 if (localMF_grow[idx_MF]==-1) {
  // do nothing
 } else
  amrex::Error("localMF_grow invalid getState_localMF");

 if (ngrow<0)
  amrex::Error("ngrow invalid");

 localMF[idx_MF]=getState(ngrow,scomp,ncomp,time);
 localMF_grow[idx_MF]=ngrow;

} //getState_localMF


void NavierStokes::getStateTensor_localMF(int idx_MF,int ngrow,
  int scomp,int ncomp,Real time) {

 if (localMF_grow[idx_MF]==-1) {
  // do nothing
 } else
  amrex::Error("localMF_grow invalid getStateTensor_localMF");

 if (ngrow<0)
  amrex::Error("ngrow invalid");

 localMF[idx_MF]=getStateTensor(ngrow,scomp,ncomp,time);
 localMF_grow[idx_MF]=ngrow;

} //subroutine getStateTensor_localMF



void NavierStokes::maskfiner_localMF(int idx_MF,int ngrow,
  Real tag,int clearbdry) {

 delete_localMF_if_exist(idx_MF,1);
 if ((ngrow<0)||(ngrow>ngrow_distance))
  amrex::Error("ngrow invalid");

 localMF[idx_MF]=maskfiner(ngrow,tag,clearbdry);
 localMF_grow[idx_MF]=ngrow;
}  // subroutine maskfiner_localMF

void NavierStokes::getStateVISC_ALL(int idx,int ngrow) {

 if (level!=0)
  amrex::Error("level invalid getStateVISC_ALL");

 int finest_level=parent->finestLevel();
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.getStateVISC(idx,ngrow);
  int scomp=0;
  int ncomp=ns_level.localMF[idx]->nComp();
  ns_level.avgDown_localMF(idx,scomp,ncomp,0);
 }

} // subroutine getStateVISC_ALL 


void NavierStokes::delete_localMF(int idx_MF,int ncomp) {

 for (int scomp=idx_MF;scomp<idx_MF+ncomp;scomp++) {
  if (localMF_grow[scomp]>=0) {
   // do nothing
  } else {
   std::cout << "level= " << level << '\n';
   std::cout << "idx_MF= " << idx_MF << '\n';
   amrex::Error("forgot to allocate the localMF variable before delete");
  }
  delete localMF[scomp];
  localMF_grow[scomp]=-1;
  localMF[scomp]=0;
 }  // scomp
 ParallelDescriptor::Barrier();

} // subroutine delete_localMF

void NavierStokes::delete_localMF_if_exist(int idx_MF,int ncomp) {

 if (localMF_grow[idx_MF]>=0) {
  delete_localMF(idx_MF,ncomp);
 } else if (localMF_grow[idx_MF]==-1) {
  // do nothing
 } else
  amrex::Error("localMF_grow[idx_MF] invalid");

} // subroutine delete_localMF_if_exist


// In the valid region: 
//   mask=tag if not covered by level+1.
// if clear_phys_boundary=0 then
//  mask=tag outside the domain.
//  mask=tag at coarse-fine ghost cell.
// if clear_phys_boundary=1 then
//  mask=1-tag outside the domain.
//  mask=tag at coarse-fine ghost cell.
// if clear_phys_boundary=2 then
//  mask=1-tag outside the domain.
//  mask=1-tag at coarse-fine ghost cell.
//  mask=tag at fine-fine uncovered ghost cell
// if clear_phys_boundary=3 then
//  mask=tag in valid region
//  mask=1-tag outside the domain.
//  mask=1-tag at coarse-fine ghost cell.
//  mask=tag at fine-fine ghost cell
MultiFab* NavierStokes::maskfiner(int ngrow,Real tag,int clear_phys_boundary) {

 if ((ngrow<0)||(ngrow>ngrow_distance))
  amrex::Error("ngrow invalid");
 
 int finest_level = parent->finestLevel();
 MultiFab* mask=new MultiFab(grids,dmap,1,ngrow,
  MFInfo().SetTag("mask"),FArrayBoxFactory());

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(mask->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*mask,false); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  const Box& tilegrid = mfi.tilebox();

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  FArrayBox& fab = (*mask)[mfi];
  if (clear_phys_boundary==0) {
   fab.setVal(tag);
  } else if (clear_phys_boundary==1) {
   fab.setVal(1.0-tag);
   Box d_box=geom.Domain();
   d_box &= fab.box();
   fab.setVal(tag,d_box,0);     // fab=1-tag outside domain
  } else if ((clear_phys_boundary==2)||
             (clear_phys_boundary==3)) {
   fab.setVal(1.0-tag);
   Box d_box=geom.Domain();
   d_box &= grids[mfi.index()];
   fab.setVal(tag,d_box,0); // fab=1-tag outside domain and coarse/fine border
  } else
   amrex::Error("clear_phys_boundary invalid");
 } // mfi
} //omp
 ns_reconcile_d_num(126);

 if ((clear_phys_boundary==0)||
     (clear_phys_boundary==1)||
     (clear_phys_boundary==2)) {

  if (level<finest_level) {
   NavierStokes& ns_fine=getLevel(level+1);
   const BoxArray& f_box=ns_fine.grids;

   if (thread_class::nthreads<1)
    amrex::Error("thread_class::nthreads invalid");
   thread_class::init_d_numPts(mask->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
   for (MFIter mfi(*mask,false); mfi.isValid(); ++mfi) {
    BL_ASSERT(grids[mfi.index()] == mfi.validbox());
    const Box& tilegrid = mfi.tilebox();

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    FArrayBox& fab = (*mask)[mfi];
    for (int j = 0; j < f_box.size(); j++) {
     Box c_box = amrex::coarsen(f_box[j],2);
     c_box &= grids[mfi.index()];
     if (c_box.ok())
      fab.setVal(1.0-tag,c_box,0);
    } // j
   } // mfi
} //omp
   ns_reconcile_d_num(127);
  } // level<finest_level
 } else if (clear_phys_boundary==3) {
  // do nothing
 } else
  amrex::Error("clear_phys_boundary invalid");

// if clear_phys_boundary!=3,
// mask=1-tag at a ghost cell that is not outside the domain
// and covered by a finer cell.
 mask->FillBoundary(geom.periodicity());  
 return mask;
}

// if spectral_override==0, then always low order average down.
void NavierStokes::avgDown_localMF(int idxMF,int scomp,int ncomp,
  int spectral_override) {

 int finest_level=parent->finestLevel();
 if (level<finest_level) {
  NavierStokes& ns_fine=getLevel(level+1);
  debug_ngrow(idxMF,0,2000);
  ns_fine.debug_ngrow(idxMF,0,2500);
  MultiFab& S_crse=*localMF[idxMF];
  MultiFab& S_fine=*ns_fine.localMF[idxMF];
  avgDown(S_crse,S_fine,scomp,ncomp,spectral_override);
 }

} // avgDown_localMF


void NavierStokes::avgDown_tag_localMF(int idxMF) {

 int finest_level=parent->finestLevel();
 if (level<finest_level) {
  NavierStokes& ns_fine=getLevel(level+1);
  debug_ngrow(idxMF,0,2000);
  ns_fine.debug_ngrow(idxMF,0,2500);
  MultiFab& S_crse=*localMF[idxMF];
  MultiFab& S_fine=*ns_fine.localMF[idxMF];
  level_avgDown_tag(S_crse,S_fine);
 }

} // subroutine avgDown_tag_localMF


void NavierStokes::avgDownBURNING_localMF(
  int burnvel_MF,int TSAT_MF) {

 int finest_level=parent->finestLevel();
 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;
 int ncomp_per_burning=AMREX_SPACEDIM;
 int ncomp_per_tsat=2;
 int nburning=nten*(ncomp_per_burning+1);
 int ntsat=nten*(ncomp_per_tsat+1);

 if (level<finest_level) {
  NavierStokes& ns_fine=getLevel(level+1);

  debug_ngrow(burnvel_MF,0,2000);
  ns_fine.debug_ngrow(burnvel_MF,0,2500);
  MultiFab& S_crseB=*localMF[burnvel_MF];
  MultiFab& S_fineB=*ns_fine.localMF[burnvel_MF];
  if ((S_crseB.nComp()==nburning)&&
      (S_fineB.nComp()==nburning)) {
   int velflag=1;
   level_avgDownBURNING(S_crseB,S_fineB,velflag);
  } else
   amrex::Error("S_crseB or S_fineB invalid nComp");

  debug_ngrow(TSAT_MF,0,2000);
  ns_fine.debug_ngrow(TSAT_MF,0,2500);
  MultiFab& S_crseT=*localMF[TSAT_MF];
  MultiFab& S_fineT=*ns_fine.localMF[TSAT_MF];
  if ((S_crseT.nComp()==ntsat)&&
      (S_fineT.nComp()==ntsat)) {
   int velflag=0;
   level_avgDownBURNING(S_crseT,S_fineT,velflag);
  } else
   amrex::Error("S_crseT or S_fineT invalid nComp");

 } // level<finest_level

} // subroutine avgDownBURNING_localMF


void NavierStokes::avgDownCURV_localMF(int idxMF) {

 int finest_level=parent->finestLevel();
 if (level<finest_level) {
  NavierStokes& ns_fine=getLevel(level+1);
  debug_ngrow(idxMF,0,2000);
  ns_fine.debug_ngrow(idxMF,0,2500);
  MultiFab& S_crse=*localMF[idxMF];
  MultiFab& S_fine=*ns_fine.localMF[idxMF];
  level_avgDownCURV(S_crse,S_fine);
 }

} // subroutine avgDownCURV_localMF


// flux variables: average down in the tangential direction 
// to the box face, copy in
// the normal direction.  Since the blocking factor is >=2, it is
// impossible to have a grid box with size of 1 cell width.
void NavierStokes::avgDown_and_Copy_localMF(
  int idx_den_MF,
  int idx_vel_MF,
  int idx_flux_MF,
  int operation_flag) {

 int finest_level=parent->finestLevel();
 int nmat=num_materials;

 int ncomp_den=0;
 int ncomp_vel=0;
 int ncomp_flux=0;
 int scomp_flux=0;
 int ncomp_flux_use=0;

 if (operation_flag==7) {  // advection 
  ncomp_den=nmat*num_state_material;
  ncomp_vel=AMREX_SPACEDIM;
  ncomp_flux=nfluxSEM;
  scomp_flux=0;
  ncomp_flux_use=nfluxSEM;
 } else if (operation_flag==1) { // P^cell->MAC
  ncomp_den=1;
  ncomp_vel=1;
  ncomp_flux=1;
  scomp_flux=0;
  ncomp_flux_use=1;
 } else if ((operation_flag==3)|| // u cell to MAC
            (operation_flag==5)|| // uMAC=uMAC+beta * diff_reg
	    (operation_flag==10)||
	    (operation_flag==11)) {
  ncomp_den=AMREX_SPACEDIM;
  ncomp_vel=AMREX_SPACEDIM;
  ncomp_flux=1;
  scomp_flux=0;
  ncomp_flux_use=1;
 } else if (operation_flag==9) {
  ncomp_den=nmat*num_state_material;
  ncomp_vel=nmat*num_state_material;
  ncomp_flux=1;
  scomp_flux=0;
  ncomp_flux_use=1;
 } else if (operation_flag==0) {
  ncomp_den=1;
  ncomp_vel=1;
  ncomp_flux=1;
  scomp_flux=0;
  ncomp_flux_use=1;
 } else if (operation_flag==8) {  // viscous flux
  ncomp_den=AMREX_SPACEDIM;
  ncomp_vel=AMREX_SPACEDIM;
  ncomp_flux=AMREX_SPACEDIM;
  scomp_flux=0;
  ncomp_flux_use=AMREX_SPACEDIM;
 } else
  amrex::Error("operation_flag invalid");

 if (localMF[idx_den_MF]->nGrow()!=1)
  amrex::Error("localMF[idx_den_MF]->nGrow()!=1");
 if (localMF[idx_vel_MF]->nGrow()!=1)
  amrex::Error("localMF[idx_vel_MF]->nGrow()!=1");
 if (localMF[idx_den_MF]->nComp()!=ncomp_den)
  amrex::Error("localMF[idx_den_MF]->nComp()!=ncomp_den");
 if (localMF[idx_vel_MF]->nComp()!=ncomp_vel)
  amrex::Error("localMF[idx_vel_MF]->nComp()!=ncomp_vel");
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  if (localMF[idx_flux_MF+dir]->nGrow()!=0)
   amrex::Error("localMF[idx_flux_MF]->nGrow()!=0");
  if (localMF[idx_flux_MF+dir]->nComp()!=ncomp_flux)
   amrex::Error("localMF[idx_flux_MF]->nComp()!=ncomp_flux");
 } // dir=0..sdim-1

 if ((level>=0)&&(level<finest_level)) {
  int f_level=level+1;
  NavierStokes& fine_lev=getLevel(f_level);
  const BoxArray& fgrids=fine_lev.grids;
  const DistributionMapping& fdmap=fine_lev.dmap;

  if (grids!=localMF[idx_den_MF]->boxArray())
   amrex::Error("grids!=localMF[idx_den_MF]->boxArray()");
  if (fgrids!=fine_lev.localMF[idx_den_MF]->boxArray())
   amrex::Error("fgrids!=fine_lev.localMF[idx_den_MF]->boxArray()");
  if (grids!=localMF[idx_vel_MF]->boxArray())
   amrex::Error("grids!=localMF[idx_vel_MF]->boxArray()");
  if (fgrids!=fine_lev.localMF[idx_vel_MF]->boxArray())
   amrex::Error("fgrids!=fine_lev.localMF[idx_vel_MF]->boxArray()");
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   if (localMF[idx_flux_MF+dir]->nComp()!=
       fine_lev.localMF[idx_flux_MF+dir]->nComp())
    amrex::Error("flux ncomp do not match");
  }

  BoxArray crse_S_fine_BA(fgrids.size());
  for (int i = 0; i < fgrids.size(); ++i) {
   crse_S_fine_BA.set(i,amrex::coarsen(fgrids[i],2));
  }
  resize_metrics(1);
  debug_ngrow(VOLUME_MF,0,700);
  fine_lev.resize_metrics(1);
  fine_lev.debug_ngrow(VOLUME_MF,0,700);

  debug_ngrow(LEVELPC_MF,1,37);
  if (localMF[LEVELPC_MF]->nComp()!=nmat*(1+AMREX_SPACEDIM))
   amrex::Error("localMF[LEVELPC_MF]->nComp()!=nmat*(1+AMREX_SPACEDIM)");

  fine_lev.debug_ngrow(LEVELPC_MF,1,37);
  if (fine_lev.localMF[LEVELPC_MF]->nComp()!=nmat*(1+AMREX_SPACEDIM))
   amrex::Error("fine_lev.localMF[LEVELPC_MF]->nComp() invalid");

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

   MultiFab& S_crse_MAC=*localMF[idx_flux_MF+dir];
   MultiFab& S_fine_MAC=*fine_lev.localMF[idx_flux_MF+dir];

   if (fine_lev.localMF[AREA_MF+dir]->boxArray()!=S_fine_MAC.boxArray()) 
    amrex::Error("fine idx_flux_MF+dir invalid boxArray");
   if (localMF[AREA_MF+dir]->boxArray()!=S_crse_MAC.boxArray())
    amrex::Error("coarseidx_flux_MF+dir invalid boxArray");

   const BoxArray& fgridsMAC=S_fine_MAC.boxArray();

   if (fgrids.size()!=fgridsMAC.size())
    amrex::Error("fgrids.size()!=fgridsMAC.size()");

   BoxArray crse_S_fine_BA_MAC(fgridsMAC.size());
   for (int i = 0; i < fgridsMAC.size(); ++i) {
    crse_S_fine_BA_MAC.set(i,amrex::coarsen(fgridsMAC[i],2));
   }
   DistributionMapping crse_dmap=fdmap;
   MultiFab crse_S_fine_MAC(crse_S_fine_BA_MAC,crse_dmap,ncomp_flux_use,0,
     MFInfo().SetTag("crse_S_fine_MAC"),FArrayBoxFactory());
   crse_S_fine_MAC.setVal(1.0e+40);

   ParallelDescriptor::Barrier();

   const Real* dx = geom.CellSize();
   const Real* dxf = fine_lev.geom.CellSize();
   const Real* prob_lo   = geom.ProbLo();
 
   if (thread_class::nthreads<1)
    amrex::Error("thread_class::nthreads invalid");
   thread_class::init_d_numPts(
      fine_lev.localMF[idx_den_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
   for (MFIter mfi(*fine_lev.localMF[idx_den_MF],false); 
		   mfi.isValid(); ++mfi) {
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

    FArrayBox& fine_fab=S_fine_MAC[gridno];
    const Box& fgridMAC=fine_fab.box();
    const int* flo=fgridMAC.loVect();
    const int* fhi=fgridMAC.hiVect();
    const Real* f_dat=fine_fab.dataPtr();

    FArrayBox& fine_den_fab=(*fine_lev.localMF[idx_den_MF])[gridno];
    const Box& fgrid_den=fine_den_fab.box();
    const int* fdenlo=fgrid_den.loVect();
    const int* fdenhi=fgrid_den.hiVect();
    const Real* f_den_dat=fine_den_fab.dataPtr();

    FArrayBox& fine_vel_fab=(*fine_lev.localMF[idx_vel_MF])[gridno];
    const Box& fgrid_vel=fine_vel_fab.box();
    if (fgrid_vel!=fgrid_den)
     amrex::Error("fgrid_vel!=fgrid_den");
    const Real* f_vel_dat=fine_vel_fab.dataPtr();

    FArrayBox& coarse_fab=crse_S_fine_MAC[gridno];
    const Box& cgridMAC=coarse_fab.box();
    const int* clo=cgridMAC.loVect();
    const int* chi=cgridMAC.hiVect();
    const Real* c_dat=coarse_fab.dataPtr();

    int bfact_c=parent->Space_blockingFactor(level);
    int bfact_f=parent->Space_blockingFactor(f_level);

    const Box& fine_fabgrid = fine_lev.grids[gridno];
    const int* fine_fablo=fine_fabgrid.loVect();
    const int* fine_fabhi=fine_fabgrid.hiVect();

    FArrayBox& maskfab=(*fine_lev.localMF[MASKSEM_MF])[mfi];

    FArrayBox& fine_LS_fab=(*fine_lev.localMF[LEVELPC_MF])[mfi];

     // in: NAVIERSTOKES_3D.F90
    FORT_AVGDOWN_COPY( 
     &enable_spectral,
     &finest_level,
     &operation_flag,
     &dir,
     prob_lo,
     dxf,
     &level,&f_level,
     &bfact_c,&bfact_f,     
     xlo_fine,dx,
     &ncomp_flux_use,
     &ncomp_den,
     &ncomp_vel,
     c_dat,ARLIM(clo),ARLIM(chi),
     f_dat,ARLIM(flo),ARLIM(fhi),
     f_den_dat,ARLIM(fdenlo),ARLIM(fdenhi),
     f_vel_dat,ARLIM(fdenlo),ARLIM(fdenhi),
     maskfab.dataPtr(),
     ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
     fine_LS_fab.dataPtr(),
     ARLIM(fine_LS_fab.loVect()),ARLIM(fine_LS_fab.hiVect()),
     ovlo,ovhi,  // ovlo,ovhi = coarsened(fine_fablo,fine_fabhi)
     fine_fablo,fine_fabhi);

   }// mfi
} //omp
   ns_reconcile_d_num(128);

   S_crse_MAC.copy(crse_S_fine_MAC,0,scomp_flux,ncomp_flux_use);
   ParallelDescriptor::Barrier();

   const Box& domain = geom.Domain();
   if (geom.isPeriodic(dir)) {
    IntVect pshift=IntVect::TheZeroVector();
    pshift[dir]=domain.length(dir);
    crse_S_fine_MAC.shift(pshift);

    ParallelDescriptor::Barrier();
    S_crse_MAC.copy(crse_S_fine_MAC,0,scomp_flux,ncomp_flux_use);
    ParallelDescriptor::Barrier();

    pshift[dir]=-2*domain.length(dir);
    crse_S_fine_MAC.shift(pshift);

    S_crse_MAC.copy(crse_S_fine_MAC,0,scomp_flux,ncomp_flux_use);
    ParallelDescriptor::Barrier();
   }  // isPeriodic(dir)

  } // dir=0..sdim-1

 } else if (level==finest_level) {
  // do nothing
 } else
  amrex::Error("level invalid10");

} // subroutine avgDown_and_Copy_localMF


// flux variables: interp in the tangential direction 
// to the box face, copy in
// the normal direction.  
void NavierStokes::interp_and_Copy_localMF(
  int idx_den_MF,
  int idx_vel_MF,
  int idx_flux_MF,
  int operation_flag) {

 int finest_level=parent->finestLevel();
 int nmat=num_materials;

 int ncomp_den=0;
 int ncomp_vel=0;
 int ncomp_flux=0;
 int scomp_flux=0;
 int ncomp_flux_use=0;


 if (operation_flag==8) {  // viscosity
  ncomp_den=AMREX_SPACEDIM;
  ncomp_vel=AMREX_SPACEDIM;
  ncomp_flux=AMREX_SPACEDIM;
  scomp_flux=0;
  ncomp_flux_use=AMREX_SPACEDIM;
 } else if (operation_flag==7) {  // advection 
  ncomp_den=nmat*num_state_material;
  ncomp_vel=AMREX_SPACEDIM;
  ncomp_flux=nfluxSEM;
  scomp_flux=0;
  ncomp_flux_use=nfluxSEM;
 } else if (operation_flag==1) { // P^cell->MAC
  ncomp_den=1;
  ncomp_vel=1;
  ncomp_flux=1;
  scomp_flux=0;
  ncomp_flux_use=1;
 } else if ((operation_flag==3)|| // u cell to MAC
            (operation_flag==5)|| // uMAC=uMAC+beta * diff_reg
	    (operation_flag==10)||
	    (operation_flag==11)) {
  ncomp_den=AMREX_SPACEDIM;
  ncomp_vel=AMREX_SPACEDIM;
  ncomp_flux=1;
  scomp_flux=0;
  ncomp_flux_use=1;
 } else if (operation_flag==9) {
  ncomp_den=nmat*num_state_material;
  ncomp_vel=nmat*num_state_material;
  ncomp_flux=1;
  scomp_flux=0;
  ncomp_flux_use=1;
 } else if (operation_flag==0) {
  ncomp_den=1;
  ncomp_vel=1;
  ncomp_flux=1;
  scomp_flux=0;
  ncomp_flux_use=1;
 } else
  amrex::Error("operation_flag invalid");

 if (localMF[idx_den_MF]->nGrow()!=1)
  amrex::Error("localMF[idx_den_MF]->nGrow()!=1");
 if (localMF[idx_vel_MF]->nGrow()!=1)
  amrex::Error("localMF[idx_vel_MF]->nGrow()!=1");
 if (localMF[idx_den_MF]->nComp()!=ncomp_den)
  amrex::Error("localMF[idx_den_MF]->nComp()!=ncomp_den");
 if (localMF[idx_vel_MF]->nComp()!=ncomp_vel)
  amrex::Error("localMF[idx_vel_MF]->nComp()!=ncomp_vel");
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  if (localMF[idx_flux_MF+dir]->nGrow()!=0)
   amrex::Error("localMF[idx_flux_MF]->nGrow()!=0");
  if (localMF[idx_flux_MF+dir]->nComp()!=ncomp_flux)
   amrex::Error("localMF[idx_flux_MF]->nComp()!=ncomp_flux");
 } // dir=0..sdim-1

 if ((level>=1)&&(level<=finest_level)) {
  int c_level=level-1;
  NavierStokes& coarse_lev=getLevel(c_level);
  const BoxArray& cgrids=coarse_lev.grids;
  const DistributionMapping& fdmap=dmap;

  if (grids!=localMF[idx_den_MF]->boxArray())
   amrex::Error("grids!=localMF[idx_den_MF]->boxArray()");
  if (cgrids!=coarse_lev.localMF[idx_den_MF]->boxArray())
   amrex::Error("cgrids!=coarse_lev.localMF[idx_den_MF]->boxArray()");
  if (grids!=localMF[idx_vel_MF]->boxArray())
   amrex::Error("grids!=localMF[idx_vel_MF]->boxArray()");
  if (cgrids!=coarse_lev.localMF[idx_vel_MF]->boxArray())
   amrex::Error("cgrids!=coarse_lev.localMF[idx_vel_MF]->boxArray()");

  BoxArray crse_S_fine_BA(grids.size());
  for (int i = 0; i < grids.size(); ++i) {
   crse_S_fine_BA.set(i,amrex::coarsen(grids[i],2));
  }
  resize_metrics(1);
  debug_ngrow(VOLUME_MF,0,700);
  coarse_lev.resize_metrics(1);
  coarse_lev.debug_ngrow(VOLUME_MF,0,700);

  resize_mask_nbr(1);
   // (1) =1 interior  =1 fine-fine ghost in domain  =0 otherwise
   // (2) =1 interior  =0 otherwise
  debug_ngrow(MASK_NBR_MF,1,253); 

  debug_ngrow(LEVELPC_MF,1,37);
  if (localMF[LEVELPC_MF]->nComp()!=nmat*(1+AMREX_SPACEDIM))
   amrex::Error("localMF[LEVELPC_MF]->nComp()!=nmat*(1+AMREX_SPACEDIM)");

  coarse_lev.debug_ngrow(LEVELPC_MF,1,37);
  if (coarse_lev.localMF[LEVELPC_MF]->nComp()!=nmat*(1+AMREX_SPACEDIM))
   amrex::Error("coarse_lev.localMF[LEVELPC_MF]->nComp() invalid");

  MultiFab* coarse_den_fine;
  MultiFab* coarse_vel_fine;
  MultiFab* coarse_mask_sem_fine;
  MultiFab* coarse_LS_fine;

  DistributionMapping crse_dmap=fdmap;
  coarse_mask_sem_fine=new MultiFab(crse_S_fine_BA,crse_dmap,1,1,
	MFInfo().SetTag("coarse_mask_sem_fine"),FArrayBoxFactory());
  coarse_LS_fine=new MultiFab(crse_S_fine_BA,crse_dmap,nmat,1,
	MFInfo().SetTag("coarse_LS_fine"),FArrayBoxFactory());
   // FabArray.H     
   // scomp,dcomp,ncomp,s_nghost,d_nghost
  coarse_mask_sem_fine->copy(*coarse_lev.localMF[MASKSEM_MF],0,0,
    1,1,1,geom.periodicity());
  coarse_LS_fine->copy(*coarse_lev.localMF[LEVELPC_MF],0,0,
    nmat,1,1,geom.periodicity());

  coarse_den_fine=new MultiFab(crse_S_fine_BA,crse_dmap,ncomp_den,1,
		  MFInfo().SetTag("coarse_den_fine"),FArrayBoxFactory());
  coarse_vel_fine=new MultiFab(crse_S_fine_BA,crse_dmap,ncomp_vel,1,
		  MFInfo().SetTag("coarse_vel_fine"),FArrayBoxFactory());
    // FabArray.H     
    // scomp,dcomp,ncomp,s_nghost,d_nghost
  coarse_den_fine->copy(*coarse_lev.localMF[idx_den_MF],0,0,
    ncomp_den,1,1,geom.periodicity());
  coarse_vel_fine->copy(*coarse_lev.localMF[idx_vel_MF],0,0,
    ncomp_vel,1,1,geom.periodicity());

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

   MultiFab& S_fine_MAC=*localMF[idx_flux_MF+dir];
   const BoxArray& fgridsMAC=S_fine_MAC.boxArray();

   if (localMF[AREA_MF+dir]->boxArray()!=fgridsMAC) 
    amrex::Error("fine idx_flux_MF+dir invalid boxArray");

   if (grids.size()!=fgridsMAC.size())
    amrex::Error("grids.size()!=fgridsMAC.size()");

   ParallelDescriptor::Barrier();

   const Real* dx = geom.CellSize();
   const Real* dxc = coarse_lev.geom.CellSize();
   const Real* prob_lo   = geom.ProbLo();
 
   if (thread_class::nthreads<1)
    amrex::Error("thread_class::nthreads invalid");
   thread_class::init_d_numPts(localMF[idx_den_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
   for (MFIter mfi(*localMF[idx_den_MF],false); mfi.isValid(); ++mfi) {
    BL_ASSERT(grids[mfi.index()] == mfi.validbox());
    const int gridno = mfi.index();

    const Box& tilegrid = mfi.tilebox();
    const Box& fabgrid = grids[gridno];
    const int* tilelo=tilegrid.loVect();
    const int* tilehi=tilegrid.hiVect();
    const int* fablo=fabgrid.loVect();
    const int* fabhi=fabgrid.hiVect();

    const Real* xlo = grid_loc[gridno].lo();

    const Box& ovgrid = crse_S_fine_BA[gridno];
    const int* ovlo=ovgrid.loVect();
    const int* ovhi=ovgrid.hiVect();

    FArrayBox& fine_fab=S_fine_MAC[gridno];
    const Box& fgridMAC=fine_fab.box();
    const int* flo=fgridMAC.loVect();
    const int* fhi=fgridMAC.hiVect();
    const Real* f_dat=fine_fab.dataPtr(scomp_flux);

    FArrayBox& coarse_den_fab=(*coarse_den_fine)[gridno];
    const Box& cgrid_den=coarse_den_fab.box();
    const int* cdenlo=cgrid_den.loVect();
    const int* cdenhi=cgrid_den.hiVect();
    const Real* c_den_dat=coarse_den_fab.dataPtr();

    FArrayBox& coarse_vel_fab=(*coarse_vel_fine)[gridno];
    const Box& cgrid_vel=coarse_vel_fab.box();
    if (cgrid_vel!=cgrid_den)
     amrex::Error("cgrid_vel!=cgrid_den");
    const Real* c_vel_dat=coarse_vel_fab.dataPtr();

    int bfact_f=parent->Space_blockingFactor(level);
    int bfact_c=parent->Space_blockingFactor(c_level);

    FArrayBox& cLS_fab=(*coarse_LS_fine)[mfi];

    FArrayBox& maskfab=(*localMF[MASKSEM_MF])[mfi];
    FArrayBox& cmaskfab=(*coarse_mask_sem_fine)[mfi];
    FArrayBox& mnbrfab=(*localMF[MASK_NBR_MF])[mfi];
    Vector<int> velbc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     // in: NAVIERSTOKES_3D.F90
    FORT_INTERP_COPY( 
     &enable_spectral,
     dxc,dx,
     &finest_level,
     &operation_flag,
     tilelo,tilehi,
     fablo,fabhi,
     &dir,
     prob_lo,
     &c_level,
     &level,
     &bfact_c,
     &bfact_f,     
     xlo,
     &ncomp_flux_use,
     &ncomp_den,
     &ncomp_vel,
     f_dat,ARLIM(flo),ARLIM(fhi),
     c_den_dat,ARLIM(cdenlo),ARLIM(cdenhi),
     c_vel_dat,ARLIM(cdenlo),ARLIM(cdenhi),
     mnbrfab.dataPtr(),
     ARLIM(mnbrfab.loVect()),ARLIM(mnbrfab.hiVect()),
     maskfab.dataPtr(),
     ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
     cmaskfab.dataPtr(),
     ARLIM(cmaskfab.loVect()),ARLIM(cmaskfab.hiVect()),
     cLS_fab.dataPtr(),
     ARLIM(cLS_fab.loVect()),ARLIM(cLS_fab.hiVect()),
     velbc.dataPtr(),
     ovlo,ovhi);  // ovlo,ovhi = crse_S_fine_BA[gridno]

   }// mfi
} //omp
   ns_reconcile_d_num(129);
  } // dir=0..sdim-1

  delete coarse_den_fine;
  delete coarse_vel_fine;

  delete coarse_mask_sem_fine;
  delete coarse_LS_fine;

 } else if (level==0) {
  // do nothing
 } else
  amrex::Error("level invalid11");

} // subroutine interp_and_Copy_localMF

void NavierStokes::sync_flux_var(int dir,int flux_MF,int ncomp_flux) {

 int finest_level=parent->finestLevel();

 if ((dir>=0)&&(dir<AMREX_SPACEDIM)) {
  
  if (localMF[AREA_MF+dir]->boxArray()==
      localMF[flux_MF+dir]->boxArray()) {
   if (localMF[flux_MF+dir]->nComp()==ncomp_flux) {

    debug_ngrow(flux_MF+dir,0,700);
    resize_mask_nbr(1);
     // (1) =1 interior  =1 fine-fine ghost in domain  =0 otherwise
     // (2) =1 interior  =0 otherwise
    debug_ngrow(MASK_NBR_MF,1,253); 
     // mask=tag if not covered by level+1 or outside the domain.
    resize_maskfiner(1,MASKCOEF_MF);
    MultiFab* flux_hold_mf=new MultiFab(grids,dmap,ncomp_flux,1,
		 MFInfo().SetTag("flux_hold_mf"),FArrayBoxFactory());
    flux_hold_mf->setVal(1.0e+30);

    for (int sync_iter=0;sync_iter<=1;sync_iter++) {

     if (thread_class::nthreads<1)
      amrex::Error("thread_class::nthreads invalid");
     thread_class::init_d_numPts(localMF[MASKCOEF_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
     for (MFIter mfi(*localMF[MASKCOEF_MF],false); mfi.isValid(); ++mfi) {
      BL_ASSERT(grids[mfi.index()] == mfi.validbox());
      const int gridno = mfi.index();

      const Box& tilegrid = mfi.tilebox();
      const Box& fabgrid = grids[gridno];
      const int* tilelo=tilegrid.loVect();
      const int* tilehi=tilegrid.hiVect();
      const int* fablo=fabgrid.loVect();
      const int* fabhi=fabgrid.hiVect();

      FArrayBox& fluxtarg=(*localMF[flux_MF+dir])[mfi];
      FArrayBox& fluxhold=(*flux_hold_mf)[mfi];

      FArrayBox& maskcov=(*localMF[MASKCOEF_MF])[mfi];
      FArrayBox& masknbr=(*localMF[MASK_NBR_MF])[mfi];
      Vector<int> presbc=getBCArray(State_Type,gridno,AMREX_SPACEDIM,1);

      int tid_current=ns_thread();
      if ((tid_current<0)||(tid_current>=thread_class::nthreads))
       amrex::Error("tid_current invalid");
      thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

      FORT_FILLBDRY_FLUX( 
       &sync_iter,
       &level,
       &finest_level,
       tilelo,tilehi,
       fablo,fabhi,
       &dir,
       &ncomp_flux,
       fluxtarg.dataPtr(),
       ARLIM(fluxtarg.loVect()),ARLIM(fluxtarg.hiVect()),
       fluxhold.dataPtr(),
       ARLIM(fluxhold.loVect()),ARLIM(fluxhold.hiVect()),
       maskcov.dataPtr(),
       ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
       masknbr.dataPtr(),
       ARLIM(masknbr.loVect()),ARLIM(masknbr.hiVect()),
       presbc.dataPtr());

     }// mfi
} //omp
     ns_reconcile_d_num(130);

     if (sync_iter==0) {
      flux_hold_mf->FillBoundary(geom.periodicity()); 
     } else if (sync_iter==1) {
      // do nothing
     } else
      amrex::Error("sync_iter invalid");

    } // sync_iter=0..1

    delete flux_hold_mf;

   } else
    amrex::Error("localMF[flux_MF+dir]->nComp() invalid");
  } else
   amrex::Error("localMF[flux_MF+dir]->boxArray() invalid");
 } else
  amrex::Error("dir invalid");
    
} // subroutine sync_flux_var

void NavierStokes::interp_flux_localMF(
  int coarse_flux_MF,
  int fine_flux_MF) {

 int finest_level=parent->finestLevel();

 int ncomp_flux=nfluxSEM;

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  if (localMF[AREA_MF+dir]->boxArray()!=
      localMF[coarse_flux_MF+dir]->boxArray())
   amrex::Error("coarse_flux_MF invalid box array");
  if (localMF[coarse_flux_MF+dir]->nGrow()!=0)
   amrex::Error("localMF[coarse_flux_MF]->nGrow()!=0");
  if (localMF[coarse_flux_MF+dir]->nComp()!=ncomp_flux)
   amrex::Error("localMF[coarse_flux_MF]->nComp()!=ncomp_flux");
  if (localMF[fine_flux_MF+dir]->nGrow()!=0)
   amrex::Error("localMF[fine_flux_MF]->nGrow()!=0");
  if (localMF[fine_flux_MF+dir]->nComp()!=ncomp_flux)
   amrex::Error("localMF[fine_flux_MF]->nComp()!=ncomp_flux");
  if (localMF[AREA_MF+dir]->boxArray()!=
      localMF[fine_flux_MF+dir]->boxArray())
   amrex::Error("fine_flux_MF invalid box array");
 } // dir=0..sdim-1

 if ((level>=1)&&(level<=finest_level)) {
  int c_level=level-1;
  NavierStokes& coarse_lev=getLevel(c_level);
  const DistributionMapping& fdmap=dmap;

  BoxArray crse_S_fine_BA(grids.size());
  for (int i = 0; i < grids.size(); ++i) {
   crse_S_fine_BA.set(i,amrex::coarsen(grids[i],2));
  }
  resize_metrics(1);
  debug_ngrow(VOLUME_MF,0,700);
  coarse_lev.resize_metrics(1);
  coarse_lev.debug_ngrow(VOLUME_MF,0,700);

  resize_mask_nbr(1);
   // (1) =1 interior  =1 fine-fine ghost in domain  =0 otherwise
   // (2) =1 interior  =0 otherwise
  debug_ngrow(MASK_NBR_MF,1,253); 

  MultiFab* coarse_mask_sem_fine;
  DistributionMapping crse_dmap=fdmap;
  coarse_mask_sem_fine=new MultiFab(crse_S_fine_BA,crse_dmap,1,1,
	MFInfo().SetTag("coarse_mask_sem_fine"),FArrayBoxFactory());
   // FabArray.H     
   // scomp,dcomp,ncomp,s_nghost,d_nghost
  coarse_mask_sem_fine->copy(*coarse_lev.localMF[MASKSEM_MF],0,0,
    1,1,1,geom.periodicity());

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

   MultiFab& S_fine_MAC=*localMF[fine_flux_MF+dir];
   const BoxArray& fgridsMAC=S_fine_MAC.boxArray();

   if (localMF[AREA_MF+dir]->boxArray()!=fgridsMAC) 
    amrex::Error("fine fine_flux_MF+dir invalid boxArray");

   if (grids.size()!=fgridsMAC.size())
    amrex::Error("grids.size()!=fgridsMAC.size()");

   BoxArray crse_S_fine_BA_MAC(fgridsMAC.size());
   for (int i = 0; i < fgridsMAC.size(); ++i) {
    crse_S_fine_BA_MAC.set(i,amrex::coarsen(fgridsMAC[i],2));
   }

   MultiFab crse_S_fine_MAC(crse_S_fine_BA_MAC,crse_dmap,ncomp_flux,0,
     MFInfo().SetTag("crse_S_fine"),FArrayBoxFactory());

     // overwrite covered fluxes with valid neighbor fluxes.
   coarse_lev.sync_flux_var(dir,coarse_flux_MF,ncomp_flux);

   crse_S_fine_MAC.copy(*coarse_lev.localMF[coarse_flux_MF+dir],0,0,
    ncomp_flux,0,0,geom.periodicity());

   ParallelDescriptor::Barrier();

   const Real* dx = geom.CellSize();
   const Real* dxc = coarse_lev.geom.CellSize();
   const Real* prob_lo   = geom.ProbLo();
 
   if (thread_class::nthreads<1)
    amrex::Error("thread_class::nthreads invalid");
   thread_class::init_d_numPts(localMF[MASKCOEF_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
   for (MFIter mfi(*localMF[MASKSEM_MF],false); mfi.isValid(); ++mfi) {
    BL_ASSERT(grids[mfi.index()] == mfi.validbox());
    const int gridno = mfi.index();

    const Box& tilegrid = mfi.tilebox();
    const Box& fabgrid = grids[gridno];
    const int* tilelo=tilegrid.loVect();
    const int* tilehi=tilegrid.hiVect();
    const int* fablo=fabgrid.loVect();
    const int* fabhi=fabgrid.hiVect();

    const Real* xlo = grid_loc[gridno].lo();

    const Box& ovgrid = crse_S_fine_BA[gridno];
    const int* ovlo=ovgrid.loVect();
    const int* ovhi=ovgrid.hiVect();

    FArrayBox& fine_fab=S_fine_MAC[gridno];
    const Box& fgridMAC=fine_fab.box();
    const int* flo=fgridMAC.loVect();
    const int* fhi=fgridMAC.hiVect();
    const Real* f_dat=fine_fab.dataPtr();

    FArrayBox& coarse_fab=crse_S_fine_MAC[gridno];
    const Box& cgridMAC=coarse_fab.box();
    const int* clo=cgridMAC.loVect();
    const int* chi=cgridMAC.hiVect();
    const Real* c_dat=coarse_fab.dataPtr();

    int bfact_f=parent->Space_blockingFactor(level);
    int bfact_c=parent->Space_blockingFactor(c_level);

    FArrayBox& maskfab=(*localMF[MASKSEM_MF])[mfi];
    FArrayBox& cmaskfab=(*coarse_mask_sem_fine)[mfi];
    FArrayBox& mnbrfab=(*localMF[MASK_NBR_MF])[mfi];
    Vector<int> velbc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    FORT_INTERP_FLUX( 
     &enable_spectral,
     dxc,dx,
     &finest_level,
     tilelo,tilehi,
     fablo,fabhi,
     &dir,
     prob_lo,
     &c_level,
     &level,
     &bfact_c,
     &bfact_f,     
     xlo,
     &ncomp_flux,
     f_dat,ARLIM(flo),ARLIM(fhi),
     c_dat,ARLIM(clo),ARLIM(chi),
     mnbrfab.dataPtr(),
     ARLIM(mnbrfab.loVect()),ARLIM(mnbrfab.hiVect()),
     maskfab.dataPtr(),
     ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
     cmaskfab.dataPtr(),
     ARLIM(cmaskfab.loVect()),ARLIM(cmaskfab.hiVect()),
     velbc.dataPtr(),
     ovlo,ovhi);  // ovlo,ovhi = crse_S_fine_BA[gridno]

   }// mfi
} //omp
   ns_reconcile_d_num(131);
  } // dir=0..sdim-1

  delete coarse_mask_sem_fine;

 } else if (level==0) {
  // do nothing
 } else
  amrex::Error("level invalid11");

} // subroutine interp_flux_localMF

// interpolate from level+1 to level.
// spectral_override==0 => always do low order average down.
void NavierStokes::avgDownEdge_localMF(int idxMF,int scomp,int ncomp,
  int start_dir,int ndir,int spectral_override,int caller_id) {

 if (1==0) {
  std::cout << "avgDownEdge_localMF caller_id= " << caller_id << '\n';
 }

 int finest_level=parent->finestLevel();
 if (level<finest_level) {
  NavierStokes& ns_fine=getLevel(level+1);
  for (int dir=start_dir;dir<start_dir+ndir;dir++) {
   debug_ngrow(idxMF+dir,0,1000+dir);
   ns_fine.debug_ngrow(idxMF+dir,0,1500+dir);
   MultiFab& S_crse=*localMF[idxMF+dir];
   MultiFab& S_fine=*ns_fine.localMF[idxMF+dir];
   int caller_id_alt=4;
   avgDownEdge(dir,S_crse,S_fine,scomp,ncomp,spectral_override,caller_id_alt);
  } // dir
 }

} // avgDownEdge_localMF

// input: XD at MAC locations, levelset function(s), normal(s)
// output: get_new_data(Umac_Type+dir,slab_step+1)  dir=0..sdim-1
// for viscosity:
// the ghost Tau is (I-nn)T^interior (I-nn) + nn T^exterior nn
//   (what about derivatives normal to the face for viscosity?)
//   (jump condition aware extrapolation?)
void NavierStokes::MAC_GRID_ELASTIC_FORCE(int im_elastic) {

 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

 if ((im_elastic>=0)&&(im_elastic<nmat)) {
  if ((particles_flag==1)||
      (particles_flag==0)) { 
   if (ns_is_rigid(im_elastic)==0) {
    if ((elastic_time[im_elastic]>0.0)&&
        (elastic_viscosity[im_elastic]>0.0)) {
     if (store_elastic_data[im_elastic]==1) {
      // do nothing
     } else
      amrex::Error("expecting store_elastic_data[im_elastic]==1");
    } else
     amrex::Error("expecting elastic_time>0 and elastic_viscosity>0");
   } else
    amrex::Error("expecting ns_is_rigid(im_elastic)==0)"); 
  } else
   amrex::Error("particls_flag invalid");
 } else
  amrex::Error("im_elastic invalid");

 int partid=0;
 while ((im_elastic_map[partid]!=im_elastic)&&
        (partid<im_elastic_map.size())) {
  partid++;
 }
 if (partid<im_elastic_map.size()) {
  // do nothing
 } else
  amrex::Error("partid invalid");

 resize_levelsetLO(2,LEVELPC_MF);
 debug_ngrow(LEVELPC_MF,2,103);
 if (localMF[LEVELPC_MF]->nComp()!=nmat*(AMREX_SPACEDIM+1))
  amrex::Error("(localMF[LEVELPC_MF]->nComp()!=nmat*(AMREX_SPACEDIM+1))");

  // https://ccse.lbl.gov/
  // https://github.com/AMReX-Codes/amrex/
 bool use_tiling=ns_tiling;

  // density, temperature
 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int finest_level=parent->finestLevel();

  // areas and volumes associated with rectangular grid cells and faces.
 resize_metrics(1);
  // if MASKCOEF_MF == 1 then the cell is not covered by a finer level.
 resize_maskfiner(1,MASKCOEF_MF);
  // mask_nbr=1 at fine-fine boundary conditions
 resize_mask_nbr(1);

 debug_ngrow(VOLUME_MF,1,100);
 debug_ngrow(MASKCOEF_MF,1,253); // maskcoef=1 if not covered by finer level.
 debug_ngrow(MASK_NBR_MF,1,253); // mask_nbr=1 at fine-fine bc.
 debug_ngrow(CELL_VISC_MATERIAL_MF,1,3);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  debug_ngrow(FACE_VAR_MF+dir,0,101); // faceden_index has the MAC density

  MultiFab& UMAC_new=get_new_data(Umac_Type+dir,slab_step+1);
  if (UMAC_new.nGrow()==0) {
   // do nothing
  } else
   amrex::Error("UMAC_new invalid ngrow");

  if (UMAC_new.nComp()!=1)
   amrex::Error("UMAC_new.nComp()!=1");

  if (localMF[AREA_MF+dir]->boxArray()!=UMAC_new.boxArray())
   amrex::Error("localMF[AREA_MF+dir]->boxArray()!=UMAC_new.boxArray()");
 } // dir=0..sdim-1

 VOF_Recon_resize(2,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,2,118);

 const Box& domain = geom.Domain();
 const int* domlo = domain.loVect();
 const int* domhi = domain.hiVect();

 const Real* dx = geom.CellSize();

  // outer loop: each force component u_t = F_elastic/density  u,v,w
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

  MultiFab& UMAC_new=get_new_data(Umac_Type+dir,slab_step+1);
  if (UMAC_new.nGrow()==0) {
   // do nothing
  } else
   amrex::Error("UMAC_new invalid ngrow");

  if (UMAC_new.nComp()!=1)
   amrex::Error("UMAC_new.nComp()!=1");

  if (localMF[AREA_MF+dir]->boxArray()!=UMAC_new.boxArray())
   amrex::Error("localMF[AREA_MF+dir]->boxArray()!=UMAC_new.boxArray()");

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(localMF[LEVELPC_MF]->boxArray().d_numPts());
 
#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(*localMF[LEVELPC_MF],use_tiling); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();
    // Gaussian quadrature point positioning
   int bfact=parent->Space_blockingFactor(level);

   const Real* xlo = grid_loc[gridno].lo();

    // fab = Fortran Array Block

   int grid_type_CC=-1;
   FArrayBox& MAC_CCfab=(*localMF[MAC_ELASTIC_FLUX_CC_MF])[mfi];

   int grid_type_XY=3;
   FArrayBox& MAC_XYfab=(*localMF[MAC_ELASTIC_FLUX_XY_MF])[mfi];

#if (AMREX_SPACEDIM==2)
   int grid_type_XZ=3;
   FArrayBox& MAC_XZfab=(*localMF[MAC_ELASTIC_FLUX_XY_MF])[mfi];
   int grid_type_YZ=3;
   FArrayBox& MAC_YZfab=(*localMF[MAC_ELASTIC_FLUX_XY_MF])[mfi];
#elif (AMREX_SPACEDIM==3)
   int grid_type_XZ=4;
   FArrayBox& MAC_XZfab=(*localMF[MAC_ELASTIC_FLUX_XZ_MF])[mfi];
   int grid_type_YZ=5;
   FArrayBox& MAC_YZfab=(*localMF[MAC_ELASTIC_FLUX_YZ_MF])[mfi];
#else
   amrex::Error("dimension bust");
#endif

   FArrayBox& xface=(*localMF[FACE_VAR_MF+dir])[mfi];

    // output
   FArrayBox& UMACNEWfab=UMAC_new[mfi];

   // mask=1.0 at interior fine bc ghost cells
   FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi];
   // maskcoef=1 if not covered by finer level or outside domain
   FArrayBox& maskcoef=(*localMF[MASKCOEF_MF])[mfi]; 

   FArrayBox& levelpcfab=(*localMF[LEVELPC_MF])[mfi];

   FArrayBox& viscfab=(*localMF[CELL_VISC_MATERIAL_MF])[mfi];
   int ncomp_visc=viscfab.nComp();

   Vector<int> velbc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);

   int rzflag=0;
   if (geom.IsRZ())
    rzflag=1;
   else if (geom.IsCartesian())
    rzflag=0;
   else if (geom.IsCYLINDRICAL())
    rzflag=3;
   else
    amrex::Error("CoordSys bust 21");

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   if (viscoelastic_model[im_elastic]==2) {

    // declared in: GODUNOV_3D.F90
    fort_mac_elastic_force(
     &im_elastic,
     &partid,
     &dir, // dir=0,1,..sdim-1  
     &ncomp_visc, 
     &visc_coef,
     &facevisc_index,
     &faceden_index,
     &massface_index,
     &vofface_index,
     &ncphys,
     velbc.dataPtr(),
     &dt_slab,
     &cur_time_slab,
     xlo,dx,
     &grid_type_CC,
     MAC_CCfab.dataPtr(),
     ARLIM(MAC_CCfab.loVect()),ARLIM(MAC_CCfab.hiVect()),
     &grid_type_XY,
     MAC_XYfab.dataPtr(),
     ARLIM(MAC_XYfab.loVect()),ARLIM(MAC_XYfab.hiVect()),
     &grid_type_XZ,
     MAC_XZfab.dataPtr(),
     ARLIM(MAC_XZfab.loVect()),ARLIM(MAC_XZfab.hiVect()),
     &grid_type_YZ,
     MAC_YZfab.dataPtr(),
     ARLIM(MAC_YZfab.loVect()),ARLIM(MAC_YZfab.hiVect()),
     viscfab.dataPtr(),
     ARLIM(viscfab.loVect()),ARLIM(viscfab.hiVect()),
     maskfab.dataPtr(), // mask=1.0 at interior fine bc ghost cells
     ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
     //maskcoef=1 if not covered by finer level or outside
     maskcoef.dataPtr(), 
     ARLIM(maskcoef.loVect()),ARLIM(maskcoef.hiVect()),
     levelpcfab.dataPtr(),
     ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),
     xface.dataPtr(),
     ARLIM(xface.loVect()),ARLIM(xface.hiVect()), 
     UMACNEWfab.dataPtr(),
     ARLIM(UMACNEWfab.loVect()),ARLIM(UMACNEWfab.hiVect()), 
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,
     &level,&finest_level,
     &rzflag,
     domlo,domhi,
     &nmat,
     &nten);

   } else if ((viscoelastic_model[im_elastic]==1)||
              (viscoelastic_model[im_elastic]==0)||
	      (viscoelastic_model[im_elastic]==3)) { //incremental
    // do nothing
   } else
    amrex::Error("viscoelastic_model[im_elastic] invalid");

  } // mfi
} // omp
  ns_reconcile_d_num(132);
 } // dir = 0..sdim-1

 make_MAC_velocity_consistent();

} // end subroutine MAC_GRID_ELASTIC_FORCE

// PEDGE_MF allocated in allocate_pressure_work_vars
void NavierStokes::apply_cell_pressure_gradient(
 int project_option,
 int energyflag,
 int idx_pres,    //nsolve=1
 int idx_umac,    //nsolve=1
 int idx_gpcell,  //sdim 
 int idx_divup) { //nsolve=1

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid apply_cell_pressure_gradient");

 if ((project_option==0)||
     (project_option==1)||
     (project_option==11)) {  //FSI_material_exists (last project)
  // do nothing
 } else
  amrex::Error("project_option invalid20");

 int num_colors=0;
 Vector<Real> blob_array;
 blob_array.resize(1);
 int blob_array_size=blob_array.size();
 
 bool use_tiling=ns_tiling;

 int nmat=num_materials;
 int nsolve=1;

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int finest_level=parent->finestLevel();

 resize_metrics(1);
 resize_maskfiner(1,MASKCOEF_MF);
 resize_mask_nbr(1);

 debug_ngrow(VOLUME_MF,0,100);
 debug_ngrow(MASKCOEF_MF,1,253); // maskcoef=1 if not covered by finer level.
 debug_ngrow(MASK_NBR_MF,1,253); // mask_nbr=1 at fine-fine bc.

  // grad p = div (pI) approx 
  //  (1/|omega|) \integral_\partial\Omega pI dot n DA=
  //  (1/|omega|) \integral_\partial\Omega p n DA
  // each face needs a left and right pressure in order to account for
  // the contribution from internal faces.
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  debug_ngrow(FACE_VAR_MF+dir,0,101);
  debug_ngrow(PEDGE_MF+dir,0,101);
  if (localMF[PEDGE_MF+dir]->nComp()!=2+nsolve)
   amrex::Error("pedge_mf invalid ncomp");
  if (localMF[AREA_MF+dir]->boxArray()!=
      localMF[PEDGE_MF+dir]->boxArray())
   amrex::Error("PEDGE boxarray does not match");
  setVal_localMF(PEDGE_MF+dir,1.0e+40,2,1,0);
 } // dir=0..sdim-1

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

 resize_levelsetLO(2,LEVELPC_MF);

 debug_ngrow(LEVELPC_MF,2,103);
 if (localMF[LEVELPC_MF]->nComp()!=nmat*(AMREX_SPACEDIM+1))
  amrex::Error("(localMF[LEVELPC_MF]->nComp()!=nmat*(AMREX_SPACEDIM+1))");

 VOF_Recon_resize(1,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,1,118);

 int fluxvel_index=0;
 int fluxden_index=AMREX_SPACEDIM;

 const Box& domain = geom.Domain();
 const int* domlo = domain.loVect();
 const int* domhi = domain.hiVect();

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

 const Real* dx = geom.CellSize();

 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

 MultiFab* presmf=localMF[idx_pres];
 if (presmf->nComp()!=nsolve)
  amrex::Error("presmf->nComp() invalid");

  // old cell velocity before application of pressure gradient.
 MultiFab* ustar;

 if (enable_spectral!=0) {
  std::cout << "non-conservative: \n";
  std::cout << "divup (1) rho div u , (2) p div u \n";
  std::cout << "conservative: \n";
  std::cout << "divup (1) 0 , (2) div (up) \n";
  amrex::Error("upgrade space time spectral element");
 }

 MultiFab* divup;
 if ((energyflag==0)|| //do not update the energy
     (energyflag==1)) {//update the energy
  ustar=getState(1,0,AMREX_SPACEDIM,cur_time_slab);
  divup=new MultiFab(grids,dmap,nsolve,0,
   MFInfo().SetTag("divup"),FArrayBoxFactory());

  //Spectral deferred correction:
  //get grad p,div(up) instead of \pm dt grad p/rho, -dt div(up)/rho
 } else if (energyflag==2) { 
  debug_ngrow(idx_gpcell,0,101);
  debug_ngrow(idx_divup,0,101);
  if (localMF[idx_gpcell]->nComp()!=AMREX_SPACEDIM)
   amrex::Error("idx_gpcell has invalid ncomp");
  if (localMF[idx_divup]->nComp()!=nsolve)
   amrex::Error("idx_divup has invalid ncomp");

  ustar=localMF[idx_gpcell];
  divup=localMF[idx_divup];
 } else
  amrex::Error("energyflag invalid");

 Vector<int> scomp;
 Vector<int> ncomp;
 int ncomp_check;
 int state_index;
  //num_materials_combine=1
 get_mm_scomp_solver(
  1,
  project_option,
  state_index,
  scomp,ncomp,ncomp_check);

 if (ncomp_check!=nsolve)
  amrex::Error("invalid ncomp");

 int scomp_den=(AMREX_SPACEDIM+1);
 int nden=nmat*num_state_material;

 int nvof=nmat*ngeom_raw;
 int nstate=S_new.nComp();
 if (nstate!=scomp_den+nden+nvof+1)
  amrex::Error("invalid ncomp in cell pressure gradient routine");

 //interpolate pressure from cell to MAC grid.
 int operation_flag_interp_pres=1; 
 // flux register is initialized to zero.
 allocate_flux_register(operation_flag_interp_pres);
 if (localMF[SEM_FLUXREG_MF]->nComp()!=AMREX_SPACEDIM)
  amrex::Error("localMF[SEM_FLUXREG_MF]->nComp() invalid1");

 if (level<finest_level) {
  avgDown_and_Copy_localMF(
    idx_pres,
    idx_pres,
    AMRSYNC_PEDGE_MF,
    operation_flag_interp_pres);
 } else if (level==finest_level) {
  // do nothing
 } else
  amrex::Error("level invalid12");

 if ((level>=1)&&(level<=finest_level)) {
  interp_and_Copy_localMF(
    idx_pres,
    idx_pres,
    AMRSYNC_PEDGE_MF,
    operation_flag_interp_pres);
 } else if (level==0) {
   // do nothing
 } else
  amrex::Error("level invalid13");

 for (int spectral_loop=0;spectral_loop<end_spectral_loop();spectral_loop++) {
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
 for (int tileloop=0;tileloop<=1;tileloop++) {
 
  // interpolate pressure to MAC grid.
  // set flag for whether to update cell velocity conservatively
  // or non-conservatively.
  // Modify MAC velocity with solid velocity or ice velocity.
  
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
   int bfact_c=bfact;
   int bfact_f=bfact;
   if (level>0)
    bfact_c=parent->Space_blockingFactor(level-1);
   if (level<finest_level)
    bfact_f=parent->Space_blockingFactor(level+1);

   const Real* xlo = grid_loc[gridno].lo();

   FArrayBox& xvel=(*localMF[idx_umac+dir])[mfi];
   FArrayBox& xface=(*localMF[FACE_VAR_MF+dir])[mfi];

   FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi];  

   FArrayBox& xp=(*localMF[PEDGE_MF+dir])[mfi];
   FArrayBox& xgp=(*localMF[AMRSYNC_PEDGE_MF+dir])[mfi];

   // mask=1.0 at interior fine bc ghost cells
   FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi];
   // maskcoef=1 if not covered by finer level or outside domain
   FArrayBox& maskcoef=(*localMF[MASKCOEF_MF])[mfi]; 

   FArrayBox& maskSEMfab=(*localMF[MASKSEM_MF])[mfi];
   FArrayBox& presfab=(*presmf)[mfi];
 
   FArrayBox& solfab=(*localMF[FSI_GHOST_MAC_MF+dir])[mfi];
   FArrayBox& levelpcfab=(*localMF[LEVELPC_MF])[mfi];

   FArrayBox& semfluxfab=(*localMF[SEM_FLUXREG_MF])[mfi];
   int ncfluxreg=semfluxfab.nComp();

   Vector<int> presbc;
   getBCArray_list(presbc,state_index,gridno,scomp,ncomp);
   if (presbc.size()!=nsolve*AMREX_SPACEDIM*2)
    amrex::Error("presbc.size() invalid");
   Vector<int> velbc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);

   Real beta=0.0;

   int rzflag=0;
   if (geom.IsRZ())
    rzflag=1;
   else if (geom.IsCartesian())
    rzflag=0;
   else if (geom.IsCYLINDRICAL())
    rzflag=3;
   else
    amrex::Error("CoordSys bust 21");

   int local_energyflag=0;
   int local_enable_spectral=enable_spectral;
   int simple_AMR_BC_flag=0;
   int ncomp_xp=2+nsolve;
   int ncomp_xgp=1;
   int ncomp_mgoni=presfab.nComp();

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // in apply_cell_pressure_gradient: p^CELL -> p^MAC 
   // AMR transfer data is in the 3rd component of xp.
   fort_cell_to_mac(
    &ncomp_mgoni, 
    &ncomp_xp, 
    &ncomp_xgp, 
    &simple_AMR_BC_flag,
    &nsolve,
    &tileloop,
    &dir,
    &operation_flag_interp_pres, //1
    &local_energyflag,
    &beta,
    &visc_coef,
    &interp_vel_increment_from_cell,
    filter_velocity.dataPtr(),
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
    &ignore_div_up,
    &pforce_index,
    &faceden_index,
    &icemask_index,
    &massface_index,
    &vofface_index,
    &ncphys,
    override_density.dataPtr(),
    constant_density_all_time.dataPtr(),
    presbc.dataPtr(),
    velbc.dataPtr(),
    &slab_step,
    &dt_slab,
    &cur_time_slab,
    xlo,dx,
    &spectral_loop,
    &ncfluxreg,
     // semfluxfab(i,j,k)  (as an array4)
     //   semfluxfab.loVect(1)<=i<=semfluxfab.hiVect(1)
     //   semfluxfab.loVect(2)<=j<=semfluxfab.hiVect(2)
     //   semfluxfab.loVect(3)<=k<=semfluxfab.hiVect(3)
    semfluxfab.dataPtr(),
     //ARLIM(x) is a macro
     //the compiler expands this macro as:
     //x(1),(2),x(3)
     //here (in 3 dimensions),
     //semfluxfab.loVect(1),semfluxfab.loVect(2),semfluxfab.loVect(3)
    ARLIM(semfluxfab.loVect()),ARLIM(semfluxfab.hiVect()),
    maskfab.dataPtr(), // mask=1.0 at interior fine bc ghost cells
    ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
    maskcoef.dataPtr(), // maskcoef=1 if not covered by finer level or outside
    ARLIM(maskcoef.loVect()),ARLIM(maskcoef.hiVect()),
    maskSEMfab.dataPtr(),
    ARLIM(maskSEMfab.loVect()),ARLIM(maskSEMfab.hiVect()),
    levelpcfab.dataPtr(),
    ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),
    solfab.dataPtr(),
    ARLIM(solfab.loVect()),ARLIM(solfab.hiVect()),
    xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()), //xcut
    xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()), 
    reconfab.dataPtr(),ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
    xgp.dataPtr(),ARLIM(xgp.loVect()),ARLIM(xgp.hiVect()),//holds AMRSYNC_PEDGE
    xp.dataPtr(),ARLIM(xp.loVect()),ARLIM(xp.hiVect()), 
    xvel.dataPtr(),ARLIM(xvel.loVect()),ARLIM(xvel.hiVect()), 
    presfab.dataPtr(),ARLIM(presfab.loVect()),ARLIM(presfab.hiVect()), //vel
    presfab.dataPtr(),ARLIM(presfab.loVect()),ARLIM(presfab.hiVect()), 
    presfab.dataPtr(),ARLIM(presfab.loVect()),ARLIM(presfab.hiVect()), //den
    presfab.dataPtr(),ARLIM(presfab.loVect()),ARLIM(presfab.hiVect()), //mgoni
    presfab.dataPtr(),ARLIM(presfab.loVect()),ARLIM(presfab.hiVect()), //color
    presfab.dataPtr(),ARLIM(presfab.loVect()),ARLIM(presfab.hiVect()), //type
    tilelo,tilehi,
    fablo,fabhi,
    &bfact,&bfact_c,&bfact_f,
    &level,&finest_level,
    &rzflag,
    domlo,domhi,
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
    &project_option,
    &SEM_upwind,
    &SEM_advection_algorithm);
  } // mfi
} // omp
  ns_reconcile_d_num(132);
 } // tileloop
 } // dir

 synchronize_flux_register(operation_flag_interp_pres,spectral_loop);
 } // spectral_loop


  // 0=use_face_pres
  // 1=grid flag (coarse/fine boundary?)
  // 2=face pressure
 if (nsolve!=1)
  amrex::Error("nsolve!=1");
 int pface_comp=2;
 int ncomp_edge=nsolve;
 int caller_id=5;
 avgDownEdge_localMF(PEDGE_MF,pface_comp,ncomp_edge,0,AMREX_SPACEDIM,1,caller_id);

  // isweep=1 calculate cell velocity from mass weighted average of face
  // velocity.
  // isweep=2 calculate cell pressure gradient, update cell velocity,
  //  update density (if non conservative), update energy.
 for (int isweep=1;isweep<=2;isweep++) {

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

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();
   int bfact=parent->Space_blockingFactor(level);

   const Real* xlo = grid_loc[gridno].lo();

   FArrayBox& xvel=(*localMF[idx_umac])[mfi];
   FArrayBox& yvel=(*localMF[idx_umac+1])[mfi];
   FArrayBox& zvel=(*localMF[idx_umac+AMREX_SPACEDIM-1])[mfi];

   FArrayBox& xface=(*localMF[FACE_VAR_MF])[mfi];
   FArrayBox& yface=(*localMF[FACE_VAR_MF+1])[mfi];
   FArrayBox& zface=(*localMF[FACE_VAR_MF+AMREX_SPACEDIM-1])[mfi];
   FArrayBox& ax = (*localMF[AREA_MF])[mfi];
   FArrayBox& ay = (*localMF[AREA_MF+1])[mfi];
   FArrayBox& az = (*localMF[AREA_MF+AMREX_SPACEDIM-1])[mfi];
   FArrayBox& vol=(*localMF[VOLUME_MF])[mfi];
   FArrayBox& solxfab=(*localMF[FSI_GHOST_MAC_MF])[mfi];
   FArrayBox& solyfab=(*localMF[FSI_GHOST_MAC_MF+1])[mfi];
   FArrayBox& solzfab=(*localMF[FSI_GHOST_MAC_MF+AMREX_SPACEDIM-1])[mfi];

   FArrayBox& levelpcfab=(*localMF[LEVELPC_MF])[mfi];

   FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi];
   FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi];//1=fine/fine bc
   FArrayBox& maskcoef=(*localMF[MASKCOEF_MF])[mfi];// 1=not covered.
   FArrayBox& maskSEMfab=(*localMF[MASKSEM_MF])[mfi];
   FArrayBox& presfab=(*presmf)[mfi];
   FArrayBox& xp=(*localMF[PEDGE_MF])[mfi];
   FArrayBox& yp=(*localMF[PEDGE_MF+1])[mfi];
   FArrayBox& zp=(*localMF[PEDGE_MF+AMREX_SPACEDIM-1])[mfi];
   FArrayBox& Snewfab=S_new[mfi]; // veldest
   FArrayBox& ustarfab=(*ustar)[mfi];
   FArrayBox& divupfab=(*divup)[mfi];

   Vector<int> presbc;
   getBCArray_list(presbc,state_index,gridno,scomp,ncomp);
   if (presbc.size()!=nsolve*AMREX_SPACEDIM*2)
    amrex::Error("presbc.size() invalid");
   Vector<int> velbc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);

   int operation_flag_interp_macvel;
   if (isweep==1)
    operation_flag_interp_macvel=2; // mac velocity -> cell velocity
   else if (isweep==2)
    operation_flag_interp_macvel=3; // (grad p)_CELL, div(up)
   else
    amrex::Error("operation_flag_interp_macvel invalid1");

   int homflag=0; // default
 
   int ok_to_call=0;
   if ((energyflag==0)||
       (energyflag==1)) {
    ok_to_call=1;

     // just find (grad p)_CELL, div(up) for space-time algorithm
   } else if (energyflag==2) { 
    if (isweep==1) { // disable mac velocity->cell velocity step here.
     ok_to_call=0;
    } else if (isweep==2) { // (grad p)_{CELL}, div(up) only
     ok_to_call=1;
    } else
     amrex::Error("isweep invalid");
   } else
    amrex::Error("energyflag invalid");
 
   // in apply_cell_pressure_gradient 

   if (ok_to_call==1) {
    int local_enable_spectral=enable_spectral;
    int use_VOF_weight=1;

    int ncomp_denold=presfab.nComp();
    int ncomp_veldest=Snewfab.nComp();
    int ncomp_dendest=Snewfab.nComp()-scomp_den;

    fort_mac_to_cell(
     &ns_time_order,
     &divu_outer_sweeps,
     &num_divu_outer_sweeps,
     // 2 (mac_vel->cell_vel or 3 (cell gradp, cell energy)
     &operation_flag_interp_macvel, 
     &energyflag,
     temperature_primitive_variable.dataPtr(),
     constant_density_all_time.dataPtr(),
     &nmat,
     &nparts,
     &nparts_def,
     im_solid_map_ptr,
     added_weight.dataPtr(),
     &nten,
     &level, 
     &finest_level,
     &project_option,
     &local_enable_spectral,
     &fluxvel_index,
     &fluxden_index,
     &facevel_index,
     &facecut_index,
     &icefacecut_index,
     &curv_index,
     &conservative_tension_force,
     &conservative_div_uu,
     filter_velocity.dataPtr(),
     &ignore_div_up,
     &pforce_index,
     &faceden_index,
     &icemask_index,
     &massface_index,
     &vofface_index,
     &ncphys,
     velbc.dataPtr(),
     presbc.dataPtr(), 
     &cur_time_slab, 
     &slab_step,
     &dt_slab,
     xlo,dx,
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,
     xp.dataPtr(),ARLIM(xp.loVect()),ARLIM(xp.hiVect()),
     yp.dataPtr(),ARLIM(yp.loVect()),ARLIM(yp.hiVect()),
     zp.dataPtr(),ARLIM(zp.loVect()),ARLIM(zp.hiVect()),
     xvel.dataPtr(),ARLIM(xvel.loVect()),ARLIM(xvel.hiVect()),
     yvel.dataPtr(),ARLIM(yvel.loVect()),ARLIM(yvel.hiVect()),
     zvel.dataPtr(),ARLIM(zvel.loVect()),ARLIM(zvel.hiVect()),
     xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()),
     yface.dataPtr(),ARLIM(yface.loVect()),ARLIM(yface.hiVect()),
     zface.dataPtr(),ARLIM(zface.loVect()),ARLIM(zface.hiVect()),
     ax.dataPtr(),ARLIM(ax.loVect()),ARLIM(ax.hiVect()),
     ay.dataPtr(),ARLIM(ay.loVect()),ARLIM(ay.hiVect()),
     az.dataPtr(),ARLIM(az.loVect()),ARLIM(az.hiVect()),
     vol.dataPtr(),ARLIM(vol.loVect()),ARLIM(vol.hiVect()),
     divupfab.dataPtr(), // rhs
     ARLIM(divupfab.loVect()),ARLIM(divupfab.hiVect()),
     Snewfab.dataPtr(),
     ARLIM(Snewfab.loVect()),ARLIM(Snewfab.hiVect()), // veldest
      // scomp_den=(AMREX_SPACEDIM+1);
     Snewfab.dataPtr(scomp_den),
     ARLIM(Snewfab.loVect()),ARLIM(Snewfab.hiVect()), // dendest
     maskfab.dataPtr(), // 1=fine/fine  0=coarse/fine
     ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
     maskcoef.dataPtr(), // 1=not covered  0=covered
     ARLIM(maskcoef.loVect()),ARLIM(maskcoef.hiVect()),
     maskSEMfab.dataPtr(), 
     ARLIM(maskSEMfab.loVect()),ARLIM(maskSEMfab.hiVect()),
     levelpcfab.dataPtr(), //levelPC
     ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),
     solxfab.dataPtr(),ARLIM(solxfab.loVect()),ARLIM(solxfab.hiVect()),
     solyfab.dataPtr(),ARLIM(solyfab.loVect()),ARLIM(solyfab.hiVect()),
     solzfab.dataPtr(),ARLIM(solzfab.loVect()),ARLIM(solzfab.hiVect()),
     levelpcfab.dataPtr(),
     ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),//cterm
     presfab.dataPtr(), 
     ARLIM(presfab.loVect()),ARLIM(presfab.hiVect()),//pold
     presfab.dataPtr(),
     ARLIM(presfab.loVect()),ARLIM(presfab.hiVect()),//denold
     ustarfab.dataPtr(),ARLIM(ustarfab.loVect()),ARLIM(ustarfab.hiVect()),
     reconfab.dataPtr(),ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
     levelpcfab.dataPtr(),
     ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),//mdot
     levelpcfab.dataPtr(),
     ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),//maskdivres
     levelpcfab.dataPtr(),
     ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),//maskres
     &SDC_outer_sweeps,
     &homflag,
     &use_VOF_weight,
     &nsolve,
     &ncomp_denold,
     &ncomp_veldest,
     &ncomp_dendest,
     &SEM_advection_algorithm);

   } else if (ok_to_call==0) {
    // do nothing
   } else
    amrex::Error("ok_to_call invalid");

  }   // mfi
} // omp
  ns_reconcile_d_num(133);

 }  // isweep=1,2

 if ((energyflag==0)||
     (energyflag==1)) {
  save_to_macvel_state(idx_umac);
  delete divup; // div(up) is discarded.
  delete ustar;
 } else if (energyflag==2) { // (grad p)_CELL, div(up) for space-time
  // do nothing
 } else
  amrex::Error("energyflag invalid");

} // subroutine apply_cell_pressure_gradient

void NavierStokes::save_to_macvel_state(int idx_umac) {

 int nsolve=1;
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  debug_ngrow(idx_umac+dir,0,111);
  if (localMF[idx_umac+dir]->nComp()==1) {
   MultiFab& Umac_new=get_new_data(Umac_Type+dir,slab_step+1);
   MultiFab::Copy(Umac_new,*localMF[idx_umac+dir],0,0,nsolve,0);
  } else
   amrex::Error("localMF[idx_umac+dir]->nComp() invalid");
 }  // dir=0..sdim-1

} // save_to_macvel_state

void NavierStokes::grid_type_to_box_type_cpp(int grid_type,
		int* box_type) {

 for (int local_dir=0;local_dir<AMREX_SPACEDIM;local_dir++) {
  box_type[local_dir]=0;
 }
 if (grid_type==-1) {
  // do nothing
 } else if ((grid_type>=0)&&(grid_type<AMREX_SPACEDIM)) {
  box_type[grid_type]=1;
 } else if (grid_type==3) {
  box_type[0]=1;
  box_type[1]=1;
 } else if ((grid_type==4)&&(AMREX_SPACEDIM==3)) {
  box_type[0]=1;
  box_type[AMREX_SPACEDIM-1]=1;
 } else if ((grid_type==5)&&(AMREX_SPACEDIM==3)) {
  box_type[1]=1;
  box_type[AMREX_SPACEDIM-1]=1;
 } else
  amrex::Error("grid_type invalid grid_type_to_box_type_cpp");

} // end subroutine grid_type_to_box_type_cpp

void NavierStokes::get_iten_cpp(int im1,int im2,int& iten,int nmat) {

 int im=-1;
 int im_opp=-1;

 if (nmat<1)
  amrex::Error("nmat invalid");

 if ((im1<1)||(im1>nmat)||
     (im2<1)||(im2>nmat)||
     (im1==im2)) {
  std::cout << "im1,im2 mismatch im1,im2=" << im1 << ' ' << im2 << '\n';
  std::cout << "nmat=" << nmat << '\n';
  amrex::Error("get_iten_cpp problem "); 
 }

 if (im1<im2) {
  im=im1;im_opp=im2;
 } else if (im2<im1) {
  im=im2;im_opp=im1;
 } else {
  amrex::Error("im1 or im2 bust");
 }
  
 if (im==1) {
  iten=im_opp-1;
 } else if (im==2) {
  iten=nmat-1+im_opp-2;
 } else if (im==3) {
  iten=2*nmat-3+im_opp-3;
 } else {
  amrex::Error("im1 or im2 not supported yet");
 }

}  // get_iten_cpp


void NavierStokes::get_inverse_iten_cpp(int& im1,int& im2,int iten,int nmat) {

 if (iten<1) {
  std::cout << "iten= " << iten << '\n';
  amrex::Error("iten invalid in get_inverse_iten_cpp");
 }
 if (nmat<1)
  amrex::Error("nmat invalid");

 im1=0;
 im2=0;
 for (int im=1;im<=nmat;im++) {
  for (int im_opp=im+1;im_opp<=nmat;im_opp++) {
   int iten_test;
   get_iten_cpp(im,im_opp,iten_test,nmat); 
   if (iten==iten_test) {
    im1=im;
    im2=im_opp;
   }
  }
 }
 if ((im1<1)||(im1>nmat)||
     (im2<1)||(im2>nmat)||
     (im1==im2)) {
  std::cout << "im1,im2 mismatch im1,im2=" << im1 << ' ' << im2 << '\n';
  std::cout << "nmat=" << nmat << '\n';
  amrex::Error("get_inverse_iten_cpp problem "); 
 }

}  // get_inverse_iten_cpp

void NavierStokes::make_MAC_velocity_consistentALL() {

 int finest_level=parent->finestLevel();

 if (level==0) {
  // do nothing
 } else
  amrex::Error("expecting level==0 in make_MAC_velocity_consistentALL");

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.make_MAC_velocity_consistent();
 } 

} // end subroutine make_MAC_velocity_consistentALL()

void NavierStokes::make_MAC_velocity_consistent() {

 int finest_level = parent->finestLevel();

 // spectral_override==0 => always do low order average down.
 // spectral_override==1 => order derived from "enable_spectral"
 int spectral_override=1;

  // avgDown all the MAC components.
 if (level<finest_level)
  avgDownMacState(Umac_Type,spectral_override); 

 int nsolve=1;

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   // ngrow,dir,scomp,ncomp,time
  MultiFab* tempmac=getStateMAC(
     Umac_Type,0,dir,0,nsolve,cur_time_slab);
  MultiFab& Umac_new=get_new_data(Umac_Type+dir,slab_step+1);
  MultiFab::Copy(Umac_new,*tempmac,0,0,nsolve,0);
  delete tempmac;
 } // dir=0..sdim-1

}  // subroutine make_MAC_velocity_consistent()

// ucell_new=ucell+old + force_cell
//
//   POTENTIALLY UNSTABLE:
// if (filter_velocity[im]==0) and (interp_option==2) then
//  umac_new=umac_old + INTERP_CELL_TO_MAC(Force_cell)
//   VERY DISSIPATIVE:
// if (filter_velocity[im]==1) and (interp_option==2) then
//  umac_new=INTERP_CELL_TO_MAC(ucell_new)
void NavierStokes::increment_face_velocityALL(
 int interp_option,
 int project_option,
 int idx_velcell,
 Real beta,
 Vector<blobclass> blobdata) {

 if (idx_velcell==DELTA_CELL_VEL_MF)
  amrex::Error("DELTA_CELL_VEL_MF reserved");
 if (idx_velcell==CURRENT_CELL_VEL_MF)
  amrex::Error("CURRENT_CELL_VEL_MF reserved");

 int finest_level=parent->finestLevel();

 getStateALL(1,cur_time_slab,0,AMREX_SPACEDIM,DELTA_CELL_VEL_MF);
 getStateALL(1,cur_time_slab,0,AMREX_SPACEDIM,CURRENT_CELL_VEL_MF);

 make_MAC_velocity_consistentALL();

  // interp_option=4 
  //  unew^{f} = 
  //   (i) unew^{f} in incompressible non-solid regions
  //   (ii) u^{f,save} + (unew^{c}-u^{c,save})^{c->f} in spectral regions or
  //        compressible regions.
  //   (iii) usolid in solid regions
 if (interp_option==4) {
  minusALL(1,AMREX_SPACEDIM,DELTA_CELL_VEL_MF,ADVECT_REGISTER_MF);
 }

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.increment_face_velocity(
    interp_option,project_option,
    idx_velcell,
    beta,blobdata); 
   // avgDownMacState, getStateMAC to fill EXT_DIR BC.
  ns_level.make_MAC_velocity_consistent();
  ParallelDescriptor::Barrier();
 }  // ilev=finest_level ... level

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
  delete_array(AMRSYNC_VEL_MF+dir);

 delete_array(DELTA_CELL_VEL_MF);
 delete_array(CURRENT_CELL_VEL_MF);

} // end subroutine increment_face_velocityALL

// interp_option=0 unew^{f} = unew^{c->f} 
// interp_option=1 unew^{f} = unew^{f} in fluid  (=usolid in solid)
// interp_option=2 unew^{f} = unew^{f} + beta * diffuse_register^{c->f}
// interp_option=3 unew^{f} = 
//    unew^{c,f -> f} in fluid  (=usolid in solid)
//    (hybrid of interp_option==0 and 1)
// interp_option=4 unew^{f} = 
//   (i) unew^{f} in incompressible non-solid regions
//   (ii) u^{f,save} + (unew^{c}-u^{c,save})^{c->f} in spectral regions or
//        compressible regions.
//   (iii) usolid in solid regions
// called from: post_init_state, do_the_advance, multiphase_project
// (when project_option==0,1,10,11,13), APPLY_REGISTERS, INCREMENT_REGISTERS
// called from NavierStokes::increment_face_velocityALL
void NavierStokes::increment_face_velocity(
 int interp_option,
 int project_option,
 int idx_velcell,
 Real beta,
 Vector<blobclass> blobdata) {

 int finest_level = parent->finestLevel();
 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

 MultiFab* levelcolor=nullptr;
 MultiFab* leveltype=nullptr;

 int num_colors=blobdata.size();
 int blob_array_size=num_colors*num_elements_blobclass;

 Vector<Real> blob_array;
 blob_array.resize(1);

 if (num_colors==0) {

  blob_array_size=blob_array.size();

 } else if (num_colors>0) {

  blob_array.resize(blob_array_size);

  int counter=0;

  for (int i=0;i<num_colors;i++) {
   copy_from_blobdata(i,counter,blob_array,blobdata);
  }  // i=0..num_colors-1
  if (counter!=blob_array_size)
   amrex::Error("counter invalid");

 } else
  amrex::Error("num_colors invalid");

 int operation_flag=3;

 int primary_vel_data=-1;
 int secondary_vel_data=-1;

 if (interp_option==0) { // unew^{f} = unew^{c->f}

  if (idx_velcell==-1) {
   primary_vel_data=CURRENT_CELL_VEL_MF; 
   secondary_vel_data=CURRENT_CELL_VEL_MF; 
  } else
   amrex::Error("idx_velcell invalid");

  if ((project_option==0)||
      (project_option==1)) {
   // do nothing
  } else
   amrex::Error("project_option invalid21");

  if (beta!=0.0)
   amrex::Error("beta invalid");

  operation_flag=3;

 } else if (interp_option==1) { // unew^{f}=unew^{f}

  if (idx_velcell==-1) {
   primary_vel_data=CURRENT_CELL_VEL_MF; 
   secondary_vel_data=CURRENT_CELL_VEL_MF; 
  } else
   amrex::Error("idx_velcell invalid");

  if (project_option==11) {  //FSI_material_exists last project

   if (num_colors>=1) {
    // do nothing
   } else
    amrex::Error("blobdata or num_colors invalid");

  } else if ((project_option==0)||
             (project_option==1)) {
   // do nothing
  } else
   amrex::Error("project_option invalid22");
  
  if (beta!=0.0)
   amrex::Error("beta invalid");

  operation_flag=4;

 } else if (interp_option==2) {//unew^{f}=unew^{f}+beta*diffuse_register^{c->f}

  if (idx_velcell>=0) {

   if ((idx_velcell==CURRENT_CELL_VEL_MF)||
       (idx_velcell==DELTA_CELL_VEL_MF))
    amrex::Error("idx_velcell collision");

   primary_vel_data=idx_velcell;  // increment
   secondary_vel_data=CURRENT_CELL_VEL_MF; 
  } else
   amrex::Error("idx_velcell invalid");

  if (project_option==3) {  // viscosity
   operation_flag=5;
  } else
   amrex::Error("project_option invalid23");

  if ((beta!=1.0)&&(beta!=-1.0))
   amrex::Error("beta invalid");

 } else if (interp_option==3) {//unew^{c,f -> f} in fluid  (=usolid in solid)

  amrex::Error("interp_option==3 not used");

  // interp_option=4 unew^{f} = 
  //   (i) unew^{f} in incompressible non-solid regions
  //   (ii) u^{f,save} + (unew^{c}-u^{c,save})^{c->f} in spectral regions or
  //        compressible regions.
  //   (iii) usolid in solid regions
 } else if (interp_option==4) {

  debug_ngrow(ADVECT_REGISTER_MF,1,2000);
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
   debug_ngrow(ADVECT_REGISTER_FACE_MF+dir,0,111);

  if ((project_option==0)||
      (project_option==1)) {
   operation_flag=11;
  } else
   amrex::Error("project_option invalid25");

  if (idx_velcell==-1) {
   primary_vel_data=DELTA_CELL_VEL_MF; 
   secondary_vel_data=CURRENT_CELL_VEL_MF; 
  } else
   amrex::Error("idx_velcell invalid");
   
 } else
  amrex::Error("interp_option invalid");

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

 int nsolve=1;

 bool use_tiling=ns_tiling;

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 resize_levelsetLO(2,LEVELPC_MF);
 debug_ngrow(LEVELPC_MF,2,110);
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  debug_ngrow(FACE_VAR_MF+dir,0,111);
 }

 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
  debug_ngrow(FSI_GHOST_MAC_MF+data_dir,0,112);
 }
 VOF_Recon_resize(1,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,1,120);

 resize_maskfiner(1,MASKCOEF_MF);
 resize_mask_nbr(1);

 int fluxvel_index=0;
 int fluxden_index=AMREX_SPACEDIM;

 const Real* dx = geom.CellSize();

 const Box& domain = geom.Domain();
 const int* domlo = domain.loVect();
 const int* domhi = domain.hiVect();

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
  debug_ngrow(FACE_VAR_MF+dir,0,122);

 allocate_flux_register(operation_flag);
 if (localMF[SEM_FLUXREG_MF]->nComp()!=AMREX_SPACEDIM)
  amrex::Error("localMF[SEM_FLUXREG_MF]->nComp() invalid2");

 if (num_colors==0) {
   // placeholders
  levelcolor=localMF[CURRENT_CELL_VEL_MF];
  leveltype=localMF[CURRENT_CELL_VEL_MF];
 } else if (num_colors>0) {
  levelcolor=localMF[COLOR_MF];
  leveltype=localMF[TYPE_MF];
 } else
  amrex::Error("num_colors invalid");

 if (levelcolor->nGrow()<1)
  amrex::Error("levelcolor->nGrow()<1");
 if (leveltype->nGrow()<1)
  amrex::Error("leveltyoe->nGrow()<1");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  new_localMF(AMRSYNC_VEL_MF+dir,1,0,dir);
  setVal_localMF(AMRSYNC_VEL_MF+dir,1.0e+40,0,1,0);
  if (localMF[AREA_MF+dir]->boxArray()!=
      localMF[AMRSYNC_VEL_MF+dir]->boxArray())
   amrex::Error("AMRSYNC_VEL boxarray does not match");
 }

 if ((interp_option==0)||
     (interp_option==2)||
     (interp_option==3)||
     (interp_option==4)) {

  if (level<finest_level) {
   avgDown_and_Copy_localMF(
     primary_vel_data,  // idx_den_MF
     primary_vel_data,  // idx_vel_MF
     AMRSYNC_VEL_MF,
     operation_flag);
  } else if (level==finest_level) {
   // do nothing
  } else
   amrex::Error("level invalid14");

  if ((level>=1)&&(level<=finest_level)) {
   interp_and_Copy_localMF(
     primary_vel_data,  // idx_den_MF
     primary_vel_data,  // idx_vel_MF
     AMRSYNC_VEL_MF,
     operation_flag);
  } else if (level==0) {
   // do nothing
  } else
   amrex::Error("level invalid15");

  // unew^{f} = unew^{f} in fluid (=usolid in solid)
 } else if (interp_option==1) {
  // do nothing
 } else
  amrex::Error("interp_option invalid");

 for (int spectral_loop=0;spectral_loop<end_spectral_loop();spectral_loop++) {

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

   MultiFab& Umac_new=get_new_data(Umac_Type+dir,slab_step+1);

   MultiFab* Umac_old;
   MultiFab* U_old;

   if (interp_option==4) { // operation_flag==11
    Umac_old=localMF[ADVECT_REGISTER_FACE_MF+dir];
    U_old=localMF[ADVECT_REGISTER_MF];
   } else if ((interp_option==0)|| //operation_flag==4
              (interp_option==1)|| //operation_flag==4
              (interp_option==2)|| //operation_flag==5
              (interp_option==3)) {//operation_flag==10

    int ncomp_MAC=Umac_new.nComp();
    Umac_old=getStateMAC(Umac_Type,0,dir,0,ncomp_MAC,cur_time_slab); 
    if (Umac_old->boxArray()==Umac_new.boxArray()) {
     // do nothing
    } else
     amrex::Error("Umac_old->boxArray() invalid");

    U_old=localMF[CURRENT_CELL_VEL_MF];
   } else
    amrex::Error("interp_option invalid");

    // if low order, then nothing done when tileloop==1
   for (int tileloop=0;tileloop<=1;tileloop++) {

    if (thread_class::nthreads<1)
     amrex::Error("thread_class::nthreads invalid");
    thread_class::init_d_numPts(
    localMF[CURRENT_CELL_VEL_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
    for (MFIter mfi(*localMF[CURRENT_CELL_VEL_MF],use_tiling); 
  	            mfi.isValid(); ++mfi) {
      BL_ASSERT(grids[mfi.index()] == mfi.validbox());
      int gridno=mfi.index();
      const Box& tilegrid = mfi.tilebox();
      const Box& fabgrid = grids[gridno];
      const int* tilelo=tilegrid.loVect();
      const int* tilehi=tilegrid.hiVect();
      const int* fablo=fabgrid.loVect();
      const int* fabhi=fabgrid.hiVect();
      int bfact=parent->Space_blockingFactor(level);
      int bfact_c=bfact;
      int bfact_f=bfact;
      if (level>0)
       bfact_c=parent->Space_blockingFactor(level-1);
      if (level<finest_level)
       bfact_f=parent->Space_blockingFactor(level+1);

      const Real* xlo = grid_loc[gridno].lo();
    
      FArrayBox& xface=(*localMF[FACE_VAR_MF+dir])[mfi];  

      FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi];  

      FArrayBox& xvel=Umac_new[mfi];

      FArrayBox& xgp=(*Umac_old)[mfi];

      FArrayBox& xp=(*localMF[AMRSYNC_VEL_MF+dir])[mfi];

      FArrayBox& pres=(*U_old)[mfi];

       // FSI_GHOST_MAC_MF is initialized in 
       //  init_FSI_GHOST_MAC_MF_ALL(caller_id)
      FArrayBox& solfab=(*localMF[FSI_GHOST_MAC_MF+dir])[mfi];

      FArrayBox& primary_velfab=(*localMF[primary_vel_data])[mfi];
      FArrayBox& secondary_velfab=(*localMF[secondary_vel_data])[mfi];

      FArrayBox& lsfab=(*localMF[LEVELPC_MF])[mfi];

      FArrayBox& colorfab=(*levelcolor)[mfi];
      FArrayBox& typefab=(*leveltype)[mfi];

      // mask=tag if not covered by level+1 or outside the domain.
      FArrayBox& maskcoeffab=(*localMF[MASKCOEF_MF])[mfi];
      FArrayBox& maskSEMfab=(*localMF[MASKSEM_MF])[mfi];

      FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi];
      FArrayBox& semfluxfab=(*localMF[SEM_FLUXREG_MF])[mfi];
      int ncfluxreg=semfluxfab.nComp();

      Vector<int> velbc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);

      int rzflag=0;
      if (geom.IsRZ())
       rzflag=1;
      else if (geom.IsCartesian())
       rzflag=0;
      else if (geom.IsCYLINDRICAL())
       rzflag=3;
      else
       amrex::Error("CoordSys bust 20");

      int energyflag=0;
      int local_enable_spectral=enable_spectral;
      int simple_AMR_BC_flag=0;
      int ncomp_xp=1;
      int ncomp_xgp=1;
      int ncomp_mgoni=AMREX_SPACEDIM;

      int tid_current=ns_thread();
      if ((tid_current<0)||(tid_current>=thread_class::nthreads))
       amrex::Error("tid_current invalid");
      thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

      // in increment_face_velocity
      fort_cell_to_mac(
       &ncomp_mgoni,
       &ncomp_xp,
       &ncomp_xgp,
       &simple_AMR_BC_flag,
       &nsolve,
       &tileloop,
       &dir,
       &operation_flag, // 3,4,5, or 10
       &energyflag,
       &beta,
       &visc_coef,
       &interp_vel_increment_from_cell,
       filter_velocity.dataPtr(),
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
       &ignore_div_up,
       &pforce_index,
       &faceden_index,
       &icemask_index,
       &massface_index,
       &vofface_index,
       &ncphys,
       override_density.dataPtr(),
       constant_density_all_time.dataPtr(),
       velbc.dataPtr(),  // presbc
       velbc.dataPtr(),  
       &slab_step,
       &dt_slab,
       &cur_time_slab, 
       xlo,dx,
       &spectral_loop,
       &ncfluxreg,
       semfluxfab.dataPtr(),
       ARLIM(semfluxfab.loVect()),ARLIM(semfluxfab.hiVect()),
       maskfab.dataPtr(), // mask=1.0 at interior fine bc ghost cells
       ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
        // mask=tag if not covered by level+1 or outside the domain.
       maskcoeffab.dataPtr(),
       ARLIM(maskcoeffab.loVect()),ARLIM(maskcoeffab.hiVect()),
       maskSEMfab.dataPtr(),
       ARLIM(maskSEMfab.loVect()),ARLIM(maskSEMfab.hiVect()),
       lsfab.dataPtr(),
       ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
       solfab.dataPtr(),
       ARLIM(solfab.loVect()),ARLIM(solfab.hiVect()), //FSI_GHOST_MAC_MF
       xface.dataPtr(),
       ARLIM(xface.loVect()),ARLIM(xface.hiVect()), //xcut
       xface.dataPtr(),
       ARLIM(xface.loVect()),ARLIM(xface.hiVect()),
       reconfab.dataPtr(),ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
       xgp.dataPtr(),ARLIM(xgp.loVect()),ARLIM(xgp.hiVect()), //holds Umac_old
       xp.dataPtr(),ARLIM(xp.loVect()),ARLIM(xp.hiVect()), //xp(holds AMRSYNC)
       xvel.dataPtr(),ARLIM(xvel.loVect()),ARLIM(xvel.hiVect()), 
       primary_velfab.dataPtr(),
       ARLIM(primary_velfab.loVect()),ARLIM(primary_velfab.hiVect()),
       pres.dataPtr(dir), // pres holds U_old
       ARLIM(pres.loVect()),ARLIM(pres.hiVect()),
       primary_velfab.dataPtr(), // den
       ARLIM(primary_velfab.loVect()),ARLIM(primary_velfab.hiVect()),
       secondary_velfab.dataPtr(), // mgoni
       ARLIM(secondary_velfab.loVect()),ARLIM(secondary_velfab.hiVect()),
       colorfab.dataPtr(),
       ARLIM(colorfab.loVect()),ARLIM(colorfab.hiVect()),
       typefab.dataPtr(),
       ARLIM(typefab.loVect()),ARLIM(typefab.hiVect()),
       tilelo,tilehi,
       fablo,fabhi,
       &bfact,&bfact_c,&bfact_f, 
       &level,&finest_level,
       &rzflag,
       domlo,domhi, 
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
       &project_option,
       &SEM_upwind,
       &SEM_advection_algorithm);
    } // mfi
} // omp
    ns_reconcile_d_num(134);
   } // tileloop

   if (interp_option==4) {
    // do nothing
   } else if ((interp_option==0)|| 
              (interp_option==1)||
              (interp_option==2)||
              (interp_option==3)) {
    delete Umac_old;
   } else
    amrex::Error("interp_option invalid");

  } // dir=0..sdim-1

  synchronize_flux_register(operation_flag,spectral_loop);
 } // spectral_loop

} // subroutine increment_face_velocity

// update the faceden_index component of FACE_VAR_MF
// called from make_physics_vars
void NavierStokes::density_TO_MAC(int project_option) {

 int operation_flag=9;

 int num_colors=0;
 Vector<Real> blob_array;
 blob_array.resize(1);
 int blob_array_size=blob_array.size();

 int finest_level = parent->finestLevel();
 int bfact=parent->Space_blockingFactor(level);
 int bfact_c=bfact;
 int bfact_f=bfact;
 if (level>0)
  bfact_c=parent->Space_blockingFactor(level-1);
 if (level<finest_level)
  bfact_f=parent->Space_blockingFactor(level+1);

 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

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

 bool use_tiling=ns_tiling;

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 resize_levelsetLO(2,LEVELPC_MF);
 debug_ngrow(LEVELPC_MF,2,110);
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  debug_ngrow(FACE_VAR_MF+dir,0,111);
   // allocated in make_physics_vars 
   // deallocated in make_physics_varsALL
  debug_ngrow(AMRSYNC_VEL_MF+dir,0,111); 
 }
 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
  debug_ngrow(FSI_GHOST_MAC_MF+data_dir,0,112);
 }
 VOF_Recon_resize(1,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,1,124);

 resize_maskfiner(1,MASKCOEF_MF);
 resize_mask_nbr(1);

 int fluxvel_index=0;
 int fluxden_index=AMREX_SPACEDIM;

 const Real* dx = geom.CellSize();

 const Box& domain = geom.Domain();
 const int* domlo = domain.loVect();
 const int* domhi = domain.hiVect();

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
  debug_ngrow(FACE_VAR_MF+dir,0,126);

 debug_ngrow(DEN_RECON_MF,1,127);
 if (localMF[DEN_RECON_MF]->nComp()!=nmat*num_state_material)
  amrex::Error("localMF[DEN_RECON_MF]->nComp() invalid");

 if ((projection_enable_spectral==1)||
     (projection_enable_spectral==2)) { //in:density_TO_MAC

  if (bfact>=2) {

   if (level<finest_level) {
    avgDown_and_Copy_localMF(
     DEN_RECON_MF,
     DEN_RECON_MF,
     AMRSYNC_VEL_MF,
     operation_flag);
   } else if (level==finest_level) {
    // do nothing
   } else
    amrex::Error("level invalid16");

   if ((level>=1)&&(level<=finest_level)) {
    interp_and_Copy_localMF(
     DEN_RECON_MF,
     DEN_RECON_MF,
     AMRSYNC_VEL_MF,
     operation_flag);
   } else if (level==0) {
    // do nothing
   } else
    amrex::Error("level invalid17");

   allocate_flux_register(operation_flag);
   if (localMF[SEM_FLUXREG_MF]->nComp()!=AMREX_SPACEDIM)
    amrex::Error("localMF[SEM_FLUXREG_MF]->nComp() invalid3");

   for (int spectral_loop=0;spectral_loop<end_spectral_loop();
        spectral_loop++) {
    for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
 
     for (int tileloop=0;tileloop<=1;tileloop++) {

      if (thread_class::nthreads<1)
       amrex::Error("thread_class::nthreads invalid");
      thread_class::init_d_numPts(
       localMF[SLOPE_RECON_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
      for (MFIter mfi(*localMF[SLOPE_RECON_MF],use_tiling); 
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
 
       FArrayBox& xface=(*localMF[FACE_VAR_MF+dir])[mfi];  
       FArrayBox& xp=(*localMF[AMRSYNC_VEL_MF+dir])[mfi];  

       FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi];  

       FArrayBox& xvel=(*localMF[FACE_VAR_MF+dir])[mfi];  
 
       FArrayBox& sol=(*localMF[FSI_GHOST_MAC_MF+dir])[mfi];
       FArrayBox& cellvelfab=(*localMF[DEN_RECON_MF])[mfi];

       FArrayBox& lsfab=(*localMF[LEVELPC_MF])[mfi];

       // mask=tag if not covered by level+1 or outside the domain.
       FArrayBox& maskcoeffab=(*localMF[MASKCOEF_MF])[mfi];
       FArrayBox& maskSEMfab=(*localMF[MASKSEM_MF])[mfi];

       FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi];
       FArrayBox& semfluxfab=(*localMF[SEM_FLUXREG_MF])[mfi];
       int ncfluxreg=semfluxfab.nComp();

       Vector<int> velbc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);

       int scomp_den=(AMREX_SPACEDIM+1);
       Vector<int> denbc=getBCArray(State_Type,gridno,scomp_den,1);

       Real beta=0.0;
 
       int rzflag=0;
       if (geom.IsRZ())
        rzflag=1;
       else if (geom.IsCartesian())
        rzflag=0;
       else if (geom.IsCYLINDRICAL())
        rzflag=3;
       else
        amrex::Error("CoordSys bust 20");

       int energyflag=0;
       int local_enable_spectral=projection_enable_spectral;
       int simple_AMR_BC_flag=0;
       int ncomp_xp=1;
       int ncomp_xgp=1;
       int ncomp_mgoni=cellvelfab.nComp();
       int nsolve=1;

       int tid_current=ns_thread();
       if ((tid_current<0)||(tid_current>=thread_class::nthreads))
        amrex::Error("tid_current invalid");
       thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

        // in density_TO_MAC
       fort_cell_to_mac(
        &ncomp_mgoni,
        &ncomp_xp,
        &ncomp_xgp,
        &simple_AMR_BC_flag,
        &nsolve,
        &tileloop,
        &dir,
        &operation_flag, // 9
        &energyflag,
        &beta,
        &visc_coef,
        &interp_vel_increment_from_cell,
        filter_velocity.dataPtr(),
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
        &ignore_div_up,
        &pforce_index,
        &faceden_index,
        &icemask_index,
        &massface_index,
        &vofface_index,
        &ncphys,
        override_density.dataPtr(),
        constant_density_all_time.dataPtr(),
        denbc.dataPtr(),  // presbc
        velbc.dataPtr(),  
        &slab_step,
        &dt_slab,
        &cur_time_slab, 
        xlo,dx,
        &spectral_loop,
        &ncfluxreg,
        semfluxfab.dataPtr(),
        ARLIM(semfluxfab.loVect()),ARLIM(semfluxfab.hiVect()),
        maskfab.dataPtr(), // mask=1.0 at interior fine bc ghost cells
        ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
         // mask=tag if not covered by level+1 or outside the domain.
        maskcoeffab.dataPtr(),
        ARLIM(maskcoeffab.loVect()),ARLIM(maskcoeffab.hiVect()),
        maskSEMfab.dataPtr(),
        ARLIM(maskSEMfab.loVect()),ARLIM(maskSEMfab.hiVect()),
        lsfab.dataPtr(),
        ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
        sol.dataPtr(),ARLIM(sol.loVect()),ARLIM(sol.hiVect()),
        xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()), //xcut
        xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()),
        reconfab.dataPtr(),ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
        xvel.dataPtr(faceden_index),
        ARLIM(xvel.loVect()),ARLIM(xvel.hiVect()), //xgp
        xp.dataPtr(),ARLIM(xp.loVect()),ARLIM(xp.hiVect()), //xp
        xvel.dataPtr(faceden_index),               //xvel=1/rho (destination)
        ARLIM(xvel.loVect()),ARLIM(xvel.hiVect()), 
        cellvelfab.dataPtr(), // contains the density
        ARLIM(cellvelfab.loVect()),ARLIM(cellvelfab.hiVect()),
        cellvelfab.dataPtr(), // pres
        ARLIM(cellvelfab.loVect()),ARLIM(cellvelfab.hiVect()),
        cellvelfab.dataPtr(), // den (source)
        ARLIM(cellvelfab.loVect()),ARLIM(cellvelfab.hiVect()),
        cellvelfab.dataPtr(), // mgoni
        ARLIM(cellvelfab.loVect()),ARLIM(cellvelfab.hiVect()),
        cellvelfab.dataPtr(), // color
        ARLIM(cellvelfab.loVect()),ARLIM(cellvelfab.hiVect()),
        cellvelfab.dataPtr(), // type
        ARLIM(cellvelfab.loVect()),ARLIM(cellvelfab.hiVect()),
        tilelo,tilehi,
        fablo,fabhi,
        &bfact,&bfact_c,&bfact_f, 
        &level,
        &finest_level,
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
        &project_option,
        &SEM_upwind,
        &SEM_advection_algorithm);
      } // mfi
} // omp
      ns_reconcile_d_num(136);
     } // tileloop
    } // dir
    ParallelDescriptor::Barrier();

    synchronize_flux_register(operation_flag,spectral_loop);
   } // spectral_loop=0..end_spectral_loop()-1

  } else if (bfact==1) {
   // do nothing
  } else
   amrex::Error("bfact invalid");
   
 } else if ((projection_enable_spectral==0)||
            (projection_enable_spectral==3)) { //in:density_TO_MAC
  // do nothing
 } else
  amrex::Error("projection_enable_spectral invalid");

} // subroutine density_TO_MAC

// vel_or_disp=0 => interpolate mac velocity
// vel_or_disp=1 => interpolate mac displacement
// dest_idx==-1 => destination is the state data.
// dest_idx>=0  => destination is localMF[dest_idx]
void NavierStokes::VELMAC_TO_CELLALL(
  int use_VOF_weight,
  int vel_or_disp,
  int dest_idx) {

 int finest_level=parent->finestLevel();

 if (level==0) {
  // do nothing
 } else
  amrex::Error("expecting level==0 in VELMAC_TO_CELLALL");

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.VELMAC_TO_CELL(use_VOF_weight,vel_or_disp,dest_idx);
 }
 if (dest_idx==-1) {
  // do nothing
 } else if (dest_idx>=0) {
  Vector<int> scompBC_map;
  scompBC_map.resize(AMREX_SPACEDIM);

  if (vel_or_disp==0) { // velocity
   for (int dir=0;dir<AMREX_SPACEDIM;dir++)
    scompBC_map[dir]=dir;
   GetStateFromLocalALL(dest_idx,localMF[dest_idx]->nGrow(),0,
     AMREX_SPACEDIM,State_Type,scompBC_map);
  } else if (vel_or_disp==1) { // displacement

   if (NUM_TENSOR_TYPE==2*AMREX_SPACEDIM) {
    // do nothing
   } else
    amrex::Error("NUM_TENSOR_TYPE invalid");

   for (int dir=0;dir<AMREX_SPACEDIM;dir++)
    scompBC_map[dir]=NUM_TENSOR_TYPE+dir;

   PCINTERP_fill_bordersALL(dest_idx,localMF[dest_idx]->nGrow(),0,
     AMREX_SPACEDIM,Tensor_Type,scompBC_map);
  } else
   amrex::Error("vel_or_disp invalid");
 } else
  amrex::Error("dest_idx invalid");

} // end subroutine VELMAC_TO_CELLALL

// vel_or_disp=0 => interpolate mac velocity
// vel_or_disp=1 => interpolate mac displacement
// dest_idx==-1 => destination is the state data.
// dest_idx>=0  => destination is localMF[dest_idx]
void NavierStokes::VELMAC_TO_CELL(
  int use_VOF_weight,
  int vel_or_disp,
  int dest_idx) {
 
 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();
 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid");

 int nsolve=1;

 int fluxvel_index=0;
 int fluxden_index=AMREX_SPACEDIM;

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
  debug_ngrow(FACE_VAR_MF+dir,0,1111);

 debug_ngrow(VOLUME_MF,0,751);
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  if (localMF[AREA_MF+dir]->boxArray()!=
      localMF[FACE_VAR_MF+dir]->boxArray())
   amrex::Error("boxarrays do not match");
 }

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
  debug_ngrow(FSI_GHOST_MAC_MF+data_dir,0,112);
 }
 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
  if (localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0)
   amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0");
  if (localMF[FSI_GHOST_MAC_MF+data_dir]->nComp()!=nparts_def*AMREX_SPACEDIM)
   amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nComp() bad");
 }

 resize_levelsetLO(2,LEVELPC_MF);
 debug_ngrow(LEVELPC_MF,2,103);
 if (localMF[LEVELPC_MF]->nComp()!=nmat*(AMREX_SPACEDIM+1))
  amrex::Error("(localMF[LEVELPC_MF]->nComp()!=nmat*(AMREX_SPACEDIM+1))");

 VOF_Recon_resize(1,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,1,128);

 const Real* dx = geom.CellSize();

 resize_maskfiner(1,MASKCOEF_MF);
 resize_mask_nbr(1);
 resize_metrics(1);
 debug_ngrow(MASKCOEF_MF,1,253); // maskcoef=1 if not covered by finer level.
 debug_ngrow(MASK_NBR_MF,1,253); // mask_nbr=1 at fine-fine bc.

 int operation_flag=2;
 int MAC_state_idx=Umac_Type;

 int local_enable_spectral=enable_spectral;

 MultiFab* face_velocity[AMREX_SPACEDIM];
 MultiFab* dest_velocity=nullptr;

 if (vel_or_disp==0) { //velocity
  MAC_state_idx=Umac_Type;
  operation_flag=2;
 } else if (vel_or_disp==1) { //displacement
  MAC_state_idx=XDmac_Type;
  operation_flag=7;
  local_enable_spectral=0;
 } else 
  amrex::Error("vel_or_disp invalid");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  face_velocity[dir]=getStateMAC(
    MAC_state_idx,0,dir,0,nsolve,cur_time_slab);
 }

 if (dest_idx==-1) {
  if (vel_or_disp==0) {
   MultiFab& S_new=get_new_data(State_Type,slab_step+1);
   dest_velocity=&S_new;
  } else if (vel_or_disp==1) {
   amrex::Error("no state cell centered disp available");
  } else
   amrex::Error("vel_or_disp invalid");
 } else if (dest_idx>=0) {
  dest_velocity=localMF[dest_idx];
 } else
  amrex::Error("dest_idx invalid");

 if (dest_velocity->nComp()>=AMREX_SPACEDIM) {
  // do nothing
 } else 
  amrex::Error("dest_velocity->nComp() invalid");

 if (dest_velocity->nGrow()>=1) {
  // do nothing
 } else 
  amrex::Error("dest_velocity->nGrow() invalid");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  debug_ngrow(FACE_VAR_MF+dir,0,129);

  MultiFab& Umac_new=get_new_data(MAC_state_idx+dir,slab_step+1);
  int ncmac=Umac_new.nComp();

  if (ncmac!=nsolve) {
   std::cout << "nmat = " << nmat << '\n';
   std::cout << "ncmac = " << ncmac << '\n';
   amrex::Error("ncmac incorrect MAC_to_CELL");
  }
 }

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(localMF[SLOPE_RECON_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*localMF[SLOPE_RECON_MF],use_tiling); mfi.isValid(); ++mfi) {
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

   // mask=1.0 at interior fine bc ghost cells
  FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi];
   //1=not cov  0=cov
  FArrayBox& maskcoef = (*localMF[MASKCOEF_MF])[mfi];

  FArrayBox& maskSEMfab=(*localMF[MASKSEM_MF])[mfi];

  FArrayBox& xface=(*localMF[FACE_VAR_MF])[mfi];  
  FArrayBox& yface=(*localMF[FACE_VAR_MF+1])[mfi];  
  FArrayBox& zface=(*localMF[FACE_VAR_MF+AMREX_SPACEDIM-1])[mfi];  

  FArrayBox& veldest=(*dest_velocity)[mfi];
  FArrayBox& solxfab=(*localMF[FSI_GHOST_MAC_MF])[mfi];
  FArrayBox& solyfab=(*localMF[FSI_GHOST_MAC_MF+1])[mfi];
  FArrayBox& solzfab=(*localMF[FSI_GHOST_MAC_MF+AMREX_SPACEDIM-1])[mfi];
  FArrayBox& voffab=(*localMF[SLOPE_RECON_MF])[mfi];
  FArrayBox& levelpcfab=(*localMF[LEVELPC_MF])[mfi];
  FArrayBox& xvel=(*face_velocity[0])[mfi];
  FArrayBox& yvel=(*face_velocity[1])[mfi];
  FArrayBox& zvel=(*face_velocity[AMREX_SPACEDIM-1])[mfi];

  FArrayBox& volfab=(*localMF[VOLUME_MF])[mfi];
  FArrayBox& areax=(*localMF[AREA_MF])[mfi];
  FArrayBox& areay=(*localMF[AREA_MF+1])[mfi];
  FArrayBox& areaz=(*localMF[AREA_MF+AMREX_SPACEDIM-1])[mfi];

  Vector<int> velbc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);

  int energyflag=0;
  int project_option=0;
  int homflag=0; // default

  int ncomp_denold=volfab.nComp();
  int ncomp_veldest=veldest.nComp();
  int ncomp_dendest=veldest.nComp();

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // we are in VELMAC_TO_CELL
   // fort_mac_to_cell is declared in LEVELSET_3D.F90
  fort_mac_to_cell(
   &ns_time_order,
   &divu_outer_sweeps,
   &num_divu_outer_sweeps,
   &operation_flag, // operation_flag=2 (mac_vel -> cell_vel) or 7 (disp)
   &energyflag,
   temperature_primitive_variable.dataPtr(),
   constant_density_all_time.dataPtr(),
   &nmat,
   &nparts,
   &nparts_def,
   im_solid_map_ptr,
   added_weight.dataPtr(),
   &nten,
   &level,
   &finest_level,
   &project_option,
   &local_enable_spectral, //0 if interp displacement to CELLS.
   &fluxvel_index,
   &fluxden_index,
   &facevel_index,
   &facecut_index,
   &icefacecut_index,
   &curv_index,
   &conservative_tension_force,
   &conservative_div_uu,
   filter_velocity.dataPtr(),
   &ignore_div_up,
   &pforce_index,
   &faceden_index,
   &icemask_index,
   &massface_index,
   &vofface_index,
   &ncphys,
   velbc.dataPtr(),
   velbc.dataPtr(), // presbc
   &cur_time_slab,
   &slab_step,
   &dt_slab, 
   xlo,dx,
   tilelo,tilehi,
   fablo,fabhi,
   &bfact,
   xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()), // xp
   yface.dataPtr(),ARLIM(yface.loVect()),ARLIM(yface.hiVect()), // yp
   zface.dataPtr(),ARLIM(zface.loVect()),ARLIM(zface.hiVect()), // zp
   xvel.dataPtr(),ARLIM(xvel.loVect()),ARLIM(xvel.hiVect()), 
   yvel.dataPtr(),ARLIM(yvel.loVect()),ARLIM(yvel.hiVect()), 
   zvel.dataPtr(),ARLIM(zvel.loVect()),ARLIM(zvel.hiVect()), 
   xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()), 
   yface.dataPtr(),ARLIM(yface.loVect()),ARLIM(yface.hiVect()),
   zface.dataPtr(),ARLIM(zface.loVect()),ARLIM(zface.hiVect()),
   areax.dataPtr(),ARLIM(areax.loVect()),ARLIM(areax.hiVect()), //ax
   areay.dataPtr(),ARLIM(areay.loVect()),ARLIM(areay.hiVect()), //ay
   areaz.dataPtr(),ARLIM(areaz.loVect()),ARLIM(areaz.hiVect()), //az
   volfab.dataPtr(),ARLIM(volfab.loVect()),ARLIM(volfab.hiVect()),
   voffab.dataPtr(),ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()), // rhs
   veldest.dataPtr(),ARLIM(veldest.loVect()),ARLIM(veldest.hiVect()), 
   veldest.dataPtr(),ARLIM(veldest.loVect()),ARLIM(veldest.hiVect()), //dendest
   maskfab.dataPtr(), // 1=fine/fine  0=coarse/fine
   ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
   maskcoef.dataPtr(), // 1=not covered  0=covered
   ARLIM(maskcoef.loVect()),ARLIM(maskcoef.hiVect()),
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
   levelpcfab.dataPtr(),
   ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),//cterm
   levelpcfab.dataPtr(),
   ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),//pold
   volfab.dataPtr(),
   ARLIM(volfab.loVect()),ARLIM(volfab.hiVect()), // denold
   veldest.dataPtr(),ARLIM(veldest.loVect()),ARLIM(veldest.hiVect()), //ustar
   voffab.dataPtr(),ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),//recon
   voffab.dataPtr(),ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),//mdot
   voffab.dataPtr(),ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),//maskdivres
   voffab.dataPtr(),ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),//maskres
   &SDC_outer_sweeps,
   &homflag,
   &use_VOF_weight,
   &nsolve,
   &ncomp_denold,
   &ncomp_veldest,
   &ncomp_dendest,
   &SEM_advection_algorithm);
 }   // mfi
} // omp
 ns_reconcile_d_num(137);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
  delete face_velocity[dir]; 

} // subroutine VELMAC_TO_CELL


//if temperature_primitive_var==0,
// add beta * (1/cv) * (u dot u/2) to temp
void NavierStokes::increment_KE(Real beta) {
 
 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();
 int nmat=num_materials;

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 if ((beta!=-1.0)&&(beta!=1.0))
  amrex::Error("beta invalid");

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,1,253); // maskcoef=1 if not covered by finer level.

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

 int ncomp_state=(AMREX_SPACEDIM+1)+nmat*num_state_material+
  nmat*ngeom_raw+1;
 if (ncomp_state!=S_new.nComp()) 
  amrex::Error("ncomp_state!=S_new.nComp()"); 

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(S_new.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(S_new,use_tiling); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  int gridno=mfi.index();
  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();

   //1=not cov  0=cov
  FArrayBox& maskcoef = (*localMF[MASKCOEF_MF])[mfi];
  FArrayBox& sfab=S_new[mfi];

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  FORT_INC_TEMP(
   &beta,
   temperature_primitive_variable.dataPtr(),
   &nmat,
   &level,
   &finest_level,
   &ncomp_state,
   tilelo,tilehi,
   fablo,fabhi,
   sfab.dataPtr(),ARLIM(sfab.loVect()),ARLIM(sfab.hiVect()), 
   maskcoef.dataPtr(), // 1=not covered  0=covered
   ARLIM(maskcoef.loVect()),ARLIM(maskcoef.hiVect()));
 }   // mfi
} // omp
 ns_reconcile_d_num(138);

} // subroutine increment_KE


// do_alloc=1 => allocate variable
// do_alloc=0 => variable already allocated
void NavierStokes::init_gradu_tensorALL(
 int im_tensor,  // =-1 if input is velocity, >=0 if input is displacement. 
 int idx_vel, //source velocity or displacement; allocated if do_alloc==1,
              //deleted if do_alloc==1.
 int do_alloc,
 int idx_cell,
 int idx_face,
 int idx_elastic_flux, //xflux,yflux,zflux
 int simple_AMR_BC_flag_viscosity) {

 int finest_level=parent->finestLevel();

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid init_gradu_tensorALL");

 int nsolve=AMREX_SPACEDIM;

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);

  if (do_alloc==1) {

   if (im_tensor==-1) {
     // ngrow,scomp,ncomp
    ns_level.getState_localMF(idx_vel,1,0,nsolve,cur_time_slab);
    if (idx_elastic_flux==-1) {
     // do nothing
    } else
     amrex::Error("idx_elastic_flux invalid");
   } else if ((im_tensor>=0)&&(im_tensor<num_materials)) {
    if (idx_elastic_flux>=0) {
     // do nothing
    } else
     amrex::Error("idx_elastic_flux invalid");
    if (enable_spectral==0) {
     int elastic_partid=-1;
     for (int i=0;i<im_elastic_map.size();i++) {
      if (im_elastic_map[i]==im_tensor)
       elastic_partid=i;
     }

     if ((elastic_partid>=0)&&(elastic_partid<im_elastic_map.size())) {
      int scomp=num_materials_viscoelastic*NUM_TENSOR_TYPE;
      ns_level.getStateTensor_localMF(idx_vel,1,scomp,
		  AMREX_SPACEDIM,cur_time_slab);
     } else
      amrex::Error("elastic_partid invalid");
    } else
     amrex::Error("enable_spectral invalid");
   } else
    amrex::Error("im_tensor invalid");
      
  } else if (do_alloc==0) {
   if ((im_tensor==-1)&&(idx_elastic_flux==-1)) {
    // do nothing
   } else
    amrex::Error("im_tensor and idx_elastic_flux = -1 if do_alloc==0");
  } else
   amrex::Error("do_alloc invalid");

  ns_level.debug_ngrow(idx_vel,1,945);
  if (ns_level.localMF[idx_vel]->nComp()<nsolve)
   amrex::Error("ns_level.localMF[idx_vel]->nComp() invalid");

//ux,vx,wx,uy,vy,wy,uz,vz,wz
  int homflag=0;
  ns_level.init_gradu_tensor(
    im_tensor, // im_tensor=-1 => velocity   0<=im_tensor<=nmat-1 => xdisp
    homflag,
    idx_vel,
    idx_cell,
    idx_face,
    idx_elastic_flux,
    simple_AMR_BC_flag_viscosity);
 } // ilev=finest_level ... level

   // grad U.
 int ntensor=AMREX_SPACEDIM*AMREX_SPACEDIM; 

//ux,vx,wx,uy,vy,wy,uz,vz,wz
 int irow=0;
 int icol=0;
 for (int i=0;i<ntensor;i++) {
  Vector<int> scompBC_map;
   // desc_lstGHOST.setComponent(Tensor_Type, ...
   // "set_tensor_bc", tensor_pc_interp 
   // FORT_EXTRAPFILL
   // (i.e. the coarse/fine BC and physical BC will be low order)

  int scomp_extrap=0;
  if ((irow==0)&&(icol==0)) { //11
   scomp_extrap=0;
  } else if ((irow==1)&&(icol==0)) { //12
   scomp_extrap=1;
  } else if ((irow==2)&&(icol==0)&&(AMREX_SPACEDIM==3)) { //13
   scomp_extrap=4;
  } else if ((irow==0)&&(icol==1)) { //21
   scomp_extrap=1;
  } else if ((irow==1)&&(icol==1)) { //22
   scomp_extrap=2;
  } else if ((irow==2)&&(icol==1)&&(AMREX_SPACEDIM==3)) { //23
   scomp_extrap=5;
  } else if ((irow==0)&&(icol==2)&&(AMREX_SPACEDIM==3)) { //31
   scomp_extrap=4;
  } else if ((irow==1)&&(icol==2)&&(AMREX_SPACEDIM==3)) { //32
   scomp_extrap=5;
  } else if ((irow==2)&&(icol==2)&&(AMREX_SPACEDIM==3)) { //33
   scomp_extrap=3;
  } else
   amrex::Error("irow or icol invalid");

  scompBC_map.resize(1);
  scompBC_map[0]=scomp_extrap;
   // idx,ngrow,scomp,ncomp,index,scompBC_map
  if (im_tensor==-1) { //input is velocity
   PCINTERP_fill_bordersALL(idx_cell,1,i,1,Tensor_Type,scompBC_map);

   // input is displacement
  } else if ((im_tensor>=0)&&(im_tensor<num_materials)) {
   PCINTERP_fill_bordersALL(idx_cell,1,i,1,Tensor_Type,scompBC_map);
  } else
   amrex::Error("im_tensor invalid");

   // 00,10,20,01,11,21,02,12,22
  icol++;
  if (icol>=AMREX_SPACEDIM) {
   irow++;
   icol=0;
  }
 } // i=0..ntensor-1

 if ((irow==AMREX_SPACEDIM)&&(icol==0)) {
  // do nothing
 } else
  amrex::Error("irow or icol invalid");

 if (do_alloc==1) {
  delete_array(idx_vel);
 } else if (do_alloc==0) {
  // do nothing
 } else
  amrex::Error("do_alloc invalid");

 delete_array(LSTENSOR_MF);
 delete_array(MASKSOLIDTENSOR_MF);

} // subroutine init_gradu_tensorALL

//called from init_gradu_tensor
//ux,vx,wx,uy,vy,wy,uz,vz,wz
// itensor_iter=0 face grad U
// itensor_iter=1 cell grad U
void NavierStokes::doit_gradu_tensor(
 int im_tensor, //-1, or, 0..nmat-1
 int homflag,int idx_vel,
 int idx_cell,int idx_face,int spectral_loop,int itensor_iter,
 MultiFab* mask3,
 int simple_AMR_BC_flag_viscosity) {

 int finest_level = parent->finestLevel();

 int elastic_partid=-1;

 if (im_tensor==-1) {
  // check nothing (gradient of velocity)
 } else if ((im_tensor>=0)&&(im_tensor<num_materials)) {
  // sanity checks for gradient of displacement vector.
  if (homflag==0) {
   // do nothing
  } else
   amrex::Error("expecting homflag=0 for grad displacement vector");
  if (enable_spectral==0) {
   elastic_partid=-1;
   for (int i=0;i<im_elastic_map.size();i++) {
    if (im_elastic_map[i]==im_tensor)
     elastic_partid=i;
   }
   if ((elastic_partid>=0)&&(elastic_partid<im_elastic_map.size())) {
    // do nothing
   } else
    amrex::Error("elastic_partid invalid");
  } else
   amrex::Error("expecting enable_spectral==0 for gradient of xdisp");
 } else
  amrex::Error("im_tensor invalid");

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid doit_gradu_tensor");

 int bfact=parent->Space_blockingFactor(level);
 int bfact_c=bfact;
 int bfact_f=bfact;
 if (level>0)
  bfact_c=parent->Space_blockingFactor(level-1);
 if (level<finest_level)
  bfact_f=parent->Space_blockingFactor(level+1);
 
 bool use_tiling=ns_tiling;

 int nmat=num_materials;

 debug_ngrow(idx_vel,1,845);
 if (localMF[idx_vel]->nComp()!=AMREX_SPACEDIM)
  amrex::Error("localMF[idx_vel] ncomp invalid");

 int ntensor=AMREX_SPACEDIM*AMREX_SPACEDIM; 

 debug_ngrow(LSTENSOR_MF,1,845);
 debug_ngrow(MASKSOLIDTENSOR_MF,1,845);
 if (localMF[LSTENSOR_MF]->nComp()!=AMREX_SPACEDIM)
  amrex::Error("lstensor has invalid ncomp");
 if (localMF[MASKSOLIDTENSOR_MF]->nComp()!=AMREX_SPACEDIM)
  amrex::Error("masksolidtensor has invalid ncomp");

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

 resize_maskfiner(1,MASKCOEF_MF);
 resize_mask_nbr(1);
 resize_metrics(1);
 VOF_Recon_resize(2,SLOPE_RECON_MF);
 resize_levelsetLO(2,LEVELPC_MF);

 debug_ngrow(MASKCOEF_MF,1,845);
 debug_ngrow(VOLUME_MF,1,845);
 debug_ngrow(SLOPE_RECON_MF,2,841);
 debug_ngrow(MASKSEM_MF,1,841);
 debug_ngrow(LEVELPC_MF,2,120);
 if (localMF[LEVELPC_MF]->nComp()!=nmat*(1+AMREX_SPACEDIM))
  amrex::Error("levelpc mf has incorrect ncomp");
 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
  debug_ngrow(FSI_GHOST_MAC_MF+data_dir,0,112);
 }
 if (localMF[SLOPE_RECON_MF]->nComp()!=nmat*ngeom_recon)
  amrex::Error("slope recon mf has incorrect ncomp");

 int rzflag=0;
 if (geom.IsRZ())
  rzflag=1;
 else if (geom.IsCartesian())
  rzflag=0;
 else if (geom.IsCYLINDRICAL())
  rzflag=3;
 else
  amrex::Error("CoordSys bust 20");

 const Real* dx = geom.CellSize();

 MultiFab* sem_flux_mf;

 if (itensor_iter==0) {  // compute grad U

  sem_flux_mf=localMF[SEM_FLUXREG_MF];
  if (sem_flux_mf->nComp()!=ntensor)
   amrex::Error("sem fluxreg mf has invalid ncomp");

  // interp grad U from MAC grid to CELL grid. (sem_flux_mf not used)
 } else if (itensor_iter==1) {  
  sem_flux_mf=localMF[idx_vel];
 } else
  amrex::Error("itensor_iter invalid");

 for (int dir=1;dir<=AMREX_SPACEDIM;dir++) {

  MultiFab* amrsync_vel_mf;
  if (simple_AMR_BC_flag_viscosity==0) {
   amrsync_vel_mf=localMF[AMRSYNC_PRES_MF+dir-1];
   if (amrsync_vel_mf->nComp()!=AMREX_SPACEDIM)
    amrex::Error("amrsync_vel_mf->nComp() invalid29");
   if (amrsync_vel_mf->boxArray()!=
       localMF[AREA_MF+dir-1]->boxArray())
    amrex::Error("amrsync_vel_mf: boxarrays do not match");
  } else if (simple_AMR_BC_flag_viscosity==1) {
   amrsync_vel_mf=localMF[idx_vel];
  } else
   amrex::Error("simple_AMR_BC_flag_viscosity invalid");

  // tileloop=0 low order
  // tileloop=1 high order
  for (int tileloop=0;tileloop<=1;tileloop++) {

   if (thread_class::nthreads<1)
    amrex::Error("thread_class::nthreads invalid");
   thread_class::init_d_numPts(localMF[idx_vel]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
   for (MFIter mfi(*localMF[idx_vel],use_tiling); mfi.isValid(); ++mfi) {
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

    FArrayBox& velfab=(*localMF[idx_vel])[mfi];

    FArrayBox& solidxfab=(*localMF[FSI_GHOST_MAC_MF])[mfi];
    FArrayBox& solidyfab=(*localMF[FSI_GHOST_MAC_MF+1])[mfi];
    FArrayBox& solidzfab=(*localMF[FSI_GHOST_MAC_MF+AMREX_SPACEDIM-1])[mfi];

    FArrayBox& levelpcfab=(*localMF[LEVELPC_MF])[mfi];
    FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi];  

    // mask0=tag if not covered by level+1 or outside the domain.
    FArrayBox& mask0fab=(*localMF[MASKCOEF_MF])[mfi];

    // mask3=tag at exterior fine/fine border.
    // mask3=1-tag at other exterior boundaries.
    FArrayBox& mask3fab=(*mask3)[mfi];
    FArrayBox& maskSEMfab=(*localMF[MASKSEM_MF])[mfi];

    FArrayBox& tensor_data=(*localMF[idx_face])[mfi];
    FArrayBox& cell_tensor_data=(*localMF[idx_cell])[mfi];
    FArrayBox& mask_tensor_data=(*localMF[MASKSOLIDTENSOR_MF])[mfi];
    FArrayBox& faceLS=(*localMF[LSTENSOR_MF])[mfi];

    FArrayBox& amrsyncfab=(*amrsync_vel_mf)[mfi];

    FArrayBox& semfluxfab=(*sem_flux_mf)[mfi];
    int ncfluxreg=semfluxfab.nComp();
    int local_enable_spectral=enable_spectral;

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // fort_face_gradients is declared in GODUNOV_3D.F90
    fort_face_gradients(
     &im_tensor,
     &elastic_partid,
     im_elastic_map.dataPtr(),
     &ns_time_order,
     &divu_outer_sweeps,
     &num_divu_outer_sweeps,
     &SDC_outer_sweeps,
     &tileloop,
     &dir,
     &slab_step,
     &itensor_iter,
     &cur_time_slab,
     temperature_primitive_variable.dataPtr(),
     &local_enable_spectral,
     velbc.dataPtr(),
     &spectral_loop,
     &ncfluxreg,
     semfluxfab.dataPtr(),
     ARLIM(semfluxfab.loVect()),ARLIM(semfluxfab.hiVect()),
     amrsyncfab.dataPtr(),
     ARLIM(amrsyncfab.loVect()),ARLIM(amrsyncfab.hiVect()),
      // mask0=1 if not covered by finer level or outside domain.
     mask0fab.dataPtr(),ARLIM(mask0fab.loVect()),ARLIM(mask0fab.hiVect()),
      // fine/fine bc ?
     mask3fab.dataPtr(),ARLIM(mask3fab.loVect()),ARLIM(mask3fab.hiVect()),
     maskSEMfab.dataPtr(),
     ARLIM(maskSEMfab.loVect()),ARLIM(maskSEMfab.hiVect()),
     faceLS.dataPtr(),
     ARLIM(faceLS.loVect()),ARLIM(faceLS.hiVect()),
     mask_tensor_data.dataPtr(),
     ARLIM(mask_tensor_data.loVect()),ARLIM(mask_tensor_data.hiVect()),
     tensor_data.dataPtr(),
     ARLIM(tensor_data.loVect()),ARLIM(tensor_data.hiVect()),
     cell_tensor_data.dataPtr(),
     ARLIM(cell_tensor_data.loVect()),ARLIM(cell_tensor_data.hiVect()),
     velfab.dataPtr(),ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
     solidxfab.dataPtr(),
     ARLIM(solidxfab.loVect()),ARLIM(solidxfab.hiVect()),
     solidyfab.dataPtr(),
     ARLIM(solidyfab.loVect()),ARLIM(solidyfab.hiVect()),
     solidzfab.dataPtr(),
     ARLIM(solidzfab.loVect()),ARLIM(solidzfab.hiVect()),
     levelpcfab.dataPtr(),
     ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),
     reconfab.dataPtr(),
     ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
     xlo,dx,
     &rzflag,
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,&bfact_c,&bfact_f,
     &level,
     &finest_level,
     &nmat,
     &nparts,
     &nparts_def,
     im_solid_map_ptr,
     &homflag,
     &ntensor,
     &SEM_upwind,
     &SEM_advection_algorithm,
     &simple_AMR_BC_flag_viscosity);
   } // mfi
} // omp
   ns_reconcile_d_num(139);
  } // tileloop=0..1
 } // dir

 int datatype=0;

 int im=0;

 for (int sc=0;sc<ntensor;sc++) {

   int dir=0;
   if ((sc>=0)&&(sc<AMREX_SPACEDIM)) { // ux,vx,wx
    dir=0;
   } else if (sc<2*AMREX_SPACEDIM) {  // uy,vy,wy
    dir=1;
   } else if ((sc<AMREX_SPACEDIM*AMREX_SPACEDIM)&&
              (AMREX_SPACEDIM==3)) { // uz,vz,wz
    dir=AMREX_SPACEDIM-1;
   } else
    amrex::Error("sc invalid");

   int sc_mat=im*ntensor+sc;

   if (itensor_iter==1) { // cell grad U
    localMF[idx_cell]->FillBoundary(sc_mat,1,geom.periodicity()); 
   } else if (itensor_iter==0) { // face grad U
    FillBoundaryTENSOR(localMF[idx_face],sc_mat,dir);
   } else
    amrex::Error("itensor_iter invalid");

 } // sc

 if (itensor_iter==0) {  // face grad U

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   FillBoundaryTENSOR(localMF[LSTENSOR_MF],dir,dir);
   FillBoundaryTENSOR(localMF[MASKSOLIDTENSOR_MF],dir,dir);
  } // dir

 } else if (itensor_iter==1) { //cell grad U
  // do nothing
 } else
  amrex::Error("itensor_iter invalid");
 
 if (check_nan==1) {

  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "in doit_gradu_tensor \n";
   std::cout << "homflag= " << homflag << '\n';
   std::cout << "spectral_loop= " << spectral_loop << '\n';
   std::cout << "itensor_iter= " << itensor_iter << '\n';
  }
  std::fflush(NULL);
  check_for_NAN(localMF[idx_vel],100);

  if (itensor_iter==1) {  // cell grad U
   datatype=2;
   check_for_NAN_TENSOR(datatype,localMF[idx_cell],101);
  } else if (itensor_iter==0) { // face grad U
   datatype=1;
   check_for_NAN_TENSOR(datatype,localMF[idx_face],102);
   for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
    check_for_NAN_TENSOR_base(datatype,localMF[LSTENSOR_MF],dir,dir,103);
    check_for_NAN_TENSOR_base(datatype,localMF[MASKSOLIDTENSOR_MF],dir,dir,104);
   }
  } else
   amrex::Error("itensor_iter invalid");

 } else if (check_nan==0) {
  // do nothing
 } else {
  amrex::Error("check_nan invalid");
 } 

} // subroutine doit_gradu_tensor

void NavierStokes::FillBoundaryTENSOR(
 MultiFab* mf,int sc,int dir) {

 if ((sc<0)||(sc>=mf->nComp()))
  amrex::Error("sc out of range");

 IndexType mactyp=IndexType::TheUMACType();
 if (dir==0) {
  mactyp=IndexType::TheUMACType();
 } else if (dir==1) {
  mactyp=IndexType::TheVMACType();
 } else if ((dir==2)&&(AMREX_SPACEDIM==3)) {
  mactyp=IndexType::TheWMACType();
 } else
  amrex::Error("dir invalid FillBoundaryTENSOR");

 MultiFab* macmf=new MultiFab(state[Umac_Type+dir].boxArray(),dmap,1,1,
    MFInfo().SetTag("macmf"),FArrayBoxFactory());

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(mf->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*mf,false); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  const Box& tilegrid = mfi.tilebox();

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  FArrayBox& cellfab=(*mf)[mfi];
  FArrayBox& macfab=(*macmf)[mfi];

  Box localbox(cellfab.box());
  localbox.growLo(dir,-1);
  FArrayBox localfab;
  localfab.resize(localbox,1);
  localfab.copy(cellfab,localbox,sc,localbox,0,1);
  localfab.SetBoxType(mactyp);

  macfab.setVal(0.0);
  macfab.copy(localfab,localfab.box(),0,localfab.box(),0,1);
 } // mfi
} // omp
 ns_reconcile_d_num(140);

 macmf->FillBoundary(geom.periodicity()); 

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(mf->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*mf,false); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  const Box& tilegrid = mfi.tilebox();

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  FArrayBox& cellfab=(*mf)[mfi];
  FArrayBox& macfab=(*macmf)[mfi];

  Box localbox(cellfab.box());
  localbox.growLo(dir,-1);
  FArrayBox localfab;
  localfab.resize(localbox,1);
  localfab.SetBoxType(mactyp);
  localfab.copy(macfab,localfab.box(),0,localfab.box(),0,1);
  localfab.SetBoxType(IndexType::TheCellType());

  cellfab.copy(localfab,localfab.box(),0,localfab.box(),sc,1);
 } // mfi
} // omp
 ns_reconcile_d_num(141);

 delete macmf;

} // subroutine FillBoundaryTENSOR

// called from apply_pressure_grad, init_gradu_tensorALL
void NavierStokes::init_gradu_tensor(
 int im_tensor, //-1, or, 0 .. nmat-1
 int homflag,
 int idx_vel,
 int idx_cell,
 int idx_face,
 int idx_elastic_flux, //xflux,yflux,zflux
 int simple_AMR_BC_flag_viscosity) {

 if (im_tensor==-1) {
  // check nothing (gradient of velocity)
 } else if ((im_tensor>=0)&&(im_tensor<num_materials)) {
  // sanity checks for gradient of displacement vector.
  if (homflag==0) {
   // do nothing
  } else
   amrex::Error("expecting homflag=0 for grad displacement vector");
  if (enable_spectral==0) {
   int elastic_partid=-1;
   for (int i=0;i<im_elastic_map.size();i++) {
    if (im_elastic_map[i]==im_tensor)
     elastic_partid=i;
   }
   if ((elastic_partid>=0)&&(elastic_partid<im_elastic_map.size())) {
    // do nothing
   } else
    amrex::Error("elastic_partid invalid");
  } else
   amrex::Error("enable_spectral invalid");
 } else
  amrex::Error("im_tensor invalid");

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid init_gradu_tensor");

  // mask3=tag at exterior fine/fine border.
  // mask3=1-tag at other exterior boundaries.
 int clear_phys_boundary=3;
 Real tag=1.0;
 MultiFab* mask3=maskfiner(1,tag,clear_phys_boundary);  
 int operation_flag=6; // evaluate tensor values

 if ((localMF_grow[idx_face]>=0)||
     (localMF_grow[idx_cell]>=0)||
     (localMF_grow[LSTENSOR_MF]>=0)||
     (localMF_grow[MASKSOLIDTENSOR_MF]>=0))
  amrex::Error("tensor scratch variables not previously deleted");

 debug_ngrow(idx_vel,1,845);
 if (localMF[idx_vel]->nComp()!=AMREX_SPACEDIM)
  amrex::Error("localMF[idx_vel]->nComp() invalid in init_gradu_tensor");

 int ntensor=AMREX_SPACEDIM*AMREX_SPACEDIM; 

 int nmasksolid=AMREX_SPACEDIM;

 new_localMF(idx_face,ntensor,1,-1);
 new_localMF(idx_cell,ntensor,1,-1);

  // x,y,z 
 new_localMF(MASKSOLIDTENSOR_MF,nmasksolid,1,-1);
  // x,y,z 
 new_localMF(LSTENSOR_MF,AMREX_SPACEDIM,1,-1);

 int itensor_iter=0; // tensor face (face grad U)

  // flux register is initialized to zero.
 allocate_flux_register(operation_flag);
 if (localMF[SEM_FLUXREG_MF]->nComp()!=ntensor)
  amrex::Error("localMF[SEM_FLUXREG_MF]->nComp() invalid4");

  // spectral_loop==0 
  //   low order face fluxes
  //   high order face fluxes in localMF[SEM_FLUXREG_MF]
  // spectral_loop==1
  //   average of flux values shared at face.
 int spectral_loop=0;
 for (spectral_loop=0;spectral_loop<end_spectral_loop();spectral_loop++) {
  doit_gradu_tensor(
   im_tensor,
   homflag,
   idx_vel,idx_cell,idx_face,spectral_loop,
   itensor_iter,mask3,
   simple_AMR_BC_flag_viscosity);
  synchronize_flux_register(operation_flag,spectral_loop);
 }

  // interpolate grad U from MAC grid to CELL grid.
 itensor_iter=1;  // tensor cell
 spectral_loop=0;
 doit_gradu_tensor(
   im_tensor,
   homflag,
   idx_vel,
   idx_cell,
   idx_face,
   spectral_loop,
   itensor_iter,mask3,
   simple_AMR_BC_flag_viscosity);

 if (im_tensor==-1) {

  // do nothing

 } else if ((im_tensor>=0)&&(im_tensor<num_materials)) {

  bool use_tiling=ns_tiling;

  int rzflag=0;
  if (geom.IsRZ())
   rzflag=1;
  else if (geom.IsCartesian())
   rzflag=0;
  else if (geom.IsCYLINDRICAL())
   rzflag=3;
  else
   amrex::Error("CoordSys bust 20");


  const Real* dx = geom.CellSize();
  int bfact=parent->Space_blockingFactor(level);

  int nmat=num_materials;
  int nden=nmat*num_state_material;

  resize_levelsetLO(2,LEVELPC_MF);
  debug_ngrow(LEVELPC_MF,2,110);

  for (int dir=1;dir<=AMREX_SPACEDIM;dir++) {

   debug_ngrow(idx_elastic_flux+dir-1,0,845);
 
   MultiFab* umac_new=&get_new_data(Umac_Type+dir-1,slab_step+1);
   MultiFab* grad_xdisp_mf=localMF[idx_elastic_flux+dir-1];
   if (grad_xdisp_mf->boxArray()==umac_new->boxArray()) {
    // do nothing
   } else
    amrex::Error("grad_xdisp_mf->boxArray()!=umac_new->boxArray()");

   if (grad_xdisp_mf->nComp()==AMREX_SPACEDIM*AMREX_SPACEDIM) {
    // do nothing
   } else
    amrex::Error("grad_xdisp_mf->nComp()!=AMREX_SPACEDIM^{2}"); 

   if (enable_spectral==0) {
    // do nothing
   } else
    amrex::Error("enable_spectral invalid");

   if (thread_class::nthreads<1)
    amrex::Error("thread_class::nthreads invalid");
   thread_class::init_d_numPts(localMF[LEVELPC_MF]->boxArray().d_numPts());
  
#ifdef _OPENMP
#pragma omp parallel
#endif
{
   for (MFIter mfi(*localMF[LEVELPC_MF],use_tiling); mfi.isValid(); ++mfi) {
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

    FArrayBox& velfab=(*localMF[idx_vel])[mfi];
    FArrayBox& levelpcfab=(*localMF[LEVELPC_MF])[mfi];

    FArrayBox& xflux=(*grad_xdisp_mf)[mfi];

    FArrayBox& xface=(*localMF[FACE_VAR_MF+dir-1])[mfi];

    FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi];  
 
    FArrayBox& tensor_data=(*localMF[idx_face])[mfi];
    FArrayBox& cell_tensor_data=(*localMF[idx_cell])[mfi];
    FArrayBox& mask_tensor_data=(*localMF[MASKSOLIDTENSOR_MF])[mfi];
    FArrayBox& faceLS=(*localMF[LSTENSOR_MF])[mfi];

    FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi];
    // maskcoef=tag if not covered by level+1 or outside the domain.
    FArrayBox& maskcoef_fab=(*localMF[MASKCOEF_MF])[mfi];

    FArrayBox& viscfab=(*localMF[CELL_VISC_MATERIAL_MF])[mfi];
    int ncomp_visc=viscfab.nComp();

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // declared in: GODUNOV_3D.F90
    //  visc_coef * viscface * (grad X + grad X^T)
    fort_crossterm_elastic(
     &ncomp_visc,
     &im_tensor, // 0..nmat-1
     &dir,
     viscfab.dataPtr(),ARLIM(viscfab.loVect()),ARLIM(viscfab.hiVect()),
     maskfab.dataPtr(), // 1=fine/fine  0=coarse/fine
     ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
     maskcoef_fab.dataPtr(), // maskcoef=tag if not cov by level+1 or outside.
     ARLIM(maskcoef_fab.loVect()),ARLIM(maskcoef_fab.hiVect()),
     faceLS.dataPtr(),
     ARLIM(faceLS.loVect()),ARLIM(faceLS.hiVect()),
     mask_tensor_data.dataPtr(),
     ARLIM(mask_tensor_data.loVect()),ARLIM(mask_tensor_data.hiVect()),
     tensor_data.dataPtr(),
     ARLIM(tensor_data.loVect()),ARLIM(tensor_data.hiVect()),
     cell_tensor_data.dataPtr(),
     ARLIM(cell_tensor_data.loVect()),ARLIM(cell_tensor_data.hiVect()),
     xlo,dx,
     &dt_slab,
     &cur_time_slab,
     velfab.dataPtr(),ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
     levelpcfab.dataPtr(),
     ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),
     xflux.dataPtr(),ARLIM(xflux.loVect()),ARLIM(xflux.hiVect()),
     xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()),
     reconfab.dataPtr(),ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
     &facevisc_index,
     &vofface_index,
     &massface_index,
     &ncphys,
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,
     &level,
     &rzflag,
     velbc.dataPtr(),
     &visc_coef,
     &nmat,
     &nden,
     &ntensor);
   } // mfi
} // omp
   ns_reconcile_d_num(142);

  } // dir=1..sdim

  int spectral_override=0;  // always do low order average down
  int ncomp_edge=-1;
  int caller_id=1;
   // idxMF,scomp,ncomp,start_dir,ndir,spectral_override,caller_id
  avgDownEdge_localMF(
    idx_elastic_flux,
    0,ncomp_edge,0,AMREX_SPACEDIM,spectral_override,caller_id); 

  for (int dir=1;dir<=AMREX_SPACEDIM;dir++) {
   localMF[idx_elastic_flux+dir-1]->
      FillBoundary(0,AMREX_SPACEDIM*AMREX_SPACEDIM,geom.periodicity());
  }

 } else
  amrex::Error("im_tensor invalid");

 delete mask3; 

} // subroutine init_gradu_tensor

  
// if projection:
// - dt*(grad p)*face_weight  
// face_weight=0 at embedded solid faces and on 
// the domain boundary where pressure has a Neumann BC.
void NavierStokes::apply_pressure_grad(
  int simple_AMR_BC_flag,
  int simple_AMR_BC_flag_viscosity,
  int homflag,
  int energyflag,
  int gp_mf,
  int pboth_mf,
  int project_option,int nsolve) {

 int finest_level = parent->finestLevel();

 int nmat=num_materials;

 if (project_option_momeqn(project_option)==1) {
  // do nothing
 } else if (project_option_momeqn(project_option)==0) {
  // do nothing
 } else
  amrex::Error("project_option_momeqn(project_option) invalid26");

 int bfact=parent->Space_blockingFactor(level);
 int bfact_c=bfact;
 int bfact_f=bfact;
 if (level>0)
  bfact_c=parent->Space_blockingFactor(level-1);
 if (level<finest_level)
  bfact_f=parent->Space_blockingFactor(level+1);
 
 bool use_tiling=ns_tiling;

 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 if ((energyflag!=0)&&(energyflag!=2))
  amrex::Error("energyflag invalid");

 debug_ngrow(pboth_mf,1,845);

 if (localMF[pboth_mf]->nComp()!=nsolve) {
  std::cout << "nsolve=" << nsolve << '\n';
  std::cout << "project_option= " << project_option << '\n';
  std::cout << "pboth ngrow= " << localMF[pboth_mf]->nGrow() << '\n';
  std::cout << "pboth ncomp= " << localMF[pboth_mf]->nComp() << '\n';
  amrex::Error("nsolve invalid28");
 }

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

 resize_maskfiner(1,MASKCOEF_MF);
 resize_mask_nbr(1);
 resize_metrics(1);
 VOF_Recon_resize(1,SLOPE_RECON_MF);

 debug_ngrow(MASKCOEF_MF,1,845);
 debug_ngrow(VOLUME_MF,1,845);
 debug_ngrow(SLOPE_RECON_MF,1,841);
 debug_ngrow(MASKSEM_MF,1,841);
 if (localMF[LEVELPC_MF]->nComp()!=nmat*(1+AMREX_SPACEDIM))
  amrex::Error("levelpc mf has incorrect ncomp");

 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
  if (localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0)
   amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0");
  if (localMF[FSI_GHOST_MAC_MF+data_dir]->nComp()!=nparts_def*AMREX_SPACEDIM)
   amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nComp() bad");
 }

 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
  debug_ngrow(FSI_GHOST_MAC_MF+data_dir,0,112);
 }
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
  debug_ngrow(FACE_VAR_MF+dir,0,122);

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

 const Real* dx = geom.CellSize();

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

  if (localMF[AMRSYNC_PRES_MF+dir]->nComp()!=nsolve)
   amrex::Error("localMF[AMRSYNC_PRES_MF+dir]->nComp() invalid29");
  if (localMF[AMRSYNC_PRES_MF+dir]->boxArray()!=
      localMF[AREA_MF+dir]->boxArray())
   amrex::Error("AMRSYNC_PRES boxarrays do not match");

  if (localMF[gp_mf+dir]->nComp()!=nsolve)
   amrex::Error("localMF[gp_mf+dir]->nComp() invalid29");
  if (localMF[AREA_MF+dir]->boxArray()!=
      localMF[gp_mf+dir]->boxArray())
   amrex::Error("gp_mf boxarrays do not match");
 } // dir=0..sdim-1

  // viscosity 
 if (project_option==3) {

  if (nsolve!=AMREX_SPACEDIM)
   amrex::Error("nsolve invalid30");

  int operation_flag=8;

  if (simple_AMR_BC_flag_viscosity==0) {

    // AMRSYNC_PRES_MF allocated in NavierStokes::allocate_pressure_work_vars
    // AMRSYNC_PRES_MF deleted in NavierStokes::remove_pressure_work_vars()
   if (level<finest_level) {
    avgDown_and_Copy_localMF( // avgdown in tan dir, copy in normal dir.
     pboth_mf,
     pboth_mf,
     AMRSYNC_PRES_MF,
     operation_flag);
   } else if (level==finest_level) {
    // do nothing
   } else
    amrex::Error("level invalid18");

   if ((level>=1)&&(level<=finest_level)) {
    interp_and_Copy_localMF(
     pboth_mf,
     pboth_mf,
     AMRSYNC_PRES_MF,
     operation_flag);
   } else if (level==0) {
    // do nothing
   } else
    amrex::Error("level invalid19");

  } else if (simple_AMR_BC_flag_viscosity==1) {
   // do nothing
  } else
   amrex::Error("simple_AMR_BC_flag_viscosity invalid");

  int im_tensor=-1;
  int idx_elastic_flux=-1;
  init_gradu_tensor(
    im_tensor,
    homflag,
    pboth_mf,
    LOCAL_CELLTENSOR_MF,
    LOCAL_FACETENSOR_MF,
    idx_elastic_flux,
    simple_AMR_BC_flag_viscosity);

  show_norm2(localMF[pboth_mf],0,localMF[pboth_mf]->nComp(),20);
  show_norm2(localMF[LOCAL_CELLTENSOR_MF],0,
     localMF[LOCAL_CELLTENSOR_MF]->nComp(),21);
  show_norm2(localMF[LOCAL_FACETENSOR_MF],0,
     localMF[LOCAL_FACETENSOR_MF]->nComp(),21);

  int nden=nmat*num_state_material;

  int ntensor=AMREX_SPACEDIM*AMREX_SPACEDIM; 

  allocate_flux_register(operation_flag);
  if (localMF[SEM_FLUXREG_MF]->nComp()!=ntensor)
   amrex::Error("localMF[SEM_FLUXREG_MF]->nComp() invalid5");

  resize_levelsetLO(2,LEVELPC_MF);
  debug_ngrow(LEVELPC_MF,2,110);

   // 1. spectral_loop==0 tileloop==0  low order grad U+grad U^T
   // 2. spectral_loop==0 tileloop==1  high order grad U+grad U^T
   // 3. spectral_loop==0 tileloop==2  fluxes+=divu  fluxes*=(-dt)*mu
   // 4. spectral_loop==0 tileloop==3  assign high order fluxes to SEM_FLUX
   // 5. spectral_loop==1 tileloop==3  avg high order fluxes 
   //
   // spectral_loop==1
   //   tileloop=0,1,2 => do nothing
  for (int spectral_loop=0;spectral_loop<end_spectral_loop();spectral_loop++) {
  for (int dir=1;dir<=AMREX_SPACEDIM;dir++) {
  for (int tileloop=0;tileloop<=3;tileloop++) {
 
   if (thread_class::nthreads<1)
    amrex::Error("thread_class::nthreads invalid");
   thread_class::init_d_numPts(localMF[pboth_mf]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
   for (MFIter mfi(*localMF[pboth_mf],use_tiling); mfi.isValid(); ++mfi) {
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

    FArrayBox& velfab=(*localMF[pboth_mf])[mfi];
    FArrayBox& levelpcfab=(*localMF[LEVELPC_MF])[mfi];
    FArrayBox& maskSEMfab=(*localMF[MASKSEM_MF])[mfi];

    FArrayBox& xflux=(*localMF[gp_mf+dir-1])[mfi];

    FArrayBox& xface=(*localMF[FACE_VAR_MF+dir-1])[mfi];

    FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi];  
 
    FArrayBox& tensor_data=(*localMF[LOCAL_FACETENSOR_MF])[mfi];
    FArrayBox& cell_tensor_data=(*localMF[LOCAL_CELLTENSOR_MF])[mfi];
    FArrayBox& mask_tensor_data=(*localMF[MASKSOLIDTENSOR_MF])[mfi];
    FArrayBox& faceLS=(*localMF[LSTENSOR_MF])[mfi];

    FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi];
    // maskcoef=tag if not covered by level+1 or outside the domain.
    FArrayBox& maskcoef_fab=(*localMF[MASKCOEF_MF])[mfi];

    FArrayBox& semfluxfab=(*localMF[SEM_FLUXREG_MF])[mfi];
    int ncfluxreg=semfluxfab.nComp();
    int local_enable_spectral=enable_spectral;

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // in: GODUNOV_3D.F90
    // -dt * visc_coef * viscface * (grad U + grad U^T)
    fort_crossterm(
     &nsolve,
     &tileloop,
     &dir,
     &operation_flag, // 8
     &local_enable_spectral,
     &spectral_loop,
     &ncfluxreg,
     semfluxfab.dataPtr(),
     ARLIM(semfluxfab.loVect()),ARLIM(semfluxfab.hiVect()),
     maskfab.dataPtr(), // 1=fine/fine  0=coarse/fine
     ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
     maskcoef_fab.dataPtr(), // maskcoef=tag if not cov by level+1 or outside.
     ARLIM(maskcoef_fab.loVect()),ARLIM(maskcoef_fab.hiVect()),
     faceLS.dataPtr(),
     ARLIM(faceLS.loVect()),ARLIM(faceLS.hiVect()),
     mask_tensor_data.dataPtr(),
     ARLIM(mask_tensor_data.loVect()),ARLIM(mask_tensor_data.hiVect()),
     tensor_data.dataPtr(),
     ARLIM(tensor_data.loVect()),ARLIM(tensor_data.hiVect()),
     cell_tensor_data.dataPtr(),
     ARLIM(cell_tensor_data.loVect()),ARLIM(cell_tensor_data.hiVect()),
     maskSEMfab.dataPtr(),
     ARLIM(maskSEMfab.loVect()),ARLIM(maskSEMfab.hiVect()),
     xlo,dx,
     &dt_slab,
     &cur_time_slab,
     velfab.dataPtr(),ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
     levelpcfab.dataPtr(),
     ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),
     xflux.dataPtr(),ARLIM(xflux.loVect()),ARLIM(xflux.hiVect()),
     xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()),
     reconfab.dataPtr(),ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
     &facevisc_index,
     &vofface_index,
     &massface_index,
     &ncphys,
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,
     &level,
     &rzflag,
     velbc.dataPtr(),
     &visc_coef,
     &nmat,
     &nden,
     &ntensor,
     &constant_viscosity,
     &homflag);
   } // mfi
} // omp
   ns_reconcile_d_num(142);

   int debug_id=30+spectral_loop*4+tileloop+10*dir;
   show_norm2(localMF[gp_mf+dir-1],0,
     localMF[gp_mf+dir-1]->nComp(),debug_id);

  } // tileloop
  } // dir

  synchronize_flux_register(operation_flag,spectral_loop);
  } // spectral_loop

  if (nsolve!=AMREX_SPACEDIM)
   amrex::Error("nsolve invalid31");

  delete_localMF(LOCAL_FACETENSOR_MF,1);
  delete_localMF(LOCAL_CELLTENSOR_MF,1);
  delete_localMF(MASKSOLIDTENSOR_MF,1);
  delete_localMF(LSTENSOR_MF,1);

  if (homflag==0) {
   // inhomogeneous Neumann BC for viscosity force
   viscous_boundary_fluxes(
    project_option,
    localMF[gp_mf],
    localMF[gp_mf+1],
    localMF[gp_mf+AMREX_SPACEDIM-1],
    nsolve);
  } else if (homflag==1) {
   // do nothing
  } else
   amrex::Error("homflag invalid in apply_pressure_grad");

  if (check_nan==1) {
   if (ParallelDescriptor::IOProcessor()) {
    std::cout << "in apply_pressure_grad \n";
    std::cout << "homflag= " << homflag << '\n';
    std::cout << "energyflag= " << energyflag << '\n';
   }
   std::fflush(NULL);
   check_for_NAN(localMF[pboth_mf],200);
   for (int dir2=0;dir2<AMREX_SPACEDIM;dir2++) {
    check_for_NAN(localMF[gp_mf+dir2],201);
   }
  } else if (check_nan==0) {
   // do nothing
  } else {
   amrex::Error("check_nan invalid");
  }

 } else if (project_option_is_valid(project_option)==1) {

  int num_colors=0;
  Vector<Real> blob_array;
  blob_array.resize(1);
  int blob_array_size=blob_array.size();

  int operation_flag=0;

  if (nsolve!=1)
   amrex::Error("nsolve invalid32");

  Vector<int> scomp;
  Vector<int> ncomp;
  int ncomp_check;
  int state_index;
   //num_materials_combine=1
  get_mm_scomp_solver(
   1,
   project_option,
   state_index,
   scomp,ncomp,ncomp_check);

  Vector<int> dombcpres(2*AMREX_SPACEDIM);
  const BCRec& descbc = get_desc_lst()[state_index].getBC(scomp[0]);
  const int* b_rec=descbc.vect();
  for (int m=0;m<2*AMREX_SPACEDIM;m++)
   dombcpres[m]=b_rec[m];

  resize_levelsetLO(2,LEVELPC_MF);

  int fluxvel_index=0;
  int fluxden_index=AMREX_SPACEDIM;

  allocate_flux_register(operation_flag);
  if (localMF[SEM_FLUXREG_MF]->nComp()!=AMREX_SPACEDIM)
   amrex::Error("localMF[SEM_FLUXREG_MF]->nComp() invalid6");

  if (simple_AMR_BC_flag==0) {

   if (level<finest_level) {
    avgDown_and_Copy_localMF(
     pboth_mf,
     pboth_mf,
     AMRSYNC_PRES_MF,
     operation_flag);
   } else if (level==finest_level) {
    // do nothing
   } else
    amrex::Error("level invalid18");

   if ((level>=1)&&(level<=finest_level)) {
    interp_and_Copy_localMF(
     pboth_mf,
     pboth_mf,
     AMRSYNC_PRES_MF,
     operation_flag);
   } else if (level==0) {
    // do nothing
   } else
    amrex::Error("level invalid19");

  } else if (simple_AMR_BC_flag==1) {
   // do nothing
  } else
   amrex::Error("simple_AMR_BC_flag invalid");

   //spectral_loop==0 (find gradients only from element data and
   //  immediate neighbors.
   //spectral_loop==1 (resolve flux differences at element boundaries)
  for (int spectral_loop=0;spectral_loop<end_spectral_loop();spectral_loop++) {
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
	  //tileloop==0 low  tileloop==1 SEM
  for (int tileloop=0;tileloop<=1;tileloop++) {

   if (thread_class::nthreads<1)
    amrex::Error("thread_class::nthreads invalid");
   thread_class::init_d_numPts(localMF[LEVELPC_MF]->boxArray().d_numPts());
 
#ifdef _OPENMP
#pragma omp parallel
#endif
{
   for (MFIter mfi(*localMF[LEVELPC_MF],use_tiling); mfi.isValid(); ++mfi) {
    BL_ASSERT(grids[mfi.index()] == mfi.validbox());
    const int gridno = mfi.index();
    const Box& tilegrid = mfi.tilebox();
    const Box& fabgrid = grids[gridno];
    const int* tilelo=tilegrid.loVect();
    const int* tilehi=tilegrid.hiVect();
    const int* fablo=fabgrid.loVect();
    const int* fabhi=fabgrid.hiVect();

    const Real* xlo = grid_loc[gridno].lo();

    FArrayBox& xp=(*localMF[AMRSYNC_PRES_MF+dir])[mfi];
    FArrayBox& xgp=(*localMF[gp_mf+dir])[mfi];
    FArrayBox& xcut=(*localMF[FACE_WEIGHT_MF+dir])[mfi]; // A/rho
    FArrayBox& xface=(*localMF[FACE_VAR_MF+dir])[mfi];

    FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi];  

    FArrayBox& maskSEMfab=(*localMF[MASKSEM_MF])[mfi];
    FArrayBox& presfab=(*localMF[pboth_mf])[mfi]; // in: apply_pressure_grad

    FArrayBox& solfab=(*localMF[FSI_GHOST_MAC_MF+dir])[mfi];
    FArrayBox& levelpcfab=(*localMF[LEVELPC_MF])[mfi];

    Vector<int> presbc;
    getBCArray_list(presbc,state_index,gridno,scomp,ncomp);
    if (presbc.size()!=nsolve*AMREX_SPACEDIM*2)
     amrex::Error("presbc.size() invalid");
    Vector<int> velbc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);
  
    Real beta=0.0;

    // mask=tag if not covered by level+1 or outside the domain.
    FArrayBox& maskcoeffab=(*localMF[MASKCOEF_MF])[mfi];
    FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi];
    FArrayBox& semfluxfab=(*localMF[SEM_FLUXREG_MF])[mfi];
    int ncfluxreg=semfluxfab.nComp();
    if (ncfluxreg!=AMREX_SPACEDIM)
     amrex::Error("ncfluxreg invalid");

    int local_enable_spectral=enable_spectral;
    int ncomp_xp=nsolve;
    int ncomp_xgp=nsolve;
    int ncomp_mgoni=presfab.nComp();

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // -grad p * FACE_WEIGHT * dt
    // in: apply_pressure_grad
    fort_cell_to_mac(
     &ncomp_mgoni,
     &ncomp_xp,
     &ncomp_xgp,
     &simple_AMR_BC_flag,
     &nsolve,
     &tileloop,
     &dir,
     &operation_flag, // 0
     &energyflag,
     &beta,
     &visc_coef,
     &interp_vel_increment_from_cell,
     filter_velocity.dataPtr(),
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
     &ignore_div_up,
     &pforce_index,
     &faceden_index, 
     &icemask_index,
     &massface_index,
     &vofface_index,
     &ncphys,
     override_density.dataPtr(),
     constant_density_all_time.dataPtr(),
     presbc.dataPtr(),
     velbc.dataPtr(),
     &slab_step,
     &dt_slab,
     &cur_time_slab,
     xlo,dx,
     &spectral_loop,
     &ncfluxreg,
     semfluxfab.dataPtr(),
     ARLIM(semfluxfab.loVect()),ARLIM(semfluxfab.hiVect()),
     maskfab.dataPtr(), // 1=fine/fine  0=coarse/fine
     ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
     maskcoeffab.dataPtr(), // 1=not covered or outside dom. 0=covered
     ARLIM(maskcoeffab.loVect()),ARLIM(maskcoeffab.hiVect()),
     maskSEMfab.dataPtr(),
     ARLIM(maskSEMfab.loVect()),ARLIM(maskSEMfab.hiVect()),
     levelpcfab.dataPtr(),
     ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),
     solfab.dataPtr(),
     ARLIM(solfab.loVect()),ARLIM(solfab.hiVect()),
     xcut.dataPtr(),ARLIM(xcut.loVect()),ARLIM(xcut.hiVect()),
     xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()),
     reconfab.dataPtr(),ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
     xgp.dataPtr(),ARLIM(xgp.loVect()),ARLIM(xgp.hiVect()),
     xp.dataPtr(),ARLIM(xp.loVect()),ARLIM(xp.hiVect()), //holds AMRSYNC_PRES
     xgp.dataPtr(),ARLIM(xgp.loVect()),ARLIM(xgp.hiVect()), // xvel
     presfab.dataPtr(), //vel
     ARLIM(presfab.loVect()),ARLIM(presfab.hiVect()),
     presfab.dataPtr(),
     ARLIM(presfab.loVect()),ARLIM(presfab.hiVect()),
     presfab.dataPtr(), //den
     ARLIM(presfab.loVect()),ARLIM(presfab.hiVect()),
     presfab.dataPtr(), //mgoni
     ARLIM(presfab.loVect()),ARLIM(presfab.hiVect()),
     presfab.dataPtr(), //color
     ARLIM(presfab.loVect()),ARLIM(presfab.hiVect()),
     presfab.dataPtr(), //type
     ARLIM(presfab.loVect()),ARLIM(presfab.hiVect()),
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
     &project_option,
     &SEM_upwind,
     &SEM_advection_algorithm);

   }  // mfi
} // omp

   ns_reconcile_d_num(143);

  } // tileloop
  } // dir

  synchronize_flux_register(operation_flag,spectral_loop);
  } // spectral_loop

  if (project_option==2) { // thermal conduction
   if (homflag==0) {
    // inhomogeneous Neumann BC for thermal conduction
    viscous_boundary_fluxes(
     project_option,
     localMF[gp_mf],
     localMF[gp_mf+1],
     localMF[gp_mf+AMREX_SPACEDIM-1],
     nsolve);
   } else if (homflag==1) {
    // do nothing
   } else
    amrex::Error("homflag invalid in apply_pressure_grad 2");
  }

 } else
  amrex::Error("project_option invalid27");

} // subroutine apply_pressure_grad


void NavierStokes::make_physics_varsALL(int project_option,
  int post_restart_flag,int caller_id) {

 if (level!=0)
  amrex::Error("level invalid make_physics_varsALL");

 if ((project_option==0)||
     (project_option==1)) {
  // do nothing
 } else
  amrex::Error("project_option invalid make_physics_varsALL");

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid");

 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;
 int nhistory=nten*2;

 int finest_level = parent->finestLevel();

 debug_ngrow(SLOPE_RECON_MF,0,90);
 if (localMF[SLOPE_RECON_MF]->nComp()==nmat*ngeom_recon) {
  // do nothing
 } else
  amrex::Error("localMF[SLOPE_RECON_MF]->nComp() invalid");

 if (1==0) {
  int do_face_decomp=0;
  for (int ilev=level;ilev<=finest_level;ilev++) {
   NavierStokes& ns_level=getLevel(ilev);
   int tessellate=3;
   ns_level.makeFaceFrac(tessellate,1,FACEFRAC_MM_MF,do_face_decomp);
   ns_level.ProcessFaceFrac(tessellate,FACEFRAC_MM_MF,FACEFRAC_SOLVE_MM_MF,0);
   ns_level.makeCellFrac(tessellate,0,CELLFRAC_MM_MF);
  } // ilev
 } else if (1==1) {
  // do nothing
 } else
  amrex::Error("bust");

  // in: NavierStokes::make_physics_varsALL
  // piecewise constant interpolation.
 allocate_levelsetLO_ALL(1,LEVELPC_MF);

   // create DIST_CURV_MF 

 curv_min.resize(thread_class::nthreads);
 curv_max.resize(thread_class::nthreads);

 for (int tid=0;tid<thread_class::nthreads;tid++) {
  curv_min[tid]=1.0e+99;
  curv_max[tid]=-1.0e+99;
 } // tid

 allocate_array(1,nhistory,-1,HISTORY_MF);

 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.makeStateCurv(project_option,post_restart_flag);
 }
  // filenames: "ANGLE_UTAN<stuff>.plt"  (cell centered)
  // if GNBC is used, then the "ghost normal" in the substrate
  // is extrapolated from the interior.
  // if GNBC is not used, then the "ghost normal" is either
  // (a) corresponding to the static angle condition or
  // (b) some dynamic angle condition depending on utan.
  // "angle" = static angle if use_DCA=-1, dynamic angle otherwise.
 if (1==0) {
  writeSanityCheckData(
    "ANGLE_UTAN",
    "in: make_physics_varsALL, after makeStateCurv: HISTORY_MF, angle, utan", 
    caller_id,
    localMF[HISTORY_MF]->nComp(),
    HISTORY_MF,
    -1, // State_Type==-1
    -1); // data_dir==-1
 }
 delete_array(HISTORY_MF);

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.avgDownCURV_localMF(DIST_CURV_MF);
 }

 int save_enable_spectral=enable_spectral;
 override_enable_spectral(viscous_enable_spectral);

 // allocate and delete HOLD_VELOCITY_DATA_MF in init_gradu_tensorALL:
 int simple_AMR_BC_flag_viscosity=1;
 int do_alloc=1; 
 int im_tensor=-1;
 int idx_elastic_flux=-1;
 init_gradu_tensorALL(
   im_tensor,
   HOLD_VELOCITY_DATA_MF,
   do_alloc,
   CELLTENSOR_MF,FACETENSOR_MF,
   idx_elastic_flux,
   simple_AMR_BC_flag_viscosity);

 override_enable_spectral(save_enable_spectral);

  //ngrow=1
 getStateVISC_ALL(CELL_VISC_MATERIAL_MF,1);

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.getStateDen_localMF(DEN_RECON_MF,1,cur_time_slab);
  ns_level.getStateMOM_DEN(MOM_DEN_MF,1,cur_time_slab);
 }

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);

  ns_level.make_physics_vars(project_option);
  ns_level.level_init_icemask();

  int bfact=parent->Space_blockingFactor(ilev);

   // average down from ilev+1 to ilev.
  
    // idxMF,scomp,ncomp,start_dir,ndir
  ns_level.avgDownEdge_localMF(FACE_VAR_MF,facecut_index,1,0,
		  AMREX_SPACEDIM,0,6);
  ns_level.avgDownEdge_localMF(FACE_VAR_MF,icefacecut_index,1,0,
		  AMREX_SPACEDIM,0,7);

  int spectral_override=0;
  if ((projection_enable_spectral==1)||
      (projection_enable_spectral==2)) {
   if (bfact>=2) {
    spectral_override=1;
   } else if (bfact==1) {
    spectral_override=0;
   } else {
    amrex::Error("bfact invalid");
   }
  } else if ((projection_enable_spectral==0)||
             (projection_enable_spectral==3)) {
   spectral_override=0;
  } else
   amrex::Error("projection_enable_spectral invalid");

  ns_level.avgDownEdge_localMF(FACE_VAR_MF,faceden_index,1,0,AMREX_SPACEDIM,
   spectral_override,8);

  ns_level.avgDownEdge_localMF(FACE_VAR_MF,facevisc_index,1,0,AMREX_SPACEDIM,
   0,9);
  ns_level.avgDownEdge_localMF(FACE_VAR_MF,faceheat_index,1,0,AMREX_SPACEDIM,
   0,10);
  if (num_species_var>0) {
   ns_level.avgDownEdge_localMF(FACE_VAR_MF,facespecies_index,
      num_species_var,0,AMREX_SPACEDIM,0,11);
  }

   // spectral_override==0
  ns_level.avgDownEdge_localMF(FACE_VAR_MF,smoothing_index,1,0,AMREX_SPACEDIM,
   0,12);

 }  // ilev=finest_level ... level

 if (1==0) {
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
    "in: make_physics_varsALL, FACE_VAR_MF",//faceden_index=2 facevisc_index=6
    caller_id,
    localMF[FACE_VAR_MF+dir]->nComp(),
    FACE_VAR_MF+dir,
    -1, // State_Type==-1
    dir);
  } // dir=0..sdim-1
 }

 delete_array(DEN_RECON_MF);
 delete_array(MOM_DEN_MF);
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
  delete_array(AMRSYNC_VEL_MF+dir);
 delete_array(CELLTENSOR_MF);
 delete_array(FACETENSOR_MF);

} // subroutine make_physics_varsALL

// called from: prelim_alloc() and make_physics_vars
void NavierStokes::allocate_physics_vars() {

 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid allocate_physics_vars");

 int nmat=num_materials;

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  new_localMF_if_not_exist(FACE_VAR_MF+dir,ncphys,0,dir);
 }

  // ncomp,ngrow,dir
 if (localMF_grow[SWEPT_CROSSING_MF]==-1) {
  new_localMF(SWEPT_CROSSING_MF,nmat,0,-1); 
  setVal_localMF(SWEPT_CROSSING_MF,1.0,0,nmat,0); //val,scomp,ncomp,ngrow
 } else if (localMF_grow[SWEPT_CROSSING_MF]>=0) {
  // do nothing
 } else
  amrex::Error("localMF_grow[SWEPT_CROSSING_MF] invalid");

 new_localMF_if_not_exist(CELL_DEDT_MF,nmat+1,1,-1); // ncomp,ngrow,dir
 new_localMF_if_not_exist(CELL_DEN_MF,nmat+1,1,-1); // ncomp,ngrow,dir
  // coeff_avg,padvect_avg 
 new_localMF_if_not_exist(CELL_SOUND_MF,2,0,-1); // ncomp,ngrow,dir

  // tessellating volume fractions.
 new_localMF_if_not_exist(CELL_VOF_MF,nmat,1,-1); // ncomp,ngrow,dir
 new_localMF_if_not_exist(CELL_VISC_MF,nmat+1,1,-1); // ncomp,ngrow,dir

} // allocate_physics_vars

void NavierStokes::allocate_levelsetLO_ALL(int ngrow,int idx) {

 int finest_level=parent->finestLevel();
 if (level!=0)
  amrex::Error("level invalid allocate_levelsetALL");

 if ((ngrow<0)||(ngrow>ngrow_distance))
  amrex::Error("ngrow invalid");

 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.allocate_levelsetLO(ngrow,idx);
 }

}  // allocate_levelsetLO_ALL

void NavierStokes::allocate_levelsetLO(int ngrow,int idx) {

 int nmat=num_materials;

 if ((ngrow<0)||(ngrow>ngrow_distance))
  amrex::Error("ngrow invalid");

 int Interp_LO=1;
 override_LS_HO(Interp_LO);  // do not use normals for coarse/fine interp.

 delete_localMF_if_exist(idx,1);
 getStateDist_localMF(idx,ngrow,cur_time_slab,17);
 if (localMF[idx]->nComp()!=nmat*(AMREX_SPACEDIM+1))
  amrex::Error("localMF[idx]->nComp()!=nmat*(AMREX_SPACEDIM+1)");
 debug_ngrow(idx,ngrow,90);

 Interp_LO=0;
 override_LS_HO(Interp_LO);  // use normals

} // subroutine allocate_levelsetLO


void NavierStokes::resize_levelsetLO(int ngrow,int idx) {

 if ((ngrow<0)||(ngrow>ngrow_distance))
  amrex::Error("ngrow invalid");

 if (localMF_grow[idx]>=0) {
  // do nothing
 } else {
  std::cout << "in resize_levelsetLO  ngrow= " << ngrow << 
   " idx= " << idx << '\n';
  std::cout << "localMF_grow[idx] = " << localMF_grow[idx] << '\n';
  amrex::Error("localMF_grow[idx]<0");
 }

 if (localMF[idx]->nGrow()==ngrow) {
  // do nothing
 } else {
  allocate_levelsetLO(ngrow,idx);
 }

} // subroutine resize_levelsetLO

// called from make_physics_varsALL
// density vars get 1/rho
void NavierStokes::make_physics_vars(int project_option) {
 
 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();
 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

  // height function curvature
  // finite difference curvature
  // pforce
  // marangoni force (sdim)
  // dir * side (dir=1..sdim, side=-1 or 1)
  // im3
  // x nten
 int num_curv=nten*(AMREX_SPACEDIM+5); 

 if ((project_option==0)||
     (project_option==1)) {
  // do nothing
 } else
  amrex::Error("project_option invalid make_physics_vars");

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
 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
  debug_ngrow(FSI_GHOST_MAC_MF+data_dir,0,112);
 }

 if (localMF[DIST_CURV_MF]->nComp()!=num_curv)
  amrex::Error("localMF[DIST_CURV_MF]->nComp() invalid");

 VOF_Recon_resize(1,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,1,240);

 if (localMF[SLOPE_RECON_MF]->nComp()!=nmat*ngeom_recon)
  amrex::Error("slope_recon_mf has incorrect ncomp");

 debug_ngrow(DIST_CURV_MF,1,241);

 resize_metrics(1);
 debug_ngrow(VOLUME_MF,1,243);
 resize_maskfiner(1,MASKCOEF_MF);
 resize_mask_nbr(1);

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int scomp_pres=AMREX_SPACEDIM;

  // in: make_physics_vars
 allocate_physics_vars();

 MultiFab& tempmf=get_new_data(State_Type,slab_step+1);

 const Real* dx = geom.CellSize();

 debug_ngrow(DEN_RECON_MF,1,6001);
 if (localMF[DEN_RECON_MF]->nComp()!=nmat*num_state_material)
  amrex::Error("den_recon has invalid ncomp in make_physics_vars");
 debug_ngrow(MOM_DEN_MF,1,6001);
 if (localMF[MOM_DEN_MF]->nComp()!=nmat)
  amrex::Error("mom_den has invalid ncomp in make_physics_vars");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  new_localMF(AMRSYNC_VEL_MF+dir,1,0,dir);
  setVal_localMF(AMRSYNC_VEL_MF+dir,1.0e+40,0,1,0);
 }

 Real problo[AMREX_SPACEDIM];
 Real probhi[AMREX_SPACEDIM];
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  problo[dir]=geom.ProbLo(dir);
  probhi[dir]=geom.ProbHi(dir);
 }
  
 localMF[CELL_SOUND_MF]->setVal(0.0,0,2,0);

 MultiFab* vofC=new MultiFab(grids,dmap,nmat,1,
  MFInfo().SetTag("vofC"),FArrayBoxFactory());

 for (int im=0;im<nmat;im++) {
  int scomp=im*ngeom_recon;
  MultiFab::Copy(*vofC,*localMF[SLOPE_RECON_MF],scomp,im,1,1);
 }

  // (dir-1)*2*nmat + (side-1)*nmat + im
 int nrefine_vof=2*nmat*AMREX_SPACEDIM;
 int ngrow_refine=1;
 MultiFab* vofF=new MultiFab(grids,dmap,nrefine_vof,ngrow_refine,
	MFInfo().SetTag("vofF"),FArrayBoxFactory());
 MultiFab* massF=new MultiFab(grids,dmap,nrefine_vof,ngrow_refine,
	MFInfo().SetTag("massF"),FArrayBoxFactory());

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(tempmf.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(tempmf,use_tiling); mfi.isValid(); ++mfi) {
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

  FArrayBox& slopefab=(*localMF[SLOPE_RECON_MF])[mfi];
  FArrayBox& denstatefab=(*localMF[DEN_RECON_MF])[mfi];
  FArrayBox& mom_denfab=(*localMF[MOM_DEN_MF])[mfi];

  FArrayBox& vofFfab=(*vofF)[mfi];
  FArrayBox& massFfab=(*massF)[mfi];

  int tessellate=3;

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  int use_mom_den=1;

    // in: LEVELSET_3D.F90
    // centroid in absolute coordinates.
  fort_build_semirefinevof(
   &tid_current,
   &tessellate,  // =3
   &ngrow_refine,
   &nrefine_vof,
   &nten,
   spec_material_id_AMBIENT.dataPtr(),
   mass_fraction_id.dataPtr(),
   cavitation_vapor_density.dataPtr(),
   override_density.dataPtr(),
   constant_density_all_time.dataPtr(),
   &use_mom_den,
   xlo,dx,
   slopefab.dataPtr(),
   ARLIM(slopefab.loVect()),ARLIM(slopefab.hiVect()),
   denstatefab.dataPtr(),
   ARLIM(denstatefab.loVect()),ARLIM(denstatefab.hiVect()),
   mom_denfab.dataPtr(),
   ARLIM(mom_denfab.loVect()),ARLIM(mom_denfab.hiVect()),
   vofFfab.dataPtr(),ARLIM(vofFfab.loVect()),ARLIM(vofFfab.hiVect()),
   massFfab.dataPtr(),ARLIM(massFfab.loVect()),ARLIM(massFfab.hiVect()),
   tilelo,tilehi,
   fablo,fabhi,
   &bfact,
   &nmat,
   &level,&finest_level);
 }  // mfi
} // omp
 ns_reconcile_d_num(144);

 resize_levelsetLO(2,LEVELPC_MF);

 int ngrow_visc=1;
 MultiFab* modvisc=new MultiFab(grids,dmap,nmat,ngrow_visc,
	MFInfo().SetTag("modvisc"),FArrayBoxFactory());

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(tempmf.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(tempmf,use_tiling); mfi.isValid(); ++mfi) {
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

  FArrayBox& slopefab=(*localMF[SLOPE_RECON_MF])[mfi];
  FArrayBox& denstatefab=(*localMF[DEN_RECON_MF])[mfi];
  FArrayBox& viscstatefab=(*localMF[CELL_VISC_MATERIAL_MF])[mfi];

  FArrayBox& levelpcfab=(*localMF[LEVELPC_MF])[mfi];

  FArrayBox& modviscfab=(*modvisc)[mfi];

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // in: LEVELSET_3D.F90
   // visc_coef passed as a parameter so that Guibo can
   // calculate the dynamic contact angle condition. 
  FORT_BUILD_MODVISC(
   &ngrow_visc,
   &cur_time_slab,
   problo,probhi,
   &visc_coef,
   &nten,
   xlo,dx,
   slopefab.dataPtr(),
   ARLIM(slopefab.loVect()),ARLIM(slopefab.hiVect()),
   denstatefab.dataPtr(), // denstate unused for now.
   ARLIM(denstatefab.loVect()),ARLIM(denstatefab.hiVect()),
   viscstatefab.dataPtr(),
   ARLIM(viscstatefab.loVect()),ARLIM(viscstatefab.hiVect()),
   levelpcfab.dataPtr(),
   ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),
   modviscfab.dataPtr(),ARLIM(modviscfab.loVect()),ARLIM(modviscfab.hiVect()),
   tilelo,tilehi,
   fablo,fabhi,
   &bfact,
   &nmat,
   &level,&finest_level);
 }  // mfi
} // omp
 ns_reconcile_d_num(145);

 Vector< Real > local_curv_min;
 Vector< Real > local_curv_max;
 local_curv_min.resize(thread_class::nthreads);
 local_curv_max.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  local_curv_min[tid]=1.0e+99;
  local_curv_max[tid]=-1.0e+99;
 } // tid

 debug_ngrow(MASKCOEF_MF,1,6001);

 debug_ngrow(MASK_NBR_MF,1,90);
 if (localMF[MASK_NBR_MF]->nComp()!=4)
  amrex::Error("invalid ncomp for mask nbr");
 
  // isweep==0  face variables
  // isweep==1  cenden,cenvof,cenDeDT,cenvisc
 for (int isweep=0;isweep<2;isweep++) {

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(tempmf.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(tempmf,use_tiling); mfi.isValid(); ++mfi) {
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

// face_frac=0 if presbc<> interior or exterior dirichlet.
   Vector<int> presbc=getBCArray(State_Type,gridno,scomp_pres,1);
   Vector<int> velbc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);

   // mask=tag if not covered by level+1 or outside the domain.
   FArrayBox& maskcov=(*localMF[MASKCOEF_MF])[mfi];
   // mask_nbr:
   // (1) =1 interior  =1 fine-fine ghost in domain  =0 otherwise
   // (2) =1 interior  =0 otherwise
   // (3) =1 interior+ngrow-1  =0 otherwise
   // (4) =1 interior+ngrow    =0 otherwise
   FArrayBox& masknbr=(*localMF[MASK_NBR_MF])[mfi];

   FArrayBox& slopefab=(*localMF[SLOPE_RECON_MF])[mfi];
   FArrayBox& curvfab=(*localMF[DIST_CURV_MF])[mfi];

   FArrayBox& denstatefab=(*localMF[DEN_RECON_MF])[mfi];
   FArrayBox& mom_denfab=(*localMF[MOM_DEN_MF])[mfi];

   FArrayBox& viscstatefab=(*localMF[CELL_VISC_MATERIAL_MF])[mfi];

   FArrayBox& levelpcfab=(*localMF[LEVELPC_MF])[mfi];

   FArrayBox& solxfab=(*localMF[FSI_GHOST_MAC_MF])[mfi];
   FArrayBox& solyfab=(*localMF[FSI_GHOST_MAC_MF+1])[mfi];
   FArrayBox& solzfab=(*localMF[FSI_GHOST_MAC_MF+AMREX_SPACEDIM-1])[mfi];
 
   FArrayBox& xfacevar=(*localMF[FACE_VAR_MF])[mfi];
   FArrayBox& yfacevar=(*localMF[FACE_VAR_MF+1])[mfi];
   FArrayBox& zfacevar=(*localMF[FACE_VAR_MF+AMREX_SPACEDIM-1])[mfi];

   // stores 1/(rho cv)   (cv=DeDT)
   FArrayBox& cDeDTfab=(*localMF[CELL_DEDT_MF])[mfi];
   FArrayBox& cdenfab=(*localMF[CELL_DEN_MF])[mfi];  // 1/rho

    // CELL_VOF_MF has the tessellating volume fractions.
   FArrayBox& cvoffab=(*localMF[CELL_VOF_MF])[mfi];  
   FArrayBox& cviscfab=(*localMF[CELL_VISC_MF])[mfi];

   FArrayBox& volfab=(*localMF[VOLUME_MF])[mfi];
   FArrayBox& vofCfab=(*vofC)[mfi];

   FArrayBox& vofFfab=(*vofF)[mfi];
   FArrayBox& massFfab=(*massF)[mfi];
   FArrayBox& modviscfab=(*modvisc)[mfi];

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // in: LEVELSET_3D.F90
   fort_init_physics_vars(
    &tid_current,
    &FD_curv_interp, 
    &local_curv_min[tid_current],
    &local_curv_max[tid_current],
    &isweep,
    &nrefine_vof,
    &smoothing_length_scale,
    denconst_interface.dataPtr(),
    viscconst_interface.dataPtr(),
    heatviscconst_interface.dataPtr(),
    speciesviscconst_interface.dataPtr(),
    &curv_index,
    &pforce_index,
    &faceden_index,
    &facecut_index,
    &icefacecut_index,
    &icemask_index,
    &facevisc_index,
    &faceheat_index,
    &facevel_index,
    &facespecies_index,
    &smoothing_index,
    &massface_index,
    &vofface_index,
    &ncphys,
    latent_heat.dataPtr(),
    freezing_model.dataPtr(),
    distribute_from_target.dataPtr(),
    &solidheat_flag,
    microlayer_size.dataPtr(), 
    microlayer_substrate.dataPtr(), 
    microlayer_temperature_substrate.dataPtr(), 
    spec_material_id_AMBIENT.dataPtr(),
    mass_fraction_id.dataPtr(),
    cavitation_vapor_density.dataPtr(),
    override_density.dataPtr(),
    constant_density_all_time.dataPtr(),
    &cur_time_slab,
    &dt_slab, //calling FORT_INIT_PHYSICS_VARS
    &project_option,
    problo,probhi,
    &visc_coef,
    &nten,
    xlo,dx,
    maskcov.dataPtr(),ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
    masknbr.dataPtr(),ARLIM(masknbr.loVect()),ARLIM(masknbr.hiVect()),
    xfacevar.dataPtr(),ARLIM(xfacevar.loVect()),ARLIM(xfacevar.hiVect()),
    yfacevar.dataPtr(),ARLIM(yfacevar.loVect()),ARLIM(yfacevar.hiVect()),
    zfacevar.dataPtr(),ARLIM(zfacevar.loVect()),ARLIM(zfacevar.hiVect()),
    curvfab.dataPtr(),
    ARLIM(curvfab.loVect()),ARLIM(curvfab.hiVect()),
    slopefab.dataPtr(),
    ARLIM(slopefab.loVect()),ARLIM(slopefab.hiVect()),
    denstatefab.dataPtr(),
    ARLIM(denstatefab.loVect()),ARLIM(denstatefab.hiVect()),
    mom_denfab.dataPtr(),
    ARLIM(mom_denfab.loVect()),ARLIM(mom_denfab.hiVect()),
    viscstatefab.dataPtr(),
    ARLIM(viscstatefab.loVect()),ARLIM(viscstatefab.hiVect()),
    solxfab.dataPtr(),
    ARLIM(solxfab.loVect()),ARLIM(solxfab.hiVect()),
    solyfab.dataPtr(),
    ARLIM(solyfab.loVect()),ARLIM(solyfab.hiVect()),
    solzfab.dataPtr(),
    ARLIM(solzfab.loVect()),ARLIM(solzfab.hiVect()),
    cDeDTfab.dataPtr(),
    ARLIM(cDeDTfab.loVect()),ARLIM(cDeDTfab.hiVect()),
    cdenfab.dataPtr(),ARLIM(cdenfab.loVect()),ARLIM(cdenfab.hiVect()),
    cvoffab.dataPtr(),ARLIM(cvoffab.loVect()),ARLIM(cvoffab.hiVect()),
    cviscfab.dataPtr(),ARLIM(cviscfab.loVect()),ARLIM(cviscfab.hiVect()),
    volfab.dataPtr(),ARLIM(volfab.loVect()),ARLIM(volfab.hiVect()),
    levelpcfab.dataPtr(),
    ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),
    vofCfab.dataPtr(),ARLIM(vofCfab.loVect()),ARLIM(vofCfab.hiVect()),
    vofFfab.dataPtr(),ARLIM(vofFfab.loVect()),ARLIM(vofFfab.hiVect()),
    massFfab.dataPtr(),ARLIM(massFfab.loVect()),ARLIM(massFfab.hiVect()),
    modviscfab.dataPtr(),ARLIM(modviscfab.loVect()),ARLIM(modviscfab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    presbc.dataPtr(), 
    velbc.dataPtr(), 
    &nmat,
    &nparts,
    &nparts_def,
    im_solid_map_ptr,
    &num_curv,
    &level,
    &finest_level);
  }  // mfi
} // omp
  ns_reconcile_d_num(146);
 } // isweep

  // correction and sanity check if spectral element method.
 density_TO_MAC(project_option);

 for (int tid=1;tid<thread_class::nthreads;tid++) {
  if (local_curv_min[tid]<local_curv_min[0])
   local_curv_min[0]=local_curv_min[tid];
  if (local_curv_max[tid]>local_curv_max[0])
   local_curv_max[0]=local_curv_max[tid];
 } // tid

 ParallelDescriptor::ReduceRealMin(local_curv_min[0]);
 ParallelDescriptor::ReduceRealMax(local_curv_max[0]);

 if ((fab_verbose==1)||(fab_verbose==3)) {

   std::cout << "make_physics_vars \n";
   std::cout << "c++ level,finest_level " << level << ' ' <<
     finest_level << '\n';
   std::cout << "c++ ngrow,csten " << ngrow_distance << ' ' <<
     curv_stencil_height << ' ' << '\n';

   std::cout << "local_curv_min= " << local_curv_min[0] << '\n';
   std::cout << "local_curv_max= " << local_curv_max[0] << '\n';

   for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

    if (thread_class::nthreads<1)
     amrex::Error("thread_class::nthreads invalid");
    thread_class::init_d_numPts(tempmf.boxArray().d_numPts());

    for (MFIter mfi(tempmf,false); mfi.isValid(); ++mfi) {
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
     std::cout << "output of face_var_mf dir= " << dir << '\n';
     int interior_only=0;
     FArrayBox& curvfab=(*localMF[FACE_VAR_MF+dir])[mfi];
     tecplot_debug(curvfab,xlo,fablo,fabhi,dx,dir,0,curv_index,
      1,interior_only);
    } // mfi
    ns_reconcile_d_num(147);

   } // dir=0..sdim-1

 } // ((fab_verbose==1)||(fab_verbose==3))

 delete vofC;
 delete vofF;
 delete massF;
 delete modvisc;
 
} // end subroutine make_physics_vars

void NavierStokes::solid_temperature() {
 
 bool use_tiling=ns_tiling;

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int nmat=num_materials;

 MultiFab &S_new = get_new_data(State_Type,slab_step+1);
 int nstate=(AMREX_SPACEDIM+1)+
  nmat*(num_state_material+ngeom_raw)+1;
 if (S_new.nComp()!=nstate) 
  amrex::Error("S_new.nComp()!=nstate");

 MultiFab &LS_new = get_new_data(LS_Type,slab_step+1);
 if (LS_new.nComp()!=nmat*(AMREX_SPACEDIM+1)) 
  amrex::Error("LS_new.nComp()!=nmat*(AMREX_SPACEDIM+1)");

 int solid_exists=0;
 for (int im=0;im<nmat;im++) {
  if (ns_is_rigid(im)==1) {
   solid_exists=1;
  } else if (ns_is_rigid(im)==0) {
   // do nothing
  } else
   amrex::Error("ns_is_rigid invalid");
 } // im

 if (solid_exists==0) {
  // do nothing
 } else if (solid_exists==1) {
  const Real* dx = geom.CellSize();

  if (solidheat_flag==0) { // diffuse in the solid

   // do nothing

  } else if ((solidheat_flag==1)||  // dirichlet at solid/fluid
             (solidheat_flag==2)) { // neumann at solid/fluid

   int dcomp=(AMREX_SPACEDIM+1);
   int nden=nmat*num_state_material;

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

    FArrayBox& snewfab=S_new[mfi];
    FArrayBox& lsnewfab=LS_new[mfi];

    const Real* xlo = grid_loc[gridno].lo();

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     // if (FSI_flag(im)==1,2,4),
     //  T(im)=TSOLID
     // else if (FSI_flag(im)==0,3,5,6,7),
     //  T(im)=TSOLID if in the solid.
     // 
     // note: if FSI_flag==2, then solid temperature copied to itself since
     //  solid temp initialized in another routine.
     // (i.e. no change here)
    fort_initsolidtemp(
     &nmat,
     &nden,
     &cur_time_slab,
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,
     snewfab.dataPtr(dcomp),
     ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()),
     lsnewfab.dataPtr(),
     ARLIM(lsnewfab.loVect()),ARLIM(lsnewfab.hiVect()),
     dx,xlo);  
 
   } // mfi
} // omp
   ns_reconcile_d_num(148);
  } else
   amrex::Error("solidheat_flag invalid");

 } else
  amrex::Error("solid_exists invalid solid_temperature");

} // solid_temperature

// adds gravity force and surface tension 
// to cell and face velocity (all components).
void NavierStokes::increment_potential_forceALL() {

 if (level!=0)
  amrex::Error("level invalid increment_potential_forceALL");

 int finest_level=parent->finestLevel();
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.increment_potential_force();
 }

} // subroutine increment_potential_forceALL()

void NavierStokes::increment_potential_force() {

 int finest_level=parent->finestLevel();
 
 bool use_tiling=ns_tiling;

 int nmat=num_materials;

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 int nstate=S_new.nComp();
 if (nstate!=(AMREX_SPACEDIM+1)+
     nmat*(num_state_material+ngeom_raw)+1)
  amrex::Error("nstate invalid");
 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 MultiFab& LS_new=get_new_data(LS_Type,slab_step+1);
 if (LS_new.nComp()!=nmat*(AMREX_SPACEDIM+1))
  amrex::Error("LS_new.nComp()!=nmat*(AMREX_SPACEDIM+1)");

 debug_ngrow(POTENTIAL_FORCE_CELL_MF,0,260);
 if (localMF[POTENTIAL_FORCE_CELL_MF]->nComp()!=AMREX_SPACEDIM)
  amrex::Error("localMF[POTENTIAL_FORCE_CELL_MF]->nComp() invalid");

 Real gravity_normalized=std::abs(gravity);
 if (invert_gravity==1)
  gravity_normalized=-gravity_normalized;
 else if (invert_gravity==0) {
  // do nothing
 } else
  amrex::Error("invert_gravity invalid");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

  debug_ngrow(POTENTIAL_FORCE_EDGE_MF+dir,0,261);
  if (localMF[POTENTIAL_FORCE_EDGE_MF+dir]->nComp()!=1) {
   std::cout << "ncomp=" << 
    localMF[POTENTIAL_FORCE_EDGE_MF+dir]->nComp() << 
    " dir= " << dir << '\n';
   amrex::Error("localMF[POTENTIAL_FORCE_EDGE_MF+dir]->nComp() invalid");
  }

  VOF_Recon_resize(1,SLOPE_RECON_MF);
  debug_ngrow(SLOPE_RECON_MF,1,240);

  const Real* dx = geom.CellSize();
  MultiFab& Umac_new=get_new_data(Umac_Type+dir,slab_step+1);

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
   FArrayBox& macfab=Umac_new[mfi];
   FArrayBox& facegrav=(*localMF[POTENTIAL_FORCE_EDGE_MF+dir])[mfi];
   FArrayBox& cellgrav=(*localMF[POTENTIAL_FORCE_CELL_MF])[mfi];
   FArrayBox& xfacefab=(*localMF[FACE_VAR_MF+dir])[mfi];
   FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi];
   FArrayBox& lsfab=LS_new[mfi];

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // in: NAVIERSTOKES_3D.F90
    // (gravity and surface tension)
    // u+=facegrav 
    // u+=cellgrav 
   FORT_ADDGRAVITY(
     &dt_slab,
     &cur_time_slab,
     &gravity_potential_form,
     &gravity_normalized,
     &gravity_dir,
     &angular_velocity,
     denconst_gravity.dataPtr(),
     &level,
     &finest_level,
     &facecut_index,
     &icefacecut_index,
     &vofface_index,
     &ncphys,
     &nmat,
     &nstate,
     tilelo,tilehi,
     fablo,fabhi,&bfact,
     xlo,dx,&dir,
     xfacefab.dataPtr(), 
     ARLIM(xfacefab.loVect()),ARLIM(xfacefab.hiVect()),
     reconfab.dataPtr(),
     ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
     lsfab.dataPtr(),ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
     snewfab.dataPtr(),
     ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()),
     macfab.dataPtr(),
     ARLIM(macfab.loVect()),ARLIM(macfab.hiVect()),
     cellgrav.dataPtr(),
     ARLIM(cellgrav.loVect()),ARLIM(cellgrav.hiVect()),
     facegrav.dataPtr(),
     ARLIM(facegrav.loVect()),ARLIM(facegrav.hiVect()) );
  } // mfi
} // omp
  ns_reconcile_d_num(149);
 }  // dir

} // increment_potential_force

// called from multiphase_project when 
// project_option==0 
void NavierStokes::deallocate_potential_forceALL() {

 int finest_level=parent->finestLevel();
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.delete_localMF(POTENTIAL_FORCE_EDGE_MF,AMREX_SPACEDIM);
  ns_level.delete_localMF(POTENTIAL_EDGE_MF,AMREX_SPACEDIM);
  ns_level.delete_localMF(POTENTIAL_FORCE_CELL_MF,1);
  ns_level.delete_localMF(AMRSYNC_PRES_MF,AMREX_SPACEDIM);
 }
} // deallocate_potential_forceALL

// called from multiphase_project when 
// project_option==0 
void NavierStokes::process_potential_forceALL() {

 int finest_level=parent->finestLevel();
 if (level!=0)
  amrex::Error("level!=0");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 if ((SDC_outer_sweeps>=0)&&(SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid process_potential_forceALL");

 allocate_array(1,2,-1,HYDROSTATIC_PRESDEN_MF);

 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   if (ns_level.localMF_grow[AMRSYNC_PRES_MF+dir]==-1) {
     // deallocated in deallocate_potential_forceALL
     // HYDROSTATIC_PRESSURE and HYDROSTATIC_DENSITY
    ns_level.new_localMF(AMRSYNC_PRES_MF+dir,2,0,dir);
    ns_level.setVal_localMF(AMRSYNC_PRES_MF+dir,1.0e+40,0,2,0);
    if (ns_level.localMF[AREA_MF+dir]->boxArray()!=
        ns_level.localMF[AMRSYNC_PRES_MF+dir]->boxArray())
     amrex::Error("AMRSYNC_PRES boxarray does not match");
   } else
    amrex::Error("ns_level.localMF_grow[AMRSYNC_PRES_MF+dir] bad");
  } // dir=0..sdim-1
 } // ilev=level..finest_level
 

 // must go from coarsest to finest level since
 // hydrostatic pressure and density are interpolated (pcinterp)
 // from coarser levels in order to fill a layer of 
 // ghost cells.
 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.init_gravity_potential();
 }

  // must go from coarsest to finest in order
  // to interpolate AMRSYNC_PRES info.
 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.process_potential_force_face();
 }

  // must go from finest to coarsest in order
  // to average down face pressures.
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.process_potential_force_cell();
 }

 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.delete_localMF(HYDROSTATIC_PRESDEN_MF,1);
 }

}  // end subroutine process_potential_forceALL

void NavierStokes::init_gravity_potential() {

 
 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid init_gravity_potential");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int nmat=num_materials;
 int pcomp=AMREX_SPACEDIM;

 MultiFab* dendata=getStateDen(1,cur_time_slab);

 const Real* dx = geom.CellSize();
 const Box& domain = geom.Domain();
 const int* domlo = domain.loVect();
 const int* domhi = domain.hiVect();

 Vector<int> dombcpres(2*AMREX_SPACEDIM);
 const BCRec& descbc = get_desc_lst()[State_Type].getBC(pcomp);
 const int* b_rec=descbc.vect();
 for (int m=0;m<2*AMREX_SPACEDIM;m++)
  dombcpres[m]=b_rec[m];

 Real gravity_normalized=std::abs(gravity);

 if (gravity_potential_form==1) {
  // do nothing
 } else if (gravity_potential_form==0) {
  gravity_normalized=0.0;
 } else
  amrex::Error("gravity_potential_form invalid");

 if (invert_gravity==1)
  gravity_normalized=-gravity_normalized;
 else if (invert_gravity==0) {
  // do nothing
 } else
  amrex::Error("invert_gravity invalid");

 int bfact=parent->Space_blockingFactor(level);
 int ngrow=1;

  // isweep=0 => interior cells updated, coarse_lev.avgDown_localMF,
  //   PCINTERP_fill_borders
  // isweep=1 => exterior cells outside domain are updated:
  //   REFLECT_EVEN BC if wall, EXT_DIR BC on the wall if
  //   outflow.
 for (int isweep=0;isweep<=1;isweep++) {

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(dendata->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(*dendata,use_tiling); mfi.isValid(); ++mfi) {
    BL_ASSERT(grids[mfi.index()] == mfi.validbox());
    const int gridno = mfi.index();
    const Box& tilegrid = mfi.tilebox();
    const Box& fabgrid = grids[gridno];
    const int* tilelo=tilegrid.loVect();
    const int* tilehi=tilegrid.hiVect();
    const int* fablo=fabgrid.loVect();
    const int* fabhi=fabgrid.hiVect();
      // bfact defined above.

    const Real* xlo = grid_loc[gridno].lo();
    FArrayBox& presdenfab=(*localMF[HYDROSTATIC_PRESDEN_MF])[mfi];
    FArrayBox& statefab=(*dendata)[mfi];

    Vector<int> presbc=getBCArray(State_Type,gridno,pcomp,1);

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     // isweep=0 => interior cells updated
     // isweep=1 => exterior cells outside domain are updated:
     //   REFLECT_EVEN BC if wall, EXT_DIR BC on the wall if
     //   outflow.
     // in: NAVIERSTOKES_3D.F90
    FORT_INITPOTENTIAL(
     &nmat,
     &ngrow,
     &override_density[0], 
     presdenfab.dataPtr(),
     ARLIM(presdenfab.loVect()),ARLIM(presdenfab.hiVect()),
     statefab.dataPtr(),
     ARLIM(statefab.loVect()),ARLIM(statefab.hiVect()),
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,
     &level,
     presbc.dataPtr(),
     dombcpres.dataPtr(),
     domlo,domhi,
     xlo,dx,
     &dt_slab,
     &gravity_normalized,
     &gravity_dir,
     &angular_velocity,
     &isweep);
  } // mfi
} // omp
  ns_reconcile_d_num(150);

  if (isweep==0) {

   if (level>0) {
    NavierStokes& coarse_lev = getLevel(level-1);
    int scomp=0;
    int ncomp=2; // HYDROSTATIC_PRESSURE and HYDROSTATIC_DENSITY
    coarse_lev.avgDown_localMF(HYDROSTATIC_PRESDEN_MF,scomp,ncomp,1);
   }

    // set_extrap_bc(bc,phys_bc)
    // FORT_EXTRAPFILL
    // pc_interp or sem_interp
   Vector<int> scompBC_map;
   scompBC_map.resize(2);
   scompBC_map[0]=0;
   scompBC_map[1]=0;

   int extrap_enable_spectral=projection_enable_spectral;
   override_enable_spectralGHOST(0,1,extrap_enable_spectral);
   PCINTERP_fill_borders(HYDROSTATIC_PRESDEN_MF,ngrow,0,2,
     State_Type,scompBC_map);
   
   extrap_enable_spectral=0;
   override_enable_spectralGHOST(0,1,extrap_enable_spectral);

  } else if (isweep==1) {
   // do nothing
  } else
   amrex::Error("isweep invalid");
 
 } // for isweep = 0..1

 delete dendata;

}  // init_gravity_potential

// called from: NavierStokes::process_potential_forceALL()
// u^cell = u^cell + cellgravforce - grad^cell p
// u^face = u^face + facegravforce - grad^face p
// reflecting boundary conditions on ppot should be identical to the
// reflecting boundary conditions on p so that
// if facegrav = grad^face ppot, then this implies that
//  div (grad^face ppot - grad^face p)=0 only if ppot=p.
void NavierStokes::process_potential_force_face() {

 int finest_level=parent->finestLevel();

 int operation_flag=2;
 
 bool use_tiling=ns_tiling;

 if ((SDC_outer_sweeps>=0)&&(SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid process_potential_force_face");

 int num_colors=0;
 Vector<Real> blob_array;
 blob_array.resize(1);
 int blob_array_size=blob_array.size();

 int nmat=num_materials;
 int nsolve=1;

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

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 resize_maskfiner(1,MASKCOEF_MF);
 resize_mask_nbr(1);
 VOF_Recon_resize(1,SLOPE_RECON_MF);
 resize_metrics(1);
 resize_levelsetLO(2,LEVELPC_MF);

 debug_ngrow(MASKCOEF_MF,1,253); // maskcoef=1 if not covered by finer level.
 debug_ngrow(MASK_NBR_MF,1,253); // mask_nbr=1 at fine-fine bc.
 debug_ngrow(SLOPE_RECON_MF,1,130);

 int fluxvel_index=0;
 int fluxden_index=AMREX_SPACEDIM;

 int pcomp=AMREX_SPACEDIM;

 MultiFab &S_new = get_new_data(State_Type,slab_step+1);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  new_localMF(POTENTIAL_FORCE_EDGE_MF+dir,1,0,dir);//grad ppot/rhopot
  new_localMF(POTENTIAL_EDGE_MF+dir,3,0,dir); //ppot,Ften-,Ften+
 }
 new_localMF(POTENTIAL_FORCE_CELL_MF,AMREX_SPACEDIM,0,-1);

 const Real* dx = geom.CellSize();
 const Box& domain = geom.Domain();
 const int* domlo = domain.loVect();
 const int* domhi = domain.hiVect();

 Vector<int> dombcpres(2*AMREX_SPACEDIM);
 const BCRec& descbc = get_desc_lst()[State_Type].getBC(pcomp);
 const int* b_rec=descbc.vect();
 for (int m=0;m<2*AMREX_SPACEDIM;m++)
  dombcpres[m]=b_rec[m];

 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

 MultiFab* dendata=getStateDen(1,cur_time_slab);

  // gpx/rhox,px,gpy/rhoy,py,gpz/rhoz,pz
 allocate_flux_register(operation_flag);
 if (localMF[SEM_FLUXREG_MF]->nComp()!=2*AMREX_SPACEDIM)
  amrex::Error("localMF[SEM_FLUXREG_MF]->nComp() invalid7");

 if (projection_enable_spectral!=enable_spectral)
  amrex::Error("projection_enable_spectral!=enable_spectral");

  // enable_spectral=0,3 => end_spectral_loop()=1
  // enable_spectral=1,2 => end_spectral_loop()=(bfact==1 ? 1:2)
 for (int spectral_loop=0;spectral_loop<end_spectral_loop();spectral_loop++) {
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
 for (int tileloop=0;tileloop<=1;tileloop++) {

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
   int bfact_c=bfact;
   int bfact_f=bfact;
   if (level>0)
    bfact_c=parent->Space_blockingFactor(level-1);
   if (level<finest_level)
    bfact_f=parent->Space_blockingFactor(level+1);

   const Real* xlo = grid_loc[gridno].lo();

   FArrayBox& xgp=(*localMF[POTENTIAL_FORCE_EDGE_MF+dir])[mfi];
   if (xgp.nComp()!=1)
    amrex::Error("xgp.nComp() invalid");
   FArrayBox& xp=(*localMF[POTENTIAL_EDGE_MF+dir])[mfi];
   if (xp.nComp()!=3)
    amrex::Error("xp.nComp() invalid");

   FArrayBox& xface=(*localMF[FACE_VAR_MF+dir])[mfi];
 
   FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi];  

   FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi];

   // mask=tag if not covered by level+1 or outside the domain.
   FArrayBox& maskcoef=(*localMF[MASKCOEF_MF])[mfi];

   FArrayBox& maskSEMfab=(*localMF[MASKSEM_MF])[mfi];
   FArrayBox& presdenfab=(*localMF[HYDROSTATIC_PRESDEN_MF])[mfi];
   FArrayBox& mgonifab=(*dendata)[mfi];

   FArrayBox& solfab=(*localMF[FSI_GHOST_MAC_MF+dir])[mfi];
   FArrayBox& levelpcfab=(*localMF[LEVELPC_MF])[mfi];

   Vector<int> presbc=getBCArray(State_Type,gridno,pcomp,1);
   Vector<int> velbc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);

   Real beta=0.0;

   FArrayBox& semfluxfab=(*localMF[SEM_FLUXREG_MF])[mfi];
   int ncfluxreg=semfluxfab.nComp();
   if (ncfluxreg!=2*AMREX_SPACEDIM) 
    amrex::Error("ncfluxreg invalid");

   int rzflag=0;
   if (geom.IsRZ())
    rzflag=1;
   else if (geom.IsCartesian())
    rzflag=0;
   else if (geom.IsCYLINDRICAL())
    rzflag=3;
   else
    amrex::Error("CoordSys bust 21");

   int local_energyflag=0;
   int local_project_option=0;
   int local_enable_spectral=enable_spectral;

   int simple_AMR_BC_flag=1;

   int ncomp_xp=3;
   int ncomp_xgp=1;
   int ncomp_mgoni=mgonifab.nComp();

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // process_potential_force_face 
   fort_cell_to_mac( 
    &ncomp_mgoni,
    &ncomp_xp,
    &ncomp_xgp,
    &simple_AMR_BC_flag,
    &nsolve,
    &tileloop,
    &dir,
    &operation_flag, //2
    &local_energyflag,
    &beta,
    &visc_coef,
    &interp_vel_increment_from_cell,
    filter_velocity.dataPtr(),
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
    &ignore_div_up,
    &pforce_index,
    &faceden_index,
    &icemask_index,
    &massface_index,
    &vofface_index,
    &ncphys,
    override_density.dataPtr(),
    constant_density_all_time.dataPtr(),
    presbc.dataPtr(),
    velbc.dataPtr(),
    &slab_step,
    &dt_slab,
    &cur_time_slab,
    xlo,dx,
    &spectral_loop,
    &ncfluxreg,
    semfluxfab.dataPtr(),
    ARLIM(semfluxfab.loVect()),ARLIM(semfluxfab.hiVect()),
    maskfab.dataPtr(), // maskfab=1 at fine/fine bc
    ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
    maskcoef.dataPtr(), // maskcoef=1 if not covered or outside domain.
    ARLIM(maskcoef.loVect()),ARLIM(maskcoef.hiVect()),
    maskSEMfab.dataPtr(),
    ARLIM(maskSEMfab.loVect()),ARLIM(maskSEMfab.hiVect()),
    levelpcfab.dataPtr(),
    ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),
    solfab.dataPtr(),
    ARLIM(solfab.loVect()),ARLIM(solfab.hiVect()),
    xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()), //xcut
    xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()), 
    reconfab.dataPtr(),ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
    xgp.dataPtr(),ARLIM(xgp.loVect()),ARLIM(xgp.hiVect()), 
    xp.dataPtr(),ARLIM(xp.loVect()),ARLIM(xp.hiVect()), 
    xgp.dataPtr(),ARLIM(xgp.loVect()),ARLIM(xgp.hiVect()), //xvel
    presdenfab.dataPtr(),
    ARLIM(presdenfab.loVect()),ARLIM(presdenfab.hiVect()), //vel
    presdenfab.dataPtr(0),  // HYDROSTATIC_PRESSURE
    ARLIM(presdenfab.loVect()),ARLIM(presdenfab.hiVect()), 
    presdenfab.dataPtr(1),  // HYDROSTATIC_DENSITY
    ARLIM(presdenfab.loVect()),ARLIM(presdenfab.hiVect()), 
    mgonifab.dataPtr(), // mgoni
    ARLIM(mgonifab.loVect()),ARLIM(mgonifab.hiVect()), 
    mgonifab.dataPtr(), // color
    ARLIM(mgonifab.loVect()),ARLIM(mgonifab.hiVect()),
    mgonifab.dataPtr(), // type
    ARLIM(mgonifab.loVect()),ARLIM(mgonifab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,
    &bfact,&bfact_c,&bfact_f,
    &level,&finest_level,
    &rzflag,
    domlo,domhi,
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
    &local_project_option,
    &SEM_upwind,
    &SEM_advection_algorithm);
  } // mfi
} // omp
  ns_reconcile_d_num(151);
 } // tileloop
 } // dir
 synchronize_flux_register(operation_flag,spectral_loop);
 } // spectral_loop

 delete dendata;

}  // subroutine process_potential_force_face


// called from: NavierStokes::process_potential_forceALL()
// u^cell = u^cell + cellgravforce - grad^cell p
// u^face = u^face + facegravforce - grad^face p
// reflecting boundary conditions on ppot should be identical to the
// reflecting boundary conditions on p so that
// if facegrav = grad^face ppot, then this implies that
//  div (grad^face ppot - grad^face p)=0 only if ppot=p.
void NavierStokes::process_potential_force_cell() {

 int finest_level=parent->finestLevel();

 int operation_flag=4;
 
 bool use_tiling=ns_tiling;

 if ((SDC_outer_sweeps>=0)&&(SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid process_potential_force_cell");

 int nmat=num_materials;
 int nsolve=1;

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

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 debug_ngrow(MASKCOEF_MF,1,253); // maskcoef=1 if not covered by finer level.
 debug_ngrow(MASK_NBR_MF,1,253); // mask_nbr=1 at fine-fine bc.
 debug_ngrow(SLOPE_RECON_MF,1,132);

 debug_ngrow(VOLUME_MF,0,751);
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  if (localMF[AREA_MF+dir]->boxArray()!=
      localMF[FACE_VAR_MF+dir]->boxArray())
   amrex::Error("boxarrays do not match");
 }

 int fluxvel_index=0;
 int fluxden_index=AMREX_SPACEDIM;

 int pcomp=AMREX_SPACEDIM;

 MultiFab &S_new = get_new_data(State_Type,slab_step+1);

 const Real* dx = geom.CellSize();

 Vector<int> dombcpres(2*AMREX_SPACEDIM);
 const BCRec& descbc = get_desc_lst()[State_Type].getBC(pcomp);
 const int* b_rec=descbc.vect();
 for (int m=0;m<2*AMREX_SPACEDIM;m++)
  dombcpres[m]=b_rec[m];

 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

 if (projection_enable_spectral!=enable_spectral)
  amrex::Error("projection_enable_spectral!=enable_spectral");

 int ncomp_edge=-1;
 int scomp=0;
  // if level<finest_level then avgdown from level+1 to level.
 avgDownEdge_localMF(POTENTIAL_FORCE_EDGE_MF,scomp,ncomp_edge,
   0,AMREX_SPACEDIM,1,12);
 avgDownEdge_localMF(POTENTIAL_EDGE_MF,scomp,ncomp_edge,
   0,AMREX_SPACEDIM,1,13);

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

  FArrayBox& xface=(*localMF[FACE_VAR_MF])[mfi];
  FArrayBox& yface=(*localMF[FACE_VAR_MF+1])[mfi];
  FArrayBox& zface=(*localMF[FACE_VAR_MF+AMREX_SPACEDIM-1])[mfi];

  FArrayBox& maskSEMfab=(*localMF[MASKSEM_MF])[mfi];
  FArrayBox& presdenfab=(*localMF[HYDROSTATIC_PRESDEN_MF])[mfi];

  FArrayBox& solxfab=(*localMF[FSI_GHOST_MAC_MF])[mfi];
  FArrayBox& solyfab=(*localMF[FSI_GHOST_MAC_MF+1])[mfi];
  FArrayBox& solzfab=(*localMF[FSI_GHOST_MAC_MF+AMREX_SPACEDIM-1])[mfi];

  FArrayBox& levelpcfab=(*localMF[LEVELPC_MF])[mfi];

  Vector<int> presbc=getBCArray(State_Type,gridno,pcomp,1);
  Vector<int> velbc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);

  FArrayBox& volfab=(*localMF[VOLUME_MF])[mfi];
  FArrayBox& areax=(*localMF[AREA_MF])[mfi];
  FArrayBox& areay=(*localMF[AREA_MF+1])[mfi];
  FArrayBox& areaz=(*localMF[AREA_MF+AMREX_SPACEDIM-1])[mfi];

  FArrayBox& xp=(*localMF[POTENTIAL_EDGE_MF])[mfi];
  FArrayBox& yp=(*localMF[POTENTIAL_EDGE_MF+1])[mfi];
  FArrayBox& zp=(*localMF[POTENTIAL_EDGE_MF+AMREX_SPACEDIM-1])[mfi];

  FArrayBox& gcell=(*localMF[POTENTIAL_FORCE_CELL_MF])[mfi];

  FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi];

  FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi]; // 1=fine/fine
  FArrayBox& maskcoef=(*localMF[MASKCOEF_MF])[mfi]; // 1=not cov

  int energyflag=0;
  int homflag=0; // default
  int local_project_option=0;
  int local_enable_spectral=enable_spectral;
  int use_VOF_weight=0;

  int ncomp_denold=1;
  int ncomp_veldest=gcell.nComp();
  int ncomp_dendest=gcell.nComp();

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // process_potential_force_cell 
  fort_mac_to_cell(
   &ns_time_order,
   &divu_outer_sweeps,
   &num_divu_outer_sweeps,
   &operation_flag, // 4 (gravity and surface tension force at cell)
   &energyflag,
   temperature_primitive_variable.dataPtr(),
   constant_density_all_time.dataPtr(),
   &nmat,
   &nparts,
   &nparts_def,
   im_solid_map_ptr,
   added_weight.dataPtr(),
   &nten,
   &level, 
   &finest_level,
   &local_project_option,
   &local_enable_spectral,
   &fluxvel_index,
   &fluxden_index,
   &facevel_index,
   &facecut_index,
   &icefacecut_index,
   &curv_index,
   &conservative_tension_force,
   &conservative_div_uu,
   filter_velocity.dataPtr(),
   &ignore_div_up,
   &pforce_index,
   &faceden_index,
   &icemask_index,
   &massface_index,
   &vofface_index,
   &ncphys,
   velbc.dataPtr(),
   presbc.dataPtr(), 
   &cur_time_slab, 
   &slab_step,
   &dt_slab, //calling fort_mac_to_cell
   xlo,dx,
   tilelo,tilehi,
   fablo,fabhi,
   &bfact,
   xp.dataPtr(),ARLIM(xp.loVect()),ARLIM(xp.hiVect()),
   yp.dataPtr(),ARLIM(yp.loVect()),ARLIM(yp.hiVect()),
   zp.dataPtr(),ARLIM(zp.loVect()),ARLIM(zp.hiVect()),
   xp.dataPtr(),ARLIM(xp.loVect()),ARLIM(xp.hiVect()), //xvel
   yp.dataPtr(),ARLIM(yp.loVect()),ARLIM(yp.hiVect()), //yvel
   zp.dataPtr(),ARLIM(zp.loVect()),ARLIM(zp.hiVect()), //zvel
   xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()),
   yface.dataPtr(),ARLIM(yface.loVect()),ARLIM(yface.hiVect()),
   zface.dataPtr(),ARLIM(zface.loVect()),ARLIM(zface.hiVect()),
   areax.dataPtr(),ARLIM(areax.loVect()),ARLIM(areax.hiVect()), //ax
   areay.dataPtr(),ARLIM(areay.loVect()),ARLIM(areay.hiVect()), //ay
   areaz.dataPtr(),ARLIM(areaz.loVect()),ARLIM(areaz.hiVect()), //az
   volfab.dataPtr(),ARLIM(volfab.loVect()),ARLIM(volfab.hiVect()),
   volfab.dataPtr(),ARLIM(volfab.loVect()),ARLIM(volfab.hiVect()), //rhs
   gcell.dataPtr(),ARLIM(gcell.loVect()),ARLIM(gcell.hiVect()), // veldest
   gcell.dataPtr(),ARLIM(gcell.loVect()),ARLIM(gcell.hiVect()), // dendest
   maskfab.dataPtr(), // 1=fine/fine  0=coarse/fine
   ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
   maskcoef.dataPtr(), // 1=not covered  0=covered
   ARLIM(maskcoef.loVect()),ARLIM(maskcoef.hiVect()),
   maskSEMfab.dataPtr(), 
   ARLIM(maskSEMfab.loVect()),ARLIM(maskSEMfab.hiVect()),
   levelpcfab.dataPtr(), //levelPC
   ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),
   solxfab.dataPtr(),ARLIM(solxfab.loVect()),ARLIM(solxfab.hiVect()),
   solyfab.dataPtr(),ARLIM(solyfab.loVect()),ARLIM(solyfab.hiVect()),
   solzfab.dataPtr(),ARLIM(solzfab.loVect()),ARLIM(solzfab.hiVect()),
   levelpcfab.dataPtr(),
   ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),//cterm
   presdenfab.dataPtr(0),  // HYDROSTATIC_PRESSURE
   ARLIM(presdenfab.loVect()),ARLIM(presdenfab.hiVect()), 
   presdenfab.dataPtr(1),  // HYDROSTATIC_DENSITY (denold)
   ARLIM(presdenfab.loVect()),ARLIM(presdenfab.hiVect()), 
   gcell.dataPtr(),ARLIM(gcell.loVect()),ARLIM(gcell.hiVect()), // ustar
   reconfab.dataPtr(),ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
   levelpcfab.dataPtr(),
   ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),//mdot
   levelpcfab.dataPtr(),
   ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),//maskdivres
   levelpcfab.dataPtr(),
   ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),//maskres
   &SDC_outer_sweeps,
   &homflag,
   &use_VOF_weight,
   &nsolve,
   &ncomp_denold,
   &ncomp_veldest,
   &ncomp_dendest,
   &SEM_advection_algorithm);

 } // mfi
} // omp
 ns_reconcile_d_num(152);

  // avgdown from level+1 to level.
 avgDown_localMF(POTENTIAL_FORCE_CELL_MF,0,AMREX_SPACEDIM,1);

}  // subroutine process_potential_force_cell


// typically ngrow=4 (see call in NavierStokes3.cpp)
void NavierStokes::metrics_dataALL(int ngrow) {

 int finest_level=parent->finestLevel();
 for (int ilev=0;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.metrics_data(ngrow);
 }
}


void NavierStokes::metrics_data_min_max_ALL(int caller_id) {

 std::fflush(NULL);

 int finest_level=parent->finestLevel();
 for (int ilev=0;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.metrics_data_min_max(caller_id);
 }
 std::fflush(NULL);
}


void NavierStokes::metrics_data(int ngrow) {
 
 bool use_tiling=ns_tiling;

 const Real* dx = geom.CellSize();

 delete_localMF_if_exist(VOLUME_MF,1);
 delete_localMF_if_exist(AREA_MF,AMREX_SPACEDIM);

 if (ngrow<0)
  amrex::Error("ngrow too small");

 new_localMF(VOLUME_MF,1,ngrow,-1);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  new_localMF(AREA_MF+dir,1,ngrow,dir);
 }

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 const BoxArray mfBA=localMF[VOLUME_MF]->boxArray();

 if (1==0) {
  for (int pid=0;pid<amrex::ParallelDescriptor::NProcs();pid++) {
   amrex::ParallelDescriptor::Barrier();
   if (amrex::ParallelDescriptor::MyProc()==pid) {
    std::fflush(NULL);
    std::cout << mfBA << '\n';
    std::cout << "NProcs()= " << 
     amrex::ParallelDescriptor::NProcs() << '\n';
    std::cout << "PROC= " << amrex::ParallelDescriptor::MyProc() << 
     " thread_class::nthreads= " << 
     thread_class::nthreads << '\n';
    std::fflush(NULL);
   }
  }  // pid=0..NProcs-1
 }

 double local_d_numPts=mfBA.d_numPts();
 thread_class::init_d_numPts(local_d_numPts);

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*localMF[VOLUME_MF],use_tiling); mfi.isValid();++mfi) {
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

  FArrayBox& volfab=(*localMF[VOLUME_MF])[mfi];
  FArrayBox& areax=(*localMF[AREA_MF])[mfi];
  FArrayBox& areay=(*localMF[AREA_MF+1])[mfi];
  FArrayBox& areaz=(*localMF[AREA_MF+AMREX_SPACEDIM-1])[mfi];

  int rzflag=0;
  if (geom.IsRZ())
   rzflag=1;
  else if (geom.IsCartesian())
   rzflag=0;
  else if (geom.IsCYLINDRICAL())
   rzflag=3;
  else
   amrex::Error("CoordSys bust 21");

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  FORT_METRICS(
   xlo,dx,
   areax.dataPtr(),ARLIM(areax.loVect()),ARLIM(areax.hiVect()),
   areay.dataPtr(),ARLIM(areay.loVect()),ARLIM(areay.hiVect()),
   areaz.dataPtr(),ARLIM(areaz.loVect()),ARLIM(areaz.hiVect()),
   volfab.dataPtr(),ARLIM(volfab.loVect()),ARLIM(volfab.hiVect()),
   tilelo,tilehi,
   fablo,fabhi,
   &bfact,
   &level,
   &ngrow,&rzflag);
 }  // mfi
}  // omp
 ns_reconcile_d_num(153);

} // subroutine metrics_data


void NavierStokes::metrics_data_min_max(int caller_id) {

 int ngrow=localMF_grow[VOLUME_MF];
 std::cout << "metrics_data_min_max caller_id " << caller_id <<'\n';
 std::cout << "metrics_data_min_max ngrow,level " << ngrow << ' ' << 
	 level <<'\n';

 std::cout << "volume_mf min= " << localMF[VOLUME_MF]->min(0,ngrow) << '\n'; 
 std::cout << "volume_mf max= " << localMF[VOLUME_MF]->max(0,ngrow) << '\n'; 
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  std::cout << "area_mf dir,min= " << dir << ' ' <<
	  localMF[AREA_MF+dir]->min(0,ngrow) << '\n'; 
  std::cout << "area_mf dir,max= " << dir << ' ' <<
	  localMF[AREA_MF+dir]->max(0,ngrow) << '\n'; 
 }

} // subroutine metrics_data_min_max


void NavierStokes::prescribe_solid_geometryALL(Real time,
  int renormalize_only,int local_truncate) {

 if (level!=0)
  amrex::Error("level should be 0 in prescribe_solid_geometryALL");
 int finest_level=parent->finestLevel();
 int nmat=num_materials;

 if (local_truncate==1) {
  Vector<Real> delta_mass_all;
  delta_mass_all.resize(nmat);

  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   if (ilev<finest_level) {
    ns_level.avgDown(LS_Type,0,nmat,0);
    ns_level.MOFavgDown();
   }
   ns_level.truncate_VOF(delta_mass_all);
   if (verbose>0) {
    if (ParallelDescriptor::IOProcessor()) {
     for (int im=0;im<nmat;im++) {
      std::cout << "truncate statistics: im,delta_mass " << im << ' ' <<
       delta_mass_all[im] << '\n';
     } // im
    } // IOProc?
   } // verbose>0?

  } // ilev

 } else if (local_truncate==0) {
  // do nothing
 } else
  amrex::Error("local_truncate invalid");

 if (renormalize_only==0) {
  if (std::abs(time-cur_time_slab)>1.0e-8)
   amrex::Error("prescribe solid at the new time");

  init_FSI_GHOST_MAC_MF_ALL(3);
 
 } else if (renormalize_only==1) {
  // do nothing
 } else
  amrex::Error("renormalize_only invalid");

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  if (ilev<finest_level) {
   ns_level.avgDown(LS_Type,0,nmat,0);
   ns_level.MOFavgDown();
  }
  ns_level.prescribe_solid_geometry(time,renormalize_only);
 }

} // end subroutine prescribe_solid_geometryALL


// 1. renormalize variables
// 2. extend from F>0 fluid regions into F=0 regions
// 3. if renormalize_only==0, 
//    a. init F,X,LS for the solid materials.
//    b. init U,T in the solid regions.
//    c. extrapolate F,X,LS from fluid regions into solid regions.
//   
// called from:
//  NavierStokes::prepare_post_process -> called from post_init_state 
//  NavierStokes::prepare_post_process -> called from post_restart
//    renormalize_only=0 when called from prepare_post_process.
//  NavierStokes::nonlinear_advection
//  NavierStokes::advance
//    renormalize_only=0 when called from advance.
//
void NavierStokes::prescribe_solid_geometry(Real time,int renormalize_only) {
 
 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "begin subroutine prescribe_solid_geometry() level= " <<
     level << '\n';
   std::fflush(NULL);
  }
 }
 bool use_tiling=ns_tiling;

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int finest_level = parent->finestLevel();

 MultiFab &S_new = get_new_data(State_Type,slab_step+1);
 MultiFab &LS_new = get_new_data(LS_Type,slab_step+1);
 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;
 int scomp_mofvars=(AMREX_SPACEDIM+1)+nmat*num_state_material;
 int dencomp=(AMREX_SPACEDIM+1);

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,1,6001);

 const Real* dx = geom.CellSize();

  // renormalize_only==1:
  //   project so that sum F_m_fluid=1
  // renormalize_only==0:
  //   correct F_m according to prescribed solid interface.
  //   project so that sum F_m_fluid=1 

 // nparts x (velocity + LS + temperature + flag)
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

 Vector<int> num_LS_extrap;
 num_LS_extrap.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  num_LS_extrap[tid]=0;
 }

 int num_LS_extrap_iter=1;

 if (renormalize_only==1) {
  // do nothing
 } else if (renormalize_only==0) {
  num_LS_extrap_iter=3;
 } else
  amrex::Error("renormalize_only invalid");	 

 for (int LS_extrap_iter=0;LS_extrap_iter<num_LS_extrap_iter; 
      LS_extrap_iter++) {

  if ((num_LS_extrap[0]==0)&&(LS_extrap_iter>=1)) {
	  // do nothing
  } else if ((num_LS_extrap[0]>=1)||(LS_extrap_iter==0)) {

   for (int tid=0;tid<thread_class::nthreads;tid++) {
    num_LS_extrap[tid]=0;
   }

   MultiFab* veldata=getState(1,0,(AMREX_SPACEDIM+1),time); 
   MultiFab* mofdata=getState(1,scomp_mofvars,nmat*ngeom_raw,time);
   MultiFab* dendata=getStateDen(1,time);
   MultiFab* lsdata=getStateDist(ngrow_distance,time,18);

   for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
    if (localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0)
     amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0");
    if (localMF[FSI_GHOST_MAC_MF+data_dir]->nComp()!=nparts_def*AMREX_SPACEDIM)
     amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nComp() bad");
   }

   if (veldata->nComp()!=(AMREX_SPACEDIM+1))
    amrex::Error("veldata incorrect ncomp");
   if (dendata->nComp()!=nmat*num_state_material)
    amrex::Error("dendata incorrect ncomp");
   if (mofdata->nComp()!=nmat*ngeom_raw)
    amrex::Error("mofdata incorrect ncomp");
   if (lsdata->nComp()!=nmat*(1+AMREX_SPACEDIM))
    amrex::Error("lsdata incorrect ncomp");
   if (lsdata->nGrow()!=ngrow_distance)
    amrex::Error("lsdata->nGrow()!=ngrow_distance");
   if (ngrow_distance!=4)
    amrex::Error("ngrow_distance invalid");

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

    FArrayBox& vofnew=S_new[mfi];
    FArrayBox& velnew=S_new[mfi];
    FArrayBox& dennew=S_new[mfi];
    FArrayBox& lsnew=LS_new[mfi];

    FArrayBox& solxfab=(*localMF[FSI_GHOST_MAC_MF])[mfi];
    FArrayBox& solyfab=(*localMF[FSI_GHOST_MAC_MF+1])[mfi];
    FArrayBox& solzfab=(*localMF[FSI_GHOST_MAC_MF+AMREX_SPACEDIM-1])[mfi];

     // mask=tag if not covered by level+1 or outside the domain.
    FArrayBox& maskcov=(*localMF[MASKCOEF_MF])[mfi];
    FArrayBox& lsfab=(*lsdata)[mfi];
    FArrayBox& moffab=(*mofdata)[mfi];
    FArrayBox& denfab=(*dendata)[mfi];
    FArrayBox& velfab=(*veldata)[mfi];

    const Real* xlo = grid_loc[gridno].lo();

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     // in: LEVELSET_3D.F90
    FORT_RENORMALIZE_PRESCRIBE(
      &tid_current,
      &level,&finest_level,
      &time,
      tilelo,tilehi,
      fablo,fabhi,&bfact,
      vofnew.dataPtr(scomp_mofvars),
      ARLIM(vofnew.loVect()),ARLIM(vofnew.hiVect()),
      solxfab.dataPtr(),
      ARLIM(solxfab.loVect()),ARLIM(solxfab.hiVect()),
      solyfab.dataPtr(),
      ARLIM(solyfab.loVect()),ARLIM(solyfab.hiVect()),
      solzfab.dataPtr(),
      ARLIM(solzfab.loVect()),ARLIM(solzfab.hiVect()),
      maskcov.dataPtr(),
      ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
      lsfab.dataPtr(),
      ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
      moffab.dataPtr(),
      ARLIM(moffab.loVect()),ARLIM(moffab.hiVect()),
      denfab.dataPtr(),
      ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
      velfab.dataPtr(),
      ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
      velnew.dataPtr(),
      ARLIM(velnew.loVect()),ARLIM(velnew.hiVect()),
      dennew.dataPtr(dencomp),
      ARLIM(dennew.loVect()),ARLIM(dennew.hiVect()),
      lsnew.dataPtr(),
      ARLIM(lsnew.loVect()),ARLIM(lsnew.hiVect()),
      xlo,dx,
      &cur_time_slab,
      &nmat,
      &nten,
      &nparts,
      &nparts_def,
      im_solid_map_ptr,
      &renormalize_only, 
      &solidheat_flag,
      &num_LS_extrap[tid_current],
      &num_LS_extrap_iter,
      &LS_extrap_iter,
      &ngrow_distance,
      constant_density_all_time.dataPtr());

   }  // mfi
} // omp
   ns_reconcile_d_num(154);

   for (int tid=1;tid<thread_class::nthreads;tid++) {
    num_LS_extrap[0]+=num_LS_extrap[tid];
   } // tid
   ParallelDescriptor::ReduceIntSum(num_LS_extrap[0]);

   delete veldata;
   delete mofdata;
   delete dendata;
   delete lsdata;

  } else
   amrex::Error("num_LS_extrap or LS_extrap_iter invalid");

 } // LS_extrap_iter=0..num_LS_extrap_iter-1

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "end subroutine prescribe_solid_geometry() level= " <<
     level << '\n';
   std::fflush(NULL);
  }
 }

}  // end subroutine prescribe_solid_geometry()

void NavierStokes::move_particles(
  AmrParticleContainer<N_EXTRA_REAL,0,0,0>& localPC_no_nbr,
  NeighborParticleContainer<N_EXTRA_REAL,0>& localPC_nbr) {

 bool use_tiling=ns_tiling;
 int max_level = parent->maxLevel();
 int finest_level=parent->finestLevel();

 if ((level>=0)&&(level<=finest_level)) {
  // do nothing
 } else 
  amrex::Error("level invalid");

 if ((level>=0)&&(level<=max_level)) {
  // do nothing
 } else 
  amrex::Error("level invalid");

 if ((level>=0)&&(level<=max_level_for_use)) {
  // do nothing
 } else 
  amrex::Error("level invalid");

 int nmat=num_materials;
 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 if (divu_outer_sweeps==0)
  vel_time_slab=prev_time_slab;
 else if (divu_outer_sweeps>0)
  vel_time_slab=cur_time_slab;
 else
  amrex::Error("divu_outer_sweeps invalid move_particles");

 int nnbr=particle_interaction_ngrow;
 if (nnbr>=1) {
  // do nothing
 } else
  amrex::Error("nnbr invalid");

 if (particles_flag==1) {

  const Real* dx = geom.CellSize();
  const Box& domain = geom.Domain();
  const int* domlo = domain.loVect();
  const int* domhi = domain.hiVect();

  int scomp_mofvars=AMREX_SPACEDIM+1+nmat*num_state_material;
  Vector<int> dombc(2*AMREX_SPACEDIM);
  const BCRec& descbc = get_desc_lst()[State_Type].getBC(scomp_mofvars);
  const int* b_rec=descbc.vect();
  for (int m=0;m<2*AMREX_SPACEDIM;m++)
   dombc[m]=b_rec[m];

   // level set function(s) prior to CLSMOF advection.
  MultiFab* LSmf=getStateDist(2,cur_time_slab,7);  
  if (LSmf->nComp()!=nmat*(1+AMREX_SPACEDIM))
   amrex::Error("LSmf invalid ncomp");
  if (LSmf->nGrow()!=2)
   amrex::Error("LSmf->nGrow()!=2");

  MultiFab* mac_velocity[AMREX_SPACEDIM];
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   mac_velocity[dir]=getStateMAC(Umac_Type,2,dir,0,1,vel_time_slab);
  }

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

   FArrayBox& lsfab=(*LSmf)[mfi];
   FArrayBox& xvelfab=(*mac_velocity[0])[mfi];
   FArrayBox& yvelfab=(*mac_velocity[1])[mfi];
   FArrayBox& zvelfab=(*mac_velocity[AMREX_SPACEDIM-1])[mfi];

   const Real* xlo = grid_loc[gridno].lo();

   auto& particles = localPC_no_nbr.GetParticles(level)
     [std::make_pair(mfi.index(),mfi.LocalTileIndex())];
   auto& particles_AoS = particles.GetArrayOfStructs();
   int Np=particles_AoS.size();

   auto& particles_NBR = localPC_nbr.GetParticles(level)
     [std::make_pair(mfi.index(),mfi.LocalTileIndex())];
   auto& particles_AoS_NBR = particles_NBR.GetArrayOfStructs();
   int Np_NBR=particles_AoS_NBR.size();

     // ParticleVector&
   auto& neighbors_local = 
     localPC_nbr.GetNeighbors(level,mfi.index(),mfi.LocalTileIndex());
   int Nn=neighbors_local.size();

    // component 1: number of particles linked to the cell.
    // component 2: the link to the list of particles.
    // Declare cell_particle_count local to the MFIter loop
    // so that this routine is thread safe.
   Box tilebox_grow=grow(tilegrid,nnbr);
   BaseFab<int> cell_particle_count(tilebox_grow,2);
   cell_particle_count.setVal(0);

    // The link index will start at 1.
   Vector< int > particle_link_data;
    // i_particle_link_1,i1,j1,k1,   (child link, parent link)
    // i_particle_link_2,i2,j2,k2,  ...
   particle_link_data.resize(Np_NBR*(1+AMREX_SPACEDIM));

   for (int i_link=0;i_link<Np_NBR*(1+AMREX_SPACEDIM);i_link++)
    particle_link_data[i_link]=0;

   int single_particle_size=AMREX_SPACEDIM+N_EXTRA_REAL;

   int dcomp=AMREX_SPACEDIM+1;
   Vector<int> denbc=getBCArray(State_Type,gridno,dcomp,
      nmat*num_state_material);
   Vector<int> velbc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     // declared in: LEVELSET_3D.F90
   fort_move_particle_container( 
     &tid_current,
     &single_particle_size,
     &particle_volume,
     &particle_relaxation_time_to_fluid,
     &particle_interaction_ngrow,
     &nmat,
     tilelo,tilehi,
     fablo,fabhi,&bfact,
     &level,
     &finest_level,
     xlo,dx,
     particles_AoS.data(),
     Np,  // pass by value
     particles_AoS_NBR.data(),
     Np_NBR,  // pass by value
     neighbors_local.data(),
     Nn,       //pass by value
     particle_link_data.dataPtr(),
     cell_particle_count.dataPtr(),
     ARLIM(cell_particle_count.loVect()),
     ARLIM(cell_particle_count.hiVect()),
     &dt_slab,
     &vel_time_slab,
     xvelfab.dataPtr(),
     ARLIM(xvelfab.loVect()),ARLIM(xvelfab.hiVect()),
     yvelfab.dataPtr(),
     ARLIM(yvelfab.loVect()),ARLIM(yvelfab.hiVect()),
     zvelfab.dataPtr(),
     ARLIM(zvelfab.loVect()),ARLIM(zvelfab.hiVect()),
     lsfab.dataPtr(),
     ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
     velbc.dataPtr(),
     denbc.dataPtr(),
     dombc.dataPtr(),
     domlo,domhi);

  }  // mfi
} // omp
  ns_reconcile_d_num(154);

  delete LSmf;
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   delete mac_velocity[dir];
  }

 } else
  amrex::Error("expecting NS_ncomp_particles>0");

} // end subroutine move_particles

// called from NavierStokes::prescribe_solid_geometryALL
void NavierStokes::truncate_VOF(Vector<Real>& delta_mass_all) {

 
 bool use_tiling=ns_tiling;

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int finest_level = parent->finestLevel();

 MultiFab &S_new = get_new_data(State_Type,slab_step+1);
 MultiFab& LS_new = get_new_data(LS_Type,slab_step+1);
 int nmat=num_materials;
 int scomp_mofvars=(AMREX_SPACEDIM+1)+
  nmat*num_state_material;
 int nc=scomp_mofvars+nmat*ngeom_raw+1;
 if (nc!=S_new.nComp())
  amrex::Error("nc invalid in truncate_VOF");
 if (LS_new.nComp()!=nmat*(AMREX_SPACEDIM+1))
  amrex::Error("LS_new ncomp invalid");

 if (delta_mass_all.size()!=nmat)
  amrex::Error("delta_mass_all has invalid size");

 Vector< Vector<Real> > local_delta_mass;
 local_delta_mass.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  local_delta_mass[tid].resize(nmat); 
  for (int im=0;im<nmat;im++)
   local_delta_mass[tid][im]=0.0;
 } // tid

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,1,80);

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

   FArrayBox& vofnew=S_new[mfi];
   FArrayBox& lsfab=LS_new[mfi];
    // mask=tag if not covered by level+1 or outside the domain.
   FArrayBox& maskcov=(*localMF[MASKCOEF_MF])[mfi];

   const Real* xlo = grid_loc[gridno].lo();

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

// trunate_volume_fractions:
// default: tessellating fluid => default==1
//          non-tesselating or tesselating solid => default==0
//  in: LEVELSET_3D.F90
   FORT_PURGEFLOTSAM(
     local_delta_mass[tid_current].dataPtr(),
     truncate_volume_fractions.dataPtr(),
     &truncate_thickness, 
     &level,&finest_level,
     &cur_time_slab,
     tilelo,tilehi,
     fablo,fabhi,&bfact,
     maskcov.dataPtr(),
     ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
     vofnew.dataPtr(scomp_mofvars),
     ARLIM(vofnew.loVect()),ARLIM(vofnew.hiVect()),
     lsfab.dataPtr(),
     ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
     xlo,dx,&nmat);

 }  // mfi
} // omp
 ns_reconcile_d_num(155);

 for (int tid=1;tid<thread_class::nthreads;tid++) {
  for (int im=0;im<nmat;im++) {
   local_delta_mass[0][im]+=local_delta_mass[tid][im];
  }
 } // tid

 ParallelDescriptor::Barrier();
 for (int im=0;im<nmat;im++) {
  ParallelDescriptor::ReduceRealSum(local_delta_mass[0][im]);
 }
 for (int im=0;im<nmat;im++) {
  delta_mass_all[im]+=local_delta_mass[0][im];
 }

}  // truncate_VOF()




void NavierStokes::output_triangles() {

 int finest_level = parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid output_triangles");
 NavierStokes& ns_level0=getLevel(0);

 bool use_tiling=ns_tiling;

  // mask=tag if not covered by level+1 or outside the domain.
 Real tag=1.0;
 int clearbdry=0;
 MultiFab* maskplot=maskfiner(1,tag,clearbdry);

   // vof,refcentroid,order,slope,intercept x nmat

 VOF_Recon_resize(1,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,1,610);

 const Real* dx = geom.CellSize();

 int nmat=num_materials;

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(localMF[SLOPE_RECON_MF]->boxArray().d_numPts());

// cannot do openmp here until each thread has its own file number.
// So, we remove the ifdef OPENMP clause ("omp parallel"), but we
// still need tiling==true if configured as such so that the
// Particle tiling structure is consistent with the Eulerian tiling
// structure.
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

   const Real* xlo = grid_loc[gridno].lo();

   FArrayBox& maskfab=(*maskplot)[mfi];
   FArrayBox& reconfab=(*localMF[SLOPE_RECON_MF])[mfi];

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // in: MARCHING_TETRA_3D.F90
   fort_isogrid(
    &tid_current,
    &visual_tessellate_vfrac,  // =0,1, or 3
    reconfab.dataPtr(),
    ARLIM(reconfab.loVect()),ARLIM(reconfab.hiVect()),
    xlo,dx,
    maskfab.dataPtr(),
    ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    &level,&gridno,&nmat);

   if (particles_flag==0) {
    // do nothing
   } else if (particles_flag==1) {

    int ipart=0;
    AmrParticleContainer<N_EXTRA_REAL,0,0,0>& localPC=
     ns_level0.get_new_dataPC(State_Type,slab_step+1,ipart);

    auto& particles = localPC.GetParticles(level)
      [std::make_pair(mfi.index(),mfi.LocalTileIndex())];
    auto& particles_AoS = particles.GetArrayOfStructs();

    int Np=particles_AoS.size();

     // declared in: NAVIERSTOKES_3D.F90
    fort_particle_grid(
      &tid_current,
      xlo,dx,
      particles_AoS.data(),
      Np,       //pass by value
      tilelo,tilehi,
      fablo,fabhi,&bfact,
      &level,&gridno,&ipart);
   } else
    amrex::Error("particles_flag invalid");

 }  // mfi
 ns_reconcile_d_num(156);

 delete maskplot;

}  // output_triangles

// dir=-1,0,1,2, or 3.
void NavierStokes::tecplot_debug(FArrayBox& fabdata,
 const Real* xlo,const int* fablo,const int* fabhi,
 const Real* dx,int dir,int id,int comp,int ncomp,
 int interior_only) {

 if (fabdata.nComp()<comp+ncomp)
  amrex::Error("ncomp invalid in tecplot_debug");
 int nsteps=parent->levelSteps(0);

 const int* growlo=fabdata.loVect();
 const int* growhi=fabdata.hiVect();
 int bfact=parent->Space_blockingFactor(level);
 fort_tecplotfab(
  &cur_time_slab,
  fabdata.dataPtr(comp),ARLIM(fabdata.loVect()),ARLIM(fabdata.hiVect()),
  growlo,growhi,
  fablo,fabhi,
  &bfact,
  xlo,dx,
  &dir,&ncomp,&interior_only, 
  &nsteps);

 std::cout << "OUTPUT is fabdata??????.plt \n";
 std::cout << "xlo= " << xlo[0] << ' ' << xlo[1] << ' ' <<
  xlo[AMREX_SPACEDIM-1] << '\n';
 std::cout << "dx= " << dx[0] << ' ' << dx[1] << ' ' <<
  dx[AMREX_SPACEDIM-1] << '\n';
 std::cout << "cur_time_slab= " << cur_time_slab << ' ' << "nsteps= " 
   << nsteps << '\n';
 std::cout << "comp= " << comp << '\n';
 std::cout << "ncomp= " << ncomp << '\n';
 std::cout << "dir= " << dir << " id= " << id << '\n';
 int fake_input;
 if (1==1) {
  std::cout << "PRESS <1> THEN <ENTER> TO CONTINUE \n";
  std::cin >> fake_input;
 }
}



// datatype=0 normal
// datatype=1 tensor face
// datatype=2 tensor cell
void NavierStokes::check_for_NAN_TENSOR_base(int datatype,MultiFab* mf,
  int sc,int dir,int id) {

 int finest_level=parent->finestLevel();
 int ncomp=mf->nComp();
 int ngrow=mf->nGrow();
 const BoxArray mfBA=mf->boxArray();
 int ngrid=mfBA.size();
 const Box& domain = geom.Domain();

 if (ngrow!=1)
  amrex::Error("ngrow invalid");

 std::fflush(NULL);

 if (mf->contains_nanTENSOR(datatype,sc,dir)==true) {
  std::cout << "id= " << id << '\n';
  std::cout << "sc= " << sc << '\n';
  std::cout << "dir= " << dir << '\n';
  std::cout << "ncomp= " << ncomp << '\n';
  std::cout << "ngrow= " << ngrow << '\n';
  std::cout << "level= " << level << '\n';
  std::cout << "finest_level= " << finest_level << '\n';
  std::cout << "domain= " << domain << '\n';
  std::cout << "ngrid= " << ngrid << '\n';
  std::cout << "mfBA= " << mfBA << '\n';
  amrex::Error("mf contains nan");
 }
 if (mf->contains_infTENSOR(datatype,sc,dir)==true) {
  std::cout << "id= " << id << '\n';
  std::cout << "sc= " << sc << '\n';
  std::cout << "dir= " << dir << '\n';
  std::cout << "ncomp= " << ncomp << '\n';
  std::cout << "ngrow= " << ngrow << '\n';
  std::cout << "level= " << level << '\n';
  std::cout << "finest_level= " << finest_level << '\n';
  std::cout << "domain= " << domain << '\n';
  std::cout << "ngrid= " << ngrid << '\n';
  std::cout << "mfBA= " << mfBA << '\n';
  amrex::Error("mf contains inf");
 }
 int force_check=1;
 int ncomp_tensor=1;
 int ngrow_tensor=0;
 Real warning_cutoff=1.0e+99;
 aggressive_debug(
  datatype,
  force_check,
  mf,
  sc,ncomp_tensor,
  ngrow_tensor,
  dir,id,
  warning_cutoff);

 std::fflush(NULL);

} // subroutine check_for_NAN_TENSOR_base




// datatype=0 normal
// datatype=1 face grad U
// datatype=2 cell grad U
void NavierStokes::check_for_NAN_TENSOR(int datatype,MultiFab* mf,int id) {

 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid check_for_NAN_TENSOR");

 int ntensor=AMREX_SPACEDIM*AMREX_SPACEDIM; 

 int ncomp=mf->nComp();
 int ngrow=mf->nGrow();

 if (ncomp!=ntensor)
  amrex::Error("ncomp invalid");
 if (ngrow!=1)
  amrex::Error("ngrow invalid");

 std::fflush(NULL);

 for (int sc=0;sc<ntensor;sc++) {

   int dir=0;
   if ((sc>=0)&&(sc<AMREX_SPACEDIM)) // ux,vx,wx
    dir=0;
   else if (sc<2*AMREX_SPACEDIM)  // uy,vy,wy
    dir=1;
   else if ((sc<AMREX_SPACEDIM*AMREX_SPACEDIM)&&
	    (AMREX_SPACEDIM==3)) // uz,vz,wz
    dir=AMREX_SPACEDIM-1;
   else
    amrex::Error("sc invalid");

   check_for_NAN_TENSOR_base(datatype,mf,sc,dir,id);
 }  // sc

 
 std::fflush(NULL);

} // subroutine check_for_NAN_TENSOR


void NavierStokes::check_for_NAN(MultiFab* mf,int id) {

 int finest_level=parent->finestLevel();
 int scomp=0;
 int ncomp=mf->nComp();
 int ngrow=mf->nGrow();
 const BoxArray mfBA=mf->boxArray();
 int ngrid=mfBA.size();
 const Box& domain = geom.Domain();

 for (int sc=0;sc<ncomp;sc++) {
  for (int ng=0;ng<=ngrow;ng++) {
   if (mf->contains_nan(sc,1,ng)==true) { // scomp,ncomp,ngrow
    std::cout << "id= " << id << '\n';
    std::cout << "sc= " << sc << '\n';
    std::cout << "ng= " << ng << '\n';
    std::cout << "ncomp= " << ncomp << '\n';
    std::cout << "ngrow= " << ngrow << '\n';
    std::cout << "level= " << level << '\n';
    std::cout << "finest_level= " << finest_level << '\n';
    std::cout << "domain= " << domain << '\n';
    std::cout << "ngrid= " << ngrid << '\n';
    std::cout << "mfBA= " << mfBA << '\n';
    amrex::Error("mf contains nan");
   }
   ParallelDescriptor::Barrier();

   if (mf->contains_inf(sc,1,ng)==true) {
    std::cout << "id= " << id << '\n';
    std::cout << "sc= " << sc << '\n';
    std::cout << "ng= " << ng << '\n';
    std::cout << "ncomp= " << ncomp << '\n';
    std::cout << "ngrow= " << ngrow << '\n';
    std::cout << "level= " << level << '\n';
    std::cout << "finest_level= " << finest_level << '\n';
    std::cout << "domain= " << domain << '\n';
    std::cout << "ngrid= " << ngrid << '\n';
    std::cout << "mfBA= " << mfBA << '\n';
    amrex::Error("mf contains inf");
   }
   ParallelDescriptor::Barrier();
  }  // ng=0..ngrow
 }  // sc=0..ncomp-1

 std::fflush(NULL);

 int datatype=0;
 int force_check=1;
 int dir=-1;
 Real warning_cutoff=1.0e+99;
 aggressive_debug(
  datatype,
  force_check,
  mf,
  scomp,ncomp,
  ngrow,
  dir,id,
  warning_cutoff);  

 std::fflush(NULL);
} // subroutine check_for_NAN

void NavierStokes::output_zones(
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
   MultiFab* magtracemf,
   int& grids_per_level,
   BoxArray& cgrids_minusBA,
   Real* slice_data,
   int do_plot,int do_slice) {

 const Real* dx = geom.CellSize();
 const Real* prob_lo = geom.ProbLo();
 const Real* prob_hi = geom.ProbHi();

 int nmat=num_materials;
 int nden=nmat*num_state_material;

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

  // x,u,p,den,T,Y1..Yn,mag vort,LS
 if (visual_ncomp==2*AMREX_SPACEDIM+3+num_species_var+1+nmat) {
  // do nothing
 } else
  amrex::Error("visual_ncomp invalid");

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

 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
  if (localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0)
   amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0");
  if (localMF[FSI_GHOST_MAC_MF+data_dir]->nComp()!=nparts_def*AMREX_SPACEDIM)
   amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nComp() bad");
 }

 VOF_Recon_resize(1,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,1,615);

 debug_ngrow(MASKSEM_MF,1,28);

 int rzflag=0;
 if (geom.IsRZ())
  rzflag=1;
 else if (geom.IsCartesian())
  rzflag=0;
 else if (geom.IsCYLINDRICAL())
  rzflag=3;
 else
  amrex::Error("CoordSys bust 22");

 int finest_level = parent->finestLevel();
 int tecplot_finest_level=finest_level;
 if (tecplot_max_level<tecplot_finest_level)
  tecplot_finest_level=tecplot_max_level;

 if (level<=tecplot_finest_level) {

  BoxArray cgrids(grids);
  BoxList cgrids_list(cgrids);

  if (level<tecplot_finest_level) {
   NavierStokes& ns_finer=getLevel(level+1);
   BoxArray fgrids(ns_finer.grids);
   fgrids.coarsen(2);
   const Box& domain=geom.Domain();
   BoxList fgrids_list(fgrids);
   BoxList fgrids_complement=amrex::complementIn(domain,fgrids_list);
   cgrids_list.intersect(fgrids_complement);
  } // level<tecplot_finest_level

  BoxArray cgrids_minusBA_temp(cgrids_list);
  grids_per_level=cgrids_minusBA_temp.size();
  cgrids_minusBA=cgrids_minusBA_temp;

  if (level==tecplot_finest_level) {
   if (grids.size()!=grids_per_level)
    amrex::Error("grids_per_level incorrect");
  }

  NavierStokes& ns_finest=getLevel(tecplot_finest_level);
  const Real* dxfinest = ns_finest.geom.CellSize();

   // x,y,z,xvel,yvel,zvel,PMG,PEOS,div,den,Temp,KE
   // (value of material with LS>0)
  int nslice=0;
  int nstate_slice=2*AMREX_SPACEDIM+6; 
 
  if ((slice_dir>=0)&&(slice_dir<AMREX_SPACEDIM)) {
   const Box& domain_finest = ns_finest.geom.Domain();
   const int* domlo_finest = domain_finest.loVect();
   const int* domhi_finest = domain_finest.hiVect();
   nslice=domhi_finest[slice_dir]-domlo_finest[slice_dir]+3;
  } else
   amrex::Error("slice_dir invalid");

  if (grids_per_level>0) {

   DistributionMapping cgrids_minus_map(cgrids_minusBA);

   MultiFab* maskSEM_minus=new MultiFab(cgrids_minusBA,cgrids_minus_map,1,0,
    MFInfo().SetTag("maskSEM_minus"),FArrayBoxFactory());

   MultiFab* velmfminus=new MultiFab(cgrids_minusBA,cgrids_minus_map,
    (AMREX_SPACEDIM+1),1,
    MFInfo().SetTag("velmfminus"),FArrayBoxFactory());

   MultiFab* vofmfminus=new MultiFab(cgrids_minusBA,cgrids_minus_map,
     nmat*ngeom_recon,1,
     MFInfo().SetTag("vofminus"),FArrayBoxFactory());

   MultiFab* presmfminus=new MultiFab(cgrids_minusBA,cgrids_minus_map,
     1,1,
     MFInfo().SetTag("presmfminus"),FArrayBoxFactory());

   MultiFab* divmfminus=new MultiFab(cgrids_minusBA,cgrids_minus_map,
     1,1,
     MFInfo().SetTag("divmfminus"),FArrayBoxFactory());

   MultiFab* div_data_minus=new MultiFab(cgrids_minusBA,cgrids_minus_map,
     1,1,
     MFInfo().SetTag("div_data_minus"),FArrayBoxFactory());

   MultiFab* denmfminus=new MultiFab(cgrids_minusBA,cgrids_minus_map,
     nden,1,
     MFInfo().SetTag("denmfminus"),FArrayBoxFactory());

   MultiFab* mom_denmfminus=new MultiFab(cgrids_minusBA,cgrids_minus_map,
     nmat,1,
     MFInfo().SetTag("mom_denmfminus"),FArrayBoxFactory());

   int elastic_ncomp=viscoelasticmf->nComp();
   if (elastic_ncomp==
       num_materials_viscoelastic*NUM_TENSOR_TYPE+AMREX_SPACEDIM) {
    // do nothing
   } else
    amrex::Error("elastic_ncomp invalid");

   MultiFab* viscoelasticmfminus=
    new MultiFab(cgrids_minusBA,cgrids_minus_map,
     elastic_ncomp,1,
     MFInfo().SetTag("viscoelasticmfminus"),FArrayBoxFactory());

   MultiFab* lsdistmfminus=
    new MultiFab(cgrids_minusBA,cgrids_minus_map,
	nmat*(AMREX_SPACEDIM+1),1,
        MFInfo().SetTag("lsdistmfminus"),FArrayBoxFactory());

   MultiFab* viscmfminus=new MultiFab(cgrids_minusBA,cgrids_minus_map,
     nmat,1,
     MFInfo().SetTag("viscmfminus"),FArrayBoxFactory());

   MultiFab* magtracemfminus=new MultiFab(cgrids_minusBA,cgrids_minus_map,
     5*nmat,1,
     MFInfo().SetTag("magtracemfminus"),FArrayBoxFactory());

   ParallelDescriptor::Barrier();

     // FabArray.H     
     // scomp,dcomp,ncomp,s_nghost,d_nghost
   maskSEM_minus->ParallelCopy(*localMF[MASKSEM_MF],0,0,
    1,0,0,geom.periodicity());

   check_for_NAN(localMF[MASKSEM_MF],1);
   check_for_NAN(maskSEM_minus,11);

     // FabArray.H     
     // scomp,dcomp,ncomp,s_nghost,d_nghost
   velmfminus->ParallelCopy(*velmf,0,0,
    (AMREX_SPACEDIM+1),1,1,geom.periodicity());

   check_for_NAN(velmf,1);
   check_for_NAN(velmfminus,11);
 
   vofmfminus->ParallelCopy(*localMF[SLOPE_RECON_MF],0,0,
    nmat*ngeom_recon,1,1,geom.periodicity());

   check_for_NAN(localMF[SLOPE_RECON_MF],2); // id==2
   check_for_NAN(vofmfminus,12);

   // scomp,dcomp,ncomp,sgrow,dgrow,period,op
   presmfminus->ParallelCopy(*presmf,0,0,1,
		   1,1,geom.periodicity()); 

   check_for_NAN(presmf,3);
   check_for_NAN(presmfminus,13);

   // scomp,dcomp,ncomp,sgrow,dgrow,period,op
   divmfminus->ParallelCopy(*divmf,0,0,1,
		   1,1,geom.periodicity()); 

   check_for_NAN(divmf,4);
   check_for_NAN(divmfminus,14);

   // scomp,dcomp,ncomp,sgrow,dgrow,period,op
   div_data_minus->ParallelCopy(*div_data,0,0,1,
		   1,1,geom.periodicity()); 

   check_for_NAN(div_data,5);
   check_for_NAN(div_data_minus,15);

   // scomp,dcomp,ncomp,sgrow,dgrow,period,op
   denmfminus->ParallelCopy(*denmf,0,0,nden,
		   1,1,geom.periodicity()); 

   // scomp,dcomp,ncomp,sgrow,dgrow,period,op
   mom_denmfminus->ParallelCopy(*mom_denmf,0,0,nmat,
		   1,1,geom.periodicity()); 

   check_for_NAN(denmf,6);
   check_for_NAN(denmfminus,16);
   check_for_NAN(mom_denmf,6);
   check_for_NAN(mom_denmfminus,16);

    // scomp,dcomp,ncomp,sgrow,dgrow,period,op
   viscoelasticmfminus->ParallelCopy(*viscoelasticmf,0,0,
     elastic_ncomp,
     1,1,geom.periodicity()); 
   check_for_NAN(viscoelasticmf,6);
   check_for_NAN(viscoelasticmfminus,16);

   // scomp,dcomp,ncomp,sgrow,dgrow,period,op
   lsdistmfminus->ParallelCopy(*lsdistmf,0,0,nmat*(1+AMREX_SPACEDIM),
     1,1,geom.periodicity()); 

   check_for_NAN(lsdistmf,7);
   check_for_NAN(lsdistmfminus,18);

   if (1==0) {
    int visc_ncomp=viscmf->nComp();
    int visc_ngrow=viscmf->nGrow();
    const BoxArray mfBA=viscmf->boxArray();
    std::cout << "viscmf ncomp= " << visc_ncomp << '\n';
    std::cout << "viscmf ngrow= " << visc_ngrow << '\n';
    std::cout << "viscmf BA= " << mfBA << '\n';
   }

   // scomp,dcomp,ncomp,sgrow,dgrow,period,op
   viscmfminus->ParallelCopy(*viscmf,0,0,nmat,
		   1,1,geom.periodicity()); 

   check_for_NAN(viscmf,9);
   check_for_NAN(viscmfminus,19);

   // scomp,dcomp,ncomp,sgrow,dgrow,period,op
   magtracemfminus->ParallelCopy(*magtracemf,0,0,5*nmat,
		   1,1,geom.periodicity()); 

   check_for_NAN(magtracemf,10);
   check_for_NAN(magtracemfminus,20);
 
   ParallelDescriptor::Barrier();

   int bfact=parent->Space_blockingFactor(level);

   if (thread_class::nthreads<1)
    amrex::Error("thread_class::nthreads invalid");
   thread_class::init_d_numPts(velmfminus->boxArray().d_numPts());

// cannot do openmp here until each thread has its own
// file handle.  Also, the update to visual_fab_output is not thread safe.
// MFIter (const FabArrayBase& fabarray,unsigned char flags_=0) ,
// is no tiling.
   for (MFIter mfi(*velmfminus,false); mfi.isValid(); ++mfi) {

    if (cgrids_minusBA[mfi.index()] != mfi.validbox())
     amrex::Error("cgrids_minusBA[mfi.index()] != mfi.validbox()");

    const Box& tilegrid = mfi.tilebox();

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    const int gridno = mfi.index();

    // cgrids_minusBA=grids at finest level.
    const Box& fabgrid = cgrids_minusBA[gridno];

    const int* lo=fabgrid.loVect();
    const int* hi=fabgrid.hiVect();

    for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

     if (lo[dir]%bfact!=0)
      amrex::Error("lo not divisible by bfact");
     if ((hi[dir]+1)%bfact!=0)
      amrex::Error("hi+1 not divisible by bfact");
 
    } // dir

    if (level==tecplot_finest_level) {
     if (cgrids_minusBA[gridno]!=grids[gridno])
      amrex::Error("box mismatch1");
    } // level=tecplot_finest_level

    FArrayBox& maskSEMfab=(*maskSEM_minus)[mfi];
    FArrayBox& velfab=(*velmfminus)[mfi];
    FArrayBox& voffab=(*vofmfminus)[mfi];
    FArrayBox& presfab=(*presmfminus)[mfi];
    FArrayBox& divfab=(*divmfminus)[mfi];
    FArrayBox& div_data_fab=(*div_data_minus)[mfi];
    FArrayBox& denfab=(*denmfminus)[mfi];
    FArrayBox& mom_denfab=(*mom_denmfminus)[mfi];
    FArrayBox& elasticfab=(*viscoelasticmfminus)[mfi];
    FArrayBox& lsdistfab=(*lsdistmfminus)[mfi];
    FArrayBox& viscfab=(*viscmfminus)[mfi];
    FArrayBox& magtracefab=(*magtracemfminus)[mfi];

      // in: NAVIERSTOKES_3D.F90
    fort_cellgrid(
     &tid_current,
     &bfact,
     visual_fab_output.dataPtr(),
     ARLIM(visual_fab_output.loVect()),
     ARLIM(visual_fab_output.hiVect()),
     visual_domain.loVect(),
     visual_domain.hiVect(),
     &visual_ncomp,
     maskSEMfab.dataPtr(),
     ARLIM(maskSEMfab.loVect()),ARLIM(maskSEMfab.hiVect()),
     velfab.dataPtr(),ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
     voffab.dataPtr(),ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),
     presfab.dataPtr(),ARLIM(presfab.loVect()),ARLIM(presfab.hiVect()),
     divfab.dataPtr(),ARLIM(divfab.loVect()),ARLIM(divfab.hiVect()),
     div_data_fab.dataPtr(),
     ARLIM(div_data_fab.loVect()),ARLIM(div_data_fab.hiVect()),
     denfab.dataPtr(),
     ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
     mom_denfab.dataPtr(),
     ARLIM(mom_denfab.loVect()),ARLIM(mom_denfab.hiVect()),
     elasticfab.dataPtr(),
     ARLIM(elasticfab.loVect()),ARLIM(elasticfab.hiVect()),
     lsdistfab.dataPtr(),ARLIM(lsdistfab.loVect()),ARLIM(lsdistfab.hiVect()),
     viscfab.dataPtr(),ARLIM(viscfab.loVect()),ARLIM(viscfab.hiVect()),
     magtracefab.dataPtr(),
     ARLIM(magtracefab.loVect()),ARLIM(magtracefab.hiVect()),
     prob_lo,
     prob_hi,
     dx,
     lo,hi,
     &level,
     &finest_level,
     &gridno,
     &visual_tessellate_vfrac,  // = 0,1,3
     &visual_option,
     &rzflag,
     &nmat,
     &nparts,
     &nparts_def,
     im_solid_map_ptr,
     &elastic_ncomp,
     slice_data,
     &nslice,
     &nstate_slice,&slice_dir,
     xslice.dataPtr(),
     dxfinest,
     &do_plot,&do_slice);
   }  // mfi
   ns_reconcile_d_num(157);

   delete viscoelasticmfminus;

   delete maskSEM_minus;
   delete velmfminus;
   delete vofmfminus;
   delete presmfminus;
   delete divmfminus;
   delete div_data_minus;
   delete denmfminus;
   delete mom_denmfminus;
   delete lsdistmfminus;
   delete viscmfminus;
   delete magtracemfminus;

  } else if (grids_per_level==0) {
  
   // do nothing

  } else 
   amrex::Error("grids_per_level is corrupt");

 } else if (level<=finest_level) {

  // do nothing

 } else {
  amrex::Error("level invalid output_zones");
 }

}  // subroutine output_zones


// data_dir=-1 cell centered data
// data_dir=0..sdim-1 face centered data.
// data_dir=sdim node data
void NavierStokes::Sanity_output_zones(
   int data_id,
   int data_dir,
   MultiFab* datamf,
   int ncomp,
   int& grids_per_level,
   BoxArray& cgrids_minusBA) {

 const Real* dx = geom.CellSize();
 const Real* prob_lo = geom.ProbLo();
 const Real* prob_hi = geom.ProbHi();

 int rzflag=0;
 if (geom.IsRZ())
  rzflag=1;
 else if (geom.IsCartesian())
  rzflag=0;
 else if (geom.IsCYLINDRICAL())
  rzflag=3;
 else
  amrex::Error("CoordSys bust 22");

 if (ParallelDescriptor::IOProcessor()) {
  std::cout << "in: Sanity_output_zones data_id=" << data_id <<
    " data_dir=" << data_dir << " ncomp=" << ncomp << '\n';
 }

 int finest_level = parent->finestLevel();
 int tecplot_finest_level=finest_level;
 if (tecplot_max_level<tecplot_finest_level)
  tecplot_finest_level=tecplot_max_level;

 if (level<=tecplot_finest_level) {

  BoxArray cgrids(grids);
  BoxList cgrids_list(cgrids);

  if (level<tecplot_finest_level) {
   NavierStokes& ns_finer=getLevel(level+1);
   BoxArray fgrids(ns_finer.grids);
   fgrids.coarsen(2);
   const Box& domain=geom.Domain();
   BoxList fgrids_list(fgrids);
   BoxList fgrids_complement=amrex::complementIn(domain,fgrids_list);
   cgrids_list.intersect(fgrids_complement);
  } // level<tecplot_finest_level

  BoxArray cgrids_minusBA_temp(cgrids_list);
  grids_per_level=cgrids_minusBA_temp.size();
  cgrids_minusBA=cgrids_minusBA_temp;

  if (level==tecplot_finest_level) {
   if (grids.size()!=grids_per_level)
    amrex::Error("grids_per_level incorrect");
  }

  if (grids_per_level>0) {

   DistributionMapping cgrids_minus_map(cgrids_minusBA);

   BoxArray minus_boxes(cgrids_minusBA);

   if (data_dir==-1) {
    // do nothing
   } else if ((data_dir>=0)&&(data_dir<AMREX_SPACEDIM)) {
    minus_boxes.surroundingNodes(data_dir);
   } else if (data_dir==AMREX_SPACEDIM) {
    minus_boxes.convert(IndexType::TheNodeType()); 
   } else
    amrex::Error("data_dir invalid sanity_output_zones");

   MultiFab* datamfminus=new MultiFab(minus_boxes,cgrids_minus_map,
    ncomp,0,
    MFInfo().SetTag("datamfminus"),FArrayBoxFactory());

   ParallelDescriptor::Barrier();

   debug_ixType_raw(datamf,data_dir,data_id);

     // FabArray.H     
     // scomp,dcomp,ncomp,s_nghost,d_nghost
   datamfminus->ParallelCopy(*datamf,0,0,
    ncomp,0,0,geom.periodicity());

   if (data_dir==-1) {
    check_for_NAN(datamf,1);
    check_for_NAN(datamfminus,11);
   } else if ((data_dir>=0)&&(data_dir<=AMREX_SPACEDIM)) {
    check_for_NAN(datamf,1);
    check_for_NAN(datamfminus,11);
   } else
    amrex::Error("data_dir invalid");
 
   ParallelDescriptor::Barrier();

   int bfact=parent->Space_blockingFactor(level);

   if (thread_class::nthreads<1)
    amrex::Error("thread_class::nthreads invalid");
   thread_class::init_d_numPts(minus_boxes.d_numPts());

// cannot do openmp here until each thread has its own
// file handle.  
// MFIter (const FabArrayBase& fabarray,unsigned char flags_=0) ,
// is no tiling.
   for (MFIter mfi(*datamfminus,false); mfi.isValid(); ++mfi) {

    if (minus_boxes[mfi.index()] != mfi.validbox())
     amrex::Error("minus_boxes[mfi.index()] != mfi.validbox()");

    const Box& tilegrid = mfi.tilebox();

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    const int gridno = mfi.index();

    // cgrids_minusBA=grids at finest level.
    const Box& fabgrid = cgrids_minusBA[gridno];

    const int* lo=fabgrid.loVect();
    const int* hi=fabgrid.hiVect();

    for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

     if (lo[dir]%bfact!=0)
      amrex::Error("lo not divisible by bfact");
     if ((hi[dir]+1)%bfact!=0)
      amrex::Error("hi+1 not divisible by bfact");
 
    } // dir

    if (level==tecplot_finest_level) {
     if (cgrids_minusBA[gridno]!=grids[gridno])
      amrex::Error("box mismatch1");
    } // level=tecplot_finest_level

    FArrayBox& datafab=(*datamfminus)[mfi];

     // in: NAVIERSTOKES_3D.F90
    fort_cellgrid_sanity(
     &tid_current,
     &data_dir,
     &bfact,
     &ncomp,
     datafab.dataPtr(),ARLIM(datafab.loVect()),ARLIM(datafab.hiVect()),
     prob_lo,
     prob_hi,
     dx,
     lo,hi,
     &level,
     &finest_level,
     &gridno,
     &rzflag);
   }  // mfi
   ns_reconcile_d_num(157);

   delete datamfminus;

  } else if (grids_per_level==0) {
  
   // do nothing

  } else 
   amrex::Error("grids_per_level is corrupt");

 } else if (level<=finest_level) {

  // do nothing

 } else {
  amrex::Error("level invalid Sanity_output_zones");
 }

}  // subroutine Sanity_output_zones


// spectral_override==0 => always low order
void NavierStokes::avgDown_localMF_ALL(int idxMF,int scomp,int ncomp,
  int spectral_override) {

 int finest_level = parent->finestLevel();
 if (level!=0)
  amrex::Error("only call with level=0");

 for (int i=finest_level-1;i>=level;i--) {
  NavierStokes& ns_level=getLevel(i);
  ns_level.avgDown_localMF(idxMF,scomp,ncomp,spectral_override);
 }
}  // subroutine avgDown_localMF_ALL

// spectral_override==0 => always low order
void NavierStokes::avgDownALL(int stateidx,int startcomp,int numcomp,
  int spectral_override) {

 int finest_level = parent->finestLevel();
 if (level!=0)
  amrex::Error("only call with level=0");

 for (int i=finest_level-1;i>=level;i--) {
  NavierStokes& ns_level=getLevel(i);
  ns_level.avgDown(stateidx,startcomp,numcomp,spectral_override);
 }
}  // subroutine avgDownALL

void NavierStokes::avgDownError_ALL() {
   
 int finest_level = parent->finestLevel();
 if (level!=0)
  amrex::Error("only call with level=0");

 for (int i=finest_level-1;i>=level;i--) {
  NavierStokes& ns_level=getLevel(i);
  ns_level.avgDownError();
 }
}  // subroutine avgDownError_ALL


void NavierStokes::copybc(Vector<int> dest,Vector<int> source,
                          int scomp,int dcomp,int ncomp) {

 int bcstatesize=2*AMREX_SPACEDIM;
 if ((dcomp+ncomp)*bcstatesize>dest.size())
  amrex::Error("dest bc too small");
 if ((scomp+ncomp)*bcstatesize>source.size()) 
  amrex::Error("source bc too small");

 for (int dir=0;dir<ncomp;dir++) {
  for (int ibc=0;ibc<bcstatesize;ibc++) {
   dest[(dir+dcomp)*bcstatesize+ibc]=source[(dir+scomp)*bcstatesize+ibc];
  }
 }

}  // copybc 


void NavierStokes::getStateALL(int ngrow,Real time,int scomp,
   int ncomp,int idx_localMF) {

 int finest_level=parent->finestLevel();
 if (level!=0)
  amrex::Error("level!=0 in getStateALL");
 if (level<=finest_level) {
  // do nothing
 } else
  amrex::Error("level or finest_level invalid");

 for (int i=finest_level;i>=level;i--) { 
  NavierStokes& ns_level=getLevel(i);
  ns_level.getState_localMF(idx_localMF,ngrow,scomp,ncomp,time);
 }

} // end subroutine getStateALL


void NavierStokes::negateALL(int idx_localMF) {

 int finest_level = parent->finestLevel();
 if (level!=0)
  amrex::Error("level!=0 in negateALL");
 int ngrow=localMF[idx_localMF]->nGrow();
 int ncomp=localMF[idx_localMF]->nComp();
 if (ncomp!=1)
  amrex::Error("ncomp invalid negateALL");
 for (int ilev=0;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.localMF[idx_localMF]->negate(0,1,ngrow);
 }
} // end subroutine negateALL

void NavierStokes::zeroALL(int ngrow,int ncomp,int idx_localMF) {

 int finest_level = parent->finestLevel();
 if (level!=0)
  amrex::Error("level!=0 in zeroALL");

 for (int ilev=0;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);

  if (ncomp>ns_level.localMF[idx_localMF]->nComp())
   amrex::Error("idx_localMF too small ncomp");
  if (ngrow>ns_level.localMF[idx_localMF]->nGrow()) 
   amrex::Error("idx_localMF too small ngrow");

  ns_level.localMF[idx_localMF]->setVal(0.0,0,ncomp,ngrow);
 }
} // end subroutine zeroALL

void NavierStokes::allocate_array(int ngrow,int ncomp,int grid_type,
  int idx_localMF) {

 int finest_level = parent->finestLevel();

 if (finest_level>=0) {

  if (level==0) {

   if ((grid_type>=-1)&&(grid_type<=5)) { //node,cell,x,y,z,...

    for (int i=finest_level;i>=level;i--) {
     NavierStokes& ns_level=getLevel(i);
      // initializes localMF[idx_localMF] to 0.0
     ns_level.new_localMF(idx_localMF,ncomp,ngrow,grid_type); 
    } // i

   } else
    amrex::Error("grid_type invalid");

  } else
   amrex::Error("level invalid");

 } else
  amrex::Error("finest_level invalid");

}  // subroutine allocate_array

void NavierStokes::allocate_independent_var(int nsolve,int idx) {

 if ((nsolve!=1)&&(nsolve!=AMREX_SPACEDIM))
  amrex::Error("nsolve invalid36");

 int finest_level = parent->finestLevel();

 if (level!=0)
  amrex::Error("level!=0 in allocate_independent_var");

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.new_localMF(idx,nsolve,1,-1); // sets values to 0.0
 } // ilev

}  // allocate_independent_var

void NavierStokes::allocate_rhs_var(int nsolve,int idx) {

 if ((nsolve!=1)&&(nsolve!=AMREX_SPACEDIM))
  amrex::Error("nsolve invalid36");

 int finest_level = parent->finestLevel();

 if (level!=0)
  amrex::Error("level!=0 in allocate_rhs_var");

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.new_localMF(idx,nsolve,0,-1); // sets values to 0.0
 } // ilev

}  // allocate_rhs_var

void NavierStokes::setVal_localMF(int idx,Real dataval,
 int scomp,int ncomp,int ngrow) {

 localMF[idx]->setVal(dataval,scomp,ncomp,ngrow);

} // setVal_localMF

void NavierStokes::setVal_array(int ngrow,int ncomp,Real dataval,
  int idx_localMF) {

 int finest_level = parent->finestLevel();

 if (finest_level>=0) {

  if (level==0) {

   for (int i=finest_level;i>=level;i--) {
    NavierStokes& ns_level=getLevel(i);
    ns_level.localMF[idx_localMF]->setVal(dataval,0,ncomp,ngrow);
   }

  } else
   amrex::Error("level invalid");

 } else
  amrex::Error("finest_level invalid");

} // subroutine setVal_array 


void NavierStokes::mult_array(int ngrow,int ncomp,Real dataval,
  int idx_localMF) {

 int finest_level = parent->finestLevel();
 if (level!=0)
  amrex::Error("level!=0 in setVal_array");

 for (int i=finest_level;i>=level;i--) {
  NavierStokes& ns_level=getLevel(i);
  ns_level.localMF[idx_localMF]->mult(dataval,0,ncomp,ngrow);
 }
}  // end subroutine mult_array

void NavierStokes::copyALL(int ngrow,int ncomp,
	int scomp,int dcomp,int idx_dest,int idx_source) {

 int finest_level = parent->finestLevel();
 if (level!=0)
  amrex::Error("level!=0 in copyALL");

 for (int i=finest_level;i>=level;i--) {
  NavierStokes& ns_level=getLevel(i);
  ns_level.Copy_localMF(idx_dest,idx_source,scomp,dcomp,ncomp,ngrow);
 }
} // end subroutine copyALL

void NavierStokes::minusALL(int ngrow,int ncomp,int idx_dest,int idx_source) {

 int finest_level = parent->finestLevel();
 if (level!=0)
  amrex::Error("level!=0 in minusALL");

 for (int i=finest_level;i>=level;i--) {
  NavierStokes& ns_level=getLevel(i);
  ns_level.minus_localMF(idx_dest,idx_source,ncomp,ngrow);
 }
} // end subroutine minusALL


void NavierStokes::delete_array(int idx_localMF) {

 int finest_level = parent->finestLevel();

 if (finest_level>=0) {

  if (level==0) {

   for (int i=level;i<=finest_level;i++) {
    NavierStokes& ns_level=getLevel(i);
    ns_level.delete_localMF(idx_localMF,1);
   }

  } else
   amrex::Error("level invalid");

 } else
  amrex::Error("finest_level invalid");

}  // subroutine delete_array

void NavierStokes::VOF_Recon_ALL(int ngrow,Real time,
  int update_flag,int init_vof_prev_time,
  int dest_mf) {

 if (level!=0)
  amrex::Error("level must be 0 ");
 if (verbose>0)
  if (ParallelDescriptor::IOProcessor()) 
   std::cout << "Start: VOF_Recon_ALL: time= " <<
    time << " update_flag= " << update_flag << 
    " init_vof_prev_time= " << init_vof_prev_time << '\n';

 if ((ngrow<1)||(ngrow>ngrow_distance))
  amrex::Error("ngrow invalid");

 int finest_level=parent->finestLevel();

  // go from coarsest to finest so that dest_mf
  // can have proper BC.
 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.VOF_Recon(ngrow,time,update_flag,
   init_vof_prev_time,dest_mf);
 }

 if (update_flag==0) {
  // do nothing
 } else if (update_flag==1) {
  avgDownError_ALL();
 } else
  amrex::Error("update_flag invalid");


} // subroutine VOF_Recon_ALL

void NavierStokes::VOF_Recon_resize(int ngrow,int dest_mf) {

 int nmat=num_materials;

 if ((ngrow<1)||(ngrow>ngrow_distance))
  amrex::Error("ngrow invalid");
 if (localMF[dest_mf]->nComp()!=nmat*ngeom_recon)
  amrex::Error("localMF[dest_mf]->nComp()!=nmat*ngeom_recon");

 if (localMF[dest_mf]->nGrow()==ngrow) {
  // do nothing
 } else if (localMF[dest_mf]->nGrow()>=0) {
  MultiFab* slopes_mf=new MultiFab(grids,dmap,nmat*ngeom_recon,0,
	  MFInfo().SetTag("slope_m"),FArrayBoxFactory());
  MultiFab::Copy(*slopes_mf,*localMF[dest_mf],0,0,nmat*ngeom_recon,0);
  delete_localMF(dest_mf,1);
  new_localMF(dest_mf,nmat*ngeom_recon,ngrow,-1); //sets values to 0.0
  MultiFab::Copy(*localMF[dest_mf],*slopes_mf,0,0,nmat*ngeom_recon,0);
 
  Vector<int> scompBC_map;
  scompBC_map.resize(nmat*ngeom_recon);
  for (int i=0;i<nmat*ngeom_recon;i++)
   scompBC_map[i]=i+1+AMREX_SPACEDIM;
  PCINTERP_fill_borders(dest_mf,ngrow,0,nmat*ngeom_recon,
    State_Type,scompBC_map);
 
  delete slopes_mf;
 } else
  amrex::Error("localMF[dest_mf]->nGrow() invalid");

} // subroutine VOF_Recon_resize

// vof,ref centroid,order,slope,intercept  x nmat
// update_flag=0 do not update the error
// update_flag=1 update S_new  (update error)
// 1. get MOF data with 1 ghost cell (so that CMOF can be chosen)
// 2. reconstruct interior cells only.
// 3. do extended filpatch; MOF used for coarse/fine and ext_dir cells.
void NavierStokes::VOF_Recon(int ngrow,Real time,
  int update_flag,int init_vof_prev_time,
  int dest_mf) {

 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();
 int max_level = parent->maxLevel();
 int nsteps=parent->levelSteps(0);

 if (ngrow_distance!=4)
  amrex::Error("ngrow_distance invalid");

 if ((ngrow<1)||(ngrow>ngrow_distance))
  amrex::Error("ngrow invalid");

 int bfact=parent->Space_blockingFactor(level);

 if ((verbose>0)&&(1==0)) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "beginning of VOF_Recon: ngrow,ngrow_distance,level: " << 
    ngrow << ' ' << ngrow_distance << ' ' << level << '\n';
   std::cout << "time,update_flag,init_vof_prev_time " <<
    time << ' ' << update_flag << ' ' << 
    init_vof_prev_time << '\n';
  } 
 }

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

  // 1000 used in linInterpFillFab
 Real teps=(cur_time_slab-prev_time_slab)/10000.0;  
 if ((time>prev_time_slab+teps)&&(time<cur_time_slab-teps))
  amrex::Error("cannot interp slope data in time");
 if ((time<prev_time_slab-teps)||(time>cur_time_slab+teps))
  amrex::Error("cannot extrapolate slope data in time");

 int nmat=num_materials;
 int nten=( (nmat-1)*(nmat-1)+nmat-1 )/2;

 int scomp_mofvars=(AMREX_SPACEDIM+1)+
  nmat*num_state_material;

 Vector< Vector<int> > total_calls;
 Vector< Vector<int> > total_iterations;
 total_calls.resize(thread_class::nthreads);
 total_iterations.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  total_calls[tid].resize(nmat);
  total_iterations[tid].resize(nmat);
  for (int im=0;im<nmat;im++) {
   total_calls[tid][im]=0;
   total_iterations[tid][im]=0;
  }
 } // tid

 double start_recon = ParallelDescriptor::second();

 MultiFab &S_new = get_new_data(State_Type,slab_step+1);

 MultiFab* lsdata=getStateDist(1,time,19);
 if (lsdata->nComp()!=nmat*(1+AMREX_SPACEDIM))
  amrex::Error("lsdata invalid ncomp");

 resize_mask_nbr(1);
 debug_ngrow(MASK_NBR_MF,1,90);
 if (localMF[MASK_NBR_MF]->nComp()!=4)
  amrex::Error("invalid ncomp for mask nbr");

 delete_localMF_if_exist(dest_mf,1);
 new_localMF(dest_mf,nmat*ngeom_recon,ngrow,-1);  // sets values to 0.0

 delete_localMF_if_exist(VOF_RECON_MF,1);
 getState_localMF(VOF_RECON_MF,1,scomp_mofvars,nmat*ngeom_raw,time);
 for (int im=0;im<nmat;im++) {
  int ibase_raw=im*ngeom_raw;
  int ibase_recon=im*ngeom_recon;
  Copy_localMF(dest_mf,VOF_RECON_MF,ibase_raw,ibase_recon,ngeom_raw,1);
  if (init_vof_prev_time==1) {
   int ngrow_save=localMF[VOF_PREV_TIME_MF]->nGrow();
   if (ngrow_save!=1)
    amrex::Error("vof prev time has invalid ngrow");
   Copy_localMF(VOF_PREV_TIME_MF,VOF_RECON_MF,ibase_raw,im,1,ngrow_save); 
  } else if (init_vof_prev_time==0) {
   // do nothing
  } else
   amrex::Error("init_vof_prev_time invalid");
 } // im=0..nmat-1

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,1,662);

 const Real* dx = geom.CellSize();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(localMF[dest_mf]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*localMF[dest_mf],use_tiling); mfi.isValid(); ++mfi) {
   BL_ASSERT(grids[mfi.index()] == mfi.validbox());
   const int gridno = mfi.index();
   const Box& tilegrid = mfi.tilebox();
   const Box& fabgrid = grids[gridno];
   const int* tilelo=tilegrid.loVect();
   const int* tilehi=tilegrid.hiVect();
   const int* fablo=fabgrid.loVect();
   const int* fabhi=fabgrid.hiVect();

   const Real* xlo = grid_loc[gridno].lo();

   FArrayBox& vfab=(*localMF[VOF_RECON_MF])[mfi];
   FArrayBox& lsfab=(*lsdata)[mfi];
   FArrayBox& mffab=(*localMF[dest_mf])[mfi];

   FArrayBox& snewfab=S_new[mfi];

   // mask_nbr:
   // (1) =1 interior  =1 fine-fine ghost in domain  =0 otherwise
   // (2) =1 interior  =0 otherwise
   // (3) =1 interior+ngrow-1  =0 otherwise
   // (4) =1 interior+ngrow    =0 otherwise
   FArrayBox& masknbr=(*localMF[MASK_NBR_MF])[mfi];

   Vector<int> vofbc=getBCArray(State_Type,gridno,scomp_mofvars,1);

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // in: PLIC_3D.F90
   fort_sloperecon(
    &tid_current,
    &gridno,
    &level,
    &finest_level,
    &max_level,
    &ngrow,
    vofbc.dataPtr(),
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    xlo,dx,
    masknbr.dataPtr(),
    ARLIM(masknbr.loVect()),ARLIM(masknbr.hiVect()),
    snewfab.dataPtr(scomp_mofvars),
    ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()),
    vfab.dataPtr(),ARLIM(vfab.loVect()),ARLIM(vfab.hiVect()),
    lsfab.dataPtr(),ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
    mffab.dataPtr(),ARLIM(mffab.loVect()),ARLIM(mffab.hiVect()),
    &nsteps,
    &time,
    &nmat,&nten,
    latent_heat.dataPtr(),
    &update_flag,
    total_calls[tid_current].dataPtr(),
    total_iterations[tid_current].dataPtr(),
    &continuous_mof, 
    &force_cmof_at_triple_junctions, 
    &partial_cmof_stencil_at_walls, 
    radius_cutoff.dataPtr());
 }  // mfi
} // omp
 ns_reconcile_d_num(158);

 for (int tid=1;tid<thread_class::nthreads;tid++) {
  for (int im=0;im<nmat;im++) {
   total_calls[0][im]+=total_calls[tid][im];
   total_iterations[0][im]+=total_iterations[tid][im];
  }
 } // tid

 for (int im=0;im<nmat;im++) {
   ParallelDescriptor::ReduceIntSum(total_calls[0][im]);
   ParallelDescriptor::ReduceIntSum(total_iterations[0][im]);
 }

 Vector<int> scompBC_map;
 scompBC_map.resize(nmat*ngeom_recon);
 for (int i=0;i<nmat*ngeom_recon;i++)
  scompBC_map[i]=i+1+AMREX_SPACEDIM;

 if (1==0) {
  for (int im=0;im<nmat;im++) {
   int vofcomp=im*ngeom_recon;
   std::cout << "lev,im,calls,iter " << level << ' ' << 
    im << ' ' << total_calls[0][im] << ' ' <<
    total_iterations[0][im] << '\n';
   Real vfracnrm=localMF[dest_mf]->norm1(vofcomp,0);
   std::cout << "im= " << im << " vfracnrm= " << vfracnrm << '\n';
  }
 }

 PCINTERP_fill_borders(dest_mf,ngrow,0,nmat*ngeom_recon,
   State_Type,scompBC_map);

 double end_recon = ParallelDescriptor::second();
 double cputime=end_recon-start_recon;

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   for (int im=0;im<nmat;im++) {
    Real itercall=0.0;
    Real r_calls=total_calls[0][im];
    Real r_iters=total_iterations[0][im];
    if (total_calls[0][im]>0) 
     itercall=r_iters/r_calls;
    std::cout << "lev,im,calls,iter,iter/call,cpu time " << level << ' ' << 
     im << ' ' << total_calls[0][im] << ' ' <<
     total_iterations[0][im] << ' ' << itercall << ' ' << cputime << '\n';
   } // im=0..nmat-1
  }
 } else if (verbose<=0) {
  // do nothing
 } else
  amrex::Error("verbose invalid");

 delete lsdata;

}  // subroutine VOF_Recon

// called from: 
//  NavierStokes::sum_integrated_quantities
//  NavierStokes::prepare_post_process
//  NavierStokes::advance
void NavierStokes::build_masksemALL() {

 if (level!=0)
  amrex::Error("level invalid in build_masksemALL");

  // mask_sweep==0: init mask within elements
  // mask_sweep==1: update mask within elements by using
  //                the mask from neighboring elements.
 int finest_level = parent->finestLevel();
 for (int mask_sweep=0;mask_sweep<2;mask_sweep++) {

   // traverse from coarsest to finest so that
   // coarse data mask will be available for filling in
   // ghost values.
  for (int ilev=level;ilev<=finest_level;ilev++) {
   NavierStokes& ns_level=getLevel(ilev);
   ns_level.build_masksem(mask_sweep);
  } // ilev

 } //mask_sweep

} // subroutine build_masksemALL

void NavierStokes::build_masksem(int mask_sweep) {

 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();
 int bfact=parent->Space_blockingFactor(level);
 int bfact_fine=bfact;
 if (level<finest_level) 
  bfact_fine=parent->Space_blockingFactor(level+1);
 if (bfact_fine>bfact)
  amrex::Error("bfact_fine invalid");

 int nmat=num_materials;

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,1,80);
 resize_mask_nbr(1);
 debug_ngrow(MASK_NBR_MF,1,90);
 if (localMF[MASK_NBR_MF]->nComp()!=4)
  amrex::Error("invalid ncomp for mask nbr");

 MultiFab* old_mask;

 if (mask_sweep==0) { // init mask within elements w/o nbr. info.

  if (localMF_grow[MASKSEM_MF]==1) {

   // do nothing

  } else if (localMF_grow[MASKSEM_MF]==-1) {

   // ncomp=1 ngrow=1
   new_localMF(MASKSEM_MF,1,1,-1);

  } else
   amrex::Error("localMF_grow invalid");

   // ncomp=1 ngrow=1
  localMF[MASKSEM_MF]->setVal(0.0,0,1,1);
  old_mask=localMF[MASKSEM_MF];

 } else if (mask_sweep==1) { // update mask using nbr. info.

  debug_ngrow(MASKSEM_MF,1,80);
  old_mask=new MultiFab(grids,dmap,1,1,
	MFInfo().SetTag("old_mask"),FArrayBoxFactory());
  MultiFab::Copy(*old_mask,*localMF[MASKSEM_MF],0,0,1,1);

 } else
  amrex::Error("mask_sweep invalid");

 Real total_cells_level=localMF[MASKCOEF_MF]->norm1();

 Vector< Vector< Real >  > spectral_cells_level;
 spectral_cells_level.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  spectral_cells_level[tid].resize(nmat);
  for (int im=0;im<nmat;im++)
   spectral_cells_level[tid][im]=0.0;
 }

 int scomp_mofvars=(AMREX_SPACEDIM+1)+
   nmat*num_state_material;

 MultiFab* vofmat=new MultiFab(grids,dmap,nmat,0,
  MFInfo().SetTag("vofmat"),FArrayBoxFactory());
 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

 for (int im=0;im<nmat;im++) {
  int scomp=scomp_mofvars+im*ngeom_raw;
  MultiFab::Copy(*vofmat,S_new,scomp,im,1,0);
 }
 const Box& domain = geom.Domain();
 const int* domlo = domain.loVect();
 const int* domhi = domain.hiVect();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(vofmat->boxArray().d_numPts());


#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*vofmat,use_tiling); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];
  const int* tilelo=tilegrid.loVect();
  const int* tilehi=tilegrid.hiVect();
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();

   // mask=tag if not covered by level+1 or outside the domain.
  FArrayBox& maskcov=(*localMF[MASKCOEF_MF])[mfi];
   // (1) =1 interior  =1 fine-fine ghost in domain  =0 otherwise
   // (2) =1 interior  =0 otherwise
  FArrayBox& masknbr=(*localMF[MASK_NBR_MF])[mfi];

  FArrayBox& maskfab=(*localMF[MASKSEM_MF])[mfi];
  FArrayBox& oldfab=(*old_mask)[mfi];
  FArrayBox& voffab=(*vofmat)[mfi];

  Vector<int> vofbc=getBCArray(State_Type,gridno,scomp_mofvars,1);

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // in: GODUNOV_3D.F90
  FORT_BUILD_MASKSEM( 
   spectral_cells_level[tid_current].dataPtr(),
   &mask_sweep,
   &level,
   &finest_level,
   &cur_time_slab,
   &enable_spectral,
   domlo,domhi,
   vofbc.dataPtr(),
   maskcov.dataPtr(),ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
   masknbr.dataPtr(),ARLIM(masknbr.loVect()),ARLIM(masknbr.hiVect()),
   maskfab.dataPtr(),ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
   oldfab.dataPtr(),ARLIM(oldfab.loVect()),ARLIM(oldfab.hiVect()),
   voffab.dataPtr(),ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),
   tilelo,tilehi,
   fablo,fabhi,
   &bfact,
   &bfact_fine,
   &nmat);
 }  // mfi
} // omp
 ns_reconcile_d_num(159);

 for (int tid=1;tid<thread_class::nthreads;tid++) {
  for (int im=0;im<nmat;im++)
   spectral_cells_level[0][im]+=spectral_cells_level[tid][im];
 }
 for (int im=0;im<nmat;im++)
  ParallelDescriptor::ReduceRealSum(spectral_cells_level[0][im]);

 Vector<int> scompBC_map;
 scompBC_map.resize(1);
 scompBC_map[0]=1+AMREX_SPACEDIM+nmat*ngeom_recon;
  // std::string maskextrap_str="maskSEMextrap"
  // FORT_EXTRAPFILL
  // mask_sem_interp (MASKINTERPPC)
  // idx,ngrow,scomp,ncomp,index,scompBC_map
 PCINTERP_fill_borders(MASKSEM_MF,1,0,1,State_Type,scompBC_map);

 delete vofmat;
  
 if (mask_sweep==1)
  delete old_mask;

 std::fflush(NULL);

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   if (mask_sweep==1) {
    std::cout << "level,finest_level " <<
     level << ' ' << finest_level << '\n';
    std::cout << "total_cells_lev " << total_cells_level << '\n';
    Real spectral_cells_sum=0.0;
    for (int im=0;im<nmat;im++) {
     std::cout << "im, spectral_cells_lev " << im << ' ' <<
      spectral_cells_level[0][im] << '\n';
     spectral_cells_sum+=spectral_cells_level[0][im];
    }
    if (total_cells_level>0.0) {
     Real frac=spectral_cells_sum/total_cells_level;
     std::cout << "frac " << frac << '\n';
    }
   } // mask_sweep==1
  } // IOProc
 } // verbose>0

 std::fflush(NULL);

 if (1==0) {
  int gridno=0;
  const Box& fabgrid = grids[gridno];
  const int* fablo=fabgrid.loVect();
  const int* fabhi=fabgrid.hiVect();
  const Real* xlo = grid_loc[gridno].lo();
  int interior_only=1;
  MultiFab& mf_out=*localMF[MASKSEM_MF];
  FArrayBox& mffab=mf_out[0];
  const Real* dx = geom.CellSize();
  int scomp=0;
  int ncomp=1;
  std::cout << " checking out masksem_mf \n";
  tecplot_debug(mffab,xlo,fablo,fabhi,dx,-1,0,scomp,ncomp,interior_only);
 }

} // subroutine build_masksem()



// If incompressible material then copy existing state pressure.
// If compressible material then P=P(rho,e).
MultiFab* NavierStokes::derive_EOS_pressure(Vector<int> local_material_type) {

 int finest_level=parent->finestLevel();
 
 bool use_tiling=ns_tiling;

 int nmat=num_materials;

 debug_ngrow(LEVELPC_MF,1,662);
 if (localMF[LEVELPC_MF]->nComp()!=nmat*(AMREX_SPACEDIM+1)) {
  std::cout << "in: derive_EOS_pressure\n";
  amrex::Error("localMF[LEVELPC_MF]->nComp()!=nmat*(AMREX_SPACEDIM+1)");
 }
 resize_levelsetLO(2,LEVELPC_MF);
 debug_ngrow(LEVELPC_MF,2,662);
 if (localMF[LEVELPC_MF]->nComp()!=nmat*(AMREX_SPACEDIM+1))
  amrex::Error("localMF[LEVELPC_MF]->nComp()!=nmat*(AMREX_SPACEDIM+1)");

 VOF_Recon_resize(1,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,1,663);

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 MultiFab* mf=getStatePres(1,cur_time_slab); 

 if (mf->nComp()!=1)
  amrex::Error("mf->nComp() invalid");

 MultiFab* denmf=getStateDen(1,cur_time_slab);  // nmat * den,temp,...
 int nden=denmf->nComp();
 
 const Real* dx = geom.CellSize();
 
 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(denmf->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*denmf,use_tiling); mfi.isValid(); ++mfi) {
 
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
  FArrayBox& levelpcfab=(*localMF[LEVELPC_MF])[mfi];
  FArrayBox& denfab=(*denmf)[mfi];
  FArrayBox& presfab=(*mf)[mfi];

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();
  
   // declared in: NAVIERSTOKES_3D.F90 
  FORT_EOS_PRESSURE(
   &level,
   &finest_level,
   local_material_type.dataPtr(),
   xlo,dx,
   presfab.dataPtr(),
   ARLIM(presfab.loVect()),ARLIM(presfab.hiVect()),
   voffab.dataPtr(),
   ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),
   levelpcfab.dataPtr(),
   ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),
   denfab.dataPtr(),ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
   tilelo,tilehi,
   fablo,fabhi,
   &bfact,
   &nmat,
   &nden);
 
 } // mfi
} // omp
 ns_reconcile_d_num(160);
 
 delete denmf;

 return mf;
} // subroutine derive_EOS_pressure()

// in NavierStokes::multiphase_project when:
// project_option==0 
//  and homflag=0,
//  the following commands are given:
//  for ilev=finest ... coarsest,
//   ns_level.init_pressure_error_indicator();
//  avgDownError_ALL();
// during regridding, the following routine checks the error:
//  NavierStokes::errorEst  (calls VFRACERROR)
void NavierStokes::init_pressure_error_indicator() {

 bool use_tiling=ns_tiling;

 int finest_level = parent->finestLevel();

 int nmat=num_materials;
 int ntensor=AMREX_SPACEDIM*AMREX_SPACEDIM;

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,1,6001);

 debug_ngrow(CELLTENSOR_MF,1,9);
 if (localMF[CELLTENSOR_MF]->nComp()!=ntensor)
  amrex::Error("localMF[CELLTENSOR_MF]->nComp() invalid");

 VOF_Recon_resize(1,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,1,664);
 if (localMF[SLOPE_RECON_MF]->nComp()!=nmat*ngeom_recon)
  amrex::Error("localMF[SLOPERECON_MF]->nComp() invalid");

 resize_levelsetLO(2,LEVELPC_MF);
 debug_ngrow(LEVELPC_MF,2,665);

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 MultiFab* denmf=getStateDen(1,cur_time_slab); 
 if (denmf->nComp()!=nmat*num_state_material)
  amrex::Error("denmf incorrect ncomp");

 MultiFab* presmf=getState(1,AMREX_SPACEDIM,1,cur_time_slab);

 MultiFab* vortmf=new MultiFab(grids,dmap,1,0,
	MFInfo().SetTag("vortmf"),FArrayBoxFactory());
 const Real* dx = geom.CellSize();
 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 int scomp_mofvars=(AMREX_SPACEDIM+1)+
   nmat*num_state_material;
 int scomp_error=scomp_mofvars+nmat*ngeom_raw;
 if (scomp_error!=S_new.nComp()-1)
  amrex::Error("scomp_error invalid");

 MultiFab* velmf=getState(1,0,AMREX_SPACEDIM,cur_time_slab);

 if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(velmf->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*velmf,use_tiling); mfi.isValid(); ++mfi) {
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
   Vector<int> velbc=
     getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);
   FArrayBox& velfab=(*velmf)[mfi];
   FArrayBox& vortfab=(*vortmf)[mfi];

   FArrayBox& cellten=(*localMF[CELLTENSOR_MF])[mfi];

   if (cellten.nComp()!=ntensor)
    amrex::Error("cellten invalid ncomp");

   int tencomp=0;
   int iproject=0;
   int onlyscalar=2; // magnitude of vorticity
   int ngrow_zero=0;

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   fort_getshear(
    &ntensor,
    cellten.dataPtr(tencomp),ARLIM(cellten.loVect()),ARLIM(cellten.hiVect()),
    velfab.dataPtr(),ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
    dx,xlo,
    vortfab.dataPtr(),
    ARLIM(vortfab.loVect()),ARLIM(vortfab.hiVect()),
    &iproject,
    &onlyscalar,
    &cur_time_slab,
    tilelo,tilehi,
    fablo,fabhi,
    &bfact,
    &level,
    velbc.dataPtr(),
    &ngrow_zero,
    &nmat);
 } // mfi
} // omp
 ns_reconcile_d_num(161);

 delete velmf;

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(denmf->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*denmf,use_tiling); mfi.isValid(); ++mfi) {

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

  FArrayBox& voffab=(*localMF[SLOPE_RECON_MF])[mfi];
  FArrayBox& denfab=(*denmf)[mfi];
  FArrayBox& presfab=(*presmf)[mfi];
  FArrayBox& vortfab=(*vortmf)[mfi];
  FArrayBox& errnew=S_new[mfi];

  if (denfab.nComp()!=nmat*num_state_material)
   amrex::Error("denfab.nComp() invalid");
  if (voffab.nComp()!=nmat*ngeom_recon)
   amrex::Error("voffab.nComp() invalid");
  if (presfab.nComp()!=1)
   amrex::Error("presfab.nComp() invalid");
  if (vortfab.nComp()!=1)
   amrex::Error("vortfab.nComp() invalid");

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  FORT_PRESSURE_INDICATOR(
   &pressure_error_flag,
   vorterr.dataPtr(),
   pressure_error_cutoff.dataPtr(),
   temperature_error_cutoff.dataPtr(),
   xlo,dx,
   errnew.dataPtr(scomp_error),
   ARLIM(errnew.loVect()),ARLIM(errnew.hiVect()),
   voffab.dataPtr(),ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),
   denfab.dataPtr(),ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
   vortfab.dataPtr(),ARLIM(vortfab.loVect()),ARLIM(vortfab.hiVect()),
   presfab.dataPtr(),ARLIM(presfab.loVect()),ARLIM(presfab.hiVect()),
   maskcov.dataPtr(),ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
   tilelo,tilehi,
   fablo,fabhi, 
   &bfact,
   &level,
   &finest_level,
   &nmat);

 } // mfi
} // omp
 ns_reconcile_d_num(162);

 delete presmf;
 delete vortmf;
 delete denmf;

} // subroutine init_pressure_error_indicator

// if project_option==0:
// 1. calculates p(rho^n+1,e_advect) and puts it in 2nd component
//    of cell_sound.
//    Equation of state to be used depends on vofPC
// 2. calculates 1/(rho c^2 dt^2)  and puts it in 1st component of cell_sound.
//    Equation of state to be used depends on vofPC
// 3. in incompressible regions, p=0
//
// div_hold=(pnew-pold)/(rho c^2 dt) + dt mdot/vol
// if project_option==11
// 1. div_hold*vol/dt is put in localMF[DIFFUSIONRHS_MF] in incomp parts
// 2. div_hold/(csound_hold*dt) is put in the 2nd component of cell_sound where
//    there are compressible materials.
//
// velocity scale: V
// time scale is : 1/V
// pressure scale: V^2
// scale for "cell_sound" is 1
//
void NavierStokes::init_advective_pressure(int project_option) {
 
 int finest_level=parent->finestLevel();
 if ((level>=0)&&(level<=finest_level)) {
  // do nothing
 } else
  amrex::Error("level invalid");

 bool use_tiling=ns_tiling;

 resize_metrics(1);
 debug_ngrow(VOLUME_MF,0,700);

 debug_ngrow(FACE_VAR_MF,0,660);

 VOF_Recon_resize(1,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,1,661);

 debug_ngrow(MASKCOEF_MF,1,6001);

 debug_ngrow(DIFFUSIONRHS_MF,0,661);

 if (localMF[DIFFUSIONRHS_MF]->nComp()!=1)
  amrex::Error("localMF[DIFFUSIONRHS_MF]->nComp() invalid");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int nmat=num_materials;

 for (int im=0;im<nmat;im++) {
  if ((compressible_dt_factor[im]>=1.0)&&
      (compressible_dt_factor[im]<=1.0e+20)) {
   // do nothing
  } else
   amrex::Error("compressible_dt_factor[im] invalid");
 }

 MultiFab* denmf=getStateDen(1,cur_time_slab);  // nmat x den,temp, ...
 int nden=denmf->nComp();

 const Real* dx = geom.CellSize();

 int nsolve=1;

 Vector<int> scomp;
 Vector<int> ncomp;
 int ncomp_check;
 int state_index;
  //num_materials_combine=1
 get_mm_scomp_solver(
  1,
  project_option,
  state_index,
  scomp,ncomp,ncomp_check);

 if (ncomp_check!=nsolve)
  amrex::Error("ncomp_check invalid");

 if (scomp.size()!=1)
  amrex::Error("scomp.size() invalid");
 if (ncomp[0]!=1)
  amrex::Error("ncomp[0] invalid");
  
 if (project_option==0) {
  if (state_index!=State_Type)
   amrex::Error("state_index invalid");
 } else if (project_option==11) { //FSI_material_exists last project
  if (state_index!=DIV_Type)
   amrex::Error("state_index invalid");
 } else
  amrex::Error("project_option invalid28");

 MultiFab& S_new=get_new_data(state_index,slab_step+1);
 MultiFab& LS_new=get_new_data(LS_Type,slab_step+1);

  // CELL_SOUND_MF
  // coeff_avg,padvect_avg 
 if (project_option==0) {
  // do nothing
 } else if (project_option==11) { //FSI_material_exists last project
   // dst,src,scomp,dcomp,ncomp,ngrow
  int sc=scomp[0];
  int dc=1; // copy 1st component of DIV_TYPE contents 
            // into 2nd component of 
            // localMF[CELL_SOUND_MF]
            // DIV_Type=-(pnew-pold)/(rho c^2 dt) + dt mdot/vol
  MultiFab::Copy(*localMF[CELL_SOUND_MF],S_new,sc,dc,1,0);
 } else
  amrex::Error("project_option invalid29");
  
 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(denmf->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*denmf,use_tiling); mfi.isValid(); ++mfi) {

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

  FArrayBox& denfab=(*denmf)[mfi];

  FArrayBox& volumefab=(*localMF[VOLUME_MF])[mfi];

   // tessellating volume fractions.
  FArrayBox& voffab=(*localMF[CELL_VOF_MF])[mfi];
  FArrayBox& csoundfab=(*localMF[CELL_SOUND_MF])[mfi];
  FArrayBox& mdotfab=(*localMF[DIFFUSIONRHS_MF])[mfi];
  FArrayBox& lsnewfab=LS_new[mfi];
 
  if (lsnewfab.nComp()!=nmat*(AMREX_SPACEDIM+1))
   amrex::Error("lsnewfab.nComp()!=nmat*(AMREX_SPACEDIM+1)");
 
  if (mdotfab.nComp()!=1)
   amrex::Error("mdotfab.nComp() invalid");
  if (csoundfab.nComp()!=2)
   amrex::Error("csoundfab.nComp() invalid");

  // mask=tag if not covered by level+1 or outside the domain.
  FArrayBox& maskcov=(*localMF[MASKCOEF_MF])[mfi];

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // in: NAVIERSTOKES_3D.F90
  FORT_ADVECTIVE_PRESSURE(
   &level,
   &finest_level,
   xlo,dx,
   &dt_slab,
   maskcov.dataPtr(),
   ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
   volumefab.dataPtr(),
   ARLIM(volumefab.loVect()),ARLIM(volumefab.hiVect()),
   lsnewfab.dataPtr(),
   ARLIM(lsnewfab.loVect()),ARLIM(lsnewfab.hiVect()),
   csoundfab.dataPtr(),
   ARLIM(csoundfab.loVect()),ARLIM(csoundfab.hiVect()),
   voffab.dataPtr(),ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),
   denfab.dataPtr(),
   ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
   mdotfab.dataPtr(),ARLIM(mdotfab.loVect()),ARLIM(mdotfab.hiVect()),
   tilelo,tilehi,
   fablo,fabhi,&bfact,
   &nmat,&nden,
   compressible_dt_factor.dataPtr(),
   &pressure_select_criterion, 
   &project_option);

 } // mfi
} // omp
 ns_reconcile_d_num(163);

  // CELL_SOUND_MF
  // coeff_avg,padvect_avg
  // dst,src,scomp,dcomp,ncomp,ngrow
 int sc=1; // padvect_avg
 int dc=scomp[0];
 MultiFab::Copy(S_new,*localMF[CELL_SOUND_MF],sc,dc,1,0);

 delete denmf;

} // subroutine init_advective_pressure


MultiFab* NavierStokes::getStateDen(int ngrow,Real time) { 

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int nmat=num_materials;
 int ncomp=nmat*num_state_material;
 int scomp=(AMREX_SPACEDIM+1);
 MultiFab* mf=getState(ngrow,scomp,ncomp,time);
 return mf;

}  // subroutine getStateDen


MultiFab* NavierStokes::getStatePres(int ngrow,Real time) {

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int scomp_pres=AMREX_SPACEDIM;
 MultiFab* mf=getState(ngrow,scomp_pres,1,time);
 return mf;

}  // subroutine getStatePres

void NavierStokes::scale_variablesALL() {

 if (level!=0)
  amrex::Error("level invalid scale_variablesALL");

 int finest_level = parent->finestLevel();
  // in: PROB.F90
 FORT_SETFORTSCALES(&projection_pressure_scale,
   &projection_velocity_scale);

 dt_slab*=projection_velocity_scale;

 int scale_flag=0;
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.scale_variables(scale_flag);
 }

} // end subroutine NavierStokes::scale_variablesALL()


void NavierStokes::unscale_variablesALL() {

 if (level!=0)
  amrex::Error("level invalid unscale_variablesALL");

 int finest_level = parent->finestLevel();
 Real dummy_scale=1.0;
 FORT_SETFORTSCALES(&dummy_scale,&dummy_scale);

 dt_slab/=projection_velocity_scale;

 int scale_flag=1;
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.scale_variables(scale_flag);
 }

} // end subroutine NavierStokes::unscale_variablesALL()

//scale Unew,Umac_new,P,mdot,solid_vars,even components of CELL_SOUND (padvect)
void NavierStokes::scale_variables(int scale_flag) {

 int finest_level = parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid scale_variables");

 int nmat=num_materials;

 int nparts=im_solid_map.size();
 if ((nparts<0)||(nparts>nmat))
  amrex::Error("nparts invalid");

 int nparts_def=nparts;
 if (nparts==0) {
  nparts_def=1;
 } else if ((nparts>=1)&&(nparts<=nmat)) {
  // do nothing
 } else
  amrex::Error("nparts invalid");

 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
  if (localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0)
   amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0");
  if (localMF[FSI_GHOST_MAC_MF+data_dir]->nComp()!=nparts_def*AMREX_SPACEDIM)
   amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nComp() bad");
 }
 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
  debug_ngrow(FSI_GHOST_MAC_MF+data_dir,0,112);
 }

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  debug_ngrow(FACE_VAR_MF+dir,0,660);
 }

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");
 
 Real pres_factor=1.0/projection_pressure_scale;
 Real vel_factor=1.0/projection_velocity_scale;
 if (scale_flag==0) {
  pres_factor=1.0/projection_pressure_scale;
  vel_factor=1.0/projection_velocity_scale;
 } else if (scale_flag==1) {
  pres_factor=projection_pressure_scale;
  vel_factor=projection_velocity_scale;
 } else
  amrex::Error("scale_flag invalid");

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 MultiFab& DIV_new=get_new_data(DIV_Type,slab_step+1);

 int nsolve=AMREX_SPACEDIM;

 S_new.mult(vel_factor,0,nsolve,0);

 nsolve=1;

 int scomp_pres=AMREX_SPACEDIM;
 S_new.mult(pres_factor,scomp_pres,nsolve,0);

  // DIV_new contains -dt (pnew-padv)/(rho c^2 dt^2) + MDOT_MF dt/vol
 DIV_new.mult(vel_factor,0,nsolve,0);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  MultiFab& Umac_new=get_new_data(Umac_Type+dir,slab_step+1);
  MultiFab& Umac_old=get_new_data(Umac_Type+dir,slab_step);
  int ncmac=Umac_new.nComp();
  if (ncmac!=nsolve) {
   std::cout << "nmat = " << nmat << '\n';
   std::cout << "ncmac = " << ncmac << '\n';
   amrex::Error("ncmac invalid scale_variables");
  }
  Umac_new.mult(vel_factor,0,nsolve,0);
  Umac_old.mult(vel_factor,0,nsolve,0);
  localMF[FACE_VAR_MF+dir]->mult(vel_factor,facevel_index,1,0);
  localMF[FSI_GHOST_MAC_MF+dir]->mult(vel_factor,0,nparts_def*AMREX_SPACEDIM,0);
 } // dir

 localMF[DIFFUSIONRHS_MF]->mult(pres_factor,0,nsolve,0);

 // coeff_avg,padvect_avg 
 localMF[CELL_SOUND_MF]->mult(pres_factor,1,1,0);

 if ((nparts>=1)&&(nparts<=nmat)) {
  MultiFab& Solid_new=get_new_data(Solid_State_Type,slab_step+1);
  Solid_new.mult(vel_factor,0,nparts*AMREX_SPACEDIM,0);
 } else if (nparts==0) {
  // do nothing
 } else
  amrex::Error("nparts invalid");

}  // subroutine scale_variables


// 1. viscosity coefficient - 1..nmat
// 2. viscoelastic coefficient - 1..nmat
// 3. relaxation time - 1..nmat
// the viscous and viscoelastic forces should both be multiplied by
// visc_coef.  
void NavierStokes::getStateVISC(int idx,int ngrow) {

 delete_localMF_if_exist(idx,1);

 if ((ngrow==0)||(ngrow==1)) {
  // do nothing
 } else
  amrex::Error("ngrow invalid in getStateVISC");
 
 int finest_level=parent->finestLevel();

 bool use_tiling=ns_tiling;

 int nmat=num_materials;
 int ntensor=AMREX_SPACEDIM*AMREX_SPACEDIM;

  // init_gradu_tensorALL is called in NavierStokes::make_physics_varsALL
  // prior to this routine being called.
  // Also, init_gradu_tensorALL is called in NavierStokes::writeTECPLOT_File
  // prior to this routine being called.
 debug_ngrow(CELLTENSOR_MF,1,9);
 if (localMF[CELLTENSOR_MF]->nComp()!=ntensor)
  amrex::Error("localMF[CELLTENSOR_MF]->nComp() invalid");

 VOF_Recon_resize(2,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,2,670);
 if (localMF[SLOPE_RECON_MF]->nComp()!=nmat*ngeom_recon)
  amrex::Error("localMF[SLOPE_RECON_MF]->nComp() invalid");

 int ncomp_visc=3*nmat;

 Vector<int> shear_thinning_fluid(nmat);

 for (int im=0;im<nmat;im++) {

  shear_thinning_fluid[im]=0;

  if (ns_is_rigid(im)==1) {
   // do nothing
  } else if (ns_is_rigid(im)==0) {

   if ((probtype==2)&&(axis_dir>0)&&(im==0))
    shear_thinning_fluid[im]=1;
   if (Carreau_beta[im]!=0.0)
    shear_thinning_fluid[im]=1;

   if (shear_thinning_fluid[im]==1) {
    // do nothing
   } else if (shear_thinning_fluid[im]==0) {
    // do nothing
   } else
    amrex::Error("shear_thinning_fluid invalid");

   if (viscosity_state_model[im]==0) {
    // do nothing
   } else if (viscosity_state_model[im]>0) { // viscosity depends on T
    // do nothing
   } else
    amrex::Error("viscosity_state_model invalid");

   if (les_model[im]==0) {
    // do nothing
   } else if (les_model[im]==1) {
    // do nothing
   } else
    amrex::Error("les_model invalid");

   if ((elastic_time[im]>0.0)&&
       (elastic_viscosity[im]>0.0)) {
    // do nothing
   } else if ((elastic_time[im]==0.0)||
              (elastic_viscosity[im]==0.0)) {
    if (viscoelastic_model[im]!=0)
     amrex::Error("viscoelastic_model[im]!=0");
   } else
    amrex::Error("elastic_time/elastic_viscosity getStateVISC");

  } else
   amrex::Error("ns_is_rigid invalid");

 } // im=0..nmat-1

 new_localMF(idx,ncomp_visc,ngrow,-1); // sets values to 0.0

 MultiFab* vel=getState(ngrow+1,0,AMREX_SPACEDIM,cur_time_slab);

 MultiFab* EOSdata=getStateDen(ngrow,cur_time_slab);

 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=nmat)) {
  // do nothing
 } else 
  amrex::Error("num_materials_viscoelastic invalid");

 MultiFab* tensor=getStateTensor(ngrow,0,
     num_materials_viscoelastic*NUM_TENSOR_TYPE,cur_time_slab);

 for (int im=0;im<nmat;im++) {

  const Real* dx = geom.CellSize();

  MultiFab* gammadot_mf=new MultiFab(grids,dmap,1,ngrow,
	MFInfo().SetTag("gammadot_mf"),FArrayBoxFactory());
  gammadot_mf->setVal(0.0);

  int scomp_tensor=0;

  if (store_elastic_data[im]==1) {

   int partid=0;
   while ((im_elastic_map[partid]!=im)&&
          (partid<im_elastic_map.size())) {
    partid++;
   }
   if (partid<im_elastic_map.size()) {
    scomp_tensor=partid*NUM_TENSOR_TYPE;
   } else
    amrex::Error("partid could not be found: getStateVISC");
     
  } else if (store_elastic_data[im]==0) {

   if (viscoelastic_model[im]!=0)
    amrex::Error("viscoelastic_model[im]!=0");

  } else
   amrex::Error("store_elastic_data invalid getStateVISC");

  if (ns_is_rigid(im)==1) {
   // no need to call getshear
  } else if (ns_is_rigid(im)==0) {

   if (shear_thinning_fluid[im]==1) {
    int iproject=0;
    int onlyscalar=1;  // mag(trace gradu) 

    if (thread_class::nthreads<1)
     amrex::Error("thread_class::nthreads invalid");
    thread_class::init_d_numPts(vel->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
    for (MFIter mfi(*vel,use_tiling); mfi.isValid(); ++mfi) {
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

     FArrayBox& gammadot=(*gammadot_mf)[mfi];

     FArrayBox& cellten=(*localMF[CELLTENSOR_MF])[mfi];
     if (cellten.nComp()!=ntensor)
      amrex::Error("cellten invalid ncomp");

     Vector<int> velbc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);

     int tid_current=ns_thread();
     if ((tid_current<0)||(tid_current>=thread_class::nthreads))
      amrex::Error("tid_current invalid");
     thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

      // since onlyscalar==1, this routine calculates mag(trace gradu)=
      //  sqrt(2 D:D)
      // since D is symmetric, D:D=trace(D^2) 
      // is invariant to coordinate transformations.
      // if levelrz==1, gradu(3,3)=u/|r|
     fort_getshear(
      &ntensor,
      cellten.dataPtr(),
      ARLIM(cellten.loVect()),ARLIM(cellten.hiVect()),
      (*vel)[mfi].dataPtr(),
      ARLIM((*vel)[mfi].loVect()),ARLIM((*vel)[mfi].hiVect()),
      dx,xlo,
      gammadot.dataPtr(),
      ARLIM(gammadot.loVect()),ARLIM(gammadot.hiVect()),
      &iproject,&onlyscalar,
      &cur_time_slab,
      tilelo,tilehi,
      fablo,fabhi,
      &bfact,
      &level,
      velbc.dataPtr(),
      &ngrow,&nmat);

    } //mfi
} // omp
    ns_reconcile_d_num(164);

   } else if (shear_thinning_fluid[im]==0) {
     // do nothing
   } else 
    amrex::Error("shear_thinning_fluid invalid");
  } else 
   amrex::Error("ns_is_rigid invalid");

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(vel->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(*vel,use_tiling); mfi.isValid(); ++mfi) {
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

   FArrayBox& gammadot=(*gammadot_mf)[mfi];

   FArrayBox& viscfab=(*localMF[idx])[mfi];

   FArrayBox& velfab=(*vel)[mfi];
   FArrayBox& eosfab=(*EOSdata)[mfi];
   FArrayBox& tensorfab=(*tensor)[mfi];

   Vector<int> velbc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);

   int fortran_im=im+1;

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // declared in: DERIVE_3D.F90
   fort_derviscosity(
      &level,
      &finest_level,
      &visc_coef,
      &fortran_im,
      &nmat,
      &dt_slab,
      &viscconst[im],
      &shear_thinning_fluid[im],
      &Carreau_alpha[im],
      &Carreau_beta[im],
      &Carreau_n[im],
      &Carreau_mu_inf[im],
      &concentration[im],
      &elastic_time[im],
      &viscosity_state_model[im],
      &viscoelastic_model[im],
      &elastic_viscosity[im],
      &elastic_regularization[im],
      &etaL[im],&etaP[im],&etaS[im],
      &polymer_factor[im],
      viscfab.dataPtr(),
      ARLIM(viscfab.loVect()),ARLIM(viscfab.hiVect()),
      velfab.dataPtr(),
      ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
      eosfab.dataPtr(),
      ARLIM(eosfab.loVect()),ARLIM(eosfab.hiVect()),
      tensorfab.dataPtr(scomp_tensor),
      ARLIM(tensorfab.loVect()),ARLIM(tensorfab.hiVect()),
      gammadot.dataPtr(),
      ARLIM(gammadot.loVect()),ARLIM(gammadot.hiVect()),
      tilelo,tilehi,
      fablo,fabhi,&bfact,
      &cur_time_slab,
      dx,xlo,
      velbc.dataPtr(),&ngrow,
      &ncomp_visc);
  } //mfi
} // omp
  ns_reconcile_d_num(165);

  if (les_model[im]==1) {

   if (thread_class::nthreads<1)
    amrex::Error("thread_class::nthreads invalid");
   thread_class::init_d_numPts(vel->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
   for (MFIter mfi(*vel,use_tiling); mfi.isValid(); ++mfi) {
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

    FArrayBox& viscfab=(*localMF[idx])[mfi];

    FArrayBox& cellten=(*localMF[CELLTENSOR_MF])[mfi];
    if (cellten.nComp()!=ntensor)
     amrex::Error("cellten invalid ncomp");

    FArrayBox& velfab=(*vel)[mfi];
    FArrayBox& eosfab=(*EOSdata)[mfi];

    Vector<int> velbc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);

    int fortran_im=im+1;

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

      // in: DERIVE_3D.F90
    FORT_DERTURBVISC(
      &level,
      &fortran_im,
      &nmat,
      &dt_slab,
      &ntensor,  // for declaring cellten
      eosfab.dataPtr(),
      ARLIM(eosfab.loVect()),ARLIM(eosfab.hiVect()),
      voffab.dataPtr(),
      ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),
      velfab.dataPtr(),
      ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
      viscfab.dataPtr(),
      ARLIM(viscfab.loVect()),ARLIM(viscfab.hiVect()),
      cellten.dataPtr(),
      ARLIM(cellten.loVect()),ARLIM(cellten.hiVect()),
      tilelo,tilehi,
      fablo,fabhi,&bfact,
      &cur_time_slab,
      dx,xlo,
      &ngrow,
      &ncomp_visc);
   } //mfi
} // omp
   ns_reconcile_d_num(166);

  } else if (les_model[im]==0) {
   // do nothing
  } else
   amrex::Error("les_model invalid");

  delete gammadot_mf;
 } // im=0..nmat-1

 delete vel;
 delete EOSdata;
 delete tensor;

}  // subroutine getStateVISC


// if FENE-CR+Carreau,
// liquid viscosity=etaS+etaP ( 1+ (beta gamma_dot)^alpha )^((n-1)/alpha)
//
// for each material, there are 4 components:
// 1. \dot{gamma}
// 2. Tr(A) if viscoelastic
//    \dot{gamma} o.t.
// 3. Tr(A) (liquid viscosity - etaS)/etaP  if FENE-CR+Carreau
//    Tr(A) if FENE-CR
//    \dot{gamma} o.t.
// 4. (3) * f(A)  if viscoelastic
//    \dot{gamma} o.t.
// 5. Last component is the vorticity magnitude.

void NavierStokes::getState_tracemag_ALL(int idx) {

 if (level!=0)
  amrex::Error("level invalid getState_tracemag_ALL");

 int finest_level=parent->finestLevel();
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.getState_tracemag(idx);
  int scomp=0;
  int ncomp=ns_level.localMF[idx]->nComp();
  ns_level.avgDown_localMF(idx,scomp,ncomp,0);
 }

} // getState_tracemag_ALL 

//ngrow=1
void NavierStokes::getState_tracemag(int idx) { 
 
 bool use_tiling=ns_tiling;

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int finest_level=parent->finestLevel();
 int nmat=num_materials;
 int ntensor=AMREX_SPACEDIM*AMREX_SPACEDIM;

 if (localMF_grow[idx]>=0)
  amrex::Error("local magtrace not previously deleted");

 int ntrace=5*nmat;
  //ngrow=1
 new_localMF(idx,ntrace,1,-1);

 VOF_Recon_resize(2,SLOPE_RECON_MF);
 debug_ngrow(SLOPE_RECON_MF,2,680);
 if (localMF[SLOPE_RECON_MF]->nComp()!=nmat*ngeom_recon)
  amrex::Error("localMF[SLOPE_RECON_MF]->nComp() invalid");

 debug_ngrow(CELLTENSOR_MF,1,9);
 if (localMF[CELLTENSOR_MF]->nComp()!=ntensor)
  amrex::Error("localMF[CELLTENSOR_MF]->nComp() invalid");

  //ngrow=1
 MultiFab* den_data=getStateDen(1,cur_time_slab);
  //ngrow=2
 MultiFab* vel_data=getState(2,0,AMREX_SPACEDIM,cur_time_slab);

 debug_ngrow(CELL_VISC_MATERIAL_MF,1,9);
 int ncomp_visc=localMF[CELL_VISC_MATERIAL_MF]->nComp();
 if (ncomp_visc!=3*nmat)
  amrex::Error("visc_data invalid ncomp");

 int ncomp_den=den_data->nComp();

 const Real* dx = geom.CellSize();

 for (int im=0;im<nmat;im++) {

  MultiFab* tensor=den_data;
  int allocate_tensor=0;
  if (ns_is_rigid(im)==0) {
   if (elastic_viscosity[im]>0.0) {
    int partid=0;
    while ((im_elastic_map[partid]!=im)&&(partid<im_elastic_map.size())) {
     partid++;
    }
    if (partid<im_elastic_map.size()) {
     int scomp_tensor=partid*NUM_TENSOR_TYPE;
      //ngrow=1
     tensor=getStateTensor(1,scomp_tensor,NUM_TENSOR_TYPE,cur_time_slab);
     allocate_tensor=1;
    } else
     amrex::Error("partid could not be found: getState_tracemag");
   } else if (elastic_viscosity[im]==0.0) {
    // do nothing
   } else
    amrex::Error("elastic_viscosity invalid");
  } else if (ns_is_rigid(im)==1) {
   // do nothing
  } else
   amrex::Error("ns_is_rigid(im) invalid");

  int ncomp_tensor=tensor->nComp();
    
  int idest=5*im;

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(den_data->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(*den_data,use_tiling); mfi.isValid(); ++mfi) {

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
   if (cellten.nComp()!=ntensor)
    amrex::Error("cellten invalid ncomp");

   FArrayBox& destfab=(*localMF[idx])[mfi];

   FArrayBox& velfab=(*vel_data)[mfi];

   Vector<int> bc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   int iproject=0;
   int onlyscalar=1;  // mag(trace gradu)
   int ngrow_getshear=1;
    // declared in: DERIVE_3D.F90
   fort_getshear(
    &ntensor,
    cellten.dataPtr(),
    ARLIM(cellten.loVect()),ARLIM(cellten.hiVect()),
    velfab.dataPtr(),
    ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
    dx,xlo,
    destfab.dataPtr(idest),
    ARLIM(destfab.loVect()),ARLIM(destfab.hiVect()),
    &iproject,&onlyscalar,
    &cur_time_slab,
    tilelo,tilehi,
    fablo,fabhi,
    &bfact,
    &level,
    bc.dataPtr(),
    &ngrow_getshear,
    &nmat);

  } //mfi
} // omp
  ns_reconcile_d_num(167);

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(den_data->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(*den_data,use_tiling); mfi.isValid(); ++mfi) {

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
   if (cellten.nComp()!=ntensor)
    amrex::Error("cellten invalid ncomp");

   FArrayBox& destfab=(*localMF[idx])[mfi];

   FArrayBox& viscfab=(*localMF[CELL_VISC_MATERIAL_MF])[mfi];

   FArrayBox& velfab=(*vel_data)[mfi];
   FArrayBox& denfab=(*den_data)[mfi];
   FArrayBox& tenfab=(*tensor)[mfi];

   Vector<int> bc=getBCArray(State_Type,gridno,0,AMREX_SPACEDIM);

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // 1. gamma dot
    // 2. Tr(A) if viscoelastic, gamma dot otherwise
    // 3. Tr(A) (liquid viscosity - etaS)/etaP  if FENE-CR+Carreau
    //    Tr(A) if FENE-CR
    //    \dot{gamma} o.t.
    // 4. (3) * f(A)  if viscoelastic
    //    \dot{gamma} o.t.
    // 5. magnitude of vorticity

   int ngrow_magtrace=1;

    //declared in: DERIVE_3D.F90
   fort_dermagtrace(
    &level,
    &finest_level,  
    &im, //im=0..nmat-1
    &ntensor,
    cellten.dataPtr(),
    ARLIM(cellten.loVect()),ARLIM(cellten.hiVect()),
    destfab.dataPtr(idest),
    ARLIM(destfab.loVect()),ARLIM(destfab.hiVect()),
    denfab.dataPtr(),
    ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
    tenfab.dataPtr(),
    ARLIM(tenfab.loVect()),ARLIM(tenfab.hiVect()),
    velfab.dataPtr(),
    ARLIM(velfab.loVect()),ARLIM(velfab.hiVect()),
    viscfab.dataPtr(),
    ARLIM(viscfab.loVect()),ARLIM(viscfab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    &ngrow_magtrace,
    dx,xlo,
    &cur_time_slab,
    bc.dataPtr(),
    &ncomp_den,
    &ncomp_tensor,
    &ncomp_visc,
    &ntrace, 
    &nmat,
    polymer_factor.dataPtr(),
    etaS.dataPtr(),
    etaP.dataPtr(),
    Carreau_beta.dataPtr(),
    elastic_time.dataPtr(),
    viscoelastic_model.dataPtr(),
    elastic_viscosity.dataPtr());
  } //mfi
} // omp
  ns_reconcile_d_num(168);

  if (allocate_tensor==0) {
   // do nothing
  } else if (allocate_tensor==1) {
   delete tensor;
  } else
   amrex::Error("allocate_tensor invalid");

 } // im=0..nmat-1

 delete den_data;
 delete vel_data;

}  // getState_tracemag

}/* namespace amrex */
