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
#include <INTEGRATED_QUANTITY.H>
#include <DRAG_COMP.H>
#include <GLOBALUTIL_F.H>
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
#include <INDEX_TYPE_MACROS.H>

#define GEOM_GROW   1
#define bogus_value 1.0e+20

#define BOGUS (1.0E+6)

namespace amrex{

void NavierStokes::Mult_localMF(int idx_dest,int idx_source,
  int scomp,int dcomp,int ncomp,int ngrow) {

 std::string local_caller_string="Mult_localMF";

 debug_ngrow(idx_dest,ngrow,local_caller_string); 
 debug_ngrow(idx_source,ngrow,local_caller_string); 
 if (localMF[idx_dest]->boxArray()!=
     localMF[idx_source]->boxArray())
  amrex::Error("Mult_localMF: boxarrays do not match");

 MultiFab::Multiply(*localMF[idx_dest],*localMF[idx_source],
    scomp,dcomp,ncomp,ngrow);

}  // subroutine Mult_localMF

void NavierStokes::Plus_localMF(int idx_dest,Real val,
  int dcomp,int ncomp,int ngrow) {

 std::string local_caller_string="Plus_localMF";

 debug_ngrow(idx_dest,ngrow,local_caller_string); 
 if (localMF[idx_dest]->boxArray()!=grids)
  amrex::Error("Plus_localMF: boxarrays do not match");

 localMF[idx_dest]->plus(val,dcomp,ncomp,ngrow);

}  // subroutine Plus_localMF


void NavierStokes::Copy_localMF(int idx_dest,int idx_source,
  int scomp,int dcomp,int ncomp,int ngrow) {

 std::string local_caller_string="Copy_localMF";

 debug_ngrow(idx_dest,ngrow,local_caller_string); 
 debug_ngrow(idx_source,ngrow,local_caller_string); 
 if (localMF[idx_dest]->boxArray()!=
     localMF[idx_source]->boxArray())
  amrex::Error("Copy_localMF: boxarrays do not match");

 MultiFab::Copy(*localMF[idx_dest],*localMF[idx_source],scomp,
   dcomp,ncomp,ngrow);

}  // copy_LocalMF


void NavierStokes::minus_localMF(int idx_dest,int idx_source,
  int scomp,int ncomp,int ngrow) {

 std::string local_caller_string="minus_localMF";

 debug_ngrow(idx_dest,ngrow,local_caller_string); 
 debug_ngrow(idx_source,ngrow,local_caller_string); 
 if (localMF[idx_dest]->boxArray()!=
     localMF[idx_source]->boxArray())
  amrex::Error("minus_localMF: boxarrays do not match");

 localMF[idx_dest]->minus(*localMF[idx_source],scomp,ncomp,ngrow);

}  // end subroutine minus_LocalMF

//grid_type=-1 ... 5
void NavierStokes::new_localMF(int idx_MF,int ncomp,int ngrow,int grid_type) {

 std::string local_caller_string="new_localMF";

 if (1==0) {
  std::fflush(NULL);
  int proc=ParallelDescriptor::MyProc();
  std::cout << "in new_localMF idx_MF= " << idx_MF << '\n';
  std::cout << "in new_localMF ncomp= " << ncomp << '\n';
  std::cout << "in new_localMF ngrow= " << ngrow << '\n';
  std::cout << "in new_localMF grid_type= " << grid_type << '\n';
  std::cout << "in new_localMF level= " << level << '\n';
  std::cout << "in new_localMF proc= " << proc << '\n';
  std::fflush(NULL);
  ParallelDescriptor::Barrier();
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

 if (1==0) {
  std::fflush(NULL);
  int proc=ParallelDescriptor::MyProc();
  std::cout << "in new_localMF2 idx_MF= " << idx_MF << '\n';
  std::cout << "in new_localMF2 ncomp= " << ncomp << '\n';
  std::cout << "in new_localMF2 ngrow= " << ngrow << '\n';
  std::cout << "in new_localMF2 grid_type= " << grid_type << '\n';
  std::cout << "in new_localMF2 level= " << level << '\n';
  std::cout << "in new_localMF2 proc= " << proc << '\n';
  std::fflush(NULL);
  ParallelDescriptor::Barrier();
 }


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

 ParallelDescriptor::Barrier();

 if (1==0) {
  std::fflush(NULL);
  int proc=ParallelDescriptor::MyProc();
  std::cout << "in new_localMF3 idx_MF= " << idx_MF << '\n';
  std::cout << "in new_localMF3 ncomp= " << ncomp << '\n';
  std::cout << "in new_localMF3 ngrow= " << ngrow << '\n';
  std::cout << "in new_localMF3 grid_type= " << grid_type << '\n';
  std::cout << "in new_localMF3 level= " << level << '\n';
  std::cout << "in new_localMF3 proc= " << proc << '\n';
  std::fflush(NULL);
  ParallelDescriptor::Barrier();
 }
 localMF[idx_MF]=new MultiFab(edge_boxes,dmap,ncomp,ngrow,
	MFInfo().SetTag("localMF[idx_MF]"),FArrayBoxFactory());

 if (1==0) {
  std::fflush(NULL);
  int proc=ParallelDescriptor::MyProc();
  std::cout << "in new_localMF3 idx_MF= " << idx_MF << '\n';
  std::cout << "in new_localMF3 ncomp= " << ncomp << '\n';
  std::cout << "in new_localMF3 ngrow= " << ngrow << '\n';
  std::cout << "in new_localMF3 grid_type= " << grid_type << '\n';
  std::cout << "in new_localMF3 level= " << level << '\n';
  std::cout << "in new_localMF3 proc= " << proc << '\n';
  std::fflush(NULL);
  ParallelDescriptor::Barrier();
 }

 ParallelDescriptor::Barrier();
 localMF[idx_MF]->setVal(0.0,0,ncomp,ngrow);
 localMF_grow[idx_MF]=ngrow;

 if (1==0) {
  std::fflush(NULL);
  int proc=ParallelDescriptor::MyProc();
  std::cout << "in new_localMF4 idx_MF= " << idx_MF << '\n';
  std::cout << "in new_localMF4 ncomp= " << ncomp << '\n';
  std::cout << "in new_localMF4 ngrow= " << ngrow << '\n';
  std::cout << "in new_localMF4 grid_type= " << grid_type << '\n';
  std::cout << "in new_localMF4 level= " << level << '\n';
  std::cout << "in new_localMF4 proc= " << proc << '\n';
  std::fflush(NULL);
  ParallelDescriptor::Barrier();
 }

 ParallelDescriptor::Barrier();
  //declared in: NavierStokes.cpp
 debug_ixType(idx_MF,grid_type,local_caller_string);

 if (1==0) {
  std::fflush(NULL);
  int proc=ParallelDescriptor::MyProc();
  std::cout << "in new_localMF5 idx_MF= " << idx_MF << '\n';
  std::cout << "in new_localMF5 ncomp= " << ncomp << '\n';
  std::cout << "in new_localMF5 ngrow= " << ngrow << '\n';
  std::cout << "in new_localMF5 grid_type= " << grid_type << '\n';
  std::cout << "in new_localMF5 level= " << level << '\n';
  std::cout << "in new_localMF5 proc= " << proc << '\n';
  std::fflush(NULL);
  ParallelDescriptor::Barrier();
 }

} //end subroutine new_localMF

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


//dir=0..sdim-1
void NavierStokes::getStateMAC_localMF(
  int idx_MF,int ngrow,int dir,Real time) {

 if ((dir>=0)&&(dir<AMREX_SPACEDIM)) {
  // do nothing
 } else
  amrex::Error("dir invalid");

 if (localMF_grow[idx_MF]==-1) {
  // do nothing
 } else {
  std::cout << "idx_MF= " << idx_MF << " ngrow = " << ngrow << 
   " dir= " << dir << " time= " << time << '\n';
  amrex::Error("localMF_grow invalid getstateMAC_localMF");
 }
 if (ngrow<0)
  amrex::Error("ngrow invalid");

   // NavierStokes::getStateMAC declared in NavierStokes.cpp
 localMF[idx_MF]=getStateMAC(ngrow,dir,time);
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
  Real time,
  const std::string& caller_string) {

 if ((ngrow<0)||(ngrow>ngrow_distance))
  amrex::Error("ngrow invalid");

 if (localMF_grow[idx_MF]==-1) {
  // do nothing
 } else {
  std::cout << "idx_MF,time " << idx_MF << ' ' << time << '\n';
  amrex::Error("localMF_grow invalid getStateDist_localMF ");
 }

 localMF[idx_MF]=getStateDist(ngrow,time,caller_string);
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

} //end subroutine getStateTensor_localMF


void NavierStokes::getStateRefineDensity_localMF(int idx_MF,int ngrow,
  int scomp,int ncomp,Real time) {

 if (localMF_grow[idx_MF]==-1) {
  // do nothing
 } else
  amrex::Error("localMF_grow invalid getStateRefineDensity_localMF");

 if (ngrow<0)
  amrex::Error("ngrow invalid");

 localMF[idx_MF]=getStateRefineDensity(ngrow,scomp,ncomp,time);
 localMF_grow[idx_MF]=ngrow;

} //end subroutine getStateRefineDensity_localMF

void NavierStokes::maskfiner_localMF(int idx_MF,int ngrow,
  Real tag,int clearbdry) {

 delete_localMF_if_exist(idx_MF,1);
 if ((ngrow<0)||(ngrow>ngrow_distance))
  amrex::Error("ngrow invalid");

 localMF[idx_MF]=maskfiner(ngrow,tag,clearbdry);
 localMF_grow[idx_MF]=ngrow;
}  // end subroutine maskfiner_localMF

//called from:
//NavierStokes::writeTECPLOT_File
//NavierStokes::init_gradu_tensor_and_material_visc_ALL
void NavierStokes::getStateVISC_ALL(const std::string& caller_string) {

 if (level!=0)
  amrex::Error("level invalid getStateVISC_ALL");

 std::string local_caller_string="getStateVISC_ALL";
 local_caller_string=caller_string+local_caller_string;

 if (divu_outer_sweeps+1==num_divu_outer_sweeps) {
  //do nothing
 } else if (divu_outer_sweeps+1<num_divu_outer_sweeps) {

  if ((slab_step>=0)&&(slab_step<ns_time_order)) {

   if (pattern_test(local_caller_string,"writeTECPLOT_File")==1) {
    //do nothing
   } else if (pattern_test(local_caller_string,"volWgtSumALL")==1) {
    //do nothing
   } else if (pattern_test(local_caller_string,"prepare_post_process")==1) {
    //do nothing
   } else if (pattern_test(local_caller_string,"nucleation_code_segment")==1) {
    //do nothing
   } else if (pattern_test(local_caller_string,"do_the_advance")==1) {
    //do nothing
   } else
    amrex::Error("local_caller_string invalid");

  } else if ((slab_step==-1)||(slab_step==ns_time_order)) {
   //do nothing
  } else
   amrex::Error("slab_step invalid");

 } else
  amrex::Error("divu_outer_sweeps invalid");

 int finest_level=parent->finestLevel();
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.getStateVISC(local_caller_string);
  int scomp=0;
  int ncomp=ns_level.localMF[CELL_VISC_MATERIAL_MF]->nComp();
  ns_level.avgDown_localMF(CELL_VISC_MATERIAL_MF,scomp,ncomp,
    LOW_ORDER_AVGDOWN);
 }

} // end subroutine getStateVISC_ALL 

//called from: NavierStokes::init_FSI_GHOST_MAC_MF_ALL
//CELL_CONDUCTIVITY_MATERIAL_MF is deleted in ::Geometry_cleanup()
void NavierStokes::getStateCONDUCTIVITY_ALL() {

 if (level!=0)
  amrex::Error("level invalid getStateCONDUCTIVITY_ALL");

 int finest_level=parent->finestLevel();

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.getStateCONDUCTIVITY();
  int scomp=0;
  int ncomp=ns_level.localMF[CELL_CONDUCTIVITY_MATERIAL_MF]->nComp();
   // spectral_override==1 => order derived from "enable_spectral"
   // spectral_override==0 => always low order.
  ns_level.avgDown_localMF(CELL_CONDUCTIVITY_MATERIAL_MF,
    scomp,ncomp,LOW_ORDER_AVGDOWN);
 }

} // end subroutine getStateCONDUCTIVITY_ALL 


void NavierStokes::delete_localMF(int idx_MF,int ncomp) {

 for (int scomp=idx_MF;scomp<idx_MF+ncomp;scomp++) {
  if (localMF_grow[scomp]>=0) {
   // do nothing
  } else {
   std::cout << "level= " << level << '\n';
   std::cout << "idx_MF= " << idx_MF << '\n';
   std::cout << "ncomp= " << ncomp << '\n';
   amrex::Error("forgot to allocate the localMF variable before delete");
  }
  delete localMF[scomp];
  ParallelDescriptor::Barrier();
  localMF_grow[scomp]=-1;
  localMF[scomp]=0;
  ParallelDescriptor::Barrier();
 }  // scomp

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
  Array4<Real> const& fab_array=fab.array();

  Box d_box_interior=fab.box(); //tag (done second)
  Box d_box_exterior=fab.box(); //1.0-tag (done first)

  if (clear_phys_boundary==0) {
   d_box_interior=fab.box(); //tag
  } else if (clear_phys_boundary==1) {
   d_box_interior=geom.Domain();
   d_box_interior &= fab.box();
  } else if ((clear_phys_boundary==2)||
             (clear_phys_boundary==3)) {
   d_box_interior=geom.Domain();
   d_box_interior &= grids[mfi.index()];
  } else
   amrex::Error("clear_phys_boundary invalid");

  const Dim3 lo3=amrex::lbound(d_box_exterior);
  const Dim3 hi3=amrex::ubound(d_box_exterior);
  for (int z=lo3.z;z<=hi3.z;++z) {
  for (int y=lo3.y;y<=hi3.y;++y) {
  for (int x=lo3.x;x<=hi3.x;++x) {
   fab_array(x,y,z,0)=1.0-tag;
  }
  }
  }

  const Dim3 lo3b=amrex::lbound(d_box_interior);
  const Dim3 hi3b=amrex::ubound(d_box_interior);
  for (int z=lo3b.z;z<=hi3b.z;++z) {
  for (int y=lo3b.y;y<=hi3b.y;++y) {
  for (int x=lo3b.x;x<=hi3b.x;++x) {
   fab_array(x,y,z,0)=tag;
  }
  }
  }

 } // mfi
} //omp
 ns_reconcile_d_num(LOOP_MASKFINER,"maskfiner");

 if ((clear_phys_boundary==0)|| //mask=tag outside domain and coarse/fine ghost
     (clear_phys_boundary==1)|| //mask=1-tag outside domain,mask=tag c/f ghost
     (clear_phys_boundary==2)) {//" ",mask=1-tag c/f ghost,mask=tag f/f ghost

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
     if (c_box.ok()) {
      Array4<Real> const& fab_array=fab.array();
      const Dim3 lo3=amrex::lbound(c_box);
      const Dim3 hi3=amrex::ubound(c_box);
      for (int z=lo3.z;z<=hi3.z;++z) {
      for (int y=lo3.y;y<=hi3.y;++y) {
      for (int x=lo3.x;x<=hi3.x;++x) {
       fab_array(x,y,z,0)=1.0-tag;
      }
      }
      }
     }

    } // j
   } // mfi
} //omp
   ns_reconcile_d_num(LOOP_MASKFINER_UNTAG,"maskfiner");
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
} // end subroutine maskfiner

MultiFab* NavierStokes::maskfiner_noghosts(Real tag) {

 int ngrow=0;

 int finest_level = parent->finestLevel();
 MultiFab* mask=new MultiFab(grids,dmap,1,ngrow,
  MFInfo().SetTag("mask"),FArrayBoxFactory());

  // sanity check for mpi and open mp
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

   // sanity check
  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  FArrayBox& fab = (*mask)[mfi];
  Array4<Real> const& fab_array=fab.array();

  Box d_box_interior=fab.box(); 

  const Dim3 lo3b=amrex::lbound(d_box_interior);
  const Dim3 hi3b=amrex::ubound(d_box_interior);
  for (int z=lo3b.z;z<=hi3b.z;++z) {
  for (int y=lo3b.y;y<=hi3b.y;++y) {
  for (int x=lo3b.x;x<=hi3b.x;++x) {
   fab_array(x,y,z,0)=tag;
  }
  }
  }

 } // mfi
} //omp

  //sanity check
 ns_reconcile_d_num(LOOP_MASKFINER,"maskfiner");

 if (level<finest_level) {

  NavierStokes& ns_fine=getLevel(level+1);
  const BoxArray& f_box=ns_fine.grids;

   //sanity check
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
    if (c_box.ok()) {
     Array4<Real> const& fab_array=fab.array();
     const Dim3 lo3=amrex::lbound(c_box);
     const Dim3 hi3=amrex::ubound(c_box);
     for (int z=lo3.z;z<=hi3.z;++z) {
     for (int y=lo3.y;y<=hi3.y;++y) {
     for (int x=lo3.x;x<=hi3.x;++x) {
      fab_array(x,y,z,0)=1.0-tag;
     }
     }
     }
    }

   } // j
  } // mfi
} //omp
   //sanity check
  ns_reconcile_d_num(LOOP_MASKFINER_UNTAG,"maskfiner");
 } // level<finest_level

 return mask;
} // end subroutine maskfiner_noghosts



// if spectral_override==0, then always low order average down.
void NavierStokes::avgDown_localMF(int idxMF,int scomp,int ncomp,
  int spectral_override) {

 std::string local_caller_string="avgDown_localMF";

 int finest_level=parent->finestLevel();
 if (level<finest_level) {
  NavierStokes& ns_fine=getLevel(level+1);
  debug_ngrow(idxMF,0,local_caller_string);
  ns_fine.debug_ngrow(idxMF,0,local_caller_string);
  MultiFab& S_crse=*localMF[idxMF];
  MultiFab& S_fine=*ns_fine.localMF[idxMF];
  avgDown(S_crse,S_fine,scomp,ncomp,spectral_override);
 }

} // end subroutine avgDown_localMF


void NavierStokes::avgDown_tag_localMF(int idxMF) {

 std::string local_caller_string="avgDown_tag_localMF";

 int finest_level=parent->finestLevel();
 if (level<finest_level) {
  NavierStokes& ns_fine=getLevel(level+1);
  debug_ngrow(idxMF,0,local_caller_string);
  ns_fine.debug_ngrow(idxMF,0,local_caller_string);
  MultiFab& S_crse=*localMF[idxMF];
  MultiFab& S_fine=*ns_fine.localMF[idxMF];
  level_avgDown_tag(S_crse,S_fine);
 }

} // end subroutine avgDown_tag_localMF


void NavierStokes::avgDownBURNING_localMF(
  int burnvel_MF,int TSAT_MF) {

 std::string local_caller_string="avgDownBURNING_localMF";

 int finest_level=parent->finestLevel();
 int nburning=EXTRAP_NCOMP_BURNING;
 int ntsat=EXTRAP_NCOMP_TSAT;

 if (level<finest_level) {
  NavierStokes& ns_fine=getLevel(level+1);

  debug_ngrow(burnvel_MF,0,local_caller_string);
  ns_fine.debug_ngrow(burnvel_MF,0,local_caller_string);
  MultiFab& S_crseB=*localMF[burnvel_MF];
  MultiFab& S_fineB=*ns_fine.localMF[burnvel_MF];
  if ((S_crseB.nComp()==nburning)&&
      (S_fineB.nComp()==nburning)) {
   int velflag=1;
   level_avgDownBURNING(S_crseB,S_fineB,velflag);
  } else
   amrex::Error("S_crseB or S_fineB invalid nComp");

  debug_ngrow(TSAT_MF,0,local_caller_string);
  ns_fine.debug_ngrow(TSAT_MF,0,local_caller_string);
  MultiFab& S_crseT=*localMF[TSAT_MF];
  MultiFab& S_fineT=*ns_fine.localMF[TSAT_MF];
  if ((S_crseT.nComp()==ntsat)&&
      (S_fineT.nComp()==ntsat)) {
   int velflag=0;
   level_avgDownBURNING(S_crseT,S_fineT,velflag);
  } else
   amrex::Error("S_crseT or S_fineT invalid nComp");

 } // level<finest_level

} // end subroutine avgDownBURNING_localMF

void NavierStokes::avgDownDRAG_MF() {

 std::string local_caller_string="avgDownDRAG_MF";

  //ngrow_make_distance=3
  //ngrow_distance=4
 debug_ngrow(DRAG_MF,ngrow_make_distance,local_caller_string);
 debug_ixType(DRAG_MF,-1,local_caller_string);
 if (localMF[DRAG_MF]->nComp()==N_DRAG) {
  // do nothing
 } else 
  amrex::Error("DRAG_MF invalid ncomp");

 int finest_level=parent->finestLevel();

 if (level<finest_level) {
  NavierStokes& ns_fine=getLevel(level+1);

  debug_ngrow(DRAG_MF,ngrow_make_distance,local_caller_string);
  ns_fine.debug_ngrow(DRAG_MF,ngrow_make_distance,local_caller_string);
  MultiFab& S_crseD=*localMF[DRAG_MF];
  MultiFab& S_fineD=*ns_fine.localMF[DRAG_MF];
  if ((S_crseD.nComp()==N_DRAG)&&
      (S_fineD.nComp()==N_DRAG)) {
   level_avgDownDRAG(S_crseD,S_fineD);
  } else
   amrex::Error("S_crseD or S_fineD invalid nComp");

 } // level<finest_level

} // end subroutine avgDownDRAG_MF

void NavierStokes::avgDownCURV_localMF(int idxMF) {

 std::string local_caller_string="avgDownCURV_localMF";

 int finest_level=parent->finestLevel();
 if (level<finest_level) {
  NavierStokes& ns_fine=getLevel(level+1);
  debug_ngrow(idxMF,0,local_caller_string);
  ns_fine.debug_ngrow(idxMF,0,local_caller_string);
  MultiFab& S_crse=*localMF[idxMF];
  MultiFab& S_fine=*ns_fine.localMF[idxMF];
  level_avgDownCURV(S_crse,S_fine);
 }

} // end subroutine avgDownCURV_localMF


// flux variables: average down in the tangential direction 
// to the box face, copy in
// the normal direction.  Since the blocking factor is >=2, it is
// impossible to have a grid box with size of 1 cell width.
void NavierStokes::avgDown_and_Copy_localMF(
  int idx_den_MF,
  int idx_vel_MF,
  int idx_flux_MF,
  int operation_flag) {

 std::string local_caller_string="avgDown_and_Copy_localMF";
 int finest_level=parent->finestLevel();

 int ncomp_den=0;
 int ncomp_vel=0;
 int ncomp_flux=0;
 int scomp_flux=0;

 if (operation_flag==OP_ISCHEME_MAC) {  // advection 
  ncomp_den=num_materials*num_state_material;
  ncomp_vel=AMREX_SPACEDIM;
  ncomp_flux=NFLUXSEM;
  scomp_flux=0;
 } else if ((operation_flag==OP_UNEW_CELL_TO_MAC)|| 
            (operation_flag==OP_UMAC_PLUS_VISC_CELL_TO_MAC)||
	    (operation_flag==OP_U_SEM_CELL_MAC_TO_MAC)) {
  ncomp_den=AMREX_SPACEDIM;
  ncomp_vel=AMREX_SPACEDIM;
  ncomp_flux=1;
  scomp_flux=0;
 } else if (operation_flag==OP_PRESGRAD_MAC) { //grad p
  ncomp_den=1;
  ncomp_vel=1;
  ncomp_flux=1;
  scomp_flux=0;
 } else if (operation_flag==OP_UGRAD_COUPLING_MAC) {  
  ncomp_den=AMREX_SPACEDIM;
  ncomp_vel=AMREX_SPACEDIM;
  ncomp_flux=AMREX_SPACEDIM;
  scomp_flux=0;
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
  debug_ngrow(VOLUME_MF,0,local_caller_string);
  fine_lev.resize_metrics(1);
  fine_lev.debug_ngrow(VOLUME_MF,0,local_caller_string);

  debug_ngrow(LEVELPC_MF,1,local_caller_string);
  if (localMF[LEVELPC_MF]->nComp()!=num_materials*(1+AMREX_SPACEDIM))
   amrex::Error("localMF[LEVELPC_MF]->nComp()!=num_materials*(1+AMREX_SPACEDIM)");

  fine_lev.debug_ngrow(LEVELPC_MF,1,local_caller_string);
  if (fine_lev.localMF[LEVELPC_MF]->nComp()!=num_materials*(1+AMREX_SPACEDIM))
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
   MultiFab crse_S_fine_MAC(crse_S_fine_BA_MAC,crse_dmap,ncomp_flux,0,
     MFInfo().SetTag("crse_S_fine_MAC"),FArrayBoxFactory());
   crse_S_fine_MAC.setVal(1.0e+30);

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

     // declared in: NAVIERSTOKES_3D.F90
    fort_avgdown_copy( 
     &enable_spectral,
     &finest_level,
     &operation_flag,
     &dir,
     prob_lo,
     dxf,
     &level,&f_level,
     &bfact_c,&bfact_f,     
     xlo_fine,dx,
     &ncomp_flux,
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
   ns_reconcile_d_num(LOOP_AVGDOWN_COPY,"avgDown_and_Copy_localMF");

   S_crse_MAC.ParallelCopy(crse_S_fine_MAC,0,scomp_flux,ncomp_flux);
   ParallelDescriptor::Barrier();

   const Box& domain = geom.Domain();
   if (geom.isPeriodic(dir)) {
    IntVect pshift=IntVect::TheZeroVector();
    pshift[dir]=domain.length(dir);
    crse_S_fine_MAC.shift(pshift);

    ParallelDescriptor::Barrier();
    S_crse_MAC.ParallelCopy(crse_S_fine_MAC,0,scomp_flux,ncomp_flux);
    ParallelDescriptor::Barrier();

    pshift[dir]=-2*domain.length(dir);
    crse_S_fine_MAC.shift(pshift);

    S_crse_MAC.ParallelCopy(crse_S_fine_MAC,0,scomp_flux,ncomp_flux);
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

 std::string local_caller_string="interp_and_Copy_localMF";

 int finest_level=parent->finestLevel();

 int ncomp_den=0;
 int ncomp_vel=0;
 int ncomp_flux=0;
 int scomp_flux=0;

 if (operation_flag==OP_UGRAD_COUPLING_MAC) {  // viscosity
  ncomp_den=AMREX_SPACEDIM;
  ncomp_vel=AMREX_SPACEDIM;
  ncomp_flux=AMREX_SPACEDIM;
  scomp_flux=0;
 } else if (operation_flag==OP_ISCHEME_MAC) {  // advection 
  ncomp_den=num_materials*num_state_material;
  ncomp_vel=AMREX_SPACEDIM;
  ncomp_flux=NFLUXSEM;
  scomp_flux=0;
 } else if ((operation_flag==OP_UNEW_CELL_TO_MAC)||
            (operation_flag==OP_UMAC_PLUS_VISC_CELL_TO_MAC)|| 
	    (operation_flag==OP_U_SEM_CELL_MAC_TO_MAC)) {
  ncomp_den=AMREX_SPACEDIM;
  ncomp_vel=AMREX_SPACEDIM;
  ncomp_flux=1;
  scomp_flux=0;
 } else if (operation_flag==OP_PRESGRAD_MAC) {
  ncomp_den=1;
  ncomp_vel=1;
  ncomp_flux=1;
  scomp_flux=0;
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
  debug_ngrow(VOLUME_MF,0,local_caller_string);
  coarse_lev.resize_metrics(1);
  coarse_lev.debug_ngrow(VOLUME_MF,0,local_caller_string);

  resize_mask_nbr(1);
   // (1) =1 interior  =1 fine-fine ghost in domain  =0 otherwise
   // (2) =1 interior  =0 otherwise
  debug_ngrow(MASK_NBR_MF,1,local_caller_string); 

  debug_ngrow(LEVELPC_MF,1,local_caller_string);
  if (localMF[LEVELPC_MF]->nComp()!=num_materials*(1+AMREX_SPACEDIM))
   amrex::Error("localMF[LEVELPC_MF]->nComp()!=num_materials*(1+AMREX_SPACEDIM)");

  coarse_lev.debug_ngrow(LEVELPC_MF,1,local_caller_string);
  if (coarse_lev.localMF[LEVELPC_MF]->nComp()!=num_materials*(1+AMREX_SPACEDIM))
   amrex::Error("coarse_lev.localMF[LEVELPC_MF]->nComp() invalid");

  MultiFab* coarse_den_fine;
  MultiFab* coarse_vel_fine;
  MultiFab* coarse_mask_sem_fine;
  MultiFab* coarse_LS_fine;

  DistributionMapping crse_dmap=fdmap;
  coarse_mask_sem_fine=new MultiFab(crse_S_fine_BA,crse_dmap,1,1,
	MFInfo().SetTag("coarse_mask_sem_fine"),FArrayBoxFactory());
  coarse_LS_fine=new MultiFab(crse_S_fine_BA,crse_dmap,num_materials,1,
	MFInfo().SetTag("coarse_LS_fine"),FArrayBoxFactory());
   // FabArray.H     
   // scomp,dcomp,ncomp,s_nghost,d_nghost
  coarse_mask_sem_fine->ParallelCopy(*coarse_lev.localMF[MASKSEM_MF],0,0,
    1,1,1,geom.periodicity());
  coarse_LS_fine->ParallelCopy(*coarse_lev.localMF[LEVELPC_MF],0,0,
    num_materials,1,1,geom.periodicity());

  coarse_den_fine=new MultiFab(crse_S_fine_BA,crse_dmap,ncomp_den,1,
		  MFInfo().SetTag("coarse_den_fine"),FArrayBoxFactory());
  coarse_vel_fine=new MultiFab(crse_S_fine_BA,crse_dmap,ncomp_vel,1,
		  MFInfo().SetTag("coarse_vel_fine"),FArrayBoxFactory());
    // FabArray.H     
    // scomp,dcomp,ncomp,s_nghost,d_nghost
  coarse_den_fine->ParallelCopy(*coarse_lev.localMF[idx_den_MF],0,0,
    ncomp_den,1,1,geom.periodicity());
  coarse_vel_fine->ParallelCopy(*coarse_lev.localMF[idx_vel_MF],0,0,
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
    Vector<int> velbc=getBCArray(State_Type,gridno,
      STATECOMP_VEL,STATE_NCOMP_VEL);

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     // declared in: NAVIERSTOKES_3D.F90
    fort_interp_copy( 
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
     &ncomp_flux,
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
   ns_reconcile_d_num(LOOP_INTERP_COPY,"interp_and_Copy_localMF");
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

 std::string local_caller_string="sync_flux_var";

 int finest_level=parent->finestLevel();

 if ((dir>=0)&&(dir<AMREX_SPACEDIM)) {
  
  if (localMF[AREA_MF+dir]->boxArray()==
      localMF[flux_MF+dir]->boxArray()) {
   if (localMF[flux_MF+dir]->nComp()==ncomp_flux) {

    debug_ngrow(flux_MF+dir,0,local_caller_string);
    resize_mask_nbr(1);
     // (1) =1 interior  =1 fine-fine ghost in domain  =0 otherwise
     // (2) =1 interior  =0 otherwise
    debug_ngrow(MASK_NBR_MF,1,local_caller_string); 
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
      Vector<int> presbc=getBCArray(State_Type,gridno,STATECOMP_PRES,1);

      int tid_current=ns_thread();
      if ((tid_current<0)||(tid_current>=thread_class::nthreads))
       amrex::Error("tid_current invalid");
      thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

       // fort_fillbdry_flux is declared in: NAVIERSTOKES_3D.F90
      fort_fillbdry_flux( 
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
     ns_reconcile_d_num(LOOP_FILLBDRY_FLUX,"sync_flux_var");

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

 std::string local_caller_string="interp_flux_localMF";

 int finest_level=parent->finestLevel();

 int ncomp_flux=NFLUXSEM;

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
  debug_ngrow(VOLUME_MF,0,local_caller_string);
  coarse_lev.resize_metrics(1);
  coarse_lev.debug_ngrow(VOLUME_MF,0,local_caller_string);

  resize_mask_nbr(1);
   // (1) =1 interior  =1 fine-fine ghost in domain  =0 otherwise
   // (2) =1 interior  =0 otherwise
  debug_ngrow(MASK_NBR_MF,1,local_caller_string); 

  MultiFab* coarse_mask_sem_fine;
  DistributionMapping crse_dmap=fdmap;
  coarse_mask_sem_fine=new MultiFab(crse_S_fine_BA,crse_dmap,1,1,
	MFInfo().SetTag("coarse_mask_sem_fine"),FArrayBoxFactory());
   // FabArray.H     
   // scomp,dcomp,ncomp,s_nghost,d_nghost
  coarse_mask_sem_fine->ParallelCopy(*coarse_lev.localMF[MASKSEM_MF],0,0,
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

   crse_S_fine_MAC.ParallelCopy(*coarse_lev.localMF[coarse_flux_MF+dir],0,0,
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
    Vector<int> velbc=getBCArray(State_Type,gridno,
      STATECOMP_VEL,STATE_NCOMP_VEL);

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     //declared in: NAVIERSTOKES_3D.F90
    fort_interp_flux( 
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
   ns_reconcile_d_num(LOOP_INTERP_FLUX,"interp_flux_localMF");
  } // dir=0..sdim-1

  delete coarse_mask_sem_fine;

 } else if (level==0) {
  // do nothing
 } else
  amrex::Error("level invalid11");

} // subroutine interp_flux_localMF

// interpolate from level+1 to level.
// spectral_override==LOW_ORDER_AVGDOWN => always do low order average down.
void NavierStokes::avgDownEdge_localMF(
  int idxMF,int scomp,int ncomp,
  int start_grid_type,int n_grid_type,
  int spectral_override,
  const std::string& caller_string) {

 if (1==0) {
  std::cout << "avgDownEdge_localMF caller_string= " << 
   caller_string << " start_grid_type= " << start_grid_type <<
   " n_grid_type= " << n_grid_type << '\n';
 }

 std::string local_caller_string="avgDownEdge_localMF";
 local_caller_string=caller_string+local_caller_string;

 int finest_level=parent->finestLevel();
 if (level<finest_level) {
  NavierStokes& ns_fine=getLevel(level+1);

  int local_grid_type=start_grid_type;

  if (n_grid_type==-1) {
   debug_ngrow(idxMF,0,local_caller_string);
   ns_fine.debug_ngrow(idxMF,0,local_caller_string);
   MultiFab& S_crse=*localMF[idxMF];
   MultiFab& S_fine=*ns_fine.localMF[idxMF];
   avgDownEdge(local_grid_type,
	S_crse,S_fine,
	scomp,ncomp,
	spectral_override,
        local_caller_string);
  } else if ((n_grid_type==1)||(n_grid_type==AMREX_SPACEDIM)) {

   for (local_grid_type=start_grid_type;
        local_grid_type<start_grid_type+n_grid_type;
        local_grid_type++) {
    debug_ngrow(idxMF+local_grid_type,0,local_caller_string);
    ns_fine.debug_ngrow(idxMF+local_grid_type,0,local_caller_string);
    MultiFab& S_crse=*localMF[idxMF+local_grid_type];
    MultiFab& S_fine=*ns_fine.localMF[idxMF+local_grid_type];
    avgDownEdge(local_grid_type,
	S_crse,S_fine,
	scomp,ncomp,
	spectral_override,
        local_caller_string);
   } // local_grid_type=start_grid_type...start_grid_type+n_grid_type-1
  } else
   amrex::Error("n_grid_type invalid");
 } else if (level==finest_level) {
  // do nothing
 } else {
  amrex::Error("level invalid");
 }

} // avgDownEdge_localMF

void NavierStokes::CELL_GRID_ELASTIC_FORCE(int im_viscoelastic,
  int elastic_force_mac_grid) {

 std::string local_caller_string="CELL_GRID_ELASTIC_FORCE";

 if ((im_viscoelastic>=0)&&(im_viscoelastic<num_materials)) {
   if (ns_is_rigid(im_viscoelastic)==0) {
    if ((elastic_time[im_viscoelastic]>0.0)&&
        (elastic_viscosity[im_viscoelastic]>0.0)) {
     if (store_elastic_data[im_viscoelastic]==1) {
      // do nothing
     } else
      amrex::Error("expecting store_elastic_data[im_viscoelastic]==1");
    } else
     amrex::Error("expecting elastic_time>0 and elastic_viscosity>0");
   } else
    amrex::Error("expecting ns_is_rigid(im_viscoelastic)==0)"); 
 } else
  amrex::Error("im_viscoelastic invalid");

 int partid=0;
 while ((im_viscoelastic_map[partid]!=im_viscoelastic)&&
        (partid<im_viscoelastic_map.size())) {
  partid++;
 }
 if ((partid>=0)&&(partid<im_viscoelastic_map.size())) {
  // do nothing
 } else
  amrex::Error("partid invalid");

 resize_levelset(2,LEVELPC_MF);
 debug_ngrow(LEVELPC_MF,2,local_caller_string);
 if (localMF[LEVELPC_MF]->nComp()!=num_materials*(AMREX_SPACEDIM+1))
  amrex::Error("(localMF[LEVELPC_MF]->nComp() invalid");

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

 debug_ngrow(VOLUME_MF,1,local_caller_string);
   // maskcoef=1 if not covered by finer level.
 debug_ngrow(MASKCOEF_MF,1,local_caller_string); 
 debug_ngrow(MASK_NBR_MF,1,local_caller_string); // mask_nbr=1 at fine-fine bc.
 debug_ngrow(CELL_VISC_MATERIAL_MF,1,local_caller_string);

 debug_ngrow(CELL_DEN_MF,1,local_caller_string);
 if (localMF[CELL_DEN_MF]->nComp()!=1)
  amrex::Error("localMF[CELL_DEN_MF]->nComp() invalid");

 const Box& domain = geom.Domain();
 const int* domlo = domain.loVect();
 const int* domhi = domain.hiVect();

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

 int nstate=STATE_NCOMP;
 if (nstate!=S_new.nComp())
  amrex::Error("nstate invalid");

 const Real* dx = geom.CellSize();

 int dir=0;

 MultiFab* save_mac_new[AMREX_SPACEDIM];
 MultiFab* save_mac_old[AMREX_SPACEDIM];
 for (dir=0;dir<AMREX_SPACEDIM;dir++) {
  save_mac_new[dir]=getStateMAC(0,dir,cur_time_slab);
  save_mac_old[dir]=getStateMAC(0,dir,cur_time_slab);
 }

  // outer loop: each force component u_t = F_elastic/density  u,v,w
 for (int force_dir=0;force_dir<AMREX_SPACEDIM;force_dir++) {

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
    
   FArrayBox& tensorfab=(*localMF[VISCOTEN_MF])[mfi];

   FArrayBox& rhoinversefab=(*localMF[CELL_DEN_MF])[mfi];

    //FACECOMP_FACEDEN component has 1/rho
   FArrayBox& xfacevar=(*localMF[FACE_VAR_MF+force_dir])[mfi];

    // output
   FArrayBox& SNEWfab=S_new[mfi];

   FArrayBox& xvel_new=(*save_mac_new[force_dir])[mfi];
   FArrayBox& xvel_old=(*save_mac_old[force_dir])[mfi];

   // mask=1.0 at interior fine bc ghost cells
   FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi];
   // maskcoef=1 if not covered by finer level or outside domain
   FArrayBox& maskcoef=(*localMF[MASKCOEF_MF])[mfi]; 

   FArrayBox& levelpcfab=(*localMF[LEVELPC_MF])[mfi];

   FArrayBox& viscfab=(*localMF[CELL_VISC_MATERIAL_MF])[mfi];
   int ncomp_visc=viscfab.nComp();

   Vector<int> velbc=getBCArray(State_Type,gridno,
     STATECOMP_VEL,STATE_NCOMP_VEL);

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // declared in: GODUNOV_3D.F90
   fort_elastic_force(
     &elastic_force_mac_grid,
     &im_viscoelastic, // 0..num_materials-1
     &partid, //0..num_materials_viscoelastic-1
     &force_dir, // force_dir=0,1,..sdim-1  
     &ncomp_visc, 
     &visc_coef,
     velbc.dataPtr(),
     &dt_slab, //fort_elastic_force
     &cur_time_slab,
     xlo,dx,
     viscfab.dataPtr(),
     ARLIM(viscfab.loVect()),ARLIM(viscfab.hiVect()),
     maskfab.dataPtr(), // mask=1.0 at interior fine bc ghost cells
     ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
     //maskcoef=1 if not covered by finer level or outside
     maskcoef.dataPtr(), 
     ARLIM(maskcoef.loVect()),ARLIM(maskcoef.hiVect()),
     levelpcfab.dataPtr(),
     ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),
     rhoinversefab.dataPtr(),
     ARLIM(rhoinversefab.loVect()),ARLIM(rhoinversefab.hiVect()),
     xfacevar.dataPtr(FACECOMP_FACEDEN), // 1/rho
     ARLIM(xfacevar.loVect()),ARLIM(xfacevar.hiVect()),
     tensorfab.dataPtr(),
     ARLIM(tensorfab.loVect()),ARLIM(tensorfab.hiVect()),
     SNEWfab.dataPtr(STATECOMP_VEL+force_dir), //force_dir=0 ... sdim-1
     ARLIM(SNEWfab.loVect()),ARLIM(SNEWfab.hiVect()), 
     xvel_new.dataPtr(), 
     ARLIM(xvel_new.loVect()),ARLIM(xvel_new.hiVect()), 
     xvel_old.dataPtr(), 
     ARLIM(xvel_old.loVect()),ARLIM(xvel_old.hiVect()), 
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,
     &level,
     &finest_level,
     &NS_geometry_coord,
     domlo,domhi);

  } // mfi
} // omp
  ns_reconcile_d_num(LOOP_ELASTIC_FORCE,"CELL_GRID_ELASTIC_FORCE");
 } // force_dir = 0..sdim-1

 if (elastic_force_mac_grid==1) {

  for (dir=0;dir<AMREX_SPACEDIM;dir++) {
   MultiFab& Umac_new=get_new_data(Umac_Type+dir,slab_step+1);
   MultiFab::Copy(Umac_new,*save_mac_new[dir],0,0,1,0);
  }

 } else if (elastic_force_mac_grid==0) {
  //do nothing
 } else
  amrex::Error("elastic_force_mac_grid invalid");

 for (dir=0;dir<AMREX_SPACEDIM;dir++) {
  delete save_mac_new[dir];
  delete save_mac_old[dir];
 }

} // end subroutine CELL_GRID_ELASTIC_FORCE

// PEDGE_MF allocated in allocate_pressure_work_vars
void NavierStokes::init_divup_cell_vel_cell(
 int project_option,
 int energyflag,
 int idx_pres,    //nsolve=1
 int idx_umac) {  //nsolve=1

 std::string local_caller_string="init_divup_cell_vel_cell";

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid init_divup_cell_vel_cell");

 if ((project_option==SOLVETYPE_PRES)||
     (project_option==SOLVETYPE_INITPROJ)) {  
  // do nothing
 } else
  amrex::Error("project_option invalid20 init_divup_cell_vel_cell");

 int num_colors=0;
 Vector<Real> blob_array;
 blob_array.resize(1);
 int blob_array_size=blob_array.size();
 
 bool use_tiling=ns_tiling;

 int nsolve=1;

 if (num_state_base!=2)
  amrex::Error("num_state_base!=2;NavierStokes::init_divup_cell_vel_cell");

 int finest_level=parent->finestLevel();

 resize_metrics(1);
 resize_maskfiner(1,MASKCOEF_MF);
 resize_mask_nbr(1);

 debug_ngrow(VOLUME_MF,0,local_caller_string);
 // maskcoef=1 if not covered by finer level.
 debug_ngrow(MASKCOEF_MF,1,local_caller_string); 
 // mask_nbr=1 at fine-fine bc.
 debug_ngrow(MASK_NBR_MF,1,local_caller_string); 

 int nparts=im_solid_map.size();
 if ((nparts<0)||(nparts>num_materials))
  amrex::Error("nparts invalid; NavierStokes::init_divup_cell_vel_cell");
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
  amrex::Error("nparts invalid; NavierStokes::init_divup_cell_vel_cell");

 resize_levelset(2,LEVELPC_MF);

 debug_ngrow(LEVELPC_MF,2,local_caller_string);
 if (localMF[LEVELPC_MF]->nComp()!=num_materials*(AMREX_SPACEDIM+1))
  amrex::Error(
    "localMF[LEVELPC_MF]->nComp()!=num_materials*(AMREX_SPACEDIM+1)");

 const Box& domain = geom.Domain();
 const int* domlo = domain.loVect();
 const int* domhi = domain.hiVect();

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

 const Real* dx = geom.CellSize();

 MultiFab* presmf=localMF[idx_pres];
 if (presmf->nComp()!=nsolve)
  amrex::Error("presmf->nComp()!=nsolve");

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
  amrex::Error("ncomp_check!=nsolve");

 int nstate=S_new.nComp();
 if (nstate!=STATE_NCOMP)
  amrex::Error("nstate!=STATE_NCOMP");

 if ((project_option==SOLVETYPE_PRES)||
     (project_option==SOLVETYPE_INITPROJ)) {  
  //do nothing
 } else
  amrex::Error("project_option invalid: init_divup_cell_vel_cell");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  debug_ngrow(FACE_VAR_MF+dir,0,local_caller_string);
  debug_ngrow(PEDGE_MF+dir,0,local_caller_string);
   // 0=use_face_pres  1= (2nd component) pface
  if (localMF[PEDGE_MF+dir]->nComp()!=NCOMP_PEDGE)
   amrex::Error("localMF[PEDGE_MF+dir]->nComp()!=NCOMP_PEDGE");
  if (localMF[AREA_MF+dir]->boxArray()!=
      localMF[PEDGE_MF+dir]->boxArray())
   amrex::Error("PEDGE boxarray does not match");
   // 0=use_face_pres=VALID_PEDGE
   // 1=face pressure=PRESSURE_PEDGE
   //scomp,ncomp,ngrow
   //pface=1.0e+30 initially.
  setVal_localMF(PEDGE_MF+dir,1.0e+30,PRESSURE_PEDGE,1,0);
 } // dir=0..sdim-1

 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
  if (localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0)
   amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0");
  if (localMF[FSI_GHOST_MAC_MF+data_dir]->nComp()!=nparts_def*AMREX_SPACEDIM)
   amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nComp() bad");
 }

  // old cell velocity before application of pressure gradient.
 MultiFab* ustar;

 MultiFab* divup;
 if ((energyflag==SUB_OP_THERMAL_DIVUP_NULL)||//do not update the temperature
     (energyflag==SUB_OP_THERMAL_DIVUP_OK)) {//update the temperature(comp)
  ustar=getState(1,STATECOMP_VEL,STATE_NCOMP_VEL,cur_time_slab);
  divup=new MultiFab(grids,dmap,nsolve,0,
   MFInfo().SetTag("divup"),FArrayBoxFactory());
 } else
  amrex::Error("energyflag invalid: init_divup_cell_vel_cell");

 //interpolate pressure from cell to MAC grid.
 int operation_flag_interp_pres=OP_PRES_CELL_TO_MAC; 
 int spectral_loop=0;
 int tileloop=0;

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
 
  // interpolate pressure to MAC grid for div(up) term.
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

   FArrayBox& xp=(*localMF[PEDGE_MF+dir])[mfi];
   if (xp.nComp()==NCOMP_PEDGE) {
    //do nothing
   } else
    amrex::Error("xp.nComp() invalid");

   // mask=1.0 at interior fine bc ghost cells
   FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi];
   // maskcoef=1 if not covered by finer level or outside domain
   FArrayBox& maskcoef=(*localMF[MASKCOEF_MF])[mfi]; 

   FArrayBox& maskSEMfab=(*localMF[MASKSEM_MF])[mfi];
   FArrayBox& presfab=(*presmf)[mfi];
 
   FArrayBox* solfab=&(*localMF[FSI_GHOST_MAC_MF+dir])[mfi];

   FArrayBox& levelpcfab=(*localMF[LEVELPC_MF])[mfi];

   int ncfluxreg=AMREX_SPACEDIM; //placeholder

   Vector<int> presbc;
   getBCArray_list(presbc,state_index,gridno,scomp,ncomp);
   if (presbc.size()!=nsolve*AMREX_SPACEDIM*2)
    amrex::Error("presbc.size() invalid");
   Vector<int> velbc=getBCArray(State_Type,gridno,
     STATECOMP_VEL,STATE_NCOMP_VEL);

   Real beta=0.0;

   int local_energyflag=SUB_OP_DEFAULT;
   int local_enable_spectral=0;
   int simple_AMR_BC_flag=0;
   int ncomp_xp=NCOMP_PEDGE;  //0=VALID_PEDGE  1=PRESSURE_PEDGE
   int ncomp_xgp=1;
   int ncomp_mgoni=presfab.nComp();

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();
  
   int ncphys_proxy=FACECOMP_NCOMP;

   // present routine: init_divup_cell_vel_cell; p^CELL -> p^MAC 
   fort_cell_to_mac(
    &ncomp_mgoni, 
    &ncomp_xp, // =2=NCOMP_PEDGE
    &ncomp_xgp, 
    &simple_AMR_BC_flag,
    &nsolve,
    &tileloop,
    &dir,
    &operation_flag_interp_pres, //OP_PRES_CELL_TO_MAC
    &local_energyflag,
    &beta,
    &visc_coef,
    &local_enable_spectral,
    &ncphys_proxy,
    constant_density_all_time.dataPtr(),
    presbc.dataPtr(),
    velbc.dataPtr(),
    &slab_step,
    &dt_slab,
    &cur_time_slab,
    xlo,dx,
    &spectral_loop,
    &ncfluxreg, //=AMREX_SPACEDIM (placeholder)
    levelpcfab.dataPtr(), //semflux placeholder
    ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),
    maskfab.dataPtr(), // mask=1.0 at interior fine bc ghost cells
    ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
    maskcoef.dataPtr(), // maskcoef=1 if not covered by finer level or outside
    ARLIM(maskcoef.loVect()),ARLIM(maskcoef.hiVect()),
    maskSEMfab.dataPtr(),
    ARLIM(maskSEMfab.loVect()),ARLIM(maskSEMfab.hiVect()),
    levelpcfab.dataPtr(),
    ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),
    solfab->dataPtr(),
    ARLIM(solfab->loVect()),ARLIM(solfab->hiVect()),
    xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()), //xcut
    xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()), 
    xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()),//xgp 
    xp.dataPtr(), //xp(1)="use_face_pres"=1 ok, xp(2)=face pressure
    ARLIM(xp.loVect()),ARLIM(xp.hiVect()), 
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
    &NS_geometry_coord,
    domlo,domhi,
    &nparts,
    &nparts_def,
    im_solid_map_ptr,
    blob_array.dataPtr(),
    &blob_array_size,
    &num_colors,
    &project_option);
  } // mfi
} // omp
  ns_reconcile_d_num(LOOP_PRES_CELL_TO_MAC,"init_divup_cell_vel_cell");
 } // dir=0..sdim-1

   // 0=use_face_pres=VALID_PEDGE (1==ok!)
   // 1=face pressure=PRESSURE_PEDGE
 if (nsolve!=1)
  amrex::Error("nsolve!=1 NavierStokes::init_divup_cell_vel_cell");
 int ncomp_edge_avgdown=1;
 int spectral_override=LOW_ORDER_AVGDOWN; // always low order.
 avgDownEdge_localMF(PEDGE_MF,PRESSURE_PEDGE,ncomp_edge_avgdown,
     0,AMREX_SPACEDIM,spectral_override,
     local_caller_string);

  // isweep=1 calculate cell velocity from mass weighted average of face
  //          velocity.
  // isweep=2 update cell velocity and
  //          update temperature (if compressible material and 
  //          energyflag=SUB_OP_THERMAL_DIVUP_OK)
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

   FArrayBox* solxfab=nullptr;
   FArrayBox* solyfab=nullptr;
   FArrayBox* solzfab=nullptr;

   if ((project_option==SOLVETYPE_PRES)||
       (project_option==SOLVETYPE_INITPROJ)) {
    solxfab=&(*localMF[FSI_GHOST_MAC_MF])[mfi];
    solyfab=&(*localMF[FSI_GHOST_MAC_MF+1])[mfi];
    solzfab=&(*localMF[FSI_GHOST_MAC_MF+AMREX_SPACEDIM-1])[mfi];
   } else {
    amrex::Error("project_option invalid 2169");
   }

   FArrayBox& levelpcfab=(*localMF[LEVELPC_MF])[mfi];

   FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi];//1=fine/fine bc
   FArrayBox& maskcoef=(*localMF[MASKCOEF_MF])[mfi];// 1=not covered.
   FArrayBox& maskSEMfab=(*localMF[MASKSEM_MF])[mfi];
   FArrayBox& presfab=(*presmf)[mfi];

   FArrayBox& xp=(*localMF[PEDGE_MF])[mfi];
   FArrayBox& yp=(*localMF[PEDGE_MF+1])[mfi];
   FArrayBox& zp=(*localMF[PEDGE_MF+AMREX_SPACEDIM-1])[mfi];

   if ((xp.nComp()==NCOMP_PEDGE)&&
       (yp.nComp()==NCOMP_PEDGE)&&
       (zp.nComp()==NCOMP_PEDGE)) {
    //do nothing
   } else
    amrex::Error("xp,yp, or zp invalid nComp()");

   FArrayBox& Snewfab=S_new[mfi]; // veldest
   FArrayBox& ustarfab=(*ustar)[mfi];
   FArrayBox& divupfab=(*divup)[mfi];

   Vector<int> presbc;
   getBCArray_list(presbc,state_index,gridno,scomp,ncomp);
   if (presbc.size()!=nsolve*AMREX_SPACEDIM*2)
    amrex::Error("presbc.size() invalid");
   Vector<int> velbc=getBCArray(State_Type,gridno,
     STATECOMP_VEL,STATE_NCOMP_VEL);

   int local_enable_spectral=enable_spectral;
   int operation_flag_interp_macvel=-1;
   if (isweep==1) {
    operation_flag_interp_macvel=OP_VEL_MAC_TO_CELL; 
    local_enable_spectral=enable_spectral;
   } else if (isweep==2) {
    operation_flag_interp_macvel=OP_VEL_DIVUP_TO_CELL;
    local_enable_spectral=0;
   } else {
    operation_flag_interp_macvel=-1; 
    amrex::Error("isweep, operation_flag_interp_macvel invalid1");
   }

   int homflag=0; // default
  
   // in init_divup_cell_vel_cell 

   int ncomp_denold=presfab.nComp();
   int ncomp_veldest=Snewfab.nComp();
   int ncomp_dendest=Snewfab.nComp()-STATECOMP_STATES;

   if (operation_flag_interp_macvel==OP_VEL_MAC_TO_CELL) {
    //do nothing
   } else if (operation_flag_interp_macvel==OP_VEL_DIVUP_TO_CELL) { 
    //do nothing
   } else
    amrex::Error("operation_flag_interp_macvel invalid");

   int ncphys_proxy=FACECOMP_NCOMP;

   fort_mac_to_cell(
    &ns_time_order,
    &divu_outer_sweeps,
    &num_divu_outer_sweeps,
    // OP_VEL_MAC_TO_CELL (mac_vel->cell_vel) or 
    // OP_VEL_DIVUP_TO_CELL ( div(up) low order only)
    &operation_flag_interp_macvel, 
    &energyflag,
    constant_density_all_time.dataPtr(),
    &nparts,
    &nparts_def,
    im_solid_map_ptr,
    &level, 
    &finest_level,
    &project_option,
    &local_enable_spectral,
    &ncphys_proxy,
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
    Snewfab.dataPtr(STATECOMP_STATES),
    ARLIM(Snewfab.loVect()),ARLIM(Snewfab.hiVect()), // dendest
    maskfab.dataPtr(), // 1=fine/fine  0=coarse/fine
    ARLIM(maskfab.loVect()),ARLIM(maskfab.hiVect()),
    maskcoef.dataPtr(), // 1=not covered  0=covered
    ARLIM(maskcoef.loVect()),ARLIM(maskcoef.hiVect()),
    maskSEMfab.dataPtr(), 
    ARLIM(maskSEMfab.loVect()),ARLIM(maskSEMfab.hiVect()),
    levelpcfab.dataPtr(), //levelPC
    ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),
    solxfab->dataPtr(),ARLIM(solxfab->loVect()),ARLIM(solxfab->hiVect()),
    solyfab->dataPtr(),ARLIM(solyfab->loVect()),ARLIM(solyfab->hiVect()),
    solzfab->dataPtr(),ARLIM(solzfab->loVect()),ARLIM(solzfab->hiVect()),
    levelpcfab.dataPtr(),
    ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),//cterm
    presfab.dataPtr(), 
    ARLIM(presfab.loVect()),ARLIM(presfab.hiVect()),//pold
    presfab.dataPtr(),
    ARLIM(presfab.loVect()),ARLIM(presfab.hiVect()),//denold
    ustarfab.dataPtr(),
    ARLIM(ustarfab.loVect()),ARLIM(ustarfab.hiVect()),
    levelpcfab.dataPtr(),
    ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),//mdot
    levelpcfab.dataPtr(),
    ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),//maskdivres
    levelpcfab.dataPtr(),
    ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),//maskres
    &SDC_outer_sweeps,
    &homflag,
    &nsolve,
    &ncomp_denold,
    &ncomp_veldest,
    &ncomp_dendest);

  }   // mfi
} // omp
  ns_reconcile_d_num(LOOP_VEL_MAC_TO_CELL,"init_divup_cell_vel_cell");

 }  // isweep=1,2

 if ((energyflag==SUB_OP_THERMAL_DIVUP_NULL)||
     (energyflag==SUB_OP_THERMAL_DIVUP_OK)) {
  save_to_macvel_state(idx_umac);
  delete divup; // div(up) is discarded.
  delete ustar;
 } else
  amrex::Error("energyflag invalid");

} // end subroutine init_divup_cell_vel_cell

void NavierStokes::save_to_macvel_state(int idx_umac) {

 std::string local_caller_string="save_to_macvel_state";

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  debug_ngrow(idx_umac+dir,0,local_caller_string);
  if (localMF[idx_umac+dir]->nComp()==1) {
   MultiFab& Umac_new=get_new_data(Umac_Type+dir,slab_step+1);
   MultiFab::Copy(Umac_new,*localMF[idx_umac+dir],0,0,1,0);
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

// im1=1,...,num_materials
// im2=1,...,num_materials
// iten=1,...,num_interfaces
void NavierStokes::get_iten_cpp(int im1,int im2,int& iten) {

 int im=-1;
 int im_opp=-1;

 if ((im1<1)||(im1>num_materials)||
     (im2<1)||(im2>num_materials)||
     (im1==im2)) {
  std::cout << "im1,im2 mismatch(C) im1,im2=" << im1 << ' ' << im2 << '\n';
  std::cout << "num_materials=" << num_materials << '\n';
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
  iten=num_materials-1+im_opp-2;
 } else if (im==3) {
  iten=2*num_materials-3+im_opp-3;
 } else {
  amrex::Error("im1 or im2 not supported yet");
 }

}  // get_iten_cpp


void NavierStokes::get_inverse_iten_cpp(int& im1,int& im2,int iten) {

 if (iten<1) {
  std::cout << "iten= " << iten << '\n';
  amrex::Error("iten invalid in get_inverse_iten_cpp");
 }
 im1=0;
 im2=0;
 for (int im=1;im<=num_materials;im++) {
  for (int im_opp=im+1;im_opp<=num_materials;im_opp++) {
   int iten_test;
   get_iten_cpp(im,im_opp,iten_test); 
   if (iten==iten_test) {
    im1=im;
    im2=im_opp;
   }
  }
 }
 if ((im1<1)||(im1>num_materials)||
     (im2<1)||(im2>num_materials)||
     (im1==im2)) {
  std::cout << "im1,im2 mismatch(D) im1,im2=" << im1 << ' ' << im2 << '\n';
  std::cout << "num_materials=" << num_materials << '\n';
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
 int spectral_override=SPECTRAL_ORDER_AVGDOWN;

  // avgDown all the MAC components.
  // Umac_Type
 if (level<finest_level)
  avgDownMacState(spectral_override); 

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   // ngrow,dir,time
   // Umac_Type
  MultiFab* tempmac=getStateMAC(0,dir,cur_time_slab);
  MultiFab& Umac_new=get_new_data(Umac_Type+dir,slab_step+1);
   // scomp,dcomp,ncomp,ngrow
  MultiFab::Copy(Umac_new,*tempmac,0,0,1,0);
  delete tempmac;
 } // dir=0..sdim-1

}  // subroutine make_MAC_velocity_consistent()

// ucell_new=ucell+old + force_cell
//
//called from: post_init_state, advance_MAC_velocity,
//  do_the_advance, multiphase_project,
//  INCREMENT_REGISTERS_ALL
void NavierStokes::increment_face_velocityALL(
 int operation_flag,
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

  //  unew^{f} = 
  //   (i) unew^{f} in non-solid regions
  //   (ii) u^{f,save} + (unew^{c}-u^{c,save})^{c->f} in spectral regions
  //   (iii) usolid in solid regions
 if (operation_flag==OP_U_SEM_CELL_MAC_TO_MAC) {

  if (enable_spectral==0) {
   //do nothing
  } else if (enable_spectral==1) {
  
   //ngrow,scomp,ncomp
   //u^{c,n+1}-u^{c,n}
   minusALL(1,0,AMREX_SPACEDIM,DELTA_CELL_VEL_MF,ADVECT_REGISTER_MF);

  } else
   amrex::Error("enable_spectral invalid");

 }

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.increment_face_velocity(
    operation_flag,
    project_option,
    idx_velcell,
    beta,
    blobdata); 
   // avgDownMacState, getStateMAC to fill EXT_DIR BC.
  ns_level.make_MAC_velocity_consistent();
  ParallelDescriptor::Barrier();
 }  // ilev=finest_level ... level

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
  delete_array(AMRSYNC_VEL_MF+dir);

 delete_array(DELTA_CELL_VEL_MF);
 delete_array(CURRENT_CELL_VEL_MF);

} // end subroutine increment_face_velocityALL

// OP_UNEW_CELL_TO_MAC unew^{f} = unew^{c->f} 
// OP_UNEW_USOL_MAC_TO_MAC unew^{f} = unew^{f} in fluid  (=usolid in solid)
// OP_UMAC_PLUS_VISC_CELL_TO_MAC 
//   unew^{f} = unew^{f} + beta * diffuse_register^{c->f}
// OP_U_SEM_CELL_MAC_TO_MAC unew^{f} = 
//   (i) unew^{f} in non-solid regions
//   (ii) u^{f,save} + (unew^{c}-u^{c,save})^{c->f} in spectral regions 
//   (iii) usolid in solid regions
// called from: post_init_state, do_the_advance, multiphase_project
// (when project_option==SOLVETYPE_PRES,SOLVETYPE_INITPROJ), 
// APPLY_REGISTERS, INCREMENT_REGISTERS
// called from NavierStokes::increment_face_velocityALL
void NavierStokes::increment_face_velocity(
 int operation_flag,
 int project_option,
 int idx_velcell,
 Real beta,
 Vector<blobclass> blobdata) {

 std::string local_caller_string="increment_face_velocity";

 int finest_level = parent->finestLevel();

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

 int primary_vel_data=-1;
 int secondary_vel_data=-1;

 if (operation_flag==OP_UNEW_CELL_TO_MAC) { // unew^{f} = unew^{c->f}

  if (idx_velcell==-1) {
   primary_vel_data=CURRENT_CELL_VEL_MF; 
   secondary_vel_data=CURRENT_CELL_VEL_MF; 
  } else
   amrex::Error("idx_velcell invalid increment_face_velocity");

  if (project_option==SOLVETYPE_INITPROJ) {
   // do nothing
  } else
   amrex::Error("project_option invalid21 increment_face_velocity");

  if (beta==0.0) {
   // do nothing
  } else
   amrex::Error("beta invalid");

 } else if (operation_flag==OP_UNEW_USOL_MAC_TO_MAC) { // unew^{f}=unew^{f}

  if (idx_velcell==-1) {
   primary_vel_data=CURRENT_CELL_VEL_MF; 
   secondary_vel_data=CURRENT_CELL_VEL_MF; 
  } else
   amrex::Error("idx_velcell invalid");

  if (project_option==SOLVETYPE_PRES) {
   // do nothing
  } else
   amrex::Error("project_option invalid22 increment_face_velocity");
  
  if (beta==0.0) {
   //do nothing
  } else
   amrex::Error("beta invalid");

  //unew^{f}=unew^{f}+beta*diffuse_register^{c->f}
 } else if (operation_flag==OP_UMAC_PLUS_VISC_CELL_TO_MAC) {

  if (idx_velcell>=0) {

   if ((idx_velcell==CURRENT_CELL_VEL_MF)||
       (idx_velcell==DELTA_CELL_VEL_MF))
    amrex::Error("idx_velcell collision");

   if (idx_velcell==REGISTER_CURRENT_MF) {
    // do nothing
   } else
    amrex::Error("expecting idx_velcell==REGISTER_CURRENT_MF");

   primary_vel_data=idx_velcell;  // " ucell^{n+1}-ucell^{n} "
   secondary_vel_data=CURRENT_CELL_VEL_MF;  //ucell^{n+1}

  } else
   amrex::Error("idx_velcell invalid");

  if (project_option==SOLVETYPE_VISC) {  
   //do nothing
  } else
   amrex::Error("expecting project_option==SOLVETYPE_VISC incr_face_vel");

  if ((beta==1.0)||(beta==-1.0)) {
   //do nothing
  } else
   amrex::Error("beta invalid");

  // unew^{f} = 
  //   (i) unew^{f} in non-solid regions
  //   (ii) u^{f,save} + (unew^{c}-u^{c,save})^{c->f} in spectral regions 
  //   (iii) usolid in solid regions
 } else if (operation_flag==OP_U_SEM_CELL_MAC_TO_MAC) {

  if (enable_spectral==0) {
   //do nothing
  } else if (enable_spectral==1) {

    //u^{c,n}
   debug_ngrow(ADVECT_REGISTER_MF,1,local_caller_string);

    //u^{f,n}
   for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
    debug_ngrow(ADVECT_REGISTER_FACE_MF+dir,1,local_caller_string);

  } else
   amrex::Error("enable_spectral invalid");

  if ((project_option==SOLVETYPE_PRES)|| //is_zalesak()==FALSE
      (project_option==SOLVETYPE_INITPROJ)) {//is_zalesak()==TRUE
   // do nothing
  } else
   amrex::Error("expecting project_option==SOLVETYPE_(PRES|INITPROJ)");

  if (idx_velcell==-1) {
   primary_vel_data=DELTA_CELL_VEL_MF; 
   secondary_vel_data=CURRENT_CELL_VEL_MF; 
  } else
   amrex::Error("idx_velcell invalid");
   
 } else
  amrex::Error("operation_flag invalid: increment_face_velocity");

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

 int nsolve=1;

 bool use_tiling=ns_tiling;

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 resize_levelset(2,LEVELPC_MF);
 debug_ngrow(LEVELPC_MF,2,local_caller_string);
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  debug_ngrow(FACE_VAR_MF+dir,0,local_caller_string);
 }

 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
  debug_ngrow(FSI_GHOST_MAC_MF+data_dir,0,local_caller_string);
 }

 resize_maskfiner(1,MASKCOEF_MF);
 resize_mask_nbr(1);

 const Real* dx = geom.CellSize();

 const Box& domain = geom.Domain();
 const int* domlo = domain.loVect();
 const int* domhi = domain.hiVect();

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
  debug_ngrow(FACE_VAR_MF+dir,0,local_caller_string);

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
  amrex::Error("leveltype->nGrow()<1");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  new_localMF(AMRSYNC_VEL_MF+dir,NCOMP_AMRSYNC_VEL_MF,0,dir);
   //val,scomp,ncomp,ngrow
  setVal_localMF(AMRSYNC_VEL_MF+dir,1.0e+30,0,NCOMP_AMRSYNC_VEL_MF,0);
  if (localMF[AREA_MF+dir]->boxArray()!=
      localMF[AMRSYNC_VEL_MF+dir]->boxArray())
   amrex::Error("AMRSYNC_VEL boxarray does not match");
 }

 if ((operation_flag==OP_UNEW_CELL_TO_MAC)||
     (operation_flag==OP_UMAC_PLUS_VISC_CELL_TO_MAC)||
     (operation_flag==OP_U_SEM_CELL_MAC_TO_MAC)) {

   //primary_vel_data=DELTA_CELL_VEL_MF if OP_U_SEM_CELL_MAC_TO_MAC
   //secondary_vel_data=CURRENT_CELL_VEL_MF if OP_U_SEM_CELL_MAC_TO_MAC

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
 } else if (operation_flag==OP_UNEW_USOL_MAC_TO_MAC) {
  // do nothing
 } else
  amrex::Error("operation_flag invalid");

 for (int spectral_loop=0;spectral_loop<end_spectral_loop();spectral_loop++) {

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

   MultiFab& Umac_new=get_new_data(Umac_Type+dir,slab_step+1);

   MultiFab* Umac_old;
   MultiFab* U_old;

   if (operation_flag==OP_U_SEM_CELL_MAC_TO_MAC) { 

    if (enable_spectral==0) {

      //u^{mac,n+1}
     Umac_old=&Umac_new;
      //u^{cell,n+1}
     U_old=localMF[CURRENT_CELL_VEL_MF];

    } else if (enable_spectral==1) {

      //u^{mac,n}
     Umac_old=localMF[ADVECT_REGISTER_FACE_MF+dir];
      //u^{cell,n}
     U_old=localMF[ADVECT_REGISTER_MF];

    } else
     amrex::Error("enable_spectral invalid");

   } else if ((operation_flag==OP_UNEW_CELL_TO_MAC)|| 
              (operation_flag==OP_UNEW_USOL_MAC_TO_MAC)|| 
              (operation_flag==OP_UMAC_PLUS_VISC_CELL_TO_MAC)) {

     //Umac_Type
    Umac_old=getStateMAC(0,dir,cur_time_slab); 
    if (Umac_old->boxArray()==Umac_new.boxArray()) {
     // do nothing
    } else
     amrex::Error("Umac_old->boxArray() invalid");

     //if operation_flag==OP_UMAC_PLUS_VISC_CELL_TO_MAC then
     //this variable is actually the latest velocity u^{n+1}
    U_old=localMF[CURRENT_CELL_VEL_MF];
   } else
    amrex::Error("operation_flag invalid");

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

      FArrayBox& xvel=Umac_new[mfi];

      FArrayBox& xgp=(*Umac_old)[mfi];

      FArrayBox& xp=(*localMF[AMRSYNC_VEL_MF+dir])[mfi];
      if (xp.nComp()==NCOMP_AMRSYNC_VEL_MF) {
       //do nothing
      } else
       amrex::Error("xp.nComp() invalid");

      FArrayBox& pres=(*U_old)[mfi];

       // FSI_GHOST_MAC_MF is initialized in 
       //  init_FSI_GHOST_MAC_MF_ALL(caller_string)
      FArrayBox& solfab=(*localMF[FSI_GHOST_MAC_MF+dir])[mfi];

      //primary_vel_data=DELTA_CELL_VEL_MF if OP_U_SEM_CELL_MAC_TO_MAC
      //secondary_vel_data=CURRENT_CELL_VEL_MF if OP_U_SEM_CELL_MAC_TO_MAC

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

      Vector<int> velbc=getBCArray(State_Type,gridno,
        STATECOMP_VEL,STATE_NCOMP_VEL);

      int energyflag=SUB_OP_DEFAULT;
      int local_enable_spectral=enable_spectral;
      int simple_AMR_BC_flag=0;
      int ncomp_xp=NCOMP_AMRSYNC_VEL_MF;
      int ncomp_xgp=1;
      int ncomp_mgoni=AMREX_SPACEDIM;

      int tid_current=ns_thread();
      if ((tid_current<0)||(tid_current>=thread_class::nthreads))
       amrex::Error("tid_current invalid");
      thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

      int ncphys_proxy=FACECOMP_NCOMP;

      // in increment_face_velocity
      // fort_cell_to_mac is declared in: LEVELSET_3D.F90
      fort_cell_to_mac(
       &ncomp_mgoni,
       &ncomp_xp, //=NCOMP_AMRSYNC_VEL_MF
       &ncomp_xgp,
       &simple_AMR_BC_flag,
       &nsolve,
       &tileloop,
       &dir,
        //OP_UNEW_CELL_TO_MAC,OP_UNEW_USOL_MAC_TO_MAC,
	//OP_UMAC_PLUS_VISC_CELL_TO_MAC, or
	//OP_U_SEM_CELL_MAC_TO_MAC
       &operation_flag, 
       &energyflag,
       &beta,
       &visc_coef,
       &local_enable_spectral,
       &ncphys_proxy,
       constant_density_all_time.dataPtr(),
       velbc.dataPtr(),  // presbc
       velbc.dataPtr(),  
       &slab_step,
       &dt_slab, //unused for OP_UNEW_...,OP_UMAC..,OP_U_COMP 
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
       xgp.dataPtr(),
       ARLIM(xgp.loVect()),ARLIM(xgp.hiVect()), //holds Umac_old
       xp.dataPtr(),ARLIM(xp.loVect()),ARLIM(xp.hiVect()), //xp(holds AMRSYNC)
       xvel.dataPtr(), //Umac_new
       ARLIM(xvel.loVect()),ARLIM(xvel.hiVect()), 
       primary_velfab.dataPtr(), // vel
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
       &NS_geometry_coord,
       domlo,domhi, 
       &nparts,
       &nparts_def,
       im_solid_map_ptr,
       blob_array.dataPtr(),
       &blob_array_size,
       &num_colors,
       &project_option);
    } // mfi
} // omp
    ns_reconcile_d_num(LOOP_VEL_CELL_TO_MAC,"increment_face_velocity");
   } // tileloop

   if (operation_flag==OP_U_SEM_CELL_MAC_TO_MAC) { 
    // do nothing
   } else if ((operation_flag==OP_UNEW_CELL_TO_MAC)|| 
              (operation_flag==OP_UNEW_USOL_MAC_TO_MAC)|| 
              (operation_flag==OP_UMAC_PLUS_VISC_CELL_TO_MAC)) {
    delete Umac_old;
   } else
    amrex::Error("operation_flag invalid");

  } // dir=0..sdim-1

  synchronize_flux_register(operation_flag,spectral_loop);
 } // spectral_loop

} // end subroutine increment_face_velocity

void NavierStokes::project_to_rigid_velocityALL() {

 int finest_level=parent->finestLevel();

 Vector<blobclass> blobdata;
 Vector< Vector<Real> > mdot_data;
 Vector< Vector<Real> > mdot_comp_data;
 Vector< Vector<Real> > mdot_data_redistribute;
 Vector< Vector<Real> > mdot_comp_data_redistribute;
 Vector<int> type_flag;

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "BEGIN: color_variable, project_to_rigid_velocityALL\n";
  }
 }

 int color_count=0;
 int coarsest_level=0;

 int idx_mdot=-1; //idx_mdot==-1 => do not collect auxiliary data.

 int tessellate=1;
 int operation_flag=OP_GATHER_MDOT;
  //calling from: void NavierStokes::project_to_rigid_velocityALL() 
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
   std::cout << "END: color_variable, project_to_rigid_velocityALL\n";
  }
 }

 make_MAC_velocity_consistentALL();

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.project_to_rigid_velocity(blobdata); 
   // avgDownMacState, getStateMAC to fill EXT_DIR BC.
  ns_level.make_MAC_velocity_consistent();
  ParallelDescriptor::Barrier();
 }  // ilev=finest_level ... level

 delete_array(TYPE_MF);
 delete_array(COLOR_MF);

 // spectral_override==1 => order derived from "enable_spectral"
 avgDownALL(State_Type,STATECOMP_VEL,STATE_NCOMP_VEL,SPECTRAL_ORDER_AVGDOWN);

} // end subroutine project_to_rigid_velocityALL

void NavierStokes::project_to_rigid_velocity(Vector<blobclass> blobdata) {

 std::string local_caller_string="project_to_rigid_velocity";

 int finest_level = parent->finestLevel();

 int num_colors=blobdata.size();
 if (num_colors>0) {
  //do nothing
 } else
  amrex::Error("num_colors invalid");

 int blob_array_size=num_colors*num_elements_blobclass;

 Vector<Real> blob_array(blob_array_size);

 int counter=0;

 for (int i=0;i<num_colors;i++) {
  copy_from_blobdata(i,counter,blob_array,blobdata);
 }  // i=0..num_colors-1
 if (counter!=blob_array_size)
  amrex::Error("counter invalid");

 bool use_tiling=ns_tiling;

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 resize_levelset(2,LEVELPC_MF);
 debug_ngrow(LEVELPC_MF,2,local_caller_string);

 resize_maskfiner(1,MASKCOEF_MF);
 resize_mask_nbr(1);

 const Real* dx = geom.CellSize();

 const Box& domain = geom.Domain();
 const int* domlo = domain.loVect();
 const int* domhi = domain.hiVect();

 if (localMF[COLOR_MF]->nGrow()<1)
  amrex::Error("localMF[COLOR_MF]->nGrow()<1");
 if (localMF[TYPE_MF]->nGrow()<1)
  amrex::Error("localMF[TYPE_MF]->nGrow()<1");

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

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

   FArrayBox& colorfab=(*localMF[COLOR_MF])[mfi];
   FArrayBox& typefab=(*localMF[TYPE_MF])[mfi];

   // mask=tag if not covered by level+1 or outside the domain.
   FArrayBox& maskcoeffab=(*localMF[MASKCOEF_MF])[mfi];

   Vector<int> velbc=getBCArray(State_Type,gridno,
     STATECOMP_VEL,STATE_NCOMP_VEL);

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // fort_project_to_rigid_velocity is declared in: LEVELSET_3D.F90
   fort_project_to_rigid_velocity(
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
     colorfab.dataPtr(),
     ARLIM(colorfab.loVect()),ARLIM(colorfab.hiVect()),
     typefab.dataPtr(),
     ARLIM(typefab.loVect()),ARLIM(typefab.hiVect()),
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,
     &level,&finest_level,
     &NS_geometry_coord,
     domlo,domhi, 
     blob_array.dataPtr(),
     &blob_array_size,
     &num_colors);
  } // mfi
} // omp
  ns_reconcile_d_num(LOOP_PROJECT_TO_RIGID_VEL,"project_to_rigid_velocity");

 } // dir=0..sdim-1

} // subroutine project_to_rigid_velocity


// dest_idx==-1 => destination is the state data.
// dest_idx>=0  => destination is localMF[dest_idx]
void NavierStokes::VELMAC_TO_CELLALL(int dest_idx) {

 int finest_level=parent->finestLevel();

 if (level==0) {
  // do nothing
 } else
  amrex::Error("expecting level==0 in VELMAC_TO_CELLALL");

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.VELMAC_TO_CELL(dest_idx);
 }

 if (dest_idx==-1) {
  // do nothing (update State_Type)
 } else if (dest_idx>=0) {
  Vector<int> scompBC_map;
  scompBC_map.resize(AMREX_SPACEDIM);

  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   scompBC_map[dir]=STATECOMP_VEL+dir;
  }
   //scomp=0
  GetStateFromLocalALL(dest_idx,localMF[dest_idx]->nGrow(),0,
    STATE_NCOMP_VEL,State_Type,scompBC_map);
	  
 } else
  amrex::Error("dest_idx invalid");

} // end subroutine VELMAC_TO_CELLALL

// dest_idx==-1 => destination is the state data.
// dest_idx>=0  => destination is localMF[dest_idx]
void NavierStokes::VELMAC_TO_CELL(int dest_idx) {
 
 std::string local_caller_string="VELMAC_TO_CELL";

 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid");

 int nsolve=1;

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
  debug_ngrow(FACE_VAR_MF+dir,0,local_caller_string);

 debug_ngrow(VOLUME_MF,0,local_caller_string);
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  if (localMF[AREA_MF+dir]->boxArray()!=
      localMF[FACE_VAR_MF+dir]->boxArray())
   amrex::Error("boxarrays do not match");
 }

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
  debug_ngrow(FSI_GHOST_MAC_MF+data_dir,0,local_caller_string);
 }
 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
  if (localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0)
   amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0");
  if (localMF[FSI_GHOST_MAC_MF+data_dir]->nComp()!=nparts_def*AMREX_SPACEDIM)
   amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nComp() bad");
 }

 resize_levelset(2,LEVELPC_MF);
 debug_ngrow(LEVELPC_MF,2,local_caller_string);
 if (localMF[LEVELPC_MF]->nComp()!=num_materials*(AMREX_SPACEDIM+1))
  amrex::Error("(localMF[LEVELPC_MF]->nComp()!=num_materials*(AMREX_SPACEDIM+1))");

 const Real* dx = geom.CellSize();

 resize_maskfiner(1,MASKCOEF_MF);
 resize_mask_nbr(1);
 resize_metrics(1);
 debug_ngrow(MASKCOEF_MF,1,local_caller_string); // maskcoef=1 if not covered by finer level.
 debug_ngrow(MASK_NBR_MF,1,local_caller_string); // mask_nbr=1 at fine-fine bc.

 int operation_flag=-1; 

 int local_enable_spectral=enable_spectral;

 MultiFab* face_velocity[AMREX_SPACEDIM];
 MultiFab* save_face_velocity[AMREX_SPACEDIM];
 MultiFab* dest_velocity=nullptr;

 operation_flag=OP_VEL_MAC_TO_CELL;

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   //ngrow=0
  face_velocity[dir]=getStateMAC(0,dir,cur_time_slab);
  save_face_velocity[dir]=face_velocity[dir];
 }

 if (dest_idx==-1) {
  MultiFab& S_new=get_new_data(State_Type,slab_step+1);
  dest_velocity=&S_new;
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
  debug_ngrow(FACE_VAR_MF+dir,0,local_caller_string);

  MultiFab& Umac_new=get_new_data(Umac_Type+dir,slab_step+1);
  int ncmac=Umac_new.nComp();

  if (ncmac!=nsolve) {
   std::cout << "num_materials = " << num_materials << '\n';
   std::cout << "ncmac = " << ncmac << '\n';
   amrex::Error("ncmac incorrect MAC_to_CELL");
  }
 }

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(localMF[LEVELPC_MF]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*localMF[LEVELPC_MF],use_tiling); mfi.isValid(); ++mfi) {
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
  FArrayBox& levelpcfab=(*localMF[LEVELPC_MF])[mfi];

  FArrayBox& xvel=(*face_velocity[0])[mfi];
  FArrayBox& yvel=(*face_velocity[1])[mfi];
  FArrayBox& zvel=(*face_velocity[AMREX_SPACEDIM-1])[mfi];

  if ((xvel.nComp()==1)&&
      (yvel.nComp()==1)&&
      (zvel.nComp()==1)) {
   //do nothing
  } else
   amrex::Error("xvel,yvel, or zvel invalid nComp()");

  FArrayBox& xvel_save=(*save_face_velocity[0])[mfi];
  FArrayBox& yvel_save=(*save_face_velocity[1])[mfi];
  FArrayBox& zvel_save=(*save_face_velocity[AMREX_SPACEDIM-1])[mfi];

  if ((xvel_save.nComp()==1)&&
      (yvel_save.nComp()==1)&&
      (zvel_save.nComp()==1)) {
   //do nothing
  } else
   amrex::Error("xvel_save,yvel_save, or zvel_save invalid nComp()");

  FArrayBox& volfab=(*localMF[VOLUME_MF])[mfi];
  FArrayBox& areax=(*localMF[AREA_MF])[mfi];
  FArrayBox& areay=(*localMF[AREA_MF+1])[mfi];
  FArrayBox& areaz=(*localMF[AREA_MF+AMREX_SPACEDIM-1])[mfi];

  Vector<int> velbc=getBCArray(State_Type,gridno,
    STATECOMP_VEL,STATE_NCOMP_VEL);

  int energyflag=SUB_OP_DEFAULT;
  int project_option=SOLVETYPE_PRES;
  int homflag=0; // default

  int ncomp_denold=volfab.nComp();
  int ncomp_veldest=veldest.nComp();
  int ncomp_dendest=veldest.nComp();

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

  int ncphys_proxy=FACECOMP_NCOMP;

   // we are in VELMAC_TO_CELL
   // fort_mac_to_cell is declared in LEVELSET_3D.F90
  fort_mac_to_cell(
   &ns_time_order,
   &divu_outer_sweeps,
   &num_divu_outer_sweeps,
   // operation_flag=OP_VEL_MAC_TO_CELL
   &operation_flag, 
   &energyflag,
   constant_density_all_time.dataPtr(),
   &nparts,
   &nparts_def,
   im_solid_map_ptr,
   &level,
   &finest_level,
   &project_option, //SOLVETYPE_PRES
   &local_enable_spectral, 
   &ncphys_proxy,
   velbc.dataPtr(),
   velbc.dataPtr(), // presbc
   &cur_time_slab,
   &slab_step,
   &dt_slab, //unused OP_VEL_MAC_TO_CELL
   xlo,dx,
   tilelo,tilehi,
   fablo,fabhi,
   &bfact,
   xvel_save.dataPtr(),
   ARLIM(xvel_save.loVect()),ARLIM(xvel_save.hiVect()), // xp
   yvel_save.dataPtr(),
   ARLIM(yvel_save.loVect()),ARLIM(yvel_save.hiVect()), // yp
   zvel_save.dataPtr(),
   ARLIM(zvel_save.loVect()),ARLIM(zvel_save.hiVect()), // zp
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
   levelpcfab.dataPtr(),
   ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()), // rhs
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
   levelpcfab.dataPtr(),
   ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),//mdot
   levelpcfab.dataPtr(),
   ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),//maskdivres
   levelpcfab.dataPtr(),
   ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),//maskres
   &SDC_outer_sweeps,
   &homflag,
   &nsolve,
   &ncomp_denold,
   &ncomp_veldest,
   &ncomp_dendest);
 }   // mfi
} // omp
 ns_reconcile_d_num(LOOP_VELMAC_TO_CELL_GENERAL,"VELMAC_TO_CELL");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
  delete face_velocity[dir]; 

} // end subroutine VELMAC_TO_CELL


// do_alloc=1 => allocate variable
// do_alloc=0 => variable already allocated
void NavierStokes::init_gradu_tensorALL(
 int idx_vel, //source velocity. 
              //allocated if do_alloc==1.
              //deleted if do_alloc==1.
 int do_alloc,
 int idx_cell,
 int idx_face,
 int simple_AMR_BC_flag_viscosity) {

 std::string local_caller_string="init_gradu_tensorALL";

 int finest_level=parent->finestLevel();

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid init_gradu_tensorALL");

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);

  if (do_alloc==1) {

    // ngrow,scomp,ncomp
   ns_level.getState_localMF(idx_vel,1,STATECOMP_VEL,STATE_NCOMP_VEL,
		cur_time_slab);
      
  } else if (do_alloc==0) {

   // do nothing

  } else
   amrex::Error("do_alloc invalid");

  ns_level.debug_ngrow(idx_vel,1,local_caller_string);
  if (ns_level.localMF[idx_vel]->nComp()<STATE_NCOMP_VEL)
   amrex::Error("ns_level.localMF[idx_vel]->nComp() invalid");

//ux,vx,wx,uy,vy,wy,uz,vz,wz
  int homflag=0;
  ns_level.init_gradu_tensor(
    homflag,
    idx_vel,
    idx_cell,
    idx_face,
    simple_AMR_BC_flag_viscosity);
 } // ilev=finest_level ... level

//ux,vx,wx,uy,vy,wy,uz,vz,wz
 int irow=0;
 int icol=0;
 for (int i=0;i<AMREX_SPACEDIM_SQR;i++) {
  Vector<int> scompBC_map;
   // desc_lstGHOST.setComponent(State_Type, ...
   // "set_tensor_bc", tensor_pc_interp 
   // fort_extrapfill
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
  scompBC_map[0]=EXTRAPCOMP_ELASTIC+scomp_extrap;
   // idx,ngrow,scomp,ncomp,index,scompBC_map
  PCINTERP_fill_bordersALL(idx_cell,1,i,1,State_Type,scompBC_map);

   // 00,10,20,01,11,21,02,12,22
  icol++;
  if (icol>=AMREX_SPACEDIM) {
   irow++;
   icol=0;
  }
 } // i=0..AMREX_SPACEDIM_SQR-1

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
 int homflag,
 int idx_vel,
 int idx_cell,
 int idx_face,
 int spectral_loop,
 int itensor_iter,
 MultiFab* mask3,
 int simple_AMR_BC_flag_viscosity) {

 std::string local_caller_string="doit_gradu_tensor";

 int finest_level = parent->finestLevel();

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


 debug_ngrow(idx_vel,1,local_caller_string);
 if (localMF[idx_vel]->nComp()!=AMREX_SPACEDIM)
  amrex::Error("localMF[idx_vel] ncomp invalid");

 debug_ngrow(LSTENSOR_MF,1,local_caller_string);
 debug_ngrow(MASKSOLIDTENSOR_MF,1,local_caller_string);
 if (localMF[LSTENSOR_MF]->nComp()!=AMREX_SPACEDIM)
  amrex::Error("lstensor has invalid ncomp");
 if (localMF[MASKSOLIDTENSOR_MF]->nComp()!=AMREX_SPACEDIM)
  amrex::Error("masksolidtensor has invalid ncomp");

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

 resize_maskfiner(1,MASKCOEF_MF);
 resize_mask_nbr(1);
 resize_metrics(1);
 resize_levelset(2,LEVELPC_MF);

 debug_ngrow(MASKCOEF_MF,1,local_caller_string);
 debug_ngrow(VOLUME_MF,1,local_caller_string);
 debug_ngrow(MASKSEM_MF,1,local_caller_string);
 debug_ngrow(LEVELPC_MF,2,local_caller_string);
 if (localMF[LEVELPC_MF]->nComp()!=num_materials*(1+AMREX_SPACEDIM))
  amrex::Error("levelpc mf has incorrect ncomp");
 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
  debug_ngrow(FSI_GHOST_MAC_MF+data_dir,0,local_caller_string);
 }

 const Real* dx = geom.CellSize();

 MultiFab* sem_flux_mf;

 if (itensor_iter==0) {  // compute grad U

  sem_flux_mf=localMF[SEM_FLUXREG_MF];
  if (sem_flux_mf->nComp()!=AMREX_SPACEDIM_SQR)
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

    Vector<int> velbc=getBCArray(State_Type,gridno,
      STATECOMP_VEL,STATE_NCOMP_VEL);

    FArrayBox& velfab=(*localMF[idx_vel])[mfi];

    FArrayBox& solidxfab=(*localMF[FSI_GHOST_MAC_MF])[mfi];
    FArrayBox& solidyfab=(*localMF[FSI_GHOST_MAC_MF+1])[mfi];
    FArrayBox& solidzfab=(*localMF[FSI_GHOST_MAC_MF+AMREX_SPACEDIM-1])[mfi];

    FArrayBox& levelpcfab=(*localMF[LEVELPC_MF])[mfi];

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

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // fort_face_gradients is declared in GODUNOV_3D.F90
    fort_face_gradients(
     &ns_time_order,
     &divu_outer_sweeps,
     &num_divu_outer_sweeps,
     &SDC_outer_sweeps,
     &tileloop,
     &dir,
     &slab_step,
     &itensor_iter,
     &cur_time_slab,
     &enable_spectral,
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
     xlo,dx,
     &NS_geometry_coord,
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,&bfact_c,&bfact_f,
     &level,
     &finest_level,
     &nparts,
     &nparts_def,
     im_solid_map_ptr,
     &homflag,
     &simple_AMR_BC_flag_viscosity);
   } // mfi
} // omp
   ns_reconcile_d_num(LOOP_FACE_GRADIENTS,"doit_gradu_tensor");
  } // tileloop=0..1
 } // dir

 int datatype=0;

 int im=0;

 for (int sc=0;sc<AMREX_SPACEDIM_SQR;sc++) {

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

   int sc_mat=im*AMREX_SPACEDIM_SQR+sc;

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
  check_for_NAN(localMF[idx_vel]);

  if (itensor_iter==1) {  // cell grad U
   datatype=2;
   check_for_NAN_TENSOR(datatype,localMF[idx_cell]);
  } else if (itensor_iter==0) { // face grad U
   datatype=1;
   check_for_NAN_TENSOR(datatype,localMF[idx_face]);
   for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
    check_for_NAN_TENSOR_base(datatype,localMF[LSTENSOR_MF],dir,dir);
    check_for_NAN_TENSOR_base(datatype,localMF[MASKSOLIDTENSOR_MF],dir,dir);
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

 IndexType mactyp=TheUMACType;
 if (dir==0) {
  mactyp=TheUMACType;
 } else if (dir==1) {
  mactyp=TheVMACType;
 } else if ((dir==2)&&(AMREX_SPACEDIM==3)) {
  mactyp=TheWMACType;
 } else
  amrex::Error("dir invalid FillBoundaryTENSOR");

  //nc=1 ng=1
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

  Array4<Real> const& cellfab_array=cellfab.array();
  Array4<Real> const& macfab_array=macfab.array();
  Array4<Real> const& localfab_array=localfab.array();

  const Dim3 lo3=amrex::lbound(localbox);
  const Dim3 hi3=amrex::ubound(localbox);
  for (int z=lo3.z;z<=hi3.z;++z) {
  for (int y=lo3.y;y<=hi3.y;++y) {
  for (int x=lo3.x;x<=hi3.x;++x) {
   localfab_array(x,y,z,0)=cellfab_array(x,y,z,sc);
  }
  }
  }

  localfab.SetBoxType(mactyp);

  Box macbox(macfab.box());
  const Dim3 lo3b=amrex::lbound(macbox);
  const Dim3 hi3b=amrex::ubound(macbox);
  for (int z=lo3b.z;z<=hi3b.z;++z) {
  for (int y=lo3b.y;y<=hi3b.y;++y) {
  for (int x=lo3b.x;x<=hi3b.x;++x) {
   macfab_array(x,y,z,0)=0.0;
  }
  }
  }

  const Dim3 lo3c=amrex::lbound(localfab.box());
  const Dim3 hi3c=amrex::ubound(localfab.box());
  for (int z=lo3c.z;z<=hi3c.z;++z) {
  for (int y=lo3c.y;y<=hi3c.y;++y) {
  for (int x=lo3c.x;x<=hi3c.x;++x) {
   macfab_array(x,y,z,0)=localfab_array(x,y,z,0);
  }
  }
  }

 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_FILL_BOUNDARY_TENSOR,"FillBoundaryTENSOR");

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

  Array4<Real> const& cellfab_array=cellfab.array();
  Array4<Real> const& macfab_array=macfab.array();
  Array4<Real> const& localfab_array=localfab.array();

  localfab.SetBoxType(mactyp);

  const Dim3 lo3=amrex::lbound(localfab.box());
  const Dim3 hi3=amrex::ubound(localfab.box());

  for (int z=lo3.z;z<=hi3.z;++z) {
  for (int y=lo3.y;y<=hi3.y;++y) {
  for (int x=lo3.x;x<=hi3.x;++x) {
   localfab_array(x,y,z,0)=macfab_array(x,y,z,0);
  }
  }
  }

  localfab.SetBoxType(IndexType::TheCellType());

  const Dim3 lo3b=amrex::lbound(localfab.box());
  const Dim3 hi3b=amrex::ubound(localfab.box());

  for (int z=lo3b.z;z<=hi3b.z;++z) {
  for (int y=lo3b.y;y<=hi3b.y;++y) {
  for (int x=lo3b.x;x<=hi3b.x;++x) {
   cellfab_array(x,y,z,sc)=localfab_array(x,y,z,0);
  }
  } 
  }

 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_OVERRIDDE_FAB_TYPE,"FillBoundaryTENSOR");

 delete macmf;

} // subroutine FillBoundaryTENSOR

// called from apply_pressure_grad, init_gradu_tensorALL
void NavierStokes::init_gradu_tensor(
 int homflag,
 int idx_vel,
 int idx_cell,
 int idx_face,
 int simple_AMR_BC_flag_viscosity) {

 std::string local_caller_string="init_gradu_tensor";

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
 int operation_flag=OP_UGRAD_MAC; // evaluate tensor values

 if ((localMF_grow[idx_face]>=0)||
     (localMF_grow[idx_cell]>=0)||
     (localMF_grow[LSTENSOR_MF]>=0)||
     (localMF_grow[MASKSOLIDTENSOR_MF]>=0))
  amrex::Error("tensor scratch variables not previously deleted");

 debug_ngrow(idx_vel,1,local_caller_string);
 if (localMF[idx_vel]->nComp()!=AMREX_SPACEDIM)
  amrex::Error("localMF[idx_vel]->nComp() invalid in init_gradu_tensor");

 int nmasksolid=AMREX_SPACEDIM;

 new_localMF(idx_face,AMREX_SPACEDIM_SQR,1,-1);
 new_localMF(idx_cell,AMREX_SPACEDIM_SQR,1,-1);

  // x,y,z 
 new_localMF(MASKSOLIDTENSOR_MF,nmasksolid,1,-1);
  // x,y,z 
 new_localMF(LSTENSOR_MF,AMREX_SPACEDIM,1,-1);

 int itensor_iter=0; // tensor face (face grad U)

  // flux register is initialized to zero.
 allocate_flux_register(operation_flag);
 if (localMF[SEM_FLUXREG_MF]->nComp()!=AMREX_SPACEDIM_SQR)
  amrex::Error("localMF[SEM_FLUXREG_MF]->nComp() invalid4");

  // spectral_loop==0 
  //   low order face fluxes
  //   high order face fluxes in localMF[SEM_FLUXREG_MF]
  // spectral_loop==1
  //   average of flux values shared at face.
 int spectral_loop=0;
 for (spectral_loop=0;spectral_loop<end_spectral_loop();spectral_loop++) {
   // NavierStokes::doit_gradu_tensor is declared in NavierStokes2.cpp
  doit_gradu_tensor(
   homflag,
   idx_vel,
   idx_cell,
   idx_face,
   spectral_loop,
   itensor_iter,
   mask3,
   simple_AMR_BC_flag_viscosity);
  synchronize_flux_register(operation_flag,spectral_loop);
 }

  // interpolate grad U from MAC grid to CELL grid.
 itensor_iter=1;  // tensor cell
 spectral_loop=0;
 doit_gradu_tensor(
   homflag,
   idx_vel,
   idx_cell,
   idx_face,
   spectral_loop,
   itensor_iter,mask3,
   simple_AMR_BC_flag_viscosity);

 delete mask3; 

} // end subroutine init_gradu_tensor

  
// if projection (energyflag=SUB_OP_FOR_MAIN):
//  gp_mf = - dt*(grad p)*face_weight  
// if called from update_SEM_forces with project_option
// equal to SOLVETYPE_HEAT or SOLVETYPE_VISC, then
//  gp_mf=-k grad T or -2 mu D respectively  (dt=1)
// if called from update_SEM_forces with project_option=SOLVETYPE_PRES,
//  gp_mf=grad P (energyflag=SUB_OP_FOR_SDC)
// face_weight=0 at embedded solid faces and on 
// the domain boundary where pressure has a Neumann BC.
void NavierStokes::apply_pressure_grad(
  int simple_AMR_BC_flag,
  int simple_AMR_BC_flag_viscosity,
  int homflag,
  int energyflag,
  int gp_mf,
  int pboth_mf,
  int project_option,
  int nsolve,
  Real dt_pressure_grad) {

 std::string local_caller_string="apply_pressure_grad";

 int finest_level = parent->finestLevel();

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

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 if ((energyflag!=SUB_OP_FOR_MAIN)&&
     (energyflag!=SUB_OP_FOR_SDC))
  amrex::Error("energyflag invalid");

 debug_ngrow(pboth_mf,1,local_caller_string);

 if (localMF[pboth_mf]->nComp()!=nsolve) {
  std::cout << "nsolve=" << nsolve << '\n';
  print_project_option(project_option);
  std::cout << "pboth ngrow= " << localMF[pboth_mf]->nGrow() << '\n';
  std::cout << "pboth ncomp= " << localMF[pboth_mf]->nComp() << '\n';
  amrex::Error("nsolve invalid28");
 }

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

 resize_maskfiner(1,MASKCOEF_MF);
 resize_mask_nbr(1);
 resize_metrics(1);

 debug_ngrow(MASKCOEF_MF,1,local_caller_string);
 debug_ngrow(VOLUME_MF,1,local_caller_string);
 if (localMF[LEVELPC_MF]->nComp()!=num_materials*(1+AMREX_SPACEDIM))
  amrex::Error("levelpc mf has incorrect ncomp");

 int local_fsi_ghost_ncomp=nparts_def*AMREX_SPACEDIM;
 int local_fsi_ghost_ngrow=0;
 int local_amrsync_pres_ncomp=nsolve;
 int local_sem_fluxreg_ncomp=AMREX_SPACEDIM*nsolve;

 if (project_option_is_valid(project_option)==1) {
  // do nothing
 } else
  amrex::Error("project_option invalid 4319");

 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
  if (localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=
      local_fsi_ghost_ngrow) {
   amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow() bad");
  }
  if (localMF[FSI_GHOST_MAC_MF+data_dir]->nComp()!=local_fsi_ghost_ncomp) {
   amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nComp() bad");
  }
 }

 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
  debug_ngrow(FSI_GHOST_MAC_MF+data_dir,
              local_fsi_ghost_ngrow,local_caller_string);
 }

 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) 
  debug_ngrow(FACE_VAR_MF+data_dir,0,local_caller_string);

 const Box& domain = geom.Domain();
 const int* domlo = domain.loVect();
 const int* domhi = domain.hiVect();

 const Real* dx = geom.CellSize();

 debug_ngrow(MASKSEM_MF,1,local_caller_string);

 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {

  if (localMF[AMRSYNC_PRES_MF+data_dir]->nComp()!=
      local_amrsync_pres_ncomp)
   amrex::Error("localMF[AMRSYNC_PRES_MF+data_dir]->nComp() invalid29");
  if (localMF[AMRSYNC_PRES_MF+data_dir]->boxArray()!=
      localMF[AREA_MF+data_dir]->boxArray())
   amrex::Error("AMRSYNC_PRES_MF boxarrays do not match");

  if (localMF[gp_mf+data_dir]->nComp()!=nsolve)
   amrex::Error("localMF[gp_mf+data_dir]->nComp() invalid29");
  if (localMF[AREA_MF+data_dir]->boxArray()!=
      localMF[gp_mf+data_dir]->boxArray())
   amrex::Error("gp_mf boxarrays do not match");
 } // data_dir=0..sdim-1

 if (project_option==SOLVETYPE_VISC) {

  if (nsolve!=AMREX_SPACEDIM)
   amrex::Error("expecting nsolve==SDIM; project_option==SOLVETYPE_VISC");

  int operation_flag=OP_UGRAD_COUPLING_MAC;

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

   // NavierStokes::init_gradu_tensor is declared in NavierStokes2.cpp
  init_gradu_tensor(
    homflag,
    pboth_mf,
    LOCAL_CELLTENSOR_MF, //should not be allocated at this point
    LOCAL_FACETENSOR_MF, //should not be allocated at this point
    simple_AMR_BC_flag_viscosity);

  show_norm2(localMF[pboth_mf],0,localMF[pboth_mf]->nComp(),20);
  show_norm2(localMF[LOCAL_CELLTENSOR_MF],0,
     localMF[LOCAL_CELLTENSOR_MF]->nComp(),21);
  show_norm2(localMF[LOCAL_FACETENSOR_MF],0,
     localMF[LOCAL_FACETENSOR_MF]->nComp(),21);

  int nden=num_materials*num_state_material;

  if (local_sem_fluxreg_ncomp==AMREX_SPACEDIM_SQR) {
   // do nothing
  } else
   amrex::Error("local_sem_fluxreg_ncomp invalid");

  allocate_flux_register(operation_flag);
  if (localMF[SEM_FLUXREG_MF]->nComp()!=AMREX_SPACEDIM_SQR)
   amrex::Error("localMF[SEM_FLUXREG_MF]->nComp() invalid5");

  resize_levelset(2,LEVELPC_MF);
  debug_ngrow(LEVELPC_MF,2,local_caller_string);

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
 
    const Real* xlo = grid_loc[gridno].lo();

    Vector<int> velbc=getBCArray(State_Type,gridno,
       STATECOMP_VEL,STATE_NCOMP_VEL);

    FArrayBox& velfab=(*localMF[pboth_mf])[mfi];
    FArrayBox& levelpcfab=(*localMF[LEVELPC_MF])[mfi];
    FArrayBox& maskSEMfab=(*localMF[MASKSEM_MF])[mfi];

    FArrayBox& xflux=(*localMF[gp_mf+dir-1])[mfi];

    FArrayBox& xface=(*localMF[FACE_VAR_MF+dir-1])[mfi];

    FArrayBox& tensor_data=(*localMF[LOCAL_FACETENSOR_MF])[mfi];
    FArrayBox& cell_tensor_data=(*localMF[LOCAL_CELLTENSOR_MF])[mfi];
    FArrayBox& mask_tensor_data=(*localMF[MASKSOLIDTENSOR_MF])[mfi];
    FArrayBox& faceLS=(*localMF[LSTENSOR_MF])[mfi];

    FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi];
    // maskcoef=tag if not covered by level+1 or outside the domain.
    FArrayBox& maskcoef_fab=(*localMF[MASKCOEF_MF])[mfi];

    FArrayBox& semfluxfab=(*localMF[SEM_FLUXREG_MF])[mfi];
    int ncfluxreg=semfluxfab.nComp();

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // declared in: GODUNOV_3D.F90
    // -dt * visc_coef * viscface * (grad U + grad U^T)
    fort_crossterm(
     &nsolve,
     &tileloop,
     &dir,
     &operation_flag, // OP_UGRAD_COUPLING_MAC
     &enable_spectral,
     &spectral_loop,
     &ncfluxreg,
     BL_TO_FORTRAN_ANYD(semfluxfab),
     BL_TO_FORTRAN_ANYD(maskfab), // 1=fine/fine  0=coarse/fine
     BL_TO_FORTRAN_ANYD(maskcoef_fab),//maskcoef=tag if not cov or outside.
     BL_TO_FORTRAN_ANYD(faceLS),
     BL_TO_FORTRAN_ANYD(mask_tensor_data),
     BL_TO_FORTRAN_ANYD(tensor_data),
     BL_TO_FORTRAN_ANYD(cell_tensor_data),
     BL_TO_FORTRAN_ANYD(maskSEMfab),
     xlo,dx,
     &dt_pressure_grad,
     &cur_time_slab,
     BL_TO_FORTRAN_ANYD(velfab),
     BL_TO_FORTRAN_ANYD(levelpcfab),
     BL_TO_FORTRAN_ANYD(xflux),
     BL_TO_FORTRAN_ANYD(xface),
     BL_TO_FORTRAN_BOX(tilegrid),
     BL_TO_FORTRAN_BOX(fabgrid),
     &bfact,
     &level,
     &NS_geometry_coord,
     velbc.dataPtr(),
     &visc_coef,
     &nden,
     &uncoupled_viscosity,
     &homflag);
   } // mfi
} // omp
   ns_reconcile_d_num(LOOP_CROSSTERM,"apply_pressure_grad");

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
   check_for_NAN(localMF[pboth_mf]);
   for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
    check_for_NAN(localMF[gp_mf+data_dir]);
   }
  } else if (check_nan==0) {
   // do nothing
  } else {
   amrex::Error("check_nan invalid");
  }

  //above: SOLVETYPE_VISC
  //below: everything else.
 } else if (project_option_is_valid(project_option)==1) {

  int num_colors=0;
  Vector<Real> blob_array;
  blob_array.resize(1);
  int blob_array_size=blob_array.size();

  int operation_flag=OP_PRESGRAD_MAC;

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

  resize_levelset(2,LEVELPC_MF);

  if (project_option_is_valid(project_option)==1) {

   allocate_flux_register(operation_flag);

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

  } else
   amrex::Error("project_option invalid 4657");

  if (localMF[SEM_FLUXREG_MF]->nComp()!=local_sem_fluxreg_ncomp)
   amrex::Error("localMF[SEM_FLUXREG_MF]->nComp() invalid6");

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
    if (xp.nComp()==nsolve) {
     if (nsolve==1) {
      //do nothing
     } else
      amrex::Error("expecting nsolve==1 OP_PRESGRAD_MAC");
    } else
     amrex::Error("expecting xp.nComp()==nsolve");

    FArrayBox& xgp=(*localMF[gp_mf+dir])[mfi];
    FArrayBox& xcut=(*localMF[FACE_WEIGHT_MF+dir])[mfi]; // A/rho
    FArrayBox& xface=(*localMF[FACE_VAR_MF+dir])[mfi];

    FArrayBox& maskSEMfab=(*localMF[MASKSEM_MF])[mfi];
    FArrayBox& presfab=(*localMF[pboth_mf])[mfi]; // in: apply_pressure_grad

    FArrayBox* solfab;
    solfab=&(*localMF[FSI_GHOST_MAC_MF+dir])[mfi];

    FArrayBox& levelpcfab=(*localMF[LEVELPC_MF])[mfi];

    Vector<int> presbc;
    getBCArray_list(presbc,state_index,gridno,scomp,ncomp);
    if (presbc.size()!=nsolve*AMREX_SPACEDIM*2)
     amrex::Error("presbc.size() invalid");
    Vector<int> velbc=getBCArray(State_Type,gridno,
      STATECOMP_VEL,STATE_NCOMP_VEL);
  
    Real beta=0.0;

    // mask=tag if not covered by level+1 or outside the domain.
    FArrayBox& maskcoeffab=(*localMF[MASKCOEF_MF])[mfi];
    FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi];

    FArrayBox& semfluxfab=(*localMF[SEM_FLUXREG_MF])[mfi];
    int ncfluxreg=semfluxfab.nComp();
    if (ncfluxreg!=local_sem_fluxreg_ncomp)
     amrex::Error("ncfluxreg invalid");

    int ncomp_xp=nsolve;
    int ncomp_xgp=nsolve;
    int ncomp_mgoni=presfab.nComp();

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    int ncphys_proxy=FACECOMP_NCOMP;

    // -grad p * FACE_WEIGHT * dt
    // fort_cell_to_mac called from: apply_pressure_grad
    // fort_cell_to_mac is declared in: LEVELSET_3D.F90
    fort_cell_to_mac(
     &ncomp_mgoni,
     &ncomp_xp,
     &ncomp_xgp,
     &simple_AMR_BC_flag,
     &nsolve,
     &tileloop,
     &dir,
     &operation_flag, // OP_PRESGRAD_MAC
     &energyflag,
     &beta,
     &visc_coef,
     &enable_spectral,
     &ncphys_proxy,
     constant_density_all_time.dataPtr(),
     presbc.dataPtr(),
     velbc.dataPtr(),
     &slab_step,
     &dt_pressure_grad,
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
     solfab->dataPtr(),
     ARLIM(solfab->loVect()),ARLIM(solfab->hiVect()),
      //xcut holds FACE_WEIGHT_MF
     xcut.dataPtr(),ARLIM(xcut.loVect()),ARLIM(xcut.hiVect()),
      //xface holds FACE_VAR_MF
     xface.dataPtr(),ARLIM(xface.loVect()),ARLIM(xface.hiVect()),
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
     &NS_geometry_coord,
     domlo,domhi,
     &nparts,
     &nparts_def,
     im_solid_map_ptr,
     blob_array.dataPtr(),
     &blob_array_size,
     &num_colors,
     &project_option);

   }  // mfi
} // omp

   ns_reconcile_d_num(LOOP_PRESGRAD_MAC,"apply_pressure_grad");

  } // tileloop
  } // dir

  if (project_option_is_valid(project_option)==1) {
   synchronize_flux_register(operation_flag,spectral_loop);
  } else
   amrex::Error("project_option invalid 4827");

  } // spectral_loop

  if (project_option==SOLVETYPE_HEAT) { 
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
  amrex::Error("project_option invalid27: apply_pressure_grad");

} // end subroutine apply_pressure_grad

void NavierStokes::init_gradu_tensor_and_material_visc_ALL(
  const std::string& caller_string) {

 std::string local_caller_string="init_gradu_tensor_and_material_visc_ALL";
 local_caller_string=caller_string+local_caller_string;

 // allocate and delete HOLD_VELOCITY_DATA_MF in init_gradu_tensorALL:
 // (since do_alloc==1)
 int simple_AMR_BC_flag_viscosity=1;
 int do_alloc=1; 
 init_gradu_tensorALL(
   HOLD_VELOCITY_DATA_MF, //alloc and delete since do_alloc==1
   do_alloc,
   CELLTENSOR_MF,
   FACETENSOR_MF,
   simple_AMR_BC_flag_viscosity);

  //localMF[CELL_VISC_MATERIAL_MF] is deleted in ::Geometry_cleanup()
  //ngrow=1
  //we are in:init_gradu_tensor_and_material_visc_ALL
 getStateVISC_ALL(local_caller_string); 

} // end subroutine init_gradu_tensor_and_material_visc_ALL

// called from:
//  NavierStokes::volWgtSumALL
//  NavierStokes::prepare_post_process
//  NavierStokes::do_the_advance
//  NavierStokes::nucleation_code_segment
//
void NavierStokes::make_physics_varsALL(int project_option,
        const std::string& caller_string) {

 std::string local_caller_string="make_physics_varsALL";
 local_caller_string=caller_string+local_caller_string;

 if (level!=0)
  amrex::Error("level invalid make_physics_varsALL");

 if ((project_option==SOLVETYPE_PRES)||
     (project_option==SOLVETYPE_INITPROJ)) {
  // do nothing
 } else
  amrex::Error("project_option invalid make_physics_varsALL");

 if ((SDC_outer_sweeps>=0)&&
     (SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid");

 int nhistory=num_interfaces*2;

 int finest_level = parent->finestLevel();

 debug_ngrow(SLOPE_RECON_MF,0,local_caller_string);
 if (localMF[SLOPE_RECON_MF]->nComp()==num_materials*ngeom_recon) {
  // do nothing
 } else
  amrex::Error("localMF[SLOPE_RECON_MF]->nComp() invalid");

  // in: NavierStokes::make_physics_varsALL
  // piecewise constant interpolation.
 allocate_levelset_ALL(1,LEVELPC_MF);

   // create DIST_CURV_MF 

 curv_min.resize(thread_class::nthreads);
 curv_max.resize(thread_class::nthreads);

 for (int tid=0;tid<thread_class::nthreads;tid++) {
  curv_min[tid]=1.0e+30;
  curv_max[tid]=-1.0e+30;
 } // tid

 allocate_array(1,nhistory,-1,HISTORY_MF);

 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.makeStateCurv(project_option,local_caller_string);
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
    local_caller_string,
    HISTORY_MF, //tower_mf_id
    localMF[HISTORY_MF]->nComp(),
    HISTORY_MF,
    -1, // State_Type==-1
    -1, // data_dir==-1
    parent->levelSteps(0)); 
 }
 delete_array(HISTORY_MF);

  //DIST_CURV_MF is initialized and filled in makeStateCurv.
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.avgDownCURV_localMF(DIST_CURV_MF);
 }

  //localMF[CELL_VISC_MATERIAL_MF] is deleted in ::Geometry_cleanup()
  //responsibility of caller to issue commands,
  // delete_array(CELLTENSOR_MF);
  // delete_array(FACETENSOR_MF);
  //
 init_gradu_tensor_and_material_visc_ALL(local_caller_string);

 debug_ngrow(CELL_VISC_MATERIAL_MF,1,local_caller_string);
 int ncomp_visc=localMF[CELL_VISC_MATERIAL_MF]->nComp();
 if (ncomp_visc!=3*num_materials)
  amrex::Error("visc_data invalid ncomp");

 debug_ngrow(CELL_CONDUCTIVITY_MATERIAL_MF,1,local_caller_string);
 int ncomp_conductivity=localMF[CELL_CONDUCTIVITY_MATERIAL_MF]->nComp();
 if (ncomp_conductivity!=num_materials)
  amrex::Error("conductivity_data invalid ncomp");

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.getStateDen_localMF(DEN_RECON_MF,1,cur_time_slab);
  ns_level.getStateMOM_DEN(MOM_DEN_MF,1,cur_time_slab);
  if ((num_materials_compressible>=1)&&
      (num_materials_compressible<=num_materials)) {
   ns_level.getStateRefineDensity_localMF(
     REFINE_DENSITY_RECON_MF,
     1,0,
     NUM_CELL_REFINE_DENSITY,
     cur_time_slab);
  } else if (num_materials_compressible==0) {
   //do nothing
  } else
   amrex::Error("num_materials_compressble invalid");
 } // for (int ilev=finest_level;ilev>=level;ilev--)

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);

  ns_level.make_physics_vars(project_option,local_caller_string);
   //NavierStokes::level_init_elasticmask_and_elasticmaskpart is 
   //declared in NavierStokes.cpp
  ns_level.level_init_elasticmask_and_elasticmaskpart();

   // average down from ilev+1 to ilev.
  
    // idxMF,scomp,ncomp,start_dir,ndir
  ns_level.avgDownEdge_localMF(FACE_VAR_MF,FACECOMP_FACECUT,1,0,
    AMREX_SPACEDIM,LOW_ORDER_AVGDOWN,local_caller_string);
  ns_level.avgDownEdge_localMF(FACE_VAR_MF,FACECOMP_ELASTICMASK,1,0,
    AMREX_SPACEDIM,LOW_ORDER_AVGDOWN,local_caller_string);
  ns_level.avgDownEdge_localMF(FACE_VAR_MF,FACECOMP_ELASTICMASKPART,1,0,
    AMREX_SPACEDIM,LOW_ORDER_AVGDOWN,local_caller_string);

   // spectral_override==0 => always low order.
  ns_level.avgDownEdge_localMF(FACE_VAR_MF,FACECOMP_FACEDEN,1,0,
	  AMREX_SPACEDIM,LOW_ORDER_AVGDOWN,local_caller_string);

  ns_level.avgDownEdge_localMF(FACE_VAR_MF,FACECOMP_FACEVISC,1,0,
          AMREX_SPACEDIM,LOW_ORDER_AVGDOWN,local_caller_string);
  ns_level.avgDownEdge_localMF(FACE_VAR_MF,FACECOMP_FACEHEAT,1,0,
	  AMREX_SPACEDIM,LOW_ORDER_AVGDOWN,local_caller_string);
  if (num_species_var>0) {
   ns_level.avgDownEdge_localMF(FACE_VAR_MF,FACECOMP_FACESPEC,
      num_species_var,0,AMREX_SPACEDIM,
      LOW_ORDER_AVGDOWN,local_caller_string);
  }

 }  // ilev=finest_level ... level

 if (1==0) {
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
    // filenames: "FACE_VAR<stuff>.plt" (MAC data)
   writeSanityCheckData(
    "FACE_VAR",
    "in: make_physics_varsALL, FACE_VAR_MF",
    local_caller_string,
    FACE_VAR_MF+dir, //tower_mf_id
    localMF[FACE_VAR_MF+dir]->nComp(),
    FACE_VAR_MF+dir,
    -1, // State_Type==-1
    dir, 
    parent->levelSteps(0));
  } // dir=0..sdim-1
 }

 delete_array(DEN_RECON_MF);
 delete_array(MOM_DEN_MF);
 if ((num_materials_compressible>=1)&&
     (num_materials_compressible<=num_materials)) {
  delete_array(REFINE_DENSITY_RECON_MF);
 } else if (num_materials_compressible==0) {
  // do nothing
 } else
  amrex::Error("num_materials_compressible invalid:make_physics_varsALL");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) 
  delete_array(AMRSYNC_VEL_MF+dir);

 //it is the responsibility of the caller to issue the following commands:
 //delete_array(CELLTENSOR_MF);
 //delete_array(FACETENSOR_MF);
 
} // end subroutine make_physics_varsALL

// called from: prelim_alloc() and make_physics_vars
void NavierStokes::allocate_physics_vars() {

 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid allocate_physics_vars");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  new_localMF_if_not_exist(FACE_VAR_MF+dir,FACECOMP_NCOMP,0,dir);
  new_localMF_if_not_exist(FACE_DEN_HOLD_MF+dir,1,0,dir);
  new_localMF_if_not_exist(FACE_VISC_HOLD_MF+dir,1,0,dir);
 }

  // ncomp,ngrow,dir
 if (localMF_grow[SWEPT_CROSSING_MF]==-1) {
  new_localMF(SWEPT_CROSSING_MF,num_materials,0,-1); 
    //val,scomp,ncomp,ngrow
  setVal_localMF(SWEPT_CROSSING_MF,1.0,0,num_materials,0);
 } else if (localMF_grow[SWEPT_CROSSING_MF]>=0) {
  // do nothing
 } else
  amrex::Error("localMF_grow[SWEPT_CROSSING_MF] invalid");

  //CELL_DEDT_MF contains 1/(rho CV)
  //CELL_DEDT_MF is passed as a parameter to:
  //fort_scalarcoeff
 new_localMF_if_not_exist(CELL_DEDT_MF,1,1,-1); // ncomp,ngrow,dir

  //CELL_DEN_MF contains 1/rho
 new_localMF_if_not_exist(CELL_DEN_MF,1,1,-1); // ncomp,ngrow,dir
 new_localMF_if_not_exist(CELL_DEN_HOLD_MF,1,1,-1); // ncomp,ngrow,dir

  // coeff_avg,padvect_avg 
 new_localMF_if_not_exist(CELL_SOUND_MF,2,0,-1); // ncomp,ngrow,dir

  // tessellating volume fractions.
 new_localMF_if_not_exist(CELL_VOF_MF,num_materials,1,-1); // ncomp,ngrow,dir
 new_localMF_if_not_exist(CELL_VISC_MF,1,1,-1); // ncomp,ngrow,dir
 new_localMF_if_not_exist(CELL_VISC_HOLD_MF,1,1,-1); // ncomp,ngrow,dir

} // allocate_physics_vars

void NavierStokes::allocate_levelset_ALL(int ngrow,int idx) {

 int finest_level=parent->finestLevel();
 if (level!=0)
  amrex::Error("level invalid allocate_levelsetALL");

 if ((ngrow<0)||(ngrow>ngrow_distance))
  amrex::Error("ngrow invalid");

 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.allocate_levelset(ngrow,idx);
 }

}  // allocate_levelset_ALL

void NavierStokes::allocate_levelset(int ngrow,int idx) {

 std::string local_caller_string="allocate_levelset";

 if ((ngrow<0)||(ngrow>ngrow_distance))
  amrex::Error("ngrow invalid");

 delete_localMF_if_exist(idx,1);
 getStateDist_localMF(idx,ngrow,cur_time_slab,local_caller_string);
 if (localMF[idx]->nComp()!=num_materials*(AMREX_SPACEDIM+1))
  amrex::Error("localMF[idx]->nComp()!=num_materials*(AMREX_SPACEDIM+1)");
 debug_ngrow(idx,ngrow,local_caller_string);

} // subroutine allocate_levelset


void NavierStokes::resize_levelset(int ngrow,int idx) {

 if ((ngrow<0)||(ngrow>ngrow_distance))
  amrex::Error("ngrow invalid");

 if (localMF_grow[idx]>=0) {
  // do nothing
 } else {
  std::cout << "in resize_levelset  ngrow= " << ngrow << 
   " idx= " << idx << '\n';
  std::cout << "localMF_grow[idx] = " << localMF_grow[idx] << '\n';
  amrex::Error("localMF_grow[idx]<0");
 }

 if (localMF[idx]->nGrow()==ngrow) {
  // do nothing
 } else {
  allocate_levelset(ngrow,idx);
 }

} // subroutine resize_levelset

// called from make_physics_varsALL
// density vars get 1/rho
void NavierStokes::make_physics_vars(int project_option,
  const std::string& caller_string) {
 
 std::string local_caller_string="make_physics_vars";
 local_caller_string=caller_string+local_caller_string;

 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();

 int num_curv=num_interfaces*CURVCOMP_NCOMP; 

 if ((project_option==SOLVETYPE_PRES)||
     (project_option==SOLVETYPE_INITPROJ)) {
  // do nothing
 } else
  amrex::Error("project_option invalid make_physics_vars");

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
 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
  debug_ngrow(FSI_GHOST_MAC_MF+data_dir,0,local_caller_string);
 }

 if (localMF[DIST_CURV_MF]->nComp()!=num_curv)
  amrex::Error("localMF[DIST_CURV_MF]->nComp() invalid");

 VOF_Recon_resize(1); //output:SLOPE_RECON_MF
 debug_ngrow(SLOPE_RECON_MF,1,local_caller_string);

 if (localMF[SLOPE_RECON_MF]->nComp()!=num_materials*ngeom_recon)
  amrex::Error("slope_recon_mf has incorrect ncomp");

 debug_ngrow(DIST_CURV_MF,1,local_caller_string);

 resize_metrics(1);
 debug_ngrow(VOLUME_MF,1,local_caller_string);
 resize_maskfiner(1,MASKCOEF_MF);
 resize_mask_nbr(1);

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

  // in: make_physics_vars
 allocate_physics_vars();

 MultiFab& tempmf=get_new_data(State_Type,slab_step+1);

 const Real* dx = geom.CellSize();

 debug_ngrow(DEN_RECON_MF,1,local_caller_string);
 if (localMF[DEN_RECON_MF]->nComp()!=num_materials*num_state_material)
  amrex::Error("den_recon has invalid ncomp in make_physics_vars");
 debug_ngrow(MOM_DEN_MF,1,local_caller_string);
 if (localMF[MOM_DEN_MF]->nComp()!=num_materials)
  amrex::Error("mom_den has invalid ncomp in make_physics_vars");

 int REFINE_DENSITY_RECON_MF_local=-1;
 if ((num_materials_compressible>=1)&&
     (num_materials_compressible<=num_materials)) {
  REFINE_DENSITY_RECON_MF_local=REFINE_DENSITY_RECON_MF;
  debug_ngrow(REFINE_DENSITY_RECON_MF,1,local_caller_string);
  if (localMF[REFINE_DENSITY_RECON_MF]->nComp()!=NUM_CELL_REFINE_DENSITY)
   amrex::Error("REFINE_DENSITY_RECON_MF invalid ncomp in make_physics_vars");
 } else if (num_materials_compressible==0) {
  REFINE_DENSITY_RECON_MF_local=DEN_RECON_MF;
 } else
  amrex::Error("num_materials_compressble invalid");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  new_localMF(AMRSYNC_VEL_MF+dir,1,0,dir);
  setVal_localMF(AMRSYNC_VEL_MF+dir,1.0e+30,0,1,0);
 }

 Real problo[AMREX_SPACEDIM];
 Real probhi[AMREX_SPACEDIM];
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  problo[dir]=geom.ProbLo(dir);
  probhi[dir]=geom.ProbHi(dir);
 }
  
 localMF[CELL_SOUND_MF]->setVal(0.0,0,2,0);

 MultiFab* vofC=new MultiFab(grids,dmap,num_materials,1,
  MFInfo().SetTag("vofC"),FArrayBoxFactory());

 for (int im=0;im<num_materials;im++) {
  int scomp=im*ngeom_recon;
  MultiFab::Copy(*vofC,*localMF[SLOPE_RECON_MF],scomp,im,1,1);
 }

  // (dir-1)*2*num_materials + (side-1)*num_materials + im
 int nrefine_vof=2*num_materials*AMREX_SPACEDIM;
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
  FArrayBox& refinedenfab=(*localMF[REFINE_DENSITY_RECON_MF_local])[mfi];

  FArrayBox& vofFfab=(*vofF)[mfi];
  FArrayBox& massFfab=(*massF)[mfi];

  int tessellate=3;

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // tessellate=3
   //  If rigid materials dominate the cell, then that cell is considered
   //  to only have the one dominant rigid material (raster cell).  
   //  In the non-raster cells, the solids have no volume.
   //
   // declared in: LEVELSET_3D.F90
  fort_build_semirefinevof(
   &tid_current,
   &tessellate,  // =3
   &ngrow_refine,
   &nrefine_vof,
   spec_material_id_AMBIENT.dataPtr(),
   mass_fraction_id.dataPtr(),
   cavitation_vapor_density.dataPtr(),
   xlo,dx,
   slopefab.dataPtr(),
   ARLIM(slopefab.loVect()),ARLIM(slopefab.hiVect()),
   denstatefab.dataPtr(),
   ARLIM(denstatefab.loVect()),
   ARLIM(denstatefab.hiVect()),
   refinedenfab.dataPtr(),
   ARLIM(refinedenfab.loVect()),
   ARLIM(refinedenfab.hiVect()),
   mom_denfab.dataPtr(),
   ARLIM(mom_denfab.loVect()),ARLIM(mom_denfab.hiVect()),
   vofFfab.dataPtr(),ARLIM(vofFfab.loVect()),ARLIM(vofFfab.hiVect()),
   massFfab.dataPtr(),ARLIM(massFfab.loVect()),ARLIM(massFfab.hiVect()),
   tilelo,tilehi,
   fablo,fabhi,
   &bfact,
   &level,&finest_level);
 }  // mfi
} // omp
 ns_reconcile_d_num(LOOP_SEMIREFINEVOF,"make_physics_vars");

 resize_levelset(2,LEVELPC_MF);

 Vector< Real > local_curv_min;
 Vector< Real > local_curv_max;
 local_curv_min.resize(thread_class::nthreads);
 local_curv_max.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  local_curv_min[tid]=1.0e+30;
  local_curv_max[tid]=-1.0e+30;
 } // tid

 debug_ngrow(MASKCOEF_MF,1,local_caller_string);

 debug_ngrow(MASK_NBR_MF,1,local_caller_string);
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
   Vector<int> presbc=getBCArray(State_Type,gridno,STATECOMP_PRES,1);
   Vector<int> velbc=getBCArray(State_Type,gridno,
      STATECOMP_VEL,STATE_NCOMP_VEL);

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
   if (curvfab.nComp()==num_interfaces*CURVCOMP_NCOMP) {
    //do nothing
   } else
    amrex::Error("(curvfab.nComp()!=num_interfaces*CURVCOMP_NCOMP)");

   FArrayBox& denstatefab=(*localMF[DEN_RECON_MF])[mfi];
   FArrayBox& mom_denfab=(*localMF[MOM_DEN_MF])[mfi];

   FArrayBox& viscstatefab=(*localMF[CELL_VISC_MATERIAL_MF])[mfi];
   if (viscstatefab.nComp()!=3*num_materials)
    amrex::Error("viscstatefab.nComp()!=3*num_materials");

   FArrayBox& conductivity_fab=(*localMF[CELL_CONDUCTIVITY_MATERIAL_MF])[mfi];
   if (conductivity_fab.nComp()!=num_materials)
    amrex::Error("conductivity_fab.nComp()!=num_materials");

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

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();
 
    // in: LEVELSET_3D.F90
   fort_init_physics_vars(
    local_caller_string.c_str(),
    local_caller_string.size(),
    &tid_current,
    &FD_curv_interp, 
    &local_curv_min[tid_current],
    &local_curv_max[tid_current],
    &isweep,
    &nrefine_vof,
    denconst_interface.dataPtr(),
    denconst_interface_min.dataPtr(),
    viscconst_interface.dataPtr(),
    heatviscconst_interface.dataPtr(),
    speciesviscconst_interface.dataPtr(),
    freezing_model.dataPtr(),
    distribute_from_target.dataPtr(),
    &solidheat_flag,
    microlayer_size.dataPtr(), 
    microlayer_substrate.dataPtr(), 
    microlayer_temperature_substrate.dataPtr(), 
    spec_material_id_AMBIENT.dataPtr(),
    mass_fraction_id.dataPtr(),
    cavitation_vapor_density.dataPtr(),
    constant_density_all_time.dataPtr(),
    &cur_time_slab,
    &dt_slab, //calling fort_init_physics_vars
    &project_option,
    problo,probhi,
    &visc_coef,
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
    viscstatefab.dataPtr(), //3*num_materials components
    ARLIM(viscstatefab.loVect()),ARLIM(viscstatefab.hiVect()),
    conductivity_fab.dataPtr(), //num_materials components
    ARLIM(conductivity_fab.loVect()),ARLIM(conductivity_fab.hiVect()),
    solxfab.dataPtr(),
    ARLIM(solxfab.loVect()),ARLIM(solxfab.hiVect()),
    solyfab.dataPtr(),
    ARLIM(solyfab.loVect()),ARLIM(solyfab.hiVect()),
    solzfab.dataPtr(),
    ARLIM(solzfab.loVect()),ARLIM(solzfab.hiVect()),
    cDeDTfab.dataPtr(),
    ARLIM(cDeDTfab.loVect()),ARLIM(cDeDTfab.hiVect()),
    cdenfab.dataPtr(),
    ARLIM(cdenfab.loVect()),ARLIM(cdenfab.hiVect()),
    cvoffab.dataPtr(),ARLIM(cvoffab.loVect()),ARLIM(cvoffab.hiVect()),
    cviscfab.dataPtr(),ARLIM(cviscfab.loVect()),ARLIM(cviscfab.hiVect()),
    volfab.dataPtr(),ARLIM(volfab.loVect()),ARLIM(volfab.hiVect()),
    levelpcfab.dataPtr(),
    ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),
    vofCfab.dataPtr(),ARLIM(vofCfab.loVect()),ARLIM(vofCfab.hiVect()),
    vofFfab.dataPtr(),ARLIM(vofFfab.loVect()),ARLIM(vofFfab.hiVect()),
    massFfab.dataPtr(),ARLIM(massFfab.loVect()),ARLIM(massFfab.hiVect()),
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    presbc.dataPtr(), 
    velbc.dataPtr(), 
    &nparts,
    &nparts_def,
    im_solid_map_ptr,
    &num_curv, //num_interfaces * CURVCOMP_NCOMP
    &level,
    &finest_level);
  }  // mfi
} // omp
  ns_reconcile_d_num(LOOP_INIT_PHYSICS_VARS,"make_physics_vars");
 } // isweep

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
   std::cout << "c++ ngrow_distance " << ngrow_distance << '\n';

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
     tecplot_debug(curvfab,xlo,fablo,fabhi,dx,dir,0,FACECOMP_CURV,
      1,interior_only);
    } // mfi
    ns_reconcile_d_num(LOOP_TECPLOT_DEBUG_CURV_POSTPROC,"make_physics_vars");

   } // dir=0..sdim-1

 } // ((fab_verbose==1)||(fab_verbose==3))

 MultiFab::Copy(*localMF[CELL_VISC_HOLD_MF],*localMF[CELL_VISC_MF],0,0,1,1);
 MultiFab::Copy(*localMF[CELL_DEN_HOLD_MF],*localMF[CELL_DEN_MF],0,0,1,1);
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  MultiFab::Copy(
    *localMF[FACE_DEN_HOLD_MF+dir],
    *localMF[FACE_VAR_MF+dir],
    FACECOMP_FACEDEN,
    0,
    1,0);
  MultiFab::Copy(
   *localMF[FACE_VISC_HOLD_MF+dir],
   *localMF[FACE_VAR_MF+dir],
   FACECOMP_FACEVISC,
   0,
   1,0);
 }

 delete vofC;
 delete vofF;
 delete massF;
 
} // end subroutine make_physics_vars

//called from: NavierStokes::veldiffuseALL
void NavierStokes::solid_temperature() {
 
 bool use_tiling=ns_tiling;

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 MultiFab &S_new = get_new_data(State_Type,slab_step+1);
 int nstate=STATE_NCOMP;
 if (S_new.nComp()!=nstate) 
  amrex::Error("S_new.nComp()!=nstate");

 MultiFab &LS_new = get_new_data(LS_Type,slab_step+1);
 if (LS_new.nComp()!=num_materials*(AMREX_SPACEDIM+1)) 
  amrex::Error("LS_new.nComp()!=num_materials*(AMREX_SPACEDIM+1)");

 int solid_exists=0;
 for (int im=0;im<num_materials;im++) {
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

   int nden=num_materials*num_state_material;

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

     // fort_initsolidtemp declared in PROB.F90
    fort_initsolidtemp(
     &nden,
     &cur_time_slab,
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,
     snewfab.dataPtr(STATECOMP_STATES),
     ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()),
     lsnewfab.dataPtr(),
     ARLIM(lsnewfab.loVect()),ARLIM(lsnewfab.hiVect()),
     dx,xlo);  
 
   } // mfi
} // omp
   ns_reconcile_d_num(LOOP_INITSOLIDTEMP,"solid_temperature");
  } else
   amrex::Error("solidheat_flag invalid");

 } else
  amrex::Error("solid_exists invalid solid_temperature");

} // end subroutine solid_temperature

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

 std::string local_caller_string="increment_potential_force";

 int finest_level=parent->finestLevel();
 
 bool use_tiling=ns_tiling;

 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 int nstate=S_new.nComp();
 if (nstate!=STATE_NCOMP)
  amrex::Error("nstate invalid");
 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 MultiFab& LS_new=get_new_data(LS_Type,slab_step+1);
 if (LS_new.nComp()!=num_materials*(AMREX_SPACEDIM+1))
  amrex::Error("LS_new.nComp()!=num_materials*(AMREX_SPACEDIM+1)");

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

  debug_ngrow(POTENTIAL_FORCE_EDGE_MF+dir,0,local_caller_string);
  if (localMF[POTENTIAL_FORCE_EDGE_MF+dir]->nComp()!=1) {
   std::cout << "ncomp=" << 
    localMF[POTENTIAL_FORCE_EDGE_MF+dir]->nComp() << 
    " dir= " << dir << '\n';
   amrex::Error("localMF[POTENTIAL_FORCE_EDGE_MF+dir]->nComp() invalid");
  }

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
   FArrayBox& macfab=Umac_new[mfi];
   FArrayBox& facegrav=(*localMF[POTENTIAL_FORCE_EDGE_MF+dir])[mfi];
   FArrayBox& xfacefab=(*localMF[FACE_VAR_MF+dir])[mfi];
   FArrayBox& lsfab=LS_new[mfi];

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // declared in: NAVIERSTOKES_3D.F90
    // (gravity and surface tension)
    // u+=facegrav 
   fort_addgravity(
     &dt_slab,
     &cur_time_slab,
     &level,
     &finest_level,
     &nstate,
     tilelo,tilehi,
     fablo,fabhi,&bfact,
     xlo,dx,&dir,
     xfacefab.dataPtr(), 
     ARLIM(xfacefab.loVect()),ARLIM(xfacefab.hiVect()),
     lsfab.dataPtr(),ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
     macfab.dataPtr(),
     ARLIM(macfab.loVect()),ARLIM(macfab.hiVect()),
     facegrav.dataPtr(),
     ARLIM(facegrav.loVect()),ARLIM(facegrav.hiVect()) );
  } // mfi
} // omp
  ns_reconcile_d_num(LOOP_ADDGRAVITY,"increment_potential_force");
 }  // dir=0..sdim-1

} // increment_potential_force

// called from multiphase_project when 
// project_option==SOLVETYPE_PRES
void NavierStokes::deallocate_potential_forceALL() {

 int finest_level=parent->finestLevel();
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.delete_localMF(POTENTIAL_FORCE_EDGE_MF,AMREX_SPACEDIM);
  ns_level.delete_localMF(AMRSYNC_PRES_MF,AMREX_SPACEDIM);
 }
} // deallocate_potential_forceALL

// called from multiphase_project when 
// project_option==SOLVETYPE_PRES
void NavierStokes::process_potential_forceALL(
 int potgrad_surface_tension_mask,int project_option) {

 std::string local_caller_string="process_potential_forceALL";

 int finest_level=parent->finestLevel();
 if (level!=0)
  amrex::Error("level!=0");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 if (project_option==SOLVETYPE_PRES) {
  //do nothing
 } else
  amrex::Error("project_option invalid5794");

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
    ns_level.setVal_localMF(AMRSYNC_PRES_MF+dir,1.0e+30,0,2,0);
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
  ns_level.process_potential_force_face(potgrad_surface_tension_mask,
     project_option);
 }

  // must go from finest to coarsest in order
  // to average down face forces.
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  int ncomp_edge=-1;
  int scomp=0;
  // if level<finest_level then avgdown from level+1 to level.
  ns_level.avgDownEdge_localMF(POTENTIAL_FORCE_EDGE_MF,scomp,ncomp_edge,
   0,AMREX_SPACEDIM,1,local_caller_string);
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

 MultiFab* dendata=getStateDen(1,cur_time_slab);

 const Real* dx = geom.CellSize();
 const Box& domain = geom.Domain();
 const int* domlo = domain.loVect();
 const int* domhi = domain.hiVect();

 Vector<int> dombcpres(2*AMREX_SPACEDIM);
 const BCRec& descbc = get_desc_lst()[State_Type].getBC(STATECOMP_PRES);
 const int* b_rec=descbc.vect();
 for (int m=0;m<2*AMREX_SPACEDIM;m++)
  dombcpres[m]=b_rec[m];

 int bfact=parent->Space_blockingFactor(level);

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

    Vector<int> presbc=getBCArray(State_Type,gridno,STATECOMP_PRES,1);

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     // isweep=0 => interior cells updated
     // isweep=1 => exterior cells outside domain are updated:
     //   REFLECT_EVEN BC if wall, EXT_DIR BC on the ghost cell if
     //   outflow.
     // if angular_velocity==0.0:
     // \vec{g} = \frac{ \nabla rho_{0} (\vec{g} dot \vec{x})}{rho_{0}}
     // fort_init_potential is declared in: NAVIERSTOKES_3D.F90
    fort_init_potential(
     &cur_time_slab,
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
     &dt_slab, //fort_init_potential
     &angular_velocity, //fort_init_potential
     &centrifugal_force_factor, //fort_init_potential
     &isweep);
  } // mfi
} // omp
  ns_reconcile_d_num(LOOP_INIT_POTENTIAL,"init_gravity_potential");

  if (isweep==0) {

   if (level>0) {
    NavierStokes& coarse_lev = getLevel(level-1);
    int scomp=0;
    int ncomp=2; // HYDROSTATIC_PRESSURE and HYDROSTATIC_DENSITY
    coarse_lev.avgDown_localMF(HYDROSTATIC_PRESDEN_MF,scomp,ncomp,1);
   }

    // set_extrap_bc(bc,phys_bc)
    // fort_extrapfill
    // pc_interp or sem_interp
   Vector<int> scompBC_map;
   scompBC_map.resize(2);
   scompBC_map[0]=0;
   scompBC_map[1]=0;

   int extrap_enable_spectral=enable_spectral;
   override_enable_spectralGHOST(0,1,extrap_enable_spectral);
    //ngrow=1  scomp=0  ncomp=2
   PCINTERP_fill_borders(HYDROSTATIC_PRESDEN_MF,1,0,2,
     State_Type,scompBC_map);
   
   extrap_enable_spectral=0;
   override_enable_spectralGHOST(0,1,extrap_enable_spectral);

  } else if (isweep==1) {
   // do nothing
  } else
   amrex::Error("isweep invalid");
 
 } // for isweep = 0..1

 delete dendata;

}  // end subroutine init_gravity_potential

// called from: NavierStokes::process_potential_forceALL()
// u^face = u^face + facegravforce - grad^face p
// reflecting boundary conditions on ppot should be identical to the
// reflecting boundary conditions on p so that
// if facegrav = grad^face ppot, then this implies that
//  div (grad^face ppot - grad^face p)=0 only if ppot=p.
void NavierStokes::process_potential_force_face(
  int potgrad_surface_tension_mask,
  int project_option) {

 std::string local_caller_string="process_potential_force_face";

 int finest_level=parent->finestLevel();

 int operation_flag=OP_POTGRAD_TO_MAC;
 
 bool use_tiling=ns_tiling;

 if ((SDC_outer_sweeps>=0)&&(SDC_outer_sweeps<ns_time_order)) {
  // do nothing
 } else
  amrex::Error("SDC_outer_sweeps invalid process_potential_force_face");

 if (project_option==SOLVETYPE_PRES) {
  //do nothing
 } else
  amrex::Error("project_option invalid 6010 process_pot_force_face ");

 int num_colors=0;
 Vector<Real> blob_array;
 blob_array.resize(1);
 int blob_array_size=blob_array.size();

 int nsolve=1;

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

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 resize_maskfiner(1,MASKCOEF_MF);
 resize_mask_nbr(1);
 resize_metrics(1);
 resize_levelset(2,LEVELPC_MF);

 // maskcoef=1 if not covered by finer level.
 debug_ngrow(MASKCOEF_MF,1,local_caller_string); 
 // mask_nbr=1 at fine-fine bc.
 debug_ngrow(MASK_NBR_MF,1,local_caller_string); 

 MultiFab &S_new = get_new_data(State_Type,slab_step+1);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  new_localMF(POTENTIAL_FORCE_EDGE_MF+dir,1,0,dir);//grad ppot/rhopot
 }

 const Real* dx = geom.CellSize();
 const Box& domain = geom.Domain();
 const int* domlo = domain.loVect();
 const int* domhi = domain.hiVect();

 Vector<int> dombcpres(2*AMREX_SPACEDIM);
 const BCRec& descbc = get_desc_lst()[State_Type].getBC(STATECOMP_PRES);
 const int* b_rec=descbc.vect();
 for (int m=0;m<2*AMREX_SPACEDIM;m++)
  dombcpres[m]=b_rec[m];

 MultiFab* dendata=getStateDen(1,cur_time_slab);

  // gpx/rhox,gpy/rhoy,gpz/rhoz
 allocate_flux_register(operation_flag);
 if (localMF[SEM_FLUXREG_MF]->nComp()!=AMREX_SPACEDIM)
  amrex::Error("localMF[SEM_FLUXREG_MF]->nComp() invalid7");

  // enable_spectral=0 => end_spectral_loop()=1
  // enable_spectral=1 => end_spectral_loop()=(bfact==1 ? 1:2)
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
   if (xgp.nComp()==1) {
    //do nothing
   } else
    amrex::Error("xgp.nComp() invalid");

   FArrayBox& xp=(*localMF[POTENTIAL_FORCE_EDGE_MF+dir])[mfi];
   if (xp.nComp()==1) {
    //do nothing
   } else
    amrex::Error("xp.nComp() invalid");

   FArrayBox& xface=(*localMF[FACE_VAR_MF+dir])[mfi];
 
   FArrayBox& maskfab=(*localMF[MASK_NBR_MF])[mfi];

   // mask=tag if not covered by level+1 or outside the domain.
   FArrayBox& maskcoef=(*localMF[MASKCOEF_MF])[mfi];

   FArrayBox& maskSEMfab=(*localMF[MASKSEM_MF])[mfi];
   FArrayBox& presdenfab=(*localMF[HYDROSTATIC_PRESDEN_MF])[mfi];
   FArrayBox& mgonifab=(*dendata)[mfi];
   if (mgonifab.nComp()==num_materials*num_state_material) {
    // do nothing
   } else
    amrex::Error("mgonifab.nComp()!=num_materials*num_state_material");

   FArrayBox& solfab=(*localMF[FSI_GHOST_MAC_MF+dir])[mfi];
   FArrayBox& levelpcfab=(*localMF[LEVELPC_MF])[mfi];

   Vector<int> presbc=getBCArray(State_Type,gridno,STATECOMP_PRES,1);
   Vector<int> velbc=getBCArray(State_Type,gridno, 
                                STATECOMP_VEL,STATE_NCOMP_VEL);

   Real beta=0.0;

   FArrayBox& semfluxfab=(*localMF[SEM_FLUXREG_MF])[mfi];
   int ncfluxreg=semfluxfab.nComp();
   if (ncfluxreg!=AMREX_SPACEDIM) 
    amrex::Error("ncfluxreg invalid");

   int local_energyflag=SUB_OP_FORCE_MASK_BASE+potgrad_surface_tension_mask;
   int local_enable_spectral=enable_spectral;

   int simple_AMR_BC_flag=1;

   int ncomp_xp=1;
   int ncomp_xgp=1;
   int ncomp_mgoni=mgonifab.nComp();

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();
   int ncphys_proxy=FACECOMP_NCOMP;

   // process_potential_force_face 
   fort_cell_to_mac( 
    &ncomp_mgoni,
    &ncomp_xp,
    &ncomp_xgp,
    &simple_AMR_BC_flag,
    &nsolve,
    &tileloop,
    &dir,
    &operation_flag, //OP_POTGRAD_TO_MAC
    &local_energyflag,//SUB_OP_FORCE_MASK_BASE+potgrad_surface_tension_mask
    &beta,
    &visc_coef,
    &local_enable_spectral,
    &ncphys_proxy,
    constant_density_all_time.dataPtr(),
    presbc.dataPtr(),
    velbc.dataPtr(),
    &slab_step,
    &dt_slab, // fort_cell_to_mac,process_potential_force_face
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
    &NS_geometry_coord,
    domlo,domhi,
    &nparts,
    &nparts_def,
    im_solid_map_ptr,
    blob_array.dataPtr(),
    &blob_array_size,
    &num_colors,
    &project_option);
  } // mfi
} // omp
  ns_reconcile_d_num(LOOP_POTGRAD_TO_MAC,"process_potential_force_face");
 } // tileloop
 } // dir
 synchronize_flux_register(operation_flag,spectral_loop);
 } // spectral_loop

 delete dendata;

}  // subroutine process_potential_force_face


// typically ngrow=4 (see call in NavierStokes3.cpp)
void NavierStokes::metrics_dataALL(int ngrow) {

 if (level==0) {
  // do nothing
 } else
  amrex::Error("level invalid in metrics_dataALL");

 int finest_level=parent->finestLevel();

 for (int ilev=0;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.metrics_data(ngrow);
 }
}


void NavierStokes::metrics_data_min_max_ALL(
  const std::string& caller_string) {

 std::fflush(NULL);

 int finest_level=parent->finestLevel();
 for (int ilev=0;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.metrics_data_min_max(caller_string);
 }
 std::fflush(NULL);
}


void NavierStokes::metrics_data(int ngrow) {
 
 int finest_level=parent->finestLevel();
 if ((level>=0)&&(level<=finest_level)) {
  // do nothing
 } else
  amrex::Error("level invalid in metrics_data");

 if (1==0) {
  std::fflush(NULL);
  int proc=ParallelDescriptor::MyProc();
  std::cout << "in metrics_data on processor=" << proc << '\n';
  std::cout << "level= " << level << '\n';
  std::cout << "finest_level= " << finest_level << '\n';
  std::cout << "ngrow= " << ngrow << '\n';
  std::fflush(NULL);
  ParallelDescriptor::Barrier();
 }
  
 bool use_tiling=ns_tiling;

 const Real* dx = geom.CellSize();

 delete_localMF_if_exist(VOLUME_MF,1);
 delete_localMF_if_exist(AREA_MF,AMREX_SPACEDIM);

 if (1==0) {
  std::fflush(NULL);
  int proc=ParallelDescriptor::MyProc();
  std::cout << "in metrics_data2 on processor=" << proc << '\n';
  std::cout << "level= " << level << '\n';
  std::cout << "finest_level= " << finest_level << '\n';
  std::cout << "ngrow= " << ngrow << '\n';
  std::fflush(NULL);
  ParallelDescriptor::Barrier();
 }

 if (ngrow<0)
  amrex::Error("ngrow too small");

 if (1==0) {
  std::fflush(NULL);
  int proc=ParallelDescriptor::MyProc();
  std::cout << "in metrics_data22 on processor=" << proc << '\n';
  std::cout << "level= " << level << '\n';
  std::cout << "finest_level= " << finest_level << '\n';
  std::cout << "ngrow= " << ngrow << '\n';
  std::fflush(NULL);
  ParallelDescriptor::Barrier();
 }

 new_localMF(VOLUME_MF,1,ngrow,-1);

 if (1==0) {
  std::fflush(NULL);
  int proc=ParallelDescriptor::MyProc();
  std::cout << "in metrics_data3 on processor=" << proc << '\n';
  std::cout << "level= " << level << '\n';
  std::cout << "finest_level= " << finest_level << '\n';
  std::cout << "ngrow= " << ngrow << '\n';
  std::fflush(NULL);
  ParallelDescriptor::Barrier();
 }
 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  new_localMF(AREA_MF+dir,1,ngrow,dir);
 }

 if (1==0) {
  std::fflush(NULL);
  int proc=ParallelDescriptor::MyProc();
  std::cout << "in metrics_data4 on processor=" << proc << '\n';
  std::cout << "level= " << level << '\n';
  std::cout << "finest_level= " << finest_level << '\n';
  std::cout << "ngrow= " << ngrow << '\n';
  std::fflush(NULL);
  ParallelDescriptor::Barrier();
 }
 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 const BoxArray mfBA=localMF[VOLUME_MF]->boxArray();

 if (1==0) {
  std::fflush(NULL);
  int proc=ParallelDescriptor::MyProc();
  std::cout << "in metrics_data5 on processor=" << proc << '\n';
  std::cout << "level= " << level << '\n';
  std::cout << "finest_level= " << finest_level << '\n';
  std::cout << "ngrow= " << ngrow << '\n';
  std::fflush(NULL);
  ParallelDescriptor::Barrier();
 }
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

 if (1==0) {
  std::fflush(NULL);
  int proc=ParallelDescriptor::MyProc();
  std::cout << "in metrics_data6 on processor=" << proc << '\n';
  std::cout << "level= " << level << '\n';
  std::cout << "finest_level= " << finest_level << '\n';
  std::cout << "ngrow= " << ngrow << '\n';
  std::fflush(NULL);
  ParallelDescriptor::Barrier();
 }
#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*localMF[VOLUME_MF],use_tiling); mfi.isValid();++mfi) {

  if (1==0) {
   std::fflush(NULL);
   int proc=ParallelDescriptor::MyProc();
   std::cout << "in metrics_data7 on processor=" << proc << '\n';
   std::cout << "level= " << level << '\n';
   std::cout << "finest_level= " << finest_level << '\n';
   std::cout << "ngrow= " << ngrow << '\n';
   std::fflush(NULL);
  }
  if (grids[mfi.index()] == mfi.validbox()) {
   // do nothing
  } else
   amrex::Error("grids[mfi.index()] == mfi.validbox() failed");

  const int gridno = mfi.index();

  if (1==0) {
   std::fflush(NULL);
   int proc=ParallelDescriptor::MyProc();
   std::cout << "prior to fort_metrics on processor=" << proc << '\n';
   std::cout << "gridno= " << gridno << '\n';
   std::fflush(NULL);
  }

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

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   //declared in: NAVIERSTOKES_3D.F90
  fort_metrics(
   xlo,dx,
   areax.dataPtr(),ARLIM(areax.loVect()),ARLIM(areax.hiVect()),
   areay.dataPtr(),ARLIM(areay.loVect()),ARLIM(areay.hiVect()),
   areaz.dataPtr(),ARLIM(areaz.loVect()),ARLIM(areaz.hiVect()),
   volfab.dataPtr(),ARLIM(volfab.loVect()),ARLIM(volfab.hiVect()),
   tilelo,tilehi,
   fablo,fabhi,
   &bfact,
   &level,
   &finest_level,
   &ngrow,
   &NS_geometry_coord);
 }  // mfi
}  // omp
 ns_reconcile_d_num(LOOP_METRICS,"metrics_data");

} // subroutine metrics_data


void NavierStokes::metrics_data_min_max(const std::string& caller_string) {

 int ngrow=localMF_grow[VOLUME_MF];
 std::cout << "metrics_data_min_max caller_string " << caller_string <<'\n';
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
  int renormalize_only,int local_truncate,
  const std::string& caller_string,
  int update_particles) {

 if (level!=0)
  amrex::Error("level should be 0 in prescribe_solid_geometryALL");

 int finest_level=parent->finestLevel();

 std::string local_caller_string="prescribe_solid_geometryALL";
 local_caller_string=caller_string+local_caller_string;

 if (renormalize_only==1) {
  if (update_particles==0) {
   //do nothing
  } else
   amrex::Error("expecting update_particles=0");
 } else if (renormalize_only==0) {
  if ((update_particles==0)||(update_particles==1)) {
   //do nothing
  } else
   amrex::Error("expecting update_particles=0 or 1");
 } else
  amrex::Error("expecting renormalize_only=0 or 1");

 if (local_truncate==1) {
  Vector<Real> delta_mass_all;
  delta_mass_all.resize(num_materials);

  for (int ilev=finest_level;ilev>=level;ilev--) {
   NavierStokes& ns_level=getLevel(ilev);
   if (ilev<finest_level) {
    ns_level.avgDown(LS_Type,0,num_materials,0);
    ns_level.MOFavgDown();
   }
   ns_level.truncate_VOF(delta_mass_all);
   if (verbose>0) {
    if (ParallelDescriptor::IOProcessor()) {
     for (int im=0;im<num_materials;im++) {
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

  if (std::abs(time-cur_time_slab)>CPP_EPS_8_5)
   amrex::Error("prescribe solid at the new time");

   //init_FSI_GHOST_MAC_MF_ALL is declared in NavierStokes.cpp
  init_FSI_GHOST_MAC_MF_ALL(renormalize_only,local_caller_string);
 
  interface_touch_flag=1; //prescribe_solid_geometryALL
			  
 } else if (renormalize_only==1) {
  // do nothing
 } else
  amrex::Error("renormalize_only invalid");

 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  if (ilev<finest_level) {
   ns_level.avgDown(LS_Type,0,num_materials,0);
   ns_level.MOFavgDown();
  }
  ns_level.prescribe_solid_geometry(time,renormalize_only);
 } //ilev=finest_level downto level
 
 if (update_particles==1) {

#ifdef AMREX_PARTICLES

   //calling from NavierStokes::prescribe_solid_geometryALL
  if ((slab_step>=0)&&(slab_step<ns_time_order)) {
   init_particle_containerALL(OP_PARTICLE_ADD,local_caller_string);
  } else
   amrex::Error("slab_step invalid");

  My_ParticleContainer& localPC_DIST=newDataPC(slab_step+1);
  int lev_min_DIST=0;
  int lev_max_DIST=-1;
  int nGrow_Redistribute_DIST=0;
  int local_Redistribute_DIST=0; 
  localPC_DIST.Redistribute(lev_min_DIST,lev_max_DIST,
    nGrow_Redistribute_DIST,local_Redistribute_DIST);

#endif

 } else if (update_particles==0) {
  //do nothing
 } else
  amrex::Error("update_particles invalid prescribe_solid_geometryALL");

} // end subroutine prescribe_solid_geometryALL

#ifdef AMREX_PARTICLES

void NavierStokes::move_particles(
  int splitting_dir,
  My_ParticleContainer& localPC,
  const std::string& caller_string) {

 if ((slab_step>=0)&&(slab_step<ns_time_order)) {
  //do nothing
 } else
  amrex::Error("slab_step invalid");

 std::string local_caller_string="move_particles";
 local_caller_string=caller_string+local_caller_string;

#if (NS_profile_solver==1)
 BLProfiler bprof(local_caller_string);
#endif

 if (pattern_test(local_caller_string,"do_the_advance")==1) {
  //do nothing
 } else {
  std::cout << local_caller_string << '\n';
  amrex::Error("caller is invalid in move_particles");
 }

 int phase_change_displacement=0;

 if (pattern_test(local_caller_string,"nonlinear_advection")==1) {
  //do nothing
 } else if (pattern_test(local_caller_string,"phase_change_code_segment")==1) {
  phase_change_displacement=1;
 } else {
  std::cout << local_caller_string << '\n';
  amrex::Error("caller is invalid in move_particles");
 }

 bool use_tiling=ns_tiling;

 if (1==1) {
  use_tiling=false;
 }

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

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");
 
 MultiFab& S_new = get_new_data(State_Type,slab_step+1);

 const Real* dx = geom.CellSize();
 const Box& domain = geom.Domain();
 const int* domlo = domain.loVect();
 const int* domhi = domain.hiVect();

 Vector<int> dombc(2*AMREX_SPACEDIM);
 const BCRec& descbc = get_desc_lst()[State_Type].getBC(STATECOMP_MOF);
 const int* b_rec=descbc.vect();
 for (int m=0;m<2*AMREX_SPACEDIM;m++)
  dombc[m]=b_rec[m];

 MultiFab* lsmf=getStateDist(1,cur_time_slab,local_caller_string); 

 MultiFab* mac_velocity[AMREX_SPACEDIM];
 MultiFab* burning_velocity;

 int burning_velocity_ncomp=1;

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  if (phase_change_displacement==0) {

   if ((splitting_dir>=0)&&(splitting_dir<AMREX_SPACEDIM)) {
    //do nothing
   } else
    amrex::Error("splitting_dir invalid");

   mac_velocity[dir]=localMF[RAW_MAC_VELOCITY_MF+dir];
   
   burning_velocity=lsmf;
   burning_velocity_ncomp=1;

  } else if (phase_change_displacement==1) {

   if (splitting_dir==-1) {
    //do nothing
   } else
    amrex::Error("splitting_dir invalid");

   mac_velocity[dir]=lsmf;
   burning_velocity=localMF[BURNING_VELOCITY_MF];
   burning_velocity_ncomp=EXTRAP_NCOMP_BURNING;

   if (localMF[BURNING_VELOCITY_MF]->nComp()!=EXTRAP_NCOMP_BURNING)
    amrex::Error("localMF[BURNING_VELOCITY_MF] incorrect ncomp");
   if (localMF[BURNING_VELOCITY_MF]->nGrow()!=ngrow_distance)
    amrex::Error("localMF[BURNING_VELOCITY_MF] incorrect ngrow");
   debug_ixType(BURNING_VELOCITY_MF,-1,local_caller_string);

  } else
   amrex::Error("phase_change_displacement invalid");
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
  int bfact=parent->Space_blockingFactor(level);

  FArrayBox& xvelfab=(*mac_velocity[0])[mfi];
  FArrayBox& yvelfab=(*mac_velocity[1])[mfi];
  FArrayBox& zvelfab=(*mac_velocity[AMREX_SPACEDIM-1])[mfi];

  FArrayBox& lsfab=(*lsmf)[mfi];

  FArrayBox& burningfab=(*burning_velocity)[mfi];

  const Real* xlo = grid_loc[gridno].lo();

    // this is an object with a pointer to AoS data
  auto& particles = localPC.GetParticles(level)
   [std::make_pair(mfi.index(),mfi.LocalTileIndex())];

  auto& particles_AoS = particles.GetArrayOfStructs();
  int Np=particles_AoS.size();

  Vector<int> denbc=getBCArray(State_Type,gridno,STATECOMP_STATES,
    num_materials*num_state_material);
  Vector<int> velbc=getBCArray(State_Type,gridno,
    STATECOMP_VEL,STATE_NCOMP_VEL);

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // declared in: LEVELSET_3D.F90
  fort_move_particle_container( 
   &splitting_dir,
   &phase_change_displacement,
   &burning_velocity_ncomp,
   &tid_current,
   tilelo,tilehi,
   fablo,fabhi,
   &bfact,
   &level,
   &finest_level,
   xlo,dx,
   particles_AoS.data(),
   Np,  // pass by value
   &dt_slab, //move_particle_container
   &vel_time_slab,
   xvelfab.dataPtr(),
   ARLIM(xvelfab.loVect()),ARLIM(xvelfab.hiVect()),
   yvelfab.dataPtr(),
   ARLIM(yvelfab.loVect()),ARLIM(yvelfab.hiVect()),
   zvelfab.dataPtr(),
   ARLIM(zvelfab.loVect()),ARLIM(zvelfab.hiVect()),
   burningfab.dataPtr(),
   ARLIM(burningfab.loVect()),
   ARLIM(burningfab.hiVect()),
   lsfab.dataPtr(),
   ARLIM(lsfab.loVect()),
   ARLIM(lsfab.hiVect()),
   velbc.dataPtr(),
   denbc.dataPtr(),
   dombc.dataPtr(),
   domlo,domhi);

 }  // mfi
} // omp
 ns_reconcile_d_num(LOOP_MOVE_PARTICLE_CONTAINER,"move_particles");

 using MyParIter=My_ParticleContainer::ParIterType;
 for (MyParIter pti(localPC,level);pti.isValid();++pti) {
  auto& particles=pti.GetArrayOfStructs();
  for (auto& p : particles) {
   int primary_material_id=p.idata(N_EXTRA_INT_PRIMARY_MATERIAL_ID);
   int secondary_material_id=p.idata(N_EXTRA_INT_SECONDARY_MATERIAL_ID);
   if ((primary_material_id>=1)&&
       (primary_material_id<=num_materials)&&
       (secondary_material_id>=1)&&
       (secondary_material_id<=num_materials)) {
    //do nothing
   } else if ((primary_material_id==-1)&&
              (secondary_material_id==-1)) {
    p.id()=-p.id();
   } else
    amrex::Error("primary_material_id or secondary_material_id invalid");
  }
 }

 delete lsmf;

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  if (phase_change_displacement==0) {
   //do nothing
  } else if (phase_change_displacement==1) {
   //do nothing
  } else
   amrex::Error("phase_change_displacement invalid");
 }

#if (NS_profile_solver==1)
 bprof.stop();
#endif

} // end subroutine move_particles

#endif

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
 
 std::string local_caller_string="prescribe_solid_geometry";

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

 if (NUM_CELL_REFINE_DENSITY==
     num_materials_compressible*ENUM_NUM_REFINE_DENSITY_TYPE) {
  // do nothing
 } else
  amrex::Error("NUM_CELL_REFINE_DENSITY invalid");

 int Refine_Density_Type_local=-1;
 if ((num_materials_compressible>=1)&&
     (num_materials_compressible<=num_materials)) {
  Refine_Density_Type_local=Refine_Density_Type;
 } else if (num_materials_compressible==0) {
  Refine_Density_Type_local=State_Type;
 } else
  amrex::Error("num_materials_compressble invalid");

 MultiFab& Refine_Density_new=
    get_new_data(Refine_Density_Type_local,slab_step+1);

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,1,local_caller_string);

 const Real* dx = geom.CellSize();

  // renormalize_only==1:
  //   project so that sum F_m_fluid=1
  // renormalize_only==0:
  //   correct F_m according to prescribed solid interface.
  //   project so that sum F_m_fluid=1 

 // nparts x (velocity + LS + temperature + flag)
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

   MultiFab* veldata=getState(1,STATECOMP_VEL, 
		STATE_NCOMP_VEL+STATE_NCOMP_PRES,time); 
   MultiFab* mofdata=getState(1,STATECOMP_MOF,num_materials*ngeom_raw,time);
   MultiFab* dendata=getStateDen(1,time);
   MultiFab* lsdata=getStateDist(ngrow_distance,time,local_caller_string);

   for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
    if (localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0)
     amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0");
    if (localMF[FSI_GHOST_MAC_MF+data_dir]->nComp()!=nparts_def*AMREX_SPACEDIM)
     amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nComp() bad");
   }

   if (veldata->nComp()!=(STATE_NCOMP_VEL+STATE_NCOMP_PRES))
    amrex::Error("veldata incorrect ncomp");
   if (dendata->nComp()!=num_materials*num_state_material)
    amrex::Error("dendata incorrect ncomp");
   if (mofdata->nComp()!=num_materials*ngeom_raw)
    amrex::Error("mofdata incorrect ncomp");
   if (lsdata->nComp()!=num_materials*(1+AMREX_SPACEDIM))
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

    FArrayBox& refinedennew=Refine_Density_new[mfi];

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
    fort_renormalize_prescribe(
      &tid_current,
      &level,&finest_level,
      &time,
      tilelo,tilehi,
      fablo,fabhi,&bfact,
      vofnew.dataPtr(STATECOMP_MOF),
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
      dennew.dataPtr(STATECOMP_STATES),
      ARLIM(dennew.loVect()),ARLIM(dennew.hiVect()),
      refinedennew.dataPtr(),
      ARLIM(refinedennew.loVect()),ARLIM(refinedennew.hiVect()),
      lsnew.dataPtr(),
      ARLIM(lsnew.loVect()),ARLIM(lsnew.hiVect()),
      xlo,dx,
      &cur_time_slab,
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
   ns_reconcile_d_num(LOOP_RENORMALIZE_PRESCRIBE,"prescribe_solid_geometry");

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


// called from NavierStokes::prescribe_solid_geometryALL
void NavierStokes::truncate_VOF(Vector<Real>& delta_mass_all) {

 std::string local_caller_string="truncate_VOF";
 
 bool use_tiling=ns_tiling;

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int finest_level = parent->finestLevel();

 MultiFab &S_new = get_new_data(State_Type,slab_step+1);
 int nc=STATE_NCOMP;
 if (nc!=S_new.nComp())
  amrex::Error("nc invalid in truncate_VOF");

 if (delta_mass_all.size()!=num_materials)
  amrex::Error("delta_mass_all has invalid size");

 Vector< Vector<Real> > local_delta_mass;
 local_delta_mass.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  local_delta_mass[tid].resize(num_materials); 
  for (int im=0;im<num_materials;im++)
   local_delta_mass[tid][im]=0.0;
 } // tid

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,1,local_caller_string);

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
   fort_purgeflotsam(
     local_delta_mass[tid_current].dataPtr(),
     &level,&finest_level,
     &cur_time_slab,
     tilelo,tilehi,
     fablo,fabhi,&bfact,
     maskcov.dataPtr(),
     ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
     vofnew.dataPtr(STATECOMP_MOF),
     ARLIM(vofnew.loVect()),ARLIM(vofnew.hiVect()),
     xlo,dx);

 }  // mfi
} // omp
 ns_reconcile_d_num(LOOP_PURGEFLOTSAM,"truncate_VOF");

 for (int tid=1;tid<thread_class::nthreads;tid++) {
  for (int im=0;im<num_materials;im++) {
   local_delta_mass[0][im]+=local_delta_mass[tid][im];
  }
 } // tid

 ParallelDescriptor::Barrier();
 for (int im=0;im<num_materials;im++) {
  ParallelDescriptor::ReduceRealSum(local_delta_mass[0][im]);
 }
 for (int im=0;im<num_materials;im++) {
  delta_mass_all[im]+=local_delta_mass[0][im];
 }

}  // truncate_VOF()

void NavierStokes::output_triangles() {

 std::string local_caller_string="output_triangles";

 int finest_level = parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid output_triangles");

// bool use_tiling=ns_tiling;
 bool use_tiling=false;

  // mask=tag if not covered by level+1 or outside the domain.
 Real tag=1.0;
 int clearbdry=0;
 MultiFab* maskplot=maskfiner(1,tag,clearbdry);

   // vof,refcentroid,order,slope,intercept x num_materials

 VOF_Recon_resize(1); //output:SLOPE_RECON_MF
 debug_ngrow(SLOPE_RECON_MF,1,local_caller_string);

 const Real* dx = geom.CellSize();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(localMF[SLOPE_RECON_MF]->boxArray().d_numPts());

// no open mp and no tiling.
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
   &level,&gridno);

#ifdef AMREX_PARTICLES

  NavierStokes& ns_level0=getLevel(0);
  My_ParticleContainer& localPC=ns_level0.newDataPC(slab_step+1);

  auto& particles_grid_tile = localPC.GetParticles(level)
   [std::make_pair(mfi.index(),mfi.LocalTileIndex())];
  auto& particles_AoS = particles_grid_tile.GetArrayOfStructs();
  unsigned int Np=particles_AoS.size();

  // declared in: NAVIERSTOKES_3D.F90
  // output to:
  // ./temptecplot/tempPARCON_pos<level><gridno>
  fort_particle_grid(
   &tid_current,
   xlo,dx,
   particles_AoS.data(),
   Np,       //pass by value
   tilelo,tilehi,
   fablo,fabhi,
   &bfact,
   &level,
   &gridno);

#endif

 }  // mfi
 ns_reconcile_d_num(LOOP_ISOGRID,"output_triangles");

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
} // end subroutine NavierStokes::tecplot_debug

// datatype=0 standard case
// datatype=1 face grad U
// datatype=2 cell grad U
Box NavierStokes::growntileboxTENSOR(
 MultiFab* mf,int datatype,int ng,int dir,
 const Box& tilebox_box,const Box& vbx_box) {

 Box bx = tilebox_box;

 if (datatype==0) {

  if (ng < -100) ng = mf->nGrow();
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
   if (bx.smallEnd(d) == vbx_box.smallEnd(d)) {
    bx.growLo(d, ng);
   }
   if (bx.bigEnd(d) == vbx_box.bigEnd(d)) {
    bx.growHi(d, ng);
   }
  } // d

 } else if ((datatype==1)||(datatype==2)) {

  if ((dir<0)||(dir>=AMREX_SPACEDIM))
   amrex::Error("dir invalid");
  if (ng!=0)
   amrex::Error("ng invalid");
  
  for (int d=0; d<AMREX_SPACEDIM; ++d) {

   if (!mf->ixType().cellCentered(d))
    amrex::Error("tensor box should be cell centered");

   if (d!=dir) {
    if (bx.smallEnd(d) == vbx_box.smallEnd(d)) {
     bx.growLo(d,1);
    }
    if (bx.bigEnd(d) == vbx_box.bigEnd(d)) {
     bx.growHi(d, 1);
    }
   } else if (d==dir) {
    if (datatype==1) {
     if (bx.bigEnd(d) == vbx_box.bigEnd(d)) {
      bx.growHi(d, 1);
     }
    } else if (datatype==2) {
     // do nothing
    } else
     amrex::Error("datatype invalid");
   }
  } // d

 } else
  amrex::Error("datatype invalid");

 return bx;

} // end subroutine Box NavierStokes::growntileboxTENSOR


bool NavierStokes::contains_nanTENSOR(MultiFab* mf,
  int datatype,int scomp,int dir) {

 if ((scomp<0)||(scomp>=mf->nComp()))
  amrex::Error("scomp invalid");
 if ((dir<0)||(dir>=AMREX_SPACEDIM))
  amrex::Error("dir invalid");

 bool r = false;

#ifdef _OPENMP
#pragma omp parallel reduction(|:r)
#endif
 for (MFIter mfi(*mf,true); mfi.isValid(); ++mfi) {
  BL_ASSERT(grids[mfi.index()] == mfi.validbox());
  const int gridno = mfi.index();
  const Box& tilegrid = mfi.tilebox();
  const Box& fabgrid = grids[gridno];
  int ng=0;

  const Box& bx = growntileboxTENSOR(mf,datatype,ng,dir,
    tilegrid,fabgrid);

  FArrayBox& fab_mf=(*mf)[mfi];
  Array4<Real> const& fab_mf_array=fab_mf.array();
  const Dim3 lo3=amrex::lbound(bx);
  const Dim3 hi3=amrex::ubound(bx);
  for (int z=lo3.z;z<=hi3.z;++z) {
  for (int y=lo3.y;y<=hi3.y;++y) {
  for (int x=lo3.x;x<=hi3.x;++x) {
   Real test_value=fab_mf_array(x,y,z,scomp);
   if (test_value==test_value) {
    //do nothing
   } else {
    r=true;
   }
  }
  }
  }
 } // mfi

 ParallelDescriptor::ReduceBoolOr(r);

 return r;
} // subroutine contains_nanTENSOR


// datatype=0 normal
// datatype=1 tensor face
// datatype=2 tensor cell
void NavierStokes::check_for_NAN_TENSOR_base(int datatype,MultiFab* mf,
  int sc,int dir) {

 int finest_level=parent->finestLevel();
 int ncomp=mf->nComp();
 int ngrow=mf->nGrow();
 const BoxArray mfBA=mf->boxArray();
 int ngrid=mfBA.size();
 const Box& domain = geom.Domain();

 if (ngrow!=1)
  amrex::Error("ngrow invalid");

 std::fflush(NULL);

 if (contains_nanTENSOR(mf,datatype,sc,dir)==true) {
  std::cout << "sc= " << sc << '\n';
  std::cout << "dir= " << dir << '\n';
  std::cout << "ncomp= " << ncomp << '\n';
  std::cout << "ngrow= " << ngrow << '\n';
  std::cout << "level= " << level << '\n';
  std::cout << "finest_level= " << finest_level << '\n';
  std::cout << "domain= " << domain << '\n';
  std::cout << "ngrid= " << ngrid << '\n';
  std::cout << "mfBA= " << mfBA << '\n';
  amrex::Error("mf contains nan ::check_for_NAN_TENSOR_base");
 }
 int force_check=1;
 int ncomp_tensor=1;
 int ngrow_tensor=0;
 Real warning_cutoff=1.0e+30;
 aggressive_debug(
  datatype,
  force_check,
  mf,
  sc,ncomp_tensor,
  ngrow_tensor,
  dir,
  warning_cutoff);

 std::fflush(NULL);

} // subroutine check_for_NAN_TENSOR_base



// datatype=0 normal
// datatype=1 face grad U
// datatype=2 cell grad U
void NavierStokes::check_for_NAN_TENSOR(int datatype,MultiFab* mf) {

 int finest_level=parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid check_for_NAN_TENSOR");

 int ncomp=mf->nComp();
 int ngrow=mf->nGrow();

 if (ncomp!=AMREX_SPACEDIM_SQR)
  amrex::Error("ncomp invalid in check_for_NAN_TENSOR");
 if (ngrow!=1)
  amrex::Error("ngrow invalid");

 std::fflush(NULL);

 for (int sc=0;sc<AMREX_SPACEDIM_SQR;sc++) {

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

   check_for_NAN_TENSOR_base(datatype,mf,sc,dir);
 }  // sc=0...AMREX_SPACEDIM_SQR-1

 
 std::fflush(NULL);

} // subroutine check_for_NAN_TENSOR


void NavierStokes::check_for_NAN(MultiFab* mf) {

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
    std::cout << "sc= " << sc << '\n';
    std::cout << "ng= " << ng << '\n';
    std::cout << "ncomp= " << ncomp << '\n';
    std::cout << "ngrow= " << ngrow << '\n';
    std::cout << "level= " << level << '\n';
    std::cout << "finest_level= " << finest_level << '\n';
    std::cout << "domain= " << domain << '\n';
    std::cout << "ngrid= " << ngrid << '\n';
    std::cout << "mfBA= " << mfBA << '\n';
    amrex::Error("mf contains nan  ::check_for_NAN");
   }
   ParallelDescriptor::Barrier();

   if (mf->contains_inf(sc,1,ng)==true) {
    std::cout << "sc= " << sc << '\n';
    std::cout << "ng= " << ng << '\n';
    std::cout << "ncomp= " << ncomp << '\n';
    std::cout << "ngrow= " << ngrow << '\n';
    std::cout << "level= " << level << '\n';
    std::cout << "finest_level= " << finest_level << '\n';
    std::cout << "domain= " << domain << '\n';
    std::cout << "ngrid= " << ngrid << '\n';
    std::cout << "mfBA= " << mfBA << '\n';
    amrex::Error("mf contains inf ::check_for_NAN");
   }
   ParallelDescriptor::Barrier();
  }  // ng=0..ngrow
 }  // sc=0..ncomp-1

 std::fflush(NULL);

 int datatype=0;
 int force_check=1;
 int dir=-1;
 Real warning_cutoff=1.0e+30;
 aggressive_debug(
  datatype,
  force_check,
  mf,
  scomp,ncomp,
  ngrow,
  dir,
  warning_cutoff);  

 std::fflush(NULL);
} // end subroutine check_for_NAN

//plot_grid_type==0 data interpolated to nodes.
//plot_grid_type==1 data lives at the cells.
void NavierStokes::output_zones(
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
   MultiFab* refine_density_mf,
   MultiFab* lsdistmf,
   MultiFab* viscmf,
   MultiFab* conductmf,
   MultiFab* magtracemf,
   MultiFab* elasticforcemf,
    //ux,vx,wx,uy,vy,wy,uz,vz,wz
   MultiFab* gradvelocitymf,
   int& grids_per_level,
   BoxArray& cgrids_minusBA,
   Real* slice_data,
   int do_plot,int do_slice) {

 std::string local_caller_string="output_zones";

 int plot_sdim_macro=AMREX_SPACEDIM;

 const Real* dx = geom.CellSize();
 const Real* prob_lo = geom.ProbLo();
 const Real* prob_hi = geom.ProbHi();

 int nden=num_materials*num_state_material;

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

  // x,u,p,den,T,Y1..Yn,mag vort,LS
 if (visual_ncomp==VISUALCOMP_NCOMP) {
  // do nothing
 } else
  amrex::Error("visual_ncomp invalid");

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

 for (int data_dir=0;data_dir<AMREX_SPACEDIM;data_dir++) {
  if (localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0)
   amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nGrow()!=0");
  if (localMF[FSI_GHOST_MAC_MF+data_dir]->nComp()!=nparts_def*AMREX_SPACEDIM)
   amrex::Error("localMF[FSI_GHOST_MAC_MF+data_dir]->nComp() bad");
 }

  //VOF_Recon_resize is declared in: NavierStokes2.cpp
 VOF_Recon_resize(1); //output:SLOPE_RECON_MF
 debug_ngrow(SLOPE_RECON_MF,1,local_caller_string);

 debug_ngrow(MASKSEM_MF,1,local_caller_string);

 int finest_level = parent->finestLevel();
 int tecplot_finest_level=finest_level;
 if (tecplot_max_level<tecplot_finest_level)
  tecplot_finest_level=tecplot_max_level;

 NavierStokes& ns_finest=getLevel(tecplot_finest_level);
 const Real* dxfinest = ns_finest.geom.CellSize();

 // x,y,z,xvel,yvel,zvel,PMG,PEOS,div,den,Temp,KE
 // (value of material with LS>0)
 int nslice=0;
 int nstate_slice=SLICECOMP_NCOMP; 

 if ((slice_dir>=0)&&(slice_dir<AMREX_SPACEDIM)) {
  const Box& domain_finest = ns_finest.geom.Domain();
  const int* domlo_finest = domain_finest.loVect();
  const int* domhi_finest = domain_finest.hiVect();
  nslice=domhi_finest[slice_dir]-domlo_finest[slice_dir]+3;
 } else
  amrex::Error("slice_dir invalid");

 int elastic_ncomp=0;
 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=num_materials)) {
  elastic_ncomp=viscoelasticmf->nComp();
 } else if (num_materials_viscoelastic==0) {
  // do nothing
 } else
  amrex::Error("num_materials_viscoelastic invalid");

 if (elastic_ncomp==num_materials_viscoelastic*ENUM_NUM_TENSOR_TYPE_REFINE) {
  // do nothing
 } else
  amrex::Error("elastic_ncomp invalid");

 int refine_density_ncomp=0;
 if ((num_materials_compressible>=1)&&
     (num_materials_compressible<=num_materials)) {
  refine_density_ncomp=refine_density_mf->nComp();
 } else if (num_materials_compressible==0) {
  // do nothing
 } else
  amrex::Error("num_materials_compressible invalid");

 if (refine_density_ncomp==
     num_materials_compressible*ENUM_NUM_REFINE_DENSITY_TYPE) {
  // do nothing
 } else
  amrex::Error("refine_density_ncomp invalid");


 check_for_NAN(localMF[MASKSEM_MF]);
 check_for_NAN(velmf);
 check_for_NAN(localMF[SLOPE_RECON_MF]);
 check_for_NAN(presmf);
 check_for_NAN(divmf);
 check_for_NAN(div_data);
 check_for_NAN(denmf);
 check_for_NAN(mom_denmf);

 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=num_materials)) {
  check_for_NAN(viscoelasticmf);
 } else if (num_materials_viscoelastic==0) {
  // do nothing
 } else
  amrex::Error("num_materials_viscoelastic invalid");

 if ((num_materials_compressible>=1)&&
     (num_materials_compressible<=num_materials)) {
  check_for_NAN(refine_density_mf);
 } else if (num_materials_compressible==0) {
  // do nothing
 } else
  amrex::Error("num_materials_compressible invalid");


 check_for_NAN(lsdistmf);
 check_for_NAN(viscmf);
 check_for_NAN(conductmf);
 check_for_NAN(magtracemf);
 check_for_NAN(elasticforcemf);
    //ux,vx,wx,uy,vy,wy,uz,vz,wz
 check_for_NAN(gradvelocitymf);

 int bfact=parent->Space_blockingFactor(level);

 if (level<=tecplot_finest_level) {

   //plot_grid_type==0 data interpolated to nodes.
  if (plot_grid_type==0) {

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

    MultiFab* maskSEM_minus=new MultiFab(cgrids_minusBA,cgrids_minus_map,1,0,
     MFInfo().SetTag("maskSEM_minus"),FArrayBoxFactory());

    MultiFab* velmfminus=new MultiFab(cgrids_minusBA,cgrids_minus_map,
     STATE_NCOMP_VEL+STATE_NCOMP_PRES,1,
     MFInfo().SetTag("velmfminus"),FArrayBoxFactory());

    MultiFab* vofmfminus=new MultiFab(cgrids_minusBA,cgrids_minus_map,
     num_materials*ngeom_recon,1,
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
     num_materials,1,
     MFInfo().SetTag("mom_denmfminus"),FArrayBoxFactory());

    MultiFab* viscoelasticmfminus;

    if ((num_materials_viscoelastic>=1)&&
        (num_materials_viscoelastic<=num_materials)) {
     viscoelasticmfminus=
      new MultiFab(cgrids_minusBA,cgrids_minus_map,
       elastic_ncomp,1,
       MFInfo().SetTag("viscoelasticmfminus"),FArrayBoxFactory());
    } else if (num_materials_viscoelastic==0) {
     viscoelasticmfminus=mom_denmfminus; //placeholder
    } else
     amrex::Error("num_materials_viscoelastic invalid");

    MultiFab* refine_density_mfminus;

    if ((num_materials_compressible>=1)&&
        (num_materials_compressible<=num_materials)) {
     refine_density_mfminus=
      new MultiFab(cgrids_minusBA,cgrids_minus_map,
       refine_density_ncomp,1,
       MFInfo().SetTag("refine_density_mfminus"),FArrayBoxFactory());
    } else if (num_materials_compressible==0) {
     refine_density_mfminus=mom_denmfminus; //placeholder
    } else
     amrex::Error("num_materials_compressible invalid");

    MultiFab* lsdistmfminus=
     new MultiFab(cgrids_minusBA,cgrids_minus_map,
 	 num_materials*(AMREX_SPACEDIM+1),1,
         MFInfo().SetTag("lsdistmfminus"),FArrayBoxFactory());

    MultiFab* viscmfminus=new MultiFab(cgrids_minusBA,cgrids_minus_map,
     num_materials,1,
     MFInfo().SetTag("viscmfminus"),FArrayBoxFactory());

    MultiFab* conductmfminus=new MultiFab(cgrids_minusBA,cgrids_minus_map,
     num_materials,1,
     MFInfo().SetTag("conductmfminus"),FArrayBoxFactory());

    MultiFab* magtracemfminus=new MultiFab(cgrids_minusBA,cgrids_minus_map,
     5*num_materials,1,
     MFInfo().SetTag("magtracemfminus"),FArrayBoxFactory());

    MultiFab* elasticforcemfminus=new MultiFab(cgrids_minusBA,cgrids_minus_map,
     AMREX_SPACEDIM,1,
     MFInfo().SetTag("elasticforcemfminus"),FArrayBoxFactory());

    //ux,vx,wx,uy,vy,wy,uz,vz,wz
    MultiFab* gradvelocitymfminus=new MultiFab(cgrids_minusBA,cgrids_minus_map,
     AMREX_SPACEDIM_SQR,1,
     MFInfo().SetTag("gradvelocitymfminus"),FArrayBoxFactory());

    MultiFab* towermfminus=new MultiFab(cgrids_minusBA,cgrids_minus_map,
     PLOTCOMP_NCOMP,1,
     MFInfo().SetTag("towermfminus"),FArrayBoxFactory());

    ParallelDescriptor::Barrier();

     // FabArray.H     
     // scomp,dcomp,ncomp,s_nghost,d_nghost
    maskSEM_minus->ParallelCopy(*localMF[MASKSEM_MF],0,0,
     1,0,0,geom.periodicity());

    check_for_NAN(maskSEM_minus);

     // FabArray.H     
     // scomp,dcomp,ncomp,s_nghost,d_nghost
    velmfminus->ParallelCopy(*velmf,0,0,
     STATE_NCOMP_VEL+STATE_NCOMP_PRES,1,1,geom.periodicity());

    check_for_NAN(velmfminus);
 
    vofmfminus->ParallelCopy(*localMF[SLOPE_RECON_MF],0,0,
     num_materials*ngeom_recon,1,1,geom.periodicity());

    check_for_NAN(vofmfminus);

    // scomp,dcomp,ncomp,sgrow,dgrow,period,op
    presmfminus->ParallelCopy(*presmf,0,0,1,
		   1,1,geom.periodicity()); 

    check_for_NAN(presmfminus);

    // scomp,dcomp,ncomp,sgrow,dgrow,period,op
    divmfminus->ParallelCopy(*divmf,0,0,1,
		   1,1,geom.periodicity()); 

    check_for_NAN(divmfminus);

    // scomp,dcomp,ncomp,sgrow,dgrow,period,op
    div_data_minus->ParallelCopy(*div_data,0,0,1,
		   1,1,geom.periodicity()); 

    check_for_NAN(div_data_minus);

    // scomp,dcomp,ncomp,sgrow,dgrow,period,op
    denmfminus->ParallelCopy(*denmf,0,0,nden,
		   1,1,geom.periodicity()); 

    // scomp,dcomp,ncomp,sgrow,dgrow,period,op
    mom_denmfminus->ParallelCopy(*mom_denmf,0,0,num_materials,
		   1,1,geom.periodicity()); 

    check_for_NAN(denmfminus);
    check_for_NAN(mom_denmfminus);

    if ((num_materials_viscoelastic>=1)&&
        (num_materials_viscoelastic<=num_materials)) {
     // scomp,dcomp,ncomp,sgrow,dgrow,period,op
     viscoelasticmfminus->ParallelCopy(*viscoelasticmf,0,0,
      elastic_ncomp,
      1,1,geom.periodicity()); 
     check_for_NAN(viscoelasticmfminus);
    } else if (num_materials_viscoelastic==0) {
     // do nothing
    } else
     amrex::Error("num_materials_viscoelastic invalid");


    if ((num_materials_compressible>=1)&&
        (num_materials_compressible<=num_materials)) {
     // scomp,dcomp,ncomp,sgrow,dgrow,period,op
     refine_density_mfminus->ParallelCopy(*refine_density_mf,0,0,
      refine_density_ncomp,
      1,1,geom.periodicity()); 
     check_for_NAN(refine_density_mfminus);
    } else if (num_materials_compressible==0) {
     // do nothing
    } else
     amrex::Error("num_materials_compressible invalid");

    // scomp,dcomp,ncomp,sgrow,dgrow,period,op
    lsdistmfminus->ParallelCopy(*lsdistmf,0,0,num_materials*(1+AMREX_SPACEDIM),
     1,1,geom.periodicity()); 

    check_for_NAN(lsdistmfminus);

    // scomp,dcomp,ncomp,sgrow,dgrow,period,op
    viscmfminus->ParallelCopy(*viscmf,0,0,num_materials,
		   1,1,geom.periodicity()); 

    check_for_NAN(viscmfminus);

    // scomp,dcomp,ncomp,sgrow,dgrow,period,op
    conductmfminus->ParallelCopy(*conductmf,0,0,num_materials,
                   1,1,geom.periodicity());

    check_for_NAN(conductmfminus);

    // scomp,dcomp,ncomp,sgrow,dgrow,period,op
    magtracemfminus->ParallelCopy(*magtracemf,0,0,5*num_materials,
		   1,1,geom.periodicity()); 

    check_for_NAN(magtracemfminus);
 
    ParallelDescriptor::Barrier();

    // scomp,dcomp,ncomp,sgrow,dgrow,period,op
    elasticforcemfminus->ParallelCopy(*elasticforcemf,0,0,AMREX_SPACEDIM,
		   1,1,geom.periodicity()); 

    check_for_NAN(elasticforcemfminus);
 
    ParallelDescriptor::Barrier();

    // scomp,dcomp,ncomp,sgrow,dgrow,period,op
    //ux,vx,wx,uy,vy,wy,uz,vz,wz
    gradvelocitymfminus->ParallelCopy(*gradvelocitymf,0,0,AMREX_SPACEDIM_SQR,
		   1,1,geom.periodicity()); 

    check_for_NAN(gradvelocitymfminus);
 
    ParallelDescriptor::Barrier();

    towermfminus->setVal(0.0,0,PLOTCOMP_NCOMP,1);

    check_for_NAN(towermfminus);

    ParallelDescriptor::Barrier();

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
     FArrayBox& refine_densityfab=(*refine_density_mfminus)[mfi];
     FArrayBox& lsdistfab=(*lsdistmfminus)[mfi];
     FArrayBox& viscfab=(*viscmfminus)[mfi];
     FArrayBox& conductfab=(*conductmfminus)[mfi];
     FArrayBox& magtracefab=(*magtracemfminus)[mfi];
     FArrayBox& elasticforcefab=(*elasticforcemfminus)[mfi];
      //ux,vx,wx,uy,vy,wy,uz,vz,wz
     FArrayBox& gradvelocityfab=(*gradvelocitymfminus)[mfi];

     FArrayBox& towerfab=(*towermfminus)[mfi];
     int ncomp_tower=towerfab.nComp();
     if (ncomp_tower==PLOTCOMP_NCOMP) {
      // do nothing
     } else
      amrex::Error("ncomp_tower!=PLOTCOMP_NCOMP");

       // declared in: NAVIERSTOKES_3D.F90
     fort_cellgrid(
      &plot_grid_type,
      &ncomp_tower,
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
      ARLIM(elasticfab.loVect()),
      ARLIM(elasticfab.hiVect()),
      refine_densityfab.dataPtr(),
      ARLIM(refine_densityfab.loVect()),
      ARLIM(refine_densityfab.hiVect()),
      lsdistfab.dataPtr(),ARLIM(lsdistfab.loVect()),ARLIM(lsdistfab.hiVect()),
      viscfab.dataPtr(),ARLIM(viscfab.loVect()),ARLIM(viscfab.hiVect()),
      conductfab.dataPtr(),
      ARLIM(conductfab.loVect()),ARLIM(conductfab.hiVect()),
      magtracefab.dataPtr(),
      ARLIM(magtracefab.loVect()),ARLIM(magtracefab.hiVect()),
      elasticforcefab.dataPtr(),
      ARLIM(elasticforcefab.loVect()),ARLIM(elasticforcefab.hiVect()),
      gradvelocityfab.dataPtr(),
      ARLIM(gradvelocityfab.loVect()),ARLIM(gradvelocityfab.hiVect()),
      towerfab.dataPtr(), //towerfab
      ARLIM(towerfab.loVect()),ARLIM(towerfab.hiVect()),
      prob_lo,
      prob_hi,
      dx,
      lo,hi, //tilelo,tilehi
      lo,hi,
      &level,
      &finest_level,
      &gridno,
      &visual_tessellate_vfrac,  // = 0,1,3
      &NS_geometry_coord,
      &nparts,
      &nparts_def,
      im_solid_map_ptr,
      &elastic_ncomp,
      &refine_density_ncomp,
      slice_data,
      &nslice,
      &nstate_slice,&slice_dir,
      xslice.dataPtr(),
      dxfinest,
      &do_plot,
      &do_slice,
      &visual_nddata_format);
    }  // mfi
    ns_reconcile_d_num(LOOP_CELLGRID,"output_zones");

    if ((num_materials_viscoelastic>=1)&&
        (num_materials_viscoelastic<=num_materials)) {
     delete viscoelasticmfminus;
    } else if (num_materials_viscoelastic==0) {
     // do nothing
    } else
     amrex::Error("num_materials_viscoelastic invalid");


    if ((num_materials_compressible>=1)&&
        (num_materials_compressible<=num_materials)) {
     delete refine_density_mfminus;
    } else if (num_materials_compressible==0) {
     // do nothing
    } else
     amrex::Error("num_materials_compressible invalid");

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
    delete conductmfminus;
    delete magtracemfminus;
    delete elasticforcemfminus;
    delete gradvelocitymfminus;

    delete towermfminus;

   } else if (grids_per_level==0) {
  
    // do nothing

   } else 
    amrex::Error("grids_per_level is corrupt");

   //plot_grid_type==1 data lives at the cells.
   //data is output "as is" to "BOXLIB" format.
  } else if (plot_grid_type==1) {

   BoxArray cgrids(grids);
   BoxList cgrids_list(cgrids);

   BoxArray cgrids_minusBA_temp(cgrids_list);
   int grids_per_level_local=cgrids_minusBA_temp.size();
   BoxArray cgrids_minusBA_local;
   cgrids_minusBA_local=cgrids_minusBA_temp;

   if (grids.size()!=grids_per_level_local) {
    amrex::Error("grids_per_level_local incorrect");
   }

   if (grids_per_level_local>0) {

    DistributionMapping cgrids_minus_map(cgrids_minusBA_local);

    ParallelDescriptor::Barrier();

    MultiFab* maskSEM_minus=localMF[MASKSEM_MF];
    check_for_NAN(maskSEM_minus);

    MultiFab* velmfminus=velmf;
    check_for_NAN(velmfminus);
 
    MultiFab* vofmfminus=localMF[SLOPE_RECON_MF];
    check_for_NAN(vofmfminus);

    MultiFab* presmfminus=presmf;
    check_for_NAN(presmfminus);

    MultiFab* divmfminus=divmf;
    check_for_NAN(divmfminus);

    MultiFab* div_data_minus=div_data;
    check_for_NAN(div_data_minus);

    MultiFab* denmfminus=denmf;
    MultiFab* mom_denmfminus=mom_denmf;

    check_for_NAN(denmfminus);
    check_for_NAN(mom_denmfminus);

    MultiFab* viscoelasticmfminus=viscoelasticmf;
    if ((num_materials_viscoelastic>=1)&&
        (num_materials_viscoelastic<=num_materials)) {
     check_for_NAN(viscoelasticmfminus);
    } else if (num_materials_viscoelastic==0) {
     // do nothing
    } else
     amrex::Error("num_materials_viscoelastic invalid:writeTECPLOT_File");


    MultiFab* refine_density_mfminus=refine_density_mf;
    if ((num_materials_compressible>=1)&&
        (num_materials_compressible<=num_materials)) {
     check_for_NAN(refine_density_mfminus);
    } else if (num_materials_compressible==0) {
     // do nothing
    } else
     amrex::Error("num_materials_compressible invalid:writeTECPLOT_File");

    MultiFab* lsdistmfminus=lsdistmf;
    check_for_NAN(lsdistmfminus);

    MultiFab* viscmfminus=viscmf;
    check_for_NAN(viscmfminus);

    MultiFab* conductmfminus=conductmf;
    check_for_NAN(conductmfminus);

    MultiFab* magtracemfminus=magtracemf;
    check_for_NAN(magtracemfminus);
 
    MultiFab* elasticforcemfminus=elasticforcemf;
    check_for_NAN(elasticforcemfminus);

    MultiFab* gradvelocitymfminus=gradvelocitymf;
    check_for_NAN(gradvelocitymfminus);

    MultiFab* towermf=localMF[MULTIFAB_TOWER_PLT_MF];
    check_for_NAN(towermf);
 
    ParallelDescriptor::Barrier();

    if (thread_class::nthreads<1)
     amrex::Error("thread_class::nthreads invalid");
    thread_class::init_d_numPts(velmfminus->boxArray().d_numPts());

    bool use_tiling=ns_tiling;
#ifdef _OPENMP
#pragma omp parallel
#endif
{
    for (MFIter mfi(*velmfminus,use_tiling); mfi.isValid(); ++mfi) {

     if (cgrids_minusBA_local[mfi.index()] != mfi.validbox())
      amrex::Error("cgrids_minusBA_local[mfi.index()] != mfi.validbox()");

     const Box& tilegrid = mfi.tilebox();

     int tid_current=ns_thread();
     if ((tid_current<0)||(tid_current>=thread_class::nthreads))
      amrex::Error("tid_current invalid");
     thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

     const int gridno = mfi.index();

     const int* tilelo=tilegrid.loVect();
     const int* tilehi=tilegrid.hiVect();

     // cgrids_minusBA_local=grids
     const Box& fabgrid = cgrids_minusBA_local[gridno];

     const int* lo=fabgrid.loVect();
     const int* hi=fabgrid.hiVect();

     for (int dir=0;dir<AMREX_SPACEDIM;dir++) {

      if (lo[dir]%bfact!=0)
       amrex::Error("lo not divisible by bfact");
      if ((hi[dir]+1)%bfact!=0)
       amrex::Error("hi+1 not divisible by bfact");
  
     } // dir

     if (cgrids_minusBA_local[gridno]!=grids[gridno])
      amrex::Error("box mismatch1");

     FArrayBox& maskSEMfab=(*maskSEM_minus)[mfi];
     FArrayBox& velfab=(*velmfminus)[mfi];
     FArrayBox& voffab=(*vofmfminus)[mfi];
     FArrayBox& presfab=(*presmfminus)[mfi];
     FArrayBox& divfab=(*divmfminus)[mfi];
     FArrayBox& div_data_fab=(*div_data_minus)[mfi];
     FArrayBox& denfab=(*denmfminus)[mfi];
     FArrayBox& mom_denfab=(*mom_denmfminus)[mfi];
     FArrayBox& elasticfab=(*viscoelasticmfminus)[mfi];
     FArrayBox& refine_densityfab=(*refine_density_mfminus)[mfi];
     FArrayBox& lsdistfab=(*lsdistmfminus)[mfi];
     FArrayBox& viscfab=(*viscmfminus)[mfi];
     FArrayBox& conductfab=(*conductmfminus)[mfi];
     FArrayBox& magtracefab=(*magtracemfminus)[mfi];
     FArrayBox& elasticforcefab=(*elasticforcemfminus)[mfi];
     FArrayBox& gradvelocityfab=(*gradvelocitymfminus)[mfi];

     FArrayBox& towerfab=(*towermf)[mfi];
     int ncomp_tower=towerfab.nComp();
     if (ncomp_tower==PLOTCOMP_NCOMP) {
      // do nothing
     } else
      amrex::Error("ncomp_tower!=PLOTCOMP_NCOMP");

       // declared in: NAVIERSTOKES_3D.F90
     fort_cellgrid(
      &plot_grid_type,
      &ncomp_tower,
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
      ARLIM(elasticfab.loVect()),
      ARLIM(elasticfab.hiVect()),
      refine_densityfab.dataPtr(),
      ARLIM(refine_densityfab.loVect()),
      ARLIM(refine_densityfab.hiVect()),
      lsdistfab.dataPtr(),ARLIM(lsdistfab.loVect()),ARLIM(lsdistfab.hiVect()),
      viscfab.dataPtr(),ARLIM(viscfab.loVect()),ARLIM(viscfab.hiVect()),
      conductfab.dataPtr(),
      ARLIM(conductfab.loVect()),ARLIM(conductfab.hiVect()),
      magtracefab.dataPtr(),
      ARLIM(magtracefab.loVect()),ARLIM(magtracefab.hiVect()),
      elasticforcefab.dataPtr(),
      ARLIM(elasticforcefab.loVect()),ARLIM(elasticforcefab.hiVect()),
      gradvelocityfab.dataPtr(),
      ARLIM(gradvelocityfab.loVect()),ARLIM(gradvelocityfab.hiVect()),
      towerfab.dataPtr(), //towerfab
      ARLIM(towerfab.loVect()),ARLIM(towerfab.hiVect()),
      prob_lo,
      prob_hi,
      dx,
      tilelo,tilehi,
      lo,hi,
      &level,
      &finest_level,
      &gridno,
      &visual_tessellate_vfrac,  // = 0,1,3
      &NS_geometry_coord,
      &nparts,
      &nparts_def,
      im_solid_map_ptr,
      &elastic_ncomp,
      &refine_density_ncomp,
      slice_data,
      &nslice,
      &nstate_slice,&slice_dir,
      xslice.dataPtr(),
      dxfinest,
      &do_plot,
      &do_slice,
      &visual_nddata_format);
    }  // mfi
} //omp
     // declared in: NavierStokes.cpp
    ns_reconcile_d_num(LOOP_CELLGRID_AS_IS,"output_zones");
    ParallelDescriptor::Barrier();

   } else if (grids_per_level_local==0) {
  
    amrex::Error("grids_per_level_local cannot be 0");

   } else 
    amrex::Error("grids_per_level_local is corrupt");

  } else
   amrex::Error("plot_grid_type invalid");

 } else if (level<=finest_level) {

  // do nothing

 } else {
  amrex::Error("level invalid output_zones");
 }

}  // end subroutine output_zones


// data_dir=-1 cell centered data
// data_dir=0..sdim-1 face centered data.
// data_dir=3  X,Y node
// data_dir=4  X,Z node
// data_dir=5  Y,Z node
void NavierStokes::Sanity_output_zones(
   const std::string& information_string,
   int tower_mf_id,
   int data_dir,
   MultiFab* datamf,
   int ncomp,
   int& grids_per_level,
   BoxArray& cgrids_minusBA) {

 std::string local_caller_string="Sanity_output_zones";

 if (tower_mf_id>=0) {
  //do nothing
 } else if (tower_mf_id-GET_NEW_DATA_OFFSET>=0) {
   //do nothing
 } else
  amrex::Error("tower_mf_id out of range");

 const Real* dx = geom.CellSize();
 const Real* prob_lo = geom.ProbLo();
 const Real* prob_hi = geom.ProbHi();

 if (ParallelDescriptor::IOProcessor()) {
  std::cout << "in: Sanity_output_zones information_string=" << 
    information_string <<
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

   int box_type[AMREX_SPACEDIM];
   grid_type_to_box_type_cpp(data_dir,box_type);

   for (int local_dir=0;local_dir<AMREX_SPACEDIM;local_dir++) {
    if (box_type[local_dir]==0) {
     // do nothing
    } else if (box_type[local_dir]==1) {
     minus_boxes.surroundingNodes(local_dir);
    } else
     amrex::Error("box_type[local_dir] invalid");
   }

   MultiFab* datamfminus=new MultiFab(minus_boxes,cgrids_minus_map,
    ncomp,0,
    MFInfo().SetTag("datamfminus"),FArrayBoxFactory());

   ParallelDescriptor::Barrier();

   debug_ixType_raw(datamf,data_dir,local_caller_string);

     // FabArray.H     
     // scomp,dcomp,ncomp,s_nghost,d_nghost
   datamfminus->ParallelCopy(*datamf,0,0,
    ncomp,0,0,geom.periodicity());

   if ((data_dir>=-1)&&(data_dir<=5)) {
    check_for_NAN(datamf);
    check_for_NAN(datamfminus);
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

     // declared in: TECPLOTUTIL.F90
    fort_cellgrid_sanity(
     &tid_current,
     &tower_mf_id,
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
     &NS_geometry_coord);
   }  // mfi
   ns_reconcile_d_num(LOOP_CELLGRID_SANITY,"Sanity_output_zones");

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

}  // end subroutine Sanity_output_zones


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
}  // end subroutine avgDownError_ALL


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
 } // for (int ilev=0;ilev<=finest_level;ilev++) 
} // end subroutine zeroALL

// initializes localMF[idx_localMF] to 0.0
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

} // end subroutine setVal_localMF

int NavierStokes::get_new_data_Type(int mfab_id) {

 if (GET_NEW_DATA_OFFSET<-100) {
  //do nothing
 } else
  amrex::Error("expecting GET_NEW_DATA_OFFSET<-100");

 if (mfab_id>=0) {
  return -1;
 } else if ((mfab_id>=GET_NEW_DATA_OFFSET)&&
   	    (mfab_id<GET_NEW_DATA_OFFSET+NUM_STATE_TYPE)&&
	    (mfab_id<0)) {
  return mfab_id-GET_NEW_DATA_OFFSET;
 } else {
  amrex::Error("mfab_id invalid");
  return GET_NEW_DATA_INVALID;
 }

} // end subroutine get_new_data_Type

void NavierStokes::Copy_array(int mfab_dest,int mfab_source,
	int scomp,int dcomp,int ncomp,int ngrow) {

 int finest_level = parent->finestLevel();

 if ((level==0)&&(finest_level>=0)) {
  // do nothing
 } else
  amrex::Error("expecting level==0 and finest_level>=0");

 for (int ilev=level;ilev<=finest_level;ilev++) {
  NavierStokes& ns_level=getLevel(ilev);
  MultiFab* mfab_dest_ptr=nullptr;
  MultiFab* mfab_source_ptr=nullptr;

  if (get_new_data_Type(mfab_dest)>=0) {
   mfab_dest_ptr=
    &ns_level.get_new_data(get_new_data_Type(mfab_dest),slab_step+1);
  } else if (get_new_data_Type(mfab_dest)==-1) {
   mfab_dest_ptr=ns_level.localMF[mfab_dest];
  } else
   amrex::Error("mfab_dest invalid");

  if (get_new_data_Type(mfab_source)>=0) {
   mfab_source_ptr=
    &ns_level.get_new_data(get_new_data_Type(mfab_source),slab_step+1);
  } else if (get_new_data_Type(mfab_source)==-1) {
   mfab_source_ptr=ns_level.localMF[mfab_source];
  } else
   amrex::Error("mfab_source invalid");

  MultiFab::Copy(*mfab_dest_ptr,*mfab_source_ptr,scomp,dcomp,ncomp,ngrow);
 } // for (int ilev=level;ilev<=finest_level;ilev++) 

} // end subroutine Copy_array

void NavierStokes::setVal_array(int ngrow,int scomp,int ncomp,Real dataval,
  int idx_localMF) {

 int finest_level = parent->finestLevel();

 if (finest_level>=0) {

  if (level==0) {

   for (int i=finest_level;i>=level;i--) {
    NavierStokes& ns_level=getLevel(i);
    MultiFab* mfab_dest_ptr=nullptr;

    if (get_new_data_Type(idx_localMF)>=0) {
     mfab_dest_ptr=
      &ns_level.get_new_data(get_new_data_Type(idx_localMF),slab_step+1);
    } else if (get_new_data_Type(idx_localMF)==-1) {
     mfab_dest_ptr=ns_level.localMF[idx_localMF];
    } else
     amrex::Error("mfab_dest invalid");

    mfab_dest_ptr->setVal(dataval,scomp,ncomp,ngrow);
   }

  } else
   amrex::Error("level invalid");

 } else
  amrex::Error("finest_level invalid");

} // end subroutine setVal_array 


void NavierStokes::mult_array(int ngrow,int ncomp,Real dataval,
  int idx_localMF) {

 int finest_level = parent->finestLevel();
 if (level!=0)
  amrex::Error("level!=0 in mult_array");

 for (int i=finest_level;i>=level;i--) {
  NavierStokes& ns_level=getLevel(i);
  ns_level.localMF[idx_localMF]->mult(dataval,0,ncomp,ngrow);
 }
}  // end subroutine mult_array

//dest=dest-source
void NavierStokes::minusALL(int ngrow,int scomp,int ncomp,
		int idx_dest,int idx_source) {

 int finest_level = parent->finestLevel();
 if (level!=0)
  amrex::Error("level!=0 in minusALL");

 for (int i=finest_level;i>=level;i--) {
  NavierStokes& ns_level=getLevel(i);

  MultiFab* mfab_dest_ptr=nullptr;
  MultiFab* mfab_source_ptr=nullptr;

  if (get_new_data_Type(idx_dest)>=0) {
   mfab_dest_ptr=
    &ns_level.get_new_data(get_new_data_Type(idx_dest),slab_step+1);
  } else if (get_new_data_Type(idx_dest)==-1) {
   mfab_dest_ptr=ns_level.localMF[idx_dest];
  } else
   amrex::Error("idx_dest invalid");

  if (get_new_data_Type(idx_source)>=0) {
   mfab_source_ptr=
    &ns_level.get_new_data(get_new_data_Type(idx_source),slab_step+1);
  } else if (get_new_data_Type(idx_source)==-1) {
   mfab_source_ptr=ns_level.localMF[idx_source];
  } else
   amrex::Error("idx_source invalid");

   // dest=dest-source
   // amrex-master/Src/Base/AMReX_MultiFab.H
  mfab_dest_ptr->minus(*mfab_source_ptr,scomp,ncomp,ngrow);
 }
} // end subroutine minusALL


//dest=dest+source
void NavierStokes::plusALL(int ngrow,int scomp,int ncomp,
		int idx_dest,int idx_source) {

 int finest_level = parent->finestLevel();
 if (level!=0)
  amrex::Error("level!=0 in plusALL");

 for (int i=finest_level;i>=level;i--) {
  NavierStokes& ns_level=getLevel(i);

  MultiFab* mfab_dest_ptr=nullptr;
  MultiFab* mfab_source_ptr=nullptr;

  if (get_new_data_Type(idx_dest)>=0) {
   mfab_dest_ptr=
    &ns_level.get_new_data(get_new_data_Type(idx_dest),slab_step+1);
  } else if (get_new_data_Type(idx_dest)==-1) {
   mfab_dest_ptr=ns_level.localMF[idx_dest];
  } else
   amrex::Error("idx_dest invalid");

  if (get_new_data_Type(idx_source)>=0) {
   mfab_source_ptr=
    &ns_level.get_new_data(get_new_data_Type(idx_source),slab_step+1);
  } else if (get_new_data_Type(idx_source)==-1) {
   mfab_source_ptr=ns_level.localMF[idx_source];
  } else
   amrex::Error("idx_source invalid");

   // dest=dest+source
   // amrex-master/Src/Base/AMReX_MultiFab.H
  mfab_dest_ptr->plus(*mfab_source_ptr,scomp,ncomp,ngrow);
 }
} // end subroutine plusALL



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

}  // end subroutine delete_array

// update_flag=
//  RECON_UPDATE_(NULL|STATE_ERR|STATE_CENTROID|STATE_ERR_AND_CENTROID)
void NavierStokes::VOF_Recon_ALL(
  const std::string& caller_string,
  Real time,
  int update_flag,
  int init_vof_prev_time) {

 if (level!=0)
  amrex::Error("level must be 0 ");

 if (interface_touch_flag==1) {

  std::string local_caller_string="VOF_Recon_ALL";
  local_caller_string=caller_string+local_caller_string;

  int local_update_flag=update_flag;

  if (update_centroid_after_recon==1) {
   // do nothing
  } else if (update_centroid_after_recon==0) {

   if (local_update_flag==RECON_UPDATE_NULL) {
    // do nothing
   } else if (local_update_flag==RECON_UPDATE_STATE_ERR) {
    // do nothing
   } else if (local_update_flag==RECON_UPDATE_STATE_CENTROID) {
    local_update_flag=RECON_UPDATE_NULL;
   } else if (local_update_flag==RECON_UPDATE_STATE_ERR_AND_CENTROID) {
    local_update_flag=RECON_UPDATE_STATE_ERR;
   } else
    amrex::Error("local_update_flag invalid");

  } else
   amrex::Error("update_centroid_after_recon invalid");

  if (verbose>0) {
   if (ParallelDescriptor::IOProcessor()) {
    std::cout << "Start: VOF_Recon_ALL: time= " <<
     time << " local_update_flag= " << local_update_flag << 
     " init_vof_prev_time= " << init_vof_prev_time << '\n';
   }
  }

  int finest_level=parent->finestLevel();

#if (NS_profile_solver==1)
  BLProfiler bprof(local_caller_string);
#endif

  int number_centroid=0;
  Real delta_centroid=0.0;

  for (int ilev=level;ilev<=finest_level;ilev++) {
   NavierStokes& ns_level=getLevel(ilev);
   ns_level.delete_localMF_if_exist(SLOPE_RECON_MF,1);
   int ngrow=1;
    // sets values to 0.0
   ns_level.new_localMF(SLOPE_RECON_MF,num_materials*ngeom_recon,ngrow,-1);  

   ns_level.delete_localMF_if_exist(VOF_RECON_MF,1);

   ns_level.getState_localMF(VOF_RECON_MF,2,STATECOMP_MOF,
     num_materials*ngeom_raw,time);

   if (ngeom_raw==AMREX_SPACEDIM+1) {
    //do nothing
   } else
    amrex::Error("ngeom_raw invalid");

   if (ngeom_recon==2*AMREX_SPACEDIM+3) {
    //do nothing
   } else
    amrex::Error("ngeom_recon invalid");

    //SLOPE_RECON_MF: num_materials x (vof,cenref,order,slope,intercept)
    //VOF_RECON_MF  : num_materials x (vof,cenref)
   for (int im=0;im<num_materials;im++) {
    int ibase_raw=im*ngeom_raw;
    int ibase_recon=im*ngeom_recon;
    ns_level.Copy_localMF(SLOPE_RECON_MF,VOF_RECON_MF,
      ibase_raw,ibase_recon,ngeom_raw,1);
    if (init_vof_prev_time==1) {
     int ngrow_save=ns_level.localMF[VOF_PREV_TIME_MF]->nGrow();
     if (ngrow_save!=2)
      amrex::Error("vof prev time has invalid ngrow");
     ns_level.Copy_localMF(VOF_PREV_TIME_MF,VOF_RECON_MF,
	ibase_raw,im,1,ngrow_save); 
    } else if (init_vof_prev_time==0) {
     // do nothing
    } else
     amrex::Error("init_vof_prev_time invalid");
   } // im=0..num_materials-1
  } // for (int ilev=level;ilev<=finest_level;ilev++)

#ifdef AMREX_PARTICLES
  int save_slab_step=slab_step;
  if (slab_step==-1) {
   slab_step=0;
  } if ((slab_step>=0)&&(slab_step<ns_time_order)) {
   //do nothing
  } else if (slab_step==ns_time_order) {
   slab_step=ns_time_order-1;
  } else
   amrex::Error("slab_step invalid");

   //calling from NavierStokes::VOF_Recon_ALL
  init_particle_containerALL(OP_PARTICLE_SLOPES,local_caller_string);

  slab_step=save_slab_step;
#endif

  // go from coarsest to finest so that SLOPE_RECON_MF
  // can have proper BC.
  for (int ilev=level;ilev<=finest_level;ilev++) {
   NavierStokes& ns_level=getLevel(ilev);
   int number_centroid_level=0;
   Real delta_centroid_level=0.0;
   ns_level.VOF_Recon(
    time,
    local_update_flag,
    init_vof_prev_time,
    delta_centroid_level,
    number_centroid_level);
   delta_centroid+=delta_centroid_level;
   number_centroid+=number_centroid_level;
  } // for (int ilev=level;ilev<=finest_level;ilev++) 

  Real single_centroid_diff=0.0;
  if (number_centroid==0) {
   //do nothing
  } else if (number_centroid>0) {
   single_centroid_diff=delta_centroid/number_centroid;
  } else
   amrex::Error("number_centroid invalid");

  if (local_update_flag==RECON_UPDATE_NULL) {

   //do nothing

  } else if (local_update_flag==RECON_UPDATE_STATE_ERR) {

   avgDownError_ALL(); //updates the new data.

  } else if (local_update_flag==RECON_UPDATE_STATE_CENTROID) {

   for (int ilev=finest_level;ilev>=level;ilev--) {
    NavierStokes& ns_level=getLevel(ilev);
    if (ilev<finest_level) {
     ns_level.MOFavgDown();
    }
   } // ilev=finest_level ... level

  } else if (local_update_flag==RECON_UPDATE_STATE_ERR_AND_CENTROID) {

   avgDownError_ALL(); //updates the new data.
		      
   for (int ilev=finest_level;ilev>=level;ilev--) {
    NavierStokes& ns_level=getLevel(ilev);
    if (ilev<finest_level) {
     ns_level.MOFavgDown();
    }
   } // ilev=finest_level ... level

  } else
   amrex::Error("local_update_flag invalid");

  if (verbose>0) {
   if (ParallelDescriptor::IOProcessor()) {
    std::cout << "continuous_mof= " << continuous_mof << '\n';
    std::cout << "number_centroid= " << number_centroid << '\n';
    std::cout << "single_centroid_diff= " << single_centroid_diff << '\n';
   } //IOProc?
  } else if (verbose==0) {
   //do nothing
  } else
   amrex::Error("verbose invalid");

#if (NS_profile_solver==1)
  bprof.stop();
#endif

 } else if (interface_touch_flag==0) {
  //do nothing
 } else
  amrex::Error("interface_touch_flag invalid");

 interface_touch_flag=0;

} // end subroutine VOF_Recon_ALL

void NavierStokes::VOF_Recon_resize(int ngrow) {

 if ((ngrow<1)||(ngrow>ngrow_distance))
  amrex::Error("ngrow invalid");
 if (localMF[SLOPE_RECON_MF]->nComp()!=num_materials*ngeom_recon)
  amrex::Error("localMF[SLOPE_RECON_MF]->nComp()!=num_materials*ngeom_recon");

 if (localMF[SLOPE_RECON_MF]->nGrow()==ngrow) {
  // do nothing
 } else if (localMF[SLOPE_RECON_MF]->nGrow()>=0) {
  MultiFab* slopes_mf=new MultiFab(grids,dmap,num_materials*ngeom_recon,0,
	  MFInfo().SetTag("slope_m"),FArrayBoxFactory());
  MultiFab::Copy(*slopes_mf,*localMF[SLOPE_RECON_MF],0,0,
		 num_materials*ngeom_recon,0);
  delete_localMF(SLOPE_RECON_MF,1);
  //sets values to 0.0
  new_localMF(SLOPE_RECON_MF,num_materials*ngeom_recon,ngrow,-1); 
  //sets values to 0.0
  MultiFab::Copy(*localMF[SLOPE_RECON_MF],*slopes_mf,
		 0,0,num_materials*ngeom_recon,0);
 
  Vector<int> scompBC_map;
  scompBC_map.resize(num_materials*ngeom_recon);
  for (int i=0;i<num_materials*ngeom_recon;i++)
   scompBC_map[i]=i+1+AMREX_SPACEDIM;
   //scomp=0
  PCINTERP_fill_borders(SLOPE_RECON_MF,ngrow,0,num_materials*ngeom_recon,
    State_Type,scompBC_map);
 
  delete slopes_mf;
 } else
  amrex::Error("localMF[SLOPE_RECON_MF]->nGrow() invalid");

} // end subroutine VOF_Recon_resize

// VOF_Recon() is called from VOF_Recon_ALL()
// vof,ref centroid,order,slope,intercept  x num_materials
// update_flag=
//  RECON_UPDATE_(NULL|STATE_ERR|STATE_CENTROID|STATE_ERR_AND_CENTROID)
// 1. get MOF data with 1 ghost cell (so that CMOF can be chosen)
// 2. reconstruct interior cells only.
// 3. do extended filpatch; MOF used for coarse/fine and ext_dir cells.
void NavierStokes::VOF_Recon(Real time,
  int update_flag,int init_vof_prev_time,
  Real& delta_centroid_level,
  int& number_centroid_level) {

 std::string local_caller_string="VOF_Recon";
 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();
 int max_level = parent->maxLevel();
 int nsteps=parent->levelSteps(0);

 if (ngrow_distance!=4)
  amrex::Error("ngrow_distance invalid");

 int bfact=parent->Space_blockingFactor(level);

 if ((verbose>0)&&(1==0)) {
  if (ParallelDescriptor::IOProcessor()) {
   std::cout << "beginning of VOF_Recon: ngrow_distance,level: " << 
    ngrow_distance << ' ' << level << '\n';
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

 Vector< Real > delta_centroid_per_core;
 Vector< int > number_centroid_per_core;
 delta_centroid_per_core.resize(thread_class::nthreads);
 number_centroid_per_core.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  delta_centroid_per_core[tid]=0.0;
  number_centroid_per_core[tid]=0;
 }

 Vector< Vector<int> > total_calls;
 Vector< Vector<int> > total_iterations;
 Vector< Vector<Real> > total_errors;
 total_calls.resize(thread_class::nthreads);
 total_iterations.resize(thread_class::nthreads);
 total_errors.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  total_calls[tid].resize(num_materials);
  total_iterations[tid].resize(num_materials);
  total_errors[tid].resize(num_materials);
  for (int im=0;im<num_materials;im++) {
   total_calls[tid][im]=0;
   total_iterations[tid][im]=0;
   total_errors[tid][im]=0.0;
  }
 } // tid=0..thread_class::nthreads-1

 double start_recon = ParallelDescriptor::second();

 MultiFab &S_new = get_new_data(State_Type,slab_step+1);

 MultiFab* lsdata=getStateDist(1,time,local_caller_string);
 if (lsdata->nComp()!=num_materials*(1+AMREX_SPACEDIM))
  amrex::Error("lsdata invalid ncomp");

 debug_ngrow(SLOPE_RECON_MF,0,local_caller_string);
 if (localMF[SLOPE_RECON_MF]->nComp()!=num_materials*ngeom_recon)
  amrex::Error("invalid ncomp for SLOPE_RECON_MF");

 debug_ngrow(VOF_RECON_MF,0,local_caller_string);
 if (localMF[VOF_RECON_MF]->nComp()!=num_materials*ngeom_raw)
  amrex::Error("invalid ncomp for VOF_RECON_MF");

 resize_mask_nbr(1);
 debug_ngrow(MASK_NBR_MF,1,local_caller_string);
 if (localMF[MASK_NBR_MF]->nComp()!=4)
  amrex::Error("invalid ncomp for mask nbr");

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,1,local_caller_string);

 const Real* dx = geom.CellSize();

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

   FArrayBox& vfab=(*localMF[VOF_RECON_MF])[mfi];
   FArrayBox& lsfab=(*lsdata)[mfi];
   FArrayBox& mffab=(*localMF[SLOPE_RECON_MF])[mfi];

   FArrayBox& snewfab=S_new[mfi];

   // mask_nbr:
   // (1) =1 interior  =1 fine-fine ghost in domain  =0 otherwise
   // (2) =1 interior  =0 otherwise
   // (3) =1 interior+ngrow-1  =0 otherwise
   // (4) =1 interior+ngrow    =0 otherwise
   FArrayBox& masknbr=(*localMF[MASK_NBR_MF])[mfi];

    // maskcov=tag if not covered by level+1 or outside the domain.
   FArrayBox& maskcov=(*localMF[MASKCOEF_MF])[mfi];

   Vector<int> vofbc=getBCArray(State_Type,gridno,STATECOMP_MOF,1);

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
    vofbc.dataPtr(),
    tilelo,tilehi,
    fablo,fabhi,&bfact,
    xlo,dx,
    maskcov.dataPtr(),
    ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
    masknbr.dataPtr(),
    ARLIM(masknbr.loVect()),ARLIM(masknbr.hiVect()),
    snewfab.dataPtr(STATECOMP_MOF),
    ARLIM(snewfab.loVect()),ARLIM(snewfab.hiVect()),
    vfab.dataPtr(),ARLIM(vfab.loVect()),ARLIM(vfab.hiVect()),
    lsfab.dataPtr(),ARLIM(lsfab.loVect()),ARLIM(lsfab.hiVect()),
    mffab.dataPtr(),ARLIM(mffab.loVect()),ARLIM(mffab.hiVect()),
    &nsteps,
    &time,
    &update_flag,
    &number_centroid_per_core[tid_current],
    &delta_centroid_per_core[tid_current],
    total_calls[tid_current].dataPtr(),
    total_iterations[tid_current].dataPtr(),
    total_errors[tid_current].dataPtr(),
    &continuous_mof,  //fort_sloperecon
    &partial_cmof_stencil_at_walls);
 }  // mfi
} // omp
 ns_reconcile_d_num(LOOP_SLOPE_RECON,"VOF_Recon");

 for (int tid=1;tid<thread_class::nthreads;tid++) {
  number_centroid_per_core[0]+=number_centroid_per_core[tid];
  delta_centroid_per_core[0]+=delta_centroid_per_core[tid];
  for (int im=0;im<num_materials;im++) {
   total_calls[0][im]+=total_calls[tid][im];
   total_iterations[0][im]+=total_iterations[tid][im];
   total_errors[0][im]+=total_errors[tid][im];
  }
 } // tid

 ParallelDescriptor::ReduceIntSum(number_centroid_per_core[0]);
 ParallelDescriptor::ReduceRealSum(delta_centroid_per_core[0]);

 number_centroid_level=number_centroid_per_core[0];
 delta_centroid_level=delta_centroid_per_core[0];

 for (int im=0;im<num_materials;im++) {
   ParallelDescriptor::ReduceIntSum(total_calls[0][im]);
   ParallelDescriptor::ReduceIntSum(total_iterations[0][im]);
   ParallelDescriptor::ReduceRealSum(total_errors[0][im]);
 }

 Vector<int> scompBC_map;
 scompBC_map.resize(num_materials*ngeom_recon);
 for (int i=0;i<num_materials*ngeom_recon;i++)
  scompBC_map[i]=i+1+AMREX_SPACEDIM;

 if (1==0) {
  for (int im=0;im<num_materials;im++) {
   int vofcomp=im*ngeom_recon;
   std::cout << "lev,im,calls,iter,err " << level << ' ' << 
    im << ' ' << total_calls[0][im] << ' ' <<
    total_iterations[0][im] << ' ' <<
    total_errors[0][im] << '\n';
   Real vfracnrm=localMF[SLOPE_RECON_MF]->norm1(vofcomp,0);
   std::cout << "im= " << im << " vfracnrm= " << vfracnrm << '\n';
  }
 }

  //scomp=0
 int ngrow=1;
 PCINTERP_fill_borders(SLOPE_RECON_MF,ngrow,0,num_materials*ngeom_recon,
   State_Type,scompBC_map);

 double end_recon = ParallelDescriptor::second();
 double cputime=end_recon-start_recon;

 if (verbose>0) {
  if (ParallelDescriptor::IOProcessor()) {
   for (int im=0;im<num_materials;im++) {
    Real iter_per_call=0.0;
    Real error_per_call=0.0;
    Real r_calls=total_calls[0][im];
    Real r_iters=total_iterations[0][im];
    Real r_error=total_errors[0][im];
    if (total_calls[0][im]>0) {
     iter_per_call=r_iters/r_calls;
     error_per_call=r_error/r_calls;
    }
    std::cout << 
     "lev,im,calls,iter,iter/call,error x refvfrac/call,cpu time " << 
     level << ' ' << 
     im << ' ' << total_calls[0][im] << ' ' <<
     total_iterations[0][im] << ' ' << iter_per_call << ' ' << 
     error_per_call << ' ' << cputime << '\n';
   } // im=0..num_materials-1
  }
 } else if (verbose<=0) {
  // do nothing
 } else
  amrex::Error("verbose invalid");

 delete lsdata;

}  // end subroutine VOF_Recon



void NavierStokes::MOF_training() {

 int finest_level=parent->finestLevel();
 int max_level = parent->maxLevel();

 if (finest_level<=max_level) {
  //do nothing
 } else
  amrex::Error("expecting finest_level<=max_level");

 if (level==0) {
  //do nothing
 } else
  amrex::Error("expecting MOF_training to be called from level=0");

 int bfact=parent->Space_blockingFactor(max_level);

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 const Real* dx = geom.CellSize();
 const Box& domain = geom.Domain();

 Box domain_max_level = domain;
 Real dx_max_level[AMREX_SPACEDIM];

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  dx_max_level[dir]=dx[dir];
 }
 for (int lev=1;lev<=max_level;lev++) {
  domain_max_level.refine(2);
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   dx_max_level[dir]=dx_max_level[dir]/2.0;
  }
 }
 const int* domlo = domain_max_level.loVect();
 const int* domhi = domain_max_level.hiVect();

 int cpp_training_lo[AMREX_SPACEDIM];
 int cpp_training_hi[AMREX_SPACEDIM];
 int op_training=0;  // 0=alloc
 int cpp_i=0;
 int cpp_j=0;
 int cpp_k=0;
 int local_continuous_mof=STANDARD_MOF;

 int tid_current=ns_thread();
 if ((tid_current<0)||(tid_current>=thread_class::nthreads))
  amrex::Error("tid_current invalid");

  // fort_MOF_DT_training is declared in: PLIC_3D.F90
  // "local_continuous_mof" varied internally.
 fort_MOF_DT_training(
   &tid_current,
   &mof_decision_tree_learning,
   &finest_level,
   &max_level,
   &bfact,
   domlo,domhi,
   dx_max_level);

 if (mof_machine_learning>0) {

  // fort_MOF_training is declared in: PLIC_3D.F90
  // op_training=0 => alloc
  // op_training=1 => generate data and do python processing
  // op_training=2 => read either NN, DT, or RF network data
  fort_MOF_training(
   &tid_current,
   &mof_machine_learning,
   &op_training, // =0 ("allocate")
   cpp_training_lo,
   cpp_training_hi,
   &cpp_i,&cpp_j,&cpp_k,
   &finest_level,
   &max_level,
   &bfact,
   domlo,domhi,
   dx_max_level,
   &local_continuous_mof); // only used if op_training=1,2

  ParallelDescriptor::Barrier();

  for (cpp_i=cpp_training_lo[0];cpp_i<=cpp_training_hi[0];cpp_i++) {
  for (cpp_j=cpp_training_lo[1];cpp_j<=cpp_training_hi[1];cpp_j++) {
#if (AMREX_SPACEDIM==3)
  for (cpp_k=cpp_training_lo[2];cpp_k<=cpp_training_hi[2];cpp_k++) {
#endif
  for (local_continuous_mof=STANDARD_MOF;
       local_continuous_mof<=CMOF_X; 
       local_continuous_mof++) {
   op_training=1;  // generate data and do python processing.
   ParallelDescriptor::Barrier();
   if (ParallelDescriptor::IOProcessor()) {
    fort_MOF_training(
     &tid_current,
     &mof_machine_learning,
     &op_training,
     cpp_training_lo,
     cpp_training_hi,
     &cpp_i,&cpp_j,&cpp_k,
     &finest_level,
     &max_level,
     &bfact,
     domlo,domhi,
     dx_max_level,
     &local_continuous_mof);
   }
   ParallelDescriptor::Barrier();
   op_training=2;  // read either NN, DT, or RF network data
   fort_MOF_training(
    &tid_current,
    &mof_machine_learning,
    &op_training,
    cpp_training_lo,
    cpp_training_hi,
    &cpp_i,&cpp_j,&cpp_k,
    &finest_level,
    &max_level,
    &bfact,
    domlo,domhi,
    dx_max_level,
    &local_continuous_mof);
   ParallelDescriptor::Barrier();
  } //local_continuous_mof
#if (AMREX_SPACEDIM==3)
  } //cpp_k
#endif
  } //cpp_j
  } //cpp_i

 } else if (mof_machine_learning==0) {
  // do nothing
 } else
  amrex::Error("mof_machine_learning invalid");

}  // end subroutine MOF_training



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

} // end subroutine build_masksemALL

void NavierStokes::build_masksem(int mask_sweep) {

 std::string local_caller_string="build_masksem";

 bool use_tiling=ns_tiling;

 int finest_level=parent->finestLevel();
 int bfact=parent->Space_blockingFactor(level);
 int bfact_fine=bfact;
 if (level<finest_level) 
  bfact_fine=parent->Space_blockingFactor(level+1);
 if (bfact_fine>bfact)
  amrex::Error("bfact_fine invalid");

 const Box& domain = geom.Domain();
 const int* domlo = domain.loVect();
 const int* domhi = domain.hiVect();

 const Real* dx = geom.CellSize();

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,1,local_caller_string);
 resize_mask_nbr(1);
 debug_ngrow(MASK_NBR_MF,1,local_caller_string);
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

  debug_ngrow(MASKSEM_MF,1,local_caller_string);
  old_mask=new MultiFab(grids,dmap,1,1,
	MFInfo().SetTag("old_mask"),FArrayBoxFactory());
  MultiFab::Copy(*old_mask,*localMF[MASKSEM_MF],0,0,1,1);

 } else
  amrex::Error("mask_sweep invalid");

 Real total_cells_level=localMF[MASKCOEF_MF]->norm1();

 Vector< Vector< Real >  > spectral_cells_level;
 spectral_cells_level.resize(thread_class::nthreads);
 for (int tid=0;tid<thread_class::nthreads;tid++) {
  spectral_cells_level[tid].resize(num_materials);
  for (int im=0;im<num_materials;im++)
   spectral_cells_level[tid][im]=0.0;
 }

 MultiFab* vofmat=new MultiFab(grids,dmap,num_materials,0,
  MFInfo().SetTag("vofmat"),FArrayBoxFactory());
 MultiFab& S_new=get_new_data(State_Type,slab_step+1);

 for (int im=0;im<num_materials;im++) {
  int scomp=STATECOMP_MOF+im*ngeom_raw;
  MultiFab::Copy(*vofmat,S_new,scomp,im,1,0);
 }

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

  Vector<int> vofbc=getBCArray(State_Type,gridno,STATECOMP_MOF,1);

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   // declared in: GODUNOV_3D.F90
  fort_build_masksem( 
   dx,
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
   &bfact_fine);
 }  // mfi
} // omp
 ns_reconcile_d_num(LOOP_BUILD_MASKSEM,"build_masksem");

 for (int tid=1;tid<thread_class::nthreads;tid++) {
  for (int im=0;im<num_materials;im++)
   spectral_cells_level[0][im]+=spectral_cells_level[tid][im];
 }
 for (int im=0;im<num_materials;im++)
  ParallelDescriptor::ReduceRealSum(spectral_cells_level[0][im]);

 Vector<int> scompBC_map;
 scompBC_map.resize(1);
 scompBC_map[0]=1+AMREX_SPACEDIM+num_materials*ngeom_recon;
  // std::string maskextrap_str="maskSEMextrap"
  // fort_extrapfill
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
    for (int im=0;im<num_materials;im++) {
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
  int scomp=0;
  int ncomp=1;
  std::cout << " checking out masksem_mf \n";
  tecplot_debug(mffab,xlo,fablo,fabhi,dx,-1,0,scomp,ncomp,interior_only);
 }

} // emd subroutine build_masksem()



// If incompressible material then copy existing state pressure.
// If compressible material then P=P(rho,e).
MultiFab* NavierStokes::derive_EOS_pressure(Vector<int> local_material_type) {

 std::string local_caller_string="derive_EOS_pressure";
 int finest_level=parent->finestLevel();
 
 bool use_tiling=ns_tiling;

 debug_ngrow(LEVELPC_MF,1,local_caller_string);
 if (localMF[LEVELPC_MF]->nComp()!=num_materials*(AMREX_SPACEDIM+1)) {
  std::cout << "in: derive_EOS_pressure\n";
  amrex::Error("localMF[LEVELPC_MF]->nComp()!=num_materials*(AMREX_SPACEDIM+1)");
 }
 resize_levelset(2,LEVELPC_MF);
 debug_ngrow(LEVELPC_MF,2,local_caller_string);
 if (localMF[LEVELPC_MF]->nComp()!=num_materials*(AMREX_SPACEDIM+1))
  amrex::Error("localMF[LEVELPC_MF]->nComp()!=num_materials*(AMREX_SPACEDIM+1)");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 MultiFab* mf=getStatePres(1,cur_time_slab); 

 if (mf->nComp()!=1)
  amrex::Error("mf->nComp() invalid");

 MultiFab* denmf=getStateDen(1,cur_time_slab);  // num_materials * den,temp,...
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
 
  FArrayBox& levelpcfab=(*localMF[LEVELPC_MF])[mfi];
  FArrayBox& denfab=(*denmf)[mfi];
  FArrayBox& presfab=(*mf)[mfi];

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();
  
   // declared in: NAVIERSTOKES_3D.F90 
  fort_eos_pressure(
   &level,
   &finest_level,
   local_material_type.dataPtr(),
   xlo,dx,
   presfab.dataPtr(),
   ARLIM(presfab.loVect()),ARLIM(presfab.hiVect()),
   levelpcfab.dataPtr(),
   ARLIM(levelpcfab.loVect()),ARLIM(levelpcfab.hiVect()),
   denfab.dataPtr(),ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
   tilelo,tilehi,
   fablo,fabhi,
   &bfact,
   &nden);
 
 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_EOS_PRESSURE,"derive_EOS_pressure");
 
 delete denmf;

 return mf;
} // subroutine derive_EOS_pressure()

void NavierStokes::level_getshear(
	MultiFab* shear_output_mf,
	MultiFab* vel_mf,
	int only_scalar,
	int destcomp,
	int ngrow) {

 std::string local_caller_string="level_getshear";

 int finest_level = parent->finestLevel();
 if ((level>=0)&&(level<=finest_level)) {
  // do nothing
 } else
  amrex::Error("level invalid");

 if ((ngrow==0)||(ngrow==1)) {
  // do nothing
 } else
  amrex::Error("ngrow invalid");

 bool use_tiling=ns_tiling;

 debug_ngrow(CELLTENSOR_MF,1,local_caller_string);
 if (localMF[CELLTENSOR_MF]->nComp()!=AMREX_SPACEDIM_SQR)
  amrex::Error("localMF[CELLTENSOR_MF]->nComp() invalid");

 if (vel_mf->nGrow()>=ngrow) {
  // do nothing
 } else
  amrex::Error("vel_mf->nGrow() invalid");

 if (vel_mf->nComp()>=STATE_NCOMP_VEL) {
  // do nothing
 } else
  amrex::Error("vel_mf->nComp() invalid");

 if (shear_output_mf->nGrow()>=ngrow) {
  // do nothing
 } else
  amrex::Error("shear_output_mf->nGrow() invalid");

 if (shear_output_mf->nComp()>=destcomp+1) {
  // do nothing
 } else
  amrex::Error("shear_output_mf->nComp() invalid");

 const Real* dx = geom.CellSize();

 if (thread_class::nthreads<1)
  amrex::Error("thread_class::nthreads invalid");
 thread_class::init_d_numPts(vel_mf->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
 for (MFIter mfi(*vel_mf,use_tiling); mfi.isValid(); ++mfi) {
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

  FArrayBox& shear_output_fab=(*shear_output_mf)[mfi];

  FArrayBox& cellten=(*localMF[CELLTENSOR_MF])[mfi];
  if (cellten.nComp()!=AMREX_SPACEDIM_SQR)
   amrex::Error("cellten invalid ncomp");

  Vector<int> velbc=getBCArray(State_Type,gridno,
    STATECOMP_VEL,STATE_NCOMP_VEL);

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   //fort_getshear is declared in: DERIVE_3D.F90
  fort_getshear(
    cellten.dataPtr(),
    ARLIM(cellten.loVect()),ARLIM(cellten.hiVect()),
    (*vel_mf)[mfi].dataPtr(),
    ARLIM((*vel_mf)[mfi].loVect()),ARLIM((*vel_mf)[mfi].hiVect()),
    dx,xlo,
    shear_output_fab.dataPtr(destcomp),
    ARLIM(shear_output_fab.loVect()),ARLIM(shear_output_fab.hiVect()),
    &only_scalar,
    &cur_time_slab,
    tilelo,tilehi,
    fablo,fabhi,
    &bfact,
    &level,
    velbc.dataPtr(),
    &ngrow);

 } //mfi
} // omp
 ns_reconcile_d_num(LOOP_GETSHEAR,"level_getshear");

} // end subroutine level_getshear

// in NavierStokes::multiphase_project when:
// project_option==SOLVETYPE_PRES 
//  and homflag=0,
//  the following commands are given:
//  for ilev=finest ... coarsest,
//   ns_level.init_pressure_error_indicator();
//  avgDownError_ALL();
// during regridding, the following routine checks the error:
//  NavierStokes::errorEst  (calls fort_vfracerror)
void NavierStokes::init_pressure_error_indicator() {

 std::string local_caller_string="init_pressure_error_indicator";

 bool use_tiling=ns_tiling;

 int finest_level = parent->finestLevel();

 resize_maskfiner(1,MASKCOEF_MF);
 debug_ngrow(MASKCOEF_MF,1,local_caller_string);

 debug_ngrow(CELLTENSOR_MF,1,local_caller_string);
 if (localMF[CELLTENSOR_MF]->nComp()!=AMREX_SPACEDIM_SQR)
  amrex::Error("localMF[CELLTENSOR_MF]->nComp() invalid");

 resize_levelset(2,LEVELPC_MF);
 debug_ngrow(LEVELPC_MF,2,local_caller_string);

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 MultiFab* denmf=getStateDen(1,cur_time_slab); 
 if (denmf->nComp()!=num_materials*num_state_material)
  amrex::Error("denmf incorrect ncomp");

 MultiFab* presmf=getState(1,AMREX_SPACEDIM,1,cur_time_slab);

 MultiFab* vortmf=new MultiFab(grids,dmap,1,0,
	MFInfo().SetTag("vortmf"),FArrayBoxFactory());
 const Real* dx = geom.CellSize();
 MultiFab& S_new=get_new_data(State_Type,slab_step+1);
 int scomp_error=STATECOMP_ERR;
 if (STATE_NCOMP!=S_new.nComp())
  amrex::Error("STATE_NCOMP!=S_new.nComp()");

 MultiFab* velmf=getState(1,STATECOMP_VEL,STATE_NCOMP_VEL,cur_time_slab);

 int only_scalar=2; // magnitude of vorticity
 int destcomp=0;
 int ngrow_zero=0;
 level_getshear(vortmf,velmf,only_scalar,destcomp,ngrow_zero);

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

  FArrayBox& LSfab=(*localMF[LEVELPC_MF])[mfi];
  FArrayBox& denfab=(*denmf)[mfi];
  FArrayBox& presfab=(*presmf)[mfi];
  FArrayBox& vortfab=(*vortmf)[mfi];
  FArrayBox& errnew=S_new[mfi];

  if (denfab.nComp()!=num_materials*num_state_material)
   amrex::Error("denfab.nComp() invalid");
  if (LSfab.nComp()!=num_materials*(AMREX_SPACEDIM+1))
   amrex::Error("LSfab.nComp() invalid");
  if (presfab.nComp()!=1)
   amrex::Error("presfab.nComp() invalid");
  if (vortfab.nComp()!=1)
   amrex::Error("vortfab.nComp() invalid");

  int tid_current=ns_thread();
  if ((tid_current<0)||(tid_current>=thread_class::nthreads))
   amrex::Error("tid_current invalid");
  thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   //fort_pressure_indicator is declared in: NAVIERSTOKES_3D.F90
  fort_pressure_indicator(
   &pressure_error_flag,
   vorterr.dataPtr(),
   pressure_error_cutoff.dataPtr(),
   temperature_error_cutoff.dataPtr(),
   xlo,dx,
   errnew.dataPtr(scomp_error),
   ARLIM(errnew.loVect()),ARLIM(errnew.hiVect()),
   LSfab.dataPtr(),ARLIM(LSfab.loVect()),ARLIM(LSfab.hiVect()),
   denfab.dataPtr(),ARLIM(denfab.loVect()),ARLIM(denfab.hiVect()),
   vortfab.dataPtr(),ARLIM(vortfab.loVect()),ARLIM(vortfab.hiVect()),
   presfab.dataPtr(),ARLIM(presfab.loVect()),ARLIM(presfab.hiVect()),
   maskcov.dataPtr(),ARLIM(maskcov.loVect()),ARLIM(maskcov.hiVect()),
   tilelo,tilehi,
   fablo,fabhi, 
   &bfact,
   &level,
   &finest_level);

 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_PRESSURE_INDICATOR,"init_pressure_error_indicator");

 delete presmf;
 delete vortmf;
 delete denmf;

} // subroutine init_pressure_error_indicator

// if project_option==SOLVETYPE_PRES:
// 1. calculates p(rho^n+1,e_advect) and puts it in 2nd component
//    of cell_sound.
//    Equation of state to be used depends on vofPC
// 2. calculates 1/(rho c^2 dt^2)  and puts it in 1st component of cell_sound.
//    Equation of state to be used depends on vofPC
// 3. in incompressible regions, p=0
//
// div_hold=(pnew-pold)/(rho c^2 dt) + dt mdot/vol
//
// velocity scale: V
// time scale is : 1/V
// pressure scale: V^2
// scale for "cell_sound" is 1
//
// init_advective_pressure is called from: NavierStokes::multiphase_project
void NavierStokes::init_advective_pressure(int project_option) {
 
 std::string local_caller_string="init_advective_pressure";

 int finest_level=parent->finestLevel();
 if ((level>=0)&&(level<=finest_level)) {
  // do nothing
 } else
  amrex::Error("level invalid");

 bool use_tiling=ns_tiling;

 resize_metrics(1);
 debug_ngrow(VOLUME_MF,0,local_caller_string);

 debug_ngrow(FACE_VAR_MF,0,local_caller_string);

 debug_ngrow(MASKCOEF_MF,1,local_caller_string);

 debug_ngrow(DIFFUSIONRHS_MF,0,local_caller_string);

 if (localMF[DIFFUSIONRHS_MF]->nComp()!=1)
  amrex::Error("localMF[DIFFUSIONRHS_MF]->nComp() invalid");

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 for (int im=0;im<num_materials;im++) {
  if ((compressible_dt_factor[im]>=1.0)&&
      (compressible_dt_factor[im]<=1.0e+20)) {
   // do nothing
  } else
   amrex::Error("compressible_dt_factor[im] invalid");
 }

 MultiFab* denmf=getStateDen(1,cur_time_slab);  // num_materials x den,temp, ...
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
  
 if (project_option==SOLVETYPE_PRES) {
  if (state_index!=State_Type)
   amrex::Error("state_index invalid");
 } else
  amrex::Error("project_option invalid28 init_advective_pressure");

 MultiFab& S_new=get_new_data(state_index,slab_step+1);
 MultiFab& LS_new=get_new_data(LS_Type,slab_step+1);

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
 
  if (lsnewfab.nComp()!=num_materials*(AMREX_SPACEDIM+1))
   amrex::Error("lsnewfab.nComp()!=num_materials*(AMREX_SPACEDIM+1)");
 
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

    // declared in: NAVIERSTOKES_3D.F90
  fort_advective_pressure(
   &level,
   &finest_level,
   xlo,dx,
   &dt_slab, //fort_advective_pressure
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
   &nden,
   compressible_dt_factor.dataPtr(),
   &pressure_select_criterion, 
   &project_option);

 } // mfi
} // omp
 ns_reconcile_d_num(LOOP_ADVECTIVE_PRESSURE,"init_advective_pressure");

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

 int ncomp=num_materials*num_state_material;
 int scomp=STATECOMP_STATES;
 MultiFab* mf=getState(ngrow,scomp,ncomp,time);
 return mf;

}  // subroutine getStateDen


MultiFab* NavierStokes::getStatePres(int ngrow,Real time) {

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 MultiFab* mf=getState(ngrow,STATECOMP_PRES,1,time);
 return mf;

}  // subroutine getStatePres

void NavierStokes::scale_variablesALL() {

 if (level!=0)
  amrex::Error("level invalid scale_variablesALL");

 int finest_level = parent->finestLevel();
  // in: PROB.F90
 fort_setfortscales(&projection_pressure_scale,
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
 fort_setfortscales(&dummy_scale,&dummy_scale);

 dt_slab/=projection_velocity_scale;

 int scale_flag=1;
 for (int ilev=finest_level;ilev>=level;ilev--) {
  NavierStokes& ns_level=getLevel(ilev);
  ns_level.scale_variables(scale_flag);
 }

} // end subroutine NavierStokes::unscale_variablesALL()

//scale Unew,Umac_new,P,mdot,solid_vars,even components of CELL_SOUND (padvect)
void NavierStokes::scale_variables(int scale_flag) {

 std::string local_caller_string="scale_variables";
 int finest_level = parent->finestLevel();
 if ((level<0)||(level>finest_level))
  amrex::Error("level invalid scale_variables");

 int nparts=im_solid_map.size();
 if ((nparts<0)||(nparts>num_materials))
  amrex::Error("nparts invalid");

 int nparts_def=nparts;
 if (nparts==0) {
  nparts_def=1;
 } else if ((nparts>=1)&&(nparts<=num_materials)) {
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
  debug_ngrow(FSI_GHOST_MAC_MF+data_dir,0,local_caller_string);
 }

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  debug_ngrow(FACE_VAR_MF+dir,0,local_caller_string);
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

 S_new.mult(pres_factor,STATECOMP_PRES,nsolve,0);

  // DIV_new contains -dt (pnew-padv)/(rho c^2 dt^2) + MDOT_MF dt/vol
 DIV_new.mult(vel_factor,0,nsolve,0);

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  MultiFab& Umac_new=get_new_data(Umac_Type+dir,slab_step+1);
  MultiFab& Umac_old=get_new_data(Umac_Type+dir,slab_step);
  int ncmac=Umac_new.nComp();
  if (ncmac!=nsolve) {
   std::cout << "num_materials = " << num_materials << '\n';
   std::cout << "ncmac = " << ncmac << '\n';
   amrex::Error("ncmac invalid scale_variables");
  }
  Umac_new.mult(vel_factor,0,nsolve,0);
  Umac_old.mult(vel_factor,0,nsolve,0);
  localMF[FACE_VAR_MF+dir]->mult(vel_factor,FACECOMP_FACEVEL,1,0);
  localMF[FSI_GHOST_MAC_MF+dir]->mult(vel_factor,0,nparts_def*AMREX_SPACEDIM,0);
 } // dir

 localMF[DIFFUSIONRHS_MF]->mult(pres_factor,0,nsolve,0);

 // coeff_avg,padvect_avg 
 localMF[CELL_SOUND_MF]->mult(pres_factor,1,1,0);

 if ((nparts>=1)&&(nparts<=num_materials)) {
  MultiFab& Solid_new=get_new_data(Solid_State_Type,slab_step+1);
  Solid_new.mult(vel_factor,0,nparts*AMREX_SPACEDIM,0);
 } else if (nparts==0) {
  // do nothing
 } else
  amrex::Error("nparts invalid");

}  // end subroutine scale_variables


// 1. viscosity coefficient - 1..num_materials
// 2. viscoelastic coefficient - 1..num_materials
// 3. relaxation time - 1..num_materials
// the viscous and viscoelastic forces should both be multiplied by
// visc_coef.  
// "getStateVISC" is called by "getStateVISC_ALL"
void NavierStokes::getStateVISC(const std::string& caller_string) {

 std::string local_caller_string="getStateVISC";
 local_caller_string=caller_string+local_caller_string;

 int ngrow=1;

 delete_localMF_if_exist(CELL_VISC_MATERIAL_MF,1);

 int finest_level=parent->finestLevel();

 bool use_tiling=ns_tiling;

  // init_gradu_tensorALL is called in NavierStokes::make_physics_varsALL
  // prior to this routine being called.
  // Also, init_gradu_tensorALL is called in NavierStokes::writeTECPLOT_File
  // prior to this routine being called.
 debug_ngrow(CELLTENSOR_MF,1,local_caller_string);
 if (localMF[CELLTENSOR_MF]->nComp()!=AMREX_SPACEDIM_SQR)
  amrex::Error("localMF[CELLTENSOR_MF]->nComp() invalid");

 VOF_Recon_resize(2); //output:SLOPE_RECON_MF
 debug_ngrow(SLOPE_RECON_MF,2,local_caller_string);
 if (localMF[SLOPE_RECON_MF]->nComp()!=num_materials*ngeom_recon)
  amrex::Error("localMF[SLOPE_RECON_MF]->nComp() invalid");

  // viscosity:                    1,...,nmat
  // viscoelastic:                 nmat+1,...,2*nmat
  // viscoelastic relaxation time: nmat+2,...,3*nmat
 int ncomp_visc=3*num_materials;

 for (int im=0;im<num_materials;im++) {

  if (ns_is_rigid(im)==1) {
   // do nothing
  } else if (ns_is_rigid(im)==0) {

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

 } // im=0..num_materials-1

  // viscosity:                    1,...,nmat
  // viscoelastic:                 nmat+1,...,2*nmat
  // viscoelastic relaxation time: nmat+2,...,3*nmat
 new_localMF(CELL_VISC_MATERIAL_MF,ncomp_visc,ngrow,-1);//sets values to 0.0

 MultiFab* vel=getState(ngrow+1,STATECOMP_VEL,STATE_NCOMP_VEL,cur_time_slab);

 MultiFab* EOSdata=getStateDen(ngrow,cur_time_slab);

 MultiFab* tensor=vel;

 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=num_materials)) {

  tensor=getStateTensor(ngrow,0,
    num_materials_viscoelastic*ENUM_NUM_TENSOR_TYPE_REFINE,cur_time_slab);

 } else if (num_materials_viscoelastic==0) {
	 // do nothing
 } else
  amrex::Error("num_materials_viscoelastic invalid");

 for (int im=0;im<num_materials;im++) {

  const Real* dx = geom.CellSize();
  NavierStokes& ns_level0=getLevel(0);
  const Real* dx_coarsest = ns_level0.geom.CellSize();

  MultiFab* gammadot_mf=new MultiFab(grids,dmap,1,ngrow,
	MFInfo().SetTag("gammadot_mf"),FArrayBoxFactory());
  gammadot_mf->setVal(0.0);

  int scomp_tensor=0;

  if (store_elastic_data[im]==1) {

   int partid=0;
   while ((im_viscoelastic_map[partid]!=im)&&
          (partid<im_viscoelastic_map.size())) {
    partid++;
   }
   if ((partid>=0)&&(partid<im_viscoelastic_map.size())) {
    scomp_tensor=partid*ENUM_NUM_TENSOR_TYPE_REFINE;
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

    // since only_scalar==1, this routine calculates: sqrt(2 D:D)
    // since D is symmetric, D:D=trace(D^2) 
    // is invariant to coordinate transformations.
    // if levelrz==COORDSYS_RZ, gradu(3,3)=u/|r|
    int only_scalar=1;  // sqrt(2 D:D)
    int destcomp=0;
    level_getshear(gammadot_mf,vel,only_scalar,destcomp,ngrow);

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

   FArrayBox& viscfab=(*localMF[CELL_VISC_MATERIAL_MF])[mfi];

   FArrayBox& velfab=(*vel)[mfi];
   FArrayBox& eosfab=(*EOSdata)[mfi];
   FArrayBox& tensorfab=(*tensor)[mfi];

   Vector<int> velbc=getBCArray(State_Type,gridno,
     STATECOMP_VEL,STATE_NCOMP_VEL);

   int fortran_im=im+1;

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    // fort_derviscosity is declared in: DERIVE_3D.F90
   fort_derviscosity(
      &level,
      &finest_level,
      &visc_coef,
      &fortran_im,
      &dt_slab, //used for viscoelastic coefficient
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
  ns_reconcile_d_num(LOOP_DERVISCOSITY,"getStateVISC");

  if ((les_model[im]==1)||
      (viscconst_eddy_wall[im]>0.0)||
      (viscconst_eddy_bulk[im]>0.0)) {

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

    // viscosity:                    1,...,nmat
    // viscoelastic:                 nmat+1,...,2*nmat
    // viscoelastic relaxation time: nmat+2,...,3*nmat
    FArrayBox& viscfab=(*localMF[CELL_VISC_MATERIAL_MF])[mfi];

    FArrayBox& cellten=(*localMF[CELLTENSOR_MF])[mfi];
    if (cellten.nComp()!=AMREX_SPACEDIM_SQR)
     amrex::Error("cellten invalid ncomp");

    FArrayBox& velfab=(*vel)[mfi];
    FArrayBox& eosfab=(*EOSdata)[mfi];

    Vector<int> velbc=getBCArray(State_Type,gridno,
      STATECOMP_VEL,STATE_NCOMP_VEL);

    int fortran_im=im+1;

    int tid_current=ns_thread();
    if ((tid_current<0)||(tid_current>=thread_class::nthreads))
     amrex::Error("tid_current invalid");
    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

      // declared in: DERIVE_3D.F90
      // WALE model, "viscconst_eddy_wall", "viscconst_eddy_bulk",
      // effects are added to "viscfab" (CELL_VISC_MATERIAL_MF)
    fort_derturbvisc(
      &les_model[im],
      &level,
      &fortran_im,
      &dt_slab, //derturbvisc
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
      dx,
      xlo,
      dx_coarsest,
      &ngrow,
      &ncomp_visc);
   } //mfi
} // omp
   ns_reconcile_d_num(LOOP_DERTURBVISC,"getStateVISC");

  } else if ((les_model[im]==0)&&
             (viscconst_eddy_wall[im]==0.0)&&
	     (viscconst_eddy_bulk[im]==0.0)) {
   // do nothing
  } else
   amrex::Error("les_model,viscconst_eddy_wall, or bulk invalid");

  delete gammadot_mf;

  Plus_localMF(CELL_VISC_MATERIAL_MF,
    viscconst_artificial[im],im,1,ngrow);
 } // im=0..num_materials-1

 delete vel;
 delete EOSdata;

 if ((num_materials_viscoelastic>=1)&&
     (num_materials_viscoelastic<=num_materials)) {
  delete tensor;
 } else if (num_materials_viscoelastic==0) {
  // do nothing
 } else
  amrex::Error("num_materials_viscoelastic invalid");

}  // end subroutine getStateVISC


//CELL_CONDUCTIVITY_MATERIAL_MF is deleted in ::Geometry_cleanup()
void NavierStokes::getStateCONDUCTIVITY() {

 std::string local_caller_string="getStateCONDUCTIVITY";

 int ngrow=1;

 delete_localMF_if_exist(CELL_CONDUCTIVITY_MATERIAL_MF,1);

 int finest_level=parent->finestLevel();

 bool use_tiling=ns_tiling;

 VOF_Recon_resize(ngrow+2); //output:SLOPE_RECON_MF
 debug_ngrow(SLOPE_RECON_MF,ngrow+2,local_caller_string);
 if (localMF[SLOPE_RECON_MF]->nComp()!=num_materials*ngeom_recon)
  amrex::Error("localMF[SLOPE_RECON_MF]->nComp() invalid");

 //sets values to 0.0
 new_localMF(CELL_CONDUCTIVITY_MATERIAL_MF,num_materials,ngrow,-1);

 MultiFab* EOSdata=getStateDen(ngrow+2,cur_time_slab);
 const Real* dx = geom.CellSize();

 for (int im=0;im<num_materials;im++) {

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

   FArrayBox& conductivity_fab=(*localMF[CELL_CONDUCTIVITY_MATERIAL_MF])[mfi];

   FArrayBox& eosfab=(*EOSdata)[mfi];

   FArrayBox& voffab=(*localMF[SLOPE_RECON_MF])[mfi];

   int fortran_im=im+1;

   int tid_current=ns_thread();
   if ((tid_current<0)||(tid_current>=thread_class::nthreads))
    amrex::Error("tid_current invalid");
   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

   int local_sumdata_size=NS_sumdata.size();

    // declared in: DERIVE_3D.F90
   fort_derconductivity(
     &local_sumdata_size,
     NS_sumdata.dataPtr(),
     &ncomp_sum_int_user1,
     &ncomp_sum_int_user2,
     &level,
     &finest_level,
     &fortran_im,
     &dt_slab, //fort_derconductivity
     conductivity_fab.dataPtr(),
     ARLIM(conductivity_fab.loVect()),
     ARLIM(conductivity_fab.hiVect()),
     eosfab.dataPtr(),
     ARLIM(eosfab.loVect()),ARLIM(eosfab.hiVect()),
     voffab.dataPtr(),
     ARLIM(voffab.loVect()),ARLIM(voffab.hiVect()),
     tilelo,tilehi,
     fablo,fabhi,
     &bfact,
     &cur_time_slab,
     dx,xlo,
     &ngrow);
  } //mfi
} // omp
  ns_reconcile_d_num(LOOP_DERCONDUCTIVITY,"getStateCONDUCTIVITY");

 } // im=0..num_materials-1

 delete EOSdata;

}  // end subroutine getStateCONDUCTIVITY



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
 
 std::string local_caller_string="getState_tracemag";

 bool use_tiling=ns_tiling;

 if (num_state_base!=2)
  amrex::Error("num_state_base invalid");

 int finest_level=parent->finestLevel();

 if (localMF_grow[idx]>=0)
  amrex::Error("local magtrace not previously deleted");

 int ntrace=5*num_materials;
  //ngrow=1
 new_localMF(idx,ntrace,1,-1);

 debug_ngrow(CELLTENSOR_MF,1,local_caller_string);
 if (localMF[CELLTENSOR_MF]->nComp()!=AMREX_SPACEDIM_SQR)
  amrex::Error("localMF[CELLTENSOR_MF]->nComp() invalid");

  //ngrow=1
 MultiFab* den_data=getStateDen(1,cur_time_slab);
  //ngrow=2
 MultiFab* vel_data=getState(2,0,AMREX_SPACEDIM,cur_time_slab);

 debug_ngrow(CELL_VISC_MATERIAL_MF,1,local_caller_string);
 int ncomp_visc=localMF[CELL_VISC_MATERIAL_MF]->nComp();
 if (ncomp_visc!=3*num_materials)
  amrex::Error("visc_data invalid ncomp");

 int ncomp_den=den_data->nComp();

 const Real* dx = geom.CellSize();

 for (int im=0;im<num_materials;im++) {

  MultiFab* tensor=den_data;
  int allocate_tensor=0;
  if (ns_is_rigid(im)==0) {
   if (elastic_viscosity[im]>0.0) {
    int partid=0;
    while ((im_viscoelastic_map[partid]!=im)&&
	   (partid<im_viscoelastic_map.size())) {
     partid++;
    }
    if ((partid>=0)&&(partid<im_viscoelastic_map.size())) {
     int scomp_tensor=partid*ENUM_NUM_TENSOR_TYPE_REFINE;
      //ngrow=1
     tensor=getStateTensor(1,scomp_tensor,
        ENUM_NUM_TENSOR_TYPE_REFINE,cur_time_slab);
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

  int idest=5*im;
  int only_scalar=1;  // sqrt(2 D:D)
  int ngrow_getshear=1;
  level_getshear(localMF[idx],vel_data,only_scalar,
     idest,ngrow_getshear);

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
   if (cellten.nComp()!=AMREX_SPACEDIM_SQR)
    amrex::Error("cellten invalid ncomp");

   FArrayBox& destfab=(*localMF[idx])[mfi];

   FArrayBox& viscfab=(*localMF[CELL_VISC_MATERIAL_MF])[mfi];

   FArrayBox& velfab=(*vel_data)[mfi];
   FArrayBox& denfab=(*den_data)[mfi];
   FArrayBox& tenfab=(*tensor)[mfi];

   Vector<int> bc=getBCArray(State_Type,gridno,
     STATECOMP_VEL,STATE_NCOMP_VEL);

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

    //fort_dermagtrace is declared in: DERIVE_3D.F90
   fort_dermagtrace(
    &level,
    &finest_level,  
    &im, //im=0..num_materials-1
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
    &ncomp_visc,
    &ntrace, 
    polymer_factor.dataPtr(),
    etaS.dataPtr(),
    etaP.dataPtr(),
    Carreau_beta.dataPtr(),
    elastic_time.dataPtr(),
    viscoelastic_model.dataPtr(),
    elastic_viscosity.dataPtr());
  } //mfi
} // omp
  ns_reconcile_d_num(LOOP_DERMAGTRACE,"getState_tracemag");

  if (allocate_tensor==0) {
   // do nothing
  } else if (allocate_tensor==1) {
   delete tensor;
  } else
   amrex::Error("allocate_tensor invalid");

 } // im=0..num_materials-1

 delete den_data;
 delete vel_data;

}  // getState_tracemag

}/* namespace amrex */
