
#include <local_thread_class.H>
#include <AMReX_Utility.H>
#include <FillPatchUtil.H>
#include <cmath>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <EXTRAP_COMP.H>

namespace amrex
{

void FillPatchSingleLevel (
  int level,
  MultiFab& mf_target, 
  Real time, 
  const MultiFab& smf, 
  int scomp, //absolute within the state smf.
  int dcomp, 
  int ncomp,
  const Geometry& geom, 
  PhysBCFunctBaseSUSSMAN& physbcf,
  Vector<int> scompBC_map,
  int bfact,
  int debug_fillpatch) {

 BL_PROFILE("FillPatchSingleLevel");

 if (scomp+ncomp <= smf.nComp()) {
  // do nothing
 } else {
  std::cout << "scomp or ncomp invalid scomp=" <<
   scomp << " ncomp= " << ncomp << "smf.nComp()= " <<
   smf.nComp() << '\n';
  amrex::Error("scomp+ncomp <= smf.nComp() failed");
 }

 if (dcomp+ncomp <= mf_target.nComp()) {
  // do nothing
 } else {
  std::cout << "dcomp or ncomp invalid dcomp=" <<
   dcomp << " ncomp= " << ncomp << "mf_target.nComp()= " <<
   mf_target.nComp() << '\n';
  amrex::Error("dcomp+ncomp <= mf_target.nComp() failed");
 }

 if (scompBC_map.size()!=ncomp)
  amrex::Error("scompBC_map has invalid size");

 if (debug_fillpatch==1) {

  for (int pid=0;pid<amrex::ParallelDescriptor::NProcs();pid++) {
   amrex::ParallelDescriptor::Barrier();
   if (amrex::ParallelDescriptor::MyProc()==pid) {
    std::fflush(NULL);
    std::cout << "FillPatchSingleLevel(1) PROC= " << 
     amrex::ParallelDescriptor::MyProc() << '\n';
    std::cout << "scomp=" << scomp << " dcomp= " << dcomp << " ncomp= " <<
     ncomp << '\n';
    std::cout << "smf.boxArray() " << smf.boxArray() << '\n';
    std::cout << "smf.DistributionMap() " << smf.DistributionMap() << '\n';
    std::cout << "mf_target.boxArray() " << mf_target.boxArray() << '\n';
    std::cout << "mf_target.DistributionMap() " << 
          mf_target.DistributionMap() << '\n';
    Periodicity my_period=geom.periodicity();
    for (int dir=0;dir<AMREX_SPACEDIM;dir++)
     std::cout << "dir=" << dir << " my_period.isPeriodic(dir) " << 
       my_period.isPeriodic(dir) << '\n';
    std::fflush(NULL);
   }
  }  // pid=0..NProcs-1

  std::fflush(NULL);
  for (int n=0;n<ncomp;n++) {
   Real smf_norm=smf.norm1(scomp+n,0);
   Real mf_norm=mf_target.norm1(dcomp+n,0);
   if (ParallelDescriptor::IOProcessor()) {
    std::cout << "n= " << n << " BCmap(n)= " << scompBC_map[n] <<
     " smf_norm= " << smf_norm << " mf_norm= " << mf_norm << '\n';
   }
  }
  std::fflush(NULL);

 } else if (debug_fillpatch==0) {
  // do nothing
 } else {
  amrex::Error("debug_fillpatch invalid");
 }

  // Ranks which do not own any destination data will not wait for the
  // other ranks in the ParallelCopy command. 
  // Barrier statements are inserted here just in case the
  // asynchronous feature of the ParallelCopy command might lead
  // to problems in future operations.
 ParallelDescriptor::Barrier();

  // src,src_comp,dest_comp,num_comp,src_nghost,dst_nghost,period
 mf_target.ParallelCopy(smf, scomp, dcomp, ncomp, IntVect{0}, 
    mf_target.nGrowVect(), geom.periodicity());

 ParallelDescriptor::Barrier();

 if (debug_fillpatch==1) {
  for (int pid=0;pid<amrex::ParallelDescriptor::NProcs();pid++) {
   amrex::ParallelDescriptor::Barrier();
   if (amrex::ParallelDescriptor::MyProc()==pid) {
    std::fflush(NULL);
    std::cout << "FillPatchSingleLevel(2) PROC= " << 
     amrex::ParallelDescriptor::MyProc() << '\n';
    std::cout << "scomp=" << scomp << " dcomp= " << dcomp << " ncomp= " <<
     ncomp << '\n';
    std::cout << "smf.boxArray() " << smf.boxArray() << '\n';
    std::cout << "smf.DistributionMap() " << smf.DistributionMap() << '\n';
    std::cout << "mf_target.boxArray() " << mf_target.boxArray() << '\n';
    std::cout << "mf_target.DistributionMap() " << 
      mf_target.DistributionMap() << '\n';
    Periodicity my_period=geom.periodicity();
    for (int dir=0;dir<AMREX_SPACEDIM;dir++)
     std::cout << "dir=" << dir << " my_period.isPeriodic(dir) " << 
       my_period.isPeriodic(dir) << '\n';
    std::fflush(NULL);
   }
  } // pid=0..NProcs-1
  std::fflush(NULL);
  for (int n=0;n<ncomp;n++) {
   Real smf_norm=smf.norm1(scomp+n,0);
   Real mf_norm=mf_target.norm1(dcomp+n,0);
   if (ParallelDescriptor::IOProcessor()) {
    std::cout << "n= " << n << " BCmap(n)= " << scompBC_map[n] <<
     " smf_norm= " << smf_norm << " mf_norm= " << mf_norm << '\n';
   }
  }
  std::fflush(NULL);
 } else if (debug_fillpatch==0) {
  // do nothing
 } else {
  amrex::Error("debug_fillpatch invalid");
 }

  //StateDataPhysBCFunct::FillBoundary is declared in: StateData.cpp
 physbcf.FillBoundary(level,mf_target,time,dcomp,scompBC_map,ncomp,bfact);
}

//FillPatchTwoLevels is called from AmrLevel.cpp:
// AmrLevel::FillPatch
// AmrLevel::InterpBordersGHOST
// AmrLevel::InterpBorders
void FillPatchTwoLevels (
 MultiFab& mf_target,   // target
 Real time,
 MultiFab& cmf,  // coarse
 MultiFab& fmf,  // fine
 int scomp, 
 int dcomp, 
 int ncomp,
 const Geometry& cgeom, 
 const Geometry& fgeom, 
 PhysBCFunctBaseSUSSMAN& cbc, 
 PhysBCFunctBaseSUSSMAN& fbc,
 Interpolater* mapper, 
 const Vector<BCRec>& global_bcs,
 Vector<int> scompBC_map,
 int levelc,int levelf,
 int bfactc,int bfactf,
 int grid_type,
 int debug_fillpatch) {

 BL_PROFILE("FillPatchTwoLevels");

 if ((levelc<0)||(levelc!=levelf-1))
  amrex::Error("levelc or levelf invalid");
 if (scompBC_map.size()!=ncomp)
  amrex::Error("scompBC_map has invalid size");
 if (global_bcs.size()<ncomp) {
  std::cout << "time= " << time << '\n';
  std::cout << "scomp,dcomp,ncomp= " << scomp << ' ' << dcomp <<
   ' ' << ncomp << '\n';
  std::cout << "levelc, levelf, bfactc, bfactf = " << levelc << ' ' << 
   levelf << ' ' << bfactc << ' ' << bfactf << '\n';
  std::cout << "global_bcs.size() " << global_bcs.size() << '\n';
  amrex::Error("global_bcs has invalid size");
 }

#ifdef AMREX_USE_EB
 amrex::Error("AMREX_USE_EB should never be defined");
#else
 EB2::IndexSpace const* index_space=nullptr;   
#endif

 int ngrow = mf_target.nGrow();
 const IntVect& ngrow_vec=mf_target.nGrowVect();	   
 
 int do_the_interp=0;

 if (ngrow>0) {
  do_the_interp=1;
 } else if (ngrow==0) {
  // do nothing
 } else {
  amrex::Error("ngrow invalid");
 }
 if (do_the_interp==0) {
  if (mf_target.boxArray()!=fmf.boxArray()) {
   do_the_interp=1;
  }
  if (do_the_interp==0) {
   if (mf_target.DistributionMap()!=fmf.DistributionMap()) {
    do_the_interp=1;
   }
  } else if (do_the_interp==1) {
   // do nothing
  } else
   amrex::Error("do_the_interp invalid");
 } else if (do_the_interp==1) {
  // do nothing
 } else
  amrex::Error("do_the_interp invalid");

 if (debug_fillpatch==1) {
  for (int pid=0;pid<amrex::ParallelDescriptor::NProcs();pid++) {
   amrex::ParallelDescriptor::Barrier();
   if (amrex::ParallelDescriptor::MyProc()==pid) {
    std::fflush(NULL);
    std::cout << "FillPatchTwoLevels(1) PROC= " << 
     amrex::ParallelDescriptor::MyProc() << '\n';
    std::cout << "scomp=" << scomp << " dcomp= " << dcomp << 
     " ncomp=" << ncomp << " levelc= " << levelc <<
     " levelf= " << levelf << " do_the_interp= " << do_the_interp << '\n';
    std::cout << "mf_target.boxArray() " << mf_target.boxArray() << '\n';
    std::cout << "mf_target.DistributionMap() " << 
      mf_target.DistributionMap() << '\n';
    std::cout << "fmf.boxArray() " << fmf.boxArray() << '\n';
    std::cout << "fmf.DistributionMap() " << fmf.DistributionMap() << '\n';
    std::cout << "cmf.boxArray() " << cmf.boxArray() << '\n';
    std::cout << "cmf.DistributionMap() " << cmf.DistributionMap() << '\n';
    std::fflush(NULL);
   }
  }  // pid=0..NProcs-1
 } else if (debug_fillpatch==0) {
  // do nothing
 } else {
  amrex::Error("debug_fillpatch invalid");
 }

 if (do_the_interp==1) {

//*mapper is derived from class Interpolater.
//Inside of Interpolater:
// virtual InterpolaterBoxCoarsener BoxCoarsener (int bfactc,int bfactf,
//   int grid_type);
  const InterpolaterBoxCoarsener& coarsener = 
    mapper->BoxCoarsener(bfactc,bfactf,grid_type);
        
  Box fdomain = fgeom.Domain();
  fdomain.convert(mf_target.boxArray().ixType());

    // find coarsen( mf_target intersect complement(fmf) within fdomain_g ).
    // note: fdomain_g is local to TheFPinfo and 
    // fdomain_g=fdomain.grow(ngrow) in the periodic directions.
    // the valid region of fpc encompass the mf_target grow regions, but the
    // boxes of fpc are disjoint.
    // "TheFPinfo" is declared in "AMReX_FabArrayBase.H"
  const FabArrayBase::FPinfo& fpc = 
   FabArrayBase::TheFPinfo(fmf, mf_target, 
      ngrow_vec, 
      coarsener,//type InterpolatorBoxCoarsener, derived from BoxConverter
      fgeom,
      cgeom,
      index_space);

  bool empty_flag=fpc.ba_crse_patch.empty();

  if (empty_flag==true) {
   // do nothing, no data on the coarse level is needed to
   // fill the fine level target.
  } else if (empty_flag==false) {

     //ngrow=0
   MultiFab mf_crse_patch(
     fpc.ba_crse_patch,
     fpc.dm_patch,
     ncomp,
     0,
     MFInfo().SetTag("mf_crse_patch"),
     *fpc.fact_crse_patch);

   mf_crse_patch.setDomainBndry(std::numeric_limits<Real>::quiet_NaN(),cgeom);
   
    // This data will be interpolated next.       
   FillPatchSingleLevel(
    levelc,
    mf_crse_patch, // target data (coarse)
    time, 
    cmf,  // coarse data source
    scomp, 
    0,   // dst_comp
    ncomp, 
    cgeom, 
    cbc,
    scompBC_map,
    bfactc,
    debug_fillpatch);

   MultiFab mf_fine_patch(fpc.ba_fine_patch,fpc.dm_patch,ncomp,0,
		  MFInfo().SetTag("mf_fine_patch"),
		  *fpc.fact_fine_patch);

   Vector< BCRec > local_bcs;
   local_bcs.resize(ncomp);
   for (int i=0;i<ncomp;i++)
    local_bcs[i]=global_bcs[scompBC_map[i]]; 
    
   bool cellcen = fpc.ba_crse_patch.ixType().cellCentered();
   if ((cellcen!=true)&&(cellcen!=false))
    amrex::Error("cellcen bust");

   if (thread_class::nthreads<1)
    amrex::Error("thread_class::nthreads invalid");
   thread_class::init_d_numPts(mf_fine_patch.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel if (cellcen)
#endif
{
   for (MFIter mfi(mf_fine_patch,false); mfi.isValid(); ++mfi) {
    const Box& tilegrid=mfi.tilebox();

    int tid_current=0;
#ifdef _OPENMP
    tid_current = omp_get_thread_num();
#endif
    if ((tid_current>=0)&&(tid_current<thread_class::nthreads)) {
     // do nothing
    } else
     amrex::Error("tid_current invalid");

    thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

    FArrayBox& sfab=mf_crse_patch[mfi];
    FArrayBox& dfab=mf_fine_patch[mfi];
    const Box& dbx = dfab.box();

    Vector<BCRec> bcr(ncomp);
    int src_comp_bcs=0;
    int dest_comp_bcr=0;
    amrex::setBC(dbx,fdomain,src_comp_bcs,dest_comp_bcr,ncomp,
     local_bcs,bcr);
  	   
     // parameter "dfab" used to correspond to mf_target.
    mapper->interp(time,
               sfab, // source
               0,
  	       dfab, //dest; dest BoxArray is disjoint from the existing fine.
  	       0,    //was dcomp
  	       ncomp,
  	       dbx,
  	       cgeom,
  	       fgeom,
  	       bcr,
               levelc,levelf,
               bfactc,bfactf,
               grid_type);
   } // mfi
} // omp
   thread_class::sync_tile_d_numPts();
   ParallelDescriptor::ReduceRealSum(thread_class::tile_d_numPts[0]);
   thread_class::reconcile_d_numPts(LOOP_MAPPER_INTERP,"FillPatchTwoLevels");

   ParallelDescriptor::Barrier();
    // src,src_comp,dest_comp,num_comp,src_nghost,dst_nghost,period
   mf_target.ParallelCopy(mf_fine_patch, 0, dcomp, ncomp, IntVect{0}, 
      ngrow_vec,fgeom.periodicity());
   ParallelDescriptor::Barrier();

  } else {
   amrex::Error("empty_flag invalid");
  }
 } else if (do_the_interp==0) {
  // do nothing, mf and fmf have same boxarray's and dmaps and
  // ngrow==0
 } else
  amrex::Error("do_the_interp invalid");

 if (debug_fillpatch==1) {
  for (int pid=0;pid<amrex::ParallelDescriptor::NProcs();pid++) {
   amrex::ParallelDescriptor::Barrier();
   if (amrex::ParallelDescriptor::MyProc()==pid) {
    std::fflush(NULL);
    std::cout << "FillPatchTwoLevels(2) PROC= " << 
     amrex::ParallelDescriptor::MyProc() << '\n';
    std::cout << "scomp=" << scomp << " dcomp= " << dcomp << 
     " ncomp=" << ncomp << " levelc= " << levelc <<
     " levelf= " << levelf << '\n';
    std::fflush(NULL);
   }
  }  // pid=0..NProcs()-1
 } else if (debug_fillpatch==0) {
  // do nothing
 } else {
  amrex::Error("debug_fillpatch invalid");
 }

 FillPatchSingleLevel(
   levelf,
   mf_target,  //mf_target init with coarse data in regions not covered by fmf
   time, 
   fmf, 
   scomp, 
   dcomp, 
   ncomp, 
   fgeom, 
   fbc,
   scompBC_map,
   bfactf,
   debug_fillpatch);

 if (debug_fillpatch==1) {
  for (int pid=0;pid<amrex::ParallelDescriptor::NProcs();pid++) {
   amrex::ParallelDescriptor::Barrier();
   if (amrex::ParallelDescriptor::MyProc()==pid) {
    std::fflush(NULL);
    std::cout << "FillPatchTwoLevels(3) PROC= " << 
     amrex::ParallelDescriptor::MyProc() << '\n';
    std::cout << "scomp=" << scomp << " dcomp= " << dcomp << 
     " ncomp=" << ncomp << " levelc= " << levelc <<
     " levelf= " << levelf << '\n';
    std::fflush(NULL);
   }
  } // pid=0..NProcs()-1
 } else if (debug_fillpatch==0) {
  // do nothing
 } else {
  amrex::Error("debug_fillpatch invalid");
 }
}  //FillPatchTwoLevels


//SUSSMAN
void FillPatchTower (
 int ngrow_root,
 int called_from_regrid,
 int finest_top_level,
 int top_level,
 MultiFab& mf_target, //target
 Real time,
 Vector<MultiFab*> tower_data,
 int scomp, 
 int dcomp, 
 int ncomp,
 Vector<const Geometry*> tower_geom,
 Vector<PhysBCFunctBaseSUSSMAN*> tower_physbc,
 Interpolater* mapper, 
 const Vector<BCRec>& global_bcs,
 Vector<int> scompBC_map,
 Vector<int> tower_bfact,
 int grid_type,
 int debug_fillpatch) {

 BL_PROFILE("FillPatchTower");

 if ((top_level>=0)&&(top_level<=finest_top_level)) {
  //do nothing
 } else
  amrex::Error("top_level invalid");

 if ((tower_data.size()==finest_top_level+1)&&
     (tower_geom.size()==finest_top_level+1)&&
     (tower_physbc.size()==finest_top_level+1)&&
     (tower_bfact.size()==finest_top_level+1)) {
  //do nothing
 } else
  amrex::Error("tower size() invalid");

 if (scompBC_map.size()!=ncomp)
  amrex::Error("scompBC_map has invalid size");
 if (global_bcs.size()<ncomp) {
  std::cout << "time= " << time << '\n';
  std::cout << "scomp,dcomp,ncomp= " << scomp << ' ' << dcomp <<
   ' ' << ncomp << '\n';
  std::cout << "top_level=" << top_level <<'\n';

  std::cout << "global_bcs.size() " << global_bcs.size() << '\n';
  amrex::Error("global_bcs has invalid size");
 }

#ifdef AMREX_USE_EB
 amrex::Error("AMREX_USE_EB should never be defined");
#else
 EB2::IndexSpace const* index_space=nullptr;   
#endif

 if (top_level==0) {

  FillPatchSingleLevel(
   top_level,
   mf_target,
   time, 
   *tower_data[top_level],
   scomp, 
   dcomp, 
   ncomp,
   *tower_geom[top_level],
   *tower_physbc[top_level],
   scompBC_map,
   tower_bfact[top_level],
   debug_fillpatch);

 } else if (top_level>=1) {

  IntVect ngrow_vec(AMREX_SPACEDIM);
  for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
   ngrow_vec[dir]=ngrow_root;
  }
 
  int do_the_interp=0;

  if (ngrow_root>=0) {
   //do nothing
  } else
   amrex::Error("ngrow_root invalid");

  if (ngrow_root>0) {
   do_the_interp=1;
  } else if (ngrow_root==0) {
   // do nothing
  } else {
   amrex::Error("ngrow_root invalid");
  }

  if (finest_top_level>top_level) {
   do_the_interp=1;
  } else if (finest_top_level==top_level) {
   //do nothing
  } else {
   amrex::Error("finest_top_level invalid");
  }
 
  if (do_the_interp==0) {
   if (mf_target.boxArray()!=
       tower_data[top_level]->boxArray()) {
    do_the_interp=1;
   }
   if (do_the_interp==0) {
    if (mf_target.DistributionMap()!=
        tower_data[top_level]->DistributionMap()) {
     do_the_interp=1;
    }
   } else if (do_the_interp==1) {
    // do nothing
   } else
    amrex::Error("do_the_interp invalid");
  } else if (do_the_interp==1) {
   // do nothing
  } else
   amrex::Error("do_the_interp invalid");

  if (do_the_interp==1) {

//*mapper is derived from class Interpolater.
//Inside of Interpolater:
//virtual InterpolaterBoxCoarsener BoxCoarsener (int bfactc,int bfactf,
//  int grid_type);
//
//Box
//InterpolaterBoxCoarsener::doit (const Box& fine) const
//{
//    return mapper->CoarseBox(fine,bfactc,bfactf,grid_type);
//}
//
   const InterpolaterBoxCoarsener& coarsener = 
     mapper->BoxCoarsener(
       tower_bfact[top_level-1],
       tower_bfact[top_level],
       grid_type);
        
   Box fdomain = tower_geom[top_level]->Domain();
   fdomain.convert(mf_target.boxArray().ixType());

    // find coarsen( mf_target intersect complement(fmf.ba) within fdomain_g ).
    // note: fdomain_g is local to TheFPinfo and 
    // fdomain_g=fdomain.grow(ngrow) in the periodic directions.
    // the valid region of fpc encompass the mf_target grow regions, but the
    // boxes of fpc are disjoint.
    // "TheFPinfo" is declared in "AMReX_FabArrayBase.H" and
    // "AMReX_FabArrayBase.cpp"
    // in: FabArrayBase::FPinfo::FPinfo,
    /*
    for (int i = iboxlo; i <= iboxhi; ++i) {
        Box bx = dstba_simplified[i];
        bx.grow(m_dstng);
        bx &= m_dstdomain;    
        BoxList const& leftover = srcba_simplified.complementIn(bx);
        if (leftover.isNotEmpty()) {
            bl.join(leftover);
        }
    }
    */
 
   const FabArrayBase::FPinfo& fpc = 
    FabArrayBase::TheFPinfo(
      *tower_data[top_level], //source
      mf_target, //dest
      ngrow_vec, //ngrow_vec=IntVect filled with ngrow_root
      coarsener,//type InterpolatorBoxCoarsener, derived from BoxConverter
      *tower_geom[top_level],
      *tower_geom[top_level-1],
      index_space);

   bool empty_flag=fpc.ba_crse_patch.empty();

   if (empty_flag==true) {
    // do nothing, no data on the coarse level is needed to
    // fill the fine level target.
   } else if (empty_flag==false) {

    int ba_crse_patch_ok=1;

    for (int gridno=0;gridno<fpc.ba_crse_patch.size();gridno++) {
     const Box& cbox=fpc.ba_crse_patch[gridno];
     for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      if (cbox.length(i)==0) {
       ba_crse_patch_ok=0;
       std::cout << "ngrow_root= " << ngrow_root << '\n';
       std::cout << "finest_top_level=" << finest_top_level << '\n';
       std::cout << "top_level=" << top_level << '\n';
       std::cout << "i=" << i << '\n';
       std::cout << fpc.ba_crse_patch << '\n';
       amrex::Error("cbox.length error");
      }
     }
    }
    
    if (ba_crse_patch_ok==1) {
 
     MultiFab mf_crse_patch(
       fpc.ba_crse_patch,
       fpc.dm_patch,
       ncomp,
       ngrow_root,
       MFInfo().SetTag("mf_crse_patch"),
       *fpc.fact_crse_patch);

     mf_crse_patch.setDomainBndry(std::numeric_limits<Real>::quiet_NaN(),
       *tower_geom[top_level-1]);
   
     FillPatchTower(
      ngrow_root,
      called_from_regrid,
      finest_top_level,
      top_level-1,
      mf_crse_patch, // target data (coarse)
      time, 
      tower_data, 
      scomp, 
      0,   // dst_comp
      ncomp, 
      tower_geom, 
      tower_physbc,
      mapper,
      global_bcs,
      scompBC_map,
      tower_bfact,
      grid_type,
      debug_fillpatch);

     MultiFab mf_fine_patch(
       fpc.ba_fine_patch,
       fpc.dm_patch,
       ncomp,
       0,
       MFInfo().SetTag("mf_fine_patch"),
       *fpc.fact_fine_patch);

     Vector< BCRec > local_bcs;
     local_bcs.resize(ncomp);
     for (int i=0;i<ncomp;i++)
      local_bcs[i]=global_bcs[scompBC_map[i]]; 
     
     bool cellcen = fpc.ba_crse_patch.ixType().cellCentered();
     if ((cellcen!=true)&&(cellcen!=false))
      amrex::Error("cellcen bust");

     if (thread_class::nthreads<1)
      amrex::Error("thread_class::nthreads invalid");
     thread_class::init_d_numPts(mf_fine_patch.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel if (cellcen)
#endif
{
     for (MFIter mfi(mf_fine_patch,false); mfi.isValid(); ++mfi) {
      const Box& tilegrid=mfi.tilebox();

      int tid_current=0;
#ifdef _OPENMP
      tid_current = omp_get_thread_num();
#endif
      if ((tid_current>=0)&&(tid_current<thread_class::nthreads)) {
       // do nothing
      } else
       amrex::Error("tid_current invalid");

      thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

      FArrayBox& sfab=mf_crse_patch[mfi];
      FArrayBox& dfab=mf_fine_patch[mfi];
      const Box& dbx = dfab.box();

      Vector<BCRec> bcr(ncomp);
      int src_comp_bcs=0;
      int dest_comp_bcr=0;
      amrex::setBC(dbx,fdomain,src_comp_bcs,dest_comp_bcr,ncomp,
       local_bcs,bcr);
   	   
       // parameter "dfab" used to correspond to mf_target.
      mapper->interp(time,
        sfab, // source
        0,
        dfab, //dest; dest BoxArray is disjoint from the existing fine.
        0,    //was dcomp
        ncomp,
        dbx,
        *tower_geom[top_level-1],
        *tower_geom[top_level],
        bcr,
        top_level-1,
        top_level,
        tower_bfact[top_level-1],
        tower_bfact[top_level],
        grid_type);
     } // mfi
} // omp
     thread_class::sync_tile_d_numPts();
     ParallelDescriptor::ReduceRealSum(thread_class::tile_d_numPts[0]);
     thread_class::reconcile_d_numPts(LOOP_MAPPER_INTERP,"FillPatchTower");

     ParallelDescriptor::Barrier();
     // src,src_comp,dest_comp,num_comp,src_nghost,dst_nghost,period
     mf_target.ParallelCopy(mf_fine_patch, 0, dcomp, ncomp, IntVect{0}, 
       ngrow_vec,tower_geom[top_level]->periodicity());
     ParallelDescriptor::Barrier();

    } else if (ba_crse_patch_ok==0) {

     //do nothing

    } else
     amrex::Error("ba_crse_patch_ok invalid");

   } else {
    amrex::Error("empty_flag invalid");
   }
  } else if (do_the_interp==0) {
   // do nothing, mf and fmf have same boxarray's and dmaps and
   // ngrow_root==0
  } else
   amrex::Error("do_the_interp invalid");

  FillPatchSingleLevel(
   top_level,
   mf_target,  //mf_target init with coarse data in regions not covered by fmf
   time, 
   *tower_data[top_level], 
   scomp, 
   dcomp, 
   ncomp, 
   *tower_geom[top_level], 
   *tower_physbc[top_level],
   scompBC_map,
   tower_bfact[top_level],
   debug_fillpatch);

 } else
  amrex::Error("top_level invalid");

}  //FillPatchTower



}  // namespace amrex
