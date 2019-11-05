
#include <sstream>

#include <unistd.h>
#include <memory>
#include <limits>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_ParmParse.H>
#include <AmrLevel.H>
#include <FillPatchUtil.H>

DescriptorList AmrLevel::desc_lst;
DescriptorList AmrLevel::desc_lstGHOST;

namespace amrex {

void
AmrLevel::manual_tags_placement (TagBoxArray&    tags,
                                 Vector<int>& bf_lev)
{}

AmrLevel::AmrLevel () noexcept
{
   parent = 0;
   level = -1;
}

// constructor
AmrLevel::AmrLevel (Amr&            papa,
                    int             lev,
                    const Geometry& level_geom,
                    const BoxArray& ba,
                    const DistributionMapping& dm,
                    Real            time)
    :
    geom(level_geom),
    grids(ba),
    dmap(dm)
{

    if (grids.size()<1) {
     std::cout << "lev= " << lev << " time= " << time << '\n';
     amrex::Error("AmrLevel: grids.size<1");
    }

    level  = lev;
    parent = &papa;


    state.resize(desc_lst.size());

    Real dt=parent->getDt();
    if (dt<=0.0) {
     std::cout << "dt= " << dt << '\n';
     amrex::Error("dt invalid");
    }

    int level_MAX_NUM_SLAB=parent->get_MAX_NUM_SLAB();
    int level_slab_dt_type=parent->get_slab_dt_type();
    if (level_MAX_NUM_SLAB<33)
     amrex::Error("level_MAX_NUM_SLAB too small");
    if ((level_slab_dt_type!=0)&&(level_slab_dt_type!=1))
     amrex::Error("level_slab_dt_type invalid");

    for (int i = 0; i < state.size(); i++)
    {
        int time_order=parent->Time_blockingFactor();
        state[i].define(level,
                        geom.Domain(),
                        grids,
                        dm,
                        desc_lst[i],
                        desc_lstGHOST[i],
                        time,
                        dt,
                        time_order,
                        level_slab_dt_type,
                        level_MAX_NUM_SLAB);
    }

    finishConstructor();
}

void
AmrLevel::restart (Amr&          papa,
                   std::istream& is) {

    parent = &papa;

    is >> level;
    is >> geom;

    grids.readFrom(is);

    int nstate;
    is >> nstate;

    int ndesc = desc_lst.size();

    int ndesc_checkpoint=0;
    for (int j=0;j<ndesc;j++) {
     if (desc_lst[j].store_in_checkpoint()==true)
      ndesc_checkpoint++;
    } 

    if (nstate!=ndesc_checkpoint)
     amrex::Error("nstate and ndesc_checkpoint do not match");

    dmap.define(grids);

    parent->SetBoxArray(level, grids);
    parent->SetDistributionMap(level, dmap);
    int level_MAX_NUM_SLAB=parent->get_MAX_NUM_SLAB();
    int level_slab_dt_type=parent->get_slab_dt_type();

    state.resize(ndesc);
    for (int i = 0; i < ndesc; i++)
    {
     if (desc_lst[i].store_in_checkpoint()==true) {
      int time_order = parent->Time_blockingFactor();
      state[i].restart(
        time_order,
        level_slab_dt_type,
        level_MAX_NUM_SLAB,
        level,
        is, 
        geom.Domain(),
        grids,
        dmap,
        desc_lst[i],
        desc_lstGHOST[i],
        papa.theRestartFile());
     } else if (desc_lst[i].store_in_checkpoint()==false) {
      Real time=parent->cumTime();
      Real dt=parent->getDt();
      int time_order = parent->Time_blockingFactor();
      state[i].define(
        level,
        geom.Domain(),
        grids,
        dmap,
        desc_lst[i],
        desc_lstGHOST[i],
        time, 
        dt,
        time_order,
        level_slab_dt_type,
        level_MAX_NUM_SLAB);
     } else
      amrex::Error("store_in_checkpoint invalid"); 
    }

    finishConstructor();
}

void
AmrLevel::finishConstructor () {

    //
    // Set physical locations of grids.
    //
    grid_loc.resize(grids.size());

    for (int i = 0; i < grid_loc.size(); i++)
    { 
        grid_loc[i] = RealBox(grids[i],geom.CellSize(),geom.ProbLo());
    }

}

void
AmrLevel::setTimeLevel (Real time,Real& dt)
{
    for (int k = 0; k < desc_lst.size(); k++)
    {
        state[k].setTimeLevel(time,dt);
    }
}

bool
AmrLevel::isStateVariable (const std::string& name,
                           int&           typ,
                           int&            n)
{
    for (typ = 0; typ < desc_lst.size(); typ++)
    {
        const StateDescriptor& desc = desc_lst[typ];

        for (n = 0; n < desc.nComp(); n++)
        {
            if (desc.name(n) == name)
                return true;
        }
    }
    return false;
}

long
AmrLevel::countCells () const
{
    long cnt = 0;
    for (int i = 0, N = grids.size(); i < N; i++)
    {
        cnt += grids[i].numPts();
    }
    return cnt;
}

// dir="ckfile" (directory)
// os=HeaderFile 
void
AmrLevel::checkPoint (const std::string& dir,
                      std::ostream&      os)
{
    int ndesc = desc_lst.size();
    int i;
    int ndesc_checkpoint=0;
    for (int j=0;j<ndesc;j++) {
     if (desc_lst[j].store_in_checkpoint()==true)
      ndesc_checkpoint++;
    } 
    //
    // Build directory to hold the MultiFabs in the StateData at this level.
    // The directory is relative the the directory containing the Header file.
    //
    // SUSSMAN
    std::string Level_string = amrex::Concatenate("Level_", level, 1);
    //
    // Now for the full pathname of that directory.
    //
    std::string FullPath = dir;
    if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/')
    {
        FullPath += '/';
    }
     //SUSSMAN
    FullPath += Level_string;
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

    if (ParallelDescriptor::IOProcessor())
    {
        os << level << '\n' << geom  << '\n';
        grids.writeOn(os);
        os << ndesc_checkpoint << '\n';
    }
    //
    // Output state data.
    //

    for (i = 0; i < ndesc; i++)
    {
        //
        // Now build the full relative pathname of the StateData.
        // The name is relative to the Header file containing this name.
        // It's the name that gets written into the Header.
        //
        // There is only one MultiFab written out at each level in HyperCLaw.
        //
        if (desc_lst[i].store_in_checkpoint()==true) {
         std::string PathNameInHdr = 
           amrex::Concatenate(Level_string + "/SD_", i, 1);
         std::string FullPathName  = 
           amrex::Concatenate(FullPath + "/SD_", i, 1);

          // os=HeaderFile 
         state[i].checkPoint(PathNameInHdr, FullPathName, os);
        }
    }
}

AmrLevel::~AmrLevel ()
{
    parent = 0;
}

void
AmrLevel::FillPatch (AmrLevel & old,
                     MultiFab& mf,
                     int       dcomp,
                     Real      time,
                     int       index,
                     int       scomp,
                     int       ncomp)
{
 BL_PROFILE("AmrLevel::FillPatch()");

 //
 // Must fill this region on crse level and interpolate.
 //
 if (level<0)
  amrex::Error("level invalid in FillPatch");
 BL_ASSERT(ncomp <= (mf.nComp()-dcomp));
 BL_ASSERT(0 <= index && index < desc_lst.size());

 Vector<int> scompBC_map;
 scompBC_map.resize(ncomp);
 for (int i=0;i<ncomp;i++)
  scompBC_map[i]=scomp+i;

 int                     DComp   = dcomp;
 const StateDescriptor&  desc    = old.desc_lst[index];
 int bfact_fine=parent->Space_blockingFactor(level);

 desc.check_inRange(scompBC_map,ncomp);

 std::vector< std::pair<int,int> > ranges = 
   desc.sameInterps(scompBC_map,ncomp);

 for (int i = 0; i < ranges.size(); i++) {
  const int     scomp_range = ranges[i].first;
  const int     ncomp_range = ranges[i].second;
  Interpolater* mapper = desc.interp(scompBC_map[scomp_range]);

  Real nudge_time;
  int best_index;
  StateData& fstatedata = old.state[index];
  fstatedata.get_time_index(time,nudge_time,best_index);

  MultiFab& fmf=fstatedata.newData(best_index);
  const Geometry& fgeom = old.geom;
  StateDataPhysBCFunct fbc(fstatedata,fgeom);

  Vector<int> local_scompBC_map;
  local_scompBC_map.resize(ncomp_range);
  for (int i=0;i<ncomp_range;i++)
   local_scompBC_map[i]=scompBC_map[scomp_range+i];

  int scomp_local=scomp_range+scomp;

  if (level==0) {

   amrex::FillPatchSingleLevel(
    level,
    mf,
    nudge_time,
    fmf,
    scomp_local,
    DComp,
    ncomp_range,
    fgeom,
    fbc,
    local_scompBC_map,
    bfact_fine);
  } else if (level>0) {

   AmrLevel&               clev    = parent->getLevel(level-1);
   const Geometry&         cgeom   = clev.geom;
   int bfact_coarse=parent->Space_blockingFactor(level-1);
   StateData& cstatedata = clev.state[index];
   MultiFab& cmf=cstatedata.newData(best_index);
   StateDataPhysBCFunct cbc(cstatedata,cgeom);

   amrex::FillPatchTwoLevels(
    mf,
    nudge_time,
    cmf,
    fmf,
    scomp_local,
    DComp,
    ncomp_range,
    cgeom,
    fgeom,
    cbc,
    fbc,
    mapper,
    desc.getBCs(),  // global_bcs
    local_scompBC_map,
    level-1,level,
    bfact_coarse,bfact_fine);  
     
  } else 
   amrex::Error("level invalid");

  DComp += ncomp_range;

 } // i=0..ranges.size()-1

}   // FillPatch


void
AmrLevel::FillCoarsePatchGHOST (
                           MultiFab& cmf, // source (scomp..scomp+ncomp-1)
                           MultiFab& mf,  // dest (scomp..scomp+ncomp-1)
                           Real      time,
                           int       index,
                           int       scomp, //cmf: scomp..scomp+ncomp-1
                           Vector<int> scompBC_map, // 0..ncomp-1
                           int       ncomp)
{
 BL_PROFILE("AmrLevel::FillCoarsePatchGHOST()");

 if (level<=0)
  amrex::Error("level invalid in FillCoarsePatchGHOST");

 if (ncomp+scomp>mf.nComp())
  amrex::Error("ncomp+scomp>mf.nComp()");
 if (ncomp+scomp>cmf.nComp())
  amrex::Error("ncomp+scomp>cmf.nComp()");
 if ((index<0)||(index>=desc_lst.size()))
  amrex::Error("(index<0)||(index>=desc_lst.size())");
 if ((index<0)||(index>=desc_lstGHOST.size()))
  amrex::Error("(index<0)||(index>=desc_lstGHOST.size())");
 if (scompBC_map.size()!=ncomp)
  amrex::Error("scompBC_map has invalid size");

 int ngrow=mf.nGrow();
 const BoxArray& mf_BA = mf.boxArray();
 if (ngrow<0)
  amrex::Error("ngrow<0 in FillCoarsePatchGHOST");

 DistributionMapping dm=mf.DistributionMap();

 const BoxArray& cmf_BA=cmf.boxArray();
 DistributionMapping cdm=cmf.DistributionMap();

  // because of the proper nesting requirement, no ghost
  // values are needed for cmf_part.
  // cmf_part: 0..ncomp-1
 MultiFab* cmf_part=new MultiFab(cmf_BA,cdm,ncomp,0,Fab_allocate);
 MultiFab::Copy(*cmf_part,cmf,scomp,0,ncomp,0);

 int                     DComp   = scomp;
 const StateDescriptor&  descGHOST = desc_lstGHOST[index];
 int bfact_fine=parent->Space_blockingFactor(level);
 int bfact_coarse=parent->Space_blockingFactor(level-1);

  // scompBC_map: 0...ncomp-1
 descGHOST.check_inRange(scompBC_map, ncomp);
 std::vector< std::pair<int,int> > ranges = 
   descGHOST.sameInterps(scompBC_map,ncomp);

 Real nudge_time;
 int best_index;
 StateData& fstatedata = state[index];
 fstatedata.get_time_index(time,nudge_time,best_index);

 AmrLevel&               clev    = parent->getLevel(level-1);
 const Geometry&         cgeom   = clev.geom;
 StateData& cstatedata = clev.state[index];
 Real c_nudge_time;
 int c_best_index;
 cstatedata.get_time_index(time,c_nudge_time,c_best_index);
 StateDataPhysBCFunctGHOST physbc_coarse(cstatedata,cgeom);

 const Box& pdomain = state[index].getDomain();

 for (int i = 0; i < ranges.size(); i++) {

  const int     scomp_range  = ranges[i].first;
  const int     ncomp_range  = ranges[i].second;
  if (i==0) {
   if (scomp_range!=0)
    amrex::Error("scomp_range!=0");
  }
  Interpolater* mapper = descGHOST.interp(scompBC_map[scomp_range]);

  BoxArray crseBA(mf_BA.size());
  for (int j = 0, N_CBA = crseBA.size(); j < N_CBA; ++j) {
   BL_ASSERT(mf_BA[j].ixType() == descGHOST.getType());
   const Box& bx = mf_BA[j];
   crseBA.set(j,mapper->CoarseBox(bx,bfact_coarse,bfact_fine));
  }

   // ghost cells do not have to be initialized.
   // call InterpBordersGHOST after FillCoarsePatchGHOST.
  MultiFab crseMF(crseBA,dm,ncomp_range,0,Fab_allocate);

  int scomp_data=scomp_range;
  int dcomp_data=scomp+scomp_data;

  if (dcomp_data!=DComp)
   amrex::Error("dcomp_data!=DComp");

  Vector<int> local_scompBC_map;
  local_scompBC_map.resize(ncomp_range);
  for (int i=0;i<ncomp_range;i++)
   local_scompBC_map[i]=scompBC_map[scomp_range+i];

  int scomp_local=scomp_range;

  const MultiFab& cmf_part_mf=*cmf_part;

  amrex::FillPatchSingleLevel(
    level-1,
    crseMF, // data to be filled 0..ncomp_range-1
    c_nudge_time,
    cmf_part_mf, // level-1 data; scomp_range..scomp_range+ncomp_range-1
    scomp_local, // = scomp_range
    0, // dstcomp
    ncomp_range,
    cgeom,
    physbc_coarse,
    local_scompBC_map,
    bfact_coarse);

  Vector< BCRec > local_bcs;
  const Vector< BCRec> & global_bcs=descGHOST.getBCs();
  local_bcs.resize(ncomp_range);
  for (int i=0;i<ncomp_range;i++)
   local_bcs[i]=global_bcs[local_scompBC_map[i]]; 

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {

   const Box& dbx = mfi.validbox();
	    
   Vector<BCRec> bcr(ncomp_range);
   int src_comp_bcs=0;
   int dest_comp_bcr=0;
     // source: local_bcs  dest: bcr
   amrex::setBC(dbx,pdomain,src_comp_bcs,dest_comp_bcr,ncomp_range,
     local_bcs,bcr);
	   
   mapper->interp(nudge_time,
                  crseMF[mfi],
                  0,  // crse_comp
		  mf[mfi],
		  DComp,
		  ncomp_range,
		  dbx,
		  cgeom,
		  geom,
		  bcr, // not used.
                  level-1,level,
		  bfact_coarse,bfact_fine);
  }  // mfi

  DComp += ncomp_range;

 } // i=0..ranges.size()-1

 delete cmf_part;

}   // FillCoarsePatchGHOST

void
AmrLevel::InterpBordersGHOST (
                     MultiFab& cmf,
                     MultiFab& mf,
                     Real      time,
                     int       index,
                     int       scomp, // source comp wrt mf
                     Vector<int> scompBC_map,
                     int       ncomp)
{
 BL_PROFILE("AmrLevel::InterpBordersGHOST()");

 if (level<0)
  amrex::Error("level invalid in InterpBordersGHOST");

 BL_ASSERT(ncomp <= (mf.nComp()-scomp));
 BL_ASSERT(0 <= index && index < desc_lst.size());
 BL_ASSERT(0 <= index && index < desc_lstGHOST.size());
 if (scompBC_map.size()!=ncomp)
  amrex::Error("scompBC_map has invalid size");

 int ngrow=mf.nGrow();
 const BoxArray& mf_BA = mf.boxArray();

 if (ngrow<=0)
  amrex::Error("ngrow<=0 in InterpBordersGHOST");

 DistributionMapping dm=mf.DistributionMap();

 MultiFab fmf(mf_BA,dm,ncomp,0,Fab_allocate);

  // dstmf,srcmf,srccomp,dstcomp,ncomp,ngrow
 MultiFab::Copy(fmf,mf,scomp,0,ncomp,0);

 MultiFab* cmf_part;

 if (level>0) {
  const BoxArray& cmf_BA=cmf.boxArray();
  DistributionMapping cdm=cmf.DistributionMap();
  cmf_part=new MultiFab(cmf_BA,cdm,ncomp,0,Fab_allocate);
  MultiFab::Copy(*cmf_part,cmf,scomp,0,ncomp,0);
 }  // level>0
 
 int                     DComp   = scomp;
 const StateDescriptor&  descGHOST = desc_lstGHOST[index];
 int bfact_fine=parent->Space_blockingFactor(level);

 descGHOST.check_inRange(scompBC_map, ncomp);
 std::vector< std::pair<int,int> > ranges = 
   descGHOST.sameInterps(scompBC_map,ncomp);

 for (int i = 0; i < ranges.size(); i++) {
  const int     scomp_range  = ranges[i].first;
  const int     ncomp_range  = ranges[i].second;
  Interpolater* mapper = descGHOST.interp(scompBC_map[scomp_range]);

  Real nudge_time;
  int best_index;
  StateData& fstatedata = state[index];
  fstatedata.get_time_index(time,nudge_time,best_index);

  const Geometry& fgeom = geom;
  StateDataPhysBCFunctGHOST fbc(fstatedata,fgeom);

  int scomp_data=scomp_range;
  int dcomp_data=scomp+scomp_data;

  if (dcomp_data!=DComp)
   amrex::Error("dcomp_data!=DComp");

  Vector<int> local_scompBC_map;
  local_scompBC_map.resize(ncomp_range);
  for (int i=0;i<ncomp_range;i++)
   local_scompBC_map[i]=scompBC_map[scomp_range+i];

  if (level==0) {
   amrex::FillPatchSingleLevel(
    level,
    mf,
    nudge_time,
    fmf,
    scomp_data,
    dcomp_data,
    ncomp_range,
    fgeom,
    fbc,
    local_scompBC_map,
    bfact_fine);
  } else if (level>0) {

   AmrLevel&               clev    = parent->getLevel(level-1);
   const Geometry&         cgeom   = clev.geom;
   int bfact_coarse=parent->Space_blockingFactor(level-1);
   StateData& cstatedata = clev.state[index];
   StateDataPhysBCFunctGHOST cbc(cstatedata,cgeom);

   amrex::FillPatchTwoLevels(
    mf,
    nudge_time,
    *cmf_part,
    fmf,
    scomp_data,
    dcomp_data,
    ncomp_range,
    cgeom,
    fgeom,
    cbc,
    fbc,
    mapper,
    descGHOST.getBCs(), // global_bcs
    local_scompBC_map,
    level-1,level,
    bfact_coarse,bfact_fine);  
     
  } else 
   amrex::Error("level invalid");

  DComp += ncomp_range;

 } // i=0..ranges.size()-1

 if (level>0) 
  delete cmf_part;

}   // InterpBordersGHOST


void
AmrLevel::InterpBorders (
                     MultiFab& cmf,
                     MultiFab& mf,
                     Real      time,
                     int       index,
                     int       scomp,
                     Vector<int> scompBC_map,
                     int       ncomp)
{
 BL_PROFILE("AmrLevel::InterpBorders()");

 if (level<0)
  amrex::Error("level invalid in InterpBorders");

 BL_ASSERT(ncomp <= (mf.nComp()-scomp));
 BL_ASSERT(0 <= index && index < desc_lst.size());
 BL_ASSERT(0 <= index && index < desc_lstGHOST.size());
 if (scompBC_map.size()!=ncomp)
  amrex::Error("scompBC_map has invalid size");

 int ngrow=mf.nGrow();
 const BoxArray& mf_BA = mf.boxArray();

 if (ngrow<=0)
  amrex::Error("ngrow<=0 in InterpBorders");

 DistributionMapping dm=mf.DistributionMap();

 MultiFab fmf(mf_BA,dm,ncomp,0,Fab_allocate);

  // dstmf,srcmf,srccomp,dstcomp,ncomp,ngrow
 MultiFab::Copy(fmf,mf,scomp,0,ncomp,0);

 MultiFab* cmf_part;

 if (level>0) {
  const BoxArray& cmf_BA=cmf.boxArray();
  DistributionMapping cdm=cmf.DistributionMap();
  cmf_part=new MultiFab(cmf_BA,cdm,ncomp,0,Fab_allocate);
  MultiFab::Copy(*cmf_part,cmf,scomp,0,ncomp,0);
 }  // level>0
 
 int                     DComp   = scomp;
 const StateDescriptor&  desc    = desc_lst[index];
 int bfact_fine=parent->Space_blockingFactor(level);

 desc.check_inRange(scompBC_map, ncomp);
 std::vector< std::pair<int,int> > ranges = 
   desc.sameInterps(scompBC_map,ncomp);

 for (int i = 0; i < ranges.size(); i++) {
  const int     scomp_range  = ranges[i].first;
  const int     ncomp_range  = ranges[i].second;
  Interpolater* mapper = desc.interp(scompBC_map[scomp_range]);

  Real nudge_time;
  int best_index;
  StateData& fstatedata = state[index];
  fstatedata.get_time_index(time,nudge_time,best_index);

  const Geometry& fgeom = geom;
  StateDataPhysBCFunct fbc(fstatedata,fgeom);

  int scomp_data=scomp_range;
  int dcomp_data=scomp+scomp_data;

  if (dcomp_data!=DComp)
   amrex::Error("dcomp_data!=DComp");

  Vector<int> local_scompBC_map;
  local_scompBC_map.resize(ncomp_range);
  for (int i=0;i<ncomp_range;i++)
   local_scompBC_map[i]=scompBC_map[scomp_range+i];

  if (level==0) {
   amrex::FillPatchSingleLevel(
    level,
    mf,
    nudge_time,
    fmf,
    scomp_data,
    dcomp_data,
    ncomp_range,
    fgeom,
    fbc,
    local_scompBC_map,
    bfact_fine);
  } else if (level>0) {

   AmrLevel&               clev    = parent->getLevel(level-1);
   const Geometry&         cgeom   = clev.geom;
   int bfact_coarse=parent->Space_blockingFactor(level-1);
   StateData& cstatedata = clev.state[index];
   StateDataPhysBCFunct cbc(cstatedata,cgeom);

   amrex::FillPatchTwoLevels(
    mf,
    nudge_time,
    *cmf_part,
    fmf,
    scomp_data,
    dcomp_data,
    ncomp_range,
    cgeom,
    fgeom,
    cbc,
    fbc,
    mapper,
    desc.getBCs(), // global_bcs
    local_scompBC_map,
    level-1,level,
    bfact_coarse,bfact_fine);  
     
  } else 
   amrex::Error("level invalid");

  DComp += ncomp_range;
 } // i=0..ranges.size()-1

 if (level>0) 
  delete cmf_part;

}   // InterpBorders


void
AmrLevel::FillCoarsePatch (MultiFab& mf,
                           int       dcomp, // update mf: dcomp..dcomp+ncomp-1
                           Real      time,
                           int       index,
                           int       scomp, // scompBC_map[i]=scomp+i
                           int       ncomp)
{
 BL_PROFILE("AmrLevel::FillCoarsePatch()");

 //
 // Must fill this region on crse level and interpolate.
 //
 if (level<=0)
  amrex::Error("level invalid in FillCoarsePatch");
 if (ncomp>(mf.nComp()-dcomp))
  amrex::Error("ncomp>(mf.nComp()-dcomp)");
 if ((index<0)||(index>=desc_lst.size()))
  amrex::Error("(index<0)||(index>=desc_lst.size())");

 Vector<int> scompBC_map;
 scompBC_map.resize(ncomp);
 for (int i=0;i<ncomp;i++)
  scompBC_map[i]=scomp+i;

 AmrLevel&               clev    = parent->getLevel(level-1);
 const Geometry&         cgeom   = clev.geom;
 int bfact_coarse=parent->Space_blockingFactor(level-1);
 StateData& cstatedata = clev.state[index];

 Real nudge_time;
 int best_index;
 cstatedata.get_time_index(time,nudge_time,best_index);

 MultiFab& cmf=cstatedata.newData(best_index);
 StateDataPhysBCFunct physbc_coarse(cstatedata,cgeom);

 int                     DComp   = dcomp;
 const StateDescriptor&  desc    = desc_lst[index];
 const Box&              pdomain = state[index].getDomain();
 const BoxArray&         mf_BA   = mf.boxArray();
 DistributionMapping dm=mf.DistributionMap();

 int bfact_fine=parent->Space_blockingFactor(level);

 desc.check_inRange(scompBC_map, ncomp);

 std::vector< std::pair<int,int> > ranges = 
   desc.sameInterps(scompBC_map,ncomp);

 for (int i = 0; i < ranges.size(); i++) {
  const int     scomp_range  = ranges[i].first;
  const int     ncomp_range  = ranges[i].second;
  Interpolater* mapper = desc.interp(scompBC_map[scomp_range]);

  BoxArray crseBA(mf_BA.size());
   
  for (int j = 0, N = crseBA.size(); j < N; ++j) {
   BL_ASSERT(mf_BA[j].ixType() == desc.getType());
   const Box& bx = mf_BA[j];
   crseBA.set(j,mapper->CoarseBox(bx,bfact_coarse,bfact_fine));
  }

    // ngrow=0
  MultiFab crseMF(crseBA,dm,ncomp_range,0,Fab_allocate);

  Vector<int> local_scompBC_map;
  local_scompBC_map.resize(ncomp_range);
  for (int i=0;i<ncomp_range;i++)
   local_scompBC_map[i]=scompBC_map[scomp_range+i];

  int scomp_local=scomp_range+scomp;

  amrex::FillPatchSingleLevel(
    level-1,
    crseMF, // data to be filled 0..ncomp_range-1
    nudge_time,
    cmf,  // level-1 data (smf); scomp_local...scomp_local+ncomp_range-1
    scomp_local,
    0, // dstcomp
    ncomp_range,
    cgeom,
    physbc_coarse,
    local_scompBC_map,
    bfact_coarse);

  Vector< BCRec > local_bcs;
  const Vector< BCRec> & global_bcs=desc.getBCs();

  local_bcs.resize(ncomp_range);
  for (int i=0;i<ncomp_range;i++)
   local_bcs[i]=global_bcs[local_scompBC_map[i]]; 

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(mf); mfi.isValid(); ++mfi) {

   const Box& dbx = mfi.validbox();
	    
   Vector<BCRec> bcr(ncomp_range);
   int src_comp_bcs=0;
   int dest_comp_bcr=0;
     // source: local_bcs  dest: bcr
   amrex::setBC(dbx,pdomain,src_comp_bcs,dest_comp_bcr,ncomp_range,
     local_bcs,bcr);
	   
   mapper->interp(nudge_time,
                  crseMF[mfi],
                  0,  // crse_comp
		  mf[mfi],
		  DComp,
		  ncomp_range,
		  dbx,
		  cgeom,
		  geom,
		  bcr,
                  level-1,level,
		  bfact_coarse,bfact_fine);
  }  // mfi

  StateDataPhysBCFunct physbc_fine(state[index],geom);
  physbc_fine.FillBoundary(level,mf,nudge_time,DComp,
    local_scompBC_map,ncomp_range,bfact_fine);

  DComp += ncomp_range;
 } // i=0..ranges.size()-1

}   // FillCoarsePatch

Vector<int>
AmrLevel::getBCArray (int State_Type,
                      int gridno,
                      int strt_comp,
                      int ncomp)
{
    Vector<int> bc(2*BL_SPACEDIM*ncomp);

    for (int n = 0; n < ncomp; n++)
    {
        const int* b_rec = state[State_Type].getBC(strt_comp+n,gridno).vect();
        for (int m = 0; m < 2*BL_SPACEDIM; m++)
            bc[2*BL_SPACEDIM*n + m] = b_rec[m];
    }

    return bc;
} // getBCArray

Vector<int>
AmrLevel::getBCArrayGHOST (int State_Type,
                      int gridno,
                      int strt_comp,
                      int ncomp)
{
 Vector<int> bc(2*BL_SPACEDIM*ncomp);

 for (int n = 0; n < ncomp; n++) {
  const int* b_rec = state[State_Type].getBCGHOST(strt_comp+n,gridno).vect();
  for (int m = 0; m < 2*BL_SPACEDIM; m++)
   bc[2*BL_SPACEDIM*n + m] = b_rec[m];
 } // n

 return bc;
} // getBCArrayGHOST

void
AmrLevel::setPlotVariables ()
{
    ParmParse pp("amr");

    if (pp.contains("plot_vars"))
    {
        std::string nm;
      
        int nPltVars = pp.countval("plot_vars");
      
        for (int i = 0; i < nPltVars; i++)
        {
            pp.get("plot_vars", nm, i);

            if (nm == "ALL") 
                parent->fillStatePlotVarList();
            else if (nm == "NONE")
                parent->clearStatePlotVarList();
            else
                parent->addStatePlotVar(nm);
        }
    }
    else 
    {
        //
        // The default is to add them all.
        //
        parent->fillStatePlotVarList();
    }
  
} // subroutine setPlotVariables

} // namespace amrex

