
#include <sstream>

#include <unistd.h>
#include <memory>
#include <limits>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_ParmParse.H>
#include <AmrLevel.H>
#include <FillPatchUtil.H>

namespace amrex {

DescriptorList AmrLevel::desc_lst;
DescriptorList AmrLevel::desc_lstGHOST;

void FSI_container_class::initData_FSI(
  int num_nodes_init,
  int num_elements_init) {

if ((num_nodes_init>=0)&&(num_elements_init>=0)) {
 num_nodes=num_nodes_init;
 num_elements=num_elements_init;

 node_list.resize(num_nodes*3);
 element_list.resize(num_elements*4);
 displacement_list.resize(num_nodes*3);
 velocity_halftime_list.resize(num_nodes*3);
 velocity_list.resize(num_nodes*3);
 force_list.resize(num_nodes*3);
 mass_list.resize(num_nodes);
 temperature_list.resize(num_nodes);
} else
 amrex::Error("num_nodes_init or num_elements_init invalid");

} // end subroutine initData_FSI()

void FSI_container_class::copyFrom_FSI(const FSI_container_class& source_FSI) {

 initData_FSI(source_FSI.num_nodes,source_FSI.num_elements);

 for (int ielem=0;ielem<4*num_elements;ielem++) {
  element_list[ielem]=source_FSI.element_list[ielem];
 } //ielem=0..4*num_elements-1

 for (int inode=0;inode<num_nodes;inode++) {
  mass_list[inode]=source_FSI.mass_list[inode];
  temperature_list[inode]=source_FSI.temperature_list[inode];
 }

 for (int inode=0;inode<3*num_nodes;inode++) {
  node_list[inode]=source_FSI.node_list[inode];
  displacement_list[inode]=source_FSI.displacement_list[inode];
  velocity_list[inode]=source_FSI.velocity_list[inode];
  velocity_halftime_list[inode]=source_FSI.velocity_halftime_list[inode];
  force_list[inode]=source_FSI.force_list[inode];
 } //inode=0..3*num_nodes-1

} // end subroutine copyFrom_FSI

void FSI_container_class::clear_FSI() {

 initData_FSI(0,0);

} // end subroutine FSI_container_class::clear_FSI() 



void
AmrLevel::manual_tags_placement (TagBoxArray&    tags,
             const Vector<IntVect>& bf_lev)
{}

const BoxArray& AmrLevel::getAreaNotToTag () noexcept
{
    return m_AreaNotToTag;
}


AmrLevel::AmrLevel () noexcept
{
   parent = 0;
   level = -1;
}

// constructor
AmrLevel::AmrLevel (AmrCore&        papa,
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

    int max_level=parent->maxLevel();

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

    int time_order=parent->Time_blockingFactor();
    int local_particles_flag=parent->global_AMR_particles_flag;
    int local_nmat=parent->global_AMR_num_materials;
    int num_SoA_var=parent->global_AMR_num_SoA_var;

    if (level==0) {

     AmrLevel0_new_dataPC.resize(level_MAX_NUM_SLAB);
     new_data_FSI.resize(level_MAX_NUM_SLAB);

     for (int i=0;i<=time_order;i++) {

      new_data_FSI[i].resize(local_nmat);
      for (int j=0;j<local_nmat;j++) {
       new_data_FSI[i][j].initData_FSI(0,0);
      }

      if (local_particles_flag==1) {

       if (num_SoA_var>=0) {
        AmrLevel0_new_dataPC[i]=
          new AmrParticleContainer<N_EXTRA_REAL,N_EXTRA_INT,0,0>(parent);
        //add Structure of Array component(s)
        //amrex-master/Src/Particle/AMReX_Particles.H
        //void AddRealComp (T communicate=true)
        for (int ns=0;ns<num_SoA_var;ns++)
         AmrLevel0_new_dataPC[i]->AddRealComp(true);
       } else
        amrex::Error("num_SoA_var invalid");

      } else if (local_particles_flag==0) {
       // do nothing
      } else
       amrex::Error("local_particles_flag invalid");

     }// for (int i=0;i<=time_order;i++) 
    } else if (level>0) {
     // do nothing
    } else
     amrex::Error("level invalid");

    for (int icomp = 0; icomp < state.size(); icomp++)
    {
        state[icomp].define(
			papa,
			level,
			max_level,
                        geom.Domain(),
                        grids,
                        dm,
                        desc_lst[icomp],
                        desc_lstGHOST[icomp],
                        time,
                        dt,
                        time_order,
                        level_slab_dt_type,
                        level_MAX_NUM_SLAB);
    }

    finishConstructor();
}



void
AmrLevel::restart (AmrCore&      papa,
                   std::istream& is,
		   int old_finest_level,
		   int new_finest_level) {

    parent = &papa;

    int max_level=parent->maxLevel();

    is >> level;
    is >> geom;

    grids.readFrom(is);

    int nstate;
    is >> nstate;

    int ndesc = desc_lst.size();

    if (nstate!=ndesc)
     amrex::Error("nstate and ndesc do not match");

    dmap.define(grids);

    parent->SetBoxArray(level, grids);
    parent->SetDistributionMap(level, dmap);
    int level_MAX_NUM_SLAB=parent->get_MAX_NUM_SLAB();
    int level_slab_dt_type=parent->get_slab_dt_type();

      // e.g. chkfile=./chk<nsteps>
    std::string chkfile=papa.theRestartFile();
    std::string Level_string = amrex::Concatenate("Level_", level, 1);
    std::string FullPath = chkfile;
    if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/') {
     FullPath += '/';
    }
    FullPath += Level_string;

    if (level==0) {

     // SUSSMAN restart the particle data
     std::string FullPathName=FullPath;
     int time_order=parent->Time_blockingFactor();
     int local_particles_flag=parent->global_AMR_particles_flag;
     int local_nmat=parent->global_AMR_num_materials;

     AmrLevel0_new_dataPC.resize(level_MAX_NUM_SLAB);
     new_data_FSI.resize(level_MAX_NUM_SLAB);

     for (int i=0;i<=time_order;i++) {

       //TODO: restart the FSI data.
      new_data_FSI[i].resize(local_nmat);
      for (int j=0;j<local_nmat;j++) {
       new_data_FSI[i][j].initData_FSI(0,0);
      }

      if (local_particles_flag==0) {
       // do nothing
      } else if (local_particles_flag==1) {

       //The actual particle data will be read from the restart file
       //later (in NavierStokes.cpp) when "post_restart" is called.
       //Particle data needs the AMR hierarchy to already be built before
       //being read from the restart file and redistributed appropriately.
       
      } else
       amrex::Error("local_particles_flag invalid");

     }//for (int i=0;i<=time_order;i++) 

    } else if (level>0) {
     // do nothing
    } else
     amrex::Error("level invalid");

    state.resize(ndesc);
    for (int icomp = 0; icomp < ndesc; icomp++)
    {
     int time_order = parent->Time_blockingFactor();
     state[icomp].restart(
        papa,
        time_order,
        level_slab_dt_type,
        level_MAX_NUM_SLAB,
        level,
	max_level,
        is, 
        geom.Domain(),
        grids,
        dmap,
        desc_lst[icomp],
        desc_lstGHOST[icomp],
        papa.theRestartFile());
    }

    finishConstructor();
}  // end subroutine AmrLevel::restart

void
AmrLevel::finishConstructor () {

    //
    // Set physical locations of grids.
    //
    grid_loc.resize(grids.size());

    for (int igrid = 0; igrid < grid_loc.size(); igrid++)
    { 
        grid_loc[igrid] = RealBox(grids[igrid],geom.CellSize(),geom.ProbLo());
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

Long
AmrLevel::countCells () const
{
    Long cnt = 0;
    for (int igrid = 0, N = grids.size(); igrid < N; igrid++)
    {
        cnt += grids[igrid].numPts();
    }
    return cnt;
}

// dir="ckfile" (directory)
// os=HeaderFile 
void
AmrLevel::checkPoint (const std::string& dir,
                      std::ostream&      os)
{
    int max_level=parent->maxLevel();

    int ndesc = desc_lst.size();

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
    if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/') {
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
        os << ndesc << '\n';
    }

    if (level==0) {

     // SUSSMAN checkpoint the particle data
     std::string FullPathName=FullPath;
     int time_order=parent->Time_blockingFactor();
     int local_particles_flag=parent->global_AMR_particles_flag;

     for (int i=0;i<=time_order;i++) {

      if (local_particles_flag==0) {
       // do nothing
      } else if (local_particles_flag==1) {

       std::string Part_name="FusionPart";
       std::stringstream i_string_stream(std::stringstream::in |
           std::stringstream::out);
       i_string_stream << i;
       std::string i_string=i_string_stream.str();
       Part_name+=i_string;
       AmrLevel0_new_dataPC[i]->Checkpoint(FullPathName,Part_name);

      } else
       amrex::Error("local_particles_flag invalid");

     }//for (int i=0;i<=time_order;i++) 

    } else if (level>0) {
     // do nothing
    } else
     amrex::Error("level invalid");

    //
    // Output state data.
    //

    for (int icomp = 0; icomp < ndesc; icomp++)
    {
        //
        // Now build the full relative pathname of the StateData.
        // The name is relative to the Header file containing this name.
        // It's the name that gets written into the Header.
        //
        // There is only one MultiFab written out at each level in HyperCLaw.
        //
        std::string PathNameInHdr = 
          amrex::Concatenate(Level_string + "/SD_", icomp, 1);
        std::string FullPathName  = 
          amrex::Concatenate(FullPath + "/SD_", icomp, 1);

         // os=HeaderFile 
        state[icomp].checkPoint(PathNameInHdr, FullPathName, os,
			level,max_level);
        
    }
} // end subroutine AmrLevel::checkPoint

AmrLevel::~AmrLevel ()
{

    if (level==0) {
     int time_order=parent->Time_blockingFactor();
     int local_particles_flag=parent->global_AMR_particles_flag;
     int local_nmat=parent->global_AMR_num_materials;

     for (int i=0;i<=time_order;i++) {

      for (int j=0;j<local_nmat;j++) {
       new_data_FSI[i][j].clear_FSI();
      }
      new_data_FSI[i].resize(0);

      if (local_particles_flag==1) {
       delete AmrLevel0_new_dataPC[i];
      } else if (local_particles_flag==0) {
       // do nothing
      } else
       amrex::Error("local_particles_flag invalid");

     } // for (int i=0;i<=time_order;i++) 
     AmrLevel0_new_dataPC.resize(0);
     new_data_FSI.resize(0);
    } else if (level>0) {
     // do nothing
    } else
     amrex::Error("level invalid");

    parent = 0;
}

AmrParticleContainer<N_EXTRA_REAL,N_EXTRA_INT,0,0>& 
   AmrLevel::newDataPC (int slab_index)
{

if (level==0) {
 int time_order=parent->Time_blockingFactor();
 int local_particles_flag=parent->global_AMR_particles_flag;
 if (local_particles_flag==1) {

  int project_slab_index=slab_index;
  if (project_slab_index==-1)
   project_slab_index=0;
  if (project_slab_index==time_order+1)
   project_slab_index=time_order;
  if ((project_slab_index<0)||
      (project_slab_index>time_order)) {
   std::cout << "time_order= " << time_order << '\n';
   std::cout << "project_slab_index= " << project_slab_index << '\n';
   amrex::Error("project_slab_index invalid1");
  }

  return *AmrLevel0_new_dataPC[project_slab_index];

 } else {
  amrex::Error("local_particles_flag invalid"); 
  return *AmrLevel0_new_dataPC[0];
 }
} else {
 amrex::Error("level invalid"); 
 return *AmrLevel0_new_dataPC[0];
}

}

const AmrParticleContainer<N_EXTRA_REAL,N_EXTRA_INT,0,0>&
AmrLevel::newDataPC (int slab_index) const
{

if (level==0) {
 int time_order=parent->Time_blockingFactor();
 int local_particles_flag=parent->global_AMR_particles_flag;
 if (local_particles_flag==1) {

  int project_slab_index=slab_index;
  if (project_slab_index==-1)
   project_slab_index=0;
  if (project_slab_index==time_order+1)
   project_slab_index=time_order;
  if ((project_slab_index<0)||
      (project_slab_index>time_order)) {
   std::cout << "time_order= " << time_order << '\n';
   std::cout << "project_slab_index= " << project_slab_index << '\n';
   amrex::Error("project_slab_index invalid2");
  }

  return *AmrLevel0_new_dataPC[project_slab_index];
 } else {
  amrex::Error("local_particles_flag invalid"); 
  return *AmrLevel0_new_dataPC[0];
 }
} else {
 amrex::Error("level invalid"); 
 return *AmrLevel0_new_dataPC[0];
}

}

void
AmrLevel::CopyNewToOldPC(int lev_max) {

if (level==0) {

 int time_order=parent->Time_blockingFactor();
 int local_particles_flag=parent->global_AMR_particles_flag;
 int local_nmat=parent->global_AMR_num_materials;

 for (int i=0;i<time_order;i++) {

  for (int j=0;j<local_nmat;j++) {
   new_data_FSI[i][j].copyFrom_FSI(new_data_FSI[time_order][j]);
  }

  if (local_particles_flag==0) {
   // do nothing
  } else if (local_particles_flag==1) {

   //amrex-master/Src/Particle/AMReX_Particles.H
   //void copyParticles (const ParticleContainerType& other,bool local=false);
   bool local=true;  
   AmrLevel0_new_dataPC[i]->clearParticles();
    //make sure hierarchy is initialized.
   AmrLevel0_new_dataPC[i]->Redistribute();

   AmrLevel0_new_dataPC[i]->copyParticles(
        *AmrLevel0_new_dataPC[time_order],local);

   int lev_min=0;
   int nGrow_Redistribute=0;
   int local_Redistribute=0;
   //amrex-master/Src/Particle/AMReX_Particles.H
   //void Redistribute (int lev_min = 0, int lev_max = -1, int nGrow = 0, int local=0);
   AmrLevel0_new_dataPC[i]->
      Redistribute(lev_min,lev_max,nGrow_Redistribute,local_Redistribute);
  } else 
   amrex::Error("local_particles_flag invalid");
 } //i=0;i<time_order;i++
} else
 amrex::Error("level invalid in CopyNewToOldPC");

} // end subroutine CopyNewToOldPC()


void
AmrLevel::CopyOldToNewPC(int lev_max) {

if (level==0) {

 int time_order=parent->Time_blockingFactor();
 int local_particles_flag=parent->global_AMR_particles_flag;
 int local_nmat=parent->global_AMR_num_materials;

 for (int i=1;i<=time_order;i++) {

  for (int j=0;j<local_nmat;j++) {
   new_data_FSI[i][j].copyFrom_FSI(new_data_FSI[0][j]);
  }

  if (local_particles_flag==0) {
   // do nothing
  } else if (local_particles_flag==1) {

   //amrex-master/Src/Particle/AMReX_Particles.H
   //void copyParticles (const ParticleContainerType& other,bool local=false);
   bool local=true;
   AmrLevel0_new_dataPC[i]->clearParticles();
    //make sure hierarchy is initialized.
   AmrLevel0_new_dataPC[i]->Redistribute();

   AmrLevel0_new_dataPC[i]->copyParticles(
      *AmrLevel0_new_dataPC[0],local);

   int lev_min=0;
   int nGrow_Redistribute=0;
   int local_Redistribute=0;
   AmrLevel0_new_dataPC[i]->
     Redistribute(lev_min,lev_max,nGrow_Redistribute,local_Redistribute);
  } else 
   amrex::Error("local_particles_flag invalid");
 } // i=1..time_order
} else
 amrex::Error("level invalid in CopyOldToNewPC");

} // end subroutine CopyOldToNewPC()



void
AmrLevel::FillPatch (AmrLevel & old,
                     MultiFab& mf,  // data to be filled
                     int       dcomp,
                     Real      time,
                     int       index,
                     int       scomp,
                     int       ncomp,
		     int debug_fillpatch)
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
 for (int icomp=0;icomp<ncomp;icomp++)
  scompBC_map[icomp]=scomp+icomp;

 int                     DComp   = dcomp;
 const StateDescriptor&  desc    = old.desc_lst[index];
 IndexType desc_typ(desc.getType());
 int desc_grid_type=-1;
 StateData::get_grid_type(desc_typ,desc_grid_type);

 int bfact_fine=parent->Space_blockingFactor(level);

 desc.check_inRange(scompBC_map,ncomp);

 std::vector< std::pair<int,int> > ranges = 
   desc.sameInterps(scompBC_map,ncomp);

 for (unsigned int igroup = 0; igroup < ranges.size(); igroup++) {
  const int     scomp_range = ranges[igroup].first;
  const int     ncomp_range = ranges[igroup].second;
  Interpolater* mapper = desc.interp(scompBC_map[scomp_range]);

  Real nudge_time;
  int best_index;
  StateData& fstatedata = old.state[index];
   //nudge_time=time_array[best_index]
   //0<=best_index<=bfact_time_order
  fstatedata.get_time_index(time,nudge_time,best_index);

  MultiFab& fmf=fstatedata.newData(best_index);

  if (fmf.DistributionMap()==old.DistributionMap()) {
   // do nothing
  } else {
   amrex::Error("fmf.DistributionMap()!=old.DistributionMap()");
  }
  if (fmf.boxArray().CellEqual(old.boxArray())) {
   // do nothing
  } else {
   amrex::Error("fmf.boxArray().CellEqual(old.boxArray()) failed");
  }

  const Geometry& fgeom = old.geom;
  StateDataPhysBCFunct fbc(fstatedata,fgeom);

  Vector<int> local_scompBC_map;
  local_scompBC_map.resize(ncomp_range);
  for (int isub=0;isub<ncomp_range;isub++)
   local_scompBC_map[isub]=scompBC_map[scomp_range+isub];

  int scomp_local=scomp_range+scomp;

  if (level==0) {

    // This routine is NOT in the "Base" directory, instead it is in the
    // present directory in the file: FillPatchUtil.cpp.
    // We are fooled into thinking this is a amrex routine since the
    // code in FillPatchUtil.cpp is hidden within a "namespace amrex"
    // block.
   amrex::FillPatchSingleLevel(
    level,
    mf,  // data to be filled
    nudge_time,
    fmf, // old data at the current level.
    scomp_local,
    DComp,
    ncomp_range,
    fgeom,
    fbc,
    local_scompBC_map,
    bfact_fine,
    debug_fillpatch);
  } else if (level>0) {

   AmrLevel&               clev    = parent->getLevel(level-1);
   const Geometry&         cgeom   = clev.geom;
   int bfact_coarse=parent->Space_blockingFactor(level-1);
   StateData& cstatedata = clev.state[index];
   MultiFab& cmf=cstatedata.newData(best_index);

   if (cmf.DistributionMap()==clev.DistributionMap()) {
    // do nothing
   } else {
    amrex::Error("cmf.DistributionMap()!=clev.DistributionMap()");
   }
   if (cmf.boxArray().CellEqual(clev.boxArray())) {
    // do nothing
   } else {
    amrex::Error("cmf.boxArray().CellEqual(clev.boxArray()) failed");
   }
   StateDataPhysBCFunct cbc(cstatedata,cgeom);

    // This routine is NOT in the "Base" directory, instead it is in the
    // present directory in the file: FillPatchUtil.cpp.
    // We are fooled into thinking this is a amrex routine since the
    // code in FillPatchUtil.cpp is hidden within a "namespace amrex"
    // block.
   amrex::FillPatchTwoLevels(
    mf,   // data to be filled
    nudge_time,
    cmf,  // new data at the previous level
    fmf,  // old data at the current level
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
    bfact_coarse,bfact_fine,
    desc_grid_type,
    debug_fillpatch);  
     
  } else 
   amrex::Error("level invalid");

  DComp += ncomp_range;

 } // igroup=0..ranges.size()-1

}   // FillPatch


void
AmrLevel::FillCoarsePatchGHOST (
                           MultiFab& cmf, // source (scomp..scomp+ncomp-1)
                           MultiFab& mf,  // dest (scomp..scomp+ncomp-1)
                           Real      time,
                           int       index,
                           int       scomp, //cmf: scomp..scomp+ncomp-1
                           Vector<int> scompBC_map, // 0..ncomp-1
                           int       ncomp,
			   int debug_fillpatch)
{
 BL_PROFILE("AmrLevel::FillCoarsePatchGHOST()");

 if (level<=0)
  amrex::Error("level invalid in FillCoarsePatchGHOST");

 if (ncomp+scomp>mf.nComp())
  amrex::Error("ncomp+scomp>mf.nComp()");
 if (ncomp+scomp>cmf.nComp())
  amrex::Error("ncomp+scomp>cmf.nComp()");

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
 MultiFab* cmf_part=new MultiFab(cmf_BA,cdm,ncomp,0,
   MFInfo().SetTag("cmf_part"),FArrayBoxFactory());
 MultiFab::Copy(*cmf_part,cmf,scomp,0,ncomp,0);

 int                     DComp   = scomp;
 const StateDescriptor&  descGHOST = desc_lstGHOST[index];
 IndexType desc_typ(descGHOST.getType());
 int desc_grid_type=-1;
 StateData::get_grid_type(desc_typ,desc_grid_type);

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

 for (unsigned int igroup = 0; igroup < ranges.size(); igroup++) {

  const int     scomp_range  = ranges[igroup].first;
  const int     ncomp_range  = ranges[igroup].second;
  if (igroup==0) {
   if (scomp_range!=0)
    amrex::Error("scomp_range!=0");
  }
  Interpolater* mapper = descGHOST.interp(scompBC_map[scomp_range]);

  BoxArray crseBA(mf_BA.size());
  for (int j = 0, N_CBA = crseBA.size(); j < N_CBA; ++j) {
   BL_ASSERT(mf_BA[j].ixType() == descGHOST.getType());
   const Box& bx = mf_BA[j];
   crseBA.set(j,mapper->CoarseBox(bx,bfact_coarse,bfact_fine,desc_grid_type));
  }

   // ghost cells do not have to be initialized.
   // call InterpBordersGHOST after FillCoarsePatchGHOST.
  MultiFab crseMF(crseBA,dm,ncomp_range,0,
	MFInfo().SetTag("crseMF"),FArrayBoxFactory());

  int scomp_data=scomp_range;
  int dcomp_data=scomp+scomp_data;

  if (dcomp_data!=DComp)
   amrex::Error("dcomp_data!=DComp");

  Vector<int> local_scompBC_map;
  local_scompBC_map.resize(ncomp_range);
  for (int isub=0;isub<ncomp_range;isub++)
   local_scompBC_map[isub]=scompBC_map[scomp_range+isub];

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
    bfact_coarse, 
    debug_fillpatch);

  Vector< BCRec > local_bcs;
  const Vector< BCRec> & global_bcs=descGHOST.getBCs();
  local_bcs.resize(ncomp_range);
  for (int isub=0;isub<ncomp_range;isub++)
   local_bcs[isub]=global_bcs[local_scompBC_map[isub]]; 

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(mf.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(mf,false); mfi.isValid(); ++mfi) {

   const Box& tilegrid=mfi.tilebox();
   const Box& dbx = mfi.validbox();

   int tid_current=0;
#ifdef _OPENMP
   tid_current = omp_get_thread_num();
#endif
   if ((tid_current>=0)&&(tid_current<thread_class::nthreads)) {
    // do nothing
   } else
    amrex::Error("tid_current invalid");

   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

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
		  bfact_coarse,bfact_fine,
                  desc_grid_type);
  }  // mfi
} // omp
  thread_class::sync_tile_d_numPts();
  ParallelDescriptor::ReduceRealSum(thread_class::tile_d_numPts[0]);
  thread_class::reconcile_d_numPts(20);

  DComp += ncomp_range;

 } // igroup=0..ranges.size()-1

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
                     int       ncomp,
		     int debug_fillpatch)
{
 BL_PROFILE("AmrLevel::InterpBordersGHOST()");

 if (level<0)
  amrex::Error("level invalid in InterpBordersGHOST");

 BL_ASSERT(ncomp <= (mf.nComp()-scomp));

 BL_ASSERT(0 <= index && index < desc_lstGHOST.size());

 if (scompBC_map.size()!=ncomp)
  amrex::Error("scompBC_map has invalid size");

 int ngrow=mf.nGrow();
 const BoxArray& mf_BA = mf.boxArray();

 if (ngrow<=0)
  amrex::Error("ngrow<=0 in InterpBordersGHOST");

 DistributionMapping dm=mf.DistributionMap();

 MultiFab fmf(mf_BA,dm,ncomp,0,
   MFInfo().SetTag("fmf"),FArrayBoxFactory());

  // dstmf,srcmf,srccomp,dstcomp,ncomp,ngrow
 MultiFab::Copy(fmf,mf,scomp,0,ncomp,0);

 MultiFab* cmf_part;

 if (level>0) {
  const BoxArray& cmf_BA=cmf.boxArray();
  DistributionMapping cdm=cmf.DistributionMap();
  cmf_part=new MultiFab(cmf_BA,cdm,ncomp,0,
   MFInfo().SetTag("cmf_part"),FArrayBoxFactory());
  MultiFab::Copy(*cmf_part,cmf,scomp,0,ncomp,0);
 }  // level>0
 
 int                     DComp   = scomp;
 const StateDescriptor&  descGHOST = desc_lstGHOST[index];
 IndexType desc_typ(descGHOST.getType());
 int desc_grid_type=-1;
 StateData::get_grid_type(desc_typ,desc_grid_type);

 int bfact_fine=parent->Space_blockingFactor(level);

 descGHOST.check_inRange(scompBC_map, ncomp);
 std::vector< std::pair<int,int> > ranges = 
   descGHOST.sameInterps(scompBC_map,ncomp);

 for (unsigned int igroup = 0; igroup < ranges.size(); igroup++) {
  const int     scomp_range  = ranges[igroup].first;
  const int     ncomp_range  = ranges[igroup].second;
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
  for (int isub=0;isub<ncomp_range;isub++)
   local_scompBC_map[isub]=scompBC_map[scomp_range+isub];

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
    bfact_fine,
    debug_fillpatch);
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
    bfact_coarse,bfact_fine,
    desc_grid_type,
    debug_fillpatch);  
     
  } else 
   amrex::Error("level invalid");

  DComp += ncomp_range;

 } // igroup=0..ranges.size()-1

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
                     int       ncomp,
		     int debug_fillpatch)
{
 BL_PROFILE("AmrLevel::InterpBorders()");

 if (level<0)
  amrex::Error("level invalid in InterpBorders");

 BL_ASSERT(ncomp <= (mf.nComp()-scomp));

 BL_ASSERT(0 <= index && index < desc_lst.size());

 if (scompBC_map.size()!=ncomp)
  amrex::Error("scompBC_map has invalid size");

 int ngrow=mf.nGrow();
 const BoxArray& mf_BA = mf.boxArray();

 if (ngrow<=0)
  amrex::Error("ngrow<=0 in InterpBorders");

 DistributionMapping dm=mf.DistributionMap();

 MultiFab fmf(mf_BA,dm,ncomp,0,
   MFInfo().SetTag("fmf"),FArrayBoxFactory());

  // dstmf,srcmf,srccomp,dstcomp,ncomp,ngrow
 MultiFab::Copy(fmf,mf,scomp,0,ncomp,0);

 MultiFab* cmf_part;

 if (level>0) {
  const BoxArray& cmf_BA=cmf.boxArray();
  DistributionMapping cdm=cmf.DistributionMap();
  cmf_part=new MultiFab(cmf_BA,cdm,ncomp,0,
    MFInfo().SetTag("cmf_part"),FArrayBoxFactory());
  MultiFab::Copy(*cmf_part,cmf,scomp,0,ncomp,0);
 }  // level>0
 
 int                     DComp   = scomp;
 const StateDescriptor&  desc    = desc_lst[index];
 IndexType desc_typ(desc.getType());
 int desc_grid_type=-1;
 StateData::get_grid_type(desc_typ,desc_grid_type);

 int bfact_fine=parent->Space_blockingFactor(level);

 desc.check_inRange(scompBC_map, ncomp);
 std::vector< std::pair<int,int> > ranges = 
   desc.sameInterps(scompBC_map,ncomp);

 for (unsigned int igroup = 0; igroup < ranges.size(); igroup++) {
  const int     scomp_range  = ranges[igroup].first;
  const int     ncomp_range  = ranges[igroup].second;
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
  for (int isub=0;isub<ncomp_range;isub++)
   local_scompBC_map[isub]=scompBC_map[scomp_range+isub];

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
    bfact_fine,
    debug_fillpatch);
  } else if (level>0) {

   AmrLevel&               clev    = parent->getLevel(level-1);
   const Geometry&         cgeom   = clev.geom;
   int bfact_coarse=parent->Space_blockingFactor(level-1);
   StateData& cstatedata = clev.state[index];
   StateDataPhysBCFunct cbc(cstatedata,cgeom);

    // declared in: FillPatchUtil.cpp
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
    bfact_coarse,bfact_fine,
    desc_grid_type,
    debug_fillpatch);  
     
  } else 
   amrex::Error("level invalid");

  DComp += ncomp_range;
 } // igroup=0..ranges.size()-1

 if (level>0) 
  delete cmf_part;

}   // InterpBorders


void
AmrLevel::FillCoarsePatch (MultiFab& mf,
                           int       dcomp, // update mf: dcomp..dcomp+ncomp-1
                           Real      time,
                           int       index,
                           int       scomp, // scompBC_map[i]=scomp+i
                           int       ncomp,
			   int debug_fillpatch)
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
 for (int icomp=0;icomp<ncomp;icomp++)
  scompBC_map[icomp]=scomp+icomp;

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
 IndexType desc_typ(desc.getType());
 int desc_grid_type=-1;
 StateData::get_grid_type(desc_typ,desc_grid_type);

 const Box&              pdomain = state[index].getDomain();
 const BoxArray&         mf_BA   = mf.boxArray();
 DistributionMapping dm=mf.DistributionMap();

 int bfact_fine=parent->Space_blockingFactor(level);

 desc.check_inRange(scompBC_map, ncomp);

 std::vector< std::pair<int,int> > ranges = 
   desc.sameInterps(scompBC_map,ncomp);

 for (unsigned int igroup = 0; igroup < ranges.size(); igroup++) {
  const int     scomp_range  = ranges[igroup].first;
  const int     ncomp_range  = ranges[igroup].second;
  Interpolater* mapper = desc.interp(scompBC_map[scomp_range]);

  BoxArray crseBA(mf_BA.size());
   
  for (int j = 0, N = crseBA.size(); j < N; ++j) {
   BL_ASSERT(mf_BA[j].ixType() == desc.getType());
   const Box& bx = mf_BA[j];
   crseBA.set(j,mapper->CoarseBox(bx,bfact_coarse,bfact_fine,desc_grid_type));
  }

    // ngrow=0
  MultiFab crseMF(crseBA,dm,ncomp_range,0,
    MFInfo().SetTag("crseMF"),FArrayBoxFactory());

  Vector<int> local_scompBC_map;
  local_scompBC_map.resize(ncomp_range);
  for (int isub=0;isub<ncomp_range;isub++)
   local_scompBC_map[isub]=scompBC_map[scomp_range+isub];

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
    bfact_coarse,
    debug_fillpatch);

  Vector< BCRec > local_bcs;
  const Vector< BCRec> & global_bcs=desc.getBCs();

  local_bcs.resize(ncomp_range);
  for (int isub=0;isub<ncomp_range;isub++)
   local_bcs[isub]=global_bcs[local_scompBC_map[isub]]; 

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(mf.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(mf,false); mfi.isValid(); ++mfi) {

   const Box& tilegrid=mfi.tilebox();
   const Box& dbx = mfi.validbox();
	    
   int tid_current=0;
#ifdef _OPENMP
   tid_current = omp_get_thread_num();
#endif
   if ((tid_current>=0)&&(tid_current<thread_class::nthreads)) {
    // do nothing
   } else
    amrex::Error("tid_current invalid");

   thread_class::tile_d_numPts[tid_current]+=tilegrid.d_numPts();

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
		  bfact_coarse,bfact_fine,
                  desc_grid_type);
  }  // mfi
} // omp
  thread_class::sync_tile_d_numPts();
  ParallelDescriptor::ReduceRealSum(thread_class::tile_d_numPts[0]);
  thread_class::reconcile_d_numPts(21);

  StateDataPhysBCFunct physbc_fine(state[index],geom);
  physbc_fine.FillBoundary(level,mf,nudge_time,DComp,
    local_scompBC_map,ncomp_range,bfact_fine);

  DComp += ncomp_range;
 } // igroup=0..ranges.size()-1

}   // FillCoarsePatch

Vector<int>
AmrLevel::getBCArray (int State_Type,
                      int gridno,
                      int strt_comp,
                      int ncomp)
{
    Vector<int> bc(2*AMREX_SPACEDIM*ncomp);

    for (int n = 0; n < ncomp; n++)
    {
        const int* b_rec = state[State_Type].getBC(strt_comp+n,gridno).vect();
        for (int m = 0; m < 2*AMREX_SPACEDIM; m++)
            bc[2*AMREX_SPACEDIM*n + m] = b_rec[m];
    }

    return bc;
} // getBCArray

Vector<int>
AmrLevel::getBCArrayGHOST (int State_Type,
                      int gridno,
                      int strt_comp,
                      int ncomp)
{
 Vector<int> bc(2*AMREX_SPACEDIM*ncomp);

 for (int n = 0; n < ncomp; n++) {
  const int* b_rec = state[State_Type].getBCGHOST(strt_comp+n,gridno).vect();
  for (int m = 0; m < 2*AMREX_SPACEDIM; m++)
   bc[2*AMREX_SPACEDIM*n + m] = b_rec[m];
 } // n

 return bc;
} // getBCArrayGHOST


} // namespace amrex

