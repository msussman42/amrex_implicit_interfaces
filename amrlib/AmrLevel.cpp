
#include <sstream>

#include <unistd.h>
#include <memory>
#include <limits>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_ParmParse.H>
#include <AmrLevel.H>
#include <FillPatchUtil.H>

#include <EXTRAP_COMP.H>

namespace amrex {

DescriptorList AmrLevel::desc_lst;
DescriptorList AmrLevel::desc_lstGHOST;

std::ofstream FSI_container_class::CTML_checkpoint_file;

void FSI_container_class::open_checkpoint(const std::string& FullPath) {

// use std::ios::in for restarting  (std::ofstream::in ok too?)
 if (ParallelDescriptor::IOProcessor()) {
  std::string CTML_FullPathName  = FullPath+"/CTML";
  CTML_checkpoint_file.open(CTML_FullPathName.c_str(),
    std::ios::out|std::ios::trunc|std::ios::binary);
  if (!CTML_checkpoint_file.good())
   amrex::FileOpenFailed(CTML_FullPathName);
  int old_prec=CTML_checkpoint_file.precision(15);
 }

} //end subroutine open_checkpoint

void FSI_container_class::close_checkpoint() {

 if (ParallelDescriptor::IOProcessor()) {
  CTML_checkpoint_file.close();
 }

}

void FSI_container_class::checkpoint(int check_id) {

 if (ParallelDescriptor::IOProcessor()) {
  CTML_checkpoint_file << check_id << '\n';

  CTML_checkpoint_file << CTML_num_solids << '\n';
  CTML_checkpoint_file << FSI_num_scalars << '\n';
  for (int dir=0;dir<3;dir++) {
   CTML_checkpoint_file << max_num_nodes[dir] << '\n';
  }
  CTML_checkpoint_file << max_num_elements << '\n';
  CTML_checkpoint_file << structured_flag << '\n';
  CTML_checkpoint_file << structure_dim << '\n';
  CTML_checkpoint_file << structure_topology << '\n';
  CTML_checkpoint_file << ngrow_node << '\n';
  CTML_checkpoint_file << node_list.size() << '\n';
  for (int i=0;i<node_list.size();i++) {
   CTML_checkpoint_file << node_list[i] << '\n';
   CTML_checkpoint_file << prev_node_list[i] << '\n';
   CTML_checkpoint_file << velocity_list[i] << '\n';
   CTML_checkpoint_file << prev_velocity_list[i] << '\n';
   CTML_checkpoint_file << init_node_list[i] << '\n';
  }
  CTML_checkpoint_file << element_list.size() << '\n';
  for (int i=0;i<element_list.size();i++) {
   CTML_checkpoint_file << element_list[i] << '\n';
  }
  CTML_checkpoint_file << mass_list.size() << '\n';
  for (int i=0;i<mass_list.size();i++) {
   CTML_checkpoint_file << mass_list[i] << '\n';
   CTML_checkpoint_file << temp_list[i] << '\n';
  }
  CTML_checkpoint_file << scalar_list.size() << '\n';
  for (int i=0;i<scalar_list.size();i++) {
   CTML_checkpoint_file << scalar_list[i] << '\n';
   CTML_checkpoint_file << prev_scalar_list[i] << '\n';
  }
 }

} //end subroutine checkpoint


void FSI_container_class::restart(int check_id,std::istream& is) {

 int local_check_id;
 is >> local_check_id;
 if (local_check_id==check_id) {
  //do nothing
 } else
  amrex::Error("local_check_id invalid");

 int local_node_list_size;
 int local_element_list_size;
 int local_mass_list_size;
 int local_scalar_list_size;

 is >> CTML_num_solids;
 is >> FSI_num_scalars;
 for (int dir=0;dir<3;dir++) {
  is >> max_num_nodes[dir];
 }
 is >> max_num_elements;
 is >> structured_flag;
 is >> structure_dim;
 is >> structure_topology;
 is >> ngrow_node;

 is >> local_node_list_size;

 node_list.resize(local_node_list_size);
 init_node_list.resize(local_node_list_size);
 prev_node_list.resize(local_node_list_size);
 velocity_list.resize(local_node_list_size);
 prev_velocity_list.resize(local_node_list_size);

 for (int i=0;i<node_list.size();i++) {
  is >> node_list[i];
  is >> prev_node_list[i];
  is >> velocity_list[i];
  is >> prev_velocity_list[i];
  is >> init_node_list[i];
 }

 is >> local_element_list_size;

 element_list.resize(local_element_list_size);
 for (int i=0;i<element_list.size();i++) {
  is >> element_list[i];
 }

 is >> local_mass_list_size;

 mass_list.resize(local_mass_list_size);
 temp_list.resize(local_mass_list_size);

 for (int i=0;i<mass_list.size();i++) {
  is >> mass_list[i];
  is >> temp_list[i];
 }

 is >> local_scalar_list_size;
 scalar_list.resize(local_scalar_list_size);
 prev_scalar_list.resize(local_scalar_list_size);
 for (int i=0;i<scalar_list.size();i++) {
  is >> scalar_list[i];
  is >> prev_scalar_list[i];
 }

} //end subroutine restart


void FSI_container_class::initData_FSI(
  const int CTML_num_solids_init,
  const int max_num_nodes_init[3],
  const int max_num_elements_init,
  const int FSI_num_scalars_init,
  const int structured_flag_init, 
  const int structure_dim_init,
  const int structure_topology_init,
  const int ngrow_node_init) {

if ((CTML_num_solids_init>=0)&&
    (ngrow_node_init>=0)&&
    (max_num_nodes_init[0]>=0)&&
    (max_num_nodes_init[1]>=0)&&
    (max_num_nodes_init[2]>=0)&&
    (max_num_elements_init>=0)&&
    (FSI_num_scalars_init>=0) {

 CTML_num_solids=CTML_num_solids_init;
 FSI_num_scalars=FSI_num_scalars_init;
 for (int dir=0;dir<3;dir++) {
  max_num_nodes[dir]=max_num_nodes_init[dir];
 }
 max_num_elements=max_num_elements_init;
 structured_flag=structured_flag_init;
 structure_dim=structure_dim_init;
 structure_topology=structure_topology_init;
 ngrow_node=ngrow_node_init;

 int max_num_nodes_grow=max_num_nodes[0];

 if (structured_flag==0) {
  //do nothing
 } else if (structured_flag==1) {
  if (structure_topology==0) { //filament

   if (AMREX_SPACEDIM==2) {
    //do nothing
   } else
    amrex::Error("filament only allowed in 2D");

   if (max_num_nodes[0]==0) {
    max_num_nodes_grow=0;
   } else if (max_num_nodes[0]>0) {
     max_num_nodes_grow=max_num_nodes[0]+2*ngrow_node;
   } else
    amrex::Error("max_num_nodes[0] invalid");

  } else if (structure_topology==1) { //sheet

   if (AMREX_SPACEDIM==3) {
    //do nothing
   } else
    amrex::Error("sheet only allowed in 3D");

   if (max_num_nodes[0]==0) {
    max_num_nodes_grow=0;
   } else if (max_num_nodes[0]>0) {
    max_num_nodes_grow=
     (max_num_nodes[0]+2*ngrow_node)*
     (max_num_nodes[1]+2*ngrow_node);
   } else
    amrex::Error("max_num_nodes[0] invalid");
  } else if (structure_topology==2) { //volumetric
   if (max_num_nodes[0]==0) {
    max_num_nodes_grow=0;
   } else if (max_num_nodes[0]>0) {
    if (AMREX_SPACEDIM==2) {
     max_num_nodes_grow=
      (max_num_nodes[0]+2*ngrow_node)*
      (max_num_nodes[1]+2*ngrow_node);
    } else if (AMREX_SPACEDIM==3) {
     max_num_nodes_grow=
      (max_num_nodes[0]+2*ngrow_node)*
      (max_num_nodes[1]+2*ngrow_node)*
      (max_num_nodes[2]+2*ngrow_node);
    } else
     amrex::Error("AMREX_SPACEDIM invalid");
   } else
    amrex::Error("max_num_nodes[0] invalid");
  } else
   amrex::Error("structure_topology invalid");
 } else
  amrex::Error("structured_flag invalid");

 prev_node_list.resize(CTML_num_solids*max_num_nodes_grow*3);
 node_list.resize(CTML_num_solids*max_num_nodes_grow*3);
 prev_velocity_list.resize(CTML_num_solids*max_num_nodes_grow*3);
 velocity_list.resize(CTML_num_solids*max_num_nodes_grow*3);
 element_list.resize(CTML_num_solids*max_num_elements*4);
 init_node_list.resize(CTML_num_solids*max_num_nodes_grow*3);
 mass_list.resize(CTML_num_solids*max_num_nodes_grow);
 temp_list.resize(CTML_num_solids*max_num_nodes_grow);
 scalar_list.resize(CTML_num_solids*max_num_nodes_grow*FSI_num_scalars);
 prev_scalar_list.resize(CTML_num_solids*max_num_nodes_grow*FSI_num_scalars);
} else {
 amrex::Error("num_solids/nodes/scalars or elements invalid");
}

} // end subroutine initData_FSI()

void FSI_container_class::FSI_flatten(Vector< Real >& flattened_data) {

 int node_list_size=node_list.size();

 if (node_list_size==velocity_list.size()) {
  //do nothing
 } else
  amrex::Error("expecting node_list_size==velocity_list.size()");

 int element_list_size=element_list.size();

 if (node_list_size==init_node_list.size()) {
  //do nothing
 } else
  amrex::Error("expecting node_list_size==init_node_list.size()");

 int mass_list_size=mass_list.size();

 if (mass_list_size==temp_list.size()) {
  //do nothing
 } else
  amrex::Error("expecting mass_list_size==temp_list.size()");

 int scalar_list_size=scalar_list.size();

 if (scalar_list_size==prev_scalar_list.size()) {
  //do nothing
 } else
  amrex::Error("expecting scalar_list_size==prev_scalar_list.size()");

 flattened_data.resize(FSIcontain_size);
 flattened_data[FSIcontain_num_solids]=CTML_num_solids; 
 flattened_data[FSIcontain_num_scalars]=FSI_num_scalars; 
 for (int dir=0;dir<3;dir++) {
  flattened_data[FSIcontain_max_num_nodes+dir]=max_num_nodes[dir]; 
 }
 flattened_data[FSIcontain_max_num_elements]=max_num_elements; 
 flattened_data[FSIcontain_structured_flag]=structured_flag; 
 flattened_data[FSIcontain_structure_dim]=structure_dim; 
 flattened_data[FSIcontain_structure_topology]=structure_topology; 
 flattened_data[FSIcontain_ngrow_node]=ngrow_node; 

 for (int i=0;i<node_list_size;i++) {
  flattened_data[FSIcontain_node_list+i]=node_list[i];
  flattened_data[FSIcontain_prev_node_list+i]=prev_node_list[i];
  flattened_data[FSIcontain_velocity_list+i]=velocity_list[i];
  flattened_data[FSIcontain_prev_velocity_list+i]=prev_velocity_list[i];
  flattened_data[FSIcontain_init_node_list+i]=init_node_list[i];
 } 
 for (int i=0;i<element_list_size;i++) {
  flattened_data[FSIcontain_element_list+i]=element_list[i];
 } 
 for (int i=0;i<mass_list_size;i++) {
  flattened_data[FSIcontain_mass_list+i]=mass_list[i];
  flattened_data[FSIcontain_temp_list+i]=temp_list[i];
 } 
 for (int i=0;i<scalar_list_size;i++) {
  flattened_data[FSIcontain_scalar_list+i]=scalar_list[i];
  flattened_data[FSIcontain_prev_scalar_list+i]=prev_scalar_list[i];
 } 

} //end subroutine FSI_flatten

void FSI_container_class::FSI_unflatten(Vector< Real > flattened_data) {

 CTML_num_solids=(int) flattened_data[FSIcontain_num_solids]; 
 FSI_num_scalars=(int) flattened_data[FSIcontain_num_scalars]; 
 for (int dir=0;dir<3;dir++) {
  max_num_nodes[dir]=(int) flattened_data[FSIcontain_max_num_nodes+dir];
 }
 max_num_elements=(int) flattened_data[FSIcontain_max_num_elements];

 structured_flag=(int) flattened_data[FSIcontain_structured_flag];
 structure_dim=(int) flattened_data[FSIcontain_structure_dim];
 structure_topology=(int) flattened_data[FSIcontain_structure_topology];
 ngrow_node=(int) flattened_data[FSIcontain_ngrow_node];

 initData_FSI(
  CTML_num_solids,
  max_num_nodes,
  max_num_elements,
  FSI_num_scalars,
  structured_flag, 
  structure_dim,
  structure_topology,
  ngrow_node);

 int node_list_size=node_list.size();
 int element_list_size=element_list.size();
 int mass_list_size=mass_list.size();
 int scalar_list_size=scalar_list.size();

 if (node_list_size==velocity_list.size()) {
  //do nothing
 } else
  amrex::Error("expecting node_list_size==velocity_list.size()");

 if (node_list_size==init_node_list.size()) {
  //do nothing
 } else
  amrex::Error("expecting node_list_size==init_node_list.size()");

 if (mass_list_size==temp_list.size()) {
  //do nothing
 } else
  amrex::Error("expecting mass_list_size==temp_list.size()");

 if (scalar_list_size==prev_scalar_list.size()) {
  //do nothing
 } else
  amrex::Error("expecting scalar_list_size==prev_scalar_list.size()");

 if (flattened_data.size()==FSIcontain_size) {
  //do nothing
 } else
  amrex::Error("flattened_data.size()==FSIcontain_size failed");

 for (int i=0;i<node_list_size;i++) {
  node_list[i]=flattened_data[FSIcontain_node_list+i];
  prev_node_list[i]=flattened_data[FSIcontain_prev_node_list+i];
  velocity_list[i]=flattened_data[FSIcontain_velocity_list+i];
  prev_velocity_list[i]=flattened_data[FSIcontain_prev_velocity_list+i];
  init_node_list[i]=flattened_data[FSIcontain_init_node_list+i];
 } 
 for (int i=0;i<element_list_size;i++) {
  element_list[i]=(int) flattened_data[FSIcontain_element_list+i];
 } 
 for (int i=0;i<mass_list_size;i++) {
  mass_list[i]=flattened_data[FSIcontain_mass_list+i];
  temp_list[i]=flattened_data[FSIcontain_temp_list+i];
 } 
 for (int i=0;i<scalar_list_size;i++) {
  scalar_list[i]=flattened_data[FSIcontain_scalar_list+i];
  prev_scalar_list[i]=flattened_data[FSIcontain_prev_scalar_list+i];
 } 

} //end subroutine FSI_unflatten


void FSI_container_class::copyFrom_FSI(const FSI_container_class& source_FSI) {

 initData_FSI(
   source_FSI.CTML_num_solids,
   source_FSI.max_num_nodes,
   source_FSI.max_num_elements,
   source_FSI.FSI_num_scalars,
   source_FSI.structured_flag,
   source_FSI.structure_dim,
   source_FSI.structure_topology,
   source_FSI.ngrow_node);

 for (int ielem=0;ielem<source_FSI.element_list.size();ielem++) {
  element_list[ielem]=source_FSI.element_list[ielem];
 } 

 if (source_FSI.node_list.size()==
     source_FSI.prev_node_list.size()) {
  //do nothing
 } else
  amrex::Error("node_list or prev_node_list size invalid");

 if (source_FSI.node_list.size()==
     source_FSI.init_node_list.size()) {
  //do nothing
 } else
  amrex::Error("node_list or init_node_list size invalid");

 if (source_FSI.node_list.size()==
     source_FSI.velocity_list.size()) {
  //do nothing
 } else
  amrex::Error("node_list or velocity_list size invalid");

 if (source_FSI.node_list.size()==
     source_FSI.prev_velocity_list.size()) {
  //do nothing
 } else
  amrex::Error("node_list or prev_velocity_list size invalid");

 for (int inode=0;inode<source_FSI.node_list.size();inode++) {
  node_list[inode]=source_FSI.node_list[inode];
  prev_node_list[inode]=source_FSI.prev_node_list[inode];
  init_node_list[inode]=source_FSI.init_node_list[inode];
  velocity_list[inode]=source_FSI.velocity_list[inode];
  prev_velocity_list[inode]=source_FSI.prev_velocity_list[inode];
 } 

 if (source_FSI.mass_list.size()==
     source_FSI.temp_list.size()) {
  //do nothing
 } else
  amrex::Error("mass_list or temp_list size invalid");

 for (int inode=0;inode<source_FSI.mass_list.size();inode++) {
  mass_list[inode]=source_FSI.mass_list[inode];
  temp_list[inode]=source_FSI.temp_list[inode];
 }

 if (source_FSI.scalar_list.size()==
     source_FSI.prev_scalar_list.size()) {
  //do nothing
 } else
  amrex::Error("scalar_list or prev_scalar_list size invalid");

 for (int inode=0;inode<source_FSI.scalar_list.size();inode++) {
  scalar_list[inode]=source_FSI.scalar_list[inode];
  prev_scalar_list[inode]=source_FSI.prev_scalar_list[inode];
 }

} // end subroutine copyFrom_FSI

void FSI_container_class::clear_FSI() {

 int local_num_solids=0;
 int local_num_scalars=0;

 int local_num_nodes[3];
 for (int dir=0;dir<3;dir++) {
  local_num_nodes[dir]=0;
 }
 int local_num_elements=0;

 initData_FSI(local_num_solids,local_num_nodes,
	local_num_elements,local_num_scalars);

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

    parent = &papa;

    if (grids.size()<1) {
     std::cout << "lev= " << lev << " time= " << time << '\n';
     amrex::Error("AmrLevel: grids.size<1");
    }

    level  = lev;

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
    int local_nmat=parent->global_AMR_num_materials;

    if (level==0) {

     new_data_FSI.resize(level_MAX_NUM_SLAB);

     for (int i=0;i<=time_order;i++) {

      new_data_FSI[i].clear_FSI();

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

 if ((old_finest_level>=0)&&
     (new_finest_level>=0)&&
     (old_finest_level<=max_level)&&
     (new_finest_level<=max_level)) {
  //do nothing
 } else
  amrex::Error("old_finest_level or new_finest_level invalid");

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

 //SUSSMAN: load CTML FSI checkpoint data
 if (level==0) {

  int time_order=parent->Time_blockingFactor();
  int local_nmat=parent->global_AMR_num_materials;

  new_data_FSI.resize(level_MAX_NUM_SLAB);

  std::string CTML_FullPathName  = FullPath+"/CTML";

  Vector<char> fileCharPtr;
   //we assume that all of the CTML data fits on each "core"
  ParallelDescriptor::ReadAndBcastFile(CTML_FullPathName, fileCharPtr);
  std::string fileCharPtrString(fileCharPtr.dataPtr());
  std::istringstream CTML_is(fileCharPtrString, std::istringstream::in);

  for (int i=0;i<=time_order;i++) {

   new_data_FSI[i].restart(i,CTML_is);

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

      //SUSSMAN: output CTML FSI checkpoint data
    if (level==0) {
     int time_order=parent->Time_blockingFactor();
     new_data_FSI[0].open_checkpoint(FullPath);
     for (int i=0;i<=time_order;i++) {
      new_data_FSI[i].checkpoint(i);
     }
     new_data_FSI[0].close_checkpoint();

     ParallelDescriptor::Barrier("AmrLevel::checkPoint");
    }
} // end subroutine AmrLevel::checkPoint

AmrLevel::~AmrLevel ()
{

    if (level==0) {
     int time_order=parent->Time_blockingFactor();
     int local_nmat=parent->global_AMR_num_materials;

     for (int i=0;i<=time_order;i++) {

      new_data_FSI[i].clear_FSI();

     } // for (int i=0;i<=time_order;i++) 
     new_data_FSI.resize(0);
    } else if (level>0) {
     // do nothing
    } else
     amrex::Error("level invalid");

    parent = 0;
}

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
  thread_class::reconcile_d_numPts(LOOP_MAPPER_INTERP_GHOST,
      "FillCoarsePatchGHOST");

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

 BL_ASSERT(ncomp<=(mf.nComp()-scomp));

 BL_ASSERT((0<=index)&&(index<desc_lstGHOST.size()));

 if (scompBC_map.size()!=ncomp)
  amrex::Error("scompBC_map has invalid size");

 int ngrow=mf.nGrow();
 const BoxArray& mf_BA = mf.boxArray();

 if (ngrow<0)
  amrex::Error("ngrow<0 in InterpBordersGHOST");

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
 
 int DComp = scomp;
 const StateDescriptor& descGHOST = desc_lstGHOST[index];
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

 if (ngrow<0)
  amrex::Error("ngrow<0 in InterpBorders");

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

 int                     DComp   = dcomp;
 const StateDescriptor&  desc    = desc_lst[index];
 IndexType desc_typ(desc.getType());
 int desc_grid_type=-1;
 StateData::get_grid_type(desc_typ,desc_grid_type);

 int ngrow=mf.nGrow();

 if (ngrow==0) {
  if ((desc_grid_type==-1)||
      ((desc_grid_type>=0)&&
       (desc_grid_type<AMREX_SPACEDIM))) {
   //do nothing
  } else {
   amrex::Error("desc_grid_type invalid (ngrow=0)");
  }
 } else if (ngrow==1) {
  if (desc_grid_type==-1) {
   //do nothing
  } else {
   amrex::Error("desc_grid_type invalid (ngrow=1)");
  }
 } else
  amrex::Error("expecting ngrow=0 or 1 in AmrLevel::FillCoarsePatch");

 const Box& domain = geom.Domain();
 const int* domlo = domain.loVect();
 const int* domhi = domain.hiVect();

 AmrLevel&               clev    = parent->getLevel(level-1);
 const Geometry&         cgeom   = clev.geom;
 int bfact_coarse=parent->Space_blockingFactor(level-1);
 StateData& cstatedata = clev.state[index];

 Real nudge_time;
 int best_index;
 cstatedata.get_time_index(time,nudge_time,best_index);

 MultiFab& cmf=cstatedata.newData(best_index);
 StateDataPhysBCFunct physbc_coarse(cstatedata,cgeom);

  //pdomain will be different depending on whether the state variable
  //is cell centered or staggared.
 const Box& pdomain = state[index].getDomain();
 const int* pdomlo = pdomain.loVect();
 const int* pdomhi = pdomain.hiVect();

 for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
  if (pdomlo[dir]==domlo[dir]) {
   // do nothing
  } else
   amrex::Error("pdomlo<>domlo");
  if (dir==desc_grid_type) {
   if (pdomhi[dir]==domhi[dir]+1) {
    // do nothing
   } else
    amrex::Error("pdomhi<>domhi+1");
  } else if (dir!=desc_grid_type) {
   if (pdomhi[dir]==domhi[dir]) {
    // do nothing
   } else
    amrex::Error("pdomhi<>domhi");
  } else
   amrex::Error("dir,desc_grid_type breakdown");
 } //dir=0..sdim-1

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

   Box grow_bx(bx);
   const int* bx_lo=bx.loVect();
   const int* bx_hi=bx.hiVect();

   if (ngrow==1) {
    for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
     if (bx_lo[dir]>domlo[dir]) 
      grow_bx.growLo(dir,ngrow);
     if (bx_hi[dir]<domhi[dir]) 
      grow_bx.growHi(dir,ngrow);
    } //dir=0..sdim-1
    Box grow_bx_test=grow_bx & domain;
    if (grow_bx_test==grow_bx) {
     //do nothing
    } else
     amrex::Error("grow_bx_test!=grow_bx");
   } else if (ngrow==0) {
    //do nothing
   } else {
    amrex::Error("ngrow invalid");
   }

   crseBA.set(j,mapper->CoarseBox(grow_bx,bfact_coarse,bfact_fine,
     desc_grid_type));
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

   const Box& mf_local_box=mf[mfi].box();

   Box dbx_test=dbx;
   for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
    dbx_test.grow(dir,ngrow);
   }

   if (dbx_test==mf_local_box) {
    //do nothing
   } else
    amrex::Error("dbx_test==mf_local_box failed");

   Box grow_bx(dbx);
   const int* bx_lo=dbx.loVect();
   const int* bx_hi=dbx.hiVect();

   if (ngrow==1) {
    for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
     if (bx_lo[dir]>domlo[dir]) 
      grow_bx.growLo(dir,ngrow);
     if (bx_hi[dir]<domhi[dir]) 
      grow_bx.growHi(dir,ngrow);
    } //dir=0..sdim-1
    Box grow_bx_test=grow_bx & domain;
    if (grow_bx_test==grow_bx) {
     //do nothing
    } else
     amrex::Error("grow_bx_test!=grow_bx");
   } else if (ngrow==0) {
    //do nothing
   } else {
    amrex::Error("ngrow invalid");
   }

   if (1==0) {
    std::cout << "ngrow= " << ngrow << '\n';
    std::cout << "dbx= " << dbx << '\n';
    std::cout << "grow_bx= " << grow_bx << '\n';
    std::cout << "DComp= " << DComp << '\n';
    std::cout << "ncomp_range= " << ncomp_range << '\n';
    std::cout << "mfi.index()= " << mfi.index() << '\n';
   }
   mf[mfi].setVal(1.0e+20,grow_bx,DComp,ncomp_range);

   mapper->interp(nudge_time,
                  crseMF[mfi],
                  0,  // crse_comp
		  mf[mfi],
		  DComp,
		  ncomp_range,
		  grow_bx,
		  cgeom,
		  geom,
		  bcr,
                  level-1,level,
		  bfact_coarse,bfact_fine,
                  desc_grid_type);

    //p=0 (infinity norm)
   Real test_norm=mf[mfi].norm(grow_bx,0,DComp,ncomp_range);
   if (test_norm<1.0e+19) {
    //do nothing
   } else {
    amrex::Error("test_norm<1.0e+19 failed");
   }

  }  // mfi
} // omp
  thread_class::sync_tile_d_numPts();
  ParallelDescriptor::ReduceRealSum(thread_class::tile_d_numPts[0]);
  thread_class::reconcile_d_numPts(LOOP_FILLCOARSEPATCH,"FillCoarsePatch");

  mf.FillBoundary(DComp,ncomp_range,geom.periodicity());

  StateDataPhysBCFunct physbc_fine(state[index],geom);
  physbc_fine.FillBoundary(level,mf,nudge_time,DComp,
    local_scompBC_map,ncomp_range,bfact_fine);

  DComp += ncomp_range;
 } // igroup=0..ranges.size()-1

 if (DComp==dcomp+ncomp) {
  //do nothing
 } else
  amrex::Error("expecting DComp==dcomp+ncomp");

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

