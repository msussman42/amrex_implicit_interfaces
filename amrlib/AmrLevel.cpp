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

dynamic_blobclass_array AmrLevel::blob_history_class;

std::ofstream FSI_container_class::CTML_checkpoint_file;

std::ofstream dynamic_blobclass_array::blob_checkpoint_file;

void snapshot_blobclass::sort_axis() {

 for (int outer=0;outer<AMREX_SPACEDIM-1;outer++) {
  for (int inner=AMREX_SPACEDIM-outer-1;inner>=1;inner--) {
   if (blob_axis_len[inner]>blob_axis_len[inner-1]) {
    Real data_save=blob_axis_len[inner-1];
    blob_axis_len[inner-1]=blob_axis_len[inner];
    blob_axis_len[inner]=data_save;
    for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
     int s1=AMREX_SPACEDIM*(inner-1)+dir;
     data_save=blob_axis_evec[s1];
     int s2=AMREX_SPACEDIM*(inner)+dir;
     blob_axis_evec[s1]=blob_axis_evec[s2];  
     blob_axis_evec[s2]=data_save; 
    } //dir=0..sdim-1 
   } //if (blob_axis_len[inner]>blob_axis_len[inner-1]) 
  } //inner
 } //outer=0..sdim-2
} // snapshot_blobclass::sort_axis()

void dynamic_blobclass_array::open_checkpoint(const std::string& FullPath) {

// use std::ios::in for restarting  (std::ofstream::in ok too?)
 if (ParallelDescriptor::IOProcessor()) {
  std::string blob_FullPathName  = FullPath+"/blob";
  blob_checkpoint_file.open(blob_FullPathName.c_str(),
    std::ios::out|std::ios::trunc|std::ios::binary);
  if (!blob_checkpoint_file.good())
   amrex::FileOpenFailed(blob_FullPathName);
  int old_prec=blob_checkpoint_file.precision(15);
 }

} //end subroutine open_checkpoint

void dynamic_blobclass_array::close_checkpoint() {

 if (ParallelDescriptor::IOProcessor()) {
  blob_checkpoint_file.close();
 }

}

void dynamic_blobclass_array::checkpoint(int check_id) {

 if (ParallelDescriptor::IOProcessor()) {
  blob_checkpoint_file << check_id << '\n';

  blob_checkpoint_file << blob_history.size() << '\n';

  blob_checkpoint_file << start_time << '\n';
  blob_checkpoint_file << end_time << '\n';
  blob_checkpoint_file << start_step << '\n';
  blob_checkpoint_file << end_step << '\n';

  for (int i=0;i<blob_history.size();i++) {

   blob_checkpoint_file << blob_history[i].im << '\n';
   blob_checkpoint_file << blob_history[i].start_time << '\n';
   blob_checkpoint_file << blob_history[i].end_time << '\n';
   blob_checkpoint_file << blob_history[i].start_step << '\n';
   blob_checkpoint_file << blob_history[i].end_step << '\n';

   blob_checkpoint_file << blob_history[i].snapshots.size() << '\n';
   for (int j=0;j<blob_history[i].snapshots.size();j++) {
    blob_checkpoint_file << blob_history[i].snapshots[j].blob_volume << '\n';
    blob_checkpoint_file << blob_history[i].snapshots[j].blob_time << '\n';
    blob_checkpoint_file << blob_history[i].snapshots[j].blob_step << '\n';
    for (int k=0;k<AMREX_SPACEDIM;k++) {
     blob_checkpoint_file << 
	 blob_history[i].snapshots[j].blob_center[k] << '\n';
     blob_checkpoint_file << 
	 blob_history[i].snapshots[j].blob_axis_len[k] << '\n';
    } //k=0..sdim-1
    for (int k=0;k<AMREX_SPACEDIM*AMREX_SPACEDIM;k++) {
     blob_checkpoint_file << 
       blob_history[i].snapshots[j].blob_axis_evec[k] << '\n';
    }
   } // j=0;j<snapshots.size()
  } // i=0;i<blob_history_size

 }

} //end subroutine dynamic_blobclass_array::checkpoint

void dynamic_blobclass_array::restart(int check_id,std::istream& is) {

 int local_check_id;
 is >> local_check_id;
 if (local_check_id==check_id) {
  //do nothing
 } else
  amrex::Error("local_check_id invalid");

 int local_blob_history_size;

 is >> local_blob_history_size;

 blob_history.resize(local_blob_history_size);

 is >> start_time;
 is >> end_time;
 is >> start_step;
 is >> end_step;

 for (int i=0;i<blob_history.size();i++) {
  is >> blob_history[i].im;
  is >> blob_history[i].start_time;
  is >> blob_history[i].end_time;
  is >> blob_history[i].start_step;
  is >> blob_history[i].end_step;

  int local_snapshots_size;
  is >> local_snapshots_size;
  blob_history[i].snapshots.resize(local_snapshots_size);
  for (int j=0;j<blob_history[i].snapshots.size();j++) {
   is >> blob_history[i].snapshots[j].blob_volume;
   is >> blob_history[i].snapshots[j].blob_time;
   is >> blob_history[i].snapshots[j].blob_step;

   for (int k=0;k<AMREX_SPACEDIM;k++) {
    is >> blob_history[i].snapshots[j].blob_center[k];
    is >> blob_history[i].snapshots[j].blob_axis_len[k];
   } //k=0..sdim-1
   for (int k=0;k<AMREX_SPACEDIM*AMREX_SPACEDIM;k++) {
    is >> blob_history[i].snapshots[j].blob_axis_evec[k];
   }
  } // j=0;j<snapshots.size()
 } // i=0;i<blob_history_size

} //end subroutine dynamic_blobclass_array::restart


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
    (FSI_num_scalars_init>=0)) {

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

 int max_num_nodes_total_grow=max_num_nodes[0];

 if (structured_flag==0) {
  max_num_nodes_total_grow=max_num_nodes[0];
 } else if (structured_flag==1) {
  if (structure_topology==0) { //filament

   if (AMREX_SPACEDIM==2) {
    //do nothing
   } else
    amrex::Error("filament only allowed in 2D");

   if (max_num_nodes[0]==0) {
    max_num_nodes_total_grow=0;
   } else if (max_num_nodes[0]>0) {
     max_num_nodes_total_grow=max_num_nodes[0]+2*ngrow_node;
   } else
    amrex::Error("max_num_nodes[0] invalid");

  } else if (structure_topology==1) { //sheet

   if (AMREX_SPACEDIM==3) {
    //do nothing
   } else
    amrex::Error("sheet only allowed in 3D");

   if (max_num_nodes[0]==0) {
    max_num_nodes_total_grow=0;
   } else if (max_num_nodes[0]>0) {
    max_num_nodes_total_grow=
     (max_num_nodes[0]+2*ngrow_node)*
     (max_num_nodes[1]+2*ngrow_node);
   } else
    amrex::Error("max_num_nodes[0] invalid");
  } else if (structure_topology==2) { //volumetric
   if (max_num_nodes[0]==0) {
    max_num_nodes_total_grow=0;
   } else if (max_num_nodes[0]>0) {
    if (AMREX_SPACEDIM==2) {
     max_num_nodes_total_grow=
      (max_num_nodes[0]+2*ngrow_node)*
      (max_num_nodes[1]+2*ngrow_node);
    } else if (AMREX_SPACEDIM==3) {
     max_num_nodes_total_grow=
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

 prev_node_list.resize(CTML_num_solids*max_num_nodes_total_grow*3);
 node_list.resize(CTML_num_solids*max_num_nodes_total_grow*3);
 prev_velocity_list.resize(CTML_num_solids*max_num_nodes_total_grow*3);
 velocity_list.resize(CTML_num_solids*max_num_nodes_total_grow*3);
 element_list.resize(CTML_num_solids*max_num_elements*4);
 init_node_list.resize(CTML_num_solids*max_num_nodes_total_grow*3);
 mass_list.resize(CTML_num_solids*max_num_nodes_total_grow);
 temp_list.resize(CTML_num_solids*max_num_nodes_total_grow);
 scalar_list.resize(CTML_num_solids*max_num_nodes_total_grow*FSI_num_scalars);
 prev_scalar_list.resize(
   CTML_num_solids*max_num_nodes_total_grow*FSI_num_scalars);
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
    if (dt>0.0) {
     //do nothing
    } else {
     std::cout << "dt= " << dt << '\n';
     amrex::Error("dt invalid");
    }

    int level_MAX_NUM_SLAB=parent->get_MAX_NUM_SLAB();
    int level_slab_dt_type=parent->get_slab_dt_type();
    if (level_MAX_NUM_SLAB<33+3) //3 to account for LSA
     amrex::Error("level_MAX_NUM_SLAB too small");
    if ((level_slab_dt_type!=0)&&(level_slab_dt_type!=1))
     amrex::Error("level_slab_dt_type invalid");

    int time_order=parent->Time_blockingFactor();
     
     //initialized in: void AmrCore::InitAmr ()
    int local_nmat=parent->global_AMR_num_materials;
    int local_num_species=parent->global_AMR_num_species_var;
    if (local_num_species>=0) {
     //do nothing
    } else
     amrex::Error("local_num_species invalid");

    if (level==0) {

#ifdef AMREX_PARTICLES
     AmrLevel0_new_dataPC.resize(level_MAX_NUM_SLAB);
#endif

      //in the constructor
     new_data_FSI.resize(level_MAX_NUM_SLAB);

     for (int i=0;i<=time_order+parent->LSA_extra_data;i++) {

      new_data_FSI[i].clear_FSI();

#ifdef AMREX_PARTICLES
      using My_ParticleContainer =
        AmrParticleContainer<N_EXTRA_REAL,N_EXTRA_INT,0,0>;
      AmrLevel0_new_dataPC[i] = std::make_unique<My_ParticleContainer>(parent);
#endif

     }// for (int i=0;i<=time_order+parent->LSA_extra_data;i++) 

    } else if ((level>0)&&(level<=max_level)) {
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

 if (level==0) {

#ifdef AMREX_PARTICLES
  AmrLevel0_new_dataPC.resize(level_MAX_NUM_SLAB);
#endif

  //SUSSMAN: load CTML FSI checkpoint data
  new_data_FSI.resize(level_MAX_NUM_SLAB);

  ParmParse ppns("ns");

  int local_num_materials=0;
  ppns.get("num_materials",local_num_materials);
  if (local_num_materials>0) {
   //do nothing
  } else
   amrex::Error("local_num_materials invalid");

  Vector< int > FSI_flag;
  FSI_flag.resize(local_num_materials);

  int query_status=ppns.queryarr("FSI_flag",FSI_flag,0,local_num_materials);

  if (query_status==1) {

   int ctml_count=0;
   for (int i=0;i<FSI_flag.size();i++) {
    if (FSI_flag[i]==FSI_SHOELE_CTML) {
     ctml_count++;
    }
   }

   if (ctml_count>0) {

    int time_order=parent->Time_blockingFactor();
    int local_nmat=parent->global_AMR_num_materials;

    std::string CTML_FullPathName  = FullPath+"/CTML";

    Vector<char> fileCharPtr;
     //we assume that all of the CTML data fits on each "core"
    ParallelDescriptor::ReadAndBcastFile(CTML_FullPathName, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream CTML_is(fileCharPtrString, std::istringstream::in);

    for (int i=0;i<=time_order+parent->LSA_extra_data;i++) {

     new_data_FSI[i].restart(i,CTML_is);

    }//for (int i=0;i<=time_order+parent->LSA_extra_data;i++) 

   } else if (ctml_count==0) {
    //do nothing
   } else
    amrex::Error("ctml_count invalid");

  } else if (query_status==0) {
   //do nothing
  } else
   amrex::Error("query_status invalid");

  std::string blob_FullPathName  = FullPath+"/blob";

  Vector<char> blob_fileCharPtr;
  //we assume that all of the blob history data fits on each "core"
  ParallelDescriptor::ReadAndBcastFile(blob_FullPathName, blob_fileCharPtr);
  std::string blob_fileCharPtrString(blob_fileCharPtr.dataPtr());
  std::istringstream blob_is(blob_fileCharPtrString, std::istringstream::in);

  blob_history_class.restart(0,blob_is);

 } else if ((level>0)&&(level<=max_level)) {
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

    int time_order=parent->Time_blockingFactor();

    if (level==0) {

     for (int i=0;i<=time_order+parent->LSA_extra_data;i++) {

#ifdef AMREX_PARTICLES
      std::string FullPathNamePart=FullPath;
      std::string Part_name="FusionPart";
      std::stringstream i_string_stream(std::stringstream::in |
          std::stringstream::out);
      i_string_stream << i;
      std::string i_string=i_string_stream.str();
      Part_name+=i_string;
      AmrLevel0_new_dataPC[i]->Checkpoint(FullPathNamePart,Part_name);
#endif

     }//for (int i=0;i<=time_order+parent->LSA_extra_data;i++) 

     //SUSSMAN: output CTML FSI checkpoint data
     ParmParse ppns("ns");

     int local_num_materials=0;
     ppns.get("num_materials",local_num_materials);
     if (local_num_materials>0) {
      //do nothing
     } else
      amrex::Error("local_num_materials invalid");

     Vector< int > FSI_flag;
     FSI_flag.resize(local_num_materials);

     int query_status=ppns.queryarr("FSI_flag",FSI_flag,0,local_num_materials);

     if (query_status==1) {

      int ctml_count=0;
      for (int i=0;i<FSI_flag.size();i++) {
       if (FSI_flag[i]==FSI_SHOELE_CTML) {
        ctml_count++;
       }
      }

      if (ctml_count>0) {

       new_data_FSI[0].open_checkpoint(FullPath);
       for (int i=0;i<=time_order+parent->LSA_extra_data;i++) {
        new_data_FSI[i].checkpoint(i);
       }
       new_data_FSI[0].close_checkpoint();

      } else if (ctml_count==0) {
       //do nothing
      } else
       amrex::Error("ctml_count invalid");

     } else if (query_status==0) {
      //do nothing
     } else
      amrex::Error("query_status invalid");

     ParallelDescriptor::Barrier("AmrLevel::checkPoint(ctml)");

     blob_history_class.open_checkpoint(FullPath);
     blob_history_class.checkpoint(0);
     blob_history_class.close_checkpoint();

     ParallelDescriptor::Barrier("AmrLevel::checkPoint(blob)");

    } else if ((level>=1)&&(level<=max_level)) {

     //do nothing
     
    } else
     amrex::Error("level invalid");

} // end subroutine AmrLevel::checkPoint

AmrLevel::~AmrLevel ()
{

    if (level==0) {
     int time_order=parent->Time_blockingFactor();
     int local_nmat=parent->global_AMR_num_materials;

     for (int i=0;i<=time_order+parent->LSA_extra_data;i++) {

      new_data_FSI[i].clear_FSI();
#ifdef AMREX_PARTICLES
      AmrLevel0_new_dataPC[i].reset();
#endif

     } // for (int i=0;i<=time_order+parent->LSA_extra_data;i++) 

#ifdef AMREX_PARTICLES
     AmrLevel0_new_dataPC.resize(0);
#endif

     new_data_FSI.resize(0);
    } else if (level>0) {
     // do nothing
    } else
     amrex::Error("level invalid");

    parent = 0;
}

#ifdef AMREX_PARTICLES
AmrParticleContainer<N_EXTRA_REAL,N_EXTRA_INT,0,0>& 
   AmrLevel::newDataPC (int slab_index)
{

 AMREX_ALWAYS_ASSERT(level==0);

 int time_order=parent->Time_blockingFactor();
 if ((slab_index<0)||
     (slab_index>time_order+parent->LSA_extra_data)) {
  std::cout << "parent->LSA_extra_data= " << parent->LSA_extra_data << '\n';
  std::cout << "time_order= " << time_order << '\n';
  std::cout << "slab_index= " << slab_index << '\n';
  amrex::Error("slab_index invalid1");
 }

 return *AmrLevel0_new_dataPC[slab_index];

}

const AmrParticleContainer<N_EXTRA_REAL,N_EXTRA_INT,0,0>&
AmrLevel::newDataPC (int slab_index) const
{

 AMREX_ALWAYS_ASSERT(level==0);
 int time_order=parent->Time_blockingFactor();

 if ((slab_index<0)||
     (slab_index>time_order+parent->LSA_extra_data)) {
  std::cout << "parent->LSA_extra_data= " << parent->LSA_extra_data << '\n';
  std::cout << "time_order= " << time_order << '\n';
  std::cout << "slab_index= " << slab_index << '\n';
  amrex::Error("slab_index invalid2");
 }

 return *AmrLevel0_new_dataPC[slab_index];

}

void
AmrLevel::CopyNewToOldPC(int lev_max) {

 AMREX_ALWAYS_ASSERT(level==0);

 int time_order=parent->Time_blockingFactor();

 for (int i=0;i<time_order;i++) {

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
  bool remove_negative=true;

  //amrex-master/Src/Particle/AMReX_Particles.H
  //void Redistribute(int lev_min=0,int lev_max=-1,
  //  int nGrow=0,int local=0,bool remove_negative=true);
  AmrLevel0_new_dataPC[i]->
     Redistribute(lev_min,lev_max,
	nGrow_Redistribute,local_Redistribute,remove_negative);
 } //i=0;i<time_order;i++

} // end subroutine CopyNewToOldPC()

void
AmrLevel::CopyOldToNewPC(int lev_max) {

 AMREX_ALWAYS_ASSERT(level==0);

 int time_order=parent->Time_blockingFactor();

 for (int i=1;i<=time_order;i++) {

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
  bool remove_negative=true;

  AmrLevel0_new_dataPC[i]->
    Redistribute(lev_min,lev_max,
	nGrow_Redistribute,local_Redistribute,remove_negative);
 } // i=1..time_order

} // end subroutine CopyOldToNewPC()

#endif

// The algorithm:
// 1. robust_coarsen:
//   (i) for all Boxes in BoxArray:
//       if box_lo is odd then box_lo--
//       if box_hi+1 is odd then box_hi++
//   (ii) return coarsen(BoxArray)
// 2. find coarsest level needed:
//     i) l=l^max
//     ii) Omega^l=Omega^TBF  TBF=To Be Filled
//     iii) Sigma=Omega^l-Omega^{l,Given}-
//          Omega^{l,physical BC)
//     iv)  if Sigma=null then only 1 level fill needed. bottom_level=l
//     v.a) Omega^{l-1}=robust_coarsen(Omega^{l})
//     v.b) Sigma=Omega^{l-1}-Omega^{l-1,Given}-
//                 coarsened(Omega^{l,Given})-
//                 Omega^{l-1,physical BC}
//     v.c) if Sigma=null then bottom_level=l-1
//     v.d) if Sigma<>null then l=l-1 go back to step v.a
// 3. if l^bottom=l^max then call fill single level
// 4. if l^bottom<l^max then
// 4.a)  for l=l^bottom+1,l^max
// 4.b)   call two level fill using l-1, and l:
//           inputs: Omega^l
//                   Data^{l-1}
//                   old Data^{l}
//           output: Data^{l}    
// 4.c) return Data^{TBF}=Data^{l^max}
void
AmrLevel::FillPatch (int called_from_regrid,
                     AmrLevel & old,
                     MultiFab& mf_to_be_filled,
                     int       dcomp,
                     Real      time,
                     int       index,
                     int       scomp, //absolute within the state (old).
                     int       ncomp,
		     int debug_fillpatch)
{
 BL_PROFILE("AmrLevel::FillPatch()");

 //
 // Must fill this region on crse level and interpolate.
 //
 if (level<0)
  amrex::Error("level invalid in FillPatch");

 BL_ASSERT(ncomp <= (mf_to_be_filled.nComp()-dcomp));
 BL_ASSERT(0 <= index && index < desc_lst.size());

 Vector<int> scompBC_map;
 scompBC_map.resize(ncomp);

 for (int icomp=0;icomp<ncomp;icomp++) {
  scompBC_map[icomp]=scomp+icomp;
 }

 int                     DComp   = dcomp;
 const StateDescriptor&  desc    = old.desc_lst[index];
 IndexType desc_typ(desc.getType());
 int desc_grid_type=-1;
 StateData::get_grid_type(desc_typ,desc_grid_type);

  //scompBC_map[0..ncomp-1] in 0..StateDescriptor.ncomp-1 ?
 desc.check_inRange(scompBC_map,ncomp);

  //scompBC_map[0...ncomp-1] in [scomp...scomp+ncomp-1]
  //ranges in [0...ncomp-1]
 std::vector< std::pair<int,int> > ranges = 
   desc.sameInterps(scompBC_map,ncomp);

 int ncomp_sanity_check=0;

 Vector<StateDataPhysBCFunct*> tower_physbc;
 Vector<PhysBCFunctBaseSUSSMAN*> tower_physbc_base;
 Vector<MultiFab*> tower_data;
 Vector<StateData*> tower_state_data;
 Vector<const Geometry*> tower_geom;
 Vector<Real> tower_nudge_time;
 Vector<int> tower_best_index;
 Vector<int> tower_bfact;

 tower_physbc.resize(level+1);
 tower_data.resize(level+1);
 tower_state_data.resize(level+1);
 tower_geom.resize(level+1);
 tower_nudge_time.resize(level+1);
 tower_best_index.resize(level+1);
 tower_bfact.resize(level+1);

 for (int ilev=0;ilev<=level;ilev++) {
  tower_physbc[ilev]=nullptr;
  tower_data[ilev]=nullptr;
  tower_state_data[ilev]=nullptr;
  tower_geom[ilev]=nullptr;
  tower_nudge_time[ilev]=0.0;
  tower_best_index[ilev]=0;
  tower_bfact[ilev]=0;
 }

 tower_state_data[level]=&old.state[index];
 tower_geom[level]=&old.geom;
  //nudge_time=time_array[best_index]
  //0<=best_index<=bfact_time_order
 tower_state_data[level]->get_time_index(time,
  tower_nudge_time[level],
  tower_best_index[level]);
 tower_data[level]=
   &(tower_state_data[level]->newData(tower_best_index[level]));
 tower_bfact[level]=parent->Space_blockingFactor(level);

 tower_physbc[level]=new StateDataPhysBCFunct(
  *tower_state_data[level],
  *tower_geom[level]);

 if (tower_data[level]->DistributionMap()==old.DistributionMap()) {
  // do nothing
 } else {
  amrex::Error("tower_data->DistributionMap()!=old.DistributionMap()");
 }
 if (tower_data[level]->boxArray().CellEqual(old.boxArray())) {
  // do nothing
 } else {
  amrex::Error("tower_data->boxArray().CellEqual(old.boxArray()) failed");
 }

 for (int ilev=0;ilev<level;ilev++) {
  AmrLevel& clev = parent->getLevel(ilev);
  tower_state_data[ilev] = &(clev.state[index]);
  tower_geom[ilev]=&(clev.geom);
  tower_state_data[ilev]->get_time_index(time,
   tower_nudge_time[ilev],
   tower_best_index[ilev]);

  tower_data[ilev]=&(tower_state_data[ilev]->newData(tower_best_index[ilev]));
  tower_bfact[ilev]=parent->Space_blockingFactor(ilev);

  if (tower_data[ilev]->DistributionMap()==clev.DistributionMap()) {
   // do nothing
  } else {
   amrex::Error("tower_data->DistributionMap()!=clev.DistributionMap()");
  }
  if (tower_data[ilev]->boxArray().CellEqual(clev.boxArray())) {
   // do nothing
  } else {
   amrex::Error("tower_data->boxArray().CellEqual(clev.boxArray()) failed");
  }

  tower_physbc[ilev]=new StateDataPhysBCFunct(
   *tower_state_data[ilev],
   *tower_geom[ilev]);
 }

 for (int ilev=0;ilev<=level;ilev++) {
  tower_physbc_base.push_back(tower_physbc[ilev]);
 }

 for (int ilev=0;ilev<=level;ilev++) {
  if (tower_data[ilev]->nComp()>=scomp+ncomp) {
   //do nothing
  } else {
   std::cout << "ilev=" << ilev << '\n';
   std::cout << "scomp=" << scomp << '\n';
   std::cout << "ncomp=" << ncomp << '\n';
   std::cout << "tower_data[ilev]->nComp()=" << 
     tower_data[ilev]->nComp() << '\n';
   amrex::Error("tower_data nComp invalid");
  }
 }

 for (unsigned int igroup = 0; igroup < ranges.size(); igroup++) {
  const int     scomp_range = ranges[igroup].first;
  const int     ncomp_range = ranges[igroup].second;
  Interpolater* mapper = desc.interp(scompBC_map[scomp_range]);

  if (scomp_range==ncomp_sanity_check) {
   //do nothing
  } else
   amrex::Error("scomp_range failed sanity check");

  ncomp_sanity_check+=ncomp_range;

  Vector<int> local_scompBC_map;
  local_scompBC_map.resize(ncomp_range);
  for (int isub=0;isub<ncomp_range;isub++)
   local_scompBC_map[isub]=scompBC_map[scomp_range+isub];

  int scomp_local=scomp_range+scomp;

    // This routine is NOT in the "Base" directory, instead it is in the
    // present directory in the file: FillPatchUtil.cpp.
    // We are fooled into thinking this is a amrex routine since the
    // code in FillPatchUtil.cpp is hidden within a "namespace amrex"
    // block.
  int ngrow_root=mf_to_be_filled.nGrow();
  amrex::FillPatchTower( //calling from FillPatch.
    ngrow_root,
    called_from_regrid,
    level,
    level,
    mf_to_be_filled, 
    tower_nudge_time[level],
    tower_data,
    scomp_local, //absolute within the state (tower_data[level])
    DComp,
    ncomp_range,
    tower_geom,
    tower_physbc_base,
    mapper,
    desc.getBCs(),  // global_bcs
    local_scompBC_map,
    tower_bfact,
    desc_grid_type,
    debug_fillpatch);

  DComp += ncomp_range;

 } // igroup=0..ranges.size()-1

 if (ncomp==ncomp_sanity_check) {
  //do nothing
 } else
  amrex::Error("ncomp sanity check failed");

 tower_data.clear(); //removes pointers from vector, but doesn't delete data
 tower_state_data.clear();
 tower_geom.clear();

 for (int ilev=0;ilev<=level;ilev++) {
  delete tower_physbc[ilev];
 }
 tower_physbc.clear();
 tower_physbc_base.clear();

}   // end subroutine AmrLevel::FillPatch


void
AmrLevel::FillCoarsePatchGHOST (
    Vector<MultiFab*> tower_data, //scomp..scomp+ncomp-1
    int level_in,
    int ngrow_in,
    Real      time,
    int       index,
    int       scomp, //tower_data: scomp..scomp+ncomp-1
    Vector<int> scompBC_map, // 0..ncomp-1
    int       ncomp,
    int debug_fillpatch)
{
 BL_PROFILE("AmrLevel::FillCoarsePatchGHOST()");

 if (level!=level_in)
  amrex::Error("level <> level_in");

 if (level<=0)
  amrex::Error("level invalid in FillCoarsePatchGHOST");

 if (tower_data[level]->nGrow()==ngrow_in) {
  //do nothing
 } else {
  std::cout << "scomp=" << scomp << '\n';
  std::cout << "ncomp=" << ncomp << '\n';
  std::cout << "ngrow_in=" << ngrow_in << '\n';
  std::cout << "tower_data[level]->nGrow()=" << 
    tower_data[level]->nGrow() << '\n';
  for (int ilev2=0;ilev2<=level;ilev2++) {
   std::cout << "tower_data[ilev2]->nGrow()=" << 
    tower_data[ilev2]->nGrow() << '\n';
   std::cout << "tower_data[ilev2]->nComp()=" << 
     tower_data[ilev2]->nComp() << '\n';
  }
  amrex::Error("tower_data nGrow invalid FillCoarsePatchGHOST");
 }

 for (int ilev=0;ilev<=level;ilev++) {
  if (tower_data[ilev]->nComp()>=scomp+ncomp) {
   //do nothing
  } else {
   std::cout << "ilev=" << ilev << '\n';
   std::cout << "scomp=" << scomp << '\n';
   std::cout << "ncomp=" << ncomp << '\n';
   std::cout << "tower_data[ilev]->nComp()=" << 
     tower_data[ilev]->nComp() << '\n';
   amrex::Error("tower_data nComp invalid");
  }
  if (tower_data[ilev]->nGrow()>=0) {
   //do nothing
  } else {
   std::cout << "ilev=" << ilev << '\n';
   std::cout << "ngrow_in=" << ngrow_in << '\n';
   std::cout << "tower_data[ilev]->nGrow()=" << 
     tower_data[ilev]->nGrow() << '\n';
   amrex::Error("tower_data nGrow invalid FillCoarsePatchGHOST");
  }
 }

 if ((index<0)||(index>=desc_lstGHOST.size()))
  amrex::Error("(index<0)||(index>=desc_lstGHOST.size())");

 if (scompBC_map.size()!=ncomp)
  amrex::Error("scompBC_map has invalid size");

 const BoxArray& mf_BA = tower_data[level]->boxArray();
 if (ngrow_in<0)
  amrex::Error("ngrow_in<0 in FillCoarsePatchGHOST");

 DistributionMapping dm=tower_data[level]->DistributionMap();

 const Box& domain = geom.Domain();
 const int* domlo = domain.loVect();
 const int* domhi = domain.hiVect();

 Vector<StateDataPhysBCFunctGHOST*> tower_physbc;
 Vector<PhysBCFunctBaseSUSSMAN*> tower_physbc_base;
 Vector<StateData*> tower_state_data;
 Vector<const Geometry*> tower_geom;
 Vector<Real> tower_nudge_time;
 Vector<int> tower_best_index;
 Vector<int> tower_bfact;

 tower_physbc.resize(level+1);
 tower_state_data.resize(level+1);
 tower_geom.resize(level+1);
 tower_nudge_time.resize(level+1);
 tower_best_index.resize(level+1);
 tower_bfact.resize(level+1);

 for (int ilev=0;ilev<=level;ilev++) {
  tower_physbc[ilev]=nullptr;
  tower_state_data[ilev]=nullptr;
  tower_geom[ilev]=nullptr;
  tower_nudge_time[ilev]=0.0;
  tower_best_index[ilev]=0;
  tower_bfact[ilev]=0;
 }

 tower_state_data[level]=&state[index];
 tower_geom[level]=&geom;
 tower_state_data[level]->get_time_index(time,
   tower_nudge_time[level],
   tower_best_index[level]);
 tower_bfact[level]=parent->Space_blockingFactor(level);

 tower_physbc[level]=new StateDataPhysBCFunctGHOST(
   *tower_state_data[level],
   *tower_geom[level]);

 for (int ilev=0;ilev<level;ilev++) {
  AmrLevel& clev = parent->getLevel(ilev);
  tower_state_data[ilev] = &(clev.state[index]);
  tower_geom[ilev]=&(clev.geom);
  tower_state_data[ilev]->get_time_index(time,
   tower_nudge_time[ilev],
   tower_best_index[ilev]);

  tower_bfact[ilev]=parent->Space_blockingFactor(ilev);

  tower_physbc[ilev]=new StateDataPhysBCFunctGHOST(
   *tower_state_data[ilev],
   *tower_geom[ilev]);
 }

 for (int ilev=0;ilev<=level;ilev++) {
  tower_physbc_base.push_back(tower_physbc[ilev]);
 }

 int DComp = scomp;

 const StateDescriptor&  descGHOST = desc_lstGHOST[index];
 IndexType desc_typ(descGHOST.getType());
 int desc_grid_type=-1;
 StateData::get_grid_type(desc_typ,desc_grid_type);

  //pdomain will be different depending on whether the state variable
  //is cell centered or staggared.
 const Box& pdomain = tower_state_data[level]->getDomain();
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

  // scompBC_map[0...ncomp-1] in 0 ... StateDescriptor.ncomp-1 ?
 descGHOST.check_inRange(scompBC_map, ncomp);

 std::vector< std::pair<int,int> > ranges = 
   descGHOST.sameInterps(scompBC_map,ncomp);

 int ncomp_sanity_check=0;

 for (unsigned int igroup = 0; igroup < ranges.size(); igroup++) {

  const int     scomp_range  = ranges[igroup].first;
  const int     ncomp_range  = ranges[igroup].second;
  if (igroup==0) {
   if (scomp_range!=0)
    amrex::Error("scomp_range!=0");
  }
  Interpolater* mapper = descGHOST.interp(scompBC_map[scomp_range]);

  if (scomp_range==ncomp_sanity_check) {
   //do nothing
  } else
   amrex::Error("scomp_range failed sanity check");

  ncomp_sanity_check+=ncomp_range;

  BoxArray crseBA(mf_BA.size());
  for (int j = 0, N_CBA = crseBA.size(); j < N_CBA; ++j) {
   BL_ASSERT(mf_BA[j].ixType() == descGHOST.getType());
   const Box& bx = mf_BA[j];

   Box grow_bx(bx);
   const int* bx_lo=bx.loVect();
   const int* bx_hi=bx.hiVect();

   for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
    if (bx_lo[dir]>domlo[dir]) 
     grow_bx.growLo(dir,ngrow_in);
    if (bx_hi[dir]<domhi[dir]) 
     grow_bx.growHi(dir,ngrow_in);
   } //dir=0..sdim-1
   Box grow_bx_test=grow_bx & domain;

   crseBA.set(j,mapper->CoarseBox(grow_bx_test,
     tower_bfact[level-1],
     tower_bfact[level],
     desc_grid_type));
  }

    // ngrow_in=0
    // This is data_to_be_filled.
  MultiFab crseMF(
    crseBA,
    dm,
    ncomp_range,
    0,
    MFInfo().SetTag("crseMF"),
    FArrayBoxFactory());

  int dcomp_data=scomp+scomp_range;

  if (dcomp_data!=DComp)
   amrex::Error("dcomp_data!=DComp");

  Vector<int> local_scompBC_map;
  local_scompBC_map.resize(ncomp_range);
  for (int isub=0;isub<ncomp_range;isub++)
   local_scompBC_map[isub]=scompBC_map[scomp_range+isub];

  int scomp_local=scomp+scomp_range;

  Vector< BCRec > local_bcs;
  const Vector< BCRec> & global_bcs=descGHOST.getBCs();

  local_bcs.resize(ncomp_range);
  for (int isub=0;isub<ncomp_range;isub++)
   local_bcs[isub]=global_bcs[local_scompBC_map[isub]]; 

  int called_from_regrid=1; //disable the sanity check for number of levels
  int ngrow_root=crseMF.nGrow();

  amrex::FillPatchTower( //calling from FillCoarsePatchGHOST
    ngrow_root,
    called_from_regrid, 
    level,
    level-1,
    crseMF, // data to be filled 0..ncomp_range-1
    tower_nudge_time[level-1],
    //level-1 data; scomp_local...scomp_local+ncomp_range-1
    tower_data, 
    scomp_local,
    0, // dstcomp
    ncomp_range,
    tower_geom,
    tower_physbc_base,
    mapper,
    global_bcs,
    local_scompBC_map,
    tower_bfact,
    desc_grid_type,
    debug_fillpatch);

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(tower_data[level]->boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(*tower_data[level],false); mfi.isValid(); ++mfi) {

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
	  
   const Box& mf_local_box=(*tower_data[level])[mfi].box();

   Box dbx_test=dbx;
   for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
    dbx_test.grow(dir,ngrow_in);
   }

   if (dbx_test==mf_local_box) {
    //do nothing
   } else
    amrex::Error("dbx_test==mf_local_box failed");

   Box grow_bx(dbx);
   const int* bx_lo=dbx.loVect();
   const int* bx_hi=dbx.hiVect();

   for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
    if (bx_lo[dir]>domlo[dir]) 
     grow_bx.growLo(dir,ngrow_in);
    if (bx_hi[dir]<domhi[dir]) 
     grow_bx.growHi(dir,ngrow_in);
   } //dir=0..sdim-1
   Box grow_bx_test=grow_bx & domain;

   Array4<Real> const& mf_array=(*tower_data[level])[mfi].array();
   const Dim3 lo3=amrex::lbound(grow_bx_test);
   const Dim3 hi3=amrex::ubound(grow_bx_test);
   for (int n=0;n<ncomp_range;++n) {
   for (int z=lo3.z;z<=hi3.z;++z) {
   for (int y=lo3.y;y<=hi3.y;++y) {
   for (int x=lo3.x;x<=hi3.x;++x) {
    mf_array(x,y,z,n+DComp)=1.0e+20;
   }
   }
   }
   }

   mapper->interp(tower_nudge_time[level-1],
     crseMF[mfi],
     0,  // crse_comp
     (*tower_data[level])[mfi],
     DComp,
     ncomp_range,
     grow_bx_test,
     *tower_geom[level-1],
     *tower_geom[level],
     bcr, // not used.
     level-1,level,
     tower_bfact[level-1],
     tower_bfact[level],
     desc_grid_type);

   for (int n=0;n<ncomp_range;++n) {
   for (int z=lo3.z;z<=hi3.z;++z) {
   for (int y=lo3.y;y<=hi3.y;++y) {
   for (int x=lo3.x;x<=hi3.x;++x) {
    Real test_norm=mf_array(x,y,z,n+DComp);
    if (test_norm<1.0e+19) {
     //do nothing
    } else {
     amrex::Error("test_norm<1.0e+19 failed");
   }
   }
   }
   }
   }

  }  // mfi
} // omp
  thread_class::sync_tile_d_numPts();
  ParallelDescriptor::ReduceRealSum(thread_class::tile_d_numPts[0]);
  thread_class::reconcile_d_numPts(LOOP_MAPPER_INTERP_GHOST,
      "FillCoarsePatchGHOST");

  tower_data[level]->FillBoundary(DComp,ncomp_range,geom.periodicity());

  tower_physbc[level]->FillBoundary(level,
    *tower_data[level],
    tower_nudge_time[level],DComp,
    local_scompBC_map,ncomp_range,
    tower_bfact[level]);

  DComp += ncomp_range;

 } // igroup=0..ranges.size()-1

 if (ncomp==ncomp_sanity_check) {
  //do nothing
 } else
  amrex::Error("ncomp sanity check failed");

 tower_state_data.clear();
 tower_geom.clear();

 for (int ilev=0;ilev<=level;ilev++) {
  delete tower_physbc[ilev];
 }
 tower_physbc.clear();
 tower_physbc_base.clear();

}   // end subroutine FillCoarsePatchGHOST

//Use virtual (aka GHOST) state data information for filling physical
//boundaries. i.e. The virtual (GHOST) state data holds no actual data, only
//holds pointers for filling the virtual data.
//called from NavierStokes3.cpp
void
AmrLevel::InterpBordersGHOST (
  Vector<MultiFab*> tower_data,
  int level_in,
  int ngrow_in,
  Real      time,
  int       index,
  int       scomp, // source comp wrt tower_data 
  Vector<int> scompBC_map, //scompBC_map[0..ncomp-1]
  int       ncomp,
  int debug_fillpatch)
{
 BL_PROFILE("AmrLevel::InterpBordersGHOST()");

 if (level!=level_in)
  amrex::Error("level <> level_in");

 if (level<0)
  amrex::Error("level invalid in InterpBordersGHOST");


 if (tower_data[level]->nGrow()==ngrow_in) {
  //do nothing
 } else {
  std::cout << "scomp=" << scomp << '\n';
  std::cout << "ncomp=" << ncomp << '\n';
  std::cout << "ngrow_in=" << ngrow_in << '\n';
  std::cout << "tower_data[level]->nGrow()=" << 
    tower_data[level]->nGrow() << '\n';
  for (int ilev2=0;ilev2<=level;ilev2++) {
   std::cout << "tower_data[ilev2]->nGrow()=" << 
    tower_data[ilev2]->nGrow() << '\n';
   std::cout << "tower_data[ilev2]->nComp()=" << 
     tower_data[ilev2]->nComp() << '\n';
  }
  amrex::Error("tower_data nGrow invalid InterpBordersGHOST");
 }

 for (int ilev=0;ilev<=level;ilev++) {
  if (tower_data[ilev]->nComp()>=scomp+ncomp) {
   //do nothing
  } else {
   std::cout << "ilev=" << ilev << '\n';
   std::cout << "scomp=" << scomp << '\n';
   std::cout << "ncomp=" << ncomp << '\n';
   std::cout << "tower_data[ilev]->nComp()=" << 
     tower_data[ilev]->nComp() << '\n';
   amrex::Error("tower_data nComp invalid");
  }
  if (tower_data[ilev]->nGrow()>=0) {
   //do nothing
  } else {
   std::cout << "ilev=" << ilev << '\n';
   std::cout << "scomp=" << scomp << '\n';
   std::cout << "ncomp=" << ncomp << '\n';
   std::cout << "ngrow_in=" << ngrow_in << '\n';
   std::cout << "tower_data[ilev]->nGrow()=" << 
     tower_data[ilev]->nGrow() << '\n';
   for (int ilev2=0;ilev2<=level;ilev2++) {
    std::cout << "tower_data[ilev2]->nGrow()=" << 
     tower_data[ilev2]->nGrow() << '\n';
    std::cout << "tower_data[ilev2]->nComp()=" << 
     tower_data[ilev2]->nComp() << '\n';
   }
   amrex::Error("tower_data nGrow invalid InterpBordersGHOST");
  }
 }

 BL_ASSERT((0<=index)&&(index<desc_lstGHOST.size()));

 if (scompBC_map.size()!=ncomp)
  amrex::Error("scompBC_map has invalid size");

 const BoxArray& mf_BA = tower_data[level]->boxArray();

 if (ngrow_in<0)
  amrex::Error("ngrow_in<0 in InterpBordersGHOST");

 DistributionMapping dm=tower_data[level]->DistributionMap();

 Vector<StateDataPhysBCFunctGHOST*> tower_physbc;
 Vector<PhysBCFunctBaseSUSSMAN*> tower_physbc_base;
 Vector<StateData*> tower_state_data;
 Vector<const Geometry*> tower_geom;
 Vector<Real> tower_nudge_time;
 Vector<int> tower_best_index;
 Vector<int> tower_bfact;

 tower_physbc.resize(level+1);
 tower_state_data.resize(level+1);
 tower_geom.resize(level+1);
 tower_nudge_time.resize(level+1);
 tower_best_index.resize(level+1);
 tower_bfact.resize(level+1);

 for (int ilev=0;ilev<=level;ilev++) {
  tower_physbc[ilev]=nullptr;
  tower_state_data[ilev]=nullptr;
  tower_geom[ilev]=nullptr;
  tower_nudge_time[ilev]=0.0;
  tower_best_index[ilev]=0;
  tower_bfact[ilev]=0;
 }

 tower_state_data[level]=&state[index];
 tower_geom[level]=&geom;
  //nudge_time=time_array[best_index]
  //0<=best_index<=bfact_time_order
 tower_state_data[level]->get_time_index(time,
  tower_nudge_time[level],
  tower_best_index[level]);
 tower_bfact[level]=parent->Space_blockingFactor(level);

 tower_physbc[level]=new StateDataPhysBCFunctGHOST(
  *tower_state_data[level],
  *tower_geom[level]);

 if (tower_data[level]->DistributionMap()==DistributionMap()) {
  // do nothing
 } else {
  amrex::Error("tower_data->DistributionMap()!=DistributionMap()");
 }
 if (tower_data[level]->boxArray().CellEqual(boxArray())) {
  // do nothing
 } else {
  amrex::Error("tower_data->boxArray().CellEqual(boxArray()) failed");
 }

 for (int ilev=0;ilev<level;ilev++) {
  AmrLevel& clev = parent->getLevel(ilev);
  tower_state_data[ilev] = &(clev.state[index]);
  tower_geom[ilev]=&(clev.geom);
  tower_state_data[ilev]->get_time_index(time,
   tower_nudge_time[ilev],
   tower_best_index[ilev]);

  tower_bfact[ilev]=parent->Space_blockingFactor(ilev);

  if (tower_data[ilev]->DistributionMap()==clev.DistributionMap()) {
   // do nothing
  } else {
   amrex::Error("tower_data->DistributionMap()!=clev.DistributionMap()");
  }
  if (tower_data[ilev]->boxArray().CellEqual(clev.boxArray())) {
   // do nothing
  } else {
   amrex::Error("tower_data->boxArray().CellEqual(clev.boxArray()) failed");
  }

  tower_physbc[ilev]=new StateDataPhysBCFunctGHOST(
   *tower_state_data[ilev],
   *tower_geom[ilev]);
 }

 for (int ilev=0;ilev<=level;ilev++) {
  tower_physbc_base.push_back(tower_physbc[ilev]);
 }

 MultiFab mf_to_be_filled(mf_BA,dm,ncomp,ngrow_in,
   MFInfo().SetTag("mf_to_be_filled"),FArrayBoxFactory());

  // dstmf,srcmf,srccomp,dstcomp,ncomp,ngrow
 MultiFab::Copy(mf_to_be_filled,*tower_data[level],scomp,0,ncomp,0);

 int DComp = 0;

 const StateDescriptor& descGHOST = desc_lstGHOST[index];
 IndexType desc_typ(descGHOST.getType());
 int desc_grid_type=-1;
 StateData::get_grid_type(desc_typ,desc_grid_type);

  //scompBC_map[0...ncomp-1] in 0...StateDescriptor.ncomp-1 ?
 descGHOST.check_inRange(scompBC_map, ncomp);

  //scompBC_map[0..ncomp-1]
  //ranges in [0..ncomp-1]
 std::vector< std::pair<int,int> > ranges = 
   descGHOST.sameInterps(scompBC_map,ncomp);

 for (unsigned int igroup = 0; igroup < ranges.size(); igroup++) {
  const int     scomp_range  = ranges[igroup].first;
  const int     ncomp_range  = ranges[igroup].second;
  Interpolater* mapper = descGHOST.interp(scompBC_map[scomp_range]);

  int dcomp_data=scomp_range;
  int scomp_data=scomp+scomp_range;

  if (dcomp_data!=DComp)
   amrex::Error("dcomp_data!=DComp");

  Vector<int> local_scompBC_map;
  local_scompBC_map.resize(ncomp_range);
  for (int isub=0;isub<ncomp_range;isub++)
   local_scompBC_map[isub]=scompBC_map[scomp_range+isub];

  int called_from_regrid=0;
  int ngrow_root=mf_to_be_filled.nGrow();

  amrex::FillPatchTower( //calling from InterpBordersGHOST
    ngrow_root,
    called_from_regrid,
    level,
    level,
    mf_to_be_filled,
    tower_nudge_time[level],
    tower_data,
    scomp_data,
    dcomp_data,
    ncomp_range,
    tower_geom,
    tower_physbc_base,
    mapper,
    descGHOST.getBCs(), // global_bcs
    local_scompBC_map,
    tower_bfact,
    desc_grid_type,
    debug_fillpatch);

  DComp += ncomp_range;

 } // igroup=0..ranges.size()-1

 tower_state_data.clear();
 tower_geom.clear();
 for (int ilev=0;ilev<=level;ilev++) {
  delete tower_physbc[ilev];
 }
 tower_physbc.clear();
 tower_physbc_base.clear();

  // dstmf,srcmf,srccomp,dstcomp,ncomp,ngrow
 MultiFab::Copy(*tower_data[level],mf_to_be_filled,0,scomp,ncomp,ngrow_in);

}   // end subroutine InterpBordersGHOST

//called from MacProj.cpp, NavierStokes3.cpp
void
AmrLevel::InterpBorders (
  Vector<MultiFab*> tower_data,
  int level_in,
  int ngrow_in,
  Real      time,
  int       index,
  int       scomp, //source comp wrt tower_data
  Vector<int> scompBC_map, //scompBC_map[0..ncomp-1]
  int       ncomp,
  int debug_fillpatch)
{
 BL_PROFILE("AmrLevel::InterpBorders()");

 if (level!=level_in)
  amrex::Error("level <> level_in");

 if (level<0)
  amrex::Error("level invalid in InterpBorders");

 if (tower_data[level]->nGrow()==ngrow_in) {
  //do nothing
 } else {
  std::cout << "scomp=" << scomp << '\n';
  std::cout << "ncomp=" << ncomp << '\n';
  std::cout << "ngrow_in=" << ngrow_in << '\n';
  std::cout << "tower_data[level]->nGrow()=" << 
    tower_data[level]->nGrow() << '\n';
  for (int ilev2=0;ilev2<=level;ilev2++) {
   std::cout << "tower_data[ilev2]->nGrow()=" << 
    tower_data[ilev2]->nGrow() << '\n';
   std::cout << "tower_data[ilev2]->nComp()=" << 
     tower_data[ilev2]->nComp() << '\n';
  }
  amrex::Error("tower_data nGrow invalid InterpBorders");
 }

 for (int ilev=0;ilev<=level;ilev++) {
  if (tower_data[ilev]->nComp()>=scomp+ncomp) {
   //do nothing
  } else {
   std::cout << "ilev=" << ilev << '\n';
   std::cout << "scomp=" << scomp << '\n';
   std::cout << "ncomp=" << ncomp << '\n';
   std::cout << "tower_data[ilev]->nComp()=" << 
     tower_data[ilev]->nComp() << '\n';
   amrex::Error("tower_data nComp invalid");
  }
  if (tower_data[ilev]->nGrow()>=0) {
   //do nothing
  } else {
   std::cout << "ilev=" << ilev << '\n';
   std::cout << "ngrow_in=" << ngrow_in << '\n';
   std::cout << "tower_data[ilev]->nGrow()=" << 
     tower_data[ilev]->nGrow() << '\n';
   amrex::Error("tower_data nGrow invalid InterpBorders");
  }
 }

 BL_ASSERT(0 <= index && index < desc_lst.size());

 if (scompBC_map.size()!=ncomp)
  amrex::Error("scompBC_map has invalid size");

 const BoxArray& mf_BA = tower_data[level]->boxArray();

 if (ngrow_in<0)
  amrex::Error("ngrow_in<0 in InterpBorders");

 DistributionMapping dm=tower_data[level]->DistributionMap();

 Vector<StateDataPhysBCFunct*> tower_physbc;
 Vector<PhysBCFunctBaseSUSSMAN*> tower_physbc_base;
 Vector<StateData*> tower_state_data;
 Vector<const Geometry*> tower_geom;
 Vector<Real> tower_nudge_time;
 Vector<int> tower_best_index;
 Vector<int> tower_bfact;

 tower_physbc.resize(level+1);
 tower_state_data.resize(level+1);
 tower_geom.resize(level+1);
 tower_nudge_time.resize(level+1);
 tower_best_index.resize(level+1);
 tower_bfact.resize(level+1);

 for (int ilev=0;ilev<=level;ilev++) {
  tower_physbc[ilev]=nullptr;
  tower_state_data[ilev]=nullptr;
  tower_geom[ilev]=nullptr;
  tower_nudge_time[ilev]=0.0;
  tower_best_index[ilev]=0;
  tower_bfact[ilev]=0;
 }

 tower_state_data[level]=&state[index];
 tower_geom[level]=&geom;
  //nudge_time=time_array[best_index]
  //0<=best_index<=bfact_time_order
 tower_state_data[level]->get_time_index(time,
  tower_nudge_time[level],
  tower_best_index[level]);
 tower_bfact[level]=parent->Space_blockingFactor(level);

 tower_physbc[level]=new StateDataPhysBCFunct(
  *tower_state_data[level],
  *tower_geom[level]);

 if (tower_data[level]->DistributionMap()==DistributionMap()) {
  // do nothing
 } else {
  amrex::Error("tower_data->DistributionMap()!=DistributionMap()");
 }
 if (tower_data[level]->boxArray().CellEqual(boxArray())) {
  // do nothing
 } else {
  amrex::Error("tower_data->boxArray().CellEqual(boxArray()) failed");
 }

 for (int ilev=0;ilev<level;ilev++) {
  AmrLevel& clev = parent->getLevel(ilev);
  tower_state_data[ilev] = &(clev.state[index]);
  tower_geom[ilev]=&(clev.geom);
  tower_state_data[ilev]->get_time_index(time,
   tower_nudge_time[ilev],
   tower_best_index[ilev]);

  tower_bfact[ilev]=parent->Space_blockingFactor(ilev);

  if (tower_data[ilev]->DistributionMap()==clev.DistributionMap()) {
   // do nothing
  } else {
   amrex::Error("tower_data->DistributionMap()!=clev.DistributionMap()");
  }
  if (tower_data[ilev]->boxArray().CellEqual(clev.boxArray())) {
   // do nothing
  } else {
   amrex::Error("tower_data->boxArray().CellEqual(clev.boxArray()) failed");
  }

  tower_physbc[ilev]=new StateDataPhysBCFunct(
   *tower_state_data[ilev],
   *tower_geom[ilev]);
 }

 for (int ilev=0;ilev<=level;ilev++) {
  tower_physbc_base.push_back(tower_physbc[ilev]);
 }

 MultiFab mf_to_be_filled(mf_BA,dm,ncomp,ngrow_in,
   MFInfo().SetTag("mf_to_be_filled"),FArrayBoxFactory());

  // dstmf,srcmf,srccomp,dstcomp,ncomp,ngrow
 MultiFab::Copy(mf_to_be_filled,*tower_data[level],scomp,0,ncomp,0);

 int DComp = 0;

 const StateDescriptor&  desc    = desc_lst[index];
 IndexType desc_typ(desc.getType());
 int desc_grid_type=-1;
 StateData::get_grid_type(desc_typ,desc_grid_type);

  //scompBC_map[0..ncomp-1] in 0...StateDescriptor.ncomp-1 ?
 desc.check_inRange(scompBC_map, ncomp);

  //scompBC_map[0..ncomp-1]
  //ranges in [0..ncomp-1]
 std::vector< std::pair<int,int> > ranges = 
   desc.sameInterps(scompBC_map,ncomp);

 for (unsigned int igroup = 0; igroup < ranges.size(); igroup++) {
  const int     scomp_range  = ranges[igroup].first;
  const int     ncomp_range  = ranges[igroup].second;
  Interpolater* mapper = desc.interp(scompBC_map[scomp_range]);

  int scomp_data=scomp+scomp_range;
  int dcomp_data=scomp_range;

  if (dcomp_data!=DComp)
   amrex::Error("dcomp_data!=DComp");

  Vector<int> local_scompBC_map;
  local_scompBC_map.resize(ncomp_range);
  for (int isub=0;isub<ncomp_range;isub++)
   local_scompBC_map[isub]=scompBC_map[scomp_range+isub];

  int called_from_regrid=0;
  int ngrow_root=mf_to_be_filled.nGrow();

  amrex::FillPatchTower( //calling from InterpBorders
    ngrow_root,
    called_from_regrid,
    level,
    level,
    mf_to_be_filled,
    tower_nudge_time[level],
    tower_data,
    scomp_data,
    dcomp_data,
    ncomp_range,
    tower_geom,
    tower_physbc_base,
    mapper,
    desc.getBCs(), // global_bcs
    local_scompBC_map,
    tower_bfact,
    desc_grid_type,
    debug_fillpatch);

  DComp += ncomp_range;
 } // igroup=0..ranges.size()-1

 tower_state_data.clear();
 tower_geom.clear();
 for (int ilev=0;ilev<=level;ilev++) {
  delete tower_physbc[ilev];
 }
 tower_physbc.clear();
 tower_physbc_base.clear();

  // dstmf,srcmf,srccomp,dstcomp,ncomp,ngrow
 MultiFab::Copy(*tower_data[level],mf_to_be_filled,0,scomp,ncomp,ngrow_in);

}   // end subroutine InterpBorders


void
AmrLevel::FillCoarsePatch (MultiFab& mf_to_be_filled,
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
 if (ncomp>(mf_to_be_filled.nComp()-dcomp))
  amrex::Error("ncomp>(mf_to_be_filled.nComp()-dcomp)");
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

 int ngrow=mf_to_be_filled.nGrow();

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

 Vector<StateDataPhysBCFunct*> tower_physbc;
 Vector<PhysBCFunctBaseSUSSMAN*> tower_physbc_base;
 Vector<MultiFab*> tower_data;
 Vector<StateData*> tower_state_data;
 Vector<const Geometry*> tower_geom;
 Vector<Real> tower_nudge_time;
 Vector<int> tower_best_index;
 Vector<int> tower_bfact;

 tower_physbc.resize(level+1);
 tower_data.resize(level+1);
 tower_state_data.resize(level+1);
 tower_geom.resize(level+1);
 tower_nudge_time.resize(level+1);
 tower_best_index.resize(level+1);
 tower_bfact.resize(level+1);

 for (int ilev=0;ilev<=level;ilev++) {
  tower_physbc[ilev]=nullptr;
  tower_data[ilev]=nullptr;
  tower_state_data[ilev]=nullptr;
  tower_geom[ilev]=nullptr;
  tower_nudge_time[ilev]=0.0;
  tower_best_index[ilev]=0;
  tower_bfact[ilev]=0;
 }

 tower_state_data[level]=&state[index];
 tower_geom[level]=&geom;
 tower_state_data[level]->get_time_index(time,
   tower_nudge_time[level],
   tower_best_index[level]);
 tower_data[level]=&mf_to_be_filled;
 tower_bfact[level]=parent->Space_blockingFactor(level);

 tower_physbc[level]=new StateDataPhysBCFunct(
   *tower_state_data[level],
   *tower_geom[level]);

 for (int ilev=0;ilev<level;ilev++) {
  AmrLevel& clev = parent->getLevel(ilev);
  tower_state_data[ilev] = &(clev.state[index]);
  tower_geom[ilev]=&(clev.geom);
  tower_state_data[ilev]->get_time_index(time,
   tower_nudge_time[ilev],
   tower_best_index[ilev]);

  tower_data[ilev]=&(tower_state_data[ilev]->newData(tower_best_index[ilev]));
  tower_bfact[ilev]=parent->Space_blockingFactor(ilev);

  tower_physbc[ilev]=new StateDataPhysBCFunct(
   *tower_state_data[ilev],
   *tower_geom[ilev]);
 }

 for (int ilev=0;ilev<=level;ilev++) {
  tower_physbc_base.push_back(tower_physbc[ilev]);
 }

 for (int ilev=0;ilev<=level;ilev++) {
  if (tower_data[ilev]->nComp()>=scomp+ncomp) {
   //do nothing
  } else {
   std::cout << "ilev=" << ilev << '\n';
   std::cout << "scomp=" << scomp << '\n';
   std::cout << "ncomp=" << ncomp << '\n';
   std::cout << "tower_data[ilev]->nComp()=" << 
     tower_data[ilev]->nComp() << '\n';
   amrex::Error("tower_data nComp invalid");
  }
 }

  //pdomain will be different depending on whether the state variable
  //is cell centered or staggared.
 const Box& pdomain = tower_state_data[level]->getDomain();
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

 const BoxArray& mf_BA = tower_data[level]->boxArray();
 DistributionMapping dm=tower_data[level]->DistributionMap();

  //scompBC_map[0..ncomp-1] in 0...StateDescriptor.ncomp-1 ?
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

   crseBA.set(j,mapper->CoarseBox(grow_bx,
     tower_bfact[level-1],
     tower_bfact[level],
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

  Vector< BCRec > local_bcs;
  const Vector< BCRec> & global_bcs=desc.getBCs();

  local_bcs.resize(ncomp_range);
  for (int isub=0;isub<ncomp_range;isub++)
   local_bcs[isub]=global_bcs[local_scompBC_map[isub]]; 

  int called_from_regrid=1;  //force sanity check to be disabled
  int ngrow_root=crseMF.nGrow();

  amrex::FillPatchTower( //calling from FillCoarsePatch
    ngrow_root,
    called_from_regrid,
    level,
    level-1,
    crseMF, // data to be filled 0..ncomp_range-1
    tower_nudge_time[level-1],
    //level-1 data (smf); scomp_local...scomp_local+ncomp_range-1
    tower_data, 
    scomp_local,
    0, // dstcomp
    ncomp_range,
    tower_geom,
    tower_physbc_base,
    mapper,
    global_bcs,
    local_scompBC_map,
    tower_bfact,
    desc_grid_type,
    debug_fillpatch);

  if (thread_class::nthreads<1)
   amrex::Error("thread_class::nthreads invalid");
  thread_class::init_d_numPts(mf_to_be_filled.boxArray().d_numPts());

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  for (MFIter mfi(*(tower_data[level]),false); mfi.isValid(); ++mfi) {

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

   const Box& mf_local_box=(*tower_data[level])[mfi].box();

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

   Array4<Real> const& mf_array=(*tower_data[level])[mfi].array();
   const Dim3 lo3=amrex::lbound(grow_bx);
   const Dim3 hi3=amrex::ubound(grow_bx);
   for (int n=0;n<ncomp_range;++n) {
   for (int z=lo3.z;z<=hi3.z;++z) {
   for (int y=lo3.y;y<=hi3.y;++y) {
   for (int x=lo3.x;x<=hi3.x;++x) {
    mf_array(x,y,z,n+DComp)=1.0e+20;
   }
   }
   }
   }

   mapper->interp(tower_nudge_time[level-1],
                  crseMF[mfi],
                  0,  // crse_comp
		  (*tower_data[level])[mfi],
		  DComp,
		  ncomp_range,
		  grow_bx,
                  *tower_geom[level-1],
                  *tower_geom[level],
		  bcr,
                  level-1,level,
		  tower_bfact[level-1],
                  tower_bfact[level],
                  desc_grid_type);

   for (int n=0;n<ncomp_range;++n) {
   for (int z=lo3.z;z<=hi3.z;++z) {
   for (int y=lo3.y;y<=hi3.y;++y) {
   for (int x=lo3.x;x<=hi3.x;++x) {
    Real test_norm=mf_array(x,y,z,n+DComp);
    if (test_norm<1.0e+19) {
     //do nothing
    } else {
     amrex::Error("test_norm<1.0e+19 failed");
   }
   }
   }
   }
   }

  }  // mfi
} // omp
  thread_class::sync_tile_d_numPts();
  ParallelDescriptor::ReduceRealSum(thread_class::tile_d_numPts[0]);
  thread_class::reconcile_d_numPts(LOOP_FILLCOARSEPATCH,"FillCoarsePatch");

  tower_data[level]->FillBoundary(DComp,ncomp_range,geom.periodicity());

  tower_physbc[level]->FillBoundary(level,
    *tower_data[level],
    tower_nudge_time[level],DComp,
    local_scompBC_map,ncomp_range,
    tower_bfact[level]);

  DComp += ncomp_range;
 } // igroup=0..ranges.size()-1

 if (DComp==dcomp+ncomp) {
  //do nothing
 } else
  amrex::Error("expecting DComp==dcomp+ncomp");

 tower_data.clear(); //removes pointers from vector, but doesn't delete data
 tower_state_data.clear();
 tower_geom.clear();

 for (int ilev=0;ilev<=level;ilev++) {
  delete tower_physbc[ilev];
 }
 tower_physbc.clear();
 tower_physbc_base.clear();

}   // end subroutine FillCoarsePatch

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

