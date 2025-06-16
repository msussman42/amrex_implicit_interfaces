
#ifndef AMREX_AmrLevel_H_
#define AMREX_AmrLevel_H_ 

#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_RealBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Array.H>
#include <AMReX_FabArray.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <AMReX_VisMF.H>
#include <Interpolater.H>
#include <AMReX_AmrCore.H>
#include <StateDescriptor.H>
#include <StateData.H>

#ifdef AMREX_PARTICLES
#include <AMReX_Particles.H>
#include <AMReX_AmrParticles.H>
#include <AMReX_NeighborParticles.H>
#endif

#include <N_EXTRA_REAL.H>

#include <memory>
#include <map>

namespace amrex {

class TagBox;
class TagBoxArray;

//
// Virtual base class for managing individual levels.
//
// AmrLevel functions both as a container for state data on a level
// and also manages the advancement of data in time.
//

//new_data_FSI only lives at level=0
class FSI_container_class
{
public:

int CTML_num_solids;
int FSI_num_scalars;
int max_num_nodes[3];
int max_num_elements;
int structured_flag; //0=unstructured 1=structured
int structure_dim; //AMREX_SPACEDIM
int structure_topology;  //0=filament (2d) 1=sheet (3d) 2=volumetric
int ngrow_node;
//node[0..num_solids-1][0..max_num_nodes_total_grow-1][0..2]
Vector< Real > node_list; 
Vector< Real > prev_node_list; 
//velocity[0..num_solids-1][0..max_num_nodes_total_grow-1][0..2]
Vector< Real > velocity_list; 
Vector< Real > prev_velocity_list; 
//element[0..num_solids-1][0..max_num_elements-1][0..3]
Vector< int > element_list;
//init_node_list[0..num_solids-1][0..max_num_nodes_total_grow-1][0..2]
Vector< Real > init_node_list; 
//mass[0..num_solids-1][0..max_num_nodes_total_grow-1]
Vector< Real > mass_list; 
//temp[0..num_solids-1][0..max_num_nodes_total_grow-1] (temperature)
Vector< Real > temp_list; 
//scalar[0..num_solids-1][0..max_num_nodes_total_grow-1][0..FSI_num_scalars-1]
Vector< Real > scalar_list; 
Vector< Real > prev_scalar_list; 

void FSI_flatten(Vector< Real >& flattened_data);
void FSI_unflatten(Vector< Real > flattened_data);

static std::ofstream CTML_checkpoint_file;

void open_checkpoint(const std::string& FullPath);
void close_checkpoint();
void checkpoint(int check_id);
void restart(int check_id,std::istream& is);

void initData_FSI(
  const int CTML_num_solids_init,
  const int max_num_nodes_init[3],
  const int max_num_elements_init,
  const int FSI_num_scalars_init,
  const int structured_flag_init=1, //structured
  const int structure_dim_init=AMREX_SPACEDIM,
  const int structure_topology_init=AMREX_SPACEDIM-2,//filament 2d,sheet 3d
  const int ngrow_node_init=2); //StructureCodeShoele

void copyFrom_FSI(const FSI_container_class& source_FSI);
void clear_FSI();

FSI_container_class() { 

  CTML_num_solids=0;
  max_num_nodes[0]=0;
  max_num_nodes[1]=0;
  max_num_nodes[2]=0;
  max_num_elements=0;
  FSI_num_scalars=0;
  structured_flag=1; //structured
  structure_dim=AMREX_SPACEDIM;
  structure_topology=AMREX_SPACEDIM-2; //filament 2d, sheet 3d
  ngrow_node=2;

}

};


class snapshot_blobclass
{
public:
 
 Real blob_volume;
 Real blob_center[AMREX_SPACEDIM];
 Real blob_axis_len[AMREX_SPACEDIM];
 //columns of "evecs" are the eigenvectors.
 //evecs(i,j) flattened index is s=n*(j-1)+i
 Real blob_axis_evec[AMREX_SPACEDIM*AMREX_SPACEDIM];
 Real blob_time;
 int blob_step;

 void sort_axis();

};

class dynamic_blobclass
{
public:

 int im;
 Real start_time;
 Real end_time;
 int start_step;
 int end_step;
 Vector<snapshot_blobclass> snapshots;

};

class dynamic_blobclass_array
{
public:

 static std::ofstream blob_checkpoint_file;

 void open_checkpoint(const std::string& FullPath);
 void close_checkpoint();
 void checkpoint(int check_id);
 void restart(int check_id,std::istream& is);

 Real start_time;
 Real end_time;
 int start_step;
 int end_step;

 Vector<dynamic_blobclass> blob_history;

};

class blobclass
{
public:
   // X=(x,y,z)  U=(u,v,w)
   // velocity=sum alpha_i phi_i(X)  i=1..6 in 3D (i=1..4 in 2D)
   // E=sum || u_k - sum alpha_i phi_i(X_k) ||^2 w_k
   // least squares matrix 
   // Aij=sum phi_i(X_k) dot phi_j(X_k) w_k
   // Bi=sum phi_i(X_k) dot u_k w_k
   // ux,uy,uz fluid_interior,fluid,solid
  Real blob_matrix[3*(2*AMREX_SPACEDIM)*(2*AMREX_SPACEDIM)]; 
  Real blob_RHS[3*(2*AMREX_SPACEDIM)]; 
   // alphax,alphay,alphaz,alphaxy alphaxz alphayz  
   //  fluid interior,fluid,solid
   // u=alphaxy(y-y0)+alphaxz(z-z0)+alphax  
   // v=-alphaxy(x-x0)+alphayz(z-z0)+alphay
   // w=-alphayz(y-y0)-alphaxz(x-x0)+alphaz
   // phix=(1 0 0) phiy=(0 1 0) phiz=(0 0 1)
   // phixy=( (y-y0)  -(x-x0)   0   )
   // phixz=( (z-z0)     0  -(x-x0) )
   // phiyz=(   0    (z-z0) -(y-y0) )
  Real blob_velocity[3*(2*AMREX_SPACEDIM)]; 
   // translational momentum integral rho u dv
   // angular momentum integral (x-x0) cross rho u dv
   // | i    j    k  |
   // |x-x0 y-y0 z-z0|
   // | u    v    w  |  = ((y-y0)w-(z-z0)v)ihat+((z-z0)u-(x-x0)w))jhat+
   //                     ((x-x0)v-(y-y0)u)khat
   // first group: sum phi_i(x_k) dot u_k wmom_k
   // second group: sum phi_i(x_k) dot phi_i(x_k) wmom_k
   // note:
   // momentum_x=alphax sum_k w_k=alphax sum phi_1 dot phi_1 wmom_k
   // momentum_y and momentum_z are analogous.
   // momentum_khat=-alphaxy sum ((xk-x0)^2 wk + (yk-y0)^2) wmomk=
   //               -alphaxy sum phixy dot phixy wmomk
   // momentum_jhat=alphaxz sum (zk-z0)^2 + (xk-x0)^2 wmomk=
   //               alphaxz sum phixz dot phixz wmomk
   // momentum_ihat=-alphayz sum (zk-z0)^2 + (yk-y0)^2 wmomk=
   //               -alphayz sum phiyz dot phiyz wmomk
  Real blob_integral_momentum[2*(2*AMREX_SPACEDIM)];
  Real blob_energy;
   // fluid_interior,fluid,solid
  Real blob_mass_for_velocity[3];
  Real blob_volume;
  Real blob_cell_count;
  Real blob_cellvol_count;
  Real blob_mass;
  Real blob_pressure;
  Real blob_second_moment[6]; //xx,xy,x(y|z),yy,y(y|z),(yy|zz)
  Real blob_center_integral[AMREX_SPACEDIM];
  Real blob_center_actual[AMREX_SPACEDIM];
  Real blob_perim;  // surface area of the blob
   // surface area of region separating the blob from a given material im.
  Vector<Real> blob_perim_mat;  
   // perimeter of contact line between im, im1, im2
  Vector< Vector<Real> > blob_triple_perim;
  int im;  // 1<=im<=num_materials
};



class AmrLevel
{

public:
    //
    // The destructor.
    //
    virtual ~AmrLevel ();
    //
    //
    virtual void writePlotFile (
     int do_plot,int do_slice,
     int SDC_outer_sweeps_in,
     int slab_step_in,
     int divu_outer_sweeps_in) = 0;

    //
    // Write current state to checkpoint file.
    //
    virtual void checkPoint (const std::string& dir,
                             std::ostream&      os);
    //
    // Restart from a checkpoint file.
    //
    virtual void restart (AmrCore&      papa,
                          std::istream& is,
			  int old_finest_level,
			  int new_finest_level);
    //
    // Is name a state variable?
    //
    static bool isStateVariable (const std::string& name,
                                int&               state_indx,
                                int&               ncomp);

    //
    // Compute the initial time step.  This is a pure virtual function
    // and hence MUST be implemented by derived classes.
    //
    virtual void computeInitialDt (int finest_level,
        Real& dt,Real stop_time)=0;
    //
    // Compute the next time step.  This is a pure virtual function
    // and hence MUST be implemented by derived classes.
    //
    virtual void computeNewDt (int finest_level,
       Real& dt,Real stop_time)=0;
    //
    // Do an integration step on this level.  Returns maximum safe
    // time step.  This is a pure virtual function and hence MUST
    // be implemented by derived classes.
    //
    virtual Real advance (Real time,Real dt) = 0;
    //
    // Contains operations to be done after a timestep.  This is a
    // pure virtual function and hence MUST be implemented by derived
    // classes.
    //
    virtual  void post_timestep (Real stop_time) = 0;
    //
    // Operations to be done after restart.  This is a pure virtual
    // function and hence MUST be implemented by derived classes.
    virtual void post_restart () = 0;

    virtual void post_regrid (int lbase,
         int start_level,int new_finest,
         int initialInit_flag,Real time) = 0;
    //
    // Operations to be done after initialization.
    // This is a pure virtual function and hence MUST be
    // implemented by derived classes.
    //
    virtual  void post_init (Real stop_time) = 0;
    //
    // Is it ok to continue the calculation?
    // This is a pure virtual function and hence MUST be
    // implemented by derived classes.
    //
    virtual  int okToContinue () = 0;
    //
    // Init grid data at problem start-up.
    // This is a pure virtual function and hence MUST be
    // implemented by derived classes.
    //
    virtual void initData () = 0;
    //
    // Set the time levels of state data.
    //
    virtual void setTimeLevel (Real time,Real& dt);

    virtual void init (
      AmrLevel& old,
      const BoxArray& ba_in,
      const DistributionMapping& dmap_in) = 0;

    virtual void init (
      const BoxArray& ba_in,
      const DistributionMapping& dmap_in) = 0;

    int Level () const;
    //
    // List of grids at this level.
    //
    const BoxArray& boxArray () const { return grids; }
    const DistributionMapping& DistributionMap() const { return dmap; }

    //
    // Number of grids at this level.
    //
    int numGrids () const;
    //
    // Returns the indices defining physical domain.
    //
    const Box& Domain () const;
    //
    // Timestep n at this level.
    //
    int nStep () const;
    //
    // Returns the geometry object.
    //
    const Geometry& Geom () const;
    //
    // Returns number of cells on level.
    //
    Long countCells () const;
    //
    // Error estimation for regridding. This is a pure virtual
    // function and hence MUST be implemented by derived classes.
    //
    virtual void errorEst (TagBoxArray& tb,
                           int          clearval,
                           int          tagval,
			   Real         time,
			   int          n_error_buf = 0,
                           int          ngrow = 0) = 0;
    //
    // Interpolate from coarse level to the valid area in dest.
    //

    void FillPatch (AmrLevel & old,
                    MultiFab& mf_to_be_filled,
                    int       dcomp,
                    Real      time,
                    int       index,
                    int       scomp,
                    int       ncomp,
		    int debug_fillpatch);

    void InterpBordersGHOST (
		    Vector<MultiFab*> tower_data,
		    int level_in,
		    int ngrow_in,
                    Real      time,
                    int       index,
                    int       scomp,
                    Vector<int> scompBC_map,
                    int       ncomp,
		    int debug_fillpatch);

    void FillCoarsePatchGHOST (
		    Vector<MultiFab*> tower_data,
		    int level_in,
		    int ngrow_in,
                    Real      time,
                    int       index,
                    int       scomp,
                    Vector<int> scompBC_map,
                    int       ncomp,
		    int debug_fillpatch);

    void InterpBorders (
		    Vector<MultiFab*> tower_data,
		    int level_in,
		    int ngrow_in,
                    Real      time,
                    int       index,
                    int       scomp,
                    Vector<int> scompBC_map,
                    int       ncomp,
		    int debug_fillpatch);

    void FillCoarsePatch (MultiFab& mf,
                          int       dcomp,
                          Real      time,
                          int       index,
                          int       scomp,
                          int       ncomp,
			  int debug_fillpatch);


    //
    // State data object.
    //
    StateData& get_state_data (int state_indx);
    //
    MultiFab& get_new_data (int state_indx,int slab_index);
    const MultiFab& get_new_data (int state_indx,int slab_index) const;

    //
    // Returns list of Descriptors.
    //
    static const DescriptorList& get_desc_lst ();
    static const DescriptorList& get_desc_lstGHOST ();
    //
    // Boundary condition access function.
    //
    Vector<int> getBCArray (int State_Type,
                           int gridno,
                           int scomp,
                           int ncomp);
    Vector<int> getBCArrayGHOST (int State_Type,
                           int gridno,
                           int scomp,
                           int ncomp);
    // 
    // Called in grid_places after other tagging routines to modify
    // the list of tagged points.  Default implementation does nothing.
    //
    virtual void manual_tags_placement (TagBoxArray&    tags,
             const Vector<IntVect>& bf_lev);
    const BoxArray& getAreaNotToTag() noexcept;

     // new_data_FSI[0..level_MAX_NUM_SLAB-1]
    Vector< FSI_container_class > new_data_FSI;

#ifdef AMREX_PARTICLES
    //! This function can be called from the parent
    virtual void particle_redistribute(
      int /*lbase*/ = 0, 
      bool /*a_init*/ = false) {;}

    void CopyNewToOldPC(int lev_max);
    void CopyOldToNewPC(int lev_max);

    AmrParticleContainer<N_EXTRA_REAL,N_EXTRA_INT,0,0>& 
       newDataPC (int slab_index);

    const AmrParticleContainer<N_EXTRA_REAL,N_EXTRA_INT,0,0>&
       newDataPC (int slab_index) const;

     // 
     // ARRAY OF STRUCTURES: x,y,z, particle id, cpu are base 
     // STRUCTURE OF ARRAYS: no base
     //
    using PLS_ParticleContainerType=
        AmrParticleContainer<N_EXTRA_REAL,N_EXTRA_INT,0,0>;

//  Vector<PLS_ParticleContainerType*> AmrLevel0_new_dataPC;
    Vector<std::unique_ptr<PLS_ParticleContainerType> > AmrLevel0_new_dataPC;
#endif

    static dynamic_blobclass_array blob_history_class;

protected:
    //
    // The constructors -- for derived classes.
    //
    AmrLevel () noexcept;

    AmrLevel (AmrCore&        papa,
              int             lev,
              const Geometry& level_geom,
              const BoxArray& bl,
              const DistributionMapping& dm,
              Real            time);
    //
    // Common code used by all constructors.
    //
    void finishConstructor (); 
    //
    // The Data.
    //
    int level;                        // AMR level (0 is coarsest).
    Geometry geom;                    // Geom at this level.
    BoxArray grids;                   // Cell-centered locations of grids.
    DistributionMapping dmap;         // Distribution of grids among processes
    Vector<RealBox> grid_loc;          // Physical locations of each grid.
    AmrCore* parent;                   // Pointer to parent AMR structure.
    static DescriptorList desc_lst;   // List of state variables.
    static DescriptorList desc_lstGHOST; 
    Vector<StateData> state;           // Array of state data.

    BoxArray m_AreaNotToTag;

private:
    //
    // Disallowed.
    //
    AmrLevel (const AmrLevel&);
    AmrLevel& operator = (const AmrLevel&);
};

inline
int
AmrLevel::Level () const
{
    return level;
}

inline
int
AmrLevel::numGrids () const
{
    return grids.size();
}

inline
const Box&
AmrLevel::Domain () const
{
    return geom.Domain();
}

inline
int
AmrLevel::nStep () const
{
    return parent->levelSteps(level);
}

inline
const Geometry&
AmrLevel::Geom () const
{
    return geom;
}

inline
StateData&
AmrLevel::get_state_data (int state_indx)
{
    return state[state_indx];
}

inline
MultiFab&
AmrLevel::get_new_data (int state_indx,int slab_index)
{
    return state[state_indx].newData(slab_index);
}

inline
const MultiFab&
AmrLevel::get_new_data (int state_indx,int slab_index) const
{
    return state[state_indx].newData(slab_index);
}

inline
const DescriptorList&
AmrLevel::get_desc_lst ()
{
    return desc_lst;
}
inline
const DescriptorList&
AmrLevel::get_desc_lstGHOST ()
{
    return desc_lstGHOST;
}

} // namespace amrex

#endif /*AMREX_AmrLevel_H_*/
