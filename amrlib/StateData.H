
#ifndef AMREX_StateData_H_
#define AMREX_StateData_H_ 

#include <memory>

#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MFCopyDescriptor.H>
#include <AMReX_BCRec.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_VisMF.H>
#include <AMReX_DistributionMapping.H>
#include <PhysBCFunctSUSSMAN.H>
#include <AMReX_Geometry.H>
#include <AMReX_RealBox.H>
#include <AMReX_Particles.H>

#include <StateDescriptor.H>

// x,y,z is not extra; original (x,y,z), time added, #steps, are extra.
#define N_EXTRA_REAL_BDRY (AMREX_SPACEDIM+2)
#define N_EXTRA_REAL_BULK (AMREX_SPACEDIM+2)
// x,y,z is not extra; mass, volume, energy
#define N_EXTRA_REAL_SPRAY (3)
#define N_EXTRA_SOA_ijk (AMREX_SPACEDIM)
namespace amrex {

//
// Current and previous level-time data.
//
// StateData holds state data on a level for the current and previous time step.
//

class StateData
{
public:
    //
    // The default constructor.
    //
    StateData ();
    //
    // Constructor that properly initializes data members.
    //
    StateData (int level,
               const Box& p_domain,
               const BoxArray&        grds,
               const DistributionMapping& dm,
               const StateDescriptor* d,
               const StateDescriptor* dGHOST,
               Real                   cur_time,
               Real                   dt,
               int                    time_order,
               int                    slab_dt_type, // 0=SEM 1=evenly spaced
               int                    MAX_NUM_SLAB);
    //
    // The destructor.
    //
    ~StateData ();

    //
    // Initializes data members if you used default constructor.
    //
    void define (int level,
                 const Box& p_domain,
                 const BoxArray&        grds,
                 const DistributionMapping& dm,
                 const StateDescriptor& d, 
                 const StateDescriptor& dGHOST, 
                 Real                   cur_time,
                 Real                   dt,
                 int                    time_order,
                 int                    slab_dt_type, // 0=SEM 1=evenly spaced
                 int                    MAX_NUM_SLAB);

    //
    // Sets up the boundary conditions.
    //
    void buildBC ();
    //
    //

    void CopyNewToOld ();
    void CopyOldToNew ();


    // Sets time of old and new data.
    // dt will be modified if t_new>0 and t_new-dt<0
    void setTimeLevel (Real t_new,Real& dt);

    //
    // Set physical bndry values
    //
    void FillBoundary (int level,
                       FArrayBox& dest,
                       Real       time,
                       const Real*    dx,
                       const RealBox& prob_domain,
                       int            dcomp,
                       Vector<int>     scompBC_map,
                       int            ncomp,
                       int            bfact);


    void FillBoundaryGHOST (int level,
                       FArrayBox& dest,
                       Real       time,
                       const Real*    dx,
                       const RealBox& prob_domain,
                       int            dcomp,
                       Vector<int>     scompBC_map,
                       int            ncomp,
                       int            bfact);

    //
    // Write the state data to a checkpoint file.
    //
    void checkPoint (const std::string& name,
                     const std::string& fullpathname,
                     std::ostream&      os);
    //
    // Initializes state data from a checkpoint file.
    //
    void restart (int time_order,
                  int slab_dt_type, // 0=SEM 1=evenly spaced
                  int MAX_NUM_SLAB,
                  int level,
                  std::istream& is,
                  const Box&             p_domain,
                  const BoxArray&        grds,
                  const DistributionMapping& dm,
                  const StateDescriptor& d,
                  const StateDescriptor& dGHOST,
                  const std::string&     restart_file);
    //
    // Returns the StateDescriptor.
    //
    const StateDescriptor* descriptor () const;
    const StateDescriptor* descriptorGHOST () const;
    //
    // Returns the valid domain.
    //
    const Box& getDomain () const;

    const BoxArray& boxArray () const { return grids; }
    const DistributionMapping& DistributionMap () const { return dmap; }

    Real slabTime (int slab_index) const;
    //
    // Returns the new data.
    //
    MultiFab& newData (int slab_index);
    //
    // Returns the new data.
    //
    const MultiFab& newData (int slab_index) const;
    //
    // Returns an array of boundary condition records for the 
    // specified component.  Each element in the array corresponds
    // to a grid in the state.
    //
    Vector<BCRec>& getBCs (int comp);
    Vector<BCRec>& getBCsGHOST (int comp);
    //
    // Returns boundary conditions of specified component on the
    //         specified grid.
    //
    const BCRec& getBC (int comp, int i) const;
    const BCRec& getBCGHOST (int comp, int i) const;
    //
    // Prints out the time interval.
    //
    void printTimeInterval (std::ostream& os) const;
    //
    // True if there is any new data available.
    //
    bool hasNewData (int slab_index) const;

    int get_bfact_time_order() const;

    void get_time_bounding_box(Real time,Real& nudge_time,
      int& start_index);
    void get_time_index(Real time,Real& nudge_time,
      int& best_index);

private:

    int StateData_MAX_NUM_SLAB;
    int StateData_slab_dt_type;

    int bfact_time_order;

    //
    // Pointer to data descriptor.
    //
    const StateDescriptor* desc;
    const StateDescriptor* descGHOST;
    //
    // Problem domain.
    //
    Box domain;
    //
    // Grids defined at this level.
    //
    BoxArray grids;
    DistributionMapping dmap;

    // MAX_NUM_SLAB entries.
    Vector<Real> time_array;
    //
    // Pointer to new-time data.
    // MAX_NUM_SLAB entries.
    Vector< MultiFab* > new_data;
     // ARRAY OF STRUCTURES: x,y,z, particle id, cpu are base 
     // STRUCTURE OF ARRAYS: no base
     // for each material, there are two particle containers: (1) interface
     // particles, (2) interior particles.
     // interface particles:
     // x,y,z is not extra; original (x,y,z) and time added are extra.
     // interior particles:
     // x,y,z is not extra; original (x,y,z) and time added are extra.
     // last two parameters are for structure of arrays.
     // NStructReal=number of extra Real variables 
     // (not including particle position)
     //
     // If ns.viscoelastic_model[im]==2 (elastic material) and
     //    ns.particleLS_flag[im]==2 (particles inside and on interface),
     //    then define "new_dataPC_BDRY" and "new_dataPC_BULK"
    Vector< ParticleContainer<N_EXTRA_REAL_BDRY,0,0,0 >* > new_dataPC_BDRY;
    Vector< ParticleContainer<N_EXTRA_REAL_BULK,0,0,0 >* > new_dataPC_BULK;
    Vector< ParticleContainer<N_EXTRA_REAL_SPRAY,0,0,0 >* > new_dataPC_SPRAY;
    Vector< ParticleContainer<0,0,0,N_EXTRA_SOA_ijk >* > 
	    structure_of_arrayLS_inner;
    Vector< ParticleContainer<0,0,0,N_EXTRA_SOA_ijk >* > 
	    structure_of_arrayLS_outer;
    //
    // Boundary conditions for the grids at this level.
    //
    Vector< Vector<BCRec> > bc;
    Vector< Vector<BCRec> > bcGHOST;

}; // StateData

class StateDataPhysBCFunct
    : public PhysBCFunctBaseSUSSMAN
{
public:

    StateDataPhysBCFunct (
      StateData& sd, 
      const Geometry& geom_);

    virtual void FillBoundary (
      int level,
      MultiFab& mf, 
      Real time,
      int dcomp, 
      Vector<int> scompBC_map, 
      int ncomp, 
      int bfact);

private:

    StateData* statedata;
    const Geometry& geom;

}; // StateDataPhysBCFunct


class StateDataPhysBCFunctGHOST
    : public PhysBCFunctBaseSUSSMAN
{
public:

    StateDataPhysBCFunctGHOST (
      StateData& sd, 
      const Geometry& geom_);

    virtual void FillBoundary (
      int level,
      MultiFab& mf, 
      Real time,
      int dcomp, 
      Vector<int> scompBC_map, 
      int ncomp, 
      int bfact);

private:

    StateData* statedata;
    const Geometry& geom;

}; // StateDataPhysBCFunctGHOST

} //namespace amrex

#endif /*AMREX_StateData_H_*/
