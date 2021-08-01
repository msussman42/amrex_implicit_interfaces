
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
#include <AMReX_AmrParticles.H>
#include <AMReX_NeighborParticles.H>

#include <StateDescriptor.H>
#include <Amr.H>

#include <N_EXTRA_REAL.H>

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
    StateData (Amr& papa,
	       int level,
	       int max_level,
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
    void define (Amr& papa,
		 int level,
		 int max_level,
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

    void CopyNewToOld (int level,int max_level);
    void CopyOldToNew (int level,int max_level);


    // Sets time of old and new data.
    // dt will be modified if t_new>0 and t_new-dt<0
    void setTimeLevel (Real t_new,Real& dt);

    static void get_grid_type(IndexType local_typ,int& grid_type);

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
                     std::ostream&      os,
		     int level,
		     int max_level);
    //
    // Initializes state data from a checkpoint file.
    //
    void restart (Amr& papa,
		  int time_order,
                  int slab_dt_type, // 0=SEM 1=evenly spaced
                  int MAX_NUM_SLAB,
                  int level,
                  int max_level,
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

      // lives only on level==0
    AmrParticleContainer<N_EXTRA_REAL,0,0,0>&
       newDataPC (int slab_index,int sub_index);

    //
    // Returns the new data.
    //
    const MultiFab& newData (int slab_index) const;

      // lives only on level==0
    const AmrParticleContainer<N_EXTRA_REAL,0,0,0>&
       newDataPC (int slab_index,int sub_index) const;

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

    int get_bfact_time_order() const;

    void get_time_bounding_box(Real time,Real& nudge_time,
      int& start_index);
    void get_time_index(Real time,Real& nudge_time,
      int& best_index);

private:

    int StateData_level;

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
     // ONLY LEVEL==0 STATEDATA PARTICLES GET INITIALIZED: THEY HOLD
     // PARTICLES THAT APPEAR ON ALL THE LEVELS.
     // 
     // ARRAY OF STRUCTURES: x,y,z, particle id, cpu are base 
     // STRUCTURE OF ARRAYS: no base
     //
     // If ns.viscoelastic_model[im]==2 (elastic material) and
     //    ns.particleLS_flag[im]==1 (particle levelset==0.0 for
     //    interface particles),
     //    then define "new_dataPC"; 
     //
    using PLS_ParticleContainerType=
	    AmrParticleContainer<N_EXTRA_REAL,0,0,0>;

     // lives only on level==0
    Vector< 
	Vector<PLS_ParticleContainerType*> 
	  > new_dataPC;
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

//PhysBCFunctBaseSussman declared in:
// amrex_implicit_interfaces/amrex-master/Src/Base/
//    PhysBCFunctSUSSMAN.[H|cpp]
//
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
