
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
#include <Derive.H>
#include <Interpolater.H>
#include <Amr.H>
#include <StateDescriptor.H>
#include <StateData.H>

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

class AmrLevel
{

public:
    //
    // What time are we at?
    //
    enum TimeLevel { AmrOldTime,
                     AmrHalfTime,
                     AmrNewTime,
                     Amr1QtrTime,
                     Amr3QtrTime,
                     AmrOtherTime };
    //
    // The destructor.
    //
    virtual ~AmrLevel ();
    //
    // A string written as the first item in writePlotFile() at
    // level zero. This MUST be defined by each derived class.
    // It is so we can distinguish between different types of
    // plot files.  This is a pure virtual function and hence MUST
    // be implemented by derived classes.
    //
    virtual std::string thePlotFileType () const = 0;
    //
    // Write plot file stuff to specified directory.  This is a
    // pure virtual function and hence MUST be implemented by
    // derived classes.
    //
    virtual void writePlotFile (
     const std::string& dir,
     std::ostream& os,
     int do_plot,int do_slice,
     int SDC_outer_sweeps_in,int slab_step_in) = 0;

    //
    // Write current state to checkpoint file.
    //
    virtual void checkPoint (const std::string& dir,
                             std::ostream&      os);
    //
    // Restart from a checkpoint file.
    //
    virtual void restart (Amr&          papa,
                          std::istream& is);
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
       Real& dt,Real stop_time,int post_regrid_flag)=0;
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
    virtual  void post_restart () = 0;

    virtual void post_regrid (int lbase,
                              int new_finest,
                              Real time) = 0;
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

    virtual void init (AmrLevel &old,
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
    long countCells () const;
    //
    // Error estimation for regridding. This is a pure virtual
    // function and hence MUST be implemented by derived classes.
    //
    virtual void errorEst (TagBoxArray& tb,
                           int          clearval,
                           int          tagval,
			   int          n_error_buf = 0,
                           int          ngrow = 0) = 0;
    //
    // Interpolate from coarse level to the valid area in dest.
    //

    void FillPatch (AmrLevel & old,
                    MultiFab& mf,
                    int       dcomp,
                    Real      time,
                    int       index,
                    int       scomp,
                    int       ncomp);

    void InterpBordersGHOST (
                    MultiFab& cmf,
                    MultiFab& mf,
                    Real      time,
                    int       index,
                    int       scomp,
                    Vector<int> scompBC_map,
                    int       ncomp);

    void FillCoarsePatchGHOST (
                    MultiFab& cmf,
                    MultiFab& mf,
                    Real      time,
                    int       index,
                    int       scomp,
                    Vector<int> scompBC_map,
                    int       ncomp);

    void InterpBorders (
                    MultiFab& cmf,
                    MultiFab& mf,
                    Real      time,
                    int       index,
                    int       scomp,
                    Vector<int> scompBC_map,
                    int       ncomp);

    void FillCoarsePatch (MultiFab& mf,
                          int       dcomp,
                          Real      time,
                          int       index,
                          int       scomp,
                          int       ncomp);

    virtual std::unique_ptr<MultiFab> derive (const std::string& name,
                              int ngrow);

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
    // Returns list of derived variables.
    //
    static DeriveList& get_derive_lst ();
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
                                        Vector<int>& bf_lev);
    //
    // Modify list of variables to be plotted
    // 
    virtual void setPlotVariables ();


#ifdef AMREX_PARTICLES
    //! This function can be called from the parent 
    virtual void particle_redistribute (int lbase = 0, bool init = false) {;}
#endif

    //
protected:
    //
    // The constructors -- for derived classes.
    //
    AmrLevel ();

    AmrLevel (Amr&            papa,
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
    Amr* parent;                      // Pointer to parent AMR structure.
    static DeriveList derive_lst;     // List of derived quantities.
    static DescriptorList desc_lst;   // List of state variables.
    static DescriptorList desc_lstGHOST; 
    Vector<StateData> state;           // Array of state data.

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
