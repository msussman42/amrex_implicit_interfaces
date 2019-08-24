/*
** (c) 1996-2000 The Regents of the University of California (through
** E.O. Lawrence Berkeley National Laboratory), subject to approval by
** the U.S. Department of Energy.  Your use of this software is under
** license -- the license agreement is attached and included in the
** directory as license.txt or you may contact Berkeley Lab's Technology
** Transfer Department at TTD@lbl.gov.  NOTICE OF U.S. GOVERNMENT RIGHTS.
** The Software was developed under funding from the U.S. Government
** which consequently retains certain rights as follows: the
** U.S. Government has been granted for itself and others acting on its
** behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
** Software to reproduce, prepare derivative works, and perform publicly
** and display publicly.  Beginning five (5) years after the date
** permission to assert copyright is obtained from the U.S. Department of
** Energy, and subject to any subsequent five (5) year renewals, the
** U.S. Government is granted for itself and others acting on its behalf
** a paid-up, nonexclusive, irrevocable, worldwide license in the
** Software to reproduce, prepare derivative works, distribute copies to
** the public, perform publicly and display publicly, and to permit
** others to do so.
*/

#ifndef _AmrLevel_H_
#define _AmrLevel_H_ 
//
// $Id: AmrLevel.H,v 1.52 2002/03/27 17:10:20 lijewski Exp $
//
#include <REAL.H>
#include <Box.H>
#include <RealBox.H>
#include <Geometry.H>
#include <FArrayBox.H>
#include <BoxArray.H>
#include <Array.H>
#include <FabArray.H>
#include <MultiFab.H>
#include <Derive.H>
#include <BCRec.H>
#include <Amr.H>
#include <StateDescriptor.H>
#include <StateData.H>
#include <VisMF.H>

class TagBox;
class TagBoxArray;

//
//@Man:
//@Memo: Virtual base class for managing individual levels.
/*@Doc:

  AmrLevel functions both as a container for state data on a level
  and also manages the advancement of data in time.
*/

class AmrLevel
{

public:

    //
    //@ManDoc: What time are we at?
    //
    enum TimeLevel { AmrOldTime, AmrHalfTime, AmrNewTime, AmrOtherTime };
    //
    //@ManDoc: The destructor.
    //
    virtual ~AmrLevel ();

    virtual std::string thePlotFileType () const = 0;

    virtual void writePlotFile (const std::string& dir,
      std::ostream& os,VisMF::How how = VisMF::OneFilePerCPU) = 0;

    //
    //@ManDoc: Write current state to checkpoint file.
    //
    virtual void checkPoint (const std::string& dir,
                             std::ostream&  os,
                             VisMF::How     how = VisMF::OneFilePerCPU);
    //
    //@ManDoc: Restart from a checkpoint file.
    //
    virtual void restart (Amr&          papa,
                          std::istream& is,
			  bool          bReadSpecial = false);
    //
    //@ManDoc: Is name a state variable?
    //
    static int isStateVariable (const std::string& name,
                                int&               state_indx,
                                int&               ncomp);

    /*@ManDoc: Compute the initial time step.  This is a pure virtual function
               and hence MUST be implemented by derived classes.
    */
    virtual void computeInitialDt (Array<Real>&          dt_level,
                                   Real                  stop_time) = 0;

    /*@ManDoc: Compute the next time step.  This is a pure virtual function
               and hence MUST be implemented by derived classes.
    */

    virtual void computeNewDt (Array<Real>&          dt_level,
                               Real                  stop_time) = 0;

    virtual void StableFillPatch(MultiFab* oldfine_data,MultiFab& S_new,
     int ngrow,int is_old,int scomp,int dcomp,int ncomp,
     const BoxArray& finegrids_coarse) = 0;
    virtual MultiFab* LinInterpData(int is_old,int scomp,int ncomp) = 0;



    /*@ManDoc: Do an integration step on this level.  Returns maximum safe
               time step.  This is a pure virtual function and hence MUST
               be implemented by derived classes.
    */
    virtual Real advance (Real time,
                          Real dt,
                          int  iteration,
                          int  ncycle,
			  int  sub_cycle) = 0;

    /*@ManDoc: Contains operations to be done after a timestep.  This is a
               pure virtual function and hence MUST be implemented by derived
               classes.
    */
    virtual  void post_timestep (int iteration,int sub_cycle,
     Real stop_time) = 0;

    /*@ManDoc: Contains operations to be done only after a full coarse
               timestep.  The default implementation does nothing.
    */
    virtual void postCoarseTimeStep (Real time);

    /*@ManDoc: Operations to be done after restart.  This is a pure virtual
               function and hence MUST be implemented by derived classes.
    */
    virtual  void post_restart () = 0;

    /*@ManDoc: Operations to be done after initialization.
               This is a pure virtual function and hence MUST be
               implemented by derived classes.
    */
    virtual  void post_init (Real stop_time) = 0;

    /*@ManDoc: Is it ok to continue the calculation?
               This is a pure virtual function and hence MUST be
               implemented by derived classes.
    */
    virtual  int okToContinue () = 0;

    /*@ManDoc: Should I regrid with this level as base level?
               This test is only evaluated when level\_count >= regrid\_int
               as well. Defaults to true.
    */
    virtual  int okToRegrid ();

    /*@ManDoc: Init grid data at problem start-up.
               This is a pure virtual function and hence MUST be
               implemented by derived classes.
    */
    virtual void initData () = 0;
    //
    //@ManDoc: Set the time levels of state data.
    //
    virtual void setTimeLevel (Real time,Real dt);
    //

    /*@ManDoc: Init data on this level from another AmrLevel (during regrid).
               This is a pure virtual function and hence MUST be
               implemented by derived classes.
    */
    virtual void init (AmrLevel &old) = 0;

    /*@ManDoc: Init data on this level after regridding if old AmrLevel
               did not previously exist. This is a pure virtual function
               and hence MUST be implemented by derived classes.
    */
    virtual void init () = 0;
    //
    //@ManDoc: Returns this AmrLevel.
    //
    int Level () const;

    int get_fine_ratio(int dir) const;

    //
    //@ManDoc: List of grids at this level.
    //
    const BoxArray& boxArray () const;
    //
    //@ManDoc: Number of grids at this level.
    //
    int numGrids () const;
    //
    //@ManDoc: Physical locations of grids at this level.
    //
    const Array<RealBox>& gridLocations () const;
    //
    //@ManDoc: Returns the indices defining physical domain.
    //
    const Box& Domain () const;
    //
    //@ManDoc: Timestep n at this level.
    //
    long nStep () const;
    //
    //@ManDoc: Returns the geometry object.
    //
    const Geometry& Geom () const;
    //
    //@ManDoc: Returns number of cells on level.
    //
    long countCells () const;

    /*@ManDoc: Error estimation for regridding. This is a pure virtual
               function and hence MUST be implemented by derived classes.
    */
    virtual void errorEst (TagBoxArray& tb, int clearval, int tagval,
                           int is_old)=0;

    /*@ManDoc: Returns a MultiFab containing the derived data for this level.
               The user is responsible for deleting this pointer when done
               with it.  If ngrow>0 the MultiFab is built on the appropriately
               grown BoxArray.
    */
    int derive_ncomp(const std::string& name);

    MultiFab* derive (const std::string& name,int is_old,int ngrow);

    /*@ManDoc: This version of derive() fills the dcomp'th component of mf
               with the derived quantity.
    */
    void derive (const std::string& name,int is_old,MultiFab& mf,
                 int dcomp);
    //
    //@ManDoc: State data object.
    //
    StateData& get_state_data (int state_indx);
    //
    //@ManDoc: State data at old time.
    //
    MultiFab& get_old_data (int state_indx);
    //
    //@ManDoc: State data at old time.
    //
    const MultiFab& get_old_data (int state_indx) const;
    //
    //@ManDoc: State data at new time.
    //
    MultiFab& get_new_data (int state_indx);
    //
    //@ManDoc: State data at new time.
    //
    const MultiFab& get_new_data (int state_indx) const;
    //
    //@ManDoc: Returns list of Descriptors.
    //
    static const DescriptorList& get_desc_lst ();
    //
    //@ManDoc: Returns list of derived variables.
    //
    static DeriveList& get_derive_lst ();
    //
    //@ManDoc: Boundary condition access function.
    //
    Array<int> getBCArray (int state_idx,
                           int gridno,
                           int scomp,
                           int ncomp);
    virtual void setPlotVariables ();
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
    Array<RealBox> grid_loc;          // Physical locations of each grid.
    Amr* parent;                      // Pointer to parent AMR structure.
    IntVect crse_ratio;               // Refinement ratio to coarser level.
    IntVect fine_ratio;               // Refinement ratio to finer level.
    static DeriveList derive_lst;     // List of derived quantities.
    static DescriptorList desc_lst;   // List of state variables.
    Array<StateData> state;           // Array of state data.

private:
    //
    // Disallowed.
    //
    AmrLevel (const AmrLevel&);
    AmrLevel& operator = (const AmrLevel&);
};


#endif /*_AmrLevel_H_*/
