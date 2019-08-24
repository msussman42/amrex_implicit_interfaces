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

#ifndef _Amr_H_
#define _Amr_H_ 

//
// $Id: Amr.H,v 1.30 2001/08/01 21:50:44 lijewski Exp $
//

#include <fstream>

#include <Box.H>
#include <Geometry.H>
#include <BoxArray.H>
#include <Array.H>
#include <PArray.H>

class AmrLevel;
class LevelBld;

//
//@Man:
//@Memo: Manages hierarchy of levels for time-dependent AMR computations.
/*@Doc:

  The Amr class is designed to manage parts of the computation  which do
  not belong on a single level, like establishing and updating the hierarchy
  of levels, global timestepping, and managing the different AmrLevels
*/

class Amr
{
public:

    int is_first_step_after_restart;

    //
    //@ManDoc: The constructor.
    //
    Amr ();
    //
    //@ManDoc: The destructor.
    //
    ~Amr ();
    //
    //@ManDoc: Init data after construction. Must be called before timestepping.
    //
    void init (Real strt_time, Real stop_time);
    //
    //@ManDoc: Set the timestep on each level.
    //
    void setDtLevel (const Array<Real>& dt_lev);
    void set_ncycle_finest(int new_ncycle_finest);
    //
    //@ManDoc: Largest allowable grid.
    //
    int maxGridSize () const;
    //
    //@ManDoc: Subcycle in time?
    //
    int subCycle () const;
    //
    //@ManDoc: Maximum allowable level.
    //
    int maxLevel () const;
    //
    //@ManDoc: Finest currently existing level.
    //
    int finestLevel () const;
    //
    //@ManDoc: Physical time.
    //
    Real cumTime () const;
    //
    //@ManDoc: Time step at specified level.
    //
    Real dtLevel (int level) const;
    //
    //@ManDoc: Array of time steps at all levels.
    //
    const Array<Real>& dtLevel () const;
    //
    //@ManDoc: Geometry object.
    //
    const Geometry& Geom (int level) const;
    //
    //@ManDoc: Number of subcycled time steps.
    //
    int nCycle (int level) const;
    int nCycleSave (int level) const;
    int bFactor(int level) const;
    //
    //@ManDoc: Refinement ratio.
    //
    IntVect refRatio (int coarse_level) const;
    //
    //@ManDoc: Maximum refinement ratio in any direction.
    //
    int MaxRefRatio (int level) const;
    //
    //@ManDoc: Refinement ratios between all levels.
    //
    const Array<IntVect>& refRatio () const;
    //
    //@ManDoc: Number of time steps at specified level.
    //
    long levelSteps (int lev) const;
    //
    //@ManDoc: Interval between regridding.
    //
    long regridInt (int lev) const;
    //
    //@ManDoc: Number of time steps between checkpoint files.
    //
    int checkInt () const;
    //
    //@ManDoc: Time between checkpoint files.
    //
    Real checkPer() const;
    //
    //@ManDoc: Number of time steps between plot files.
    //
    int plotInt () const;
    //
    //@ManDoc: Time between plot files.
    //
    Real plotPer () const;

    void output_checks(Real cumtime,int& check_test,int& plot_test,
      long& levelstepszero,int endofstep);

    static const std::list<std::string>& statePlotVars ();

    static bool isStatePlotVar (const std::string& name);

    static void addStatePlotVar (const std::string& name);

    static void deleteStatePlotVar (const std::string& name);

    static void clearStatePlotVarList ();

    static void fillStatePlotVarList ();

    static const std::list<std::string>& derivePlotVars ();

    static bool isDerivePlotVar (const std::string& name);

    static void addDerivePlotVar (const std::string& name);

    static void deleteDerivePlotVar (const std::string& name);

    static void clearDerivePlotVarList ();

    static void fillDerivePlotVarList ();
    //
    //@ManDoc: Number of buffer cells in error estimator.
    //
    int nErrorBuf (int lev) const;
    //
    //@ManDoc: Minimum allowable grid efficiency.
    //
    Real gridEff () const;
    //
    //@ManDoc: Number of cells to define proper nesting 
    //
    int nProper () const;
    //
    //@ManDoc: List of grids at each level.
    //
    const BoxArray& boxArray (int lev) const;
    //
    //@ManDoc: AmrLevel lev.
    //
    AmrLevel&  getLevel (int lev);
    //
    //@ManDoc: PArray of AmrLevels.
    //
    PArray<AmrLevel>& getAmrLevels ();
    //
    //@ManDoc: Total number of cells.
    //
    long cellCount ();
    //
    //@ManDoc: Number of cells at given level.
    //
    long cellCount (int lev);
    //
    //@ManDoc: Total number of grids.
    //
    int numGrids ();
    //
    //@ManDoc: Number of grids at given level.
    //
    int numGrids (int lev);
    //
    //@ManDoc: More work to be done?
    //
    virtual int okToContinue ();
    //
    //@ManDoc: Do a complete integration cycle.
    //
    void coarseTimeStep (Real stop_time);
    //
    //@ManDoc: Retrieve derived data. User is responsible for deleting pointer.
    //
    MultiFab* derive (const std::string& name,
      int is_old,int lev,int ngrow);

    //
    //@ManDoc: Name of the restart file.
    //
    const std::string& theRestartFile () const;
    //
    //@ManDoc: The ith datalog file.  Do with it what you want.
    //
    std::ostream& DataLog (int i);

protected:

    //
    // Used by grid_places().
    //
    static void ProjPeriodic (BoxDomain&      bd,
                              const Geometry& geom);
    //
    // Write the plot file to be used for visualization.
    //
    void writePlotFile (const std::string& root, long num);
    //
    // Initialize grid hierarchy.
    //
    void initialInit (Real strt_time, Real stop_time);
    //
    // Check for valid input.
    //
    void checkInput ();
    //
    // Write current state into a chk* file.
    //
    void checkPoint ();
    //
    // Restart from a checkpoint file.
    //
    void restart (const std::string& filename);
    //
    // Define and initialize coarsest level.
    //
    void defBaseLevel (Real start_time);
    //
    // Define and initialize refined levels.
    //
    void bldFineLevels (Real start_time);
    //
    // Rebuild grid hierarchy finer than lbase.
    //
    void regrid (int  lbase,
                 Real time);
    //
    // Define new grid locations (called from regrid) and put into new_grids.
    //
    void grid_places (int              lbase,
                      Real             time,
                      int&             new_finest,
                      Array<BoxArray>& new_grids);
    //
    // Do a single timestep on level L.
    //
    void timeStep (int  level,
                   Real time,
                   int  iteration,
                   int  ncycle,
                   Real prev_dt,Real stop_time);

    void printGridInfo (std::ostream& os,
                        int           min_lev,
                        int           max_lev);

    void setRecordGridInfo (const std::string&);

    void setRecordRunInfo (const std::string&);

    void setRecordDataInfo (int i, const std::string&);
    //
    // The data ...
    //
    std::string      grids_file;   // Grids file that will bypass regridding.
    int              max_level;    // Maximum allowed level.
    int              finest_level; // Current finest level.
    Array<IntVect>   ref_ratio;    // Refinement ratios [0:finest_level-1]
    PArray<AmrLevel> amr_level;    // Array of levels
    Array<Geometry>  geom;         // Index extent of domain at all levels.
    Real             cumtime;      // Physical time variable.
    Array<Real>      dt_level;     // Timestep at this level.
    Array<long>      level_steps;  // Number of time steps at this level.
    Array<long>      level_count;
    Array<long>      level_cells;
    Array<int>       n_cycle;
    Array<int>       n_cycle_save;
    bool             isPeriodic[BL_SPACEDIM];  // Domain periodic?
    Array<long>      regrid_int;      // Interval between regridding.
    Array<int>       n_error_buf;     // Buffer cells around each tagged cell.
    Array<int>       blocking_factor; // Blocking factor in grid generation.
    Real             grid_eff;        // Grid efficiency.
    int              max_grid_size;   // Maximum allowable grid size.
    int              n_proper;        // # cells required for proper nesting.
    long             last_checkpoint; // Step number of previous checkpoint.
    int              check_int;       // How often checkpoint (# time steps).
    Real             check_per;       // How often checkpoint (units of time).
    std::string      check_file_root; // Root name of checkpoint file.
    long             last_plotfile;   // Step number of previous plotfile.
    int              plot_int;        // How often plotfile (# of time steps)
    Real             plot_per;        // How often plotfile (in units of time)
    std::string      plot_file_root;  // Root name of plotfile.
    int              verbose;
    int              record_grid_info;
    int              record_run_info;
    std::ofstream    gridlog;
    std::ofstream    runlog;
    PArray<std::ofstream> datalog;
    int              sub_cycle;
    std::string      restart_file;
    std::string      probin_file;
    LevelBld*        levelbld;
    //
    // The static data ...
    //
    static std::list<std::string> state_plot_vars;  // State Vars to dump to plotfile 
    static std::list<std::string> derive_plot_vars; // Derived Vars to dump to plotfile 
    static bool                   first_plotfile;

private:
    //
    // Disallowed
    //
    Amr (const Amr& rhs);
    Amr& operator= (const Amr& rhs);
};

#endif /*_Amr_H_*/
