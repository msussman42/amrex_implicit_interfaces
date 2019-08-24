
#ifndef _StateData_H_
#define _StateData_H_ 

#include <Box.H>
#include <BoxArray.H>
#include <MultiFab.H>
#include <BCRec.H>
#include <Array.H>
#include <VisMF.H>
#include <DistributionMapping.H>
#include <PhysBCFunct.H>
#include <Geometry.H>
#include <RealBox.H>
#include <StateDescriptor.H>

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
                       Array<int>     scompBC_map,
                       int            ncomp,
                       int            bfact);


    void FillBoundaryGHOST (int level,
                       FArrayBox& dest,
                       Real       time,
                       const Real*    dx,
                       const RealBox& prob_domain,
                       int            dcomp,
                       Array<int>     scompBC_map,
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
    Array<BCRec>& getBCs (int comp);
    Array<BCRec>& getBCsGHOST (int comp);
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
    Array<Real> time_array;
    //
    // Pointer to new-time data.
    // MAX_NUM_SLAB entries.
    Array<MultiFab*> new_data;
    //
    // Boundary conditions for the grids at this level.
    //
    Array< Array<BCRec> > bc;
    Array< Array<BCRec> > bcGHOST;

}; // StateData

class StateDataPhysBCFunct
    : public PhysBCFunctBase
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
      Array<int> scompBC_map, 
      int ncomp, 
      int bfact);

private:

    StateData* statedata;
    const Geometry& geom;

}; // StateDataPhysBCFunct


class StateDataPhysBCFunctGHOST
    : public PhysBCFunctBase
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
      Array<int> scompBC_map, 
      int ncomp, 
      int bfact);

private:

    StateData* statedata;
    const Geometry& geom;

}; // StateDataPhysBCFunctGHOST



#endif /*_StateData_H_*/
