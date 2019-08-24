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

#ifndef _StateData_H_
#define _StateData_H_ 
//
// $Id: StateData.H,v 1.19 2001/08/01 21:50:46 lijewski Exp $
//
#include <Box.H>
#include <BoxArray.H>
#include <MultiFab.H>
#include <BCRec.H>
#include <Array.H>
#include <VisMF.H>

class RealBox;
class StateDescriptor;

namespace BoxLib
{
    void readBoxArray (BoxArray&, std::istream&, bool);
}
//
//@Man:
//@Memo: Current and previous level-time data.
/*@Doc:

  StateData holds state data on a level for the current and previous time step.
*/

class StateData
{
public:
    //
    //@ManDoc: The default constructor.
    //
    StateData ();
    //
    //@ManDoc: Constructor that properly initializes data members.
    //
    StateData (const Box&             p_domain,
               const BoxArray&        grds,
               const StateDescriptor* d,
               Real                   cur_time,
               Real                   dt);
    //
    //@ManDoc: The destructor.
    //
    ~StateData ();
    //
    //@ManDoc: Initializes data members if you used default constructor.
    //
    void define (const Box&             p_domain,
                 const BoxArray&        grds,
                 const StateDescriptor& d, 
                 Real                   cur_time,
                 Real                   dt);
    //
    //@ManDoc: Sets up the boundary conditions.
    //
    void buildBC ();
    //
    //@ManDoc: Old data becomes new data and new time is incremented by dt.
    //
    void CopyNewToOld (Real dt);
    //
    //@ManDoc: Sets time of old and new data.
    //
    void setTimeLevel (Real t_new,Real dt);
    //
    //@ManDoc: Sets time of old data.
    //
    void setOldTimeLevel (Real t_old);
    //
    //@ManDoc: Sets time of new data.
    //
    void setNewTimeLevel (Real t_new);

    void checkPoint (const std::string& name,
                     const std::string& fullpathname,
                     std::ostream&      os,
                     VisMF::How         how);
    //
    //@ManDoc: Initializes state data from a checkpoint file.
    //
    void restart (std::istream&          is,
                  const StateDescriptor& d,
                  const std::string&     restart_file,
		  bool                   bReadSpecial = false);
    //
    //@ManDoc: Returns the StateDescriptor.
    //
    const StateDescriptor* descriptor () const;
    //
    //@ManDoc: Returns the valid domain.
    //
    const Box& getDomain () const;

    //
    //@ManDoc: Returns the BoxArray.
    //
    const BoxArray& boxArray () const;
    //
    //@ManDoc: Returns the current time.
    //
    Real curTime () const;
    //
    //@ManDoc: Returns the previous time.
    //
    Real prevTime () const;
    //
    //@ManDoc: Returns the new data.
    //
    MultiFab& newData ();
    //
    //@ManDoc: Returns the new data.
    //
    const MultiFab& newData () const;
    //
    //@ManDoc: Returns the old data.
    //
    MultiFab& oldData ();
    //
    //@ManDoc: Returns the old data.
    //
    const MultiFab& oldData () const;
    //
    //@ManDoc: Returns the FAB of new data at grid index `i'.
    //
    FArrayBox& newGrid (int i);
    //
    //@ManDoc: Returns the FAB of old data at grid index `i'.
    //
    FArrayBox& oldGrid (int i);
    //
    //@ManDoc: Returns an array of boundary condition records for the 
    //         specified component.  Each element in the array corresponds
    //         to a grid in the state.
    //
    Array<BCRec>& getBCs (int comp);
    //
    //@ManDoc: Returns boundary conditions of specified component on the
    //         specified grid.
    //
    const BCRec& getBC (int comp, int i) const;
    //
    //@ManDoc: Prints out the time interval.
    //
    void printTimeInterval (std::ostream& os) const;
    //
    bool hasData () const;
    //

protected:

    struct TimeInterval
    {
        Real start;
        Real stop;
    };
    //
    // Pointer to data descriptor.
    //
    const StateDescriptor* desc;
    //
    // Problem domain.
    //
    Box domain;
    //
    // Grids defined at this level.
    //
    BoxArray grids;
    //
    // Time variable assoc with new data.
    //
    TimeInterval new_time;
    //
    // Time variable assoc with old data.
    //
    TimeInterval old_time;
    //
    // Pointer to new-time data.
    //
    MultiFab* new_data;
    //
    // Pointer to previous time data.
    //
    MultiFab* old_data;
    //
    // Boundary conditions for the grids at this level.
    //
    Array< Array<BCRec> > bc;
};

#endif /*_StateData_H_*/
