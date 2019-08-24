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

#ifndef _StateDescriptor_H_
#define _StateDescriptor_H_ 
//
// $Id: StateDescriptor.H,v 1.19 2001/08/01 21:50:46 lijewski Exp $
//
#include <utility>
#include <vector>

#include <Array.H>
#include <PArray.H>
#include <REAL.H>
#include <Box.H>
#include <BCRec.H>
#include <ArrayLim.H>

//
//@Man:
//@Memo: Attributes of StateData.
/*@Doc:

  Attributes of StateData.
*/

class StateDescriptor
{
    friend class DescriptorList;

public:
    //
    //@ManDoc: Defines the temporal centeredness of data.
    //
    enum TimeCenter { Point = 0 , Interval };
    //
    //@ManDoc: The default constructor.
    //
    StateDescriptor ();
    //
    //@ManDoc: Constructor that sets all data members.
    //
    StateDescriptor (IndexType     btyp,
                     TimeCenter    ttyp,
                     int           ident,
                     int           nextra,
                     int           num_comp);
    //
    //@ManDoc: The destructor.
    //
    ~StateDescriptor ();
    //
    //@ManDoc: Define the data members if constructed with default construtor.
    //
    void define (IndexType     btyp,
                 TimeCenter    ttyp,
                 int           ident,
                 int           nextra,
                 int           num_comp);
    //
    //@ManDoc: Sets details of `comp' component.
    //
    void setComponent (int comp,
                       const std::string&   nm,
                       const BCRec& bc,
                       int interp_type);

    //
    //@ManDoc: Output names of components.
    //
    void dumpNames (std::ostream& os,
                    int           start_comp,
                    int           num_comp) const;
    //
    //@ManDoc: Returns the IndexType.
    //
    IndexType getType () const;
    //
    //@ManDoc: Returns StateDescriptor::TimeCenter.
    //
    StateDescriptor::TimeCenter timeType () const;
    //
    //@ManDoc: Returns number of components.
    //
    int nComp () const;
    //
    //@ManDoc: Returns the grow factor.
    //
    int nExtra () const;
    //
    //@ManDoc: Returns the name of specified component.
    //
    const std::string& name (int i) const;
    //
    int get_interp_type(int i) const;

    //@ManDoc: Returns the BCRec of specified component.
    //
    const BCRec& getBC (int i) const;
    //
    //@ManDoc: Returns all BCRecs.
    //
    const Array<BCRec>& getBCs () const;
    //
    //@ManDoc: Is sc>=0 \&\& sc+nc<=ncomp ?
    //
    int inRange (int sc, int nc) const;
    //
protected:

    IndexType          type;     // Cell centered, node centered ...
    TimeCenter         t_type;   // Temporal centering
    int                id;       // Unique id
    int                ncomp;    // Number of components
    int                ngrow;    // Grow factor
    Array<std::string> names;    // Printable names of components
    Array<int>         interp_type_array;
    Array<BCRec>       bc;       // Array of bndry types for entire level
    
};

//
//@Man:
//@Memo: A List of StateDescriptors.
/*@Doc:

  A container class for StateDescriptors.
*/

class DescriptorList
{
public:
    //
    //@ManDoc: The constructor.
    //
    DescriptorList ();
    //
    //@ManDoc: Set the list to its default state.
    //
    void clear ();
    //
    //@ManDoc: Returns number of elements in the list.
    //
    int size () const;
    //
    //@ManDoc: Adds new StateDescriptor at index `indx' to list.
    //
    void addDescriptor (int                         indx,
                        IndexType                   typ,
                        StateDescriptor::TimeCenter ttyp,
                        int                         nextra,
                        int                         num_comp);
    //
    //@ManDoc: Calls setComponent() on StateDescriptor at index `indx'.
    //
    void setComponent (int                               indx,
                       int                               comp,
                       const std::string&                nm,
                       const BCRec&                      bc,
                       int interp_type);
    //
    //@ManDoc: Returns StateDescriptor at index `k'.
    //
    const StateDescriptor& operator[] (int k) const;

private:
    //
    // These are disallowed.
    //
    DescriptorList (const DescriptorList&);
    DescriptorList& operator= (const DescriptorList&);

    PArray<StateDescriptor> desc;
};

#endif /*_StateDescriptor_H_*/
