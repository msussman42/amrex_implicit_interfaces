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

#ifndef _Derive_H_
#define _Derive_H_
//
// $Id: Derive.H,v 1.15 2001/08/01 21:50:45 lijewski Exp $
//
#include <list>
#include <string>

#include <ArrayLim.H>
#include <REAL.H>
#include <Box.H>

class DescriptorList;

extern "C"
{
    //
    // Type of extern "C" function called by DeriveRec to 
    // compute derived quantity.
    //
    typedef void (*DeriveFunc)(
      Real* data, ARLIM_P(dlo), ARLIM_P(dhi),
      const int* nvar, const Real* compdat,
      ARLIM_P(compdat_lo), ARLIM_P(compdat_hi),
      const int* ncomp,
      const int* lo, const int* hi,const int* ngrow,
      const int* domain_lo, const int* domain_hi,
      const Real* delta, const Real* xlo,
      const Real* time, const Real* dt, 
      const int* bcrec,const int* level,const int* grid_no) ;
}

//
//@Man:
//@Memo: Derived Type Record
/*@Doc:

  Computes quantities derived from state data.

  DeriveRec is designed to compute quantities which can be derived
  from the state data contained in AmrLevel and its derivatives. Some
  examples might be kinetic energy, vorticity, concentration gradients ...
*/

class DeriveRec
{
   friend class DeriveList;

public:
    //
    //@ManDoc: A pointer to function taking and returning a Box. 
    //
    typedef Box (*DeriveBoxMap)(const Box&);
    //
    //@ManDoc: The destructor.
    //
    ~DeriveRec ();
    //
    //@ManDoc: The name of the derived type.
    //
    const std::string& name () const;
    //
    //@ManDoc: The names of components
    //
    const std::string& variableName (int comp) const;
    //
    //@ManDoc: The IndexType of the derived type.
    //
    IndexType deriveType () const;
    //
    //@ManDoc: The DeriveFunc used to calculate the derived type.
    //
    DeriveFunc derFunc () const;
    //
    //@ManDoc: Maps state data box to derived data box.
    //
    DeriveBoxMap boxMap () const;
    //
    //@ManDoc: Number of components in the derived type.
    //
    int numDerive () const;
    //
    //@ManDoc: Number of different chunks of state data needed for derived type.
    //
    int numRange () const;
    //
    //@ManDoc: Total number of state variables needed for derived type.
    //
    int numState () const;
    //
    //@ManDoc: The boundary conditions.
    //
    const int* getBC () const;

    /*@ManDoc: Sets state\_indx, src\_comp and num\_comp for the kth
               range (or chunk) of state data components needed to calculate
               the derived quantity.
    */
    void getRange (int  k,
                   int& state_indx,
                   int& src_comp,
                   int& num_comp) const;
protected:
    //
    // Constructor.
    //
    DeriveRec (const std::string& name,
               IndexType      result_type,
               int            nvar_derive,
               DeriveFunc     der_func,
               DeriveBoxMap   box_map);
    //
    // Constructor.
    //
    DeriveRec (const std::string&  name,
               IndexType           result_type,
               int                 nvar_derive,
	       Array<std::string>& var_names,
               DeriveFunc          der_func,
               DeriveBoxMap        box_map);

    void addRange (const DescriptorList& d_list,
                   int                   state_indx,
                   int                   src_comp,
                   int                   num_comp);

    void buildBC (const DescriptorList& d_list);
    //
    // An element of a linked list to point to state quantities in AmrLevels.
    //
    struct StateRange
    {
        int         typ;
        int         sc;
        int         nc;
        StateRange* next;
    };
    //
    // Name of derived quantity.
    //
    std::string derive_name;
    //    
    // Names of derived variables
    //
    Array<std::string> variable_names;  
    //
    // Type of derived quantity.
    //
    IndexType der_type;
    //
    // Number of components in derived quantity.
    //
    int n_derive;
    //
    // Function that computes derived quantity from state variables.
    //
    DeriveFunc func;
    //
    // Box mapper that specifies constituent region given derived region.
    //
    DeriveBoxMap bx_map;
    //
    // Total number of state variables.
    //
    int n_state;
    //
    // Number of state ranges.
    //
    int nsr;
    //
    // List of state data subranges required to derive given quantity.
    //
    StateRange* rng;
    //
    // Array of bndry types.
    //
    int* bcr;
};

//
//@Man:
//@Memo: A list of DeriveRecs.
/*@Doc:

  DeriveList manages and provides access to the list of DeriveRecs.
*/

class DeriveList
{
public:
    //
    //@ManDoc: The default constructor.
    //
    DeriveList ();
    //
    //@ManDoc: Determines whether quantity identified by <name> is in the registry.
    //
    bool canDerive (const std::string& name) const;
    //
    //@ManDoc: Access the particular record in registry.
    //
    const DeriveRec* get (const std::string& name) const;
    //
    //@ManDoc: Adds another entry to the registry.
    //
    void add (const std::string&      name,
              IndexType               result_type,
              int                     nvar_derive,
              DeriveFunc              der_func,
              DeriveRec::DeriveBoxMap box_map);
    //
    //@ManDoc: Adds another entry to the registry.
    //
    void add (const std::string&      name,
              IndexType               result_type,
              int                     nvar_derive,
	      Array<std::string>&     var_names,
              DeriveFunc              der_func,
              DeriveRec::DeriveBoxMap box_map);
    //
    //@ManDoc: Adds another StateRange to the DeriveRec identified by <name>.
    //
    void addComponent (const std::string&    name,
                       const DescriptorList& d_list,
                       int                   state_indx,
                       int                   start_comp,
                       int                   ncomp);
 
    std::list<DeriveRec>& dlist ();

private:
    //
    // Disallowed.
    //
    DeriveList (const DeriveList&);
    DeriveList& operator= (const DeriveList&);

    std::list<DeriveRec> lst;
};

#endif /*_Derive_H_*/
