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

#ifndef _BCREC_H_
#define _BCREC_H_ 
//
// $Id: BCRec.H,v 1.9 2001/08/01 21:50:44 lijewski Exp $
//
#include <Box.H>
#include <BC_TYPES.H>
//
//@Man:
//@Memo: Boundary Condition Records.
/*@Doc:

  Necessary information and functions for computing boundary conditions.
*/

class BCRec
{
public:
    //
    //@ManDoc: The default constructor -- sets to invalid boundary types.
    //
    BCRec ();
    //
    //@ManDoc: The constructor.
    //
    BCRec (D_DECL(int loX, int loY, int loZ),
           D_DECL(int hiX, int hiY, int hiZ));
    //
    //@ManDoc: Another constructor.
    //
    BCRec (const int* lo,
           const int* hi);

    /*@ManDoc: Yet another constructor.  Inherits bndry types from bc\_domain
               when bx lies on edge of domain otherwise gets interior Dirchlet.
    */
    BCRec (const Box&   bx,
           const Box&   domain,
           const BCRec& bc_domain);
    //
    //@ManDoc: Explicitly set lo bndry value.
    //
    void setLo (int dir,
                int bc_val);
    //
    //@ManDoc: Explicitly set hi bndry value.
    //
    void setHi (int dir,
                int bc_val);
    //
    //@ManDoc: Return bndry values (used in calls to FORTRAN).
    //
    const int* vect () const;
    //
    //@ManDoc: Return low-end boundary data.
    //
    const int* lo () const;
    //
    //@ManDoc: Return high-end boundary data.
    //
    const int* hi () const;
    //
    //@ManDoc: Return low-end boundary data in direction <dir>.
    //
    int lo (int dir) const;
    //
    //@ManDoc: Return high-end boundary data in direction <dir>.
    //
    int hi (int dir) const;
    //
    //@ManDoc: ASCII write to ostream.
    //
    friend std::ostream& operator << (std::ostream&, const BCRec&);

protected:
    //
    // Array of integer values describing boundary conditions.
    //
    int bc[2*BL_SPACEDIM];
};

namespace BoxLib
{
    //
    // Function for setting a BC.
    //
    void setBC (const Box&   bx, 
                const Box&   domain,
                const BCRec& bc_dom, 
                BCRec&       bcr);
    //
    // Function for setting array of BCs.
    //
    void setBC (const Box&          bx, 
                const Box&          domain, 
                int                 src_comp, 
                int                 dest_comp, 
                int                 ncomp,
                const Array<BCRec>& bc_dom, 
                Array<BCRec>&       bcr);
}

#endif /*_BCREC_H_*/
