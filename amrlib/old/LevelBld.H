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


#ifndef _LEVELBLD_H_
#define _LEVELBLD_H_

//
// $Id: LevelBld.H,v 1.6 2000/10/02 20:48:42 lijewski Exp $
//

#include <REAL.H>
#include <Box.H>
#include <Geometry.H>
#include <Amr.H>
#include <AmrLevel.H>

extern LevelBld* getLevelBld ();

//
//@Man:
//@Memo: Builds problem-specific AmrLevels.
/*@Doc:

  Abstract base class specifying an interface for building problem-specific
  AmrLevels.  
*/

class LevelBld
{
public:

    /*@ManDoc: Perform any problem-dependent setup such as physical
               boundary condition and derived quantities.
               This is a pure virtual function and hence MUST
               be implemented by derived classes.
    */
    virtual void variableSetUp () = 0;

    /*@ManDoc: Perform any problem-dependent cleanup.
               This is a pure virtual function and hence MUST
               be implemented by derived classes.
    */
    virtual void variableCleanUp () = 0;

    /*@ManDoc: This is a virtual constructor for types derived
               from AmrLevel.  The derived type is initialized
               with the default constructor.
               This is a pure virtual function and hence MUST
               be implemented by derived classes.
    */
    virtual AmrLevel* operator() () = 0;

    /*@ManDoc: This is a virtual constructor for types derived
               from AmrLevel.  The derived type is initialized
               with the five specified variables.
               This is a pure virtual function and hence MUST
               be implemented by derived classes.
    */
    virtual AmrLevel* operator() (Amr&            papa,
                                  int             lev,
                                  const Geometry& geom_lev,
                                  const BoxArray& ba,
                                  Real            time) = 0;
};

#endif /*_LEVELBLD_H_*/
