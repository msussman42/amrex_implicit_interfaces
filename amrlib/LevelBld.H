
#ifndef _LEVELBLD_H_
#define _LEVELBLD_H_

#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <Amr.H>
#include <AmrLevel.H>

extern LevelBld* getLevelBld ();

//
// Builds problem-specific AmrLevels.
//
// Abstract base class specifying an interface for building problem-specific
// AmrLevels.  
//

class LevelBld
{
public:

    virtual ~LevelBld () {}
    virtual void variableSetUp () = 0;
    virtual void variableCleanUp () = 0;
    virtual AmrLevel* operator() () = 0;
    virtual AmrLevel* operator() (Amr&            papa,
                                  int             lev,
                                  const Geometry& geom_lev,
                                  const BoxArray& ba,
                                  const DistributionMapping& dm,
                                  Real            time) = 0;
};

#endif /*_LEVELBLD_H_*/
