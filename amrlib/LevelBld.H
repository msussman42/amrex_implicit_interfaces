
#ifndef _LEVELBLD_H_
#define _LEVELBLD_H_

#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_AmrCore.H>
#include <AmrLevel.H>

namespace amrex {


class LevelBld
{
public:

    virtual ~LevelBld () {}
    virtual void variableSetUp () = 0;
    virtual void variableCleanUp () = 0;
    virtual AmrLevel* operator() () = 0;
    virtual AmrLevel* operator() (AmrCore&        papa,
                                  int             lev,
                                  const Geometry& geom_lev,
                                  const BoxArray& ba,
                                  const DistributionMapping& dm,
                                  Real            time) = 0;
};

}

extern AMREX_ATTRIBUTE_WEAK amrex::LevelBld* getLevelBld ();

#endif /*_LEVELBLD_H_*/
