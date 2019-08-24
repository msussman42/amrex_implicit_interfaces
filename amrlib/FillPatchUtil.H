#ifndef _FillPatchUtil_H_
#define _FillPatchUtil_H_

#include <MultiFab.H>
#include <Geometry.H>
#include <PhysBCFunct.H>
#include <Interpolater.H>

namespace BoxLib
{

    void FillPatchSingleLevel (
      int level,
      MultiFab& mf, 
      Real time, 
      const MultiFab& smf,
      int scomp, 
      int dcomp, 
      int ncomp,
      const Geometry& geom, 
      PhysBCFunctBase& physbcf,
      Array<int> scompBC_map,
      int bfact);

    void FillPatchTwoLevels (
      MultiFab& mf, 
      Real time,
      MultiFab& cmf, 
      MultiFab& fmf, 
      int scomp, int dcomp, int ncomp,
      const Geometry& cgeom, const Geometry& fgeom, 
      PhysBCFunctBase& cbc, PhysBCFunctBase& fbc,
      Interpolater* mapper, const Array<BCRec>& bcs,
      Array<int> scompBC_map,
      int levelc,int levelf,
      int bfactc,int bfactf);

}

#endif
