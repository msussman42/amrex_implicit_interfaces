#ifndef AMREX_FillPatchUtil_H_
#define AMREX_FillPatchUtil_H_

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <PhysBCFunctSUSSMAN.H>
#include <AMReX_Array.H>
#include <Interpolater.H>

namespace amrex
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
      PhysBCFunctBaseSUSSMAN& physbcf,
      Vector<int> scompBC_map,
      int bfact, 
      int debug_fillpatch);

    void FillPatchTwoLevels (
      MultiFab& mf, 
      Real time,
      MultiFab& cmf, 
      MultiFab& fmf, 
      int scomp, int dcomp, int ncomp,
      const Geometry& cgeom, const Geometry& fgeom, 
      PhysBCFunctBaseSUSSMAN& cbc, PhysBCFunctBaseSUSSMAN& fbc,
      Interpolater* mapper, const Vector<BCRec>& bcs,
      Vector<int> scompBC_map,
      int levelc,int levelf,
      int bfactc,int bfactf, 
      int grid_type,
      int debug_fillpatch);

} // namespace amrex

#endif
