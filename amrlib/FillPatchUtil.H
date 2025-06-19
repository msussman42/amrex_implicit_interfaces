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
      PhysBCFunctBaseSUSSMAN& cbc, 
      PhysBCFunctBaseSUSSMAN& fbc,
      Interpolater* mapper, 
      const Vector<BCRec>& global_bcs,
      Vector<int> scompBC_map,
      int levelc,int levelf,
      int bfactc,int bfactf, 
      int grid_type,
      int debug_fillpatch);

//SUSSMAN
    void FillPatchTower (
      int ngrow_root,
      int called_from_regrid,
      int finest_top_level,
      int top_level,
      MultiFab& mf_target, //target
      Real time,
      Vector<MultiFab*> tower_data,
      int scomp, 
      int dcomp, 
      int ncomp,
      Vector<const Geometry*> tower_geom,
      Vector<PhysBCFunctBaseSUSSMAN*> tower_physbc,
      Interpolater* mapper, 
      const Vector<BCRec>& bcs,
      Vector<int> scompBC_map,
      Vector<int> tower_bfact,
      int grid_type,
      int debug_fillpatch);

} // namespace amrex

#endif
