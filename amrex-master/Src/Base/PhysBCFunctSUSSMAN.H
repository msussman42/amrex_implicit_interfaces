#ifndef _PhysBCFunct_H_
#define _PhysBCFunct_H_

#include <AMReX_BCRec.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>

namespace amrex {

class PhysBCFunctBaseSUSSMAN {
public:

    virtual void FillBoundary (
      int level,
      MultiFab& mf, 
      Real time,
      int dcomp, 
      Vector<int> scompBC_map,
      int ncomp, 
      int bfact) = 0;

    virtual ~PhysBCFunctBaseSUSSMAN () = 0;
};

inline PhysBCFunctBaseSUSSMAN::~PhysBCFunctBaseSUSSMAN () {}

class PhysBCFunctSUSSMAN
    : public PhysBCFunctBaseSUSSMAN
{
public:

    PhysBCFunctSUSSMAN () {}

    virtual ~PhysBCFunctSUSSMAN () {}

    PhysBCFunctSUSSMAN (
      const Geometry& geom);

    void define (
      const Geometry& geom);

    virtual void FillBoundary (
      int level,
      MultiFab& mf, 
      Real time,
      int dcomp,
      Vector<int> scompBC_map,
      int ncomp, 
      int bfact) override;

    using PhysBCFunctBaseSUSSMAN::FillBoundary;
private:
    Geometry m_geom;
};

}
#endif
