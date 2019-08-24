#ifndef _PhysBCFunct_H_
#define _PhysBCFunct_H_

#include <BCRec.H>
#include <Geometry.H>
#include <MultiFab.H>
#include <ArrayLim.H>
#include <Pointers.H>

class PhysBCFunctBase {
public:

    virtual void FillBoundary (
      int level,
      MultiFab& mf, 
      Real time,
      int dcomp, 
      Array<int> scompBC_map,
      int ncomp, 
      int bfact) = 0;

    virtual ~PhysBCFunctBase () = 0;
};

inline PhysBCFunctBase::~PhysBCFunctBase () {}

class PhysBCFunct
    : public PhysBCFunctBase
{
public:

    PhysBCFunct () {}

    virtual ~PhysBCFunct () {}

    PhysBCFunct (
      const Geometry& geom);

    void define (
      const Geometry& geom);

    virtual void FillBoundary (
      int level,
      MultiFab& mf, 
      Real time,
      int dcomp,
      Array<int> scompBC_map,
      int ncomp, 
      int bfact) override;

    using PhysBCFunctBase::FillBoundary;
private:
    Geometry m_geom;
};


#endif
