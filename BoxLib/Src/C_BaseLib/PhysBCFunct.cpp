
#include <PhysBCFunct.H>

PhysBCFunct::PhysBCFunct (
  const Geometry& geom) : m_geom(geom)
{ }

void
PhysBCFunct::define (const Geometry& geom)
{
    m_geom = geom;
}

void
PhysBCFunct::FillBoundary (
  int level,
  MultiFab& mf, 
  Real time,
  int dcomp, 
  Array<int> scompBC_map,
  int ncomp,
  int bfact) {

 BoxLib::Error("this routine should never be called");

}
