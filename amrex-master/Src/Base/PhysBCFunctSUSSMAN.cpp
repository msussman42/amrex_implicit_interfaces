
#include <PhysBCFunctSUSSMAN.H>

PhysBCFunctSUSSMAN::PhysBCFunctSUSSMAN (
  const Geometry& geom) : m_geom(geom)
{ }

void
PhysBCFunctSUSSMAN::define (const Geometry& geom)
{
    m_geom = geom;
}

void
PhysBCFunctSUSSMAN::FillBoundary (
  int level,
  MultiFab& mf, 
  Real time,
  int dcomp, 
  Array<int> scompBC_map,
  int ncomp,
  int bfact) {

 BoxLib::Error("this routine should never be called");

}
