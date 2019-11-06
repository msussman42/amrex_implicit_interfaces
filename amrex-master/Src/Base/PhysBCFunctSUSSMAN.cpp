
#include <PhysBCFunctSUSSMAN.H>

namespace amrex {

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
  Vector<int> scompBC_map,
  int ncomp,
  int bfact) {

 amrex::Error("this routine should never be called");

}

} // namespace amrex
