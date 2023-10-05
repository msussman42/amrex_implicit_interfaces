
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>

using namespace amrex;

//GeometryData is declared in:
//Src/Base/AMReX_Geometry.H
//BCRec is declared in:
//Src/Base/AMReX_BCRec.H
struct CnsFillExtDir
{
    AMREX_GPU_DEVICE
    void operator() (const IntVect& iv, Array4<Real> const& dest,
                     const int dcomp, const int num_comp,
                     GeometryData const& geom, const Real time,
                     const BCRec* bcr, const int bcomp,
                     const int scomp) const
        {
            // do something for external Dirichlet (BCType::ext_dir)
	    const Real* prob_lo=geom.ProbLo();
	    const Real* prob_hi=geom.ProbHi();
	    const Real* dx=geom.CellSize();
	    Real x[3];
	    int i=iv[0];
	    int j=iv[1];
	    int k=iv[2];
	    x[0] = prob_lo[0] + (i+Real(0.5))*dx[0];
            x[1] = prob_lo[1] + (j+Real(0.5))*dx[1];
            x[2] = prob_lo[2] + (k+Real(0.5))*dx[2];

	    if (x[0]<prob_lo[0]) {
             for (int nc=scomp;nc<scomp+num_comp;nc++) {
	      const int* lo_bc=bcr[nc-scomp+bcomp].lo();
  	      if (lo_bc[0]==BCType::ext_dir) {
		      AMREX_ASSERT(nc==0);
		      dest(i,j,k,nc-scomp+dcomp)=Real(0.0);
	      }
	     }
	    }

	    if (x[0]>prob_hi[0]) {
             for (int nc=scomp;nc<scomp+num_comp;nc++) {
  	      const int* hi_bc=bcr[nc-scomp+bcomp].hi();
  	      if (hi_bc[0]==BCType::ext_dir) {
		      AMREX_ASSERT(nc==0);
		      dest(i,j,k,nc-scomp+dcomp)=Real(1.0)+time;
	      }
	     }
	    }

        }
};

// bx                  : Cells outside physical domain and inside bx are filled.
// data, dcomp, numcomp: Fill numcomp components of data starting from dcomp.
// bcr, bcomp          : bcr[bcomp] specifies BC for component dcomp and so on.
// scomp               : component index for dcomp as in the desciptor set up in CNS::variableSetUp.

void cns_bcfill (Box const& bx, FArrayBox& data,
                 const int dcomp, const int numcomp,
                 Geometry const& geom, const Real time,
                 const Vector<BCRec>& bcr, const int bcomp,
                 const int scomp)
{
    GpuBndryFuncFab<CnsFillExtDir> gpu_bndry_func(CnsFillExtDir{});
    gpu_bndry_func(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);
}
