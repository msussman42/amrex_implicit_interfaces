
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

	    int i=iv[0];
	    int j=0;
	    int k=0;
#if (AMREX_SPACEDIM==2)||(AMREX_SPACEDIM==3)
	    j=iv[1];
#if (AMREX_SPACEDIM==3)
	    k=iv[2];
#endif
#endif

	    dest(i,j,k,dcomp)=((i<0) ? Real(0.0) : Real(1.0)+time);

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
