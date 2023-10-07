#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include "CNS_parm.H"
#include "cns_prob_parm.H"

AMREX_GPU_DEVICE
inline
void
cns_initdata (int i, int j, int k, 
  amrex::Array4<amrex::Real> const& state,
  amrex::GeometryData const& /*geomdata*/,
  amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
  amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo,
  Parm const& /*parm*/, 
  ProbParm const& /*prob_parm*/)
{
    using amrex::Real;

    AMREX_ASSERT(state.nComp()==AMREX_SPACEDIM+1);
    Real x[AMREX_SPACEDIM];
    x[0] = prob_lo[0] + (i+Real(0.5))*dx[0];
#if ((AMREX_SPACEDIM==2)||(AMREX_SPACEDIM==3))
    x[1] = prob_lo[1] + (j+Real(0.5))*dx[1];
#if (AMREX_SPACEDIM==3)
    x[2] = prob_lo[2] + (k+Real(0.5))*dx[2];
#endif
#endif
     //which is optimal on a GPU vs CPU? single or double precision?
     //AMD? INTEL? NVIDEA?
     //number of cores?
     //architecture?
     //GPU memory?
     //max_grid_size
     //better to have the compiler take care of "if" statements or manually?
    state(i,j,k,0) = x[0];
    state(i,j,k,1) = Real(1.0);
    for (int dir=2;dir<=AMREX_SPACEDIM;dir++) {
     state(i,j,k,dir) = Real(0.0);
    }
}

#endif
