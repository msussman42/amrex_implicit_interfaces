#ifndef CNS_PROB_H_
#define CNS_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include "CNS_index_macros.H"
#include "CNS_parm.H"
#include "cns_prob_parm.H"

AMREX_GPU_DEVICE
inline
void
cns_initdata (int i, int j, int k, amrex::Array4<amrex::Real> const& state,
              amrex::GeometryData const& geomdata, Parm const& parm, ProbParm const& prob_parm)
{
    using amrex::Real;

    const Real* prob_lo = geomdata.ProbLo();
    const Real* dx      = geomdata.CellSize();

    AMREX_ALWAYS_ASSERT(prob_parm.num_state_variables==4);
    Real x[3];
    x[0] = prob_lo[0] + i*dx[0];
    x[1] = prob_lo[1] + j*dx[1];
    x[2] = prob_lo[2] + k*dx[2];
     //which is optimal on a GPU vs CPU? single or double precision?
     //AMD? INTEL? NVIDEA?
     //number of cores?
     //architecture?
     //GPU memory?
     //max_grid_size
     //better to have the compiler take care of "if" statements or manually?
    state(i,j,k,0) = x[0];
    state(i,j,k,1) = Real(1.0);
    state(i,j,k,2) = Real(0.0);
    state(i,j,k,3) = Real(0.0);
}

#endif
