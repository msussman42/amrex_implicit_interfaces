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


AMREX_GPU_DEVICE
inline
void
cns_compute_errors (int i, int j, int k, 
  amrex::Array4<amrex::Real> const& snewfab,
  amrex::Array4<amrex::Real> const& errorfab,
  amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
  amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo,
  const amrex::Real time)
{
   using amrex::Real;

   AMREX_ASSERT(snewfab.nComp()==AMREX_SPACEDIM+1);
   AMREX_ASSERT(errorfab.nComp()==2);


   Real volume=dx[0];
#if (AMREX_SPACEDIM==2)
   volume=volume*dx[1];
#endif
#if (AMREX_SPACEDIM==3)
   volume=volume*dx[1]*dx[2];
#endif
   Real x = prob_lo[0] + (i+Real(0.5))*dx[0];

   errorfab(i,j,k,0)=volume;

   Real p_ref=Real(0.0);
   if (time<Real(1.0)-x) {
    p_ref=x;
   } else if (time<Real(1.0)+x) {
    p_ref=Real(2.0)*x+time-Real(1.0);
   } else if (time<Real(3.0)-x) {
    p_ref=Real(3.0)*x;
   }

   Real p_act=snewfab(i,j,k,0);
   Real local_error=std::abs(p_ref-p_act)*volume;
   errorfab(i,j,k,1)=local_error;
}


#endif
