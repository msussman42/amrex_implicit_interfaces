#ifndef CNS_K_H_
#define CNS_K_H_

#include <AMReX_FArrayBox.H>
#include <limits>
#include <cmath>

#include "cns_prob.H"
#include "CNS_parm.H"

AMREX_GPU_HOST_DEVICE
inline
amrex::Real
cns_estdt (amrex::Box const& bx, 
	   amrex::Array4<Real const> const& /*state*/,
           amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
           Parm const& /*parm*/,
	   ProbParm const& prob_parm) noexcept
{
    using amrex::Real;

    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
    Real dt = std::numeric_limits<Real>::max();
#else
    Real dt = Real(1.e37);
#endif

    for (int dir=0;dir<AMREX_SPACEDIM;dir++) {
     AMREX_ALWAYS_ASSERT(dx[dir]>0.0);
    }
    AMREX_ALWAYS_ASSERT(prob_parm.sound_speed>0.0);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {

              Real dtx = dx[0]/(amrex::Math::abs(prob_parm.sound_speed));
#if (AMREX_SPACEDIM==1)
              dt = amrex::min(dt,dtx);
#elif ((AMREX_SPACEDIM==2)||(AMREX_SPACEDIM==3))
              Real dty = dx[1]/(amrex::Math::abs(prob_parm.sound_speed));
#if (AMREX_SPACEDIM==2)
              dt = amrex::min(dt,amrex::min(dtx,dty));
#elif (AMREX_SPACEDIM==3)
              Real dtz = dx[AMREX_SPACEDIM-1]/
		  (amrex::Math::abs(prob_parm.sound_speed));
              dt = amrex::min(dt,amrex::min(dtx,amrex::min(dty,dtz)));
#else
	      amrex::Error("AMREX_SPACEDIM problem");
#endif
#else
	      amrex::Error("AMREX_SPACEDIM problem");
#endif


            }
        }
    }
    AMREX_ALWAYS_ASSERT(dt>0.0);

    return dt;
}


#endif
