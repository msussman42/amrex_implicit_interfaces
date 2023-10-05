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

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {

              AMREX_ASSERT(dx[0]>0.0);
              AMREX_ASSERT(dx[1]>0.0);
              AMREX_ASSERT(dx[2]>0.0);
              AMREX_ASSERT(prob_parm.sound_speed>0.0);

              Real dtx = dx[0]/(amrex::Math::abs(prob_parm.sound_speed));
              Real dty = dx[1]/(amrex::Math::abs(prob_parm.sound_speed));
              Real dtz = dx[2]/(amrex::Math::abs(prob_parm.sound_speed));
              dt = amrex::min(dt,amrex::min(dtx,amrex::min(dty,dtz)));

	      AMREX_ASSERT(dt>0.0);
            }
        }
    }

    return dt;
}


#endif
