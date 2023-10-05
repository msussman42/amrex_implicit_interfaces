#ifndef CNS_TAGGING_H_
#define CNS_TAGGING_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_TagBox.H>
#include <cmath>

AMREX_GPU_HOST_DEVICE
inline
void
cns_tag_perror (int i, int j, int k,
                amrex::Array4<char> const& tag,
                amrex::Array4<amrex::Real const> const& p,
                amrex::Real pgrad_threshold, char tagval) noexcept
{
    amrex::Real ax = amrex::Math::abs(p(i+1,j,k) - p(i,j,k));
    amrex::Real ay = amrex::Math::abs(p(i,j+1,k) - p(i,j,k));
    amrex::Real az = amrex::Math::abs(p(i,j,k+1) - p(i,j,k));
    ax = amrex::max(ax,amrex::Math::abs(p(i,j,k) - p(i-1,j,k)));
    ay = amrex::max(ay,amrex::Math::abs(p(i,j,k) - p(i,j-1,k)));
    az = amrex::max(az,amrex::Math::abs(p(i,j,k) - p(i,j,k-1)));
    if (amrex::max(ax,ay,az) >= pgrad_threshold) {
        tag(i,j,k) = tagval;
    }
}

#endif
