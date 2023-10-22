#ifndef AMREX_AmrParticles_H_
#define AMREX_AmrParticles_H_
#include <AMReX_Config.H>

#include <AMReX_Particles.H>
#include <AMReX_TracerParticles.H>
#include <AMReX_AmrParGDB.H>

namespace amrex {

template <typename T_ParticleType, int NArrayReal=0, int NArrayInt=0,
          template<class> class Allocator=DefaultAllocator, class CellAssignor=DefaultAssignor>
class AmrParticleContainer_impl // NOLINT(cppcoreguidelines-virtual-class-destructor)
    : public ParticleContainer_impl<T_ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>
{

public:

    using ParticleType = T_ParticleType;

    AmrParticleContainer_impl ()
        : ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>()
    {
    }

    AmrParticleContainer_impl (AmrCore* amr_core)
        : ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>(amr_core->GetParGDB())
    {
    }

    AmrParticleContainer_impl (const Vector<Geometry>            & geom,
                               const Vector<DistributionMapping> & dmap,
                               const Vector<BoxArray>            & ba,
                               const Vector<int>                 & rr)
        : ParticleContainer_impl<ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>(geom, dmap, ba, rr)
    {
    }

    ~AmrParticleContainer_impl () override = default;

    AmrParticleContainer_impl ( const AmrParticleContainer_impl &) = delete;
    AmrParticleContainer_impl& operator= ( const AmrParticleContainer_impl & ) = delete;

    AmrParticleContainer_impl ( AmrParticleContainer_impl && ) noexcept = default;
    AmrParticleContainer_impl& operator= ( AmrParticleContainer_impl && ) noexcept = default;
};

template <int T_NStructReal, int T_NStructInt=0, int T_NArrayReal=0, int T_NArrayInt=0,
          template<class> class Allocator=DefaultAllocator, class CellAssignor=DefaultAssignor>
using AmrParticleContainer = AmrParticleContainer_impl<Particle<T_NStructReal, T_NStructInt>, T_NArrayReal, T_NArrayInt, Allocator, CellAssignor>;

class AmrTracerParticleContainer
    : public TracerParticleContainer
{
public:

    AmrTracerParticleContainer (AmrCore* amr_core)
        : TracerParticleContainer(amr_core->GetParGDB())
    {
    }

    ~AmrTracerParticleContainer () override = default;

    AmrTracerParticleContainer ( const AmrTracerParticleContainer &) = delete;
    AmrTracerParticleContainer& operator= ( const AmrTracerParticleContainer & ) = delete;

    AmrTracerParticleContainer ( AmrTracerParticleContainer && ) noexcept = default;
    AmrTracerParticleContainer& operator= ( AmrTracerParticleContainer && ) noexcept = default;
};

}

#endif
