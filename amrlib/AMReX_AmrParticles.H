#ifndef AMREX_AmrParticles_H_
#define AMREX_AmrParticles_H_
#include <AMReX_Config.H>

#include <AMReX_Particles.H>
//SUSSMAN
//#include <AMReX_TracerParticles.H>
#include <AMReX_AmrParGDB.H>

namespace amrex {

//SUSSMAN: do not include the grid transfer routines.
//
template <int NStructReal, int NStructInt=0, int NArrayReal=0, int NArrayInt=0,
          template<class> class Allocator=DefaultAllocator>
class AmrParticleContainer
    : public ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
{

public:

    typedef Particle<NStructReal, NStructInt> ParticleType;

    AmrParticleContainer ()
        : ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>()
    {
    }

    AmrParticleContainer (AmrCore* amr_core)
        : ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>(amr_core->GetParGDB())
    {
    }

    AmrParticleContainer (const Vector<Geometry>            & geom,
                          const Vector<DistributionMapping> & dmap,
                          const Vector<BoxArray>            & ba,
                          const Vector<int>                 & rr)
        : ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>(geom, dmap, ba, rr)
    {
    }

    ~AmrParticleContainer () {}

    AmrParticleContainer ( const AmrParticleContainer &) = delete;
    AmrParticleContainer& operator= ( const AmrParticleContainer & ) = delete;

    AmrParticleContainer ( AmrParticleContainer && ) = default;
    AmrParticleContainer& operator= ( AmrParticleContainer && ) = default;
};

}

#endif
