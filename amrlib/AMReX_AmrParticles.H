#ifndef AMREX_AmrParticles_H_
#define AMREX_AmrParticles_H_

#include <AMReX_Particles.H>
#include <AMReX_AmrParGDB.H>

namespace amrex {

template <int NStructReal, int NStructInt=0, int NArrayReal=0, int NArrayInt=0>
class AmrParticleContainer
        : public ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
{

public:

    typedef Particle<NStructReal, NStructInt> ParticleType;

      //SUSSMAN    
    AmrParticleContainer (Amr* amr_core)
        : ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>(amr_core->GetParGDB())
    {
    }

    AmrParticleContainer (const Vector<Geometry>            & geom, 
                          const Vector<DistributionMapping> & dmap,
                          const Vector<BoxArray>            & ba,
                          const Vector<int>                 & rr)
        : ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>(geom, dmap, ba, rr)
    {
    }
    
    ~AmrParticleContainer () {}
    
};


}

#endif
