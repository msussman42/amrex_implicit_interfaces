#ifndef Particles_F_H_
#define Particles_F_H_

#include <BLFort.H>

#ifdef USE_PARTICLES

extern "C"
{

    void BoxLib_particle_set_position (Real*, int ns, int np,
                                      D_DECL(const Real* x,
                                             const Real* y,
                                             const Real* z));

    void BoxLib_particle_get_position (const Real*, int ns, int np,
                                      D_DECL(Real* x,
                                             Real* y,
                                             Real* z));

    void BoxLib_deposit_cic(const Real*, int ns, int np,
                           Real* rho, const int* lo, const int* hi,
                           const Real* plo, const Real* dx);

    void BoxLib_interpolate_cic(const Real*, int ns, int np,
                               Real* rho, const int* lo, const int* hi, int ncomp,
                               const Real* plo, const Real* dx);

BL_FORT_PROC_DECL(PART_SUMMASSDOWN,part_summassdown)
    (BL_FORT_FAB_ARG(crse_fab),
     const BL_FORT_FAB_ARG(fine_fab),
     const int* ovlo, const int* ovhi, const int* rat);

BL_FORT_PROC_DECL(PART_SUMMASSUP,part_summassup)
    (const BL_FORT_FAB_ARG(crse_fab),
     BL_FORT_FAB_ARG(fine_fab),
     const int* ovlo, const int* ovhi, const int* rat);

BL_FORT_PROC_DECL(PART_SUMDENSUP,part_sumdensup)
    (const BL_FORT_FAB_ARG(crse_fab),
     BL_FORT_FAB_ARG(fine_fab),
     const int* ovlo, const int* ovhi, const int* rat);

#if (BL_SPACEDIM==2)
    BL_FORT_PROC_DECL(BICUBICINTERP,bicubicinterp)
    (const Real* x, const Real* y, const Real* fab,
     const int* lo, const int* hi,
     const Real* plo, const Real* phi, const Real* dx, Real* z);
#endif

};

#endif

#endif /*_Particles_F_H_*/
