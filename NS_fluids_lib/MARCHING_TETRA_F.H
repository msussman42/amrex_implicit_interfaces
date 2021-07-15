#ifndef _MARCHING_TETRA_F_H_
#define _MARCHING_TETRA_F_H_

#ifdef BL_LANG_FORT
#else

#include <AMReX_ArrayLim.H>

namespace amrex{

extern "C" {


    void fort_isogrid(
     int* tid,
     int* visual_tessellate_vfrac,
     Real* recon,AMREX_ARLIM_P(reconlo),AMREX_ARLIM_P(reconhi),
     const Real* xlo,const Real* dx,
     Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,int* bfact,
     const int* level,const int* gridno,int* nmat);

    void fort_combinetriangles(
     int* grids_per_level,int* finest_level,
     int* nsteps,int* imaterial,int* arrdim,
     Real* time,int* plotint);


    void fort_isogridsingle(
     Real* ls,AMREX_ARLIM_P(lslo),AMREX_ARLIM_P(lshi),
     const Real* xlo,const Real* dx,
     Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
     const int* lo,const int* hi,
     const int* level,const int* gridno);

    void fort_combinetrianglessingle(
     int* grids_per_level,int* finest_level,
     int* nsteps,int* arrdim);

}
}/* namespace amrex */
    
#endif
#endif /*_MARCHING_TETRA_F_H_*/
