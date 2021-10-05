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
     int* bfact,
     const int* level,
     const int* gridno);

     //sweep==0 => find closest point map
     //sweep==1 => set "nice" values far from the interface.
    void fort_closest_point_map(
     int* tid, //thread id
     int* sweep,  //sweep=0 or 1
     int* ngrow,
     Real* minLS,
     Real* maxLS,
     Real* ls,AMREX_ARLIM_P(lslo),AMREX_ARLIM_P(lshi),
     Real* ls_out,
     AMREX_ARLIM_P(ls_outlo),AMREX_ARLIM_P(ls_outhi),
     Real* ls_grad_out,
     AMREX_ARLIM_P(ls_grad_outlo),
     AMREX_ARLIM_P(ls_grad_outhi),
     const Real* xlo,const Real* dx,
     Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,
     int* bfact,
     const int* level);


    void fort_combinetrianglessingle(
     int* grids_per_level,int* finest_level,
     int* nsteps,int* arrdim);

}
}/* namespace amrex */
    
#endif
#endif /*_MARCHING_TETRA_F_H_*/
