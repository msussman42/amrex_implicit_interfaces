#ifndef N_EXTRA_REAL_H_
#define N_EXTRA_REAL_H_

#ifdef BL_LANG_FORT
! x,y,z is not extra; (x0,y0,z0,r,u,v,w,den,T,insert time) is extra. 
! (r can be a distance, or a radius)
! Structure of Array components can be added real time.
#define N_EXTRA_REAL (2*AMREX_SPACEDIM+4)
#else
// x,y,z is not extra; (x0,y0,z0,r,u,v,w,den,T,insert time) is extra. 
// (r can be a distance, or a radius)
// Structure of Array components can be added real time.
#define N_EXTRA_REAL (2*AMREX_SPACEDIM+4)
#  endif

#endif
