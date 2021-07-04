#ifndef N_EXTRA_REAL_H_
#define N_EXTRA_REAL_H_

#ifdef BL_LANG_FORT
! x,y,z is core (base); 
! (x0,y0,z0,r,u,v,w,den,T,insert time,vorticity,type) is extra. 
! (r is a radius)
! Structure of Array components can be added real time.
#else
// x,y,z is core (base); 
// (x0,y0,z0,r,u,v,w,den,T,insert time,vorticity,type) is extra. 
// (r is a radius)
// Structure of Array components can be added real time.
#endif

#define N_EXTRA_REAL (3*AMREX_SPACEDIM+5)

#endif
