#ifndef N_EXTRA_REAL_H_
#define N_EXTRA_REAL_H_

#if (0)
x,y,z is core (base); 
(x0,y0,z0,u,v,w,den,T,insert time) is extra
Structure of Array components can be added real time.
#endif

#define N_EXTRA_REAL (2*AMREX_SPACEDIM+3)
#define N_EXTRA_REAL_X0 (0)
#define N_EXTRA_REAL_U (AMREX_SPACEDIM)
#define N_EXTRA_REAL_DEN (2*AMREX_SPACEDIM)
#define N_EXTRA_REAL_T   (2*AMREX_SPACEDIM+1)
#define N_EXTRA_REAL_INSERT_TIME (2*AMREX_SPACEDIM+2)

#endif
