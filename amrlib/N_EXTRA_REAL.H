#ifndef N_EXTRA_REAL_H_
#define N_EXTRA_REAL_H_

#ifdef BL_LANG_FORT
! x,y,z is core (base); 
! (x0,y0,z0,r,u,v,w,den,T,
!  insert time,type_molecule,type_atom) is extra. 
! (r is a radius)
! Structure of Array components can be added real time.
#else
// x,y,z is core (base); 
// (x0,y0,z0,r,u,v,w,den,T,
//  insert time,type_molecule,type_atom) is extra. 
// (r is a radius)
// Structure of Array components can be added real time.
#endif

#define N_EXTRA_REAL (2*AMREX_SPACEDIM+6)
#define N_EXTRA_REAL_r (AMREX_SPACEDIM)
#define N_EXTRA_REAL_u (AMREX_SPACEDIM+1)
#define N_EXTRA_REAL_den (2*AMREX_SPACEDIM+1)
#define N_EXTRA_REAL_T   (2*AMREX_SPACEDIM+2)
#define N_EXTRA_REAL_INSERT_TIME (2*AMREX_SPACEDIM+3)
#define N_EXTRA_REAL_type_molecule (2*AMREX_SPACEDIM+4)
#define N_EXTRA_REAL_type_atom (2*AMREX_SPACEDIM+5)

#endif
