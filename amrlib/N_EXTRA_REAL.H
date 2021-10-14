#ifndef N_EXTRA_REAL_H_
#define N_EXTRA_REAL_H_

#ifdef BL_LANG_FORT
! x,y,z is core (base); 
! (x0,y0,z0,r,u,v,w,pres,den,T,stress,
!  insert time,type_molecule,type_atom) is extra. 
! (r is a radius)
! Structure of Array components can be added real time.
#else
// x,y,z is core (base); 
// (x0,y0,z0,r,u,v,w,pres,den,T,stress,
//  insert time,type_molecule,type_atom) is extra. 
// (r is a radius)
// Structure of Array components can be added real time.
#endif

#define N_EXTRA_REAL (4*AMREX_SPACEDIM+7)
#define N_EXTRA_REAL_INSERT_TIME (4*AMREX_SPACEDIM+4)

#endif
