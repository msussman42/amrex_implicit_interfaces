#ifndef N_EXTRA_REAL_H_
#define N_EXTRA_REAL_H_

#if (0)
x,y,z is core (base); insert time (REAL), x0,y0,z0 (REAL) 
and material_id (int) are extra
Structure of Array components can be added real time.
#endif

#define N_EXTRA_REAL (1+AMREX_SPACEDIM)
#define N_EXTRA_REAL_INSERT_TIME (0)
#define N_EXTRA_REAL_X0 (N_EXTRA_REAL_INSERT_TIME+1)
#define N_EXTRA_REAL_Y0 (N_EXTRA_REAL_X0+1)
#define N_EXTRA_REAL_Z0 (N_EXTRA_REAL_X0+AMREX_SPACEDIM-1)

#define N_EXTRA_INT (1)
#define N_EXTRA_INT_MATERIAL_ID (0)

#endif
