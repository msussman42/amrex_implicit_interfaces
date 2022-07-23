#ifndef AMREX_ARRAYLIM_SUSSMAN_H_
#define AMREX_ARRAYLIM_SUSSMAN_H_

#include <AMReX_REAL.H>

#ifdef __cplusplus

#else

#if !defined(BL_LANG_FORT)

/* C stuff */

#else

#if __STDC__==1 || defined(__INTEL_COMPILER)


#if (AMREX_SPACEDIM == 2)

#define DIMS3D(a) a##_lo, a##_hi
#define DIMDEC3D(a) a##_lo(3), a##_hi(3)
#define DIMV3D(a) a##_lo(1):a##_hi(1), a##_lo(2):a##_hi(2), a##_lo(3):a##_hi(3)
#define ARG3D_L1(a) a##_lo(1)
#define ARG3D_L2(a) a##_lo(2)
#define ARG3D_L3(a) a##_lo(3)
#define ARG3D_H1(a) a##_hi(1)
#define ARG3D_H2(a) a##_hi(2)
#define ARG3D_H3(a) a##_hi(3)

#endif

#if (AMREX_SPACEDIM == 3)

#define DIMS3D(a) a##_lo, a##_hi
#define DIMDEC3D(a) a##_lo(3), a##_hi(3)
#define DIMV3D(a) a##_lo(1):a##_hi(1), a##_lo(2):a##_hi(2), a##_lo(3):a##_hi(3)
#define ARG3D_L1(a) a##_lo(1)
#define ARG3D_L2(a) a##_lo(2)
#define ARG3D_L3(a) a##_lo(3)
#define ARG3D_H1(a) a##_hi(1)
#define ARG3D_H2(a) a##_hi(2)
#define ARG3D_H3(a) a##_hi(3)

#endif

#else

#if (AMREX_SPACEDIM == 2)

#define DIMS3D(a) a/**/_l1, a/**/_l2, a/**/_l3, a/**/_h1, a/**/_h2, a/**/_h3
#define DIMDEC3D(a) a/**/_l1, a/**/_l2, a/**/_l3, a/**/_h1, a/**/_h2, a/**/_h3
#define DIMV3D(a) a/**/_l1:a/**/_h1, a/**/_l2:a/**/_h2, a/**/_l3:a/**/_h3
#define ARG3D_L1(a) a/**/_l1
#define ARG3D_L2(a) a/**/_l2
#define ARG3D_L3(a) a/**/_l3
#define ARG3D_H1(a) a/**/_h1
#define ARG3D_H2(a) a/**/_h2
#define ARG3D_H3(a) a/**/_h3

#endif

#if (AMREX_SPACEDIM == 3)

#define DIMS3D(a) a/**/_l1, a/**/_l2, a/**/_l3, a/**/_h1, a/**/_h2, a/**/_h3
#define DIMDEC3D(a) a/**/_l1, a/**/_l2, a/**/_l3, a/**/_h1, a/**/_h2, a/**/_h3
#define DIMV3D(a) a/**/_l1:a/**/_h1, a/**/_l2:a/**/_h2, a/**/_l3:a/**/_h3
#define ARG3D_L1(a) a/**/_l1
#define ARG3D_L2(a) a/**/_l2
#define ARG3D_L3(a) a/**/_l3
#define ARG3D_H1(a) a/**/_h1
#define ARG3D_H2(a) a/**/_h2
#define ARG3D_H3(a) a/**/_h3

#endif

#endif /*__STDC__==1*/

#endif /*else of !defined(BL_LANG_FORT)*/

#endif /*else of __cplusplus*/

#endif /*BL_ARRAYLIM_H*/
