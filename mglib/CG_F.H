/*
** (c) 1996-2000 The Regents of the University of California (through
** E.O. Lawrence Berkeley National Laboratory), subject to approval by
** the U.S. Department of Energy.  Your use of this software is under
** license -- the license agreement is attached and included in the
** directory as license.txt or you may contact Berkeley Lab's Technology
** Transfer Department at TTD@lbl.gov.  NOTICE OF U.S. GOVERNMENT RIGHTS.
** The Software was developed under funding from the U.S. Government
** which consequently retains certain rights as follows: the
** U.S. Government has been granted for itself and others acting on its
** behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
** Software to reproduce, prepare derivative works, and perform publicly
** and display publicly.  Beginning five (5) years after the date
** permission to assert copyright is obtained from the U.S. Department of
** Energy, and subject to any subsequent five (5) year renewals, the
** U.S. Government is granted for itself and others acting on its behalf
** a paid-up, nonexclusive, irrevocable, worldwide license in the
** Software to reproduce, prepare derivative works, distribute copies to
** the public, perform publicly and display publicly, and to permit
** others to do so.
*/

#ifndef _CG_F_H_
#define _CG_F_H_

/*
** $Id: CG_F.H,v 1.3 2000/08/24 16:02:46 car Exp $
*/

#include <REAL.H>

#ifdef BL_LANG_FORT
#define FORT_CGUPDATE  cgupdate3dgen
#define FORT_CGADVCP   cgadvcp3dgen
#define FORT_CGXDOTY   cgxdoty3dgen
#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_CGUPDATE  CGUPDATE3DGEN
#    define FORT_CGADVCP   CGADVCP3DGEN
#    define FORT_CGXDOTY   CGXDOTY3DGEN
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_CGUPDATE  cgupdate3dgen
#    define FORT_CGADVCP   cgadvcp3dgen
#    define FORT_CGXDOTY   cgxdoty3dgen
#  else
#    define FORT_CGUPDATE  cgupdate3dgen_
#    define FORT_CGADVCP   cgadvcp3dgen_
#    define FORT_CGXDOTY   cgxdoty3dgen_
#endif
#endif

#include <ArrayLim.H>

extern "C"
{
    void FORT_CGUPDATE (
        Real* phi,  ARLIM_P(phi_lo),  ARLIM_P(phi_hi),
        const Real *aa,
        const Real* tmpy, ARLIM_P(tmpy_lo), ARLIM_P(tmpy_hi),
        const Real* tmpp, ARLIM_P(tmpp_lo), ARLIM_P(tmpp_hi),
        const int* tilelo, const int* tilehi,
        const int* fablo,const int* fabhi,
        int* bfact,int* bfact_top);

    void FORT_CGADVCP (
        Real* tmpp,ARLIM_P(tmpp_lo), ARLIM_P(tmpp_hi),
        const Real* tmpr, ARLIM_P(tmpr_lo), ARLIM_P(tmpr_hi),
        const Real* tmpy, ARLIM_P(tmpy_lo), ARLIM_P(tmpy_hi),
        const Real *b,
        const int* tilelo, const int* tilehi,
        const int* fablo,const int* fabhi,
        int* bfact,int* bfact_top);

    void FORT_CGXDOTY (
      int* ncomp,
      Real *pw,
      const Real* phi,  ARLIM_P(phi_lo),  ARLIM_P(phi_hi),
      const Real* tmpw, ARLIM_P(tmpw_lo), ARLIM_P(tmpw_hi),
      const int* tilelo, const int* tilehi,
      const int* fablo,const int* fabhi,
      int* bfact,int* bfact_top);


}
#endif

#endif /*_CG_F_H_*/
