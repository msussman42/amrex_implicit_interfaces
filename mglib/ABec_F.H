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

#ifndef _AABBEC_F_H_
#define _AABBEC_F_H_

/*
** $Id: ABec_F.H,v 1.4 2000/08/24 20:28:28 car Exp $
*/

#include <REAL.H>

#ifdef BL_LANG_FORT
#define FORT_GSRB          gsrb
#define FORT_ADOTX         adotx
#define FORT_DIAGSUM       diagsum
#else
#  ifdef BL_FORT_USE_UPPERCASE
#define FORT_GSRB     GSRB
#define FORT_ADOTX    ADOTX
#define FORT_DIAGSUM  DIAGSUM
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#define FORT_GSRB     gsrb
#define FORT_ADOTX    adotx
#define FORT_DIAGSUM  diagsum
#  else
#define FORT_GSRB     gsrb_
#define FORT_ADOTX    adotx_
#define FORT_DIAGSUM  diagsum_
#endif
#endif

#include <ArrayLim.H>

extern "C"
{

    void FORT_GSRB (
        int* isweep,
	int* num_sweeps,
        Real* offdiag_coeff, 
        int* check_for_singular, 
        Real* diag_regularization, 
        Real* phi,ARLIM_P(philo), ARLIM_P(phihi),
        const Real* rhs, ARLIM_P(rhslo), ARLIM_P(rhshi),
        const Real* diagnonsing, 
	ARLIM_P(diagnonsinglo),  ARLIM_P(diagnonsinghi),
        const Real* diagsing, 
        const Real* bXleft,
        const Real* bXright,
        const Real* bYleft,
        const Real* bYright,
        const Real* bZleft,
        const Real* bZright,
        const Real* icbx,
        const Real* icby,
        const Real* icbz,
        const Real* icdiag,
        const Real* icdiagrb,
        const Real* mask,
        const Real* ax,
        const Real* solnsave,
        const Real* rhssave,
        const Real* redsoln,
        const Real* blacksoln,
        const int *tilelo, const int *tilehi,
        const int *fablo, const int *fabhi,
        int* bfact,int* bfact_top,
        const int* smooth_type);

    void FORT_ADOTX(
        Real* offdiag_coeff, 
        int* check_for_singular, 
        Real* diag_regularization, 
        Real *y, ARLIM_P(ylo), ARLIM_P(yhi),
        const Real *x, ARLIM_P(xlo), ARLIM_P(xhi),
        const Real *diagnonsing,
	ARLIM_P(diagnonsinglo),ARLIM_P(diagnonsinghi),
        const Real *diagsing,
        const Real* bXleft,
        const Real* bXright,
        const Real* bYleft,
        const Real* bYright,
        const Real* bZleft,
        const Real* bZright,
        const int *tilelo, const int *tilehi,
        const int *fablo, const int *fabhi,
        int* bfact,int* bfact_top);
   
    void FORT_DIAGSUM(
        Real *y      , ARLIM_P(y_lo), ARLIM_P(y_hi),
        const Real* a , ARLIM_P(a_lo),  ARLIM_P(a_hi),
        const Real* bX, ARLIM_P(bX_lo), ARLIM_P(bX_hi),
        const Real* bY, ARLIM_P(bY_lo), ARLIM_P(bY_hi),
        const Real* bZ, ARLIM_P(bZ_lo), ARLIM_P(bZ_hi),
        const int *tilelo, const int *tilehi,
        const int *fablo, const int *fabhi,
        int* bfact,int* bfact_top);

}
#endif

#endif /*_AABBEC_F_H_*/
