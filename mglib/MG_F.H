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

#ifndef _MG_F_H_
#define _MG_F_H_

/*
** $Id: MG_F.H,v 1.2 1999/08/25 18:04:10 propp Exp $
*/

#include <REAL.H>

#ifdef BL_LANG_FORT
#define FORT_AVERAGE   average_for_mg
#define FORT_INTERP    interp_for_mg
#else
#  ifdef BL_FORT_USE_UPPERCASE
#define FORT_AVERAGE   AVERAGE_FOR_MG
#define FORT_INTERP    INTERP_FOR_MG
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#define FORT_AVERAGE   average_for_mg
#define FORT_INTERP    interp_for_mg
#  else
#define FORT_AVERAGE   average_for_mg_
#define FORT_INTERP    interp_for_mg_
#endif
#endif

#include <ArrayLim.H>

extern "C"
{
    void FORT_AVERAGE (
        Real* crse,       ARLIM_P(crse_lo), ARLIM_P(crse_hi),
        const Real* fine, ARLIM_P(fine_lo), ARLIM_P(fine_hi),
        const int *tlo, const int *thi,
        int* iaverage, 
        int* bfact_coarse,
        int* bfact_fine,
        int* bfact_top);

    void FORT_INTERP (
        int* bfact,int* bfact_f,int* bfact_top,
        Real* fine,       ARLIM_P(fine_lo), ARLIM_P(fine_hi),
        const Real* crse, ARLIM_P(crse_lo), ARLIM_P(crse_hi),
        const int *tlo, const int *thi);
}
#endif

#endif /*_MG_F_H_*/
