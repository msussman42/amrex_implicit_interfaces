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

#ifndef _MAKESLICE_F_H_
#define _MAKESLICE_F_H_

/*
** $Id: MAKESLICE_F.H,v 1.6 1999/04/05 17:39:21 car Exp $
*/

#ifdef BL_LANG_FORT

#if (BL_SPACEDIM==3)
#    define FORT_MAKESLICE      makeslice
#    define FORT_MAKEMASKSLICE  makemaskslice
#endif
#endif

#ifdef BL_LANG_CC

#  ifdef BL_FORT_USE_UPPERCASE
#if (BL_SPACEDIM==3)
#    define FORT_MAKESLICE      MAKESLICE
#    define FORT_MAKEMASKSLICE  MAKEMASKSLICE
#endif

#  else
#  ifdef BL_FORT_USE_LOWERCASE
#if (BL_SPACEDIM==3)
#    define FORT_MAKESLICE      makeslice
#    define FORT_MAKEMASKSLICE  makemaskslice
#endif
#  else
#if (BL_SPACEDIM==3)
#    define FORT_MAKESLICE      makeslice_
#    define FORT_MAKEMASKSLICE  makemaskslice_
#endif
#  endif
#  endif

extern "C"
{
#if (BL_SPACEDIM==3)
    void FORT_MAKESLICE(const Real* data1,ARLIM_P(data1l), ARLIM_P(data1h),
                        Real* data2,ARLIM_P(data2l), ARLIM_P(data2h),
                        Real* data_min, Real* data_max,Real* slice_val,
                        int* slice,const Real* dx);

    void FORT_MAKEMASKSLICE(const int* data1,ARLIM_P(data1l), ARLIM_P(data1h),
                            int* data2,ARLIM_P(data2l), ARLIM_P(data2h),
                            Real* data_min, Real* data_max,Real* slice_val,
                            int* slice,const Real* dx);
#endif
};
#endif

#endif /*_MAKESLICE_F_H_*/
