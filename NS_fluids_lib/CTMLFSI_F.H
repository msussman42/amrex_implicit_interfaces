#ifndef _CTMLFSI_F_H_
#define _CTMLFSI_F_H_ 

#ifdef BL_LANG_FORT

#else

#include <AMReX_ArrayLim.H>

namespace amrex{

extern "C" {
			
 void fort_ctmltransferforce( 
  const int* tilelo, const int* tilehi,
  const int* fablo, const int* fabhi, 
  Real* velfab,AMREX_ARLIM_P(velfablo),AMREX_ARLIM_P(velfabhi),
  Real* forcefab,AMREX_ARLIM_P(forcefablo),AMREX_ARLIM_P(forcefabhi));

};

}/* namespace amrex */

#endif

#endif
