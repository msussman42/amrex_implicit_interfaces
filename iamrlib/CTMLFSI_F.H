#ifndef _CTMLFSI_F_H_
#define _CTMLFSI_F_H_ 

#ifdef BL_LANG_FORT
#    define FORT_CTMLTRANSFERFORCE ctmltransferforce
#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_CTMLTRANSFERFORCE CTMLTRANSFERFORCE
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_CTMLTRANSFERFORCE ctmltransferforce
#  else
#    define FORT_CTMLTRANSFERFORCE ctmltransferforce_
#  endif
#endif

#include <ArrayLim.H>

extern "C" {
			
 void FORT_CTMLTRANSFERFORCE( 
  const int* tilelo, const int* tilehi,
  const int* fablo, const int* fabhi, 
  Real* velfab, ARLIM_P(velfablo), ARLIM_P(velfabhi),
  Real* forcefab, ARLIM_P(forcefablo), ARLIM_P(forcefabhi));

};
#endif

#endif
