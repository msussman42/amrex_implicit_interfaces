#ifndef _MOF_REDIST_F_H_
#define _MOF_REDIST_F_H_ 

#ifdef BL_LANG_FORT
#    define FORT_FD_NORMAL fd_normal
#    define FORT_LEVELSTRIP levelstrip
#    define FORT_STENINIT steninit
#    define FORT_FACEINITTEST faceinittest
#    define FORT_FACEINIT faceinit
#    define FORT_CORRECT_UNINIT correct_uninit
#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_FD_NORMAL FD_NORMAL
#    define FORT_LEVELSTRIP LEVELSTRIP
#    define FORT_STENINIT STENINIT
#    define FORT_FACEINITTEST FACEINITTEST
#    define FORT_FACEINIT FACEINIT
#    define FORT_CORRECT_UNINIT CORRECT_UNINIT
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_FD_NORMAL fd_normal
#    define FORT_LEVELSTRIP levelstrip
#    define FORT_STENINIT steninit
#    define FORT_FACEINITTEST faceinittest
#    define FORT_FACEINIT faceinit
#    define FORT_CORRECT_UNINIT correct_uninit
#  else
#    define FORT_FD_NORMAL fd_normal_
#    define FORT_LEVELSTRIP levelstrip_
#    define FORT_STENINIT steninit_
#    define FORT_FACEINITTEST faceinittest_
#    define FORT_FACEINIT faceinit_
#    define FORT_CORRECT_UNINIT correct_uninit_
#  endif
#  endif

#include <ArrayLim.H>

  extern "C" {


     void FORT_FD_NORMAL(
       int* level,
       int* finest_level,
       Real* lsnew,ARLIM_P(lsnewlo),ARLIM_P(lsnewhi),
       Real* lsnrmFD,ARLIM_P(lsnrmFDlo),ARLIM_P(lsnrmFDhi),
       const int* tilelo,const int* tilehi,
       const int* fablo,const int* fabhi,
       int* bfact,
       const Real* xlo,const Real* dx,
       int* nmat);

     void FORT_LEVELSTRIP(
       int* nprocessed,
       Real* minLS,
       Real* maxLS,
       Real* max_problen,
       int* level,
       int* finest_level,
       Real* latent_heat,
       Real* maskfab,ARLIM_P(maskfablo),ARLIM_P(maskfabhi),
       Real* facefab,ARLIM_P(facefablo),ARLIM_P(facefabhi),
       Real* facetest,ARLIM_P(facetestlo),ARLIM_P(facetesthi),
       Real* stenfab,ARLIM_P(stenfablo),ARLIM_P(stenfabhi),
       Real* vofrecon,ARLIM_P(vofreconlo),ARLIM_P(vofreconhi),
       Real* origdist,ARLIM_P(origdistlo),ARLIM_P(origdisthi),
       Real* newfab,ARLIM_P(newfablo),ARLIM_P(newfabhi),
       Real* touchfab,ARLIM_P(touchfablo),ARLIM_P(touchfabhi),
       Real* crsetouch,ARLIM_P(crsetouchlo),ARLIM_P(crsetouchhi),
       Real* crsedist,ARLIM_P(crsedistlo),ARLIM_P(crsedisthi),
       const int* tilelo,const int* tilehi,
       const int* fablo,const int* fabhi,
       int* bfact,
       const int* bc,const int* rz_flag,
       const Real* xlo,const Real* dx,
       Real* time,
       int* ngrowLS,
       int* nmat,int* nten,int* nstar,int* nface);

     void FORT_CORRECT_UNINIT(
       Real* minLS,
       Real* maxLS,
       Real* max_problen,
       int* level,
       int* finest_level,
       Real* newfab,ARLIM_P(newfablo),ARLIM_P(newfabhi),
       Real* touchfab,ARLIM_P(touchfablo),ARLIM_P(touchfabhi),
       const int* tilelo,const int* tilehi,
       const int* fablo,const int* fabhi,
       int* bfact,
       const Real* xlo,const Real* dx,
       Real* time,
       int* nmat);

     void FORT_STENINIT(
       int* level,int* finest_level,
       Real* stenfab,ARLIM_P(stenfablo),ARLIM_P(stenfabhi),
       Real* mask,ARLIM_P(masklo),ARLIM_P(maskhi),
       Real* vofrecon,ARLIM_P(vofreconlo),ARLIM_P(vofreconhi),
       const int* tilelo,const int* tilehi,
       const int* fablo,const int* fabhi,
       int* bfact,
       const int* rz_flag,
       const Real* xlo,const Real* dx,
       Real* time,int* ngrowLS,
       int* nmat,int* nstar);


     void FORT_FACEINIT(
       int* tid,
       int* tessellate,
       int* nten,
       int* level,
       int* finest_level,
       Real* facefab,ARLIM_P(facefablo),ARLIM_P(facefabhi),
       Real* maskfab,ARLIM_P(maskfablo),ARLIM_P(maskfabhi),
       Real* vofrecon,ARLIM_P(vofreconlo),ARLIM_P(vofreconhi),
       const int* tilelo,const int* tilehi,
       const int* fablo,const int* fabhi,
       int* bfact,
       const int* rz_flag,
       const Real* xlo,const Real* dx,
       Real* time,
       int* ngrow,
       int* nmat,
       int* nface,
       int* nface_decomp);

     void FORT_FACEINITTEST(
       int* tid,
       int* tessellate,
       int* level,
       int* finest_level,
       Real* facefab,ARLIM_P(facefablo),ARLIM_P(facefabhi),
       Real* facetest,ARLIM_P(facetestlo),ARLIM_P(facetesthi),
       Real* maskfab,ARLIM_P(maskfablo),ARLIM_P(maskfabhi),
       Real* vofrecon,ARLIM_P(vofreconlo),ARLIM_P(vofreconhi),
       const int* tilelo,const int* tilehi,
       const int* fablo,const int* fabhi,
       int* bfact,
       const int* rz_flag,
       const Real* xlo,const Real* dx,
       Real* time,
       int* ngrow,
       int* nmat,
       int* nface);

  };
#endif

#endif

