#ifndef _MOF_REDIST_F_H_
#define _MOF_REDIST_F_H_ 

#ifdef BL_LANG_FORT

#else

#include <AMReX_ArrayLim.H>

namespace amrex{

  extern "C" {


     void fort_fd_normal(
       int* level,
       int* finest_level,
       Real* lsnew,AMREX_ARLIM_P(lsnewlo),AMREX_ARLIM_P(lsnewhi),
       Real* lsnrmFD,AMREX_ARLIM_P(lsnrmFDlo),AMREX_ARLIM_P(lsnrmFDhi),
       const int* tilelo,const int* tilehi,
       const int* fablo,const int* fabhi,
       int* bfact,
       const Real* xlo,const Real* dx);


     void fort_fd_node_normal(
       int* level,
       int* finest_level,
       Real* lsnew,AMREX_ARLIM_P(lsnewlo),AMREX_ARLIM_P(lsnewhi),
       Real* lsnrmFD,AMREX_ARLIM_P(lsnrmFDlo),AMREX_ARLIM_P(lsnrmFDhi),
       const int* tilelo,const int* tilehi,
       const int* fablo,const int* fabhi,
       int* bfact,
       const Real* xlo,const Real* dx,
       int* n_normal);

     void fort_node_to_cell(
       int* tid_current,
       int* level,
       int* finest_level,
       int* height_function_flag,
       Real* F_new,AMREX_ARLIM_P(F_new_lo),AMREX_ARLIM_P(F_new_hi),
       Real* lsnew,AMREX_ARLIM_P(lsnewlo),AMREX_ARLIM_P(lsnewhi),
       Real* lsnrmFD,AMREX_ARLIM_P(lsnrmFDlo),AMREX_ARLIM_P(lsnrmFDhi),
       Real* curvcell,AMREX_ARLIM_P(curvcelllo),AMREX_ARLIM_P(curvcellhi),
       const int* tilelo,const int* tilehi,
       const int* fablo,const int* fabhi,
       int* bfact,
       const Real* xlo,const Real* dx,
       int* n_normal);


     void fort_levelstrip(
       int* nprocessed,
       Real* minLS,
       Real* maxLS,
       Real* max_problen,
       int* level,
       int* finest_level,
       Real* maskfab,AMREX_ARLIM_P(maskfablo),AMREX_ARLIM_P(maskfabhi),
       Real* stenfab,AMREX_ARLIM_P(stenfablo),AMREX_ARLIM_P(stenfabhi),
       Real* vofrecon,AMREX_ARLIM_P(vofreconlo),AMREX_ARLIM_P(vofreconhi),
       Real* newfab,AMREX_ARLIM_P(newfablo),AMREX_ARLIM_P(newfabhi),
       Real* touchfab,AMREX_ARLIM_P(touchfablo),AMREX_ARLIM_P(touchfabhi),
       Real* crsetouch,AMREX_ARLIM_P(crsetouchlo),AMREX_ARLIM_P(crsetouchhi),
       Real* crsedist,AMREX_ARLIM_P(crsedistlo),AMREX_ARLIM_P(crsedisthi),
       const int* tilelo,const int* tilehi,
       const int* fablo,const int* fabhi,
       int* bfact,
       const int* bc,const int* rz_flag,
       const Real* xlo,const Real* dx,
       Real* time);

     void fort_correct_uninit(
       Real* minLS,
       Real* maxLS,
       Real* max_problen,
       int* level,
       int* finest_level,
       Real* newfab,AMREX_ARLIM_P(newfablo),AMREX_ARLIM_P(newfabhi),
       Real* touchfab,AMREX_ARLIM_P(touchfablo),AMREX_ARLIM_P(touchfabhi),
       const int* tilelo,const int* tilehi,
       const int* fablo,const int* fabhi,
       int* bfact,
       const Real* xlo,const Real* dx,
       Real* time);

     void fort_steninit(
       int* level,int* finest_level,
       Real* stenfab,AMREX_ARLIM_P(stenfablo),AMREX_ARLIM_P(stenfabhi),
       Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
       Real* vofrecon,AMREX_ARLIM_P(vofreconlo),AMREX_ARLIM_P(vofreconhi),
       const int* tilelo,const int* tilehi,
       const int* fablo,const int* fabhi,
       int* bfact,
       const int* rz_flag,
       const Real* xlo,const Real* dx,
       Real* time);


     void fort_faceinit(
       int* tid,
       int* tessellate,
       int* level,
       int* finest_level,
       Real* facefab,AMREX_ARLIM_P(facefablo),AMREX_ARLIM_P(facefabhi),
       Real* maskfab,AMREX_ARLIM_P(maskfablo),AMREX_ARLIM_P(maskfabhi),
       Real* vofrecon,AMREX_ARLIM_P(vofreconlo),AMREX_ARLIM_P(vofreconhi),
       const int* tilelo,const int* tilehi,
       const int* fablo,const int* fabhi,
       int* bfact,
       const int* rz_flag,
       const Real* xlo,const Real* dx,
       Real* time,
       int* ngrow,
       int* nface);

     void fort_faceprocess(
       int* ngrow_source,
       int* ngrow_dest,
       int* tid,
       int* dir,
       int* tessellate,
       int* level,
       int* finest_level,
       Real* dstfab,AMREX_ARLIM_P(dstfablo),AMREX_ARLIM_P(dstfabhi),
       Real* facefab,AMREX_ARLIM_P(facefablo),AMREX_ARLIM_P(facefabhi),
       Real* vofrecon,AMREX_ARLIM_P(vofreconlo),AMREX_ARLIM_P(vofreconhi),
       const int* tilelo,const int* tilehi,
       const int* fablo,const int* fabhi,
       int* bfact,
       const int* rz_flag,
       const Real* xlo,const Real* dx,
       Real* time,
       int* nface_src,int* nface_dst);

  }

}/* namespace amrex */

#endif

#endif

