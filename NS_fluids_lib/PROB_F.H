#ifndef _PROB_F_H_
#define _PROB_F_H_

#ifdef BL_LANG_FORT
#    define FORT_FLUSH_FORTRAN flush_fortran
#    define FORT_SET_PERIODIC_VAR set_periodic_var
#    define FORT_SETFORTSCALES setfortscales
#    define FORT_FORCEVELOCITY forcevelocity
#    define FORT_INITVELOCITY initvelocity
#    define FORT_INITDATA_ALLOC initdata_alloc
#    define FORT_INITDATA    initdata
#    define FORT_ADDNOISE    addnoise
#    define FORT_INITGRIDMAP initgridmap
#    define FORT_INITDATASOLID initdatasolid
#    define FORT_INITSOLIDTEMP initsolidtemp
#    define FORT_VFRACERROR  vfracerror
#    define FORT_OVERRIDELSBC overridelsbc
#    define FORT_OVERRIDEPBC overridepbc
#    define FORT_VISCFLUXFILL viscfluxfill
#    define FORT_UMACFILL    umacfill
#    define FORT_VMACFILL    vmacfill
#    define FORT_WMACFILL    wmacfill
#    define FORT_XDMACFILL   xdmacfill
#    define FORT_YDMACFILL   ydmacfill
#    define FORT_ZDMACFILL   zdmacfill
#    define FORT_X_EXTRAPFILL x_extrapfill
#    define FORT_Y_EXTRAPFILL y_extrapfill
#    define FORT_Z_EXTRAPFILL z_extrapfill
#    define FORT_VELFILL     velfill
#    define FORT_GROUP_VELFILL group_velfill
#    define FORT_SOLVFILL    solvfill
#    define FORT_GROUP_SOLVFILL group_solvfill
#    define FORT_PRESSUREFILL  pressurefill
#    define FORT_STATEFILL   statefill
#    define FORT_GROUP_STATEFILL group_statefill
#    define FORT_TENSORFILL   tensorfill
#    define FORT_GROUP_TENSORFILL group_tensorfill
#    define FORT_XDISPLACEFILL   xdisplacefill
#    define FORT_SCALARFILL  scalarfill
#    define FORT_EXTRAPFILL  extrapfill
#    define FORT_GROUP_EXTRAPFILL group_extrapfill
#    define FORT_MOFFILL     moffill
#    define FORT_EXTMOFFILL  extmoffill
#    define FORT_GROUP_MOFFILL group_moffill
#    define FORT_GROUP_EXTMOFFILL group_extmoffill
#    define FORT_LS_HO_FILL ls_ho_fill
#    define FORT_GROUP_LS_HO_FILL group_ls_ho_fill
#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_FLUSH_FORTRAN FLUSH_FORTRAN
#    define FORT_SET_PERIODIC_VAR SET_PERIODIC_VAR
#    define FORT_SETFORTSCALES SETFORTSCALES
#    define FORT_FORCEVELOCITY FORCEVELOCITY
#    define FORT_INITVELOCITY INITVELOCITY
#    define FORT_INITGRIDMAP INITGRIDMAP
#    define FORT_INITDATA_ALLOC INITDATA_ALLOC
#    define FORT_INITDATA    INITDATA
#    define FORT_ADDNOISE    ADDNOISE
#    define FORT_SUPERELASTICOVERRIDE SUPERELASTICOVERRIDE
#    define FORT_INITDATASOLID INITDATASOLID
#    define FORT_INITSOLIDTEMP INITSOLIDTEMP
#    define FORT_VFRACERROR  VFRACERROR
#    define FORT_OVERRIDELSBC OVERRIDELSBC
#    define FORT_OVERRIDEPBC OVERRIDEPBC
#    define FORT_VISCFLUXFILL VISCFLUXFILL
#    define FORT_UMACFILL    UMACFILL
#    define FORT_VMACFILL    VMACFILL
#    define FORT_WMACFILL    WMACFILL
#    define FORT_XDMACFILL   XDMACFILL
#    define FORT_YDMACFILL   YDMACFILL
#    define FORT_ZDMACFILL   ZDMACFILL
#    define FORT_X_EXTRAPFILL X_EXTRAPFILL
#    define FORT_Y_EXTRAPFILL Y_EXTRAPFILL
#    define FORT_Z_EXTRAPFILL Z_EXTRAPFILL
#    define FORT_VELFILL     VELFILL
#    define FORT_GROUP_VELFILL GROUP_VELFILL
#    define FORT_SOLVFILL    SOLVFILL
#    define FORT_GROUP_SOLVFILL GROUP_SOLVFILL
#    define FORT_PRESSUREFILL  PRESSUREFILL
#    define FORT_STATEFILL  STATEFILL 
#    define FORT_GROUP_STATEFILL GROUP_STATEFILL
#    define FORT_TENSORFILL   TENSORFILL
#    define FORT_GROUP_TENSORFILL GROUP_TENSORFILL
#    define FORT_XDISPLACEFILL   XDISPLACEFILL
#    define FORT_SCALARFILL  SCALARFILL
#    define FORT_EXTRAPFILL  EXTRAPFILL
#    define FORT_GROUP_EXTRAPFILL GROUP_EXTRAPFILL
#    define FORT_MOFFILL     MOFFILL
#    define FORT_EXTMOFFILL  EXTMOFFILL
#    define FORT_GROUP_MOFFILL GROUP_MOFFILL
#    define FORT_GROUP_EXTMOFFILL GROUP_EXTMOFFILL
#    define FORT_LS_HO_FILL LS_HO_FILL
#    define FORT_GROUP_LS_HO_FILL GROUP_LS_HO_FILL
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_FLUSH_FORTRAN flush_fortran
#    define FORT_SET_PERIODIC_VAR set_periodic_var
#    define FORT_SETFORTSCALES setfortscales
#    define FORT_FORCEVELOCITY forcevelocity
#    define FORT_INITVELOCITY initvelocity
#    define FORT_INITGRIDMAP initgridmap
#    define FORT_INITDATA_ALLOC initdata_alloc
#    define FORT_INITDATA    initdata
#    define FORT_ADDNOISE    addnoise
#    define FORT_INITDATASOLID initdatasolid
#    define FORT_INITSOLIDTEMP initsolidtemp
#    define FORT_VFRACERROR  vfracerror
#    define FORT_OVERRIDELSBC overridelsbc
#    define FORT_OVERRIDEPBC overridepbc
#    define FORT_VISCFLUXFILL viscfluxfill
#    define FORT_UMACFILL    umacfill
#    define FORT_VMACFILL    vmacfill
#    define FORT_WMACFILL    wmacfill
#    define FORT_XDMACFILL   xdmacfill
#    define FORT_YDMACFILL   ydmacfill
#    define FORT_ZDMACFILL   zdmacfill
#    define FORT_X_EXTRAPFILL x_extrapfill
#    define FORT_Y_EXTRAPFILL y_extrapfill
#    define FORT_Z_EXTRAPFILL z_extrapfill
#    define FORT_VELFILL     velfill
#    define FORT_GROUP_VELFILL group_velfill
#    define FORT_SOLVFILL    solvfill
#    define FORT_GROUP_SOLVFILL group_solvfill
#    define FORT_PRESSUREFILL  pressurefill
#    define FORT_STATEFILL  statefill
#    define FORT_GROUP_STATEFILL group_statefill
#    define FORT_TENSORFILL   tensorfill
#    define FORT_GROUP_TENSORFILL group_tensorfill
#    define FORT_XDISPLACEFILL   xdisplacefill
#    define FORT_SCALARFILL  scalarfill
#    define FORT_EXTRAPFILL  extrapfill
#    define FORT_GROUP_EXTRAPFILL group_extrapfill
#    define FORT_MOFFILL     moffill
#    define FORT_EXTMOFFILL  extmoffill
#    define FORT_GROUP_MOFFILL group_moffill
#    define FORT_GROUP_EXTMOFFILL group_extmoffill
#    define FORT_LS_HO_FILL ls_ho_fill
#    define FORT_GROUP_LS_HO_FILL group_ls_ho_fill
#  else
#    define FORT_FLUSH_FORTRAN flush_fortran_
#    define FORT_SET_PERIODIC_VAR set_periodic_var_
#    define FORT_SETFORTSCALES setfortscales_
#    define FORT_FORCEVELOCITY forcevelocity_
#    define FORT_INITVELOCITY initvelocity_
#    define FORT_INITGRIDMAP initgridmap_
#    define FORT_INITDATA_ALLOC initdata_alloc_
#    define FORT_INITDATA    initdata_
#    define FORT_ADDNOISE    addnoise_
#    define FORT_INITDATASOLID initdatasolid_
#    define FORT_INITSOLIDTEMP initsolidtemp_
#    define FORT_VFRACERROR  vfracerror_
#    define FORT_OVERRIDELSBC overridelsbc_
#    define FORT_OVERRIDEPBC overridepbc_
#    define FORT_VISCFLUXFILL viscfluxfill_
#    define FORT_UMACFILL    umacfill_
#    define FORT_VMACFILL    vmacfill_
#    define FORT_WMACFILL    wmacfill_
#    define FORT_XDMACFILL   xdmacfill_
#    define FORT_YDMACFILL   ydmacfill_
#    define FORT_ZDMACFILL   zdmacfill_
#    define FORT_X_EXTRAPFILL x_extrapfill_
#    define FORT_Y_EXTRAPFILL y_extrapfill_
#    define FORT_Z_EXTRAPFILL z_extrapfill_
#    define FORT_VELFILL     velfill_
#    define FORT_GROUP_VELFILL group_velfill_
#    define FORT_SOLVFILL    solvfill_
#    define FORT_GROUP_SOLVFILL group_solvfill_
#    define FORT_PRESSUREFILL  pressurefill_
#    define FORT_STATEFILL  statefill_
#    define FORT_GROUP_STATEFILL group_statefill_
#    define FORT_TENSORFILL   tensorfill_
#    define FORT_GROUP_TENSORFILL group_tensorfill_
#    define FORT_XDISPLACEFILL   xdisplacefill_
#    define FORT_SCALARFILL  scalarfill_
#    define FORT_EXTRAPFILL  extrapfill_
#    define FORT_GROUP_EXTRAPFILL group_extrapfill_
#    define FORT_MOFFILL     moffill_
#    define FORT_EXTMOFFILL  extmoffill_
#    define FORT_GROUP_MOFFILL group_moffill_
#    define FORT_GROUP_EXTMOFFILL group_extmoffill_
#    define FORT_LS_HO_FILL ls_ho_fill_
#    define FORT_GROUP_LS_HO_FILL group_ls_ho_fill_
#  endif
#  endif

#include "AMReX_ArrayLim.H"

namespace amrex{

  extern "C" {
     void FORT_FLUSH_FORTRAN();

     void FORT_SET_PERIODIC_VAR(int* periodic_flag);

     void FORT_SETFORTSCALES(Real* pressure_scale,
       Real* velocity_scale);


     void FORT_OVERRIDELSBC(int* homflag);
     void FORT_OVERRIDEPBC(int* homflag,int* project_option);

     void FORT_INITVELOCITY (
      const int* level,const Real* time,
      const int* tilelo,const int* tilehi,
      const int* fablo,const int* fabhi,
      int* bfact,
      Real* vel,AMREX_ARLIM_P(vello),AMREX_ARLIM_P(velhi),
      const Real* dx,const Real* xlo, const Real* xhi,
      Real* Re,Real* We,Real* RGASRWATER,int* use_lsa);

     void FORT_INITDATA_ALLOC(
      int* nmat,int* nten,int* nc,
      Real* latent_heat,
      int* freezing_model,
      int* distribute_from_target,
      Real* saturation_temp,
      const Real* dx);

     void FORT_INITDATA (
      int* tid,
      int* adapt_quad_depth,
      int* level,int* max_level,
      Real* time,
      const int* tilelo,const int* tilehi,
      const int* fablo,const int* fabhi,
      int* bfact,
      int* nc,
      int* nmat,int* nten,
      Real* latent_heat,
      Real* saturation_temp,
      int* radius_cutoff,
      Real* scal,AMREX_ARLIM_P(scallo),AMREX_ARLIM_P(scalhi),
      Real* LS,AMREX_ARLIM_P(LSlo),AMREX_ARLIM_P(LShi),
      const Real* dx,const Real* xlo, const Real* xhi);

     void FORT_ADDNOISE(
      int* dir,
      Real* angular_velocity,
      int* perturbation_mode,
      Real* perturbation_eps_temp,
      Real* perturbation_eps_vel,
      int* nstate,
      int* nmat,
      const Real* xlo,const Real* dx,
      Real* Snew,AMREX_ARLIM_P(Snewlo),AMREX_ARLIM_P(Snewhi), 
      Real* LSnew,AMREX_ARLIM_P(LSnewlo),AMREX_ARLIM_P(LSnewhi), 
      Real* MAC,AMREX_ARLIM_P(MAClo),AMREX_ARLIM_P(MAChi), 
      const int* tilelo,const int* tilehi,
      const int* fablo,const int* fabhi,
      int* bfact,
      int* level,
      int* finest_level);

     void FORT_INITGRIDMAP (
      const int* max_level,
      int* bfact_space_level,
      int* bfact_grid_level,
      const int* domlo,const int* domhi,
      const Real* dx,
      Real* problo, Real* probhi);


     void FORT_INITDATASOLID (
      int* nmat,
      int* nparts,
      int* nFSI_sub,
      int* nFSI,
      int* ngrowFSI,
      int* im_solid_map,
      Real* time,
      const int* tilelo,const int* tilehi,
      const int* fablo,const int* fabhi,
      int* bfact,
      Real* solid,AMREX_ARLIM_P(solidlo),AMREX_ARLIM_P(solidhi),
      const Real* dx,const Real* xlo, const Real* xhi);

     void FORT_INITSOLIDTEMP (
      int* nmat,
      int* nden,
      Real* time,
      const int* tilelo,const int* tilehi,
      const int* fablo,const int* fabhi,
      int* bfact,
      Real* snew,AMREX_ARLIM_P(snewlo),AMREX_ARLIM_P(snewhi),
      Real* lsnew,AMREX_ARLIM_P(lsnewlo),AMREX_ARLIM_P(lsnewhi),
      const Real* dx,const Real* xlo);

     void FORT_VFRACERROR  (
      int* tag, AMREX_ARLIM_P(tlo), AMREX_ARLIM_P(thi),
      const int* tagval, const int* clearval,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* tilelo, const int* tilehi, 
      const int* fablo, const int* fabhi,int* bfact, 
      const int* nvar,
      const int* domain_lo, const int* domain_hi,
      const Real* dx, const Real* xlo,
      const Real* prob_lo, 
      const Real* time,
      const int* level,
      const int* max_level,
      const int* max_level_two_materials,
      int* nblocks,
      Real* xblocks,Real* yblocks,Real* zblocks,
      Real* rxblocks,Real* ryblocks,Real* rzblocks,
      int* ncoarseblocks,
      Real* xcoarseblocks,Real* ycoarseblocks,Real* zcoarseblocks,
      Real* rxcoarseblocks,Real* rycoarseblocks,Real* rzcoarseblocks);

     void FORT_FORCEVELOCITY(
      int* nsolveMM_FACE,
      const Real* prob_lo,const Real* prob_hi,
      Real* vel, AMREX_ARLIM_P(vel_lo), AMREX_ARLIM_P(vel_hi),
      Real* velmac, AMREX_ARLIM_P(velmac_lo), AMREX_ARLIM_P(velmac_hi),
      const int* dir,
      const Real* xlo,const Real* dx,
      const int* tilelo, const int* tilehi,
      const int* fablo, const int* fabhi,
      int* bfact,
      Real* time,
      int* presbc,
      Real* outflow_velocity_buffer_size);

      // grd_lo=problo+dx*(datalo-domlo)
      // x=grd_lo+(i-datalo+1/2)dx=problo+dx*(i-domlo+1/2)

     void FORT_STATEFILL (
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void FORT_TENSORFILL (
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void FORT_XDISPLACEFILL (
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void FORT_PRESSUREFILL (
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact); 

     void FORT_GROUP_STATEFILL (
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void FORT_GROUP_TENSORFILL (
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void FORT_SCALARFILL (
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void FORT_EXTRAPFILL (
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void FORT_GROUP_EXTRAPFILL (
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void FORT_LS_HO_FILL (
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void FORT_GROUP_LS_HO_FILL (
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void FORT_MOFFILL (
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void FORT_EXTMOFFILL (
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void FORT_GROUP_MOFFILL (
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void FORT_GROUP_EXTMOFFILL (
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void FORT_VELFILL (
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void FORT_GROUP_VELFILL (
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);


     void FORT_SOLVFILL (
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void FORT_GROUP_SOLVFILL (
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void FORT_VISCFLUXFILL(
      Real* macrolayer_size,
      int* microlayer_substrate,
      Real* microlayer_temperature_substrate,
      Real* latent_heat,
      int* freezing_model,
      Real* saturation_temp,
      int* nsolveMM_FACE,
      int* nsolve,
      int* nsolveMM,
      int* dir,
      const Real* xlo,const Real* dx,
      int* velbc,
      int* tempbc,
      int* domtempbc,
      Real* LS, 
      AMREX_ARLIM_P(LSlo), AMREX_ARLIM_P(LShi),
      Real* area, 
      AMREX_ARLIM_P(arealo), AMREX_ARLIM_P(areahi),
      Real* xflux, 
      AMREX_ARLIM_P(xfluxlo), AMREX_ARLIM_P(xfluxhi),
      const int* tilelo,const int* tilehi,
      const int* fablo,const int* fabhi,
      int* bfact,
      const int* domlo,const int* domhi,
      Real* dt,
      int* nmat,
      int* nten,
      int* solidheat_flag,
      int* project_option,
      Real* time);

     void FORT_UMACFILL (
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);


     void FORT_VMACFILL (
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void FORT_WMACFILL (
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

  }

}/* namespace amrex */

#endif
#endif 
