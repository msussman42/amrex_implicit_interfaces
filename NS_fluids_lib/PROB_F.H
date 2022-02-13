#ifndef _PROB_F_H_
#define _PROB_F_H_

#ifdef BL_LANG_FORT
#    define FORT_PRESSUREFILL  pressurefill
#    define FORT_GROUP_STATEFILL group_statefill
#    define FORT_TENSORFILL   tensorfill
#    define FORT_GROUP_TENSORFILL group_tensorfill
#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_PRESSUREFILL  PRESSUREFILL
#    define FORT_GROUP_STATEFILL GROUP_STATEFILL
#    define FORT_TENSORFILL   TENSORFILL
#    define FORT_GROUP_TENSORFILL GROUP_TENSORFILL
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_PRESSUREFILL  pressurefill
#    define FORT_GROUP_STATEFILL group_statefill
#    define FORT_TENSORFILL   tensorfill
#    define FORT_GROUP_TENSORFILL group_tensorfill
#  else
#    define FORT_PRESSUREFILL  pressurefill_
#    define FORT_GROUP_STATEFILL group_statefill_
#    define FORT_TENSORFILL   tensorfill_
#    define FORT_GROUP_TENSORFILL group_tensorfill_
#  endif
#  endif

#include "AMReX_ArrayLim.H"

namespace amrex{

  extern "C" {
     void fort_flush_fortran();

     void fort_set_periodic_var(int* periodic_flag);

     void fort_setfortscales(Real* pressure_scale,
       Real* velocity_scale);


     void fort_overridelsbc(int* homflag);
     void fort_overridepbc(int* homflag,int* project_option);

     void fort_initvelocity(
      const int* level,const Real* time,
      const int* tilelo,const int* tilehi,
      const int* fablo,const int* fabhi,
      int* bfact,
      Real* vel,AMREX_ARLIM_P(vello),AMREX_ARLIM_P(velhi),
      const Real* dx,const Real* xlo, const Real* xhi,
      Real* Re,Real* We,Real* RGASRWATER);

     void fort_initdata_alloc(
      int* nmat,int* nten,int* nc,
      Real* latent_heat,
      int* freezing_model,
      int* distribute_from_target,
      Real* saturation_temp,
      const Real* dx);

     void fort_initdata(
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
      Real* scal,AMREX_ARLIM_P(scallo),AMREX_ARLIM_P(scalhi),
      Real* LS,AMREX_ARLIM_P(LSlo),AMREX_ARLIM_P(LShi),
      const Real* dx,const Real* xlo, const Real* xhi);

     void fort_addnoise(
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

     void fort_initgridmap(
      const int* max_level,
      int* bfact_space_level,
      int* bfact_grid_level,
      const int* domlo,const int* domhi,
      const Real* dx,
      Real* problo, Real* probhi);

     void fort_init_regions_list(
	int* constant_density_all_time,
	int* num_materials,
	int* nthreads);
  
     void fort_delete_regions_list(int* ioproc);

     void fort_initdatasolid(
      int* nmat,
      int* nparts,
      int* nFSI,
      int* ngrow_make_distance_in,
      int* im_solid_map,
      Real* time,
      const int* tilelo,const int* tilehi,
      const int* fablo,const int* fabhi,
      int* bfact,
      Real* solid,AMREX_ARLIM_P(solidlo),AMREX_ARLIM_P(solidhi),
      const Real* dx,const Real* xlo, const Real* xhi);

     void fort_initsolidtemp(
      int* nmat,
      int* nden,
      Real* time,
      const int* tilelo,const int* tilehi,
      const int* fablo,const int* fabhi,
      int* bfact,
      Real* snew,AMREX_ARLIM_P(snewlo),AMREX_ARLIM_P(snewhi),
      Real* lsnew,AMREX_ARLIM_P(lsnewlo),AMREX_ARLIM_P(lsnewhi),
      const Real* dx,const Real* xlo);

     void fort_vfracerror(
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
      const int* max_level_for_use,
      int* nblocks,
      Real* xblocks,Real* yblocks,Real* zblocks,
      Real* rxblocks,Real* ryblocks,Real* rzblocks,
      int* ncoarseblocks,
      Real* xcoarseblocks,Real* ycoarseblocks,Real* zcoarseblocks,
      Real* rxcoarseblocks,Real* rycoarseblocks,Real* rzcoarseblocks);

     void fort_forcevelocity(
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

     void fort_statefill(
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void FORT_TENSORFILL (
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void FORT_PRESSUREFILL (
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact); 

     void FORT_GROUP_STATEFILL (
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void FORT_GROUP_TENSORFILL (
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_scalarfill(
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_extrapfill (
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_group_extrapfill (
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_ls_ho_fill(
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_group_ls_ho_fill(
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_moffill(
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_extmoffill(
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_group_moffill(
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_group_extmoffill(
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_velfill(
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_group_velfill(
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);


     void fort_solvfill(
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_group_solvfill(
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_viscfluxfill(
      Real* macrolayer_size,
      int* microlayer_substrate,
      Real* microlayer_temperature_substrate,
      Real* latent_heat,
      int* freezing_model,
      Real* saturation_temp,
      int* nsolve,
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

     void fort_umacfill(
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_xdmacfill(
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_x_extrapfill (
      int* grid_type,
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
