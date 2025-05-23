#ifndef _PROB_F_H_
#define _PROB_F_H_

#ifdef BL_LANG_FORT

#else

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
      const Real* dx,const Real* xlo, const Real* xhi);

     void fort_initdata_alloc(
      int* nc,
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
      Real* saturation_temp,
      Real* scal,
      AMREX_ARLIM_P(scallo),
      AMREX_ARLIM_P(scalhi),
      Real* refineden,
      AMREX_ARLIM_P(refinedenlo),
      AMREX_ARLIM_P(refinedenhi),
      Real* LS,
      AMREX_ARLIM_P(LSlo),
      AMREX_ARLIM_P(LShi),
      const Real* dx,const Real* xlo, const Real* xhi,
      Real* centroid_noise_factor);

     void fort_addnoise(
      int* dir,
      Real* angular_velocity_vector, //parameter to fort_addnoise
      int* perturbation_mode,
      Real* perturbation_eps_temp,
      Real* perturbation_eps_vel,
      int* nstate,
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
      int* verbose,
      int* ioproc,
      const int* max_level,
      int* bfact_space_level,
      int* bfact_grid_level,
      const int* domlo,const int* domhi,
      const Real* dx,
      Real* problo, Real* probhi);

     void fort_init_regions_list(
	int* constant_density_all_time,
	int* nthreads);
  
     void fort_delete_regions_list(int* ioproc);

     void fort_initdatasolid(
      int* nparts,
      int* im_solid_map,
      Real* time,
      const int* tilelo,const int* tilehi,
      const int* fablo,const int* fabhi,
      int* bfact,
      Real* solid,AMREX_ARLIM_P(solidlo),AMREX_ARLIM_P(solidhi),
      Real* LSnew,AMREX_ARLIM_P(LSnewlo),AMREX_ARLIM_P(LSnewhi),
      Real* snew,AMREX_ARLIM_P(snewlo),AMREX_ARLIM_P(snewhi),
      const Real* dx,const Real* xlo, const Real* xhi);

     void fort_initsolidtemp(
      int* nden,
      Real* time,
      const int* tilelo,const int* tilehi,
      const int* fablo,const int* fabhi,
      int* bfact,
      Real* snew,AMREX_ARLIM_P(snewlo),AMREX_ARLIM_P(snewhi),
      Real* lsnew,AMREX_ARLIM_P(lsnewlo),AMREX_ARLIM_P(lsnewhi),
      const Real* dx,const Real* xlo);

     void fort_vfracerror(
      int* tid_current,
      int* error_set_count,
      int* tag, 
      AMREX_ARLIM_P(tlo), AMREX_ARLIM_P(thi),
      const int* tagval, const int* clearval,
      Real* errfab, 
      AMREX_ARLIM_P(err_lo), AMREX_ARLIM_P(err_hi),
      Real* Snew,AMREX_ARLIM_P(Snewlo),AMREX_ARLIM_P(Snewhi), 
      Real* LSnew,AMREX_ARLIM_P(LSnewlo),AMREX_ARLIM_P(LSnewhi), 
      const int* tilelo, const int* tilehi, 
      const int* fablo, const int* fabhi,
      int* bfact, 
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
      int* tid_in,
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_refine_densityfill(
      int* tid_in,
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_pressurefill(
      int* tid_in,
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact); 

     void fort_group_statefill(
      int* tid_in,
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);


     void fort_group_refine_densityfill(
      int* tid_in,
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_group_tensorfill(
      int* tid_in,
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_scalarfill(
      int* tid_in,
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_extrapfill (
      int* tid_in,
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_group_extrapfill (
      int* tid_in,
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_ls_fill(
      int* tid_in,
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_group_ls_fill(
      int* tid_in,
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_moffill(
      int* tid_in,
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_extmoffill(
      int* tid_in,
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_group_moffill(
      int* tid_in,
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_group_extmoffill(
      int* tid_in,
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_velfill(
      int* tid_in,
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_group_velfill(
      int* tid_in,
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);


     void fort_solvfill(
      int* tid_in,
      int* grid_type,
      int* level,
      Real* data, AMREX_ARLIM_P(data_lo), AMREX_ARLIM_P(data_hi),
      const int* dom_lo, const int* dom_hi,
      const Real* dx, const Real* grd_lo,
      const Real* time, const int* bc,
      const int* scomp,int* ncomp,int* bfact);

     void fort_group_solvfill(
      int* tid_in,
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
      int* solidheat_flag,
      int* project_option,
      Real* time);

     void fort_umacfill(
      int* tid_in,
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
