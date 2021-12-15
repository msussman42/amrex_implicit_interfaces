#ifndef _NAVIERSTOKES_F_H_
#define _NAVIERSTOKES_F_H_

/*
** $Id: NAVIERSTOKES_F.H,v 1.10 2002/03/26 20:51:08 lijewski Exp $
*/

#ifdef BL_LANG_FORT
#    define FORT_PRESSURE_INDICATOR pressure_indicator
#    define FORT_ZALESAK_CELL zalesak_cell
#    define FORT_EDGEAVGDOWN edgeavgdown
#    define FORT_FLUIDSOLIDCOR    fluidsolidcor
#    define FORT_MOFAVGDOWN  mofavgdown
#    define FORT_ERRORAVGDOWN erroravgdown
#    define FORT_AVGDOWN     avgdown
#    define FORT_AVGDOWN_LOW avgdown_low
#    define FORT_AVGDOWN_TAG avgdown_tag
#    define FORT_AVGDOWN_CURV avgdown_curv
#    define FORT_DIAGINV     diaginv
#    define FORT_OUTPUTSLICE outputslice
#    define FORT_IO_COMPARE io_compare
#    define FORT_METRICS metrics
#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_PRESSURE_INDICATOR PRESSURE_INDICATOR
#    define FORT_ZALESAK_CELL ZALESAK_CELL
#    define FORT_EDGEAVGDOWN EDGEAVGDOWN
#    define FORT_VOFFLUX_CORRECT VOFFLUX_CORRECT
#    define FORT_FLUIDSOLIDCOR    FLUIDSOLIDCOR
#    define FORT_MOFAVGDOWN  MOFAVGDOWN
#    define FORT_ERRORAVGDOWN ERRORAVGDOWN
#    define FORT_AVGDOWN     AVGDOWN
#    define FORT_AVGDOWN_LOW AVGDOWN_LOW
#    define FORT_AVGDOWN_TAG AVGDOWN_TAG
#    define FORT_AVGDOWN_CURV AVGDOWN_CURV
#    define FORT_DIAGINV     DIAGINV
#    define FORT_OUTPUTSLICE OUTPUTSLICE
#    define FORT_IO_COMPARE IO_COMPARE
#    define FORT_METRICS METRICS
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_PRESSURE_INDICATOR pressure_indicator
#    define FORT_ZALESAK_CELL zalesak_cell
#    define FORT_EDGEAVGDOWN edgeavgdown
#    define FORT_FLUIDSOLIDCOR    fluidsolidcor
#    define FORT_MOFAVGDOWN  mofavgdown
#    define FORT_ERRORAVGDOWN erroravgdown
#    define FORT_AVGDOWN     avgdown
#    define FORT_AVGDOWN_LOW avgdown_low
#    define FORT_AVGDOWN_TAG avgdown_tag
#    define FORT_AVGDOWN_CURV avgdown_curv
#    define FORT_DIAGINV     diaginv
#    define FORT_OUTPUTSLICE outputslice
#    define FORT_IO_COMPARE io_compare
#    define FORT_METRICS metrics
#  else
#    define FORT_PRESSURE_INDICATOR pressure_indicator_
#    define FORT_ZALESAK_CELL zalesak_cell_
#    define FORT_EDGEAVGDOWN edgeavgdown_
#    define FORT_FLUIDSOLIDCOR    fluidsolidcor_
#    define FORT_MOFAVGDOWN  mofavgdown_
#    define FORT_ERRORAVGDOWN erroravgdown_
#    define FORT_AVGDOWN     avgdown_
#    define FORT_AVGDOWN_LOW avgdown_low_
#    define FORT_AVGDOWN_TAG avgdown_tag_
#    define FORT_AVGDOWN_CURV avgdown_curv_
#    define FORT_DIAGINV     diaginv_
#    define FORT_OUTPUTSLICE outputslice_
#    define FORT_IO_COMPARE io_compare_
#    define FORT_METRICS metrics_
#  endif
#  endif

#include <AMReX_ArrayLim.H>

namespace amrex{

extern "C" {

    void fort_advective_pressure(
     int* level,
     int* finest_level,
     const Real* xlo,const Real* dx,
     Real* dt,
     Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi),
     Real* vol,AMREX_ARLIM_P(vollo),AMREX_ARLIM_P(volhi),
     Real* lsnew,AMREX_ARLIM_P(lsnewlo),AMREX_ARLIM_P(lsnewhi),
     Real* csound,AMREX_ARLIM_P(csoundlo),AMREX_ARLIM_P(csoundhi),
     Real* cvof,AMREX_ARLIM_P(cvoflo),AMREX_ARLIM_P(cvofhi),
     Real* den,AMREX_ARLIM_P(denlo),AMREX_ARLIM_P(denhi),
     Real* mdot,AMREX_ARLIM_P(mdotlo),AMREX_ARLIM_P(mdothi),
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,
     int* bfact,
     int* nmat,int* nden,
     Real* compressible_dt_factor,
     int* pressure_select_criterion,
     int* project_option);

    void fort_update_div(
     const Real* xlo,const Real* dx,
     Real* dt,
     Real* vol,AMREX_ARLIM_P(vollo),AMREX_ARLIM_P(volhi),
     Real* csound,AMREX_ARLIM_P(csoundlo),AMREX_ARLIM_P(csoundhi),
     Real* mdot,AMREX_ARLIM_P(mdotlo),AMREX_ARLIM_P(mdothi),
     Real* pnew,AMREX_ARLIM_P(pnewlo),AMREX_ARLIM_P(pnewhi),
     Real* divnew,AMREX_ARLIM_P(divnewlo),AMREX_ARLIM_P(divnewhi),
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,
     int* bfact,
     int* nmat);


    void fort_eos_pressure(
     int* level,
     int* finest_level,
     int* local_material_type,
     const Real* xlo,const Real* dx,
     Real* pres,AMREX_ARLIM_P(preslo),AMREX_ARLIM_P(preshi),
     Real* recon,AMREX_ARLIM_P(reconlo),AMREX_ARLIM_P(reconhi),
     Real* levelpc,AMREX_ARLIM_P(levelpclo),AMREX_ARLIM_P(levelpchi),
     Real* den,AMREX_ARLIM_P(denlo),AMREX_ARLIM_P(denhi),
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,
     int* bfact,
     int* nmat,
     int* nden);


    void FORT_PRESSURE_INDICATOR(
     int* pressure_error_flag,
     Real* vorterr,
     Real* pressure_error_cutoff,
     Real* temperature_error_cutoff,
     const Real* xlo,const Real* dx,
     Real* errnew,AMREX_ARLIM_P(errnewlo),AMREX_ARLIM_P(errnewhi),
     Real* slrecon,AMREX_ARLIM_P(slreconlo),AMREX_ARLIM_P(slreconhi),
     Real* den,AMREX_ARLIM_P(denlo),AMREX_ARLIM_P(denhi),
     Real* vort,AMREX_ARLIM_P(vortlo),AMREX_ARLIM_P(vorthi),
     Real* pres,AMREX_ARLIM_P(preslo),AMREX_ARLIM_P(preshi),
     Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi),
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,
     int* bfact,
     int* level,
     int* finest_level,
     int* nmat);

    void fort_coflow(
     Real* time,
     const int* fdomlo,
     const int* fdomhi,
     int* Z_dir,  //0..sdim-1
     int* R_dir,  //0..sdim-1
     int* num_cells,
     Real* coflow_Z, 
     Real* coflow_R_of_Z);

    void fort_cellgrid(
     int* tid,
     int* bfact,
     Real* fabout,AMREX_ARLIM_P(faboutlo),AMREX_ARLIM_P(fabouthi),
     const int* vislo,const int* vishi,
     int* visual_ncomp,
     Real* maskSEM,AMREX_ARLIM_P(maskSEMlo),AMREX_ARLIM_P(maskSEMhi),
     Real* vel,AMREX_ARLIM_P(vello),AMREX_ARLIM_P(velhi),
     Real* vof,AMREX_ARLIM_P(voflo),AMREX_ARLIM_P(vofhi),
     Real* pres,AMREX_ARLIM_P(preslo),AMREX_ARLIM_P(preshi),
     Real* div,AMREX_ARLIM_P(divlo),AMREX_ARLIM_P(divhi),
     Real* divdat,AMREX_ARLIM_P(divdatlo),AMREX_ARLIM_P(divdathi),
     Real* den,
     AMREX_ARLIM_P(denlo),AMREX_ARLIM_P(denhi),
     Real* mom_den,
     AMREX_ARLIM_P(mom_denlo),AMREX_ARLIM_P(mom_denhi),
     Real* elastic,AMREX_ARLIM_P(elasticlo),AMREX_ARLIM_P(elastichi),
     Real* lsdist,AMREX_ARLIM_P(lsdistlo),AMREX_ARLIM_P(lsdisthi),
     Real* visc,AMREX_ARLIM_P(visclo),AMREX_ARLIM_P(vischi),
     Real* conduct,AMREX_ARLIM_P(conductlo),AMREX_ARLIM_P(conducthi),
     Real* trace,AMREX_ARLIM_P(tracelo),AMREX_ARLIM_P(tracehi),
     Real* elasticforce,
     AMREX_ARLIM_P(elasticforcelo),
     AMREX_ARLIM_P(elasticforcehi),
     const Real* problo,
     const Real* probhi,
     const Real* dx,
     const int* lo,const int* hi,
     const int* level,
     const int* finest_level,
     const int* gridno,
     int* visual_tessellate_vfrac,
     int* rz_flag,
     int* nmat,
     int* nparts,
     int* nparts_def,
     int* im_solid_map,
     int* elastic_ncomp,
     Real* slice_data,int* nslice,
     int* nstate_slice,int* slice_dir,
     Real* xslice,const Real* dxfinest,
     int* do_plot,int* do_slice);
 
    void fort_memstatus(int* procnum);

    void fort_combinezones(
     int* total_number_grids,
     int* grids_per_level_array,
     int* levels_array,
     int* bfact_array,
     int* gridno_array,
     int* gridlo_array,
     int* gridhi_array,
     int* finest_level,
     int* nsteps,
     int* num_levels,
     Real* time,
     int* visual_revolve,
     int* plotint,
     int* nmat,
     int* nparts, 
     int* nparts_def,
     int* im_solid_map);


    void FORT_OUTPUTSLICE(Real* time,int* nsteps,int* sliceint,
     Real* slice_data,int* nslice,int* nstate_slice);

    void FORT_IO_COMPARE(
     int* nmat,
     int* nsteps,
     int* do_input,
     int* visual_compare,
     Real* time,
     Real* fabinput,
     AMREX_ARLIM_P(fabinput_lo),AMREX_ARLIM_P(fabinput_hi),
     Real* faboutput,
     AMREX_ARLIM_P(faboutput_lo),AMREX_ARLIM_P(faboutput_hi),
     const int* vislo,const int* vishi,
     int* visual_ncomp);
 
    void FORT_ZALESAK_CELL(
       const Real* xlo,const Real* dx,
       Real* u,const int* domlo,const int* domhi,
       const int* tilelo,const int* tilehi,
       const int* fablo,const int* fabhi,
       int* bfact,
       int* level,
      AMREX_ARLIM_P(ulo),AMREX_ARLIM_P(uhi),const Real* time);

    void fort_init_potential(
      int* nmat,
      Real* presden,AMREX_ARLIM_P(presdenlo),AMREX_ARLIM_P(presdenhi),
      Real* state,AMREX_ARLIM_P(statelo),AMREX_ARLIM_P(statehi),
      const int* tilelo,const int* tilehi,
      const int* fablo,const int* fabhi,
      int* bfact,
      int* level,
      int* presbc,int* dombc,
      const int* domlo,const int* domhi,
      const Real* xlo,const Real* dx,Real* dt,
      Real* gravity_normalized,
      int* gravity_dir,
      Real* angular_velocity,
      int* isweep);

    void FORT_METRICS(
         const Real* xlo,const Real* dx,
         Real* areax,AMREX_ARLIM_P(areaxlo),AMREX_ARLIM_P(areaxhi),
         Real* areay,AMREX_ARLIM_P(areaylo),AMREX_ARLIM_P(areayhi),
         Real* areaz,AMREX_ARLIM_P(areazlo),AMREX_ARLIM_P(areazhi),
         Real* vol,AMREX_ARLIM_P(vollo),AMREX_ARLIM_P(volhi),
         const int* tilelo,const int* tilehi,
         const int* fablo,const int* fabhi,
         int* bfact,
         int* level,
         int* ngrow,int* rzflag);

    void FORT_EDGEAVGDOWN(
       int* enable_spectral,
       int* finest_level,
       int* spectral_override,
       const Real* problo,
       const Real* dxf,
       int* level_c,int* level_f,
       int* bfact_c,int* bfact_f,
       const Real* xlo_fine,const Real* dx,
       const int* grid_type,
       Real* crse,AMREX_ARLIM_P(clo),AMREX_ARLIM_P(chi),
       const Real* fine,AMREX_ARLIM_P(flo),AMREX_ARLIM_P(fhi),
       const Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
       const int* loc, const int* hic, 
       const int* lof,const int* hif,
       int* ncomp);
    

    void FORT_AVGDOWN(
      int* enable_spectral,
      int* finest_level,
      int* spectral_override,
      const Real* problo,
      const Real* dxf,
      int* level_c,int* level_f,
      int* bfact_c,int* bfact_f,
      const Real* xlo_fine,const Real* dx,
      int* ncomp,
      const Real* crse,AMREX_ARLIM_P(clo),AMREX_ARLIM_P(chi),
      const Real* fine,AMREX_ARLIM_P(fine_lo),AMREX_ARLIM_P(fine_hi),
      const Real* mask,AMREX_ARLIM_P(mask_lo),AMREX_ARLIM_P(mask_hi),
      const int* lo, const int* hi,
      const int* lof, const int* hif);


    void fort_avgdown_copy(
      int* enable_spectral,
      int* finest_level,
      int* operation_flag,
      int* dir,
      const Real* problo,
      const Real* dxf,
      int* level_c,int* level_f,
      int* bfact_c,int* bfact_f,
      const Real* xlo_fine,const Real* dx,
      int* ncomp_flux,
      int* ncomp_den,
      int* ncomp_vel,
      const Real* crse,AMREX_ARLIM_P(clo),AMREX_ARLIM_P(chi),
      const Real* fine,AMREX_ARLIM_P(fine_lo),AMREX_ARLIM_P(fine_hi),
      const Real* den_fine,
      AMREX_ARLIM_P(den_fine_lo),AMREX_ARLIM_P(den_fine_hi),
      const Real* vel_fine,
      AMREX_ARLIM_P(vel_fine_lo),AMREX_ARLIM_P(vel_fine_hi),
      const Real* mask,AMREX_ARLIM_P(mask_lo),AMREX_ARLIM_P(mask_hi),
      const Real* fine_LS,AMREX_ARLIM_P(fine_LS_lo),AMREX_ARLIM_P(fine_LS_hi),
      const int* lo, const int* hi,
      const int* lof, const int* hif);


    void fort_interp_copy(
      int* enable_spectral,
      const Real* dxc,
      const Real* dx,
      int* finest_level,
      int* operation_flag,
      const int* tilelo,const int* tilehi,
      const int* fablo,const int* fabhi,
      int* dir,
      const Real* problo,
      int* level_c,int* level,
      int* bfact_c,int* bfact,
      const Real* xlo,
      int* ncomp_flux,
      int* ncomp_den,
      int* ncomp_vel,
      const Real* fine,AMREX_ARLIM_P(fine_lo),AMREX_ARLIM_P(fine_hi),
      const Real* den_crse,
      AMREX_ARLIM_P(den_crse_lo),AMREX_ARLIM_P(den_crse_hi),
      const Real* vel_crse,
      AMREX_ARLIM_P(vel_crse_lo),AMREX_ARLIM_P(vel_crse_hi),
      const Real* masknbr,AMREX_ARLIM_P(masknbr_lo),AMREX_ARLIM_P(masknbr_hi),
      const Real* masksem,AMREX_ARLIM_P(masksem_lo),AMREX_ARLIM_P(masksem_hi),
      const Real* cmasksem,
      AMREX_ARLIM_P(cmasksem_lo),AMREX_ARLIM_P(cmasksem_hi),
      const Real* cLS,AMREX_ARLIM_P(cLS_lo),AMREX_ARLIM_P(cLS_hi),
      int* velbc,
      const int* loc, const int* hic);

    void fort_interp_flux(
      int* enable_spectral,
      const Real* dxc,
      const Real* dx,
      int* finest_level,
      const int* tilelo,const int* tilehi,
      const int* fablo,const int* fabhi,
      int* dir,
      const Real* problo,
      int* level_c,int* level,
      int* bfact_c,int* bfact,
      const Real* xlo,
      int* ncomp_flux,
      const Real* fine,AMREX_ARLIM_P(fine_lo),AMREX_ARLIM_P(fine_hi),
      const Real* crse,AMREX_ARLIM_P(crse_lo),AMREX_ARLIM_P(crse_hi),
      const Real* masknbr,AMREX_ARLIM_P(masknbr_lo),AMREX_ARLIM_P(masknbr_hi),
      const Real* masksem,AMREX_ARLIM_P(masksem_lo),AMREX_ARLIM_P(masksem_hi),
      const Real* cmasksem,
      AMREX_ARLIM_P(cmasksem_lo),AMREX_ARLIM_P(cmasksem_hi),
      int* velbc,
      const int* loc, const int* hic);


    void fort_fillbdry_flux(
      int* sync_iter,
      int* level,
      int* finest_level,
      const int* tilelo,const int* tilehi,
      const int* fablo,const int* fabhi,
      int* dir,
      int* ncomp_flux,
      const Real* fluxtarg,
      AMREX_ARLIM_P(fluxtarg_lo),AMREX_ARLIM_P(fluxtarg_hi),
      const Real* fluxhold,
      AMREX_ARLIM_P(fluxhold_lo),AMREX_ARLIM_P(fluxhold_hi),
      const Real* maskcov,AMREX_ARLIM_P(maskcov_lo),AMREX_ARLIM_P(maskcov_hi),
      const Real* masknbr,AMREX_ARLIM_P(masknbr_lo),AMREX_ARLIM_P(masknbr_hi),
      int* presbc);

    void FORT_AVGDOWN_LOW(
      const Real* problo,
      const Real* dxf,
      int* level_c,int* level_f,
      int* bfact_c,int* bfact_f,
      const Real* xlo_fine,const Real* dx,
      int* ncomp,
      const Real* crse,AMREX_ARLIM_P(clo),AMREX_ARLIM_P(chi),
      const Real* fine,AMREX_ARLIM_P(fine_lo),AMREX_ARLIM_P(fine_hi),
      const int* lo, const int* hi,
      const int* lof, const int* hif);

    void FORT_AVGDOWN_TAG(
      const Real* problo,
      const Real* dxf,
      int* level_c,int* level_f,
      int* bfact_c,int* bfact_f,
      const Real* xlo_fine,const Real* dx,
      int* ncomp,
      const Real* crse,AMREX_ARLIM_P(clo),AMREX_ARLIM_P(chi),
      const Real* fine,AMREX_ARLIM_P(fine_lo),AMREX_ARLIM_P(fine_hi),
      const int* lo, const int* hi,
      const int* lof, const int* hif);

    void fort_avgdown_burning(
      int* velflag,
      const Real* problo,
      const Real* dxf,
      int* level_c,int* level_f,
      int* bfact_c,int* bfact_f,
      const Real* xlo_fine,const Real* dx,
      int* ncomp,
      int* nmat,
      int* nten,
      const Real* crse,AMREX_ARLIM_P(clo),AMREX_ARLIM_P(chi),
      const Real* fine,AMREX_ARLIM_P(fine_lo),AMREX_ARLIM_P(fine_hi),
      const int* lo, const int* hi,
      const int* lof, const int* hif);


    void fort_avgdown_drag(
      const Real* problo,
      const Real* dxf,
      int* level_c,int* level_f,
      int* bfact_c,int* bfact_f,
      const Real* xlo_fine,const Real* dx,
      int* ncomp,
      int* nmat,
      const Real* crse,AMREX_ARLIM_P(clo),AMREX_ARLIM_P(chi),
      const Real* fine,AMREX_ARLIM_P(fine_lo),AMREX_ARLIM_P(fine_hi),
      const int* lo, const int* hi,
      const int* lof, const int* hif);


    void FORT_AVGDOWN_CURV(
      const Real* problo,
      const Real* dxf,
      int* level_c,int* level_f,
      int* bfact_c,int* bfact_f,
      const Real* xlo_fine,const Real* dx,
      int* ncomp,int* nmat,int* nten,
      const Real* crse,AMREX_ARLIM_P(clo),AMREX_ARLIM_P(chi),
      const Real* fine,AMREX_ARLIM_P(fine_lo),AMREX_ARLIM_P(fine_hi),
      const int* lo, const int* hi,
      const int* lof, const int* hif);

    void FORT_MOFAVGDOWN(
      Real* time,
      const Real* problo,
      const Real* dxc,
      const Real* dxf,
      int* bfact_c,int* bfact_f,
      const Real* crse,AMREX_ARLIM_P(clo),AMREX_ARLIM_P(chi),
      const Real* fine,AMREX_ARLIM_P(fine_lo),AMREX_ARLIM_P(fine_hi),
      const int* lo, const int* hi,
      int* nmat);

    void FORT_ERRORAVGDOWN(
      const Real* problo,
      const Real* dxf,
      int* bfact_c,int* bfact_f,
      const Real* crse,AMREX_ARLIM_P(clo),AMREX_ARLIM_P(chi),
      const Real* fine,AMREX_ARLIM_P(fine_lo),AMREX_ARLIM_P(fine_hi),
      const int* lo, const int* hi);


    void fort_addgravity(
      Real* dt,
      Real* cur_time,
      Real* gravity_normalized,
      int* gravity_dir,
      Real* angular_velocity,
      int* level,
      int* finest_level,
      int* facecut_index,
      int* icefacecut_index,
      int* vofface_index,
      int* ncphys,
      int* nmat,
      int* nstate,
      const int* tilelo,const int* tilehi,
      const int* fablo,const int* fabhi,int* bfact,
      const Real* xlo,const Real* dx,int* dir,
      Real* xface,AMREX_ARLIM_P(xfacelo),AMREX_ARLIM_P(xfacehi),
      Real* recon,AMREX_ARLIM_P(reconlo),AMREX_ARLIM_P(reconhi),
      Real* lsnew,AMREX_ARLIM_P(lsnewlo),AMREX_ARLIM_P(lsnewhi),
      Real* macnew,AMREX_ARLIM_P(macnewlo),AMREX_ARLIM_P(macnewhi),
      Real* facegrav,AMREX_ARLIM_P(facegravlo),AMREX_ARLIM_P(facegravhi) );


    void FORT_FLUIDSOLIDCOR(
      int* level,
      int* finest_level,
      int* velcomp,
      int* nsolve,
      int* facecut_index,
      int* icefacecut_index,
      int* ncphys,
      int* project_option,
      const int* tilelo,const int* tilehi,
      const int* fablo,const int* fabhi,int* bfact,
      int* presbc,
      Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi), 
      Real* xface,AMREX_ARLIM_P(xfacelo),AMREX_ARLIM_P(xfacehi), 
      Real* yface,AMREX_ARLIM_P(yfacelo),AMREX_ARLIM_P(yfacehi), 
      Real* zface,AMREX_ARLIM_P(zfacelo),AMREX_ARLIM_P(zfacehi), 
      Real* xgp,AMREX_ARLIM_P(xgplo),AMREX_ARLIM_P(xgphi), 
      Real* ygp,AMREX_ARLIM_P(ygplo),AMREX_ARLIM_P(ygphi), 
      Real* zgp,AMREX_ARLIM_P(zgplo),AMREX_ARLIM_P(zgphi), 
      Real* xsrc,AMREX_ARLIM_P(xsrclo),AMREX_ARLIM_P(xsrchi), 
      Real* ysrc,AMREX_ARLIM_P(ysrclo),AMREX_ARLIM_P(ysrchi), 
      Real* zsrc,AMREX_ARLIM_P(zsrclo),AMREX_ARLIM_P(zsrchi), 
      Real* xdest,AMREX_ARLIM_P(xdestlo),AMREX_ARLIM_P(xdesthi), 
      Real* ydest,AMREX_ARLIM_P(ydestlo),AMREX_ARLIM_P(ydesthi), 
      Real* zdest,AMREX_ARLIM_P(zdestlo),AMREX_ARLIM_P(zdesthi), 
      const Real* xlo,const Real* dx,Real* cur_time,
      int* nmat);

    void fort_summass(
     int* tid,
     int* level,
     int* finest_level,
     int* ncomp_sum_int_user1,
     int* ncomp_sum_int_user2,
     int* adapt_quad_depth,
     int* slice_dir,Real* xslice,
     Real* problo,Real* probhi,
     const Real* xlo,const Real* dx,
     const Real* cellten,AMREX_ARLIM_P(celltenlo),AMREX_ARLIM_P(celltenhi),
     const Real* ls,AMREX_ARLIM_P(lslo),AMREX_ARLIM_P(lshi),
     const Real* masksem,AMREX_ARLIM_P(masksemlo),AMREX_ARLIM_P(masksemhi),
     const Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
     const Real* drag,AMREX_ARLIM_P(draglo),AMREX_ARLIM_P(draghi),
     const Real* slopes,AMREX_ARLIM_P(slopeslo),AMREX_ARLIM_P(slopeshi),
     const Real* den,AMREX_ARLIM_P(denlo),AMREX_ARLIM_P(denhi),
     const Real* vel,AMREX_ARLIM_P(vello),AMREX_ARLIM_P(velhi),
     const Real* visco,AMREX_ARLIM_P(viscolo),AMREX_ARLIM_P(viscohi),
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,
     int* bfact,
     Real* time,
     Real* local_result,
     Real* result_all_levels,
     int* sumdata_type,
     int* sumdata_sweep,
     int* resultsize,
     int* num_cells,
     Real* local_coflow_Z,
     Real* local_coflow_R_of_Z,
     int* Z_dir,
     int* R_dir,
     int* nmat,
     int* ntensor,
     int* den_Ncomp,
     int* isweep); // =0 or 1

     void fort_regionsum(
       int* tid_current,
       int* isweep, //isweep=0 or 1
       int* constant_density_all_time,
       Real* cur_time,
       Real* dt,
       const Real* dx,
       const Real* xlo,
       int* nmat,
       int* nstate,
       Real* snew,AMREX_ARLIM_P(snewlo),AMREX_ARLIM_P(snewhi),
       Real* umacnew,
       AMREX_ARLIM_P(umacnewlo),AMREX_ARLIM_P(umacnewhi),
       Real* vmacnew,
       AMREX_ARLIM_P(vmacnewlo),AMREX_ARLIM_P(vmacnewhi),
       Real* wmacnew,
       AMREX_ARLIM_P(wmacnewlo),AMREX_ARLIM_P(wmacnewhi),
       Real* mdot,
       AMREX_ARLIM_P(mdotlo),AMREX_ARLIM_P(mdothi),
       Real* DEN,AMREX_ARLIM_P(DENlo),AMREX_ARLIM_P(DENhi),
       Real* VOF,AMREX_ARLIM_P(VOFlo),AMREX_ARLIM_P(VOFhi),
       Real* volumefab,
       AMREX_ARLIM_P(volumefablo),AMREX_ARLIM_P(volumefabhi),
       Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
       const int* tilelo, const int* tilehi,
       const int* fablo, const int* fabhi, 
       int* bfact,
       int* level,
       int* finest_level);

    void fort_get_number_regions(int* cpp_number_regions);
    void fort_get_region_data(
      int* isweep,
      Real* cpp_energy_per_kelvin,
      Real* cpp_mass,
      Real* cpp_energy,
      Real* cpp_volume,
      Real* cpp_volume_raster,
      Real* cpp_mass_after,
      Real* cpp_energy_after,
      Real* cpp_volume_after);

    void fort_put_region_data(
      int* isweep,
      Real* cpp_energy_per_kelvin,
      Real* cpp_mass,
      Real* cpp_energy,
      Real* cpp_volume,
      Real* cpp_volume_raster,
      Real* cpp_mass_after,
      Real* cpp_energy_after,
      Real* cpp_volume_after);

    void fort_reduce_sum_regions(int* isweep);

    void fort_sumdot(Real* mass1,
      const Real* dat1,AMREX_ARLIM_P(dat1_lo),AMREX_ARLIM_P(dat1_hi),
      const Real* dat2,AMREX_ARLIM_P(dat2_lo),AMREX_ARLIM_P(dat2_hi),
      const Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
      const int* tilelo,const int* tilehi,
      const int* fablo,const int* fabhi,
      int* bfact,
      int* debug_dot_product,
      int* levelno,const int* gridno,
      int* nsolve);

    void fort_sumdot_ones_size(
      Real* fab_sum,
      int* fab_flag,
      const Real* ones_fab,
      AMREX_ARLIM_P(ones_fablo),AMREX_ARLIM_P(ones_fabhi),
      const Real* type_fab,
      AMREX_ARLIM_P(type_fablo),AMREX_ARLIM_P(type_fabhi),
      const Real* color_fab,
      AMREX_ARLIM_P(color_fablo),AMREX_ARLIM_P(color_fabhi),
      const Real* alpha_fab,
      AMREX_ARLIM_P(alpha_fablo),AMREX_ARLIM_P(alpha_fabhi),
      const Real* mask_fab,
      AMREX_ARLIM_P(mask_fablo),AMREX_ARLIM_P(mask_fabhi),
      const int* tilelo,const int* tilehi,
      const int* fablo,const int* fabhi,
      int* bfact,
      int* levelno,const int* gridno,
      int* nsolve,
      int* presbc,
      int* type_flag,
      int* color_count, 
      int* project_option);


    void fort_sumdot_ones(
      Real* fab_sum,
      int* fab_flag,
      const Real* data_fab,
      AMREX_ARLIM_P(data_fablo),AMREX_ARLIM_P(data_fabhi),
      const Real* ones_fab,
      AMREX_ARLIM_P(ones_fablo),AMREX_ARLIM_P(ones_fabhi),
      const Real* type_fab,
      AMREX_ARLIM_P(type_fablo),AMREX_ARLIM_P(type_fabhi),
      const Real* color_fab,
      AMREX_ARLIM_P(color_fablo),AMREX_ARLIM_P(color_fabhi),
      const Real* alpha_fab,
      AMREX_ARLIM_P(alpha_fablo),AMREX_ARLIM_P(alpha_fabhi),
      const Real* mask_fab,
      AMREX_ARLIM_P(mask_fablo),AMREX_ARLIM_P(mask_fabhi),
      const int* tilelo,const int* tilehi,
      const int* fablo,const int* fabhi,
      int* bfact,
      int* levelno,const int* gridno,
      int* nsolve,
      int* presbc,
      int* type_flag,
      int* color_count, 
      int* project_option);

    void fort_fabcom_ones(
      Real* beta,
      int* singular_patch_flag,
      const Real* data_fab,
      AMREX_ARLIM_P(data_fablo),AMREX_ARLIM_P(data_fabhi),
      const Real* ones_fab,
      AMREX_ARLIM_P(ones_fablo),AMREX_ARLIM_P(ones_fabhi),
      const Real* type_fab,
      AMREX_ARLIM_P(type_fablo),AMREX_ARLIM_P(type_fabhi),
      const Real* color_fab,
      AMREX_ARLIM_P(color_fablo),AMREX_ARLIM_P(color_fabhi),
      const Real* alpha_fab,
      AMREX_ARLIM_P(alpha_fablo),AMREX_ARLIM_P(alpha_fabhi),
      const Real* mask_fab,
      AMREX_ARLIM_P(mask_fablo),AMREX_ARLIM_P(mask_fabhi),
      const int* tilelo,const int* tilehi,
      const int* fablo,const int* fabhi,
      int* bfact,
      int* levelno,const int* gridno,
      int* nsolve,
      int* presbc,
      int* type_flag,
      int* color_count, 
      int* project_option);


    void fort_fabcom(
     const Real* fabx,AMREX_ARLIM_P(fabx_lo),AMREX_ARLIM_P(fabx_hi),
     const Real* faby,AMREX_ARLIM_P(faby_lo),AMREX_ARLIM_P(faby_hi),
     const Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
     Real* fabz,AMREX_ARLIM_P(fabz_lo),AMREX_ARLIM_P(fabz_hi),
     const Real* beta,
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,int* bfact,
     int* nsolve);

    void FORT_DIAGINV(
      const Real* diagfab, 
      AMREX_ARLIM_P(diagfablo),AMREX_ARLIM_P(diagfabhi),
      const Real* resid,AMREX_ARLIM_P(residlo),AMREX_ARLIM_P(residhi),
      const Real* xnew,AMREX_ARLIM_P(xnewlo),AMREX_ARLIM_P(xnewhi),
      const Real* xold,AMREX_ARLIM_P(xoldlo),AMREX_ARLIM_P(xoldhi),
      const Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
      const int* tilelo,const int* tilehi,
      const int* fablo,const int* fabhi,int* bfact);

    void fort_particle_grid(
       int* tid,
       const Real* xlo,const Real* dx,
       const void* particles,
       int Np, 
       const int* tilelo,const int* tilehi,
       const int* fablo,const int* fabhi,
       int* bfact,
       int* level,
       const int* gridno,
       int* ipart);

    void fort_combine_particles(
     int* grids_per_level,int* finest_level,
     int* nsteps,int* ipart,
     int* NS_ncomp_particles,
     int* arrdim,
     Real* time,int* plotint);

}
}/* namespace amrex */
    
#endif
#endif /*_NAVIERSTOKES_F_H_*/
