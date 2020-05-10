#ifndef _NAVIERSTOKES_F_H_
#define _NAVIERSTOKES_F_H_

/*
** $Id: NAVIERSTOKES_F.H,v 1.10 2002/03/26 20:51:08 lijewski Exp $
*/

#ifdef BL_LANG_FORT
#    define FORT_PRESSURE_INDICATOR pressure_indicator
#    define FORT_EOS_PRESSURE eos_pressure
#    define FORT_ADVECTIVE_PRESSURE advective_pressure
#    define FORT_UPDATE_DIV update_div
#    define FORT_COFLOW coflow
#    define FORT_ZALESAKNODE zalesaknode
#    define FORT_EDGEAVGDOWN edgeavgdown
#    define FORT_VOFFLUX_CORRECT vofflux_correct
#    define FORT_FLUIDSOLIDCOR    fluidsolidcor
#    define FORT_ADDGRAVITY  addgravity
#    define FORT_INITPOTENTIAL initpotential
#    define FORT_MOFAVGDOWN  mofavgdown
#    define FORT_ERRORAVGDOWN erroravgdown
#    define FORT_AVGDOWN     avgdown
#    define FORT_AVGDOWN_COPY avgdown_copy
#    define FORT_INTERP_COPY interp_copy
#    define FORT_INTERP_FLUX interp_flux
#    define FORT_FILLBDRY_FLUX fillbdry_flux
#    define FORT_AVGDOWN_LOW avgdown_low
#    define FORT_AVGDOWN_TAG avgdown_tag
#    define FORT_AVGDOWN_BURNING avgdown_burning
#    define FORT_AVGDOWN_CURV avgdown_curv
#    define FORT_SUMMASS     summass
#    define FORT_SUMDOT      sumdot
#    define FORT_DIAGINV     diaginv
#    define FORT_FABCOM      fabcom
#    define FORT_ISOGRID     isogrid
#    define FORT_ISOGRIDSINGLE     isogridsingle
#    define FORT_COMBINETRIANGLESSINGLE combinetrianglessingle
#    define FORT_CELLGRID    cellgrid
#    define FORT_CELLGRID_SANITY cellgrid_sanity
#    define FORT_COMBINETRIANGLES combinetriangles
#    define FORT_COMBINEZONES combinezones
#    define FORT_COMBINEZONES_SANITY combinezones_sanity
#    define FORT_OUTPUTSLICE outputslice
#    define FORT_IO_COMPARE io_compare
#    define FORT_MEMSTATUS memstatus
#    define FORT_METRICS metrics
#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_PRESSURE_INDICATOR PRESSURE_INDICATOR
#    define FORT_EOS_PRESSURE EOS_PRESSURE
#    define FORT_ADVECTIVE_PRESSURE ADVECTIVE_PRESSURE
#    define FORT_UPDATE_DIV UPDATE_DIV
#    define FORT_COFLOW COFLOW
#    define FORT_ZALESAKNODE ZALESAKNODE
#    define FORT_EDGEAVGDOWN EDGEAVGDOWN
#    define FORT_VOFFLUX_CORRECT VOFFLUX_CORRECT
#    define FORT_FLUIDSOLIDCOR    FLUIDSOLIDCOR
#    define FORT_ADDGRAVITY  ADDGRAVITY
#    define FORT_INITPOTENTIAL INITPOTENTIAL
#    define FORT_MOFAVGDOWN  MOFAVGDOWN
#    define FORT_ERRORAVGDOWN ERRORAVGDOWN
#    define FORT_AVGDOWN     AVGDOWN
#    define FORT_AVGDOWN_COPY AVGDOWN_COPY
#    define FORT_INTERP_COPY INTERP_COPY
#    define FORT_INTERP_FLUX INTERP_FLUX
#    define FORT_FILLBDRY_FLUX FILLBDRY_FLUX
#    define FORT_AVGDOWN_LOW AVGDOWN_LOW
#    define FORT_AVGDOWN_TAG AVGDOWN_TAG
#    define FORT_AVGDOWN_BURNING AVGDOWN_BURNING
#    define FORT_AVGDOWN_CURV AVGDOWN_CURV
#    define FORT_SUMMASS     SUMMASS
#    define FORT_SUMDOT      SUMDOT
#    define FORT_DIAGINV     DIAGINV
#    define FORT_FABCOM      FABCOM
#    define FORT_ISOGRID     ISOGRID
#    define FORT_ISOGRIDSINGLE ISOGRIDSINGLE
#    define FORT_COMBINETRIANGLESSINGLE COMBINETRIANGLESSINGLE
#    define FORT_CELLGRID    CELLGRID
#    define FORT_CELLGRID_SANITY CELLGRID_SANITY
#    define FORT_COMBINETRIANGLES COMBINETRIANGLES
#    define FORT_MEMSTATUS MEMSTATUS
#    define FORT_COMBINEZONES COMBINEZONES
#    define FORT_COMBINEZONES_SANITY COMBINEZONES_SANITY
#    define FORT_OUTPUTSLICE OUTPUTSLICE
#    define FORT_IO_COMPARE IO_COMPARE
#    define FORT_METRICS METRICS
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_PRESSURE_INDICATOR pressure_indicator
#    define FORT_EOS_PRESSURE eos_pressure
#    define FORT_ADVECTIVE_PRESSURE advective_pressure
#    define FORT_UPDATE_DIV update_div
#    define FORT_COFLOW coflow
#    define FORT_ZALESAKNODE zalesaknode
#    define FORT_EDGEAVGDOWN edgeavgdown
#    define FORT_VOFFLUX_CORRECT vofflux_correct
#    define FORT_FLUIDSOLIDCOR    fluidsolidcor
#    define FORT_ADDGRAVITY  addgravity
#    define FORT_INITPOTENTIAL initpotential
#    define FORT_MOFAVGDOWN  mofavgdown
#    define FORT_ERRORAVGDOWN erroravgdown
#    define FORT_AVGDOWN     avgdown
#    define FORT_AVGDOWN_COPY avgdown_copy
#    define FORT_INTERP_COPY interp_copy
#    define FORT_INTERP_FLUX interp_flux
#    define FORT_FILLBDRY_FLUX fillbdry_flux
#    define FORT_AVGDOWN_LOW avgdown_low
#    define FORT_AVGDOWN_TAG avgdown_tag
#    define FORT_AVGDOWN_BURNING avgdown_burning
#    define FORT_AVGDOWN_CURV avgdown_curv
#    define FORT_SUMMASS     summass
#    define FORT_SUMDOT      sumdot
#    define FORT_DIAGINV     diaginv
#    define FORT_FABCOM      fabcom
#    define FORT_ISOGRID     isogrid
#    define FORT_ISOGRIDSINGLE     isogridsingle
#    define FORT_COMBINETRIANGLESSINGLE combinetrianglessingle
#    define FORT_CELLGRID    cellgrid
#    define FORT_CELLGRID_SANITY cellgrid_sanity
#    define FORT_COMBINETRIANGLES combinetriangles
#    define FORT_MEMSTATUS memstatus
#    define FORT_COMBINEZONES combinezones
#    define FORT_COMBINEZONES_SANITY combinezones_sanity
#    define FORT_OUTPUTSLICE outputslice
#    define FORT_IO_COMPARE io_compare
#    define FORT_METRICS metrics
#  else
#    define FORT_PRESSURE_INDICATOR pressure_indicator_
#    define FORT_EOS_PRESSURE eos_pressure_
#    define FORT_ADVECTIVE_PRESSURE advective_pressure_
#    define FORT_UPDATE_DIV update_div_
#    define FORT_COFLOW coflow_
#    define FORT_ZALESAKNODE zalesaknode_
#    define FORT_EDGEAVGDOWN edgeavgdown_
#    define FORT_VOFFLUX_CORRECT vofflux_correct_
#    define FORT_FLUIDSOLIDCOR    fluidsolidcor_
#    define FORT_INITPOTENTIAL initpotential_
#    define FORT_ADDGRAVITY  addgravity_
#    define FORT_MOFAVGDOWN  mofavgdown_
#    define FORT_ERRORAVGDOWN erroravgdown_
#    define FORT_AVGDOWN     avgdown_
#    define FORT_AVGDOWN_COPY avgdown_copy_
#    define FORT_INTERP_COPY interp_copy_
#    define FORT_INTERP_FLUX interp_flux_
#    define FORT_FILLBDRY_FLUX fillbdry_flux_
#    define FORT_AVGDOWN_LOW avgdown_low_
#    define FORT_AVGDOWN_TAG avgdown_tag_
#    define FORT_AVGDOWN_BURNING avgdown_burning_
#    define FORT_AVGDOWN_CURV avgdown_curv_
#    define FORT_SUMMASS     summass_
#    define FORT_SUMDOT      sumdot_
#    define FORT_DIAGINV     diaginv_
#    define FORT_FABCOM      fabcom_
#    define FORT_ISOGRID     isogrid_
#    define FORT_ISOGRIDSINGLE     isogridsingle_
#    define FORT_COMBINETRIANGLESSINGLE combinetrianglessingle_
#    define FORT_CELLGRID    cellgrid_
#    define FORT_CELLGRID_SANITY cellgrid_sanity_
#    define FORT_COMBINETRIANGLES combinetriangles_
#    define FORT_MEMSTATUS memstatus_
#    define FORT_COMBINEZONES combinezones_
#    define FORT_COMBINEZONES_SANITY combinezones_sanity_
#    define FORT_OUTPUTSLICE outputslice_
#    define FORT_IO_COMPARE io_compare_
#    define FORT_METRICS metrics_
#  endif
#  endif

#include <AMReX_ArrayLim.H>

namespace amrex{

extern "C" {

    void FORT_ADVECTIVE_PRESSURE(
     int* level,
     int* finest_level,
     const Real* xlo,const Real* dx,
     Real* dt,
     int* make_interface_incomp,
     Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi),
     Real* lsnew,AMREX_ARLIM_P(lsnewlo),AMREX_ARLIM_P(lsnewhi),
     Real* csound,AMREX_ARLIM_P(csoundlo),AMREX_ARLIM_P(csoundhi),
     Real* cvof,AMREX_ARLIM_P(cvoflo),AMREX_ARLIM_P(cvofhi),
     Real* den,AMREX_ARLIM_P(denlo),AMREX_ARLIM_P(denhi),
     Real* mdot,AMREX_ARLIM_P(mdotlo),AMREX_ARLIM_P(mdothi),
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,
     int* bfact,
     int* nmat,int* nden,
     Real* pgrad_dt_factor,
     int* pressure_select_criterion,
     int* project_option);

    void FORT_UPDATE_DIV(
     const Real* xlo,const Real* dx,
     Real* dt,
     Real* csound,AMREX_ARLIM_P(csoundlo),AMREX_ARLIM_P(csoundhi),
     Real* mdot,AMREX_ARLIM_P(mdotlo),AMREX_ARLIM_P(mdothi),
     Real* pnew,AMREX_ARLIM_P(pnewlo),AMREX_ARLIM_P(pnewhi),
     Real* divnew,AMREX_ARLIM_P(divnewlo),AMREX_ARLIM_P(divnewhi),
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,
     int* bfact,
     int* nmat);


    void FORT_EOS_PRESSURE(
     int* level,
     int* finest_level,
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

    void FORT_COFLOW(
     Real* time,int* js,int* je,int* NN,Real* ZZ,Real* FF,
     int* dirx,int* diry,int* cut_flag);

    void FORT_CELLGRID(
     int* tid,
     int* bfact,
     Real* fabout,AMREX_ARLIM_P(faboutlo),AMREX_ARLIM_P(fabouthi),
     const int* vislo,const int* vishi,
     int* vis_ncomp,
     Real* maskSEM,AMREX_ARLIM_P(maskSEMlo),AMREX_ARLIM_P(maskSEMhi),
     Real* vel,AMREX_ARLIM_P(vello),AMREX_ARLIM_P(velhi),
     Real* velsolid,AMREX_ARLIM_P(velsolidlo),AMREX_ARLIM_P(velsolidhi),
     Real* vof,AMREX_ARLIM_P(voflo),AMREX_ARLIM_P(vofhi),
     Real* pres,AMREX_ARLIM_P(preslo),AMREX_ARLIM_P(preshi),
     Real* div,AMREX_ARLIM_P(divlo),AMREX_ARLIM_P(divhi),
     Real* divdat,AMREX_ARLIM_P(divdatlo),AMREX_ARLIM_P(divdathi),
     Real* den,AMREX_ARLIM_P(denlo),AMREX_ARLIM_P(denhi),
     Real* elastic,AMREX_ARLIM_P(elasticlo),AMREX_ARLIM_P(elastichi),
     Real* lsdist,AMREX_ARLIM_P(lsdistlo),AMREX_ARLIM_P(lsdisthi),
     Real* visc,AMREX_ARLIM_P(visclo),AMREX_ARLIM_P(vischi),
     Real* trace,AMREX_ARLIM_P(tracelo),AMREX_ARLIM_P(tracehi),
     const Real* problo,
     const Real* probhi,
     const Real* dx,
     const int* lo,const int* hi,
     const int* level,
     const int* finest_level,
     const int* gridno,
     int* visual_tessellate_vfrac,
     int* visual_option,
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


    void FORT_CELLGRID_SANITY(
     int* tid,
     int* data_dir,
     int* bfact,
     int* ncomp,
     Real* datafab,AMREX_ARLIM_P(datafablo),AMREX_ARLIM_P(datafabhi),
     const Real* problo,
     const Real* probhi,
     const Real* dx,
     const int* lo,const int* hi,
     const int* level,
     const int* finest_level,
     const int* gridno,
     int* rz_flag);

    void FORT_ISOGRID(
     int* tid,
     int* visual_tessellate_vfrac,
     Real* recon,AMREX_ARLIM_P(reconlo),AMREX_ARLIM_P(reconhi),
     const Real* xlo,const Real* dx,
     Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,int* bfact,
     const int* level,const int* gridno,int* nmat);

    void FORT_COMBINETRIANGLES(
     int* grids_per_level,int* finest_level,
     int* nsteps,int* imaterial,int* arrdim,
     Real* time,int* plotint);


    void FORT_ISOGRIDSINGLE(
     Real* ls,AMREX_ARLIM_P(lslo),AMREX_ARLIM_P(lshi),
     const Real* xlo,const Real* dx,
     Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
     const int* lo,const int* hi,
     const int* level,const int* gridno);

    void FORT_COMBINETRIANGLESSINGLE(
     int* grids_per_level,int* finest_level,
     int* nsteps,int* arrdim);

 
    void FORT_MEMSTATUS(int* procnum);

    void FORT_COMBINEZONES(
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
     int* visual_option,
     int* visual_revolve,
     int* plotint,
     int* nmat,
     int* nparts, 
     int* nparts_def,
     int* im_solid_map);


    void FORT_COMBINEZONES_SANITY(
     char* root_char_array,
     int* n_root,
     int* data_dir,
     int* total_number_grids,
     int* grids_per_level_array,
     int* levels_array,
     int* bfact_array,
     int* gridno_array,
     int* gridlo_array,
     int* gridhi_array,
     int* finest_level,
     int* SDC_outer_sweeps,
     int* slab_step,
     int* data_id,
     int* nsteps,
     int* num_levels,
     Real* time,
     int* visual_option,
     int* visual_revolve,
     int* ncomp);


    void FORT_OUTPUTSLICE(Real* time,int* nsteps,int* sliceint,
     Real* slice_data,int* nslice,int* nstate_slice,
     int* visual_option);

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
     int* vis_ncomp);
 
    void FORT_ZALESAKNODE(
       const Real* xlo,const Real* dx,
       Real* u,const int* domlo,const int* domhi,
       const int* tilelo,const int* tilehi,
       const int* fablo,const int* fabhi,
       int* bfact,
       int* level,
      AMREX_ARLIM_P(ulo),AMREX_ARLIM_P(uhi),const Real* time);

    void FORT_INITPOTENTIAL(
      int* nmat,
      int* ngrow,
      int* override_density,
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
       const int* dir,
       Real* crse,AMREX_ARLIM_P(clo),AMREX_ARLIM_P(chi),
       const Real* fine,AMREX_ARLIM_P(flo),AMREX_ARLIM_P(fhi),
       const Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
       const int* loc, const int* hic, 
       const int* lof,const int* hif,
       int* ncomp);
    

    void FORT_VOFFLUX_CORRECT(
       int* finest_level,
       const Real* problo,
       const Real* dxf,
       int* level_c,int* level_f,
       int* bfact_c,int* bfact_f,
       const Real* xlo_fine,const Real* dx,
       const int* dir,
       Real* crse,AMREX_ARLIM_P(clo),AMREX_ARLIM_P(chi),
       const Real* fine,AMREX_ARLIM_P(flo),AMREX_ARLIM_P(fhi),
       const int* loc, const int* hic, 
       const int* lof,const int* hif,
       int* nmat);

 
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


    void FORT_AVGDOWN_COPY(
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
      const Real* den_fine,AMREX_ARLIM_P(den_fine_lo),AMREX_ARLIM_P(den_fine_hi),
      const Real* vel_fine,AMREX_ARLIM_P(vel_fine_lo),AMREX_ARLIM_P(vel_fine_hi),
      const Real* mask,AMREX_ARLIM_P(mask_lo),AMREX_ARLIM_P(mask_hi),
      const Real* fine_LS,AMREX_ARLIM_P(fine_LS_lo),AMREX_ARLIM_P(fine_LS_hi),
      const int* lo, const int* hi,
      const int* lof, const int* hif);


    void FORT_INTERP_COPY(
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
      const Real* den_crse,AMREX_ARLIM_P(den_crse_lo),AMREX_ARLIM_P(den_crse_hi),
      const Real* vel_crse,AMREX_ARLIM_P(vel_crse_lo),AMREX_ARLIM_P(vel_crse_hi),
      const Real* masknbr,AMREX_ARLIM_P(masknbr_lo),AMREX_ARLIM_P(masknbr_hi),
      const Real* masksem,AMREX_ARLIM_P(masksem_lo),AMREX_ARLIM_P(masksem_hi),
      const Real* cmasksem,AMREX_ARLIM_P(cmasksem_lo),AMREX_ARLIM_P(cmasksem_hi),
      const Real* cLS,AMREX_ARLIM_P(cLS_lo),AMREX_ARLIM_P(cLS_hi),
      int* velbc,
      const int* loc, const int* hic);

    void FORT_INTERP_FLUX(
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
      const Real* cmasksem,AMREX_ARLIM_P(cmasksem_lo),AMREX_ARLIM_P(cmasksem_hi),
      int* velbc,
      const int* loc, const int* hic);


    void FORT_FILLBDRY_FLUX(
      int* sync_iter,
      int* level,
      int* finest_level,
      const int* tilelo,const int* tilehi,
      const int* fablo,const int* fabhi,
      int* dir,
      int* ncomp_flux,
      const Real* fluxtarg,AMREX_ARLIM_P(fluxtarg_lo),AMREX_ARLIM_P(fluxtarg_hi),
      const Real* fluxhold,AMREX_ARLIM_P(fluxhold_lo),AMREX_ARLIM_P(fluxhold_hi),
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

    void FORT_AVGDOWN_BURNING(
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


    void FORT_ADDGRAVITY(
      Real* denconst_gravity,
      int* nsolveMM_FACE,
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
      Real* snew,AMREX_ARLIM_P(snewlo),AMREX_ARLIM_P(snewhi),
      Real* macnew,AMREX_ARLIM_P(macnewlo),AMREX_ARLIM_P(macnewhi),
      Real* cellgrav,AMREX_ARLIM_P(cellgravlo),AMREX_ARLIM_P(cellgravhi),
      Real* facegrav,AMREX_ARLIM_P(facegravlo),AMREX_ARLIM_P(facegravhi) );


    void FORT_FLUIDSOLIDCOR(
      int* level,
      int* finest_level,
      int* velcomp,
      int* nsolve,
      int* nsolveMM,
      int* nsolveMM_face,
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
      Real* alt_xface,AMREX_ARLIM_P(alt_xfacelo),AMREX_ARLIM_P(alt_xfacehi), 
      Real* alt_yface,AMREX_ARLIM_P(alt_yfacelo),AMREX_ARLIM_P(alt_yfacehi), 
      Real* alt_zface,AMREX_ARLIM_P(alt_zfacelo),AMREX_ARLIM_P(alt_zfacehi), 
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

    void FORT_SUMMASS (
     int* tid,
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
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,
     int* bfact,
     Real* time,
     Real* local_result,
     Real* result_all_levels,
     int* sumdata_type,
     int* sumdata_sweep,
     int* resultsize,
     int* NN,Real* ZZ,Real* FF,
     int* dirx,int* diry,int* cut_flag,
     int* nmat,
     int* ntensorMM,
     int* den_Ncomp,
     int* isweep);


    void FORT_SUMDOT(Real* mass1,
      const Real* dat1,AMREX_ARLIM_P(dat1_lo),AMREX_ARLIM_P(dat1_hi),
      const Real* dat2,AMREX_ARLIM_P(dat2_lo),AMREX_ARLIM_P(dat2_hi),
      const Real* dotmask,AMREX_ARLIM_P(dotmasklo),AMREX_ARLIM_P(dotmaskhi),
      const Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
      const int* tilelo,const int* tilehi,
      const int* fablo,const int* fabhi,
      int* bfact,
      int* debug_dot_product,
      int* levelno,const int* gridno,
      int* nsolve,
      int* nsolveMM,
      int* num_materials_face);

    void FORT_FABCOM(
     const Real* fabx,AMREX_ARLIM_P(fabx_lo),AMREX_ARLIM_P(fabx_hi),
     const Real* faby,AMREX_ARLIM_P(faby_lo),AMREX_ARLIM_P(faby_hi),
     const Real* dotmask,AMREX_ARLIM_P(dotmasklo),AMREX_ARLIM_P(dotmaskhi),
     const Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
     Real* fabz,AMREX_ARLIM_P(fabz_lo),AMREX_ARLIM_P(fabz_hi),
     const Real* beta,
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,int* bfact,
     int* nsolve,
     int* nsolveMM,
     int* num_materials_face);

    void FORT_DIAGINV(
      const Real* diagfab, 
      AMREX_ARLIM_P(diagfablo),AMREX_ARLIM_P(diagfabhi),
      const Real* resid,AMREX_ARLIM_P(residlo),AMREX_ARLIM_P(residhi),
      const Real* xnew,AMREX_ARLIM_P(xnewlo),AMREX_ARLIM_P(xnewhi),
      const Real* xold,AMREX_ARLIM_P(xoldlo),AMREX_ARLIM_P(xoldhi),
      const Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
      const int* tilelo,const int* tilehi,
      const int* fablo,const int* fabhi,int* bfact);

}
}/* namespace amrex */
    
#endif
#endif /*_NAVIERSTOKES_F_H_*/
