#ifndef _LEVEL_F_H_
#define _LEVEL_F_H_ 

#ifdef BL_LANG_FORT
#    define FORT_CELLFACEINIT cellfaceinit
#    define FORT_NEUMANNEXTRAP neumannextrap
#    define FORT_PURGEFLOTSAM purgeflotsam
#    define FORT_INITRECALESCE initrecalesce
#    define FORT_GETTYPEFAB gettypefab
#    define FORT_LEVELRECOLOR levelrecolor
#    define FORT_COLORFILL colorfill
#    define FORT_GRIDRECOLOR gridrecolor
#    define FORT_LEVELCOLORINIT levelcolorinit
#    define FORT_AVGDOWNCOLOR avgdowncolor
#    define FORT_COPYFINECOARSECOLOR copyfinecoarsecolor

#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_CELLFACEINIT CELLFACEINIT
#    define FORT_NEUMANNEXTRAP NEUMANNEXTRAP
#    define FORT_PURGEFLOTSAM PURGEFLOTSAM
#    define FORT_INITSUPERELASTIC INITSUPERELASTIC
#    define FORT_INITRECALESCE INITRECALESCE
#    define FORT_GETTYPEFAB GETTYPEFAB
#    define FORT_LEVELRECOLOR LEVELRECOLOR
#    define FORT_COLORFILL COLORFILL
#    define FORT_GRIDRECOLOR GRIDRECOLOR
#    define FORT_LEVELCOLORINIT LEVELCOLORINIT
#    define FORT_AVGDOWNCOLOR AVGDOWNCOLOR
#    define FORT_COPYFINECOARSECOLOR COPYFINECOARSECOLOR
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_CELLFACEINIT cellfaceinit
#    define FORT_NEUMANNEXTRAP neumannextrap
#    define FORT_PURGEFLOTSAM purgeflotsam
#    define FORT_INITRECALESCE initrecalesce
#    define FORT_GETTYPEFAB gettypefab
#    define FORT_LEVELRECOLOR levelrecolor
#    define FORT_COLORFILL colorfill
#    define FORT_GRIDRECOLOR gridrecolor
#    define FORT_LEVELCOLORINIT levelcolorinit
#    define FORT_AVGDOWNCOLOR avgdowncolor
#    define FORT_COPYFINECOARSECOLOR copyfinecoarsecolor
#  else
#    define FORT_CELLFACEINIT cellfaceinit_
#    define FORT_NEUMANNEXTRAP neumannextrap_
#    define FORT_PURGEFLOTSAM purgeflotsam_
#    define FORT_INITRECALESCE initrecalesce_
#    define FORT_GETTYPEFAB gettypefab_
#    define FORT_LEVELRECOLOR levelrecolor_
#    define FORT_COLORFILL colorfill_
#    define FORT_GRIDRECOLOR gridrecolor_
#    define FORT_LEVELCOLORINIT levelcolorinit_
#    define FORT_AVGDOWNCOLOR avgdowncolor_
#    define FORT_COPYFINECOARSECOLOR copyfinecoarsecolor_
#  endif
#  endif

#include <AMReX_ArrayLim.H>

namespace amrex{

  extern "C" {

    void FORT_GETTYPEFAB(
     const Real* source_fab,
     AMREX_ARLIM_P(source_fablo),AMREX_ARLIM_P(source_fabhi),
     const Real* typefab,AMREX_ARLIM_P(typefablo),AMREX_ARLIM_P(typefabhi),
     const Real* xlo,
     const Real* dx,
     const int* tilelo, const int* tilehi,
     const int* fablo, const int* fabhi,int* bfact,
     int* type_flag,
     int* nmat,
     int* ncomp_type,
     int* ncomp_source,
     int* zero_diag_flag);

     void fort_getcolorsum(
       int* tid_current,
       int* operation_flag,
       int* sweep_num,
       int* tessellate,
       int* distribute_mdot_evenly,
       int* constant_volume_mdot,
       Real* latent_heat,
       int* distribute_from_target,
       int* constant_density_all_time,
       Real* dt,
       const Real* dx,
       const Real* xlo,
       int* nmat,
       int* nten,
       int* nstate,
       Real* snew,AMREX_ARLIM_P(snewlo),AMREX_ARLIM_P(snewhi),
       Real* mdot,
       AMREX_ARLIM_P(mdotlo),AMREX_ARLIM_P(mdothi),
       Real* mdot_comp,
       AMREX_ARLIM_P(mdot_comp_lo),AMREX_ARLIM_P(mdot_comp_hi),
       Real* LS,AMREX_ARLIM_P(LSlo),AMREX_ARLIM_P(LShi),
       Real* VEL,AMREX_ARLIM_P(VELlo),AMREX_ARLIM_P(VELhi),
       Real* DEN,AMREX_ARLIM_P(DENlo),AMREX_ARLIM_P(DENhi),
       Real* VOF,AMREX_ARLIM_P(VOFlo),AMREX_ARLIM_P(VOFhi),
       Real* facefab,AMREX_ARLIM_P(facefablo),AMREX_ARLIM_P(facefabhi),
       Real* xface,AMREX_ARLIM_P(xfacelo),AMREX_ARLIM_P(xfacehi),
       Real* yface,AMREX_ARLIM_P(yfacelo),AMREX_ARLIM_P(yfacehi),
       Real* zface,AMREX_ARLIM_P(zfacelo),AMREX_ARLIM_P(zfacehi),
       Real* areax,AMREX_ARLIM_P(areaxlo),AMREX_ARLIM_P(areaxhi),
       Real* areay,AMREX_ARLIM_P(areaylo),AMREX_ARLIM_P(areayhi),
       Real* areaz,AMREX_ARLIM_P(areazlo),AMREX_ARLIM_P(areazhi),
       Real* cellfab,AMREX_ARLIM_P(cellfablo),AMREX_ARLIM_P(cellfabhi),
       Real* typefab,AMREX_ARLIM_P(typefablo),AMREX_ARLIM_P(typefabhi),
       Real* color,AMREX_ARLIM_P(colorlo),AMREX_ARLIM_P(colorhi),
       Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
       const int* tilelo, const int* tilehi,
       const int* fablo, const int* fabhi, 
       int* bfact,
       int* level,
       int* finest_level,
       int* rzflag,
       int* num_colors,
       Real* cum_blobarray,
       Real* cum_mdot_array,
       Real* cum_mdot_comp_array,
       Real* level_blobarray,
       int* level_blobtypearray,
       Real* level_mdot_array,
       Real* level_mdot_comp_array,
       Real* level_mdot_redistribute_array,
       Real* level_mdot_comp_redistribute_array,
       int* arraysize,
       int* mdot_arraysize,
       int* ncomp_mdot_alloc,
       int* ncomp_mdot,
       int* levelbc,
       int* velbc,
       int* material_type_lowmach,
       int* material_type_visual,
       int* nface,
       int* nface_dst,
       int* ncellfrac);

     void fort_get_lowmach_divu(
       int* tid_current,
       int* sweep_num,
       int* constant_density_all_time,
       Real* dt,
       const Real* dx,
       const Real* xlo,
       int* nmat,
       int* nten,
       Real* mdot_local, // mdot_local initialized to 0
       AMREX_ARLIM_P(mdot_locallo),AMREX_ARLIM_P(mdot_localhi),
       Real* mdot_global, // mdot_global incremented 
       AMREX_ARLIM_P(mdot_globallo),AMREX_ARLIM_P(mdot_globalhi),
       Real* LS,AMREX_ARLIM_P(LSlo),AMREX_ARLIM_P(LShi),
       Real* DEN,AMREX_ARLIM_P(DENlo),AMREX_ARLIM_P(DENhi),
       Real* DTDt,AMREX_ARLIM_P(DTDtlo),AMREX_ARLIM_P(DTDthi),
       Real* VOF,AMREX_ARLIM_P(VOFlo),AMREX_ARLIM_P(VOFhi),
       Real* typefab,AMREX_ARLIM_P(typefablo),AMREX_ARLIM_P(typefabhi),
       Real* color,AMREX_ARLIM_P(colorlo),AMREX_ARLIM_P(colorhi),
       Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
       const int* tilelo, const int* tilehi,
       const int* fablo, const int* fabhi, 
       int* bfact,
       int* level,
       int* finest_level,
       int* rzflag,
       int* num_colors,
       Real* level_blobarray,
       Real* level_mdot_array,
       int* arraysize,
       int* mdot_arraysize,
       int* material_type_lowmach);


    void FORT_LEVELRECOLOR(
     Real* color,AMREX_ARLIM_P(colorlo),AMREX_ARLIM_P(colorhi),
     const Real* xlo,const Real* dx,
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,int* bfact,
     int* domaincolormap,
     int* max_colors_level,int* level,int* base_level,int* arrsize);

    void FORT_GRIDRECOLOR(
     Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
     Real* color,AMREX_ARLIM_P(colorlo),AMREX_ARLIM_P(colorhi),
     const Real* xlo,const Real* dx,
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,
     int* bfact,
     int* levelcolormap,
     int* max_colors_grid,
     int* number_grids,int* arrsize);

    void FORT_COLORFILL(
     const Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
     const Real* typefab,AMREX_ARLIM_P(typefablo),AMREX_ARLIM_P(typefabhi),
     const Real* color,AMREX_ARLIM_P(colorlo),AMREX_ARLIM_P(colorhi),
     const int* ijk,AMREX_ARLIM_P(ijklo),AMREX_ARLIM_P(ijkhi),
     const int* lo,const int* hi,int* ipass,int* number_grids,
     int* color_per_grid,const int* gridno,int* max_colors_grid,
     int* typedim);

    void FORT_LEVELCOLORINIT(
     Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
     Real* color,AMREX_ARLIM_P(colorlo),AMREX_ARLIM_P(colorhi),
     const Real* xlo,const Real* dx,
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,
     int* bfact,
     int* level_color,
     int* max_colors_level,
     int* arrsize,int* check_corners);

    void FORT_AVGDOWNCOLOR(
     const Real* prob_lo,
     const Real* dxf,
     int* bfact_f,int* bfact,
     const Real* xlo_fine,
     const Real* dx,
     const Real* crse,AMREX_ARLIM_P(crselo),AMREX_ARLIM_P(crsehi),
     const Real* fine,AMREX_ARLIM_P(fine_lo),AMREX_ARLIM_P(fine_hi),
     const Real* typef,AMREX_ARLIM_P(typef_lo),AMREX_ARLIM_P(typef_hi),
     const Real* typec,AMREX_ARLIM_P(typec_lo),AMREX_ARLIM_P(typec_hi),
     const int* clo,const int* chi);

    void FORT_COPYFINECOARSECOLOR(
     const Real* prob_lo,
     const Real* dxf,
     int* bfact_f,int* bfact,
     const Real* xlo_fine,
     const Real* dx,
     const Real* crse,AMREX_ARLIM_P(crselo),AMREX_ARLIM_P(crsehi),
     const Real* fine,AMREX_ARLIM_P(fine_lo),AMREX_ARLIM_P(fine_hi),
     const Real* typef,AMREX_ARLIM_P(typef_lo),AMREX_ARLIM_P(typef_hi),
     const Real* maskf,AMREX_ARLIM_P(maskf_lo),AMREX_ARLIM_P(maskf_hi),
     const int* clo,const int* chi,const int* flo,const int* fhi,
     int* zero_diag_flag);
    

    void FORT_NEUMANNEXTRAP(
     const int* lo,const int* hi,
     const Real* ls,AMREX_ARLIM_P(lslo),AMREX_ARLIM_P(lshi),
     int* bc,int* dombc,
     const int* domlo,const int* domhi);


     void fort_renormalize_prescribe(
       int* tid,
       int* level,int* finest_level,
       Real* solid_time,
       const int* tilelo,const int* tilehi,
       const int* fablo,const int* fabhi,int* bfact,
       Real* vofnew,AMREX_ARLIM_P(vofnewlo),AMREX_ARLIM_P(vofnewhi),
       Real* solxfab,AMREX_ARLIM_P(solxfablo),AMREX_ARLIM_P(solxfabhi),
       Real* solyfab,AMREX_ARLIM_P(solyfablo),AMREX_ARLIM_P(solyfabhi),
       Real* solzfab,AMREX_ARLIM_P(solzfablo),AMREX_ARLIM_P(solzfabhi),
       Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi),
       Real* LS,AMREX_ARLIM_P(LSlo),AMREX_ARLIM_P(LShi),
       Real* mof,AMREX_ARLIM_P(moflo),AMREX_ARLIM_P(mofhi),
       Real* den,AMREX_ARLIM_P(denlo),AMREX_ARLIM_P(denhi),
       Real* vel,AMREX_ARLIM_P(vello),AMREX_ARLIM_P(velhi),
       Real* velnew,AMREX_ARLIM_P(velnewlo),AMREX_ARLIM_P(velnewhi),
       Real* dennew,AMREX_ARLIM_P(dennewlo),AMREX_ARLIM_P(dennewhi),
       Real* lsnew,AMREX_ARLIM_P(lsnewlo),AMREX_ARLIM_P(lsnewhi),
       const Real* xlo,const Real* dx, 
       Real* time,
       int* nmat,
       int* nten,
       int* nparts,
       int* nparts_def,
       int* im_solid_map,
       int* renormalize_only,
       int* solidheat_flag,
       int* num_LS_extrap,
       int* num_LS_extrap_iter,
       int* LS_extrap_iter,
       int* ngrow_distance,
       int* constant_density_all_time);


     void FORT_PURGEFLOTSAM(
       Real* delta_mass,
       int* truncate_volume_fractions,
       Real* truncate_thickness,
       int* level,int* finest_level,
       Real* time,
       const int* tilelo,const int* tilehi,
       const int* fablo,const int* fabhi,int* bfact,
       Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi),
       Real* vofnew,AMREX_ARLIM_P(vofnewlo),AMREX_ARLIM_P(vofnewhi),
       Real* LS,AMREX_ARLIM_P(LSlo),AMREX_ARLIM_P(LShi),
       const Real* xlo,const Real* dx, 
       int* nmat);

     void FORT_INITRECALESCE(
       int* recalesce_material,
       Real* recalesce_state_old,
       int* recalesce_num_state,int* nmat);

      void fort_cell_to_mac(
        int* ncomp_mgoni,
        int* ncomp_xp,
        int* ncomp_xgp,
        int* simple_AMR_BC_flag,
        int* nsolve,
        int* tileloop,
        int* dir,
        int* operation_flag,
        int* energyflag,
        Real* beta,Real* visc_coef,
        int* enable_spectral,
        int* fluxvel_index,
        int* fluxden_index,
        int* facevel_index,
        int* facecut_index,
        int* icefacecut_index,
        int* curv_index,
        int* ignore_div_up, 
        int* pforce_index,
        int* faceden_index,
        int* icemask_index,
        int* massface_index,
        int* vofface_index,
        int* ncphys,
        int* constant_density_all_time,
        int* presbc,
        int* velbc,
        int* slab_step,
        Real* dt,
        Real* time,
        const Real* xlo,const Real* dx,
        int* spectral_loop,
        int* ncfluxreg,
        Real* semflux,AMREX_ARLIM_P(semfluxlo),AMREX_ARLIM_P(semfluxhi),
        Real* maskfine,AMREX_ARLIM_P(maskfinelo),AMREX_ARLIM_P(maskfinehi),
        Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi),
        Real* maskSEM,AMREX_ARLIM_P(maskSEMlo),AMREX_ARLIM_P(maskSEMhi),
        Real* levelPC,AMREX_ARLIM_P(levelPClo),AMREX_ARLIM_P(levelPChi),
        Real* sol,AMREX_ARLIM_P(sollo),AMREX_ARLIM_P(solhi),
        Real* xcut,AMREX_ARLIM_P(xcutlo),AMREX_ARLIM_P(xcuthi), 
        Real* xface,AMREX_ARLIM_P(xfacelo),AMREX_ARLIM_P(xfacehi), 
        Real* recon,AMREX_ARLIM_P(reconlo),AMREX_ARLIM_P(reconhi), 
        Real* xgp,AMREX_ARLIM_P(xgplo),AMREX_ARLIM_P(xgphi), 
        Real* xp,AMREX_ARLIM_P(xplo),AMREX_ARLIM_P(xphi), 
        Real* xvel,AMREX_ARLIM_P(xvello),AMREX_ARLIM_P(xvelhi), 
	Real* vel,AMREX_ARLIM_P(vello),AMREX_ARLIM_P(velhi),
	Real* pres,AMREX_ARLIM_P(preslo),AMREX_ARLIM_P(preshi),
	Real* den,AMREX_ARLIM_P(denlo),AMREX_ARLIM_P(denhi),
	Real* mgoni,AMREX_ARLIM_P(mgonilo),AMREX_ARLIM_P(mgonihi),
	Real* colorfab,AMREX_ARLIM_P(colorfablo),AMREX_ARLIM_P(colorfabhi),
	Real* typefab,AMREX_ARLIM_P(typefablo),AMREX_ARLIM_P(typefabhi),
        const int* tilelo,const int* tilehi,
        const int* fablo,const int* fabhi,
        int* bfact,int* bfact_c,int* bfact_f,
        int* level,int* finest_level,
        const int* rz_flag,
        const int* domlo,const int* domhi,
        int* nmat,
        int* nparts,
        int* nparts_def,
        int* im_solid_map,
        Real* added_weight,
        Real* blob_array,
        int* blob_array_size,
        int* num_colors,
        int* nten,
        int* project_option,
        int* SEM_upwind);


      void fort_inc_temp(
        Real* beta,
        int* nmat,
        int* level,
        int* finest_level,
        int* ncomp_state,
        const int* tilelo,const int* tilehi,
        const int* fablo,const int* fabhi,
        Real* state,AMREX_ARLIM_P(statelo),AMREX_ARLIM_P(statehi), 
         //1==not covered
        Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi)); 


      void fort_mac_to_cell(
        int* ns_time_order,
        int* divu_outer_sweeps,
        int* num_divu_outer_sweeps,
        int* operation_flag,
        int* energyflag,
        int* constant_density_all_time,
        int* nmat,
        int* nparts,
        int* nparts_def,
        int* im_solid_map,
        Real* added_weight,
        int* nten,
        int* level,
        int* finest_level,
        int* project_option,
        int* enable_spectral,
        int* fluxvel_index,
        int* fluxden_index,
        int* facevel_index,
        int* facecut_index,
        int* icefacecut_index,
        int* curv_index,
        int* ignore_div_up, 
        int* pforce_index,
        int* faceden_index,
        int* icemask_index,
        int* massface_index,
        int* vofface_index,
        int* ncphys,
        int* velbc,
        int* presbc,
        Real* time,
        int* slab_step,
        Real* dt,
        const Real* xlo,const Real* dx,
        const int* tilelo,const int* tilehi,
        const int* fablo,const int* fabhi,
        int* bfact,
        Real* xp,AMREX_ARLIM_P(xplo),AMREX_ARLIM_P(xphi), 
        Real* yp,AMREX_ARLIM_P(yplo),AMREX_ARLIM_P(yphi), 
        Real* zp,AMREX_ARLIM_P(zplo),AMREX_ARLIM_P(zphi), 
        Real* xvel,AMREX_ARLIM_P(xvello),AMREX_ARLIM_P(xvelhi), 
        Real* yvel,AMREX_ARLIM_P(yvello),AMREX_ARLIM_P(yvelhi), 
        Real* zvel,AMREX_ARLIM_P(zvello),AMREX_ARLIM_P(zvelhi), 
        Real* xface,AMREX_ARLIM_P(xfacelo),AMREX_ARLIM_P(xfacehi), 
        Real* yface,AMREX_ARLIM_P(yfacelo),AMREX_ARLIM_P(yfacehi), 
        Real* zface,AMREX_ARLIM_P(zfacelo),AMREX_ARLIM_P(zfacehi), 
        Real* ax,AMREX_ARLIM_P(axlo),AMREX_ARLIM_P(axhi), 
        Real* ay,AMREX_ARLIM_P(aylo),AMREX_ARLIM_P(ayhi), 
        Real* az,AMREX_ARLIM_P(azlo),AMREX_ARLIM_P(azhi), 
        Real* vol,AMREX_ARLIM_P(vollo),AMREX_ARLIM_P(volhi), 
        Real* rhs,AMREX_ARLIM_P(rhslo),AMREX_ARLIM_P(rhshi), 
        Real* veldest,AMREX_ARLIM_P(veldestlo),AMREX_ARLIM_P(veldesthi), 
        Real* dendest,AMREX_ARLIM_P(dendestlo),AMREX_ARLIM_P(dendesthi), 
        Real* maskfine,
        AMREX_ARLIM_P(maskfinelo),AMREX_ARLIM_P(maskfinehi),//1=fine/fine
        Real* maskcov,
        AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi),//1=not covered
        Real* maskSEM,AMREX_ARLIM_P(maskSEMlo),AMREX_ARLIM_P(maskSEMhi),
        Real* levelPC,AMREX_ARLIM_P(levelPClo),AMREX_ARLIM_P(levelPChi),
        Real* solxfab,
	AMREX_ARLIM_P(solxfablo),AMREX_ARLIM_P(solxfabhi),
        Real* solyfab,
	AMREX_ARLIM_P(solyfablo),AMREX_ARLIM_P(solyfabhi),
        Real* solzfab,
	AMREX_ARLIM_P(solzfablo),AMREX_ARLIM_P(solzfabhi),
        Real* cterm,AMREX_ARLIM_P(ctermlo),AMREX_ARLIM_P(ctermhi),
        Real* pold,AMREX_ARLIM_P(poldlo),AMREX_ARLIM_P(poldhi),
        Real* denold,AMREX_ARLIM_P(denoldlo),AMREX_ARLIM_P(denoldhi),
        Real* ustar,AMREX_ARLIM_P(ustarlo),AMREX_ARLIM_P(ustarhi),
        Real* recon,AMREX_ARLIM_P(reconlo),AMREX_ARLIM_P(reconhi),
        Real* mdot,AMREX_ARLIM_P(mdotlo),AMREX_ARLIM_P(mdothi),
        Real* maskdivres,
        AMREX_ARLIM_P(maskdivreslo),AMREX_ARLIM_P(maskdivreshi),
        Real* maskres,
        AMREX_ARLIM_P(maskreslo),AMREX_ARLIM_P(maskreshi),
        int* SDC_outer_sweeps,
        int* homflag,
        int* nsolve,
        int* ncomp_denold,
        int* ncomp_veldest,
        int* ncomp_dendest);


      void fort_init_physics_vars(
       int* tid,
       int* FD_curv_interp,
       Real* curv_min,
       Real* curv_max,
       int* isweep,
       int* nrefine_vof,
       Real* smoothing_length_scale,
       Real* den_interface,
       Real* visc_interface,
       Real* heatvisc_interface,
       Real* speciesvisc_interface,
       int* curv_index,
       int* pforce_index,
       int* faceden_index,
       int* facecut_index,
       int* icefacecut_index,
       int* icemask_index,
       int* facevisc_index,
       int* faceheat_index,
       int* facevel_index,
       int* facespecies_index,
       int* smoothing_index,
       int* massface_index,
       int* vofface_index,
       int* ncphys,
       Real* latent_heat,
       int* freezing_model,
       int* distribute_from_target,
       int* solidheat_flag,
       Real* microlayer_size,
       int* microlayer_substrate,
       Real* microlayer_temperature_substrate,
       int* spec_material_id,
       int* mass_fraction_id,
       Real* cavitation_vapor_density,
       int* constant_density_all_time,
       Real* time,
       Real* dt_slab,
       int* project_option,
       Real* problo,Real* probhi,
       Real* visc_coef,
       int* nten,
       const Real* xlo,const Real* dx,
       Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi),
       Real* masknbr,AMREX_ARLIM_P(masknbrlo),AMREX_ARLIM_P(masknbrhi),
       Real* xface,AMREX_ARLIM_P(xfacelo),AMREX_ARLIM_P(xfacehi),
       Real* yface,AMREX_ARLIM_P(yfacelo),AMREX_ARLIM_P(yfacehi),
       Real* zface,AMREX_ARLIM_P(zfacelo),AMREX_ARLIM_P(zfacehi),
       Real* curv,AMREX_ARLIM_P(curvlo),AMREX_ARLIM_P(curvhi), 
       Real* slope,AMREX_ARLIM_P(slopelo),AMREX_ARLIM_P(slopehi), 
       Real* denstate,
       AMREX_ARLIM_P(denstatelo),AMREX_ARLIM_P(denstatehi), 
       Real* mom_den,
       AMREX_ARLIM_P(mom_denlo),AMREX_ARLIM_P(mom_denhi), 
       Real* viscstate,AMREX_ARLIM_P(viscstatelo),AMREX_ARLIM_P(viscstatehi), 
       Real* conductstate,
       AMREX_ARLIM_P(conductstatelo),
       AMREX_ARLIM_P(conductstatehi), 
       Real* solxfab,
       AMREX_ARLIM_P(solxfablo),AMREX_ARLIM_P(solxfabhi), 
       Real* solyfab,
       AMREX_ARLIM_P(solyfablo),AMREX_ARLIM_P(solyfabhi), 
       Real* solzfab,
       AMREX_ARLIM_P(solzfablo),AMREX_ARLIM_P(solzfabhi), 
       Real* cenDeDT,AMREX_ARLIM_P(cenDeDTlo),AMREX_ARLIM_P(cenDeDThi),
       Real* cenden,AMREX_ARLIM_P(cendenlo),AMREX_ARLIM_P(cendenhi),
       Real* cenvof,AMREX_ARLIM_P(cenvoflo),AMREX_ARLIM_P(cenvofhi),
       Real* cenvisc,AMREX_ARLIM_P(cenvisclo),AMREX_ARLIM_P(cenvischi),
       Real* vol,AMREX_ARLIM_P(vollo),AMREX_ARLIM_P(volhi),
       Real* levelpc,AMREX_ARLIM_P(levelpclo),AMREX_ARLIM_P(levelpchi),
       Real* vofC,AMREX_ARLIM_P(vofClo),AMREX_ARLIM_P(vofChi),
       Real* vofF,AMREX_ARLIM_P(vofFlo),AMREX_ARLIM_P(vofFhi),
       Real* massF,AMREX_ARLIM_P(massFlo),AMREX_ARLIM_P(massFhi),
       const int* tilelo,const int* tilehi,
       const int* fablo,const int* fabhi,
       int* bfact,
       const int* presbc,
       const int* velbc,
       int* nmat,
       int* nparts,
       int* nparts_def,
       int* im_solid_map,
       int* num_curv,
       int* level,
       int* finest_level);


      void fort_build_semirefinevof(
        int* tid,
        int* tessellate,
        int* ngrow_refine,
        int* nrefine_vof,
        int* nten,
        int* spec_material_id_AMBIENT,
        int* mass_fraction_id,
        Real* cavitation_vapor_density,
        const Real* xlo,const Real* dx,
        Real* slope,AMREX_ARLIM_P(slopelo),AMREX_ARLIM_P(slopehi), 
        Real* denstate,
        AMREX_ARLIM_P(denstatelo),AMREX_ARLIM_P(denstatehi), 
        Real* mom_den,
        AMREX_ARLIM_P(mom_denlo),AMREX_ARLIM_P(mom_denhi), 
        Real* vofF,AMREX_ARLIM_P(vofFlo),AMREX_ARLIM_P(vofFhi),
        Real* massF,AMREX_ARLIM_P(massFlo),AMREX_ARLIM_P(massFhi),
        const int* tilelo,const int* tilehi,
        const int* fablo,const int* fabhi,
        int* bfact,
        int* nmat,
        int* level,int* finest_level);

      void fort_buildfacewt(
        int* facewt_iter,
        int* level,
        int* finest_level,
        int* nsolve,
        int* local_face_index,
        int* facecut_index,
        int* icefacecut_index,
        int* ncphys,
        int* nmat,
        const Real* xlo,const Real* dx,
        Real* offdiagcheck,
        AMREX_ARLIM_P(offdiagchecklo),AMREX_ARLIM_P(offdiagcheckhi), 
        Real* recon,AMREX_ARLIM_P(reconlo),AMREX_ARLIM_P(reconhi), 
        Real* cenden,AMREX_ARLIM_P(cendenlo),AMREX_ARLIM_P(cendenhi), 
        Real* cenvisc,AMREX_ARLIM_P(cenvisclo),AMREX_ARLIM_P(cenvischi), 
        Real* xfwt,AMREX_ARLIM_P(xfwtlo),AMREX_ARLIM_P(xfwthi), 
        Real* yfwt,AMREX_ARLIM_P(yfwtlo),AMREX_ARLIM_P(yfwthi), 
        Real* zfwt,AMREX_ARLIM_P(zfwtlo),AMREX_ARLIM_P(zfwthi), 
        Real* xface,AMREX_ARLIM_P(xfacelo),AMREX_ARLIM_P(xfacehi), 
        Real* yface,AMREX_ARLIM_P(yfacelo),AMREX_ARLIM_P(yfacehi), 
        Real* zface,AMREX_ARLIM_P(zfacelo),AMREX_ARLIM_P(zfacehi), 
        Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi), 
        const int* tilelo,const int* tilehi,
        const int* fablo,const int* fabhi,int* bfact,
        Real* min_face_wt,
        Real* max_face_wt,
        int* presbc,
        Real* visc_coef,
        int* uncoupled_viscosity,
        int* project_option);


     void fort_l1_distance(
       int* level,
       int* finest_level,
       Real* maskfab,AMREX_ARLIM_P(maskfablo),AMREX_ARLIM_P(maskfabhi),
       Real* lsfab,AMREX_ARLIM_P(lsfablo),AMREX_ARLIM_P(lsfabhi),
       Real* l1lsfab,AMREX_ARLIM_P(l1lsfablo),AMREX_ARLIM_P(l1lsfabhi),
       const int* tilelo,const int* tilehi,
       const int* fablo,const int* fabhi,
       int* bfact,
       const Real* xlo,const Real* dx,
       int* ngrowLS,
       int* nmat);


     void FORT_CELLFACEINIT(
       int* tid,
       int* tessellate,
       int* nten,
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
       Real* time,int* ngrow,
       int* nmat,int* ncellfrac);

     void fort_curvstrip(
       int* post_restart_flag,
       int* vof_height_function,
       int* level,
       int* finest_level,
       Real* curv_min,
       Real* curv_max,
       int* nhistory,
       Real* history_dat,
       AMREX_ARLIM_P(history_datlo),AMREX_ARLIM_P(history_dathi),
       Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi),
       Real* vol,AMREX_ARLIM_P(vollo),AMREX_ARLIM_P(volhi),
       Real* areax,AMREX_ARLIM_P(areaxlo),AMREX_ARLIM_P(areaxhi),
       Real* areay,AMREX_ARLIM_P(areaylo),AMREX_ARLIM_P(areayhi),
       Real* areaz,AMREX_ARLIM_P(areazlo),AMREX_ARLIM_P(areazhi),
       Real* maskfab,AMREX_ARLIM_P(maskfablo),AMREX_ARLIM_P(maskfabhi),
       Real* LSPC,AMREX_ARLIM_P(LSPClo),AMREX_ARLIM_P(LSPChi),
       Real* recon,AMREX_ARLIM_P(reconlo),AMREX_ARLIM_P(reconhi),
       Real* curvfab,AMREX_ARLIM_P(curvfablo),AMREX_ARLIM_P(curvfabhi),
       Real* velfab,AMREX_ARLIM_P(velfablo),AMREX_ARLIM_P(velfabhi),
       Real* denfab,AMREX_ARLIM_P(denfablo),AMREX_ARLIM_P(denfabhi),
       const int* tilelo,const int* tilehi,
       const int* fablo,const int* fabhi,
       int* bfact,
       int* bfact_grid,
       const int* rz_flag,
       const Real* xlo,const Real* dx,
       Real* time,
       Real* visc_coef,
       int* nmat,int* nten,
       int* num_curv,
       int* ngrow_distance,
       int* curv_sten);

     void fort_init_particle_container(
        Real* particles_weight_XD,
        Real* particles_weight_VEL,
        int* tid,
        int* single_particle_size,
        int* isweep,
        int* append_flag,
        int* particle_nsubdivide,
        int* particle_max_per_nsubdivide,
        int* particle_min_per_nsubdivide,
        int* nmat,
        const int* tilelo,const int* tilehi,
        const int* fablo,const int* fabhi,
        int* bfact,
        int* level,
        int* finest_level,
	Real* cur_time_slab,
        const Real* xlo,const Real* dx,
        const void* particles,
        int Np, 
        Real* new_particles, // size is "new_Pdata_size"
        int* new_Pdata_size, 
        int* Np_append, 
	int* particle_link_data,
	int* particle_delete_flag,
	int* cell_particle_count,
	AMREX_ARLIM_P(cell_particle_countlo),
	AMREX_ARLIM_P(cell_particle_counthi),
        Real* velfab,
	AMREX_ARLIM_P(velfablo),
	AMREX_ARLIM_P(velfabhi),
        Real* xd,
	AMREX_ARLIM_P(xdlo),
	AMREX_ARLIM_P(xdhi),
        Real* yd,
	AMREX_ARLIM_P(ydlo),
	AMREX_ARLIM_P(ydhi),
        Real* zd,
	AMREX_ARLIM_P(zdlo),
	AMREX_ARLIM_P(zdhi),
        Real* LS,
	AMREX_ARLIM_P(LSlo),
	AMREX_ARLIM_P(LShi),
        Real* mfiner,
	AMREX_ARLIM_P(mfinerlo),
	AMREX_ARLIM_P(mfinerhi) );


     void fort_move_particle_container(
        int* tid,
        int* single_particle_size,
        Real* particle_volume,
        Real* particle_relaxation_time_to_fluid,
        int* particle_interaction_ngrow,
        int* nmat,
        const int* tilelo,const int* tilehi,
        const int* fablo,const int* fabhi,
        int* bfact,
        int* level,
        int* finest_level,
        const Real* xlo,const Real* dx,
        const void* particles,
        int Np, 
        const void* particles_NBR,
        int Np_NBR, 
        const void* nbr_particles,
        int Np_NBR_only, 
	int* particle_link_data,
	int* cell_particle_count,
	AMREX_ARLIM_P(cell_particle_countlo),
	AMREX_ARLIM_P(cell_particle_counthi),
	Real* dt,
	Real* vel_time_slab,
        Real* umac,AMREX_ARLIM_P(umaclo),AMREX_ARLIM_P(umachi),
        Real* vmac,AMREX_ARLIM_P(vmaclo),AMREX_ARLIM_P(vmachi),
        Real* wmac,AMREX_ARLIM_P(wmaclo),AMREX_ARLIM_P(wmachi),
        Real* LS,AMREX_ARLIM_P(LSlo),AMREX_ARLIM_P(LShi),
	int* velbc,
	int* denbc,
	int* dombc,
	const int* domlo,
	const int* domhi);
  };

}/* namespace amrex */

#endif

#endif

