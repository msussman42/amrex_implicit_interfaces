#ifndef _LEVEL_F_H_
#define _LEVEL_F_H_ 

#ifdef BL_LANG_FORT
#    define FORT_CURVSTRIP curvstrip
#    define FORT_L1_DISTANCE l1_distance
#    define FORT_DOTMASK_BUILD dotmask_build
#    define FORT_CELLFACEINIT cellfaceinit
#    define FORT_NEUMANNEXTRAP neumannextrap
#    define FORT_BUILDFACEWT buildfacewt
#    define FORT_INIT_PHYSICS_VARS init_physics_vars
#    define FORT_BUILD_SEMIREFINEVOF build_semirefinevof
#    define FORT_BUILD_MODVISC build_modvisc
#    define FORT_BUILD_MASK_UNSPLIT build_mask_unsplit
#    define FORT_CELL_TO_MAC cell_to_mac
#    define FORT_MAC_TO_CELL mac_to_cell
#    define FORT_INC_TEMP inc_temp
#    define FORT_RENORMALIZE_PRESCRIBE renormalize_prescribe
#    define FORT_PURGEFLOTSAM purgeflotsam
#    define FORT_INITRECALESCE initrecalesce
#    define FORT_GETTYPEFAB gettypefab
#    define FORT_GETCOLORSUM getcolorsum
#    define FORT_LEVELRECOLOR levelrecolor
#    define FORT_COLORFILL colorfill
#    define FORT_GRIDRECOLOR gridrecolor
#    define FORT_LEVELCOLORINIT levelcolorinit
#    define FORT_AVGDOWNCOLOR avgdowncolor
#    define FORT_COPYFINECOARSECOLOR copyfinecoarsecolor

#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_CURVSTRIP CURVSTRIP
#    define FORT_L1_DISTANCE L1_DISTANCE
#    define FORT_DOTMASK_BUILD DOTMASK_BUILD
#    define FORT_CELLFACEINIT CELLFACEINIT
#    define FORT_NEUMANNEXTRAP NEUMANNEXTRAP
#    define FORT_BUILDFACEWT BUILDFACEWT
#    define FORT_INIT_PHYSICS_VARS INIT_PHYSICS_VARS
#    define FORT_BUILD_SEMIREFINEVOF BUILD_SEMIREFINEVOF
#    define FORT_BUILD_MODVISC BUILD_MODVISC
#    define FORT_BUILD_MASK_UNSPLIT BUILD_MASK_UNSPLIT
#    define FORT_CELL_TO_MAC CELL_TO_MAC
#    define FORT_MAC_TO_CELL MAC_TO_CELL
#    define FORT_INC_TEMP INC_TEMP
#    define FORT_RENORMALIZE_PRESCRIBE RENORMALIZE_PRESCRIBE
#    define FORT_PURGEFLOTSAM PURGEFLOTSAM
#    define FORT_INITSUPERELASTIC INITSUPERELASTIC
#    define FORT_INITRECALESCE INITRECALESCE
#    define FORT_GETTYPEFAB GETTYPEFAB
#    define FORT_GETCOLORSUM GETCOLORSUM
#    define FORT_LEVELRECOLOR LEVELRECOLOR
#    define FORT_COLORFILL COLORFILL
#    define FORT_GRIDRECOLOR GRIDRECOLOR
#    define FORT_LEVELCOLORINIT LEVELCOLORINIT
#    define FORT_AVGDOWNCOLOR AVGDOWNCOLOR
#    define FORT_COPYFINECOARSECOLOR COPYFINECOARSECOLOR
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_CURVSTRIP curvstrip
#    define FORT_L1_DISTANCE l1_distance
#    define FORT_DOTMASK_BUILD dotmask_build
#    define FORT_CELLFACEINIT cellfaceinit
#    define FORT_NEUMANNEXTRAP neumannextrap
#    define FORT_BUILDFACEWT buildfacewt
#    define FORT_INIT_PHYSICS_VARS init_physics_vars
#    define FORT_BUILD_SEMIREFINEVOF build_semirefinevof
#    define FORT_BUILD_MODVISC build_modvisc
#    define FORT_BUILD_MASK_UNSPLIT build_mask_unsplit
#    define FORT_CELL_TO_MAC cell_to_mac
#    define FORT_MAC_TO_CELL mac_to_cell
#    define FORT_INC_TEMP inc_temp
#    define FORT_RENORMALIZE_PRESCRIBE renormalize_prescribe
#    define FORT_PURGEFLOTSAM purgeflotsam
#    define FORT_INITRECALESCE initrecalesce
#    define FORT_GETTYPEFAB gettypefab
#    define FORT_GETCOLORSUM getcolorsum
#    define FORT_LEVELRECOLOR levelrecolor
#    define FORT_COLORFILL colorfill
#    define FORT_GRIDRECOLOR gridrecolor
#    define FORT_LEVELCOLORINIT levelcolorinit
#    define FORT_AVGDOWNCOLOR avgdowncolor
#    define FORT_COPYFINECOARSECOLOR copyfinecoarsecolor
#  else
#    define FORT_CURVSTRIP curvstrip_
#    define FORT_L1_DISTANCE l1_distance_
#    define FORT_DOTMASK_BUILD dotmask_build_
#    define FORT_CELLFACEINIT cellfaceinit_
#    define FORT_NEUMANNEXTRAP neumannextrap_
#    define FORT_BUILDFACEWT buildfacewt_
#    define FORT_INIT_PHYSICS_VARS init_physics_vars_
#    define FORT_BUILD_SEMIREFINEVOF build_semirefinevof_
#    define FORT_BUILD_MODVISC build_modvisc_
#    define FORT_BUILD_MASK_UNSPLIT build_mask_unsplit_
#    define FORT_CELL_TO_MAC cell_to_mac_
#    define FORT_MAC_TO_CELL mac_to_cell_
#    define FORT_INC_TEMP inc_temp_
#    define FORT_RENORMALIZE_PRESCRIBE renormalize_prescribe_
#    define FORT_PURGEFLOTSAM purgeflotsam_
#    define FORT_INITRECALESCE initrecalesce_
#    define FORT_GETTYPEFAB gettypefab_
#    define FORT_GETCOLORSUM getcolorsum_
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
     const Real* LS,AMREX_ARLIM_P(LSlo),AMREX_ARLIM_P(LShi),
     const Real* typefab,AMREX_ARLIM_P(typefablo),AMREX_ARLIM_P(typefabhi),
     const Real* xlo,
     const Real* dx,
     const int* tilelo, const int* tilehi,
     const int* fablo, const int* fabhi,int* bfact,
     int* type_flag,int* nmat);

     void FORT_GETCOLORSUM(
       int* operation_flag,
       int* sweep_num,
       int* tessellate,
       int* distribute_mdot_evenly,
       Real* dt,
       const Real* dx,
       const Real* xlo,
       int* nmat,
       int* nstate,
       Real* snew,AMREX_ARLIM_P(snewlo),AMREX_ARLIM_P(snewhi),
       Real* mdot,AMREX_ARLIM_P(mdotlo),AMREX_ARLIM_P(mdothi),
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
       Real* level_blobarray,
       int* level_blobtypearray,
       Real* level_mdot_array,
       Real* level_mdot_redistribute_array,
       int* arraysize,
       int* mdot_arraysize,
       int* num_elements_blobclass,
       int* ncomp_mdot_alloc,
       int* ncomp_mdot,
       int* levelbc,
       int* velbc,
       int* nface,int* nface_dst,int* ncellfrac);

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
     const int* clo,const int* chi,const int* flo,const int* fhi);
    

    void FORT_NEUMANNEXTRAP(
     const int* lo,const int* hi,
     const Real* ls,AMREX_ARLIM_P(lslo),AMREX_ARLIM_P(lshi),
     int* bc,int* dombc,
     const int* domlo,const int* domhi);


     void FORT_RENORMALIZE_PRESCRIBE(
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
       int* ngrow_distance);


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

      void FORT_CELL_TO_MAC(
        int* ncomp_xp,
        int* ncomp_xgp,
        int* simple_AMR_BC_flag,
        int* nsolveMM_FACE,
        int* num_materials_face,
        int* tileloop,
        int* dir,
        int* operation_flag,
        int* energyflag,
        Real* beta,Real* visc_coef,
        int* face_flag,
        int* interp_vel_increment_from_cell,
        int* temperature_primitive_variable,
        int* enable_spectral,
        int* fluxvel_index,
        int* fluxden_index,
        int* facevel_index,
        int* facecut_index,
        int* icefacecut_index,
        int* curv_index,
        int* conservative_tension_force,
        int* conservative_div_uu,
        int* interp_presgrad_increment_from_face, 
        int* ignore_div_up, 
        int* pforce_index,
        int* faceden_index,
        int* icemask_index,
        int* massface_index,
        int* vofface_index,
        int* ncphys,
        int* override_density,
        int* solvability_projection,
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
        Real* xfacemm,AMREX_ARLIM_P(xfacemmlo),AMREX_ARLIM_P(xfacemmhi), 
        Real* xcellmm,AMREX_ARLIM_P(xcellmmlo),AMREX_ARLIM_P(xcellmmhi), 
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
        int* num_elements_blobclass,
        int* num_colors,
        int* nten,
        int* nfacefrac,
        int* ncellfrac,
        int* project_option,
        int* SEM_upwind,
        int* SEM_advection_algorithm);

      void FORT_MAC_TO_CELL(
        int* nsolveMM_FACE,
        int* num_materials_face,
        int* ns_time_order,
        int* divu_outer_sweeps,
        int* num_divu_outer_sweeps,
        int* operation_flag,
        int* energyflag,
        int* temperature_primitive_variable,
        int* nmat,
        int* nparts,
        int* nparts_def,
        int* im_solid_map,
        Real* added_weight,
        int* nten,
        int* level,
        int* finest_level,
        int* face_flag,
        int* local_solvability_projection,
        int* project_option,
        int* enable_spectral,
        int* fluxvel_index,
        int* fluxden_index,
        int* facevel_index,
        int* facecut_index,
        int* icefacecut_index,
        int* curv_index,
        int* conservative_tension_force,
        int* conservative_div_uu,
        int* interp_presgrad_increment_from_face, 
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
        Real* maskfine,AMREX_ARLIM_P(maskfinelo),AMREX_ARLIM_P(maskfinehi),//1=fine/fine
        Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi),//1=not covered
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
        int* use_VOF_weight,
        int* nsolve,
        int* ncomp_denold,
        int* ncomp_veldest,
        int* ncomp_dendest,
        int* SEM_advection_algorithm);


      void FORT_INC_TEMP(
        Real* beta,
        int* temperature_primitive_variable,
        int* nmat,
        int* level,
        int* finest_level,
        int* ncomp_state,
        const int* tilelo,const int* tilehi,
        const int* fablo,const int* fabhi,
        Real* state,AMREX_ARLIM_P(statelo),AMREX_ARLIM_P(statehi), 
        Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi)); //1=not covered

      void FORT_INIT_PHYSICS_VARS(
       int* tid,
       int* FD_curv_interp,
       Real* curv_min,
       Real* curv_max,
       int* isweep,
       int* nrefine_vof,
       Real* den_interface,
       Real* visc_interface,
       Real* heatvisc_interface,
       Real* speciesvisc_interface,
       int* diffusionface_flag,
       int* temperatureface_flag,
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
       Real* species_evaporation_density,
       Real* cavitation_vapor_density,
       int* override_density,
       Real* time,
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
       Real* denstate,AMREX_ARLIM_P(denstatelo),AMREX_ARLIM_P(denstatehi), 
       Real* viscstate,AMREX_ARLIM_P(viscstatelo),AMREX_ARLIM_P(viscstatehi), 
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
       Real* modvisc,AMREX_ARLIM_P(modvisclo),AMREX_ARLIM_P(modvischi),
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


      void FORT_BUILD_SEMIREFINEVOF(
        int* tid,
        int* tessellate,
        int* ngrow_refine,
        int* nrefine_vof,
        int* nrefine_cen,
        int* nten,
        int* spec_material_id_AMBIENT,
        int* mass_fraction_id,
        Real* species_evaporation_density,
        Real* cavitation_vapor_density,
        int* override_density,
        const Real* xlo,const Real* dx,
        Real* slope,AMREX_ARLIM_P(slopelo),AMREX_ARLIM_P(slopehi), 
        Real* denstate,AMREX_ARLIM_P(denstatelo),AMREX_ARLIM_P(denstatehi), 
        Real* vofF,AMREX_ARLIM_P(vofFlo),AMREX_ARLIM_P(vofFhi),
        Real* cenF,AMREX_ARLIM_P(cenFlo),AMREX_ARLIM_P(cenFhi),
        Real* massF,AMREX_ARLIM_P(massFlo),AMREX_ARLIM_P(massFhi),
        const int* tilelo,const int* tilehi,
        const int* fablo,const int* fabhi,
        int* bfact,
        int* nmat,
        int* level,int* finest_level);

      void FORT_BUILD_MASK_UNSPLIT(
        int* unsplit_flag,
        const Real* xlo,const Real* dx,
        Real* maskunsplit,AMREX_ARLIM_P(maskunsplitlo),AMREX_ARLIM_P(maskunsplithi), 
        Real* vofls0,AMREX_ARLIM_P(vofls0lo),AMREX_ARLIM_P(vofls0hi), 
        const int* tilelo,const int* tilehi,
        const int* fablo,const int* fabhi,
        int* bfact,
        int* nmat,
        int* level,int* finest_level);

      void FORT_BUILD_MODVISC(
        int* ngrow_visc,
        Real* time,
        Real* problo,Real* probhi,
        Real* visc_coef,
        int* nten,
        const Real* xlo,const Real* dx,
        Real* slope,AMREX_ARLIM_P(slopelo),AMREX_ARLIM_P(slopehi), 
        Real* denstate,AMREX_ARLIM_P(denstatelo),AMREX_ARLIM_P(denstatehi), 
        Real* viscstate,AMREX_ARLIM_P(viscstatelo),AMREX_ARLIM_P(viscstatehi), 
        Real* levelpc,AMREX_ARLIM_P(levelpclo),AMREX_ARLIM_P(levelpchi),
        Real* modvisc,AMREX_ARLIM_P(modvisclo),AMREX_ARLIM_P(modvischi),
        const int* tilelo,const int* tilehi,
        const int* fablo,const int* fabhi,
        int* bfact,
        int* nmat,
        int* level,int* finest_level);


      void FORT_BUILDFACEWT(
        int* facewt_iter,
        int* num_materials_face,
        int* level,
        int* finest_level,
        int* nsolve,
        int* nsolveMM,
        int* nsolveMM_FACE,
        int* nfacefrac,
        int* ncellfrac,
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
        Real* cellmm,AMREX_ARLIM_P(cellmmlo),AMREX_ARLIM_P(cellmmhi), 
        Real* xfacemm,AMREX_ARLIM_P(xfacemmlo),AMREX_ARLIM_P(xfacemmhi), 
        Real* yfacemm,AMREX_ARLIM_P(yfacemmlo),AMREX_ARLIM_P(yfacemmhi), 
        Real* zfacemm,AMREX_ARLIM_P(zfacemmlo),AMREX_ARLIM_P(zfacemmhi), 
        Real* xfwt,AMREX_ARLIM_P(xfwtlo),AMREX_ARLIM_P(xfwthi), 
        Real* yfwt,AMREX_ARLIM_P(yfwtlo),AMREX_ARLIM_P(yfwthi), 
        Real* zfwt,AMREX_ARLIM_P(zfwtlo),AMREX_ARLIM_P(zfwthi), 
        Real* xface,AMREX_ARLIM_P(xfacelo),AMREX_ARLIM_P(xfacehi), 
        Real* yface,AMREX_ARLIM_P(yfacelo),AMREX_ARLIM_P(yfacehi), 
        Real* zface,AMREX_ARLIM_P(zfacelo),AMREX_ARLIM_P(zfacehi), 
        Real* alt_xface,AMREX_ARLIM_P(alt_xfacelo),AMREX_ARLIM_P(alt_xfacehi), 
        Real* alt_yface,AMREX_ARLIM_P(alt_yfacelo),AMREX_ARLIM_P(alt_yfacehi), 
        Real* alt_zface,AMREX_ARLIM_P(alt_zfacelo),AMREX_ARLIM_P(alt_zfacehi), 
        Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi), 
        const int* tilelo,const int* tilehi,
        const int* fablo,const int* fabhi,int* bfact,
        Real* min_face_wt,
        Real* max_face_wt,
        int* singular_possible,
        int* solvability_projection,
        int* solvability_level_flag,
        int* presbc,
        Real* visc_coef,
        int* constant_viscosity,
        int* project_option);


     void FORT_L1_DISTANCE(
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


     void FORT_DOTMASK_BUILD(
       int* num_materials_face,
       int* level,
       int* finest_level,
       Real* dotmask,AMREX_ARLIM_P(dotmasklo),AMREX_ARLIM_P(dotmaskhi),
       Real* maskfab,AMREX_ARLIM_P(maskfablo),AMREX_ARLIM_P(maskfabhi),
       Real* vofrecon,AMREX_ARLIM_P(vofreconlo),AMREX_ARLIM_P(vofreconhi),
       const int* tilelo,const int* tilehi,
       const int* fablo,const int* fabhi,
       int* bfact,
       const Real* xlo,const Real* dx,
       Real* time,
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

     void FORT_CURVSTRIP(
       int* post_restart_flag,
       int* conservative_tension_force,
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
       Real* LSHO,AMREX_ARLIM_P(LSHOlo),AMREX_ARLIM_P(LSHOhi),
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


     void fort_assimilate_lvlset_from_particles(
       Real* particles_weight,
       int* tid,
       int* im_PLS,
       int* isweep,
       int* level,
       int* finest_level,
       Real* solid_time,
       const int* tilelo,const int* tilehi,
       const int* fablo,const int* fabhi,int* bfact,
       Real* vofnew,AMREX_ARLIM_P(vofnewlo),AMREX_ARLIM_P(vofnewhi),
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
       int* ngrow_distance,
       const void* particles,
       const void* nbr_particles,
       int Np, 
       int Nn, 
       int* matrix_points,
       int* RHS_points,
       int* ncomp_accumulate,
       Real* matrixfab,
       AMREX_ARLIM_P(matrixfablo),
       AMREX_ARLIM_P(matrixfabhi));


     void fort_init_particle_container(
        Real* particles_weight,
        int* tid,
        int* single_particle_size,
        int* isweep,
        int* append_flag,
        int* particle_nsubdivide,
        int* particle_max_per_nsubdivide,
        int* particleLS_flag,
        int* imPLS,
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
        Real* xdisplace,
	AMREX_ARLIM_P(xdisplacelo),
	AMREX_ARLIM_P(xdisplacehi),
        Real* LS,AMREX_ARLIM_P(LSlo),AMREX_ARLIM_P(LShi));


     void fort_move_particle_container(
        int* tid,
        int* imPLS,
        int* nmat,
        const int* tilelo,const int* tilehi,
        const int* fablo,const int* fabhi,
        int* bfact,
        int* level,
        int* finest_level,
        const Real* xlo,const Real* dx,
        const void* particles,
        int Np, 
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

