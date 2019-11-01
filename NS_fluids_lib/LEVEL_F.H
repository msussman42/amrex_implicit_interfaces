#ifndef _LEVEL_F_H_
#define _LEVEL_F_H_ 

#ifdef BL_LANG_FORT
#    define FORT_CURVSTRIP curvstrip
#    define FORT_L1_DISTANCE l1_distance
#    define FORT_DOTMASK_BUILD dotmask_build
#    define FORT_FACEPROCESS faceprocess
#    define FORT_CELLFACEINIT cellfaceinit
#    define FORT_NEUMANNEXTRAP neumannextrap
#    define FORT_BUILDFACEWT buildfacewt
#    define FORT_INIT_PHYSICS_VARS init_physics_vars
#    define FORT_BUILD_SEMIREFINEVOF build_semirefinevof
#    define FORT_BUILD_MODVISC build_modvisc
#    define FORT_BUILD_MASK_UNSPLIT build_mask_unsplit
#    define FORT_CELL_TO_MAC cell_to_mac
#    define FORT_GFMUPDATE gfmupdate
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
#    define FORT_FACEPROCESS FACEPROCESS
#    define FORT_CELLFACEINIT CELLFACEINIT
#    define FORT_NEUMANNEXTRAP NEUMANNEXTRAP
#    define FORT_BUILDFACEWT BUILDFACEWT
#    define FORT_INIT_PHYSICS_VARS INIT_PHYSICS_VARS
#    define FORT_BUILD_SEMIREFINEVOF BUILD_SEMIREFINEVOF
#    define FORT_BUILD_MODVISC BUILD_MODVISC
#    define FORT_BUILD_MASK_UNSPLIT BUILD_MASK_UNSPLIT
#    define FORT_CELL_TO_MAC CELL_TO_MAC
#    define FORT_GFMUPDATE GFMUPDATE
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
#    define FORT_FACEPROCESS faceprocess
#    define FORT_CELLFACEINIT cellfaceinit
#    define FORT_NEUMANNEXTRAP neumannextrap
#    define FORT_BUILDFACEWT buildfacewt
#    define FORT_INIT_PHYSICS_VARS init_physics_vars
#    define FORT_BUILD_SEMIREFINEVOF build_semirefinevof
#    define FORT_BUILD_MODVISC build_modvisc
#    define FORT_BUILD_MASK_UNSPLIT build_mask_unsplit
#    define FORT_CELL_TO_MAC cell_to_mac
#    define FORT_GFMUPDATE gfmupdate
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
#    define FORT_FACEPROCESS faceprocess_
#    define FORT_CELLFACEINIT cellfaceinit_
#    define FORT_NEUMANNEXTRAP neumannextrap_
#    define FORT_BUILDFACEWT buildfacewt_
#    define FORT_INIT_PHYSICS_VARS init_physics_vars_
#    define FORT_BUILD_SEMIREFINEVOF build_semirefinevof_
#    define FORT_BUILD_MODVISC build_modvisc_
#    define FORT_BUILD_MASK_UNSPLIT build_mask_unsplit_
#    define FORT_CELL_TO_MAC cell_to_mac_
#    define FORT_GFMUPDATE gfmupdate_
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

#include <ArrayLim.H>

  extern "C" {

    void FORT_GETTYPEFAB(
     const Real* LS,ARLIM_P(LSlo),ARLIM_P(LShi),
     const Real* typefab,ARLIM_P(typefablo),ARLIM_P(typefabhi),
     const Real* xlo,
     const Real* dx,
     const int* tilelo, const int* tilehi,
     const int* fablo, const int* fabhi,int* bfact,
     int* type_flag,int* nmat);

     void FORT_GETCOLORSUM(
       int* sweep_num,
       const Real* dx,
       const Real* xlo,
       int* nmat,
       Real* LS, ARLIM_P(LSlo),ARLIM_P(LShi),
       Real* VEL, ARLIM_P(VELlo),ARLIM_P(VELhi),
       Real* DEN, ARLIM_P(DENlo),ARLIM_P(DENhi),
       Real* VOF, ARLIM_P(VOFlo),ARLIM_P(VOFhi),
       Real* facefab, ARLIM_P(facefablo),ARLIM_P(facefabhi),
       Real* xface, ARLIM_P(xfacelo),ARLIM_P(xfacehi),
       Real* yface, ARLIM_P(yfacelo),ARLIM_P(yfacehi),
       Real* zface, ARLIM_P(zfacelo),ARLIM_P(zfacehi),
       Real* areax, ARLIM_P(areaxlo),ARLIM_P(areaxhi),
       Real* areay, ARLIM_P(areaylo),ARLIM_P(areayhi),
       Real* areaz, ARLIM_P(areazlo),ARLIM_P(areazhi),
       Real* cellfab, ARLIM_P(cellfablo),ARLIM_P(cellfabhi),
       Real* typefab, ARLIM_P(typefablo),ARLIM_P(typefabhi),
       Real* color, ARLIM_P(colorlo),ARLIM_P(colorhi),
       Real* mask, ARLIM_P(masklo),ARLIM_P(maskhi),
       const int* tilelo, const int* tilehi,
       const int* fablo, const int* fabhi, 
       int* bfact,
       int* level,
       int* finest_level,
       int* rzflag,
       int* num_colors,
       Real* cum_blobarray,
       Real* level_blobarray,
       int* level_blobtypearray,
       int* arraysize,int* num_elements_blobclass,
       int* levelbc,
       int* velbc,
       int* nface,int* nface_dst,int* ncellfrac);

    void FORT_LEVELRECOLOR(
     Real* color,ARLIM_P(colorlo),ARLIM_P(colorhi),
     const Real* xlo,const Real* dx,
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,int* bfact,
     int* domaincolormap,
     int* max_colors_level,int* level,int* base_level,int* arrsize);

    void FORT_GRIDRECOLOR(
     Real* mask,ARLIM_P(masklo),ARLIM_P(maskhi),
     Real* color,ARLIM_P(colorlo),ARLIM_P(colorhi),
     const Real* xlo,const Real* dx,
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,
     int* bfact,
     int* levelcolormap,
     int* max_colors_grid,
     int* number_grids,int* arrsize);

    void FORT_COLORFILL(
     const Real* mask,ARLIM_P(masklo),ARLIM_P(maskhi),
     const Real* typefab,ARLIM_P(typefablo),ARLIM_P(typefabhi),
     const Real* color,ARLIM_P(colorlo),ARLIM_P(colorhi),
     const int* ijk,ARLIM_P(ijklo),ARLIM_P(ijkhi),
     const int* lo,const int* hi,int* ipass,int* number_grids,
     int* color_per_grid,const int* gridno,int* max_colors_grid,
     int* typedim);

    void FORT_LEVELCOLORINIT(
     Real* mask,ARLIM_P(masklo),ARLIM_P(maskhi),
     Real* color,ARLIM_P(colorlo),ARLIM_P(colorhi),
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
     const Real* crse,ARLIM_P(crselo),ARLIM_P(crsehi),
     const Real* fine,ARLIM_P(fine_lo),ARLIM_P(fine_hi),
     const Real* typef,ARLIM_P(typef_lo),ARLIM_P(typef_hi),
     const Real* typec,ARLIM_P(typec_lo),ARLIM_P(typec_hi),
     const int* clo,const int* chi);

    void FORT_COPYFINECOARSECOLOR(
     const Real* prob_lo,
     const Real* dxf,
     int* bfact_f,int* bfact,
     const Real* xlo_fine,
     const Real* dx,
     const Real* crse,ARLIM_P(crselo),ARLIM_P(crsehi),
     const Real* fine,ARLIM_P(fine_lo),ARLIM_P(fine_hi),
     const Real* typef,ARLIM_P(typef_lo),ARLIM_P(typef_hi),
     const Real* maskf,ARLIM_P(maskf_lo),ARLIM_P(maskf_hi),
     const int* clo,const int* chi,const int* flo,const int* fhi);
    

    void FORT_NEUMANNEXTRAP(
     const int* lo,const int* hi,
     const Real* ls,ARLIM_P(lslo),ARLIM_P(lshi),
     int* bc,int* dombc,
     const int* domlo,const int* domhi);


     void FORT_RENORMALIZE_PRESCRIBE(
       int* tid,
       int* level,int* finest_level,
       Real* solid_time,
       const int* tilelo,const int* tilehi,
       const int* fablo,const int* fabhi,int* bfact,
       Real* vofnew,ARLIM_P(vofnewlo),ARLIM_P(vofnewhi),
       Real* solid,ARLIM_P(solidlo),ARLIM_P(solidhi),
       Real* maskcov,ARLIM_P(maskcovlo),ARLIM_P(maskcovhi),
       Real* LS,ARLIM_P(LSlo),ARLIM_P(LShi),
       Real* mof,ARLIM_P(moflo),ARLIM_P(mofhi),
       Real* den,ARLIM_P(denlo),ARLIM_P(denhi),
       Real* vel,ARLIM_P(vello),ARLIM_P(velhi),
       Real* velnew,ARLIM_P(velnewlo),ARLIM_P(velnewhi),
       Real* dennew,ARLIM_P(dennewlo),ARLIM_P(dennewhi),
       Real* lsnew,ARLIM_P(lsnewlo),ARLIM_P(lsnewhi),
       const Real* xlo,const Real* dx, 
       Real* time,
       int* nmat,
       int* nparts,
       int* nparts_def,
       int* im_solid_map,
       int* renormalize_only,
       int* solidheat_flag);


     void FORT_PURGEFLOTSAM(
       Real* delta_mass,
       int* truncate_volume_fractions,
       Real* truncate_thickness,
       int* level,int* finest_level,
       Real* time,
       const int* tilelo,const int* tilehi,
       const int* fablo,const int* fabhi,int* bfact,
       Real* maskcov,ARLIM_P(maskcovlo),ARLIM_P(maskcovhi),
       Real* vofnew,ARLIM_P(vofnewlo),ARLIM_P(vofnewhi),
       Real* LS,ARLIM_P(LSlo),ARLIM_P(LShi),
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
        int* pforce_index,
        int* faceden_index,
        int* icemask_index,
        int* massface_index,
        int* vofface_index,
        int* ncphys,
        int* make_interface_incomp,
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
        Real* semflux,ARLIM_P(semfluxlo),ARLIM_P(semfluxhi),
        Real* maskfine,ARLIM_P(maskfinelo),ARLIM_P(maskfinehi),
        Real* maskcov,ARLIM_P(maskcovlo),ARLIM_P(maskcovhi),
        Real* maskSEM,ARLIM_P(maskSEMlo),ARLIM_P(maskSEMhi),
        Real* levelPC,ARLIM_P(levelPClo),ARLIM_P(levelPChi),
        Real* sol,ARLIM_P(sollo),ARLIM_P(solhi),
        Real* xcut,ARLIM_P(xcutlo),ARLIM_P(xcuthi), 
        Real* xface,ARLIM_P(xfacelo),ARLIM_P(xfacehi), 
        Real* xfacemm,ARLIM_P(xfacemmlo),ARLIM_P(xfacemmhi), 
        Real* xcellmm,ARLIM_P(xcellmmlo),ARLIM_P(xcellmmhi), 
        Real* recon,ARLIM_P(reconlo),ARLIM_P(reconhi), 
        Real* xgp,ARLIM_P(xgplo),ARLIM_P(xgphi), 
        Real* xp,ARLIM_P(xplo),ARLIM_P(xphi), 
        Real* xvel,ARLIM_P(xvello),ARLIM_P(xvelhi), 
	Real* vel,ARLIM_P(vello),ARLIM_P(velhi),
	Real* pres,ARLIM_P(preslo),ARLIM_P(preshi),
	Real* den,ARLIM_P(denlo),ARLIM_P(denhi),
	Real* mgoni,ARLIM_P(mgonilo),ARLIM_P(mgonihi),
	Real* colorfab,ARLIM_P(colorfablo),ARLIM_P(colorfabhi),
	Real* typefab,ARLIM_P(typefablo),ARLIM_P(typefabhi),
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
        Real* prescribed_solid_scale,
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

      void FORT_GFMUPDATE(
        int* tid,
        int* dirmac,
	Real* gfm_error,
        int* velbc,
        int* slab_step,
        Real* dt,
        Real* time,
        const Real* xlo,const Real* dx,
        Real* masknbr,ARLIM_P(masknbrlo),ARLIM_P(masknbrhi),
        Real* maskcov,ARLIM_P(maskcovlo),ARLIM_P(maskcovhi),
        Real* levelPC,ARLIM_P(levelPClo),ARLIM_P(levelPChi),
        Real* sol,ARLIM_P(sollo),ARLIM_P(solhi),
        Real* recon,ARLIM_P(reconlo),ARLIM_P(reconhi), 
        Real* xvel,ARLIM_P(xvello),ARLIM_P(xvelhi), 
	Real* xsrc,ARLIM_P(xsrclo),ARLIM_P(xsrchi),
	Real* ysrc,ARLIM_P(ysrclo),ARLIM_P(ysrchi),
	Real* zsrc,ARLIM_P(zsrclo),ARLIM_P(zsrchi),
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
        Real* prescribed_solid_scale,
        int* prescribed_solid_method,
        int* project_option);

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
        Real* prescribed_solid_scale,
        Real* added_weight,
        int* nten,
        int* level,
        int* finest_level,
        int* face_flag,
        int* make_interface_incomp,
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
        Real* xp,ARLIM_P(xplo),ARLIM_P(xphi), 
        Real* yp,ARLIM_P(yplo),ARLIM_P(yphi), 
        Real* zp,ARLIM_P(zplo),ARLIM_P(zphi), 
        Real* xvel,ARLIM_P(xvello),ARLIM_P(xvelhi), 
        Real* yvel,ARLIM_P(yvello),ARLIM_P(yvelhi), 
        Real* zvel,ARLIM_P(zvello),ARLIM_P(zvelhi), 
        Real* xface,ARLIM_P(xfacelo),ARLIM_P(xfacehi), 
        Real* yface,ARLIM_P(yfacelo),ARLIM_P(yfacehi), 
        Real* zface,ARLIM_P(zfacelo),ARLIM_P(zfacehi), 
        Real* ax,ARLIM_P(axlo),ARLIM_P(axhi), 
        Real* ay,ARLIM_P(aylo),ARLIM_P(ayhi), 
        Real* az,ARLIM_P(azlo),ARLIM_P(azhi), 
        Real* vol,ARLIM_P(vollo),ARLIM_P(volhi), 
        Real* rhs,ARLIM_P(rhslo),ARLIM_P(rhshi), 
        Real* veldest,ARLIM_P(veldestlo),ARLIM_P(veldesthi), 
        Real* dendest,ARLIM_P(dendestlo),ARLIM_P(dendesthi), 
        Real* maskfine,ARLIM_P(maskfinelo),ARLIM_P(maskfinehi),//1=fine/fine
        Real* maskcov,ARLIM_P(maskcovlo),ARLIM_P(maskcovhi),//1=not covered
        Real* maskSEM,ARLIM_P(maskSEMlo),ARLIM_P(maskSEMhi),
        Real* levelPC,ARLIM_P(levelPClo),ARLIM_P(levelPChi),
        Real* sol,ARLIM_P(sollo),ARLIM_P(solhi),
        Real* cterm,ARLIM_P(ctermlo),ARLIM_P(ctermhi),
        Real* pold,ARLIM_P(poldlo),ARLIM_P(poldhi),
        Real* denold,ARLIM_P(denoldlo),ARLIM_P(denoldhi),
        Real* ustar,ARLIM_P(ustarlo),ARLIM_P(ustarhi),
        Real* recon,ARLIM_P(reconlo),ARLIM_P(reconhi),
        Real* mdot,ARLIM_P(mdotlo),ARLIM_P(mdothi),
        Real* maskdivres,
        ARLIM_P(maskdivreslo),ARLIM_P(maskdivreshi),
        Real* maskres,
        ARLIM_P(maskreslo),ARLIM_P(maskreshi),
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
        Real* state,ARLIM_P(statelo),ARLIM_P(statehi), 
        Real* maskcov,ARLIM_P(maskcovlo),ARLIM_P(maskcovhi)); //1=not covered

      void FORT_INIT_PHYSICS_VARS(
       int* tid,
       int* FD_curv_select,
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
       int* cavitation_species,
       int* override_density,
       Real* time,
       int* project_option,
       Real* problo,Real* probhi,
       Real* visc_coef,
       int* nten,
       const Real* xlo,const Real* dx,
       Real* maskcov,ARLIM_P(maskcovlo),ARLIM_P(maskcovhi),
       Real* masknbr,ARLIM_P(masknbrlo),ARLIM_P(masknbrhi),
       Real* xface,ARLIM_P(xfacelo),ARLIM_P(xfacehi),
       Real* yface,ARLIM_P(yfacelo),ARLIM_P(yfacehi),
       Real* zface,ARLIM_P(zfacelo),ARLIM_P(zfacehi),
       Real* curv,ARLIM_P(curvlo),ARLIM_P(curvhi), 
       Real* slope,ARLIM_P(slopelo),ARLIM_P(slopehi), 
       Real* denstate,ARLIM_P(denstatelo),ARLIM_P(denstatehi), 
       Real* viscstate,ARLIM_P(viscstatelo),ARLIM_P(viscstatehi), 
       Real* sol,ARLIM_P(sollo),ARLIM_P(solhi), 
       Real* cenDeDT,ARLIM_P(cenDeDTlo),ARLIM_P(cenDeDThi),
       Real* cenden,ARLIM_P(cendenlo),ARLIM_P(cendenhi),
       Real* cenvof,ARLIM_P(cenvoflo),ARLIM_P(cenvofhi),
       Real* cenvisc,ARLIM_P(cenvisclo),ARLIM_P(cenvischi),
       Real* vol,ARLIM_P(vollo),ARLIM_P(volhi),
       Real* levelpc,ARLIM_P(levelpclo),ARLIM_P(levelpchi),
       Real* vofC,ARLIM_P(vofClo),ARLIM_P(vofChi),
       Real* vofF,ARLIM_P(vofFlo),ARLIM_P(vofFhi),
       Real* massF,ARLIM_P(massFlo),ARLIM_P(massFhi),
       Real* modvisc,ARLIM_P(modvisclo),ARLIM_P(modvischi),
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
       Real* prescribed_solid_scale,
       int* level,
       int* finest_level);


      void FORT_BUILD_SEMIREFINEVOF(
        int* tid,
        int* tessellate,
        int* ngrow_refine,
        int* nrefine_vof,
        int* nrefine_cen,
        int* nten,
        int* spec_material_id,
        int* mass_fraction_id,
        Real* species_evaporation_density,
        Real* cavitation_vapor_density,
        int* cavitation_species,
        int* override_density,
        const Real* xlo,const Real* dx,
        Real* slope,ARLIM_P(slopelo),ARLIM_P(slopehi), 
        Real* denstate,ARLIM_P(denstatelo),ARLIM_P(denstatehi), 
        Real* vofF,ARLIM_P(vofFlo),ARLIM_P(vofFhi),
        Real* cenF,ARLIM_P(cenFlo),ARLIM_P(cenFhi),
        Real* massF,ARLIM_P(massFlo),ARLIM_P(massFhi),
        const int* tilelo,const int* tilehi,
        const int* fablo,const int* fabhi,
        int* bfact,
        int* nmat,
        int* level,int* finest_level);

      void FORT_BUILD_MASK_UNSPLIT(
        int* unsplit_flag,
        int* make_interface_incomp,
        const Real* xlo,const Real* dx,
        Real* maskunsplit,ARLIM_P(maskunsplitlo),ARLIM_P(maskunsplithi), 
        Real* vofls0,ARLIM_P(vofls0lo),ARLIM_P(vofls0hi), 
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
        Real* slope,ARLIM_P(slopelo),ARLIM_P(slopehi), 
        Real* denstate,ARLIM_P(denstatelo),ARLIM_P(denstatehi), 
        Real* viscstate,ARLIM_P(viscstatelo),ARLIM_P(viscstatehi), 
        Real* levelpc,ARLIM_P(levelpclo),ARLIM_P(levelpchi),
        Real* modvisc,ARLIM_P(modvisclo),ARLIM_P(modvischi),
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
	ARLIM_P(offdiagchecklo),ARLIM_P(offdiagcheckhi), 
        Real* recon,ARLIM_P(reconlo),ARLIM_P(reconhi), 
        Real* cenden,ARLIM_P(cendenlo),ARLIM_P(cendenhi), 
        Real* cenvisc,ARLIM_P(cenvisclo),ARLIM_P(cenvischi), 
        Real* cellmm,ARLIM_P(cellmmlo),ARLIM_P(cellmmhi), 
        Real* xfacemm,ARLIM_P(xfacemmlo),ARLIM_P(xfacemmhi), 
        Real* yfacemm,ARLIM_P(yfacemmlo),ARLIM_P(yfacemmhi), 
        Real* zfacemm,ARLIM_P(zfacemmlo),ARLIM_P(zfacemmhi), 
        Real* xfwt,ARLIM_P(xfwtlo),ARLIM_P(xfwthi), 
        Real* yfwt,ARLIM_P(yfwtlo),ARLIM_P(yfwthi), 
        Real* zfwt,ARLIM_P(zfwtlo),ARLIM_P(zfwthi), 
        Real* xface,ARLIM_P(xfacelo),ARLIM_P(xfacehi), 
        Real* yface,ARLIM_P(yfacelo),ARLIM_P(yfacehi), 
        Real* zface,ARLIM_P(zfacelo),ARLIM_P(zfacehi), 
        Real* alt_xface,ARLIM_P(alt_xfacelo),ARLIM_P(alt_xfacehi), 
        Real* alt_yface,ARLIM_P(alt_yfacelo),ARLIM_P(alt_yfacehi), 
        Real* alt_zface,ARLIM_P(alt_zfacelo),ARLIM_P(alt_zfacehi), 
        Real* mask,ARLIM_P(masklo),ARLIM_P(maskhi), 
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
        Real* prescribed_solid_scale,
        int* project_option);


     void FORT_L1_DISTANCE(
       int* level,
       int* finest_level,
       Real* maskfab,ARLIM_P(maskfablo),ARLIM_P(maskfabhi),
       Real* lsfab,ARLIM_P(lsfablo),ARLIM_P(lsfabhi),
       Real* l1lsfab,ARLIM_P(l1lsfablo),ARLIM_P(l1lsfabhi),
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
       Real* dotmask,ARLIM_P(dotmasklo),ARLIM_P(dotmaskhi),
       Real* maskfab,ARLIM_P(maskfablo),ARLIM_P(maskfabhi),
       Real* vofrecon,ARLIM_P(vofreconlo),ARLIM_P(vofreconhi),
       const int* tilelo,const int* tilehi,
       const int* fablo,const int* fabhi,
       int* bfact,
       const Real* xlo,const Real* dx,
       Real* time,
       int* nmat);

     void FORT_FACEPROCESS(
       int* tid,
       int* dir,
       int* tessellate,
       int* level,
       int* finest_level,
       Real* dstfab,ARLIM_P(dstfablo),ARLIM_P(dstfabhi),
       Real* facefab,ARLIM_P(facefablo),ARLIM_P(facefabhi),
       Real* vofrecon,ARLIM_P(vofreconlo),ARLIM_P(vofreconhi),
       const int* tilelo,const int* tilehi,
       const int* fablo,const int* fabhi,
       int* bfact,
       const int* rz_flag,
       const Real* xlo,const Real* dx,
       Real* time,
       int* nmat,int* nface_src,int* nface_dst);

     void FORT_CELLFACEINIT(
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
       Real* time,int* ngrow,
       int* nmat,int* ncellfrac);

     void FORT_CURVSTRIP(
       int* post_restart_flag,
       int* conservative_tension_force,
       int* level,
       int* finest_level,
       Real* curv_min,
       Real* curv_max,
       Real* maskcov, ARLIM_P(maskcovlo),ARLIM_P(maskcovhi),
       Real* vol, ARLIM_P(vollo),ARLIM_P(volhi),
       Real* areax, ARLIM_P(areaxlo),ARLIM_P(areaxhi),
       Real* areay, ARLIM_P(areaylo),ARLIM_P(areayhi),
       Real* areaz, ARLIM_P(areazlo),ARLIM_P(areazhi),
       Real* maskfab,ARLIM_P(maskfablo),ARLIM_P(maskfabhi),
       Real* LSPC,ARLIM_P(LSPClo),ARLIM_P(LSPChi),
       Real* LSHO,ARLIM_P(LSHOlo),ARLIM_P(LSHOhi),
       Real* curvfab,ARLIM_P(curvfablo),ARLIM_P(curvfabhi),
       Real* velfab,ARLIM_P(velfablo),ARLIM_P(velfabhi),
       Real* denfab,ARLIM_P(denfablo),ARLIM_P(denfabhi),
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

  };
#endif

#endif

