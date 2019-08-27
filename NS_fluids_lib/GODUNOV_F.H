#ifndef _GODUNOV_F_H_
#define _GODUNOV_F_H_

#ifdef BL_LANG_FORT
#    define FORT_SOD_SANITY   sod_sanity
#    define FORT_AGGRESSIVE   aggressive
#    define FORT_DENCOR       dencor
#    define FORT_STEFANSOLVER stefansolver
#    define FORT_HEATSOURCE_FACE heatsource_face
#    define FORT_WALLFUNCTION wallfunction
#    define FORT_INITJUMPTERM initjumpterm
#    define FORT_INIT_ICEMASK init_icemask
#    define FORT_RENORM_MDOT renorm_mdot
#    define FORT_TAGEXPANSION tagexpansion
#    define FORT_DISTRIBUTEEXPANSION distributeexpansion
#    define FORT_CLEAREXPANSION clearexpansion
#    define FORT_VFRAC_SPLIT  vfrac_split
#    define FORT_VFRAC_UNSPLIT vfrac_unsplit
#    define FORT_ESTDT        estdt
#    define FORT_CROSSTERM    crossterm
#    define FORT_FACE_GRADIENTS face_gradients
#    define FORT_COMBINEVEL   combinevel
#    define FORT_COMBINEVELFACE combinevelface
#    define FORT_VISCTENSORHEAT visctensorheat
#    define FORT_VELADVANCE  veladvance
#    define FORT_VELMAC_OVERRIDE velmac_override
#    define FORT_UPDATETENSOR updatetensor
#    define FORT_FIX_HOOP_TENSOR fix_hoop_tensor
#    define FORT_EXTRAPTENSOR extraptensor
#    define FORT_BUILD_MOMENT build_moment
#    define FORT_COPY_VEL_ON_SIGN copy_vel_on_sign
#    define FORT_MAKETENSOR maketensor
#    define FORT_TENSORFORCE tensorforce
#    define FORT_TENSORHEAT tensorheat
#    define FORT_HEATSOURCE heatsource
#    define FORT_MARANGONIFORCE marangoniforce
#    define FORT_SEMDELTAFORCE semdeltaforce
#    define FORT_SEMDELTAFORCE_FACE semdeltaforce_face
#    define FORT_UPDATESEMFORCE updatesemforce
#    define FORT_UPDATESEMFORCE_FACE updatesemforce_face
#    define FORT_SDC_TIME_QUAD sdc_time_quad
#    define FORT_SDC_TIME_QUAD_FACE sdc_time_quad_face
#    define FORT_BUILD_MASKSEM build_masksem
#    define FORT_BUILD_SLOPES build_slopes
#    define FORT_BUILD_SLOPES_FACE build_slopes_face
#    define FORT_BUILD_CONSERVE build_conserve
#    define FORT_BUILD_MACVOF build_macvof
#    define FORT_BUILD_NEWMAC build_newmac
#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_SOD_SANITY   SOD_SANITY
#    define FORT_AGGRESSIVE   AGGRESSIVE
#    define FORT_DENCOR       DENCOR
#    define FORT_STEFANSOLVER STEFANSOLVER
#    define FORT_HEATSOURCE_FACE HEATSOURCE_FACE
#    define FORT_WALLFUNCTION WALLFUNCTION
#    define FORT_INITJUMPTERM INITJUMPTERM
#    define FORT_INIT_ICEMASK INIT_ICEMASK
#    define FORT_RENORM_MDOT RENORM_MDOT
#    define FORT_TAGEXPANSION TAGEXPANSION
#    define FORT_DISTRIBUTEEXPANSION DISTRIBUTEEXPANSION
#    define FORT_CLEAREXPANSION CLEAREXPANSION
#    define FORT_VFRAC_SPLIT  VFRAC_SPLIT
#    define FORT_VFRAC_UNSPLIT VFRAC_UNSPLIT
#    define FORT_ESTDT        ESTDT
#    define FORT_CROSSTERM    CROSSTERM
#    define FORT_FACE_GRADIENTS FACE_GRADIENTS
#    define FORT_COMBINEVEL   COMBINEVEL
#    define FORT_COMBINEVELFACE COMBINEVELFACE
#    define FORT_VISCTENSORHEAT     VISCTENSORHEAT
#    define FORT_VELADVANCE  VELADVANCE
#    define FORT_VELMAC_OVERRIDE VELMAC_OVERRIDE
#    define FORT_UPDATETENSOR UPDATETENSOR
#    define FORT_FIX_HOOP_TENSOR FIX_HOOP_TENSOR
#    define FORT_EXTRAPTENSOR EXTRAPTENSOR
#    define FORT_BUILD_MOMENT BUILD_MOMENT
#    define FORT_COPY_VEL_ON_SIGN COPY_VEL_ON_SIGN
#    define FORT_MAKETENSOR MAKETENSOR
#    define FORT_TENSORFORCE TENSORFORCE
#    define FORT_TENSORHEAT TENSORHEAT
#    define FORT_MARANGONIFORCE MARANGONIFORCE
#    define FORT_HEATSOURCE HEATSOURCE
#    define FORT_SEMDELTAFORCE SEMDELTAFORCE
#    define FORT_SEMDELTAFORCE_FACE SEMDELTAFORCE_FACE
#    define FORT_UPDATESEMFORCE UPDATESEMFORCE
#    define FORT_UPDATESEMFORCE_FACE UPDATESEMFORCE_FACE
#    define FORT_SDC_TIME_QUAD SDC_TIME_QUAD
#    define FORT_SDC_TIME_QUAD_FACE SDC_TIME_QUAD_FACE
#    define FORT_BUILD_MASKSEM BUILD_MASKSEM
#    define FORT_BUILD_SLOPES BUILD_SLOPES
#    define FORT_BUILD_SLOPES_FACE BUILD_SLOPES_FACE
#    define FORT_BUILD_CONSERVE BUILD_CONSERVE
#    define FORT_BUILD_MACVOF BUILD_MACVOF
#    define FORT_BUILD_NEWMAC BUILD_NEWMAC
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_SOD_SANITY   sod_sanity
#    define FORT_AGGRESSIVE   aggressive
#    define FORT_DENCOR       dencor
#    define FORT_STEFANSOLVER stefansolver
#    define FORT_HEATSOURCE_FACE heatsource_face
#    define FORT_WALLFUNCTION wallfunction
#    define FORT_INITJUMPTERM initjumpterm
#    define FORT_INIT_ICEMASK init_icemask
#    define FORT_RENORM_MDOT renorm_mdot
#    define FORT_TAGEXPANSION tagexpansion
#    define FORT_DISTRIBUTEEXPANSION distributeexpansion
#    define FORT_CLEAREXPANSION clearexpansion
#    define FORT_VFRAC_SPLIT  vfrac_split
#    define FORT_VFRAC_UNSPLIT vfrac_unsplit
#    define FORT_ESTDT        estdt
#    define FORT_CROSSTERM    crossterm
#    define FORT_FACE_GRADIENTS face_gradients
#    define FORT_COMBINEVEL   combinevel
#    define FORT_COMBINEVELFACE combinevelface
#    define FORT_VISCTENSORHEAT visctensorheat
#    define FORT_VELADVANCE  veladvance
#    define FORT_VELMAC_OVERRIDE velmac_override
#    define FORT_UPDATETENSOR updatetensor
#    define FORT_FIX_HOOP_TENSOR fix_hoop_tensor
#    define FORT_EXTRAPTENSOR extraptensor
#    define FORT_BUILD_MOMENT build_moment
#    define FORT_COPY_VEL_ON_SIGN copy_vel_on_sign
#    define FORT_MAKETENSOR maketensor
#    define FORT_TENSORFORCE tensorforce
#    define FORT_TENSORHEAT tensorheat
#    define FORT_MARANGONIFORCE marangoniforce
#    define FORT_HEATSOURCE heatsource
#    define FORT_SEMDELTAFORCE semdeltaforce
#    define FORT_SEMDELTAFORCE_FACE semdeltaforce_face
#    define FORT_UPDATESEMFORCE updatesemforce
#    define FORT_UPDATESEMFORCE_FACE updatesemforce_face
#    define FORT_SDC_TIME_QUAD sdc_time_quad
#    define FORT_SDC_TIME_QUAD_FACE sdc_time_quad_face
#    define FORT_BUILD_MASKSEM build_masksem
#    define FORT_BUILD_SLOPES build_slopes
#    define FORT_BUILD_SLOPES_FACE build_slopes_face
#    define FORT_BUILD_CONSERVE build_conserve
#    define FORT_BUILD_MACVOF build_macvof
#    define FORT_BUILD_NEWMAC build_newmac
#  else
#    define FORT_SOD_SANITY   sod_sanity_
#    define FORT_AGGRESSIVE   aggressive_
#    define FORT_DENCOR       dencor_
#    define FORT_STEFANSOLVER stefansolver_
#    define FORT_HEATSOURCE_FACE heatsource_face_
#    define FORT_WALLFUNCTION wallfunction_
#    define FORT_INITJUMPTERM initjumpterm_
#    define FORT_INIT_ICEMASK init_icemask_
#    define FORT_RENORM_MDOT renorm_mdot_
#    define FORT_TAGEXPANSION tagexpansion_
#    define FORT_DISTRIBUTEEXPANSION distributeexpansion_
#    define FORT_CLEAREXPANSION clearexpansion_
#    define FORT_VFRAC_SPLIT  vfrac_split_
#    define FORT_VFRAC_UNSPLIT vfrac_unsplit_
#    define FORT_ESTDT        estdt_
#    define FORT_CROSSTERM    crossterm_
#    define FORT_FACE_GRADIENTS face_gradients_
#    define FORT_COMBINEVEL   combinevel_
#    define FORT_COMBINEVELFACE combinevelface_
#    define FORT_VISCTENSORHEAT visctensorheat_
#    define FORT_VELADVANCE  veladvance_
#    define FORT_VELMAC_OVERRIDE velmac_override_
#    define FORT_UPDATETENSOR updatetensor_
#    define FORT_FIX_HOOP_TENSOR fix_hoop_tensor_
#    define FORT_EXTRAPTENSOR extraptensor_
#    define FORT_BUILD_MOMENT build_moment_
#    define FORT_COPY_VEL_ON_SIGN copy_vel_on_sign_
#    define FORT_MAKETENSOR maketensor_
#    define FORT_TENSORFORCE tensorforce_
#    define FORT_TENSORHEAT tensorheat_
#    define FORT_MARANGONIFORCE marangoniforce_
#    define FORT_HEATSOURCE heatsource_
#    define FORT_SEMDELTAFORCE semdeltaforce_
#    define FORT_SEMDELTAFORCE_FACE semdeltaforce_face_
#    define FORT_UPDATESEMFORCE updatesemforce_
#    define FORT_UPDATESEMFORCE_FACE updatesemforce_face_
#    define FORT_SDC_TIME_QUAD sdc_time_quad_
#    define FORT_SDC_TIME_QUAD_FACE sdc_time_quad_face_
#    define FORT_BUILD_MASKSEM build_masksem_
#    define FORT_BUILD_SLOPES build_slopes_
#    define FORT_BUILD_SLOPES_FACE build_slopes_face_
#    define FORT_BUILD_CONSERVE build_conserve_
#    define FORT_BUILD_MACVOF build_macvof_
#    define FORT_BUILD_NEWMAC build_newmac_
#  endif
#  endif

#include <ArrayLim.H>

extern "C"
{
   void FORT_SOD_SANITY(
    int* id,int* nc,
    const int* lo,const int* hi,
    Real* snew,ARLIM_P(snewlo),ARLIM_P(snewhi));

   void FORT_AGGRESSIVE(
    int* datatype,
    Real* warning_cutoff,
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,
    const int* growlo,const int* growhi,
    int* bfact,
    const Real* dx,
    int* scomp,
    int* ncomp,
    int* ndefined,
    int* ngrow,int* dir,int* id,
    int* verbose,
    int* force_check,
    const int* gridno,
    int* ngrid,int* level,int* finest_level,
    Real* mf,ARLIM_P(mflo),ARLIM_P(mfhi));

   void FORT_TENSORFORCE(
    int* elasticface_flag,
    int* massface_index,
    int* vofface_index,
    int* ncphys,
    int* nstate,
    const Real* xlo,const Real* dx,
    Real* xface,ARLIM_P(xfacelo),ARLIM_P(xfacehi),
    Real* yface,ARLIM_P(yfacelo),ARLIM_P(yfacehi),
    Real* zface,ARLIM_P(zfacelo),ARLIM_P(zfacehi),
    Real* lsfab,ARLIM_P(lsfablo),ARLIM_P(lsfabhi),
    Real* rhoinverse,
    ARLIM_P(rhoinverselo),ARLIM_P(rhoinversehi),
    Real* velnew,ARLIM_P(velnewlo),ARLIM_P(velnewhi),
    Real* tenfab,
    ARLIM_P(tenfablo),ARLIM_P(tenfabhi),
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,int* bfact,
    int* level,
    Real* dt,
    int* irz,
    int* im,int* nmat,int* nden);

   void FORT_TENSORHEAT(
    int* elasticface_flag,
    int* massface_index,
    int* vofface_index,
    int* ncphys,
    int* ntensor,
    int* ntensorMM,
    int* nstate,
    const Real* xlo,const Real* dx,
    Real* xface,ARLIM_P(xfacelo),ARLIM_P(xfacehi),
    Real* yface,ARLIM_P(yfacelo),ARLIM_P(yfacehi),
    Real* zface,ARLIM_P(zfacelo),ARLIM_P(zfacehi),
    Real* lsfab,ARLIM_P(lsfablo),ARLIM_P(lsfabhi),
    Real* DeDTinverse,
    ARLIM_P(DeDTinverselo),ARLIM_P(DeDTinversehi),
    Real* vischeat,ARLIM_P(vischeatlo),ARLIM_P(vischeathi),
    Real* tenfab,ARLIM_P(tenfablo),ARLIM_P(tenfabhi),
    Real* gradufab,ARLIM_P(gradufablo),ARLIM_P(gradufabhi),
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,int* bfact,
    int* level,
    Real* dt,
    int* irz,
    int* im,int* nmat,int* nden);

   void FORT_MARANGONIFORCE(
    int* conservative_tension_force,
    int* isweep,
    int* nstate,
    int* nten,
    int* num_curv,
    const Real* xlo,const Real* dx,
    int* facecut_index,
    int* icefacecut_index,
    int* curv_index,
    int* pforce_index,
    int* faceden_index,
    int* icemask_index,
    int* massface_index,
    int* vofface_index,
    int* ncphys,
    Real* xface,ARLIM_P(xfacelo),ARLIM_P(xfacehi),
    Real* yface,ARLIM_P(yfacelo),ARLIM_P(yfacehi),
    Real* zface,ARLIM_P(zfacelo),ARLIM_P(zfacehi),
    Real* xp,ARLIM_P(xplo),ARLIM_P(xphi),
    Real* yp,ARLIM_P(yplo),ARLIM_P(yphi),
    Real* zp,ARLIM_P(zplo),ARLIM_P(zphi),
    Real* maskcov,ARLIM_P(maskcovlo),ARLIM_P(maskcovhi),
    Real* masknbr,ARLIM_P(masknbrlo),ARLIM_P(masknbrhi),
    Real* vol,ARLIM_P(vollo),ARLIM_P(volhi),
    Real* areax,ARLIM_P(areaxlo),ARLIM_P(areaxhi),
    Real* areay,ARLIM_P(areaylo),ARLIM_P(areayhi),
    Real* areaz,ARLIM_P(areazlo),ARLIM_P(areazhi),
    Real* xflux,ARLIM_P(xfluxlo),ARLIM_P(xfluxhi),
    Real* yflux,ARLIM_P(yfluxlo),ARLIM_P(yfluxhi),
    Real* zflux,ARLIM_P(zfluxlo),ARLIM_P(zfluxhi),
    Real* vel,ARLIM_P(vello),ARLIM_P(velhi),
    Real* den,ARLIM_P(denlo),ARLIM_P(denhi),
    Real* ls,ARLIM_P(lslo),ARLIM_P(lshi),
    Real* lsho,ARLIM_P(lsholo),ARLIM_P(lshohi),
    Real* rhoinverse,
    ARLIM_P(rhoinverselo),ARLIM_P(rhoinversehi),
    Real* vof,ARLIM_P(voflo),ARLIM_P(vofhi),
    Real* curv,ARLIM_P(curvlo),ARLIM_P(curvhi),
    Real* velnew,ARLIM_P(velnewlo),ARLIM_P(velnewhi),
    Real* umacnew,ARLIM_P(umacnewlo),ARLIM_P(umacnewhi),
    Real* vmacnew,ARLIM_P(vmacnewlo),ARLIM_P(vmacnewhi),
    Real* wmacnew,ARLIM_P(wmacnewlo),ARLIM_P(wmacnewhi),
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,
    int* bfact,
    int* bfact_grid,
    int* level,
    int* finest_level,
    Real* dt,
    Real* cur_time,
    Real* visc_coef,
    int* solvability_projection,
    const int* presbc,
    const int* velbc,
    const int* vofbc,
    int* nmat,
    int* nden);

   void FORT_SEMDELTAFORCE(
    int* nstate,
    int* nfluxSEM,
    int* nstate_SDC,
    int* nmat,
    int* project_option,
    const Real* xlo,const Real* dx,
    Real* deltafab,ARLIM_P(deltafablo),ARLIM_P(deltafabhi),
    Real* maskSEMfab,ARLIM_P(maskSEMfablo),ARLIM_P(maskSEMfabhi),
    Real* rhoinverse,
    ARLIM_P(rhoinverselo),ARLIM_P(rhoinversehi),
    Real* DeDTinverse,
    ARLIM_P(DeDTinverselo),ARLIM_P(DeDTinversehi),
    Real* velnew,ARLIM_P(velnewlo),ARLIM_P(velnewhi),
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,int* bfact,
    int* level,
    Real* dt);


   void FORT_SEMDELTAFORCE_FACE(
    int* dir,
    int* faceden_index,
    int* ncphys,
    const Real* xlo,const Real* dx,
    Real* deltafab,ARLIM_P(deltafablo),ARLIM_P(deltafabhi),
    Real* maskSEMfab,ARLIM_P(maskSEMfablo),ARLIM_P(maskSEMfabhi),
    Real* xface,ARLIM_P(xfacelo),ARLIM_P(xfacehi),
    Real* xmac,ARLIM_P(xmaclo),ARLIM_P(xmachi),
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,int* bfact,
    int* level,
    Real* dt);



   void FORT_UPDATESEMFORCE(
    int* ns_time_order,
    int* slab_step,
    int* nsolve,
    int* update_spectral,
    int* update_stable,
    int* nstate,
    int* nfluxSEM,
    int* nstate_SDC,
    int* nmat,
    int* project_option,
    const Real* xlo,const Real* dx,
    Real* gpfab,ARLIM_P(gpfablo),ARLIM_P(gpfabhi),
    Real* divfab,ARLIM_P(divfablo),ARLIM_P(divfabhi),
    Real* hoopfab,ARLIM_P(hoopfablo),ARLIM_P(hoopfabhi),
    Real* HOfab,ARLIM_P(HOfablo),ARLIM_P(HOfabhi),
    Real* LOfab,ARLIM_P(LOfablo),ARLIM_P(LOfabhi),
    Real* maskSEMfab,ARLIM_P(maskSEMfablo),ARLIM_P(maskSEMfabhi),
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,int* bfact,
    int* level,
    Real* dt);


   void FORT_UPDATESEMFORCE_FACE(
    int* project_option,
    int* num_materials_face,
    int* nsolveMM_FACE,
    int* ns_time_order,
    int* dir,
    int* slab_step,
    int* update_spectral,
    int* update_stable,
    int* nmat,
    const Real* xlo,const Real* dx,
    Real* gpfab,ARLIM_P(gpfablo),ARLIM_P(gpfabhi),
    Real* HOfab,ARLIM_P(HOfablo),ARLIM_P(HOfabhi),
    Real* LOfab,ARLIM_P(LOfablo),ARLIM_P(LOfabhi),
    Real* maskSEMfab,ARLIM_P(maskSEMfablo),ARLIM_P(maskSEMfabhi),
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,int* bfact,
    int* level,
    Real* dt);



   void FORT_SDC_TIME_QUAD(
    int* HOncomp,int* LOncomp,int* delta_ncomp,
    int* nstate,
    int* nfluxSEM,
    int* nstate_SDC,
    int* nmat,
    const Real* xlo,const Real* dx,
    Real* delta,ARLIM_P(deltalo),ARLIM_P(deltahi),
    Real* HOfab,ARLIM_P(HOfablo),ARLIM_P(HOfabhi),
    Real* LOfab,ARLIM_P(LOfablo),ARLIM_P(LOfabhi),
    Real* maskSEM,ARLIM_P(maskSEMlo),ARLIM_P(maskSEMhi),
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,int* bfact,
    int* level,
    int* finest_level,
    Real* dt);


   void FORT_SDC_TIME_QUAD_FACE(
    int* dir,
    int* HOncomp,
    int* LOncomp,
    int* delta_ncomp,
    int* nstate,
    int* nfluxSEM,
    int* nstate_SDC,
    int* nmat,
    const Real* xlo,const Real* dx,
    Real* delta,ARLIM_P(deltalo),ARLIM_P(deltahi),
    Real* HOfab,ARLIM_P(HOfablo),ARLIM_P(HOfabhi),
    Real* LOfab,ARLIM_P(LOfablo),ARLIM_P(LOfabhi),
    Real* maskSEM,ARLIM_P(maskSEMlo),ARLIM_P(maskSEMhi),
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,int* bfact,
    int* level,
    int* finest_level,
    Real* dt);


   void FORT_BUILD_MASKSEM(
    Real* spectral_cells_level,
    int* mask_sweep,
    int* level,
    int* finest_level,
    int* enable_spectral,
    const int* domlo,const int* domhi,
    int* vofbc,
    Real* maskcov,ARLIM_P(maskcovlo),ARLIM_P(maskcovhi),
    Real* masknbr,ARLIM_P(masknbrlo),ARLIM_P(masknbrhi),
    Real* mask,ARLIM_P(masklo),ARLIM_P(maskhi),
    Real* oldmask,ARLIM_P(oldmasklo),ARLIM_P(oldmaskhi),
    Real* vfrac,ARLIM_P(vfraclo),ARLIM_P(vfrachi),
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,
    int* bfact,
    int* bfact_fine,
    int* nmat);

   void FORT_BUILD_SLOPES(
    Real* masknbr,ARLIM_P(masknbrlo),ARLIM_P(masknbrhi),
    Real* recon,ARLIM_P(reconlo),ARLIM_P(reconhi),
    Real* data,ARLIM_P(datalo),ARLIM_P(datahi),
    Real* slopes,ARLIM_P(slopeslo),ARLIM_P(slopeshi),
    int* nc_conserve,
    int* nmat,
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,int* bfact,
    int* level,
    int* finest_level,
    int* velbc,
    const Real* xlo,const Real* dx,
    int* normdir,
    int* ngrow,
    int* advection_order,
    int* density_advection_order,
    int* slope_limiter_option);


   void FORT_BUILD_SLOPES_FACE(
    Real* masknbr,ARLIM_P(masknbrlo),ARLIM_P(masknbrhi),
    Real* vof,ARLIM_P(voflo),ARLIM_P(vofhi),
    Real* data,ARLIM_P(datalo),ARLIM_P(datahi),
    Real* slopes,ARLIM_P(slopeslo),ARLIM_P(slopeshi),
    int* nmat,
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,int* bfact,
    int* level,
    int* finest_level,
    int* velbc,
    const Real* xlo,const Real* dx,
    int* normdir,
    int* slopedir,
    int* ngrow,
    int* advection_order,
    int* slope_limiter_type);


   void FORT_BUILD_CONSERVE(
    int* iden_base,
    int* override_density,
    int* temperature_primitive_variable,
    Real* conserve,ARLIM_P(conservelo),ARLIM_P(conservehi),
    Real* den,ARLIM_P(denlo),ARLIM_P(denhi),
    Real* vel,ARLIM_P(vello),ARLIM_P(velhi),
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,int* bfact,
    int* nmat,int* ngrow,
    int* normdir,
    int* nc_conserve,
    int* nc_den);


   void FORT_BUILD_MACVOF(
    int* unsplit_advection,
    int* nsolveMM_FACE,
    int* level,
    int* finest_level,
    int* normdir,
    int* nrefine_vof,
    int* nrefine_cen,
    Real* vofF,ARLIM_P(vofFlo),ARLIM_P(vofFhi),
    Real* cenF,ARLIM_P(cenFlo),ARLIM_P(cenFhi),
    Real* xold,ARLIM_P(xoldlo),ARLIM_P(xoldhi),
    Real* xvof,ARLIM_P(xvoflo),ARLIM_P(xvofhi),
    Real* xvel,ARLIM_P(xvello),ARLIM_P(xvelhi), // xvelleft,xvelright
    Real* xvelslope,ARLIM_P(xvelslopelo),ARLIM_P(xvelslopehi), //xvelslope,xcen
    const Real* xlo,const Real* dx,
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,
    int* bfact,
    int* nmat,int* ngrow,
    int* ngrowmac,int* veldir);


   void FORT_BUILD_NEWMAC(
    int* nsolveMM_FACE,
    int* normdir,
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,
    int* bfact,
    Real* xmomside,ARLIM_P(xmomsidelo),ARLIM_P(xmomsidehi),
    Real* ymomside,ARLIM_P(ymomsidelo),ARLIM_P(ymomsidehi),
    Real* zmomside,ARLIM_P(zmomsidelo),ARLIM_P(zmomsidehi),
    Real* xmassside,ARLIM_P(xmasssidelo),ARLIM_P(xmasssidehi),
    Real* ymassside,ARLIM_P(ymasssidelo),ARLIM_P(ymasssidehi),
    Real* zmassside,ARLIM_P(zmasssidelo),ARLIM_P(zmasssidehi),
    Real* xmac,ARLIM_P(xmaclo),ARLIM_P(xmachi),
    Real* ymac,ARLIM_P(ymaclo),ARLIM_P(ymachi),
    Real* zmac,ARLIM_P(zmaclo),ARLIM_P(zmachi),
    Real* mask,ARLIM_P(masklo),ARLIM_P(maskhi),
    Real* maskunsplit,
    ARLIM_P(maskunsplitlo),ARLIM_P(maskunsplithi),
    const Real* xlo,const Real* dx,
    int* nmat,
    int* level,
    int* finest_level);

 
   void FORT_MAKETENSOR(
    int* ncomp_visc,int* im,
    const Real* xlo,const Real* dx,
    Real* visc,ARLIM_P(visclo),ARLIM_P(vischi),
    Real* tensor,ARLIM_P(tlo),ARLIM_P(thi),
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,int* bfact,
    Real* elastic_viscosity,
    Real* etaS,
    Real* elastic_time,
    int* viscoelastic_model,
    Real* polymer_factor,
    int* irz,int* ngrow,int* nmat);


   void FORT_COPY_VEL_ON_SIGN(
    int* im_part,
    int* nparts,
    int* partid,
    int* ngrowFSI,
    int* nFSI,
    int* nFSI_sub,
    const Real* xlo,const Real* dx,
    Real* snew,ARLIM_P(snewlo),ARLIM_P(snewhi),
    Real* fsi,ARLIM_P(fsilo),ARLIM_P(fsihi),
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,int* bfact,
    int* nmat,int* nstate);

   void FORT_BUILD_MOMENT(
    int* level,
    int* finest_level,
    int* nFSI,
    int* nFSI_sub,
    int* nparts,
    int* ngrowFSI,
    int* im_solid_map,
    const Real* xlo,const Real* dx,
    Real* snew,ARLIM_P(snewlo),ARLIM_P(snewhi),
    Real* lsnew,ARLIM_P(lsnewlo),ARLIM_P(lsnewhi),
    Real* fsi,ARLIM_P(fsilo),ARLIM_P(fsihi),
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,int* bfact,
    int* nmat,int* nstate);

 
   void FORT_UPDATETENSOR(
    int* level,
    int* finest_level,
    int* nmat,int* im,int* ncomp_visc,
    Real* vof,ARLIM_P(voflo),ARLIM_P(vofhi),
    Real* visc,ARLIM_P(visclo),ARLIM_P(vischi),
    Real* tendata,ARLIM_P(tendatalo),ARLIM_P(tendatahi),
    const Real* dx,const Real* xlo,
    Real* vel,ARLIM_P(vello),ARLIM_P(velhi),
    Real* tensor_new,ARLIM_P(tnewlo),ARLIM_P(tnewhi),
    Real* tensor_old,ARLIM_P(toldlo),ARLIM_P(toldhi),
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,
    int* bfact,
    Real* dt,
    Real* elastic_time,
    int* viscoelastic_model,
    Real* polymer_factor,
    int* irz,const int* bcrec,
    int* transposegradu);


   void FORT_FIX_HOOP_TENSOR(
    int* level,
    int* finest_level,
    int* nmat,int* im,
    Real* vof,ARLIM_P(voflo),ARLIM_P(vofhi),
    const Real* dx,const Real* xlo,
    Real* tensor_new,ARLIM_P(tnewlo),ARLIM_P(tnewhi),
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,
    int* bfact,
    int* irz);


   void FORT_EXTRAPTENSOR(
    int* level,
    int* finest_level,
    int* nmat,
    int* im,
    int* ngrow_extrap,
    Real* vof,ARLIM_P(voflo),ARLIM_P(vofhi),
    const Real* dx,const Real* xlo,
    Real* tensor_new,ARLIM_P(tnewlo),ARLIM_P(tnewhi),
    Real* tensor_old,ARLIM_P(toldlo),ARLIM_P(toldhi),
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,int* bfact);

 
   void FORT_ESTDT(
      int* nsolveMM_FACE,
      int* enable_spectral,
      int* microlayer_substrate,
      Real* microlayer_angle,
      Real* microlayer_size,
      Real* macrolayer_size,
      Real* latent_heat,
      Real* reaction_rate,
      int* freezing_model,
      int* distribute_from_target,
      Real* saturation_temp,
      int* mass_fraction_id,
      Real* species_evaporation_density,
      const Real* Umac,  ARLIM_P(umaclo), ARLIM_P(umachi),
      const Real* Ucell,  ARLIM_P(ucelllo), ARLIM_P(ucellhi),
      const Real* sol,  ARLIM_P(sollo), ARLIM_P(solhi),
      Real* den,ARLIM_P(denlo),ARLIM_P(denhi),
      Real* vof,ARLIM_P(voflo),ARLIM_P(vofhi),
      Real* dist,ARLIM_P(distlo),ARLIM_P(disthi),
      const Real* xlo,const Real* dx,
      const int* tilelo, const int* tilehi,
      const int* fablo, const int* fabhi,
      int* bfact,
      Real* cap_wave_speed,
      Real* u_max,
      Real* u_max_estdt,
      Real* dt_min,
      int* rzflag,
      Real* Uref,Real* Lref,
      int* nten,
      int* use_lsa,
      Real* denconst,
      Real* denconst_gravity,
      Real* visc_coef,
      Real* ns_gravity,
      int* terminal_velocity_dt,
      int* dir,
      int* nmat,
      int* nparts,
      int* nparts_def,
      int* im_solid_map,
      int* material_type,
      Real* time,int* shock_timestep,
      Real* cfl,
      int* EILE_flag,
      int* level,
      int* finest_level);


   void FORT_WALLFUNCTION(
    int* im_solid_map,
    int* level,
    int* finest_level,
    int* ngrow_law_of_wall,
    int* ngrow_distance,
    int* nmat,int* nparts,int* nden,
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,
    int* bfact,
    const Real* xlo,const Real* dx,
    Real* dt,
    Real* LS,ARLIM_P(LSlo),ARLIM_P(LShi), 
    Real* state,ARLIM_P(statelo),ARLIM_P(statehi), 
    Real* ufluid,ARLIM_P(ufluidlo),ARLIM_P(ufluidhi),
    Real* usolid,ARLIM_P(usolidlo),ARLIM_P(usolidhi),
    Real* ughost,ARLIM_P(ughostlo),ARLIM_P(ughosthi));


   void FORT_STEFANSOLVER(
    int* solidheat_flag,
    Real* microlayer_size,
    int* microlayer_substrate,
    Real* microlayer_temperature_substrate,
    int* adjust_temperature,
    int* nmat,int* nten,int* nstate,
    Real* latent_heat,
    int* freezing_model,
    int* distribute_from_target,
    Real* saturation_temp,
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,
    int* bfact,
    int* level,
    int* finest_level,
    int* nfacefrac,
    int* ncellfrac,
    const Real* xlo,const Real* dx,
    Real* dt,
    Real* cellmm,ARLIM_P(cellmmlo),ARLIM_P(cellmmhi), 
    Real* xfacemm,ARLIM_P(xfacemmlo),ARLIM_P(xfacemmhi), 
    Real* yfacemm,ARLIM_P(yfacemmlo),ARLIM_P(yfacemmhi), 
    Real* zfacemm,ARLIM_P(zfacemmlo),ARLIM_P(zfacemmhi), 
    Real* swept,ARLIM_P(sweptlo),ARLIM_P(swepthi), 
    Real* LS,ARLIM_P(LSlo),ARLIM_P(LShi), 
    Real* thermal,ARLIM_P(thermallo),ARLIM_P(thermalhi), 
    Real* Snew,ARLIM_P(Snewlo),ARLIM_P(Snewhi), 
    Real* DeDT,ARLIM_P(DeDTlo),ARLIM_P(DeDThi),
    Real* den,ARLIM_P(denlo),ARLIM_P(denhi),
    Real* coeff,ARLIM_P(coefflo),ARLIM_P(coeffhi),
    Real* vol,ARLIM_P(vollo),ARLIM_P(volhi),
    Real* heatx,ARLIM_P(heatxlo),ARLIM_P(heatxhi),
    Real* heaty,ARLIM_P(heatylo),ARLIM_P(heatyhi),
    Real* heatz,ARLIM_P(heatzlo),ARLIM_P(heatzhi),
    Real* areax,ARLIM_P(areaxlo),ARLIM_P(areaxhi),
    Real* areay,ARLIM_P(areaylo),ARLIM_P(areayhi),
    Real* areaz,ARLIM_P(areazlo),ARLIM_P(areazhi) ); 


   void FORT_HEATSOURCE_FACE(
    int* nmat,int* nten,int* nstate,
    Real* latent_heat,
    Real* saturation_temp,
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,
    int* bfact,
    const Real* xlo,const Real* dx,
    Real* dt,
    Real* LS,ARLIM_P(LSlo),ARLIM_P(LShi), 
    Real* Snew,ARLIM_P(Snewlo),ARLIM_P(Snewhi), 
    Real* DeDT,ARLIM_P(DeDTlo),ARLIM_P(DeDThi),
    Real* den,ARLIM_P(denlo),ARLIM_P(denhi),
    Real* vol,ARLIM_P(vollo),ARLIM_P(volhi),
    Real* heatx,ARLIM_P(heatxlo),ARLIM_P(heatxhi),
    Real* heaty,ARLIM_P(heatylo),ARLIM_P(heatyhi),
    Real* heatz,ARLIM_P(heatzlo),ARLIM_P(heatzhi),
    Real* areax,ARLIM_P(areaxlo),ARLIM_P(areaxhi),
    Real* areay,ARLIM_P(areaylo),ARLIM_P(areayhi),
    Real* areaz,ARLIM_P(areazlo),ARLIM_P(areazhi) ); 


   void FORT_HEATSOURCE(
    int* nstate,
    int* nmat,
    int* nden,
    const Real* xlo,const Real* dx,
    Real* temperature_source,
    Real* temperature_source_cen,
    Real* temperature_source_rad,
    Real* rhoinverse,
    ARLIM_P(rhoinverselo),ARLIM_P(rhoinversehi), 
    Real* DeDTinverse,
    ARLIM_P(DeDTinverselo),ARLIM_P(DeDTinversehi), 
    Real* Tnew,ARLIM_P(Tnewlo),ARLIM_P(Tnewhi), 
    Real* lsfab,ARLIM_P(lsfablo),ARLIM_P(lsfabhi), 
    Real* recon,ARLIM_P(reconlo),ARLIM_P(reconhi), 
    Real* vol,ARLIM_P(vollo),ARLIM_P(volhi),
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,
    int* bfact,
    int* level,
    int* finest_level,
    Real* dt,Real* time);

   void FORT_INITJUMPTERM(
    Real* mdotplus,
    Real* mdotminus,
    Real* mdotcount,
    int* ngrow_expansion,
    Real* time,
    int* level,int* finest_level,
    int* nmat,int* nten,
    Real* latent_heat,
    Real* saturation_temp,
    int* freezing_model,
    int* distribute_from_target,
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,
    int* bfact,
    const Real* xlo,const Real* dx,Real* dt,
    Real* maskcov,ARLIM_P(maskcovlo),ARLIM_P(maskcovhi),
    Real* jump,ARLIM_P(jumplo),ARLIM_P(jumphi),
    Real* mdot,ARLIM_P(mdotlo),ARLIM_P(mdothi),
    Real* LSnew,ARLIM_P(LSnewlo),ARLIM_P(LSnewhi));


   void FORT_INIT_ICEMASK(
    Real* time,
    int* facecut_index,
    int* icefacecut_index,
    int* icemask_index,
    int* massface_index,
    int* vofface_index,
    int* ncphys,
    int* level,int* finest_level,
    int* nmat,int* nten,
    Real* latent_heat,
    Real* saturation_temp,
    int* freezing_model,
    int* distribute_from_target,
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,
    int* bfact,
    const Real* xlo,const Real* dx,Real* dt,
    Real* maskcov,ARLIM_P(maskcovlo),ARLIM_P(maskcovhi),
    Real* xface,ARLIM_P(xfacelo),ARLIM_P(xfacehi),
    Real* yface,ARLIM_P(yfacelo),ARLIM_P(yfacehi),
    Real* zface,ARLIM_P(zfacelo),ARLIM_P(zfacehi),
    Real* LSnew,ARLIM_P(LSnewlo),ARLIM_P(LSnewhi),
    Real* recon,ARLIM_P(reconlo),ARLIM_P(reconhi));


   void FORT_RENORM_MDOT(
    Real* mdotplus,Real* mdotminus,Real* mdotcount,
    int* level,int* finest_level,
    int* nmat,int* nten,
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,
    int* bfact,
    const Real* xlo,const Real* dx,Real* dt,
    Real* maskcov,ARLIM_P(maskcovlo),ARLIM_P(maskcovhi),
    Real* mdot,ARLIM_P(mdotlo),ARLIM_P(mdothi));


   void FORT_TAGEXPANSION(
    Real* latent_heat,
    int* freezing_model,
    int* distribute_from_target,
    int* ngrow_expansion,
    Real* time,
    int* vofbc,
    Real* expect_mdot_sign,
    Real* mdot_sum,
    int* im_source,
    int* im_dest,
    int* indexEXP,
    int* level,int* finest_level,
    int* nmat,int* nten,
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,
    int* bfact,
    const Real* xlo,const Real* dx,Real* dt,
    Real* maskcov,ARLIM_P(maskcovlo),ARLIM_P(maskcovhi),
    Real* tag,ARLIM_P(tablo),ARLIM_P(taghi),
    Real* jump,ARLIM_P(jumplo),ARLIM_P(jumphi),
    Real* LSnew,ARLIM_P(LSnewlo),ARLIM_P(LSnewhi),
    Real* recon,ARLIM_P(reconlo),ARLIM_P(reconhi));


   void FORT_DISTRIBUTEEXPANSION(
    int* ngrow_expansion,
    int* im_source,
    int* im_dest,
    int* indexEXP,
    int* level,int* finest_level,
    int* nmat,int* nten,
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,
    int* bfact,
    const Real* xlo,const Real* dx,Real* dt,
    Real* maskcov,ARLIM_P(maskcovlo),ARLIM_P(maskcovhi),
    Real* LS,ARLIM_P(LSlo),ARLIM_P(LShi),
    Real* tag,ARLIM_P(tablo),ARLIM_P(taghi),
    Real* jump,ARLIM_P(jumplo),ARLIM_P(jumphi));

   void FORT_CLEAREXPANSION(
    int* ngrow_expansion,
    Real* mdot_sum,
    Real* mdot_lost,
    int* im_source,
    int* im_dest,
    int* indexEXP,
    int* level,int* finest_level,
    int* nmat,int* nten,
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,
    int* bfact,
    const Real* xlo,const Real* dx,Real* dt,
    Real* maskcov,ARLIM_P(maskcovlo),ARLIM_P(maskcovhi),
    Real* tag,ARLIM_P(taglo),ARLIM_P(taghi),
    Real* jump,ARLIM_P(jumplo),ARLIM_P(jumphi));


   void FORT_VFRAC_SPLIT(
     int* nsolveMM_FACE,
     int* nprocessed,
     int* tid,
     int* make_interface_incomp,
     Real* added_weight,
     Real* density_floor,
     Real* density_ceiling,
     int* solidheat_flag,
     int* temperature_primitive_variable,
     int* dencomp,int* mofcomp,int* errcomp,
     Real* latent_heat,
     int* freezing_model,
     int* distribute_from_target,
     int* nten,
     int* face_flag,
     int* override_density,
     int* velbc,
     int* EILE_flag,
     int* dir_counter,
     int* normdir,
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,
     int* bfact,
     int* bfact_f,
     Real* dt,
     Real* time,
     Real* passive_veltime,
     Real* vofflux,ARLIM_P(voffluxlo),ARLIM_P(voffluxhi),
      // original data
     Real* LS,ARLIM_P(LSlo),ARLIM_P(LShi),
     Real* den,ARLIM_P(denlo),ARLIM_P(denhi),
     Real* tensor,ARLIM_P(tensorlo),ARLIM_P(tensorhi),
     Real* velfab,ARLIM_P(velfablo),ARLIM_P(velfabhi),
      // slope data
     Real* vofslope,ARLIM_P(vofslopelo),ARLIM_P(vofslopehi), 
      // this is the result
     Real* snew,ARLIM_P(snewlo),ARLIM_P(snewhi), 
     Real* tennew,ARLIM_P(tennewlo),ARLIM_P(tennewhi), 
     Real* LSnew,ARLIM_P(LSnewlo),ARLIM_P(LSnewhi), 
      // other vars.
     Real* ucell,ARLIM_P(ucelllo),ARLIM_P(ucellhi),
     Real* vofls0,ARLIM_P(vofls0lo),ARLIM_P(vofls0hi),
     Real* conscor,ARLIM_P(conscorlo),ARLIM_P(conscorhi),
     Real* mask,ARLIM_P(masklo),ARLIM_P(maskhi),
     Real* masknbr,ARLIM_P(masknbrlo),ARLIM_P(masknbrhi),
     Real* unode,ARLIM_P(unodelo),ARLIM_P(unodehi),
     const Real* xlo,const Real* dx,
      // local variables
     Real* conserve,ARLIM_P(conservelo),ARLIM_P(conservehi),
     Real* xvel,ARLIM_P(xvello),ARLIM_P(xvelhi),
     Real* yvel,ARLIM_P(yvello),ARLIM_P(yvelhi),
     Real* zvel,ARLIM_P(zvello),ARLIM_P(zvelhi),
     Real* xvelslope,ARLIM_P(xvelslopelo),ARLIM_P(xvelslopehi),
     Real* yvelslope,ARLIM_P(yvelslopelo),ARLIM_P(yvelslopehi),
     Real* zvelslope,ARLIM_P(zvelslopelo),ARLIM_P(zvelslopehi),
     Real* momslope,ARLIM_P(momslopelo),ARLIM_P(momslopehi),
     Real* xmomside,ARLIM_P(xmomsidelo),ARLIM_P(xmomsidehi),
     Real* ymomside,ARLIM_P(ymomsidelo),ARLIM_P(ymomsidehi),
     Real* zmomside,ARLIM_P(zmomsidelo),ARLIM_P(zmomsidehi),
     Real* xmassside,ARLIM_P(xmasssidelo),ARLIM_P(xmasssidehi),
     Real* ymassside,ARLIM_P(ymasssidelo),ARLIM_P(ymasssidehi),
     Real* zmassside,ARLIM_P(zmasssidelo),ARLIM_P(zmasssidehi),
     int* ngrow,
     int* ngrow_mac_old,
     int* nc_conserve,
     int* iden_base,
     int* nmat,
     int* map_forward,
     int* recon_ncomp,
     int* den_recon_ncomp,
     int* ncomp_state,
     int* ntensor,
     int* nc_bucket,
     int* nrefine_vof,
     int* verbose,
     const int* gridno,int* ngrid,
     int* level,
     int* finest_level,
     int* dombc,
     const int* domlo,const int* domhi);

   void FORT_VFRAC_UNSPLIT(
     int* unsplit_flag,
     int* nsolveMM_FACE,
     int* nprocessed,
     int* tid,
     int* make_interface_incomp,
     Real* added_weight,
     Real* density_floor,
     Real* density_ceiling,
     int* solidheat_flag,
     int* temperature_primitive_variable,
     int* dencomp,int* mofcomp,int* errcomp,
     Real* latent_heat,
     int* freezing_model,
     int* distribute_from_target,
     int* nten,
     int* face_flag,
     int* override_density,
     int* velbc,
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,
     int* bfact,
     int* bfact_f,
     Real* dt,
     Real* time,
     Real* passive_veltime,
      // original data
     Real* LS,ARLIM_P(LSlo),ARLIM_P(LShi),
     Real* den,ARLIM_P(denlo),ARLIM_P(denhi),
     Real* tensor,ARLIM_P(tensorlo),ARLIM_P(tensorhi),
     Real* velfab,ARLIM_P(velfablo),ARLIM_P(velfabhi),
      // slope data
     Real* vofslope,ARLIM_P(vofslopelo),ARLIM_P(vofslopehi), 
      // this is the result
     Real* snew,ARLIM_P(snewlo),ARLIM_P(snewhi), 
     Real* tennew,ARLIM_P(tennewlo),ARLIM_P(tennewhi), 
     Real* LSnew,ARLIM_P(LSnewlo),ARLIM_P(LSnewhi), 
      // other vars.
     Real* ucell,ARLIM_P(ucelllo),ARLIM_P(ucellhi),
     Real* vofls0,ARLIM_P(vofls0lo),ARLIM_P(vofls0hi),
     Real* mask,ARLIM_P(masklo),ARLIM_P(maskhi),
     Real* masknbr,ARLIM_P(masknbrlo),ARLIM_P(masknbrhi),
     Real* maskunsplit,
     ARLIM_P(maskunsplitlo),ARLIM_P(maskunsplithi),
     Real* unode,ARLIM_P(unodelo),ARLIM_P(unodehi),
     Real* vnode,ARLIM_P(vnodelo),ARLIM_P(vnodehi),
     Real* wnode,ARLIM_P(wnodelo),ARLIM_P(wnodehi),
     const Real* xlo,const Real* dx,
      // local variables
     Real* conserve,ARLIM_P(conservelo),ARLIM_P(conservehi),
     Real* xvel,ARLIM_P(xvello),ARLIM_P(xvelhi),
     Real* yvel,ARLIM_P(yvello),ARLIM_P(yvelhi),
     Real* zvel,ARLIM_P(zvello),ARLIM_P(zvelhi),
     Real* xmomside,ARLIM_P(xmomsidelo),ARLIM_P(xmomsidehi),
     Real* ymomside,ARLIM_P(ymomsidelo),ARLIM_P(ymomsidehi),
     Real* zmomside,ARLIM_P(zmomsidelo),ARLIM_P(zmomsidehi),
     Real* xmassside,ARLIM_P(xmasssidelo),ARLIM_P(xmasssidehi),
     Real* ymassside,ARLIM_P(ymasssidelo),ARLIM_P(ymasssidehi),
     Real* zmassside,ARLIM_P(zmasssidelo),ARLIM_P(zmasssidehi),
     int* ngrow,
     int* ngrow_mac_old,
     int* nc_conserve,
     int* iden_base,
     int* nmat,
     int* recon_ncomp,
     int* den_recon_ncomp,
     int* ncomp_state,
     int* ntensor,
     int* nc_bucket,
     int* nrefine_vof,
     int* verbose,
     const int* gridno,int* ngrid,
     int* level,
     int* finest_level,
     int* dombc,
     const int* domlo,const int* domhi);

   void FORT_DENCOR(
     int* presbc,
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,
     int* bfact,
     Real* dt,
     Real* mask,ARLIM_P(masklo),ARLIM_P(maskhi),
     Real* masknbr,ARLIM_P(masknbrlo),ARLIM_P(masknbrhi),
     Real* vol,ARLIM_P(vollo),ARLIM_P(volhi),
     Real* dennew,ARLIM_P(dennewlo),ARLIM_P(dennewhi), 
     const Real* xlo,const Real* dx,
     Real* gravity_normalized,
     Real* DrhoDT,int* override_density,
     int* nmat,int* level,int* finest_level);


   void FORT_COMBINEVEL(
     int* tid,
     int* hflag,
     int* num_materials_combine,
     Real* latent_heat,
     int* freezing_model,
     int* distribute_from_target,
     Real* saturation_temp,
     int* hydrate_flag, // scalar
     int* nmat,
     int* nparts,
     int* nparts_def,
     int* im_solid_map,
     int* nten,
     int* nsolve,
     int* nsolveMM,
     int* project_option,
     int* combine_idx,
     int* combine_flag,
     int* nstate_main,
     int* ncomp_cell,
     int* scomp,
     int* ncomp,
     int* scomp_size,
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,
     int* bfact,
     int* level,
     int* finest_level,
     Real* maskcov,ARLIM_P(maskcovlo),ARLIM_P(maskcovhi),
     Real* sol,ARLIM_P(sollo),ARLIM_P(solhi),
     Real* LSNEW,ARLIM_P(LSNEWlo),ARLIM_P(LSNEWhi),
     Real* LS,ARLIM_P(LSlo),ARLIM_P(LShi),
     Real* vof,ARLIM_P(voflo),ARLIM_P(vofhi),
     Real* cellfab,ARLIM_P(cellfablo),ARLIM_P(cellfabhi),
     Real* newcell,ARLIM_P(newcelllo),ARLIM_P(newcellhi),
     Real* state,ARLIM_P(statelo),ARLIM_P(statehi),
     int* velbc,
     int* listbc,
     const Real* xlo,const Real* dx,
     Real* cur_time,
     Real* prescribed_solid_scale);


   void FORT_COMBINEVELFACE(
     int* tid,
     int* num_materials_combine,
     int* hflag,
     int* facecut_index,
     int* icefacecut_index,
     int* vofface_index,
     int* massface_index,
     int* ncphys,
     int* nmat,
     int* nparts,
     int* nparts_def,
     int* im_solid_map,
     int* nten,
     int* nsolve,
     int* nsolveMM,
     int* nsolveMM_FACE,
     int* project_option,
     int* combine_idx,
     int* combine_flag,
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,
     int* bfact,
     int* level,
     int* finest_level,
     int* velbc,
     Real* vof,ARLIM_P(voflo),ARLIM_P(vofhi),
     Real* mac,ARLIM_P(maclo),ARLIM_P(machi),
     Real* xface,ARLIM_P(xfacelo),ARLIM_P(xfacehi),
     Real* LS,ARLIM_P(LSlo),ARLIM_P(LShi),
     Real* sol,ARLIM_P(sollo),ARLIM_P(solhi),
     const Real* xlo,const Real* dx,
     int* dir,
     Real* cur_time,
     Real* prescribed_solid_scale);

   void FORT_CROSSTERM(
     int* nsolveMM_FACE,
     int* tileloop,
     int* dir,
     int* operation_flag, // 8
     int* enable_spectral,
     int* spectral_loop,
     int* ncfluxreg,
     Real* semflux,ARLIM_P(semfluxlo),ARLIM_P(semfluxhi),
     Real* mask,ARLIM_P(masklo),ARLIM_P(maskhi),
     Real* maskcov,ARLIM_P(maskcovlo),ARLIM_P(maskcovhi),
     Real* faceLS,ARLIM_P(faceLSlo),ARLIM_P(faceLShi),
     Real* mdata,ARLIM_P(mdatalo),ARLIM_P(mdatahi),
     Real* tdata,ARLIM_P(tdatalo),ARLIM_P(tdatahi),
     Real* c_tdata,ARLIM_P(c_tdatalo),ARLIM_P(c_tdatahi),
     Real* maskSEM,ARLIM_P(maskSEMlo),ARLIM_P(maskSEMhi),
     const Real* xlo,const Real* dx,
     Real* dt,
     Real* vel,ARLIM_P(vello),ARLIM_P(velhi),
     Real* levelpc,
     ARLIM_P(levelpclo),ARLIM_P(levelpchi),
     Real* xflux,ARLIM_P(xfluxlo),ARLIM_P(xfluxhi),
     Real* xface,ARLIM_P(xfacelo),ARLIM_P(xfacehi),
     Real* xfacemm,ARLIM_P(xfacemmlo),ARLIM_P(xfacemmhi), 
     Real* xcellmm,ARLIM_P(xcellmmlo),ARLIM_P(xcellmmhi), 
     Real* recon,ARLIM_P(reconlo),ARLIM_P(reconhi), 
     int* facevisc_index,
     int* vofface_index,
     int* massface_index,
     int* ncphys,
     const int* tilelo,const int* tilehi, 
     const int* fablo,const int* fabhi, 
     int* bfact,
     int* level,
     int* rzflag,
     int* velbc,
     Real* visc_coef, 
     int* nmat,
     int* nden,
     int* nfacefrac,
     int* ncellfrac,
     int* ntensor,
     int* ntensorMM,
     int* use_StewartLay,
     int* constant_viscosity,
     int* homflag);


   void FORT_FACE_GRADIENTS(
     int* ns_time_order,
     int* divu_outer_sweeps,
     int* num_divu_outer_sweeps,
     int* SDC_outer_sweeps,
     int* tileloop,
     int* dir,
     int* slab_step,
     int* itensor_iter,
     Real* time,
     int* temperature_primitive_variable,
     int* face_flag,
     int* enable_spectral,
     int* velbc,
     int* spectral_loop,
     int* ncfluxreg,
     Real* semflux,ARLIM_P(semfluxlo),ARLIM_P(semfluxhi),
     Real* mask0,ARLIM_P(mask0lo),ARLIM_P(mask0hi),
     Real* mask3,ARLIM_P(mask3lo),ARLIM_P(mask3hi),
     Real* maskSEM,ARLIM_P(maskSEMlo),ARLIM_P(maskSEMhi),
     Real* faceLS,ARLIM_P(faceLSlo),ARLIM_P(faceLShi),
     Real* mdata,ARLIM_P(mdatalo),ARLIM_P(mdatahi),
     Real* tdata,ARLIM_P(tdatalo),ARLIM_P(tdatahi),
     Real* c_tdata,ARLIM_P(c_tdatalo),ARLIM_P(c_tdatahi),
     Real* vel,ARLIM_P(vello),ARLIM_P(velhi),
     Real* solid,ARLIM_P(solidlo),ARLIM_P(solidlhi),
     Real* levelpc,
     ARLIM_P(levelpclo),ARLIM_P(levelpchi),
     Real* recon,
     ARLIM_P(reconlo),ARLIM_P(reconhi),
     const Real* xlo,
     const Real* dx,
     int* rzflag,
     const int* tilelo,const int* tilehi, 
     const int* fablo,const int* fabhi, 
     int* bfact,int* bfact_c,int* bfact_f,
     int* level,
     int* finest_level,
     int* nmat,
     int* nparts,
     int* nparts_def,
     int* im_solid_map,
     int* homflag,
     int* ntensor,
     int* ntensorMM,
     int* SEM_upwind,
     int* SEM_advection_algorithm);


   void FORT_VISCTENSORHEAT(
     int* ntensor,
     int* ntensorMM,
     int* nsolve,
     int* nsolveMM,
     int* nsolveMM_FACE,
     int* nstate,
     const Real* xlo,const Real* dx,
     Real* lsfab,ARLIM_P(lsfablo),ARLIM_P(lsfabhi),
     Real* DeDTinverse,ARLIM_P(DeDTinverselo),ARLIM_P(DeDTinversehi),
     Real* heatnew,ARLIM_P(heatnewlo),ARLIM_P(heatnewhi),
     Real* xstress,ARLIM_P(xstresslo),ARLIM_P(xstresshi),
     Real* ystress,ARLIM_P(ystresslo),ARLIM_P(ystresshi),
     Real* zstress,ARLIM_P(zstresslo),ARLIM_P(zstresshi),
     Real* gradu,ARLIM_P(gradulo),ARLIM_P(graduhi),
     const int* tilelo,const int* tilehi, 
     const int* fablo,const int* fabhi, 
     int* bfact,
     int* level,
     const Real* dt,int* irz,
     int* nmat,int* nden);

   void FORT_VELADVANCE(
     int* nmat,
     int* nparts,
     int* nparts_def,
     int* im_solid_map,
     int* nsolve,
     int* nsolveMM,
     int* nstate,
     const Real* xlo,const Real* dx,
     Real* solid,ARLIM_P(solidlo),ARLIM_P(solidhi),
     Real* snew,ARLIM_P(snewlo),ARLIM_P(snewhi),
     Real* lsnew,ARLIM_P(lsnewlo),ARLIM_P(lsnewhi),
     Real* du,ARLIM_P(dulo),ARLIM_P(duhi),
     Real* advect,ARLIM_P(advectlo),ARLIM_P(advecthi),
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,
     int* bfact);


   void FORT_VELMAC_OVERRIDE(
     int* isweep,
     int* unsplit_displacement,
     int* nsolveMM_FACE,
     int* nmat,
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,
     int* bfact,
     int* velbc,
     Real* dt,Real* time,
     Real* passive_veltime,
     Real* vel_time,int* normdir,
     Real* unodetemp,ARLIM_P(unodetemplo),ARLIM_P(unodetemphi),
     Real* unode,ARLIM_P(unodelo),ARLIM_P(unodehi),
     Real* ucell,ARLIM_P(ucelllo),ARLIM_P(ucellhi),
     const Real* xlo,const Real* dx,
     int* mac_grow,
     int* map_forward,
     int* level,
     int* finest_level,
     int* SDC_outer_sweeps,
     int* ns_time_order,
     int* divu_outer_sweeps,
     int* num_divu_outer_sweeps);


}
#endif

#endif 
