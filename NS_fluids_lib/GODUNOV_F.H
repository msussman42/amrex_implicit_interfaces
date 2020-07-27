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

#include <AMReX_ArrayLim.H>
#include <AMReX_BLFort.H>

namespace amrex{

extern "C"
{
   void FORT_SOD_SANITY(
    int* id,int* nc,
    const int* lo,const int* hi,
    Real* snew,AMREX_ARLIM_P(snewlo),AMREX_ARLIM_P(snewhi));

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
    Real* mf,AMREX_ARLIM_P(mflo),AMREX_ARLIM_P(mfhi));

   void FORT_TENSORFORCE(
    int* elasticface_flag,
    int* massface_index,
    int* vofface_index,
    int* ncphys,
    int* nstate,
    const Real* xlo,const Real* dx,
    Real* xface,AMREX_ARLIM_P(xfacelo),AMREX_ARLIM_P(xfacehi),
    Real* yface,AMREX_ARLIM_P(yfacelo),AMREX_ARLIM_P(yfacehi),
    Real* zface,AMREX_ARLIM_P(zfacelo),AMREX_ARLIM_P(zfacehi),
    Real* lsfab,AMREX_ARLIM_P(lsfablo),AMREX_ARLIM_P(lsfabhi),
    Real* rhoinverse,
   AMREX_ARLIM_P(rhoinverselo),AMREX_ARLIM_P(rhoinversehi),
    Real* velnew,AMREX_ARLIM_P(velnewlo),AMREX_ARLIM_P(velnewhi),
    Real* tenfab,
   AMREX_ARLIM_P(tenfablo),AMREX_ARLIM_P(tenfabhi),
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
    Real* xface,AMREX_ARLIM_P(xfacelo),AMREX_ARLIM_P(xfacehi),
    Real* yface,AMREX_ARLIM_P(yfacelo),AMREX_ARLIM_P(yfacehi),
    Real* zface,AMREX_ARLIM_P(zfacelo),AMREX_ARLIM_P(zfacehi),
    Real* lsfab,AMREX_ARLIM_P(lsfablo),AMREX_ARLIM_P(lsfabhi),
    Real* DeDTinverse,
   AMREX_ARLIM_P(DeDTinverselo),AMREX_ARLIM_P(DeDTinversehi),
    Real* vischeat,AMREX_ARLIM_P(vischeatlo),AMREX_ARLIM_P(vischeathi),
    Real* tenfab,AMREX_ARLIM_P(tenfablo),AMREX_ARLIM_P(tenfabhi),
    Real* gradufab,AMREX_ARLIM_P(gradufablo),AMREX_ARLIM_P(gradufabhi),
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
    Real* xface,AMREX_ARLIM_P(xfacelo),AMREX_ARLIM_P(xfacehi),
    Real* yface,AMREX_ARLIM_P(yfacelo),AMREX_ARLIM_P(yfacehi),
    Real* zface,AMREX_ARLIM_P(zfacelo),AMREX_ARLIM_P(zfacehi),
    Real* xp,AMREX_ARLIM_P(xplo),AMREX_ARLIM_P(xphi),
    Real* yp,AMREX_ARLIM_P(yplo),AMREX_ARLIM_P(yphi),
    Real* zp,AMREX_ARLIM_P(zplo),AMREX_ARLIM_P(zphi),
    Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi),
    Real* masknbr,AMREX_ARLIM_P(masknbrlo),AMREX_ARLIM_P(masknbrhi),
    Real* vol,AMREX_ARLIM_P(vollo),AMREX_ARLIM_P(volhi),
    Real* areax,AMREX_ARLIM_P(areaxlo),AMREX_ARLIM_P(areaxhi),
    Real* areay,AMREX_ARLIM_P(areaylo),AMREX_ARLIM_P(areayhi),
    Real* areaz,AMREX_ARLIM_P(areazlo),AMREX_ARLIM_P(areazhi),
    Real* xflux,AMREX_ARLIM_P(xfluxlo),AMREX_ARLIM_P(xfluxhi),
    Real* yflux,AMREX_ARLIM_P(yfluxlo),AMREX_ARLIM_P(yfluxhi),
    Real* zflux,AMREX_ARLIM_P(zfluxlo),AMREX_ARLIM_P(zfluxhi),
    Real* vel,AMREX_ARLIM_P(vello),AMREX_ARLIM_P(velhi),
    Real* den,AMREX_ARLIM_P(denlo),AMREX_ARLIM_P(denhi),
    Real* ls,AMREX_ARLIM_P(lslo),AMREX_ARLIM_P(lshi),
    Real* lsho,AMREX_ARLIM_P(lsholo),AMREX_ARLIM_P(lshohi),
    Real* rhoinverse,
   AMREX_ARLIM_P(rhoinverselo),AMREX_ARLIM_P(rhoinversehi),
    Real* vof,AMREX_ARLIM_P(voflo),AMREX_ARLIM_P(vofhi),
    Real* curv,AMREX_ARLIM_P(curvlo),AMREX_ARLIM_P(curvhi),
    Real* velnew,AMREX_ARLIM_P(velnewlo),AMREX_ARLIM_P(velnewhi),
    Real* umacnew,AMREX_ARLIM_P(umacnewlo),AMREX_ARLIM_P(umacnewhi),
    Real* vmacnew,AMREX_ARLIM_P(vmacnewlo),AMREX_ARLIM_P(vmacnewhi),
    Real* wmacnew,AMREX_ARLIM_P(wmacnewlo),AMREX_ARLIM_P(wmacnewhi),
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
    Real* deltafab,AMREX_ARLIM_P(deltafablo),AMREX_ARLIM_P(deltafabhi),
    Real* maskSEMfab,AMREX_ARLIM_P(maskSEMfablo),AMREX_ARLIM_P(maskSEMfabhi),
    Real* rhoinverse,
   AMREX_ARLIM_P(rhoinverselo),AMREX_ARLIM_P(rhoinversehi),
    Real* DeDTinverse,
   AMREX_ARLIM_P(DeDTinverselo),AMREX_ARLIM_P(DeDTinversehi),
    Real* velnew,AMREX_ARLIM_P(velnewlo),AMREX_ARLIM_P(velnewhi),
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,int* bfact,
    int* level,
    Real* dt);


   void FORT_SEMDELTAFORCE_FACE(
    int* dir,
    int* faceden_index,
    int* ncphys,
    const Real* xlo,const Real* dx,
    Real* deltafab,AMREX_ARLIM_P(deltafablo),AMREX_ARLIM_P(deltafabhi),
    Real* maskSEMfab,AMREX_ARLIM_P(maskSEMfablo),AMREX_ARLIM_P(maskSEMfabhi),
    Real* xface,AMREX_ARLIM_P(xfacelo),AMREX_ARLIM_P(xfacehi),
    Real* xmac,AMREX_ARLIM_P(xmaclo),AMREX_ARLIM_P(xmachi),
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
    Real* gpfab,AMREX_ARLIM_P(gpfablo),AMREX_ARLIM_P(gpfabhi),
    Real* divfab,AMREX_ARLIM_P(divfablo),AMREX_ARLIM_P(divfabhi),
    Real* hoopfab,AMREX_ARLIM_P(hoopfablo),AMREX_ARLIM_P(hoopfabhi),
    Real* HOfab,AMREX_ARLIM_P(HOfablo),AMREX_ARLIM_P(HOfabhi),
    Real* LOfab,AMREX_ARLIM_P(LOfablo),AMREX_ARLIM_P(LOfabhi),
    Real* maskSEMfab,AMREX_ARLIM_P(maskSEMfablo),AMREX_ARLIM_P(maskSEMfabhi),
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
    Real* gpfab,AMREX_ARLIM_P(gpfablo),AMREX_ARLIM_P(gpfabhi),
    Real* HOfab,AMREX_ARLIM_P(HOfablo),AMREX_ARLIM_P(HOfabhi),
    Real* LOfab,AMREX_ARLIM_P(LOfablo),AMREX_ARLIM_P(LOfabhi),
    Real* maskSEMfab,AMREX_ARLIM_P(maskSEMfablo),AMREX_ARLIM_P(maskSEMfabhi),
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
    Real* delta,AMREX_ARLIM_P(deltalo),AMREX_ARLIM_P(deltahi),
    Real* HOfab,AMREX_ARLIM_P(HOfablo),AMREX_ARLIM_P(HOfabhi),
    Real* LOfab,AMREX_ARLIM_P(LOfablo),AMREX_ARLIM_P(LOfabhi),
    Real* maskSEM,AMREX_ARLIM_P(maskSEMlo),AMREX_ARLIM_P(maskSEMhi),
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
    Real* delta,AMREX_ARLIM_P(deltalo),AMREX_ARLIM_P(deltahi),
    Real* HOfab,AMREX_ARLIM_P(HOfablo),AMREX_ARLIM_P(HOfabhi),
    Real* LOfab,AMREX_ARLIM_P(LOfablo),AMREX_ARLIM_P(LOfabhi),
    Real* maskSEM,AMREX_ARLIM_P(maskSEMlo),AMREX_ARLIM_P(maskSEMhi),
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
    Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi),
    Real* masknbr,AMREX_ARLIM_P(masknbrlo),AMREX_ARLIM_P(masknbrhi),
    Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
    Real* oldmask,AMREX_ARLIM_P(oldmasklo),AMREX_ARLIM_P(oldmaskhi),
    Real* vfrac,AMREX_ARLIM_P(vfraclo),AMREX_ARLIM_P(vfrachi),
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,
    int* bfact,
    int* bfact_fine,
    int* nmat);

   void FORT_BUILD_SLOPES(
    Real* masknbr,AMREX_ARLIM_P(masknbrlo),AMREX_ARLIM_P(masknbrhi),
    Real* recon,AMREX_ARLIM_P(reconlo),AMREX_ARLIM_P(reconhi),
    Real* data,AMREX_ARLIM_P(datalo),AMREX_ARLIM_P(datahi),
    Real* slopes,AMREX_ARLIM_P(slopeslo),AMREX_ARLIM_P(slopeshi),
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
    Real* masknbr,AMREX_ARLIM_P(masknbrlo),AMREX_ARLIM_P(masknbrhi),
    Real* vof,AMREX_ARLIM_P(voflo),AMREX_ARLIM_P(vofhi),
    Real* data,AMREX_ARLIM_P(datalo),AMREX_ARLIM_P(datahi),
    Real* slopes,AMREX_ARLIM_P(slopeslo),AMREX_ARLIM_P(slopeshi),
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
    Real* conserve,AMREX_ARLIM_P(conservelo),AMREX_ARLIM_P(conservehi),
    Real* den,AMREX_ARLIM_P(denlo),AMREX_ARLIM_P(denhi),
    Real* vel,AMREX_ARLIM_P(vello),AMREX_ARLIM_P(velhi),
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
    Real* vofF,AMREX_ARLIM_P(vofFlo),AMREX_ARLIM_P(vofFhi),
    Real* cenF,AMREX_ARLIM_P(cenFlo),AMREX_ARLIM_P(cenFhi),
    Real* xold,AMREX_ARLIM_P(xoldlo),AMREX_ARLIM_P(xoldhi),
    Real* xvof,AMREX_ARLIM_P(xvoflo),AMREX_ARLIM_P(xvofhi),
    Real* xvel,AMREX_ARLIM_P(xvello),AMREX_ARLIM_P(xvelhi), // xvelleft,xvelright
    Real* xvelslope,AMREX_ARLIM_P(xvelslopelo),AMREX_ARLIM_P(xvelslopehi), //xvelslope,xcen
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
    Real* xmomside,AMREX_ARLIM_P(xmomsidelo),AMREX_ARLIM_P(xmomsidehi),
    Real* ymomside,AMREX_ARLIM_P(ymomsidelo),AMREX_ARLIM_P(ymomsidehi),
    Real* zmomside,AMREX_ARLIM_P(zmomsidelo),AMREX_ARLIM_P(zmomsidehi),
    Real* xmassside,AMREX_ARLIM_P(xmasssidelo),AMREX_ARLIM_P(xmasssidehi),
    Real* ymassside,AMREX_ARLIM_P(ymasssidelo),AMREX_ARLIM_P(ymasssidehi),
    Real* zmassside,AMREX_ARLIM_P(zmasssidelo),AMREX_ARLIM_P(zmasssidehi),
    Real* xmac,AMREX_ARLIM_P(xmaclo),AMREX_ARLIM_P(xmachi),
    Real* ymac,AMREX_ARLIM_P(ymaclo),AMREX_ARLIM_P(ymachi),
    Real* zmac,AMREX_ARLIM_P(zmaclo),AMREX_ARLIM_P(zmachi),
    Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
    Real* maskunsplit,
   AMREX_ARLIM_P(maskunsplitlo),AMREX_ARLIM_P(maskunsplithi),
    const Real* xlo,const Real* dx,
    int* nmat,
    int* level,
    int* finest_level);

 
   void FORT_MAKETENSOR(
    int* ncomp_visc,int* im,
    const Real* xlo,const Real* dx,
    Real* visc,AMREX_ARLIM_P(visclo),AMREX_ARLIM_P(vischi),
    Real* tensor,AMREX_ARLIM_P(tlo),AMREX_ARLIM_P(thi),
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
    Real* snew,AMREX_ARLIM_P(snewlo),AMREX_ARLIM_P(snewhi),
    Real* fsi,AMREX_ARLIM_P(fsilo),AMREX_ARLIM_P(fsihi),
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
    Real* snew,AMREX_ARLIM_P(snewlo),AMREX_ARLIM_P(snewhi),
    Real* lsnew,AMREX_ARLIM_P(lsnewlo),AMREX_ARLIM_P(lsnewhi),
    Real* fsi,AMREX_ARLIM_P(fsilo),AMREX_ARLIM_P(fsihi),
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,int* bfact,
    int* nmat,int* nstate);

 
   void FORT_UPDATETENSOR(
    int* level,
    int* finest_level,
    int* nmat,int* im,int* ncomp_visc,
    Real* vof,AMREX_ARLIM_P(voflo),AMREX_ARLIM_P(vofhi),
    Real* visc,AMREX_ARLIM_P(visclo),AMREX_ARLIM_P(vischi),
    Real* tendata,AMREX_ARLIM_P(tendatalo),AMREX_ARLIM_P(tendatahi),
    const Real* dx,const Real* xlo,
    Real* vel,AMREX_ARLIM_P(vello),AMREX_ARLIM_P(velhi),
    Real* tensor_new,AMREX_ARLIM_P(tnewlo),AMREX_ARLIM_P(tnewhi),
    Real* tensor_old,AMREX_ARLIM_P(toldlo),AMREX_ARLIM_P(toldhi),
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
    Real* vof,AMREX_ARLIM_P(voflo),AMREX_ARLIM_P(vofhi),
    const Real* dx,const Real* xlo,
    Real* tensor_new,AMREX_ARLIM_P(tnewlo),AMREX_ARLIM_P(tnewhi),
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
    Real* vof,AMREX_ARLIM_P(voflo),AMREX_ARLIM_P(vofhi),
    const Real* dx,const Real* xlo,
    Real* tensor_new,AMREX_ARLIM_P(tnewlo),AMREX_ARLIM_P(tnewhi),
    Real* tensor_old,AMREX_ARLIM_P(toldlo),AMREX_ARLIM_P(toldhi),
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,int* bfact);

 
   void FORT_ESTDT(
      int* nsolveMM_FACE,
      int* enable_spectral,
      Real* elastic_time,
      int* microlayer_substrate,
      Real* microlayer_angle,
      Real* microlayer_size,
      Real* macrolayer_size,
      Real* latent_heat,
      Real* reaction_rate,
      int* freezing_model,
      int* Tanasawa_or_Schrage,
      int* distribute_from_target,
      Real* saturation_temp,
      int* mass_fraction_id,
      Real* molar_mass,
      Real* species_molar_mass,
      Real* species_evaporation_density,
      const Real* Umac, AMREX_ARLIM_P(umaclo),AMREX_ARLIM_P(umachi),
      const Real* Ucell, AMREX_ARLIM_P(ucelllo),AMREX_ARLIM_P(ucellhi),
      const Real* sol, AMREX_ARLIM_P(sollo),AMREX_ARLIM_P(solhi),
      Real* den,AMREX_ARLIM_P(denlo),AMREX_ARLIM_P(denhi),
      Real* vof,AMREX_ARLIM_P(voflo),AMREX_ARLIM_P(vofhi),
      Real* dist,AMREX_ARLIM_P(distlo),AMREX_ARLIM_P(disthi),
      const Real* xlo,const Real* dx,
      const int* tilelo, const int* tilehi,
      const int* fablo, const int* fabhi,
      int* bfact,
      Real* cap_wave_speed,
      Real* u_max,
      Real* u_max_estdt,
      Real* u_max_cap_wave,
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
    int* data_dir,
    int* law_of_the_wall,
    int* im_solid_map,
    int* level,
    int* finest_level,
    int* ngrow_law_of_wall,
    int* ngrow_distance,
    int* nmat,
    int* nparts,
    int* nparts_ghost,
    int* nden,
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,
    int* bfact,
    const Real* xlo,const Real* dx,
    Real* dt,
    Real* time,
    Real* LSCP,AMREX_ARLIM_P(LSCPlo),AMREX_ARLIM_P(LSCPhi), 
    Real* LSFD,AMREX_ARLIM_P(LSFDlo),AMREX_ARLIM_P(LSFDhi), 
    Real* state,AMREX_ARLIM_P(statelo),AMREX_ARLIM_P(statehi), 
    Real* ufluid,AMREX_ARLIM_P(ufluidlo),AMREX_ARLIM_P(ufluidhi),
    Real* usolid,AMREX_ARLIM_P(usolidlo),AMREX_ARLIM_P(usolidhi),
    Real* ughost,AMREX_ARLIM_P(ughostlo),AMREX_ARLIM_P(ughosthi),
    Real* history_dat,
    AMREX_ARLIM_P(history_datlo),AMREX_ARLIM_P(history_dathi),
    int* nhistory,
    Real* visc_coef);


   void FORT_STEFANSOLVER(
    int* project_option,
    int* solidheat_flag,
    Real* microlayer_size,
    int* microlayer_substrate,
    Real* microlayer_temperature_substrate,
    int* adjust_temperature,
    int* nmat,
    int* nten,
    int* nstate,
    int* ntsat,
    int* nden,
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
    Real* STATEFAB,AMREX_ARLIM_P(STATEFABlo),AMREX_ARLIM_P(STATEFABhi), 
    Real* TSATFAB,AMREX_ARLIM_P(TSATFABlo),AMREX_ARLIM_P(TSATFABhi), 
    Real* cellmm,AMREX_ARLIM_P(cellmmlo),AMREX_ARLIM_P(cellmmhi), 
    Real* xfacemm,AMREX_ARLIM_P(xfacemmlo),AMREX_ARLIM_P(xfacemmhi), 
    Real* yfacemm,AMREX_ARLIM_P(yfacemmlo),AMREX_ARLIM_P(yfacemmhi), 
    Real* zfacemm,AMREX_ARLIM_P(zfacemmlo),AMREX_ARLIM_P(zfacemmhi), 
    Real* swept,AMREX_ARLIM_P(sweptlo),AMREX_ARLIM_P(swepthi), 
    Real* LS,AMREX_ARLIM_P(LSlo),AMREX_ARLIM_P(LShi), 
    Real* thermal,AMREX_ARLIM_P(thermallo),AMREX_ARLIM_P(thermalhi), 
    Real* Snew,AMREX_ARLIM_P(Snewlo),AMREX_ARLIM_P(Snewhi), 
    Real* DeDT,AMREX_ARLIM_P(DeDTlo),AMREX_ARLIM_P(DeDThi),
    Real* den,AMREX_ARLIM_P(denlo),AMREX_ARLIM_P(denhi),
    Real* coeff,AMREX_ARLIM_P(coefflo),AMREX_ARLIM_P(coeffhi),
    Real* vol,AMREX_ARLIM_P(vollo),AMREX_ARLIM_P(volhi),
    Real* heatx,AMREX_ARLIM_P(heatxlo),AMREX_ARLIM_P(heatxhi),
    Real* heaty,AMREX_ARLIM_P(heatylo),AMREX_ARLIM_P(heatyhi),
    Real* heatz,AMREX_ARLIM_P(heatzlo),AMREX_ARLIM_P(heatzhi),
    Real* areax,AMREX_ARLIM_P(areaxlo),AMREX_ARLIM_P(areaxhi),
    Real* areay,AMREX_ARLIM_P(areaylo),AMREX_ARLIM_P(areayhi),
    Real* areaz,AMREX_ARLIM_P(areazlo),AMREX_ARLIM_P(areazhi) ); 


   void FORT_HEATSOURCE_FACE(
    int* nmat,int* nten,int* nstate,
    Real* latent_heat,
    Real* saturation_temp,
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,
    int* bfact,
    const Real* xlo,const Real* dx,
    Real* dt,
    Real* time,
    int* level,
    int* finest_level,
    Real* LS,AMREX_ARLIM_P(LSlo),AMREX_ARLIM_P(LShi), 
    Real* Snew,AMREX_ARLIM_P(Snewlo),AMREX_ARLIM_P(Snewhi), 
    Real* DeDT,AMREX_ARLIM_P(DeDTlo),AMREX_ARLIM_P(DeDThi),
    Real* den,AMREX_ARLIM_P(denlo),AMREX_ARLIM_P(denhi),
    Real* vol,AMREX_ARLIM_P(vollo),AMREX_ARLIM_P(volhi),
    Real* heatx,AMREX_ARLIM_P(heatxlo),AMREX_ARLIM_P(heatxhi),
    Real* heaty,AMREX_ARLIM_P(heatylo),AMREX_ARLIM_P(heatyhi),
    Real* heatz,AMREX_ARLIM_P(heatzlo),AMREX_ARLIM_P(heatzhi),
    Real* areax,AMREX_ARLIM_P(areaxlo),AMREX_ARLIM_P(areaxhi),
    Real* areay,AMREX_ARLIM_P(areaylo),AMREX_ARLIM_P(areayhi),
    Real* areaz,AMREX_ARLIM_P(areazlo),AMREX_ARLIM_P(areazhi) ); 


   void FORT_HEATSOURCE(
    int* nstate,
    int* nmat,
    int* nden,
    const Real* xlo,const Real* dx,
    Real* temperature_source,
    Real* temperature_source_cen,
    Real* temperature_source_rad,
    Real* rhoinverse,
   AMREX_ARLIM_P(rhoinverselo),AMREX_ARLIM_P(rhoinversehi), 
    Real* DeDTinverse,
   AMREX_ARLIM_P(DeDTinverselo),AMREX_ARLIM_P(DeDTinversehi), 
    Real* Tnew,AMREX_ARLIM_P(Tnewlo),AMREX_ARLIM_P(Tnewhi), 
    Real* lsfab,AMREX_ARLIM_P(lsfablo),AMREX_ARLIM_P(lsfabhi), 
    Real* recon,AMREX_ARLIM_P(reconlo),AMREX_ARLIM_P(reconhi), 
    Real* vol,AMREX_ARLIM_P(vollo),AMREX_ARLIM_P(volhi),
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
    Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi),
    Real* jump,AMREX_ARLIM_P(jumplo),AMREX_ARLIM_P(jumphi),
    Real* mdot,AMREX_ARLIM_P(mdotlo),AMREX_ARLIM_P(mdothi),
    Real* LSnew,AMREX_ARLIM_P(LSnewlo),AMREX_ARLIM_P(LSnewhi));


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
    Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi),
    Real* xface,AMREX_ARLIM_P(xfacelo),AMREX_ARLIM_P(xfacehi),
    Real* yface,AMREX_ARLIM_P(yfacelo),AMREX_ARLIM_P(yfacehi),
    Real* zface,AMREX_ARLIM_P(zfacelo),AMREX_ARLIM_P(zfacehi),
    Real* LSnew,AMREX_ARLIM_P(LSnewlo),AMREX_ARLIM_P(LSnewhi),
    Real* recon,AMREX_ARLIM_P(reconlo),AMREX_ARLIM_P(reconhi));


   void FORT_RENORM_MDOT(
    Real* mdotplus,Real* mdotminus,Real* mdotcount,
    int* level,int* finest_level,
    int* nmat,int* nten,
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,
    int* bfact,
    const Real* xlo,const Real* dx,Real* dt,
    Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi),
    Real* mdot,AMREX_ARLIM_P(mdotlo),AMREX_ARLIM_P(mdothi));


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
    Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi),
    Real* tag,AMREX_ARLIM_P(tablo),AMREX_ARLIM_P(taghi),
    Real* jump,AMREX_ARLIM_P(jumplo),AMREX_ARLIM_P(jumphi),
    Real* LSnew,AMREX_ARLIM_P(LSnewlo),AMREX_ARLIM_P(LSnewhi),
    Real* recon,AMREX_ARLIM_P(reconlo),AMREX_ARLIM_P(reconhi));


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
    Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi),
    Real* LS,AMREX_ARLIM_P(LSlo),AMREX_ARLIM_P(LShi),
    Real* tag,AMREX_ARLIM_P(tablo),AMREX_ARLIM_P(taghi),
    Real* jump,AMREX_ARLIM_P(jumplo),AMREX_ARLIM_P(jumphi));

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
    Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi),
    Real* tag,AMREX_ARLIM_P(taglo),AMREX_ARLIM_P(taghi),
    Real* jump,AMREX_ARLIM_P(jumplo),AMREX_ARLIM_P(jumphi));


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
     int* VOF_reflux,
     int* dir_counter,
     int* normdir,
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,
     int* bfact,
     int* bfact_f,
     Real* dt,
     Real* time,
     Real* passive_veltime,
     Real* vofflux,AMREX_ARLIM_P(voffluxlo),AMREX_ARLIM_P(voffluxhi),
      // original data
     Real* LS,AMREX_ARLIM_P(LSlo),AMREX_ARLIM_P(LShi),
     Real* den,AMREX_ARLIM_P(denlo),AMREX_ARLIM_P(denhi),
     Real* tensor,AMREX_ARLIM_P(tensorlo),AMREX_ARLIM_P(tensorhi),
     Real* velfab,AMREX_ARLIM_P(velfablo),AMREX_ARLIM_P(velfabhi),
      // slope data
     Real* vofslope,AMREX_ARLIM_P(vofslopelo),AMREX_ARLIM_P(vofslopehi), 
      // this is the result
     Real* snew,AMREX_ARLIM_P(snewlo),AMREX_ARLIM_P(snewhi), 
     Real* tennew,AMREX_ARLIM_P(tennewlo),AMREX_ARLIM_P(tennewhi), 
     Real* LSnew,AMREX_ARLIM_P(LSnewlo),AMREX_ARLIM_P(LSnewhi), 
      // other vars.
     Real* ucell,AMREX_ARLIM_P(ucelllo),AMREX_ARLIM_P(ucellhi),
     Real* vofls0,AMREX_ARLIM_P(vofls0lo),AMREX_ARLIM_P(vofls0hi),
     Real* conscor,AMREX_ARLIM_P(conscorlo),AMREX_ARLIM_P(conscorhi),
     Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
     Real* masknbr,AMREX_ARLIM_P(masknbrlo),AMREX_ARLIM_P(masknbrhi),
     Real* unode,AMREX_ARLIM_P(unodelo),AMREX_ARLIM_P(unodehi),
     const Real* xlo,const Real* dx,
      // local variables
     Real* conserve,AMREX_ARLIM_P(conservelo),AMREX_ARLIM_P(conservehi),
     Real* xvel,AMREX_ARLIM_P(xvello),AMREX_ARLIM_P(xvelhi),
     Real* yvel,AMREX_ARLIM_P(yvello),AMREX_ARLIM_P(yvelhi),
     Real* zvel,AMREX_ARLIM_P(zvello),AMREX_ARLIM_P(zvelhi),
     Real* xvelslope,AMREX_ARLIM_P(xvelslopelo),AMREX_ARLIM_P(xvelslopehi),
     Real* yvelslope,AMREX_ARLIM_P(yvelslopelo),AMREX_ARLIM_P(yvelslopehi),
     Real* zvelslope,AMREX_ARLIM_P(zvelslopelo),AMREX_ARLIM_P(zvelslopehi),
     Real* momslope,AMREX_ARLIM_P(momslopelo),AMREX_ARLIM_P(momslopehi),
     Real* xmomside,AMREX_ARLIM_P(xmomsidelo),AMREX_ARLIM_P(xmomsidehi),
     Real* ymomside,AMREX_ARLIM_P(ymomsidelo),AMREX_ARLIM_P(ymomsidehi),
     Real* zmomside,AMREX_ARLIM_P(zmomsidelo),AMREX_ARLIM_P(zmomsidehi),
     Real* xmassside,AMREX_ARLIM_P(xmasssidelo),AMREX_ARLIM_P(xmasssidehi),
     Real* ymassside,AMREX_ARLIM_P(ymasssidelo),AMREX_ARLIM_P(ymasssidehi),
     Real* zmassside,AMREX_ARLIM_P(zmasssidelo),AMREX_ARLIM_P(zmasssidehi),
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
     Real* LS,AMREX_ARLIM_P(LSlo),AMREX_ARLIM_P(LShi),
     Real* den,AMREX_ARLIM_P(denlo),AMREX_ARLIM_P(denhi),
     Real* tensor,AMREX_ARLIM_P(tensorlo),AMREX_ARLIM_P(tensorhi),
     Real* velfab,AMREX_ARLIM_P(velfablo),AMREX_ARLIM_P(velfabhi),
      // slope data
     Real* vofslope,AMREX_ARLIM_P(vofslopelo),AMREX_ARLIM_P(vofslopehi), 
      // this is the result
     Real* snew,AMREX_ARLIM_P(snewlo),AMREX_ARLIM_P(snewhi), 
     Real* tennew,AMREX_ARLIM_P(tennewlo),AMREX_ARLIM_P(tennewhi), 
     Real* LSnew,AMREX_ARLIM_P(LSnewlo),AMREX_ARLIM_P(LSnewhi), 
      // other vars.
     Real* ucell,AMREX_ARLIM_P(ucelllo),AMREX_ARLIM_P(ucellhi),
     Real* vofls0,AMREX_ARLIM_P(vofls0lo),AMREX_ARLIM_P(vofls0hi),
     Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
     Real* masknbr,AMREX_ARLIM_P(masknbrlo),AMREX_ARLIM_P(masknbrhi),
     Real* maskunsplit,
    AMREX_ARLIM_P(maskunsplitlo),AMREX_ARLIM_P(maskunsplithi),
     Real* unode,AMREX_ARLIM_P(unodelo),AMREX_ARLIM_P(unodehi),
     Real* vnode,AMREX_ARLIM_P(vnodelo),AMREX_ARLIM_P(vnodehi),
     Real* wnode,AMREX_ARLIM_P(wnodelo),AMREX_ARLIM_P(wnodehi),
     const Real* xlo,const Real* dx,
      // local variables
     Real* conserve,AMREX_ARLIM_P(conservelo),AMREX_ARLIM_P(conservehi),
     Real* xvel,AMREX_ARLIM_P(xvello),AMREX_ARLIM_P(xvelhi),
     Real* yvel,AMREX_ARLIM_P(yvello),AMREX_ARLIM_P(yvelhi),
     Real* zvel,AMREX_ARLIM_P(zvello),AMREX_ARLIM_P(zvelhi),
     Real* xmomside,AMREX_ARLIM_P(xmomsidelo),AMREX_ARLIM_P(xmomsidehi),
     Real* ymomside,AMREX_ARLIM_P(ymomsidelo),AMREX_ARLIM_P(ymomsidehi),
     Real* zmomside,AMREX_ARLIM_P(zmomsidelo),AMREX_ARLIM_P(zmomsidehi),
     Real* xmassside,AMREX_ARLIM_P(xmasssidelo),AMREX_ARLIM_P(xmasssidehi),
     Real* ymassside,AMREX_ARLIM_P(ymasssidelo),AMREX_ARLIM_P(ymasssidehi),
     Real* zmassside,AMREX_ARLIM_P(zmasssidelo),AMREX_ARLIM_P(zmasssidehi),
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
     int* spec_material_id_AMBIENT,
     Real* species_evaporation_density,
     int* presbc,
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,
     int* bfact,
     Real* dt,
     Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
     Real* masknbr,AMREX_ARLIM_P(masknbrlo),AMREX_ARLIM_P(masknbrhi),
     Real* vol,AMREX_ARLIM_P(vollo),AMREX_ARLIM_P(volhi),
     Real* dennew,AMREX_ARLIM_P(dennewlo),AMREX_ARLIM_P(dennewhi), 
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
     Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi),
     Real* solxfab,
     AMREX_ARLIM_P(solxfablo),AMREX_ARLIM_P(solxfabhi),
     Real* solyfab,
     AMREX_ARLIM_P(solyfablo),AMREX_ARLIM_P(solyfabhi),
     Real* solzfab,
     AMREX_ARLIM_P(solzfablo),AMREX_ARLIM_P(solzfabhi),
     Real* LSNEW,AMREX_ARLIM_P(LSNEWlo),AMREX_ARLIM_P(LSNEWhi),
     Real* LS,AMREX_ARLIM_P(LSlo),AMREX_ARLIM_P(LShi),
     Real* vof,AMREX_ARLIM_P(voflo),AMREX_ARLIM_P(vofhi),
     Real* cellfab,AMREX_ARLIM_P(cellfablo),AMREX_ARLIM_P(cellfabhi),
     Real* newcell,AMREX_ARLIM_P(newcelllo),AMREX_ARLIM_P(newcellhi),
     Real* state,AMREX_ARLIM_P(statelo),AMREX_ARLIM_P(statehi),
     int* velbc,
     int* listbc,
     const Real* xlo,const Real* dx,
     Real* cur_time);


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
     Real* vof,AMREX_ARLIM_P(voflo),AMREX_ARLIM_P(vofhi),
     Real* mac,AMREX_ARLIM_P(maclo),AMREX_ARLIM_P(machi),
     Real* xface,AMREX_ARLIM_P(xfacelo),AMREX_ARLIM_P(xfacehi),
     Real* LS,AMREX_ARLIM_P(LSlo),AMREX_ARLIM_P(LShi),
     Real* solfab,
     AMREX_ARLIM_P(solfablo),AMREX_ARLIM_P(solfabhi),
     const Real* xlo,const Real* dx,
     int* dir,
     Real* cur_time);

   void FORT_CROSSTERM(
     int* nsolveMM_FACE,
     int* tileloop,
     int* dir,
     int* operation_flag, // 8
     int* enable_spectral,
     int* spectral_loop,
     int* ncfluxreg,
     Real* semflux,AMREX_ARLIM_P(semfluxlo),AMREX_ARLIM_P(semfluxhi),
     Real* mask,AMREX_ARLIM_P(masklo),AMREX_ARLIM_P(maskhi),
     Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi),
     Real* faceLS,AMREX_ARLIM_P(faceLSlo),AMREX_ARLIM_P(faceLShi),
     Real* mdata,AMREX_ARLIM_P(mdatalo),AMREX_ARLIM_P(mdatahi),
     Real* tdata,AMREX_ARLIM_P(tdatalo),AMREX_ARLIM_P(tdatahi),
     Real* c_tdata,AMREX_ARLIM_P(c_tdatalo),AMREX_ARLIM_P(c_tdatahi),
     Real* maskSEM,AMREX_ARLIM_P(maskSEMlo),AMREX_ARLIM_P(maskSEMhi),
     const Real* xlo,const Real* dx,
     Real* dt,
     Real* vel,AMREX_ARLIM_P(vello),AMREX_ARLIM_P(velhi),
     Real* levelpc,
    AMREX_ARLIM_P(levelpclo),AMREX_ARLIM_P(levelpchi),
     Real* xflux,AMREX_ARLIM_P(xfluxlo),AMREX_ARLIM_P(xfluxhi),
     Real* xface,AMREX_ARLIM_P(xfacelo),AMREX_ARLIM_P(xfacehi),
     Real* xfacemm,AMREX_ARLIM_P(xfacemmlo),AMREX_ARLIM_P(xfacemmhi), 
     Real* xcellmm,AMREX_ARLIM_P(xcellmmlo),AMREX_ARLIM_P(xcellmmhi), 
     Real* recon,AMREX_ARLIM_P(reconlo),AMREX_ARLIM_P(reconhi), 
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
     Real* semflux,AMREX_ARLIM_P(semfluxlo),AMREX_ARLIM_P(semfluxhi),
     Real* amrsync,AMREX_ARLIM_P(amrsynclo),AMREX_ARLIM_P(amrsynchi),
     Real* mask0,AMREX_ARLIM_P(mask0lo),AMREX_ARLIM_P(mask0hi),
     Real* mask3,AMREX_ARLIM_P(mask3lo),AMREX_ARLIM_P(mask3hi),
     Real* maskSEM,AMREX_ARLIM_P(maskSEMlo),AMREX_ARLIM_P(maskSEMhi),
     Real* faceLS,AMREX_ARLIM_P(faceLSlo),AMREX_ARLIM_P(faceLShi),
     Real* mdata,AMREX_ARLIM_P(mdatalo),AMREX_ARLIM_P(mdatahi),
     Real* tdata,AMREX_ARLIM_P(tdatalo),AMREX_ARLIM_P(tdatahi),
     Real* c_tdata,AMREX_ARLIM_P(c_tdatalo),AMREX_ARLIM_P(c_tdatahi),
     Real* vel,AMREX_ARLIM_P(vello),AMREX_ARLIM_P(velhi),
     Real* solidx,
     AMREX_ARLIM_P(solidxlo),AMREX_ARLIM_P(solidxhi),
     Real* solidy,
     AMREX_ARLIM_P(solidylo),AMREX_ARLIM_P(solidyhi),
     Real* solidz,
     AMREX_ARLIM_P(solidzlo),AMREX_ARLIM_P(solidzhi),
     Real* levelpc,
     AMREX_ARLIM_P(levelpclo),AMREX_ARLIM_P(levelpchi),
     Real* recon,
     AMREX_ARLIM_P(reconlo),AMREX_ARLIM_P(reconhi),
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
     int* SEM_advection_algorithm,
     int* simple_AMR_BC_flag_viscosity);


   void FORT_VISCTENSORHEAT(
     int* ntensor,
     int* ntensorMM,
     int* nsolve,
     int* nsolveMM,
     int* nsolveMM_FACE,
     int* nstate,
     const Real* xlo,const Real* dx,
     Real* lsfab,AMREX_ARLIM_P(lsfablo),AMREX_ARLIM_P(lsfabhi),
     Real* DeDTinverse,AMREX_ARLIM_P(DeDTinverselo),AMREX_ARLIM_P(DeDTinversehi),
     Real* heatnew,AMREX_ARLIM_P(heatnewlo),AMREX_ARLIM_P(heatnewhi),
     Real* xstress,AMREX_ARLIM_P(xstresslo),AMREX_ARLIM_P(xstresshi),
     Real* ystress,AMREX_ARLIM_P(ystresslo),AMREX_ARLIM_P(ystresshi),
     Real* zstress,AMREX_ARLIM_P(zstresslo),AMREX_ARLIM_P(zstresshi),
     Real* gradu,AMREX_ARLIM_P(gradulo),AMREX_ARLIM_P(graduhi),
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
     Real* solxfab,
     AMREX_ARLIM_P(solxfablo),AMREX_ARLIM_P(solxfabhi),
     Real* solyfab,
     AMREX_ARLIM_P(solyfablo),AMREX_ARLIM_P(solyfabhi),
     Real* solzfab,
     AMREX_ARLIM_P(solzfablo),AMREX_ARLIM_P(solzfabhi),
     Real* snew,AMREX_ARLIM_P(snewlo),AMREX_ARLIM_P(snewhi),
     Real* lsnew,AMREX_ARLIM_P(lsnewlo),AMREX_ARLIM_P(lsnewhi),
     Real* du,AMREX_ARLIM_P(dulo),AMREX_ARLIM_P(duhi),
     Real* advect,AMREX_ARLIM_P(advectlo),AMREX_ARLIM_P(advecthi),
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
     Real* vel_time,
     int* dir_absolute_direct_split,
     int* normdir,
     Real* unodetemp,AMREX_ARLIM_P(unodetemplo),AMREX_ARLIM_P(unodetemphi),
     Real* unode,AMREX_ARLIM_P(unodelo),AMREX_ARLIM_P(unodehi),
     Real* ucell,AMREX_ARLIM_P(ucelllo),AMREX_ARLIM_P(ucellhi),
     const Real* xlo,const Real* dx,
     int* mac_grow,
     int* map_forward,
     int* level,
     int* finest_level,
     int* SDC_outer_sweeps,
     int* ns_time_order,
     int* divu_outer_sweeps,
     int* num_divu_outer_sweeps);


   void fort_assimilate_tensor_from_particles(
     int* tid,
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,
     int* bfact,
     int* level,
     int* finest_level,
     const Real* xlo,const Real* dx,
     const amrex::Real* particles,
     const amrex::Real* Nparticles,
     int Np, // must be pass by value just in case "long" on the c++ side.
     int Nn, // must be pass by value just in case "long" on the c++ side.
     int* nTNEW,
     int* matrix_points,
     int* RHS_points,
     int* ncomp_accumulate,
     int* ipart_type,
     Real* TNEWfab,
     AMREX_ARLIM_P(TNEWfablo),
     AMREX_ARLIM_P(TNEWfabhi),
     Real* matrixfab,
     AMREX_ARLIM_P(matrixfablo),
     AMREX_ARLIM_P(matrixfabhi));

}

}/* namespace amrex */
#endif

#endif 
