#ifndef _MASS_TRANSFER_F_H_
#define _MASS_TRANSFER_F_H_

#ifdef BL_LANG_FORT
#    define FORT_CONVERTMATERIAL convertmaterial
#    define FORT_NODEDISPLACE nodedisplace
#    define FORT_RATEMASSCHANGE ratemasschange
#    define FORT_EXTEND_BURNING_VEL extend_burning_vel
#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_CONVERTMATERIAL CONVERTMATERIAL
#    define FORT_NODEDISPLACE NODEDISPLACE
#    define FORT_RATEMASSCHANGE RATEMASSCHANGE
#    define FORT_EXTEND_BURNING_VEL EXTEND_BURNING_VEL
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_CONVERTMATERIAL convertmaterial
#    define FORT_NODEDISPLACE nodedisplace
#    define FORT_RATEMASSCHANGE ratemasschange
#    define FORT_EXTEND_BURNING_VEL extend_burning_vel
#  else
#    define FORT_CONVERTMATERIAL convertmaterial_
#    define FORT_NODEDISPLACE nodedisplace_
#    define FORT_RATEMASSCHANGE ratemasschange_
#    define FORT_EXTEND_BURNING_VEL extend_burning_vel_
#  endif
#  endif

#include <AMReX_ArrayLim.H>

namespace amrex{

extern "C"
{

   void FORT_CONVERTMATERIAL(
    int* tid,
    int* isweep,
    int* solvability_projection,
    int* ngrow_expansion,
    int* level,int* finest_level,
    int* normal_probe_size,
    int* nmat,
    int* nten,
    int* nden,
    int* nstate,
    int* ntsat,
    Real* density_floor_expansion,
    Real* density_ceiling_expansion,
    Real* latent_heat,
    Real* saturation_temp,
    int* freezing_model,
    int* mass_fraction_id,
    Real* species_evaporation_density,
    int* distribute_from_target,
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,
    int* bfact,
    int* vofbc,
    const Real* xlo,const Real* dx,
    Real* dt,
    Real* delta_mass,
    Real* DVOF,
    Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi), 
    Real* deltaVOF,AMREX_ARLIM_P(deltaVOFlo),AMREX_ARLIM_P(deltaVOFhi), 
    Real* nodevel,AMREX_ARLIM_P(nodevello),AMREX_ARLIM_P(nodevelhi), 
    Real* jumpfab,AMREX_ARLIM_P(jumpfablo),AMREX_ARLIM_P(jumpfabhi),
    Real* tsatfab,AMREX_ARLIM_P(tsatfablo),AMREX_ARLIM_P(tsatfabhi),
    Real* LSold,AMREX_ARLIM_P(LSoldlo),AMREX_ARLIM_P(LSoldhi),
    Real* LSnew,AMREX_ARLIM_P(LSnewlo),AMREX_ARLIM_P(LSnewhi),      
    Real* recon,AMREX_ARLIM_P(reconlo),AMREX_ARLIM_P(reconhi),      
    Real* snew,AMREX_ARLIM_P(snewlo),AMREX_ARLIM_P(snewhi),      
    Real* EOS,AMREX_ARLIM_P(EOSlo),AMREX_ARLIM_P(EOShi),
    Real* swept,AMREX_ARLIM_P(sweptlo),AMREX_ARLIM_P(swepthi)); 


   void FORT_NODEDISPLACE(
     int* nmat,
     int* nten,
     int* nburning,
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,
     int* bfact,
     int* velbc,
     Real* dt,
     Real* unode,AMREX_ARLIM_P(unodelo),AMREX_ARLIM_P(unodehi),
     Real* ucell,AMREX_ARLIM_P(ucelllo),AMREX_ARLIM_P(ucellhi),
     const Real* xlo,const Real* dx,
     int* level,int* finest_level);

   void FORT_RATEMASSCHANGE(
    int* tid,
    int* nucleation_flag,
    int* stefan_flag,
    int* level,
    int* finest_level,
    int* normal_probe_size,
    int* ngrow_distance,
    int* nstate,
    int* nmat,
    int* nten,
    int* nburning,
    int* ntsat,
    int* nden,
    Real* density_floor_expansion,
    Real* density_ceiling_expansion,
    int* custom_nucleation_model,
    int* do_the_nucleate,
    Real* nucleate_pos,
    const int* nucleate_pos_size,
    Real* nucleation_temp,
    Real* nucleation_pressure,
    Real* nucleation_pmg,
    Real* nucleation_mach,
    Real* cavitation_pressure,
    Real* cavitation_vapor_density,
    Real* cavitation_tension,
    int* microlayer_substrate,
    Real* microlayer_angle,
    Real* microlayer_size,
    Real* macrolayer_size,
    Real* max_contact_line_size,
    Real* latent_heat,
    int* use_exact_temperature,
    Real* reaction_rate,
    Real* saturation_temp,
    Real* saturation_temp_curv,
    Real* saturation_temp_vel,
    Real* saturation_temp_min,
    Real* saturation_temp_max,
    int* freezing_model,
    int* Tanasawa_or_Schrage,
    int* distribute_from_target,
    int* mass_fraction_id,
    Real* species_evaporation_density,
    Real* molar_mass,
    Real* species_molar_mass,
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,
    int* bfact,
    const Real* xlo,const Real* dx,
    Real* time,
    Real* dt,
    int* arraysize,
    Real* blob_array,
    int* num_elements_blobclass,
    int* color_count,
    Real* colorfab,AMREX_ARLIM_P(colorfablo),AMREX_ARLIM_P(colorfabhi), 
    Real* typefab,AMREX_ARLIM_P(typefablo),AMREX_ARLIM_P(typefabhi), 
    Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi), 
    Real* burnvel,AMREX_ARLIM_P(burnvello),AMREX_ARLIM_P(burnvelhi), 
    Real* Tsatfab,AMREX_ARLIM_P(Tsatfablo),AMREX_ARLIM_P(Tsatfabhi), 
    Real* LS,AMREX_ARLIM_P(LSlo),AMREX_ARLIM_P(LShi), 
    Real* LSnew,AMREX_ARLIM_P(LSnewlo),AMREX_ARLIM_P(LSnewhi), 
    Real* Snew,AMREX_ARLIM_P(Snewlo),AMREX_ARLIM_P(Snewhi), 
    Real* LS_slopes_FD,
    AMREX_ARLIM_P(LS_slopes_FDlo),
    AMREX_ARLIM_P(LS_slopes_FDhi), 
    Real* EOS,AMREX_ARLIM_P(EOSlo),AMREX_ARLIM_P(EOShi),
    Real* recon,AMREX_ARLIM_P(reconlo),AMREX_ARLIM_P(reconhi),
    Real* pres,AMREX_ARLIM_P(preslo),AMREX_ARLIM_P(preshi),
    Real* pres_eos,AMREX_ARLIM_P(pres_eoslo),AMREX_ARLIM_P(pres_eoshi),
    Real* curvfab,
    AMREX_ARLIM_P(curvfablo),
    AMREX_ARLIM_P(curvfabhi) ); 


   void FORT_EXTEND_BURNING_VEL(
    int* velflag,
    int* level,
    int* finest_level,
    const Real* xlo,const Real* dx,
    int* nmat,
    int* nten,
    int* nburning,
    int* ngrow,
    Real* latent_heat,
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,
    int* bfact,
    Real* burnvel,AMREX_ARLIM_P(burnvello),AMREX_ARLIM_P(burnvelhi), 
    Real* LS,AMREX_ARLIM_P(LSlo),AMREX_ARLIM_P(LShi));

}

}/* namespace amrex */

#endif

#endif 
