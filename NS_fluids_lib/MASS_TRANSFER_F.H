#ifndef _MASS_TRANSFER_F_H_
#define _MASS_TRANSFER_F_H_

#ifdef BL_LANG_FORT

#else

#include <AMReX_ArrayLim.H>

namespace amrex{

extern "C"
{

   void fort_convertmaterial(
    int* tid,
    int* im_outer,
    int* im_opp_outer,
    int* ngrow_expansion,
    int* level,int* finest_level,
    int* nmat,
    int* nten,
    int* nden,
    int* nstate,
    int* ntsat,
    int* supermesh_flag,
    Real* latent_heat,
    Real* saturation_temp,
    int* freezing_model,
    int* Tanasawa_or_Schrage_or_Kassemi,
    int* mass_fraction_id,
    int* distribute_from_target,
    int* constant_density_all_time,
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,
    int* bfact,
    Real* min_stefan_velocity_for_dt,
    int* vofbc,
    const Real* xlo,const Real* dx,
    Real* dt,
    Real* delta_mass,
    Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi), 
    Real* conductstate,
    AMREX_ARLIM_P(conductstatelo),
    AMREX_ARLIM_P(conductstatehi), 
    Real* nodevel,AMREX_ARLIM_P(nodevello),AMREX_ARLIM_P(nodevelhi), 
    Real* jumpfab,AMREX_ARLIM_P(jumpfablo),AMREX_ARLIM_P(jumpfabhi),
    Real* tsatfab,AMREX_ARLIM_P(tsatfablo),AMREX_ARLIM_P(tsatfabhi),
    Real* LSold,AMREX_ARLIM_P(LSoldlo),AMREX_ARLIM_P(LSoldhi),
    Real* LSnew,AMREX_ARLIM_P(LSnewlo),AMREX_ARLIM_P(LSnewhi),      
    Real* recon,AMREX_ARLIM_P(reconlo),AMREX_ARLIM_P(reconhi),      
    Real* snew,AMREX_ARLIM_P(snewlo),AMREX_ARLIM_P(snewhi),      
    Real* EOS,AMREX_ARLIM_P(EOSlo),AMREX_ARLIM_P(EOShi),
    Real* swept,AMREX_ARLIM_P(sweptlo),AMREX_ARLIM_P(swepthi)); 


   void fort_nodedisplace(
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
     Real* oldLS,AMREX_ARLIM_P(oldLSlo),AMREX_ARLIM_P(oldLShi),
     const Real* xlo,const Real* dx,
     int* level,int* finest_level);

   void fort_ratemasschange(
    int* tid,
    int* nucleation_flag,
    int* stefan_flag,
    int* level,
    int* finest_level,
    int* ngrow_distance,
    int* nstate,
    int* nmat,
    int* nten,
    int* nburning,
    int* ntsat,
    int* nden,
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
    Real* R_Palmore_Desjardins,
    Real* smoothing_length_scale,
    Real* latent_heat,
    int* use_exact_temperature,
    Real* reaction_rate,
    Real* hardwire_Y_gamma,
    Real* hardwire_T_gamma,
    Real* accommodation_coefficient,
    Real* reference_pressure,
    Real* saturation_temp,
    Real* saturation_temp_curv,
    Real* saturation_temp_vel,
    Real* saturation_temp_min,
    Real* saturation_temp_max,
    int* freezing_model,
    int* Tanasawa_or_Schrage_or_Kassemi,
    int* interface_mass_transfer_model,
    int* distribute_from_target,
    int* mass_fraction_id,
    int* constant_density_all_time,
    int* material_type_evap, 
    Real* molar_mass,
    Real* species_molar_mass,
    int* use_supermesh,
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,
    int* bfact,
    const Real* xlo,const Real* dx,
    Real* time,
    Real* dt,
    int* arraysize,
    Real* blob_array,
    int* color_count,
    Real* colorfab,AMREX_ARLIM_P(colorfablo),AMREX_ARLIM_P(colorfabhi), 
    Real* typefab,AMREX_ARLIM_P(typefablo),AMREX_ARLIM_P(typefabhi), 
    Real* maskcov,AMREX_ARLIM_P(maskcovlo),AMREX_ARLIM_P(maskcovhi), 
    Real* conductstate,
    AMREX_ARLIM_P(conductstatelo),
    AMREX_ARLIM_P(conductstatehi), 
    Real* burnvel,AMREX_ARLIM_P(burnvello),AMREX_ARLIM_P(burnvelhi), 
    Real* Tsatfab,AMREX_ARLIM_P(Tsatfablo),AMREX_ARLIM_P(Tsatfabhi), 
    Real* smoothfab,AMREX_ARLIM_P(smoothfablo),AMREX_ARLIM_P(smoothfabhi), 
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


   void fort_extend_burning_vel(
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

   void fort_extend_drag(
    int* level,
    int* finest_level,
    const Real* xlo,const Real* dx,
    int* nmat,
    int* nten,
    int* ncomp,
    int* ngrow,
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,
    int* bfact,
    Real* drag,AMREX_ARLIM_P(draglo),AMREX_ARLIM_P(draghi), 
    Real* LS,AMREX_ARLIM_P(LSlo),AMREX_ARLIM_P(LShi));

}

}/* namespace amrex */

#endif

#endif 
