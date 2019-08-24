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

#include <ArrayLim.H>

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
    Real* maskcov,ARLIM_P(maskcovlo),ARLIM_P(maskcovhi), 
    Real* deltaVOF,ARLIM_P(deltaVOFlo),ARLIM_P(deltaVOFhi), 
    Real* nodevel,ARLIM_P(nodevello),ARLIM_P(nodevelhi), 
    Real* jump,ARLIM_P(jumplo),ARLIM_P(jumphi),
    Real* LSold,ARLIM_P(LSoldlo),ARLIM_P(LSoldhi),
    Real* LSnew,ARLIM_P(LSnewlo),ARLIM_P(LSnewhi),      
    Real* recon,ARLIM_P(reconlo),ARLIM_P(reconhi),      
    Real* snew,ARLIM_P(snewlo),ARLIM_P(snewhi),      
    Real* EOS,ARLIM_P(EOSlo),ARLIM_P(EOShi),
    Real* swept,ARLIM_P(sweptlo),ARLIM_P(swepthi)); 


   void FORT_NODEDISPLACE(
     int* nmat,
     int* nten,
     int* nburning,
     const int* tilelo,const int* tilehi,
     const int* fablo,const int* fabhi,
     int* bfact,
     int* velbc,
     Real* dt,
     Real* unode,ARLIM_P(unodelo),ARLIM_P(unodehi),
     Real* ucell,ARLIM_P(ucelllo),ARLIM_P(ucellhi),
     const Real* xlo,const Real* dx,
     int* level,int* finest_level);

   void FORT_RATEMASSCHANGE(
    int* stefan_flag,
    int* level,
    int* finest_level,
    int* normal_probe_size,
    int* ngrow_distance,
    int* nmat,
    int* nten,
    int* nburning,
    int* nden,
    Real* density_floor_expansion,
    Real* density_ceiling_expansion,
    int* microlayer_substrate,
    Real* microlayer_angle,
    Real* microlayer_size,
    Real* macrolayer_size,
    Real* max_contact_line_size,
    Real* latent_heat,
    int* use_exact_temperature,
    Real* reaction_rate,
    Real* saturation_temp,
    int* freezing_model,
    int* distribute_from_target,
    int* mass_fraction_id,
    Real* species_evaporation_density,
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
    Real* colorfab,ARLIM_P(colorfablo),ARLIM_P(colorfabhi), 
    Real* typefab,ARLIM_P(typefablo),ARLIM_P(typefabhi), 
    Real* maskcov,ARLIM_P(maskcovlo),ARLIM_P(maskcovhi), 
    Real* burnvel,ARLIM_P(burnvello),ARLIM_P(burnvelhi), 
    Real* LS,ARLIM_P(LSlo),ARLIM_P(LShi), 
    Real* LSnew,ARLIM_P(LSnewlo),ARLIM_P(LSnewhi), 
    Real* LS_slopes_FD,
    ARLIM_P(LS_slopes_FDlo),
    ARLIM_P(LS_slopes_FDhi), 
    Real* EOS,ARLIM_P(EOSlo),ARLIM_P(EOShi),
    Real* recon,ARLIM_P(reconlo),ARLIM_P(reconhi),
    Real* pres,ARLIM_P(preslo),ARLIM_P(preshi)); 


   void FORT_EXTEND_BURNING_VEL(
    int* level,int* finest_level,
    const Real* xlo,const Real* dx,
    int* nmat,
    int* nten,
    int* nburning,
    int* ngrow,
    Real* latent_heat,
    const int* tilelo,const int* tilehi,
    const int* fablo,const int* fabhi,
    int* bfact,
    Real* burnvel,ARLIM_P(burnvello),ARLIM_P(burnvelhi), 
    Real* LS,ARLIM_P(LSlo),ARLIM_P(LShi));

}
#endif

#endif 
