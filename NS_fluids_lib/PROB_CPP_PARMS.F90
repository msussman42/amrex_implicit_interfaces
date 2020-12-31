#undef BL_LANG_CC
#define BL_LANG_FORT
        
#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"

#include "PROB_CPP_PARMS_F.H"

#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else  
print *,"dimension bust"
stop
#endif



      subroutine FORT_OVERRIDE( &
        ccmax_level, &
        ccbfact_space_order, &
        ccbfact_time_order, &
        ccprescribe_temperature_outflow, &
        ccsolidheat_flag, &
        rz_flag, &
        ccFSI_flag, &
        ccZEYU_DCA_SELECT, &
        ccinvert_solid_levelset, &
        ccdenfact, &
        ccvelfact, &
        ccn_sites, &
        ccnucleation_init_time, &
        ccpos_sites, &
        ccxblob,ccyblob,cczblob,ccradblob, &
        ccxblob2,ccyblob2,cczblob2,ccradblob2, &
        ccxblob3,ccyblob3,cczblob3,ccradblob3, &
        ccxblob4,ccyblob4,cczblob4,ccradblob4, &
        ccxblob5,ccyblob5,cczblob5,ccradblob5, &
        ccxblob6,ccyblob6,cczblob6,ccradblob6, &
        ccxblob7,ccyblob7,cczblob7,ccradblob7, &
        ccxblob8,ccyblob8,cczblob8,ccradblob8, &
        ccxblob9,ccyblob9,cczblob9,ccradblob9, &
        ccxblob10,ccyblob10,cczblob10,ccradblob10, &
        ccxactive,ccyactive,cczactive, &
        ccractivex, &
        ccractivey, &
        ccractivez, &
        ccprobtype, &
        ccadv_dir, &
        ccadv_vel, &
        ccaxis_dir, &
        ccrgasinlet, &
        ccvinletgas, &
        cctwall, &
        ccadvbot, &
        ccinflow_pressure, &
        ccoutflow_pressure, &
        ccperiod_time, &
        ccproblox,ccprobloy,ccprobloz, &
        ccprobhix,ccprobhiy,ccprobhiz, &
        ccnum_species_var, &
        ccnum_materials_viscoelastic, &
        ccnum_state_material, &
        ccnum_state_base, &
        ccngeom_raw, &
        ccngeom_recon, &
        ccnum_materials, &
        ccnum_materials_vel, &
        ccnum_materials_scalar_solve, &
        ccmaterial_type, &
        ccnten, &
        ccdrhodt, &
        ccdrhodz, &
        cctempconst, &
        ccinitial_temperature, &
        cctempcutoff, &
        cctempcutoffmax, &
        ccstiffPINF, &
        ccstiffCP, &
        ccstiffCV, &
        ccstiffGAMMA, &
        ccdenconst, &
        ccden_floor, &
        ccden_ceiling, &
        cccavdenconst, &
        ccviscconst, &
        ccviscconst_eddy, &
        ccviscosity_state_model, &
        ccelastic_viscosity, &
        ccelastic_time, &
        ccviscoelastic_model, &
        cclame_coefficient, &
        cclinear_elastic_model, &
        ccshear_modulus, &
        ccstore_elastic_data, &
        ccheatviscconst, &
        ccprerecalesce_heatviscconst, &
        ccprerecalesce_viscconst, &
        ccprerecalesce_stiffCP, &
        ccprerecalesce_stiffCV, &
        ccspeciesconst, &
        ccspeciesviscconst, &
        cclatent_heat, &
        ccsaturation_temp, &
        ccmolar_mass, &
        ccspecies_molar_mass, &
        cctension, &
        cctension_slope, &
        cctension_T0, &
        cctension_min, &
        ccprefreeze_tension, &
        ccMUSHY_THICK, &
        ccgravity, &
        ccgravity_dir, &
        ccinvert_gravity, &
        ccstop_time, &
        ioproc)
      use LegendreNodes
      use probf90_module
      use global_utility_module
      use geometry_intersect_module
      use hydrateReactor_module
      use unimaterialChannel_module
      use shockdrop
      use USERDEF_module
      use CAV3D_module
      use HELIX_module
      use TSPRAY_module
      use CAV2Dstep_module
      use ZEYU_droplet_impact_module
      use STUB_module
      use GENERAL_PHASE_CHANGE_module
      use ICE_ON_SUBSTRATE_module
      use SIMPLE_PALMORE_DESJARDINS_module
      use SIMPLE_KASSEMI_module
      use DROP_IN_SHEAR_module
      use MITSUHIRO_MELTING_module
      use CRYOGENIC_TANK1_module
      use CRYOGENIC_TANK2_module
      use CRYOGENIC_TANK_MK_module
      use flexible_plate_impact_module
      use CONE3D_module
      use WAVY_Channel_module
      use rigid_FSI_module
      
      IMPLICIT NONE
      
      INTEGER_T, intent(in) :: ccmax_level
      INTEGER_T, intent(in) :: ccbfact_space_order(0:ccmax_level)
      INTEGER_T, intent(in) :: ccbfact_time_order
      INTEGER_T, intent(in) :: ccnum_materials
      
      INTEGER_T, intent(in) :: ccnten
      REAL_T, intent(in) :: ccMUSHY_THICK
      REAL_T, intent(in) :: ccgravity
      INTEGER_T, intent(in) :: ccgravity_dir
      INTEGER_T, intent(in) :: ccinvert_gravity
      INTEGER_T, intent(in) :: ccFSI_flag(ccnum_materials)
      INTEGER_T, intent(in) :: ccZEYU_DCA_SELECT
      INTEGER_T, intent(in) :: ccinvert_solid_levelset
      INTEGER_T, intent(in) :: ccprescribe_temperature_outflow
      INTEGER_T, intent(in) :: ccsolidheat_flag
      INTEGER_T, intent(in) :: rz_flag,ioproc
      INTEGER_T, intent(in) :: ccprobtype,ccadv_dir,ccaxis_dir
      
      REAL_T, intent(in) :: ccdenfact,ccvelfact
      REAL_T, intent(in) :: ccxblob,ccyblob,cczblob,ccradblob
      REAL_T, intent(in) :: ccxblob2,ccyblob2,cczblob2,ccradblob2
      REAL_T, intent(in) :: ccxblob3,ccyblob3,cczblob3,ccradblob3
      REAL_T, intent(in) :: ccxblob4,ccyblob4,cczblob4,ccradblob4
      REAL_T, intent(in) :: ccxblob5,ccyblob5,cczblob5,ccradblob5
      REAL_T, intent(in) :: ccxblob6,ccyblob6,cczblob6,ccradblob6
      REAL_T, intent(in) :: ccxblob7,ccyblob7,cczblob7,ccradblob7
      REAL_T, intent(in) :: ccxblob8,ccyblob8,cczblob8,ccradblob8
      REAL_T, intent(in) :: ccxblob9,ccyblob9,cczblob9,ccradblob9
      REAL_T, intent(in) :: ccxblob10,ccyblob10,cczblob10,ccradblob10
      REAL_T, intent(in) :: ccxactive,ccyactive,cczactive
      REAL_T, intent(in) :: ccractivex,ccractivey,ccractivez
      REAL_T, intent(in) :: ccadv_vel,ccrgasinlet
      REAL_T, intent(in) :: ccvinletgas,cctwall
      REAL_T, intent(in) :: ccadvbot
      REAL_T, intent(in) :: ccinflow_pressure
      REAL_T, intent(in) :: ccoutflow_pressure
      REAL_T, intent(in) :: ccperiod_time
      REAL_T, intent(in) :: ccproblox,ccprobloy,ccprobloz
      REAL_T, intent(in) :: ccprobhix,ccprobhiy,ccprobhiz
      REAL_T, intent(in) :: ccstop_time
      
      INTEGER_T, intent(in) :: ccnum_species_var
      
      INTEGER_T, intent(in) :: ccnum_materials_viscoelastic
      INTEGER_T :: nelastic
      
      INTEGER_T, intent(in) :: ccnum_state_material
      INTEGER_T, intent(in) :: ccnum_state_base
      INTEGER_T, intent(in) :: ccngeom_raw
      INTEGER_T, intent(in) :: ccngeom_recon
      INTEGER_T, intent(in) :: ccnum_materials_vel
      INTEGER_T, intent(in) :: ccnum_materials_scalar_solve
      
      INTEGER_T, intent(in) :: ccmaterial_type(ccnum_materials)
      REAL_T, intent(in) :: ccdrhodt(ccnum_materials)
      REAL_T, intent(in) :: ccdrhodz(ccnum_materials)
      REAL_T, intent(in) :: cctempconst(ccnum_materials)
      REAL_T, intent(in) :: ccinitial_temperature(ccnum_materials)
      REAL_T, intent(in) :: cctempcutoff(ccnum_materials)
      REAL_T, intent(in) :: cctempcutoffmax(ccnum_materials)
      REAL_T, intent(in) :: ccstiffPINF(ccnum_materials)
      REAL_T, intent(in) :: ccstiffCP(ccnum_materials)
      REAL_T, intent(in) :: ccstiffCV(ccnum_materials)
      REAL_T, intent(in) :: ccstiffGAMMA(ccnum_materials)
      REAL_T, intent(in) :: ccdenconst(ccnum_materials)
      REAL_T, intent(in) :: ccden_floor(ccnum_materials)
      REAL_T, intent(in) :: ccden_ceiling(ccnum_materials)
      REAL_T, intent(in) :: cccavdenconst(ccnum_materials)
      REAL_T, intent(in) :: ccviscconst(ccnum_materials)
      REAL_T, intent(in) :: ccviscconst_eddy(ccnum_materials)
      INTEGER_T, intent(in) :: ccviscosity_state_model(ccnum_materials)
      REAL_T, intent(in) :: ccelastic_viscosity(ccnum_materials)
      REAL_T, intent(in) :: ccelastic_time(ccnum_materials)
      INTEGER_T, intent(in) :: ccviscoelastic_model(ccnum_materials)
      REAL_T, intent(in) :: cclame_coefficient(ccnum_materials)
      INTEGER_T, intent(in) :: cclinear_elastic_model(ccnum_materials)
      REAL_T, intent(in) :: ccshear_modulus(ccnum_materials)
      INTEGER_T, intent(in) :: ccstore_elastic_data(ccnum_materials)
      REAL_T, intent(in) :: ccheatviscconst(ccnum_materials)
      REAL_T, intent(in) :: ccprerecalesce_heatviscconst(ccnum_materials)
      REAL_T, intent(in) :: ccprerecalesce_viscconst(ccnum_materials)
      REAL_T, intent(in) :: ccprerecalesce_stiffCP(ccnum_materials)
      REAL_T, intent(in) :: ccprerecalesce_stiffCV(ccnum_materials)
      REAL_T, intent(in) :: &
        ccspeciesconst((ccnum_species_var+1)*ccnum_materials)
      REAL_T, intent(in) :: &
        ccspeciesviscconst((ccnum_species_var+1)*ccnum_materials)
      REAL_T, intent(in) :: cclatent_heat(2*ccnten)
      REAL_T, intent(in) :: ccsaturation_temp(2*ccnten)
      REAL_T, intent(in) :: ccmolar_mass(ccnum_materials)
      REAL_T, intent(in) :: ccspecies_molar_mass(ccnum_species_var+1)
      REAL_T, intent(in) :: cctension(ccnten)
      REAL_T, intent(in) :: cctension_slope(ccnten)
      REAL_T, intent(in) :: cctension_T0(ccnten)
      REAL_T, intent(in) :: cctension_min(ccnten)
      REAL_T, intent(in) :: ccprefreeze_tension(ccnten)
      
      INTEGER_T, intent(in) :: ccn_sites
      REAL_T, intent(in) :: ccnucleation_init_time
      REAL_T, intent(in) :: ccpos_sites(400)
      
      character*12 namestr1
      character*13 namestr2
      INTEGER_T i
      
      INTEGER_T im,iten
      INTEGER_T nten
      INTEGER_T level,bfactmax
      REAL_T :: massfrac_parm(ccnum_species_var+1)
      
      probtype=ccprobtype
      num_materials=ccnum_materials
      if (num_materials.lt.MAX_NUM_MATERIALS) then
       ! do nothing
      else
       print *,"increase MAX_NUM_MATERIALS; aborting"
       stop
      endif
      
      ! USER DEFINED (used by "is_in_probtype_list")
      ! IN ORDER TO ADD A NEW TEST PROBLEM:
      ! 1. increment probtype_list_size by 1.
      ! 2. set used_probtypes( <new value for probtype_list_size> ) = 
      !    new_probtype
      ! 3. a) link template subroutine names to new module names 
      !       (see probtype.eq.2002 below for an example)
      !    b) add appropriate "use" command to the beginning of this routine. 
      ! 4. create new module file (e.g. by copying an existing module file)
      ! 5. update Make.package accordingly (2 places)
      ! 6. create inputs file
      probtype_list_size=10
      used_probtypes(1)=2000 ! flexible_plate_impact
      used_probtypes(2)=421  ! CRYOGENIC_TANK1
      used_probtypes(3)=414  ! MITSUHIRO_MELTING
      used_probtypes(4)=2001 ! ICE_ON_SUBSTRATE
      used_probtypes(5)=2002 ! 1D TEST FROM PALMORE and Desjardins
      used_probtypes(6)=55   ! GENERAL_PHASE_CHANGE
      used_probtypes(7)=422  ! CRYOGENIC_TANK2
      used_probtypes(8)=423  ! CRYOGENIC_TANK_MK
      used_probtypes(9)=424  ! DROP_IN_SHEAR
      used_probtypes(10)=2003 ! 1D TEST, Kassemi model for PD test 
      
      SUB_INIT_MODULE=>INIT_STUB_MODULE
      SUB_LS=>STUB_LS
      SUB_VEL=>STUB_VEL
      SUB_EOS=>EOS_STUB
      SUB_SOUNDSQR=>SOUNDSQR_STUB
      SUB_INTERNAL=>INTERNAL_STUB
      SUB_TEMPERATURE=>TEMPERATURE_STUB
      SUB_PRES=>STUB_PRES
      SUB_STATE=>STUB_STATE
      SUB_LS_BC=>STUB_LS_BC
      SUB_VEL_BC=>STUB_VEL_BC
      SUB_PRES_BC=>STUB_PRES_BC
      SUB_STATE_BC=>STUB_STATE_BC
      SUB_HEATSOURCE=>STUB_HEATSOURCE
      SUB_EB_heat_source=>STUB_EB_heat_source
      SUB_microcell_heat_coeff=>STUB_microcell_heat_coeff
      SUB_velfreestream=>STUB_velfreestream
      SUB_nucleation=>STUB_nucleation
      SUB_CFL_HELPER=>STUB_CFL_HELPER
      SUB_hydro_pressure_density=>STUB_hydro_pressure_density
      SUB_correct_pres_rho_hydrostatic=>STUB_correct_pres_rho_hydrostatic
      SUB_ASSIMILATE=>STUB_ASSIMILATE
      SUB_SUMINT=>STUB_SUMINT
      SUB_clamped_LS_no_scale=>STUB_CLAMPED_LS
      
      if (probtype.eq.421) then
       SUB_INIT_MODULE=>INIT_CRYOGENIC_TANK1_MODULE
       SUB_LS=>CRYOGENIC_TANK1_LS
       SUB_VEL=>CRYOGENIC_TANK1_VEL
       SUB_EOS=>EOS_CRYOGENIC_TANK1
       SUB_SOUNDSQR=>SOUNDSQR_CRYOGENIC_TANK1
       SUB_INTERNAL=>INTERNAL_CRYOGENIC_TANK1
       SUB_TEMPERATURE=>TEMPERATURE_CRYOGENIC_TANK1
       SUB_PRES=>CRYOGENIC_TANK1_PRES
       SUB_STATE=>CRYOGENIC_TANK1_STATE
       SUB_LS_BC=>CRYOGENIC_TANK1_LS_BC
       SUB_VEL_BC=>CRYOGENIC_TANK1_VEL_BC
       SUB_PRES_BC=>CRYOGENIC_TANK1_PRES_BC
       SUB_STATE_BC=>CRYOGENIC_TANK1_STATE_BC
       SUB_HEATSOURCE=>CRYOGENIC_TANK1_HEATSOURCE
      
      else if (probtype.eq.422) then
       SUB_INIT_MODULE=>INIT_CRYOGENIC_TANK2_MODULE
       SUB_LS=>CRYOGENIC_TANK2_LS
       SUB_VEL=>CRYOGENIC_TANK2_VEL
       SUB_EOS=>EOS_CRYOGENIC_TANK2
       SUB_SOUNDSQR=>SOUNDSQR_CRYOGENIC_TANK2
       SUB_INTERNAL=>INTERNAL_CRYOGENIC_TANK2
       SUB_TEMPERATURE=>TEMPERATURE_CRYOGENIC_TANK2
       SUB_PRES=>CRYOGENIC_TANK2_PRES
       SUB_STATE=>CRYOGENIC_TANK2_STATE
       SUB_LS_BC=>CRYOGENIC_TANK2_LS_BC
       SUB_VEL_BC=>CRYOGENIC_TANK2_VEL_BC
       SUB_PRES_BC=>CRYOGENIC_TANK2_PRES_BC
       SUB_STATE_BC=>CRYOGENIC_TANK2_STATE_BC
       SUB_HEATSOURCE=>CRYOGENIC_TANK2_HEATSOURCE
      
      else if (probtype.eq.423) then

       SUB_INIT_MODULE=>INIT_CRYOGENIC_TANK_MK_MODULE
       SUB_LS=>CRYOGENIC_TANK_MK_LS
       SUB_VEL=>CRYOGENIC_TANK_MK_VEL
       SUB_EOS=>EOS_CRYOGENIC_TANK_MK
       SUB_SOUNDSQR=>SOUNDSQR_CRYOGENIC_TANK_MK
       SUB_INTERNAL=>INTERNAL_CRYOGENIC_TANK_MK
       SUB_TEMPERATURE=>TEMPERATURE_CRYOGENIC_TANK_MK
       SUB_PRES=>CRYOGENIC_TANK_MK_PRES
       SUB_STATE=>CRYOGENIC_TANK_MK_STATE
       SUB_LS_BC=>CRYOGENIC_TANK_MK_LS_BC
       SUB_VEL_BC=>CRYOGENIC_TANK_MK_VEL_BC
       SUB_PRES_BC=>CRYOGENIC_TANK_MK_PRES_BC
       SUB_STATE_BC=>CRYOGENIC_TANK_MK_STATE_BC
       SUB_correct_pres_rho_hydrostatic=> &
              CRYOGENIC_TANK_MK_correct_pres_rho_hydrostatic
       SUB_HEATSOURCE=>CRYOGENIC_TANK_MK_HEATSOURCE

      else if (probtype.eq.424) then

       SUB_INIT_MODULE=>INIT_DROP_IN_SHEAR_MODULE
       SUB_CFL_HELPER=>DROP_IN_SHEAR_CFL_HELPER
       SUB_LS=>DROP_IN_SHEAR_LS
       SUB_VEL=>DROP_IN_SHEAR_VEL
       SUB_PRES=>DROP_IN_SHEAR_PRES
       SUB_STATE=>DROP_IN_SHEAR_STATE
       SUB_LS_BC=>DROP_IN_SHEAR_LS_BC
       SUB_VEL_BC=>DROP_IN_SHEAR_VEL_BC
       SUB_PRES_BC=>DROP_IN_SHEAR_PRES_BC
       SUB_STATE_BC=>DROP_IN_SHEAR_STATE_BC
       SUB_HEATSOURCE=>DROP_IN_SHEAR_HEATSOURCE

      else if (probtype.eq.414) then
       SUB_INIT_MODULE=>INIT_MITSUHIRO_MELTING_MODULE
       SUB_LS=>MITSUHIRO_MELTING_LS
       SUB_VEL=>MITSUHIRO_MELTING_VEL
       SUB_PRES=>MITSUHIRO_MELTING_PRES
       SUB_STATE=>MITSUHIRO_MELTING_STATE
       SUB_LS_BC=>MITSUHIRO_MELTING_LS_BC
       SUB_VEL_BC=>MITSUHIRO_MELTING_VEL_BC
       SUB_PRES_BC=>MITSUHIRO_MELTING_PRES_BC
       SUB_STATE_BC=>MITSUHIRO_MELTING_STATE_BC
       SUB_HEATSOURCE=>MITSUHIRO_MELTING_HEATSOURCE
      else if (probtype.eq.2001) then
       SUB_INIT_MODULE=>INIT_ICE_ON_SUBSTRATE_MODULE
       SUB_LS=>ICE_ON_SUBSTRATE_LS
       SUB_VEL=>ICE_ON_SUBSTRATE_VEL
       SUB_PRES=>ICE_ON_SUBSTRATE_PRES
       SUB_STATE=>ICE_ON_SUBSTRATE_STATE
       SUB_LS_BC=>ICE_ON_SUBSTRATE_LS_BC
       SUB_VEL_BC=>ICE_ON_SUBSTRATE_VEL_BC
       SUB_PRES_BC=>ICE_ON_SUBSTRATE_PRES_BC
       SUB_STATE_BC=>ICE_ON_SUBSTRATE_STATE_BC
       SUB_HEATSOURCE=>ICE_ON_SUBSTRATE_HEATSOURCE
      else if (probtype.eq.2002) then
       SUB_INIT_MODULE=>INIT_SIMPLE_PALMORE_DESJARDINS_MODULE
       SUB_LS=>SIMPLE_PALMORE_DESJARDINS_LS
       SUB_VEL=>SIMPLE_PALMORE_DESJARDINS_VEL
       SUB_PRES=>SIMPLE_PALMORE_DESJARDINS_PRES
       SUB_STATE=>SIMPLE_PALMORE_DESJARDINS_STATE
       SUB_LS_BC=>SIMPLE_PALMORE_DESJARDINS_LS_BC
       SUB_VEL_BC=>SIMPLE_PALMORE_DESJARDINS_VEL_BC
       SUB_PRES_BC=>SIMPLE_PALMORE_DESJARDINS_PRES_BC
       SUB_STATE_BC=>SIMPLE_PALMORE_DESJARDINS_STATE_BC
       SUB_HEATSOURCE=>SIMPLE_PALMORE_DESJARDINS_HEATSOURCE
       SUB_SUMINT=>SIMPLE_PALMORE_DESJARDINS_SUMINT ! compare with analytical
       SUB_ASSIMILATE=>SIMPLE_PALMORE_DESJARDINS_ASSIMILATE
      else if (probtype.eq.2003) then
       SUB_INIT_MODULE=>INIT_SIMPLE_KASSEMI_MODULE
       SUB_LS=>SIMPLE_KASSEMI_LS
       SUB_VEL=>SIMPLE_KASSEMI_VEL
       SUB_PRES=>SIMPLE_KASSEMI_PRES
       SUB_STATE=>SIMPLE_KASSEMI_STATE
       SUB_LS_BC=>SIMPLE_KASSEMI_LS_BC
       SUB_VEL_BC=>SIMPLE_KASSEMI_VEL_BC

       SUB_EOS=>EOS_KASSEMI_MK
       SUB_SOUNDSQR=>SOUNDSQR_KASSEMI_MK
       SUB_INTERNAL=>INTERNAL_KASSEMI_MK
       SUB_TEMPERATURE=>TEMPERATURE_KASSEMI_MK

       SUB_PRES_BC=>SIMPLE_KASSEMI_PRES_BC
       SUB_STATE_BC=>SIMPLE_KASSEMI_STATE_BC
       SUB_HEATSOURCE=>SIMPLE_KASSEMI_HEATSOURCE
       SUB_SUMINT=>SIMPLE_KASSEMI_SUMINT ! compare with analytical
       SUB_ASSIMILATE=>SIMPLE_KASSEMI_ASSIMILATE
      else if (probtype.eq.2000) then
       SUB_INIT_MODULE=>INIT_flexible_plate_impact_MODULE
       SUB_clamped_LS_no_scale=>flexible_plate_clamped_LS
       SUB_LS=>flexible_plate_impact_LS
       SUB_VEL=>flexible_plate_impact_VEL
       SUB_PRES=>flexible_plate_impact_PRES
       SUB_STATE=>flexible_plate_impact_STATE
       SUB_LS_BC=>flexible_plate_impact_LS_BC
       SUB_VEL_BC=>flexible_plate_impact_VEL_BC
       SUB_PRES_BC=>flexible_plate_impact_PRES_BC
       SUB_STATE_BC=>flexible_plate_impact_STATE_BC
       SUB_HEATSOURCE=>flexible_plate_impact_HEATSOURCE
      else if (probtype.eq.55) then
       SUB_INIT_MODULE=>INIT_GENERAL_PHASE_CHANGE_MODULE
       SUB_LS=>GENERAL_PHASE_CHANGE_LS
       SUB_VEL=>GENERAL_PHASE_CHANGE_VEL
       SUB_PRES=>GENERAL_PHASE_CHANGE_PRES
       SUB_STATE=>GENERAL_PHASE_CHANGE_STATE
       SUB_LS_BC=>GENERAL_PHASE_CHANGE_LS_BC
       SUB_VEL_BC=>GENERAL_PHASE_CHANGE_VEL_BC
       SUB_PRES_BC=>GENERAL_PHASE_CHANGE_PRES_BC
       SUB_STATE_BC=>GENERAL_PHASE_CHANGE_STATE_BC
       SUB_HEATSOURCE=>GENERAL_PHASE_CHANGE_HEATSOURCE
       SUB_EB_heat_source=>GENERAL_PHASE_CHANGE_EB_heat_source
       SUB_microcell_heat_coeff=>GENERAL_PHASE_CHANGE_microcell_heat_coeff
       SUB_velfreestream=>GENERAL_PHASE_CHANGE_velfreestream
       SUB_nucleation=>GENERAL_PHASE_CHANGE_nucleation
       SUB_CFL_HELPER=>GENERAL_PHASE_CHANGE_CFL_HELPER
       SUB_SUMINT=>GENERAL_PHASE_CHANGE_SUMINT ! Nusseltt number
       SUB_hydro_pressure_density=>GENERAL_PHASE_CHANGE_hydro_pressure_density
      else
       ! assign null routines here that would cause the program to abort
       ! if called.  In otherwords, these are routines THAT MUST BE DEFINED
       ! and CANNOT depend on the STUB routines.
       SUB_INIT_MODULE=>NULL()
       SUB_LS=>NULL()
!       SUB_clamped_LS_no_scale=>NULL()
       SUB_VEL=>NULL()
       SUB_EOS=>NULL()
       SUB_SOUNDSQR=>NULL()
       SUB_INTERNAL=>NULL()
       SUB_TEMPERATURE=>NULL()
       SUB_PRES=>NULL()
       SUB_STATE=>NULL()
       SUB_LS_BC=>NULL()
       SUB_VEL_BC=>NULL()
       SUB_PRES_BC=>NULL()
       SUB_STATE_BC=>NULL()
       SUB_HEATSOURCE=>NULL()
       SUB_EB_heat_source=>NULL()
       SUB_microcell_heat_coeff=>NULL()
       SUB_velfreestream=>NULL()
       SUB_nucleation=>NULL()
       SUB_CFL_HELPER=>NULL()
       SUB_hydro_pressure_density=>NULL()
       SUB_correct_pres_rho_hydrostatic=>NULL()
      endif
      
      global_pressure_scale=one
      global_velocity_scale=one
      
      fort_max_level=ccmax_level
      do level=0,fort_max_level
       bfact_space_order(level)=ccbfact_space_order(level)
      enddo
      bfact_time_order=ccbfact_time_order
      
      if (fort_max_level.lt.0) then
       print *,"fort_max_level invalid"
       stop
      endif
      do level=0,fort_max_level
       if (bfact_space_order(level).lt.1) then
        print *,"bfact_space_order invalid: ",bfact_space_order(level)
        stop
       endif
      enddo
      if (bfact_time_order.lt.1) then
       print *,"bfact_time_order invalid"
       stop
      endif
      bfactmax=16
      call sanity_check(bfactmax+2)
      call init_cache(bfactmax+2)
      
      if (ioproc.eq.1) then
       print *,"LEGEND FOR SPECTRAL ELEMENT NODES"
       print *,"0 Legendre Nodes"
       print *,"1 Clenshaw Curtis Nodes (Chebychev)"
       print *,"SPTYPE= ",SPTYPE
       print *,"TMTYPE= ",TMTYPE
      else if (ioproc.eq.0) then
       ! do nothing
      else
       print *,"ioproc invalid"
       stop
      endif
       
      problox=ccproblox
      probloy=ccprobloy
      probloz=ccprobloz
      probhix=ccprobhix
      probhiy=ccprobhiy
      probhiz=ccprobhiz
      
      if (SDIM.eq.2) then
       probhiz=probhiy
       probloz=probloy
      endif
      
      problenx=probhix-problox
      probleny=probhiy-probloy
      problenz=probhiz-probloz
      
      if ((problenx.le.zero).or. &
          (probleny.le.zero).or. &
          (problenz.le.zero)) then
       print *,"problenx or probleny or problenz invalid"
       stop
      endif
      
      fort_stop_time=ccstop_time
      
      prescribe_temperature_outflow= &
        ccprescribe_temperature_outflow
      
      fort_solidheat_flag=ccsolidheat_flag
      
      levelrz=rz_flag
      denfact=ccdenfact
      velfact=ccvelfact
      
      n_sites=ccn_sites
      nucleation_init_time=ccnucleation_init_time
      
      if (nucleation_init_time.lt.zero) then
       print *,"nucleation_init_time.lt.zero"
       stop
      endif
      
      if ((n_sites.lt.0).or.(4*n_sites.gt.400)) then
       print *,"n_sites invalid"
       stop
      endif
      do i=1,4*n_sites
       pos_sites(i)=ccpos_sites(i)
      enddo
      
      xblob=ccxblob
      yblob=ccyblob
      zblob=cczblob
      radblob=ccradblob
      
      xblob2=ccxblob2
      yblob2=ccyblob2
      zblob2=cczblob2
      radblob2=ccradblob2
      
      xblob3=ccxblob3
      yblob3=ccyblob3
      zblob3=cczblob3
      radblob3=ccradblob3
      
      xblob4=ccxblob4
      yblob4=ccyblob4
      zblob4=cczblob4
      radblob4=ccradblob4
      
      xblob5=ccxblob5
      yblob5=ccyblob5
      zblob5=cczblob5
      radblob5=ccradblob5
      
      xblob6=ccxblob6
      yblob6=ccyblob6
      zblob6=cczblob6
      radblob6=ccradblob6
      
      xblob7=ccxblob7
      yblob7=ccyblob7
      zblob7=cczblob7
      radblob7=ccradblob7
      
      xblob8=ccxblob8
      yblob8=ccyblob8
      zblob8=cczblob8
      radblob8=ccradblob8
      
      xblob9=ccxblob9
      yblob9=ccyblob9
      zblob9=cczblob9
      radblob9=ccradblob9
      
      xblob10=ccxblob10
      yblob10=ccyblob10
      zblob10=cczblob10
      radblob10=ccradblob10
      
      xblobarr(1)=xblob 
      yblobarr(1)=yblob 
      zblobarr(1)=zblob 
      radblobarr(1)=radblob 
      
      xblobarr(2)=xblob2 
      yblobarr(2)=yblob2 
      zblobarr(2)=zblob2
      radblobarr(2)=radblob2 
      
      xblobarr(3)=xblob3 
      yblobarr(3)=yblob3 
      zblobarr(3)=zblob3
      radblobarr(3)=radblob3 
      
      xblobarr(4)=xblob4 
      yblobarr(4)=yblob4 
      zblobarr(4)=zblob4
      radblobarr(4)=radblob4 
      
      xblobarr(5)=xblob5 
      yblobarr(5)=yblob5 
      zblobarr(5)=zblob5
      radblobarr(5)=radblob5 
      
      xblobarr(6)=xblob6 
      yblobarr(6)=yblob6 
      zblobarr(6)=zblob6
      radblobarr(6)=radblob6 
      
      xblobarr(7)=xblob7 
      yblobarr(7)=yblob7 
      zblobarr(7)=zblob7
      radblobarr(7)=radblob7 
      
      xblobarr(8)=xblob8 
      yblobarr(8)=yblob8 
      zblobarr(8)=zblob8
      radblobarr(8)=radblob8 
      
      xblobarr(9)=xblob9 
      yblobarr(9)=yblob9 
      zblobarr(9)=zblob9
      radblobarr(9)=radblob9 
      
      xblobarr(10)=xblob10 
      yblobarr(10)=yblob10 
      zblobarr(10)=zblob10
      radblobarr(10)=radblob10 
      
      xactive=ccxactive
      yactive=ccyactive
      zactive=cczactive
      ractivex=ccractivex
      ractivey=ccractivey
      ractivez=ccractivez
      
      adv_dir=ccadv_dir
      adv_vel=ccadv_vel
      axis_dir=ccaxis_dir
      rgasinlet=ccrgasinlet
      vinletgas=ccvinletgas
      twall=cctwall
      advbot=ccadvbot
      inflow_pressure=ccinflow_pressure
      outflow_pressure=ccoutflow_pressure
      period_time=ccperiod_time
      
      do im=1,num_materials
      
       fort_material_type(im)=ccmaterial_type(im)
       FSI_flag(im)=ccFSI_flag(im)
       if (fort_material_type(im).eq.0) then
        ! do nothing
       else if (fort_material_type(im).eq.999) then
        if ((FSI_flag(im).ne.1).and. &
            (FSI_flag(im).ne.2).and. &
            (FSI_flag(im).ne.4)) then
         print *,"FSI_flag invalid"
         stop
        endif
       else if ((fort_material_type(im).gt.0).and. &
                (fort_material_type(im).le.MAX_NUM_EOS)) then
        ! do nothing
       else
        print *,"fort_material_type invalid"
        stop
       endif
      
      enddo ! im=1..num_materials
      
      fort_ZEYU_DCA_SELECT=ccZEYU_DCA_SELECT
      
      invert_solid_levelset=ccinvert_solid_levelset
      if ((invert_solid_levelset.ne.0).and. &
          (invert_solid_levelset.ne.1)) then
       print *,"invert_solid_levelset invalid in override"
       stop
      endif
      
      num_species_var=ccnum_species_var
      if (num_species_var.lt.MAX_NUM_SPECIES) then
       ! do nothing
      else
       print *,"num_species_var too large, increase MAX_NUM_SPECIES"
       stop
      endif
      
      num_materials_viscoelastic=ccnum_materials_viscoelastic
      num_materials_vel=ccnum_materials_vel
      num_materials_scalar_solve=ccnum_materials_scalar_solve
      
      num_state_base=ccnum_state_base
      if (num_state_base.ne.2) then
       print *,"num_state_base invalid 60"
       stop
      endif
      num_state_material=num_state_base  ! den,T
      num_state_material=num_state_material+num_species_var
      
      if (num_state_material.ne.ccnum_state_material) then
       print *,"ccnum_state_material invalid"
       stop
      endif
      
      if ((num_species_var.lt.0).or. &
          (num_materials_viscoelastic.lt.0).or. &
          (num_materials_viscoelastic.gt.num_materials).or. &
          (num_materials.lt.1).or. &
          (num_materials.gt.MAX_NUM_MATERIALS).or. &
          (num_materials.gt.100).or. &
          (num_materials_vel.ne.1).or. &
          ((num_materials_scalar_solve.ne.1).and. &
           (num_materials_scalar_solve.ne.num_materials)) ) then
       print *,"material parameters illegal"
       stop
      endif
      
      ngeom_raw=ccngeom_raw
      ngeom_recon=ccngeom_recon
      
      if (ngeom_recon.ne.2*SDIM+3) then
       print *,"ngeom_recon invalid"
       stop
      endif
      if (ngeom_raw.ne.SDIM+1) then
       print *,"ngeom_raw invalid"
       stop
      endif
      
      if (ioproc.eq.1) then
      
       print *,"fort_max_level= ",fort_max_level
       do level=0,fort_max_level
        print *,"level,bfact_space_order ",level,bfact_space_order(level)
       enddo
       print *,"bfact_time_order ",bfact_time_order
      
       print *,"fort material parameters"
       print *,"numspec,num_mat_visc,MAX_NUM_MATERIALS,num_materials ", &
        num_species_var,num_materials_viscoelastic, &
        MAX_NUM_MATERIALS,num_materials
       print *,"num_materials_vel ",num_materials_vel
       print *,"num_materials_scalar_solve ",num_materials_scalar_solve
       print *,"MAX_NUM_EOS ",MAX_NUM_EOS
       print *,"ngeom_raw ",ngeom_raw
       print *,"ngeom_recon ",ngeom_recon
       print *,"fort: num_state_material ",num_state_material
       print *,"fort: num_state_base ",num_state_base
      
      else if (ioproc.eq.0) then
       ! do nothing
      else
       print *,"ioproc invalid"
       stop
      endif
      
      do im=1,num_species_var*num_materials
       fort_speciesconst(im)=ccspeciesconst(im)
       fort_speciesviscconst(im)=ccspeciesviscconst(im)
      enddo
      
      do im=1,num_materials
      
       fort_drhodt(im)=ccdrhodt(im)
       fort_drhodz(im)=ccdrhodz(im)
       fort_tempconst(im)=cctempconst(im)
       fort_initial_temperature(im)=ccinitial_temperature(im)
       fort_tempcutoff(im)=cctempcutoff(im) ! default 1.0E-8
       fort_tempcutoffmax(im)=cctempcutoffmax(im) ! default 1.0D+99
       fort_stiffPINF(im)=ccstiffPINF(im)
       fort_stiffCP(im)=ccstiffCP(im)
       fort_stiffCV(im)=ccstiffCV(im)
       fort_stiffGAMMA(im)=ccstiffGAMMA(im)
       fort_denconst(im)=ccdenconst(im)
       fort_density_floor(im)=ccden_floor(im)
       fort_density_ceiling(im)=ccden_ceiling(im)
       fort_cavdenconst(im)=cccavdenconst(im)
       fort_viscconst(im)=ccviscconst(im)
       fort_viscconst_eddy(im)=ccviscconst_eddy(im)
       fort_viscosity_state_model(im)= &
         ccviscosity_state_model(im)
       fort_elastic_viscosity(im)=ccelastic_viscosity(im)
       fort_elastic_time(im)=ccelastic_time(im)
       fort_viscoelastic_model(im)=ccviscoelastic_model(im)
       fort_lame_coefficient(im)=cclame_coefficient(im)
       fort_linear_elastic_model(im)=cclinear_elastic_model(im)
       fort_shear_modulus(im)=ccshear_modulus(im)
       fort_store_elastic_data(im)=ccstore_elastic_data(im)
       fort_heatviscconst(im)=ccheatviscconst(im)
       fort_prerecalesce_heatviscconst(im)=ccprerecalesce_heatviscconst(im)
       fort_prerecalesce_viscconst(im)=ccprerecalesce_viscconst(im)
       fort_prerecalesce_stiffCP(im)=ccprerecalesce_stiffCP(im)
       fort_prerecalesce_stiffCV(im)=ccprerecalesce_stiffCV(im)
      
       fort_im_elastic_map(im)=-1
      
      enddo ! im=1..num_materials
      
      nelastic=0
      do im=1,num_materials
       if (fort_store_elastic_data(im).eq.0) then
        ! do nothing
       else if (fort_store_elastic_data(im).eq.1) then
        nelastic=nelastic+1
        fort_im_elastic_map(nelastic)=im-1
       else
        print *,"fort_store_elastic_data(im) invalid"
        stop
       endif
      enddo ! im=1..num_materials
      
      if (nelastic.ne.num_materials_viscoelastic) then
       print *,"nelastic.ne.num_materials_viscoelastic"
       print *,"nelastic ",nelastic
       print *,"num_materials_viscoelastic ",num_materials_viscoelastic
       do im=1,num_materials
        print *,"im,fort_store_elastic_data(im) ",im, &
                fort_store_elastic_data(im)
       enddo
       stop
      endif
      
      nten=( (num_materials-1)*(num_materials-1)+ &
         num_materials-1 )/2
      if (nten.ne.ccnten) then
       print *,"nten or ccnten invalid"
       print *,"nten=",nten
       print *,"ccnten=",ccnten
       stop
      endif
      if (nten.ge.MAX_NUM_INTERFACES) then
       print *,"too many surface tension coefficients, increase "
       print *,"MAX_NUM_INTERFACES"
       print *,"nten= ",nten
       stop
      endif
      do iten=1,nten
       fort_latent_heat(iten)=cclatent_heat(iten)
       fort_latent_heat(nten+iten)=cclatent_heat(nten+iten)
       fort_saturation_temp(iten)=ccsaturation_temp(iten)
       fort_saturation_temp(nten+iten)=ccsaturation_temp(nten+iten)
      
       fort_tension(iten)=cctension(iten)
       fort_tension_slope(iten)=cctension_slope(iten)
       fort_tension_T0(iten)=cctension_T0(iten)
       fort_tension_min(iten)=cctension_min(iten)
       fort_prefreeze_tension(iten)=ccprefreeze_tension(iten)
      enddo
      
      do im=1,num_materials
       fort_molar_mass(im)=ccmolar_mass(im)
       if (fort_molar_mass(im).gt.zero) then
        ! do nothing
       else
        print *,"fort_molar_mass invalid"
        stop
       endif
      enddo
      do im=1,num_species_var
       fort_species_molar_mass(im)=ccspecies_molar_mass(im)
       if (fort_species_molar_mass(im).gt.zero) then
        ! do nothing
       else
        print *,"fort_species_molar_mass invalid"
        stop
       endif
      enddo
      
      if (ioproc.eq.1) then
       print *,"n_sites= ",n_sites
       do i=1,4*n_sites
        print *,"i,pos_sites=",i,pos_sites(i)
       enddo 
       do im=1,num_materials
        print *,"im,mat type ",im,fort_material_type(im)
        print *,"im,fort_molar_mass ",im,fort_molar_mass(im)
        print *,"im,drhodt ",im,fort_drhodt(im)
        print *,"im,drhodz ",im,fort_drhodz(im)
        print *,"im,temp ",im,fort_tempconst(im)
        print *,"im,initial_temp ",im,fort_initial_temperature(im)
        print *,"im,tempcutoff ",im,fort_tempcutoff(im)
        print *,"im,tempcutoffmax ",im,fort_tempcutoffmax(im)
        print *,"im,stiffPINF ",im,fort_stiffPINF(im)
        print *,"im,stiffCP ",im,fort_stiffCP(im)
        print *,"im,stiffCV ",im,fort_stiffCV(im)
        print *,"im,stiffGAMMA ",im,fort_stiffGAMMA(im)
        print *,"im,den ",im,fort_denconst(im)
        print *,"im,den_floor ",im,fort_density_floor(im)
        print *,"im,den_ceiling ",im,fort_density_ceiling(im)
        print *,"im,cavden ",im,fort_cavdenconst(im)
        print *,"im,visc ",im,fort_viscconst(im)
        print *,"im,visc_eddy ",im,fort_viscconst_eddy(im)
        print *,"im,viscosity_state_model ",im, &
         fort_viscosity_state_model(im)
        print *,"im,fort_elastic_viscosity ",im,fort_elastic_viscosity(im)
        print *,"im,fort_elastic_time ",im,fort_elastic_time(im)
        print *,"im,fort_viscoelastic_model ",im,fort_viscoelastic_model(im)
        print *,"im,fort_lame_coefficient ",im,fort_lame_coefficient(im)
        print *,"im,fort_linear_elastic_model ",im, &
                fort_linear_elastic_model(im)
        print *,"im,fort_shear_modulus ",im,fort_shear_modulus(im)
        print *,"im,fort_store_elastic_data ",im,fort_store_elastic_data(im)
        print *,"im,fort_im_elastic_map ",im,fort_im_elastic_map(im)
        print *,"im,heatvisc ",im,fort_heatviscconst(im)
        print *,"im,prerecalesce_heatvisc ",im, &
         fort_prerecalesce_heatviscconst(im)
        print *,"im,prerecalesce_visc ",im, &
         fort_prerecalesce_viscconst(im)
        print *,"im,prerecalesce_cp ",im, &
         fort_prerecalesce_stiffCP(im)
        print *,"im,prerecalesce_cv ",im, &
         fort_prerecalesce_stiffCV(im)
       enddo ! im
       do im=1,num_species_var*num_materials
        print *,"im,species ",im,fort_speciesconst(im)
        print *,"im,speciesvisc ",im,fort_speciesviscconst(im)
       enddo
       do im=1,num_species_var
        print *,"im,fort_species_molar_mass ",im,fort_species_molar_mass(im)
       enddo
      
       do iten=1,nten
        print *,"iten,tension ",iten,fort_tension(iten)
        print *,"iten,tension_slope ",iten,fort_tension_slope(iten)
        print *,"iten,tension_T0 ",iten,fort_tension_T0(iten)
        print *,"iten,tension_min ",iten,fort_tension_min(iten)
        print *,"iten,prefreeze_tension ",iten,fort_prefreeze_tension(iten)
       enddo
      endif
      
      FORT_MUSHY_THICK=ccMUSHY_THICK
      
      gravity=ccgravity
      gravity_dir=ccgravity_dir
      invert_gravity=ccinvert_gravity
      
       ! in: GLOBALUTIL.F90
      call init_density_at_depth()
      
      pres_homflag=0
      vel_homflag=0
      temp_homflag=0
      species_homflag=0
      ls_homflag=1  ! default 90 degree contact angle.
      
      inflow_count=0
      outflow_count=0
      last_inflow_index=1
      last_outflow_index=1
      
      if (is_in_probtype_list().eq.1) then
      
       call SUB_INIT_MODULE()
      
      else if ((probtype.eq.110).and.(SDIM.eq.2)) then
      
       print *,"opening InflowBC.dat and OutflowBC.dat"
       namestr1='InflowBC.dat' 
       namestr2='OutflowBC.dat' 
      
       open(unit=11,file=namestr1)
       read(11,*) inflow_count
       print *,"inflow_count= ",inflow_count
       do i=1,inflow_count
        read(11,*) inflow_time(i),inflow_velocity(i), &
          inflow_elevation(i)
       enddo
       close(11)
      
       open(unit=12,file=namestr2)
       read(12,*) outflow_count
       print *,"outflow_count= ",outflow_count
       do i=1,outflow_count
        read(12,*) outflow_time(i),outflow_velocity(i), &
          outflow_elevation(i)
       enddo
       close(12)
      
       call shallow_water_solve()
      
       ! above: probtype==110
      else if ((probtype.eq.1).and. &
               ((axis_dir.eq.150).or. &
                (axis_dir.eq.151))) then
      
       call shockdrop_init()
      
      else if (probtype.eq.411) then
      
       call INIT_CAV3D_MODULE()
      
      else if (probtype.eq.401) then
      
       call INIT_HELIX_MODULE()
      
      else if (probtype.eq.402) then
      
       call INIT_TSPRAY_MODULE()
      else if (probtype.eq.412) then
      
       call INIT_CAV2Dstep_MODULE()
      
      else if (probtype.eq.413) then
      
       call INIT_ZEYU_droplet_impact_MODULE()
      
      else if (probtype.eq.533) then
      
       call INIT_rigid_FSI_MODULE()
      
      else if (probtype.eq.311) then
      
       call INIT_USERDEF_MODULE()
      
      else if (probtype.eq.222) then
      
       call INIT_CONE3D_MODULE()
      
      else if (probtype.eq.915) then
      
       call INIT_WAVY_MODULE()
      
      else if (probtype.eq.199) then ! hydrate
       if (num_materials.ne.3) then
        print *,"num_materials invalid probtype=199"
        stop
       endif
       if (num_species_var.ne.1) then
        print *,"num_species_var invalid"
        stop
       endif
        ! 1=water 2=methane 3=hydrate
       call INIT_HYDRATE_MODULE( &
        fort_tempconst(1),fort_tempconst(3), &
        fort_tempconst(2), &
        fort_denconst(1),fort_denconst(3), &
        fort_denconst(2), &
        outflow_pressure, &
        probhiy,probhix, &
        yblob,yblob2) 
      
       ! in: subroutine FORT_OVERRIDE
      else if (probtype.eq.220) then
       if (num_materials.ne.3) then
        print *,"num_materials invalid probtype=220"
        stop
       endif
       if (num_species_var.ne.0) then
        print *,"num_species_var invalid probtype=220"
        stop
       endif
        ! 1=Material 2=ghost material 3=fiber
       call UNIMAT_INIT_MODULE( &
        fort_denconst(1), &
        fort_tempconst(1), &
        outflow_pressure, &
        velfact, &
        probhix, &
        probhiy) 
      
      else
       ! do nothing 
      endif
      
       ! this loop occurs after user defined initialization.
      do im=1,num_materials
       call init_massfrac_parm(fort_denconst(im),massfrac_parm,im)
       call INTERNAL_material(fort_denconst(im),massfrac_parm, &
           fort_tempconst(im), &
           fort_energyconst(im),fort_material_type(im),im)
       call INTERNAL_material(fort_denconst(im),massfrac_parm, &
           fort_tempcutoff(im), &
           fort_energycutoff(im),fort_material_type(im),im)
      enddo
      
      if (ioproc.eq.1) then
       do im=1,num_materials
        print *,"im,energy ",im,fort_energyconst(im)
        print *,"im,energycutoff ",im,fort_energycutoff(im)
       enddo
      else if (ioproc.ne.0) then
       print *,"ioproc invalid"
       stop
      endif
      
      if (ioproc.eq.1) then
       print *,"prescribe_temperature_outflow (fortran)= ", &
        prescribe_temperature_outflow
       print *,"fort_solidheat_flag (fortran)= ", &
         fort_solidheat_flag
       print *,"density_at_depth ",density_at_depth
      
       print *,"fort: problox,y,z,hix,y,z ",problox,probloy,probloz, &
        probhix,probhiy,probhiz
       print *,"fort: problenx,y,z ",problenx,probleny,problenz
      
       if ((FORT_MUSHY_THICK.ge.one).and. &
           (FORT_MUSHY_THICK.le.four)) then
        ! do nothing
       else
        print *,"MUSHY_THICK invalid"
        stop
       endif
      
       if ((gravity_dir.lt.1).or.(gravity_dir.gt.SDIM).or. &
           (invert_gravity.lt.0).or.(invert_gravity.gt.1).or. &
           (ls_homflag.lt.0).or.(ls_homflag.gt.1).or. &
           (pres_homflag.lt.0).or.(pres_homflag.gt.1).or. &
           (vel_homflag.lt.0).or.(vel_homflag.gt.1).or. &
           (temp_homflag.lt.0).or.(temp_homflag.gt.1).or. &
           (species_homflag.lt.0).or.(species_homflag.gt.1).or. &
           (fort_solidheat_flag.lt.0).or. &
           (fort_solidheat_flag.gt.2).or. &
           (prescribe_temperature_outflow.lt.0).or. &
           (prescribe_temperature_outflow.gt.3)) then
        print *,"parameters invalid"
        stop
       endif
      
       print *,"fort:MUSHY_THICK ",FORT_MUSHY_THICK
      
       print *,"fort:gravity,gravity_dir,invert_gravity ",gravity, &
        gravity_dir,invert_gravity
       print *,"fort:pres,vel,temp,spec,ls homflag ",pres_homflag, &
        vel_homflag,temp_homflag,species_homflag,ls_homflag
      
       print *,"fort: stop_time ",fort_stop_time
      
       print *,"fort: rz ",levelrz
       do im=1,num_materials
        print *,"fort: im,FSI_flag ",im,FSI_flag(im)
       enddo
       print *,"fort: invert_solid_levelset ",invert_solid_levelset
       print *,"fort: denfact,velfact,xblob,yblob,zblob ", &
        denfact,velfact,xblob,yblob,zblob
       print *,"fort: radblob,probtype,adv_dir,adv_vel,axis_dir ", &
        radblob,probtype,adv_dir,adv_vel,axis_dir
       print *,"fort: xblob2,yblob2,zblob2,radblob2", &
        xblob2,yblob2,zblob2,radblob2
       print *,"fort: xblob3,yblob3,zblob3,radblob3", &
        xblob3,yblob3,zblob3,radblob3
       print *,"fort: xblob4,yblob4,zblob4,radblob4", &
        xblob4,yblob4,zblob4,radblob4
       print *,"fort: xblob5,yblob5,zblob5,radblob5", &
        xblob5,yblob5,zblob5,radblob5
       print *,"fort: xblob6,yblob6,zblob6,radblob6", &
        xblob6,yblob6,zblob6,radblob6
       print *,"fort: xblob7,yblob7,zblob7,radblob7", &
        xblob7,yblob7,zblob7,radblob7
       print *,"fort: xblob8,yblob8,zblob8,radblob8", &
        xblob8,yblob8,zblob8,radblob8
       print *,"fort: xblob9,yblob9,zblob9,radblob9", &
        xblob9,yblob9,zblob9,radblob9
       print *,"fort: xblob10,yblob10,zblob10,radblob10", &
        xblob10,yblob10,zblob10,radblob10
      
       print *,"fort: xactive,yactive,zactive,ractivexyz", &
        xactive,yactive,zactive,ractivex,ractivey,ractivez
       print *,"fort:vinletgas",vinletgas
       print *,"fort:rgasinlet",rgasinlet
       print *,"fort:advbot",advbot
       print *,"fort: inflow_pressure ",inflow_pressure
       print *,"fort: outflow_pressure ",outflow_pressure
       print *,"fort:period_time",period_time
       print *,"fort:twall",twall
       print *,"fort:fort_ZEYU_DCA_SELECT",fort_ZEYU_DCA_SELECT
      
       print *,"fort:end of override routine"
      
      else if (ioproc.ne.0) then
       print *,"ioproc invalid"
       stop
      endif
      
      do im=1,num_materials
       if (fort_material_type(im).eq.0) then
        ! do nothing
       else if (fort_material_type(im).eq.999) then
        ! do nothing
       else if ((fort_material_type(im).gt.0).and. &
                (fort_material_type(im).le.MAX_NUM_EOS)) then
        call debug_EOS(im)
       else
        print *,"fort_material_type invalid"
        stop
       endif
      enddo ! im=1..num_materials
      
      return
      end subroutine FORT_OVERRIDE
