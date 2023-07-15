#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif
        
#include "AMReX_FORT_INTEGER.H"
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

#include "EXTRAP_COMP.H"

      subroutine fort_blb_init( &
       blb_matrix_in, &
       blb_rhs_in, &
       blb_vel_in, &
       blb_int_mom_in, &
       blb_energy_in, &
       blb_mass_vel_in, &
       blb_vol_in, &
       blb_cen_int_in, &
       blb_cen_act_in, &
       blb_perim_in, &
       blb_perim_mat_in, &
       blb_triple_perim_in, &
       blb_cell_cnt_in, &
       blb_cellvol_cnt_in, &
       blb_mass_in, &
       blb_pres_in, &
       num_elements_blobclass_in) &
      bind(c,name='fort_blb_init')

      use probcommon_module
      IMPLICIT NONE

      INTEGER_T, INTENT(in) :: blb_matrix_in
      INTEGER_T, INTENT(in) :: blb_rhs_in
      INTEGER_T, INTENT(in) :: blb_vel_in
      INTEGER_T, INTENT(in) :: blb_int_mom_in
      INTEGER_T, INTENT(in) :: blb_energy_in
      INTEGER_T, INTENT(in) :: blb_mass_vel_in
      INTEGER_T, INTENT(in) :: blb_vol_in
      INTEGER_T, INTENT(in) :: blb_cen_int_in
      INTEGER_T, INTENT(in) :: blb_cen_act_in
      INTEGER_T, INTENT(in) :: blb_perim_in
      INTEGER_T, INTENT(in) :: blb_perim_mat_in
      INTEGER_T, INTENT(in) :: blb_triple_perim_in
      INTEGER_T, INTENT(in) :: blb_cell_cnt_in
      INTEGER_T, INTENT(in) :: blb_cellvol_cnt_in
      INTEGER_T, INTENT(in) :: blb_mass_in
      INTEGER_T, INTENT(in) :: blb_pres_in
      INTEGER_T, INTENT(in) :: num_elements_blobclass_in
    
      BLB_MATRIX=blb_matrix_in
      BLB_RHS=blb_rhs_in
      BLB_VEL=blb_vel_in
      BLB_INT_MOM=blb_int_mom_in
      BLB_ENERGY=blb_energy_in
      BLB_MASS_VEL=blb_mass_vel_in
      BLB_VOL=blb_vol_in
      BLB_CEN_INT=blb_cen_int_in
      BLB_CEN_ACT=blb_cen_act_in
      BLB_PERIM=blb_perim_in
      BLB_PERIM_MAT=blb_perim_mat_in
      BLB_TRIPLE_PERIM=blb_triple_perim_in
      BLB_CELL_CNT=blb_cell_cnt_in
      BLB_CELLVOL_CNT=blb_cellvol_cnt_in
      BLB_MASS=blb_mass_in
      BLB_PRES=blb_pres_in
      num_elements_blobclass=num_elements_blobclass_in

      if ((BLB_MATRIX.eq.0).and. &
          (BLB_RHS.eq.BLB_MATRIX+ &
             3*(2*AMREX_SPACEDIM)*(2*AMREX_SPACEDIM)).and. &
          (BLB_VEL.eq.BLB_RHS+3*(2*AMREX_SPACEDIM)).and. &
          (BLB_INT_MOM.eq.BLB_VEL+3*(2*AMREX_SPACEDIM)).and. &
          (BLB_ENERGY.eq.BLB_INT_MOM+2*(2*AMREX_SPACEDIM)).and. &
          (BLB_MASS_VEL.eq.BLB_ENERGY+1).and. &
          (BLB_VOL.eq.BLB_MASS_VEL+3).and. &
          (BLB_CEN_INT.eq.BLB_VOL+1).and. &
          (BLB_CEN_ACT.eq.BLB_CEN_INT+AMREX_SPACEDIM).and. &
          (BLB_PERIM.eq.BLB_CEN_ACT+AMREX_SPACEDIM).and. &
          (BLB_PERIM_MAT.eq.BLB_PERIM+1).and. &
          (BLB_TRIPLE_PERIM.eq.BLB_PERIM_MAT+num_materials).and. &
          (BLB_CELL_CNT.eq.BLB_TRIPLE_PERIM+num_materials*num_materials).and. &
          (BLB_CELLVOL_CNT.eq.BLB_CELL_CNT+1).and. &
          (BLB_MASS.eq.BLB_CELLVOL_CNT+1).and. &
          (BLB_PRES.eq.BLB_MASS+1).and. &
          (num_elements_blobclass.eq.BLB_PRES+1)) then
          ! do nothing
      else
       print *,"BLB parameters invalid"
       print *,"BLB_MATRIX ",BLB_MATRIX
       print *,"BLB_RHS ",BLB_RHS
       print *,"BLB_VEL ",BLB_VEL
       print *,"BLB_INT_MOM ",BLB_INT_MOM
       print *,"BLB_ENERGY ",BLB_ENERGY
       print *,"BLB_MASS_VEL ",BLB_MASS_VEL
       print *,"BLB_VOL ",BLB_VOL
       print *,"BLB_CEN_INT ",BLB_CEN_INT
       print *,"BLB_CEN_ACT ",BLB_CEN_ACT
       print *,"BLB_PERIM ",BLB_PERIM
       print *,"BLB_PERIM_MAT ",BLB_PERIM_MAT
       print *,"BLB_TRIPLE_PERIM ",BLB_TRIPLE_PERIM
       print *,"BLB_CELL_CNT ",BLB_CELL_CNT
       print *,"BLB_CELLVOL_CNT ",BLB_CELLVOL_CNT
       print *,"BLB_MASS ",BLB_MASS
       print *,"BLB_PRES ",BLB_PRES
       print *,"num_elements_blobclass ",num_elements_blobclass
       print *,"num_materials=",num_materials
       print *,"AMREX_SPACEDIM=",AMREX_SPACEDIM
       stop
      endif

      end subroutine fort_blb_init

       !called from NavierStokes.cpp: void fortran_deallocate_parameters
      subroutine fort_deallocate_module( &
        ) &
      bind(c,name='fort_deallocate_module')

      use probcommon_module
      use global_utility_module
      IMPLICIT NONE

      if (is_in_probtype_list().eq.1) then
       call SUB_DEALLOCATE_MODULE()
      else if (is_in_probtype_list().eq.0) then
       ! do nothing
      else
       print *,"is_in_probtype_list invalid"
       stop
      endif

      return
      end subroutine fort_deallocate_module

       ! ns.mof_ordering overrides this
       ! called from: fortran_parameters()  (before "pp.queryAdd")
       ! called from: NavierStokes::read_params()  (before "pp.queryAdd")
      subroutine fort_mof_ordering_override( &
        mof_ordering_local, &
        mof_error_ordering_local, &
        FSI_flag_temp) &
      bind(c,name='fort_mof_ordering_override')

      use probcommon_module
      IMPLICIT NONE

      INTEGER_T, INTENT(in) :: mof_error_ordering_local
      INTEGER_T, INTENT(in) :: FSI_flag_temp(num_materials)
      INTEGER_T, INTENT(out) :: mof_ordering_local(num_materials)
      INTEGER_T :: im
      INTEGER_T :: local_FSI_flag

      if ((num_materials.lt.1).or.(num_materials.gt.9999)) then
       print *,"num_materials invalid"
       stop
      endif

      do im=1,num_materials
       mof_ordering_local(im)=0

       local_FSI_flag=FSI_flag_temp(im)

       if ((local_FSI_flag.eq.FSI_FLUID).or. &
           (local_FSI_flag.eq.FSI_FLUID_NODES_INIT)) then
        ! do nothing, tessellating
       else if (local_FSI_flag.eq.FSI_PRESCRIBED_PROBF90) then
        mof_ordering_local(im)=1  ! non-tessellating
       else if (local_FSI_flag.eq.FSI_PRESCRIBED_NODES) then
        mof_ordering_local(im)=1  ! non-tessellating
       else if ((local_FSI_flag.eq.FSI_ICE_PROBF90).or. &
                (local_FSI_flag.eq.FSI_ICE_STATIC).or. &
                (local_FSI_flag.eq.FSI_ICE_NODES_INIT)) then
        ! do nothing, tessellating
   
        ! FSI elastic link w/Kourosh (sci_clsvof.F90)
       else if (local_FSI_flag.eq.FSI_SHOELE_CTML) then
        mof_ordering_local(im)=1 ! non-tessellating
       else if (local_FSI_flag.eq.FSI_RIGID_NOTPRESCRIBED) then
        mof_ordering_local(im)=1 ! tessellating
       else
        print *,"local_FSI_flag invalid in fort_mof_ordering_override"
        print *,"im,local_FSI_flag ",im,local_FSI_flag
        stop
       endif
      enddo !im=1,num_materials

      ! default: centroid farthest from uncaptured centroid.
      if (mof_error_ordering_local.eq.0) then

       do im=1,num_materials

        local_FSI_flag=FSI_flag_temp(im)

        if ((local_FSI_flag.eq.FSI_FLUID).or. &
            (local_FSI_flag.eq.FSI_FLUID_NODES_INIT)) then
         mof_ordering_local(im)=num_materials
        else if (local_FSI_flag.eq.FSI_PRESCRIBED_PROBF90) then 
         mof_ordering_local(im)=1 ! non-tessellating
        else if (local_FSI_flag.eq.FSI_PRESCRIBED_NODES) then 
         mof_ordering_local(im)=1 ! non-tessellating
        else if ((local_FSI_flag.eq.FSI_ICE_PROBF90).or. &
                 (local_FSI_flag.eq.FSI_ICE_STATIC).or. &
                 (local_FSI_flag.eq.FSI_ICE_NODES_INIT)) then 
         mof_ordering_local(im)=num_materials ! tessellating

        else if (local_FSI_flag.eq.FSI_SHOELE_CTML) then  
         mof_ordering_local(im)=1  ! non-tessellating
        else if (local_FSI_flag.eq.FSI_RIGID_NOTPRESCRIBED) then 
         mof_ordering_local(im)=1  ! tessellating
        else
         print *,"local_FSI_flag invalid in fort_mof_ordering_override"
         print *,"im,local_FSI_flag ",im,local_FSI_flag
         stop
        endif

       enddo !im=1,num_materials

       ! impinge jets unlike material
       if ((probtype.eq.530).and.(AMREX_SPACEDIM.eq.3)) then
        if (axis_dir.eq.1) then
         mof_ordering_local(2)=num_materials+1! make gas have low priority
        else if (axis_dir.ne.0) then
         print *,"axis_dir invalid probtype=530"
         stop
        endif
       endif 

       ! 2d colliding droplets, boiling, freezing problems
       if ((probtype.eq.55).and.(AMREX_SPACEDIM.eq.2)) then
        if (radblob7.gt.zero) then
         mof_ordering_local(2)=num_materials+1! make gas have low priority
        else if (radblob7.le.zero) then
         ! do nothing
        else
         print *,"radblob7 is NaN"
         stop
        endif
        if (axis_dir.eq.0) then
         ! do nothing
        else if (axis_dir.eq.1) then
         ! 0=water 1=gas 2=ice 3=cold plate
         mof_ordering_local(3)=1
        else if (axis_dir.eq.5) then
         ! 0=water 1=gas 2=ice 3=cold plate
         if (num_materials.eq.4) then
          mof_ordering_local(1)=3
          mof_ordering_local(2)=1
          mof_ordering_local(3)=2
          mof_ordering_local(4)=1
         else
          print *,"expecting num_materials==4 axis_dir==5 probtype==55"
          stop
         endif
         ! 0=water 1=vapor 2=hot plate or
         ! 0=water 1=vapor 2=gas 3=hot plate 
        else if (axis_dir.eq.6) then  ! nucleate boiling incompressible
         mof_ordering_local(num_materials)=1
         ! 0=water 1=vapor 2=hot plate or
         ! 0=water 1=vapor 2=gas 3=hot plate 
        else if (axis_dir.eq.7) then  ! nucleate boiling compressible
         mof_ordering_local(num_materials)=1
        else
         print *,"axis_dir invalid probtype.eq.55"
         stop
        endif
       endif

       if (probtype.eq.540) then
        if ((radblob4.gt.zero).and.(radblob3.gt.zero)) then
         print *,"conflict of parametrs for 540"
         stop
        else if (radblob3.gt.zero) then  
         mof_ordering_local(2)=num_materials+1!make gas have low priority
        else if (radblob4.gt.zero) then
         mof_ordering_local(3)=num_materials+1!make filament gas low priority
        endif
       endif

       if (probtype.eq.202) then  ! liquidlens
        mof_ordering_local(2)=1 ! make (circle) material 2 have high priority
       endif 

       if ((probtype.eq.17).and. &
           (num_materials.eq.3).and. &
           (1.eq.0)) then! droplet impact 3mat
        mof_ordering_local(2)=1 ! make gas material 2 have high priority
       endif

      else if (mof_error_ordering_local.eq.1) then

       ! mof_ordering_local already init above.
  
      else
       print *,"mof_error_ordering_local invalid"
       stop
      endif

      end subroutine fort_mof_ordering_override

       ! fortran_parameters is called from main.cpp prior to:
       !  1. AmrCore* amrptr = new AmrCore();
       !  2. amrptr->init(strt_time,stop_time);
       ! fort_override is called from 
       !  NavierStokes.cpp: void fortran_parameters
      subroutine fort_override( &
        cc_int_size, &
        ccmax_level, &
        ccn_cell, &
        ccbfact_space_order, &
        ccbfact_time_order, &
        ccprescribe_temperature_outflow, &
        ccsolidheat_flag, &
        rz_flag, &
        ccFSI_flag, &
        ccnum_local_aux_grids, &
        ccZEYU_DCA_SELECT, &
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
        ccmaterial_type, &
        ccnten, &
        ccDrhoDT, &
        cctempconst, &
        ccinitial_temperature, &
        cctempcutoff, &
        cctempcutoffmax, &
        ccstiffPINF, &
        ccR_Palmore_Desjardins, &
        ccstiffCP, &
        ccstiffCV, &
        ccstiffGAMMA, &
        ccdenconst, &
        ccden_floor, &
        ccden_ceiling, &
        cccavdenconst, &
        ccviscconst, &
        ccviscconst_eddy_wall, &
        ccviscconst_eddy_bulk, &
        ccheatviscconst_eddy_wall, &
        ccheatviscconst_eddy_bulk, &
        ccthermal_microlayer_size, &
        ccshear_microlayer_size, &
        ccbuoyancy_microlayer_size, &
        ccphasechange_microlayer_size, &
        ccviscosity_state_model, &
        ccelastic_viscosity, &
        ccelastic_time, &
        ccviscoelastic_model, &
        ccstore_elastic_data, &
        ccheatflux_factor, &
        ccheatviscconst, &
        ccprerecalesce_heatviscconst, &
        ccprerecalesce_viscconst, &
        ccprerecalesce_stiffCP, &
        ccprerecalesce_stiffCV, &
        ccspeciesconst, &
        ccspeciesviscconst, &
        cclatent_heat, &
        cclatent_heat_slope, &
        cclatent_heat_T0, &
        cclatent_heat_min, &
        ccsaturation_temp, &
        ccreference_pressure, &
        ccmolar_mass, &
        ccspecies_molar_mass, &
        cctension, &
        cctension_init, &
        cctension_slope, &
        cctension_T0, &
        cctension_min, &
        ccprefreeze_tension, &
        ccgravity_vector, &
        ccstop_time, &
        ccCarreau_alpha, &
        ccCarreau_beta, &
        ccCarreau_n, &
        ccCarreau_mu_inf, &
        ccshear_thinning_fluid, &
        ccpolymer_factor, &
        ccconcentration, &
        ccetaL, &
        ccetaS, &
        ccetaP, &
        ccvisc_coef, &
        ccangular_velocity, &
        ccgrid_stretching_parameter, &
        ioproc) &
      bind(c,name='fort_override')

      use probcommon_module
      use LegendreNodes
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
      use ROTATING_ANNULUS_module
      use CAVITY_PHASE_CHANGE_module
      use ICE_ON_SUBSTRATE_module
      use SIMPLE_PALMORE_DESJARDINS_module
      use SIMPLE_KASSEMI_module
      use DROP_IN_SHEAR_module
      use MITSUHIRO_MELTING_module
      use AHMED_ICE_RESISTANT_module
      use FABRIC_DROP_MODULE
      use CRYOGENIC_TANK1_module
      use CRYOGENIC_TANK2_module
      use CRYOGENIC_TANK_MK_module
      use flexible_plate_impact_module
      use CONE3D_module
      use YAOHONG_INKJET_module
      use WAVY_Channel_module
      use rigid_FSI_module
      use passive_advect_module ! probtype=28,29,31
      
      IMPLICIT NONE
      
      INTEGER_T, INTENT(in) :: cc_int_size
      INTEGER_T, INTENT(in) :: ccmax_level
      INTEGER_T, INTENT(in) :: ccn_cell(SDIM)
      INTEGER_T, INTENT(in) :: ccbfact_space_order(0:ccmax_level)
      INTEGER_T, INTENT(in) :: ccbfact_time_order
      INTEGER_T, INTENT(in) :: ccnum_materials
      
      INTEGER_T, INTENT(in) :: ccnten
      REAL_T, INTENT(in) :: ccgravity_vector(SDIM)
      INTEGER_T, INTENT(in) :: ccFSI_flag(ccnum_materials)
      INTEGER_T, INTENT(in) :: ccnum_local_aux_grids
      INTEGER_T, INTENT(in) :: ccZEYU_DCA_SELECT
      INTEGER_T, INTENT(in) :: ccprescribe_temperature_outflow
      INTEGER_T, INTENT(in) :: ccsolidheat_flag
      INTEGER_T, INTENT(in) :: rz_flag
      INTEGER_T, INTENT(in) :: ioproc
      INTEGER_T, INTENT(in) :: ccprobtype,ccadv_dir,ccaxis_dir
      
      REAL_T, INTENT(in) :: ccdenfact,ccvelfact
      REAL_T, INTENT(in) :: ccxblob,ccyblob,cczblob,ccradblob
      REAL_T, INTENT(in) :: ccxblob2,ccyblob2,cczblob2,ccradblob2
      REAL_T, INTENT(in) :: ccxblob3,ccyblob3,cczblob3,ccradblob3
      REAL_T, INTENT(in) :: ccxblob4,ccyblob4,cczblob4,ccradblob4
      REAL_T, INTENT(in) :: ccxblob5,ccyblob5,cczblob5,ccradblob5
      REAL_T, INTENT(in) :: ccxblob6,ccyblob6,cczblob6,ccradblob6
      REAL_T, INTENT(in) :: ccxblob7,ccyblob7,cczblob7,ccradblob7
      REAL_T, INTENT(in) :: ccxblob8,ccyblob8,cczblob8,ccradblob8
      REAL_T, INTENT(in) :: ccxblob9,ccyblob9,cczblob9,ccradblob9
      REAL_T, INTENT(in) :: ccxblob10,ccyblob10,cczblob10,ccradblob10
      REAL_T, INTENT(in) :: ccxactive,ccyactive,cczactive
      REAL_T, INTENT(in) :: ccractivex,ccractivey,ccractivez
      REAL_T, INTENT(in) :: ccadv_vel,ccrgasinlet
      REAL_T, INTENT(in) :: ccvinletgas,cctwall
      REAL_T, INTENT(in) :: ccadvbot
      REAL_T, INTENT(in) :: ccinflow_pressure
      REAL_T, INTENT(in) :: ccoutflow_pressure
      REAL_T, INTENT(in) :: ccperiod_time
      REAL_T, INTENT(in) :: ccproblox,ccprobloy,ccprobloz
      REAL_T, INTENT(in) :: ccprobhix,ccprobhiy,ccprobhiz
      REAL_T, INTENT(in) :: ccstop_time
      
      INTEGER_T, INTENT(in) :: ccnum_species_var
      
      INTEGER_T, INTENT(in) :: ccnum_materials_viscoelastic
      INTEGER_T :: nelastic
      
      INTEGER_T, INTENT(in) :: ccnum_state_material
      INTEGER_T, INTENT(in) :: ccnum_state_base
      INTEGER_T, INTENT(in) :: ccngeom_raw
      INTEGER_T, INTENT(in) :: ccngeom_recon
      
      INTEGER_T, INTENT(in) :: ccmaterial_type(ccnum_materials)
      REAL_T, INTENT(in) :: ccDrhoDT(ccnum_materials)
      REAL_T, INTENT(in) :: cctempconst(ccnum_materials)
      REAL_T, INTENT(in) :: ccinitial_temperature(ccnum_materials)
      REAL_T, INTENT(in) :: cctempcutoff(ccnum_materials)
      REAL_T, INTENT(in) :: cctempcutoffmax(ccnum_materials)
      REAL_T, INTENT(in) :: ccstiffPINF(ccnum_materials)
      REAL_T, INTENT(in) :: ccR_Palmore_Desjardins
      REAL_T, INTENT(in) :: ccstiffCP(ccnum_materials)
      REAL_T, INTENT(in) :: ccstiffCV(ccnum_materials)
      REAL_T, INTENT(in) :: ccstiffGAMMA(ccnum_materials)
      REAL_T, INTENT(in) :: ccdenconst(ccnum_materials)
      REAL_T, INTENT(in) :: ccden_floor(ccnum_materials)
      REAL_T, INTENT(in) :: ccden_ceiling(ccnum_materials)
      REAL_T, INTENT(in) :: cccavdenconst(ccnum_materials)
      REAL_T, INTENT(in) :: ccviscconst(ccnum_materials)

      REAL_T, INTENT(in) :: ccviscconst_eddy_wall(ccnum_materials)
      REAL_T, INTENT(in) :: ccviscconst_eddy_bulk(ccnum_materials)
      REAL_T, INTENT(in) :: ccheatviscconst_eddy_wall(ccnum_materials)
      REAL_T, INTENT(in) :: ccheatviscconst_eddy_bulk(ccnum_materials)

      REAL_T, INTENT(in) :: ccthermal_microlayer_size(ccnum_materials)
      REAL_T, INTENT(in) :: ccshear_microlayer_size(ccnum_materials)
      REAL_T, INTENT(in) :: ccbuoyancy_microlayer_size(ccnum_materials)
      REAL_T, INTENT(in) :: ccphasechange_microlayer_size(ccnum_materials)

      INTEGER_T, INTENT(in) :: ccviscosity_state_model(ccnum_materials)
      REAL_T, INTENT(in) :: ccelastic_viscosity(ccnum_materials)
      REAL_T, INTENT(in) :: ccelastic_time(ccnum_materials)
      INTEGER_T, INTENT(in) :: ccviscoelastic_model(ccnum_materials)
      INTEGER_T, INTENT(in) :: ccstore_elastic_data(ccnum_materials)
      REAL_T, INTENT(in) :: ccheatflux_factor(ccnum_materials)
      REAL_T, INTENT(in) :: ccheatviscconst(ccnum_materials)
      REAL_T, INTENT(in) :: ccprerecalesce_heatviscconst(ccnum_materials)
      REAL_T, INTENT(in) :: ccprerecalesce_viscconst(ccnum_materials)
      REAL_T, INTENT(in) :: ccprerecalesce_stiffCP(ccnum_materials)
      REAL_T, INTENT(in) :: ccprerecalesce_stiffCV(ccnum_materials)
      REAL_T, INTENT(in) :: &
        ccspeciesconst((ccnum_species_var+1)*ccnum_materials)
      REAL_T, INTENT(in) :: &
        ccspeciesviscconst((ccnum_species_var+1)*ccnum_materials)
      REAL_T, INTENT(in) :: cclatent_heat(2*ccnten)
      REAL_T, INTENT(in) :: cclatent_heat_slope(2*ccnten)
      REAL_T, INTENT(in) :: cclatent_heat_T0(2*ccnten)
      REAL_T, INTENT(in) :: cclatent_heat_min(2*ccnten)
      REAL_T, INTENT(in) :: ccsaturation_temp(2*ccnten)
      REAL_T, INTENT(in) :: ccreference_pressure(2*ccnten)
      REAL_T, INTENT(in) :: ccmolar_mass(ccnum_materials)
      REAL_T, INTENT(in) :: ccspecies_molar_mass(ccnum_species_var+1)
      REAL_T, INTENT(in) :: cctension(ccnten)
      REAL_T, INTENT(in) :: cctension_init(ccnten)
      REAL_T, INTENT(in) :: cctension_slope(ccnten)
      REAL_T, INTENT(in) :: cctension_T0(ccnten)
      REAL_T, INTENT(in) :: cctension_min(ccnten)
      REAL_T, INTENT(in) :: ccprefreeze_tension(ccnten)
      
      INTEGER_T, INTENT(in) :: ccn_sites
      REAL_T, INTENT(in) :: ccnucleation_init_time
      REAL_T, INTENT(in) :: ccpos_sites(5000)
     
      REAL_T, INTENT(in) :: ccCarreau_alpha(ccnum_materials)
      REAL_T, INTENT(in) :: ccCarreau_beta(ccnum_materials)
      REAL_T, INTENT(in) :: ccCarreau_n(ccnum_materials)
      REAL_T, INTENT(in) :: ccCarreau_mu_inf(ccnum_materials)

      INTEGER_T, INTENT(in) :: ccshear_thinning_fluid(ccnum_materials)

      REAL_T, INTENT(in) :: ccpolymer_factor(ccnum_materials)
      REAL_T, INTENT(in) :: ccconcentration(ccnum_materials)
      REAL_T, INTENT(in) :: ccetaL(ccnum_materials)
      REAL_T, INTENT(in) :: ccetaS(ccnum_materials)
      REAL_T, INTENT(in) :: ccetaP(ccnum_materials)
      REAL_T, INTENT(in) :: ccgrid_stretching_parameter(SDIM)

      REAL_T, INTENT(in) :: ccvisc_coef
      REAL_T, INTENT(in) :: ccangular_velocity


      character*12 namestr1
      character*13 namestr2
      INTEGER_T i
      INTEGER_T local_dir
      
      INTEGER_T im,iten
      INTEGER_T level,bfactmax
      REAL_T :: massfrac_parm(ccnum_species_var+1)
      INTEGER_T :: fort_double_size,fort_int_size
      
      probtype=ccprobtype
      num_materials=ccnum_materials
      if (num_materials.lt.MAX_NUM_MATERIALS) then
       ! do nothing
      else
       print *,"increase MAX_NUM_MATERIALS; aborting"
       stop
      endif
      num_interfaces=( (num_materials-1)*(num_materials-1)+ &
             (num_materials-1) )/2
      if (num_interfaces.eq.ccnten) then
       ! do nothing
      else
       print *,"num_interfaces<>ccnten"
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
      probtype_list_size=21
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
      used_probtypes(11)=222 ! CONE3D_module 
      used_probtypes(12)=2011 ! YAOHONG_INKJET 
      used_probtypes(13)=425  ! AHMED_ICE_RESISTANT
      used_probtypes(14)=7001 ! FABRIC_DROP
      used_probtypes(15)=915 ! wavy channel or Tomas
      used_probtypes(16)=411 ! CAV3D_module
      used_probtypes(17)=28  ! passive_advect_module
      used_probtypes(18)=29  ! passive_advect_module
      used_probtypes(19)=31  ! passive_advect_module
      used_probtypes(20)=710 ! CAVITY_PHASE_CHANGE
      used_probtypes(21)=82  ! Differentially Heated Rotating Annulus: 
                             ! ROTATING_ANNULUS
      
      SUB_INIT_MODULE=>INIT_STUB_MODULE
      SUB_DEALLOCATE_MODULE=>DEALLOCATE_STUB_MODULE
      SUB_LS=>STUB_LS
      SUB_OVERRIDE_TAGFLAG=>STUB_OVERRIDE_TAGFLAG
      SUB_AUX_DATA=>STUB_AUX_DATA
      SUB_OVERRIDE_FSI_SIGN_LS_VEL_TEMP=>STUB_OVERRIDE_FSI_SIGN_LS_VEL_TEMP
      SUB_GET_OUTSIDE_POINT=>STUB_GET_OUTSIDE_POINT
      SUB_BOUNDING_BOX_AUX=>STUB_BOUNDING_BOX_AUX
      SUB_VEL=>STUB_VEL
      SUB_EOS=>EOS_STUB

      SUB_VARIABLE_SURFACE_TENSION=>STUB_VARIABLE_SURFACE_TENSION
      SUB_VARIABLE_LATENT_HEAT=>STUB_VARIABLE_LATENT_HEAT

      SUB_UNITLESS_EXPANSION_FACTOR=>STUB_UNITLESS_EXPANSION_FACTOR
      SUB_INTERNAL_GRAVITY_WAVE_FLAG=>STUB_INTERNAL_GRAVITY_WAVE_FLAG

      SUB_dVdT=>dVdT_STUB
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
      SUB_ICE_SUBSTRATE_DISTANCE=>STUB_ICE_SUBSTRATE_DISTANCE
      SUB_CFL_HELPER=>STUB_CFL_HELPER
      SUB_correct_pres_rho_hydrostatic=>STUB_correct_pres_rho_hydrostatic
      SUB_ASSIMILATE=>STUB_ASSIMILATE
      SUB_SUMINT=>STUB_SUMINT
      SUB_check_vel_rigid=>STUB_check_vel_rigid
      SUB_clamped_LS_no_scale=>STUB_clamped_LS

      SUB_wallfunc=>STUB_wallfunc

      SUB_MAPPING_WEIGHT_COEFF=>STUB_MAPPING_WEIGHT_COEFF

      SUB_INIT_REGIONS_LIST=>STUB_INIT_REGIONS_LIST
      SUB_CHARFN_REGION=>STUB_CHARFN_REGION
      SUB_DELETE_REGIONS_LIST=>STUB_DELETE_REGIONS_LIST
      
      SUB_THERMAL_K=>STUB_THERMAL_K
      SUB_INTERFACE_TEMPERATURE=>STUB_INTERFACE_TEMPERATURE
      SUB_MDOT=>STUB_MDOT
      SUB_K_EFFECTIVE=>STUB_K_EFFECTIVE

      SUB_reference_wavelen=>STUB_reference_wavelen

      SUB_OPEN_CASFILE=>STUB_OPEN_CASFILE
      SUB_OPEN_AUXFILE=>STUB_OPEN_AUXFILE
      SUB_ORDER_NODES=>STUB_ORDER_NODES
      SUB_FSI_SLICE=>STUB_FSI_SLICE

      SUB_T0_Boussinesq=>STUB_T0_Boussinesq
      SUB_V0_Coriolis=>STUB_V0_Coriolis

      if (probtype.eq.421) then
       SUB_INIT_MODULE=>INIT_CRYOGENIC_TANK1_MODULE
       SUB_LS=>CRYOGENIC_TANK1_LS
       SUB_VEL=>CRYOGENIC_TANK1_VEL
       SUB_EOS=>EOS_CRYOGENIC_TANK1
       SUB_dVdT=>dVdT_CRYOGENIC_TANK1
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
       SUB_dVdT=>dVdT_CRYOGENIC_TANK2
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
       SUB_dVdT=>dVdT_CRYOGENIC_TANK_MK
       SUB_SOUNDSQR=>SOUNDSQR_CRYOGENIC_TANK_MK
       SUB_INTERNAL=>INTERNAL_CRYOGENIC_TANK_MK
       SUB_TEMPERATURE=>TEMPERATURE_CRYOGENIC_TANK_MK
       SUB_PRES=>CRYOGENIC_TANK_MK_PRES
       SUB_STATE=>CRYOGENIC_TANK_MK_STATE
       SUB_LS_BC=>CRYOGENIC_TANK_MK_LS_BC
       SUB_VEL_BC=>CRYOGENIC_TANK_MK_VEL_BC
       SUB_PRES_BC=>CRYOGENIC_TANK_MK_PRES_BC
       SUB_STATE_BC=>CRYOGENIC_TANK_MK_STATE_BC
       SUB_HEATSOURCE=>CRYOGENIC_TANK_MK_HEATSOURCE

       SUB_SUMINT=>CRYOGENIC_TANK_MK_SUMINT
       SUB_INIT_REGIONS_LIST=>CRYOGENIC_TANK_MK_INIT_REGIONS_LIST
       SUB_CHARFN_REGION=>CRYOGENIC_TANK_MK_CHARFN_REGION
       SUB_THERMAL_K=>CRYOGENIC_TANK_MK_THERMAL_K

       SUB_wallfunc=>CRYOGENIC_TANK_MK_wallfunc
       SUB_K_EFFECTIVE=>CRYOGENIC_TANK_MK_K_EFFECTIVE

       SUB_OPEN_CASFILE=>CRYOGENIC_TANK_MK_OPEN_CASFILE
       SUB_OPEN_AUXFILE=>CRYOGENIC_TANK_MK_OPEN_AUXFILE

       SUB_BOUNDING_BOX_AUX=>CRYOGENIC_TANK_MK_BOUNDING_BOX_AUX
       SUB_OVERRIDE_FSI_SIGN_LS_VEL_TEMP=> &
         CRYOGENIC_TANK_MK_OVERRIDE_FSI_SIGN_LS_VEL_TEMP
       SUB_GET_OUTSIDE_POINT=>CRYOGENIC_TANK_MK_GET_OUTSIDE_POINT

       SUB_MAPPING_WEIGHT_COEFF=>CRYOGENIC_TANK_MK_MAPPING_WEIGHT_COEFF

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
       SUB_ASSIMILATE=>DROP_IN_SHEAR_ASSIMILATE

      else if ((probtype.eq.28).or. &
               (probtype.eq.29).or. &
               (probtype.eq.31)) then

       SUB_INIT_MODULE=>INIT_passive_advect_MODULE
       SUB_LS=>passive_advect_LS
       SUB_VEL=>passive_advect_VEL
       SUB_PRES=>passive_advect_PRES
       SUB_STATE=>passive_advect_STATE
       SUB_LS_BC=>passive_advect_LS_BC
       SUB_VEL_BC=>passive_advect_VEL_BC
       SUB_PRES_BC=>passive_advect_PRES_BC
       SUB_STATE_BC=>passive_advect_STATE_BC
       SUB_clamped_LS_no_scale=>passive_advect_clamped_LS
       SUB_CFL_HELPER=>passive_advect_CFL_HELPER

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

      else if (probtype.eq.915) then

       SUB_INIT_MODULE=>INIT_WAVY_MODULE
       SUB_LS=>WAVY_INIT_LS
       SUB_VEL=>WAVY_INIT_VEL
       SUB_PRES=>WAVY_INIT_PRES
       SUB_STATE=>WAVY_INIT_STATE
       SUB_LS_BC=>WAVY_LS_BC

       SUB_AUX_DATA=>WAVY_AUX_DATA
       SUB_OVERRIDE_TAGFLAG=>WAVY_OVERRIDE_TAGFLAG
       SUB_BOUNDING_BOX_AUX=>WAVY_BOUNDING_BOX_AUX

       SUB_VEL_BC=>WAVY_VEL_BC
       SUB_PRES_BC=>WAVY_PRES_BC
       SUB_STATE_BC=>WAVY_STATE_BC
       SUB_HEATSOURCE=>WAVY_HEATSOURCE

      else if (probtype.eq.7001) then

       SUB_INIT_MODULE=>INIT_FABRIC_DROP_MODULE
       SUB_LS=>FABRIC_DROP_LS
       SUB_VEL=>FABRIC_DROP_VEL
       SUB_PRES=>FABRIC_DROP_PRES
       SUB_STATE=>FABRIC_DROP_STATE
       SUB_LS_BC=>FABRIC_DROP_LS_BC
       SUB_VEL_BC=>FABRIC_DROP_VEL_BC
       SUB_PRES_BC=>FABRIC_DROP_PRES_BC
       SUB_STATE_BC=>FABRIC_DROP_STATE_BC
       SUB_HEATSOURCE=>FABRIC_DROP_HEATSOURCE

      else if (probtype.eq.411) then ! cav3D.F90
       SUB_INIT_MODULE=>INIT_CAV3D_MODULE
       SUB_LS=>CAV3D_LS
       SUB_VEL=>CAV3D_VEL
       SUB_PRES=>CAV3D_PRES
       SUB_STATE=>CAV3D_STATE
       SUB_LS_BC=>CAV3D_LS_BC
       SUB_VEL_BC=>CAV3D_VEL_BC
       SUB_PRES_BC=>CAV3D_PRES_BC
       SUB_STATE_BC=>CAV3D_STATE_BC
       SUB_HEATSOURCE=>CAV3D_HEATSOURCE

       SUB_OPEN_CASFILE=>OPEN_CAV3D_CASFILE
       SUB_ORDER_NODES=>CAV3D_ORDER_NODES
       SUB_FSI_SLICE=>CAV3D_SLICE

      else if (probtype.eq.425) then
       SUB_INIT_MODULE=>INIT_AHMED_ICE_RESISTANT_MODULE
       SUB_LS=>AHMED_ICE_RESISTANT_LS
       SUB_VEL=>AHMED_ICE_RESISTANT_VEL
       SUB_PRES=>AHMED_ICE_RESISTANT_PRES
       SUB_STATE=>AHMED_ICE_RESISTANT_STATE
       SUB_LS_BC=>AHMED_ICE_RESISTANT_LS_BC
       SUB_VEL_BC=>AHMED_ICE_RESISTANT_VEL_BC
       SUB_PRES_BC=>AHMED_ICE_RESISTANT_PRES_BC
       SUB_STATE_BC=>AHMED_ICE_RESISTANT_STATE_BC
       SUB_HEATSOURCE=>AHMED_ICE_RESISTANT_HEATSOURCE
       SUB_ASSIMILATE=>AHMED_ICE_RESISTANT_ASSIMILATE

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
       SUB_dVdT=>dVdT_KASSEMI_MK
       SUB_SOUNDSQR=>SOUNDSQR_KASSEMI_MK
       SUB_INTERNAL=>INTERNAL_KASSEMI_MK
       SUB_TEMPERATURE=>TEMPERATURE_KASSEMI_MK

       SUB_PRES_BC=>SIMPLE_KASSEMI_PRES_BC
       SUB_STATE_BC=>SIMPLE_KASSEMI_STATE_BC
       SUB_HEATSOURCE=>SIMPLE_KASSEMI_HEATSOURCE
       SUB_SUMINT=>SIMPLE_KASSEMI_SUMINT ! compare with analytical
       SUB_ASSIMILATE=>SIMPLE_KASSEMI_ASSIMILATE
      else if (probtype.eq.222) then
       SUB_INIT_MODULE=>INIT_CONE3D_MODULE
       SUB_LS=>CONE3D_LS
       SUB_VEL=>CONE3D_VEL
       SUB_PRES=>CONE3D_PRES
       SUB_STATE=>CONE3D_STATE
       SUB_LS_BC=>CONE3D_LS_BC
       SUB_VEL_BC=>CONE3D_VEL_BC
       SUB_PRES_BC=>CONE3D_PRES_BC
       SUB_STATE_BC=>CONE3D_STATE_BC
      else if (probtype.eq.2011) then
       SUB_INIT_MODULE=>INIT_YAOHONG_INKJET_MODULE
       SUB_LS=>YAOHONG_INKJET_LS
       SUB_VEL=>YAOHONG_INKJET_VEL
       SUB_PRES=>YAOHONG_INKJET_PRES
       SUB_STATE=>YAOHONG_INKJET_STATE
       SUB_LS_BC=>YAOHONG_INKJET_LS_BC
       SUB_VEL_BC=>YAOHONG_INKJET_VEL_BC
       SUB_PRES_BC=>YAOHONG_INKJET_PRES_BC
       SUB_STATE_BC=>YAOHONG_INKJET_STATE_BC
      else if (probtype.eq.2000) then
       SUB_INIT_MODULE=>INIT_flexible_plate_impact_MODULE
       SUB_check_vel_rigid=>flexible_plate_check_vel_rigid
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
       SUB_ASSIMILATE=>flexible_plate_impact_ASSIMILATE
      else if (probtype.eq.710) then
       SUB_INIT_MODULE=>INIT_CAVITY_PHASE_CHANGE_MODULE
       SUB_LS=>CAVITY_PHASE_CHANGE_LS
       SUB_VEL=>CAVITY_PHASE_CHANGE_VEL
       SUB_PRES=>CAVITY_PHASE_CHANGE_PRES
       SUB_STATE=>CAVITY_PHASE_CHANGE_STATE
       SUB_LS_BC=>CAVITY_PHASE_CHANGE_LS_BC
       SUB_VEL_BC=>CAVITY_PHASE_CHANGE_VEL_BC
       SUB_PRES_BC=>CAVITY_PHASE_CHANGE_PRES_BC
       SUB_STATE_BC=>CAVITY_PHASE_CHANGE_STATE_BC
       SUB_velfreestream=>CAVITY_PHASE_CHANGE_velfreestream
       SUB_nucleation=>Satomodel_nucleation
       SUB_INIT_REGIONS_LIST=>CAVITY_BOILING_INIT_REGIONS_LIST
       SUB_CHARFN_REGION=>CAVITY_BOILING_CHARFN_REGION
       SUB_THERMAL_K=>CAVITY_BOILING_THERMAL_K
      else if (probtype.eq.55) then
       SUB_INIT_MODULE=>INIT_GENERAL_PHASE_CHANGE_MODULE
       SUB_check_vel_rigid=>GENERAL_PHASE_CHANGE_check_vel_rigid
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
       SUB_ICE_SUBSTRATE_DISTANCE=>GENERAL_PHASE_CHANGE_ICE_SUBSTRATE_DISTANCE
       SUB_CFL_HELPER=>GENERAL_PHASE_CHANGE_CFL_HELPER
       SUB_SUMINT=>GENERAL_PHASE_CHANGE_SUMINT ! Nusseltt number
      else if (probtype.eq.82) then
       SUB_INIT_MODULE=>INIT_ROTATING_ANNULUS_MODULE
       SUB_LS=>ROTATING_ANNULUS_LS
       SUB_VEL=>ROTATING_ANNULUS_VEL
       SUB_PRES=>ROTATING_ANNULUS_PRES
       SUB_STATE=>ROTATING_ANNULUS_STATE
       SUB_LS_BC=>ROTATING_ANNULUS_LS_BC
       SUB_VEL_BC=>ROTATING_ANNULUS_VEL_BC
       SUB_PRES_BC=>ROTATING_ANNULUS_PRES_BC
       SUB_STATE_BC=>ROTATING_ANNULUS_STATE_BC
       SUB_INTERNAL_GRAVITY_WAVE_FLAG=> &
         ROTATING_ANNULUS_INTERNAL_GRAVITY_WAVE_FLAG
       SUB_MAPPING_WEIGHT_COEFF=>ROTATING_ANNULUS_MAPPING_WEIGHT_COEFF
       SUB_T0_Boussinesq=>ROTATING_ANNULUS_T0_Boussinesq
       SUB_V0_Coriolis=>ROTATING_ANNULUS_V0_Coriolis
      else
       ! assign null routines here that would cause the program to abort
       ! if called.  In otherwords, these are routines, that if called,
       ! MUST BE DEFINED and CANNOT depend on the STUB routines.
       ! In otherwords, these routines must be uniquely defined for each
       ! user defined problem.
       SUB_INIT_MODULE=>NULL() ! always called
       SUB_LS=>NULL() ! always called
!       SUB_clamped_LS_no_scale=>NULL()
       SUB_VEL=>NULL() ! always called
       SUB_EOS=>NULL()
       SUB_dVdT=>NULL()
       SUB_SOUNDSQR=>NULL()
       SUB_INTERNAL=>NULL()
       SUB_TEMPERATURE=>NULL()
       SUB_PRES=>NULL() ! always called
       SUB_STATE=>NULL() ! always called
       SUB_LS_BC=>NULL() ! always called
       SUB_VEL_BC=>NULL() ! always called
       SUB_PRES_BC=>NULL() ! always called
       SUB_STATE_BC=>NULL() ! always called
       SUB_HEATSOURCE=>NULL()
       SUB_EB_heat_source=>NULL()
       SUB_microcell_heat_coeff=>NULL()
       SUB_velfreestream=>NULL()
       SUB_nucleation=>NULL()
       SUB_CFL_HELPER=>NULL()
       SUB_correct_pres_rho_hydrostatic=>NULL()
      endif
     
      fort_double_size=SIZEOF(global_pressure_scale)
      fort_int_size=SIZEOF(local_dir)
     
      if ((fort_double_size.eq.8).and. &
          (fort_int_size.eq.cc_int_size)) then
       print *,"fort_double_size=",fort_double_size     
       print *,"fort_int_size=",fort_int_size     
       print *,"cc_int_size=",cc_int_size     
      else
       print *,"fort_double_size or fort_int_size invalid"
       print *,"fort_double_size=",fort_double_size     
       print *,"fort_int_size=",fort_int_size     
       print *,"cc_int_size=",cc_int_size     
       stop
      endif

      global_pressure_scale=one
      global_velocity_scale=one
     
      do local_dir=1,SDIM
       fort_n_cell(local_dir)=ccn_cell(local_dir)
       if ((fort_n_cell(local_dir).ge.4).and. &
           (4*(fort_n_cell(local_dir)/4).eq. &
            fort_n_cell(local_dir))) then
        ! do nothing
       else
        print *,"fort_n_cell must be divisible by 4"
        stop
       endif
      enddo

      if (SDIM.eq.2) then
       fort_n_cell(3)=fort_n_cell(SDIM)
      else if (SDIM.eq.3) then
       ! do nothing
      else
       print *,"dimension bust"
       stop
      endif
  
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
       ! declared in: GLOBALUTIL.F90
       ! define GG,GL weights for interpolation and derivatives
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
      
      if ((problenx.gt.zero).and. &
          (probleny.gt.zero).and. &
          (problenz.gt.zero)) then
       !do nothing
      else
       print *,"problenx or probleny or problenz invalid"
       stop
      endif
     
      problo_array(1)=problox 
      problo_array(2)=probloy 
      problo_array(3)=probloz 
      probhi_array(1)=probhix 
      probhi_array(2)=probhiy 
      probhi_array(3)=probhiz 
      do local_dir=1,3
       problen_array(local_dir)= &
          probhi_array(local_dir)-problo_array(local_dir)
       if (problen_array(local_dir).gt.zero) then
        ! do nothing
       else
        print *,"problen_array(local_dir) invalid"
        stop
       endif
      enddo ! local_dir=1,3

      fort_stop_time=ccstop_time
      
      prescribe_temperature_outflow= &
        ccprescribe_temperature_outflow
      
      fort_solidheat_flag=ccsolidheat_flag
      
      levelrz=rz_flag
      denfact=ccdenfact
      velfact=ccvelfact
      
      n_sites=ccn_sites
      nucleation_init_time=ccnucleation_init_time
      
      if (nucleation_init_time.ge.zero) then
       ! do nothing
      else
       print *,"nucleation_init_time.lt.zero"
       stop
      endif
      
      if ((n_sites.lt.0).or.(4*n_sites.gt.5000)) then
       print *,"n_sites invalid(PROB_CPP_PARMS.F90), n_sites=",n_sites
       print *,"pos_sites allocated 1..5000"
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
     
      fort_num_local_aux_grids=ccnum_local_aux_grids

      if (fort_num_local_aux_grids.gt.0) then
       if (aux_data_allocated.eq.0) then
        ALLOCATE(contain_aux(fort_num_local_aux_grids))
       else
        print *,"aux_data_allocated invalid"
        stop
       endif
      else if (fort_num_local_aux_grids.eq.0) then
       ! do nothing
      else
       print *,"fort_num_local_aux_grids invalid"
       stop
      endif

      do im=1,num_materials
      
       fort_material_type(im)=ccmaterial_type(im)

       FSI_flag(im)=ccFSI_flag(im)

       if (fort_material_type(im).eq.0) then
        ! do nothing
       else if (fort_material_type(im).eq.999) then
        if ((FSI_flag(im).eq.FSI_PRESCRIBED_PROBF90).or. &
            (FSI_flag(im).eq.FSI_PRESCRIBED_NODES).or. &
            (FSI_flag(im).eq.FSI_SHOELE_CTML)) then 
         !do nothing
        else
         print *,"FSI_flag invalid in fort_override"
         print *,"im=",im
         print *,"FSI_flag(im)=",FSI_flag(im)
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
      
      num_species_var=ccnum_species_var
      if (num_species_var.lt.MAX_NUM_SPECIES) then
       ! do nothing
      else
       print *,"num_species_var too large, increase MAX_NUM_SPECIES"
       stop
      endif
      
      num_materials_viscoelastic=ccnum_materials_viscoelastic
      
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
          (num_interfaces.lt.1).or. &
          (num_materials.gt.MAX_NUM_MATERIALS).or. &
          (num_interfaces.gt.100).or. &
          (num_materials.gt.100)) then
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
       print *,"num_materials= ",num_materials
       print *,"num_interfaces= ",num_interfaces
       print *,"numspec,num_mat_visc,MAX_NUM_MATERIALS ", &
        num_species_var,num_materials_viscoelastic, &
        MAX_NUM_MATERIALS
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

      fort_R_Palmore_Desjardins=ccR_Palmore_Desjardins
      
      do im=1,num_materials
      
       fort_DrhoDT(im)=ccDrhoDT(im)
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

       fort_viscconst_eddy_wall(im)=ccviscconst_eddy_wall(im)
       fort_viscconst_eddy_bulk(im)=ccviscconst_eddy_bulk(im)
       fort_heatviscconst_eddy_wall(im)=ccheatviscconst_eddy_wall(im)
       fort_heatviscconst_eddy_bulk(im)=ccheatviscconst_eddy_bulk(im)

       fort_thermal_microlayer_size(im)=ccthermal_microlayer_size(im)
       fort_shear_microlayer_size(im)=ccshear_microlayer_size(im)
       fort_buoyancy_microlayer_size(im)=ccbuoyancy_microlayer_size(im)
       fort_phasechange_microlayer_size(im)=ccphasechange_microlayer_size(im)

       if (fort_thermal_microlayer_size(im).le.zero) then
        print *,"fort_thermal_microlayer_size(im) must be positive"
        stop
       endif
       if (fort_shear_microlayer_size(im).le.zero) then
        print *,"fort_shear_microlayer_size(im) must be positive"
        stop
       endif
       if (fort_buoyancy_microlayer_size(im).le.zero) then
        print *,"fort_buoyancy_microlayer_size(im) must be positive"
        stop
       endif
       if (fort_phasechange_microlayer_size(im).le.zero) then
        print *,"fort_phasechange_microlayer_size(im) must be positive"
        stop
       endif

       fort_viscosity_state_model(im)= &
         ccviscosity_state_model(im)
       fort_elastic_viscosity(im)=ccelastic_viscosity(im)
       fort_elastic_time(im)=ccelastic_time(im)
       fort_viscoelastic_model(im)=ccviscoelastic_model(im)
       fort_store_elastic_data(im)=ccstore_elastic_data(im)
       fort_heatflux_factor(im)=ccheatflux_factor(im)
       fort_heatviscconst(im)=ccheatviscconst(im)
       fort_prerecalesce_heatviscconst(im)=ccprerecalesce_heatviscconst(im)
       fort_prerecalesce_viscconst(im)=ccprerecalesce_viscconst(im)
       fort_prerecalesce_stiffCP(im)=ccprerecalesce_stiffCP(im)
       fort_prerecalesce_stiffCV(im)=ccprerecalesce_stiffCV(im)
      
       fort_im_elastic_map(im)=-1
     
       fort_Carreau_alpha(im)=ccCarreau_alpha(im)
       fort_Carreau_beta(im)=ccCarreau_beta(im)
       fort_Carreau_n(im)=ccCarreau_n(im)
       fort_Carreau_mu_inf(im)=ccCarreau_mu_inf(im)

       fort_shear_thinning_fluid(im)=ccshear_thinning_fluid(im)

       fort_polymer_factor(im)=ccpolymer_factor(im)
       fort_concentration(im)=ccconcentration(im)
       fort_etaL(im)=ccetaL(im)
       fort_etaS(im)=ccetaS(im)
       fort_etaP(im)=ccetaP(im)

      enddo ! im=1..num_materials

      do local_dir=1,SDIM
       fort_grid_stretching_parameter(local_dir)= &
          ccgrid_stretching_parameter(local_dir)
      enddo

      fort_visc_coef=ccvisc_coef

      fort_angular_velocity=ccangular_velocity
      
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
      
      if (num_interfaces.ge.MAX_NUM_INTERFACES) then
       print *,"too many surface tension coefficients, increase "
       print *,"MAX_NUM_INTERFACES"
       print *,"num_interfaces= ",num_interfaces
       stop
      endif

      do iten=1,num_interfaces
       fort_latent_heat(iten)=cclatent_heat(iten)
       fort_latent_heat(num_interfaces+iten)=cclatent_heat(num_interfaces+iten)
       fort_latent_heat_slope(iten)=cclatent_heat_slope(iten)
       fort_latent_heat_slope(num_interfaces+iten)= &
               cclatent_heat_slope(num_interfaces+iten)
       fort_latent_heat_T0(iten)=cclatent_heat_T0(iten)
       fort_latent_heat_T0(num_interfaces+iten)= &
               cclatent_heat_T0(num_interfaces+iten)
       fort_latent_heat_min(iten)=cclatent_heat_min(iten)
       fort_latent_heat_min(num_interfaces+iten)= &
               cclatent_heat_min(num_interfaces+iten)

       fort_saturation_temp(iten)=ccsaturation_temp(iten)
       fort_saturation_temp(num_interfaces+iten)= &
               ccsaturation_temp(num_interfaces+iten)
      
       fort_reference_pressure(iten)=ccreference_pressure(iten)
       fort_reference_pressure(num_interfaces+iten)= &
               ccreference_pressure(num_interfaces+iten)

       fort_tension(iten)=cctension(iten)
       fort_tension_init(iten)=cctension_init(iten)
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
       if (fort_heatflux_factor(im).ge.zero) then
        ! do nothing
       else
        print *,"fort_heatflux_factor invalid"
        stop
       endif
      enddo ! do im=1,num_materials

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

       print *,"fort_R_Palmore_Desjardins ", &
         fort_R_Palmore_Desjardins

       do im=1,num_materials
        print *,"im,mat type ",im,fort_material_type(im)
        print *,"im,fort_molar_mass ",im,fort_molar_mass(im)
        print *,"im,DrhoDT ",im,fort_DrhoDT(im)
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

        print *,"im,viscconst_eddy_wall ",im,fort_viscconst_eddy_wall(im)
        print *,"im,viscconst_eddy_bulk ",im,fort_viscconst_eddy_bulk(im)
        print *,"im,heatviscconst_eddy_wall ",im, &
                fort_heatviscconst_eddy_wall(im)
        print *,"im,heatviscconst_eddy_bulk ",im, &
                fort_heatviscconst_eddy_bulk(im)

        print *,"im,thermal_microlayer_size ",im, &
                fort_thermal_microlayer_size(im)
        print *,"im,shear_microlayer_size ",im, &
                fort_shear_microlayer_size(im)
        print *,"im,buoyancy_microlayer_size ",im, &
                fort_buoyancy_microlayer_size(im)
        print *,"im,phasechange_microlayer_size ",im, &
                fort_phasechange_microlayer_size(im)

        print *,"im,viscosity_state_model ",im, &
         fort_viscosity_state_model(im)
        print *,"im,fort_elastic_viscosity ",im,fort_elastic_viscosity(im)
        print *,"im,fort_elastic_time ",im,fort_elastic_time(im)
        print *,"im,fort_viscoelastic_model ",im,fort_viscoelastic_model(im)
        print *,"im,fort_store_elastic_data ",im,fort_store_elastic_data(im)
        print *,"im,fort_im_elastic_map ",im,fort_im_elastic_map(im)

        print *,"im,fort_Carreau_alpha ",im,fort_Carreau_alpha(im)
        print *,"im,fort_Carreau_beta ",im,fort_Carreau_beta(im)
        print *,"im,fort_Carreau_n ",im,fort_Carreau_n(im)
        print *,"im,fort_Carreau_mu_inf ",im,fort_Carreau_mu_inf(im)

        print *,"im,fort_shear_thinning_fluid ", &
                im,fort_shear_thinning_fluid(im)
      
        print *,"im,fort_polymer_factor ",im,fort_polymer_factor(im)
        print *,"im,fort_concentration ",im,fort_concentration(im)

        print *,"im,fort_etaL ",im,fort_etaL(im)
        print *,"im,fort_etaS ",im,fort_etaS(im)
        print *,"im,fort_etaP ",im,fort_etaP(im)

        print *,"im,fort_heatflux_factor ",im,fort_heatflux_factor(im)
       
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

       do local_dir=1,SDIM
        print *,"local_dir,fort_grid_stretching_parameter ",local_dir, &
           fort_grid_stretching_parameter(local_dir)
       enddo

       do im=1,num_species_var*num_materials
        print *,"im,species ",im,fort_speciesconst(im)
        print *,"im,speciesvisc ",im,fort_speciesviscconst(im)
       enddo
       do im=1,num_species_var
        print *,"im,fort_species_molar_mass ",im,fort_species_molar_mass(im)
       enddo
     
       print *,"fort_visc_coef= ",fort_visc_coef

       print *,"fort_angular_velocity= ",fort_angular_velocity

       do iten=1,num_interfaces
        print *,"iten,tension ",iten,fort_tension(iten)
        print *,"iten,tension_init ",iten,fort_tension_init(iten)
        print *,"iten,tension_slope ",iten,fort_tension_slope(iten)
        if (fort_tension_slope(iten).le.zero) then
         ! do nothing
        else
         print *,"fort_tension_slope must be non-positive"
         stop
        endif
        print *,"iten,tension_T0 ",iten,fort_tension_T0(iten)
        print *,"iten,tension_min ",iten,fort_tension_min(iten)
        print *,"iten,prefreeze_tension ",iten,fort_prefreeze_tension(iten)
        print *,"iten,fort_saturation_temp ",iten,fort_saturation_temp(iten)
        print *,"iten,fort_reference_pressure ",iten, &
          fort_reference_pressure(iten)
        print *,"iten+num_interfaces,fort_saturation_temp ", &
          iten+num_interfaces, &
          fort_saturation_temp(iten+num_interfaces)
        print *,"iten+num_interfaces,fort_reference_pressure ", &
          iten+num_interfaces, &
          fort_reference_pressure(iten+num_interfaces)
       enddo ! iten=1..num_interfaces
      endif
     
      if (ioproc.eq.1) then
       do iten=1,2*num_interfaces
        print *,"iten,fort_latent_heat ",iten,fort_latent_heat(iten)
        print *,"iten,fort_latent_heat_slope ",iten,fort_latent_heat_slope(iten)
        print *,"iten,fort_latent_heat_T0 ",iten,fort_latent_heat_T0(iten)
        print *,"iten,fort_latent_heat_min ",iten,fort_latent_heat_min(iten)
        if (fort_latent_heat_slope(iten).le.zero) then
         ! do nothing
        else
         print *,"fort_latent_heat_slope must be non-positive"
         stop
        endif
       enddo ! iten=1..2*num_interfaces
      endif

      do local_dir=1,SDIM
       gravity_vector(local_dir)=ccgravity_vector(local_dir)
      enddo
      
       ! declared in: GLOBALUTIL.F90
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
      
       ! above: probtype.eq.110
      else if ((probtype.eq.1).and. &
               ((axis_dir.eq.150).or. &
                (axis_dir.eq.151))) then
      
       call shockdrop_init()
      
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
      
       ! in: subroutine fort_override
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
     
       ! initialize grid mapping variables here.
       ! mapping_n_cell=n_cell * 2^max_level
       !  index,dir
       ! REAL_T, allocatable, dimension(:,:) :: mapping_comp_to_phys
       ! REAL_T, allocatable, dimension(:,:) :: mapping_phys_to_comp
       ! INTEGER_T :: mapping_n_cell(0:2)
       ! INTEGER_T :: mapping_allocated=0

      mapping_allocated=1
      mapping_n_cell_max=0
      do local_dir=0,SDIM-1
       mapping_n_cell(local_dir)=fort_n_cell(local_dir+1)
       do level=1,fort_max_level
        mapping_n_cell(local_dir)=2*mapping_n_cell(local_dir)
       enddo
       mapping_n_cell_max=max(mapping_n_cell_max,mapping_n_cell(local_dir))
      enddo ! local_dir=0,SDIM-1
      allocate(mapping_comp_to_phys(0:mapping_n_cell_max,0:2))
      allocate(mapping_phys_to_comp(0:mapping_n_cell_max,0:2))
      do local_dir=0,SDIM-1
       call single_dimension_grid_mapping(local_dir)
      enddo

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
      
       if ((ls_homflag.lt.0).or.(ls_homflag.gt.1).or. &
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
      
       print *,"fort:gravity(1,..,sdim) ", &
        gravity_vector(1), &
        gravity_vector(2), &
        gravity_vector(SDIM)
       print *,"fort:pres,vel,temp,spec,ls homflag ",pres_homflag, &
        vel_homflag,temp_homflag,species_homflag,ls_homflag
      
       print *,"fort: stop_time ",fort_stop_time
      
       print *,"fort: rz ",levelrz
       do im=1,num_materials
        print *,"fort: im,FSI_flag ",im,FSI_flag(im)
       enddo

       print *,"fort: fort_num_local_aux_grids= ",fort_num_local_aux_grids

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
      end subroutine fort_override


      subroutine fort_override_MAIN_GLOBALS( &
        cc_int_size, &
        ccnum_species_var, &
        ccnum_materials_viscoelastic, &
        ccnum_state_material, &
        ccnum_state_base, &
        ccngeom_raw, &
        ccngeom_recon, &
        ccnum_materials, &
        ccnten, &
        ioproc) &
      bind(c,name='fort_override_MAIN_GLOBALS')

      use probcommon_module
      
      IMPLICIT NONE
      
      INTEGER_T, INTENT(in) :: cc_int_size
      INTEGER_T, INTENT(in) :: ccnum_materials
      INTEGER_T, INTENT(in) :: ccnten
      INTEGER_T, INTENT(in) :: ccnum_species_var
      INTEGER_T, INTENT(in) :: ccnum_materials_viscoelastic
      INTEGER_T, INTENT(in) :: ccnum_state_material
      INTEGER_T, INTENT(in) :: ccnum_state_base
      INTEGER_T, INTENT(in) :: ccngeom_raw
      INTEGER_T, INTENT(in) :: ccngeom_recon
      INTEGER_T, INTENT(in) :: ioproc

      INTEGER_T :: local_var_int
      REAL_T :: local_var_double
      INTEGER_T :: fort_double_size,fort_int_size
      
      num_materials=ccnum_materials
      if (num_materials.lt.MAX_NUM_MATERIALS) then
       ! do nothing
      else
       print *,"increase MAX_NUM_MATERIALS; aborting"
       stop
      endif
      num_interfaces=( (num_materials-1)*(num_materials-1)+ &
             (num_materials-1) )/2
      if (num_interfaces.eq.ccnten) then
       ! do nothing
      else
       print *,"num_interfaces<>ccnten"
       stop
      endif

      fort_double_size=SIZEOF(local_var_double)
      fort_int_size=SIZEOF(local_var_int)
     
      if ((fort_double_size.eq.8).and. &
          (fort_int_size.eq.cc_int_size)) then
       print *,"fort_override_MAIN_GLOBALS"
       print *,"fort_double_size=",fort_double_size     
       print *,"fort_int_size=",fort_int_size     
       print *,"cc_int_size=",cc_int_size     
      else
       print *,"fort_override_MAIN_GLOBALS"
       print *,"fort_double_size or fort_int_size invalid"
       print *,"fort_double_size=",fort_double_size     
       print *,"fort_int_size=",fort_int_size     
       print *,"cc_int_size=",cc_int_size     
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
          (num_interfaces.lt.1).or. &
          (num_materials.gt.MAX_NUM_MATERIALS).or. &
          (num_interfaces.gt.100).or. &
          (num_materials.gt.100)) then
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
       print *,"fort_override_MAIN_GLOBALS"
      
       print *,"num_materials= ",num_materials
       print *,"num_interfaces= ",num_interfaces
       print *,"numspec,num_mat_visc,MAX_NUM_MATERIALS ", &
        num_species_var,num_materials_viscoelastic, &
        MAX_NUM_MATERIALS
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
      
      return
      end subroutine fort_override_MAIN_GLOBALS
