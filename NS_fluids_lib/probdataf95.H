      integer fort_max_level
      integer fort_n_cell(3)
      integer bfact_space_order(0:100)
      integer bfact_time_order
 
      real(amrex_real) adv_vel
      real(amrex_real) rgasinlet
      real(amrex_real) vinletgas
      real(amrex_real) denfact
      real(amrex_real) xblobarr(10)
      real(amrex_real) yblobarr(10)
      real(amrex_real) zblobarr(10)
      real(amrex_real) radblobarr(10)
      real(amrex_real) xblob,yblob,zblob,radblob
      real(amrex_real) xblob2,yblob2,zblob2,radblob2
      real(amrex_real) xblob3,yblob3,zblob3,radblob3
      real(amrex_real) xblob4,yblob4,zblob4,radblob4
      real(amrex_real) xblob5,yblob5,zblob5,radblob5
      real(amrex_real) xblob6,yblob6,zblob6,radblob6
      real(amrex_real) xblob7,yblob7,zblob7,radblob7
      real(amrex_real) xblob8,yblob8,zblob8,radblob8
      real(amrex_real) xblob9,yblob9,zblob9,radblob9
      real(amrex_real) xblob10,yblob10,zblob10,radblob10
      real(amrex_real) xactive,yactive,zactive
      real(amrex_real) ractivex,ractivey,ractivez
      real(amrex_real) velfact
      real(amrex_real) advbot
      real(amrex_real) inflow_pressure,outflow_pressure
      real(amrex_real) period_time
      real(amrex_real) twall
      real(amrex_real) problox,probloy,probloz
      real(amrex_real) probhix,probhiy,probhiz
      real(amrex_real) problenx,probleny,problenz
      real(amrex_real) problo_array(3)
      real(amrex_real) probhi_array(3)
      real(amrex_real) problen_array(3)
      integer fort_is_periodic(3)
      real(amrex_real) fort_stop_time
      real(amrex_real) fort_energyconst(MAX_NUM_MATERIALS)
      real(amrex_real) fort_energycutoff(MAX_NUM_MATERIALS)
      real(amrex_real) fort_DrhoDT(MAX_NUM_MATERIALS)
      real(amrex_real) fort_tempconst(MAX_NUM_MATERIALS)
      real(amrex_real) fort_initial_temperature(MAX_NUM_MATERIALS)
      real(amrex_real) fort_tempcutoff(MAX_NUM_MATERIALS)
      real(amrex_real) fort_tempcutoffmax(MAX_NUM_MATERIALS)
      real(amrex_real) fort_denconst(MAX_NUM_MATERIALS)
      real(amrex_real) fort_density_floor(MAX_NUM_MATERIALS)
      real(amrex_real) fort_density_ceiling(MAX_NUM_MATERIALS)
      real(amrex_real) fort_cavdenconst(MAX_NUM_MATERIALS)
      real(amrex_real) fort_viscconst(MAX_NUM_MATERIALS)

      real(amrex_real) fort_viscconst_eddy_wall(MAX_NUM_MATERIALS)
      real(amrex_real) fort_viscconst_eddy_bulk(MAX_NUM_MATERIALS)
      real(amrex_real) fort_heatviscconst_eddy_wall(MAX_NUM_MATERIALS)
      real(amrex_real) fort_heatviscconst_eddy_bulk(MAX_NUM_MATERIALS)

      real(amrex_real) fort_thermal_microlayer_size(MAX_NUM_MATERIALS)
      real(amrex_real) fort_shear_microlayer_size(MAX_NUM_MATERIALS)
      real(amrex_real) fort_buoyancy_microlayer_size(MAX_NUM_MATERIALS)
      real(amrex_real) fort_phasechange_microlayer_size(MAX_NUM_MATERIALS)

      integer fort_viscosity_state_model(MAX_NUM_MATERIALS)
      real(amrex_real) fort_heatflux_factor(MAX_NUM_MATERIALS)
      real(amrex_real) fort_heatviscconst(MAX_NUM_MATERIALS)

      real(amrex_real) fort_prerecalesce_heatviscconst(MAX_NUM_MATERIALS)
      real(amrex_real) fort_prerecalesce_viscconst(MAX_NUM_MATERIALS)
      real(amrex_real) fort_prerecalesce_stiffCP(MAX_NUM_MATERIALS)
      real(amrex_real) fort_prerecalesce_stiffCV(MAX_NUM_MATERIALS)
      
      real(amrex_real) fort_speciesconst(MAX_NUM_SPECIES*MAX_NUM_MATERIALS)
      real(amrex_real) fort_speciesviscconst(MAX_NUM_SPECIES*MAX_NUM_MATERIALS)

      real(amrex_real) fort_tension(MAX_NUM_INTERFACES)
      real(amrex_real) fort_tension_init(MAX_NUM_INTERFACES)
      real(amrex_real) fort_tension_slope(MAX_NUM_INTERFACES)
      real(amrex_real) fort_tension_T0(MAX_NUM_INTERFACES)
      real(amrex_real) fort_tension_min(MAX_NUM_INTERFACES)
      real(amrex_real) fort_prefreeze_tension(MAX_NUM_INTERFACES)

      real(amrex_real) fort_latent_heat(2*MAX_NUM_INTERFACES)
      real(amrex_real) fort_latent_heat_slope(2*MAX_NUM_INTERFACES)
      real(amrex_real) fort_latent_heat_T0(2*MAX_NUM_INTERFACES)
      real(amrex_real) fort_latent_heat_min(2*MAX_NUM_INTERFACES)

      real(amrex_real) fort_saturation_temp(2*MAX_NUM_INTERFACES)
      real(amrex_real) fort_reference_pressure(2*MAX_NUM_INTERFACES)

      real(amrex_real) fort_molar_mass(MAX_NUM_MATERIALS)
      real(amrex_real) fort_species_molar_mass(MAX_NUM_SPECIES)

      real(amrex_real) fort_alpha(2*MAX_NUM_INTERFACES)
      real(amrex_real) fort_beta(2*MAX_NUM_INTERFACES)
      real(amrex_real) fort_expansion_factor(2*MAX_NUM_INTERFACES)
      real(amrex_real) fort_time_radblob(2*MAX_NUM_INTERFACES)
      real(amrex_real) fort_stefan_number(2*MAX_NUM_INTERFACES)
      real(amrex_real) fort_jacob_number(2*MAX_NUM_INTERFACES)

      real(amrex_real) fort_stiffPINF(MAX_NUM_MATERIALS)

      real(amrex_real) fort_R_Palmore_Desjardins

      real(amrex_real) fort_stiffCP(MAX_NUM_MATERIALS)
      real(amrex_real) fort_stiffCV(MAX_NUM_MATERIALS)
      real(amrex_real) fort_stiffGAMMA(MAX_NUM_MATERIALS)

      real(amrex_real) density_at_depth
      real(amrex_real) global_pressure_scale
      real(amrex_real) global_velocity_scale
      real(amrex_real) gravity_vector(AMREX_SPACEDIM)

      integer num_species_var
      integer num_materials_viscoelastic
      integer num_materials_compressible
      integer fort_im_elastic_map(MAX_NUM_MATERIALS)
      integer fort_im_refine_density_map(MAX_NUM_MATERIALS)
      integer fort_viscoelastic_model(MAX_NUM_MATERIALS)
      real(amrex_real) fort_elastic_viscosity(MAX_NUM_MATERIALS)
      real(amrex_real) fort_elastic_time(MAX_NUM_MATERIALS)
      integer fort_store_elastic_data(MAX_NUM_MATERIALS)
      integer fort_store_refine_density_data(MAX_NUM_MATERIALS)

      real(amrex_real) fort_Carreau_alpha(MAX_NUM_MATERIALS)
      real(amrex_real) fort_Carreau_beta(MAX_NUM_MATERIALS)
      real(amrex_real) fort_Carreau_n(MAX_NUM_MATERIALS)
      real(amrex_real) fort_Carreau_mu_inf(MAX_NUM_MATERIALS)
      integer fort_shear_thinning_fluid(MAX_NUM_MATERIALS)
      real(amrex_real) fort_polymer_factor(MAX_NUM_MATERIALS)
      real(amrex_real) fort_concentration(MAX_NUM_MATERIALS)
      real(amrex_real) fort_etaL(MAX_NUM_MATERIALS)
      real(amrex_real) fort_etaS(MAX_NUM_MATERIALS)
      real(amrex_real) fort_etaP(MAX_NUM_MATERIALS)

      real(amrex_real) fort_visc_coef
      real(amrex_real) fort_angular_velocity

      real(amrex_real) fort_grid_stretching_parameter(SDIM)

      integer num_state_material
      integer num_state_base
      integer ngeom_raw
      integer ngeom_recon
      integer num_materials
      integer num_interfaces
      integer probtype, adv_dir, axis_dir
      integer levelrz
      integer prescribe_temperature_outflow
      integer fort_solidheat_flag
      integer FSI_flag(MAX_NUM_MATERIALS)
      integer fort_material_type(MAX_NUM_MATERIALS)
      integer fort_material_type_interface(MAX_NUM_INTERFACES)
      integer fort_conserve_total_energy
      integer fort_ZEYU_DCA_SELECT
      integer pres_homflag
      integer vel_homflag
      integer temp_homflag
      integer species_homflag
      integer ls_homflag
      integer fort_finest_level
      integer n_sites
      real(amrex_real) nucleation_init_time
      real(amrex_real) pos_sites(5000)

