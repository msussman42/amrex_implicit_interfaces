      INTEGER_T fort_max_level
      INTEGER_T fort_n_cell(3)
      INTEGER_T bfact_space_order(0:100)
      INTEGER_T bfact_time_order
 
      REAL_T adv_vel
      REAL_T rgasinlet
      REAL_T vinletgas
      REAL_T denfact
      REAL_T xblobarr(10)
      REAL_T yblobarr(10)
      REAL_T zblobarr(10)
      REAL_T radblobarr(10)
      REAL_T xblob,yblob,zblob,radblob
      REAL_T xblob2,yblob2,zblob2,radblob2
      REAL_T xblob3,yblob3,zblob3,radblob3
      REAL_T xblob4,yblob4,zblob4,radblob4
      REAL_T xblob5,yblob5,zblob5,radblob5
      REAL_T xblob6,yblob6,zblob6,radblob6
      REAL_T xblob7,yblob7,zblob7,radblob7
      REAL_T xblob8,yblob8,zblob8,radblob8
      REAL_T xblob9,yblob9,zblob9,radblob9
      REAL_T xblob10,yblob10,zblob10,radblob10
      REAL_T xactive,yactive,zactive
      REAL_T ractivex,ractivey,ractivez
      REAL_T velfact
      REAL_T advbot
      REAL_T inflow_pressure,outflow_pressure
      REAL_T period_time
      REAL_T twall
      REAL_T problox,probloy,probloz
      REAL_T probhix,probhiy,probhiz
      REAL_T problenx,probleny,problenz
      REAL_T problo_array(3)
      REAL_T probhi_array(3)
      REAL_T problen_array(3)
      INTEGER_T fort_is_periodic(3)
      REAL_T fort_stop_time
      REAL_T fort_energyconst(MAX_NUM_MATERIALS)
      REAL_T fort_energycutoff(MAX_NUM_MATERIALS)
      REAL_T fort_DrhoDT(MAX_NUM_MATERIALS)
      REAL_T fort_tempconst(MAX_NUM_MATERIALS)
      REAL_T fort_initial_temperature(MAX_NUM_MATERIALS)
      REAL_T fort_tempcutoff(MAX_NUM_MATERIALS)
      REAL_T fort_tempcutoffmax(MAX_NUM_MATERIALS)
      REAL_T fort_denconst(MAX_NUM_MATERIALS)
      REAL_T fort_density_floor(MAX_NUM_MATERIALS)
      REAL_T fort_density_ceiling(MAX_NUM_MATERIALS)
      REAL_T fort_cavdenconst(MAX_NUM_MATERIALS)
      REAL_T fort_viscconst(MAX_NUM_MATERIALS)

      REAL_T fort_viscconst_eddy_wall(MAX_NUM_MATERIALS)
      REAL_T fort_viscconst_eddy_bulk(MAX_NUM_MATERIALS)
      REAL_T fort_heatviscconst_eddy_wall(MAX_NUM_MATERIALS)
      REAL_T fort_heatviscconst_eddy_bulk(MAX_NUM_MATERIALS)

      REAL_T fort_thermal_microlayer_size(MAX_NUM_MATERIALS)
      REAL_T fort_shear_microlayer_size(MAX_NUM_MATERIALS)
      REAL_T fort_buoyancy_microlayer_size(MAX_NUM_MATERIALS)
      REAL_T fort_phasechange_microlayer_size(MAX_NUM_MATERIALS)

      INTEGER_T fort_viscosity_state_model(MAX_NUM_MATERIALS)
      REAL_T fort_heatflux_factor(MAX_NUM_MATERIALS)
      REAL_T fort_heatviscconst(MAX_NUM_MATERIALS)

      REAL_T fort_prerecalesce_heatviscconst(MAX_NUM_MATERIALS)
      REAL_T fort_prerecalesce_viscconst(MAX_NUM_MATERIALS)
      REAL_T fort_prerecalesce_stiffCP(MAX_NUM_MATERIALS)
      REAL_T fort_prerecalesce_stiffCV(MAX_NUM_MATERIALS)
      
      REAL_T fort_speciesconst(MAX_NUM_SPECIES*MAX_NUM_MATERIALS)
      REAL_T fort_speciesviscconst(MAX_NUM_SPECIES*MAX_NUM_MATERIALS)

      REAL_T fort_static_surface_tension_duration

      REAL_T fort_static_tension(MAX_NUM_INTERFACES)
      REAL_T fort_tension(MAX_NUM_INTERFACES)
      REAL_T fort_tension_init(MAX_NUM_INTERFACES)
      REAL_T fort_tension_slope(MAX_NUM_INTERFACES)
      REAL_T fort_tension_T0(MAX_NUM_INTERFACES)
      REAL_T fort_tension_min(MAX_NUM_INTERFACES)
      REAL_T fort_prefreeze_tension(MAX_NUM_INTERFACES)

      REAL_T fort_latent_heat(2*MAX_NUM_INTERFACES)
      REAL_T fort_latent_heat_slope(2*MAX_NUM_INTERFACES)
      REAL_T fort_latent_heat_T0(2*MAX_NUM_INTERFACES)
      REAL_T fort_latent_heat_min(2*MAX_NUM_INTERFACES)

      REAL_T fort_saturation_temp(2*MAX_NUM_INTERFACES)
      REAL_T fort_reference_pressure(2*MAX_NUM_INTERFACES)

      REAL_T fort_molar_mass(MAX_NUM_MATERIALS)
      REAL_T fort_species_molar_mass(MAX_NUM_SPECIES)

      REAL_T fort_alpha(2*MAX_NUM_INTERFACES)
      REAL_T fort_beta(2*MAX_NUM_INTERFACES)
      REAL_T fort_expansion_factor(2*MAX_NUM_INTERFACES)
      REAL_T fort_time_radblob(2*MAX_NUM_INTERFACES)
      REAL_T fort_stefan_number(2*MAX_NUM_INTERFACES)
      REAL_T fort_jacob_number(2*MAX_NUM_INTERFACES)

      REAL_T fort_stiffPINF(MAX_NUM_MATERIALS)

      REAL_T fort_R_Palmore_Desjardins

      REAL_T fort_stiffCP(MAX_NUM_MATERIALS)
      REAL_T fort_stiffCV(MAX_NUM_MATERIALS)
      REAL_T fort_stiffGAMMA(MAX_NUM_MATERIALS)

      REAL_T density_at_depth
      REAL_T global_pressure_scale
      REAL_T global_velocity_scale
      REAL_T gravity_vector(AMREX_SPACEDIM)
      REAL_T gravity_boussinesq_vector(AMREX_SPACEDIM)

      INTEGER_T num_species_var
      INTEGER_T num_materials_viscoelastic
      INTEGER_T fort_im_elastic_map(MAX_NUM_MATERIALS)
      INTEGER_T fort_viscoelastic_model(MAX_NUM_MATERIALS)
      REAL_T fort_elastic_viscosity(MAX_NUM_MATERIALS)
      REAL_T fort_elastic_time(MAX_NUM_MATERIALS)
      INTEGER_T fort_store_elastic_data(MAX_NUM_MATERIALS)

      REAL_T fort_Carreau_alpha(MAX_NUM_MATERIALS)
      REAL_T fort_Carreau_beta(MAX_NUM_MATERIALS)
      REAL_T fort_Carreau_n(MAX_NUM_MATERIALS)
      REAL_T fort_Carreau_mu_inf(MAX_NUM_MATERIALS)
      INTEGER_T fort_shear_thinning_fluid(MAX_NUM_MATERIALS)
      REAL_T fort_polymer_factor(MAX_NUM_MATERIALS)
      REAL_T fort_concentration(MAX_NUM_MATERIALS)
      REAL_T fort_etaL(MAX_NUM_MATERIALS)
      REAL_T fort_etaS(MAX_NUM_MATERIALS)
      REAL_T fort_etaP(MAX_NUM_MATERIALS)

      REAL_T fort_visc_coef

      INTEGER_T num_state_material
      INTEGER_T num_state_base
      INTEGER_T ngeom_raw
      INTEGER_T ngeom_recon
      INTEGER_T num_materials
      INTEGER_T num_interfaces
      INTEGER_T probtype, adv_dir, axis_dir
      INTEGER_T levelrz
      INTEGER_T prescribe_temperature_outflow
      INTEGER_T fort_solidheat_flag
      INTEGER_T FSI_flag(MAX_NUM_MATERIALS)
      INTEGER_T fort_material_type(MAX_NUM_MATERIALS)
      INTEGER_T fort_ZEYU_DCA_SELECT
      INTEGER_T pres_homflag
      INTEGER_T vel_homflag
      INTEGER_T temp_homflag
      INTEGER_T species_homflag
      INTEGER_T ls_homflag
      INTEGER_T fort_finest_level
      INTEGER_T n_sites
      REAL_T nucleation_init_time
      REAL_T pos_sites(5000)

