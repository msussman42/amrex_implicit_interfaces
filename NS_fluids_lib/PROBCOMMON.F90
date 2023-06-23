#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_FORT_INTEGER.H"
#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_ArrayLim.H"

#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else  
print *,"dimension bust"
stop
#endif


module probcommon_module_types

      use NeuralNetwork
      use DecisionTree
      use RandomForest

 
      type law_of_wall_parm_type
      INTEGER_T :: level
      INTEGER_T :: finest_level
      INTEGER_T :: bfact
      REAL_T :: visc_coef
      REAL_T :: time
      REAL_T :: dt
      REAL_T, pointer :: usolid_raster(:)
      REAL_T, pointer :: n_raster(:) ! points to solid
      REAL_T, pointer :: x_image_raster(:)
      REAL_T, pointer :: x_probe_raster(:)
      REAL_T, pointer :: x_projection_raster(:)
      REAL_T, pointer :: dx(:)
      REAL_T :: dxmin
      REAL_T, pointer :: xlo(:)
      INTEGER_T, pointer :: fablo(:)
      INTEGER_T, pointer :: fabhi(:)
      REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: LSCP
      REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: LSFD
      REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: state ! nden comp.
      REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: ufluid
      REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: usolid
      end type law_of_wall_parm_type

      type user_defined_sum_int_type
       INTEGER_T ncomp_sum_int_user1
       INTEGER_T ncomp_sum_int_user2
       INTEGER_T ncomp_sum_int_user12
       REAL_T, pointer :: problo(:)      
       REAL_T, pointer :: probhi(:) 
       INTEGER_T :: igrid,jgrid,kgrid
       REAL_T :: volgrid
       INTEGER_T :: nhalf
       INTEGER_T :: bfact
       INTEGER_T :: den_ncomp
       INTEGER_T :: level
       INTEGER_T :: finest_level
       INTEGER_T, pointer :: tilelo(:)
       INTEGER_T, pointer :: tilehi(:)
       INTEGER_T, pointer :: fablo(:)
       INTEGER_T, pointer :: fabhi(:)
       REAL_T, pointer :: xlo(:)
       REAL_T, pointer :: dx(:)
       REAL_T, pointer :: xsten(:,:)
       REAL_T :: time
       REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: cellten
       REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: lsfab
        ! 1..num_materials*ngeom_recon
       REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: slopes
        ! num_materials * num_state_material
        ! density1,temperature1,species1_1,...,species_N_1
        ! density2,temperature2,species1_2,...,species_N_2
        ! density3,temperature3,species1_3,...,species_N_3
       REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: den
       REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: vel
       REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: visco
      end type user_defined_sum_int_type

      type nucleation_parm_type_input
       INTEGER_T :: tid
       INTEGER_T :: local_freezing_model
       REAL_T :: LL
       INTEGER_T :: i,j,k
       INTEGER_T :: im_source
       INTEGER_T :: im_dest
       REAL_T :: dxmaxLS
       INTEGER_T :: bfact
       INTEGER_T :: level
       INTEGER_T :: finest_level
       REAL_T, pointer :: dx(:)
       REAL_T, pointer :: xlo(:)
       INTEGER_T :: nstate
       INTEGER_T, pointer :: fablo(:)
       INTEGER_T, pointer :: fabhi(:)
       REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: EOS
       REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: Snew
       REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: LSnew
       REAL_T, pointer, dimension(D_DECL(:,:,:)) :: pres
       REAL_T, pointer, dimension(D_DECL(:,:,:)) :: pres_eos
       INTEGER_T :: custom_nucleation_model
       INTEGER_T :: do_the_nucleate
       INTEGER_T :: nucleate_pos_size
       REAL_T, pointer :: nucleate_pos(:)
       REAL_T, pointer :: nucleation_temp(:)
       REAL_T, pointer :: nucleation_pressure(:)
       REAL_T, pointer :: nucleation_pmg(:)
       REAL_T, pointer :: nucleation_mach(:)
       REAL_T, pointer :: cavitation_pressure(:)
       REAL_T, pointer :: cavitation_vapor_density(:)
       REAL_T, pointer :: cavitation_tension(:)
       REAL_T :: local_TSAT
       REAL_T :: prev_time
       REAL_T :: cur_time
       REAL_T :: dt
      end type nucleation_parm_type_input


      type assimilate_parm_type
      INTEGER_T :: level
      INTEGER_T :: finest_level
      INTEGER_T :: bfact
      INTEGER_T :: nparts_ghost
      INTEGER_T :: nparts
      INTEGER_T, pointer :: im_solid_map(:)
      INTEGER_T :: nstate
      INTEGER_T :: nhalf
      REAL_T :: cur_time
      REAL_T :: dt
      REAL_T, pointer :: dx(:)
      REAL_T, pointer :: xsten(:,:)
      REAL_T :: dxmin
      REAL_T, pointer :: xlo(:)
      INTEGER_T, pointer :: fablo(:)
      INTEGER_T, pointer :: fabhi(:)
      REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: ughostx
      REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: ughosty
      REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: ughostz
      end type assimilate_parm_type

      type assimilate_out_parm_type
      REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: state ! nstate comp.
      REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: LS_state
      REAL_T, pointer, dimension(D_DECL(:,:,:)) :: macx
      REAL_T, pointer, dimension(D_DECL(:,:,:)) :: macy
      REAL_T, pointer, dimension(D_DECL(:,:,:)) :: macz
      end type assimilate_out_parm_type

       ! used by deriv_from_grid_util
      type deriv_from_grid_parm_type
      INTEGER_T :: scomp
      INTEGER_T :: ncomp
      INTEGER_T :: level
      INTEGER_T :: finest_level
      INTEGER_T :: bfact
      INTEGER_T :: index_flux(SDIM)  !flux point where derivative needed.
      INTEGER_T :: dir_deriv !dir_deriv=1..sdim or dir_deriv=-1
      INTEGER_T :: box_type_flux(SDIM) !0=CELL 1=NODE
      INTEGER_T :: box_type_data(SDIM)
      INTEGER_T :: grid_type_flux  ! -1..5
      INTEGER_T :: grid_type_data  ! -1..5
      REAL_T :: dx(SDIM)
      REAL_T :: xlo(SDIM)
      INTEGER_T :: fablo(SDIM)
      INTEGER_T :: fabhi(SDIM)
!      REAL_T, INTENT(in), pointer, dimension(D_DECL(:,:,:),:) :: disp_dataptr
      end type deriv_from_grid_parm_type

       ! used by single_deriv_from_grid_util
      type single_deriv_from_grid_parm_type
      INTEGER_T :: level
      INTEGER_T :: finest_level
      INTEGER_T :: bfact
      INTEGER_T :: index_flux(SDIM)  !flux point where derivative needed.
      INTEGER_T :: dir_deriv !dir_deriv=1..sdim or dir_deriv=-1
      INTEGER_T :: box_type_flux(SDIM) !0=CELL 1=NODE
      INTEGER_T :: box_type_data(SDIM)
      INTEGER_T :: grid_type_flux  ! -1..5
      INTEGER_T :: grid_type_data  ! -1..5
      REAL_T :: dx(SDIM)
      REAL_T :: xlo(SDIM)
      INTEGER_T :: fablo(SDIM)
      INTEGER_T :: fabhi(SDIM)
!      REAL_T, INTENT(in), pointer, dimension(D_DECL(:,:,:)) :: disp_dataptr
      end type single_deriv_from_grid_parm_type

       ! used by interp_from_grid_util
      type interp_from_grid_parm_type
      INTEGER_T :: scomp
      INTEGER_T :: ncomp
      INTEGER_T :: level
      INTEGER_T :: finest_level
      INTEGER_T :: bfact
      REAL_T :: xtarget(SDIM)
      REAL_T :: dx(SDIM)
      REAL_T :: xlo(SDIM)
      INTEGER_T :: fablo(SDIM)
      INTEGER_T :: fabhi(SDIM)
!      REAL_T, INTENT(in), pointer, dimension(D_DECL(:,:,:),:) :: stateptr 
      end type interp_from_grid_parm_type


       ! used by single_interp_from_grid_util
      type single_interp_from_grid_parm_type
      INTEGER_T :: level
      INTEGER_T :: finest_level
      INTEGER_T :: bfact
      REAL_T :: xtarget(SDIM)
      REAL_T :: dx(SDIM)
      REAL_T :: xlo(SDIM)
      INTEGER_T :: fablo(SDIM)
      INTEGER_T :: fabhi(SDIM)
!      REAL_T, INTENT(in), pointer, dimension(D_DECL(:,:,:)) :: stateptr 
      end type single_interp_from_grid_parm_type


      type interp_from_grid_out_parm_type
      REAL_T, pointer :: data_interp(:)
      end type interp_from_grid_out_parm_type

     contains

end module probcommon_module_types

module probcommon_module
use probcommon_module_types

implicit none

! fort_cavdenconst added: December 25, 2018
! fort_elastic_viscosity added: August 22, 2018
! fort_im_elastic_map added: August 22, 2018
! num_materials_elastic deleted: August 21, 2018
! fort_viscconst_eddy_wall added: July 15, 2018
! fort_initial_temperature added: April 10, 2018
! nucleation_init_time added: May 5, 2018
! fort_density_floor and fort_density_ceiling added: January 8, 2019
! fort_ZEYU_DCA_SELECT added: January 20, 2020
! fort_molar_mass added: September 5, 2020.
! fort_max_num_materials and fort_max_num_eos deleted: September 6, 2020.
! fort_stiffCV added September 29.
! fort_store_elastic_data added October 3.
! fort_solidheat_flag added October 13.
! fort_lame_coefficient added October 21.
! fort_shear_modulus added October 28.
! fort_linear_elastic_model added November 4.
! fort_elastic_time added December 30
! fort_viscoelastic_model added December 30
! fort_R_Palmore_Desjardins added January 30
! fort_n_cell added March 11.
! pos_sites has 5000 components, May 8.
! num_materials_vel deleted, June 7.
! num_materials_scalar_solve deleted, June 9.
! fort_reference_pressure added, July 8.
! fort_drhodz deleted, Aug 2.
! fort_drhodt renamed to fort_DrhoDT, Aug 2.
!
! Added November 4,2021:
! fort_Carreau_alpha,
! fort_Carreau_beta,
! fort_Carreau_n,
! fort_Carreau_mu_inf,
! fort_shear_thinning_fluid,
! fort_polymer_factor,
! fort_concentration,
! fort_etaL,
! fort_etaS,
! fort_etaP,
! fort_visc_coef
! Added November 24,2021:
!  fort_viscconst_eddy_wall 
!  fort_viscconst_eddy_bulk
!  fort_heatviscconst_eddy_wall 
!  fort_heatviscconst_eddy_bulk
!  fort_thermal_microlayer_size 
!  fort_shear_microlayer_size 
!  fort_buoyancy_microlayer_size 
!  fort_phasechange_microlayer_size 
! Added November 30,2021:
!  fort_heatflux_factor
! Added December 16,2021:
!  num_interfaces
! Added February 12,2022:
!  problo_array,probhi_array,problen_array
! Added April 14, 2022:
!  fort_damping_coefficient
! added June 17, 2022:
!  fort_latent_heat_slope, fort_latent_heat_T0, fort_latent_heat_min
! deleted August 10, 2022:
!   fort_lame_coefficient 
!   fort_shear_modulus 
!   fort_linear_elastic_model 
! added August 27, 2022:
!   gravity_vector(SDIM)
! deleted August 27, 2022:
!   gravity,gravity_dir,invert_gravity
! added November 3, 2022: 
!   tension_init
! added December 3, 2022: 
!   gravity_boussinesq_vector
! deleted December 18, 2022:
!   fort_damping_coefficient
! added January 22, 2023:
!   fort_static_surface_tension_duration
! added March 3, 2023:
!   fort_static_tension
! deleted March 9, 2023:
!   VISCSTATIC
! deleted March 10, 2023:
!   fort_static_surface_tension_duration
! added April 19, 2023:
!   fort_grid_stretching_parameter
! deleted May 11, 2023: 
!   gravity_boussinesq_vector
! added June 8, deleted June 16 2023:
!   fort_growth_angle_primary_mat

      INTEGER_T, PARAMETER :: MOF_TRAINING_NDIM_DECISIONS=AMREX_SPACEDIM
      INTEGER_T, PARAMETER :: MOF_TRAINING_NDIM_CLASSIFY=AMREX_SPACEDIM-1

      Type training_model_type
        Type(Neural_Network) :: NN_ZHOUTENG_LOCAL
        Type(Decision_Tree) :: DT_ZHOUTENG_LOCAL
        Type(Random_Forest) :: RF_ZHOUTENG_LOCAL
      end Type training_model_type

      Type(training_model_type), allocatable, dimension(D_DECL(:,:,:),:) :: &
        training_array
      INTEGER_T :: training_finest_level=-1
      INTEGER_T :: training_lo(SDIM)
      INTEGER_T :: training_hi(SDIM)

      !https://www.stat.cmu.edu/~cshalizi/350-2006/lecture-10.pdf
      !S=sum_{c \in leaves(T)} n_{c} V_{c}
      !V_{c}=(1/n_{c})sum_{i \in C} (y_{i}-m_{c})^{2}
      Type branch_type
       INTEGER_T :: ndata
       INTEGER_T :: parent_id
       INTEGER_T :: parent_level
       INTEGER_T :: current_id
       INTEGER_T :: current_level
       INTEGER_T :: splittingrule
       INTEGER_T :: median_index
       REAL_T :: median_value
       REAL_T, pointer :: data_decisions(:,:) !datanum, data_idx
       REAL_T, pointer :: data_classify(:,:) !datanum, data_idx
       INTEGER_T :: child1_id
       INTEGER_T :: child_level
       INTEGER_T :: child2_id
      end Type branch_type

      Type level_branch_type
        INTEGER_T :: nbranches
        Type(branch_type), pointer :: branch_list(:)
      end Type level_branch_type

      Type tree_type
       INTEGER_T :: max_number_tree_levels
       INTEGER_T :: number_tree_levels
       INTEGER_T, pointer :: nbranches_level(:) 
       Type(level_branch_type), pointer :: branch_list_level(:)
      end Type tree_type

      Type(tree_type), allocatable, dimension(D_DECL(:,:,:),:) :: &
          decision_tree_array
      INTEGER_T :: decision_tree_finest_level=-1
      INTEGER_T :: decision_tree_lo(3)
      INTEGER_T :: decision_tree_hi(3)

      INTEGER_T, PARAMETER :: MAX_NUM_MATERIALS=10
      !num_interfaces=( (num_materials-1)*(num_materials-1)+num_materials-1 )/2
      INTEGER_T, PARAMETER :: MAX_NUM_INTERFACES=55
      INTEGER_T, PARAMETER :: MAX_NUM_SPECIES=10
       ! for user definable EOS, the EOS is defined based on "probtype"
       ! not "material_type"; but it is essential that
       ! 1<=material_type<=MAX_NUM_EOS
      INTEGER_T, PARAMETER :: MAX_NUM_EOS=998

      REAL_T, PARAMETER :: CLAMPED_EVERYWHERE_LS=99999.0D0
      REAL_T, PARAMETER :: CLAMPED_NO_WHERE_LS=-99999.0D0

#include "probdataf95.H"

      INTEGER_T, PARAMETER :: DEBUG_EVAPORATION=0
      INTEGER_T, PARAMETER :: EVAPORATION_iter_max=50
      REAL_T, PARAMETER :: EVAPORATION_TOL=1.0D-10

      INTEGER_T, PARAMETER :: OLD_DODECANE=1

      INTEGER_T, PARAMETER :: DEBUG_DYNAMIC_CONTACT_ANGLE=0

      REAL_T, PARAMETER :: GNBC_RADIUS=2.0d0

      INTEGER_T, PARAMETER :: ngrow_make_distance=3
      INTEGER_T, PARAMETER :: ngrow_distance=4

      REAL_T, PARAMETER :: GAMMA_SIMPLE_PARMS=1.4

         ! R=CP-CV
         ! CP=1.007D+7 Specific heat at constant pressure cgs ergs/(Kelvin g)
         ! GAMMA=CP/CV=1.39861
      REAL_T, PARAMETER :: R_AIR_PARMS=0.287D+7  ! ergs/(Kelvin g)
      REAL_T, PARAMETER :: CV_AIR_PARMS=0.72D+7  ! ergs/(Kelvin g)
      REAL_T, PARAMETER :: PCAV_TAIT=220.2726D0  ! cgs dyne/cm^2
      REAL_T, PARAMETER :: PCAV_TAIT_VACUUM=220.2726D0  ! dyne/cm^2
      REAL_T, PARAMETER :: A_TAIT=1.0D+6  ! dyne/cm^2
      REAL_T, PARAMETER :: B_TAIT=3.31D+9  ! dyne/cm^2
      REAL_T, PARAMETER :: RHOBAR_TAIT=1.0D0  ! g/cm^3
      REAL_T, PARAMETER :: GAMMA_TAIT=7.15D0

      REAL_T, PARAMETER :: P0_tillotson=1.0D+6 ! dyne/cm^2
      REAL_T, PARAMETER :: a_hydro_tillotson=0.7  ! dimensionless
      REAL_T, PARAMETER :: b_hydro_tillotson=0.15 ! dimensionless
      REAL_T, PARAMETER :: rho_IV_tillotson=0.958 ! g/cm^3
       ! 1 joule = 1 kg (m/s)^2=1 kg m^2/s^2
       ! 1 erg=1 g (cm/s)^2=1 g cm^2/s^2
       ! 
       ! MKS units of rho e: 1 joule/m^3=1 kg/(m s^2)
       ! CGS units of rho e: 1 erg/cm^3=1 g/(cm s^2)
       ! CGS units of e: 1 erg/cm^3  / (g/cm^3) = 1 erg/g
       ! units of pressure: 1 N/m^2=1 kg m/s^2/m^2=1kg/(m s^2)
       ! 1 bar=10^5 Pa = 10^5 N/m^2 = 10^6 dyne/cm^2
       ! 1 Pa=1 N/m^2 = 1 kg m/s^2 /m^2 =1kg/(m s^2)=1000 g/(m s^2)=
       !      10 g/(cm s^2)=10 dyne/cm^2
       ! 1 kbar=10^8 Pa=10^9 dyne/cm^2
      REAL_T, PARAMETER :: A_strain_tillotson=21.8D+9 ! dyne/cm^2
      REAL_T, PARAMETER :: B_strain_tillotson=132.5D+9 ! dyne/cm^2
      REAL_T, PARAMETER :: E0_tillotson=0.07D+12 ! erg/g
      REAL_T, PARAMETER :: E_IV_tillotson=0.00419D+12 ! erg/g
      REAL_T, PARAMETER :: E_CV_tillotson=0.025D+12   ! erg/g
      REAL_T, PARAMETER :: rho_cav_tillotson=0.995 ! g/cm^3
      REAL_T, PARAMETER :: P_cav_tillotson=5.0D+4 ! dyne/cm^2
      REAL_T, PARAMETER :: T_cav_tillotson=305.9 ! degrees Kelvin
      REAL_T, PARAMETER :: alpha_tillotson=10.0 
      REAL_T, PARAMETER :: beta_tillotson=5.0 

      REAL_T, PARAMETER :: TEMPERATURE_FLOOR=1.0D-20  ! default: 1.0D-20
      INTEGER_T, PARAMETER :: visual_RT_transform=1

      INTEGER_T, PARAMETER :: bubbleInPackedColumn=1001
      REAL_T, PARAMETER :: MASK_FINEST_TOL=1.0D-3

      REAL_T, PARAMETER :: OVERFLOW_CUTOFF=1.0D+20

      ! For inputs.curvature_converge with axis_dir=210 (sanity check),1.0D-12
      REAL_T, PARAMETER :: FACETOL_DVOL=1.0D-3

      REAL_T, PARAMETER :: VOFTOL_REDIST=1.0D-3

      REAL_T, PARAMETER :: FACETOL_REDIST=1.0D-2
      REAL_T, PARAMETER :: FACETOL_SANITY=1.0D-3
       ! Default: LS_CURV_TOL=1.0D-2 
       ! For inputs.curvature_converge with axis_dir=210 (sanity check),1.0D-12
      REAL_T, PARAMETER :: LS_CURV_TOL=1.0D-2
      REAL_T, PARAMETER :: LSTOL=1.0D-2
      REAL_T, PARAMETER :: VOFTOL_SLOPES=1.0D-2

      REAL_T, PARAMETER :: VOFTOL=1.0D-8

      REAL_T, PARAMETER :: VOFTOL_AREAFRAC=1.0D-1
       ! Default: VOFTOL_MULTI_VOLUME=1.0D-12
      REAL_T, PARAMETER :: VOFTOL_MULTI_VOLUME=1.0D-12
      REAL_T, PARAMETER :: VOFTOL_MULTI_VOLUME_SANITY=1.0D-8
       ! used for checking for 0 magnitude
      REAL_T, PARAMETER :: MLSVOFTOL=1.0D-14
       ! used for checking if centroid in box.
      REAL_T, PARAMETER :: CENTOL=1.0D-13
      REAL_T, PARAMETER :: EVAP_BISECTION_TOL=1.0D-12
       ! Default: INTERCEPT_TOL=1.0D-12
      REAL_T, PARAMETER :: INTERCEPT_TOL=1.0D-12
      REAL_T, PARAMETER :: FRAC_PAIR_TOL=1.0D-12

      REAL_T, PARAMETER :: TANGENT_EPS=1.0D-2

       !SEM_IMAGE_BC_ALG=1 => use "reflect_odd" for dirichlet,
       !                      use "reflect even" for Neumann
      INTEGER_T, PARAMETER :: SEM_IMAGE_BC_ALG=0

      INTEGER_T, PARAMETER :: POLYGON_LIST_MAX=1000

       ! variables for comparing with Fred Stern's experiment
      INTEGER_T, PARAMETER :: SHALLOW_M=100
      INTEGER_T, PARAMETER :: SHALLOW_N=1000
      REAL_T, PARAMETER :: SHALLOW_TIME=12.0

      INTEGER_T :: BLB_MATRIX=-32767
      INTEGER_T :: BLB_RHS=-32767
      INTEGER_T :: BLB_VEL=-32767
      INTEGER_T :: BLB_INT_MOM=-32767
      INTEGER_T :: BLB_ENERGY=-32767
      INTEGER_T :: BLB_MASS_VEL=-32767
      INTEGER_T :: BLB_VOL=-32767
      INTEGER_T :: BLB_CEN_INT=-32767
      INTEGER_T :: BLB_CEN_ACT=-32767
      INTEGER_T :: BLB_PERIM=-32767
      INTEGER_T :: BLB_PERIM_MAT=-32767
      INTEGER_T :: BLB_TRIPLE_PERIM=-32767
      INTEGER_T :: BLB_CELL_CNT=-32767
      INTEGER_T :: BLB_CELLVOL_CNT=-32767
      INTEGER_T :: BLB_MASS=-32767
      INTEGER_T :: BLB_PRES=-32767
      INTEGER_T :: num_elements_blobclass=-32767

      REAL_T shallow_water_data(0:SHALLOW_M,0:SHALLOW_N,2)
      REAL_T inflow_time(3000)
      REAL_T inflow_elevation(3000)
      REAL_T inflow_velocity(3000)
      REAL_T outflow_time(3000)
      REAL_T outflow_elevation(3000)
      REAL_T outflow_velocity(3000)
      INTEGER_T inflow_count,outflow_count
      INTEGER_T last_inflow_index,last_outflow_index

       ! variables from "rfiledata"
      REAL_T zstatic(0:300),rstatic(0:300) 

       ! variables from "pressure_bcs"
      real*8  dt_pressure_bcs
      real*8  time_pressure_bcs(0:100) ,  pressbc_pressure_bcs(0:100,1:3)
      INTEGER_T selectpress
       ! variables from "vel_bcs"
      real*8  timehist_velbc(0:100), &
        zpos_velbc(1:50),velbc_velbc(0:100,1:50), &
        period_velbc,rigidwall_velbc
      INTEGER_T itime_velbc,ipos_velbc

        ! level,index,dir
      REAL_T, allocatable, dimension(:,:,:) :: grid_cache
      INTEGER_T cache_index_low,cache_index_high,cache_max_level
      INTEGER_T :: grid_cache_allocated=0

       ! solve: div_X (1/w(x)) grad_X x=0  x(Xlo)=Xlo   x(Xhi)=Xhi
       ! index,dir
      REAL_T, allocatable, dimension(:,:) :: mapping_comp_to_phys
      REAL_T, allocatable, dimension(:,:) :: mapping_phys_to_comp
      INTEGER_T :: mapping_n_cell(0:2)
      INTEGER_T :: mapping_n_cell_max=0
      INTEGER_T :: mapping_allocated=0

      INTEGER_T :: number_of_source_regions=0
      INTEGER_T :: number_of_threads_regions=0

      type region_info_type
       INTEGER_T :: region_material_id
       REAL_T :: region_dt
       REAL_T :: region_mass_flux   ! e.g. kg/s
       REAL_T :: region_volume_flux ! e.g. m^3/s
       REAL_T :: region_temperature_prescribe ! e.g. degrees Kelvin
       REAL_T :: region_velocity_prescribe(SDIM) ! e.g. degrees Kelvin
       REAL_T :: region_energy_flux ! e.g. J/s = Watts
       REAL_T :: region_volume_raster ! e.g. m^3
       REAL_T :: region_volume        ! e.g. m^3
       REAL_T :: region_mass          ! e.g. kg
       REAL_T :: region_energy        ! e.g. J
       REAL_T :: region_energy_per_kelvin ! e.g. J/K
       REAL_T :: region_volume_after
       REAL_T :: region_mass_after
       REAL_T :: region_energy_after
      end type region_info_type

       ! first index: 1..number_of_source_regions
       ! second index: 0..number_of_threads_regions
      type(region_info_type), allocatable :: regions_list(:,:)
      
       ! interface particle container class
       ! refined data 
      type elem_contain_type
       INTEGER_T :: numElems
       INTEGER_T, pointer :: ElemData(:)
      end type elem_contain_type

       ! interface particle container class
       ! original coarse data
      type node_contain_type
       INTEGER_T :: numNodes
       INTEGER_T, pointer :: NodeData(:)
      end type node_contain_type

      type level_contain_type
        ! tid,partid,tilenum
       type(elem_contain_type), pointer :: level_elem_data(:,:,:)
       type(node_contain_type), pointer :: level_node_data(:,:,:)
        ! tid,tilenum,dir
       INTEGER_T, pointer :: tilelo3D(:,:,:)
       INTEGER_T, pointer :: tilehi3D(:,:,:)
       REAL_T, pointer :: xlo3D(:,:,:)
        ! tid,tilenum
       INTEGER_T, pointer :: gridno3D(:,:)
        ! tid
       INTEGER_T, pointer :: num_tiles_on_thread3D_proc(:)
       INTEGER_T :: max_num_tiles_on_thread3D_proc
       INTEGER_T :: num_grids_on_level
       INTEGER_T :: num_grids_on_level_proc
      end type level_contain_type

       ! level
      type(level_contain_type), dimension(:), allocatable :: contain_elem

      INTEGER_T :: container_allocated=0

       ! level
      INTEGER_T, allocatable, dimension(:) :: level_container_allocated

      type aux_contain_type
       INTEGER_T :: lo3D(3)
       INTEGER_T :: hi3D(3)
       REAL_T :: xlo3D(3)
       REAL_T :: xhi3D(3)
       REAL_T :: dx3D
       INTEGER_T :: aux_ncells_max_side
       REAL_T, dimension(:,:,:,:), pointer :: LS3D ! level set data
      end type aux_contain_type

      REAL_T, dimension(:,:,:,:), target, allocatable :: aux_xdata3D
      REAL_T, dimension(:,:,:,:), target, allocatable :: aux_FSIdata3D
      REAL_T, dimension(:,:,:,:), target, allocatable :: aux_masknbr3D

      INTEGER_T :: fort_num_local_aux_grids=0
      INTEGER_T :: aux_data_allocated=0
      type(aux_contain_type), dimension(:), allocatable :: contain_aux

      INTEGER_T :: used_probtypes(1000)
      INTEGER_T :: probtype_list_size

      ABSTRACT INTERFACE

      subroutine TEMPLATE_wallfunc( &
        dir, & ! =1,2,3
        data_dir, & ! =0,1,2
        dxmin, &
        x_projection_raster, &
        dx, &
        n_raster, & ! points to solid
        u, & !INTENT(in) uimage_raster_solid_frame(dir)
        uimage_tngt_mag, & !INTENT(in) 
        wall_model_velocity, & ! INTENT(in)
        dist_probe, & ! INTENT(in)
        dist_fluid, & ! INTENT(in)
        temperature_image, & !INTENT(in) 
        temperature_wall, & ! INTENT(in)      
        temperature_wall_max, & ! INTENT(in)      
        viscosity_molecular, & ! INTENT(in)      
        viscosity_eddy_wall, & ! INTENT(in)      
        y, & !INTENT(in) distance from image to wall
        ughost_tngt, & ! INTENT(out)
        im_fluid, &  ! INTENT(in)
        critical_length) ! INTENT(in) used for sanity check
      INTEGER_T, INTENT(in) :: dir ! 1,2,3
      INTEGER_T, INTENT(in) :: data_dir ! 0,1,2
      REAL_T, INTENT(in) :: dxmin
      REAL_T, INTENT(in), pointer :: x_projection_raster(:)
      REAL_T, INTENT(in), pointer :: dx(:)
      REAL_T, INTENT(in), pointer :: n_raster(:) ! points to solid
      INTEGER_T, INTENT(in) :: im_fluid
      REAL_T, INTENT(in) :: u !uimage_raster_solid_frame(dir)
      REAL_T, INTENT(in) :: uimage_tngt_mag
      REAL_T, INTENT(in) :: wall_model_velocity
      REAL_T, INTENT(in) :: dist_probe
      REAL_T, INTENT(in) :: dist_fluid
      REAL_T, INTENT(in) :: temperature_image
      REAL_T, INTENT(in) :: temperature_wall
      REAL_T, INTENT(in) :: temperature_wall_max
      REAL_T, INTENT(in) :: viscosity_molecular
      REAL_T, INTENT(in) :: viscosity_eddy_wall
      REAL_T, INTENT(in) :: y !delta_r
      REAL_T, INTENT(in) :: critical_length
      REAL_T, INTENT(out) :: ughost_tngt  ! dir direction
      end subroutine TEMPLATE_wallfunc

       ! returns (1/w) where w>>1 in "trouble" regions
      subroutine TEMPLATE_MAPPING_WEIGHT_COEFF(dir,wt,phys_x)
      INTEGER_T, INTENT(in) :: dir
      REAL_T, INTENT(out) :: wt
      REAL_T, INTENT(in) :: phys_x
      end subroutine TEMPLATE_MAPPING_WEIGHT_COEFF

      subroutine TEMPLATE_INIT_REGIONS_LIST(constant_density_all_time, &
          num_materials_in,num_threads_in)
      INTEGER_T, INTENT(in) :: num_materials_in
      INTEGER_T, INTENT(in) :: num_threads_in
      INTEGER_T, INTENT(in) :: constant_density_all_time(num_materials_in)
      end subroutine TEMPLATE_INIT_REGIONS_LIST

      subroutine TEMPLATE_CHARFN_REGION(region_id,x,cur_time,charfn_out)
      INTEGER_T, INTENT(in) :: region_id
      REAL_T, INTENT(in) :: x(SDIM)
      REAL_T, INTENT(in) :: cur_time
      REAL_T, INTENT(out) :: charfn_out
      end subroutine TEMPLATE_CHARFN_REGION

      subroutine TEMPLATE_T0_Boussinesq(x,dx,cur_time,T0)
      REAL_T, INTENT(in) :: x(SDIM)
      REAL_T, INTENT(in) :: dx(SDIM)
      REAL_T, INTENT(in) :: cur_time
      REAL_T, INTENT(out) :: T0
      end subroutine TEMPLATE_T0_Boussinesq

      subroutine TEMPLATE_V0_Coriolis(x,dx,cur_time,T0)
      REAL_T, INTENT(in) :: x(SDIM)
      REAL_T, INTENT(in) :: dx(SDIM)
      REAL_T, INTENT(in) :: cur_time
      REAL_T, INTENT(out) :: T0
      end subroutine TEMPLATE_V0_Coriolis

      subroutine TEMPLATE_THERMAL_K(x,dx,cur_time, &
        density, &
        temperature, &
        thermal_k, &
        im, &
        near_interface, &
        im_solid, &
        temperature_wall, &
        temperature_wall_max, &
        temperature_probe, &
        nrm) ! nrm points from solid to fluid
      INTEGER_T, INTENT(in) :: im
      INTEGER_T, INTENT(in) :: im_solid
      INTEGER_T, INTENT(in) :: near_interface
      REAL_T, INTENT(in) :: x(SDIM)
      REAL_T, INTENT(in) :: dx(SDIM)
      REAL_T, INTENT(in) :: cur_time
      REAL_T, INTENT(in) :: density
      REAL_T, INTENT(in) :: temperature
      REAL_T, INTENT(in) :: temperature_wall
      REAL_T, INTENT(in) :: temperature_wall_max
      REAL_T, INTENT(in) :: temperature_probe
      REAL_T, INTENT(in) :: nrm(SDIM) ! nrm points from solid to fluid
      REAL_T, INTENT(inout) :: thermal_k
      end subroutine TEMPLATE_THERMAL_K

      subroutine TEMPLATE_INTERFACE_TEMPERATURE( &
        interface_mass_transfer_model, &
        probe_constrain, &
        ireverse, &
        iten, &        
        xI, &        
        cur_time, &        
        prev_time, &        
        dt, &        
        TI, &
        YI, &
        user_override_TI_YI, &
        molar_mass, & ! index: 1..num_materials
        species_molar_mass, & ! index: 1..num_species_var
        ksrc_predict, &
        kdst_predict, &
        ksrc_physical, &
        kdst_physical, &
        T_probe_src, &
        T_probe_dst, &
        LL, &
        dxprobe_src, &
        dxprobe_dst, &
        num_materials_in, &
        num_species_var_in)
      INTEGER_T, INTENT(in) :: interface_mass_transfer_model
      INTEGER_T, INTENT(in) :: num_materials_in
      INTEGER_T, INTENT(in) :: num_species_var_in
      INTEGER_T, INTENT(in) :: probe_constrain
      INTEGER_T, INTENT(in) :: ireverse
      INTEGER_T, INTENT(in) :: iten
      REAL_T, INTENT(in) :: xI(SDIM)
      REAL_T, INTENT(in) :: cur_time
      REAL_T, INTENT(in) :: prev_time
      REAL_T, INTENT(in) :: dt
      REAL_T, INTENT(inout) :: TI
      REAL_T, INTENT(inout) :: YI
      INTEGER_T, INTENT(inout) :: user_override_TI_YI
      REAL_T, INTENT(in) :: molar_mass(num_materials_in)
      REAL_T, INTENT(in) :: species_molar_mass(num_species_var_in)
      REAL_T, INTENT(in) :: ksrc_predict
      REAL_T, INTENT(in) :: kdst_predict
      REAL_T, INTENT(in) :: ksrc_physical
      REAL_T, INTENT(in) :: kdst_physical
      REAL_T, INTENT(in) :: T_probe_src
      REAL_T, INTENT(in) :: T_probe_dst
      REAL_T, INTENT(in) :: LL
      REAL_T, INTENT(in) :: dxprobe_src
      REAL_T, INTENT(in) :: dxprobe_dst
      end subroutine TEMPLATE_INTERFACE_TEMPERATURE

      subroutine TEMPLATE_MDOT( &
        num_materials_in, &
        num_species_var_in, &
        interface_mass_transfer_model, &
        xI, & 
        ispec, &
        molar_mass, & ! 1..num_materials
        species_molar_mass, & ! 1..num_species_var+1
        im_source, &
        im_dest, &
        mdot, & ! INTENT(out)
        mdot_override, & ! INTENT(inout)
        ksrc_derived, &
        kdst_derived, &
        ksrc_physical, &
        kdst_physical, &
        T_probe_src, &
        T_probe_dst, &
        TI, &
        LL, &
        dxprobe_src, &
        dxprobe_dst)
      INTEGER_T, INTENT(in) :: interface_mass_transfer_model
      INTEGER_T, INTENT(in) :: num_materials_in
      INTEGER_T, INTENT(in) :: num_species_var_in
      INTEGER_T, INTENT(in) :: ispec
      INTEGER_T, INTENT(in) :: im_source
      INTEGER_T, INTENT(in) :: im_dest
      REAL_T, INTENT(in) :: xI(SDIM)
      REAL_T, INTENT(in) :: TI
      REAL_T, INTENT(in) :: molar_mass(num_materials_in)
      REAL_T, INTENT(in) :: species_molar_mass(num_species_var_in+1)
      REAL_T, INTENT(out) :: mdot
      INTEGER_T, INTENT(inout) :: mdot_override
      REAL_T, INTENT(in) :: ksrc_derived
      REAL_T, INTENT(in) :: kdst_derived
      REAL_T, INTENT(in) :: ksrc_physical
      REAL_T, INTENT(in) :: kdst_physical
      REAL_T, INTENT(in) :: T_probe_src
      REAL_T, INTENT(in) :: T_probe_dst
      REAL_T, INTENT(in) :: LL
      REAL_T, INTENT(in) :: dxprobe_src
      REAL_T, INTENT(in) :: dxprobe_dst
      end subroutine TEMPLATE_MDOT


      subroutine TEMPLATE_K_EFFECTIVE( &
        interface_mass_transfer_model, &
        ireverse, &
        iten, &        
        molar_mass, & ! index: 1..num_materials
        species_molar_mass, & ! index: 1..num_species_var
        k_model_predict, &
        k_model_correct, &
        k_physical_base, &
        T_probe_src, &
        T_probe_dst, &
        dxprobe_src, &
        dxprobe_dst, &
        LL, &
        num_materials_in, &
        num_species_var_in)
      INTEGER_T, INTENT(in) :: interface_mass_transfer_model
      INTEGER_T, INTENT(in) :: num_materials_in
      INTEGER_T, INTENT(in) :: num_species_var_in
      INTEGER_T, INTENT(in) :: ireverse
      INTEGER_T, INTENT(in) :: iten
      REAL_T, INTENT(in) :: molar_mass(num_materials_in)
      REAL_T, INTENT(in) :: species_molar_mass(num_species_var_in)
      REAL_T, INTENT(in) :: k_model_predict(2) ! src,dst
      REAL_T, INTENT(inout) :: k_model_correct(2) ! src,dst
      REAL_T, INTENT(in) :: k_physical_base(2) ! src, dst
      REAL_T, INTENT(in) :: T_probe_src
      REAL_T, INTENT(in) :: T_probe_dst
      REAL_T, INTENT(in) :: LL
      REAL_T, INTENT(in) :: dxprobe_src
      REAL_T, INTENT(in) :: dxprobe_dst
      end subroutine TEMPLATE_K_EFFECTIVE


      subroutine TEMPLATE_reference_wavelen(wavelen)
      REAL_T, INTENT(inout) :: wavelen
      end subroutine TEMPLATE_reference_wavelen

      subroutine TEMPLATE_DELETE_REGIONS_LIST()
      end subroutine TEMPLATE_DELETE_REGIONS_LIST

      subroutine TEMPLATE_INIT_MODULE()
      end subroutine TEMPLATE_INIT_MODULE

      subroutine TEMPLATE_DEALLOCATE_MODULE()
      end subroutine TEMPLATE_DEALLOCATE_MODULE

      subroutine TEMPLATE_correct_pres_rho_hydrostatic( &
        i,j,k,level, &
        angular_velocity, &
        centrifugal_force_factor, &
        dt, &
        rho_hydrostatic, &
        pres_hydrostatic, &
        state_ptr)
      INTEGER_T, INTENT(in) :: i,j,k,level
      REAL_T, INTENT(in) :: angular_velocity
      REAL_T, INTENT(in) :: centrifugal_force_factor
      REAL_T, INTENT(in) :: dt
      REAL_T, INTENT(inout) :: rho_hydrostatic
      REAL_T, INTENT(inout) :: pres_hydrostatic
      REAL_T, INTENT(in),pointer :: state_ptr(D_DECL(:,:,:),:)
      end subroutine TEMPLATE_correct_pres_rho_hydrostatic
        
      subroutine TEMPLATE_CFL_HELPER(time,dir,uu,dx)
      INTEGER_T, INTENT(in) :: dir
      REAL_T, INTENT(in) :: time
      REAL_T, INTENT(inout) :: uu
      REAL_T, INTENT(in) :: dx(SDIM)
      end subroutine TEMPLATE_CFL_HELPER


      subroutine TEMPLATE_SUMINT(GRID_DATA_IN,increment_out1, &
            increment_out2,nsum1,nsum2,isweep)
      use probcommon_module_types

      INTEGER_T, INTENT(in) :: nsum1,nsum2,isweep
      type(user_defined_sum_int_type), INTENT(in) :: GRID_DATA_IN
      REAL_T, INTENT(inout) :: increment_out1(nsum1)
      REAL_T, INTENT(inout) :: increment_out2(nsum2)
      end subroutine TEMPLATE_SUMINT
      
      subroutine TEMPLATE_LS(x,t,LS,nmat)
      INTEGER_T, INTENT(in) :: nmat
      REAL_T, INTENT(in) :: x(SDIM)
      REAL_T, INTENT(in) :: t
      REAL_T, INTENT(out) :: LS(nmat)
      end subroutine TEMPLATE_LS

      subroutine TEMPLATE_OVERRIDE_TAGFLAG(xsten,nhalf,time,rflag,tagflag)
      INTEGER_T, INTENT(in) :: nhalf
      REAL_T, INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      REAL_T, INTENT(in) :: time
      REAL_T, INTENT(inout) :: rflag
      INTEGER_T, INTENT(inout) :: tagflag
      end subroutine TEMPLATE_OVERRIDE_TAGFLAG

      subroutine TEMPLATE_AUX_DATA(auxcomp,x,LS)
      INTEGER_T, INTENT(in) :: auxcomp
      REAL_T, INTENT(in) :: x(3)
      REAL_T, INTENT(out) :: LS
      end subroutine TEMPLATE_AUX_DATA

      subroutine TEMPLATE_OVERRIDE_FSI_SIGN_LS_VEL_TEMP( &
        exterior_BB, &
        interior_BB, &
        xcell,time,LS,VEL,TEMP,MASK,lev77,im_part,part_id)
      REAL_T, INTENT(in) :: exterior_BB(3,2)
      REAL_T, INTENT(in) :: interior_BB(3,2)
      REAL_T, INTENT(in) :: xcell(3)
      REAL_T, INTENT(in) :: time
      REAL_T, INTENT(out) :: LS
      REAL_T, INTENT(out) :: VEL(3)
      REAL_T, INTENT(out) :: TEMP
      INTEGER_T, INTENT(out) :: MASK
      INTEGER_T, INTENT(in) :: lev77 !lev77=-1 for aux, >=0 otherwise.
      INTEGER_T, INTENT(in) :: im_part
      INTEGER_T, INTENT(in) :: part_id
      end subroutine TEMPLATE_OVERRIDE_FSI_SIGN_LS_VEL_TEMP

      subroutine TEMPLATE_GET_OUTSIDE_POINT( &
        exterior_BB, &
        xcell,time,x_outside,im_part,part_id)
      REAL_T, INTENT(in) :: exterior_BB(3,2)
      REAL_T, INTENT(in) :: xcell(3)
      REAL_T, INTENT(in) :: time
      REAL_T, INTENT(out) :: x_outside(3)
      INTEGER_T, INTENT(in) :: im_part
      INTEGER_T, INTENT(in) :: part_id
      end subroutine TEMPLATE_GET_OUTSIDE_POINT

      subroutine TEMPLATE_BOUNDING_BOX_AUX(auxcomp, &
          minnode,maxnode,LS_FROM_SUBROUTINE,aux_ncells_max_side)
      INTEGER_T, INTENT(in) :: auxcomp
      REAL_T, INTENT(inout) :: minnode(3)
      REAL_T, INTENT(inout) :: maxnode(3)
      INTEGER_T, INTENT(out) :: LS_FROM_SUBROUTINE
      INTEGER_T, INTENT(out) :: aux_ncells_max_side
      end subroutine TEMPLATE_BOUNDING_BOX_AUX


      subroutine TEMPLATE_check_vel_rigid(x,t,vel,dir)
       REAL_T, INTENT(in) :: x(SDIM)
       REAL_T, INTENT(in) :: t
       REAL_T, INTENT(in) :: vel
       INTEGER_T, INTENT(in) :: dir
      end subroutine TEMPLATE_check_vel_rigid

      subroutine TEMPLATE_clamped_LS(x,t,LS,vel,temperature,prescribed_flag,dx)
       REAL_T, INTENT(in) :: x(SDIM)
       REAL_T, INTENT(in) :: dx(SDIM)
       REAL_T, INTENT(in) :: t
       REAL_T, INTENT(out) :: LS
       REAL_T, INTENT(out) :: vel(SDIM)
       REAL_T, INTENT(out) :: temperature
       INTEGER_T, INTENT(out) :: prescribed_flag
      end subroutine TEMPLATE_clamped_LS

      subroutine TEMPLATE_VEL(x,t,LS,VEL,velsolid_flag,dx,nmat)
      INTEGER_T, INTENT(in) :: nmat
      REAL_T, INTENT(in) :: x(SDIM)
      REAL_T, INTENT(in) :: t
      REAL_T, INTENT(in) :: dx(SDIM)
      REAL_T, INTENT(in) :: LS(nmat)
      REAL_T, INTENT(out) :: VEL(SDIM)
      INTEGER_T, INTENT(in) :: velsolid_flag
      end subroutine TEMPLATE_VEL

      subroutine TEMPLATE_EOS(rho,massfrac_var, &
        internal_energy,pressure, &
        imattype,im,num_species_var_in)
      INTEGER_T, INTENT(in) :: imattype,im,num_species_var_in
      REAL_T, INTENT(in) :: rho
      REAL_T, INTENT(in) :: massfrac_var(num_species_var_in+1)
      REAL_T, INTENT(in) :: internal_energy
      REAL_T, INTENT(out) :: pressure
      end subroutine TEMPLATE_EOS

      subroutine TEMPLATE_UNITLESS_EXPANSION_FACTOR( &
        im,temperature,temperature_base,expansion_factor)
      INTEGER_T, INTENT(in) :: im
      REAL_T, INTENT(in) :: temperature
      REAL_T, INTENT(in) :: temperature_base
      REAL_T, INTENT(out) :: expansion_factor
      end subroutine TEMPLATE_UNITLESS_EXPANSION_FACTOR

      subroutine TEMPLATE_INTERNAL_GRAVITY_WAVE_FLAG( &
        internal_wave_exists)
      INTEGER_T, INTENT(out) :: internal_wave_exists
      end subroutine TEMPLATE_INTERNAL_GRAVITY_WAVE_FLAG

      subroutine TEMPLATE_dVdT(dVdT,massfrac_var, &
        pressure,temperature, &
        imattype,im,num_species_var_in)
      INTEGER_T, INTENT(in) :: imattype,im,num_species_var_in
      REAL_T, INTENT(in) :: pressure,temperature
      REAL_T, INTENT(in) :: massfrac_var(num_species_var_in+1)
      REAL_T, INTENT(out) :: dVdT
      end subroutine TEMPLATE_dVdT


      subroutine TEMPLATE_SOUNDSQR(rho,massfrac_var, &
        internal_energy,soundsqr, &
        imattype,im,num_species_var_in)
      INTEGER_T, INTENT(in) :: imattype,im,num_species_var_in
      REAL_T, INTENT(in) :: rho
      REAL_T, INTENT(in) :: massfrac_var(num_species_var_in+1)
      REAL_T, INTENT(in) :: internal_energy
      REAL_T, INTENT(out) :: soundsqr
      end subroutine TEMPLATE_SOUNDSQR

      subroutine TEMPLATE_INTERNAL(rho,massfrac_var, &
        temperature,local_internal_energy, &
        imattype,im,num_species_var_in)
      INTEGER_T, INTENT(in) :: imattype,im,num_species_var_in
      REAL_T, INTENT(in) :: rho
      REAL_T, INTENT(in) :: massfrac_var(num_species_var_in+1)
      REAL_T, INTENT(in) :: temperature 
      REAL_T, INTENT(out) :: local_internal_energy
      end subroutine TEMPLATE_INTERNAL

      subroutine TEMPLATE_TEMPERATURE(rho,massfrac_var, &
        temperature,internal_energy, &
        imattype,im,num_species_var_in)
      INTEGER_T, INTENT(in) :: imattype,im,num_species_var_in
      REAL_T, INTENT(in) :: rho
      REAL_T, INTENT(in) :: massfrac_var(num_species_var_in+1)
      REAL_T, INTENT(out) :: temperature 
      REAL_T, INTENT(in) :: internal_energy
      end subroutine TEMPLATE_TEMPERATURE

      subroutine TEMPLATE_PRES(x,t,LS,PRES,nmat)
      INTEGER_T, INTENT(in) :: nmat
      REAL_T, INTENT(in) :: x(SDIM)
      REAL_T, INTENT(in) :: t
      REAL_T, INTENT(in) :: LS(nmat)
      REAL_T, INTENT(out) :: PRES
      end subroutine TEMPLATE_PRES

      subroutine TEMPLATE_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
      INTEGER_T, INTENT(in) :: bcflag
      INTEGER_T, INTENT(in) :: nmat
      INTEGER_T, INTENT(in) :: nstate_mat
      REAL_T, INTENT(in) :: x(SDIM)
      REAL_T, INTENT(in) :: t
      REAL_T, INTENT(in) :: LS(nmat)
      REAL_T, INTENT(out) :: STATE(nmat*nstate_mat)
      end subroutine TEMPLATE_STATE

      subroutine TEMPLATE_LS_BC(xwall,xghost,t,LS, &
       LS_in,dir,side,dx,nmat)
      INTEGER_T, INTENT(in) :: nmat
      REAL_T, INTENT(in) :: xwall
      REAL_T, INTENT(in) :: xghost(SDIM)
      REAL_T, INTENT(in) :: t
      REAL_T, INTENT(inout) :: LS(nmat)
      REAL_T, INTENT(in) :: LS_in(nmat)
      INTEGER_T, INTENT(in) :: dir,side
      REAL_T, INTENT(in) :: dx(SDIM)
      end subroutine TEMPLATE_LS_BC

      subroutine TEMPLATE_VEL_BC(xwall,xghost,t,LS, &
        VEL,VEL_in,veldir,dir,side,dx,nmat)
      INTEGER_T, INTENT(in) :: nmat
      REAL_T, INTENT(in) :: xwall
      REAL_T, INTENT(in) :: xghost(SDIM)
      REAL_T, INTENT(in) :: t
      REAL_T, INTENT(in) :: LS(nmat)
      REAL_T, INTENT(inout) :: VEL
      REAL_T, INTENT(in) :: VEL_in
      INTEGER_T, INTENT(in) :: veldir,dir,side
      REAL_T, INTENT(in) :: dx(SDIM)
      end subroutine TEMPLATE_VEL_BC

      subroutine TEMPLATE_PRES_BC(xwall,xghost,t,LS, &
        PRES,PRES_in,dir,side,dx,nmat)
      INTEGER_T, INTENT(in) :: nmat
      REAL_T, INTENT(in) :: xwall
      REAL_T, INTENT(in) :: xghost(SDIM)
      REAL_T, INTENT(in) :: t
      REAL_T, INTENT(in) :: LS(nmat)
      REAL_T, INTENT(inout) :: PRES
      REAL_T, INTENT(in) :: PRES_in
      INTEGER_T, INTENT(in) :: dir,side
      REAL_T, INTENT(in) :: dx(SDIM)
      end subroutine TEMPLATE_PRES_BC

      subroutine TEMPLATE_STATE_BC(xwall,xghost,t,LS, &
       STATE,STATE_merge,STATE_in,im,istate,dir,side,dx, &
       nmat)
      INTEGER_T, INTENT(in) :: nmat
      REAL_T, INTENT(in) :: xwall
      REAL_T, INTENT(in) :: xghost(SDIM)
      REAL_T, INTENT(in) :: t
      REAL_T, INTENT(in) :: LS(nmat)
      REAL_T, INTENT(inout) :: STATE
      REAL_T, INTENT(inout) :: STATE_merge
      REAL_T, INTENT(in) :: STATE_in
      INTEGER_T, INTENT(in) :: dir,side
      REAL_T, INTENT(in) :: dx(SDIM)
      INTEGER_T, INTENT(in) :: istate,im
      end subroutine TEMPLATE_STATE_BC

      subroutine TEMPLATE_HEATSOURCE( &
        im,VFRAC, &
        time, &
        x, &
        xsten, & ! xsten(-nhalf:nhalf,SDIM)
        nhalf, &
        temp, &
        heat_source,den,CV,dt, &
        nmat)
      INTEGER_T, INTENT(in) :: nmat
      INTEGER_T, INTENT(in) :: im
      REAL_T, INTENT(in) :: VFRAC(nmat)
      REAL_T, INTENT(in) :: time
      INTEGER_T, INTENT(in) :: nhalf
      REAL_T, INTENT(in) :: x(SDIM)
      REAL_T, INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      REAL_T, INTENT(in) :: temp(nmat)
      REAL_T, INTENT(in) :: den(nmat)
      REAL_T, INTENT(in) :: CV(nmat)
      REAL_T, INTENT(in) :: dt
      REAL_T, INTENT(out) :: heat_source
      end subroutine TEMPLATE_HEATSOURCE

      subroutine TEMPLATE_EB_heat_source(time,dt,xsten,nhalf, &
        heat_flux,heat_dir,heat_side)
      INTEGER_T, INTENT(in) :: nhalf
      REAL_T, dimension(-nhalf:nhalf,SDIM), INTENT(in) :: xsten
      REAL_T, INTENT(in) :: time
      REAL_T, INTENT(in) :: dt
      REAL_T, INTENT(out) :: heat_flux
      INTEGER_T, INTENT(out) :: heat_dir
      INTEGER_T, INTENT(out) :: heat_side
      end subroutine TEMPLATE_EB_heat_source

      subroutine TEMPLATE_velfreestream(problen,local_buffer)
      REAL_T, INTENT(inout) :: local_buffer(2*SDIM)
      REAL_T, INTENT(in)    :: problen(SDIM)
      end subroutine TEMPLATE_velfreestream

      subroutine TEMPLATE_nucleation(nucleate_in,xsten,nhalf,make_seed)
      use probcommon_module_types
      INTEGER_T, INTENT(in) :: nhalf
      REAL_T, dimension(-nhalf:nhalf,SDIM), INTENT(in) :: xsten
      INTEGER_T, INTENT(inout) :: make_seed
      type(nucleation_parm_type_input), INTENT(in) :: nucleate_in
      end subroutine TEMPLATE_nucleation

      subroutine TEMPLATE_ICE_SUBSTRATE_DISTANCE(xtarget,dist)
      REAL_T, INTENT(in) :: xtarget(SDIM)
      REAL_T, INTENT(out) :: dist
      end subroutine TEMPLATE_ICE_SUBSTRATE_DISTANCE


      subroutine TEMPLATE_microcell_heat_coeff(heatcoeff,dx,veldir)
      REAL_T, INTENT(in) :: dx(SDIM)
      INTEGER_T, INTENT(in) :: veldir
      REAL_T, INTENT(inout) :: heatcoeff
      end subroutine TEMPLATE_microcell_heat_coeff

      subroutine TEMPLATE_ASSIMILATE(assimilate_in,assimilate_out, &
         i,j,k,cell_flag)
      use probcommon_module_types
      type(assimilate_parm_type), INTENT(in) :: assimilate_in
      type(assimilate_out_parm_type), INTENT(inout) :: assimilate_out
      INTEGER_T, INTENT(in) :: i,j,k,cell_flag
      end subroutine TEMPLATE_ASSIMILATE

      subroutine TEMPLATE_FSI_SLICE(xmap3D,xslice3D,problo3D,probhi3D,dx_slice)
      REAL_T, INTENT(in) :: dx_slice
      INTEGER_T, INTENT(inout) :: xmap3D(3)
      REAL_T, INTENT(inout) :: xslice3D(3)
      REAL_T, INTENT(out) :: problo3D(3)
      REAL_T, INTENT(out) :: probhi3D(3)
      end subroutine TEMPLATE_FSI_SLICE

      subroutine TEMPLATE_OPEN_CASFILE(part_id,unit_id,file_format)
      INTEGER_T, INTENT(in) :: part_id
      INTEGER_T, INTENT(in) :: unit_id
      INTEGER_T, INTENT(out) :: file_format
      end subroutine TEMPLATE_OPEN_CASFILE

      subroutine TEMPLATE_OPEN_AUXFILE(part_id,unit_id,file_format)
      INTEGER_T, INTENT(in) :: part_id
      INTEGER_T, INTENT(in) :: unit_id
      INTEGER_T, INTENT(out) :: file_format
      end subroutine TEMPLATE_OPEN_AUXFILE

      subroutine TEMPLATE_ORDER_NODES(nodes,nodemap)
      REAL_T, INTENT(in) :: nodes(3,3) ! dir,nodenum
      INTEGER_T, INTENT(inout) :: nodemap(3)
      end subroutine TEMPLATE_ORDER_NODES

      subroutine TEMPLATE_VARIABLE_SURFACE_TENSION( &
        xpos, &
        time, &
        iten, &
        temperature, &
        tension)
      INTEGER_T, INTENT(in) :: iten
      REAL_T, INTENT(in) :: time,temperature
      REAL_T, INTENT(in) :: xpos(SDIM)
      REAL_T, INTENT(inout) :: tension
      end subroutine TEMPLATE_VARIABLE_SURFACE_TENSION

      subroutine TEMPLATE_VARIABLE_LATENT_HEAT( &
        iten, &
        temperature, &
        latent_heat)
      INTEGER_T, INTENT(in) :: iten
      REAL_T, INTENT(in) :: temperature
      REAL_T, INTENT(inout) :: latent_heat ! always positive
      end subroutine TEMPLATE_VARIABLE_LATENT_HEAT


      END INTERFACE

      PROCEDURE(TEMPLATE_INIT_MODULE), POINTER :: SUB_INIT_MODULE
      PROCEDURE(TEMPLATE_DEALLOCATE_MODULE), POINTER :: SUB_DEALLOCATE_MODULE
      PROCEDURE(TEMPLATE_correct_pres_rho_hydrostatic), POINTER :: &
              SUB_correct_pres_rho_hydrostatic
      PROCEDURE(TEMPLATE_CFL_HELPER), POINTER :: SUB_CFL_HELPER
      PROCEDURE(TEMPLATE_SUMINT), POINTER :: SUB_SUMINT
      PROCEDURE(TEMPLATE_LS), POINTER :: SUB_LS
      PROCEDURE(TEMPLATE_OVERRIDE_TAGFLAG), POINTER :: SUB_OVERRIDE_TAGFLAG
      PROCEDURE(TEMPLATE_AUX_DATA), POINTER :: SUB_AUX_DATA
      PROCEDURE(TEMPLATE_OVERRIDE_FSI_SIGN_LS_VEL_TEMP), POINTER :: &
        SUB_OVERRIDE_FSI_SIGN_LS_VEL_TEMP
      PROCEDURE(TEMPLATE_GET_OUTSIDE_POINT), POINTER :: &
        SUB_GET_OUTSIDE_POINT
      PROCEDURE(TEMPLATE_BOUNDING_BOX_AUX), POINTER :: SUB_BOUNDING_BOX_AUX
      PROCEDURE(TEMPLATE_check_vel_rigid), POINTER :: SUB_check_vel_rigid
      PROCEDURE(TEMPLATE_clamped_LS), POINTER :: SUB_clamped_LS_no_scale
      PROCEDURE(TEMPLATE_VEL), POINTER :: SUB_VEL
      PROCEDURE(TEMPLATE_EOS), POINTER :: SUB_EOS
      PROCEDURE(TEMPLATE_UNITLESS_EXPANSION_FACTOR), POINTER :: &
              SUB_UNITLESS_EXPANSION_FACTOR
      PROCEDURE(TEMPLATE_INTERNAL_GRAVITY_WAVE_FLAG), POINTER :: &
              SUB_INTERNAL_GRAVITY_WAVE_FLAG
      PROCEDURE(TEMPLATE_dVdT), POINTER :: SUB_dVdT
      PROCEDURE(TEMPLATE_SOUNDSQR), POINTER :: SUB_SOUNDSQR
      PROCEDURE(TEMPLATE_INTERNAL), POINTER :: SUB_INTERNAL
      PROCEDURE(TEMPLATE_TEMPERATURE), POINTER :: SUB_TEMPERATURE
      PROCEDURE(TEMPLATE_PRES), POINTER :: SUB_PRES
      PROCEDURE(TEMPLATE_STATE), POINTER :: SUB_STATE
      PROCEDURE(TEMPLATE_LS_BC), POINTER :: SUB_LS_BC
      PROCEDURE(TEMPLATE_VEL_BC), POINTER :: SUB_VEL_BC
      PROCEDURE(TEMPLATE_PRES_BC), POINTER :: SUB_PRES_BC
      PROCEDURE(TEMPLATE_STATE_BC), POINTER :: SUB_STATE_BC
      PROCEDURE(TEMPLATE_HEATSOURCE), POINTER :: SUB_HEATSOURCE
      PROCEDURE(TEMPLATE_EB_heat_source), POINTER :: SUB_EB_heat_source
      PROCEDURE(TEMPLATE_velfreestream), POINTER :: SUB_velfreestream
      PROCEDURE(TEMPLATE_nucleation), POINTER :: SUB_nucleation
      PROCEDURE(TEMPLATE_ICE_SUBSTRATE_DISTANCE), POINTER :: &
              SUB_ICE_SUBSTRATE_DISTANCE
      PROCEDURE(TEMPLATE_microcell_heat_coeff), POINTER :: &
              SUB_microcell_heat_coeff
      PROCEDURE(TEMPLATE_ASSIMILATE), POINTER :: SUB_ASSIMILATE

      PROCEDURE(TEMPLATE_FSI_SLICE), POINTER :: SUB_FSI_SLICE
      PROCEDURE(TEMPLATE_OPEN_CASFILE), POINTER :: SUB_OPEN_CASFILE
      PROCEDURE(TEMPLATE_OPEN_AUXFILE), POINTER :: SUB_OPEN_AUXFILE
      PROCEDURE(TEMPLATE_ORDER_NODES), POINTER :: SUB_ORDER_NODES

      PROCEDURE(TEMPLATE_wallfunc), POINTER :: SUB_wallfunc

      PROCEDURE(TEMPLATE_MAPPING_WEIGHT_COEFF), POINTER :: &
              SUB_MAPPING_WEIGHT_COEFF

      PROCEDURE(TEMPLATE_INIT_REGIONS_LIST), POINTER :: SUB_INIT_REGIONS_LIST
      PROCEDURE(TEMPLATE_CHARFN_REGION), POINTER :: SUB_CHARFN_REGION
      PROCEDURE(TEMPLATE_DELETE_REGIONS_LIST), POINTER ::  &
              SUB_DELETE_REGIONS_LIST

      PROCEDURE(TEMPLATE_T0_Boussinesq), POINTER :: SUB_T0_Boussinesq
      PROCEDURE(TEMPLATE_V0_Coriolis), POINTER :: SUB_V0_Coriolis

      PROCEDURE(TEMPLATE_THERMAL_K), POINTER :: SUB_THERMAL_K
      PROCEDURE(TEMPLATE_INTERFACE_TEMPERATURE), POINTER :: &
              SUB_INTERFACE_TEMPERATURE
      PROCEDURE(TEMPLATE_MDOT), POINTER :: &
              SUB_MDOT
      PROCEDURE(TEMPLATE_K_EFFECTIVE), POINTER :: &
              SUB_K_EFFECTIVE

      PROCEDURE(TEMPLATE_reference_wavelen), POINTER :: SUB_reference_wavelen

      PROCEDURE(TEMPLATE_VARIABLE_SURFACE_TENSION), POINTER :: &
              SUB_VARIABLE_SURFACE_TENSION
      PROCEDURE(TEMPLATE_VARIABLE_LATENT_HEAT), POINTER :: &
              SUB_VARIABLE_LATENT_HEAT

contains

      subroutine EOS_tait_ADIABATIC_rhohydro(rho,pressure)
      IMPLICIT NONE

      REAL_T, INTENT(in) :: rho
      REAL_T, INTENT(out) :: pressure
      REAL_T A,B,rhobar,pcav


      A=A_TAIT   ! dyne/cm^2
      B=B_TAIT  ! dyne/cm^2
      rhobar=fort_denconst(1) ! g/cm^3

      if (rhobar.lt.0.001) then
       print *,"rhobar invalid in eos tait adiabatic rhohydro"
       stop
      endif

      pcav=PCAV_TAIT

      if (rho.gt.zero) then
       ! do nothing
      else
       print *,"rho invalid"
       stop
      endif

      pressure=B*( (rho/rhobar)**GAMMA_TAIT - one ) + A

      if (pressure.lt.pcav) then
       pressure=pcav
      endif

      return
      end subroutine EOS_tait_ADIABATIC_rhohydro


      subroutine EOS_tait_ADIABATIC(rho,pressure)
      IMPLICIT NONE

      REAL_T rho,pressure
      REAL_T A,B,rhobar,GAMMA,pcav

      A=A_TAIT   ! dyne/cm^2
      B=B_TAIT  ! dyne/cm^2
      rhobar=RHOBAR_TAIT ! g/cm^3
      GAMMA=GAMMA_TAIT
      pcav=PCAV_TAIT 

      if (rho.le.zero) then
       print *,"rho invalid"
       stop
      endif

      pressure=B*( (rho/rhobar)**GAMMA - one ) + A

      if (pressure.lt.pcav) then
       pressure=pcav
      endif

      return
      end subroutine EOS_tait_ADIABATIC

end module probcommon_module

