#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

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
      use amrex_fort_module, only : amrex_real

      use NeuralNetwork
      use DecisionTree
      use RandomForest

 
      type law_of_wall_parm_type
      integer :: level
      integer :: finest_level
      integer :: bfact
      real(amrex_real) :: visc_coef
      real(amrex_real) :: time
      real(amrex_real) :: dt
      real(amrex_real), pointer :: usolid_raster(:)
      real(amrex_real), pointer :: n_raster(:) ! points to solid
      real(amrex_real), pointer :: x_image_raster(:)
      real(amrex_real), pointer :: x_probe_raster(:)
      real(amrex_real), pointer :: x_projection_raster(:)
      real(amrex_real), pointer :: dx(:)
      real(amrex_real) :: dxmin
      real(amrex_real), pointer :: xlo(:)
      integer, pointer :: fablo(:)
      integer, pointer :: fabhi(:)
      real(amrex_real), pointer, dimension(D_DECL(:,:,:),:) :: LSCP
      real(amrex_real), pointer, dimension(D_DECL(:,:,:),:) :: LSFD
      real(amrex_real), pointer, dimension(D_DECL(:,:,:),:) :: state ! nden comp.
      real(amrex_real), pointer, dimension(D_DECL(:,:,:),:) :: ufluid
      real(amrex_real), pointer, dimension(D_DECL(:,:,:),:) :: usolid
      end type law_of_wall_parm_type

      type user_defined_sum_int_type
       integer ncomp_sum_int_user1
       integer ncomp_sum_int_user2
       integer ncomp_sum_int_user12
       real(amrex_real), pointer :: problo(:)      
       real(amrex_real), pointer :: probhi(:) 
       integer :: igrid,jgrid,kgrid
       real(amrex_real) :: volgrid
       integer :: nhalf
       integer :: bfact
       integer :: den_ncomp
       integer :: level
       integer :: finest_level
       integer, pointer :: tilelo(:)
       integer, pointer :: tilehi(:)
       integer, pointer :: fablo(:)
       integer, pointer :: fabhi(:)
       real(amrex_real), pointer :: xlo(:)
       real(amrex_real), pointer :: dx(:)
       real(amrex_real), pointer :: xsten(:,:)
       real(amrex_real) :: time
       real(amrex_real), pointer, dimension(D_DECL(:,:,:),:) :: cellten
       real(amrex_real), pointer, dimension(D_DECL(:,:,:),:) :: lsfab
        ! 1..num_materials*ngeom_recon
       real(amrex_real), pointer, dimension(D_DECL(:,:,:),:) :: slopes
        ! num_materials * num_state_material
        ! density1,temperature1,species1_1,...,species_N_1
        ! density2,temperature2,species1_2,...,species_N_2
        ! density3,temperature3,species1_3,...,species_N_3
       real(amrex_real), pointer, dimension(D_DECL(:,:,:),:) :: den
       real(amrex_real), pointer, dimension(D_DECL(:,:,:),:) :: vel
       real(amrex_real), pointer, dimension(D_DECL(:,:,:),:) :: visco
      end type user_defined_sum_int_type

      type nucleation_parm_type_input
       integer :: tid
       integer :: local_freezing_model
       real(amrex_real) :: LL
       integer :: i,j,k
       integer :: im_source
       integer :: im_dest
       real(amrex_real) :: dxmaxLS
       integer :: bfact
       integer :: level
       integer :: finest_level
       real(amrex_real), pointer :: dx(:)
       real(amrex_real), pointer :: xlo(:)
       integer :: nstate
       integer, pointer :: fablo(:)
       integer, pointer :: fabhi(:)
       real(amrex_real), pointer, dimension(D_DECL(:,:,:),:) :: EOS
       real(amrex_real), pointer, dimension(D_DECL(:,:,:),:) :: Snew
       real(amrex_real), pointer, dimension(D_DECL(:,:,:),:) :: LSnew
       real(amrex_real), pointer, dimension(D_DECL(:,:,:)) :: pres
       real(amrex_real), pointer, dimension(D_DECL(:,:,:)) :: pres_eos
       integer :: custom_nucleation_model
       integer :: do_the_nucleate
       integer :: nucleate_pos_size
       real(amrex_real), pointer :: nucleate_pos(:)
       real(amrex_real), pointer :: nucleation_temp(:)
       real(amrex_real), pointer :: nucleation_pressure(:)
       real(amrex_real), pointer :: nucleation_pmg(:)
       real(amrex_real), pointer :: nucleation_mach(:)
       real(amrex_real), pointer :: cavitation_pressure(:)
       real(amrex_real), pointer :: cavitation_vapor_density(:)
       real(amrex_real), pointer :: cavitation_tension(:)
       real(amrex_real) :: local_TSAT
       real(amrex_real) :: prev_time
       real(amrex_real) :: cur_time
       real(amrex_real) :: dt
      end type nucleation_parm_type_input

      type user_defined_force_parm_type_input
       integer :: i,j,k
       real(amrex_real), pointer :: dx(:)
       real(amrex_real), pointer :: xlo(:)
       integer, pointer :: fablo(:)
       integer, pointer :: fabhi(:)
       real(amrex_real), pointer, dimension(D_DECL(:,:,:)) :: thermal
       real(amrex_real), pointer, dimension(D_DECL(:,:,:),:) :: uold
       real(amrex_real), pointer, dimension(D_DECL(:,:,:),:) :: lsnew
       real(amrex_real), pointer, dimension(D_DECL(:,:,:)) :: one_over_den 
       real(amrex_real) :: cur_time
       real(amrex_real) :: dt
      end type user_defined_force_parm_type_input

      type assimilate_parm_type
      integer :: level
      integer :: finest_level
      integer :: bfact
      integer :: nparts_ghost
      integer :: nparts
      integer, pointer :: im_solid_map(:)
      integer :: nstate
      integer :: nhalf
      real(amrex_real) :: cur_time
      real(amrex_real) :: dt
      real(amrex_real), pointer :: dx(:)
      real(amrex_real), pointer :: xsten(:,:)
      real(amrex_real) :: dxmin
      real(amrex_real), pointer :: xlo(:)
      integer, pointer :: fablo(:)
      integer, pointer :: fabhi(:)
      real(amrex_real), pointer, dimension(D_DECL(:,:,:),:) :: ughostx
      real(amrex_real), pointer, dimension(D_DECL(:,:,:),:) :: ughosty
      real(amrex_real), pointer, dimension(D_DECL(:,:,:),:) :: ughostz
      end type assimilate_parm_type

      type assimilate_out_parm_type
!nstate components
      real(amrex_real), pointer, dimension(D_DECL(:,:,:),:) :: state
      real(amrex_real), pointer, dimension(D_DECL(:,:,:),:) :: LS_state
      real(amrex_real), pointer, dimension(D_DECL(:,:,:)) :: macx
      real(amrex_real), pointer, dimension(D_DECL(:,:,:)) :: macy
      real(amrex_real), pointer, dimension(D_DECL(:,:,:)) :: macz
      end type assimilate_out_parm_type

       ! used by deriv_from_grid_util
      type deriv_from_grid_parm_type
      integer :: scomp
      integer :: ncomp
      integer :: level
      integer :: finest_level
      integer :: bfact
      integer :: index_flux(SDIM)  !flux point where derivative needed.
      integer :: dir_deriv !dir_deriv=1..sdim or dir_deriv=-1
      integer :: box_type_flux(SDIM) !0=CELL 1=NODE
      integer :: box_type_data(SDIM)
      integer :: grid_type_flux  ! -1..5
      integer :: grid_type_data  ! -1..5
      real(amrex_real) :: dx(SDIM)
      real(amrex_real) :: xlo(SDIM)
      integer :: fablo(SDIM)
      integer :: fabhi(SDIM)
!      real(amrex_real), INTENT(in), pointer, dimension(D_DECL(:,:,:),:) :: disp_dataptr
      end type deriv_from_grid_parm_type

       ! used by single_deriv_from_grid_util
      type single_deriv_from_grid_parm_type
      integer :: level
      integer :: finest_level
      integer :: bfact
      integer :: index_flux(SDIM)  !flux point where derivative needed.
      integer :: dir_deriv !dir_deriv=1..sdim or dir_deriv=-1
      integer :: box_type_flux(SDIM) !0=CELL 1=NODE
      integer :: box_type_data(SDIM)
      integer :: grid_type_flux  ! -1..5
      integer :: grid_type_data  ! -1..5
      real(amrex_real) :: dx(SDIM)
      real(amrex_real) :: xlo(SDIM)
      integer :: fablo(SDIM)
      integer :: fabhi(SDIM)
      end type single_deriv_from_grid_parm_type

       ! used by interp_from_grid_util
      type interp_from_grid_parm_type
      integer :: scomp
      integer :: ncomp
      integer :: level
      integer :: finest_level
      integer :: bfact
      real(amrex_real) :: xtarget(SDIM)
      real(amrex_real) :: dx(SDIM)
      real(amrex_real) :: xlo(SDIM)
      integer :: fablo(SDIM)
      integer :: fabhi(SDIM)
      end type interp_from_grid_parm_type


       ! used by single_interp_from_grid_util
      type single_interp_from_grid_parm_type
      integer :: level
      integer :: finest_level
      integer :: bfact
      real(amrex_real) :: xtarget(SDIM)
      real(amrex_real) :: dx(SDIM)
      real(amrex_real) :: xlo(SDIM)
      integer :: fablo(SDIM)
      integer :: fabhi(SDIM)
      end type single_interp_from_grid_parm_type


      type interp_from_grid_out_parm_type
      real(amrex_real), pointer :: data_interp(:)
      end type interp_from_grid_out_parm_type

     contains

end module probcommon_module_types

module probcommon_module
use amrex_fort_module, only : amrex_real
use probcommon_module_types

implicit none

! fort_material_conservation_form added June 15, 2024
! fort_im_viscoelastic_map added July 16,2024
! ngrow_distance,ngrow_make_distance added Feb 27,2025
! fort_yield_stress added March 4, 2025
! fort_mechanical_to_thermal added July 13, 2025
! fort_stiff_sound_speed added March 30, 2025
! Johnson Cook Softening parameters added April 08, 2025

      integer, PARAMETER :: MOF_TRAINING_NDIM_DECISIONS=AMREX_SPACEDIM
      integer, PARAMETER :: MOF_TRAINING_NDIM_CLASSIFY=AMREX_SPACEDIM-1

      Type training_model_type
        Type(Neural_Network) :: NN_ZHOUTENG_LOCAL
        Type(Decision_Tree) :: DT_ZHOUTENG_LOCAL
        Type(Random_Forest) :: RF_ZHOUTENG_LOCAL
      end Type training_model_type

      Type(training_model_type), allocatable, dimension(D_DECL(:,:,:),:) :: &
        training_array
      integer :: training_max_level=-1
      integer :: training_lo(SDIM)
      integer :: training_hi(SDIM)

      !https://www.stat.cmu.edu/~cshalizi/350-2006/lecture-10.pdf
      !S=sum_{c \in leaves(T)} n_{c} V_{c}
      !V_{c}=(1/n_{c})sum_{i \in C} (y_{i}-m_{c})^{2}
      Type branch_type
       integer :: ndata
       integer :: parent_id
       integer :: parent_level
       integer :: current_id
       integer :: current_level
       integer :: splittingrule
       integer :: median_index
       real(amrex_real) :: median_value
       real(amrex_real), pointer :: data_decisions(:,:) !datanum, data_idx
       real(amrex_real), pointer :: data_classify(:,:) !datanum, data_idx
       integer :: child1_id
       integer :: child_level
       integer :: child2_id
      end Type branch_type

      Type level_branch_type
        integer :: nbranches
        Type(branch_type), pointer :: branch_list(:)
      end Type level_branch_type

      Type tree_type
       integer :: max_number_tree_levels
       integer :: number_tree_levels
       integer, pointer :: nbranches_level(:) 
       Type(level_branch_type), pointer :: branch_list_level(:)
      end Type tree_type

      Type(tree_type), allocatable, dimension(D_DECL(:,:,:),:) :: &
          decision_tree_array
      integer :: decision_tree_max_level=-1
      integer :: decision_tree_lo(3)
      integer :: decision_tree_hi(3)

      integer, PARAMETER :: MAX_NUM_MATERIALS=10
      !num_interfaces=( (num_materials-1)*(num_materials-1)+num_materials-1 )/2
      integer, PARAMETER :: MAX_NUM_INTERFACES=55
      integer, PARAMETER :: MAX_NUM_SPECIES=10
       ! for user definable EOS, the EOS is defined based on "probtype"
       ! not "material_type"; but it is essential that
       ! 1<=material_type<=MAX_NUM_EOS
      integer, PARAMETER :: MAX_NUM_EOS=998

      real(amrex_real), PARAMETER :: CLAMPED_EVERYWHERE_LS=99999.0D0
      real(amrex_real), PARAMETER :: CLAMPED_NO_WHERE_LS=-99999.0D0

#include "probdataf95.H"

      integer, PARAMETER :: DEBUG_EVAPORATION=0
      integer, PARAMETER :: EVAPORATION_iter_max=50

      integer, PARAMETER :: OLD_DODECANE=1

      integer, PARAMETER :: DEBUG_DYNAMIC_CONTACT_ANGLE=0

      real(amrex_real), PARAMETER :: GNBC_RADIUS=2.0d0

      real(amrex_real), PARAMETER :: GAMMA_SIMPLE_PARMS=1.4d0

         ! R=CP-CV
         ! CP=1.007D+7 Specific heat at constant pressure cgs ergs/(Kelvin g)
         ! GAMMA=CP/CV=1.39861
         ! for air, molar_mass=28.9 g/mol =>
         !  R=(2.87E+6 erg/(K g))*(28.9 g/mol)=8.2943E+7 ergs/(Kelvin mol)
      real(amrex_real), PARAMETER :: R_AIR_PARMS=0.287D+7  ! ergs/(Kelvin g)
      real(amrex_real), PARAMETER :: CV_AIR_PARMS=0.72D+7  ! ergs/(Kelvin g)
      real(amrex_real), PARAMETER :: PCAV_TAIT=220.2726D0  ! cgs dyne/cm^2
      real(amrex_real), PARAMETER :: PCAV_ELASTIC=220.2726D0  ! cgs dyne/cm^2
      real(amrex_real), PARAMETER :: PCAV_TAIT_VACUUM=220.2726D0  ! dyne/cm^2
      real(amrex_real), PARAMETER :: A_TAIT=1.0D+6  ! dyne/cm^2
      real(amrex_real), PARAMETER :: B_TAIT=3.31D+9  ! dyne/cm^2
      real(amrex_real), PARAMETER :: RHOBAR_TAIT=1.0D0  ! g/cm^3
      real(amrex_real), PARAMETER :: GAMMA_TAIT=7.15D0

      real(amrex_real), PARAMETER :: A_GALINSTAN=1.0D+6  ! dyne/cm^2
      real(amrex_real), PARAMETER :: B_GALINSTAN=6.76D+10  ! dyne/cm^2
      real(amrex_real), PARAMETER :: GAMMA_GALINSTAN=7.15D0

      real(amrex_real), PARAMETER :: A_ELASTIC=1.0D+6  ! dyne/cm^2
      real(amrex_real), PARAMETER :: B_ELASTIC=6.76D+10  ! dyne/cm^2
      real(amrex_real), PARAMETER :: GAMMA_ELASTIC=7.15D0

      real(amrex_real), PARAMETER :: P0_tillotson=1.0D+6 ! dyne/cm^2
      real(amrex_real), PARAMETER :: a_hydro_tillotson=0.7d0  ! dimensionless
      real(amrex_real), PARAMETER :: b_hydro_tillotson=0.15d0 ! dimensionless
      real(amrex_real), PARAMETER :: rho_IV_tillotson=0.958d0 ! g/cm^3
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
      real(amrex_real), PARAMETER :: A_strain_tillotson=21.8D+9 ! dyne/cm^2
      real(amrex_real), PARAMETER :: B_strain_tillotson=132.5D+9 ! dyne/cm^2
      real(amrex_real), PARAMETER :: E0_tillotson=0.07D+12 ! erg/g
      real(amrex_real), PARAMETER :: E_IV_tillotson=0.00419D+12 ! erg/g
      real(amrex_real), PARAMETER :: E_CV_tillotson=0.025D+12   ! erg/g
      real(amrex_real), PARAMETER :: rho_cav_tillotson=0.995d0 ! g/cm^3
      real(amrex_real), PARAMETER :: P_cav_tillotson=5.0D+4 ! dyne/cm^2
      real(amrex_real), PARAMETER :: T_cav_tillotson=305.9d0 ! degrees Kelvin
      real(amrex_real), PARAMETER :: alpha_tillotson=10.0d0 
      real(amrex_real), PARAMETER :: beta_tillotson=5.0d0 

      real(amrex_real), PARAMETER :: P_cav_mie_gruneisen=1.0D+6 ! dyne/cm^2

      real(amrex_real), PARAMETER :: omega_wardlaw_tillotson=0.28d0
      real(amrex_real), PARAMETER :: A_wardlaw_tillotson=2.2D+10 ! dyne/cm^2
      real(amrex_real), PARAMETER :: B_wardlaw_tillotson=9.94D+10 ! dyne/cm^2
      real(amrex_real), PARAMETER :: C_wardlaw_tillotson=1.457D+11 ! dyne/cm^2
      real(amrex_real), PARAMETER :: &
         rho_cav_wardlaw_tillotson=0.99995681d0 !g/cm^3
      real(amrex_real), PARAMETER :: &
         sound_cav_wardlaw_tillotson=148295.1d0 !cm/s (293 Kelvin)
      integer, PARAMETER :: visual_RT_transform=1
      integer, PARAMETER :: bubbleInPackedColumn=1001

      real(amrex_real), PARAMETER :: room_temperature=293.0d0

      real(amrex_real), PARAMETER :: incomp_thickness=2.0d0

#ifdef BL_USE_FLOAT

      real(amrex_real), PARAMETER :: EVAPORATION_TOL=BL_REAL_E(1.0,-6)
      real(amrex_real), PARAMETER :: TEMPERATURE_FLOOR=BL_REAL_E(1.0,-20)
      real(amrex_real), PARAMETER :: MASK_FINEST_TOL=BL_REAL_E(1.0,-3)
      real(amrex_real), PARAMETER :: OVERFLOW_CUTOFF=BL_REAL_E(1.0,+20)
      real(amrex_real), PARAMETER :: FACETOL_REDIST=BL_REAL_E(1.0,-2)
      ! Default: LS_CURV_TOL=BL_REAL_E(1.0,-2)
      !inputs.curvature_converge with axis_dir=210 (sanity check),
      !BL_REAL_E(1.0,-12)
      real(amrex_real), PARAMETER :: LS_CURV_TOL=BL_REAL_E(1.0,-2)
      real(amrex_real), PARAMETER :: VOFTOL_SLOPES=BL_REAL_E(1.0,-2)
      real(amrex_real), PARAMETER :: VOFTOL_AREAFRAC=BL_REAL_E(1.0,-1)
      real(amrex_real), PARAMETER :: EVAP_BISECTION_TOL=BL_REAL_E(1.0,-5)
      real(amrex_real), PARAMETER :: INTERCEPT_TOL=BL_REAL_E(1.0,-6)
      real(amrex_real), PARAMETER :: FRAC_PAIR_TOL=BL_REAL_E(1.0,-5)

      real(amrex_real), PARAMETER :: TANGENT_EPS=BL_REAL_E(1.0,-2)

      real(amrex_real), PARAMETER :: EPS30=BL_REAL_E(1.0,-30)
      real(amrex_real), PARAMETER :: EPS20=BL_REAL_E(1.0,-6)
      real(amrex_real), PARAMETER :: EPS15=BL_REAL_E(1.0,-6)
      real(amrex_real), PARAMETER :: EPS14=BL_REAL_E(1.0,-6)
      real(amrex_real), PARAMETER :: EPS13=BL_REAL_E(1.0,-6)
      real(amrex_real), PARAMETER :: EPS12=BL_REAL_E(1.0,-6)
      real(amrex_real), PARAMETER :: EPS11=BL_REAL_E(1.0,-6)
      real(amrex_real), PARAMETER :: EPS10=BL_REAL_E(1.0,-6)
      real(amrex_real), PARAMETER :: EPS9=BL_REAL_E(1.0,-6)
      real(amrex_real), PARAMETER :: EPS8=BL_REAL_E(1.0,-6)
      real(amrex_real), PARAMETER :: EPS7=BL_REAL_E(1.0,-6)
      real(amrex_real), PARAMETER :: EPS6=BL_REAL_E(1.0,-6)
      real(amrex_real), PARAMETER :: EPS5=BL_REAL_E(1.0,-5)
      real(amrex_real), PARAMETER :: EPS4=BL_REAL_E(1.0,-4)
      real(amrex_real), PARAMETER :: EPS3=BL_REAL_E(1.0,-3)
      real(amrex_real), PARAMETER :: EPS2=BL_REAL_E(1.0,-2)
      real(amrex_real), PARAMETER :: EPS1=BL_REAL_E(1.0,-1)

      real(amrex_real), PARAMETER :: EPS_3_2=BL_REAL_E(1.0,-2)
      real(amrex_real), PARAMETER :: EPS_8_2=BL_REAL_E(1.0,-2)
      real(amrex_real), PARAMETER :: EPS_8_3=BL_REAL_E(1.0,-3)
      real(amrex_real), PARAMETER :: EPS_8_4=BL_REAL_E(1.0,-4)
      real(amrex_real), PARAMETER :: EPS_10_3=BL_REAL_E(1.0,-3)
      real(amrex_real), PARAMETER :: EPS_10_4=BL_REAL_E(1.0,-4)
      real(amrex_real), PARAMETER :: EPS_10_5=BL_REAL_E(1.0,-5)
      real(amrex_real), PARAMETER :: EPS_12_2=BL_REAL_E(1.0,-2)
      real(amrex_real), PARAMETER :: EPS_12_4=BL_REAL_E(1.0,-4)
      real(amrex_real), PARAMETER :: EPS_12_6=BL_REAL_E(1.0,-6)
      real(amrex_real), PARAMETER :: EPS_14_7=BL_REAL_E(1.0,-7)
      real(amrex_real), PARAMETER :: EPS_13_5=BL_REAL_E(1.0,-5)
#else

      real(amrex_real), PARAMETER :: EVAPORATION_TOL=BL_REAL_E(1.0,-10)
      real(amrex_real), PARAMETER :: TEMPERATURE_FLOOR=BL_REAL_E(1.0,-20)
      real(amrex_real), PARAMETER :: MASK_FINEST_TOL=BL_REAL_E(1.0,-3)
      real(amrex_real), PARAMETER :: OVERFLOW_CUTOFF=BL_REAL_E(1.0,+20)
      real(amrex_real), PARAMETER :: FACETOL_REDIST=BL_REAL_E(1.0,-2)
      !Default: LS_CURV_TOL=BL_REAL_E(1.0,-2) 
      !inputs.curvature_converge with axis_dir=210 (sanity check),
      !BL_REAL_E(1.0,-12)
      real(amrex_real), PARAMETER :: LS_CURV_TOL=BL_REAL_E(1.0,-2)
      real(amrex_real), PARAMETER :: VOFTOL_SLOPES=BL_REAL_E(1.0,-2)
      real(amrex_real), PARAMETER :: VOFTOL_AREAFRAC=BL_REAL_E(1.0,-1)
      real(amrex_real), PARAMETER :: EVAP_BISECTION_TOL=BL_REAL_E(1.0,-12)
      real(amrex_real), PARAMETER :: INTERCEPT_TOL=BL_REAL_E(1.0,-12)
      real(amrex_real), PARAMETER :: FRAC_PAIR_TOL=BL_REAL_E(1.0,-12)

      real(amrex_real), PARAMETER :: TANGENT_EPS=BL_REAL_E(1.0,-2)

      real(amrex_real), PARAMETER :: EPS30=BL_REAL_E(1.0,-30)
      real(amrex_real), PARAMETER :: EPS20=BL_REAL_E(1.0,-14)
      real(amrex_real), PARAMETER :: EPS15=BL_REAL_E(1.0,-14)
      real(amrex_real), PARAMETER :: EPS14=BL_REAL_E(1.0,-14)
      real(amrex_real), PARAMETER :: EPS13=BL_REAL_E(1.0,-13)
      real(amrex_real), PARAMETER :: EPS12=BL_REAL_E(1.0,-12)
      real(amrex_real), PARAMETER :: EPS11=BL_REAL_E(1.0,-11)
      real(amrex_real), PARAMETER :: EPS10=BL_REAL_E(1.0,-10)
      real(amrex_real), PARAMETER :: EPS9=BL_REAL_E(1.0,-9)
      real(amrex_real), PARAMETER :: EPS8=BL_REAL_E(1.0,-8)
      real(amrex_real), PARAMETER :: EPS7=BL_REAL_E(1.0,-7)
      real(amrex_real), PARAMETER :: EPS6=BL_REAL_E(1.0,-6)
      real(amrex_real), PARAMETER :: EPS5=BL_REAL_E(1.0,-5)
      real(amrex_real), PARAMETER :: EPS4=BL_REAL_E(1.0,-4)
      real(amrex_real), PARAMETER :: EPS3=BL_REAL_E(1.0,-3)
      real(amrex_real), PARAMETER :: EPS2=BL_REAL_E(1.0,-2)
      real(amrex_real), PARAMETER :: EPS1=BL_REAL_E(1.0,-1)

      real(amrex_real), PARAMETER :: EPS_3_2=BL_REAL_E(1.0,-3)
      real(amrex_real), PARAMETER :: EPS_8_2=BL_REAL_E(1.0,-8)
      real(amrex_real), PARAMETER :: EPS_8_3=BL_REAL_E(1.0,-8)
      real(amrex_real), PARAMETER :: EPS_8_4=BL_REAL_E(1.0,-8)
      real(amrex_real), PARAMETER :: EPS_10_3=BL_REAL_E(1.0,-10)
      real(amrex_real), PARAMETER :: EPS_10_4=BL_REAL_E(1.0,-10)
      real(amrex_real), PARAMETER :: EPS_10_5=BL_REAL_E(1.0,-10)
      real(amrex_real), PARAMETER :: EPS_12_2=BL_REAL_E(1.0,-12)
      real(amrex_real), PARAMETER :: EPS_12_4=BL_REAL_E(1.0,-12)
      real(amrex_real), PARAMETER :: EPS_12_6=BL_REAL_E(1.0,-12)
      real(amrex_real), PARAMETER :: EPS_14_7=BL_REAL_E(1.0,-14)
      real(amrex_real), PARAMETER :: EPS_13_5=BL_REAL_E(1.0,-13)
#endif
      real(amrex_real), PARAMETER :: EPS_UNCAPTURED=EPS_8_4
      real(amrex_real), PARAMETER :: VOFTOL=EPS_8_4
      real(amrex_real), PARAMETER :: EPS_FULL_WEAK=EPS2

      integer, PARAMETER :: tecplot_real=8
      integer, PARAMETER :: tecplot_real_short=4

      real(amrex_real), PARAMETER :: FSI_PRESSURE_FORCE_ONLY=1

       !SEM_IMAGE_BC_ALG=1 => use "reflect_odd" for dirichlet,
       !                      use "reflect even" for Neumann
      integer, PARAMETER :: SEM_IMAGE_BC_ALG=0

      integer, PARAMETER :: POLYGON_LIST_MAX=1000

       ! variables for comparing with Fred Stern's experiment
      integer, PARAMETER :: SHALLOW_M=100
      integer, PARAMETER :: SHALLOW_N=1000
      real(amrex_real), PARAMETER :: SHALLOW_TIME=12.0

      integer :: BLB_MATRIX=-32767
      integer :: BLB_RHS=-32767
      integer :: BLB_VEL=-32767
      integer :: BLB_INT_MOM=-32767
      integer :: BLB_ENERGY=-32767
      integer :: BLB_MASS_VEL=-32767
      integer :: BLB_VOL=-32767
      integer :: BLB_CEN_INT=-32767
      integer :: BLB_CEN_ACT=-32767
      integer :: BLB_PERIM=-32767
      integer :: BLB_PERIM_MAT=-32767
      integer :: BLB_TRIPLE_PERIM=-32767
      integer :: BLB_CELL_CNT=-32767
      integer :: BLB_CELLVOL_CNT=-32767
      integer :: BLB_MASS=-32767
      integer :: BLB_PRES=-32767
      integer :: BLB_SECONDMOMENT=-32767
      integer :: num_elements_blobclass=-32767

      real(amrex_real) shallow_water_data(0:SHALLOW_M,0:SHALLOW_N,2)
      real(amrex_real) inflow_time(3000)
      real(amrex_real) inflow_elevation(3000)
      real(amrex_real) inflow_velocity(3000)
      real(amrex_real) outflow_time(3000)
      real(amrex_real) outflow_elevation(3000)
      real(amrex_real) outflow_velocity(3000)
      integer inflow_count,outflow_count
      integer last_inflow_index,last_outflow_index

       ! variables from "rfiledata"
      real(amrex_real) zstatic(0:300),rstatic(0:300) 

       ! variables from "pressure_bcs"
      real(amrex_real)  dt_pressure_bcs
      real(amrex_real)  time_pressure_bcs(0:100) 
      real(amrex_real)  pressbc_pressure_bcs(0:100,1:3)
      integer selectpress
       ! variables from "vel_bcs"
      real(amrex_real)  timehist_velbc(0:100), &
        zpos_velbc(1:50),velbc_velbc(0:100,1:50), &
        period_velbc,rigidwall_velbc
      integer itime_velbc,ipos_velbc

        ! level,index,dir
      real(amrex_real), allocatable, dimension(:,:,:) :: grid_cache
      integer cache_index_low,cache_index_high,cache_max_level
      integer :: grid_cache_allocated=0

       ! solve: div_X (1/w(x)) grad_X x=0  x(Xlo)=Xlo   x(Xhi)=Xhi
       ! index,dir
      real(amrex_real), allocatable, dimension(:,:) :: mapping_comp_to_phys
      real(amrex_real), allocatable, dimension(:,:) :: mapping_phys_to_comp
      integer :: use_identity_mapping=1
      integer :: mapping_n_cell(0:2)
      integer :: mapping_n_cell_max=0
      integer :: mapping_allocated=0

      integer :: number_of_source_regions=0
      integer :: number_of_threads_regions=0

      type region_info_type
       integer :: region_material_id
       real(amrex_real) :: region_dt
       real(amrex_real) :: region_mass_flux   ! e.g. kg/s
       real(amrex_real) :: region_volume_flux ! e.g. m^3/s
       real(amrex_real) :: region_temperature_prescribe ! e.g. degrees Kelvin
       real(amrex_real) :: region_velocity_prescribe(SDIM) ! e.g. degrees Kelvin
       real(amrex_real) :: region_energy_flux ! e.g. J/s = Watts
       real(amrex_real) :: region_volume_raster ! e.g. m^3
       real(amrex_real) :: region_volume        ! e.g. m^3
       real(amrex_real) :: region_mass          ! e.g. kg
       real(amrex_real) :: region_energy        ! e.g. J
       real(amrex_real) :: region_energy_per_kelvin ! e.g. J/K
       real(amrex_real) :: region_volume_after
       real(amrex_real) :: region_mass_after
       real(amrex_real) :: region_energy_after
      end type region_info_type

       ! first index: 1..number_of_source_regions
       ! second index: 0..number_of_threads_regions
      type(region_info_type), allocatable :: regions_list(:,:)
      
       ! interface particle container class
       ! refined data 
      type elem_contain_type
       integer :: numElems
       integer, pointer :: ElemData(:)
      end type elem_contain_type

       ! interface particle container class
       ! original coarse data
      type node_contain_type
       integer :: numNodes
       integer, pointer :: NodeData(:)
      end type node_contain_type

      type level_contain_type
        ! tid,partid,tilenum
       type(elem_contain_type), pointer :: level_elem_data(:,:,:)
       type(node_contain_type), pointer :: level_node_data(:,:,:)
        ! tid,tilenum,dir
       integer, pointer :: tilelo3D(:,:,:)
       integer, pointer :: tilehi3D(:,:,:)
       real(amrex_real), pointer :: xlo3D(:,:,:)
        ! tid,tilenum
       integer, pointer :: gridno3D(:,:)
        ! tid
       integer, pointer :: num_tiles_on_thread3D_proc(:)
       integer :: max_num_tiles_on_thread3D_proc
       integer :: num_grids_on_level
       integer :: num_grids_on_level_proc
      end type level_contain_type

       ! level
      type(level_contain_type), dimension(:), allocatable :: contain_elem

      integer :: container_allocated=0

       ! level
      integer, allocatable, dimension(:) :: level_container_allocated

      type aux_contain_type
       integer :: lo3D(3)
       integer :: hi3D(3)
       real(amrex_real) :: xlo3D(3)
       real(amrex_real) :: xhi3D(3)
       real(amrex_real) :: dx3D
       integer :: aux_ncells_max_side
       real(amrex_real), dimension(:,:,:,:), pointer :: LS3D ! level set data
      end type aux_contain_type

      real(amrex_real), dimension(:,:,:,:), target, allocatable :: aux_xdata3D
      real(amrex_real), dimension(:,:,:,:), target, allocatable :: aux_FSIdata3D
      real(amrex_real), dimension(:,:,:,:), target, allocatable :: aux_masknbr3D

      integer :: fort_num_local_aux_grids=0
      integer :: aux_data_allocated=0
      type(aux_contain_type), dimension(:), allocatable :: contain_aux

      integer :: used_probtypes(1000)
      integer :: probtype_list_size

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
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: dir ! 1,2,3
      integer, INTENT(in) :: data_dir ! 0,1,2
      real(amrex_real), INTENT(in) :: dxmin
      real(amrex_real), INTENT(in), pointer :: x_projection_raster(:)
      real(amrex_real), INTENT(in), pointer :: dx(:)
      real(amrex_real), INTENT(in), pointer :: n_raster(:) ! points to solid
      integer, INTENT(in) :: im_fluid
      real(amrex_real), INTENT(in) :: u !uimage_raster_solid_frame(dir)
      real(amrex_real), INTENT(in) :: uimage_tngt_mag
      real(amrex_real), INTENT(in) :: wall_model_velocity
      real(amrex_real), INTENT(in) :: dist_probe
      real(amrex_real), INTENT(in) :: dist_fluid
      real(amrex_real), INTENT(in) :: temperature_image
      real(amrex_real), INTENT(in) :: temperature_wall
      real(amrex_real), INTENT(in) :: temperature_wall_max
      real(amrex_real), INTENT(in) :: viscosity_molecular
      real(amrex_real), INTENT(in) :: viscosity_eddy_wall
      real(amrex_real), INTENT(in) :: y !delta_r
      real(amrex_real), INTENT(in) :: critical_length
      real(amrex_real), INTENT(out) :: ughost_tngt  ! dir direction
      end subroutine TEMPLATE_wallfunc

       ! returns (1/w) where w>>1 in "trouble" regions
      subroutine TEMPLATE_MAPPING_WEIGHT_COEFF(dir,wt,phys_x)
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: dir
      real(amrex_real), INTENT(out) :: wt
      real(amrex_real), INTENT(in) :: phys_x
      end subroutine TEMPLATE_MAPPING_WEIGHT_COEFF

      subroutine TEMPLATE_INIT_REGIONS_LIST(constant_density_all_time, &
          num_materials_in,num_threads_in)
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: num_materials_in
      integer, INTENT(in) :: num_threads_in
      integer, INTENT(in) :: constant_density_all_time(num_materials_in)
      end subroutine TEMPLATE_INIT_REGIONS_LIST

      subroutine TEMPLATE_CHARFN_REGION(region_id,x,cur_time,charfn_out)
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: region_id
      real(amrex_real), INTENT(in) :: x(SDIM)
      real(amrex_real), INTENT(in) :: cur_time
      real(amrex_real), INTENT(out) :: charfn_out
      end subroutine TEMPLATE_CHARFN_REGION

      subroutine TEMPLATE_T0_Boussinesq(x,dx,cur_time,im,T0)
      use amrex_fort_module, only : amrex_real
      real(amrex_real), INTENT(in) :: x(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: cur_time
      integer, INTENT(in) :: im
      real(amrex_real), INTENT(out) :: T0
      end subroutine TEMPLATE_T0_Boussinesq

      subroutine TEMPLATE_V0_Coriolis(x,dx,cur_time,V0)
      use amrex_fort_module, only : amrex_real
      real(amrex_real), INTENT(in) :: x(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: cur_time
      real(amrex_real), INTENT(out) :: V0(SDIM)
      end subroutine TEMPLATE_V0_Coriolis

      subroutine TEMPLATE_angular_velocity_vector( &
       x,cur_time, &
       angular_velocity_vector, &
       angular_velocity_vector_custom, &
       angular_velocity_vector_dot, &
       lever_arm_in, &
       lever_arm_custom)
      use amrex_fort_module, only : amrex_real
      real(amrex_real), INTENT(in) :: x(SDIM)
      real(amrex_real), INTENT(in) :: cur_time
      real(amrex_real), INTENT(in) :: angular_velocity_vector(3)
      real(amrex_real), INTENT(in) :: lever_arm_in(SDIM)
      real(amrex_real), INTENT(out) :: angular_velocity_vector_custom(3)
      real(amrex_real), INTENT(out) :: lever_arm_custom(SDIM)
      real(amrex_real), INTENT(out) :: angular_velocity_vector_dot(3)
      end subroutine TEMPLATE_angular_velocity_vector

      subroutine TEMPLATE_gravity_vector(x,cur_time, &
       gravity_vector_in,gravity_vector_out)
      use amrex_fort_module, only : amrex_real
      real(amrex_real), INTENT(in) :: x(SDIM)
      real(amrex_real), INTENT(in) :: cur_time
      real(amrex_real), INTENT(in) :: gravity_vector_in(SDIM)
      real(amrex_real), INTENT(out) :: gravity_vector_out(SDIM)
      end subroutine TEMPLATE_gravity_vector

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
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: im
      integer, INTENT(in) :: im_solid
      integer, INTENT(in) :: near_interface
      real(amrex_real), INTENT(in) :: x(SDIM)
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: cur_time
      real(amrex_real), INTENT(in) :: density
      real(amrex_real), INTENT(in) :: temperature
      real(amrex_real), INTENT(in) :: temperature_wall
      real(amrex_real), INTENT(in) :: temperature_wall_max
      real(amrex_real), INTENT(in) :: temperature_probe
      real(amrex_real), INTENT(in) :: nrm(SDIM) ! nrm points from solid to fluid
      real(amrex_real), INTENT(inout) :: thermal_k
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
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: interface_mass_transfer_model
      integer, INTENT(in) :: num_materials_in
      integer, INTENT(in) :: num_species_var_in
      integer, INTENT(in) :: probe_constrain
      integer, INTENT(in) :: ireverse
      integer, INTENT(in) :: iten
      real(amrex_real), INTENT(in) :: xI(SDIM)
      real(amrex_real), INTENT(in) :: cur_time
      real(amrex_real), INTENT(in) :: prev_time
      real(amrex_real), INTENT(in) :: dt
      real(amrex_real), INTENT(inout) :: TI
      real(amrex_real), INTENT(inout) :: YI
      integer, INTENT(inout) :: user_override_TI_YI
      real(amrex_real), INTENT(in) :: molar_mass(num_materials_in)
      real(amrex_real), INTENT(in) :: species_molar_mass(num_species_var_in)
      real(amrex_real), INTENT(in) :: ksrc_predict
      real(amrex_real), INTENT(in) :: kdst_predict
      real(amrex_real), INTENT(in) :: ksrc_physical
      real(amrex_real), INTENT(in) :: kdst_physical
      real(amrex_real), INTENT(in) :: T_probe_src
      real(amrex_real), INTENT(in) :: T_probe_dst
      real(amrex_real), INTENT(in) :: LL
      real(amrex_real), INTENT(in) :: dxprobe_src
      real(amrex_real), INTENT(in) :: dxprobe_dst
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
        interp_valid_flag_src, &
        interp_valid_flag_dst, &
        T_probe_src, &
        T_probe_dst, &
        TI, &
        LL, &
        dxprobe_src, &
        dxprobe_dst)
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: interface_mass_transfer_model
      integer, INTENT(in) :: num_materials_in
      integer, INTENT(in) :: num_species_var_in
      integer, INTENT(in) :: ispec
      integer, INTENT(in) :: im_source
      integer, INTENT(in) :: im_dest
      real(amrex_real), INTENT(in) :: xI(SDIM)
      real(amrex_real), INTENT(in) :: TI
      real(amrex_real), INTENT(in) :: molar_mass(num_materials_in)
      real(amrex_real), INTENT(in) :: species_molar_mass(num_species_var_in+1)
      real(amrex_real), INTENT(out) :: mdot
      integer, INTENT(inout) :: mdot_override
      real(amrex_real), INTENT(in) :: ksrc_derived
      real(amrex_real), INTENT(in) :: kdst_derived
      real(amrex_real), INTENT(in) :: ksrc_physical
      real(amrex_real), INTENT(in) :: kdst_physical
      integer, INTENT(in) :: interp_valid_flag_src
      integer, INTENT(in) :: interp_valid_flag_dst
      real(amrex_real), INTENT(in) :: T_probe_src
      real(amrex_real), INTENT(in) :: T_probe_dst
      real(amrex_real), INTENT(in) :: LL
      real(amrex_real), INTENT(in) :: dxprobe_src
      real(amrex_real), INTENT(in) :: dxprobe_dst
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
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: interface_mass_transfer_model
      integer, INTENT(in) :: num_materials_in
      integer, INTENT(in) :: num_species_var_in
      integer, INTENT(in) :: ireverse
      integer, INTENT(in) :: iten
      real(amrex_real), INTENT(in) :: molar_mass(num_materials_in)
      real(amrex_real), INTENT(in) :: species_molar_mass(num_species_var_in)
      real(amrex_real), INTENT(in) :: k_model_predict(2) ! src,dst
      real(amrex_real), INTENT(inout) :: k_model_correct(2) ! src,dst
      real(amrex_real), INTENT(in) :: k_physical_base(2) ! src, dst
      real(amrex_real), INTENT(in) :: T_probe_src
      real(amrex_real), INTENT(in) :: T_probe_dst
      real(amrex_real), INTENT(in) :: LL
      real(amrex_real), INTENT(in) :: dxprobe_src
      real(amrex_real), INTENT(in) :: dxprobe_dst
      end subroutine TEMPLATE_K_EFFECTIVE


      subroutine TEMPLATE_reference_wavelen(wavelen)
      use amrex_fort_module, only : amrex_real
      real(amrex_real), INTENT(inout) :: wavelen
      end subroutine TEMPLATE_reference_wavelen

      subroutine TEMPLATE_DELETE_REGIONS_LIST()
      end subroutine TEMPLATE_DELETE_REGIONS_LIST

      subroutine TEMPLATE_INIT_MODULE()
      end subroutine TEMPLATE_INIT_MODULE

      subroutine TEMPLATE_DEALLOCATE_MODULE()
      end subroutine TEMPLATE_DEALLOCATE_MODULE

      subroutine TEMPLATE_correct_pres_rho_hydrostatic( &
        i,j,k,level, &
        angular_velocity_vector, &
        lever_arm_in, &
        centrifugal_force_factor, &
        dt, &
        rho_hydrostatic, &
        pres_hydrostatic, &
        state_ptr)
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: i,j,k,level
      real(amrex_real), INTENT(in) :: angular_velocity_vector(3)
      real(amrex_real), INTENT(in) :: lever_arm_in(SDIM)
      real(amrex_real), INTENT(in) :: centrifugal_force_factor
      real(amrex_real), INTENT(in) :: dt
      real(amrex_real), INTENT(inout) :: rho_hydrostatic
      real(amrex_real), INTENT(inout) :: pres_hydrostatic
      real(amrex_real), INTENT(in),pointer :: state_ptr(D_DECL(:,:,:),:)
      end subroutine TEMPLATE_correct_pres_rho_hydrostatic
        
      subroutine TEMPLATE_CFL_HELPER(time,dir,uu,dx)
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: dir
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(inout) :: uu
      real(amrex_real), INTENT(in) :: dx(SDIM)
      end subroutine TEMPLATE_CFL_HELPER

      subroutine TEMPLATE_USER_DEFINED_FORCE(xpoint,output_force, &
                      force_input)
      use amrex_fort_module, only : amrex_real
      use probcommon_module_types

      type(user_defined_force_parm_type_input), INTENT(in) :: force_input
      real(amrex_real), INTENT(in) :: xpoint(AMREX_SPACEDIM)
      real(amrex_real), INTENT(out) :: output_force(AMREX_SPACEDIM)
      end subroutine TEMPLATE_USER_DEFINED_FORCE

      subroutine TEMPLATE_SUMINT(GRID_DATA_IN,increment_out1, &
            increment_out2,nsum1,nsum2,isweep)
      use amrex_fort_module, only : amrex_real
      use probcommon_module_types

      integer, INTENT(in) :: nsum1,nsum2,isweep
      type(user_defined_sum_int_type), INTENT(in) :: GRID_DATA_IN
      real(amrex_real), INTENT(inout) :: increment_out1(nsum1)
      real(amrex_real), INTENT(inout) :: increment_out2(nsum2)
      end subroutine TEMPLATE_SUMINT
      
      subroutine TEMPLATE_LS(x,t,LS,nmat)
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: nmat
      real(amrex_real), INTENT(in) :: x(SDIM)
      real(amrex_real), INTENT(in) :: t
      real(amrex_real), INTENT(out) :: LS(nmat)
      end subroutine TEMPLATE_LS

      subroutine TEMPLATE_OVERRIDE_TAGFLAG( &
        i,j,k, &
        level,max_level, &
        snew_ptr,lsnew_ptr, &
        xsten,nhalf,time, &
        rflag,tagflag)
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: i,j,k
      integer, INTENT(in) :: level,max_level
      integer, INTENT(in) :: nhalf
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(inout) :: rflag
      integer, INTENT(inout) :: tagflag
      real(amrex_real), INTENT(in),pointer :: snew_ptr(D_DECL(:,:,:),:)
      real(amrex_real), INTENT(in),pointer :: lsnew_ptr(D_DECL(:,:,:),:)
      end subroutine TEMPLATE_OVERRIDE_TAGFLAG

      subroutine TEMPLATE_AUX_DATA(auxcomp,x,LS)
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: auxcomp
      real(amrex_real), INTENT(in) :: x(3)
      real(amrex_real), INTENT(out) :: LS
      end subroutine TEMPLATE_AUX_DATA

      subroutine TEMPLATE_OVERRIDE_FSI_SIGN_LS_VEL_TEMP( &
        exterior_BB, &
        interior_BB, &
        xcell,time,LS,VEL,TEMP,MASK,lev77,im_part,part_id)
      use amrex_fort_module, only : amrex_real
      real(amrex_real), INTENT(in) :: exterior_BB(3,2)
      real(amrex_real), INTENT(in) :: interior_BB(3,2)
      real(amrex_real), INTENT(in) :: xcell(3)
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(out) :: LS
      real(amrex_real), INTENT(out) :: VEL(3)
      real(amrex_real), INTENT(out) :: TEMP
      integer, INTENT(out) :: MASK
      integer, INTENT(in) :: lev77 !lev77=-1 for aux, >=0 otherwise.
      integer, INTENT(in) :: im_part
      integer, INTENT(in) :: part_id
      end subroutine TEMPLATE_OVERRIDE_FSI_SIGN_LS_VEL_TEMP

      subroutine TEMPLATE_GET_OUTSIDE_POINT( &
        exterior_BB, &
        xcell,time,x_outside,im_part,part_id)
      use amrex_fort_module, only : amrex_real
      real(amrex_real), INTENT(in) :: exterior_BB(3,2)
      real(amrex_real), INTENT(in) :: xcell(3)
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(out) :: x_outside(3)
      integer, INTENT(in) :: im_part
      integer, INTENT(in) :: part_id
      end subroutine TEMPLATE_GET_OUTSIDE_POINT

      subroutine TEMPLATE_BOUNDING_BOX_AUX(auxcomp, &
          minnode,maxnode,LS_FROM_SUBROUTINE,aux_ncells_max_side)
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: auxcomp
      real(amrex_real), INTENT(inout) :: minnode(3)
      real(amrex_real), INTENT(inout) :: maxnode(3)
      integer, INTENT(out) :: LS_FROM_SUBROUTINE
      integer, INTENT(out) :: aux_ncells_max_side
      end subroutine TEMPLATE_BOUNDING_BOX_AUX


      subroutine TEMPLATE_check_vel_rigid(x,t,vel,dir)
       use amrex_fort_module, only : amrex_real
       real(amrex_real), INTENT(in) :: x(SDIM)
       real(amrex_real), INTENT(in) :: t
       real(amrex_real), INTENT(in) :: vel
       integer, INTENT(in) :: dir
      end subroutine TEMPLATE_check_vel_rigid

      subroutine TEMPLATE_verification_flag(verification_flag)
      integer, INTENT(out) :: verification_flag
      end subroutine TEMPLATE_verification_flag

      subroutine TEMPLATE_clamped_LS(x,t,LS,vel,temperature,prescribed_flag,dx)
       use amrex_fort_module, only : amrex_real
       real(amrex_real), INTENT(in) :: x(SDIM)
       real(amrex_real), INTENT(in) :: dx(SDIM)
       real(amrex_real), INTENT(in) :: t
       real(amrex_real), INTENT(out) :: LS
       real(amrex_real), INTENT(out) :: vel(SDIM)
       real(amrex_real), INTENT(out) :: temperature
       integer, INTENT(out) :: prescribed_flag
      end subroutine TEMPLATE_clamped_LS

      subroutine TEMPLATE_VEL(x,t,LS,VEL,velsolid_flag,dx,nmat)
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: nmat
      real(amrex_real), INTENT(in) :: x(SDIM)
      real(amrex_real), INTENT(in) :: t
      real(amrex_real), INTENT(in) :: dx(SDIM)
      real(amrex_real), INTENT(in) :: LS(nmat)
      real(amrex_real), INTENT(out) :: VEL(SDIM)
      integer, INTENT(in) :: velsolid_flag
      end subroutine TEMPLATE_VEL

      subroutine TEMPLATE_EOS(rho,massfrac_var, &
        internal_energy,pressure, &
        imattype,im,num_species_var_in)
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: imattype,im,num_species_var_in
      real(amrex_real), INTENT(in) :: rho
      real(amrex_real), INTENT(in) :: massfrac_var(num_species_var_in+1)
      real(amrex_real), INTENT(in) :: internal_energy
      real(amrex_real), INTENT(out) :: pressure
      end subroutine TEMPLATE_EOS

      subroutine TEMPLATE_UNITLESS_EXPANSION_FACTOR( &
        im,temperature,temperature_base,expansion_factor)
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: im
      real(amrex_real), INTENT(in) :: temperature
      real(amrex_real), INTENT(in) :: temperature_base
      real(amrex_real), INTENT(out) :: expansion_factor
      end subroutine TEMPLATE_UNITLESS_EXPANSION_FACTOR

      subroutine TEMPLATE_INTERNAL_GRAVITY_WAVE_FLAG( &
        internal_wave_exists)
      integer, INTENT(out) :: internal_wave_exists
      end subroutine TEMPLATE_INTERNAL_GRAVITY_WAVE_FLAG

      subroutine TEMPLATE_dVdT(dVdT,massfrac_var, &
        pressure,temperature, &
        imattype,im,num_species_var_in)
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: imattype,im,num_species_var_in
      real(amrex_real), INTENT(in) :: pressure,temperature
      real(amrex_real), INTENT(in) :: massfrac_var(num_species_var_in+1)
      real(amrex_real), INTENT(out) :: dVdT
      end subroutine TEMPLATE_dVdT


      subroutine TEMPLATE_SOUNDSQR(rho,massfrac_var, &
        internal_energy,soundsqr, &
        imattype,im,num_species_var_in)
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: imattype,im,num_species_var_in
      real(amrex_real), INTENT(in) :: rho
      real(amrex_real), INTENT(in) :: massfrac_var(num_species_var_in+1)
      real(amrex_real), INTENT(in) :: internal_energy
      real(amrex_real), INTENT(out) :: soundsqr
      end subroutine TEMPLATE_SOUNDSQR

      subroutine TEMPLATE_INTERNAL(rho,massfrac_var, &
        temperature,local_internal_energy, &
        imattype,im,num_species_var_in)
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: imattype,im,num_species_var_in
      real(amrex_real), INTENT(in) :: rho
      real(amrex_real), INTENT(in) :: massfrac_var(num_species_var_in+1)
      real(amrex_real), INTENT(in) :: temperature 
      real(amrex_real), INTENT(out) :: local_internal_energy
      end subroutine TEMPLATE_INTERNAL

      subroutine TEMPLATE_TEMPERATURE(rho,massfrac_var, &
        temperature,internal_energy, &
        imattype,im,num_species_var_in)
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: imattype,im,num_species_var_in
      real(amrex_real), INTENT(in) :: rho
      real(amrex_real), INTENT(in) :: massfrac_var(num_species_var_in+1)
      real(amrex_real), INTENT(out) :: temperature 
      real(amrex_real), INTENT(in) :: internal_energy
      end subroutine TEMPLATE_TEMPERATURE

      subroutine TEMPLATE_PRES(x,t,LS,PRES,nmat)
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: nmat
      real(amrex_real), INTENT(in) :: x(SDIM)
      real(amrex_real), INTENT(in) :: t
      real(amrex_real), INTENT(in) :: LS(nmat)
      real(amrex_real), INTENT(out) :: PRES
      end subroutine TEMPLATE_PRES

      subroutine TEMPLATE_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: bcflag
      integer, INTENT(in) :: nmat
      integer, INTENT(in) :: nstate_mat
      real(amrex_real), INTENT(in) :: x(SDIM)
      real(amrex_real), INTENT(in) :: t
      real(amrex_real), INTENT(in) :: LS(nmat)
      real(amrex_real), INTENT(out) :: STATE(nmat*nstate_mat)
      end subroutine TEMPLATE_STATE

      subroutine TEMPLATE_LS_BC(xwall,xghost,t,LS, &
       LS_in,dir,side,dx,nmat)
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: nmat
      real(amrex_real), INTENT(in) :: xwall
      real(amrex_real), INTENT(in) :: xghost(SDIM)
      real(amrex_real), INTENT(in) :: t
      real(amrex_real), INTENT(inout) :: LS(nmat)
      real(amrex_real), INTENT(in) :: LS_in(nmat)
      integer, INTENT(in) :: dir,side
      real(amrex_real), INTENT(in) :: dx(SDIM)
      end subroutine TEMPLATE_LS_BC

      subroutine TEMPLATE_VEL_BC(xwall,xghost,t,LS, &
        VEL,VEL_in,veldir,dir,side,dx,nmat)
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: nmat
      real(amrex_real), INTENT(in) :: xwall
      real(amrex_real), INTENT(in) :: xghost(SDIM)
      real(amrex_real), INTENT(in) :: t
      real(amrex_real), INTENT(in) :: LS(nmat)
      real(amrex_real), INTENT(inout) :: VEL
      real(amrex_real), INTENT(in) :: VEL_in
      integer, INTENT(in) :: veldir,dir,side
      real(amrex_real), INTENT(in) :: dx(SDIM)
      end subroutine TEMPLATE_VEL_BC

      subroutine TEMPLATE_PRES_BC(xwall,xghost,t,LS, &
        PRES,PRES_in,dir,side,dx,nmat)
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: nmat
      real(amrex_real), INTENT(in) :: xwall
      real(amrex_real), INTENT(in) :: xghost(SDIM)
      real(amrex_real), INTENT(in) :: t
      real(amrex_real), INTENT(in) :: LS(nmat)
      real(amrex_real), INTENT(inout) :: PRES
      real(amrex_real), INTENT(in) :: PRES_in
      integer, INTENT(in) :: dir,side
      real(amrex_real), INTENT(in) :: dx(SDIM)
      end subroutine TEMPLATE_PRES_BC

      subroutine TEMPLATE_STATE_BC(xwall,xghost,t,LS, &
       STATE,STATE_merge,STATE_in,im,istate,dir,side,dx, &
       nmat)
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: nmat
      real(amrex_real), INTENT(in) :: xwall
      real(amrex_real), INTENT(in) :: xghost(SDIM)
      real(amrex_real), INTENT(in) :: t
      real(amrex_real), INTENT(in) :: LS(nmat)
      real(amrex_real), INTENT(inout) :: STATE
      real(amrex_real), INTENT(inout) :: STATE_merge
      real(amrex_real), INTENT(in) :: STATE_in
      integer, INTENT(in) :: dir,side
      real(amrex_real), INTENT(in) :: dx(SDIM)
      integer, INTENT(in) :: istate,im
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
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: nmat
      integer, INTENT(in) :: im
      real(amrex_real), INTENT(in) :: VFRAC(nmat)
      real(amrex_real), INTENT(in) :: time
      integer, INTENT(in) :: nhalf
      real(amrex_real), INTENT(in) :: x(SDIM)
      real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
      real(amrex_real), INTENT(in) :: temp(nmat)
      real(amrex_real), INTENT(in) :: den(nmat)
      real(amrex_real), INTENT(in) :: CV(nmat)
      real(amrex_real), INTENT(in) :: dt
      real(amrex_real), INTENT(out) :: heat_source
      end subroutine TEMPLATE_HEATSOURCE

      subroutine TEMPLATE_EB_heat_source(time,dt,xsten,nhalf, &
        heat_flux,heat_dir,heat_side)
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: nhalf
      real(amrex_real), dimension(-nhalf:nhalf,SDIM), INTENT(in) :: xsten
      real(amrex_real), INTENT(in) :: time
      real(amrex_real), INTENT(in) :: dt
      real(amrex_real), INTENT(out) :: heat_flux
      integer, INTENT(out) :: heat_dir
      integer, INTENT(out) :: heat_side
      end subroutine TEMPLATE_EB_heat_source

      subroutine TEMPLATE_velfreestream(problen,local_buffer)
      use amrex_fort_module, only : amrex_real
      real(amrex_real), INTENT(inout) :: local_buffer(2*SDIM)
      real(amrex_real), INTENT(in)    :: problen(SDIM)
      end subroutine TEMPLATE_velfreestream

      subroutine TEMPLATE_nucleation(nucleate_in,xsten,nhalf,make_seed)
      use amrex_fort_module, only : amrex_real
      use probcommon_module_types
      integer, INTENT(in) :: nhalf
      real(amrex_real), dimension(-nhalf:nhalf,SDIM), INTENT(in) :: xsten
      integer, INTENT(inout) :: make_seed
      type(nucleation_parm_type_input), INTENT(in) :: nucleate_in
      end subroutine TEMPLATE_nucleation

      subroutine TEMPLATE_ICE_SUBSTRATE_DISTANCE(xtarget,dist)
      use amrex_fort_module, only : amrex_real
      real(amrex_real), INTENT(in) :: xtarget(SDIM)
      real(amrex_real), INTENT(out) :: dist
      end subroutine TEMPLATE_ICE_SUBSTRATE_DISTANCE


      subroutine TEMPLATE_microcell_heat_coeff(heatcoeff,dx,veldir)
      use amrex_fort_module, only : amrex_real
      real(amrex_real), INTENT(in) :: dx(SDIM)
      integer, INTENT(in) :: veldir
      real(amrex_real), INTENT(inout) :: heatcoeff
      end subroutine TEMPLATE_microcell_heat_coeff

      subroutine TEMPLATE_ASSIMILATE(assimilate_in,assimilate_out, &
         i,j,k,cell_flag)
      use amrex_fort_module, only : amrex_real
      use probcommon_module_types
      type(assimilate_parm_type), INTENT(in) :: assimilate_in
      type(assimilate_out_parm_type), INTENT(inout) :: assimilate_out
      integer, INTENT(in) :: i,j,k,cell_flag
      end subroutine TEMPLATE_ASSIMILATE

      subroutine TEMPLATE_INIT_EVAL( &
        i,j,k,dir, &
        xpoint,cur_time, &
        scomp_size, &
        ncomp_size, &
        State_Type, &
        LS_Type, &
        DIV_Type, &
        Solid_State_Type, &
        Tensor_Type, &
        Refine_Density_Type, &
        ncomp_total, &
        scomp_array, &
        ncomp_array, &
        local_cell_evec, &
        local_velx, &
        local_vely, &
        local_velz)
      use amrex_fort_module, only : amrex_real

      integer, INTENT(in) :: i,j,k,dir
      real(amrex_real), INTENT(in) :: xpoint(SDIM)
      real(amrex_real), INTENT(in) :: cur_time
      integer, INTENT(in) :: scomp_size
      integer, INTENT(in) :: ncomp_size
      integer, INTENT(in) :: State_Type
      integer, INTENT(in) :: LS_Type
      integer, INTENT(in) :: DIV_Type
      integer, INTENT(in) :: Solid_State_Type
      integer, INTENT(in) :: Tensor_Type
      integer, INTENT(in) :: Refine_Density_Type
      integer, INTENT(in) :: ncomp_total
      integer, INTENT(in) :: scomp_array(scomp_size)
      integer, INTENT(in) :: ncomp_array(ncomp_size)
      real(amrex_real), INTENT(inout) :: local_cell_evec(ncomp_total)
      real(amrex_real), INTENT(inout) :: local_velx
      real(amrex_real), INTENT(inout) :: local_vely
      real(amrex_real), INTENT(inout) :: local_velz
      end subroutine TEMPLATE_INIT_EVAL


      subroutine TEMPLATE_FSI_SLICE(xmap3D,xslice3D,problo3D,probhi3D,dx_slice)
      use amrex_fort_module, only : amrex_real
      real(amrex_real), INTENT(in) :: dx_slice
      integer, INTENT(inout) :: xmap3D(3)
      real(amrex_real), INTENT(inout) :: xslice3D(3)
      real(amrex_real), INTENT(out) :: problo3D(3)
      real(amrex_real), INTENT(out) :: probhi3D(3)
      end subroutine TEMPLATE_FSI_SLICE

      subroutine TEMPLATE_OPEN_CASFILE(part_id,unit_id,file_format)
      integer, INTENT(in) :: part_id
      integer, INTENT(in) :: unit_id
      integer, INTENT(out) :: file_format
      end subroutine TEMPLATE_OPEN_CASFILE

      subroutine TEMPLATE_OPEN_AUXFILE(part_id,unit_id,file_format)
      integer, INTENT(in) :: part_id
      integer, INTENT(in) :: unit_id
      integer, INTENT(out) :: file_format
      end subroutine TEMPLATE_OPEN_AUXFILE

      subroutine TEMPLATE_ORDER_NODES(nodes,nodemap)
      use amrex_fort_module, only : amrex_real
      real(amrex_real), INTENT(in) :: nodes(3,3) ! dir,nodenum
      integer, INTENT(inout) :: nodemap(3)
      end subroutine TEMPLATE_ORDER_NODES

      subroutine TEMPLATE_VARIABLE_SURFACE_TENSION( &
        xpos, &
        time, &
        iten, &
        temperature, &
        tension)
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: iten
      real(amrex_real), INTENT(in) :: time,temperature
      real(amrex_real), INTENT(in) :: xpos(SDIM)
      real(amrex_real), INTENT(inout) :: tension
      end subroutine TEMPLATE_VARIABLE_SURFACE_TENSION

      subroutine TEMPLATE_VARIABLE_LATENT_HEAT( &
        iten, &
        temperature, &
        latent_heat)
      use amrex_fort_module, only : amrex_real
      integer, INTENT(in) :: iten
      real(amrex_real), INTENT(in) :: temperature
      real(amrex_real), INTENT(inout) :: latent_heat ! always positive
      end subroutine TEMPLATE_VARIABLE_LATENT_HEAT


      END INTERFACE

      PROCEDURE(TEMPLATE_INIT_MODULE), POINTER :: SUB_INIT_MODULE
      PROCEDURE(TEMPLATE_DEALLOCATE_MODULE), POINTER :: SUB_DEALLOCATE_MODULE
      PROCEDURE(TEMPLATE_correct_pres_rho_hydrostatic), POINTER :: &
              SUB_correct_pres_rho_hydrostatic
      PROCEDURE(TEMPLATE_CFL_HELPER), POINTER :: SUB_CFL_HELPER
      PROCEDURE(TEMPLATE_SUMINT), POINTER :: SUB_SUMINT
      PROCEDURE(TEMPLATE_USER_DEFINED_FORCE), POINTER :: SUB_USER_DEFINED_FORCE
      PROCEDURE(TEMPLATE_LS), POINTER :: SUB_LS
      PROCEDURE(TEMPLATE_OVERRIDE_TAGFLAG), POINTER :: SUB_OVERRIDE_TAGFLAG
      PROCEDURE(TEMPLATE_AUX_DATA), POINTER :: SUB_AUX_DATA
      PROCEDURE(TEMPLATE_OVERRIDE_FSI_SIGN_LS_VEL_TEMP), POINTER :: &
        SUB_OVERRIDE_FSI_SIGN_LS_VEL_TEMP
      PROCEDURE(TEMPLATE_GET_OUTSIDE_POINT), POINTER :: &
        SUB_GET_OUTSIDE_POINT
      PROCEDURE(TEMPLATE_BOUNDING_BOX_AUX), POINTER :: SUB_BOUNDING_BOX_AUX
      PROCEDURE(TEMPLATE_check_vel_rigid), POINTER :: SUB_check_vel_rigid
      PROCEDURE(TEMPLATE_verification_flag), POINTER :: SUB_verification_flag
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
      PROCEDURE(TEMPLATE_INIT_EVAL), POINTER :: SUB_INIT_EVAL

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
      PROCEDURE(TEMPLATE_angular_velocity_vector), POINTER ::  &
              SUB_angular_velocity_vector
      PROCEDURE(TEMPLATE_gravity_vector), POINTER :: SUB_gravity_vector

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

      real(amrex_real), INTENT(in) :: rho
      real(amrex_real), INTENT(out) :: pressure
      real(amrex_real) A,B,rhobar,pcav


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

      real(amrex_real) rho,pressure
      real(amrex_real) A,B,rhobar,GAMMA,pcav

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

