#undef BL_LANG_CC
#define BL_LANG_FORT

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

      type user_defined_sum_int_type
       INTEGER_T ncomp_sum_int_user
       REAL_T, pointer :: problo(:)      
       REAL_T, pointer :: probhi(:) 
       INTEGER_T :: igrid,jgrid,kgrid
       REAL_T :: volgrid
       INTEGER_T :: nhalf
       INTEGER_T :: nmat
       INTEGER_T :: bfact
       INTEGER_T :: ntensorMM
       INTEGER_T :: den_ncomp
       INTEGER_T, pointer :: tilelo(:)
       INTEGER_T, pointer :: tilehi(:)
       INTEGER_T, pointer :: fablo(:)
       INTEGER_T, pointer :: fabhi(:)
       REAL_T, pointer :: xlo(:)
       REAL_T, pointer :: dx(:)
       REAL_T, pointer :: xsten(:,:)
       REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: cellten
       REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: lsfab
       REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: slopes
       REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: den
       REAL_T, pointer, dimension(D_DECL(:,:,:),:) :: vel
      end type user_defined_sum_int_type

     contains

end module probcommon_module_types

module probcommon_module
use probcommon_module_types

implicit none

! fort_cavdenconst added: December 25, 2018
! fort_elastic_viscosity added: August 22, 2018
! fort_im_elastic_map added: August 22, 2018
! num_materials_elastic deleted: August 21, 2018
! fort_viscconst_eddy added: July 15, 2018
! fort_initial_temperature added: April 10, 2018
! nucleation_init_time added: May 5, 2018
! fort_density_floor and fort_density_ceiling added: January 8, 2019
! fort_ZEYU_DCA_SELECT added: January 20, 2020

#include "probdataf95.H"

      INTEGER_T, PARAMETER :: OLD_DODECANE=1

      INTEGER_T, PARAMETER :: DEBUG_DYNAMIC_CONTACT_ANGLE=1

      REAL_T, PARAMETER :: GNBC_RADIUS=2.0d0

      INTEGER_T, PARAMETER :: ngrow_make_distance=3

      INTEGER_T, PARAMETER :: FORT_NUM_TENSOR_TYPE=2*SDIM

      REAL_T, PARAMETER :: GAMMA_SIMPLE_PARMS=1.4

         ! R=CP-CV
         ! CP=1.007D+7 Specific heat at constant pressure cgs ergs/(Kelvin g)
         ! GAMMA=CP/CV=1.39861
      REAL_T, PARAMETER :: R_AIR_PARMS=0.287D+7  ! ergs/(Kelvin g)
      REAL_T, PARAMETER :: R_Palmore_Desjardins=0.287D+7  ! ergs/(Kelvin g)
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
      INTEGER_T, PARAMETER :: DO_SANITY_CHECK=0
      INTEGER_T, PARAMETER :: COARSE_FINE_VELAVG=1
      REAL_T, PARAMETER :: MASK_FINEST_TOL=1.0D-3

      REAL_T, PARAMETER :: ICEFACECUT_EPS=1.0D-5

       ! Default: FACETOL_DVOL=1.0D-3 
       ! Default: FACETOL_DVOL=1.0D-6 (prototype code)
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
       ! Default: VOFTOL=1.0D-8
       ! Default: VOFTOL=1.0D-10 (prototype code)
      REAL_T, PARAMETER :: VOFTOL=1.0D-8
      REAL_T, PARAMETER :: VOFTOL_AREAFRAC=1.0D-1
       ! Default: VOFTOL_MULTI_VOLUME=1.0D-12
      REAL_T, PARAMETER :: VOFTOL_MULTI_VOLUME=1.0D-12
      REAL_T, PARAMETER :: VOFTOL_MULTI_VOLUME_SANITY=1.0D-8
       ! used for checking for 0 magnitude
      REAL_T, PARAMETER :: MLSVOFTOL=1.0D-14
       ! used for checking if centroid in box.
      REAL_T, PARAMETER :: CENTOL=1.0D-13
       ! Default: INTERCEPT_TOL=1.0D-12
      REAL_T, PARAMETER :: INTERCEPT_TOL=1.0D-12
      REAL_T, PARAMETER :: FRAC_PAIR_TOL=1.0D-12

      REAL_T, PARAMETER :: TANGENT_EPS=1.0D-2

      INTEGER_T, PARAMETER :: SEM_IMAGE_BC_ALG=1

      INTEGER_T, PARAMETER :: POLYGON_LIST_MAX=1000

       ! variables for comparing with Fred Stern's experiment
      INTEGER_T, PARAMETER :: SHALLOW_M=100
      INTEGER_T, PARAMETER :: SHALLOW_N=1000
      REAL_T, PARAMETER :: SHALLOW_TIME=12.0

      REAL_T shallow_water_data(0:SHALLOW_M,0:SHALLOW_N,2)
      REAL_T inflow_time(3000)
      REAL_T inflow_elevation(3000)
      REAL_T inflow_velocity(3000)
      REAL_T outflow_time(3000)
      REAL_T outflow_elevation(3000)
      REAL_T outflow_velocity(3000)
      INTEGER_T inflow_count,outflow_count
      INTEGER_T last_inflow_index,last_outflow_index

      INTEGER_T, PARAMETER :: recalesce_num_state=6
      INTEGER_T recalesce_material(100)
      REAL_T recalesce_state_old(recalesce_num_state*100)

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

      INTEGER_T :: used_probtypes(1000)
      INTEGER_T :: probtype_list_size

      ABSTRACT INTERFACE

      subroutine TEMPLATE_INIT_MODULE()
      end subroutine TEMPLATE_INIT_MODULE

      subroutine TEMPLATE_hydro_pressure_density( &
                        xpos,rho,pres)
      REAL_T, intent(in) :: xpos(SDIM)
      REAL_T, intent(out) :: rho
      REAL_T, intent(out) :: pres
      end subroutine TEMPLATE_hydro_pressure_density

      subroutine TEMPLATE_CFL_HELPER(time,dir,uu,dx)
      INTEGER_T, intent(in) :: dir
      REAL_T, intent(in) :: time
      REAL_T, intent(inout) :: uu
      REAL_T, intent(in) :: dx(SDIM)
      end subroutine TEMPLATE_CFL_HELPER


      subroutine TEMPLATE_SUMINT(GRID_DATA_IN,increment_out,nsum)
      use probcommon_module_types

      INTEGER_T, intent(in) :: nsum
      type(user_defined_sum_int_type), intent(in) :: GRID_DATA_IN
      REAL_T, intent(out) :: increment_out(nsum)
      end subroutine TEMPLATE_SUMINT
      
      subroutine TEMPLATE_LS(x,t,LS,nmat)
      INTEGER_T, intent(in) :: nmat
      REAL_T, intent(in) :: x(SDIM)
      REAL_T, intent(in) :: t
      REAL_T, intent(out) :: LS(nmat)
      end subroutine TEMPLATE_LS

      subroutine TEMPLATE_VEL(x,t,LS,VEL,velsolid_flag,dx,nmat)
      INTEGER_T, intent(in) :: nmat
      REAL_T, intent(in) :: x(SDIM)
      REAL_T, intent(in) :: t
      REAL_T, intent(in) :: dx(SDIM)
      REAL_T, intent(in) :: LS(nmat)
      REAL_T, intent(out) :: VEL(SDIM)
      INTEGER_T, intent(in) :: velsolid_flag
      end subroutine TEMPLATE_VEL

      subroutine TEMPLATE_EOS(rho,internal_energy,pressure, &
        imattype,im)
      INTEGER_T, intent(in) :: imattype,im
      REAL_T, intent(in) :: rho
      REAL_T, intent(in) :: internal_energy
      REAL_T, intent(out) :: pressure
      end subroutine TEMPLATE_EOS

      subroutine TEMPLATE_SOUNDSQR(rho,internal_energy,soundsqr, &
        imattype,im)
      INTEGER_T, intent(in) :: imattype,im
      REAL_T, intent(in) :: rho
      REAL_T, intent(in) :: internal_energy
      REAL_T, intent(out) :: soundsqr
      end subroutine TEMPLATE_SOUNDSQR

      subroutine TEMPLATE_INTERNAL(rho,temperature,local_internal_energy, &
        imattype,im)
      INTEGER_T, intent(in) :: imattype,im
      REAL_T, intent(in) :: rho
      REAL_T, intent(in) :: temperature 
      REAL_T, intent(out) :: local_internal_energy
      end subroutine TEMPLATE_INTERNAL

      subroutine TEMPLATE_TEMPERATURE(rho,temperature,internal_energy, &
        imattype,im)
      INTEGER_T, intent(in) :: imattype,im
      REAL_T, intent(in) :: rho
      REAL_T, intent(out) :: temperature 
      REAL_T, intent(in) :: internal_energy
      end subroutine TEMPLATE_TEMPERATURE

      subroutine TEMPLATE_PRES(x,t,LS,PRES,nmat)
      INTEGER_T, intent(in) :: nmat
      REAL_T, intent(in) :: x(SDIM)
      REAL_T, intent(in) :: t
      REAL_T, intent(in) :: LS(nmat)
      REAL_T, intent(out) :: PRES
      end subroutine TEMPLATE_PRES

      subroutine TEMPLATE_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
      INTEGER_T, intent(in) :: bcflag
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: nstate_mat
      REAL_T, intent(in) :: x(SDIM)
      REAL_T, intent(in) :: t
      REAL_T, intent(in) :: LS(nmat)
      REAL_T, intent(out) :: STATE(nmat*nstate_mat)
      end subroutine TEMPLATE_STATE

      subroutine TEMPLATE_LS_BC(xwall,xghost,t,LS, &
       LS_in,dir,side,dx,nmat)
      INTEGER_T, intent(in) :: nmat
      REAL_T, intent(in) :: xwall
      REAL_T, intent(in) :: xghost(SDIM)
      REAL_T, intent(in) :: t
      REAL_T, intent(inout) :: LS(nmat)
      REAL_T, intent(in) :: LS_in(nmat)
      INTEGER_T, intent(in) :: dir,side
      REAL_T, intent(in) :: dx(SDIM)
      end subroutine TEMPLATE_LS_BC

      subroutine TEMPLATE_VEL_BC(xwall,xghost,t,LS, &
        VEL,VEL_in,veldir,dir,side,dx,nmat)
      INTEGER_T, intent(in) :: nmat
      REAL_T, intent(in) :: xwall
      REAL_T, intent(in) :: xghost(SDIM)
      REAL_T, intent(in) :: t
      REAL_T, intent(in) :: LS(nmat)
      REAL_T, intent(inout) :: VEL
      REAL_T, intent(in) :: VEL_in
      INTEGER_T, intent(in) :: veldir,dir,side
      REAL_T, intent(in) :: dx(SDIM)
      end subroutine TEMPLATE_VEL_BC

      subroutine TEMPLATE_PRES_BC(xwall,xghost,t,LS, &
        PRES,PRES_in,dir,side,dx,nmat)
      INTEGER_T, intent(in) :: nmat
      REAL_T, intent(in) :: xwall
      REAL_T, intent(in) :: xghost(SDIM)
      REAL_T, intent(in) :: t
      REAL_T, intent(in) :: LS(nmat)
      REAL_T, intent(inout) :: PRES
      REAL_T, intent(in) :: PRES_in
      INTEGER_T, intent(in) :: dir,side
      REAL_T, intent(in) :: dx(SDIM)
      end subroutine TEMPLATE_PRES_BC

      subroutine TEMPLATE_STATE_BC(xwall,xghost,t,LS, &
       STATE,STATE_merge,STATE_in,im,istate,dir,side,dx, &
       nmat)
      INTEGER_T, intent(in) :: nmat
      REAL_T, intent(in) :: xwall
      REAL_T, intent(in) :: xghost(SDIM)
      REAL_T, intent(in) :: t
      REAL_T, intent(in) :: LS(nmat)
      REAL_T, intent(inout) :: STATE
      REAL_T, intent(inout) :: STATE_merge
      REAL_T, intent(in) :: STATE_in
      INTEGER_T, intent(in) :: dir,side
      REAL_T, intent(in) :: dx(SDIM)
      INTEGER_T, intent(in) :: istate,im
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
      INTEGER_T, intent(in) :: nmat
      INTEGER_T, intent(in) :: im
      REAL_T, intent(in) :: VFRAC(nmat)
      REAL_T, intent(in) :: time
      INTEGER_T, intent(in) :: nhalf
      REAL_T, intent(in) :: x(SDIM)
      REAL_T, intent(in) :: xsten(-nhalf:nhalf,SDIM)
      REAL_T, intent(in) :: temp(nmat)
      REAL_T, intent(in) :: den(nmat)
      REAL_T, intent(in) :: CV(nmat)
      REAL_T, intent(in) :: dt
      REAL_T, intent(out) :: heat_source
      end subroutine TEMPLATE_HEATSOURCE

      subroutine TEMPLATE_EB_heat_source(time,dt,xsten,nhalf, &
        heat_flux,heat_dir,heat_side)
      INTEGER_T, intent(in) :: nhalf
      REAL_T, dimension(-nhalf:nhalf,SDIM), intent(in) :: xsten
      REAL_T, intent(in) :: time
      REAL_T, intent(in) :: dt
      REAL_T, intent(out) :: heat_flux
      INTEGER_T, intent(out) :: heat_dir
      INTEGER_T, intent(out) :: heat_side
      end subroutine TEMPLATE_EB_heat_source

      END INTERFACE

      PROCEDURE(TEMPLATE_INIT_MODULE), POINTER :: SUB_INIT_MODULE
      PROCEDURE(TEMPLATE_hydro_pressure_density), POINTER :: &
              SUB_hydro_pressure_density
      PROCEDURE(TEMPLATE_CFL_HELPER), POINTER :: SUB_CFL_HELPER
      PROCEDURE(TEMPLATE_SUMINT), POINTER :: SUB_SUMINT
      PROCEDURE(TEMPLATE_LS), POINTER :: SUB_LS
      PROCEDURE(TEMPLATE_VEL), POINTER :: SUB_VEL
      PROCEDURE(TEMPLATE_EOS), POINTER :: SUB_EOS
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

contains

end module probcommon_module

