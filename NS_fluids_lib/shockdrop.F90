#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"
#include "AMReX_ArrayLim.H"
#include "EXTRAP_COMP.H"

! MARCO USES inputs3d.nwave probtype=395
#define shockdrop_PROB_TYPE 3001

#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else
print *,"dimension bust"
stop
#endif

MODULE shockdrop
use amrex_fort_module, only : amrex_real
implicit none 

 !see run2d/NormalShockWaveNASA.F90
 !see run2d/inputs.shockdrop

real(amrex_real) shockdrop_R
real(amrex_real) shockdrop_cp
real(amrex_real) shockdrop_cv

real(amrex_real) shockdrop_We
real(amrex_real) shockdrop_Oh
real(amrex_real) shockdrop_Re
real(amrex_real) shockdrop_inertial_time_scale

 !upstream supersonic relative to shock
real(amrex_real) shockdrop_P0
real(amrex_real) shockdrop_T0
real(amrex_real) shockdrop_DEN0
real(amrex_real) shockdrop_M0
real(amrex_real) shockdrop_VEL0
real(amrex_real) shockdrop_init_VEL0
real(amrex_real) shockdrop_C0
real(amrex_real) shockdrop_EE0

 !downstream subsonic relative to shock
real(amrex_real) shockdrop_P1
real(amrex_real) shockdrop_T1
real(amrex_real) shockdrop_DEN1
real(amrex_real) shockdrop_M1
real(amrex_real) shockdrop_VEL1
real(amrex_real) shockdrop_C1
real(amrex_real) shockdrop_EE1

real(amrex_real) shockdrop_gamma

integer, parameter :: n_data=94
real(amrex_real) :: vel_data(n_data) ! m/s
real(amrex_real) :: time_data(n_data) !ms

CONTAINS

!stateIndex=1 U
!stateIndex=2 V
!stateIndex=3 rho
!stateIndex=4 e
!stateIndex=5 T
!stateIndex=6 p
subroutine N_wave_solution(time,x,y,z,T_o,Mn1,shockAngle,stateIndex,preState,postState,dState_dt)
use probcommon_module
use global_utility_module
IMPLICIT NONE

      real(amrex_real), intent(in) :: time
      real(amrex_real), intent(in) :: x,y,z
      real(amrex_real) deltaP
      real(amrex_real) P1,T1,Gamma_constant,GammaP1,c1,e1,rho1
      real(amrex_real) M1
      real(amrex_real) Db,Rb,Lb,Fw
      real(amrex_real) U2,V2,rho2,e2,P2,T2
      real(amrex_real) vn,vt
      real(amrex_real), intent(out) :: shockAngle
      real(amrex_real), intent(out) :: T_o
      real(amrex_real), intent(out) :: Mn1
      real(amrex_real), intent(out) :: preState,postState,dState_dt
      integer, intent(in) :: stateIndex

      real(amrex_real) factor1

      real(amrex_real) dRho_dt,dP_dt,dE_dt,dU_dt,dV_dt,dT_dt

      !Bullet dimensions
      Db = 0.42d0 * 2.8d0 !Diameter
      Rb = Db/2.0d0 !Radius
      Lb = 1.128d0 / 2.0d0 !ogival length of bullet

      !Shape factor
      Fw = Rb / Lb ** 0.25

      !Farfield Mach
      M1 = radblob2

      !Pre-shock conditions
      P1 = inflow_pressure
      T1 = fort_tempconst(2)
      Gamma_constant = 1.4d0
      GammaP1 = Gamma_constant + 1.0d0
      rho1=fort_denconst(2)
      call INTERNAL_air(rho1,T1,e1)
      call SOUNDSQR_air(rho1,e1,c1)
      c1=sqrt(c1)

      deltaP = P1 * 2.0**(0.25) * Gamma_constant*GammaP1**(-0.5)*(M1**2 - 1.)**(1./8.) * Fw * (0.31)**(-.75)

      P2 = P1 + deltaP

      Mn1 = (((P2/P1) - 1.d0) * GammaP1/(2.0d0*Gamma_constant) + 1.0d0) ** 0.5
      shockAngle = asin(Mn1/M1)

       ! ANALYTICAL STUDY OF SONIC BOOM FROM SUPERSONIC PROJECTILES
       ! Gottlieb and Ritzel
       ! Marco Arienti's formulation used a more arcane formulation from
       ! this paper:
       ! Ballistic Wave from Projectiles and Vehicles of Simple Geometry
       ! Varnier and Le Pape
       ! It is more intuitive to have:
       ! T_o= 2^(1/4) ....
       ! and
       ! Fw=Diameter_b/Lb ** 0.25
      T_o = 2.0**(5.0/4.0) * GammaP1**0.5 / c1 * M1 * (M1 ** 2 - 1.0)**(-3.0/8.0) * Fw * (0.31)**0.25
      T_o = T_o / zblob4

      rho1 = fort_denconst(2)
      T1 = fort_tempconst(2)
      call INTERNAL_air(rho1,T1,e1)
      call postshock_air(Mn1,rho1,T1,rho2,T2,P2)
      call INTERNAL_air(rho2,T2,e2)
      call vnpostshock_air(T1,Mn1,vn)
      vt=c1*radblob2*cos(shockAngle)
      V2 = 0.0d0
      U2 = Mn1 * c1 - vn
      e2=e2+0.5d0*(U2*U2+V2*V2)

      factor1 = radblob4

      dP_dt = 2.0d0 * (P1 - P2)  / T_o  * factor1
      dRho_dt = 2.0d0 * (rho1 - rho2) / T_o  * factor1
      dE_dt = 2.0d0 * (e1 - e2) / T_o  * factor1
      dU_dt = 2.0d0 * (0.0d0 - U2) / T_o  * factor1
      dV_dt = 2.0d0 * (0.0d0 - V2) / T_o  * factor1
      dT_dt = 2.0d0 * (T1 - T2) / T_o  * factor1


      if(stateIndex.eq.1)then
         preState = 0.0d0
         postState = U2 + 0.0d0
         dState_dt = dU_dt
      else if(stateIndex.eq.2)then
         preState = 0.0d0
         postState = V2 + 0.0d0
         dState_dt = dV_dt
      else if(stateIndex.eq.3)then
         preState = rho1
         postState = rho2
         dState_dt = dRho_dt
      else if(stateIndex.eq.4)then
         preState = e1
         postState = e2
         dState_dt = dE_dt
      else if(stateIndex.eq.5)then
         preState = T1
         postState = T2
         dState_dt = dT_dt
      else if(stateIndex.eq.6)then
         preState = inflow_pressure
         postState = p2
         dState_dt = dP_dt
      else
         print*,'invalided stateIndex n wave'
         stop
      end if

return
end subroutine N_wave_solution


!This routine is for PROB 394 -- don't call outside without care
subroutine N_wave_solution1(time,x,y,z,shockAngle,Mn1)
use probcommon_module
use global_utility_module
IMPLICIT NONE

      real(amrex_real), intent(in) :: time
      real(amrex_real), intent(in) :: x,y,z
      real(amrex_real) deltaP
      real(amrex_real), intent(out) :: shockAngle
      real(amrex_real), intent(out) :: Mn1

      real(amrex_real) P1,T1,Gamma_constant,GammaP1,c1
      real(amrex_real) M1
      real(amrex_real) Db,Lb,Rb,Fw

      real(amrex_real) rho1,e1
      real(amrex_real) P2

      !Bullet dimensions
      Db = 0.42d0 * 2.8d0 !Diameter
      Rb = Db/2.0d0 !Radius
      Lb = 1.128d0 / 2.0d0 !ogival length of bullet

      !Shape factor
      Fw = Rb / Lb ** 0.25

      !Farfield Mach
      M1 = radblob2

      !Pre-shock conditions
      P1 = inflow_pressure

      T1 = fort_tempconst(2)
      rho1 = fort_denconst(2)
      call INTERNAL_air(rho1,T1,e1)
      Gamma_constant = 1.4d0
      GammaP1 = Gamma_constant + 1.0d0
      call SOUNDSQR_air(rho1,e1,c1)
      c1=sqrt(c1)

      deltaP = P1 * 2.0**(0.25) * Gamma_constant*GammaP1**(-0.5)*(M1**2 - 1.)**(1./8.) * Fw * (0.31)**(-.75)

      P2 = P1 + deltaP

      Mn1 = (((P2/P1) - 1.d0) * GammaP1/(2.0d0*Gamma_constant) + 1.0d0) ** 0.5
      shockAngle = asin(Mn1/M1)

return
end subroutine N_wave_solution1

subroutine N_wave_solution2(time,x,y,z,T_o,Mn1,shockAngle,stateIndex,preState,postState,dState_dt)
use probcommon_module
use global_utility_module
IMPLICIT NONE
      
      real(amrex_real), intent(in) :: time    
      real(amrex_real), intent(in) :: x,y,z
      real(amrex_real) deltaP
      real(amrex_real) P1,T1,Gamma_constant,GammaP1,c1,e1,rho1
      real(amrex_real) M1
      real(amrex_real) Db,Rb,Lb,Fw
      real(amrex_real) U2,V2,rho2,e2,P2,T2
      real(amrex_real) vn,vt
      real(amrex_real), intent(out) :: shockAngle
      real(amrex_real), intent(out) :: T_o
      real(amrex_real), intent(out) :: Mn1
      real(amrex_real), intent(out) :: preState,postState,dState_dt
      integer, intent(in) :: stateIndex

      real(amrex_real) factor1
             
      real(amrex_real) dRho_dt,dP_dt,dE_dt,dU_dt,dV_dt,dT_dt

      !Bullet dimensions
      Db = 0.42d0 * 1.8d0 !Diameter
      Rb = Db/2.0d0 !Radius
      Lb = 1.128d0 / 2.0d0 !ogival length of bullet

      !Shape factor
      Fw = Rb / Lb ** 0.25

      !Farfield Mach
      M1 = radblob2

      !Pre-shock conditions
      P1 = inflow_pressure
      T1 = fort_tempconst(2)
      Gamma_constant = 1.4d0
      GammaP1 = Gamma_constant + 1.0d0
      rho1=fort_denconst(2)
      call INTERNAL_air(rho1,T1,e1)
      call SOUNDSQR_air(rho1,e1,c1)
      c1=sqrt(c1)

      deltaP = P1 * 2.0**(0.25) * Gamma_constant*GammaP1**(-0.5)*(M1**2 - 1.)**(1./8.) * Fw * (0.5)**(-.75)

      P2 = P1 + deltaP

      Mn1 = (((P2/P1) - 1.d0) * GammaP1/(2.0d0*Gamma_constant) + 1.0d0) ** 0.5
      shockAngle = asin(Mn1/M1)

      T_o = 2.0**(5.0/4.0) * GammaP1**0.5 / c1 * M1 * (M1 ** 2 - 1.0)**(-3.0/8.0) * Fw * (0.5)**0.25
      T_o = T_o / 4.0d0

      rho1 = fort_denconst(2)
      T1 = fort_tempconst(2)
      call INTERNAL_air(rho1,T1,e1)
      call postshock_air(Mn1,rho1,T1,rho2,T2,P2)
      call INTERNAL_air(rho2,T2,e2)
      call vnpostshock_air(T1,Mn1,vn)
      vt=c1*radblob2*cos(shockAngle)
      U2=c1*radblob2-vn*sin(shockAngle)-vt*cos(shockAngle)
      V2=-vn*cos(shockAngle)+vt*sin(shockAngle)
      e2=e2+0.5d0*(U2*U2+V2*V2)

      factor1 = 1.0d0

      dP_dt = 2.0d0 * (P1 - P2)  / T_o  * factor1
      dRho_dt = 2.0d0 * (rho1 - rho2) / T_o  * factor1
      dE_dt = 2.0d0 * (e1 - e2) / T_o  * factor1
      dU_dt = 2.0d0 * (0.0d0 - U2) / T_o  * factor1
      dV_dt = 2.0d0 * (0.0d0 - V2) / T_o  * factor1
      dT_dt = 2.0d0 * (T1 - T2) / T_o  * factor1


      if(stateIndex.eq.1)then
         preState = 0.0d0
         postState = U2 + 0.0d0
         dState_dt = dU_dt
      else if(stateIndex.eq.2)then
         preState = 0.0d0
         postState = V2 + 0.0d0
         dState_dt = dV_dt
      else if(stateIndex.eq.3)then
         preState = rho1
         postState = rho2
         dState_dt = dRho_dt
      else if(stateIndex.eq.4)then
         preState = e1
         postState = e2
         dState_dt = dE_dt
      else if(stateIndex.eq.5)then
         preState = T1
         postState = T2
         dState_dt = dT_dt
      else if(stateIndex.eq.6)then
         preState = inflow_pressure
         postState = p2
         dState_dt = dP_dt
      else
         print*,'invalided stateIndex n wave'
         stop
      end if


return
end subroutine N_wave_solution2


subroutine recompute_globals(time)
use probcommon_module
use global_utility_module
IMPLICIT NONE
real(amrex_real), intent(in) :: time
real(amrex_real) test_pres
real(amrex_real) water_pressure
integer :: idata

!some remarks:
! in a frame of reference wrt the shock:
! shockspeed=0.0
! M0 C0 = V0  (supersonic)
! M1 C1 = V1  (subsonic)
! in a frame of reference wrt the drop:
! V0_drop = 0
! shockspeed=-V0
! V1_downstream=V1-V0 (opposite side of the drop wrt shock)
 !material_type=5 EOS_air
 !see PROBCOMMON.F90
shockdrop_R=R_AIR_PARMS !ergs/(Kelvin g) (0.287D+7)
shockdrop_cv=CV_AIR_PARMS !ergs/(Kelvin g) (0.72D+7)
shockdrop_cp=shockdrop_cv+shockdrop_R !ergs/(Kelvin g)

! shockdrop_M0=1.4
! shockdrop_M0=3.0
! shockdrop_M0=1.17
! shockdrop_M0=1.0017
if ((axis_dir.eq.150).or. &
    (axis_dir.eq.151).or. &
    (axis_dir.eq.152).or. &
    (axis_dir.eq.154)) then
 shockdrop_M0=vinletgas
else if (axis_dir.eq.153) then
 shockdrop_M0=radblob2 !far field Mach
else
 print *,"axis_dir invalid: ",axis_dir
 stop
endif

if (shockdrop_M0.gt.one) then
 !do nothing
else
 print *,"shockdrop_M0 invalid: ",shockdrop_M0
 stop
endif

if ((axis_dir.eq.150).or. &
    (axis_dir.eq.151).or. &
    (axis_dir.eq.152).or. &
    (axis_dir.eq.154)) then
 shockdrop_T0=278.0d0 !tempconst(2)
 shockdrop_DEN0=0.00125335272d0 !denconst(2)
else if (axis_dir.eq.153) then
 shockdrop_T0=300.0d0 !tempconst(2)
 shockdrop_DEN0=0.00102 !denconst(2)
else
 print *,"axis_dir invalid: ",axis_dir
 stop
endif

if (abs(shockdrop_DEN0-fort_denconst(2)).le.EPS6) then
 !do nothing
else
 print *,"shockdrop_DEN0 ",shockdrop_DEN0
 print *,"fort_denconst(2) ",fort_denconst(2)
 print *,"mismatch"
 stop
endif

if (abs(shockdrop_T0-fort_tempconst(2)).le.EPS6) then
 !do nothing
else
 print *,"shockdrop_T0 ",shockdrop_T0
 print *,"fort_tempconst(2) ",fort_tempconst(2)
 print *,"mismatch"
 stop
endif
if (abs(shockdrop_T0-fort_tempconst(1)).le.EPS6) then
 !do nothing
else
 print *,"shockdrop_T0 ",shockdrop_T0
 print *,"fort_tempconst(1) ",fort_tempconst(1)
 print *,"mismatch"
 stop
endif


 ! this value must be consistent with material_type=5 parameters
shockdrop_gamma=shockdrop_cp/shockdrop_cv

 ! if compressible liquid, this value should be the same
 ! as the equilibrium pressure of the liquid drop.
shockdrop_P0=(shockdrop_gamma-1.0d0)* &
 shockdrop_DEN0*shockdrop_cv*shockdrop_T0

shockdrop_C0=(shockdrop_gamma*shockdrop_P0/shockdrop_DEN0)**half
shockdrop_VEL0=shockdrop_M0*shockdrop_C0

shockdrop_P1=shockdrop_P0* &
 (two*shockdrop_gamma*(shockdrop_M0**2)-shockdrop_gamma+one)/ &
 (shockdrop_gamma+one)

shockdrop_T1=shockdrop_T0* &
 ( two*shockdrop_gamma*(shockdrop_M0**2)-shockdrop_gamma+one )* &
 ( (shockdrop_gamma-one)*(shockdrop_M0**2)+two )/ &
 ( ((shockdrop_gamma+one)*shockdrop_M0)**2 )

shockdrop_DEN1=shockdrop_DEN0* &
 ((shockdrop_gamma+one)*(shockdrop_M0**2))/ &
 ((shockdrop_gamma-one)*(shockdrop_M0**2)+two)

shockdrop_C1=sqrt(shockdrop_gamma*shockdrop_P1/shockdrop_DEN1)
shockdrop_M1= &
 ((shockdrop_gamma-one)*(shockdrop_M0**2)+two)/ &
 (two*shockdrop_gamma*(shockdrop_M0**2)-shockdrop_gamma+one)
shockdrop_M1=sqrt(shockdrop_M1)
shockdrop_VEL1=shockdrop_M1*shockdrop_C1

!shock_vel0-shock_vel1=M0 * C0 - M1 * C1

 !general_hydrostic_pressure is declared in GLOBALUTIL.F90
call general_hydrostatic_pressure(test_pres)

if (probtype.eq.shockdrop_PROB_TYPE) then

 if ((axis_dir.eq.150).or. &
     (axis_dir.eq.151).or. &
     (axis_dir.eq.152).or. &
     (axis_dir.eq.154)) then
  if (abs(test_pres-shockdrop_P0)/test_pres.gt.EPS3) then
   print *,"shockdrop_P0 inconsistent w/general_hydrostatic_pressure"
   print *,"test_pres=",test_pres
   print *,"shockdrop_P0=",shockdrop_P0
   print *,"shockdrop_gamma=",shockdrop_gamma
   print *,"shockdrop_cv=",shockdrop_cv
   print *,"shockdrop_DEN0=",shockdrop_DEN0
   print *,"shockdrop_T0=",shockdrop_T0
   print *,"fort_denconst(1)=",fort_denconst(1)
   print *,"fort_tempconst(1)=",fort_tempconst(1)
   call EOS_tait_rho(fort_denconst(1),fort_tempconst(1),water_pressure)
   print *,"pressure from EOS_tait_rho: ",water_pressure
   print *,"gravity_vector=",gravity_vector
   stop
  endif
 else if (axis_dir.eq.153) then
  !check nothing
 else
  print *,"axis_dir invalid: ",axis_dir
  stop
 endif

else
 print *,"probtype invalid:",probtype
 stop
endif

if (fort_material_type(2).ne.5) then
 print *,"only material_type=5 supported for gas for this problem"
 print *,"fort_material_type(2)=",fort_material_type(2)
 stop
endif

! in shock frame of reference, the upstream velocity is -|V0| and the
! downstream velocity is -|V1|.  In upstream frame of reference the 
! downstream velocity is |V0|-|V1|

if (shockdrop_VEL0.gt.shockdrop_VEL1) then
 shockdrop_We=shockdrop_DEN1* &
   ((shockdrop_VEL0-shockdrop_VEL1)**2)*two*radblob/ &
   fort_tension(1)
 shockdrop_Oh=fort_viscconst(1)/ &
   sqrt(fort_denconst(1)*fort_tension(1)*two*radblob)
 shockdrop_Re=shockdrop_DEN1*(shockdrop_VEL0-shockdrop_VEL1)* &
   two*radblob/fort_viscconst(2)
 shockdrop_inertial_time_scale=two*radblob* &
    sqrt(fort_denconst(1)/shockdrop_DEN1)/ &
    (shockdrop_VEL0-shockdrop_VEL1) 
else
 print *,"shockdrop_VEL0 or shockdrop_VEL1 invalid: ", &
   shockdrop_VEL0,shockdrop_VEL1
 stop
endif

if (axis_dir.eq.150) then
 !do nothing (shock drop)
else if (axis_dir.eq.151) then
 !do nothing (shock column)
else if (axis_dir.eq.152) then
 !shock cylinder
 idata=1
 do while ((time_data(idata)*1.0D-3.lt.time).and.(idata.lt.n_data)) 
  idata=idata+1
  if ((time_data(idata).ge.time_data(idata-1)).and. &
      (time_data(idata-1).ge.zero)) then
   !do nothing
  else
   print *,"time_data invalid"
   stop
  endif
 enddo
 if (vel_data(idata).ge.zero) then
  !do nothing
 else
  print *,"vel_data invalid"
  stop
 endif
 shockdrop_VEL0=shockdrop_VEL1+vel_data(idata)*1.0D+2
else if (axis_dir.eq.154) then
 !do nothing (Arienti shock cylinder but not compare with experiments)
else if (axis_dir.eq.153) then
 !do nothing (Arienti shock sphere)
else
 print *,"axis_dir invalid: ",axis_dir
 stop
endif

return
end subroutine recompute_globals

! units are cgs
subroutine shockdrop_init()
use probcommon_module
use global_utility_module
IMPLICIT NONE
real(amrex_real) time
integer :: test_n_data
integer :: idata

time=0.0d0

if (axis_dir.eq.150) then
 !do nothing (shock drop)
else if (axis_dir.eq.151) then
 !do nothing (shock column)
else if (axis_dir.eq.152) then
 !shock cylinder
 print *,"reading GuildenbecherDataRaw"
 open(unit=2, file='GuildenbecherDataRaw')
 read(2,*) test_n_data
 print *,"test_n_data=",test_n_data
 if (test_n_data.eq.n_data) then
  !do nothing
 else
  print *,"test_n_data invalid: ",test_n_data,n_data
  stop
 endif
  !time_data ms
  !vel_data  m/s
 do idata=1,n_data
  read(2,*) time_data(idata),vel_data(idata)
 enddo
 close(2)
else if (axis_dir.eq.154) then
 !do nothing
else if (axis_dir.eq.153) then
 !do nothing (shock sphere ARIENTI)
else
 print *,"axis_dir invalid: ",axis_dir
 stop
endif

shockdrop_init_VEL0=shockdrop_VEL0

call recompute_globals(time) !modifies shockdrop_VEL0

print *,"shockdrop: init VEL0 ",shockdrop_init_VEL0
print *,"shockdrop: VEL0 ",shockdrop_VEL0

print *,"shockdrop: upstream den,approaching SPEED,T,M,C ", &
 shockdrop_DEN0,shockdrop_VEL0,shockdrop_T0,shockdrop_M0,shockdrop_C0
print *,"shockdrop: downstream den,SPEED,T,M,C ", &
 shockdrop_DEN1,shockdrop_VEL1,shockdrop_T1,shockdrop_M1,shockdrop_C1

! in shock frame of reference, the upstream velocity is -|V0| and the
! downstream velocity is -|V1|.  In upstream frame of reference the 
! downstream velocity is |V0|-|V1|

if (shockdrop_VEL0.gt.shockdrop_VEL1) then
 print *,"shockdrop_DownStreamVelocity=", &
    shockdrop_VEL0-shockdrop_VEL1
 print *,"shockdrop_We=",shockdrop_We
 print *,"shockdrop_Oh=",shockdrop_Oh
 print *,"shockdrop_Re=",shockdrop_Re
 print *,"shockdrop_inertial_time_scale=",shockdrop_inertial_time_scale
else
 print *,"shockdrop_VEL0 or shockdrop_VEL1 invalid: ", &
   shockdrop_VEL0,shockdrop_VEL1
 stop
endif

return
end subroutine shockdrop_init

subroutine setup_stream(dir_stream,dir_transverse)
use probcommon_module

integer, intent(out) :: dir_stream,dir_transverse

 if (SDIM.eq.3) then
  dir_stream=1
  dir_transverse=2
 else if (SDIM.eq.2) then
  if ((axis_dir.eq.150).or. &
      (axis_dir.eq.151).or. &
      (axis_dir.eq.152).or. &
      (axis_dir.eq.154)) then
   dir_stream=1
   dir_transverse=2
  else if (axis_dir.eq.153) then !Arienti shock sphere
   dir_stream=2
   dir_transverse=1
  else
   print *,"axis_dir invalid: ",axis_dir
   stop
  endif
 else
  print *,"dimension bust: ",SDIM
  stop
 endif

return
end subroutine setup_stream

subroutine shockdrop_velocity(x,t,ls,vel, &
  velsolid_flag,dx,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: dx(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: ls(nmat)
real(amrex_real) :: ls_local
real(amrex_real), INTENT(out) :: vel(SDIM)
integer dir
integer, INTENT(in) :: velsolid_flag
real(amrex_real) :: To,Mn,shockAngle,preState,x_vel,dudt,alpha,vn
real(amrex_real) :: rho1,e1,c1sqr,c1,u1,u2
real(amrex_real) :: Tcross,Tcross2,Tcross3
integer :: dir_stream,dir_transverse

call setup_stream(dir_stream,dir_transverse)

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
if ((velsolid_flag.eq.0).or. &
    (velsolid_flag.eq.1)) then
 ! do nothing
else 
 print *,"velsolid_flag invalid: ",velsolid_flag
 stop
endif
do dir=1,SDIM
 if (dx(dir).gt.zero) then
  ! do nothing
 else
  print *,"dx invalid: ",dx(dir)
  stop
 endif
enddo

do dir=1,SDIM
 vel(dir)=zero
enddo

call shockdrop_dropLS(x(dir_stream),x(dir_transverse),x(SDIM),ls_local)

if (ls_local.ge.zero) then
 ! do nothing in y/z direction (drop is upstream from shock and
 ! stationary in the "upstream frame of reference")
 ! shock velocity > 0
 ! shock is approaching with speed: shockdrop_VEL0
 if ((axis_dir.eq.150).or. &
     (axis_dir.eq.151).or. &
     (axis_dir.eq.152).or. &
     (axis_dir.eq.154)) then
  vel(dir_stream)=advbot ! velocity in the "periodic" direction
 else if (axis_dir.eq.153) then !Arienti shock sphere
  !do nothing
 else
  print *,"axis_dir invalid: ",axis_dir
  stop
 endif

else
 call shockdrop_shockLS(x(dir_stream),x(dir_transverse),x(SDIM),ls_local)
 ! in shock frame of reference:
 ! upstream: v=-shockdrop_VEL0
 ! downstream: v=-shockdrop_VEL1
 if ((axis_dir.eq.150).or. &
     (axis_dir.eq.151).or. &
     (axis_dir.eq.152).or. &
     (axis_dir.eq.154)) then

  if (ls_local.ge.zero) then  ! upstream (above the shock)
   vel(SDIM)=zero
  else if (ls_local.le.zero) then
   vel(SDIM)=shockdrop_VEL0-shockdrop_VEL1
  else
   print *,"ls_local invalid: ",ls_local
   stop
  endif

 else if (axis_dir.eq.153) then !Arienti shock sphere

  if (t.eq.zero) then
   call N_wave_solution(t,x(dir_stream),x(dir_transverse),x(SDIM),To,Mn,shockAngle,1,preState,x_vel,dudt)
   alpha=0.5d0*(1d0+0.5d0*datan(200.5*ls_local/dx(dir_stream))/datan(1d0))
   call vnpostshock_air(fort_tempconst(2),Mn,vn)
   rho1=fort_denconst(2)
   call INTERNAL_air(rho1,fort_tempconst(2),e1)
   call SOUNDSQR_air(rho1,e1,c1sqr)
   c1=sqrt(c1sqr)
   vel(dir_stream)=(1d0-alpha)*(Mn*c1 - vn)
   vel(dir_transverse)=zero
   if (SDIM.eq.3) then
    vel(SDIM)=zero
   endif
  else if (t.gt.zero) then

   call N_wave_solution(t,x(dir_stream),x(dir_transverse),x(SDIM),To,Mn,shockAngle,1,u1,u2,dudt)
   rho1=fort_denconst(2)
   call INTERNAL_air(rho1,fort_tempconst(2),e1)
   call SOUNDSQR_air(rho1,e1,c1sqr)
   Tcross = (x(dir_stream) - (xblob2) + 0.5d0*dx(dir_stream))/(sqrt(c1sqr)*Mn)
   Tcross2 = Tcross + To
   Tcross3 = Tcross2 + To * xblob4
      
   if(t.lt.Tcross) then
    ! pre-shock conditions
    vel(dir_stream) = u1
   else if (t.lt.Tcross2)then
    ! post-shock conditions
    vel(dir_stream) = u2 + dudt * (t-Tcross)
   else if (t.lt.Tcross3)then
    ! re-compression
    u2 = u2 + dudt * (Tcross2-Tcross)
    vel(dir_stream) = u2 - dudt * (t-Tcross2)/(2.0*xblob4)*yblob4
   else
    vel(dir_stream) = u1
   endif

  else
   print *,"t out of range: ",t
   stop
  endif

 else
  print *,"axis_dir invalid: ",axis_dir
  stop
 endif

endif 

return
end subroutine shockdrop_velocity

 !gets mapped to SUB_CFL_HELPER
 !SUB_CFL_HELPER is called from check_user_defined_velbc (PROB.F90)
 !check_user_defined_velbc is called from fort_estdt (GODUNOV_3D.F90)
subroutine shockdrop_maxvelocity(time,dir,vel,dx)
use probcommon_module
IMPLICIT NONE
integer, INTENT(in) :: dir
real(amrex_real), INTENT(in) :: time
real(amrex_real), INTENT(inout) ::  vel
real(amrex_real) :: vel_local
real(amrex_real), INTENT(in) :: dx(SDIM)

call recompute_globals(time)

vel_local=abs(abs(shockdrop_VEL0)-abs(shockdrop_VEL1))
vel_local=max(abs(advbot),abs(vel_local))
vel=max(abs(vel),abs(vel_local))

return
end subroutine shockdrop_maxvelocity


subroutine shockdrop_LS(x,t,LS,nmat)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(out) :: LS(nmat)
integer :: dir_stream,dir_transverse

call setup_stream(dir_stream,dir_transverse)

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

call shockdrop_dropLS(x(dir_stream),x(dir_transverse),x(SDIM),LS(1))
LS(2)=-LS(1)

return
end subroutine shockdrop_LS

subroutine shockdrop_pressure(x,t,ls,pres,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: ls(nmat)
real(amrex_real) :: ls_local
real(amrex_real), INTENT(out) :: pres
real(amrex_real) :: To,Mn,shockAngle,p1,p2,dpdt
real(amrex_real) :: rho1,e1,c1sqr
real(amrex_real) :: Tcross,Tcross2,Tcross3
integer :: dir_stream,dir_transverse

call setup_stream(dir_stream,dir_transverse)

if (num_materials.eq.nmat) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

call shockdrop_dropLS(x(dir_stream),x(dir_transverse),x(SDIM),ls_local)
if (ls_local.ge.zero) then ! liquid
 pres=shockdrop_P0
else
 call shockdrop_shockLS(x(dir_stream),x(dir_transverse),x(SDIM),ls_local)

 if ((axis_dir.eq.150).or. &
     (axis_dir.eq.151).or. &
     (axis_dir.eq.152).or. &
     (axis_dir.eq.154)) then

  if (ls_local.ge.zero) then  ! upstream (above the approaching shock)
   pres=shockdrop_P0
  else
   pres=shockdrop_P1
  endif

 else if (axis_dir.eq.153) then !Arienti shock sphere

    !stateIndex=6 p
  call N_wave_solution(t,x(dir_stream),x(dir_transverse),x(SDIM),To,Mn,shockAngle,6,p1,p2,dpdt)
  rho1=fort_denconst(2)
  call INTERNAL_air(rho1,fort_tempconst(2),e1)
  call SOUNDSQR_air(rho1,e1,c1sqr)
  Tcross = (x(dir_stream) - (xblob2))/(sqrt(c1sqr)*Mn)
  Tcross2 = Tcross + To
  Tcross3 = Tcross2 + To * xblob4

  if (t.lt.Tcross) then
   ! pre-shock conditions
   pres = inflow_pressure
  else if (t.lt.Tcross2) then
   ! post-shock conditions
   pres = p2 + dpdt * (t-Tcross)
  else if (t.lt.Tcross3) then
   ! re-compression
   p2  = p2 + dpdt * (Tcross2-Tcross)
   pres = p2 - dpdt * (t-Tcross2) / (2.0*xblob4)*yblob4
  else
   pres = inflow_pressure
  endif

 else
  print *,"axis_dir invalid: ",axis_dir
  stop
 endif
endif 

return
end subroutine shockdrop_pressure


subroutine shockdrop_gas_density(t,x,y,z,den)
use probcommon_module
use global_utility_module
IMPLICIT NONE
real(amrex_real), intent(in) :: t,x,y,z
real(amrex_real), intent(out) :: den
real(amrex_real) LS
real(amrex_real) :: To,Mn,shockAngle,vn,alpha,preState,postState,dState_dt
real(amrex_real) :: rho1,e1,c1sqr,c1,T1,T2,rho2,e2,vt,U2,V2,P2
real(amrex_real) :: Tcross,Tcross2,Tcross3
real(amrex_real) :: dx_local(SDIM)
integer :: dir,ilev
integer :: dir_stream,dir_transverse

call setup_stream(dir_stream,dir_transverse)

if (fort_finest_level.ge.0) then
 !do nothing
else
 print *,"fort_finest_level invalid: ",fort_finest_level
 stop
endif

do dir=1,SDIM
 if (fort_n_cell(dir).gt.0) then
  if (problen_array(dir).gt.zero) then
   dx_local(dir)=problen_array(dir)/fort_n_cell(dir)
   do ilev=1,fort_finest_level
    dx_local(dir)=dx_local(dir)/two
   enddo
   if (dx_local(dir).gt.zero) then
    !do nothing
   else
    print *,"dx_local invalid: ",dir,dx_local(dir)
    stop
   endif
  else
   print *,"problen_array invalid: ",dir,problen_array(dir)
   stop
  endif
 else
  print *,"fort_n_cell invalid: ",dir,fort_n_cell(dir)
  stop
 endif
enddo  !dir=1..sdim

call shockdrop_dropLS(x,y,z,LS)
if (LS.ge.zero) then ! liquid
 den=shockdrop_DEN0 !fort_denconst(2)
else if (LS.le.zero) then !gas
 call shockdrop_shockLS(x,y,z,LS)

 if ((axis_dir.eq.150).or. &
     (axis_dir.eq.151).or. &
     (axis_dir.eq.152).or. &
     (axis_dir.eq.154)) then

  if (LS.ge.zero) then  ! upstream (above the shock)
   den=shockdrop_DEN0
  else
   den=shockdrop_DEN1
  endif

 else if (axis_dir.eq.153) then !Arienti shock sphere

  if (t.eq.zero) then
   call N_wave_solution1(t,x,y,z,shockAngle,Mn)
   rho1=fort_denconst(2)
   T1=fort_tempconst(2)
   call INTERNAL_air(rho1,T1,e1)
   call SOUNDSQR_air(rho1,e1,c1sqr)
   c1=sqrt(c1sqr)
   call postshock_air(Mn,rho1,T1,rho2,T2,P2)
   call INTERNAL_air(rho2,T2,e2)
   call vnpostshock_air(T1,Mn,vn)
   vt=c1*radblob2*cos(shockAngle)
   U2= Mn * c1  - vn
   V2= 0.0d0
   alpha=0.5d0*(1d0+0.5d0*datan(200.5*LS/dx_local(dir_stream))/datan(1d0))
   den=rho1*alpha+rho2*(1d0-alpha)

  else if (t.gt.zero) then

   call N_wave_solution(t,x,y,z,To,Mn,shockAngle,3, &
    preState,postState,dState_dt)
   rho1=fort_denconst(2)
   call INTERNAL_air(rho1,fort_tempconst(2),e1)
   call SOUNDSQR_air(rho1,e1,c1sqr)
   Tcross = (x - (xblob2))/(sqrt(c1sqr)*Mn)
   Tcross2 = Tcross + To
   Tcross3 = Tcross2 + To * xblob4
   if(t.lt.Tcross) then
    ! pre-shock conditions
    den = preState
   else if (t.lt.Tcross2)then
    ! post-shock conditions
    den = postState + dState_dt * (t-Tcross)
   else if (t.lt.Tcross3)then
    ! re-compression
    den = postState + dState_dt * (Tcross2-Tcross)
    den = den - dState_dt * (t-Tcross2) / (2.0*xblob4)*yblob4
   else
    den = preState
   endif

  else
   print *,"t invalid: ",t
   stop
  endif

 else
  print *,"axis_dir invalid: ",axis_dir
  stop
 endif

else
 print *,"LS invalid: ",LS
 stop
endif 

return
end subroutine shockdrop_gas_density


subroutine shockdrop_gas_temperature(t,x,y,z,temp)
use probcommon_module
use global_utility_module
IMPLICIT NONE
real(amrex_real), intent(in) :: t,x,y,z
real(amrex_real), intent(out) :: temp
real(amrex_real) LS
real(amrex_real) :: To,Mn,shockAngle,vn,alpha,preState,postState,dState_dt
real(amrex_real) :: rho1,e1,c1sqr,c1,T1,T2,rho2,e2,vt,U2,V2,P2
real(amrex_real) :: Tcross,Tcross2,Tcross3
real(amrex_real) :: den
real(amrex_real) :: dx_local(SDIM)
integer :: dir,ilev
integer :: dir_stream,dir_transverse

call setup_stream(dir_stream,dir_transverse)

if (fort_finest_level.ge.0) then
 !do nothing
else
 print *,"fort_finest_level invalid: ",fort_finest_level
 stop
endif

do dir=1,SDIM
 if (fort_n_cell(dir).gt.0) then
  if (problen_array(dir).gt.zero) then
   dx_local(dir)=problen_array(dir)/fort_n_cell(dir)
   do ilev=1,fort_finest_level
    dx_local(dir)=dx_local(dir)/two
   enddo
   if (dx_local(dir).gt.zero) then
    !do nothing
   else
    print *,"dx_local invalid: ",dir,dx_local(dir)
    stop
   endif
  else
   print *,"problen_array invalid: ",dir,problen_array(dir)
   stop
  endif
 else
  print *,"fort_n_cell invalid: ",dir,fort_n_cell(dir)
  stop
 endif
enddo  !dir=1..sdim

call shockdrop_dropLS(x,y,z,LS)

if (LS.ge.zero) then ! liquid

 temp=shockdrop_T0

else if (LS.le.zero) then !gas

 call shockdrop_shockLS(x,y,z,LS)

 if ((axis_dir.eq.150).or. &
     (axis_dir.eq.151).or. &
     (axis_dir.eq.152).or. &
     (axis_dir.eq.154)) then

  if (LS.ge.zero) then  ! upstream (above the shock)
   temp=shockdrop_T0
  else
   temp=shockdrop_T1
  endif

 else if (axis_dir.eq.153) then !Arienti shock sphere

  if (t.eq.zero) then

   call N_wave_solution1(t,x,y,z,shockAngle,Mn)
   rho1=fort_denconst(2)
   T1=fort_tempconst(2)
   call INTERNAL_air(rho1,T1,e1)
   call SOUNDSQR_air(rho1,e1,c1sqr)
   c1=sqrt(c1sqr)
   call postshock_air(Mn,rho1,T1,rho2,T2,P2)
   call INTERNAL_air(rho2,T2,e2)
   call vnpostshock_air(T1,Mn,vn)
   vt=c1*radblob2*cos(shockAngle)
   U2= Mn * c1  - vn
   V2= 0.0d0
   alpha=0.5d0*(1d0+0.5d0*datan(200.5*LS/dx_local(dir_stream))/datan(1d0))
   den=rho1*alpha+rho2*(1d0-alpha)
   call TEMPERATURE_air(den,temp,e1*alpha+e2*(1d0-alpha))

  else if (t.gt.zero) then

   call N_wave_solution(t,x,y,z,To,Mn,shockAngle,5, &
    preState,postState,dState_dt)
   rho1=fort_denconst(2)
   call INTERNAL_air(rho1,fort_tempconst(2),e1)
   call SOUNDSQR_air(rho1,e1,c1sqr)
   Tcross = (x - (xblob2))/(sqrt(c1sqr)*Mn)
   Tcross2 = Tcross + To
   Tcross3 = Tcross2 + To * xblob4
   if(t.lt.Tcross) then
    ! pre-shock conditions
    temp = preState
   else if (t.lt.Tcross2)then
    ! post-shock conditions
    temp = postState + dState_dt * (t-Tcross)
   else if (t.lt.Tcross3)then
    ! re-compression
    temp = postState + dState_dt * (Tcross2-Tcross)
    temp = temp - dState_dt * (t-Tcross2) / (2.0*xblob4)*yblob4
   else
    temp = preState
   endif

  else
   print *,"t invalid: ",t
   stop
  endif

 else
  print *,"axis_dir invalid: ",axis_dir
  stop
 endif

else
 print *,"LS invalid: ",LS
 stop
endif 

return
end subroutine shockdrop_gas_temperature


! probtype="shockdrop_PROB_TYPE" in the inputs file
! axis_dir=150 shock drop
! axis_dir=151 shock column
! axis_dir=152,154 shock cylinder
! axis_dir=153 Arienti shock sphere
! LS>0 upstream of the shock  z>zblob2
! LS<0 downstream of the shock z<zblob2
SUBROUTINE shockdrop_shockLS(x,y,z,LS)
use probcommon_module
IMPLICIT NONE
real(amrex_real),INTENT(in) :: x,y,z
real(amrex_real),INTENT(out) :: LS

if ((axis_dir.eq.150).or. &
    (axis_dir.eq.151).or. &
    (axis_dir.eq.152).or. &
    (axis_dir.eq.154)) then

! downstream vel>0    shock-->   upstream vel=0
! in shock frame of reference:
! (subsonic)downstream<-   shock vel=0  (supersonic)upstream <--
 LS=z-zblob2
else if (axis_dir.eq.153) then !Arienti shock sphere
 LS=x-xblob2
else
 print *,"axis_dir invalid: ",axis_dir
 stop
endif

return
END SUBROUTINE shockdrop_shockLS

! probtype="shockdrop_PROB_TYPE" in the inputs file
! axis_dir=150 shock drop
! axis_dir=151 shock column
! axis_dir=152,154 shock cylinder
! axis_dir=153 Arienti shock sphere
! LS>0 in the drop
subroutine shockdrop_dropLS(x,y,z,LS)
use probcommon_module
IMPLICIT NONE
real(amrex_real),INTENT(in) :: x,y,z
real(amrex_real),INTENT(out) :: LS
real(amrex_real) mag

if (axis_dir.eq.150) then ! shock drop
 mag=(x-xblob)**2+(y-yblob)**2
 if (SDIM.eq.3) then
  mag=mag+(z-zblob)**2
 endif
 LS=radblob-sqrt(mag) 
else if (axis_dir.eq.151) then ! shock column
 if (SDIM.eq.2) then
  mag=(y-yblob)**2
  LS=radblob-sqrt(mag) 
 else if (SDIM.eq.3) then
  mag=(z-zblob)**2
  LS=radblob-sqrt(mag) 
 else 
  print *,"dimension bust"
  stop
 endif
else if ((axis_dir.eq.152).or. &
         (axis_dir.eq.154)) then ! shock cylinder
 if (SDIM.eq.2) then
  mag=(x-xblob)**2+(y-yblob)**2
  LS=radblob-sqrt(mag) 
 else if (SDIM.eq.3) then
  mag=(y-yblob)**2+(z-zblob)**2
  LS=radblob-sqrt(mag) 
 else 
  print *,"dimension bust"
  stop
 endif
else if (axis_dir.eq.153) then ! Arienti shock sphere
 if (SDIM.eq.2) then
  mag=(x-xblob3)**2+(y-yblob3)**2
  LS=radblob3-sqrt(mag) 
 else if (SDIM.eq.3) then
  mag=(x-xblob3)**2+(y-yblob3)**2+(z-zblob3)**2
  LS=radblob3-sqrt(mag) 
 else 
  print *,"dimension bust: ",SDIM
  stop
 endif
else
 print *,"axis_dir invalid in shockdrop_dropLS: ",axis_dir
 stop
endif

return
END SUBROUTINE shockdrop_dropLS


subroutine shockdrop_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: bcflag !0=called from initialize  1=called from bc
integer, INTENT(in) :: nmat
integer, INTENT(in) :: nstate_mat
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: LS(nmat)
real(amrex_real), INTENT(out) :: STATE(nmat*nstate_mat)
integer im,ibase,n
integer :: dir_stream,dir_transverse

call setup_stream(dir_stream,dir_transverse)

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
if (nstate_mat.eq.num_state_material) then
 ! do nothing
else
 print *,"nstate_mat invalid"
 stop
endif

if ((num_materials.eq.2).and. &
    (num_state_material.ge.2).and. &
    (probtype.eq.shockdrop_PROB_TYPE)) then
 do im=1,num_materials
  ibase=(im-1)*num_state_material
  STATE(ibase+ENUM_DENVAR+1)=fort_denconst(im)

  if (im.eq.1) then
   !do nothing
  else if (im.eq.2) then
   call shockdrop_gas_density(t,x(dir_stream),x(dir_transverse),x(SDIM),STATE(ibase+ENUM_DENVAR+1))
  else
   print *,"im out of range"
   stop
  endif

  if (t.eq.zero) then
   STATE(ibase+ENUM_TEMPERATUREVAR+1)=fort_initial_temperature(im)
  else if (t.gt.zero) then
   STATE(ibase+ENUM_TEMPERATUREVAR+1)=fort_tempconst(im)
  else
   print *,"t invalid"
   stop
  endif

  if (im.eq.1) then
   !do nothing
  else if (im.eq.2) then
   call shockdrop_gas_temperature(t,x(dir_stream),x(dir_transverse),x(SDIM), &
     STATE(ibase+ENUM_TEMPERATUREVAR+1))
  else
   print *,"im out of range: ",im
   stop
  endif

  do n=1,num_species_var
   STATE(ibase+ENUM_SPECIESVAR+n)=fort_speciesconst((n-1)*num_materials+im)
  enddo
 enddo ! im=1..num_materials
else
 print *,"num_materials,num_state_material, or probtype invalid"
 stop
endif
 
return
end subroutine shockdrop_STATE

 ! dir=1..sdim  side=1..2
subroutine shockdrop_LS_BC(xwall,xghost,t,LS, &
   LS_in,dir,side,dx,nmat)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: xwall
real(amrex_real), INTENT(in) :: xghost(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(inout) :: LS(nmat)
real(amrex_real), INTENT(in) :: LS_in(nmat)
integer, INTENT(in) :: dir,side
real(amrex_real), INTENT(in) :: dx(SDIM)

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then
 call shockdrop_LS(xghost,t,LS,nmat)
else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine shockdrop_LS_BC


 ! dir=1..sdim  side=1..2 veldir=1..sdim
subroutine shockdrop_VEL_BC(xwall,xghost,t,LS, &
   VEL,VEL_in,veldir,dir,side,dx,nmat)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: xwall
real(amrex_real), INTENT(in) :: xghost(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: LS(nmat)
real(amrex_real), INTENT(inout) :: VEL
real(amrex_real), INTENT(in) :: VEL_in
integer, INTENT(in) :: veldir,dir,side
real(amrex_real), INTENT(in) :: dx(SDIM)
real(amrex_real) local_VEL(SDIM)
integer velsolid_flag

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
velsolid_flag=0
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2).and. &
    (veldir.ge.1).and.(veldir.le.SDIM)) then

 call shockdrop_velocity(xghost,t,LS,local_VEL,velsolid_flag,dx,nmat)
 VEL=local_VEL(veldir)

else
 print *,"dir,side, or veldir invalid"
 stop
endif

return
end subroutine shockdrop_VEL_BC

 ! dir=1..sdim  side=1..2
subroutine shockdrop_PRES_BC(xwall,xghost,t,LS, &
   PRES,PRES_in,dir,side,dx,nmat)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: xwall
real(amrex_real), INTENT(in) :: xghost(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: LS(nmat)
real(amrex_real), INTENT(inout) :: PRES
real(amrex_real), INTENT(in) :: PRES_in
integer, INTENT(in) :: dir,side
real(amrex_real), INTENT(in) :: dx(SDIM)

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif


if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then

 call shockdrop_pressure(xghost,t,LS,PRES,nmat)

else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine shockdrop_PRES_BC

 ! dir=1..sdim  side=1..2
subroutine shockdrop_STATE_BC(xwall,xghost,t,LS, &
   STATE,STATE_merge,STATE_in,im,istate,dir,side,dx,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: xwall
real(amrex_real), INTENT(in) :: xghost(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: LS(nmat)
real(amrex_real) local_STATE(nmat*num_state_material)
real(amrex_real), INTENT(inout) :: STATE
real(amrex_real), INTENT(inout) :: STATE_merge
real(amrex_real), INTENT(in) :: STATE_in
integer, INTENT(in) :: dir,side
real(amrex_real), INTENT(in) :: dx(SDIM)
integer, INTENT(in) :: istate,im
integer ibase,im_crit
integer local_bcflag

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
local_bcflag=1


if ((istate.ge.1).and. &
    (istate.le.num_state_material).and. &
    (im.ge.1).and. &
    (im.le.num_materials)) then
 call shockdrop_STATE(xghost,t,LS,local_STATE, &
         local_bcflag,nmat,num_state_material)
 ibase=(im-1)*num_state_material
 STATE=local_STATE(ibase+istate)
 call get_primary_material(LS,im_crit)
 ibase=(im_crit-1)*num_state_material
 STATE_merge=local_STATE(ibase+istate)
else
 print *,"istate invalid"
 stop
endif

return
end subroutine shockdrop_STATE_BC

subroutine shockdrop_OVERRIDE_TAGFLAG( &
  i,j,k, &
  level,max_level, &
  snew_ptr,lsnew_ptr, &
  xsten,nhalf,time, &
  rflag,tagflag)
use amrex_fort_module, only : amrex_real
use probcommon_module
use global_utility_module
IMPLICIT NONE
integer, INTENT(in) :: i,j,k
integer, INTENT(in) :: level,max_level
integer, INTENT(in) :: nhalf
real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
real(amrex_real), INTENT(in) :: time
real(amrex_real), INTENT(inout) :: rflag
integer, INTENT(inout) :: tagflag
real(amrex_real), INTENT(in),pointer :: snew_ptr(D_DECL(:,:,:),:)
real(amrex_real), INTENT(in),pointer :: lsnew_ptr(D_DECL(:,:,:),:)
real(amrex_real), dimension(3) :: local_x
real(amrex_real), dimension(SDIM) :: local_delta
real(amrex_real) :: F_LIQUID,LS_LIQUID,LS_SHOCK,P_diff,mag
integer :: dir
integer :: ii,jj,kk
integer :: dir_stream,dir_transverse

call setup_stream(dir_stream,dir_transverse)

if (nhalf.lt.3) then
 print *,"nhalf invalid shock drop override tagflag"
 stop
endif
if ((level.ge.0).and.(level.lt.max_level)) then
 ! do nothing
else
 print *,"level and/or max_level invalid"
 print *,"level=",level
 print *,"max_level=",max_level
 stop
endif
do dir=1,SDIM
 local_x(dir)=xsten(0,dir)
enddo
local_x(3)=xsten(0,SDIM)
do dir=1,SDIM
 local_delta(dir)=xsten(1,dir)-xsten(-1,dir)
 if (local_delta(dir).gt.zero) then
  ! do nothing
 else
  print *,"local_delta invalid shockdrop_override_tagflag"
  stop
 endif
enddo !dir=1..sdim

if ((num_materials.ge.2).and. &
    (probtype.eq.shockdrop_PROB_TYPE)) then

 if ((axis_dir.ge.150).and.(axis_dir.le.154)) then

  ii=0
  jj=0
  kk=0
  if (AMREX_SPACEDIM.eq.2) then
   jj=1
  else if (AMREX_SPACEDIM.eq.3) then
   kk=1
  else
   print *,"dimension problem"
   stop
  endif

  rflag=0.0d0
  tagflag=0
  call shockdrop_shockLS(local_x(dir_stream),local_x(dir_transverse),local_x(SDIM),LS_SHOCK)
  LS_LIQUID=lsnew_ptr(D_DECL(i,j,k),1)
  F_LIQUID=snew_ptr(D_DECL(i,j,k),STATECOMP_MOF+1)
  if (time.eq.zero) then
   if (abs(LS_SHOCK).le.two*local_delta(SDIM)) then
    p_diff=abs((shockdrop_P1-shockdrop_P0)/shockdrop_P0)
   else if (abs(LS_SHOCK).ge.two*local_delta(SDIM)) then
    P_diff=zero
   else
    print *,"LS_SHOCK invalid: ",LS_SHOCK
    stop
   endif
  else if (time.gt.zero) then
   P_diff=abs((snew_ptr(D_DECL(i+ii,j+jj,k+kk),STATECOMP_PRES+1)- &
               snew_ptr(D_DECL(i,j,k),STATECOMP_PRES+1))/shockdrop_P0)
  else
   print *,"time invalid: ",time
   stop
  endif

  if (axis_dir.eq.150) then
   !do nothing (shock drop)
  else if (axis_dir.eq.151) then
   !do nothing (shock column)
  else if (axis_dir.eq.152) then
   !shock cylinder
   P_diff=zero
  else if (axis_dir.eq.154) then
   !shock cylinder
   P_diff=zero
  else if (axis_dir.eq.153) then
   !shock sphere (Arienti)
   P_diff=zero
  else
   print *,"axis_dir invalid: ",axis_dir
   stop
  endif

  mag=(local_x(dir_transverse)-yblob)**2
  if (AMREX_SPACEDIM.eq.3) then
   mag=mag+(local_x(SDIM)-zblob)**2
  endif
  mag=sqrt(mag)

  if ((axis_dir.eq.150).or. &
      (axis_dir.eq.151).or. &
      (axis_dir.eq.152).or. &
      (axis_dir.eq.154)) then

   if (mag.gt.100.d0*radblob) then
    !do nothing
   else if (mag.ge.zero) then
    if ((LS_LIQUID.ge.zero).or.(F_LIQUID.ge.0.1d0)) then
     rflag=1.0d0
     tagflag=1
    endif
    if (P_diff.ge.0.01d0) then
     rflag=1.0d0
     tagflag=1
    endif
   else
    print *,"mag invalid shockdrop.F90: ",mag
    stop
   endif

  else if (axis_dir.eq.153) then

   if ((LS_LIQUID.ge.zero).or.(F_LIQUID.ge.0.1d0)) then
    rflag=1.0d0
    tagflag=1
   endif

  else
   print *,"axis_dir invalid shockdrop.F90: ",axis_dir
   stop
  endif

 else
  print *,"axis_dir invalid shockdrop.F90: ",axis_dir
  stop
 endif

else
 print *,"num_materials or probtype invalid"
 print *,"num_materials: ",num_materials
 print *,"probtype: ",probtype
 stop
endif

end subroutine shockdrop_OVERRIDE_TAGFLAG

END MODULE shockdrop
