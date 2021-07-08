#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"

#include "AMReX_ArrayLim.H"


#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else
print *,"dimension bust"
stop
#endif

! probtype==423 (see run2d/inputs.CRYOGENIC_TANK_MK)
module CRYOGENIC_TANK_MK_module

implicit none

INTEGER_T, PARAMETER :: TANK_MK_MATERIAL_TYPE=24

!! MIDDLE OF THE TANK IS AT Z=0 
! Tank inner radius
REAL_T :: TANK_MK_RADIUS
! Tank inner height (inner wall height)
REAL_T :: TANK_MK_HEIGHT
! Location of liquid-gas interface in respect to
! middle of tank at z=0
REAL_T :: TANK_MK_INTERFACE_LOCATION
! Tank speherical end radius
REAL_T :: TANK_MK_END_RADIUS
! Tank speherical end curvature center (0,C_z)
REAL_T :: TANK_MK_END_CENTER
! Heater flux
REAL_T :: TANK_MK_HEATER_FLUX
! Heater location in dim=2 direction
REAL_T :: TANK_MK_HEATER_LOW
REAL_T :: TANK_MK_HEATER_HIGH
REAL_T :: TANK_MK_HEATER_R
REAL_T :: TANK_MK_HEATER_R_LOW

! Flat or spherical interface
REAL_T :: TANK_MK_INTERFACE_RADIUS
REAL_T :: TANK_MK_BUBBLE_X
REAL_T :: TANK_MK_BUBBLE_Y
REAL_T :: TANK_MK_BUBBLE_Z

! Initial mixture pressure at the hight point of the tank
REAL_T :: TANK_MK_INITIAL_PRESSURE
! Universal gas constant [J/(mol K)]
REAL_T :: TANK_MK_R_UNIV

REAL_T :: TANK_MK_GAS_GAMMA
REAL_T :: TANK_MK_GAS_CP
REAL_T :: TANK_MK_GAS_CV

contains

 ! do any initial preparation needed
 subroutine INIT_CRYOGENIC_TANK_MK_MODULE()
  use probcommon_module
  implicit none
  
  
  
  TANK_MK_RADIUS             = xblob
  TANK_MK_HEIGHT             = yblob
  TANK_MK_INTERFACE_LOCATION = zblob

  TANK_MK_INTERFACE_RADIUS = radblob2
  TANK_MK_BUBBLE_X         = xblob2
  TANK_MK_BUBBLE_Y         = yblob2
  TANK_MK_BUBBLE_Z         = zblob2

  TANK_MK_HEATER_FLUX      = xblob3
  TANK_MK_HEATER_LOW       = yblob3
  TANK_MK_HEATER_HIGH      = zblob3
  TANK_MK_HEATER_R         = radblob3
  TANK_MK_HEATER_R_LOW     = radblob4

  TANK_MK_END_RADIUS       = xblob4
  TANK_MK_END_CENTER       = yblob4

  
  

  ! ASSUMING IDEA GAS => The gas heat cpacities should satisfy this
  ! R_spc = C_{p,spc}-C_{v,spc}
  ! to have ideal mixture gas as well =>
  ! Only C_p or C_v can be picked from table and the other one
  ! calculated from equation above.
  ! Here we pick C_{v,spc} from input file.
  ! C_{p,spc} = C_{v,spc} + R_spc

  TANK_MK_R_UNIV = 8.31446261815324D0

  TANK_MK_GAS_CV = fort_stiffCV(2) ![J/(kg K)]
!  TANK_MK_GAS_CP = fort_stiffCP(2) ![J∕(kg·K)]
  TANK_MK_GAS_CP = TANK_MK_GAS_CV + TANK_MK_R_UNIV/fort_molar_mass(2)  ! [J∕(kg·K)]
  TANK_MK_GAS_GAMMA = TANK_MK_GAS_CP / TANK_MK_GAS_CV

  ! Initial pressure based on the given density and temeprature
  ! at the highest point of the tank
  if(fort_material_type(2).eq.0) then
   ! incompressible
   TANK_MK_INITIAL_PRESSURE = outflow_pressure
  elseif(fort_material_type(2).eq.TANK_MK_MATERIAL_TYPE) then
   ! compressible
   ! P = rho R_sp T = rho (gamma-1) U
   TANK_MK_INITIAL_PRESSURE = &
    fort_denconst(2)*(TANK_MK_GAS_GAMMA-one)*&
    fort_initial_temperature(2)*TANK_MK_GAS_CV
  else
   print *,"material type invalid for pressure setup!"
   stop
  endif
  
  return
 end subroutine INIT_CRYOGENIC_TANK_MK_MODULE


 ! fluids tessellate the domain, solids are immersed. 
 ! fluid interfaces are extended into solids.
 ! material 1 is liquid
 ! material 2 is gas
 ! material 3 is solid

 ! MK = Validation of two-phase CFD models for propellant tank 
 !    self-pressurization: Crossing fluid types, scales, and gravity levels 
 !  Mohammad Kassemi , Olga Kartuzova, Sonya Hylton

 subroutine CRYOGENIC_TANK_MK_LS(x,t,LS,nmat)
  use probcommon_module
  IMPLICIT NONE

  INTEGER_T, intent(in) :: nmat
  REAL_T, intent(in) :: x(SDIM)
  REAL_T, intent(in) :: t
  REAL_T, intent(out) :: LS(nmat)

  if (nmat.eq.num_materials) then
   ! do nothing
  else
   print *,"nmat invalid"
   stop
  endif

   ! material 1= liquid  (e.g. Freon 113)
   ! material 2= vapor  
   ! material 3= tank geometry (e.g. acrylic)
  if ((num_materials.eq.3).and.(probtype.eq.423)) then
   ! liquid
   if (TANK_MK_INTERFACE_RADIUS.eq.zero) then
    LS(1)=TANK_MK_INTERFACE_LOCATION-x(SDIM)
   else if (TANK_MK_INTERFACE_RADIUS.gt.zero) then
    if (SDIM.eq.2) then
     LS(1)=sqrt((x(1)-TANK_MK_BUBBLE_X)**2+&
                (x(2)-TANK_MK_BUBBLE_Y)**2)&
               -TANK_MK_INTERFACE_RADIUS
    else if (SDIM.eq.3) then 
     LS(1)=sqrt((x(1)-TANK_MK_BUBBLE_X)**2+&
                (x(2)-TANK_MK_BUBBLE_Y)**2+&
                (x(SDIM)-TANK_MK_BUBBLE_Z)**2)&
               -TANK_MK_INTERFACE_RADIUS
    else
     print *,"dimension bust"
     stop
    endif
   else
    print *,"radblob2 invalid"
    stop
   endif 
   LS(2)=-LS(1)
   ! Solid
   LS(3)=SOLID_TOP_HALF_DIST(x)
  else
   print *,"num_materials ", num_materials
   print *,"probtype ", probtype
   print *,"num_materials or probtype invalid"
   stop
  endif
  ! print*,"X= ",x," LS= ", LS
  return
 end subroutine CRYOGENIC_TANK_MK_LS

 ! if SOLID VELOCITY requested everywhere (including outside of the solid),
 ! then velsolid==1
 subroutine CRYOGENIC_TANK_MK_VEL(x,t,LS,VEL,velsolid_flag,dx,nmat)
  use probcommon_module
  IMPLICIT NONE

  INTEGER_T, intent(in) :: nmat
  REAL_T, intent(in) :: x(SDIM)
  REAL_T, intent(in) :: t
  REAL_T, intent(in) :: dx(SDIM)
  REAL_T, intent(in) :: LS(nmat)
  REAL_T, intent(out) :: VEL(SDIM)
  INTEGER_T, intent(in) :: velsolid_flag
  INTEGER_T dir

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
   print *,"velsolid_flag invalid"
   stop
  endif

  if((t.eq.zero).or.(velsolid_flag.eq.0)) then
   do dir=1,SDIM
    VEL(dir)=zero
   enddo
  else
   ! do nothing
  endif

  return 
 end subroutine CRYOGENIC_TANK_MK_VEL

REAL_T function SOLID_TOP_HALF_DIST(P)
 ! Returns the signed distance function to the
 ! cylindrical tank with spherical ends.
 ! The tank is symmetrical to x(SDIM)=0;
 ! The axis of cylinder is along dim=SDIM direction
 ! Inside the tank < 0
 ! Outside the tank > 0
 implicit none

 REAL_T, intent(in), dimension(SDIM) :: P
 REAL_T R,Z,FRZ,D1,D2
 
 if (SDIM.eq.2) then
  R=abs(P(1))
  Z=abs(P(2))
 elseif (SDIM.eq.3) then
  R=abs(sqrt(P(1)**2+P(SDIM)**2))
  Z=abs(P(2))
 else
  print *,"Dimension bust at DIST_FINITE_CYLINDER"
  stop
 endif

  ! R=Z=0 is the dead center of the tank; computation domain goes from
  ! -ZTOTAL/2, ... ,ZTOTAL/2
 if (TANK_MK_END_RADIUS.eq.0.0d0) then ! rectangular tank

  if (Z.le.TANK_MK_HEIGHT/2.0) then
   if (R.ge.TANK_MK_RADIUS) then
    SOLID_TOP_HALF_DIST=R-TANK_MK_RADIUS
   else
    D1 = TANK_MK_RADIUS - R
    D2 = TANK_MK_HEIGHT/2.0-Z
    if (D1.lt.D2) then
     SOLID_TOP_HALF_DIST=-D1
    else
     SOLID_TOP_HALF_DIST=-D2
    endif
   endif
  else 
   if (R.le.TANK_MK_RADIUS) then
    SOLID_TOP_HALF_DIST=Z-TANK_MK_HEIGHT/2.0
   else
    SOLID_TOP_HALF_DIST= &
     sqrt( (Z-TANK_MK_HEIGHT/2.0)**2+ &
           (R-TANK_MK_RADIUS)**2 )
   endif
  endif
    
 else if (TANK_MK_END_RADIUS.gt.0.0d0) then
 
  ! Equation of the line passing through
  ! (0,TANK_MK_END_CENTER)
  !          and
  ! (TANK_MK_RADIUS,TANK_MK_HEIGHT/2)
  FRZ=R*(TANK_MK_END_CENTER-TANK_MK_HEIGHT/2.0d0)/TANK_MK_RADIUS &
     + Z - TANK_MK_END_CENTER

  if (FRZ.le.zero) then
   ! Below line
   if (R.le.TANK_MK_RADIUS) then
    ! In the tank
    D1 = TANK_MK_RADIUS - R
    D2 = sqrt((R-TANK_MK_RADIUS)**2 + (Z-TANK_MK_HEIGHT/2)**2)
    SOLID_TOP_HALF_DIST = -min(D1,D2)
   elseif (R.gt.TANK_MK_RADIUS) then
    ! Out of the tank
    if(Z.le.TANK_MK_HEIGHT/2) then
     ! Below the cap base line
     SOLID_TOP_HALF_DIST = R - TANK_MK_RADIUS
    elseif (Z.gt.TANK_MK_HEIGHT/2) then
     ! Above the cap baseline
     SOLID_TOP_HALF_DIST = &
      sqrt((R-TANK_MK_RADIUS)**2 + (Z-TANK_MK_HEIGHT/2)**2)
    else
     print *,"Z invalid!"
     stop
    endif ! Z
   else
    print *,"Z invalid!"
    stop
   endif ! R
  elseif (FRZ.gt.zero) then
   ! Above line
   ! Distance to curvature center
   D2 = sqrt(R**2 + (Z-TANK_MK_END_CENTER)**2)
   if (D2.le.TANK_MK_END_RADIUS) then
    ! Inside tank
    if (Z.le.TANK_MK_HEIGHT/2) then
     ! Below the cap baseline
     D1 = TANK_MK_RADIUS - R
    elseif (Z.gt.TANK_MK_HEIGHT/2) then
     ! Above the cap baseline
     D1 = sqrt((R-TANK_MK_RADIUS)**2 + (Z-TANK_MK_HEIGHT/2)**2)
    else
     print *,"Z invalid!"
     stop
    endif ! Z
    SOLID_TOP_HALF_DIST = -min(TANK_MK_END_RADIUS-D2,D1)
   elseif (D2.gt.TANK_MK_END_RADIUS) then
    SOLID_TOP_HALF_DIST = D2-TANK_MK_END_RADIUS
   else
    print *,"Distance to curvature center invalid!"
    stop
   end if ! D2
  
  else
   print *,"Line equation invalid!"
   stop
  end if ! FRZ
 else 
  print *,"TANK_MK_END_RADIUS invalid"
  stop
 endif

end function SOLID_TOP_HALF_DIST

!***********************************************
! compressible material functions for 
!   (ns.material_type = TANK_MK_MATERIAL_TYPE)
! C_spc => Specific heat capacity [J(kg K)]
! C_m   =>    Moalr heat capacity [J(mol K)]
!
! U = C_{v,spc} T
! [U] = J/kg= J/(kg K)  K
! R_spc = C_{p,scp}-C_{v,scp}
! R_unv = C_{p,m} - C_{v,m}
! gamma = C_{p,scp}/C_{v,scp}
!
! p = rho R_spc T 
!   = rho R_spc x U/C_{v,spc}
!   = rho (C_{p,scp}-C_{v,scp})/C_{v,spc} U
!   = rhp (gamma-1) U
!
! a = sqrt(gamma R_sp T) = sqrt(gamma p/rho)
!

subroutine EOS_CRYOGENIC_TANK_MK(rho,massfrac_var, &
  internal_energy,pressure, &
  imattype,im,num_species_var_in)
 use probcommon_module
 use global_utility_module
 IMPLICIT NONE
 INTEGER_T, intent(in) :: imattype,im,num_species_var_in
 REAL_T, intent(in) :: rho
 REAL_T, intent(in) :: massfrac_var(num_species_var_in+1)
 REAL_T, intent(in) :: internal_energy
 REAL_T, intent(out) :: pressure
 INTEGER_T :: dummy_input

 if (num_species_var_in.eq.num_species_var) then
  if (im.eq.2) then
   if (imattype.eq.TANK_MK_MATERIAL_TYPE) then
    ! p = rho (gamme-1) U
    pressure=rho * (TANK_MK_GAS_GAMMA-one) * internal_energy
   else
    print *,"imattype= ",imattype
    print *,"imattype invalid EOS_CRYOGENIC_TANK_MK"
    print *,"break point and gdb: "
    print *,"(1) compile with the -g option"
    print *,"(2) break CRYOGENIC_TANK_MK.F90:350"
    print *,"By pressing <CTRL C> during this read statement, the"
    print *,"gdb debugger will produce a stacktrace."
    print *,"type 0 then <enter> to exit the program"
    read *,dummy_input
    error stop
   endif
  else
   call EOS_material_CORE(rho,massfrac_var, &
         internal_energy,pressure,imattype,im)
  endif
 else
  print *,"num_species_var_in invalid"
  stop
 endif

 return
end subroutine EOS_CRYOGENIC_TANK_MK

subroutine dVdT_CRYOGENIC_TANK_MK(dVdT,massfrac_var, &
   pressure,temperature, &
   imattype,im,num_species_var_in)
use probcommon_module
use global_utility_module
IMPLICIT NONE
INTEGER_T, intent(in) :: imattype,im,num_species_var_in
REAL_T, intent(in) :: pressure
REAL_T, intent(in) :: massfrac_var(num_species_var_in+1)
REAL_T, intent(in) :: temperature
REAL_T, intent(out) :: dVdT
INTEGER_T :: dummy_input

 if (pressure.gt.zero) then
  ! do nothing
 else
  print *,"pressure invalid"
  stop
 endif
 if (temperature.gt.zero) then
  ! do nothing
 else
  print *,"temperature invalid"
  stop
 endif

 if (num_species_var_in.eq.num_species_var) then
  if (im.eq.2) then
   if (imattype.eq.TANK_MK_MATERIAL_TYPE) then
    ! p = rho (gamma-1) e  rho=1/V
    ! V = (gamma-1) e/p=(gamma-1)Cv T/p
    ! dVdT=(gamma-1)Cv/p
    dVdT=(TANK_MK_GAS_GAMMA-one) * TANK_MK_GAS_CV/pressure
   else
    print *,"imattype= ",imattype
    print *,"imattype invalid dVdT_CRYOGENIC_TANK_MK"
    print *,"break point and gdb: "
    print *,"(1) compile with the -g option"
    print *,"(2) break CRYOGENIC_TANK_MK.F90:350"
    print *,"By pressing <CTRL C> during this read statement, the"
    print *,"gdb debugger will produce a stacktrace."
    print *,"type 0 then <enter> to exit the program"
    read *,dummy_input
    error stop
   endif
  else
   call dVdT_material_CORE(dVdT,massfrac_var, &
         pressure,temperature, &
         imattype,im)
  endif
 else
  print *,"num_species_var_in invalid"
  stop
 endif

 return
end subroutine dVdT_CRYOGENIC_TANK_MK


subroutine SOUNDSQR_CRYOGENIC_TANK_MK(rho,massfrac_var, &
  internal_energy,soundsqr, &
  imattype,im,num_species_var_in)
 use probcommon_module
 use global_utility_module
 IMPLICIT NONE
 INTEGER_T, intent(in) :: imattype,im,num_species_var_in
 REAL_T, intent(in) :: rho
 REAL_T, intent(in) :: massfrac_var(num_species_var_in+1)
 REAL_T, intent(in) :: internal_energy
 REAL_T, intent(out) :: soundsqr
 REAL_T pressure

 if (num_species_var_in.eq.num_species_var) then
  if (im.eq.2) then
   if (imattype.eq.TANK_MK_MATERIAL_TYPE) then
     ! a = sqrt(gamma R_sp T) = sqrt(gamma p/rho)
    call EOS_CRYOGENIC_TANK_MK(rho,massfrac_var, &
     internal_energy,pressure,imattype,im,num_species_var_in)
    if (rho.gt.zero) then
     soundsqr=TANK_MK_GAS_GAMMA*pressure/rho
    else
     print *,"rho invalid"
     stop
    endif
   else
    print *,"imattype= ",imattype
    print *,"imattype invalid SOUNDSQR CRYOGENIC TANK_MK"
    stop
   endif
  else
   call SOUNDSQR_material_CORE(rho,massfrac_var, &
    internal_energy,soundsqr, &
    imattype,im)
  endif
 else
  print *,"num_species_var_in invalid"
  stop
 endif

 return
end subroutine SOUNDSQR_CRYOGENIC_TANK_MK

subroutine INTERNAL_CRYOGENIC_TANK_MK(rho,massfrac_var, &
  temperature,local_internal_energy, &
  imattype,im,num_species_var_in)
 use probcommon_module
 use global_utility_module
 IMPLICIT NONE
 INTEGER_T, intent(in) :: imattype,im,num_species_var_in
 REAL_T, intent(in) :: rho
 REAL_T, intent(in) :: massfrac_var(num_species_var_in+1)
 REAL_T, intent(in) :: temperature 
 REAL_T, intent(out) :: local_internal_energy

 if (num_species_var_in.eq.num_species_var) then
  if (im.eq.2) then
   if ((imattype.eq.TANK_MK_MATERIAL_TYPE).or. &
       (imattype.eq.0)) then 
    ! U_mix = C_{v,spc} T
    local_internal_energy=TANK_MK_GAS_CV*temperature
   else
    print *,"imattype= ",imattype
    print *,"imattype invalid INTERNAL CRYOGENIC TANK_MK"
    stop
   endif
  else
   call INTERNAL_material_CORE(rho,massfrac_var, &
    temperature,local_internal_energy, &
    imattype,im)
  endif
 else
  print *,"num_species_var_in invalid"
  stop
 endif

 return
end subroutine INTERNAL_CRYOGENIC_TANK_MK

subroutine TEMPERATURE_CRYOGENIC_TANK_MK(rho,massfrac_var, &
  temperature,internal_energy, &
  imattype,im,num_species_var_in)
 use probcommon_module
 use global_utility_module
 IMPLICIT NONE
 INTEGER_T, intent(in) :: imattype,im,num_species_var_in
 REAL_T, intent(in) :: rho
 REAL_T, intent(in) :: massfrac_var(num_species_var_in+1)
 REAL_T, intent(out) :: temperature 
 REAL_T, intent(in) :: internal_energy

 if (num_species_var_in.eq.num_species_var) then
  if (im.eq.2) then
   if ((imattype.eq.TANK_MK_MATERIAL_TYPE).or. &
       (imattype.eq.0)) then 
    ! T = U / C_{v,spc}
    if (TANK_MK_GAS_CV.gt.zero) then
     temperature=internal_energy/TANK_MK_GAS_CV
    else
     print *,"TANK_MK_GAS_CV invalid 1"
     stop
    endif
   else
    print *,"imattype= ",imattype
    print *,"imattype invalid TEMPERATURE_CRYOGENIC_TANK_MK"
    stop
   endif
  else
   call TEMPERATURE_material_CORE(rho,massfrac_var, &
     temperature,internal_energy, &
     imattype,im)
  endif
 else
  print *,"num_species_var_in invalid"
  stop
 endif

 return
end subroutine TEMPERATURE_CRYOGENIC_TANK_MK

!***********************************************
! called by the boundary condition routine
! might be called at initialization, so put a placeholder pressure here.
subroutine CRYOGENIC_TANK_MK_PRES_UTIL(x,PRES,rho_hyd)
use probcommon_module
IMPLICIT NONE

REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(out) :: PRES
REAL_T, intent(out) :: rho_hyd
INTEGER_T simple_hyd_p

simple_hyd_p=1

!PRES=TANK_MK_INITIAL_GAS_PRESSURE
if(fort_material_type(2).eq.0) then

 if (simple_hyd_p.eq.0) then
  rho_hyd=fort_denconst(2)
  ! incompressible gas
  ! Flat open top x_2: TANK_MK_HEIGHT/two
  ! Known pressure(P_1) at top (outflow_pressure)
  ! P_2=P_1 + rho*g*(z_1-z_2)  [g>0]
  if (x(2).ge.TANK_MK_INTERFACE_LOCATION) then
   PRES=TANK_MK_INITIAL_PRESSURE+&
       fort_denconst(2)*(TANK_MK_HEIGHT/two-x(2))*(abs(gravity)) 
  elseif (x(2).lt.TANK_MK_INTERFACE_LOCATION) then
   PRES=TANK_MK_INITIAL_PRESSURE+&
       fort_denconst(2)*(TANK_MK_HEIGHT/two-TANK_MK_INTERFACE_LOCATION)* &
       (abs(gravity))+ &
       fort_denconst(1)*(TANK_MK_INTERFACE_LOCATION-x(2))*(abs(gravity))
  else
   print *,"x(2) is invalid in CRYOGENIC_TANK_MK_PRES!"
   stop
  endif

 else if (simple_hyd_p.eq.1) then
  rho_hyd=fort_denconst(1)
  PRES=-abs(gravity)*rho_hyd*(x(2)-probhiy-probhiy)
 else
  print *,"simple_hyd_p invalid"
  stop
 endif
elseif (fort_material_type(2).eq.TANK_MK_MATERIAL_TYPE) then
 ! compressible gas
 ! Known pressure(P_1) at top (based on given density and temperature)
 ! P_2=P_1 * exp(g*(z_1-z_2)/(R_sp*T_0))  [g>0]
 rho_hyd=fort_denconst(2)
 if (x(2).ge.TANK_MK_INTERFACE_LOCATION) then
  PRES=TANK_MK_INITIAL_PRESSURE*&
       exp((TANK_MK_END_CENTER+TANK_MK_END_RADIUS-x(2))*abs(gravity)/&
           (TANK_MK_R_UNIV/fort_molar_mass(2)*fort_initial_temperature(2)))
 elseif (x(2).lt.TANK_MK_INTERFACE_LOCATION) then
  PRES=TANK_MK_INITIAL_PRESSURE*&
       exp((TANK_MK_END_CENTER+TANK_MK_END_RADIUS-TANK_MK_INTERFACE_LOCATION)*&
            abs(gravity)/&
           (TANK_MK_R_UNIV/fort_molar_mass(2)*fort_initial_temperature(2)))+&
       fort_denconst(1)*(TANK_MK_INTERFACE_LOCATION-x(2))*(abs(gravity))
 else
  print *,"x(2) is invalid in CRYOGENIC_TANK_MK_PRES!"
  stop
 endif
else
 print  *,"invalid material type in pressure setup!"
 stop
endif

return 
end subroutine CRYOGENIC_TANK_MK_PRES_UTIL


subroutine CRYOGENIC_TANK_MK_PRES(x,t,LS,PRES,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(nmat)
REAL_T, intent(out) :: PRES
REAL_T :: rho_hyd

if (num_materials.eq.nmat) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

call CRYOGENIC_TANK_MK_PRES_UTIL(x,PRES,rho_hyd)

return 
end subroutine CRYOGENIC_TANK_MK_PRES



subroutine CRYOGENIC_TANK_MK_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: bcflag !0=called from initialize  1=called from bc
INTEGER_T, intent(in) :: nmat
INTEGER_T, intent(in) :: nstate_mat
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(nmat)
REAL_T, intent(out) :: STATE(nmat*nstate_mat)
INTEGER_T im,ibase,n
REAL_T pressure

 ! num_state_material=2 (default)  density and temperature
 ! num_state_material>2 if scalar (species) variables added.

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

if ((num_materials.eq.3).and. &
    (num_state_material.ge.2).and. &
    (probtype.eq.423)) then
 do im=1,num_materials
  ibase=(im-1)*num_state_material
  ! density
  if(im.eq.2) then
   if(fort_material_type(2).eq.0) then
    ! incompressible
    STATE(ibase+1)=fort_denconst(im)
   elseif(fort_material_type(2).eq.TANK_MK_MATERIAL_TYPE) then
    ! compressible
    ! rho =P/(R_sp T)
    call CRYOGENIC_TANK_MK_PRES(x,t,LS,pressure,nmat)
    STATE(ibase+1) = pressure/&
     (TANK_MK_R_UNIV/fort_molar_mass(2)*fort_initial_temperature(2))
   else
    print *,"material type invalid for density setup!"
    stop
   endif ! material_type
  else if ((im.eq.1).or.(im.eq.3)) then ! liquid or tank walls
   STATE(ibase+1)=fort_denconst(im)
  else
   print *,"im invalid"
   stop
  endif
  ! temperature
  if (t.eq.zero) then
   STATE(ibase+2)=fort_initial_temperature(im)
  else if (t.gt.zero) then
   STATE(ibase+2)=fort_tempconst(im)
  else
   print *,"t invalid"
   stop
  endif
  ! species
  do n=1,num_species_var
   STATE(ibase+2+n)=fort_speciesconst((n-1)*num_materials+im)
  enddo
  if (t.eq.zero) then
   if (im.eq.2) then
    den=STATE(ibase+1)
    temperature=STATE(ibase+2)
    call init_massfrac_parm(den,massfrac_parm,im)
    do n=1,num_species_var
     massfrac_parm(n)=STATE(ibase+2+n)
    enddo
    call INTERNAL_CRYOGENIC_TANK_MK(den,massfrac_parm, &
     temperature,internal_energy,TANK_MK_MATERIAL_TYPE,im,num_species_var)
    call EOS_CRYOGENIC_TANK_MK(den,massfrac_parm, &
     internal_energy,pressure,TANK_MK_MATERIAL_TYPE,im,num_species_var)
    if (abs(TANK_MK_R_UNIV-fort_R_Palmore_Desjardins).le. &
            VOFTOL*TANK_MK_R_UNIV) then
     ! do nothing
    else
     print *,"mismatch between TANK_MK_R_UNIV and fort_R_Palmore_Desjardins"
     stop
    endif
    call Pgamma_Clausius_Clapyron(Pgamma, &
            fort_reference_pressure(1), &
            temperature, &
            fort_saturation_temp(1), &
            fort_latent_heat(1), &
            TANK_MK_R_UNIV,fort_molar_mass(2))
   else if ((im.eq.1).or.(im.eq.3)) then
    ! do nothing
   else
    print *,"im invalid"
    stop
   endif
  else if (t.gt.zero) then
   ! do nothing
  else
   print *,"t invalid"
   stop
  endif
 enddo ! im=1..num_materials
else
 print *,"num_materials,num_state_material, or probtype invalid"
 stop
endif
 
return
end subroutine CRYOGENIC_TANK_MK_STATE

 ! dir=1..sdim  side=1..2
subroutine CRYOGENIC_TANK_MK_LS_BC(xwall,xghost,t,LS, &
   LS_in,dir,side,dx,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: xwall
REAL_T, intent(in) :: xghost(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(inout) :: LS(nmat)
REAL_T, intent(in) :: LS_in(nmat)
INTEGER_T, intent(in) :: dir,side
REAL_T, intent(in) :: dx(SDIM)

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then
 call CRYOGENIC_TANK_MK_LS(xghost,t,LS,nmat)
else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine CRYOGENIC_TANK_MK_LS_BC


 ! dir=1..sdim  side=1..2 veldir=1..sdim
subroutine CRYOGENIC_TANK_MK_VEL_BC(xwall,xghost,t,LS, &
   VEL,VEL_in,veldir,dir,side,dx,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: xwall
REAL_T, intent(in) :: xghost(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(nmat)
REAL_T, intent(inout) :: VEL
REAL_T, intent(in) :: VEL_in
INTEGER_T, intent(in) :: veldir,dir,side
REAL_T, intent(in) :: dx(SDIM)
REAL_T local_VEL(SDIM)
INTEGER_T velsolid_flag

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

 call CRYOGENIC_TANK_MK_VEL(xghost,t,LS,local_VEL,velsolid_flag,dx,nmat)
 VEL=local_VEL(veldir)

else
 print *,"dir,side, or veldir invalid"
 stop
endif

return
end subroutine CRYOGENIC_TANK_MK_VEL_BC

 ! dir=1..sdim  side=1..2
subroutine CRYOGENIC_TANK_MK_PRES_BC(xwall,xghost,t,LS, &
   PRES,PRES_in,dir,side,dx,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: xwall
REAL_T, intent(in) :: xghost(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(nmat)
REAL_T, intent(inout) :: PRES
REAL_T, intent(in) :: PRES_in
INTEGER_T, intent(in) :: dir,side
REAL_T, intent(in) :: dx(SDIM)

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then

 call CRYOGENIC_TANK_MK_PRES(xghost,t,LS,PRES,nmat)

else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine CRYOGENIC_TANK_MK_PRES_BC

 ! dir=1..sdim  side=1..2
subroutine CRYOGENIC_TANK_MK_STATE_BC(xwall,xghost,t,LS, &
   STATE,STATE_merge,STATE_in,im,istate,dir,side,dx,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: xwall
REAL_T, intent(in) :: xghost(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(nmat)
REAL_T local_STATE(nmat*num_state_material)
REAL_T, intent(inout) :: STATE
REAL_T, intent(inout) :: STATE_merge
REAL_T, intent(in) :: STATE_in
INTEGER_T, intent(in) :: dir,side
REAL_T, intent(in) :: dx(SDIM)
INTEGER_T, intent(in) :: istate,im
INTEGER_T ibase,im_crit,im_loop
INTEGER_T local_bcflag

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
 call CRYOGENIC_TANK_MK_STATE(xghost,t,LS,local_STATE, &
         local_bcflag,nmat,num_state_material)
 ibase=(im-1)*num_state_material
 STATE=local_STATE(ibase+istate)
 im_crit=1
 do im_loop=2,num_materials
  if (LS(im_loop).gt.LS(im_crit)) then
   im_crit=im_loop
  endif
 enddo
 ibase=(im_crit-1)*num_state_material
 STATE_merge=local_STATE(ibase+istate)
else
 print *,"istate invalid"
 stop
endif

return
end subroutine CRYOGENIC_TANK_MK_STATE_BC

subroutine CRYOGENIC_TANK_MK_HEATSOURCE( &
     im,VFRAC, &
     time, &
     x, &
     xsten, & ! xsten(-nhalf:nhalf,SDIM)
     nhalf, &
     temp, &
     heat_source, & ! unit of tildeQ not Q
     den,CV,dt, &
     nmat)
use probcommon_module
IMPLICIT NONE

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

INTEGER_T dir
REAL_T local_dx(SDIM)
REAL_T flux_magnitude

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

do dir=1,SDIM
 local_dx(dir)=xsten(1,dir)-xsten(-1,dir)
 if (local_dx(dir).gt.zero) then
  ! do nothing
 else
  print *,"local_dx invalid"
  stop
 endif
enddo

flux_magnitude=zero

if ((num_materials.eq.3).and.(probtype.eq.423)) then
 heat_source=zero
else
 print *,"num_materials ", num_materials
 print *,"probtype ", probtype
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine CRYOGENIC_TANK_MK_HEATSOURCE

! called from subroutine general_hydrostatic_pressure_density which
! is declared in PROB.F90
subroutine CRYOGENIC_TANK_MK_correct_pres_rho_hydrostatic( &
   pres_hydrostatic,rho_hydrostatic, &
   xpos, &
   gravity_normalized, &
   gravity_dir_parm)
use probcommon_module
IMPLICIT NONE

REAL_T, intent(inout) :: rho_hydrostatic
REAL_T, intent(inout) :: pres_hydrostatic
REAL_T, intent(in) :: xpos(SDIM)
REAL_T, intent(in) :: gravity_normalized ! usually |g| (point down case)
INTEGER_T, intent(in) :: gravity_dir_parm

if ((gravity_dir_parm.ge.1).and. &
    (gravity_dir_parm.le.SDIM)) then
 ! do nothing
else
 print *,"gravity_dir_parm invalid"
 stop
endif

if ((num_materials.eq.3).and. &
    (num_state_material.ge.2).and. &
    (probtype.eq.423)) then

 call CRYOGENIC_TANK_MK_PRES_UTIL(xpos,pres_hydrostatic,rho_hydrostatic)

 if (xpos(2).ge.TANK_MK_INTERFACE_LOCATION) then
  if(fort_material_type(2).eq.0) then
   ! incompressible
  elseif(fort_material_type(2).eq.TANK_MK_MATERIAL_TYPE) then
    ! compressible
    ! rho =P/(R_sp T)
   rho_hydrostatic = pres_hydrostatic/&
     (TANK_MK_R_UNIV/fort_molar_mass(2)*fort_initial_temperature(2))
  else
   print *,"material type invalid for pres den hydrostatic!"
   stop
  endif ! material_type

 else if (xpos(2).le.TANK_MK_INTERFACE_LOCATION) then
  if(fort_material_type(2).eq.0) then
   ! incompressible
  elseif(fort_material_type(2).eq.TANK_MK_MATERIAL_TYPE) then
   rho_hydrostatic=fort_denconst(1)
  else
   print *,"material type invalid for pres den hydrostatic!"
   stop
  endif ! material_type
 else
  print *,"xpos(2) invalid"
  stop
 endif

else
 print *,"num_materials,num_state_material, or probtype invalid"
 stop
endif

end subroutine CRYOGENIC_TANK_MK_correct_pres_rho_hydrostatic

subroutine CRYOGENIC_TANK_MK_SUMINT(GRID_DATA_IN,increment_out1, &
                increment_out2,nsum1,nsum2,isweep)
use probcommon_module_types
use probcommon_module
use global_utility_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nsum1,nsum2,isweep
type(user_defined_sum_int_type), intent(in) :: GRID_DATA_IN
REAL_T, intent(inout) :: increment_out1(nsum1)
REAL_T, intent(inout) :: increment_out2(nsum2)

REAL_T massfrac_parm(num_species_var+1)
REAL_T T1_probe(SDIM)
INTEGER_T im
INTEGER_T dir
INTEGER_T dencomp,local_ispec
REAL_T den,temperature,internal_energy,pressure
REAL_T support_r
REAL_T dx_coarsest
REAL_T charfn
REAL_T volgrid
REAL_T denom

INTEGER_T :: i,j,k

i=GRID_DATA_IN%igrid
j=GRID_DATA_IN%jgrid
k=GRID_DATA_IN%kgrid
if ((num_materials.eq.3).and.(probtype.eq.423)) then

 if ((nsum1.eq.1).and.(nsum2.eq.2)) then
  ! integral of region surrounding T1 in Figure 3 of Barsi and Kassemi, 2013
  ! T1: r=0.0  Z=0.2921 relative to bottom of tank cylindrical section.
  ! bottom of cylindrical section: -0.2032
  ! so T1_probe_z=-0.2032+0.2921=0.0889
  T1_probe(1)=zero
  T1_probe(2)=0.0889
  im=2 ! vapor
  dencomp=(im-1)*num_state_material+1
  den=GRID_DATA_IN%den(D_DECL(i,j,k),dencomp)
  temperature=GRID_DATA_IN%den(D_DECL(i,j,k),dencomp+1)
  call init_massfrac_parm(den,massfrac_parm,im)
  do local_ispec=1,num_species_var
   massfrac_parm(local_ispec)= &
       GRID_DATA_IN%den(D_DECL(i,j,k),dencomp+1+local_ispec)
  enddo
  call INTERNAL_CRYOGENIC_TANK_MK(den,massfrac_parm, &
    temperature,internal_energy,TANK_MK_MATERIAL_TYPE,im,num_species_var)
  call EOS_CRYOGENIC_TANK_MK(den,massfrac_parm, &
     internal_energy,pressure,TANK_MK_MATERIAL_TYPE,im,num_species_var)
  support_r=zero
  do dir=1,SDIM
   support_r=support_r+(GRID_DATA_IN%xsten(0,dir)-T1_probe(dir))**2
  enddo
  support_r=sqrt(support_r) 
  dx_coarsest=TANK_MK_HEIGHT/32.0d0
  if (support_r.le.dx_coarsest) then
   charfn=one
  else if (support_r.gt.dx_coarsest) then
   charfn=zero
  else
   print *,"support_r invalid"
   stop
  endif
  volgrid=GRID_DATA_IN%volgrid
  if (isweep.eq.0) then
   increment_out1(1)=charfn*volgrid
  else if (isweep.eq.1) then
   denom=increment_out1(1)
   if (denom.gt.zero) then
    increment_out2(1)=charfn*volgrid*pressure/denom
    increment_out2(2)=charfn*volgrid*temperature/denom
   else
    print *,"expecting denom>0.0"
    stop
   endif
  else
   print *,"isweep invalid"
   stop
  endif
 else
  print *,"nsum1 or nsum2 invalid"
  stop
 endif

else
 print *,"num_materials ", num_materials
 print *,"probtype ", probtype
 print *,"num_materials or probtype invalid"
 stop
endif

end subroutine CRYOGENIC_TANK_MK_SUMINT

subroutine CRYOGENIC_TANK_MK_INIT_REGIONS_LIST(constant_density_all_time, &
      num_materials_in,num_threads_in)
use probcommon_module
use geometry_intersect_module

IMPLICIT NONE

INTEGER_T, intent(in) :: num_materials_in
INTEGER_T, intent(in) :: num_threads_in
INTEGER_T, intent(in) :: constant_density_all_time(num_materials_in)
INTEGER_T :: im,iregion

 if (num_materials_in.eq.num_materials) then
  ! do nothing
 else
  print *,"num_materials_in invalid"
  stop
 endif
 if (num_threads_in.eq.geom_nthreads) then
  ! do nothing
 else
  print *,"num_threads_in invalid"
  stop
 endif
 do im=1,num_materials
  if ((constant_density_all_time(im).eq.0).or. &
      (constant_density_all_time(im).eq.1)) then
   ! do nothing
  else
   print *,"constant_density_all_time(im) invalid"
   stop
  endif
 enddo ! im=1..num_materials

 number_of_source_regions=1
 number_of_threads_regions=num_threads_in
 allocate(regions_list(1:number_of_source_regions, &
                       0:number_of_threads_regions))

 do iregion=1,number_of_source_regions
  regions_list(iregion,0)%region_material_id=0
  regions_list(iregion,0)%region_dt=zero
  regions_list(iregion,0)%region_mass_flux=zero
  regions_list(iregion,0)%region_volume_flux=zero
  regions_list(iregion,0)%region_energy_flux=zero
  regions_list(iregion,0)%region_volume_raster=zero 
  regions_list(iregion,0)%region_volume=zero 
  regions_list(iregion,0)%region_mass=zero 
  regions_list(iregion,0)%region_energy=zero 
  regions_list(iregion,0)%region_energy_per_kelvin=zero 
  regions_list(iregion,0)%region_volume_after=zero 
  regions_list(iregion,0)%region_mass_after=zero 
  regions_list(iregion,0)%region_energy_after=zero 
 enddo ! iregion=1,number_of_source_regions

 regions_list(1,0)%region_material_id=1
 regions_list(1,0)%region_energy_flux=TANK_MK_HEATER_FLUX ! Watts=J/s

end subroutine CRYOGENIC_TANK_MK_INIT_REGIONS_LIST

subroutine CRYOGENIC_TANK_MK_CHARFN_REGION(region_id,x,cur_time,charfn_out)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: region_id
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: cur_time
REAL_T, intent(out) :: charfn_out

if ((num_materials.eq.3).and.(probtype.eq.423)) then
 if (region_id.eq.1) then
  if ((abs(x(1)).le.TANK_MK_HEATER_R).and.&
      (abs(x(1)).ge.TANK_MK_HEATER_R_LOW).and.&
      (x(2).ge.TANK_MK_HEATER_LOW).and.&
      (x(2).le.TANK_MK_HEATER_HIGH)) then
   charfn_out=one
  else if ((abs(x(1)).gt.TANK_MK_HEATER_R).or. &
           (abs(x(1)).lt.TANK_MK_HEATER_R_LOW).or. &
           (x(2).lt.TANK_MK_HEATER_LOW).or. &
           (x(2).gt.TANK_MK_HEATER_HIGH)) then
   charfn_out=zero
  else
   print *,"position bust"
   stop
  endif

 else
  print *,"region_id invalid"
  stop
 endif
else
 print *,"num_materials ", num_materials
 print *,"probtype ", probtype
 print *,"num_materials or probtype invalid"
 stop
endif

end subroutine CRYOGENIC_TANK_MK_CHARFN_REGION

end module CRYOGENIC_TANK_MK_module
