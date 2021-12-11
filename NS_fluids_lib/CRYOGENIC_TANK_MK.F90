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
REAL_T :: TANK_MK_HEATER_WATTS
! Heater location in dim=2 direction
REAL_T :: TANK_MK_HEATER_WALL_MODEL
REAL_T :: TANK_MK_HEATER_LOW
REAL_T :: TANK_MK_HEATER_HIGH
REAL_T :: TANK_MK_HEATER_R
REAL_T :: TANK_MK_HEATER_R_LOW

REAL_T :: TANK_MK_HEATER_THICK
REAL_T :: TANK_MK_HEATER_TOTAL_ANGLE

REAL_T :: TANK_MK_INSULATE_R
REAL_T :: TANK_MK_INSULATE_R_HIGH
REAL_T :: TANK_MK_INSULATE_THICK

REAL_T :: TANK_MK_NOZZLE_RAD
REAL_T :: TANK_MK_NOZZLE_HT
REAL_T :: TANK_MK_NOZZLE_THICK_OUTLET
REAL_T :: TANK_MK_NOZZLE_BASE

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

  TANK_MK_END_RADIUS       = xblob4
  TANK_MK_END_CENTER       = yblob4

  TANK_MK_INTERFACE_RADIUS = radblob2
  TANK_MK_BUBBLE_X         = xblob2
  TANK_MK_BUBBLE_Y         = yblob2
  TANK_MK_BUBBLE_Z         = zblob2

  TANK_MK_HEATER_WATTS      = xblob3

   ! see Barsi and Kassemi 2013, Journal of Thermal Science and Engineering
   ! Applications.
   ! ZBOT
  if (axis_dir.eq.0) then ! volume=pi(.09^2-0.08^2)*.02=pi(.01)(.17)(0.02)
   TANK_MK_HEATER_WALL_MODEL = 0.1683d0
   TANK_MK_HEATER_LOW       = -0.1683d0
   TANK_MK_HEATER_HIGH      = TANK_MK_HEATER_LOW+0.0254d0
   TANK_MK_HEATER_R_LOW     = 0.1016d0
   TANK_MK_HEATER_R         = TANK_MK_HEATER_R_LOW+0.027d0

   TANK_MK_INSULATE_R = xblob+0.027d0 !xblob is the tank cavity radius
   TANK_MK_INSULATE_THICK = 0.0508
   TANK_MK_INSULATE_R_HIGH = yblob/2.0d0 !yblob is height of cylindrical part

   ! TPCE
  else if (axis_dir.eq.1) then ! heater on top
   TANK_MK_HEATER_WALL_MODEL = 0.0
   TANK_MK_HEATER_THICK     = 0.01d0
   TANK_MK_HEATER_TOTAL_ANGLE = 30.0d0*Pi/180.0d0

   TANK_MK_NOZZLE_RAD=0.005D0  !dx coarsest =0.0015625, "1cm diameter."
   TANK_MK_NOZZLE_HT=0.064D0
   TANK_MK_NOZZLE_THICK_OUTLET=0.005d0
   TANK_MK_NOZZLE_BASE=-half*TANK_MK_HEIGHT-TANK_MK_END_RADIUS

  else
   print *,"axis_dir invalid"
   stop
  endif

  ! ASSUMING IDEAL GAS => The gas heat cpacities should satisfy this
  ! R_spc = C_{p,spc}-C_{v,spc}
  ! to have ideal mixture gas as well =>
  ! Only C_p or C_v can be picked from table and the other one
  ! calculated from equation above.
  ! Here we pick C_{v,spc} from input file.
  ! C_{p,spc} = C_{v,spc} + R_spc

  TANK_MK_R_UNIV = 8.31446261815324D0

  TANK_MK_GAS_CV = fort_stiffCV(2) ![J/(kg K)]
!  TANK_MK_GAS_CP = fort_stiffCP(2) ![J∕(kg·K)]
  TANK_MK_GAS_CP=TANK_MK_GAS_CV+TANK_MK_R_UNIV/fort_molar_mass(2) ![J∕(kg·K)]
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

  ! LS>0 in the nozzle 
 subroutine CRYOGENIC_TANK_MK_LS_NOZZLE(x,LS)
  use probcommon_module
  use global_utility_module
  IMPLICIT NONE

  REAL_T, intent(in) :: x(SDIM)
  REAL_T, intent(out) :: LS
  REAL_T :: xlo,xhi,ylo,yhi

  if ((TANK_MK_NOZZLE_RAD.gt.0.0d0).and. &
      (TANK_MK_NOZZLE_BASE.lt.0.0d0).and. &
      (TANK_MK_NOZZLE_HT.gt.0.0d0).and. &
      (TANK_MK_NOZZLE_THICK_OUTLET.gt.0.0d0)) then
   xlo=-TANK_MK_NOZZLE_RAD
   xhi=TANK_MK_NOZZLE_RAD
   ! "TANK_MK_HEIGHT" term insures nozzle is flush against the
   ! bottom of the tank.
   ylo=TANK_MK_NOZZLE_BASE-TANK_MK_HEIGHT
   yhi=TANK_MK_NOZZLE_BASE+TANK_MK_NOZZLE_HT
   call squaredist(x(1),x(2),xlo,xhi,ylo,yhi,LS)
   LS=-LS !now, LS>0 in the nozzle.
  else
   print *,"CRYOGENIC_TANK_MK_LS_NOZZLE parameter problem"
   stop
  endif

  return
 end subroutine CRYOGENIC_TANK_MK_LS_NOZZLE

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
  use global_utility_module
  IMPLICIT NONE

  INTEGER_T, intent(in) :: nmat
  REAL_T, intent(in) :: x(SDIM)
  REAL_T, intent(in) :: t
  REAL_T, intent(out) :: LS(nmat)
  REAL_T :: nozzle_dist,LS_A
  INTEGER_T :: caller_id

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
   if (TANK_MK_INTERFACE_RADIUS.eq.0.0d0) then
    LS(1)=TANK_MK_INTERFACE_LOCATION-x(SDIM)
   else if (TANK_MK_INTERFACE_RADIUS.gt.0.0d0) then
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

   if (axis_dir.eq.0) then
    ! do nothing
   else if (axis_dir.eq.1) then ! TPCE
    call CRYOGENIC_TANK_MK_LS_NOZZLE(x,nozzle_dist)
    if (nozzle_dist.gt.LS(3)) then ! nozzle_dist>0 in the nozzle
     LS(3)=nozzle_dist
    endif
    caller_id=0
    call CRYOGENIC_TANK_MK_LS_HEATER_A(x,LS_A,caller_id) !LS_A>0 in heater
    if (LS_A.gt.LS(3)) then
     LS(3)=LS_A
    endif
   else
    print *,"axis_dir invalid"
    stop
   endif

  else
   print *,"num_materials ", num_materials
   print *,"probtype ", probtype
   print *,"num_materials or probtype invalid"
   stop
  endif
  ! print*,"X= ",x," LS= ", LS
  return
 end subroutine CRYOGENIC_TANK_MK_LS

   !LS>0 in heater A region.
   !caller_id=0 if called from geometry routine
   !caller_id=1 if called from heater source routine.
 subroutine CRYOGENIC_TANK_MK_LS_HEATER_A(x,LS,caller_id)
  use probcommon_module
  use global_utility_module
  IMPLICIT NONE

  REAL_T, intent(in) :: x(SDIM)
  REAL_T, intent(out) :: LS
  INTEGER_T, intent(in) :: caller_id
  REAL_T zdiff
  REAL_T angle_end_center
  REAL_T shell_R
  REAL_T shell_center
  REAL_T r_crit
  REAL_T z_crit

  zdiff=x(2)-TANK_MK_END_CENTER

  if (zdiff.le.0.0d0) then
   LS=-TANK_MK_END_RADIUS
  else if (zdiff.gt.0.0d0) then
   angle_end_center=atan(abs(x(1))/zdiff)
   if ((angle_end_center.ge.0.0d0).and. &
       (angle_end_center.lt.0.5d0*Pi)) then
    shell_R=sqrt(x(1)**2+zdiff**2)
    shell_center=TANK_MK_END_RADIUS-0.5d0*TANK_MK_HEATER_THICK

    if (caller_id.eq.1) then
     LS=0.5d0*TANK_MK_HEATER_THICK-abs(shell_R-shell_center)
    else if (caller_id.eq.0) then
     LS=shell_R-TANK_MK_END_RADIUS+TANK_MK_HEATER_THICK
    else
     print *,"caller_id invalid"
     stop
    endif

    if ((TANK_MK_HEATER_TOTAL_ANGLE.gt.0.0d0).and. &
        (TANK_MK_HEATER_TOTAL_ANGLE.le.0.5d0*Pi)) then
     if (angle_end_center.le.TANK_MK_HEATER_TOTAL_ANGLE) then
      ! do nothing
     else if (angle_end_center.ge.TANK_MK_HEATER_TOTAL_ANGLE) then
      r_crit=sin(TANK_MK_HEATER_TOTAL_ANGLE)*(TANK_MK_END_RADIUS- &
            0.5d0*TANK_MK_HEATER_THICK)
      z_crit=cos(TANK_MK_HEATER_TOTAL_ANGLE)*(TANK_MK_END_RADIUS- &
            0.5d0*TANK_MK_HEATER_THICK)
      if ((r_crit.gt.0.0d0).and.(z_crit.gt.0.0d0).and. &
          (r_crit.le.TANK_MK_END_RADIUS).and. &
          (z_crit.le.TANK_MK_END_RADIUS)) then
       LS=-sqrt( (abs(x(1))-r_crit)**2+(zdiff-z_crit)**2 )
      else
       print *,"r_crit or z_crit invalid"
       stop
      endif
     else
      print *,"angle_end_center is NaN"
      stop
     endif 
    else
     print *,"TANK_MK_HEATER_TOTAL_ANGLE invalid"
     stop
    endif
   else
    print *,"angle_end_center invalid"
    stop
   endif
  else
   print *,"zdiff is NaN"
   stop
  endif

  return
 end subroutine CRYOGENIC_TANK_MK_LS_HEATER_A

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

  if((t.eq.0.0d0).or.(velsolid_flag.eq.0)) then
   do dir=1,SDIM
    VEL(dir)=0.0d0
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

  if (FRZ.le.0.0d0) then
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
  elseif (FRZ.gt.0.0d0) then
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
! C_m   =>    Molar heat capacity [J(mol K)]
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
    print *,"(breakpoint) break point and gdb: "
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

 if (pressure.gt.0.0d0) then
  ! do nothing
 else
  print *,"pressure invalid dVdT_CRYOGENIC_TANK_MK:",pressure
  stop
 endif
 if (temperature.gt.0.0d0) then
  ! do nothing
 else
  print *,"temperature invalid dVdT_CRYOGENIC_TANK_MK:",temperature
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
    print *,"(breakpoint) break point and gdb: "
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
    if (rho.gt.0.0d0) then
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
    if (TANK_MK_GAS_CV.gt.0.0d0) then
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
use global_utility_module
IMPLICIT NONE

INTEGER_T, intent(in) :: bcflag !0=called from initialize  1=called from bc
INTEGER_T, intent(in) :: nmat
INTEGER_T, intent(in) :: nstate_mat
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(nmat)
REAL_T, intent(out) :: STATE(nmat*nstate_mat)
INTEGER_T im,ibase,n
REAL_T den,temperature,internal_energy,pressure,Pgamma
REAL_T massfrac_parm(num_species_var+1)

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
   print *,"im invalid 758 ",im
   stop
  endif
  ! temperature
  if (t.eq.0.0d0) then
   STATE(ibase+2)=fort_initial_temperature(im)
  else if (t.gt.0.0d0) then
   STATE(ibase+2)=fort_tempconst(im)
  else
   print *,"t invalid"
   stop
  endif
  ! species
  do n=1,num_species_var
   STATE(ibase+2+n)=fort_speciesconst((n-1)*num_materials+im)
  enddo
  if (t.eq.0.0d0) then
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
    if (abs(Pgamma-pressure).le.VOFTOL*pressure) then
     ! do nothing
    else
     print *,"mismatch between Pgamma and Pgas"
     print *,"Pgamma=",Pgamma
     print *,"Pgas=",pressure
     print *,"reference pressure=",fort_reference_pressure(1)
     print *,"modify reference pressure to be: ", &
        fort_reference_pressure(1)*pressure/Pgamma
     stop
    endif
             
   else if ((im.eq.1).or.(im.eq.3)) then
    ! do nothing
   else
    print *,"im invalid 814 ",im
    stop
   endif
  else if (t.gt.0.0d0) then
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
 if (local_dx(dir).gt.0.0d0) then
  ! do nothing
 else
  print *,"local_dx invalid"
  stop
 endif
enddo

flux_magnitude=0.0d0

if ((num_materials.eq.3).and.(probtype.eq.423)) then
 heat_source=0.0d0
else
 print *,"num_materials ", num_materials
 print *,"probtype ", probtype
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine CRYOGENIC_TANK_MK_HEATSOURCE

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
REAL_T T1_probe(SDIM)  !ZBOT
REAL_T T4_probe(SDIM)  !TPCE
INTEGER_T im
INTEGER_T dir
INTEGER_T dencomp,local_ispec
REAL_T den,temperature,internal_energy,pressure
REAL_T support_r
REAL_T dx_coarsest
REAL_T charfn
REAL_T volgrid
REAL_T denom

INTEGER_T :: level,finest_level

INTEGER_T :: i,j,k
INTEGER_T :: ilev

i=GRID_DATA_IN%igrid
j=GRID_DATA_IN%jgrid
k=GRID_DATA_IN%kgrid
level=GRID_DATA_IN%level
finest_level=GRID_DATA_IN%finest_level

if ((level.le.finest_level).and.(level.ge.0)) then
 ! do nothing
else
 print *,"level invalid"
 stop
endif

if ((num_materials.eq.3).and.(probtype.eq.423)) then

 if ((nsum1.eq.1).and.(nsum2.eq.2)) then
  ! integral of region surrounding T1 in Figure 3 of Barsi and Kassemi, 2013
  ! T1: r=0.0  Z=0.2921 relative to bottom of tank cylindrical section.
  ! (note the very bottom of the tank corresponds to z=-0.2032-0.0254)
  ! bottom of cylindrical section: -0.2032
  ! so T1_probe_z=-0.2032+0.2921=0.0889
  T1_probe(1)=0.0d0
  T1_probe(2)=0.0889d0
  T4_probe(1)=0.0d0
  T4_probe(2)=half*TANK_MK_HEIGHT+TANK_MK_END_RADIUS-0.025d0
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
  support_r=0.0d0
  do dir=1,SDIM
   if (axis_dir.eq.0) then
    support_r=support_r+(GRID_DATA_IN%xsten(0,dir)-T1_probe(dir))**2
   else if (axis_dir.eq.1) then
    support_r=support_r+(GRID_DATA_IN%xsten(0,dir)-T4_probe(dir))**2
   else
    print *,"axis_dir invalid"
    stop
   endif
  enddo
  support_r=sqrt(support_r) 
  dx_coarsest=GRID_DATA_IN%dx(SDIM)
  do ilev=0,level-1
   dx_coarsest=2.0d0*dx_coarsest
  enddo

  if (support_r.le.2.0d0*dx_coarsest) then
   charfn=one
  else if (support_r.gt.2.0d0*dx_coarsest) then
   charfn=0.0d0
  else
   print *,"support_r invalid"
   stop
  endif
  volgrid=GRID_DATA_IN%volgrid
  if (isweep.eq.0) then
   increment_out1(1)=charfn*volgrid
   if (1.eq.0) then
    print *,"nsum1,nsum2 ",nsum1,nsum2
    print *,"charfn,volgrid,pressure,temperature ", &
       charfn,volgrid,pressure,temperature
    print *,"i,j,k ",i,j,k
   endif

  else if (isweep.eq.1) then
   denom=increment_out1(1)
   if (denom.gt.0.0d0) then
    increment_out2(1)=charfn*volgrid*pressure/denom
    increment_out2(2)=charfn*volgrid*temperature/denom
   else
    print *,"expecting denom>0.0:",denom
    print *,"nsum1,nsum2 ",nsum1,nsum2
    print *,"charfn,volgrid,pressure,temperature ", &
       charfn,volgrid,pressure,temperature
    print *,"i,j,k ",i,j,k
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

IMPLICIT NONE

INTEGER_T, intent(in) :: num_materials_in
INTEGER_T, intent(in) :: num_threads_in
INTEGER_T, intent(in) :: constant_density_all_time(num_materials_in)
INTEGER_T :: im,iregion,dir

 if (num_materials_in.eq.num_materials) then
  ! do nothing
 else
  print *,"num_materials_in invalid"
  stop
 endif
 if (num_threads_in.ge.1) then
  ! do nothing
 else
  print *,"num_threads_in invalid: ",num_threads_in
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

 if (axis_dir.eq.0) then
  number_of_source_regions=1 ! side heater
 else if (axis_dir.eq.1) then
  number_of_source_regions=3 ! heater A, inflow, outflow
 else
  print *,"axis_dir invalid"
  stop
 endif

 number_of_threads_regions=num_threads_in
 allocate(regions_list(1:number_of_source_regions, &
                       0:number_of_threads_regions))

 do iregion=1,number_of_source_regions
  regions_list(iregion,0)%region_material_id=0
  regions_list(iregion,0)%region_dt=0.0d0
  regions_list(iregion,0)%region_mass_flux=0.0d0
  regions_list(iregion,0)%region_volume_flux=0.0d0
  regions_list(iregion,0)%region_temperature_prescribe=0.0d0
  do dir=1,SDIM
   regions_list(iregion,0)%region_velocity_prescribe(dir)=0.0d0
  enddo
  regions_list(iregion,0)%region_energy_flux=0.0d0
  regions_list(iregion,0)%region_volume_raster=0.0d0 
  regions_list(iregion,0)%region_volume=0.0d0 
  regions_list(iregion,0)%region_mass=0.0d0 
  regions_list(iregion,0)%region_energy=0.0d0 
  regions_list(iregion,0)%region_energy_per_kelvin=0.0d0 
  regions_list(iregion,0)%region_volume_after=0.0d0 
  regions_list(iregion,0)%region_mass_after=0.0d0 
  regions_list(iregion,0)%region_energy_after=0.0d0 
 enddo ! iregion=1,number_of_source_regions

 if (axis_dir.eq.0) then
  regions_list(1,0)%region_material_id=3
  regions_list(1,0)%region_energy_flux=TANK_MK_HEATER_WATTS ! Watts=J/s
 else if (axis_dir.eq.1) then
  regions_list(1,0)%region_material_id=3
  regions_list(1,0)%region_energy_flux=TANK_MK_HEATER_WATTS ! Watts=J/s
   ! inflow
  regions_list(2,0)%region_material_id=1
  regions_list(2,0)%region_volume_flux=xblob5
  regions_list(2,0)%region_mass_flux=xblob5*fort_denconst(1)
  regions_list(2,0)%region_temperature_prescribe=xblob6
  if (TANK_MK_NOZZLE_RAD.gt.0.0d0) then
   regions_list(2,0)%region_velocity_prescribe(SDIM)= &
      xblob5/(Pi*(TANK_MK_NOZZLE_RAD**2.0d0))
  else
   print *,"TANK_MK_NOZZLE_RAD invalid"
   stop
  endif
   ! outflow
  regions_list(3,0)%region_material_id=1
  regions_list(3,0)%region_volume_flux=-xblob5
  regions_list(3,0)%region_mass_flux=-xblob5*fort_denconst(1)
 else
  print *,"axis_dir invalid"
  stop
 endif
 
end subroutine CRYOGENIC_TANK_MK_INIT_REGIONS_LIST

subroutine CRYOGENIC_TANK_MK_CHARFN_REGION(region_id,x,cur_time,charfn_out)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: region_id
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: cur_time
REAL_T, intent(out) :: charfn_out
REAL_T :: TANK_MK_R_WIDTH
REAL_T :: shell_R,shell_center,LS_SHELL,LS_A,LS_nozzle,zdiff
INTEGER_T :: caller_id

if ((num_materials.eq.3).and.(probtype.eq.423)) then

 if (axis_dir.eq.0) then ! ZBOT
  TANK_MK_R_WIDTH=TANK_MK_HEATER_R-TANK_MK_HEATER_R_LOW
  if (TANK_MK_R_WIDTH.gt.0.0d0) then
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
     charfn_out=0.0d0
    else
     print *,"position bust"
     stop
    endif
   else
    print *,"region_id invalid"
    stop
   endif
  else
   print *,"TANK_MK_R_WIDTH invalid"
   stop
  endif

 else if (axis_dir.eq.1) then !TPCE

  if (region_id.eq.1) then
   caller_id=1
   call CRYOGENIC_TANK_MK_LS_HEATER_A(x,LS_A,caller_id) !LS_A>0 in heater
   if (LS_A.ge.0.0d0) then
    charfn_out=one
   else if (LS_A.le.0.0d0) then
    charfn_out=0.0d0
   else
    print *,"LS_A invalid"
    stop
   endif
  else if (region_id.eq.2) then ! inflow
   if ((TANK_MK_NOZZLE_RAD.gt.0.0d0).and. &
       (TANK_MK_NOZZLE_BASE.lt.0.0d0).and. &
       (TANK_MK_NOZZLE_HT.gt.0.0d0).and. &
       (TANK_MK_NOZZLE_THICK_OUTLET.gt.0.0d0)) then
    if ((abs(x(1)).le.TANK_MK_NOZZLE_RAD).and. &
        (x(2).gt.TANK_MK_NOZZLE_BASE+TANK_MK_NOZZLE_HT).and. &
        (x(2).le.TANK_MK_NOZZLE_BASE+TANK_MK_NOZZLE_HT+ &
                 TANK_MK_NOZZLE_THICK_OUTLET)) then
     charfn_out=one
    else
     charfn_out=0.0d0
    endif
   else
    print *,"region_id==2 parameter problem"
    stop
   endif
  else if (region_id.eq.3) then ! outflow
   call CRYOGENIC_TANK_MK_LS_NOZZLE(x,LS_nozzle)
   if (LS_nozzle.ge.0.0d0) then
    charfn_out=0.0d0
   else if (LS_nozzle.le.0.0d0) then
    if ((TANK_MK_END_CENTER.gt.0.0d0).and. &
        (TANK_MK_HEATER_THICK.gt.0.0d0).and. &
        (TANK_MK_END_RADIUS.gt.0.0d0)) then
     zdiff=x(2)+TANK_MK_END_CENTER
     if (zdiff.ge.0.0d0) then
      charfn_out=0.0d0
     else if (zdiff.le.-TANK_MK_END_RADIUS) then
      charfn_out=0.0d0
     else if ((zdiff.lt.0.0d0).and. &
              (zdiff.gt.-TANK_MK_END_RADIUS)) then
      shell_R=sqrt(x(1)**2+zdiff**2)
      shell_center=TANK_MK_END_RADIUS-0.5d0*TANK_MK_HEATER_THICK
      LS_SHELL=0.5d0*TANK_MK_HEATER_THICK-abs(shell_R-shell_center)
      if (LS_SHELL.ge.0.0d0) then
       charfn_out=one
      else if (LS_SHELL.le.0.0d0) then
       charfn_out=0.0d0
      else
       print *,"LS_SHELL invalid"
       stop
      endif
     else
      print *,"zdiff is NaN"
      stop
     endif
    else
     print *,"parameter bust"
     stop
    endif
   else
    print *,"LS_nozzle is NaN"
    stop
   endif

  else
   print *,"region_id invalid"
   stop
  endif

 else
  print *,"axis_dir invalid"
  stop
 endif

else
 print *,"num_materials ", num_materials
 print *,"probtype ", probtype
 print *,"num_materials or probtype invalid"
 stop
endif

end subroutine CRYOGENIC_TANK_MK_CHARFN_REGION


subroutine CRYOGENIC_TANK_MK_THERMAL_K(x,dx,cur_time, &
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
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: im
INTEGER_T, intent(in) :: im_solid
INTEGER_T, intent(in) :: near_interface
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: dx(SDIM)
REAL_T, intent(in) :: cur_time
REAL_T, intent(in) :: density
REAL_T, intent(in) :: temperature
REAL_T, intent(in) :: temperature_wall
REAL_T, intent(in) :: temperature_wall_max
REAL_T, intent(in) :: temperature_probe
REAL_T, intent(in) :: nrm(SDIM) ! nrm points from solid to fluid
REAL_T, intent(inout) :: thermal_k
REAL_T :: Ra,Gr,Pr,psi,alpha
REAL_T :: mu_w,rho_w,nu,thermal_conductivity,Cp,xi,R,thermal_diffusivity
REAL_T :: gravity_local,expansion_coefficient
INTEGER_T :: turb_flag

REAL_T :: LS_A
INTEGER_T :: caller_id

if (probtype.eq.423) then
 ! do nothing
else
 print *,"probtype invalid"
 stop
endif

if ((im_solid.ge.0).and.(im_solid.le.num_materials)) then
 ! do nothing
else
 print *,"im_solid invalid"
 stop
endif

if ((im.ge.1).and.(im.le.num_materials)) then

 if (fort_thermal_microlayer_size(im).gt.zero) then
  !do nothing
 else
  print *,"thermal_microlayer_size(im) invalid"
  stop
 endif  

 if (axis_dir.eq.0) then ! ZBOT

  if (im.eq.2) then ! vapor
   ! do nothing
  else if ((im.eq.1).or.(im.eq.3)) then ! liquid or solid
   if ((abs(x(1)).le.TANK_MK_HEATER_R).and.&
       (abs(x(1)).ge.TANK_MK_HEATER_R_LOW-dx(1)).and.&
       (x(2).ge.TANK_MK_HEATER_LOW).and.&
       (x(2).le.TANK_MK_HEATER_HIGH)) then
    thermal_k=fort_heatviscconst(im)* &
      max(one,dx(1)/fort_thermal_microlayer_size(im))
   else if ((abs(x(2)).ge.TANK_MK_INSULATE_THICK+TANK_MK_HEIGHT/2.0d0).or. &
            (abs(x(1)).ge.TANK_MK_INSULATE_R_HIGH)) then
    if (im.eq.3) then
     thermal_k=0.0d0
    else if (im.eq.1) then
     ! do nothing
    else
     print *,"im invalid"
     stop
    endif
   else if ((abs(x(2)).le.TANK_MK_HEIGHT/2.0d0).and. &
            (abs(x(1)).ge.TANK_MK_INSULATE_R)) then
    if (im.eq.3) then
     thermal_k=0.0d0
    else if (im.eq.1) then
     ! do nothing
    else
     print *,"im invalid"
     stop
    endif
   else if (im.eq.3) then
    ! do nothing
   else if (im.eq.1) then
    mu_w=fort_viscconst(im) 
    rho_w=fort_denconst(im)
    if (rho_w.gt.zero) then
     ! do nothing
    else
     print *,"rho_w invalid"
     stop
    endif
    nu=mu_w/rho_w

    thermal_conductivity=fort_heatviscconst(im)
    Cp=fort_stiffCP(im)
    if (thermal_conductivity.gt.zero) then
     ! do nothing
    else
     print *,"expecting thermal_conductivity to be positive"
     stop
    endif
    if (Cp.gt.zero) then
     ! do nothing
    else
     print *,"Cp invalid"
     stop
    endif
    xi=x(SDIM)-TANK_MK_HEATER_LOW
    R=x(1)
    if ((xi.gt.0.0d0).and. &
        (xi.le.TANK_MK_HEATER_WALL_MODEL).and. &
        (nrm(1).eq.-one).and. &
        (near_interface.eq.1).and. &
        (temperature_wall_max.gt.temperature_probe)) then
     thermal_diffusivity=thermal_conductivity/(rho_w*Cp)
     gravity_local=abs(gravity)
     expansion_coefficient=abs(fort_DrhoDT(im)) ! units: 1/temperature
     Gr=gravity_local*expansion_coefficient* &
       (temperature_wall_max-temperature_probe)*(xi**3.0)/(nu*nu)
     Pr=nu/thermal_diffusivity
     Ra=Gr*Pr

     if(Ra.lt.1.0e+9)then
      turb_flag=0
     else if(Ra.ge.1.0e+9)then
      turb_flag=1
     else
      print *,"Ra number invalid"
      stop
     endif
     psi=(1.0d0+(0.492/Pr)**(9.0d0/16.0d0))**(-16.0d0/9.0d0)
     if (turb_flag.eq.0) then
      alpha=(thermal_conductivity/xi)* &
            (0.68d0+0.503d0*((Ra*psi)**0.25d0))
     else if (turb_flag.eq.1) then
      alpha=(thermal_conductivity/xi)* &
            (0.15d0*((Ra*psi)**(1.0d0/3.0d0)))
     else
      print *,"turb_flag invalid"
      stop
     endif 
     thermal_k=max(thermal_k,alpha*dx(1))
    else if ((xi.le.0.0d0).or. &
             (xi.ge.TANK_MK_HEATER_WALL_MODEL).or. &
             (nrm(1).ne.-one).or. &
             (near_interface.eq.0).or. &
             (temperature_wall_max.le.temperature_probe)) then
     ! do nothing
    else
     print *,"xi,nrm,temp_wall, or temp_probe invalid"
     stop
    endif
   else
    print *,"im must be 1 or 3"
    stop
   endif

  else
   print *,"im invalid"
   stop
  endif

 else if (axis_dir.eq.1) then !TPCE
  if (im.eq.2) then ! vapor
   ! do nothing
  else if ((im.eq.1).or.(im.eq.3)) then ! liquid or solid
   caller_id=1
   call CRYOGENIC_TANK_MK_LS_HEATER_A(x,LS_A,caller_id)
   if (LS_A.gt.-dx(1)) then
    thermal_k=fort_heatviscconst(im)* &
      max(one,dx(1)/fort_thermal_microlayer_size(im))
   else if (LS_A.le.-dx(1)) then
    if (im.eq.3) then
     thermal_k=0.0d0
    else if (im.eq.1) then
     ! do nothing
    else
     print *,"im invalid"
     stop
    endif
   else
    print *,"LS_A is NaN"
    stop
   endif
  else
   print *,"im invalid"
   stop
  endif

 else
  print *,"axis_dir out of range"
  stop
 endif

else 
 print *,"im invalid in CRYOGENIC_TANK_MK_THERMAL_K"
 stop
endif

end subroutine CRYOGENIC_TANK_MK_THERMAL_K

! Temperature Stratification in a Cryogenic Fuel Tank
! JOURNAL OF THERMOPHYSICS AND HEAT TRANSFER
! Vol. 27, No. 1, January–March 2013
! Daigle, Smelyanskiy, Boschee, Foygel
!   also
! Upper Stage Tank Thermodynamic Modeling Using SINDA/FLUINT
! 42nd AIAA/ASME/SAE/ASEE Joint Propulsion Conference & Exhibit
! 9 - 12 July 2006, Sacramento, California
! AIAA 2006-5051
! Paul Schallhorn, D. Michael Campbell, Sukhdeep Chase,
! Jorge Piquero, Cindy Fortenberry, Xiaoyi Li, Lisa Grob
!
! This routine only called when "law_of_the_wall==1" associated with 
! im_fluid.
subroutine wallfunc_thermocorrelation( &
  dir, & ! =1,2,3
  data_dir, & ! =0,1,2
  dxmin, &
  x_projection_raster, &
  dx, &
  n_raster, & ! points to solid
  u, & !intent(in) uimage_raster_solid_frame(dir)
  uimage_tngt_mag, & !intent(in) 
  wall_model_velocity, & ! intent(in)
  dist_probe, & ! intent(in)
  dist_fluid, & ! intent(in)
  temperature_image, & !intent(in) 
  temperature_wall, & ! intent(in)      
  temperature_wall_max, & ! intent(in)      
  viscosity_molecular, & ! intent(in)      
  viscosity_eddy_wall, & ! intent(in)      
  y, & !intent(in) distance from image to wall
  ughost_tngt, & ! intent(out)
  im_fluid, &  ! intent(in)
  critical_length) ! intent(in) used for sanity check
use probcommon_module
implicit none
INTEGER_T, intent(in) :: dir ! 1,2,3
INTEGER_T, intent(in) :: data_dir ! 0,1,2
REAL_T, intent(in) :: dxmin
REAL_T, intent(in), pointer :: x_projection_raster(:)
REAL_T, intent(in), pointer :: dx(:)
REAL_T, intent(in), pointer :: n_raster(:) ! points to solid
INTEGER_T, intent(in) :: im_fluid
REAL_T, intent(in) :: u !uimage_raster_solid_frame(dir)
REAL_T, intent(in) :: uimage_tngt_mag
REAL_T, intent(in) :: wall_model_velocity
REAL_T, intent(in) :: dist_probe
REAL_T, intent(in) :: dist_fluid
REAL_T, intent(in) :: temperature_image
REAL_T, intent(in) :: temperature_wall
REAL_T, intent(in) :: temperature_wall_max
REAL_T, intent(in) :: viscosity_molecular
REAL_T, intent(in) :: viscosity_eddy_wall
REAL_T, intent(in) :: y !delta_r
REAL_T, intent(in) :: critical_length
REAL_T, intent(out) :: ughost_tngt  ! dir direction

REAL_T :: rho_w !wall density
REAL_T :: mu_w  !mu_w: wall molecular viscosity
REAL_T :: thermal_conductivity
REAL_T :: thermal_diffusivity
REAL_T :: gravity_local
REAL_T :: expansion_coefficient
REAL_T :: Cp
REAL_T :: Jtemp,Jtemp_no_area,dtemp,vtemp
REAL_T,parameter :: local_pi=4.0*atan(1.0d0)
INTEGER_T :: turb_flag  ! 1 for turb  0 for laminar
REAL_T :: Ra,Gr,Pr
REAL_T :: nu
REAL_T :: xi
REAL_T :: R
REAL_T :: macro_scale_thickness

if ((im_fluid.lt.1).or.(im_fluid.gt.num_materials)) then
 print *,"im_fluid invalid in wallfunc_thermocorrelation"
 stop
endif

mu_w=fort_viscconst(im_fluid) 
rho_w=fort_denconst(im_fluid)

if (dx(1).gt.0.0d0) then
 ! do nothing
else
 print *,"dx(1) invalid"
 stop
endif

if (y.gt.zero) then
 ! do nothing
else
 print *,"y should be positive"
 stop
endif
if (dxmin.gt.zero) then
 ! do nothing
else
 print *,"dxmin should be positive"
 stop
endif

if (rho_w.gt.zero) then
 ! do nothing
else
 print *,"rho_w invalid"
 stop
endif

if (viscosity_eddy_wall.eq.fort_viscconst_eddy_wall(im_fluid)) then
 ! do nothing
else
 print *,"viscosity_eddy_wall.eq.fort_viscconst_eddy_wall(im_fluid) == false"
 stop
endif

if (mu_w.eq.viscosity_molecular) then
 ! do nothing
else
 print *,"mu_w.eq.viscosity_molecular == false"
 stop
endif

nu=mu_w/rho_w

thermal_conductivity=fort_heatviscconst(im_fluid)
Cp=fort_stiffCP(im_fluid)

if (thermal_conductivity.gt.zero) then
 ! do nothing
else
 print *,"expecting thermal_conductivity to be positive"
 stop
endif

if (Cp.gt.zero) then
 ! do nothing
else
 print *,"Cp invalid"
 stop
endif

xi=x_projection_raster(SDIM)-TANK_MK_HEATER_LOW
R=x_projection_raster(1)

if (1.eq.0) then
 print *,"xi=",xi
 print *,"n_raster(1)=",n_raster(1)
 print *,"temperature_wall ",temperature_wall
 print *,"temperature_wall_max ",temperature_wall_max
 print *,"temperature_image ",temperature_image
endif

if (dir.eq.data_dir+1) then
 if (abs(n_raster(dir)).eq.one) then
  ! do nothing
 else
  print *,"abs(n_raster(dir)) invalid"
  stop
 endif
else if (dir.ne.data_dir+1) then
 if (n_raster(dir).eq.zero) then
  ! do nothing
 else
  print *,"n_raster(dir) invalid"
  stop
 endif
else
 print *,"dir or data_dir corrupt"
 stop
endif

if ((xi.gt.0.0d0).and. &
    (n_raster(1).eq.one)) then

 thermal_diffusivity=thermal_conductivity/(rho_w*Cp)
 gravity_local=abs(gravity)
 expansion_coefficient=abs(fort_DrhoDT(im_fluid)) ! units: 1/temperature
 Pr=nu/thermal_diffusivity

 if (temperature_wall_max.gt.temperature_image) then

  Gr=gravity_local*expansion_coefficient* &
    (temperature_wall_max-temperature_image)*(xi**3.0)/(nu*nu)
  Ra=Gr*Pr

  if(Ra.lt.1.0e+9)then
   turb_flag=0
  else if(Ra.ge.1.0e+9)then
   turb_flag=1
  else
   print *,"Ra number invalid"
   stop
  endif

  vtemp=1.185d0*(nu/xi)*(Gr/(1.0d0+0.494*(Pr**(2.0d0/3.0d0))))**0.5d0

   ! Jtemp is mass flux through horizontal face of length dtemp.
   ! Jtemp has units of mass/sec=rho * velocity * area_face
   ! We assume that J=rho * u_tan * 2 pi R * dr
   ! u_tan=J/(2 pi R rho dr)
  if(turb_flag.eq.1)then
   dtemp=xi*0.565*((1.0d0+0.494d0*Pr**(2.0d0/3.0d0))/Gr)**(0.1d0)/ &
        (Pr**(8.0d0/15.0d0))
   Jtemp_no_area=0.1436d0*rho_w*vtemp
  else if(turb_flag.eq.0)then
   dtemp=xi*3.93*((0.952+Pr)/(Gr*(Pr**2.0d0)))**0.25d0
   Jtemp_no_area=0.0833d0*rho_w*vtemp
  else
   print *,"wrong turb_flag"
   stop
  endif
  Jtemp=2.0d0*local_pi*R*Jtemp_no_area*dtemp

  ! it is known that converged solutions can be obtained on a 128x512 grid.
  ! on a 16x64 grid, choose a thickness associated to the finer grid.
  macro_scale_thickness=dx(1)/16.0d0
  if ((macro_scale_thickness.lt.dtemp).or.(1.eq.1)) then
   macro_scale_thickness=dtemp
  endif
! ughost_tngt=(Jtemp_no_area/rho_w)*dtemp/macro_scale_thickness
  ughost_tngt=vtemp ! vtemp is the maximum tangential velocity in the
                    ! boundary layer region.

 else if (temperature_wall_max.le.temperature_image) then
  Gr=zero
  Ra=zero
  turb_flag=0
  vtemp=zero
  dtemp=zero
  Jtemp_no_area=zero
  Jtemp=zero
  macro_scale_thickness=zero
  ughost_tngt=zero
 else
  print *,"temperature_wall_max or temperature_image is NaN"
  stop
 endif

 if (wall_model_velocity.eq.zero) then
  ! do nothing
 else if (wall_model_velocity.ne.zero) then
  ughost_tngt=wall_model_velocity
 else
  print *,"wall_model_velocity is NaN"
  stop
 endif

  ! do not prescribe a wall velocity if near (or in) another fluid.
  !
 if ((dist_probe.lt.dx(SDIM)).or. &
     (dist_fluid.lt.dx(SDIM))) then
  ughost_tngt=0.0d0
 else if ((dist_probe.ge.dx(SDIM)).and. &
          (dist_fluid.ge.dx(SDIM))) then
  ! do nothing
 else
  print *,"dist_probe or dist_fluid is NaN"
  stop
 endif

 if (1.eq.0) then
  print *,"xi=",xi
  print *,"Gr,Pr,Ra,vtemp,dtemp ",Gr,Pr,Ra,vtemp,dtemp
  print *,"Jtemp=",Jtemp
  print *,"R=",R
  print *,"rho_w=",rho_w
  print *,"dx(1)=",dx(1)
  print *,"ughost_tngt=",ughost_tngt
 endif

else if ((xi.le.0.0d0).or. &
         (n_raster(1).ne.one)) then
 ughost_tngt=0.0d0
else
 print *,"xi,n_raster,twall or timage became corrupt"
 stop
endif

end subroutine wallfunc_thermocorrelation


subroutine CRYOGENIC_TANK_MK_wallfunc( &
  dir, & ! =1,2,3
  data_dir, & ! =0,1,2
  dxmin, &
  x_projection_raster, &
  dx, &
  n_raster, & ! points to solid
  u, & !intent(in) uimage_raster_solid_frame(dir)
  uimage_tngt_mag, & !intent(in) 
  wall_model_velocity, & ! intent(in)
  dist_probe, & ! intent(in)
  dist_fluid, & ! intent(in)
  temperature_image, & !intent(in) 
  temperature_wall, & ! intent(in)      
  temperature_wall_max, & ! intent(in)      
  viscosity_molecular, & ! intent(in)      
  viscosity_eddy_wall, & ! intent(in)      
  y, & !intent(in) distance from image to wall
  ughost_tngt, & ! intent(out)
  im_fluid, &  ! intent(in)
  critical_length) ! intent(in) used for sanity check
use probcommon_module
use global_utility_module
implicit none
INTEGER_T, intent(in) :: dir ! 1,2,3
INTEGER_T, intent(in) :: data_dir ! 0,1,2
REAL_T, intent(in) :: dxmin
REAL_T, intent(in), pointer :: x_projection_raster(:)
REAL_T, intent(in), pointer :: dx(:)
REAL_T, intent(in), pointer :: n_raster(:) ! points to solid
INTEGER_T, intent(in) :: im_fluid
REAL_T, intent(in) :: u !uimage_raster_solid_frame(dir)
REAL_T, intent(in) :: uimage_tngt_mag
REAL_T, intent(in) :: wall_model_velocity
REAL_T, intent(in) :: dist_probe
REAL_T, intent(in) :: dist_fluid
REAL_T, intent(in) :: temperature_image
REAL_T, intent(in) :: temperature_wall
REAL_T, intent(in) :: temperature_wall_max
REAL_T, intent(in) :: viscosity_molecular
REAL_T, intent(in) :: viscosity_eddy_wall
REAL_T, intent(in) :: y !delta_r
REAL_T, intent(in) :: critical_length
REAL_T, intent(out) :: ughost_tngt  ! dir direction

 if (1.eq.0) then
  ! remark: "subroutine wallfunc_newtonsmethod" is 
  ! declared in GLOBALUTIL.F90
  call wallfunc_newtonsmethod( &
   dir, & ! =1,2,3
   data_dir, & ! =0,1,2
   dxmin, &
   x_projection_raster, &
   dx, &
   n_raster, & ! points to solid
   u, & !intent(in) uimage_raster_solid_frame(dir)
   uimage_tngt_mag, & !intent(in) 
   wall_model_velocity, & ! intent(in)
   dist_probe, & ! intent(in)
   dist_fluid, & ! intent(in)
   temperature_image, & !intent(in) 
   temperature_wall, & ! intent(in)      
   temperature_wall_max, & ! intent(in)      
   viscosity_molecular, & ! intent(in)      
   viscosity_eddy_wall, & ! intent(in)      
   y, & !intent(in) distance from image to wall
   ughost_tngt, & ! intent(out)
   im_fluid, &  ! intent(in)
   critical_length) ! intent(in) used for sanity check
 else
  call wallfunc_thermocorrelation( &
   dir, & ! =1,2,3
   data_dir, & ! =0,1,2
   dxmin, &
   x_projection_raster, &
   dx, &
   n_raster, & ! points to solid
   u, & !intent(in) uimage_raster_solid_frame(dir)
   uimage_tngt_mag, & !intent(in) 
   wall_model_velocity, & ! intent(in)
   dist_probe, & ! intent(in)
   dist_fluid, & ! intent(in)
   temperature_image, & !intent(in) 
   temperature_wall, & ! intent(in)      
   temperature_wall_max, & ! intent(in)      
   viscosity_molecular, & ! intent(in)      
   viscosity_eddy_wall, & ! intent(in)      
   y, & !intent(in) distance from image to wall
   ughost_tngt, & ! intent(out)
   im_fluid, &  ! intent(in)
   critical_length) ! intent(in) used for sanity check
 endif

end subroutine CRYOGENIC_TANK_MK_wallfunc

subroutine CRYOGENIC_TANK_MK_K_EFFECTIVE( &
  interface_mass_transfer_model, &
  ireverse, &
  iten, &        
  molar_mass, & ! index: 1..nmat
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
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: interface_mass_transfer_model
INTEGER_T, intent(in) :: num_materials_in
INTEGER_T, intent(in) :: num_species_var_in
INTEGER_T, intent(in) :: ireverse
INTEGER_T, intent(in) :: iten
REAL_T, intent(in) :: molar_mass(num_materials_in)
REAL_T, intent(in) :: species_molar_mass(num_species_var_in)
REAL_T, intent(in) :: k_model_predict(2) ! src,dst
REAL_T, intent(inout) :: k_model_correct(2) ! src,dst
REAL_T, intent(in) :: k_physical_base(2) ! src, dst
REAL_T, intent(in) :: T_probe_src
REAL_T, intent(in) :: T_probe_dst
REAL_T, intent(in) :: LL
REAL_T, intent(in) :: dxprobe_src
REAL_T, intent(in) :: dxprobe_dst

REAL_T :: RA1,RA2

if (interface_mass_transfer_model.eq.1) then

 RA1=(1.7069d+9)*abs(302.0d0-T_probe_src)
 RA2=(5.8395e+8)*abs(302.0d0-T_probe_dst)

 if (RA1.gt.1.0d+4.and.RA1.le.1.0d+7)then
  k_model_correct(1)=dxprobe_src*k_model_predict(1)/0.1016d0* &
           0.54d0*(RA1**(1.0d0/4.0d0))
 else if (RA1.ge.1.0d+7.and.RA1.le.1.0d+11)then
  k_model_correct(1)=dxprobe_src*k_model_predict(1)/0.1016d0* &
           0.15d0*(RA1**(1.0d0/3.0d0))
 else
   print *,"invalid Ra1 number"
   stop
 endif

 if (RA2.gt.1.0d+4.and.RA2.le.1.0d+7)then
  k_model_correct(2)=dxprobe_dst*k_model_predict(2)/0.1016d0* &
           0.54d0*(RA2**(1.0d0/4.0d0))
 else if (RA2.ge.1.0d+7.and.RA2.le.1.0d+11)then
  k_model_correct(2)=dxprobe_dst*k_model_predict(2)/0.1016d0* &
           0.15d0*(RA2**(1.0d0/3.0d0))
 else
   print *,"invalid Ra2 number"
   stop
 endif

else if (interface_mass_transfer_model.eq.0) then
 ! do nothing
else
 print *,"interface_mass_transfer_model invalid"
 stop
endif

end subroutine CRYOGENIC_TANK_MK_K_EFFECTIVE


end module CRYOGENIC_TANK_MK_module
