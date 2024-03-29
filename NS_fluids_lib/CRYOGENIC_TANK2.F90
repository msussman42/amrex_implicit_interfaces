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


#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else
print *,"dimension bust"
stop
#endif

! probtype==422 (see run2d/inputs.CRYOGENIC_TANK2)
module CRYOGENIC_TANK2_module
use amrex_fort_module, only : amrex_real

implicit none                   
! Tank outter radius
real(amrex_real) :: TANK2_RADIUS         
! Tank outher height
real(amrex_real) :: TANK2_HEIGHT         
! Tank wall thickness
real(amrex_real) :: TANK2_THICKNESS      
! Location of liquid-gas interface in respect to z=0
real(amrex_real) :: TANK2_LIQUID_HEIGHT  

! Initial mixture pressure
real(amrex_real) :: TANK2_INITIAL_PRESSURE
! Universal gas constant [J/(mol K)]
real(amrex_real) :: TANK2_R_UNIV

real(amrex_real) :: TANK2_GAS_GAMMA
real(amrex_real) :: TANK2_GAS_CP
real(amrex_real) :: TANK2_GAS_CV

real(amrex_real) :: TANK2_HEATER_FLUX
! Heater location in dim=2 direction
real(amrex_real) :: TANK2_HEATER_LOW
real(amrex_real) :: TANK2_HEATER_HIGH
real(amrex_real) :: TANK2_HEATER_R

contains

 ! do any initial preparation needed
 subroutine INIT_CRYOGENIC_TANK2_MODULE()
  use probcommon_module
  implicit none
  TANK2_RADIUS = xblob
  TANK2_HEIGHT = yblob
  TANK2_THICKNESS = radblob
  TANK2_LIQUID_HEIGHT = zblob

  TANK2_HEATER_FLUX      = xblob3
  TANK2_HEATER_LOW       = yblob3
  TANK2_HEATER_HIGH      = zblob3
  TANK2_HEATER_R         = radblob3

  ! ASSUMING IDEA GAS => The gas heat cpacities should satisfy this
  ! R_spc = C_{p,spc}-C_{v,spc}
  ! to have ideal mixture gas as well =>
  ! Only C_p or C_v can be picked from table and the other one
  ! calculated from equation above.
  ! Here we pick C_{v,spc} from input file.
  ! C_{p,spc} = C_{v,spc} + R_spc

  TANK2_R_UNIV = 8.31446261815324D0

  TANK2_GAS_CV = fort_stiffCV(2) ![J/(kg K)]
!  TANK2_GAS_CP = fort_stiffCP(2) ![J∕(kg·K)]
  TANK2_GAS_CP = TANK2_GAS_CV + TANK2_R_UNIV/fort_molar_mass(2)  ! [J∕(kg·K)]
  TANK2_GAS_GAMMA = TANK2_GAS_CP / TANK2_GAS_CV

  ! Initial pressure based on the given density and pressure
  ! P = rho R_sp T = rho (gamma-1) U
  TANK2_INITIAL_PRESSURE = &
   fort_denconst(2)*(TANK2_GAS_GAMMA-one)*&
   fort_initial_temperature(2)*TANK2_GAS_CV 
  
  return
 end subroutine INIT_CRYOGENIC_TANK2_MODULE


 ! fluids tessellate the domain, solids are immersed. 
 ! fluid interfaces are extended into solids.
 ! material 1 is liquid
 ! material 2 is gas
 ! material 3 is solid

 subroutine CRYOGENIC_TANK2_LS(x,t,LS,nmat)
  use probcommon_module
  IMPLICIT NONE

  integer, INTENT(in) :: nmat
  real(amrex_real), INTENT(in) :: x(SDIM)
  real(amrex_real), INTENT(in) :: t
  real(amrex_real), INTENT(out) :: LS(nmat)
  real(amrex_real) ls_o,ls_i
  real(amrex_real), PARAMETER :: stub_zero=zero

  if (nmat.eq.num_materials) then
   ! do nothing
  else
   print *,"nmat invalid"
   stop
  endif

  if ((num_materials.eq.3).and.(probtype.eq.422)) then
   ! liquid
   LS(1)=TANK2_LIQUID_HEIGHT-x(2)

   if (radblob2.eq.zero) then
    ! do nothing
   else if (radblob2.gt.zero) then
    if (SDIM.eq.2) then
     LS(1)=sqrt((x(1)-xblob2)**2+(x(2)-yblob2)**2)-radblob2
    else if (SDIM.eq.3) then 
     LS(1)=sqrt((x(1)-xblob2)**2+(x(2)-yblob2)**2+ &
                (x(SDIM)-zblob2)**2)-radblob2
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
   ls_o = DIST_FINITE_CYLINDER(x,TANK2_RADIUS,stub_zero,TANK2_HEIGHT)
   ls_i = DIST_FINITE_CYLINDER(x,TANK2_RADIUS-TANK2_THICKNESS,&
    TANK2_THICKNESS,TANK2_HEIGHT-TANK2_THICKNESS)
   if((ls_o.ge.zero).and.(ls_i.ge.zero)) then
    ! outside of tank
    LS(3) = -min(ls_o,ls_i)
   else if((ls_o.le.zero).and.(ls_i.le.zero)) then
    ! inside of tank cavity
    LS(3) = -min(-ls_o,-ls_i)
   else if((ls_o.lt.zero).and.(ls_i.gt.zero)) then
    ! inside of tank wall
    LS(3) = min(-ls_o,ls_i)
   else
    print *,"tank level set calculation failed!"
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
 end subroutine CRYOGENIC_TANK2_LS

 ! if SOLID VELOCITY requested everywhere (including outside of the solid),
 ! then velsolid==1
 subroutine CRYOGENIC_TANK2_VEL(x,t,LS,VEL,velsolid_flag,dx,nmat)
  use probcommon_module
  IMPLICIT NONE

  integer, INTENT(in) :: nmat
  real(amrex_real), INTENT(in) :: x(SDIM)
  real(amrex_real), INTENT(in) :: t
  real(amrex_real), INTENT(in) :: dx(SDIM)
  real(amrex_real), INTENT(in) :: LS(nmat)
  real(amrex_real), INTENT(out) :: VEL(SDIM)
  integer, INTENT(in) :: velsolid_flag
  integer dir

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
 end subroutine CRYOGENIC_TANK2_VEL

real(amrex_real) function DIST_FINITE_CYLINDER(P,R_cyl,H_bot,H_top)
 ! Returns the signed distance function to the cylinder
 ! surfaces (including top and bottom)
 ! The axis of cylinder is along SDIM=2 direction
 ! Cylinder radus is R_cyl
 ! Bottom and top faces are at H_bot and H_top
 ! Inside the cylinder < 0
 ! Outside the cylinder > 0
 implicit none

 real(amrex_real), INTENT(in), dimension(SDIM) :: P
 real(amrex_real), INTENT(in) :: R_cyl
 real(amrex_real), INTENT(in) :: H_bot
 real(amrex_real), INTENT(in) :: H_top
 
 real(amrex_real) x,y,z,r
 real(amrex_real) dist_cyl, dist_end
 
 x=P(1)
 y=P(2)
 if (SDIM.eq.2) then
  z=zero
 else if(SDIM.eq.3) then
  z=P(SDIM)
 else
  print *,"dimension bust at DIST_FINITE_CYLINDER"
 endif

 r = sqrt(x**2+z**2)
 if((H_bot.le.y).and.(y.le.H_top)) then
  ! between top and bottom
  if(r.ge.R_cyl) then
   ! outside
   DIST_FINITE_CYLINDER = r-R_cyl
  else if (r.lt.R_cyl) then
   ! inside
   dist_cyl = R_cyl-r
   dist_end = min(H_top-y,y-H_bot)
   DIST_FINITE_CYLINDER = -min(dist_cyl,dist_end)
  else
   print *,"r=",r
   print *,"invalid r value at DIST_FINITE_CYLINDER (1)"
   stop
  endif

 else if (y.gt.H_top) then
  ! higher than top
  if(r.le.R_cyl) then
   ! inside infinite cylinder
   DIST_FINITE_CYLINDER = y-H_top
  else if (r.gt.R_cyl) then
   ! outside infinite cylinder
   ! distance to the edge of the top
   DIST_FINITE_CYLINDER = &
    sqrt((r-R_cyl)**2 + (y-H_top)**2)
  else
   print *,"r=",r
   print *,"invalid r value at DIST_FINITE_CYLINDER (2)"
   stop
  endif

 else if (y.lt.H_bot) then
  ! lower than bottom
  if(r.le.R_cyl) then
   ! inside infinite cylinder
   DIST_FINITE_CYLINDER = H_bot-y
  else if (r.gt.R_cyl) then
   ! outside infinite cylinder
   ! distance to the edge of the bottom
   DIST_FINITE_CYLINDER = &
    sqrt((r-R_cyl)**2 + (H_bot-y)**2)
  else
   print *,"r=",r
   print *,"invalid r value at DIST_FINITE_CYLINDER (3)"
   stop
  endif
 else
  print *,"invalid y value at DIST_FINITE_CYLINDER"
  stop
 endif
end function DIST_FINITE_CYLINDER

!***********************************************
! compressible material functions for (ns.material_type = 24)
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

subroutine EOS_CRYOGENIC_TANK2(rho,massfrac_var, &
  internal_energy,pressure, &
  imattype,im,num_species_var_in)
 use probcommon_module
 use global_utility_module
 IMPLICIT NONE
 integer, INTENT(in) :: imattype,im,num_species_var_in
 real(amrex_real), INTENT(in) :: rho
 real(amrex_real), INTENT(in) :: massfrac_var(num_species_var_in+1)
 real(amrex_real), INTENT(in) :: internal_energy
 real(amrex_real), INTENT(out) :: pressure

 if (num_species_var_in.eq.num_species_var) then
  if (im.eq.2) then
   if (imattype.eq.24) then
    ! p = rho (gamme-1) U
    pressure=rho * (TANK2_GAS_GAMMA-one) * internal_energy
   else
    print *,"imattype= ",imattype
    print *,"imattype invalid EOS_CRYOGENIC_TANK2"
    stop
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
end subroutine EOS_CRYOGENIC_TANK2


subroutine dVdT_CRYOGENIC_TANK2(dVdT,massfrac_var, &
   pressure,temperature, &
   imattype,im,num_species_var_in)
use probcommon_module
use global_utility_module
IMPLICIT NONE
integer, INTENT(in) :: imattype,im,num_species_var_in
real(amrex_real), INTENT(in) :: pressure
real(amrex_real), INTENT(in) :: massfrac_var(num_species_var_in+1)
real(amrex_real), INTENT(in) :: temperature
real(amrex_real), INTENT(out) :: dVdT
integer :: dummy_input

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
   if (imattype.eq.24) then
    ! p = rho (gamma-1) e  rho=1/V
    ! V = (gamma-1) e/p=(gamma-1)Cv T/p
    ! dVdT=(gamma-1)Cv/p
    dVdT=(TANK2_GAS_GAMMA-one) * TANK2_GAS_CV/pressure
   else
    print *,"imattype= ",imattype
    print *,"imattype invalid dVdT_CRYOGENIC_TANK2"
    print *,"break point and gdb: "
    print *,"(1) compile with the -g option"
    print *,"(2) break CRYOGENIC_TANK2.F90:350"
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
end subroutine dVdT_CRYOGENIC_TANK2


subroutine SOUNDSQR_CRYOGENIC_TANK2(rho,massfrac_var, &
  internal_energy,soundsqr, &
  imattype,im,num_species_var_in)
 use probcommon_module
 use global_utility_module
 IMPLICIT NONE
 integer, INTENT(in) :: imattype,im,num_species_var_in
 real(amrex_real), INTENT(in) :: rho
 real(amrex_real), INTENT(in) :: massfrac_var(num_species_var_in+1)
 real(amrex_real), INTENT(in) :: internal_energy
 real(amrex_real), INTENT(out) :: soundsqr
 real(amrex_real) pressure

 if (num_species_var_in.eq.num_species_var) then
  if (im.eq.2) then
   if (imattype.eq.24) then
     ! a = sqrt(gamma R_sp T) = sqrt(gamma p/rho)
    call EOS_CRYOGENIC_TANK2(rho,massfrac_var, &
     internal_energy,pressure,imattype,im,num_species_var_in)
    if (rho.gt.zero) then
     soundsqr=TANK2_GAS_GAMMA*pressure/rho
    else
     print *,"rho invalid"
     stop
    endif
   else
    print *,"imattype= ",imattype
    print *,"imattype invalid SOUNDSQR CRYOGENIC TANK2"
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
end subroutine SOUNDSQR_CRYOGENIC_TANK2

subroutine INTERNAL_CRYOGENIC_TANK2(rho,massfrac_var, &
  temperature,local_internal_energy, &
  imattype,im,num_species_var_in)
 use probcommon_module
 use global_utility_module
 IMPLICIT NONE
 integer, INTENT(in) :: imattype,im,num_species_var_in
 real(amrex_real), INTENT(in) :: rho
 real(amrex_real), INTENT(in) :: massfrac_var(num_species_var_in+1)
 real(amrex_real), INTENT(in) :: temperature 
 real(amrex_real), INTENT(out) :: local_internal_energy

 if (num_species_var_in.eq.num_species_var) then
  if (im.eq.2) then
   if ((imattype.eq.24).or.(imattype.eq.0)) then 
    ! U_mix = C_{v,spc} T
    local_internal_energy=TANK2_GAS_CV*temperature
   else
    print *,"imattype= ",imattype
    print *,"imattype invalid INTERNAL CRYOGENIC TANK2"
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
end subroutine INTERNAL_CRYOGENIC_TANK2

subroutine TEMPERATURE_CRYOGENIC_TANK2(rho,massfrac_var, &
  temperature,internal_energy, &
  imattype,im,num_species_var_in)
 use probcommon_module
 use global_utility_module
 IMPLICIT NONE
 integer, INTENT(in) :: imattype,im,num_species_var_in
 real(amrex_real), INTENT(in) :: rho
 real(amrex_real), INTENT(in) :: massfrac_var(num_species_var_in+1)
 real(amrex_real), INTENT(out) :: temperature 
 real(amrex_real), INTENT(in) :: internal_energy

 if (num_species_var_in.eq.num_species_var) then
  if (im.eq.2) then
   if ((imattype.eq.24).or.(imattype.eq.0)) then 
    ! T = U / C_{v,spc}
    if (TANK2_GAS_CV.gt.zero) then
     temperature=internal_energy/TANK2_GAS_CV
    else
     print *,"TANK2_GAS_CV invalid 1"
     stop
    endif
   else
    print *,"imattype= ",imattype
    print *,"imattype invalid TEMPERATURE_CRYOGENIC_TANK2"
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
end subroutine TEMPERATURE_CRYOGENIC_TANK2

!***********************************************
! called by the boundary condition routine
! might be called at initialization, so put a placeholder pressure here.
subroutine CRYOGENIC_TANK2_PRES(x,t,LS,PRES,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: LS(nmat)
real(amrex_real), INTENT(out) :: PRES
integer gravity_dir

if (num_materials.eq.nmat) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

!PRES=TANK2_INITIAL_GAS_PRESSURE 

 call fort_derive_gravity_dir(gravity_vector,gravity_dir)

 if (x(2).ge.TANK2_LIQUID_HEIGHT) then
  PRES=TANK2_INITIAL_PRESSURE
 elseif (x(2).lt.TANK2_LIQUID_HEIGHT) then
  PRES=TANK2_INITIAL_PRESSURE +&
    fort_denconst(1)*(TANK2_LIQUID_HEIGHT-x(2))* &
    (abs(gravity_vector(gravity_dir)))
 else
  print *,"x(2) is invalid in CRYOGENIC_TANK2_PRES!"
  stop
 endif

return 
end subroutine CRYOGENIC_TANK2_PRES


subroutine CRYOGENIC_TANK2_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
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
    (probtype.eq.422)) then
 do im=1,num_materials
  ibase=(im-1)*num_state_material
  STATE(ibase+ENUM_DENVAR+1)=fort_denconst(im)
  if (t.eq.zero) then
   STATE(ibase+ENUM_TEMPERATUREVAR+1)=fort_initial_temperature(im)
  else if (t.gt.zero) then
   STATE(ibase+ENUM_TEMPERATUREVAR+1)=fort_tempconst(im)
  else
   print *,"t invalid"
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
end subroutine CRYOGENIC_TANK2_STATE

 ! dir=1..sdim  side=1..2
subroutine CRYOGENIC_TANK2_LS_BC(xwall,xghost,t,LS, &
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
 call CRYOGENIC_TANK2_LS(xghost,t,LS,nmat)
else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine CRYOGENIC_TANK2_LS_BC


 ! dir=1..sdim  side=1..2 veldir=1..sdim
subroutine CRYOGENIC_TANK2_VEL_BC(xwall,xghost,t,LS, &
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

 call CRYOGENIC_TANK2_VEL(xghost,t,LS,local_VEL,velsolid_flag,dx,nmat)
 VEL=local_VEL(veldir)

else
 print *,"dir,side, or veldir invalid"
 stop
endif

return
end subroutine CRYOGENIC_TANK2_VEL_BC

 ! dir=1..sdim  side=1..2
subroutine CRYOGENIC_TANK2_PRES_BC(xwall,xghost,t,LS, &
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

 call CRYOGENIC_TANK2_PRES(xghost,t,LS,PRES,nmat)

else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine CRYOGENIC_TANK2_PRES_BC

 ! dir=1..sdim  side=1..2
subroutine CRYOGENIC_TANK2_STATE_BC(xwall,xghost,t,LS, &
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
 call CRYOGENIC_TANK2_STATE(xghost,t,LS,local_STATE, &
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
end subroutine CRYOGENIC_TANK2_STATE_BC


! suppose inhomogeneous flux condition: -k grad T = q
! 1. T_t - div k grad T = 0
! 3. T_t - (1/V) sum (A_{i} (k grad T)_i dot n_i) =0  n_i=outward facing normal
! 4. at the right wall:
!     (k grad T)_{right} = -q_{right}
! 5. T_{t} - (1/V) sum_{except right} (A_{i} (k grad T)_{i} dot n_{i} =
!      (1/V)A_{right} (-q_{right} dot n_right)
!
! TANK_MK_HEATER_FLUX \equiv -q dot n
! MEHDI VAHAB HEAT SOURCE
! T^new=T^* + dt * (tildeQ)/(rho cv)    
! dt=seconds  rho=kg/m^3   cv=Joules/(kg K)
! second * J/(m^3 s)  * (m^3/kg)  *  (K kg/J) = degrees Kelvin
! tildeQ units: J/(m^3 s)
! Q = k grad T = W/(m K) K/m = W/m^2= J/(m^2 s)
! tildeQ=Q * area/volume
! called from: GODUNOV_3D.F90, subroutine fort_heatsource
! in fort_heatsource:
! T_local(im)=T_local(im)+ &
!   dt*DeDTinverse(D_DECL(i,j,k),1)*heat_source_total  im=1..nmat
!      (1/V)A_{right} (-q_{right} dot n_right)  => tildeQ corresponds
! to q * A/V

subroutine CRYOGENIC_TANK2_HEATSOURCE( &
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

integer dir
real(amrex_real) local_dx(SDIM)
real(amrex_real) flux_magnitude
real(amrex_real) denom

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

if ((num_materials.eq.3).and.(probtype.eq.422)) then
 heat_source=zero
 if (levelrz.eq.COORDSYS_RZ) then
  if ((abs(xsten(-1,1)).le.TANK2_HEATER_R).and.&
   (abs(xsten(1,1)).ge.TANK2_HEATER_R).and.&
   (xsten(1,2).ge.TANK2_HEATER_LOW).and.&
   (xsten(-1,2).le.TANK2_HEATER_HIGH)) then
   ! area=2 pi rf dz
   ! vol =2 pi rc dr dz
   ! area/vol=rf/(rc dr)
   ! input file value in J/(s.m^2) (flux into the face)
   ! Transforming to J/(s.m^3) (flux into the control volume)
   flux_magnitude=TANK2_HEATER_FLUX
   denom=xsten(0,1)*local_dx(1)
   if (denom.gt.zero) then
    flux_magnitude=flux_magnitude*xsten(1,1)/denom
   else
    print *,"denom invalid 3"
    stop
   endif
  endif
 else if (levelrz.eq.COORDSYS_CARTESIAN) then
  if ((abs(xsten(-1,1)).le.TANK2_HEATER_R).and.&
   (abs(xsten(1,1)).ge.(-TANK2_HEATER_R)).and.&
   (xsten(1,2).ge.TANK2_HEATER_LOW).and.&
   (xsten(-1,2).le.TANK2_HEATER_HIGH)) then

   flux_magnitude=TANK2_HEATER_FLUX/local_dx(2)
  endif
 else
  print *,"levelrz invalid"
  stop
 endif
 heat_source=heat_source+flux_magnitude
else
 print *,"num_materials ", num_materials
 print *,"probtype ", probtype
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine CRYOGENIC_TANK2_HEATSOURCE

end module CRYOGENIC_TANK2_module
