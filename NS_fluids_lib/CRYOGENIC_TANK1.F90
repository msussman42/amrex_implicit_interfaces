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

! probtype==421 (see run2d/inputs.CRYOGENIC_TANK1)
module CRYOGENIC_TANK1_module

implicit none                   
! Tank outter radius
REAL_T :: TANK1_RADIUS         
! Tank outher height
REAL_T :: TANK1_HEIGHT         
! Tank wall thickness
REAL_T :: TANK1_THICKNESS      
! Location of liquid-gas interface in respect to z=0
REAL_T :: TANK1_LIQUID_HEIGHT  
! Initial vapor pressure
REAL_T :: TANK1_INITIAL_VAPOR_PRESSURE

! Vapor thermodygamma =c_p/c_v
REAL_T :: TANK1_VAPOR_GAMMA
REAL_T :: TANK1_VAPOR_CP
REAL_T :: TANK1_VAPOR_CV
contains

 ! do any initial preparation needed
 subroutine INIT_CRYOGENIC_TANK1_MODULE()
  use probcommon_module
  implicit none
  TANK1_RADIUS = xblob
  TANK1_HEIGHT = yblob
  TANK1_THICKNESS = radblob
  TANK1_LIQUID_HEIGHT = zblob

  TANK1_VAPOR_GAMMA =  1.666666667D0
  TANK1_VAPOR_CV = 6.490D3 ! [J∕(kg·K)]
  TANK1_VAPOR_CP =  TANK1_VAPOR_CV * TANK1_VAPOR_GAMMA ! [J∕(kg·K)]
  TANK1_INITIAL_VAPOR_PRESSURE = fort_denconst(2)*(TANK1_VAPOR_GAMMA-one)*fort_tempconst(2)
  
  return
 end subroutine INIT_CRYOGENIC_TANK1_MODULE


 ! fluids tessellate the domain, solids are immersed. 
 ! fluid interfaces are extended into solids.
 ! material 1 is liquid
 ! material 2 is gas
 ! material 3 is solid

 subroutine CRYOGENIC_TANK1_LS(x,t,LS)
  use probcommon_module
  IMPLICIT NONE

  REAL_T, intent(in) :: x(SDIM)
  REAL_T, intent(in) :: t
  REAL_T, intent(out) :: LS(num_materials)
  REAL_T ls_o,ls_i

  if ((num_materials.eq.3).and.(probtype.eq.421)) then
   ! liquid
   LS(1)=TANK1_LIQUID_HEIGHT-x(2)
   LS(2)=-LS(1)

   ! Solid
   ls_o = DIST_FINITE_CYLINDER(x,TANK1_RADIUS,zero,TANK1_HEIGHT)
   ls_i = DIST_FINITE_CYLINDER(x,TANK1_RADIUS-TANK1_THICKNESS,&
    TANK1_THICKNESS,TANK1_HEIGHT-TANK1_THICKNESS)
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
 end subroutine CRYOGENIC_TANK1_LS

 ! if SOLID VELOCITY requested everywhere (including outside of the solid),
 ! then velsolid==1
 subroutine CRYOGENIC_TANK1_VEL(x,t,LS,VEL,velsolid_flag,dx)
  use probcommon_module
  IMPLICIT NONE

  REAL_T, intent(in) :: x(SDIM)
  REAL_T, intent(in) :: t
  REAL_T, intent(in) :: dx(SDIM)
  REAL_T, intent(in) :: LS(num_materials)
  REAL_T, intent(out) :: VEL(SDIM)
  INTEGER_T, intent(in) :: velsolid_flag
  INTEGER_T dir

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
 end subroutine CRYOGENIC_TANK1_VEL

REAL_T function DIST_FINITE_CYLINDER(P,R_cyl,H_bot,H_top)
 ! Returns the signed distance function to the cylinder
 ! surfaces (including top and bottom)
 ! The axis of cylinder is along SDIM=2 direction
 ! Cylinder radus is R_cyl
 ! Bottom and top faces are at H_bot and H_top
 ! Inside the cylinder < 0
 ! Outside the cylinder > 0
 implicit none

 REAL_T, intent(in), dimension(SDIM) :: P
 REAL_T, intent(in) :: R_cyl
 REAL_T, intent(in) :: H_bot
 REAL_T, intent(in) :: H_top
 
 REAL_T x,y,z,r
 REAL_T dist_cyl, dist_end
 
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
subroutine EOS_CRYOGENIC_TANK1(rho,internal_energy,pressure)
 IMPLICIT NONE
 REAL_T, intent(in) :: rho
 REAL_T, intent(in) :: internal_energy
 REAL_T, intent(out) :: pressure

 pressure=rho*(TANK1_VAPOR_GAMMA-1)*internal_energy

 return
end subroutine EOS_CRYOGENIC_TANK1

subroutine SOUNDSQR_CRYOGENIC_TANK1(rho,internal_energy,soundsqr)
 IMPLICIT NONE
 REAL_T, intent(in) :: rho
 REAL_T, intent(in) :: internal_energy
 REAL_T, intent(out) :: soundsqr
 REAL_T pressure

 call EOS_CRYOGENIC_TANK1(rho,internal_energy,pressure)
 soundsqr=TANK1_VAPOR_GAMMA*pressure/rho

 return
end subroutine SOUNDSQR_CRYOGENIC_TANK1

subroutine INTERNAL_CRYOGENIC_TANK1(rho,temperature,local_internal_energy)
 IMPLICIT NONE
 REAL_T, intent(in) :: rho
 REAL_T, intent(in) :: temperature 
 REAL_T, intent(out) :: local_internal_energy
 
 local_internal_energy=TANK1_VAPOR_CV*temperature

 return
end subroutine INTERNAL_CRYOGENIC_TANK1

subroutine TEMPERATURE_CRYOGENIC_TANK1(rho,temperature,internal_energy)
 IMPLICIT NONE
 REAL_T, intent(in) :: rho
 REAL_T, intent(out) :: temperature 
 REAL_T, intent(in) :: internal_energy

 temperature=internal_energy/TANK1_VAPOR_CV

 return
end subroutine TEMPERATURE_CRYOGENIC_TANK1

!***********************************************
! called by the boundary condition routine
! might be called at initialization, so put a placeholder pressure here.
subroutine CRYOGENIC_TANK1_PRES(x,t,LS,PRES)
use probcommon_module
IMPLICIT NONE

REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(num_materials)
REAL_T, intent(out) :: PRES

PRES=TANK1_INITIAL_VAPOR_PRESSURE 

return 
end subroutine CRYOGENIC_TANK1_PRES


subroutine CRYOGENIC_TANK1_STATE(x,t,LS,STATE)
use probcommon_module
IMPLICIT NONE

REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(num_materials)
REAL_T, intent(out) :: STATE(num_materials*num_state_material)
INTEGER_T im,ibase,n

 ! num_state_material=2 (default)  density and temperature
 ! num_state_material>2 if scalar (species) variables added.
if ((num_materials.eq.3).and. &
    (num_state_material.ge.2).and. &
    (probtype.eq.421)) then
 do im=1,num_materials
  ibase=(im-1)*num_state_material
  STATE(ibase+1)=fort_denconst(im)
  if (t.eq.zero) then
   STATE(ibase+2)=fort_initial_temperature(im)
  else if (t.gt.zero) then
   STATE(ibase+2)=fort_tempconst(im)
  else
   print *,"t invalid"
   stop
  endif
  do n=1,num_species_var
   STATE(ibase+2+n)=fort_speciesconst((n-1)*num_materials+im)
  enddo
 enddo ! im=1..num_materials
else
 print *,"num_materials,num_state_material, or probtype invalid"
 stop
endif
 
return
end subroutine CRYOGENIC_TANK1_STATE

 ! dir=1..sdim  side=1..2
subroutine CRYOGENIC_TANK1_LS_BC(xwall,xghost,t,LS, &
   LS_in,dir,side,dx)
use probcommon_module
IMPLICIT NONE

REAL_T xwall
REAL_T xghost(SDIM)
REAL_T t
REAL_T LS(num_materials)
REAL_T LS_in(num_materials)
INTEGER_T dir,side
REAL_T dx(SDIM)

if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then
 call CRYOGENIC_TANK1_LS(xghost,t,LS)
else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine CRYOGENIC_TANK1_LS_BC


 ! dir=1..sdim  side=1..2 veldir=1..sdim
subroutine CRYOGENIC_TANK1_VEL_BC(xwall,xghost,t,LS, &
   VEL,VEL_in,veldir,dir,side,dx)
use probcommon_module
IMPLICIT NONE

REAL_T xwall
REAL_T xghost(SDIM)
REAL_T t
REAL_T LS(num_materials)
REAL_T VEL
REAL_T VEL_in
INTEGER_T veldir,dir,side
REAL_T dx(SDIM)
REAL_T local_VEL(SDIM)
INTEGER_T velsolid_flag

velsolid_flag=0
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2).and. &
    (veldir.ge.1).and.(veldir.le.SDIM)) then

 call CRYOGENIC_TANK1_VEL(xghost,t,LS,local_VEL,velsolid_flag,dx)
 VEL=local_VEL(veldir)

else
 print *,"dir,side, or veldir invalid"
 stop
endif

return
end subroutine CRYOGENIC_TANK1_VEL_BC

 ! dir=1..sdim  side=1..2
subroutine CRYOGENIC_TANK1_PRES_BC(xwall,xghost,t,LS, &
   PRES,PRES_in,dir,side,dx)
use probcommon_module
IMPLICIT NONE

REAL_T xwall
REAL_T xghost(SDIM)
REAL_T t
REAL_T LS(num_materials)
REAL_T PRES
REAL_T PRES_in
INTEGER_T dir,side
REAL_T dx(SDIM)

if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then

 call CRYOGENIC_TANK1_PRES(xghost,t,LS,PRES)

else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine CRYOGENIC_TANK1_PRES_BC

 ! dir=1..sdim  side=1..2
subroutine CRYOGENIC_TANK1_STATE_BC(xwall,xghost,t,LS, &
   STATE,STATE_merge,STATE_in,im,istate,dir,side,dx)
use probcommon_module
IMPLICIT NONE

REAL_T xwall
REAL_T xghost(SDIM)
REAL_T t
REAL_T LS(num_materials)
REAL_T local_STATE(num_materials*num_state_material)
REAL_T STATE
REAL_T STATE_merge
REAL_T STATE_in
INTEGER_T dir,side
REAL_T dx(SDIM)
INTEGER_T istate,im
INTEGER_T ibase,im_crit,im_loop

if ((istate.ge.1).and. &
    (istate.le.num_state_material).and. &
    (im.ge.1).and. &
    (im.le.num_materials)) then
 call CRYOGENIC_TANK1_STATE(xghost,t,LS,local_STATE)
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
end subroutine CRYOGENIC_TANK1_STATE_BC

! suppose inhomogeneous flux condition: k grad T dot n= q
! n outward facing normal
! 1. set "k"=0 on the boundary
! 2. T_t + div k grad T = F
! 3. T_t + ((k grad T)_right - (k grad T)_left)/dx = F 
! 4. at the left wall:
!    T_t + (k grad T)_right/dx = F - q/dx 
subroutine CRYOGENIC_TANK1_HEATSOURCE(im,VFRAC,time,x,temp, &
     heat_source,den,CV,dt)
use probcommon_module
IMPLICIT NONE

INTEGER_T im
REAL_T VFRAC(num_materials)
REAL_T time
REAL_T x(SDIM)
REAL_T, intent(in) :: temp(num_materials)
REAL_T den(num_materials)
REAL_T CV(num_materials)
REAL_T dt
REAL_T heat_source

if ((num_materials.eq.3).and.(probtype.eq.421)) then
 heat_source=zero
else
 print *,"num_materials ", num_materials
 print *,"probtype ", probtype
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine CRYOGENIC_TANK1_HEATSOURCE

end module CRYOGENIC_TANK1_module
