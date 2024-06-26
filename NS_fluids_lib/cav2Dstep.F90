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

! probtype==412 (see run2d/inputs.CAVTEST2D_STEP)
module CAV2Dstep_module
use amrex_fort_module, only : amrex_real

implicit none                   

contains

  ! do any initial preparation needed
subroutine INIT_CAV2Dstep_MODULE()
IMPLICIT NONE

return
end subroutine INIT_CAV2Dstep_MODULE

! Phi>0 in the solid
subroutine nozzle2d_step(x,y,Phi) !return nozzle signed distance function Phi
use probcommon_module

 !  1 ___________________________________ topwall
 !     -->                      -->
 !     -->             __________________ stepheight
 !     -->            |
 !     -->            |
 !  0 ________________|                  botwall
 !   -L              0 stepoffset      L
implicit none
real(amrex_real), INTENT(in) :: x, y !spatial coordinates
real(amrex_real), INTENT(out) :: Phi !LS dist, Phi>0 in the solid
 
integer :: insideflag
real(amrex_real) :: topwall, botwall, stepheight, stepoffset
real(amrex_real) :: temp
 
!geometry params, requires: botwall<stepheight<topwall, -L<stepoffset<L
topwall = one !y
botwall = zero !y
stepheight = half !y
stepoffset = zero !x

Phi = zero
insideflag = 0
temp = 0

!!inside
!topwall (horiz)
if (y.ge.topwall) then
 Phi = y-topwall
 insideflag = 1
endif
!botwall,wide (horiz)
if (y.le.botwall) then 
 insideflag = 1
 Phi = botwall-y
endif
!step (vert)
if (y.le.stepheight .and. x.ge.stepoffset) then
 Phi = x-stepoffset
 insideflag = 1
endif
!corner
if (y.le.botwall .and. x.ge.stepoffset) then
 temp = sqrt((y-botwall)**2+(x-stepoffset)**2)
 Phi = MAX(Phi,temp)
 insideflag = 1
endif
!botwall (narrow)
if (y.le.stepheight .and. x.ge.stepoffset) then
 Phi = MIN(Phi, stepheight-y)
 insideflag = 1
endif

if (insideflag==0) then
 !!outside
 !topwall (horiz)
 if (y.lt.topwall) then
  Phi = y-topwall
 endif
 !botwall,wide (horiz)
 if (x.lt.stepoffset .and. y.gt.botwall) then
  temp = botwall-y
 endif
 !botwall,narrow (horiz)
 if (x.ge.stepoffset .and. y.gt.stepheight) then
  temp = stepheight-y
 endif
 Phi = MAX(Phi,temp)
 !step (vert)
 if (y.gt.botwall .and. y.lt.stepheight .and. x.lt.stepoffset) then
  temp = x-stepoffset
  Phi = MAX(Phi,temp)
 endif
 !corner
 if (y.ge.stepheight .and. y.lt.topwall .and. x.le.stepoffset) then
  temp = 0-sqrt((y-stepheight)**2+(x-stepoffset)**2)
  Phi = MAX(Phi,temp)
 endif
endif !end outside

end subroutine nozzle2d_step


 ! fluids tessellate the domain, solids are immersed. 
subroutine CAV2Dstep_LS(x,t,LS)
use probcommon_module
IMPLICIT NONE

real(amrex_real) x(SDIM)
real(amrex_real) t
integer im
real(amrex_real) LS(num_materials)

if ((num_materials.eq.3).and.(probtype.eq.412)) then
 do im=1,num_materials
  if (im.eq.1) then !liquid
   LS(im)=99999.0
  else if (im.eq.2) then ! ambient
   LS(im)=-99999.0
  else if (im.eq.num_materials) then ! geometry 
   call nozzle2d_step(x(1),x(SDIM),LS(im))
  else
   print *,"im invalid"
   stop
  endif
 enddo ! im=1..num_materials
else
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine CAV2Dstep_LS

subroutine CAV2Dstep_VEL(x,t,LS,VEL,velsolid_flag)
use probcommon_module
IMPLICIT NONE

real(amrex_real) x(SDIM)
real(amrex_real) t
real(amrex_real) LS(num_materials)
real(amrex_real) VEL(SDIM)
integer dir
integer velsolid_flag

if ((velsolid_flag.eq.0).or. &
    (velsolid_flag.eq.1)) then
 ! do nothing
else 
 print *,"velsolid_flag invalid"
 stop
endif

! flow enters wide part (yhi) and exits narrow part (ylo)
if ((adv_dir.ge.1).and.(adv_dir.le.SDIM)) then

 if (adv_dir.eq.1) then

  if ((LS(2).ge.zero).or. &
      (velsolid_flag.eq.1)) then
   ! in solid
   do dir=1,SDIM
    VEL(dir)=zero
   enddo
  else if ((LS(2).le.zero).and. &
           (velsolid_flag.eq.0)) then
   VEL(1)=adv_vel
   do dir=2,SDIM
    VEL(dir)=zero
   enddo
  else
   print *,"LS(2) bust"
   stop
  endif

 else
  print *,"expecting adv_dir = 1 in CAV2Dstep_VEL"
  stop
 endif

else
 print *,"adv_dir invalid in CAV2Dstep_VEL"
 stop
endif

return 
end subroutine CAV2Dstep_VEL

subroutine CAV2Dstep_PRES(x,t,LS,PRES)
use probcommon_module
IMPLICIT NONE

real(amrex_real) x(SDIM)
real(amrex_real) t
real(amrex_real) LS(num_materials)
real(amrex_real) PRES

PRES=outflow_pressure

return 
end subroutine CAV2Dstep_PRES


subroutine CAV2Dstep_STATE(x,t,LS,STATE)
use probcommon_module
IMPLICIT NONE

real(amrex_real) x(SDIM)
real(amrex_real) t
real(amrex_real) LS(num_materials)
real(amrex_real) STATE(num_materials*num_state_material)
integer im,ibase,n

if ((num_materials.eq.3).and. &
    (num_state_material.ge.2).and. &
    (probtype.eq.412)) then
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
end subroutine CAV2Dstep_STATE

 ! dir=1..sdim  side=1..2
subroutine CAV2Dstep_LS_BC(xwall,xghost,t,LS, &
   LS_in,dir,side,dx)
use probcommon_module
IMPLICIT NONE

real(amrex_real) xwall
real(amrex_real) xghost(SDIM)
real(amrex_real) t
real(amrex_real) LS(num_materials)
real(amrex_real) LS_in(num_materials)
integer dir,side
real(amrex_real) dx(SDIM)

if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then
 call CAV2Dstep_LS(xghost,t,LS)
else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine CAV2Dstep_LS_BC


 ! dir=1..sdim  side=1..2 veldir=1..sdim
subroutine CAV2Dstep_VEL_BC(xwall,xghost,t,LS, &
   VEL,VEL_in,veldir,dir,side,dx)
use probcommon_module
IMPLICIT NONE

real(amrex_real) xwall
real(amrex_real) xghost(SDIM)
real(amrex_real) t
real(amrex_real) LS(num_materials)
real(amrex_real) VEL
real(amrex_real) VEL_in
integer veldir,dir,side
real(amrex_real) dx(SDIM)
real(amrex_real) local_VEL(SDIM)
integer velsolid_flag

velsolid_flag=0
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2).and. &
    (veldir.ge.1).and.(veldir.le.SDIM)) then

 call CAV2Dstep_VEL(xghost,t,LS,local_VEL,velsolid_flag)
 VEL=local_VEL(veldir)

else
 print *,"dir,side, or veldir invalid"
 stop
endif

return
end subroutine CAV2Dstep_VEL_BC

 ! dir=1..sdim  side=1..2
subroutine CAV2Dstep_PRES_BC(xwall,xghost,t,LS, &
   PRES,PRES_in,dir,side,dx)
use probcommon_module
IMPLICIT NONE

real(amrex_real) xwall
real(amrex_real) xghost(SDIM)
real(amrex_real) t
real(amrex_real) LS(num_materials)
real(amrex_real) PRES
real(amrex_real) PRES_in
integer dir,side
real(amrex_real) dx(SDIM)

if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then

 call CAV2Dstep_PRES(xghost,t,LS,PRES)

else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine CAV2Dstep_PRES_BC

 ! dir=1..sdim  side=1..2
subroutine CAV2Dstep_STATE_BC(xwall,xghost,t,LS, &
   STATE,STATE_merge,STATE_in,im,istate,dir,side,dx)
use probcommon_module
use global_utility_module
IMPLICIT NONE

real(amrex_real) xwall
real(amrex_real) xghost(SDIM)
real(amrex_real) t
real(amrex_real) LS(num_materials)
real(amrex_real) local_STATE(num_materials*num_state_material)
real(amrex_real) STATE
real(amrex_real) STATE_merge
real(amrex_real) STATE_in
integer dir,side
real(amrex_real) dx(SDIM)
integer istate,im
integer ibase,im_crit

if ((istate.ge.1).and. &
    (istate.le.num_state_material).and. &
    (im.ge.1).and. &
    (im.le.num_materials)) then
 call CAV2Dstep_STATE(xghost,t,LS,local_STATE)
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
end subroutine CAV2Dstep_STATE_BC

subroutine CAV2Dstep_HEATSOURCE(im,VFRAC,time,x,temp, &
     heat_source,den,CV,dt)
use probcommon_module
IMPLICIT NONE

integer im
real(amrex_real) VFRAC(num_materials)
real(amrex_real) time
real(amrex_real) x(SDIM)
real(amrex_real) temp(num_materials)
real(amrex_real) den(num_materials)
real(amrex_real) CV(num_materials)
real(amrex_real) dt
real(amrex_real) heat_source

if ((num_materials.eq.3).and.(probtype.eq.412)) then
 heat_source=zero
else
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine CAV2Dstep_HEATSOURCE

end module CAV2Dstep_module
