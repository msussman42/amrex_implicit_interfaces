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

! probtype==2000 (see run2d/inputs.flexible_plate_impact)
module flexible_plate_impact_module

implicit none 

REAL_T :: DEF_VAPOR_GAMMA

contains

  ! do any initial preparation needed
subroutine INIT_flexible_plate_impact_MODULE()
IMPLICIT NONE

  DEF_VAPOR_GAMMA =  1.666666667D0

return
end subroutine INIT_flexible_plate_impact_MODULE

! Phi>0 in the plate
subroutine flexible_substrateLS(x,Phi) 
use probcommon_module
use global_utility_module
implicit none
REAL_T, INTENT(in), dimension(SDIM) :: x !spatial coordinates
REAL_T, INTENT(out) :: Phi !LS dist, Phi>0 in the substrate

REAL_T substrate_height

if ((radblob2.gt.zero).and. &
    (radblob3.gt.zero).and. &
    (radblob4.gt.zero)) then
 ! do nothing
else
 print *,"radblob2, radblob3, or radblob4 invalid"
 stop
endif

if (SDIM.eq.2) then
 call squaredist(x(1),x(2), & ! substrate_height<0 in object
   xblob2-radblob2, &
   xblob2+radblob2, &
   yblob2-radblob3, &
   yblob2+radblob3, &
   substrate_height)
else if (SDIM.eq.3) then
 call cubedist( &
   xblob2-radblob2, &
   xblob2+radblob2, &
   yblob2-radblob3, &
   yblob2+radblob3, &
   zblob2-radblob4, &
   zblob2+radblob4, &
   x(1),x(2),x(SDIM), &
   substrate_height) ! substrate_height<0 in object
else
 print *,"dimension bust"
 stop
endif

 Phi=-substrate_height

end subroutine flexible_substrateLS


 ! fluids tessellate the domain, solids are immersed. 
subroutine flexible_plate_impact_LS(x,t,LS,nmat)
use probcommon_module
IMPLICIT NONE

  INTEGER_T, INTENT(in) :: nmat
  REAL_T, INTENT(in) :: x(SDIM)
  REAL_T, INTENT(in) :: t
  REAL_T, INTENT(out) :: LS(nmat)

  if (nmat.eq.num_materials) then
   ! do nothing
  else
   print *,"nmat invalid"
   stop
  endif

if ((num_materials.eq.3).and.(probtype.eq.2000)) then
 ! liquid
 if (SDIM.eq.3) then
  LS(1)=radblob-sqrt((x(1)-xblob)**2+(x(2)-yblob)**2+(x(SDIM)-zblob)**2)
 else if (SDIM.eq.2) then
  LS(1)=radblob-sqrt((x(1)-xblob)**2+(x(2)-yblob)**2)
 else
  print *,"dimension bust"
  stop
 endif
 LS(2)=-LS(1)

 call flexible_substrateLS(x,LS(3))

 if (LS(2).ge.zero) then
  ! intersection of the complements of the liquid and plate
  LS(2)=min(-LS(1),-LS(3)) 
 else if (LS(2).le.zero) then
  ! do nothing
 else
  print *,"LS(2) invalid"
  stop
 endif

else
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine flexible_plate_impact_LS

subroutine flexible_plate_check_vel_rigid(x,t,vel,dir)
use probcommon_module
IMPLICIT NONE

REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: vel
INTEGER_T, INTENT(in) :: dir

if (t.ge.0.0d0) then
 ! do nothing
else
 print *,"t invalid"
 stop
endif
if ((dir.ge.1).and.(dir.le.SDIM)) then
 ! do nothing
else
 print *,"dir invalid"
 stop
endif

if (probtype.eq.2000) then
 if (vel.eq.0.0d0) then
  ! do nothing
 else
  print *,"flexible_plate_check_vel_rigid: vel not expected"
  print *,"t,dir,vel ",t,dir,vel
  stop
 endif
else
 print *,"probtype invalid"
 stop
endif

return
end subroutine flexible_plate_check_vel_rigid

subroutine flexible_plate_clamped_LS(x,t,LS,vel, &
     temperature,prescribed_flag,dx)
use probcommon_module
use global_utility_module
IMPLICIT NONE

  REAL_T, INTENT(in) :: x(SDIM)
  REAL_T, INTENT(in) :: dx(SDIM)
  REAL_T, INTENT(in) :: t
  REAL_T, INTENT(out) :: LS
  REAL_T, INTENT(out) :: vel(SDIM)
  REAL_T, INTENT(out) :: temperature
  INTEGER_T, INTENT(out) :: prescribed_flag
  REAL_T :: LS_left,LS_right
  INTEGER_T :: dir


if (probtype.eq.2000) then

 prescribed_flag=0

 if (SDIM.eq.2) then
  call squaredist(x(1),x(2), & ! LS_right<0 in object
   xblob2+radblob2-radblob5, &
   xblob2+radblob2, &
   yblob2-radblob3, &
   yblob2+radblob3, &
   LS_right)
  call squaredist(x(1),x(2), & ! LS_left<0 in object
   xblob2-radblob2, &
   xblob2-radblob2+radblob5, &
   yblob2-radblob3, &
   yblob2+radblob3, &
   LS_left)
  LS_left=-LS_left
  LS_right=-LS_right
  LS=max(LS_left,LS_right)
  do dir=1,SDIM
   vel(dir)=zero
  enddo
  temperature=293.0d0  ! room temperature
 else
  print *,"this code only for 2d for now"
  stop
 endif
else
 print *,"probtype invalid"
 stop
endif

return
end subroutine flexible_plate_clamped_LS



subroutine flexible_plate_impact_VEL(x,t,LS,VEL,velsolid_flag,dx,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: dx(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(nmat)
REAL_T, INTENT(out) :: VEL(SDIM)
INTEGER_T dir
INTEGER_T, INTENT(in) :: velsolid_flag

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

do dir=1,SDIM
 if (dx(dir).gt.zero) then
  ! do nothing
 else
  print *,"dx invalid"
  stop
 endif
enddo

if (adv_dir.eq.SDIM) then

  ! material 3 is the substrate
  ! velsolid_flag==1 if initializing the solid velocity.
 if ((LS(3).ge.zero).or.(velsolid_flag.eq.1)) then
  ! in solid
  do dir=1,SDIM
   VEL(dir)=zero
  enddo
 else if ((LS(3).le.zero).and. &
          (velsolid_flag.eq.0)) then

     ! material 1 is the drop
  if ((LS(1).ge.zero).or. &
      (LS(1).ge.-two*dx(1))) then
   ! in drop
   do dir=1,SDIM
    VEL(dir)=zero
   enddo
   VEL(SDIM)=-abs(advbot)
  else if (LS(2).ge.zero) then ! material 2 is the gas

   do dir=1,SDIM
    VEL(dir)=zero
   enddo

  else
   print *,"LS bust"
   stop
  endif
 else
  print *,"LS(3) or velsolid bust"
  stop
 endif

else
 print *,"expecting adv_dir = SDIM in flexible_plate_impact_VEL"
 stop
endif


return 
end subroutine flexible_plate_impact_VEL



! These next routines only used for compressible materials.
!***********************************************
! compressible material functions for (ns.material_type = 24)
subroutine EOS_flexible_plate_impact(rho,internal_energy,pressure, &
  imattype,im)
 IMPLICIT NONE
 INTEGER_T, INTENT(in) :: imattype,im
 REAL_T, INTENT(in) :: rho
 REAL_T, INTENT(in) :: internal_energy
 REAL_T, INTENT(out) :: pressure

 if (imattype.eq.24) then
  pressure=rho*(DEF_VAPOR_GAMMA-1.0D0)*internal_energy
 else
  print *,"imattype invalid EOS_flexible_plate_impact"
  stop
 endif

 return
end subroutine EOS_flexible_plate_impact

subroutine SOUNDSQR_flexible_plate_impact(rho,internal_energy,soundsqr, &
  imattype,im)
 IMPLICIT NONE
 INTEGER_T, INTENT(in) :: imattype,im
 REAL_T, INTENT(in) :: rho
 REAL_T, INTENT(in) :: internal_energy
 REAL_T, INTENT(out) :: soundsqr
 REAL_T pressure

 if (imattype.eq.24) then
  call EOS_flexible_plate_impact(rho,internal_energy,pressure,imattype,im)
  soundsqr=DEF_VAPOR_GAMMA*pressure/rho
 else
  print *,"imattype invalid SOUNDSQR_flexible_plate_impact"
  stop
 endif

 return
end subroutine SOUNDSQR_flexible_plate_impact

subroutine INTERNAL_flexible_plate_impact(rho,temperature, &
  local_internal_energy, &
  imattype,im)
 use global_utility_module
 IMPLICIT NONE
 INTEGER_T, INTENT(in) :: imattype,im
 REAL_T, INTENT(in) :: rho
 REAL_T, INTENT(in) :: temperature 
 REAL_T, INTENT(out) :: local_internal_energy

 call INTERNAL_default(rho,temperature,local_internal_energy, &
        imattype,im)

 return
end subroutine INTERNAL_flexible_plate_impact

subroutine TEMPERATURE_flexible_plate_impact(rho,temperature,internal_energy, &
  imattype,im)
 use global_utility_module
 IMPLICIT NONE
 INTEGER_T, INTENT(in) :: imattype,im
 REAL_T, INTENT(in) :: rho
 REAL_T, INTENT(out) :: temperature 
 REAL_T, INTENT(in) :: internal_energy

 call TEMPERATURE_default(rho,temperature,internal_energy, &
        imattype,im)

 return
end subroutine TEMPERATURE_flexible_plate_impact

! This routine will not effect the simulation since
! all of the domain BC will be no-slip.
subroutine flexible_plate_impact_PRES(x,t,LS,PRES,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(nmat)
REAL_T, INTENT(out) :: PRES

if (num_materials.eq.nmat) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
PRES=outflow_pressure

return 
end subroutine flexible_plate_impact_PRES


subroutine flexible_plate_impact_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: bcflag !0=called from initialize  1=called from bc
INTEGER_T, INTENT(in) :: nmat
INTEGER_T, INTENT(in) :: nstate_mat
REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(nmat)
REAL_T, INTENT(out) :: STATE(nmat*nstate_mat)
INTEGER_T im,ibase,n

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
    (probtype.eq.2000)) then
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
end subroutine flexible_plate_impact_STATE

 ! dir=1..sdim  side=1..2
subroutine flexible_plate_impact_LS_BC(xwall,xghost,t,LS, &
   LS_in,dir,side,dx,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
REAL_T, INTENT(in) :: xwall
REAL_T, INTENT(in) :: xghost(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(inout) :: LS(nmat)
REAL_T, INTENT(in) :: LS_in(nmat)
INTEGER_T, INTENT(in) :: dir,side
REAL_T, INTENT(in) :: dx(SDIM)

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then
 call flexible_plate_impact_LS(xghost,t,LS,nmat)
else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine flexible_plate_impact_LS_BC


 ! dir=1..sdim  side=1..2 veldir=1..sdim
subroutine flexible_plate_impact_VEL_BC(xwall,xghost,t,LS, &
   VEL,VEL_in,veldir,dir,side,dx,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
REAL_T, INTENT(in) :: xwall
REAL_T, INTENT(in) :: xghost(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(nmat)
REAL_T, INTENT(inout) :: VEL
REAL_T, INTENT(in) :: VEL_in
INTEGER_T, INTENT(in) :: veldir,dir,side
REAL_T, INTENT(in) :: dx(SDIM)
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

 call flexible_plate_impact_VEL(xghost,t,LS,local_VEL,velsolid_flag,dx,nmat)
 VEL=local_VEL(veldir)

else
 print *,"dir,side, or veldir invalid"
 stop
endif

return
end subroutine flexible_plate_impact_VEL_BC

 ! dir=1..sdim  side=1..2
subroutine flexible_plate_impact_PRES_BC(xwall,xghost,t,LS, &
   PRES,PRES_in,dir,side,dx,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
REAL_T, INTENT(in) :: xwall
REAL_T, INTENT(in) :: xghost(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(nmat)
REAL_T, INTENT(inout) :: PRES
REAL_T, INTENT(in) :: PRES_in
INTEGER_T, INTENT(in) :: dir,side
REAL_T, INTENT(in) :: dx(SDIM)

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif


if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then

 call flexible_plate_impact_PRES(xghost,t,LS,PRES,nmat)

else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine flexible_plate_impact_PRES_BC

 ! dir=1..sdim  side=1..2
subroutine flexible_plate_impact_STATE_BC(xwall,xghost,t,LS, &
   STATE,STATE_merge,STATE_in,im,istate,dir,side,dx,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
REAL_T, INTENT(in) :: xwall
REAL_T, INTENT(in) :: xghost(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(nmat)
REAL_T local_STATE(nmat*num_state_material)
REAL_T, INTENT(inout) :: STATE
REAL_T, INTENT(inout) :: STATE_merge
REAL_T, INTENT(in) :: STATE_in
INTEGER_T, INTENT(in) :: dir,side
REAL_T, INTENT(in) :: dx(SDIM)
INTEGER_T, INTENT(in) :: istate,im
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
 call flexible_plate_impact_STATE(xghost,t,LS,local_STATE, &
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
end subroutine flexible_plate_impact_STATE_BC

subroutine flexible_plate_impact_HEATSOURCE( &
     im,VFRAC, &
     time, &
     x, &
     xsten, & ! xsten(-nhalf:nhalf,SDIM)
     nhalf, &
     temp, &
     heat_source,den,CV,dt, &
     nmat)
use probcommon_module
IMPLICIT NONE

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

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

if ((num_materials.eq.3).and.(probtype.eq.2000)) then
 heat_source=zero
else
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine flexible_plate_impact_HEATSOURCE

end module flexible_plate_impact_module
