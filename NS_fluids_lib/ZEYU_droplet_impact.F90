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

! probtype==413 (see run2d/inputs.ZEYU_droplet_impact)
module ZEYU_droplet_impact_module

implicit none                   

contains

  ! do any initial preparation needed
subroutine INIT_ZEYU_droplet_impact_MODULE()
IMPLICIT NONE

return
end subroutine INIT_ZEYU_droplet_impact_MODULE

! Phi>0 in the solid
subroutine substrateLS(x,Phi) 
use probcommon_module
implicit none
REAL_T, intent(in), dimension(SDIM) :: x !spatial coordinates
REAL_T, intent(out) :: Phi !LS dist, Phi>0 in the substrate

REAL_T substrate_height

if (abs(zblob2-yblob2).le.1.0D-14) then
 substrate_height=zblob2  ! substrate thickness
else
 print *,"zblob2 or yblob2 invalid"
 stop
endif

if (abs(x(SDIM)).le.1.0D-20) then
 Phi=substrate_height-x(SDIM)
else
 print *,"x(SDIM) invalid"
 stop
endif

end subroutine substrateLS


 ! fluids tessellate the domain, solids are immersed. 
subroutine ZEYU_droplet_impact_LS(x,t,LS)
use probcommon_module
IMPLICIT NONE

REAL_T x(SDIM)
REAL_T t
INTEGER_T im
REAL_T LS(num_materials)

if ((num_materials.eq.3).and.(probtype.eq.413)) then
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

 call substrateLS(x,LS(3))
else
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine ZEYU_droplet_impact_LS

subroutine ZEYU_droplet_impact_LS_VEL(x,t,LS,VEL,velsolid_flag,dx)
use probcommon_module
IMPLICIT NONE

REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(num_materials)
REAL_T, intent(out) :: VEL(SDIM)
INTEGER_T dir
INTEGER_T, intent(in) :: velsolid_flag

if ((velsolid_flag.eq.0).or. &
    (velsolid_flag.eq.1)) then
 ! do nothing
else 
 print *,"velsolid_flag invalid"
 stop
endif

if (adv_dir.eq.SDIM) then

  ! material 3 is the substrate
 if ((LS(3).ge.zero).or.(velsolid_flag.eq.1)) then
  ! in solid
  do dir=1,SDIM
   VEL(dir)=zero
  enddo
 else if ((LS(3).le.zero).and. &
          (velsolid_flag.eq.0)) then

     ! material 1 is the drop
  if ((LS(1).ge.zero).or. &
      (LS(1).ge.-dx(1)) then
   ! in drop
   do dir=1,SDIM
    VEL(dir)=-abs(advbot)
   enddo
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
 print *,"expecting adv_dir = SDIM in ZEYU_droplet_impact_LS_VEL"
 stop
endif


return 
end subroutine ZEYU_droplet_impact_LS_VEL

subroutine CAV2Dstep_PRES(x,t,LS,PRES)
use probcommon_module
IMPLICIT NONE

REAL_T x(SDIM)
REAL_T t
REAL_T LS(num_materials)
REAL_T PRES

PRES=outflow_pressure

return 
end subroutine CAV2Dstep_PRES


subroutine CAV2Dstep_STATE(x,t,LS,STATE)
use probcommon_module
IMPLICIT NONE

REAL_T x(SDIM)
REAL_T t
REAL_T LS(num_materials)
REAL_T STATE(num_materials*num_state_material)
INTEGER_T im,ibase,n

if ((num_materials.eq.2).and. &
    (num_state_material.ge.2).and. &
    (probtype.eq.412)) then
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
end subroutine CAV2Dstep_STATE

 ! dir=1..sdim  side=1..2
subroutine CAV2Dstep_LS_BC(xwall,xghost,t,LS, &
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
 call CAV2Dstep_STATE(xghost,t,LS,local_STATE)
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
end subroutine CAV2Dstep_STATE_BC

subroutine CAV2Dstep_HEATSOURCE(im,VFRAC,time,x,temp, &
     heat_source,den,CV,dt)
use probcommon_module
IMPLICIT NONE

INTEGER_T im
REAL_T VFRAC(num_materials)
REAL_T time
REAL_T x(SDIM)
REAL_T temp(num_materials)
REAL_T den(num_materials)
REAL_T CV(num_materials)
REAL_T dt
REAL_T heat_source

if ((num_materials.eq.2).and.(probtype.eq.412)) then
 heat_source=zero
else
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine CAV2Dstep_HEATSOURCE

end module CAV2Dstep_module
