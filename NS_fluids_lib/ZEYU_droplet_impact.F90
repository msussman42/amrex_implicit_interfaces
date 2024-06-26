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

! probtype==413 (see run2d/inputs.ZEYU_droplet_impact)
module ZEYU_droplet_impact_module
use amrex_fort_module, only : amrex_real

implicit none                   

contains

  ! do any initial preparation needed
subroutine INIT_ZEYU_droplet_impact_MODULE()
IMPLICIT NONE

return
end subroutine INIT_ZEYU_droplet_impact_MODULE

! Phi>0 in the solid
subroutine ZEYU_substrateLS(x,Phi) 
use probcommon_module
implicit none
real(amrex_real), INTENT(in), dimension(SDIM) :: x !spatial coordinates
real(amrex_real), INTENT(out) :: Phi !LS dist, Phi>0 in the substrate

real(amrex_real) substrate_height

if (abs(zblob2-yblob2).le.EPS14) then
 substrate_height=zblob2  ! substrate thickness
else
 print *,"zblob2 or yblob2 invalid"
 stop
endif

if (abs(x(SDIM)).le.1.0D+20) then
 Phi=substrate_height-x(SDIM)
else
 print *,"x(SDIM) invalid"
 stop
endif

end subroutine ZEYU_substrateLS


 ! fluids tessellate the domain, solids are immersed. 
subroutine ZEYU_droplet_impact_LS(x,t,LS)
use probcommon_module
IMPLICIT NONE

real(amrex_real) x(SDIM)
real(amrex_real) t
real(amrex_real) LS(num_materials)

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

 call ZEYU_substrateLS(x,LS(3))
else
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine ZEYU_droplet_impact_LS

subroutine ZEYU_droplet_impact_LS_VEL(x,t,LS,VEL,velsolid_flag,dx)
use probcommon_module
IMPLICIT NONE

real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: dx(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: LS(num_materials)
real(amrex_real), INTENT(out) :: VEL(SDIM)
integer dir
integer, INTENT(in) :: velsolid_flag

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
 print *,"expecting adv_dir = SDIM in ZEYU_droplet_impact_LS_VEL"
 stop
endif


return 
end subroutine ZEYU_droplet_impact_LS_VEL

subroutine ZEYU_droplet_impact_PRES(x,t,LS,PRES)
use probcommon_module
IMPLICIT NONE

real(amrex_real) x(SDIM)
real(amrex_real) t
real(amrex_real) LS(num_materials)
real(amrex_real) PRES

PRES=outflow_pressure

return 
end subroutine ZEYU_droplet_impact_PRES


subroutine ZEYU_droplet_impact_STATE(x,t,LS,STATE)
use probcommon_module
IMPLICIT NONE

real(amrex_real) x(SDIM)
real(amrex_real) t
real(amrex_real) LS(num_materials)
real(amrex_real) STATE(num_materials*num_state_material)
integer im,ibase,n

if ((num_materials.eq.3).and. &
    (num_state_material.ge.2).and. &
    (probtype.eq.413)) then
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
end subroutine ZEYU_droplet_impact_STATE

 ! dir=1..sdim  side=1..2
subroutine ZEYU_droplet_impact_LS_BC(xwall,xghost,t,LS, &
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
 call ZEYU_droplet_impact_LS(xghost,t,LS)
else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine ZEYU_droplet_impact_LS_BC


 ! dir=1..sdim  side=1..2 veldir=1..sdim
subroutine ZEYU_droplet_impact_VEL_BC(xwall,xghost,t,LS, &
   VEL,VEL_in,veldir,dir,side,dx)
use probcommon_module
IMPLICIT NONE

real(amrex_real), INTENT(in) :: xwall
real(amrex_real), INTENT(in) :: xghost(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: LS(num_materials)
real(amrex_real), INTENT(out) :: VEL
real(amrex_real), INTENT(in) :: VEL_in
integer, INTENT(in) :: veldir,dir,side
real(amrex_real), INTENT(in) :: dx(SDIM)
real(amrex_real) local_VEL(SDIM)
integer velsolid_flag

velsolid_flag=0
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2).and. &
    (veldir.ge.1).and.(veldir.le.SDIM)) then

 call ZEYU_droplet_impact_LS_VEL(xghost,t,LS,local_VEL,velsolid_flag,dx)
 VEL=local_VEL(veldir)

else
 print *,"dir,side, or veldir invalid"
 stop
endif

return
end subroutine ZEYU_droplet_impact_VEL_BC

 ! dir=1..sdim  side=1..2
subroutine ZEYU_droplet_impact_PRES_BC(xwall,xghost,t,LS, &
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

 call ZEYU_droplet_impact_PRES(xghost,t,LS,PRES)

else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine ZEYU_droplet_impact_PRES_BC

 ! dir=1..sdim  side=1..2
subroutine ZEYU_droplet_impact_STATE_BC(xwall,xghost,t,LS, &
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
 call ZEYU_droplet_impact_STATE(xghost,t,LS,local_STATE)
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
end subroutine ZEYU_droplet_impact_STATE_BC

subroutine ZEYU_droplet_impact_HEATSOURCE(im,VFRAC,time,x,temp, &
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

if ((num_materials.eq.2).and.(probtype.eq.413)) then
 heat_source=zero
else
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine ZEYU_droplet_impact_HEATSOURCE

end module ZEYU_droplet_impact_module
