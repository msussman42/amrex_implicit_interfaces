#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_FORT_INTEGER.H"
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

! probtype==82
module ROTATING_ANNULUS_module

implicit none                   

contains

  ! do any initial preparation needed
subroutine INIT_ROTATING_ANNULUS_MODULE()
IMPLICIT NONE

return
end subroutine INIT_ROTATING_ANNULUS_MODULE


subroutine ROTATING_ANNULUS_INTERNAL_GRAVITY_WAVE_FLAG(internal_wave_exists)
use probcommon_module
use global_utility_module
IMPLICIT NONE
INTEGER_T, INTENT(out) :: internal_wave_exists

internal_wave_exists=1

return
end subroutine ROTATING_ANNULUS_INTERNAL_GRAVITY_WAVE_FLAG

 ! fluids tessellate the domain, solids are immersed. 
subroutine ROTATING_ANNULUS_LS(x,t,LS,nmat)
use probcommon_module
use global_utility_module
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
  if (num_materials.eq.2) then 
   ! do nothing
  else
   print *,"expecting num_materials=2 for ROTATING_ANNULUS"
   stop
  endif

  LS(1)=99999.90
  LS(2)=-LS(1)

return
end subroutine ROTATING_ANNULUS_LS

! initial velocity is zero
subroutine ROTATING_ANNULUS_VEL(x,t,LS,VEL,velsolid_flag,dx,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: dx(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(nmat)
REAL_T, INTENT(out) :: VEL(SDIM)
INTEGER_T dir
INTEGER_T, INTENT(in) :: velsolid_flag
REAL_T :: temp
REAL_T :: xmid,zmid

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

do dir=1,SDIM
 VEL(dir)=zero
enddo

if (probtype.eq.82) then
 ! do nothing
else
 print *,"expecting probtype==82"
 stop
endif

return 
end subroutine ROTATING_ANNULUS_VEL


! this routine used as a default when
! pressure boundary conditions are prescribed.
! For the case when only top wall is 
! "outflow" (outflow in quotes since ice shrinks when
! melting), and flow is incompressible, ok to make the top wall pressure zero.
subroutine ROTATING_ANNULUS_PRES(x,t,LS,PRES,nmat)
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
if (probtype.eq.82) then
 ! do nothing
else
 print *,"expecting probtype==82"
 stop
endif
PRES=zero

return 
end subroutine ROTATING_ANNULUS_PRES

! fort_tempconst(1) is the inner wall temperature
! twall is the outer wall temperature
subroutine ROTATING_ANNULUS_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: bcflag !0=called from initialize  1=called from bc
INTEGER_T, INTENT(in) :: nmat
INTEGER_T, INTENT(in) :: nstate_mat
REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(nmat)
REAL_T, INTENT(out) :: STATE(nmat*nstate_mat)
INTEGER_T im,ibase,n
REAL_T t_sloping

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

if (num_state_material.eq.2) then
 ! do nothing
else
 print *,"expecting num_state_material.eq.2"
 stop
endif

if (levelrz.eq.3) then
 ! do nothing
else
 print *,"expecting levelrz=3"
 stop
endif

if ((problenx.eq.problen_array(1)).and. &
    (problenx.gt.zero).and. &
    (problox.gt.zero).and. &
    (x(1).gt.zero)) then
 ! do nothing
else
 print *,"problenx or problox or x(1) invalid"
 stop
endif

if (twall.ge.fort_tempconst(1)) then
 ! do nothing
else
 print *,"TA or TB invalid"
 stop
endif

if (probtype.eq.82) then

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

  if (x(1).lt.problox) then
   t_sloping=fort_tempconst(1)
  else if (x(1).gt.probhix) then
   t_sloping=twall
  else if ((x(1).ge.problox).and. &
           (x(1).le.probhix)) then
   t_sloping=fort_tempconst(1)+ &
      (twall-fort_tempconst(1))*(x(1)-problox)/problenx
  else
   print *,"x(1) invalid"
   stop
  endif

  if (t_sloping.gt.zero) then
   STATE(ibase+ENUM_TEMPERATUREVAR+1)=t_sloping
  else
   print *,"need t_sloping>0"
   stop
  endif

   ! initial species in inputs?
  do n=1,num_species_var
   STATE(ibase+ENUM_SPECIESVAR+n)=fort_speciesconst((n-1)*num_materials+im)
  enddo

 enddo ! im=1..num_materials

else
 print *,"num_materials,num_state_material, or probtype invalid"
 print *,"aborting ROTATING_ANNULUS_STATE"
 stop
endif
 
return
end subroutine ROTATING_ANNULUS_STATE

 ! dir=1..sdim  side=1..2
subroutine ROTATING_ANNULUS_LS_BC(xwall,xghost,t,LS, &
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
REAL_T, INTENT(in) ::  dx(SDIM)

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then
 call ROTATING_ANNULUS_LS(xghost,t,LS,nmat)
else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine ROTATING_ANNULUS_LS_BC


 ! dir=1..sdim  side=1..2 veldir=1..sdim
subroutine ROTATING_ANNULUS_VEL_BC(xwall,xghost,t,LS, &
   VEL,VEL_in,veldir,dir,side,dx,nmat)
use probcommon_module
use global_utility_module
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
REAL_T temp
INTEGER_T for_dt

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
if (probtype.eq.82) then
 ! do nothing
else
 print *,"expecting probtype==82"
 stop
endif
velsolid_flag=0
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2).and. &
    (veldir.ge.1).and.(veldir.le.SDIM)) then

 call ROTATING_ANNULUS_VEL(xghost,t,LS,local_VEL,velsolid_flag,dx,nmat)

else
 print *,"dir,side, or veldir invalid"
 stop
endif

return
end subroutine ROTATING_ANNULUS_VEL_BC

 ! dir=1..sdim  side=1..2
subroutine ROTATING_ANNULUS_PRES_BC(xwall,xghost,t,LS, &
   PRES,PRES_in,dir,side,dx,nmat)
use probcommon_module
use global_utility_module
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
if (probtype.eq.82) then
 ! do nothing
else
 print *,"expecting probtype==82"
 stop
endif

if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then

 call ROTATING_ANNULUS_PRES(xghost,t,LS,PRES,nmat)

else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine ROTATING_ANNULUS_PRES_BC

 ! dir=1..sdim  side=1..2
subroutine ROTATING_ANNULUS_STATE_BC(xwall,xghost,t,LS, &
   STATE,STATE_merge,STATE_in,im,istate,dir,side,dx,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
REAL_T, INTENT(in) :: xwall
REAL_T, INTENT(in) :: xghost(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(nmat)
REAL_T :: local_STATE(nmat*num_state_material)
REAL_T, INTENT(inout) :: STATE
REAL_T, INTENT(inout) :: STATE_merge
REAL_T, INTENT(in) :: STATE_in
INTEGER_T, INTENT(in) :: dir,side
REAL_T, INTENT(in) :: dx(SDIM)
INTEGER_T, INTENT(in) :: istate,im
INTEGER_T ibase,im_crit
INTEGER_T local_bcflag

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

if (num_materials.eq.2) then
 ! do nothing
else
 print *,"expecting num_materials.eq.2"
 stop
endif

if (levelrz.eq.3) then
 ! do nothing
else
 print *,"expecting levelrz=3"
 stop
endif

if ((problenx.eq.problen_array(1)).and. &
    (problenx.gt.zero).and. &
    (problox.gt.zero).and. &
    (xghost(1).gt.zero)) then
 ! do nothing
else
 print *,"problenx or problox or xghost(1) invalid"
 stop
endif
if (dir.eq.1) then
 if (xwall.gt.zero) then
  ! do nothing
 else
  print *,"xwall invalid"
  stop
 endif
else if ((dir.eq.2).or.(dir.eq.3)) then
 ! do nothing
else
 print *,"dir invalid"
 stop
endif

if (num_state_material.eq.2) then
 ! do nothing
else
 print *,"expecting num_state_material.eq.2"
 stop
endif

if (probtype.eq.82) then
 ! do nothing
else
 print *,"expecting probtype.eq.82"
 stop
endif

! fort_tempconst(1) is the inner wall temperature
! twall is the outer wall temperature
local_bcflag=1

if ((istate.ge.1).and. &
    (istate.le.num_state_material).and. &
    (im.ge.1).and. &
    (im.le.num_materials)) then
 call ROTATING_ANNULUS_STATE(xghost,t,LS,local_STATE,local_bcflag, &
         nmat,num_state_material)
 ibase=(im-1)*num_state_material
 STATE=local_STATE(ibase+istate)
 call get_primary_material(LS,im_crit)
 if (im_crit.eq.1) then
  ! do nothing
 else
  print *,"expecting im_crit=1"
  stop
 endif
 ibase=(im_crit-1)*num_state_material
 STATE_merge=local_STATE(ibase+istate)
else
 print *,"istate invalid"
 stop
endif

return
end subroutine ROTATING_ANNULUS_STATE_BC

end module ROTATING_ANNULUS_module
