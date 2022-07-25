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

! probtype==415 (see run2d/inputs.shock_solid_sphere) 
!  use probtype==401 (run3d/inputs.HELIX) as a guide.
module MEHDI_SHOCK_SPHERE

implicit none                   

contains

  ! do any initial preparation needed
subroutine INIT_MEHDI_SHOCK_SPHERE()
IMPLICIT NONE

return
end subroutine INIT_MEHDI_SHOCK_SPHERE

subroutine MEHDI_SS_LS(x,t,LS)
use probcommon_module
use global_utility_module
IMPLICIT NONE

REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(out) :: LS(num_materials)

if ((num_materials.eq.3).and.(probtype.eq.415)) then

  do im=1,num_materials
   if (im.eq.1) then !liquid
    LS(im)=-99999.0
   else if (im.eq.2) then !gas
    LS(im)=99999.0
   else if (im.eq.3) then ! geometry (placeholder)
    LS(im)=-99999.0
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
end subroutine MEHDI_SS_LS

! initial velocity is zero
subroutine MITSUHIRO_LS_VEL(x,t,LS,VEL,velsolid_flag,dx)
use probcommon_module
IMPLICIT NONE

REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: dx(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(num_materials)
REAL_T, INTENT(out) :: VEL(SDIM)
INTEGER_T dir
INTEGER_T, INTENT(in) :: velsolid_flag

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

return 
end subroutine MITSUHIRO_LS_VEL

! this routine used if pressure boundary conditions are prescribed,
! since only top wall is "outflow" (outflow in quotes since ice shrinks when
! melting), and flow is incompressible, ok to make the top wall pressure zero.
subroutine MITSUHIRO_PRES(x,t,LS,PRES)
use probcommon_module
IMPLICIT NONE

REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(num_materials)
REAL_T, INTENT(out) :: PRES

PRES=zero

return 
end subroutine MITSUHIRO_PRES



subroutine MITSUHIRO_STATE(x,t,LS,STATE,bcflag)
use probcommon_module
use global_utility_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: bcflag !0=called from initialize  1=called from bc
REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(num_materials)
REAL_T, INTENT(out) :: STATE(num_materials*num_state_material)
INTEGER_T im,ibase,n

if ((num_materials.eq.4).and. &
    (num_state_material.ge.3).and. & ! density, temperature, vapor spec
    (probtype.eq.414)) then
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
   ! always assume Dirichlet boundary condition at zlo for temperature.
  call outside_temperature(t,x(1),x(2),x(SDIM), &
     STATE(ibase+ENUM_TEMPERATUREVAR+1),im,bcflag)

   ! initial species in inputs?
  do n=1,num_species_var
   STATE(ibase+ENUM_SPECIESVAR+n)=fort_speciesconst((n-1)*num_materials+im)
  enddo
 enddo ! im=1..num_materials
else
 print *,"num_materials,num_state_material, or probtype invalid"
 stop
endif
 
return
end subroutine MITSUHIRO_STATE

 ! dir=1..sdim  side=1..2
subroutine MITSUHIRO_LS_BC(xwall,xghost,t,LS, &
   LS_in,dir,side,dx)
use probcommon_module
IMPLICIT NONE

REAL_T, INTENT(in) :: xwall
REAL_T, INTENT(in) :: xghost(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(out) :: LS(num_materials)
REAL_T, INTENT(in) :: LS_in(num_materials)
INTEGER_T, INTENT(in) :: dir,side
REAL_T, INTENT(in) ::  dx(SDIM)

if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then
 call MITSUHIRO_LS(xghost,t,LS)
else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine MITSUHIRO_LS_BC


 ! dir=1..sdim  side=1..2 veldir=1..sdim
subroutine MITSUHIRO_VEL_BC(xwall,xghost,t,LS, &
   VEL,VEL_in,veldir,dir,side,dx)
use probcommon_module
IMPLICIT NONE

REAL_T, INTENT(in) :: xwall
REAL_T, INTENT(in) :: xghost(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(num_materials)
REAL_T, INTENT(out) :: VEL
REAL_T, INTENT(in) :: VEL_in
INTEGER_T, INTENT(in) :: veldir,dir,side
REAL_T, INTENT(in) :: dx(SDIM)
REAL_T local_VEL(SDIM)
INTEGER_T velsolid_flag

velsolid_flag=0
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2).and. &
    (veldir.ge.1).and.(veldir.le.SDIM)) then

 call MITSUHIRO_LS_VEL(xghost,t,LS,local_VEL,velsolid_flag,dx)
 VEL=local_VEL(veldir)

else
 print *,"dir,side, or veldir invalid"
 stop
endif

return
end subroutine MITSUHIRO_VEL_BC

 ! dir=1..sdim  side=1..2
subroutine MITSUHIRO_PRES_BC(xwall,xghost,t,LS, &
   PRES,PRES_in,dir,side,dx)
use probcommon_module
IMPLICIT NONE

REAL_T, INTENT(in) :: xwall
REAL_T, INTENT(in) :: xghost(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(num_materials)
REAL_T, INTENT(out) :: PRES
REAL_T, INTENT(in) :: PRES_in
INTEGER_T, INTENT(in) :: dir,side
REAL_T, INTENT(in) :: dx(SDIM)

if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then

 call MITSUHIRO_PRES(xghost,t,LS,PRES)

else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine MITSUHIRO_PRES_BC

 ! dir=1..sdim  side=1..2
subroutine MITSUHIRO_STATE_BC(xwall,xghost,t,LS, &
   STATE,STATE_merge,STATE_in,im,istate,dir,side,dx)
use probcommon_module
IMPLICIT NONE

REAL_T, INTENT(in) :: xwall
REAL_T, INTENT(in) :: xghost(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(num_materials)
REAL_T :: local_STATE(num_materials*num_state_material)
REAL_T, INTENT(out) :: STATE
REAL_T, INTENT(out) :: STATE_merge
REAL_T, INTENT(in) :: STATE_in
INTEGER_T, INTENT(in) :: dir,side
REAL_T, INTENT(in) :: dx(SDIM)
INTEGER_T istate,im
INTEGER_T ibase,im_crit,im_loop
INTEGER_T local_bcflag

local_bcflag=1

if ((istate.ge.1).and. &
    (istate.le.num_state_material).and. &
    (im.ge.1).and. &
    (im.le.num_materials)) then
 call MITSUHIRO_STATE(xghost,t,LS,local_STATE,local_bcflag)
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
end subroutine MITSUHIRO_STATE_BC

subroutine MITSUHIRO_HEATSOURCE(im,VFRAC,time,x,temp, &
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

if ((num_materials.eq.4).and.(probtype.eq.414)) then
 heat_source=zero
else
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine MITSUHIRO_HEATSOURCE

end module MITSUHIRO_MELTING_module
