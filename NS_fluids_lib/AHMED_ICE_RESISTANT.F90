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

! probtype==425 (see run2d/inputs.AHMED_ICE_RESISTANT or
! run3d/inputs.AHMED_ICE_RESISTANT3d)
module AHMED_ICE_RESISTANT_module

implicit none                   

REAL_T :: DEF_VAPOR_GAMMA

contains

  ! do any initial preparation needed
subroutine INIT_AHMED_ICE_RESISTANT_MODULE()
IMPLICIT NONE

  DEF_VAPOR_GAMMA =  1.666666667D0

return
end subroutine INIT_AHMED_ICE_RESISTANT_MODULE

! Phi>0 in the solid
subroutine AHMED_substrateLS(x,Phi) 
use probcommon_module
use global_utility_module
implicit none
REAL_T, intent(in), dimension(SDIM) :: x !spatial coordinates
REAL_T, intent(out) :: Phi !LS dist, Phi>0 in the substrate
REAL_T :: yhalf,xshift
REAL_T :: local_time
INTEGER_T :: im
INTEGER_T :: expected_nmat

if (axis_dir.eq.0) then
 expected_nmat=4
else if (axis_dir.eq.1) then
 expected_nmat=3
else if (axis_dir.eq.2) then
 expected_nmat=4
else
 print *,"axis_dir invalid"
 stop
endif

if ((num_materials.eq.expected_nmat).and.(probtype.eq.425)) then
 local_time=zero
 im=num_materials
 if (SDIM.eq.2) then
  yhalf=0.2d0
  xshift=x(1)+0.2d0
 else if (SDIM.eq.3) then
  yhalf=x(2)
  xshift=x(1)
 else
  print *,"dimension bust"
  stop
 endif
  ! patterned_substrates is declared in GLOBALUTIL.F90.
  ! Phi<0 in the substrate, Phi>0 in the fluid
 call patterned_substrates(xshift,yhalf,x(SDIM),Phi,local_time,im)
 Phi=-Phi
else
 print *,"num_materials or probtype invalid"
 stop
endif

end subroutine AHMED_substrateLS

subroutine AHMED_ICE_RESISTANT_LS(x,t,LS,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(out) :: LS(nmat)
REAL_T :: ice_vertical
REAL_T :: substrate_height
INTEGER_T :: expected_nmat

  if (nmat.eq.num_materials) then
   ! do nothing
  else
   print *,"nmat invalid"
   stop
  endif

if (axis_dir.eq.0) then
 expected_nmat=4
else if (axis_dir.eq.1) then
 expected_nmat=3
else if (axis_dir.eq.2) then
 expected_nmat=4
else
 print *,"axis_dir invalid"
 stop
endif

if ((num_materials.eq.expected_nmat).and.(probtype.eq.425)) then
  ! fluids tessellate the domain, substrate is embedded.

  ! water is material 1
 if (SDIM.eq.2) then
  LS(1)=radblob-sqrt((x(1)-xblob)**2+(x(2)-yblob)**2)
 else if (SDIM.eq.3) then
  LS(1)=radblob-sqrt((x(1)-xblob)**2+(x(2)-yblob)**2+ &
                     (x(SDIM)-zblob)**2)
 else
  print *,"dimension bust"
  stop
 endif

 if ((axis_dir.eq.0).or.(axis_dir.eq.2)) then
  ! oil is material 3
  ! zblob2 is the altitude of the oil layer
  LS(3)=zblob2-x(SDIM)
  ! air
  LS(2)=-(max(LS(1),LS(3)))
 else if (axis_dir.eq.1) then
  LS(2)=-LS(1)
 else
  print *,"axis_dir invalid"
  stop
 endif

 call AHMED_substrateLS(x,LS(num_materials))
else
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine AHMED_ICE_RESISTANT_LS

! initial velocity is zero
subroutine AHMED_ICE_RESISTANT_VEL(x,t,LS,VEL,velsolid_flag,dx,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: dx(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(nmat)
REAL_T, intent(out) :: VEL(SDIM)
INTEGER_T dir
INTEGER_T, intent(in) :: velsolid_flag

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

  if (axis_dir.eq.0) then

   if (LS(1).ge.-dx(1)) then
    VEL(SDIM)=-abs(adv_vel)
   else if (LS(1).le.-dx(1)) then
    ! do nothing
   else
    print *,"LS(1) invalid"
    stop
   endif

  else if ((axis_dir.eq.1).or. &
           (axis_dir.eq.2)) then
   VEL(1)=adv_vel

   if (LS(1).ge.-dx(1)) then
    VEL(SDIM)=-abs(adv_vel)
   else if (LS(1).le.-dx(1)) then
    VEL(1)=xblob3
    if (SDIM.eq.2) then
     VEL(SDIM)=yblob3
    else if (SDIM.eq.3) then
     VEL(SDIM)=zblob3
    else
     print *,"SDIM invalid"
     stop
    endif
   else
    print *,"LS(1) invalid"
    stop
   endif

  else
   print *,"axis_dir invalid"
   stop
  endif

return 
end subroutine AHMED_ICE_RESISTANT_VEL


! this routine used if pressure boundary conditions are prescribed,
! since only top wall is "outflow" (outflow in quotes since ice shrinks when
! melting), and flow is incompressible, ok to make the top wall pressure zero.
subroutine AHMED_ICE_RESISTANT_PRES(x,t,LS,PRES,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(nmat)
REAL_T, intent(out) :: PRES

if (num_materials.eq.nmat) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
PRES=zero

return 
end subroutine AHMED_ICE_RESISTANT_PRES



subroutine AHMED_ICE_RESISTANT_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
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
INTEGER_T :: expected_nmat

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

if (axis_dir.eq.0) then
 expected_nmat=4
else if (axis_dir.eq.1) then
 expected_nmat=3
else if (axis_dir.eq.2) then
 expected_nmat=4
else
 print *,"axis_dir invalid"
 stop
endif

if ((num_materials.eq.expected_nmat).and. &
    (num_state_material.ge.2).and. & ! density, temperature, vapor spec
    (probtype.eq.425)) then
 do im=1,num_materials
  ibase=(im-1)*num_state_material
  STATE(ibase+1)=fort_denconst(im) ! density prescribed in the inputs file.
  if (t.eq.zero) then
   STATE(ibase+2)=fort_initial_temperature(im) !initial temperature in inputs
  else if (t.gt.zero) then
   STATE(ibase+2)=fort_tempconst(im)
  else
   print *,"t invalid"
   stop
  endif
   ! initial species in inputs?
  do n=1,num_species_var
   STATE(ibase+2+n)=fort_speciesconst((n-1)*num_materials+im)
  enddo
 enddo ! im=1..num_materials
else
 print *,"num_materials,num_state_material, or probtype invalid"
 stop
endif
 
return
end subroutine AHMED_ICE_RESISTANT_STATE

 ! dir=1..sdim  side=1..2
subroutine AHMED_ICE_RESISTANT_LS_BC(xwall,xghost,t,LS, &
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
REAL_T, intent(in) ::  dx(SDIM)

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then
 call AHMED_ICE_RESISTANT_LS(xghost,t,LS,nmat)
else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine AHMED_ICE_RESISTANT_LS_BC


 ! dir=1..sdim  side=1..2 veldir=1..sdim
subroutine AHMED_ICE_RESISTANT_VEL_BC(xwall,xghost,t,LS, &
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

 call AHMED_ICE_RESISTANT_VEL(xghost,t,LS,local_VEL,velsolid_flag,dx,nmat)
 VEL=local_VEL(veldir)

else
 print *,"dir,side, or veldir invalid"
 stop
endif

return
end subroutine AHMED_ICE_RESISTANT_VEL_BC

 ! dir=1..sdim  side=1..2
subroutine AHMED_ICE_RESISTANT_PRES_BC(xwall,xghost,t,LS, &
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

 call AHMED_ICE_RESISTANT_PRES(xghost,t,LS,PRES,nmat)

else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine AHMED_ICE_RESISTANT_PRES_BC

 ! dir=1..sdim  side=1..2
subroutine AHMED_ICE_RESISTANT_STATE_BC(xwall,xghost,t,LS, &
   STATE,STATE_merge,STATE_in,im,istate,dir,side,dx,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: xwall
REAL_T, intent(in) :: xghost(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(nmat)
REAL_T :: local_STATE(nmat*num_state_material)
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
 call AHMED_ICE_RESISTANT_STATE(xghost,t,LS,local_STATE,local_bcflag, &
         nmat,num_state_material)
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
end subroutine AHMED_ICE_RESISTANT_STATE_BC

subroutine AHMED_ICE_RESISTANT_HEATSOURCE(im,VFRAC,time,x, &
     xsten,nhalf,temp, &
     heat_source,den,CV,dt,nmat)
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
INTEGER_T :: expected_nmat

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

if (axis_dir.eq.0) then
 expected_nmat=4
else if (axis_dir.eq.1) then
 expected_nmat=3
else if (axis_dir.eq.2) then
 expected_nmat=4
else
 print *,"axis_dir invalid"
 stop
endif

if ((num_materials.eq.expected_nmat).and. &
    (probtype.eq.425)) then
 heat_source=zero
else
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine AHMED_ICE_RESISTANT_HEATSOURCE

end module AHMED_ICE_RESISTANT_module
