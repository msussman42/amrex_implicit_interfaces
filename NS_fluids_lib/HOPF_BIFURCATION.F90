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

! probtype==820
! driven cavity problem
module HOPF_BIFURCATION_module
use amrex_fort_module, only : amrex_real

implicit none 

contains

  ! do any initial preparation needed
subroutine INIT_HOPF_BIFURCATION_MODULE()
IMPLICIT NONE

return
end subroutine INIT_HOPF_BIFURCATION_MODULE

! referring to:
! Hopf Bifurcation of the
! Unsteady Regularized Driven Cavity Flow
!  JIE SHEN
! Density is 1 everywhere, and there are no Boussinesq type approximations,
! so we set internal_wave_exists=0
subroutine HOPF_BIFURCATION_INTERNAL_GRAVITY_WAVE_FLAG(internal_wave_exists)
use probcommon_module
use global_utility_module
IMPLICIT NONE
integer, INTENT(out) :: internal_wave_exists

internal_wave_exists=0

return
end subroutine HOPF_BIFURCATION_INTERNAL_GRAVITY_WAVE_FLAG

 ! fluids tessellate the domain, solids are immersed. 
subroutine HOPF_BIFURCATION_LS(x,t,LS,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(out) :: LS(nmat)

  if (nmat.eq.num_materials) then
   ! do nothing
  else
   print *,"nmat invalid"
   stop
  endif
  if (num_materials.eq.2) then 
   ! do nothing
  else
   print *,"expecting num_materials=2 for HOPF_BIFURCATION"
   stop
  endif

  LS(1)=99999.90
  LS(2)=-LS(1)

return
end subroutine HOPF_BIFURCATION_LS

! initial velocity is zero
subroutine HOPF_BIFURCATION_VEL(x,t,LS,VEL,velsolid_flag,dx,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: dx(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: LS(nmat)
real(amrex_real), INTENT(out) :: VEL(SDIM)
integer dir
integer, INTENT(in) :: velsolid_flag

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

if (levelrz.eq.COORDSYS_CARTESIAN) then
 !do nothing
else
 print *,"expecting levelrz=COORDSYS_CARTESIAN "
 stop
endif

do dir=1,SDIM
 VEL(dir)=zero
enddo

!u=(16 x^2 (1-x)^2 , 0) upper lid
!u=(0,0) all other walls.
!u=(0,0) t=0 interior of domain.
if (probtype.eq.820) then
 if ((problenx.eq.one).and. &
     (probleny.eq.one).and. &
     (problox.eq.zero).and. &
     (probloy.eq.zero).and. &
     (probhix.eq.one).and. &
     (probhiy.eq.one)) then
  if (x(1).le.zero) then
   !do nothing
  else if (x(1).ge.one) then
   !do nothing
  else if (x(2).lt.one) then
   ! do nothing
  else if (x(2).ge.one) then
   !u=(16 x^2 (1-x)^2 , 0) upper lid
   VEL(1)=16.0d0*(x(1)**2)*((one-x(1))**2)
  else
   print *,"x() is NaN"
   stop
  endif 
 else
  print *,"problen, etc invalid"
  stop
 endif
else
 print *,"expecting probtype==820"
 stop
endif

return 
end subroutine HOPF_BIFURCATION_VEL


! this routine used as a default when
! pressure boundary conditions are prescribed.
subroutine HOPF_BIFURCATION_PRES(x,t,LS,PRES,nmat)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: LS(nmat)
real(amrex_real), INTENT(out) :: PRES

if (num_materials.eq.nmat) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
if (probtype.eq.820) then
 ! do nothing
else
 print *,"expecting probtype==820"
 stop
endif
PRES=zero

return 
end subroutine HOPF_BIFURCATION_PRES

subroutine HOPF_BIFURCATION_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: bcflag !0=called from initialize  1=called from bc
integer, INTENT(in) :: nmat
integer, INTENT(in) :: nstate_mat
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: LS(nmat)
real(amrex_real), INTENT(out) :: STATE(nmat*nstate_mat)
integer im,ibase,n
real(amrex_real) dT_dr_local
real(amrex_real) T0
real(amrex_real) dx_local(SDIM)
integer local_dir

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

if ((problenx.eq.problen_array(1)).and. &
    (probleny.eq.problen_array(2)).and. &
    (problenx.gt.zero).and. &
    (probleny.gt.zero)) then
 ! do nothing
else
 print *,"problenx or probleny invalid"
 stop
endif

if (probtype.eq.820) then

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

  do local_dir=1,SDIM
   dx_local(local_dir)=one
  enddo

  do n=1,num_species_var
   STATE(ibase+ENUM_SPECIESVAR+n)=fort_speciesconst((n-1)*num_materials+im)
  enddo

 enddo ! im=1..num_materials

else
 print *,"probtype invalid"
 print *,"aborting HOPF_BIFURCATION_STATE"
 stop
endif
 
return
end subroutine HOPF_BIFURCATION_STATE

 ! dir=1..sdim  side=1..2
subroutine HOPF_BIFURCATION_LS_BC(xwall,xghost,t,LS, &
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
real(amrex_real), INTENT(in) ::  dx(SDIM)

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then
 call HOPF_BIFURCATION_LS(xghost,t,LS,nmat)
else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine HOPF_BIFURCATION_LS_BC


 ! dir=1..sdim  side=1..2 veldir=1..sdim
subroutine HOPF_BIFURCATION_VEL_BC(xwall,xghost,t,LS, &
   VEL,VEL_in,veldir,dir,side,dx,nmat)
use probcommon_module
use global_utility_module
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
if (probtype.eq.820) then
 ! do nothing
else
 print *,"expecting probtype==820"
 stop
endif

velsolid_flag=0

if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2).and. &
    (veldir.ge.1).and.(veldir.le.SDIM)) then

 call HOPF_BIFURCATION_VEL(xghost,t,LS,local_VEL,velsolid_flag,dx,nmat)
 VEL=local_VEL(veldir)

else
 print *,"dir,side, or veldir invalid"
 stop
endif

return
end subroutine HOPF_BIFURCATION_VEL_BC

 ! dir=1..sdim  side=1..2
subroutine HOPF_BIFURCATION_PRES_BC(xwall,xghost,t,LS, &
   PRES,PRES_in,dir,side,dx,nmat)
use probcommon_module
use global_utility_module
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
if (probtype.eq.820) then
 ! do nothing
else
 print *,"expecting probtype==820"
 stop
endif

if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then

 call HOPF_BIFURCATION_PRES(xghost,t,LS,PRES,nmat)

else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine HOPF_BIFURCATION_PRES_BC

 ! dir=1..sdim  side=1..2
subroutine HOPF_BIFURCATION_STATE_BC(xwall,xghost,t,LS, &
   STATE,STATE_merge,STATE_in,im,istate,dir,side,dx,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: xwall
real(amrex_real), INTENT(in) :: xghost(SDIM)
real(amrex_real) :: xwall_vec(SDIM)
integer :: radial_dir ! 1 or 2
integer :: local_dir
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: LS(nmat)
real(amrex_real) :: local_STATE(nmat*num_state_material)
real(amrex_real), INTENT(inout) :: STATE
real(amrex_real), INTENT(inout) :: STATE_merge
real(amrex_real), INTENT(in) :: STATE_in
integer, INTENT(in) :: dir,side
real(amrex_real), INTENT(in) :: dx(SDIM)
integer, INTENT(in) :: istate,im
integer ibase,im_crit
integer, PARAMETER :: local_bcflag=1

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

if ((problenx.eq.problen_array(1)).and. &
    (probleny.eq.problen_array(2))) then
 ! do nothing
else
 print *,"problenx or probleny"
 stop
endif
if (dir.eq.1) then
 ! do nothing
else if (dir.eq.2) then
 ! do nothing
else if ((dir.eq.3).and.(SDIM.eq.3)) then
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

if (probtype.eq.820) then
 ! do nothing
else
 print *,"expecting probtype.eq.820"
 stop
endif

call get_primary_material(LS,im_crit)
if (im_crit.eq.1) then
 ! do nothing
else
 print *,"expecting im_crit=1"
 stop
endif

if ((istate.ge.1).and. &
    (istate.le.num_state_material).and. &
    (im.ge.1).and. &
    (im.le.num_materials)) then

 call HOPF_BIFURCATION_STATE(xghost,t,LS,local_STATE,local_bcflag, &
         nmat,num_state_material)
 ibase=(im-1)*num_state_material
 STATE=local_STATE(ibase+istate)
 ibase=(im_crit-1)*num_state_material
 STATE_merge=local_STATE(ibase+istate)

else
 print *,"istate invalid"
 stop
endif

return
end subroutine HOPF_BIFURCATION_STATE_BC

end module HOPF_BIFURCATION_module
