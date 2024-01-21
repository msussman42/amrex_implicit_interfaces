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

! probtype==82
module ROTATING_ANNULUS_module
use amrex_fort_module, only : amrex_real

implicit none 

real(amrex_real), PARAMETER :: dT_dr=0.0d0
real(amrex_real), PARAMETER :: dT_dz=0.0d0

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
integer, INTENT(out) :: internal_wave_exists

internal_wave_exists=1

return
end subroutine ROTATING_ANNULUS_INTERNAL_GRAVITY_WAVE_FLAG

 ! fluids tessellate the domain, solids are immersed. 
subroutine ROTATING_ANNULUS_LS(x,t,LS,nmat)
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

if (levelrz.eq.COORDSYS_CYLINDRICAL) then
 !do nothing
else if (levelrz.eq.COORDSYS_CARTESIAN) then
 !do nothing
else
 print *,"expecting levelrz=COORDSYS_CYLINDRICAL or "
 print *,"levelrz=COORDSYS_CARTESIAN "
 stop
endif

do dir=1,SDIM
 VEL(dir)=zero
enddo

if (probtype.eq.82) then
 call ROTATING_ANNULUS_V0_Coriolis(x,dx,t,VEL)
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

integer, INTENT(in) :: bcflag !0=called from initialize  1=called from bc
integer, INTENT(in) :: nmat
integer, INTENT(in) :: nstate_mat
real(amrex_real), INTENT(in) :: x(SDIM)
integer :: radial_dir
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: LS(nmat)
real(amrex_real), INTENT(out) :: STATE(nmat*nstate_mat)
integer im,ibase,n
real(amrex_real) dT_dr_local
real(amrex_real) T0
real(amrex_real) dx_local(SDIM)
integer local_dir
integer, PARAMETER :: im_liquid=1

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

if (levelrz.eq.COORDSYS_CYLINDRICAL) then
 radial_dir=1
else if (levelrz.eq.COORDSYS_CARTESIAN) then
 radial_dir=2
else
 print *,"expecting levelrz=COORDSYS_CYLINDRICAL or "
 print *,"levelrz=COORDSYS_CARTESIAN "
 stop
endif

if ((problenx.eq.problen_array(1)).and. &
    (probleny.eq.problen_array(2)).and. &
    (problenx.gt.zero).and. &
    (probleny.gt.zero).and. &
    (problo_array(radial_dir).gt.zero).and. &
    (x(radial_dir).gt.zero)) then
 ! do nothing
else
 print *,"problenx or problo_array or x(radial_dir) invalid"
 stop
endif

dT_dr_local=(twall-fort_tempconst(1))/problen_array(radial_dir)

if (twall.ge.fort_tempconst(1)) then
 if (abs(dT_dr-dT_dr_local).le.EPS2) then
  ! do nothig
 else
  print *,"TA (fort_tempconst(1) or TB (twall) invalid"
  print *,"(abs(dT_dr-dT_dr_local).le.EPS2) failed"
  print *,"twall=",twall
  print *,"fort_tempconst(1)=",fort_tempconst(1)
  stop
 endif
else
 print *,"TA (fort_tempconst(1) or TB (twall) invalid"
 print *,"twall=",twall
 print *,"fort_tempconst(1)=",fort_tempconst(1)
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

  do local_dir=1,SDIM
   dx_local(local_dir)=one
  enddo

  call ROTATING_ANNULUS_T0_Boussinesq(x,dx_local,t,im_liquid,T0)

  STATE(ibase+ENUM_TEMPERATUREVAR+1)=T0

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
 VEL=local_VEL(veldir)

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

do local_dir=1,SDIM
 xwall_vec(local_dir)=xghost(local_dir)
enddo
xwall_vec(dir)=xwall

if (levelrz.eq.COORDSYS_CYLINDRICAL) then
 radial_dir=1
else if (levelrz.eq.COORDSYS_CARTESIAN) then
 radial_dir=2
else
 print *,"expecting levelrz=COORDSYS_CYLINDRICAL or "
 print *,"levelrz=COORDSYS_CARTESIAN "
 stop
endif

if ((problenx.eq.problen_array(1)).and. &
    (probleny.eq.problen_array(2)).and. &
    (problenx.gt.zero).and. &
    (probleny.gt.zero).and. &
    (problo_array(radial_dir).gt.zero).and. &
    (xghost(radial_dir).gt.zero)) then
 ! do nothing
else
 print *,"problenx or problo_array or xghost(radial_dir) invalid"
 stop
endif
if (xwall_vec(radial_dir).gt.zero) then
 ! do nothing
else
 print *,"xwall_vec(radial_dir) invalid"
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

if (probtype.eq.82) then
 ! do nothing
else
 print *,"expecting probtype.eq.82"
 stop
endif

call get_primary_material(LS,im_crit)
if (im_crit.eq.1) then
 ! do nothing
else
 print *,"expecting im_crit=1"
 stop
endif

! fort_tempconst(1) is the inner wall temperature
! twall is the outer wall temperature

if ((istate.ge.1).and. &
    (istate.le.num_state_material).and. &
    (im.ge.1).and. &
    (im.le.num_materials)) then

 call ROTATING_ANNULUS_STATE(xghost,t,LS,local_STATE,local_bcflag, &
         nmat,num_state_material)
 ibase=(im-1)*num_state_material
 STATE=local_STATE(ibase+istate)
 ibase=(im_crit-1)*num_state_material
 STATE_merge=local_STATE(ibase+istate)

 if (istate.eq.1+ENUM_TEMPERATUREVAR) then
  call ROTATING_ANNULUS_STATE(xwall_vec,t,LS,local_STATE,local_bcflag, &
         nmat,num_state_material)
  ibase=(im-1)*num_state_material
  STATE=local_STATE(ibase+istate)
  ibase=(im_crit-1)*num_state_material
  STATE_merge=local_STATE(ibase+istate)

 else if (istate.eq.1+ENUM_DENVAR) then
  ! do nothing
 else if (istate.gt.1+ENUM_TEMPERATUREVAR) then
  ! do nothing
 else
  print *,"istate invalid"
  stop
 endif

else
 print *,"istate invalid"
 stop
endif

return
end subroutine ROTATING_ANNULUS_STATE_BC

! returns (1/w) where w>>1 in "trouble" regions
subroutine ROTATING_ANNULUS_MAPPING_WEIGHT_COEFF(dir,wt,phys_x)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: dir !0,1,2
real(amrex_real), INTENT(out) :: wt
real(amrex_real), INTENT(in) :: phys_x
real(amrex_real) :: scaling
real(amrex_real) :: mid_x
integer :: radial_dir ! 1 or 2
integer :: azimuthal_dir ! 1 or 2

if (SDIM.eq.3) then
 ! do nothing
else
 print *,"expecting SDIM==3"
 stop
endif

if (levelrz.eq.COORDSYS_CYLINDRICAL) then
 radial_dir=1
 azimuthal_dir=2
else if (levelrz.eq.COORDSYS_CARTESIAN) then
 radial_dir=2
 azimuthal_dir=1
else
 print *,"expecting levelrz=COORDSYS_CYLINDRICAL or "
 print *,"levelrz=COORDSYS_CARTESIAN "
 stop
endif

if ((dir.ge.0).and.(dir.lt.SDIM)) then
 ! do nothing
else
 print *,"dir invalid: ",dir
 stop
endif

if ((phys_x.ge.zero).or.(phys_x.le.zero)) then
 ! do nothing
else
 print *,"phys_x is NaN: ",phys_x
 stop
endif

if (1.eq.1) then
 wt=one
else if (1.eq.0) then
 if (dir+1.eq.radial_dir) then
  scaling=problen_array(radial_dir)
  mid_x=half*(problo_array(radial_dir)+probhi_array(radial_dir))

  if (phys_x.le.mid_x) then
   wt=one+(scaling/(problo_array(radial_dir)-phys_x))**2
  else if (phys_x.ge.mid_x) then
   wt=one+(scaling/(probhi_array(radial_dir)-phys_x))**2
  else
   print *,"phys_x bust"
   stop
  endif
 else if (dir+1.eq.azimuthal_dir) then
  wt=one
 else if ((dir.eq.SDIM).and.(SDIM.eq.3)) then
  scaling=(probhiz-probloz)
  if (phys_x.le.half*(probloz+probhiz)) then
   wt=one+(scaling/(probloz-phys_x))**2
  else if (phys_x.ge.half*(probloz+probhiz)) then
   wt=one+(scaling/(probhiz-phys_x))**2
  else
   print *,"phys_x bust"
   stop
  endif
 else 
  print *,"dir invalid: ",dir
  stop
 endif
else
 print *,"incorrect option"
 stop
endif

wt=one/wt

return
end subroutine ROTATING_ANNULUS_MAPPING_WEIGHT_COEFF

subroutine ROTATING_ANNULUS_T0_Boussinesq(x,dx,cur_time,im,T0)
use probcommon_module
use global_utility_module
IMPLICIT NONE

real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: dx(SDIM)
real(amrex_real), INTENT(in) :: cur_time
integer, INTENT(in) :: im
real(amrex_real), INTENT(out) :: T0

 if (cur_time.ge.0.0d0) then
  ! do nothing
 else
  print *,"cur_time invalid: ",cur_time
  stop
 endif
 if ((im.ge.1).and.(im.le.num_materials)) then
  ! do nothing
 else
  print *,"im invalid"
  stop
 endif
 if (SDIM.eq.3) then
  ! do nothing
 else
  print *,"expecting sdim=3"
  stop
 endif

 T0=fort_tempconst(im)
 if (levelrz.eq.COORDSYS_CYLINDRICAL) then
  T0=T0+dT_dr*x(1)
 else if (levelrz.eq.COORDSYS_CARTESIAN) then
  T0=T0+dT_dr*x(2)
 else
  print *,"expecting levelrz=COORDSYS_CYLINDRICAL or "
  print *,"levelrz=COORDSYS_CARTESIAN "
  stop
 endif
 T0=T0+dT_dz*x(SDIM)
        
end subroutine ROTATING_ANNULUS_T0_Boussinesq

subroutine ROTATING_ANNULUS_V0_Coriolis(x,dx,cur_time,V0)
use probcommon_module
use global_utility_module
IMPLICIT NONE

real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: dx(SDIM)
real(amrex_real), INTENT(in) :: cur_time
real(amrex_real), INTENT(out) :: V0(SDIM)

integer :: dir
real(amrex_real) :: dV_dz

 if (cur_time.ge.0.0d0) then
  ! do nothing
 else
  print *,"cur_time invalid: ",cur_time
  stop
 endif
 if (SDIM.eq.3) then
  ! do nothing
 else
  print *,"expecting sdim=3"
  stop
 endif
 if (fort_angular_velocity.gt.0.0d0) then
  ! do nothing
 else
  print *,"expecting fort_angular_velocity>0.0d0"
  stop
 endif

 do dir=1,SDIM
  V0(dir)=0.0d0
 enddo

  ! rho=rho0*(1+fort_DrhoDT(im)*(T-T0))
 dV_dz=dT_dr*abs(gravity_vector(SDIM)*fort_DrhoDT(1))/ &
             (two*fort_angular_velocity)

 if (levelrz.eq.COORDSYS_CYLINDRICAL) then
  V0(2)=dV_dz*x(SDIM)
 else if (levelrz.eq.COORDSYS_CARTESIAN) then
  V0(1)=-dV_dz*x(SDIM)
 else
  print *,"expecting levelrz=COORDSYS_CYLINDRICAL or "
  print *,"levelrz=COORDSYS_CARTESIAN "
  stop
 endif

end subroutine ROTATING_ANNULUS_V0_Coriolis



end module ROTATING_ANNULUS_module
