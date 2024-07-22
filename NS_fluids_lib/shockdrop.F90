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

#define shockdrop_PROB_TYPE 3001

#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else
print *,"dimension bust"
stop
#endif

MODULE shockdrop
use amrex_fort_module, only : amrex_real
implicit none 

 !see run2d/NormalShockWaveNASA.F90
 !see run2d/inputs.shockdrop

real(amrex_real) shockdrop_R
real(amrex_real) shockdrop_cp
real(amrex_real) shockdrop_cv

real(amrex_real) shockdrop_We
real(amrex_real) shockdrop_Oh
real(amrex_real) shockdrop_Re
real(amrex_real) shockdrop_inertial_time_scale

 !upstream supersonic relative to shock
real(amrex_real) shockdrop_P0
real(amrex_real) shockdrop_T0
real(amrex_real) shockdrop_DEN0
real(amrex_real) shockdrop_M0
real(amrex_real) shockdrop_VEL0
real(amrex_real) shockdrop_C0
real(amrex_real) shockdrop_EE0

 !downstream subsonic relative to shock
real(amrex_real) shockdrop_P1
real(amrex_real) shockdrop_T1
real(amrex_real) shockdrop_DEN1
real(amrex_real) shockdrop_M1
real(amrex_real) shockdrop_VEL1
real(amrex_real) shockdrop_C1
real(amrex_real) shockdrop_EE1

real(amrex_real) shockdrop_gamma

CONTAINS

! units are cgs
subroutine shockdrop_init()
use probcommon_module
use global_utility_module
IMPLICIT NONE
real(amrex_real) test_pres

 !material_type=5 EOS_air
 !see PROBCOMMON.F90
shockdrop_R=R_AIR_PARMS !ergs/(Kelvin g)
shockdrop_cv=CV_AIR_PARMS !ergs/(Kelvin g)
shockdrop_cp=shockdrop_cv+shockdrop_R !ergs/(Kelvin g)

! shockdrop_M0=1.4
! shockdrop_M0=3.0
! shockdrop_M0=1.17
shockdrop_M0=vinletgas

if (shockdrop_M0.gt.one) then
 !do nothing
else
 print *,"shockdrop_M0 invalid: ",shockdrop_M0
 stop
endif

shockdrop_T0=278.0d0 !tempconst(2)
shockdrop_DEN0=0.00125335272d0 !denconst(2)

if (abs(shockdrop_DEN0-fort_denconst(2)).le.EPS6) then
 !do nothing
else
 print *,"shockdrop_DEN0 ",shockdrop_DEN0
 print *,"fort_denconst(2) ",fort_denconst(2)
 print *,"mismatch"
 stop
endif

if (abs(shockdrop_T0-fort_tempconst(2)).le.EPS6) then
 !do nothing
else
 print *,"shockdrop_T0 ",shockdrop_T0
 print *,"fort_tempconst(2) ",fort_tempconst(2)
 print *,"mismatch"
 stop
endif
if (abs(shockdrop_T0-fort_tempconst(1)).le.EPS6) then
 !do nothing
else
 print *,"shockdrop_T0 ",shockdrop_T0
 print *,"fort_tempconst(1) ",fort_tempconst(1)
 print *,"mismatch"
 stop
endif


 ! this value must be consistent with material_type=5 parameters
shockdrop_gamma=shockdrop_cp/shockdrop_cv

 ! if compressible liquid, this value should be the same
 ! as the equilibrium pressure of the liquid drop.
shockdrop_P0=(shockdrop_gamma-1.0d0)* &
 shockdrop_DEN0*shockdrop_cv*shockdrop_T0

shockdrop_C0=(shockdrop_gamma*shockdrop_P0/shockdrop_DEN0)**half
shockdrop_VEL0=shockdrop_M0*shockdrop_C0

shockdrop_P1=shockdrop_P0* &
 (two*shockdrop_gamma*(shockdrop_M0**2)-shockdrop_gamma+one)/ &
 (shockdrop_gamma+one)

shockdrop_T1=shockdrop_T0* &
 ( two*shockdrop_gamma*(shockdrop_M0**2)-shockdrop_gamma+one )* &
 ( (shockdrop_gamma-one)*(shockdrop_M0**2)+two )/ &
 ( ((shockdrop_gamma+one)*shockdrop_M0)**2 )

shockdrop_DEN1=shockdrop_DEN0* &
 ((shockdrop_gamma+one)*(shockdrop_M0**2))/ &
 ((shockdrop_gamma-one)*(shockdrop_M0**2)+two)

shockdrop_C1=sqrt(shockdrop_gamma*shockdrop_P1/shockdrop_DEN1)
shockdrop_M1= &
 ((shockdrop_gamma-one)*(shockdrop_M0**2)+two)/ &
 (two*shockdrop_gamma*(shockdrop_M0**2)-shockdrop_gamma+one)
shockdrop_M1=sqrt(shockdrop_M1)
shockdrop_VEL1=shockdrop_M1*shockdrop_C1

!shock_vel0-shock_vel1=M0 * C0 - M1 * C1

call general_hydrostatic_pressure(test_pres)
if (abs(test_pres-shockdrop_P0)/test_pres.gt.1.0E-8) then
 print *,"shockdrop_P0 inconsistent w/general_hydrostatic_pressure"
 print *,"test_pres=",test_pres
 print *,"shockdrop_P0=",shockdrop_P0
 stop
endif
if (fort_material_type(2).ne.5) then
 print *,"only material_type=5 supported for gas for this problem"
 print *,"fort_material_type(2)=",fort_material_type(2)
 stop
endif

print *,"shockdrop: upstream den,approaching SPEED,T,M,C ", &
 shockdrop_DEN0,shockdrop_VEL0,shockdrop_T0,shockdrop_M0,shockdrop_C0
print *,"shockdrop: downstream den,SPEED,T,M,C ", &
 shockdrop_DEN1,shockdrop_VEL1,shockdrop_T1,shockdrop_M1,shockdrop_C1

! in shock frame of reference, the upstream velocity is -|V0| and the
! downstream velocity is -|V1|.  In upstream frame of reference the 
! downstream velocity is |V0|-|V1|

if (shockdrop_VEL0.gt.shockdrop_VEL1) then
 shockdrop_We=shockdrop_DEN1* &
   ((shockdrop_VEL0-shockdrop_VEL1)**2)*two*radblob/ &
   fort_tension(1)
 shockdrop_Oh=fort_viscconst(1)/ &
   sqrt(fort_denconst(1)*fort_tension(1)*two*radblob)
 shockdrop_Re=shockdrop_DEN1*(shockdrop_VEL0-shockdrop_VEL1)* &
   two*radblob/fort_viscconst(2)
 shockdrop_inertial_time_scale=two*radblob* &
    sqrt(fort_denconst(1)/shockdrop_DEN1)/ &
    (shockdrop_VEL0-shockdrop_VEL1) 
 print *,"shockdrop_We=",shockdrop_We
 print *,"shockdrop_Oh=",shockdrop_Oh
 print *,"shockdrop_Re=",shockdrop_Re
 print *,"shockdrop_inertial_time_scale=",shockdrop_inertial_time_scale
else
 print *,"shockdrop_VEL0 or shockdrop_VEL1 invalid: ", &
   shockdrop_VEL0,shockdrop_VEL1
 stop
endif

return
end subroutine shockdrop_init

subroutine shockdrop_velocity(x,t,ls,vel, &
  velsolid_flag,dx,nmat)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: dx(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: ls(nmat)
real(amrex_real) :: ls_local
real(amrex_real), INTENT(out) :: vel(SDIM)
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
 print *,"velsolid_flag invalid: ",velsolid_flag
 stop
endif
do dir=1,SDIM
 if (dx(dir).gt.zero) then
  ! do nothing
 else
  print *,"dx invalid: ",dx(dir)
  stop
 endif
enddo

do dir=1,SDIM
 vel(dir)=zero
enddo

call shockdrop_dropLS(x(1),x(2),x(SDIM),ls_local)
if (ls_local.ge.zero) then
 ! do nothing in y/z direction (drop is upstream from shock and
 ! stationary in the "upstream frame of reference")
 ! shock velocity > 0
 ! shock is approaching with speed: shockdrop_VEL0
 vel(1)=advbot ! velocity in the "periodic" direction
else
 call shockdrop_shockLS(x(1),x(2),x(SDIM),ls_local)
 ! in shock frame of reference:
 ! upstream: v=-shockdrop_VEL0
 ! downstream: v=-shockdrop_VEL1
 if (ls_local.ge.zero) then  ! upstream (above the shock)
  vel(SDIM)=zero
 else
  vel(SDIM)=shockdrop_VEL0-shockdrop_VEL1
 endif
endif 

return
end subroutine shockdrop_velocity

 !gets mapped to SUB_CFL_HELPER
subroutine shockdrop_maxvelocity(time,dir,vel,dx)
use probcommon_module
IMPLICIT NONE
integer, INTENT(in) :: dir
real(amrex_real), INTENT(in) :: time
real(amrex_real), INTENT(inout) ::  vel
real(amrex_real) :: vel_local
real(amrex_real), INTENT(in) :: dx(SDIM)

if ((dir.lt.0).or.(dir.ge.SDIM)) then
 print *,"dir invalid: ",dir
 stop
endif
vel_local=max(abs(shockdrop_VEL0),abs(shockdrop_VEL1))
vel_local=max(abs(advbot),abs(vel_local))
vel=max(abs(vel),abs(vel_local))

return
end subroutine shockdrop_maxvelocity


subroutine shockdrop_LS(x,t,LS,nmat)
use probcommon_module
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

call shockdrop_dropLS(x(1),x(2),x(SDIM),LS(1))
LS(2)=-LS(1)

return
end subroutine shockdrop_LS

subroutine shockdrop_pressure(x,t,ls,pres,nmat)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: ls(nmat)
real(amrex_real) :: ls_local
real(amrex_real), INTENT(out) :: pres


if (num_materials.eq.nmat) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

call shockdrop_dropLS(x(1),x(2),x(SDIM),ls_local)
if (ls_local.ge.zero) then ! liquid
 pres=shockdrop_P0
else
 call shockdrop_shockLS(x(1),x(2),x(SDIM),ls_local)
 if (ls_local.ge.zero) then  ! upstream (above the approaching shock)
  pres=shockdrop_P0
 else
  pres=shockdrop_P1
 endif
endif 

return
end subroutine shockdrop_pressure


subroutine shockdrop_gas_density(x,y,z,den)
use probcommon_module
IMPLICIT NONE
real(amrex_real), intent(in) :: x,y,z
real(amrex_real), intent(out) :: den
real(amrex_real) LS

if (SDIM.eq.2) then
 if (abs(z-y).le.1.0E-6) then
  !do nothing
 else
  print *,"expecting z=y"
  stop
 endif
endif

call shockdrop_dropLS(x,y,z,LS)
if (LS.ge.zero) then ! liquid
 den=shockdrop_DEN0
else
 call shockdrop_shockLS(x,y,z,LS)
 if (LS.ge.zero) then  ! upstream (above the shock)
  den=shockdrop_DEN0
 else
  den=shockdrop_DEN1
 endif
endif 

return
end subroutine shockdrop_gas_density


subroutine shockdrop_gas_temperature(x,y,z,temp)
use probcommon_module
IMPLICIT NONE
real(amrex_real), intent(in) :: x,y,z
real(amrex_real), intent(out) :: temp
real(amrex_real) LS

if (SDIM.eq.2) then
 if (abs(z-y).le.1.0E-6) then
  !do nothing
 else
  print *,"expecting z=y"
  stop
 endif
endif

call shockdrop_dropLS(x,y,z,LS)
if (LS.ge.zero) then ! liquid
 temp=shockdrop_T0
else
 call shockdrop_shockLS(x,y,z,LS)
 if (LS.ge.zero) then  ! upstream (above the shock)
  temp=shockdrop_T0
 else
  temp=shockdrop_T1
 endif
endif 

return
end subroutine shockdrop_gas_temperature


! probtype=1 in the inputs file
! axis_dir=150 shock drop
! axis_dir=151 shock column
! axis_dir=152 shock cylinder
! LS>0 upstream of the shock  z>zblob2
! LS<0 downstream of the shock z<zblob2
SUBROUTINE shockdrop_shockLS(x,y,z,LS)
use probcommon_module
IMPLICIT NONE
real(amrex_real),INTENT(in) :: x,y,z
real(amrex_real),INTENT(out) :: LS

if (SDIM.eq.2) then
 if (abs(z-y).le.1.0E-6) then
  !do nothing
 else
  print *,"z<>y error"
  stop
 endif
endif

if ((axis_dir.eq.150).or. &
    (axis_dir.eq.151).or. &
    (axis_dir.eq.152)) then
 LS=z-zblob2
else
 print *,"axis_dir invalid"
 stop
endif

return
END SUBROUTINE shockdrop_shockLS

! probtype=1 in the inputs file
! axis_dir=150 shock drop
! axis_dir=151 shock column
! axis_dir=152 shock cylinder
! LS>0 in the drop
subroutine shockdrop_dropLS(x,y,z,LS)
use probcommon_module
IMPLICIT NONE
real(amrex_real),INTENT(in) :: x,y,z
real(amrex_real),INTENT(out) :: LS
real(amrex_real) mag

if (SDIM.eq.2) then
 if (abs(z-y).le.1.0E-6) then
  !do nothing
 else
  print *,"z<>y error"
  stop
 endif
endif

if (axis_dir.eq.150) then ! shock drop
 mag=(x-xblob)**2+(y-yblob)**2
 if (SDIM.eq.3) then
  mag=mag+(z-zblob)**2
 endif
 LS=radblob-sqrt(mag) 
else if (axis_dir.eq.151) then ! shock column
 if (SDIM.eq.2) then
  mag=(y-yblob)**2
  LS=radblob-sqrt(mag) 
 else if (SDIM.eq.3) then
  mag=(z-zblob)**2
  LS=radblob-sqrt(mag) 
 else 
  print *,"dimension bust"
  stop
 endif
else if (axis_dir.eq.152) then ! shock cylinder
 if (SDIM.eq.2) then
  mag=(y-yblob)**2
  LS=radblob-sqrt(mag) 
 else if (SDIM.eq.3) then
  mag=(y-yblob)**2+(z-zblob)**2
  LS=radblob-sqrt(mag) 
 else 
  print *,"dimension bust"
  stop
 endif
else
 print *,"axis_dir invalid in shockdrop_dropLS"
 stop
endif

return
END SUBROUTINE shockdrop_dropLS


subroutine shockdrop_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: bcflag !0=called from initialize  1=called from bc
integer, INTENT(in) :: nmat
integer, INTENT(in) :: nstate_mat
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: LS(nmat)
real(amrex_real), INTENT(out) :: STATE(nmat*nstate_mat)
integer im,ibase,n

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

if ((num_materials.eq.2).and. &
    (num_state_material.ge.2).and. &
    (probtype.eq.shockdrop_PROB_TYPE)) then
 do im=1,num_materials
  ibase=(im-1)*num_state_material
  STATE(ibase+ENUM_DENVAR+1)=fort_denconst(im)

  if (im.eq.1) then
   !do nothing
  else if (im.eq.2) then
   call shockdrop_gas_density(x(1),x(2),x(SDIM),STATE(ibase+ENUM_DENVAR+1))
  else
   print *,"im out of range"
   stop
  endif

  if (t.eq.zero) then
   STATE(ibase+ENUM_TEMPERATUREVAR+1)=fort_initial_temperature(im)
  else if (t.gt.zero) then
   STATE(ibase+ENUM_TEMPERATUREVAR+1)=fort_tempconst(im)
  else
   print *,"t invalid"
   stop
  endif

  if (im.eq.1) then
   !do nothing
  else if (im.eq.2) then
   call shockdrop_gas_temperature(x(1),x(2),x(SDIM), &
     STATE(ibase+ENUM_TEMPERATUREVAR+1))
  else
   print *,"im out of range"
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
end subroutine shockdrop_STATE

 ! dir=1..sdim  side=1..2
subroutine shockdrop_LS_BC(xwall,xghost,t,LS, &
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
real(amrex_real), INTENT(in) :: dx(SDIM)

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then
 call shockdrop_LS(xghost,t,LS,nmat)
else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine shockdrop_LS_BC


 ! dir=1..sdim  side=1..2 veldir=1..sdim
subroutine shockdrop_VEL_BC(xwall,xghost,t,LS, &
   VEL,VEL_in,veldir,dir,side,dx,nmat)
use probcommon_module
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
velsolid_flag=0
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2).and. &
    (veldir.ge.1).and.(veldir.le.SDIM)) then

 call shockdrop_velocity(xghost,t,LS,local_VEL,velsolid_flag,dx,nmat)
 VEL=local_VEL(veldir)

else
 print *,"dir,side, or veldir invalid"
 stop
endif

return
end subroutine shockdrop_VEL_BC

 ! dir=1..sdim  side=1..2
subroutine shockdrop_PRES_BC(xwall,xghost,t,LS, &
   PRES,PRES_in,dir,side,dx,nmat)
use probcommon_module
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


if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then

 call shockdrop_pressure(xghost,t,LS,PRES,nmat)

else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine shockdrop_PRES_BC

 ! dir=1..sdim  side=1..2
subroutine shockdrop_STATE_BC(xwall,xghost,t,LS, &
   STATE,STATE_merge,STATE_in,im,istate,dir,side,dx,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: xwall
real(amrex_real), INTENT(in) :: xghost(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: LS(nmat)
real(amrex_real) local_STATE(nmat*num_state_material)
real(amrex_real), INTENT(inout) :: STATE
real(amrex_real), INTENT(inout) :: STATE_merge
real(amrex_real), INTENT(in) :: STATE_in
integer, INTENT(in) :: dir,side
real(amrex_real), INTENT(in) :: dx(SDIM)
integer, INTENT(in) :: istate,im
integer ibase,im_crit
integer local_bcflag

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
 call shockdrop_STATE(xghost,t,LS,local_STATE, &
         local_bcflag,nmat,num_state_material)
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
end subroutine shockdrop_STATE_BC

subroutine shockdrop_OVERRIDE_TAGFLAG( &
  i,j,k, &
  level,max_level, &
  snew_ptr,lsnew_ptr, &
  xsten,nhalf,time, &
  rflag,tagflag)
use amrex_fort_module, only : amrex_real
use probcommon_module
use global_utility_module
IMPLICIT NONE
integer, INTENT(in) :: i,j,k
integer, INTENT(in) :: level,max_level
integer, INTENT(in) :: nhalf
real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
real(amrex_real), INTENT(in) :: time
real(amrex_real), INTENT(inout) :: rflag
integer, INTENT(inout) :: tagflag
real(amrex_real), INTENT(in),pointer :: snew_ptr(D_DECL(:,:,:),:)
real(amrex_real), INTENT(in),pointer :: lsnew_ptr(D_DECL(:,:,:),:)
real(amrex_real), dimension(3) :: local_x
real(amrex_real), dimension(SDIM) :: local_delta
real(amrex_real) :: F_LIQUID,LS_LIQUID,LS_SHOCK,P_diff,mag
integer :: dir
integer :: ii,jj,kk

if (nhalf.lt.3) then
 print *,"nhalf invalid shock drop override tagflag"
 stop
endif
if ((level.ge.0).and.(level.lt.max_level)) then
 ! do nothing
else
 print *,"level and/or max_level invalid"
 print *,"level=",level
 print *,"max_level=",max_level
 stop
endif
do dir=1,SDIM
 local_x(dir)=xsten(0,dir)
enddo
local_x(3)=xsten(0,SDIM)
do dir=1,SDIM
 local_delta(dir)=xsten(1,dir)-xsten(-1,dir)
 if (local_delta(dir).gt.zero) then
  ! do nothing
 else
  print *,"local_delta invalid shockdrop_override_tagflag"
  stop
 endif
enddo !dir=1..sdim

if ((num_materials.ge.2).and. &
    (probtype.eq.shockdrop_PROB_TYPE)) then

 if ((axis_dir.ge.150).and.(axis_dir.le.152)) then

  ii=0
  jj=0
  kk=0
  if (AMREX_SPACEDIM.eq.2) then
   jj=1
  else if (AMREX_SPACEDIM.eq.3) then
   kk=1
  else
   print *,"dimension problem"
   stop
  endif

  rflag=0.0d0
  tagflag=0
  call shockdrop_shockLS(local_x(1),local_x(2),local_x(SDIM),LS_SHOCK)
  LS_LIQUID=lsnew_ptr(D_DECL(i,j,k),1)
  F_LIQUID=snew_ptr(D_DECL(i,j,k),STATECOMP_MOF+1)
  if (time.eq.zero) then
   if (abs(LS_SHOCK).le.two*local_delta(SDIM)) then
    p_diff=abs((shockdrop_P1-shockdrop_P0)/shockdrop_P0)
   else if (abs(LS_SHOCK).ge.two*local_delta(SDIM)) then
    P_diff=zero
   else
    print *,"LS_SHOCK invalid: ",LS_SHOCK
    stop
   endif
  else if (time.gt.zero) then
   P_diff=abs((snew_ptr(D_DECL(i+ii,j+jj,k+kk),STATECOMP_PRES+1)- &
               snew_ptr(D_DECL(i,j,k),STATECOMP_PRES+1))/shockdrop_P0)
  else
   print *,"time invalid: ",time
   stop
  endif

  mag=(local_x(1)-xblob)**2+(local_x(2)-yblob)**2
  if (AMREX_SPACEDIM.eq.3) then
   mag=mag+(local_x(SDIM)-zblob)**2
  endif
  mag=sqrt(mag)
  if (mag.gt.five*radblob) then
   !do nothing
  else if (mag.ge.zero) then
   if ((LS_LIQUID.ge.zero).or.(F_LIQUID.ge.0.1d0)) then
    rflag=1.0d0
    tagflag=1
   endif
   if (P_diff.ge.0.01d0) then
    rflag=1.0d0
    tagflag=1
   endif
  else
   print *,"mag invalid: ",mag
   stop
  endif

 else
  print *,"axis_dir invalid: ",axis_dir
  stop
 endif

else
 print *,"num_materials or probtype invalid"
 print *,"num_materials: ",num_materials
 print *,"probtype: ",probtype
 stop
endif

end subroutine shockdrop_OVERRIDE_TAGFLAG

END MODULE shockdrop
