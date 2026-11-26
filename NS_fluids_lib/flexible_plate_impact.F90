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
use amrex_fort_module, only : amrex_real

implicit none 

real(amrex_real) :: DEF_VAPOR_GAMMA

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
real(amrex_real), INTENT(in), dimension(SDIM) :: x !spatial coordinates
real(amrex_real), INTENT(out) :: Phi !LS dist, Phi>0 in the substrate
integer, parameter :: im_solid=3

 if (num_materials.eq.im_solid) then

  if (FSI_flag(im_solid).eq.FSI_SHOELE_CTML) then
   !CTML takes care of this.
   Phi=-99999.0
  else if ((FSI_flag(im_solid).eq.FSI_EULERIAN_ELASTIC).or. &
           (FSI_flag(im_solid).eq.FSI_PRESCRIBED_PROBF90).or. &
           (FSI_flag(im_solid).eq.FSI_ICE_EULERIAN_ELASTIC).or. &
           (FSI_flag(im_solid).eq.FSI_ICE_STATIC)) then

   if (FSI_flag(im_solid).eq.FSI_EULERIAN_ELASTIC) then
    if (is_elastic(im_solid).eq.1) then
     !do nothing
    else
     print *,"expecting is_elastic(im_solid)=1"
     stop
    endif
   else if (FSI_flag(im_solid).eq.FSI_ICE_EULERIAN_ELASTIC) then
    if (is_elastic(im_solid).eq.1) then
     !do nothing
    else
     print *,"expecting is_elastic(im_solid)=1"
     stop
    endif
   else if (FSI_flag(im_solid).eq.FSI_ICE_STATIC) then
    if (is_elastic(im_solid).eq.1) then
     !do nothing
    else
     print *,"expecting is_elastic(im_solid)=1"
     stop
    endif
   else if (FSI_flag(im_solid).eq.FSI_PRESCRIBED_PROBF90) then
    if (is_elastic(im_solid).eq.0) then
     !do nothing
    else
     print *,"expecting is_elastic(im_solid)=0"
     stop
    endif
    if (is_rigid(im_solid).eq.1) then
     !do nothing
    else
     print *,"expecting is_rigid(im_solid)=1"
     stop
    endif
   else
    print *,"FSI_flag invalid: ",im_solid,FSI_flag(im_solid)
    stop
   endif

   if (AMREX_SPACEDIM.eq.2) then
      !squaredist returns Phi<0 in the square.
    call squaredist(x(1),x(2), &
      xblob2-radblob2, &
      xblob2+radblob2, &
      yblob2-radblob3, &
      yblob2+radblob3, &
      Phi)
    Phi=-Phi
   else if (AMREX_SPACEDIM.eq.3) then
      !cylinderdist returns Phi<0 in the cylinder
    call cylinderdist( &
      x(1),x(2),x(SDIM), &
      xblob2,yblob2, &
      radblob2, &
      zblob2-radblob3, &
      zblob2+radblob3, &
      Phi)
    Phi=-Phi
   else
    print *,"dimension bust"
    stop
   endif
  else
   print *,"FSI_flag invalid: ",im_solid,FSI_flag(im_solid)
   stop
  endif
 else
  print *,"num_materials invalid: ",num_materials
  stop
 endif

end subroutine flexible_substrateLS


 ! fluids tessellate the domain, solids are immersed. 
subroutine flexible_plate_impact_LS(x,t,LS,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

  integer, INTENT(in) :: nmat
  real(amrex_real), INTENT(in) :: x(SDIM)
  real(amrex_real), INTENT(in) :: t
  real(amrex_real), INTENT(out) :: LS(nmat)
  integer, parameter :: im_solid=3

  if (nmat.eq.num_materials) then
   ! do nothing
  else
   print *,"nmat invalid"
   stop
  endif

if ((num_materials.eq.im_solid).and. &
    (probtype.eq.2000)) then

 if (FSI_flag(im_solid).eq.FSI_SHOELE_CTML) then
  !do nothing
 else if ((FSI_flag(im_solid).eq.FSI_EULERIAN_ELASTIC).or. &
          (FSI_flag(im_solid).eq.FSI_PRESCRIBED_PROBF90).or. &
          (FSI_flag(im_solid).eq.FSI_ICE_EULERIAN_ELASTIC).or. &
          (FSI_flag(im_solid).eq.FSI_ICE_STATIC)) then

  if (FSI_flag(im_solid).eq.FSI_EULERIAN_ELASTIC) then
   if (is_elastic(im_solid).eq.1) then
    !do nothing
   else
    print *,"expecting is_elastic(im_solid)=1"
    stop
   endif
  else if (FSI_flag(im_solid).eq.FSI_ICE_EULERIAN_ELASTIC) then
   if (is_elastic(im_solid).eq.1) then
    !do nothing
   else
    print *,"expecting is_elastic(im_solid)=1"
    stop
   endif
  else if (FSI_flag(im_solid).eq.FSI_ICE_STATIC) then
   if (is_elastic(im_solid).eq.1) then
    !do nothing
   else
    print *,"expecting is_elastic(im_solid)=1"
    stop
   endif
  else if (FSI_flag(im_solid).eq.FSI_PRESCRIBED_PROBF90) then
   if (is_elastic(im_solid).eq.0) then
    !do nothing
   else
    print *,"expecting is_elastic(im_solid)=0"
    stop
   endif
   if (is_rigid(im_solid).eq.1) then
    !do nothing
   else
    print *,"expecting is_rigid(im_solid)=1"
    stop
   endif
  else
   print *,"FSI_flag invalid: ",im_solid,FSI_flag(im_solid)
   stop
  endif
 else
  print *,"FSI_flag invalid: ",im_solid,FSI_flag(im_solid)
  stop
 endif

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

 call flexible_substrateLS(x,LS(im_solid))

 if (FSI_flag(im_solid).eq.FSI_SHOELE_CTML) then
  !do nothing
 else if (FSI_flag(im_solid).eq.FSI_EULERIAN_ELASTIC) then
   if (is_elastic(im_solid).eq.1) then
    !do nothing
   else
    print *,"expecting is_elastic(im_solid)=1"
    stop
   endif
 else if (FSI_flag(im_solid).eq.FSI_PRESCRIBED_PROBF90) then
   if (is_elastic(im_solid).eq.0) then
    !do nothing
   else
    print *,"expecting is_elastic(im_solid)=0"
    stop
   endif
   if (is_rigid(im_solid).eq.1) then
    !do nothing
   else
    print *,"expecting is_rigid(im_solid)=1"
    stop
   endif
 else if (FSI_flag(im_solid).eq.FSI_ICE_EULERIAN_ELASTIC) then
   if (is_elastic(im_solid).eq.1) then
    !do nothing
   else
    print *,"expecting is_elastic(im_solid)=1"
    stop
   endif
 else if (FSI_flag(im_solid).eq.FSI_ICE_STATIC) then
   if (is_elastic(im_solid).eq.1) then
    !do nothing
   else
    print *,"expecting is_elastic(im_solid)=1"
    stop
   endif
 else
  print *,"FSI_flag(im_solid) invalid ",im_solid,FSI_flag
  stop
 endif

 if (num_materials.eq.3) then
  !do nothing
 else
  print *,"num_materials invalid ",num_materials
  stop
 endif

else
 print *,"num_materials or probtype invalid ",num_materials,probtype
 stop
endif

return
end subroutine flexible_plate_impact_LS

subroutine flexible_plate_check_vel_rigid(x,t,vel,dir)
use probcommon_module
IMPLICIT NONE

real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: vel
integer, INTENT(in) :: dir

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
 ! do nothing
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

  real(amrex_real), INTENT(in) :: x(SDIM)
  real(amrex_real), INTENT(in) :: dx(SDIM)
  real(amrex_real), INTENT(in) :: t
  real(amrex_real), INTENT(out) :: LS
  real(amrex_real), INTENT(out) :: vel(SDIM)
  real(amrex_real), INTENT(out) :: temperature
  integer, INTENT(out) :: prescribed_flag
  integer :: dir
  integer, parameter :: im_solid=3
  real(amrex_real) :: radeps,LS_A,LS_B

if (probtype.eq.2000) then

 prescribed_flag=0 !prescribed_flag=1 if "zalesak's" problem
 do dir=1,SDIM
  vel(dir)=zero
 enddo
 LS=-99999.0d0
 LS_A=-99999.0d0
 LS_B=-99999.0d0

 temperature=fort_tempconst(num_materials)

 if (num_materials.eq.im_solid) then

  if (FSI_flag(im_solid).eq.FSI_SHOELE_CTML) then
   LS=-99999.0d0
  else if ((FSI_flag(im_solid).eq.FSI_EULERIAN_ELASTIC).or. &
           (FSI_flag(im_solid).eq.FSI_PRESCRIBED_PROBF90).or. &
           (FSI_flag(im_solid).eq.FSI_ICE_EULERIAN_ELASTIC).or. &
           (FSI_flag(im_solid).eq.FSI_ICE_STATIC)) then

   radeps=radblob2/10.0d0

   if (AMREX_SPACEDIM.eq.2) then
      !squaredist returns Phi<0 in the square.
    call squaredist(x(1),x(2), &
      xblob2-radblob2-radeps, &
      xblob2-radblob2+radeps, &
      yblob2-radblob3, &
      yblob2+radblob3, &
      LS_A)
    LS_A=-LS_A
   
      !squaredist returns Phi<0 in the square.
    call squaredist(x(1),x(2), &
      xblob2+radblob2-radeps, &
      xblob2+radblob2+radeps, &
      yblob2-radblob3, &
      yblob2+radblob3, &
      LS_B)
    LS_B=-LS_B

   else if (AMREX_SPACEDIM.eq.3) then

      !annulusdist returns Phi<0 in the annulus
    call annulusdist( &
      x(1),x(2),x(SDIM), &
      xblob2,yblob2, &
      radblob2, &
      radeps, & 
      zblob2-radblob3, &
      zblob2+radblob3, &
      LS_A)
    LS_A=-LS_A
    LS_B=LS_A

   else
    print *,"dimension bust"
    stop
   endif

   if ((LS_A.ge.zero).or.(LS_B.ge.zero)) then
    LS=99999.0d0
   else if ((LS_A.le.zero).and.(LS_B.le.zero)) then
    LS=-99999.0d0
   else
    print *,"LS_A or LS_B invalid: ",LS_A,LS_B
    stop
   endif

  else
   print *,"FSI_flag invalid: ",im_solid,FSI_flag(im_solid)
   stop
  endif

  if (num_materials.eq.3) then
   !do nothing
  else
   print *,"num_materials invalid ",num_materials
   stop
  endif 

 else
  print *,"num_materials invalid: ",num_materials
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
 print *,"nmat invalid ",nmat
 stop
endif

if ((velsolid_flag.eq.0).or. &
    (velsolid_flag.eq.1)) then
 ! do nothing
else 
 print *,"velsolid_flag invalid ",velsolid_flag
 stop
endif

do dir=1,SDIM
 if (dx(dir).gt.zero) then
  ! do nothing
 else
  print *,"dx invalid ",dx
  stop
 endif
enddo

if (adv_dir.eq.SDIM) then

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

  else if (LS(num_materials).ge.zero) then ! material 3 is the ice or rigid

   do dir=1,SDIM
    VEL(dir)=zero
   enddo

  else
   print *,"LS bust in flexible_plate_impact_VEL LS=",LS
   stop
  endif

else
 print *,"expecting adv_dir = SDIM in flexible_plate_impact_VEL adv_dir=", &
         adv_dir
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
 integer, INTENT(in) :: imattype,im
 real(amrex_real), INTENT(in) :: rho
 real(amrex_real), INTENT(in) :: internal_energy
 real(amrex_real), INTENT(out) :: pressure

 print *,"EOS_flexible_plate_impact should not be used; the default is fine"
 stop

 return
end subroutine EOS_flexible_plate_impact

subroutine SOUNDSQR_flexible_plate_impact(rho,internal_energy,soundsqr, &
  imattype,im)
 IMPLICIT NONE
 integer, INTENT(in) :: imattype,im
 real(amrex_real), INTENT(in) :: rho
 real(amrex_real), INTENT(in) :: internal_energy
 real(amrex_real), INTENT(out) :: soundsqr

 print *,"SOUNDSQR_flexible_plate_impact shouldn't be used; the default is ok"
 stop

 return
end subroutine SOUNDSQR_flexible_plate_impact

subroutine INTERNAL_flexible_plate_impact(rho,temperature, &
  local_internal_energy, &
  imattype,im)
 use global_utility_module
 IMPLICIT NONE
 integer, INTENT(in) :: imattype,im
 real(amrex_real), INTENT(in) :: rho
 real(amrex_real), INTENT(in) :: temperature 
 real(amrex_real), INTENT(out) :: local_internal_energy

 print *,"INTERNAL_flexible_plate_impact shouldn't be used; the default is ok"
 stop

 return
end subroutine INTERNAL_flexible_plate_impact

subroutine TEMPERATURE_flexible_plate_impact(rho,temperature,internal_energy, &
  imattype,im)
 use global_utility_module
 IMPLICIT NONE
 integer, INTENT(in) :: imattype,im
 real(amrex_real), INTENT(in) :: rho
 real(amrex_real), INTENT(out) :: temperature 
 real(amrex_real), INTENT(in) :: internal_energy

 print *,"TEMPERATURE_flexible_plate_impact shouldn't be used; default is ok"
 stop

 return
end subroutine TEMPERATURE_flexible_plate_impact

! This routine will not effect the simulation since
! all of the domain BC will be no-slip.
subroutine flexible_plate_impact_PRES(x,t,LS,PRES,nmat)
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
PRES=A_ELASTIC

return 
end subroutine flexible_plate_impact_PRES


subroutine flexible_plate_impact_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
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
 print *,"num_materials,num_state_material, or probtype invalid ", &
   num_materials,num_state_material,probtype
 stop
endif
 
return
end subroutine flexible_plate_impact_STATE

 ! dir=1..sdim  side=1..2
subroutine flexible_plate_impact_LS_BC(xwall,xghost,t,LS, &
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
 call flexible_plate_impact_LS(xghost,t,LS,nmat)
else
 print *,"dir or side invalid ",dir,side
 stop
endif

return
end subroutine flexible_plate_impact_LS_BC


 ! dir=1..sdim  side=1..2 veldir=1..sdim
subroutine flexible_plate_impact_VEL_BC(xwall,xghost,t,LS, &
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

 call flexible_plate_impact_VEL(xghost,t,LS,local_VEL,velsolid_flag,dx,nmat)
 VEL=local_VEL(veldir)

else
 print *,"dir,side, or veldir invalid ",dir,side,veldir
 stop
endif

return
end subroutine flexible_plate_impact_VEL_BC

 ! dir=1..sdim  side=1..2
subroutine flexible_plate_impact_PRES_BC(xwall,xghost,t,LS, &
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
 print *,"nmat invalid ",nmat
 stop
endif


if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then

 call flexible_plate_impact_PRES(xghost,t,LS,PRES,nmat)

else
 print *,"dir or side invalid ",dir,side
 stop
endif

return
end subroutine flexible_plate_impact_PRES_BC

 ! dir=1..sdim  side=1..2
subroutine flexible_plate_impact_STATE_BC(xwall,xghost,t,LS, &
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
 call flexible_plate_impact_STATE(xghost,t,LS,local_STATE, &
         local_bcflag,nmat,num_state_material)
 ibase=(im-1)*num_state_material
 STATE=local_STATE(ibase+istate)
 call get_primary_material(dx,LS,im_crit)
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

integer, INTENT(in) :: nmat
integer, INTENT(in) :: im
real(amrex_real), INTENT(in) :: VFRAC(nmat)
real(amrex_real), INTENT(in) :: time
integer, INTENT(in) :: nhalf
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
real(amrex_real), INTENT(in) :: temp(nmat)
real(amrex_real), INTENT(in) :: den(nmat)
real(amrex_real), INTENT(in) :: CV(nmat)
real(amrex_real), INTENT(in) :: dt
real(amrex_real), INTENT(out) :: heat_source

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

if ((num_materials.eq.3).and. &
    (probtype.eq.2000)) then
 heat_source=zero
else
 print *,"num_materials or probtype invalid ",num_materials,probtype
 stop
endif

return
end subroutine flexible_plate_impact_HEATSOURCE


subroutine flexible_plate_impact_ASSIMILATE( &
     assimilate_in,assimilate_out, &
     i,j,k,cell_flag)
use probcommon_module
use global_utility_module
use geometry_intersect_module
IMPLICIT NONE

type(assimilate_parm_type), INTENT(in) :: assimilate_in
type(assimilate_out_parm_type), INTENT(inout) :: assimilate_out
integer, INTENT(in) :: i,j,k,cell_flag

integer :: nstate,nstate_test
real(amrex_real) :: xcrit(SDIM)
integer :: dir
integer :: im
real(amrex_real) lnode(4*(SDIM-1))
integer :: i1,j1,k1,k1lo,k1hi,inode
real(amrex_real) :: xdata(SDIM)
real(amrex_real) :: volcell
real(amrex_real) :: cencell(SDIM)
real(amrex_real) :: LS_clamped
real(amrex_real) :: temperature_clamped
real(amrex_real) :: vel_clamped(SDIM)
integer :: prescribed_flag
real(amrex_real) :: LS_flexible(num_materials)
real(amrex_real) :: VFRAC_flexible(num_materials)
real(amrex_real) :: vfrac_override
real(amrex_real) :: vfrac_sum
real(amrex_real) :: centroid_override(SDIM)
real(amrex_real) :: facearea_temp
real(amrex_real) :: lateral_dist
integer :: modify_flag
integer :: vfrac_comp
integer :: tcomp
integer, parameter :: im_solid=3

nstate=assimilate_in%nstate

nstate_test=STATE_NCOMP
if (nstate.eq.nstate_test) then
 ! do nothing
else
 print *,"nstate invalid"
 print *,"nstate=",nstate
 print *,"nstate_test=",nstate_test
 stop
endif

if ((num_materials.eq.3).and. &
    (num_state_material.ge.2).and. & 
    (probtype.eq.2000)) then

 do dir=1,SDIM
  xcrit(dir)=assimilate_in%xsten(0,dir)
 enddo
 
 if (assimilate_in%nhalf.ge.2) then
  ! do nothing
 else
  print *,"(assimilate_in%nhalf.ge.2) violated"
  stop
 endif

 if (SDIM.eq.2) then
  if (abs(yblob2-zblob2).le.EPS11) then
   !do nothing
  else
   print *,"expecting yblob2==zblob2 ",yblob2,zblob2
   stop
  endif
 endif

 if (cell_flag.eq.-1) then
  if (num_materials.eq.3) then

   if (radblob6.gt.zero) then

    if (radblob6.ge.radblob3) then
     print *,"cannot have radblob6>=radblob3 ",radblob3,radblob6
     stop
    endif

    if (radblob6.gt.assimilate_in%dx(1)) then
     if (radblob6.le.four*assimilate_in%dx(1)) then

      if (assimilate_out%LS_state(D_DECL(i,j,k),im_solid).ge.zero) then

       if (xcrit(SDIM).le.zblob2+radblob3-radblob6) then

        if (SDIM.eq.2) then
         lateral_dist=abs(xcrit(1)-xblob2)
        else if (SDIM.eq.3) then
         lateral_dist=sqrt((xcrit(1)-xblob2)**2+(xcrit(2)-yblob2)**2)
        else
         print *,"dimension bust"
         stop
        endif

        if (lateral_dist.le.radblob2) then

         modify_flag=0

         if (xcrit(SDIM).le.zblob2) then
          modify_flag=1
         else if (lateral_dist.ge.half*radblob2) then
          modify_flag=1
         else if &
            (assimilate_out%LS_state(D_DECL(i,j,k),im_solid).ge. &
             (one-EPS2)*radblob6) then
          modify_flag=1
         endif

         if (FSI_flag(im_solid).eq.FSI_ICE_STATIC) then
          if (modify_flag.eq.1) then
           !do nothing
          else
           print *,"expecting modify_flag==1 ",modify_flag
           print *,"FSI_flag= ",FSI_flag
           stop
          endif
         else if (FSI_flag(im_solid).eq.FSI_ICE_EULERIAN_ELASTIC) then
          !do nothing
         else
          print *,"FSI_flag(im_solid) invalid ",im_solid,FSI_flag(im_solid)
          stop
         endif

         if (modify_flag.eq.1) then
          do im=1,num_materials
           tcomp= &
             STATECOMP_STATES+(im-1)*num_state_material+ENUM_TEMPERATUREVAR+1
           assimilate_out%state(D_DECL(i,j,k),tcomp)=fort_tempconst(im_solid)
          enddo
         else if (modify_flag.eq.0) then
          !do nothing
         else
          print *,"modify_flag invalid ",modify_flag
          stop
         endif

        else
         print *,"lateral_dist too big: ",lateral_dist
         print *,"radblob2=",radblob2
         stop
        endif

       else if (xcrit(SDIM).gt.zblob2+radblob3-radblob6) then
        !do nothing
       else
        print *,"parameter bust in flexible_plate_impact_ASSIMILATE"
        print *,"xcrit=",xcrit
        print *,"zblob2= ",zblob2
        print *,"radblob3 ",radblob3
        print *,"radblob6 ",radblob6
        print *,"LS= ",assimilate_out%LS_state(D_DECL(i,j,k),im_solid)
        stop
       endif 

      else if (assimilate_out%LS_state(D_DECL(i,j,k),im_solid).lt.zero) then
       !do nothing
      else
       print *,"parameter bust in flexible_plate_impact_ASSIMILATE"
       print *,"LS= ",assimilate_out%LS_state(D_DECL(i,j,k),im_solid)
       stop
      endif

     else
      print *,"radblob6 too big ",radblob6
      print *,"dx= ",assimilate_in%dx
      stop
     endif
    else
     print *,"radblob6 too small ",radblob6
     print *,"dx= ",assimilate_in%dx
     stop
    endif
   else if (radblob6.eq.zero) then
    !do nothing
   else
    print *,"radblob6 invalid ",radblob6
    stop
   endif

  else
   print *,"num_materials invalid ",num_materials
   stop
  endif
 else if ((cell_flag.ge.0).and.(cell_flag.lt.SDIM)) then
  !do nothing
 else
  print *,"cell_flag invalid ",cell_flag
  stop
 endif

 call flexible_plate_clamped_LS(xcrit,assimilate_in%cur_time, &
   LS_clamped,vel_clamped,temperature_clamped, &
   prescribed_flag,assimilate_in%dx)

 if (LS_clamped.ge.zero) then

  if (cell_flag.eq.0) then ! MAC GRID X
   assimilate_out%macx(D_DECL(i,j,k))=vel_clamped(1)
  else if (cell_flag.eq.1) then ! MAC GRID Y
   assimilate_out%macy(D_DECL(i,j,k))=vel_clamped(2)
  else if ((cell_flag.eq.2).and.(SDIM.eq.3)) then ! MAC GRID Z
   assimilate_out%macz(D_DECL(i,j,k))=vel_clamped(SDIM)
  else if (cell_flag.eq.-1) then
   do dir=1,SDIM
    assimilate_out%state(D_DECL(i,j,k),STATECOMP_VEL+dir)=vel_clamped(dir)
   enddo

   call flexible_plate_impact_LS(xcrit,assimilate_in%cur_time, &
      LS_flexible,num_materials)

   do im=1,num_materials
    assimilate_out%LS_state(D_DECL(i,j,k),im)=LS_flexible(im)
   enddo

   k1lo=0
   k1hi=0
   if (SDIM.eq.3) then
    k1lo=-1
    k1hi=1
   else if (SDIM.eq.2) then
    ! do nothing
   else
    print *,"dimension bust"
    stop
   endif

   do im=1,num_materials

    inode=1
    do k1=k1lo,k1hi,2
    do j1=-1,1,2
    do i1=-1,1,2
     xdata(1)=assimilate_in%xsten(i1,1)
     xdata(2)=assimilate_in%xsten(j1,2)
     if (SDIM.eq.3) then
      xdata(SDIM)=assimilate_in%xsten(k1,SDIM)
     endif
     call flexible_plate_impact_LS(xdata,assimilate_in%cur_time, &
       LS_flexible,num_materials)
     lnode(inode)=LS_flexible(im)
     inode=inode+1
    enddo !i1=-1,1,2
    enddo !j1=-1,1,2
    enddo !k1=-1,1,2
     ! getvolumeNODE is declared in MOF.F90
    call getvolumeNODE( &
     volcell, & !intent(out)
     assimilate_in%bfact, &
     assimilate_in%dx, &
     assimilate_in%xsten, &
     assimilate_in%nhalf, &
     lnode, &
     vfrac_override, &
     facearea_temp, &
     centroid_override, &
     VOFTOL_MATERIAL, &
     SDIM)
    call CISBOX( &
     assimilate_in%xsten, &
     assimilate_in%nhalf, &
     assimilate_in%xlo, &
     assimilate_in%dx, &
     i,j,k, &
     assimilate_in%bfact, &
     assimilate_in%level, &
     volcell,cencell,SDIM)

    vfrac_comp=(im-1)*ngeom_raw+1
    VFRAC_flexible(im)=vfrac_override
    assimilate_out%state(D_DECL(i,j,k),STATECOMP_MOF+vfrac_comp)= &
      vfrac_override
    do dir=1,SDIM
     assimilate_out%state(D_DECL(i,j,k),STATECOMP_MOF+vfrac_comp+dir)= &
      centroid_override(dir)-cencell(dir)
    enddo

   enddo !do im=1,num_materials

   vfrac_sum=zero
   do im=1,num_materials
    if ((is_rigid(im).eq.0).and.(is_elastic(im).eq.0)) then 
     vfrac_sum=vfrac_sum+VFRAC_flexible(im)
    else if ((is_rigid(im).eq.1).or.(is_elastic(im).eq.1)) then 
     !do nothing
    else
     print *,"is_rigid(im) invalid ",im,is_rigid(im)
     print *,"or is_elastic(im) invalid ",im,is_elastic(im)
     stop
    endif
   enddo !im=1,num_materials

   if (vfrac_sum.gt.VOFTOL_MATERIAL) then
    do im=1,num_materials
     vfrac_comp=(im-1)*ngeom_raw+1
     if ((is_rigid(im).eq.0).and.(is_elastic(im).eq.0)) then 
      assimilate_out%state(D_DECL(i,j,k),STATECOMP_MOF+vfrac_comp)= &
       VFRAC_flexible(im)/vfrac_sum
     else if ((is_rigid(im).eq.1).or.(is_elastic(im).eq.1)) then 
      !do nothing
     else
      print *,"is_rigid(im) invalid ",im,is_rigid(im)
      print *,"or is_elastic(im) invalid ",im,is_elastic(im)
      stop
     endif
    enddo !do im=1,num_materials
   else
    print *,"vfrac_sum invalid ",vfrac_sum
    stop
   endif 

  else 
   print *,"cell_flag invalid ",cell_flag
   stop
  endif

 else if (LS_clamped.lt.zero) then
  ! do nothing
 else
  print *,"LS_clamped invalid ",LS_clamped
  stop
 endif

else
 print *,"num_materials,num_state_material, or probtype invalid ", &
    num_materials,num_state_material,probtype
 stop
endif

return
end subroutine flexible_plate_impact_ASSIMILATE


subroutine flexible_plate_impact_OVERRIDE_TAGFLAG( &
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
real(amrex_real) :: F_LIQUID,LS_LIQUID
real(amrex_real) :: F_PLATE,LS_PLATE
integer :: dir
integer :: vofcomp

if (nhalf.lt.3) then
 print *,"nhalf invalid flexible plate tagflag"
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
  print *,"local_delta invalid flexible_plate_impact_override_tagflag"
  stop
 endif
enddo !dir=1..sdim

if ((num_materials.eq.3).and. &
    (probtype.eq.2000)) then

 rflag=0.0d0
 tagflag=0
 LS_LIQUID=lsnew_ptr(D_DECL(i,j,k),1)
 vofcomp=1
 F_LIQUID=snew_ptr(D_DECL(i,j,k),STATECOMP_MOF+vofcomp)
 LS_PLATE=lsnew_ptr(D_DECL(i,j,k),3)
 vofcomp=2*ngeom_raw+1
 F_PLATE=snew_ptr(D_DECL(i,j,k),STATECOMP_MOF+vofcomp)
 if ((LS_LIQUID.ge.zero).or. &
     (F_LIQUID.ge.0.1d0).or. &
     (LS_PLATE.ge.zero).or. &
     (F_PLATE.ge.0.1d0)) then
  rflag=1.0d0
  tagflag=1
 endif

else
 print *,"num_materials or probtype invalid"
 print *,"num_materials: ",num_materials
 print *,"probtype: ",probtype
 stop
endif

end subroutine flexible_plate_impact_OVERRIDE_TAGFLAG


subroutine flexible_plate_impact_VARIABLE_SURFACE_TENSION( &
  xpos, &
  time, &
  iten, &
  temperature, &
  tension)
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: iten
real(amrex_real), INTENT(in) :: time,temperature
real(amrex_real), INTENT(in) :: xpos(SDIM)
real(amrex_real), INTENT(inout) :: tension
real(amrex_real) :: theta

 if ((iten.ge.1).and.(iten.le.num_interfaces)) then
  ! do nothing
 else
  print *,"iten invalid: ",iten
  stop
 endif
 if (temperature.gt.0.0d0) then
  ! do nothing
 else
  print *,"temperature invalid: ",temperature
  stop
 endif
 if (time.ge.0.0d0) then
  ! do nothing
 else
  print *,"time invalid: ",time
  stop
 endif
 if (tension.ge.0.0d0) then
  ! do nothing
 else
  print *,"tension invalid: ",tension
  stop
 endif

 if (probtype.eq.2000) then
   if (radblob9.eq.zero) then
    !do nothing
   else if (radblob9.gt.zero) then
    if (time.eq.zero) then
     tension=fort_tension_init(iten)
    else if (time.ge.radblob9) then
     tension=fort_tension(iten)
    else if ((time.ge.zero).and. &
             (time.le.radblob9)) then
     theta=time/radblob9
     tension=(one-theta)*fort_tension_init(iten)+ &
             theta*fort_tension(iten)
    else
     print *,"time invalid flex_plate_imp_var_surface_tension: ",time
     stop
    endif
   else
    print *,"radblob9 invalid: ",radblob9
    stop
   endif
 else
  print *,"unexpected probtype(flexible_plate_impact SURFTEN): ",probtype
  stop
 endif

end subroutine flexible_plate_impact_VARIABLE_SURFACE_TENSION



end module flexible_plate_impact_module
