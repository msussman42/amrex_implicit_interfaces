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

module STUB_module
use amrex_fort_module, only : amrex_real

implicit none                   

real(amrex_real) :: DEF_VAPOR_GAMMA

contains

! do any initial preparation needed
subroutine INIT_STUB_MODULE()
use probcommon_module
IMPLICIT NONE

 number_of_source_regions=0

 print *,"INIT_STUB_MODULE should not be called"
 stop

return
end subroutine INIT_STUB_MODULE


! do any final deallocation needed
subroutine DEALLOCATE_STUB_MODULE()
use probcommon_module
IMPLICIT NONE

 print *,"DEALLOCATE_STUB_MODULE being called"

return
end subroutine DEALLOCATE_STUB_MODULE


subroutine STUB_CFL_HELPER(time,dir,uu,dx)
IMPLICIT NONE
integer, INTENT(in) :: dir
real(amrex_real), INTENT(in) :: time
real(amrex_real), INTENT(inout) :: uu
real(amrex_real), INTENT(in) :: dx(SDIM)

if ((dir.lt.0).or.(dir.ge.SDIM)) then
 print *,"dir invalid"
 stop
endif

if (dir.eq.0) then
 ! do nothing
else if (dir.eq.1) then
 ! do nothing
else if ((dir.eq.2).and.(SDIM.eq.3)) then
 ! do nothing
else
 print *,"dir invalid STUB_CFL_HELPER"
 stop
endif

return
end subroutine STUB_CFL_HELPER

subroutine STUB_OVERRIDE_TAGFLAG( &
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

 if (nhalf.lt.3) then
  print *,"nhalf invalid stub override tagflag"
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

end subroutine STUB_OVERRIDE_TAGFLAG

subroutine STUB_OVERRIDE_FSI_SIGN_LS_VEL_TEMP( &
 exterior_BB, &
 interior_BB, &
 xcell,time,LS,VEL,TEMP,MASK,lev77,im_part,part_id)
use probcommon_module
use global_utility_module
real(amrex_real), INTENT(in) :: exterior_BB(3,2)
real(amrex_real), INTENT(in) :: interior_BB(3,2)
real(amrex_real), INTENT(in) :: xcell(3)
real(amrex_real), INTENT(in) :: time
real(amrex_real), INTENT(out) :: LS
real(amrex_real), INTENT(out) :: VEL(3)
real(amrex_real), INTENT(out) :: TEMP
integer, INTENT(out) :: MASK
integer, INTENT(in) :: lev77 !lev77=-1 for aux, >=0 otherwise.
integer, INTENT(in) :: im_part
integer, INTENT(in) :: part_id

 if ((lev77.eq.-1).or. &
     (lev77.ge.1)) then
  ! do nothing
 else 
  print *,"lev77 invalid"
  stop
 endif

 MASK=FSI_NOTHING_VALID
 
end subroutine STUB_OVERRIDE_FSI_SIGN_LS_VEL_TEMP

subroutine STUB_GET_OUTSIDE_POINT( &
 exterior_BB, &
 xcell,time,x_outside,im_part,part_id)
use probcommon_module
use global_utility_module
real(amrex_real), INTENT(in) :: exterior_BB(3,2)
real(amrex_real), INTENT(in) :: xcell(3)
real(amrex_real), INTENT(in) :: time
real(amrex_real), INTENT(out) :: x_outside(3)
integer, INTENT(in) :: im_part
integer, INTENT(in) :: part_id
integer :: dir
real(amrex_real) :: BB_len

 do dir=1,3
  BB_len=exterior_BB(dir,2)-exterior_BB(dir,1)
  if (BB_len.gt.zero) then
   x_outside(dir)=exterior_BB(dir,2)+0.05d0*BB_len
  else
   print *,"BB_len invalid"
   stop
  endif
 enddo ! dir=1..3
 
end subroutine STUB_GET_OUTSIDE_POINT

subroutine STUB_AUX_DATA(auxcomp,x,LS)
use probcommon_module
use global_utility_module
IMPLICIT NONE
integer, INTENT(in) :: auxcomp
real(amrex_real), INTENT(in) :: x(3)
real(amrex_real), INTENT(out) :: LS

 print *,"STUB_AUX_DATA should not be called"
 stop

end subroutine STUB_AUX_DATA

subroutine STUB_BOUNDING_BOX_AUX(auxcomp, &
    minnode,maxnode,LS_FROM_SUBROUTINE,aux_ncells_max_side)
use probcommon_module
use global_utility_module
IMPLICIT NONE
integer, INTENT(in) :: auxcomp
real(amrex_real), INTENT(inout) :: minnode(3)
real(amrex_real), INTENT(inout) :: maxnode(3)
integer, INTENT(out) :: LS_FROM_SUBROUTINE
integer, INTENT(out) :: aux_ncells_max_side

 LS_FROM_SUBROUTINE=0
 aux_ncells_max_side=-1
 if ((auxcomp.ge.1).and.(auxcomp.le.fort_num_local_aux_grids)) then
  LS_FROM_SUBROUTINE=0
 else
  print *,"auxcomp invalid in STUB_BOUNDING_BOX_AUX"
  stop
 endif

end subroutine STUB_BOUNDING_BOX_AUX


 ! fluids tessellate the domain, solids are immersed. 
subroutine STUB_LS(x,t,LS,nmat)
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

  print *,"STUB_LS should not be called"
  stop

return
end subroutine STUB_LS

subroutine STUB_check_vel_rigid(x,t,vel,dir)
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

 if (vel.eq.0.0d0) then
  ! do nothing
 else if (vel.gt.0.0d0) then
  ! do nothing
 else if (vel.lt.0.0d0) then
  ! do nothing
 else
  print *,"STUB_check_vel_rigid: vel not expected"
  stop
 endif

return
end subroutine STUB_check_vel_rigid

subroutine STUB_verification_flag(verification_flag)
IMPLICIT NONE

integer, INTENT(out) :: verification_flag

 verification_flag=0

return
end subroutine STUB_verification_flag

subroutine STUB_clamped_LS(x,t,LS,vel,temperature,prescribed_flag,dx)
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
integer dir

 LS=CLAMPED_NO_WHERE_LS
 do dir=1,SDIM
  vel(dir)=zero
 enddo
 temperature=room_temperature
 prescribed_flag=0

return
end subroutine STUB_clamped_LS


subroutine STUB_clamped_LS_jetting_or_cav(x,t,LS,vel, &
               temperature,prescribed_flag,dx)
use probcommon_module
use global_utility_module
use global_distance_module
IMPLICIT NONE

real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: dx(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(out) :: LS
real(amrex_real), INTENT(out) :: vel(SDIM)
real(amrex_real), INTENT(out) :: temperature
integer, INTENT(out) :: prescribed_flag
integer dir
integer, parameter :: for_clamped=1
integer :: solid_id !=1 or 2
integer :: backing_id !=3 or 2

 if ((probtype.eq.42).or. &
     ((probtype.eq.46).and.(axis_dir.eq.10)).or. &
     ((probtype.eq.46).and.(axis_dir.eq.11))) then

  if (probtype.eq.42) then
   backing_id=3
  else if (probtype.eq.46) then
   backing_id=2
  else
   print *,"probtype invalid"
   stop
  endif

  if (num_materials.ge.3) then
   !do nothing
  else
   print *,"expecting num_materials.ge.3: ",num_materials
   stop
  endif

  LS=CLAMPED_NO_WHERE_LS
  do dir=1,SDIM
   vel(dir)=zero
  enddo
  temperature=room_temperature
  prescribed_flag=0 !prescribed_flag=1 if "zalesak's" problem

  if ((FSI_flag(backing_id).eq.FSI_EULERIAN_ELASTIC).or. &
      (FSI_flag(backing_id).eq.FSI_RIGID_NOTPRESCRIBED)) then
   LS=-99999.0d0
   temperature=293.0d0

    !LS<0 in the clamped regions
    !substrate thickness=radblob2
    !          distplate=yblob2
    !substrate length: -xblob2 < r < xblob2
    !substrate vertical: yblob+distplate < z < yblob+distplate+radblob2
    ! (biofilm is on the bottom)
    ! 
   solid_id=1
   call jetting_plate_dist(x(1),x(2),x(SDIM),LS,solid_id,for_clamped)
   LS=-LS
   if (LS.ge.zero) then
    LS=99999.0d0
   else if (LS.lt.zero) then
    LS=-99999.0d0
   else
    print *,"LS is NAN STUB_PROCS.F90 ",LS,solid_id
    stop
   endif
   if (LS.lt.zero) then
    if (num_materials.eq.3) then
     !do nothing
    else if (num_materials.eq.4) then
     solid_id=2
     call jetting_plate_dist(x(1),x(2),x(SDIM),LS,solid_id,for_clamped)
     LS=-LS
     if (LS.ge.zero) then
      LS=99999.0d0
     else if (LS.lt.zero) then
      LS=-99999.0d0
     else
      print *,"LS is NAN STUB_PROCS.F90 ",LS,solid_id
      stop
     endif
    else
     print *,"num_materials invalid: ",num_materials
     stop
    endif
   else if (LS.ge.zero) then
    !do nothing
   else
    print *,"LS invalid: ",LS
    stop
   endif
    
  else if ((FSI_flag(backing_id).eq.FSI_PRESCRIBED_NODES).or. &
           (FSI_flag(backing_id).eq.FSI_SHOELE_CTML).or. &
           (FSI_flag(backing_id).eq.FSI_PRESCRIBED_PROBF90)) then
   !do nothing
  else
   print *,"FSI_flag(backing_id) invalid: ",FSI_flag(backing_id)
   stop
  endif
 else
  print *,"expecting probtype.eq.42 or .eq.46 ",probtype
  print *,"axis_dir ",axis_dir
  stop
 endif

return
end subroutine STUB_clamped_LS_jetting_or_cav

subroutine STUB_clamped_hydrobulge(x,t,LS,vel, &
               temperature,prescribed_flag,dx)
use probcommon_module
use global_utility_module
use global_distance_module
IMPLICIT NONE

real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: dx(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(out) :: LS
real(amrex_real), INTENT(out) :: vel(SDIM)
real(amrex_real), INTENT(out) :: temperature
integer, INTENT(out) :: prescribed_flag
real(amrex_real) :: mag
integer :: dir

 if ((probtype.eq.36).and.(axis_dir.eq.310)) then

  if (num_materials.eq.4) then
   !do nothing
  else
   print *,"expecting num_materials.eq.4: ",num_materials
   stop
  endif

  LS=CLAMPED_NO_WHERE_LS ! -99999.0
  do dir=1,SDIM
   vel(dir)=zero
  enddo
  temperature=room_temperature
  prescribed_flag=0 !prescribed_flag=1 if "zalesak's" problem

  if ((FSI_flag(3).eq.FSI_EULERIAN_ELASTIC).or. &
      (FSI_flag(3).eq.FSI_RIGID_NOTPRESCRIBED)) then
   LS=-99999.0d0
   temperature=293.0d0

   if (SDIM.eq.2) then
    mag=abs(x(1))
   else if (SDIM.eq.3) then
    mag=sqrt(x(1)**2+x(2)**2)
   else
    print *,"dimension bust"
    stop
   endif
   if (abs(x(SDIM)).lt.half*radblob4) then
    !do nothing
   else if (abs(x(SDIM)).ge.half*radblob4) then
    if ((mag.ge.half*radblob3-radblob2).and. &
        (mag.le.half*radblob3+radblob2)) then
     LS=99999.0d0
    else if ((mag.le.half*radblob3-radblob2).or. &
             (mag.ge.half*radblob3+radblob2)) then
     !do nothing
    else
     print *,"mag invalid (STUB_clamped_hydrobulge) ",mag
     stop
    endif
   else
    print *,"abs(x(SDIM)) invalid: ",x(SDIM)
    stop
   endif

  else
   print *,"FSI_flag(3) invalid: ",FSI_flag(3)
   print *,"probtype (STUB_clamped_hydrobulge) = ",probtype
   stop
  endif
 else
  print *,"expecting probtype.eq.36(310) ",probtype
  print *,"axis_dir ",axis_dir
  stop
 endif

return
end subroutine STUB_clamped_hydrobulge



! initial velocity is zero
subroutine STUB_VEL(x,t,LS,VEL,velsolid_flag,dx,nmat)
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

print *,"STUB_VEL should not be called"
stop

return 
end subroutine STUB_VEL

subroutine STUB_INTERNAL_GRAVITY_WAVE_FLAG(internal_wave_exists)
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(out) :: internal_wave_exists

internal_wave_exists=0

end subroutine STUB_INTERNAL_GRAVITY_WAVE_FLAG

! density(T) = density_base * (1+expansion_factor(T))
! remark: expansion_factor(T)=density(T)/density_base - 1
subroutine STUB_UNITLESS_EXPANSION_FACTOR( &
  im,temperature,temperature_base,expansion_factor)
 use probcommon_module
 use global_utility_module
 IMPLICIT NONE

 integer, INTENT(in) :: im
 real(amrex_real), INTENT(in) :: temperature
 real(amrex_real), INTENT(in) :: temperature_base
 real(amrex_real), INTENT(out) :: expansion_factor

 if ((im.ge.1).and.(im.le.num_materials)) then
  if (temperature.gt.zero) then
   if (fort_DrhoDT(im).le.zero) then
    if (temperature_base.gt.zero) then
     expansion_factor=fort_DrhoDT(im)*(temperature-temperature_base)
    else
     print *,"temperature_base must be positive"
     stop
    endif
   else
    print *,"fort_DrhoDT must be nonpositive"
    stop
   endif
  else
   print *,"temperature must be positive"
   stop
  endif

 else
  print *,"im invalid"
  stop
 endif

 return
end subroutine STUB_UNITLESS_EXPANSION_FACTOR

subroutine EOS_STUB(rho,massfrac_var, &
  internal_energy,pressure, &
  imattype,im,num_species_var_in)
 use probcommon_module
 use global_utility_module
 IMPLICIT NONE
 integer, INTENT(in) :: imattype,im,num_species_var_in
 real(amrex_real), INTENT(in) :: rho
 real(amrex_real), INTENT(in) :: massfrac_var(num_species_var_in+1)
 real(amrex_real), INTENT(in) :: internal_energy
 real(amrex_real), INTENT(out) :: pressure

 if (num_species_var_in.eq.num_species_var) then
  call EOS_material_CORE(rho,massfrac_var, &
         internal_energy,pressure,imattype,im)
 else
  print *,"num_species_var_in invalid"
  stop
 endif

 return
end subroutine EOS_STUB


subroutine dVdT_STUB(dVdT,massfrac_var, &
  pressure,temperature, &
  imattype,im,num_species_var_in)
 use probcommon_module
 use global_utility_module
 IMPLICIT NONE
 integer, INTENT(in) :: imattype,im,num_species_var_in
 real(amrex_real), INTENT(in) :: pressure,temperature
 real(amrex_real), INTENT(in) :: massfrac_var(num_species_var_in+1)
 real(amrex_real), INTENT(out) :: dVdT

 if (num_species_var_in.eq.num_species_var) then
  call dVdT_material_CORE(dVdT,massfrac_var, &
         pressure,temperature,imattype,im)
 else
  print *,"num_species_var_in invalid"
  stop
 endif

 return
end subroutine dVdT_STUB



subroutine SOUNDSQR_STUB(rho,massfrac_var, &
  internal_energy,soundsqr, &
  imattype,im,num_species_var_in)
 use probcommon_module
 use global_utility_module
 IMPLICIT NONE
 integer, INTENT(in) :: imattype,im,num_species_var_in
 real(amrex_real), INTENT(in) :: rho
 real(amrex_real), INTENT(in) :: massfrac_var(num_species_var_in+1)
 real(amrex_real), INTENT(in) :: internal_energy
 real(amrex_real), INTENT(out) :: soundsqr

 if (num_species_var_in.eq.num_species_var) then
  call SOUNDSQR_material_CORE(rho,massfrac_var, &
   internal_energy,soundsqr, &
   imattype,im)
 else
  print *,"num_species_var_in invalid"
  stop
 endif

 return
end subroutine SOUNDSQR_STUB


subroutine INTERNAL_STUB(rho,massfrac_var, &
  temperature, &
  local_internal_energy, &
  imattype,im,num_species_var_in)
 use probcommon_module
 use global_utility_module
 IMPLICIT NONE
 integer, INTENT(in) :: imattype,im,num_species_var_in
 real(amrex_real), INTENT(in) :: rho
 real(amrex_real), INTENT(in) :: massfrac_var(num_species_var_in+1)
 real(amrex_real), INTENT(in) :: temperature 
 real(amrex_real), INTENT(out) :: local_internal_energy

 if (num_species_var_in.eq.num_species_var) then
  call INTERNAL_material_CORE(rho,massfrac_var, &
   temperature,local_internal_energy, &
   imattype,im)
 else
  print *,"num_species_var_in invalid"
  stop
 endif

 return
end subroutine INTERNAL_STUB

subroutine TEMPERATURE_STUB(rho,massfrac_var, &
  temperature,internal_energy, &
  imattype,im,num_species_var_in)
 use probcommon_module
 use global_utility_module
 IMPLICIT NONE
 integer, INTENT(in) :: imattype,im,num_species_var_in
 real(amrex_real), INTENT(in) :: rho
 real(amrex_real), INTENT(in) :: massfrac_var(num_species_var_in+1)
 real(amrex_real), INTENT(out) :: temperature 
 real(amrex_real), INTENT(in) :: internal_energy

 if (num_species_var_in.eq.num_species_var) then
  call TEMPERATURE_material_CORE(rho,massfrac_var, &
     temperature,internal_energy, &
     imattype,im)
 else
  print *,"num_species_var_in invalid"
  stop
 endif

 return
end subroutine TEMPERATURE_STUB


! this routine used if pressure boundary conditions are prescribed,
! since only top wall is "outflow" (outflow in quotes since ice shrinks when
! melting), and flow is incompressible, ok to make the top wall pressure zero.
subroutine STUB_PRES(x,t,LS,PRES,nmat)
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
PRES=zero

return 
end subroutine STUB_PRES



subroutine STUB_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
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

print *,"STUB_STATE should not be called"
stop

return
end subroutine STUB_STATE

 ! dir=1..sdim  side=1..2
subroutine STUB_LS_BC(xwall,xghost,t,LS, &
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
 call STUB_LS(xghost,t,LS,nmat)
else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine STUB_LS_BC


 ! dir=1..sdim  side=1..2 veldir=1..sdim
subroutine STUB_VEL_BC(xwall,xghost,t,LS, &
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
if (probtype.eq.-1) then
 ! do nothing
else
 print *,"expecting probtype==-1"
 stop
endif
velsolid_flag=0
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2).and. &
    (veldir.ge.1).and.(veldir.le.SDIM)) then

 call STUB_VEL(xghost,t,LS,local_VEL,velsolid_flag,dx,nmat)
 VEL=local_VEL(veldir)

else
 print *,"dir,side, or veldir invalid"
 stop
endif

return
end subroutine STUB_VEL_BC

 ! dir=1..sdim  side=1..2
subroutine STUB_PRES_BC(xwall,xghost,t,LS, &
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

 call STUB_PRES(xghost,t,LS,PRES,nmat)

else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine STUB_PRES_BC

 ! dir=1..sdim  side=1..2
subroutine STUB_STATE_BC(xwall,xghost,t,LS, &
   STATE,STATE_merge,STATE_in,im,istate,dir,side,dx,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: xwall
real(amrex_real), INTENT(in) :: xghost(SDIM)
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
 call STUB_STATE(xghost,t,LS,local_STATE,local_bcflag, &
         nmat,num_state_material)
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
end subroutine STUB_STATE_BC

! MITSUHIRO: THIS ROUTINE MUST BE CUSTOMIZED.
! water=material 1
! gas=material 2
! ice=material 3
! substrate=material 4
! if (VFRAC(3)>1/2) then
!  set MITSUHIRO_CUSTOM_TEMPERATURE accordingly
!  note: xsten(0,1),xsten(0,2),xsten(0,SDIM) is the coordinate at which a
!  heat source might be prescribed.
!  x(1)=xsten(0,1)
!  x(2)=xsten(0,2)
!  x(SDIM)=xsten(0,SDIM)
!
! this routine called from PROB.F90: subroutine get_local_heat_source
! get_local_heat_source is called from GODUNOV_3D.F90: FORT_HEATSOURCE
subroutine STUB_HEATSOURCE(im,VFRAC,time,x, &
     xsten,nhalf,temp, &
     heat_source,den,CV,dt,nmat)
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

heat_source=zero
 

return
end subroutine STUB_HEATSOURCE


 ! MEHDI VAHAB HEAT SOURCE
 ! called from: GODUNOV_3D.F90, 
 !  HEATSOURCE_FACE when center cell is a solid
 ! and the opposite cell is a fluid cell.
 ! heat_dir=1,2,3
 ! heat_side=1,2
subroutine STUB_EB_heat_source(time,dt,xsten,nhalf, &
      heat_flux,heat_dir,heat_side)
IMPLICIT NONE

integer, INTENT(in) :: nhalf
real(amrex_real), dimension(-nhalf:nhalf,SDIM), INTENT(in) :: xsten
real(amrex_real), INTENT(in) :: time
real(amrex_real), INTENT(in) :: dt
real(amrex_real), INTENT(out) :: heat_flux
integer, INTENT(out) :: heat_dir
integer, INTENT(out) :: heat_side

if (time.lt.zero) then
 print *,"time invalid"
 stop
endif
if (dt.le.zero) then
 print *,"dt invalid"
 stop
endif

heat_flux=zero

return
end subroutine STUB_EB_heat_source

  ! only called at faces with an adjoining solid cell and
  ! an adjoining fluid cell.
subroutine STUB_microcell_heat_coeff(heatcoeff,dx,veldir)
IMPLICIT NONE

real(amrex_real), INTENT(in) :: dx(SDIM)
integer, INTENT(in) :: veldir
real(amrex_real), INTENT(inout) :: heatcoeff


return
end subroutine STUB_microcell_heat_coeff

subroutine STUB_velfreestream(problen,local_buffer)
IMPLICIT NONE

real(amrex_real), INTENT(inout) :: local_buffer(2*SDIM)
real(amrex_real), INTENT(in)    :: problen(SDIM)

return
end subroutine STUB_velfreestream

! this routine called from PROB.F90
subroutine STUB_nucleation(nucleate_in,xsten,nhalf,make_seed)
use probcommon_module_types
IMPLICIT NONE
integer, INTENT(in) :: nhalf
real(amrex_real), dimension(-nhalf:nhalf,SDIM), INTENT(in) :: xsten
integer, INTENT(inout) :: make_seed
type(nucleation_parm_type_input), INTENT(in) :: nucleate_in

! i=nucleate_in%i
! j=nucleate_in%j
! k=nucleate_in%k
! 1<=im<=num_materials
! temperature_component=(im-1)*num_state_material+ENUM_TEMPERATUREVAR+1
! temperature=nucleate_in%EOS(D_DECL(i,j,k),temperature_component)
make_seed=0

return
end subroutine STUB_nucleation

subroutine STUB_ICE_SUBSTRATE_DISTANCE( &
                xtarget,dist)
use probcommon_module
use global_utility_module
IMPLICIT NONE

real(amrex_real), INTENT(in) :: xtarget(SDIM)
real(amrex_real), INTENT(out) :: dist

 dist=-9999.0d0

end subroutine STUB_ICE_SUBSTRATE_DISTANCE


subroutine STUB_correct_pres_rho_hydrostatic( &
  i,j,k,level, &
  angular_velocity, &!INTENT(in) STUB_correct_pres_rho_hydrostatic
  centrifugal_force_factor, &!INTENT(in) STUB_correct_pres_rho_hydrostatic
  dt, &
  rho_hydrostatic, &
  pres_hydrostatic, &
  state_ptr)
IMPLICIT NONE

integer, INTENT(in) :: i,j,k,level
real(amrex_real), INTENT(in) :: angular_velocity
real(amrex_real), INTENT(in) :: centrifugal_force_factor
real(amrex_real), INTENT(in) :: dt
real(amrex_real), INTENT(inout) :: rho_hydrostatic
real(amrex_real), INTENT(inout) :: pres_hydrostatic
real(amrex_real), INTENT(in),pointer :: state_ptr(D_DECL(:,:,:),:)


 if (dt.gt.zero) then
  ! do nothing
 else
  print *,"dt must be positive"
  stop
 endif
 if (angular_velocity.ge.zero) then
  ! do nothing
 else
  print *,"angular_velocity should be nonneg (counter clockwise): ", &
     angular_velocity
  stop
 endif
 if ((centrifugal_force_factor.ge.zero).and. &
     (centrifugal_force_factor.le.one)) then
  ! do nothing
 else
  print *,"expecting 0<=centrifugal_force_factor<=1"
  stop
 endif


end subroutine STUB_correct_pres_rho_hydrostatic

subroutine STUB_FSI_SLICE(xmap3D,xslice3D,problo3D,probhi3D,dx_slice)
use probcommon_module
use global_utility_module
IMPLICIT NONE
real(amrex_real), INTENT(in) :: dx_slice
integer, INTENT(inout) :: xmap3D(3)
real(amrex_real), INTENT(inout) :: xslice3D(3)
real(amrex_real), INTENT(out) :: problo3D(3)
real(amrex_real), INTENT(out) :: probhi3D(3)


  !CTML_FSI_flagF() is declared in GLOBALUTIL.F90
 if (CTML_FSI_flagF().eq.1) then 
  xmap3D(1)=1
  xmap3D(2)=2
  xmap3D(3)=0
  xslice3D(3)=zero
  problo3D(3)=-half*dx_slice
  probhi3D(3)=half*dx_slice
 else if (CTML_FSI_flagF().eq.0) then

   ! 537 is 6 hole injector
  if ((probtype.eq.538).or. &
      (probtype.eq.537).or. &
      (probtype.eq.541)) then
   xmap3D(3)=2
   xmap3D(1)=1
   xmap3D(2)=0
   xslice3D(2)=zero
   problo3D(2)=problo_array(1)
   probhi3D(2)=probhi_array(1)

    ! injector C
   if (probtype.eq.541) then
    if (problo_array(1).ne.zero) then
     print *,"problo_array(1).ne.zero"
     stop
    endif
    xmap3D(1)=1
    xmap3D(2)=2
    xmap3D(3)=0
    xslice3D(1)=zero
    xslice3D(2)=zero
    xslice3D(3)=zero
    problo3D(3)=-1e-3
    probhi3D(3)=1e-3
   endif

  else if (probtype.eq.701) then  ! flapping wing
   xmap3D(1)=1
   xmap3D(3)=2
   xmap3D(2)=0
   xslice3D(2)=0.05
   problo3D(2)=-0.1
   probhi3D(2)=0.2
  else if(probtype.eq.539) then ! the surface is 3D
   xmap3D(1)=1
   xmap3D(2)=2
   xmap3D(3)=0
   xslice3D(3)=0.0
   problo3D(3)=-0.014
   probhi3D(3)=0.014
  else if (probtype.eq.9) then ! ship wave
   xmap3D(1)=1
   xmap3D(3)=2
   xmap3D(2)=0
   xslice3D(2)=0.0
   problo3D(2)=0.0
   probhi3D(2)=0.25
  else if (probtype.eq.5700) then
   xmap3D(1)=1
   xmap3D(2)=2
   xmap3D(3)=0
   xslice3D(3)=0.31
   problo3D(3)=0.0
   probhi3D(3)=0.62
  else if ((probtype.eq.400).or. &
           (probtype.eq.406).or. & ! fractal
           (probtype.eq.404)) then ! gingerbread man or Xue
   xmap3D(1)=1
   xmap3D(2)=2
   xmap3D(3)=0
   xslice3D(3)=zero
   problo3D(3)=-half*dx_slice
   probhi3D(3)=half*dx_slice
  else if (probtype.eq.401) then ! helix
   print *,"this geometry has no 2D analogue"
   stop
  else if (probtype.eq.415) then ! shock sphere
   xmap3D(1)=1
   xmap3D(2)=2
   xmap3D(3)=0
   xslice3D(3)=zero
   problo3D(3)=-half*dx_slice
   probhi3D(3)=half*dx_slice
  endif

 else
  print *,"CTML_FSI_flagF invalid"
  stop
 endif 

end subroutine STUB_FSI_SLICE

subroutine STUB_OPEN_CASFILE(part_id,unit_id,file_format)
IMPLICIT NONE

integer, INTENT(in) :: part_id
integer, INTENT(in) :: unit_id
integer, INTENT(out) :: file_format

 print *,"need to define a routine for SUB_OPEN_CASFILE"
 stop

return
end subroutine STUB_OPEN_CASFILE


subroutine STUB_OPEN_AUXFILE(part_id,unit_id,file_format)
IMPLICIT NONE

integer, INTENT(in) :: part_id
integer, INTENT(in) :: unit_id
integer, INTENT(out) :: file_format

 print *,"need to define a routine for SUB_OPEN_AUXFILE"
 stop

return
end subroutine STUB_OPEN_AUXFILE


subroutine STUB_ORDER_NODES(nodes,nodemap)
IMPLICIT NONE

real(amrex_real), INTENT(in) :: nodes(3,3) ! dir,nodenum
integer, INTENT(inout) :: nodemap(3)

 print *,"need to define a routine for SUB_ORDER_NODES"
 stop

return
end subroutine STUB_ORDER_NODES


subroutine STUB_ASSIMILATE( &
  assimilate_in,assimilate_out,i,j,k,cell_flag)
use probcommon_module
IMPLICIT NONE

type(assimilate_parm_type), INTENT(in) :: assimilate_in
type(assimilate_out_parm_type), INTENT(inout) :: assimilate_out
integer, INTENT(in) :: i,j,k,cell_flag

return
end subroutine STUB_ASSIMILATE


subroutine STUB_SUMINT(GRID_DATA_IN,increment_out1, &
                increment_out2,nsum1,nsum2,isweep)
use probcommon_module_types
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: nsum1,nsum2,isweep
type(user_defined_sum_int_type), INTENT(in) :: GRID_DATA_IN
real(amrex_real), INTENT(inout) :: increment_out1(nsum1)
real(amrex_real), INTENT(inout) :: increment_out2(nsum2)
integer :: i,j,k

i=GRID_DATA_IN%igrid
j=GRID_DATA_IN%jgrid
k=GRID_DATA_IN%kgrid

end subroutine STUB_SUMINT

subroutine STUB_HYDROBULGE_SUMINT(GRID_DATA_IN,increment_out1, &
                increment_out2,nsum1,nsum2,isweep)
use probcommon_module_types
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: nsum1,nsum2,isweep
type(user_defined_sum_int_type), INTENT(in) :: GRID_DATA_IN
real(amrex_real), INTENT(inout) :: increment_out1(nsum1)
real(amrex_real), INTENT(inout) :: increment_out2(nsum2)

real(amrex_real) massfrac_parm(num_species_var+1)
integer im
integer im_primary
integer dir
integer dencomp,local_ispec
real(amrex_real) den,temperature,internal_energy,pressure
real(amrex_real) support_r
real(amrex_real) dx_this_level
real(amrex_real) volgrid
real(amrex_real) denom
real(amrex_real) LS(num_materials)

integer :: level,finest_level

integer :: i,j,k

i=GRID_DATA_IN%igrid
j=GRID_DATA_IN%jgrid
k=GRID_DATA_IN%kgrid
level=GRID_DATA_IN%level
finest_level=GRID_DATA_IN%finest_level

if ((level.le.finest_level).and.(level.ge.0)) then
 ! do nothing
else
 print *,"level invalid"
 stop
endif

if ((num_materials.eq.4).and.(probtype.eq.36).and. &
    (axis_dir.eq.310)) then

 if ((nsum1.eq.1).and.(nsum2.eq.1)) then

  if (isweep.eq.0) then
   increment_out1(1)=zero
  else if (isweep.eq.1) then
   increment_out2(1)=zero
  else
   print *,"isweep invalid"
   stop
  endif
   
  im=1 ! liquid
  dencomp=(im-1)*num_state_material+1+ENUM_DENVAR
  den=GRID_DATA_IN%den(D_DECL(i,j,k),dencomp)
  temperature=GRID_DATA_IN%den(D_DECL(i,j,k),dencomp+1)
  call init_massfrac_parm(den,massfrac_parm,im)
  do local_ispec=1,num_species_var
   massfrac_parm(local_ispec)= &
       GRID_DATA_IN%den(D_DECL(i,j,k),dencomp+1+local_ispec)
  enddo
  call INTERNAL_material_CORE(den,massfrac_parm, &
    temperature,internal_energy,fort_material_type(im), &
    im)
  call EOS_material_CORE(den,massfrac_parm, &
     internal_energy,pressure,fort_material_type(im), &
     im)

  dx_this_level=GRID_DATA_IN%dx(SDIM)

  do im=1,num_materials
   LS(im)=GRID_DATA_IN%lsfab(D_DECL(i,j,k),im)
  enddo
  call get_primary_material(LS,im_primary)

   !liquid,jwl,aluminum side walls,gas
  if (im_primary.eq.1) then
   if ((abs(LS(im_primary)).le.dx_this_level).and. &
       (abs(LS(3)).le.dx_this_level)) then
  
    support_r=0.0d0
    do dir=1,SDIM
     if (dir.ne.1) then
      support_r=support_r+(GRID_DATA_IN%xsten(0,dir))**2
     endif
    enddo
    support_r=sqrt(support_r) 
    if (support_r.le.dx_this_level) then

     volgrid=GRID_DATA_IN%volgrid
     if (isweep.eq.0) then
      increment_out1(1)=volgrid
     else if (isweep.eq.1) then
      denom=increment_out1(1)
      if (denom.gt.0.0d0) then
       increment_out2(1)=volgrid*pressure/denom
      else
       print *,"expecting denom>0.0:",denom
       print *,"nsum1,nsum2 ",nsum1,nsum2
       print *,"charfn,volgrid,pressure,temperature ", &
         volgrid,pressure,temperature
       print *,"i,j,k ",i,j,k
       stop
      endif
     else
      print *,"isweep invalid"
      stop
     endif
    else if (support_r.ge.dx_this_level) then
     !do nothing
    else
     print *,"support_r invalid"
     stop
    endif
   else if ((abs(LS(im_primary)).ge.dx_this_level).or. &
            (abs(LS(3)).ge.dx_this_level)) then
    !do nothing
   else
    print *,"LS(im_primary) invalid"
    stop
   endif
  else if ((im_primary.ge.1).and.(im_primary.le.num_materials)) then
   !do nothing
  else
   print *,"im_primary invalid"
   stop
  endif
 else
  print *,"nsum1 or nsum2 invalid"
  stop
 endif

else
 print *,"num_materials ", num_materials
 print *,"probtype ", probtype
 print *,"axis_dir ", axis_dir
 print *,"num_materials or probtype or axis_dir invalid"
 stop
endif

end subroutine STUB_HYDROBULGE_SUMINT


subroutine STUB_USER_DEFINED_FORCE(xpoint,output_force,force_input)
use probcommon_module_types
use probcommon_module
IMPLICIT NONE

type(user_defined_force_parm_type_input), INTENT(in) :: force_input
real(amrex_real), INTENT(in) :: xpoint(AMREX_SPACEDIM)
real(amrex_real), INTENT(out) :: output_force(AMREX_SPACEDIM)

integer :: dir
integer :: i,j,k
real(amrex_real) :: one_over_den

i=force_input%i
j=force_input%j
k=force_input%k

one_over_den=force_input%one_over_den(D_DECL(i,j,k))
if (one_over_den.gt.zero) then
 ! do nothing
else
 print *,"one_over_den invalid in STUB_USER_DEFINED_FORCE"
 stop
endif

do dir=1,AMREX_SPACEDIM
 output_force(dir)=zero
enddo

end subroutine STUB_USER_DEFINED_FORCE


subroutine STUB_wallfunc( &
  dir, & ! =1,2,3
  data_dir, & ! =0,1,2
  dxmin, &
  x_projection_raster, &
  dx, &
  n_raster, & ! points to solid
  u, & !INTENT(in) uimage_raster_solid_frame(dir)
  uimage_tngt_mag, & !INTENT(in) 
  wall_model_velocity, & ! INTENT(in)
  dist_probe, & ! INTENT(in)
  dist_fluid, & ! INTENT(in)
  temperature_image, & !INTENT(in) 
  temperature_wall, & ! INTENT(in)      
  temperature_wall_max, & ! INTENT(in)      
  viscosity_molecular, & ! INTENT(in)      
  viscosity_eddy_wall, & ! INTENT(in)      
  y, & !INTENT(in) distance from image to wall
  ughost_tngt, & ! INTENT(out)
  im_fluid, &  ! INTENT(in)
  critical_length) ! INTENT(in) used for sanity check
use probcommon_module
use global_utility_module
implicit none
integer, INTENT(in) :: dir ! 1,2,3
integer, INTENT(in) :: data_dir ! 0,1,2
real(amrex_real), INTENT(in) :: dxmin
real(amrex_real), INTENT(in), pointer :: x_projection_raster(:)
real(amrex_real), INTENT(in), pointer :: dx(:)
real(amrex_real), INTENT(in), pointer :: n_raster(:) ! points to solid
integer, INTENT(in) :: im_fluid
real(amrex_real), INTENT(in) :: u !uimage_raster_solid_frame(dir)
real(amrex_real), INTENT(in) :: uimage_tngt_mag
real(amrex_real), INTENT(in) :: wall_model_velocity
real(amrex_real), INTENT(in) :: dist_probe
real(amrex_real), INTENT(in) :: dist_fluid
real(amrex_real), INTENT(in) :: temperature_image
real(amrex_real), INTENT(in) :: temperature_wall
real(amrex_real), INTENT(in) :: temperature_wall_max
real(amrex_real), INTENT(in) :: viscosity_molecular
real(amrex_real), INTENT(in) :: viscosity_eddy_wall
real(amrex_real), INTENT(in) :: y !delta_r
real(amrex_real), INTENT(in) :: critical_length
real(amrex_real), INTENT(out) :: ughost_tngt  ! dir direction

 call wallfunc_newtonsmethod( &
  dir, & ! =1,2,3
  data_dir, & ! =0,1,2
  dxmin, &
  x_projection_raster, &
  dx, &
  n_raster, & ! points to solid
  u, & !INTENT(in) uimage_raster_solid_frame(dir)
  uimage_tngt_mag, & !INTENT(in)
  wall_model_velocity, & ! INTENT(in)
  dist_probe, & ! INTENT(in)
  dist_fluid, & ! INTENT(in)
  temperature_image, & !INTENT(in) 
  temperature_wall, & ! INTENT(in)      
  temperature_wall_max, & ! INTENT(in)      
  viscosity_molecular, & ! INTENT(in)      
  viscosity_eddy_wall, & ! INTENT(in)      
  y, & !INTENT(in) distance from image to wall
  ughost_tngt, & ! INTENT(out)
  im_fluid, &  ! INTENT(in)
  critical_length) ! INTENT(in) used for sanity check

end subroutine STUB_wallfunc

! returns (1/w) where w>>1 in "trouble" regions
subroutine STUB_MAPPING_WEIGHT_COEFF(dir,wt,phys_x)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: dir
real(amrex_real), INTENT(out) :: wt
real(amrex_real), INTENT(in) :: phys_x

if ((dir.ge.0).and.(dir.lt.SDIM)) then
 ! do nothing
else
 print *,"dir invalid"
 stop
endif
if ((phys_x.ge.zero).or.(phys_x.le.zero)) then
 ! do nothing
else
 print *,"phys_x is NaN"
 stop
endif

wt=one

return
end subroutine STUB_MAPPING_WEIGHT_COEFF


subroutine STUB_INIT_REGIONS_LIST(constant_density_all_time, &
      num_materials_in,num_threads_in)
use probcommon_module

IMPLICIT NONE

integer, INTENT(in) :: num_materials_in
integer, INTENT(in) :: num_threads_in
integer, INTENT(in) :: constant_density_all_time(num_materials_in)
integer :: im

 if (num_materials_in.eq.num_materials) then
  ! do nothing
 else
  print *,"num_materials_in invalid"
  stop
 endif
 if (num_threads_in.ge.1) then
  ! do nothing
 else
  print *,"num_threads_in invalid: ",num_threads_in
  stop
 endif
 do im=1,num_materials
  if ((constant_density_all_time(im).eq.0).or. &
      (constant_density_all_time(im).eq.1)) then
   ! do nothing
  else
   print *,"constant_density_all_time(im) invalid"
   stop
  endif
 enddo ! im=1..num_materials

end subroutine STUB_INIT_REGIONS_LIST

subroutine STUB_CHARFN_REGION(region_id,x,cur_time,charfn_out)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: region_id
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: cur_time
real(amrex_real), INTENT(out) :: charfn_out

 ! 1<=region_id<=number_of_source_regions

 print *,"STUB: this routine should not be called if no regions"
 stop

end subroutine STUB_CHARFN_REGION

subroutine STUB_DELETE_REGIONS_LIST()
use probcommon_module
IMPLICIT NONE

end subroutine STUB_DELETE_REGIONS_LIST

subroutine STUB_THERMAL_K(x,dx,cur_time, &
  density, &
  temperature, &
  thermal_k, &
  im, &
  near_interface, &
  im_solid, &
  temperature_wall, &
  temperature_wall_max, &
  temperature_probe, &
  nrm) ! nrm points from solid to fluid
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: im
integer, INTENT(in) :: im_solid
integer, INTENT(in) :: near_interface
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: dx(SDIM)
real(amrex_real), INTENT(in) :: cur_time
real(amrex_real), INTENT(in) :: density
real(amrex_real), INTENT(in) :: temperature
real(amrex_real), INTENT(in) :: temperature_wall
real(amrex_real), INTENT(in) :: temperature_wall_max
real(amrex_real), INTENT(in) :: temperature_probe
real(amrex_real), INTENT(in) :: nrm(SDIM) ! nrm points from solid to fluid
real(amrex_real), INTENT(inout) :: thermal_k

if ((im.ge.1).and.(im.le.num_materials)) then
 ! do nothing
else 
 print *,"im invalid in STUB_THERMAL_K"
 stop
endif
if ((im_solid.ge.0).and.(im_solid.le.num_materials)) then
 ! do nothing
else
 print *,"im_solid invalid STUB_THERMAL_K"
 stop
endif

end subroutine STUB_THERMAL_K

subroutine STUB_INTERFACE_TEMPERATURE( &
  interface_mass_transfer_model, &
  probe_constrain, &
  ireverse, &
  iten, &        
  xI, &        
  cur_time, &        
  prev_time, &        
  dt, &        
  TI, &
  YI, &
  user_override_TI_YI, &
  molar_mass, & ! index: 1..nmat
  species_molar_mass, & ! index: 1..num_species_var
  ksrc_predict, &
  kdst_predict, &
  ksrc_physical, &
  kdst_physical, &
  T_probe_src, &
  T_probe_dst, &
  LL, &
  dxprobe_src, &
  dxprobe_dst, &
  num_materials_in, &
  num_species_var_in)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: interface_mass_transfer_model
integer, INTENT(in) :: num_materials_in
integer, INTENT(in) :: num_species_var_in
integer, INTENT(in) :: probe_constrain
integer, INTENT(in) :: ireverse
integer, INTENT(in) :: iten
real(amrex_real), INTENT(in) :: xI(SDIM)
real(amrex_real), INTENT(in) :: cur_time
real(amrex_real), INTENT(in) :: prev_time
real(amrex_real), INTENT(in) :: dt
real(amrex_real), INTENT(inout) :: TI
real(amrex_real), INTENT(inout) :: YI
integer, INTENT(inout) :: user_override_TI_YI
real(amrex_real), INTENT(in) :: molar_mass(num_materials_in)
real(amrex_real), INTENT(in) :: species_molar_mass(num_species_var_in)
real(amrex_real), INTENT(in) :: ksrc_predict
real(amrex_real), INTENT(in) :: kdst_predict
real(amrex_real), INTENT(in) :: ksrc_physical
real(amrex_real), INTENT(in) :: kdst_physical
real(amrex_real), INTENT(in) :: T_probe_src
real(amrex_real), INTENT(in) :: T_probe_dst
real(amrex_real), INTENT(in) :: LL
real(amrex_real), INTENT(in) :: dxprobe_src
real(amrex_real), INTENT(in) :: dxprobe_dst

end subroutine STUB_INTERFACE_TEMPERATURE

subroutine STUB_MDOT( &
  num_materials_in, &
  num_species_var_in, &
  interface_mass_transfer_model, &
  xI, & 
  ispec, &
  molar_mass, & ! 1..nmat
  species_molar_mass, & ! 1..num_species_var+1
  im_source, &
  im_dest, &
  mdot, & ! INTENT(out)
  mdot_override, & ! INTENT(inout)
  ksrc_derived, &
  kdst_derived, &
  ksrc_physical, &
  kdst_physical, &
  interp_valid_flag_src, &
  interp_valid_flag_dst, &
  T_probe_src, &
  T_probe_dst, &
  TI, &
  LL, &
  dxprobe_src, &
  dxprobe_dst)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: interface_mass_transfer_model
integer, INTENT(in) :: num_materials_in
integer, INTENT(in) :: num_species_var_in
integer, INTENT(in) :: ispec
integer, INTENT(in) :: im_source
integer, INTENT(in) :: im_dest
real(amrex_real), INTENT(in) :: xI(SDIM)
real(amrex_real), INTENT(in) :: TI
real(amrex_real), INTENT(in) :: molar_mass(num_materials_in)
real(amrex_real), INTENT(in) :: species_molar_mass(num_species_var_in+1)
real(amrex_real), INTENT(out) :: mdot
integer, INTENT(inout) :: mdot_override
real(amrex_real), INTENT(in) :: ksrc_derived
real(amrex_real), INTENT(in) :: kdst_derived
real(amrex_real), INTENT(in) :: ksrc_physical
real(amrex_real), INTENT(in) :: kdst_physical
integer, INTENT(in) :: interp_valid_flag_src
integer, INTENT(in) :: interp_valid_flag_dst
real(amrex_real), INTENT(in) :: T_probe_src
real(amrex_real), INTENT(in) :: T_probe_dst
real(amrex_real), INTENT(in) :: LL
real(amrex_real), INTENT(in) :: dxprobe_src
real(amrex_real), INTENT(in) :: dxprobe_dst
real(amrex_real) DTsrc,DTdst,mdotsrc,mdotdst,mdotsum

mdot_override=0

if (interface_mass_transfer_model.eq.0) then
 ! do nothing
else if (interface_mass_transfer_model.eq.999) then
 mdot_override=1
 if (interp_valid_flag_src.eq.1) then
  DTsrc=T_probe_src-TI
 else if (interp_valid_flag_src.eq.0) then
  DTsrc=zero
 else
  print *,"interp_valid_flag_src invalid"
  stop
 endif
 if (interp_valid_flag_dst.eq.1) then
  DTdst=T_probe_dst-TI
 else if (interp_valid_flag_dst.eq.0) then
  DTdst=zero
 else
  print *,"interp_valid_flag_dst invalid"
  stop
 endif

 mdotsrc=ksrc_derived*DTsrc/(LL*dxprobe_src)
 mdotdst=kdst_derived*DTdst/(LL*dxprobe_dst)
 mdotsum=mdotsrc+mdotdst
 if (mdotsum.gt.zero) then
  ! do nothing
 else if (mdotsum.le.zero) then
  mdotsum=zero
 else
  print *,"mdotsum invalid in STUB_MDOT"
  stop
 endif
 mdot=mdotsum
else if (interface_mass_transfer_model.gt.0) then
 ! do nothing
else
 print *,"interface_mass_transfer_model invalid"
 stop
endif

end subroutine STUB_MDOT

subroutine STUB_K_EFFECTIVE( &
  interface_mass_transfer_model, &
  ireverse, &
  iten, &        
  molar_mass, & ! index: 1..nmat
  species_molar_mass, & ! index: 1..num_species_var
  k_model_predict, &
  k_model_correct, &
  k_physical_base, &
  T_probe_src, &
  T_probe_dst, &
  dxprobe_src, &
  dxprobe_dst, &
  LL, &
  num_materials_in, &
  num_species_var_in)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: interface_mass_transfer_model
integer, INTENT(in) :: num_materials_in
integer, INTENT(in) :: num_species_var_in
integer, INTENT(in) :: ireverse
integer, INTENT(in) :: iten
real(amrex_real), INTENT(in) :: molar_mass(num_materials_in)
real(amrex_real), INTENT(in) :: species_molar_mass(num_species_var_in)
real(amrex_real), INTENT(in) :: k_model_predict(2) ! src,dst
real(amrex_real), INTENT(inout) :: k_model_correct(2) ! src,dst
real(amrex_real), INTENT(in) :: k_physical_base(2) ! src, dst
real(amrex_real), INTENT(in) :: T_probe_src
real(amrex_real), INTENT(in) :: T_probe_dst
real(amrex_real), INTENT(in) :: LL
real(amrex_real), INTENT(in) :: dxprobe_src
real(amrex_real), INTENT(in) :: dxprobe_dst

end subroutine STUB_K_EFFECTIVE

subroutine STUB_reference_wavelen(wavelen)
use probcommon_module
use global_utility_module
IMPLICIT NONE
real(amrex_real), INTENT(inout) :: wavelen
real(amrex_real) :: default_wavelen
integer :: dir_local
integer :: gravity_dir

 call fort_derive_gravity_dir(gravity_vector,gravity_dir)
 if ((gravity_dir.ge.1).and.(gravity_dir.le.SDIM)) then
  ! do nothing
 else
  print *,"gravity_dir invalid"
  stop
 endif

 default_wavelen=zero
 do dir_local=1,SDIM
  if (gravity_dir.ne.dir_local) then
   default_wavelen=default_wavelen+problen_array(dir_local)**2
  endif
 enddo !dir_local=1..sdim
 default_wavelen=sqrt(default_wavelen)

 if ((wavelen.gt.zero).and. &
     (wavelen.le.default_wavelen*(one+EPS2))) then
  ! do nothing
 else
  print *,"input wavelen out of range"
  print *,"wavelen (input) = ",wavelen
  print *,"default_wavelen (diagonal distance) = ",default_wavelen
  stop
 endif

end subroutine STUB_reference_wavelen

subroutine STUB_VARIABLE_SURFACE_TENSION( &
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

 if ((iten.ge.1).and.(iten.le.num_interfaces)) then
  ! do nothing
 else
  print *,"iten invalid"
  stop
 endif
 if (temperature.gt.0.0d0) then
  ! do nothing
 else
  print *,"temperature invalid"
  stop
 endif
 if (tension.ge.0.0d0) then
  ! do nothing
 else
  print *,"tension invalid"
  stop
 endif

end subroutine STUB_VARIABLE_SURFACE_TENSION

subroutine STUB_VARIABLE_LATENT_HEAT( &
  iten, &
  temperature, &
  latent_heat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: iten
real(amrex_real), INTENT(in) :: temperature
real(amrex_real), INTENT(inout) :: latent_heat ! always positive

 if ((iten.ge.1).and.(iten.le.2*num_interfaces)) then
  ! do nothing
 else
  print *,"iten invalid"
  stop
 endif
 if (temperature.gt.0.0d0) then
  ! do nothing
 else
  print *,"temperature invalid"
  stop
 endif
 if (latent_heat.gt.0.0d0) then
  ! do nothing
 else
  print *,"latent_heat invalid"
  stop
 endif

end subroutine STUB_VARIABLE_LATENT_HEAT

subroutine STUB_T0_Boussinesq(x,dx,cur_time,im,T0)
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

 T0=fort_tempconst(im)

end subroutine STUB_T0_Boussinesq

subroutine STUB_V0_Coriolis(x,dx,cur_time,V0)
use probcommon_module
use global_utility_module
IMPLICIT NONE

real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: dx(SDIM)
real(amrex_real), INTENT(in) :: cur_time
real(amrex_real), INTENT(out) :: V0(SDIM)

integer dir

 if (cur_time.ge.0.0d0) then
  ! do nothing
 else
  print *,"cur_time invalid: ",cur_time
  stop
 endif

 do dir=1,SDIM
  V0(dir)=0.0d0
 enddo

end subroutine STUB_V0_Coriolis

subroutine STUB_angular_velocity(x,cur_time, &
   angular_velocity,angular_velocity_custom, &
   angular_velocity_dot,lever_arm)
use probcommon_module
use global_utility_module
IMPLICIT NONE

real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: cur_time
real(amrex_real), INTENT(in) :: angular_velocity
real(amrex_real), INTENT(out) :: angular_velocity_custom
real(amrex_real), INTENT(out) :: angular_velocity_dot
real(amrex_real), INTENT(out) :: lever_arm

 if (cur_time.ge.0.0d0) then
  ! do nothing
 else
  print *,"cur_time invalid: ",cur_time
  stop
 endif

 angular_velocity_custom=angular_velocity
 angular_velocity_dot=zero
 lever_arm=zero

end subroutine STUB_angular_velocity


subroutine STUB_gravity_vector(x,cur_time, &
   gravity_vector_in, &
   gravity_vector_out)
use probcommon_module
use global_utility_module
IMPLICIT NONE

real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: cur_time
real(amrex_real), INTENT(in) :: gravity_vector_in(SDIM)
real(amrex_real), INTENT(out) :: gravity_vector_out(SDIM)
integer :: dir

 if (cur_time.ge.0.0d0) then
  ! do nothing
 else
  print *,"cur_time invalid STUB_gravity_vector: ",cur_time
  stop
 endif

 do dir=1,SDIM
  gravity_vector_out(dir)=gravity_vector_in(dir)
 enddo

end subroutine STUB_gravity_vector

end module STUB_module
