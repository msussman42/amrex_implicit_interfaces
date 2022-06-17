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

implicit none                   

REAL_T :: DEF_VAPOR_GAMMA

contains

  ! do any initial preparation needed
subroutine INIT_STUB_MODULE()
IMPLICIT NONE

print *,"INIT_STUB_MODULE should not be called"
stop

return
end subroutine INIT_STUB_MODULE


subroutine STUB_CFL_HELPER(time,dir,uu,dx)
IMPLICIT NONE
INTEGER_T, intent(in) :: dir
REAL_T, intent(in) :: time
REAL_T, intent(inout) :: uu
REAL_T, intent(in) :: dx(SDIM)

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

subroutine STUB_OVERRIDE_TAGFLAG(xsten,nhalf,time,rflag,tagflag)
use probcommon_module
use global_utility_module
IMPLICIT NONE
INTEGER_T, intent(in) :: nhalf
REAL_T, intent(in) :: xsten(-nhalf:nhalf,SDIM)
REAL_T, intent(in) :: time
REAL_T, intent(inout) :: rflag
INTEGER_T, intent(inout) :: tagflag

 if (nhalf.lt.1) then
  print *,"nhalf invalid stub override tagflag"
  stop
 endif

end subroutine STUB_OVERRIDE_TAGFLAG

subroutine STUB_OVERRIDE_FSI_SIGN_LS_VEL_TEMP( &
 xcell,time,LS,VEL,TEMP,MASK,lev77,im_part,part_id)
use probcommon_module
use global_utility_module
REAL_T, intent(in) :: xcell(3)
REAL_T, intent(in) :: time
REAL_T, intent(out) :: LS
REAL_T, intent(out) :: VEL(3)
REAL_T, intent(out) :: TEMP
INTEGER_T, intent(out) :: MASK
INTEGER_T, intent(in) :: lev77 !lev77=-1 for aux, >=0 otherwise.
INTEGER_T, intent(in) :: im_part
INTEGER_T, intent(in) :: part_id

 if ((lev77.eq.-1).or. &
     (lev77.ge.1)) then
  ! do nothing
 else 
  print *,"lev77 invalid"
  stop
 endif

 MASK=FSI_NOTHING_VALID
 
end subroutine STUB_OVERRIDE_FSI_SIGN_LS_VEL_TEMP

subroutine STUB_AUX_DATA(auxcomp,x,LS)
use probcommon_module
use global_utility_module
IMPLICIT NONE
INTEGER_T, intent(in) :: auxcomp
REAL_T, intent(in) :: x(3)
REAL_T, intent(out) :: LS

 print *,"STUB_AUX_DATA should not be called"
 stop

end subroutine STUB_AUX_DATA

subroutine STUB_BOUNDING_BOX_AUX(auxcomp, &
    minnode,maxnode,LS_FROM_SUBROUTINE,aux_ncells_max_side)
use probcommon_module
use global_utility_module
IMPLICIT NONE
INTEGER_T, intent(in) :: auxcomp
REAL_T, intent(inout) :: minnode(3)
REAL_T, intent(inout) :: maxnode(3)
INTEGER_T, intent(out) :: LS_FROM_SUBROUTINE
INTEGER_T, intent(out) :: aux_ncells_max_side

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

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(out) :: LS(nmat)

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

REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: vel
INTEGER_T, intent(in) :: dir

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

subroutine STUB_clamped_LS(x,t,LS,vel,temperature,dx)
use probcommon_module
use global_utility_module
IMPLICIT NONE

REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: dx(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(out) :: LS
REAL_T, intent(out) :: vel(SDIM)
REAL_T, intent(out) :: temperature
INTEGER_T dir

 LS=CLAMPED_NO_WHERE_LS
 do dir=1,SDIM
  vel(dir)=zero
 enddo
 temperature=293.0d0

return
end subroutine STUB_clamped_LS



! initial velocity is zero
subroutine STUB_VEL(x,t,LS,VEL,velsolid_flag,dx,nmat)
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

print *,"STUB_VEL should not be called"
stop

return 
end subroutine STUB_VEL

! density(T) = density_base * (1+expansion_factor(T))
! remark: expansion_factor(T)=density(T)/density_base - 1
subroutine STUB_UNITLESS_EXPANSION_FACTOR( &
  im,temperature,temperature_base,expansion_factor)
 use probcommon_module
 use global_utility_module
 IMPLICIT NONE

 INTEGER_T, intent(in) :: im
 REAL_T, intent(in) :: temperature
 REAL_T, intent(in) :: temperature_base
 REAL_T, intent(out) :: expansion_factor

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
 INTEGER_T, intent(in) :: imattype,im,num_species_var_in
 REAL_T, intent(in) :: rho
 REAL_T, intent(in) :: massfrac_var(num_species_var_in+1)
 REAL_T, intent(in) :: internal_energy
 REAL_T, intent(out) :: pressure

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
 INTEGER_T, intent(in) :: imattype,im,num_species_var_in
 REAL_T, intent(in) :: pressure,temperature
 REAL_T, intent(in) :: massfrac_var(num_species_var_in+1)
 REAL_T, intent(out) :: dVdT

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
 INTEGER_T, intent(in) :: imattype,im,num_species_var_in
 REAL_T, intent(in) :: rho
 REAL_T, intent(in) :: massfrac_var(num_species_var_in+1)
 REAL_T, intent(in) :: internal_energy
 REAL_T, intent(out) :: soundsqr

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
 INTEGER_T, intent(in) :: imattype,im,num_species_var_in
 REAL_T, intent(in) :: rho
 REAL_T, intent(in) :: massfrac_var(num_species_var_in+1)
 REAL_T, intent(in) :: temperature 
 REAL_T, intent(out) :: local_internal_energy

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
 INTEGER_T, intent(in) :: imattype,im,num_species_var_in
 REAL_T, intent(in) :: rho
 REAL_T, intent(in) :: massfrac_var(num_species_var_in+1)
 REAL_T, intent(out) :: temperature 
 REAL_T, intent(in) :: internal_energy

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
end subroutine STUB_PRES



subroutine STUB_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
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
if (probtype.eq.55) then
 ! do nothing
else
 print *,"expecting probtype==55"
 stop
endif
velsolid_flag=0
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2).and. &
    (veldir.ge.1).and.(veldir.le.SDIM)) then

 call STUB_VEL(xghost,t,LS,local_VEL,velsolid_flag,dx,nmat)

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
 call STUB_STATE(xghost,t,LS,local_STATE,local_bcflag, &
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

INTEGER_T, intent(in) :: nhalf
REAL_T, dimension(-nhalf:nhalf,SDIM), intent(in) :: xsten
REAL_T, intent(in) :: time
REAL_T, intent(in) :: dt
REAL_T, intent(out) :: heat_flux
INTEGER_T, intent(out) :: heat_dir
INTEGER_T, intent(out) :: heat_side

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

REAL_T, intent(in) :: dx(SDIM)
INTEGER_T, intent(in) :: veldir
REAL_T, intent(inout) :: heatcoeff


return
end subroutine STUB_microcell_heat_coeff

subroutine STUB_velfreestream(problen,local_buffer)
IMPLICIT NONE

REAL_T, intent(inout) :: local_buffer(2*SDIM)
REAL_T, intent(in)    :: problen(SDIM)

return
end subroutine STUB_velfreestream


subroutine STUB_nucleation(nucleate_in,xsten,nhalf,make_seed)
use probcommon_module_types
IMPLICIT NONE
INTEGER_T, intent(in) :: nhalf
REAL_T, dimension(-nhalf:nhalf,SDIM), intent(in) :: xsten
INTEGER_T, intent(inout) :: make_seed
type(nucleation_parm_type_input), intent(in) :: nucleate_in

make_seed=0

return
end subroutine STUB_nucleation

subroutine STUB_ICE_SUBSTRATE_DISTANCE( &
                xtarget,dist)
use probcommon_module
use global_utility_module
IMPLICIT NONE

REAL_T, intent(in) :: xtarget(SDIM)
REAL_T, intent(out) :: dist

 dist=-9999.0d0

end subroutine STUB_ICE_SUBSTRATE_DISTANCE


subroutine STUB_correct_pres_rho_hydrostatic( &
  i,j,k,level, &
  gravity_normalized, &
  gravity_dir_parm, &
  angular_velocity, &
  dt, &
  rho_hydrostatic, &
  pres_hydrostatic, &
  state_ptr)
IMPLICIT NONE

INTEGER_T, intent(in) :: i,j,k,level
INTEGER_T, intent(in) :: gravity_dir_parm
REAL_T, intent(in) :: angular_velocity
REAL_T, intent(in) :: gravity_normalized
REAL_T, intent(in) :: dt
REAL_T, intent(inout) :: rho_hydrostatic
REAL_T, intent(inout) :: pres_hydrostatic
REAL_T, intent(in),pointer :: state_ptr(D_DECL(:,:,:),:)

end subroutine STUB_correct_pres_rho_hydrostatic

subroutine STUB_FSI_SLICE(xmap3D,xslice3D,problo3D,probhi3D,dx_slice)
use probcommon_module
use global_utility_module
IMPLICIT NONE
REAL_T, intent(in) :: dx_slice
INTEGER_T, intent(inout) :: xmap3D(3)
REAL_T, intent(inout) :: xslice3D(3)
REAL_T, intent(out) :: problo3D(3)
REAL_T, intent(out) :: probhi3D(3)


  !CTML_FSI_flagF(num_materials) is declared in GLOBALUTIL.F90
 if (CTML_FSI_flagF(num_materials).eq.1) then ! FSI_flag==4 or 8
  xmap3D(1)=1
  xmap3D(2)=2
  xmap3D(3)=0
  xslice3D(3)=zero
  problo3D(3)=-half*dx_slice
  probhi3D(3)=half*dx_slice
 else if (CTML_FSI_flagF(num_materials).eq.0) then

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

INTEGER_T, intent(in) :: part_id
INTEGER_T, intent(in) :: unit_id
INTEGER_T, intent(out) :: file_format

 print *,"need to define a routine for SUB_OPEN_CASFILE"
 stop

return
end subroutine STUB_OPEN_CASFILE


subroutine STUB_OPEN_AUXFILE(part_id,unit_id,file_format)
IMPLICIT NONE

INTEGER_T, intent(in) :: part_id
INTEGER_T, intent(in) :: unit_id
INTEGER_T, intent(out) :: file_format

 print *,"need to define a routine for SUB_OPEN_AUXFILE"
 stop

return
end subroutine STUB_OPEN_AUXFILE


subroutine STUB_ORDER_NODES(nodes,nodemap)
IMPLICIT NONE

REAL_T, intent(in) :: nodes(3,3) ! dir,nodenum
INTEGER_T, intent(inout) :: nodemap(3)

 print *,"need to define a routine for SUB_ORDER_NODES"
 stop

return
end subroutine STUB_ORDER_NODES


subroutine STUB_ASSIMILATE( &
  assimilate_in,assimilate_out,i,j,k,cell_flag)
use probcommon_module
IMPLICIT NONE

type(assimilate_parm_type), intent(in) :: assimilate_in
type(assimilate_out_parm_type), intent(inout) :: assimilate_out
INTEGER_T, intent(in) :: i,j,k,cell_flag

return
end subroutine STUB_ASSIMILATE


subroutine STUB_SUMINT(GRID_DATA_IN,increment_out1, &
                increment_out2,nsum1,nsum2,isweep)
use probcommon_module_types
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nsum1,nsum2,isweep
type(user_defined_sum_int_type), intent(in) :: GRID_DATA_IN
REAL_T, intent(inout) :: increment_out1(nsum1)
REAL_T, intent(inout) :: increment_out2(nsum2)
INTEGER_T :: i,j,k

i=GRID_DATA_IN%igrid
j=GRID_DATA_IN%jgrid
k=GRID_DATA_IN%kgrid

end subroutine STUB_SUMINT

subroutine STUB_wallfunc( &
  dir, & ! =1,2,3
  data_dir, & ! =0,1,2
  dxmin, &
  x_projection_raster, &
  dx, &
  n_raster, & ! points to solid
  u, & !intent(in) uimage_raster_solid_frame(dir)
  uimage_tngt_mag, & !intent(in) 
  wall_model_velocity, & ! intent(in)
  dist_probe, & ! intent(in)
  dist_fluid, & ! intent(in)
  temperature_image, & !intent(in) 
  temperature_wall, & ! intent(in)      
  temperature_wall_max, & ! intent(in)      
  viscosity_molecular, & ! intent(in)      
  viscosity_eddy_wall, & ! intent(in)      
  y, & !intent(in) distance from image to wall
  ughost_tngt, & ! intent(out)
  im_fluid, &  ! intent(in)
  critical_length) ! intent(in) used for sanity check
use probcommon_module
use global_utility_module
implicit none
INTEGER_T, intent(in) :: dir ! 1,2,3
INTEGER_T, intent(in) :: data_dir ! 0,1,2
REAL_T, intent(in) :: dxmin
REAL_T, intent(in), pointer :: x_projection_raster(:)
REAL_T, intent(in), pointer :: dx(:)
REAL_T, intent(in), pointer :: n_raster(:) ! points to solid
INTEGER_T, intent(in) :: im_fluid
REAL_T, intent(in) :: u !uimage_raster_solid_frame(dir)
REAL_T, intent(in) :: uimage_tngt_mag
REAL_T, intent(in) :: wall_model_velocity
REAL_T, intent(in) :: dist_probe
REAL_T, intent(in) :: dist_fluid
REAL_T, intent(in) :: temperature_image
REAL_T, intent(in) :: temperature_wall
REAL_T, intent(in) :: temperature_wall_max
REAL_T, intent(in) :: viscosity_molecular
REAL_T, intent(in) :: viscosity_eddy_wall
REAL_T, intent(in) :: y !delta_r
REAL_T, intent(in) :: critical_length
REAL_T, intent(out) :: ughost_tngt  ! dir direction

 call wallfunc_newtonsmethod( &
  dir, & ! =1,2,3
  data_dir, & ! =0,1,2
  dxmin, &
  x_projection_raster, &
  dx, &
  n_raster, & ! points to solid
  u, & !intent(in) uimage_raster_solid_frame(dir)
  uimage_tngt_mag, & !intent(in)
  wall_model_velocity, & ! intent(in)
  dist_probe, & ! intent(in)
  dist_fluid, & ! intent(in)
  temperature_image, & !intent(in) 
  temperature_wall, & ! intent(in)      
  temperature_wall_max, & ! intent(in)      
  viscosity_molecular, & ! intent(in)      
  viscosity_eddy_wall, & ! intent(in)      
  y, & !intent(in) distance from image to wall
  ughost_tngt, & ! intent(out)
  im_fluid, &  ! intent(in)
  critical_length) ! intent(in) used for sanity check

end subroutine STUB_wallfunc

subroutine STUB_INIT_REGIONS_LIST(constant_density_all_time, &
      num_materials_in,num_threads_in)
use probcommon_module

IMPLICIT NONE

INTEGER_T, intent(in) :: num_materials_in
INTEGER_T, intent(in) :: num_threads_in
INTEGER_T, intent(in) :: constant_density_all_time(num_materials_in)
INTEGER_T :: im

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

 number_of_source_regions=0

end subroutine STUB_INIT_REGIONS_LIST

subroutine STUB_CHARFN_REGION(region_id,x,cur_time,charfn_out)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: region_id
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: cur_time
REAL_T, intent(out) :: charfn_out

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

INTEGER_T, intent(in) :: im
INTEGER_T, intent(in) :: im_solid
INTEGER_T, intent(in) :: near_interface
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: dx(SDIM)
REAL_T, intent(in) :: cur_time
REAL_T, intent(in) :: density
REAL_T, intent(in) :: temperature
REAL_T, intent(in) :: temperature_wall
REAL_T, intent(in) :: temperature_wall_max
REAL_T, intent(in) :: temperature_probe
REAL_T, intent(in) :: nrm(SDIM) ! nrm points from solid to fluid
REAL_T, intent(inout) :: thermal_k

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
  probe_ok_gradient_src, &
  probe_ok_gradient_dst, &
  LL, &
  dxprobe_src, &
  dxprobe_dst, &
  num_materials_in, &
  num_species_var_in)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: interface_mass_transfer_model
INTEGER_T, intent(in) :: num_materials_in
INTEGER_T, intent(in) :: num_species_var_in
INTEGER_T, intent(in) :: probe_constrain
INTEGER_T, intent(in) :: ireverse
INTEGER_T, intent(in) :: iten
REAL_T, intent(in) :: xI(SDIM)
REAL_T, intent(in) :: cur_time
REAL_T, intent(in) :: prev_time
REAL_T, intent(in) :: dt
REAL_T, intent(inout) :: TI
REAL_T, intent(inout) :: YI
INTEGER_T, intent(inout) :: user_override_TI_YI
REAL_T, intent(in) :: molar_mass(num_materials_in)
REAL_T, intent(in) :: species_molar_mass(num_species_var_in)
REAL_T, intent(in) :: ksrc_predict
REAL_T, intent(in) :: kdst_predict
REAL_T, intent(in) :: ksrc_physical
REAL_T, intent(in) :: kdst_physical
REAL_T, intent(in) :: T_probe_src
REAL_T, intent(in) :: T_probe_dst
INTEGER_T, intent(in) :: probe_ok_gradient_src
INTEGER_T, intent(in) :: probe_ok_gradient_dst
REAL_T, intent(in) :: LL
REAL_T, intent(in) :: dxprobe_src
REAL_T, intent(in) :: dxprobe_dst

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
  mdot, & ! intent(out)
  mdot_override, & ! intent(inout)
  ksrc_derived, &
  kdst_derived, &
  ksrc_physical, &
  kdst_physical, &
  T_probe_src, &
  T_probe_dst, &
  probe_ok_gradient_src, &
  probe_ok_gradient_dst, &
  TI, &
  LL, &
  dxprobe_src, &
  dxprobe_dst)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: interface_mass_transfer_model
INTEGER_T, intent(in) :: num_materials_in
INTEGER_T, intent(in) :: num_species_var_in
INTEGER_T, intent(in) :: ispec
INTEGER_T, intent(in) :: im_source
INTEGER_T, intent(in) :: im_dest
REAL_T, intent(in) :: xI(SDIM)
REAL_T, intent(in) :: TI
REAL_T, intent(in) :: molar_mass(num_materials_in)
REAL_T, intent(in) :: species_molar_mass(num_species_var_in+1)
REAL_T, intent(out) :: mdot
INTEGER_T, intent(inout) :: mdot_override
REAL_T, intent(in) :: ksrc_derived
REAL_T, intent(in) :: kdst_derived
REAL_T, intent(in) :: ksrc_physical
REAL_T, intent(in) :: kdst_physical
REAL_T, intent(in) :: T_probe_src
REAL_T, intent(in) :: T_probe_dst
INTEGER_T, intent(in) :: probe_ok_gradient_src
INTEGER_T, intent(in) :: probe_ok_gradient_dst
REAL_T, intent(in) :: LL
REAL_T, intent(in) :: dxprobe_src
REAL_T, intent(in) :: dxprobe_dst
REAL_T DTsrc,DTdst,mdotsrc,mdotdst,mdotsum

mdot_override=0

if (interface_mass_transfer_model.eq.0) then
 ! do nothing
else if (interface_mass_transfer_model.eq.999) then
 mdot_override=1
 DTsrc=T_probe_src-TI
 DTdst=T_probe_dst-TI
 if (probe_ok_gradient_src.eq.1) then
  ! do nothing
 else if (probe_ok_gradient_src.eq.0) then
  DTsrc=zero
 else
  print *,"probe_ok_gradient_src invalid"
  stop
 endif
 if (probe_ok_gradient_dst.eq.1) then
  ! do nothing
 else if (probe_ok_gradient_dst.eq.0) then
  DTdst=zero
 else
  print *,"probe_ok_gradient_dst invalid"
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
  probe_ok_gradient_src, &
  probe_ok_gradient_dst, &
  dxprobe_src, &
  dxprobe_dst, &
  LL, &
  num_materials_in, &
  num_species_var_in)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: interface_mass_transfer_model
INTEGER_T, intent(in) :: num_materials_in
INTEGER_T, intent(in) :: num_species_var_in
INTEGER_T, intent(in) :: ireverse
INTEGER_T, intent(in) :: iten
REAL_T, intent(in) :: molar_mass(num_materials_in)
REAL_T, intent(in) :: species_molar_mass(num_species_var_in)
REAL_T, intent(in) :: k_model_predict(2) ! src,dst
REAL_T, intent(inout) :: k_model_correct(2) ! src,dst
REAL_T, intent(in) :: k_physical_base(2) ! src, dst
REAL_T, intent(in) :: T_probe_src
REAL_T, intent(in) :: T_probe_dst
INTEGER_T, intent(in) :: probe_ok_gradient_src
INTEGER_T, intent(in) :: probe_ok_gradient_dst
REAL_T, intent(in) :: LL
REAL_T, intent(in) :: dxprobe_src
REAL_T, intent(in) :: dxprobe_dst

end subroutine STUB_K_EFFECTIVE

subroutine STUB_reference_wavelen(wavelen)
use probcommon_module
IMPLICIT NONE
REAL_T, intent(inout) :: wavelen
REAL_T :: default_wavelen
INTEGER_T :: dir_local

 default_wavelen=zero
 do dir_local=1,SDIM
  if (gravity_dir.ne.dir_local) then
   default_wavelen=max(default_wavelen,problen_array(dir_local))
  endif
 enddo

 if ((wavelen.gt.zero).and. &
     (wavelen.le.default_wavelen*(one+VOFTOL))) then
  ! do nothing
 else
  print *,"input wavelen out of range"
  print *,"wavelen (input) = ",wavelen
  print *,"default_wavelen (max wavelen) = ",default_wavelen
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

INTEGER_T, intent(in) :: iten
REAL_T, intent(in) :: time,temperature
REAL_T, intent(in) :: xpos(SDIM)
REAL_T, intent(inout) :: tension

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

INTEGER_T, intent(in) :: iten
REAL_T, intent(in) :: temperature
REAL_T, intent(inout) :: latent_heat ! always positive

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

end module STUB_module
