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

! probtype==55 
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



subroutine STUB_hydro_pressure_density( &
  xpos,rho,pres,from_boundary_hydrostatic)
use probcommon_module
IMPLICIT NONE

REAL_T, intent(in) :: xpos(SDIM)
REAL_T, intent(inout) :: rho
REAL_T, intent(inout) :: pres
INTEGER_T, intent(in) :: from_boundary_hydrostatic

pres=zero

return
end subroutine STUB_hydro_pressure_density

subroutine STUB_ASSIMILATE( &
  assimilate_in,assimilate_out,i,j,k,cell_flag,data_dir)
use probcommon_module
IMPLICIT NONE

type(assimilate_parm_type), intent(in) :: assimilate_in
type(assimilate_out_parm_type), intent(inout) :: assimilate_out
INTEGER_T, intent(in) :: i,j,k,cell_flag,data_dir

return
end subroutine STUB_ASSIMILATE


subroutine STUB_SUMINT(GRID_DATA_IN,increment_out,nsum)
use probcommon_module_types
use probcommon_module

INTEGER_T, intent(in) :: nsum
type(user_defined_sum_int_type), intent(in) :: GRID_DATA_IN
REAL_T, intent(out) :: increment_out(nsum)
INTEGER_T :: i,j,k

i=GRID_DATA_IN%igrid
j=GRID_DATA_IN%jgrid
k=GRID_DATA_IN%kgrid

end subroutine STUB_SUMINT


end module STUB_module
