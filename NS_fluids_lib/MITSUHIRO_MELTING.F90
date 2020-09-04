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

! probtype==414 (see run2d/inputs.MITSUHIRO_block_ice_melt)
module MITSUHIRO_MELTING_module

implicit none                   

REAL_T :: DEF_VAPOR_GAMMA

contains

  ! do any initial preparation needed
subroutine INIT_MITSUHIRO_MELTING_MODULE()
IMPLICIT NONE

  DEF_VAPOR_GAMMA =  1.666666667D0

return
end subroutine INIT_MITSUHIRO_MELTING_MODULE

! Phi>0 in the solid
subroutine MITSUHIRO_substrateLS(x,Phi) 
use probcommon_module
implicit none
REAL_T, intent(in), dimension(SDIM) :: x !spatial coordinates
REAL_T, intent(out) :: Phi !LS dist, Phi>0 in the substrate

REAL_T substrate_height

if (abs(zblob2-yblob2).le.1.0D-14) then
 substrate_height=zblob2  ! substrate thickness
else
 print *,"zblob2 or yblob2 invalid"
 stop
endif

if (abs(x(SDIM)).le.1.0D+20) then
 Phi=substrate_height-x(SDIM)
else
 print *,"x(SDIM) invalid"
 stop
endif

end subroutine MITSUHIRO_substrateLS

! ice + water thickness = radblob
! water thickness = radblob3  usually radblob3 << radlob (initial water
! is "seed" for starting the melting process)
!   ---------
!   | ice   |
!   ---------
!   |water  |
!-------------------
!   substrate

 ! fluids tessellate the domain, solids are immersed. 
subroutine MITSUHIRO_MELTING_LS(x,t,LS,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(out) :: LS(nmat)
REAL_T :: ice_vertical
REAL_T :: substrate_height

  if (nmat.eq.num_materials) then
   ! do nothing
  else
   print *,"nmat invalid"
   stop
  endif

if (abs(zblob2-yblob2).le.1.0D-14) then
 substrate_height=zblob2  ! substrate thickness
else
 print *,"zblob2 or yblob2 invalid"
 stop
endif

if ((num_materials.eq.4).and.(probtype.eq.414)) then

 if (SDIM.eq.2) then
  ice_vertical=yblob
 else if (SDIM.eq.3) then
  ice_vertical=zblob
 else
  print *,"dimension bust"
  stop
 endif
 if (abs(substrate_height-(ice_vertical-half*radblob)).gt.VOFTOL) then
  print *,"bottom of original water+ice block must coincide w/substrate"
  print *,"sdim= ",SDIM
  print *,"yblob=",yblob
  print *,"zblob=",zblob
  print *,"ice_vertical=",ice_vertical
  print *,"substrate_height=",substrate_height
  print *,"radblob=",radblob
  stop
 endif

  !dist<0 inside the square
  !water below the ice
 if (SDIM.eq.2) then
  call squaredist(x(1),x(2),xblob-half*radblob,xblob+half*radblob, &
     -substrate_height-radblob3,substrate_height+radblob3,LS(1))
  LS(1)=-LS(1) ! water
  !ice (is above water)
  call squaredist(x(1),x(2),xblob-half*radblob,xblob+half*radblob, &
    substrate_height+radblob3,substrate_height+radblob,LS(3))
  LS(3)=-LS(3)

  !air; dist<0 inside the square
  call squaredist(x(1),x(2),xblob-half*radblob,xblob+half*radblob, &
    -substrate_height-radblob,substrate_height+radblob,LS(2))

 else if (SDIM.eq.3) then

  !dist<0 inside the square
  !water below the ice
  call cubedist(xblob-half*radblob,xblob+half*radblob, &
     yblob-half*radblob,yblob+half*radblob, &
     -substrate_height-radblob3,substrate_height+radblob3, &
     x(1),x(2),x(SDIM),LS(1))
  LS(1)=-LS(1)
  !ice
  call cubedist(xblob-half*radblob,xblob+half*radblob, &
     yblob-half*radblob,yblob+half*radblob, &
     substrate_height+radblob3,substrate_height+radblob, &
     x(1),x(2),x(SDIM),LS(3))
  LS(3)=-LS(3)

  !air; dist<0 inside the square
  !air everywhere not ice or water.
  ! important: fluids tessellate domain, solids (i.e. substrate)
  ! are embedded.
  call cubedist(xblob-half*radblob,xblob+half*radblob, &
     yblob-half*radblob,yblob+half*radblob, &
     -substrate_height-radblob,substrate_height+radblob, &
     x(1),x(2),x(SDIM),LS(2))  ! air

 else
  print *,"dimension bust"
  stop
 endif

 call MITSUHIRO_substrateLS(x,LS(4))
else
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine MITSUHIRO_MELTING_LS

! initial velocity is zero
subroutine MITSUHIRO_MELTING_VEL(x,t,LS,VEL,velsolid_flag,dx,nmat)
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

return 
end subroutine MITSUHIRO_MELTING_VEL


! These next routines only used for compressible materials.
!***********************************************
! compressible material functions for (ns.material_type = 24)
subroutine EOS_MITSUHIRO_MELTING(rho,massfrac_var, &
  internal_energy,pressure, &
  imattype,im)
 use global_utility_module
 IMPLICIT NONE
 INTEGER_T, intent(in) :: imattype,im
 REAL_T, intent(in) :: rho
 REAL_T, intent(in) :: massfrac_var(num_species_var+1)
 REAL_T, intent(in) :: internal_energy
 REAL_T, intent(out) :: pressure

 call EOS_material_CORE(rho,massfrac_var, &
         internal_energy,pressure,imattype,im)

 return
end subroutine EOS_MITSUHIRO_MELTING

subroutine SOUNDSQR_MITSUHIRO_MELTING(rho,massfrac_var, &
  internal_energy,soundsqr, &
  imattype,im)
 use global_utility_module
 IMPLICIT NONE
 INTEGER_T, intent(in) :: imattype,im
 REAL_T, intent(in) :: rho
 REAL_T, intent(in) :: massfrac_var(num_species_var+1)
 REAL_T, intent(in) :: internal_energy
 REAL_T, intent(out) :: soundsqr

 call SOUNDSQR_material_CORE(rho,massfrac_var, &
   internal_energy,soundsqr, &
   imattype,im)

 return
end subroutine SOUNDSQR_MITSUHIRO_MELTING

subroutine INTERNAL_MITSUHIRO_MELTING(rho,massfrac_var, &
  temperature,local_internal_energy, &
  imattype,im)
 use global_utility_module
 IMPLICIT NONE
 INTEGER_T, intent(in) :: imattype,im
 REAL_T, intent(in) :: rho
 REAL_T, intent(in) :: massfrac_var(num_species_var+1)
 REAL_T, intent(in) :: temperature 
 REAL_T, intent(out) :: local_internal_energy

 call INTERNAL_material_CORE(rho,massfrac_var, &
   temperature,local_internal_energy, &
   imattype,im)

 return
end subroutine INTERNAL_MITSUHIRO_MELTING

subroutine TEMPERATURE_MITSUHIRO_MELTING(rho,massfrac_var, &
  temperature,internal_energy, &
  imattype,im)
 use global_utility_module
 IMPLICIT NONE
 INTEGER_T, intent(in) :: imattype,im
 REAL_T, intent(in) :: rho
 REAL_T, intent(in) :: massfrac_var(num_species_var+1)
 REAL_T, intent(out) :: temperature 
 REAL_T, intent(in) :: internal_energy

 call TEMPERATURE_material_CORE(rho,massfrac_var, &
     temperature,internal_energy, &
     imattype,im)

 return
end subroutine TEMPERATURE_MITSUHIRO_MELTING


! this routine used if pressure boundary conditions are prescribed,
! since only top wall is "outflow" (outflow in quotes since ice shrinks when
! melting), and flow is incompressible, ok to make the top wall pressure zero.
subroutine MITSUHIRO_MELTING_PRES(x,t,LS,PRES,nmat)
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
end subroutine MITSUHIRO_MELTING_PRES



subroutine MITSUHIRO_MELTING_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
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
if ((num_materials.eq.4).and. &
    (num_state_material.ge.3).and. & ! density, temperature, vapor spec
    (probtype.eq.414)) then
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
   ! always assume Dirichlet boundary condition at zlo for temperature.
  call outside_temperature(t,x(1),x(2),x(SDIM),STATE(ibase+2),im,bcflag)

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
end subroutine MITSUHIRO_MELTING_STATE

 ! dir=1..sdim  side=1..2
subroutine MITSUHIRO_MELTING_LS_BC(xwall,xghost,t,LS, &
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
 call MITSUHIRO_MELTING_LS(xghost,t,LS,nmat)
else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine MITSUHIRO_MELTING_LS_BC


 ! dir=1..sdim  side=1..2 veldir=1..sdim
subroutine MITSUHIRO_MELTING_VEL_BC(xwall,xghost,t,LS, &
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

 call MITSUHIRO_MELTING_VEL(xghost,t,LS,local_VEL,velsolid_flag,dx,nmat)
 VEL=local_VEL(veldir)

else
 print *,"dir,side, or veldir invalid"
 stop
endif

return
end subroutine MITSUHIRO_MELTING_VEL_BC

 ! dir=1..sdim  side=1..2
subroutine MITSUHIRO_MELTING_PRES_BC(xwall,xghost,t,LS, &
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

 call MITSUHIRO_MELTING_PRES(xghost,t,LS,PRES,nmat)

else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine MITSUHIRO_MELTING_PRES_BC

 ! dir=1..sdim  side=1..2
subroutine MITSUHIRO_MELTING_STATE_BC(xwall,xghost,t,LS, &
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
 call MITSUHIRO_MELTING_STATE(xghost,t,LS,local_STATE,local_bcflag, &
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
end subroutine MITSUHIRO_MELTING_STATE_BC

subroutine MITSUHIRO_MELTING_HEATSOURCE(im,VFRAC,time,x, &
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

if ((num_materials.eq.4).and.(probtype.eq.414)) then
 heat_source=zero
else
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine MITSUHIRO_MELTING_HEATSOURCE

end module MITSUHIRO_MELTING_module
