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

! probtype==414 (see run2d/inputs.MITSUHIRO_block_ice_melt)
module MITSUHIRO_MELTING_module
use amrex_fort_module, only : amrex_real

implicit none                   

real(amrex_real) :: DEF_VAPOR_GAMMA

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
real(amrex_real), INTENT(in), dimension(SDIM) :: x !spatial coordinates
real(amrex_real), INTENT(out) :: Phi !LS dist, Phi>0 in the substrate

real(amrex_real) substrate_height

if (SDIM.eq.2) then
 if (abs(zblob2-yblob2).le.EPS14) then
  substrate_height=zblob2  ! substrate thickness
 else
  print *,"zblob2 or yblob2 invalid (they should be the same) 2D"
  print *,"zblob2,yblob2   = ",zblob2,yblob2
  stop
 endif
else if (SDIM.eq.3) then
 substrate_height=zblob2  ! substrate thickness
else
 print *,"dimension bust"
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

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(out) :: LS(nmat)
real(amrex_real) :: ice_vertical
real(amrex_real) :: substrate_height
real(amrex_real) :: box_xlo,box_xhi,box_ylo,box_yhi,box_zlo,box_zhi

  if (nmat.eq.num_materials) then
   ! do nothing
  else
   print *,"nmat invalid"
   stop
  endif

if (SDIM.eq.2) then
 if (abs(zblob2-yblob2).le.EPS14) then
  substrate_height=zblob2  ! substrate thickness
 else
  print *,"zblob2 or yblob2 invalid"
  stop
 endif
else if (SDIM.eq.3) then
 substrate_height=zblob2  ! substrate thickness
else
 print *,"dimension bust"
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
 if (abs(substrate_height-(ice_vertical-half*radblob)).gt.EPS2) then
  print *,"bottom of original water+ice block must coincide w/substrate"
  print *,"probtype=",probtype
  print *,"sdim= ",SDIM
  print *,"yblob=",yblob
  print *,"zblob=",zblob
  print *,"ice_vertical (yblob if 2d, zblob if 3d)=",ice_vertical
  print *,"substrate_height (zblob2) =",substrate_height
  print *,"radblob=",radblob
  print *,"radblob/2=",half*radblob
  print *,"ice_vertical-radblob/2 ",ice_vertical-half*radblob
  stop
 endif

  !dist<0 inside the square
  !water below the ice
 box_xlo=xblob-half*radblob
 box_xhi=xblob+half*radblob
 if (SDIM.eq.2) then
  box_ylo=-substrate_height-radblob3
  box_yhi=substrate_height+radblob3
  call squaredist(x(1),x(2),box_xlo,box_xhi, &
     box_ylo,box_yhi,LS(1))
  LS(1)=-LS(1) ! water
  !ice (is above water)
  box_ylo=substrate_height+radblob3
  box_yhi=substrate_height+radblob
  call squaredist(x(1),x(2),box_xlo,box_xhi, &
    box_ylo,box_yhi,LS(3))
  LS(3)=-LS(3)

  !air; dist<0 inside the square
  box_ylo=-substrate_height-radblob
  box_yhi=substrate_height+radblob
  call squaredist(x(1),x(2),box_xlo,box_xhi, &
    box_ylo,box_yhi,LS(2))

 else if (SDIM.eq.3) then
  box_ylo=yblob-half*radblob
  box_yhi=yblob+half*radblob

  !dist<0 inside the square
  !water below the ice
  box_zlo=-substrate_height-radblob3
  box_zhi=substrate_height+radblob3
  call cubedist(box_xlo,box_xhi, &
     box_ylo,box_yhi, &
     box_zlo,box_zhi, &
     x(1),x(2),x(SDIM),LS(1))
  LS(1)=-LS(1)
  !ice
  box_zlo=substrate_height+radblob3
  box_zhi=substrate_height+radblob
  call cubedist(box_xlo,box_xhi, &
     box_ylo,box_yhi, &
     box_zlo,box_zhi, &
     x(1),x(2),x(SDIM),LS(3))
  LS(3)=-LS(3)

  !air; dist<0 inside the square
  !air everywhere not ice or water.
  ! important: fluids tessellate domain, solids (i.e. substrate)
  ! are embedded.
  box_zlo=-substrate_height-radblob
  box_zhi=substrate_height+radblob
  call cubedist(box_xlo,box_xhi, &
     box_ylo,box_yhi, &
     box_zlo,box_zhi, &
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

return 
end subroutine MITSUHIRO_MELTING_VEL


! this routine used if pressure boundary conditions are prescribed,
! since only top wall is "outflow" (outflow in quotes since ice shrinks when
! melting), and flow is incompressible, ok to make the top wall pressure zero.
subroutine MITSUHIRO_MELTING_PRES(x,t,LS,PRES,nmat)
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
end subroutine MITSUHIRO_MELTING_PRES



subroutine MITSUHIRO_MELTING_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
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
  STATE(ibase+ENUM_DENVAR+1)=fort_denconst(im) ! density prescribed in the inputs file.
  if (t.eq.zero) then
   STATE(ibase+ENUM_TEMPERATUREVAR+1)=fort_initial_temperature(im) !initial temperature in inputs
  else if (t.gt.zero) then
   STATE(ibase+ENUM_TEMPERATUREVAR+1)=fort_tempconst(im)
  else
   print *,"t invalid"
   stop
  endif
   ! always assume Dirichlet boundary condition at zlo for temperature.
  call outside_temperature(t,x(1),x(2),x(SDIM),STATE(ibase+ENUM_TEMPERATUREVAR+1),im,bcflag)

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
end subroutine MITSUHIRO_MELTING_STATE

 ! dir=1..sdim  side=1..2
subroutine MITSUHIRO_MELTING_LS_BC(xwall,xghost,t,LS, &
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
 call MITSUHIRO_MELTING_STATE(xghost,t,LS,local_STATE,local_bcflag, &
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
end subroutine MITSUHIRO_MELTING_STATE_BC

subroutine MITSUHIRO_MELTING_HEATSOURCE(im,VFRAC,time,x, &
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

if ((num_materials.eq.4).and.(probtype.eq.414)) then
 heat_source=zero
else
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine MITSUHIRO_MELTING_HEATSOURCE

end module MITSUHIRO_MELTING_module
