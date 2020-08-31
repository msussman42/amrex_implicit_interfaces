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
module GENERAL_PHASE_CHANGE_module

implicit none                   

REAL_T :: DEF_VAPOR_GAMMA

contains

  ! do any initial preparation needed
subroutine INIT_GENERAL_PHASE_CHANGE_MODULE()
IMPLICIT NONE

  DEF_VAPOR_GAMMA =  1.666666667D0

return
end subroutine INIT_GENERAL_PHASE_CHANGE_MODULE


 ! this is velocity boundary condition at the top of the domain.  
subroutine acoustic_pulse_bc(time,vel_pulse,xsten,nhalf,for_dt)
IMPLICIT NONE

INTEGER_T, intent(in) :: for_dt
REAL_T, intent(in) :: time
REAL_T, intent(out) :: vel_pulse
INTEGER_T, intent(in) :: nhalf
REAL_T, intent(in) :: xsten(-nhalf:nhalf,SDIM)
REAL_T x,y,z

x=xsten(0,1)
y=xsten(0,2)
z=xsten(0,SDIM)

if ((time.ge.zero).and.(time.le.1.0e+20)) then
 ! do nothing
else
 print *,"time invalid"
 stop
endif

if (abs(x)+abs(y)+abs(z).le.1.0e+20) then
 ! do nothing
else
 print *,"x,y, or z invalid"
 stop
endif
if ((for_dt.eq.0).or.(for_dt.eq.1)) then
 ! do nothing
else
 print *,"for_dt invalid"
 stop
endif

if (probtype.eq.55) then

 if (axis_dir.eq.7) then

  vel_pulse=zero
 
  if (1.eq.1) then 

   if (for_dt.eq.1) then
    vel_pulse=200.0
   else if (for_dt.eq.0) then 
    if ((time.ge.5.0e-8).and.(time.le.5.0e-7)) then
     vel_pulse=-200.0
    endif
   else
    print *,"for_dt invalid"
    stop
   endif

  endif

 else
  print *,"unexpected axis_dir"
  stop
 endif

else
 print *,"unexpected probtype"
 stop
endif

return
end subroutine acoustic_pulse_bc 

 ! density_at_depth previously initialized by:
 ! init_density_at_depth() 
 ! called from: FORT_DENCOR, general_hydrostatic_pressure_density,
 ! boundary_hydrostatic, EOS_air_rho2, EOS_air_rho2_ADIABAT,
 ! SOUNDSQR_air_rho2, EOS_error_ind, presBDRYCOND, FORT_INITDATA 
subroutine GENERAL_PHASE_CHANGE_hydro_pressure_density( &
  xpos,rho,pres)
IMPLICIT NONE

REAL_T, intent(in) :: xpos(SDIM)
REAL_T, intent(out) :: rho
REAL_T, intent(out) :: pres
REAL_T denfree,zfree
REAL_T den_top,z_top
REAL_T z_at_depth

 if (probtype.eq.55) then
  if ((axis_dir.ge.0).and.(axis_dir.le.5)) then
   ! do nothing (compressible drop) ??
   FIX ME (boundary_hydrostatic)
  endif

  if (axis_dir.eq.7) then

   if (SDIM.eq.2) then
    z_at_depth=probloy
    z_top=probhiy
   else if (SDIM.eq.3) then
    z_at_depth=probloz
    z_top=probhiz
   else
    print *,"dimension bust GENERAL_PHASE_CHANGE_hydro_pressure_density"
    stop
   endif

   if (z_at_depth.ne.zero) then
    print *,"z_at_depth must be 0 for compressible boiling problem"
    stop
   endif

   den_top=fort_denconst(1)

   if (xpos(SDIM).gt.z_top) then
    rho=den_top ! atmos pressure at top of domain
   else
    rho= &
     ((density_at_depth-den_top)/ &
      (z_at_depth-z_top))*(xpos(SDIM)-z_top)+den_top
   endif
   call EOS_tait_ADIABATIC_rhohydro(rho,pres)


  else if (fort_material_type(1).eq.13) then

   denfree=fort_denconst(1)
   if (SDIM.eq.2) then
    zfree=probhiy
    z_at_depth=probloy
   else if (SDIM.eq.3) then
    zfree=probhiz
    z_at_depth=probloz
   else
    print *,"dimension bust"
    stop
   endif

   ! density_at_depth is found so that
   ! (p(density_at_depth)-p(rho_0))/(rho_0 (z_at_depth-zfree))=g
   !
   if (xpos(SDIM).gt.zfree) then
    rho=denfree
   else
    rho= &
     ((density_at_depth-denfree)/ &
      (z_at_depth-zfree))*(xpos(SDIM)-zfree)+denfree
   endif
   call EOS_tait_ADIABATIC_rhohydro(rho,pres)
  else
   print *,"axis_dir invalid GENERAL_PHASE_CHANGE_hydro_pressure_density"
   stop
  endif
 else
  print *,"probtype invalid GENERAL_PHASE_CHANGE_hydro_pressure_density"
  stop
 endif

return
end subroutine GENERAL_PHASE_CHANGE_hydro_pressure_density

subroutine GENERAL_PHASE_CHANGE_CFL_HELPER(time,dir,uu,dx)
INTEGER_T, intent(in) :: dir
REAL_T, intent(in) :: time
REAL_T, intent(inout) :: uu
REAL_T, intent(in) :: dx(SDIM)

INTEGER_T dir2
INTEGER_T side
REAL_T utest,uscale
REAL_T xsten_dummy(-1:1,SDIM)
INTEGER_T for_dt
INTEGER_T nhalf

nhalf=1

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
 print *,"dir invalid GENERAL_PHASE_CHANGE_CFL_HELPER"
 stop
endif

if (probtype.eq.55) then
 if (axis_dir.eq.7) then
  do dir2=1,SDIM
  do side=-nhalf,nhalf
   xsten_dummy(side,dir2)=dx(dir2)*half*side
  enddo
  enddo
  for_dt=1
  call acoustic_pulse_bc(time,utest,xsten_dummy,nhalf,for_dt)
  uu=max(abs(uu),abs(utest))
 endif
else
 print *,"unexpected probtype"
 stop
endif

return
end subroutine GENERAL_PHASE_CHANGE_CFL_HELPER


! Phi>0 in the solid
subroutine BOTTOM_substrateLS(x,Phi) 
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

end subroutine BOTTOM_substrateLS

! ice + water thickness = radblob
! water thickness = radblob3  usually radblob3 << radlob (initial water
! is "seed" for starting the melting process)
!   ---------
!   | water |
!   ---------
!   | ice   |
!-------------------
!   substrate

 ! fluids tessellate the domain, solids are immersed. 
subroutine GENERAL_PHASE_CHANGE_LS(x,t,LS,nmat)
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

if ((num_materials.eq.4).and.(probtype.eq.2001)) then

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

  ! water thickness=radblob3
  !dist<0 inside the square
  !water above the ice
 if (SDIM.eq.2) then
  call squaredist(x(1),x(2),xblob-half*radblob,xblob+half*radblob, &
     substrate_height+radblob-radblob3,substrate_height+radblob,LS(1))
  LS(1)=-LS(1) ! water
  !ice (is below water)
  call squaredist(x(1),x(2),xblob-half*radblob,xblob+half*radblob, &
    -substrate_height-radblob,substrate_height+radblob-radblob3,LS(3))
  LS(3)=-LS(3)

  !air; dist<0 inside the square
  call squaredist(x(1),x(2),xblob-half*radblob,xblob+half*radblob, &
    -substrate_height-radblob,substrate_height+radblob,LS(2))

 else if (SDIM.eq.3) then

  !dist<0 inside the square
  !water above the ice
  call cubedist(xblob-half*radblob,xblob+half*radblob, &
     yblob-half*radblob,yblob+half*radblob, &
     substrate_height+radblob-radblob3,substrate_height+radblob, &
     x(1),x(2),x(SDIM),LS(1))
  LS(1)=-LS(1)
  !ice
  call cubedist(xblob-half*radblob,xblob+half*radblob, &
     yblob-half*radblob,yblob+half*radblob, &
     -substrate_height-radblob,substrate_height+radblob-radblob3, &
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

 call BOTTOM_substrateLS(x,LS(4))
else
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine GENERAL_PHASE_CHANGE_LS

! initial velocity is zero
subroutine GENERAL_PHASE_CHANGE_VEL(x,t,LS,VEL,velsolid_flag,dx,nmat)
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
end subroutine GENERAL_PHASE_CHANGE_VEL


! These next routines only used for compressible materials.
!***********************************************
! compressible material functions for (ns.material_type = 24)
subroutine EOS_GENERAL_PHASE_CHANGE(rho,internal_energy,pressure, &
  imattype,im)
 IMPLICIT NONE
 INTEGER_T, intent(in) :: imattype,im
 REAL_T, intent(in) :: rho
 REAL_T, intent(in) :: internal_energy
 REAL_T, intent(out) :: pressure

 if (imattype.eq.24) then
  pressure=rho*(DEF_VAPOR_GAMMA-1.0D0)*internal_energy
 else
  print *,"imattype invalid EOS_GENERAL_PHASE_CHANGE"
  stop
 endif

 return
end subroutine EOS_GENERAL_PHASE_CHANGE

subroutine SOUNDSQR_GENERAL_PHASE_CHANGE(rho,internal_energy,soundsqr, &
  imattype,im)
 IMPLICIT NONE
 INTEGER_T, intent(in) :: imattype,im
 REAL_T, intent(in) :: rho
 REAL_T, intent(in) :: internal_energy
 REAL_T, intent(out) :: soundsqr
 REAL_T pressure

 if (imattype.eq.24) then
  call EOS_GENERAL_PHASE_CHANGE(rho,internal_energy,pressure,imattype,im)
  soundsqr=DEF_VAPOR_GAMMA*pressure/rho
 else
  print *,"imattype invalid SOUNDSQR_GENERAL_PHASE_CHANGE"
  stop
 endif

 return
end subroutine SOUNDSQR_GENERAL_PHASE_CHANGE


subroutine INTERNAL_GENERAL_PHASE_CHANGE(rho,temperature,local_internal_energy, &
  imattype,im)
 use global_utility_module
 IMPLICIT NONE
 INTEGER_T, intent(in) :: imattype,im
 REAL_T, intent(in) :: rho
 REAL_T, intent(in) :: temperature 
 REAL_T, intent(out) :: local_internal_energy

 call INTERNAL_default(rho,temperature,local_internal_energy, &
        imattype,im)

 return
end subroutine INTERNAL_GENERAL_PHASE_CHANGE

subroutine TEMPERATURE_GENERAL_PHASE_CHANGE(rho,temperature,internal_energy, &
  imattype,im)
 use global_utility_module
 IMPLICIT NONE
 INTEGER_T, intent(in) :: imattype,im
 REAL_T, intent(in) :: rho
 REAL_T, intent(out) :: temperature 
 REAL_T, intent(in) :: internal_energy

 call TEMPERATURE_default(rho,temperature,internal_energy, &
        imattype,im)

 return
end subroutine TEMPERATURE_GENERAL_PHASE_CHANGE


! this routine used if pressure boundary conditions are prescribed,
! since only top wall is "outflow" (outflow in quotes since ice shrinks when
! melting), and flow is incompressible, ok to make the top wall pressure zero.
subroutine GENERAL_PHASE_CHANGE_PRES(x,t,LS,PRES,nmat)
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
end subroutine GENERAL_PHASE_CHANGE_PRES



subroutine GENERAL_PHASE_CHANGE_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
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
    (num_state_material.eq.2).and. & ! density, temperature
    (probtype.eq.2001)) then
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
end subroutine GENERAL_PHASE_CHANGE_STATE

 ! dir=1..sdim  side=1..2
subroutine GENERAL_PHASE_CHANGE_LS_BC(xwall,xghost,t,LS, &
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
 call GENERAL_PHASE_CHANGE_LS(xghost,t,LS,nmat)
else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine GENERAL_PHASE_CHANGE_LS_BC


 ! dir=1..sdim  side=1..2 veldir=1..sdim
subroutine GENERAL_PHASE_CHANGE_VEL_BC(xwall,xghost,t,LS, &
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

 call GENERAL_PHASE_CHANGE_VEL(xghost,t,LS,local_VEL,velsolid_flag,dx,nmat)
 VEL=local_VEL(veldir)

else
 print *,"dir,side, or veldir invalid"
 stop
endif

return
end subroutine GENERAL_PHASE_CHANGE_VEL_BC

 ! dir=1..sdim  side=1..2
subroutine GENERAL_PHASE_CHANGE_PRES_BC(xwall,xghost,t,LS, &
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

 call GENERAL_PHASE_CHANGE_PRES(xghost,t,LS,PRES,nmat)

else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine GENERAL_PHASE_CHANGE_PRES_BC

 ! dir=1..sdim  side=1..2
subroutine GENERAL_PHASE_CHANGE_STATE_BC(xwall,xghost,t,LS, &
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
 call GENERAL_PHASE_CHANGE_STATE(xghost,t,LS,local_STATE,local_bcflag, &
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
end subroutine GENERAL_PHASE_CHANGE_STATE_BC

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
subroutine GENERAL_PHASE_CHANGE_HEATSOURCE(im,VFRAC,time,x, &
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
REAL_T :: MITSUHIRO_CUSTOM_TEMPERATURE

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
 
 ! set a hot temperature in the liquid
if ((num_materials.eq.4).and.(probtype.eq.2001)) then

 heat_source=zero
 if (VFRAC(1).ge.half) then ! in the liquid
  MITSUHIRO_CUSTOM_TEMPERATURE=fort_tempconst(1)
  heat_source=(MITSUHIRO_CUSTOM_TEMPERATURE-temp(1))* &
               den(1)*CV(1)/dt
 else if (VFRAC(1).le.half) then
  ! do nothing
 else
  print *,"VFRAC(1) invalid"
  stop
 endif 

else
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine GENERAL_PHASE_CHANGE_HEATSOURCE

end module GENERAL_PHASE_CHANGE_module
