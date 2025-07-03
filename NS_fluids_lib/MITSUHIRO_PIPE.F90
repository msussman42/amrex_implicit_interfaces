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

! probtype==41 (see run2d/inputs.MITSUHIRO_block_ice_melt)
module MITSUHIRO_PIPE_module
use amrex_fort_module, only : amrex_real

implicit none                   

real(amrex_real) :: DEF_VAPOR_GAMMA

contains

subroutine inletpipedist(x,y,z,dist)
use global_utility_module
use global_distance_module

IMPLICIT NONE

integer im
real(amrex_real), INTENT(out) :: dist(num_materials)
real(amrex_real), INTENT(in) :: x,y,z
real(amrex_real) ht,rr,initial_time
integer im_solid_pipe

if (SDIM.eq.2) then
 if (abs(z-y).le.EPS2) then
  ! do nothing
 else
  print *,"z<>y in inletpipedist"
  stop
 endif
endif

im_solid_pipe=im_solid_primary()

initial_time=zero

if (probtype.eq.41) then
 ! do nothing
else
 print *,"probtype invalid for inlet pipedist"
 stop
endif

if ((axis_dir.lt.0).or.(axis_dir.gt.5)) then
 print *,"axis dir invalid for pipe problem"
 stop
endif

if (axis_dir.eq.0) then

 if (SDIM.eq.2) then
  dist(1)=xblob+radblob*cos(two*Pi*y/yblob)-x
 else if (SDIM.eq.3) then
  rr=sqrt(y**2+z**2)
  dist(1)=rr-zblob-radblob*sin(two*Pi*x/xblob)
 else
  print *,"dimension bust"
  stop
 endif

else if ((axis_dir.eq.1).or. &
         (axis_dir.eq.2).or. &
         (axis_dir.eq.3)) then

 if (SDIM.eq.3) then
  dist(1)=z-zblob-radblob*sin(two*Pi*x/xblob) 
 else if (SDIM.eq.2) then
  dist(1)=abs(x)-xblob-radblob*sin(two*Pi*y/yblob)
  if (axis_dir.eq.2) then ! gas in the middle
   ht=radblob2+radblob*sin(two*Pi*y/yblob)
   if (x.gt.xblob) then
    dist(1)=x-(xblob+ht)
   else
    dist(1)=(xblob-ht)-x
   endif
  endif
  if (axis_dir.eq.3) then ! blowout problem
   dist(1)=y
  endif
 else
  print *,"dimension bust"
  stop
 endif

else if (axis_dir.eq.4) then

 dist(1)=abs(x)-xblob-radblob*sin(two*Pi*y/yblob)
 dist(1)=-dist(1)

else if (axis_dir.eq.5) then

 dist(1)=-99999.0

else 
 print *,"axis dir invalid"
 stop
endif

dist(2)=-dist(1)
do im=3,num_materials
 dist(im)=-99999.0
enddo

  ! axis_dir=4 is the comparison with Linear Stability Analysis
if (axis_dir.ne.4) then
 if ((im_solid_pipe.lt.1).or. &
     (im_solid_pipe.gt.num_materials)) then
  print *,"im_solid_pipe invalid 2"
  stop
 endif
   ! in inlet_pipe_dist; positive in solid
 call materialdistsolid(x,y,z,dist(im_solid_pipe),initial_time, &
  im_solid_pipe)
endif

return
end subroutine inletpipedist

subroutine get_pipe_velocity(xsten,nhalf,dx,bfact,vel,time)
use global_utility_module

IMPLICIT NONE

integer, INTENT(in) :: nhalf,bfact
real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
real(amrex_real) x,y,z,r
real(amrex_real), INTENT(in) :: time
real(amrex_real), INTENT(in) :: dx(SDIM)
real(amrex_real) cenbc(num_materials,SDIM)
real(amrex_real), INTENT(out) :: vel(SDIM)
real(amrex_real) VOF(num_materials)
integer dir2
integer im_solid_pipe
real(amrex_real) x_vel,y_vel,z_vel

im_solid_pipe=im_solid_primary()

if (nhalf.lt.3) then
 print *,"nhalf invalid get pipe velocity"
 stop
endif
if (bfact.lt.1) then
 print *,"bfact invalid200"
 stop
endif

if (probtype.eq.41) then
 ! do nothing
else 
 print *,"probtype invalid in get pipe velocity"
 stop
endif

x=xsten(0,1)
y=xsten(0,2)
z=xsten(0,SDIM)

if (SDIM.eq.2) then
 if (abs(z-y).le.EPS2) then
  ! do nothing
 else
  print *,"z<>y line 4073"
  stop
 endif
endif

if ((axis_dir.lt.0).or.(axis_dir.gt.5)) then
 print *,"get_pipe_velocity: axis dir invalid for pipe problem"
 stop
endif

do dir2=1,SDIM
 vel(dir2)=zero
enddo

call get_pipe_vfrac(xsten,nhalf,dx,bfact,VOF,cenbc)  

  ! axis_dir=4 is LSA comparison
if (axis_dir.eq.4) then

 if (VOF(1).gt.zero) then
  if (time.gt.zero) then
   vel(SDIM)=vinletgas
  else if (time.eq.zero) then
   vel(SDIM)=yblob4
  else
   print *,"time invalid in get pipe velocity"
   print *,"time= ",time
   print *,"vof1,vof2 ",VOF(1),VOF(2)
   print *,"axis_dir ",axis_dir
   stop
  endif
 else if ((VOF(1).eq.zero).and.(VOF(2).gt.zero)) then
  if (time.gt.zero) then
   vel(SDIM)=advbot
  else if (time.eq.zero) then
   vel(SDIM)=xblob4
  else
   print *,"time invalid in get pipe velocity"
   print *,"time= ",time
   print *,"vof1,vof2 ",VOF(1),VOF(2)
   print *,"axis_dir ",axis_dir
   stop
  endif
 endif
else if (axis_dir.eq.3) then
 if ((im_solid_pipe.lt.1).or. &
     (im_solid_pipe.gt.num_materials)) then
  print *,"im_solid_pipe invalid 2.9"
  stop
 endif
 if (VOF(im_solid_pipe).gt.zero) then
  vel(SDIM)=zero
 else if (z.ge.VOFTOL*dx(SDIM)) then !y=z if 2D
  if (time.gt.zero) then
   vel(SDIM)=vinletgas
  else if (time.eq.zero) then
   vel(SDIM)=yblob4
  else
   print *,"time invalid in get pipe velocity"
   print *,"time= ",time
   print *,"vof1,vof2 ",VOF(1),VOF(2)
   print *,"axis_dir ",axis_dir
   stop
  endif
 else if (z.le.VOFTOL*dx(SDIM)) then ! y=z if 2D
  if (time.gt.zero) then
   vel(SDIM)=advbot
  else if (time.eq.zero) then
   vel(SDIM)=xblob4
  else
   print *,"time invalid in get pipe velocity"
   print *,"time= ",time
   print *,"vof1,vof2 ",VOF(1),VOF(2)
   print *,"axis_dir ",axis_dir
   stop
  endif
 endif

else if ((axis_dir.eq.0).or. &
         (axis_dir.eq.1).or. &
         (axis_dir.eq.2)) then

 if ((im_solid_pipe.lt.1).or. &
     (im_solid_pipe.gt.num_materials)) then
  print *,"im_solid_pipe invalid 3"
  stop
 endif
 if (VOF(im_solid_pipe).gt.zero) then
  vel(SDIM)=zero
 else if (VOF(1).gt.zero) then
  if (time.gt.zero) then
   vel(SDIM)=vinletgas
  else if (time.eq.zero) then
   vel(SDIM)=yblob4
  else
   print *,"time invalid in get pipe velocity"
   print *,"time= ",time
   print *,"vof1,vof2 ",VOF(1),VOF(2)
   print *,"axis_dir ",axis_dir
   stop
  endif
 else if ((VOF(1).eq.zero).and.(VOF(2).gt.zero)) then
  if (time.gt.zero) then
   vel(SDIM)=advbot
  else if (time.eq.zero) then
   vel(SDIM)=xblob4
  else
   print *,"time invalid in get pipe velocity"
   print *,"time= ",time
   print *,"vof1,vof2 ",VOF(1),VOF(2)
   print *,"axis_dir ",axis_dir
   stop
  endif
 endif

 ! above: axis_dir=0,1,2 (and above that 3,4)

else if (axis_dir.eq.5) then

 if (SDIM.eq.3) then
  ! do nothing
 else if (SDIM.eq.2) then
  z=zero
 else
  print *,"dimension bust"
  stop
 endif
 r = sqrt(y*y+z*z)

 if (r.gt.radblob) then
  x_vel=zero
 else 
  x_vel=adv_vel*r*(radblob-r)
 endif
 y_vel=zero
 z_vel=zero
 if((z-0.3)**2+(y-0.3)**2+(x-0.5)**2.le.0.4)then
   x_vel = .3d0
   y_vel = .4d0
   z_vel = -.2d0
 else if((z+0.3)**2+(y-1.7)**2+(x-5.5)**2.le.0.4)then
   x_vel = -.3d0
   y_vel = -.4d0
   z_vel = .2d0
 else if((z-0.9)**2+(y)**2+(x-1.5)**2.le.1.)then
   x_vel = -x_vel
   y_vel = -y_vel
   z_vel = -z_vel 
 else if((z+0.9)**2+(y)**2+(x-3.5)**2.le.1.)then
   x_vel = -y_vel
   y_vel = -x_vel
   z_vel = -z_vel 
 else if((z)**2+(y)**2+(x-2.5)**2.le.0.1)then
   x_vel = -z_vel
   y_vel = -y_vel
   z_vel = -x_vel 
 endif
 vel(1)=x_vel
 vel(2)=y_vel
 if (SDIM.eq.3) then
  vel(SDIM)=z_vel
 endif
else
 print *,"axis_dir invalid"
 stop
endif

return
end subroutine get_pipe_velocity





  ! do any initial preparation needed
subroutine INIT_MITSUHIRO_PIPE_MODULE()
IMPLICIT NONE

  DEF_VAPOR_GAMMA =  1.666666667D0

return
end subroutine INIT_MITSUHIRO_PIPE_MODULE

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

! ice + water thickness = total height = radblob
! radius of ice block = radblob4 (if radblob4==0, then default to radblob/2)
! water thickness = radblob3  usually radblob3 << radlob (initial water
! is "seed" for starting the melting process)
!   ---------
!   | ice   |
!   ---------
!   |water  |
!-------------------
!   substrate

 ! fluids tessellate the domain, solids are immersed. 
subroutine MITSUHIRO_PIPE_LS(x,t,LS,nmat)
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
   print *,"nmat invalid: ",nmat
   stop
  endif

if (probtype.eq.41) then

 if (axis_dir.eq.5) then
  if (SDIM.eq.2) then
   dist=radblob-sqrt((y-yblob)**2)
  else if (SDIM.eq.3) then
   dist=radblob-sqrt((y-yblob)**2+(z-zblob)**2) 
  else
   print *,"dimension bust"
   stop
  endif 
 else if ((axis_dir.eq.0).and.(SDIM.eq.3)) then
! x is free stream direction
  if (levelrz.eq.COORDSYS_CARTESIAN) then 
   dist=zblob2-sqrt(y**2+z**2)
  else
   print *,"levelrz invalid MITSUHIRO_PIPE"
   stop
  endif

! pipe problem  soliddist: dist>0 fluid
 else if (SDIM.eq.2) then  

   ! axis_dir=4 comparison with LSA
  if (axis_dir.eq.4) then
   dist=99999.0
  else
   pipexlo=problox
   pipexhi=probhix
   if ((axis_dir.eq.1).or.(axis_dir.eq.2)) then
    pipexlo=zero
    pipexhi=two*radblob3
   endif
 
   if (x.lt.pipexlo) then
    dist=pipexlo-x
   else if (x.gt.pipexhi) then
    dist=x-pipexhi
   else
    dist=-min( x-pipexlo, pipexhi-x )
   endif
   dist=-dist
  endif  ! axis_dir<> 4

! pipe - vapordist 2D or 3D
       else if (probtype.eq.41) then
        call inletpipedist(x,y,z,distbatch)
        dist=distbatch(1)


else
 print *,"probtype invalid: ",probtype
 stop
endif

return
end subroutine MITSUHIRO_PIPE_LS

! initial velocity is zero
subroutine MITSUHIRO_PIPE_VEL(x,t,LS,VEL,velsolid_flag,dx,nmat)
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
 print *,"velsolid_flag invalid: ",velsolid_flag
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
end subroutine MITSUHIRO_PIPE_VEL


! this routine used if pressure boundary conditions are prescribed,
! since only top wall is "outflow" (outflow in quotes since ice shrinks when
! melting), and flow is incompressible, ok to make the top wall pressure zero.
subroutine MITSUHIRO_PIPE_PRES(x,t,LS,PRES,nmat)
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
end subroutine MITSUHIRO_PIPE_PRES



subroutine MITSUHIRO_PIPE_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
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
    (probtype.eq.41)) then
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
end subroutine MITSUHIRO_PIPE_STATE

 ! dir=1..sdim  side=1..2
subroutine MITSUHIRO_PIPE_LS_BC(xwall,xghost,t,LS, &
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
 call MITSUHIRO_PIPE_LS(xghost,t,LS,nmat)
else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine MITSUHIRO_PIPE_LS_BC


 ! dir=1..sdim  side=1..2 veldir=1..sdim
subroutine MITSUHIRO_PIPE_VEL_BC(xwall,xghost,t,LS, &
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

 call MITSUHIRO_PIPE_VEL(xghost,t,LS,local_VEL,velsolid_flag,dx,nmat)
 VEL=local_VEL(veldir)

else
 print *,"dir,side, or veldir invalid"
 stop
endif

return
end subroutine MITSUHIRO_PIPE_VEL_BC

 ! dir=1..sdim  side=1..2
subroutine MITSUHIRO_PIPE_PRES_BC(xwall,xghost,t,LS, &
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

 call MITSUHIRO_PIPE_PRES(xghost,t,LS,PRES,nmat)

else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine MITSUHIRO_PIPE_PRES_BC

 ! dir=1..sdim  side=1..2
subroutine MITSUHIRO_PIPE_STATE_BC(xwall,xghost,t,LS, &
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
 call MITSUHIRO_PIPE_STATE(xghost,t,LS,local_STATE,local_bcflag, &
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
end subroutine MITSUHIRO_PIPE_STATE_BC

subroutine MITSUHIRO_PIPE_HEATSOURCE(im,VFRAC,time,x, &
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

if ((num_materials.eq.4).and.(probtype.eq.41)) then
 heat_source=zero
else
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine MITSUHIRO_PIPE_HEATSOURCE

subroutine MITSUHIRO_PIPE_VARIABLE_SURFACE_TENSION( &
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

 if (probtype.eq.41) then
   if (radblob9.eq.zero) then
    !do nothing
   else if (radblob9.gt.zero) then
    if (time.eq.zero) then
     tension=fort_tension_init(iten)
    else if ((time.ge.zero).and. &
             (time.le.radblob9)) then
     tension=fort_tension_init(iten)
    else if (time.ge.two*radblob9) then
     tension=fort_tension(iten)
    else if ((time.ge.radblob9).and. &
             (time.le.two*radblob9)) then
     theta=(time-radblob9)/radblob9
     tension=(one-theta)*fort_tension_init(iten)+ &
             theta*fort_tension(iten)
    else
     print *,"time invalid: ",time
     stop
    endif
   else
    print *,"radblob9 invalid: ",radblob9
    stop
   endif
 else
  print *,"unexpected probtype: ",probtype
  stop
 endif

end subroutine MITSUHIRO_PIPE_VARIABLE_SURFACE_TENSION

end module MITSUHIRO_PIPE_module
