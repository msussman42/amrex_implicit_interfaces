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

! probtype==41 or probtype==3 (see run2d/inputs.growthrate.LSA)
module MITSUHIRO_PIPE_module
use amrex_fort_module, only : amrex_real

implicit none                   

real(amrex_real) :: DEF_VAPOR_GAMMA

contains

subroutine pipedist(x,y,z,dist,nmat)
use probcommon_module
use global_utility_module

IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(out) :: dist(nmat)
real(amrex_real), INTENT(in) :: x,y,z
real(amrex_real) ht,rr
real(amrex_real) pipexlo,pipexhi
integer im_solid_pipe

if (SDIM.eq.2) then
 if (abs(z-y).le.EPS2) then
  ! do nothing
 else
  print *,"z<>y in pipedist (x,y,z): ",x,y,z
  stop
 endif
else if (SDIM.eq.3) then
 !checking nothing
else
 print *,"dimension bust"
 stop
endif

im_solid_pipe=im_solid_primary()

if ((im_solid_pipe.eq.0).or. &
    (im_solid_pipe.eq.nmat)) then
 !do nothing
else
 print *,"expecting im_solid_pipe=0 or nmat: ",im_solid_pipe,nmat
 stop
endif

if ((probtype.eq.41).or.(probtype.eq.3)) then
 ! do nothing
else
 print *,"probtype invalid for pipedist: ",probtype
 stop
endif

if ((axis_dir.lt.0).or.(axis_dir.gt.5)) then
 print *,"axis dir invalid for pipe problem: ",axis_dir
 stop
endif

if (probtype.eq.3) then

 dist(1)=xblob+radblob*cos(two*Pi*y/yblob)-x

else if (axis_dir.eq.0) then

 print *,"this code is obsolete"
 stop

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

 dist(1)=-abs(x)+xblob+radblob*sin(two*Pi*y/yblob)

else if (axis_dir.eq.5) then ! all gas

 dist(1)=-99999.0

else 
 print *,"axis dir invalid: ",axis_dir
 stop
endif

dist(2)=-dist(1)

if (nmat.eq.2) then
 !do nothing
else if (nmat.eq.3) then
 dist(nmat)=-99999.0
else
 print *,"expecting num_materials=2 or 3 pipe_dist: ",nmat
 stop
endif

  ! axis_dir=4 is the comparison with Linear Stability Analysis
if (axis_dir.eq.4) then
 if (nmat.eq.2) then
  !do nothing
 else
  print *,"expecting nmat==2: ",nmat
  stop
 endif
else if ((axis_dir.eq.0).or. &
         (axis_dir.eq.1).or. & 
         (axis_dir.eq.2).or. & 
         (axis_dir.eq.3).or. & 
         (axis_dir.eq.5)) then

 if (im_solid_pipe.eq.3) then
  !do nothing
 else
  print *,"im_solid_pipe invalid 2: ",im_solid_pipe
  stop
 endif

 if (axis_dir.eq.5) then

  if (SDIM.eq.2) then
   dist(im_solid_pipe)=-radblob+sqrt((y-yblob)**2)
  else if (SDIM.eq.3) then
   dist(im_solid_pipe)=-radblob+sqrt((y-yblob)**2+(z-zblob)**2) 
  else
   print *,"dimension bust"
   stop
  endif 

 else if ((axis_dir.eq.0).and.(SDIM.eq.3)) then

! x is free stream direction
  if (levelrz.eq.COORDSYS_CARTESIAN) then 
   dist(im_solid_pipe)=-zblob2+sqrt(y**2+z**2)
  else
   print *,"levelrz invalid MITSUHIRO_PIPE"
   stop
  endif

 else   

  pipexlo=problox
  pipexhi=probhix
  if ((axis_dir.eq.1).or. &
      (axis_dir.eq.2).or. &
      (axis_dir.eq.3)) then
   pipexlo=zero
   pipexhi=two*radblob3
  endif
 
  if (x.lt.pipexlo) then
   dist(im_solid_pipe)=pipexlo-x
  else if (x.gt.pipexhi) then
   dist(im_solid_pipe)=x-pipexhi
  else
   dist(im_solid_pipe)=-min( x-pipexlo, pipexhi-x )
  endif

 endif

else
 print *,"axis_dir invalid in pipedist: ",axis_dir
 stop
endif

return
end subroutine pipedist

subroutine get_pipe_velocity(x,y,z,dx,vel,time,nmat)
use probcommon_module
use global_utility_module

IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: x,y,z
real(amrex_real) r,zlocal
real(amrex_real), INTENT(in) :: time
real(amrex_real), INTENT(in) :: dx(SDIM)
real(amrex_real), INTENT(out) :: vel(SDIM)
real(amrex_real) dist(nmat)
integer dir2
real(amrex_real) x_vel,y_vel,z_vel

if ((probtype.eq.41).or.(probtype.eq.3)) then
 ! do nothing
else 
 print *,"probtype invalid in get pipe velocity"
 stop
endif

if (SDIM.eq.2) then
 if (abs(z-y).le.EPS2) then
  ! do nothing
 else
  print *,"z<>y line 4073: ",x,y,z
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

call pipedist(x,y,z,dist,nmat)

  ! axis_dir=4 is LSA comparison
if (axis_dir.eq.4) then

 if (dist(1).ge.zero) then
  if (time.gt.zero) then
   vel(SDIM)=vinletgas
  else if (time.eq.zero) then
   vel(SDIM)=yblob4
  else
   print *,"time invalid in get pipe velocity"
   print *,"time= ",time
   print *,"dist= ",dist
   print *,"axis_dir ",axis_dir
   stop
  endif
 else if (dist(1).le.zero) then
  if (time.gt.zero) then
   vel(SDIM)=advbot
  else if (time.eq.zero) then
   vel(SDIM)=xblob4
  else
   print *,"time invalid in get pipe velocity"
   print *,"time= ",time
   print *,"dist= ",dist
   print *,"axis_dir ",axis_dir
   stop
  endif
 else
  print *,"dist(1) corrupt: ",dist(1)
  stop
 endif
else if (axis_dir.eq.3) then
 if (nmat.eq.3) then
  !do nothing
 else
  print *,"expecting nmat==3: ",nmat
  stop
 endif
 if (dist(nmat).ge.zero) then
  vel(SDIM)=zero
 else if (z.ge.EPS2*dx(SDIM)) then !y=z if 2D
  if (time.gt.zero) then
   vel(SDIM)=vinletgas
  else if (time.eq.zero) then
   vel(SDIM)=yblob4
  else
   print *,"time invalid in get pipe velocity"
   print *,"time= ",time
   print *,"dist= ",dist
   print *,"axis_dir ",axis_dir
   stop
  endif
 else if (z.le.EPS2*dx(SDIM)) then ! y=z if 2D
  if (time.gt.zero) then
   vel(SDIM)=advbot
  else if (time.eq.zero) then
   vel(SDIM)=xblob4
  else
   print *,"time invalid in get pipe velocity"
   print *,"time= ",time
   print *,"dist= ",dist
   print *,"axis_dir ",axis_dir
   stop
  endif
 else
  print *,"z corrupt: ",z
  stop
 endif

else if ((axis_dir.eq.0).or. &
         (axis_dir.eq.1).or. &
         (axis_dir.eq.2)) then

 if (nmat.eq.3) then
  !do nothing
 else
  print *,"expecting nmat==3: ",nmat
  stop
 endif
 if (dist(nmat).ge.zero) then
  vel(SDIM)=zero
 else if (dist(1).ge.zero) then
  if (time.gt.zero) then
   vel(SDIM)=vinletgas
  else if (time.eq.zero) then
   vel(SDIM)=yblob4
  else
   print *,"time invalid in get pipe velocity"
   print *,"time= ",time
   print *,"dist= ",dist
   print *,"axis_dir ",axis_dir
   stop
  endif
 else if (dist(1).le.zero) then
  if (time.gt.zero) then
   vel(SDIM)=advbot
  else if (time.eq.zero) then
   vel(SDIM)=xblob4
  else
   print *,"time invalid in get pipe velocity"
   print *,"time= ",time
   print *,"dist= ",dist
   print *,"axis_dir ",axis_dir
   stop
  endif
 else
  print *,"dist(1) corrupt: ",dist(1)
  stop
 endif

 ! above: axis_dir=0,1,2 (and above that 3,4)

else if (axis_dir.eq.5) then

 if (SDIM.eq.3) then
  r = sqrt(y*y+z*z)
  zlocal=z
 else if (SDIM.eq.2) then
  r = sqrt(y*y)
  zlocal=zero
 else
  print *,"dimension bust"
  stop
 endif

 if (r.gt.radblob) then
  x_vel=zero
 else 
  x_vel=adv_vel*r*(radblob-r)
 endif
 y_vel=zero
 z_vel=zero
 if((zlocal-0.3)**2+(y-0.3)**2+(x-0.5)**2.le.0.4)then
   x_vel = .3d0
   y_vel = .4d0
   z_vel = -.2d0
 else if((zlocal+0.3)**2+(y-1.7)**2+(x-5.5)**2.le.0.4)then
   x_vel = -.3d0
   y_vel = -.4d0
   z_vel = .2d0
 else if((zlocal-0.9)**2+(y)**2+(x-1.5)**2.le.1.)then
   x_vel = -x_vel
   y_vel = -y_vel
   z_vel = -z_vel 
 else if((zlocal+0.9)**2+(y)**2+(x-3.5)**2.le.1.)then
   x_vel = -y_vel
   y_vel = -x_vel
   z_vel = -z_vel 
 else if((zlocal)**2+(y)**2+(x-2.5)**2.le.0.1)then
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
 print *,"axis_dir invalid get_pipe_velocity axis_dir:",axis_dir
 stop
endif

return
end subroutine get_pipe_velocity


! dir=velocity component
subroutine MITSUHIRO_PIPE_CFL_HELPER(time,dir,uu,dx)
use probcommon_module
implicit none
integer, INTENT(in) :: dir
real(amrex_real), INTENT(in) :: time
real(amrex_real), INTENT(inout) :: uu
real(amrex_real), INTENT(in) :: dx(SDIM)

real(amrex_real) utest


if ((dir.lt.0).or.(dir.ge.SDIM)) then
 print *,"dir invalid: ",dir
 stop
endif

if (dir.eq.0) then
 ! do nothing
else if (dir.eq.1) then
 ! do nothing
else if ((dir.eq.2).and.(SDIM.eq.3)) then
 ! do nothing
else
 print *,"dir invalid MITSUHIRO_PIPE_CFL_HELPER: ",dir
 stop
endif

if ((probtype.eq.41).or.(probtype.eq.3)) then
 if (dir.eq.adv_dir-1) then
  utest=abs(adv_vel)
  uu=max(abs(uu),abs(utest))
  utest=abs(advbot)
  uu=max(abs(uu),abs(utest))
  utest=abs(vinletgas)
  uu=max(abs(uu),abs(utest))
 endif
else
 print *,"expecting probtype=41,3:",probtype
 stop
endif

return
end subroutine MITSUHIRO_PIPE_CFL_HELPER


  ! do any initial preparation needed
subroutine INIT_MITSUHIRO_PIPE_MODULE()
IMPLICIT NONE

  DEF_VAPOR_GAMMA =  1.666666667D0

return
end subroutine INIT_MITSUHIRO_PIPE_MODULE

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

if ((probtype.eq.41).or.(probtype.eq.3)) then

 call pipedist(x(1),x(2),x(SDIM),LS,nmat)

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

if (velsolid_flag.eq.1) then
 !do nothing
else if (velsolid_flag.eq.0) then

 call get_pipe_velocity(x(1),x(2),x(SDIM),dx,VEL,t,nmat)

else
 print *,"velsolid_flag in MITSUHIRO_PIPE_VEL invalid: ",velsolid_flag
 stop
endif

return 
end subroutine MITSUHIRO_PIPE_VEL


! this routine used if pressure boundary conditions are prescribed,
! since only top wall is "outflow" (outflow in quotes since ice shrinks when
! melting), and flow is incompressible, ok to make the top wall pressure zero.
subroutine MITSUHIRO_PIPE_PRES(x,t,LS,PRES,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: LS(nmat)
real(amrex_real), INTENT(out) :: PRES
real(amrex_real) :: pipexlo,pipexhi
integer :: gravity_dir

call fort_derive_gravity_dir(gravity_vector,gravity_dir)

if (num_materials.eq.nmat) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
PRES=zero

if ((probtype.eq.41).or.(probtype.eq.3)) then ! presBDRYCOND 2D yhi

 pipexlo=problox
 pipexhi=probhix
 if ((axis_dir.eq.1).or.(axis_dir.eq.2).or.(axis_dir.eq.3)) then
  pipexlo=zero
  pipexhi=two*radblob3
 endif

 if (axis_dir.eq.0) then
  ! do nothing
 else if (axis_dir.eq.1) then

  if (x(1).lt.pipexlo) then
   PRES=zero
  else if (x(1).lt.xblob) then
   PRES=fort_denconst(2)*abs(gravity_vector(gravity_dir))*(x(1)-pipexlo)
  else
   PRES= &
    fort_denconst(2)*abs(gravity_vector(gravity_dir))* &
    (xblob-pipexlo)+ &
    fort_denconst(1)*abs(gravity_vector(gravity_dir))*(x(1)-xblob) 
  endif

 else if (axis_dir.eq.2) then

  if (x(1).lt.pipexlo) then
   PRES=zero
  else if (x(1).lt.xblob-radblob2) then
   PRES=fort_denconst(1)*abs(gravity_vector(gravity_dir))*(x(1)-pipexlo)
  else if (x(1).lt.xblob+radblob2) then
   PRES=fort_denconst(1)*abs(gravity_vector(gravity_dir))*(xblob-radblob2-pipexlo)+ &
       fort_denconst(2)*abs(gravity_vector(gravity_dir))*(x(1)-xblob+radblob2) 
  else
   PRES=fort_denconst(1)*abs(gravity_vector(gravity_dir))*(xblob-radblob2-pipexlo)+ &
       fort_denconst(2)*abs(gravity_vector(gravity_dir))*(two*radblob2)+ & 
       fort_denconst(1)*abs(gravity_vector(gravity_dir))*(x(1)-xblob-radblob2) 
  endif

 else if (axis_dir.eq.3) then

  if (x(1).lt.pipexlo) then
   PRES=zero
  else 
   PRES=fort_denconst(1)*abs(gravity_vector(gravity_dir))*(x(1)-pipexlo)
  endif

 else if (axis_dir.eq.4) then
  ! do nothing 
 else
  print *,"axis_dir invalid for pipe problem"
  stop
 endif  

else
 print *,"probtype invalid for pipe problem"
 stop
endif  

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
if ((num_materials.ge.2).and. &
    (num_state_material.ge.2).and. & ! density, temperature, vapor spec
    ((probtype.eq.41).or. &
     (probtype.eq.3))) then
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

end module MITSUHIRO_PIPE_module
