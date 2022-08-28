#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif
 
#include "AMReX_FORT_INTEGER.H"
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

! probtype==915
! axis_dir=0 is run2d/inputs.wavychannel
! axis_dir=1 is run3d/inputs.viscoHelix


!probtype==915 axis_dir=1 run3d/inputs.viscoHelix
!-----------------------------------------------
!yblob2=#tube radius
!yblob3=pitch angle
!yblob4=w rotation
!yblob7=pulling velocity
!yblob5=R
module WAVY_Channel_module

implicit none                   

contains

subroutine INIT_WAVY_MODULE() 
use probcommon_module
IMPLICIT NONE

return
end subroutine INIT_WAVY_MODULE

REAL_T function DIST_FINITE_CYLHEAD(P,R_cyl,H_bot,H_top)
 ! Returns the signed distance function to the cylinder
 ! surfaces (including top and bottom)
 ! The axis of cylinder is along SDIM=2 direction
 ! Cylinder radus is R_cyl
 ! Bottom and top faces are at H_bot and H_top
 ! Inside the cylinder < 0
 ! Outside the cylinder > 0
 implicit none

 REAL_T, INTENT(in), dimension(SDIM) :: P
 REAL_T, INTENT(in) :: R_cyl
 REAL_T, INTENT(in) :: H_bot
 REAL_T, INTENT(in) :: H_top
 
 REAL_T x,y,z,r
 REAL_T dist_cyl, dist_end
 
 x=P(1)
 y=P(2)
 z=P(SDIM)

 r = sqrt(x**2+y**2)
 if((H_bot.le.z).and.(z.le.H_top)) then
  ! between top and bottom
  if(r.ge.R_cyl) then
   ! outside
   DIST_FINITE_CYLHEAD = r-R_cyl
  else if (r.lt.R_cyl) then
   ! inside
   dist_cyl = R_cyl-r
   dist_end = min(H_top-z,z-H_bot)
   DIST_FINITE_CYLHEAD = -min(dist_cyl,dist_end)
  else
   print *,"r=",r
   print *,"invalid r value at DIST_FINITE_CYLHEAD (1)"
   stop
  endif

 else if (z.gt.H_top) then
  ! higher than top
  if(r.le.R_cyl) then
   ! inside infinite cylinder
   DIST_FINITE_CYLHEAD = z-H_top
  else if (r.gt.R_cyl) then
   ! outside infinite cylinder
   ! distance to the edge of the top
   DIST_FINITE_CYLHEAD = &
    sqrt((r-R_cyl)**2 + (z-H_top)**2)
  else
   print *,"r=",r
   print *,"invalid r value at DIST_FINITE_CYLHEAD (2)"
   stop
  endif

 else if (z.lt.H_bot) then
  ! lower than bottom
  if(r.le.R_cyl) then
   ! inside infinite cylinder
   DIST_FINITE_CYLHEAD = H_bot-z
  else if (r.gt.R_cyl) then
   ! outside infinite cylinder
   ! distance to the edge of the bottom
   DIST_FINITE_CYLHEAD = &
    sqrt((r-R_cyl)**2 + (H_bot-z)**2)
  else
   print *,"r=",r
   print *,"invalid r value at DIST_FINITE_CYLHEAD (3)"
   stop
  endif
 else
  print *,"invalid z value at DIST_FINITE_CYLHEAD"
  stop
 endif
end function DIST_FINITE_CYLHEAD

subroutine WAVY_INIT_LS_core(x,t,LS,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(out) :: LS(nmat)
INTEGER_T im
REAL_T tt,ss,wvel
REAL_T thickness,Radius,LambdaWave,piin
REAL_T HelixLength, scalelength, cylinderHeight,LStmp1,LStmp2,BodyStart

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

ss=t
wvel=yblob7
piin=4.0d0*atan(1.0d0) 

if (radblob4.gt.zero) then
 scalelength=radblob4
 thickness=yblob2/radblob4
 Radius = yblob5/radblob4
 HelixLength=radblob2/radblob4
 LambdaWave=(Radius-thickness/2.0)/tan(yblob3)
 cylinderHeight=radblob3/radblob4
 BodyStart=radblob5/radblob4
else if (radblob4.eq.zero) then
 scalelength=zero
 thickness=zero
 Radius = zero
 HelixLength=zero
 LambdaWave=zero
 cylinderHeight=zero
 BodyStart=zero
else
 print *,"radblob4 invalid"
 stop
endif

if (axis_dir.eq.0) then

 ! 12 x 12 domain
 if (SDIM.eq.2) then
  if ((num_materials.eq.3).and.(probtype.eq.915)) then
   do im=1,num_materials 
    if (im.eq.1) then !vapor
     LS(im)=x(SDIM)-7.0
    else if (im.eq.2) then ! liquid
     LS(im)=7.0 - x(SDIM)
    else if (im.eq.3) then ! wavy substrate
     ! does not have to be an exact distance.  Main thing
     ! is that the function is smooth and the zero LS is correct.
     ! also, since RZ, use cos(x(1) pi/3) since this is symmetric
     ! at r=0.
     LS(im)=three+sin(two*Pi*t)+cos(x(1)*Pi/three)-x(SDIM)
    else
     print *,"im invalid"
     stop
    endif
   enddo ! im=1..num_materials
  else
   print *,"num_materials or probtype invalid"
   stop
  endif
 else
  print *,"not setup for this dimension yet dim=",SDIM
  stop
 endif
else if (axis_dir.eq.1) then

 if (SDIM.eq.3) then
  if ((num_materials.eq.3).and.(probtype.eq.915)) then
   do im=1,num_materials 
    if (im.eq.1) then !water
     LS(im)=-1000.0
    elseif(im.eq.2) then
     LS(im)=1000.0
    elseif (im.eq.3) then ! solid helix
     call GET_ROOT(x(1),x(2),x(SDIM), &
      tt,ss,wvel,yblob5,yblob3,yblob4, &
      Radius-thickness/2.0d0, &
      LambdaWave,helixlength,scalelength,BodyStart)
     LS(im)=-(sqrt((x(2)-yblob7)**2+(x(SDIM)-yblob8)**2)-radblob7) !circle or cylinder
!    LS(im) = 0.5- sqrt( (cos(tt+t)-x(1))**2 + (sin(tt+t)-x(2))**2 + (tt-x(3))**2 )  !without z velocity
! print*, 'tt=',tt
!print*, 'Radius*cos(tt+yblob4*t)',Radius*cos(tt+yblob4*t)
!print*, 'Radius*sin(tt+yblob4*t)',Radius*sin(tt+yblob4*t)
!print*, 'LambdaWave*tt',LambdaWave*tt
     LStmp1  = thickness/2.0d0-  &
          sqrt( ((Radius-thickness/2.0d0)*cos(tt+yblob4*t)-x(1))**2 &
              + ((Radius-thickness/2.0d0)*sin(tt+yblob4*t)-x(2))**2 &
              + (LambdaWave*tt+wvel*t+BodyStart-x(SDIM))**2 ) 
     LStmp2  =-DIST_FINITE_CYLHEAD(x(1:SDIM),Radius, &
      HelixLength+BodyStart+wvel*t, &
      HelixLength+cylinderHeight+BodyStart+wvel*t)   !minus or plus
     LS(im) =max(LStmp1,LStmp2)  !check later
    else
     print *,"im invalid, im=",im
     stop
    endif
   enddo ! im=1..num_materials
  else
   print *,"num_materials or probtype invalid",num_materials,probtype
   stop
  endif
 else
  print *,"not setup for this dimension yet dim=",SDIM
  stop
 endif
else
 print *,"axis_dir invalid in WAVY_INIT_LS_core"
 stop
endif

return
end subroutine WAVY_INIT_LS_core

subroutine WAVY_INIT_LS(x,t,LS,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(out) :: LS(nmat)
REAL_T :: x3D(3)
REAL_T :: x3D_foot(3)
INTEGER_T :: dir
INTEGER_T im
INTEGER_T auxcomp

if (t.ge.0.0d0) then
 ! do nothing
else
 print *,"t invalid"
 stop
endif

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

do dir=1,SDIM
 x3D(dir)=x(dir)
enddo
x3D(3)=x(SDIM)

if (axis_dir.eq.0) then
 call WAVY_INIT_LS_core(x,t,LS,nmat)
else if (axis_dir.eq.1) then

 if (SDIM.eq.3) then
  if ((num_materials.eq.3).and.(probtype.eq.915)) then
   do im=1,num_materials 
    if (im.eq.1) then !water
     LS(im)=-1000.0
    elseif(im.eq.2) then
     LS(im)=1000.0
    elseif (im.eq.3) then ! solid helix
     auxcomp=1
     x3D_foot(1)=x3D(1)*cos(yblob4*t)+x3D(2)*sin(yblob4*t)
     x3D_foot(2)=x3D(2)*cos(yblob4*t)-x3D(1)*sin(yblob4*t)
     x3D_foot(3)=x3D(3)
     call interp_from_aux_grid(auxcomp,x3D_foot,LS(im))
    else
     print *,"im invalid, im=",im
     stop
    endif
   enddo ! im=1..num_materials
  else
   print *,"num_materials or probtype invalid",num_materials,probtype
   stop
  endif
 else
  print *,"not setup for this dimension yet dim=",SDIM
  stop
 endif
else
 print *,"axis_dir invalid in WAVY_INIT_LS"
 stop
endif

return
end subroutine WAVY_INIT_LS

subroutine WAVY_BOUNDING_BOX_AUX(auxcomp, &
    minnode,maxnode,LS_FROM_SUBROUTINE,aux_ncells_max_side)
use probcommon_module
use global_utility_module
IMPLICIT NONE
INTEGER_T, INTENT(in) :: auxcomp
REAL_T, INTENT(inout) :: minnode(3)
REAL_T, INTENT(inout) :: maxnode(3)
INTEGER_T, INTENT(out) :: LS_FROM_SUBROUTINE
INTEGER_T, INTENT(out) :: aux_ncells_max_side

 if (auxcomp.eq.1) then
  if (axis_dir.eq.1) then
   if (SDIM.eq.3) then
    if (num_materials.eq.3) then
     LS_FROM_SUBROUTINE=1
     aux_ncells_max_side=100
     minnode(1)=-0.15
     maxnode(1)=0.15
     minnode(2)=-0.15
     maxnode(2)=0.15
     minnode(3)=3.925
     maxnode(3)=6.875
    else
     print *,"num_materials invalid in WAVY_BOUNDING_BOX_AUX"
     stop
    endif
   else
    print *,"sdim invalid in WAVY_BOUNDING_BOX_AUX"
    stop
   endif
  else
   print *,"axis_dir invalid in WAVY_BOUNDING_BOX_AUX"
   stop
  endif

 else
  print *,"auxcomp invalid in WAVY_BOUNDING_BOX_AUX"
  stop
 endif


end subroutine WAVY_BOUNDING_BOX_AUX

subroutine WAVY_AUX_DATA(auxcomp,x,LS)
use probcommon_module
use global_utility_module
IMPLICIT NONE
INTEGER_T, INTENT(in) :: auxcomp
REAL_T, INTENT(in) :: x(3)
REAL_T, INTENT(out) :: LS
REAL_T :: local_LS(num_materials)
REAL_T :: local_time

 if (auxcomp.eq.1) then
  if (axis_dir.eq.1) then
   if (SDIM.eq.3) then
    if (num_materials.eq.3) then
     local_time=0.0d0
     call WAVY_INIT_LS_core(x,local_time,local_LS,num_materials)
     LS=local_LS(num_materials)
    else
     print *,"num_materials invalid in WAVY_AUX_DATA"
     stop
    endif
   else
    print *,"sdim invalid in WAVY_AUX_DATA"
    stop
   endif
  else
   print *,"axis_dir invalid in WAVY_AUX_DATA"
   stop
  endif

 else
  print *,"auxcomp invalid in WAVY_AUX_DATA"
  stop
 endif

end subroutine WAVY_AUX_DATA

subroutine WAVY_OVERRIDE_TAGFLAG(xsten,nhalf,time,rflag,tagflag)
use probcommon_module
use global_utility_module
IMPLICIT NONE
INTEGER_T, INTENT(in) :: nhalf
REAL_T, INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
REAL_T, INTENT(in) :: time
REAL_T, INTENT(inout) :: rflag
INTEGER_T, INTENT(inout) :: tagflag
REAL_T, dimension(3) :: local_x
REAL_T, dimension(SDIM) :: local_delta
INTEGER_T :: dir
INTEGER_T :: auxcomp
REAL_T :: LS

if (nhalf.lt.1) then
 print *,"nhalf invalid wavy override tagflag"
 stop
endif
do dir=1,SDIM
 local_x(dir)=xsten(0,dir)
enddo
local_x(3)=xsten(0,SDIM)
auxcomp=1
call WAVY_AUX_DATA(auxcomp,local_x,LS)
do dir=1,SDIM
 local_delta(dir)=xsten(1,dir)-xsten(-1,dir)
 if (local_delta(dir).gt.zero) then
  ! do nothing
 else
  print *,"local_delta invalid wavy_override_tagflag"
  stop
 endif
enddo !dir=1..sdim
if (abs(LS).le.local_delta(1)) then 
 rflag=1.0d0
 tagflag=1
else if (abs(LS).gt.local_delta(1)) then
 ! do nothing
else
 print *,"bust in wavy_override_tagflag"
 stop
endif 
 
end subroutine WAVY_OVERRIDE_TAGFLAG

subroutine WAVY_INIT_VEL(x,t,LS,VEL,velsolid_flag,dx,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: dx(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(nmat)
REAL_T, INTENT(out) :: VEL(SDIM)
INTEGER_T dir
INTEGER_T, INTENT(in) :: velsolid_flag
INTEGER_T im_solid_Tomas
INTEGER_T im_solid_wavy

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

if (axis_dir.eq.0) then ! wavy_channel
 do dir=1,SDIM
  VEL(dir)=zero
 enddo
 im_solid_wavy=3
 if ((LS(im_solid_wavy).ge.zero).or. &
     (velsolid_flag.eq.1)) then
  VEL(SDIM)=two*Pi*cos(two*Pi*t)
 else if ((LS(im_solid_wavy).le.zero).and. &
          (velsolid_flag.eq.0)) then
  ! do nothing
 else
  print *,"LS(im_solid_wavy) bust"
  stop
 endif
else if (axis_dir.eq.1) then ! Tomas Solano

 if (SDIM.eq.3) then

  do dir=1,SDIM
   VEL(dir)=zero
  enddo
  im_solid_Tomas=3
  if ((LS(im_solid_Tomas).ge.zero).or. &
      (velsolid_flag.eq.1)) then
   if (1.eq.0) then
    VEL(2)=-x(SDIM)
    VEL(SDIM)=x(2)
    VEL(1)=0.0d0
   else
    VEL(1)=-x(2)*yblob4
    VEL(2)=x(1)*yblob4
    VEL(SDIM)=yblob7
   endif

  else if ((LS(im_solid_Tomas).lt.zero).and. &
           (velsolid_flag.eq.0)) then

   if (1.eq.1) then
    VEL(adv_dir) = adv_vel
   else
    VEL(adv_dir) = zero
   endif

  else
   print *,"LS(im_solid_tomas) bust"
   stop
  endif
 else
  print *,"expecting 3D if axis_dir.eq.1"
  stop
 endif

else
 print *,"axis_dir invalid"
 stop
endif

return 
end subroutine WAVY_INIT_VEL

subroutine WAVY_INIT_PRES(x,t,LS,PRES,nmat)
use probcommon_module
use global_utility_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(nmat)
REAL_T, INTENT(out) :: PRES
REAL_T :: gravity_dz
INTEGER_T :: gravity_dir

call fort_derive_gravity_dir(gravity_vector,gravity_dir)

if (num_materials.eq.nmat) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

if (gravity_dir.eq.1) then
 gravity_dz=x(gravity_dir)-probhix
else if (gravity_dir.eq.2) then
 gravity_dz=x(gravity_dir)-probhiy
else if ((gravity_dir.eq.SDIM).and.(SDIM.eq.3)) then
 gravity_dz=x(gravity_dir)-probhiz
else
 print *,"dimension bust"
 stop
endif

if (axis_dir.eq.0) then ! wavy channel
 PRES=-fort_denconst(1)*abs(gravity_vector(gravity_dir))*gravity_dz
else if (axis_dir.eq.1) then ! Tomas
 PRES=zero
else
 print *,"axis_dir invalid"
 stop
endif

return 
end subroutine WAVY_INIT_PRES

subroutine WAVY_INIT_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: bcflag !0=called from initialize  1=called from bc
INTEGER_T, INTENT(in) :: nmat
INTEGER_T, INTENT(in) :: nstate_mat
REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(nmat)
REAL_T, INTENT(out) :: STATE(nmat*nstate_mat)
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

if (axis_dir.eq.0) then ! wavy channel

 if ((num_materials.eq.3).and. &
     (num_state_material.eq.2).and. &
     (num_species_var.eq.0).and. &
     (probtype.eq.915)) then
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
  print *,"num_materials,num_state_material, or probtype invalid"
  stop
 endif

else if (axis_dir.eq.1) then ! Tomas
 
 if ((num_materials.eq.3).and. &
     (num_state_material.eq.2).and. &
     (num_species_var.eq.0).and. &
     (probtype.eq.915)) then
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
  enddo ! im=1..num_materials
 else
  print *,"num_materials,num_state_material,num_species_var,",num_materials,num_state_material,num_species_var
  print *,"or probtype invalid"
  stop
 endif
else
 print *,"axis_dir invalid in WAVY_INIT_STATE"
 stop
endif

return
end subroutine WAVY_INIT_STATE

 ! dir=1..sdim  side=1..2
subroutine WAVY_LS_BC(xwall,xghost,t,LS, &
   LS_in,dir,side,dx,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
REAL_T, INTENT(in) :: xwall
REAL_T, INTENT(in) :: xghost(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(inout) :: LS(nmat)
REAL_T, INTENT(in) :: LS_in(nmat)
INTEGER_T, INTENT(in) :: dir,side
REAL_T, INTENT(in) :: dx(SDIM)

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif
if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then
 call WAVY_INIT_LS(xghost,t,LS,nmat)
else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine WAVY_LS_BC

 ! dir=1..sdim  side=1..2 veldir=1..sdim
subroutine WAVY_VEL_BC(xwall,xghost,t,LS, &
   VEL,VEL_in,veldir,dir,side,dx,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
REAL_T, INTENT(in) :: xwall
REAL_T, INTENT(in) :: xghost(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(nmat)
REAL_T, INTENT(inout) :: VEL
REAL_T, INTENT(in) :: VEL_in
INTEGER_T, INTENT(in) :: veldir,dir,side
REAL_T, INTENT(in) :: dx(SDIM)
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

 call WAVY_INIT_VEL(xghost,t,LS,local_VEL,velsolid_flag,dx,nmat)
 VEL=local_VEL(veldir)

else
 print *,"dir,side, or veldir invalid"
 stop
endif

return
end subroutine WAVY_VEL_BC

 ! dir=1..sdim  side=1..2
subroutine WAVY_PRES_BC(xwall,xghost,t,LS, &
   PRES,PRES_in,dir,side,dx,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
REAL_T, INTENT(in) :: xwall
REAL_T, INTENT(in) :: xghost(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(nmat)
REAL_T, INTENT(inout) :: PRES
REAL_T, INTENT(in) :: PRES_in
INTEGER_T, INTENT(in) :: dir,side
REAL_T, INTENT(in) :: dx(SDIM)

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then

 call WAVY_INIT_PRES(xghost,t,LS,PRES,nmat)

else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine WAVY_PRES_BC

 ! dir=1..sdim  side=1..2
subroutine WAVY_STATE_BC(xwall,xghost,t,LS, &
   STATE,STATE_merge,STATE_in,im,istate,dir,side,dx,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
REAL_T, INTENT(in) :: xwall
REAL_T, INTENT(in) :: xghost(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(nmat)
REAL_T local_STATE(nmat*num_state_material)
REAL_T, INTENT(inout) :: STATE
REAL_T, INTENT(inout) :: STATE_merge
REAL_T, INTENT(in) :: STATE_in
INTEGER_T, INTENT(in) :: dir,side
REAL_T, INTENT(in) :: dx(SDIM)
INTEGER_T, INTENT(in) :: istate,im
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
 call WAVY_INIT_STATE(xghost,t,LS,local_STATE, &
         local_bcflag,nmat,num_state_material)
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
end subroutine WAVY_STATE_BC


subroutine WAVY_HEATSOURCE( &
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

INTEGER_T, INTENT(in) :: nmat
INTEGER_T, INTENT(in) :: im
REAL_T, INTENT(in) :: VFRAC(nmat)
REAL_T, INTENT(in) :: time
INTEGER_T, INTENT(in) :: nhalf
REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
REAL_T, INTENT(in) :: temp(nmat)
REAL_T, INTENT(in) :: den(nmat)
REAL_T, INTENT(in) :: CV(nmat)
REAL_T, INTENT(in) :: dt
REAL_T, INTENT(out) :: heat_source

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

if ((num_materials.eq.3).and.(probtype.eq.915)) then
 heat_source=zero
else
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine WAVY_HEATSOURCE


end module WAVY_Channel_module


! Written by: Tomas Solano
! Written date : 01/07/2020

! This program finds the root of a function using bisection or interpolation

! This algorithm was originated by T. J. Dekker.  An Algol 60 version,
! with some improvements as given in MATLAB's function fzero, is given 
! by R. P. Brent in "Algorithms for Minimization Without Derivatives", 
! Prentice-Hall, 1973. Also in Forsythe, Malcolm and Moler, "Computer Methods
! for Mathematical Computations", Prentice-Hall, 1976.

!call GET_ROOT(x(1),x(2),x(3),tt,ss,wvel,yblob5,yblob3,yblob4)

SUBROUTINE GET_ROOT(px,py,pz,z,ss,wvel,yblob5,yblob3,yblob4,Radius,LambdaWave,helixlength,scalelength,BodyStart)
  
  REAL(KIND=8):: z,ss,wvel,yblob5,yblob3,yblob4
  REAL(KIND=8) :: Radius,LambdaWave,helixlength,scalelength,BodyStart
INTERFACE
  FUNCTION F_Tomas(x,px,py,pz,ss,wvel,yblob5,yblob3,yblob4,Radius,LambdaWave,helixlength,scalelength) RESULT(g)
    IMPLICIT NONE
    REAL(KIND=8),INTENT(IN):: x,px,py,pz,ss,wvel,yblob5,yblob3,yblob4,Radius,LambdaWave,helixlength,scalelength
    REAL(KIND=8):: g
  END Function F_Tomas
END INTERFACE    
  
  INTEGER:: errCode
  INTEGER,PARAMETER:: MAXITER=25
  INTEGER:: neval   ! not used
  REAL(KIND=8):: xZero,fZero
  REAL(KIND=8),PARAMETER:: a = 0.0
!   REAL(KIND=8),PARAMETER:: b = 10.0
  REAL(KIND=8),PARAMETER:: TOL = 1E-10
  REAL(KIND=8), INTENT(IN) :: px,py,pz
  
!   CALL BrentZeroDouble(a,b,F,px,py,pz,tol,MAXITER,neval,errCode,xZero,fZero)
  CALL BrentZeroDouble(a,F_Tomas,px,py,pz,ss,wvel,yblob5,yblob3,yblob4,Radius,LambdaWave,helixlength,scalelength,BodyStart,tol,MAXITER,neval,errCode,xZero,fZero)

!BrentZeroDouble(x0,F_Tomas,px,py,pz,ss,wvel,yblob5,yblob3,yblob4,tol,maxIter,neval,errCode,xZero,fZero)


  z=xZero
  
END SUBROUTINE GET_ROOT
!======================================================================================!

!======================================================================================!
FUNCTION F_Tomas(x,px,py,pz,ss,wvel,yblob5,yblob3,yblob4,Radius,LambdaWave,helixlength,scalelength) RESULT (FReturn)        !Function of which to find zero

  REAL(KIND=8), INTENT(IN) :: x,px,py,pz,ss,wvel,yblob5,yblob3,yblob4,Radius,LambdaWave,helixlength,scalelength
  REAL(KIND=8) :: FReturn

!  FReturn=2*px*sin(x+ss)-2*py*cos(x+ss)+2*x-2*pz
  FReturn=2*Radius*px*sin(x+yblob4*ss)-2*Radius*py*cos(x+yblob4*ss)+2**LambdaWave**2*x+2*LambdaWave*wvel*ss - 2*LambdaWave*pz
  RETURN
END Function F_Tomas   
!======================================================================================!
!        LS(im) = yblob2/2.0/radblob4- sqrt( Radius*cos(tt+yblob4*t)-x(1))**2 + (Radius*sin(tt+yblob4*t)-x(2))**2 + (LambdaWave*tt+wvel*t-x(3))**2 ) 
!======================================================================================!
! SUBROUTINE
! BrentZeroDouble(ax,bx,F,px,py,pz,tol,maxIter,neval,errCode,xZero,fZero)
SUBROUTINE BrentZeroDouble(x0,F_Tomas,px,py,pz,ss,wvel,yblob5,yblob3,yblob4,Radius,LambdaWave,helixlength,scalelength,BodyStart,tol,maxIter,neval,errCode,xZero,fZero)
! PURPOSE - Compute a zero of F in the interval (ax,bx)

  ! REAL(KIND=8),INTENT(IN):: ax,bx   ! left and right enKIND=8oints of interval
  REAL(KIND=8), INTENT(IN) :: x0,px,py,pz,ss,wvel,yblob5,yblob3,yblob4,Radius,LambdaWave,helixlength,scalelength,BodyStart !initial guess, point in question,time
  REAL(KIND=8),INTENT(IN):: tol     ! desired interval of uncertainity 
  INTEGER,INTENT(IN):: maxIter   ! max number of iterations allowed. 25 is good
  INTEGER,INTENT(OUT):: neval
  INTEGER,INTENT(OUT):: errCode   ! =0 is OK; =1 too many iterations
                                  ! =2 if F(ax) and F(bx) have the same sign
  REAL(KIND=8),INTENT(OUT):: xZero,fZero ! the last and best value of the zero                                   
      
INTERFACE
  FUNCTION F_Tomas(x,px,py,pz,ss,wvel,yblob5,yblob3,yblob4,Radius,LambdaWave,helixlength,scalelength) RESULT(g)
    IMPLICIT NONE
    REAL(KIND=8),INTENT(IN):: x,px,py,pz,ss,wvel,yblob5,yblob3,yblob4,Radius,LambdaWave,helixlength,scalelength
    REAL(KIND=8):: g
  END Function F_Tomas
END INTERFACE    
  
  REAL(KIND=8):: a,b,c,d,e,eps
  REAL(KIND=8):: fa,fb,fc,tol1
  INTEGER:: kIter,i,i_min,itest,itestp,ntestpmax
!  INTEGER:: method   ! =0 bisection; =1 linear; =2 inverse quadratic
  REAL(KIND=8):: xm,p,q,r,s,dist,min_dist,Ltotal,s1start,s1end,s_start,s_end
  REAL(KIND=8),PARAMETER:: ZERO=0.0, ONE=1.0, TWO=2.0, THREE=3.0, HALF=0.5, twosqrt = sqrt(2.)
!----------------------------------------------------------------------------
  eps=EPSILON(x0)
  tol1=ONE+eps 
  !------------------------------------------------------------------------------- 
  !The following section was added from MATLAB's fzero function to find better
  !init point
!  a=pz
   a=0.5
   min_dist=1e9
   
!======================================================================================!
!        LS(im) = yblob2/2.0/radblob4- sqrt( Radius*cos(tt+yblob4*t)-x(1))**2 + (Radius*sin(tt+yblob4*t)-x(2))**2 + (LambdaWave*tt+wvel*t-x(3))**2 ) 
!======================================================================================!
! LambdaWave=(Radius-thickness/2.0)/tan(yblob3)

   Ltotal=helixlength/LambdaWave
   DO i=-0,200
     s = (dble(i)/200.0)*Ltotal
!	 s=s1*cos(yblob3)
     dist = (Radius*cos(s+yblob4*ss)-px)**2 + (Radius*sin(s+yblob4*ss)-py)**2 + (LambdaWave*s+wvel*ss+BodyStart-pz)**2
     IF(dist.lt.min_dist) THEN
       min_dist=dist
       a = s
       i_min=i
     ENDIF
   ENDDO

   DO i=max(i_min-3,0)*10,min(i_min+3,200)*10
     s = (dble(i)/2000.0)*Ltotal
!	 s=s1*cos(yblob3)
     dist = (Radius*cos(s+yblob4*ss)-px)**2 + (Radius*sin(s+yblob4*ss)-py)**2 + (LambdaWave*s+wvel*ss+BodyStart-pz)**2
     IF(dist.lt.min_dist) THEN
       min_dist=dist
       a = s
       i_min=i
     ENDIF
   ENDDO


   b= a
   fb=min_dist
   
   xZero=b
   fZero=fb
   errCode=0   ! SUCCESS! The proper way to leave
   RETURN
   
   
   if (i_min >75) then
        a= (70.0/80.0)*Ltotal
		b= (80.0/80.0)*Ltotal
		fa=	F_Tomas(a,px,py,pz,ss,wvel,yblob5,yblob3,yblob4,Radius,LambdaWave,helixlength,scalelength)
        fb=	F_Tomas(b,px,py,pz,ss,wvel,yblob5,yblob3,yblob4,Radius,LambdaWave,helixlength,scalelength)
   else if (i_min< 5) then
        a= (0)*Ltotal
		b= (10.0/80.0)*Ltotal
		fa=	F_Tomas(a,px,py,pz,ss,wvel,yblob5,yblob3,yblob4,Radius,LambdaWave,helixlength,scalelength)
        fb=	F_Tomas(b,px,py,pz,ss,wvel,yblob5,yblob3,yblob4,Radius,LambdaWave,helixlength,scalelength)
	else
		
   itest=0
   itestp=0
   ntestpmax=20;
   print*, i_min
   Do While (itest .eq. 0)
   s1start=(max(0.0,dble(i_min-1-itestp))/40.0)*Ltotal
   s1end=(min(40.0,dble(i_min+1+itestp))/40.0)*Ltotal 
	s_start=s1start !*cos(yblob3)   
	s_end=s1end  !*cos(yblob3)   

   fa=	F_Tomas(s_start,px,py,pz,ss,wvel,yblob5,yblob3,yblob4,Radius,LambdaWave,helixlength,scalelength)
   a=s_start
   fb=	F_Tomas(s_end,px,py,pz,ss,wvel,yblob5,yblob3,yblob4,Radius,LambdaWave,helixlength,scalelength)
   b=s_end   
   neval=2   
!        fa = F_Tomas(a,px,py,pz,ss,wvel,yblob5,yblob3,yblob4)
!        neval=neval+1   
        IF (abs(fa) .lt. eps)THEN
          xZero = a
          fZero = fa
          RETURN
        ENDIF
        IF (abs(fb) .lt. eps)THEN
          xZero = b
          fZero = fb
          RETURN
        ENDIF		
        IF(fa*fb > ZERO) THEN
		  itestp=itestp+1
		  IF(itestp .gt.  ntestpmax) then
		     stop "Error in bracketing"
		   Endif
		ELSE
		  itest=1
          Exit	
        Endif
  ENDDO		
  endif
		     
  c=a
  fc = fa
  d=b-a
  e=d
  !-----------------------------------------------------------------------------  
  DO kIter=1,maxIter
    IF ( (fb>0 .AND. fc>0) .OR. (fb<0 .AND. fc<0) ) THEN
      c=a  ! we insist that b and c straddle the zero
      fc=fa
      d=b-a
      e=d
    END IF

    IF (ABS(fc) < ABS(fb)) THEN
      a=b    ! we insist that b be the better guess of b and c
      b=c
      c=a
      fa=fb
      fb=fc
      fc=fa
    END IF

    tol1=TWO*eps*ABS(b)+HALF*tol   ! convergence test
    xm=HALF*(c-b)
    IF (ABS(xm) <= tol1 .OR. abs(fb) .lt. eps) THEN
      xZero=b
      fZero=fb
      errCode=0   ! SUCCESS! The proper way to leave
      RETURN
    END IF


    IF (ABS(e) < tol1 .OR. ABS(fa) <= ABS(fb) ) THEN
      d=xm   ! bisection
      e=d
!      method=0
    ELSE
      IF (abs(a-c).lt. eps) THEN
        s=fb/fa   ! linear interpolation
        p=TWO*xm*s
        q=ONE-s
!        method=1
      ELSE
        q=fa/fc   ! inverse quadratic interpolation
        r=fb/fc
        s=fb/fa
        p=s*(TWO*xm*q*(q-r)-(b-a)*(r-ONE))
        q=(q-ONE)*(r-ONE)*(s-ONE)
!        method=2
      END IF
      IF (p > ZERO) q=-q   ! adjust signs
      p=ABS(p)
      IF (p+p >= (THREE*xm*q-ABS(tol1*q)) .OR. p+p >= ABS(e*q) ) THEN
        d=xm   ! don't interpolate. Use bisection
        e=d
 !       method=-1
      ELSE
        e=d   ! OK, use interpolation
        d=p/q
      END IF
    END IF  
   
    a=b   ! complete step. a becomes the previous iteration
    fa=fb
   IF (ABS(d) > tol1) THEN
     b=b+d
   ELSE  
     b=b+SIGN(tol1,xm)
   END IF
     
   fb=F_Tomas(b,px,py,pz,ss,wvel,yblob5,yblob3,yblob4,Radius,LambdaWave,helixlength,scalelength)   ! the newest and best value (we hope)
   neval=neval+1   ! keep count of the function evaluations
  END DO
  
! The loop should never terminate. If it does, return the last iteration
!  and set errCode to 1
!  xZero=b
  
  xZero=b
  fZero=fb
  errCode=1
  stop "Error in finding minimum"
  RETURN  

END Subroutine BrentZeroDouble   
