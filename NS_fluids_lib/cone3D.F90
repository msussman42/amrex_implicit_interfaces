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

!!$ in_vel = (xblob10,yblob10,zblob10)
!!$ shear for input velocity
!!$ if t<xblob9
!!$  if yblob8 == 0.0 
!!$     u = u * (1 + yblob9 * y * exp(-t/radblob9))
!!$  else
!!$   if -yblob8 < y < yblob8
!!$     u = u * (1 + yblob9 * exp(-t/radblob9) * sin(pi*y/yblob8) * 
!!$             sin(radblob8*t))
!!$                                  
!!$ CONE => ns.axis_dir=0                                  
!!$ (before: in_vel=adv_vel) 
!!$                                  
!!$                  /
!!$  in_vel         /
!!$ ---->          /    SOLID
!!$               /
!!$         APEX / __________BISECT__\
!!$              \  |                /
!!$               \ / THETA_C
!!$                \
!!$                 \
!!$         GAS      \
!!$                                  
!!$                                 
!!$ To set the cone position/orientation change                      
!!$ the APEX(xyzblob) BISECT(xyzblob2) vectors, and THETA_C(radblob) angle
!!$                                 
!!$                                 
!!$
!!$ Cylinder => ns.axis_dir=1                                
!!$                    
!!$   in_vel          
!!$ ---->          _----_   
!!$               /      \
!!$              / SOLID  \
!!$              \        /
!!$               \_    _/   
!!$                 ----
!!$                  
!!$         GAS       
 
!!$ To set the cylinder position/orientation change                      
!!$ MEHDI_CENTER(xyzblob3), MEHDI_AXIS(xyzblob4), and MEHDI_RADIUS(radblob3)
!!$                                 
 
!!$ Sphere => ns.axis_dir=2
!!$                    
!!$   in_vel          
!!$ ---->          _----_   
!!$               /      \
!!$              / SOLID  \
!!$              \        /
!!$               \_    _/   
!!$                 ----
!!$                  
!!$         GAS       
!!$                                  
!!$ To set the sphere position/orientation change                      
!!$ the MEHDI_CENTER(xyzblob3), and MEHDI_RADIUS(radblob3)
 
#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else
 print *,"dimension bust"
 stop
#endif

 ! probtype==222
 module CONE3D_module

  implicit none
  REAL_T :: MEHDI_IN_VEL(SDIM) ! Inflow velocity  
  REAL_T :: APEX(SDIM)   ! Apex of the cone
  REAL_T :: BISECT(SDIM) ! Vector goes from apex toward base and 
                         ! bisect the cone angle
  REAL_T :: THETA_C      ! Half of the cone angle in radian 
  ! Cylinder and sphere
  REAL_T :: MEHDI_CENTER(SDIM)  ! A point on axis of cylinder
                                ! or center of sphere 
  REAL_T :: MEHDI_AXIS(SDIM)    ! Cylinder axis vector
  REAL_T :: MEHDI_RADIUS        ! Cylinder or sphere radius

 contains

  ! initial preparations
  subroutine INIT_CONE3D_MODULE()
   use probcommon_module
   IMPLICIT NONE

   INTEGER_T :: MEHDI_DIR
   REAL_T :: dotprod

   MEHDI_IN_VEL(1) = xblob10
   MEHDI_IN_VEL(2) = yblob10
   if (SDIM.eq.3) then
    MEHDI_IN_VEL(SDIM) = zblob10
   endif

   ! Cone
   if (axis_dir.eq.0) then
   
    if (adv_vel.ge.zero) then
     do MEHDI_DIR=1,SDIM
      MEHDI_IN_VEL(MEHDI_DIR)=zero
     enddo
     MEHDI_IN_VEL(1)=adv_vel
    else
     print *,"adv_vel invalid"
     stop
    endif

    ! APEX => xblob,yblob,zblob
    APEX(1) = xblob
    APEX(2) = yblob
    if (SDIM.eq.3) then
     APEX(SDIM) = zblob
    endif

    ! BISECT => xblob2,yblob2,zblob2
    BISECT(1) = xblob2
    BISECT(2) = yblob2
    if (SDIM.eq.3) then
     BISECT(SDIM) = zblob2
    endif
    dotprod=zero
    do MEHDI_DIR=1,SDIM
     dotprod=dotprod+BISECT(MEHDI_DIR)**2
    enddo
    dotprod=sqrt(dotprod)
    if (dotprod.gt.zero) then
     do MEHDI_DIR=1,SDIM
      BISECT(MEHDI_DIR) = BISECT(MEHDI_DIR) / dotprod 
     enddo
    else
     print *,"dotprod invalid"
     stop
    endif

    ! THETA_C => radblob
    THETA_C = radblob

    ! Cylinder
   else if (axis_dir.eq.1) then
    if(SDIM.eq.3)then
     MEHDI_CENTER(1) = xblob3
     MEHDI_CENTER(2) = yblob3
     MEHDI_CENTER(SDIM) = zblob3
     
     MEHDI_RADIUS = radblob3

     MEHDI_AXIS(1) = xblob4
     MEHDI_AXIS(2) = yblob4
     MEHDI_AXIS(SDIM) = zblob4

     dotprod=zero
     do MEHDI_DIR=1,SDIM
      dotprod=dotprod+MEHDI_AXIS(MEHDI_DIR)**2
     enddo
     dotprod=sqrt(dotprod)
     if (dotprod.gt.zero) then
      do MEHDI_DIR=1,SDIM
       MEHDI_AXIS(MEHDI_DIR) = MEHDI_AXIS(MEHDI_DIR) / dotprod 
      enddo
     else
      print *,"dotprod invalid"
      stop
     endif
    else
     print *,"cylinder(axis_dir=1) is only avaialble in 3D!"
     stop
    end if
    
    ! Sphere
   else if (axis_dir.eq.2) then
    MEHDI_CENTER(1) = xblob3
    MEHDI_CENTER(2) = yblob3
    if(SDIM.eq.3)then
     MEHDI_CENTER(SDIM) = zblob3
    endif

    MEHDI_RADIUS = radblob3
   else
    print *,"axis_dir invalid in cone3D.F90"
    stop
   end if
   if(1==1) then
    print *,"MEHDI_RADIUS", MEHDI_RADIUS
    do MEHDI_DIR = 1,SDIM
     print *,"MEHDI_AXIS", MEHDI_AXIS(MEHDI_DIR)
     print *,"MEHDI_CENTER", MEHDI_CENTER(MEHDI_DIR)
    enddo
    print *,"axis_dir", axis_dir
   endif


   return
  end subroutine INIT_CONE3D_MODULE

  ! fluids tessellate the domain, solids are immersed. 
  subroutine CONE3D_LS(x,t,LS,nmat)
  use probcommon_module
  IMPLICIT NONE

  INTEGER_T, INTENT(in) :: nmat
  REAL_T, INTENT(in) :: x(SDIM)
  REAL_T, INTENT(in) :: t
  REAL_T, INTENT(out) :: LS(nmat)
  INTEGER_T im

  if (nmat.eq.num_materials) then
   ! do nothing
  else
   print *,"nmat invalid"
   stop
  endif

  if ((num_materials.eq.2).and.(probtype.eq.222)) then
    do im=1,num_materials
     if (im.eq.1) then !air
      LS(im)=9999.0
     else if (im.eq.2) then ! solid

      if (axis_dir.eq.0) then
       LS(im) = DIST_CONE(x)
      else if (axis_dir.eq.1) then
       LS(im) = DIST_CYLINDER_MEHDI(x)
      else if (axis_dir.eq.2) then
       LS(im) = DIST_SPHERE_MEHDI(x)
      else
       print *,"axis_dir invalid"
       stop
      endif

!!$      print *,"x:", x
!!$      print *,"t:", t
!!$      print *,"im:", im
!!$      print *,"LS(im):", LS(im)
     else
      print *,"im invalid"
      stop
     endif
    enddo ! im=1..num_materials
  else
    print *,"num_materials or probtype invalid"
    stop
  endif

  return
  end subroutine CONE3D_LS

  !****************************************************
  function DIST_CONE(P)
   ! Returns the signed distance function to the cone surface
   ! Inside the cone > 0
   ! Outside the cone < 0
   IMPLICIT NONE

   REAL_T M(SDIM)
   REAL_T P(SDIM)
   REAL_T AP(SDIM)
   REAL_T AM(SDIM)
   REAL_T PM(SDIM)
   REAL_T angle
   REAL_T DIST_CONE

   REAL_T :: PI=3.141592653589793

   AP = P-APEX

#ifdef BL_USE_DOUBLE
   angle = DACOS(DOT_PRODUCT(AP,BISECT)/NORM2(AP))
#else
   angle = ACOS(DOT_PRODUCT(AP,BISECT)/NORM2(AP))
#endif

   if (angle.ge.(PI/2+THETA_C)) then
    ! There is no normal to the cone surface so return
    ! distance to the cone apex
    DIST_CONE = -NORM2(AP)
   else
    ! distance to cone surface
    ! M is the intersection of the normal to the cone surface from P  
    ! and the bisect semiline 
    ! d = (P-Apex) . Bisect 
    ! M = Apex + d Biset
    ! LS = |AM| sin(theta_c) - |PM| sin(theta_c)
    M = APEX + DOT_PRODUCT(AP,BISECT) * BISECT
    AM = M - APEX
    PM = M - P
#ifdef BL_USE_DOUBLE
    DIST_CONE = NORM2(AM)*DSIN(THETA_C) - NORM2(PM)*DCOS(THETA_C) 
#else
    DIST_CONE = NORM2(AM)*SIN(THETA_C) - NORM2(PM)*COS(THETA_C) 
#endif
   end if

   return

  end function DIST_CONE


   function DIST_CYLINDER_MEHDI(P)
   ! Returns the signed distance function to the cylinder surface
   ! Inside the cylinder > 0
   ! Outside the cylinder < 0
   IMPLICIT NONE

   REAL_T P(SDIM)
   REAL_T A(SDIM)
   REAL_T AP(SDIM)
   REAL_T CP(SDIM)
   REAL_T DIST_CYLINDER_MEHDI
   REAL_T local_dotprod
   INTEGER_T dir

   local_dotprod=zero
   do dir=1,SDIM
    CP(dir) = P(dir)-MEHDI_CENTER(dir)
    local_dotprod=local_dotprod+CP(dir)*MEHDI_AXIS(dir)
   enddo
   do dir=1,SDIM
    A(dir) = MEHDI_CENTER(dir) + local_dotprod * MEHDI_AXIS(dir)
   enddo
   local_dotprod=zero
   do dir=1,SDIM
    AP(dir) = P(dir)-A(dir)
    local_dotprod=local_dotprod+AP(dir)**2
   enddo
   local_dotprod=sqrt(local_dotprod)
   DIST_CYLINDER_MEHDI = MEHDI_RADIUS - local_dotprod
    
   return
   end function DIST_CYLINDER_MEHDI
 
   function DIST_SPHERE_MEHDI(P)
    ! Returns the signed distance function to the sphere surface
    ! Inside the cylinder > 0
    ! Outside the cylinder < 0
    IMPLICIT NONE
 
    REAL_T P(SDIM)
    REAL_T CP(SDIM)
    REAL_T DIST_SPHERE_MEHDI
    REAL_T local_dotprod
    INTEGER_T dir

    local_dotprod=zero
    do dir=1,SDIM
     CP(dir) = P(dir)-MEHDI_CENTER(dir)
     local_dotprod=local_dotprod+CP(dir)**2
    enddo
    local_dotprod=sqrt(local_dotprod)

    DIST_SPHERE_MEHDI = MEHDI_RADIUS - local_dotprod
    
    return
   end function DIST_SPHERE_MEHDI


  !****************************************************
subroutine CONE3D_VEL(x,t,LS,VEL,velsolid_flag,dx,nmat)
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

REAL_T local_PI

if ((velsolid_flag.eq.0).or. &
    (velsolid_flag.eq.1)) then
 ! do nothing
else 
 print *,"velsolid_flag invalid"
 stop
endif

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

local_PI=4.0d0*atan(one)

if ((LS(2).ge.zero).or. &
    (velsolid_flag.eq.1)) then
 ! in solid
 do dir=1,SDIM
  VEL(dir)=zero
 enddo
else if ((LS(2).le.zero).and. &
         (velsolid_flag.eq.0)) then

 ! boundary values
 do dir=1,SDIM
  VEL(dir)=MEHDI_IN_VEL(dir)
 enddo
 if (xblob9.gt.zero) then
  if ((t.ge.zero).and.(t.le.xblob9)) then
   if (yblob8.eq.zero) then
    if (radblob9.gt.zero) then
     VEL(1) = VEL(1)*(one+yblob9*x(2)*exp(-t/radblob9))
    else
     print *,"radblob9 invalid"
     stop
    endif
   else if (yblob8.gt.zero) then
    if (abs(x(2)).le.yblob8) then
     if (radblob9.gt.zero) then
      VEL(1) = VEL(1)*(one+yblob9*exp(-t/radblob9)* &
         sin(local_PI*x(2)/yblob8)*cos(2*local_PI*t/radblob8))
     else
      print *,"radblob9 invalid"
      stop
     endif
    else if (abs(x(2)).ge.yblob8) then
     ! do nothing
    else
     print *,"x(2) invalid"
     stop
    endif
   else
    print *,"yblob8 invalid"
    stop
   endif
  else if (t.ge.xblob9) then
   ! do nothing
  else
   print *,"t invalid"
   stop
  endif
 else if (xblob9.eq.zero) then
  ! do nothing
 else
  print *,"xblob9 invalid"
  stop
 endif

else
 print *,"LS(2) bust"
 print *,"t=",t
 print *,"num_materials=",num_materials
 print *,"LS(1),LS(2) ",LS(1),LS(2)
 print *,"velsolid_flag ",velsolid_flag
 stop
endif

return 
end subroutine CONE3D_VEL

  !****************************************************
subroutine CONE3D_PRES(x,t,LS,PRES,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(nmat)
REAL_T, INTENT(out) :: PRES

REAL_T gravity_dz

if (SDIM.eq.2) then
 gravity_dz=x(SDIM)-probhiy
else if (SDIM.eq.3) then
 if (1.eq.0) then
  gravity_dz=x(SDIM)-probhiz
 else if (1.eq.1) then
  gravity_dz=x(2)-probhiy
 else
  print *,"do not know what the gravity orientation is"
  stop
 endif
else
 print *,"dimension bust"
 stop
endif

PRES=-fort_denconst(1)*abs(gravity_vector(SDIM))*gravity_dz

return 
end subroutine CONE3D_PRES

!****************************************************
subroutine CONE3D_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
   use probcommon_module
   IMPLICIT NONE

INTEGER_T, INTENT(in) :: bcflag !0=called from initialize  1=called from bc
INTEGER_T, INTENT(in) :: nmat
INTEGER_T, INTENT(in) :: nstate_mat
REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(nmat)
REAL_T, INTENT(out) :: STATE(nmat*nstate_mat)
INTEGER_T im,ibase

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

if ((num_materials.eq.2).and. &
    (num_state_material.eq.2).and. &
    (probtype.eq.222)) then
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
 print *,"num_materials,num_state_material, or probtype invalid"
 stop
endif

return
end subroutine CONE3D_STATE

  !****************************************************
  ! dir=1..sdim  side=1..2
subroutine CONE3D_LS_BC(xwall,xghost,t,LS, &
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
 call CONE3D_LS(xghost,t,LS,nmat)
else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine CONE3D_LS_BC

  !****************************************************
  ! dir=1..sdim  side=1..2 veldir=1..sdim
subroutine CONE3D_VEL_BC(xwall,xghost,t,LS, &
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

 call CONE3D_VEL(xghost,t,LS,local_VEL,velsolid_flag,dx,nmat)
 VEL=local_VEL(veldir)

else
 print *,"dir,side, or veldir invalid"
 stop
endif

return
end subroutine CONE3D_VEL_BC

!****************************************************
! dir=1..sdim  side=1..2
subroutine CONE3D_PRES_BC(xwall,xghost,t,LS, &
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

 call CONE3D_PRES(xghost,t,LS,PRES,nmat)

else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine CONE3D_PRES_BC


function is_CONE3D_overlay(nmat,im)
use probcommon_module
IMPLICIT NONE

INTEGER_T is_CONE3D_overlay
INTEGER_T, INTENT(in) :: nmat,im

if (nmat.eq.num_materials) then
 if (num_materials.eq.2) then
  if ((im.ge.1).and.(im.le.nmat)) then
   if (im.eq.2) then 
    is_CONE3D_overlay=1
   else
    is_CONE3D_overlay=0
   endif
  else
   print *,"im invalid in is_CONE3D_overlay"
   stop
  endif
 else
  print *,"num_materials invalid in is_CONE3D_overlay"
  stop
 endif
else
 print *,"nmat invalid in is_CONE3D_overlay"
 stop
endif
 
return
end function is_CONE3D_overlay


!****************************************************
! dir=1..sdim  side=1..2
subroutine CONE3D_STATE_BC(xwall,xghost,t,LS, &
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
 call CONE3D_STATE(xghost,t,LS,local_STATE, &
   local_bcflag,nmat,num_state_material)
 ibase=(im-1)*num_state_material
 STATE=local_STATE(ibase+istate)
 im_crit=1
 do im_loop=2,num_materials
  if (LS(im_loop).gt.LS(im_crit)) then
   im_crit=im_loop
  endif
 enddo

 do im_loop=1,num_materials
  if (is_CONE3D_overlay(num_materials,im_loop).eq.1) then
   if (LS(im_loop).ge.zero) then
    im_crit=im_loop
   else if (LS(im_loop).le.zero) then
    ! do nothing
   else
    print *,"LS(im_loop) invalid"
    stop
   endif
  else if (is_CONE3D_overlay(num_materials,im_loop).eq.0) then
   ! do nothing
  else
   print *,"is_CONE3D_overlay(num_materials,im_loop) invalid"
   stop
  endif
 enddo ! im_loop=1,num_materials

 ibase=(im_crit-1)*num_state_material
 STATE_merge=local_STATE(ibase+istate)
else
 print *,"istate invalid"
 stop
endif

return
end subroutine CONE3D_STATE_BC

end module CONE3D_module
