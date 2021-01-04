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

 ! probtype==2011
 module YAOHONG_INKJET_module

 implicit none

! time vs. pressure, time unit(ms), pressure unit (kPa)
 INTEGER_T, parameter :: N_pressure=500
 REAL_T :: t_pressure(N_pressure), press_nozzle(N_pressure)

 contains



subroutine pressure_input()
implicit none

! time vs. pressure, time unit(ms), pressure unit (kPa)

INTEGER_T  :: i

print *,"opening press_YAOHONG.in"

open(30,file='press_YAOHONG.in')
do i=1,N_pressure
   read(30,*) t_pressure(i), press_nozzle(i)
end do
close(30)
end subroutine pressure_input


  ! initial preparations
  subroutine INIT_YAOHONG_INKJET_MODULE()
   use probcommon_module
   IMPLICIT NONE

   call pressure_input()

   return
  end subroutine INIT_YAOHONG_INKJET_MODULE


subroutine press_interp(t_input,t, press, N, p_output)                        
implicit none
INTEGER_T   :: N
REAL_T         :: t(N), press(N)
REAL_T         :: t_input, p_output
INTEGER_T              :: m  
m=floor(t_input*5)
if (m<N) then
   p_output=press(m)+(t_input-t(m))/(t(m+1)-t(m))*(press(m+1)-press(m))
else
   write(*,*) "press_interp: input time is out of initial set up"
end if
end subroutine press_interp

subroutine LS_geometry(x,y,ls)
implicit none

! setup the computational domain with level set function
! 11.5x1.5 (mm) domain
! ls for nozzel geometry

REAL_T   :: x, y, ls
REAL_T   :: d1,d2
REAL_T   :: w1=0.3d0, w2=1.0d0, w=1.5d0 
REAL_T   :: l1=4.0d0,  l2 =2.5d0, l3=4.5d0, l4=0.5d0

if (y<=l1.and.y>=0) then
   if (x<=w1) then
      ls=sqrt((x-w1)**2+(l1-y)**2)
   else
      ls=l1-y
   end if
else if (y-l1<=l2.and.y>=0) then
   if (x<=w1) then
      ls= w1-x
   else 
      d1=x-w1
      d2=y-l1
      ls=-min(d1,d2)
   end if
else if (y-l1-l2<l3) then
   ls=-(45.0*x-7.0*y+32)/sqrt(45.0**2+7.0**2)
else if (y-l1-l2-l3<=l4.and.y>=0) then
   if (x<w2) then
      d1=-(45.0*x-7.0*y+32)/sqrt(45.0**2+7.0**2)
      d2=w2-x
      ls=min(d1,d2)
   else
      ls=w2-x
   end if
else
   write(*,*) "ls geometry, (x,y) is out of the domain",x,y
end if
end subroutine LS_geometry

subroutine LS_air_water(x,y,ls)
implicit none
REAL_T   :: x, y, ls
REAL_T   :: d1,d2
REAL_T   :: w1=0.3d0, w2=1.0d0, w=1.5d0 
REAL_T   :: l1=4.0,  l2 =2.5, l3=4.5, l4=0.5, l=11.5

if (y<=l1.and.y>=0) then
   if (x<=w1) then
      ls=y-l1
   else
      ls=-sqrt((x-w1)**2+(l1-y)**2)
   end if
else if (y-l1<=l2.and.y>=0) then
   if (x<=w1) then
      d1=w1-x
      d2=y-l1
      ls= min(d1,d2)
   else
      ls=w1-x
   end if
else if (y-l1-l2<=l3.and.y>=0) then
   if (x<=w1) then
      d1=y-l1
      d2=-(45.0*x-7.0*y+32)/sqrt(45.0**2+7.0**2)
      ls= min(d1,d2)
   else
      ls=-(45.0*x-7.0*y+32)/sqrt(45.0**2+7.0**2)
   end if
else if (y-l1-l2-l3<=l4.and.y>=0) then
   if (x<=w2) then
      d1=-(45.0*x-7.0*y+32)/sqrt(45.0**2+7.0**2)
      d2=w2-x
      ls=min(d1,d2)
   else
      ls=w2-x
   end if
else
   write(*,*) "ls air_water, (x,y) is out of the domain",x,y
end if
end subroutine LS_air_water




  ! fluids tessellate the domain, solids are immersed. 
  subroutine YAOHONG_INKJET_LS(x,t,LS,nmat)
  use probcommon_module
  IMPLICIT NONE

  INTEGER_T, intent(in) :: nmat
  REAL_T, intent(in) :: x(SDIM)
  REAL_T, intent(in) :: t
  REAL_T, intent(out) :: LS(nmat)
  INTEGER_T im

  if (nmat.eq.num_materials) then
   ! do nothing
  else
   print *,"nmat invalid"
   stop
  endif

  if ((num_materials.eq.3).and.(probtype.eq.2011)) then
    do im=1,num_materials
     if (im.eq.1) then
      call LS_air_water(x(1),x(2),LS(1))
     else if (im.eq.2) then
      call LS_air_water(x(1),x(2),LS(2))
      LS(2)=-LS(2)
     else if (im.eq.3) then
      call LS_geometry(x(1),x(2),LS(3))
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
  end subroutine YAOHONG_INKJET_LS


  !****************************************************
subroutine YAOHONG_INKJET_VEL(x,t,LS,VEL,velsolid_flag,dx,nmat)
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
  VEL(dir)=zero
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
end subroutine YAOHONG_INKJET_VEL

  !****************************************************
subroutine YAOHONG_INKJET_PRES(x,t,LS,PRES,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(nmat)
REAL_T, intent(out) :: PRES

REAL_T gravity_dz


PRES=zero

return 
end subroutine YAOHONG_INKJET_PRES

!****************************************************
subroutine YAOHONG_INKJET_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
   use probcommon_module
   IMPLICIT NONE

INTEGER_T, intent(in) :: bcflag !0=called from initialize  1=called from bc
INTEGER_T, intent(in) :: nmat
INTEGER_T, intent(in) :: nstate_mat
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(nmat)
REAL_T, intent(out) :: STATE(nmat*nstate_mat)
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

if ((num_materials.eq.3).and. &
    (num_state_material.eq.2).and. &
    (probtype.eq.2011)) then
 do im=1,num_materials

  ibase=(im-1)*num_state_material
  STATE(ibase+1)=fort_denconst(im)

  if (t.eq.zero) then
   STATE(ibase+2)=fort_initial_temperature(im)
  else if (t.gt.zero) then
   STATE(ibase+2)=fort_tempconst(im)
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
end subroutine YAOHONG_INKJET_STATE

  !****************************************************
  ! dir=1..sdim  side=1..2
subroutine YAOHONG_INKJET_LS_BC(xwall,xghost,t,LS, &
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
REAL_T, intent(in) :: dx(SDIM)

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif

if ((dir.ge.1).and.(dir.le.SDIM).and. &
    (side.ge.1).and.(side.le.2)) then
 call YAOHONG_INKJET_LS(xghost,t,LS,nmat)
else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine YAOHONG_INKJET_LS_BC

  !****************************************************
  ! dir=1..sdim  side=1..2 veldir=1..sdim
subroutine YAOHONG_INKJET_VEL_BC(xwall,xghost,t,LS, &
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

 call YAOHONG_INKJET_VEL(xghost,t,LS,local_VEL,velsolid_flag,dx,nmat)
 VEL=local_VEL(veldir)

else
 print *,"dir,side, or veldir invalid"
 stop
endif

return
end subroutine YAOHONG_INKJET_VEL_BC

!****************************************************
! dir=1..sdim  side=1..2
subroutine YAOHONG_INKJET_PRES_BC(xwall,xghost,t,LS, &
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

 call YAOHONG_INKJET_PRES(xghost,t,LS,PRES,nmat)

else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine YAOHONG_INKJET_PRES_BC


function is_YAOHONG_INKJET_overlay(nmat,im)
use probcommon_module
IMPLICIT NONE

INTEGER_T is_YAOHONG_INKJET_overlay
INTEGER_T, intent(in) :: nmat,im

if (nmat.eq.num_materials) then
 if (num_materials.eq.3) then
  if ((im.ge.1).and.(im.le.nmat)) then
   if (im.eq.3) then 
    is_YAOHONG_INKJET_overlay=1
   else
    is_YAOHONG_INKJET_overlay=0
   endif
  else
   print *,"im invalid in is_YAOHONG_INKJET_overlay"
   stop
  endif
 else
  print *,"num_materials invalid in is_YAOHONG_INKJET_overlay"
  stop
 endif
else
 print *,"nmat invalid in is_YAOHONG_INKJET_overlay"
 stop
endif
 
return
end function is_YAOHONG_INKJET_overlay


!****************************************************
! dir=1..sdim  side=1..2
subroutine YAOHONG_INKJET_STATE_BC(xwall,xghost,t,LS, &
 STATE,STATE_merge,STATE_in,im,istate,dir,side,dx,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: xwall
REAL_T, intent(in) :: xghost(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(in) :: LS(nmat)
REAL_T local_STATE(nmat*num_state_material)
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
 call YAOHONG_INKJET_STATE(xghost,t,LS,local_STATE, &
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
  if (is_YAOHONG_INKJET_overlay(num_materials,im_loop).eq.1) then
   if (LS(im_loop).ge.zero) then
    im_crit=im_loop
   else if (LS(im_loop).le.zero) then
    ! do nothing
   else
    print *,"LS(im_loop) invalid"
    stop
   endif
  else if (is_YAOHONG_INKJET_overlay(num_materials,im_loop).eq.0) then
   ! do nothing
  else
   print *,"is_YAOHONG_INKJET_overlay(num_materials,im_loop) invalid"
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
end subroutine YAOHONG_INKJET_STATE_BC

end module YAOHONG_INKJET_module
