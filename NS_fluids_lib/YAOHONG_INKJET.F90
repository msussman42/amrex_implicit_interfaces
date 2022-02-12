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

 ! probtype==2011
 module YAOHONG_INKJET_module

 implicit none

! time vs. pressure, time unit(ms), pressure unit (kPa)
 INTEGER_T, parameter :: N_pressure=500
 REAL_T :: t_pressure(N_pressure)
 REAL_T :: press_nozzle(N_pressure)

 contains

! 11.5x1.5 mm -> 0.0115 x 0.0015 m  axisymmetric
! surface tension=0.064 N/m
! static angle: 110 degrees
! viscosity liquid: 0.017 Pa s
! viscosity air   : 0.000018 Pa s
! density liquid  : 1180 kg/m^3
! density air     : 1.169 kg/m^3

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
INTEGER_T, intent(in) :: N
REAL_T, intent(in) :: t(N), press(N)
REAL_T, intent(in) :: t_input
REAL_T, intent(out) :: p_output
INTEGER_T           :: m  
REAL_T :: t_ms

! change the unit to ms
! delta t=1/5 ms
! M delta_t = T = 100ms
! M=T/delta t= 100/(1/5)=500
! m delta t = t
! m=t/delta t = 5 t 

if (t_input.ge.0.0d0) then
 ! do nothing
else
 print *,"t_input cannot be negative"
 stop
endif

t_ms=t_input*1000.0

m=floor(t_ms*5.0d0)
m=m+1
if (m.lt.N) then
 p_output=press(m)+(t_ms-t(m))/(t(m+1)-t(m))*(press(m+1)-press(m))
 p_output=p_output*1000.0d0 ! change the units to Pa 
else
 print *,"press_interp: input time is out of initial set up"
 stop
end if
end subroutine press_interp

! ls>0 in fluid, ls<0 in solid
subroutine LS_geometry(x_i,y_i,ls)
implicit none

! setup the computational domain with level set function
! 11.5x1.5 (mm) domain
! ls for nozzel geometry

REAL_T, intent(in) :: x_i,y_i
REAL_T, intent(out) :: ls
REAL_T   :: x, y
REAL_T   :: d1,d2
REAL_T   :: w1=0.3d0, w2=1.0d0, w=1.5d0 
REAL_T   :: l1=4.0d0,  l2 =2.5d0, l3=4.5d0, l4=0.5d0

x=x_i*1000.0d0  ! convert from MKS to mm
y=y_i*1000.0d0  ! convert from MKS to mm

if (SDIM.eq.2) then
 x=abs(x)
else
 print *,"not ready for 3d"
 stop
endif

if (y.le.l1) then
   if (x.le.w1) then
      ls=sqrt((x-w1)**2+(l1-y)**2)
   else if (x.le.w) then
      d1=l1-y
      d2=w-x
      ls=min(d1,d2)
   else if (x.ge.w) then
      ls=w-x
   else
      print *,"x invalid"
      stop
   end if
else if ((y.le.l1+l2).and.(y.ge.l1)) then
   if (x.le.w1) then
      ls= w1-x
   else if ((x.ge.w1).and.(x.le.w)) then
      d1=x-w1
      d2=y-l1
      ls=-min(d1,d2)
   else if (x.ge.w) then
      d1=x-w1
      d2=sqrt((x-w)**2+(y-l1)**2)
      ls=-min(d1,d2)
   else
    print *,"x invalid"
    stop
   end if
else if ((y.le.l1+l2+l3).and.(y.ge.l1+l2)) then
   ls=-(45.0d0*x-7.0d0*y+32.0d0)/sqrt(45.0d0**2+7.0d0**2)
else if ((y.le.l1+l2+l3+l4).and.(y.ge.l1+l1+l2)) then
   if (x.le.w2) then
      d1=-(45.0d0*x-7.0d0*y+32.0d0)/sqrt(45.0d0**2+7.0d0**2)
      d2=w2-x
      ls=min(d1,d2)
   else
      ls=w2-x
   end if
else if (y.ge.l1+l2+l3+l4) then
   ls=w2-x
else
   write(*,*) "ls geometry, (x,y) is out of the domain",x,y
   stop
end if
ls=ls*0.001d0  ! change the unit to meters
end subroutine LS_geometry

! ls>0 in the water, ls<0 in the air
subroutine LS_air_water(x_i,y_i,ls)
implicit none
REAL_T, intent(in) :: x_i,y_i
REAL_T, intent(out) :: ls
REAL_T :: x, y
REAL_T :: d1,d2,e2
REAL_T :: w1=0.3d0, w2=1.0d0, w=1.5d0 
REAL_T :: l1=4.0d0, l2 =2.5d0, l3=4.5d0, l4=0.5d0, l=11.5d0

x=x_i*1000.0d0
y=y_i*1000.0d0

if (SDIM.eq.2) then
 x=abs(x)
else
 print *,"not ready for 3d"
 stop
endif

e2=0.001d0

if (x.le.w1) then
 ls=y-l1
else if (x.ge.w1) then
 ls=y-l1-(x-w1)
else
 print *,"x invalid"
 stop
end if

!   write(*,*) x,y,ls
ls=ls*e2 ! change the unit to meter
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
      LS(3)=-LS(3)
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

if ((num_materials.eq.3).and.(probtype.eq.2011)) then
 if ((LS(nmat).ge.zero).or. &
     (velsolid_flag.eq.1)) then
  ! in solid
  do dir=1,SDIM
   VEL(dir)=zero
  enddo
 else if ((LS(nmat).le.zero).and. &
          (velsolid_flag.eq.0)) then

  ! boundary values
  do dir=1,SDIM
   VEL(dir)=zero
  enddo

 else
  print *,"LS(nmat) bust"
  print *,"t=",t
  print *,"num_materials=",num_materials
  print *,"LS(1),LS(2),LS(3) ",LS(1),LS(2),LS(3)
  print *,"velsolid_flag ",velsolid_flag
  stop
 endif

else
 print *,"num_materials or probtype invalid"
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
REAL_T :: zhi

if (t.ge.zero) then
 ! do nothing
else
 print *,"t invalid"
 stop
endif

if (SDIM.eq.3) then
 zhi=probhiz
else if (SDIM.eq.2) then
 zhi=probhiy
else
 print *,"dimension bust"
 stop
endif
call press_interp(t,t_pressure,press_nozzle, &
        N_pressure,PRES)
if (x(SDIM).ge.half*zhi) then
 ! do nothing
else if (x(SDIM).le.half*zhi) then
 PRES=zero
else
 print *,"x(SDIM) or zhi invalid"
 stop
endif

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

 if (1.eq.0) then
  if (dir.eq.SDIM) then
   print *,"dir,side,x,y,z,t,PRES ", &
    dir,side,xghost(1),xghost(2),xghost(SDIM),t,PRES
  endif
 endif

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
