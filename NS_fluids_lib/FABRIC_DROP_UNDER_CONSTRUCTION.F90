#undef BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "AMReX_REAL.H"
#include "AMReX_CONSTANTS.H"
#include "AMReX_SPACE.H"
#include "AMReX_BC_TYPES.H"

#include "AMReX_ArrayLim.H"

! #define DEBUG_FABRIC_DROP
#define FABRIC_DROP_PROB_TYPE 7001


#if (AMREX_SPACEDIM==3)
#define SDIM 3
#elif (AMREX_SPACEDIM==2)
#define SDIM 2
#else
print *,"dimension bust"
stop
#endif

! probtype==FABRIC_DROP_PROB_TYPE (see run2d/inputs.FABRIC_DROP)
module FABRIC_DROP_MODULE
implicit none 
REAL(KIND=8),PARAMETER        :: pi=4.0d0*atan(1.0d0)
REAL(KIND=8),PARAMETER        :: a_wavy=0.4d0
REAL(KIND=8),PARAMETER        :: r_1=0.14d0  ! radius of wavy thread
REAL(KIND=8),PARAMETER        :: r_2=0.14d0  ! radius of flat thread
INTEGER,PARAMETER             :: N1=4   ! wavy threads
INTEGER,PARAMETER             :: N2=4   ! Straight threads
integer,parameter             :: P=192  ! number of partition points
REAL(KIND=8),PARAMETER        :: omega=pi/(3.0d0/4.0d0)

contains

subroutine l2normd(s,x1,x2, x1x2norm)
implicit none

integer,intent(in)       :: s
real(kind=8),intent(in)  :: x1(s),x2(s)
real(kind=8)             :: x1x2norm

integer                  :: i
real(kind=8),allocatable :: diff(:)

x1x2norm = 0.0d0
allocate(diff(s))
do i = 1,s
 diff(i) = x1(i)-x2(i)
enddo

do i = 1,s
 x1x2norm = x1x2norm + diff(i)**2.0d0
enddo

x1x2norm = sqrt(x1x2norm)

deallocate(diff)

end subroutine l2normd



subroutine find_xc(la,lb,flag,xz,xc)
! flag =1 sinecurve   =2 cosinecurve
implicit none

real(kind=8),intent(in) :: la,lb
integer,intent(in)      :: flag
real(kind=8),intent(in) :: xz(2)
integer                 :: k
real(kind=8)            :: spl(2,P+1)
real(kind=8)            :: dist,dtemp
real(kind=8),intent(out):: xc(3)
real(kind=8)            :: spltemp(2)

   do k=1,P+1
     spl(1,k)=la+(k-1)*(lb-la)/real(P,8)
    if(flag.eq.0)then
     spl(2,k)=a_wavy*sin(omega*(spl(1,k)+1.5d0)) 
    elseif(flag.eq.1)then
     spl(2,k)=-a_wavy*sin(omega*(spl(1,k)+1.5d0))
    else
     print *,"flag invalid"
     stop
    endif
   enddo
   dist=1.0e+8
   xc=0.0d0
   do k=1,P+1
    spltemp(1)=spl(1,k)
    spltemp(2)=spl(2,k)
    call l2normd(2,xz,spltemp,dtemp)
    if(dtemp .lt. dist)then
     dist=dtemp
     xc(1)=spl(1,k)
     xc(3)=spl(2,k)
    endif
   enddo  

end subroutine find_xc
! do any initial preparation needed
subroutine INIT_FABRIC_DROP_MODULE()
use probcommon_module
IMPLICIT NONE
REAL_T R
REAL_T internal_x(3)
INTEGER_T ithread
INTEGER_T i,j,k
INTEGER_T N, N_max


#ifdef DEBUG_FABRIC_DROP
INTEGER_T inode
#endif



if ((probtype.ne.FABRIC_DROP_PROB_TYPE).or.(SDIM.ne.3)) then
 print *,"probtype or SDIM invalid!"
 stop
endif



print *,"... done!"


return
end subroutine INIT_FABRIC_DROP_MODULE


subroutine FABRIC_DROP_LS(x,t,LS,nmat)
use probcommon_module
IMPLICIT NONE


INTEGER_T, intent(in) :: nmat
REAL_T, intent(in) :: x(SDIM)
REAL_T, intent(in) :: t
REAL_T, intent(out) :: LS(nmat)
real(kind=8)         :: xlo,xhi,ylo,yhi,zlo,zhi
real(kind=8)         :: spl(2,P+1)
real(kind=8)         :: hN1,hN2,xctemp
real(kind=8)         :: xy(3),xc(3),xz(2),xh(3)
integer              :: i,j,k,l,ky,kx
real(kind=8)         :: lstemp1,lstemp2
real(kind=8)         :: a,b
real(kind=8)         :: dtemp,dist
real(kind=8)         :: ynf(N1+2)  ! wavy thread
real(kind=8)         :: xnf(N2+2)  ! straight thread
integer              :: fnflag


if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif


LS(1)=radblob-sqrt((x(1)-xblob)**2+(x(2)-yblob)**2+(x(SDIM)-zblob)**2)
LS(2)=-LS(1)


xlo=problox
xhi=probhix
ylo=probloy
yhi=probhiy
zlo=probloz
zhi=probhiz

hN1=(yhi-ylo)/real(N1,8)
hN2=(xhi-xlo)/real(N2,8)

!hN1=(prob_hi(2)-prob_lo(2))/real(N1,8)
!hN2=(prob_hi(1)-prob_lo(1))/real(N2,8)

do i=1,N1+2
 ynf(i)=ylo-0.5d0*hN1+(i-1)*hN1 
enddo
do i=1,N2+2
 xnf(i)=xlo-0.5d0*hN2+(i-1)*hN2
enddo

!do i=1,N1+2
! ynf(i)=prob_lo(2)-0.5d0*hN1+(i-1)*hN1 
!enddo
!do i=1,N2+2
! xnf(i)=prob_lo(1)-0.5d0*hN2+(i-1)*hN2
!enddo

LS(3)=0.0d0

   xy(1)=x(1)
   xy(2)=x(2)
   xy(3)=x(3)
   xz(1)=x(1)
   xz(2)=x(3)

  do ky=1,N1+2
   if( x(2).ge.ynf(ky)-0.5d0*hN1 .and. x(2).lt. ynf(ky)+0.5d0*hN1)then
    do kx=1,N2+2
     if( x(1).ge.xnf(kx) .and. x(1).lt.xnf(kx)+hN2)then
      a=xnf(kx)
      b=xnf(kx)+hN2      
      exit
     endif
    enddo
    xctemp=ynf(ky)
    fnflag=mod(ky,2)
    exit
   endif
  enddo
 
  call find_xc(a,b,fnflag,xz,xc)
  xc(2)=xctemp
  call l2normd(3,xy,xc,lstemp1)
  lstemp1=lstemp1-r_1


    do kx=1,N2+2
     if( x(1).ge.xnf(kx)-0.5d0*hN2 .and. x(1).lt.xnf(kx)+0.5d0*hN2)then
      xh(1)=xnf(kx)
      xh(2)=xy(2)
      xh(3)=0.0d0     
      exit
     endif
    enddo
   call l2normd(3,xy,xh,lstemp2)
   lstemp2=lstemp2-r_2

 if(lstemp1.lt.0.0d0)then
  LS(3)=lstemp1
 elseif(lstemp2.lt.0.0d0)then
  LS(3)=lstemp2
 else
  LS(3)=min(lstemp1,lstemp2)
 endif


return
end subroutine FABRIC_DROP_LS



subroutine FABRIC_DROP_VEL(x,t,LS,VEL,velsolid_flag,dx,nmat)
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

if (adv_dir.eq.SDIM) then
  ! material 3 is the fabric
  ! velsolid_flag==1 if initializing the solid velocity.
 if ((LS(3).ge.zero).or.(velsolid_flag.eq.1)) then
  ! in solid
  do dir=1,SDIM
   VEL(dir)=zero
  enddo
 else if ((LS(3).le.zero).and. &
          (velsolid_flag.eq.0)) then
  ! material 1 is the drop
  if ((LS(1).ge.zero).or. &
      (LS(1).ge.-two*dx(1))) then
   ! in drop
   do dir=1,SDIM
    VEL(dir)=zero
   enddo
   VEL(SDIM)=-abs(advbot)
   ! material 2 is the gas
  else if (LS(2).ge.zero) then 
   do dir=1,SDIM
    VEL(dir)=zero
   enddo
  else
   print *,"LS bust"
   stop
  endif
 else
  print *,"LS(3) or velsolid bust"
  stop
 endif

else
 print *,"expecting adv_dir = SDIM in FABRIC_DROP"
 stop
endif

return 
end subroutine FABRIC_DROP_VEL






! ! These next routines only used for compressible materials.
! !***********************************************
! ! compressible material functions for (ns.material_type = 24)
! subroutine EOS_FABRIC_DROP(rho,internal_energy,pressure, &
!   imattype,im)
!  IMPLICIT NONE
!  INTEGER_T, intent(in) :: imattype,im
!  REAL_T, intent(in) :: rho
!  REAL_T, intent(in) :: internal_energy
!  REAL_T, intent(out) :: pressure

!  if (imattype.eq.24) then
!   pressure=zero
!   print *,"Not expecting compressible material! "
!  else
!   print *,"imattype invalid EOS_FABRIC_DROP"
!   stop
!  endif

!  return
! end subroutine EOS_FABRIC_DROP

! subroutine SOUNDSQR_FABRIC_DROP(rho,internal_energy,soundsqr, &
!   imattype,im)
!  IMPLICIT NONE
!  INTEGER_T, intent(in) :: imattype,im
!  REAL_T, intent(in) :: rho
!  REAL_T, intent(in) :: internal_energy
!  REAL_T, intent(out) :: soundsqr
!  REAL_T pressure

!  if (imattype.eq.24) then
!   call EOS_FABRIC_DROP(rho,internal_energy,pressure,imattype,im)
!   soundsqr=zero
!   print *,"Not expecting compressible material! "
!   stop

!  else
!   print *,"imattype invalid SOUNDSQR_FABRIC_DROP"
!   stop
!  endif

!  return
! end subroutine SOUNDSQR_FABRIC_DROP

! subroutine INTERNAL_FABRIC_DROP(rho,temperature, &
!   local_internal_energy, &
!   imattype,im)
!  use global_utility_module
!  IMPLICIT NONE
!  INTEGER_T, intent(in) :: imattype,im
!  REAL_T, intent(in) :: rho
!  REAL_T, intent(in) :: temperature 
!  REAL_T, intent(out) :: local_internal_energy

!  call INTERNAL_default(rho,temperature,local_internal_energy, &
!         imattype,im)

!  return
! end subroutine INTERNAL_FABRIC_DROP

! subroutine TEMPERATURE_FABRIC_DROP(rho,temperature,internal_energy, &
!   imattype,im)
!  use global_utility_module
!  IMPLICIT NONE
!  INTEGER_T, intent(in) :: imattype,im
!  REAL_T, intent(in) :: rho
!  REAL_T, intent(out) :: temperature 
!  REAL_T, intent(in) :: internal_energy

!  call TEMPERATURE_default(rho,temperature,internal_energy, &
!         imattype,im)

!  return
! end subroutine TEMPERATURE_FABRIC_DROP

! This routine will not effect the simulation since
! all of the domain BC will be no-slip.


subroutine FABRIC_DROP_PRES(x,t,LS,PRES,nmat)
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
PRES=outflow_pressure

return 
end subroutine FABRIC_DROP_PRES


subroutine FABRIC_DROP_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
use probcommon_module
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

if ((num_materials.eq.3).and. &
    (num_state_material.ge.2).and. &
    (probtype.eq.FABRIC_DROP_PROB_TYPE)) then
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
  do n=1,num_species_var
   STATE(ibase+2+n)=fort_speciesconst((n-1)*num_materials+im)
  enddo
 enddo ! im=1..num_materials
else
 print *,"num_materials,num_state_material, or probtype invalid"
 stop
endif
 
return
end subroutine FABRIC_DROP_STATE

 ! dir=1..sdim  side=1..2
subroutine FABRIC_DROP_LS_BC(xwall,xghost,t,LS, &
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
 call FABRIC_DROP_LS(xghost,t,LS,nmat)
else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine FABRIC_DROP_LS_BC


 ! dir=1..sdim  side=1..2 veldir=1..sdim
subroutine FABRIC_DROP_VEL_BC(xwall,xghost,t,LS, &
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

 call FABRIC_DROP_VEL(xghost,t,LS,local_VEL,velsolid_flag,dx,nmat)
 VEL=local_VEL(veldir)

else
 print *,"dir,side, or veldir invalid"
 stop
endif

return
end subroutine FABRIC_DROP_VEL_BC

 ! dir=1..sdim  side=1..2
subroutine FABRIC_DROP_PRES_BC(xwall,xghost,t,LS, &
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

 call FABRIC_DROP_PRES(xghost,t,LS,PRES,nmat)

else
 print *,"dir or side invalid"
 stop
endif

return
end subroutine FABRIC_DROP_PRES_BC

 ! dir=1..sdim  side=1..2
subroutine FABRIC_DROP_STATE_BC(xwall,xghost,t,LS, &
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
 call FABRIC_DROP_STATE(xghost,t,LS,local_STATE, &
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
end subroutine FABRIC_DROP_STATE_BC

subroutine FABRIC_DROP_HEATSOURCE( &
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

if ((num_materials.eq.3).and.(probtype.eq.FABRIC_DROP_PROB_TYPE)) then
 heat_source=zero
else
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine FABRIC_DROP_HEATSOURCE

end module FABRIC_DROP_MODULE
