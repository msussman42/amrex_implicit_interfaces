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
integer,parameter             :: P=96  ! number of partition points
REAL(KIND=8),PARAMETER        :: omega=pi/(3.0d0/4.0d0)

REAL_T, allocatable, dimension(:,:,:) :: internal_thread_ls
REAL_T internal_dx(3)

contains

subroutine l2normd(s,x1,x2, x1x2norm)
implicit none

integer,INTENT(in)       :: s
real(kind=8),INTENT(in)  :: x1(s),x2(s)
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

subroutine dist_point_to_lined(sd,p1,p2,x,pout,dist)
implicit none
! represent the line in parametric form,(v = f(s))
! v^x = x1 + (x2-x1)s
! v^y = y1 + (y2-y1)s
! v^z = z1 + (z2-z1)s
! 
! if the closest point on the line to point x  is outside p1 -- p2, 
! --------> return the distance from x either to p1 or p2 which is shorter.
! otherwise
! --------> return the distance from x to the cloest point

integer,INTENT(in)           :: sd
real(kind=8),INTENT(in)      ::  p1(sd),p2(sd),x(sd)
real(kind=8),INTENT(out)     ::  dist
real(kind=8),INTENT(out)     ::  pout(sd)

real(kind=8)                 :: diff10,diff21,diffx
real(kind=8),allocatable     :: x10(:), x21(:)
integer                      :: i
real(kind=8)                 :: s


dist = 0.0d0

allocate(x10(sd),x21(sd))
do i = 1,sd
 x10(i) = p1(i) - x(i)
 x21(i) = p2(i) - p1(i)
enddo

if (maxval(abs(x21)) .lt. 10d-8)then
 print *,"p1 and p2 are coincide with each other",p1,p2
 stop
endif

call l2normd(sd, p1, x, diff10)
call l2normd(sd, p2, p1, diff21)

!print *,"diff10",diff10
!print *,"diff21",diff21

s = -1.0d0*(dot_product(x10,x21))/(diff21**2.0d0)

!write(13,*) "s=",s

if(s .gt. 1.0d0)then
 call l2normd(sd, p2, x,dist)
 do i=1,sd
  pout(i)=p1(i)
 enddo
elseif(s .lt. 0.0d0)then
 call l2normd(sd,p1,x,dist)
 do i=1,sd
  pout(i)=p2(i)
 enddo
else
! if(abs((diff10**2.0d0 )*(diff21**2.0d0) - & 
!        (dot_product(x10,x21))**2.0d0) .lt. 1.0e-10)then
!   dist= 0.0d0
! else
  dist = sqrt(((diff10**2.0d0 )*(diff21**2.0d0) - & 
        (dot_product(x10,x21))**2.0d0)/ &
        (diff21**2.0d0))
  do i=1,sd
   pout(i)=p1(i)+s*(p2(i) - p1(i))
  enddo
! endif
endif

deallocate(x10,x21) 

end subroutine dist_point_to_lined


subroutine find_xc(la,lb,flag,xz,xc)
! flag =1 sinecurve   =2 cosinecurve
implicit none

real(kind=8),INTENT(in) :: la,lb
integer,INTENT(in)      :: flag
real(kind=8),INTENT(in) :: xz(2)
integer                 :: k
real(kind=8),allocatable         :: spl(:,:)
real(kind=8)            :: dist,dtemp
real(kind=8),INTENT(out):: xc(3)
real(kind=8)            :: spltemp(2),pout(2)


allocate(spl(2,P+1))
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

!   do k=1,P+1
!    spltemp(1)=spl(1,k)
!    spltemp(2)=spl(2,k)
!    call l2normd(2,xz,spltemp,dtemp)
!    if(dtemp .lt. dist)then
!     dist=dtemp
!     xc(1)=spl(1,k)
!     xc(3)=spl(2,k)
!    endif
!   enddo  
   do k=1,P
    call dist_point_to_lined(2,spl(:,k),spl(:,k+1),xz,pout,dtemp)
    if(dtemp .lt. dist)then
     dist=dtemp
     xc(1)=pout(1)
     xc(3)=pout(2)
    endif
   enddo 

 deallocate(spl)
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

! Evaluate level set function on internal fine mesh
allocate(internal_thread_ls(IDNINT(xblob5),IDNINT(yblob5),IDNINT(zblob5)))

internal_dx(1)=(xblob4-xblob3)/xblob5
internal_dx(2)=(yblob4-yblob3)/yblob5
internal_dx(3)=(zblob4-zblob3)/zblob5

print *,"Internal dx for thread LS calculation:", internal_dx
print *,"Calculating thread level set on internal mesh..."

do k=1,IDNINT(zblob5)
#ifdef DEBUG_FABRIC_DROP
 print *,"k=", k
#endif
 do j=1,IDNINT(yblob5)
  do i=1,IDNINT(xblob5)
   internal_x(1)=xblob3+(i-half)*internal_dx(1)
   internal_x(2)=yblob3+(j-half)*internal_dx(2)
   internal_x(3)=zblob3+(k-half)*internal_dx(3)

!   internal_thread_ls(i,j,k) = DIST_THREADS(internal_x)
   call FABRIC_LS(internal_x,internal_thread_ls(i,j,k))
   enddo
 enddo
enddo



print *,"... done!"


return
end subroutine INIT_FABRIC_DROP_MODULE

subroutine FABRIC_DROP_LS(x,t,LS,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(out) :: LS(nmat)

REAL_T x_d(3),x_0(3),x_p(3)
INTEGER_T ind(3)
REAL_T c_000, c_001, c_010, c_011, c_100, c_101, c_110, c_111
REAL_T c_00, c_01, c_10, c_11
REAL_T c_0, c_1, c
INTEGER_T dir

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif


LS(1)=radblob-sqrt((x(1)-xblob)**2+(x(2)-yblob)**2+(x(SDIM)-zblob)**2)
LS(2)=-LS(1)

! Direct thread level set calculation 
! LS(3)= DIST_THREADS(x)

! Approximate level set calcualtion from internal fine mesh
if ((x(SDIM).lt.(zblob3+internal_dx(SDIM))).or. &
    (x(SDIM).gt.(zblob4-internal_dx(SDIM)))) then
 ! use large level set value
 LS(3)=-abs(radblob3)
else
 ! assuming [xblob3,yblob3]:[xblob4,yblob4] is the repeating pattern
 x_p(1)=modulo(x(1),xblob4)
 x_p(2)=modulo(x(2),yblob4)
 if (SDIM.eq.3) then
  x_p(SDIM)=x(SDIM)
 endif

 ind(1)=IDINT((x_p(1)-xblob3)/internal_dx(1))
 ind(2)=IDINT((x_p(2)-yblob3)/internal_dx(2))
 if (SDIM.eq.3) then
  ind(SDIM)=IDINT((x_p(SDIM)-zblob3)/internal_dx(SDIM))
 endif
 do dir=1,SDIM
  if (ind(dir).lt.1) then
   ind(dir)=1
  endif

  if (ind(dir)+1.gt.UBOUND(internal_thread_ls,dir)) then
   ind(dir)=UBOUND(internal_thread_ls,dir)-1
  endif
 enddo ! dir=1..sdim

 ! x is in the cube {ind, ind+1}
 ! trilinear interpolation
 ! (https://en.wikipedia.org/wiki/Trilinear_interpolation)
 x_0(1)=xblob3+(ind(1)-half)*internal_dx(1)
 x_0(2)=yblob3+(ind(2)-half)*internal_dx(2)
 if (SDIM.eq.3) then
  x_0(SDIM)=zblob3+(ind(SDIM)-half)*internal_dx(SDIM)
 endif

 x_d=(x_p-x_0)/internal_dx ! component-wise division

 c_000=internal_thread_ls(ind(1)+0, ind(2)+0, ind(3)+0)
 c_100=internal_thread_ls(ind(1)+1, ind(2)+0, ind(3)+0)
 c_010=internal_thread_ls(ind(1)+0, ind(2)+1, ind(3)+0)
 c_110=internal_thread_ls(ind(1)+1, ind(2)+1, ind(3)+0)
 c_001=internal_thread_ls(ind(1)+0, ind(2)+0, ind(3)+1)
 c_101=internal_thread_ls(ind(1)+1, ind(2)+0, ind(3)+1)
 c_011=internal_thread_ls(ind(1)+0, ind(2)+1, ind(3)+1)
 c_111=internal_thread_ls(ind(1)+1, ind(2)+1, ind(3)+1)

 c_00=c_000*(one-x_d(1)) + c_100*x_d(1)
 c_01=c_001*(one-x_d(1)) + c_101*x_d(1)
 c_10=c_010*(one-x_d(1)) + c_110*x_d(1)
 c_11=c_011*(one-x_d(1)) + c_111*x_d(1)

 c_0=c_00*(one-x_d(2)) + c_10*x_d(2)
 c_1=c_01*(one-x_d(2)) + c_11*x_d(2)

 if (SDIM.eq.3) then
  c=c_0*(one-x_d(SDIM)) + c_1*x_d(SDIM)
 endif

 LS(3)=c

endif

! print *,"X, LS", x, LS

return
end subroutine FABRIC_DROP_LS





subroutine FABRIC_LS(x,LS)
use probcommon_module
IMPLICIT NONE


REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(out) :: LS
real(kind=8)         :: xlo,xhi,ylo,yhi,zlo,zhi
real(kind=8)         :: hN1,hN2,xctemp
real(kind=8)         :: xy(3),xc(3),xz(2),xh(3)
integer              :: i,j,k,l,ky,kx
real(kind=8)         :: lstemp1,lstemp2
real(kind=8)         :: a,b
real(kind=8)         :: dtemp,dist
real(kind=8),allocatable :: ynf(:)  ! wavy thread
real(kind=8),allocatable :: xnf(:)  ! straight thread
integer              :: fnflag



allocate(ynf(N1+2),xnf(N2+2))

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

LS=0.0d0

   xy(1)=x(1)
   xy(2)=x(2)
   if (SDIM.eq.3) then
    xy(SDIM)=x(SDIM)
   endif
   xz(1)=x(1)
   if (SDIM.eq.3) then
    xz(2)=x(SDIM)
   endif

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
  lstemp1=r_1-lstemp1


    do kx=1,N2+2
     if( x(1).ge.xnf(kx)-0.5d0*hN2 .and. x(1).lt.xnf(kx)+0.5d0*hN2)then
      xh(1)=xnf(kx)
      xh(2)=xy(2)
      xh(3)=0.0d0     
      exit
     endif
    enddo
   call l2normd(3,xy,xh,lstemp2)
   lstemp2=r_2-lstemp2

 if(lstemp1.ge.0.0d0.and.lstemp2.ge.0.0d0)then
   print *,"lstemp1 and lstemp2 can not both be positive"
   stop
 endif

 LS=min(abs(lstemp1),abs(lstemp2))
 if(lstemp1.ge.0.0d0 .or. lstemp2.ge.0.0d0)then
  ! DO NOTHING
 else
  LS=-1.0d0*LS
 endif

deallocate(xnf,ynf)

return
end subroutine FABRIC_LS



subroutine FABRIC_DROP_VEL(x,t,LS,VEL,velsolid_flag,dx,nmat)
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
!  INTEGER_T, INTENT(in) :: imattype,im
!  REAL_T, INTENT(in) :: rho
!  REAL_T, INTENT(in) :: internal_energy
!  REAL_T, INTENT(out) :: pressure

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
!  INTEGER_T, INTENT(in) :: imattype,im
!  REAL_T, INTENT(in) :: rho
!  REAL_T, INTENT(in) :: internal_energy
!  REAL_T, INTENT(out) :: soundsqr
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
!  INTEGER_T, INTENT(in) :: imattype,im
!  REAL_T, INTENT(in) :: rho
!  REAL_T, INTENT(in) :: temperature 
!  REAL_T, INTENT(out) :: local_internal_energy

!  call INTERNAL_default(rho,temperature,local_internal_energy, &
!         imattype,im)

!  return
! end subroutine INTERNAL_FABRIC_DROP

! subroutine TEMPERATURE_FABRIC_DROP(rho,temperature,internal_energy, &
!   imattype,im)
!  use global_utility_module
!  IMPLICIT NONE
!  INTEGER_T, INTENT(in) :: imattype,im
!  REAL_T, INTENT(in) :: rho
!  REAL_T, INTENT(out) :: temperature 
!  REAL_T, INTENT(in) :: internal_energy

!  call TEMPERATURE_default(rho,temperature,internal_energy, &
!         imattype,im)

!  return
! end subroutine TEMPERATURE_FABRIC_DROP

! This routine will not effect the simulation since
! all of the domain BC will be no-slip.


subroutine FABRIC_DROP_PRES(x,t,LS,PRES,nmat)
use probcommon_module
IMPLICIT NONE

INTEGER_T, INTENT(in) :: nmat
REAL_T, INTENT(in) :: x(SDIM)
REAL_T, INTENT(in) :: t
REAL_T, INTENT(in) :: LS(nmat)
REAL_T, INTENT(out) :: PRES

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

if ((num_materials.eq.3).and. &
    (num_state_material.ge.2).and. &
    (probtype.eq.FABRIC_DROP_PROB_TYPE)) then
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
 
return
end subroutine FABRIC_DROP_STATE

 ! dir=1..sdim  side=1..2
subroutine FABRIC_DROP_LS_BC(xwall,xghost,t,LS, &
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
use global_utility_module
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
INTEGER_T ibase,im_crit
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
 call get_primary_material(LS,im_crit)
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

if ((num_materials.eq.3).and.(probtype.eq.FABRIC_DROP_PROB_TYPE)) then
 heat_source=zero
else
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine FABRIC_DROP_HEATSOURCE

end module FABRIC_DROP_MODULE
