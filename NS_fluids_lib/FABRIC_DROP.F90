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

! #define DEBUG_FABRIC_DROP
! This configuration file is authored by:
! Mehdi Vahab, Yang Liu, and Kourosh Shoele
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
use amrex_fort_module, only : amrex_real
implicit none 
real(amrex_real),PARAMETER        :: pi=4.0d0*atan(1.0d0)
real(amrex_real),PARAMETER        :: a_wavy=0.28d0
REAL(KIND=8),PARAMETER            :: r_1=0.10d0  ! radius of wavy thread
REAL(KIND=8),PARAMETER            :: r_2=0.10d0  ! radius of flat thread
!real(amrex_real),PARAMETER        :: a_wavy=0.4d0
!real(amrex_real),PARAMETER        :: r_1=0.14d0  ! radius of wavy thread
!real(amrex_real),PARAMETER        :: r_2=0.14d0  ! radius of flat thread

!THIS IS -1.5<x<1.5  -1.5<y<1.5
!integer,PARAMETER             :: N1=4   ! wavy threads
!integer,PARAMETER             :: N2=4   ! Straight threads
!real(amrex_real),PARAMETER    :: omega=2.0d0*pi/1.5d0 ! 2 pi / 1.5

!THIS IS -3.0<x<3.0  -3.0<y<3.0
!8 * .75 = 6
!.75-2*r_1=.55
!xblob4=yblob4=1.5
integer,PARAMETER             :: N1=8   ! wavy threads
integer,PARAMETER             :: N2=8   ! Straight threads
real(amrex_real),PARAMETER    :: omega=2.0d0*pi/1.5d0 ! 2 pi / 1.5

!THIS IS -3.0<x<3.0  -3.0<y<3.0
!12 * .5 = 6
!.5-2*r_1=.3
!xblob4=yblob4=1.0
!integer,PARAMETER             :: N1=12   ! wavy threads
!integer,PARAMETER             :: N2=12   ! Straight threads
!real(amrex_real),PARAMETER    :: omega=2.0d0*pi/1.0d0 !2 pi / 1.0

integer,parameter             :: P=96  ! number of partition points

real(amrex_real), allocatable, dimension(:,:,:) :: internal_thread_ls
real(amrex_real) internal_dx(3)

contains

subroutine l2normd(s,x1,x2, x1x2norm)
implicit none

integer,INTENT(in)       :: s
real(amrex_real),INTENT(in)  :: x1(s),x2(s)
real(amrex_real)             :: x1x2norm

integer                  :: i
real(amrex_real),allocatable :: diff(:)

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
use probcommon_module
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
real(amrex_real),INTENT(in)      ::  p1(sd),p2(sd),x(sd)
real(amrex_real),INTENT(out)     ::  dist
real(amrex_real),INTENT(out)     ::  pout(sd)

real(amrex_real)                 :: diff10,diff21
real(amrex_real),allocatable     :: x10(:), x21(:)
integer                      :: i
real(amrex_real)                 :: s


dist = 0.0d0

allocate(x10(sd),x21(sd))
do i = 1,sd
 x10(i) = p1(i) - x(i)
 x21(i) = p2(i) - p1(i)
enddo

if (maxval(abs(x21)) .lt. EPS_8_4)then
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

real(amrex_real),INTENT(in) :: la,lb
integer,INTENT(in)      :: flag
real(amrex_real),INTENT(in) :: xz(2)
integer                 :: k
real(amrex_real),allocatable         :: spl(:,:)
real(amrex_real)            :: dist,dtemp
real(amrex_real),INTENT(out):: xc(3)
real(amrex_real)            :: pout(2)


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
real(amrex_real) internal_x(3)
integer i,j,k
real(8) :: xblob5_dbl,yblob5_dbl,zblob5_dbl

if ((probtype.ne.FABRIC_DROP_PROB_TYPE).or.(SDIM.ne.3)) then
 print *,"probtype or SDIM invalid!"
 stop
endif

if (axis_dir.eq.0) then

! Evaluate level set function on internal fine mesh
xblob5_dbl=xblob5
yblob5_dbl=yblob5
zblob5_dbl=zblob5
allocate(internal_thread_ls(IDNINT(xblob5_dbl),IDNINT(yblob5_dbl),IDNINT(zblob5_dbl)))

if ((xblob3.eq.0.0d0).and.(yblob3.eq.0.0d0)) then
 !do nothing
else
 print *,"expecting xblob3 and yblob3 = 0.0d0 "
 stop
endif

if (N1.eq.N2) then
 !do nothing
else
 print *,"expecting N1==N2"
 stop
endif

if (abs(N1*xblob4*0.5d0-(probhix-problox)).le.EPS2) then
 !do nothing
else
 print *,"N1 or xblob4 invalid: ",N1,xblob4
 stop
endif
if (abs(N2*yblob4*0.5d0-(probhiy-probloy)).le.EPS2) then
 !do nothing
else
 print *,"N2 or yblob4 invalid: ",N2,yblob4
 stop
endif

if ((xblob4.eq.1.5d0).and.(yblob4.eq.1.5d0)) then
 if (abs(omega-2.0d0*pi/xblob4).le.1.0D-2) then
  !do nothing
 else
  print *,"expecting omega= 2 pi/xblob4 ",omega
  print *,"xblob4 ",xblob4
  stop
 endif
else if ((xblob4.eq.1.0d0).and.(yblob4.eq.1.0d0)) then
 if (abs(omega-2.0d0*pi/xblob4).le.1.0D-2) then
  !do nothing
 else
  print *,"expecting omega= 2 pi/xblob4 ",omega
  print *,"xblob4 ",xblob4
  stop
 endif
else
 print *,"expecting xblob4=yblob4=1.5 or 1.0:",xblob4,yblob4
 stop
endif

internal_dx(1)=(xblob4-xblob3)/xblob5
internal_dx(2)=(yblob4-yblob3)/yblob5
internal_dx(3)=(zblob4-zblob3)/zblob5

print *,"Internal dx for thread LS calculation:", internal_dx
print *,"Calculating thread level set on internal mesh..."

do k=1,IDNINT(zblob5_dbl)
do j=1,IDNINT(yblob5_dbl)
do i=1,IDNINT(xblob5_dbl)
 internal_x(1)=xblob3+(i-half)*internal_dx(1)
 internal_x(2)=yblob3+(j-half)*internal_dx(2)
 internal_x(3)=zblob3+(k-half)*internal_dx(3)

 call FABRIC_LS(internal_x,internal_thread_ls(i,j,k))
enddo
enddo
enddo

else if (axis_dir.eq.1) then !mesh
 !do nothing
else if (axis_dir.eq.2) then !cylinder
 !do nothing
else
 print *,"axis_dir invalid"
 stop
endif

return
end subroutine INIT_FABRIC_DROP_MODULE

subroutine FABRIC_DROP_LS(x,t,LS,nmat)
use probcommon_module
IMPLICIT NONE

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(out) :: LS(nmat)

real(amrex_real) x_d(3),x_0(3),x_p(3)
integer ind(3)
real(amrex_real) c_000, c_001, c_010, c_011, c_100, c_101, c_110, c_111
real(amrex_real) c_00, c_01, c_10, c_11
real(amrex_real) c_0, c_1, c
real(8) :: input_to_IDINT
integer dir

if (nmat.eq.num_materials) then
 ! do nothing
else
 print *,"nmat invalid"
 stop
endif


LS(1)=radblob-sqrt((x(1)-xblob)**2+(x(2)-yblob)**2+(x(SDIM)-zblob)**2)
LS(2)=-LS(1)

if (axis_dir.eq.0) then !woven fabric

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

 input_to_IDINT=(x_p(1)-xblob3)/internal_dx(1)
 ind(1)=IDINT(input_to_IDINT)
 input_to_IDINT=(x_p(2)-yblob3)/internal_dx(2)
 ind(2)=IDINT(input_to_IDINT)
 if (SDIM.eq.3) then
  input_to_IDINT=(x_p(SDIM)-zblob3)/internal_dx(SDIM)
  ind(SDIM)=IDINT(input_to_IDINT)
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

else if (axis_dir.eq.1) then !mesh

 if(abs(x(SDIM)).ge.yblob6*2.0d0)then
  LS(3)=-yblob6*2.0d0
 else
  call MESH_LS(x,LS(3))
 endif
elseif(axis_dir.eq.2)then  ! single cylinder
 call CYL_LS(x,LS(3))   
else
 print *,"axis_dir invalid"
 stop
endif
! print *,"X, LS", x, LS

return
end subroutine FABRIC_DROP_LS


subroutine FABRIC_LS(x,LS)
use probcommon_module
IMPLICIT NONE


real(amrex_real), INTENT(in) :: x(SDIM)
real(amrex_real), INTENT(out) :: LS
real(amrex_real)         :: xlo,xhi,ylo,yhi,zlo,zhi
real(amrex_real)         :: hN1,hN2,xctemp
real(amrex_real)         :: xy(3),xc(3),xz(2),xh(3)
integer              :: i,ky,kx
real(amrex_real)         :: lstemp1,lstemp2
real(amrex_real)         :: a,b
real(amrex_real),allocatable :: ynf(:)  ! wavy thread
real(amrex_real),allocatable :: xnf(:)  ! straight thread
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
  enddo ! do ky=1,N1+2
 
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

subroutine MESH_LS(x_in,LS)
use probcommon_module
IMPLICIT NONE

real(amrex_real), INTENT(in) :: x_in(SDIM)
real(amrex_real), INTENT(out) :: LS
real(amrex_real)         :: lstemp1,lstemp2

real(amrex_real)         :: tempsx,tempsy
real(amrex_real)         :: centerx,centery,centerz
integer        :: TN1,TN2,i
real(amrex_real)         :: tx(SDIM),ty(SDIM)
real(amrex_real)         :: x(SDIM)

if(axis_dir.ne.1)then
 print *,"set axis_dir=1 for this validation test"
 stop
endif

do i=1,SDIM
 x(i)=x_in(i)
enddo

! center of the mesh, center of the very middle hole
centerx=0.0d0
centery=0.0d0
centerz=0.0d0
! xblob6=2b length of one side of the hole  
! yblob6=d thread diameter  
tempsx=xblob6+yblob6
tempsy=xblob6+yblob6
!print *,"tempsx",tempsx,"tempsy",tempsy

TN1 =floor((x(1)-centerx)/tempsx)
tx(1)=(real(TN1,8)+0.5d0)*tempsx
tx(2)=x(2)
tx(SDIM)=centerz
lstemp1=0.5d0*yblob6-sqrt((tx(1)-x(1))**2+(tx(SDIM)-x(SDIM))**2)
!lstemp1=0.5d0*yblob6-sqrt((tx(1)-x(1))**2.0d0+(tx(3)-x(3))**2.0d0)
LS=lstemp1

if(1.eq.1)then
TN2 =floor((x(2)-centery)/tempsy)
ty(2)=(real(TN2,8)+0.5d0)*tempsy
ty(1)=x(1)
ty(SDIM)=centerz
lstemp2=0.5d0*yblob6-sqrt((ty(2)-x(2))**2+(ty(SDIM)-x(SDIM))**2)


LS=0.0d0
if(lstemp1.lt.0.0d0.and.lstemp2.lt.0.0d0)then
 ls=max(lstemp1,lstemp2)
elseif(lstemp1.ge.0.0d0.and.lstemp2.lt.0.0d0)then
 ls=lstemp1
elseif(lstemp2.ge.0.0d0.and.lstemp1.lt.0.0d0)then
 ls=lstemp2
elseif(lstemp1.ge.0.0d0.and.lstemp2.ge.0.0d0)then
 ls=min(lstemp1,lstemp2)
else
 print *,"check thread ls setup" 
 print *,"lstemp1=",lstemp1
 print *,"lstemp2=",lstemp2
 stop
endif
endif

return
end subroutine MESH_LS

subroutine CYL_LS(x_in,LS)
use probcommon_module
IMPLICIT NONE

real(amrex_real), INTENT(in) :: x_in(SDIM)
real(amrex_real), INTENT(out) :: LS
real(amrex_real)         :: lstemp1

real(amrex_real)         :: centerx,centery,centerz
integer        :: i
real(amrex_real)         :: tx(SDIM)
real(amrex_real)         :: x(SDIM)

if(axis_dir.ne.2)then
 print *,"set aixs_dir=1 for this single CYL test"
 stop
endif

do i=1,SDIM
 x(i)=x_in(i)
enddo

! center of the mesh, center of the very middle hole
centerx=0.0d0
centery=0.0d0
centerz=0.0d0
! xblob6=2b length of one side of the hole  
! yblob6=d thread diameter  

! simple cylinder validation test
tx(1)=x(1)
tx(2)=centery
tx(SDIM)=centerz
lstemp1=0.5d0*yblob6-sqrt((tx(2)-x(2))**2+(tx(SDIM)-x(SDIM))**2)

return
end subroutine CYL_LS



subroutine FABRIC_DROP_VEL(x,t,LS,VEL,velsolid_flag,dx,nmat)
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
  print *,"dx invalid: ",dx(dir)
  stop
 endif
enddo

if (adv_dir.eq.SDIM) then
  ! material 3 is the fabric
  ! velsolid_flag==1 if initializing the solid velocity.
 if ((LS(3).ge.zero).or. &
     (velsolid_flag.eq.1)) then
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
   print *,"LS bust: ",LS(1),LS(2)
   stop
  endif
 else
  print *,"LS(3) or velsolid_flag bust: ",LS(3),velsolid_flag
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
!  integer, INTENT(in) :: imattype,im
!  real(amrex_real), INTENT(in) :: rho
!  real(amrex_real), INTENT(in) :: internal_energy
!  real(amrex_real), INTENT(out) :: pressure

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
!  integer, INTENT(in) :: imattype,im
!  real(amrex_real), INTENT(in) :: rho
!  real(amrex_real), INTENT(in) :: internal_energy
!  real(amrex_real), INTENT(out) :: soundsqr
!  real(amrex_real) pressure

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
!  integer, INTENT(in) :: imattype,im
!  real(amrex_real), INTENT(in) :: rho
!  real(amrex_real), INTENT(in) :: temperature 
!  real(amrex_real), INTENT(out) :: local_internal_energy

!  call INTERNAL_default(rho,temperature,local_internal_energy, &
!         imattype,im)

!  return
! end subroutine INTERNAL_FABRIC_DROP

! subroutine TEMPERATURE_FABRIC_DROP(rho,temperature,internal_energy, &
!   imattype,im)
!  use global_utility_module
!  IMPLICIT NONE
!  integer, INTENT(in) :: imattype,im
!  real(amrex_real), INTENT(in) :: rho
!  real(amrex_real), INTENT(out) :: temperature 
!  real(amrex_real), INTENT(in) :: internal_energy

!  call TEMPERATURE_default(rho,temperature,internal_energy, &
!         imattype,im)

!  return
! end subroutine TEMPERATURE_FABRIC_DROP

! This routine will not effect the simulation since
! all of the domain BC will be no-slip.


subroutine FABRIC_DROP_PRES(x,t,LS,PRES,nmat)
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
PRES=outflow_pressure

return 
end subroutine FABRIC_DROP_PRES

subroutine FABRIC_DROP_STATE(x,t,LS,STATE,bcflag,nmat,nstate_mat)
use probcommon_module
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

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: xwall
real(amrex_real), INTENT(in) :: xghost(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(inout) :: LS(nmat)
real(amrex_real), INTENT(in) :: LS_in(nmat)
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

integer, INTENT(in) :: nmat
real(amrex_real), INTENT(in) :: xwall
real(amrex_real), INTENT(in) :: xghost(SDIM)
real(amrex_real), INTENT(in) :: t
real(amrex_real), INTENT(in) :: LS(nmat)
real(amrex_real) local_STATE(nmat*num_state_material)
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

if ((num_materials.eq.3).and.(probtype.eq.FABRIC_DROP_PROB_TYPE)) then
 heat_source=zero
else
 print *,"num_materials or probtype invalid"
 stop
endif

return
end subroutine FABRIC_DROP_HEATSOURCE


subroutine FABRIC_DROP_SUMINT(GRID_DATA_IN,increment_out1, &
       increment_out2,nsum1,nsum2,isweep)
use probcommon_module_types
use probcommon_module
use global_utility_module
IMPLICIT NONE

integer, INTENT(in) :: nsum1,nsum2,isweep
type(user_defined_sum_int_type), INTENT(in) :: GRID_DATA_IN
real(amrex_real), INTENT(inout) :: increment_out1(nsum1)
real(amrex_real), INTENT(inout) :: increment_out2(nsum2)

integer :: level,finest_level

integer :: i,j,k,dir
real(amrex_real) :: FABRIC_DROP_LOW,FABRIC_DROP_HIGH
real(amrex_real) :: xval,yval,zval
real(amrex_real) :: LS_SOLID,F_LIQUID,volgrid,denom
real(amrex_real) :: xcentroid,ycentroid

i=GRID_DATA_IN%igrid
j=GRID_DATA_IN%jgrid
k=GRID_DATA_IN%kgrid
level=GRID_DATA_IN%level
finest_level=GRID_DATA_IN%finest_level

if ((level.le.finest_level).and.(level.ge.0)) then
 ! do nothing
else
 print *,"level invalid: ",level
 stop
endif

if (SDIM.eq.3) then
 !do nothing
else
 print *,"expecting 3d"
 stop
endif

if ((num_materials.eq.3).and. &
    (probtype.eq.FABRIC_DROP_PROB_TYPE)) then

  ! zero moment + x,y moment = 3 integrations (3D)
  ! above, in, and below the fabric: 9 integrations
  ! r^2 moment (3 integrations)
 if ((nsum1.eq.9).and.(nsum2.eq.3)) then

  if (axis_dir.eq.0) then
   FABRIC_DROP_LOW=-(a_wavy+r_1-0.005d0)
   FABRIC_DROP_HIGH=(a_wavy+r_1-0.005d0)
  else if ((axis_dir.eq.1).or. &
           (axis_dir.eq.2)) then
   FABRIC_DROP_LOW=-half*yblob6
   FABRIC_DROP_HIGH=half*yblob6
  else
   print *,"axis_dir invalid"
   stop
  endif

  if (isweep.eq.0) then
   do dir=1,9
    increment_out1(dir)=zero
   enddo
  else if (isweep.eq.1) then
   do dir=1,3
    increment_out2(dir)=zero
   enddo
  else
   print *,"isweep invalid: ",isweep
   stop
  endif

  xval=GRID_DATA_IN%xsten(0,1)
  yval=GRID_DATA_IN%xsten(0,2)
  zval=GRID_DATA_IN%xsten(0,SDIM)

  LS_SOLID=GRID_DATA_IN%lsfab(D_DECL(i,j,k),3)
  F_LIQUID=GRID_DATA_IN%slopes(D_DECL(i,j,k),1)
  volgrid=GRID_DATA_IN%volgrid

  if (LS_SOLID.ge.zero) then
   !do nothing
  else if (LS_SOLID.lt.zero) then

   if (zval.gt.FABRIC_DROP_HIGH) then
    dir=0
   else if (zval.ge.FABRIC_DROP_LOW) then
    dir=1
   else if (zval.lt.FABRIC_DROP_LOW) then
    dir=2
   else
    print *,"zval invalid: ",zval
    stop
   endif

   if (isweep.eq.0) then
    increment_out1(3*dir+1)=F_LIQUID*volgrid
    increment_out1(3*dir+2)=xval*F_LIQUID*volgrid
    increment_out1(3*dir+3)=yval*F_LIQUID*volgrid
   else if (isweep.eq.1) then
    denom=increment_out1(3*dir+1)
    xcentroid=zero
    ycentroid=zero
    if (denom.gt.0.0d0) then
     xcentroid=increment_out1(3*dir+2)/denom
     ycentroid=increment_out1(3*dir+3)/denom
     increment_out2(dir+1)= &
       ((xval-xcentroid)**2+(yval-ycentroid)**2)*volgrid*F_LIQUID
    else if (denom.eq.0.0d0) then
     ! do nothing
    else
     print *,"denom invalid: ",denom
     stop
    endif
   else
    print *,"isweep invalid: ",isweep
    stop
   endif
  else
   print *,"LS_SOLID invalid: ",LS_SOLID
   stop
  endif

 else
  print *,"nsum1 or nsum2 invalid: ",nsum1,nsum2
  stop
 endif

else
 print *,"num_materials ", num_materials
 print *,"probtype ", probtype
 print *,"num_materials or probtype invalid"
 stop
endif

end subroutine FABRIC_DROP_SUMINT

subroutine FABRIC_DROP_OVERRIDE_TAGFLAG( &
  i,j,k, &
  level,max_level, &
  snew_ptr,lsnew_ptr, &
  xsten,nhalf,time, &
  rflag,tagflag)
use amrex_fort_module, only : amrex_real
use probcommon_module
use global_utility_module
IMPLICIT NONE
integer, INTENT(in) :: i,j,k
integer, INTENT(in) :: level,max_level
integer, INTENT(in) :: nhalf
real(amrex_real), INTENT(in) :: xsten(-nhalf:nhalf,SDIM)
real(amrex_real), INTENT(in) :: time
real(amrex_real), INTENT(inout) :: rflag
integer, INTENT(inout) :: tagflag
real(amrex_real), INTENT(in),pointer :: snew_ptr(D_DECL(:,:,:),:)
real(amrex_real), INTENT(in),pointer :: lsnew_ptr(D_DECL(:,:,:),:)
real(amrex_real), dimension(3) :: local_x
real(amrex_real), dimension(SDIM) :: local_delta
real(amrex_real) :: F_LIQUID,LS_LIQUID,LS_FABRIC
integer :: dir

if (nhalf.lt.3) then
 print *,"nhalf invalid fabric drop override tagflag"
 stop
endif
if ((level.ge.0).and.(level.lt.max_level)) then
 ! do nothing
else
 print *,"level and/or max_level invalid"
 print *,"level=",level
 print *,"max_level=",max_level
 stop
endif
do dir=1,SDIM
 local_x(dir)=xsten(0,dir)
enddo
local_x(3)=xsten(0,SDIM)
do dir=1,SDIM
 local_delta(dir)=xsten(1,dir)-xsten(-1,dir)
 if (local_delta(dir).gt.zero) then
  ! do nothing
 else
  print *,"local_delta invalid fabric_drop_override_tagflag"
  stop
 endif
enddo !dir=1..sdim

if ((num_materials.eq.3).and. &
    (probtype.eq.FABRIC_DROP_PROB_TYPE)) then

 if ((axis_dir.ge.0).and.(axis_dir.le.2)) then

  rflag=0.0d0
  tagflag=0
  LS_LIQUID=lsnew_ptr(D_DECL(i,j,k),1)
  LS_FABRIC=lsnew_ptr(D_DECL(i,j,k),3)
  F_LIQUID=snew_ptr(D_DECL(i,j,k),STATECOMP_MOF+1)

  if (LS_FABRIC.ge.zero) then

   !do nothing

  else if (LS_FABRIC.lt.zero) then

   if (abs(LS_LIQUID).le.local_delta(1)) then
    rflag=1.0d0
    tagflag=1
   else if (abs(LS_LIQUID).ge.local_delta(1)) then
    !do nothing
   else
    print *,"LS_LIQUID invalid: ",LS_LIQUID
    print *,"local_delta(1) :",local_delta(1)
    stop
   endif
 
   if ((F_LIQUID.ge.0.1d0).and.(F_LIQUID.le.0.9d0)) then
    rflag=1.0d0
    tagflag=1
   else if ((F_LIQUID.ge.-EPS1).and.(F_LIQUID.le.0.1d0)) then
    ! do nothing
   else if ((F_LIQUID.ge.0.9d0).and.(F_LIQUID.le.1.0d0+EPS1)) then
    ! do nothing
   else
    print *,"F_LIQUID invalid: ",F_LIQUID
    stop
   endif
           
  else
   print *,"LS_FABRIC invalid: ",LS_FABRIC
   stop
  endif
 else
  print *,"axis_dir invalid: ",axis_dir
  stop
 endif
else
 print *,"num_materials or probtype invalid"
 print *,"num_materials: ",num_materials
 print *,"probtype: ",probtype
 stop
endif

end subroutine FABRIC_DROP_OVERRIDE_TAGFLAG

end module FABRIC_DROP_MODULE
