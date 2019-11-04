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

MODULE marangoni


CONTAINS

 ! >>>>>>>>>>>>>>>>>>position .vs. temperature
 ! heat pipe
subroutine position_Temp(flag,delta,delta2,x,y,z,Tout)  
!  See Kundan et al 2015 Langmuir  Figure 2.
!  For the heating side input 2.4W,  with temperature 450K
!  If flag = 0 ,     _ _ _ _ _ _ _ _ _
!                   | | | | | | | | | |
!                   | | | | | | | | | |       set delta to anything
!                   | | | | | | | | | |                     
!                   |_|_|_|_|_|_|_|_|_|
!                   
!  If flag = 1,     _ _ _ _ _ _ _ _ _
!                  |   ___________   | 
!                  |  |  room     |  |thickness of surrounding part is "delta"
!                  |  |___temp____|  |
!                  |_ __ _ _ _ _ _ _ |
!
!
! 0<y<3

implicit none

integer               :: flag
REAL_T                :: delta
REAL_T                :: delta2
REAL_T,intent(in)     :: x,y,z
REAL_T,parameter      :: p1 = -0.94381
REAL_T,parameter      :: p2 = 10.921
REAL_T,parameter      :: p3 = -43.593
REAL_T,parameter      :: p4 = 51.946
REAL_T,parameter      :: p5 = 73.773
REAL_T,parameter      :: p6 = -216.18
REAL_T,parameter      :: p7 = 465.15
REAL_T,parameter      :: TLO = 317.8415
REAL_T,parameter      :: THI = 465.15
REAL_T                :: Tout,T_expand

if (SDIM.eq.2) then
 if (abs(y-z).ge.1.0E-8) then
  print *,"expecting y=z"
  stop
 endif
endif

if (y.le.0.0d0) then
 T_expand=THI
else if (y.ge.3.0d0) then
 T_expand=TLO
else if ((y.ge.0.0d0).and.(y.le.3.0d0)) then

 T_expand=p1*(y**6.0d0) + p2*(y**5.0d0) + p3*(y**4.0d0) &           
         +p4*(y**3.0d0) + p5*(y**2.0d0) + p6*y + p7         
else
 print *,"y invalid"
 stop
endif

if (flag .eq. 0) then
 Tout=T_expand
elseif(flag .eq. 1) then
 if((y.gt.(3.0d0-delta2)).or. &
    (y.lt.delta2).or.  &
    (x.gt.(0.15d0-delta2)).or. &
    (x.lt.(-0.15d0+delta2))) then
  Tout=T_expand
 else 
  Tout = 293.0d0  
 endif  
else
 write(*,*) "Flag is not 0 or 1"
 stop
endif

end subroutine position_Temp



!  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>temperature .vs. surface tension

subroutine Temp_ST(T,r)
! temperature .vs.  surface tension
! input a temperature, output the corresponding surface tension
! unit :  temperature T:  K
!         Surface Tension r: dyne/cm   = 10e-3(N/m)

implicit none

REAL_T ,intent(in)   :: T
REAL_T               :: r


r =  (7.0d0-3.135d0)/(-40.0d0)*(T-383.0d0) + 7.0d0

end subroutine Temp_ST



!>>>>>>>>>>>>>>>>>>>>>>> distance function of the initial bubble
subroutine dist_long_bubble(delta,delta2,x,y,z,dist)
! distance function for the initial bubble
! output:  signed distance to the interface 
!---------------------------------------------------
!     two semicircle on the two ends
!     c1 and c2 are the center of the semicircle
!    
!                    +
!     a1 _______________________________ b1
!      /                                \
!   a2| c1           -               c2  | b2
!      \ _______________________________/                               
!     a3                                b3 
! 

implicit none

REAL_T,intent(in)        :: x,y,z
REAL_T                   :: delta   ! liquid film thickness
REAL_T                   :: delta2  ! init dist from top
REAL_T                   :: dist
REAL_T,dimension(SDIM)   :: c1,c2
REAL_T                   :: r       ! radius of the semicircle
REAL_T                   :: dist1
REAL_T,dimension(SDIM)   :: xvec

 if (SDIM.eq.2) then
  if (abs(y-z).ge.1.0E-8) then
   print *,"expecting y=z"
   stop
  endif
 endif
 xvec(1)=x
 xvec(2)=y
 if (SDIM.eq.3) then
  xvec(SDIM)=z
 endif

 r  = 0.15d0 - delta 

 c1(2) = 0.3d0 + r
 c1(1) = 0.0d0
 c2(2) = 3.0d0 - r - delta2
 c2(1) = 0.0d0


IF(y .lt. c1(2))then
    call cal_dist(xvec,c1,dist1)
    dist = dist1 - r
ELSEIF(y .gt. c2(2))THEN
    call cal_dist(xvec,c2,dist1)
    dist = dist1 - r
ELSE
    dist = abs(x) - r
ENDIF
  

end subroutine dist_long_bubble



subroutine cal_dist(x,y,dist)
! calculate the distance between two points in 2D plane

implicit none

REAL_T,intent(in)    :: x(SDIM),y(SDIM)
REAL_T,intent(out)   :: dist

dist = sqrt((x(1)-y(1))**2.0d0 + (x(2)-y(2))**2.0d0) 

end subroutine cal_dist


END MODULE marangoni
