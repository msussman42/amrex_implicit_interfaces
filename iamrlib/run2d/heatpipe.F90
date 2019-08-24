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
REAL*8                :: delta
REAL*8                :: delta2
REAL*8,intent(in)     :: x,y,z
REAL*8,parameter      :: p1 = -0.94381
REAL*8,parameter      :: p2 = 10.921
REAL*8,parameter      :: p3 = -43.593
REAL*8,parameter      :: p4 = 51.946
REAL*8,parameter      :: p5 = 73.773
REAL*8,parameter      :: p6 = -216.18
REAL*8,parameter      :: p7 = 465.15
REAL*8,parameter      :: TLO = 317.8415
REAL*8,parameter      :: THI = 465.15
REAL*8                :: Tout,T_expand

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

END MODULE marangoni

program test_marangoni
use marangoni

real*8 x,y,z,Tout,delta
integer flag

flag=0
delta=0.0
x=0.0
y=0.0
z=0.0
y=1.0E-2
call position_Temp(flag,delta,delta,x,y,z,Tout)
print *,"Tout=",Tout

end program
