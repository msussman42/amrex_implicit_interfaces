      program main
      IMPLICIT NONE

      real*8 mypi
      real*8 rad,theta,sin1,sin2,d
      real*8 rtheta,e0

      mypi=4.d0*atan(1.d0)
      theta=mypi/3.d0  ! 60 degrees
      rtheta=sqrt(mypi/(2.d0*theta-sin(2.d0*theta)))
      d=2.d0*sin(mypi-theta)*rtheta
      e0=(1.d0+cos(mypi-theta))*rtheta

      print *,"theta= ",theta
      print *,"Rtheta=",rtheta
      print *,"d= ",d
      print *,"e0= ",e0
  
      end

