      program main
      IMPLICIT NONE

      real*8 sigma12,sigma13,sigma23
      real*8 mypi
      real*8 rad,theta1,theta2,sin1,sin2,d
      real*8 rtheta,e0
      real*8 gamma1,gamma2,gamma3
      integer select

      select=2

      if (select.eq.1) then
       sigma12=0.044444444444444d0
       sigma23=sigma12
       sigma13=0.055555555555555d0
      else if (select.eq.2) then
       sigma12=0.044444444444444d0
       sigma23=sigma12
       sigma13=0.066666666666666d0
      endif

      gamma1=0.5d0*(sigma12-sigma23+sigma13)
      gamma2=0.5d0*(sigma12+sigma23-sigma13)
      gamma3=0.5d0*(-sigma12+sigma23+sigma13)
       
      rad=0.15d0
      mypi=4.d0*atan(1.d0)
      theta1=acos(-sigma13/(2.d0*sigma23))
      theta2=2.d0*mypi-2.d0*theta1
      sin1=sin(theta2)
      sin2=sin(mypi-theta1)
      d=(theta2-sin1)/(4.d0*mypi*(rad**2)*(sin2**2))
      d=1.d0/sqrt(d)
      print *,"gamma1,gamma2,g1/g2 (scale) ", 
     &  gamma1/sigma12,gamma2/sigma12,gamma1/gamma2
      print *,"gamma1,gamma3,g1/g3 (scale) ", 
     &  gamma1/sigma13,gamma3/sigma13,gamma1/gamma3
      print *,"gamma2,gamma3,g2/g3 (scale) ", 
     &  gamma2/sigma23,gamma3/sigma23,gamma2/gamma3
      print *,"gamma1,gamma2,gamma3 ",gamma1,gamma2,gamma3
      print *,"sigma12,sigma23,sigma13 ",sigma12,sigma23,sigma13
      print *,"d= ",d

      theta1=acos(sigma13/(2.d0*sigma23))
      rad=0.15d0
      rtheta=sqrt(mypi/(2.d0*theta1-sin(2.d0*theta1)))*rad
      d=2.d0*sin(mypi-theta1)*rtheta
      e0=2.d0*(1.d0+cos(mypi-theta1))*rtheta
      print *,"using formulas from Arienti and Sussman"
      print *,"theta= ",theta1
      print *,"theta in degrees= ",theta1*180.0/mypi
      print *,"R0=    ",rad
      print *,"Rtheta=",rtheta
      print *,"d=     ",d
      print *,"e0=    ",e0
  
      end

