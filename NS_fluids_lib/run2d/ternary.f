      program main
      IMPLICIT NONE

      real*8 sigma12,sigma13,sigma23
      real*8 mypi
      real*8 rad,theta1,theta2,sin1,sin2,d
      real*8 term1,term2
      real*8 rtheta,e0
      real*8 gamma1,gamma2,gamma3
      real*8 pos1(3),pos2(3),distpos1pos2
      real*8 virtual_radius
      real*8 minor_axis
      integer dir
      integer select

      select=6

c sin theta1/sigma23 = sin theta2/sigma13 = sin theta 3/sigma12
      if (select.eq.1) then
       sigma12=0.044444444444444d0
       sigma23=sigma12
       sigma13=0.055555555555555d0
      else if (select.eq.2) then
       sigma12=0.044444444444444d0
       sigma23=sigma12
       sigma13=0.066666666666666d0
      else if (select.eq.3) then
       sigma12=0.044444444444444d0
       sigma23=sigma12
       sigma13=0.044444444444444d0
      else if (select.eq.4) then
       sigma12=0.044444444444444d0
       sigma23=sigma12
       sigma13=0.033333333333333d0
      else if (select.eq.5) then
       sigma12=0.7555d0
       sigma13=1.0d0
       sigma23=0.5495d0
      else if (select.eq.6) then
       sigma12=1.0d0
       sigma23=sigma12
       sigma13=1.25d0
      else
       print *,"select invalid: ",select
       stop
      endif

      gamma1=0.5d0*(sigma12-sigma23+sigma13)
      gamma2=0.5d0*(sigma12+sigma23-sigma13)
      gamma3=0.5d0*(-sigma12+sigma23+sigma13)
      mypi=4.d0*atan(1.d0)
      
      if ((select.ge.1).and.(select.le.4)) then 
       rad=0.15d0
       theta1=acos(-sigma13/(2.d0*sigma23))
       theta2=2.d0*mypi-2.d0*theta1
       sin1=sin(theta2)
c      sin2=sin(0.5d0*theta2)
       sin2=sin(mypi-theta1)
       d=(theta2-sin1)/(4.d0*mypi*(rad**2)*(sin2**2))
       d=1.d0/sqrt(d)
       virtual_radius=d/(2.0d0*sin2)
       minor_axis=2.0d0*virtual_radius*(1.0d0-cos(0.5d0*theta2))
       print *,"gamma1,gamma2,g1/g2 (scale) ", 
     &  gamma1/sigma12,gamma2/sigma12,gamma1/gamma2
       print *,"gamma1,gamma3,g1/g3 (scale) ", 
     &  gamma1/sigma13,gamma3/sigma13,gamma1/gamma3
       print *,"gamma2,gamma3,g2/g3 (scale) ", 
     &  gamma2/sigma23,gamma3/sigma23,gamma2/gamma3
       print *,"gamma1,gamma2,gamma3 ",gamma1,gamma2,gamma3
       print *,"sigma12,sigma23,sigma13 ",sigma12,sigma23,sigma13
       print *,"d= ",d
       print *,"virtual_radius= ",virtual_radius
       print *,"minor_axis= ",minor_axis

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

       pos1(1)=0.2745d0
       pos1(2)=0.5012d0
       pos1(3)=0.0d0
       pos2(1)=0.7256d0
       pos2(2)=0.508d0
       pos2(3)=0.0d0
       distpos1pos2=0.0d0
       do dir=1,3
        distpos1pos2=distpos1pos2+(pos1(dir)-pos2(dir))**2
       enddo
       distpos1pos2=sqrt(distpos1pos2)
       print *,"pos1 ",pos1(1),pos1(2),pos1(3)
       print *,"pos2 ",pos2(1),pos2(2),pos2(3)
       print *,"distpos1pos2 ",distpos1pos2

      else if ((select.ge.5).and.(select.le.6)) then
       rad=1.0d0
       theta1=acos((sigma12**2+sigma13**2-sigma23**2)/
     &    (2.0d0*sigma12*sigma13))
       theta2=acos((sigma23**2+sigma13**2-sigma12**2)/
     &    (2.0d0*sigma23*sigma13))
       term1=(theta1/sin(theta1)-cos(theta1))/sin(theta1)
       term2=(theta2/sin(theta2)-cos(theta2))/sin(theta2)
       d=2.0d0*sqrt(mypi/(term1+term2))
       minor_axis=0.5d0*d*( (1.0d0-cos(theta1))/sin(theta1)+
     &   (1.0d0-cos(theta2))/sin(theta2) )
       print *,"d= ",d
       print *,"minor_axis= ",minor_axis
      else
       print *,"select invalid: ",select
       stop
      endif

      end
