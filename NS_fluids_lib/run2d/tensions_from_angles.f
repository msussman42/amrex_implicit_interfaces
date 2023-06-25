      program main
      IMPLICIT NONE

      real*8 sigma12,sigma13,sigma23
      real*8 theta1,theta2,theta3
      real*8 mypi
      integer select

      select=2

      if (select.eq.1) then
       sigma12=0.831d0   ! (silicone) liquid/gas
       theta1=60.0d0 ! (silicone) liquid
       theta2=168.0d0 ! gas
       theta3=132.0d0 ! ice
      else if (select.eq.2) then
       sigma12=1740.0d0   ! liquid/gas
       theta1=90.0d0+18.0d0 ! (silicone) liquid
       theta2=180.0d0-18.0d0 ! gas
       theta3=90.0d0 ! ice
      else
       print *,"select invalid"
       stop
      endif

      mypi=4.d0*atan(1.d0)

      if (abs(360.0d0-theta1-theta2-theta3).gt.1.0d-8) then
       print *,"must have: theta1+theta2+theta3=360"
       stop
      endif
      print *,"theta1,theta2,theta3 ",theta1,theta2,theta3
      theta1=theta1*mypi/180.0d0
      theta2=theta2*mypi/180.0d0
      theta3=theta3*mypi/180.0d0
      sigma13=sigma12*sin(theta2)/sin(theta3)
      sigma23=sigma12*sin(theta1)/sin(theta3)
      print *,"sigma12,sigma13,sigma23 ",sigma12,sigma13,sigma23
      end

