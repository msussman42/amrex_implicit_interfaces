      program main
      IMPLICIT NONE

      real*8 sigma12,sigma13,sigma23
      real*8 theta1,theta2,theta3
      real*8 mypi
      integer select

      select=3

      if (select.eq.1) then
       sigma12=0.831d0   ! liquid/gas
       theta1=96.0d0 ! liquid
       theta2=168.0d0 ! gas
       theta3=96.0d0 ! solid
      else if (select.eq.2) then
       sigma12=0.831d0   ! liquid/gas
       theta1=107.5d0 ! liquid
       theta2=145.0d0 ! gas
       theta3=107.5d0 ! solid
      else if (select.eq.3) then
       sigma12=0.831d0   ! (silicone) liquid/gas
       theta1=60.0d0 ! (silicone) liquid
       theta2=168.0d0 ! gas
       theta3=132.0d0 ! solid
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
