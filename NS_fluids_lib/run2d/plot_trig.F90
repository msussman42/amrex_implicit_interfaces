      program main
      IMPLICIT NONE
      integer nplot,i
      real*8 x_low
      real*8 x_high
      real*8 x_local
      real*8 S2
      real*8 fx
      real*8 local_pi

      local_pi=4.0d0*atan(1.0d0)

      nplot=100
      x_low=-local_pi
      x_high=local_pi

      do i=0,nplot
       x_local=x_low+(x_high-x_low)*i/nplot
       S2=0.25d0*local_pi-(2.0d0/local_pi)*cos(x_local)+sin(x_local)
       if (x_local.le.0.0d0) then
        fx=0.0d0
       else
        fx=x_local
       endif
       print *,"x,S2(x),f(x) ",x_local,S2,fx
      enddo

      return
      end

