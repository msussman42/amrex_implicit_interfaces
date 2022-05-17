      program main
      IMPLICIT NONE

      real*8 sigma12,sigma13,sigma23
      real*8 mypi
      real*8 theta_L
      integer select

      select=2

      if (select.eq.1) then
       sigma12=0.831d0   ! liquid/gas
       sigma23=0.831d0  ! gas/solid
       sigma13=0.1341d0 ! liquid/solid
      else if (select.eq.2) then
       sigma12=0.831d0   ! liquid/gas
       sigma23=0.831d0  ! gas/solid
       sigma13=0.243d0 ! liquid/solid
      else
       print *,"select invalid"
       stop
      endif

      mypi=4.d0*atan(1.d0)

c sigma_GS-sigma_LS=sigma_LG * cos(theta_L)
      theta_L=dacos((sigma23-sigma13)/sigma12)
      theta_L=theta_L*180.0d0/mypi
      print *,"sigma12=",sigma12
      print *,"sigma23=",sigma23
      print *,"sigma13=",sigma13
      print *,"theta_L= (degrees) ",theta_L
      end

