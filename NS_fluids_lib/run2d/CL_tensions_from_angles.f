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
       theta_L=45.0d0 ! degrees
      else if (select.eq.2) then
       sigma12=0.831d0   ! liquid/gas
       sigma23=0.831d0  ! gas/solid
       theta_L=33.0d0 ! degrees
      else
       print *,"select invalid"
       stop
      endif

      mypi=4.d0*atan(1.d0)

c sigma_GS-sigma_LS=sigma_LG * cos(theta_L)
      sigma13=sigma23-sigma12*dcos(mypi*theta_L/180.0d0)
      print *,"sigma12=",sigma12
      print *,"sigma23=",sigma23
      print *,"sigma13=",sigma13
      print *,"theta_L= (degrees) ",theta_L
      end

