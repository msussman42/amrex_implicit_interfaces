      program main
      IMPLICIT NONE
      integer nplot,i
      real*8 temp_low
      real*8 temp_high
      real*8 temp_local
      real*8 sigma

      nplot=100
      temp_low=240.0
      temp_high=290.0 !room temperature

      do i=0,nplot
       temp_local=temp_low+(temp_high-temp_low)*i/nplot
       sigma=235.8*((1.0-temp_local/647.096)**1.256)* &
               (1.0-0.625*(1.0-temp_local/647.096))
       print *,"temp,sigma ",temp_local,sigma
      enddo

      return
      end

