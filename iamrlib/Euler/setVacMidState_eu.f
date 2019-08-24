      subroutine setVacMidState(rhoA,velA,preA,sndA,dir, 
     &     rhoB,velB,preB,sndB)
      implicit none  

      include 'parameters.h'

c     Input:
      real*8 rhoA,velA,preA,sndA
      integer dir
c     Output:
      real*8 rhoB,velB,preB,sndB 
c     Auxiliary.
      real*8 veljump,veljump2vel

      include 'mat_constants.h'

      rhoB = 0.d0
      preB = 0.d0
      sndB = 0.d0

      veljump = -rc1*sndA
      velB = veljump2vel(velA,veljump,dir)

      end
