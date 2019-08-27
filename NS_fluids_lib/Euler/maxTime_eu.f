      function maxTime(vel,snd,dx)
      implicit none

      include 'mat_constants.h'

c     Input:
      real*8 vel,snd,dx
c     Output:
      real*8 maxTime

      maxTime = dx/( dabs(vel)+snd )

      end
