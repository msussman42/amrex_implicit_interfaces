      subroutine setRegMidState(rhoA,velA,preA,sndA,vel,pre, 
     &     rhoB,velB,preB,sndB)
      implicit none

      include 'wave_constants.h'  

c     Input:
      real*8 rhoA,velA,preA,sndA, vel,pre
c     Output:
      real*8 rhoB,velB,preB,sndB
c     Auxiliary:
      real*8 onWaveCurve, thermSt2snd
      external pre2densExpan,pre2densCompr
      real*8 pre2densExpan,pre2densCompr
c-----
           
      velB = vel
      preB = pre
      rhoB = onWaveCurve(preB,preA,rhoA,AHEAD,
     &              pre2densExpan,pre2densCompr)
      sndB = thermSt2snd(rhoB,preB)

      end
