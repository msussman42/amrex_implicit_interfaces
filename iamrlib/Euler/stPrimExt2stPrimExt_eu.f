      subroutine stPrimExt2stPrimExt(rhoIn,velIn,preIn,sndIn, 
     &     rhoOut,velOut,preOut,sndOut)
      implicit none

c     Input:
      real*8 rhoIn,velIn,preIn,sndIn
c     Output:
      real*8 rhoOut,velOut,preOut,sndOut
c-----
           
      rhoOut = rhoIn
      velOut = velIn
      preOut = preIn
      sndOut = sndIn

      end
