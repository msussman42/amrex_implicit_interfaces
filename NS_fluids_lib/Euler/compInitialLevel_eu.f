      subroutine compInitialLevel(time,nIS,xID,uIS,cIS)
      implicit none

      include 'initial.h'
      include 'wave_constants.h'

c     Input:
      real*8 time
      integer nIS
      real*8 xID(0:N_INIT_DISC)
c     Input/Output:
      real*8 uIS(0:N_INIT_STATE,1:N_EXT_VARS)
c     Output: 
      real*8 cIS(0:N_INIT_STATE,1:N_EQUATIONS)
c     Auxiliary:      
      integer i
      real*8 thermSt2snd

      do i=1,nIS
c        Speed of Sound         
         uIS(i,ISND) = thermSt2snd(uIS(i,IRHO),uIS(i,IPRE))

c        Conservative Variables
         call stPrim2stCons(uIS(i,IRHO),uIS(i,IVEL),uIS(i,IPRE), 
     &        cIS(i,IRHO),cIS(i,IMOM),cIS(i,IENE))
      enddo
      end
