      subroutine setExactPars(time,nIS,cIS,uIS,xID,dx)
      implicit none
c
c------------------------------------------------
c     An isolated shock moves at constant speed
c------------------------------------------------
c
      include 'parameters.h'
      include 'initial.h'
      include 'wave_constants.h'
      include 'exact.h'

c     Input:
      real*8 time
      integer nIS
      real*8 cIS(0:N_INIT_STATE,1:N_EQUATIONS)
      real*8 uIS(0:N_INIT_STATE,1:N_EXT_VARS)
      real*8 xID(0:N_INIT_DISC)
      real*8 dx
c     Auxiliary:
      integer i,iVar
      integer dir
      real*8 primExtSt2W

      lCSt=1
      rCSt=2

      timeInit = time
      
      nStateInit = nIS

      do i=1,nIS+1

c        i-th Discontinuity 
         xInit(i) = xID(i)

         if( i .eq. nIS+1 ) goto 10

c        i-th State
c           Extended Primitive Variables
         do iVar=1,N_EXT_VARS
            uInit(i,iVar) = uIS(i,iVar)
         enddo
c           Conservative Variables
         do iVar=1,N_EQUATIONS
            cInit(i,iVar) = cIS(i,iVar)
         enddo 

      enddo
 10   continue
      
      if ( uInit(1,IPRE) .lt. uInit(2,IPRE) ) then
         !left-facing wave
         dir = -1 
      else
         !rigth-facing wave
         dir = +1 
      endif
      wInit(2) = primExtSt2W(
     &     uInit(1,IRHO),uInit(1,IVEL),uInit(1,IPRE),uInit(1,ISND),
     &     uInit(2,IRHO),uInit(2,IVEL),uInit(2,IPRE),uInit(2,ISND),dir)
      wInit(1) = 0.d0
      wInit(3) = 0.d0

      dxInit = dx

      end
