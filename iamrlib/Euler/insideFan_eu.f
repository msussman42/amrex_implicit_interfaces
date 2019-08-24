      subroutine insideFan(rhoC,velC,preC,sndC,dir,centSide,DxDt, 
     &     rho,vel,pre,snd)
      implicit none
c
c
c     Computes state W = {rho,vel,pre; ... snd}; 
c         state W lies INSIDE a FAN 
c             --- ON an EXPANSION portion of the wave-curve;
c         the wave-curve is centered at state C.
c
c     C is the state at which the wave-curve is centered.
c     With respect to the wave, 
c     state C could be either Ahead or Behind.
c
c     Whether C is Ahead or Behind 
c     depends on the value of centSide
c             centSide =  -1 state C is Behind the wave
c                         +1 state C  is Ahead of the wave
c
c-----------------------------------------------------------------------
c                              INPUT 
c-----------------------------------------------------------------------
c
c    STATE C (AHEAD of or BEHIND THE FAN)
c    --------------------
c  rhoC      density
c  velC      velocity
c  preC      pressure
c  sndC      sound speed
c
c    DxDt --- slope of the characteristic inside the fan on which 
c             the STATE is to be calculated
c
c    dir =1 for right-facing wave, =-1 for left-facing wave
c
c-----------------------------------------------------------------------
c                              OUTPUT 
c-----------------------------------------------------------------------
c
c    STATE inside the FAN on the characteristic with the slope DxDt
c    ---------------------------------------------------------------
c  rho    density                         
c  vel    velocity 
c  pre    pressure
c  snd    sound speed
c

      include 'mat_constants.h'

c     Input:
      real*8  rhoC,velC,preC,sndC
      integer dir,centSide
      real*8 DxDt
c     Output:
      real*8 rho, vel, pre, snd
c     Auxiliary:
      real*8 compSnd, compPower

      call slopeInFan2velNsnd(velC,sndC,dir,DxDt, vel,snd)

      compSnd = snd/sndC
      pre=preC*( compPower(compSnd,rc2) )
      rho=rhoC*(pre/preC)*( compPower(compSnd,-2.d0) ) 

      end
