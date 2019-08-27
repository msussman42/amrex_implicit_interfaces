      subroutine slopeInFan2velNsnd(velC,sndC,dir,DxDt,vel,snd)
      implicit none
c 
c-----------------------------------------------------------------------
c                              INPUT 
c-----------------------------------------------------------------------
c
c    STATE C (AHEAD of or BEHIND THE FAN)
c  velC      velocity
c  sndC      sound speed
c
c     dir =1 for right-facing wave, =-1 for left-facing wave
c
c     DxDt --- slope of the ray on which the STATE is to be calculated
c
c-----------------------------------------------------------------------
c                              OUTPUT 
c-----------------------------------------------------------------------
c
c    STATE inside the FAN on the characteristic with the slope DxDt
c    ---------------------------------------------------------------
c  vel    velocity 
c  snd    sound speed
c
c-----------------------------------------------------------------------
c
      include 'mat_constants.h'

c     Input:
      real*8  velC,sndC
      integer dir
      real*8 DxDt
c     Output:
      real*8 vel,snd

      vel=gk*( velC+rc1*(DxDt-dir*sndC) )
      snd=(DxDt-vel)*dir

      end
